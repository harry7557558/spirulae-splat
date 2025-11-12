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

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_9)
{
    float _S87 = (*right_8).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_9.x;
    float sum_14 = _S87 + (*right_8).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_9.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_9.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S88 = (*right_8).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_9.x;
    float sum_16 = _S88 + (*right_8).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_9.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_9.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S89 = (*right_8).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_9.x;
    float sum_18 = _S89 + (*right_8).primal_0.rows[int(2)].y * dOut_9.y;
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

inline __device__ bool generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_4, float3  * ray_o_0, float3  * ray_d_0)
{
    float2  _S90 = uv_4;
    *ray_o_0 = - mul_7(t_1, R_2);
    bool _S91 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S90);
    if(!_S91)
    {
        return false;
    }
    float3  raydir_1;
    if(is_fisheye_2)
    {
        float2  _S92 = _S90;
        float theta_2 = length_0(_S90);
        float _S93;
        if(theta_2 < 0.00100000004749745f)
        {
            _S93 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S93 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S94 = make_float3 ((_S92 * make_float2 (_S93)).x, (_S92 * make_float2 (_S93)).y, (F32_cos((theta_2))));
        raydir_1 = _S94;
    }
    else
    {
        raydir_1 = make_float3 (_S90.x, _S90.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_1, R_2));
    return true;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S95;
    bool _S96;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S97, Matrix<float, 3, 3>  _S98)
{
    return mul_7(_S97, _S98);
}

inline __device__ float s_primal_ctx_sin_0(float _S99)
{
    return (F32_sin((_S99)));
}

inline __device__ float s_primal_ctx_cos_0(float _S100)
{
    return (F32_cos((_S100)));
}

inline __device__ bool s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_5, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S95 = make_float2 (0.0f);
    _s_diff_ctx_0->_S96 = false;
    float3  _S101 = make_float3 (0.0f);
    float3  _S102 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    float2  _S103 = uv_5;
    bool _S104 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S103);
    _s_diff_ctx_0->_S95 = _S103;
    _s_diff_ctx_0->_S96 = _S104;
    float2  _S105 = _S103;
    float3  raydir_2;
    bool _S106;
    if(!!_S104)
    {
        if(is_fisheye_3)
        {
            float _S107 = length_0(_S105);
            float _S108;
            if(_S107 < 0.00100000004749745f)
            {
                _S108 = 1.0f - _S107 * _S107 / 6.0f;
            }
            else
            {
                _S108 = s_primal_ctx_sin_0(_S107) / _S107;
            }
            float3  _S109 = make_float3 ((_S105 * make_float2 (_S108)).x, (_S105 * make_float2 (_S108)).y, s_primal_ctx_cos_0(_S107));
            raydir_2 = _S109;
        }
        else
        {
            raydir_2 = make_float3 (_S105.x, _S105.y, 1.0f);
        }
        float3  _S110 = normalize_0(s_primal_ctx_mul_0(raydir_2, dpR_0));
        _S106 = true;
        raydir_2 = _S110;
    }
    else
    {
        _S106 = false;
        raydir_2 = _S101;
    }
    *dpray_o_0 = _S102;
    *dpray_d_0 = raydir_2;
    return _S106;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S111, float _S112)
{
    _d_sqrt_0(_S111, _S112);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float _s_dOut_0)
{
    float _S113 = (*dpx_5).primal_0.x;
    float _S114 = (*dpx_5).primal_0.y;
    float _S115 = (*dpx_5).primal_0.z;
    DiffPair_float_0 _S116;
    (&_S116)->primal_0 = _S113 * _S113 + _S114 * _S114 + _S115 * _S115;
    (&_S116)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S116, _s_dOut_0);
    float _S117 = (*dpx_5).primal_0.z * _S116.differential_0;
    float _S118 = _S117 + _S117;
    float _S119 = (*dpx_5).primal_0.y * _S116.differential_0;
    float _S120 = _S119 + _S119;
    float _S121 = (*dpx_5).primal_0.x * _S116.differential_0;
    float _S122 = _S121 + _S121;
    float3  _S123 = make_float3 (0.0f);
    *&((&_S123)->z) = _S118;
    *&((&_S123)->y) = _S120;
    *&((&_S123)->x) = _S122;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S123;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S124, float _S125)
{
    s_bwd_prop_length_impl_0(_S124, _S125);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  _s_dOut_1)
{
    float _S126 = length_1((*dpx_6).primal_0);
    float3  _S127 = (*dpx_6).primal_0 * _s_dOut_1;
    float3  _S128 = make_float3 (1.0f / _S126) * _s_dOut_1;
    float _S129 = - ((_S127.x + _S127.y + _S127.z) / (_S126 * _S126));
    float3  _S130 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S131;
    (&_S131)->primal_0 = (*dpx_6).primal_0;
    (&_S131)->differential_0 = _S130;
    s_bwd_length_impl_0(&_S131, _S129);
    float3  _S132 = _S128 + _S131.differential_0;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S132;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S133, float3  _S134)
{
    s_bwd_prop_normalize_impl_0(_S133, _S134);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S135, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S136, float3  _S137)
{
    _d_mul_1(_S135, _S136, _S137);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_6, bool is_fisheye_4, FixedArray<float, 10>  * dist_coeffs_6, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S138 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S139 = *dpt_1;
    float3  _S140 = make_float3 (0.0f);
    bool _S141 = !!_s_diff_ctx_1->_S96;
    float3  raydir_3;
    float3  _S142;
    if(_S141)
    {
        if(is_fisheye_4)
        {
            float _S143 = length_0(_s_diff_ctx_1->_S95);
            float _S144;
            if(_S143 < 0.00100000004749745f)
            {
                _S144 = 1.0f - _S143 * _S143 / 6.0f;
            }
            else
            {
                _S144 = s_primal_ctx_sin_0(_S143) / _S143;
            }
            float3  _S145 = make_float3 ((_s_diff_ctx_1->_S95 * make_float2 (_S144)).x, (_s_diff_ctx_1->_S95 * make_float2 (_S144)).y, s_primal_ctx_cos_0(_S143));
            raydir_3 = _S145;
        }
        else
        {
            raydir_3 = make_float3 (_s_diff_ctx_1->_S95.x, _s_diff_ctx_1->_S95.y, 1.0f);
        }
        float3  _S146 = raydir_3;
        raydir_3 = s_primal_ctx_mul_0(raydir_3, _S138.primal_0);
        _S142 = _S146;
    }
    else
    {
        raydir_3 = _S140;
        _S142 = _S140;
    }
    Matrix<float, 3, 3>  _S147 = makeMatrix<float, 3, 3> (0.0f);
    Matrix<float, 3, 3>  _S148;
    if(_S141)
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S149;
        (&_S149)->primal_0 = raydir_3;
        (&_S149)->differential_0 = _S140;
        s_bwd_normalize_impl_0(&_S149, dpray_d_1);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S150;
        (&_S150)->primal_0 = _S142;
        (&_S150)->differential_0 = _S140;
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S151;
        (&_S151)->primal_0 = _S138.primal_0;
        (&_S151)->differential_0 = _S147;
        s_bwd_prop_mul_0(&_S150, &_S151, _S149.differential_0);
        _S148 = _S151.differential_0;
    }
    else
    {
        _S148 = _S147;
    }
    float3  _S152 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S153;
    (&_S153)->primal_0 = _S139.primal_0;
    (&_S153)->differential_0 = _S140;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S154;
    (&_S154)->primal_0 = _S138.primal_0;
    (&_S154)->differential_0 = _S147;
    s_bwd_prop_mul_0(&_S153, &_S154, _S152);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S153.differential_0;
    Matrix<float, 3, 3>  _S155 = _S154.differential_0 + _S148;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S155;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S156, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S157, float2  _S158, bool _S159, FixedArray<float, 10>  * _S160, float3  _S161, float3  _S162)
{
    float3  _S163;
    float3  _S164;
    s_bwd_prop_generate_ray_Intermediates_0 _S165;
    bool _S166 = s_primal_ctx_generate_ray_0((*_S156).primal_0, (*_S157).primal_0, _S158, _S159, _S160, &_S163, &_S164, &_S165);
    s_bwd_prop_generate_ray_Intermediates_0 _S167 = _S165;
    s_bwd_prop_generate_ray_0(_S156, _S157, _S158, _S159, _S160, _S161, _S162, &_S167);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_7, bool is_fisheye_5, FixedArray<float, 10>  * dist_coeffs_7, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S168 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S168;
    float3  _S169 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S169;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_7, is_fisheye_5, dist_coeffs_7, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S170 = float(width_0);
    float _S171 = float(height_0);
    float _S172 = 0.30000001192092896f * (0.5f * _S170 / fx_0);
    float _S173 = 0.30000001192092896f * (0.5f * _S171 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_2 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S170 - cx_0) / fx_0 + _S172), ((F32_max((- (cx_0 / fx_0 + _S172)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S171 - cy_0) / fy_0 + _S173), ((F32_max((- (cy_0 / fy_0 + _S173)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_2, cov3d_0), transpose_1(J_2));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_10)
{
    DiffPair_float_0 _S174 = *dpx_7;
    bool _S175;
    if(((*dpx_7).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S175 = ((*dpx_7).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S175 = false;
    }
    float _S176;
    if(_S175)
    {
        _S176 = dOut_10;
    }
    else
    {
        _S176 = 0.0f;
    }
    dpx_7->primal_0 = _S174.primal_0;
    dpx_7->differential_0 = _S176;
    DiffPair_float_0 _S177 = *dpMin_0;
    if((_S174.primal_0) < ((*dpMin_0).primal_0))
    {
        _S176 = dOut_10;
    }
    else
    {
        _S176 = 0.0f;
    }
    dpMin_0->primal_0 = _S177.primal_0;
    dpMin_0->differential_0 = _S176;
    DiffPair_float_0 _S178 = *dpMax_0;
    if(((*dpx_7).primal_0) > ((*dpMax_0).primal_0))
    {
        _S176 = dOut_10;
    }
    else
    {
        _S176 = 0.0f;
    }
    dpMax_0->primal_0 = _S178.primal_0;
    dpMax_0->differential_0 = _S176;
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

inline __device__ bool persp_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_8, uint width_1, uint height_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float2  * _S179;
    float2  * _S180;
    float2  * _S181;
    bool _S182;
    float2  * _S183;
    float2  * _S184;
    float2  * _S185;
    bool _S186;
    float2  * _S187;
    bool _S188;
    int2  _S189 = make_int2 (int(0));
    float2  _S190 = make_float2 ((float)_S189.x, (float)_S189.y);
    *mean2d_1 = _S190;
    *cov2d_1 = makeMatrix<float, 2, 2> (0.0f);
    float fx_1 = intrins_0.x;
    float fy_1 = intrins_0.y;
    float cx_1 = intrins_0.z;
    float cy_1 = intrins_0.w;
    float _S191 = float(width_1);
    float _S192 = float(height_1);
    float _S193 = 0.30000001192092896f * (0.5f * _S191 / fx_1) * fx_1;
    float lim_x_pos_0 = _S191 + _S193;
    float lim_x_neg_0 = cx_1 + _S193;
    float _S194 = 0.30000001192092896f * (0.5f * _S192 / fy_1) * fy_1;
    float lim_y_pos_0 = _S192 + _S194;
    float lim_y_neg_0 = cy_1 + _S194;
    FixedArray<float2 , 7>  proj_points_0;
    for(;;)
    {
        bool _S195;
        _S179 = &proj_points_0[int(0)];
        for(;;)
        {
            float _S196 = sigmas_0->p_0[int(0)].z;
            proj_points_0[int(0)] = float2 {sigmas_0->p_0[int(0)].x, sigmas_0->p_0[int(0)].y} / make_float2 (_S196);
            if(_S196 < 0.0f)
            {
                _S195 = true;
            }
            else
            {
                bool _S197 = is_valid_distortion(proj_points_0[int(0)], dist_coeffs_8);
                _S195 = !_S197;
            }
            if(_S195)
            {
                break;
            }
            float u_4 = proj_points_0[int(0)].x;
            float v_4 = proj_points_0[int(0)].y;
            float r2_4 = u_4 * u_4 + v_4 * v_4;
            float2  _S198 = proj_points_0[int(0)] * make_float2 (1.0f + r2_4 * ((*dist_coeffs_8)[int(0)] + r2_4 * ((*dist_coeffs_8)[int(1)] + r2_4 * ((*dist_coeffs_8)[int(2)] + r2_4 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_4 * v_4 + (*dist_coeffs_8)[int(5)] * (r2_4 + 2.0f * u_4 * u_4) + (*dist_coeffs_8)[int(6)] * r2_4, 2.0f * (*dist_coeffs_8)[int(5)] * u_4 * v_4 + (*dist_coeffs_8)[int(4)] * (r2_4 + 2.0f * v_4 * v_4) + (*dist_coeffs_8)[int(7)] * r2_4);
            float2  _S199 = _S198 + make_float2 ((*dist_coeffs_8)[int(8)] * _S198.x + (*dist_coeffs_8)[int(9)] * _S198.y, 0.0f);
            proj_points_0[int(0)] = make_float2 (fx_1 * _S199.x + cx_1, fy_1 * _S199.y + cy_1);
            break;
        }
        bool all_valid_0 = true & (!_S195);
        _S180 = &proj_points_0[int(1)];
        for(;;)
        {
            float _S200 = sigmas_0->p_0[int(1)].z;
            proj_points_0[int(1)] = float2 {sigmas_0->p_0[int(1)].x, sigmas_0->p_0[int(1)].y} / make_float2 (_S200);
            if(_S200 < 0.0f)
            {
                _S195 = true;
            }
            else
            {
                bool _S201 = is_valid_distortion(proj_points_0[int(1)], dist_coeffs_8);
                _S195 = !_S201;
            }
            if(_S195)
            {
                break;
            }
            float u_5 = proj_points_0[int(1)].x;
            float v_5 = proj_points_0[int(1)].y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S202 = proj_points_0[int(1)] * make_float2 (1.0f + r2_5 * ((*dist_coeffs_8)[int(0)] + r2_5 * ((*dist_coeffs_8)[int(1)] + r2_5 * ((*dist_coeffs_8)[int(2)] + r2_5 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_5 * v_5 + (*dist_coeffs_8)[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + (*dist_coeffs_8)[int(6)] * r2_5, 2.0f * (*dist_coeffs_8)[int(5)] * u_5 * v_5 + (*dist_coeffs_8)[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + (*dist_coeffs_8)[int(7)] * r2_5);
            float2  _S203 = _S202 + make_float2 ((*dist_coeffs_8)[int(8)] * _S202.x + (*dist_coeffs_8)[int(9)] * _S202.y, 0.0f);
            proj_points_0[int(1)] = make_float2 (fx_1 * _S203.x + cx_1, fy_1 * _S203.y + cy_1);
            break;
        }
        bool all_valid_1 = all_valid_0 & (!_S195);
        for(;;)
        {
            _S181 = &proj_points_0[int(2)];
            for(;;)
            {
                float _S204 = sigmas_0->p_0[int(2)].z;
                proj_points_0[int(2)] = float2 {sigmas_0->p_0[int(2)].x, sigmas_0->p_0[int(2)].y} / make_float2 (_S204);
                if(_S204 < 0.0f)
                {
                    _S195 = true;
                }
                else
                {
                    bool _S205 = is_valid_distortion(proj_points_0[int(2)], dist_coeffs_8);
                    _S195 = !_S205;
                }
                if(_S195)
                {
                    break;
                }
                float u_6 = proj_points_0[int(2)].x;
                float v_6 = proj_points_0[int(2)].y;
                float r2_6 = u_6 * u_6 + v_6 * v_6;
                float2  _S206 = proj_points_0[int(2)] * make_float2 (1.0f + r2_6 * ((*dist_coeffs_8)[int(0)] + r2_6 * ((*dist_coeffs_8)[int(1)] + r2_6 * ((*dist_coeffs_8)[int(2)] + r2_6 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_6 * v_6 + (*dist_coeffs_8)[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + (*dist_coeffs_8)[int(6)] * r2_6, 2.0f * (*dist_coeffs_8)[int(5)] * u_6 * v_6 + (*dist_coeffs_8)[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + (*dist_coeffs_8)[int(7)] * r2_6);
                float2  _S207 = _S206 + make_float2 ((*dist_coeffs_8)[int(8)] * _S206.x + (*dist_coeffs_8)[int(9)] * _S206.y, 0.0f);
                proj_points_0[int(2)] = make_float2 (fx_1 * _S207.x + cx_1, fy_1 * _S207.y + cy_1);
                break;
            }
            _S182 = all_valid_1 & (!_S195);
            break;
        }
        _S183 = &proj_points_0[int(3)];
        for(;;)
        {
            float _S208 = sigmas_0->p_0[int(3)].z;
            proj_points_0[int(3)] = float2 {sigmas_0->p_0[int(3)].x, sigmas_0->p_0[int(3)].y} / make_float2 (_S208);
            if(_S208 < 0.0f)
            {
                _S195 = true;
            }
            else
            {
                bool _S209 = is_valid_distortion(proj_points_0[int(3)], dist_coeffs_8);
                _S195 = !_S209;
            }
            if(_S195)
            {
                break;
            }
            float u_7 = proj_points_0[int(3)].x;
            float v_7 = proj_points_0[int(3)].y;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float2  _S210 = proj_points_0[int(3)] * make_float2 (1.0f + r2_7 * ((*dist_coeffs_8)[int(0)] + r2_7 * ((*dist_coeffs_8)[int(1)] + r2_7 * ((*dist_coeffs_8)[int(2)] + r2_7 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_7 * v_7 + (*dist_coeffs_8)[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + (*dist_coeffs_8)[int(6)] * r2_7, 2.0f * (*dist_coeffs_8)[int(5)] * u_7 * v_7 + (*dist_coeffs_8)[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + (*dist_coeffs_8)[int(7)] * r2_7);
            float2  _S211 = _S210 + make_float2 ((*dist_coeffs_8)[int(8)] * _S210.x + (*dist_coeffs_8)[int(9)] * _S210.y, 0.0f);
            proj_points_0[int(3)] = make_float2 (fx_1 * _S211.x + cx_1, fy_1 * _S211.y + cy_1);
            break;
        }
        bool all_valid_2 = _S182 & (!_S195);
        _S184 = &proj_points_0[int(4)];
        for(;;)
        {
            float _S212 = sigmas_0->p_0[int(4)].z;
            proj_points_0[int(4)] = float2 {sigmas_0->p_0[int(4)].x, sigmas_0->p_0[int(4)].y} / make_float2 (_S212);
            if(_S212 < 0.0f)
            {
                _S195 = true;
            }
            else
            {
                bool _S213 = is_valid_distortion(proj_points_0[int(4)], dist_coeffs_8);
                _S195 = !_S213;
            }
            if(_S195)
            {
                break;
            }
            float u_8 = proj_points_0[int(4)].x;
            float v_8 = proj_points_0[int(4)].y;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float2  _S214 = proj_points_0[int(4)] * make_float2 (1.0f + r2_8 * ((*dist_coeffs_8)[int(0)] + r2_8 * ((*dist_coeffs_8)[int(1)] + r2_8 * ((*dist_coeffs_8)[int(2)] + r2_8 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_8 * v_8 + (*dist_coeffs_8)[int(5)] * (r2_8 + 2.0f * u_8 * u_8) + (*dist_coeffs_8)[int(6)] * r2_8, 2.0f * (*dist_coeffs_8)[int(5)] * u_8 * v_8 + (*dist_coeffs_8)[int(4)] * (r2_8 + 2.0f * v_8 * v_8) + (*dist_coeffs_8)[int(7)] * r2_8);
            float2  _S215 = _S214 + make_float2 ((*dist_coeffs_8)[int(8)] * _S214.x + (*dist_coeffs_8)[int(9)] * _S214.y, 0.0f);
            proj_points_0[int(4)] = make_float2 (fx_1 * _S215.x + cx_1, fy_1 * _S215.y + cy_1);
            break;
        }
        bool all_valid_3 = all_valid_2 & (!_S195);
        for(;;)
        {
            _S185 = &proj_points_0[int(5)];
            for(;;)
            {
                float _S216 = sigmas_0->p_0[int(5)].z;
                proj_points_0[int(5)] = float2 {sigmas_0->p_0[int(5)].x, sigmas_0->p_0[int(5)].y} / make_float2 (_S216);
                if(_S216 < 0.0f)
                {
                    _S195 = true;
                }
                else
                {
                    bool _S217 = is_valid_distortion(proj_points_0[int(5)], dist_coeffs_8);
                    _S195 = !_S217;
                }
                if(_S195)
                {
                    break;
                }
                float u_9 = proj_points_0[int(5)].x;
                float v_9 = proj_points_0[int(5)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S218 = proj_points_0[int(5)] * make_float2 (1.0f + r2_9 * ((*dist_coeffs_8)[int(0)] + r2_9 * ((*dist_coeffs_8)[int(1)] + r2_9 * ((*dist_coeffs_8)[int(2)] + r2_9 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_9 * v_9 + (*dist_coeffs_8)[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + (*dist_coeffs_8)[int(6)] * r2_9, 2.0f * (*dist_coeffs_8)[int(5)] * u_9 * v_9 + (*dist_coeffs_8)[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + (*dist_coeffs_8)[int(7)] * r2_9);
                float2  _S219 = _S218 + make_float2 ((*dist_coeffs_8)[int(8)] * _S218.x + (*dist_coeffs_8)[int(9)] * _S218.y, 0.0f);
                proj_points_0[int(5)] = make_float2 (fx_1 * _S219.x + cx_1, fy_1 * _S219.y + cy_1);
                break;
            }
            _S186 = all_valid_3 & (!_S195);
            break;
        }
        _S187 = &proj_points_0[int(6)];
        for(;;)
        {
            float _S220 = sigmas_0->p_0[int(6)].z;
            proj_points_0[int(6)] = float2 {sigmas_0->p_0[int(6)].x, sigmas_0->p_0[int(6)].y} / make_float2 (_S220);
            if(_S220 < 0.0f)
            {
                _S195 = true;
            }
            else
            {
                bool _S221 = is_valid_distortion(proj_points_0[int(6)], dist_coeffs_8);
                _S195 = !_S221;
            }
            if(_S195)
            {
                break;
            }
            float u_10 = proj_points_0[int(6)].x;
            float v_10 = proj_points_0[int(6)].y;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float2  _S222 = proj_points_0[int(6)] * make_float2 (1.0f + r2_10 * ((*dist_coeffs_8)[int(0)] + r2_10 * ((*dist_coeffs_8)[int(1)] + r2_10 * ((*dist_coeffs_8)[int(2)] + r2_10 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_10 * v_10 + (*dist_coeffs_8)[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + (*dist_coeffs_8)[int(6)] * r2_10, 2.0f * (*dist_coeffs_8)[int(5)] * u_10 * v_10 + (*dist_coeffs_8)[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + (*dist_coeffs_8)[int(7)] * r2_10);
            float2  _S223 = _S222 + make_float2 ((*dist_coeffs_8)[int(8)] * _S222.x + (*dist_coeffs_8)[int(9)] * _S222.y, 0.0f);
            proj_points_0[int(6)] = make_float2 (fx_1 * _S223.x + cx_1, fy_1 * _S223.y + cy_1);
            break;
        }
        _S188 = _S186 & (!_S195);
        break;
    }
    if(!_S188)
    {
        return false;
    }
    float2  _S224 = *mean2d_1 + make_float2 (sigmas_0->w_mean_0[int(0)]) * *_S179;
    *mean2d_1 = _S224;
    float2  _S225 = _S224 + make_float2 (sigmas_0->w_mean_0[int(1)]) * *_S180;
    *mean2d_1 = _S225;
    float2  _S226 = _S225 + make_float2 (sigmas_0->w_mean_0[int(2)]) * *_S181;
    *mean2d_1 = _S226;
    float2  _S227 = _S226 + make_float2 (sigmas_0->w_mean_0[int(3)]) * *_S183;
    *mean2d_1 = _S227;
    float2  _S228 = _S227 + make_float2 (sigmas_0->w_mean_0[int(4)]) * *_S184;
    *mean2d_1 = _S228;
    float2  _S229 = _S228 + make_float2 (sigmas_0->w_mean_0[int(5)]) * *_S185;
    *mean2d_1 = _S229;
    float2  _S230 = _S229 + make_float2 (sigmas_0->w_mean_0[int(6)]) * *_S187;
    *mean2d_1 = _S230;
    float _S231 = - lim_x_neg_0;
    float _S232 = - lim_y_neg_0;
    float2  _S233 = make_float2 (clamp_0(_S230.x, _S231, lim_x_pos_0), clamp_0(_S230.y, _S232, lim_y_pos_0));
    float2  d_0 = make_float2 (clamp_0((*_S179).x, _S231, lim_x_pos_0), clamp_0((*_S179).y, _S232, lim_y_pos_0)) - _S233;
    float _S234 = d_0.x;
    float _S235 = d_0.y;
    float _S236 = _S234 * _S235;
    Matrix<float, 2, 2>  _S237 = *cov2d_1 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S234 * _S234, _S236, _S236, _S235 * _S235);
    *cov2d_1 = _S237;
    float2  d_1 = make_float2 (clamp_0((*_S180).x, _S231, lim_x_pos_0), clamp_0((*_S180).y, _S232, lim_y_pos_0)) - _S233;
    float _S238 = d_1.x;
    float _S239 = d_1.y;
    float _S240 = _S238 * _S239;
    Matrix<float, 2, 2>  _S241 = _S237 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S238 * _S238, _S240, _S240, _S239 * _S239);
    *cov2d_1 = _S241;
    float2  d_2 = make_float2 (clamp_0((*_S181).x, _S231, lim_x_pos_0), clamp_0((*_S181).y, _S232, lim_y_pos_0)) - _S233;
    float _S242 = d_2.x;
    float _S243 = d_2.y;
    float _S244 = _S242 * _S243;
    Matrix<float, 2, 2>  _S245 = _S241 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S242 * _S242, _S244, _S244, _S243 * _S243);
    *cov2d_1 = _S245;
    float2  d_3 = make_float2 (clamp_0((*_S183).x, _S231, lim_x_pos_0), clamp_0((*_S183).y, _S232, lim_y_pos_0)) - _S233;
    float _S246 = d_3.x;
    float _S247 = d_3.y;
    float _S248 = _S246 * _S247;
    Matrix<float, 2, 2>  _S249 = _S245 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S246 * _S246, _S248, _S248, _S247 * _S247);
    *cov2d_1 = _S249;
    float2  d_4 = make_float2 (clamp_0((*_S184).x, _S231, lim_x_pos_0), clamp_0((*_S184).y, _S232, lim_y_pos_0)) - _S233;
    float _S250 = d_4.x;
    float _S251 = d_4.y;
    float _S252 = _S250 * _S251;
    Matrix<float, 2, 2>  _S253 = _S249 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S250 * _S250, _S252, _S252, _S251 * _S251);
    *cov2d_1 = _S253;
    float2  d_5 = make_float2 (clamp_0((*_S185).x, _S231, lim_x_pos_0), clamp_0((*_S185).y, _S232, lim_y_pos_0)) - _S233;
    float _S254 = d_5.x;
    float _S255 = d_5.y;
    float _S256 = _S254 * _S255;
    Matrix<float, 2, 2>  _S257 = _S253 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S254 * _S254, _S256, _S256, _S255 * _S255);
    *cov2d_1 = _S257;
    float2  d_6 = make_float2 (clamp_0((*_S187).x, _S231, lim_x_pos_0), clamp_0((*_S187).y, _S232, lim_y_pos_0)) - _S233;
    float _S258 = d_6.x;
    float _S259 = d_6.y;
    float _S260 = _S258 * _S259;
    *cov2d_1 = _S257 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S258 * _S258, _S260, _S260, _S259 * _S259);
    return true;
}

inline __device__ bool persp_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_9, uint width_2, uint height_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    int2  _S261 = make_int2 (int(0));
    float2  _S262 = make_float2 ((float)_S261.x, (float)_S261.y);
    *mean2d_2 = _S262;
    *cov2d_2 = makeMatrix<float, 2, 2> (0.0f);
    float fx_2 = intrins_1.x;
    float fy_2 = intrins_1.y;
    float cx_2 = intrins_1.z;
    float cy_2 = intrins_1.w;
    float _S263 = float(width_2);
    float _S264 = float(height_2);
    float _S265 = 0.30000001192092896f * (0.5f * _S263 / fx_2) * fx_2;
    float lim_x_pos_1 = _S263 + _S265;
    float lim_x_neg_1 = cx_2 + _S265;
    float _S266 = 0.30000001192092896f * (0.5f * _S264 / fy_2) * fy_2;
    float lim_y_pos_1 = _S264 + _S266;
    float lim_y_neg_1 = cy_2 + _S266;
    FixedArray<float2 , 7>  proj_points_1;
    float2  _S267 = float2 {sigmas_1->p_0[int(0)].x, sigmas_1->p_0[int(0)].y} / make_float2 (sigmas_1->p_0[int(0)].z);
    float u_11 = _S267.x;
    float v_11 = _S267.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float _S268 = 2.0f * (*dist_coeffs_9)[int(4)];
    float _S269 = 2.0f * (*dist_coeffs_9)[int(5)];
    float2  _S270 = _S267 * make_float2 (1.0f + r2_11 * ((*dist_coeffs_9)[int(0)] + r2_11 * ((*dist_coeffs_9)[int(1)] + r2_11 * ((*dist_coeffs_9)[int(2)] + r2_11 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_11 * v_11 + (*dist_coeffs_9)[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + (*dist_coeffs_9)[int(6)] * r2_11, _S269 * u_11 * v_11 + (*dist_coeffs_9)[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + (*dist_coeffs_9)[int(7)] * r2_11);
    float2  _S271 = _S270 + make_float2 ((*dist_coeffs_9)[int(8)] * _S270.x + (*dist_coeffs_9)[int(9)] * _S270.y, 0.0f);
    float _S272 = fx_2 * _S271.x + cx_2;
    float _S273 = fy_2 * _S271.y + cy_2;
    float2  _S274 = make_float2 (_S272, _S273);
    proj_points_1[int(0)] = _S274;
    float2  _S275 = float2 {sigmas_1->p_0[int(1)].x, sigmas_1->p_0[int(1)].y} / make_float2 (sigmas_1->p_0[int(1)].z);
    float u_12 = _S275.x;
    float v_12 = _S275.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float2  _S276 = _S275 * make_float2 (1.0f + r2_12 * ((*dist_coeffs_9)[int(0)] + r2_12 * ((*dist_coeffs_9)[int(1)] + r2_12 * ((*dist_coeffs_9)[int(2)] + r2_12 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_12 * v_12 + (*dist_coeffs_9)[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + (*dist_coeffs_9)[int(6)] * r2_12, _S269 * u_12 * v_12 + (*dist_coeffs_9)[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + (*dist_coeffs_9)[int(7)] * r2_12);
    float2  _S277 = _S276 + make_float2 ((*dist_coeffs_9)[int(8)] * _S276.x + (*dist_coeffs_9)[int(9)] * _S276.y, 0.0f);
    float _S278 = fx_2 * _S277.x + cx_2;
    float _S279 = fy_2 * _S277.y + cy_2;
    float2  _S280 = make_float2 (_S278, _S279);
    proj_points_1[int(1)] = _S280;
    float2  _S281 = float2 {sigmas_1->p_0[int(2)].x, sigmas_1->p_0[int(2)].y} / make_float2 (sigmas_1->p_0[int(2)].z);
    float u_13 = _S281.x;
    float v_13 = _S281.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S282 = _S281 * make_float2 (1.0f + r2_13 * ((*dist_coeffs_9)[int(0)] + r2_13 * ((*dist_coeffs_9)[int(1)] + r2_13 * ((*dist_coeffs_9)[int(2)] + r2_13 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_13 * v_13 + (*dist_coeffs_9)[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + (*dist_coeffs_9)[int(6)] * r2_13, _S269 * u_13 * v_13 + (*dist_coeffs_9)[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + (*dist_coeffs_9)[int(7)] * r2_13);
    float2  _S283 = _S282 + make_float2 ((*dist_coeffs_9)[int(8)] * _S282.x + (*dist_coeffs_9)[int(9)] * _S282.y, 0.0f);
    float _S284 = fx_2 * _S283.x + cx_2;
    float _S285 = fy_2 * _S283.y + cy_2;
    float2  _S286 = make_float2 (_S284, _S285);
    proj_points_1[int(2)] = _S286;
    float2  _S287 = float2 {sigmas_1->p_0[int(3)].x, sigmas_1->p_0[int(3)].y} / make_float2 (sigmas_1->p_0[int(3)].z);
    float u_14 = _S287.x;
    float v_14 = _S287.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S288 = _S287 * make_float2 (1.0f + r2_14 * ((*dist_coeffs_9)[int(0)] + r2_14 * ((*dist_coeffs_9)[int(1)] + r2_14 * ((*dist_coeffs_9)[int(2)] + r2_14 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_14 * v_14 + (*dist_coeffs_9)[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + (*dist_coeffs_9)[int(6)] * r2_14, _S269 * u_14 * v_14 + (*dist_coeffs_9)[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + (*dist_coeffs_9)[int(7)] * r2_14);
    float2  _S289 = _S288 + make_float2 ((*dist_coeffs_9)[int(8)] * _S288.x + (*dist_coeffs_9)[int(9)] * _S288.y, 0.0f);
    float _S290 = fx_2 * _S289.x + cx_2;
    float _S291 = fy_2 * _S289.y + cy_2;
    float2  _S292 = make_float2 (_S290, _S291);
    proj_points_1[int(3)] = _S292;
    float2  _S293 = float2 {sigmas_1->p_0[int(4)].x, sigmas_1->p_0[int(4)].y} / make_float2 (sigmas_1->p_0[int(4)].z);
    float u_15 = _S293.x;
    float v_15 = _S293.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S294 = _S293 * make_float2 (1.0f + r2_15 * ((*dist_coeffs_9)[int(0)] + r2_15 * ((*dist_coeffs_9)[int(1)] + r2_15 * ((*dist_coeffs_9)[int(2)] + r2_15 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_15 * v_15 + (*dist_coeffs_9)[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + (*dist_coeffs_9)[int(6)] * r2_15, _S269 * u_15 * v_15 + (*dist_coeffs_9)[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + (*dist_coeffs_9)[int(7)] * r2_15);
    float2  _S295 = _S294 + make_float2 ((*dist_coeffs_9)[int(8)] * _S294.x + (*dist_coeffs_9)[int(9)] * _S294.y, 0.0f);
    float _S296 = fx_2 * _S295.x + cx_2;
    float _S297 = fy_2 * _S295.y + cy_2;
    float2  _S298 = make_float2 (_S296, _S297);
    proj_points_1[int(4)] = _S298;
    float2  _S299 = float2 {sigmas_1->p_0[int(5)].x, sigmas_1->p_0[int(5)].y} / make_float2 (sigmas_1->p_0[int(5)].z);
    float u_16 = _S299.x;
    float v_16 = _S299.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S300 = _S299 * make_float2 (1.0f + r2_16 * ((*dist_coeffs_9)[int(0)] + r2_16 * ((*dist_coeffs_9)[int(1)] + r2_16 * ((*dist_coeffs_9)[int(2)] + r2_16 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_16 * v_16 + (*dist_coeffs_9)[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + (*dist_coeffs_9)[int(6)] * r2_16, _S269 * u_16 * v_16 + (*dist_coeffs_9)[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + (*dist_coeffs_9)[int(7)] * r2_16);
    float2  _S301 = _S300 + make_float2 ((*dist_coeffs_9)[int(8)] * _S300.x + (*dist_coeffs_9)[int(9)] * _S300.y, 0.0f);
    float _S302 = fx_2 * _S301.x + cx_2;
    float _S303 = fy_2 * _S301.y + cy_2;
    float2  _S304 = make_float2 (_S302, _S303);
    proj_points_1[int(5)] = _S304;
    float2  _S305 = float2 {sigmas_1->p_0[int(6)].x, sigmas_1->p_0[int(6)].y} / make_float2 (sigmas_1->p_0[int(6)].z);
    float u_17 = _S305.x;
    float v_17 = _S305.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float2  _S306 = _S305 * make_float2 (1.0f + r2_17 * ((*dist_coeffs_9)[int(0)] + r2_17 * ((*dist_coeffs_9)[int(1)] + r2_17 * ((*dist_coeffs_9)[int(2)] + r2_17 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S268 * u_17 * v_17 + (*dist_coeffs_9)[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + (*dist_coeffs_9)[int(6)] * r2_17, _S269 * u_17 * v_17 + (*dist_coeffs_9)[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + (*dist_coeffs_9)[int(7)] * r2_17);
    float2  _S307 = _S306 + make_float2 ((*dist_coeffs_9)[int(8)] * _S306.x + (*dist_coeffs_9)[int(9)] * _S306.y, 0.0f);
    float _S308 = fx_2 * _S307.x + cx_2;
    float _S309 = fy_2 * _S307.y + cy_2;
    float2  _S310 = make_float2 (_S308, _S309);
    proj_points_1[int(6)] = _S310;
    float2  _S311 = *mean2d_2 + make_float2 (sigmas_1->w_mean_0[int(0)]) * _S274;
    *mean2d_2 = _S311;
    float2  _S312 = _S311 + make_float2 (sigmas_1->w_mean_0[int(1)]) * _S280;
    *mean2d_2 = _S312;
    float2  _S313 = _S312 + make_float2 (sigmas_1->w_mean_0[int(2)]) * _S286;
    *mean2d_2 = _S313;
    float2  _S314 = _S313 + make_float2 (sigmas_1->w_mean_0[int(3)]) * _S292;
    *mean2d_2 = _S314;
    float2  _S315 = _S314 + make_float2 (sigmas_1->w_mean_0[int(4)]) * _S298;
    *mean2d_2 = _S315;
    float2  _S316 = _S315 + make_float2 (sigmas_1->w_mean_0[int(5)]) * _S304;
    *mean2d_2 = _S316;
    float2  _S317 = _S316 + make_float2 (sigmas_1->w_mean_0[int(6)]) * _S310;
    *mean2d_2 = _S317;
    float _S318 = - lim_x_neg_1;
    float _S319 = - lim_y_neg_1;
    float2  _S320 = make_float2 (clamp_0(_S317.x, _S318, lim_x_pos_1), clamp_0(_S317.y, _S319, lim_y_pos_1));
    float2  d_7 = make_float2 (clamp_0(_S272, _S318, lim_x_pos_1), clamp_0(_S273, _S319, lim_y_pos_1)) - _S320;
    float _S321 = d_7.x;
    float _S322 = d_7.y;
    float _S323 = _S321 * _S322;
    Matrix<float, 2, 2>  _S324 = *cov2d_2 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S321 * _S321, _S323, _S323, _S322 * _S322);
    *cov2d_2 = _S324;
    float2  d_8 = make_float2 (clamp_0(_S278, _S318, lim_x_pos_1), clamp_0(_S279, _S319, lim_y_pos_1)) - _S320;
    float _S325 = d_8.x;
    float _S326 = d_8.y;
    float _S327 = _S325 * _S326;
    Matrix<float, 2, 2>  _S328 = _S324 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S325 * _S325, _S327, _S327, _S326 * _S326);
    *cov2d_2 = _S328;
    float2  d_9 = make_float2 (clamp_0(_S284, _S318, lim_x_pos_1), clamp_0(_S285, _S319, lim_y_pos_1)) - _S320;
    float _S329 = d_9.x;
    float _S330 = d_9.y;
    float _S331 = _S329 * _S330;
    Matrix<float, 2, 2>  _S332 = _S328 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S329 * _S329, _S331, _S331, _S330 * _S330);
    *cov2d_2 = _S332;
    float2  d_10 = make_float2 (clamp_0(_S290, _S318, lim_x_pos_1), clamp_0(_S291, _S319, lim_y_pos_1)) - _S320;
    float _S333 = d_10.x;
    float _S334 = d_10.y;
    float _S335 = _S333 * _S334;
    Matrix<float, 2, 2>  _S336 = _S332 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S333 * _S333, _S335, _S335, _S334 * _S334);
    *cov2d_2 = _S336;
    float2  d_11 = make_float2 (clamp_0(_S296, _S318, lim_x_pos_1), clamp_0(_S297, _S319, lim_y_pos_1)) - _S320;
    float _S337 = d_11.x;
    float _S338 = d_11.y;
    float _S339 = _S337 * _S338;
    Matrix<float, 2, 2>  _S340 = _S336 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S337 * _S337, _S339, _S339, _S338 * _S338);
    *cov2d_2 = _S340;
    float2  d_12 = make_float2 (clamp_0(_S302, _S318, lim_x_pos_1), clamp_0(_S303, _S319, lim_y_pos_1)) - _S320;
    float _S341 = d_12.x;
    float _S342 = d_12.y;
    float _S343 = _S341 * _S342;
    Matrix<float, 2, 2>  _S344 = _S340 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S341 * _S341, _S343, _S343, _S342 * _S342);
    *cov2d_2 = _S344;
    float2  d_13 = make_float2 (clamp_0(_S308, _S318, lim_x_pos_1), clamp_0(_S309, _S319, lim_y_pos_1)) - _S320;
    float _S345 = d_13.x;
    float _S346 = d_13.y;
    float _S347 = _S345 * _S346;
    *cov2d_2 = _S344 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S345 * _S345, _S347, _S347, _S346 * _S346);
    return true;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_8, float dOut_11)
{
    DiffPair_float_0 _S348 = *dpx_8;
    float _S349 = - (*dpy_4).primal_0 / ((*dpx_8).primal_0 * (*dpx_8).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S349;
    float _S350 = _S348.primal_0 / (_S348.primal_0 * _S348.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S350;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S351, float _S352)
{
    return (F32_atan2((_S351), (_S352)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S353, DiffPair_float_0 * _S354, float _S355)
{
    _d_atan2_0(_S353, _S354, _S355);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_9, float _s_dOut_2)
{
    float _S356 = (*dpx_9).primal_0.x;
    float _S357 = (*dpx_9).primal_0.y;
    DiffPair_float_0 _S358;
    (&_S358)->primal_0 = _S356 * _S356 + _S357 * _S357;
    (&_S358)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S358, _s_dOut_2);
    float _S359 = (*dpx_9).primal_0.y * _S358.differential_0;
    float _S360 = _S359 + _S359;
    float _S361 = (*dpx_9).primal_0.x * _S358.differential_0;
    float _S362 = _S361 + _S361;
    float2  _S363 = make_float2 (0.0f);
    *&((&_S363)->y) = _S360;
    *&((&_S363)->x) = _S362;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S363;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S364, float _S365)
{
    s_bwd_prop_length_impl_1(_S364, _S365);
    return;
}

inline __device__ bool fisheye_proj_3dgs_0(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_10, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    float k_0;
    float2  _S366;
    float _S367;
    float _S368;
    bool _S369;
    for(;;)
    {
        float2  _S370 = float2 {mean3d_1.x, mean3d_1.y};
        _S366 = _S370;
        float r_8 = length_0(_S370);
        _S367 = r_8;
        float _S371 = mean3d_1.z;
        _S368 = _S371;
        float theta_3 = (F32_atan2((r_8), (_S371)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S371;
        }
        else
        {
            k_0 = theta_3 / r_8;
        }
        float2  _S372 = _S370 * make_float2 (k_0);
        *mean2d_3 = _S372;
        bool _S373 = is_valid_distortion(_S372, dist_coeffs_10);
        bool _S374 = !_S373;
        _S369 = _S374;
        if(_S374)
        {
            break;
        }
        float u_18 = (*mean2d_3).x;
        float v_18 = (*mean2d_3).y;
        float r2_18 = u_18 * u_18 + v_18 * v_18;
        float2  _S375 = *mean2d_3 * make_float2 (1.0f + r2_18 * ((*dist_coeffs_10)[int(0)] + r2_18 * ((*dist_coeffs_10)[int(1)] + r2_18 * ((*dist_coeffs_10)[int(2)] + r2_18 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_18 * v_18 + (*dist_coeffs_10)[int(5)] * (r2_18 + 2.0f * u_18 * u_18) + (*dist_coeffs_10)[int(6)] * r2_18, 2.0f * (*dist_coeffs_10)[int(5)] * u_18 * v_18 + (*dist_coeffs_10)[int(4)] * (r2_18 + 2.0f * v_18 * v_18) + (*dist_coeffs_10)[int(7)] * r2_18);
        float2  _S376 = _S375 + make_float2 ((*dist_coeffs_10)[int(8)] * _S375.x + (*dist_coeffs_10)[int(9)] * _S375.y, 0.0f);
        *mean2d_3 = make_float2 (intrins_2.x * _S376.x + intrins_2.z, intrins_2.y * _S376.y + intrins_2.w);
        break;
    }
    if(!!_S369)
    {
        return false;
    }
    Matrix<float, 2, 3>  J_3;
    float2  _S377 = make_float2 (0.0f);
    float2  seed_4 = _S377;
    *&((&seed_4)->x) = 1.0f;
    float2  _S378 = seed_4;
    float _S379 = s_primal_ctx_atan2_0(_S367, _S368);
    bool _S380 = _S379 < 0.00100000004749745f;
    float _S381;
    float _S382;
    float _S383;
    if(_S380)
    {
        float _S384 = 1.0f - _S379 * _S379 / 3.0f;
        float _S385 = _S368 * _S368;
        k_0 = _S384 / _S368;
        _S381 = 0.0f;
        _S382 = _S385;
        _S383 = _S384;
    }
    else
    {
        float _S386 = _S367 * _S367;
        k_0 = _S379 / _S367;
        _S381 = _S386;
        _S382 = 0.0f;
        _S383 = 0.0f;
    }
    float2  _S387 = make_float2 (k_0);
    float2  _S388 = _S366 * make_float2 (k_0);
    float u_19 = _S388.x;
    float v_19 = _S388.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S389 = (*dist_coeffs_10)[int(2)] + r2_19 * (*dist_coeffs_10)[int(3)];
    float _S390 = (*dist_coeffs_10)[int(1)] + r2_19 * _S389;
    float _S391 = (*dist_coeffs_10)[int(0)] + r2_19 * _S390;
    float _S392 = 2.0f * (*dist_coeffs_10)[int(4)];
    float _S393 = 2.0f * (*dist_coeffs_10)[int(5)];
    float fx_3 = intrins_2.x;
    float fy_3 = intrins_2.y;
    float _S394 = fx_3 * _S378.x;
    float2  _S395 = make_float2 (_S394, fy_3 * _S378.y) + make_float2 ((*dist_coeffs_10)[int(8)] * _S394, (*dist_coeffs_10)[int(9)] * _S394);
    float2  _S396 = _S388 * _S395;
    float _S397 = (*dist_coeffs_10)[int(4)] * _S395.y;
    float _S398 = (*dist_coeffs_10)[int(5)] * _S395.x;
    float _S399 = _S396.x + _S396.y;
    float _S400 = r2_19 * _S399;
    float _S401 = r2_19 * _S400;
    float _S402 = (*dist_coeffs_10)[int(7)] * _S395.y + _S397 + (*dist_coeffs_10)[int(6)] * _S395.x + _S398 + _S391 * _S399 + _S390 * _S400 + _S389 * _S401 + (*dist_coeffs_10)[int(3)] * (r2_19 * _S401);
    float _S403 = v_19 * _S402;
    float _S404 = u_19 * _S402;
    float2  _S405 = make_float2 (1.0f + r2_19 * _S391) * _S395 + make_float2 (_S393 * (v_19 * _S395.y) + 2.0f * u_19 * _S398 + 2.0f * (u_19 * _S398) + _S392 * (v_19 * _S395.x) + _S404 + _S404, 2.0f * v_19 * _S397 + 2.0f * (v_19 * _S397) + _S393 * u_19 * _S395.y + _S392 * u_19 * _S395.x + _S403 + _S403);
    float2  _S406 = _S366 * _S405;
    float2  _S407 = _S387 * _S405;
    float _S408 = _S406.x + _S406.y;
    if(_S380)
    {
        float _S409 = _S408 / _S382;
        float _S410 = _S383 * - _S409;
        float _S411 = _S379 * (0.3333333432674408f * - (_S368 * _S409));
        k_0 = _S411 + _S411;
        _S381 = _S410;
        _S382 = 0.0f;
    }
    else
    {
        float _S412 = _S408 / _S381;
        float _S413 = _S379 * - _S412;
        k_0 = _S367 * _S412;
        _S381 = 0.0f;
        _S382 = _S413;
    }
    float _S414 = _S367;
    DiffPair_float_0 _S415;
    (&_S415)->primal_0 = _S367;
    (&_S415)->differential_0 = 0.0f;
    float _S416 = _S368;
    DiffPair_float_0 _S417;
    (&_S417)->primal_0 = _S368;
    (&_S417)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S415, &_S417, k_0);
    float _S418 = _S417.differential_0 + _S381;
    float _S419 = _S415.differential_0 + _S382;
    float2  _S420 = make_float2 (0.0f);
    float2  _S421 = _S366;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S422;
    (&_S422)->primal_0 = _S366;
    (&_S422)->differential_0 = _S420;
    s_bwd_length_impl_1(&_S422, _S419);
    float2  _S423 = _S422.differential_0 + _S407;
    float3  _S424 = make_float3 (_S423.x, _S423.y, _S418);
    J_3[int(0)] = _S424;
    float2  seed_5 = _S377;
    *&((&seed_5)->y) = 1.0f;
    float2  _S425 = seed_5;
    if(_S380)
    {
        float _S426 = 1.0f - _S379 * _S379 / 3.0f;
        float _S427 = _S368 * _S368;
        k_0 = _S426 / _S368;
        _S381 = 0.0f;
        _S382 = _S427;
        _S383 = _S426;
    }
    else
    {
        float _S428 = _S367 * _S367;
        k_0 = _S379 / _S367;
        _S381 = _S428;
        _S382 = 0.0f;
        _S383 = 0.0f;
    }
    float2  _S429 = make_float2 (k_0);
    float2  _S430 = _S366 * make_float2 (k_0);
    float u_20 = _S430.x;
    float v_20 = _S430.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S431 = (*dist_coeffs_10)[int(2)] + r2_20 * (*dist_coeffs_10)[int(3)];
    float _S432 = (*dist_coeffs_10)[int(1)] + r2_20 * _S431;
    float _S433 = (*dist_coeffs_10)[int(0)] + r2_20 * _S432;
    float _S434 = fx_3 * _S425.x;
    float2  _S435 = make_float2 (_S434, fy_3 * _S425.y) + make_float2 ((*dist_coeffs_10)[int(8)] * _S434, (*dist_coeffs_10)[int(9)] * _S434);
    float2  _S436 = _S430 * _S435;
    float _S437 = (*dist_coeffs_10)[int(4)] * _S435.y;
    float _S438 = (*dist_coeffs_10)[int(5)] * _S435.x;
    float _S439 = _S436.x + _S436.y;
    float _S440 = r2_20 * _S439;
    float _S441 = r2_20 * _S440;
    float _S442 = (*dist_coeffs_10)[int(7)] * _S435.y + _S437 + (*dist_coeffs_10)[int(6)] * _S435.x + _S438 + _S433 * _S439 + _S432 * _S440 + _S431 * _S441 + (*dist_coeffs_10)[int(3)] * (r2_20 * _S441);
    float _S443 = v_20 * _S442;
    float _S444 = u_20 * _S442;
    float2  _S445 = make_float2 (1.0f + r2_20 * _S433) * _S435 + make_float2 (_S393 * (v_20 * _S435.y) + 2.0f * u_20 * _S438 + 2.0f * (u_20 * _S438) + _S392 * (v_20 * _S435.x) + _S444 + _S444, 2.0f * v_20 * _S437 + 2.0f * (v_20 * _S437) + _S393 * u_20 * _S435.y + _S392 * u_20 * _S435.x + _S443 + _S443);
    float2  _S446 = _S366 * _S445;
    float2  _S447 = _S429 * _S445;
    float _S448 = _S446.x + _S446.y;
    if(_S380)
    {
        float _S449 = _S448 / _S382;
        float _S450 = _S383 * - _S449;
        float _S451 = _S379 * (0.3333333432674408f * - (_S368 * _S449));
        k_0 = _S451 + _S451;
        _S381 = _S450;
        _S382 = 0.0f;
    }
    else
    {
        float _S452 = _S448 / _S381;
        float _S453 = _S379 * - _S452;
        k_0 = _S367 * _S452;
        _S381 = 0.0f;
        _S382 = _S453;
    }
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S414;
    (&_S454)->differential_0 = 0.0f;
    DiffPair_float_0 _S455;
    (&_S455)->primal_0 = _S416;
    (&_S455)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S454, &_S455, k_0);
    float _S456 = _S455.differential_0 + _S381;
    float _S457 = _S454.differential_0 + _S382;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S458;
    (&_S458)->primal_0 = _S421;
    (&_S458)->differential_0 = _S420;
    s_bwd_length_impl_1(&_S458, _S457);
    float2  _S459 = _S458.differential_0 + _S447;
    float3  _S460 = make_float3 (_S459.x, _S459.y, _S456);
    J_3[int(1)] = _S460;
    *cov2d_3 = mul_6(mul_5(J_3, cov3d_1), transpose_1(J_3));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_1(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_11, Matrix<float, 2, 2>  * cov2d_4, float2  * mean2d_4)
{
    float2  _S461 = float2 {mean3d_2.x, mean3d_2.y};
    float r_9 = length_0(_S461);
    float _S462 = mean3d_2.z;
    float theta_4 = (F32_atan2((r_9), (_S462)));
    float k_1;
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S462;
    }
    else
    {
        k_1 = theta_4 / r_9;
    }
    float2  _S463 = _S461 * make_float2 (k_1);
    float u_21 = _S463.x;
    float v_21 = _S463.y;
    float r2_21 = u_21 * u_21 + v_21 * v_21;
    float _S464 = 2.0f * (*dist_coeffs_11)[int(4)];
    float _S465 = 2.0f * (*dist_coeffs_11)[int(5)];
    float2  _S466 = _S463 * make_float2 (1.0f + r2_21 * ((*dist_coeffs_11)[int(0)] + r2_21 * ((*dist_coeffs_11)[int(1)] + r2_21 * ((*dist_coeffs_11)[int(2)] + r2_21 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S464 * u_21 * v_21 + (*dist_coeffs_11)[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + (*dist_coeffs_11)[int(6)] * r2_21, _S465 * u_21 * v_21 + (*dist_coeffs_11)[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + (*dist_coeffs_11)[int(7)] * r2_21);
    float2  _S467 = _S466 + make_float2 ((*dist_coeffs_11)[int(8)] * _S466.x + (*dist_coeffs_11)[int(9)] * _S466.y, 0.0f);
    float fx_4 = intrins_3.x;
    float fy_4 = intrins_3.y;
    *mean2d_4 = make_float2 (fx_4 * _S467.x + intrins_3.z, fy_4 * _S467.y + intrins_3.w);
    Matrix<float, 2, 3>  J_4;
    float2  _S468 = make_float2 (0.0f);
    float2  seed_6 = _S468;
    *&((&seed_6)->x) = 1.0f;
    float2  _S469 = seed_6;
    float _S470 = s_primal_ctx_atan2_0(r_9, _S462);
    bool _S471 = _S470 < 0.00100000004749745f;
    float _S472;
    float _S473;
    float _S474;
    if(_S471)
    {
        float _S475 = 1.0f - _S470 * _S470 / 3.0f;
        float _S476 = _S462 * _S462;
        k_1 = _S475 / _S462;
        _S472 = 0.0f;
        _S473 = _S476;
        _S474 = _S475;
    }
    else
    {
        float _S477 = r_9 * r_9;
        k_1 = _S470 / r_9;
        _S472 = _S477;
        _S473 = 0.0f;
        _S474 = 0.0f;
    }
    float2  _S478 = make_float2 (k_1);
    float2  _S479 = _S461 * make_float2 (k_1);
    float u_22 = _S479.x;
    float v_22 = _S479.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S480 = (*dist_coeffs_11)[int(2)] + r2_22 * (*dist_coeffs_11)[int(3)];
    float _S481 = (*dist_coeffs_11)[int(1)] + r2_22 * _S480;
    float _S482 = (*dist_coeffs_11)[int(0)] + r2_22 * _S481;
    float _S483 = fx_4 * _S469.x;
    float2  _S484 = make_float2 (_S483, fy_4 * _S469.y) + make_float2 ((*dist_coeffs_11)[int(8)] * _S483, (*dist_coeffs_11)[int(9)] * _S483);
    float2  _S485 = _S479 * _S484;
    float _S486 = (*dist_coeffs_11)[int(4)] * _S484.y;
    float _S487 = (*dist_coeffs_11)[int(5)] * _S484.x;
    float _S488 = _S485.x + _S485.y;
    float _S489 = r2_22 * _S488;
    float _S490 = r2_22 * _S489;
    float _S491 = (*dist_coeffs_11)[int(7)] * _S484.y + _S486 + (*dist_coeffs_11)[int(6)] * _S484.x + _S487 + _S482 * _S488 + _S481 * _S489 + _S480 * _S490 + (*dist_coeffs_11)[int(3)] * (r2_22 * _S490);
    float _S492 = v_22 * _S491;
    float _S493 = u_22 * _S491;
    float2  _S494 = make_float2 (1.0f + r2_22 * _S482) * _S484 + make_float2 (_S465 * (v_22 * _S484.y) + 2.0f * u_22 * _S487 + 2.0f * (u_22 * _S487) + _S464 * (v_22 * _S484.x) + _S493 + _S493, 2.0f * v_22 * _S486 + 2.0f * (v_22 * _S486) + _S465 * u_22 * _S484.y + _S464 * u_22 * _S484.x + _S492 + _S492);
    float2  _S495 = _S461 * _S494;
    float2  _S496 = _S478 * _S494;
    float _S497 = _S495.x + _S495.y;
    if(_S471)
    {
        float _S498 = _S497 / _S473;
        float _S499 = _S474 * - _S498;
        float _S500 = _S470 * (0.3333333432674408f * - (_S462 * _S498));
        k_1 = _S500 + _S500;
        _S472 = _S499;
        _S473 = 0.0f;
    }
    else
    {
        float _S501 = _S497 / _S472;
        float _S502 = _S470 * - _S501;
        k_1 = r_9 * _S501;
        _S472 = 0.0f;
        _S473 = _S502;
    }
    DiffPair_float_0 _S503;
    (&_S503)->primal_0 = r_9;
    (&_S503)->differential_0 = 0.0f;
    DiffPair_float_0 _S504;
    (&_S504)->primal_0 = _S462;
    (&_S504)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S503, &_S504, k_1);
    float _S505 = _S504.differential_0 + _S472;
    float _S506 = _S503.differential_0 + _S473;
    float2  _S507 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S508;
    (&_S508)->primal_0 = _S461;
    (&_S508)->differential_0 = _S507;
    s_bwd_length_impl_1(&_S508, _S506);
    float2  _S509 = _S508.differential_0 + _S496;
    float3  _S510 = make_float3 (_S509.x, _S509.y, _S505);
    J_4[int(0)] = _S510;
    float2  seed_7 = _S468;
    *&((&seed_7)->y) = 1.0f;
    float2  _S511 = seed_7;
    if(_S471)
    {
        float _S512 = 1.0f - _S470 * _S470 / 3.0f;
        float _S513 = _S462 * _S462;
        k_1 = _S512 / _S462;
        _S472 = 0.0f;
        _S473 = _S513;
        _S474 = _S512;
    }
    else
    {
        float _S514 = r_9 * r_9;
        k_1 = _S470 / r_9;
        _S472 = _S514;
        _S473 = 0.0f;
        _S474 = 0.0f;
    }
    float2  _S515 = make_float2 (k_1);
    float2  _S516 = _S461 * make_float2 (k_1);
    float u_23 = _S516.x;
    float v_23 = _S516.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float _S517 = (*dist_coeffs_11)[int(2)] + r2_23 * (*dist_coeffs_11)[int(3)];
    float _S518 = (*dist_coeffs_11)[int(1)] + r2_23 * _S517;
    float _S519 = (*dist_coeffs_11)[int(0)] + r2_23 * _S518;
    float _S520 = fx_4 * _S511.x;
    float2  _S521 = make_float2 (_S520, fy_4 * _S511.y) + make_float2 ((*dist_coeffs_11)[int(8)] * _S520, (*dist_coeffs_11)[int(9)] * _S520);
    float2  _S522 = _S516 * _S521;
    float _S523 = (*dist_coeffs_11)[int(4)] * _S521.y;
    float _S524 = (*dist_coeffs_11)[int(5)] * _S521.x;
    float _S525 = _S522.x + _S522.y;
    float _S526 = r2_23 * _S525;
    float _S527 = r2_23 * _S526;
    float _S528 = (*dist_coeffs_11)[int(7)] * _S521.y + _S523 + (*dist_coeffs_11)[int(6)] * _S521.x + _S524 + _S519 * _S525 + _S518 * _S526 + _S517 * _S527 + (*dist_coeffs_11)[int(3)] * (r2_23 * _S527);
    float _S529 = v_23 * _S528;
    float _S530 = u_23 * _S528;
    float2  _S531 = make_float2 (1.0f + r2_23 * _S519) * _S521 + make_float2 (_S465 * (v_23 * _S521.y) + 2.0f * u_23 * _S524 + 2.0f * (u_23 * _S524) + _S464 * (v_23 * _S521.x) + _S530 + _S530, 2.0f * v_23 * _S523 + 2.0f * (v_23 * _S523) + _S465 * u_23 * _S521.y + _S464 * u_23 * _S521.x + _S529 + _S529);
    float2  _S532 = _S461 * _S531;
    float2  _S533 = _S515 * _S531;
    float _S534 = _S532.x + _S532.y;
    if(_S471)
    {
        float _S535 = _S534 / _S473;
        float _S536 = _S474 * - _S535;
        float _S537 = _S470 * (0.3333333432674408f * - (_S462 * _S535));
        k_1 = _S537 + _S537;
        _S472 = _S536;
        _S473 = 0.0f;
    }
    else
    {
        float _S538 = _S534 / _S472;
        float _S539 = _S470 * - _S538;
        k_1 = r_9 * _S538;
        _S472 = 0.0f;
        _S473 = _S539;
    }
    DiffPair_float_0 _S540;
    (&_S540)->primal_0 = r_9;
    (&_S540)->differential_0 = 0.0f;
    DiffPair_float_0 _S541;
    (&_S541)->primal_0 = _S462;
    (&_S541)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S540, &_S541, k_1);
    float _S542 = _S541.differential_0 + _S472;
    float _S543 = _S540.differential_0 + _S473;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S544;
    (&_S544)->primal_0 = _S461;
    (&_S544)->differential_0 = _S507;
    s_bwd_length_impl_1(&_S544, _S543);
    float2  _S545 = _S544.differential_0 + _S533;
    float3  _S546 = make_float3 (_S545.x, _S545.y, _S542);
    J_4[int(1)] = _S546;
    *cov2d_4 = mul_6(mul_5(J_4, cov3d_2), transpose_1(J_4));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_2, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_12, Matrix<float, 2, 2>  * cov2d_5, float2  * mean2d_5)
{
    float2  * _S547;
    bool _S548;
    float2  * _S549;
    bool _S550;
    float2  * _S551;
    bool _S552;
    bool _S553;
    float2  * _S554;
    bool _S555;
    float2  * _S556;
    bool _S557;
    float2  * _S558;
    bool _S559;
    bool _S560;
    float2  * _S561;
    bool _S562;
    bool _S563;
    int2  _S564 = make_int2 (int(0));
    float2  _S565 = make_float2 ((float)_S564.x, (float)_S564.y);
    *mean2d_5 = _S565;
    *cov2d_5 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_2;
    for(;;)
    {
        float k_2;
        _S547 = &proj_points_2[int(0)];
        for(;;)
        {
            float2  _S566 = float2 {sigmas_2->p_0[int(0)].x, sigmas_2->p_0[int(0)].y};
            float r_10 = length_0(_S566);
            float _S567 = sigmas_2->p_0[int(0)].z;
            float theta_5 = (F32_atan2((r_10), (_S567)));
            if(theta_5 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_5 * theta_5 / 3.0f) / _S567;
            }
            else
            {
                k_2 = theta_5 / r_10;
            }
            float2  _S568 = _S566 * make_float2 (k_2);
            proj_points_2[int(0)] = _S568;
            bool _S569 = is_valid_distortion(_S568, dist_coeffs_12);
            bool _S570 = !_S569;
            _S548 = _S570;
            if(_S570)
            {
                break;
            }
            float u_24 = proj_points_2[int(0)].x;
            float v_24 = proj_points_2[int(0)].y;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float2  _S571 = proj_points_2[int(0)] * make_float2 (1.0f + r2_24 * ((*dist_coeffs_12)[int(0)] + r2_24 * ((*dist_coeffs_12)[int(1)] + r2_24 * ((*dist_coeffs_12)[int(2)] + r2_24 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_24 * v_24 + (*dist_coeffs_12)[int(5)] * (r2_24 + 2.0f * u_24 * u_24) + (*dist_coeffs_12)[int(6)] * r2_24, 2.0f * (*dist_coeffs_12)[int(5)] * u_24 * v_24 + (*dist_coeffs_12)[int(4)] * (r2_24 + 2.0f * v_24 * v_24) + (*dist_coeffs_12)[int(7)] * r2_24);
            float2  _S572 = _S571 + make_float2 ((*dist_coeffs_12)[int(8)] * _S571.x + (*dist_coeffs_12)[int(9)] * _S571.y, 0.0f);
            proj_points_2[int(0)] = make_float2 (intrins_4.x * _S572.x + intrins_4.z, intrins_4.y * _S572.y + intrins_4.w);
            break;
        }
        bool all_valid_4 = true & (!_S548);
        _S549 = &proj_points_2[int(1)];
        for(;;)
        {
            float2  _S573 = float2 {sigmas_2->p_0[int(1)].x, sigmas_2->p_0[int(1)].y};
            float r_11 = length_0(_S573);
            float _S574 = sigmas_2->p_0[int(1)].z;
            float theta_6 = (F32_atan2((r_11), (_S574)));
            if(theta_6 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_6 * theta_6 / 3.0f) / _S574;
            }
            else
            {
                k_2 = theta_6 / r_11;
            }
            float2  _S575 = _S573 * make_float2 (k_2);
            proj_points_2[int(1)] = _S575;
            bool _S576 = is_valid_distortion(_S575, dist_coeffs_12);
            bool _S577 = !_S576;
            _S550 = _S577;
            if(_S577)
            {
                break;
            }
            float u_25 = proj_points_2[int(1)].x;
            float v_25 = proj_points_2[int(1)].y;
            float r2_25 = u_25 * u_25 + v_25 * v_25;
            float2  _S578 = proj_points_2[int(1)] * make_float2 (1.0f + r2_25 * ((*dist_coeffs_12)[int(0)] + r2_25 * ((*dist_coeffs_12)[int(1)] + r2_25 * ((*dist_coeffs_12)[int(2)] + r2_25 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_25 * v_25 + (*dist_coeffs_12)[int(5)] * (r2_25 + 2.0f * u_25 * u_25) + (*dist_coeffs_12)[int(6)] * r2_25, 2.0f * (*dist_coeffs_12)[int(5)] * u_25 * v_25 + (*dist_coeffs_12)[int(4)] * (r2_25 + 2.0f * v_25 * v_25) + (*dist_coeffs_12)[int(7)] * r2_25);
            float2  _S579 = _S578 + make_float2 ((*dist_coeffs_12)[int(8)] * _S578.x + (*dist_coeffs_12)[int(9)] * _S578.y, 0.0f);
            proj_points_2[int(1)] = make_float2 (intrins_4.x * _S579.x + intrins_4.z, intrins_4.y * _S579.y + intrins_4.w);
            break;
        }
        bool all_valid_5 = all_valid_4 & (!_S550);
        for(;;)
        {
            _S551 = &proj_points_2[int(2)];
            for(;;)
            {
                float2  _S580 = float2 {sigmas_2->p_0[int(2)].x, sigmas_2->p_0[int(2)].y};
                float r_12 = length_0(_S580);
                float _S581 = sigmas_2->p_0[int(2)].z;
                float theta_7 = (F32_atan2((r_12), (_S581)));
                if(theta_7 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_7 * theta_7 / 3.0f) / _S581;
                }
                else
                {
                    k_2 = theta_7 / r_12;
                }
                float2  _S582 = _S580 * make_float2 (k_2);
                proj_points_2[int(2)] = _S582;
                bool _S583 = is_valid_distortion(_S582, dist_coeffs_12);
                bool _S584 = !_S583;
                _S552 = _S584;
                if(_S584)
                {
                    break;
                }
                float u_26 = proj_points_2[int(2)].x;
                float v_26 = proj_points_2[int(2)].y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S585 = proj_points_2[int(2)] * make_float2 (1.0f + r2_26 * ((*dist_coeffs_12)[int(0)] + r2_26 * ((*dist_coeffs_12)[int(1)] + r2_26 * ((*dist_coeffs_12)[int(2)] + r2_26 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_26 * v_26 + (*dist_coeffs_12)[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + (*dist_coeffs_12)[int(6)] * r2_26, 2.0f * (*dist_coeffs_12)[int(5)] * u_26 * v_26 + (*dist_coeffs_12)[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + (*dist_coeffs_12)[int(7)] * r2_26);
                float2  _S586 = _S585 + make_float2 ((*dist_coeffs_12)[int(8)] * _S585.x + (*dist_coeffs_12)[int(9)] * _S585.y, 0.0f);
                proj_points_2[int(2)] = make_float2 (intrins_4.x * _S586.x + intrins_4.z, intrins_4.y * _S586.y + intrins_4.w);
                break;
            }
            _S553 = all_valid_5 & (!_S552);
            break;
        }
        _S554 = &proj_points_2[int(3)];
        for(;;)
        {
            float2  _S587 = float2 {sigmas_2->p_0[int(3)].x, sigmas_2->p_0[int(3)].y};
            float r_13 = length_0(_S587);
            float _S588 = sigmas_2->p_0[int(3)].z;
            float theta_8 = (F32_atan2((r_13), (_S588)));
            if(theta_8 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_8 * theta_8 / 3.0f) / _S588;
            }
            else
            {
                k_2 = theta_8 / r_13;
            }
            float2  _S589 = _S587 * make_float2 (k_2);
            proj_points_2[int(3)] = _S589;
            bool _S590 = is_valid_distortion(_S589, dist_coeffs_12);
            bool _S591 = !_S590;
            _S555 = _S591;
            if(_S591)
            {
                break;
            }
            float u_27 = proj_points_2[int(3)].x;
            float v_27 = proj_points_2[int(3)].y;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float2  _S592 = proj_points_2[int(3)] * make_float2 (1.0f + r2_27 * ((*dist_coeffs_12)[int(0)] + r2_27 * ((*dist_coeffs_12)[int(1)] + r2_27 * ((*dist_coeffs_12)[int(2)] + r2_27 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_27 * v_27 + (*dist_coeffs_12)[int(5)] * (r2_27 + 2.0f * u_27 * u_27) + (*dist_coeffs_12)[int(6)] * r2_27, 2.0f * (*dist_coeffs_12)[int(5)] * u_27 * v_27 + (*dist_coeffs_12)[int(4)] * (r2_27 + 2.0f * v_27 * v_27) + (*dist_coeffs_12)[int(7)] * r2_27);
            float2  _S593 = _S592 + make_float2 ((*dist_coeffs_12)[int(8)] * _S592.x + (*dist_coeffs_12)[int(9)] * _S592.y, 0.0f);
            proj_points_2[int(3)] = make_float2 (intrins_4.x * _S593.x + intrins_4.z, intrins_4.y * _S593.y + intrins_4.w);
            break;
        }
        bool all_valid_6 = _S553 & (!_S555);
        _S556 = &proj_points_2[int(4)];
        for(;;)
        {
            float2  _S594 = float2 {sigmas_2->p_0[int(4)].x, sigmas_2->p_0[int(4)].y};
            float r_14 = length_0(_S594);
            float _S595 = sigmas_2->p_0[int(4)].z;
            float theta_9 = (F32_atan2((r_14), (_S595)));
            if(theta_9 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_9 * theta_9 / 3.0f) / _S595;
            }
            else
            {
                k_2 = theta_9 / r_14;
            }
            float2  _S596 = _S594 * make_float2 (k_2);
            proj_points_2[int(4)] = _S596;
            bool _S597 = is_valid_distortion(_S596, dist_coeffs_12);
            bool _S598 = !_S597;
            _S557 = _S598;
            if(_S598)
            {
                break;
            }
            float u_28 = proj_points_2[int(4)].x;
            float v_28 = proj_points_2[int(4)].y;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float2  _S599 = proj_points_2[int(4)] * make_float2 (1.0f + r2_28 * ((*dist_coeffs_12)[int(0)] + r2_28 * ((*dist_coeffs_12)[int(1)] + r2_28 * ((*dist_coeffs_12)[int(2)] + r2_28 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_28 * v_28 + (*dist_coeffs_12)[int(5)] * (r2_28 + 2.0f * u_28 * u_28) + (*dist_coeffs_12)[int(6)] * r2_28, 2.0f * (*dist_coeffs_12)[int(5)] * u_28 * v_28 + (*dist_coeffs_12)[int(4)] * (r2_28 + 2.0f * v_28 * v_28) + (*dist_coeffs_12)[int(7)] * r2_28);
            float2  _S600 = _S599 + make_float2 ((*dist_coeffs_12)[int(8)] * _S599.x + (*dist_coeffs_12)[int(9)] * _S599.y, 0.0f);
            proj_points_2[int(4)] = make_float2 (intrins_4.x * _S600.x + intrins_4.z, intrins_4.y * _S600.y + intrins_4.w);
            break;
        }
        bool all_valid_7 = all_valid_6 & (!_S557);
        for(;;)
        {
            _S558 = &proj_points_2[int(5)];
            for(;;)
            {
                float2  _S601 = float2 {sigmas_2->p_0[int(5)].x, sigmas_2->p_0[int(5)].y};
                float r_15 = length_0(_S601);
                float _S602 = sigmas_2->p_0[int(5)].z;
                float theta_10 = (F32_atan2((r_15), (_S602)));
                if(theta_10 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_10 * theta_10 / 3.0f) / _S602;
                }
                else
                {
                    k_2 = theta_10 / r_15;
                }
                float2  _S603 = _S601 * make_float2 (k_2);
                proj_points_2[int(5)] = _S603;
                bool _S604 = is_valid_distortion(_S603, dist_coeffs_12);
                bool _S605 = !_S604;
                _S559 = _S605;
                if(_S605)
                {
                    break;
                }
                float u_29 = proj_points_2[int(5)].x;
                float v_29 = proj_points_2[int(5)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S606 = proj_points_2[int(5)] * make_float2 (1.0f + r2_29 * ((*dist_coeffs_12)[int(0)] + r2_29 * ((*dist_coeffs_12)[int(1)] + r2_29 * ((*dist_coeffs_12)[int(2)] + r2_29 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_29 * v_29 + (*dist_coeffs_12)[int(5)] * (r2_29 + 2.0f * u_29 * u_29) + (*dist_coeffs_12)[int(6)] * r2_29, 2.0f * (*dist_coeffs_12)[int(5)] * u_29 * v_29 + (*dist_coeffs_12)[int(4)] * (r2_29 + 2.0f * v_29 * v_29) + (*dist_coeffs_12)[int(7)] * r2_29);
                float2  _S607 = _S606 + make_float2 ((*dist_coeffs_12)[int(8)] * _S606.x + (*dist_coeffs_12)[int(9)] * _S606.y, 0.0f);
                proj_points_2[int(5)] = make_float2 (intrins_4.x * _S607.x + intrins_4.z, intrins_4.y * _S607.y + intrins_4.w);
                break;
            }
            _S560 = all_valid_7 & (!_S559);
            break;
        }
        _S561 = &proj_points_2[int(6)];
        for(;;)
        {
            float2  _S608 = float2 {sigmas_2->p_0[int(6)].x, sigmas_2->p_0[int(6)].y};
            float r_16 = length_0(_S608);
            float _S609 = sigmas_2->p_0[int(6)].z;
            float theta_11 = (F32_atan2((r_16), (_S609)));
            if(theta_11 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_11 * theta_11 / 3.0f) / _S609;
            }
            else
            {
                k_2 = theta_11 / r_16;
            }
            float2  _S610 = _S608 * make_float2 (k_2);
            proj_points_2[int(6)] = _S610;
            bool _S611 = is_valid_distortion(_S610, dist_coeffs_12);
            bool _S612 = !_S611;
            _S562 = _S612;
            if(_S612)
            {
                break;
            }
            float u_30 = proj_points_2[int(6)].x;
            float v_30 = proj_points_2[int(6)].y;
            float r2_30 = u_30 * u_30 + v_30 * v_30;
            float2  _S613 = proj_points_2[int(6)] * make_float2 (1.0f + r2_30 * ((*dist_coeffs_12)[int(0)] + r2_30 * ((*dist_coeffs_12)[int(1)] + r2_30 * ((*dist_coeffs_12)[int(2)] + r2_30 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_30 * v_30 + (*dist_coeffs_12)[int(5)] * (r2_30 + 2.0f * u_30 * u_30) + (*dist_coeffs_12)[int(6)] * r2_30, 2.0f * (*dist_coeffs_12)[int(5)] * u_30 * v_30 + (*dist_coeffs_12)[int(4)] * (r2_30 + 2.0f * v_30 * v_30) + (*dist_coeffs_12)[int(7)] * r2_30);
            float2  _S614 = _S613 + make_float2 ((*dist_coeffs_12)[int(8)] * _S613.x + (*dist_coeffs_12)[int(9)] * _S613.y, 0.0f);
            proj_points_2[int(6)] = make_float2 (intrins_4.x * _S614.x + intrins_4.z, intrins_4.y * _S614.y + intrins_4.w);
            break;
        }
        _S563 = _S560 & (!_S562);
        break;
    }
    if(!_S563)
    {
        return false;
    }
    float2  _S615 = *mean2d_5 + make_float2 (sigmas_2->w_mean_0[int(0)]) * *_S547;
    *mean2d_5 = _S615;
    float2  _S616 = _S615 + make_float2 (sigmas_2->w_mean_0[int(1)]) * *_S549;
    *mean2d_5 = _S616;
    float2  _S617 = _S616 + make_float2 (sigmas_2->w_mean_0[int(2)]) * *_S551;
    *mean2d_5 = _S617;
    float2  _S618 = _S617 + make_float2 (sigmas_2->w_mean_0[int(3)]) * *_S554;
    *mean2d_5 = _S618;
    float2  _S619 = _S618 + make_float2 (sigmas_2->w_mean_0[int(4)]) * *_S556;
    *mean2d_5 = _S619;
    float2  _S620 = _S619 + make_float2 (sigmas_2->w_mean_0[int(5)]) * *_S558;
    *mean2d_5 = _S620;
    float2  _S621 = _S620 + make_float2 (sigmas_2->w_mean_0[int(6)]) * *_S561;
    *mean2d_5 = _S621;
    float2  d_14 = *_S547 - _S621;
    float _S622 = d_14.x;
    float _S623 = d_14.y;
    float _S624 = _S622 * _S623;
    Matrix<float, 2, 2>  _S625 = *cov2d_5 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S622 * _S622, _S624, _S624, _S623 * _S623);
    *cov2d_5 = _S625;
    float2  d_15 = *_S549 - *mean2d_5;
    float _S626 = d_15.x;
    float _S627 = d_15.y;
    float _S628 = _S626 * _S627;
    Matrix<float, 2, 2>  _S629 = _S625 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S626 * _S626, _S628, _S628, _S627 * _S627);
    *cov2d_5 = _S629;
    float2  d_16 = *_S551 - *mean2d_5;
    float _S630 = d_16.x;
    float _S631 = d_16.y;
    float _S632 = _S630 * _S631;
    Matrix<float, 2, 2>  _S633 = _S629 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S630 * _S630, _S632, _S632, _S631 * _S631);
    *cov2d_5 = _S633;
    float2  d_17 = *_S554 - *mean2d_5;
    float _S634 = d_17.x;
    float _S635 = d_17.y;
    float _S636 = _S634 * _S635;
    Matrix<float, 2, 2>  _S637 = _S633 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S634 * _S634, _S636, _S636, _S635 * _S635);
    *cov2d_5 = _S637;
    float2  d_18 = *_S556 - *mean2d_5;
    float _S638 = d_18.x;
    float _S639 = d_18.y;
    float _S640 = _S638 * _S639;
    Matrix<float, 2, 2>  _S641 = _S637 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S638 * _S638, _S640, _S640, _S639 * _S639);
    *cov2d_5 = _S641;
    float2  d_19 = *_S558 - *mean2d_5;
    float _S642 = d_19.x;
    float _S643 = d_19.y;
    float _S644 = _S642 * _S643;
    Matrix<float, 2, 2>  _S645 = _S641 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S642 * _S642, _S644, _S644, _S643 * _S643);
    *cov2d_5 = _S645;
    float2  d_20 = *_S561 - *mean2d_5;
    float _S646 = d_20.x;
    float _S647 = d_20.y;
    float _S648 = _S646 * _S647;
    *cov2d_5 = _S645 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S646 * _S646, _S648, _S648, _S647 * _S647);
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_3, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_13, Matrix<float, 2, 2>  * cov2d_6, float2  * mean2d_6)
{
    int2  _S649 = make_int2 (int(0));
    float2  _S650 = make_float2 ((float)_S649.x, (float)_S649.y);
    *mean2d_6 = _S650;
    *cov2d_6 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_3;
    float2  _S651 = float2 {sigmas_3->p_0[int(0)].x, sigmas_3->p_0[int(0)].y};
    float r_17 = length_0(_S651);
    float _S652 = sigmas_3->p_0[int(0)].z;
    float theta_12 = (F32_atan2((r_17), (_S652)));
    float k_3;
    if(theta_12 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_12 * theta_12 / 3.0f) / _S652;
    }
    else
    {
        k_3 = theta_12 / r_17;
    }
    float2  _S653 = _S651 * make_float2 (k_3);
    float u_31 = _S653.x;
    float v_31 = _S653.y;
    float r2_31 = u_31 * u_31 + v_31 * v_31;
    float _S654 = 2.0f * (*dist_coeffs_13)[int(4)];
    float _S655 = 2.0f * (*dist_coeffs_13)[int(5)];
    float2  _S656 = _S653 * make_float2 (1.0f + r2_31 * ((*dist_coeffs_13)[int(0)] + r2_31 * ((*dist_coeffs_13)[int(1)] + r2_31 * ((*dist_coeffs_13)[int(2)] + r2_31 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_31 * v_31 + (*dist_coeffs_13)[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + (*dist_coeffs_13)[int(6)] * r2_31, _S655 * u_31 * v_31 + (*dist_coeffs_13)[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + (*dist_coeffs_13)[int(7)] * r2_31);
    float2  _S657 = _S656 + make_float2 ((*dist_coeffs_13)[int(8)] * _S656.x + (*dist_coeffs_13)[int(9)] * _S656.y, 0.0f);
    float fx_5 = intrins_5.x;
    float fy_5 = intrins_5.y;
    float cx_3 = intrins_5.z;
    float cy_3 = intrins_5.w;
    proj_points_3[int(0)] = make_float2 (fx_5 * _S657.x + cx_3, fy_5 * _S657.y + cy_3);
    float2  _S658 = float2 {sigmas_3->p_0[int(1)].x, sigmas_3->p_0[int(1)].y};
    float r_18 = length_0(_S658);
    float _S659 = sigmas_3->p_0[int(1)].z;
    float theta_13 = (F32_atan2((r_18), (_S659)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_13 * theta_13 / 3.0f) / _S659;
    }
    else
    {
        k_3 = theta_13 / r_18;
    }
    float2  _S660 = _S658 * make_float2 (k_3);
    float u_32 = _S660.x;
    float v_32 = _S660.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float2  _S661 = _S660 * make_float2 (1.0f + r2_32 * ((*dist_coeffs_13)[int(0)] + r2_32 * ((*dist_coeffs_13)[int(1)] + r2_32 * ((*dist_coeffs_13)[int(2)] + r2_32 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_32 * v_32 + (*dist_coeffs_13)[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + (*dist_coeffs_13)[int(6)] * r2_32, _S655 * u_32 * v_32 + (*dist_coeffs_13)[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + (*dist_coeffs_13)[int(7)] * r2_32);
    float2  _S662 = _S661 + make_float2 ((*dist_coeffs_13)[int(8)] * _S661.x + (*dist_coeffs_13)[int(9)] * _S661.y, 0.0f);
    proj_points_3[int(1)] = make_float2 (fx_5 * _S662.x + cx_3, fy_5 * _S662.y + cy_3);
    float2  _S663 = float2 {sigmas_3->p_0[int(2)].x, sigmas_3->p_0[int(2)].y};
    float r_19 = length_0(_S663);
    float _S664 = sigmas_3->p_0[int(2)].z;
    float theta_14 = (F32_atan2((r_19), (_S664)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_14 * theta_14 / 3.0f) / _S664;
    }
    else
    {
        k_3 = theta_14 / r_19;
    }
    float2  _S665 = _S663 * make_float2 (k_3);
    float u_33 = _S665.x;
    float v_33 = _S665.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S666 = _S665 * make_float2 (1.0f + r2_33 * ((*dist_coeffs_13)[int(0)] + r2_33 * ((*dist_coeffs_13)[int(1)] + r2_33 * ((*dist_coeffs_13)[int(2)] + r2_33 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_33 * v_33 + (*dist_coeffs_13)[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + (*dist_coeffs_13)[int(6)] * r2_33, _S655 * u_33 * v_33 + (*dist_coeffs_13)[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + (*dist_coeffs_13)[int(7)] * r2_33);
    float2  _S667 = _S666 + make_float2 ((*dist_coeffs_13)[int(8)] * _S666.x + (*dist_coeffs_13)[int(9)] * _S666.y, 0.0f);
    proj_points_3[int(2)] = make_float2 (fx_5 * _S667.x + cx_3, fy_5 * _S667.y + cy_3);
    float2  _S668 = float2 {sigmas_3->p_0[int(3)].x, sigmas_3->p_0[int(3)].y};
    float r_20 = length_0(_S668);
    float _S669 = sigmas_3->p_0[int(3)].z;
    float theta_15 = (F32_atan2((r_20), (_S669)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_15 * theta_15 / 3.0f) / _S669;
    }
    else
    {
        k_3 = theta_15 / r_20;
    }
    float2  _S670 = _S668 * make_float2 (k_3);
    float u_34 = _S670.x;
    float v_34 = _S670.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S671 = _S670 * make_float2 (1.0f + r2_34 * ((*dist_coeffs_13)[int(0)] + r2_34 * ((*dist_coeffs_13)[int(1)] + r2_34 * ((*dist_coeffs_13)[int(2)] + r2_34 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_34 * v_34 + (*dist_coeffs_13)[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + (*dist_coeffs_13)[int(6)] * r2_34, _S655 * u_34 * v_34 + (*dist_coeffs_13)[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + (*dist_coeffs_13)[int(7)] * r2_34);
    float2  _S672 = _S671 + make_float2 ((*dist_coeffs_13)[int(8)] * _S671.x + (*dist_coeffs_13)[int(9)] * _S671.y, 0.0f);
    proj_points_3[int(3)] = make_float2 (fx_5 * _S672.x + cx_3, fy_5 * _S672.y + cy_3);
    float2  _S673 = float2 {sigmas_3->p_0[int(4)].x, sigmas_3->p_0[int(4)].y};
    float r_21 = length_0(_S673);
    float _S674 = sigmas_3->p_0[int(4)].z;
    float theta_16 = (F32_atan2((r_21), (_S674)));
    if(theta_16 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_16 * theta_16 / 3.0f) / _S674;
    }
    else
    {
        k_3 = theta_16 / r_21;
    }
    float2  _S675 = _S673 * make_float2 (k_3);
    float u_35 = _S675.x;
    float v_35 = _S675.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S676 = _S675 * make_float2 (1.0f + r2_35 * ((*dist_coeffs_13)[int(0)] + r2_35 * ((*dist_coeffs_13)[int(1)] + r2_35 * ((*dist_coeffs_13)[int(2)] + r2_35 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_35 * v_35 + (*dist_coeffs_13)[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + (*dist_coeffs_13)[int(6)] * r2_35, _S655 * u_35 * v_35 + (*dist_coeffs_13)[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + (*dist_coeffs_13)[int(7)] * r2_35);
    float2  _S677 = _S676 + make_float2 ((*dist_coeffs_13)[int(8)] * _S676.x + (*dist_coeffs_13)[int(9)] * _S676.y, 0.0f);
    proj_points_3[int(4)] = make_float2 (fx_5 * _S677.x + cx_3, fy_5 * _S677.y + cy_3);
    float2  _S678 = float2 {sigmas_3->p_0[int(5)].x, sigmas_3->p_0[int(5)].y};
    float r_22 = length_0(_S678);
    float _S679 = sigmas_3->p_0[int(5)].z;
    float theta_17 = (F32_atan2((r_22), (_S679)));
    if(theta_17 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_17 * theta_17 / 3.0f) / _S679;
    }
    else
    {
        k_3 = theta_17 / r_22;
    }
    float2  _S680 = _S678 * make_float2 (k_3);
    float u_36 = _S680.x;
    float v_36 = _S680.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S681 = _S680 * make_float2 (1.0f + r2_36 * ((*dist_coeffs_13)[int(0)] + r2_36 * ((*dist_coeffs_13)[int(1)] + r2_36 * ((*dist_coeffs_13)[int(2)] + r2_36 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_36 * v_36 + (*dist_coeffs_13)[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + (*dist_coeffs_13)[int(6)] * r2_36, _S655 * u_36 * v_36 + (*dist_coeffs_13)[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + (*dist_coeffs_13)[int(7)] * r2_36);
    float2  _S682 = _S681 + make_float2 ((*dist_coeffs_13)[int(8)] * _S681.x + (*dist_coeffs_13)[int(9)] * _S681.y, 0.0f);
    proj_points_3[int(5)] = make_float2 (fx_5 * _S682.x + cx_3, fy_5 * _S682.y + cy_3);
    float2  _S683 = float2 {sigmas_3->p_0[int(6)].x, sigmas_3->p_0[int(6)].y};
    float r_23 = length_0(_S683);
    float _S684 = sigmas_3->p_0[int(6)].z;
    float theta_18 = (F32_atan2((r_23), (_S684)));
    if(theta_18 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_18 * theta_18 / 3.0f) / _S684;
    }
    else
    {
        k_3 = theta_18 / r_23;
    }
    float2  _S685 = _S683 * make_float2 (k_3);
    float u_37 = _S685.x;
    float v_37 = _S685.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S686 = _S685 * make_float2 (1.0f + r2_37 * ((*dist_coeffs_13)[int(0)] + r2_37 * ((*dist_coeffs_13)[int(1)] + r2_37 * ((*dist_coeffs_13)[int(2)] + r2_37 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S654 * u_37 * v_37 + (*dist_coeffs_13)[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + (*dist_coeffs_13)[int(6)] * r2_37, _S655 * u_37 * v_37 + (*dist_coeffs_13)[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + (*dist_coeffs_13)[int(7)] * r2_37);
    float2  _S687 = _S686 + make_float2 ((*dist_coeffs_13)[int(8)] * _S686.x + (*dist_coeffs_13)[int(9)] * _S686.y, 0.0f);
    float2  _S688 = make_float2 (fx_5 * _S687.x + cx_3, fy_5 * _S687.y + cy_3);
    proj_points_3[int(6)] = _S688;
    float2  _S689 = *mean2d_6 + make_float2 (sigmas_3->w_mean_0[int(0)]) * proj_points_3[int(0)];
    *mean2d_6 = _S689;
    float2  _S690 = _S689 + make_float2 (sigmas_3->w_mean_0[int(1)]) * proj_points_3[int(1)];
    *mean2d_6 = _S690;
    float2  _S691 = _S690 + make_float2 (sigmas_3->w_mean_0[int(2)]) * proj_points_3[int(2)];
    *mean2d_6 = _S691;
    float2  _S692 = _S691 + make_float2 (sigmas_3->w_mean_0[int(3)]) * proj_points_3[int(3)];
    *mean2d_6 = _S692;
    float2  _S693 = _S692 + make_float2 (sigmas_3->w_mean_0[int(4)]) * proj_points_3[int(4)];
    *mean2d_6 = _S693;
    float2  _S694 = _S693 + make_float2 (sigmas_3->w_mean_0[int(5)]) * proj_points_3[int(5)];
    *mean2d_6 = _S694;
    float2  _S695 = _S694 + make_float2 (sigmas_3->w_mean_0[int(6)]) * _S688;
    *mean2d_6 = _S695;
    float2  d_21 = proj_points_3[int(0)] - _S695;
    float _S696 = d_21.x;
    float _S697 = d_21.y;
    float _S698 = _S696 * _S697;
    Matrix<float, 2, 2>  _S699 = *cov2d_6 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S696 * _S696, _S698, _S698, _S697 * _S697);
    *cov2d_6 = _S699;
    float2  d_22 = proj_points_3[int(1)] - *mean2d_6;
    float _S700 = d_22.x;
    float _S701 = d_22.y;
    float _S702 = _S700 * _S701;
    Matrix<float, 2, 2>  _S703 = _S699 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S700 * _S700, _S702, _S702, _S701 * _S701);
    *cov2d_6 = _S703;
    float2  d_23 = proj_points_3[int(2)] - *mean2d_6;
    float _S704 = d_23.x;
    float _S705 = d_23.y;
    float _S706 = _S704 * _S705;
    Matrix<float, 2, 2>  _S707 = _S703 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S704 * _S704, _S706, _S706, _S705 * _S705);
    *cov2d_6 = _S707;
    float2  d_24 = proj_points_3[int(3)] - *mean2d_6;
    float _S708 = d_24.x;
    float _S709 = d_24.y;
    float _S710 = _S708 * _S709;
    Matrix<float, 2, 2>  _S711 = _S707 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S708 * _S708, _S710, _S710, _S709 * _S709);
    *cov2d_6 = _S711;
    float2  d_25 = proj_points_3[int(4)] - *mean2d_6;
    float _S712 = d_25.x;
    float _S713 = d_25.y;
    float _S714 = _S712 * _S713;
    Matrix<float, 2, 2>  _S715 = _S711 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S712 * _S712, _S714, _S714, _S713 * _S713);
    *cov2d_6 = _S715;
    float2  d_26 = proj_points_3[int(5)] - *mean2d_6;
    float _S716 = d_26.x;
    float _S717 = d_26.y;
    float _S718 = _S716 * _S717;
    Matrix<float, 2, 2>  _S719 = _S715 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S716 * _S716, _S718, _S718, _S717 * _S717);
    *cov2d_6 = _S719;
    float2  d_27 = _S688 - *mean2d_6;
    float _S720 = d_27.x;
    float _S721 = d_27.y;
    float _S722 = _S720 * _S721;
    *cov2d_6 = _S719 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S720 * _S720, _S722, _S722, _S721 * _S721);
    return true;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_6, float fy_6, float cx_4, float cy_4, Matrix<float, 2, 2>  * cov2d_7, float2  * mean2d_7)
{
    Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
    *cov2d_7 = mul_6(mul_5(J_5, cov3d_3), transpose_1(J_5));
    *mean2d_7 = make_float2 (fx_6 * mean3d_3.x + cx_4, fy_6 * mean3d_3.y + cy_4);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    float _S723 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S723;
    float _S724 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S724;
    float det_blur_0 = _S723 * _S724 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_10, float dOut_12)
{
    float _S725 = (F32_exp(((*dpx_10).primal_0))) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S725;
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
    float3  _S726 = exp_0((*dpx_11).primal_0) * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S726;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_12, float dOut_14)
{
    float _S727 = 1.0f / (*dpx_12).primal_0 * dOut_14;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S727;
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

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_7, float fy_7, float cx_5, float cy_5, FixedArray<float, 10>  * dist_coeffs_14, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_8, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S728 = mean_c_0.z;
        bool _S729;
        if(_S728 < near_plane_0)
        {
            _S729 = true;
        }
        else
        {
            _S729 = _S728 > far_plane_0;
        }
        if(_S729)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        float3  _S730 = exp_0(scale_2);
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
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S730.x, 0.0f, 0.0f, 0.0f, _S730.y, 0.0f, 0.0f, 0.0f, _S730.z));
        Matrix<float, 3, 3>  _S731 = transpose_0(R_4);
        float _S732 = float(image_width_0);
        float _S733 = float(image_height_0);
        float _S734 = 0.30000001192092896f * (0.5f * _S732 / fx_7);
        float _S735 = 0.30000001192092896f * (0.5f * _S733 / fy_7);
        float rz_1 = 1.0f / mean_c_0.z;
        float rz2_1 = rz_1 * rz_1;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_7 * rz_1, 0.0f, - fx_7 * (mean_c_0.z * (F32_min(((_S732 - cx_5) / fx_7 + _S734), ((F32_max((- (cx_5 / fx_7 + _S734)), (mean_c_0.x * rz_1))))))) * rz2_1, 0.0f, fy_7 * rz_1, - fy_7 * (mean_c_0.z * (F32_min(((_S733 - cy_5) / fy_7 + _S735), ((F32_max((- (cy_5 / fy_7 + _S735)), (mean_c_0.y * rz_1))))))) * rz2_1);
        covar2d_0 = mul_6(mul_5(J_6, mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S731)), transpose_1(J_6));
        *mean2d_8 = make_float2 (fx_7 * mean_c_0.x * rz_1 + cx_5, fy_7 * mean_c_0.y * rz_1 + cy_5);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S736 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S736;
        float _S737 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S737;
        float det_blur_1 = _S736 * _S737 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S738 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S729 = true;
        }
        else
        {
            _S729 = xmin_0 >= _S732;
        }
        if(_S729)
        {
            _S729 = true;
        }
        else
        {
            _S729 = ymax_0 <= 0.0f;
        }
        if(_S729)
        {
            _S729 = true;
        }
        else
        {
            _S729 = ymin_0 >= _S733;
        }
        if(_S729)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S738.rows[int(0)].x, _S738.rows[int(0)].y, _S738.rows[int(1)].y);
        float3  _S739 = mean_0 - - mul_0(_S731, t_3);
        float3  _S740 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S740;
        float _S741 = _S739.x;
        float _S742 = _S739.y;
        float _S743 = _S739.z;
        float norm_0 = (F32_sqrt((_S741 * _S741 + _S742 * _S742 + _S743 * _S743)));
        float x_16 = _S741 / norm_0;
        float y_3 = _S742 / norm_0;
        float z_0 = _S743 / norm_0;
        float3  _S744 = _S740 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_16) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S744;
        float z2_4 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_16 * x_16 - y_3 * y_3;
        float fS1_0 = 2.0f * x_16 * y_3;
        float3  _S745 = _S744 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_16) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S745;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S745 + (make_float3 (-0.59004360437393188f * (x_16 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_16) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_16 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_8, float fy_8, float cx_6, float cy_6, FixedArray<float, 10>  * dist_coeffs_15, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_9, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S746 = length_1(mean_c_1);
        bool _S747;
        if(_S746 < near_plane_1)
        {
            _S747 = true;
        }
        else
        {
            _S747 = _S746 > far_plane_1;
        }
        if(_S747)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_6 = make_float4 (fx_8, fy_8, cx_6, cy_6);
        Matrix<float, 2, 2>  covar2d_1;
        float3  _S748 = exp_0(scale_3);
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
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S748.x, 0.0f, 0.0f, 0.0f, _S748.y, 0.0f, 0.0f, 0.0f, _S748.z));
        Matrix<float, 3, 3>  _S749 = transpose_0(R_5);
        bool _S750 = fisheye_proj_3dgs_0(mean_c_1, mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S749), intrins_6, dist_coeffs_15, &covar2d_1, mean2d_9);
        if(!(true & _S750))
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S751 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S751;
        float _S752 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S752;
        float det_blur_2 = _S751 * _S752 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S753 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S747 = true;
        }
        else
        {
            _S747 = xmin_1 >= float(image_width_1);
        }
        if(_S747)
        {
            _S747 = true;
        }
        else
        {
            _S747 = ymax_1 <= 0.0f;
        }
        if(_S747)
        {
            _S747 = true;
        }
        else
        {
            _S747 = ymin_1 >= float(image_height_1);
        }
        if(_S747)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S753.rows[int(0)].x, _S753.rows[int(0)].y, _S753.rows[int(1)].y);
        float3  _S754 = mean_1 - - mul_0(_S749, t_4);
        float3  _S755 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S755;
        float _S756 = _S754.x;
        float _S757 = _S754.y;
        float _S758 = _S754.z;
        float norm_1 = (F32_sqrt((_S756 * _S756 + _S757 * _S757 + _S758 * _S758)));
        float x_18 = _S756 / norm_1;
        float y_4 = _S757 / norm_1;
        float z_1 = _S758 / norm_1;
        float3  _S759 = _S755 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_18) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S759;
        float z2_6 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_18 * x_18 - y_4 * y_4;
        float fS1_1 = 2.0f * x_18 * y_4;
        float3  _S760 = _S759 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_18) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S760;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S760 + (make_float3 (-0.59004360437393188f * (x_18 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_18) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_18 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_9, float fy_9, float cx_7, float cy_7, FixedArray<float, 10>  * dist_coeffs_16, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_10, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S761 = mean_c_2.z;
        bool _S762;
        if(_S761 < near_plane_2)
        {
            _S762 = true;
        }
        else
        {
            _S762 = _S761 > far_plane_2;
        }
        if(_S762)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        float3  _S763 = exp_0(scale_4);
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
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S763.x, 0.0f, 0.0f, 0.0f, _S763.y, 0.0f, 0.0f, 0.0f, _S763.z));
        Matrix<float, 3, 3>  _S764 = transpose_0(R_6);
        Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
        covar2d_2 = mul_6(mul_5(J_7, mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S764)), transpose_1(J_7));
        *mean2d_10 = make_float2 (fx_9 * mean_c_2.x + cx_7, fy_9 * mean_c_2.y + cy_7);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S765 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S765;
        float _S766 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S766;
        float det_blur_3 = _S765 * _S766 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S767 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S762 = true;
        }
        else
        {
            _S762 = xmin_2 >= float(image_width_2);
        }
        if(_S762)
        {
            _S762 = true;
        }
        else
        {
            _S762 = ymax_2 <= 0.0f;
        }
        if(_S762)
        {
            _S762 = true;
        }
        else
        {
            _S762 = ymin_2 >= float(image_height_2);
        }
        if(_S762)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S767.rows[int(0)].x, _S767.rows[int(0)].y, _S767.rows[int(1)].y);
        float3  _S768 = mean_2 - - mul_0(_S764, t_5);
        float3  _S769 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S769;
        float _S770 = _S768.x;
        float _S771 = _S768.y;
        float _S772 = _S768.z;
        float norm_2 = (F32_sqrt((_S770 * _S770 + _S771 * _S771 + _S772 * _S772)));
        float x_20 = _S770 / norm_2;
        float y_5 = _S771 / norm_2;
        float z_2 = _S772 / norm_2;
        float3  _S773 = _S769 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_20) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S773;
        float z2_8 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_20 * x_20 - y_5 * y_5;
        float fS1_2 = 2.0f * x_20 * y_5;
        float3  _S774 = _S773 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_20) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S774;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S774 + (make_float3 (-0.59004360437393188f * (x_20 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_20) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_10, float fy_10, float cx_8, float cy_8, FixedArray<float, 10>  * dist_coeffs_17, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_11, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
        float _S775 = mean_c_3.z;
        bool _S776;
        if(_S775 < near_plane_3)
        {
            _S776 = true;
        }
        else
        {
            _S776 = _S775 > far_plane_3;
        }
        if(_S776)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_7 = make_float4 (fx_10, fy_10, cx_8, cy_8);
        Matrix<float, 2, 2>  covar2d_3;
        float3  _S777 = exp_0(scale_5);
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
        Matrix<float, 3, 3>  _S778 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_3;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S779 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S779 * _S777.x) * _S778.rows[0U];
        float3  _S780 = mean_3 + delta_0;
        float3  _S781 = mean_3 - delta_0;
        float3  delta_1 = make_float3 (_S779 * _S777.y) * _S778.rows[1U];
        float3  _S782 = mean_3 + delta_1;
        float3  _S783 = mean_3 - delta_1;
        float3  delta_2 = make_float3 (_S779 * _S777.z) * _S778.rows[2U];
        float3  _S784 = mean_3 + delta_2;
        float3  _S785 = mean_3 - delta_2;
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
        (&ret_0)->p_0[1U] = mul_0(R_7, _S780) + t_6;
        (&ret_0)->p_0[2U] = mul_0(R_7, _S782) + t_6;
        (&ret_0)->p_0[3U] = mul_0(R_7, _S784) + t_6;
        (&ret_0)->p_0[4U] = mul_0(R_7, _S781) + t_6;
        (&ret_0)->p_0[5U] = mul_0(R_7, _S783) + t_6;
        (&ret_0)->p_0[6U] = mul_0(R_7, _S785) + t_6;
        SigmaPoints_0 _S786 = ret_0;
        bool _S787 = persp_proj_3dgs_ut_0(&_S786, intrins_7, dist_coeffs_17, image_width_3, image_height_3, &covar2d_3, mean2d_11);
        if(!(true & _S787))
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S788 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S788;
        float _S789 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S789;
        float det_blur_4 = _S788 * _S789 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
            _S776 = true;
        }
        else
        {
            _S776 = xmin_3 >= float(image_width_3);
        }
        if(_S776)
        {
            _S776 = true;
        }
        else
        {
            _S776 = ymax_3 <= 0.0f;
        }
        if(_S776)
        {
            _S776 = true;
        }
        else
        {
            _S776 = ymin_3 >= float(image_height_3);
        }
        if(_S776)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_5);
        float3  _S790 = mean_3 - - mul_0(transpose_0(R_7), t_6);
        float3  _S791 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S791;
        float _S792 = _S790.x;
        float _S793 = _S790.y;
        float _S794 = _S790.z;
        float norm_3 = (F32_sqrt((_S792 * _S792 + _S793 * _S793 + _S794 * _S794)));
        float x_22 = _S792 / norm_3;
        float y_6 = _S793 / norm_3;
        float z_3 = _S794 / norm_3;
        float3  _S795 = _S791 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_22) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S795;
        float z2_10 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_22 * x_22 - y_6 * y_6;
        float fS1_3 = 2.0f * x_22 * y_6;
        float3  _S796 = _S795 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_22) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S796;
        float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S796 + (make_float3 (-0.59004360437393188f * (x_22 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_22) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_22 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_11, float fy_11, float cx_9, float cy_9, FixedArray<float, 10>  * dist_coeffs_18, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_12, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
        float _S797 = length_1(mean_c_4);
        bool _S798;
        if(_S797 < near_plane_4)
        {
            _S798 = true;
        }
        else
        {
            _S798 = _S797 > far_plane_4;
        }
        if(_S798)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_8 = make_float4 (fx_11, fy_11, cx_9, cy_9);
        Matrix<float, 2, 2>  covar2d_4;
        float3  _S799 = exp_0(scale_6);
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
        Matrix<float, 3, 3>  _S800 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_4;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S801 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S801 * _S799.x) * _S800.rows[0U];
        float3  _S802 = mean_4 + delta_3;
        float3  _S803 = mean_4 - delta_3;
        float3  delta_4 = make_float3 (_S801 * _S799.y) * _S800.rows[1U];
        float3  _S804 = mean_4 + delta_4;
        float3  _S805 = mean_4 - delta_4;
        float3  delta_5 = make_float3 (_S801 * _S799.z) * _S800.rows[2U];
        float3  _S806 = mean_4 + delta_5;
        float3  _S807 = mean_4 - delta_5;
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
        (&ret_1)->p_0[1U] = mul_0(R_8, _S802) + t_7;
        (&ret_1)->p_0[2U] = mul_0(R_8, _S804) + t_7;
        (&ret_1)->p_0[3U] = mul_0(R_8, _S806) + t_7;
        (&ret_1)->p_0[4U] = mul_0(R_8, _S803) + t_7;
        (&ret_1)->p_0[5U] = mul_0(R_8, _S805) + t_7;
        (&ret_1)->p_0[6U] = mul_0(R_8, _S807) + t_7;
        SigmaPoints_0 _S808 = ret_1;
        bool _S809 = fisheye_proj_3dgs_ut_0(&_S808, intrins_8, dist_coeffs_18, &covar2d_4, mean2d_12);
        if(!(true & _S809))
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S810 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S810;
        float _S811 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S811;
        float det_blur_5 = _S810 * _S811 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
            _S798 = true;
        }
        else
        {
            _S798 = xmin_4 >= float(image_width_4);
        }
        if(_S798)
        {
            _S798 = true;
        }
        else
        {
            _S798 = ymax_4 <= 0.0f;
        }
        if(_S798)
        {
            _S798 = true;
        }
        else
        {
            _S798 = ymin_4 >= float(image_height_4);
        }
        if(_S798)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_6);
        float3  _S812 = mean_4 - - mul_0(transpose_0(R_8), t_7);
        float3  _S813 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S813;
        float _S814 = _S812.x;
        float _S815 = _S812.y;
        float _S816 = _S812.z;
        float norm_4 = (F32_sqrt((_S814 * _S814 + _S815 * _S815 + _S816 * _S816)));
        float x_24 = _S814 / norm_4;
        float y_7 = _S815 / norm_4;
        float z_4 = _S816 / norm_4;
        float3  _S817 = _S813 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_24) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S817;
        float z2_12 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_24 * x_24 - y_7 * y_7;
        float fS1_4 = 2.0f * x_24 * y_7;
        float3  _S818 = _S817 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_24) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S818;
        float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S818 + (make_float3 (-0.59004360437393188f * (x_24 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_24) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_24 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_12, float fy_12, float cx_10, float cy_10, FixedArray<float, 10>  * dist_coeffs_19, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_13, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S819 = exp_0(scale_7);
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
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S819.x, 0.0f, 0.0f, 0.0f, _S819.y, 0.0f, 0.0f, 0.0f, _S819.z));
    Matrix<float, 3, 3>  _S820 = transpose_0(R_9);
    float _S821 = float(image_width_5);
    float _S822 = float(image_height_5);
    float _S823 = 0.30000001192092896f * (0.5f * _S821 / fx_12);
    float _S824 = 0.30000001192092896f * (0.5f * _S822 / fy_12);
    float rz_2 = 1.0f / mean_c_5.z;
    float rz2_2 = rz_2 * rz_2;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_12 * rz_2, 0.0f, - fx_12 * (mean_c_5.z * (F32_min(((_S821 - cx_10) / fx_12 + _S823), ((F32_max((- (cx_10 / fx_12 + _S823)), (mean_c_5.x * rz_2))))))) * rz2_2, 0.0f, fy_12 * rz_2, - fy_12 * (mean_c_5.z * (F32_min(((_S822 - cy_10) / fy_12 + _S824), ((F32_max((- (cy_10 / fy_12 + _S824)), (mean_c_5.y * rz_2))))))) * rz2_2);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_8, mul_4(mul_4(R_9, mul_4(M_5, transpose_0(M_5))), _S820)), transpose_1(J_8));
    *mean2d_13 = make_float2 (fx_12 * mean_c_5.x * rz_2 + cx_10, fy_12 * mean_c_5.y * rz_2 + cy_10);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S825 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S825;
    float _S826 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S826;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S825 * _S826 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S827 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
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
    *conic_5 = make_float3 (_S827.rows[int(0)].x, _S827.rows[int(0)].y, _S827.rows[int(1)].y);
    float3  _S828 = mean_5 - - mul_0(_S820, t_8);
    float3  _S829 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S829;
    float _S830 = _S828.x;
    float _S831 = _S828.y;
    float _S832 = _S828.z;
    float norm_5 = (F32_sqrt((_S830 * _S830 + _S831 * _S831 + _S832 * _S832)));
    float x_26 = _S830 / norm_5;
    float y_8 = _S831 / norm_5;
    float z_5 = _S832 / norm_5;
    float3  _S833 = _S829 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_26) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S833;
    float z2_14 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_26 * x_26 - y_8 * y_8;
    float fS1_5 = 2.0f * x_26 * y_8;
    float3  _S834 = _S833 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_26) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S834;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S834 + (make_float3 (-0.59004360437393188f * (x_26 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_26) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_26 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_13, float fy_13, float cx_11, float cy_11, FixedArray<float, 10>  * dist_coeffs_20, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_14, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_6) + t_9;
    float3  _S835 = exp_0(scale_8);
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
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S835.x, 0.0f, 0.0f, 0.0f, _S835.y, 0.0f, 0.0f, 0.0f, _S835.z));
    Matrix<float, 3, 3>  _S836 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    bool _S837 = fisheye_proj_3dgs_1(mean_c_6, mul_4(mul_4(R_10, mul_4(M_6, transpose_0(M_6))), _S836), make_float4 (fx_13, fy_13, cx_11, cy_11), dist_coeffs_20, &covar2d_6, mean2d_14);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S838 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S838;
    float _S839 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S839;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S838 * _S839 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S840 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
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
    *conic_6 = make_float3 (_S840.rows[int(0)].x, _S840.rows[int(0)].y, _S840.rows[int(1)].y);
    float3  _S841 = mean_6 - - mul_0(_S836, t_9);
    float3  _S842 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S842;
    float _S843 = _S841.x;
    float _S844 = _S841.y;
    float _S845 = _S841.z;
    float norm_6 = (F32_sqrt((_S843 * _S843 + _S844 * _S844 + _S845 * _S845)));
    float x_28 = _S843 / norm_6;
    float y_9 = _S844 / norm_6;
    float z_6 = _S845 / norm_6;
    float3  _S846 = _S842 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_28) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S846;
    float z2_16 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_28 * x_28 - y_9 * y_9;
    float fS1_6 = 2.0f * x_28 * y_9;
    float3  _S847 = _S846 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_16 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_28) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S847;
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S847 + (make_float3 (-0.59004360437393188f * (x_28 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_16 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_28) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_28 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_14, float fy_14, float cx_12, float cy_12, FixedArray<float, 10>  * dist_coeffs_21, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_15, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_7) + t_10;
    float3  _S848 = exp_0(scale_9);
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
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S848.x, 0.0f, 0.0f, 0.0f, _S848.y, 0.0f, 0.0f, 0.0f, _S848.z));
    Matrix<float, 3, 3>  _S849 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_14, 0.0f, 0.0f, 0.0f, fy_14, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_9, mul_4(mul_4(R_11, mul_4(M_7, transpose_0(M_7))), _S849)), transpose_1(J_9));
    *mean2d_15 = make_float2 (fx_14 * mean_c_7.x + cx_12, fy_14 * mean_c_7.y + cy_12);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S850 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S850;
    float _S851 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S851;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S850 * _S851 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S852 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
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
    *conic_7 = make_float3 (_S852.rows[int(0)].x, _S852.rows[int(0)].y, _S852.rows[int(1)].y);
    float3  _S853 = mean_7 - - mul_0(_S849, t_10);
    float3  _S854 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S854;
    float _S855 = _S853.x;
    float _S856 = _S853.y;
    float _S857 = _S853.z;
    float norm_7 = (F32_sqrt((_S855 * _S855 + _S856 * _S856 + _S857 * _S857)));
    float x_30 = _S855 / norm_7;
    float y_10 = _S856 / norm_7;
    float z_7 = _S857 / norm_7;
    float3  _S858 = _S854 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_30) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S858;
    float z2_18 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_30 * x_30 - y_10 * y_10;
    float fS1_7 = 2.0f * x_30 * y_10;
    float3  _S859 = _S858 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_18 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_30) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S859;
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S859 + (make_float3 (-0.59004360437393188f * (x_30 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_18 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_30) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_30 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_15, float fy_15, float cx_13, float cy_13, FixedArray<float, 10>  * dist_coeffs_22, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_16, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_8) + t_11;
    float4  intrins_9 = make_float4 (fx_15, fy_15, cx_13, cy_13);
    float3  _S860 = exp_0(scale_10);
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
    Matrix<float, 3, 3>  _S861 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))));
    SigmaPoints_0 ret_2;
    (&ret_2)->p_0[int(0)] = mean_8;
    (&ret_2)->w_mean_0[int(0)] = 0.0f;
    (&ret_2)->w_cov_0[int(0)] = 2.0f;
    float _S862 = (F32_sqrt((3.0f)));
    float3  delta_6 = make_float3 (_S862 * _S860.x) * _S861.rows[0U];
    float3  _S863 = mean_8 + delta_6;
    float3  _S864 = mean_8 - delta_6;
    float3  delta_7 = make_float3 (_S862 * _S860.y) * _S861.rows[1U];
    float3  _S865 = mean_8 + delta_7;
    float3  _S866 = mean_8 - delta_7;
    float3  delta_8 = make_float3 (_S862 * _S860.z) * _S861.rows[2U];
    float3  _S867 = mean_8 + delta_8;
    float3  _S868 = mean_8 - delta_8;
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
    (&ret_2)->p_0[1U] = mul_0(R_12, _S863) + t_11;
    (&ret_2)->p_0[2U] = mul_0(R_12, _S865) + t_11;
    (&ret_2)->p_0[3U] = mul_0(R_12, _S867) + t_11;
    (&ret_2)->p_0[4U] = mul_0(R_12, _S864) + t_11;
    (&ret_2)->p_0[5U] = mul_0(R_12, _S866) + t_11;
    (&ret_2)->p_0[6U] = mul_0(R_12, _S868) + t_11;
    SigmaPoints_0 _S869 = ret_2;
    Matrix<float, 2, 2>  covar2d_8;
    bool _S870 = persp_proj_3dgs_ut_1(&_S869, intrins_9, dist_coeffs_22, image_width_8, image_height_8, &covar2d_8, mean2d_16);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S871 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S871;
    float _S872 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S872;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S871 * _S872 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
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
    float3  _S873 = mean_8 - - mul_0(transpose_0(R_12), t_11);
    float3  _S874 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S874;
    float _S875 = _S873.x;
    float _S876 = _S873.y;
    float _S877 = _S873.z;
    float norm_8 = (F32_sqrt((_S875 * _S875 + _S876 * _S876 + _S877 * _S877)));
    float x_32 = _S875 / norm_8;
    float y_11 = _S876 / norm_8;
    float z_8 = _S877 / norm_8;
    float3  _S878 = _S874 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_32) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S878;
    float z2_20 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_32 * x_32 - y_11 * y_11;
    float fS1_8 = 2.0f * x_32 * y_11;
    float3  _S879 = _S878 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_20 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_32) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S879;
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S879 + (make_float3 (-0.59004360437393188f * (x_32 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_20 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_32) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_32 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_16, float fy_16, float cx_14, float cy_14, FixedArray<float, 10>  * dist_coeffs_23, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_17, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_9) + t_12;
    float4  intrins_10 = make_float4 (fx_16, fy_16, cx_14, cy_14);
    float3  _S880 = exp_0(scale_11);
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
    Matrix<float, 3, 3>  _S881 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))));
    SigmaPoints_0 ret_3;
    (&ret_3)->p_0[int(0)] = mean_9;
    (&ret_3)->w_mean_0[int(0)] = 0.0f;
    (&ret_3)->w_cov_0[int(0)] = 2.0f;
    float _S882 = (F32_sqrt((3.0f)));
    float3  delta_9 = make_float3 (_S882 * _S880.x) * _S881.rows[0U];
    float3  _S883 = mean_9 + delta_9;
    float3  _S884 = mean_9 - delta_9;
    float3  delta_10 = make_float3 (_S882 * _S880.y) * _S881.rows[1U];
    float3  _S885 = mean_9 + delta_10;
    float3  _S886 = mean_9 - delta_10;
    float3  delta_11 = make_float3 (_S882 * _S880.z) * _S881.rows[2U];
    float3  _S887 = mean_9 + delta_11;
    float3  _S888 = mean_9 - delta_11;
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
    (&ret_3)->p_0[1U] = mul_0(R_13, _S883) + t_12;
    (&ret_3)->p_0[2U] = mul_0(R_13, _S885) + t_12;
    (&ret_3)->p_0[3U] = mul_0(R_13, _S887) + t_12;
    (&ret_3)->p_0[4U] = mul_0(R_13, _S884) + t_12;
    (&ret_3)->p_0[5U] = mul_0(R_13, _S886) + t_12;
    (&ret_3)->p_0[6U] = mul_0(R_13, _S888) + t_12;
    SigmaPoints_0 _S889 = ret_3;
    Matrix<float, 2, 2>  covar2d_9;
    bool _S890 = fisheye_proj_3dgs_ut_1(&_S889, intrins_10, dist_coeffs_23, &covar2d_9, mean2d_17);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S891 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S891;
    float _S892 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S892;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S891 * _S892 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
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
    float3  _S893 = mean_9 - - mul_0(transpose_0(R_13), t_12);
    float3  _S894 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S894;
    float _S895 = _S893.x;
    float _S896 = _S893.y;
    float _S897 = _S893.z;
    float norm_9 = (F32_sqrt((_S895 * _S895 + _S896 * _S896 + _S897 * _S897)));
    float x_34 = _S895 / norm_9;
    float y_12 = _S896 / norm_9;
    float z_9 = _S897 / norm_9;
    float3  _S898 = _S894 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_34) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S898;
    float z2_22 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_34 * x_34 - y_12 * y_12;
    float fS1_9 = 2.0f * x_34 * y_12;
    float3  _S899 = _S898 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_34) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S899;
    float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S899 + (make_float3 (-0.59004360437393188f * (x_34 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_34) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_34 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S900, float3  _S901)
{
    return mul_0(_S900, _S901);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S902)
{
    return exp_0(_S902);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S903, Matrix<float, 3, 3>  _S904)
{
    return mul_4(_S903, _S904);
}

inline __device__ float s_primal_ctx_max_0(float _S905, float _S906)
{
    return (F32_max((_S905), (_S906)));
}

inline __device__ float s_primal_ctx_min_0(float _S907, float _S908)
{
    return (F32_min((_S907), (_S908)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S909, Matrix<float, 3, 3>  _S910)
{
    return mul_5(_S909, _S910);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S911, Matrix<float, 3, 2>  _S912)
{
    return mul_6(_S911, _S912);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S913)
{
    return (F32_sqrt((_S913)));
}

inline __device__ float s_primal_ctx_exp_1(float _S914)
{
    return (F32_exp((_S914)));
}

inline __device__ float s_primal_ctx_log_0(float _S915)
{
    return (F32_log((_S915)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S916, float3  _S917)
{
    return dot_0(_S916, _S917);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S918, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S919, float3  _S920)
{
    _d_max_vector_0(_S918, _S919, _S920);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S921, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S922, float3  _S923)
{
    _d_mul_0(_S921, _S922, _S923);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S924, float _S925)
{
    _d_log_0(_S924, _S925);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S926, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S927, float _S928)
{
    _d_dot_0(_S926, _S927, _S928);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S929, DiffPair_float_0 * _S930, float _S931)
{
    _d_min_0(_S929, _S930, _S931);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S932, float _S933)
{
    _d_exp_0(_S932, _S933);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S934, DiffPair_float_0 * _S935, float _S936)
{
    _d_max_0(_S934, _S935, _S936);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S937, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S938, Matrix<float, 2, 2>  _S939)
{
    mul_3(_S937, _S938, _S939);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S940, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S941, Matrix<float, 2, 3>  _S942)
{
    mul_2(_S940, _S941, _S942);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S943, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S944, Matrix<float, 3, 3>  _S945)
{
    mul_1(_S943, _S944, _S945);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S946, float3  _S947)
{
    _d_exp_vector_0(_S946, _S947);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_17, float fy_17, float cx_15, float cy_15, FixedArray<float, 10>  * dist_coeffs_24, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S948 = s_primal_ctx_exp_0(scale_12);
    float _S949 = quat_13.y;
    float x2_13 = _S949 * _S949;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S950 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S948.x, 0.0f, 0.0f, 0.0f, _S948.y, 0.0f, 0.0f, 0.0f, _S948.z);
    Matrix<float, 3, 3>  _S951 = s_primal_ctx_mul_2(_S950, S_0);
    Matrix<float, 3, 3>  _S952 = transpose_0(_S951);
    Matrix<float, 3, 3>  _S953 = s_primal_ctx_mul_2(_S951, _S952);
    Matrix<float, 3, 3>  _S954 = s_primal_ctx_mul_2(R_14, _S953);
    Matrix<float, 3, 3>  _S955 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S956 = s_primal_ctx_mul_2(_S954, _S955);
    float _S957 = float(image_width_10);
    float _S958 = float(image_height_10);
    float _S959 = 0.30000001192092896f * (0.5f * _S957 / fx_17);
    float lim_x_pos_2 = (_S957 - cx_15) / fx_17 + _S959;
    float _S960 = 0.30000001192092896f * (0.5f * _S958 / fy_17);
    float lim_y_pos_2 = (_S958 - cy_15) / fy_17 + _S960;
    float rz_3 = 1.0f / mean_c_10.z;
    float _S961 = mean_c_10.z * mean_c_10.z;
    float rz2_3 = rz_3 * rz_3;
    float _S962 = - (cx_15 / fx_17 + _S959);
    float _S963 = mean_c_10.x * rz_3;
    float _S964 = s_primal_ctx_max_0(_S962, _S963);
    float _S965 = s_primal_ctx_min_0(lim_x_pos_2, _S964);
    float _S966 = - (cy_15 / fy_17 + _S960);
    float _S967 = mean_c_10.y * rz_3;
    float _S968 = s_primal_ctx_max_0(_S966, _S967);
    float _S969 = s_primal_ctx_min_0(lim_y_pos_2, _S968);
    float _S970 = - fx_17;
    float _S971 = _S970 * (mean_c_10.z * _S965);
    float _S972 = - fy_17;
    float _S973 = _S972 * (mean_c_10.z * _S969);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (fx_17 * rz_3, 0.0f, _S971 * rz2_3, 0.0f, fy_17 * rz_3, _S973 * rz2_3);
    Matrix<float, 2, 3>  _S974 = s_primal_ctx_mul_3(J_10, _S956);
    Matrix<float, 3, 2>  _S975 = transpose_1(J_10);
    Matrix<float, 2, 2>  _S976 = s_primal_ctx_mul_4(_S974, _S975);
    float _S977 = fx_17 * mean_c_10.x;
    float _S978 = fy_17 * mean_c_10.y;
    float _S979 = _S976.rows[int(0)].y * _S976.rows[int(1)].x;
    float det_orig_11 = _S976.rows[int(0)].x * _S976.rows[int(1)].y - _S979;
    float _S980 = _S976.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S981 = _S976;
    *&(((&_S981)->rows + (int(0)))->x) = _S980;
    float _S982 = _S976.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S981)->rows + (int(1)))->y) = _S982;
    Matrix<float, 2, 2>  _S983 = _S981;
    Matrix<float, 2, 2>  _S984 = _S981;
    float det_blur_6 = _S980 * _S982 - _S979;
    float _S985 = det_orig_11 / det_blur_6;
    float _S986 = det_blur_6 * det_blur_6;
    float _S987 = s_primal_ctx_max_0(0.0f, _S985);
    float _S988 = s_primal_ctx_sqrt_0(_S987);
    float invdet_7 = 1.0f / det_blur_6;
    float _S989 = - _S976.rows[int(0)].y;
    float _S990 = - _S976.rows[int(1)].x;
    float _S991 = - in_opacity_10;
    float _S992 = 1.0f + s_primal_ctx_exp_1(_S991);
    float _S993 = 1.0f / _S992;
    float _S994 = _S992 * _S992;
    float _S995;
    if(antialiased_10)
    {
        _S995 = _S993 * _S988;
    }
    else
    {
        _S995 = _S993;
    }
    float _S996 = _S995 / 0.00392156885936856f;
    float _S997 = 2.0f * s_primal_ctx_log_0(_S996);
    float _S998 = s_primal_ctx_sqrt_0(_S997);
    float _S999 = _S983.rows[int(0)].x;
    float _S1000 = _S984.rows[int(1)].y;
    float _S1001 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S1002 = mean_10 - - s_primal_ctx_mul_1(_S955, t_13);
    float _S1003 = _S1002.x;
    float _S1004 = _S1002.y;
    float _S1005 = _S1002.z;
    float _S1006 = _S1003 * _S1003 + _S1004 * _S1004 + _S1005 * _S1005;
    float _S1007 = s_primal_ctx_sqrt_0(_S1006);
    float x_35 = _S1003 / _S1007;
    float3  _S1008 = make_float3 (x_35);
    float _S1009 = _S1007 * _S1007;
    float y_13 = _S1004 / _S1007;
    float z_10 = _S1005 / _S1007;
    float3  _S1010 = make_float3 (z_10);
    float _S1011 = - y_13;
    float3  _S1012 = make_float3 (_S1011);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_35 * x_35 - y_13 * y_13;
    float _S1013 = 2.0f * x_35;
    float fS1_10 = _S1013 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S1014 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_35;
    float3  _S1015 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S1016 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S1017 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S1018 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S1019 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S1019;
    float3  _S1020 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_35;
    float3  _S1021 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S1022 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S1023 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S1024 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_35 * fC1_10 - y_13 * fS1_10);
    float3  _S1025 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_35 * fS1_10 + y_13 * fC1_10);
    float3  _S1026 = make_float3 (pSH9_0);
    float3  _S1027 = make_float3 (0.0f);
    float3  _S1028 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1029;
    (&_S1029)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1011) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_35) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S1029)->differential_0 = _S1028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1030;
    (&_S1030)->primal_0 = _S1027;
    (&_S1030)->differential_0 = _S1028;
    s_bwd_prop_max_0(&_S1029, &_S1030, v_rgb_0);
    float3  _S1031 = _S1025 * _S1029.differential_0;
    float3  _S1032 = (*sh_coeffs_10)[int(15)] * _S1029.differential_0;
    float3  _S1033 = _S1023 * _S1029.differential_0;
    float3  _S1034 = (*sh_coeffs_10)[int(14)] * _S1029.differential_0;
    float3  _S1035 = _S1021 * _S1029.differential_0;
    float3  _S1036 = (*sh_coeffs_10)[int(13)] * _S1029.differential_0;
    float3  _S1037 = _S1020 * _S1029.differential_0;
    float3  _S1038 = (*sh_coeffs_10)[int(12)] * _S1029.differential_0;
    float3  _S1039 = _S1022 * _S1029.differential_0;
    float3  _S1040 = (*sh_coeffs_10)[int(11)] * _S1029.differential_0;
    float3  _S1041 = _S1024 * _S1029.differential_0;
    float3  _S1042 = (*sh_coeffs_10)[int(10)] * _S1029.differential_0;
    float3  _S1043 = _S1026 * _S1029.differential_0;
    float3  _S1044 = (*sh_coeffs_10)[int(9)] * _S1029.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S1044.x + _S1044.y + _S1044.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S1032.x + _S1032.y + _S1032.z);
    float _S1045 = _S1042.x + _S1042.y + _S1042.z;
    float _S1046 = _S1034.x + _S1034.y + _S1034.z;
    float _S1047 = _S1040.x + _S1040.y + _S1040.z;
    float _S1048 = _S1036.x + _S1036.y + _S1036.z;
    float _S1049 = _S1038.x + _S1038.y + _S1038.z;
    float _S1050 = - s_diff_fC2_T_0;
    float3  _S1051 = _S1017 * _S1029.differential_0;
    float3  _S1052 = (*sh_coeffs_10)[int(8)] * _S1029.differential_0;
    float3  _S1053 = _S1015 * _S1029.differential_0;
    float3  _S1054 = (*sh_coeffs_10)[int(7)] * _S1029.differential_0;
    float3  _S1055 = _S1014 * _S1029.differential_0;
    float3  _S1056 = (*sh_coeffs_10)[int(6)] * _S1029.differential_0;
    float3  _S1057 = _S1016 * _S1029.differential_0;
    float3  _S1058 = (*sh_coeffs_10)[int(5)] * _S1029.differential_0;
    float3  _S1059 = _S1018 * _S1029.differential_0;
    float3  _S1060 = (*sh_coeffs_10)[int(4)] * _S1029.differential_0;
    float _S1061 = _S1058.x + _S1058.y + _S1058.z;
    float _S1062 = _S1054.x + _S1054.y + _S1054.z;
    float _S1063 = fTmp1B_10 * _S1045 + x_35 * s_diff_fS2_T_0 + y_13 * _S1050 + 0.54627424478530884f * (_S1060.x + _S1060.y + _S1060.z);
    float _S1064 = fTmp1B_10 * _S1046 + y_13 * s_diff_fS2_T_0 + x_35 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S1052.x + _S1052.y + _S1052.z);
    float _S1065 = y_13 * - _S1064;
    float _S1066 = x_35 * _S1064;
    float _S1067 = z_10 * (1.86588168144226074f * (z_10 * _S1049) + -2.28522896766662598f * (y_13 * _S1047 + x_35 * _S1048) + 0.94617468118667603f * (_S1056.x + _S1056.y + _S1056.z));
    float3  _S1068 = make_float3 (0.48860251903533936f) * _S1029.differential_0;
    float3  _S1069 = - _S1068;
    float3  _S1070 = _S1008 * _S1069;
    float3  _S1071 = (*sh_coeffs_10)[int(3)] * _S1069;
    float3  _S1072 = _S1010 * _S1068;
    float3  _S1073 = (*sh_coeffs_10)[int(2)] * _S1068;
    float3  _S1074 = _S1012 * _S1068;
    float3  _S1075 = (*sh_coeffs_10)[int(1)] * _S1068;
    float _S1076 = (_S1019 * _S1049 + 1.44530570507049561f * (fS1_10 * _S1045 + fC1_10 * _S1046) + -1.09254848957061768f * (y_13 * _S1061 + x_35 * _S1062) + _S1067 + _S1067 + _S1073.x + _S1073.y + _S1073.z) / _S1009;
    float _S1077 = _S1007 * _S1076;
    float _S1078 = (fTmp0C_10 * _S1047 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S1050 + fTmp0B_10 * _S1061 + _S1013 * _S1063 + _S1065 + _S1065 + - (_S1075.x + _S1075.y + _S1075.z)) / _S1009;
    float _S1079 = _S1007 * _S1078;
    float _S1080 = (fTmp0C_10 * _S1048 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S1062 + 2.0f * (y_13 * _S1063) + _S1066 + _S1066 + _S1071.x + _S1071.y + _S1071.z) / _S1009;
    float _S1081 = _S1007 * _S1080;
    float _S1082 = _S1005 * - _S1076 + _S1004 * - _S1078 + _S1003 * - _S1080;
    DiffPair_float_0 _S1083;
    (&_S1083)->primal_0 = _S1006;
    (&_S1083)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1083, _S1082);
    float _S1084 = _S1005 * _S1083.differential_0;
    float _S1085 = _S1004 * _S1083.differential_0;
    float _S1086 = _S1003 * _S1083.differential_0;
    float3  _S1087 = make_float3 (0.282094806432724f) * _S1029.differential_0;
    float3  _S1088 = make_float3 (_S1081 + _S1086 + _S1086, _S1079 + _S1085 + _S1085, _S1077 + _S1084 + _S1084);
    float3  _S1089 = - - _S1088;
    Matrix<float, 3, 3>  _S1090 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1091;
    (&_S1091)->primal_0 = _S955;
    (&_S1091)->differential_0 = _S1090;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1092;
    (&_S1092)->primal_0 = t_13;
    (&_S1092)->differential_0 = _S1028;
    s_bwd_prop_mul_1(&_S1091, &_S1092, _S1089);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1093 = _S1091;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1094 = _S1092;
    float2  _S1095 = make_float2 (0.0f);
    float2  _S1096 = _S1095;
    *&((&_S1096)->y) = v_conic_0.z;
    float2  _S1097 = _S1095;
    *&((&_S1097)->y) = v_conic_0.y;
    *&((&_S1097)->x) = v_conic_0.x;
    float _S1098 = 0.5f * v_depth_0;
    DiffPair_float_0 _S1099;
    (&_S1099)->primal_0 = _S1001;
    (&_S1099)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1099, _S1098);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1100;
    (&_S1100)->primal_0 = mean_c_10;
    (&_S1100)->differential_0 = _S1028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1101;
    (&_S1101)->primal_0 = mean_c_10;
    (&_S1101)->differential_0 = _S1028;
    s_bwd_prop_dot_0(&_S1100, &_S1101, _S1099.differential_0);
    DiffPair_float_0 _S1102;
    (&_S1102)->primal_0 = _S1000;
    (&_S1102)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1102, 0.0f);
    DiffPair_float_0 _S1103;
    (&_S1103)->primal_0 = _S999;
    (&_S1103)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1103, 0.0f);
    DiffPair_float_0 _S1104;
    (&_S1104)->primal_0 = 3.32999992370605469f;
    (&_S1104)->differential_0 = 0.0f;
    DiffPair_float_0 _S1105;
    (&_S1105)->primal_0 = _S998;
    (&_S1105)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1104, &_S1105, 0.0f);
    DiffPair_float_0 _S1106;
    (&_S1106)->primal_0 = _S997;
    (&_S1106)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1106, _S1105.differential_0);
    float _S1107 = 2.0f * _S1106.differential_0;
    DiffPair_float_0 _S1108;
    (&_S1108)->primal_0 = _S996;
    (&_S1108)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1108, _S1107);
    float _S1109 = v_opacity_0 + 254.9999847412109375f * _S1108.differential_0;
    Matrix<float, 2, 2>  _S1110 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1111 = _S1110;
    _S1111[int(1)] = _S1096;
    _S1111[int(0)] = _S1097;
    Matrix<float, 2, 2>  _S1112 = _S1111;
    FixedArray<float3 , 16>  _S1113;
    _S1113[int(0)] = _S1028;
    _S1113[int(1)] = _S1028;
    _S1113[int(2)] = _S1028;
    _S1113[int(3)] = _S1028;
    _S1113[int(4)] = _S1028;
    _S1113[int(5)] = _S1028;
    _S1113[int(6)] = _S1028;
    _S1113[int(7)] = _S1028;
    _S1113[int(8)] = _S1028;
    _S1113[int(9)] = _S1028;
    _S1113[int(10)] = _S1028;
    _S1113[int(11)] = _S1028;
    _S1113[int(12)] = _S1028;
    _S1113[int(13)] = _S1028;
    _S1113[int(14)] = _S1028;
    _S1113[int(15)] = _S1028;
    _S1113[int(7)] = _S1053;
    _S1113[int(0)] = _S1087;
    _S1113[int(1)] = _S1074;
    _S1113[int(2)] = _S1072;
    _S1113[int(3)] = _S1070;
    _S1113[int(4)] = _S1059;
    _S1113[int(5)] = _S1057;
    _S1113[int(6)] = _S1055;
    _S1113[int(15)] = _S1031;
    _S1113[int(8)] = _S1051;
    _S1113[int(9)] = _S1043;
    _S1113[int(10)] = _S1041;
    _S1113[int(11)] = _S1039;
    _S1113[int(12)] = _S1037;
    _S1113[int(13)] = _S1035;
    _S1113[int(14)] = _S1033;
    float3  _S1114 = _S1113[int(0)];
    float3  _S1115 = _S1113[int(1)];
    float3  _S1116 = _S1113[int(2)];
    float3  _S1117 = _S1113[int(3)];
    float3  _S1118 = _S1113[int(4)];
    float3  _S1119 = _S1113[int(5)];
    float3  _S1120 = _S1113[int(6)];
    float3  _S1121 = _S1113[int(7)];
    float3  _S1122 = _S1113[int(8)];
    float3  _S1123 = _S1113[int(9)];
    float3  _S1124 = _S1113[int(10)];
    float3  _S1125 = _S1113[int(11)];
    float3  _S1126 = _S1113[int(12)];
    float3  _S1127 = _S1113[int(13)];
    float3  _S1128 = _S1113[int(14)];
    float3  _S1129 = _S1113[int(15)];
    float3  _S1130 = _S1101.differential_0 + _S1100.differential_0;
    float2  _S1131 = make_float2 (0.0f, _S1102.differential_0);
    float2  _S1132 = make_float2 (_S1103.differential_0, 0.0f);
    float _S1133;
    if(antialiased_10)
    {
        float _S1134 = _S993 * _S1109;
        _S995 = _S988 * _S1109;
        _S1133 = _S1134;
    }
    else
    {
        _S995 = _S1109;
        _S1133 = 0.0f;
    }
    float _S1135 = - (_S995 / _S994);
    DiffPair_float_0 _S1136;
    (&_S1136)->primal_0 = _S991;
    (&_S1136)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1136, _S1135);
    float _S1137 = - _S1136.differential_0;
    float _S1138 = invdet_7 * _S1112.rows[int(1)].y;
    float _S1139 = - (invdet_7 * _S1112.rows[int(1)].x);
    float _S1140 = - (invdet_7 * _S1112.rows[int(0)].y);
    float _S1141 = invdet_7 * _S1112.rows[int(0)].x;
    float _S1142 = - ((_S980 * _S1112.rows[int(1)].y + _S990 * _S1112.rows[int(1)].x + _S989 * _S1112.rows[int(0)].y + _S982 * _S1112.rows[int(0)].x) / _S986);
    DiffPair_float_0 _S1143;
    (&_S1143)->primal_0 = _S987;
    (&_S1143)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1143, _S1133);
    DiffPair_float_0 _S1144;
    (&_S1144)->primal_0 = 0.0f;
    (&_S1144)->differential_0 = 0.0f;
    DiffPair_float_0 _S1145;
    (&_S1145)->primal_0 = _S985;
    (&_S1145)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1144, &_S1145, _S1143.differential_0);
    float _S1146 = _S1145.differential_0 / _S986;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1146;
    float _S1147 = _S1142 + det_orig_11 * - _S1146;
    float _S1148 = - _S1147;
    float _S1149 = _S980 * _S1147;
    float _S1150 = _S982 * _S1147;
    Matrix<float, 2, 2>  _S1151 = _S1110;
    _S1151[int(1)] = _S1131;
    _S1151[int(0)] = _S1132;
    _S981 = _S1151;
    *&(((&_S981)->rows + (int(1)))->y) = 0.0f;
    float _S1152 = _S1141 + _S1149 + _S1151.rows[int(1)].y;
    *&(((&_S981)->rows + (int(0)))->x) = 0.0f;
    float _S1153 = _S1138 + _S1150 + _S1151.rows[int(0)].x;
    float _S1154 = _S1148 + - s_diff_det_orig_T_0;
    float _S1155 = _S1139 + _S976.rows[int(0)].y * _S1154;
    float _S1156 = _S1140 + _S976.rows[int(1)].x * _S1154;
    float _S1157 = _S976.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1158 = _S1152 + _S976.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1159 = _S1095;
    *&((&_S1159)->x) = _S1155;
    *&((&_S1159)->y) = _S1158;
    float _S1160 = _S1153 + _S1157;
    float2  _S1161 = _S1095;
    *&((&_S1161)->y) = _S1156;
    *&((&_S1161)->x) = _S1160;
    float _S1162 = _S978 * v_mean2d_0.y;
    float _S1163 = fy_17 * (rz_3 * v_mean2d_0.y);
    float _S1164 = _S977 * v_mean2d_0.x;
    float _S1165 = fx_17 * (rz_3 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S1166 = _S1110;
    _S1166[int(1)] = _S1159;
    _S1166[int(0)] = _S1161;
    Matrix<float, 2, 2>  _S1167 = _S981 + _S1166;
    Matrix<float, 2, 3>  _S1168 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1169;
    (&_S1169)->primal_0 = _S974;
    (&_S1169)->differential_0 = _S1168;
    Matrix<float, 3, 2>  _S1170 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1171;
    (&_S1171)->primal_0 = _S975;
    (&_S1171)->differential_0 = _S1170;
    s_bwd_prop_mul_2(&_S1169, &_S1171, _S1167);
    Matrix<float, 2, 3>  _S1172 = transpose_2(_S1171.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1173;
    (&_S1173)->primal_0 = J_10;
    (&_S1173)->differential_0 = _S1168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1174;
    (&_S1174)->primal_0 = _S956;
    (&_S1174)->differential_0 = _S1090;
    s_bwd_prop_mul_3(&_S1173, &_S1174, _S1169.differential_0);
    Matrix<float, 2, 3>  _S1175 = _S1172 + _S1173.differential_0;
    float _S1176 = _S973 * _S1175.rows[int(1)].z;
    float s_diff_ty_T_0 = _S972 * (rz2_3 * _S1175.rows[int(1)].z);
    float _S1177 = fy_17 * _S1175.rows[int(1)].y;
    float _S1178 = _S971 * _S1175.rows[int(0)].z;
    float s_diff_tx_T_0 = _S970 * (rz2_3 * _S1175.rows[int(0)].z);
    float _S1179 = fx_17 * _S1175.rows[int(0)].x;
    float _S1180 = mean_c_10.z * s_diff_ty_T_0;
    float _S1181 = _S969 * s_diff_ty_T_0;
    DiffPair_float_0 _S1182;
    (&_S1182)->primal_0 = lim_y_pos_2;
    (&_S1182)->differential_0 = 0.0f;
    DiffPair_float_0 _S1183;
    (&_S1183)->primal_0 = _S968;
    (&_S1183)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1182, &_S1183, _S1180);
    DiffPair_float_0 _S1184;
    (&_S1184)->primal_0 = _S966;
    (&_S1184)->differential_0 = 0.0f;
    DiffPair_float_0 _S1185;
    (&_S1185)->primal_0 = _S967;
    (&_S1185)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1184, &_S1185, _S1183.differential_0);
    float _S1186 = mean_c_10.y * _S1185.differential_0;
    float _S1187 = rz_3 * _S1185.differential_0;
    float _S1188 = mean_c_10.z * s_diff_tx_T_0;
    float _S1189 = _S965 * s_diff_tx_T_0;
    DiffPair_float_0 _S1190;
    (&_S1190)->primal_0 = lim_x_pos_2;
    (&_S1190)->differential_0 = 0.0f;
    DiffPair_float_0 _S1191;
    (&_S1191)->primal_0 = _S964;
    (&_S1191)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1190, &_S1191, _S1188);
    DiffPair_float_0 _S1192;
    (&_S1192)->primal_0 = _S962;
    (&_S1192)->differential_0 = 0.0f;
    DiffPair_float_0 _S1193;
    (&_S1193)->primal_0 = _S963;
    (&_S1193)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1192, &_S1193, _S1191.differential_0);
    float _S1194 = rz_3 * (_S1176 + _S1178);
    float _S1195 = _S1181 + _S1189 + - ((_S1162 + _S1164 + _S1177 + _S1179 + _S1186 + mean_c_10.x * _S1193.differential_0 + _S1194 + _S1194) / _S961);
    float _S1196 = _S1163 + _S1187;
    float _S1197 = _S1165 + rz_3 * _S1193.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1198;
    (&_S1198)->primal_0 = _S954;
    (&_S1198)->differential_0 = _S1090;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1199;
    (&_S1199)->primal_0 = _S955;
    (&_S1199)->differential_0 = _S1090;
    s_bwd_prop_mul_4(&_S1198, &_S1199, _S1174.differential_0);
    Matrix<float, 3, 3>  _S1200 = transpose_0(_S1199.differential_0 + _S1093.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1201;
    (&_S1201)->primal_0 = R_14;
    (&_S1201)->differential_0 = _S1090;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1202;
    (&_S1202)->primal_0 = _S953;
    (&_S1202)->differential_0 = _S1090;
    s_bwd_prop_mul_4(&_S1201, &_S1202, _S1198.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1203;
    (&_S1203)->primal_0 = _S951;
    (&_S1203)->differential_0 = _S1090;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1204;
    (&_S1204)->primal_0 = _S952;
    (&_S1204)->differential_0 = _S1090;
    s_bwd_prop_mul_4(&_S1203, &_S1204, _S1202.differential_0);
    Matrix<float, 3, 3>  _S1205 = _S1203.differential_0 + transpose_0(_S1204.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1206;
    (&_S1206)->primal_0 = _S950;
    (&_S1206)->differential_0 = _S1090;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1207;
    (&_S1207)->primal_0 = S_0;
    (&_S1207)->differential_0 = _S1090;
    s_bwd_prop_mul_4(&_S1206, &_S1207, _S1205);
    Matrix<float, 3, 3>  _S1208 = transpose_0(_S1206.differential_0);
    float _S1209 = 2.0f * - _S1208.rows[int(2)].z;
    float _S1210 = 2.0f * _S1208.rows[int(2)].y;
    float _S1211 = 2.0f * _S1208.rows[int(2)].x;
    float _S1212 = 2.0f * _S1208.rows[int(1)].z;
    float _S1213 = 2.0f * - _S1208.rows[int(1)].y;
    float _S1214 = 2.0f * _S1208.rows[int(1)].x;
    float _S1215 = 2.0f * _S1208.rows[int(0)].z;
    float _S1216 = 2.0f * _S1208.rows[int(0)].y;
    float _S1217 = 2.0f * - _S1208.rows[int(0)].x;
    float _S1218 = - _S1214 + _S1216;
    float _S1219 = _S1211 + - _S1215;
    float _S1220 = - _S1210 + _S1212;
    float _S1221 = _S1210 + _S1212;
    float _S1222 = _S1211 + _S1215;
    float _S1223 = _S1214 + _S1216;
    float _S1224 = quat_13.w * (_S1213 + _S1217);
    float _S1225 = quat_13.z * (_S1209 + _S1217);
    float _S1226 = quat_13.y * (_S1209 + _S1213);
    float _S1227 = quat_13.x * _S1218 + quat_13.z * _S1221 + quat_13.y * _S1222 + _S1224 + _S1224;
    float _S1228 = quat_13.x * _S1219 + quat_13.w * _S1221 + quat_13.y * _S1223 + _S1225 + _S1225;
    float _S1229 = quat_13.x * _S1220 + quat_13.w * _S1222 + quat_13.z * _S1223 + _S1226 + _S1226;
    float _S1230 = quat_13.w * _S1218 + quat_13.z * _S1219 + quat_13.y * _S1220;
    float3  _S1231 = _S1028;
    *&((&_S1231)->z) = _S1207.differential_0.rows[int(2)].z;
    *&((&_S1231)->y) = _S1207.differential_0.rows[int(1)].y;
    *&((&_S1231)->x) = _S1207.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1232;
    (&_S1232)->primal_0 = scale_12;
    (&_S1232)->differential_0 = _S1028;
    s_bwd_prop_exp_1(&_S1232, _S1231);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1233 = _S1232;
    float3  _S1234 = _S1028;
    *&((&_S1234)->z) = _S1195;
    *&((&_S1234)->y) = _S1196;
    *&((&_S1234)->x) = _S1197;
    float3  _S1235 = _S1130 + _S1234;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1236;
    (&_S1236)->primal_0 = R_14;
    (&_S1236)->differential_0 = _S1090;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1237;
    (&_S1237)->primal_0 = mean_10;
    (&_S1237)->differential_0 = _S1028;
    s_bwd_prop_mul_1(&_S1236, &_S1237, _S1235);
    float3  _S1238 = _S1235 + _S1094.differential_0;
    Matrix<float, 3, 3>  _S1239 = _S1200 + _S1201.differential_0 + _S1236.differential_0;
    float4  _S1240 = make_float4 (0.0f);
    *&((&_S1240)->w) = _S1227;
    *&((&_S1240)->z) = _S1228;
    *&((&_S1240)->y) = _S1229;
    *&((&_S1240)->x) = _S1230;
    float4  _S1241 = _S1240;
    float3  _S1242 = _S1237.differential_0 + _S1088;
    *v_mean_0 = _S1242;
    *v_quat_0 = _S1241;
    *v_scale_0 = _S1233.differential_0;
    *v_in_opacity_0 = _S1137;
    (*v_sh_coeffs_0)[int(0)] = _S1114;
    (*v_sh_coeffs_0)[int(1)] = _S1115;
    (*v_sh_coeffs_0)[int(2)] = _S1116;
    (*v_sh_coeffs_0)[int(3)] = _S1117;
    (*v_sh_coeffs_0)[int(4)] = _S1118;
    (*v_sh_coeffs_0)[int(5)] = _S1119;
    (*v_sh_coeffs_0)[int(6)] = _S1120;
    (*v_sh_coeffs_0)[int(7)] = _S1121;
    (*v_sh_coeffs_0)[int(8)] = _S1122;
    (*v_sh_coeffs_0)[int(9)] = _S1123;
    (*v_sh_coeffs_0)[int(10)] = _S1124;
    (*v_sh_coeffs_0)[int(11)] = _S1125;
    (*v_sh_coeffs_0)[int(12)] = _S1126;
    (*v_sh_coeffs_0)[int(13)] = _S1127;
    (*v_sh_coeffs_0)[int(14)] = _S1128;
    (*v_sh_coeffs_0)[int(15)] = _S1129;
    *v_R_1 = _S1239;
    *v_t_1 = _S1238;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S1243;
    DiffPair_float_0 _S1244;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S1245;
    DiffPair_float_0 _S1246;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1247;
    DiffPair_float_0 _S1248;
    DiffPair_float_0 _S1249;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1250;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S1251;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1252;
};

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S1253, float _S1254)
{
    return s_primal_ctx_atan2_0(_S1253, _S1254);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S1255;
    DiffPair_float_0 _S1256;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S1257 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S1255 = _S1257;
    _s_diff_ctx_2->_S1256 = _S1257;
    (&_s_diff_ctx_2->_S1255)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1255)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S1256)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1256)->differential_0 = 0.0f;
    DiffPair_float_0 _S1258 = *dpdpy_0;
    _s_diff_ctx_2->_S1255 = *dpdpy_0;
    DiffPair_float_0 _S1259 = *dpdpx_0;
    _s_diff_ctx_2->_S1256 = *dpdpx_0;
    float _S1260 = _S1259.primal_0 * _S1259.primal_0 + _S1258.primal_0 * _S1258.primal_0;
    float _S1261 = - _S1258.primal_0 / _S1260 * dpdOut_0;
    float _S1262 = _S1259.primal_0 / _S1260 * dpdOut_0;
    dpdpy_0->primal_0 = _S1258.primal_0;
    dpdpy_0->differential_0 = _S1262;
    dpdpx_0->primal_0 = _S1259.primal_0;
    dpdpx_0->differential_0 = _S1261;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S1263, DiffPair_float_0 * _S1264, float _S1265, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S1266 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S1243 = _S1266;
    _s_diff_ctx_3->_S1244 = _S1266;
    (&_s_diff_ctx_3->_S1243)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1243)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S1244)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1244)->differential_0 = 0.0f;
    DiffPair_float_0 _S1267 = *_S1263;
    _s_diff_ctx_3->_S1243 = *_S1263;
    DiffPair_float_0 _S1268 = *_S1264;
    _s_diff_ctx_3->_S1244 = *_S1264;
    DiffPair_float_0 _S1269 = _S1267;
    DiffPair_float_0 _S1270 = _S1268;
    s_bwd_prop_d_atan2_Intermediates_0 _S1271;
    (&_S1271)->_S1255 = _S1266;
    (&_S1271)->_S1256 = _S1266;
    s_primal_ctx_d_atan2_0(&_S1269, &_S1270, _S1265, &_S1271);
    *_S1263 = _S1269;
    *_S1264 = _S1270;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1272;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1273;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1274;
    DiffPair_float_0 _S1275;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1276;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1277;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S1278 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S1277 = _S1278;
    (&_s_diff_ctx_4->_S1277)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1277)->differential_0 = 0.0f;
    DiffPair_float_0 _S1279 = *dpdpx_1;
    _s_diff_ctx_4->_S1277 = *dpdpx_1;
    float _S1280 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S1279.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S1279.primal_0;
    dpdpx_1->differential_0 = _S1280;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1281, float _S1282, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S1283 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S1273 = _S1283;
    (&_s_diff_ctx_5->_S1273)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S1273)->differential_0 = 0.0f;
    DiffPair_float_0 _S1284 = *_S1281;
    _s_diff_ctx_5->_S1273 = *_S1281;
    DiffPair_float_0 _S1285 = _S1284;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1286;
    (&_S1286)->_S1277 = _S1283;
    s_primal_ctx_d_sqrt_0(&_S1285, _S1282, &_S1286);
    *_S1281 = _S1285;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S1287 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1288 = { _S1287, _S1287 };
    DiffPair_float_0 _S1289 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1290 = { _S1289 };
    _s_diff_ctx_6->_S1274 = _S1288;
    _s_diff_ctx_6->_S1275 = _S1289;
    _s_diff_ctx_6->_S1276 = _S1290;
    (&_s_diff_ctx_6->_S1274)->primal_0 = _S1287;
    (&_s_diff_ctx_6->_S1274)->differential_0 = _S1287;
    (&_s_diff_ctx_6->_S1275)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1275)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1291 = *dpdpx_2;
    _s_diff_ctx_6->_S1274 = *dpdpx_2;
    float _S1292 = _S1291.primal_0.x;
    float _S1293 = _S1291.primal_0.y;
    DiffPair_float_0 _S1294;
    (&_S1294)->primal_0 = _S1292 * _S1292 + _S1293 * _S1293;
    (&_S1294)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S1294, dp_s_dOut_0, &_s_diff_ctx_6->_S1276);
    _s_diff_ctx_6->_S1275 = _S1294;
    float _S1295 = _S1291.primal_0.y * _S1294.differential_0;
    float _S1296 = _S1295 + _S1295;
    float _S1297 = _S1291.primal_0.x * _S1294.differential_0;
    float _S1298 = _S1297 + _S1297;
    float2  _S1299 = _S1287;
    *&((&_S1299)->y) = _S1296;
    *&((&_S1299)->x) = _S1298;
    dpdpx_2->primal_0 = _S1291.primal_0;
    dpdpx_2->differential_0 = _S1299;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1300, float _S1301, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S1302 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1303 = { _S1302, _S1302 };
    _s_diff_ctx_7->_S1272 = _S1303;
    (&_s_diff_ctx_7->_S1272)->primal_0 = _S1302;
    (&_s_diff_ctx_7->_S1272)->differential_0 = _S1302;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1304 = *_S1300;
    _s_diff_ctx_7->_S1272 = *_S1300;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1305 = _S1304;
    DiffPair_float_0 _S1306 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1307 = { _S1306 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1308;
    (&_S1308)->_S1274 = _S1303;
    (&_S1308)->_S1275 = _S1306;
    (&_S1308)->_S1276 = _S1307;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1305, _S1301, &_S1308);
    *_S1300 = _S1305;
    return;
}

inline __device__ bool s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float4  dpintrins_0, FixedArray<float, 10>  * dpdist_coeffs_0, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S1309 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1310 = { _S1309, _S1309 };
    _s_diff_ctx_8->_S1245 = _S1309;
    _s_diff_ctx_8->_S1246 = _S1309;
    _s_diff_ctx_8->_S1247 = _S1310;
    _s_diff_ctx_8->_S1248 = _S1309;
    _s_diff_ctx_8->_S1249 = _S1309;
    _s_diff_ctx_8->_S1250 = _S1310;
    (&_s_diff_ctx_8->_S1245)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1245)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1246)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1246)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1248)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1248)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1249)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1249)->differential_0 = 0.0f;
    float2  _S1311 = make_float2 (0.0f);
    float2  _S1312 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S1313 = length_0(_S1312);
    float _S1314 = dpmean3d_0.z;
    float _S1315 = s_primal_ctx_atan2_0(_S1313, _S1314);
    float k_4;
    if(_S1315 < 0.00100000004749745f)
    {
        k_4 = (1.0f - _S1315 * _S1315 / 3.0f) / _S1314;
    }
    else
    {
        k_4 = _S1315 / _S1313;
    }
    float2  _S1316 = _S1312 * make_float2 (k_4);
    float u_38 = _S1316.x;
    float v_38 = _S1316.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float _S1317 = 2.0f * (*dpdist_coeffs_0)[int(4)];
    float _S1318 = 2.0f * (*dpdist_coeffs_0)[int(5)];
    float2  _S1319 = _S1316 * make_float2 (1.0f + r2_38 * ((*dpdist_coeffs_0)[int(0)] + r2_38 * ((*dpdist_coeffs_0)[int(1)] + r2_38 * ((*dpdist_coeffs_0)[int(2)] + r2_38 * (*dpdist_coeffs_0)[int(3)])))) + make_float2 (_S1317 * u_38 * v_38 + (*dpdist_coeffs_0)[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + (*dpdist_coeffs_0)[int(6)] * r2_38, _S1318 * u_38 * v_38 + (*dpdist_coeffs_0)[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + (*dpdist_coeffs_0)[int(7)] * r2_38);
    float2  _S1320 = _S1319 + make_float2 ((*dpdist_coeffs_0)[int(8)] * _S1319.x + (*dpdist_coeffs_0)[int(9)] * _S1319.y, 0.0f);
    float fx_18 = dpintrins_0.x;
    float fy_18 = dpintrins_0.y;
    float2  _S1321 = make_float2 (fx_18 * _S1320.x + dpintrins_0.z, fy_18 * _S1320.y + dpintrins_0.w);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (0.0f);
    float _S1322 = s_primal_ctx_s_primal_ctx_atan2_0(_S1313, _S1314);
    bool _S1323 = _S1322 < 0.00100000004749745f;
    float _S1324;
    float _S1325;
    float _S1326;
    if(_S1323)
    {
        float _S1327 = 1.0f - _S1322 * _S1322 / 3.0f;
        float _S1328 = _S1314 * _S1314;
        k_4 = _S1327 / _S1314;
        _S1324 = 0.0f;
        _S1325 = _S1328;
        _S1326 = _S1327;
    }
    else
    {
        float _S1329 = _S1313 * _S1313;
        k_4 = _S1322 / _S1313;
        _S1324 = _S1329;
        _S1325 = 0.0f;
        _S1326 = 0.0f;
    }
    float2  _S1330 = make_float2 (k_4);
    float2  _S1331 = _S1312 * make_float2 (k_4);
    float u_39 = _S1331.x;
    float v_39 = _S1331.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float _S1332 = (*dpdist_coeffs_0)[int(2)] + r2_39 * (*dpdist_coeffs_0)[int(3)];
    float _S1333 = (*dpdist_coeffs_0)[int(1)] + r2_39 * _S1332;
    float _S1334 = (*dpdist_coeffs_0)[int(0)] + r2_39 * _S1333;
    float2  _S1335 = make_float2 (fx_18, 0.0f) + make_float2 ((*dpdist_coeffs_0)[int(8)] * fx_18, (*dpdist_coeffs_0)[int(9)] * fx_18);
    float2  _S1336 = _S1331 * _S1335;
    float _S1337 = (*dpdist_coeffs_0)[int(4)] * _S1335.y;
    float _S1338 = (*dpdist_coeffs_0)[int(5)] * _S1335.x;
    float _S1339 = _S1336.x + _S1336.y;
    float _S1340 = r2_39 * _S1339;
    float _S1341 = r2_39 * _S1340;
    float _S1342 = (*dpdist_coeffs_0)[int(7)] * _S1335.y + _S1337 + (*dpdist_coeffs_0)[int(6)] * _S1335.x + _S1338 + _S1334 * _S1339 + _S1333 * _S1340 + _S1332 * _S1341 + (*dpdist_coeffs_0)[int(3)] * (r2_39 * _S1341);
    float _S1343 = v_39 * _S1342;
    float _S1344 = u_39 * _S1342;
    float2  _S1345 = make_float2 (1.0f + r2_39 * _S1334) * _S1335 + make_float2 (_S1318 * (v_39 * _S1335.y) + 2.0f * u_39 * _S1338 + 2.0f * (u_39 * _S1338) + _S1317 * (v_39 * _S1335.x) + _S1344 + _S1344, 2.0f * v_39 * _S1337 + 2.0f * (v_39 * _S1337) + _S1318 * u_39 * _S1335.y + _S1317 * u_39 * _S1335.x + _S1343 + _S1343);
    float2  _S1346 = _S1312 * _S1345;
    float2  _S1347 = _S1330 * _S1345;
    float _S1348 = _S1346.x + _S1346.y;
    if(_S1323)
    {
        float _S1349 = _S1348 / _S1325;
        float _S1350 = _S1326 * - _S1349;
        float _S1351 = _S1322 * (0.3333333432674408f * - (_S1314 * _S1349));
        k_4 = _S1351 + _S1351;
        _S1324 = _S1350;
        _S1325 = 0.0f;
    }
    else
    {
        float _S1352 = _S1348 / _S1324;
        float _S1353 = _S1322 * - _S1352;
        k_4 = _S1313 * _S1352;
        _S1324 = 0.0f;
        _S1325 = _S1353;
    }
    DiffPair_float_0 _S1354;
    (&_S1354)->primal_0 = _S1313;
    (&_S1354)->differential_0 = 0.0f;
    DiffPair_float_0 _S1355;
    (&_S1355)->primal_0 = _S1314;
    (&_S1355)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1354, &_S1355, k_4, &_s_diff_ctx_8->_S1247);
    _s_diff_ctx_8->_S1245 = _S1354;
    _s_diff_ctx_8->_S1246 = _S1355;
    float _S1356 = _S1355.differential_0 + _S1324;
    float _S1357 = _S1354.differential_0 + _S1325;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1358;
    (&_S1358)->primal_0 = _S1312;
    (&_S1358)->differential_0 = _S1311;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1359 = { _S1311, _S1311 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1360;
    (&_S1360)->_S1272 = _S1359;
    s_primal_ctx_s_bwd_length_impl_0(&_S1358, _S1357, &_S1360);
    float2  _S1361 = _S1358.differential_0 + _S1347;
    float3  _S1362 = make_float3 (_S1361.x, _S1361.y, _S1356);
    Matrix<float, 2, 3>  _S1363 = J_11;
    _S1363[int(0)] = _S1362;
    if(_S1323)
    {
        float _S1364 = 1.0f - _S1322 * _S1322 / 3.0f;
        float _S1365 = _S1314 * _S1314;
        k_4 = _S1364 / _S1314;
        _S1324 = 0.0f;
        _S1325 = _S1365;
        _S1326 = _S1364;
    }
    else
    {
        float _S1366 = _S1313 * _S1313;
        k_4 = _S1322 / _S1313;
        _S1324 = _S1366;
        _S1325 = 0.0f;
        _S1326 = 0.0f;
    }
    float2  _S1367 = make_float2 (k_4);
    float2  _S1368 = _S1312 * make_float2 (k_4);
    float u_40 = _S1368.x;
    float v_40 = _S1368.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S1369 = (*dpdist_coeffs_0)[int(2)] + r2_40 * (*dpdist_coeffs_0)[int(3)];
    float _S1370 = (*dpdist_coeffs_0)[int(1)] + r2_40 * _S1369;
    float _S1371 = (*dpdist_coeffs_0)[int(0)] + r2_40 * _S1370;
    float2  _S1372 = make_float2 (0.0f, fy_18);
    float2  _S1373 = _S1368 * _S1372;
    float _S1374 = (*dpdist_coeffs_0)[int(4)] * fy_18;
    float _S1375 = _S1373.x + _S1373.y;
    float _S1376 = r2_40 * _S1375;
    float _S1377 = r2_40 * _S1376;
    float _S1378 = (*dpdist_coeffs_0)[int(7)] * fy_18 + _S1374 + _S1371 * _S1375 + _S1370 * _S1376 + _S1369 * _S1377 + (*dpdist_coeffs_0)[int(3)] * (r2_40 * _S1377);
    float _S1379 = v_40 * _S1378;
    float _S1380 = u_40 * _S1378;
    float2  _S1381 = make_float2 (1.0f + r2_40 * _S1371) * _S1372 + make_float2 (_S1318 * (v_40 * fy_18) + _S1380 + _S1380, 2.0f * v_40 * _S1374 + 2.0f * (v_40 * _S1374) + _S1318 * u_40 * fy_18 + _S1379 + _S1379);
    float2  _S1382 = _S1312 * _S1381;
    float2  _S1383 = _S1367 * _S1381;
    float _S1384 = _S1382.x + _S1382.y;
    if(_S1323)
    {
        float _S1385 = _S1384 / _S1325;
        float _S1386 = _S1326 * - _S1385;
        float _S1387 = _S1322 * (0.3333333432674408f * - (_S1314 * _S1385));
        k_4 = _S1387 + _S1387;
        _S1324 = _S1386;
        _S1325 = 0.0f;
    }
    else
    {
        float _S1388 = _S1384 / _S1324;
        float _S1389 = _S1322 * - _S1388;
        k_4 = _S1313 * _S1388;
        _S1324 = 0.0f;
        _S1325 = _S1389;
    }
    DiffPair_float_0 _S1390;
    (&_S1390)->primal_0 = _S1313;
    (&_S1390)->differential_0 = 0.0f;
    DiffPair_float_0 _S1391;
    (&_S1391)->primal_0 = _S1314;
    (&_S1391)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1390, &_S1391, k_4, &_s_diff_ctx_8->_S1250);
    _s_diff_ctx_8->_S1248 = _S1390;
    _s_diff_ctx_8->_S1249 = _S1391;
    float _S1392 = _S1391.differential_0 + _S1324;
    float _S1393 = _S1390.differential_0 + _S1325;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1394;
    (&_S1394)->primal_0 = _S1312;
    (&_S1394)->differential_0 = _S1311;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1395;
    (&_S1395)->_S1272 = _S1359;
    s_primal_ctx_s_bwd_length_impl_0(&_S1394, _S1393, &_S1395);
    float2  _S1396 = _S1394.differential_0 + _S1383;
    float3  _S1397 = make_float3 (_S1396.x, _S1396.y, _S1392);
    _S1363[int(1)] = _S1397;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S1363, dpcov3d_0), transpose_1(_S1363));
    *dpmean2d_0 = _S1321;
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
    DiffPair_1 _S1398 = *dpdpx_3;
    float _S1399 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S1277)->primal_0);
    float _S1400 = s_primal_ctx_sqrt_0(_S1399);
    float _S1401 = 0.5f / _S1400 * (*dpdpx_3).differential_0.differential_0;
    float _S1402 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S1400 * _S1400));
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = _S1399;
    (&_S1403)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1403, _S1402);
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = 1.00000001168609742e-07f;
    (&_S1404)->differential_0 = 0.0f;
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = (&_s_diff_ctx_9->_S1277)->primal_0;
    (&_S1405)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1404, &_S1405, _S1403.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S1405.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S1401;
    dpdpx_3->primal_0 = _S1398.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S1406, DiffPair_float_0 * _S1407, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S1408 = *_S1406;
    DiffPair_float_0 _S1409 = _s_diff_ctx_10->_S1273;
    DiffPair_float_0 _S1410 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S1411;
    (&_S1411)->_S1277 = _S1410;
    s_primal_ctx_d_sqrt_0(&_S1409, (*_S1407).primal_0, &_S1411);
    DiffPair_float_0 _S1412 = { (*_S1406).differential_0.primal_0, (*_S1406).differential_0.differential_0 };
    DiffPair_1 _S1413;
    (&_S1413)->primal_0 = _s_diff_ctx_10->_S1273;
    (&_S1413)->differential_0 = _S1412;
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = (*_S1407).primal_0;
    (&_S1414)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1415 = _S1411;
    s_bwd_prop_d_sqrt_0(&_S1413, &_S1414, &_S1415);
    DiffPair_float_0 _S1416 = { _S1413.differential_0.primal_0, _S1413.differential_0.differential_0 };
    _S1407->primal_0 = (*_S1407).primal_0;
    _S1407->differential_0 = _S1414.differential_0;
    _S1406->primal_0 = _S1408.primal_0;
    _S1406->differential_0 = _S1416;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S1417, float _s_dOut_3)
{
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = (*_S1417).primal_0;
    (&_S1418)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1418, _s_dOut_3);
    _S1417->primal_0 = (*_S1417).primal_0;
    _S1417->differential_0 = _S1418.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S1419 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->y);
    DiffPair_float_0 _S1420 = { len_0, 0.0f };
    float2  _S1421 = make_float2 (0.0f);
    float _S1422 = (*dpdpx_5).differential_0.differential_0.x;
    float _S1423 = _S1422 + _S1422;
    float _S1424 = (&_s_diff_ctx_11->_S1275)->differential_0 * _S1423;
    float _S1425 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S1426 = (&_s_diff_ctx_11->_S1275)->differential_0 * _S1425;
    DiffPair_float_0 _S1427 = { 0.0f, *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->x) * _S1423 + *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->y) * _S1425 };
    DiffPair_1 _S1428;
    (&_S1428)->primal_0 = _S1420;
    (&_S1428)->differential_0 = _S1427;
    DiffPair_float_0 _S1429;
    (&_S1429)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S1429)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S1428, &_S1429, &_s_diff_ctx_11->_S1276);
    DiffPair_float_0 _S1430;
    (&_S1430)->primal_0 = len_0;
    (&_S1430)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1430, 0.0f);
    float _S1431 = _S1428.differential_0.primal_0 + _S1430.differential_0;
    float _S1432 = *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->y) * _S1431;
    float _S1433 = _S1426 + _S1432 + _S1432;
    float _S1434 = *&((&(&_s_diff_ctx_11->_S1274)->primal_0)->x) * _S1431;
    float _S1435 = _S1424 + _S1434 + _S1434;
    float2  _S1436 = _S1421;
    *&((&_S1436)->y) = _S1433;
    *&((&_S1436)->x) = _S1435;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S1419.differential_0.primal_0 + _S1436, _S1421 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S1429.differential_0;
    dpdpx_5->primal_0 = _S1419.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S1437 = (*dpdpx_7).primal_0.x;
    float _S1438 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S1439;
    (&_S1439)->primal_0 = _S1437 * _S1437 + _S1438 * _S1438;
    (&_S1439)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1439, _s_dOut_4);
    float _S1440 = (*dpdpx_7).primal_0.y * _S1439.differential_0;
    float _S1441 = _S1440 + _S1440;
    float _S1442 = (*dpdpx_7).primal_0.x * _S1439.differential_0;
    float _S1443 = _S1442 + _S1442;
    float2  _S1444 = make_float2 (0.0f);
    *&((&_S1444)->y) = _S1441;
    *&((&_S1444)->x) = _S1443;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S1444;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S1445, DiffPair_float_0 * _S1446, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S1447 = *_S1445;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1448 = _s_diff_ctx_12->_S1272;
    float2  _S1449 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1450 = { _S1449, _S1449 };
    DiffPair_float_0 _S1451 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1452 = { _S1451 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1453;
    (&_S1453)->_S1274 = _S1450;
    (&_S1453)->_S1275 = _S1451;
    (&_S1453)->_S1276 = _S1452;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1448, (*_S1446).primal_0, &_S1453);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1454 = { (*_S1445).differential_0.primal_0, (*_S1445).differential_0.differential_0 };
    DiffPair_0 _S1455;
    (&_S1455)->primal_0 = _s_diff_ctx_12->_S1272;
    (&_S1455)->differential_0 = _S1454;
    DiffPair_float_0 _S1456;
    (&_S1456)->primal_0 = (*_S1446).primal_0;
    (&_S1456)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1457 = _S1453;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S1455, &_S1456, &_S1457);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1458;
    (&_S1458)->primal_0 = (&_s_diff_ctx_12->_S1272)->primal_0;
    (&_S1458)->differential_0 = _S1449;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S1458, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1459 = { _S1455.differential_0.primal_0 + _S1458.differential_0, _S1455.differential_0.differential_0 };
    _S1446->primal_0 = (*_S1446).primal_0;
    _S1446->differential_0 = _S1456.differential_0;
    _S1445->primal_0 = _S1447.primal_0;
    _S1445->differential_0 = _S1459;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S1460 = *dpdpy_1;
    DiffPair_1 _S1461 = *dpdpx_8;
    float _S1462 = - (&_s_diff_ctx_13->_S1255)->primal_0;
    float _S1463 = (&_s_diff_ctx_13->_S1256)->primal_0 * (&_s_diff_ctx_13->_S1256)->primal_0 + (&_s_diff_ctx_13->_S1255)->primal_0 * (&_s_diff_ctx_13->_S1255)->primal_0;
    float _S1464 = _S1463 * _S1463;
    float _S1465 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S1464;
    float _S1466 = (&_s_diff_ctx_13->_S1256)->primal_0 * - _S1465;
    float _S1467 = (&_s_diff_ctx_13->_S1255)->primal_0 * _S1466;
    float _S1468 = (&_s_diff_ctx_13->_S1256)->primal_0 * _S1466;
    float _S1469 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S1464;
    float _S1470 = _S1462 * - _S1469;
    float _S1471 = (&_s_diff_ctx_13->_S1255)->primal_0 * _S1470;
    float _S1472 = (&_s_diff_ctx_13->_S1256)->primal_0 * _S1470;
    DiffPair_float_0 dpdpx_9 = { _S1472 + _S1472 + ((*dpdpx_8).differential_0.primal_0 + (_S1468 + _S1468 + _S1463 * _S1465)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S1467 + _S1467 + (*dpdpy_1).differential_0.primal_0 + _S1471 + _S1471 + - (_S1463 * _S1469), 0.0f };
    float _S1473 = (&_s_diff_ctx_13->_S1256)->primal_0 / _S1463 * (*dpdpy_1).differential_0.differential_0 + _S1462 / _S1463 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S1473;
    dpdpy_1->primal_0 = _S1460.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S1461.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S1474, DiffPair_1 * _S1475, DiffPair_float_0 * _S1476, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S1477 = *_S1474;
    DiffPair_1 _S1478 = *_S1475;
    DiffPair_float_0 _S1479 = _s_diff_ctx_14->_S1243;
    DiffPair_float_0 _S1480 = _s_diff_ctx_14->_S1244;
    DiffPair_float_0 _S1481 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S1482;
    (&_S1482)->_S1255 = _S1481;
    (&_S1482)->_S1256 = _S1481;
    s_primal_ctx_d_atan2_0(&_S1479, &_S1480, (*_S1476).primal_0, &_S1482);
    DiffPair_float_0 _S1483 = { (*_S1475).differential_0.primal_0, (*_S1475).differential_0.differential_0 };
    DiffPair_float_0 _S1484 = { (*_S1474).differential_0.primal_0, (*_S1474).differential_0.differential_0 };
    DiffPair_1 _S1485;
    (&_S1485)->primal_0 = _s_diff_ctx_14->_S1243;
    (&_S1485)->differential_0 = _S1484;
    DiffPair_1 _S1486;
    (&_S1486)->primal_0 = _s_diff_ctx_14->_S1244;
    (&_S1486)->differential_0 = _S1483;
    DiffPair_float_0 _S1487;
    (&_S1487)->primal_0 = (*_S1476).primal_0;
    (&_S1487)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S1488 = _S1482;
    s_bwd_prop_d_atan2_0(&_S1485, &_S1486, &_S1487, &_S1488);
    DiffPair_float_0 _S1489 = { _S1486.differential_0.primal_0, _S1486.differential_0.differential_0 };
    DiffPair_float_0 _S1490 = { _S1485.differential_0.primal_0, _S1485.differential_0.differential_0 };
    _S1476->primal_0 = (*_S1476).primal_0;
    _S1476->differential_0 = _S1487.differential_0;
    _S1474->primal_0 = _S1477.primal_0;
    _S1474->differential_0 = _S1490;
    _S1475->primal_0 = _S1478.primal_0;
    _S1475->differential_0 = _S1489;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S1491, DiffPair_float_0 * _S1492, float _s_dOut_5)
{
    DiffPair_float_0 _S1493;
    (&_S1493)->primal_0 = (*_S1491).primal_0;
    (&_S1493)->differential_0 = 0.0f;
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = (*_S1492).primal_0;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1493, &_S1494, _s_dOut_5);
    _S1492->primal_0 = (*_S1492).primal_0;
    _S1492->differential_0 = _S1494.differential_0;
    _S1491->primal_0 = (*_S1491).primal_0;
    _S1491->differential_0 = _S1493.differential_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpintrins_1, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_1, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1495 = *dpcov3d_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1496 = *dpintrins_1;
    FixedArray<float, 10>  _S1497 = dpdist_coeffs_1->primal_0;
    float2  _S1498 = make_float2 (0.0f);
    float2  _S1499 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S1500 = length_0(_S1499);
    float _S1501 = (*dpmean3d_1).primal_0.z;
    float _S1502 = s_primal_ctx_atan2_0(_S1500, _S1501);
    bool _S1503 = _S1502 < 0.00100000004749745f;
    float k_5;
    float _S1504;
    float _S1505;
    float _S1506;
    if(_S1503)
    {
        float _S1507 = 1.0f - _S1502 * _S1502 / 3.0f;
        float _S1508 = _S1501 * _S1501;
        k_5 = _S1507 / _S1501;
        _S1504 = 0.0f;
        _S1505 = _S1508;
        _S1506 = _S1507;
    }
    else
    {
        float _S1509 = _S1500 * _S1500;
        k_5 = _S1502 / _S1500;
        _S1504 = _S1509;
        _S1505 = 0.0f;
        _S1506 = 0.0f;
    }
    float2  _S1510 = make_float2 (k_5);
    float2  _S1511 = _S1499 * make_float2 (k_5);
    float u_41 = _S1511.x;
    float v_41 = _S1511.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float _S1512 = _S1497[int(2)] + r2_41 * _S1497[int(3)];
    float _S1513 = _S1497[int(1)] + r2_41 * _S1512;
    float _S1514 = _S1497[int(0)] + r2_41 * _S1513;
    float radial_1 = 1.0f + r2_41 * _S1514;
    float2  _S1515 = make_float2 (radial_1);
    float _S1516 = 2.0f * _S1497[int(4)];
    float _S1517 = _S1516 * u_41;
    float _S1518 = 2.0f * u_41;
    float _S1519 = r2_41 + _S1518 * u_41;
    float _S1520 = 2.0f * _S1497[int(5)];
    float _S1521 = _S1520 * u_41;
    float _S1522 = 2.0f * v_41;
    float _S1523 = r2_41 + _S1522 * v_41;
    float2  _S1524 = _S1511 * make_float2 (radial_1) + make_float2 (_S1517 * v_41 + _S1497[int(5)] * _S1519 + _S1497[int(6)] * r2_41, _S1521 * v_41 + _S1497[int(4)] * _S1523 + _S1497[int(7)] * r2_41);
    float _S1525 = _S1524.x;
    float _S1526 = _S1524.y;
    float2  _S1527 = _S1524 + make_float2 (_S1497[int(8)] * _S1525 + _S1497[int(9)] * _S1526, 0.0f);
    float fx_19 = _S1496.primal_0.x;
    float fy_19 = _S1496.primal_0.y;
    float _S1528 = _S1527.x;
    float _S1529 = _S1527.y;
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float _S1530 = s_primal_ctx_s_primal_ctx_atan2_0(_S1500, _S1501);
    bool _S1531 = _S1530 < 0.00100000004749745f;
    float _S1532;
    float _S1533;
    float _S1534;
    if(_S1531)
    {
        float _S1535 = 1.0f - _S1530 * _S1530 / 3.0f;
        float _S1536 = _S1501 * _S1501;
        k_5 = _S1535 / _S1501;
        _S1532 = 0.0f;
        _S1533 = _S1536;
        _S1534 = _S1535;
    }
    else
    {
        float _S1537 = _S1500 * _S1500;
        k_5 = _S1530 / _S1500;
        _S1532 = _S1537;
        _S1533 = 0.0f;
        _S1534 = 0.0f;
    }
    float2  _S1538 = make_float2 (k_5);
    float2  _S1539 = _S1499 * make_float2 (k_5);
    float u_42 = _S1539.x;
    float v_42 = _S1539.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float _S1540 = _S1497[int(2)] + r2_42 * _S1497[int(3)];
    float _S1541 = _S1497[int(1)] + r2_42 * _S1540;
    float _S1542 = _S1497[int(0)] + r2_42 * _S1541;
    float2  _S1543 = make_float2 (1.0f + r2_42 * _S1542);
    float _S1544 = _S1516 * u_42;
    float _S1545 = 2.0f * u_42;
    float _S1546 = _S1520 * u_42;
    float _S1547 = 2.0f * v_42;
    float2  _S1548 = make_float2 (fx_19, 0.0f) + make_float2 (_S1497[int(8)] * fx_19, _S1497[int(9)] * fx_19);
    float2  _S1549 = _S1539 * _S1548;
    float _S1550 = _S1497[int(4)] * _S1548.y;
    float _S1551 = v_42 * _S1548.y;
    float _S1552 = _S1497[int(5)] * _S1548.x;
    float _S1553 = v_42 * _S1548.x;
    float _S1554 = _S1549.x + _S1549.y;
    float _S1555 = r2_42 * _S1554;
    float _S1556 = r2_42 * _S1555;
    float _S1557 = r2_42 * _S1556;
    float _S1558 = _S1497[int(7)] * _S1548.y + _S1550 + _S1497[int(6)] * _S1548.x + _S1552 + _S1542 * _S1554 + _S1541 * _S1555 + _S1540 * _S1556 + _S1497[int(3)] * _S1557;
    float _S1559 = v_42 * _S1558;
    float _S1560 = u_42 * _S1558;
    float2  _S1561 = _S1543 * _S1548 + make_float2 (_S1520 * _S1551 + _S1545 * _S1552 + 2.0f * (u_42 * _S1552) + _S1516 * _S1553 + _S1560 + _S1560, _S1547 * _S1550 + 2.0f * (v_42 * _S1550) + _S1546 * _S1548.y + _S1544 * _S1548.x + _S1559 + _S1559);
    float2  _S1562 = _S1499 * _S1561;
    float2  _S1563 = _S1538 * _S1561;
    float _S1564 = _S1562.x + _S1562.y;
    float k_6;
    float _S1565;
    float _S1566;
    float _S1567;
    float _S1568;
    float _S1569;
    float _S1570;
    float _S1571;
    float _S1572;
    if(_S1531)
    {
        float _S1573 = _S1564 / _S1533;
        float _S1574 = _S1533 * _S1533;
        float _S1575 = - _S1573;
        float _S1576 = _S1534 * _S1575;
        float _S1577 = 0.3333333432674408f * - (_S1501 * _S1573);
        float _S1578 = _S1530 * _S1577;
        k_5 = _S1578 + _S1578;
        k_6 = _S1576;
        _S1565 = 0.0f;
        _S1566 = 0.0f;
        _S1567 = 0.0f;
        _S1568 = 0.0f;
        _S1569 = _S1577;
        _S1570 = _S1573;
        _S1571 = _S1575;
        _S1572 = _S1574;
    }
    else
    {
        float _S1579 = _S1564 / _S1532;
        float _S1580 = _S1532 * _S1532;
        float _S1581 = - _S1579;
        float _S1582 = _S1530 * _S1581;
        k_5 = _S1500 * _S1579;
        k_6 = 0.0f;
        _S1565 = _S1582;
        _S1566 = _S1579;
        _S1567 = _S1581;
        _S1568 = _S1580;
        _S1569 = 0.0f;
        _S1570 = 0.0f;
        _S1571 = 0.0f;
        _S1572 = 0.0f;
    }
    DiffPair_float_0 _S1583 = { _S1500, 0.0f };
    DiffPair_float_0 _S1584 = { _S1501, 0.0f };
    float _S1585 = (&_s_diff_ctx_15->_S1246)->differential_0 + k_6;
    float _S1586 = (&_s_diff_ctx_15->_S1245)->differential_0 + _S1565;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1587 = { _S1499, _S1498 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1588;
    (&_S1588)->primal_0 = _S1499;
    (&_S1588)->differential_0 = _S1498;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1589 = { _S1498, _S1498 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1590;
    (&_S1590)->_S1272 = _S1589;
    s_primal_ctx_s_bwd_length_impl_0(&_S1588, _S1586, &_S1590);
    float2  _S1591 = _S1588.differential_0 + _S1563;
    float3  _S1592 = make_float3 (_S1591.x, _S1591.y, _S1585);
    Matrix<float, 2, 3>  _S1593 = J_12;
    _S1593[int(0)] = _S1592;
    float _S1594;
    float _S1595;
    if(_S1531)
    {
        float _S1596 = 1.0f - _S1530 * _S1530 / 3.0f;
        float _S1597 = _S1501 * _S1501;
        k_6 = _S1596 / _S1501;
        _S1565 = 0.0f;
        _S1594 = _S1597;
        _S1595 = _S1596;
    }
    else
    {
        float _S1598 = _S1500 * _S1500;
        k_6 = _S1530 / _S1500;
        _S1565 = _S1598;
        _S1594 = 0.0f;
        _S1595 = 0.0f;
    }
    float2  _S1599 = make_float2 (k_6);
    float2  _S1600 = _S1499 * make_float2 (k_6);
    float u_43 = _S1600.x;
    float v_43 = _S1600.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float _S1601 = _S1497[int(2)] + r2_43 * _S1497[int(3)];
    float _S1602 = _S1497[int(1)] + r2_43 * _S1601;
    float _S1603 = _S1497[int(0)] + r2_43 * _S1602;
    float2  _S1604 = make_float2 (1.0f + r2_43 * _S1603);
    float _S1605 = _S1516 * u_43;
    float _S1606 = 2.0f * u_43;
    float _S1607 = _S1520 * u_43;
    float _S1608 = 2.0f * v_43;
    float2  _S1609 = make_float2 (0.0f, fy_19);
    float2  _S1610 = _S1600 * _S1609;
    float _S1611 = _S1497[int(4)] * fy_19;
    float _S1612 = v_43 * fy_19;
    float _S1613 = _S1610.x + _S1610.y;
    float _S1614 = r2_43 * _S1613;
    float _S1615 = r2_43 * _S1614;
    float _S1616 = r2_43 * _S1615;
    float _S1617 = _S1497[int(7)] * fy_19 + _S1611 + _S1603 * _S1613 + _S1602 * _S1614 + _S1601 * _S1615 + _S1497[int(3)] * _S1616;
    float _S1618 = v_43 * _S1617;
    float _S1619 = u_43 * _S1617;
    float2  _S1620 = _S1604 * _S1609 + make_float2 (_S1520 * _S1612 + _S1619 + _S1619, _S1608 * _S1611 + 2.0f * (v_43 * _S1611) + _S1607 * fy_19 + _S1618 + _S1618);
    float2  _S1621 = _S1499 * _S1620;
    float2  _S1622 = _S1599 * _S1620;
    float _S1623 = _S1621.x + _S1621.y;
    float _S1624;
    float _S1625;
    float _S1626;
    float _S1627;
    float _S1628;
    float _S1629;
    float _S1630;
    float _S1631;
    float _S1632;
    if(_S1531)
    {
        float _S1633 = _S1623 / _S1594;
        float _S1634 = _S1594 * _S1594;
        float _S1635 = - _S1633;
        float _S1636 = _S1595 * _S1635;
        float _S1637 = 0.3333333432674408f * - (_S1501 * _S1633);
        float _S1638 = _S1530 * _S1637;
        k_6 = _S1638 + _S1638;
        _S1624 = _S1636;
        _S1625 = 0.0f;
        _S1626 = 0.0f;
        _S1627 = 0.0f;
        _S1628 = 0.0f;
        _S1629 = _S1637;
        _S1630 = _S1633;
        _S1631 = _S1635;
        _S1632 = _S1634;
    }
    else
    {
        float _S1639 = _S1623 / _S1565;
        float _S1640 = _S1565 * _S1565;
        float _S1641 = - _S1639;
        float _S1642 = _S1530 * _S1641;
        k_6 = _S1500 * _S1639;
        _S1624 = 0.0f;
        _S1625 = _S1642;
        _S1626 = _S1639;
        _S1627 = _S1641;
        _S1628 = _S1640;
        _S1629 = 0.0f;
        _S1630 = 0.0f;
        _S1631 = 0.0f;
        _S1632 = 0.0f;
    }
    float _S1643 = (&_s_diff_ctx_15->_S1249)->differential_0 + _S1624;
    float _S1644 = (&_s_diff_ctx_15->_S1248)->differential_0 + _S1625;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1645;
    (&_S1645)->primal_0 = _S1499;
    (&_S1645)->differential_0 = _S1498;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1646;
    (&_S1646)->_S1272 = _S1589;
    s_primal_ctx_s_bwd_length_impl_0(&_S1645, _S1644, &_S1646);
    float2  _S1647 = _S1645.differential_0 + _S1622;
    float3  _S1648 = make_float3 (_S1647.x, _S1647.y, _S1643);
    _S1593[int(1)] = _S1648;
    Matrix<float, 3, 2>  _S1649 = transpose_1(_S1593);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1650;
    (&_S1650)->primal_0 = s_primal_ctx_mul_3(_S1593, _S1495.primal_0);
    (&_S1650)->differential_0 = J_12;
    Matrix<float, 3, 2>  _S1651 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1652;
    (&_S1652)->primal_0 = _S1649;
    (&_S1652)->differential_0 = _S1651;
    s_bwd_prop_mul_2(&_S1650, &_S1652, dpcov2d_1);
    Matrix<float, 2, 3>  _S1653 = transpose_2(_S1652.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1654;
    (&_S1654)->primal_0 = _S1593;
    (&_S1654)->differential_0 = J_12;
    Matrix<float, 3, 3>  _S1655 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1656;
    (&_S1656)->primal_0 = _S1495.primal_0;
    (&_S1656)->differential_0 = _S1655;
    s_bwd_prop_mul_3(&_S1654, &_S1656, _S1650.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1657 = _S1656;
    Matrix<float, 2, 3>  _S1658 = _S1653 + _S1654.differential_0;
    float2  _S1659 = _S1498;
    *&((&_S1659)->y) = _S1658.rows[int(1)].y;
    *&((&_S1659)->x) = _S1658.rows[int(1)].x;
    float2  _S1660 = _S1659;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1661 = { _S1498, _S1659 };
    DiffPair_0 _S1662;
    (&_S1662)->primal_0 = _S1587;
    (&_S1662)->differential_0 = _S1661;
    DiffPair_float_0 _S1663;
    (&_S1663)->primal_0 = _S1644;
    (&_S1663)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1664 = _S1646;
    s_bwd_prop_s_bwd_length_impl_0(&_S1662, &_S1663, &_S1664);
    DiffPair_0 _S1665 = _S1662;
    DiffPair_float_0 _S1666 = _S1663;
    DiffPair_float_0 _S1667 = { 0.0f, _S1658.rows[int(1)].z };
    DiffPair_float_0 _S1668 = { 0.0f, _S1663.differential_0 };
    DiffPair_1 _S1669;
    (&_S1669)->primal_0 = _S1583;
    (&_S1669)->differential_0 = _S1668;
    DiffPair_1 _S1670;
    (&_S1670)->primal_0 = _S1584;
    (&_S1670)->differential_0 = _S1667;
    DiffPair_float_0 _S1671;
    (&_S1671)->primal_0 = k_6;
    (&_S1671)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1669, &_S1670, &_S1671, &_s_diff_ctx_15->_S1250);
    DiffPair_1 _S1672 = _S1669;
    DiffPair_1 _S1673 = _S1670;
    DiffPair_float_0 _S1674 = _S1671;
    if(_S1531)
    {
        float _S1675 = _S1674.differential_0 + _S1674.differential_0;
        float _S1676 = _S1629 * _S1675;
        float _S1677 = - (0.3333333432674408f * (_S1530 * _S1675));
        float _S1678 = _S1631 * _S1658.rows[int(1)].z;
        float _S1679 = (_S1501 * _S1677 + - (_S1595 * _S1658.rows[int(1)].z)) / _S1632;
        float _S1680 = _S1623 * - _S1679;
        float _S1681 = _S1630 * _S1677 + _S1673.differential_0.primal_0;
        k_6 = _S1594 * _S1679;
        _S1624 = 0.0f;
        _S1625 = _S1680;
        _S1626 = _S1678;
        _S1627 = _S1672.differential_0.primal_0;
        _S1628 = _S1676;
        _S1629 = _S1681;
    }
    else
    {
        float _S1682 = _S1627 * _S1666.differential_0;
        float _S1683 = (_S1500 * _S1674.differential_0 + - (_S1530 * _S1666.differential_0)) / _S1628;
        float _S1684 = _S1623 * - _S1683;
        float _S1685 = _S1626 * _S1674.differential_0 + _S1672.differential_0.primal_0;
        k_6 = _S1565 * _S1683;
        _S1624 = _S1684;
        _S1625 = 0.0f;
        _S1626 = 0.0f;
        _S1627 = _S1685;
        _S1628 = _S1682;
        _S1629 = _S1673.differential_0.primal_0;
    }
    float2  _S1686 = _S1599 * _S1660;
    float2  _S1687 = _S1620 * _S1660;
    float2  _S1688 = _S1498;
    *&((&_S1688)->y) = k_6;
    *&((&_S1688)->x) = k_6;
    float2  _S1689 = _S1620 * _S1688;
    float2  _S1690 = _S1686 + _S1499 * _S1688;
    float _S1691 = _S1690.x;
    float _S1692 = _S1691 + _S1691;
    float _S1693 = _S1617 * _S1692;
    float _S1694 = _S1690.y + _S1690.y;
    float _S1695 = _S1617 * _S1694;
    float _S1696 = u_43 * _S1692 + v_43 * _S1694;
    float _S1697 = _S1497[int(3)] * _S1696;
    float _S1698 = _S1616 * _S1696;
    float _S1699 = _S1615 * _S1696;
    float _S1700 = _S1615 * _S1697;
    float _S1701 = _S1614 * _S1696;
    float _S1702 = _S1601 * _S1696 + r2_43 * _S1697;
    float _S1703 = _S1614 * _S1702;
    float _S1704 = _S1613 * _S1696;
    float _S1705 = _S1602 * _S1696 + r2_43 * _S1702;
    float _S1706 = _S1613 * _S1705;
    float _S1707 = _S1603 * _S1696 + r2_43 * _S1705;
    float _S1708 = v_43 * (_S1516 * _S1690.x);
    float _S1709 = _S1605 * _S1690.y;
    float _S1710 = _S1497[int(5)] * (_S1696 + u_43 * (2.0f * _S1690.x) + _S1606 * _S1690.x);
    float _S1711 = _S1497[int(6)] * _S1696;
    float _S1712 = _S1520 * _S1690.x;
    float _S1713 = _S1612 * _S1690.x;
    float _S1714 = v_43 * _S1712;
    float _S1715 = fy_19 * _S1712;
    float _S1716 = _S1607 * _S1690.y;
    float _S1717 = fy_19 * _S1690.y;
    float _S1718 = 2.0f * _S1690.y;
    float _S1719 = _S1611 * _S1718;
    float _S1720 = _S1611 * _S1690.y;
    float _S1721 = _S1696 + v_43 * _S1718 + _S1608 * _S1690.y;
    float _S1722 = _S1497[int(4)] * _S1721;
    float _S1723 = fy_19 * _S1721;
    float _S1724 = _S1497[int(7)] * _S1696;
    float _S1725 = fy_19 * _S1696;
    float2  _S1726 = _S1604 * _S1690;
    float2  _S1727 = _S1609 * _S1690;
    float2  _S1728 = _S1498;
    *&((&_S1728)->y) = _S1707;
    *&((&_S1728)->x) = _S1707;
    float2  _S1729 = _S1609 * _S1728;
    float _S1730 = _S1714 + _S1716 + _S1722 + _S1724;
    float _S1731 = _S1708 + _S1709 + _S1710 + _S1711;
    float2  _S1732 = _S1726 + _S1600 * _S1728;
    float2  _S1733 = _S1498;
    *&((&_S1733)->y) = _S1730;
    *&((&_S1733)->x) = _S1731;
    float2  _S1734 = _S1732 + _S1733;
    float _S1735 = _S1727.x + _S1727.y;
    float _S1736 = _S1704 + r2_43 * _S1735;
    float _S1737 = _S1701 + r2_43 * _S1736;
    float _S1738 = _S1699 + r2_43 * _S1737;
    float _S1739 = _S1700 + _S1703 + _S1706 + _S1603 * _S1735 + _S1602 * _S1736 + _S1601 * _S1737 + _S1497[int(3)] * _S1738;
    float _S1740 = v_43 * _S1739;
    float _S1741 = u_43 * _S1739;
    float2  _S1742 = _S1689 + _S1665.differential_0.primal_0;
    float2  _S1743 = _S1729 + make_float2 (_S1693 + _S1520 * _S1717 + _S1741 + _S1741, _S1695 + _S1715 + _S1719 + 2.0f * _S1720 + _S1740 + _S1740);
    float _S1744 = _S1713 + u_43 * _S1717;
    float _S1745 = _S1698 + r2_43 * _S1738;
    float2  _S1746 = _S1499 * _S1743;
    float _S1747 = _S1687.x + _S1687.y + _S1746.x + _S1746.y;
    float2  _S1748 = _S1599 * _S1743 + _S1742;
    if(_S1531)
    {
        float _S1749 = _S1501 * _S1625;
        float _S1750 = _S1747 / _S1594;
        float _S1751 = _S1530 * (0.3333333432674408f * - (_S1626 + _S1501 * _S1750));
        float _S1752 = _S1749 + _S1749 + _S1595 * - _S1750 + _S1629;
        k_6 = _S1751 + _S1751 + _S1628;
        _S1565 = _S1752;
        _S1594 = _S1627;
    }
    else
    {
        float _S1753 = _S1500 * _S1624;
        float _S1754 = _S1747 / _S1565;
        float _S1755 = _S1753 + _S1753 + _S1530 * - _S1754 + _S1627;
        k_6 = _S1500 * _S1754 + _S1628;
        _S1565 = _S1629;
        _S1594 = _S1755;
    }
    DiffPair_float_0 _S1756;
    (&_S1756)->primal_0 = _S1500;
    (&_S1756)->differential_0 = 0.0f;
    DiffPair_float_0 _S1757;
    (&_S1757)->primal_0 = _S1501;
    (&_S1757)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1756, &_S1757, k_6);
    float _S1758 = _S1757.differential_0 + _S1565;
    float _S1759 = _S1756.differential_0 + _S1594;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1760;
    (&_S1760)->primal_0 = _S1499;
    (&_S1760)->differential_0 = _S1498;
    s_bwd_length_impl_1(&_S1760, _S1759);
    float2  _S1761 = _S1760.differential_0 + _S1748;
    float3  _S1762 = make_float3 (_S1761.x, _S1761.y, _S1758);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1763;
    (&_S1763)->primal_0 = _S1499;
    (&_S1763)->differential_0 = _S1498;
    s_bwd_length_impl_1(&_S1763, 0.0f);
    float3  _S1764 = _S1762 + make_float3 (_S1763.differential_0.x, _S1763.differential_0.y, 0.0f);
    float2  _S1765 = _S1498;
    *&((&_S1765)->y) = _S1658.rows[int(0)].y;
    *&((&_S1765)->x) = _S1658.rows[int(0)].x;
    float2  _S1766 = _S1765;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1767 = { _S1498, _S1765 };
    DiffPair_0 _S1768;
    (&_S1768)->primal_0 = _S1587;
    (&_S1768)->differential_0 = _S1767;
    DiffPair_float_0 _S1769;
    (&_S1769)->primal_0 = _S1586;
    (&_S1769)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1770 = _S1590;
    s_bwd_prop_s_bwd_length_impl_0(&_S1768, &_S1769, &_S1770);
    DiffPair_0 _S1771 = _S1768;
    DiffPair_float_0 _S1772 = _S1769;
    DiffPair_float_0 _S1773 = { 0.0f, _S1658.rows[int(0)].z };
    DiffPair_float_0 _S1774 = { 0.0f, _S1769.differential_0 };
    DiffPair_1 _S1775;
    (&_S1775)->primal_0 = _S1583;
    (&_S1775)->differential_0 = _S1774;
    DiffPair_1 _S1776;
    (&_S1776)->primal_0 = _S1584;
    (&_S1776)->differential_0 = _S1773;
    DiffPair_float_0 _S1777;
    (&_S1777)->primal_0 = k_5;
    (&_S1777)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1775, &_S1776, &_S1777, &_s_diff_ctx_15->_S1247);
    DiffPair_1 _S1778 = _S1775;
    DiffPair_1 _S1779 = _S1776;
    DiffPair_float_0 _S1780 = _S1777;
    if(_S1531)
    {
        float _S1781 = _S1780.differential_0 + _S1780.differential_0;
        float _S1782 = _S1569 * _S1781;
        float _S1783 = - (0.3333333432674408f * (_S1530 * _S1781));
        float _S1784 = _S1571 * _S1658.rows[int(0)].z;
        float _S1785 = (_S1501 * _S1783 + - (_S1534 * _S1658.rows[int(0)].z)) / _S1572;
        float _S1786 = _S1564 * - _S1785;
        float _S1787 = _S1570 * _S1783 + _S1779.differential_0.primal_0;
        k_5 = _S1533 * _S1785;
        k_6 = 0.0f;
        _S1565 = _S1786;
        _S1566 = _S1784;
        _S1567 = _S1778.differential_0.primal_0;
        _S1568 = _S1782;
        _S1569 = _S1787;
    }
    else
    {
        float _S1788 = _S1567 * _S1772.differential_0;
        float _S1789 = (_S1500 * _S1780.differential_0 + - (_S1530 * _S1772.differential_0)) / _S1568;
        float _S1790 = _S1564 * - _S1789;
        float _S1791 = _S1566 * _S1780.differential_0 + _S1778.differential_0.primal_0;
        k_5 = _S1532 * _S1789;
        k_6 = _S1790;
        _S1565 = 0.0f;
        _S1566 = 0.0f;
        _S1567 = _S1791;
        _S1568 = _S1788;
        _S1569 = _S1779.differential_0.primal_0;
    }
    float2  _S1792 = _S1538 * _S1766;
    float2  _S1793 = _S1561 * _S1766;
    float2  _S1794 = _S1498;
    *&((&_S1794)->y) = k_5;
    *&((&_S1794)->x) = k_5;
    float2  _S1795 = _S1561 * _S1794;
    float2  _S1796 = _S1792 + _S1499 * _S1794;
    float _S1797 = _S1796.x;
    float _S1798 = _S1797 + _S1797;
    float _S1799 = _S1558 * _S1798;
    float _S1800 = _S1796.y + _S1796.y;
    float _S1801 = _S1558 * _S1800;
    float _S1802 = u_42 * _S1798 + v_42 * _S1800;
    float _S1803 = _S1497[int(3)] * _S1802;
    float _S1804 = _S1557 * _S1802;
    float _S1805 = _S1556 * _S1802;
    float _S1806 = _S1556 * _S1803;
    float _S1807 = _S1555 * _S1802;
    float _S1808 = _S1540 * _S1802 + r2_42 * _S1803;
    float _S1809 = _S1555 * _S1808;
    float _S1810 = _S1554 * _S1802;
    float _S1811 = _S1541 * _S1802 + r2_42 * _S1808;
    float _S1812 = _S1554 * _S1811;
    float _S1813 = _S1542 * _S1802 + r2_42 * _S1811;
    float _S1814 = _S1516 * _S1796.x;
    float _S1815 = _S1553 * _S1796.x;
    float _S1816 = v_42 * _S1814;
    float _S1817 = _S1548.x * _S1814;
    float _S1818 = _S1544 * _S1796.y;
    float _S1819 = _S1548.x * _S1796.y;
    float _S1820 = 2.0f * _S1796.x;
    float _S1821 = _S1552 * _S1820;
    float _S1822 = _S1552 * _S1796.x;
    float _S1823 = _S1802 + u_42 * _S1820 + _S1545 * _S1796.x;
    float _S1824 = _S1497[int(5)] * _S1823;
    float _S1825 = _S1548.x * _S1823;
    float _S1826 = _S1497[int(6)] * _S1802;
    float _S1827 = _S1548.x * _S1802;
    float _S1828 = _S1520 * _S1796.x;
    float _S1829 = _S1551 * _S1796.x;
    float _S1830 = v_42 * _S1828;
    float _S1831 = _S1548.y * _S1828;
    float _S1832 = _S1546 * _S1796.y;
    float _S1833 = _S1548.y * _S1796.y;
    float _S1834 = 2.0f * _S1796.y;
    float _S1835 = _S1550 * _S1834;
    float _S1836 = _S1550 * _S1796.y;
    float _S1837 = _S1802 + v_42 * _S1834 + _S1547 * _S1796.y;
    float _S1838 = _S1497[int(4)] * _S1837;
    float _S1839 = _S1548.y * _S1837;
    float _S1840 = _S1497[int(7)] * _S1802;
    float _S1841 = _S1548.y * _S1802;
    float2  _S1842 = _S1543 * _S1796;
    float2  _S1843 = _S1548 * _S1796;
    float2  _S1844 = _S1498;
    *&((&_S1844)->y) = _S1813;
    *&((&_S1844)->x) = _S1813;
    float2  _S1845 = _S1548 * _S1844;
    float _S1846 = _S1830 + _S1832 + _S1838 + _S1840;
    float _S1847 = _S1816 + _S1818 + _S1824 + _S1826;
    float2  _S1848 = _S1842 + _S1539 * _S1844;
    float2  _S1849 = _S1498;
    *&((&_S1849)->y) = _S1846;
    *&((&_S1849)->x) = _S1847;
    float2  _S1850 = _S1848 + _S1849;
    float _S1851 = fx_19 * _S1850.x;
    float _S1852 = fx_19 * _S1850.y;
    float _S1853 = _S1843.x + _S1843.y;
    float _S1854 = _S1810 + r2_42 * _S1853;
    float _S1855 = _S1807 + r2_42 * _S1854;
    float _S1856 = _S1805 + r2_42 * _S1855;
    float _S1857 = _S1806 + _S1809 + _S1812 + _S1542 * _S1853 + _S1541 * _S1854 + _S1540 * _S1855 + _S1497[int(3)] * _S1856;
    float _S1858 = v_42 * _S1857;
    float _S1859 = u_42 * _S1857;
    float2  _S1860 = _S1795 + _S1771.differential_0.primal_0;
    float _S1861 = _S1841 + _S1725;
    float _S1862 = _S1815 + u_42 * _S1819;
    float _S1863 = _S1497[int(8)] * _S1850.x + _S1497[int(9)] * _S1850.y + _S1850.x;
    float _S1864 = _S1855 + _S1737;
    float _S1865 = _S1854 + _S1736;
    float2  _S1866 = _S1845 + make_float2 (_S1799 + _S1821 + _S1520 * _S1833 + 2.0f * _S1822 + _S1516 * _S1819 + _S1859 + _S1859, _S1801 + _S1817 + _S1831 + _S1835 + 2.0f * _S1836 + _S1858 + _S1858);
    float _S1867 = _S1829 + u_42 * _S1833 + _S1744;
    float _S1868 = _S1804 + r2_42 * _S1856 + _S1745;
    float _S1869 = _S1856 + _S1738;
    float _S1870 = _S1839 + _S1723;
    float2  _S1871 = _S1499 * _S1866;
    float _S1872 = _S1793.x + _S1793.y + _S1871.x + _S1871.y;
    float2  _S1873 = _S1538 * _S1866 + _S1860;
    if(_S1531)
    {
        float _S1874 = _S1501 * _S1565;
        float _S1875 = _S1872 / _S1533;
        float _S1876 = _S1530 * (0.3333333432674408f * - (_S1566 + _S1501 * _S1875));
        float _S1877 = _S1874 + _S1874 + _S1534 * - _S1875 + _S1569;
        k_5 = _S1876 + _S1876 + _S1568;
        _S1532 = _S1877;
        _S1533 = _S1567;
    }
    else
    {
        float _S1878 = _S1500 * k_6;
        float _S1879 = _S1872 / _S1532;
        float _S1880 = _S1878 + _S1878 + _S1530 * - _S1879 + _S1567;
        k_5 = _S1500 * _S1879 + _S1568;
        _S1532 = _S1569;
        _S1533 = _S1880;
    }
    DiffPair_float_0 _S1881;
    (&_S1881)->primal_0 = _S1500;
    (&_S1881)->differential_0 = 0.0f;
    DiffPair_float_0 _S1882;
    (&_S1882)->primal_0 = _S1501;
    (&_S1882)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1881, &_S1882, k_5);
    float _S1883 = _S1882.differential_0 + _S1532;
    float _S1884 = _S1881.differential_0 + _S1533;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1885;
    (&_S1885)->primal_0 = _S1499;
    (&_S1885)->differential_0 = _S1498;
    s_bwd_length_impl_1(&_S1885, _S1884);
    float2  _S1886 = _S1885.differential_0 + _S1873;
    float3  _S1887 = make_float3 (_S1886.x, _S1886.y, _S1883);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1888;
    (&_S1888)->primal_0 = _S1499;
    (&_S1888)->differential_0 = _S1498;
    s_bwd_length_impl_1(&_S1888, 0.0f);
    float _S1889 = fx_19 * dpmean2d_1.x;
    float2  _S1890 = make_float2 (_S1889, fy_19 * dpmean2d_1.y) + make_float2 (_S1497[int(8)] * _S1889, _S1497[int(9)] * _S1889);
    float2  _S1891 = _S1511 * _S1890;
    float2  _S1892 = _S1515 * _S1890;
    float _S1893 = _S1497[int(4)] * _S1890.y;
    float _S1894 = v_41 * _S1890.y;
    float _S1895 = _S1497[int(5)] * _S1890.x;
    float _S1896 = v_41 * _S1890.x;
    float _S1897 = _S1891.x + _S1891.y;
    float _S1898 = r2_41 * _S1897;
    float _S1899 = r2_41 * _S1898;
    float _S1900 = r2_41 * _S1899;
    float _S1901 = _S1497[int(7)] * _S1890.y + _S1893 + _S1497[int(6)] * _S1890.x + _S1895 + _S1514 * _S1897 + _S1513 * _S1898 + _S1512 * _S1899 + _S1497[int(3)] * _S1900;
    float _S1902 = v_41 * _S1901;
    float _S1903 = u_41 * _S1901;
    float _S1904 = _S1522 * _S1893 + 2.0f * (v_41 * _S1893) + _S1521 * _S1890.y + _S1517 * _S1890.x + _S1902 + _S1902;
    float _S1905 = _S1520 * _S1894 + _S1518 * _S1895 + 2.0f * (u_41 * _S1895) + _S1516 * _S1896 + _S1903 + _S1903;
    float _S1906 = _S1526 * _S1889 + _S1852;
    float _S1907 = _S1525 * _S1889 + _S1851;
    float _S1908 = r2_41 * _S1890.y + _S1861;
    float _S1909 = r2_41 * _S1890.x + _S1827;
    float _S1910 = 2.0f * (u_41 * _S1894 + _S1867) + _S1519 * _S1890.x + _S1825;
    float _S1911 = _S1523 * _S1890.y + 2.0f * (u_41 * _S1896 + _S1862) + _S1870;
    float _S1912 = r2_41 * _S1900 + _S1868;
    float _S1913 = _S1900 + _S1869;
    float _S1914 = _S1899 + _S1864;
    float _S1915 = _S1898 + _S1865;
    float3  _S1916 = _S1887 + make_float3 (_S1888.differential_0.x, _S1888.differential_0.y, 0.0f) + _S1764;
    float4  _S1917 = make_float4 (_S1528 * dpmean2d_1.x + _S1863, _S1529 * dpmean2d_1.y + _S1734.y, dpmean2d_1.x, dpmean2d_1.y);
    FixedArray<float, 10>  _S1918;
    _S1918[int(0)] = 0.0f;
    _S1918[int(1)] = 0.0f;
    _S1918[int(2)] = 0.0f;
    _S1918[int(3)] = 0.0f;
    _S1918[int(4)] = 0.0f;
    _S1918[int(5)] = 0.0f;
    _S1918[int(6)] = 0.0f;
    _S1918[int(7)] = 0.0f;
    _S1918[int(8)] = 0.0f;
    _S1918[int(9)] = 0.0f;
    _S1918[int(9)] = _S1906;
    _S1918[int(8)] = _S1907;
    _S1918[int(7)] = _S1908;
    _S1918[int(6)] = _S1909;
    _S1918[int(5)] = _S1910;
    _S1918[int(4)] = _S1911;
    _S1918[int(3)] = _S1912;
    _S1918[int(2)] = _S1913;
    _S1918[int(1)] = _S1914;
    _S1918[int(0)] = _S1915;
    FixedArray<float, 10>  _S1919 = {
        _S1918[int(0)], _S1918[int(1)], _S1918[int(2)], _S1918[int(3)], _S1918[int(4)], _S1918[int(5)], _S1918[int(6)], _S1918[int(7)], _S1918[int(8)], _S1918[int(9)]
    };
    float2  _S1920 = _S1892 + make_float2 (_S1905, _S1904);
    float2  _S1921 = _S1499 * _S1920;
    float2  _S1922 = _S1510 * _S1920;
    float _S1923 = _S1921.x + _S1921.y;
    if(_S1503)
    {
        float _S1924 = _S1923 / _S1505;
        float _S1925 = _S1506 * - _S1924;
        float _S1926 = _S1502 * (0.3333333432674408f * - (_S1501 * _S1924));
        k_5 = _S1926 + _S1926;
        _S1504 = _S1925;
        _S1505 = 0.0f;
    }
    else
    {
        float _S1927 = _S1923 / _S1504;
        float _S1928 = _S1502 * - _S1927;
        k_5 = _S1500 * _S1927;
        _S1504 = 0.0f;
        _S1505 = _S1928;
    }
    DiffPair_float_0 _S1929;
    (&_S1929)->primal_0 = _S1500;
    (&_S1929)->differential_0 = 0.0f;
    DiffPair_float_0 _S1930;
    (&_S1930)->primal_0 = _S1501;
    (&_S1930)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1929, &_S1930, k_5);
    float _S1931 = _S1930.differential_0 + _S1504;
    float _S1932 = _S1929.differential_0 + _S1505;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1933;
    (&_S1933)->primal_0 = _S1499;
    (&_S1933)->differential_0 = _S1498;
    s_bwd_length_impl_1(&_S1933, _S1932);
    float2  _S1934 = _S1933.differential_0 + _S1922;
    dpdist_coeffs_1->primal_0 = dpdist_coeffs_1->primal_0;
    dpdist_coeffs_1->differential_0 = _S1919;
    dpintrins_1->primal_0 = (*dpintrins_1).primal_0;
    dpintrins_1->differential_0 = _S1917;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1657.differential_0;
    float3  _S1935 = _S1916 + make_float3 (_S1934.x, _S1934.y, _S1931);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1935;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_20, float fy_20, float cx_16, float cy_16, FixedArray<float, 10>  * dist_coeffs_25, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1936 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1937 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1938 = { _S1937, _S1937 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1939 = { _S1937, _S1937, _S1938, _S1937, _S1937, _S1938 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1940;
    (&_S1940)->_S1251 = _S1936;
    (&_S1940)->_S1252 = _S1939;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float4  intrins_11 = make_float4 (fx_20, fy_20, cx_16, cy_16);
    float3  _S1941 = s_primal_ctx_exp_0(scale_13);
    float _S1942 = quat_14.y;
    float x2_14 = _S1942 * _S1942;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1943 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1941.x, 0.0f, 0.0f, 0.0f, _S1941.y, 0.0f, 0.0f, 0.0f, _S1941.z);
    Matrix<float, 3, 3>  _S1944 = s_primal_ctx_mul_2(_S1943, S_1);
    Matrix<float, 3, 3>  _S1945 = transpose_0(_S1944);
    Matrix<float, 3, 3>  _S1946 = s_primal_ctx_mul_2(_S1944, _S1945);
    Matrix<float, 3, 3>  _S1947 = s_primal_ctx_mul_2(R_15, _S1946);
    Matrix<float, 3, 3>  _S1948 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1949 = s_primal_ctx_mul_2(_S1947, _S1948);
    Matrix<float, 2, 2>  _S1950 = _S1936;
    float2  _S1951 = make_float2 (0.0f);
    float2  _S1952 = _S1951;
    bool _S1953 = s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1949, intrins_11, dist_coeffs_25, &_S1950, &_S1952, &(&_S1940)->_S1252);
    (&_S1940)->_S1251 = _S1950;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1954 = _S1940;
    float _S1955 = _S1940._S1251.rows[int(0)].y * _S1940._S1251.rows[int(1)].x;
    float det_orig_12 = _S1940._S1251.rows[int(0)].x * _S1940._S1251.rows[int(1)].y - _S1955;
    float _S1956 = _S1940._S1251.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1957 = _S1940._S1251;
    *&(((&_S1957)->rows + (int(0)))->x) = _S1956;
    float _S1958 = _S1940._S1251.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1957)->rows + (int(1)))->y) = _S1958;
    Matrix<float, 2, 2>  _S1959 = _S1957;
    Matrix<float, 2, 2>  _S1960 = _S1957;
    float det_blur_7 = _S1956 * _S1958 - _S1955;
    float _S1961 = det_orig_12 / det_blur_7;
    float _S1962 = det_blur_7 * det_blur_7;
    float _S1963 = s_primal_ctx_max_0(0.0f, _S1961);
    float _S1964 = s_primal_ctx_sqrt_0(_S1963);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1965 = - _S1940._S1251.rows[int(0)].y;
    float _S1966 = - _S1940._S1251.rows[int(1)].x;
    float _S1967 = - in_opacity_11;
    float _S1968 = 1.0f + s_primal_ctx_exp_1(_S1967);
    float _S1969 = 1.0f / _S1968;
    float _S1970 = _S1968 * _S1968;
    float _S1971;
    if(antialiased_11)
    {
        _S1971 = _S1969 * _S1964;
    }
    else
    {
        _S1971 = _S1969;
    }
    float _S1972 = _S1971 / 0.00392156885936856f;
    float _S1973 = 2.0f * s_primal_ctx_log_0(_S1972);
    float _S1974 = s_primal_ctx_sqrt_0(_S1973);
    float _S1975 = _S1959.rows[int(0)].x;
    float _S1976 = _S1960.rows[int(1)].y;
    float _S1977 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1978 = mean_11 - - s_primal_ctx_mul_1(_S1948, t_14);
    float _S1979 = _S1978.x;
    float _S1980 = _S1978.y;
    float _S1981 = _S1978.z;
    float _S1982 = _S1979 * _S1979 + _S1980 * _S1980 + _S1981 * _S1981;
    float _S1983 = s_primal_ctx_sqrt_0(_S1982);
    float x_36 = _S1979 / _S1983;
    float3  _S1984 = make_float3 (x_36);
    float _S1985 = _S1983 * _S1983;
    float y_14 = _S1980 / _S1983;
    float z_11 = _S1981 / _S1983;
    float3  _S1986 = make_float3 (z_11);
    float _S1987 = - y_14;
    float3  _S1988 = make_float3 (_S1987);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_36 * x_36 - y_14 * y_14;
    float _S1989 = 2.0f * x_36;
    float fS1_11 = _S1989 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S1990 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_36;
    float3  _S1991 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1992 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1993 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1994 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1995 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1995;
    float3  _S1996 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_36;
    float3  _S1997 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1998 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1999 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S2000 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_36 * fC1_11 - y_14 * fS1_11);
    float3  _S2001 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_36 * fS1_11 + y_14 * fC1_11);
    float3  _S2002 = make_float3 (pSH9_1);
    float3  _S2003 = make_float3 (0.0f);
    float3  _S2004 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2005;
    (&_S2005)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1987) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_36) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S2005)->differential_0 = _S2004;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2006;
    (&_S2006)->primal_0 = _S2003;
    (&_S2006)->differential_0 = _S2004;
    s_bwd_prop_max_0(&_S2005, &_S2006, v_rgb_1);
    float3  _S2007 = _S2001 * _S2005.differential_0;
    float3  _S2008 = (*sh_coeffs_11)[int(15)] * _S2005.differential_0;
    float3  _S2009 = _S1999 * _S2005.differential_0;
    float3  _S2010 = (*sh_coeffs_11)[int(14)] * _S2005.differential_0;
    float3  _S2011 = _S1997 * _S2005.differential_0;
    float3  _S2012 = (*sh_coeffs_11)[int(13)] * _S2005.differential_0;
    float3  _S2013 = _S1996 * _S2005.differential_0;
    float3  _S2014 = (*sh_coeffs_11)[int(12)] * _S2005.differential_0;
    float3  _S2015 = _S1998 * _S2005.differential_0;
    float3  _S2016 = (*sh_coeffs_11)[int(11)] * _S2005.differential_0;
    float3  _S2017 = _S2000 * _S2005.differential_0;
    float3  _S2018 = (*sh_coeffs_11)[int(10)] * _S2005.differential_0;
    float3  _S2019 = _S2002 * _S2005.differential_0;
    float3  _S2020 = (*sh_coeffs_11)[int(9)] * _S2005.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S2020.x + _S2020.y + _S2020.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S2008.x + _S2008.y + _S2008.z);
    float _S2021 = _S2018.x + _S2018.y + _S2018.z;
    float _S2022 = _S2010.x + _S2010.y + _S2010.z;
    float _S2023 = _S2016.x + _S2016.y + _S2016.z;
    float _S2024 = _S2012.x + _S2012.y + _S2012.z;
    float _S2025 = _S2014.x + _S2014.y + _S2014.z;
    float _S2026 = - s_diff_fC2_T_1;
    float3  _S2027 = _S1993 * _S2005.differential_0;
    float3  _S2028 = (*sh_coeffs_11)[int(8)] * _S2005.differential_0;
    float3  _S2029 = _S1991 * _S2005.differential_0;
    float3  _S2030 = (*sh_coeffs_11)[int(7)] * _S2005.differential_0;
    float3  _S2031 = _S1990 * _S2005.differential_0;
    float3  _S2032 = (*sh_coeffs_11)[int(6)] * _S2005.differential_0;
    float3  _S2033 = _S1992 * _S2005.differential_0;
    float3  _S2034 = (*sh_coeffs_11)[int(5)] * _S2005.differential_0;
    float3  _S2035 = _S1994 * _S2005.differential_0;
    float3  _S2036 = (*sh_coeffs_11)[int(4)] * _S2005.differential_0;
    float _S2037 = _S2034.x + _S2034.y + _S2034.z;
    float _S2038 = _S2030.x + _S2030.y + _S2030.z;
    float _S2039 = fTmp1B_11 * _S2021 + x_36 * s_diff_fS2_T_1 + y_14 * _S2026 + 0.54627424478530884f * (_S2036.x + _S2036.y + _S2036.z);
    float _S2040 = fTmp1B_11 * _S2022 + y_14 * s_diff_fS2_T_1 + x_36 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S2028.x + _S2028.y + _S2028.z);
    float _S2041 = y_14 * - _S2040;
    float _S2042 = x_36 * _S2040;
    float _S2043 = z_11 * (1.86588168144226074f * (z_11 * _S2025) + -2.28522896766662598f * (y_14 * _S2023 + x_36 * _S2024) + 0.94617468118667603f * (_S2032.x + _S2032.y + _S2032.z));
    float3  _S2044 = make_float3 (0.48860251903533936f) * _S2005.differential_0;
    float3  _S2045 = - _S2044;
    float3  _S2046 = _S1984 * _S2045;
    float3  _S2047 = (*sh_coeffs_11)[int(3)] * _S2045;
    float3  _S2048 = _S1986 * _S2044;
    float3  _S2049 = (*sh_coeffs_11)[int(2)] * _S2044;
    float3  _S2050 = _S1988 * _S2044;
    float3  _S2051 = (*sh_coeffs_11)[int(1)] * _S2044;
    float _S2052 = (_S1995 * _S2025 + 1.44530570507049561f * (fS1_11 * _S2021 + fC1_11 * _S2022) + -1.09254848957061768f * (y_14 * _S2037 + x_36 * _S2038) + _S2043 + _S2043 + _S2049.x + _S2049.y + _S2049.z) / _S1985;
    float _S2053 = _S1983 * _S2052;
    float _S2054 = (fTmp0C_11 * _S2023 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S2026 + fTmp0B_11 * _S2037 + _S1989 * _S2039 + _S2041 + _S2041 + - (_S2051.x + _S2051.y + _S2051.z)) / _S1985;
    float _S2055 = _S1983 * _S2054;
    float _S2056 = (fTmp0C_11 * _S2024 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S2038 + 2.0f * (y_14 * _S2039) + _S2042 + _S2042 + _S2047.x + _S2047.y + _S2047.z) / _S1985;
    float _S2057 = _S1983 * _S2056;
    float _S2058 = _S1981 * - _S2052 + _S1980 * - _S2054 + _S1979 * - _S2056;
    DiffPair_float_0 _S2059;
    (&_S2059)->primal_0 = _S1982;
    (&_S2059)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2059, _S2058);
    float _S2060 = _S1981 * _S2059.differential_0;
    float _S2061 = _S1980 * _S2059.differential_0;
    float _S2062 = _S1979 * _S2059.differential_0;
    float3  _S2063 = make_float3 (0.282094806432724f) * _S2005.differential_0;
    float3  _S2064 = make_float3 (_S2057 + _S2062 + _S2062, _S2055 + _S2061 + _S2061, _S2053 + _S2060 + _S2060);
    float3  _S2065 = - - _S2064;
    Matrix<float, 3, 3>  _S2066 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2067;
    (&_S2067)->primal_0 = _S1948;
    (&_S2067)->differential_0 = _S2066;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2068;
    (&_S2068)->primal_0 = t_14;
    (&_S2068)->differential_0 = _S2004;
    s_bwd_prop_mul_1(&_S2067, &_S2068, _S2065);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2069 = _S2067;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2070 = _S2068;
    float2  _S2071 = _S1951;
    *&((&_S2071)->y) = v_conic_1.z;
    float2  _S2072 = _S1951;
    *&((&_S2072)->y) = v_conic_1.y;
    *&((&_S2072)->x) = v_conic_1.x;
    float _S2073 = 0.5f * v_depth_1;
    DiffPair_float_0 _S2074;
    (&_S2074)->primal_0 = _S1977;
    (&_S2074)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2074, _S2073);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2075;
    (&_S2075)->primal_0 = mean_c_11;
    (&_S2075)->differential_0 = _S2004;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2076;
    (&_S2076)->primal_0 = mean_c_11;
    (&_S2076)->differential_0 = _S2004;
    s_bwd_prop_dot_0(&_S2075, &_S2076, _S2074.differential_0);
    DiffPair_float_0 _S2077;
    (&_S2077)->primal_0 = _S1976;
    (&_S2077)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2077, 0.0f);
    DiffPair_float_0 _S2078;
    (&_S2078)->primal_0 = _S1975;
    (&_S2078)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2078, 0.0f);
    DiffPair_float_0 _S2079;
    (&_S2079)->primal_0 = 3.32999992370605469f;
    (&_S2079)->differential_0 = 0.0f;
    DiffPair_float_0 _S2080;
    (&_S2080)->primal_0 = _S1974;
    (&_S2080)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2079, &_S2080, 0.0f);
    DiffPair_float_0 _S2081;
    (&_S2081)->primal_0 = _S1973;
    (&_S2081)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2081, _S2080.differential_0);
    float _S2082 = 2.0f * _S2081.differential_0;
    DiffPair_float_0 _S2083;
    (&_S2083)->primal_0 = _S1972;
    (&_S2083)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2083, _S2082);
    float2  _S2084 = make_float2 (_S2078.differential_0, 0.0f);
    float _S2085 = v_opacity_1 + 254.9999847412109375f * _S2083.differential_0;
    Matrix<float, 2, 2>  _S2086 = _S1936;
    _S2086[int(1)] = _S2071;
    _S2086[int(0)] = _S2072;
    Matrix<float, 2, 2>  _S2087 = _S2086;
    FixedArray<float3 , 16>  _S2088;
    _S2088[int(0)] = _S2004;
    _S2088[int(1)] = _S2004;
    _S2088[int(2)] = _S2004;
    _S2088[int(3)] = _S2004;
    _S2088[int(4)] = _S2004;
    _S2088[int(5)] = _S2004;
    _S2088[int(6)] = _S2004;
    _S2088[int(7)] = _S2004;
    _S2088[int(8)] = _S2004;
    _S2088[int(9)] = _S2004;
    _S2088[int(10)] = _S2004;
    _S2088[int(11)] = _S2004;
    _S2088[int(12)] = _S2004;
    _S2088[int(13)] = _S2004;
    _S2088[int(14)] = _S2004;
    _S2088[int(15)] = _S2004;
    _S2088[int(7)] = _S2029;
    _S2088[int(0)] = _S2063;
    _S2088[int(1)] = _S2050;
    _S2088[int(2)] = _S2048;
    _S2088[int(3)] = _S2046;
    _S2088[int(4)] = _S2035;
    _S2088[int(5)] = _S2033;
    _S2088[int(6)] = _S2031;
    _S2088[int(15)] = _S2007;
    _S2088[int(8)] = _S2027;
    _S2088[int(9)] = _S2019;
    _S2088[int(10)] = _S2017;
    _S2088[int(11)] = _S2015;
    _S2088[int(12)] = _S2013;
    _S2088[int(13)] = _S2011;
    _S2088[int(14)] = _S2009;
    float3  _S2089 = _S2088[int(0)];
    float3  _S2090 = _S2088[int(1)];
    float3  _S2091 = _S2088[int(2)];
    float3  _S2092 = _S2088[int(3)];
    float3  _S2093 = _S2088[int(4)];
    float3  _S2094 = _S2088[int(5)];
    float3  _S2095 = _S2088[int(6)];
    float3  _S2096 = _S2088[int(7)];
    float3  _S2097 = _S2088[int(8)];
    float3  _S2098 = _S2088[int(9)];
    float3  _S2099 = _S2088[int(10)];
    float3  _S2100 = _S2088[int(11)];
    float3  _S2101 = _S2088[int(12)];
    float3  _S2102 = _S2088[int(13)];
    float3  _S2103 = _S2088[int(14)];
    float3  _S2104 = _S2088[int(15)];
    float3  _S2105 = _S2076.differential_0 + _S2075.differential_0;
    float2  _S2106 = make_float2 (0.0f, _S2077.differential_0);
    float _S2107;
    if(antialiased_11)
    {
        float _S2108 = _S1969 * _S2085;
        _S1971 = _S1964 * _S2085;
        _S2107 = _S2108;
    }
    else
    {
        _S1971 = _S2085;
        _S2107 = 0.0f;
    }
    float _S2109 = - (_S1971 / _S1970);
    DiffPair_float_0 _S2110;
    (&_S2110)->primal_0 = _S1967;
    (&_S2110)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2110, _S2109);
    float _S2111 = - _S2110.differential_0;
    float _S2112 = invdet_8 * _S2087.rows[int(1)].y;
    float _S2113 = - (invdet_8 * _S2087.rows[int(1)].x);
    float _S2114 = - (invdet_8 * _S2087.rows[int(0)].y);
    float _S2115 = invdet_8 * _S2087.rows[int(0)].x;
    float _S2116 = - ((_S1956 * _S2087.rows[int(1)].y + _S1966 * _S2087.rows[int(1)].x + _S1965 * _S2087.rows[int(0)].y + _S1958 * _S2087.rows[int(0)].x) / _S1962);
    DiffPair_float_0 _S2117;
    (&_S2117)->primal_0 = _S1963;
    (&_S2117)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2117, _S2107);
    DiffPair_float_0 _S2118;
    (&_S2118)->primal_0 = 0.0f;
    (&_S2118)->differential_0 = 0.0f;
    DiffPair_float_0 _S2119;
    (&_S2119)->primal_0 = _S1961;
    (&_S2119)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2118, &_S2119, _S2117.differential_0);
    float _S2120 = _S2119.differential_0 / _S1962;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2120;
    float _S2121 = _S2116 + det_orig_12 * - _S2120;
    float _S2122 = - _S2121;
    float _S2123 = _S1956 * _S2121;
    float _S2124 = _S1958 * _S2121;
    Matrix<float, 2, 2>  _S2125 = _S1936;
    _S2125[int(1)] = _S2106;
    _S2125[int(0)] = _S2084;
    _S1957 = _S2125;
    *&(((&_S1957)->rows + (int(1)))->y) = 0.0f;
    float _S2126 = _S2115 + _S2123 + _S2125.rows[int(1)].y;
    *&(((&_S1957)->rows + (int(0)))->x) = 0.0f;
    float _S2127 = _S2112 + _S2124 + _S2125.rows[int(0)].x;
    float _S2128 = _S2122 + - s_diff_det_orig_T_1;
    float _S2129 = _S2113 + _S1954._S1251.rows[int(0)].y * _S2128;
    float _S2130 = _S2114 + _S1954._S1251.rows[int(1)].x * _S2128;
    float _S2131 = _S1954._S1251.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2132 = _S2126 + _S1954._S1251.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2133 = _S1951;
    *&((&_S2133)->x) = _S2129;
    *&((&_S2133)->y) = _S2132;
    float _S2134 = _S2127 + _S2131;
    float2  _S2135 = _S1951;
    *&((&_S2135)->y) = _S2130;
    *&((&_S2135)->x) = _S2134;
    Matrix<float, 2, 2>  _S2136 = _S1936;
    _S2136[int(1)] = _S2133;
    _S2136[int(0)] = _S2135;
    Matrix<float, 2, 2>  _S2137 = _S1957 + _S2136;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2138;
    (&_S2138)->primal_0 = mean_c_11;
    (&_S2138)->differential_0 = _S2004;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2139;
    (&_S2139)->primal_0 = _S1949;
    (&_S2139)->differential_0 = _S2066;
    float4  _S2140 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2141;
    (&_S2141)->primal_0 = intrins_11;
    (&_S2141)->differential_0 = _S2140;
    FixedArray<float, 10>  _S2142 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2143;
    (&_S2143)->primal_0 = *dist_coeffs_25;
    (&_S2143)->differential_0 = _S2142;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2144 = _S1954._S1252;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2138, &_S2139, &_S2141, &_S2143, _S2137, v_mean2d_1, &_S2144);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2145;
    (&_S2145)->primal_0 = _S1947;
    (&_S2145)->differential_0 = _S2066;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2146;
    (&_S2146)->primal_0 = _S1948;
    (&_S2146)->differential_0 = _S2066;
    s_bwd_prop_mul_4(&_S2145, &_S2146, _S2139.differential_0);
    Matrix<float, 3, 3>  _S2147 = transpose_0(_S2146.differential_0 + _S2069.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2148;
    (&_S2148)->primal_0 = R_15;
    (&_S2148)->differential_0 = _S2066;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2149;
    (&_S2149)->primal_0 = _S1946;
    (&_S2149)->differential_0 = _S2066;
    s_bwd_prop_mul_4(&_S2148, &_S2149, _S2145.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2150;
    (&_S2150)->primal_0 = _S1944;
    (&_S2150)->differential_0 = _S2066;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2151;
    (&_S2151)->primal_0 = _S1945;
    (&_S2151)->differential_0 = _S2066;
    s_bwd_prop_mul_4(&_S2150, &_S2151, _S2149.differential_0);
    Matrix<float, 3, 3>  _S2152 = _S2150.differential_0 + transpose_0(_S2151.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2153;
    (&_S2153)->primal_0 = _S1943;
    (&_S2153)->differential_0 = _S2066;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2154;
    (&_S2154)->primal_0 = S_1;
    (&_S2154)->differential_0 = _S2066;
    s_bwd_prop_mul_4(&_S2153, &_S2154, _S2152);
    Matrix<float, 3, 3>  _S2155 = transpose_0(_S2153.differential_0);
    float _S2156 = 2.0f * - _S2155.rows[int(2)].z;
    float _S2157 = 2.0f * _S2155.rows[int(2)].y;
    float _S2158 = 2.0f * _S2155.rows[int(2)].x;
    float _S2159 = 2.0f * _S2155.rows[int(1)].z;
    float _S2160 = 2.0f * - _S2155.rows[int(1)].y;
    float _S2161 = 2.0f * _S2155.rows[int(1)].x;
    float _S2162 = 2.0f * _S2155.rows[int(0)].z;
    float _S2163 = 2.0f * _S2155.rows[int(0)].y;
    float _S2164 = 2.0f * - _S2155.rows[int(0)].x;
    float _S2165 = - _S2161 + _S2163;
    float _S2166 = _S2158 + - _S2162;
    float _S2167 = - _S2157 + _S2159;
    float _S2168 = _S2157 + _S2159;
    float _S2169 = _S2158 + _S2162;
    float _S2170 = _S2161 + _S2163;
    float _S2171 = quat_14.w * (_S2160 + _S2164);
    float _S2172 = quat_14.z * (_S2156 + _S2164);
    float _S2173 = quat_14.y * (_S2156 + _S2160);
    float _S2174 = quat_14.x * _S2165 + quat_14.z * _S2168 + quat_14.y * _S2169 + _S2171 + _S2171;
    float _S2175 = quat_14.x * _S2166 + quat_14.w * _S2168 + quat_14.y * _S2170 + _S2172 + _S2172;
    float _S2176 = quat_14.x * _S2167 + quat_14.w * _S2169 + quat_14.z * _S2170 + _S2173 + _S2173;
    float _S2177 = quat_14.w * _S2165 + quat_14.z * _S2166 + quat_14.y * _S2167;
    float3  _S2178 = _S2004;
    *&((&_S2178)->z) = _S2154.differential_0.rows[int(2)].z;
    *&((&_S2178)->y) = _S2154.differential_0.rows[int(1)].y;
    *&((&_S2178)->x) = _S2154.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2179;
    (&_S2179)->primal_0 = scale_13;
    (&_S2179)->differential_0 = _S2004;
    s_bwd_prop_exp_1(&_S2179, _S2178);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2180 = _S2179;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2181;
    (&_S2181)->primal_0 = mean_c_11;
    (&_S2181)->differential_0 = _S2004;
    s_bwd_length_impl_0(&_S2181, 0.0f);
    float3  _S2182 = _S2138.differential_0 + _S2181.differential_0 + _S2105;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2183;
    (&_S2183)->primal_0 = R_15;
    (&_S2183)->differential_0 = _S2066;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2184;
    (&_S2184)->primal_0 = mean_11;
    (&_S2184)->differential_0 = _S2004;
    s_bwd_prop_mul_1(&_S2183, &_S2184, _S2182);
    float3  _S2185 = _S2182 + _S2070.differential_0;
    Matrix<float, 3, 3>  _S2186 = _S2147 + _S2148.differential_0 + _S2183.differential_0;
    float4  _S2187 = _S2140;
    *&((&_S2187)->w) = _S2174;
    *&((&_S2187)->z) = _S2175;
    *&((&_S2187)->y) = _S2176;
    *&((&_S2187)->x) = _S2177;
    float4  _S2188 = _S2187;
    float3  _S2189 = _S2184.differential_0 + _S2064;
    *v_mean_1 = _S2189;
    *v_quat_1 = _S2188;
    *v_scale_1 = _S2180.differential_0;
    *v_in_opacity_1 = _S2111;
    (*v_sh_coeffs_1)[int(0)] = _S2089;
    (*v_sh_coeffs_1)[int(1)] = _S2090;
    (*v_sh_coeffs_1)[int(2)] = _S2091;
    (*v_sh_coeffs_1)[int(3)] = _S2092;
    (*v_sh_coeffs_1)[int(4)] = _S2093;
    (*v_sh_coeffs_1)[int(5)] = _S2094;
    (*v_sh_coeffs_1)[int(6)] = _S2095;
    (*v_sh_coeffs_1)[int(7)] = _S2096;
    (*v_sh_coeffs_1)[int(8)] = _S2097;
    (*v_sh_coeffs_1)[int(9)] = _S2098;
    (*v_sh_coeffs_1)[int(10)] = _S2099;
    (*v_sh_coeffs_1)[int(11)] = _S2100;
    (*v_sh_coeffs_1)[int(12)] = _S2101;
    (*v_sh_coeffs_1)[int(13)] = _S2102;
    (*v_sh_coeffs_1)[int(14)] = _S2103;
    (*v_sh_coeffs_1)[int(15)] = _S2104;
    *v_R_2 = _S2186;
    *v_t_2 = _S2185;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_21, float fy_21, float cx_17, float cy_17, FixedArray<float, 10>  * dist_coeffs_26, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S2190 = s_primal_ctx_exp_0(scale_14);
    float _S2191 = quat_15.y;
    float x2_15 = _S2191 * _S2191;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S2192 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2190.x, 0.0f, 0.0f, 0.0f, _S2190.y, 0.0f, 0.0f, 0.0f, _S2190.z);
    Matrix<float, 3, 3>  _S2193 = s_primal_ctx_mul_2(_S2192, S_2);
    Matrix<float, 3, 3>  _S2194 = transpose_0(_S2193);
    Matrix<float, 3, 3>  _S2195 = s_primal_ctx_mul_2(_S2193, _S2194);
    Matrix<float, 3, 3>  _S2196 = s_primal_ctx_mul_2(R_16, _S2195);
    Matrix<float, 3, 3>  _S2197 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S2198 = s_primal_ctx_mul_2(_S2196, _S2197);
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (fx_21, 0.0f, 0.0f, 0.0f, fy_21, 0.0f);
    Matrix<float, 2, 3>  _S2199 = s_primal_ctx_mul_3(J_13, _S2198);
    Matrix<float, 3, 2>  _S2200 = transpose_1(J_13);
    Matrix<float, 2, 2>  _S2201 = s_primal_ctx_mul_4(_S2199, _S2200);
    float _S2202 = _S2201.rows[int(0)].y * _S2201.rows[int(1)].x;
    float det_orig_13 = _S2201.rows[int(0)].x * _S2201.rows[int(1)].y - _S2202;
    float _S2203 = _S2201.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2204 = _S2201;
    *&(((&_S2204)->rows + (int(0)))->x) = _S2203;
    float _S2205 = _S2201.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2204)->rows + (int(1)))->y) = _S2205;
    Matrix<float, 2, 2>  _S2206 = _S2204;
    Matrix<float, 2, 2>  _S2207 = _S2204;
    float det_blur_8 = _S2203 * _S2205 - _S2202;
    float _S2208 = det_orig_13 / det_blur_8;
    float _S2209 = det_blur_8 * det_blur_8;
    float _S2210 = s_primal_ctx_max_0(0.0f, _S2208);
    float _S2211 = s_primal_ctx_sqrt_0(_S2210);
    float invdet_9 = 1.0f / det_blur_8;
    float _S2212 = - _S2201.rows[int(0)].y;
    float _S2213 = - _S2201.rows[int(1)].x;
    float _S2214 = - in_opacity_12;
    float _S2215 = 1.0f + s_primal_ctx_exp_1(_S2214);
    float _S2216 = 1.0f / _S2215;
    float _S2217 = _S2215 * _S2215;
    float _S2218;
    if(antialiased_12)
    {
        _S2218 = _S2216 * _S2211;
    }
    else
    {
        _S2218 = _S2216;
    }
    float _S2219 = _S2218 / 0.00392156885936856f;
    float _S2220 = 2.0f * s_primal_ctx_log_0(_S2219);
    float _S2221 = s_primal_ctx_sqrt_0(_S2220);
    float _S2222 = _S2206.rows[int(0)].x;
    float _S2223 = _S2207.rows[int(1)].y;
    float _S2224 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S2225 = mean_12 - - s_primal_ctx_mul_1(_S2197, t_15);
    float _S2226 = _S2225.x;
    float _S2227 = _S2225.y;
    float _S2228 = _S2225.z;
    float _S2229 = _S2226 * _S2226 + _S2227 * _S2227 + _S2228 * _S2228;
    float _S2230 = s_primal_ctx_sqrt_0(_S2229);
    float x_37 = _S2226 / _S2230;
    float3  _S2231 = make_float3 (x_37);
    float _S2232 = _S2230 * _S2230;
    float y_15 = _S2227 / _S2230;
    float z_12 = _S2228 / _S2230;
    float3  _S2233 = make_float3 (z_12);
    float _S2234 = - y_15;
    float3  _S2235 = make_float3 (_S2234);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_37 * x_37 - y_15 * y_15;
    float _S2236 = 2.0f * x_37;
    float fS1_12 = _S2236 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S2237 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_37;
    float3  _S2238 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S2239 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S2240 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S2241 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S2242 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S2242;
    float3  _S2243 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_37;
    float3  _S2244 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S2245 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S2246 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S2247 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_37 * fC1_12 - y_15 * fS1_12);
    float3  _S2248 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_37 * fS1_12 + y_15 * fC1_12);
    float3  _S2249 = make_float3 (pSH9_2);
    float3  _S2250 = make_float3 (0.0f);
    float3  _S2251 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2252;
    (&_S2252)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2234) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_37) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S2252)->differential_0 = _S2251;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2253;
    (&_S2253)->primal_0 = _S2250;
    (&_S2253)->differential_0 = _S2251;
    s_bwd_prop_max_0(&_S2252, &_S2253, v_rgb_2);
    float3  _S2254 = _S2248 * _S2252.differential_0;
    float3  _S2255 = (*sh_coeffs_12)[int(15)] * _S2252.differential_0;
    float3  _S2256 = _S2246 * _S2252.differential_0;
    float3  _S2257 = (*sh_coeffs_12)[int(14)] * _S2252.differential_0;
    float3  _S2258 = _S2244 * _S2252.differential_0;
    float3  _S2259 = (*sh_coeffs_12)[int(13)] * _S2252.differential_0;
    float3  _S2260 = _S2243 * _S2252.differential_0;
    float3  _S2261 = (*sh_coeffs_12)[int(12)] * _S2252.differential_0;
    float3  _S2262 = _S2245 * _S2252.differential_0;
    float3  _S2263 = (*sh_coeffs_12)[int(11)] * _S2252.differential_0;
    float3  _S2264 = _S2247 * _S2252.differential_0;
    float3  _S2265 = (*sh_coeffs_12)[int(10)] * _S2252.differential_0;
    float3  _S2266 = _S2249 * _S2252.differential_0;
    float3  _S2267 = (*sh_coeffs_12)[int(9)] * _S2252.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S2267.x + _S2267.y + _S2267.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S2255.x + _S2255.y + _S2255.z);
    float _S2268 = _S2265.x + _S2265.y + _S2265.z;
    float _S2269 = _S2257.x + _S2257.y + _S2257.z;
    float _S2270 = _S2263.x + _S2263.y + _S2263.z;
    float _S2271 = _S2259.x + _S2259.y + _S2259.z;
    float _S2272 = _S2261.x + _S2261.y + _S2261.z;
    float _S2273 = - s_diff_fC2_T_2;
    float3  _S2274 = _S2240 * _S2252.differential_0;
    float3  _S2275 = (*sh_coeffs_12)[int(8)] * _S2252.differential_0;
    float3  _S2276 = _S2238 * _S2252.differential_0;
    float3  _S2277 = (*sh_coeffs_12)[int(7)] * _S2252.differential_0;
    float3  _S2278 = _S2237 * _S2252.differential_0;
    float3  _S2279 = (*sh_coeffs_12)[int(6)] * _S2252.differential_0;
    float3  _S2280 = _S2239 * _S2252.differential_0;
    float3  _S2281 = (*sh_coeffs_12)[int(5)] * _S2252.differential_0;
    float3  _S2282 = _S2241 * _S2252.differential_0;
    float3  _S2283 = (*sh_coeffs_12)[int(4)] * _S2252.differential_0;
    float _S2284 = _S2281.x + _S2281.y + _S2281.z;
    float _S2285 = _S2277.x + _S2277.y + _S2277.z;
    float _S2286 = fTmp1B_12 * _S2268 + x_37 * s_diff_fS2_T_2 + y_15 * _S2273 + 0.54627424478530884f * (_S2283.x + _S2283.y + _S2283.z);
    float _S2287 = fTmp1B_12 * _S2269 + y_15 * s_diff_fS2_T_2 + x_37 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S2275.x + _S2275.y + _S2275.z);
    float _S2288 = y_15 * - _S2287;
    float _S2289 = x_37 * _S2287;
    float _S2290 = z_12 * (1.86588168144226074f * (z_12 * _S2272) + -2.28522896766662598f * (y_15 * _S2270 + x_37 * _S2271) + 0.94617468118667603f * (_S2279.x + _S2279.y + _S2279.z));
    float3  _S2291 = make_float3 (0.48860251903533936f) * _S2252.differential_0;
    float3  _S2292 = - _S2291;
    float3  _S2293 = _S2231 * _S2292;
    float3  _S2294 = (*sh_coeffs_12)[int(3)] * _S2292;
    float3  _S2295 = _S2233 * _S2291;
    float3  _S2296 = (*sh_coeffs_12)[int(2)] * _S2291;
    float3  _S2297 = _S2235 * _S2291;
    float3  _S2298 = (*sh_coeffs_12)[int(1)] * _S2291;
    float _S2299 = (_S2242 * _S2272 + 1.44530570507049561f * (fS1_12 * _S2268 + fC1_12 * _S2269) + -1.09254848957061768f * (y_15 * _S2284 + x_37 * _S2285) + _S2290 + _S2290 + _S2296.x + _S2296.y + _S2296.z) / _S2232;
    float _S2300 = _S2230 * _S2299;
    float _S2301 = (fTmp0C_12 * _S2270 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S2273 + fTmp0B_12 * _S2284 + _S2236 * _S2286 + _S2288 + _S2288 + - (_S2298.x + _S2298.y + _S2298.z)) / _S2232;
    float _S2302 = _S2230 * _S2301;
    float _S2303 = (fTmp0C_12 * _S2271 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S2285 + 2.0f * (y_15 * _S2286) + _S2289 + _S2289 + _S2294.x + _S2294.y + _S2294.z) / _S2232;
    float _S2304 = _S2230 * _S2303;
    float _S2305 = _S2228 * - _S2299 + _S2227 * - _S2301 + _S2226 * - _S2303;
    DiffPair_float_0 _S2306;
    (&_S2306)->primal_0 = _S2229;
    (&_S2306)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2306, _S2305);
    float _S2307 = _S2228 * _S2306.differential_0;
    float _S2308 = _S2227 * _S2306.differential_0;
    float _S2309 = _S2226 * _S2306.differential_0;
    float3  _S2310 = make_float3 (0.282094806432724f) * _S2252.differential_0;
    float3  _S2311 = make_float3 (_S2304 + _S2309 + _S2309, _S2302 + _S2308 + _S2308, _S2300 + _S2307 + _S2307);
    float3  _S2312 = - - _S2311;
    Matrix<float, 3, 3>  _S2313 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2314;
    (&_S2314)->primal_0 = _S2197;
    (&_S2314)->differential_0 = _S2313;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2315;
    (&_S2315)->primal_0 = t_15;
    (&_S2315)->differential_0 = _S2251;
    s_bwd_prop_mul_1(&_S2314, &_S2315, _S2312);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2316 = _S2314;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2317 = _S2315;
    float2  _S2318 = make_float2 (0.0f);
    float2  _S2319 = _S2318;
    *&((&_S2319)->y) = v_conic_2.z;
    float2  _S2320 = _S2318;
    *&((&_S2320)->y) = v_conic_2.y;
    *&((&_S2320)->x) = v_conic_2.x;
    float _S2321 = 0.5f * v_depth_2;
    DiffPair_float_0 _S2322;
    (&_S2322)->primal_0 = _S2224;
    (&_S2322)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2322, _S2321);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2323;
    (&_S2323)->primal_0 = mean_c_12;
    (&_S2323)->differential_0 = _S2251;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2324;
    (&_S2324)->primal_0 = mean_c_12;
    (&_S2324)->differential_0 = _S2251;
    s_bwd_prop_dot_0(&_S2323, &_S2324, _S2322.differential_0);
    DiffPair_float_0 _S2325;
    (&_S2325)->primal_0 = _S2223;
    (&_S2325)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2325, 0.0f);
    DiffPair_float_0 _S2326;
    (&_S2326)->primal_0 = _S2222;
    (&_S2326)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2326, 0.0f);
    DiffPair_float_0 _S2327;
    (&_S2327)->primal_0 = 3.32999992370605469f;
    (&_S2327)->differential_0 = 0.0f;
    DiffPair_float_0 _S2328;
    (&_S2328)->primal_0 = _S2221;
    (&_S2328)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2327, &_S2328, 0.0f);
    DiffPair_float_0 _S2329;
    (&_S2329)->primal_0 = _S2220;
    (&_S2329)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2329, _S2328.differential_0);
    float _S2330 = 2.0f * _S2329.differential_0;
    DiffPair_float_0 _S2331;
    (&_S2331)->primal_0 = _S2219;
    (&_S2331)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2331, _S2330);
    float _S2332 = v_opacity_2 + 254.9999847412109375f * _S2331.differential_0;
    Matrix<float, 2, 2>  _S2333 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2334 = _S2333;
    _S2334[int(1)] = _S2319;
    _S2334[int(0)] = _S2320;
    Matrix<float, 2, 2>  _S2335 = _S2334;
    FixedArray<float3 , 16>  _S2336;
    _S2336[int(0)] = _S2251;
    _S2336[int(1)] = _S2251;
    _S2336[int(2)] = _S2251;
    _S2336[int(3)] = _S2251;
    _S2336[int(4)] = _S2251;
    _S2336[int(5)] = _S2251;
    _S2336[int(6)] = _S2251;
    _S2336[int(7)] = _S2251;
    _S2336[int(8)] = _S2251;
    _S2336[int(9)] = _S2251;
    _S2336[int(10)] = _S2251;
    _S2336[int(11)] = _S2251;
    _S2336[int(12)] = _S2251;
    _S2336[int(13)] = _S2251;
    _S2336[int(14)] = _S2251;
    _S2336[int(15)] = _S2251;
    _S2336[int(7)] = _S2276;
    _S2336[int(0)] = _S2310;
    _S2336[int(1)] = _S2297;
    _S2336[int(2)] = _S2295;
    _S2336[int(3)] = _S2293;
    _S2336[int(4)] = _S2282;
    _S2336[int(5)] = _S2280;
    _S2336[int(6)] = _S2278;
    _S2336[int(15)] = _S2254;
    _S2336[int(8)] = _S2274;
    _S2336[int(9)] = _S2266;
    _S2336[int(10)] = _S2264;
    _S2336[int(11)] = _S2262;
    _S2336[int(12)] = _S2260;
    _S2336[int(13)] = _S2258;
    _S2336[int(14)] = _S2256;
    float3  _S2337 = _S2336[int(0)];
    float3  _S2338 = _S2336[int(1)];
    float3  _S2339 = _S2336[int(2)];
    float3  _S2340 = _S2336[int(3)];
    float3  _S2341 = _S2336[int(4)];
    float3  _S2342 = _S2336[int(5)];
    float3  _S2343 = _S2336[int(6)];
    float3  _S2344 = _S2336[int(7)];
    float3  _S2345 = _S2336[int(8)];
    float3  _S2346 = _S2336[int(9)];
    float3  _S2347 = _S2336[int(10)];
    float3  _S2348 = _S2336[int(11)];
    float3  _S2349 = _S2336[int(12)];
    float3  _S2350 = _S2336[int(13)];
    float3  _S2351 = _S2336[int(14)];
    float3  _S2352 = _S2336[int(15)];
    float3  _S2353 = _S2324.differential_0 + _S2323.differential_0;
    float2  _S2354 = make_float2 (0.0f, _S2325.differential_0);
    float2  _S2355 = make_float2 (_S2326.differential_0, 0.0f);
    float _S2356;
    if(antialiased_12)
    {
        float _S2357 = _S2216 * _S2332;
        _S2218 = _S2211 * _S2332;
        _S2356 = _S2357;
    }
    else
    {
        _S2218 = _S2332;
        _S2356 = 0.0f;
    }
    float _S2358 = - (_S2218 / _S2217);
    DiffPair_float_0 _S2359;
    (&_S2359)->primal_0 = _S2214;
    (&_S2359)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2359, _S2358);
    float _S2360 = - _S2359.differential_0;
    float _S2361 = invdet_9 * _S2335.rows[int(1)].y;
    float _S2362 = - (invdet_9 * _S2335.rows[int(1)].x);
    float _S2363 = - (invdet_9 * _S2335.rows[int(0)].y);
    float _S2364 = invdet_9 * _S2335.rows[int(0)].x;
    float _S2365 = - ((_S2203 * _S2335.rows[int(1)].y + _S2213 * _S2335.rows[int(1)].x + _S2212 * _S2335.rows[int(0)].y + _S2205 * _S2335.rows[int(0)].x) / _S2209);
    DiffPair_float_0 _S2366;
    (&_S2366)->primal_0 = _S2210;
    (&_S2366)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2366, _S2356);
    DiffPair_float_0 _S2367;
    (&_S2367)->primal_0 = 0.0f;
    (&_S2367)->differential_0 = 0.0f;
    DiffPair_float_0 _S2368;
    (&_S2368)->primal_0 = _S2208;
    (&_S2368)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2367, &_S2368, _S2366.differential_0);
    float _S2369 = _S2368.differential_0 / _S2209;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2369;
    float _S2370 = _S2365 + det_orig_13 * - _S2369;
    float _S2371 = - _S2370;
    float _S2372 = _S2203 * _S2370;
    float _S2373 = _S2205 * _S2370;
    Matrix<float, 2, 2>  _S2374 = _S2333;
    _S2374[int(1)] = _S2354;
    _S2374[int(0)] = _S2355;
    _S2204 = _S2374;
    *&(((&_S2204)->rows + (int(1)))->y) = 0.0f;
    float _S2375 = _S2364 + _S2372 + _S2374.rows[int(1)].y;
    *&(((&_S2204)->rows + (int(0)))->x) = 0.0f;
    float _S2376 = _S2361 + _S2373 + _S2374.rows[int(0)].x;
    float _S2377 = _S2371 + - s_diff_det_orig_T_2;
    float _S2378 = _S2362 + _S2201.rows[int(0)].y * _S2377;
    float _S2379 = _S2363 + _S2201.rows[int(1)].x * _S2377;
    float _S2380 = _S2201.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2381 = _S2375 + _S2201.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2382 = _S2318;
    *&((&_S2382)->x) = _S2378;
    *&((&_S2382)->y) = _S2381;
    float _S2383 = _S2376 + _S2380;
    float2  _S2384 = _S2318;
    *&((&_S2384)->y) = _S2379;
    *&((&_S2384)->x) = _S2383;
    float _S2385 = fy_21 * v_mean2d_2.y;
    float _S2386 = fx_21 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S2387 = _S2333;
    _S2387[int(1)] = _S2382;
    _S2387[int(0)] = _S2384;
    Matrix<float, 2, 2>  _S2388 = _S2204 + _S2387;
    Matrix<float, 2, 3>  _S2389 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2390;
    (&_S2390)->primal_0 = _S2199;
    (&_S2390)->differential_0 = _S2389;
    Matrix<float, 3, 2>  _S2391 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2392;
    (&_S2392)->primal_0 = _S2200;
    (&_S2392)->differential_0 = _S2391;
    s_bwd_prop_mul_2(&_S2390, &_S2392, _S2388);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2393;
    (&_S2393)->primal_0 = J_13;
    (&_S2393)->differential_0 = _S2389;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2394;
    (&_S2394)->primal_0 = _S2198;
    (&_S2394)->differential_0 = _S2313;
    s_bwd_prop_mul_3(&_S2393, &_S2394, _S2390.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2395;
    (&_S2395)->primal_0 = _S2196;
    (&_S2395)->differential_0 = _S2313;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2396;
    (&_S2396)->primal_0 = _S2197;
    (&_S2396)->differential_0 = _S2313;
    s_bwd_prop_mul_4(&_S2395, &_S2396, _S2394.differential_0);
    Matrix<float, 3, 3>  _S2397 = transpose_0(_S2396.differential_0 + _S2316.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2398;
    (&_S2398)->primal_0 = R_16;
    (&_S2398)->differential_0 = _S2313;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2399;
    (&_S2399)->primal_0 = _S2195;
    (&_S2399)->differential_0 = _S2313;
    s_bwd_prop_mul_4(&_S2398, &_S2399, _S2395.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2400;
    (&_S2400)->primal_0 = _S2193;
    (&_S2400)->differential_0 = _S2313;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2401;
    (&_S2401)->primal_0 = _S2194;
    (&_S2401)->differential_0 = _S2313;
    s_bwd_prop_mul_4(&_S2400, &_S2401, _S2399.differential_0);
    Matrix<float, 3, 3>  _S2402 = _S2400.differential_0 + transpose_0(_S2401.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2403;
    (&_S2403)->primal_0 = _S2192;
    (&_S2403)->differential_0 = _S2313;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2404;
    (&_S2404)->primal_0 = S_2;
    (&_S2404)->differential_0 = _S2313;
    s_bwd_prop_mul_4(&_S2403, &_S2404, _S2402);
    Matrix<float, 3, 3>  _S2405 = transpose_0(_S2403.differential_0);
    float _S2406 = 2.0f * - _S2405.rows[int(2)].z;
    float _S2407 = 2.0f * _S2405.rows[int(2)].y;
    float _S2408 = 2.0f * _S2405.rows[int(2)].x;
    float _S2409 = 2.0f * _S2405.rows[int(1)].z;
    float _S2410 = 2.0f * - _S2405.rows[int(1)].y;
    float _S2411 = 2.0f * _S2405.rows[int(1)].x;
    float _S2412 = 2.0f * _S2405.rows[int(0)].z;
    float _S2413 = 2.0f * _S2405.rows[int(0)].y;
    float _S2414 = 2.0f * - _S2405.rows[int(0)].x;
    float _S2415 = - _S2411 + _S2413;
    float _S2416 = _S2408 + - _S2412;
    float _S2417 = - _S2407 + _S2409;
    float _S2418 = _S2407 + _S2409;
    float _S2419 = _S2408 + _S2412;
    float _S2420 = _S2411 + _S2413;
    float _S2421 = quat_15.w * (_S2410 + _S2414);
    float _S2422 = quat_15.z * (_S2406 + _S2414);
    float _S2423 = quat_15.y * (_S2406 + _S2410);
    float _S2424 = quat_15.x * _S2415 + quat_15.z * _S2418 + quat_15.y * _S2419 + _S2421 + _S2421;
    float _S2425 = quat_15.x * _S2416 + quat_15.w * _S2418 + quat_15.y * _S2420 + _S2422 + _S2422;
    float _S2426 = quat_15.x * _S2417 + quat_15.w * _S2419 + quat_15.z * _S2420 + _S2423 + _S2423;
    float _S2427 = quat_15.w * _S2415 + quat_15.z * _S2416 + quat_15.y * _S2417;
    float3  _S2428 = _S2251;
    *&((&_S2428)->z) = _S2404.differential_0.rows[int(2)].z;
    *&((&_S2428)->y) = _S2404.differential_0.rows[int(1)].y;
    *&((&_S2428)->x) = _S2404.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2429;
    (&_S2429)->primal_0 = scale_14;
    (&_S2429)->differential_0 = _S2251;
    s_bwd_prop_exp_1(&_S2429, _S2428);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2430 = _S2429;
    float3  _S2431 = _S2251;
    *&((&_S2431)->y) = _S2385;
    *&((&_S2431)->x) = _S2386;
    float3  _S2432 = _S2353 + _S2431;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2433;
    (&_S2433)->primal_0 = R_16;
    (&_S2433)->differential_0 = _S2313;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2434;
    (&_S2434)->primal_0 = mean_12;
    (&_S2434)->differential_0 = _S2251;
    s_bwd_prop_mul_1(&_S2433, &_S2434, _S2432);
    float3  _S2435 = _S2432 + _S2317.differential_0;
    Matrix<float, 3, 3>  _S2436 = _S2397 + _S2398.differential_0 + _S2433.differential_0;
    float4  _S2437 = make_float4 (0.0f);
    *&((&_S2437)->w) = _S2424;
    *&((&_S2437)->z) = _S2425;
    *&((&_S2437)->y) = _S2426;
    *&((&_S2437)->x) = _S2427;
    float4  _S2438 = _S2437;
    float3  _S2439 = _S2434.differential_0 + _S2311;
    *v_mean_2 = _S2439;
    *v_quat_2 = _S2438;
    *v_scale_2 = _S2430.differential_0;
    *v_in_opacity_2 = _S2360;
    (*v_sh_coeffs_2)[int(0)] = _S2337;
    (*v_sh_coeffs_2)[int(1)] = _S2338;
    (*v_sh_coeffs_2)[int(2)] = _S2339;
    (*v_sh_coeffs_2)[int(3)] = _S2340;
    (*v_sh_coeffs_2)[int(4)] = _S2341;
    (*v_sh_coeffs_2)[int(5)] = _S2342;
    (*v_sh_coeffs_2)[int(6)] = _S2343;
    (*v_sh_coeffs_2)[int(7)] = _S2344;
    (*v_sh_coeffs_2)[int(8)] = _S2345;
    (*v_sh_coeffs_2)[int(9)] = _S2346;
    (*v_sh_coeffs_2)[int(10)] = _S2347;
    (*v_sh_coeffs_2)[int(11)] = _S2348;
    (*v_sh_coeffs_2)[int(12)] = _S2349;
    (*v_sh_coeffs_2)[int(13)] = _S2350;
    (*v_sh_coeffs_2)[int(14)] = _S2351;
    (*v_sh_coeffs_2)[int(15)] = _S2352;
    *v_R_3 = _S2436;
    *v_t_3 = _S2435;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2440;
};

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_22, float fy_22, float cx_18, float cy_18, FixedArray<float, 10>  * dist_coeffs_27, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 2, 2>  _S2441 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2442;
    (&_S2442)->_S2440 = _S2441;
    float2  _S2443 = make_float2 (0.0f);
    float3  _S2444 = make_float3 (0.0f);
    float4  intrins_12 = make_float4 (fx_22, fy_22, cx_18, cy_18);
    float3  _S2445 = s_primal_ctx_exp_0(scale_15);
    float _S2446 = quat_16.y;
    float x2_16 = _S2446 * _S2446;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S2447 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16))));
    FixedArray<float3 , 7>  _S2448 = {
        _S2444, _S2444, _S2444, _S2444, _S2444, _S2444, _S2444
    };
    FixedArray<float, 7>  _S2449 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2450;
    (&_S2450)->p_0 = _S2448;
    (&_S2450)->w_mean_0 = _S2449;
    (&_S2450)->w_cov_0 = _S2449;
    (&_S2450)->p_0[int(0)] = mean_13;
    SigmaPoints_0 _S2451 = _S2450;
    (&_S2451)->w_mean_0[int(0)] = 0.0f;
    (&_S2451)->w_cov_0[int(0)] = 2.0f;
    float _S2452 = s_primal_ctx_sqrt_0(3.0f);
    float _S2453 = _S2452 * _S2445.x;
    float3  delta_12 = make_float3 (_S2453) * _S2447.rows[0U];
    float3  _S2454 = mean_13 + delta_12;
    float3  _S2455 = mean_13 - delta_12;
    float _S2456 = _S2452 * _S2445.y;
    float3  delta_13 = make_float3 (_S2456) * _S2447.rows[1U];
    float3  _S2457 = mean_13 + delta_13;
    float3  _S2458 = mean_13 - delta_13;
    float _S2459 = _S2452 * _S2445.z;
    float3  delta_14 = make_float3 (_S2459) * _S2447.rows[2U];
    float3  _S2460 = mean_13 + delta_14;
    float3  _S2461 = mean_13 - delta_14;
    (&_S2451)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2451)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2451)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2451)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2451)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2451)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2451)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2462 = _S2450;
    (&_S2451)->p_0[0U] = s_primal_ctx_mul_1(R_17, _S2450.p_0[0U]) + t_16;
    (&_S2451)->p_0[1U] = s_primal_ctx_mul_1(R_17, _S2454) + t_16;
    (&_S2451)->p_0[2U] = s_primal_ctx_mul_1(R_17, _S2457) + t_16;
    (&_S2451)->p_0[3U] = s_primal_ctx_mul_1(R_17, _S2460) + t_16;
    (&_S2451)->p_0[4U] = s_primal_ctx_mul_1(R_17, _S2455) + t_16;
    (&_S2451)->p_0[5U] = s_primal_ctx_mul_1(R_17, _S2458) + t_16;
    (&_S2451)->p_0[6U] = s_primal_ctx_mul_1(R_17, _S2461) + t_16;
    float2  _S2463 = _S2443;
    Matrix<float, 2, 2>  _S2464 = _S2441;
    SigmaPoints_0 _S2465 = _S2451;
    bool _S2466 = persp_proj_3dgs_ut_1(&_S2465, intrins_12, dist_coeffs_27, image_width_13, image_height_13, &_S2464, &_S2463);
    (&_S2442)->_S2440 = _S2464;
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2467 = _S2442;
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S2468 = make_float3 (_S2453);
    float3  _S2469 = make_float3 (_S2456);
    float3  _S2470 = make_float3 (_S2459);
    float _S2471 = _S2442._S2440.rows[int(0)].y * _S2442._S2440.rows[int(1)].x;
    float det_orig_14 = _S2442._S2440.rows[int(0)].x * _S2442._S2440.rows[int(1)].y - _S2471;
    float _S2472 = _S2442._S2440.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2473 = _S2442._S2440;
    *&(((&_S2473)->rows + (int(0)))->x) = _S2472;
    float _S2474 = _S2442._S2440.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2473)->rows + (int(1)))->y) = _S2474;
    Matrix<float, 2, 2>  _S2475 = _S2473;
    Matrix<float, 2, 2>  _S2476 = _S2473;
    float det_blur_9 = _S2472 * _S2474 - _S2471;
    float _S2477 = det_orig_14 / det_blur_9;
    float _S2478 = det_blur_9 * det_blur_9;
    float _S2479 = s_primal_ctx_max_0(0.0f, _S2477);
    float _S2480 = s_primal_ctx_sqrt_0(_S2479);
    float _S2481 = - in_opacity_13;
    float _S2482 = 1.0f + s_primal_ctx_exp_1(_S2481);
    float _S2483 = 1.0f / _S2482;
    float _S2484 = _S2482 * _S2482;
    float _S2485;
    if(antialiased_13)
    {
        _S2485 = _S2483 * _S2480;
    }
    else
    {
        _S2485 = _S2483;
    }
    float _S2486 = _S2485 / 0.00392156885936856f;
    float _S2487 = 2.0f * s_primal_ctx_log_0(_S2486);
    float _S2488 = s_primal_ctx_sqrt_0(_S2487);
    float _S2489 = _S2475.rows[int(0)].x;
    float _S2490 = _S2476.rows[int(1)].y;
    float _S2491 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S2492 = - scale_15;
    Matrix<float, 3, 3>  _S2493 = transpose_0(R_17);
    float3  _S2494 = mean_13 - - s_primal_ctx_mul_1(_S2493, t_16);
    float _S2495 = _S2494.x;
    float _S2496 = _S2494.y;
    float _S2497 = _S2494.z;
    float _S2498 = _S2495 * _S2495 + _S2496 * _S2496 + _S2497 * _S2497;
    float _S2499 = s_primal_ctx_sqrt_0(_S2498);
    float x_38 = _S2495 / _S2499;
    float3  _S2500 = make_float3 (x_38);
    float _S2501 = _S2499 * _S2499;
    float y_16 = _S2496 / _S2499;
    float z_13 = _S2497 / _S2499;
    float3  _S2502 = make_float3 (z_13);
    float _S2503 = - y_16;
    float3  _S2504 = make_float3 (_S2503);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_38 * x_38 - y_16 * y_16;
    float _S2505 = 2.0f * x_38;
    float fS1_13 = _S2505 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S2506 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_38;
    float3  _S2507 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S2508 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2509 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2510 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S2511 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S2511;
    float3  _S2512 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_38;
    float3  _S2513 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S2514 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2515 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2516 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_38 * fC1_13 - y_16 * fS1_13);
    float3  _S2517 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_38 * fS1_13 + y_16 * fC1_13);
    float3  _S2518 = make_float3 (pSH9_3);
    float3  _S2519 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2520;
    (&_S2520)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2503) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_38) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S2520)->differential_0 = _S2444;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2521;
    (&_S2521)->primal_0 = _S2519;
    (&_S2521)->differential_0 = _S2444;
    s_bwd_prop_max_0(&_S2520, &_S2521, v_rgb_3);
    float3  _S2522 = _S2517 * _S2520.differential_0;
    float3  _S2523 = (*sh_coeffs_13)[int(15)] * _S2520.differential_0;
    float3  _S2524 = _S2515 * _S2520.differential_0;
    float3  _S2525 = (*sh_coeffs_13)[int(14)] * _S2520.differential_0;
    float3  _S2526 = _S2513 * _S2520.differential_0;
    float3  _S2527 = (*sh_coeffs_13)[int(13)] * _S2520.differential_0;
    float3  _S2528 = _S2512 * _S2520.differential_0;
    float3  _S2529 = (*sh_coeffs_13)[int(12)] * _S2520.differential_0;
    float3  _S2530 = _S2514 * _S2520.differential_0;
    float3  _S2531 = (*sh_coeffs_13)[int(11)] * _S2520.differential_0;
    float3  _S2532 = _S2516 * _S2520.differential_0;
    float3  _S2533 = (*sh_coeffs_13)[int(10)] * _S2520.differential_0;
    float3  _S2534 = _S2518 * _S2520.differential_0;
    float3  _S2535 = (*sh_coeffs_13)[int(9)] * _S2520.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2535.x + _S2535.y + _S2535.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2523.x + _S2523.y + _S2523.z);
    float _S2536 = _S2533.x + _S2533.y + _S2533.z;
    float _S2537 = _S2525.x + _S2525.y + _S2525.z;
    float _S2538 = _S2531.x + _S2531.y + _S2531.z;
    float _S2539 = _S2527.x + _S2527.y + _S2527.z;
    float _S2540 = _S2529.x + _S2529.y + _S2529.z;
    float _S2541 = - s_diff_fC2_T_3;
    float3  _S2542 = _S2509 * _S2520.differential_0;
    float3  _S2543 = (*sh_coeffs_13)[int(8)] * _S2520.differential_0;
    float3  _S2544 = _S2507 * _S2520.differential_0;
    float3  _S2545 = (*sh_coeffs_13)[int(7)] * _S2520.differential_0;
    float3  _S2546 = _S2506 * _S2520.differential_0;
    float3  _S2547 = (*sh_coeffs_13)[int(6)] * _S2520.differential_0;
    float3  _S2548 = _S2508 * _S2520.differential_0;
    float3  _S2549 = (*sh_coeffs_13)[int(5)] * _S2520.differential_0;
    float3  _S2550 = _S2510 * _S2520.differential_0;
    float3  _S2551 = (*sh_coeffs_13)[int(4)] * _S2520.differential_0;
    float _S2552 = _S2549.x + _S2549.y + _S2549.z;
    float _S2553 = _S2545.x + _S2545.y + _S2545.z;
    float _S2554 = fTmp1B_13 * _S2536 + x_38 * s_diff_fS2_T_3 + y_16 * _S2541 + 0.54627424478530884f * (_S2551.x + _S2551.y + _S2551.z);
    float _S2555 = fTmp1B_13 * _S2537 + y_16 * s_diff_fS2_T_3 + x_38 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2543.x + _S2543.y + _S2543.z);
    float _S2556 = y_16 * - _S2555;
    float _S2557 = x_38 * _S2555;
    float _S2558 = z_13 * (1.86588168144226074f * (z_13 * _S2540) + -2.28522896766662598f * (y_16 * _S2538 + x_38 * _S2539) + 0.94617468118667603f * (_S2547.x + _S2547.y + _S2547.z));
    float3  _S2559 = make_float3 (0.48860251903533936f) * _S2520.differential_0;
    float3  _S2560 = - _S2559;
    float3  _S2561 = _S2500 * _S2560;
    float3  _S2562 = (*sh_coeffs_13)[int(3)] * _S2560;
    float3  _S2563 = _S2502 * _S2559;
    float3  _S2564 = (*sh_coeffs_13)[int(2)] * _S2559;
    float3  _S2565 = _S2504 * _S2559;
    float3  _S2566 = (*sh_coeffs_13)[int(1)] * _S2559;
    float _S2567 = (_S2511 * _S2540 + 1.44530570507049561f * (fS1_13 * _S2536 + fC1_13 * _S2537) + -1.09254848957061768f * (y_16 * _S2552 + x_38 * _S2553) + _S2558 + _S2558 + _S2564.x + _S2564.y + _S2564.z) / _S2501;
    float _S2568 = _S2499 * _S2567;
    float _S2569 = (fTmp0C_13 * _S2538 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2541 + fTmp0B_13 * _S2552 + _S2505 * _S2554 + _S2556 + _S2556 + - (_S2566.x + _S2566.y + _S2566.z)) / _S2501;
    float _S2570 = _S2499 * _S2569;
    float _S2571 = (fTmp0C_13 * _S2539 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2553 + 2.0f * (y_16 * _S2554) + _S2557 + _S2557 + _S2562.x + _S2562.y + _S2562.z) / _S2501;
    float _S2572 = _S2499 * _S2571;
    float _S2573 = _S2497 * - _S2567 + _S2496 * - _S2569 + _S2495 * - _S2571;
    DiffPair_float_0 _S2574;
    (&_S2574)->primal_0 = _S2498;
    (&_S2574)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2574, _S2573);
    float _S2575 = _S2497 * _S2574.differential_0;
    float _S2576 = _S2496 * _S2574.differential_0;
    float _S2577 = _S2495 * _S2574.differential_0;
    float3  _S2578 = make_float3 (0.282094806432724f) * _S2520.differential_0;
    float3  _S2579 = make_float3 (_S2572 + _S2577 + _S2577, _S2570 + _S2576 + _S2576, _S2568 + _S2575 + _S2575);
    float3  _S2580 = - - _S2579;
    Matrix<float, 3, 3>  _S2581 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2582;
    (&_S2582)->primal_0 = _S2493;
    (&_S2582)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2583;
    (&_S2583)->primal_0 = t_16;
    (&_S2583)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2582, &_S2583, _S2580);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2584 = _S2583;
    Matrix<float, 3, 3>  _S2585 = transpose_0(_S2582.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2586;
    (&_S2586)->primal_0 = _S2492;
    (&_S2586)->differential_0 = _S2444;
    s_bwd_prop_exp_1(&_S2586, v_conic_3);
    float3  _S2587 = - _S2586.differential_0;
    float _S2588 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2589;
    (&_S2589)->primal_0 = _S2491;
    (&_S2589)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2589, _S2588);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2590;
    (&_S2590)->primal_0 = mean_c_13;
    (&_S2590)->differential_0 = _S2444;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2591;
    (&_S2591)->primal_0 = mean_c_13;
    (&_S2591)->differential_0 = _S2444;
    s_bwd_prop_dot_0(&_S2590, &_S2591, _S2589.differential_0);
    DiffPair_float_0 _S2592;
    (&_S2592)->primal_0 = _S2490;
    (&_S2592)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2592, 0.0f);
    DiffPair_float_0 _S2593;
    (&_S2593)->primal_0 = _S2489;
    (&_S2593)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2593, 0.0f);
    DiffPair_float_0 _S2594;
    (&_S2594)->primal_0 = 3.32999992370605469f;
    (&_S2594)->differential_0 = 0.0f;
    DiffPair_float_0 _S2595;
    (&_S2595)->primal_0 = _S2488;
    (&_S2595)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2594, &_S2595, 0.0f);
    DiffPair_float_0 _S2596;
    (&_S2596)->primal_0 = _S2487;
    (&_S2596)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2596, _S2595.differential_0);
    float _S2597 = 2.0f * _S2596.differential_0;
    DiffPair_float_0 _S2598;
    (&_S2598)->primal_0 = _S2486;
    (&_S2598)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2598, _S2597);
    float _S2599 = v_opacity_3 + 254.9999847412109375f * _S2598.differential_0;
    FixedArray<float3 , 16>  _S2600;
    _S2600[int(0)] = _S2444;
    _S2600[int(1)] = _S2444;
    _S2600[int(2)] = _S2444;
    _S2600[int(3)] = _S2444;
    _S2600[int(4)] = _S2444;
    _S2600[int(5)] = _S2444;
    _S2600[int(6)] = _S2444;
    _S2600[int(7)] = _S2444;
    _S2600[int(8)] = _S2444;
    _S2600[int(9)] = _S2444;
    _S2600[int(10)] = _S2444;
    _S2600[int(11)] = _S2444;
    _S2600[int(12)] = _S2444;
    _S2600[int(13)] = _S2444;
    _S2600[int(14)] = _S2444;
    _S2600[int(15)] = _S2444;
    _S2600[int(7)] = _S2544;
    _S2600[int(0)] = _S2578;
    _S2600[int(1)] = _S2565;
    _S2600[int(2)] = _S2563;
    _S2600[int(3)] = _S2561;
    _S2600[int(4)] = _S2550;
    _S2600[int(5)] = _S2548;
    _S2600[int(6)] = _S2546;
    _S2600[int(15)] = _S2522;
    _S2600[int(8)] = _S2542;
    _S2600[int(9)] = _S2534;
    _S2600[int(10)] = _S2532;
    _S2600[int(11)] = _S2530;
    _S2600[int(12)] = _S2528;
    _S2600[int(13)] = _S2526;
    _S2600[int(14)] = _S2524;
    float3  _S2601 = _S2600[int(0)];
    float3  _S2602 = _S2600[int(1)];
    float3  _S2603 = _S2600[int(2)];
    float3  _S2604 = _S2600[int(3)];
    float3  _S2605 = _S2600[int(4)];
    float3  _S2606 = _S2600[int(5)];
    float3  _S2607 = _S2600[int(6)];
    float3  _S2608 = _S2600[int(7)];
    float3  _S2609 = _S2600[int(8)];
    float3  _S2610 = _S2600[int(9)];
    float3  _S2611 = _S2600[int(10)];
    float3  _S2612 = _S2600[int(11)];
    float3  _S2613 = _S2600[int(12)];
    float3  _S2614 = _S2600[int(13)];
    float3  _S2615 = _S2600[int(14)];
    float3  _S2616 = _S2600[int(15)];
    float3  _S2617 = _S2591.differential_0 + _S2590.differential_0;
    float2  _S2618 = make_float2 (0.0f, _S2592.differential_0);
    float2  _S2619 = make_float2 (_S2593.differential_0, 0.0f);
    float _S2620;
    if(antialiased_13)
    {
        float _S2621 = _S2483 * _S2599;
        _S2485 = _S2480 * _S2599;
        _S2620 = _S2621;
    }
    else
    {
        _S2485 = _S2599;
        _S2620 = 0.0f;
    }
    float _S2622 = - (_S2485 / _S2484);
    DiffPair_float_0 _S2623;
    (&_S2623)->primal_0 = _S2481;
    (&_S2623)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2623, _S2622);
    float _S2624 = - _S2623.differential_0;
    DiffPair_float_0 _S2625;
    (&_S2625)->primal_0 = _S2479;
    (&_S2625)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2625, _S2620);
    DiffPair_float_0 _S2626;
    (&_S2626)->primal_0 = 0.0f;
    (&_S2626)->differential_0 = 0.0f;
    DiffPair_float_0 _S2627;
    (&_S2627)->primal_0 = _S2477;
    (&_S2627)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2626, &_S2627, _S2625.differential_0);
    float _S2628 = _S2627.differential_0 / _S2478;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2628;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2628;
    float _S2629 = - s_diff_det_blur_T_0;
    float _S2630 = _S2472 * s_diff_det_blur_T_0;
    float _S2631 = _S2474 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2632 = _S2441;
    _S2632[int(1)] = _S2618;
    _S2632[int(0)] = _S2619;
    float _S2633 = _S2631 + _S2632.rows[int(0)].x;
    float _S2634 = _S2629 + - s_diff_det_orig_T_3;
    float _S2635 = _S2467._S2440.rows[int(0)].y * _S2634;
    float _S2636 = _S2467._S2440.rows[int(1)].x * _S2634;
    float _S2637 = _S2467._S2440.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2638 = _S2630 + _S2632.rows[int(1)].y + _S2467._S2440.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2639 = _S2443;
    *&((&_S2639)->x) = _S2635;
    *&((&_S2639)->y) = _S2638;
    float _S2640 = _S2633 + _S2637;
    float2  _S2641 = _S2443;
    *&((&_S2641)->y) = _S2636;
    *&((&_S2641)->x) = _S2640;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2642;
    (&_S2642)->primal_0 = R_17;
    (&_S2642)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2643;
    (&_S2643)->primal_0 = _S2461;
    (&_S2643)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2642, &_S2643, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2644;
    (&_S2644)->primal_0 = R_17;
    (&_S2644)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2645;
    (&_S2645)->primal_0 = _S2458;
    (&_S2645)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2644, &_S2645, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2646;
    (&_S2646)->primal_0 = R_17;
    (&_S2646)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2647;
    (&_S2647)->primal_0 = _S2455;
    (&_S2647)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2646, &_S2647, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2648;
    (&_S2648)->primal_0 = R_17;
    (&_S2648)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2649;
    (&_S2649)->primal_0 = _S2460;
    (&_S2649)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2648, &_S2649, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2650;
    (&_S2650)->primal_0 = R_17;
    (&_S2650)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2651;
    (&_S2651)->primal_0 = _S2457;
    (&_S2651)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2650, &_S2651, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2652;
    (&_S2652)->primal_0 = R_17;
    (&_S2652)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2653;
    (&_S2653)->primal_0 = _S2454;
    (&_S2653)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2652, &_S2653, _S2444);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2654;
    (&_S2654)->primal_0 = R_17;
    (&_S2654)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2655;
    (&_S2655)->primal_0 = _S2462.p_0[0U];
    (&_S2655)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2654, &_S2655, _S2444);
    float3  _S2656 = - _S2643.differential_0 + _S2649.differential_0;
    float3  _S2657 = _S2470 * _S2656;
    float3  _S2658 = _S2447.rows[2U] * _S2656;
    float _S2659 = _S2452 * (_S2658.x + _S2658.y + _S2658.z);
    float3  _S2660 = - _S2645.differential_0 + _S2651.differential_0;
    float3  _S2661 = _S2469 * _S2660;
    float3  _S2662 = _S2447.rows[1U] * _S2660;
    float _S2663 = _S2452 * (_S2662.x + _S2662.y + _S2662.z);
    float3  _S2664 = - _S2647.differential_0 + _S2653.differential_0;
    float3  _S2665 = _S2468 * _S2664;
    float3  _S2666 = _S2447.rows[0U] * _S2664;
    float _S2667 = _S2452 * (_S2666.x + _S2666.y + _S2666.z);
    Matrix<float, 3, 3>  _S2668 = _S2581;
    _S2668[2U] = _S2657;
    _S2668[1U] = _S2661;
    _S2668[0U] = _S2665;
    Matrix<float, 3, 3>  _S2669 = transpose_0(transpose_0(_S2668));
    float _S2670 = 2.0f * - _S2669.rows[int(2)].z;
    float _S2671 = 2.0f * _S2669.rows[int(2)].y;
    float _S2672 = 2.0f * _S2669.rows[int(2)].x;
    float _S2673 = 2.0f * _S2669.rows[int(1)].z;
    float _S2674 = 2.0f * - _S2669.rows[int(1)].y;
    float _S2675 = 2.0f * _S2669.rows[int(1)].x;
    float _S2676 = 2.0f * _S2669.rows[int(0)].z;
    float _S2677 = 2.0f * _S2669.rows[int(0)].y;
    float _S2678 = 2.0f * - _S2669.rows[int(0)].x;
    float _S2679 = - _S2675 + _S2677;
    float _S2680 = _S2672 + - _S2676;
    float _S2681 = - _S2671 + _S2673;
    float _S2682 = _S2671 + _S2673;
    float _S2683 = _S2672 + _S2676;
    float _S2684 = _S2675 + _S2677;
    float _S2685 = quat_16.w * (_S2674 + _S2678);
    float _S2686 = quat_16.z * (_S2670 + _S2678);
    float _S2687 = quat_16.y * (_S2670 + _S2674);
    float _S2688 = quat_16.x * _S2679 + quat_16.z * _S2682 + quat_16.y * _S2683 + _S2685 + _S2685;
    float _S2689 = quat_16.x * _S2680 + quat_16.w * _S2682 + quat_16.y * _S2684 + _S2686 + _S2686;
    float _S2690 = quat_16.x * _S2681 + quat_16.w * _S2683 + quat_16.z * _S2684 + _S2687 + _S2687;
    float _S2691 = quat_16.w * _S2679 + quat_16.z * _S2680 + quat_16.y * _S2681;
    float3  _S2692 = _S2444;
    *&((&_S2692)->z) = _S2659;
    *&((&_S2692)->y) = _S2663;
    *&((&_S2692)->x) = _S2667;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2693;
    (&_S2693)->primal_0 = scale_15;
    (&_S2693)->differential_0 = _S2444;
    s_bwd_prop_exp_1(&_S2693, _S2692);
    Matrix<float, 2, 2>  _S2694 = _S2441;
    _S2694[int(1)] = _S2639;
    _S2694[int(0)] = _S2641;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2695;
    (&_S2695)->primal_0 = R_17;
    (&_S2695)->differential_0 = _S2581;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2696;
    (&_S2696)->primal_0 = mean_13;
    (&_S2696)->differential_0 = _S2444;
    s_bwd_prop_mul_1(&_S2695, &_S2696, _S2617);
    float3  _S2697 = _S2617 + _S2584.differential_0;
    Matrix<float, 3, 3>  _S2698 = _S2642.differential_0 + _S2644.differential_0 + _S2646.differential_0 + _S2648.differential_0 + _S2650.differential_0 + _S2652.differential_0 + _S2654.differential_0 + _S2695.differential_0 + _S2585;
    float3  _S2699 = _S2693.differential_0 + _S2587;
    float4  _S2700 = make_float4 (0.0f);
    *&((&_S2700)->w) = _S2688;
    *&((&_S2700)->z) = _S2689;
    *&((&_S2700)->y) = _S2690;
    *&((&_S2700)->x) = _S2691;
    float4  _S2701 = _S2700;
    float3  _S2702 = _S2643.differential_0 + _S2649.differential_0 + _S2645.differential_0 + _S2651.differential_0 + _S2647.differential_0 + _S2653.differential_0 + _S2696.differential_0 + _S2579;
    *v_mean_3 = _S2702;
    *v_quat_3 = _S2701;
    *v_scale_3 = _S2699;
    *v_in_opacity_3 = _S2624;
    (*v_sh_coeffs_3)[int(0)] = _S2601;
    (*v_sh_coeffs_3)[int(1)] = _S2602;
    (*v_sh_coeffs_3)[int(2)] = _S2603;
    (*v_sh_coeffs_3)[int(3)] = _S2604;
    (*v_sh_coeffs_3)[int(4)] = _S2605;
    (*v_sh_coeffs_3)[int(5)] = _S2606;
    (*v_sh_coeffs_3)[int(6)] = _S2607;
    (*v_sh_coeffs_3)[int(7)] = _S2608;
    (*v_sh_coeffs_3)[int(8)] = _S2609;
    (*v_sh_coeffs_3)[int(9)] = _S2610;
    (*v_sh_coeffs_3)[int(10)] = _S2611;
    (*v_sh_coeffs_3)[int(11)] = _S2612;
    (*v_sh_coeffs_3)[int(12)] = _S2613;
    (*v_sh_coeffs_3)[int(13)] = _S2614;
    (*v_sh_coeffs_3)[int(14)] = _S2615;
    (*v_sh_coeffs_3)[int(15)] = _S2616;
    *v_R_4 = _S2698;
    *v_t_4 = _S2697;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2703;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_23, float fy_23, float cx_19, float cy_19, FixedArray<float, 10>  * dist_coeffs_28, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2704 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2705;
    (&_S2705)->_S2703 = _S2704;
    float2  _S2706 = make_float2 (0.0f);
    float3  _S2707 = make_float3 (0.0f);
    float4  intrins_13 = make_float4 (fx_23, fy_23, cx_19, cy_19);
    float3  _S2708 = s_primal_ctx_exp_0(scale_16);
    float _S2709 = quat_17.y;
    float x2_17 = _S2709 * _S2709;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2710 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17))));
    FixedArray<float3 , 7>  _S2711 = {
        _S2707, _S2707, _S2707, _S2707, _S2707, _S2707, _S2707
    };
    FixedArray<float, 7>  _S2712 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2713;
    (&_S2713)->p_0 = _S2711;
    (&_S2713)->w_mean_0 = _S2712;
    (&_S2713)->w_cov_0 = _S2712;
    (&_S2713)->p_0[int(0)] = mean_14;
    SigmaPoints_0 _S2714 = _S2713;
    (&_S2714)->w_mean_0[int(0)] = 0.0f;
    (&_S2714)->w_cov_0[int(0)] = 2.0f;
    float _S2715 = s_primal_ctx_sqrt_0(3.0f);
    float _S2716 = _S2715 * _S2708.x;
    float3  delta_15 = make_float3 (_S2716) * _S2710.rows[0U];
    float3  _S2717 = mean_14 + delta_15;
    float3  _S2718 = mean_14 - delta_15;
    float _S2719 = _S2715 * _S2708.y;
    float3  delta_16 = make_float3 (_S2719) * _S2710.rows[1U];
    float3  _S2720 = mean_14 + delta_16;
    float3  _S2721 = mean_14 - delta_16;
    float _S2722 = _S2715 * _S2708.z;
    float3  delta_17 = make_float3 (_S2722) * _S2710.rows[2U];
    float3  _S2723 = mean_14 + delta_17;
    float3  _S2724 = mean_14 - delta_17;
    (&_S2714)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2714)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2714)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2714)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2714)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2714)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2714)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2725 = _S2713;
    (&_S2714)->p_0[0U] = s_primal_ctx_mul_1(R_18, _S2713.p_0[0U]) + t_17;
    (&_S2714)->p_0[1U] = s_primal_ctx_mul_1(R_18, _S2717) + t_17;
    (&_S2714)->p_0[2U] = s_primal_ctx_mul_1(R_18, _S2720) + t_17;
    (&_S2714)->p_0[3U] = s_primal_ctx_mul_1(R_18, _S2723) + t_17;
    (&_S2714)->p_0[4U] = s_primal_ctx_mul_1(R_18, _S2718) + t_17;
    (&_S2714)->p_0[5U] = s_primal_ctx_mul_1(R_18, _S2721) + t_17;
    (&_S2714)->p_0[6U] = s_primal_ctx_mul_1(R_18, _S2724) + t_17;
    float2  _S2726 = _S2706;
    Matrix<float, 2, 2>  _S2727 = _S2704;
    SigmaPoints_0 _S2728 = _S2714;
    bool _S2729 = fisheye_proj_3dgs_ut_1(&_S2728, intrins_13, dist_coeffs_28, &_S2727, &_S2726);
    (&_S2705)->_S2703 = _S2727;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2730 = _S2705;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2731 = make_float3 (_S2716);
    float3  _S2732 = make_float3 (_S2719);
    float3  _S2733 = make_float3 (_S2722);
    float _S2734 = _S2705._S2703.rows[int(0)].y * _S2705._S2703.rows[int(1)].x;
    float det_orig_15 = _S2705._S2703.rows[int(0)].x * _S2705._S2703.rows[int(1)].y - _S2734;
    float _S2735 = _S2705._S2703.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2736 = _S2705._S2703;
    *&(((&_S2736)->rows + (int(0)))->x) = _S2735;
    float _S2737 = _S2705._S2703.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2736)->rows + (int(1)))->y) = _S2737;
    Matrix<float, 2, 2>  _S2738 = _S2736;
    Matrix<float, 2, 2>  _S2739 = _S2736;
    float det_blur_10 = _S2735 * _S2737 - _S2734;
    float _S2740 = det_orig_15 / det_blur_10;
    float _S2741 = det_blur_10 * det_blur_10;
    float _S2742 = s_primal_ctx_max_0(0.0f, _S2740);
    float _S2743 = s_primal_ctx_sqrt_0(_S2742);
    float _S2744 = - in_opacity_14;
    float _S2745 = 1.0f + s_primal_ctx_exp_1(_S2744);
    float _S2746 = 1.0f / _S2745;
    float _S2747 = _S2745 * _S2745;
    float _S2748;
    if(antialiased_14)
    {
        _S2748 = _S2746 * _S2743;
    }
    else
    {
        _S2748 = _S2746;
    }
    float _S2749 = _S2748 / 0.00392156885936856f;
    float _S2750 = 2.0f * s_primal_ctx_log_0(_S2749);
    float _S2751 = s_primal_ctx_sqrt_0(_S2750);
    float _S2752 = _S2738.rows[int(0)].x;
    float _S2753 = _S2739.rows[int(1)].y;
    float _S2754 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2755 = - scale_16;
    Matrix<float, 3, 3>  _S2756 = transpose_0(R_18);
    float3  _S2757 = mean_14 - - s_primal_ctx_mul_1(_S2756, t_17);
    float _S2758 = _S2757.x;
    float _S2759 = _S2757.y;
    float _S2760 = _S2757.z;
    float _S2761 = _S2758 * _S2758 + _S2759 * _S2759 + _S2760 * _S2760;
    float _S2762 = s_primal_ctx_sqrt_0(_S2761);
    float x_39 = _S2758 / _S2762;
    float3  _S2763 = make_float3 (x_39);
    float _S2764 = _S2762 * _S2762;
    float y_17 = _S2759 / _S2762;
    float z_14 = _S2760 / _S2762;
    float3  _S2765 = make_float3 (z_14);
    float _S2766 = - y_17;
    float3  _S2767 = make_float3 (_S2766);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_39 * x_39 - y_17 * y_17;
    float _S2768 = 2.0f * x_39;
    float fS1_14 = _S2768 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2769 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_39;
    float3  _S2770 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2771 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2772 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2773 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2774 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2774;
    float3  _S2775 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_39;
    float3  _S2776 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2777 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2778 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2779 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_39 * fC1_14 - y_17 * fS1_14);
    float3  _S2780 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_39 * fS1_14 + y_17 * fC1_14);
    float3  _S2781 = make_float3 (pSH9_4);
    float3  _S2782 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2783;
    (&_S2783)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2766) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_39) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2783)->differential_0 = _S2707;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2784;
    (&_S2784)->primal_0 = _S2782;
    (&_S2784)->differential_0 = _S2707;
    s_bwd_prop_max_0(&_S2783, &_S2784, v_rgb_4);
    float3  _S2785 = _S2780 * _S2783.differential_0;
    float3  _S2786 = (*sh_coeffs_14)[int(15)] * _S2783.differential_0;
    float3  _S2787 = _S2778 * _S2783.differential_0;
    float3  _S2788 = (*sh_coeffs_14)[int(14)] * _S2783.differential_0;
    float3  _S2789 = _S2776 * _S2783.differential_0;
    float3  _S2790 = (*sh_coeffs_14)[int(13)] * _S2783.differential_0;
    float3  _S2791 = _S2775 * _S2783.differential_0;
    float3  _S2792 = (*sh_coeffs_14)[int(12)] * _S2783.differential_0;
    float3  _S2793 = _S2777 * _S2783.differential_0;
    float3  _S2794 = (*sh_coeffs_14)[int(11)] * _S2783.differential_0;
    float3  _S2795 = _S2779 * _S2783.differential_0;
    float3  _S2796 = (*sh_coeffs_14)[int(10)] * _S2783.differential_0;
    float3  _S2797 = _S2781 * _S2783.differential_0;
    float3  _S2798 = (*sh_coeffs_14)[int(9)] * _S2783.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2798.x + _S2798.y + _S2798.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2786.x + _S2786.y + _S2786.z);
    float _S2799 = _S2796.x + _S2796.y + _S2796.z;
    float _S2800 = _S2788.x + _S2788.y + _S2788.z;
    float _S2801 = _S2794.x + _S2794.y + _S2794.z;
    float _S2802 = _S2790.x + _S2790.y + _S2790.z;
    float _S2803 = _S2792.x + _S2792.y + _S2792.z;
    float _S2804 = - s_diff_fC2_T_4;
    float3  _S2805 = _S2772 * _S2783.differential_0;
    float3  _S2806 = (*sh_coeffs_14)[int(8)] * _S2783.differential_0;
    float3  _S2807 = _S2770 * _S2783.differential_0;
    float3  _S2808 = (*sh_coeffs_14)[int(7)] * _S2783.differential_0;
    float3  _S2809 = _S2769 * _S2783.differential_0;
    float3  _S2810 = (*sh_coeffs_14)[int(6)] * _S2783.differential_0;
    float3  _S2811 = _S2771 * _S2783.differential_0;
    float3  _S2812 = (*sh_coeffs_14)[int(5)] * _S2783.differential_0;
    float3  _S2813 = _S2773 * _S2783.differential_0;
    float3  _S2814 = (*sh_coeffs_14)[int(4)] * _S2783.differential_0;
    float _S2815 = _S2812.x + _S2812.y + _S2812.z;
    float _S2816 = _S2808.x + _S2808.y + _S2808.z;
    float _S2817 = fTmp1B_14 * _S2799 + x_39 * s_diff_fS2_T_4 + y_17 * _S2804 + 0.54627424478530884f * (_S2814.x + _S2814.y + _S2814.z);
    float _S2818 = fTmp1B_14 * _S2800 + y_17 * s_diff_fS2_T_4 + x_39 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2806.x + _S2806.y + _S2806.z);
    float _S2819 = y_17 * - _S2818;
    float _S2820 = x_39 * _S2818;
    float _S2821 = z_14 * (1.86588168144226074f * (z_14 * _S2803) + -2.28522896766662598f * (y_17 * _S2801 + x_39 * _S2802) + 0.94617468118667603f * (_S2810.x + _S2810.y + _S2810.z));
    float3  _S2822 = make_float3 (0.48860251903533936f) * _S2783.differential_0;
    float3  _S2823 = - _S2822;
    float3  _S2824 = _S2763 * _S2823;
    float3  _S2825 = (*sh_coeffs_14)[int(3)] * _S2823;
    float3  _S2826 = _S2765 * _S2822;
    float3  _S2827 = (*sh_coeffs_14)[int(2)] * _S2822;
    float3  _S2828 = _S2767 * _S2822;
    float3  _S2829 = (*sh_coeffs_14)[int(1)] * _S2822;
    float _S2830 = (_S2774 * _S2803 + 1.44530570507049561f * (fS1_14 * _S2799 + fC1_14 * _S2800) + -1.09254848957061768f * (y_17 * _S2815 + x_39 * _S2816) + _S2821 + _S2821 + _S2827.x + _S2827.y + _S2827.z) / _S2764;
    float _S2831 = _S2762 * _S2830;
    float _S2832 = (fTmp0C_14 * _S2801 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2804 + fTmp0B_14 * _S2815 + _S2768 * _S2817 + _S2819 + _S2819 + - (_S2829.x + _S2829.y + _S2829.z)) / _S2764;
    float _S2833 = _S2762 * _S2832;
    float _S2834 = (fTmp0C_14 * _S2802 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2816 + 2.0f * (y_17 * _S2817) + _S2820 + _S2820 + _S2825.x + _S2825.y + _S2825.z) / _S2764;
    float _S2835 = _S2762 * _S2834;
    float _S2836 = _S2760 * - _S2830 + _S2759 * - _S2832 + _S2758 * - _S2834;
    DiffPair_float_0 _S2837;
    (&_S2837)->primal_0 = _S2761;
    (&_S2837)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2837, _S2836);
    float _S2838 = _S2760 * _S2837.differential_0;
    float _S2839 = _S2759 * _S2837.differential_0;
    float _S2840 = _S2758 * _S2837.differential_0;
    float3  _S2841 = make_float3 (0.282094806432724f) * _S2783.differential_0;
    float3  _S2842 = make_float3 (_S2835 + _S2840 + _S2840, _S2833 + _S2839 + _S2839, _S2831 + _S2838 + _S2838);
    float3  _S2843 = - - _S2842;
    Matrix<float, 3, 3>  _S2844 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2845;
    (&_S2845)->primal_0 = _S2756;
    (&_S2845)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2846;
    (&_S2846)->primal_0 = t_17;
    (&_S2846)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2845, &_S2846, _S2843);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2847 = _S2846;
    Matrix<float, 3, 3>  _S2848 = transpose_0(_S2845.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2849;
    (&_S2849)->primal_0 = _S2755;
    (&_S2849)->differential_0 = _S2707;
    s_bwd_prop_exp_1(&_S2849, v_conic_4);
    float3  _S2850 = - _S2849.differential_0;
    float _S2851 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2852;
    (&_S2852)->primal_0 = _S2754;
    (&_S2852)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2852, _S2851);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2853;
    (&_S2853)->primal_0 = mean_c_14;
    (&_S2853)->differential_0 = _S2707;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2854;
    (&_S2854)->primal_0 = mean_c_14;
    (&_S2854)->differential_0 = _S2707;
    s_bwd_prop_dot_0(&_S2853, &_S2854, _S2852.differential_0);
    DiffPair_float_0 _S2855;
    (&_S2855)->primal_0 = _S2753;
    (&_S2855)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2855, 0.0f);
    DiffPair_float_0 _S2856;
    (&_S2856)->primal_0 = _S2752;
    (&_S2856)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2856, 0.0f);
    DiffPair_float_0 _S2857;
    (&_S2857)->primal_0 = 3.32999992370605469f;
    (&_S2857)->differential_0 = 0.0f;
    DiffPair_float_0 _S2858;
    (&_S2858)->primal_0 = _S2751;
    (&_S2858)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2857, &_S2858, 0.0f);
    DiffPair_float_0 _S2859;
    (&_S2859)->primal_0 = _S2750;
    (&_S2859)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2859, _S2858.differential_0);
    float _S2860 = 2.0f * _S2859.differential_0;
    DiffPair_float_0 _S2861;
    (&_S2861)->primal_0 = _S2749;
    (&_S2861)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2861, _S2860);
    float _S2862 = v_opacity_4 + 254.9999847412109375f * _S2861.differential_0;
    FixedArray<float3 , 16>  _S2863;
    _S2863[int(0)] = _S2707;
    _S2863[int(1)] = _S2707;
    _S2863[int(2)] = _S2707;
    _S2863[int(3)] = _S2707;
    _S2863[int(4)] = _S2707;
    _S2863[int(5)] = _S2707;
    _S2863[int(6)] = _S2707;
    _S2863[int(7)] = _S2707;
    _S2863[int(8)] = _S2707;
    _S2863[int(9)] = _S2707;
    _S2863[int(10)] = _S2707;
    _S2863[int(11)] = _S2707;
    _S2863[int(12)] = _S2707;
    _S2863[int(13)] = _S2707;
    _S2863[int(14)] = _S2707;
    _S2863[int(15)] = _S2707;
    _S2863[int(7)] = _S2807;
    _S2863[int(0)] = _S2841;
    _S2863[int(1)] = _S2828;
    _S2863[int(2)] = _S2826;
    _S2863[int(3)] = _S2824;
    _S2863[int(4)] = _S2813;
    _S2863[int(5)] = _S2811;
    _S2863[int(6)] = _S2809;
    _S2863[int(15)] = _S2785;
    _S2863[int(8)] = _S2805;
    _S2863[int(9)] = _S2797;
    _S2863[int(10)] = _S2795;
    _S2863[int(11)] = _S2793;
    _S2863[int(12)] = _S2791;
    _S2863[int(13)] = _S2789;
    _S2863[int(14)] = _S2787;
    float3  _S2864 = _S2863[int(0)];
    float3  _S2865 = _S2863[int(1)];
    float3  _S2866 = _S2863[int(2)];
    float3  _S2867 = _S2863[int(3)];
    float3  _S2868 = _S2863[int(4)];
    float3  _S2869 = _S2863[int(5)];
    float3  _S2870 = _S2863[int(6)];
    float3  _S2871 = _S2863[int(7)];
    float3  _S2872 = _S2863[int(8)];
    float3  _S2873 = _S2863[int(9)];
    float3  _S2874 = _S2863[int(10)];
    float3  _S2875 = _S2863[int(11)];
    float3  _S2876 = _S2863[int(12)];
    float3  _S2877 = _S2863[int(13)];
    float3  _S2878 = _S2863[int(14)];
    float3  _S2879 = _S2863[int(15)];
    float3  _S2880 = _S2854.differential_0 + _S2853.differential_0;
    float2  _S2881 = make_float2 (0.0f, _S2855.differential_0);
    float2  _S2882 = make_float2 (_S2856.differential_0, 0.0f);
    float _S2883;
    if(antialiased_14)
    {
        float _S2884 = _S2746 * _S2862;
        _S2748 = _S2743 * _S2862;
        _S2883 = _S2884;
    }
    else
    {
        _S2748 = _S2862;
        _S2883 = 0.0f;
    }
    float _S2885 = - (_S2748 / _S2747);
    DiffPair_float_0 _S2886;
    (&_S2886)->primal_0 = _S2744;
    (&_S2886)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2886, _S2885);
    float _S2887 = - _S2886.differential_0;
    DiffPair_float_0 _S2888;
    (&_S2888)->primal_0 = _S2742;
    (&_S2888)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2888, _S2883);
    DiffPair_float_0 _S2889;
    (&_S2889)->primal_0 = 0.0f;
    (&_S2889)->differential_0 = 0.0f;
    DiffPair_float_0 _S2890;
    (&_S2890)->primal_0 = _S2740;
    (&_S2890)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2889, &_S2890, _S2888.differential_0);
    float _S2891 = _S2890.differential_0 / _S2741;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2891;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2891;
    float _S2892 = - s_diff_det_blur_T_1;
    float _S2893 = _S2735 * s_diff_det_blur_T_1;
    float _S2894 = _S2737 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2895 = _S2704;
    _S2895[int(1)] = _S2881;
    _S2895[int(0)] = _S2882;
    float _S2896 = _S2894 + _S2895.rows[int(0)].x;
    float _S2897 = _S2892 + - s_diff_det_orig_T_4;
    float _S2898 = _S2730._S2703.rows[int(0)].y * _S2897;
    float _S2899 = _S2730._S2703.rows[int(1)].x * _S2897;
    float _S2900 = _S2730._S2703.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2901 = _S2893 + _S2895.rows[int(1)].y + _S2730._S2703.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2902 = _S2706;
    *&((&_S2902)->x) = _S2898;
    *&((&_S2902)->y) = _S2901;
    float _S2903 = _S2896 + _S2900;
    float2  _S2904 = _S2706;
    *&((&_S2904)->y) = _S2899;
    *&((&_S2904)->x) = _S2903;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2905;
    (&_S2905)->primal_0 = R_18;
    (&_S2905)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2906;
    (&_S2906)->primal_0 = _S2724;
    (&_S2906)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2905, &_S2906, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2907;
    (&_S2907)->primal_0 = R_18;
    (&_S2907)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2908;
    (&_S2908)->primal_0 = _S2721;
    (&_S2908)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2907, &_S2908, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2909;
    (&_S2909)->primal_0 = R_18;
    (&_S2909)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2910;
    (&_S2910)->primal_0 = _S2718;
    (&_S2910)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2909, &_S2910, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2911;
    (&_S2911)->primal_0 = R_18;
    (&_S2911)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2912;
    (&_S2912)->primal_0 = _S2723;
    (&_S2912)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2911, &_S2912, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2913;
    (&_S2913)->primal_0 = R_18;
    (&_S2913)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2914;
    (&_S2914)->primal_0 = _S2720;
    (&_S2914)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2913, &_S2914, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2915;
    (&_S2915)->primal_0 = R_18;
    (&_S2915)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2916;
    (&_S2916)->primal_0 = _S2717;
    (&_S2916)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2915, &_S2916, _S2707);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2917;
    (&_S2917)->primal_0 = R_18;
    (&_S2917)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2918;
    (&_S2918)->primal_0 = _S2725.p_0[0U];
    (&_S2918)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2917, &_S2918, _S2707);
    float3  _S2919 = - _S2906.differential_0 + _S2912.differential_0;
    float3  _S2920 = _S2733 * _S2919;
    float3  _S2921 = _S2710.rows[2U] * _S2919;
    float _S2922 = _S2715 * (_S2921.x + _S2921.y + _S2921.z);
    float3  _S2923 = - _S2908.differential_0 + _S2914.differential_0;
    float3  _S2924 = _S2732 * _S2923;
    float3  _S2925 = _S2710.rows[1U] * _S2923;
    float _S2926 = _S2715 * (_S2925.x + _S2925.y + _S2925.z);
    float3  _S2927 = - _S2910.differential_0 + _S2916.differential_0;
    float3  _S2928 = _S2731 * _S2927;
    float3  _S2929 = _S2710.rows[0U] * _S2927;
    float _S2930 = _S2715 * (_S2929.x + _S2929.y + _S2929.z);
    Matrix<float, 3, 3>  _S2931 = _S2844;
    _S2931[2U] = _S2920;
    _S2931[1U] = _S2924;
    _S2931[0U] = _S2928;
    Matrix<float, 3, 3>  _S2932 = transpose_0(transpose_0(_S2931));
    float _S2933 = 2.0f * - _S2932.rows[int(2)].z;
    float _S2934 = 2.0f * _S2932.rows[int(2)].y;
    float _S2935 = 2.0f * _S2932.rows[int(2)].x;
    float _S2936 = 2.0f * _S2932.rows[int(1)].z;
    float _S2937 = 2.0f * - _S2932.rows[int(1)].y;
    float _S2938 = 2.0f * _S2932.rows[int(1)].x;
    float _S2939 = 2.0f * _S2932.rows[int(0)].z;
    float _S2940 = 2.0f * _S2932.rows[int(0)].y;
    float _S2941 = 2.0f * - _S2932.rows[int(0)].x;
    float _S2942 = - _S2938 + _S2940;
    float _S2943 = _S2935 + - _S2939;
    float _S2944 = - _S2934 + _S2936;
    float _S2945 = _S2934 + _S2936;
    float _S2946 = _S2935 + _S2939;
    float _S2947 = _S2938 + _S2940;
    float _S2948 = quat_17.w * (_S2937 + _S2941);
    float _S2949 = quat_17.z * (_S2933 + _S2941);
    float _S2950 = quat_17.y * (_S2933 + _S2937);
    float _S2951 = quat_17.x * _S2942 + quat_17.z * _S2945 + quat_17.y * _S2946 + _S2948 + _S2948;
    float _S2952 = quat_17.x * _S2943 + quat_17.w * _S2945 + quat_17.y * _S2947 + _S2949 + _S2949;
    float _S2953 = quat_17.x * _S2944 + quat_17.w * _S2946 + quat_17.z * _S2947 + _S2950 + _S2950;
    float _S2954 = quat_17.w * _S2942 + quat_17.z * _S2943 + quat_17.y * _S2944;
    float3  _S2955 = _S2707;
    *&((&_S2955)->z) = _S2922;
    *&((&_S2955)->y) = _S2926;
    *&((&_S2955)->x) = _S2930;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2956;
    (&_S2956)->primal_0 = scale_16;
    (&_S2956)->differential_0 = _S2707;
    s_bwd_prop_exp_1(&_S2956, _S2955);
    Matrix<float, 2, 2>  _S2957 = _S2704;
    _S2957[int(1)] = _S2902;
    _S2957[int(0)] = _S2904;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2958;
    (&_S2958)->primal_0 = mean_c_14;
    (&_S2958)->differential_0 = _S2707;
    s_bwd_length_impl_0(&_S2958, 0.0f);
    float3  _S2959 = _S2958.differential_0 + _S2880;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2960;
    (&_S2960)->primal_0 = R_18;
    (&_S2960)->differential_0 = _S2844;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2961;
    (&_S2961)->primal_0 = mean_14;
    (&_S2961)->differential_0 = _S2707;
    s_bwd_prop_mul_1(&_S2960, &_S2961, _S2959);
    float3  _S2962 = _S2959 + _S2847.differential_0;
    Matrix<float, 3, 3>  _S2963 = _S2905.differential_0 + _S2907.differential_0 + _S2909.differential_0 + _S2911.differential_0 + _S2913.differential_0 + _S2915.differential_0 + _S2917.differential_0 + _S2960.differential_0 + _S2848;
    float3  _S2964 = _S2956.differential_0 + _S2850;
    float4  _S2965 = make_float4 (0.0f);
    *&((&_S2965)->w) = _S2951;
    *&((&_S2965)->z) = _S2952;
    *&((&_S2965)->y) = _S2953;
    *&((&_S2965)->x) = _S2954;
    float4  _S2966 = _S2965;
    float3  _S2967 = _S2906.differential_0 + _S2912.differential_0 + _S2908.differential_0 + _S2914.differential_0 + _S2910.differential_0 + _S2916.differential_0 + _S2961.differential_0 + _S2842;
    *v_mean_4 = _S2967;
    *v_quat_4 = _S2966;
    *v_scale_4 = _S2964;
    *v_in_opacity_4 = _S2887;
    (*v_sh_coeffs_4)[int(0)] = _S2864;
    (*v_sh_coeffs_4)[int(1)] = _S2865;
    (*v_sh_coeffs_4)[int(2)] = _S2866;
    (*v_sh_coeffs_4)[int(3)] = _S2867;
    (*v_sh_coeffs_4)[int(4)] = _S2868;
    (*v_sh_coeffs_4)[int(5)] = _S2869;
    (*v_sh_coeffs_4)[int(6)] = _S2870;
    (*v_sh_coeffs_4)[int(7)] = _S2871;
    (*v_sh_coeffs_4)[int(8)] = _S2872;
    (*v_sh_coeffs_4)[int(9)] = _S2873;
    (*v_sh_coeffs_4)[int(10)] = _S2874;
    (*v_sh_coeffs_4)[int(11)] = _S2875;
    (*v_sh_coeffs_4)[int(12)] = _S2876;
    (*v_sh_coeffs_4)[int(13)] = _S2877;
    (*v_sh_coeffs_4)[int(14)] = _S2878;
    (*v_sh_coeffs_4)[int(15)] = _S2879;
    *v_R_5 = _S2963;
    *v_t_5 = _S2962;
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
    float _S2968 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2968 * _S2968;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2969 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2970 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2971;
    (&_S2971)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2971)->differential_0 = _S2970;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2972;
    (&_S2972)->primal_0 = _S2969;
    (&_S2972)->differential_0 = _S2970;
    s_bwd_prop_mul_4(&_S2971, &_S2972, _s_dOut_6);
    Matrix<float, 3, 3>  _S2973 = transpose_0(transpose_0(_S2972.differential_0));
    float _S2974 = 2.0f * - _S2973.rows[int(2)].z;
    float _S2975 = 2.0f * _S2973.rows[int(2)].y;
    float _S2976 = 2.0f * _S2973.rows[int(2)].x;
    float _S2977 = 2.0f * _S2973.rows[int(1)].z;
    float _S2978 = 2.0f * - _S2973.rows[int(1)].y;
    float _S2979 = 2.0f * _S2973.rows[int(1)].x;
    float _S2980 = 2.0f * _S2973.rows[int(0)].z;
    float _S2981 = 2.0f * _S2973.rows[int(0)].y;
    float _S2982 = 2.0f * - _S2973.rows[int(0)].x;
    float _S2983 = - _S2979 + _S2981;
    float _S2984 = _S2976 + - _S2980;
    float _S2985 = - _S2975 + _S2977;
    float _S2986 = _S2975 + _S2977;
    float _S2987 = _S2976 + _S2980;
    float _S2988 = _S2979 + _S2981;
    float _S2989 = (*dpquat_0).primal_0.w * (_S2978 + _S2982);
    float _S2990 = (*dpquat_0).primal_0.z * (_S2974 + _S2982);
    float _S2991 = (*dpquat_0).primal_0.y * (_S2974 + _S2978);
    float _S2992 = (*dpquat_0).primal_0.x * _S2983 + (*dpquat_0).primal_0.z * _S2986 + (*dpquat_0).primal_0.y * _S2987 + _S2989 + _S2989;
    float _S2993 = (*dpquat_0).primal_0.x * _S2984 + (*dpquat_0).primal_0.w * _S2986 + (*dpquat_0).primal_0.y * _S2988 + _S2990 + _S2990;
    float _S2994 = (*dpquat_0).primal_0.x * _S2985 + (*dpquat_0).primal_0.w * _S2987 + (*dpquat_0).primal_0.z * _S2988 + _S2991 + _S2991;
    float _S2995 = (*dpquat_0).primal_0.w * _S2983 + (*dpquat_0).primal_0.z * _S2984 + (*dpquat_0).primal_0.y * _S2985;
    float3  _S2996 = make_float3 (_S2971.differential_0.rows[int(0)].x, _S2971.differential_0.rows[int(1)].y, _S2971.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2996;
    float4  _S2997 = make_float4 (0.0f);
    *&((&_S2997)->w) = _S2992;
    *&((&_S2997)->z) = _S2993;
    *&((&_S2997)->y) = _S2994;
    *&((&_S2997)->x) = _S2995;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2997;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2998, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2999, Matrix<float, 3, 3>  _S3000)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2998, _S2999, _S3000);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S3001 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S3001;
    float3  _S3002 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S3002;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S3003 = dOut_16.y;
    float _S3004 = dOut_16.z;
    float _S3005 = dOut_16.x;
    float _S3006 = (*a_0).primal_0.z * _S3003 + - (*a_0).primal_0.y * _S3004;
    float _S3007 = - (*a_0).primal_0.z * _S3005 + (*a_0).primal_0.x * _S3004;
    float _S3008 = (*a_0).primal_0.y * _S3005 + - (*a_0).primal_0.x * _S3003;
    float3  _S3009 = make_float3 (- (*b_0).primal_0.z * _S3003 + (*b_0).primal_0.y * _S3004, (*b_0).primal_0.z * _S3005 + - (*b_0).primal_0.x * _S3004, - (*b_0).primal_0.y * _S3005 + (*b_0).primal_0.x * _S3003);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S3009;
    float3  _S3010 = make_float3 (_S3006, _S3007, _S3008);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S3010;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S3011 = left_10.y;
    float _S3012 = right_10.z;
    float _S3013 = left_10.z;
    float _S3014 = right_10.y;
    float _S3015 = right_10.x;
    float _S3016 = left_10.x;
    return make_float3 (_S3011 * _S3012 - _S3013 * _S3014, _S3013 * _S3015 - _S3016 * _S3012, _S3016 * _S3014 - _S3011 * _S3015);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S3017, float3  _S3018)
{
    return cross_0(_S3017, _S3018);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3019, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3020, float3  _S3021)
{
    _d_cross_0(_S3019, _S3020, _S3021);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_7)
{
    float3  _S3022 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S3023 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S3022);
    float3  _S3024 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S3025 = s_primal_ctx_cross_0(_S3024, _S3023);
    float _S3026 = -0.5f * s_primal_ctx_dot_0(_S3025, _S3025);
    float _S3027 = s_primal_ctx_dot_0(_S3024, _S3024);
    float _S3028 = _S3026 / _S3027;
    float _S3029 = _S3027 * _S3027;
    float _S3030 = (*dpopacity_0).primal_0 * _s_dOut_7;
    float _S3031 = s_primal_ctx_exp_1(_S3028) * _s_dOut_7;
    DiffPair_float_0 _S3032;
    (&_S3032)->primal_0 = _S3028;
    (&_S3032)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3032, _S3030);
    float _S3033 = _S3032.differential_0 / _S3029;
    float _S3034 = _S3026 * - _S3033;
    float _S3035 = _S3027 * _S3033;
    float3  _S3036 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3037;
    (&_S3037)->primal_0 = _S3024;
    (&_S3037)->differential_0 = _S3036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3038;
    (&_S3038)->primal_0 = _S3024;
    (&_S3038)->differential_0 = _S3036;
    s_bwd_prop_dot_0(&_S3037, &_S3038, _S3034);
    float _S3039 = -0.5f * _S3035;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3040;
    (&_S3040)->primal_0 = _S3025;
    (&_S3040)->differential_0 = _S3036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3041;
    (&_S3041)->primal_0 = _S3025;
    (&_S3041)->differential_0 = _S3036;
    s_bwd_prop_dot_0(&_S3040, &_S3041, _S3039);
    float3  _S3042 = _S3041.differential_0 + _S3040.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3043;
    (&_S3043)->primal_0 = _S3024;
    (&_S3043)->differential_0 = _S3036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3044;
    (&_S3044)->primal_0 = _S3023;
    (&_S3044)->differential_0 = _S3036;
    s_bwd_prop_cross_0(&_S3043, &_S3044, _S3042);
    float3  _S3045 = _S3038.differential_0 + _S3037.differential_0 + _S3043.differential_0;
    Matrix<float, 3, 3>  _S3046 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3047;
    (&_S3047)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3047)->differential_0 = _S3046;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3048;
    (&_S3048)->primal_0 = (*dpray_d_2).primal_0;
    (&_S3048)->differential_0 = _S3036;
    s_bwd_prop_mul_1(&_S3047, &_S3048, _S3045);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3049;
    (&_S3049)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3049)->differential_0 = _S3046;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3050;
    (&_S3050)->primal_0 = _S3022;
    (&_S3050)->differential_0 = _S3036;
    s_bwd_prop_mul_1(&_S3049, &_S3050, _S3044.differential_0);
    float3  _S3051 = - _S3050.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S3048.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S3050.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S3031;
    Matrix<float, 3, 3>  _S3052 = _S3047.differential_0 + _S3049.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S3052;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S3051;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3053, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3054, DiffPair_float_0 * _S3055, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3056, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3057, float _S3058)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S3053, _S3054, _S3055, _S3056, _S3057, _S3058);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S3059 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S3059;
    Matrix<float, 3, 3>  _S3060 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S3060;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S3059;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S3059;
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
    float3  _S3061 = mean_17 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S3061, _S3061) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3062 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S3063 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S3064;
    (&_S3064)->primal_0 = s_primal_ctx_dot_0(_S3062, _S3062) + 9.99999997475242708e-07f;
    (&_S3064)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3064, _S3063);
    float3  _S3065 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3066;
    (&_S3066)->primal_0 = _S3062;
    (&_S3066)->differential_0 = _S3065;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3067;
    (&_S3067)->primal_0 = _S3062;
    (&_S3067)->differential_0 = _S3065;
    s_bwd_prop_dot_0(&_S3066, &_S3067, _S3064.differential_0);
    float3  _S3068 = _S3067.differential_0 + _S3066.differential_0;
    float3  _S3069 = - _S3068;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S3065;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S3069;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3070 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3070;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3068;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3071, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3072, DiffPair_float_0 * _S3073, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3074, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3075, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3076, float3  _S3077, float _S3078)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3071, _S3072, _S3073, _S3074, _S3075, _S3076, _S3077, _S3078);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3079 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S3079;
    Matrix<float, 3, 3>  _S3080 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3080;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3079;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S3079;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S3079;
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
    float _S3081 = _slang_select(((*dpx_14).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_14).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3081;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_15, float dOut_18)
{
    float _S3082 = (F32_exp2(((*dpx_15).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3082;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  dOut_19)
{
    float3  _S3083 = make_float3 (1.0f) / (*dpx_16).primal_0 * dOut_19;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S3083;
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

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_24, float fy_24, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S3084 = mean_c_15.z;
        bool _S3085;
        if(_S3084 < near_plane_10)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3084 > far_plane_10;
        }
        if(_S3085)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3086 = scale_19.x;
        float sx_0 = (F32_exp((_S3086)));
        float _S3087 = scale_19.y;
        float sy_0 = (F32_exp((_S3087)));
        float sz_0 = scale_19.z - 0.5f * (_S3086 + _S3087);
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
        Matrix<float, 3, 3>  _S3088 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S3088, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S3088, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S3088, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S3089 = vert0_c_0.z;
        float _S3090 = vert1_c_0.z;
        float _S3091 = vert2_c_0.z;
        if(_S3089 < near_plane_10)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3089 > far_plane_10;
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3090 < near_plane_10;
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3090 > far_plane_10;
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3091 < near_plane_10;
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = _S3091 > far_plane_10;
        }
        if(_S3085)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        for(;;)
        {
            *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3089);
            if(_S3089 < 0.0f)
            {
                _S3085 = true;
            }
            else
            {
                bool _S3092 = is_valid_distortion(*uv0_0, dist_coeffs_29);
                _S3085 = !_S3092;
            }
            if(_S3085)
            {
                break;
            }
            float u_44 = (*uv0_0).x;
            float v_44 = (*uv0_0).y;
            float r2_44 = u_44 * u_44 + v_44 * v_44;
            float2  _S3093 = *uv0_0 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_29)[int(0)] + r2_44 * ((*dist_coeffs_29)[int(1)] + r2_44 * ((*dist_coeffs_29)[int(2)] + r2_44 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_44 * v_44 + (*dist_coeffs_29)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_29)[int(6)] * r2_44, 2.0f * (*dist_coeffs_29)[int(5)] * u_44 * v_44 + (*dist_coeffs_29)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_29)[int(7)] * r2_44);
            float2  _S3094 = _S3093 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3093.x + (*dist_coeffs_29)[int(9)] * _S3093.y, 0.0f);
            *uv0_0 = make_float2 (fx_24 * _S3094.x + cx_20, fy_24 * _S3094.y + cy_20);
            break;
        }
        bool all_valid_8 = true & (!_S3085);
        for(;;)
        {
            *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3090);
            if(_S3090 < 0.0f)
            {
                _S3085 = true;
            }
            else
            {
                bool _S3095 = is_valid_distortion(*uv1_0, dist_coeffs_29);
                _S3085 = !_S3095;
            }
            if(_S3085)
            {
                break;
            }
            float u_45 = (*uv1_0).x;
            float v_45 = (*uv1_0).y;
            float r2_45 = u_45 * u_45 + v_45 * v_45;
            float2  _S3096 = *uv1_0 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_29)[int(0)] + r2_45 * ((*dist_coeffs_29)[int(1)] + r2_45 * ((*dist_coeffs_29)[int(2)] + r2_45 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_45 * v_45 + (*dist_coeffs_29)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_29)[int(6)] * r2_45, 2.0f * (*dist_coeffs_29)[int(5)] * u_45 * v_45 + (*dist_coeffs_29)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_29)[int(7)] * r2_45);
            float2  _S3097 = _S3096 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3096.x + (*dist_coeffs_29)[int(9)] * _S3096.y, 0.0f);
            *uv1_0 = make_float2 (fx_24 * _S3097.x + cx_20, fy_24 * _S3097.y + cy_20);
            break;
        }
        bool all_valid_9 = all_valid_8 & (!_S3085);
        for(;;)
        {
            *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3091);
            if(_S3091 < 0.0f)
            {
                _S3085 = true;
            }
            else
            {
                bool _S3098 = is_valid_distortion(*uv2_0, dist_coeffs_29);
                _S3085 = !_S3098;
            }
            if(_S3085)
            {
                break;
            }
            float u_46 = (*uv2_0).x;
            float v_46 = (*uv2_0).y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S3099 = *uv2_0 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_29)[int(0)] + r2_46 * ((*dist_coeffs_29)[int(1)] + r2_46 * ((*dist_coeffs_29)[int(2)] + r2_46 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_46 * v_46 + (*dist_coeffs_29)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_29)[int(6)] * r2_46, 2.0f * (*dist_coeffs_29)[int(5)] * u_46 * v_46 + (*dist_coeffs_29)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_29)[int(7)] * r2_46);
            float2  _S3100 = _S3099 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3099.x + (*dist_coeffs_29)[int(9)] * _S3099.y, 0.0f);
            *uv2_0 = make_float2 (fx_24 * _S3100.x + cx_20, fy_24 * _S3100.y + cy_20);
            break;
        }
        if(!(all_valid_9 & (!_S3085)))
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
            _S3085 = true;
        }
        else
        {
            _S3085 = xmin_5 >= float(image_width_15);
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = ymax_5 <= 0.0f;
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            _S3085 = ymin_5 >= float(image_height_15);
        }
        if(_S3085)
        {
            _S3085 = true;
        }
        else
        {
            if(_S3084 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3085 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3085 = false;
                }
                if(_S3085)
                {
                    _S3085 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3085 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3085 = false;
                    }
                }
            }
            else
            {
                _S3085 = false;
            }
        }
        if(_S3085)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3101 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S3102 = _S3101.x;
        float _S3103 = _S3101.y;
        float _S3104 = _S3101.z;
        float norm_10 = (F32_sqrt((_S3102 * _S3102 + _S3103 * _S3103 + _S3104 * _S3104)));
        float x_43 = _S3102 / norm_10;
        float y_18 = _S3103 / norm_10;
        float z_15 = _S3104 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_43 * x_43 - y_18 * y_18;
        float fS1_15 = 2.0f * x_43 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_43) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_43) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_43 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_43) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_43 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3105 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3105);
        float3  _S3106 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3107 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3106 + _S3107 + make_float3 (0.5f), _S3105);
        (*rgb_12)[int(2)] = max_0(_S3106 - _S3107 + make_float3 (0.5f), _S3105);
        float3  _S3108 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3108 * make_float3 (float(- (F32_sign((dot_0(_S3108, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_25, float fy_25, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    bool _S3109;
    bool _S3110;
    bool _S3111;
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S3112 = length_1(mean_c_16);
        bool _S3113;
        if(_S3112 < near_plane_11)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3112 > far_plane_11;
        }
        if(_S3113)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3114 = scale_20.x;
        float sx_1 = (F32_exp((_S3114)));
        float _S3115 = scale_20.y;
        float sy_1 = (F32_exp((_S3115)));
        float sz_1 = scale_20.z - 0.5f * (_S3114 + _S3115);
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
        Matrix<float, 3, 3>  _S3116 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S3116, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S3116, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S3116, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S3117 = length_1(vert0_c_1);
        float _S3118 = length_1(vert1_c_1);
        float _S3119 = length_1(vert2_c_1);
        if(_S3117 < near_plane_11)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3117 > far_plane_11;
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3118 < near_plane_11;
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3118 > far_plane_11;
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3119 < near_plane_11;
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = _S3119 > far_plane_11;
        }
        if(_S3113)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float k_7;
        for(;;)
        {
            float2  _S3120 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_24 = length_0(_S3120);
            float _S3121 = vert0_c_1.z;
            float theta_19 = (F32_atan2((r_24), (_S3121)));
            if(theta_19 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_19 * theta_19 / 3.0f) / _S3121;
            }
            else
            {
                k_7 = theta_19 / r_24;
            }
            float2  _S3122 = _S3120 * make_float2 (k_7);
            *uv0_1 = _S3122;
            bool _S3123 = is_valid_distortion(_S3122, dist_coeffs_30);
            bool _S3124 = !_S3123;
            _S3109 = _S3124;
            if(_S3124)
            {
                break;
            }
            float u_47 = (*uv0_1).x;
            float v_47 = (*uv0_1).y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S3125 = *uv0_1 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_30)[int(0)] + r2_47 * ((*dist_coeffs_30)[int(1)] + r2_47 * ((*dist_coeffs_30)[int(2)] + r2_47 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_47 * v_47 + (*dist_coeffs_30)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_30)[int(6)] * r2_47, 2.0f * (*dist_coeffs_30)[int(5)] * u_47 * v_47 + (*dist_coeffs_30)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_30)[int(7)] * r2_47);
            float2  _S3126 = _S3125 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3125.x + (*dist_coeffs_30)[int(9)] * _S3125.y, 0.0f);
            *uv0_1 = make_float2 (fx_25 * _S3126.x + cx_21, fy_25 * _S3126.y + cy_21);
            break;
        }
        bool all_valid_10 = true & (!_S3109);
        for(;;)
        {
            float2  _S3127 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_25 = length_0(_S3127);
            float _S3128 = vert1_c_1.z;
            float theta_20 = (F32_atan2((r_25), (_S3128)));
            if(theta_20 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_20 * theta_20 / 3.0f) / _S3128;
            }
            else
            {
                k_7 = theta_20 / r_25;
            }
            float2  _S3129 = _S3127 * make_float2 (k_7);
            *uv1_1 = _S3129;
            bool _S3130 = is_valid_distortion(_S3129, dist_coeffs_30);
            bool _S3131 = !_S3130;
            _S3110 = _S3131;
            if(_S3131)
            {
                break;
            }
            float u_48 = (*uv1_1).x;
            float v_48 = (*uv1_1).y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S3132 = *uv1_1 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_30)[int(0)] + r2_48 * ((*dist_coeffs_30)[int(1)] + r2_48 * ((*dist_coeffs_30)[int(2)] + r2_48 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_48 * v_48 + (*dist_coeffs_30)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_30)[int(6)] * r2_48, 2.0f * (*dist_coeffs_30)[int(5)] * u_48 * v_48 + (*dist_coeffs_30)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_30)[int(7)] * r2_48);
            float2  _S3133 = _S3132 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3132.x + (*dist_coeffs_30)[int(9)] * _S3132.y, 0.0f);
            *uv1_1 = make_float2 (fx_25 * _S3133.x + cx_21, fy_25 * _S3133.y + cy_21);
            break;
        }
        bool all_valid_11 = all_valid_10 & (!_S3110);
        for(;;)
        {
            float2  _S3134 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_26 = length_0(_S3134);
            float _S3135 = vert2_c_1.z;
            float theta_21 = (F32_atan2((r_26), (_S3135)));
            if(theta_21 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_21 * theta_21 / 3.0f) / _S3135;
            }
            else
            {
                k_7 = theta_21 / r_26;
            }
            float2  _S3136 = _S3134 * make_float2 (k_7);
            *uv2_1 = _S3136;
            bool _S3137 = is_valid_distortion(_S3136, dist_coeffs_30);
            bool _S3138 = !_S3137;
            _S3111 = _S3138;
            if(_S3138)
            {
                break;
            }
            float u_49 = (*uv2_1).x;
            float v_49 = (*uv2_1).y;
            float r2_49 = u_49 * u_49 + v_49 * v_49;
            float2  _S3139 = *uv2_1 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_30)[int(0)] + r2_49 * ((*dist_coeffs_30)[int(1)] + r2_49 * ((*dist_coeffs_30)[int(2)] + r2_49 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_49 * v_49 + (*dist_coeffs_30)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_30)[int(6)] * r2_49, 2.0f * (*dist_coeffs_30)[int(5)] * u_49 * v_49 + (*dist_coeffs_30)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_30)[int(7)] * r2_49);
            float2  _S3140 = _S3139 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3139.x + (*dist_coeffs_30)[int(9)] * _S3139.y, 0.0f);
            *uv2_1 = make_float2 (fx_25 * _S3140.x + cx_21, fy_25 * _S3140.y + cy_21);
            break;
        }
        if(!(all_valid_11 & (!_S3111)))
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
            _S3113 = true;
        }
        else
        {
            _S3113 = xmin_6 >= float(image_width_16);
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = ymax_6 <= 0.0f;
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            _S3113 = ymin_6 >= float(image_height_16);
        }
        if(_S3113)
        {
            _S3113 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3113 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3113 = false;
                }
                if(_S3113)
                {
                    _S3113 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3113 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3113 = false;
                    }
                }
            }
            else
            {
                _S3113 = false;
            }
        }
        if(_S3113)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3117, _S3118, _S3119) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3141 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S3142 = _S3141.x;
        float _S3143 = _S3141.y;
        float _S3144 = _S3141.z;
        float norm_11 = (F32_sqrt((_S3142 * _S3142 + _S3143 * _S3143 + _S3144 * _S3144)));
        float x_45 = _S3142 / norm_11;
        float y_19 = _S3143 / norm_11;
        float z_16 = _S3144 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_45 * x_45 - y_19 * y_19;
        float fS1_16 = 2.0f * x_45 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_45) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_45) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_45 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_45) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_45 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3145 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3145);
        float3  _S3146 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3147 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3146 + _S3147 + make_float3 (0.5f), _S3145);
        (*rgb_13)[int(2)] = max_0(_S3146 - _S3147 + make_float3 (0.5f), _S3145);
        float3  _S3148 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3148 * make_float3 (float(- (F32_sign((dot_0(_S3148, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_26, float fy_26, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
    float _S3149 = scale_21.x;
    float sx_2 = (F32_exp((_S3149)));
    float _S3150 = scale_21.y;
    float sy_2 = (F32_exp((_S3150)));
    float sz_2 = scale_21.z - 0.5f * (_S3149 + _S3150);
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
    Matrix<float, 3, 3>  _S3151 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S3151, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S3151, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S3151, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
    float2  _S3152 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_50 = _S3152.x;
    float v_50 = _S3152.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S3153 = 2.0f * (*dist_coeffs_31)[int(4)];
    float _S3154 = 2.0f * (*dist_coeffs_31)[int(5)];
    float2  _S3155 = _S3152 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_31)[int(0)] + r2_50 * ((*dist_coeffs_31)[int(1)] + r2_50 * ((*dist_coeffs_31)[int(2)] + r2_50 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3153 * u_50 * v_50 + (*dist_coeffs_31)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_31)[int(6)] * r2_50, _S3154 * u_50 * v_50 + (*dist_coeffs_31)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_31)[int(7)] * r2_50);
    float2  _S3156 = _S3155 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3155.x + (*dist_coeffs_31)[int(9)] * _S3155.y, 0.0f);
    *uv0_2 = make_float2 (fx_26 * _S3156.x + cx_22, fy_26 * _S3156.y + cy_22);
    float2  _S3157 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_51 = _S3157.x;
    float v_51 = _S3157.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S3158 = _S3157 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_31)[int(0)] + r2_51 * ((*dist_coeffs_31)[int(1)] + r2_51 * ((*dist_coeffs_31)[int(2)] + r2_51 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3153 * u_51 * v_51 + (*dist_coeffs_31)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_31)[int(6)] * r2_51, _S3154 * u_51 * v_51 + (*dist_coeffs_31)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_31)[int(7)] * r2_51);
    float2  _S3159 = _S3158 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3158.x + (*dist_coeffs_31)[int(9)] * _S3158.y, 0.0f);
    *uv1_2 = make_float2 (fx_26 * _S3159.x + cx_22, fy_26 * _S3159.y + cy_22);
    float2  _S3160 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_52 = _S3160.x;
    float v_52 = _S3160.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float2  _S3161 = _S3160 * make_float2 (1.0f + r2_52 * ((*dist_coeffs_31)[int(0)] + r2_52 * ((*dist_coeffs_31)[int(1)] + r2_52 * ((*dist_coeffs_31)[int(2)] + r2_52 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3153 * u_52 * v_52 + (*dist_coeffs_31)[int(5)] * (r2_52 + 2.0f * u_52 * u_52) + (*dist_coeffs_31)[int(6)] * r2_52, _S3154 * u_52 * v_52 + (*dist_coeffs_31)[int(4)] * (r2_52 + 2.0f * v_52 * v_52) + (*dist_coeffs_31)[int(7)] * r2_52);
    float2  _S3162 = _S3161 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3161.x + (*dist_coeffs_31)[int(9)] * _S3161.y, 0.0f);
    float _S3163 = fx_26 * _S3162.x + cx_22;
    float _S3164 = fy_26 * _S3162.y + cy_22;
    float2  _S3165 = make_float2 (_S3163, _S3164);
    *uv2_2 = _S3165;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S3165 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3165)));
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3163))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3164))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3163))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3164))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S3166 = mean_21 - - mul_0(transpose_0(R_21), t_20);
    float _S3167 = _S3166.x;
    float _S3168 = _S3166.y;
    float _S3169 = _S3166.z;
    float norm_12 = (F32_sqrt((_S3167 * _S3167 + _S3168 * _S3168 + _S3169 * _S3169)));
    float x_47 = _S3167 / norm_12;
    float y_20 = _S3168 / norm_12;
    float z_17 = _S3169 / norm_12;
    float z2_40 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_47 * x_47 - y_20 * y_20;
    float fS1_17 = 2.0f * x_47 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_47) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_47) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_47 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_47) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_47 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S3170 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3170);
    float3  _S3171 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S3172 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S3171 + _S3172 + make_float3 (0.5f), _S3170);
    (*rgb_14)[int(2)] = max_0(_S3171 - _S3172 + make_float3 (0.5f), _S3170);
    float3  _S3173 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S3173 * make_float3 (float(- (F32_sign((dot_0(_S3173, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_27, float fy_27, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S3174 = scale_22.x;
    float sx_3 = (F32_exp((_S3174)));
    float _S3175 = scale_22.y;
    float sy_3 = (F32_exp((_S3175)));
    float sz_3 = scale_22.z - 0.5f * (_S3174 + _S3175);
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
    Matrix<float, 3, 3>  _S3176 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S3176, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S3176, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S3176, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    float2  _S3177 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_27 = length_0(_S3177);
    float _S3178 = vert0_c_3.z;
    float theta_22 = (F32_atan2((r_27), (_S3178)));
    float k_8;
    if(theta_22 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_22 * theta_22 / 3.0f) / _S3178;
    }
    else
    {
        k_8 = theta_22 / r_27;
    }
    float2  _S3179 = _S3177 * make_float2 (k_8);
    float u_53 = _S3179.x;
    float v_53 = _S3179.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S3180 = 2.0f * (*dist_coeffs_32)[int(4)];
    float _S3181 = 2.0f * (*dist_coeffs_32)[int(5)];
    float2  _S3182 = _S3179 * make_float2 (1.0f + r2_53 * ((*dist_coeffs_32)[int(0)] + r2_53 * ((*dist_coeffs_32)[int(1)] + r2_53 * ((*dist_coeffs_32)[int(2)] + r2_53 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3180 * u_53 * v_53 + (*dist_coeffs_32)[int(5)] * (r2_53 + 2.0f * u_53 * u_53) + (*dist_coeffs_32)[int(6)] * r2_53, _S3181 * u_53 * v_53 + (*dist_coeffs_32)[int(4)] * (r2_53 + 2.0f * v_53 * v_53) + (*dist_coeffs_32)[int(7)] * r2_53);
    float2  _S3183 = _S3182 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3182.x + (*dist_coeffs_32)[int(9)] * _S3182.y, 0.0f);
    *uv0_3 = make_float2 (fx_27 * _S3183.x + cx_23, fy_27 * _S3183.y + cy_23);
    float2  _S3184 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_28 = length_0(_S3184);
    float _S3185 = vert1_c_3.z;
    float theta_23 = (F32_atan2((r_28), (_S3185)));
    if(theta_23 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_23 * theta_23 / 3.0f) / _S3185;
    }
    else
    {
        k_8 = theta_23 / r_28;
    }
    float2  _S3186 = _S3184 * make_float2 (k_8);
    float u_54 = _S3186.x;
    float v_54 = _S3186.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float2  _S3187 = _S3186 * make_float2 (1.0f + r2_54 * ((*dist_coeffs_32)[int(0)] + r2_54 * ((*dist_coeffs_32)[int(1)] + r2_54 * ((*dist_coeffs_32)[int(2)] + r2_54 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3180 * u_54 * v_54 + (*dist_coeffs_32)[int(5)] * (r2_54 + 2.0f * u_54 * u_54) + (*dist_coeffs_32)[int(6)] * r2_54, _S3181 * u_54 * v_54 + (*dist_coeffs_32)[int(4)] * (r2_54 + 2.0f * v_54 * v_54) + (*dist_coeffs_32)[int(7)] * r2_54);
    float2  _S3188 = _S3187 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3187.x + (*dist_coeffs_32)[int(9)] * _S3187.y, 0.0f);
    *uv1_3 = make_float2 (fx_27 * _S3188.x + cx_23, fy_27 * _S3188.y + cy_23);
    float2  _S3189 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_29 = length_0(_S3189);
    float _S3190 = vert2_c_3.z;
    float theta_24 = (F32_atan2((r_29), (_S3190)));
    if(theta_24 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_24 * theta_24 / 3.0f) / _S3190;
    }
    else
    {
        k_8 = theta_24 / r_29;
    }
    float2  _S3191 = _S3189 * make_float2 (k_8);
    float u_55 = _S3191.x;
    float v_55 = _S3191.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float2  _S3192 = _S3191 * make_float2 (1.0f + r2_55 * ((*dist_coeffs_32)[int(0)] + r2_55 * ((*dist_coeffs_32)[int(1)] + r2_55 * ((*dist_coeffs_32)[int(2)] + r2_55 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3180 * u_55 * v_55 + (*dist_coeffs_32)[int(5)] * (r2_55 + 2.0f * u_55 * u_55) + (*dist_coeffs_32)[int(6)] * r2_55, _S3181 * u_55 * v_55 + (*dist_coeffs_32)[int(4)] * (r2_55 + 2.0f * v_55 * v_55) + (*dist_coeffs_32)[int(7)] * r2_55);
    float2  _S3193 = _S3192 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3192.x + (*dist_coeffs_32)[int(9)] * _S3192.y, 0.0f);
    float _S3194 = fx_27 * _S3193.x + cx_23;
    float _S3195 = fy_27 * _S3193.y + cy_23;
    float2  _S3196 = make_float2 (_S3194, _S3195);
    *uv2_3 = _S3196;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3196 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3196)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3194))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3195))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3194))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3195))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3197 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S3198 = _S3197.x;
    float _S3199 = _S3197.y;
    float _S3200 = _S3197.z;
    float norm_13 = (F32_sqrt((_S3198 * _S3198 + _S3199 * _S3199 + _S3200 * _S3200)));
    float x_49 = _S3198 / norm_13;
    float y_21 = _S3199 / norm_13;
    float z_18 = _S3200 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_49 * x_49 - y_21 * y_21;
    float fS1_18 = 2.0f * x_49 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_49) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_49) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_49 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_49) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_49 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3201 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3201);
    float3  _S3202 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3203 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3202 + _S3203 + make_float3 (0.5f), _S3201);
    (*rgb_15)[int(2)] = max_0(_S3202 - _S3203 + make_float3 (0.5f), _S3201);
    float3  _S3204 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3204 * make_float3 (float(- (F32_sign((dot_0(_S3204, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3205, float3  _S3206)
{
    _d_log_vector_0(_S3205, _S3206);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3207, float _S3208)
{
    _d_exp2_0(_S3207, _S3208);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3209, float _S3210)
{
    _d_abs_0(_S3209, _S3210);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_28, float fy_28, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float _S3211 = scale_23.x;
    float _S3212 = s_primal_ctx_exp_1(_S3211);
    float _S3213 = scale_23.y;
    float _S3214 = s_primal_ctx_exp_1(_S3213);
    float sz_4 = scale_23.z - 0.5f * (_S3211 + _S3213);
    float _S3215 = quat_24.y;
    float x2_24 = _S3215 * _S3215;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S3216 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S3217 = make_float3 (_S3212, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S3216, _S3217) + mean_23;
    float _S3218 = -0.5f + sz_4;
    float3  _S3219 = make_float3 (_S3212 * _S3218, _S3214, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S3216, _S3219) + mean_23;
    float _S3220 = -0.5f - sz_4;
    float3  _S3221 = make_float3 (_S3212 * _S3220, - _S3214, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S3216, _S3221) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S3222 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3223 = vert0_c_4.z;
    float2  _S3224 = make_float2 (_S3223);
    float2  _S3225 = _S3222 / make_float2 (_S3223);
    float2  _S3226 = make_float2 (_S3223 * _S3223);
    float u_56 = _S3225.x;
    float v_56 = _S3225.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S3227 = (*dist_coeffs_33)[int(2)] + r2_56 * (*dist_coeffs_33)[int(3)];
    float _S3228 = (*dist_coeffs_33)[int(1)] + r2_56 * _S3227;
    float _S3229 = (*dist_coeffs_33)[int(0)] + r2_56 * _S3228;
    float radial_2 = 1.0f + r2_56 * _S3229;
    float2  _S3230 = make_float2 (radial_2);
    float _S3231 = 2.0f * (*dist_coeffs_33)[int(4)];
    float _S3232 = _S3231 * u_56;
    float _S3233 = 2.0f * u_56;
    float _S3234 = 2.0f * (*dist_coeffs_33)[int(5)];
    float _S3235 = _S3234 * u_56;
    float _S3236 = 2.0f * v_56;
    float2  _S3237 = _S3225 * make_float2 (radial_2) + make_float2 (_S3232 * v_56 + (*dist_coeffs_33)[int(5)] * (r2_56 + _S3233 * u_56) + (*dist_coeffs_33)[int(6)] * r2_56, _S3235 * v_56 + (*dist_coeffs_33)[int(4)] * (r2_56 + _S3236 * v_56) + (*dist_coeffs_33)[int(7)] * r2_56);
    float2  _S3238 = _S3237 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3237.x + (*dist_coeffs_33)[int(9)] * _S3237.y, 0.0f);
    float _S3239 = fx_28 * _S3238.x + cx_24;
    float _S3240 = fy_28 * _S3238.y + cy_24;
    float2  _S3241 = make_float2 (_S3239, _S3240);
    float2  _S3242 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3243 = vert1_c_4.z;
    float2  _S3244 = make_float2 (_S3243);
    float2  _S3245 = _S3242 / make_float2 (_S3243);
    float2  _S3246 = make_float2 (_S3243 * _S3243);
    float u_57 = _S3245.x;
    float v_57 = _S3245.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S3247 = (*dist_coeffs_33)[int(2)] + r2_57 * (*dist_coeffs_33)[int(3)];
    float _S3248 = (*dist_coeffs_33)[int(1)] + r2_57 * _S3247;
    float _S3249 = (*dist_coeffs_33)[int(0)] + r2_57 * _S3248;
    float radial_3 = 1.0f + r2_57 * _S3249;
    float2  _S3250 = make_float2 (radial_3);
    float _S3251 = _S3231 * u_57;
    float _S3252 = 2.0f * u_57;
    float _S3253 = _S3234 * u_57;
    float _S3254 = 2.0f * v_57;
    float2  _S3255 = _S3245 * make_float2 (radial_3) + make_float2 (_S3251 * v_57 + (*dist_coeffs_33)[int(5)] * (r2_57 + _S3252 * u_57) + (*dist_coeffs_33)[int(6)] * r2_57, _S3253 * v_57 + (*dist_coeffs_33)[int(4)] * (r2_57 + _S3254 * v_57) + (*dist_coeffs_33)[int(7)] * r2_57);
    float2  _S3256 = _S3255 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3255.x + (*dist_coeffs_33)[int(9)] * _S3255.y, 0.0f);
    float _S3257 = fx_28 * _S3256.x + cx_24;
    float _S3258 = fy_28 * _S3256.y + cy_24;
    float2  _S3259 = make_float2 (_S3257, _S3258);
    float2  _S3260 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3261 = vert2_c_4.z;
    float2  _S3262 = make_float2 (_S3261);
    float2  _S3263 = _S3260 / make_float2 (_S3261);
    float2  _S3264 = make_float2 (_S3261 * _S3261);
    float u_58 = _S3263.x;
    float v_58 = _S3263.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S3265 = (*dist_coeffs_33)[int(2)] + r2_58 * (*dist_coeffs_33)[int(3)];
    float _S3266 = (*dist_coeffs_33)[int(1)] + r2_58 * _S3265;
    float _S3267 = (*dist_coeffs_33)[int(0)] + r2_58 * _S3266;
    float radial_4 = 1.0f + r2_58 * _S3267;
    float2  _S3268 = make_float2 (radial_4);
    float _S3269 = _S3231 * u_58;
    float _S3270 = 2.0f * u_58;
    float _S3271 = _S3234 * u_58;
    float _S3272 = 2.0f * v_58;
    float2  _S3273 = _S3263 * make_float2 (radial_4) + make_float2 (_S3269 * v_58 + (*dist_coeffs_33)[int(5)] * (r2_58 + _S3270 * u_58) + (*dist_coeffs_33)[int(6)] * r2_58, _S3271 * v_58 + (*dist_coeffs_33)[int(4)] * (r2_58 + _S3272 * v_58) + (*dist_coeffs_33)[int(7)] * r2_58);
    float2  _S3274 = _S3273 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3273.x + (*dist_coeffs_33)[int(9)] * _S3273.y, 0.0f);
    float _S3275 = fx_28 * _S3274.x + cx_24;
    float _S3276 = fy_28 * _S3274.y + cy_24;
    float2  _S3277 = make_float2 (_S3275, _S3276);
    float2  e0_4 = _S3259 - _S3241;
    float2  e1_4 = _S3277 - _S3259;
    float2  e2_0 = _S3241 - _S3277;
    float _S3278 = e0_4.x;
    float _S3279 = e1_4.y;
    float _S3280 = e0_4.y;
    float _S3281 = e1_4.x;
    float _S3282 = _S3278 * _S3279 - _S3280 * _S3281;
    float _S3283 = 1.0f - hardness_4.y;
    float _S3284 = -1.0f / _S3283;
    float _S3285 = _S3283 * _S3283;
    float _S3286 = s_primal_ctx_max_0(_S3239, _S3257);
    float _S3287 = s_primal_ctx_min_0(_S3239, _S3257);
    float _S3288 = s_primal_ctx_max_0(_S3240, _S3258);
    float _S3289 = s_primal_ctx_min_0(_S3240, _S3258);
    float3  _S3290 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3291 = transpose_0(R_23);
    float3  _S3292 = mean_23 - - s_primal_ctx_mul_1(_S3291, t_22);
    float _S3293 = _S3292.x;
    float _S3294 = _S3292.y;
    float _S3295 = _S3292.z;
    float _S3296 = _S3293 * _S3293 + _S3294 * _S3294 + _S3295 * _S3295;
    float _S3297 = s_primal_ctx_sqrt_0(_S3296);
    float x_50 = _S3293 / _S3297;
    float3  _S3298 = make_float3 (x_50);
    float _S3299 = _S3297 * _S3297;
    float y_22 = _S3294 / _S3297;
    float z_19 = _S3295 / _S3297;
    float3  _S3300 = make_float3 (z_19);
    float _S3301 = - y_22;
    float3  _S3302 = make_float3 (_S3301);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_50 * x_50 - y_22 * y_22;
    float _S3303 = 2.0f * x_50;
    float fS1_19 = _S3303 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S3304 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_50;
    float3  _S3305 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3306 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3307 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3308 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3309 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3309;
    float3  _S3310 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_50;
    float3  _S3311 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3312 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3313 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3314 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_50 * fC1_19 - y_22 * fS1_19);
    float3  _S3315 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_50 * fS1_19 + y_22 * fC1_19);
    float3  _S3316 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3301) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_50) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3317 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3318 = make_float3 (0.0f);
    float3  _S3319 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3320 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3321 = make_float3 (_S3320);
    float3  _S3322 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3320);
    float3  _S3323 = _S3319 + _S3322 + make_float3 (0.5f);
    float3  _S3324 = _S3319 - _S3322 + make_float3 (0.5f);
    float3  _S3325 = vert1_c_4 - vert0_c_4;
    float3  _S3326 = vert2_c_4 - vert0_c_4;
    float3  _S3327 = s_primal_ctx_cross_0(_S3325, _S3326);
    float3  _S3328 = normalize_0(_S3327);
    float3  _S3329 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3328, mean_c_19)))))) * v_normal_0;
    float3  _S3330 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3331;
    (&_S3331)->primal_0 = _S3328;
    (&_S3331)->differential_0 = _S3330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3332;
    (&_S3332)->primal_0 = mean_c_19;
    (&_S3332)->differential_0 = _S3330;
    s_bwd_prop_dot_0(&_S3331, &_S3332, 0.0f);
    float3  _S3333 = _S3329 + _S3331.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3334;
    (&_S3334)->primal_0 = _S3327;
    (&_S3334)->differential_0 = _S3330;
    s_bwd_normalize_impl_0(&_S3334, _S3333);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3335;
    (&_S3335)->primal_0 = _S3325;
    (&_S3335)->differential_0 = _S3330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3336;
    (&_S3336)->primal_0 = _S3326;
    (&_S3336)->differential_0 = _S3330;
    s_bwd_prop_cross_0(&_S3335, &_S3336, _S3334.differential_0);
    float3  _S3337 = - _S3336.differential_0;
    float3  _S3338 = - _S3335.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3339;
    (&_S3339)->primal_0 = _S3324;
    (&_S3339)->differential_0 = _S3330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3340;
    (&_S3340)->primal_0 = _S3318;
    (&_S3340)->differential_0 = _S3330;
    s_bwd_prop_max_0(&_S3339, &_S3340, (*v_rgb_6)[int(2)]);
    float3  _S3341 = - _S3339.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3342;
    (&_S3342)->primal_0 = _S3323;
    (&_S3342)->differential_0 = _S3330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3343;
    (&_S3343)->primal_0 = _S3318;
    (&_S3343)->differential_0 = _S3330;
    s_bwd_prop_max_0(&_S3342, &_S3343, (*v_rgb_6)[int(1)]);
    float3  _S3344 = _S3321 * (_S3341 + _S3342.differential_0);
    float3  _S3345 = _S3339.differential_0 + _S3342.differential_0;
    float3  _S3346 = make_float3 (0.5f) * - _S3345;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3347;
    (&_S3347)->primal_0 = _S3317;
    (&_S3347)->differential_0 = _S3330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3348;
    (&_S3348)->primal_0 = _S3318;
    (&_S3348)->differential_0 = _S3330;
    s_bwd_prop_max_0(&_S3347, &_S3348, (*v_rgb_6)[int(0)]);
    float3  _S3349 = _S3346 + _S3347.differential_0;
    float3  _S3350 = _S3345 + _S3347.differential_0;
    float3  _S3351 = _S3315 * _S3350;
    float3  _S3352 = (*sh_coeffs_19)[int(15)] * _S3350;
    float3  _S3353 = _S3313 * _S3350;
    float3  _S3354 = (*sh_coeffs_19)[int(14)] * _S3350;
    float3  _S3355 = _S3311 * _S3350;
    float3  _S3356 = (*sh_coeffs_19)[int(13)] * _S3350;
    float3  _S3357 = _S3310 * _S3350;
    float3  _S3358 = (*sh_coeffs_19)[int(12)] * _S3350;
    float3  _S3359 = _S3312 * _S3350;
    float3  _S3360 = (*sh_coeffs_19)[int(11)] * _S3350;
    float3  _S3361 = _S3314 * _S3350;
    float3  _S3362 = (*sh_coeffs_19)[int(10)] * _S3350;
    float3  _S3363 = _S3316 * _S3350;
    float3  _S3364 = (*sh_coeffs_19)[int(9)] * _S3350;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3364.x + _S3364.y + _S3364.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3352.x + _S3352.y + _S3352.z);
    float _S3365 = _S3362.x + _S3362.y + _S3362.z;
    float _S3366 = _S3354.x + _S3354.y + _S3354.z;
    float _S3367 = _S3360.x + _S3360.y + _S3360.z;
    float _S3368 = _S3356.x + _S3356.y + _S3356.z;
    float _S3369 = _S3358.x + _S3358.y + _S3358.z;
    float _S3370 = - s_diff_fC2_T_5;
    float3  _S3371 = _S3307 * _S3350;
    float3  _S3372 = (*sh_coeffs_19)[int(8)] * _S3350;
    float3  _S3373 = _S3305 * _S3350;
    float3  _S3374 = (*sh_coeffs_19)[int(7)] * _S3350;
    float3  _S3375 = _S3304 * _S3350;
    float3  _S3376 = (*sh_coeffs_19)[int(6)] * _S3350;
    float3  _S3377 = _S3306 * _S3350;
    float3  _S3378 = (*sh_coeffs_19)[int(5)] * _S3350;
    float3  _S3379 = _S3308 * _S3350;
    float3  _S3380 = (*sh_coeffs_19)[int(4)] * _S3350;
    float _S3381 = _S3378.x + _S3378.y + _S3378.z;
    float _S3382 = _S3374.x + _S3374.y + _S3374.z;
    float _S3383 = fTmp1B_19 * _S3365 + x_50 * s_diff_fS2_T_5 + y_22 * _S3370 + 0.54627424478530884f * (_S3380.x + _S3380.y + _S3380.z);
    float _S3384 = fTmp1B_19 * _S3366 + y_22 * s_diff_fS2_T_5 + x_50 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3372.x + _S3372.y + _S3372.z);
    float _S3385 = y_22 * - _S3384;
    float _S3386 = x_50 * _S3384;
    float _S3387 = z_19 * (1.86588168144226074f * (z_19 * _S3369) + -2.28522896766662598f * (y_22 * _S3367 + x_50 * _S3368) + 0.94617468118667603f * (_S3376.x + _S3376.y + _S3376.z));
    float3  _S3388 = make_float3 (0.48860251903533936f) * _S3350;
    float3  _S3389 = - _S3388;
    float3  _S3390 = _S3298 * _S3389;
    float3  _S3391 = (*sh_coeffs_19)[int(3)] * _S3389;
    float3  _S3392 = _S3300 * _S3388;
    float3  _S3393 = (*sh_coeffs_19)[int(2)] * _S3388;
    float3  _S3394 = _S3302 * _S3388;
    float3  _S3395 = (*sh_coeffs_19)[int(1)] * _S3388;
    float _S3396 = (_S3309 * _S3369 + 1.44530570507049561f * (fS1_19 * _S3365 + fC1_19 * _S3366) + -1.09254848957061768f * (y_22 * _S3381 + x_50 * _S3382) + _S3387 + _S3387 + _S3393.x + _S3393.y + _S3393.z) / _S3299;
    float _S3397 = _S3297 * _S3396;
    float _S3398 = (fTmp0C_19 * _S3367 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3370 + fTmp0B_19 * _S3381 + _S3303 * _S3383 + _S3385 + _S3385 + - (_S3395.x + _S3395.y + _S3395.z)) / _S3299;
    float _S3399 = _S3297 * _S3398;
    float _S3400 = (fTmp0C_19 * _S3368 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3382 + 2.0f * (y_22 * _S3383) + _S3386 + _S3386 + _S3391.x + _S3391.y + _S3391.z) / _S3299;
    float _S3401 = _S3297 * _S3400;
    float _S3402 = _S3295 * - _S3396 + _S3294 * - _S3398 + _S3293 * - _S3400;
    DiffPair_float_0 _S3403;
    (&_S3403)->primal_0 = _S3296;
    (&_S3403)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3403, _S3402);
    float _S3404 = _S3295 * _S3403.differential_0;
    float _S3405 = _S3294 * _S3403.differential_0;
    float _S3406 = _S3293 * _S3403.differential_0;
    float3  _S3407 = make_float3 (0.282094806432724f) * _S3350;
    float3  _S3408 = make_float3 (_S3401 + _S3406 + _S3406, _S3399 + _S3405 + _S3405, _S3397 + _S3404 + _S3404);
    float3  _S3409 = - - _S3408;
    Matrix<float, 3, 3>  _S3410 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3411;
    (&_S3411)->primal_0 = _S3291;
    (&_S3411)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3412;
    (&_S3412)->primal_0 = t_22;
    (&_S3412)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3411, &_S3412, _S3409);
    Matrix<float, 3, 3>  _S3413 = transpose_0(_S3411.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3414;
    (&_S3414)->primal_0 = _S3290;
    (&_S3414)->differential_0 = _S3330;
    s_bwd_prop_log_1(&_S3414, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3415;
    (&_S3415)->primal_0 = vert2_c_4;
    (&_S3415)->differential_0 = _S3330;
    s_bwd_length_impl_0(&_S3415, _S3414.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3416;
    (&_S3416)->primal_0 = vert1_c_4;
    (&_S3416)->differential_0 = _S3330;
    s_bwd_length_impl_0(&_S3416, _S3414.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3417;
    (&_S3417)->primal_0 = vert0_c_4;
    (&_S3417)->differential_0 = _S3330;
    s_bwd_length_impl_0(&_S3417, _S3414.differential_0.x);
    DiffPair_float_0 _S3418;
    (&_S3418)->primal_0 = _S3289;
    (&_S3418)->differential_0 = 0.0f;
    DiffPair_float_0 _S3419;
    (&_S3419)->primal_0 = _S3276;
    (&_S3419)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3418, &_S3419, 0.0f);
    DiffPair_float_0 _S3420;
    (&_S3420)->primal_0 = _S3240;
    (&_S3420)->differential_0 = 0.0f;
    DiffPair_float_0 _S3421;
    (&_S3421)->primal_0 = _S3258;
    (&_S3421)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3420, &_S3421, _S3418.differential_0);
    DiffPair_float_0 _S3422;
    (&_S3422)->primal_0 = _S3288;
    (&_S3422)->differential_0 = 0.0f;
    DiffPair_float_0 _S3423;
    (&_S3423)->primal_0 = _S3276;
    (&_S3423)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3422, &_S3423, 0.0f);
    DiffPair_float_0 _S3424;
    (&_S3424)->primal_0 = _S3240;
    (&_S3424)->differential_0 = 0.0f;
    DiffPair_float_0 _S3425;
    (&_S3425)->primal_0 = _S3258;
    (&_S3425)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3424, &_S3425, _S3422.differential_0);
    DiffPair_float_0 _S3426;
    (&_S3426)->primal_0 = _S3287;
    (&_S3426)->differential_0 = 0.0f;
    DiffPair_float_0 _S3427;
    (&_S3427)->primal_0 = _S3275;
    (&_S3427)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3426, &_S3427, 0.0f);
    DiffPair_float_0 _S3428;
    (&_S3428)->primal_0 = _S3239;
    (&_S3428)->differential_0 = 0.0f;
    DiffPair_float_0 _S3429;
    (&_S3429)->primal_0 = _S3257;
    (&_S3429)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3428, &_S3429, _S3426.differential_0);
    DiffPair_float_0 _S3430;
    (&_S3430)->primal_0 = _S3286;
    (&_S3430)->differential_0 = 0.0f;
    DiffPair_float_0 _S3431;
    (&_S3431)->primal_0 = _S3275;
    (&_S3431)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3430, &_S3431, 0.0f);
    DiffPair_float_0 _S3432;
    (&_S3432)->primal_0 = _S3239;
    (&_S3432)->differential_0 = 0.0f;
    DiffPair_float_0 _S3433;
    (&_S3433)->primal_0 = _S3257;
    (&_S3433)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3432, &_S3433, _S3430.differential_0);
    DiffPair_float_0 _S3434;
    (&_S3434)->primal_0 = _S3284;
    (&_S3434)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3434, 0.0f);
    float _S3435 = - (-1.0f * - (_S3434.differential_0 / _S3285));
    float2  _S3436 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3437;
    (&_S3437)->primal_0 = e2_0;
    (&_S3437)->differential_0 = _S3436;
    s_bwd_length_impl_1(&_S3437, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3438;
    (&_S3438)->primal_0 = e1_4;
    (&_S3438)->differential_0 = _S3436;
    s_bwd_length_impl_1(&_S3438, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3439;
    (&_S3439)->primal_0 = e0_4;
    (&_S3439)->differential_0 = _S3436;
    s_bwd_length_impl_1(&_S3439, -0.0f);
    DiffPair_float_0 _S3440;
    (&_S3440)->primal_0 = _S3282;
    (&_S3440)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3440, 0.0f);
    float _S3441 = - _S3440.differential_0;
    float2  _S3442 = _S3438.differential_0 + make_float2 (_S3280 * _S3441, _S3278 * _S3440.differential_0);
    float2  _S3443 = _S3439.differential_0 + make_float2 (_S3279 * _S3440.differential_0, _S3281 * _S3441);
    float2  _S3444 = v_uv2_0 + - _S3437.differential_0 + _S3442;
    float _S3445 = fx_28 * (_S3427.differential_0 + _S3431.differential_0 + _S3444.x);
    float2  _S3446 = make_float2 (_S3445, fy_28 * (_S3419.differential_0 + _S3423.differential_0 + _S3444.y)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3445, (*dist_coeffs_33)[int(9)] * _S3445);
    float2  _S3447 = _S3263 * _S3446;
    float _S3448 = (*dist_coeffs_33)[int(4)] * _S3446.y;
    float _S3449 = (*dist_coeffs_33)[int(5)] * _S3446.x;
    float _S3450 = _S3447.x + _S3447.y;
    float _S3451 = r2_58 * _S3450;
    float _S3452 = r2_58 * _S3451;
    float _S3453 = (*dist_coeffs_33)[int(7)] * _S3446.y + _S3448 + (*dist_coeffs_33)[int(6)] * _S3446.x + _S3449 + _S3267 * _S3450 + _S3266 * _S3451 + _S3265 * _S3452 + (*dist_coeffs_33)[int(3)] * (r2_58 * _S3452);
    float _S3454 = v_58 * _S3453;
    float _S3455 = u_58 * _S3453;
    float2  _S3456 = (_S3268 * _S3446 + make_float2 (_S3234 * (v_58 * _S3446.y) + _S3270 * _S3449 + 2.0f * (u_58 * _S3449) + _S3231 * (v_58 * _S3446.x) + _S3455 + _S3455, _S3272 * _S3448 + 2.0f * (v_58 * _S3448) + _S3271 * _S3446.y + _S3269 * _S3446.x + _S3454 + _S3454)) / _S3264;
    float2  _S3457 = _S3260 * - _S3456;
    float2  _S3458 = _S3262 * _S3456;
    float2  _S3459 = v_uv1_0 + - _S3442 + _S3443;
    float _S3460 = fx_28 * (_S3429.differential_0 + _S3433.differential_0 + _S3459.x);
    float2  _S3461 = make_float2 (_S3460, fy_28 * (_S3421.differential_0 + _S3425.differential_0 + _S3459.y)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3460, (*dist_coeffs_33)[int(9)] * _S3460);
    float2  _S3462 = _S3245 * _S3461;
    float _S3463 = (*dist_coeffs_33)[int(4)] * _S3461.y;
    float _S3464 = (*dist_coeffs_33)[int(5)] * _S3461.x;
    float _S3465 = _S3462.x + _S3462.y;
    float _S3466 = r2_57 * _S3465;
    float _S3467 = r2_57 * _S3466;
    float _S3468 = (*dist_coeffs_33)[int(7)] * _S3461.y + _S3463 + (*dist_coeffs_33)[int(6)] * _S3461.x + _S3464 + _S3249 * _S3465 + _S3248 * _S3466 + _S3247 * _S3467 + (*dist_coeffs_33)[int(3)] * (r2_57 * _S3467);
    float _S3469 = v_57 * _S3468;
    float _S3470 = u_57 * _S3468;
    float2  _S3471 = (_S3250 * _S3461 + make_float2 (_S3234 * (v_57 * _S3461.y) + _S3252 * _S3464 + 2.0f * (u_57 * _S3464) + _S3231 * (v_57 * _S3461.x) + _S3470 + _S3470, _S3254 * _S3463 + 2.0f * (v_57 * _S3463) + _S3253 * _S3461.y + _S3251 * _S3461.x + _S3469 + _S3469)) / _S3246;
    float2  _S3472 = _S3242 * - _S3471;
    float2  _S3473 = _S3244 * _S3471;
    float _S3474 = _S3472.x + _S3472.y;
    float2  _S3475 = v_uv0_0 + _S3437.differential_0 + - _S3443;
    float _S3476 = fx_28 * (_S3428.differential_0 + _S3432.differential_0 + _S3475.x);
    float2  _S3477 = make_float2 (_S3476, fy_28 * (_S3420.differential_0 + _S3424.differential_0 + _S3475.y)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3476, (*dist_coeffs_33)[int(9)] * _S3476);
    float2  _S3478 = _S3225 * _S3477;
    float _S3479 = (*dist_coeffs_33)[int(4)] * _S3477.y;
    float _S3480 = (*dist_coeffs_33)[int(5)] * _S3477.x;
    float _S3481 = _S3478.x + _S3478.y;
    float _S3482 = r2_56 * _S3481;
    float _S3483 = r2_56 * _S3482;
    float _S3484 = (*dist_coeffs_33)[int(7)] * _S3477.y + _S3479 + (*dist_coeffs_33)[int(6)] * _S3477.x + _S3480 + _S3229 * _S3481 + _S3228 * _S3482 + _S3227 * _S3483 + (*dist_coeffs_33)[int(3)] * (r2_56 * _S3483);
    float _S3485 = v_56 * _S3484;
    float _S3486 = u_56 * _S3484;
    float2  _S3487 = (_S3230 * _S3477 + make_float2 (_S3234 * (v_56 * _S3477.y) + _S3233 * _S3480 + 2.0f * (u_56 * _S3480) + _S3231 * (v_56 * _S3477.x) + _S3486 + _S3486, _S3236 * _S3479 + 2.0f * (v_56 * _S3479) + _S3235 * _S3477.y + _S3232 * _S3477.x + _S3485 + _S3485)) / _S3226;
    float2  _S3488 = _S3222 * - _S3487;
    float2  _S3489 = _S3224 * _S3487;
    float _S3490 = _S3488.x + _S3488.y;
    float3  _S3491 = _S3336.differential_0 + _S3415.differential_0 + make_float3 (_S3458.x, _S3458.y, _S3457.x + _S3457.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3492;
    (&_S3492)->primal_0 = R_23;
    (&_S3492)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3493;
    (&_S3493)->primal_0 = vert2_0;
    (&_S3493)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3492, &_S3493, _S3491);
    float3  _S3494 = _S3335.differential_0 + _S3416.differential_0 + make_float3 (_S3473.x, _S3473.y, _S3474);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3495;
    (&_S3495)->primal_0 = R_23;
    (&_S3495)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3496;
    (&_S3496)->primal_0 = vert1_0;
    (&_S3496)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3495, &_S3496, _S3494);
    float3  _S3497 = _S3337 + _S3338 + _S3417.differential_0 + make_float3 (_S3489.x, _S3489.y, _S3490);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3498;
    (&_S3498)->primal_0 = R_23;
    (&_S3498)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3499;
    (&_S3499)->primal_0 = vert0_0;
    (&_S3499)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3498, &_S3499, _S3497);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3500;
    (&_S3500)->primal_0 = _S3216;
    (&_S3500)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3501;
    (&_S3501)->primal_0 = _S3221;
    (&_S3501)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3500, &_S3501, _S3493.differential_0);
    float _S3502 = - _S3501.differential_0.y;
    float _S3503 = _S3220 * _S3501.differential_0.x;
    float _S3504 = - (_S3212 * _S3501.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3505;
    (&_S3505)->primal_0 = _S3216;
    (&_S3505)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3506;
    (&_S3506)->primal_0 = _S3219;
    (&_S3506)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3505, &_S3506, _S3496.differential_0);
    float _S3507 = _S3212 * _S3506.differential_0.x;
    float _S3508 = _S3218 * _S3506.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3509;
    (&_S3509)->primal_0 = _S3216;
    (&_S3509)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3510;
    (&_S3510)->primal_0 = _S3217;
    (&_S3510)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3509, &_S3510, _S3499.differential_0);
    Matrix<float, 3, 3>  _S3511 = transpose_0(_S3500.differential_0 + _S3505.differential_0 + _S3509.differential_0);
    float _S3512 = 2.0f * - _S3511.rows[int(2)].z;
    float _S3513 = 2.0f * _S3511.rows[int(2)].y;
    float _S3514 = 2.0f * _S3511.rows[int(2)].x;
    float _S3515 = 2.0f * _S3511.rows[int(1)].z;
    float _S3516 = 2.0f * - _S3511.rows[int(1)].y;
    float _S3517 = 2.0f * _S3511.rows[int(1)].x;
    float _S3518 = 2.0f * _S3511.rows[int(0)].z;
    float _S3519 = 2.0f * _S3511.rows[int(0)].y;
    float _S3520 = 2.0f * - _S3511.rows[int(0)].x;
    float _S3521 = - _S3517 + _S3519;
    float _S3522 = _S3514 + - _S3518;
    float _S3523 = - _S3513 + _S3515;
    float _S3524 = _S3513 + _S3515;
    float _S3525 = _S3514 + _S3518;
    float _S3526 = _S3517 + _S3519;
    float _S3527 = quat_24.w * (_S3516 + _S3520);
    float _S3528 = quat_24.z * (_S3512 + _S3520);
    float _S3529 = quat_24.y * (_S3512 + _S3516);
    float _S3530 = quat_24.x * _S3521 + quat_24.z * _S3524 + quat_24.y * _S3525 + _S3527 + _S3527;
    float _S3531 = quat_24.x * _S3522 + quat_24.w * _S3524 + quat_24.y * _S3526 + _S3528 + _S3528;
    float _S3532 = quat_24.x * _S3523 + quat_24.w * _S3525 + quat_24.z * _S3526 + _S3529 + _S3529;
    float _S3533 = quat_24.w * _S3521 + quat_24.z * _S3522 + quat_24.y * _S3523;
    float _S3534 = _S3504 + _S3507;
    float _S3535 = 0.5f * - _S3534;
    float _S3536 = _S3502 + _S3506.differential_0.y;
    DiffPair_float_0 _S3537;
    (&_S3537)->primal_0 = _S3213;
    (&_S3537)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3537, _S3536);
    float _S3538 = _S3535 + _S3537.differential_0;
    float _S3539 = _S3503 + _S3508 + _S3510.differential_0.x;
    DiffPair_float_0 _S3540;
    (&_S3540)->primal_0 = _S3211;
    (&_S3540)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3540, _S3539);
    float _S3541 = _S3535 + _S3540.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3542;
    (&_S3542)->primal_0 = R_23;
    (&_S3542)->differential_0 = _S3410;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3543;
    (&_S3543)->primal_0 = mean_23;
    (&_S3543)->differential_0 = _S3330;
    s_bwd_prop_mul_1(&_S3542, &_S3543, _S3332.differential_0);
    float3  _S3544 = _S3412.differential_0 + _S3491 + _S3494 + _S3497 + _S3332.differential_0;
    Matrix<float, 3, 3>  _S3545 = _S3413 + _S3492.differential_0 + _S3495.differential_0 + _S3498.differential_0 + _S3542.differential_0;
    FixedArray<float3 , 2>  _S3546;
    _S3546[int(0)] = _S3330;
    _S3546[int(1)] = _S3330;
    _S3546[int(1)] = _S3344;
    _S3546[int(0)] = _S3349;
    FixedArray<float3 , 16>  _S3547;
    _S3547[int(0)] = _S3330;
    _S3547[int(1)] = _S3330;
    _S3547[int(2)] = _S3330;
    _S3547[int(3)] = _S3330;
    _S3547[int(4)] = _S3330;
    _S3547[int(5)] = _S3330;
    _S3547[int(6)] = _S3330;
    _S3547[int(7)] = _S3330;
    _S3547[int(8)] = _S3330;
    _S3547[int(9)] = _S3330;
    _S3547[int(10)] = _S3330;
    _S3547[int(11)] = _S3330;
    _S3547[int(12)] = _S3330;
    _S3547[int(13)] = _S3330;
    _S3547[int(14)] = _S3330;
    _S3547[int(15)] = _S3330;
    _S3547[int(15)] = _S3351;
    _S3547[int(14)] = _S3353;
    _S3547[int(13)] = _S3355;
    _S3547[int(12)] = _S3357;
    _S3547[int(11)] = _S3359;
    _S3547[int(10)] = _S3361;
    _S3547[int(9)] = _S3363;
    _S3547[int(8)] = _S3371;
    _S3547[int(7)] = _S3373;
    _S3547[int(6)] = _S3375;
    _S3547[int(5)] = _S3377;
    _S3547[int(4)] = _S3379;
    _S3547[int(3)] = _S3390;
    _S3547[int(2)] = _S3392;
    _S3547[int(1)] = _S3394;
    _S3547[int(0)] = _S3407;
    float2  _S3548 = v_out_hardness_0 + make_float2 (0.0f, _S3435);
    float3  _S3549 = make_float3 (_S3541, _S3538, _S3534);
    float4  _S3550 = make_float4 (0.0f);
    *&((&_S3550)->w) = _S3530;
    *&((&_S3550)->z) = _S3531;
    *&((&_S3550)->y) = _S3532;
    *&((&_S3550)->x) = _S3533;
    *v_mean_7 = _S3408 + _S3493.differential_0 + _S3496.differential_0 + _S3499.differential_0 + _S3543.differential_0;
    *v_quat_6 = _S3550;
    *v_scale_6 = _S3549;
    *v_hardness_0 = _S3548;
    *v_sh_coeffs_5 = _S3547;
    *v_ch_coeffs_0 = _S3546;
    *v_R_6 = _S3545;
    *v_t_6 = _S3544;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_29, float fy_29, float cx_25, float cy_25, FixedArray<float, 10>  * dist_coeffs_34, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S3551 = scale_24.x;
    float _S3552 = s_primal_ctx_exp_1(_S3551);
    float _S3553 = scale_24.y;
    float _S3554 = s_primal_ctx_exp_1(_S3553);
    float sz_5 = scale_24.z - 0.5f * (_S3551 + _S3553);
    float _S3555 = quat_25.y;
    float x2_25 = _S3555 * _S3555;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3556 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3557 = make_float3 (_S3552, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S3556, _S3557) + mean_24;
    float _S3558 = -0.5f + sz_5;
    float3  _S3559 = make_float3 (_S3552 * _S3558, _S3554, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S3556, _S3559) + mean_24;
    float _S3560 = -0.5f - sz_5;
    float3  _S3561 = make_float3 (_S3552 * _S3560, - _S3554, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S3556, _S3561) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    float _S3562 = length_1(vert0_c_5);
    float _S3563 = length_1(vert1_c_5);
    float _S3564 = length_1(vert2_c_5);
    float2  _S3565 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3566 = length_0(_S3565);
    float _S3567 = vert0_c_5.z;
    float _S3568 = s_primal_ctx_atan2_0(_S3566, _S3567);
    bool _S3569 = _S3568 < 0.00100000004749745f;
    float k_9;
    float _S3570;
    float _S3571;
    float _S3572;
    if(_S3569)
    {
        float _S3573 = 1.0f - _S3568 * _S3568 / 3.0f;
        float _S3574 = _S3567 * _S3567;
        k_9 = _S3573 / _S3567;
        _S3570 = 0.0f;
        _S3571 = _S3574;
        _S3572 = _S3573;
    }
    else
    {
        float _S3575 = _S3566 * _S3566;
        k_9 = _S3568 / _S3566;
        _S3570 = _S3575;
        _S3571 = 0.0f;
        _S3572 = 0.0f;
    }
    float2  _S3576 = make_float2 (k_9);
    float2  _S3577 = _S3565 * make_float2 (k_9);
    float u_59 = _S3577.x;
    float v_59 = _S3577.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S3578 = (*dist_coeffs_34)[int(2)] + r2_59 * (*dist_coeffs_34)[int(3)];
    float _S3579 = (*dist_coeffs_34)[int(1)] + r2_59 * _S3578;
    float _S3580 = (*dist_coeffs_34)[int(0)] + r2_59 * _S3579;
    float radial_5 = 1.0f + r2_59 * _S3580;
    float2  _S3581 = make_float2 (radial_5);
    float _S3582 = 2.0f * (*dist_coeffs_34)[int(4)];
    float _S3583 = _S3582 * u_59;
    float _S3584 = 2.0f * u_59;
    float _S3585 = 2.0f * (*dist_coeffs_34)[int(5)];
    float _S3586 = _S3585 * u_59;
    float _S3587 = 2.0f * v_59;
    float2  _S3588 = _S3577 * make_float2 (radial_5) + make_float2 (_S3583 * v_59 + (*dist_coeffs_34)[int(5)] * (r2_59 + _S3584 * u_59) + (*dist_coeffs_34)[int(6)] * r2_59, _S3586 * v_59 + (*dist_coeffs_34)[int(4)] * (r2_59 + _S3587 * v_59) + (*dist_coeffs_34)[int(7)] * r2_59);
    float2  _S3589 = _S3588 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3588.x + (*dist_coeffs_34)[int(9)] * _S3588.y, 0.0f);
    float _S3590 = fx_29 * _S3589.x + cx_25;
    float _S3591 = fy_29 * _S3589.y + cy_25;
    float2  _S3592 = make_float2 (_S3590, _S3591);
    float2  _S3593 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3594 = length_0(_S3593);
    float _S3595 = vert1_c_5.z;
    float _S3596 = s_primal_ctx_atan2_0(_S3594, _S3595);
    bool _S3597 = _S3596 < 0.00100000004749745f;
    float _S3598;
    float _S3599;
    float _S3600;
    if(_S3597)
    {
        float _S3601 = 1.0f - _S3596 * _S3596 / 3.0f;
        float _S3602 = _S3595 * _S3595;
        k_9 = _S3601 / _S3595;
        _S3598 = 0.0f;
        _S3599 = _S3602;
        _S3600 = _S3601;
    }
    else
    {
        float _S3603 = _S3594 * _S3594;
        k_9 = _S3596 / _S3594;
        _S3598 = _S3603;
        _S3599 = 0.0f;
        _S3600 = 0.0f;
    }
    float2  _S3604 = make_float2 (k_9);
    float2  _S3605 = _S3593 * make_float2 (k_9);
    float u_60 = _S3605.x;
    float v_60 = _S3605.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S3606 = (*dist_coeffs_34)[int(2)] + r2_60 * (*dist_coeffs_34)[int(3)];
    float _S3607 = (*dist_coeffs_34)[int(1)] + r2_60 * _S3606;
    float _S3608 = (*dist_coeffs_34)[int(0)] + r2_60 * _S3607;
    float radial_6 = 1.0f + r2_60 * _S3608;
    float2  _S3609 = make_float2 (radial_6);
    float _S3610 = _S3582 * u_60;
    float _S3611 = 2.0f * u_60;
    float _S3612 = _S3585 * u_60;
    float _S3613 = 2.0f * v_60;
    float2  _S3614 = _S3605 * make_float2 (radial_6) + make_float2 (_S3610 * v_60 + (*dist_coeffs_34)[int(5)] * (r2_60 + _S3611 * u_60) + (*dist_coeffs_34)[int(6)] * r2_60, _S3612 * v_60 + (*dist_coeffs_34)[int(4)] * (r2_60 + _S3613 * v_60) + (*dist_coeffs_34)[int(7)] * r2_60);
    float2  _S3615 = _S3614 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3614.x + (*dist_coeffs_34)[int(9)] * _S3614.y, 0.0f);
    float _S3616 = fx_29 * _S3615.x + cx_25;
    float _S3617 = fy_29 * _S3615.y + cy_25;
    float2  _S3618 = make_float2 (_S3616, _S3617);
    float2  _S3619 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3620 = length_0(_S3619);
    float _S3621 = vert2_c_5.z;
    float _S3622 = s_primal_ctx_atan2_0(_S3620, _S3621);
    bool _S3623 = _S3622 < 0.00100000004749745f;
    float _S3624;
    float _S3625;
    float _S3626;
    if(_S3623)
    {
        float _S3627 = 1.0f - _S3622 * _S3622 / 3.0f;
        float _S3628 = _S3621 * _S3621;
        k_9 = _S3627 / _S3621;
        _S3624 = 0.0f;
        _S3625 = _S3628;
        _S3626 = _S3627;
    }
    else
    {
        float _S3629 = _S3620 * _S3620;
        k_9 = _S3622 / _S3620;
        _S3624 = _S3629;
        _S3625 = 0.0f;
        _S3626 = 0.0f;
    }
    float2  _S3630 = make_float2 (k_9);
    float2  _S3631 = _S3619 * make_float2 (k_9);
    float u_61 = _S3631.x;
    float v_61 = _S3631.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S3632 = (*dist_coeffs_34)[int(2)] + r2_61 * (*dist_coeffs_34)[int(3)];
    float _S3633 = (*dist_coeffs_34)[int(1)] + r2_61 * _S3632;
    float _S3634 = (*dist_coeffs_34)[int(0)] + r2_61 * _S3633;
    float radial_7 = 1.0f + r2_61 * _S3634;
    float2  _S3635 = make_float2 (radial_7);
    float _S3636 = _S3582 * u_61;
    float _S3637 = 2.0f * u_61;
    float _S3638 = _S3585 * u_61;
    float _S3639 = 2.0f * v_61;
    float2  _S3640 = _S3631 * make_float2 (radial_7) + make_float2 (_S3636 * v_61 + (*dist_coeffs_34)[int(5)] * (r2_61 + _S3637 * u_61) + (*dist_coeffs_34)[int(6)] * r2_61, _S3638 * v_61 + (*dist_coeffs_34)[int(4)] * (r2_61 + _S3639 * v_61) + (*dist_coeffs_34)[int(7)] * r2_61);
    float2  _S3641 = _S3640 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3640.x + (*dist_coeffs_34)[int(9)] * _S3640.y, 0.0f);
    float _S3642 = fx_29 * _S3641.x + cx_25;
    float _S3643 = fy_29 * _S3641.y + cy_25;
    float2  _S3644 = make_float2 (_S3642, _S3643);
    float2  e0_5 = _S3618 - _S3592;
    float2  e1_5 = _S3644 - _S3618;
    float2  e2_1 = _S3592 - _S3644;
    float _S3645 = e0_5.x;
    float _S3646 = e1_5.y;
    float _S3647 = e0_5.y;
    float _S3648 = e1_5.x;
    float _S3649 = _S3645 * _S3646 - _S3647 * _S3648;
    float _S3650 = 1.0f - hardness_5.y;
    float _S3651 = -1.0f / _S3650;
    float _S3652 = _S3650 * _S3650;
    float _S3653 = s_primal_ctx_max_0(_S3590, _S3616);
    float _S3654 = s_primal_ctx_min_0(_S3590, _S3616);
    float _S3655 = s_primal_ctx_max_0(_S3591, _S3617);
    float _S3656 = s_primal_ctx_min_0(_S3591, _S3617);
    float3  _S3657 = make_float3 (_S3562, _S3563, _S3564) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3658 = transpose_0(R_24);
    float3  _S3659 = mean_24 - - s_primal_ctx_mul_1(_S3658, t_23);
    float _S3660 = _S3659.x;
    float _S3661 = _S3659.y;
    float _S3662 = _S3659.z;
    float _S3663 = _S3660 * _S3660 + _S3661 * _S3661 + _S3662 * _S3662;
    float _S3664 = s_primal_ctx_sqrt_0(_S3663);
    float x_51 = _S3660 / _S3664;
    float3  _S3665 = make_float3 (x_51);
    float _S3666 = _S3664 * _S3664;
    float y_23 = _S3661 / _S3664;
    float z_20 = _S3662 / _S3664;
    float3  _S3667 = make_float3 (z_20);
    float _S3668 = - y_23;
    float3  _S3669 = make_float3 (_S3668);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_51 * x_51 - y_23 * y_23;
    float _S3670 = 2.0f * x_51;
    float fS1_20 = _S3670 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3671 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_51;
    float3  _S3672 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3673 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3674 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3675 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3676 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3676;
    float3  _S3677 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_51;
    float3  _S3678 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3679 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3680 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3681 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_51 * fC1_20 - y_23 * fS1_20);
    float3  _S3682 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_51 * fS1_20 + y_23 * fC1_20);
    float3  _S3683 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3668) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_51) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3684 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3685 = make_float3 (0.0f);
    float3  _S3686 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3687 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3688 = make_float3 (_S3687);
    float3  _S3689 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3687);
    float3  _S3690 = _S3686 + _S3689 + make_float3 (0.5f);
    float3  _S3691 = _S3686 - _S3689 + make_float3 (0.5f);
    float3  _S3692 = vert1_c_5 - vert0_c_5;
    float3  _S3693 = vert2_c_5 - vert0_c_5;
    float3  _S3694 = s_primal_ctx_cross_0(_S3692, _S3693);
    float3  _S3695 = normalize_0(_S3694);
    float3  _S3696 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3695, mean_c_20)))))) * v_normal_1;
    float3  _S3697 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3698;
    (&_S3698)->primal_0 = _S3695;
    (&_S3698)->differential_0 = _S3697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3699;
    (&_S3699)->primal_0 = mean_c_20;
    (&_S3699)->differential_0 = _S3697;
    s_bwd_prop_dot_0(&_S3698, &_S3699, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3700 = _S3699;
    float3  _S3701 = _S3696 + _S3698.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3702;
    (&_S3702)->primal_0 = _S3694;
    (&_S3702)->differential_0 = _S3697;
    s_bwd_normalize_impl_0(&_S3702, _S3701);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3703;
    (&_S3703)->primal_0 = _S3692;
    (&_S3703)->differential_0 = _S3697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3704;
    (&_S3704)->primal_0 = _S3693;
    (&_S3704)->differential_0 = _S3697;
    s_bwd_prop_cross_0(&_S3703, &_S3704, _S3702.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3705 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3706 = _S3704;
    float3  _S3707 = - _S3704.differential_0;
    float3  _S3708 = - _S3703.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3709;
    (&_S3709)->primal_0 = _S3691;
    (&_S3709)->differential_0 = _S3697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3710;
    (&_S3710)->primal_0 = _S3685;
    (&_S3710)->differential_0 = _S3697;
    s_bwd_prop_max_0(&_S3709, &_S3710, (*v_rgb_7)[int(2)]);
    float3  _S3711 = - _S3709.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3712;
    (&_S3712)->primal_0 = _S3690;
    (&_S3712)->differential_0 = _S3697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3713;
    (&_S3713)->primal_0 = _S3685;
    (&_S3713)->differential_0 = _S3697;
    s_bwd_prop_max_0(&_S3712, &_S3713, (*v_rgb_7)[int(1)]);
    float3  _S3714 = _S3688 * (_S3711 + _S3712.differential_0);
    float3  _S3715 = _S3709.differential_0 + _S3712.differential_0;
    float3  _S3716 = make_float3 (0.5f) * - _S3715;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3717;
    (&_S3717)->primal_0 = _S3684;
    (&_S3717)->differential_0 = _S3697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3718;
    (&_S3718)->primal_0 = _S3685;
    (&_S3718)->differential_0 = _S3697;
    s_bwd_prop_max_0(&_S3717, &_S3718, (*v_rgb_7)[int(0)]);
    float3  _S3719 = _S3716 + _S3717.differential_0;
    float3  _S3720 = _S3715 + _S3717.differential_0;
    float3  _S3721 = _S3682 * _S3720;
    float3  _S3722 = (*sh_coeffs_20)[int(15)] * _S3720;
    float3  _S3723 = _S3680 * _S3720;
    float3  _S3724 = (*sh_coeffs_20)[int(14)] * _S3720;
    float3  _S3725 = _S3678 * _S3720;
    float3  _S3726 = (*sh_coeffs_20)[int(13)] * _S3720;
    float3  _S3727 = _S3677 * _S3720;
    float3  _S3728 = (*sh_coeffs_20)[int(12)] * _S3720;
    float3  _S3729 = _S3679 * _S3720;
    float3  _S3730 = (*sh_coeffs_20)[int(11)] * _S3720;
    float3  _S3731 = _S3681 * _S3720;
    float3  _S3732 = (*sh_coeffs_20)[int(10)] * _S3720;
    float3  _S3733 = _S3683 * _S3720;
    float3  _S3734 = (*sh_coeffs_20)[int(9)] * _S3720;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3734.x + _S3734.y + _S3734.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3722.x + _S3722.y + _S3722.z);
    float _S3735 = _S3732.x + _S3732.y + _S3732.z;
    float _S3736 = _S3724.x + _S3724.y + _S3724.z;
    float _S3737 = _S3730.x + _S3730.y + _S3730.z;
    float _S3738 = _S3726.x + _S3726.y + _S3726.z;
    float _S3739 = _S3728.x + _S3728.y + _S3728.z;
    float _S3740 = - s_diff_fC2_T_6;
    float3  _S3741 = _S3674 * _S3720;
    float3  _S3742 = (*sh_coeffs_20)[int(8)] * _S3720;
    float3  _S3743 = _S3672 * _S3720;
    float3  _S3744 = (*sh_coeffs_20)[int(7)] * _S3720;
    float3  _S3745 = _S3671 * _S3720;
    float3  _S3746 = (*sh_coeffs_20)[int(6)] * _S3720;
    float3  _S3747 = _S3673 * _S3720;
    float3  _S3748 = (*sh_coeffs_20)[int(5)] * _S3720;
    float3  _S3749 = _S3675 * _S3720;
    float3  _S3750 = (*sh_coeffs_20)[int(4)] * _S3720;
    float _S3751 = _S3748.x + _S3748.y + _S3748.z;
    float _S3752 = _S3744.x + _S3744.y + _S3744.z;
    float _S3753 = fTmp1B_20 * _S3735 + x_51 * s_diff_fS2_T_6 + y_23 * _S3740 + 0.54627424478530884f * (_S3750.x + _S3750.y + _S3750.z);
    float _S3754 = fTmp1B_20 * _S3736 + y_23 * s_diff_fS2_T_6 + x_51 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3742.x + _S3742.y + _S3742.z);
    float _S3755 = y_23 * - _S3754;
    float _S3756 = x_51 * _S3754;
    float _S3757 = z_20 * (1.86588168144226074f * (z_20 * _S3739) + -2.28522896766662598f * (y_23 * _S3737 + x_51 * _S3738) + 0.94617468118667603f * (_S3746.x + _S3746.y + _S3746.z));
    float3  _S3758 = make_float3 (0.48860251903533936f) * _S3720;
    float3  _S3759 = - _S3758;
    float3  _S3760 = _S3665 * _S3759;
    float3  _S3761 = (*sh_coeffs_20)[int(3)] * _S3759;
    float3  _S3762 = _S3667 * _S3758;
    float3  _S3763 = (*sh_coeffs_20)[int(2)] * _S3758;
    float3  _S3764 = _S3669 * _S3758;
    float3  _S3765 = (*sh_coeffs_20)[int(1)] * _S3758;
    float _S3766 = (_S3676 * _S3739 + 1.44530570507049561f * (fS1_20 * _S3735 + fC1_20 * _S3736) + -1.09254848957061768f * (y_23 * _S3751 + x_51 * _S3752) + _S3757 + _S3757 + _S3763.x + _S3763.y + _S3763.z) / _S3666;
    float _S3767 = _S3664 * _S3766;
    float _S3768 = (fTmp0C_20 * _S3737 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3740 + fTmp0B_20 * _S3751 + _S3670 * _S3753 + _S3755 + _S3755 + - (_S3765.x + _S3765.y + _S3765.z)) / _S3666;
    float _S3769 = _S3664 * _S3768;
    float _S3770 = (fTmp0C_20 * _S3738 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3752 + 2.0f * (y_23 * _S3753) + _S3756 + _S3756 + _S3761.x + _S3761.y + _S3761.z) / _S3666;
    float _S3771 = _S3664 * _S3770;
    float _S3772 = _S3662 * - _S3766 + _S3661 * - _S3768 + _S3660 * - _S3770;
    DiffPair_float_0 _S3773;
    (&_S3773)->primal_0 = _S3663;
    (&_S3773)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3773, _S3772);
    float _S3774 = _S3662 * _S3773.differential_0;
    float _S3775 = _S3661 * _S3773.differential_0;
    float _S3776 = _S3660 * _S3773.differential_0;
    float3  _S3777 = make_float3 (0.282094806432724f) * _S3720;
    float3  _S3778 = make_float3 (_S3771 + _S3776 + _S3776, _S3769 + _S3775 + _S3775, _S3767 + _S3774 + _S3774);
    float3  _S3779 = - - _S3778;
    Matrix<float, 3, 3>  _S3780 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3781;
    (&_S3781)->primal_0 = _S3658;
    (&_S3781)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3782;
    (&_S3782)->primal_0 = t_23;
    (&_S3782)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3781, &_S3782, _S3779);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3783 = _S3782;
    Matrix<float, 3, 3>  _S3784 = transpose_0(_S3781.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3785;
    (&_S3785)->primal_0 = _S3657;
    (&_S3785)->differential_0 = _S3697;
    s_bwd_prop_log_1(&_S3785, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3786 = _S3785;
    DiffPair_float_0 _S3787;
    (&_S3787)->primal_0 = _S3656;
    (&_S3787)->differential_0 = 0.0f;
    DiffPair_float_0 _S3788;
    (&_S3788)->primal_0 = _S3643;
    (&_S3788)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3787, &_S3788, 0.0f);
    DiffPair_float_0 _S3789;
    (&_S3789)->primal_0 = _S3591;
    (&_S3789)->differential_0 = 0.0f;
    DiffPair_float_0 _S3790;
    (&_S3790)->primal_0 = _S3617;
    (&_S3790)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3789, &_S3790, _S3787.differential_0);
    DiffPair_float_0 _S3791;
    (&_S3791)->primal_0 = _S3655;
    (&_S3791)->differential_0 = 0.0f;
    DiffPair_float_0 _S3792;
    (&_S3792)->primal_0 = _S3643;
    (&_S3792)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3791, &_S3792, 0.0f);
    DiffPair_float_0 _S3793;
    (&_S3793)->primal_0 = _S3591;
    (&_S3793)->differential_0 = 0.0f;
    DiffPair_float_0 _S3794;
    (&_S3794)->primal_0 = _S3617;
    (&_S3794)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3793, &_S3794, _S3791.differential_0);
    DiffPair_float_0 _S3795;
    (&_S3795)->primal_0 = _S3654;
    (&_S3795)->differential_0 = 0.0f;
    DiffPair_float_0 _S3796;
    (&_S3796)->primal_0 = _S3642;
    (&_S3796)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3795, &_S3796, 0.0f);
    DiffPair_float_0 _S3797;
    (&_S3797)->primal_0 = _S3590;
    (&_S3797)->differential_0 = 0.0f;
    DiffPair_float_0 _S3798;
    (&_S3798)->primal_0 = _S3616;
    (&_S3798)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3797, &_S3798, _S3795.differential_0);
    DiffPair_float_0 _S3799;
    (&_S3799)->primal_0 = _S3653;
    (&_S3799)->differential_0 = 0.0f;
    DiffPair_float_0 _S3800;
    (&_S3800)->primal_0 = _S3642;
    (&_S3800)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3799, &_S3800, 0.0f);
    DiffPair_float_0 _S3801;
    (&_S3801)->primal_0 = _S3590;
    (&_S3801)->differential_0 = 0.0f;
    DiffPair_float_0 _S3802;
    (&_S3802)->primal_0 = _S3616;
    (&_S3802)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3801, &_S3802, _S3799.differential_0);
    DiffPair_float_0 _S3803;
    (&_S3803)->primal_0 = _S3651;
    (&_S3803)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3803, 0.0f);
    float _S3804 = - (-1.0f * - (_S3803.differential_0 / _S3652));
    float2  _S3805 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3806;
    (&_S3806)->primal_0 = e2_1;
    (&_S3806)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3806, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3807;
    (&_S3807)->primal_0 = e1_5;
    (&_S3807)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3807, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3808;
    (&_S3808)->primal_0 = e0_5;
    (&_S3808)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3808, -0.0f);
    DiffPair_float_0 _S3809;
    (&_S3809)->primal_0 = _S3649;
    (&_S3809)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3809, 0.0f);
    float _S3810 = - _S3809.differential_0;
    float2  _S3811 = _S3807.differential_0 + make_float2 (_S3647 * _S3810, _S3645 * _S3809.differential_0);
    float2  _S3812 = _S3808.differential_0 + make_float2 (_S3646 * _S3809.differential_0, _S3648 * _S3810);
    float2  _S3813 = v_uv2_1 + - _S3806.differential_0 + _S3811;
    float _S3814 = fx_29 * (_S3796.differential_0 + _S3800.differential_0 + _S3813.x);
    float2  _S3815 = make_float2 (_S3814, fy_29 * (_S3788.differential_0 + _S3792.differential_0 + _S3813.y)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3814, (*dist_coeffs_34)[int(9)] * _S3814);
    float2  _S3816 = _S3631 * _S3815;
    float2  _S3817 = _S3635 * _S3815;
    float _S3818 = (*dist_coeffs_34)[int(4)] * _S3815.y;
    float _S3819 = (*dist_coeffs_34)[int(5)] * _S3815.x;
    float _S3820 = _S3816.x + _S3816.y;
    float _S3821 = r2_61 * _S3820;
    float _S3822 = r2_61 * _S3821;
    float _S3823 = (*dist_coeffs_34)[int(7)] * _S3815.y + _S3818 + (*dist_coeffs_34)[int(6)] * _S3815.x + _S3819 + _S3634 * _S3820 + _S3633 * _S3821 + _S3632 * _S3822 + (*dist_coeffs_34)[int(3)] * (r2_61 * _S3822);
    float _S3824 = v_61 * _S3823;
    float _S3825 = u_61 * _S3823;
    float _S3826 = _S3639 * _S3818 + 2.0f * (v_61 * _S3818) + _S3638 * _S3815.y + _S3636 * _S3815.x + _S3824 + _S3824;
    float _S3827 = _S3585 * (v_61 * _S3815.y) + _S3637 * _S3819 + 2.0f * (u_61 * _S3819) + _S3582 * (v_61 * _S3815.x) + _S3825 + _S3825;
    float2  _S3828 = v_uv0_1 + _S3806.differential_0 + - _S3812;
    float2  _S3829 = v_out_hardness_1 + make_float2 (0.0f, _S3804);
    float _S3830 = _S3798.differential_0 + _S3802.differential_0;
    float2  _S3831 = v_uv1_1 + - _S3811 + _S3812;
    float3  _S3832 = _S3707 + _S3708;
    FixedArray<float3 , 2>  _S3833;
    _S3833[int(0)] = _S3697;
    _S3833[int(1)] = _S3697;
    _S3833[int(1)] = _S3714;
    _S3833[int(0)] = _S3719;
    float3  _S3834 = _S3833[int(0)];
    float3  _S3835 = _S3833[int(1)];
    FixedArray<float3 , 16>  _S3836;
    _S3836[int(0)] = _S3697;
    _S3836[int(1)] = _S3697;
    _S3836[int(2)] = _S3697;
    _S3836[int(3)] = _S3697;
    _S3836[int(4)] = _S3697;
    _S3836[int(5)] = _S3697;
    _S3836[int(6)] = _S3697;
    _S3836[int(7)] = _S3697;
    _S3836[int(8)] = _S3697;
    _S3836[int(9)] = _S3697;
    _S3836[int(10)] = _S3697;
    _S3836[int(11)] = _S3697;
    _S3836[int(12)] = _S3697;
    _S3836[int(13)] = _S3697;
    _S3836[int(14)] = _S3697;
    _S3836[int(15)] = _S3697;
    _S3836[int(7)] = _S3743;
    _S3836[int(0)] = _S3777;
    _S3836[int(1)] = _S3764;
    _S3836[int(2)] = _S3762;
    _S3836[int(3)] = _S3760;
    _S3836[int(4)] = _S3749;
    _S3836[int(5)] = _S3747;
    _S3836[int(6)] = _S3745;
    _S3836[int(15)] = _S3721;
    _S3836[int(8)] = _S3741;
    _S3836[int(9)] = _S3733;
    _S3836[int(10)] = _S3731;
    _S3836[int(11)] = _S3729;
    _S3836[int(12)] = _S3727;
    _S3836[int(13)] = _S3725;
    _S3836[int(14)] = _S3723;
    float3  _S3837 = _S3836[int(0)];
    float3  _S3838 = _S3836[int(1)];
    float3  _S3839 = _S3836[int(2)];
    float3  _S3840 = _S3836[int(3)];
    float3  _S3841 = _S3836[int(4)];
    float3  _S3842 = _S3836[int(5)];
    float3  _S3843 = _S3836[int(6)];
    float3  _S3844 = _S3836[int(7)];
    float3  _S3845 = _S3836[int(8)];
    float3  _S3846 = _S3836[int(9)];
    float3  _S3847 = _S3836[int(10)];
    float3  _S3848 = _S3836[int(11)];
    float3  _S3849 = _S3836[int(12)];
    float3  _S3850 = _S3836[int(13)];
    float3  _S3851 = _S3836[int(14)];
    float3  _S3852 = _S3836[int(15)];
    float _S3853 = _S3797.differential_0 + _S3801.differential_0;
    float _S3854 = _S3789.differential_0 + _S3793.differential_0;
    float _S3855 = _S3790.differential_0 + _S3794.differential_0;
    float2  _S3856 = _S3817 + make_float2 (_S3827, _S3826);
    float2  _S3857 = _S3619 * _S3856;
    float2  _S3858 = _S3630 * _S3856;
    float _S3859 = _S3857.x + _S3857.y;
    if(_S3623)
    {
        float _S3860 = _S3859 / _S3625;
        float _S3861 = _S3626 * - _S3860;
        float _S3862 = _S3622 * (0.3333333432674408f * - (_S3621 * _S3860));
        k_9 = _S3862 + _S3862;
        _S3624 = _S3861;
        _S3625 = 0.0f;
    }
    else
    {
        float _S3863 = _S3859 / _S3624;
        float _S3864 = _S3622 * - _S3863;
        k_9 = _S3620 * _S3863;
        _S3624 = 0.0f;
        _S3625 = _S3864;
    }
    DiffPair_float_0 _S3865;
    (&_S3865)->primal_0 = _S3620;
    (&_S3865)->differential_0 = 0.0f;
    DiffPair_float_0 _S3866;
    (&_S3866)->primal_0 = _S3621;
    (&_S3866)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3865, &_S3866, k_9);
    float _S3867 = _S3866.differential_0 + _S3624;
    float _S3868 = _S3865.differential_0 + _S3625;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3869;
    (&_S3869)->primal_0 = _S3619;
    (&_S3869)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3869, _S3868);
    float2  _S3870 = _S3869.differential_0 + _S3858;
    float _S3871 = fx_29 * (_S3831.x + _S3830);
    float2  _S3872 = make_float2 (_S3871, fy_29 * (_S3831.y + _S3855)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3871, (*dist_coeffs_34)[int(9)] * _S3871);
    float2  _S3873 = _S3605 * _S3872;
    float _S3874 = (*dist_coeffs_34)[int(4)] * _S3872.y;
    float _S3875 = (*dist_coeffs_34)[int(5)] * _S3872.x;
    float _S3876 = _S3873.x + _S3873.y;
    float _S3877 = r2_60 * _S3876;
    float _S3878 = r2_60 * _S3877;
    float _S3879 = (*dist_coeffs_34)[int(7)] * _S3872.y + _S3874 + (*dist_coeffs_34)[int(6)] * _S3872.x + _S3875 + _S3608 * _S3876 + _S3607 * _S3877 + _S3606 * _S3878 + (*dist_coeffs_34)[int(3)] * (r2_60 * _S3878);
    float _S3880 = v_60 * _S3879;
    float _S3881 = u_60 * _S3879;
    float3  _S3882 = _S3706.differential_0 + make_float3 (_S3870.x, _S3870.y, _S3867);
    float2  _S3883 = _S3609 * _S3872 + make_float2 (_S3585 * (v_60 * _S3872.y) + _S3611 * _S3875 + 2.0f * (u_60 * _S3875) + _S3582 * (v_60 * _S3872.x) + _S3881 + _S3881, _S3613 * _S3874 + 2.0f * (v_60 * _S3874) + _S3612 * _S3872.y + _S3610 * _S3872.x + _S3880 + _S3880);
    float2  _S3884 = _S3593 * _S3883;
    float2  _S3885 = _S3604 * _S3883;
    float _S3886 = _S3884.x + _S3884.y;
    if(_S3597)
    {
        float _S3887 = _S3886 / _S3599;
        float _S3888 = _S3600 * - _S3887;
        float _S3889 = _S3596 * (0.3333333432674408f * - (_S3595 * _S3887));
        k_9 = _S3889 + _S3889;
        _S3598 = _S3888;
        _S3599 = 0.0f;
    }
    else
    {
        float _S3890 = _S3886 / _S3598;
        float _S3891 = _S3596 * - _S3890;
        k_9 = _S3594 * _S3890;
        _S3598 = 0.0f;
        _S3599 = _S3891;
    }
    DiffPair_float_0 _S3892;
    (&_S3892)->primal_0 = _S3594;
    (&_S3892)->differential_0 = 0.0f;
    DiffPair_float_0 _S3893;
    (&_S3893)->primal_0 = _S3595;
    (&_S3893)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3892, &_S3893, k_9);
    float _S3894 = _S3893.differential_0 + _S3598;
    float _S3895 = _S3892.differential_0 + _S3599;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3896;
    (&_S3896)->primal_0 = _S3593;
    (&_S3896)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3896, _S3895);
    float2  _S3897 = _S3896.differential_0 + _S3885;
    float _S3898 = fx_29 * (_S3828.x + _S3853);
    float2  _S3899 = make_float2 (_S3898, fy_29 * (_S3828.y + _S3854)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3898, (*dist_coeffs_34)[int(9)] * _S3898);
    float2  _S3900 = _S3577 * _S3899;
    float _S3901 = (*dist_coeffs_34)[int(4)] * _S3899.y;
    float _S3902 = (*dist_coeffs_34)[int(5)] * _S3899.x;
    float _S3903 = _S3900.x + _S3900.y;
    float _S3904 = r2_59 * _S3903;
    float _S3905 = r2_59 * _S3904;
    float _S3906 = (*dist_coeffs_34)[int(7)] * _S3899.y + _S3901 + (*dist_coeffs_34)[int(6)] * _S3899.x + _S3902 + _S3580 * _S3903 + _S3579 * _S3904 + _S3578 * _S3905 + (*dist_coeffs_34)[int(3)] * (r2_59 * _S3905);
    float _S3907 = v_59 * _S3906;
    float _S3908 = u_59 * _S3906;
    float3  _S3909 = _S3705.differential_0 + make_float3 (_S3897.x, _S3897.y, _S3894);
    float2  _S3910 = _S3581 * _S3899 + make_float2 (_S3585 * (v_59 * _S3899.y) + _S3584 * _S3902 + 2.0f * (u_59 * _S3902) + _S3582 * (v_59 * _S3899.x) + _S3908 + _S3908, _S3587 * _S3901 + 2.0f * (v_59 * _S3901) + _S3586 * _S3899.y + _S3583 * _S3899.x + _S3907 + _S3907);
    float2  _S3911 = _S3565 * _S3910;
    float2  _S3912 = _S3576 * _S3910;
    float _S3913 = _S3911.x + _S3911.y;
    if(_S3569)
    {
        float _S3914 = _S3913 / _S3571;
        float _S3915 = _S3572 * - _S3914;
        float _S3916 = _S3568 * (0.3333333432674408f * - (_S3567 * _S3914));
        k_9 = _S3916 + _S3916;
        _S3570 = _S3915;
        _S3571 = 0.0f;
    }
    else
    {
        float _S3917 = _S3913 / _S3570;
        float _S3918 = _S3568 * - _S3917;
        k_9 = _S3566 * _S3917;
        _S3570 = 0.0f;
        _S3571 = _S3918;
    }
    DiffPair_float_0 _S3919;
    (&_S3919)->primal_0 = _S3566;
    (&_S3919)->differential_0 = 0.0f;
    DiffPair_float_0 _S3920;
    (&_S3920)->primal_0 = _S3567;
    (&_S3920)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3919, &_S3920, k_9);
    float _S3921 = _S3920.differential_0 + _S3570;
    float _S3922 = _S3919.differential_0 + _S3571;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3923;
    (&_S3923)->primal_0 = _S3565;
    (&_S3923)->differential_0 = _S3805;
    s_bwd_length_impl_1(&_S3923, _S3922);
    float2  _S3924 = _S3923.differential_0 + _S3912;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3925;
    (&_S3925)->primal_0 = vert2_c_5;
    (&_S3925)->differential_0 = _S3697;
    s_bwd_length_impl_0(&_S3925, _S3786.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3926;
    (&_S3926)->primal_0 = vert1_c_5;
    (&_S3926)->differential_0 = _S3697;
    s_bwd_length_impl_0(&_S3926, _S3786.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3927;
    (&_S3927)->primal_0 = vert0_c_5;
    (&_S3927)->differential_0 = _S3697;
    s_bwd_length_impl_0(&_S3927, _S3786.differential_0.x);
    float3  _S3928 = _S3925.differential_0 + _S3882;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3929;
    (&_S3929)->primal_0 = R_24;
    (&_S3929)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3930;
    (&_S3930)->primal_0 = vert2_1;
    (&_S3930)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3929, &_S3930, _S3928);
    float3  _S3931 = _S3926.differential_0 + _S3909;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3932;
    (&_S3932)->primal_0 = R_24;
    (&_S3932)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3933;
    (&_S3933)->primal_0 = vert1_1;
    (&_S3933)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3932, &_S3933, _S3931);
    float3  _S3934 = _S3927.differential_0 + _S3832 + make_float3 (_S3924.x, _S3924.y, _S3921);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3935;
    (&_S3935)->primal_0 = R_24;
    (&_S3935)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3936;
    (&_S3936)->primal_0 = vert0_1;
    (&_S3936)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3935, &_S3936, _S3934);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3937;
    (&_S3937)->primal_0 = _S3556;
    (&_S3937)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3938;
    (&_S3938)->primal_0 = _S3561;
    (&_S3938)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3937, &_S3938, _S3930.differential_0);
    float _S3939 = - _S3938.differential_0.y;
    float _S3940 = _S3560 * _S3938.differential_0.x;
    float _S3941 = - (_S3552 * _S3938.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3942;
    (&_S3942)->primal_0 = _S3556;
    (&_S3942)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3943;
    (&_S3943)->primal_0 = _S3559;
    (&_S3943)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3942, &_S3943, _S3933.differential_0);
    float _S3944 = _S3552 * _S3943.differential_0.x;
    float _S3945 = _S3558 * _S3943.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3946;
    (&_S3946)->primal_0 = _S3556;
    (&_S3946)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3947;
    (&_S3947)->primal_0 = _S3557;
    (&_S3947)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3946, &_S3947, _S3936.differential_0);
    Matrix<float, 3, 3>  _S3948 = transpose_0(_S3937.differential_0 + _S3942.differential_0 + _S3946.differential_0);
    float _S3949 = 2.0f * - _S3948.rows[int(2)].z;
    float _S3950 = 2.0f * _S3948.rows[int(2)].y;
    float _S3951 = 2.0f * _S3948.rows[int(2)].x;
    float _S3952 = 2.0f * _S3948.rows[int(1)].z;
    float _S3953 = 2.0f * - _S3948.rows[int(1)].y;
    float _S3954 = 2.0f * _S3948.rows[int(1)].x;
    float _S3955 = 2.0f * _S3948.rows[int(0)].z;
    float _S3956 = 2.0f * _S3948.rows[int(0)].y;
    float _S3957 = 2.0f * - _S3948.rows[int(0)].x;
    float _S3958 = - _S3954 + _S3956;
    float _S3959 = _S3951 + - _S3955;
    float _S3960 = - _S3950 + _S3952;
    float _S3961 = _S3950 + _S3952;
    float _S3962 = _S3951 + _S3955;
    float _S3963 = _S3954 + _S3956;
    float _S3964 = quat_25.w * (_S3953 + _S3957);
    float _S3965 = quat_25.z * (_S3949 + _S3957);
    float _S3966 = quat_25.y * (_S3949 + _S3953);
    float _S3967 = quat_25.x * _S3958 + quat_25.z * _S3961 + quat_25.y * _S3962 + _S3964 + _S3964;
    float _S3968 = quat_25.x * _S3959 + quat_25.w * _S3961 + quat_25.y * _S3963 + _S3965 + _S3965;
    float _S3969 = quat_25.x * _S3960 + quat_25.w * _S3962 + quat_25.z * _S3963 + _S3966 + _S3966;
    float _S3970 = quat_25.w * _S3958 + quat_25.z * _S3959 + quat_25.y * _S3960;
    float _S3971 = _S3941 + _S3944;
    float _S3972 = 0.5f * - _S3971;
    float _S3973 = _S3939 + _S3943.differential_0.y;
    DiffPair_float_0 _S3974;
    (&_S3974)->primal_0 = _S3553;
    (&_S3974)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3974, _S3973);
    float _S3975 = _S3972 + _S3974.differential_0;
    float _S3976 = _S3940 + _S3945 + _S3947.differential_0.x;
    DiffPair_float_0 _S3977;
    (&_S3977)->primal_0 = _S3551;
    (&_S3977)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3977, _S3976);
    float _S3978 = _S3972 + _S3977.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3979;
    (&_S3979)->primal_0 = mean_c_20;
    (&_S3979)->differential_0 = _S3697;
    s_bwd_length_impl_0(&_S3979, 0.0f);
    float3  _S3980 = _S3979.differential_0 + _S3700.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3981;
    (&_S3981)->primal_0 = R_24;
    (&_S3981)->differential_0 = _S3780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3982;
    (&_S3982)->primal_0 = mean_24;
    (&_S3982)->differential_0 = _S3697;
    s_bwd_prop_mul_1(&_S3981, &_S3982, _S3980);
    float3  _S3983 = _S3928 + _S3931 + _S3934 + _S3980 + _S3783.differential_0;
    Matrix<float, 3, 3>  _S3984 = _S3929.differential_0 + _S3932.differential_0 + _S3935.differential_0 + _S3981.differential_0 + _S3784;
    float3  _S3985 = make_float3 (_S3978, _S3975, _S3971);
    float4  _S3986 = make_float4 (0.0f);
    *&((&_S3986)->w) = _S3967;
    *&((&_S3986)->z) = _S3968;
    *&((&_S3986)->y) = _S3969;
    *&((&_S3986)->x) = _S3970;
    float4  _S3987 = _S3986;
    float3  _S3988 = _S3930.differential_0 + _S3933.differential_0 + _S3936.differential_0 + _S3982.differential_0 + _S3778;
    *v_mean_8 = _S3988;
    *v_quat_7 = _S3987;
    *v_scale_7 = _S3985;
    *v_hardness_1 = _S3829;
    (*v_sh_coeffs_6)[int(0)] = _S3837;
    (*v_sh_coeffs_6)[int(1)] = _S3838;
    (*v_sh_coeffs_6)[int(2)] = _S3839;
    (*v_sh_coeffs_6)[int(3)] = _S3840;
    (*v_sh_coeffs_6)[int(4)] = _S3841;
    (*v_sh_coeffs_6)[int(5)] = _S3842;
    (*v_sh_coeffs_6)[int(6)] = _S3843;
    (*v_sh_coeffs_6)[int(7)] = _S3844;
    (*v_sh_coeffs_6)[int(8)] = _S3845;
    (*v_sh_coeffs_6)[int(9)] = _S3846;
    (*v_sh_coeffs_6)[int(10)] = _S3847;
    (*v_sh_coeffs_6)[int(11)] = _S3848;
    (*v_sh_coeffs_6)[int(12)] = _S3849;
    (*v_sh_coeffs_6)[int(13)] = _S3850;
    (*v_sh_coeffs_6)[int(14)] = _S3851;
    (*v_sh_coeffs_6)[int(15)] = _S3852;
    (*v_ch_coeffs_1)[int(0)] = _S3834;
    (*v_ch_coeffs_1)[int(1)] = _S3835;
    *v_R_7 = _S3984;
    *v_t_7 = _S3983;
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
        DiffPair_float_0 _S3989 = *dpx_17;
        float _S3990 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3990;
        float _S3991 = val_0 * (F32_log((_S3989.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3991;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_1)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3992 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3992))));
    float2  _S3993 = p_1 - v0_0;
    float2  _S3994 = normalize_1(e0_6);
    float2  _S3995 = p_1 - v1_0;
    float2  _S3996 = normalize_1(e1_6);
    float2  _S3997 = p_1 - v2_0;
    float2  _S3998 = normalize_1(e2_2);
    float _S3999 = hardness_6.x;
    float _S4000 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3993.x * _S3994.y - _S3993.y * _S3994.x)), (se_0 * (_S3995.x * _S3996.y - _S3995.y * _S3996.x))))), (se_0 * (_S3997.x * _S3998.y - _S3997.y * _S3998.x)))) / ((F32_abs((_S3992))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S4000))));
    float _S4001;
    if(a_1 <= 0.0f)
    {
        _S4001 = 0.0f;
    }
    else
    {
        _S4001 = (F32_min(((F32_pow((a_1), (_S4000)))), (0.99900001287460327f)));
    }
    return _S3999 * _S4001;
}

inline __device__ float s_primal_ctx_abs_0(float _S4002)
{
    return (F32_abs((_S4002)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S4003, float _S4004, float _S4005)
{
    return clamp_0(_S4003, _S4004, _S4005);
}

inline __device__ float s_primal_ctx_exp2_0(float _S4006)
{
    return (F32_exp2((_S4006)));
}

inline __device__ float s_primal_ctx_pow_0(float _S4007, float _S4008)
{
    return (F32_pow((_S4007), (_S4008)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S4009, DiffPair_float_0 * _S4010, float _S4011)
{
    _d_pow_0(_S4009, _S4010, _S4011);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S4012, DiffPair_float_0 * _S4013, DiffPair_float_0 * _S4014, float _S4015)
{
    _d_clamp_0(_S4012, _S4013, _S4014, _S4015);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_8)
{
    float _S4016 = length_0((*dpx_18).primal_0);
    float2  _S4017 = (*dpx_18).primal_0 * _s_dOut_8;
    float2  _S4018 = make_float2 (1.0f / _S4016) * _s_dOut_8;
    float _S4019 = - ((_S4017.x + _S4017.y) / (_S4016 * _S4016));
    float2  _S4020 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4021;
    (&_S4021)->primal_0 = (*dpx_18).primal_0;
    (&_S4021)->differential_0 = _S4020;
    s_bwd_length_impl_1(&_S4021, _S4019);
    float2  _S4022 = _S4018 + _S4021.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S4022;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4023, float2  _S4024)
{
    s_bwd_prop_normalize_impl_1(_S4023, _S4024);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_2, float _s_dOut_9)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S4025 = e0_7.x;
    float _S4026 = e1_7.y;
    float _S4027 = e0_7.y;
    float _S4028 = e1_7.x;
    float _S4029 = _S4025 * _S4026 - _S4027 * _S4028;
    float se_1 = float((F32_sign((_S4029))));
    float2  _S4030 = p_2 - (*dpv0_0).primal_0;
    float2  _S4031 = normalize_1(e0_7);
    float _S4032 = _S4030.x;
    float _S4033 = _S4031.y;
    float _S4034 = _S4030.y;
    float _S4035 = _S4031.x;
    float de0_0 = se_1 * (_S4032 * _S4033 - _S4034 * _S4035);
    float2  _S4036 = p_2 - (*dpv1_0).primal_0;
    float2  _S4037 = normalize_1(e1_7);
    float _S4038 = _S4036.x;
    float _S4039 = _S4037.y;
    float _S4040 = _S4036.y;
    float _S4041 = _S4037.x;
    float de1_0 = se_1 * (_S4038 * _S4039 - _S4040 * _S4041);
    float2  _S4042 = p_2 - (*dpv2_0).primal_0;
    float2  _S4043 = normalize_1(e2_3);
    float _S4044 = _S4042.x;
    float _S4045 = _S4043.y;
    float _S4046 = _S4042.y;
    float _S4047 = _S4043.x;
    float de2_0 = se_1 * (_S4044 * _S4045 - _S4046 * _S4047);
    float _S4048 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4049 = s_primal_ctx_max_0(_S4048, de2_0);
    float _S4050 = s_primal_ctx_abs_0(_S4029);
    float _S4051 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4050 / _S4051;
    float _S4052 = _S4051 * _S4051;
    float _S4053 = (*dphardness_0).primal_0.x;
    float _S4054 = (*dphardness_0).primal_0.y;
    float _S4055 = dmax_0 * dmax_0;
    float _S4056 = 1.0f + _S4049 / dmax_0;
    float _S4057 = 1.0f - s_primal_ctx_clamp_0(_S4054, 0.00499999988824129f, 0.98000001907348633f);
    float _S4058 = -1.0f / _S4057;
    float _S4059 = _S4057 * _S4057;
    float _S4060 = 1.0f - s_primal_ctx_exp2_0(_S4058);
    float a_2 = 1.0f - _S4056 * _S4060;
    bool _S4061 = a_2 <= 0.0f;
    float _S4062;
    float _S4063;
    if(_S4061)
    {
        _S4062 = 0.0f;
        _S4063 = 0.0f;
    }
    else
    {
        float _S4064 = s_primal_ctx_pow_0(a_2, _S4057);
        _S4062 = s_primal_ctx_min_0(_S4064, 0.99900001287460327f);
        _S4063 = _S4064;
    }
    float _S4065 = _S4053 * _s_dOut_9;
    float _S4066 = _S4062 * _s_dOut_9;
    if(_S4061)
    {
        _S4062 = 0.0f;
        _S4063 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4067;
        (&_S4067)->primal_0 = _S4063;
        (&_S4067)->differential_0 = 0.0f;
        DiffPair_float_0 _S4068;
        (&_S4068)->primal_0 = 0.99900001287460327f;
        (&_S4068)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4067, &_S4068, _S4065);
        DiffPair_float_0 _S4069;
        (&_S4069)->primal_0 = a_2;
        (&_S4069)->differential_0 = 0.0f;
        DiffPair_float_0 _S4070;
        (&_S4070)->primal_0 = _S4057;
        (&_S4070)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4069, &_S4070, _S4067.differential_0);
        _S4062 = _S4069.differential_0;
        _S4063 = _S4070.differential_0;
    }
    float _S4071 = - _S4062;
    float _S4072 = _S4060 * _S4071;
    float _S4073 = - (_S4056 * _S4071);
    DiffPair_float_0 _S4074;
    (&_S4074)->primal_0 = _S4058;
    (&_S4074)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4074, _S4073);
    float _S4075 = - (-1.0f * - (_S4074.differential_0 / _S4059) + _S4063);
    float _S4076 = _S4072 / _S4055;
    float s_diff_dmax_T_0 = _S4049 * - _S4076;
    float _S4077 = dmax_0 * _S4076;
    DiffPair_float_0 _S4078;
    (&_S4078)->primal_0 = _S4054;
    (&_S4078)->differential_0 = 0.0f;
    DiffPair_float_0 _S4079;
    (&_S4079)->primal_0 = 0.00499999988824129f;
    (&_S4079)->differential_0 = 0.0f;
    DiffPair_float_0 _S4080;
    (&_S4080)->primal_0 = 0.98000001907348633f;
    (&_S4080)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4078, &_S4079, &_S4080, _S4075);
    float _S4081 = s_diff_dmax_T_0 / _S4052;
    float _S4082 = _S4050 * - _S4081;
    float _S4083 = _S4051 * _S4081;
    float2  _S4084 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4085;
    (&_S4085)->primal_0 = e2_3;
    (&_S4085)->differential_0 = _S4084;
    s_bwd_length_impl_1(&_S4085, _S4082);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4086;
    (&_S4086)->primal_0 = e1_7;
    (&_S4086)->differential_0 = _S4084;
    s_bwd_length_impl_1(&_S4086, _S4082);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4087;
    (&_S4087)->primal_0 = e0_7;
    (&_S4087)->differential_0 = _S4084;
    s_bwd_length_impl_1(&_S4087, _S4082);
    DiffPair_float_0 _S4088;
    (&_S4088)->primal_0 = _S4029;
    (&_S4088)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4088, _S4083);
    DiffPair_float_0 _S4089;
    (&_S4089)->primal_0 = _S4048;
    (&_S4089)->differential_0 = 0.0f;
    DiffPair_float_0 _S4090;
    (&_S4090)->primal_0 = de2_0;
    (&_S4090)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4089, &_S4090, _S4077);
    DiffPair_float_0 _S4091;
    (&_S4091)->primal_0 = de0_0;
    (&_S4091)->differential_0 = 0.0f;
    DiffPair_float_0 _S4092;
    (&_S4092)->primal_0 = de1_0;
    (&_S4092)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4091, &_S4092, _S4089.differential_0);
    float _S4093 = se_1 * _S4090.differential_0;
    float _S4094 = - _S4093;
    float _S4095 = _S4047 * _S4094;
    float _S4096 = _S4045 * _S4093;
    float2  _S4097 = make_float2 (_S4046 * _S4094, _S4044 * _S4093);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4098;
    (&_S4098)->primal_0 = e2_3;
    (&_S4098)->differential_0 = _S4084;
    s_bwd_normalize_impl_1(&_S4098, _S4097);
    float2  _S4099 = - make_float2 (_S4096, _S4095);
    float _S4100 = se_1 * _S4092.differential_0;
    float _S4101 = - _S4100;
    float _S4102 = _S4041 * _S4101;
    float _S4103 = _S4039 * _S4100;
    float2  _S4104 = make_float2 (_S4040 * _S4101, _S4038 * _S4100);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4105;
    (&_S4105)->primal_0 = e1_7;
    (&_S4105)->differential_0 = _S4084;
    s_bwd_normalize_impl_1(&_S4105, _S4104);
    float2  _S4106 = - make_float2 (_S4103, _S4102);
    float _S4107 = se_1 * _S4091.differential_0;
    float _S4108 = - _S4107;
    float _S4109 = _S4035 * _S4108;
    float _S4110 = _S4033 * _S4107;
    float2  _S4111 = make_float2 (_S4034 * _S4108, _S4032 * _S4107);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4112;
    (&_S4112)->primal_0 = e0_7;
    (&_S4112)->differential_0 = _S4084;
    s_bwd_normalize_impl_1(&_S4112, _S4111);
    float2  _S4113 = - make_float2 (_S4110, _S4109);
    float _S4114 = - _S4088.differential_0;
    float2  _S4115 = _S4085.differential_0 + _S4098.differential_0;
    float2  _S4116 = - _S4115;
    float2  _S4117 = _S4086.differential_0 + _S4105.differential_0 + make_float2 (_S4027 * _S4114, _S4025 * _S4088.differential_0);
    float2  _S4118 = - _S4117;
    float2  _S4119 = _S4087.differential_0 + _S4112.differential_0 + make_float2 (_S4026 * _S4088.differential_0, _S4028 * _S4114);
    float2  _S4120 = - _S4119;
    float2  _S4121 = make_float2 (_S4066, _S4078.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4121;
    float2  _S4122 = _S4099 + _S4116 + _S4117;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4122;
    float2  _S4123 = _S4106 + _S4118 + _S4119;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4123;
    float2  _S4124 = _S4113 + _S4115 + _S4120;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4124;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4125, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4126, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4127, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4128, float2  _S4129, float _S4130)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4125, _S4126, _S4127, _S4128, _S4129, _S4130);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_3, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4131 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4131;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4131;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4131;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4131;
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
    float2  _S4132 = p_4 - v0_2;
    float2  _S4133 = p_4 - v1_2;
    float2  _S4134 = p_4 - v2_2;
    float _S4135 = e0_8.x;
    float _S4136 = e1_8.y;
    float _S4137 = e0_8.y;
    float _S4138 = e1_8.x;
    float _S4139 = _S4135 * _S4136 - _S4137 * _S4138;
    float se_2 = float((F32_sign((_S4139))));
    float _S4140 = hardness_8.x;
    float _S4141 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4132.x * _S4137 - _S4132.y * _S4135)), (se_2 * (_S4133.x * _S4136 - _S4133.y * _S4138))))), (se_2 * (_S4134.x * e2_4.y - _S4134.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4132 - e0_8 * make_float2 (clamp_0(dot_1(_S4132, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4133 - e1_8 * make_float2 (clamp_0(dot_1(_S4133, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4134 - e2_4 * make_float2 (clamp_0(dot_1(_S4134, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4139))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4141))));
    float _S4142;
    if(a_3 <= 0.0f)
    {
        _S4142 = 0.0f;
    }
    else
    {
        _S4142 = (F32_min(((F32_pow((a_3), (_S4141)))), (0.99900001287460327f)));
    }
    return _S4140 * _S4142;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4143, float2  _S4144)
{
    return dot_1(_S4143, _S4144);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4145, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4146, float _S4147)
{
    _d_dot_1(_S4145, _S4146, _S4147);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_5, float _s_dOut_10)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4148 = p_5 - (*dpv0_1).primal_0;
    float _S4149 = s_primal_ctx_dot_1(_S4148, e0_9);
    float _S4150 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4151 = _S4149 / _S4150;
    float _S4152 = _S4150 * _S4150;
    float _S4153 = s_primal_ctx_clamp_0(_S4151, 0.0f, 1.0f);
    float2  _S4154 = make_float2 (_S4153);
    float2  _S4155 = _S4148 - e0_9 * make_float2 (_S4153);
    float _S4156 = length_0(_S4155);
    float2  _S4157 = p_5 - (*dpv1_1).primal_0;
    float _S4158 = s_primal_ctx_dot_1(_S4157, e1_9);
    float _S4159 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4160 = _S4158 / _S4159;
    float _S4161 = _S4159 * _S4159;
    float _S4162 = s_primal_ctx_clamp_0(_S4160, 0.0f, 1.0f);
    float2  _S4163 = make_float2 (_S4162);
    float2  _S4164 = _S4157 - e1_9 * make_float2 (_S4162);
    float _S4165 = length_0(_S4164);
    float2  _S4166 = p_5 - (*dpv2_1).primal_0;
    float _S4167 = s_primal_ctx_dot_1(_S4166, e2_5);
    float _S4168 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4169 = _S4167 / _S4168;
    float _S4170 = _S4168 * _S4168;
    float _S4171 = s_primal_ctx_clamp_0(_S4169, 0.0f, 1.0f);
    float2  _S4172 = make_float2 (_S4171);
    float2  _S4173 = _S4166 - e2_5 * make_float2 (_S4171);
    float _S4174 = length_0(_S4173);
    float _S4175 = e0_9.x;
    float _S4176 = e1_9.y;
    float _S4177 = e0_9.y;
    float _S4178 = e1_9.x;
    float _S4179 = _S4175 * _S4176 - _S4177 * _S4178;
    float se_3 = float((F32_sign((_S4179))));
    float _S4180 = _S4148.x;
    float _S4181 = _S4148.y;
    float s0_0 = se_3 * (_S4180 * _S4177 - _S4181 * _S4175);
    float _S4182 = _S4157.x;
    float _S4183 = _S4157.y;
    float s1_0 = se_3 * (_S4182 * _S4176 - _S4183 * _S4178);
    float _S4184 = _S4166.x;
    float _S4185 = e2_5.y;
    float _S4186 = _S4166.y;
    float _S4187 = e2_5.x;
    float s2_0 = se_3 * (_S4184 * _S4185 - _S4186 * _S4187);
    float _S4188 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4188, s2_0)))));
    float _S4189 = s_primal_ctx_min_0(_S4156, _S4165);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4189, _S4174);
    float _S4190 = s_primal_ctx_abs_0(_S4179);
    float _S4191 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4190 / _S4191;
    float _S4192 = _S4191 * _S4191;
    float _S4193 = (*dphardness_1).primal_0.x;
    float _S4194 = (*dphardness_1).primal_0.y;
    float _S4195 = dmax_1 * dmax_1;
    float _S4196 = 1.0f + dv_0 / dmax_1;
    float _S4197 = 1.0f - s_primal_ctx_clamp_0(_S4194, 0.00499999988824129f, 0.98000001907348633f);
    float _S4198 = -1.0f / _S4197;
    float _S4199 = _S4197 * _S4197;
    float _S4200 = 1.0f - s_primal_ctx_exp2_0(_S4198);
    float a_4 = 1.0f - _S4196 * _S4200;
    bool _S4201 = a_4 <= 0.0f;
    float _S4202;
    float _S4203;
    if(_S4201)
    {
        _S4202 = 0.0f;
        _S4203 = 0.0f;
    }
    else
    {
        float _S4204 = s_primal_ctx_pow_0(a_4, _S4197);
        _S4202 = s_primal_ctx_min_0(_S4204, 0.99900001287460327f);
        _S4203 = _S4204;
    }
    float _S4205 = _S4193 * _s_dOut_10;
    float _S4206 = _S4202 * _s_dOut_10;
    if(_S4201)
    {
        _S4202 = 0.0f;
        _S4203 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4207;
        (&_S4207)->primal_0 = _S4203;
        (&_S4207)->differential_0 = 0.0f;
        DiffPair_float_0 _S4208;
        (&_S4208)->primal_0 = 0.99900001287460327f;
        (&_S4208)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4207, &_S4208, _S4205);
        DiffPair_float_0 _S4209;
        (&_S4209)->primal_0 = a_4;
        (&_S4209)->differential_0 = 0.0f;
        DiffPair_float_0 _S4210;
        (&_S4210)->primal_0 = _S4197;
        (&_S4210)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4209, &_S4210, _S4207.differential_0);
        _S4202 = _S4209.differential_0;
        _S4203 = _S4210.differential_0;
    }
    float _S4211 = - _S4202;
    float _S4212 = _S4200 * _S4211;
    float _S4213 = - (_S4196 * _S4211);
    DiffPair_float_0 _S4214;
    (&_S4214)->primal_0 = _S4198;
    (&_S4214)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4214, _S4213);
    float _S4215 = - (-1.0f * - (_S4214.differential_0 / _S4199) + _S4203);
    float _S4216 = _S4212 / _S4195;
    float s_diff_dmax_T_1 = dv_0 * - _S4216;
    float s_diff_dv_T_0 = dmax_1 * _S4216;
    DiffPair_float_0 _S4217;
    (&_S4217)->primal_0 = _S4194;
    (&_S4217)->differential_0 = 0.0f;
    DiffPair_float_0 _S4218;
    (&_S4218)->primal_0 = 0.00499999988824129f;
    (&_S4218)->differential_0 = 0.0f;
    DiffPair_float_0 _S4219;
    (&_S4219)->primal_0 = 0.98000001907348633f;
    (&_S4219)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4217, &_S4218, &_S4219, _S4215);
    float _S4220 = s_diff_dmax_T_1 / _S4192;
    float _S4221 = _S4190 * - _S4220;
    float _S4222 = _S4191 * _S4220;
    float2  _S4223 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4224;
    (&_S4224)->primal_0 = e2_5;
    (&_S4224)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4224, _S4221);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4225;
    (&_S4225)->primal_0 = e1_9;
    (&_S4225)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4225, _S4221);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4226;
    (&_S4226)->primal_0 = e0_9;
    (&_S4226)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4226, _S4221);
    DiffPair_float_0 _S4227;
    (&_S4227)->primal_0 = _S4179;
    (&_S4227)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4227, _S4222);
    float _S4228 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4229;
    (&_S4229)->primal_0 = _S4189;
    (&_S4229)->differential_0 = 0.0f;
    DiffPair_float_0 _S4230;
    (&_S4230)->primal_0 = _S4174;
    (&_S4230)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4229, &_S4230, _S4228);
    DiffPair_float_0 _S4231;
    (&_S4231)->primal_0 = _S4156;
    (&_S4231)->differential_0 = 0.0f;
    DiffPair_float_0 _S4232;
    (&_S4232)->primal_0 = _S4165;
    (&_S4232)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4231, &_S4232, _S4229.differential_0);
    DiffPair_float_0 _S4233;
    (&_S4233)->primal_0 = _S4188;
    (&_S4233)->differential_0 = 0.0f;
    DiffPair_float_0 _S4234;
    (&_S4234)->primal_0 = s2_0;
    (&_S4234)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4233, &_S4234, 0.0f);
    DiffPair_float_0 _S4235;
    (&_S4235)->primal_0 = s0_0;
    (&_S4235)->differential_0 = 0.0f;
    DiffPair_float_0 _S4236;
    (&_S4236)->primal_0 = s1_0;
    (&_S4236)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4235, &_S4236, _S4233.differential_0);
    float _S4237 = se_3 * _S4234.differential_0;
    float _S4238 = - _S4237;
    float _S4239 = _S4186 * _S4238;
    float _S4240 = _S4187 * _S4238;
    float _S4241 = _S4184 * _S4237;
    float _S4242 = _S4185 * _S4237;
    float _S4243 = se_3 * _S4236.differential_0;
    float _S4244 = - _S4243;
    float _S4245 = _S4178 * _S4244;
    float _S4246 = _S4176 * _S4243;
    float _S4247 = se_3 * _S4235.differential_0;
    float _S4248 = - _S4247;
    float _S4249 = _S4175 * _S4248;
    float _S4250 = _S4177 * _S4247;
    float _S4251 = - _S4227.differential_0;
    float _S4252 = _S4183 * _S4244 + _S4177 * _S4251;
    float _S4253 = _S4180 * _S4247 + _S4178 * _S4251;
    float _S4254 = _S4182 * _S4243 + _S4175 * _S4227.differential_0;
    float _S4255 = _S4181 * _S4248 + _S4176 * _S4227.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4256;
    (&_S4256)->primal_0 = _S4173;
    (&_S4256)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4256, _S4230.differential_0);
    float2  _S4257 = - _S4256.differential_0;
    float2  _S4258 = e2_5 * _S4257;
    float2  _S4259 = _S4172 * _S4257;
    float _S4260 = _S4258.x + _S4258.y;
    DiffPair_float_0 _S4261;
    (&_S4261)->primal_0 = _S4169;
    (&_S4261)->differential_0 = 0.0f;
    DiffPair_float_0 _S4262;
    (&_S4262)->primal_0 = 0.0f;
    (&_S4262)->differential_0 = 0.0f;
    DiffPair_float_0 _S4263;
    (&_S4263)->primal_0 = 1.0f;
    (&_S4263)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4261, &_S4262, &_S4263, _S4260);
    float _S4264 = _S4261.differential_0 / _S4170;
    float _S4265 = _S4167 * - _S4264;
    float _S4266 = _S4168 * _S4264;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4267;
    (&_S4267)->primal_0 = e2_5;
    (&_S4267)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4268;
    (&_S4268)->primal_0 = e2_5;
    (&_S4268)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4267, &_S4268, _S4265);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4269;
    (&_S4269)->primal_0 = _S4166;
    (&_S4269)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4270;
    (&_S4270)->primal_0 = e2_5;
    (&_S4270)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4269, &_S4270, _S4266);
    float2  _S4271 = - (_S4256.differential_0 + _S4269.differential_0 + make_float2 (_S4242, _S4240));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4272;
    (&_S4272)->primal_0 = _S4164;
    (&_S4272)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4272, _S4232.differential_0);
    float2  _S4273 = - _S4272.differential_0;
    float2  _S4274 = e1_9 * _S4273;
    float2  _S4275 = _S4163 * _S4273;
    float _S4276 = _S4274.x + _S4274.y;
    DiffPair_float_0 _S4277;
    (&_S4277)->primal_0 = _S4160;
    (&_S4277)->differential_0 = 0.0f;
    DiffPair_float_0 _S4278;
    (&_S4278)->primal_0 = 0.0f;
    (&_S4278)->differential_0 = 0.0f;
    DiffPair_float_0 _S4279;
    (&_S4279)->primal_0 = 1.0f;
    (&_S4279)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4277, &_S4278, &_S4279, _S4276);
    float _S4280 = _S4277.differential_0 / _S4161;
    float _S4281 = _S4158 * - _S4280;
    float _S4282 = _S4159 * _S4280;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4283;
    (&_S4283)->primal_0 = e1_9;
    (&_S4283)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4284;
    (&_S4284)->primal_0 = e1_9;
    (&_S4284)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4283, &_S4284, _S4281);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4285;
    (&_S4285)->primal_0 = _S4157;
    (&_S4285)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4286;
    (&_S4286)->primal_0 = e1_9;
    (&_S4286)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4285, &_S4286, _S4282);
    float2  _S4287 = - (_S4272.differential_0 + _S4285.differential_0 + make_float2 (_S4246, _S4245));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4288;
    (&_S4288)->primal_0 = _S4155;
    (&_S4288)->differential_0 = _S4223;
    s_bwd_length_impl_1(&_S4288, _S4231.differential_0);
    float2  _S4289 = - _S4288.differential_0;
    float2  _S4290 = e0_9 * _S4289;
    float2  _S4291 = _S4154 * _S4289;
    float _S4292 = _S4290.x + _S4290.y;
    DiffPair_float_0 _S4293;
    (&_S4293)->primal_0 = _S4151;
    (&_S4293)->differential_0 = 0.0f;
    DiffPair_float_0 _S4294;
    (&_S4294)->primal_0 = 0.0f;
    (&_S4294)->differential_0 = 0.0f;
    DiffPair_float_0 _S4295;
    (&_S4295)->primal_0 = 1.0f;
    (&_S4295)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4293, &_S4294, &_S4295, _S4292);
    float _S4296 = _S4293.differential_0 / _S4152;
    float _S4297 = _S4149 * - _S4296;
    float _S4298 = _S4150 * _S4296;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4299;
    (&_S4299)->primal_0 = e0_9;
    (&_S4299)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4300;
    (&_S4300)->primal_0 = e0_9;
    (&_S4300)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4299, &_S4300, _S4297);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4301;
    (&_S4301)->primal_0 = _S4148;
    (&_S4301)->differential_0 = _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4302;
    (&_S4302)->primal_0 = e0_9;
    (&_S4302)->differential_0 = _S4223;
    s_bwd_prop_dot_1(&_S4301, &_S4302, _S4298);
    float2  _S4303 = - (_S4288.differential_0 + _S4301.differential_0 + make_float2 (_S4250, _S4249));
    float2  _S4304 = _S4224.differential_0 + _S4259 + _S4268.differential_0 + _S4267.differential_0 + _S4270.differential_0 + make_float2 (_S4239, _S4241);
    float2  _S4305 = - _S4304;
    float2  _S4306 = _S4225.differential_0 + _S4275 + _S4284.differential_0 + _S4283.differential_0 + _S4286.differential_0 + make_float2 (_S4252, _S4254);
    float2  _S4307 = - _S4306;
    float2  _S4308 = _S4226.differential_0 + _S4291 + _S4300.differential_0 + _S4299.differential_0 + _S4302.differential_0 + make_float2 (_S4255, _S4253);
    float2  _S4309 = - _S4308;
    float2  _S4310 = make_float2 (_S4206, _S4217.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4310;
    float2  _S4311 = _S4271 + _S4305 + _S4306;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4311;
    float2  _S4312 = _S4287 + _S4307 + _S4308;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4312;
    float2  _S4313 = _S4303 + _S4304 + _S4309;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4313;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4314, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4315, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4316, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4317, float2  _S4318, float _S4319)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4314, _S4315, _S4316, _S4317, _S4318, _S4319);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_6, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4320 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4320;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4320;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4320;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4320;
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
    float _S4321 = 0.3333333432674408f * dpdepth_1;
    float3  _S4322 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4323 = make_float3 (0.0f);
    float3  _S4324 = _S4323;
    *&((&_S4324)->z) = _S4321;
    *&((&_S4324)->y) = _S4321;
    *&((&_S4324)->x) = _S4321;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4324;
    FixedArray<float3 , 3>  _S4325;
    _S4325[int(0)] = _S4323;
    _S4325[int(1)] = _S4323;
    _S4325[int(2)] = _S4323;
    _S4325[int(2)] = _S4322;
    _S4325[int(1)] = _S4322;
    _S4325[int(0)] = _S4322;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4325;
    float2  _S4326 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4326;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4326;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4326;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4327, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4328, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4329, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4330, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4331, float2  _S4332, float3  _S4333, float _S4334)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4327, _S4328, _S4329, _S4330, _S4331, _S4332, _S4333, _S4334);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_9, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4335 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4335;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4335;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4335;
    float3  _S4336 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4337 = { _S4336, _S4336, _S4336 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4337;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4336;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_9, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_30, float fy_30, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_25) + t_24;
        float _S4338 = mean_c_21.z;
        bool _S4339;
        if(_S4338 < near_plane_14)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4338 > far_plane_14;
        }
        if(_S4339)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4340 = scale_25.x;
        float sx_4 = (F32_exp((_S4340)));
        float _S4341 = scale_25.y;
        float sy_4 = (F32_exp((_S4341)));
        float sz_6 = scale_25.z - 0.5f * (_S4340 + _S4341);
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
        Matrix<float, 3, 3>  _S4342 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
        float3  vert0_2 = mul_0(_S4342, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
        float3  vert1_2 = mul_0(_S4342, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
        float3  vert2_2 = mul_0(_S4342, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
        float3  vert0_c_6 = mul_0(R_25, vert0_2) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_2) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_2) + t_24;
        float _S4343 = vert0_c_6.z;
        float _S4344 = vert1_c_6.z;
        float _S4345 = vert2_c_6.z;
        if(_S4343 < near_plane_14)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4343 > far_plane_14;
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4344 < near_plane_14;
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4344 > far_plane_14;
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4345 < near_plane_14;
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = _S4345 > far_plane_14;
        }
        if(_S4339)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_4;
        for(;;)
        {
            float2  uv0_5 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4343);
            if(_S4343 < 0.0f)
            {
                _S4339 = true;
            }
            else
            {
                bool _S4346 = is_valid_distortion(uv0_5, dist_coeffs_35);
                _S4339 = !_S4346;
            }
            if(_S4339)
            {
                uv0_4 = uv0_5;
                break;
            }
            float u_62 = uv0_5.x;
            float v_62 = uv0_5.y;
            float r2_62 = u_62 * u_62 + v_62 * v_62;
            float2  _S4347 = uv0_5 * make_float2 (1.0f + r2_62 * ((*dist_coeffs_35)[int(0)] + r2_62 * ((*dist_coeffs_35)[int(1)] + r2_62 * ((*dist_coeffs_35)[int(2)] + r2_62 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_62 * v_62 + (*dist_coeffs_35)[int(5)] * (r2_62 + 2.0f * u_62 * u_62) + (*dist_coeffs_35)[int(6)] * r2_62, 2.0f * (*dist_coeffs_35)[int(5)] * u_62 * v_62 + (*dist_coeffs_35)[int(4)] * (r2_62 + 2.0f * v_62 * v_62) + (*dist_coeffs_35)[int(7)] * r2_62);
            float2  _S4348 = _S4347 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4347.x + (*dist_coeffs_35)[int(9)] * _S4347.y, 0.0f);
            uv0_4 = make_float2 (fx_30 * _S4348.x + cx_26, fy_30 * _S4348.y + cy_26);
            break;
        }
        float2  uv1_4;
        bool all_valid_12 = true & (!_S4339);
        for(;;)
        {
            float2  uv1_5 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4344);
            if(_S4344 < 0.0f)
            {
                _S4339 = true;
            }
            else
            {
                bool _S4349 = is_valid_distortion(uv1_5, dist_coeffs_35);
                _S4339 = !_S4349;
            }
            if(_S4339)
            {
                uv1_4 = uv1_5;
                break;
            }
            float u_63 = uv1_5.x;
            float v_63 = uv1_5.y;
            float r2_63 = u_63 * u_63 + v_63 * v_63;
            float2  _S4350 = uv1_5 * make_float2 (1.0f + r2_63 * ((*dist_coeffs_35)[int(0)] + r2_63 * ((*dist_coeffs_35)[int(1)] + r2_63 * ((*dist_coeffs_35)[int(2)] + r2_63 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_63 * v_63 + (*dist_coeffs_35)[int(5)] * (r2_63 + 2.0f * u_63 * u_63) + (*dist_coeffs_35)[int(6)] * r2_63, 2.0f * (*dist_coeffs_35)[int(5)] * u_63 * v_63 + (*dist_coeffs_35)[int(4)] * (r2_63 + 2.0f * v_63 * v_63) + (*dist_coeffs_35)[int(7)] * r2_63);
            float2  _S4351 = _S4350 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4350.x + (*dist_coeffs_35)[int(9)] * _S4350.y, 0.0f);
            uv1_4 = make_float2 (fx_30 * _S4351.x + cx_26, fy_30 * _S4351.y + cy_26);
            break;
        }
        float2  uv2_4;
        bool all_valid_13 = all_valid_12 & (!_S4339);
        for(;;)
        {
            float2  uv2_5 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4345);
            if(_S4345 < 0.0f)
            {
                _S4339 = true;
            }
            else
            {
                bool _S4352 = is_valid_distortion(uv2_5, dist_coeffs_35);
                _S4339 = !_S4352;
            }
            if(_S4339)
            {
                uv2_4 = uv2_5;
                break;
            }
            float u_64 = uv2_5.x;
            float v_64 = uv2_5.y;
            float r2_64 = u_64 * u_64 + v_64 * v_64;
            float2  _S4353 = uv2_5 * make_float2 (1.0f + r2_64 * ((*dist_coeffs_35)[int(0)] + r2_64 * ((*dist_coeffs_35)[int(1)] + r2_64 * ((*dist_coeffs_35)[int(2)] + r2_64 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_64 * v_64 + (*dist_coeffs_35)[int(5)] * (r2_64 + 2.0f * u_64 * u_64) + (*dist_coeffs_35)[int(6)] * r2_64, 2.0f * (*dist_coeffs_35)[int(5)] * u_64 * v_64 + (*dist_coeffs_35)[int(4)] * (r2_64 + 2.0f * v_64 * v_64) + (*dist_coeffs_35)[int(7)] * r2_64);
            float2  _S4354 = _S4353 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4353.x + (*dist_coeffs_35)[int(9)] * _S4353.y, 0.0f);
            uv2_4 = make_float2 (fx_30 * _S4354.x + cx_26, fy_30 * _S4354.y + cy_26);
            break;
        }
        if(!(all_valid_13 & (!_S4339)))
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_10 = uv1_4 - uv0_4;
        float2  e1_10 = uv2_4 - uv1_4;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(uv0_4 - uv2_4)));
        float _S4355 = uv0_4.x;
        float _S4356 = uv1_4.x;
        float _S4357 = uv2_4.x;
        float xmax_7 = (F32_max(((F32_max((_S4355), (_S4356)))), (_S4357))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4355), (_S4356)))), (_S4357))) - offset_4;
        float _S4358 = uv0_4.y;
        float _S4359 = uv1_4.y;
        float _S4360 = uv2_4.y;
        float ymax_7 = (F32_max(((F32_max((_S4358), (_S4359)))), (_S4360))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4358), (_S4359)))), (_S4360))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = xmin_7 >= float(image_width_21);
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = ymax_7 <= 0.0f;
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            _S4339 = ymin_7 >= float(image_height_21);
        }
        if(_S4339)
        {
            _S4339 = true;
        }
        else
        {
            if(_S4338 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4339 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4339 = false;
                }
                if(_S4339)
                {
                    _S4339 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4339 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4339 = false;
                    }
                }
            }
            else
            {
                _S4339 = false;
            }
        }
        if(_S4339)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4361 = mean_25 - - mul_0(transpose_0(R_25), t_24);
        float _S4362 = _S4361.x;
        float _S4363 = _S4361.y;
        float _S4364 = _S4361.z;
        float norm_14 = (F32_sqrt((_S4362 * _S4362 + _S4363 * _S4363 + _S4364 * _S4364)));
        float x_53 = _S4362 / norm_14;
        float y_24 = _S4363 / norm_14;
        float z_21 = _S4364 / norm_14;
        float z2_48 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_53 * x_53 - y_24 * y_24;
        float fS1_21 = 2.0f * x_53 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_48 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_53) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_48 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_53) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_53 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_48 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_53) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_53 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4365 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4365);
        float3  _S4366 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4367 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4366 + _S4367 + make_float3 (0.5f), _S4365);
        (*rgbs_0)[int(2)] = max_0(_S4366 - _S4367 + make_float3 (0.5f), _S4365);
        (*verts_0)[int(0)] = vert0_2;
        (*verts_0)[int(1)] = vert1_2;
        (*verts_0)[int(2)] = vert2_2;
        float3  _S4368 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4368 * make_float3 (float(- (F32_sign((dot_0(_S4368, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_31, float fy_31, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    bool _S4369;
    bool _S4370;
    bool _S4371;
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_26) + t_25;
        float _S4372 = length_1(mean_c_22);
        bool _S4373;
        if(_S4372 < near_plane_15)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4372 > far_plane_15;
        }
        if(_S4373)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4374 = scale_26.x;
        float sx_5 = (F32_exp((_S4374)));
        float _S4375 = scale_26.y;
        float sy_5 = (F32_exp((_S4375)));
        float sz_7 = scale_26.z - 0.5f * (_S4374 + _S4375);
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
        Matrix<float, 3, 3>  _S4376 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_49), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_49), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4376, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4376, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4376, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_7 = mul_0(R_26, vert0_3) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_3) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_3) + t_25;
        float _S4377 = length_1(vert0_c_7);
        float _S4378 = length_1(vert1_c_7);
        float _S4379 = length_1(vert2_c_7);
        if(_S4377 < near_plane_15)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4377 > far_plane_15;
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4378 < near_plane_15;
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4378 > far_plane_15;
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4379 < near_plane_15;
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = _S4379 > far_plane_15;
        }
        if(_S4373)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_6;
        float k_10;
        for(;;)
        {
            float2  _S4380 = float2 {vert0_c_7.x, vert0_c_7.y};
            float r_30 = length_0(_S4380);
            float _S4381 = vert0_c_7.z;
            float theta_25 = (F32_atan2((r_30), (_S4381)));
            if(theta_25 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_25 * theta_25 / 3.0f) / _S4381;
            }
            else
            {
                k_10 = theta_25 / r_30;
            }
            float2  uv0_7 = _S4380 * make_float2 (k_10);
            bool _S4382 = is_valid_distortion(uv0_7, dist_coeffs_36);
            bool _S4383 = !_S4382;
            _S4369 = _S4383;
            if(_S4383)
            {
                uv0_6 = uv0_7;
                break;
            }
            float u_65 = uv0_7.x;
            float v_65 = uv0_7.y;
            float r2_65 = u_65 * u_65 + v_65 * v_65;
            float2  _S4384 = uv0_7 * make_float2 (1.0f + r2_65 * ((*dist_coeffs_36)[int(0)] + r2_65 * ((*dist_coeffs_36)[int(1)] + r2_65 * ((*dist_coeffs_36)[int(2)] + r2_65 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_65 * v_65 + (*dist_coeffs_36)[int(5)] * (r2_65 + 2.0f * u_65 * u_65) + (*dist_coeffs_36)[int(6)] * r2_65, 2.0f * (*dist_coeffs_36)[int(5)] * u_65 * v_65 + (*dist_coeffs_36)[int(4)] * (r2_65 + 2.0f * v_65 * v_65) + (*dist_coeffs_36)[int(7)] * r2_65);
            float2  _S4385 = _S4384 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4384.x + (*dist_coeffs_36)[int(9)] * _S4384.y, 0.0f);
            uv0_6 = make_float2 (fx_31 * _S4385.x + cx_27, fy_31 * _S4385.y + cy_27);
            break;
        }
        float2  uv1_6;
        bool all_valid_14 = true & (!_S4369);
        for(;;)
        {
            float2  _S4386 = float2 {vert1_c_7.x, vert1_c_7.y};
            float r_31 = length_0(_S4386);
            float _S4387 = vert1_c_7.z;
            float theta_26 = (F32_atan2((r_31), (_S4387)));
            if(theta_26 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_26 * theta_26 / 3.0f) / _S4387;
            }
            else
            {
                k_10 = theta_26 / r_31;
            }
            float2  uv1_7 = _S4386 * make_float2 (k_10);
            bool _S4388 = is_valid_distortion(uv1_7, dist_coeffs_36);
            bool _S4389 = !_S4388;
            _S4370 = _S4389;
            if(_S4389)
            {
                uv1_6 = uv1_7;
                break;
            }
            float u_66 = uv1_7.x;
            float v_66 = uv1_7.y;
            float r2_66 = u_66 * u_66 + v_66 * v_66;
            float2  _S4390 = uv1_7 * make_float2 (1.0f + r2_66 * ((*dist_coeffs_36)[int(0)] + r2_66 * ((*dist_coeffs_36)[int(1)] + r2_66 * ((*dist_coeffs_36)[int(2)] + r2_66 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_66 * v_66 + (*dist_coeffs_36)[int(5)] * (r2_66 + 2.0f * u_66 * u_66) + (*dist_coeffs_36)[int(6)] * r2_66, 2.0f * (*dist_coeffs_36)[int(5)] * u_66 * v_66 + (*dist_coeffs_36)[int(4)] * (r2_66 + 2.0f * v_66 * v_66) + (*dist_coeffs_36)[int(7)] * r2_66);
            float2  _S4391 = _S4390 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4390.x + (*dist_coeffs_36)[int(9)] * _S4390.y, 0.0f);
            uv1_6 = make_float2 (fx_31 * _S4391.x + cx_27, fy_31 * _S4391.y + cy_27);
            break;
        }
        float2  uv2_6;
        bool all_valid_15 = all_valid_14 & (!_S4370);
        for(;;)
        {
            float2  _S4392 = float2 {vert2_c_7.x, vert2_c_7.y};
            float r_32 = length_0(_S4392);
            float _S4393 = vert2_c_7.z;
            float theta_27 = (F32_atan2((r_32), (_S4393)));
            if(theta_27 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_27 * theta_27 / 3.0f) / _S4393;
            }
            else
            {
                k_10 = theta_27 / r_32;
            }
            float2  uv2_7 = _S4392 * make_float2 (k_10);
            bool _S4394 = is_valid_distortion(uv2_7, dist_coeffs_36);
            bool _S4395 = !_S4394;
            _S4371 = _S4395;
            if(_S4395)
            {
                uv2_6 = uv2_7;
                break;
            }
            float u_67 = uv2_7.x;
            float v_67 = uv2_7.y;
            float r2_67 = u_67 * u_67 + v_67 * v_67;
            float2  _S4396 = uv2_7 * make_float2 (1.0f + r2_67 * ((*dist_coeffs_36)[int(0)] + r2_67 * ((*dist_coeffs_36)[int(1)] + r2_67 * ((*dist_coeffs_36)[int(2)] + r2_67 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_67 * v_67 + (*dist_coeffs_36)[int(5)] * (r2_67 + 2.0f * u_67 * u_67) + (*dist_coeffs_36)[int(6)] * r2_67, 2.0f * (*dist_coeffs_36)[int(5)] * u_67 * v_67 + (*dist_coeffs_36)[int(4)] * (r2_67 + 2.0f * v_67 * v_67) + (*dist_coeffs_36)[int(7)] * r2_67);
            float2  _S4397 = _S4396 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4396.x + (*dist_coeffs_36)[int(9)] * _S4396.y, 0.0f);
            uv2_6 = make_float2 (fx_31 * _S4397.x + cx_27, fy_31 * _S4397.y + cy_27);
            break;
        }
        if(!(all_valid_15 & (!_S4371)))
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_11 = uv1_6 - uv0_6;
        float2  e1_11 = uv2_6 - uv1_6;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(uv0_6 - uv2_6)));
        float _S4398 = uv0_6.x;
        float _S4399 = uv1_6.x;
        float _S4400 = uv2_6.x;
        float xmax_8 = (F32_max(((F32_max((_S4398), (_S4399)))), (_S4400))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4398), (_S4399)))), (_S4400))) - offset_5;
        float _S4401 = uv0_6.y;
        float _S4402 = uv1_6.y;
        float _S4403 = uv2_6.y;
        float ymax_8 = (F32_max(((F32_max((_S4401), (_S4402)))), (_S4403))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4401), (_S4402)))), (_S4403))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = xmin_8 >= float(image_width_22);
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = ymax_8 <= 0.0f;
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            _S4373 = ymin_8 >= float(image_height_22);
        }
        if(_S4373)
        {
            _S4373 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4373 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4373 = false;
                }
                if(_S4373)
                {
                    _S4373 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4373 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4373 = false;
                    }
                }
            }
            else
            {
                _S4373 = false;
            }
        }
        if(_S4373)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4404 = mean_26 - - mul_0(transpose_0(R_26), t_25);
        float _S4405 = _S4404.x;
        float _S4406 = _S4404.y;
        float _S4407 = _S4404.z;
        float norm_15 = (F32_sqrt((_S4405 * _S4405 + _S4406 * _S4406 + _S4407 * _S4407)));
        float x_55 = _S4405 / norm_15;
        float y_25 = _S4406 / norm_15;
        float z_22 = _S4407 / norm_15;
        float z2_50 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_55 * x_55 - y_25 * y_25;
        float fS1_22 = 2.0f * x_55 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_50 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_55) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_50 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_55) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_55 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_50 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_55) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_55 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4408 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4408);
        float3  _S4409 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4410 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4409 + _S4410 + make_float3 (0.5f), _S4408);
        (*rgbs_1)[int(2)] = max_0(_S4409 - _S4410 + make_float3 (0.5f), _S4408);
        (*verts_1)[int(0)] = vert0_3;
        (*verts_1)[int(1)] = vert1_3;
        (*verts_1)[int(2)] = vert2_3;
        float3  _S4411 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4411 * make_float3 (float(- (F32_sign((dot_0(_S4411, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_32, float fy_32, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_37, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_27) + t_26;
    float _S4412 = scale_27.x;
    float sx_6 = (F32_exp((_S4412)));
    float _S4413 = scale_27.y;
    float sy_6 = (F32_exp((_S4413)));
    float sz_8 = scale_27.z - 0.5f * (_S4412 + _S4413);
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
    Matrix<float, 3, 3>  _S4414 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_51), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_51), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
    float3  vert0_4 = mul_0(_S4414, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
    float3  vert1_4 = mul_0(_S4414, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
    float3  vert2_4 = mul_0(_S4414, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
    float3  vert0_c_8 = mul_0(R_27, vert0_4) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_4) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_4) + t_26;
    float2  _S4415 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z);
    float u_68 = _S4415.x;
    float v_68 = _S4415.y;
    float r2_68 = u_68 * u_68 + v_68 * v_68;
    float _S4416 = 2.0f * (*dist_coeffs_37)[int(4)];
    float _S4417 = 2.0f * (*dist_coeffs_37)[int(5)];
    float2  _S4418 = _S4415 * make_float2 (1.0f + r2_68 * ((*dist_coeffs_37)[int(0)] + r2_68 * ((*dist_coeffs_37)[int(1)] + r2_68 * ((*dist_coeffs_37)[int(2)] + r2_68 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4416 * u_68 * v_68 + (*dist_coeffs_37)[int(5)] * (r2_68 + 2.0f * u_68 * u_68) + (*dist_coeffs_37)[int(6)] * r2_68, _S4417 * u_68 * v_68 + (*dist_coeffs_37)[int(4)] * (r2_68 + 2.0f * v_68 * v_68) + (*dist_coeffs_37)[int(7)] * r2_68);
    float2  _S4419 = _S4418 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4418.x + (*dist_coeffs_37)[int(9)] * _S4418.y, 0.0f);
    float _S4420 = fx_32 * _S4419.x + cx_28;
    float _S4421 = fy_32 * _S4419.y + cy_28;
    float2  uv0_8 = make_float2 (_S4420, _S4421);
    float2  _S4422 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z);
    float u_69 = _S4422.x;
    float v_69 = _S4422.y;
    float r2_69 = u_69 * u_69 + v_69 * v_69;
    float2  _S4423 = _S4422 * make_float2 (1.0f + r2_69 * ((*dist_coeffs_37)[int(0)] + r2_69 * ((*dist_coeffs_37)[int(1)] + r2_69 * ((*dist_coeffs_37)[int(2)] + r2_69 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4416 * u_69 * v_69 + (*dist_coeffs_37)[int(5)] * (r2_69 + 2.0f * u_69 * u_69) + (*dist_coeffs_37)[int(6)] * r2_69, _S4417 * u_69 * v_69 + (*dist_coeffs_37)[int(4)] * (r2_69 + 2.0f * v_69 * v_69) + (*dist_coeffs_37)[int(7)] * r2_69);
    float2  _S4424 = _S4423 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4423.x + (*dist_coeffs_37)[int(9)] * _S4423.y, 0.0f);
    float _S4425 = fx_32 * _S4424.x + cx_28;
    float _S4426 = fy_32 * _S4424.y + cy_28;
    float2  uv1_8 = make_float2 (_S4425, _S4426);
    float2  _S4427 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z);
    float u_70 = _S4427.x;
    float v_70 = _S4427.y;
    float r2_70 = u_70 * u_70 + v_70 * v_70;
    float2  _S4428 = _S4427 * make_float2 (1.0f + r2_70 * ((*dist_coeffs_37)[int(0)] + r2_70 * ((*dist_coeffs_37)[int(1)] + r2_70 * ((*dist_coeffs_37)[int(2)] + r2_70 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4416 * u_70 * v_70 + (*dist_coeffs_37)[int(5)] * (r2_70 + 2.0f * u_70 * u_70) + (*dist_coeffs_37)[int(6)] * r2_70, _S4417 * u_70 * v_70 + (*dist_coeffs_37)[int(4)] * (r2_70 + 2.0f * v_70 * v_70) + (*dist_coeffs_37)[int(7)] * r2_70);
    float2  _S4429 = _S4428 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4428.x + (*dist_coeffs_37)[int(9)] * _S4428.y, 0.0f);
    float _S4430 = fx_32 * _S4429.x + cx_28;
    float _S4431 = fy_32 * _S4429.y + cy_28;
    float2  uv2_8 = make_float2 (_S4430, _S4431);
    float2  e0_12 = uv1_8 - uv0_8;
    float2  e1_12 = uv2_8 - uv1_8;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(uv0_8 - uv2_8)));
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4420), (_S4425)))), (_S4430))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4421), (_S4426)))), (_S4431))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4420), (_S4425)))), (_S4430))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4421), (_S4426)))), (_S4431))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4432 = mean_27 - - mul_0(transpose_0(R_27), t_26);
    float _S4433 = _S4432.x;
    float _S4434 = _S4432.y;
    float _S4435 = _S4432.z;
    float norm_16 = (F32_sqrt((_S4433 * _S4433 + _S4434 * _S4434 + _S4435 * _S4435)));
    float x_57 = _S4433 / norm_16;
    float y_26 = _S4434 / norm_16;
    float z_23 = _S4435 / norm_16;
    float z2_52 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_57 * x_57 - y_26 * y_26;
    float fS1_23 = 2.0f * x_57 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_52 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_57) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_52 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_57) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_57 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_52 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_57) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_57 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S4436 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4436);
    float3  _S4437 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S4438 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S4437 + _S4438 + make_float3 (0.5f), _S4436);
    (*rgbs_2)[int(2)] = max_0(_S4437 - _S4438 + make_float3 (0.5f), _S4436);
    (*verts_2)[int(0)] = vert0_4;
    (*verts_2)[int(1)] = vert1_4;
    (*verts_2)[int(2)] = vert2_4;
    float3  _S4439 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S4439 * make_float3 (float(- (F32_sign((dot_0(_S4439, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_33, float fy_33, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_38, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_28) + t_27;
    float _S4440 = scale_28.x;
    float sx_7 = (F32_exp((_S4440)));
    float _S4441 = scale_28.y;
    float sy_7 = (F32_exp((_S4441)));
    float sz_9 = scale_28.z - 0.5f * (_S4440 + _S4441);
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
    Matrix<float, 3, 3>  _S4442 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_53), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_53), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4442, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4442, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4442, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_9 = mul_0(R_28, vert0_5) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_5) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_5) + t_27;
    float2  _S4443 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_33 = length_0(_S4443);
    float _S4444 = vert0_c_9.z;
    float theta_28 = (F32_atan2((r_33), (_S4444)));
    float k_11;
    if(theta_28 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_28 * theta_28 / 3.0f) / _S4444;
    }
    else
    {
        k_11 = theta_28 / r_33;
    }
    float2  _S4445 = _S4443 * make_float2 (k_11);
    float u_71 = _S4445.x;
    float v_71 = _S4445.y;
    float r2_71 = u_71 * u_71 + v_71 * v_71;
    float _S4446 = 2.0f * (*dist_coeffs_38)[int(4)];
    float _S4447 = 2.0f * (*dist_coeffs_38)[int(5)];
    float2  _S4448 = _S4445 * make_float2 (1.0f + r2_71 * ((*dist_coeffs_38)[int(0)] + r2_71 * ((*dist_coeffs_38)[int(1)] + r2_71 * ((*dist_coeffs_38)[int(2)] + r2_71 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4446 * u_71 * v_71 + (*dist_coeffs_38)[int(5)] * (r2_71 + 2.0f * u_71 * u_71) + (*dist_coeffs_38)[int(6)] * r2_71, _S4447 * u_71 * v_71 + (*dist_coeffs_38)[int(4)] * (r2_71 + 2.0f * v_71 * v_71) + (*dist_coeffs_38)[int(7)] * r2_71);
    float2  _S4449 = _S4448 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4448.x + (*dist_coeffs_38)[int(9)] * _S4448.y, 0.0f);
    float _S4450 = fx_33 * _S4449.x + cx_29;
    float _S4451 = fy_33 * _S4449.y + cy_29;
    float2  uv0_9 = make_float2 (_S4450, _S4451);
    float2  _S4452 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_34 = length_0(_S4452);
    float _S4453 = vert1_c_9.z;
    float theta_29 = (F32_atan2((r_34), (_S4453)));
    if(theta_29 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_29 * theta_29 / 3.0f) / _S4453;
    }
    else
    {
        k_11 = theta_29 / r_34;
    }
    float2  _S4454 = _S4452 * make_float2 (k_11);
    float u_72 = _S4454.x;
    float v_72 = _S4454.y;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float2  _S4455 = _S4454 * make_float2 (1.0f + r2_72 * ((*dist_coeffs_38)[int(0)] + r2_72 * ((*dist_coeffs_38)[int(1)] + r2_72 * ((*dist_coeffs_38)[int(2)] + r2_72 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4446 * u_72 * v_72 + (*dist_coeffs_38)[int(5)] * (r2_72 + 2.0f * u_72 * u_72) + (*dist_coeffs_38)[int(6)] * r2_72, _S4447 * u_72 * v_72 + (*dist_coeffs_38)[int(4)] * (r2_72 + 2.0f * v_72 * v_72) + (*dist_coeffs_38)[int(7)] * r2_72);
    float2  _S4456 = _S4455 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4455.x + (*dist_coeffs_38)[int(9)] * _S4455.y, 0.0f);
    float _S4457 = fx_33 * _S4456.x + cx_29;
    float _S4458 = fy_33 * _S4456.y + cy_29;
    float2  uv1_9 = make_float2 (_S4457, _S4458);
    float2  _S4459 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_35 = length_0(_S4459);
    float _S4460 = vert2_c_9.z;
    float theta_30 = (F32_atan2((r_35), (_S4460)));
    if(theta_30 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_30 * theta_30 / 3.0f) / _S4460;
    }
    else
    {
        k_11 = theta_30 / r_35;
    }
    float2  _S4461 = _S4459 * make_float2 (k_11);
    float u_73 = _S4461.x;
    float v_73 = _S4461.y;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float2  _S4462 = _S4461 * make_float2 (1.0f + r2_73 * ((*dist_coeffs_38)[int(0)] + r2_73 * ((*dist_coeffs_38)[int(1)] + r2_73 * ((*dist_coeffs_38)[int(2)] + r2_73 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4446 * u_73 * v_73 + (*dist_coeffs_38)[int(5)] * (r2_73 + 2.0f * u_73 * u_73) + (*dist_coeffs_38)[int(6)] * r2_73, _S4447 * u_73 * v_73 + (*dist_coeffs_38)[int(4)] * (r2_73 + 2.0f * v_73 * v_73) + (*dist_coeffs_38)[int(7)] * r2_73);
    float2  _S4463 = _S4462 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4462.x + (*dist_coeffs_38)[int(9)] * _S4462.y, 0.0f);
    float _S4464 = fx_33 * _S4463.x + cx_29;
    float _S4465 = fy_33 * _S4463.y + cy_29;
    float2  uv2_9 = make_float2 (_S4464, _S4465);
    float2  e0_13 = uv1_9 - uv0_9;
    float2  e1_13 = uv2_9 - uv1_9;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(uv0_9 - uv2_9)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4450), (_S4457)))), (_S4464))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4451), (_S4458)))), (_S4465))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4450), (_S4457)))), (_S4464))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4451), (_S4458)))), (_S4465))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4466 = mean_28 - - mul_0(transpose_0(R_28), t_27);
    float _S4467 = _S4466.x;
    float _S4468 = _S4466.y;
    float _S4469 = _S4466.z;
    float norm_17 = (F32_sqrt((_S4467 * _S4467 + _S4468 * _S4468 + _S4469 * _S4469)));
    float x_59 = _S4467 / norm_17;
    float y_27 = _S4468 / norm_17;
    float z_24 = _S4469 / norm_17;
    float z2_54 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_59 * x_59 - y_27 * y_27;
    float fS1_24 = 2.0f * x_59 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_54 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_59) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_54 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_59) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_59 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_54 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_59) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_59 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4470 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4470);
    float3  _S4471 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4472 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4471 + _S4472 + make_float3 (0.5f), _S4470);
    (*rgbs_3)[int(2)] = max_0(_S4471 - _S4472 + make_float3 (0.5f), _S4470);
    (*verts_3)[int(0)] = vert0_5;
    (*verts_3)[int(1)] = vert1_5;
    (*verts_3)[int(2)] = vert2_5;
    float3  _S4473 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4473 * make_float3 (float(- (F32_sign((dot_0(_S4473, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_34, float fy_34, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_39, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_29) + t_28;
    float _S4474 = scale_29.x;
    float _S4475 = s_primal_ctx_exp_1(_S4474);
    float _S4476 = scale_29.y;
    float _S4477 = s_primal_ctx_exp_1(_S4476);
    float sz_10 = scale_29.z - 0.5f * (_S4474 + _S4476);
    float _S4478 = quat_30.y;
    float x2_30 = _S4478 * _S4478;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_55 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S4479 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_55), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_55), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  _S4480 = make_float3 (_S4475, 0.0f, 0.0f);
    float3  vert0_6 = s_primal_ctx_mul_1(_S4479, _S4480) + mean_29;
    float _S4481 = -0.5f + sz_10;
    float3  _S4482 = make_float3 (_S4475 * _S4481, _S4477, 0.0f);
    float3  vert1_6 = s_primal_ctx_mul_1(_S4479, _S4482) + mean_29;
    float _S4483 = -0.5f - sz_10;
    float3  _S4484 = make_float3 (_S4475 * _S4483, - _S4477, 0.0f);
    float3  vert2_6 = s_primal_ctx_mul_1(_S4479, _S4484) + mean_29;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_6) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_6) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_6) + t_28;
    float2  _S4485 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4486 = vert0_c_10.z;
    float2  _S4487 = make_float2 (_S4486);
    float2  _S4488 = _S4485 / make_float2 (_S4486);
    float2  _S4489 = make_float2 (_S4486 * _S4486);
    float u_74 = _S4488.x;
    float v_74 = _S4488.y;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float _S4490 = (*dist_coeffs_39)[int(2)] + r2_74 * (*dist_coeffs_39)[int(3)];
    float _S4491 = (*dist_coeffs_39)[int(1)] + r2_74 * _S4490;
    float _S4492 = (*dist_coeffs_39)[int(0)] + r2_74 * _S4491;
    float radial_8 = 1.0f + r2_74 * _S4492;
    float2  _S4493 = make_float2 (radial_8);
    float _S4494 = 2.0f * (*dist_coeffs_39)[int(4)];
    float _S4495 = _S4494 * u_74;
    float _S4496 = 2.0f * u_74;
    float _S4497 = 2.0f * (*dist_coeffs_39)[int(5)];
    float _S4498 = _S4497 * u_74;
    float _S4499 = 2.0f * v_74;
    float2  _S4500 = _S4488 * make_float2 (radial_8) + make_float2 (_S4495 * v_74 + (*dist_coeffs_39)[int(5)] * (r2_74 + _S4496 * u_74) + (*dist_coeffs_39)[int(6)] * r2_74, _S4498 * v_74 + (*dist_coeffs_39)[int(4)] * (r2_74 + _S4499 * v_74) + (*dist_coeffs_39)[int(7)] * r2_74);
    float2  _S4501 = _S4500 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4500.x + (*dist_coeffs_39)[int(9)] * _S4500.y, 0.0f);
    float _S4502 = fx_34 * _S4501.x + cx_30;
    float _S4503 = fy_34 * _S4501.y + cy_30;
    float2  uv0_10 = make_float2 (_S4502, _S4503);
    float2  _S4504 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4505 = vert1_c_10.z;
    float2  _S4506 = make_float2 (_S4505);
    float2  _S4507 = _S4504 / make_float2 (_S4505);
    float2  _S4508 = make_float2 (_S4505 * _S4505);
    float u_75 = _S4507.x;
    float v_75 = _S4507.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S4509 = (*dist_coeffs_39)[int(2)] + r2_75 * (*dist_coeffs_39)[int(3)];
    float _S4510 = (*dist_coeffs_39)[int(1)] + r2_75 * _S4509;
    float _S4511 = (*dist_coeffs_39)[int(0)] + r2_75 * _S4510;
    float radial_9 = 1.0f + r2_75 * _S4511;
    float2  _S4512 = make_float2 (radial_9);
    float _S4513 = _S4494 * u_75;
    float _S4514 = 2.0f * u_75;
    float _S4515 = _S4497 * u_75;
    float _S4516 = 2.0f * v_75;
    float2  _S4517 = _S4507 * make_float2 (radial_9) + make_float2 (_S4513 * v_75 + (*dist_coeffs_39)[int(5)] * (r2_75 + _S4514 * u_75) + (*dist_coeffs_39)[int(6)] * r2_75, _S4515 * v_75 + (*dist_coeffs_39)[int(4)] * (r2_75 + _S4516 * v_75) + (*dist_coeffs_39)[int(7)] * r2_75);
    float2  _S4518 = _S4517 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4517.x + (*dist_coeffs_39)[int(9)] * _S4517.y, 0.0f);
    float _S4519 = fx_34 * _S4518.x + cx_30;
    float _S4520 = fy_34 * _S4518.y + cy_30;
    float2  uv1_10 = make_float2 (_S4519, _S4520);
    float2  _S4521 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4522 = vert2_c_10.z;
    float2  _S4523 = make_float2 (_S4522);
    float2  _S4524 = _S4521 / make_float2 (_S4522);
    float2  _S4525 = make_float2 (_S4522 * _S4522);
    float u_76 = _S4524.x;
    float v_76 = _S4524.y;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float _S4526 = (*dist_coeffs_39)[int(2)] + r2_76 * (*dist_coeffs_39)[int(3)];
    float _S4527 = (*dist_coeffs_39)[int(1)] + r2_76 * _S4526;
    float _S4528 = (*dist_coeffs_39)[int(0)] + r2_76 * _S4527;
    float radial_10 = 1.0f + r2_76 * _S4528;
    float2  _S4529 = make_float2 (radial_10);
    float _S4530 = _S4494 * u_76;
    float _S4531 = 2.0f * u_76;
    float _S4532 = _S4497 * u_76;
    float _S4533 = 2.0f * v_76;
    float2  _S4534 = _S4524 * make_float2 (radial_10) + make_float2 (_S4530 * v_76 + (*dist_coeffs_39)[int(5)] * (r2_76 + _S4531 * u_76) + (*dist_coeffs_39)[int(6)] * r2_76, _S4532 * v_76 + (*dist_coeffs_39)[int(4)] * (r2_76 + _S4533 * v_76) + (*dist_coeffs_39)[int(7)] * r2_76);
    float2  _S4535 = _S4534 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4534.x + (*dist_coeffs_39)[int(9)] * _S4534.y, 0.0f);
    float _S4536 = fx_34 * _S4535.x + cx_30;
    float _S4537 = fy_34 * _S4535.y + cy_30;
    float2  uv2_10 = make_float2 (_S4536, _S4537);
    float2  e0_14 = uv1_10 - uv0_10;
    float2  e1_14 = uv2_10 - uv1_10;
    float2  e2_6 = uv0_10 - uv2_10;
    float _S4538 = e0_14.x;
    float _S4539 = e1_14.y;
    float _S4540 = e0_14.y;
    float _S4541 = e1_14.x;
    float _S4542 = _S4538 * _S4539 - _S4540 * _S4541;
    float _S4543 = 1.0f - hardness_14.y;
    float _S4544 = -1.0f / _S4543;
    float _S4545 = _S4543 * _S4543;
    float _S4546 = s_primal_ctx_max_0(_S4502, _S4519);
    float _S4547 = s_primal_ctx_min_0(_S4502, _S4519);
    float _S4548 = s_primal_ctx_max_0(_S4503, _S4520);
    float _S4549 = s_primal_ctx_min_0(_S4503, _S4520);
    float3  _S4550 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4551 = length_1(_S4550) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4552 = transpose_0(R_29);
    float3  _S4553 = mean_29 - - s_primal_ctx_mul_1(_S4552, t_28);
    float _S4554 = _S4553.x;
    float _S4555 = _S4553.y;
    float _S4556 = _S4553.z;
    float _S4557 = _S4554 * _S4554 + _S4555 * _S4555 + _S4556 * _S4556;
    float _S4558 = s_primal_ctx_sqrt_0(_S4557);
    float x_60 = _S4554 / _S4558;
    float3  _S4559 = make_float3 (x_60);
    float _S4560 = _S4558 * _S4558;
    float y_28 = _S4555 / _S4558;
    float z_25 = _S4556 / _S4558;
    float3  _S4561 = make_float3 (z_25);
    float _S4562 = - y_28;
    float3  _S4563 = make_float3 (_S4562);
    float z2_56 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_60 * x_60 - y_28 * y_28;
    float _S4564 = 2.0f * x_60;
    float fS1_25 = _S4564 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_56 - 0.31539157032966614f;
    float3  _S4565 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_60;
    float3  _S4566 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4567 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4568 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4569 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_56 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4570 = 1.86588168144226074f * z2_56 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4570;
    float3  _S4571 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_60;
    float3  _S4572 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4573 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4574 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4575 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_60 * fC1_25 - y_28 * fS1_25);
    float3  _S4576 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_60 * fS1_25 + y_28 * fC1_25);
    float3  _S4577 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4562) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_60) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4578 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4579 = make_float3 (0.0f);
    float3  _S4580 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4581 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4582 = make_float3 (_S4581);
    float3  _S4583 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4581);
    float3  _S4584 = _S4580 + _S4583 + make_float3 (0.5f);
    float3  _S4585 = _S4580 - _S4583 + make_float3 (0.5f);
    float3  _S4586 = vert1_c_10 - vert0_c_10;
    float3  _S4587 = vert2_c_10 - vert0_c_10;
    float3  _S4588 = s_primal_ctx_cross_0(_S4586, _S4587);
    float3  _S4589 = normalize_0(_S4588);
    float3  _S4590 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4589, mean_c_25)))))) * v_normal_2;
    float3  _S4591 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4592;
    (&_S4592)->primal_0 = _S4589;
    (&_S4592)->differential_0 = _S4591;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4593;
    (&_S4593)->primal_0 = mean_c_25;
    (&_S4593)->differential_0 = _S4591;
    s_bwd_prop_dot_0(&_S4592, &_S4593, 0.0f);
    float3  _S4594 = _S4590 + _S4592.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4595;
    (&_S4595)->primal_0 = _S4588;
    (&_S4595)->differential_0 = _S4591;
    s_bwd_normalize_impl_0(&_S4595, _S4594);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4596;
    (&_S4596)->primal_0 = _S4586;
    (&_S4596)->differential_0 = _S4591;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4597;
    (&_S4597)->primal_0 = _S4587;
    (&_S4597)->differential_0 = _S4591;
    s_bwd_prop_cross_0(&_S4596, &_S4597, _S4595.differential_0);
    float3  _S4598 = - _S4597.differential_0;
    float3  _S4599 = - _S4596.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4600;
    (&_S4600)->primal_0 = _S4585;
    (&_S4600)->differential_0 = _S4591;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4601;
    (&_S4601)->primal_0 = _S4579;
    (&_S4601)->differential_0 = _S4591;
    s_bwd_prop_max_0(&_S4600, &_S4601, (*v_rgbs_0)[int(2)]);
    float3  _S4602 = - _S4600.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4603;
    (&_S4603)->primal_0 = _S4584;
    (&_S4603)->differential_0 = _S4591;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4604;
    (&_S4604)->primal_0 = _S4579;
    (&_S4604)->differential_0 = _S4591;
    s_bwd_prop_max_0(&_S4603, &_S4604, (*v_rgbs_0)[int(1)]);
    float3  _S4605 = _S4582 * (_S4602 + _S4603.differential_0);
    float3  _S4606 = _S4600.differential_0 + _S4603.differential_0;
    float3  _S4607 = make_float3 (0.5f) * - _S4606;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4608;
    (&_S4608)->primal_0 = _S4578;
    (&_S4608)->differential_0 = _S4591;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4609;
    (&_S4609)->primal_0 = _S4579;
    (&_S4609)->differential_0 = _S4591;
    s_bwd_prop_max_0(&_S4608, &_S4609, (*v_rgbs_0)[int(0)]);
    float3  _S4610 = _S4607 + _S4608.differential_0;
    float3  _S4611 = _S4606 + _S4608.differential_0;
    float3  _S4612 = _S4576 * _S4611;
    float3  _S4613 = (*sh_coeffs_25)[int(15)] * _S4611;
    float3  _S4614 = _S4574 * _S4611;
    float3  _S4615 = (*sh_coeffs_25)[int(14)] * _S4611;
    float3  _S4616 = _S4572 * _S4611;
    float3  _S4617 = (*sh_coeffs_25)[int(13)] * _S4611;
    float3  _S4618 = _S4571 * _S4611;
    float3  _S4619 = (*sh_coeffs_25)[int(12)] * _S4611;
    float3  _S4620 = _S4573 * _S4611;
    float3  _S4621 = (*sh_coeffs_25)[int(11)] * _S4611;
    float3  _S4622 = _S4575 * _S4611;
    float3  _S4623 = (*sh_coeffs_25)[int(10)] * _S4611;
    float3  _S4624 = _S4577 * _S4611;
    float3  _S4625 = (*sh_coeffs_25)[int(9)] * _S4611;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4625.x + _S4625.y + _S4625.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4613.x + _S4613.y + _S4613.z);
    float _S4626 = _S4623.x + _S4623.y + _S4623.z;
    float _S4627 = _S4615.x + _S4615.y + _S4615.z;
    float _S4628 = _S4621.x + _S4621.y + _S4621.z;
    float _S4629 = _S4617.x + _S4617.y + _S4617.z;
    float _S4630 = _S4619.x + _S4619.y + _S4619.z;
    float _S4631 = - s_diff_fC2_T_7;
    float3  _S4632 = _S4568 * _S4611;
    float3  _S4633 = (*sh_coeffs_25)[int(8)] * _S4611;
    float3  _S4634 = _S4566 * _S4611;
    float3  _S4635 = (*sh_coeffs_25)[int(7)] * _S4611;
    float3  _S4636 = _S4565 * _S4611;
    float3  _S4637 = (*sh_coeffs_25)[int(6)] * _S4611;
    float3  _S4638 = _S4567 * _S4611;
    float3  _S4639 = (*sh_coeffs_25)[int(5)] * _S4611;
    float3  _S4640 = _S4569 * _S4611;
    float3  _S4641 = (*sh_coeffs_25)[int(4)] * _S4611;
    float _S4642 = _S4639.x + _S4639.y + _S4639.z;
    float _S4643 = _S4635.x + _S4635.y + _S4635.z;
    float _S4644 = fTmp1B_25 * _S4626 + x_60 * s_diff_fS2_T_7 + y_28 * _S4631 + 0.54627424478530884f * (_S4641.x + _S4641.y + _S4641.z);
    float _S4645 = fTmp1B_25 * _S4627 + y_28 * s_diff_fS2_T_7 + x_60 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4633.x + _S4633.y + _S4633.z);
    float _S4646 = y_28 * - _S4645;
    float _S4647 = x_60 * _S4645;
    float _S4648 = z_25 * (1.86588168144226074f * (z_25 * _S4630) + -2.28522896766662598f * (y_28 * _S4628 + x_60 * _S4629) + 0.94617468118667603f * (_S4637.x + _S4637.y + _S4637.z));
    float3  _S4649 = make_float3 (0.48860251903533936f) * _S4611;
    float3  _S4650 = - _S4649;
    float3  _S4651 = _S4559 * _S4650;
    float3  _S4652 = (*sh_coeffs_25)[int(3)] * _S4650;
    float3  _S4653 = _S4561 * _S4649;
    float3  _S4654 = (*sh_coeffs_25)[int(2)] * _S4649;
    float3  _S4655 = _S4563 * _S4649;
    float3  _S4656 = (*sh_coeffs_25)[int(1)] * _S4649;
    float _S4657 = (_S4570 * _S4630 + 1.44530570507049561f * (fS1_25 * _S4626 + fC1_25 * _S4627) + -1.09254848957061768f * (y_28 * _S4642 + x_60 * _S4643) + _S4648 + _S4648 + _S4654.x + _S4654.y + _S4654.z) / _S4560;
    float _S4658 = _S4558 * _S4657;
    float _S4659 = (fTmp0C_25 * _S4628 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4631 + fTmp0B_25 * _S4642 + _S4564 * _S4644 + _S4646 + _S4646 + - (_S4656.x + _S4656.y + _S4656.z)) / _S4560;
    float _S4660 = _S4558 * _S4659;
    float _S4661 = (fTmp0C_25 * _S4629 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4643 + 2.0f * (y_28 * _S4644) + _S4647 + _S4647 + _S4652.x + _S4652.y + _S4652.z) / _S4560;
    float _S4662 = _S4558 * _S4661;
    float _S4663 = _S4556 * - _S4657 + _S4555 * - _S4659 + _S4554 * - _S4661;
    DiffPair_float_0 _S4664;
    (&_S4664)->primal_0 = _S4557;
    (&_S4664)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4664, _S4663);
    float _S4665 = _S4556 * _S4664.differential_0;
    float _S4666 = _S4555 * _S4664.differential_0;
    float _S4667 = _S4554 * _S4664.differential_0;
    float3  _S4668 = make_float3 (0.282094806432724f) * _S4611;
    float3  _S4669 = make_float3 (_S4662 + _S4667 + _S4667, _S4660 + _S4666 + _S4666, _S4658 + _S4665 + _S4665);
    float3  _S4670 = - - _S4669;
    Matrix<float, 3, 3>  _S4671 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4672;
    (&_S4672)->primal_0 = _S4552;
    (&_S4672)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4673;
    (&_S4673)->primal_0 = t_28;
    (&_S4673)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4672, &_S4673, _S4670);
    Matrix<float, 3, 3>  _S4674 = transpose_0(_S4672.differential_0);
    DiffPair_float_0 _S4675;
    (&_S4675)->primal_0 = _S4551;
    (&_S4675)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4675, v_depth_9);
    float _S4676 = 0.3333333432674408f * _S4675.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4677;
    (&_S4677)->primal_0 = _S4550;
    (&_S4677)->differential_0 = _S4591;
    s_bwd_length_impl_0(&_S4677, _S4676);
    DiffPair_float_0 _S4678;
    (&_S4678)->primal_0 = _S4549;
    (&_S4678)->differential_0 = 0.0f;
    DiffPair_float_0 _S4679;
    (&_S4679)->primal_0 = _S4537;
    (&_S4679)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4678, &_S4679, 0.0f);
    DiffPair_float_0 _S4680;
    (&_S4680)->primal_0 = _S4503;
    (&_S4680)->differential_0 = 0.0f;
    DiffPair_float_0 _S4681;
    (&_S4681)->primal_0 = _S4520;
    (&_S4681)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4680, &_S4681, _S4678.differential_0);
    DiffPair_float_0 _S4682;
    (&_S4682)->primal_0 = _S4548;
    (&_S4682)->differential_0 = 0.0f;
    DiffPair_float_0 _S4683;
    (&_S4683)->primal_0 = _S4537;
    (&_S4683)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4682, &_S4683, 0.0f);
    DiffPair_float_0 _S4684;
    (&_S4684)->primal_0 = _S4503;
    (&_S4684)->differential_0 = 0.0f;
    DiffPair_float_0 _S4685;
    (&_S4685)->primal_0 = _S4520;
    (&_S4685)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4684, &_S4685, _S4682.differential_0);
    DiffPair_float_0 _S4686;
    (&_S4686)->primal_0 = _S4547;
    (&_S4686)->differential_0 = 0.0f;
    DiffPair_float_0 _S4687;
    (&_S4687)->primal_0 = _S4536;
    (&_S4687)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4686, &_S4687, 0.0f);
    DiffPair_float_0 _S4688;
    (&_S4688)->primal_0 = _S4502;
    (&_S4688)->differential_0 = 0.0f;
    DiffPair_float_0 _S4689;
    (&_S4689)->primal_0 = _S4519;
    (&_S4689)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4688, &_S4689, _S4686.differential_0);
    DiffPair_float_0 _S4690;
    (&_S4690)->primal_0 = _S4546;
    (&_S4690)->differential_0 = 0.0f;
    DiffPair_float_0 _S4691;
    (&_S4691)->primal_0 = _S4536;
    (&_S4691)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4690, &_S4691, 0.0f);
    DiffPair_float_0 _S4692;
    (&_S4692)->primal_0 = _S4502;
    (&_S4692)->differential_0 = 0.0f;
    DiffPair_float_0 _S4693;
    (&_S4693)->primal_0 = _S4519;
    (&_S4693)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4692, &_S4693, _S4690.differential_0);
    DiffPair_float_0 _S4694;
    (&_S4694)->primal_0 = _S4544;
    (&_S4694)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4694, 0.0f);
    float _S4695 = - (-1.0f * - (_S4694.differential_0 / _S4545));
    float2  _S4696 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4697;
    (&_S4697)->primal_0 = e2_6;
    (&_S4697)->differential_0 = _S4696;
    s_bwd_length_impl_1(&_S4697, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4698;
    (&_S4698)->primal_0 = e1_14;
    (&_S4698)->differential_0 = _S4696;
    s_bwd_length_impl_1(&_S4698, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4699;
    (&_S4699)->primal_0 = e0_14;
    (&_S4699)->differential_0 = _S4696;
    s_bwd_length_impl_1(&_S4699, -0.0f);
    DiffPair_float_0 _S4700;
    (&_S4700)->primal_0 = _S4542;
    (&_S4700)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4700, 0.0f);
    float _S4701 = - _S4700.differential_0;
    float2  _S4702 = _S4698.differential_0 + make_float2 (_S4540 * _S4701, _S4538 * _S4700.differential_0);
    float2  _S4703 = _S4699.differential_0 + make_float2 (_S4539 * _S4700.differential_0, _S4541 * _S4701);
    float2  _S4704 = - _S4697.differential_0 + _S4702;
    float _S4705 = fx_34 * (_S4687.differential_0 + _S4691.differential_0 + _S4704.x);
    float2  _S4706 = make_float2 (_S4705, fy_34 * (_S4679.differential_0 + _S4683.differential_0 + _S4704.y)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S4705, (*dist_coeffs_39)[int(9)] * _S4705);
    float2  _S4707 = _S4524 * _S4706;
    float _S4708 = (*dist_coeffs_39)[int(4)] * _S4706.y;
    float _S4709 = (*dist_coeffs_39)[int(5)] * _S4706.x;
    float _S4710 = _S4707.x + _S4707.y;
    float _S4711 = r2_76 * _S4710;
    float _S4712 = r2_76 * _S4711;
    float _S4713 = (*dist_coeffs_39)[int(7)] * _S4706.y + _S4708 + (*dist_coeffs_39)[int(6)] * _S4706.x + _S4709 + _S4528 * _S4710 + _S4527 * _S4711 + _S4526 * _S4712 + (*dist_coeffs_39)[int(3)] * (r2_76 * _S4712);
    float _S4714 = v_76 * _S4713;
    float _S4715 = u_76 * _S4713;
    float2  _S4716 = (_S4529 * _S4706 + make_float2 (_S4497 * (v_76 * _S4706.y) + _S4531 * _S4709 + 2.0f * (u_76 * _S4709) + _S4494 * (v_76 * _S4706.x) + _S4715 + _S4715, _S4533 * _S4708 + 2.0f * (v_76 * _S4708) + _S4532 * _S4706.y + _S4530 * _S4706.x + _S4714 + _S4714)) / _S4525;
    float2  _S4717 = _S4521 * - _S4716;
    float2  _S4718 = _S4523 * _S4716;
    float2  _S4719 = - _S4702 + _S4703;
    float _S4720 = fx_34 * (_S4689.differential_0 + _S4693.differential_0 + _S4719.x);
    float2  _S4721 = make_float2 (_S4720, fy_34 * (_S4681.differential_0 + _S4685.differential_0 + _S4719.y)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S4720, (*dist_coeffs_39)[int(9)] * _S4720);
    float2  _S4722 = _S4507 * _S4721;
    float _S4723 = (*dist_coeffs_39)[int(4)] * _S4721.y;
    float _S4724 = (*dist_coeffs_39)[int(5)] * _S4721.x;
    float _S4725 = _S4722.x + _S4722.y;
    float _S4726 = r2_75 * _S4725;
    float _S4727 = r2_75 * _S4726;
    float _S4728 = (*dist_coeffs_39)[int(7)] * _S4721.y + _S4723 + (*dist_coeffs_39)[int(6)] * _S4721.x + _S4724 + _S4511 * _S4725 + _S4510 * _S4726 + _S4509 * _S4727 + (*dist_coeffs_39)[int(3)] * (r2_75 * _S4727);
    float _S4729 = v_75 * _S4728;
    float _S4730 = u_75 * _S4728;
    float2  _S4731 = (_S4512 * _S4721 + make_float2 (_S4497 * (v_75 * _S4721.y) + _S4514 * _S4724 + 2.0f * (u_75 * _S4724) + _S4494 * (v_75 * _S4721.x) + _S4730 + _S4730, _S4516 * _S4723 + 2.0f * (v_75 * _S4723) + _S4515 * _S4721.y + _S4513 * _S4721.x + _S4729 + _S4729)) / _S4508;
    float2  _S4732 = _S4504 * - _S4731;
    float2  _S4733 = _S4506 * _S4731;
    float _S4734 = _S4732.x + _S4732.y;
    float2  _S4735 = _S4697.differential_0 + - _S4703;
    float _S4736 = fx_34 * (_S4688.differential_0 + _S4692.differential_0 + _S4735.x);
    float2  _S4737 = make_float2 (_S4736, fy_34 * (_S4680.differential_0 + _S4684.differential_0 + _S4735.y)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S4736, (*dist_coeffs_39)[int(9)] * _S4736);
    float2  _S4738 = _S4488 * _S4737;
    float _S4739 = (*dist_coeffs_39)[int(4)] * _S4737.y;
    float _S4740 = (*dist_coeffs_39)[int(5)] * _S4737.x;
    float _S4741 = _S4738.x + _S4738.y;
    float _S4742 = r2_74 * _S4741;
    float _S4743 = r2_74 * _S4742;
    float _S4744 = (*dist_coeffs_39)[int(7)] * _S4737.y + _S4739 + (*dist_coeffs_39)[int(6)] * _S4737.x + _S4740 + _S4492 * _S4741 + _S4491 * _S4742 + _S4490 * _S4743 + (*dist_coeffs_39)[int(3)] * (r2_74 * _S4743);
    float _S4745 = v_74 * _S4744;
    float _S4746 = u_74 * _S4744;
    float2  _S4747 = (_S4493 * _S4737 + make_float2 (_S4497 * (v_74 * _S4737.y) + _S4496 * _S4740 + 2.0f * (u_74 * _S4740) + _S4494 * (v_74 * _S4737.x) + _S4746 + _S4746, _S4499 * _S4739 + 2.0f * (v_74 * _S4739) + _S4498 * _S4737.y + _S4495 * _S4737.x + _S4745 + _S4745)) / _S4489;
    float2  _S4748 = _S4485 * - _S4747;
    float2  _S4749 = _S4487 * _S4747;
    float _S4750 = _S4748.x + _S4748.y;
    float3  _S4751 = _S4597.differential_0 + _S4677.differential_0 + make_float3 (_S4718.x, _S4718.y, _S4717.x + _S4717.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4752;
    (&_S4752)->primal_0 = R_29;
    (&_S4752)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4753;
    (&_S4753)->primal_0 = vert2_6;
    (&_S4753)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4752, &_S4753, _S4751);
    float3  _S4754 = _S4596.differential_0 + _S4677.differential_0 + make_float3 (_S4733.x, _S4733.y, _S4734);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4755;
    (&_S4755)->primal_0 = R_29;
    (&_S4755)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4756;
    (&_S4756)->primal_0 = vert1_6;
    (&_S4756)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4755, &_S4756, _S4754);
    float3  _S4757 = _S4598 + _S4599 + _S4677.differential_0 + make_float3 (_S4749.x, _S4749.y, _S4750);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4758;
    (&_S4758)->primal_0 = R_29;
    (&_S4758)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4759;
    (&_S4759)->primal_0 = vert0_6;
    (&_S4759)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4758, &_S4759, _S4757);
    float3  _S4760 = (*v_verts_0)[int(2)] + _S4753.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4761;
    (&_S4761)->primal_0 = _S4479;
    (&_S4761)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4762;
    (&_S4762)->primal_0 = _S4484;
    (&_S4762)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4761, &_S4762, _S4760);
    float _S4763 = - _S4762.differential_0.y;
    float _S4764 = _S4483 * _S4762.differential_0.x;
    float _S4765 = - (_S4475 * _S4762.differential_0.x);
    float3  _S4766 = (*v_verts_0)[int(1)] + _S4756.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4767;
    (&_S4767)->primal_0 = _S4479;
    (&_S4767)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4768;
    (&_S4768)->primal_0 = _S4482;
    (&_S4768)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4767, &_S4768, _S4766);
    float _S4769 = _S4475 * _S4768.differential_0.x;
    float _S4770 = _S4481 * _S4768.differential_0.x;
    float3  _S4771 = (*v_verts_0)[int(0)] + _S4759.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4772;
    (&_S4772)->primal_0 = _S4479;
    (&_S4772)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4773;
    (&_S4773)->primal_0 = _S4480;
    (&_S4773)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4772, &_S4773, _S4771);
    Matrix<float, 3, 3>  _S4774 = transpose_0(_S4761.differential_0 + _S4767.differential_0 + _S4772.differential_0);
    float _S4775 = 2.0f * - _S4774.rows[int(2)].z;
    float _S4776 = 2.0f * _S4774.rows[int(2)].y;
    float _S4777 = 2.0f * _S4774.rows[int(2)].x;
    float _S4778 = 2.0f * _S4774.rows[int(1)].z;
    float _S4779 = 2.0f * - _S4774.rows[int(1)].y;
    float _S4780 = 2.0f * _S4774.rows[int(1)].x;
    float _S4781 = 2.0f * _S4774.rows[int(0)].z;
    float _S4782 = 2.0f * _S4774.rows[int(0)].y;
    float _S4783 = 2.0f * - _S4774.rows[int(0)].x;
    float _S4784 = - _S4780 + _S4782;
    float _S4785 = _S4777 + - _S4781;
    float _S4786 = - _S4776 + _S4778;
    float _S4787 = _S4776 + _S4778;
    float _S4788 = _S4777 + _S4781;
    float _S4789 = _S4780 + _S4782;
    float _S4790 = quat_30.w * (_S4779 + _S4783);
    float _S4791 = quat_30.z * (_S4775 + _S4783);
    float _S4792 = quat_30.y * (_S4775 + _S4779);
    float _S4793 = quat_30.x * _S4784 + quat_30.z * _S4787 + quat_30.y * _S4788 + _S4790 + _S4790;
    float _S4794 = quat_30.x * _S4785 + quat_30.w * _S4787 + quat_30.y * _S4789 + _S4791 + _S4791;
    float _S4795 = quat_30.x * _S4786 + quat_30.w * _S4788 + quat_30.z * _S4789 + _S4792 + _S4792;
    float _S4796 = quat_30.w * _S4784 + quat_30.z * _S4785 + quat_30.y * _S4786;
    float _S4797 = _S4765 + _S4769;
    float _S4798 = 0.5f * - _S4797;
    float _S4799 = _S4763 + _S4768.differential_0.y;
    DiffPair_float_0 _S4800;
    (&_S4800)->primal_0 = _S4476;
    (&_S4800)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4800, _S4799);
    float _S4801 = _S4798 + _S4800.differential_0;
    float _S4802 = _S4764 + _S4770 + _S4773.differential_0.x;
    DiffPair_float_0 _S4803;
    (&_S4803)->primal_0 = _S4474;
    (&_S4803)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4803, _S4802);
    float _S4804 = _S4798 + _S4803.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4805;
    (&_S4805)->primal_0 = R_29;
    (&_S4805)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4806;
    (&_S4806)->primal_0 = mean_29;
    (&_S4806)->differential_0 = _S4591;
    s_bwd_prop_mul_1(&_S4805, &_S4806, _S4593.differential_0);
    float3  _S4807 = _S4673.differential_0 + _S4751 + _S4754 + _S4757 + _S4593.differential_0;
    Matrix<float, 3, 3>  _S4808 = _S4674 + _S4752.differential_0 + _S4755.differential_0 + _S4758.differential_0 + _S4805.differential_0;
    FixedArray<float3 , 2>  _S4809;
    _S4809[int(0)] = _S4591;
    _S4809[int(1)] = _S4591;
    _S4809[int(1)] = _S4605;
    _S4809[int(0)] = _S4610;
    FixedArray<float3 , 16>  _S4810;
    _S4810[int(0)] = _S4591;
    _S4810[int(1)] = _S4591;
    _S4810[int(2)] = _S4591;
    _S4810[int(3)] = _S4591;
    _S4810[int(4)] = _S4591;
    _S4810[int(5)] = _S4591;
    _S4810[int(6)] = _S4591;
    _S4810[int(7)] = _S4591;
    _S4810[int(8)] = _S4591;
    _S4810[int(9)] = _S4591;
    _S4810[int(10)] = _S4591;
    _S4810[int(11)] = _S4591;
    _S4810[int(12)] = _S4591;
    _S4810[int(13)] = _S4591;
    _S4810[int(14)] = _S4591;
    _S4810[int(15)] = _S4591;
    _S4810[int(15)] = _S4612;
    _S4810[int(14)] = _S4614;
    _S4810[int(13)] = _S4616;
    _S4810[int(12)] = _S4618;
    _S4810[int(11)] = _S4620;
    _S4810[int(10)] = _S4622;
    _S4810[int(9)] = _S4624;
    _S4810[int(8)] = _S4632;
    _S4810[int(7)] = _S4634;
    _S4810[int(6)] = _S4636;
    _S4810[int(5)] = _S4638;
    _S4810[int(4)] = _S4640;
    _S4810[int(3)] = _S4651;
    _S4810[int(2)] = _S4653;
    _S4810[int(1)] = _S4655;
    _S4810[int(0)] = _S4668;
    float2  _S4811 = make_float2 (0.0f, _S4695);
    float3  _S4812 = make_float3 (_S4804, _S4801, _S4797);
    float4  _S4813 = make_float4 (0.0f);
    *&((&_S4813)->w) = _S4793;
    *&((&_S4813)->z) = _S4794;
    *&((&_S4813)->y) = _S4795;
    *&((&_S4813)->x) = _S4796;
    *v_mean_9 = _S4669 + _S4760 + _S4766 + _S4771 + _S4806.differential_0;
    *v_quat_8 = _S4813;
    *v_scale_8 = _S4812;
    *v_hardness_4 = _S4811;
    *v_sh_coeffs_7 = _S4810;
    *v_ch_coeffs_2 = _S4809;
    *v_R_8 = _S4808;
    *v_t_8 = _S4807;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_35, float fy_35, float cx_31, float cy_31, FixedArray<float, 10>  * dist_coeffs_40, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_30) + t_29;
    float _S4814 = scale_30.x;
    float _S4815 = s_primal_ctx_exp_1(_S4814);
    float _S4816 = scale_30.y;
    float _S4817 = s_primal_ctx_exp_1(_S4816);
    float sz_11 = scale_30.z - 0.5f * (_S4814 + _S4816);
    float _S4818 = quat_31.y;
    float x2_31 = _S4818 * _S4818;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_57 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4819 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_57), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_57), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4820 = make_float3 (_S4815, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4819, _S4820) + mean_30;
    float _S4821 = -0.5f + sz_11;
    float3  _S4822 = make_float3 (_S4815 * _S4821, _S4817, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4819, _S4822) + mean_30;
    float _S4823 = -0.5f - sz_11;
    float3  _S4824 = make_float3 (_S4815 * _S4823, - _S4817, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4819, _S4824) + mean_30;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_7) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_7) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_7) + t_29;
    float2  _S4825 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4826 = length_0(_S4825);
    float _S4827 = vert0_c_11.z;
    float _S4828 = s_primal_ctx_atan2_0(_S4826, _S4827);
    bool _S4829 = _S4828 < 0.00100000004749745f;
    float k_12;
    float _S4830;
    float _S4831;
    float _S4832;
    if(_S4829)
    {
        float _S4833 = 1.0f - _S4828 * _S4828 / 3.0f;
        float _S4834 = _S4827 * _S4827;
        k_12 = _S4833 / _S4827;
        _S4830 = 0.0f;
        _S4831 = _S4834;
        _S4832 = _S4833;
    }
    else
    {
        float _S4835 = _S4826 * _S4826;
        k_12 = _S4828 / _S4826;
        _S4830 = _S4835;
        _S4831 = 0.0f;
        _S4832 = 0.0f;
    }
    float2  _S4836 = make_float2 (k_12);
    float2  _S4837 = _S4825 * make_float2 (k_12);
    float u_77 = _S4837.x;
    float v_77 = _S4837.y;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float _S4838 = (*dist_coeffs_40)[int(2)] + r2_77 * (*dist_coeffs_40)[int(3)];
    float _S4839 = (*dist_coeffs_40)[int(1)] + r2_77 * _S4838;
    float _S4840 = (*dist_coeffs_40)[int(0)] + r2_77 * _S4839;
    float radial_11 = 1.0f + r2_77 * _S4840;
    float2  _S4841 = make_float2 (radial_11);
    float _S4842 = 2.0f * (*dist_coeffs_40)[int(4)];
    float _S4843 = _S4842 * u_77;
    float _S4844 = 2.0f * u_77;
    float _S4845 = 2.0f * (*dist_coeffs_40)[int(5)];
    float _S4846 = _S4845 * u_77;
    float _S4847 = 2.0f * v_77;
    float2  _S4848 = _S4837 * make_float2 (radial_11) + make_float2 (_S4843 * v_77 + (*dist_coeffs_40)[int(5)] * (r2_77 + _S4844 * u_77) + (*dist_coeffs_40)[int(6)] * r2_77, _S4846 * v_77 + (*dist_coeffs_40)[int(4)] * (r2_77 + _S4847 * v_77) + (*dist_coeffs_40)[int(7)] * r2_77);
    float2  _S4849 = _S4848 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4848.x + (*dist_coeffs_40)[int(9)] * _S4848.y, 0.0f);
    float _S4850 = fx_35 * _S4849.x + cx_31;
    float _S4851 = fy_35 * _S4849.y + cy_31;
    float2  uv0_11 = make_float2 (_S4850, _S4851);
    float2  _S4852 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4853 = length_0(_S4852);
    float _S4854 = vert1_c_11.z;
    float _S4855 = s_primal_ctx_atan2_0(_S4853, _S4854);
    bool _S4856 = _S4855 < 0.00100000004749745f;
    float _S4857;
    float _S4858;
    float _S4859;
    if(_S4856)
    {
        float _S4860 = 1.0f - _S4855 * _S4855 / 3.0f;
        float _S4861 = _S4854 * _S4854;
        k_12 = _S4860 / _S4854;
        _S4857 = 0.0f;
        _S4858 = _S4861;
        _S4859 = _S4860;
    }
    else
    {
        float _S4862 = _S4853 * _S4853;
        k_12 = _S4855 / _S4853;
        _S4857 = _S4862;
        _S4858 = 0.0f;
        _S4859 = 0.0f;
    }
    float2  _S4863 = make_float2 (k_12);
    float2  _S4864 = _S4852 * make_float2 (k_12);
    float u_78 = _S4864.x;
    float v_78 = _S4864.y;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float _S4865 = (*dist_coeffs_40)[int(2)] + r2_78 * (*dist_coeffs_40)[int(3)];
    float _S4866 = (*dist_coeffs_40)[int(1)] + r2_78 * _S4865;
    float _S4867 = (*dist_coeffs_40)[int(0)] + r2_78 * _S4866;
    float radial_12 = 1.0f + r2_78 * _S4867;
    float2  _S4868 = make_float2 (radial_12);
    float _S4869 = _S4842 * u_78;
    float _S4870 = 2.0f * u_78;
    float _S4871 = _S4845 * u_78;
    float _S4872 = 2.0f * v_78;
    float2  _S4873 = _S4864 * make_float2 (radial_12) + make_float2 (_S4869 * v_78 + (*dist_coeffs_40)[int(5)] * (r2_78 + _S4870 * u_78) + (*dist_coeffs_40)[int(6)] * r2_78, _S4871 * v_78 + (*dist_coeffs_40)[int(4)] * (r2_78 + _S4872 * v_78) + (*dist_coeffs_40)[int(7)] * r2_78);
    float2  _S4874 = _S4873 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4873.x + (*dist_coeffs_40)[int(9)] * _S4873.y, 0.0f);
    float _S4875 = fx_35 * _S4874.x + cx_31;
    float _S4876 = fy_35 * _S4874.y + cy_31;
    float2  uv1_11 = make_float2 (_S4875, _S4876);
    float2  _S4877 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4878 = length_0(_S4877);
    float _S4879 = vert2_c_11.z;
    float _S4880 = s_primal_ctx_atan2_0(_S4878, _S4879);
    bool _S4881 = _S4880 < 0.00100000004749745f;
    float _S4882;
    float _S4883;
    float _S4884;
    if(_S4881)
    {
        float _S4885 = 1.0f - _S4880 * _S4880 / 3.0f;
        float _S4886 = _S4879 * _S4879;
        k_12 = _S4885 / _S4879;
        _S4882 = 0.0f;
        _S4883 = _S4886;
        _S4884 = _S4885;
    }
    else
    {
        float _S4887 = _S4878 * _S4878;
        k_12 = _S4880 / _S4878;
        _S4882 = _S4887;
        _S4883 = 0.0f;
        _S4884 = 0.0f;
    }
    float2  _S4888 = make_float2 (k_12);
    float2  _S4889 = _S4877 * make_float2 (k_12);
    float u_79 = _S4889.x;
    float v_79 = _S4889.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S4890 = (*dist_coeffs_40)[int(2)] + r2_79 * (*dist_coeffs_40)[int(3)];
    float _S4891 = (*dist_coeffs_40)[int(1)] + r2_79 * _S4890;
    float _S4892 = (*dist_coeffs_40)[int(0)] + r2_79 * _S4891;
    float radial_13 = 1.0f + r2_79 * _S4892;
    float2  _S4893 = make_float2 (radial_13);
    float _S4894 = _S4842 * u_79;
    float _S4895 = 2.0f * u_79;
    float _S4896 = _S4845 * u_79;
    float _S4897 = 2.0f * v_79;
    float2  _S4898 = _S4889 * make_float2 (radial_13) + make_float2 (_S4894 * v_79 + (*dist_coeffs_40)[int(5)] * (r2_79 + _S4895 * u_79) + (*dist_coeffs_40)[int(6)] * r2_79, _S4896 * v_79 + (*dist_coeffs_40)[int(4)] * (r2_79 + _S4897 * v_79) + (*dist_coeffs_40)[int(7)] * r2_79);
    float2  _S4899 = _S4898 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4898.x + (*dist_coeffs_40)[int(9)] * _S4898.y, 0.0f);
    float _S4900 = fx_35 * _S4899.x + cx_31;
    float _S4901 = fy_35 * _S4899.y + cy_31;
    float2  uv2_11 = make_float2 (_S4900, _S4901);
    float2  e0_15 = uv1_11 - uv0_11;
    float2  e1_15 = uv2_11 - uv1_11;
    float2  e2_7 = uv0_11 - uv2_11;
    float _S4902 = e0_15.x;
    float _S4903 = e1_15.y;
    float _S4904 = e0_15.y;
    float _S4905 = e1_15.x;
    float _S4906 = _S4902 * _S4903 - _S4904 * _S4905;
    float _S4907 = 1.0f - hardness_15.y;
    float _S4908 = -1.0f / _S4907;
    float _S4909 = _S4907 * _S4907;
    float _S4910 = s_primal_ctx_max_0(_S4850, _S4875);
    float _S4911 = s_primal_ctx_min_0(_S4850, _S4875);
    float _S4912 = s_primal_ctx_max_0(_S4851, _S4876);
    float _S4913 = s_primal_ctx_min_0(_S4851, _S4876);
    float3  _S4914 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4915 = length_1(_S4914) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4916 = transpose_0(R_30);
    float3  _S4917 = mean_30 - - s_primal_ctx_mul_1(_S4916, t_29);
    float _S4918 = _S4917.x;
    float _S4919 = _S4917.y;
    float _S4920 = _S4917.z;
    float _S4921 = _S4918 * _S4918 + _S4919 * _S4919 + _S4920 * _S4920;
    float _S4922 = s_primal_ctx_sqrt_0(_S4921);
    float x_61 = _S4918 / _S4922;
    float3  _S4923 = make_float3 (x_61);
    float _S4924 = _S4922 * _S4922;
    float y_29 = _S4919 / _S4922;
    float z_26 = _S4920 / _S4922;
    float3  _S4925 = make_float3 (z_26);
    float _S4926 = - y_29;
    float3  _S4927 = make_float3 (_S4926);
    float z2_58 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_61 * x_61 - y_29 * y_29;
    float _S4928 = 2.0f * x_61;
    float fS1_26 = _S4928 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_58 - 0.31539157032966614f;
    float3  _S4929 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_61;
    float3  _S4930 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4931 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4932 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4933 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_58 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4934 = 1.86588168144226074f * z2_58 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4934;
    float3  _S4935 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_61;
    float3  _S4936 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4937 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4938 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4939 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_61 * fC1_26 - y_29 * fS1_26);
    float3  _S4940 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_61 * fS1_26 + y_29 * fC1_26);
    float3  _S4941 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4926) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_61) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4942 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4943 = make_float3 (0.0f);
    float3  _S4944 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4945 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4946 = make_float3 (_S4945);
    float3  _S4947 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4945);
    float3  _S4948 = _S4944 + _S4947 + make_float3 (0.5f);
    float3  _S4949 = _S4944 - _S4947 + make_float3 (0.5f);
    float3  _S4950 = vert1_c_11 - vert0_c_11;
    float3  _S4951 = vert2_c_11 - vert0_c_11;
    float3  _S4952 = s_primal_ctx_cross_0(_S4950, _S4951);
    float3  _S4953 = normalize_0(_S4952);
    float3  _S4954 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4953, mean_c_26)))))) * v_normal_3;
    float3  _S4955 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4956;
    (&_S4956)->primal_0 = _S4953;
    (&_S4956)->differential_0 = _S4955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4957;
    (&_S4957)->primal_0 = mean_c_26;
    (&_S4957)->differential_0 = _S4955;
    s_bwd_prop_dot_0(&_S4956, &_S4957, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4958 = _S4957;
    float3  _S4959 = _S4954 + _S4956.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4960;
    (&_S4960)->primal_0 = _S4952;
    (&_S4960)->differential_0 = _S4955;
    s_bwd_normalize_impl_0(&_S4960, _S4959);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4961;
    (&_S4961)->primal_0 = _S4950;
    (&_S4961)->differential_0 = _S4955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4962;
    (&_S4962)->primal_0 = _S4951;
    (&_S4962)->differential_0 = _S4955;
    s_bwd_prop_cross_0(&_S4961, &_S4962, _S4960.differential_0);
    float3  _S4963 = - _S4962.differential_0;
    float3  _S4964 = - _S4961.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4965;
    (&_S4965)->primal_0 = _S4949;
    (&_S4965)->differential_0 = _S4955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4966;
    (&_S4966)->primal_0 = _S4943;
    (&_S4966)->differential_0 = _S4955;
    s_bwd_prop_max_0(&_S4965, &_S4966, (*v_rgbs_1)[int(2)]);
    float3  _S4967 = - _S4965.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4968;
    (&_S4968)->primal_0 = _S4948;
    (&_S4968)->differential_0 = _S4955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4969;
    (&_S4969)->primal_0 = _S4943;
    (&_S4969)->differential_0 = _S4955;
    s_bwd_prop_max_0(&_S4968, &_S4969, (*v_rgbs_1)[int(1)]);
    float3  _S4970 = _S4946 * (_S4967 + _S4968.differential_0);
    float3  _S4971 = _S4965.differential_0 + _S4968.differential_0;
    float3  _S4972 = make_float3 (0.5f) * - _S4971;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4973;
    (&_S4973)->primal_0 = _S4942;
    (&_S4973)->differential_0 = _S4955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4974;
    (&_S4974)->primal_0 = _S4943;
    (&_S4974)->differential_0 = _S4955;
    s_bwd_prop_max_0(&_S4973, &_S4974, (*v_rgbs_1)[int(0)]);
    float3  _S4975 = _S4972 + _S4973.differential_0;
    float3  _S4976 = _S4971 + _S4973.differential_0;
    float3  _S4977 = _S4940 * _S4976;
    float3  _S4978 = (*sh_coeffs_26)[int(15)] * _S4976;
    float3  _S4979 = _S4938 * _S4976;
    float3  _S4980 = (*sh_coeffs_26)[int(14)] * _S4976;
    float3  _S4981 = _S4936 * _S4976;
    float3  _S4982 = (*sh_coeffs_26)[int(13)] * _S4976;
    float3  _S4983 = _S4935 * _S4976;
    float3  _S4984 = (*sh_coeffs_26)[int(12)] * _S4976;
    float3  _S4985 = _S4937 * _S4976;
    float3  _S4986 = (*sh_coeffs_26)[int(11)] * _S4976;
    float3  _S4987 = _S4939 * _S4976;
    float3  _S4988 = (*sh_coeffs_26)[int(10)] * _S4976;
    float3  _S4989 = _S4941 * _S4976;
    float3  _S4990 = (*sh_coeffs_26)[int(9)] * _S4976;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4990.x + _S4990.y + _S4990.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4978.x + _S4978.y + _S4978.z);
    float _S4991 = _S4988.x + _S4988.y + _S4988.z;
    float _S4992 = _S4980.x + _S4980.y + _S4980.z;
    float _S4993 = _S4986.x + _S4986.y + _S4986.z;
    float _S4994 = _S4982.x + _S4982.y + _S4982.z;
    float _S4995 = _S4984.x + _S4984.y + _S4984.z;
    float _S4996 = - s_diff_fC2_T_8;
    float3  _S4997 = _S4932 * _S4976;
    float3  _S4998 = (*sh_coeffs_26)[int(8)] * _S4976;
    float3  _S4999 = _S4930 * _S4976;
    float3  _S5000 = (*sh_coeffs_26)[int(7)] * _S4976;
    float3  _S5001 = _S4929 * _S4976;
    float3  _S5002 = (*sh_coeffs_26)[int(6)] * _S4976;
    float3  _S5003 = _S4931 * _S4976;
    float3  _S5004 = (*sh_coeffs_26)[int(5)] * _S4976;
    float3  _S5005 = _S4933 * _S4976;
    float3  _S5006 = (*sh_coeffs_26)[int(4)] * _S4976;
    float _S5007 = _S5004.x + _S5004.y + _S5004.z;
    float _S5008 = _S5000.x + _S5000.y + _S5000.z;
    float _S5009 = fTmp1B_26 * _S4991 + x_61 * s_diff_fS2_T_8 + y_29 * _S4996 + 0.54627424478530884f * (_S5006.x + _S5006.y + _S5006.z);
    float _S5010 = fTmp1B_26 * _S4992 + y_29 * s_diff_fS2_T_8 + x_61 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4998.x + _S4998.y + _S4998.z);
    float _S5011 = y_29 * - _S5010;
    float _S5012 = x_61 * _S5010;
    float _S5013 = z_26 * (1.86588168144226074f * (z_26 * _S4995) + -2.28522896766662598f * (y_29 * _S4993 + x_61 * _S4994) + 0.94617468118667603f * (_S5002.x + _S5002.y + _S5002.z));
    float3  _S5014 = make_float3 (0.48860251903533936f) * _S4976;
    float3  _S5015 = - _S5014;
    float3  _S5016 = _S4923 * _S5015;
    float3  _S5017 = (*sh_coeffs_26)[int(3)] * _S5015;
    float3  _S5018 = _S4925 * _S5014;
    float3  _S5019 = (*sh_coeffs_26)[int(2)] * _S5014;
    float3  _S5020 = _S4927 * _S5014;
    float3  _S5021 = (*sh_coeffs_26)[int(1)] * _S5014;
    float _S5022 = (_S4934 * _S4995 + 1.44530570507049561f * (fS1_26 * _S4991 + fC1_26 * _S4992) + -1.09254848957061768f * (y_29 * _S5007 + x_61 * _S5008) + _S5013 + _S5013 + _S5019.x + _S5019.y + _S5019.z) / _S4924;
    float _S5023 = _S4922 * _S5022;
    float _S5024 = (fTmp0C_26 * _S4993 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4996 + fTmp0B_26 * _S5007 + _S4928 * _S5009 + _S5011 + _S5011 + - (_S5021.x + _S5021.y + _S5021.z)) / _S4924;
    float _S5025 = _S4922 * _S5024;
    float _S5026 = (fTmp0C_26 * _S4994 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S5008 + 2.0f * (y_29 * _S5009) + _S5012 + _S5012 + _S5017.x + _S5017.y + _S5017.z) / _S4924;
    float _S5027 = _S4922 * _S5026;
    float _S5028 = _S4920 * - _S5022 + _S4919 * - _S5024 + _S4918 * - _S5026;
    DiffPair_float_0 _S5029;
    (&_S5029)->primal_0 = _S4921;
    (&_S5029)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5029, _S5028);
    float _S5030 = _S4920 * _S5029.differential_0;
    float _S5031 = _S4919 * _S5029.differential_0;
    float _S5032 = _S4918 * _S5029.differential_0;
    float3  _S5033 = make_float3 (0.282094806432724f) * _S4976;
    float3  _S5034 = make_float3 (_S5027 + _S5032 + _S5032, _S5025 + _S5031 + _S5031, _S5023 + _S5030 + _S5030);
    float3  _S5035 = - - _S5034;
    Matrix<float, 3, 3>  _S5036 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5037;
    (&_S5037)->primal_0 = _S4916;
    (&_S5037)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5038;
    (&_S5038)->primal_0 = t_29;
    (&_S5038)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5037, &_S5038, _S5035);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5039 = _S5038;
    Matrix<float, 3, 3>  _S5040 = transpose_0(_S5037.differential_0);
    DiffPair_float_0 _S5041;
    (&_S5041)->primal_0 = _S4915;
    (&_S5041)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5041, v_depth_10);
    float _S5042 = 0.3333333432674408f * _S5041.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5043;
    (&_S5043)->primal_0 = _S4914;
    (&_S5043)->differential_0 = _S4955;
    s_bwd_length_impl_0(&_S5043, _S5042);
    DiffPair_float_0 _S5044;
    (&_S5044)->primal_0 = _S4913;
    (&_S5044)->differential_0 = 0.0f;
    DiffPair_float_0 _S5045;
    (&_S5045)->primal_0 = _S4901;
    (&_S5045)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5044, &_S5045, 0.0f);
    DiffPair_float_0 _S5046;
    (&_S5046)->primal_0 = _S4851;
    (&_S5046)->differential_0 = 0.0f;
    DiffPair_float_0 _S5047;
    (&_S5047)->primal_0 = _S4876;
    (&_S5047)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5046, &_S5047, _S5044.differential_0);
    DiffPair_float_0 _S5048;
    (&_S5048)->primal_0 = _S4912;
    (&_S5048)->differential_0 = 0.0f;
    DiffPair_float_0 _S5049;
    (&_S5049)->primal_0 = _S4901;
    (&_S5049)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5048, &_S5049, 0.0f);
    DiffPair_float_0 _S5050;
    (&_S5050)->primal_0 = _S4851;
    (&_S5050)->differential_0 = 0.0f;
    DiffPair_float_0 _S5051;
    (&_S5051)->primal_0 = _S4876;
    (&_S5051)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5050, &_S5051, _S5048.differential_0);
    DiffPair_float_0 _S5052;
    (&_S5052)->primal_0 = _S4911;
    (&_S5052)->differential_0 = 0.0f;
    DiffPair_float_0 _S5053;
    (&_S5053)->primal_0 = _S4900;
    (&_S5053)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5052, &_S5053, 0.0f);
    DiffPair_float_0 _S5054;
    (&_S5054)->primal_0 = _S4850;
    (&_S5054)->differential_0 = 0.0f;
    DiffPair_float_0 _S5055;
    (&_S5055)->primal_0 = _S4875;
    (&_S5055)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5054, &_S5055, _S5052.differential_0);
    DiffPair_float_0 _S5056;
    (&_S5056)->primal_0 = _S4910;
    (&_S5056)->differential_0 = 0.0f;
    DiffPair_float_0 _S5057;
    (&_S5057)->primal_0 = _S4900;
    (&_S5057)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5056, &_S5057, 0.0f);
    DiffPair_float_0 _S5058;
    (&_S5058)->primal_0 = _S4850;
    (&_S5058)->differential_0 = 0.0f;
    DiffPair_float_0 _S5059;
    (&_S5059)->primal_0 = _S4875;
    (&_S5059)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5058, &_S5059, _S5056.differential_0);
    DiffPair_float_0 _S5060;
    (&_S5060)->primal_0 = _S4908;
    (&_S5060)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5060, 0.0f);
    float _S5061 = - (-1.0f * - (_S5060.differential_0 / _S4909));
    float2  _S5062 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5063;
    (&_S5063)->primal_0 = e2_7;
    (&_S5063)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5063, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5064;
    (&_S5064)->primal_0 = e1_15;
    (&_S5064)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5064, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5065;
    (&_S5065)->primal_0 = e0_15;
    (&_S5065)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5065, -0.0f);
    DiffPair_float_0 _S5066;
    (&_S5066)->primal_0 = _S4906;
    (&_S5066)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5066, 0.0f);
    float _S5067 = - _S5066.differential_0;
    float2  _S5068 = _S5064.differential_0 + make_float2 (_S4904 * _S5067, _S4902 * _S5066.differential_0);
    float2  _S5069 = - _S5068;
    float2  _S5070 = _S5065.differential_0 + make_float2 (_S4903 * _S5066.differential_0, _S4905 * _S5067);
    float2  _S5071 = - _S5070;
    float2  _S5072 = - _S5063.differential_0 + _S5068;
    float _S5073 = fx_35 * (_S5053.differential_0 + _S5057.differential_0 + _S5072.x);
    float2  _S5074 = make_float2 (_S5073, fy_35 * (_S5045.differential_0 + _S5049.differential_0 + _S5072.y)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S5073, (*dist_coeffs_40)[int(9)] * _S5073);
    float2  _S5075 = _S4889 * _S5074;
    float2  _S5076 = _S4893 * _S5074;
    float _S5077 = (*dist_coeffs_40)[int(4)] * _S5074.y;
    float _S5078 = (*dist_coeffs_40)[int(5)] * _S5074.x;
    float _S5079 = _S5075.x + _S5075.y;
    float _S5080 = r2_79 * _S5079;
    float _S5081 = r2_79 * _S5080;
    float _S5082 = (*dist_coeffs_40)[int(7)] * _S5074.y + _S5077 + (*dist_coeffs_40)[int(6)] * _S5074.x + _S5078 + _S4892 * _S5079 + _S4891 * _S5080 + _S4890 * _S5081 + (*dist_coeffs_40)[int(3)] * (r2_79 * _S5081);
    float _S5083 = v_79 * _S5082;
    float _S5084 = u_79 * _S5082;
    float _S5085 = _S4897 * _S5077 + 2.0f * (v_79 * _S5077) + _S4896 * _S5074.y + _S4894 * _S5074.x + _S5083 + _S5083;
    float _S5086 = _S4845 * (v_79 * _S5074.y) + _S4895 * _S5078 + 2.0f * (u_79 * _S5078) + _S4842 * (v_79 * _S5074.x) + _S5084 + _S5084;
    float3  _S5087 = _S4961.differential_0 + _S5043.differential_0;
    float3  _S5088 = _S4963 + _S4964 + _S5043.differential_0;
    float3  _S5089 = _S4962.differential_0 + _S5043.differential_0;
    FixedArray<float3 , 2>  _S5090;
    _S5090[int(0)] = _S4955;
    _S5090[int(1)] = _S4955;
    _S5090[int(1)] = _S4970;
    _S5090[int(0)] = _S4975;
    float3  _S5091 = _S5090[int(0)];
    float3  _S5092 = _S5090[int(1)];
    FixedArray<float3 , 16>  _S5093;
    _S5093[int(0)] = _S4955;
    _S5093[int(1)] = _S4955;
    _S5093[int(2)] = _S4955;
    _S5093[int(3)] = _S4955;
    _S5093[int(4)] = _S4955;
    _S5093[int(5)] = _S4955;
    _S5093[int(6)] = _S4955;
    _S5093[int(7)] = _S4955;
    _S5093[int(8)] = _S4955;
    _S5093[int(9)] = _S4955;
    _S5093[int(10)] = _S4955;
    _S5093[int(11)] = _S4955;
    _S5093[int(12)] = _S4955;
    _S5093[int(13)] = _S4955;
    _S5093[int(14)] = _S4955;
    _S5093[int(15)] = _S4955;
    _S5093[int(7)] = _S4999;
    _S5093[int(0)] = _S5033;
    _S5093[int(1)] = _S5020;
    _S5093[int(2)] = _S5018;
    _S5093[int(3)] = _S5016;
    _S5093[int(4)] = _S5005;
    _S5093[int(5)] = _S5003;
    _S5093[int(6)] = _S5001;
    _S5093[int(15)] = _S4977;
    _S5093[int(8)] = _S4997;
    _S5093[int(9)] = _S4989;
    _S5093[int(10)] = _S4987;
    _S5093[int(11)] = _S4985;
    _S5093[int(12)] = _S4983;
    _S5093[int(13)] = _S4981;
    _S5093[int(14)] = _S4979;
    float3  _S5094 = _S5093[int(0)];
    float3  _S5095 = _S5093[int(1)];
    float3  _S5096 = _S5093[int(2)];
    float3  _S5097 = _S5093[int(3)];
    float3  _S5098 = _S5093[int(4)];
    float3  _S5099 = _S5093[int(5)];
    float3  _S5100 = _S5093[int(6)];
    float3  _S5101 = _S5093[int(7)];
    float3  _S5102 = _S5093[int(8)];
    float3  _S5103 = _S5093[int(9)];
    float3  _S5104 = _S5093[int(10)];
    float3  _S5105 = _S5093[int(11)];
    float3  _S5106 = _S5093[int(12)];
    float3  _S5107 = _S5093[int(13)];
    float3  _S5108 = _S5093[int(14)];
    float3  _S5109 = _S5093[int(15)];
    float _S5110 = _S5054.differential_0 + _S5058.differential_0;
    float2  _S5111 = _S5063.differential_0 + _S5071;
    float _S5112 = _S5047.differential_0 + _S5051.differential_0;
    float _S5113 = _S5046.differential_0 + _S5050.differential_0;
    float2  _S5114 = _S5069 + _S5070;
    float _S5115 = _S5055.differential_0 + _S5059.differential_0;
    float2  _S5116 = make_float2 (0.0f, _S5061);
    float2  _S5117 = _S5076 + make_float2 (_S5086, _S5085);
    float2  _S5118 = _S4877 * _S5117;
    float2  _S5119 = _S4888 * _S5117;
    float _S5120 = _S5118.x + _S5118.y;
    if(_S4881)
    {
        float _S5121 = _S5120 / _S4883;
        float _S5122 = _S4884 * - _S5121;
        float _S5123 = _S4880 * (0.3333333432674408f * - (_S4879 * _S5121));
        k_12 = _S5123 + _S5123;
        _S4882 = _S5122;
        _S4883 = 0.0f;
    }
    else
    {
        float _S5124 = _S5120 / _S4882;
        float _S5125 = _S4880 * - _S5124;
        k_12 = _S4878 * _S5124;
        _S4882 = 0.0f;
        _S4883 = _S5125;
    }
    DiffPair_float_0 _S5126;
    (&_S5126)->primal_0 = _S4878;
    (&_S5126)->differential_0 = 0.0f;
    DiffPair_float_0 _S5127;
    (&_S5127)->primal_0 = _S4879;
    (&_S5127)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5126, &_S5127, k_12);
    float _S5128 = _S5127.differential_0 + _S4882;
    float _S5129 = _S5126.differential_0 + _S4883;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5130;
    (&_S5130)->primal_0 = _S4877;
    (&_S5130)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5130, _S5129);
    float2  _S5131 = _S5130.differential_0 + _S5119;
    float _S5132 = fx_35 * (_S5114.x + _S5115);
    float2  _S5133 = make_float2 (_S5132, fy_35 * (_S5114.y + _S5112)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S5132, (*dist_coeffs_40)[int(9)] * _S5132);
    float2  _S5134 = _S4864 * _S5133;
    float _S5135 = (*dist_coeffs_40)[int(4)] * _S5133.y;
    float _S5136 = (*dist_coeffs_40)[int(5)] * _S5133.x;
    float _S5137 = _S5134.x + _S5134.y;
    float _S5138 = r2_78 * _S5137;
    float _S5139 = r2_78 * _S5138;
    float _S5140 = (*dist_coeffs_40)[int(7)] * _S5133.y + _S5135 + (*dist_coeffs_40)[int(6)] * _S5133.x + _S5136 + _S4867 * _S5137 + _S4866 * _S5138 + _S4865 * _S5139 + (*dist_coeffs_40)[int(3)] * (r2_78 * _S5139);
    float _S5141 = v_78 * _S5140;
    float _S5142 = u_78 * _S5140;
    float3  _S5143 = _S5089 + make_float3 (_S5131.x, _S5131.y, _S5128);
    float2  _S5144 = _S4868 * _S5133 + make_float2 (_S4845 * (v_78 * _S5133.y) + _S4870 * _S5136 + 2.0f * (u_78 * _S5136) + _S4842 * (v_78 * _S5133.x) + _S5142 + _S5142, _S4872 * _S5135 + 2.0f * (v_78 * _S5135) + _S4871 * _S5133.y + _S4869 * _S5133.x + _S5141 + _S5141);
    float2  _S5145 = _S4852 * _S5144;
    float2  _S5146 = _S4863 * _S5144;
    float _S5147 = _S5145.x + _S5145.y;
    if(_S4856)
    {
        float _S5148 = _S5147 / _S4858;
        float _S5149 = _S4859 * - _S5148;
        float _S5150 = _S4855 * (0.3333333432674408f * - (_S4854 * _S5148));
        k_12 = _S5150 + _S5150;
        _S4857 = _S5149;
        _S4858 = 0.0f;
    }
    else
    {
        float _S5151 = _S5147 / _S4857;
        float _S5152 = _S4855 * - _S5151;
        k_12 = _S4853 * _S5151;
        _S4857 = 0.0f;
        _S4858 = _S5152;
    }
    DiffPair_float_0 _S5153;
    (&_S5153)->primal_0 = _S4853;
    (&_S5153)->differential_0 = 0.0f;
    DiffPair_float_0 _S5154;
    (&_S5154)->primal_0 = _S4854;
    (&_S5154)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5153, &_S5154, k_12);
    float _S5155 = _S5154.differential_0 + _S4857;
    float _S5156 = _S5153.differential_0 + _S4858;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5157;
    (&_S5157)->primal_0 = _S4852;
    (&_S5157)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5157, _S5156);
    float2  _S5158 = _S5157.differential_0 + _S5146;
    float _S5159 = fx_35 * (_S5111.x + _S5110);
    float2  _S5160 = make_float2 (_S5159, fy_35 * (_S5111.y + _S5113)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S5159, (*dist_coeffs_40)[int(9)] * _S5159);
    float2  _S5161 = _S4837 * _S5160;
    float _S5162 = (*dist_coeffs_40)[int(4)] * _S5160.y;
    float _S5163 = (*dist_coeffs_40)[int(5)] * _S5160.x;
    float _S5164 = _S5161.x + _S5161.y;
    float _S5165 = r2_77 * _S5164;
    float _S5166 = r2_77 * _S5165;
    float _S5167 = (*dist_coeffs_40)[int(7)] * _S5160.y + _S5162 + (*dist_coeffs_40)[int(6)] * _S5160.x + _S5163 + _S4840 * _S5164 + _S4839 * _S5165 + _S4838 * _S5166 + (*dist_coeffs_40)[int(3)] * (r2_77 * _S5166);
    float _S5168 = v_77 * _S5167;
    float _S5169 = u_77 * _S5167;
    float3  _S5170 = _S5087 + make_float3 (_S5158.x, _S5158.y, _S5155);
    float2  _S5171 = _S4841 * _S5160 + make_float2 (_S4845 * (v_77 * _S5160.y) + _S4844 * _S5163 + 2.0f * (u_77 * _S5163) + _S4842 * (v_77 * _S5160.x) + _S5169 + _S5169, _S4847 * _S5162 + 2.0f * (v_77 * _S5162) + _S4846 * _S5160.y + _S4843 * _S5160.x + _S5168 + _S5168);
    float2  _S5172 = _S4825 * _S5171;
    float2  _S5173 = _S4836 * _S5171;
    float _S5174 = _S5172.x + _S5172.y;
    if(_S4829)
    {
        float _S5175 = _S5174 / _S4831;
        float _S5176 = _S4832 * - _S5175;
        float _S5177 = _S4828 * (0.3333333432674408f * - (_S4827 * _S5175));
        k_12 = _S5177 + _S5177;
        _S4830 = _S5176;
        _S4831 = 0.0f;
    }
    else
    {
        float _S5178 = _S5174 / _S4830;
        float _S5179 = _S4828 * - _S5178;
        k_12 = _S4826 * _S5178;
        _S4830 = 0.0f;
        _S4831 = _S5179;
    }
    DiffPair_float_0 _S5180;
    (&_S5180)->primal_0 = _S4826;
    (&_S5180)->differential_0 = 0.0f;
    DiffPair_float_0 _S5181;
    (&_S5181)->primal_0 = _S4827;
    (&_S5181)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5180, &_S5181, k_12);
    float _S5182 = _S5181.differential_0 + _S4830;
    float _S5183 = _S5180.differential_0 + _S4831;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5184;
    (&_S5184)->primal_0 = _S4825;
    (&_S5184)->differential_0 = _S5062;
    s_bwd_length_impl_1(&_S5184, _S5183);
    float2  _S5185 = _S5184.differential_0 + _S5173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5186;
    (&_S5186)->primal_0 = vert2_c_11;
    (&_S5186)->differential_0 = _S4955;
    s_bwd_length_impl_0(&_S5186, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5187;
    (&_S5187)->primal_0 = vert1_c_11;
    (&_S5187)->differential_0 = _S4955;
    s_bwd_length_impl_0(&_S5187, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5188;
    (&_S5188)->primal_0 = vert0_c_11;
    (&_S5188)->differential_0 = _S4955;
    s_bwd_length_impl_0(&_S5188, 0.0f);
    float3  _S5189 = _S5186.differential_0 + _S5143;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5190;
    (&_S5190)->primal_0 = R_30;
    (&_S5190)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5191;
    (&_S5191)->primal_0 = vert2_7;
    (&_S5191)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5190, &_S5191, _S5189);
    float3  _S5192 = _S5187.differential_0 + _S5170;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5193;
    (&_S5193)->primal_0 = R_30;
    (&_S5193)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5194;
    (&_S5194)->primal_0 = vert1_7;
    (&_S5194)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5193, &_S5194, _S5192);
    float3  _S5195 = _S5188.differential_0 + _S5088 + make_float3 (_S5185.x, _S5185.y, _S5182);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5196;
    (&_S5196)->primal_0 = R_30;
    (&_S5196)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5197;
    (&_S5197)->primal_0 = vert0_7;
    (&_S5197)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5196, &_S5197, _S5195);
    float3  _S5198 = _S5191.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5199;
    (&_S5199)->primal_0 = _S4819;
    (&_S5199)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5200;
    (&_S5200)->primal_0 = _S4824;
    (&_S5200)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5199, &_S5200, _S5198);
    float _S5201 = - _S5200.differential_0.y;
    float _S5202 = _S4823 * _S5200.differential_0.x;
    float _S5203 = - (_S4815 * _S5200.differential_0.x);
    float3  _S5204 = _S5194.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5205;
    (&_S5205)->primal_0 = _S4819;
    (&_S5205)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5206;
    (&_S5206)->primal_0 = _S4822;
    (&_S5206)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5205, &_S5206, _S5204);
    float _S5207 = _S4815 * _S5206.differential_0.x;
    float _S5208 = _S4821 * _S5206.differential_0.x;
    float3  _S5209 = _S5197.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5210;
    (&_S5210)->primal_0 = _S4819;
    (&_S5210)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5211;
    (&_S5211)->primal_0 = _S4820;
    (&_S5211)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5210, &_S5211, _S5209);
    Matrix<float, 3, 3>  _S5212 = transpose_0(_S5199.differential_0 + _S5205.differential_0 + _S5210.differential_0);
    float _S5213 = 2.0f * - _S5212.rows[int(2)].z;
    float _S5214 = 2.0f * _S5212.rows[int(2)].y;
    float _S5215 = 2.0f * _S5212.rows[int(2)].x;
    float _S5216 = 2.0f * _S5212.rows[int(1)].z;
    float _S5217 = 2.0f * - _S5212.rows[int(1)].y;
    float _S5218 = 2.0f * _S5212.rows[int(1)].x;
    float _S5219 = 2.0f * _S5212.rows[int(0)].z;
    float _S5220 = 2.0f * _S5212.rows[int(0)].y;
    float _S5221 = 2.0f * - _S5212.rows[int(0)].x;
    float _S5222 = - _S5218 + _S5220;
    float _S5223 = _S5215 + - _S5219;
    float _S5224 = - _S5214 + _S5216;
    float _S5225 = _S5214 + _S5216;
    float _S5226 = _S5215 + _S5219;
    float _S5227 = _S5218 + _S5220;
    float _S5228 = quat_31.w * (_S5217 + _S5221);
    float _S5229 = quat_31.z * (_S5213 + _S5221);
    float _S5230 = quat_31.y * (_S5213 + _S5217);
    float _S5231 = quat_31.x * _S5222 + quat_31.z * _S5225 + quat_31.y * _S5226 + _S5228 + _S5228;
    float _S5232 = quat_31.x * _S5223 + quat_31.w * _S5225 + quat_31.y * _S5227 + _S5229 + _S5229;
    float _S5233 = quat_31.x * _S5224 + quat_31.w * _S5226 + quat_31.z * _S5227 + _S5230 + _S5230;
    float _S5234 = quat_31.w * _S5222 + quat_31.z * _S5223 + quat_31.y * _S5224;
    float _S5235 = _S5203 + _S5207;
    float _S5236 = 0.5f * - _S5235;
    float _S5237 = _S5201 + _S5206.differential_0.y;
    DiffPair_float_0 _S5238;
    (&_S5238)->primal_0 = _S4816;
    (&_S5238)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5238, _S5237);
    float _S5239 = _S5236 + _S5238.differential_0;
    float _S5240 = _S5202 + _S5208 + _S5211.differential_0.x;
    DiffPair_float_0 _S5241;
    (&_S5241)->primal_0 = _S4814;
    (&_S5241)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5241, _S5240);
    float _S5242 = _S5236 + _S5241.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5243;
    (&_S5243)->primal_0 = mean_c_26;
    (&_S5243)->differential_0 = _S4955;
    s_bwd_length_impl_0(&_S5243, 0.0f);
    float3  _S5244 = _S5243.differential_0 + _S4958.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5245;
    (&_S5245)->primal_0 = R_30;
    (&_S5245)->differential_0 = _S5036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5246;
    (&_S5246)->primal_0 = mean_30;
    (&_S5246)->differential_0 = _S4955;
    s_bwd_prop_mul_1(&_S5245, &_S5246, _S5244);
    float3  _S5247 = _S5189 + _S5192 + _S5195 + _S5244 + _S5039.differential_0;
    Matrix<float, 3, 3>  _S5248 = _S5190.differential_0 + _S5193.differential_0 + _S5196.differential_0 + _S5245.differential_0 + _S5040;
    float3  _S5249 = make_float3 (_S5242, _S5239, _S5235);
    float4  _S5250 = make_float4 (0.0f);
    *&((&_S5250)->w) = _S5231;
    *&((&_S5250)->z) = _S5232;
    *&((&_S5250)->y) = _S5233;
    *&((&_S5250)->x) = _S5234;
    float4  _S5251 = _S5250;
    float3  _S5252 = _S5198 + _S5204 + _S5209 + _S5246.differential_0 + _S5034;
    *v_mean_10 = _S5252;
    *v_quat_9 = _S5251;
    *v_scale_9 = _S5249;
    *v_hardness_5 = _S5116;
    (*v_sh_coeffs_8)[int(0)] = _S5094;
    (*v_sh_coeffs_8)[int(1)] = _S5095;
    (*v_sh_coeffs_8)[int(2)] = _S5096;
    (*v_sh_coeffs_8)[int(3)] = _S5097;
    (*v_sh_coeffs_8)[int(4)] = _S5098;
    (*v_sh_coeffs_8)[int(5)] = _S5099;
    (*v_sh_coeffs_8)[int(6)] = _S5100;
    (*v_sh_coeffs_8)[int(7)] = _S5101;
    (*v_sh_coeffs_8)[int(8)] = _S5102;
    (*v_sh_coeffs_8)[int(9)] = _S5103;
    (*v_sh_coeffs_8)[int(10)] = _S5104;
    (*v_sh_coeffs_8)[int(11)] = _S5105;
    (*v_sh_coeffs_8)[int(12)] = _S5106;
    (*v_sh_coeffs_8)[int(13)] = _S5107;
    (*v_sh_coeffs_8)[int(14)] = _S5108;
    (*v_sh_coeffs_8)[int(15)] = _S5109;
    (*v_ch_coeffs_3)[int(0)] = _S5091;
    (*v_ch_coeffs_3)[int(1)] = _S5092;
    *v_R_9 = _S5248;
    *v_t_9 = _S5247;
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle(FixedArray<float3 , 3>  * verts_4, float2  hardness_16, float3  ray_o_5, float3  ray_d_5)
{
    float3  v1v0_0 = (*verts_4)[int(1)] - (*verts_4)[int(0)];
    float3  v2v0_0 = (*verts_4)[int(2)] - (*verts_4)[int(0)];
    float3  rov0_0 = ray_o_5 - (*verts_4)[int(0)];
    float3  n_0 = cross_0(v1v0_0, v2v0_0);
    float3  q_2 = cross_0(rov0_0, ray_d_5);
    float d_28 = 1.0f / dot_0(ray_d_5, n_0);
    float u_80 = d_28 * dot_0(- q_2, v2v0_0);
    float v_80 = d_28 * dot_0(q_2, v1v0_0);
    float t_30 = d_28 * dot_0(- n_0, rov0_0);
    bool _S5253;
    if(u_80 >= 0.0f)
    {
        _S5253 = v_80 >= 0.0f;
    }
    else
    {
        _S5253 = false;
    }
    if(_S5253)
    {
        _S5253 = (u_80 + v_80) <= 1.0f;
    }
    else
    {
        _S5253 = false;
    }
    if(_S5253)
    {
        _S5253 = t_30 >= 0.0f;
    }
    else
    {
        _S5253 = false;
    }
    if(!_S5253)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_80), (v_80)))), ((F32_sqrt((0.5f))) * (1.0f - u_80 - v_80)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5254;
    if(opac_0 < 0.0f)
    {
        _S5254 = 0.0f;
    }
    else
    {
        _S5254 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5254;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5255 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5256 = *dpray_d_4;
    float3  v1v0_1 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_1 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_1 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5257 = s_primal_ctx_cross_0(v1v0_1, v2v0_1);
    float3  _S5258 = s_primal_ctx_cross_0(rov0_1, (*dpray_d_4).primal_0);
    float _S5259 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S5257);
    float d_29 = 1.0f / _S5259;
    float _S5260 = _S5259 * _S5259;
    float3  _S5261 = - _S5258;
    float _S5262 = s_primal_ctx_dot_0(_S5261, v2v0_1);
    float u_81 = d_29 * _S5262;
    float _S5263 = s_primal_ctx_dot_0(_S5258, v1v0_1);
    float v_81 = d_29 * _S5263;
    float3  _S5264 = - _S5257;
    float t_31 = d_29 * s_primal_ctx_dot_0(_S5264, rov0_1);
    bool _S5265;
    if(u_81 >= 0.0f)
    {
        _S5265 = v_81 >= 0.0f;
    }
    else
    {
        _S5265 = false;
    }
    if(_S5265)
    {
        _S5265 = (u_81 + v_81) <= 1.0f;
    }
    else
    {
        _S5265 = false;
    }
    if(_S5265)
    {
        _S5265 = t_31 >= 0.0f;
    }
    else
    {
        _S5265 = false;
    }
    bool _S5266 = !!_S5265;
    float _S5267;
    float _S5268;
    float _S5269;
    float _S5270;
    float _S5271;
    float _S5272;
    float _S5273;
    float _S5274;
    float _S5275;
    float _S5276;
    float _S5277;
    if(_S5266)
    {
        float _S5278 = s_primal_ctx_min_0(u_81, v_81);
        float _S5279 = s_primal_ctx_sqrt_0(0.5f);
        float _S5280 = _S5279 * (1.0f - u_81 - v_81);
        float _S5281 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5278, _S5280) * _S5281;
        float _S5282 = _S5255.primal_0.y;
        float _S5283 = 1.0f - opac_1;
        float _S5284 = 1.0f - s_primal_ctx_clamp_0(_S5282, 0.0f, 0.99989998340606689f);
        float _S5285 = 1.0f / _S5284;
        float _S5286 = _S5284 * _S5284;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5283, _S5285);
        float o_1 = _S5255.primal_0.x;
        bool _S5287 = opac_1 < 0.0f;
        if(_S5287)
        {
            _S5267 = 0.0f;
        }
        else
        {
            _S5267 = o_1 * w_1;
        }
        _S5265 = _S5287;
        _S5268 = o_1;
        _S5269 = w_1;
        _S5270 = _S5283;
        _S5271 = _S5285;
        _S5272 = _S5286;
        _S5273 = _S5282;
        _S5274 = _S5281;
        _S5275 = _S5278;
        _S5276 = _S5280;
        _S5277 = _S5279;
    }
    else
    {
        _S5265 = false;
        _S5267 = 0.0f;
        _S5268 = 0.0f;
        _S5269 = 0.0f;
        _S5270 = 0.0f;
        _S5271 = 0.0f;
        _S5272 = 0.0f;
        _S5273 = 0.0f;
        _S5274 = 0.0f;
        _S5275 = 0.0f;
        _S5276 = 0.0f;
        _S5277 = 0.0f;
    }
    float2  _S5288 = make_float2 (0.0f);
    float2  _S5289;
    if(_S5266)
    {
        if(_S5265)
        {
            _S5267 = 0.0f;
            _S5268 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5290;
            (&_S5290)->primal_0 = _S5267;
            (&_S5290)->differential_0 = 0.0f;
            DiffPair_float_0 _S5291;
            (&_S5291)->primal_0 = 0.99500000476837158f;
            (&_S5291)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5290, &_S5291, _s_dOut_11);
            float _S5292 = _S5268 * _S5290.differential_0;
            _S5267 = _S5269 * _S5290.differential_0;
            _S5268 = _S5292;
        }
        float _S5293 = - _S5268;
        DiffPair_float_0 _S5294;
        (&_S5294)->primal_0 = _S5270;
        (&_S5294)->differential_0 = 0.0f;
        DiffPair_float_0 _S5295;
        (&_S5295)->primal_0 = _S5271;
        (&_S5295)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5294, &_S5295, _S5293);
        float _S5296 = - - (_S5295.differential_0 / _S5272);
        float s_diff_opac_T_0 = - _S5294.differential_0;
        DiffPair_float_0 _S5297;
        (&_S5297)->primal_0 = _S5273;
        (&_S5297)->differential_0 = 0.0f;
        DiffPair_float_0 _S5298;
        (&_S5298)->primal_0 = 0.0f;
        (&_S5298)->differential_0 = 0.0f;
        DiffPair_float_0 _S5299;
        (&_S5299)->primal_0 = 0.99989998340606689f;
        (&_S5299)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5297, &_S5298, &_S5299, _S5296);
        float _S5300 = _S5274 * s_diff_opac_T_0;
        DiffPair_float_0 _S5301;
        (&_S5301)->primal_0 = _S5275;
        (&_S5301)->differential_0 = 0.0f;
        DiffPair_float_0 _S5302;
        (&_S5302)->primal_0 = _S5276;
        (&_S5302)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5301, &_S5302, _S5300);
        float _S5303 = - (_S5277 * _S5302.differential_0);
        DiffPair_float_0 _S5304;
        (&_S5304)->primal_0 = u_81;
        (&_S5304)->differential_0 = 0.0f;
        DiffPair_float_0 _S5305;
        (&_S5305)->primal_0 = v_81;
        (&_S5305)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5304, &_S5305, _S5301.differential_0);
        float2  _S5306 = make_float2 (_S5267, _S5297.differential_0);
        float _S5307 = _S5303 + _S5305.differential_0;
        _S5267 = _S5303 + _S5304.differential_0;
        _S5268 = _S5307;
        _S5289 = _S5306;
    }
    else
    {
        _S5267 = 0.0f;
        _S5268 = 0.0f;
        _S5289 = _S5288;
    }
    float3  _S5308 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5309;
    (&_S5309)->primal_0 = _S5264;
    (&_S5309)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5310;
    (&_S5310)->primal_0 = rov0_1;
    (&_S5310)->differential_0 = _S5308;
    s_bwd_prop_dot_0(&_S5309, &_S5310, 0.0f);
    float3  _S5311 = - _S5309.differential_0;
    float _S5312 = d_29 * _S5268;
    float _S5313 = _S5263 * _S5268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5314;
    (&_S5314)->primal_0 = _S5258;
    (&_S5314)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5315;
    (&_S5315)->primal_0 = v1v0_1;
    (&_S5315)->differential_0 = _S5308;
    s_bwd_prop_dot_0(&_S5314, &_S5315, _S5312);
    float _S5316 = d_29 * _S5267;
    float _S5317 = _S5262 * _S5267;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5318;
    (&_S5318)->primal_0 = _S5261;
    (&_S5318)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5319;
    (&_S5319)->primal_0 = v2v0_1;
    (&_S5319)->differential_0 = _S5308;
    s_bwd_prop_dot_0(&_S5318, &_S5319, _S5316);
    float3  _S5320 = - _S5318.differential_0;
    float _S5321 = - ((_S5313 + _S5317) / _S5260);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5322;
    (&_S5322)->primal_0 = _S5256.primal_0;
    (&_S5322)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5323;
    (&_S5323)->primal_0 = _S5257;
    (&_S5323)->differential_0 = _S5308;
    s_bwd_prop_dot_0(&_S5322, &_S5323, _S5321);
    float3  _S5324 = _S5314.differential_0 + _S5320;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5325;
    (&_S5325)->primal_0 = rov0_1;
    (&_S5325)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5326;
    (&_S5326)->primal_0 = _S5256.primal_0;
    (&_S5326)->differential_0 = _S5308;
    s_bwd_prop_cross_0(&_S5325, &_S5326, _S5324);
    float3  _S5327 = _S5311 + _S5323.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5328;
    (&_S5328)->primal_0 = v1v0_1;
    (&_S5328)->differential_0 = _S5308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5329;
    (&_S5329)->primal_0 = v2v0_1;
    (&_S5329)->differential_0 = _S5308;
    s_bwd_prop_cross_0(&_S5328, &_S5329, _S5327);
    float3  _S5330 = _S5310.differential_0 + _S5325.differential_0;
    float3  _S5331 = _S5319.differential_0 + _S5329.differential_0;
    float3  _S5332 = _S5315.differential_0 + _S5328.differential_0;
    float3  _S5333 = - _S5330 + - _S5331 + - _S5332;
    float3  _S5334 = _S5322.differential_0 + _S5326.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S5334;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S5330;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5289;
    FixedArray<float3 , 3>  _S5335;
    _S5335[int(0)] = _S5308;
    _S5335[int(1)] = _S5308;
    _S5335[int(2)] = _S5308;
    _S5335[int(2)] = _S5331;
    _S5335[int(0)] = _S5333;
    _S5335[int(1)] = _S5332;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5335;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5336, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5337, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5338, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5339, float _S5340)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5336, _S5337, _S5338, _S5339, _S5340);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_5, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5341 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5342 = { _S5341, _S5341, _S5341 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_5;
    (&dp_verts_0)->differential_0 = _S5342;
    float2  _S5343 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5343;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S5341;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S5341;
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
    float d_30 = 1.0f / dot_0(ray_d_8, n_2);
    float u_82 = d_30 * dot_0(- q_3, v2v0_2);
    float v_82 = d_30 * dot_0(q_3, v1v0_2);
    *depth_20 = d_30 * dot_0(- n_2, rov0_2);
    *color_13 = (*rgbs_5)[int(0)] * make_float3 (1.0f - u_82 - v_82) + (*rgbs_5)[int(1)] * make_float3 (u_82) + (*rgbs_5)[int(2)] * make_float3 (v_82);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_3 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_3 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_3 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S5344 = s_primal_ctx_cross_0(v1v0_3, v2v0_3);
    float3  _S5345 = s_primal_ctx_cross_0(rov0_3, (*dpray_d_5).primal_0);
    float _S5346 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S5344);
    float d_31 = 1.0f / _S5346;
    float _S5347 = _S5346 * _S5346;
    float3  _S5348 = - _S5345;
    float _S5349 = s_primal_ctx_dot_0(_S5348, v2v0_3);
    float u_83 = d_31 * _S5349;
    float3  _S5350 = make_float3 (u_83);
    float _S5351 = s_primal_ctx_dot_0(_S5345, v1v0_3);
    float v_83 = d_31 * _S5351;
    float3  _S5352 = make_float3 (v_83);
    float3  _S5353 = - _S5344;
    float _S5354 = s_primal_ctx_dot_0(_S5353, rov0_3);
    float _S5355 = d_31 * _S5354;
    float3  _S5356 = make_float3 (1.0f - u_83 - v_83);
    DiffPair_float_0 _S5357;
    (&_S5357)->primal_0 = s_primal_ctx_max_0(_S5355, 9.999999960041972e-13f);
    (&_S5357)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5357, dpdepth_2);
    DiffPair_float_0 _S5358;
    (&_S5358)->primal_0 = _S5355;
    (&_S5358)->differential_0 = 0.0f;
    DiffPair_float_0 _S5359;
    (&_S5359)->primal_0 = 9.999999960041972e-13f;
    (&_S5359)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5358, &_S5359, _S5357.differential_0);
    float3  _S5360 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5361 = _S5352 * dpcolor_1;
    float3  _S5362 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5363 = _S5350 * dpcolor_1;
    float3  _S5364 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5365 = _S5356 * dpcolor_1;
    float _S5366 = - (_S5364.x + _S5364.y + _S5364.z);
    float _S5367 = d_31 * _S5358.differential_0;
    float _S5368 = _S5354 * _S5358.differential_0;
    float3  _S5369 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5370;
    (&_S5370)->primal_0 = _S5353;
    (&_S5370)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5371;
    (&_S5371)->primal_0 = rov0_3;
    (&_S5371)->differential_0 = _S5369;
    s_bwd_prop_dot_0(&_S5370, &_S5371, _S5367);
    float3  _S5372 = - _S5370.differential_0;
    float _S5373 = _S5366 + _S5360.x + _S5360.y + _S5360.z;
    float _S5374 = d_31 * _S5373;
    float _S5375 = _S5351 * _S5373;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5376;
    (&_S5376)->primal_0 = _S5345;
    (&_S5376)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5377;
    (&_S5377)->primal_0 = v1v0_3;
    (&_S5377)->differential_0 = _S5369;
    s_bwd_prop_dot_0(&_S5376, &_S5377, _S5374);
    float _S5378 = _S5366 + _S5362.x + _S5362.y + _S5362.z;
    float _S5379 = d_31 * _S5378;
    float _S5380 = _S5349 * _S5378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5381;
    (&_S5381)->primal_0 = _S5348;
    (&_S5381)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5382;
    (&_S5382)->primal_0 = v2v0_3;
    (&_S5382)->differential_0 = _S5369;
    s_bwd_prop_dot_0(&_S5381, &_S5382, _S5379);
    float3  _S5383 = - _S5381.differential_0;
    float _S5384 = - ((_S5368 + _S5375 + _S5380) / _S5347);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5385;
    (&_S5385)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5385)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5386;
    (&_S5386)->primal_0 = _S5344;
    (&_S5386)->differential_0 = _S5369;
    s_bwd_prop_dot_0(&_S5385, &_S5386, _S5384);
    float3  _S5387 = _S5376.differential_0 + _S5383;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5388;
    (&_S5388)->primal_0 = rov0_3;
    (&_S5388)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5389;
    (&_S5389)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5389)->differential_0 = _S5369;
    s_bwd_prop_cross_0(&_S5388, &_S5389, _S5387);
    float3  _S5390 = _S5372 + _S5386.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5391;
    (&_S5391)->primal_0 = v1v0_3;
    (&_S5391)->differential_0 = _S5369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5392;
    (&_S5392)->primal_0 = v2v0_3;
    (&_S5392)->differential_0 = _S5369;
    s_bwd_prop_cross_0(&_S5391, &_S5392, _S5390);
    float3  _S5393 = _S5371.differential_0 + _S5388.differential_0;
    float3  _S5394 = _S5382.differential_0 + _S5392.differential_0;
    float3  _S5395 = _S5377.differential_0 + _S5391.differential_0;
    float3  _S5396 = - _S5393 + - _S5394 + - _S5395;
    float3  _S5397 = _S5385.differential_0 + _S5389.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S5397;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S5393;
    FixedArray<float3 , 3>  _S5398;
    _S5398[int(0)] = _S5369;
    _S5398[int(1)] = _S5369;
    _S5398[int(2)] = _S5369;
    _S5398[int(2)] = _S5361;
    _S5398[int(1)] = _S5363;
    _S5398[int(0)] = _S5365;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5398;
    FixedArray<float3 , 3>  _S5399;
    _S5399[int(0)] = _S5369;
    _S5399[int(1)] = _S5369;
    _S5399[int(2)] = _S5369;
    _S5399[int(2)] = _S5394;
    _S5399[int(0)] = _S5396;
    _S5399[int(1)] = _S5395;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5399;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5400, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5401, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5402, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5403, float3  _S5404, float _S5405)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5400, _S5401, _S5402, _S5403, _S5404, _S5405);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_8, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_9, float3  ray_d_9, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5406 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5407 = { _S5406, _S5406, _S5406 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_8;
    (&dp_verts_1)->differential_0 = _S5407;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5407;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_9;
    (&dp_ray_o_3)->differential_0 = _S5406;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_9;
    (&dp_ray_d_3)->differential_0 = _S5406;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

