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
    float3  grd_1 = mul_0(iscl_rot_2, ray_d_2);
    *depth_10 = (F32_log(((F32_max((- dot_0(mul_0(iscl_rot_2, ray_o_2 - mean_18), grd_1) / dot_0(grd_1, grd_1)), (9.99999997475242708e-07f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3007 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S3008 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S3007);
    float3  _S3009 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S3010 = - s_primal_ctx_dot_0(_S3008, _S3009);
    float _S3011 = s_primal_ctx_dot_0(_S3009, _S3009);
    float _S3012 = _S3010 / _S3011;
    float _S3013 = _S3011 * _S3011;
    DiffPair_float_0 _S3014;
    (&_S3014)->primal_0 = s_primal_ctx_max_0(_S3012, 9.99999997475242708e-07f);
    (&_S3014)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3014, dpdepth_0);
    DiffPair_float_0 _S3015;
    (&_S3015)->primal_0 = _S3012;
    (&_S3015)->differential_0 = 0.0f;
    DiffPair_float_0 _S3016;
    (&_S3016)->primal_0 = 9.99999997475242708e-07f;
    (&_S3016)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3015, &_S3016, _S3014.differential_0);
    float _S3017 = _S3015.differential_0 / _S3013;
    float _S3018 = _S3010 * - _S3017;
    float _S3019 = _S3011 * _S3017;
    float3  _S3020 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3021;
    (&_S3021)->primal_0 = _S3009;
    (&_S3021)->differential_0 = _S3020;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3022;
    (&_S3022)->primal_0 = _S3009;
    (&_S3022)->differential_0 = _S3020;
    s_bwd_prop_dot_0(&_S3021, &_S3022, _S3018);
    float _S3023 = - _S3019;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3024;
    (&_S3024)->primal_0 = _S3008;
    (&_S3024)->differential_0 = _S3020;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3025;
    (&_S3025)->primal_0 = _S3009;
    (&_S3025)->differential_0 = _S3020;
    s_bwd_prop_dot_0(&_S3024, &_S3025, _S3023);
    float3  _S3026 = _S3022.differential_0 + _S3021.differential_0 + _S3025.differential_0;
    Matrix<float, 3, 3>  _S3027 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3028;
    (&_S3028)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S3028)->differential_0 = _S3027;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3029;
    (&_S3029)->primal_0 = (*dpray_d_1).primal_0;
    (&_S3029)->differential_0 = _S3020;
    s_bwd_prop_mul_1(&_S3028, &_S3029, _S3026);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3030;
    (&_S3030)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S3030)->differential_0 = _S3027;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3031;
    (&_S3031)->primal_0 = _S3007;
    (&_S3031)->differential_0 = _S3020;
    s_bwd_prop_mul_1(&_S3030, &_S3031, _S3024.differential_0);
    float3  _S3032 = - _S3031.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S3029.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S3031.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3033 = _S3028.differential_0 + _S3030.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3033;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3032;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3034, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3035, DiffPair_float_0 * _S3036, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3037, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3038, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3039, float3  _S3040, float _S3041)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3034, _S3035, _S3036, _S3037, _S3038, _S3039, _S3040, _S3041);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_19, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3042 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_19;
    (&dp_mean_1)->differential_0 = _S3042;
    Matrix<float, 3, 3>  _S3043 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3043;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S3042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S3042;
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
    float _S3044 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S3044;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_18)
{
    float _S3045 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3045;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  dOut_19)
{
    float3  _S3046 = make_float3 (1.0f) / (*dpx_15).primal_0 * dOut_19;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3046;
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
        float _S3047 = mean_c_15.z;
        bool _S3048;
        if(_S3047 < near_plane_10)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3047 > far_plane_10;
        }
        if(_S3048)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3049 = scale_20.x;
        float sx_1 = (F32_exp((_S3049)));
        float _S3050 = scale_20.y;
        float sy_1 = (F32_exp((_S3050)));
        float sz_1 = scale_20.z - 0.5f * (_S3049 + _S3050);
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
        Matrix<float, 3, 3>  _S3051 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_36), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_36), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_0 = mul_0(R_21, mul_0(_S3051, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_18;
        float3  vert1_c_0 = mul_0(R_21, mul_0(_S3051, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_18;
        float3  vert2_c_0 = mul_0(R_21, mul_0(_S3051, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_18;
        float _S3052 = vert0_c_0.z;
        float _S3053 = vert1_c_0.z;
        float _S3054 = vert2_c_0.z;
        if(_S3052 < near_plane_10)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3052 > far_plane_10;
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3053 < near_plane_10;
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3053 > far_plane_10;
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3054 < near_plane_10;
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = _S3054 > far_plane_10;
        }
        if(_S3048)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        for(;;)
        {
            *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3052);
            if(_S3052 < 0.0f)
            {
                _S3048 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3055 = camera_distortion_jac_0(*uv0_0, dist_coeffs_28);
                _S3048 = !((F32_min((determinant_0(_S3055)), ((F32_min((_S3055.rows[int(0)].x), (_S3055.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3048)
            {
                break;
            }
            float u_44 = (*uv0_0).x;
            float v_44 = (*uv0_0).y;
            float r2_44 = u_44 * u_44 + v_44 * v_44;
            float2  _S3056 = *uv0_0 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_28)[int(0)] + r2_44 * ((*dist_coeffs_28)[int(1)] + r2_44 * ((*dist_coeffs_28)[int(2)] + r2_44 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_44 * v_44 + (*dist_coeffs_28)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_28)[int(6)] * r2_44, 2.0f * (*dist_coeffs_28)[int(5)] * u_44 * v_44 + (*dist_coeffs_28)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_28)[int(7)] * r2_44);
            float2  _S3057 = _S3056 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3056.x + (*dist_coeffs_28)[int(9)] * _S3056.y, 0.0f);
            *uv0_0 = make_float2 (fx_24 * _S3057.x + cx_19, fy_24 * _S3057.y + cy_19);
            break;
        }
        bool all_valid_8 = true & (!_S3048);
        for(;;)
        {
            *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3053);
            if(_S3053 < 0.0f)
            {
                _S3048 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3058 = camera_distortion_jac_0(*uv1_0, dist_coeffs_28);
                _S3048 = !((F32_min((determinant_0(_S3058)), ((F32_min((_S3058.rows[int(0)].x), (_S3058.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3048)
            {
                break;
            }
            float u_45 = (*uv1_0).x;
            float v_45 = (*uv1_0).y;
            float r2_45 = u_45 * u_45 + v_45 * v_45;
            float2  _S3059 = *uv1_0 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_28)[int(0)] + r2_45 * ((*dist_coeffs_28)[int(1)] + r2_45 * ((*dist_coeffs_28)[int(2)] + r2_45 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_45 * v_45 + (*dist_coeffs_28)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_28)[int(6)] * r2_45, 2.0f * (*dist_coeffs_28)[int(5)] * u_45 * v_45 + (*dist_coeffs_28)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_28)[int(7)] * r2_45);
            float2  _S3060 = _S3059 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3059.x + (*dist_coeffs_28)[int(9)] * _S3059.y, 0.0f);
            *uv1_0 = make_float2 (fx_24 * _S3060.x + cx_19, fy_24 * _S3060.y + cy_19);
            break;
        }
        bool all_valid_9 = all_valid_8 & (!_S3048);
        for(;;)
        {
            *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3054);
            if(_S3054 < 0.0f)
            {
                _S3048 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3061 = camera_distortion_jac_0(*uv2_0, dist_coeffs_28);
                _S3048 = !((F32_min((determinant_0(_S3061)), ((F32_min((_S3061.rows[int(0)].x), (_S3061.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3048)
            {
                break;
            }
            float u_46 = (*uv2_0).x;
            float v_46 = (*uv2_0).y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S3062 = *uv2_0 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_28)[int(0)] + r2_46 * ((*dist_coeffs_28)[int(1)] + r2_46 * ((*dist_coeffs_28)[int(2)] + r2_46 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_46 * v_46 + (*dist_coeffs_28)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_28)[int(6)] * r2_46, 2.0f * (*dist_coeffs_28)[int(5)] * u_46 * v_46 + (*dist_coeffs_28)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_28)[int(7)] * r2_46);
            float2  _S3063 = _S3062 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3062.x + (*dist_coeffs_28)[int(9)] * _S3062.y, 0.0f);
            *uv2_0 = make_float2 (fx_24 * _S3063.x + cx_19, fy_24 * _S3063.y + cy_19);
            break;
        }
        if(!(all_valid_9 & (!_S3048)))
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
            _S3048 = true;
        }
        else
        {
            _S3048 = xmin_5 >= float(image_width_15);
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = ymax_5 <= 0.0f;
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            _S3048 = ymin_5 >= float(image_height_15);
        }
        if(_S3048)
        {
            _S3048 = true;
        }
        else
        {
            if(_S3047 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3048 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3048 = false;
                }
                if(_S3048)
                {
                    _S3048 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3048 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3048 = false;
                    }
                }
            }
            else
            {
                _S3048 = false;
            }
        }
        if(_S3048)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3064 = mean_20 - - mul_0(transpose_0(R_21), t_18);
        float _S3065 = _S3064.x;
        float _S3066 = _S3064.y;
        float _S3067 = _S3064.z;
        float norm_10 = (F32_sqrt((_S3065 * _S3065 + _S3066 * _S3066 + _S3067 * _S3067)));
        float x_44 = _S3065 / norm_10;
        float y_18 = _S3066 / norm_10;
        float z_15 = _S3067 / norm_10;
        float z2_37 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_44 * x_44 - y_18 * y_18;
        float fS1_15 = 2.0f * x_44 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_37 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_44) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_37 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_44) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_37 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_44) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3068 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3068);
        float3  _S3069 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3070 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3069 + _S3070 + make_float3 (0.5f), _S3068);
        (*rgb_12)[int(2)] = max_0(_S3069 - _S3070 + make_float3 (0.5f), _S3068);
        float3  _S3071 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3071 * make_float3 (float(- (F32_sign((dot_0(_S3071, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_22, float3  t_19, float fx_25, float fy_25, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    bool _S3072;
    bool _S3073;
    bool _S3074;
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_22, mean_21) + t_19;
        float _S3075 = length_1(mean_c_16);
        bool _S3076;
        if(_S3075 < near_plane_11)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3075 > far_plane_11;
        }
        if(_S3076)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3077 = scale_21.x;
        float sx_2 = (F32_exp((_S3077)));
        float _S3078 = scale_21.y;
        float sy_2 = (F32_exp((_S3078)));
        float sz_2 = scale_21.z - 0.5f * (_S3077 + _S3078);
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
        Matrix<float, 3, 3>  _S3079 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_38), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_38), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
        float3  vert0_c_1 = mul_0(R_22, mul_0(_S3079, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_19;
        float3  vert1_c_1 = mul_0(R_22, mul_0(_S3079, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_19;
        float3  vert2_c_1 = mul_0(R_22, mul_0(_S3079, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_19;
        float _S3080 = length_1(vert0_c_1);
        float _S3081 = length_1(vert1_c_1);
        float _S3082 = length_1(vert2_c_1);
        if(_S3080 < near_plane_11)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3080 > far_plane_11;
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3081 < near_plane_11;
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3081 > far_plane_11;
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3082 < near_plane_11;
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = _S3082 > far_plane_11;
        }
        if(_S3076)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float k_7;
        for(;;)
        {
            float2  _S3083 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_24 = length_0(_S3083);
            float _S3084 = vert0_c_1.z;
            float theta_20 = (F32_atan2((r_24), (_S3084)));
            if(theta_20 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_20 * theta_20 / 3.0f) / _S3084;
            }
            else
            {
                k_7 = theta_20 / r_24;
            }
            float2  _S3085 = _S3083 * make_float2 (k_7);
            *uv0_1 = _S3085;
            Matrix<float, 2, 2>  _S3086 = camera_distortion_jac_0(_S3085, dist_coeffs_29);
            bool _S3087 = !((F32_min((determinant_0(_S3086)), ((F32_min((_S3086.rows[int(0)].x), (_S3086.rows[int(1)].y)))))) > 0.0f);
            _S3072 = _S3087;
            if(_S3087)
            {
                break;
            }
            float u_47 = (*uv0_1).x;
            float v_47 = (*uv0_1).y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S3088 = *uv0_1 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_29)[int(0)] + r2_47 * ((*dist_coeffs_29)[int(1)] + r2_47 * ((*dist_coeffs_29)[int(2)] + r2_47 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_47 * v_47 + (*dist_coeffs_29)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_29)[int(6)] * r2_47, 2.0f * (*dist_coeffs_29)[int(5)] * u_47 * v_47 + (*dist_coeffs_29)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_29)[int(7)] * r2_47);
            float2  _S3089 = _S3088 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3088.x + (*dist_coeffs_29)[int(9)] * _S3088.y, 0.0f);
            *uv0_1 = make_float2 (fx_25 * _S3089.x + cx_20, fy_25 * _S3089.y + cy_20);
            break;
        }
        bool all_valid_10 = true & (!_S3072);
        for(;;)
        {
            float2  _S3090 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_25 = length_0(_S3090);
            float _S3091 = vert1_c_1.z;
            float theta_21 = (F32_atan2((r_25), (_S3091)));
            if(theta_21 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_21 * theta_21 / 3.0f) / _S3091;
            }
            else
            {
                k_7 = theta_21 / r_25;
            }
            float2  _S3092 = _S3090 * make_float2 (k_7);
            *uv1_1 = _S3092;
            Matrix<float, 2, 2>  _S3093 = camera_distortion_jac_0(_S3092, dist_coeffs_29);
            bool _S3094 = !((F32_min((determinant_0(_S3093)), ((F32_min((_S3093.rows[int(0)].x), (_S3093.rows[int(1)].y)))))) > 0.0f);
            _S3073 = _S3094;
            if(_S3094)
            {
                break;
            }
            float u_48 = (*uv1_1).x;
            float v_48 = (*uv1_1).y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S3095 = *uv1_1 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_29)[int(0)] + r2_48 * ((*dist_coeffs_29)[int(1)] + r2_48 * ((*dist_coeffs_29)[int(2)] + r2_48 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_48 * v_48 + (*dist_coeffs_29)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_29)[int(6)] * r2_48, 2.0f * (*dist_coeffs_29)[int(5)] * u_48 * v_48 + (*dist_coeffs_29)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_29)[int(7)] * r2_48);
            float2  _S3096 = _S3095 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3095.x + (*dist_coeffs_29)[int(9)] * _S3095.y, 0.0f);
            *uv1_1 = make_float2 (fx_25 * _S3096.x + cx_20, fy_25 * _S3096.y + cy_20);
            break;
        }
        bool all_valid_11 = all_valid_10 & (!_S3073);
        for(;;)
        {
            float2  _S3097 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_26 = length_0(_S3097);
            float _S3098 = vert2_c_1.z;
            float theta_22 = (F32_atan2((r_26), (_S3098)));
            if(theta_22 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_22 * theta_22 / 3.0f) / _S3098;
            }
            else
            {
                k_7 = theta_22 / r_26;
            }
            float2  _S3099 = _S3097 * make_float2 (k_7);
            *uv2_1 = _S3099;
            Matrix<float, 2, 2>  _S3100 = camera_distortion_jac_0(_S3099, dist_coeffs_29);
            bool _S3101 = !((F32_min((determinant_0(_S3100)), ((F32_min((_S3100.rows[int(0)].x), (_S3100.rows[int(1)].y)))))) > 0.0f);
            _S3074 = _S3101;
            if(_S3101)
            {
                break;
            }
            float u_49 = (*uv2_1).x;
            float v_49 = (*uv2_1).y;
            float r2_49 = u_49 * u_49 + v_49 * v_49;
            float2  _S3102 = *uv2_1 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_29)[int(0)] + r2_49 * ((*dist_coeffs_29)[int(1)] + r2_49 * ((*dist_coeffs_29)[int(2)] + r2_49 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_49 * v_49 + (*dist_coeffs_29)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_29)[int(6)] * r2_49, 2.0f * (*dist_coeffs_29)[int(5)] * u_49 * v_49 + (*dist_coeffs_29)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_29)[int(7)] * r2_49);
            float2  _S3103 = _S3102 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3102.x + (*dist_coeffs_29)[int(9)] * _S3102.y, 0.0f);
            *uv2_1 = make_float2 (fx_25 * _S3103.x + cx_20, fy_25 * _S3103.y + cy_20);
            break;
        }
        if(!(all_valid_11 & (!_S3074)))
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
            _S3076 = true;
        }
        else
        {
            _S3076 = xmin_6 >= float(image_width_16);
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = ymax_6 <= 0.0f;
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            _S3076 = ymin_6 >= float(image_height_16);
        }
        if(_S3076)
        {
            _S3076 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3076 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3076 = false;
                }
                if(_S3076)
                {
                    _S3076 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3076 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3076 = false;
                    }
                }
            }
            else
            {
                _S3076 = false;
            }
        }
        if(_S3076)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3080, _S3081, _S3082) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3104 = mean_21 - - mul_0(transpose_0(R_22), t_19);
        float _S3105 = _S3104.x;
        float _S3106 = _S3104.y;
        float _S3107 = _S3104.z;
        float norm_11 = (F32_sqrt((_S3105 * _S3105 + _S3106 * _S3106 + _S3107 * _S3107)));
        float x_46 = _S3105 / norm_11;
        float y_19 = _S3106 / norm_11;
        float z_16 = _S3107 / norm_11;
        float z2_39 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_46 * x_46 - y_19 * y_19;
        float fS1_16 = 2.0f * x_46 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_39 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_46) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_39 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_46) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_39 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_46) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3108 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3108);
        float3  _S3109 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3110 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3109 + _S3110 + make_float3 (0.5f), _S3108);
        (*rgb_13)[int(2)] = max_0(_S3109 - _S3110 + make_float3 (0.5f), _S3108);
        float3  _S3111 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3111 * make_float3 (float(- (F32_sign((dot_0(_S3111, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_23, float3  t_20, float fx_26, float fy_26, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_23, mean_22) + t_20;
    float _S3112 = scale_22.x;
    float sx_3 = (F32_exp((_S3112)));
    float _S3113 = scale_22.y;
    float sy_3 = (F32_exp((_S3113)));
    float sz_3 = scale_22.z - 0.5f * (_S3112 + _S3113);
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
    Matrix<float, 3, 3>  _S3114 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_40), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_40), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_2 = mul_0(R_23, mul_0(_S3114, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_20;
    float3  vert1_c_2 = mul_0(R_23, mul_0(_S3114, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_20;
    float3  vert2_c_2 = mul_0(R_23, mul_0(_S3114, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_20;
    float2  _S3115 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_50 = _S3115.x;
    float v_50 = _S3115.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S3116 = 2.0f * (*dist_coeffs_30)[int(4)];
    float _S3117 = 2.0f * (*dist_coeffs_30)[int(5)];
    float2  _S3118 = _S3115 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_30)[int(0)] + r2_50 * ((*dist_coeffs_30)[int(1)] + r2_50 * ((*dist_coeffs_30)[int(2)] + r2_50 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3116 * u_50 * v_50 + (*dist_coeffs_30)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_30)[int(6)] * r2_50, _S3117 * u_50 * v_50 + (*dist_coeffs_30)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_30)[int(7)] * r2_50);
    float2  _S3119 = _S3118 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3118.x + (*dist_coeffs_30)[int(9)] * _S3118.y, 0.0f);
    *uv0_2 = make_float2 (fx_26 * _S3119.x + cx_21, fy_26 * _S3119.y + cy_21);
    float2  _S3120 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_51 = _S3120.x;
    float v_51 = _S3120.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S3121 = _S3120 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_30)[int(0)] + r2_51 * ((*dist_coeffs_30)[int(1)] + r2_51 * ((*dist_coeffs_30)[int(2)] + r2_51 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3116 * u_51 * v_51 + (*dist_coeffs_30)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_30)[int(6)] * r2_51, _S3117 * u_51 * v_51 + (*dist_coeffs_30)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_30)[int(7)] * r2_51);
    float2  _S3122 = _S3121 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3121.x + (*dist_coeffs_30)[int(9)] * _S3121.y, 0.0f);
    *uv1_2 = make_float2 (fx_26 * _S3122.x + cx_21, fy_26 * _S3122.y + cy_21);
    float2  _S3123 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_52 = _S3123.x;
    float v_52 = _S3123.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float2  _S3124 = _S3123 * make_float2 (1.0f + r2_52 * ((*dist_coeffs_30)[int(0)] + r2_52 * ((*dist_coeffs_30)[int(1)] + r2_52 * ((*dist_coeffs_30)[int(2)] + r2_52 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3116 * u_52 * v_52 + (*dist_coeffs_30)[int(5)] * (r2_52 + 2.0f * u_52 * u_52) + (*dist_coeffs_30)[int(6)] * r2_52, _S3117 * u_52 * v_52 + (*dist_coeffs_30)[int(4)] * (r2_52 + 2.0f * v_52 * v_52) + (*dist_coeffs_30)[int(7)] * r2_52);
    float2  _S3125 = _S3124 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3124.x + (*dist_coeffs_30)[int(9)] * _S3124.y, 0.0f);
    float _S3126 = fx_26 * _S3125.x + cx_21;
    float _S3127 = fy_26 * _S3125.y + cy_21;
    float2  _S3128 = make_float2 (_S3126, _S3127);
    *uv2_2 = _S3128;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S3128 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3128)));
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3126))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3127))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3126))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3127))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S3129 = mean_22 - - mul_0(transpose_0(R_23), t_20);
    float _S3130 = _S3129.x;
    float _S3131 = _S3129.y;
    float _S3132 = _S3129.z;
    float norm_12 = (F32_sqrt((_S3130 * _S3130 + _S3131 * _S3131 + _S3132 * _S3132)));
    float x_48 = _S3130 / norm_12;
    float y_20 = _S3131 / norm_12;
    float z_17 = _S3132 / norm_12;
    float z2_41 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_48 * x_48 - y_20 * y_20;
    float fS1_17 = 2.0f * x_48 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_41 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_48) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_41 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_48) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_41 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_48) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S3133 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3133);
    float3  _S3134 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S3135 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S3134 + _S3135 + make_float3 (0.5f), _S3133);
    (*rgb_14)[int(2)] = max_0(_S3134 - _S3135 + make_float3 (0.5f), _S3133);
    float3  _S3136 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S3136 * make_float3 (float(- (F32_sign((dot_0(_S3136, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_24, float3  t_21, float fx_27, float fy_27, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_24, mean_23) + t_21;
    float _S3137 = scale_23.x;
    float sx_4 = (F32_exp((_S3137)));
    float _S3138 = scale_23.y;
    float sy_4 = (F32_exp((_S3138)));
    float sz_4 = scale_23.z - 0.5f * (_S3137 + _S3138);
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
    Matrix<float, 3, 3>  _S3139 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_42), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_42), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  vert0_c_3 = mul_0(R_24, mul_0(_S3139, make_float3 (sx_4, 0.0f, 0.0f)) + mean_23) + t_21;
    float3  vert1_c_3 = mul_0(R_24, mul_0(_S3139, make_float3 (sx_4 * (-0.5f + sz_4), sy_4, 0.0f)) + mean_23) + t_21;
    float3  vert2_c_3 = mul_0(R_24, mul_0(_S3139, make_float3 (sx_4 * (-0.5f - sz_4), - sy_4, 0.0f)) + mean_23) + t_21;
    float2  _S3140 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_27 = length_0(_S3140);
    float _S3141 = vert0_c_3.z;
    float theta_23 = (F32_atan2((r_27), (_S3141)));
    float k_8;
    if(theta_23 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_23 * theta_23 / 3.0f) / _S3141;
    }
    else
    {
        k_8 = theta_23 / r_27;
    }
    float2  _S3142 = _S3140 * make_float2 (k_8);
    float u_53 = _S3142.x;
    float v_53 = _S3142.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S3143 = 2.0f * (*dist_coeffs_31)[int(4)];
    float _S3144 = 2.0f * (*dist_coeffs_31)[int(5)];
    float2  _S3145 = _S3142 * make_float2 (1.0f + r2_53 * ((*dist_coeffs_31)[int(0)] + r2_53 * ((*dist_coeffs_31)[int(1)] + r2_53 * ((*dist_coeffs_31)[int(2)] + r2_53 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3143 * u_53 * v_53 + (*dist_coeffs_31)[int(5)] * (r2_53 + 2.0f * u_53 * u_53) + (*dist_coeffs_31)[int(6)] * r2_53, _S3144 * u_53 * v_53 + (*dist_coeffs_31)[int(4)] * (r2_53 + 2.0f * v_53 * v_53) + (*dist_coeffs_31)[int(7)] * r2_53);
    float2  _S3146 = _S3145 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3145.x + (*dist_coeffs_31)[int(9)] * _S3145.y, 0.0f);
    *uv0_3 = make_float2 (fx_27 * _S3146.x + cx_22, fy_27 * _S3146.y + cy_22);
    float2  _S3147 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_28 = length_0(_S3147);
    float _S3148 = vert1_c_3.z;
    float theta_24 = (F32_atan2((r_28), (_S3148)));
    if(theta_24 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_24 * theta_24 / 3.0f) / _S3148;
    }
    else
    {
        k_8 = theta_24 / r_28;
    }
    float2  _S3149 = _S3147 * make_float2 (k_8);
    float u_54 = _S3149.x;
    float v_54 = _S3149.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float2  _S3150 = _S3149 * make_float2 (1.0f + r2_54 * ((*dist_coeffs_31)[int(0)] + r2_54 * ((*dist_coeffs_31)[int(1)] + r2_54 * ((*dist_coeffs_31)[int(2)] + r2_54 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3143 * u_54 * v_54 + (*dist_coeffs_31)[int(5)] * (r2_54 + 2.0f * u_54 * u_54) + (*dist_coeffs_31)[int(6)] * r2_54, _S3144 * u_54 * v_54 + (*dist_coeffs_31)[int(4)] * (r2_54 + 2.0f * v_54 * v_54) + (*dist_coeffs_31)[int(7)] * r2_54);
    float2  _S3151 = _S3150 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3150.x + (*dist_coeffs_31)[int(9)] * _S3150.y, 0.0f);
    *uv1_3 = make_float2 (fx_27 * _S3151.x + cx_22, fy_27 * _S3151.y + cy_22);
    float2  _S3152 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_29 = length_0(_S3152);
    float _S3153 = vert2_c_3.z;
    float theta_25 = (F32_atan2((r_29), (_S3153)));
    if(theta_25 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_25 * theta_25 / 3.0f) / _S3153;
    }
    else
    {
        k_8 = theta_25 / r_29;
    }
    float2  _S3154 = _S3152 * make_float2 (k_8);
    float u_55 = _S3154.x;
    float v_55 = _S3154.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float2  _S3155 = _S3154 * make_float2 (1.0f + r2_55 * ((*dist_coeffs_31)[int(0)] + r2_55 * ((*dist_coeffs_31)[int(1)] + r2_55 * ((*dist_coeffs_31)[int(2)] + r2_55 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3143 * u_55 * v_55 + (*dist_coeffs_31)[int(5)] * (r2_55 + 2.0f * u_55 * u_55) + (*dist_coeffs_31)[int(6)] * r2_55, _S3144 * u_55 * v_55 + (*dist_coeffs_31)[int(4)] * (r2_55 + 2.0f * v_55 * v_55) + (*dist_coeffs_31)[int(7)] * r2_55);
    float2  _S3156 = _S3155 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3155.x + (*dist_coeffs_31)[int(9)] * _S3155.y, 0.0f);
    float _S3157 = fx_27 * _S3156.x + cx_22;
    float _S3158 = fy_27 * _S3156.y + cy_22;
    float2  _S3159 = make_float2 (_S3157, _S3158);
    *uv2_3 = _S3159;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3159 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3159)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3157))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3158))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3157))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3158))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3160 = mean_23 - - mul_0(transpose_0(R_24), t_21);
    float _S3161 = _S3160.x;
    float _S3162 = _S3160.y;
    float _S3163 = _S3160.z;
    float norm_13 = (F32_sqrt((_S3161 * _S3161 + _S3162 * _S3162 + _S3163 * _S3163)));
    float x_50 = _S3161 / norm_13;
    float y_21 = _S3162 / norm_13;
    float z_18 = _S3163 / norm_13;
    float z2_43 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_50 * x_50 - y_21 * y_21;
    float fS1_18 = 2.0f * x_50 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_43 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_50) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_43 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_50) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_50 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_43 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_50) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_50 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3164 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3164);
    float3  _S3165 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3166 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3165 + _S3166 + make_float3 (0.5f), _S3164);
    (*rgb_15)[int(2)] = max_0(_S3165 - _S3166 + make_float3 (0.5f), _S3164);
    float3  _S3167 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3167 * make_float3 (float(- (F32_sign((dot_0(_S3167, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  _s_dOut_9)
{
    float _S3168 = length_1((*dpx_16).primal_0);
    float3  _S3169 = (*dpx_16).primal_0 * _s_dOut_9;
    float3  _S3170 = make_float3 (1.0f / _S3168) * _s_dOut_9;
    float _S3171 = - ((_S3169.x + _S3169.y + _S3169.z) / (_S3168 * _S3168));
    float3  _S3172 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3173;
    (&_S3173)->primal_0 = (*dpx_16).primal_0;
    (&_S3173)->differential_0 = _S3172;
    s_bwd_length_impl_1(&_S3173, _S3171);
    float3  _S3174 = _S3170 + _S3173.differential_0;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S3174;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3175, float3  _S3176)
{
    s_bwd_prop_normalize_impl_0(_S3175, _S3176);
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3177, float3  _S3178)
{
    _d_log_vector_0(_S3177, _S3178);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3179, float _S3180)
{
    _d_exp2_0(_S3179, _S3180);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3181, float _S3182)
{
    _d_abs_0(_S3181, _S3182);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_25, float3  t_22, float fx_28, float fy_28, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_0(R_25, mean_24) + t_22;
    float _S3183 = scale_24.x;
    float _S3184 = s_primal_ctx_exp_1(_S3183);
    float _S3185 = scale_24.y;
    float _S3186 = s_primal_ctx_exp_1(_S3185);
    float sz_5 = scale_24.z - 0.5f * (_S3183 + _S3185);
    float _S3187 = quat_25.y;
    float x2_25 = _S3187 * _S3187;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_44 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3188 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_44), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_44), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3189 = make_float3 (_S3184, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_0(_S3188, _S3189) + mean_24;
    float _S3190 = -0.5f + sz_5;
    float3  _S3191 = make_float3 (_S3184 * _S3190, _S3186, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_0(_S3188, _S3191) + mean_24;
    float _S3192 = -0.5f - sz_5;
    float3  _S3193 = make_float3 (_S3184 * _S3192, - _S3186, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_0(_S3188, _S3193) + mean_24;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_25, vert0_1) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_25, vert1_1) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_25, vert2_1) + t_22;
    float2  _S3194 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3195 = vert0_c_4.z;
    float2  _S3196 = make_float2 (_S3195);
    float2  _S3197 = _S3194 / make_float2 (_S3195);
    float2  _S3198 = make_float2 (_S3195 * _S3195);
    float u_56 = _S3197.x;
    float v_56 = _S3197.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S3199 = (*dist_coeffs_32)[int(2)] + r2_56 * (*dist_coeffs_32)[int(3)];
    float _S3200 = (*dist_coeffs_32)[int(1)] + r2_56 * _S3199;
    float _S3201 = (*dist_coeffs_32)[int(0)] + r2_56 * _S3200;
    float radial_1 = 1.0f + r2_56 * _S3201;
    float2  _S3202 = make_float2 (radial_1);
    float _S3203 = 2.0f * (*dist_coeffs_32)[int(4)];
    float _S3204 = _S3203 * u_56;
    float _S3205 = 2.0f * u_56;
    float _S3206 = 2.0f * (*dist_coeffs_32)[int(5)];
    float _S3207 = _S3206 * u_56;
    float _S3208 = 2.0f * v_56;
    float2  _S3209 = _S3197 * make_float2 (radial_1) + make_float2 (_S3204 * v_56 + (*dist_coeffs_32)[int(5)] * (r2_56 + _S3205 * u_56) + (*dist_coeffs_32)[int(6)] * r2_56, _S3207 * v_56 + (*dist_coeffs_32)[int(4)] * (r2_56 + _S3208 * v_56) + (*dist_coeffs_32)[int(7)] * r2_56);
    float2  _S3210 = _S3209 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3209.x + (*dist_coeffs_32)[int(9)] * _S3209.y, 0.0f);
    float _S3211 = fx_28 * _S3210.x + cx_23;
    float _S3212 = fy_28 * _S3210.y + cy_23;
    float2  _S3213 = make_float2 (_S3211, _S3212);
    float2  _S3214 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3215 = vert1_c_4.z;
    float2  _S3216 = make_float2 (_S3215);
    float2  _S3217 = _S3214 / make_float2 (_S3215);
    float2  _S3218 = make_float2 (_S3215 * _S3215);
    float u_57 = _S3217.x;
    float v_57 = _S3217.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S3219 = (*dist_coeffs_32)[int(2)] + r2_57 * (*dist_coeffs_32)[int(3)];
    float _S3220 = (*dist_coeffs_32)[int(1)] + r2_57 * _S3219;
    float _S3221 = (*dist_coeffs_32)[int(0)] + r2_57 * _S3220;
    float radial_2 = 1.0f + r2_57 * _S3221;
    float2  _S3222 = make_float2 (radial_2);
    float _S3223 = _S3203 * u_57;
    float _S3224 = 2.0f * u_57;
    float _S3225 = _S3206 * u_57;
    float _S3226 = 2.0f * v_57;
    float2  _S3227 = _S3217 * make_float2 (radial_2) + make_float2 (_S3223 * v_57 + (*dist_coeffs_32)[int(5)] * (r2_57 + _S3224 * u_57) + (*dist_coeffs_32)[int(6)] * r2_57, _S3225 * v_57 + (*dist_coeffs_32)[int(4)] * (r2_57 + _S3226 * v_57) + (*dist_coeffs_32)[int(7)] * r2_57);
    float2  _S3228 = _S3227 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3227.x + (*dist_coeffs_32)[int(9)] * _S3227.y, 0.0f);
    float _S3229 = fx_28 * _S3228.x + cx_23;
    float _S3230 = fy_28 * _S3228.y + cy_23;
    float2  _S3231 = make_float2 (_S3229, _S3230);
    float2  _S3232 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3233 = vert2_c_4.z;
    float2  _S3234 = make_float2 (_S3233);
    float2  _S3235 = _S3232 / make_float2 (_S3233);
    float2  _S3236 = make_float2 (_S3233 * _S3233);
    float u_58 = _S3235.x;
    float v_58 = _S3235.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S3237 = (*dist_coeffs_32)[int(2)] + r2_58 * (*dist_coeffs_32)[int(3)];
    float _S3238 = (*dist_coeffs_32)[int(1)] + r2_58 * _S3237;
    float _S3239 = (*dist_coeffs_32)[int(0)] + r2_58 * _S3238;
    float radial_3 = 1.0f + r2_58 * _S3239;
    float2  _S3240 = make_float2 (radial_3);
    float _S3241 = _S3203 * u_58;
    float _S3242 = 2.0f * u_58;
    float _S3243 = _S3206 * u_58;
    float _S3244 = 2.0f * v_58;
    float2  _S3245 = _S3235 * make_float2 (radial_3) + make_float2 (_S3241 * v_58 + (*dist_coeffs_32)[int(5)] * (r2_58 + _S3242 * u_58) + (*dist_coeffs_32)[int(6)] * r2_58, _S3243 * v_58 + (*dist_coeffs_32)[int(4)] * (r2_58 + _S3244 * v_58) + (*dist_coeffs_32)[int(7)] * r2_58);
    float2  _S3246 = _S3245 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3245.x + (*dist_coeffs_32)[int(9)] * _S3245.y, 0.0f);
    float _S3247 = fx_28 * _S3246.x + cx_23;
    float _S3248 = fy_28 * _S3246.y + cy_23;
    float2  _S3249 = make_float2 (_S3247, _S3248);
    float2  e0_4 = _S3231 - _S3213;
    float2  e1_4 = _S3249 - _S3231;
    float2  e2_0 = _S3213 - _S3249;
    float _S3250 = e0_4.x;
    float _S3251 = e1_4.y;
    float _S3252 = e0_4.y;
    float _S3253 = e1_4.x;
    float _S3254 = _S3250 * _S3251 - _S3252 * _S3253;
    float _S3255 = 1.0f - hardness_4.y;
    float _S3256 = -1.0f / _S3255;
    float _S3257 = _S3255 * _S3255;
    float _S3258 = s_primal_ctx_max_0(_S3211, _S3229);
    float _S3259 = s_primal_ctx_min_0(_S3211, _S3229);
    float _S3260 = s_primal_ctx_max_0(_S3212, _S3230);
    float _S3261 = s_primal_ctx_min_0(_S3212, _S3230);
    float3  _S3262 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3263 = transpose_0(R_25);
    float3  _S3264 = mean_24 - - s_primal_ctx_mul_0(_S3263, t_22);
    float _S3265 = _S3264.x;
    float _S3266 = _S3264.y;
    float _S3267 = _S3264.z;
    float _S3268 = _S3265 * _S3265 + _S3266 * _S3266 + _S3267 * _S3267;
    float _S3269 = s_primal_ctx_sqrt_0(_S3268);
    float x_51 = _S3265 / _S3269;
    float3  _S3270 = make_float3 (x_51);
    float _S3271 = _S3269 * _S3269;
    float y_22 = _S3266 / _S3269;
    float z_19 = _S3267 / _S3269;
    float3  _S3272 = make_float3 (z_19);
    float _S3273 = - y_22;
    float3  _S3274 = make_float3 (_S3273);
    float z2_45 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_51 * x_51 - y_22 * y_22;
    float _S3275 = 2.0f * x_51;
    float fS1_19 = _S3275 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_45 - 0.31539157032966614f;
    float3  _S3276 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_51;
    float3  _S3277 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3278 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3279 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3280 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_45 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3281 = 1.86588168144226074f * z2_45 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3281;
    float3  _S3282 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_51;
    float3  _S3283 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3284 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3285 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3286 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_51 * fC1_19 - y_22 * fS1_19);
    float3  _S3287 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_51 * fS1_19 + y_22 * fC1_19);
    float3  _S3288 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3273) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_51) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3289 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3290 = make_float3 (0.0f);
    float3  _S3291 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3292 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3293 = make_float3 (_S3292);
    float3  _S3294 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3292);
    float3  _S3295 = _S3291 + _S3294 + make_float3 (0.5f);
    float3  _S3296 = _S3291 - _S3294 + make_float3 (0.5f);
    float3  _S3297 = vert1_c_4 - vert0_c_4;
    float3  _S3298 = vert2_c_4 - vert0_c_4;
    float3  _S3299 = s_primal_ctx_cross_0(_S3297, _S3298);
    float3  _S3300 = normalize_0(_S3299);
    float3  _S3301 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3300, mean_c_19)))))) * v_normal_0;
    float3  _S3302 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3303;
    (&_S3303)->primal_0 = _S3300;
    (&_S3303)->differential_0 = _S3302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3304;
    (&_S3304)->primal_0 = mean_c_19;
    (&_S3304)->differential_0 = _S3302;
    s_bwd_prop_dot_0(&_S3303, &_S3304, 0.0f);
    float3  _S3305 = _S3301 + _S3303.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3306;
    (&_S3306)->primal_0 = _S3299;
    (&_S3306)->differential_0 = _S3302;
    s_bwd_normalize_impl_0(&_S3306, _S3305);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3307;
    (&_S3307)->primal_0 = _S3297;
    (&_S3307)->differential_0 = _S3302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3308;
    (&_S3308)->primal_0 = _S3298;
    (&_S3308)->differential_0 = _S3302;
    s_bwd_prop_cross_0(&_S3307, &_S3308, _S3306.differential_0);
    float3  _S3309 = - _S3308.differential_0;
    float3  _S3310 = - _S3307.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3311;
    (&_S3311)->primal_0 = _S3296;
    (&_S3311)->differential_0 = _S3302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3312;
    (&_S3312)->primal_0 = _S3290;
    (&_S3312)->differential_0 = _S3302;
    s_bwd_prop_max_0(&_S3311, &_S3312, (*v_rgb_6)[int(2)]);
    float3  _S3313 = - _S3311.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3314;
    (&_S3314)->primal_0 = _S3295;
    (&_S3314)->differential_0 = _S3302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3315;
    (&_S3315)->primal_0 = _S3290;
    (&_S3315)->differential_0 = _S3302;
    s_bwd_prop_max_0(&_S3314, &_S3315, (*v_rgb_6)[int(1)]);
    float3  _S3316 = _S3293 * (_S3313 + _S3314.differential_0);
    float3  _S3317 = _S3311.differential_0 + _S3314.differential_0;
    float3  _S3318 = make_float3 (0.5f) * - _S3317;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3319;
    (&_S3319)->primal_0 = _S3289;
    (&_S3319)->differential_0 = _S3302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3320;
    (&_S3320)->primal_0 = _S3290;
    (&_S3320)->differential_0 = _S3302;
    s_bwd_prop_max_0(&_S3319, &_S3320, (*v_rgb_6)[int(0)]);
    float3  _S3321 = _S3318 + _S3319.differential_0;
    float3  _S3322 = _S3317 + _S3319.differential_0;
    float3  _S3323 = _S3287 * _S3322;
    float3  _S3324 = (*sh_coeffs_19)[int(15)] * _S3322;
    float3  _S3325 = _S3285 * _S3322;
    float3  _S3326 = (*sh_coeffs_19)[int(14)] * _S3322;
    float3  _S3327 = _S3283 * _S3322;
    float3  _S3328 = (*sh_coeffs_19)[int(13)] * _S3322;
    float3  _S3329 = _S3282 * _S3322;
    float3  _S3330 = (*sh_coeffs_19)[int(12)] * _S3322;
    float3  _S3331 = _S3284 * _S3322;
    float3  _S3332 = (*sh_coeffs_19)[int(11)] * _S3322;
    float3  _S3333 = _S3286 * _S3322;
    float3  _S3334 = (*sh_coeffs_19)[int(10)] * _S3322;
    float3  _S3335 = _S3288 * _S3322;
    float3  _S3336 = (*sh_coeffs_19)[int(9)] * _S3322;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3336.x + _S3336.y + _S3336.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3324.x + _S3324.y + _S3324.z);
    float _S3337 = _S3334.x + _S3334.y + _S3334.z;
    float _S3338 = _S3326.x + _S3326.y + _S3326.z;
    float _S3339 = _S3332.x + _S3332.y + _S3332.z;
    float _S3340 = _S3328.x + _S3328.y + _S3328.z;
    float _S3341 = _S3330.x + _S3330.y + _S3330.z;
    float _S3342 = - s_diff_fC2_T_5;
    float3  _S3343 = _S3279 * _S3322;
    float3  _S3344 = (*sh_coeffs_19)[int(8)] * _S3322;
    float3  _S3345 = _S3277 * _S3322;
    float3  _S3346 = (*sh_coeffs_19)[int(7)] * _S3322;
    float3  _S3347 = _S3276 * _S3322;
    float3  _S3348 = (*sh_coeffs_19)[int(6)] * _S3322;
    float3  _S3349 = _S3278 * _S3322;
    float3  _S3350 = (*sh_coeffs_19)[int(5)] * _S3322;
    float3  _S3351 = _S3280 * _S3322;
    float3  _S3352 = (*sh_coeffs_19)[int(4)] * _S3322;
    float _S3353 = _S3350.x + _S3350.y + _S3350.z;
    float _S3354 = _S3346.x + _S3346.y + _S3346.z;
    float _S3355 = fTmp1B_19 * _S3337 + x_51 * s_diff_fS2_T_5 + y_22 * _S3342 + 0.54627424478530884f * (_S3352.x + _S3352.y + _S3352.z);
    float _S3356 = fTmp1B_19 * _S3338 + y_22 * s_diff_fS2_T_5 + x_51 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3344.x + _S3344.y + _S3344.z);
    float _S3357 = y_22 * - _S3356;
    float _S3358 = x_51 * _S3356;
    float _S3359 = z_19 * (1.86588168144226074f * (z_19 * _S3341) + -2.28522896766662598f * (y_22 * _S3339 + x_51 * _S3340) + 0.94617468118667603f * (_S3348.x + _S3348.y + _S3348.z));
    float3  _S3360 = make_float3 (0.48860251903533936f) * _S3322;
    float3  _S3361 = - _S3360;
    float3  _S3362 = _S3270 * _S3361;
    float3  _S3363 = (*sh_coeffs_19)[int(3)] * _S3361;
    float3  _S3364 = _S3272 * _S3360;
    float3  _S3365 = (*sh_coeffs_19)[int(2)] * _S3360;
    float3  _S3366 = _S3274 * _S3360;
    float3  _S3367 = (*sh_coeffs_19)[int(1)] * _S3360;
    float _S3368 = (_S3281 * _S3341 + 1.44530570507049561f * (fS1_19 * _S3337 + fC1_19 * _S3338) + -1.09254848957061768f * (y_22 * _S3353 + x_51 * _S3354) + _S3359 + _S3359 + _S3365.x + _S3365.y + _S3365.z) / _S3271;
    float _S3369 = _S3269 * _S3368;
    float _S3370 = (fTmp0C_19 * _S3339 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3342 + fTmp0B_19 * _S3353 + _S3275 * _S3355 + _S3357 + _S3357 + - (_S3367.x + _S3367.y + _S3367.z)) / _S3271;
    float _S3371 = _S3269 * _S3370;
    float _S3372 = (fTmp0C_19 * _S3340 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3354 + 2.0f * (y_22 * _S3355) + _S3358 + _S3358 + _S3363.x + _S3363.y + _S3363.z) / _S3271;
    float _S3373 = _S3269 * _S3372;
    float _S3374 = _S3267 * - _S3368 + _S3266 * - _S3370 + _S3265 * - _S3372;
    DiffPair_float_0 _S3375;
    (&_S3375)->primal_0 = _S3268;
    (&_S3375)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3375, _S3374);
    float _S3376 = _S3267 * _S3375.differential_0;
    float _S3377 = _S3266 * _S3375.differential_0;
    float _S3378 = _S3265 * _S3375.differential_0;
    float3  _S3379 = make_float3 (0.282094806432724f) * _S3322;
    float3  _S3380 = make_float3 (_S3373 + _S3378 + _S3378, _S3371 + _S3377 + _S3377, _S3369 + _S3376 + _S3376);
    float3  _S3381 = - - _S3380;
    Matrix<float, 3, 3>  _S3382 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3383;
    (&_S3383)->primal_0 = _S3263;
    (&_S3383)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3384;
    (&_S3384)->primal_0 = t_22;
    (&_S3384)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3383, &_S3384, _S3381);
    Matrix<float, 3, 3>  _S3385 = transpose_0(_S3383.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3386;
    (&_S3386)->primal_0 = _S3262;
    (&_S3386)->differential_0 = _S3302;
    s_bwd_prop_log_1(&_S3386, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3387;
    (&_S3387)->primal_0 = vert2_c_4;
    (&_S3387)->differential_0 = _S3302;
    s_bwd_length_impl_1(&_S3387, _S3386.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3388;
    (&_S3388)->primal_0 = vert1_c_4;
    (&_S3388)->differential_0 = _S3302;
    s_bwd_length_impl_1(&_S3388, _S3386.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3389;
    (&_S3389)->primal_0 = vert0_c_4;
    (&_S3389)->differential_0 = _S3302;
    s_bwd_length_impl_1(&_S3389, _S3386.differential_0.x);
    DiffPair_float_0 _S3390;
    (&_S3390)->primal_0 = _S3261;
    (&_S3390)->differential_0 = 0.0f;
    DiffPair_float_0 _S3391;
    (&_S3391)->primal_0 = _S3248;
    (&_S3391)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3390, &_S3391, 0.0f);
    DiffPair_float_0 _S3392;
    (&_S3392)->primal_0 = _S3212;
    (&_S3392)->differential_0 = 0.0f;
    DiffPair_float_0 _S3393;
    (&_S3393)->primal_0 = _S3230;
    (&_S3393)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3392, &_S3393, _S3390.differential_0);
    DiffPair_float_0 _S3394;
    (&_S3394)->primal_0 = _S3260;
    (&_S3394)->differential_0 = 0.0f;
    DiffPair_float_0 _S3395;
    (&_S3395)->primal_0 = _S3248;
    (&_S3395)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3394, &_S3395, 0.0f);
    DiffPair_float_0 _S3396;
    (&_S3396)->primal_0 = _S3212;
    (&_S3396)->differential_0 = 0.0f;
    DiffPair_float_0 _S3397;
    (&_S3397)->primal_0 = _S3230;
    (&_S3397)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3396, &_S3397, _S3394.differential_0);
    DiffPair_float_0 _S3398;
    (&_S3398)->primal_0 = _S3259;
    (&_S3398)->differential_0 = 0.0f;
    DiffPair_float_0 _S3399;
    (&_S3399)->primal_0 = _S3247;
    (&_S3399)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3398, &_S3399, 0.0f);
    DiffPair_float_0 _S3400;
    (&_S3400)->primal_0 = _S3211;
    (&_S3400)->differential_0 = 0.0f;
    DiffPair_float_0 _S3401;
    (&_S3401)->primal_0 = _S3229;
    (&_S3401)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3400, &_S3401, _S3398.differential_0);
    DiffPair_float_0 _S3402;
    (&_S3402)->primal_0 = _S3258;
    (&_S3402)->differential_0 = 0.0f;
    DiffPair_float_0 _S3403;
    (&_S3403)->primal_0 = _S3247;
    (&_S3403)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3402, &_S3403, 0.0f);
    DiffPair_float_0 _S3404;
    (&_S3404)->primal_0 = _S3211;
    (&_S3404)->differential_0 = 0.0f;
    DiffPair_float_0 _S3405;
    (&_S3405)->primal_0 = _S3229;
    (&_S3405)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3404, &_S3405, _S3402.differential_0);
    DiffPair_float_0 _S3406;
    (&_S3406)->primal_0 = _S3256;
    (&_S3406)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3406, 0.0f);
    float _S3407 = - (-1.0f * - (_S3406.differential_0 / _S3257));
    float2  _S3408 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3409;
    (&_S3409)->primal_0 = e2_0;
    (&_S3409)->differential_0 = _S3408;
    s_bwd_length_impl_0(&_S3409, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3410;
    (&_S3410)->primal_0 = e1_4;
    (&_S3410)->differential_0 = _S3408;
    s_bwd_length_impl_0(&_S3410, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3411;
    (&_S3411)->primal_0 = e0_4;
    (&_S3411)->differential_0 = _S3408;
    s_bwd_length_impl_0(&_S3411, -0.0f);
    DiffPair_float_0 _S3412;
    (&_S3412)->primal_0 = _S3254;
    (&_S3412)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3412, 0.0f);
    float _S3413 = - _S3412.differential_0;
    float2  _S3414 = _S3410.differential_0 + make_float2 (_S3252 * _S3413, _S3250 * _S3412.differential_0);
    float2  _S3415 = _S3411.differential_0 + make_float2 (_S3251 * _S3412.differential_0, _S3253 * _S3413);
    float2  _S3416 = v_uv2_0 + - _S3409.differential_0 + _S3414;
    float _S3417 = fx_28 * (_S3399.differential_0 + _S3403.differential_0 + _S3416.x);
    float2  _S3418 = make_float2 (_S3417, fy_28 * (_S3391.differential_0 + _S3395.differential_0 + _S3416.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3417, (*dist_coeffs_32)[int(9)] * _S3417);
    float2  _S3419 = _S3235 * _S3418;
    float _S3420 = (*dist_coeffs_32)[int(4)] * _S3418.y;
    float _S3421 = (*dist_coeffs_32)[int(5)] * _S3418.x;
    float _S3422 = _S3419.x + _S3419.y;
    float _S3423 = r2_58 * _S3422;
    float _S3424 = r2_58 * _S3423;
    float _S3425 = (*dist_coeffs_32)[int(7)] * _S3418.y + _S3420 + (*dist_coeffs_32)[int(6)] * _S3418.x + _S3421 + _S3239 * _S3422 + _S3238 * _S3423 + _S3237 * _S3424 + (*dist_coeffs_32)[int(3)] * (r2_58 * _S3424);
    float _S3426 = v_58 * _S3425;
    float _S3427 = u_58 * _S3425;
    float2  _S3428 = (_S3240 * _S3418 + make_float2 (_S3206 * (v_58 * _S3418.y) + _S3242 * _S3421 + 2.0f * (u_58 * _S3421) + _S3203 * (v_58 * _S3418.x) + _S3427 + _S3427, _S3244 * _S3420 + 2.0f * (v_58 * _S3420) + _S3243 * _S3418.y + _S3241 * _S3418.x + _S3426 + _S3426)) / _S3236;
    float2  _S3429 = _S3232 * - _S3428;
    float2  _S3430 = _S3234 * _S3428;
    float2  _S3431 = v_uv1_0 + - _S3414 + _S3415;
    float _S3432 = fx_28 * (_S3401.differential_0 + _S3405.differential_0 + _S3431.x);
    float2  _S3433 = make_float2 (_S3432, fy_28 * (_S3393.differential_0 + _S3397.differential_0 + _S3431.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3432, (*dist_coeffs_32)[int(9)] * _S3432);
    float2  _S3434 = _S3217 * _S3433;
    float _S3435 = (*dist_coeffs_32)[int(4)] * _S3433.y;
    float _S3436 = (*dist_coeffs_32)[int(5)] * _S3433.x;
    float _S3437 = _S3434.x + _S3434.y;
    float _S3438 = r2_57 * _S3437;
    float _S3439 = r2_57 * _S3438;
    float _S3440 = (*dist_coeffs_32)[int(7)] * _S3433.y + _S3435 + (*dist_coeffs_32)[int(6)] * _S3433.x + _S3436 + _S3221 * _S3437 + _S3220 * _S3438 + _S3219 * _S3439 + (*dist_coeffs_32)[int(3)] * (r2_57 * _S3439);
    float _S3441 = v_57 * _S3440;
    float _S3442 = u_57 * _S3440;
    float2  _S3443 = (_S3222 * _S3433 + make_float2 (_S3206 * (v_57 * _S3433.y) + _S3224 * _S3436 + 2.0f * (u_57 * _S3436) + _S3203 * (v_57 * _S3433.x) + _S3442 + _S3442, _S3226 * _S3435 + 2.0f * (v_57 * _S3435) + _S3225 * _S3433.y + _S3223 * _S3433.x + _S3441 + _S3441)) / _S3218;
    float2  _S3444 = _S3214 * - _S3443;
    float2  _S3445 = _S3216 * _S3443;
    float _S3446 = _S3444.x + _S3444.y;
    float2  _S3447 = v_uv0_0 + _S3409.differential_0 + - _S3415;
    float _S3448 = fx_28 * (_S3400.differential_0 + _S3404.differential_0 + _S3447.x);
    float2  _S3449 = make_float2 (_S3448, fy_28 * (_S3392.differential_0 + _S3396.differential_0 + _S3447.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3448, (*dist_coeffs_32)[int(9)] * _S3448);
    float2  _S3450 = _S3197 * _S3449;
    float _S3451 = (*dist_coeffs_32)[int(4)] * _S3449.y;
    float _S3452 = (*dist_coeffs_32)[int(5)] * _S3449.x;
    float _S3453 = _S3450.x + _S3450.y;
    float _S3454 = r2_56 * _S3453;
    float _S3455 = r2_56 * _S3454;
    float _S3456 = (*dist_coeffs_32)[int(7)] * _S3449.y + _S3451 + (*dist_coeffs_32)[int(6)] * _S3449.x + _S3452 + _S3201 * _S3453 + _S3200 * _S3454 + _S3199 * _S3455 + (*dist_coeffs_32)[int(3)] * (r2_56 * _S3455);
    float _S3457 = v_56 * _S3456;
    float _S3458 = u_56 * _S3456;
    float2  _S3459 = (_S3202 * _S3449 + make_float2 (_S3206 * (v_56 * _S3449.y) + _S3205 * _S3452 + 2.0f * (u_56 * _S3452) + _S3203 * (v_56 * _S3449.x) + _S3458 + _S3458, _S3208 * _S3451 + 2.0f * (v_56 * _S3451) + _S3207 * _S3449.y + _S3204 * _S3449.x + _S3457 + _S3457)) / _S3198;
    float2  _S3460 = _S3194 * - _S3459;
    float2  _S3461 = _S3196 * _S3459;
    float _S3462 = _S3460.x + _S3460.y;
    float3  _S3463 = _S3308.differential_0 + _S3387.differential_0 + make_float3 (_S3430.x, _S3430.y, _S3429.x + _S3429.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3464;
    (&_S3464)->primal_0 = R_25;
    (&_S3464)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3465;
    (&_S3465)->primal_0 = vert2_1;
    (&_S3465)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3464, &_S3465, _S3463);
    float3  _S3466 = _S3307.differential_0 + _S3388.differential_0 + make_float3 (_S3445.x, _S3445.y, _S3446);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3467;
    (&_S3467)->primal_0 = R_25;
    (&_S3467)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3468;
    (&_S3468)->primal_0 = vert1_1;
    (&_S3468)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3467, &_S3468, _S3466);
    float3  _S3469 = _S3309 + _S3310 + _S3389.differential_0 + make_float3 (_S3461.x, _S3461.y, _S3462);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3470;
    (&_S3470)->primal_0 = R_25;
    (&_S3470)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3471;
    (&_S3471)->primal_0 = vert0_1;
    (&_S3471)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3470, &_S3471, _S3469);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3472;
    (&_S3472)->primal_0 = _S3188;
    (&_S3472)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3473;
    (&_S3473)->primal_0 = _S3193;
    (&_S3473)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3472, &_S3473, _S3465.differential_0);
    float _S3474 = - _S3473.differential_0.y;
    float _S3475 = _S3192 * _S3473.differential_0.x;
    float _S3476 = - (_S3184 * _S3473.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3477;
    (&_S3477)->primal_0 = _S3188;
    (&_S3477)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3478;
    (&_S3478)->primal_0 = _S3191;
    (&_S3478)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3477, &_S3478, _S3468.differential_0);
    float _S3479 = _S3184 * _S3478.differential_0.x;
    float _S3480 = _S3190 * _S3478.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3481;
    (&_S3481)->primal_0 = _S3188;
    (&_S3481)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3482;
    (&_S3482)->primal_0 = _S3189;
    (&_S3482)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3481, &_S3482, _S3471.differential_0);
    Matrix<float, 3, 3>  _S3483 = transpose_0(_S3472.differential_0 + _S3477.differential_0 + _S3481.differential_0);
    float _S3484 = 2.0f * - _S3483.rows[int(2)].z;
    float _S3485 = 2.0f * _S3483.rows[int(2)].y;
    float _S3486 = 2.0f * _S3483.rows[int(2)].x;
    float _S3487 = 2.0f * _S3483.rows[int(1)].z;
    float _S3488 = 2.0f * - _S3483.rows[int(1)].y;
    float _S3489 = 2.0f * _S3483.rows[int(1)].x;
    float _S3490 = 2.0f * _S3483.rows[int(0)].z;
    float _S3491 = 2.0f * _S3483.rows[int(0)].y;
    float _S3492 = 2.0f * - _S3483.rows[int(0)].x;
    float _S3493 = - _S3489 + _S3491;
    float _S3494 = _S3486 + - _S3490;
    float _S3495 = - _S3485 + _S3487;
    float _S3496 = _S3485 + _S3487;
    float _S3497 = _S3486 + _S3490;
    float _S3498 = _S3489 + _S3491;
    float _S3499 = quat_25.w * (_S3488 + _S3492);
    float _S3500 = quat_25.z * (_S3484 + _S3492);
    float _S3501 = quat_25.y * (_S3484 + _S3488);
    float _S3502 = quat_25.x * _S3493 + quat_25.z * _S3496 + quat_25.y * _S3497 + _S3499 + _S3499;
    float _S3503 = quat_25.x * _S3494 + quat_25.w * _S3496 + quat_25.y * _S3498 + _S3500 + _S3500;
    float _S3504 = quat_25.x * _S3495 + quat_25.w * _S3497 + quat_25.z * _S3498 + _S3501 + _S3501;
    float _S3505 = quat_25.w * _S3493 + quat_25.z * _S3494 + quat_25.y * _S3495;
    float _S3506 = _S3476 + _S3479;
    float _S3507 = 0.5f * - _S3506;
    float _S3508 = _S3474 + _S3478.differential_0.y;
    DiffPair_float_0 _S3509;
    (&_S3509)->primal_0 = _S3185;
    (&_S3509)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3509, _S3508);
    float _S3510 = _S3507 + _S3509.differential_0;
    float _S3511 = _S3475 + _S3480 + _S3482.differential_0.x;
    DiffPair_float_0 _S3512;
    (&_S3512)->primal_0 = _S3183;
    (&_S3512)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3512, _S3511);
    float _S3513 = _S3507 + _S3512.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3514;
    (&_S3514)->primal_0 = R_25;
    (&_S3514)->differential_0 = _S3382;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3515;
    (&_S3515)->primal_0 = mean_24;
    (&_S3515)->differential_0 = _S3302;
    s_bwd_prop_mul_1(&_S3514, &_S3515, _S3304.differential_0);
    float3  _S3516 = _S3384.differential_0 + _S3463 + _S3466 + _S3469 + _S3304.differential_0;
    Matrix<float, 3, 3>  _S3517 = _S3385 + _S3464.differential_0 + _S3467.differential_0 + _S3470.differential_0 + _S3514.differential_0;
    FixedArray<float3 , 2>  _S3518;
    _S3518[int(0)] = _S3302;
    _S3518[int(1)] = _S3302;
    _S3518[int(1)] = _S3316;
    _S3518[int(0)] = _S3321;
    FixedArray<float3 , 16>  _S3519;
    _S3519[int(0)] = _S3302;
    _S3519[int(1)] = _S3302;
    _S3519[int(2)] = _S3302;
    _S3519[int(3)] = _S3302;
    _S3519[int(4)] = _S3302;
    _S3519[int(5)] = _S3302;
    _S3519[int(6)] = _S3302;
    _S3519[int(7)] = _S3302;
    _S3519[int(8)] = _S3302;
    _S3519[int(9)] = _S3302;
    _S3519[int(10)] = _S3302;
    _S3519[int(11)] = _S3302;
    _S3519[int(12)] = _S3302;
    _S3519[int(13)] = _S3302;
    _S3519[int(14)] = _S3302;
    _S3519[int(15)] = _S3302;
    _S3519[int(15)] = _S3323;
    _S3519[int(14)] = _S3325;
    _S3519[int(13)] = _S3327;
    _S3519[int(12)] = _S3329;
    _S3519[int(11)] = _S3331;
    _S3519[int(10)] = _S3333;
    _S3519[int(9)] = _S3335;
    _S3519[int(8)] = _S3343;
    _S3519[int(7)] = _S3345;
    _S3519[int(6)] = _S3347;
    _S3519[int(5)] = _S3349;
    _S3519[int(4)] = _S3351;
    _S3519[int(3)] = _S3362;
    _S3519[int(2)] = _S3364;
    _S3519[int(1)] = _S3366;
    _S3519[int(0)] = _S3379;
    float2  _S3520 = v_out_hardness_0 + make_float2 (0.0f, _S3407);
    float3  _S3521 = make_float3 (_S3513, _S3510, _S3506);
    float4  _S3522 = make_float4 (0.0f);
    *&((&_S3522)->w) = _S3502;
    *&((&_S3522)->z) = _S3503;
    *&((&_S3522)->y) = _S3504;
    *&((&_S3522)->x) = _S3505;
    *v_mean_7 = _S3380 + _S3465.differential_0 + _S3468.differential_0 + _S3471.differential_0 + _S3515.differential_0;
    *v_quat_6 = _S3522;
    *v_scale_6 = _S3521;
    *v_hardness_0 = _S3520;
    *v_sh_coeffs_5 = _S3519;
    *v_ch_coeffs_0 = _S3518;
    *v_R_7 = _S3517;
    *v_t_6 = _S3516;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_26, float3  t_23, float fx_29, float fy_29, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_0(R_26, mean_25) + t_23;
    float _S3523 = scale_25.x;
    float _S3524 = s_primal_ctx_exp_1(_S3523);
    float _S3525 = scale_25.y;
    float _S3526 = s_primal_ctx_exp_1(_S3525);
    float sz_6 = scale_25.z - 0.5f * (_S3523 + _S3525);
    float _S3527 = quat_26.y;
    float x2_26 = _S3527 * _S3527;
    float y2_26 = quat_26.z * quat_26.z;
    float z2_46 = quat_26.w * quat_26.w;
    float xy_26 = quat_26.y * quat_26.z;
    float xz_26 = quat_26.y * quat_26.w;
    float yz_26 = quat_26.z * quat_26.w;
    float wx_26 = quat_26.x * quat_26.y;
    float wy_26 = quat_26.x * quat_26.z;
    float wz_26 = quat_26.x * quat_26.w;
    Matrix<float, 3, 3>  _S3528 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_46), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_46), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
    float3  _S3529 = make_float3 (_S3524, 0.0f, 0.0f);
    float3  vert0_2 = s_primal_ctx_mul_0(_S3528, _S3529) + mean_25;
    float _S3530 = -0.5f + sz_6;
    float3  _S3531 = make_float3 (_S3524 * _S3530, _S3526, 0.0f);
    float3  vert1_2 = s_primal_ctx_mul_0(_S3528, _S3531) + mean_25;
    float _S3532 = -0.5f - sz_6;
    float3  _S3533 = make_float3 (_S3524 * _S3532, - _S3526, 0.0f);
    float3  vert2_2 = s_primal_ctx_mul_0(_S3528, _S3533) + mean_25;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_26, vert0_2) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_26, vert1_2) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_26, vert2_2) + t_23;
    float _S3534 = length_1(vert0_c_5);
    float _S3535 = length_1(vert1_c_5);
    float _S3536 = length_1(vert2_c_5);
    float2  _S3537 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3538 = length_0(_S3537);
    float _S3539 = vert0_c_5.z;
    float _S3540 = s_primal_ctx_atan2_0(_S3538, _S3539);
    bool _S3541 = _S3540 < 0.00100000004749745f;
    float k_9;
    float _S3542;
    float _S3543;
    float _S3544;
    if(_S3541)
    {
        float _S3545 = 1.0f - _S3540 * _S3540 / 3.0f;
        float _S3546 = _S3539 * _S3539;
        k_9 = _S3545 / _S3539;
        _S3542 = 0.0f;
        _S3543 = _S3546;
        _S3544 = _S3545;
    }
    else
    {
        float _S3547 = _S3538 * _S3538;
        k_9 = _S3540 / _S3538;
        _S3542 = _S3547;
        _S3543 = 0.0f;
        _S3544 = 0.0f;
    }
    float2  _S3548 = make_float2 (k_9);
    float2  _S3549 = _S3537 * make_float2 (k_9);
    float u_59 = _S3549.x;
    float v_59 = _S3549.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S3550 = (*dist_coeffs_33)[int(2)] + r2_59 * (*dist_coeffs_33)[int(3)];
    float _S3551 = (*dist_coeffs_33)[int(1)] + r2_59 * _S3550;
    float _S3552 = (*dist_coeffs_33)[int(0)] + r2_59 * _S3551;
    float radial_4 = 1.0f + r2_59 * _S3552;
    float2  _S3553 = make_float2 (radial_4);
    float _S3554 = 2.0f * (*dist_coeffs_33)[int(4)];
    float _S3555 = _S3554 * u_59;
    float _S3556 = 2.0f * u_59;
    float _S3557 = 2.0f * (*dist_coeffs_33)[int(5)];
    float _S3558 = _S3557 * u_59;
    float _S3559 = 2.0f * v_59;
    float2  _S3560 = _S3549 * make_float2 (radial_4) + make_float2 (_S3555 * v_59 + (*dist_coeffs_33)[int(5)] * (r2_59 + _S3556 * u_59) + (*dist_coeffs_33)[int(6)] * r2_59, _S3558 * v_59 + (*dist_coeffs_33)[int(4)] * (r2_59 + _S3559 * v_59) + (*dist_coeffs_33)[int(7)] * r2_59);
    float2  _S3561 = _S3560 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3560.x + (*dist_coeffs_33)[int(9)] * _S3560.y, 0.0f);
    float _S3562 = fx_29 * _S3561.x + cx_24;
    float _S3563 = fy_29 * _S3561.y + cy_24;
    float2  _S3564 = make_float2 (_S3562, _S3563);
    float2  _S3565 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3566 = length_0(_S3565);
    float _S3567 = vert1_c_5.z;
    float _S3568 = s_primal_ctx_atan2_0(_S3566, _S3567);
    bool _S3569 = _S3568 < 0.00100000004749745f;
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
    float u_60 = _S3577.x;
    float v_60 = _S3577.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S3578 = (*dist_coeffs_33)[int(2)] + r2_60 * (*dist_coeffs_33)[int(3)];
    float _S3579 = (*dist_coeffs_33)[int(1)] + r2_60 * _S3578;
    float _S3580 = (*dist_coeffs_33)[int(0)] + r2_60 * _S3579;
    float radial_5 = 1.0f + r2_60 * _S3580;
    float2  _S3581 = make_float2 (radial_5);
    float _S3582 = _S3554 * u_60;
    float _S3583 = 2.0f * u_60;
    float _S3584 = _S3557 * u_60;
    float _S3585 = 2.0f * v_60;
    float2  _S3586 = _S3577 * make_float2 (radial_5) + make_float2 (_S3582 * v_60 + (*dist_coeffs_33)[int(5)] * (r2_60 + _S3583 * u_60) + (*dist_coeffs_33)[int(6)] * r2_60, _S3584 * v_60 + (*dist_coeffs_33)[int(4)] * (r2_60 + _S3585 * v_60) + (*dist_coeffs_33)[int(7)] * r2_60);
    float2  _S3587 = _S3586 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3586.x + (*dist_coeffs_33)[int(9)] * _S3586.y, 0.0f);
    float _S3588 = fx_29 * _S3587.x + cx_24;
    float _S3589 = fy_29 * _S3587.y + cy_24;
    float2  _S3590 = make_float2 (_S3588, _S3589);
    float2  _S3591 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3592 = length_0(_S3591);
    float _S3593 = vert2_c_5.z;
    float _S3594 = s_primal_ctx_atan2_0(_S3592, _S3593);
    bool _S3595 = _S3594 < 0.00100000004749745f;
    float _S3596;
    float _S3597;
    float _S3598;
    if(_S3595)
    {
        float _S3599 = 1.0f - _S3594 * _S3594 / 3.0f;
        float _S3600 = _S3593 * _S3593;
        k_9 = _S3599 / _S3593;
        _S3596 = 0.0f;
        _S3597 = _S3600;
        _S3598 = _S3599;
    }
    else
    {
        float _S3601 = _S3592 * _S3592;
        k_9 = _S3594 / _S3592;
        _S3596 = _S3601;
        _S3597 = 0.0f;
        _S3598 = 0.0f;
    }
    float2  _S3602 = make_float2 (k_9);
    float2  _S3603 = _S3591 * make_float2 (k_9);
    float u_61 = _S3603.x;
    float v_61 = _S3603.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S3604 = (*dist_coeffs_33)[int(2)] + r2_61 * (*dist_coeffs_33)[int(3)];
    float _S3605 = (*dist_coeffs_33)[int(1)] + r2_61 * _S3604;
    float _S3606 = (*dist_coeffs_33)[int(0)] + r2_61 * _S3605;
    float radial_6 = 1.0f + r2_61 * _S3606;
    float2  _S3607 = make_float2 (radial_6);
    float _S3608 = _S3554 * u_61;
    float _S3609 = 2.0f * u_61;
    float _S3610 = _S3557 * u_61;
    float _S3611 = 2.0f * v_61;
    float2  _S3612 = _S3603 * make_float2 (radial_6) + make_float2 (_S3608 * v_61 + (*dist_coeffs_33)[int(5)] * (r2_61 + _S3609 * u_61) + (*dist_coeffs_33)[int(6)] * r2_61, _S3610 * v_61 + (*dist_coeffs_33)[int(4)] * (r2_61 + _S3611 * v_61) + (*dist_coeffs_33)[int(7)] * r2_61);
    float2  _S3613 = _S3612 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3612.x + (*dist_coeffs_33)[int(9)] * _S3612.y, 0.0f);
    float _S3614 = fx_29 * _S3613.x + cx_24;
    float _S3615 = fy_29 * _S3613.y + cy_24;
    float2  _S3616 = make_float2 (_S3614, _S3615);
    float2  e0_5 = _S3590 - _S3564;
    float2  e1_5 = _S3616 - _S3590;
    float2  e2_1 = _S3564 - _S3616;
    float _S3617 = e0_5.x;
    float _S3618 = e1_5.y;
    float _S3619 = e0_5.y;
    float _S3620 = e1_5.x;
    float _S3621 = _S3617 * _S3618 - _S3619 * _S3620;
    float _S3622 = 1.0f - hardness_5.y;
    float _S3623 = -1.0f / _S3622;
    float _S3624 = _S3622 * _S3622;
    float _S3625 = s_primal_ctx_max_0(_S3562, _S3588);
    float _S3626 = s_primal_ctx_min_0(_S3562, _S3588);
    float _S3627 = s_primal_ctx_max_0(_S3563, _S3589);
    float _S3628 = s_primal_ctx_min_0(_S3563, _S3589);
    float3  _S3629 = make_float3 (_S3534, _S3535, _S3536) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3630 = transpose_0(R_26);
    float3  _S3631 = mean_25 - - s_primal_ctx_mul_0(_S3630, t_23);
    float _S3632 = _S3631.x;
    float _S3633 = _S3631.y;
    float _S3634 = _S3631.z;
    float _S3635 = _S3632 * _S3632 + _S3633 * _S3633 + _S3634 * _S3634;
    float _S3636 = s_primal_ctx_sqrt_0(_S3635);
    float x_52 = _S3632 / _S3636;
    float3  _S3637 = make_float3 (x_52);
    float _S3638 = _S3636 * _S3636;
    float y_23 = _S3633 / _S3636;
    float z_20 = _S3634 / _S3636;
    float3  _S3639 = make_float3 (z_20);
    float _S3640 = - y_23;
    float3  _S3641 = make_float3 (_S3640);
    float z2_47 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_52 * x_52 - y_23 * y_23;
    float _S3642 = 2.0f * x_52;
    float fS1_20 = _S3642 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_47 - 0.31539157032966614f;
    float3  _S3643 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_52;
    float3  _S3644 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3645 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3646 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3647 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_47 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3648 = 1.86588168144226074f * z2_47 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3648;
    float3  _S3649 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_52;
    float3  _S3650 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3651 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3652 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3653 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_52 * fC1_20 - y_23 * fS1_20);
    float3  _S3654 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_52 * fS1_20 + y_23 * fC1_20);
    float3  _S3655 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3640) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_52) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3656 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3657 = make_float3 (0.0f);
    float3  _S3658 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3659 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3660 = make_float3 (_S3659);
    float3  _S3661 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3659);
    float3  _S3662 = _S3658 + _S3661 + make_float3 (0.5f);
    float3  _S3663 = _S3658 - _S3661 + make_float3 (0.5f);
    float3  _S3664 = vert1_c_5 - vert0_c_5;
    float3  _S3665 = vert2_c_5 - vert0_c_5;
    float3  _S3666 = s_primal_ctx_cross_0(_S3664, _S3665);
    float3  _S3667 = normalize_0(_S3666);
    float3  _S3668 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3667, mean_c_20)))))) * v_normal_1;
    float3  _S3669 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3670;
    (&_S3670)->primal_0 = _S3667;
    (&_S3670)->differential_0 = _S3669;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3671;
    (&_S3671)->primal_0 = mean_c_20;
    (&_S3671)->differential_0 = _S3669;
    s_bwd_prop_dot_0(&_S3670, &_S3671, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3672 = _S3671;
    float3  _S3673 = _S3668 + _S3670.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3674;
    (&_S3674)->primal_0 = _S3666;
    (&_S3674)->differential_0 = _S3669;
    s_bwd_normalize_impl_0(&_S3674, _S3673);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3675;
    (&_S3675)->primal_0 = _S3664;
    (&_S3675)->differential_0 = _S3669;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3676;
    (&_S3676)->primal_0 = _S3665;
    (&_S3676)->differential_0 = _S3669;
    s_bwd_prop_cross_0(&_S3675, &_S3676, _S3674.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3677 = _S3675;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3678 = _S3676;
    float3  _S3679 = - _S3676.differential_0;
    float3  _S3680 = - _S3675.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3681;
    (&_S3681)->primal_0 = _S3663;
    (&_S3681)->differential_0 = _S3669;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3682;
    (&_S3682)->primal_0 = _S3657;
    (&_S3682)->differential_0 = _S3669;
    s_bwd_prop_max_0(&_S3681, &_S3682, (*v_rgb_7)[int(2)]);
    float3  _S3683 = - _S3681.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3684;
    (&_S3684)->primal_0 = _S3662;
    (&_S3684)->differential_0 = _S3669;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3685;
    (&_S3685)->primal_0 = _S3657;
    (&_S3685)->differential_0 = _S3669;
    s_bwd_prop_max_0(&_S3684, &_S3685, (*v_rgb_7)[int(1)]);
    float3  _S3686 = _S3660 * (_S3683 + _S3684.differential_0);
    float3  _S3687 = _S3681.differential_0 + _S3684.differential_0;
    float3  _S3688 = make_float3 (0.5f) * - _S3687;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3689;
    (&_S3689)->primal_0 = _S3656;
    (&_S3689)->differential_0 = _S3669;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3690;
    (&_S3690)->primal_0 = _S3657;
    (&_S3690)->differential_0 = _S3669;
    s_bwd_prop_max_0(&_S3689, &_S3690, (*v_rgb_7)[int(0)]);
    float3  _S3691 = _S3688 + _S3689.differential_0;
    float3  _S3692 = _S3687 + _S3689.differential_0;
    float3  _S3693 = _S3654 * _S3692;
    float3  _S3694 = (*sh_coeffs_20)[int(15)] * _S3692;
    float3  _S3695 = _S3652 * _S3692;
    float3  _S3696 = (*sh_coeffs_20)[int(14)] * _S3692;
    float3  _S3697 = _S3650 * _S3692;
    float3  _S3698 = (*sh_coeffs_20)[int(13)] * _S3692;
    float3  _S3699 = _S3649 * _S3692;
    float3  _S3700 = (*sh_coeffs_20)[int(12)] * _S3692;
    float3  _S3701 = _S3651 * _S3692;
    float3  _S3702 = (*sh_coeffs_20)[int(11)] * _S3692;
    float3  _S3703 = _S3653 * _S3692;
    float3  _S3704 = (*sh_coeffs_20)[int(10)] * _S3692;
    float3  _S3705 = _S3655 * _S3692;
    float3  _S3706 = (*sh_coeffs_20)[int(9)] * _S3692;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3706.x + _S3706.y + _S3706.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3694.x + _S3694.y + _S3694.z);
    float _S3707 = _S3704.x + _S3704.y + _S3704.z;
    float _S3708 = _S3696.x + _S3696.y + _S3696.z;
    float _S3709 = _S3702.x + _S3702.y + _S3702.z;
    float _S3710 = _S3698.x + _S3698.y + _S3698.z;
    float _S3711 = _S3700.x + _S3700.y + _S3700.z;
    float _S3712 = - s_diff_fC2_T_6;
    float3  _S3713 = _S3646 * _S3692;
    float3  _S3714 = (*sh_coeffs_20)[int(8)] * _S3692;
    float3  _S3715 = _S3644 * _S3692;
    float3  _S3716 = (*sh_coeffs_20)[int(7)] * _S3692;
    float3  _S3717 = _S3643 * _S3692;
    float3  _S3718 = (*sh_coeffs_20)[int(6)] * _S3692;
    float3  _S3719 = _S3645 * _S3692;
    float3  _S3720 = (*sh_coeffs_20)[int(5)] * _S3692;
    float3  _S3721 = _S3647 * _S3692;
    float3  _S3722 = (*sh_coeffs_20)[int(4)] * _S3692;
    float _S3723 = _S3720.x + _S3720.y + _S3720.z;
    float _S3724 = _S3716.x + _S3716.y + _S3716.z;
    float _S3725 = fTmp1B_20 * _S3707 + x_52 * s_diff_fS2_T_6 + y_23 * _S3712 + 0.54627424478530884f * (_S3722.x + _S3722.y + _S3722.z);
    float _S3726 = fTmp1B_20 * _S3708 + y_23 * s_diff_fS2_T_6 + x_52 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3714.x + _S3714.y + _S3714.z);
    float _S3727 = y_23 * - _S3726;
    float _S3728 = x_52 * _S3726;
    float _S3729 = z_20 * (1.86588168144226074f * (z_20 * _S3711) + -2.28522896766662598f * (y_23 * _S3709 + x_52 * _S3710) + 0.94617468118667603f * (_S3718.x + _S3718.y + _S3718.z));
    float3  _S3730 = make_float3 (0.48860251903533936f) * _S3692;
    float3  _S3731 = - _S3730;
    float3  _S3732 = _S3637 * _S3731;
    float3  _S3733 = (*sh_coeffs_20)[int(3)] * _S3731;
    float3  _S3734 = _S3639 * _S3730;
    float3  _S3735 = (*sh_coeffs_20)[int(2)] * _S3730;
    float3  _S3736 = _S3641 * _S3730;
    float3  _S3737 = (*sh_coeffs_20)[int(1)] * _S3730;
    float _S3738 = (_S3648 * _S3711 + 1.44530570507049561f * (fS1_20 * _S3707 + fC1_20 * _S3708) + -1.09254848957061768f * (y_23 * _S3723 + x_52 * _S3724) + _S3729 + _S3729 + _S3735.x + _S3735.y + _S3735.z) / _S3638;
    float _S3739 = _S3636 * _S3738;
    float _S3740 = (fTmp0C_20 * _S3709 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3712 + fTmp0B_20 * _S3723 + _S3642 * _S3725 + _S3727 + _S3727 + - (_S3737.x + _S3737.y + _S3737.z)) / _S3638;
    float _S3741 = _S3636 * _S3740;
    float _S3742 = (fTmp0C_20 * _S3710 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3724 + 2.0f * (y_23 * _S3725) + _S3728 + _S3728 + _S3733.x + _S3733.y + _S3733.z) / _S3638;
    float _S3743 = _S3636 * _S3742;
    float _S3744 = _S3634 * - _S3738 + _S3633 * - _S3740 + _S3632 * - _S3742;
    DiffPair_float_0 _S3745;
    (&_S3745)->primal_0 = _S3635;
    (&_S3745)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3745, _S3744);
    float _S3746 = _S3634 * _S3745.differential_0;
    float _S3747 = _S3633 * _S3745.differential_0;
    float _S3748 = _S3632 * _S3745.differential_0;
    float3  _S3749 = make_float3 (0.282094806432724f) * _S3692;
    float3  _S3750 = make_float3 (_S3743 + _S3748 + _S3748, _S3741 + _S3747 + _S3747, _S3739 + _S3746 + _S3746);
    float3  _S3751 = - - _S3750;
    Matrix<float, 3, 3>  _S3752 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3753;
    (&_S3753)->primal_0 = _S3630;
    (&_S3753)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3754;
    (&_S3754)->primal_0 = t_23;
    (&_S3754)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3753, &_S3754, _S3751);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3755 = _S3754;
    Matrix<float, 3, 3>  _S3756 = transpose_0(_S3753.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3757;
    (&_S3757)->primal_0 = _S3629;
    (&_S3757)->differential_0 = _S3669;
    s_bwd_prop_log_1(&_S3757, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3758 = _S3757;
    DiffPair_float_0 _S3759;
    (&_S3759)->primal_0 = _S3628;
    (&_S3759)->differential_0 = 0.0f;
    DiffPair_float_0 _S3760;
    (&_S3760)->primal_0 = _S3615;
    (&_S3760)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3759, &_S3760, 0.0f);
    DiffPair_float_0 _S3761;
    (&_S3761)->primal_0 = _S3563;
    (&_S3761)->differential_0 = 0.0f;
    DiffPair_float_0 _S3762;
    (&_S3762)->primal_0 = _S3589;
    (&_S3762)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3761, &_S3762, _S3759.differential_0);
    DiffPair_float_0 _S3763;
    (&_S3763)->primal_0 = _S3627;
    (&_S3763)->differential_0 = 0.0f;
    DiffPair_float_0 _S3764;
    (&_S3764)->primal_0 = _S3615;
    (&_S3764)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3763, &_S3764, 0.0f);
    DiffPair_float_0 _S3765;
    (&_S3765)->primal_0 = _S3563;
    (&_S3765)->differential_0 = 0.0f;
    DiffPair_float_0 _S3766;
    (&_S3766)->primal_0 = _S3589;
    (&_S3766)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3765, &_S3766, _S3763.differential_0);
    DiffPair_float_0 _S3767;
    (&_S3767)->primal_0 = _S3626;
    (&_S3767)->differential_0 = 0.0f;
    DiffPair_float_0 _S3768;
    (&_S3768)->primal_0 = _S3614;
    (&_S3768)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3767, &_S3768, 0.0f);
    DiffPair_float_0 _S3769;
    (&_S3769)->primal_0 = _S3562;
    (&_S3769)->differential_0 = 0.0f;
    DiffPair_float_0 _S3770;
    (&_S3770)->primal_0 = _S3588;
    (&_S3770)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3769, &_S3770, _S3767.differential_0);
    DiffPair_float_0 _S3771;
    (&_S3771)->primal_0 = _S3625;
    (&_S3771)->differential_0 = 0.0f;
    DiffPair_float_0 _S3772;
    (&_S3772)->primal_0 = _S3614;
    (&_S3772)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3771, &_S3772, 0.0f);
    DiffPair_float_0 _S3773;
    (&_S3773)->primal_0 = _S3562;
    (&_S3773)->differential_0 = 0.0f;
    DiffPair_float_0 _S3774;
    (&_S3774)->primal_0 = _S3588;
    (&_S3774)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3773, &_S3774, _S3771.differential_0);
    DiffPair_float_0 _S3775;
    (&_S3775)->primal_0 = _S3623;
    (&_S3775)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3775, 0.0f);
    float _S3776 = - (-1.0f * - (_S3775.differential_0 / _S3624));
    float2  _S3777 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3778;
    (&_S3778)->primal_0 = e2_1;
    (&_S3778)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3778, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3779;
    (&_S3779)->primal_0 = e1_5;
    (&_S3779)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3779, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3780;
    (&_S3780)->primal_0 = e0_5;
    (&_S3780)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3780, -0.0f);
    DiffPair_float_0 _S3781;
    (&_S3781)->primal_0 = _S3621;
    (&_S3781)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3781, 0.0f);
    float _S3782 = - _S3781.differential_0;
    float2  _S3783 = _S3779.differential_0 + make_float2 (_S3619 * _S3782, _S3617 * _S3781.differential_0);
    float2  _S3784 = _S3780.differential_0 + make_float2 (_S3618 * _S3781.differential_0, _S3620 * _S3782);
    float2  _S3785 = v_uv2_1 + - _S3778.differential_0 + _S3783;
    float _S3786 = fx_29 * (_S3768.differential_0 + _S3772.differential_0 + _S3785.x);
    float2  _S3787 = make_float2 (_S3786, fy_29 * (_S3760.differential_0 + _S3764.differential_0 + _S3785.y)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3786, (*dist_coeffs_33)[int(9)] * _S3786);
    float2  _S3788 = _S3603 * _S3787;
    float2  _S3789 = _S3607 * _S3787;
    float _S3790 = (*dist_coeffs_33)[int(4)] * _S3787.y;
    float _S3791 = (*dist_coeffs_33)[int(5)] * _S3787.x;
    float _S3792 = _S3788.x + _S3788.y;
    float _S3793 = r2_61 * _S3792;
    float _S3794 = r2_61 * _S3793;
    float _S3795 = (*dist_coeffs_33)[int(7)] * _S3787.y + _S3790 + (*dist_coeffs_33)[int(6)] * _S3787.x + _S3791 + _S3606 * _S3792 + _S3605 * _S3793 + _S3604 * _S3794 + (*dist_coeffs_33)[int(3)] * (r2_61 * _S3794);
    float _S3796 = v_61 * _S3795;
    float _S3797 = u_61 * _S3795;
    float _S3798 = _S3611 * _S3790 + 2.0f * (v_61 * _S3790) + _S3610 * _S3787.y + _S3608 * _S3787.x + _S3796 + _S3796;
    float _S3799 = _S3557 * (v_61 * _S3787.y) + _S3609 * _S3791 + 2.0f * (u_61 * _S3791) + _S3554 * (v_61 * _S3787.x) + _S3797 + _S3797;
    float2  _S3800 = v_uv0_1 + _S3778.differential_0 + - _S3784;
    float2  _S3801 = v_out_hardness_1 + make_float2 (0.0f, _S3776);
    float _S3802 = _S3770.differential_0 + _S3774.differential_0;
    float2  _S3803 = v_uv1_1 + - _S3783 + _S3784;
    float3  _S3804 = _S3679 + _S3680;
    FixedArray<float3 , 2>  _S3805;
    _S3805[int(0)] = _S3669;
    _S3805[int(1)] = _S3669;
    _S3805[int(1)] = _S3686;
    _S3805[int(0)] = _S3691;
    float3  _S3806 = _S3805[int(0)];
    float3  _S3807 = _S3805[int(1)];
    FixedArray<float3 , 16>  _S3808;
    _S3808[int(0)] = _S3669;
    _S3808[int(1)] = _S3669;
    _S3808[int(2)] = _S3669;
    _S3808[int(3)] = _S3669;
    _S3808[int(4)] = _S3669;
    _S3808[int(5)] = _S3669;
    _S3808[int(6)] = _S3669;
    _S3808[int(7)] = _S3669;
    _S3808[int(8)] = _S3669;
    _S3808[int(9)] = _S3669;
    _S3808[int(10)] = _S3669;
    _S3808[int(11)] = _S3669;
    _S3808[int(12)] = _S3669;
    _S3808[int(13)] = _S3669;
    _S3808[int(14)] = _S3669;
    _S3808[int(15)] = _S3669;
    _S3808[int(7)] = _S3715;
    _S3808[int(0)] = _S3749;
    _S3808[int(1)] = _S3736;
    _S3808[int(2)] = _S3734;
    _S3808[int(3)] = _S3732;
    _S3808[int(4)] = _S3721;
    _S3808[int(5)] = _S3719;
    _S3808[int(6)] = _S3717;
    _S3808[int(15)] = _S3693;
    _S3808[int(8)] = _S3713;
    _S3808[int(9)] = _S3705;
    _S3808[int(10)] = _S3703;
    _S3808[int(11)] = _S3701;
    _S3808[int(12)] = _S3699;
    _S3808[int(13)] = _S3697;
    _S3808[int(14)] = _S3695;
    float3  _S3809 = _S3808[int(0)];
    float3  _S3810 = _S3808[int(1)];
    float3  _S3811 = _S3808[int(2)];
    float3  _S3812 = _S3808[int(3)];
    float3  _S3813 = _S3808[int(4)];
    float3  _S3814 = _S3808[int(5)];
    float3  _S3815 = _S3808[int(6)];
    float3  _S3816 = _S3808[int(7)];
    float3  _S3817 = _S3808[int(8)];
    float3  _S3818 = _S3808[int(9)];
    float3  _S3819 = _S3808[int(10)];
    float3  _S3820 = _S3808[int(11)];
    float3  _S3821 = _S3808[int(12)];
    float3  _S3822 = _S3808[int(13)];
    float3  _S3823 = _S3808[int(14)];
    float3  _S3824 = _S3808[int(15)];
    float _S3825 = _S3769.differential_0 + _S3773.differential_0;
    float _S3826 = _S3761.differential_0 + _S3765.differential_0;
    float _S3827 = _S3762.differential_0 + _S3766.differential_0;
    float2  _S3828 = _S3789 + make_float2 (_S3799, _S3798);
    float2  _S3829 = _S3591 * _S3828;
    float2  _S3830 = _S3602 * _S3828;
    float _S3831 = _S3829.x + _S3829.y;
    if(_S3595)
    {
        float _S3832 = _S3831 / _S3597;
        float _S3833 = _S3598 * - _S3832;
        float _S3834 = _S3594 * (0.3333333432674408f * - (_S3593 * _S3832));
        k_9 = _S3834 + _S3834;
        _S3596 = _S3833;
        _S3597 = 0.0f;
    }
    else
    {
        float _S3835 = _S3831 / _S3596;
        float _S3836 = _S3594 * - _S3835;
        k_9 = _S3592 * _S3835;
        _S3596 = 0.0f;
        _S3597 = _S3836;
    }
    DiffPair_float_0 _S3837;
    (&_S3837)->primal_0 = _S3592;
    (&_S3837)->differential_0 = 0.0f;
    DiffPair_float_0 _S3838;
    (&_S3838)->primal_0 = _S3593;
    (&_S3838)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3837, &_S3838, k_9);
    float _S3839 = _S3838.differential_0 + _S3596;
    float _S3840 = _S3837.differential_0 + _S3597;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3841;
    (&_S3841)->primal_0 = _S3591;
    (&_S3841)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3841, _S3840);
    float2  _S3842 = _S3841.differential_0 + _S3830;
    float _S3843 = fx_29 * (_S3803.x + _S3802);
    float2  _S3844 = make_float2 (_S3843, fy_29 * (_S3803.y + _S3827)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3843, (*dist_coeffs_33)[int(9)] * _S3843);
    float2  _S3845 = _S3577 * _S3844;
    float _S3846 = (*dist_coeffs_33)[int(4)] * _S3844.y;
    float _S3847 = (*dist_coeffs_33)[int(5)] * _S3844.x;
    float _S3848 = _S3845.x + _S3845.y;
    float _S3849 = r2_60 * _S3848;
    float _S3850 = r2_60 * _S3849;
    float _S3851 = (*dist_coeffs_33)[int(7)] * _S3844.y + _S3846 + (*dist_coeffs_33)[int(6)] * _S3844.x + _S3847 + _S3580 * _S3848 + _S3579 * _S3849 + _S3578 * _S3850 + (*dist_coeffs_33)[int(3)] * (r2_60 * _S3850);
    float _S3852 = v_60 * _S3851;
    float _S3853 = u_60 * _S3851;
    float3  _S3854 = _S3678.differential_0 + make_float3 (_S3842.x, _S3842.y, _S3839);
    float2  _S3855 = _S3581 * _S3844 + make_float2 (_S3557 * (v_60 * _S3844.y) + _S3583 * _S3847 + 2.0f * (u_60 * _S3847) + _S3554 * (v_60 * _S3844.x) + _S3853 + _S3853, _S3585 * _S3846 + 2.0f * (v_60 * _S3846) + _S3584 * _S3844.y + _S3582 * _S3844.x + _S3852 + _S3852);
    float2  _S3856 = _S3565 * _S3855;
    float2  _S3857 = _S3576 * _S3855;
    float _S3858 = _S3856.x + _S3856.y;
    if(_S3569)
    {
        float _S3859 = _S3858 / _S3571;
        float _S3860 = _S3572 * - _S3859;
        float _S3861 = _S3568 * (0.3333333432674408f * - (_S3567 * _S3859));
        k_9 = _S3861 + _S3861;
        _S3570 = _S3860;
        _S3571 = 0.0f;
    }
    else
    {
        float _S3862 = _S3858 / _S3570;
        float _S3863 = _S3568 * - _S3862;
        k_9 = _S3566 * _S3862;
        _S3570 = 0.0f;
        _S3571 = _S3863;
    }
    DiffPair_float_0 _S3864;
    (&_S3864)->primal_0 = _S3566;
    (&_S3864)->differential_0 = 0.0f;
    DiffPair_float_0 _S3865;
    (&_S3865)->primal_0 = _S3567;
    (&_S3865)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3864, &_S3865, k_9);
    float _S3866 = _S3865.differential_0 + _S3570;
    float _S3867 = _S3864.differential_0 + _S3571;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3868;
    (&_S3868)->primal_0 = _S3565;
    (&_S3868)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3868, _S3867);
    float2  _S3869 = _S3868.differential_0 + _S3857;
    float _S3870 = fx_29 * (_S3800.x + _S3825);
    float2  _S3871 = make_float2 (_S3870, fy_29 * (_S3800.y + _S3826)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3870, (*dist_coeffs_33)[int(9)] * _S3870);
    float2  _S3872 = _S3549 * _S3871;
    float _S3873 = (*dist_coeffs_33)[int(4)] * _S3871.y;
    float _S3874 = (*dist_coeffs_33)[int(5)] * _S3871.x;
    float _S3875 = _S3872.x + _S3872.y;
    float _S3876 = r2_59 * _S3875;
    float _S3877 = r2_59 * _S3876;
    float _S3878 = (*dist_coeffs_33)[int(7)] * _S3871.y + _S3873 + (*dist_coeffs_33)[int(6)] * _S3871.x + _S3874 + _S3552 * _S3875 + _S3551 * _S3876 + _S3550 * _S3877 + (*dist_coeffs_33)[int(3)] * (r2_59 * _S3877);
    float _S3879 = v_59 * _S3878;
    float _S3880 = u_59 * _S3878;
    float3  _S3881 = _S3677.differential_0 + make_float3 (_S3869.x, _S3869.y, _S3866);
    float2  _S3882 = _S3553 * _S3871 + make_float2 (_S3557 * (v_59 * _S3871.y) + _S3556 * _S3874 + 2.0f * (u_59 * _S3874) + _S3554 * (v_59 * _S3871.x) + _S3880 + _S3880, _S3559 * _S3873 + 2.0f * (v_59 * _S3873) + _S3558 * _S3871.y + _S3555 * _S3871.x + _S3879 + _S3879);
    float2  _S3883 = _S3537 * _S3882;
    float2  _S3884 = _S3548 * _S3882;
    float _S3885 = _S3883.x + _S3883.y;
    if(_S3541)
    {
        float _S3886 = _S3885 / _S3543;
        float _S3887 = _S3544 * - _S3886;
        float _S3888 = _S3540 * (0.3333333432674408f * - (_S3539 * _S3886));
        k_9 = _S3888 + _S3888;
        _S3542 = _S3887;
        _S3543 = 0.0f;
    }
    else
    {
        float _S3889 = _S3885 / _S3542;
        float _S3890 = _S3540 * - _S3889;
        k_9 = _S3538 * _S3889;
        _S3542 = 0.0f;
        _S3543 = _S3890;
    }
    DiffPair_float_0 _S3891;
    (&_S3891)->primal_0 = _S3538;
    (&_S3891)->differential_0 = 0.0f;
    DiffPair_float_0 _S3892;
    (&_S3892)->primal_0 = _S3539;
    (&_S3892)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3891, &_S3892, k_9);
    float _S3893 = _S3892.differential_0 + _S3542;
    float _S3894 = _S3891.differential_0 + _S3543;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3895;
    (&_S3895)->primal_0 = _S3537;
    (&_S3895)->differential_0 = _S3777;
    s_bwd_length_impl_0(&_S3895, _S3894);
    float2  _S3896 = _S3895.differential_0 + _S3884;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3897;
    (&_S3897)->primal_0 = vert2_c_5;
    (&_S3897)->differential_0 = _S3669;
    s_bwd_length_impl_1(&_S3897, _S3758.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3898;
    (&_S3898)->primal_0 = vert1_c_5;
    (&_S3898)->differential_0 = _S3669;
    s_bwd_length_impl_1(&_S3898, _S3758.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3899;
    (&_S3899)->primal_0 = vert0_c_5;
    (&_S3899)->differential_0 = _S3669;
    s_bwd_length_impl_1(&_S3899, _S3758.differential_0.x);
    float3  _S3900 = _S3897.differential_0 + _S3854;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3901;
    (&_S3901)->primal_0 = R_26;
    (&_S3901)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3902;
    (&_S3902)->primal_0 = vert2_2;
    (&_S3902)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3901, &_S3902, _S3900);
    float3  _S3903 = _S3898.differential_0 + _S3881;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3904;
    (&_S3904)->primal_0 = R_26;
    (&_S3904)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3905;
    (&_S3905)->primal_0 = vert1_2;
    (&_S3905)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3904, &_S3905, _S3903);
    float3  _S3906 = _S3899.differential_0 + _S3804 + make_float3 (_S3896.x, _S3896.y, _S3893);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3907;
    (&_S3907)->primal_0 = R_26;
    (&_S3907)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3908;
    (&_S3908)->primal_0 = vert0_2;
    (&_S3908)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3907, &_S3908, _S3906);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3909;
    (&_S3909)->primal_0 = _S3528;
    (&_S3909)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3910;
    (&_S3910)->primal_0 = _S3533;
    (&_S3910)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3909, &_S3910, _S3902.differential_0);
    float _S3911 = - _S3910.differential_0.y;
    float _S3912 = _S3532 * _S3910.differential_0.x;
    float _S3913 = - (_S3524 * _S3910.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3914;
    (&_S3914)->primal_0 = _S3528;
    (&_S3914)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3915;
    (&_S3915)->primal_0 = _S3531;
    (&_S3915)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3914, &_S3915, _S3905.differential_0);
    float _S3916 = _S3524 * _S3915.differential_0.x;
    float _S3917 = _S3530 * _S3915.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3918;
    (&_S3918)->primal_0 = _S3528;
    (&_S3918)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3919;
    (&_S3919)->primal_0 = _S3529;
    (&_S3919)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3918, &_S3919, _S3908.differential_0);
    Matrix<float, 3, 3>  _S3920 = transpose_0(_S3909.differential_0 + _S3914.differential_0 + _S3918.differential_0);
    float _S3921 = 2.0f * - _S3920.rows[int(2)].z;
    float _S3922 = 2.0f * _S3920.rows[int(2)].y;
    float _S3923 = 2.0f * _S3920.rows[int(2)].x;
    float _S3924 = 2.0f * _S3920.rows[int(1)].z;
    float _S3925 = 2.0f * - _S3920.rows[int(1)].y;
    float _S3926 = 2.0f * _S3920.rows[int(1)].x;
    float _S3927 = 2.0f * _S3920.rows[int(0)].z;
    float _S3928 = 2.0f * _S3920.rows[int(0)].y;
    float _S3929 = 2.0f * - _S3920.rows[int(0)].x;
    float _S3930 = - _S3926 + _S3928;
    float _S3931 = _S3923 + - _S3927;
    float _S3932 = - _S3922 + _S3924;
    float _S3933 = _S3922 + _S3924;
    float _S3934 = _S3923 + _S3927;
    float _S3935 = _S3926 + _S3928;
    float _S3936 = quat_26.w * (_S3925 + _S3929);
    float _S3937 = quat_26.z * (_S3921 + _S3929);
    float _S3938 = quat_26.y * (_S3921 + _S3925);
    float _S3939 = quat_26.x * _S3930 + quat_26.z * _S3933 + quat_26.y * _S3934 + _S3936 + _S3936;
    float _S3940 = quat_26.x * _S3931 + quat_26.w * _S3933 + quat_26.y * _S3935 + _S3937 + _S3937;
    float _S3941 = quat_26.x * _S3932 + quat_26.w * _S3934 + quat_26.z * _S3935 + _S3938 + _S3938;
    float _S3942 = quat_26.w * _S3930 + quat_26.z * _S3931 + quat_26.y * _S3932;
    float _S3943 = _S3913 + _S3916;
    float _S3944 = 0.5f * - _S3943;
    float _S3945 = _S3911 + _S3915.differential_0.y;
    DiffPair_float_0 _S3946;
    (&_S3946)->primal_0 = _S3525;
    (&_S3946)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3946, _S3945);
    float _S3947 = _S3944 + _S3946.differential_0;
    float _S3948 = _S3912 + _S3917 + _S3919.differential_0.x;
    DiffPair_float_0 _S3949;
    (&_S3949)->primal_0 = _S3523;
    (&_S3949)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3949, _S3948);
    float _S3950 = _S3944 + _S3949.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3951;
    (&_S3951)->primal_0 = mean_c_20;
    (&_S3951)->differential_0 = _S3669;
    s_bwd_length_impl_1(&_S3951, 0.0f);
    float3  _S3952 = _S3951.differential_0 + _S3672.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3953;
    (&_S3953)->primal_0 = R_26;
    (&_S3953)->differential_0 = _S3752;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3954;
    (&_S3954)->primal_0 = mean_25;
    (&_S3954)->differential_0 = _S3669;
    s_bwd_prop_mul_1(&_S3953, &_S3954, _S3952);
    float3  _S3955 = _S3900 + _S3903 + _S3906 + _S3952 + _S3755.differential_0;
    Matrix<float, 3, 3>  _S3956 = _S3901.differential_0 + _S3904.differential_0 + _S3907.differential_0 + _S3953.differential_0 + _S3756;
    float3  _S3957 = make_float3 (_S3950, _S3947, _S3943);
    float4  _S3958 = make_float4 (0.0f);
    *&((&_S3958)->w) = _S3939;
    *&((&_S3958)->z) = _S3940;
    *&((&_S3958)->y) = _S3941;
    *&((&_S3958)->x) = _S3942;
    float4  _S3959 = _S3958;
    float3  _S3960 = _S3902.differential_0 + _S3905.differential_0 + _S3908.differential_0 + _S3954.differential_0 + _S3750;
    *v_mean_8 = _S3960;
    *v_quat_7 = _S3959;
    *v_scale_7 = _S3957;
    *v_hardness_1 = _S3801;
    (*v_sh_coeffs_6)[int(0)] = _S3809;
    (*v_sh_coeffs_6)[int(1)] = _S3810;
    (*v_sh_coeffs_6)[int(2)] = _S3811;
    (*v_sh_coeffs_6)[int(3)] = _S3812;
    (*v_sh_coeffs_6)[int(4)] = _S3813;
    (*v_sh_coeffs_6)[int(5)] = _S3814;
    (*v_sh_coeffs_6)[int(6)] = _S3815;
    (*v_sh_coeffs_6)[int(7)] = _S3816;
    (*v_sh_coeffs_6)[int(8)] = _S3817;
    (*v_sh_coeffs_6)[int(9)] = _S3818;
    (*v_sh_coeffs_6)[int(10)] = _S3819;
    (*v_sh_coeffs_6)[int(11)] = _S3820;
    (*v_sh_coeffs_6)[int(12)] = _S3821;
    (*v_sh_coeffs_6)[int(13)] = _S3822;
    (*v_sh_coeffs_6)[int(14)] = _S3823;
    (*v_sh_coeffs_6)[int(15)] = _S3824;
    (*v_ch_coeffs_1)[int(0)] = _S3806;
    (*v_ch_coeffs_1)[int(1)] = _S3807;
    *v_R_8 = _S3956;
    *v_t_7 = _S3955;
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
        DiffPair_float_0 _S3961 = *dpx_17;
        float _S3962 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3962;
        float _S3963 = val_0 * (F32_log((_S3961.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3963;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_1)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3964 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3964))));
    float2  _S3965 = p_1 - v0_0;
    float2  _S3966 = normalize_1(e0_6);
    float2  _S3967 = p_1 - v1_0;
    float2  _S3968 = normalize_1(e1_6);
    float2  _S3969 = p_1 - v2_0;
    float2  _S3970 = normalize_1(e2_2);
    float _S3971 = hardness_6.x;
    float _S3972 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3965.x * _S3966.y - _S3965.y * _S3966.x)), (se_0 * (_S3967.x * _S3968.y - _S3967.y * _S3968.x))))), (se_0 * (_S3969.x * _S3970.y - _S3969.y * _S3970.x)))) / ((F32_abs((_S3964))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3972))));
    float _S3973;
    if(a_1 <= 0.0f)
    {
        _S3973 = 0.0f;
    }
    else
    {
        _S3973 = (F32_min(((F32_pow((a_1), (_S3972)))), (0.99900001287460327f)));
    }
    return _S3971 * _S3973;
}

inline __device__ float s_primal_ctx_abs_0(float _S3974)
{
    return (F32_abs((_S3974)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3975, float _S3976, float _S3977)
{
    return clamp_0(_S3975, _S3976, _S3977);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3978)
{
    return (F32_exp2((_S3978)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3979, float _S3980)
{
    return (F32_pow((_S3979), (_S3980)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3981, DiffPair_float_0 * _S3982, float _S3983)
{
    _d_pow_0(_S3981, _S3982, _S3983);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3984, DiffPair_float_0 * _S3985, DiffPair_float_0 * _S3986, float _S3987)
{
    _d_clamp_0(_S3984, _S3985, _S3986, _S3987);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_10)
{
    float _S3988 = length_0((*dpx_18).primal_0);
    float2  _S3989 = (*dpx_18).primal_0 * _s_dOut_10;
    float2  _S3990 = make_float2 (1.0f / _S3988) * _s_dOut_10;
    float _S3991 = - ((_S3989.x + _S3989.y) / (_S3988 * _S3988));
    float2  _S3992 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3993;
    (&_S3993)->primal_0 = (*dpx_18).primal_0;
    (&_S3993)->differential_0 = _S3992;
    s_bwd_length_impl_0(&_S3993, _S3991);
    float2  _S3994 = _S3990 + _S3993.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S3994;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3995, float2  _S3996)
{
    s_bwd_prop_normalize_impl_1(_S3995, _S3996);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_2, float _s_dOut_11)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S3997 = e0_7.x;
    float _S3998 = e1_7.y;
    float _S3999 = e0_7.y;
    float _S4000 = e1_7.x;
    float _S4001 = _S3997 * _S3998 - _S3999 * _S4000;
    float se_1 = float((F32_sign((_S4001))));
    float2  _S4002 = p_2 - (*dpv0_0).primal_0;
    float2  _S4003 = normalize_1(e0_7);
    float _S4004 = _S4002.x;
    float _S4005 = _S4003.y;
    float _S4006 = _S4002.y;
    float _S4007 = _S4003.x;
    float de0_0 = se_1 * (_S4004 * _S4005 - _S4006 * _S4007);
    float2  _S4008 = p_2 - (*dpv1_0).primal_0;
    float2  _S4009 = normalize_1(e1_7);
    float _S4010 = _S4008.x;
    float _S4011 = _S4009.y;
    float _S4012 = _S4008.y;
    float _S4013 = _S4009.x;
    float de1_0 = se_1 * (_S4010 * _S4011 - _S4012 * _S4013);
    float2  _S4014 = p_2 - (*dpv2_0).primal_0;
    float2  _S4015 = normalize_1(e2_3);
    float _S4016 = _S4014.x;
    float _S4017 = _S4015.y;
    float _S4018 = _S4014.y;
    float _S4019 = _S4015.x;
    float de2_0 = se_1 * (_S4016 * _S4017 - _S4018 * _S4019);
    float _S4020 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4021 = s_primal_ctx_max_0(_S4020, de2_0);
    float _S4022 = s_primal_ctx_abs_0(_S4001);
    float _S4023 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4022 / _S4023;
    float _S4024 = _S4023 * _S4023;
    float _S4025 = (*dphardness_0).primal_0.x;
    float _S4026 = (*dphardness_0).primal_0.y;
    float _S4027 = dmax_0 * dmax_0;
    float _S4028 = 1.0f + _S4021 / dmax_0;
    float _S4029 = 1.0f - s_primal_ctx_clamp_0(_S4026, 0.00499999988824129f, 0.98000001907348633f);
    float _S4030 = -1.0f / _S4029;
    float _S4031 = _S4029 * _S4029;
    float _S4032 = 1.0f - s_primal_ctx_exp2_0(_S4030);
    float a_2 = 1.0f - _S4028 * _S4032;
    bool _S4033 = a_2 <= 0.0f;
    float _S4034;
    float _S4035;
    if(_S4033)
    {
        _S4034 = 0.0f;
        _S4035 = 0.0f;
    }
    else
    {
        float _S4036 = s_primal_ctx_pow_0(a_2, _S4029);
        _S4034 = s_primal_ctx_min_0(_S4036, 0.99900001287460327f);
        _S4035 = _S4036;
    }
    float _S4037 = _S4025 * _s_dOut_11;
    float _S4038 = _S4034 * _s_dOut_11;
    if(_S4033)
    {
        _S4034 = 0.0f;
        _S4035 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4039;
        (&_S4039)->primal_0 = _S4035;
        (&_S4039)->differential_0 = 0.0f;
        DiffPair_float_0 _S4040;
        (&_S4040)->primal_0 = 0.99900001287460327f;
        (&_S4040)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4039, &_S4040, _S4037);
        DiffPair_float_0 _S4041;
        (&_S4041)->primal_0 = a_2;
        (&_S4041)->differential_0 = 0.0f;
        DiffPair_float_0 _S4042;
        (&_S4042)->primal_0 = _S4029;
        (&_S4042)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4041, &_S4042, _S4039.differential_0);
        _S4034 = _S4041.differential_0;
        _S4035 = _S4042.differential_0;
    }
    float _S4043 = - _S4034;
    float _S4044 = _S4032 * _S4043;
    float _S4045 = - (_S4028 * _S4043);
    DiffPair_float_0 _S4046;
    (&_S4046)->primal_0 = _S4030;
    (&_S4046)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4046, _S4045);
    float _S4047 = - (-1.0f * - (_S4046.differential_0 / _S4031) + _S4035);
    float _S4048 = _S4044 / _S4027;
    float s_diff_dmax_T_0 = _S4021 * - _S4048;
    float _S4049 = dmax_0 * _S4048;
    DiffPair_float_0 _S4050;
    (&_S4050)->primal_0 = _S4026;
    (&_S4050)->differential_0 = 0.0f;
    DiffPair_float_0 _S4051;
    (&_S4051)->primal_0 = 0.00499999988824129f;
    (&_S4051)->differential_0 = 0.0f;
    DiffPair_float_0 _S4052;
    (&_S4052)->primal_0 = 0.98000001907348633f;
    (&_S4052)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4050, &_S4051, &_S4052, _S4047);
    float _S4053 = s_diff_dmax_T_0 / _S4024;
    float _S4054 = _S4022 * - _S4053;
    float _S4055 = _S4023 * _S4053;
    float2  _S4056 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4057;
    (&_S4057)->primal_0 = e2_3;
    (&_S4057)->differential_0 = _S4056;
    s_bwd_length_impl_0(&_S4057, _S4054);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4058;
    (&_S4058)->primal_0 = e1_7;
    (&_S4058)->differential_0 = _S4056;
    s_bwd_length_impl_0(&_S4058, _S4054);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4059;
    (&_S4059)->primal_0 = e0_7;
    (&_S4059)->differential_0 = _S4056;
    s_bwd_length_impl_0(&_S4059, _S4054);
    DiffPair_float_0 _S4060;
    (&_S4060)->primal_0 = _S4001;
    (&_S4060)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4060, _S4055);
    DiffPair_float_0 _S4061;
    (&_S4061)->primal_0 = _S4020;
    (&_S4061)->differential_0 = 0.0f;
    DiffPair_float_0 _S4062;
    (&_S4062)->primal_0 = de2_0;
    (&_S4062)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4061, &_S4062, _S4049);
    DiffPair_float_0 _S4063;
    (&_S4063)->primal_0 = de0_0;
    (&_S4063)->differential_0 = 0.0f;
    DiffPair_float_0 _S4064;
    (&_S4064)->primal_0 = de1_0;
    (&_S4064)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4063, &_S4064, _S4061.differential_0);
    float _S4065 = se_1 * _S4062.differential_0;
    float _S4066 = - _S4065;
    float _S4067 = _S4019 * _S4066;
    float _S4068 = _S4017 * _S4065;
    float2  _S4069 = make_float2 (_S4018 * _S4066, _S4016 * _S4065);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4070;
    (&_S4070)->primal_0 = e2_3;
    (&_S4070)->differential_0 = _S4056;
    s_bwd_normalize_impl_1(&_S4070, _S4069);
    float2  _S4071 = - make_float2 (_S4068, _S4067);
    float _S4072 = se_1 * _S4064.differential_0;
    float _S4073 = - _S4072;
    float _S4074 = _S4013 * _S4073;
    float _S4075 = _S4011 * _S4072;
    float2  _S4076 = make_float2 (_S4012 * _S4073, _S4010 * _S4072);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4077;
    (&_S4077)->primal_0 = e1_7;
    (&_S4077)->differential_0 = _S4056;
    s_bwd_normalize_impl_1(&_S4077, _S4076);
    float2  _S4078 = - make_float2 (_S4075, _S4074);
    float _S4079 = se_1 * _S4063.differential_0;
    float _S4080 = - _S4079;
    float _S4081 = _S4007 * _S4080;
    float _S4082 = _S4005 * _S4079;
    float2  _S4083 = make_float2 (_S4006 * _S4080, _S4004 * _S4079);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4084;
    (&_S4084)->primal_0 = e0_7;
    (&_S4084)->differential_0 = _S4056;
    s_bwd_normalize_impl_1(&_S4084, _S4083);
    float2  _S4085 = - make_float2 (_S4082, _S4081);
    float _S4086 = - _S4060.differential_0;
    float2  _S4087 = _S4057.differential_0 + _S4070.differential_0;
    float2  _S4088 = - _S4087;
    float2  _S4089 = _S4058.differential_0 + _S4077.differential_0 + make_float2 (_S3999 * _S4086, _S3997 * _S4060.differential_0);
    float2  _S4090 = - _S4089;
    float2  _S4091 = _S4059.differential_0 + _S4084.differential_0 + make_float2 (_S3998 * _S4060.differential_0, _S4000 * _S4086);
    float2  _S4092 = - _S4091;
    float2  _S4093 = make_float2 (_S4038, _S4050.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4093;
    float2  _S4094 = _S4071 + _S4088 + _S4089;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4094;
    float2  _S4095 = _S4078 + _S4090 + _S4091;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4095;
    float2  _S4096 = _S4085 + _S4087 + _S4092;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4096;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4097, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4098, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4099, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4100, float2  _S4101, float _S4102)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4097, _S4098, _S4099, _S4100, _S4101, _S4102);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_3, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4103 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4103;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4103;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4103;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4103;
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
    float2  _S4104 = p_4 - v0_2;
    float2  _S4105 = p_4 - v1_2;
    float2  _S4106 = p_4 - v2_2;
    float _S4107 = e0_8.x;
    float _S4108 = e1_8.y;
    float _S4109 = e0_8.y;
    float _S4110 = e1_8.x;
    float _S4111 = _S4107 * _S4108 - _S4109 * _S4110;
    float se_2 = float((F32_sign((_S4111))));
    float _S4112 = hardness_8.x;
    float _S4113 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4104.x * _S4109 - _S4104.y * _S4107)), (se_2 * (_S4105.x * _S4108 - _S4105.y * _S4110))))), (se_2 * (_S4106.x * e2_4.y - _S4106.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4104 - e0_8 * make_float2 (clamp_0(dot_1(_S4104, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4105 - e1_8 * make_float2 (clamp_0(dot_1(_S4105, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4106 - e2_4 * make_float2 (clamp_0(dot_1(_S4106, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4111))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4113))));
    float _S4114;
    if(a_3 <= 0.0f)
    {
        _S4114 = 0.0f;
    }
    else
    {
        _S4114 = (F32_min(((F32_pow((a_3), (_S4113)))), (0.99900001287460327f)));
    }
    return _S4112 * _S4114;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4115, float2  _S4116)
{
    return dot_1(_S4115, _S4116);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4117, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4118, float _S4119)
{
    _d_dot_1(_S4117, _S4118, _S4119);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_5, float _s_dOut_12)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4120 = p_5 - (*dpv0_1).primal_0;
    float _S4121 = s_primal_ctx_dot_1(_S4120, e0_9);
    float _S4122 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4123 = _S4121 / _S4122;
    float _S4124 = _S4122 * _S4122;
    float _S4125 = s_primal_ctx_clamp_0(_S4123, 0.0f, 1.0f);
    float2  _S4126 = make_float2 (_S4125);
    float2  _S4127 = _S4120 - e0_9 * make_float2 (_S4125);
    float _S4128 = length_0(_S4127);
    float2  _S4129 = p_5 - (*dpv1_1).primal_0;
    float _S4130 = s_primal_ctx_dot_1(_S4129, e1_9);
    float _S4131 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4132 = _S4130 / _S4131;
    float _S4133 = _S4131 * _S4131;
    float _S4134 = s_primal_ctx_clamp_0(_S4132, 0.0f, 1.0f);
    float2  _S4135 = make_float2 (_S4134);
    float2  _S4136 = _S4129 - e1_9 * make_float2 (_S4134);
    float _S4137 = length_0(_S4136);
    float2  _S4138 = p_5 - (*dpv2_1).primal_0;
    float _S4139 = s_primal_ctx_dot_1(_S4138, e2_5);
    float _S4140 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4141 = _S4139 / _S4140;
    float _S4142 = _S4140 * _S4140;
    float _S4143 = s_primal_ctx_clamp_0(_S4141, 0.0f, 1.0f);
    float2  _S4144 = make_float2 (_S4143);
    float2  _S4145 = _S4138 - e2_5 * make_float2 (_S4143);
    float _S4146 = length_0(_S4145);
    float _S4147 = e0_9.x;
    float _S4148 = e1_9.y;
    float _S4149 = e0_9.y;
    float _S4150 = e1_9.x;
    float _S4151 = _S4147 * _S4148 - _S4149 * _S4150;
    float se_3 = float((F32_sign((_S4151))));
    float _S4152 = _S4120.x;
    float _S4153 = _S4120.y;
    float s0_0 = se_3 * (_S4152 * _S4149 - _S4153 * _S4147);
    float _S4154 = _S4129.x;
    float _S4155 = _S4129.y;
    float s1_0 = se_3 * (_S4154 * _S4148 - _S4155 * _S4150);
    float _S4156 = _S4138.x;
    float _S4157 = e2_5.y;
    float _S4158 = _S4138.y;
    float _S4159 = e2_5.x;
    float s2_0 = se_3 * (_S4156 * _S4157 - _S4158 * _S4159);
    float _S4160 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4160, s2_0)))));
    float _S4161 = s_primal_ctx_min_0(_S4128, _S4137);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4161, _S4146);
    float _S4162 = s_primal_ctx_abs_0(_S4151);
    float _S4163 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4162 / _S4163;
    float _S4164 = _S4163 * _S4163;
    float _S4165 = (*dphardness_1).primal_0.x;
    float _S4166 = (*dphardness_1).primal_0.y;
    float _S4167 = dmax_1 * dmax_1;
    float _S4168 = 1.0f + dv_0 / dmax_1;
    float _S4169 = 1.0f - s_primal_ctx_clamp_0(_S4166, 0.00499999988824129f, 0.98000001907348633f);
    float _S4170 = -1.0f / _S4169;
    float _S4171 = _S4169 * _S4169;
    float _S4172 = 1.0f - s_primal_ctx_exp2_0(_S4170);
    float a_4 = 1.0f - _S4168 * _S4172;
    bool _S4173 = a_4 <= 0.0f;
    float _S4174;
    float _S4175;
    if(_S4173)
    {
        _S4174 = 0.0f;
        _S4175 = 0.0f;
    }
    else
    {
        float _S4176 = s_primal_ctx_pow_0(a_4, _S4169);
        _S4174 = s_primal_ctx_min_0(_S4176, 0.99900001287460327f);
        _S4175 = _S4176;
    }
    float _S4177 = _S4165 * _s_dOut_12;
    float _S4178 = _S4174 * _s_dOut_12;
    if(_S4173)
    {
        _S4174 = 0.0f;
        _S4175 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4179;
        (&_S4179)->primal_0 = _S4175;
        (&_S4179)->differential_0 = 0.0f;
        DiffPair_float_0 _S4180;
        (&_S4180)->primal_0 = 0.99900001287460327f;
        (&_S4180)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4179, &_S4180, _S4177);
        DiffPair_float_0 _S4181;
        (&_S4181)->primal_0 = a_4;
        (&_S4181)->differential_0 = 0.0f;
        DiffPair_float_0 _S4182;
        (&_S4182)->primal_0 = _S4169;
        (&_S4182)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4181, &_S4182, _S4179.differential_0);
        _S4174 = _S4181.differential_0;
        _S4175 = _S4182.differential_0;
    }
    float _S4183 = - _S4174;
    float _S4184 = _S4172 * _S4183;
    float _S4185 = - (_S4168 * _S4183);
    DiffPair_float_0 _S4186;
    (&_S4186)->primal_0 = _S4170;
    (&_S4186)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4186, _S4185);
    float _S4187 = - (-1.0f * - (_S4186.differential_0 / _S4171) + _S4175);
    float _S4188 = _S4184 / _S4167;
    float s_diff_dmax_T_1 = dv_0 * - _S4188;
    float s_diff_dv_T_0 = dmax_1 * _S4188;
    DiffPair_float_0 _S4189;
    (&_S4189)->primal_0 = _S4166;
    (&_S4189)->differential_0 = 0.0f;
    DiffPair_float_0 _S4190;
    (&_S4190)->primal_0 = 0.00499999988824129f;
    (&_S4190)->differential_0 = 0.0f;
    DiffPair_float_0 _S4191;
    (&_S4191)->primal_0 = 0.98000001907348633f;
    (&_S4191)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4189, &_S4190, &_S4191, _S4187);
    float _S4192 = s_diff_dmax_T_1 / _S4164;
    float _S4193 = _S4162 * - _S4192;
    float _S4194 = _S4163 * _S4192;
    float2  _S4195 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4196;
    (&_S4196)->primal_0 = e2_5;
    (&_S4196)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4196, _S4193);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4197;
    (&_S4197)->primal_0 = e1_9;
    (&_S4197)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4197, _S4193);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4198;
    (&_S4198)->primal_0 = e0_9;
    (&_S4198)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4198, _S4193);
    DiffPair_float_0 _S4199;
    (&_S4199)->primal_0 = _S4151;
    (&_S4199)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4199, _S4194);
    float _S4200 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4201;
    (&_S4201)->primal_0 = _S4161;
    (&_S4201)->differential_0 = 0.0f;
    DiffPair_float_0 _S4202;
    (&_S4202)->primal_0 = _S4146;
    (&_S4202)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4201, &_S4202, _S4200);
    DiffPair_float_0 _S4203;
    (&_S4203)->primal_0 = _S4128;
    (&_S4203)->differential_0 = 0.0f;
    DiffPair_float_0 _S4204;
    (&_S4204)->primal_0 = _S4137;
    (&_S4204)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4203, &_S4204, _S4201.differential_0);
    DiffPair_float_0 _S4205;
    (&_S4205)->primal_0 = _S4160;
    (&_S4205)->differential_0 = 0.0f;
    DiffPair_float_0 _S4206;
    (&_S4206)->primal_0 = s2_0;
    (&_S4206)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4205, &_S4206, 0.0f);
    DiffPair_float_0 _S4207;
    (&_S4207)->primal_0 = s0_0;
    (&_S4207)->differential_0 = 0.0f;
    DiffPair_float_0 _S4208;
    (&_S4208)->primal_0 = s1_0;
    (&_S4208)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4207, &_S4208, _S4205.differential_0);
    float _S4209 = se_3 * _S4206.differential_0;
    float _S4210 = - _S4209;
    float _S4211 = _S4158 * _S4210;
    float _S4212 = _S4159 * _S4210;
    float _S4213 = _S4156 * _S4209;
    float _S4214 = _S4157 * _S4209;
    float _S4215 = se_3 * _S4208.differential_0;
    float _S4216 = - _S4215;
    float _S4217 = _S4150 * _S4216;
    float _S4218 = _S4148 * _S4215;
    float _S4219 = se_3 * _S4207.differential_0;
    float _S4220 = - _S4219;
    float _S4221 = _S4147 * _S4220;
    float _S4222 = _S4149 * _S4219;
    float _S4223 = - _S4199.differential_0;
    float _S4224 = _S4155 * _S4216 + _S4149 * _S4223;
    float _S4225 = _S4152 * _S4219 + _S4150 * _S4223;
    float _S4226 = _S4154 * _S4215 + _S4147 * _S4199.differential_0;
    float _S4227 = _S4153 * _S4220 + _S4148 * _S4199.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4228;
    (&_S4228)->primal_0 = _S4145;
    (&_S4228)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4228, _S4202.differential_0);
    float2  _S4229 = - _S4228.differential_0;
    float2  _S4230 = e2_5 * _S4229;
    float2  _S4231 = _S4144 * _S4229;
    float _S4232 = _S4230.x + _S4230.y;
    DiffPair_float_0 _S4233;
    (&_S4233)->primal_0 = _S4141;
    (&_S4233)->differential_0 = 0.0f;
    DiffPair_float_0 _S4234;
    (&_S4234)->primal_0 = 0.0f;
    (&_S4234)->differential_0 = 0.0f;
    DiffPair_float_0 _S4235;
    (&_S4235)->primal_0 = 1.0f;
    (&_S4235)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4233, &_S4234, &_S4235, _S4232);
    float _S4236 = _S4233.differential_0 / _S4142;
    float _S4237 = _S4139 * - _S4236;
    float _S4238 = _S4140 * _S4236;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4239;
    (&_S4239)->primal_0 = e2_5;
    (&_S4239)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4240;
    (&_S4240)->primal_0 = e2_5;
    (&_S4240)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4239, &_S4240, _S4237);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4241;
    (&_S4241)->primal_0 = _S4138;
    (&_S4241)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4242;
    (&_S4242)->primal_0 = e2_5;
    (&_S4242)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4241, &_S4242, _S4238);
    float2  _S4243 = - (_S4228.differential_0 + _S4241.differential_0 + make_float2 (_S4214, _S4212));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4244;
    (&_S4244)->primal_0 = _S4136;
    (&_S4244)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4244, _S4204.differential_0);
    float2  _S4245 = - _S4244.differential_0;
    float2  _S4246 = e1_9 * _S4245;
    float2  _S4247 = _S4135 * _S4245;
    float _S4248 = _S4246.x + _S4246.y;
    DiffPair_float_0 _S4249;
    (&_S4249)->primal_0 = _S4132;
    (&_S4249)->differential_0 = 0.0f;
    DiffPair_float_0 _S4250;
    (&_S4250)->primal_0 = 0.0f;
    (&_S4250)->differential_0 = 0.0f;
    DiffPair_float_0 _S4251;
    (&_S4251)->primal_0 = 1.0f;
    (&_S4251)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4249, &_S4250, &_S4251, _S4248);
    float _S4252 = _S4249.differential_0 / _S4133;
    float _S4253 = _S4130 * - _S4252;
    float _S4254 = _S4131 * _S4252;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4255;
    (&_S4255)->primal_0 = e1_9;
    (&_S4255)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4256;
    (&_S4256)->primal_0 = e1_9;
    (&_S4256)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4255, &_S4256, _S4253);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4257;
    (&_S4257)->primal_0 = _S4129;
    (&_S4257)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4258;
    (&_S4258)->primal_0 = e1_9;
    (&_S4258)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4257, &_S4258, _S4254);
    float2  _S4259 = - (_S4244.differential_0 + _S4257.differential_0 + make_float2 (_S4218, _S4217));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4260;
    (&_S4260)->primal_0 = _S4127;
    (&_S4260)->differential_0 = _S4195;
    s_bwd_length_impl_0(&_S4260, _S4203.differential_0);
    float2  _S4261 = - _S4260.differential_0;
    float2  _S4262 = e0_9 * _S4261;
    float2  _S4263 = _S4126 * _S4261;
    float _S4264 = _S4262.x + _S4262.y;
    DiffPair_float_0 _S4265;
    (&_S4265)->primal_0 = _S4123;
    (&_S4265)->differential_0 = 0.0f;
    DiffPair_float_0 _S4266;
    (&_S4266)->primal_0 = 0.0f;
    (&_S4266)->differential_0 = 0.0f;
    DiffPair_float_0 _S4267;
    (&_S4267)->primal_0 = 1.0f;
    (&_S4267)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4265, &_S4266, &_S4267, _S4264);
    float _S4268 = _S4265.differential_0 / _S4124;
    float _S4269 = _S4121 * - _S4268;
    float _S4270 = _S4122 * _S4268;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4271;
    (&_S4271)->primal_0 = e0_9;
    (&_S4271)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4272;
    (&_S4272)->primal_0 = e0_9;
    (&_S4272)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4271, &_S4272, _S4269);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4273;
    (&_S4273)->primal_0 = _S4120;
    (&_S4273)->differential_0 = _S4195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4274;
    (&_S4274)->primal_0 = e0_9;
    (&_S4274)->differential_0 = _S4195;
    s_bwd_prop_dot_1(&_S4273, &_S4274, _S4270);
    float2  _S4275 = - (_S4260.differential_0 + _S4273.differential_0 + make_float2 (_S4222, _S4221));
    float2  _S4276 = _S4196.differential_0 + _S4231 + _S4240.differential_0 + _S4239.differential_0 + _S4242.differential_0 + make_float2 (_S4211, _S4213);
    float2  _S4277 = - _S4276;
    float2  _S4278 = _S4197.differential_0 + _S4247 + _S4256.differential_0 + _S4255.differential_0 + _S4258.differential_0 + make_float2 (_S4224, _S4226);
    float2  _S4279 = - _S4278;
    float2  _S4280 = _S4198.differential_0 + _S4263 + _S4272.differential_0 + _S4271.differential_0 + _S4274.differential_0 + make_float2 (_S4227, _S4225);
    float2  _S4281 = - _S4280;
    float2  _S4282 = make_float2 (_S4178, _S4189.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4282;
    float2  _S4283 = _S4243 + _S4277 + _S4278;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4283;
    float2  _S4284 = _S4259 + _S4279 + _S4280;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4284;
    float2  _S4285 = _S4275 + _S4276 + _S4281;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4285;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4286, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4287, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4288, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4289, float2  _S4290, float _S4291)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4286, _S4287, _S4288, _S4289, _S4290, _S4291);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_6, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4292 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4292;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4292;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4292;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4292;
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
    float _S4293 = 0.3333333432674408f * dpdepth_1;
    float3  _S4294 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4295 = make_float3 (0.0f);
    float3  _S4296 = _S4295;
    *&((&_S4296)->z) = _S4293;
    *&((&_S4296)->y) = _S4293;
    *&((&_S4296)->x) = _S4293;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4296;
    FixedArray<float3 , 3>  _S4297;
    _S4297[int(0)] = _S4295;
    _S4297[int(1)] = _S4295;
    _S4297[int(2)] = _S4295;
    _S4297[int(2)] = _S4294;
    _S4297[int(1)] = _S4294;
    _S4297[int(0)] = _S4294;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4297;
    float2  _S4298 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4298;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4298;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4298;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4299, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4300, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4301, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4302, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4303, float2  _S4304, float3  _S4305, float _S4306)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4299, _S4300, _S4301, _S4302, _S4303, _S4304, _S4305, _S4306);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_9, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4307 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4307;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4307;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4307;
    float3  _S4308 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4309 = { _S4308, _S4308, _S4308 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4309;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4308;
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
        float _S4310 = mean_c_21.z;
        bool _S4311;
        if(_S4310 < near_plane_14)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4310 > far_plane_14;
        }
        if(_S4311)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4312 = scale_26.x;
        float sx_5 = (F32_exp((_S4312)));
        float _S4313 = scale_26.y;
        float sy_5 = (F32_exp((_S4313)));
        float sz_7 = scale_26.z - 0.5f * (_S4312 + _S4313);
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
        Matrix<float, 3, 3>  _S4314 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_48), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_48), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4314, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4314, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4314, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_6 = mul_0(R_27, vert0_3) + t_24;
        float3  vert1_c_6 = mul_0(R_27, vert1_3) + t_24;
        float3  vert2_c_6 = mul_0(R_27, vert2_3) + t_24;
        float _S4315 = vert0_c_6.z;
        float _S4316 = vert1_c_6.z;
        float _S4317 = vert2_c_6.z;
        if(_S4315 < near_plane_14)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4315 > far_plane_14;
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4316 < near_plane_14;
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4316 > far_plane_14;
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4317 < near_plane_14;
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = _S4317 > far_plane_14;
        }
        if(_S4311)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_4;
        for(;;)
        {
            float2  uv0_5 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4315);
            if(_S4315 < 0.0f)
            {
                _S4311 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4318 = camera_distortion_jac_0(uv0_5, dist_coeffs_34);
                _S4311 = !((F32_min((determinant_0(_S4318)), ((F32_min((_S4318.rows[int(0)].x), (_S4318.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4311)
            {
                uv0_4 = uv0_5;
                break;
            }
            float u_62 = uv0_5.x;
            float v_62 = uv0_5.y;
            float r2_62 = u_62 * u_62 + v_62 * v_62;
            float2  _S4319 = uv0_5 * make_float2 (1.0f + r2_62 * ((*dist_coeffs_34)[int(0)] + r2_62 * ((*dist_coeffs_34)[int(1)] + r2_62 * ((*dist_coeffs_34)[int(2)] + r2_62 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_62 * v_62 + (*dist_coeffs_34)[int(5)] * (r2_62 + 2.0f * u_62 * u_62) + (*dist_coeffs_34)[int(6)] * r2_62, 2.0f * (*dist_coeffs_34)[int(5)] * u_62 * v_62 + (*dist_coeffs_34)[int(4)] * (r2_62 + 2.0f * v_62 * v_62) + (*dist_coeffs_34)[int(7)] * r2_62);
            float2  _S4320 = _S4319 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4319.x + (*dist_coeffs_34)[int(9)] * _S4319.y, 0.0f);
            uv0_4 = make_float2 (fx_30 * _S4320.x + cx_25, fy_30 * _S4320.y + cy_25);
            break;
        }
        float2  uv1_4;
        bool all_valid_12 = true & (!_S4311);
        for(;;)
        {
            float2  uv1_5 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4316);
            if(_S4316 < 0.0f)
            {
                _S4311 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4321 = camera_distortion_jac_0(uv1_5, dist_coeffs_34);
                _S4311 = !((F32_min((determinant_0(_S4321)), ((F32_min((_S4321.rows[int(0)].x), (_S4321.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4311)
            {
                uv1_4 = uv1_5;
                break;
            }
            float u_63 = uv1_5.x;
            float v_63 = uv1_5.y;
            float r2_63 = u_63 * u_63 + v_63 * v_63;
            float2  _S4322 = uv1_5 * make_float2 (1.0f + r2_63 * ((*dist_coeffs_34)[int(0)] + r2_63 * ((*dist_coeffs_34)[int(1)] + r2_63 * ((*dist_coeffs_34)[int(2)] + r2_63 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_63 * v_63 + (*dist_coeffs_34)[int(5)] * (r2_63 + 2.0f * u_63 * u_63) + (*dist_coeffs_34)[int(6)] * r2_63, 2.0f * (*dist_coeffs_34)[int(5)] * u_63 * v_63 + (*dist_coeffs_34)[int(4)] * (r2_63 + 2.0f * v_63 * v_63) + (*dist_coeffs_34)[int(7)] * r2_63);
            float2  _S4323 = _S4322 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4322.x + (*dist_coeffs_34)[int(9)] * _S4322.y, 0.0f);
            uv1_4 = make_float2 (fx_30 * _S4323.x + cx_25, fy_30 * _S4323.y + cy_25);
            break;
        }
        float2  uv2_4;
        bool all_valid_13 = all_valid_12 & (!_S4311);
        for(;;)
        {
            float2  uv2_5 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4317);
            if(_S4317 < 0.0f)
            {
                _S4311 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4324 = camera_distortion_jac_0(uv2_5, dist_coeffs_34);
                _S4311 = !((F32_min((determinant_0(_S4324)), ((F32_min((_S4324.rows[int(0)].x), (_S4324.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4311)
            {
                uv2_4 = uv2_5;
                break;
            }
            float u_64 = uv2_5.x;
            float v_64 = uv2_5.y;
            float r2_64 = u_64 * u_64 + v_64 * v_64;
            float2  _S4325 = uv2_5 * make_float2 (1.0f + r2_64 * ((*dist_coeffs_34)[int(0)] + r2_64 * ((*dist_coeffs_34)[int(1)] + r2_64 * ((*dist_coeffs_34)[int(2)] + r2_64 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_64 * v_64 + (*dist_coeffs_34)[int(5)] * (r2_64 + 2.0f * u_64 * u_64) + (*dist_coeffs_34)[int(6)] * r2_64, 2.0f * (*dist_coeffs_34)[int(5)] * u_64 * v_64 + (*dist_coeffs_34)[int(4)] * (r2_64 + 2.0f * v_64 * v_64) + (*dist_coeffs_34)[int(7)] * r2_64);
            float2  _S4326 = _S4325 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4325.x + (*dist_coeffs_34)[int(9)] * _S4325.y, 0.0f);
            uv2_4 = make_float2 (fx_30 * _S4326.x + cx_25, fy_30 * _S4326.y + cy_25);
            break;
        }
        if(!(all_valid_13 & (!_S4311)))
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_10 = uv1_4 - uv0_4;
        float2  e1_10 = uv2_4 - uv1_4;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(uv0_4 - uv2_4)));
        float _S4327 = uv0_4.x;
        float _S4328 = uv1_4.x;
        float _S4329 = uv2_4.x;
        float xmax_7 = (F32_max(((F32_max((_S4327), (_S4328)))), (_S4329))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4327), (_S4328)))), (_S4329))) - offset_4;
        float _S4330 = uv0_4.y;
        float _S4331 = uv1_4.y;
        float _S4332 = uv2_4.y;
        float ymax_7 = (F32_max(((F32_max((_S4330), (_S4331)))), (_S4332))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4330), (_S4331)))), (_S4332))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = xmin_7 >= float(image_width_21);
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = ymax_7 <= 0.0f;
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            _S4311 = ymin_7 >= float(image_height_21);
        }
        if(_S4311)
        {
            _S4311 = true;
        }
        else
        {
            if(_S4310 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4311 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4311 = false;
                }
                if(_S4311)
                {
                    _S4311 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4311 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4311 = false;
                    }
                }
            }
            else
            {
                _S4311 = false;
            }
        }
        if(_S4311)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4333 = mean_26 - - mul_0(transpose_0(R_27), t_24);
        float _S4334 = _S4333.x;
        float _S4335 = _S4333.y;
        float _S4336 = _S4333.z;
        float norm_14 = (F32_sqrt((_S4334 * _S4334 + _S4335 * _S4335 + _S4336 * _S4336)));
        float x_54 = _S4334 / norm_14;
        float y_24 = _S4335 / norm_14;
        float z_21 = _S4336 / norm_14;
        float z2_49 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_54 * x_54 - y_24 * y_24;
        float fS1_21 = 2.0f * x_54 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_49 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_54) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_49 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_54) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_54 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_49 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_54) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_54 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4337 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4337);
        float3  _S4338 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4339 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4338 + _S4339 + make_float3 (0.5f), _S4337);
        (*rgbs_0)[int(2)] = max_0(_S4338 - _S4339 + make_float3 (0.5f), _S4337);
        (*verts_0)[int(0)] = vert0_3;
        (*verts_0)[int(1)] = vert1_3;
        (*verts_0)[int(2)] = vert2_3;
        float3  _S4340 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4340 * make_float3 (float(- (F32_sign((dot_0(_S4340, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_28, float3  t_25, float fx_31, float fy_31, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    bool _S4341;
    bool _S4342;
    bool _S4343;
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_28, mean_27) + t_25;
        float _S4344 = length_1(mean_c_22);
        bool _S4345;
        if(_S4344 < near_plane_15)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4344 > far_plane_15;
        }
        if(_S4345)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4346 = scale_27.x;
        float sx_6 = (F32_exp((_S4346)));
        float _S4347 = scale_27.y;
        float sy_6 = (F32_exp((_S4347)));
        float sz_8 = scale_27.z - 0.5f * (_S4346 + _S4347);
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
        Matrix<float, 3, 3>  _S4348 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_50), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_50), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
        float3  vert0_4 = mul_0(_S4348, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
        float3  vert1_4 = mul_0(_S4348, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
        float3  vert2_4 = mul_0(_S4348, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
        float3  vert0_c_7 = mul_0(R_28, vert0_4) + t_25;
        float3  vert1_c_7 = mul_0(R_28, vert1_4) + t_25;
        float3  vert2_c_7 = mul_0(R_28, vert2_4) + t_25;
        float _S4349 = length_1(vert0_c_7);
        float _S4350 = length_1(vert1_c_7);
        float _S4351 = length_1(vert2_c_7);
        if(_S4349 < near_plane_15)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4349 > far_plane_15;
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4350 < near_plane_15;
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4350 > far_plane_15;
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4351 < near_plane_15;
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = _S4351 > far_plane_15;
        }
        if(_S4345)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_6;
        float k_10;
        for(;;)
        {
            float2  _S4352 = float2 {vert0_c_7.x, vert0_c_7.y};
            float r_30 = length_0(_S4352);
            float _S4353 = vert0_c_7.z;
            float theta_26 = (F32_atan2((r_30), (_S4353)));
            if(theta_26 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_26 * theta_26 / 3.0f) / _S4353;
            }
            else
            {
                k_10 = theta_26 / r_30;
            }
            float2  uv0_7 = _S4352 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4354 = camera_distortion_jac_0(uv0_7, dist_coeffs_35);
            bool _S4355 = !((F32_min((determinant_0(_S4354)), ((F32_min((_S4354.rows[int(0)].x), (_S4354.rows[int(1)].y)))))) > 0.0f);
            _S4341 = _S4355;
            if(_S4355)
            {
                uv0_6 = uv0_7;
                break;
            }
            float u_65 = uv0_7.x;
            float v_65 = uv0_7.y;
            float r2_65 = u_65 * u_65 + v_65 * v_65;
            float2  _S4356 = uv0_7 * make_float2 (1.0f + r2_65 * ((*dist_coeffs_35)[int(0)] + r2_65 * ((*dist_coeffs_35)[int(1)] + r2_65 * ((*dist_coeffs_35)[int(2)] + r2_65 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_65 * v_65 + (*dist_coeffs_35)[int(5)] * (r2_65 + 2.0f * u_65 * u_65) + (*dist_coeffs_35)[int(6)] * r2_65, 2.0f * (*dist_coeffs_35)[int(5)] * u_65 * v_65 + (*dist_coeffs_35)[int(4)] * (r2_65 + 2.0f * v_65 * v_65) + (*dist_coeffs_35)[int(7)] * r2_65);
            float2  _S4357 = _S4356 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4356.x + (*dist_coeffs_35)[int(9)] * _S4356.y, 0.0f);
            uv0_6 = make_float2 (fx_31 * _S4357.x + cx_26, fy_31 * _S4357.y + cy_26);
            break;
        }
        float2  uv1_6;
        bool all_valid_14 = true & (!_S4341);
        for(;;)
        {
            float2  _S4358 = float2 {vert1_c_7.x, vert1_c_7.y};
            float r_31 = length_0(_S4358);
            float _S4359 = vert1_c_7.z;
            float theta_27 = (F32_atan2((r_31), (_S4359)));
            if(theta_27 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_27 * theta_27 / 3.0f) / _S4359;
            }
            else
            {
                k_10 = theta_27 / r_31;
            }
            float2  uv1_7 = _S4358 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4360 = camera_distortion_jac_0(uv1_7, dist_coeffs_35);
            bool _S4361 = !((F32_min((determinant_0(_S4360)), ((F32_min((_S4360.rows[int(0)].x), (_S4360.rows[int(1)].y)))))) > 0.0f);
            _S4342 = _S4361;
            if(_S4361)
            {
                uv1_6 = uv1_7;
                break;
            }
            float u_66 = uv1_7.x;
            float v_66 = uv1_7.y;
            float r2_66 = u_66 * u_66 + v_66 * v_66;
            float2  _S4362 = uv1_7 * make_float2 (1.0f + r2_66 * ((*dist_coeffs_35)[int(0)] + r2_66 * ((*dist_coeffs_35)[int(1)] + r2_66 * ((*dist_coeffs_35)[int(2)] + r2_66 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_66 * v_66 + (*dist_coeffs_35)[int(5)] * (r2_66 + 2.0f * u_66 * u_66) + (*dist_coeffs_35)[int(6)] * r2_66, 2.0f * (*dist_coeffs_35)[int(5)] * u_66 * v_66 + (*dist_coeffs_35)[int(4)] * (r2_66 + 2.0f * v_66 * v_66) + (*dist_coeffs_35)[int(7)] * r2_66);
            float2  _S4363 = _S4362 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4362.x + (*dist_coeffs_35)[int(9)] * _S4362.y, 0.0f);
            uv1_6 = make_float2 (fx_31 * _S4363.x + cx_26, fy_31 * _S4363.y + cy_26);
            break;
        }
        float2  uv2_6;
        bool all_valid_15 = all_valid_14 & (!_S4342);
        for(;;)
        {
            float2  _S4364 = float2 {vert2_c_7.x, vert2_c_7.y};
            float r_32 = length_0(_S4364);
            float _S4365 = vert2_c_7.z;
            float theta_28 = (F32_atan2((r_32), (_S4365)));
            if(theta_28 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_28 * theta_28 / 3.0f) / _S4365;
            }
            else
            {
                k_10 = theta_28 / r_32;
            }
            float2  uv2_7 = _S4364 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4366 = camera_distortion_jac_0(uv2_7, dist_coeffs_35);
            bool _S4367 = !((F32_min((determinant_0(_S4366)), ((F32_min((_S4366.rows[int(0)].x), (_S4366.rows[int(1)].y)))))) > 0.0f);
            _S4343 = _S4367;
            if(_S4367)
            {
                uv2_6 = uv2_7;
                break;
            }
            float u_67 = uv2_7.x;
            float v_67 = uv2_7.y;
            float r2_67 = u_67 * u_67 + v_67 * v_67;
            float2  _S4368 = uv2_7 * make_float2 (1.0f + r2_67 * ((*dist_coeffs_35)[int(0)] + r2_67 * ((*dist_coeffs_35)[int(1)] + r2_67 * ((*dist_coeffs_35)[int(2)] + r2_67 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_67 * v_67 + (*dist_coeffs_35)[int(5)] * (r2_67 + 2.0f * u_67 * u_67) + (*dist_coeffs_35)[int(6)] * r2_67, 2.0f * (*dist_coeffs_35)[int(5)] * u_67 * v_67 + (*dist_coeffs_35)[int(4)] * (r2_67 + 2.0f * v_67 * v_67) + (*dist_coeffs_35)[int(7)] * r2_67);
            float2  _S4369 = _S4368 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4368.x + (*dist_coeffs_35)[int(9)] * _S4368.y, 0.0f);
            uv2_6 = make_float2 (fx_31 * _S4369.x + cx_26, fy_31 * _S4369.y + cy_26);
            break;
        }
        if(!(all_valid_15 & (!_S4343)))
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_11 = uv1_6 - uv0_6;
        float2  e1_11 = uv2_6 - uv1_6;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(uv0_6 - uv2_6)));
        float _S4370 = uv0_6.x;
        float _S4371 = uv1_6.x;
        float _S4372 = uv2_6.x;
        float xmax_8 = (F32_max(((F32_max((_S4370), (_S4371)))), (_S4372))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4370), (_S4371)))), (_S4372))) - offset_5;
        float _S4373 = uv0_6.y;
        float _S4374 = uv1_6.y;
        float _S4375 = uv2_6.y;
        float ymax_8 = (F32_max(((F32_max((_S4373), (_S4374)))), (_S4375))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4373), (_S4374)))), (_S4375))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = xmin_8 >= float(image_width_22);
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = ymax_8 <= 0.0f;
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            _S4345 = ymin_8 >= float(image_height_22);
        }
        if(_S4345)
        {
            _S4345 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4345 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4345 = false;
                }
                if(_S4345)
                {
                    _S4345 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4345 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4345 = false;
                    }
                }
            }
            else
            {
                _S4345 = false;
            }
        }
        if(_S4345)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4376 = mean_27 - - mul_0(transpose_0(R_28), t_25);
        float _S4377 = _S4376.x;
        float _S4378 = _S4376.y;
        float _S4379 = _S4376.z;
        float norm_15 = (F32_sqrt((_S4377 * _S4377 + _S4378 * _S4378 + _S4379 * _S4379)));
        float x_56 = _S4377 / norm_15;
        float y_25 = _S4378 / norm_15;
        float z_22 = _S4379 / norm_15;
        float z2_51 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_56 * x_56 - y_25 * y_25;
        float fS1_22 = 2.0f * x_56 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_51 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_56) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_51 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_56) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_56 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_51 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_56) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_56 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4380 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4380);
        float3  _S4381 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4382 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4381 + _S4382 + make_float3 (0.5f), _S4380);
        (*rgbs_1)[int(2)] = max_0(_S4381 - _S4382 + make_float3 (0.5f), _S4380);
        (*verts_1)[int(0)] = vert0_4;
        (*verts_1)[int(1)] = vert1_4;
        (*verts_1)[int(2)] = vert2_4;
        float3  _S4383 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4383 * make_float3 (float(- (F32_sign((dot_0(_S4383, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_29, float3  t_26, float fx_32, float fy_32, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_29, mean_28) + t_26;
    float _S4384 = scale_28.x;
    float sx_7 = (F32_exp((_S4384)));
    float _S4385 = scale_28.y;
    float sy_7 = (F32_exp((_S4385)));
    float sz_9 = scale_28.z - 0.5f * (_S4384 + _S4385);
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
    Matrix<float, 3, 3>  _S4386 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_52), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_52), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4386, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4386, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4386, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_8 = mul_0(R_29, vert0_5) + t_26;
    float3  vert1_c_8 = mul_0(R_29, vert1_5) + t_26;
    float3  vert2_c_8 = mul_0(R_29, vert2_5) + t_26;
    float2  _S4387 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z);
    float u_68 = _S4387.x;
    float v_68 = _S4387.y;
    float r2_68 = u_68 * u_68 + v_68 * v_68;
    float _S4388 = 2.0f * (*dist_coeffs_36)[int(4)];
    float _S4389 = 2.0f * (*dist_coeffs_36)[int(5)];
    float2  _S4390 = _S4387 * make_float2 (1.0f + r2_68 * ((*dist_coeffs_36)[int(0)] + r2_68 * ((*dist_coeffs_36)[int(1)] + r2_68 * ((*dist_coeffs_36)[int(2)] + r2_68 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4388 * u_68 * v_68 + (*dist_coeffs_36)[int(5)] * (r2_68 + 2.0f * u_68 * u_68) + (*dist_coeffs_36)[int(6)] * r2_68, _S4389 * u_68 * v_68 + (*dist_coeffs_36)[int(4)] * (r2_68 + 2.0f * v_68 * v_68) + (*dist_coeffs_36)[int(7)] * r2_68);
    float2  _S4391 = _S4390 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4390.x + (*dist_coeffs_36)[int(9)] * _S4390.y, 0.0f);
    float _S4392 = fx_32 * _S4391.x + cx_27;
    float _S4393 = fy_32 * _S4391.y + cy_27;
    float2  uv0_8 = make_float2 (_S4392, _S4393);
    float2  _S4394 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z);
    float u_69 = _S4394.x;
    float v_69 = _S4394.y;
    float r2_69 = u_69 * u_69 + v_69 * v_69;
    float2  _S4395 = _S4394 * make_float2 (1.0f + r2_69 * ((*dist_coeffs_36)[int(0)] + r2_69 * ((*dist_coeffs_36)[int(1)] + r2_69 * ((*dist_coeffs_36)[int(2)] + r2_69 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4388 * u_69 * v_69 + (*dist_coeffs_36)[int(5)] * (r2_69 + 2.0f * u_69 * u_69) + (*dist_coeffs_36)[int(6)] * r2_69, _S4389 * u_69 * v_69 + (*dist_coeffs_36)[int(4)] * (r2_69 + 2.0f * v_69 * v_69) + (*dist_coeffs_36)[int(7)] * r2_69);
    float2  _S4396 = _S4395 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4395.x + (*dist_coeffs_36)[int(9)] * _S4395.y, 0.0f);
    float _S4397 = fx_32 * _S4396.x + cx_27;
    float _S4398 = fy_32 * _S4396.y + cy_27;
    float2  uv1_8 = make_float2 (_S4397, _S4398);
    float2  _S4399 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z);
    float u_70 = _S4399.x;
    float v_70 = _S4399.y;
    float r2_70 = u_70 * u_70 + v_70 * v_70;
    float2  _S4400 = _S4399 * make_float2 (1.0f + r2_70 * ((*dist_coeffs_36)[int(0)] + r2_70 * ((*dist_coeffs_36)[int(1)] + r2_70 * ((*dist_coeffs_36)[int(2)] + r2_70 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4388 * u_70 * v_70 + (*dist_coeffs_36)[int(5)] * (r2_70 + 2.0f * u_70 * u_70) + (*dist_coeffs_36)[int(6)] * r2_70, _S4389 * u_70 * v_70 + (*dist_coeffs_36)[int(4)] * (r2_70 + 2.0f * v_70 * v_70) + (*dist_coeffs_36)[int(7)] * r2_70);
    float2  _S4401 = _S4400 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4400.x + (*dist_coeffs_36)[int(9)] * _S4400.y, 0.0f);
    float _S4402 = fx_32 * _S4401.x + cx_27;
    float _S4403 = fy_32 * _S4401.y + cy_27;
    float2  uv2_8 = make_float2 (_S4402, _S4403);
    float2  e0_12 = uv1_8 - uv0_8;
    float2  e1_12 = uv2_8 - uv1_8;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(uv0_8 - uv2_8)));
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4392), (_S4397)))), (_S4402))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4393), (_S4398)))), (_S4403))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4392), (_S4397)))), (_S4402))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4393), (_S4398)))), (_S4403))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4404 = mean_28 - - mul_0(transpose_0(R_29), t_26);
    float _S4405 = _S4404.x;
    float _S4406 = _S4404.y;
    float _S4407 = _S4404.z;
    float norm_16 = (F32_sqrt((_S4405 * _S4405 + _S4406 * _S4406 + _S4407 * _S4407)));
    float x_58 = _S4405 / norm_16;
    float y_26 = _S4406 / norm_16;
    float z_23 = _S4407 / norm_16;
    float z2_53 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_58 * x_58 - y_26 * y_26;
    float fS1_23 = 2.0f * x_58 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_53 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_58) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_53 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_58) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_58 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_53 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_58) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_58 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S4408 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4408);
    float3  _S4409 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S4410 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S4409 + _S4410 + make_float3 (0.5f), _S4408);
    (*rgbs_2)[int(2)] = max_0(_S4409 - _S4410 + make_float3 (0.5f), _S4408);
    (*verts_2)[int(0)] = vert0_5;
    (*verts_2)[int(1)] = vert1_5;
    (*verts_2)[int(2)] = vert2_5;
    float3  _S4411 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S4411 * make_float3 (float(- (F32_sign((dot_0(_S4411, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_30, float3  t_27, float fx_33, float fy_33, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_37, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_30, mean_29) + t_27;
    float _S4412 = scale_29.x;
    float sx_8 = (F32_exp((_S4412)));
    float _S4413 = scale_29.y;
    float sy_8 = (F32_exp((_S4413)));
    float sz_10 = scale_29.z - 0.5f * (_S4412 + _S4413);
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
    Matrix<float, 3, 3>  _S4414 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_54), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_54), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  vert0_6 = mul_0(_S4414, make_float3 (sx_8, 0.0f, 0.0f)) + mean_29;
    float3  vert1_6 = mul_0(_S4414, make_float3 (sx_8 * (-0.5f + sz_10), sy_8, 0.0f)) + mean_29;
    float3  vert2_6 = mul_0(_S4414, make_float3 (sx_8 * (-0.5f - sz_10), - sy_8, 0.0f)) + mean_29;
    float3  vert0_c_9 = mul_0(R_30, vert0_6) + t_27;
    float3  vert1_c_9 = mul_0(R_30, vert1_6) + t_27;
    float3  vert2_c_9 = mul_0(R_30, vert2_6) + t_27;
    float2  _S4415 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_33 = length_0(_S4415);
    float _S4416 = vert0_c_9.z;
    float theta_29 = (F32_atan2((r_33), (_S4416)));
    float k_11;
    if(theta_29 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_29 * theta_29 / 3.0f) / _S4416;
    }
    else
    {
        k_11 = theta_29 / r_33;
    }
    float2  _S4417 = _S4415 * make_float2 (k_11);
    float u_71 = _S4417.x;
    float v_71 = _S4417.y;
    float r2_71 = u_71 * u_71 + v_71 * v_71;
    float _S4418 = 2.0f * (*dist_coeffs_37)[int(4)];
    float _S4419 = 2.0f * (*dist_coeffs_37)[int(5)];
    float2  _S4420 = _S4417 * make_float2 (1.0f + r2_71 * ((*dist_coeffs_37)[int(0)] + r2_71 * ((*dist_coeffs_37)[int(1)] + r2_71 * ((*dist_coeffs_37)[int(2)] + r2_71 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4418 * u_71 * v_71 + (*dist_coeffs_37)[int(5)] * (r2_71 + 2.0f * u_71 * u_71) + (*dist_coeffs_37)[int(6)] * r2_71, _S4419 * u_71 * v_71 + (*dist_coeffs_37)[int(4)] * (r2_71 + 2.0f * v_71 * v_71) + (*dist_coeffs_37)[int(7)] * r2_71);
    float2  _S4421 = _S4420 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4420.x + (*dist_coeffs_37)[int(9)] * _S4420.y, 0.0f);
    float _S4422 = fx_33 * _S4421.x + cx_28;
    float _S4423 = fy_33 * _S4421.y + cy_28;
    float2  uv0_9 = make_float2 (_S4422, _S4423);
    float2  _S4424 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_34 = length_0(_S4424);
    float _S4425 = vert1_c_9.z;
    float theta_30 = (F32_atan2((r_34), (_S4425)));
    if(theta_30 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_30 * theta_30 / 3.0f) / _S4425;
    }
    else
    {
        k_11 = theta_30 / r_34;
    }
    float2  _S4426 = _S4424 * make_float2 (k_11);
    float u_72 = _S4426.x;
    float v_72 = _S4426.y;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float2  _S4427 = _S4426 * make_float2 (1.0f + r2_72 * ((*dist_coeffs_37)[int(0)] + r2_72 * ((*dist_coeffs_37)[int(1)] + r2_72 * ((*dist_coeffs_37)[int(2)] + r2_72 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4418 * u_72 * v_72 + (*dist_coeffs_37)[int(5)] * (r2_72 + 2.0f * u_72 * u_72) + (*dist_coeffs_37)[int(6)] * r2_72, _S4419 * u_72 * v_72 + (*dist_coeffs_37)[int(4)] * (r2_72 + 2.0f * v_72 * v_72) + (*dist_coeffs_37)[int(7)] * r2_72);
    float2  _S4428 = _S4427 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4427.x + (*dist_coeffs_37)[int(9)] * _S4427.y, 0.0f);
    float _S4429 = fx_33 * _S4428.x + cx_28;
    float _S4430 = fy_33 * _S4428.y + cy_28;
    float2  uv1_9 = make_float2 (_S4429, _S4430);
    float2  _S4431 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_35 = length_0(_S4431);
    float _S4432 = vert2_c_9.z;
    float theta_31 = (F32_atan2((r_35), (_S4432)));
    if(theta_31 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_31 * theta_31 / 3.0f) / _S4432;
    }
    else
    {
        k_11 = theta_31 / r_35;
    }
    float2  _S4433 = _S4431 * make_float2 (k_11);
    float u_73 = _S4433.x;
    float v_73 = _S4433.y;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float2  _S4434 = _S4433 * make_float2 (1.0f + r2_73 * ((*dist_coeffs_37)[int(0)] + r2_73 * ((*dist_coeffs_37)[int(1)] + r2_73 * ((*dist_coeffs_37)[int(2)] + r2_73 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4418 * u_73 * v_73 + (*dist_coeffs_37)[int(5)] * (r2_73 + 2.0f * u_73 * u_73) + (*dist_coeffs_37)[int(6)] * r2_73, _S4419 * u_73 * v_73 + (*dist_coeffs_37)[int(4)] * (r2_73 + 2.0f * v_73 * v_73) + (*dist_coeffs_37)[int(7)] * r2_73);
    float2  _S4435 = _S4434 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4434.x + (*dist_coeffs_37)[int(9)] * _S4434.y, 0.0f);
    float _S4436 = fx_33 * _S4435.x + cx_28;
    float _S4437 = fy_33 * _S4435.y + cy_28;
    float2  uv2_9 = make_float2 (_S4436, _S4437);
    float2  e0_13 = uv1_9 - uv0_9;
    float2  e1_13 = uv2_9 - uv1_9;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(uv0_9 - uv2_9)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4422), (_S4429)))), (_S4436))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4423), (_S4430)))), (_S4437))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4422), (_S4429)))), (_S4436))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4423), (_S4430)))), (_S4437))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4438 = mean_29 - - mul_0(transpose_0(R_30), t_27);
    float _S4439 = _S4438.x;
    float _S4440 = _S4438.y;
    float _S4441 = _S4438.z;
    float norm_17 = (F32_sqrt((_S4439 * _S4439 + _S4440 * _S4440 + _S4441 * _S4441)));
    float x_60 = _S4439 / norm_17;
    float y_27 = _S4440 / norm_17;
    float z_24 = _S4441 / norm_17;
    float z2_55 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_60 * x_60 - y_27 * y_27;
    float fS1_24 = 2.0f * x_60 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_55 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_60) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_55 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_60) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_60 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_55 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_60) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_60 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4442 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4442);
    float3  _S4443 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4444 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4443 + _S4444 + make_float3 (0.5f), _S4442);
    (*rgbs_3)[int(2)] = max_0(_S4443 - _S4444 + make_float3 (0.5f), _S4442);
    (*verts_3)[int(0)] = vert0_6;
    (*verts_3)[int(1)] = vert1_6;
    (*verts_3)[int(2)] = vert2_6;
    float3  _S4445 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4445 * make_float3 (float(- (F32_sign((dot_0(_S4445, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_31, float3  t_28, float fx_34, float fy_34, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_38, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_0(R_31, mean_30) + t_28;
    float _S4446 = scale_30.x;
    float _S4447 = s_primal_ctx_exp_1(_S4446);
    float _S4448 = scale_30.y;
    float _S4449 = s_primal_ctx_exp_1(_S4448);
    float sz_11 = scale_30.z - 0.5f * (_S4446 + _S4448);
    float _S4450 = quat_31.y;
    float x2_31 = _S4450 * _S4450;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_56 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4451 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_56), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_56), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4452 = make_float3 (_S4447, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_0(_S4451, _S4452) + mean_30;
    float _S4453 = -0.5f + sz_11;
    float3  _S4454 = make_float3 (_S4447 * _S4453, _S4449, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_0(_S4451, _S4454) + mean_30;
    float _S4455 = -0.5f - sz_11;
    float3  _S4456 = make_float3 (_S4447 * _S4455, - _S4449, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_0(_S4451, _S4456) + mean_30;
    float3  vert0_c_10 = s_primal_ctx_mul_0(R_31, vert0_7) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_0(R_31, vert1_7) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_0(R_31, vert2_7) + t_28;
    float2  _S4457 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4458 = vert0_c_10.z;
    float2  _S4459 = make_float2 (_S4458);
    float2  _S4460 = _S4457 / make_float2 (_S4458);
    float2  _S4461 = make_float2 (_S4458 * _S4458);
    float u_74 = _S4460.x;
    float v_74 = _S4460.y;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float _S4462 = (*dist_coeffs_38)[int(2)] + r2_74 * (*dist_coeffs_38)[int(3)];
    float _S4463 = (*dist_coeffs_38)[int(1)] + r2_74 * _S4462;
    float _S4464 = (*dist_coeffs_38)[int(0)] + r2_74 * _S4463;
    float radial_7 = 1.0f + r2_74 * _S4464;
    float2  _S4465 = make_float2 (radial_7);
    float _S4466 = 2.0f * (*dist_coeffs_38)[int(4)];
    float _S4467 = _S4466 * u_74;
    float _S4468 = 2.0f * u_74;
    float _S4469 = 2.0f * (*dist_coeffs_38)[int(5)];
    float _S4470 = _S4469 * u_74;
    float _S4471 = 2.0f * v_74;
    float2  _S4472 = _S4460 * make_float2 (radial_7) + make_float2 (_S4467 * v_74 + (*dist_coeffs_38)[int(5)] * (r2_74 + _S4468 * u_74) + (*dist_coeffs_38)[int(6)] * r2_74, _S4470 * v_74 + (*dist_coeffs_38)[int(4)] * (r2_74 + _S4471 * v_74) + (*dist_coeffs_38)[int(7)] * r2_74);
    float2  _S4473 = _S4472 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4472.x + (*dist_coeffs_38)[int(9)] * _S4472.y, 0.0f);
    float _S4474 = fx_34 * _S4473.x + cx_29;
    float _S4475 = fy_34 * _S4473.y + cy_29;
    float2  uv0_10 = make_float2 (_S4474, _S4475);
    float2  _S4476 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4477 = vert1_c_10.z;
    float2  _S4478 = make_float2 (_S4477);
    float2  _S4479 = _S4476 / make_float2 (_S4477);
    float2  _S4480 = make_float2 (_S4477 * _S4477);
    float u_75 = _S4479.x;
    float v_75 = _S4479.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S4481 = (*dist_coeffs_38)[int(2)] + r2_75 * (*dist_coeffs_38)[int(3)];
    float _S4482 = (*dist_coeffs_38)[int(1)] + r2_75 * _S4481;
    float _S4483 = (*dist_coeffs_38)[int(0)] + r2_75 * _S4482;
    float radial_8 = 1.0f + r2_75 * _S4483;
    float2  _S4484 = make_float2 (radial_8);
    float _S4485 = _S4466 * u_75;
    float _S4486 = 2.0f * u_75;
    float _S4487 = _S4469 * u_75;
    float _S4488 = 2.0f * v_75;
    float2  _S4489 = _S4479 * make_float2 (radial_8) + make_float2 (_S4485 * v_75 + (*dist_coeffs_38)[int(5)] * (r2_75 + _S4486 * u_75) + (*dist_coeffs_38)[int(6)] * r2_75, _S4487 * v_75 + (*dist_coeffs_38)[int(4)] * (r2_75 + _S4488 * v_75) + (*dist_coeffs_38)[int(7)] * r2_75);
    float2  _S4490 = _S4489 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4489.x + (*dist_coeffs_38)[int(9)] * _S4489.y, 0.0f);
    float _S4491 = fx_34 * _S4490.x + cx_29;
    float _S4492 = fy_34 * _S4490.y + cy_29;
    float2  uv1_10 = make_float2 (_S4491, _S4492);
    float2  _S4493 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4494 = vert2_c_10.z;
    float2  _S4495 = make_float2 (_S4494);
    float2  _S4496 = _S4493 / make_float2 (_S4494);
    float2  _S4497 = make_float2 (_S4494 * _S4494);
    float u_76 = _S4496.x;
    float v_76 = _S4496.y;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float _S4498 = (*dist_coeffs_38)[int(2)] + r2_76 * (*dist_coeffs_38)[int(3)];
    float _S4499 = (*dist_coeffs_38)[int(1)] + r2_76 * _S4498;
    float _S4500 = (*dist_coeffs_38)[int(0)] + r2_76 * _S4499;
    float radial_9 = 1.0f + r2_76 * _S4500;
    float2  _S4501 = make_float2 (radial_9);
    float _S4502 = _S4466 * u_76;
    float _S4503 = 2.0f * u_76;
    float _S4504 = _S4469 * u_76;
    float _S4505 = 2.0f * v_76;
    float2  _S4506 = _S4496 * make_float2 (radial_9) + make_float2 (_S4502 * v_76 + (*dist_coeffs_38)[int(5)] * (r2_76 + _S4503 * u_76) + (*dist_coeffs_38)[int(6)] * r2_76, _S4504 * v_76 + (*dist_coeffs_38)[int(4)] * (r2_76 + _S4505 * v_76) + (*dist_coeffs_38)[int(7)] * r2_76);
    float2  _S4507 = _S4506 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4506.x + (*dist_coeffs_38)[int(9)] * _S4506.y, 0.0f);
    float _S4508 = fx_34 * _S4507.x + cx_29;
    float _S4509 = fy_34 * _S4507.y + cy_29;
    float2  uv2_10 = make_float2 (_S4508, _S4509);
    float2  e0_14 = uv1_10 - uv0_10;
    float2  e1_14 = uv2_10 - uv1_10;
    float2  e2_6 = uv0_10 - uv2_10;
    float _S4510 = e0_14.x;
    float _S4511 = e1_14.y;
    float _S4512 = e0_14.y;
    float _S4513 = e1_14.x;
    float _S4514 = _S4510 * _S4511 - _S4512 * _S4513;
    float _S4515 = 1.0f - hardness_14.y;
    float _S4516 = -1.0f / _S4515;
    float _S4517 = _S4515 * _S4515;
    float _S4518 = s_primal_ctx_max_0(_S4474, _S4491);
    float _S4519 = s_primal_ctx_min_0(_S4474, _S4491);
    float _S4520 = s_primal_ctx_max_0(_S4475, _S4492);
    float _S4521 = s_primal_ctx_min_0(_S4475, _S4492);
    float3  _S4522 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4523 = length_1(_S4522) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4524 = transpose_0(R_31);
    float3  _S4525 = mean_30 - - s_primal_ctx_mul_0(_S4524, t_28);
    float _S4526 = _S4525.x;
    float _S4527 = _S4525.y;
    float _S4528 = _S4525.z;
    float _S4529 = _S4526 * _S4526 + _S4527 * _S4527 + _S4528 * _S4528;
    float _S4530 = s_primal_ctx_sqrt_0(_S4529);
    float x_61 = _S4526 / _S4530;
    float3  _S4531 = make_float3 (x_61);
    float _S4532 = _S4530 * _S4530;
    float y_28 = _S4527 / _S4530;
    float z_25 = _S4528 / _S4530;
    float3  _S4533 = make_float3 (z_25);
    float _S4534 = - y_28;
    float3  _S4535 = make_float3 (_S4534);
    float z2_57 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_61 * x_61 - y_28 * y_28;
    float _S4536 = 2.0f * x_61;
    float fS1_25 = _S4536 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_57 - 0.31539157032966614f;
    float3  _S4537 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_61;
    float3  _S4538 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4539 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4540 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4541 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_57 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4542 = 1.86588168144226074f * z2_57 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4542;
    float3  _S4543 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_61;
    float3  _S4544 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4545 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4546 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4547 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_61 * fC1_25 - y_28 * fS1_25);
    float3  _S4548 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_61 * fS1_25 + y_28 * fC1_25);
    float3  _S4549 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4534) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_61) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4550 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4551 = make_float3 (0.0f);
    float3  _S4552 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4553 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4554 = make_float3 (_S4553);
    float3  _S4555 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4553);
    float3  _S4556 = _S4552 + _S4555 + make_float3 (0.5f);
    float3  _S4557 = _S4552 - _S4555 + make_float3 (0.5f);
    float3  _S4558 = vert1_c_10 - vert0_c_10;
    float3  _S4559 = vert2_c_10 - vert0_c_10;
    float3  _S4560 = s_primal_ctx_cross_0(_S4558, _S4559);
    float3  _S4561 = normalize_0(_S4560);
    float3  _S4562 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4561, mean_c_25)))))) * v_normal_2;
    float3  _S4563 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4564;
    (&_S4564)->primal_0 = _S4561;
    (&_S4564)->differential_0 = _S4563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4565;
    (&_S4565)->primal_0 = mean_c_25;
    (&_S4565)->differential_0 = _S4563;
    s_bwd_prop_dot_0(&_S4564, &_S4565, 0.0f);
    float3  _S4566 = _S4562 + _S4564.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4567;
    (&_S4567)->primal_0 = _S4560;
    (&_S4567)->differential_0 = _S4563;
    s_bwd_normalize_impl_0(&_S4567, _S4566);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4568;
    (&_S4568)->primal_0 = _S4558;
    (&_S4568)->differential_0 = _S4563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4569;
    (&_S4569)->primal_0 = _S4559;
    (&_S4569)->differential_0 = _S4563;
    s_bwd_prop_cross_0(&_S4568, &_S4569, _S4567.differential_0);
    float3  _S4570 = - _S4569.differential_0;
    float3  _S4571 = - _S4568.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4572;
    (&_S4572)->primal_0 = _S4557;
    (&_S4572)->differential_0 = _S4563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4573;
    (&_S4573)->primal_0 = _S4551;
    (&_S4573)->differential_0 = _S4563;
    s_bwd_prop_max_0(&_S4572, &_S4573, (*v_rgbs_0)[int(2)]);
    float3  _S4574 = - _S4572.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4575;
    (&_S4575)->primal_0 = _S4556;
    (&_S4575)->differential_0 = _S4563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4576;
    (&_S4576)->primal_0 = _S4551;
    (&_S4576)->differential_0 = _S4563;
    s_bwd_prop_max_0(&_S4575, &_S4576, (*v_rgbs_0)[int(1)]);
    float3  _S4577 = _S4554 * (_S4574 + _S4575.differential_0);
    float3  _S4578 = _S4572.differential_0 + _S4575.differential_0;
    float3  _S4579 = make_float3 (0.5f) * - _S4578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4580;
    (&_S4580)->primal_0 = _S4550;
    (&_S4580)->differential_0 = _S4563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4581;
    (&_S4581)->primal_0 = _S4551;
    (&_S4581)->differential_0 = _S4563;
    s_bwd_prop_max_0(&_S4580, &_S4581, (*v_rgbs_0)[int(0)]);
    float3  _S4582 = _S4579 + _S4580.differential_0;
    float3  _S4583 = _S4578 + _S4580.differential_0;
    float3  _S4584 = _S4548 * _S4583;
    float3  _S4585 = (*sh_coeffs_25)[int(15)] * _S4583;
    float3  _S4586 = _S4546 * _S4583;
    float3  _S4587 = (*sh_coeffs_25)[int(14)] * _S4583;
    float3  _S4588 = _S4544 * _S4583;
    float3  _S4589 = (*sh_coeffs_25)[int(13)] * _S4583;
    float3  _S4590 = _S4543 * _S4583;
    float3  _S4591 = (*sh_coeffs_25)[int(12)] * _S4583;
    float3  _S4592 = _S4545 * _S4583;
    float3  _S4593 = (*sh_coeffs_25)[int(11)] * _S4583;
    float3  _S4594 = _S4547 * _S4583;
    float3  _S4595 = (*sh_coeffs_25)[int(10)] * _S4583;
    float3  _S4596 = _S4549 * _S4583;
    float3  _S4597 = (*sh_coeffs_25)[int(9)] * _S4583;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4597.x + _S4597.y + _S4597.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4585.x + _S4585.y + _S4585.z);
    float _S4598 = _S4595.x + _S4595.y + _S4595.z;
    float _S4599 = _S4587.x + _S4587.y + _S4587.z;
    float _S4600 = _S4593.x + _S4593.y + _S4593.z;
    float _S4601 = _S4589.x + _S4589.y + _S4589.z;
    float _S4602 = _S4591.x + _S4591.y + _S4591.z;
    float _S4603 = - s_diff_fC2_T_7;
    float3  _S4604 = _S4540 * _S4583;
    float3  _S4605 = (*sh_coeffs_25)[int(8)] * _S4583;
    float3  _S4606 = _S4538 * _S4583;
    float3  _S4607 = (*sh_coeffs_25)[int(7)] * _S4583;
    float3  _S4608 = _S4537 * _S4583;
    float3  _S4609 = (*sh_coeffs_25)[int(6)] * _S4583;
    float3  _S4610 = _S4539 * _S4583;
    float3  _S4611 = (*sh_coeffs_25)[int(5)] * _S4583;
    float3  _S4612 = _S4541 * _S4583;
    float3  _S4613 = (*sh_coeffs_25)[int(4)] * _S4583;
    float _S4614 = _S4611.x + _S4611.y + _S4611.z;
    float _S4615 = _S4607.x + _S4607.y + _S4607.z;
    float _S4616 = fTmp1B_25 * _S4598 + x_61 * s_diff_fS2_T_7 + y_28 * _S4603 + 0.54627424478530884f * (_S4613.x + _S4613.y + _S4613.z);
    float _S4617 = fTmp1B_25 * _S4599 + y_28 * s_diff_fS2_T_7 + x_61 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4605.x + _S4605.y + _S4605.z);
    float _S4618 = y_28 * - _S4617;
    float _S4619 = x_61 * _S4617;
    float _S4620 = z_25 * (1.86588168144226074f * (z_25 * _S4602) + -2.28522896766662598f * (y_28 * _S4600 + x_61 * _S4601) + 0.94617468118667603f * (_S4609.x + _S4609.y + _S4609.z));
    float3  _S4621 = make_float3 (0.48860251903533936f) * _S4583;
    float3  _S4622 = - _S4621;
    float3  _S4623 = _S4531 * _S4622;
    float3  _S4624 = (*sh_coeffs_25)[int(3)] * _S4622;
    float3  _S4625 = _S4533 * _S4621;
    float3  _S4626 = (*sh_coeffs_25)[int(2)] * _S4621;
    float3  _S4627 = _S4535 * _S4621;
    float3  _S4628 = (*sh_coeffs_25)[int(1)] * _S4621;
    float _S4629 = (_S4542 * _S4602 + 1.44530570507049561f * (fS1_25 * _S4598 + fC1_25 * _S4599) + -1.09254848957061768f * (y_28 * _S4614 + x_61 * _S4615) + _S4620 + _S4620 + _S4626.x + _S4626.y + _S4626.z) / _S4532;
    float _S4630 = _S4530 * _S4629;
    float _S4631 = (fTmp0C_25 * _S4600 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4603 + fTmp0B_25 * _S4614 + _S4536 * _S4616 + _S4618 + _S4618 + - (_S4628.x + _S4628.y + _S4628.z)) / _S4532;
    float _S4632 = _S4530 * _S4631;
    float _S4633 = (fTmp0C_25 * _S4601 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4615 + 2.0f * (y_28 * _S4616) + _S4619 + _S4619 + _S4624.x + _S4624.y + _S4624.z) / _S4532;
    float _S4634 = _S4530 * _S4633;
    float _S4635 = _S4528 * - _S4629 + _S4527 * - _S4631 + _S4526 * - _S4633;
    DiffPair_float_0 _S4636;
    (&_S4636)->primal_0 = _S4529;
    (&_S4636)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4636, _S4635);
    float _S4637 = _S4528 * _S4636.differential_0;
    float _S4638 = _S4527 * _S4636.differential_0;
    float _S4639 = _S4526 * _S4636.differential_0;
    float3  _S4640 = make_float3 (0.282094806432724f) * _S4583;
    float3  _S4641 = make_float3 (_S4634 + _S4639 + _S4639, _S4632 + _S4638 + _S4638, _S4630 + _S4637 + _S4637);
    float3  _S4642 = - - _S4641;
    Matrix<float, 3, 3>  _S4643 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4644;
    (&_S4644)->primal_0 = _S4524;
    (&_S4644)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4645;
    (&_S4645)->primal_0 = t_28;
    (&_S4645)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4644, &_S4645, _S4642);
    Matrix<float, 3, 3>  _S4646 = transpose_0(_S4644.differential_0);
    DiffPair_float_0 _S4647;
    (&_S4647)->primal_0 = _S4523;
    (&_S4647)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4647, v_depth_9);
    float _S4648 = 0.3333333432674408f * _S4647.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4649;
    (&_S4649)->primal_0 = _S4522;
    (&_S4649)->differential_0 = _S4563;
    s_bwd_length_impl_1(&_S4649, _S4648);
    DiffPair_float_0 _S4650;
    (&_S4650)->primal_0 = _S4521;
    (&_S4650)->differential_0 = 0.0f;
    DiffPair_float_0 _S4651;
    (&_S4651)->primal_0 = _S4509;
    (&_S4651)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4650, &_S4651, 0.0f);
    DiffPair_float_0 _S4652;
    (&_S4652)->primal_0 = _S4475;
    (&_S4652)->differential_0 = 0.0f;
    DiffPair_float_0 _S4653;
    (&_S4653)->primal_0 = _S4492;
    (&_S4653)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4652, &_S4653, _S4650.differential_0);
    DiffPair_float_0 _S4654;
    (&_S4654)->primal_0 = _S4520;
    (&_S4654)->differential_0 = 0.0f;
    DiffPair_float_0 _S4655;
    (&_S4655)->primal_0 = _S4509;
    (&_S4655)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4654, &_S4655, 0.0f);
    DiffPair_float_0 _S4656;
    (&_S4656)->primal_0 = _S4475;
    (&_S4656)->differential_0 = 0.0f;
    DiffPair_float_0 _S4657;
    (&_S4657)->primal_0 = _S4492;
    (&_S4657)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4656, &_S4657, _S4654.differential_0);
    DiffPair_float_0 _S4658;
    (&_S4658)->primal_0 = _S4519;
    (&_S4658)->differential_0 = 0.0f;
    DiffPair_float_0 _S4659;
    (&_S4659)->primal_0 = _S4508;
    (&_S4659)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4658, &_S4659, 0.0f);
    DiffPair_float_0 _S4660;
    (&_S4660)->primal_0 = _S4474;
    (&_S4660)->differential_0 = 0.0f;
    DiffPair_float_0 _S4661;
    (&_S4661)->primal_0 = _S4491;
    (&_S4661)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4660, &_S4661, _S4658.differential_0);
    DiffPair_float_0 _S4662;
    (&_S4662)->primal_0 = _S4518;
    (&_S4662)->differential_0 = 0.0f;
    DiffPair_float_0 _S4663;
    (&_S4663)->primal_0 = _S4508;
    (&_S4663)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4662, &_S4663, 0.0f);
    DiffPair_float_0 _S4664;
    (&_S4664)->primal_0 = _S4474;
    (&_S4664)->differential_0 = 0.0f;
    DiffPair_float_0 _S4665;
    (&_S4665)->primal_0 = _S4491;
    (&_S4665)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4664, &_S4665, _S4662.differential_0);
    DiffPair_float_0 _S4666;
    (&_S4666)->primal_0 = _S4516;
    (&_S4666)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4666, 0.0f);
    float _S4667 = - (-1.0f * - (_S4666.differential_0 / _S4517));
    float2  _S4668 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4669;
    (&_S4669)->primal_0 = e2_6;
    (&_S4669)->differential_0 = _S4668;
    s_bwd_length_impl_0(&_S4669, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4670;
    (&_S4670)->primal_0 = e1_14;
    (&_S4670)->differential_0 = _S4668;
    s_bwd_length_impl_0(&_S4670, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4671;
    (&_S4671)->primal_0 = e0_14;
    (&_S4671)->differential_0 = _S4668;
    s_bwd_length_impl_0(&_S4671, -0.0f);
    DiffPair_float_0 _S4672;
    (&_S4672)->primal_0 = _S4514;
    (&_S4672)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4672, 0.0f);
    float _S4673 = - _S4672.differential_0;
    float2  _S4674 = _S4670.differential_0 + make_float2 (_S4512 * _S4673, _S4510 * _S4672.differential_0);
    float2  _S4675 = _S4671.differential_0 + make_float2 (_S4511 * _S4672.differential_0, _S4513 * _S4673);
    float2  _S4676 = - _S4669.differential_0 + _S4674;
    float _S4677 = fx_34 * (_S4659.differential_0 + _S4663.differential_0 + _S4676.x);
    float2  _S4678 = make_float2 (_S4677, fy_34 * (_S4651.differential_0 + _S4655.differential_0 + _S4676.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4677, (*dist_coeffs_38)[int(9)] * _S4677);
    float2  _S4679 = _S4496 * _S4678;
    float _S4680 = (*dist_coeffs_38)[int(4)] * _S4678.y;
    float _S4681 = (*dist_coeffs_38)[int(5)] * _S4678.x;
    float _S4682 = _S4679.x + _S4679.y;
    float _S4683 = r2_76 * _S4682;
    float _S4684 = r2_76 * _S4683;
    float _S4685 = (*dist_coeffs_38)[int(7)] * _S4678.y + _S4680 + (*dist_coeffs_38)[int(6)] * _S4678.x + _S4681 + _S4500 * _S4682 + _S4499 * _S4683 + _S4498 * _S4684 + (*dist_coeffs_38)[int(3)] * (r2_76 * _S4684);
    float _S4686 = v_76 * _S4685;
    float _S4687 = u_76 * _S4685;
    float2  _S4688 = (_S4501 * _S4678 + make_float2 (_S4469 * (v_76 * _S4678.y) + _S4503 * _S4681 + 2.0f * (u_76 * _S4681) + _S4466 * (v_76 * _S4678.x) + _S4687 + _S4687, _S4505 * _S4680 + 2.0f * (v_76 * _S4680) + _S4504 * _S4678.y + _S4502 * _S4678.x + _S4686 + _S4686)) / _S4497;
    float2  _S4689 = _S4493 * - _S4688;
    float2  _S4690 = _S4495 * _S4688;
    float2  _S4691 = - _S4674 + _S4675;
    float _S4692 = fx_34 * (_S4661.differential_0 + _S4665.differential_0 + _S4691.x);
    float2  _S4693 = make_float2 (_S4692, fy_34 * (_S4653.differential_0 + _S4657.differential_0 + _S4691.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4692, (*dist_coeffs_38)[int(9)] * _S4692);
    float2  _S4694 = _S4479 * _S4693;
    float _S4695 = (*dist_coeffs_38)[int(4)] * _S4693.y;
    float _S4696 = (*dist_coeffs_38)[int(5)] * _S4693.x;
    float _S4697 = _S4694.x + _S4694.y;
    float _S4698 = r2_75 * _S4697;
    float _S4699 = r2_75 * _S4698;
    float _S4700 = (*dist_coeffs_38)[int(7)] * _S4693.y + _S4695 + (*dist_coeffs_38)[int(6)] * _S4693.x + _S4696 + _S4483 * _S4697 + _S4482 * _S4698 + _S4481 * _S4699 + (*dist_coeffs_38)[int(3)] * (r2_75 * _S4699);
    float _S4701 = v_75 * _S4700;
    float _S4702 = u_75 * _S4700;
    float2  _S4703 = (_S4484 * _S4693 + make_float2 (_S4469 * (v_75 * _S4693.y) + _S4486 * _S4696 + 2.0f * (u_75 * _S4696) + _S4466 * (v_75 * _S4693.x) + _S4702 + _S4702, _S4488 * _S4695 + 2.0f * (v_75 * _S4695) + _S4487 * _S4693.y + _S4485 * _S4693.x + _S4701 + _S4701)) / _S4480;
    float2  _S4704 = _S4476 * - _S4703;
    float2  _S4705 = _S4478 * _S4703;
    float _S4706 = _S4704.x + _S4704.y;
    float2  _S4707 = _S4669.differential_0 + - _S4675;
    float _S4708 = fx_34 * (_S4660.differential_0 + _S4664.differential_0 + _S4707.x);
    float2  _S4709 = make_float2 (_S4708, fy_34 * (_S4652.differential_0 + _S4656.differential_0 + _S4707.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4708, (*dist_coeffs_38)[int(9)] * _S4708);
    float2  _S4710 = _S4460 * _S4709;
    float _S4711 = (*dist_coeffs_38)[int(4)] * _S4709.y;
    float _S4712 = (*dist_coeffs_38)[int(5)] * _S4709.x;
    float _S4713 = _S4710.x + _S4710.y;
    float _S4714 = r2_74 * _S4713;
    float _S4715 = r2_74 * _S4714;
    float _S4716 = (*dist_coeffs_38)[int(7)] * _S4709.y + _S4711 + (*dist_coeffs_38)[int(6)] * _S4709.x + _S4712 + _S4464 * _S4713 + _S4463 * _S4714 + _S4462 * _S4715 + (*dist_coeffs_38)[int(3)] * (r2_74 * _S4715);
    float _S4717 = v_74 * _S4716;
    float _S4718 = u_74 * _S4716;
    float2  _S4719 = (_S4465 * _S4709 + make_float2 (_S4469 * (v_74 * _S4709.y) + _S4468 * _S4712 + 2.0f * (u_74 * _S4712) + _S4466 * (v_74 * _S4709.x) + _S4718 + _S4718, _S4471 * _S4711 + 2.0f * (v_74 * _S4711) + _S4470 * _S4709.y + _S4467 * _S4709.x + _S4717 + _S4717)) / _S4461;
    float2  _S4720 = _S4457 * - _S4719;
    float2  _S4721 = _S4459 * _S4719;
    float _S4722 = _S4720.x + _S4720.y;
    float3  _S4723 = _S4569.differential_0 + _S4649.differential_0 + make_float3 (_S4690.x, _S4690.y, _S4689.x + _S4689.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4724;
    (&_S4724)->primal_0 = R_31;
    (&_S4724)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4725;
    (&_S4725)->primal_0 = vert2_7;
    (&_S4725)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4724, &_S4725, _S4723);
    float3  _S4726 = _S4568.differential_0 + _S4649.differential_0 + make_float3 (_S4705.x, _S4705.y, _S4706);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4727;
    (&_S4727)->primal_0 = R_31;
    (&_S4727)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4728;
    (&_S4728)->primal_0 = vert1_7;
    (&_S4728)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4727, &_S4728, _S4726);
    float3  _S4729 = _S4570 + _S4571 + _S4649.differential_0 + make_float3 (_S4721.x, _S4721.y, _S4722);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4730;
    (&_S4730)->primal_0 = R_31;
    (&_S4730)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4731;
    (&_S4731)->primal_0 = vert0_7;
    (&_S4731)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4730, &_S4731, _S4729);
    float3  _S4732 = (*v_verts_0)[int(2)] + _S4725.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4733;
    (&_S4733)->primal_0 = _S4451;
    (&_S4733)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4734;
    (&_S4734)->primal_0 = _S4456;
    (&_S4734)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4733, &_S4734, _S4732);
    float _S4735 = - _S4734.differential_0.y;
    float _S4736 = _S4455 * _S4734.differential_0.x;
    float _S4737 = - (_S4447 * _S4734.differential_0.x);
    float3  _S4738 = (*v_verts_0)[int(1)] + _S4728.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4739;
    (&_S4739)->primal_0 = _S4451;
    (&_S4739)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4740;
    (&_S4740)->primal_0 = _S4454;
    (&_S4740)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4739, &_S4740, _S4738);
    float _S4741 = _S4447 * _S4740.differential_0.x;
    float _S4742 = _S4453 * _S4740.differential_0.x;
    float3  _S4743 = (*v_verts_0)[int(0)] + _S4731.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4744;
    (&_S4744)->primal_0 = _S4451;
    (&_S4744)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4745;
    (&_S4745)->primal_0 = _S4452;
    (&_S4745)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4744, &_S4745, _S4743);
    Matrix<float, 3, 3>  _S4746 = transpose_0(_S4733.differential_0 + _S4739.differential_0 + _S4744.differential_0);
    float _S4747 = 2.0f * - _S4746.rows[int(2)].z;
    float _S4748 = 2.0f * _S4746.rows[int(2)].y;
    float _S4749 = 2.0f * _S4746.rows[int(2)].x;
    float _S4750 = 2.0f * _S4746.rows[int(1)].z;
    float _S4751 = 2.0f * - _S4746.rows[int(1)].y;
    float _S4752 = 2.0f * _S4746.rows[int(1)].x;
    float _S4753 = 2.0f * _S4746.rows[int(0)].z;
    float _S4754 = 2.0f * _S4746.rows[int(0)].y;
    float _S4755 = 2.0f * - _S4746.rows[int(0)].x;
    float _S4756 = - _S4752 + _S4754;
    float _S4757 = _S4749 + - _S4753;
    float _S4758 = - _S4748 + _S4750;
    float _S4759 = _S4748 + _S4750;
    float _S4760 = _S4749 + _S4753;
    float _S4761 = _S4752 + _S4754;
    float _S4762 = quat_31.w * (_S4751 + _S4755);
    float _S4763 = quat_31.z * (_S4747 + _S4755);
    float _S4764 = quat_31.y * (_S4747 + _S4751);
    float _S4765 = quat_31.x * _S4756 + quat_31.z * _S4759 + quat_31.y * _S4760 + _S4762 + _S4762;
    float _S4766 = quat_31.x * _S4757 + quat_31.w * _S4759 + quat_31.y * _S4761 + _S4763 + _S4763;
    float _S4767 = quat_31.x * _S4758 + quat_31.w * _S4760 + quat_31.z * _S4761 + _S4764 + _S4764;
    float _S4768 = quat_31.w * _S4756 + quat_31.z * _S4757 + quat_31.y * _S4758;
    float _S4769 = _S4737 + _S4741;
    float _S4770 = 0.5f * - _S4769;
    float _S4771 = _S4735 + _S4740.differential_0.y;
    DiffPair_float_0 _S4772;
    (&_S4772)->primal_0 = _S4448;
    (&_S4772)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4772, _S4771);
    float _S4773 = _S4770 + _S4772.differential_0;
    float _S4774 = _S4736 + _S4742 + _S4745.differential_0.x;
    DiffPair_float_0 _S4775;
    (&_S4775)->primal_0 = _S4446;
    (&_S4775)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4775, _S4774);
    float _S4776 = _S4770 + _S4775.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4777;
    (&_S4777)->primal_0 = R_31;
    (&_S4777)->differential_0 = _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4778;
    (&_S4778)->primal_0 = mean_30;
    (&_S4778)->differential_0 = _S4563;
    s_bwd_prop_mul_1(&_S4777, &_S4778, _S4565.differential_0);
    float3  _S4779 = _S4645.differential_0 + _S4723 + _S4726 + _S4729 + _S4565.differential_0;
    Matrix<float, 3, 3>  _S4780 = _S4646 + _S4724.differential_0 + _S4727.differential_0 + _S4730.differential_0 + _S4777.differential_0;
    FixedArray<float3 , 2>  _S4781;
    _S4781[int(0)] = _S4563;
    _S4781[int(1)] = _S4563;
    _S4781[int(1)] = _S4577;
    _S4781[int(0)] = _S4582;
    FixedArray<float3 , 16>  _S4782;
    _S4782[int(0)] = _S4563;
    _S4782[int(1)] = _S4563;
    _S4782[int(2)] = _S4563;
    _S4782[int(3)] = _S4563;
    _S4782[int(4)] = _S4563;
    _S4782[int(5)] = _S4563;
    _S4782[int(6)] = _S4563;
    _S4782[int(7)] = _S4563;
    _S4782[int(8)] = _S4563;
    _S4782[int(9)] = _S4563;
    _S4782[int(10)] = _S4563;
    _S4782[int(11)] = _S4563;
    _S4782[int(12)] = _S4563;
    _S4782[int(13)] = _S4563;
    _S4782[int(14)] = _S4563;
    _S4782[int(15)] = _S4563;
    _S4782[int(15)] = _S4584;
    _S4782[int(14)] = _S4586;
    _S4782[int(13)] = _S4588;
    _S4782[int(12)] = _S4590;
    _S4782[int(11)] = _S4592;
    _S4782[int(10)] = _S4594;
    _S4782[int(9)] = _S4596;
    _S4782[int(8)] = _S4604;
    _S4782[int(7)] = _S4606;
    _S4782[int(6)] = _S4608;
    _S4782[int(5)] = _S4610;
    _S4782[int(4)] = _S4612;
    _S4782[int(3)] = _S4623;
    _S4782[int(2)] = _S4625;
    _S4782[int(1)] = _S4627;
    _S4782[int(0)] = _S4640;
    float2  _S4783 = make_float2 (0.0f, _S4667);
    float3  _S4784 = make_float3 (_S4776, _S4773, _S4769);
    float4  _S4785 = make_float4 (0.0f);
    *&((&_S4785)->w) = _S4765;
    *&((&_S4785)->z) = _S4766;
    *&((&_S4785)->y) = _S4767;
    *&((&_S4785)->x) = _S4768;
    *v_mean_9 = _S4641 + _S4732 + _S4738 + _S4743 + _S4778.differential_0;
    *v_quat_8 = _S4785;
    *v_scale_8 = _S4784;
    *v_hardness_4 = _S4783;
    *v_sh_coeffs_7 = _S4782;
    *v_ch_coeffs_2 = _S4781;
    *v_R_9 = _S4780;
    *v_t_8 = _S4779;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_31, float4  quat_32, float3  scale_31, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_32, float3  t_29, float fx_35, float fy_35, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_39, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_0(R_32, mean_31) + t_29;
    float _S4786 = scale_31.x;
    float _S4787 = s_primal_ctx_exp_1(_S4786);
    float _S4788 = scale_31.y;
    float _S4789 = s_primal_ctx_exp_1(_S4788);
    float sz_12 = scale_31.z - 0.5f * (_S4786 + _S4788);
    float _S4790 = quat_32.y;
    float x2_32 = _S4790 * _S4790;
    float y2_32 = quat_32.z * quat_32.z;
    float z2_58 = quat_32.w * quat_32.w;
    float xy_32 = quat_32.y * quat_32.z;
    float xz_32 = quat_32.y * quat_32.w;
    float yz_32 = quat_32.z * quat_32.w;
    float wx_32 = quat_32.x * quat_32.y;
    float wy_32 = quat_32.x * quat_32.z;
    float wz_32 = quat_32.x * quat_32.w;
    Matrix<float, 3, 3>  _S4791 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_32 + z2_58), 2.0f * (xy_32 + wz_32), 2.0f * (xz_32 - wy_32), 2.0f * (xy_32 - wz_32), 1.0f - 2.0f * (x2_32 + z2_58), 2.0f * (yz_32 + wx_32), 2.0f * (xz_32 + wy_32), 2.0f * (yz_32 - wx_32), 1.0f - 2.0f * (x2_32 + y2_32)));
    float3  _S4792 = make_float3 (_S4787, 0.0f, 0.0f);
    float3  vert0_8 = s_primal_ctx_mul_0(_S4791, _S4792) + mean_31;
    float _S4793 = -0.5f + sz_12;
    float3  _S4794 = make_float3 (_S4787 * _S4793, _S4789, 0.0f);
    float3  vert1_8 = s_primal_ctx_mul_0(_S4791, _S4794) + mean_31;
    float _S4795 = -0.5f - sz_12;
    float3  _S4796 = make_float3 (_S4787 * _S4795, - _S4789, 0.0f);
    float3  vert2_8 = s_primal_ctx_mul_0(_S4791, _S4796) + mean_31;
    float3  vert0_c_11 = s_primal_ctx_mul_0(R_32, vert0_8) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_0(R_32, vert1_8) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_0(R_32, vert2_8) + t_29;
    float2  _S4797 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4798 = length_0(_S4797);
    float _S4799 = vert0_c_11.z;
    float _S4800 = s_primal_ctx_atan2_0(_S4798, _S4799);
    bool _S4801 = _S4800 < 0.00100000004749745f;
    float k_12;
    float _S4802;
    float _S4803;
    float _S4804;
    if(_S4801)
    {
        float _S4805 = 1.0f - _S4800 * _S4800 / 3.0f;
        float _S4806 = _S4799 * _S4799;
        k_12 = _S4805 / _S4799;
        _S4802 = 0.0f;
        _S4803 = _S4806;
        _S4804 = _S4805;
    }
    else
    {
        float _S4807 = _S4798 * _S4798;
        k_12 = _S4800 / _S4798;
        _S4802 = _S4807;
        _S4803 = 0.0f;
        _S4804 = 0.0f;
    }
    float2  _S4808 = make_float2 (k_12);
    float2  _S4809 = _S4797 * make_float2 (k_12);
    float u_77 = _S4809.x;
    float v_77 = _S4809.y;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float _S4810 = (*dist_coeffs_39)[int(2)] + r2_77 * (*dist_coeffs_39)[int(3)];
    float _S4811 = (*dist_coeffs_39)[int(1)] + r2_77 * _S4810;
    float _S4812 = (*dist_coeffs_39)[int(0)] + r2_77 * _S4811;
    float radial_10 = 1.0f + r2_77 * _S4812;
    float2  _S4813 = make_float2 (radial_10);
    float _S4814 = 2.0f * (*dist_coeffs_39)[int(4)];
    float _S4815 = _S4814 * u_77;
    float _S4816 = 2.0f * u_77;
    float _S4817 = 2.0f * (*dist_coeffs_39)[int(5)];
    float _S4818 = _S4817 * u_77;
    float _S4819 = 2.0f * v_77;
    float2  _S4820 = _S4809 * make_float2 (radial_10) + make_float2 (_S4815 * v_77 + (*dist_coeffs_39)[int(5)] * (r2_77 + _S4816 * u_77) + (*dist_coeffs_39)[int(6)] * r2_77, _S4818 * v_77 + (*dist_coeffs_39)[int(4)] * (r2_77 + _S4819 * v_77) + (*dist_coeffs_39)[int(7)] * r2_77);
    float2  _S4821 = _S4820 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4820.x + (*dist_coeffs_39)[int(9)] * _S4820.y, 0.0f);
    float _S4822 = fx_35 * _S4821.x + cx_30;
    float _S4823 = fy_35 * _S4821.y + cy_30;
    float2  uv0_11 = make_float2 (_S4822, _S4823);
    float2  _S4824 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4825 = length_0(_S4824);
    float _S4826 = vert1_c_11.z;
    float _S4827 = s_primal_ctx_atan2_0(_S4825, _S4826);
    bool _S4828 = _S4827 < 0.00100000004749745f;
    float _S4829;
    float _S4830;
    float _S4831;
    if(_S4828)
    {
        float _S4832 = 1.0f - _S4827 * _S4827 / 3.0f;
        float _S4833 = _S4826 * _S4826;
        k_12 = _S4832 / _S4826;
        _S4829 = 0.0f;
        _S4830 = _S4833;
        _S4831 = _S4832;
    }
    else
    {
        float _S4834 = _S4825 * _S4825;
        k_12 = _S4827 / _S4825;
        _S4829 = _S4834;
        _S4830 = 0.0f;
        _S4831 = 0.0f;
    }
    float2  _S4835 = make_float2 (k_12);
    float2  _S4836 = _S4824 * make_float2 (k_12);
    float u_78 = _S4836.x;
    float v_78 = _S4836.y;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float _S4837 = (*dist_coeffs_39)[int(2)] + r2_78 * (*dist_coeffs_39)[int(3)];
    float _S4838 = (*dist_coeffs_39)[int(1)] + r2_78 * _S4837;
    float _S4839 = (*dist_coeffs_39)[int(0)] + r2_78 * _S4838;
    float radial_11 = 1.0f + r2_78 * _S4839;
    float2  _S4840 = make_float2 (radial_11);
    float _S4841 = _S4814 * u_78;
    float _S4842 = 2.0f * u_78;
    float _S4843 = _S4817 * u_78;
    float _S4844 = 2.0f * v_78;
    float2  _S4845 = _S4836 * make_float2 (radial_11) + make_float2 (_S4841 * v_78 + (*dist_coeffs_39)[int(5)] * (r2_78 + _S4842 * u_78) + (*dist_coeffs_39)[int(6)] * r2_78, _S4843 * v_78 + (*dist_coeffs_39)[int(4)] * (r2_78 + _S4844 * v_78) + (*dist_coeffs_39)[int(7)] * r2_78);
    float2  _S4846 = _S4845 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4845.x + (*dist_coeffs_39)[int(9)] * _S4845.y, 0.0f);
    float _S4847 = fx_35 * _S4846.x + cx_30;
    float _S4848 = fy_35 * _S4846.y + cy_30;
    float2  uv1_11 = make_float2 (_S4847, _S4848);
    float2  _S4849 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4850 = length_0(_S4849);
    float _S4851 = vert2_c_11.z;
    float _S4852 = s_primal_ctx_atan2_0(_S4850, _S4851);
    bool _S4853 = _S4852 < 0.00100000004749745f;
    float _S4854;
    float _S4855;
    float _S4856;
    if(_S4853)
    {
        float _S4857 = 1.0f - _S4852 * _S4852 / 3.0f;
        float _S4858 = _S4851 * _S4851;
        k_12 = _S4857 / _S4851;
        _S4854 = 0.0f;
        _S4855 = _S4858;
        _S4856 = _S4857;
    }
    else
    {
        float _S4859 = _S4850 * _S4850;
        k_12 = _S4852 / _S4850;
        _S4854 = _S4859;
        _S4855 = 0.0f;
        _S4856 = 0.0f;
    }
    float2  _S4860 = make_float2 (k_12);
    float2  _S4861 = _S4849 * make_float2 (k_12);
    float u_79 = _S4861.x;
    float v_79 = _S4861.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S4862 = (*dist_coeffs_39)[int(2)] + r2_79 * (*dist_coeffs_39)[int(3)];
    float _S4863 = (*dist_coeffs_39)[int(1)] + r2_79 * _S4862;
    float _S4864 = (*dist_coeffs_39)[int(0)] + r2_79 * _S4863;
    float radial_12 = 1.0f + r2_79 * _S4864;
    float2  _S4865 = make_float2 (radial_12);
    float _S4866 = _S4814 * u_79;
    float _S4867 = 2.0f * u_79;
    float _S4868 = _S4817 * u_79;
    float _S4869 = 2.0f * v_79;
    float2  _S4870 = _S4861 * make_float2 (radial_12) + make_float2 (_S4866 * v_79 + (*dist_coeffs_39)[int(5)] * (r2_79 + _S4867 * u_79) + (*dist_coeffs_39)[int(6)] * r2_79, _S4868 * v_79 + (*dist_coeffs_39)[int(4)] * (r2_79 + _S4869 * v_79) + (*dist_coeffs_39)[int(7)] * r2_79);
    float2  _S4871 = _S4870 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4870.x + (*dist_coeffs_39)[int(9)] * _S4870.y, 0.0f);
    float _S4872 = fx_35 * _S4871.x + cx_30;
    float _S4873 = fy_35 * _S4871.y + cy_30;
    float2  uv2_11 = make_float2 (_S4872, _S4873);
    float2  e0_15 = uv1_11 - uv0_11;
    float2  e1_15 = uv2_11 - uv1_11;
    float2  e2_7 = uv0_11 - uv2_11;
    float _S4874 = e0_15.x;
    float _S4875 = e1_15.y;
    float _S4876 = e0_15.y;
    float _S4877 = e1_15.x;
    float _S4878 = _S4874 * _S4875 - _S4876 * _S4877;
    float _S4879 = 1.0f - hardness_15.y;
    float _S4880 = -1.0f / _S4879;
    float _S4881 = _S4879 * _S4879;
    float _S4882 = s_primal_ctx_max_0(_S4822, _S4847);
    float _S4883 = s_primal_ctx_min_0(_S4822, _S4847);
    float _S4884 = s_primal_ctx_max_0(_S4823, _S4848);
    float _S4885 = s_primal_ctx_min_0(_S4823, _S4848);
    float3  _S4886 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4887 = length_1(_S4886) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4888 = transpose_0(R_32);
    float3  _S4889 = mean_31 - - s_primal_ctx_mul_0(_S4888, t_29);
    float _S4890 = _S4889.x;
    float _S4891 = _S4889.y;
    float _S4892 = _S4889.z;
    float _S4893 = _S4890 * _S4890 + _S4891 * _S4891 + _S4892 * _S4892;
    float _S4894 = s_primal_ctx_sqrt_0(_S4893);
    float x_62 = _S4890 / _S4894;
    float3  _S4895 = make_float3 (x_62);
    float _S4896 = _S4894 * _S4894;
    float y_29 = _S4891 / _S4894;
    float z_26 = _S4892 / _S4894;
    float3  _S4897 = make_float3 (z_26);
    float _S4898 = - y_29;
    float3  _S4899 = make_float3 (_S4898);
    float z2_59 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_62 * x_62 - y_29 * y_29;
    float _S4900 = 2.0f * x_62;
    float fS1_26 = _S4900 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_59 - 0.31539157032966614f;
    float3  _S4901 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_62;
    float3  _S4902 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4903 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4904 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4905 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_59 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4906 = 1.86588168144226074f * z2_59 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4906;
    float3  _S4907 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_62;
    float3  _S4908 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4909 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4910 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4911 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_62 * fC1_26 - y_29 * fS1_26);
    float3  _S4912 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_62 * fS1_26 + y_29 * fC1_26);
    float3  _S4913 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4898) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_62) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4914 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4915 = make_float3 (0.0f);
    float3  _S4916 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4917 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4918 = make_float3 (_S4917);
    float3  _S4919 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4917);
    float3  _S4920 = _S4916 + _S4919 + make_float3 (0.5f);
    float3  _S4921 = _S4916 - _S4919 + make_float3 (0.5f);
    float3  _S4922 = vert1_c_11 - vert0_c_11;
    float3  _S4923 = vert2_c_11 - vert0_c_11;
    float3  _S4924 = s_primal_ctx_cross_0(_S4922, _S4923);
    float3  _S4925 = normalize_0(_S4924);
    float3  _S4926 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4925, mean_c_26)))))) * v_normal_3;
    float3  _S4927 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4928;
    (&_S4928)->primal_0 = _S4925;
    (&_S4928)->differential_0 = _S4927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4929;
    (&_S4929)->primal_0 = mean_c_26;
    (&_S4929)->differential_0 = _S4927;
    s_bwd_prop_dot_0(&_S4928, &_S4929, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4930 = _S4929;
    float3  _S4931 = _S4926 + _S4928.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4932;
    (&_S4932)->primal_0 = _S4924;
    (&_S4932)->differential_0 = _S4927;
    s_bwd_normalize_impl_0(&_S4932, _S4931);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4933;
    (&_S4933)->primal_0 = _S4922;
    (&_S4933)->differential_0 = _S4927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4934;
    (&_S4934)->primal_0 = _S4923;
    (&_S4934)->differential_0 = _S4927;
    s_bwd_prop_cross_0(&_S4933, &_S4934, _S4932.differential_0);
    float3  _S4935 = - _S4934.differential_0;
    float3  _S4936 = - _S4933.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4937;
    (&_S4937)->primal_0 = _S4921;
    (&_S4937)->differential_0 = _S4927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4938;
    (&_S4938)->primal_0 = _S4915;
    (&_S4938)->differential_0 = _S4927;
    s_bwd_prop_max_0(&_S4937, &_S4938, (*v_rgbs_1)[int(2)]);
    float3  _S4939 = - _S4937.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4940;
    (&_S4940)->primal_0 = _S4920;
    (&_S4940)->differential_0 = _S4927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4941;
    (&_S4941)->primal_0 = _S4915;
    (&_S4941)->differential_0 = _S4927;
    s_bwd_prop_max_0(&_S4940, &_S4941, (*v_rgbs_1)[int(1)]);
    float3  _S4942 = _S4918 * (_S4939 + _S4940.differential_0);
    float3  _S4943 = _S4937.differential_0 + _S4940.differential_0;
    float3  _S4944 = make_float3 (0.5f) * - _S4943;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4945;
    (&_S4945)->primal_0 = _S4914;
    (&_S4945)->differential_0 = _S4927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4946;
    (&_S4946)->primal_0 = _S4915;
    (&_S4946)->differential_0 = _S4927;
    s_bwd_prop_max_0(&_S4945, &_S4946, (*v_rgbs_1)[int(0)]);
    float3  _S4947 = _S4944 + _S4945.differential_0;
    float3  _S4948 = _S4943 + _S4945.differential_0;
    float3  _S4949 = _S4912 * _S4948;
    float3  _S4950 = (*sh_coeffs_26)[int(15)] * _S4948;
    float3  _S4951 = _S4910 * _S4948;
    float3  _S4952 = (*sh_coeffs_26)[int(14)] * _S4948;
    float3  _S4953 = _S4908 * _S4948;
    float3  _S4954 = (*sh_coeffs_26)[int(13)] * _S4948;
    float3  _S4955 = _S4907 * _S4948;
    float3  _S4956 = (*sh_coeffs_26)[int(12)] * _S4948;
    float3  _S4957 = _S4909 * _S4948;
    float3  _S4958 = (*sh_coeffs_26)[int(11)] * _S4948;
    float3  _S4959 = _S4911 * _S4948;
    float3  _S4960 = (*sh_coeffs_26)[int(10)] * _S4948;
    float3  _S4961 = _S4913 * _S4948;
    float3  _S4962 = (*sh_coeffs_26)[int(9)] * _S4948;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4962.x + _S4962.y + _S4962.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4950.x + _S4950.y + _S4950.z);
    float _S4963 = _S4960.x + _S4960.y + _S4960.z;
    float _S4964 = _S4952.x + _S4952.y + _S4952.z;
    float _S4965 = _S4958.x + _S4958.y + _S4958.z;
    float _S4966 = _S4954.x + _S4954.y + _S4954.z;
    float _S4967 = _S4956.x + _S4956.y + _S4956.z;
    float _S4968 = - s_diff_fC2_T_8;
    float3  _S4969 = _S4904 * _S4948;
    float3  _S4970 = (*sh_coeffs_26)[int(8)] * _S4948;
    float3  _S4971 = _S4902 * _S4948;
    float3  _S4972 = (*sh_coeffs_26)[int(7)] * _S4948;
    float3  _S4973 = _S4901 * _S4948;
    float3  _S4974 = (*sh_coeffs_26)[int(6)] * _S4948;
    float3  _S4975 = _S4903 * _S4948;
    float3  _S4976 = (*sh_coeffs_26)[int(5)] * _S4948;
    float3  _S4977 = _S4905 * _S4948;
    float3  _S4978 = (*sh_coeffs_26)[int(4)] * _S4948;
    float _S4979 = _S4976.x + _S4976.y + _S4976.z;
    float _S4980 = _S4972.x + _S4972.y + _S4972.z;
    float _S4981 = fTmp1B_26 * _S4963 + x_62 * s_diff_fS2_T_8 + y_29 * _S4968 + 0.54627424478530884f * (_S4978.x + _S4978.y + _S4978.z);
    float _S4982 = fTmp1B_26 * _S4964 + y_29 * s_diff_fS2_T_8 + x_62 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4970.x + _S4970.y + _S4970.z);
    float _S4983 = y_29 * - _S4982;
    float _S4984 = x_62 * _S4982;
    float _S4985 = z_26 * (1.86588168144226074f * (z_26 * _S4967) + -2.28522896766662598f * (y_29 * _S4965 + x_62 * _S4966) + 0.94617468118667603f * (_S4974.x + _S4974.y + _S4974.z));
    float3  _S4986 = make_float3 (0.48860251903533936f) * _S4948;
    float3  _S4987 = - _S4986;
    float3  _S4988 = _S4895 * _S4987;
    float3  _S4989 = (*sh_coeffs_26)[int(3)] * _S4987;
    float3  _S4990 = _S4897 * _S4986;
    float3  _S4991 = (*sh_coeffs_26)[int(2)] * _S4986;
    float3  _S4992 = _S4899 * _S4986;
    float3  _S4993 = (*sh_coeffs_26)[int(1)] * _S4986;
    float _S4994 = (_S4906 * _S4967 + 1.44530570507049561f * (fS1_26 * _S4963 + fC1_26 * _S4964) + -1.09254848957061768f * (y_29 * _S4979 + x_62 * _S4980) + _S4985 + _S4985 + _S4991.x + _S4991.y + _S4991.z) / _S4896;
    float _S4995 = _S4894 * _S4994;
    float _S4996 = (fTmp0C_26 * _S4965 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4968 + fTmp0B_26 * _S4979 + _S4900 * _S4981 + _S4983 + _S4983 + - (_S4993.x + _S4993.y + _S4993.z)) / _S4896;
    float _S4997 = _S4894 * _S4996;
    float _S4998 = (fTmp0C_26 * _S4966 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4980 + 2.0f * (y_29 * _S4981) + _S4984 + _S4984 + _S4989.x + _S4989.y + _S4989.z) / _S4896;
    float _S4999 = _S4894 * _S4998;
    float _S5000 = _S4892 * - _S4994 + _S4891 * - _S4996 + _S4890 * - _S4998;
    DiffPair_float_0 _S5001;
    (&_S5001)->primal_0 = _S4893;
    (&_S5001)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5001, _S5000);
    float _S5002 = _S4892 * _S5001.differential_0;
    float _S5003 = _S4891 * _S5001.differential_0;
    float _S5004 = _S4890 * _S5001.differential_0;
    float3  _S5005 = make_float3 (0.282094806432724f) * _S4948;
    float3  _S5006 = make_float3 (_S4999 + _S5004 + _S5004, _S4997 + _S5003 + _S5003, _S4995 + _S5002 + _S5002);
    float3  _S5007 = - - _S5006;
    Matrix<float, 3, 3>  _S5008 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5009;
    (&_S5009)->primal_0 = _S4888;
    (&_S5009)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5010;
    (&_S5010)->primal_0 = t_29;
    (&_S5010)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5009, &_S5010, _S5007);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5011 = _S5010;
    Matrix<float, 3, 3>  _S5012 = transpose_0(_S5009.differential_0);
    DiffPair_float_0 _S5013;
    (&_S5013)->primal_0 = _S4887;
    (&_S5013)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5013, v_depth_10);
    float _S5014 = 0.3333333432674408f * _S5013.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5015;
    (&_S5015)->primal_0 = _S4886;
    (&_S5015)->differential_0 = _S4927;
    s_bwd_length_impl_1(&_S5015, _S5014);
    DiffPair_float_0 _S5016;
    (&_S5016)->primal_0 = _S4885;
    (&_S5016)->differential_0 = 0.0f;
    DiffPair_float_0 _S5017;
    (&_S5017)->primal_0 = _S4873;
    (&_S5017)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5016, &_S5017, 0.0f);
    DiffPair_float_0 _S5018;
    (&_S5018)->primal_0 = _S4823;
    (&_S5018)->differential_0 = 0.0f;
    DiffPair_float_0 _S5019;
    (&_S5019)->primal_0 = _S4848;
    (&_S5019)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5018, &_S5019, _S5016.differential_0);
    DiffPair_float_0 _S5020;
    (&_S5020)->primal_0 = _S4884;
    (&_S5020)->differential_0 = 0.0f;
    DiffPair_float_0 _S5021;
    (&_S5021)->primal_0 = _S4873;
    (&_S5021)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5020, &_S5021, 0.0f);
    DiffPair_float_0 _S5022;
    (&_S5022)->primal_0 = _S4823;
    (&_S5022)->differential_0 = 0.0f;
    DiffPair_float_0 _S5023;
    (&_S5023)->primal_0 = _S4848;
    (&_S5023)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5022, &_S5023, _S5020.differential_0);
    DiffPair_float_0 _S5024;
    (&_S5024)->primal_0 = _S4883;
    (&_S5024)->differential_0 = 0.0f;
    DiffPair_float_0 _S5025;
    (&_S5025)->primal_0 = _S4872;
    (&_S5025)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5024, &_S5025, 0.0f);
    DiffPair_float_0 _S5026;
    (&_S5026)->primal_0 = _S4822;
    (&_S5026)->differential_0 = 0.0f;
    DiffPair_float_0 _S5027;
    (&_S5027)->primal_0 = _S4847;
    (&_S5027)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5026, &_S5027, _S5024.differential_0);
    DiffPair_float_0 _S5028;
    (&_S5028)->primal_0 = _S4882;
    (&_S5028)->differential_0 = 0.0f;
    DiffPair_float_0 _S5029;
    (&_S5029)->primal_0 = _S4872;
    (&_S5029)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5028, &_S5029, 0.0f);
    DiffPair_float_0 _S5030;
    (&_S5030)->primal_0 = _S4822;
    (&_S5030)->differential_0 = 0.0f;
    DiffPair_float_0 _S5031;
    (&_S5031)->primal_0 = _S4847;
    (&_S5031)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5030, &_S5031, _S5028.differential_0);
    DiffPair_float_0 _S5032;
    (&_S5032)->primal_0 = _S4880;
    (&_S5032)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5032, 0.0f);
    float _S5033 = - (-1.0f * - (_S5032.differential_0 / _S4881));
    float2  _S5034 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5035;
    (&_S5035)->primal_0 = e2_7;
    (&_S5035)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5035, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5036;
    (&_S5036)->primal_0 = e1_15;
    (&_S5036)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5036, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5037;
    (&_S5037)->primal_0 = e0_15;
    (&_S5037)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5037, -0.0f);
    DiffPair_float_0 _S5038;
    (&_S5038)->primal_0 = _S4878;
    (&_S5038)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5038, 0.0f);
    float _S5039 = - _S5038.differential_0;
    float2  _S5040 = _S5036.differential_0 + make_float2 (_S4876 * _S5039, _S4874 * _S5038.differential_0);
    float2  _S5041 = - _S5040;
    float2  _S5042 = _S5037.differential_0 + make_float2 (_S4875 * _S5038.differential_0, _S4877 * _S5039);
    float2  _S5043 = - _S5042;
    float2  _S5044 = - _S5035.differential_0 + _S5040;
    float _S5045 = fx_35 * (_S5025.differential_0 + _S5029.differential_0 + _S5044.x);
    float2  _S5046 = make_float2 (_S5045, fy_35 * (_S5017.differential_0 + _S5021.differential_0 + _S5044.y)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5045, (*dist_coeffs_39)[int(9)] * _S5045);
    float2  _S5047 = _S4861 * _S5046;
    float2  _S5048 = _S4865 * _S5046;
    float _S5049 = (*dist_coeffs_39)[int(4)] * _S5046.y;
    float _S5050 = (*dist_coeffs_39)[int(5)] * _S5046.x;
    float _S5051 = _S5047.x + _S5047.y;
    float _S5052 = r2_79 * _S5051;
    float _S5053 = r2_79 * _S5052;
    float _S5054 = (*dist_coeffs_39)[int(7)] * _S5046.y + _S5049 + (*dist_coeffs_39)[int(6)] * _S5046.x + _S5050 + _S4864 * _S5051 + _S4863 * _S5052 + _S4862 * _S5053 + (*dist_coeffs_39)[int(3)] * (r2_79 * _S5053);
    float _S5055 = v_79 * _S5054;
    float _S5056 = u_79 * _S5054;
    float _S5057 = _S4869 * _S5049 + 2.0f * (v_79 * _S5049) + _S4868 * _S5046.y + _S4866 * _S5046.x + _S5055 + _S5055;
    float _S5058 = _S4817 * (v_79 * _S5046.y) + _S4867 * _S5050 + 2.0f * (u_79 * _S5050) + _S4814 * (v_79 * _S5046.x) + _S5056 + _S5056;
    float3  _S5059 = _S4933.differential_0 + _S5015.differential_0;
    float3  _S5060 = _S4935 + _S4936 + _S5015.differential_0;
    float3  _S5061 = _S4934.differential_0 + _S5015.differential_0;
    FixedArray<float3 , 2>  _S5062;
    _S5062[int(0)] = _S4927;
    _S5062[int(1)] = _S4927;
    _S5062[int(1)] = _S4942;
    _S5062[int(0)] = _S4947;
    float3  _S5063 = _S5062[int(0)];
    float3  _S5064 = _S5062[int(1)];
    FixedArray<float3 , 16>  _S5065;
    _S5065[int(0)] = _S4927;
    _S5065[int(1)] = _S4927;
    _S5065[int(2)] = _S4927;
    _S5065[int(3)] = _S4927;
    _S5065[int(4)] = _S4927;
    _S5065[int(5)] = _S4927;
    _S5065[int(6)] = _S4927;
    _S5065[int(7)] = _S4927;
    _S5065[int(8)] = _S4927;
    _S5065[int(9)] = _S4927;
    _S5065[int(10)] = _S4927;
    _S5065[int(11)] = _S4927;
    _S5065[int(12)] = _S4927;
    _S5065[int(13)] = _S4927;
    _S5065[int(14)] = _S4927;
    _S5065[int(15)] = _S4927;
    _S5065[int(7)] = _S4971;
    _S5065[int(0)] = _S5005;
    _S5065[int(1)] = _S4992;
    _S5065[int(2)] = _S4990;
    _S5065[int(3)] = _S4988;
    _S5065[int(4)] = _S4977;
    _S5065[int(5)] = _S4975;
    _S5065[int(6)] = _S4973;
    _S5065[int(15)] = _S4949;
    _S5065[int(8)] = _S4969;
    _S5065[int(9)] = _S4961;
    _S5065[int(10)] = _S4959;
    _S5065[int(11)] = _S4957;
    _S5065[int(12)] = _S4955;
    _S5065[int(13)] = _S4953;
    _S5065[int(14)] = _S4951;
    float3  _S5066 = _S5065[int(0)];
    float3  _S5067 = _S5065[int(1)];
    float3  _S5068 = _S5065[int(2)];
    float3  _S5069 = _S5065[int(3)];
    float3  _S5070 = _S5065[int(4)];
    float3  _S5071 = _S5065[int(5)];
    float3  _S5072 = _S5065[int(6)];
    float3  _S5073 = _S5065[int(7)];
    float3  _S5074 = _S5065[int(8)];
    float3  _S5075 = _S5065[int(9)];
    float3  _S5076 = _S5065[int(10)];
    float3  _S5077 = _S5065[int(11)];
    float3  _S5078 = _S5065[int(12)];
    float3  _S5079 = _S5065[int(13)];
    float3  _S5080 = _S5065[int(14)];
    float3  _S5081 = _S5065[int(15)];
    float _S5082 = _S5026.differential_0 + _S5030.differential_0;
    float2  _S5083 = _S5035.differential_0 + _S5043;
    float _S5084 = _S5019.differential_0 + _S5023.differential_0;
    float _S5085 = _S5018.differential_0 + _S5022.differential_0;
    float2  _S5086 = _S5041 + _S5042;
    float _S5087 = _S5027.differential_0 + _S5031.differential_0;
    float2  _S5088 = make_float2 (0.0f, _S5033);
    float2  _S5089 = _S5048 + make_float2 (_S5058, _S5057);
    float2  _S5090 = _S4849 * _S5089;
    float2  _S5091 = _S4860 * _S5089;
    float _S5092 = _S5090.x + _S5090.y;
    if(_S4853)
    {
        float _S5093 = _S5092 / _S4855;
        float _S5094 = _S4856 * - _S5093;
        float _S5095 = _S4852 * (0.3333333432674408f * - (_S4851 * _S5093));
        k_12 = _S5095 + _S5095;
        _S4854 = _S5094;
        _S4855 = 0.0f;
    }
    else
    {
        float _S5096 = _S5092 / _S4854;
        float _S5097 = _S4852 * - _S5096;
        k_12 = _S4850 * _S5096;
        _S4854 = 0.0f;
        _S4855 = _S5097;
    }
    DiffPair_float_0 _S5098;
    (&_S5098)->primal_0 = _S4850;
    (&_S5098)->differential_0 = 0.0f;
    DiffPair_float_0 _S5099;
    (&_S5099)->primal_0 = _S4851;
    (&_S5099)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5098, &_S5099, k_12);
    float _S5100 = _S5099.differential_0 + _S4854;
    float _S5101 = _S5098.differential_0 + _S4855;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5102;
    (&_S5102)->primal_0 = _S4849;
    (&_S5102)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5102, _S5101);
    float2  _S5103 = _S5102.differential_0 + _S5091;
    float _S5104 = fx_35 * (_S5086.x + _S5087);
    float2  _S5105 = make_float2 (_S5104, fy_35 * (_S5086.y + _S5084)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5104, (*dist_coeffs_39)[int(9)] * _S5104);
    float2  _S5106 = _S4836 * _S5105;
    float _S5107 = (*dist_coeffs_39)[int(4)] * _S5105.y;
    float _S5108 = (*dist_coeffs_39)[int(5)] * _S5105.x;
    float _S5109 = _S5106.x + _S5106.y;
    float _S5110 = r2_78 * _S5109;
    float _S5111 = r2_78 * _S5110;
    float _S5112 = (*dist_coeffs_39)[int(7)] * _S5105.y + _S5107 + (*dist_coeffs_39)[int(6)] * _S5105.x + _S5108 + _S4839 * _S5109 + _S4838 * _S5110 + _S4837 * _S5111 + (*dist_coeffs_39)[int(3)] * (r2_78 * _S5111);
    float _S5113 = v_78 * _S5112;
    float _S5114 = u_78 * _S5112;
    float3  _S5115 = _S5061 + make_float3 (_S5103.x, _S5103.y, _S5100);
    float2  _S5116 = _S4840 * _S5105 + make_float2 (_S4817 * (v_78 * _S5105.y) + _S4842 * _S5108 + 2.0f * (u_78 * _S5108) + _S4814 * (v_78 * _S5105.x) + _S5114 + _S5114, _S4844 * _S5107 + 2.0f * (v_78 * _S5107) + _S4843 * _S5105.y + _S4841 * _S5105.x + _S5113 + _S5113);
    float2  _S5117 = _S4824 * _S5116;
    float2  _S5118 = _S4835 * _S5116;
    float _S5119 = _S5117.x + _S5117.y;
    if(_S4828)
    {
        float _S5120 = _S5119 / _S4830;
        float _S5121 = _S4831 * - _S5120;
        float _S5122 = _S4827 * (0.3333333432674408f * - (_S4826 * _S5120));
        k_12 = _S5122 + _S5122;
        _S4829 = _S5121;
        _S4830 = 0.0f;
    }
    else
    {
        float _S5123 = _S5119 / _S4829;
        float _S5124 = _S4827 * - _S5123;
        k_12 = _S4825 * _S5123;
        _S4829 = 0.0f;
        _S4830 = _S5124;
    }
    DiffPair_float_0 _S5125;
    (&_S5125)->primal_0 = _S4825;
    (&_S5125)->differential_0 = 0.0f;
    DiffPair_float_0 _S5126;
    (&_S5126)->primal_0 = _S4826;
    (&_S5126)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5125, &_S5126, k_12);
    float _S5127 = _S5126.differential_0 + _S4829;
    float _S5128 = _S5125.differential_0 + _S4830;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5129;
    (&_S5129)->primal_0 = _S4824;
    (&_S5129)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5129, _S5128);
    float2  _S5130 = _S5129.differential_0 + _S5118;
    float _S5131 = fx_35 * (_S5083.x + _S5082);
    float2  _S5132 = make_float2 (_S5131, fy_35 * (_S5083.y + _S5085)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5131, (*dist_coeffs_39)[int(9)] * _S5131);
    float2  _S5133 = _S4809 * _S5132;
    float _S5134 = (*dist_coeffs_39)[int(4)] * _S5132.y;
    float _S5135 = (*dist_coeffs_39)[int(5)] * _S5132.x;
    float _S5136 = _S5133.x + _S5133.y;
    float _S5137 = r2_77 * _S5136;
    float _S5138 = r2_77 * _S5137;
    float _S5139 = (*dist_coeffs_39)[int(7)] * _S5132.y + _S5134 + (*dist_coeffs_39)[int(6)] * _S5132.x + _S5135 + _S4812 * _S5136 + _S4811 * _S5137 + _S4810 * _S5138 + (*dist_coeffs_39)[int(3)] * (r2_77 * _S5138);
    float _S5140 = v_77 * _S5139;
    float _S5141 = u_77 * _S5139;
    float3  _S5142 = _S5059 + make_float3 (_S5130.x, _S5130.y, _S5127);
    float2  _S5143 = _S4813 * _S5132 + make_float2 (_S4817 * (v_77 * _S5132.y) + _S4816 * _S5135 + 2.0f * (u_77 * _S5135) + _S4814 * (v_77 * _S5132.x) + _S5141 + _S5141, _S4819 * _S5134 + 2.0f * (v_77 * _S5134) + _S4818 * _S5132.y + _S4815 * _S5132.x + _S5140 + _S5140);
    float2  _S5144 = _S4797 * _S5143;
    float2  _S5145 = _S4808 * _S5143;
    float _S5146 = _S5144.x + _S5144.y;
    if(_S4801)
    {
        float _S5147 = _S5146 / _S4803;
        float _S5148 = _S4804 * - _S5147;
        float _S5149 = _S4800 * (0.3333333432674408f * - (_S4799 * _S5147));
        k_12 = _S5149 + _S5149;
        _S4802 = _S5148;
        _S4803 = 0.0f;
    }
    else
    {
        float _S5150 = _S5146 / _S4802;
        float _S5151 = _S4800 * - _S5150;
        k_12 = _S4798 * _S5150;
        _S4802 = 0.0f;
        _S4803 = _S5151;
    }
    DiffPair_float_0 _S5152;
    (&_S5152)->primal_0 = _S4798;
    (&_S5152)->differential_0 = 0.0f;
    DiffPair_float_0 _S5153;
    (&_S5153)->primal_0 = _S4799;
    (&_S5153)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5152, &_S5153, k_12);
    float _S5154 = _S5153.differential_0 + _S4802;
    float _S5155 = _S5152.differential_0 + _S4803;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5156;
    (&_S5156)->primal_0 = _S4797;
    (&_S5156)->differential_0 = _S5034;
    s_bwd_length_impl_0(&_S5156, _S5155);
    float2  _S5157 = _S5156.differential_0 + _S5145;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5158;
    (&_S5158)->primal_0 = vert2_c_11;
    (&_S5158)->differential_0 = _S4927;
    s_bwd_length_impl_1(&_S5158, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5159;
    (&_S5159)->primal_0 = vert1_c_11;
    (&_S5159)->differential_0 = _S4927;
    s_bwd_length_impl_1(&_S5159, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5160;
    (&_S5160)->primal_0 = vert0_c_11;
    (&_S5160)->differential_0 = _S4927;
    s_bwd_length_impl_1(&_S5160, 0.0f);
    float3  _S5161 = _S5158.differential_0 + _S5115;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5162;
    (&_S5162)->primal_0 = R_32;
    (&_S5162)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5163;
    (&_S5163)->primal_0 = vert2_8;
    (&_S5163)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5162, &_S5163, _S5161);
    float3  _S5164 = _S5159.differential_0 + _S5142;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5165;
    (&_S5165)->primal_0 = R_32;
    (&_S5165)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5166;
    (&_S5166)->primal_0 = vert1_8;
    (&_S5166)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5165, &_S5166, _S5164);
    float3  _S5167 = _S5160.differential_0 + _S5060 + make_float3 (_S5157.x, _S5157.y, _S5154);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5168;
    (&_S5168)->primal_0 = R_32;
    (&_S5168)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5169;
    (&_S5169)->primal_0 = vert0_8;
    (&_S5169)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5168, &_S5169, _S5167);
    float3  _S5170 = _S5163.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5171;
    (&_S5171)->primal_0 = _S4791;
    (&_S5171)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5172;
    (&_S5172)->primal_0 = _S4796;
    (&_S5172)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5171, &_S5172, _S5170);
    float _S5173 = - _S5172.differential_0.y;
    float _S5174 = _S4795 * _S5172.differential_0.x;
    float _S5175 = - (_S4787 * _S5172.differential_0.x);
    float3  _S5176 = _S5166.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5177;
    (&_S5177)->primal_0 = _S4791;
    (&_S5177)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5178;
    (&_S5178)->primal_0 = _S4794;
    (&_S5178)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5177, &_S5178, _S5176);
    float _S5179 = _S4787 * _S5178.differential_0.x;
    float _S5180 = _S4793 * _S5178.differential_0.x;
    float3  _S5181 = _S5169.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5182;
    (&_S5182)->primal_0 = _S4791;
    (&_S5182)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5183;
    (&_S5183)->primal_0 = _S4792;
    (&_S5183)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5182, &_S5183, _S5181);
    Matrix<float, 3, 3>  _S5184 = transpose_0(_S5171.differential_0 + _S5177.differential_0 + _S5182.differential_0);
    float _S5185 = 2.0f * - _S5184.rows[int(2)].z;
    float _S5186 = 2.0f * _S5184.rows[int(2)].y;
    float _S5187 = 2.0f * _S5184.rows[int(2)].x;
    float _S5188 = 2.0f * _S5184.rows[int(1)].z;
    float _S5189 = 2.0f * - _S5184.rows[int(1)].y;
    float _S5190 = 2.0f * _S5184.rows[int(1)].x;
    float _S5191 = 2.0f * _S5184.rows[int(0)].z;
    float _S5192 = 2.0f * _S5184.rows[int(0)].y;
    float _S5193 = 2.0f * - _S5184.rows[int(0)].x;
    float _S5194 = - _S5190 + _S5192;
    float _S5195 = _S5187 + - _S5191;
    float _S5196 = - _S5186 + _S5188;
    float _S5197 = _S5186 + _S5188;
    float _S5198 = _S5187 + _S5191;
    float _S5199 = _S5190 + _S5192;
    float _S5200 = quat_32.w * (_S5189 + _S5193);
    float _S5201 = quat_32.z * (_S5185 + _S5193);
    float _S5202 = quat_32.y * (_S5185 + _S5189);
    float _S5203 = quat_32.x * _S5194 + quat_32.z * _S5197 + quat_32.y * _S5198 + _S5200 + _S5200;
    float _S5204 = quat_32.x * _S5195 + quat_32.w * _S5197 + quat_32.y * _S5199 + _S5201 + _S5201;
    float _S5205 = quat_32.x * _S5196 + quat_32.w * _S5198 + quat_32.z * _S5199 + _S5202 + _S5202;
    float _S5206 = quat_32.w * _S5194 + quat_32.z * _S5195 + quat_32.y * _S5196;
    float _S5207 = _S5175 + _S5179;
    float _S5208 = 0.5f * - _S5207;
    float _S5209 = _S5173 + _S5178.differential_0.y;
    DiffPair_float_0 _S5210;
    (&_S5210)->primal_0 = _S4788;
    (&_S5210)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5210, _S5209);
    float _S5211 = _S5208 + _S5210.differential_0;
    float _S5212 = _S5174 + _S5180 + _S5183.differential_0.x;
    DiffPair_float_0 _S5213;
    (&_S5213)->primal_0 = _S4786;
    (&_S5213)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5213, _S5212);
    float _S5214 = _S5208 + _S5213.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5215;
    (&_S5215)->primal_0 = mean_c_26;
    (&_S5215)->differential_0 = _S4927;
    s_bwd_length_impl_1(&_S5215, 0.0f);
    float3  _S5216 = _S5215.differential_0 + _S4930.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5217;
    (&_S5217)->primal_0 = R_32;
    (&_S5217)->differential_0 = _S5008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5218;
    (&_S5218)->primal_0 = mean_31;
    (&_S5218)->differential_0 = _S4927;
    s_bwd_prop_mul_1(&_S5217, &_S5218, _S5216);
    float3  _S5219 = _S5161 + _S5164 + _S5167 + _S5216 + _S5011.differential_0;
    Matrix<float, 3, 3>  _S5220 = _S5162.differential_0 + _S5165.differential_0 + _S5168.differential_0 + _S5217.differential_0 + _S5012;
    float3  _S5221 = make_float3 (_S5214, _S5211, _S5207);
    float4  _S5222 = make_float4 (0.0f);
    *&((&_S5222)->w) = _S5203;
    *&((&_S5222)->z) = _S5204;
    *&((&_S5222)->y) = _S5205;
    *&((&_S5222)->x) = _S5206;
    float4  _S5223 = _S5222;
    float3  _S5224 = _S5170 + _S5176 + _S5181 + _S5218.differential_0 + _S5006;
    *v_mean_10 = _S5224;
    *v_quat_9 = _S5223;
    *v_scale_9 = _S5221;
    *v_hardness_5 = _S5088;
    (*v_sh_coeffs_8)[int(0)] = _S5066;
    (*v_sh_coeffs_8)[int(1)] = _S5067;
    (*v_sh_coeffs_8)[int(2)] = _S5068;
    (*v_sh_coeffs_8)[int(3)] = _S5069;
    (*v_sh_coeffs_8)[int(4)] = _S5070;
    (*v_sh_coeffs_8)[int(5)] = _S5071;
    (*v_sh_coeffs_8)[int(6)] = _S5072;
    (*v_sh_coeffs_8)[int(7)] = _S5073;
    (*v_sh_coeffs_8)[int(8)] = _S5074;
    (*v_sh_coeffs_8)[int(9)] = _S5075;
    (*v_sh_coeffs_8)[int(10)] = _S5076;
    (*v_sh_coeffs_8)[int(11)] = _S5077;
    (*v_sh_coeffs_8)[int(12)] = _S5078;
    (*v_sh_coeffs_8)[int(13)] = _S5079;
    (*v_sh_coeffs_8)[int(14)] = _S5080;
    (*v_sh_coeffs_8)[int(15)] = _S5081;
    (*v_ch_coeffs_3)[int(0)] = _S5063;
    (*v_ch_coeffs_3)[int(1)] = _S5064;
    *v_R_10 = _S5220;
    *v_t_9 = _S5219;
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
    bool _S5225;
    if((*u_80) >= 0.0f)
    {
        _S5225 = (*v_80) >= 0.0f;
    }
    else
    {
        _S5225 = false;
    }
    if(_S5225)
    {
        _S5225 = (*u_80 + *v_80) <= 1.0f;
    }
    else
    {
        _S5225 = false;
    }
    if(_S5225)
    {
        _S5225 = (*t_30) >= 0.0f;
    }
    else
    {
        _S5225 = false;
    }
    return _S5225;
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
    bool _S5226;
    if(u_81 >= 0.0f)
    {
        _S5226 = v_81 >= 0.0f;
    }
    else
    {
        _S5226 = false;
    }
    if(_S5226)
    {
        _S5226 = (u_81 + v_81) <= 1.0f;
    }
    else
    {
        _S5226 = false;
    }
    if(_S5226)
    {
        _S5226 = t_31 >= 0.0f;
    }
    else
    {
        _S5226 = false;
    }
    if(!_S5226)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_81), (v_81)))), ((F32_sqrt((0.5f))) * (1.0f - u_81 - v_81)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5227;
    if(opac_0 < 0.0f)
    {
        _S5227 = 0.0f;
    }
    else
    {
        _S5227 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5227;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_13)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5228 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5229 = *dpray_d_2;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_2).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5230 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S5231 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_2).primal_0);
    float _S5232 = s_primal_ctx_dot_0((*dpray_d_2).primal_0, _S5230);
    float d_30 = 1.0f / _S5232;
    float _S5233 = _S5232 * _S5232;
    float3  _S5234 = - _S5231;
    float _S5235 = s_primal_ctx_dot_0(_S5234, v2v0_2);
    float u_82 = d_30 * _S5235;
    float _S5236 = s_primal_ctx_dot_0(_S5231, v1v0_2);
    float v_82 = d_30 * _S5236;
    float3  _S5237 = - _S5230;
    float t_32 = d_30 * s_primal_ctx_dot_0(_S5237, rov0_2);
    bool _S5238;
    if(u_82 >= 0.0f)
    {
        _S5238 = v_82 >= 0.0f;
    }
    else
    {
        _S5238 = false;
    }
    if(_S5238)
    {
        _S5238 = (u_82 + v_82) <= 1.0f;
    }
    else
    {
        _S5238 = false;
    }
    if(_S5238)
    {
        _S5238 = t_32 >= 0.0f;
    }
    else
    {
        _S5238 = false;
    }
    bool _S5239 = !!_S5238;
    float _S5240;
    float _S5241;
    float _S5242;
    float _S5243;
    float _S5244;
    float _S5245;
    float _S5246;
    float _S5247;
    float _S5248;
    float _S5249;
    float _S5250;
    if(_S5239)
    {
        float _S5251 = s_primal_ctx_min_0(u_82, v_82);
        float _S5252 = s_primal_ctx_sqrt_0(0.5f);
        float _S5253 = _S5252 * (1.0f - u_82 - v_82);
        float _S5254 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5251, _S5253) * _S5254;
        float _S5255 = _S5228.primal_0.y;
        float _S5256 = 1.0f - opac_1;
        float _S5257 = 1.0f - s_primal_ctx_clamp_0(_S5255, 0.0f, 0.99989998340606689f);
        float _S5258 = 1.0f / _S5257;
        float _S5259 = _S5257 * _S5257;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5256, _S5258);
        float o_1 = _S5228.primal_0.x;
        bool _S5260 = opac_1 < 0.0f;
        if(_S5260)
        {
            _S5240 = 0.0f;
        }
        else
        {
            _S5240 = o_1 * w_1;
        }
        _S5238 = _S5260;
        _S5241 = o_1;
        _S5242 = w_1;
        _S5243 = _S5256;
        _S5244 = _S5258;
        _S5245 = _S5259;
        _S5246 = _S5255;
        _S5247 = _S5254;
        _S5248 = _S5251;
        _S5249 = _S5253;
        _S5250 = _S5252;
    }
    else
    {
        _S5238 = false;
        _S5240 = 0.0f;
        _S5241 = 0.0f;
        _S5242 = 0.0f;
        _S5243 = 0.0f;
        _S5244 = 0.0f;
        _S5245 = 0.0f;
        _S5246 = 0.0f;
        _S5247 = 0.0f;
        _S5248 = 0.0f;
        _S5249 = 0.0f;
        _S5250 = 0.0f;
    }
    float2  _S5261 = make_float2 (0.0f);
    float2  _S5262;
    if(_S5239)
    {
        if(_S5238)
        {
            _S5240 = 0.0f;
            _S5241 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5263;
            (&_S5263)->primal_0 = _S5240;
            (&_S5263)->differential_0 = 0.0f;
            DiffPair_float_0 _S5264;
            (&_S5264)->primal_0 = 0.99500000476837158f;
            (&_S5264)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5263, &_S5264, _s_dOut_13);
            float _S5265 = _S5241 * _S5263.differential_0;
            _S5240 = _S5242 * _S5263.differential_0;
            _S5241 = _S5265;
        }
        float _S5266 = - _S5241;
        DiffPair_float_0 _S5267;
        (&_S5267)->primal_0 = _S5243;
        (&_S5267)->differential_0 = 0.0f;
        DiffPair_float_0 _S5268;
        (&_S5268)->primal_0 = _S5244;
        (&_S5268)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5267, &_S5268, _S5266);
        float _S5269 = - - (_S5268.differential_0 / _S5245);
        float s_diff_opac_T_0 = - _S5267.differential_0;
        DiffPair_float_0 _S5270;
        (&_S5270)->primal_0 = _S5246;
        (&_S5270)->differential_0 = 0.0f;
        DiffPair_float_0 _S5271;
        (&_S5271)->primal_0 = 0.0f;
        (&_S5271)->differential_0 = 0.0f;
        DiffPair_float_0 _S5272;
        (&_S5272)->primal_0 = 0.99989998340606689f;
        (&_S5272)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5270, &_S5271, &_S5272, _S5269);
        float _S5273 = _S5247 * s_diff_opac_T_0;
        DiffPair_float_0 _S5274;
        (&_S5274)->primal_0 = _S5248;
        (&_S5274)->differential_0 = 0.0f;
        DiffPair_float_0 _S5275;
        (&_S5275)->primal_0 = _S5249;
        (&_S5275)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5274, &_S5275, _S5273);
        float _S5276 = - (_S5250 * _S5275.differential_0);
        DiffPair_float_0 _S5277;
        (&_S5277)->primal_0 = u_82;
        (&_S5277)->differential_0 = 0.0f;
        DiffPair_float_0 _S5278;
        (&_S5278)->primal_0 = v_82;
        (&_S5278)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5277, &_S5278, _S5274.differential_0);
        float2  _S5279 = make_float2 (_S5240, _S5270.differential_0);
        float _S5280 = _S5276 + _S5278.differential_0;
        _S5240 = _S5276 + _S5277.differential_0;
        _S5241 = _S5280;
        _S5262 = _S5279;
    }
    else
    {
        _S5240 = 0.0f;
        _S5241 = 0.0f;
        _S5262 = _S5261;
    }
    float3  _S5281 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5282;
    (&_S5282)->primal_0 = _S5237;
    (&_S5282)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5283;
    (&_S5283)->primal_0 = rov0_2;
    (&_S5283)->differential_0 = _S5281;
    s_bwd_prop_dot_0(&_S5282, &_S5283, 0.0f);
    float3  _S5284 = - _S5282.differential_0;
    float _S5285 = d_30 * _S5241;
    float _S5286 = _S5236 * _S5241;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5287;
    (&_S5287)->primal_0 = _S5231;
    (&_S5287)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5288;
    (&_S5288)->primal_0 = v1v0_2;
    (&_S5288)->differential_0 = _S5281;
    s_bwd_prop_dot_0(&_S5287, &_S5288, _S5285);
    float _S5289 = d_30 * _S5240;
    float _S5290 = _S5235 * _S5240;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5291;
    (&_S5291)->primal_0 = _S5234;
    (&_S5291)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5292;
    (&_S5292)->primal_0 = v2v0_2;
    (&_S5292)->differential_0 = _S5281;
    s_bwd_prop_dot_0(&_S5291, &_S5292, _S5289);
    float3  _S5293 = - _S5291.differential_0;
    float _S5294 = - ((_S5286 + _S5290) / _S5233);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5295;
    (&_S5295)->primal_0 = _S5229.primal_0;
    (&_S5295)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5296;
    (&_S5296)->primal_0 = _S5230;
    (&_S5296)->differential_0 = _S5281;
    s_bwd_prop_dot_0(&_S5295, &_S5296, _S5294);
    float3  _S5297 = _S5287.differential_0 + _S5293;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5298;
    (&_S5298)->primal_0 = rov0_2;
    (&_S5298)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5299;
    (&_S5299)->primal_0 = _S5229.primal_0;
    (&_S5299)->differential_0 = _S5281;
    s_bwd_prop_cross_0(&_S5298, &_S5299, _S5297);
    float3  _S5300 = _S5284 + _S5296.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5301;
    (&_S5301)->primal_0 = v1v0_2;
    (&_S5301)->differential_0 = _S5281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5302;
    (&_S5302)->primal_0 = v2v0_2;
    (&_S5302)->differential_0 = _S5281;
    s_bwd_prop_cross_0(&_S5301, &_S5302, _S5300);
    float3  _S5303 = _S5283.differential_0 + _S5298.differential_0;
    float3  _S5304 = _S5292.differential_0 + _S5302.differential_0;
    float3  _S5305 = _S5288.differential_0 + _S5301.differential_0;
    float3  _S5306 = - _S5303 + - _S5304 + - _S5305;
    float3  _S5307 = _S5295.differential_0 + _S5299.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S5307;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S5303;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5262;
    FixedArray<float3 , 3>  _S5308;
    _S5308[int(0)] = _S5281;
    _S5308[int(1)] = _S5281;
    _S5308[int(2)] = _S5281;
    _S5308[int(2)] = _S5304;
    _S5308[int(0)] = _S5306;
    _S5308[int(1)] = _S5305;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5308;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5309, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5310, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5311, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5312, float _S5313)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5309, _S5310, _S5311, _S5312, _S5313);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_6, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5314 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5315 = { _S5314, _S5314, _S5314 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_6;
    (&dp_verts_0)->differential_0 = _S5315;
    float2  _S5316 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5316;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S5314;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S5314;
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
    float3  _S5317 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S5318 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_3).primal_0);
    float _S5319 = s_primal_ctx_dot_0((*dpray_d_3).primal_0, _S5317);
    float d_32 = 1.0f / _S5319;
    float _S5320 = _S5319 * _S5319;
    float3  _S5321 = - _S5318;
    float _S5322 = s_primal_ctx_dot_0(_S5321, v2v0_4);
    float u_84 = d_32 * _S5322;
    float3  _S5323 = make_float3 (u_84);
    float _S5324 = s_primal_ctx_dot_0(_S5318, v1v0_4);
    float v_84 = d_32 * _S5324;
    float3  _S5325 = make_float3 (v_84);
    float3  _S5326 = - _S5317;
    float _S5327 = s_primal_ctx_dot_0(_S5326, rov0_4);
    float _S5328 = d_32 * _S5327;
    float3  _S5329 = make_float3 (1.0f - u_84 - v_84);
    DiffPair_float_0 _S5330;
    (&_S5330)->primal_0 = s_primal_ctx_max_0(_S5328, 9.999999960041972e-13f);
    (&_S5330)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5330, dpdepth_2);
    DiffPair_float_0 _S5331;
    (&_S5331)->primal_0 = _S5328;
    (&_S5331)->differential_0 = 0.0f;
    DiffPair_float_0 _S5332;
    (&_S5332)->primal_0 = 9.999999960041972e-13f;
    (&_S5332)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5331, &_S5332, _S5330.differential_0);
    float3  _S5333 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5334 = _S5325 * dpcolor_1;
    float3  _S5335 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5336 = _S5323 * dpcolor_1;
    float3  _S5337 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5338 = _S5329 * dpcolor_1;
    float _S5339 = - (_S5337.x + _S5337.y + _S5337.z);
    float _S5340 = d_32 * _S5331.differential_0;
    float _S5341 = _S5327 * _S5331.differential_0;
    float3  _S5342 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5343;
    (&_S5343)->primal_0 = _S5326;
    (&_S5343)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5344;
    (&_S5344)->primal_0 = rov0_4;
    (&_S5344)->differential_0 = _S5342;
    s_bwd_prop_dot_0(&_S5343, &_S5344, _S5340);
    float3  _S5345 = - _S5343.differential_0;
    float _S5346 = _S5339 + _S5333.x + _S5333.y + _S5333.z;
    float _S5347 = d_32 * _S5346;
    float _S5348 = _S5324 * _S5346;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5349;
    (&_S5349)->primal_0 = _S5318;
    (&_S5349)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5350;
    (&_S5350)->primal_0 = v1v0_4;
    (&_S5350)->differential_0 = _S5342;
    s_bwd_prop_dot_0(&_S5349, &_S5350, _S5347);
    float _S5351 = _S5339 + _S5335.x + _S5335.y + _S5335.z;
    float _S5352 = d_32 * _S5351;
    float _S5353 = _S5322 * _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5354;
    (&_S5354)->primal_0 = _S5321;
    (&_S5354)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5355;
    (&_S5355)->primal_0 = v2v0_4;
    (&_S5355)->differential_0 = _S5342;
    s_bwd_prop_dot_0(&_S5354, &_S5355, _S5352);
    float3  _S5356 = - _S5354.differential_0;
    float _S5357 = - ((_S5341 + _S5348 + _S5353) / _S5320);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5358;
    (&_S5358)->primal_0 = (*dpray_d_3).primal_0;
    (&_S5358)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5359;
    (&_S5359)->primal_0 = _S5317;
    (&_S5359)->differential_0 = _S5342;
    s_bwd_prop_dot_0(&_S5358, &_S5359, _S5357);
    float3  _S5360 = _S5349.differential_0 + _S5356;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5361;
    (&_S5361)->primal_0 = rov0_4;
    (&_S5361)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5362;
    (&_S5362)->primal_0 = (*dpray_d_3).primal_0;
    (&_S5362)->differential_0 = _S5342;
    s_bwd_prop_cross_0(&_S5361, &_S5362, _S5360);
    float3  _S5363 = _S5345 + _S5359.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5364;
    (&_S5364)->primal_0 = v1v0_4;
    (&_S5364)->differential_0 = _S5342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5365;
    (&_S5365)->primal_0 = v2v0_4;
    (&_S5365)->differential_0 = _S5342;
    s_bwd_prop_cross_0(&_S5364, &_S5365, _S5363);
    float3  _S5366 = _S5344.differential_0 + _S5361.differential_0;
    float3  _S5367 = _S5355.differential_0 + _S5365.differential_0;
    float3  _S5368 = _S5350.differential_0 + _S5364.differential_0;
    float3  _S5369 = - _S5366 + - _S5367 + - _S5368;
    float3  _S5370 = _S5358.differential_0 + _S5362.differential_0;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S5370;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S5366;
    FixedArray<float3 , 3>  _S5371;
    _S5371[int(0)] = _S5342;
    _S5371[int(1)] = _S5342;
    _S5371[int(2)] = _S5342;
    _S5371[int(2)] = _S5334;
    _S5371[int(1)] = _S5336;
    _S5371[int(0)] = _S5338;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5371;
    FixedArray<float3 , 3>  _S5372;
    _S5372[int(0)] = _S5342;
    _S5372[int(1)] = _S5342;
    _S5372[int(2)] = _S5342;
    _S5372[int(2)] = _S5367;
    _S5372[int(0)] = _S5369;
    _S5372[int(1)] = _S5368;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5372;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5373, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5374, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5375, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5376, float3  _S5377, float _S5378)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5373, _S5374, _S5375, _S5376, _S5377, _S5378);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_9, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_9, float3  ray_d_9, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5379 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5380 = { _S5379, _S5379, _S5379 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_9;
    (&dp_verts_1)->differential_0 = _S5380;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5380;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_9;
    (&dp_ray_o_3)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_9;
    (&dp_ray_d_3)->differential_0 = _S5379;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  * densities_0, FixedArray<float3 , 16>  * sh_coeffs_27, Matrix<float, 3, 3>  R_33, float3  t_33, float fx_36, float fy_36, float cx_31, float cy_31, FixedArray<float, 10>  * dist_coeffs_40, uint image_width_27, uint image_height_27, float near_plane_18, float far_plane_18, int4  * aabb_xyxy_18, float * depth_21, float3  * rgbs_7)
{
    float2  * _S5381;
    float2  * _S5382;
    float2  * _S5383;
    float2  * _S5384;
    float2  * _S5385;
    float2  * _S5386;
    float2  * _S5387;
    float2  * _S5388;
    bool _S5389;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_0;
        float3  _S5390 = mul_0(R_33, pos_0) + t_33;
        pos_c_0[int(0)] = _S5390;
        float _S5391 = _S5390.z;
        float _S5392 = (F32_min((far_plane_18), (_S5391)));
        float _S5393 = (F32_max((near_plane_18), (_S5391)));
        float3  _S5394 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 0.0f)) + t_33;
        pos_c_0[int(1)] = _S5394;
        float _S5395 = _S5394.z;
        float _S5396 = (F32_min((_S5392), (_S5395)));
        float _S5397 = (F32_max((_S5393), (_S5395)));
        float3  _S5398 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(2)] = _S5398;
        float _S5399 = _S5398.z;
        float _S5400 = (F32_min((_S5396), (_S5399)));
        float _S5401 = (F32_max((_S5397), (_S5399)));
        float3  _S5402 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(3)] = _S5402;
        float _S5403 = _S5402.z;
        float _S5404 = (F32_min((_S5400), (_S5403)));
        float _S5405 = (F32_max((_S5401), (_S5403)));
        float3  _S5406 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(4)] = _S5406;
        float _S5407 = _S5406.z;
        float _S5408 = (F32_min((_S5404), (_S5407)));
        float _S5409 = (F32_max((_S5405), (_S5407)));
        float3  _S5410 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(5)] = _S5410;
        float _S5411 = _S5410.z;
        float _S5412 = (F32_min((_S5408), (_S5411)));
        float _S5413 = (F32_max((_S5409), (_S5411)));
        float3  _S5414 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 1.0f)) + t_33;
        pos_c_0[int(6)] = _S5414;
        float _S5415 = _S5414.z;
        float _S5416 = (F32_min((_S5412), (_S5415)));
        float _S5417 = (F32_max((_S5413), (_S5415)));
        float3  _S5418 = mul_0(R_33, pos_0 + make_float3 (size_0)) + t_33;
        pos_c_0[int(7)] = _S5418;
        float _S5419 = _S5418.z;
        float _S5420 = (F32_min((_S5416), (_S5419)));
        float _S5421 = (F32_max((_S5417), (_S5419)));
        bool _S5422;
        if(_S5420 < near_plane_18)
        {
            _S5422 = true;
        }
        else
        {
            _S5422 = _S5421 > far_plane_18;
        }
        if(_S5422)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_27 = mul_0(R_33, pos_0 + make_float3 (0.5f * size_0)) + t_33;
        FixedArray<float2 , 8>  uv_7;
        for(;;)
        {
            float3  _S5423 = pos_c_0[int(0)];
            _S5381 = &uv_7[int(0)];
            for(;;)
            {
                float _S5424 = _S5423.z;
                uv_7[int(0)] = float2 {_S5423.x, _S5423.y} / make_float2 (_S5424);
                if(_S5424 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5425 = camera_distortion_jac_0(uv_7[int(0)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5425)), ((F32_min((_S5425.rows[int(0)].x), (_S5425.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_85 = uv_7[int(0)].x;
                float v_85 = uv_7[int(0)].y;
                float r2_80 = u_85 * u_85 + v_85 * v_85;
                float2  _S5426 = uv_7[int(0)] * make_float2 (1.0f + r2_80 * ((*dist_coeffs_40)[int(0)] + r2_80 * ((*dist_coeffs_40)[int(1)] + r2_80 * ((*dist_coeffs_40)[int(2)] + r2_80 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_85 * v_85 + (*dist_coeffs_40)[int(5)] * (r2_80 + 2.0f * u_85 * u_85) + (*dist_coeffs_40)[int(6)] * r2_80, 2.0f * (*dist_coeffs_40)[int(5)] * u_85 * v_85 + (*dist_coeffs_40)[int(4)] * (r2_80 + 2.0f * v_85 * v_85) + (*dist_coeffs_40)[int(7)] * r2_80);
                float2  _S5427 = _S5426 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5426.x + (*dist_coeffs_40)[int(9)] * _S5426.y, 0.0f);
                uv_7[int(0)] = make_float2 (fx_36 * _S5427.x + cx_31, fy_36 * _S5427.y + cy_31);
                break;
            }
            bool all_valid_16 = true & (!_S5422);
            float3  _S5428 = pos_c_0[int(1)];
            _S5382 = &uv_7[int(1)];
            for(;;)
            {
                float _S5429 = _S5428.z;
                uv_7[int(1)] = float2 {_S5428.x, _S5428.y} / make_float2 (_S5429);
                if(_S5429 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5430 = camera_distortion_jac_0(uv_7[int(1)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5430)), ((F32_min((_S5430.rows[int(0)].x), (_S5430.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_86 = uv_7[int(1)].x;
                float v_86 = uv_7[int(1)].y;
                float r2_81 = u_86 * u_86 + v_86 * v_86;
                float2  _S5431 = uv_7[int(1)] * make_float2 (1.0f + r2_81 * ((*dist_coeffs_40)[int(0)] + r2_81 * ((*dist_coeffs_40)[int(1)] + r2_81 * ((*dist_coeffs_40)[int(2)] + r2_81 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_86 * v_86 + (*dist_coeffs_40)[int(5)] * (r2_81 + 2.0f * u_86 * u_86) + (*dist_coeffs_40)[int(6)] * r2_81, 2.0f * (*dist_coeffs_40)[int(5)] * u_86 * v_86 + (*dist_coeffs_40)[int(4)] * (r2_81 + 2.0f * v_86 * v_86) + (*dist_coeffs_40)[int(7)] * r2_81);
                float2  _S5432 = _S5431 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5431.x + (*dist_coeffs_40)[int(9)] * _S5431.y, 0.0f);
                uv_7[int(1)] = make_float2 (fx_36 * _S5432.x + cx_31, fy_36 * _S5432.y + cy_31);
                break;
            }
            bool all_valid_17 = all_valid_16 & (!_S5422);
            float3  _S5433 = pos_c_0[int(2)];
            _S5383 = &uv_7[int(2)];
            for(;;)
            {
                float _S5434 = _S5433.z;
                uv_7[int(2)] = float2 {_S5433.x, _S5433.y} / make_float2 (_S5434);
                if(_S5434 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5435 = camera_distortion_jac_0(uv_7[int(2)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5435)), ((F32_min((_S5435.rows[int(0)].x), (_S5435.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_87 = uv_7[int(2)].x;
                float v_87 = uv_7[int(2)].y;
                float r2_82 = u_87 * u_87 + v_87 * v_87;
                float2  _S5436 = uv_7[int(2)] * make_float2 (1.0f + r2_82 * ((*dist_coeffs_40)[int(0)] + r2_82 * ((*dist_coeffs_40)[int(1)] + r2_82 * ((*dist_coeffs_40)[int(2)] + r2_82 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_87 * v_87 + (*dist_coeffs_40)[int(5)] * (r2_82 + 2.0f * u_87 * u_87) + (*dist_coeffs_40)[int(6)] * r2_82, 2.0f * (*dist_coeffs_40)[int(5)] * u_87 * v_87 + (*dist_coeffs_40)[int(4)] * (r2_82 + 2.0f * v_87 * v_87) + (*dist_coeffs_40)[int(7)] * r2_82);
                float2  _S5437 = _S5436 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5436.x + (*dist_coeffs_40)[int(9)] * _S5436.y, 0.0f);
                uv_7[int(2)] = make_float2 (fx_36 * _S5437.x + cx_31, fy_36 * _S5437.y + cy_31);
                break;
            }
            bool all_valid_18 = all_valid_17 & (!_S5422);
            float3  _S5438 = pos_c_0[int(3)];
            _S5384 = &uv_7[int(3)];
            for(;;)
            {
                float _S5439 = _S5438.z;
                uv_7[int(3)] = float2 {_S5438.x, _S5438.y} / make_float2 (_S5439);
                if(_S5439 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5440 = camera_distortion_jac_0(uv_7[int(3)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5440)), ((F32_min((_S5440.rows[int(0)].x), (_S5440.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_88 = uv_7[int(3)].x;
                float v_88 = uv_7[int(3)].y;
                float r2_83 = u_88 * u_88 + v_88 * v_88;
                float2  _S5441 = uv_7[int(3)] * make_float2 (1.0f + r2_83 * ((*dist_coeffs_40)[int(0)] + r2_83 * ((*dist_coeffs_40)[int(1)] + r2_83 * ((*dist_coeffs_40)[int(2)] + r2_83 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_88 * v_88 + (*dist_coeffs_40)[int(5)] * (r2_83 + 2.0f * u_88 * u_88) + (*dist_coeffs_40)[int(6)] * r2_83, 2.0f * (*dist_coeffs_40)[int(5)] * u_88 * v_88 + (*dist_coeffs_40)[int(4)] * (r2_83 + 2.0f * v_88 * v_88) + (*dist_coeffs_40)[int(7)] * r2_83);
                float2  _S5442 = _S5441 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5441.x + (*dist_coeffs_40)[int(9)] * _S5441.y, 0.0f);
                uv_7[int(3)] = make_float2 (fx_36 * _S5442.x + cx_31, fy_36 * _S5442.y + cy_31);
                break;
            }
            bool all_valid_19 = all_valid_18 & (!_S5422);
            float3  _S5443 = pos_c_0[int(4)];
            _S5385 = &uv_7[int(4)];
            for(;;)
            {
                float _S5444 = _S5443.z;
                uv_7[int(4)] = float2 {_S5443.x, _S5443.y} / make_float2 (_S5444);
                if(_S5444 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5445 = camera_distortion_jac_0(uv_7[int(4)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5445)), ((F32_min((_S5445.rows[int(0)].x), (_S5445.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_89 = uv_7[int(4)].x;
                float v_89 = uv_7[int(4)].y;
                float r2_84 = u_89 * u_89 + v_89 * v_89;
                float2  _S5446 = uv_7[int(4)] * make_float2 (1.0f + r2_84 * ((*dist_coeffs_40)[int(0)] + r2_84 * ((*dist_coeffs_40)[int(1)] + r2_84 * ((*dist_coeffs_40)[int(2)] + r2_84 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_89 * v_89 + (*dist_coeffs_40)[int(5)] * (r2_84 + 2.0f * u_89 * u_89) + (*dist_coeffs_40)[int(6)] * r2_84, 2.0f * (*dist_coeffs_40)[int(5)] * u_89 * v_89 + (*dist_coeffs_40)[int(4)] * (r2_84 + 2.0f * v_89 * v_89) + (*dist_coeffs_40)[int(7)] * r2_84);
                float2  _S5447 = _S5446 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5446.x + (*dist_coeffs_40)[int(9)] * _S5446.y, 0.0f);
                uv_7[int(4)] = make_float2 (fx_36 * _S5447.x + cx_31, fy_36 * _S5447.y + cy_31);
                break;
            }
            bool all_valid_20 = all_valid_19 & (!_S5422);
            float3  _S5448 = pos_c_0[int(5)];
            _S5386 = &uv_7[int(5)];
            for(;;)
            {
                float _S5449 = _S5448.z;
                uv_7[int(5)] = float2 {_S5448.x, _S5448.y} / make_float2 (_S5449);
                if(_S5449 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5450 = camera_distortion_jac_0(uv_7[int(5)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5450)), ((F32_min((_S5450.rows[int(0)].x), (_S5450.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_90 = uv_7[int(5)].x;
                float v_90 = uv_7[int(5)].y;
                float r2_85 = u_90 * u_90 + v_90 * v_90;
                float2  _S5451 = uv_7[int(5)] * make_float2 (1.0f + r2_85 * ((*dist_coeffs_40)[int(0)] + r2_85 * ((*dist_coeffs_40)[int(1)] + r2_85 * ((*dist_coeffs_40)[int(2)] + r2_85 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_90 * v_90 + (*dist_coeffs_40)[int(5)] * (r2_85 + 2.0f * u_90 * u_90) + (*dist_coeffs_40)[int(6)] * r2_85, 2.0f * (*dist_coeffs_40)[int(5)] * u_90 * v_90 + (*dist_coeffs_40)[int(4)] * (r2_85 + 2.0f * v_90 * v_90) + (*dist_coeffs_40)[int(7)] * r2_85);
                float2  _S5452 = _S5451 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5451.x + (*dist_coeffs_40)[int(9)] * _S5451.y, 0.0f);
                uv_7[int(5)] = make_float2 (fx_36 * _S5452.x + cx_31, fy_36 * _S5452.y + cy_31);
                break;
            }
            bool all_valid_21 = all_valid_20 & (!_S5422);
            float3  _S5453 = pos_c_0[int(6)];
            _S5387 = &uv_7[int(6)];
            for(;;)
            {
                float _S5454 = _S5453.z;
                uv_7[int(6)] = float2 {_S5453.x, _S5453.y} / make_float2 (_S5454);
                if(_S5454 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5455 = camera_distortion_jac_0(uv_7[int(6)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5455)), ((F32_min((_S5455.rows[int(0)].x), (_S5455.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_91 = uv_7[int(6)].x;
                float v_91 = uv_7[int(6)].y;
                float r2_86 = u_91 * u_91 + v_91 * v_91;
                float2  _S5456 = uv_7[int(6)] * make_float2 (1.0f + r2_86 * ((*dist_coeffs_40)[int(0)] + r2_86 * ((*dist_coeffs_40)[int(1)] + r2_86 * ((*dist_coeffs_40)[int(2)] + r2_86 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_91 * v_91 + (*dist_coeffs_40)[int(5)] * (r2_86 + 2.0f * u_91 * u_91) + (*dist_coeffs_40)[int(6)] * r2_86, 2.0f * (*dist_coeffs_40)[int(5)] * u_91 * v_91 + (*dist_coeffs_40)[int(4)] * (r2_86 + 2.0f * v_91 * v_91) + (*dist_coeffs_40)[int(7)] * r2_86);
                float2  _S5457 = _S5456 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5456.x + (*dist_coeffs_40)[int(9)] * _S5456.y, 0.0f);
                uv_7[int(6)] = make_float2 (fx_36 * _S5457.x + cx_31, fy_36 * _S5457.y + cy_31);
                break;
            }
            bool all_valid_22 = all_valid_21 & (!_S5422);
            float3  _S5458 = pos_c_0[int(7)];
            _S5388 = &uv_7[int(7)];
            for(;;)
            {
                float _S5459 = _S5458.z;
                uv_7[int(7)] = float2 {_S5458.x, _S5458.y} / make_float2 (_S5459);
                if(_S5459 < 0.0f)
                {
                    _S5422 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5460 = camera_distortion_jac_0(uv_7[int(7)], dist_coeffs_40);
                    _S5422 = !((F32_min((determinant_0(_S5460)), ((F32_min((_S5460.rows[int(0)].x), (_S5460.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5422)
                {
                    break;
                }
                float u_92 = uv_7[int(7)].x;
                float v_92 = uv_7[int(7)].y;
                float r2_87 = u_92 * u_92 + v_92 * v_92;
                float2  _S5461 = uv_7[int(7)] * make_float2 (1.0f + r2_87 * ((*dist_coeffs_40)[int(0)] + r2_87 * ((*dist_coeffs_40)[int(1)] + r2_87 * ((*dist_coeffs_40)[int(2)] + r2_87 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_92 * v_92 + (*dist_coeffs_40)[int(5)] * (r2_87 + 2.0f * u_92 * u_92) + (*dist_coeffs_40)[int(6)] * r2_87, 2.0f * (*dist_coeffs_40)[int(5)] * u_92 * v_92 + (*dist_coeffs_40)[int(4)] * (r2_87 + 2.0f * v_92 * v_92) + (*dist_coeffs_40)[int(7)] * r2_87);
                float2  _S5462 = _S5461 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5461.x + (*dist_coeffs_40)[int(9)] * _S5461.y, 0.0f);
                uv_7[int(7)] = make_float2 (fx_36 * _S5462.x + cx_31, fy_36 * _S5462.y + cy_31);
                break;
            }
            _S5389 = all_valid_22 & (!_S5422);
            break;
        }
        if(!_S5389)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*densities_0)[int(0)]), ((*densities_0)[int(1)])))), ((*densities_0)[int(2)])))), ((*densities_0)[int(3)])))), ((*densities_0)[int(4)])))), ((*densities_0)[int(5)])))), ((*densities_0)[int(6)])))), ((*densities_0)[int(7)]))) * size_0 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5463 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5381).x), ((*_S5382).x)))), ((*_S5383).x)))), ((*_S5384).x)))), ((*_S5385).x)))), ((*_S5386).x)))), ((*_S5387).x)))), ((*_S5388).x)));
        float _S5464 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5381).x), ((*_S5382).x)))), ((*_S5383).x)))), ((*_S5384).x)))), ((*_S5385).x)))), ((*_S5386).x)))), ((*_S5387).x)))), ((*_S5388).x)));
        float _S5465 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5381).y), ((*_S5382).y)))), ((*_S5383).y)))), ((*_S5384).y)))), ((*_S5385).y)))), ((*_S5386).y)))), ((*_S5387).y)))), ((*_S5388).y)));
        float _S5466 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5381).y), ((*_S5382).y)))), ((*_S5383).y)))), ((*_S5384).y)))), ((*_S5385).y)))), ((*_S5386).y)))), ((*_S5387).y)))), ((*_S5388).y)));
        if(_S5463 <= 0.0f)
        {
            _S5422 = true;
        }
        else
        {
            _S5422 = _S5464 >= float(image_width_27);
        }
        if(_S5422)
        {
            _S5422 = true;
        }
        else
        {
            _S5422 = _S5465 <= 0.0f;
        }
        if(_S5422)
        {
            _S5422 = true;
        }
        else
        {
            _S5422 = _S5466 >= float(image_height_27);
        }
        if(_S5422)
        {
            _S5422 = true;
        }
        else
        {
            if(_S5420 <= 0.0f)
            {
                if(_S5464 <= 0.0f)
                {
                    _S5422 = _S5463 >= float(image_width_27);
                }
                else
                {
                    _S5422 = false;
                }
                if(_S5422)
                {
                    _S5422 = true;
                }
                else
                {
                    if(_S5466 <= 0.0f)
                    {
                        _S5422 = _S5465 >= float(image_width_27);
                    }
                    else
                    {
                        _S5422 = false;
                    }
                }
            }
            else
            {
                _S5422 = false;
            }
        }
        if(_S5422)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_18 = make_int4 (int((F32_floor((_S5464)))), int((F32_floor((_S5466)))), int((F32_ceil((_S5463)))), int((F32_ceil((_S5465)))));
        *depth_21 = 0.5f * (F32_log((dot_0(mean_c_27, mean_c_27) + 9.99999997475242708e-07f)));
        float3  _S5467 = mean_c_27 - - mul_0(transpose_0(R_33), t_33);
        float3  _S5468 = make_float3 (0.282094806432724f) * (*sh_coeffs_27)[int(0)];
        *rgbs_7 = _S5468;
        float _S5469 = _S5467.x;
        float _S5470 = _S5467.y;
        float _S5471 = _S5467.z;
        float norm_18 = (F32_sqrt((_S5469 * _S5469 + _S5470 * _S5470 + _S5471 * _S5471)));
        float x_63 = _S5469 / norm_18;
        float y_30 = _S5470 / norm_18;
        float z_27 = _S5471 / norm_18;
        float3  _S5472 = _S5468 + make_float3 (0.48860251903533936f) * (make_float3 (- y_30) * (*sh_coeffs_27)[int(1)] + make_float3 (z_27) * (*sh_coeffs_27)[int(2)] - make_float3 (x_63) * (*sh_coeffs_27)[int(3)]);
        *rgbs_7 = _S5472;
        float z2_60 = z_27 * z_27;
        float fTmp0B_27 = -1.09254848957061768f * z_27;
        float fC1_27 = x_63 * x_63 - y_30 * y_30;
        float fS1_27 = 2.0f * x_63 * y_30;
        float3  _S5473 = _S5472 + (make_float3 (0.54627424478530884f * fS1_27) * (*sh_coeffs_27)[int(4)] + make_float3 (fTmp0B_27 * y_30) * (*sh_coeffs_27)[int(5)] + make_float3 (0.94617468118667603f * z2_60 - 0.31539157032966614f) * (*sh_coeffs_27)[int(6)] + make_float3 (fTmp0B_27 * x_63) * (*sh_coeffs_27)[int(7)] + make_float3 (0.54627424478530884f * fC1_27) * (*sh_coeffs_27)[int(8)]);
        *rgbs_7 = _S5473;
        float fTmp0C_27 = -2.28522896766662598f * z2_60 + 0.4570457935333252f;
        float fTmp1B_27 = 1.44530570507049561f * z_27;
        *rgbs_7 = max_0(_S5473 + (make_float3 (-0.59004360437393188f * (x_63 * fS1_27 + y_30 * fC1_27)) * (*sh_coeffs_27)[int(9)] + make_float3 (fTmp1B_27 * fS1_27) * (*sh_coeffs_27)[int(10)] + make_float3 (fTmp0C_27 * y_30) * (*sh_coeffs_27)[int(11)] + make_float3 (z_27 * (1.86588168144226074f * z2_60 - 1.11952900886535645f)) * (*sh_coeffs_27)[int(12)] + make_float3 (fTmp0C_27 * x_63) * (*sh_coeffs_27)[int(13)] + make_float3 (fTmp1B_27 * fC1_27) * (*sh_coeffs_27)[int(14)] + make_float3 (-0.59004360437393188f * (x_63 * fC1_27 - y_30 * fS1_27)) * (*sh_coeffs_27)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  * densities_1, FixedArray<float3 , 16>  * sh_coeffs_28, Matrix<float, 3, 3>  R_34, float3  t_34, float fx_37, float fy_37, float cx_32, float cy_32, FixedArray<float, 10>  * dist_coeffs_41, uint image_width_28, uint image_height_28, float near_plane_19, float far_plane_19, int4  * aabb_xyxy_19, float * depth_22, float3  * rgbs_8)
{
    float2  * _S5474;
    bool _S5475;
    float2  * _S5476;
    bool _S5477;
    float2  * _S5478;
    bool _S5479;
    float2  * _S5480;
    bool _S5481;
    float2  * _S5482;
    bool _S5483;
    float2  * _S5484;
    bool _S5485;
    float2  * _S5486;
    bool _S5487;
    float2  * _S5488;
    bool _S5489;
    bool _S5490;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S5491 = mul_0(R_34, pos_1) + t_34;
        pos_c_1[int(0)] = _S5491;
        float _S5492 = length_1(_S5491);
        float _S5493 = (F32_min((far_plane_19), (_S5492)));
        float _S5494 = (F32_max((near_plane_19), (_S5492)));
        float3  _S5495 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_34;
        pos_c_1[int(1)] = _S5495;
        float _S5496 = length_1(_S5495);
        float _S5497 = (F32_min((_S5493), (_S5496)));
        float _S5498 = (F32_max((_S5494), (_S5496)));
        float3  _S5499 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(2)] = _S5499;
        float _S5500 = length_1(_S5499);
        float _S5501 = (F32_min((_S5497), (_S5500)));
        float _S5502 = (F32_max((_S5498), (_S5500)));
        float3  _S5503 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(3)] = _S5503;
        float _S5504 = length_1(_S5503);
        float _S5505 = (F32_min((_S5501), (_S5504)));
        float _S5506 = (F32_max((_S5502), (_S5504)));
        float3  _S5507 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(4)] = _S5507;
        float _S5508 = length_1(_S5507);
        float _S5509 = (F32_min((_S5505), (_S5508)));
        float _S5510 = (F32_max((_S5506), (_S5508)));
        float3  _S5511 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(5)] = _S5511;
        float _S5512 = length_1(_S5511);
        float _S5513 = (F32_min((_S5509), (_S5512)));
        float _S5514 = (F32_max((_S5510), (_S5512)));
        float3  _S5515 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_34;
        pos_c_1[int(6)] = _S5515;
        float _S5516 = length_1(_S5515);
        float _S5517 = (F32_min((_S5513), (_S5516)));
        float _S5518 = (F32_max((_S5514), (_S5516)));
        float3  _S5519 = mul_0(R_34, pos_1 + make_float3 (size_1)) + t_34;
        pos_c_1[int(7)] = _S5519;
        float _S5520 = length_1(_S5519);
        float _S5521 = (F32_min((_S5517), (_S5520)));
        float _S5522 = (F32_max((_S5518), (_S5520)));
        bool _S5523;
        if(_S5521 < near_plane_19)
        {
            _S5523 = true;
        }
        else
        {
            _S5523 = _S5522 > far_plane_19;
        }
        if(_S5523)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_28 = mul_0(R_34, pos_1 + make_float3 (0.5f * size_1)) + t_34;
        FixedArray<float2 , 8>  uv_8;
        for(;;)
        {
            float k_13;
            float3  _S5524 = pos_c_1[int(0)];
            _S5474 = &uv_8[int(0)];
            for(;;)
            {
                float2  _S5525 = float2 {_S5524.x, _S5524.y};
                float r_36 = length_0(_S5525);
                float _S5526 = _S5524.z;
                float theta_32 = (F32_atan2((r_36), (_S5526)));
                if(theta_32 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_32 * theta_32 / 3.0f) / _S5526;
                }
                else
                {
                    k_13 = theta_32 / r_36;
                }
                float2  _S5527 = _S5525 * make_float2 (k_13);
                uv_8[int(0)] = _S5527;
                Matrix<float, 2, 2>  _S5528 = camera_distortion_jac_0(_S5527, dist_coeffs_41);
                bool _S5529 = !((F32_min((determinant_0(_S5528)), ((F32_min((_S5528.rows[int(0)].x), (_S5528.rows[int(1)].y)))))) > 0.0f);
                _S5475 = _S5529;
                if(_S5529)
                {
                    break;
                }
                float u_93 = uv_8[int(0)].x;
                float v_93 = uv_8[int(0)].y;
                float r2_88 = u_93 * u_93 + v_93 * v_93;
                float2  _S5530 = uv_8[int(0)] * make_float2 (1.0f + r2_88 * ((*dist_coeffs_41)[int(0)] + r2_88 * ((*dist_coeffs_41)[int(1)] + r2_88 * ((*dist_coeffs_41)[int(2)] + r2_88 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_93 * v_93 + (*dist_coeffs_41)[int(5)] * (r2_88 + 2.0f * u_93 * u_93) + (*dist_coeffs_41)[int(6)] * r2_88, 2.0f * (*dist_coeffs_41)[int(5)] * u_93 * v_93 + (*dist_coeffs_41)[int(4)] * (r2_88 + 2.0f * v_93 * v_93) + (*dist_coeffs_41)[int(7)] * r2_88);
                float2  _S5531 = _S5530 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5530.x + (*dist_coeffs_41)[int(9)] * _S5530.y, 0.0f);
                uv_8[int(0)] = make_float2 (fx_37 * _S5531.x + cx_32, fy_37 * _S5531.y + cy_32);
                break;
            }
            bool all_valid_23 = true & (!_S5475);
            float3  _S5532 = pos_c_1[int(1)];
            _S5476 = &uv_8[int(1)];
            for(;;)
            {
                float2  _S5533 = float2 {_S5532.x, _S5532.y};
                float r_37 = length_0(_S5533);
                float _S5534 = _S5532.z;
                float theta_33 = (F32_atan2((r_37), (_S5534)));
                if(theta_33 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_33 * theta_33 / 3.0f) / _S5534;
                }
                else
                {
                    k_13 = theta_33 / r_37;
                }
                float2  _S5535 = _S5533 * make_float2 (k_13);
                uv_8[int(1)] = _S5535;
                Matrix<float, 2, 2>  _S5536 = camera_distortion_jac_0(_S5535, dist_coeffs_41);
                bool _S5537 = !((F32_min((determinant_0(_S5536)), ((F32_min((_S5536.rows[int(0)].x), (_S5536.rows[int(1)].y)))))) > 0.0f);
                _S5477 = _S5537;
                if(_S5537)
                {
                    break;
                }
                float u_94 = uv_8[int(1)].x;
                float v_94 = uv_8[int(1)].y;
                float r2_89 = u_94 * u_94 + v_94 * v_94;
                float2  _S5538 = uv_8[int(1)] * make_float2 (1.0f + r2_89 * ((*dist_coeffs_41)[int(0)] + r2_89 * ((*dist_coeffs_41)[int(1)] + r2_89 * ((*dist_coeffs_41)[int(2)] + r2_89 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_94 * v_94 + (*dist_coeffs_41)[int(5)] * (r2_89 + 2.0f * u_94 * u_94) + (*dist_coeffs_41)[int(6)] * r2_89, 2.0f * (*dist_coeffs_41)[int(5)] * u_94 * v_94 + (*dist_coeffs_41)[int(4)] * (r2_89 + 2.0f * v_94 * v_94) + (*dist_coeffs_41)[int(7)] * r2_89);
                float2  _S5539 = _S5538 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5538.x + (*dist_coeffs_41)[int(9)] * _S5538.y, 0.0f);
                uv_8[int(1)] = make_float2 (fx_37 * _S5539.x + cx_32, fy_37 * _S5539.y + cy_32);
                break;
            }
            bool all_valid_24 = all_valid_23 & (!_S5477);
            float3  _S5540 = pos_c_1[int(2)];
            _S5478 = &uv_8[int(2)];
            for(;;)
            {
                float2  _S5541 = float2 {_S5540.x, _S5540.y};
                float r_38 = length_0(_S5541);
                float _S5542 = _S5540.z;
                float theta_34 = (F32_atan2((r_38), (_S5542)));
                if(theta_34 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_34 * theta_34 / 3.0f) / _S5542;
                }
                else
                {
                    k_13 = theta_34 / r_38;
                }
                float2  _S5543 = _S5541 * make_float2 (k_13);
                uv_8[int(2)] = _S5543;
                Matrix<float, 2, 2>  _S5544 = camera_distortion_jac_0(_S5543, dist_coeffs_41);
                bool _S5545 = !((F32_min((determinant_0(_S5544)), ((F32_min((_S5544.rows[int(0)].x), (_S5544.rows[int(1)].y)))))) > 0.0f);
                _S5479 = _S5545;
                if(_S5545)
                {
                    break;
                }
                float u_95 = uv_8[int(2)].x;
                float v_95 = uv_8[int(2)].y;
                float r2_90 = u_95 * u_95 + v_95 * v_95;
                float2  _S5546 = uv_8[int(2)] * make_float2 (1.0f + r2_90 * ((*dist_coeffs_41)[int(0)] + r2_90 * ((*dist_coeffs_41)[int(1)] + r2_90 * ((*dist_coeffs_41)[int(2)] + r2_90 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_95 * v_95 + (*dist_coeffs_41)[int(5)] * (r2_90 + 2.0f * u_95 * u_95) + (*dist_coeffs_41)[int(6)] * r2_90, 2.0f * (*dist_coeffs_41)[int(5)] * u_95 * v_95 + (*dist_coeffs_41)[int(4)] * (r2_90 + 2.0f * v_95 * v_95) + (*dist_coeffs_41)[int(7)] * r2_90);
                float2  _S5547 = _S5546 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5546.x + (*dist_coeffs_41)[int(9)] * _S5546.y, 0.0f);
                uv_8[int(2)] = make_float2 (fx_37 * _S5547.x + cx_32, fy_37 * _S5547.y + cy_32);
                break;
            }
            bool all_valid_25 = all_valid_24 & (!_S5479);
            float3  _S5548 = pos_c_1[int(3)];
            _S5480 = &uv_8[int(3)];
            for(;;)
            {
                float2  _S5549 = float2 {_S5548.x, _S5548.y};
                float r_39 = length_0(_S5549);
                float _S5550 = _S5548.z;
                float theta_35 = (F32_atan2((r_39), (_S5550)));
                if(theta_35 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_35 * theta_35 / 3.0f) / _S5550;
                }
                else
                {
                    k_13 = theta_35 / r_39;
                }
                float2  _S5551 = _S5549 * make_float2 (k_13);
                uv_8[int(3)] = _S5551;
                Matrix<float, 2, 2>  _S5552 = camera_distortion_jac_0(_S5551, dist_coeffs_41);
                bool _S5553 = !((F32_min((determinant_0(_S5552)), ((F32_min((_S5552.rows[int(0)].x), (_S5552.rows[int(1)].y)))))) > 0.0f);
                _S5481 = _S5553;
                if(_S5553)
                {
                    break;
                }
                float u_96 = uv_8[int(3)].x;
                float v_96 = uv_8[int(3)].y;
                float r2_91 = u_96 * u_96 + v_96 * v_96;
                float2  _S5554 = uv_8[int(3)] * make_float2 (1.0f + r2_91 * ((*dist_coeffs_41)[int(0)] + r2_91 * ((*dist_coeffs_41)[int(1)] + r2_91 * ((*dist_coeffs_41)[int(2)] + r2_91 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_96 * v_96 + (*dist_coeffs_41)[int(5)] * (r2_91 + 2.0f * u_96 * u_96) + (*dist_coeffs_41)[int(6)] * r2_91, 2.0f * (*dist_coeffs_41)[int(5)] * u_96 * v_96 + (*dist_coeffs_41)[int(4)] * (r2_91 + 2.0f * v_96 * v_96) + (*dist_coeffs_41)[int(7)] * r2_91);
                float2  _S5555 = _S5554 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5554.x + (*dist_coeffs_41)[int(9)] * _S5554.y, 0.0f);
                uv_8[int(3)] = make_float2 (fx_37 * _S5555.x + cx_32, fy_37 * _S5555.y + cy_32);
                break;
            }
            bool all_valid_26 = all_valid_25 & (!_S5481);
            float3  _S5556 = pos_c_1[int(4)];
            _S5482 = &uv_8[int(4)];
            for(;;)
            {
                float2  _S5557 = float2 {_S5556.x, _S5556.y};
                float r_40 = length_0(_S5557);
                float _S5558 = _S5556.z;
                float theta_36 = (F32_atan2((r_40), (_S5558)));
                if(theta_36 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_36 * theta_36 / 3.0f) / _S5558;
                }
                else
                {
                    k_13 = theta_36 / r_40;
                }
                float2  _S5559 = _S5557 * make_float2 (k_13);
                uv_8[int(4)] = _S5559;
                Matrix<float, 2, 2>  _S5560 = camera_distortion_jac_0(_S5559, dist_coeffs_41);
                bool _S5561 = !((F32_min((determinant_0(_S5560)), ((F32_min((_S5560.rows[int(0)].x), (_S5560.rows[int(1)].y)))))) > 0.0f);
                _S5483 = _S5561;
                if(_S5561)
                {
                    break;
                }
                float u_97 = uv_8[int(4)].x;
                float v_97 = uv_8[int(4)].y;
                float r2_92 = u_97 * u_97 + v_97 * v_97;
                float2  _S5562 = uv_8[int(4)] * make_float2 (1.0f + r2_92 * ((*dist_coeffs_41)[int(0)] + r2_92 * ((*dist_coeffs_41)[int(1)] + r2_92 * ((*dist_coeffs_41)[int(2)] + r2_92 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_97 * v_97 + (*dist_coeffs_41)[int(5)] * (r2_92 + 2.0f * u_97 * u_97) + (*dist_coeffs_41)[int(6)] * r2_92, 2.0f * (*dist_coeffs_41)[int(5)] * u_97 * v_97 + (*dist_coeffs_41)[int(4)] * (r2_92 + 2.0f * v_97 * v_97) + (*dist_coeffs_41)[int(7)] * r2_92);
                float2  _S5563 = _S5562 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5562.x + (*dist_coeffs_41)[int(9)] * _S5562.y, 0.0f);
                uv_8[int(4)] = make_float2 (fx_37 * _S5563.x + cx_32, fy_37 * _S5563.y + cy_32);
                break;
            }
            bool all_valid_27 = all_valid_26 & (!_S5483);
            float3  _S5564 = pos_c_1[int(5)];
            _S5484 = &uv_8[int(5)];
            for(;;)
            {
                float2  _S5565 = float2 {_S5564.x, _S5564.y};
                float r_41 = length_0(_S5565);
                float _S5566 = _S5564.z;
                float theta_37 = (F32_atan2((r_41), (_S5566)));
                if(theta_37 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_37 * theta_37 / 3.0f) / _S5566;
                }
                else
                {
                    k_13 = theta_37 / r_41;
                }
                float2  _S5567 = _S5565 * make_float2 (k_13);
                uv_8[int(5)] = _S5567;
                Matrix<float, 2, 2>  _S5568 = camera_distortion_jac_0(_S5567, dist_coeffs_41);
                bool _S5569 = !((F32_min((determinant_0(_S5568)), ((F32_min((_S5568.rows[int(0)].x), (_S5568.rows[int(1)].y)))))) > 0.0f);
                _S5485 = _S5569;
                if(_S5569)
                {
                    break;
                }
                float u_98 = uv_8[int(5)].x;
                float v_98 = uv_8[int(5)].y;
                float r2_93 = u_98 * u_98 + v_98 * v_98;
                float2  _S5570 = uv_8[int(5)] * make_float2 (1.0f + r2_93 * ((*dist_coeffs_41)[int(0)] + r2_93 * ((*dist_coeffs_41)[int(1)] + r2_93 * ((*dist_coeffs_41)[int(2)] + r2_93 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_98 * v_98 + (*dist_coeffs_41)[int(5)] * (r2_93 + 2.0f * u_98 * u_98) + (*dist_coeffs_41)[int(6)] * r2_93, 2.0f * (*dist_coeffs_41)[int(5)] * u_98 * v_98 + (*dist_coeffs_41)[int(4)] * (r2_93 + 2.0f * v_98 * v_98) + (*dist_coeffs_41)[int(7)] * r2_93);
                float2  _S5571 = _S5570 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5570.x + (*dist_coeffs_41)[int(9)] * _S5570.y, 0.0f);
                uv_8[int(5)] = make_float2 (fx_37 * _S5571.x + cx_32, fy_37 * _S5571.y + cy_32);
                break;
            }
            bool all_valid_28 = all_valid_27 & (!_S5485);
            float3  _S5572 = pos_c_1[int(6)];
            _S5486 = &uv_8[int(6)];
            for(;;)
            {
                float2  _S5573 = float2 {_S5572.x, _S5572.y};
                float r_42 = length_0(_S5573);
                float _S5574 = _S5572.z;
                float theta_38 = (F32_atan2((r_42), (_S5574)));
                if(theta_38 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_38 * theta_38 / 3.0f) / _S5574;
                }
                else
                {
                    k_13 = theta_38 / r_42;
                }
                float2  _S5575 = _S5573 * make_float2 (k_13);
                uv_8[int(6)] = _S5575;
                Matrix<float, 2, 2>  _S5576 = camera_distortion_jac_0(_S5575, dist_coeffs_41);
                bool _S5577 = !((F32_min((determinant_0(_S5576)), ((F32_min((_S5576.rows[int(0)].x), (_S5576.rows[int(1)].y)))))) > 0.0f);
                _S5487 = _S5577;
                if(_S5577)
                {
                    break;
                }
                float u_99 = uv_8[int(6)].x;
                float v_99 = uv_8[int(6)].y;
                float r2_94 = u_99 * u_99 + v_99 * v_99;
                float2  _S5578 = uv_8[int(6)] * make_float2 (1.0f + r2_94 * ((*dist_coeffs_41)[int(0)] + r2_94 * ((*dist_coeffs_41)[int(1)] + r2_94 * ((*dist_coeffs_41)[int(2)] + r2_94 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_99 * v_99 + (*dist_coeffs_41)[int(5)] * (r2_94 + 2.0f * u_99 * u_99) + (*dist_coeffs_41)[int(6)] * r2_94, 2.0f * (*dist_coeffs_41)[int(5)] * u_99 * v_99 + (*dist_coeffs_41)[int(4)] * (r2_94 + 2.0f * v_99 * v_99) + (*dist_coeffs_41)[int(7)] * r2_94);
                float2  _S5579 = _S5578 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5578.x + (*dist_coeffs_41)[int(9)] * _S5578.y, 0.0f);
                uv_8[int(6)] = make_float2 (fx_37 * _S5579.x + cx_32, fy_37 * _S5579.y + cy_32);
                break;
            }
            bool all_valid_29 = all_valid_28 & (!_S5487);
            float3  _S5580 = pos_c_1[int(7)];
            _S5488 = &uv_8[int(7)];
            for(;;)
            {
                float2  _S5581 = float2 {_S5580.x, _S5580.y};
                float r_43 = length_0(_S5581);
                float _S5582 = _S5580.z;
                float theta_39 = (F32_atan2((r_43), (_S5582)));
                if(theta_39 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_39 * theta_39 / 3.0f) / _S5582;
                }
                else
                {
                    k_13 = theta_39 / r_43;
                }
                float2  _S5583 = _S5581 * make_float2 (k_13);
                uv_8[int(7)] = _S5583;
                Matrix<float, 2, 2>  _S5584 = camera_distortion_jac_0(_S5583, dist_coeffs_41);
                bool _S5585 = !((F32_min((determinant_0(_S5584)), ((F32_min((_S5584.rows[int(0)].x), (_S5584.rows[int(1)].y)))))) > 0.0f);
                _S5489 = _S5585;
                if(_S5585)
                {
                    break;
                }
                float u_100 = uv_8[int(7)].x;
                float v_100 = uv_8[int(7)].y;
                float r2_95 = u_100 * u_100 + v_100 * v_100;
                float2  _S5586 = uv_8[int(7)] * make_float2 (1.0f + r2_95 * ((*dist_coeffs_41)[int(0)] + r2_95 * ((*dist_coeffs_41)[int(1)] + r2_95 * ((*dist_coeffs_41)[int(2)] + r2_95 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_100 * v_100 + (*dist_coeffs_41)[int(5)] * (r2_95 + 2.0f * u_100 * u_100) + (*dist_coeffs_41)[int(6)] * r2_95, 2.0f * (*dist_coeffs_41)[int(5)] * u_100 * v_100 + (*dist_coeffs_41)[int(4)] * (r2_95 + 2.0f * v_100 * v_100) + (*dist_coeffs_41)[int(7)] * r2_95);
                float2  _S5587 = _S5586 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5586.x + (*dist_coeffs_41)[int(9)] * _S5586.y, 0.0f);
                uv_8[int(7)] = make_float2 (fx_37 * _S5587.x + cx_32, fy_37 * _S5587.y + cy_32);
                break;
            }
            _S5490 = all_valid_29 & (!_S5489);
            break;
        }
        if(!_S5490)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*densities_1)[int(0)]), ((*densities_1)[int(1)])))), ((*densities_1)[int(2)])))), ((*densities_1)[int(3)])))), ((*densities_1)[int(4)])))), ((*densities_1)[int(5)])))), ((*densities_1)[int(6)])))), ((*densities_1)[int(7)]))) * size_1 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5588 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5474).x), ((*_S5476).x)))), ((*_S5478).x)))), ((*_S5480).x)))), ((*_S5482).x)))), ((*_S5484).x)))), ((*_S5486).x)))), ((*_S5488).x)));
        float _S5589 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5474).x), ((*_S5476).x)))), ((*_S5478).x)))), ((*_S5480).x)))), ((*_S5482).x)))), ((*_S5484).x)))), ((*_S5486).x)))), ((*_S5488).x)));
        float _S5590 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5474).y), ((*_S5476).y)))), ((*_S5478).y)))), ((*_S5480).y)))), ((*_S5482).y)))), ((*_S5484).y)))), ((*_S5486).y)))), ((*_S5488).y)));
        float _S5591 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5474).y), ((*_S5476).y)))), ((*_S5478).y)))), ((*_S5480).y)))), ((*_S5482).y)))), ((*_S5484).y)))), ((*_S5486).y)))), ((*_S5488).y)));
        if(_S5588 <= 0.0f)
        {
            _S5523 = true;
        }
        else
        {
            _S5523 = _S5589 >= float(image_width_28);
        }
        if(_S5523)
        {
            _S5523 = true;
        }
        else
        {
            _S5523 = _S5590 <= 0.0f;
        }
        if(_S5523)
        {
            _S5523 = true;
        }
        else
        {
            _S5523 = _S5591 >= float(image_height_28);
        }
        if(_S5523)
        {
            _S5523 = true;
        }
        else
        {
            if(_S5521 <= 0.0f)
            {
                if(_S5589 <= 0.0f)
                {
                    _S5523 = _S5588 >= float(image_width_28);
                }
                else
                {
                    _S5523 = false;
                }
                if(_S5523)
                {
                    _S5523 = true;
                }
                else
                {
                    if(_S5591 <= 0.0f)
                    {
                        _S5523 = _S5590 >= float(image_width_28);
                    }
                    else
                    {
                        _S5523 = false;
                    }
                }
            }
            else
            {
                _S5523 = false;
            }
        }
        if(_S5523)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_19 = make_int4 (int((F32_floor((_S5589)))), int((F32_floor((_S5591)))), int((F32_ceil((_S5588)))), int((F32_ceil((_S5590)))));
        *depth_22 = 0.5f * (F32_log((dot_0(mean_c_28, mean_c_28) + 9.99999997475242708e-07f)));
        float3  _S5592 = mean_c_28 - - mul_0(transpose_0(R_34), t_34);
        float3  _S5593 = make_float3 (0.282094806432724f) * (*sh_coeffs_28)[int(0)];
        *rgbs_8 = _S5593;
        float _S5594 = _S5592.x;
        float _S5595 = _S5592.y;
        float _S5596 = _S5592.z;
        float norm_19 = (F32_sqrt((_S5594 * _S5594 + _S5595 * _S5595 + _S5596 * _S5596)));
        float x_64 = _S5594 / norm_19;
        float y_31 = _S5595 / norm_19;
        float z_28 = _S5596 / norm_19;
        float3  _S5597 = _S5593 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_28)[int(1)] + make_float3 (z_28) * (*sh_coeffs_28)[int(2)] - make_float3 (x_64) * (*sh_coeffs_28)[int(3)]);
        *rgbs_8 = _S5597;
        float z2_61 = z_28 * z_28;
        float fTmp0B_28 = -1.09254848957061768f * z_28;
        float fC1_28 = x_64 * x_64 - y_31 * y_31;
        float fS1_28 = 2.0f * x_64 * y_31;
        float3  _S5598 = _S5597 + (make_float3 (0.54627424478530884f * fS1_28) * (*sh_coeffs_28)[int(4)] + make_float3 (fTmp0B_28 * y_31) * (*sh_coeffs_28)[int(5)] + make_float3 (0.94617468118667603f * z2_61 - 0.31539157032966614f) * (*sh_coeffs_28)[int(6)] + make_float3 (fTmp0B_28 * x_64) * (*sh_coeffs_28)[int(7)] + make_float3 (0.54627424478530884f * fC1_28) * (*sh_coeffs_28)[int(8)]);
        *rgbs_8 = _S5598;
        float fTmp0C_28 = -2.28522896766662598f * z2_61 + 0.4570457935333252f;
        float fTmp1B_28 = 1.44530570507049561f * z_28;
        *rgbs_8 = max_0(_S5598 + (make_float3 (-0.59004360437393188f * (x_64 * fS1_28 + y_31 * fC1_28)) * (*sh_coeffs_28)[int(9)] + make_float3 (fTmp1B_28 * fS1_28) * (*sh_coeffs_28)[int(10)] + make_float3 (fTmp0C_28 * y_31) * (*sh_coeffs_28)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_61 - 1.11952900886535645f)) * (*sh_coeffs_28)[int(12)] + make_float3 (fTmp0C_28 * x_64) * (*sh_coeffs_28)[int(13)] + make_float3 (fTmp1B_28 * fC1_28) * (*sh_coeffs_28)[int(14)] + make_float3 (-0.59004360437393188f * (x_64 * fC1_28 - y_31 * fS1_28)) * (*sh_coeffs_28)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  * densities_2, FixedArray<float3 , 16>  * sh_coeffs_29, Matrix<float, 3, 3>  R_35, float3  t_35, float fx_38, float fy_38, float cx_33, float cy_33, FixedArray<float, 10>  * dist_coeffs_42, uint image_width_29, uint image_height_29, float near_plane_20, float far_plane_20, int4  * aabb_xyxy_20, float * depth_23, float3  * rgbs_9)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S5599 = mul_0(R_35, pos_2) + t_35;
    pos_c_2[int(0)] = _S5599;
    float3  _S5600 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_35;
    pos_c_2[int(1)] = _S5600;
    float3  _S5601 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(2)] = _S5601;
    float3  _S5602 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(3)] = _S5602;
    float3  _S5603 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(4)] = _S5603;
    float3  _S5604 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(5)] = _S5604;
    float3  _S5605 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_35;
    pos_c_2[int(6)] = _S5605;
    float3  _S5606 = mul_0(R_35, pos_2 + make_float3 (size_2)) + t_35;
    pos_c_2[int(7)] = _S5606;
    float3  mean_c_29 = mul_0(R_35, pos_2 + make_float3 (0.5f * size_2)) + t_35;
    FixedArray<float2 , 8>  uv_9;
    float2  _S5607 = float2 {_S5599.x, _S5599.y} / make_float2 (_S5599.z);
    float u_101 = _S5607.x;
    float v_101 = _S5607.y;
    float r2_96 = u_101 * u_101 + v_101 * v_101;
    float _S5608 = 2.0f * (*dist_coeffs_42)[int(4)];
    float _S5609 = 2.0f * (*dist_coeffs_42)[int(5)];
    float2  _S5610 = _S5607 * make_float2 (1.0f + r2_96 * ((*dist_coeffs_42)[int(0)] + r2_96 * ((*dist_coeffs_42)[int(1)] + r2_96 * ((*dist_coeffs_42)[int(2)] + r2_96 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_101 * v_101 + (*dist_coeffs_42)[int(5)] * (r2_96 + 2.0f * u_101 * u_101) + (*dist_coeffs_42)[int(6)] * r2_96, _S5609 * u_101 * v_101 + (*dist_coeffs_42)[int(4)] * (r2_96 + 2.0f * v_101 * v_101) + (*dist_coeffs_42)[int(7)] * r2_96);
    float2  _S5611 = _S5610 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5610.x + (*dist_coeffs_42)[int(9)] * _S5610.y, 0.0f);
    float _S5612 = fx_38 * _S5611.x + cx_33;
    float _S5613 = fy_38 * _S5611.y + cy_33;
    uv_9[int(0)] = make_float2 (_S5612, _S5613);
    float2  _S5614 = float2 {_S5600.x, _S5600.y} / make_float2 (_S5600.z);
    float u_102 = _S5614.x;
    float v_102 = _S5614.y;
    float r2_97 = u_102 * u_102 + v_102 * v_102;
    float2  _S5615 = _S5614 * make_float2 (1.0f + r2_97 * ((*dist_coeffs_42)[int(0)] + r2_97 * ((*dist_coeffs_42)[int(1)] + r2_97 * ((*dist_coeffs_42)[int(2)] + r2_97 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_102 * v_102 + (*dist_coeffs_42)[int(5)] * (r2_97 + 2.0f * u_102 * u_102) + (*dist_coeffs_42)[int(6)] * r2_97, _S5609 * u_102 * v_102 + (*dist_coeffs_42)[int(4)] * (r2_97 + 2.0f * v_102 * v_102) + (*dist_coeffs_42)[int(7)] * r2_97);
    float2  _S5616 = _S5615 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5615.x + (*dist_coeffs_42)[int(9)] * _S5615.y, 0.0f);
    float _S5617 = fx_38 * _S5616.x + cx_33;
    float _S5618 = fy_38 * _S5616.y + cy_33;
    uv_9[int(1)] = make_float2 (_S5617, _S5618);
    float2  _S5619 = float2 {_S5601.x, _S5601.y} / make_float2 (_S5601.z);
    float u_103 = _S5619.x;
    float v_103 = _S5619.y;
    float r2_98 = u_103 * u_103 + v_103 * v_103;
    float2  _S5620 = _S5619 * make_float2 (1.0f + r2_98 * ((*dist_coeffs_42)[int(0)] + r2_98 * ((*dist_coeffs_42)[int(1)] + r2_98 * ((*dist_coeffs_42)[int(2)] + r2_98 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_103 * v_103 + (*dist_coeffs_42)[int(5)] * (r2_98 + 2.0f * u_103 * u_103) + (*dist_coeffs_42)[int(6)] * r2_98, _S5609 * u_103 * v_103 + (*dist_coeffs_42)[int(4)] * (r2_98 + 2.0f * v_103 * v_103) + (*dist_coeffs_42)[int(7)] * r2_98);
    float2  _S5621 = _S5620 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5620.x + (*dist_coeffs_42)[int(9)] * _S5620.y, 0.0f);
    float _S5622 = fx_38 * _S5621.x + cx_33;
    float _S5623 = fy_38 * _S5621.y + cy_33;
    uv_9[int(2)] = make_float2 (_S5622, _S5623);
    float2  _S5624 = float2 {_S5602.x, _S5602.y} / make_float2 (_S5602.z);
    float u_104 = _S5624.x;
    float v_104 = _S5624.y;
    float r2_99 = u_104 * u_104 + v_104 * v_104;
    float2  _S5625 = _S5624 * make_float2 (1.0f + r2_99 * ((*dist_coeffs_42)[int(0)] + r2_99 * ((*dist_coeffs_42)[int(1)] + r2_99 * ((*dist_coeffs_42)[int(2)] + r2_99 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_104 * v_104 + (*dist_coeffs_42)[int(5)] * (r2_99 + 2.0f * u_104 * u_104) + (*dist_coeffs_42)[int(6)] * r2_99, _S5609 * u_104 * v_104 + (*dist_coeffs_42)[int(4)] * (r2_99 + 2.0f * v_104 * v_104) + (*dist_coeffs_42)[int(7)] * r2_99);
    float2  _S5626 = _S5625 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5625.x + (*dist_coeffs_42)[int(9)] * _S5625.y, 0.0f);
    float _S5627 = fx_38 * _S5626.x + cx_33;
    float _S5628 = fy_38 * _S5626.y + cy_33;
    uv_9[int(3)] = make_float2 (_S5627, _S5628);
    float2  _S5629 = float2 {_S5603.x, _S5603.y} / make_float2 (_S5603.z);
    float u_105 = _S5629.x;
    float v_105 = _S5629.y;
    float r2_100 = u_105 * u_105 + v_105 * v_105;
    float2  _S5630 = _S5629 * make_float2 (1.0f + r2_100 * ((*dist_coeffs_42)[int(0)] + r2_100 * ((*dist_coeffs_42)[int(1)] + r2_100 * ((*dist_coeffs_42)[int(2)] + r2_100 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_105 * v_105 + (*dist_coeffs_42)[int(5)] * (r2_100 + 2.0f * u_105 * u_105) + (*dist_coeffs_42)[int(6)] * r2_100, _S5609 * u_105 * v_105 + (*dist_coeffs_42)[int(4)] * (r2_100 + 2.0f * v_105 * v_105) + (*dist_coeffs_42)[int(7)] * r2_100);
    float2  _S5631 = _S5630 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5630.x + (*dist_coeffs_42)[int(9)] * _S5630.y, 0.0f);
    float _S5632 = fx_38 * _S5631.x + cx_33;
    float _S5633 = fy_38 * _S5631.y + cy_33;
    uv_9[int(4)] = make_float2 (_S5632, _S5633);
    float2  _S5634 = float2 {_S5604.x, _S5604.y} / make_float2 (_S5604.z);
    float u_106 = _S5634.x;
    float v_106 = _S5634.y;
    float r2_101 = u_106 * u_106 + v_106 * v_106;
    float2  _S5635 = _S5634 * make_float2 (1.0f + r2_101 * ((*dist_coeffs_42)[int(0)] + r2_101 * ((*dist_coeffs_42)[int(1)] + r2_101 * ((*dist_coeffs_42)[int(2)] + r2_101 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_106 * v_106 + (*dist_coeffs_42)[int(5)] * (r2_101 + 2.0f * u_106 * u_106) + (*dist_coeffs_42)[int(6)] * r2_101, _S5609 * u_106 * v_106 + (*dist_coeffs_42)[int(4)] * (r2_101 + 2.0f * v_106 * v_106) + (*dist_coeffs_42)[int(7)] * r2_101);
    float2  _S5636 = _S5635 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5635.x + (*dist_coeffs_42)[int(9)] * _S5635.y, 0.0f);
    float _S5637 = fx_38 * _S5636.x + cx_33;
    float _S5638 = fy_38 * _S5636.y + cy_33;
    uv_9[int(5)] = make_float2 (_S5637, _S5638);
    float2  _S5639 = float2 {_S5605.x, _S5605.y} / make_float2 (_S5605.z);
    float u_107 = _S5639.x;
    float v_107 = _S5639.y;
    float r2_102 = u_107 * u_107 + v_107 * v_107;
    float2  _S5640 = _S5639 * make_float2 (1.0f + r2_102 * ((*dist_coeffs_42)[int(0)] + r2_102 * ((*dist_coeffs_42)[int(1)] + r2_102 * ((*dist_coeffs_42)[int(2)] + r2_102 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_107 * v_107 + (*dist_coeffs_42)[int(5)] * (r2_102 + 2.0f * u_107 * u_107) + (*dist_coeffs_42)[int(6)] * r2_102, _S5609 * u_107 * v_107 + (*dist_coeffs_42)[int(4)] * (r2_102 + 2.0f * v_107 * v_107) + (*dist_coeffs_42)[int(7)] * r2_102);
    float2  _S5641 = _S5640 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5640.x + (*dist_coeffs_42)[int(9)] * _S5640.y, 0.0f);
    float _S5642 = fx_38 * _S5641.x + cx_33;
    float _S5643 = fy_38 * _S5641.y + cy_33;
    uv_9[int(6)] = make_float2 (_S5642, _S5643);
    float2  _S5644 = float2 {_S5606.x, _S5606.y} / make_float2 (_S5606.z);
    float u_108 = _S5644.x;
    float v_108 = _S5644.y;
    float r2_103 = u_108 * u_108 + v_108 * v_108;
    float2  _S5645 = _S5644 * make_float2 (1.0f + r2_103 * ((*dist_coeffs_42)[int(0)] + r2_103 * ((*dist_coeffs_42)[int(1)] + r2_103 * ((*dist_coeffs_42)[int(2)] + r2_103 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5608 * u_108 * v_108 + (*dist_coeffs_42)[int(5)] * (r2_103 + 2.0f * u_108 * u_108) + (*dist_coeffs_42)[int(6)] * r2_103, _S5609 * u_108 * v_108 + (*dist_coeffs_42)[int(4)] * (r2_103 + 2.0f * v_108 * v_108) + (*dist_coeffs_42)[int(7)] * r2_103);
    float2  _S5646 = _S5645 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5645.x + (*dist_coeffs_42)[int(9)] * _S5645.y, 0.0f);
    float _S5647 = fx_38 * _S5646.x + cx_33;
    float _S5648 = fy_38 * _S5646.y + cy_33;
    uv_9[int(7)] = make_float2 (_S5647, _S5648);
    *aabb_xyxy_20 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5612), (_S5617)))), (_S5622)))), (_S5627)))), (_S5632)))), (_S5637)))), (_S5642)))), (_S5647))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5613), (_S5618)))), (_S5623)))), (_S5628)))), (_S5633)))), (_S5638)))), (_S5643)))), (_S5648))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5612), (_S5617)))), (_S5622)))), (_S5627)))), (_S5632)))), (_S5637)))), (_S5642)))), (_S5647))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5613), (_S5618)))), (_S5623)))), (_S5628)))), (_S5633)))), (_S5638)))), (_S5643)))), (_S5648))))))));
    *depth_23 = 0.5f * (F32_log((dot_0(mean_c_29, mean_c_29) + 9.99999997475242708e-07f)));
    float3  _S5649 = mean_c_29 - - mul_0(transpose_0(R_35), t_35);
    float3  _S5650 = make_float3 (0.282094806432724f) * (*sh_coeffs_29)[int(0)];
    *rgbs_9 = _S5650;
    float _S5651 = _S5649.x;
    float _S5652 = _S5649.y;
    float _S5653 = _S5649.z;
    float norm_20 = (F32_sqrt((_S5651 * _S5651 + _S5652 * _S5652 + _S5653 * _S5653)));
    float x_65 = _S5651 / norm_20;
    float y_32 = _S5652 / norm_20;
    float z_29 = _S5653 / norm_20;
    float3  _S5654 = _S5650 + make_float3 (0.48860251903533936f) * (make_float3 (- y_32) * (*sh_coeffs_29)[int(1)] + make_float3 (z_29) * (*sh_coeffs_29)[int(2)] - make_float3 (x_65) * (*sh_coeffs_29)[int(3)]);
    *rgbs_9 = _S5654;
    float z2_62 = z_29 * z_29;
    float fTmp0B_29 = -1.09254848957061768f * z_29;
    float fC1_29 = x_65 * x_65 - y_32 * y_32;
    float fS1_29 = 2.0f * x_65 * y_32;
    float3  _S5655 = _S5654 + (make_float3 (0.54627424478530884f * fS1_29) * (*sh_coeffs_29)[int(4)] + make_float3 (fTmp0B_29 * y_32) * (*sh_coeffs_29)[int(5)] + make_float3 (0.94617468118667603f * z2_62 - 0.31539157032966614f) * (*sh_coeffs_29)[int(6)] + make_float3 (fTmp0B_29 * x_65) * (*sh_coeffs_29)[int(7)] + make_float3 (0.54627424478530884f * fC1_29) * (*sh_coeffs_29)[int(8)]);
    *rgbs_9 = _S5655;
    float fTmp0C_29 = -2.28522896766662598f * z2_62 + 0.4570457935333252f;
    float fTmp1B_29 = 1.44530570507049561f * z_29;
    *rgbs_9 = max_0(_S5655 + (make_float3 (-0.59004360437393188f * (x_65 * fS1_29 + y_32 * fC1_29)) * (*sh_coeffs_29)[int(9)] + make_float3 (fTmp1B_29 * fS1_29) * (*sh_coeffs_29)[int(10)] + make_float3 (fTmp0C_29 * y_32) * (*sh_coeffs_29)[int(11)] + make_float3 (z_29 * (1.86588168144226074f * z2_62 - 1.11952900886535645f)) * (*sh_coeffs_29)[int(12)] + make_float3 (fTmp0C_29 * x_65) * (*sh_coeffs_29)[int(13)] + make_float3 (fTmp1B_29 * fC1_29) * (*sh_coeffs_29)[int(14)] + make_float3 (-0.59004360437393188f * (x_65 * fC1_29 - y_32 * fS1_29)) * (*sh_coeffs_29)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  * densities_3, FixedArray<float3 , 16>  * sh_coeffs_30, Matrix<float, 3, 3>  R_36, float3  t_36, float fx_39, float fy_39, float cx_34, float cy_34, FixedArray<float, 10>  * dist_coeffs_43, uint image_width_30, uint image_height_30, float near_plane_21, float far_plane_21, int4  * aabb_xyxy_21, float * depth_24, float3  * rgbs_10)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S5656 = mul_0(R_36, pos_3) + t_36;
    pos_c_3[int(0)] = _S5656;
    pos_c_3[int(1)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_36;
    pos_c_3[int(2)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(3)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(4)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(5)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(6)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_36;
    pos_c_3[int(7)] = mul_0(R_36, pos_3 + make_float3 (size_3)) + t_36;
    float3  mean_c_30 = mul_0(R_36, pos_3 + make_float3 (0.5f * size_3)) + t_36;
    FixedArray<float2 , 8>  uv_10;
    float2  _S5657 = float2 {_S5656.x, _S5656.y};
    float r_44 = length_0(_S5657);
    float _S5658 = _S5656.z;
    float theta_40 = (F32_atan2((r_44), (_S5658)));
    float k_14;
    if(theta_40 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_40 * theta_40 / 3.0f) / _S5658;
    }
    else
    {
        k_14 = theta_40 / r_44;
    }
    float2  _S5659 = _S5657 * make_float2 (k_14);
    float u_109 = _S5659.x;
    float v_109 = _S5659.y;
    float r2_104 = u_109 * u_109 + v_109 * v_109;
    float _S5660 = 2.0f * (*dist_coeffs_43)[int(4)];
    float _S5661 = 2.0f * (*dist_coeffs_43)[int(5)];
    float2  _S5662 = _S5659 * make_float2 (1.0f + r2_104 * ((*dist_coeffs_43)[int(0)] + r2_104 * ((*dist_coeffs_43)[int(1)] + r2_104 * ((*dist_coeffs_43)[int(2)] + r2_104 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_109 * v_109 + (*dist_coeffs_43)[int(5)] * (r2_104 + 2.0f * u_109 * u_109) + (*dist_coeffs_43)[int(6)] * r2_104, _S5661 * u_109 * v_109 + (*dist_coeffs_43)[int(4)] * (r2_104 + 2.0f * v_109 * v_109) + (*dist_coeffs_43)[int(7)] * r2_104);
    float2  _S5663 = _S5662 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5662.x + (*dist_coeffs_43)[int(9)] * _S5662.y, 0.0f);
    uv_10[int(0)] = make_float2 (fx_39 * _S5663.x + cx_34, fy_39 * _S5663.y + cy_34);
    float2  _S5664 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_45 = length_0(_S5664);
    float _S5665 = pos_c_3[int(1)].z;
    float theta_41 = (F32_atan2((r_45), (_S5665)));
    if(theta_41 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_41 * theta_41 / 3.0f) / _S5665;
    }
    else
    {
        k_14 = theta_41 / r_45;
    }
    float2  _S5666 = _S5664 * make_float2 (k_14);
    float u_110 = _S5666.x;
    float v_110 = _S5666.y;
    float r2_105 = u_110 * u_110 + v_110 * v_110;
    float2  _S5667 = _S5666 * make_float2 (1.0f + r2_105 * ((*dist_coeffs_43)[int(0)] + r2_105 * ((*dist_coeffs_43)[int(1)] + r2_105 * ((*dist_coeffs_43)[int(2)] + r2_105 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_110 * v_110 + (*dist_coeffs_43)[int(5)] * (r2_105 + 2.0f * u_110 * u_110) + (*dist_coeffs_43)[int(6)] * r2_105, _S5661 * u_110 * v_110 + (*dist_coeffs_43)[int(4)] * (r2_105 + 2.0f * v_110 * v_110) + (*dist_coeffs_43)[int(7)] * r2_105);
    float2  _S5668 = _S5667 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5667.x + (*dist_coeffs_43)[int(9)] * _S5667.y, 0.0f);
    uv_10[int(1)] = make_float2 (fx_39 * _S5668.x + cx_34, fy_39 * _S5668.y + cy_34);
    float2  _S5669 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_46 = length_0(_S5669);
    float _S5670 = pos_c_3[int(2)].z;
    float theta_42 = (F32_atan2((r_46), (_S5670)));
    if(theta_42 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_42 * theta_42 / 3.0f) / _S5670;
    }
    else
    {
        k_14 = theta_42 / r_46;
    }
    float2  _S5671 = _S5669 * make_float2 (k_14);
    float u_111 = _S5671.x;
    float v_111 = _S5671.y;
    float r2_106 = u_111 * u_111 + v_111 * v_111;
    float2  _S5672 = _S5671 * make_float2 (1.0f + r2_106 * ((*dist_coeffs_43)[int(0)] + r2_106 * ((*dist_coeffs_43)[int(1)] + r2_106 * ((*dist_coeffs_43)[int(2)] + r2_106 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_111 * v_111 + (*dist_coeffs_43)[int(5)] * (r2_106 + 2.0f * u_111 * u_111) + (*dist_coeffs_43)[int(6)] * r2_106, _S5661 * u_111 * v_111 + (*dist_coeffs_43)[int(4)] * (r2_106 + 2.0f * v_111 * v_111) + (*dist_coeffs_43)[int(7)] * r2_106);
    float2  _S5673 = _S5672 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5672.x + (*dist_coeffs_43)[int(9)] * _S5672.y, 0.0f);
    uv_10[int(2)] = make_float2 (fx_39 * _S5673.x + cx_34, fy_39 * _S5673.y + cy_34);
    float2  _S5674 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_47 = length_0(_S5674);
    float _S5675 = pos_c_3[int(3)].z;
    float theta_43 = (F32_atan2((r_47), (_S5675)));
    if(theta_43 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_43 * theta_43 / 3.0f) / _S5675;
    }
    else
    {
        k_14 = theta_43 / r_47;
    }
    float2  _S5676 = _S5674 * make_float2 (k_14);
    float u_112 = _S5676.x;
    float v_112 = _S5676.y;
    float r2_107 = u_112 * u_112 + v_112 * v_112;
    float2  _S5677 = _S5676 * make_float2 (1.0f + r2_107 * ((*dist_coeffs_43)[int(0)] + r2_107 * ((*dist_coeffs_43)[int(1)] + r2_107 * ((*dist_coeffs_43)[int(2)] + r2_107 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_112 * v_112 + (*dist_coeffs_43)[int(5)] * (r2_107 + 2.0f * u_112 * u_112) + (*dist_coeffs_43)[int(6)] * r2_107, _S5661 * u_112 * v_112 + (*dist_coeffs_43)[int(4)] * (r2_107 + 2.0f * v_112 * v_112) + (*dist_coeffs_43)[int(7)] * r2_107);
    float2  _S5678 = _S5677 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5677.x + (*dist_coeffs_43)[int(9)] * _S5677.y, 0.0f);
    uv_10[int(3)] = make_float2 (fx_39 * _S5678.x + cx_34, fy_39 * _S5678.y + cy_34);
    float2  _S5679 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_48 = length_0(_S5679);
    float _S5680 = pos_c_3[int(4)].z;
    float theta_44 = (F32_atan2((r_48), (_S5680)));
    if(theta_44 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_44 * theta_44 / 3.0f) / _S5680;
    }
    else
    {
        k_14 = theta_44 / r_48;
    }
    float2  _S5681 = _S5679 * make_float2 (k_14);
    float u_113 = _S5681.x;
    float v_113 = _S5681.y;
    float r2_108 = u_113 * u_113 + v_113 * v_113;
    float2  _S5682 = _S5681 * make_float2 (1.0f + r2_108 * ((*dist_coeffs_43)[int(0)] + r2_108 * ((*dist_coeffs_43)[int(1)] + r2_108 * ((*dist_coeffs_43)[int(2)] + r2_108 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_113 * v_113 + (*dist_coeffs_43)[int(5)] * (r2_108 + 2.0f * u_113 * u_113) + (*dist_coeffs_43)[int(6)] * r2_108, _S5661 * u_113 * v_113 + (*dist_coeffs_43)[int(4)] * (r2_108 + 2.0f * v_113 * v_113) + (*dist_coeffs_43)[int(7)] * r2_108);
    float2  _S5683 = _S5682 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5682.x + (*dist_coeffs_43)[int(9)] * _S5682.y, 0.0f);
    uv_10[int(4)] = make_float2 (fx_39 * _S5683.x + cx_34, fy_39 * _S5683.y + cy_34);
    float2  _S5684 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_49 = length_0(_S5684);
    float _S5685 = pos_c_3[int(5)].z;
    float theta_45 = (F32_atan2((r_49), (_S5685)));
    if(theta_45 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_45 * theta_45 / 3.0f) / _S5685;
    }
    else
    {
        k_14 = theta_45 / r_49;
    }
    float2  _S5686 = _S5684 * make_float2 (k_14);
    float u_114 = _S5686.x;
    float v_114 = _S5686.y;
    float r2_109 = u_114 * u_114 + v_114 * v_114;
    float2  _S5687 = _S5686 * make_float2 (1.0f + r2_109 * ((*dist_coeffs_43)[int(0)] + r2_109 * ((*dist_coeffs_43)[int(1)] + r2_109 * ((*dist_coeffs_43)[int(2)] + r2_109 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_114 * v_114 + (*dist_coeffs_43)[int(5)] * (r2_109 + 2.0f * u_114 * u_114) + (*dist_coeffs_43)[int(6)] * r2_109, _S5661 * u_114 * v_114 + (*dist_coeffs_43)[int(4)] * (r2_109 + 2.0f * v_114 * v_114) + (*dist_coeffs_43)[int(7)] * r2_109);
    float2  _S5688 = _S5687 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5687.x + (*dist_coeffs_43)[int(9)] * _S5687.y, 0.0f);
    uv_10[int(5)] = make_float2 (fx_39 * _S5688.x + cx_34, fy_39 * _S5688.y + cy_34);
    float2  _S5689 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_50 = length_0(_S5689);
    float _S5690 = pos_c_3[int(6)].z;
    float theta_46 = (F32_atan2((r_50), (_S5690)));
    if(theta_46 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_46 * theta_46 / 3.0f) / _S5690;
    }
    else
    {
        k_14 = theta_46 / r_50;
    }
    float2  _S5691 = _S5689 * make_float2 (k_14);
    float u_115 = _S5691.x;
    float v_115 = _S5691.y;
    float r2_110 = u_115 * u_115 + v_115 * v_115;
    float2  _S5692 = _S5691 * make_float2 (1.0f + r2_110 * ((*dist_coeffs_43)[int(0)] + r2_110 * ((*dist_coeffs_43)[int(1)] + r2_110 * ((*dist_coeffs_43)[int(2)] + r2_110 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_115 * v_115 + (*dist_coeffs_43)[int(5)] * (r2_110 + 2.0f * u_115 * u_115) + (*dist_coeffs_43)[int(6)] * r2_110, _S5661 * u_115 * v_115 + (*dist_coeffs_43)[int(4)] * (r2_110 + 2.0f * v_115 * v_115) + (*dist_coeffs_43)[int(7)] * r2_110);
    float2  _S5693 = _S5692 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5692.x + (*dist_coeffs_43)[int(9)] * _S5692.y, 0.0f);
    uv_10[int(6)] = make_float2 (fx_39 * _S5693.x + cx_34, fy_39 * _S5693.y + cy_34);
    float2  _S5694 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_51 = length_0(_S5694);
    float _S5695 = pos_c_3[int(7)].z;
    float theta_47 = (F32_atan2((r_51), (_S5695)));
    if(theta_47 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_47 * theta_47 / 3.0f) / _S5695;
    }
    else
    {
        k_14 = theta_47 / r_51;
    }
    float2  _S5696 = _S5694 * make_float2 (k_14);
    float u_116 = _S5696.x;
    float v_116 = _S5696.y;
    float r2_111 = u_116 * u_116 + v_116 * v_116;
    float2  _S5697 = _S5696 * make_float2 (1.0f + r2_111 * ((*dist_coeffs_43)[int(0)] + r2_111 * ((*dist_coeffs_43)[int(1)] + r2_111 * ((*dist_coeffs_43)[int(2)] + r2_111 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5660 * u_116 * v_116 + (*dist_coeffs_43)[int(5)] * (r2_111 + 2.0f * u_116 * u_116) + (*dist_coeffs_43)[int(6)] * r2_111, _S5661 * u_116 * v_116 + (*dist_coeffs_43)[int(4)] * (r2_111 + 2.0f * v_116 * v_116) + (*dist_coeffs_43)[int(7)] * r2_111);
    float2  _S5698 = _S5697 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5697.x + (*dist_coeffs_43)[int(9)] * _S5697.y, 0.0f);
    float _S5699 = fx_39 * _S5698.x + cx_34;
    float _S5700 = fy_39 * _S5698.y + cy_34;
    uv_10[int(7)] = make_float2 (_S5699, _S5700);
    *aabb_xyxy_21 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_10[int(0)].x), (uv_10[int(1)].x)))), (uv_10[int(2)].x)))), (uv_10[int(3)].x)))), (uv_10[int(4)].x)))), (uv_10[int(5)].x)))), (uv_10[int(6)].x)))), (_S5699))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_10[int(0)].y), (uv_10[int(1)].y)))), (uv_10[int(2)].y)))), (uv_10[int(3)].y)))), (uv_10[int(4)].y)))), (uv_10[int(5)].y)))), (uv_10[int(6)].y)))), (_S5700))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_10[int(0)].x), (uv_10[int(1)].x)))), (uv_10[int(2)].x)))), (uv_10[int(3)].x)))), (uv_10[int(4)].x)))), (uv_10[int(5)].x)))), (uv_10[int(6)].x)))), (_S5699))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_10[int(0)].y), (uv_10[int(1)].y)))), (uv_10[int(2)].y)))), (uv_10[int(3)].y)))), (uv_10[int(4)].y)))), (uv_10[int(5)].y)))), (uv_10[int(6)].y)))), (_S5700))))))));
    *depth_24 = 0.5f * (F32_log((dot_0(mean_c_30, mean_c_30) + 9.99999997475242708e-07f)));
    float3  _S5701 = mean_c_30 - - mul_0(transpose_0(R_36), t_36);
    float3  _S5702 = make_float3 (0.282094806432724f) * (*sh_coeffs_30)[int(0)];
    *rgbs_10 = _S5702;
    float _S5703 = _S5701.x;
    float _S5704 = _S5701.y;
    float _S5705 = _S5701.z;
    float norm_21 = (F32_sqrt((_S5703 * _S5703 + _S5704 * _S5704 + _S5705 * _S5705)));
    float x_66 = _S5703 / norm_21;
    float y_33 = _S5704 / norm_21;
    float z_30 = _S5705 / norm_21;
    float3  _S5706 = _S5702 + make_float3 (0.48860251903533936f) * (make_float3 (- y_33) * (*sh_coeffs_30)[int(1)] + make_float3 (z_30) * (*sh_coeffs_30)[int(2)] - make_float3 (x_66) * (*sh_coeffs_30)[int(3)]);
    *rgbs_10 = _S5706;
    float z2_63 = z_30 * z_30;
    float fTmp0B_30 = -1.09254848957061768f * z_30;
    float fC1_30 = x_66 * x_66 - y_33 * y_33;
    float fS1_30 = 2.0f * x_66 * y_33;
    float3  _S5707 = _S5706 + (make_float3 (0.54627424478530884f * fS1_30) * (*sh_coeffs_30)[int(4)] + make_float3 (fTmp0B_30 * y_33) * (*sh_coeffs_30)[int(5)] + make_float3 (0.94617468118667603f * z2_63 - 0.31539157032966614f) * (*sh_coeffs_30)[int(6)] + make_float3 (fTmp0B_30 * x_66) * (*sh_coeffs_30)[int(7)] + make_float3 (0.54627424478530884f * fC1_30) * (*sh_coeffs_30)[int(8)]);
    *rgbs_10 = _S5707;
    float fTmp0C_30 = -2.28522896766662598f * z2_63 + 0.4570457935333252f;
    float fTmp1B_30 = 1.44530570507049561f * z_30;
    *rgbs_10 = max_0(_S5707 + (make_float3 (-0.59004360437393188f * (x_66 * fS1_30 + y_33 * fC1_30)) * (*sh_coeffs_30)[int(9)] + make_float3 (fTmp1B_30 * fS1_30) * (*sh_coeffs_30)[int(10)] + make_float3 (fTmp0C_30 * y_33) * (*sh_coeffs_30)[int(11)] + make_float3 (z_30 * (1.86588168144226074f * z2_63 - 1.11952900886535645f)) * (*sh_coeffs_30)[int(12)] + make_float3 (fTmp0C_30 * x_66) * (*sh_coeffs_30)[int(13)] + make_float3 (fTmp1B_30 * fC1_30) * (*sh_coeffs_30)[int(14)] + make_float3 (-0.59004360437393188f * (x_66 * fC1_30 - y_33 * fS1_30)) * (*sh_coeffs_30)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  * densities_4, FixedArray<float3 , 16>  * sh_coeffs_31, Matrix<float, 3, 3>  R_37, float3  t_37, float fx_40, float fy_40, float cx_35, float cy_35, FixedArray<float, 10>  * dist_coeffs_44, uint image_width_31, uint image_height_31, float3  v_rgb_8, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_9, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_10)
{
    float3  _S5708 = s_primal_ctx_mul_0(R_37, pos_4) + t_37;
    float _S5709 = _S5708.z;
    float2  _S5710 = make_float2 (_S5709);
    float _S5711 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S5709);
    float _S5712 = s_primal_ctx_max_0(0.0f, _S5709);
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S5713 = s_primal_ctx_mul_0(R_37, pos_i_0) + t_37;
    float _S5714 = _S5713.z;
    float2  _S5715 = make_float2 (_S5714);
    float _S5716 = s_primal_ctx_min_0(_S5711, _S5714);
    float _S5717 = s_primal_ctx_max_0(_S5712, _S5714);
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S5718 = s_primal_ctx_mul_0(R_37, pos_i_1) + t_37;
    float _S5719 = _S5718.z;
    float2  _S5720 = make_float2 (_S5719);
    float _S5721 = s_primal_ctx_min_0(_S5716, _S5719);
    float _S5722 = s_primal_ctx_max_0(_S5717, _S5719);
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S5723 = s_primal_ctx_mul_0(R_37, pos_i_2) + t_37;
    float _S5724 = _S5723.z;
    float2  _S5725 = make_float2 (_S5724);
    float _S5726 = s_primal_ctx_min_0(_S5721, _S5724);
    float _S5727 = s_primal_ctx_max_0(_S5722, _S5724);
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S5728 = s_primal_ctx_mul_0(R_37, pos_i_3) + t_37;
    float _S5729 = _S5728.z;
    float2  _S5730 = make_float2 (_S5729);
    float _S5731 = s_primal_ctx_min_0(_S5726, _S5729);
    float _S5732 = s_primal_ctx_max_0(_S5727, _S5729);
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S5733 = s_primal_ctx_mul_0(R_37, pos_i_4) + t_37;
    float _S5734 = _S5733.z;
    float2  _S5735 = make_float2 (_S5734);
    float _S5736 = s_primal_ctx_min_0(_S5731, _S5734);
    float _S5737 = s_primal_ctx_max_0(_S5732, _S5734);
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S5738 = s_primal_ctx_mul_0(R_37, pos_i_5) + t_37;
    float _S5739 = _S5738.z;
    float2  _S5740 = make_float2 (_S5739);
    float _S5741 = s_primal_ctx_min_0(_S5736, _S5739);
    float _S5742 = s_primal_ctx_max_0(_S5737, _S5739);
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S5743 = s_primal_ctx_mul_0(R_37, pos_i_6) + t_37;
    float _S5744 = _S5743.z;
    float2  _S5745 = make_float2 (_S5744);
    float3  _S5746 = pos_4 + make_float3 (0.5f * size_4);
    float3  mean_c_31 = s_primal_ctx_mul_0(R_37, _S5746) + t_37;
    float2  _S5747 = float2 {_S5708.x, _S5708.y};
    float2  _S5748 = _S5747 / make_float2 (_S5709);
    float2  _S5749 = make_float2 (_S5709 * _S5709);
    float u_117 = _S5748.x;
    float v_117 = _S5748.y;
    float r2_112 = u_117 * u_117 + v_117 * v_117;
    float _S5750 = (*dist_coeffs_44)[int(2)] + r2_112 * (*dist_coeffs_44)[int(3)];
    float _S5751 = (*dist_coeffs_44)[int(1)] + r2_112 * _S5750;
    float _S5752 = (*dist_coeffs_44)[int(0)] + r2_112 * _S5751;
    float radial_13 = 1.0f + r2_112 * _S5752;
    float2  _S5753 = make_float2 (radial_13);
    float _S5754 = 2.0f * (*dist_coeffs_44)[int(4)];
    float _S5755 = _S5754 * u_117;
    float _S5756 = 2.0f * u_117;
    float _S5757 = 2.0f * (*dist_coeffs_44)[int(5)];
    float _S5758 = _S5757 * u_117;
    float _S5759 = 2.0f * v_117;
    float2  _S5760 = _S5748 * make_float2 (radial_13) + make_float2 (_S5755 * v_117 + (*dist_coeffs_44)[int(5)] * (r2_112 + _S5756 * u_117) + (*dist_coeffs_44)[int(6)] * r2_112, _S5758 * v_117 + (*dist_coeffs_44)[int(4)] * (r2_112 + _S5759 * v_117) + (*dist_coeffs_44)[int(7)] * r2_112);
    float2  _S5761 = _S5760 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5760.x + (*dist_coeffs_44)[int(9)] * _S5760.y, 0.0f);
    float _S5762 = fx_40 * _S5761.x + cx_35;
    float _S5763 = fy_40 * _S5761.y + cy_35;
    float2  _S5764 = float2 {_S5713.x, _S5713.y};
    float2  _S5765 = _S5764 / make_float2 (_S5714);
    float2  _S5766 = make_float2 (_S5714 * _S5714);
    float u_118 = _S5765.x;
    float v_118 = _S5765.y;
    float r2_113 = u_118 * u_118 + v_118 * v_118;
    float _S5767 = (*dist_coeffs_44)[int(2)] + r2_113 * (*dist_coeffs_44)[int(3)];
    float _S5768 = (*dist_coeffs_44)[int(1)] + r2_113 * _S5767;
    float _S5769 = (*dist_coeffs_44)[int(0)] + r2_113 * _S5768;
    float radial_14 = 1.0f + r2_113 * _S5769;
    float2  _S5770 = make_float2 (radial_14);
    float _S5771 = _S5754 * u_118;
    float _S5772 = 2.0f * u_118;
    float _S5773 = _S5757 * u_118;
    float _S5774 = 2.0f * v_118;
    float2  _S5775 = _S5765 * make_float2 (radial_14) + make_float2 (_S5771 * v_118 + (*dist_coeffs_44)[int(5)] * (r2_113 + _S5772 * u_118) + (*dist_coeffs_44)[int(6)] * r2_113, _S5773 * v_118 + (*dist_coeffs_44)[int(4)] * (r2_113 + _S5774 * v_118) + (*dist_coeffs_44)[int(7)] * r2_113);
    float2  _S5776 = _S5775 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5775.x + (*dist_coeffs_44)[int(9)] * _S5775.y, 0.0f);
    float _S5777 = fx_40 * _S5776.x + cx_35;
    float _S5778 = fy_40 * _S5776.y + cy_35;
    float2  _S5779 = float2 {_S5718.x, _S5718.y};
    float2  _S5780 = _S5779 / make_float2 (_S5719);
    float2  _S5781 = make_float2 (_S5719 * _S5719);
    float u_119 = _S5780.x;
    float v_119 = _S5780.y;
    float r2_114 = u_119 * u_119 + v_119 * v_119;
    float _S5782 = (*dist_coeffs_44)[int(2)] + r2_114 * (*dist_coeffs_44)[int(3)];
    float _S5783 = (*dist_coeffs_44)[int(1)] + r2_114 * _S5782;
    float _S5784 = (*dist_coeffs_44)[int(0)] + r2_114 * _S5783;
    float radial_15 = 1.0f + r2_114 * _S5784;
    float2  _S5785 = make_float2 (radial_15);
    float _S5786 = _S5754 * u_119;
    float _S5787 = 2.0f * u_119;
    float _S5788 = _S5757 * u_119;
    float _S5789 = 2.0f * v_119;
    float2  _S5790 = _S5780 * make_float2 (radial_15) + make_float2 (_S5786 * v_119 + (*dist_coeffs_44)[int(5)] * (r2_114 + _S5787 * u_119) + (*dist_coeffs_44)[int(6)] * r2_114, _S5788 * v_119 + (*dist_coeffs_44)[int(4)] * (r2_114 + _S5789 * v_119) + (*dist_coeffs_44)[int(7)] * r2_114);
    float2  _S5791 = _S5790 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5790.x + (*dist_coeffs_44)[int(9)] * _S5790.y, 0.0f);
    float _S5792 = fx_40 * _S5791.x + cx_35;
    float _S5793 = fy_40 * _S5791.y + cy_35;
    float2  _S5794 = float2 {_S5723.x, _S5723.y};
    float2  _S5795 = _S5794 / make_float2 (_S5724);
    float2  _S5796 = make_float2 (_S5724 * _S5724);
    float u_120 = _S5795.x;
    float v_120 = _S5795.y;
    float r2_115 = u_120 * u_120 + v_120 * v_120;
    float _S5797 = (*dist_coeffs_44)[int(2)] + r2_115 * (*dist_coeffs_44)[int(3)];
    float _S5798 = (*dist_coeffs_44)[int(1)] + r2_115 * _S5797;
    float _S5799 = (*dist_coeffs_44)[int(0)] + r2_115 * _S5798;
    float radial_16 = 1.0f + r2_115 * _S5799;
    float2  _S5800 = make_float2 (radial_16);
    float _S5801 = _S5754 * u_120;
    float _S5802 = 2.0f * u_120;
    float _S5803 = _S5757 * u_120;
    float _S5804 = 2.0f * v_120;
    float2  _S5805 = _S5795 * make_float2 (radial_16) + make_float2 (_S5801 * v_120 + (*dist_coeffs_44)[int(5)] * (r2_115 + _S5802 * u_120) + (*dist_coeffs_44)[int(6)] * r2_115, _S5803 * v_120 + (*dist_coeffs_44)[int(4)] * (r2_115 + _S5804 * v_120) + (*dist_coeffs_44)[int(7)] * r2_115);
    float2  _S5806 = _S5805 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5805.x + (*dist_coeffs_44)[int(9)] * _S5805.y, 0.0f);
    float _S5807 = fx_40 * _S5806.x + cx_35;
    float _S5808 = fy_40 * _S5806.y + cy_35;
    float2  _S5809 = float2 {_S5728.x, _S5728.y};
    float2  _S5810 = _S5809 / make_float2 (_S5729);
    float2  _S5811 = make_float2 (_S5729 * _S5729);
    float u_121 = _S5810.x;
    float v_121 = _S5810.y;
    float r2_116 = u_121 * u_121 + v_121 * v_121;
    float _S5812 = (*dist_coeffs_44)[int(2)] + r2_116 * (*dist_coeffs_44)[int(3)];
    float _S5813 = (*dist_coeffs_44)[int(1)] + r2_116 * _S5812;
    float _S5814 = (*dist_coeffs_44)[int(0)] + r2_116 * _S5813;
    float radial_17 = 1.0f + r2_116 * _S5814;
    float2  _S5815 = make_float2 (radial_17);
    float _S5816 = _S5754 * u_121;
    float _S5817 = 2.0f * u_121;
    float _S5818 = _S5757 * u_121;
    float _S5819 = 2.0f * v_121;
    float2  _S5820 = _S5810 * make_float2 (radial_17) + make_float2 (_S5816 * v_121 + (*dist_coeffs_44)[int(5)] * (r2_116 + _S5817 * u_121) + (*dist_coeffs_44)[int(6)] * r2_116, _S5818 * v_121 + (*dist_coeffs_44)[int(4)] * (r2_116 + _S5819 * v_121) + (*dist_coeffs_44)[int(7)] * r2_116);
    float2  _S5821 = _S5820 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5820.x + (*dist_coeffs_44)[int(9)] * _S5820.y, 0.0f);
    float _S5822 = fx_40 * _S5821.x + cx_35;
    float _S5823 = fy_40 * _S5821.y + cy_35;
    float2  _S5824 = float2 {_S5733.x, _S5733.y};
    float2  _S5825 = _S5824 / make_float2 (_S5734);
    float2  _S5826 = make_float2 (_S5734 * _S5734);
    float u_122 = _S5825.x;
    float v_122 = _S5825.y;
    float r2_117 = u_122 * u_122 + v_122 * v_122;
    float _S5827 = (*dist_coeffs_44)[int(2)] + r2_117 * (*dist_coeffs_44)[int(3)];
    float _S5828 = (*dist_coeffs_44)[int(1)] + r2_117 * _S5827;
    float _S5829 = (*dist_coeffs_44)[int(0)] + r2_117 * _S5828;
    float radial_18 = 1.0f + r2_117 * _S5829;
    float2  _S5830 = make_float2 (radial_18);
    float _S5831 = _S5754 * u_122;
    float _S5832 = 2.0f * u_122;
    float _S5833 = _S5757 * u_122;
    float _S5834 = 2.0f * v_122;
    float2  _S5835 = _S5825 * make_float2 (radial_18) + make_float2 (_S5831 * v_122 + (*dist_coeffs_44)[int(5)] * (r2_117 + _S5832 * u_122) + (*dist_coeffs_44)[int(6)] * r2_117, _S5833 * v_122 + (*dist_coeffs_44)[int(4)] * (r2_117 + _S5834 * v_122) + (*dist_coeffs_44)[int(7)] * r2_117);
    float2  _S5836 = _S5835 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5835.x + (*dist_coeffs_44)[int(9)] * _S5835.y, 0.0f);
    float _S5837 = fx_40 * _S5836.x + cx_35;
    float _S5838 = fy_40 * _S5836.y + cy_35;
    float2  _S5839 = float2 {_S5738.x, _S5738.y};
    float2  _S5840 = _S5839 / make_float2 (_S5739);
    float2  _S5841 = make_float2 (_S5739 * _S5739);
    float u_123 = _S5840.x;
    float v_123 = _S5840.y;
    float r2_118 = u_123 * u_123 + v_123 * v_123;
    float _S5842 = (*dist_coeffs_44)[int(2)] + r2_118 * (*dist_coeffs_44)[int(3)];
    float _S5843 = (*dist_coeffs_44)[int(1)] + r2_118 * _S5842;
    float _S5844 = (*dist_coeffs_44)[int(0)] + r2_118 * _S5843;
    float radial_19 = 1.0f + r2_118 * _S5844;
    float2  _S5845 = make_float2 (radial_19);
    float _S5846 = _S5754 * u_123;
    float _S5847 = 2.0f * u_123;
    float _S5848 = _S5757 * u_123;
    float _S5849 = 2.0f * v_123;
    float2  _S5850 = _S5840 * make_float2 (radial_19) + make_float2 (_S5846 * v_123 + (*dist_coeffs_44)[int(5)] * (r2_118 + _S5847 * u_123) + (*dist_coeffs_44)[int(6)] * r2_118, _S5848 * v_123 + (*dist_coeffs_44)[int(4)] * (r2_118 + _S5849 * v_123) + (*dist_coeffs_44)[int(7)] * r2_118);
    float2  _S5851 = _S5850 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5850.x + (*dist_coeffs_44)[int(9)] * _S5850.y, 0.0f);
    float _S5852 = fx_40 * _S5851.x + cx_35;
    float _S5853 = fy_40 * _S5851.y + cy_35;
    float2  _S5854 = float2 {_S5743.x, _S5743.y};
    float2  _S5855 = _S5854 / make_float2 (_S5744);
    float2  _S5856 = make_float2 (_S5744 * _S5744);
    float u_124 = _S5855.x;
    float v_124 = _S5855.y;
    float r2_119 = u_124 * u_124 + v_124 * v_124;
    float _S5857 = (*dist_coeffs_44)[int(2)] + r2_119 * (*dist_coeffs_44)[int(3)];
    float _S5858 = (*dist_coeffs_44)[int(1)] + r2_119 * _S5857;
    float _S5859 = (*dist_coeffs_44)[int(0)] + r2_119 * _S5858;
    float radial_20 = 1.0f + r2_119 * _S5859;
    float2  _S5860 = make_float2 (radial_20);
    float _S5861 = _S5754 * u_124;
    float _S5862 = 2.0f * u_124;
    float _S5863 = _S5757 * u_124;
    float _S5864 = 2.0f * v_124;
    float2  _S5865 = _S5855 * make_float2 (radial_20) + make_float2 (_S5861 * v_124 + (*dist_coeffs_44)[int(5)] * (r2_119 + _S5862 * u_124) + (*dist_coeffs_44)[int(6)] * r2_119, _S5863 * v_124 + (*dist_coeffs_44)[int(4)] * (r2_119 + _S5864 * v_124) + (*dist_coeffs_44)[int(7)] * r2_119);
    float2  _S5866 = _S5865 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5865.x + (*dist_coeffs_44)[int(9)] * _S5865.y, 0.0f);
    float _S5867 = fx_40 * _S5866.x + cx_35;
    float _S5868 = fy_40 * _S5866.y + cy_35;
    float _S5869 = s_primal_ctx_max_0(_S5762, _S5777);
    float _S5870 = s_primal_ctx_min_0(_S5762, _S5777);
    float _S5871 = s_primal_ctx_max_0(_S5763, _S5778);
    float _S5872 = s_primal_ctx_min_0(_S5763, _S5778);
    float _S5873 = s_primal_ctx_max_0(_S5869, _S5792);
    float _S5874 = s_primal_ctx_min_0(_S5870, _S5792);
    float _S5875 = s_primal_ctx_max_0(_S5871, _S5793);
    float _S5876 = s_primal_ctx_min_0(_S5872, _S5793);
    float _S5877 = s_primal_ctx_max_0(_S5873, _S5807);
    float _S5878 = s_primal_ctx_min_0(_S5874, _S5807);
    float _S5879 = s_primal_ctx_max_0(_S5875, _S5808);
    float _S5880 = s_primal_ctx_min_0(_S5876, _S5808);
    float _S5881 = s_primal_ctx_max_0(_S5877, _S5822);
    float _S5882 = s_primal_ctx_min_0(_S5878, _S5822);
    float _S5883 = s_primal_ctx_max_0(_S5879, _S5823);
    float _S5884 = s_primal_ctx_min_0(_S5880, _S5823);
    float _S5885 = s_primal_ctx_max_0(_S5881, _S5837);
    float _S5886 = s_primal_ctx_min_0(_S5882, _S5837);
    float _S5887 = s_primal_ctx_max_0(_S5883, _S5838);
    float _S5888 = s_primal_ctx_min_0(_S5884, _S5838);
    float _S5889 = s_primal_ctx_max_0(_S5885, _S5852);
    float _S5890 = s_primal_ctx_min_0(_S5886, _S5852);
    float _S5891 = s_primal_ctx_max_0(_S5887, _S5853);
    float _S5892 = s_primal_ctx_min_0(_S5888, _S5853);
    float _S5893 = s_primal_ctx_dot_0(mean_c_31, mean_c_31) + 9.99999997475242708e-07f;
    Matrix<float, 3, 3>  _S5894 = transpose_0(R_37);
    float3  _S5895 = mean_c_31 - - s_primal_ctx_mul_0(_S5894, t_37);
    float _S5896 = _S5895.x;
    float _S5897 = _S5895.y;
    float _S5898 = _S5895.z;
    float _S5899 = _S5896 * _S5896 + _S5897 * _S5897 + _S5898 * _S5898;
    float _S5900 = s_primal_ctx_sqrt_0(_S5899);
    float x_67 = _S5896 / _S5900;
    float3  _S5901 = make_float3 (x_67);
    float _S5902 = _S5900 * _S5900;
    float y_34 = _S5897 / _S5900;
    float z_31 = _S5898 / _S5900;
    float3  _S5903 = make_float3 (z_31);
    float _S5904 = - y_34;
    float3  _S5905 = make_float3 (_S5904);
    float z2_64 = z_31 * z_31;
    float fTmp0B_31 = -1.09254848957061768f * z_31;
    float fC1_31 = x_67 * x_67 - y_34 * y_34;
    float _S5906 = 2.0f * x_67;
    float fS1_31 = _S5906 * y_34;
    float pSH6_9 = 0.94617468118667603f * z2_64 - 0.31539157032966614f;
    float3  _S5907 = make_float3 (pSH6_9);
    float pSH7_9 = fTmp0B_31 * x_67;
    float3  _S5908 = make_float3 (pSH7_9);
    float pSH5_9 = fTmp0B_31 * y_34;
    float3  _S5909 = make_float3 (pSH5_9);
    float pSH8_9 = 0.54627424478530884f * fC1_31;
    float3  _S5910 = make_float3 (pSH8_9);
    float pSH4_9 = 0.54627424478530884f * fS1_31;
    float3  _S5911 = make_float3 (pSH4_9);
    float fTmp0C_31 = -2.28522896766662598f * z2_64 + 0.4570457935333252f;
    float fTmp1B_31 = 1.44530570507049561f * z_31;
    float _S5912 = 1.86588168144226074f * z2_64 - 1.11952900886535645f;
    float pSH12_9 = z_31 * _S5912;
    float3  _S5913 = make_float3 (pSH12_9);
    float pSH13_9 = fTmp0C_31 * x_67;
    float3  _S5914 = make_float3 (pSH13_9);
    float pSH11_9 = fTmp0C_31 * y_34;
    float3  _S5915 = make_float3 (pSH11_9);
    float pSH14_9 = fTmp1B_31 * fC1_31;
    float3  _S5916 = make_float3 (pSH14_9);
    float pSH10_9 = fTmp1B_31 * fS1_31;
    float3  _S5917 = make_float3 (pSH10_9);
    float pSH15_9 = -0.59004360437393188f * (x_67 * fC1_31 - y_34 * fS1_31);
    float3  _S5918 = make_float3 (pSH15_9);
    float pSH9_9 = -0.59004360437393188f * (x_67 * fS1_31 + y_34 * fC1_31);
    float3  _S5919 = make_float3 (pSH9_9);
    float3  _S5920 = make_float3 (0.0f);
    float3  _S5921 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5922;
    (&_S5922)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_31)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S5904) * (*sh_coeffs_31)[int(1)] + make_float3 (z_31) * (*sh_coeffs_31)[int(2)] - make_float3 (x_67) * (*sh_coeffs_31)[int(3)]) + (make_float3 (pSH4_9) * (*sh_coeffs_31)[int(4)] + make_float3 (pSH5_9) * (*sh_coeffs_31)[int(5)] + make_float3 (pSH6_9) * (*sh_coeffs_31)[int(6)] + make_float3 (pSH7_9) * (*sh_coeffs_31)[int(7)] + make_float3 (pSH8_9) * (*sh_coeffs_31)[int(8)]) + (make_float3 (pSH9_9) * (*sh_coeffs_31)[int(9)] + make_float3 (pSH10_9) * (*sh_coeffs_31)[int(10)] + make_float3 (pSH11_9) * (*sh_coeffs_31)[int(11)] + make_float3 (pSH12_9) * (*sh_coeffs_31)[int(12)] + make_float3 (pSH13_9) * (*sh_coeffs_31)[int(13)] + make_float3 (pSH14_9) * (*sh_coeffs_31)[int(14)] + make_float3 (pSH15_9) * (*sh_coeffs_31)[int(15)]) + make_float3 (0.5f);
    (&_S5922)->differential_0 = _S5921;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5923;
    (&_S5923)->primal_0 = _S5920;
    (&_S5923)->differential_0 = _S5921;
    s_bwd_prop_max_0(&_S5922, &_S5923, v_rgb_8);
    float3  _S5924 = _S5918 * _S5922.differential_0;
    float3  _S5925 = (*sh_coeffs_31)[int(15)] * _S5922.differential_0;
    float3  _S5926 = _S5916 * _S5922.differential_0;
    float3  _S5927 = (*sh_coeffs_31)[int(14)] * _S5922.differential_0;
    float3  _S5928 = _S5914 * _S5922.differential_0;
    float3  _S5929 = (*sh_coeffs_31)[int(13)] * _S5922.differential_0;
    float3  _S5930 = _S5913 * _S5922.differential_0;
    float3  _S5931 = (*sh_coeffs_31)[int(12)] * _S5922.differential_0;
    float3  _S5932 = _S5915 * _S5922.differential_0;
    float3  _S5933 = (*sh_coeffs_31)[int(11)] * _S5922.differential_0;
    float3  _S5934 = _S5917 * _S5922.differential_0;
    float3  _S5935 = (*sh_coeffs_31)[int(10)] * _S5922.differential_0;
    float3  _S5936 = _S5919 * _S5922.differential_0;
    float3  _S5937 = (*sh_coeffs_31)[int(9)] * _S5922.differential_0;
    float s_diff_fS2_T_9 = -0.59004360437393188f * (_S5937.x + _S5937.y + _S5937.z);
    float s_diff_fC2_T_9 = -0.59004360437393188f * (_S5925.x + _S5925.y + _S5925.z);
    float _S5938 = _S5935.x + _S5935.y + _S5935.z;
    float _S5939 = _S5927.x + _S5927.y + _S5927.z;
    float _S5940 = _S5933.x + _S5933.y + _S5933.z;
    float _S5941 = _S5929.x + _S5929.y + _S5929.z;
    float _S5942 = _S5931.x + _S5931.y + _S5931.z;
    float _S5943 = - s_diff_fC2_T_9;
    float3  _S5944 = _S5910 * _S5922.differential_0;
    float3  _S5945 = (*sh_coeffs_31)[int(8)] * _S5922.differential_0;
    float3  _S5946 = _S5908 * _S5922.differential_0;
    float3  _S5947 = (*sh_coeffs_31)[int(7)] * _S5922.differential_0;
    float3  _S5948 = _S5907 * _S5922.differential_0;
    float3  _S5949 = (*sh_coeffs_31)[int(6)] * _S5922.differential_0;
    float3  _S5950 = _S5909 * _S5922.differential_0;
    float3  _S5951 = (*sh_coeffs_31)[int(5)] * _S5922.differential_0;
    float3  _S5952 = _S5911 * _S5922.differential_0;
    float3  _S5953 = (*sh_coeffs_31)[int(4)] * _S5922.differential_0;
    float _S5954 = _S5951.x + _S5951.y + _S5951.z;
    float _S5955 = _S5947.x + _S5947.y + _S5947.z;
    float _S5956 = fTmp1B_31 * _S5938 + x_67 * s_diff_fS2_T_9 + y_34 * _S5943 + 0.54627424478530884f * (_S5953.x + _S5953.y + _S5953.z);
    float _S5957 = fTmp1B_31 * _S5939 + y_34 * s_diff_fS2_T_9 + x_67 * s_diff_fC2_T_9 + 0.54627424478530884f * (_S5945.x + _S5945.y + _S5945.z);
    float _S5958 = y_34 * - _S5957;
    float _S5959 = x_67 * _S5957;
    float _S5960 = z_31 * (1.86588168144226074f * (z_31 * _S5942) + -2.28522896766662598f * (y_34 * _S5940 + x_67 * _S5941) + 0.94617468118667603f * (_S5949.x + _S5949.y + _S5949.z));
    float3  _S5961 = make_float3 (0.48860251903533936f) * _S5922.differential_0;
    float3  _S5962 = - _S5961;
    float3  _S5963 = _S5901 * _S5962;
    float3  _S5964 = (*sh_coeffs_31)[int(3)] * _S5962;
    float3  _S5965 = _S5903 * _S5961;
    float3  _S5966 = (*sh_coeffs_31)[int(2)] * _S5961;
    float3  _S5967 = _S5905 * _S5961;
    float3  _S5968 = (*sh_coeffs_31)[int(1)] * _S5961;
    float _S5969 = (_S5912 * _S5942 + 1.44530570507049561f * (fS1_31 * _S5938 + fC1_31 * _S5939) + -1.09254848957061768f * (y_34 * _S5954 + x_67 * _S5955) + _S5960 + _S5960 + _S5966.x + _S5966.y + _S5966.z) / _S5902;
    float _S5970 = _S5900 * _S5969;
    float _S5971 = (fTmp0C_31 * _S5940 + fC1_31 * s_diff_fS2_T_9 + fS1_31 * _S5943 + fTmp0B_31 * _S5954 + _S5906 * _S5956 + _S5958 + _S5958 + - (_S5968.x + _S5968.y + _S5968.z)) / _S5902;
    float _S5972 = _S5900 * _S5971;
    float _S5973 = (fTmp0C_31 * _S5941 + fS1_31 * s_diff_fS2_T_9 + fC1_31 * s_diff_fC2_T_9 + fTmp0B_31 * _S5955 + 2.0f * (y_34 * _S5956) + _S5959 + _S5959 + _S5964.x + _S5964.y + _S5964.z) / _S5902;
    float _S5974 = _S5900 * _S5973;
    float _S5975 = _S5898 * - _S5969 + _S5897 * - _S5971 + _S5896 * - _S5973;
    DiffPair_float_0 _S5976;
    (&_S5976)->primal_0 = _S5899;
    (&_S5976)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5976, _S5975);
    float _S5977 = _S5898 * _S5976.differential_0;
    float _S5978 = _S5897 * _S5976.differential_0;
    float _S5979 = _S5896 * _S5976.differential_0;
    float3  _S5980 = make_float3 (0.282094806432724f) * _S5922.differential_0;
    float3  _S5981 = make_float3 (_S5974 + _S5979 + _S5979, _S5972 + _S5978 + _S5978, _S5970 + _S5977 + _S5977);
    float3  _S5982 = - - _S5981;
    Matrix<float, 3, 3>  _S5983 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5984;
    (&_S5984)->primal_0 = _S5894;
    (&_S5984)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5985;
    (&_S5985)->primal_0 = t_37;
    (&_S5985)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S5984, &_S5985, _S5982);
    Matrix<float, 3, 3>  _S5986 = transpose_0(_S5984.differential_0);
    DiffPair_float_0 _S5987;
    (&_S5987)->primal_0 = _S5893;
    (&_S5987)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5987, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5988;
    (&_S5988)->primal_0 = mean_c_31;
    (&_S5988)->differential_0 = _S5921;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5989;
    (&_S5989)->primal_0 = mean_c_31;
    (&_S5989)->differential_0 = _S5921;
    s_bwd_prop_dot_0(&_S5988, &_S5989, _S5987.differential_0);
    DiffPair_float_0 _S5990;
    (&_S5990)->primal_0 = _S5892;
    (&_S5990)->differential_0 = 0.0f;
    DiffPair_float_0 _S5991;
    (&_S5991)->primal_0 = _S5868;
    (&_S5991)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5990, &_S5991, 0.0f);
    DiffPair_float_0 _S5992;
    (&_S5992)->primal_0 = _S5891;
    (&_S5992)->differential_0 = 0.0f;
    DiffPair_float_0 _S5993;
    (&_S5993)->primal_0 = _S5868;
    (&_S5993)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5992, &_S5993, 0.0f);
    DiffPair_float_0 _S5994;
    (&_S5994)->primal_0 = _S5890;
    (&_S5994)->differential_0 = 0.0f;
    DiffPair_float_0 _S5995;
    (&_S5995)->primal_0 = _S5867;
    (&_S5995)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5994, &_S5995, 0.0f);
    DiffPair_float_0 _S5996;
    (&_S5996)->primal_0 = _S5889;
    (&_S5996)->differential_0 = 0.0f;
    DiffPair_float_0 _S5997;
    (&_S5997)->primal_0 = _S5867;
    (&_S5997)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5996, &_S5997, 0.0f);
    DiffPair_float_0 _S5998;
    (&_S5998)->primal_0 = _S5888;
    (&_S5998)->differential_0 = 0.0f;
    DiffPair_float_0 _S5999;
    (&_S5999)->primal_0 = _S5853;
    (&_S5999)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5998, &_S5999, _S5990.differential_0);
    DiffPair_float_0 _S6000;
    (&_S6000)->primal_0 = _S5887;
    (&_S6000)->differential_0 = 0.0f;
    DiffPair_float_0 _S6001;
    (&_S6001)->primal_0 = _S5853;
    (&_S6001)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6000, &_S6001, _S5992.differential_0);
    DiffPair_float_0 _S6002;
    (&_S6002)->primal_0 = _S5886;
    (&_S6002)->differential_0 = 0.0f;
    DiffPair_float_0 _S6003;
    (&_S6003)->primal_0 = _S5852;
    (&_S6003)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6002, &_S6003, _S5994.differential_0);
    DiffPair_float_0 _S6004;
    (&_S6004)->primal_0 = _S5885;
    (&_S6004)->differential_0 = 0.0f;
    DiffPair_float_0 _S6005;
    (&_S6005)->primal_0 = _S5852;
    (&_S6005)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6004, &_S6005, _S5996.differential_0);
    DiffPair_float_0 _S6006;
    (&_S6006)->primal_0 = _S5884;
    (&_S6006)->differential_0 = 0.0f;
    DiffPair_float_0 _S6007;
    (&_S6007)->primal_0 = _S5838;
    (&_S6007)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6006, &_S6007, _S5998.differential_0);
    DiffPair_float_0 _S6008;
    (&_S6008)->primal_0 = _S5883;
    (&_S6008)->differential_0 = 0.0f;
    DiffPair_float_0 _S6009;
    (&_S6009)->primal_0 = _S5838;
    (&_S6009)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6008, &_S6009, _S6000.differential_0);
    DiffPair_float_0 _S6010;
    (&_S6010)->primal_0 = _S5882;
    (&_S6010)->differential_0 = 0.0f;
    DiffPair_float_0 _S6011;
    (&_S6011)->primal_0 = _S5837;
    (&_S6011)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6010, &_S6011, _S6002.differential_0);
    DiffPair_float_0 _S6012;
    (&_S6012)->primal_0 = _S5881;
    (&_S6012)->differential_0 = 0.0f;
    DiffPair_float_0 _S6013;
    (&_S6013)->primal_0 = _S5837;
    (&_S6013)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6012, &_S6013, _S6004.differential_0);
    DiffPair_float_0 _S6014;
    (&_S6014)->primal_0 = _S5880;
    (&_S6014)->differential_0 = 0.0f;
    DiffPair_float_0 _S6015;
    (&_S6015)->primal_0 = _S5823;
    (&_S6015)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6014, &_S6015, _S6006.differential_0);
    DiffPair_float_0 _S6016;
    (&_S6016)->primal_0 = _S5879;
    (&_S6016)->differential_0 = 0.0f;
    DiffPair_float_0 _S6017;
    (&_S6017)->primal_0 = _S5823;
    (&_S6017)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6016, &_S6017, _S6008.differential_0);
    DiffPair_float_0 _S6018;
    (&_S6018)->primal_0 = _S5878;
    (&_S6018)->differential_0 = 0.0f;
    DiffPair_float_0 _S6019;
    (&_S6019)->primal_0 = _S5822;
    (&_S6019)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6018, &_S6019, _S6010.differential_0);
    DiffPair_float_0 _S6020;
    (&_S6020)->primal_0 = _S5877;
    (&_S6020)->differential_0 = 0.0f;
    DiffPair_float_0 _S6021;
    (&_S6021)->primal_0 = _S5822;
    (&_S6021)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6020, &_S6021, _S6012.differential_0);
    DiffPair_float_0 _S6022;
    (&_S6022)->primal_0 = _S5876;
    (&_S6022)->differential_0 = 0.0f;
    DiffPair_float_0 _S6023;
    (&_S6023)->primal_0 = _S5808;
    (&_S6023)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6022, &_S6023, _S6014.differential_0);
    DiffPair_float_0 _S6024;
    (&_S6024)->primal_0 = _S5875;
    (&_S6024)->differential_0 = 0.0f;
    DiffPair_float_0 _S6025;
    (&_S6025)->primal_0 = _S5808;
    (&_S6025)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6024, &_S6025, _S6016.differential_0);
    DiffPair_float_0 _S6026;
    (&_S6026)->primal_0 = _S5874;
    (&_S6026)->differential_0 = 0.0f;
    DiffPair_float_0 _S6027;
    (&_S6027)->primal_0 = _S5807;
    (&_S6027)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6026, &_S6027, _S6018.differential_0);
    DiffPair_float_0 _S6028;
    (&_S6028)->primal_0 = _S5873;
    (&_S6028)->differential_0 = 0.0f;
    DiffPair_float_0 _S6029;
    (&_S6029)->primal_0 = _S5807;
    (&_S6029)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6028, &_S6029, _S6020.differential_0);
    DiffPair_float_0 _S6030;
    (&_S6030)->primal_0 = _S5872;
    (&_S6030)->differential_0 = 0.0f;
    DiffPair_float_0 _S6031;
    (&_S6031)->primal_0 = _S5793;
    (&_S6031)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6030, &_S6031, _S6022.differential_0);
    DiffPair_float_0 _S6032;
    (&_S6032)->primal_0 = _S5871;
    (&_S6032)->differential_0 = 0.0f;
    DiffPair_float_0 _S6033;
    (&_S6033)->primal_0 = _S5793;
    (&_S6033)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6032, &_S6033, _S6024.differential_0);
    DiffPair_float_0 _S6034;
    (&_S6034)->primal_0 = _S5870;
    (&_S6034)->differential_0 = 0.0f;
    DiffPair_float_0 _S6035;
    (&_S6035)->primal_0 = _S5792;
    (&_S6035)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6034, &_S6035, _S6026.differential_0);
    DiffPair_float_0 _S6036;
    (&_S6036)->primal_0 = _S5869;
    (&_S6036)->differential_0 = 0.0f;
    DiffPair_float_0 _S6037;
    (&_S6037)->primal_0 = _S5792;
    (&_S6037)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6036, &_S6037, _S6028.differential_0);
    DiffPair_float_0 _S6038;
    (&_S6038)->primal_0 = _S5763;
    (&_S6038)->differential_0 = 0.0f;
    DiffPair_float_0 _S6039;
    (&_S6039)->primal_0 = _S5778;
    (&_S6039)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6038, &_S6039, _S6030.differential_0);
    DiffPair_float_0 _S6040;
    (&_S6040)->primal_0 = _S5763;
    (&_S6040)->differential_0 = 0.0f;
    DiffPair_float_0 _S6041;
    (&_S6041)->primal_0 = _S5778;
    (&_S6041)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6040, &_S6041, _S6032.differential_0);
    DiffPair_float_0 _S6042;
    (&_S6042)->primal_0 = _S5762;
    (&_S6042)->differential_0 = 0.0f;
    DiffPair_float_0 _S6043;
    (&_S6043)->primal_0 = _S5777;
    (&_S6043)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6042, &_S6043, _S6034.differential_0);
    DiffPair_float_0 _S6044;
    (&_S6044)->primal_0 = _S5762;
    (&_S6044)->differential_0 = 0.0f;
    DiffPair_float_0 _S6045;
    (&_S6045)->primal_0 = _S5777;
    (&_S6045)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6044, &_S6045, _S6036.differential_0);
    float _S6046 = fx_40 * (_S5995.differential_0 + _S5997.differential_0);
    float2  _S6047 = make_float2 (_S6046, fy_40 * (_S5991.differential_0 + _S5993.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6046, (*dist_coeffs_44)[int(9)] * _S6046);
    float2  _S6048 = _S5855 * _S6047;
    float _S6049 = (*dist_coeffs_44)[int(4)] * _S6047.y;
    float _S6050 = (*dist_coeffs_44)[int(5)] * _S6047.x;
    float _S6051 = _S6048.x + _S6048.y;
    float _S6052 = r2_119 * _S6051;
    float _S6053 = r2_119 * _S6052;
    float _S6054 = (*dist_coeffs_44)[int(7)] * _S6047.y + _S6049 + (*dist_coeffs_44)[int(6)] * _S6047.x + _S6050 + _S5859 * _S6051 + _S5858 * _S6052 + _S5857 * _S6053 + (*dist_coeffs_44)[int(3)] * (r2_119 * _S6053);
    float _S6055 = v_124 * _S6054;
    float _S6056 = u_124 * _S6054;
    float2  _S6057 = (_S5860 * _S6047 + make_float2 (_S5757 * (v_124 * _S6047.y) + _S5862 * _S6050 + 2.0f * (u_124 * _S6050) + _S5754 * (v_124 * _S6047.x) + _S6056 + _S6056, _S5864 * _S6049 + 2.0f * (v_124 * _S6049) + _S5863 * _S6047.y + _S5861 * _S6047.x + _S6055 + _S6055)) / _S5856;
    float2  _S6058 = _S5854 * - _S6057;
    float2  _S6059 = _S5745 * _S6057;
    float _S6060 = fx_40 * (_S6003.differential_0 + _S6005.differential_0);
    float2  _S6061 = make_float2 (_S6060, fy_40 * (_S5999.differential_0 + _S6001.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6060, (*dist_coeffs_44)[int(9)] * _S6060);
    float2  _S6062 = _S5840 * _S6061;
    float _S6063 = (*dist_coeffs_44)[int(4)] * _S6061.y;
    float _S6064 = (*dist_coeffs_44)[int(5)] * _S6061.x;
    float _S6065 = _S6062.x + _S6062.y;
    float _S6066 = r2_118 * _S6065;
    float _S6067 = r2_118 * _S6066;
    float _S6068 = (*dist_coeffs_44)[int(7)] * _S6061.y + _S6063 + (*dist_coeffs_44)[int(6)] * _S6061.x + _S6064 + _S5844 * _S6065 + _S5843 * _S6066 + _S5842 * _S6067 + (*dist_coeffs_44)[int(3)] * (r2_118 * _S6067);
    float _S6069 = v_123 * _S6068;
    float _S6070 = u_123 * _S6068;
    float2  _S6071 = (_S5845 * _S6061 + make_float2 (_S5757 * (v_123 * _S6061.y) + _S5847 * _S6064 + 2.0f * (u_123 * _S6064) + _S5754 * (v_123 * _S6061.x) + _S6070 + _S6070, _S5849 * _S6063 + 2.0f * (v_123 * _S6063) + _S5848 * _S6061.y + _S5846 * _S6061.x + _S6069 + _S6069)) / _S5841;
    float2  _S6072 = _S5839 * - _S6071;
    float2  _S6073 = _S5740 * _S6071;
    float _S6074 = fx_40 * (_S6011.differential_0 + _S6013.differential_0);
    float2  _S6075 = make_float2 (_S6074, fy_40 * (_S6007.differential_0 + _S6009.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6074, (*dist_coeffs_44)[int(9)] * _S6074);
    float2  _S6076 = _S5825 * _S6075;
    float _S6077 = (*dist_coeffs_44)[int(4)] * _S6075.y;
    float _S6078 = (*dist_coeffs_44)[int(5)] * _S6075.x;
    float _S6079 = _S6076.x + _S6076.y;
    float _S6080 = r2_117 * _S6079;
    float _S6081 = r2_117 * _S6080;
    float _S6082 = (*dist_coeffs_44)[int(7)] * _S6075.y + _S6077 + (*dist_coeffs_44)[int(6)] * _S6075.x + _S6078 + _S5829 * _S6079 + _S5828 * _S6080 + _S5827 * _S6081 + (*dist_coeffs_44)[int(3)] * (r2_117 * _S6081);
    float _S6083 = v_122 * _S6082;
    float _S6084 = u_122 * _S6082;
    float2  _S6085 = (_S5830 * _S6075 + make_float2 (_S5757 * (v_122 * _S6075.y) + _S5832 * _S6078 + 2.0f * (u_122 * _S6078) + _S5754 * (v_122 * _S6075.x) + _S6084 + _S6084, _S5834 * _S6077 + 2.0f * (v_122 * _S6077) + _S5833 * _S6075.y + _S5831 * _S6075.x + _S6083 + _S6083)) / _S5826;
    float2  _S6086 = _S5824 * - _S6085;
    float2  _S6087 = _S5735 * _S6085;
    float _S6088 = fx_40 * (_S6019.differential_0 + _S6021.differential_0);
    float2  _S6089 = make_float2 (_S6088, fy_40 * (_S6015.differential_0 + _S6017.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6088, (*dist_coeffs_44)[int(9)] * _S6088);
    float2  _S6090 = _S5810 * _S6089;
    float _S6091 = (*dist_coeffs_44)[int(4)] * _S6089.y;
    float _S6092 = (*dist_coeffs_44)[int(5)] * _S6089.x;
    float _S6093 = _S6090.x + _S6090.y;
    float _S6094 = r2_116 * _S6093;
    float _S6095 = r2_116 * _S6094;
    float _S6096 = (*dist_coeffs_44)[int(7)] * _S6089.y + _S6091 + (*dist_coeffs_44)[int(6)] * _S6089.x + _S6092 + _S5814 * _S6093 + _S5813 * _S6094 + _S5812 * _S6095 + (*dist_coeffs_44)[int(3)] * (r2_116 * _S6095);
    float _S6097 = v_121 * _S6096;
    float _S6098 = u_121 * _S6096;
    float2  _S6099 = (_S5815 * _S6089 + make_float2 (_S5757 * (v_121 * _S6089.y) + _S5817 * _S6092 + 2.0f * (u_121 * _S6092) + _S5754 * (v_121 * _S6089.x) + _S6098 + _S6098, _S5819 * _S6091 + 2.0f * (v_121 * _S6091) + _S5818 * _S6089.y + _S5816 * _S6089.x + _S6097 + _S6097)) / _S5811;
    float2  _S6100 = _S5809 * - _S6099;
    float2  _S6101 = _S5730 * _S6099;
    float _S6102 = fx_40 * (_S6027.differential_0 + _S6029.differential_0);
    float2  _S6103 = make_float2 (_S6102, fy_40 * (_S6023.differential_0 + _S6025.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6102, (*dist_coeffs_44)[int(9)] * _S6102);
    float2  _S6104 = _S5795 * _S6103;
    float _S6105 = (*dist_coeffs_44)[int(4)] * _S6103.y;
    float _S6106 = (*dist_coeffs_44)[int(5)] * _S6103.x;
    float _S6107 = _S6104.x + _S6104.y;
    float _S6108 = r2_115 * _S6107;
    float _S6109 = r2_115 * _S6108;
    float _S6110 = (*dist_coeffs_44)[int(7)] * _S6103.y + _S6105 + (*dist_coeffs_44)[int(6)] * _S6103.x + _S6106 + _S5799 * _S6107 + _S5798 * _S6108 + _S5797 * _S6109 + (*dist_coeffs_44)[int(3)] * (r2_115 * _S6109);
    float _S6111 = v_120 * _S6110;
    float _S6112 = u_120 * _S6110;
    float2  _S6113 = (_S5800 * _S6103 + make_float2 (_S5757 * (v_120 * _S6103.y) + _S5802 * _S6106 + 2.0f * (u_120 * _S6106) + _S5754 * (v_120 * _S6103.x) + _S6112 + _S6112, _S5804 * _S6105 + 2.0f * (v_120 * _S6105) + _S5803 * _S6103.y + _S5801 * _S6103.x + _S6111 + _S6111)) / _S5796;
    float2  _S6114 = _S5794 * - _S6113;
    float2  _S6115 = _S5725 * _S6113;
    float _S6116 = fx_40 * (_S6035.differential_0 + _S6037.differential_0);
    float2  _S6117 = make_float2 (_S6116, fy_40 * (_S6031.differential_0 + _S6033.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6116, (*dist_coeffs_44)[int(9)] * _S6116);
    float2  _S6118 = _S5780 * _S6117;
    float _S6119 = (*dist_coeffs_44)[int(4)] * _S6117.y;
    float _S6120 = (*dist_coeffs_44)[int(5)] * _S6117.x;
    float _S6121 = _S6118.x + _S6118.y;
    float _S6122 = r2_114 * _S6121;
    float _S6123 = r2_114 * _S6122;
    float _S6124 = (*dist_coeffs_44)[int(7)] * _S6117.y + _S6119 + (*dist_coeffs_44)[int(6)] * _S6117.x + _S6120 + _S5784 * _S6121 + _S5783 * _S6122 + _S5782 * _S6123 + (*dist_coeffs_44)[int(3)] * (r2_114 * _S6123);
    float _S6125 = v_119 * _S6124;
    float _S6126 = u_119 * _S6124;
    float2  _S6127 = (_S5785 * _S6117 + make_float2 (_S5757 * (v_119 * _S6117.y) + _S5787 * _S6120 + 2.0f * (u_119 * _S6120) + _S5754 * (v_119 * _S6117.x) + _S6126 + _S6126, _S5789 * _S6119 + 2.0f * (v_119 * _S6119) + _S5788 * _S6117.y + _S5786 * _S6117.x + _S6125 + _S6125)) / _S5781;
    float2  _S6128 = _S5779 * - _S6127;
    float2  _S6129 = _S5720 * _S6127;
    float _S6130 = fx_40 * (_S6043.differential_0 + _S6045.differential_0);
    float2  _S6131 = make_float2 (_S6130, fy_40 * (_S6039.differential_0 + _S6041.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6130, (*dist_coeffs_44)[int(9)] * _S6130);
    float2  _S6132 = _S5765 * _S6131;
    float _S6133 = (*dist_coeffs_44)[int(4)] * _S6131.y;
    float _S6134 = (*dist_coeffs_44)[int(5)] * _S6131.x;
    float _S6135 = _S6132.x + _S6132.y;
    float _S6136 = r2_113 * _S6135;
    float _S6137 = r2_113 * _S6136;
    float _S6138 = (*dist_coeffs_44)[int(7)] * _S6131.y + _S6133 + (*dist_coeffs_44)[int(6)] * _S6131.x + _S6134 + _S5769 * _S6135 + _S5768 * _S6136 + _S5767 * _S6137 + (*dist_coeffs_44)[int(3)] * (r2_113 * _S6137);
    float _S6139 = v_118 * _S6138;
    float _S6140 = u_118 * _S6138;
    float2  _S6141 = (_S5770 * _S6131 + make_float2 (_S5757 * (v_118 * _S6131.y) + _S5772 * _S6134 + 2.0f * (u_118 * _S6134) + _S5754 * (v_118 * _S6131.x) + _S6140 + _S6140, _S5774 * _S6133 + 2.0f * (v_118 * _S6133) + _S5773 * _S6131.y + _S5771 * _S6131.x + _S6139 + _S6139)) / _S5766;
    float2  _S6142 = _S5764 * - _S6141;
    float2  _S6143 = _S5715 * _S6141;
    float _S6144 = fx_40 * (_S6042.differential_0 + _S6044.differential_0);
    float2  _S6145 = make_float2 (_S6144, fy_40 * (_S6038.differential_0 + _S6040.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6144, (*dist_coeffs_44)[int(9)] * _S6144);
    float2  _S6146 = _S5748 * _S6145;
    float _S6147 = (*dist_coeffs_44)[int(4)] * _S6145.y;
    float _S6148 = (*dist_coeffs_44)[int(5)] * _S6145.x;
    float _S6149 = _S6146.x + _S6146.y;
    float _S6150 = r2_112 * _S6149;
    float _S6151 = r2_112 * _S6150;
    float _S6152 = (*dist_coeffs_44)[int(7)] * _S6145.y + _S6147 + (*dist_coeffs_44)[int(6)] * _S6145.x + _S6148 + _S5752 * _S6149 + _S5751 * _S6150 + _S5750 * _S6151 + (*dist_coeffs_44)[int(3)] * (r2_112 * _S6151);
    float _S6153 = v_117 * _S6152;
    float _S6154 = u_117 * _S6152;
    float2  _S6155 = (_S5753 * _S6145 + make_float2 (_S5757 * (v_117 * _S6145.y) + _S5756 * _S6148 + 2.0f * (u_117 * _S6148) + _S5754 * (v_117 * _S6145.x) + _S6154 + _S6154, _S5759 * _S6147 + 2.0f * (v_117 * _S6147) + _S5758 * _S6145.y + _S5755 * _S6145.x + _S6153 + _S6153)) / _S5749;
    float2  _S6156 = _S5747 * - _S6155;
    float2  _S6157 = _S5710 * _S6155;
    float3  _S6158 = _S5981 + _S5989.differential_0 + _S5988.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6159;
    (&_S6159)->primal_0 = R_37;
    (&_S6159)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6160;
    (&_S6160)->primal_0 = _S5746;
    (&_S6160)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6159, &_S6160, _S6158);
    DiffPair_float_0 _S6161;
    (&_S6161)->primal_0 = _S5742;
    (&_S6161)->differential_0 = 0.0f;
    DiffPair_float_0 _S6162;
    (&_S6162)->primal_0 = _S5744;
    (&_S6162)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6161, &_S6162, 0.0f);
    DiffPair_float_0 _S6163;
    (&_S6163)->primal_0 = _S5741;
    (&_S6163)->differential_0 = 0.0f;
    DiffPair_float_0 _S6164;
    (&_S6164)->primal_0 = _S5744;
    (&_S6164)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6163, &_S6164, 0.0f);
    float3  _S6165 = make_float3 (_S6059.x, _S6059.y, _S6162.differential_0 + _S6164.differential_0 + _S6058.x + _S6058.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6166;
    (&_S6166)->primal_0 = R_37;
    (&_S6166)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6167;
    (&_S6167)->primal_0 = pos_i_6;
    (&_S6167)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6166, &_S6167, _S6165);
    DiffPair_float_0 _S6168;
    (&_S6168)->primal_0 = _S5737;
    (&_S6168)->differential_0 = 0.0f;
    DiffPair_float_0 _S6169;
    (&_S6169)->primal_0 = _S5739;
    (&_S6169)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6168, &_S6169, _S6161.differential_0);
    DiffPair_float_0 _S6170;
    (&_S6170)->primal_0 = _S5736;
    (&_S6170)->differential_0 = 0.0f;
    DiffPair_float_0 _S6171;
    (&_S6171)->primal_0 = _S5739;
    (&_S6171)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6170, &_S6171, _S6163.differential_0);
    float3  _S6172 = make_float3 (_S6073.x, _S6073.y, _S6169.differential_0 + _S6171.differential_0 + _S6072.x + _S6072.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6173;
    (&_S6173)->primal_0 = R_37;
    (&_S6173)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6174;
    (&_S6174)->primal_0 = pos_i_5;
    (&_S6174)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6173, &_S6174, _S6172);
    DiffPair_float_0 _S6175;
    (&_S6175)->primal_0 = _S5732;
    (&_S6175)->differential_0 = 0.0f;
    DiffPair_float_0 _S6176;
    (&_S6176)->primal_0 = _S5734;
    (&_S6176)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6175, &_S6176, _S6168.differential_0);
    DiffPair_float_0 _S6177;
    (&_S6177)->primal_0 = _S5731;
    (&_S6177)->differential_0 = 0.0f;
    DiffPair_float_0 _S6178;
    (&_S6178)->primal_0 = _S5734;
    (&_S6178)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6177, &_S6178, _S6170.differential_0);
    float3  _S6179 = make_float3 (_S6087.x, _S6087.y, _S6176.differential_0 + _S6178.differential_0 + _S6086.x + _S6086.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6180;
    (&_S6180)->primal_0 = R_37;
    (&_S6180)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6181;
    (&_S6181)->primal_0 = pos_i_4;
    (&_S6181)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6180, &_S6181, _S6179);
    DiffPair_float_0 _S6182;
    (&_S6182)->primal_0 = _S5727;
    (&_S6182)->differential_0 = 0.0f;
    DiffPair_float_0 _S6183;
    (&_S6183)->primal_0 = _S5729;
    (&_S6183)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6182, &_S6183, _S6175.differential_0);
    DiffPair_float_0 _S6184;
    (&_S6184)->primal_0 = _S5726;
    (&_S6184)->differential_0 = 0.0f;
    DiffPair_float_0 _S6185;
    (&_S6185)->primal_0 = _S5729;
    (&_S6185)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6184, &_S6185, _S6177.differential_0);
    float3  _S6186 = make_float3 (_S6101.x, _S6101.y, _S6183.differential_0 + _S6185.differential_0 + _S6100.x + _S6100.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6187;
    (&_S6187)->primal_0 = R_37;
    (&_S6187)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6188;
    (&_S6188)->primal_0 = pos_i_3;
    (&_S6188)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6187, &_S6188, _S6186);
    DiffPair_float_0 _S6189;
    (&_S6189)->primal_0 = _S5722;
    (&_S6189)->differential_0 = 0.0f;
    DiffPair_float_0 _S6190;
    (&_S6190)->primal_0 = _S5724;
    (&_S6190)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6189, &_S6190, _S6182.differential_0);
    DiffPair_float_0 _S6191;
    (&_S6191)->primal_0 = _S5721;
    (&_S6191)->differential_0 = 0.0f;
    DiffPair_float_0 _S6192;
    (&_S6192)->primal_0 = _S5724;
    (&_S6192)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6191, &_S6192, _S6184.differential_0);
    float3  _S6193 = make_float3 (_S6115.x, _S6115.y, _S6190.differential_0 + _S6192.differential_0 + _S6114.x + _S6114.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6194;
    (&_S6194)->primal_0 = R_37;
    (&_S6194)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6195;
    (&_S6195)->primal_0 = pos_i_2;
    (&_S6195)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6194, &_S6195, _S6193);
    DiffPair_float_0 _S6196;
    (&_S6196)->primal_0 = _S5717;
    (&_S6196)->differential_0 = 0.0f;
    DiffPair_float_0 _S6197;
    (&_S6197)->primal_0 = _S5719;
    (&_S6197)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6196, &_S6197, _S6189.differential_0);
    DiffPair_float_0 _S6198;
    (&_S6198)->primal_0 = _S5716;
    (&_S6198)->differential_0 = 0.0f;
    DiffPair_float_0 _S6199;
    (&_S6199)->primal_0 = _S5719;
    (&_S6199)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6198, &_S6199, _S6191.differential_0);
    float3  _S6200 = make_float3 (_S6129.x, _S6129.y, _S6197.differential_0 + _S6199.differential_0 + _S6128.x + _S6128.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6201;
    (&_S6201)->primal_0 = R_37;
    (&_S6201)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6202;
    (&_S6202)->primal_0 = pos_i_1;
    (&_S6202)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6201, &_S6202, _S6200);
    DiffPair_float_0 _S6203;
    (&_S6203)->primal_0 = _S5712;
    (&_S6203)->differential_0 = 0.0f;
    DiffPair_float_0 _S6204;
    (&_S6204)->primal_0 = _S5714;
    (&_S6204)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6203, &_S6204, _S6196.differential_0);
    DiffPair_float_0 _S6205;
    (&_S6205)->primal_0 = _S5711;
    (&_S6205)->differential_0 = 0.0f;
    DiffPair_float_0 _S6206;
    (&_S6206)->primal_0 = _S5714;
    (&_S6206)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6205, &_S6206, _S6198.differential_0);
    float3  _S6207 = make_float3 (_S6143.x, _S6143.y, _S6204.differential_0 + _S6206.differential_0 + _S6142.x + _S6142.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6208;
    (&_S6208)->primal_0 = R_37;
    (&_S6208)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6209;
    (&_S6209)->primal_0 = pos_i_0;
    (&_S6209)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6208, &_S6209, _S6207);
    DiffPair_float_0 _S6210;
    (&_S6210)->primal_0 = 0.0f;
    (&_S6210)->differential_0 = 0.0f;
    DiffPair_float_0 _S6211;
    (&_S6211)->primal_0 = _S5709;
    (&_S6211)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6210, &_S6211, _S6203.differential_0);
    DiffPair_float_0 _S6212;
    (&_S6212)->primal_0 = 1.00000001504746622e+30f;
    (&_S6212)->differential_0 = 0.0f;
    DiffPair_float_0 _S6213;
    (&_S6213)->primal_0 = _S5709;
    (&_S6213)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6212, &_S6213, _S6205.differential_0);
    float3  _S6214 = make_float3 (_S6157.x, _S6157.y, _S6211.differential_0 + _S6213.differential_0 + _S6156.x + _S6156.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6215;
    (&_S6215)->primal_0 = R_37;
    (&_S6215)->differential_0 = _S5983;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6216;
    (&_S6216)->primal_0 = pos_4;
    (&_S6216)->differential_0 = _S5921;
    s_bwd_prop_mul_1(&_S6215, &_S6216, _S6214);
    float3  _S6217 = _S5985.differential_0 + _S6158 + _S6165 + _S6172 + _S6179 + _S6186 + _S6193 + _S6200 + _S6207 + _S6214;
    Matrix<float, 3, 3>  _S6218 = _S5986 + _S6159.differential_0 + _S6166.differential_0 + _S6173.differential_0 + _S6180.differential_0 + _S6187.differential_0 + _S6194.differential_0 + _S6201.differential_0 + _S6208.differential_0 + _S6215.differential_0;
    FixedArray<float3 , 16>  _S6219;
    _S6219[int(0)] = _S5921;
    _S6219[int(1)] = _S5921;
    _S6219[int(2)] = _S5921;
    _S6219[int(3)] = _S5921;
    _S6219[int(4)] = _S5921;
    _S6219[int(5)] = _S5921;
    _S6219[int(6)] = _S5921;
    _S6219[int(7)] = _S5921;
    _S6219[int(8)] = _S5921;
    _S6219[int(9)] = _S5921;
    _S6219[int(10)] = _S5921;
    _S6219[int(11)] = _S5921;
    _S6219[int(12)] = _S5921;
    _S6219[int(13)] = _S5921;
    _S6219[int(14)] = _S5921;
    _S6219[int(15)] = _S5921;
    _S6219[int(15)] = _S5924;
    _S6219[int(14)] = _S5926;
    _S6219[int(13)] = _S5928;
    _S6219[int(12)] = _S5930;
    _S6219[int(11)] = _S5932;
    _S6219[int(10)] = _S5934;
    _S6219[int(9)] = _S5936;
    _S6219[int(8)] = _S5944;
    _S6219[int(7)] = _S5946;
    _S6219[int(6)] = _S5948;
    _S6219[int(5)] = _S5950;
    _S6219[int(4)] = _S5952;
    _S6219[int(3)] = _S5963;
    _S6219[int(2)] = _S5965;
    _S6219[int(1)] = _S5967;
    _S6219[int(0)] = _S5980;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_9 = _S6219;
    *v_R_11 = _S6218;
    *v_t_10 = _S6217;
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  * densities_5, FixedArray<float3 , 16>  * sh_coeffs_32, Matrix<float, 3, 3>  R_38, float3  t_38, float fx_41, float fy_41, float cx_36, float cy_36, FixedArray<float, 10>  * dist_coeffs_45, uint image_width_32, uint image_height_32, float3  v_rgb_9, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_10, Matrix<float, 3, 3>  * v_R_12, float3  * v_t_11)
{
    float3  _S6220 = s_primal_ctx_mul_0(R_38, pos_5) + t_38;
    float _S6221 = length_1(_S6220);
    float _S6222 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S6221);
    float _S6223 = s_primal_ctx_max_0(0.0f, _S6221);
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S6224 = s_primal_ctx_mul_0(R_38, pos_i_7) + t_38;
    float _S6225 = length_1(_S6224);
    float _S6226 = s_primal_ctx_min_0(_S6222, _S6225);
    float _S6227 = s_primal_ctx_max_0(_S6223, _S6225);
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S6228 = s_primal_ctx_mul_0(R_38, pos_i_8) + t_38;
    float _S6229 = length_1(_S6228);
    float _S6230 = s_primal_ctx_min_0(_S6226, _S6229);
    float _S6231 = s_primal_ctx_max_0(_S6227, _S6229);
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S6232 = s_primal_ctx_mul_0(R_38, pos_i_9) + t_38;
    float _S6233 = length_1(_S6232);
    float _S6234 = s_primal_ctx_min_0(_S6230, _S6233);
    float _S6235 = s_primal_ctx_max_0(_S6231, _S6233);
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S6236 = s_primal_ctx_mul_0(R_38, pos_i_10) + t_38;
    float _S6237 = length_1(_S6236);
    float _S6238 = s_primal_ctx_min_0(_S6234, _S6237);
    float _S6239 = s_primal_ctx_max_0(_S6235, _S6237);
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S6240 = s_primal_ctx_mul_0(R_38, pos_i_11) + t_38;
    float _S6241 = length_1(_S6240);
    float _S6242 = s_primal_ctx_min_0(_S6238, _S6241);
    float _S6243 = s_primal_ctx_max_0(_S6239, _S6241);
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S6244 = s_primal_ctx_mul_0(R_38, pos_i_12) + t_38;
    float _S6245 = length_1(_S6244);
    float _S6246 = s_primal_ctx_min_0(_S6242, _S6245);
    float _S6247 = s_primal_ctx_max_0(_S6243, _S6245);
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S6248 = s_primal_ctx_mul_0(R_38, pos_i_13) + t_38;
    float _S6249 = length_1(_S6248);
    float3  _S6250 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_32 = s_primal_ctx_mul_0(R_38, _S6250) + t_38;
    float2  _S6251 = float2 {_S6220.x, _S6220.y};
    float _S6252 = length_0(_S6251);
    float _S6253 = _S6220.z;
    float _S6254 = s_primal_ctx_atan2_0(_S6252, _S6253);
    bool _S6255 = _S6254 < 0.00100000004749745f;
    float k_15;
    float _S6256;
    float _S6257;
    float _S6258;
    if(_S6255)
    {
        float _S6259 = 1.0f - _S6254 * _S6254 / 3.0f;
        float _S6260 = _S6253 * _S6253;
        k_15 = _S6259 / _S6253;
        _S6256 = 0.0f;
        _S6257 = _S6260;
        _S6258 = _S6259;
    }
    else
    {
        float _S6261 = _S6252 * _S6252;
        k_15 = _S6254 / _S6252;
        _S6256 = _S6261;
        _S6257 = 0.0f;
        _S6258 = 0.0f;
    }
    float2  _S6262 = make_float2 (k_15);
    float2  _S6263 = _S6251 * make_float2 (k_15);
    float u_125 = _S6263.x;
    float v_125 = _S6263.y;
    float r2_120 = u_125 * u_125 + v_125 * v_125;
    float _S6264 = (*dist_coeffs_45)[int(2)] + r2_120 * (*dist_coeffs_45)[int(3)];
    float _S6265 = (*dist_coeffs_45)[int(1)] + r2_120 * _S6264;
    float _S6266 = (*dist_coeffs_45)[int(0)] + r2_120 * _S6265;
    float radial_21 = 1.0f + r2_120 * _S6266;
    float2  _S6267 = make_float2 (radial_21);
    float _S6268 = 2.0f * (*dist_coeffs_45)[int(4)];
    float _S6269 = _S6268 * u_125;
    float _S6270 = 2.0f * u_125;
    float _S6271 = 2.0f * (*dist_coeffs_45)[int(5)];
    float _S6272 = _S6271 * u_125;
    float _S6273 = 2.0f * v_125;
    float2  _S6274 = _S6263 * make_float2 (radial_21) + make_float2 (_S6269 * v_125 + (*dist_coeffs_45)[int(5)] * (r2_120 + _S6270 * u_125) + (*dist_coeffs_45)[int(6)] * r2_120, _S6272 * v_125 + (*dist_coeffs_45)[int(4)] * (r2_120 + _S6273 * v_125) + (*dist_coeffs_45)[int(7)] * r2_120);
    float2  _S6275 = _S6274 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6274.x + (*dist_coeffs_45)[int(9)] * _S6274.y, 0.0f);
    float _S6276 = fx_41 * _S6275.x + cx_36;
    float _S6277 = fy_41 * _S6275.y + cy_36;
    float2  _S6278 = float2 {_S6224.x, _S6224.y};
    float _S6279 = length_0(_S6278);
    float _S6280 = _S6224.z;
    float _S6281 = s_primal_ctx_atan2_0(_S6279, _S6280);
    bool _S6282 = _S6281 < 0.00100000004749745f;
    float _S6283;
    float _S6284;
    float _S6285;
    if(_S6282)
    {
        float _S6286 = 1.0f - _S6281 * _S6281 / 3.0f;
        float _S6287 = _S6280 * _S6280;
        k_15 = _S6286 / _S6280;
        _S6283 = 0.0f;
        _S6284 = _S6287;
        _S6285 = _S6286;
    }
    else
    {
        float _S6288 = _S6279 * _S6279;
        k_15 = _S6281 / _S6279;
        _S6283 = _S6288;
        _S6284 = 0.0f;
        _S6285 = 0.0f;
    }
    float2  _S6289 = make_float2 (k_15);
    float2  _S6290 = _S6278 * make_float2 (k_15);
    float u_126 = _S6290.x;
    float v_126 = _S6290.y;
    float r2_121 = u_126 * u_126 + v_126 * v_126;
    float _S6291 = (*dist_coeffs_45)[int(2)] + r2_121 * (*dist_coeffs_45)[int(3)];
    float _S6292 = (*dist_coeffs_45)[int(1)] + r2_121 * _S6291;
    float _S6293 = (*dist_coeffs_45)[int(0)] + r2_121 * _S6292;
    float radial_22 = 1.0f + r2_121 * _S6293;
    float2  _S6294 = make_float2 (radial_22);
    float _S6295 = _S6268 * u_126;
    float _S6296 = 2.0f * u_126;
    float _S6297 = _S6271 * u_126;
    float _S6298 = 2.0f * v_126;
    float2  _S6299 = _S6290 * make_float2 (radial_22) + make_float2 (_S6295 * v_126 + (*dist_coeffs_45)[int(5)] * (r2_121 + _S6296 * u_126) + (*dist_coeffs_45)[int(6)] * r2_121, _S6297 * v_126 + (*dist_coeffs_45)[int(4)] * (r2_121 + _S6298 * v_126) + (*dist_coeffs_45)[int(7)] * r2_121);
    float2  _S6300 = _S6299 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6299.x + (*dist_coeffs_45)[int(9)] * _S6299.y, 0.0f);
    float _S6301 = fx_41 * _S6300.x + cx_36;
    float _S6302 = fy_41 * _S6300.y + cy_36;
    float2  _S6303 = float2 {_S6228.x, _S6228.y};
    float _S6304 = length_0(_S6303);
    float _S6305 = _S6228.z;
    float _S6306 = s_primal_ctx_atan2_0(_S6304, _S6305);
    bool _S6307 = _S6306 < 0.00100000004749745f;
    float _S6308;
    float _S6309;
    float _S6310;
    if(_S6307)
    {
        float _S6311 = 1.0f - _S6306 * _S6306 / 3.0f;
        float _S6312 = _S6305 * _S6305;
        k_15 = _S6311 / _S6305;
        _S6308 = 0.0f;
        _S6309 = _S6312;
        _S6310 = _S6311;
    }
    else
    {
        float _S6313 = _S6304 * _S6304;
        k_15 = _S6306 / _S6304;
        _S6308 = _S6313;
        _S6309 = 0.0f;
        _S6310 = 0.0f;
    }
    float2  _S6314 = make_float2 (k_15);
    float2  _S6315 = _S6303 * make_float2 (k_15);
    float u_127 = _S6315.x;
    float v_127 = _S6315.y;
    float r2_122 = u_127 * u_127 + v_127 * v_127;
    float _S6316 = (*dist_coeffs_45)[int(2)] + r2_122 * (*dist_coeffs_45)[int(3)];
    float _S6317 = (*dist_coeffs_45)[int(1)] + r2_122 * _S6316;
    float _S6318 = (*dist_coeffs_45)[int(0)] + r2_122 * _S6317;
    float radial_23 = 1.0f + r2_122 * _S6318;
    float2  _S6319 = make_float2 (radial_23);
    float _S6320 = _S6268 * u_127;
    float _S6321 = 2.0f * u_127;
    float _S6322 = _S6271 * u_127;
    float _S6323 = 2.0f * v_127;
    float2  _S6324 = _S6315 * make_float2 (radial_23) + make_float2 (_S6320 * v_127 + (*dist_coeffs_45)[int(5)] * (r2_122 + _S6321 * u_127) + (*dist_coeffs_45)[int(6)] * r2_122, _S6322 * v_127 + (*dist_coeffs_45)[int(4)] * (r2_122 + _S6323 * v_127) + (*dist_coeffs_45)[int(7)] * r2_122);
    float2  _S6325 = _S6324 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6324.x + (*dist_coeffs_45)[int(9)] * _S6324.y, 0.0f);
    float _S6326 = fx_41 * _S6325.x + cx_36;
    float _S6327 = fy_41 * _S6325.y + cy_36;
    float2  _S6328 = float2 {_S6232.x, _S6232.y};
    float _S6329 = length_0(_S6328);
    float _S6330 = _S6232.z;
    float _S6331 = s_primal_ctx_atan2_0(_S6329, _S6330);
    bool _S6332 = _S6331 < 0.00100000004749745f;
    float _S6333;
    float _S6334;
    float _S6335;
    if(_S6332)
    {
        float _S6336 = 1.0f - _S6331 * _S6331 / 3.0f;
        float _S6337 = _S6330 * _S6330;
        k_15 = _S6336 / _S6330;
        _S6333 = 0.0f;
        _S6334 = _S6337;
        _S6335 = _S6336;
    }
    else
    {
        float _S6338 = _S6329 * _S6329;
        k_15 = _S6331 / _S6329;
        _S6333 = _S6338;
        _S6334 = 0.0f;
        _S6335 = 0.0f;
    }
    float2  _S6339 = make_float2 (k_15);
    float2  _S6340 = _S6328 * make_float2 (k_15);
    float u_128 = _S6340.x;
    float v_128 = _S6340.y;
    float r2_123 = u_128 * u_128 + v_128 * v_128;
    float _S6341 = (*dist_coeffs_45)[int(2)] + r2_123 * (*dist_coeffs_45)[int(3)];
    float _S6342 = (*dist_coeffs_45)[int(1)] + r2_123 * _S6341;
    float _S6343 = (*dist_coeffs_45)[int(0)] + r2_123 * _S6342;
    float radial_24 = 1.0f + r2_123 * _S6343;
    float2  _S6344 = make_float2 (radial_24);
    float _S6345 = _S6268 * u_128;
    float _S6346 = 2.0f * u_128;
    float _S6347 = _S6271 * u_128;
    float _S6348 = 2.0f * v_128;
    float2  _S6349 = _S6340 * make_float2 (radial_24) + make_float2 (_S6345 * v_128 + (*dist_coeffs_45)[int(5)] * (r2_123 + _S6346 * u_128) + (*dist_coeffs_45)[int(6)] * r2_123, _S6347 * v_128 + (*dist_coeffs_45)[int(4)] * (r2_123 + _S6348 * v_128) + (*dist_coeffs_45)[int(7)] * r2_123);
    float2  _S6350 = _S6349 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6349.x + (*dist_coeffs_45)[int(9)] * _S6349.y, 0.0f);
    float _S6351 = fx_41 * _S6350.x + cx_36;
    float _S6352 = fy_41 * _S6350.y + cy_36;
    float2  _S6353 = float2 {_S6236.x, _S6236.y};
    float _S6354 = length_0(_S6353);
    float _S6355 = _S6236.z;
    float _S6356 = s_primal_ctx_atan2_0(_S6354, _S6355);
    bool _S6357 = _S6356 < 0.00100000004749745f;
    float _S6358;
    float _S6359;
    float _S6360;
    if(_S6357)
    {
        float _S6361 = 1.0f - _S6356 * _S6356 / 3.0f;
        float _S6362 = _S6355 * _S6355;
        k_15 = _S6361 / _S6355;
        _S6358 = 0.0f;
        _S6359 = _S6362;
        _S6360 = _S6361;
    }
    else
    {
        float _S6363 = _S6354 * _S6354;
        k_15 = _S6356 / _S6354;
        _S6358 = _S6363;
        _S6359 = 0.0f;
        _S6360 = 0.0f;
    }
    float2  _S6364 = make_float2 (k_15);
    float2  _S6365 = _S6353 * make_float2 (k_15);
    float u_129 = _S6365.x;
    float v_129 = _S6365.y;
    float r2_124 = u_129 * u_129 + v_129 * v_129;
    float _S6366 = (*dist_coeffs_45)[int(2)] + r2_124 * (*dist_coeffs_45)[int(3)];
    float _S6367 = (*dist_coeffs_45)[int(1)] + r2_124 * _S6366;
    float _S6368 = (*dist_coeffs_45)[int(0)] + r2_124 * _S6367;
    float radial_25 = 1.0f + r2_124 * _S6368;
    float2  _S6369 = make_float2 (radial_25);
    float _S6370 = _S6268 * u_129;
    float _S6371 = 2.0f * u_129;
    float _S6372 = _S6271 * u_129;
    float _S6373 = 2.0f * v_129;
    float2  _S6374 = _S6365 * make_float2 (radial_25) + make_float2 (_S6370 * v_129 + (*dist_coeffs_45)[int(5)] * (r2_124 + _S6371 * u_129) + (*dist_coeffs_45)[int(6)] * r2_124, _S6372 * v_129 + (*dist_coeffs_45)[int(4)] * (r2_124 + _S6373 * v_129) + (*dist_coeffs_45)[int(7)] * r2_124);
    float2  _S6375 = _S6374 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6374.x + (*dist_coeffs_45)[int(9)] * _S6374.y, 0.0f);
    float _S6376 = fx_41 * _S6375.x + cx_36;
    float _S6377 = fy_41 * _S6375.y + cy_36;
    float2  _S6378 = float2 {_S6240.x, _S6240.y};
    float _S6379 = length_0(_S6378);
    float _S6380 = _S6240.z;
    float _S6381 = s_primal_ctx_atan2_0(_S6379, _S6380);
    bool _S6382 = _S6381 < 0.00100000004749745f;
    float _S6383;
    float _S6384;
    float _S6385;
    if(_S6382)
    {
        float _S6386 = 1.0f - _S6381 * _S6381 / 3.0f;
        float _S6387 = _S6380 * _S6380;
        k_15 = _S6386 / _S6380;
        _S6383 = 0.0f;
        _S6384 = _S6387;
        _S6385 = _S6386;
    }
    else
    {
        float _S6388 = _S6379 * _S6379;
        k_15 = _S6381 / _S6379;
        _S6383 = _S6388;
        _S6384 = 0.0f;
        _S6385 = 0.0f;
    }
    float2  _S6389 = make_float2 (k_15);
    float2  _S6390 = _S6378 * make_float2 (k_15);
    float u_130 = _S6390.x;
    float v_130 = _S6390.y;
    float r2_125 = u_130 * u_130 + v_130 * v_130;
    float _S6391 = (*dist_coeffs_45)[int(2)] + r2_125 * (*dist_coeffs_45)[int(3)];
    float _S6392 = (*dist_coeffs_45)[int(1)] + r2_125 * _S6391;
    float _S6393 = (*dist_coeffs_45)[int(0)] + r2_125 * _S6392;
    float radial_26 = 1.0f + r2_125 * _S6393;
    float2  _S6394 = make_float2 (radial_26);
    float _S6395 = _S6268 * u_130;
    float _S6396 = 2.0f * u_130;
    float _S6397 = _S6271 * u_130;
    float _S6398 = 2.0f * v_130;
    float2  _S6399 = _S6390 * make_float2 (radial_26) + make_float2 (_S6395 * v_130 + (*dist_coeffs_45)[int(5)] * (r2_125 + _S6396 * u_130) + (*dist_coeffs_45)[int(6)] * r2_125, _S6397 * v_130 + (*dist_coeffs_45)[int(4)] * (r2_125 + _S6398 * v_130) + (*dist_coeffs_45)[int(7)] * r2_125);
    float2  _S6400 = _S6399 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6399.x + (*dist_coeffs_45)[int(9)] * _S6399.y, 0.0f);
    float _S6401 = fx_41 * _S6400.x + cx_36;
    float _S6402 = fy_41 * _S6400.y + cy_36;
    float2  _S6403 = float2 {_S6244.x, _S6244.y};
    float _S6404 = length_0(_S6403);
    float _S6405 = _S6244.z;
    float _S6406 = s_primal_ctx_atan2_0(_S6404, _S6405);
    bool _S6407 = _S6406 < 0.00100000004749745f;
    float _S6408;
    float _S6409;
    float _S6410;
    if(_S6407)
    {
        float _S6411 = 1.0f - _S6406 * _S6406 / 3.0f;
        float _S6412 = _S6405 * _S6405;
        k_15 = _S6411 / _S6405;
        _S6408 = 0.0f;
        _S6409 = _S6412;
        _S6410 = _S6411;
    }
    else
    {
        float _S6413 = _S6404 * _S6404;
        k_15 = _S6406 / _S6404;
        _S6408 = _S6413;
        _S6409 = 0.0f;
        _S6410 = 0.0f;
    }
    float2  _S6414 = make_float2 (k_15);
    float2  _S6415 = _S6403 * make_float2 (k_15);
    float u_131 = _S6415.x;
    float v_131 = _S6415.y;
    float r2_126 = u_131 * u_131 + v_131 * v_131;
    float _S6416 = (*dist_coeffs_45)[int(2)] + r2_126 * (*dist_coeffs_45)[int(3)];
    float _S6417 = (*dist_coeffs_45)[int(1)] + r2_126 * _S6416;
    float _S6418 = (*dist_coeffs_45)[int(0)] + r2_126 * _S6417;
    float radial_27 = 1.0f + r2_126 * _S6418;
    float2  _S6419 = make_float2 (radial_27);
    float _S6420 = _S6268 * u_131;
    float _S6421 = 2.0f * u_131;
    float _S6422 = _S6271 * u_131;
    float _S6423 = 2.0f * v_131;
    float2  _S6424 = _S6415 * make_float2 (radial_27) + make_float2 (_S6420 * v_131 + (*dist_coeffs_45)[int(5)] * (r2_126 + _S6421 * u_131) + (*dist_coeffs_45)[int(6)] * r2_126, _S6422 * v_131 + (*dist_coeffs_45)[int(4)] * (r2_126 + _S6423 * v_131) + (*dist_coeffs_45)[int(7)] * r2_126);
    float2  _S6425 = _S6424 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6424.x + (*dist_coeffs_45)[int(9)] * _S6424.y, 0.0f);
    float _S6426 = fx_41 * _S6425.x + cx_36;
    float _S6427 = fy_41 * _S6425.y + cy_36;
    float2  _S6428 = float2 {_S6248.x, _S6248.y};
    float _S6429 = length_0(_S6428);
    float _S6430 = _S6248.z;
    float _S6431 = s_primal_ctx_atan2_0(_S6429, _S6430);
    bool _S6432 = _S6431 < 0.00100000004749745f;
    float _S6433;
    float _S6434;
    float _S6435;
    if(_S6432)
    {
        float _S6436 = 1.0f - _S6431 * _S6431 / 3.0f;
        float _S6437 = _S6430 * _S6430;
        k_15 = _S6436 / _S6430;
        _S6433 = 0.0f;
        _S6434 = _S6437;
        _S6435 = _S6436;
    }
    else
    {
        float _S6438 = _S6429 * _S6429;
        k_15 = _S6431 / _S6429;
        _S6433 = _S6438;
        _S6434 = 0.0f;
        _S6435 = 0.0f;
    }
    float2  _S6439 = make_float2 (k_15);
    float2  _S6440 = _S6428 * make_float2 (k_15);
    float u_132 = _S6440.x;
    float v_132 = _S6440.y;
    float r2_127 = u_132 * u_132 + v_132 * v_132;
    float _S6441 = (*dist_coeffs_45)[int(2)] + r2_127 * (*dist_coeffs_45)[int(3)];
    float _S6442 = (*dist_coeffs_45)[int(1)] + r2_127 * _S6441;
    float _S6443 = (*dist_coeffs_45)[int(0)] + r2_127 * _S6442;
    float radial_28 = 1.0f + r2_127 * _S6443;
    float2  _S6444 = make_float2 (radial_28);
    float _S6445 = _S6268 * u_132;
    float _S6446 = 2.0f * u_132;
    float _S6447 = _S6271 * u_132;
    float _S6448 = 2.0f * v_132;
    float2  _S6449 = _S6440 * make_float2 (radial_28) + make_float2 (_S6445 * v_132 + (*dist_coeffs_45)[int(5)] * (r2_127 + _S6446 * u_132) + (*dist_coeffs_45)[int(6)] * r2_127, _S6447 * v_132 + (*dist_coeffs_45)[int(4)] * (r2_127 + _S6448 * v_132) + (*dist_coeffs_45)[int(7)] * r2_127);
    float2  _S6450 = _S6449 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6449.x + (*dist_coeffs_45)[int(9)] * _S6449.y, 0.0f);
    float _S6451 = fx_41 * _S6450.x + cx_36;
    float _S6452 = fy_41 * _S6450.y + cy_36;
    float _S6453 = s_primal_ctx_max_0(_S6276, _S6301);
    float _S6454 = s_primal_ctx_min_0(_S6276, _S6301);
    float _S6455 = s_primal_ctx_max_0(_S6277, _S6302);
    float _S6456 = s_primal_ctx_min_0(_S6277, _S6302);
    float _S6457 = s_primal_ctx_max_0(_S6453, _S6326);
    float _S6458 = s_primal_ctx_min_0(_S6454, _S6326);
    float _S6459 = s_primal_ctx_max_0(_S6455, _S6327);
    float _S6460 = s_primal_ctx_min_0(_S6456, _S6327);
    float _S6461 = s_primal_ctx_max_0(_S6457, _S6351);
    float _S6462 = s_primal_ctx_min_0(_S6458, _S6351);
    float _S6463 = s_primal_ctx_max_0(_S6459, _S6352);
    float _S6464 = s_primal_ctx_min_0(_S6460, _S6352);
    float _S6465 = s_primal_ctx_max_0(_S6461, _S6376);
    float _S6466 = s_primal_ctx_min_0(_S6462, _S6376);
    float _S6467 = s_primal_ctx_max_0(_S6463, _S6377);
    float _S6468 = s_primal_ctx_min_0(_S6464, _S6377);
    float _S6469 = s_primal_ctx_max_0(_S6465, _S6401);
    float _S6470 = s_primal_ctx_min_0(_S6466, _S6401);
    float _S6471 = s_primal_ctx_max_0(_S6467, _S6402);
    float _S6472 = s_primal_ctx_min_0(_S6468, _S6402);
    float _S6473 = s_primal_ctx_max_0(_S6469, _S6426);
    float _S6474 = s_primal_ctx_min_0(_S6470, _S6426);
    float _S6475 = s_primal_ctx_max_0(_S6471, _S6427);
    float _S6476 = s_primal_ctx_min_0(_S6472, _S6427);
    float _S6477 = s_primal_ctx_dot_0(mean_c_32, mean_c_32) + 9.99999997475242708e-07f;
    Matrix<float, 3, 3>  _S6478 = transpose_0(R_38);
    float3  _S6479 = mean_c_32 - - s_primal_ctx_mul_0(_S6478, t_38);
    float _S6480 = _S6479.x;
    float _S6481 = _S6479.y;
    float _S6482 = _S6479.z;
    float _S6483 = _S6480 * _S6480 + _S6481 * _S6481 + _S6482 * _S6482;
    float _S6484 = s_primal_ctx_sqrt_0(_S6483);
    float x_68 = _S6480 / _S6484;
    float3  _S6485 = make_float3 (x_68);
    float _S6486 = _S6484 * _S6484;
    float y_35 = _S6481 / _S6484;
    float z_32 = _S6482 / _S6484;
    float3  _S6487 = make_float3 (z_32);
    float _S6488 = - y_35;
    float3  _S6489 = make_float3 (_S6488);
    float z2_65 = z_32 * z_32;
    float fTmp0B_32 = -1.09254848957061768f * z_32;
    float fC1_32 = x_68 * x_68 - y_35 * y_35;
    float _S6490 = 2.0f * x_68;
    float fS1_32 = _S6490 * y_35;
    float pSH6_10 = 0.94617468118667603f * z2_65 - 0.31539157032966614f;
    float3  _S6491 = make_float3 (pSH6_10);
    float pSH7_10 = fTmp0B_32 * x_68;
    float3  _S6492 = make_float3 (pSH7_10);
    float pSH5_10 = fTmp0B_32 * y_35;
    float3  _S6493 = make_float3 (pSH5_10);
    float pSH8_10 = 0.54627424478530884f * fC1_32;
    float3  _S6494 = make_float3 (pSH8_10);
    float pSH4_10 = 0.54627424478530884f * fS1_32;
    float3  _S6495 = make_float3 (pSH4_10);
    float fTmp0C_32 = -2.28522896766662598f * z2_65 + 0.4570457935333252f;
    float fTmp1B_32 = 1.44530570507049561f * z_32;
    float _S6496 = 1.86588168144226074f * z2_65 - 1.11952900886535645f;
    float pSH12_10 = z_32 * _S6496;
    float3  _S6497 = make_float3 (pSH12_10);
    float pSH13_10 = fTmp0C_32 * x_68;
    float3  _S6498 = make_float3 (pSH13_10);
    float pSH11_10 = fTmp0C_32 * y_35;
    float3  _S6499 = make_float3 (pSH11_10);
    float pSH14_10 = fTmp1B_32 * fC1_32;
    float3  _S6500 = make_float3 (pSH14_10);
    float pSH10_10 = fTmp1B_32 * fS1_32;
    float3  _S6501 = make_float3 (pSH10_10);
    float pSH15_10 = -0.59004360437393188f * (x_68 * fC1_32 - y_35 * fS1_32);
    float3  _S6502 = make_float3 (pSH15_10);
    float pSH9_10 = -0.59004360437393188f * (x_68 * fS1_32 + y_35 * fC1_32);
    float3  _S6503 = make_float3 (pSH9_10);
    float3  _S6504 = make_float3 (0.0f);
    float3  _S6505 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6506;
    (&_S6506)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_32)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S6488) * (*sh_coeffs_32)[int(1)] + make_float3 (z_32) * (*sh_coeffs_32)[int(2)] - make_float3 (x_68) * (*sh_coeffs_32)[int(3)]) + (make_float3 (pSH4_10) * (*sh_coeffs_32)[int(4)] + make_float3 (pSH5_10) * (*sh_coeffs_32)[int(5)] + make_float3 (pSH6_10) * (*sh_coeffs_32)[int(6)] + make_float3 (pSH7_10) * (*sh_coeffs_32)[int(7)] + make_float3 (pSH8_10) * (*sh_coeffs_32)[int(8)]) + (make_float3 (pSH9_10) * (*sh_coeffs_32)[int(9)] + make_float3 (pSH10_10) * (*sh_coeffs_32)[int(10)] + make_float3 (pSH11_10) * (*sh_coeffs_32)[int(11)] + make_float3 (pSH12_10) * (*sh_coeffs_32)[int(12)] + make_float3 (pSH13_10) * (*sh_coeffs_32)[int(13)] + make_float3 (pSH14_10) * (*sh_coeffs_32)[int(14)] + make_float3 (pSH15_10) * (*sh_coeffs_32)[int(15)]) + make_float3 (0.5f);
    (&_S6506)->differential_0 = _S6505;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6507;
    (&_S6507)->primal_0 = _S6504;
    (&_S6507)->differential_0 = _S6505;
    s_bwd_prop_max_0(&_S6506, &_S6507, v_rgb_9);
    float3  _S6508 = _S6502 * _S6506.differential_0;
    float3  _S6509 = (*sh_coeffs_32)[int(15)] * _S6506.differential_0;
    float3  _S6510 = _S6500 * _S6506.differential_0;
    float3  _S6511 = (*sh_coeffs_32)[int(14)] * _S6506.differential_0;
    float3  _S6512 = _S6498 * _S6506.differential_0;
    float3  _S6513 = (*sh_coeffs_32)[int(13)] * _S6506.differential_0;
    float3  _S6514 = _S6497 * _S6506.differential_0;
    float3  _S6515 = (*sh_coeffs_32)[int(12)] * _S6506.differential_0;
    float3  _S6516 = _S6499 * _S6506.differential_0;
    float3  _S6517 = (*sh_coeffs_32)[int(11)] * _S6506.differential_0;
    float3  _S6518 = _S6501 * _S6506.differential_0;
    float3  _S6519 = (*sh_coeffs_32)[int(10)] * _S6506.differential_0;
    float3  _S6520 = _S6503 * _S6506.differential_0;
    float3  _S6521 = (*sh_coeffs_32)[int(9)] * _S6506.differential_0;
    float s_diff_fS2_T_10 = -0.59004360437393188f * (_S6521.x + _S6521.y + _S6521.z);
    float s_diff_fC2_T_10 = -0.59004360437393188f * (_S6509.x + _S6509.y + _S6509.z);
    float _S6522 = _S6519.x + _S6519.y + _S6519.z;
    float _S6523 = _S6511.x + _S6511.y + _S6511.z;
    float _S6524 = _S6517.x + _S6517.y + _S6517.z;
    float _S6525 = _S6513.x + _S6513.y + _S6513.z;
    float _S6526 = _S6515.x + _S6515.y + _S6515.z;
    float _S6527 = - s_diff_fC2_T_10;
    float3  _S6528 = _S6494 * _S6506.differential_0;
    float3  _S6529 = (*sh_coeffs_32)[int(8)] * _S6506.differential_0;
    float3  _S6530 = _S6492 * _S6506.differential_0;
    float3  _S6531 = (*sh_coeffs_32)[int(7)] * _S6506.differential_0;
    float3  _S6532 = _S6491 * _S6506.differential_0;
    float3  _S6533 = (*sh_coeffs_32)[int(6)] * _S6506.differential_0;
    float3  _S6534 = _S6493 * _S6506.differential_0;
    float3  _S6535 = (*sh_coeffs_32)[int(5)] * _S6506.differential_0;
    float3  _S6536 = _S6495 * _S6506.differential_0;
    float3  _S6537 = (*sh_coeffs_32)[int(4)] * _S6506.differential_0;
    float _S6538 = _S6535.x + _S6535.y + _S6535.z;
    float _S6539 = _S6531.x + _S6531.y + _S6531.z;
    float _S6540 = fTmp1B_32 * _S6522 + x_68 * s_diff_fS2_T_10 + y_35 * _S6527 + 0.54627424478530884f * (_S6537.x + _S6537.y + _S6537.z);
    float _S6541 = fTmp1B_32 * _S6523 + y_35 * s_diff_fS2_T_10 + x_68 * s_diff_fC2_T_10 + 0.54627424478530884f * (_S6529.x + _S6529.y + _S6529.z);
    float _S6542 = y_35 * - _S6541;
    float _S6543 = x_68 * _S6541;
    float _S6544 = z_32 * (1.86588168144226074f * (z_32 * _S6526) + -2.28522896766662598f * (y_35 * _S6524 + x_68 * _S6525) + 0.94617468118667603f * (_S6533.x + _S6533.y + _S6533.z));
    float3  _S6545 = make_float3 (0.48860251903533936f) * _S6506.differential_0;
    float3  _S6546 = - _S6545;
    float3  _S6547 = _S6485 * _S6546;
    float3  _S6548 = (*sh_coeffs_32)[int(3)] * _S6546;
    float3  _S6549 = _S6487 * _S6545;
    float3  _S6550 = (*sh_coeffs_32)[int(2)] * _S6545;
    float3  _S6551 = _S6489 * _S6545;
    float3  _S6552 = (*sh_coeffs_32)[int(1)] * _S6545;
    float _S6553 = (_S6496 * _S6526 + 1.44530570507049561f * (fS1_32 * _S6522 + fC1_32 * _S6523) + -1.09254848957061768f * (y_35 * _S6538 + x_68 * _S6539) + _S6544 + _S6544 + _S6550.x + _S6550.y + _S6550.z) / _S6486;
    float _S6554 = _S6484 * _S6553;
    float _S6555 = (fTmp0C_32 * _S6524 + fC1_32 * s_diff_fS2_T_10 + fS1_32 * _S6527 + fTmp0B_32 * _S6538 + _S6490 * _S6540 + _S6542 + _S6542 + - (_S6552.x + _S6552.y + _S6552.z)) / _S6486;
    float _S6556 = _S6484 * _S6555;
    float _S6557 = (fTmp0C_32 * _S6525 + fS1_32 * s_diff_fS2_T_10 + fC1_32 * s_diff_fC2_T_10 + fTmp0B_32 * _S6539 + 2.0f * (y_35 * _S6540) + _S6543 + _S6543 + _S6548.x + _S6548.y + _S6548.z) / _S6486;
    float _S6558 = _S6484 * _S6557;
    float _S6559 = _S6482 * - _S6553 + _S6481 * - _S6555 + _S6480 * - _S6557;
    DiffPair_float_0 _S6560;
    (&_S6560)->primal_0 = _S6483;
    (&_S6560)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S6560, _S6559);
    float _S6561 = _S6482 * _S6560.differential_0;
    float _S6562 = _S6481 * _S6560.differential_0;
    float _S6563 = _S6480 * _S6560.differential_0;
    float3  _S6564 = make_float3 (0.282094806432724f) * _S6506.differential_0;
    float3  _S6565 = make_float3 (_S6558 + _S6563 + _S6563, _S6556 + _S6562 + _S6562, _S6554 + _S6561 + _S6561);
    float3  _S6566 = - - _S6565;
    Matrix<float, 3, 3>  _S6567 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6568;
    (&_S6568)->primal_0 = _S6478;
    (&_S6568)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6569;
    (&_S6569)->primal_0 = t_38;
    (&_S6569)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6568, &_S6569, _S6566);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6570 = _S6569;
    Matrix<float, 3, 3>  _S6571 = transpose_0(_S6568.differential_0);
    DiffPair_float_0 _S6572;
    (&_S6572)->primal_0 = _S6477;
    (&_S6572)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S6572, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6573;
    (&_S6573)->primal_0 = mean_c_32;
    (&_S6573)->differential_0 = _S6505;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6574;
    (&_S6574)->primal_0 = mean_c_32;
    (&_S6574)->differential_0 = _S6505;
    s_bwd_prop_dot_0(&_S6573, &_S6574, _S6572.differential_0);
    DiffPair_float_0 _S6575;
    (&_S6575)->primal_0 = _S6476;
    (&_S6575)->differential_0 = 0.0f;
    DiffPair_float_0 _S6576;
    (&_S6576)->primal_0 = _S6452;
    (&_S6576)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6575, &_S6576, 0.0f);
    DiffPair_float_0 _S6577;
    (&_S6577)->primal_0 = _S6475;
    (&_S6577)->differential_0 = 0.0f;
    DiffPair_float_0 _S6578;
    (&_S6578)->primal_0 = _S6452;
    (&_S6578)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6577, &_S6578, 0.0f);
    DiffPair_float_0 _S6579;
    (&_S6579)->primal_0 = _S6474;
    (&_S6579)->differential_0 = 0.0f;
    DiffPair_float_0 _S6580;
    (&_S6580)->primal_0 = _S6451;
    (&_S6580)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6579, &_S6580, 0.0f);
    DiffPair_float_0 _S6581;
    (&_S6581)->primal_0 = _S6473;
    (&_S6581)->differential_0 = 0.0f;
    DiffPair_float_0 _S6582;
    (&_S6582)->primal_0 = _S6451;
    (&_S6582)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6581, &_S6582, 0.0f);
    DiffPair_float_0 _S6583;
    (&_S6583)->primal_0 = _S6472;
    (&_S6583)->differential_0 = 0.0f;
    DiffPair_float_0 _S6584;
    (&_S6584)->primal_0 = _S6427;
    (&_S6584)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6583, &_S6584, _S6575.differential_0);
    DiffPair_float_0 _S6585;
    (&_S6585)->primal_0 = _S6471;
    (&_S6585)->differential_0 = 0.0f;
    DiffPair_float_0 _S6586;
    (&_S6586)->primal_0 = _S6427;
    (&_S6586)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6585, &_S6586, _S6577.differential_0);
    DiffPair_float_0 _S6587;
    (&_S6587)->primal_0 = _S6470;
    (&_S6587)->differential_0 = 0.0f;
    DiffPair_float_0 _S6588;
    (&_S6588)->primal_0 = _S6426;
    (&_S6588)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6587, &_S6588, _S6579.differential_0);
    DiffPair_float_0 _S6589;
    (&_S6589)->primal_0 = _S6469;
    (&_S6589)->differential_0 = 0.0f;
    DiffPair_float_0 _S6590;
    (&_S6590)->primal_0 = _S6426;
    (&_S6590)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6589, &_S6590, _S6581.differential_0);
    DiffPair_float_0 _S6591;
    (&_S6591)->primal_0 = _S6468;
    (&_S6591)->differential_0 = 0.0f;
    DiffPair_float_0 _S6592;
    (&_S6592)->primal_0 = _S6402;
    (&_S6592)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6591, &_S6592, _S6583.differential_0);
    DiffPair_float_0 _S6593;
    (&_S6593)->primal_0 = _S6467;
    (&_S6593)->differential_0 = 0.0f;
    DiffPair_float_0 _S6594;
    (&_S6594)->primal_0 = _S6402;
    (&_S6594)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6593, &_S6594, _S6585.differential_0);
    DiffPair_float_0 _S6595;
    (&_S6595)->primal_0 = _S6466;
    (&_S6595)->differential_0 = 0.0f;
    DiffPair_float_0 _S6596;
    (&_S6596)->primal_0 = _S6401;
    (&_S6596)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6595, &_S6596, _S6587.differential_0);
    DiffPair_float_0 _S6597;
    (&_S6597)->primal_0 = _S6465;
    (&_S6597)->differential_0 = 0.0f;
    DiffPair_float_0 _S6598;
    (&_S6598)->primal_0 = _S6401;
    (&_S6598)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6597, &_S6598, _S6589.differential_0);
    DiffPair_float_0 _S6599;
    (&_S6599)->primal_0 = _S6464;
    (&_S6599)->differential_0 = 0.0f;
    DiffPair_float_0 _S6600;
    (&_S6600)->primal_0 = _S6377;
    (&_S6600)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6599, &_S6600, _S6591.differential_0);
    DiffPair_float_0 _S6601;
    (&_S6601)->primal_0 = _S6463;
    (&_S6601)->differential_0 = 0.0f;
    DiffPair_float_0 _S6602;
    (&_S6602)->primal_0 = _S6377;
    (&_S6602)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6601, &_S6602, _S6593.differential_0);
    DiffPair_float_0 _S6603;
    (&_S6603)->primal_0 = _S6462;
    (&_S6603)->differential_0 = 0.0f;
    DiffPair_float_0 _S6604;
    (&_S6604)->primal_0 = _S6376;
    (&_S6604)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6603, &_S6604, _S6595.differential_0);
    DiffPair_float_0 _S6605;
    (&_S6605)->primal_0 = _S6461;
    (&_S6605)->differential_0 = 0.0f;
    DiffPair_float_0 _S6606;
    (&_S6606)->primal_0 = _S6376;
    (&_S6606)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6605, &_S6606, _S6597.differential_0);
    DiffPair_float_0 _S6607;
    (&_S6607)->primal_0 = _S6460;
    (&_S6607)->differential_0 = 0.0f;
    DiffPair_float_0 _S6608;
    (&_S6608)->primal_0 = _S6352;
    (&_S6608)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6607, &_S6608, _S6599.differential_0);
    DiffPair_float_0 _S6609;
    (&_S6609)->primal_0 = _S6459;
    (&_S6609)->differential_0 = 0.0f;
    DiffPair_float_0 _S6610;
    (&_S6610)->primal_0 = _S6352;
    (&_S6610)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6609, &_S6610, _S6601.differential_0);
    DiffPair_float_0 _S6611;
    (&_S6611)->primal_0 = _S6458;
    (&_S6611)->differential_0 = 0.0f;
    DiffPair_float_0 _S6612;
    (&_S6612)->primal_0 = _S6351;
    (&_S6612)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6611, &_S6612, _S6603.differential_0);
    DiffPair_float_0 _S6613;
    (&_S6613)->primal_0 = _S6457;
    (&_S6613)->differential_0 = 0.0f;
    DiffPair_float_0 _S6614;
    (&_S6614)->primal_0 = _S6351;
    (&_S6614)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6613, &_S6614, _S6605.differential_0);
    DiffPair_float_0 _S6615;
    (&_S6615)->primal_0 = _S6456;
    (&_S6615)->differential_0 = 0.0f;
    DiffPair_float_0 _S6616;
    (&_S6616)->primal_0 = _S6327;
    (&_S6616)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6615, &_S6616, _S6607.differential_0);
    DiffPair_float_0 _S6617;
    (&_S6617)->primal_0 = _S6455;
    (&_S6617)->differential_0 = 0.0f;
    DiffPair_float_0 _S6618;
    (&_S6618)->primal_0 = _S6327;
    (&_S6618)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6617, &_S6618, _S6609.differential_0);
    DiffPair_float_0 _S6619;
    (&_S6619)->primal_0 = _S6454;
    (&_S6619)->differential_0 = 0.0f;
    DiffPair_float_0 _S6620;
    (&_S6620)->primal_0 = _S6326;
    (&_S6620)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6619, &_S6620, _S6611.differential_0);
    DiffPair_float_0 _S6621;
    (&_S6621)->primal_0 = _S6453;
    (&_S6621)->differential_0 = 0.0f;
    DiffPair_float_0 _S6622;
    (&_S6622)->primal_0 = _S6326;
    (&_S6622)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6621, &_S6622, _S6613.differential_0);
    DiffPair_float_0 _S6623;
    (&_S6623)->primal_0 = _S6277;
    (&_S6623)->differential_0 = 0.0f;
    DiffPair_float_0 _S6624;
    (&_S6624)->primal_0 = _S6302;
    (&_S6624)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6623, &_S6624, _S6615.differential_0);
    DiffPair_float_0 _S6625;
    (&_S6625)->primal_0 = _S6277;
    (&_S6625)->differential_0 = 0.0f;
    DiffPair_float_0 _S6626;
    (&_S6626)->primal_0 = _S6302;
    (&_S6626)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6625, &_S6626, _S6617.differential_0);
    DiffPair_float_0 _S6627;
    (&_S6627)->primal_0 = _S6276;
    (&_S6627)->differential_0 = 0.0f;
    DiffPair_float_0 _S6628;
    (&_S6628)->primal_0 = _S6301;
    (&_S6628)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6627, &_S6628, _S6619.differential_0);
    DiffPair_float_0 _S6629;
    (&_S6629)->primal_0 = _S6276;
    (&_S6629)->differential_0 = 0.0f;
    DiffPair_float_0 _S6630;
    (&_S6630)->primal_0 = _S6301;
    (&_S6630)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6629, &_S6630, _S6621.differential_0);
    float _S6631 = fx_41 * (_S6580.differential_0 + _S6582.differential_0);
    float2  _S6632 = make_float2 (_S6631, fy_41 * (_S6576.differential_0 + _S6578.differential_0)) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6631, (*dist_coeffs_45)[int(9)] * _S6631);
    float2  _S6633 = _S6440 * _S6632;
    float2  _S6634 = _S6444 * _S6632;
    float _S6635 = (*dist_coeffs_45)[int(4)] * _S6632.y;
    float _S6636 = (*dist_coeffs_45)[int(5)] * _S6632.x;
    float _S6637 = _S6633.x + _S6633.y;
    float _S6638 = r2_127 * _S6637;
    float _S6639 = r2_127 * _S6638;
    float _S6640 = (*dist_coeffs_45)[int(7)] * _S6632.y + _S6635 + (*dist_coeffs_45)[int(6)] * _S6632.x + _S6636 + _S6443 * _S6637 + _S6442 * _S6638 + _S6441 * _S6639 + (*dist_coeffs_45)[int(3)] * (r2_127 * _S6639);
    float _S6641 = v_132 * _S6640;
    float _S6642 = u_132 * _S6640;
    float _S6643 = _S6448 * _S6635 + 2.0f * (v_132 * _S6635) + _S6447 * _S6632.y + _S6445 * _S6632.x + _S6641 + _S6641;
    float _S6644 = _S6271 * (v_132 * _S6632.y) + _S6446 * _S6636 + 2.0f * (u_132 * _S6636) + _S6268 * (v_132 * _S6632.x) + _S6642 + _S6642;
    FixedArray<float3 , 16>  _S6645;
    _S6645[int(0)] = _S6505;
    _S6645[int(1)] = _S6505;
    _S6645[int(2)] = _S6505;
    _S6645[int(3)] = _S6505;
    _S6645[int(4)] = _S6505;
    _S6645[int(5)] = _S6505;
    _S6645[int(6)] = _S6505;
    _S6645[int(7)] = _S6505;
    _S6645[int(8)] = _S6505;
    _S6645[int(9)] = _S6505;
    _S6645[int(10)] = _S6505;
    _S6645[int(11)] = _S6505;
    _S6645[int(12)] = _S6505;
    _S6645[int(13)] = _S6505;
    _S6645[int(14)] = _S6505;
    _S6645[int(15)] = _S6505;
    _S6645[int(7)] = _S6530;
    _S6645[int(0)] = _S6564;
    _S6645[int(1)] = _S6551;
    _S6645[int(2)] = _S6549;
    _S6645[int(3)] = _S6547;
    _S6645[int(4)] = _S6536;
    _S6645[int(5)] = _S6534;
    _S6645[int(6)] = _S6532;
    _S6645[int(15)] = _S6508;
    _S6645[int(8)] = _S6528;
    _S6645[int(9)] = _S6520;
    _S6645[int(10)] = _S6518;
    _S6645[int(11)] = _S6516;
    _S6645[int(12)] = _S6514;
    _S6645[int(13)] = _S6512;
    _S6645[int(14)] = _S6510;
    float3  _S6646 = _S6645[int(0)];
    float3  _S6647 = _S6645[int(1)];
    float3  _S6648 = _S6645[int(2)];
    float3  _S6649 = _S6645[int(3)];
    float3  _S6650 = _S6645[int(4)];
    float3  _S6651 = _S6645[int(5)];
    float3  _S6652 = _S6645[int(6)];
    float3  _S6653 = _S6645[int(7)];
    float3  _S6654 = _S6645[int(8)];
    float3  _S6655 = _S6645[int(9)];
    float3  _S6656 = _S6645[int(10)];
    float3  _S6657 = _S6645[int(11)];
    float3  _S6658 = _S6645[int(12)];
    float3  _S6659 = _S6645[int(13)];
    float3  _S6660 = _S6645[int(14)];
    float3  _S6661 = _S6645[int(15)];
    float3  _S6662 = _S6565 + _S6574.differential_0 + _S6573.differential_0;
    float _S6663 = _S6628.differential_0 + _S6630.differential_0;
    float _S6664 = _S6623.differential_0 + _S6625.differential_0;
    float _S6665 = _S6624.differential_0 + _S6626.differential_0;
    float _S6666 = _S6620.differential_0 + _S6622.differential_0;
    float _S6667 = _S6627.differential_0 + _S6629.differential_0;
    float _S6668 = _S6584.differential_0 + _S6586.differential_0;
    float _S6669 = _S6588.differential_0 + _S6590.differential_0;
    float _S6670 = _S6592.differential_0 + _S6594.differential_0;
    float _S6671 = _S6596.differential_0 + _S6598.differential_0;
    float _S6672 = _S6600.differential_0 + _S6602.differential_0;
    float _S6673 = _S6604.differential_0 + _S6606.differential_0;
    float _S6674 = _S6608.differential_0 + _S6610.differential_0;
    float _S6675 = _S6612.differential_0 + _S6614.differential_0;
    float _S6676 = _S6616.differential_0 + _S6618.differential_0;
    float2  _S6677 = _S6634 + make_float2 (_S6644, _S6643);
    float2  _S6678 = _S6428 * _S6677;
    float2  _S6679 = _S6439 * _S6677;
    float _S6680 = _S6678.x + _S6678.y;
    if(_S6432)
    {
        float _S6681 = _S6680 / _S6434;
        float _S6682 = _S6435 * - _S6681;
        float _S6683 = _S6431 * (0.3333333432674408f * - (_S6430 * _S6681));
        k_15 = _S6683 + _S6683;
        _S6433 = _S6682;
        _S6434 = 0.0f;
    }
    else
    {
        float _S6684 = _S6680 / _S6433;
        float _S6685 = _S6431 * - _S6684;
        k_15 = _S6429 * _S6684;
        _S6433 = 0.0f;
        _S6434 = _S6685;
    }
    DiffPair_float_0 _S6686;
    (&_S6686)->primal_0 = _S6429;
    (&_S6686)->differential_0 = 0.0f;
    DiffPair_float_0 _S6687;
    (&_S6687)->primal_0 = _S6430;
    (&_S6687)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6686, &_S6687, k_15);
    float _S6688 = _S6687.differential_0 + _S6433;
    float _S6689 = _S6686.differential_0 + _S6434;
    float2  _S6690 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6691;
    (&_S6691)->primal_0 = _S6428;
    (&_S6691)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6691, _S6689);
    float2  _S6692 = _S6691.differential_0 + _S6679;
    float _S6693 = fx_41 * _S6669;
    float2  _S6694 = make_float2 (_S6693, fy_41 * _S6668) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6693, (*dist_coeffs_45)[int(9)] * _S6693);
    float2  _S6695 = _S6415 * _S6694;
    float _S6696 = (*dist_coeffs_45)[int(4)] * _S6694.y;
    float _S6697 = (*dist_coeffs_45)[int(5)] * _S6694.x;
    float _S6698 = _S6695.x + _S6695.y;
    float _S6699 = r2_126 * _S6698;
    float _S6700 = r2_126 * _S6699;
    float _S6701 = (*dist_coeffs_45)[int(7)] * _S6694.y + _S6696 + (*dist_coeffs_45)[int(6)] * _S6694.x + _S6697 + _S6418 * _S6698 + _S6417 * _S6699 + _S6416 * _S6700 + (*dist_coeffs_45)[int(3)] * (r2_126 * _S6700);
    float _S6702 = v_131 * _S6701;
    float _S6703 = u_131 * _S6701;
    float3  _S6704 = make_float3 (_S6692.x, _S6692.y, _S6688);
    float2  _S6705 = _S6419 * _S6694 + make_float2 (_S6271 * (v_131 * _S6694.y) + _S6421 * _S6697 + 2.0f * (u_131 * _S6697) + _S6268 * (v_131 * _S6694.x) + _S6703 + _S6703, _S6423 * _S6696 + 2.0f * (v_131 * _S6696) + _S6422 * _S6694.y + _S6420 * _S6694.x + _S6702 + _S6702);
    float2  _S6706 = _S6403 * _S6705;
    float2  _S6707 = _S6414 * _S6705;
    float _S6708 = _S6706.x + _S6706.y;
    if(_S6407)
    {
        float _S6709 = _S6708 / _S6409;
        float _S6710 = _S6410 * - _S6709;
        float _S6711 = _S6406 * (0.3333333432674408f * - (_S6405 * _S6709));
        k_15 = _S6711 + _S6711;
        _S6408 = _S6710;
        _S6409 = 0.0f;
    }
    else
    {
        float _S6712 = _S6708 / _S6408;
        float _S6713 = _S6406 * - _S6712;
        k_15 = _S6404 * _S6712;
        _S6408 = 0.0f;
        _S6409 = _S6713;
    }
    DiffPair_float_0 _S6714;
    (&_S6714)->primal_0 = _S6404;
    (&_S6714)->differential_0 = 0.0f;
    DiffPair_float_0 _S6715;
    (&_S6715)->primal_0 = _S6405;
    (&_S6715)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6714, &_S6715, k_15);
    float _S6716 = _S6715.differential_0 + _S6408;
    float _S6717 = _S6714.differential_0 + _S6409;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6718;
    (&_S6718)->primal_0 = _S6403;
    (&_S6718)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6718, _S6717);
    float2  _S6719 = _S6718.differential_0 + _S6707;
    float _S6720 = fx_41 * _S6671;
    float2  _S6721 = make_float2 (_S6720, fy_41 * _S6670) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6720, (*dist_coeffs_45)[int(9)] * _S6720);
    float2  _S6722 = _S6390 * _S6721;
    float _S6723 = (*dist_coeffs_45)[int(4)] * _S6721.y;
    float _S6724 = (*dist_coeffs_45)[int(5)] * _S6721.x;
    float _S6725 = _S6722.x + _S6722.y;
    float _S6726 = r2_125 * _S6725;
    float _S6727 = r2_125 * _S6726;
    float _S6728 = (*dist_coeffs_45)[int(7)] * _S6721.y + _S6723 + (*dist_coeffs_45)[int(6)] * _S6721.x + _S6724 + _S6393 * _S6725 + _S6392 * _S6726 + _S6391 * _S6727 + (*dist_coeffs_45)[int(3)] * (r2_125 * _S6727);
    float _S6729 = v_130 * _S6728;
    float _S6730 = u_130 * _S6728;
    float3  _S6731 = make_float3 (_S6719.x, _S6719.y, _S6716);
    float2  _S6732 = _S6394 * _S6721 + make_float2 (_S6271 * (v_130 * _S6721.y) + _S6396 * _S6724 + 2.0f * (u_130 * _S6724) + _S6268 * (v_130 * _S6721.x) + _S6730 + _S6730, _S6398 * _S6723 + 2.0f * (v_130 * _S6723) + _S6397 * _S6721.y + _S6395 * _S6721.x + _S6729 + _S6729);
    float2  _S6733 = _S6378 * _S6732;
    float2  _S6734 = _S6389 * _S6732;
    float _S6735 = _S6733.x + _S6733.y;
    if(_S6382)
    {
        float _S6736 = _S6735 / _S6384;
        float _S6737 = _S6385 * - _S6736;
        float _S6738 = _S6381 * (0.3333333432674408f * - (_S6380 * _S6736));
        k_15 = _S6738 + _S6738;
        _S6383 = _S6737;
        _S6384 = 0.0f;
    }
    else
    {
        float _S6739 = _S6735 / _S6383;
        float _S6740 = _S6381 * - _S6739;
        k_15 = _S6379 * _S6739;
        _S6383 = 0.0f;
        _S6384 = _S6740;
    }
    DiffPair_float_0 _S6741;
    (&_S6741)->primal_0 = _S6379;
    (&_S6741)->differential_0 = 0.0f;
    DiffPair_float_0 _S6742;
    (&_S6742)->primal_0 = _S6380;
    (&_S6742)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6741, &_S6742, k_15);
    float _S6743 = _S6742.differential_0 + _S6383;
    float _S6744 = _S6741.differential_0 + _S6384;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6745;
    (&_S6745)->primal_0 = _S6378;
    (&_S6745)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6745, _S6744);
    float2  _S6746 = _S6745.differential_0 + _S6734;
    float _S6747 = fx_41 * _S6673;
    float2  _S6748 = make_float2 (_S6747, fy_41 * _S6672) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6747, (*dist_coeffs_45)[int(9)] * _S6747);
    float2  _S6749 = _S6365 * _S6748;
    float _S6750 = (*dist_coeffs_45)[int(4)] * _S6748.y;
    float _S6751 = (*dist_coeffs_45)[int(5)] * _S6748.x;
    float _S6752 = _S6749.x + _S6749.y;
    float _S6753 = r2_124 * _S6752;
    float _S6754 = r2_124 * _S6753;
    float _S6755 = (*dist_coeffs_45)[int(7)] * _S6748.y + _S6750 + (*dist_coeffs_45)[int(6)] * _S6748.x + _S6751 + _S6368 * _S6752 + _S6367 * _S6753 + _S6366 * _S6754 + (*dist_coeffs_45)[int(3)] * (r2_124 * _S6754);
    float _S6756 = v_129 * _S6755;
    float _S6757 = u_129 * _S6755;
    float3  _S6758 = make_float3 (_S6746.x, _S6746.y, _S6743);
    float2  _S6759 = _S6369 * _S6748 + make_float2 (_S6271 * (v_129 * _S6748.y) + _S6371 * _S6751 + 2.0f * (u_129 * _S6751) + _S6268 * (v_129 * _S6748.x) + _S6757 + _S6757, _S6373 * _S6750 + 2.0f * (v_129 * _S6750) + _S6372 * _S6748.y + _S6370 * _S6748.x + _S6756 + _S6756);
    float2  _S6760 = _S6353 * _S6759;
    float2  _S6761 = _S6364 * _S6759;
    float _S6762 = _S6760.x + _S6760.y;
    if(_S6357)
    {
        float _S6763 = _S6762 / _S6359;
        float _S6764 = _S6360 * - _S6763;
        float _S6765 = _S6356 * (0.3333333432674408f * - (_S6355 * _S6763));
        k_15 = _S6765 + _S6765;
        _S6358 = _S6764;
        _S6359 = 0.0f;
    }
    else
    {
        float _S6766 = _S6762 / _S6358;
        float _S6767 = _S6356 * - _S6766;
        k_15 = _S6354 * _S6766;
        _S6358 = 0.0f;
        _S6359 = _S6767;
    }
    DiffPair_float_0 _S6768;
    (&_S6768)->primal_0 = _S6354;
    (&_S6768)->differential_0 = 0.0f;
    DiffPair_float_0 _S6769;
    (&_S6769)->primal_0 = _S6355;
    (&_S6769)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6768, &_S6769, k_15);
    float _S6770 = _S6769.differential_0 + _S6358;
    float _S6771 = _S6768.differential_0 + _S6359;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6772;
    (&_S6772)->primal_0 = _S6353;
    (&_S6772)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6772, _S6771);
    float2  _S6773 = _S6772.differential_0 + _S6761;
    float _S6774 = fx_41 * _S6675;
    float2  _S6775 = make_float2 (_S6774, fy_41 * _S6674) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6774, (*dist_coeffs_45)[int(9)] * _S6774);
    float2  _S6776 = _S6340 * _S6775;
    float _S6777 = (*dist_coeffs_45)[int(4)] * _S6775.y;
    float _S6778 = (*dist_coeffs_45)[int(5)] * _S6775.x;
    float _S6779 = _S6776.x + _S6776.y;
    float _S6780 = r2_123 * _S6779;
    float _S6781 = r2_123 * _S6780;
    float _S6782 = (*dist_coeffs_45)[int(7)] * _S6775.y + _S6777 + (*dist_coeffs_45)[int(6)] * _S6775.x + _S6778 + _S6343 * _S6779 + _S6342 * _S6780 + _S6341 * _S6781 + (*dist_coeffs_45)[int(3)] * (r2_123 * _S6781);
    float _S6783 = v_128 * _S6782;
    float _S6784 = u_128 * _S6782;
    float3  _S6785 = make_float3 (_S6773.x, _S6773.y, _S6770);
    float2  _S6786 = _S6344 * _S6775 + make_float2 (_S6271 * (v_128 * _S6775.y) + _S6346 * _S6778 + 2.0f * (u_128 * _S6778) + _S6268 * (v_128 * _S6775.x) + _S6784 + _S6784, _S6348 * _S6777 + 2.0f * (v_128 * _S6777) + _S6347 * _S6775.y + _S6345 * _S6775.x + _S6783 + _S6783);
    float2  _S6787 = _S6328 * _S6786;
    float2  _S6788 = _S6339 * _S6786;
    float _S6789 = _S6787.x + _S6787.y;
    if(_S6332)
    {
        float _S6790 = _S6789 / _S6334;
        float _S6791 = _S6335 * - _S6790;
        float _S6792 = _S6331 * (0.3333333432674408f * - (_S6330 * _S6790));
        k_15 = _S6792 + _S6792;
        _S6333 = _S6791;
        _S6334 = 0.0f;
    }
    else
    {
        float _S6793 = _S6789 / _S6333;
        float _S6794 = _S6331 * - _S6793;
        k_15 = _S6329 * _S6793;
        _S6333 = 0.0f;
        _S6334 = _S6794;
    }
    DiffPair_float_0 _S6795;
    (&_S6795)->primal_0 = _S6329;
    (&_S6795)->differential_0 = 0.0f;
    DiffPair_float_0 _S6796;
    (&_S6796)->primal_0 = _S6330;
    (&_S6796)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6795, &_S6796, k_15);
    float _S6797 = _S6796.differential_0 + _S6333;
    float _S6798 = _S6795.differential_0 + _S6334;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6799;
    (&_S6799)->primal_0 = _S6328;
    (&_S6799)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6799, _S6798);
    float2  _S6800 = _S6799.differential_0 + _S6788;
    float _S6801 = fx_41 * _S6666;
    float2  _S6802 = make_float2 (_S6801, fy_41 * _S6676) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6801, (*dist_coeffs_45)[int(9)] * _S6801);
    float2  _S6803 = _S6315 * _S6802;
    float _S6804 = (*dist_coeffs_45)[int(4)] * _S6802.y;
    float _S6805 = (*dist_coeffs_45)[int(5)] * _S6802.x;
    float _S6806 = _S6803.x + _S6803.y;
    float _S6807 = r2_122 * _S6806;
    float _S6808 = r2_122 * _S6807;
    float _S6809 = (*dist_coeffs_45)[int(7)] * _S6802.y + _S6804 + (*dist_coeffs_45)[int(6)] * _S6802.x + _S6805 + _S6318 * _S6806 + _S6317 * _S6807 + _S6316 * _S6808 + (*dist_coeffs_45)[int(3)] * (r2_122 * _S6808);
    float _S6810 = v_127 * _S6809;
    float _S6811 = u_127 * _S6809;
    float3  _S6812 = make_float3 (_S6800.x, _S6800.y, _S6797);
    float2  _S6813 = _S6319 * _S6802 + make_float2 (_S6271 * (v_127 * _S6802.y) + _S6321 * _S6805 + 2.0f * (u_127 * _S6805) + _S6268 * (v_127 * _S6802.x) + _S6811 + _S6811, _S6323 * _S6804 + 2.0f * (v_127 * _S6804) + _S6322 * _S6802.y + _S6320 * _S6802.x + _S6810 + _S6810);
    float2  _S6814 = _S6303 * _S6813;
    float2  _S6815 = _S6314 * _S6813;
    float _S6816 = _S6814.x + _S6814.y;
    if(_S6307)
    {
        float _S6817 = _S6816 / _S6309;
        float _S6818 = _S6310 * - _S6817;
        float _S6819 = _S6306 * (0.3333333432674408f * - (_S6305 * _S6817));
        k_15 = _S6819 + _S6819;
        _S6308 = _S6818;
        _S6309 = 0.0f;
    }
    else
    {
        float _S6820 = _S6816 / _S6308;
        float _S6821 = _S6306 * - _S6820;
        k_15 = _S6304 * _S6820;
        _S6308 = 0.0f;
        _S6309 = _S6821;
    }
    DiffPair_float_0 _S6822;
    (&_S6822)->primal_0 = _S6304;
    (&_S6822)->differential_0 = 0.0f;
    DiffPair_float_0 _S6823;
    (&_S6823)->primal_0 = _S6305;
    (&_S6823)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6822, &_S6823, k_15);
    float _S6824 = _S6823.differential_0 + _S6308;
    float _S6825 = _S6822.differential_0 + _S6309;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6826;
    (&_S6826)->primal_0 = _S6303;
    (&_S6826)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6826, _S6825);
    float2  _S6827 = _S6826.differential_0 + _S6815;
    float _S6828 = fx_41 * _S6663;
    float2  _S6829 = make_float2 (_S6828, fy_41 * _S6665) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6828, (*dist_coeffs_45)[int(9)] * _S6828);
    float2  _S6830 = _S6290 * _S6829;
    float _S6831 = (*dist_coeffs_45)[int(4)] * _S6829.y;
    float _S6832 = (*dist_coeffs_45)[int(5)] * _S6829.x;
    float _S6833 = _S6830.x + _S6830.y;
    float _S6834 = r2_121 * _S6833;
    float _S6835 = r2_121 * _S6834;
    float _S6836 = (*dist_coeffs_45)[int(7)] * _S6829.y + _S6831 + (*dist_coeffs_45)[int(6)] * _S6829.x + _S6832 + _S6293 * _S6833 + _S6292 * _S6834 + _S6291 * _S6835 + (*dist_coeffs_45)[int(3)] * (r2_121 * _S6835);
    float _S6837 = v_126 * _S6836;
    float _S6838 = u_126 * _S6836;
    float3  _S6839 = make_float3 (_S6827.x, _S6827.y, _S6824);
    float2  _S6840 = _S6294 * _S6829 + make_float2 (_S6271 * (v_126 * _S6829.y) + _S6296 * _S6832 + 2.0f * (u_126 * _S6832) + _S6268 * (v_126 * _S6829.x) + _S6838 + _S6838, _S6298 * _S6831 + 2.0f * (v_126 * _S6831) + _S6297 * _S6829.y + _S6295 * _S6829.x + _S6837 + _S6837);
    float2  _S6841 = _S6278 * _S6840;
    float2  _S6842 = _S6289 * _S6840;
    float _S6843 = _S6841.x + _S6841.y;
    if(_S6282)
    {
        float _S6844 = _S6843 / _S6284;
        float _S6845 = _S6285 * - _S6844;
        float _S6846 = _S6281 * (0.3333333432674408f * - (_S6280 * _S6844));
        k_15 = _S6846 + _S6846;
        _S6283 = _S6845;
        _S6284 = 0.0f;
    }
    else
    {
        float _S6847 = _S6843 / _S6283;
        float _S6848 = _S6281 * - _S6847;
        k_15 = _S6279 * _S6847;
        _S6283 = 0.0f;
        _S6284 = _S6848;
    }
    DiffPair_float_0 _S6849;
    (&_S6849)->primal_0 = _S6279;
    (&_S6849)->differential_0 = 0.0f;
    DiffPair_float_0 _S6850;
    (&_S6850)->primal_0 = _S6280;
    (&_S6850)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6849, &_S6850, k_15);
    float _S6851 = _S6850.differential_0 + _S6283;
    float _S6852 = _S6849.differential_0 + _S6284;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6853;
    (&_S6853)->primal_0 = _S6278;
    (&_S6853)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6853, _S6852);
    float2  _S6854 = _S6853.differential_0 + _S6842;
    float _S6855 = fx_41 * _S6667;
    float2  _S6856 = make_float2 (_S6855, fy_41 * _S6664) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6855, (*dist_coeffs_45)[int(9)] * _S6855);
    float2  _S6857 = _S6263 * _S6856;
    float _S6858 = (*dist_coeffs_45)[int(4)] * _S6856.y;
    float _S6859 = (*dist_coeffs_45)[int(5)] * _S6856.x;
    float _S6860 = _S6857.x + _S6857.y;
    float _S6861 = r2_120 * _S6860;
    float _S6862 = r2_120 * _S6861;
    float _S6863 = (*dist_coeffs_45)[int(7)] * _S6856.y + _S6858 + (*dist_coeffs_45)[int(6)] * _S6856.x + _S6859 + _S6266 * _S6860 + _S6265 * _S6861 + _S6264 * _S6862 + (*dist_coeffs_45)[int(3)] * (r2_120 * _S6862);
    float _S6864 = v_125 * _S6863;
    float _S6865 = u_125 * _S6863;
    float3  _S6866 = make_float3 (_S6854.x, _S6854.y, _S6851);
    float2  _S6867 = _S6267 * _S6856 + make_float2 (_S6271 * (v_125 * _S6856.y) + _S6270 * _S6859 + 2.0f * (u_125 * _S6859) + _S6268 * (v_125 * _S6856.x) + _S6865 + _S6865, _S6273 * _S6858 + 2.0f * (v_125 * _S6858) + _S6272 * _S6856.y + _S6269 * _S6856.x + _S6864 + _S6864);
    float2  _S6868 = _S6251 * _S6867;
    float2  _S6869 = _S6262 * _S6867;
    float _S6870 = _S6868.x + _S6868.y;
    if(_S6255)
    {
        float _S6871 = _S6870 / _S6257;
        float _S6872 = _S6258 * - _S6871;
        float _S6873 = _S6254 * (0.3333333432674408f * - (_S6253 * _S6871));
        k_15 = _S6873 + _S6873;
        _S6256 = _S6872;
        _S6257 = 0.0f;
    }
    else
    {
        float _S6874 = _S6870 / _S6256;
        float _S6875 = _S6254 * - _S6874;
        k_15 = _S6252 * _S6874;
        _S6256 = 0.0f;
        _S6257 = _S6875;
    }
    DiffPair_float_0 _S6876;
    (&_S6876)->primal_0 = _S6252;
    (&_S6876)->differential_0 = 0.0f;
    DiffPair_float_0 _S6877;
    (&_S6877)->primal_0 = _S6253;
    (&_S6877)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6876, &_S6877, k_15);
    float _S6878 = _S6877.differential_0 + _S6256;
    float _S6879 = _S6876.differential_0 + _S6257;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6880;
    (&_S6880)->primal_0 = _S6251;
    (&_S6880)->differential_0 = _S6690;
    s_bwd_length_impl_0(&_S6880, _S6879);
    float2  _S6881 = _S6880.differential_0 + _S6869;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6882;
    (&_S6882)->primal_0 = R_38;
    (&_S6882)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6883;
    (&_S6883)->primal_0 = _S6250;
    (&_S6883)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6882, &_S6883, _S6662);
    DiffPair_float_0 _S6884;
    (&_S6884)->primal_0 = _S6247;
    (&_S6884)->differential_0 = 0.0f;
    DiffPair_float_0 _S6885;
    (&_S6885)->primal_0 = _S6249;
    (&_S6885)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6884, &_S6885, 0.0f);
    DiffPair_float_0 _S6886;
    (&_S6886)->primal_0 = _S6246;
    (&_S6886)->differential_0 = 0.0f;
    DiffPair_float_0 _S6887;
    (&_S6887)->primal_0 = _S6249;
    (&_S6887)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6886, &_S6887, 0.0f);
    float _S6888 = _S6885.differential_0 + _S6887.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6889;
    (&_S6889)->primal_0 = _S6248;
    (&_S6889)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6889, _S6888);
    float3  _S6890 = _S6889.differential_0 + _S6704;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6891;
    (&_S6891)->primal_0 = R_38;
    (&_S6891)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6892;
    (&_S6892)->primal_0 = pos_i_13;
    (&_S6892)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6891, &_S6892, _S6890);
    DiffPair_float_0 _S6893;
    (&_S6893)->primal_0 = _S6243;
    (&_S6893)->differential_0 = 0.0f;
    DiffPair_float_0 _S6894;
    (&_S6894)->primal_0 = _S6245;
    (&_S6894)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6893, &_S6894, _S6884.differential_0);
    DiffPair_float_0 _S6895;
    (&_S6895)->primal_0 = _S6242;
    (&_S6895)->differential_0 = 0.0f;
    DiffPair_float_0 _S6896;
    (&_S6896)->primal_0 = _S6245;
    (&_S6896)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6895, &_S6896, _S6886.differential_0);
    float _S6897 = _S6894.differential_0 + _S6896.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6898;
    (&_S6898)->primal_0 = _S6244;
    (&_S6898)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6898, _S6897);
    float3  _S6899 = _S6898.differential_0 + _S6731;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6900;
    (&_S6900)->primal_0 = R_38;
    (&_S6900)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6901;
    (&_S6901)->primal_0 = pos_i_12;
    (&_S6901)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6900, &_S6901, _S6899);
    DiffPair_float_0 _S6902;
    (&_S6902)->primal_0 = _S6239;
    (&_S6902)->differential_0 = 0.0f;
    DiffPair_float_0 _S6903;
    (&_S6903)->primal_0 = _S6241;
    (&_S6903)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6902, &_S6903, _S6893.differential_0);
    DiffPair_float_0 _S6904;
    (&_S6904)->primal_0 = _S6238;
    (&_S6904)->differential_0 = 0.0f;
    DiffPair_float_0 _S6905;
    (&_S6905)->primal_0 = _S6241;
    (&_S6905)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6904, &_S6905, _S6895.differential_0);
    float _S6906 = _S6903.differential_0 + _S6905.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6907;
    (&_S6907)->primal_0 = _S6240;
    (&_S6907)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6907, _S6906);
    float3  _S6908 = _S6907.differential_0 + _S6758;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6909;
    (&_S6909)->primal_0 = R_38;
    (&_S6909)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6910;
    (&_S6910)->primal_0 = pos_i_11;
    (&_S6910)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6909, &_S6910, _S6908);
    DiffPair_float_0 _S6911;
    (&_S6911)->primal_0 = _S6235;
    (&_S6911)->differential_0 = 0.0f;
    DiffPair_float_0 _S6912;
    (&_S6912)->primal_0 = _S6237;
    (&_S6912)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6911, &_S6912, _S6902.differential_0);
    DiffPair_float_0 _S6913;
    (&_S6913)->primal_0 = _S6234;
    (&_S6913)->differential_0 = 0.0f;
    DiffPair_float_0 _S6914;
    (&_S6914)->primal_0 = _S6237;
    (&_S6914)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6913, &_S6914, _S6904.differential_0);
    float _S6915 = _S6912.differential_0 + _S6914.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6916;
    (&_S6916)->primal_0 = _S6236;
    (&_S6916)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6916, _S6915);
    float3  _S6917 = _S6916.differential_0 + _S6785;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6918;
    (&_S6918)->primal_0 = R_38;
    (&_S6918)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6919;
    (&_S6919)->primal_0 = pos_i_10;
    (&_S6919)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6918, &_S6919, _S6917);
    DiffPair_float_0 _S6920;
    (&_S6920)->primal_0 = _S6231;
    (&_S6920)->differential_0 = 0.0f;
    DiffPair_float_0 _S6921;
    (&_S6921)->primal_0 = _S6233;
    (&_S6921)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6920, &_S6921, _S6911.differential_0);
    DiffPair_float_0 _S6922;
    (&_S6922)->primal_0 = _S6230;
    (&_S6922)->differential_0 = 0.0f;
    DiffPair_float_0 _S6923;
    (&_S6923)->primal_0 = _S6233;
    (&_S6923)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6922, &_S6923, _S6913.differential_0);
    float _S6924 = _S6921.differential_0 + _S6923.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6925;
    (&_S6925)->primal_0 = _S6232;
    (&_S6925)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6925, _S6924);
    float3  _S6926 = _S6925.differential_0 + _S6812;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6927;
    (&_S6927)->primal_0 = R_38;
    (&_S6927)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6928;
    (&_S6928)->primal_0 = pos_i_9;
    (&_S6928)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6927, &_S6928, _S6926);
    DiffPair_float_0 _S6929;
    (&_S6929)->primal_0 = _S6227;
    (&_S6929)->differential_0 = 0.0f;
    DiffPair_float_0 _S6930;
    (&_S6930)->primal_0 = _S6229;
    (&_S6930)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6929, &_S6930, _S6920.differential_0);
    DiffPair_float_0 _S6931;
    (&_S6931)->primal_0 = _S6226;
    (&_S6931)->differential_0 = 0.0f;
    DiffPair_float_0 _S6932;
    (&_S6932)->primal_0 = _S6229;
    (&_S6932)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6931, &_S6932, _S6922.differential_0);
    float _S6933 = _S6930.differential_0 + _S6932.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6934;
    (&_S6934)->primal_0 = _S6228;
    (&_S6934)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6934, _S6933);
    float3  _S6935 = _S6934.differential_0 + _S6839;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6936;
    (&_S6936)->primal_0 = R_38;
    (&_S6936)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6937;
    (&_S6937)->primal_0 = pos_i_8;
    (&_S6937)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6936, &_S6937, _S6935);
    DiffPair_float_0 _S6938;
    (&_S6938)->primal_0 = _S6223;
    (&_S6938)->differential_0 = 0.0f;
    DiffPair_float_0 _S6939;
    (&_S6939)->primal_0 = _S6225;
    (&_S6939)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6938, &_S6939, _S6929.differential_0);
    DiffPair_float_0 _S6940;
    (&_S6940)->primal_0 = _S6222;
    (&_S6940)->differential_0 = 0.0f;
    DiffPair_float_0 _S6941;
    (&_S6941)->primal_0 = _S6225;
    (&_S6941)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6940, &_S6941, _S6931.differential_0);
    float _S6942 = _S6939.differential_0 + _S6941.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6943;
    (&_S6943)->primal_0 = _S6224;
    (&_S6943)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6943, _S6942);
    float3  _S6944 = _S6943.differential_0 + _S6866;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6945;
    (&_S6945)->primal_0 = R_38;
    (&_S6945)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6946;
    (&_S6946)->primal_0 = pos_i_7;
    (&_S6946)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6945, &_S6946, _S6944);
    DiffPair_float_0 _S6947;
    (&_S6947)->primal_0 = 0.0f;
    (&_S6947)->differential_0 = 0.0f;
    DiffPair_float_0 _S6948;
    (&_S6948)->primal_0 = _S6221;
    (&_S6948)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6947, &_S6948, _S6938.differential_0);
    DiffPair_float_0 _S6949;
    (&_S6949)->primal_0 = 1.00000001504746622e+30f;
    (&_S6949)->differential_0 = 0.0f;
    DiffPair_float_0 _S6950;
    (&_S6950)->primal_0 = _S6221;
    (&_S6950)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6949, &_S6950, _S6940.differential_0);
    float _S6951 = _S6948.differential_0 + _S6950.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6952;
    (&_S6952)->primal_0 = _S6220;
    (&_S6952)->differential_0 = _S6505;
    s_bwd_length_impl_1(&_S6952, _S6951);
    float3  _S6953 = _S6952.differential_0 + make_float3 (_S6881.x, _S6881.y, _S6878);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6954;
    (&_S6954)->primal_0 = R_38;
    (&_S6954)->differential_0 = _S6567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6955;
    (&_S6955)->primal_0 = pos_5;
    (&_S6955)->differential_0 = _S6505;
    s_bwd_prop_mul_1(&_S6954, &_S6955, _S6953);
    float3  _S6956 = _S6662 + _S6890 + _S6899 + _S6908 + _S6917 + _S6926 + _S6935 + _S6944 + _S6953 + _S6570.differential_0;
    Matrix<float, 3, 3>  _S6957 = _S6882.differential_0 + _S6891.differential_0 + _S6900.differential_0 + _S6909.differential_0 + _S6918.differential_0 + _S6927.differential_0 + _S6936.differential_0 + _S6945.differential_0 + _S6954.differential_0 + _S6571;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_10)[int(0)] = _S6646;
    (*v_sh_coeffs_10)[int(1)] = _S6647;
    (*v_sh_coeffs_10)[int(2)] = _S6648;
    (*v_sh_coeffs_10)[int(3)] = _S6649;
    (*v_sh_coeffs_10)[int(4)] = _S6650;
    (*v_sh_coeffs_10)[int(5)] = _S6651;
    (*v_sh_coeffs_10)[int(6)] = _S6652;
    (*v_sh_coeffs_10)[int(7)] = _S6653;
    (*v_sh_coeffs_10)[int(8)] = _S6654;
    (*v_sh_coeffs_10)[int(9)] = _S6655;
    (*v_sh_coeffs_10)[int(10)] = _S6656;
    (*v_sh_coeffs_10)[int(11)] = _S6657;
    (*v_sh_coeffs_10)[int(12)] = _S6658;
    (*v_sh_coeffs_10)[int(13)] = _S6659;
    (*v_sh_coeffs_10)[int(14)] = _S6660;
    (*v_sh_coeffs_10)[int(15)] = _S6661;
    *v_R_12 = _S6957;
    *v_t_11 = _S6956;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_19, float3  dOut_21)
{
    float3  _S6958 = _slang_select(((*dpx_19).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_19).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_21;
    dpx_19->primal_0 = (*dpx_19).primal_0;
    dpx_19->differential_0 = _S6958;
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
    float3  _S6959 = - (m_2 * (ray_o_10 - center_0));
    float3  ta_0 = _S6959 - k_16;
    float3  tb_0 = _S6959 + k_16;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S6960 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S6960;
    return (*t0_0) < _S6960;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_20, DiffPair_float_0 * dpy_7, DiffPair_float_0 * dps_0, float dOut_22)
{
    float _S6961 = (1.0f - (*dps_0).primal_0) * dOut_22;
    dpx_20->primal_0 = (*dpx_20).primal_0;
    dpx_20->differential_0 = _S6961;
    DiffPair_float_0 _S6962 = *dpy_7;
    float _S6963 = (*dps_0).primal_0 * dOut_22;
    dpy_7->primal_0 = (*dpy_7).primal_0;
    dpy_7->differential_0 = _S6963;
    float _S6964 = (_S6962.primal_0 - (*dpx_20).primal_0) * dOut_22;
    dps_0->primal_0 = _S6962.primal_0;
    dps_0->differential_0 = _S6964;
    return;
}

inline __device__ float lerp_0(float x_70, float y_36, float s_0)
{
    return x_70 + (y_36 - x_70) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_2)
{
    float _S6965 = w_2.z;
    float _S6966 = 1.0f - _S6965;
    float _S6967 = w_2.y;
    float _S6968 = 1.0f - _S6967;
    float _S6969 = _S6966 * _S6968;
    float _S6970 = w_2.x;
    float _S6971 = 1.0f - _S6970;
    float _S6972 = _S6966 * _S6967;
    float _S6973 = _S6965 * _S6968;
    float _S6974 = _S6965 * _S6967;
    return _S6969 * _S6971 * (*densities_6)[int(0)] + _S6969 * _S6970 * (*densities_6)[int(1)] + _S6972 * _S6971 * (*densities_6)[int(2)] + _S6972 * _S6970 * (*densities_6)[int(3)] + _S6973 * _S6971 * (*densities_6)[int(4)] + _S6973 * _S6970 * (*densities_6)[int(5)] + _S6974 * _S6971 * (*densities_6)[int(6)] + _S6974 * _S6970 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  * densities_7, float3  ray_o_11, float3  ray_d_11)
{
    float _S6975 = 0.5f * size_6;
    float3  m_3 = make_float3 (1.0f) / ray_d_11;
    float3  k_17 = abs_0(m_3) * make_float3 (_S6975);
    float3  _S6976 = - (m_3 * (ray_o_11 - (pos_6 + make_float3 (_S6975))));
    float3  ta_1 = _S6976 - k_17;
    float3  tb_1 = _S6976 + k_17;
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
        float _S6977 = interp_0(densities_7, (ray_o_11 + ray_d_11 * make_float3 (lerp_0(t0_1, t1_1, (float(i_12) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6));
        float _S6978;
        if(_S6977 > 1.10000002384185791f)
        {
            _S6978 = _S6977;
        }
        else
        {
            _S6978 = (F32_exp((0.90909093618392944f * _S6977 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S6978;
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
    float _S6979;
};

inline __device__ float3  s_primal_ctx_abs_1(float3  _S6980)
{
    return abs_0(_S6980);
}

inline __device__ float s_primal_ctx_lerp_0(float _S6981, float _S6982, float _S6983)
{
    return lerp_0(_S6981, _S6982, _S6983);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S6984 = dpw_0.z;
    float _S6985 = 1.0f - _S6984;
    float _S6986 = dpw_0.y;
    float _S6987 = 1.0f - _S6986;
    float _S6988 = _S6985 * _S6987;
    float _S6989 = dpw_0.x;
    float _S6990 = 1.0f - _S6989;
    float _S6991 = _S6985 * _S6986;
    float _S6992 = _S6984 * _S6987;
    float _S6993 = _S6984 * _S6986;
    return _S6988 * _S6990 * (*dpdensities_0)[int(0)] + _S6988 * _S6989 * (*dpdensities_0)[int(1)] + _S6991 * _S6990 * (*dpdensities_0)[int(2)] + _S6991 * _S6989 * (*dpdensities_0)[int(3)] + _S6992 * _S6990 * (*dpdensities_0)[int(4)] + _S6992 * _S6989 * (*dpdensities_0)[int(5)] + _S6993 * _S6990 * (*dpdensities_0)[int(6)] + _S6993 * _S6989 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_4, float3  dpray_d_4, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_14)
{
    _s_diff_ctx_14->_S6979 = 0.0f;
    _s_diff_ctx_14->_S6979 = 0.0f;
    float _S6994 = 0.5f * size_7;
    float3  m_4 = make_float3 (1.0f) / dpray_d_4;
    float3  k_18 = s_primal_ctx_abs_1(m_4) * make_float3 (_S6994);
    float3  _S6995 = - (m_4 * (dpray_o_4 - (pos_7 + make_float3 (_S6994))));
    float3  ta_2 = _S6995 - k_18;
    float3  tb_2 = _S6995 + k_18;
    float _S6996 = s_primal_ctx_max_0(s_primal_ctx_max_0(ta_2.x, ta_2.y), s_primal_ctx_max_0(ta_2.z, 0.0f));
    float _S6997 = s_primal_ctx_min_0(s_primal_ctx_min_0(tb_2.x, tb_2.y), tb_2.z);
    float accum_2;
    if(!!(_S6996 < _S6997))
    {
        float _S6998 = - (_S6997 - _S6996);
        bool _runFlag_0 = true;
        accum_2 = 0.0f;
        int i_13 = int(0);
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_14->_S6979 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            float _S6999;
            int _S7000;
            if(i_13 < int(8))
            {
                float _S7001 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_4 + dpray_d_4 * make_float3 (s_primal_ctx_lerp_0(_S6996, _S6997, (float(i_13) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S7001 > 1.10000002384185791f)
                {
                    _S6999 = _S7001;
                }
                else
                {
                    _S6999 = s_primal_ctx_exp_1(0.90909093618392944f * _S7001 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S6999;
                _S7000 = int(2);
                _S6999 = accum_3;
            }
            else
            {
                _S7000 = int(1);
                _S6999 = 0.0f;
            }
            if(_S7000 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                int _S7002 = i_13 + int(1);
                accum_2 = _S6999;
                i_13 = _S7002;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = s_primal_ctx_min_0(1.0f - s_primal_ctx_exp_1(_S6998 / 8.0f * accum_2), 0.99900001287460327f);
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7003, float3  _S7004)
{
    _d_abs_vector_0(_S7003, _S7004);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_14)
{
    float _S7005 = (*dpw_1).primal_0.z;
    float _S7006 = 1.0f - _S7005;
    float _S7007 = (*dpw_1).primal_0.y;
    float _S7008 = 1.0f - _S7007;
    float _S7009 = _S7006 * _S7008;
    float _S7010 = (*dpw_1).primal_0.x;
    float _S7011 = 1.0f - _S7010;
    float _S7012 = _S7006 * _S7007;
    float _S7013 = _S7005 * _S7008;
    float _S7014 = _S7005 * _S7007;
    float _S7015 = _S7014 * _S7010 * _s_dOut_14;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_14;
    float _S7016 = _S7014 * _S7011 * _s_dOut_14;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_14;
    float _S7017 = _S7013 * _S7010 * _s_dOut_14;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_14;
    float _S7018 = _S7013 * _S7011 * _s_dOut_14;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_14;
    float _S7019 = _S7012 * _S7010 * _s_dOut_14;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_14;
    float _S7020 = _S7012 * _S7011 * _s_dOut_14;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_14;
    float _S7021 = _S7009 * _S7010 * _s_dOut_14;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_14;
    float _S7022 = _S7009 * _S7011 * _s_dOut_14;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_14;
    float _S7023 = _S7010 * s_diff_w7_T_0 + _S7011 * s_diff_w6_T_0;
    float _S7024 = _S7010 * s_diff_w5_T_0 + _S7011 * s_diff_w4_T_0;
    float _S7025 = _S7010 * s_diff_w3_T_0 + _S7011 * s_diff_w2_T_0;
    float _S7026 = _S7010 * s_diff_w1_T_0 + _S7011 * s_diff_w0_T_0;
    float3  _S7027 = make_float3 (_S7014 * s_diff_w7_T_0 + _S7013 * s_diff_w5_T_0 + _S7012 * s_diff_w3_T_0 + _S7009 * s_diff_w1_T_0 + - (_S7014 * s_diff_w6_T_0 + _S7013 * s_diff_w4_T_0 + _S7012 * s_diff_w2_T_0 + _S7009 * s_diff_w0_T_0), _S7005 * _S7023 + _S7006 * _S7025 + - (_S7005 * _S7024 + _S7006 * _S7026), _S7007 * _S7023 + _S7008 * _S7024 + - (_S7007 * _S7025 + _S7008 * _S7026));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S7027;
    FixedArray<float, 8>  _S7028;
    _S7028[int(0)] = 0.0f;
    _S7028[int(1)] = 0.0f;
    _S7028[int(2)] = 0.0f;
    _S7028[int(3)] = 0.0f;
    _S7028[int(4)] = 0.0f;
    _S7028[int(5)] = 0.0f;
    _S7028[int(6)] = 0.0f;
    _S7028[int(7)] = 0.0f;
    _S7028[int(7)] = _S7015;
    _S7028[int(6)] = _S7016;
    _S7028[int(5)] = _S7017;
    _S7028[int(4)] = _S7018;
    _S7028[int(3)] = _S7019;
    _S7028[int(2)] = _S7020;
    _S7028[int(1)] = _S7021;
    _S7028[int(0)] = _S7022;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S7028;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S7029, DiffPair_float_0 * _S7030, DiffPair_float_0 * _S7031, float _S7032)
{
    _d_lerp_0(_S7029, _S7030, _S7031, _S7032);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float _s_dOut_15, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_15)
{
    FixedArray<float, 8>  _S7033 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7034 = *dpray_o_5;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7035 = *dpray_d_5;
    float _S7036 = 0.5f * size_8;
    float3  _S7037 = make_float3 (_S7036);
    float3  _S7038 = make_float3 (size_8);
    float3  m_5 = make_float3 (1.0f) / (*dpray_d_5).primal_0;
    float3  _S7039 = (*dpray_d_5).primal_0 * (*dpray_d_5).primal_0;
    float3  _S7040 = (*dpray_o_5).primal_0 - (pos_8 + make_float3 (_S7036));
    float3  k_19 = s_primal_ctx_abs_1(m_5) * make_float3 (_S7036);
    float3  _S7041 = - (m_5 * _S7040);
    float3  ta_3 = _S7041 - k_19;
    float3  tb_3 = _S7041 + k_19;
    float _S7042 = ta_3.x;
    float _S7043 = ta_3.y;
    float _S7044 = s_primal_ctx_max_0(_S7042, _S7043);
    float _S7045 = ta_3.z;
    float _S7046 = s_primal_ctx_max_0(_S7045, 0.0f);
    float _S7047 = s_primal_ctx_max_0(_S7044, _S7046);
    float _S7048 = tb_3.x;
    float _S7049 = tb_3.y;
    float _S7050 = s_primal_ctx_min_0(_S7048, _S7049);
    float _S7051 = tb_3.z;
    float _S7052 = s_primal_ctx_min_0(_S7050, _S7051);
    bool _S7053 = !!(_S7047 < _S7052);
    float _S7054;
    float _S7055;
    float _S7056;
    if(_S7053)
    {
        float _S7057 = - (_S7052 - _S7047) / 8.0f;
        float _S7058 = _S7057 * _s_diff_ctx_15->_S6979;
        _S7054 = 1.0f - s_primal_ctx_exp_1(_S7058);
        _S7055 = _S7058;
        _S7056 = _S7057;
    }
    else
    {
        _S7054 = 0.0f;
        _S7055 = 0.0f;
        _S7056 = 0.0f;
    }
    float3  _S7059 = make_float3 (0.0f);
    FixedArray<float, 8>  _S7060 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    float3  _S7061;
    float3  _S7062;
    FixedArray<float, 8>  _S7063;
    if(_S7053)
    {
        DiffPair_float_0 _S7064;
        (&_S7064)->primal_0 = _S7054;
        (&_S7064)->differential_0 = 0.0f;
        DiffPair_float_0 _S7065;
        (&_S7065)->primal_0 = 0.99900001287460327f;
        (&_S7065)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S7064, &_S7065, _s_dOut_15);
        float _S7066 = - _S7064.differential_0;
        DiffPair_float_0 _S7067;
        (&_S7067)->primal_0 = _S7055;
        (&_S7067)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S7067, _S7066);
        float _S7068 = _S7056 * _S7067.differential_0;
        float _S7069 = 0.125f * (_s_diff_ctx_15->_S6979 * _S7067.differential_0);
        int _dc_0 = int(8);
        _S7054 = _S7068;
        _S7055 = 0.0f;
        _S7056 = 0.0f;
        _S7061 = _S7059;
        _S7062 = _S7059;
        _S7063[int(0)] = 0.0f;
        _S7063[int(1)] = 0.0f;
        _S7063[int(2)] = 0.0f;
        _S7063[int(3)] = 0.0f;
        _S7063[int(4)] = 0.0f;
        _S7063[int(5)] = 0.0f;
        _S7063[int(6)] = 0.0f;
        _S7063[int(7)] = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S7070 = _dc_0 < int(8);
            float _S7071;
            float _S7072;
            int _S7073;
            float3  _S7074;
            float3  _S7075;
            bool _S7076;
            if(_S7070)
            {
                float _S7077 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S7078 = s_primal_ctx_lerp_0(_S7047, _S7052, _S7077);
                float3  _S7079 = make_float3 (_S7078);
                float3  _S7080 = (_S7034.primal_0 + _S7035.primal_0 * make_float3 (_S7078) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S7081 = _S7033;
                float _S7082 = s_primal_ctx_interp_0(&_S7081, _S7080);
                bool _S7083 = _S7082 > 1.10000002384185791f;
                if(_S7083)
                {
                    _S7071 = 0.0f;
                }
                else
                {
                    _S7071 = 0.90909093618392944f * _S7082 - 0.90468984842300415f;
                }
                _S7073 = int(2);
                _S7076 = _S7083;
                _S7074 = _S7080;
                _S7075 = _S7079;
                _S7072 = _S7077;
            }
            else
            {
                _S7073 = int(1);
                _S7076 = false;
                _S7071 = 0.0f;
                _S7074 = _S7059;
                _S7075 = _S7059;
                _S7072 = 0.0f;
            }
            float _S7084;
            float _S7085;
            if(!(_S7073 != int(2)))
            {
                _S7084 = _S7054;
                _S7085 = 0.0f;
            }
            else
            {
                _S7084 = 0.0f;
                _S7085 = _S7054;
            }
            if(_S7070)
            {
                float _S7086 = _S7084 + _S7085;
                float _S7087;
                if(_S7076)
                {
                    _S7087 = _S7084;
                }
                else
                {
                    DiffPair_float_0 _S7088;
                    (&_S7088)->primal_0 = _S7071;
                    (&_S7088)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S7088, _S7084);
                    _S7087 = 0.90909093618392944f * _S7088.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S7089;
                (&_S7089)->primal_0 = _S7033;
                (&_S7089)->differential_0 = _S7060;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S7090;
                (&_S7090)->primal_0 = _S7074;
                (&_S7090)->differential_0 = _S7059;
                s_bwd_prop_interp_0(&_S7089, &_S7090, _S7087);
                float3  _S7091 = _S7090.differential_0 / _S7038;
                float3  _S7092 = _S7035.primal_0 * _S7091;
                float3  _S7093 = _S7075 * _S7091;
                float _S7094 = _S7092.x + _S7092.y + _S7092.z;
                DiffPair_float_0 _S7095;
                (&_S7095)->primal_0 = _S7047;
                (&_S7095)->differential_0 = 0.0f;
                DiffPair_float_0 _S7096;
                (&_S7096)->primal_0 = _S7052;
                (&_S7096)->differential_0 = 0.0f;
                DiffPair_float_0 _S7097;
                (&_S7097)->primal_0 = _S7072;
                (&_S7097)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S7095, &_S7096, &_S7097, _S7094);
                float _S7098 = (&_S7089)->differential_0[int(0)] + _S7063[int(0)];
                float _S7099 = (&_S7089)->differential_0[int(1)] + _S7063[int(1)];
                float _S7100 = (&_S7089)->differential_0[int(2)] + _S7063[int(2)];
                float _S7101 = (&_S7089)->differential_0[int(3)] + _S7063[int(3)];
                float _S7102 = (&_S7089)->differential_0[int(4)] + _S7063[int(4)];
                float _S7103 = (&_S7089)->differential_0[int(5)] + _S7063[int(5)];
                float _S7104 = (&_S7089)->differential_0[int(6)] + _S7063[int(6)];
                float _S7105 = (&_S7089)->differential_0[int(7)] + _S7063[int(7)];
                float3  _S7106 = _S7091 + _S7062;
                float3  _S7107 = _S7093 + _S7061;
                float _S7108 = _S7096.differential_0 + _S7055;
                float _S7109 = _S7095.differential_0 + _S7056;
                _S7054 = _S7086;
                _S7055 = _S7108;
                _S7056 = _S7109;
                _S7061 = _S7107;
                _S7062 = _S7106;
                _S7063[int(0)] = _S7098;
                _S7063[int(1)] = _S7099;
                _S7063[int(2)] = _S7100;
                _S7063[int(3)] = _S7101;
                _S7063[int(4)] = _S7102;
                _S7063[int(5)] = _S7103;
                _S7063[int(6)] = _S7104;
                _S7063[int(7)] = _S7105;
            }
            else
            {
                _S7054 = _S7085;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S7110 = - _S7069;
        float _S7111 = - _S7110 + _S7056;
        _S7054 = _S7110 + _S7055;
        _S7055 = _S7111;
    }
    else
    {
        _S7054 = 0.0f;
        _S7055 = 0.0f;
        _S7061 = _S7059;
        _S7062 = _S7059;
        _S7063[int(0)] = 0.0f;
        _S7063[int(1)] = 0.0f;
        _S7063[int(2)] = 0.0f;
        _S7063[int(3)] = 0.0f;
        _S7063[int(4)] = 0.0f;
        _S7063[int(5)] = 0.0f;
        _S7063[int(6)] = 0.0f;
        _S7063[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S7112;
    (&_S7112)->primal_0 = _S7050;
    (&_S7112)->differential_0 = 0.0f;
    DiffPair_float_0 _S7113;
    (&_S7113)->primal_0 = _S7051;
    (&_S7113)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7112, &_S7113, _S7054);
    DiffPair_float_0 _S7114;
    (&_S7114)->primal_0 = _S7048;
    (&_S7114)->differential_0 = 0.0f;
    DiffPair_float_0 _S7115;
    (&_S7115)->primal_0 = _S7049;
    (&_S7115)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7114, &_S7115, _S7112.differential_0);
    DiffPair_float_0 _S7116;
    (&_S7116)->primal_0 = _S7044;
    (&_S7116)->differential_0 = 0.0f;
    DiffPair_float_0 _S7117;
    (&_S7117)->primal_0 = _S7046;
    (&_S7117)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7116, &_S7117, _S7055);
    DiffPair_float_0 _S7118;
    (&_S7118)->primal_0 = _S7045;
    (&_S7118)->differential_0 = 0.0f;
    DiffPair_float_0 _S7119;
    (&_S7119)->primal_0 = 0.0f;
    (&_S7119)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7118, &_S7119, _S7117.differential_0);
    DiffPair_float_0 _S7120;
    (&_S7120)->primal_0 = _S7042;
    (&_S7120)->differential_0 = 0.0f;
    DiffPair_float_0 _S7121;
    (&_S7121)->primal_0 = _S7043;
    (&_S7121)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7120, &_S7121, _S7116.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S7114.differential_0, _S7115.differential_0, _S7113.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S7120.differential_0, _S7121.differential_0, _S7118.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S7122 = _S7037 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7123;
    (&_S7123)->primal_0 = m_5;
    (&_S7123)->differential_0 = _S7059;
    s_bwd_prop_abs_1(&_S7123, _S7122);
    float3  _S7124 = m_5 * s_diff_n_T_0;
    float3  _S7125 = - ((_S7123.differential_0 + _S7040 * s_diff_n_T_0) / _S7039) + _S7061;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S7125;
    float3  _S7126 = _S7124 + _S7062;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S7126;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S7063;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S7127, float _S7128, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S7129, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7130, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7131, float _S7132)
{
    FixedArray<float, 8>  _S7133 = _S7129->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S7134;
    float _S7135 = s_primal_ctx_evaluate_alpha_voxel_0(_S7127, _S7128, &_S7133, (*_S7130).primal_0, (*_S7131).primal_0, &_S7134);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S7136 = _S7134;
    s_bwd_prop_evaluate_alpha_voxel_0(_S7127, _S7128, _S7129, _S7130, _S7131, _S7132, &_S7136);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  * densities_8, float3  ray_o_12, float3  ray_d_12, float v_alpha_4, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_5, float3  * v_ray_d_5)
{
    FixedArray<float, 8>  _S7137 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = *densities_8;
    (&dp_densities_0)->differential_0 = _S7137;
    float3  _S7138 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_4;
    (&dp_ray_o_4)->primal_0 = ray_o_12;
    (&dp_ray_o_4)->differential_0 = _S7138;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_4;
    (&dp_ray_d_4)->primal_0 = ray_d_12;
    (&dp_ray_d_4)->differential_0 = _S7138;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_4, &dp_ray_d_4, v_alpha_4);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_5 = dp_ray_o_4.differential_0;
    *v_ray_d_5 = dp_ray_d_4.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  * densities_9, float3  rgb_16, float3  ray_o_13, float3  ray_d_13, float3  * out_rgb_1, float * depth_25)
{
    *out_rgb_1 = rgb_16;
    float _S7139 = 0.5f * size_10;
    float3  m_6 = make_float3 (1.0f) / ray_d_13;
    float3  k_20 = abs_0(m_6) * make_float3 (_S7139);
    float3  _S7140 = - (m_6 * (ray_o_13 - (pos_10 + make_float3 (_S7139))));
    float3  ta_4 = _S7140 - k_20;
    float3  tb_4 = _S7140 + k_20;
    float _S7141 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S7142 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
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
        float t_39 = lerp_0(_S7141, _S7142, (float(i_14) + 0.5f) / 8.0f);
        float _S7143 = interp_0(densities_9, (ray_o_13 + ray_d_13 * make_float3 (t_39) - pos_10) / make_float3 (size_10));
        float _S7144;
        if(_S7143 > 1.10000002384185791f)
        {
            _S7144 = _S7143;
        }
        else
        {
            _S7144 = (F32_exp((0.90909093618392944f * _S7143 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S7144;
        float depth_accum_1 = depth_accum_0 + t_39 * _S7144;
        i_14 = i_14 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_25 = (F32_log(((F32_max((depth_accum_0 / accum_4), (0.0f))) + 9.99999997475242708e-07f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S7145;
    float _S7146;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_1, float3  dpray_o_6, float3  dpray_d_6, float3  * dpout_rgb_1, float * dpdepth_3, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_16)
{
    _s_diff_ctx_16->_S7145 = 0.0f;
    _s_diff_ctx_16->_S7146 = 0.0f;
    float _S7147 = 0.5f * size_11;
    float3  m_7 = make_float3 (1.0f) / dpray_d_6;
    float3  k_21 = s_primal_ctx_abs_1(m_7) * make_float3 (_S7147);
    float3  _S7148 = - (m_7 * (dpray_o_6 - (pos_11 + make_float3 (_S7147))));
    float3  ta_5 = _S7148 - k_21;
    float3  tb_5 = _S7148 + k_21;
    float _S7149 = s_primal_ctx_max_0(s_primal_ctx_max_0(ta_5.x, ta_5.y), s_primal_ctx_max_0(ta_5.z, 0.0f));
    float _S7150 = s_primal_ctx_min_0(s_primal_ctx_min_0(tb_5.x, tb_5.y), tb_5.z);
    bool _runFlag_1 = true;
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int i_15 = int(0);
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_16->_S7145 = depth_accum_2;
        _s_diff_ctx_16->_S7146 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        float _S7151;
        float _S7152;
        int _S7153;
        if(i_15 < int(8))
        {
            float _S7154 = s_primal_ctx_lerp_0(_S7149, _S7150, (float(i_15) + 0.5f) / 8.0f);
            float _S7155 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_6 + dpray_d_6 * make_float3 (_S7154) - pos_11) / make_float3 (size_11));
            if(_S7155 > 1.10000002384185791f)
            {
                _S7151 = _S7155;
            }
            else
            {
                _S7151 = s_primal_ctx_exp_1(0.90909093618392944f * _S7155 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S7151;
            float depth_accum_3 = depth_accum_2 + _S7154 * _S7151;
            _S7153 = int(1);
            _S7151 = accum_7;
            _S7152 = depth_accum_3;
        }
        else
        {
            _S7153 = int(0);
            _S7151 = 0.0f;
            _S7152 = 0.0f;
        }
        if(_S7153 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            int _S7156 = i_15 + int(1);
            accum_6 = _S7151;
            depth_accum_2 = _S7152;
            i_15 = _S7156;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S7157 = s_primal_ctx_log_0(s_primal_ctx_max_0(depth_accum_2 / accum_6, 0.0f) + 9.99999997475242708e-07f);
    *dpout_rgb_1 = dprgb_1;
    *dpdepth_3 = _S7157;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_7, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_7, float3  dpout_rgb_2, float dpdepth_4, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_17)
{
    FixedArray<float, 8>  _S7158 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7159 = *dpray_o_7;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7160 = *dpray_d_7;
    float3  _S7161 = make_float3 (size_12);
    float _S7162 = 0.5f * size_12;
    float3  _S7163 = make_float3 (_S7162);
    float3  m_8 = make_float3 (1.0f) / (*dpray_d_7).primal_0;
    float3  _S7164 = (*dpray_d_7).primal_0 * (*dpray_d_7).primal_0;
    float3  _S7165 = (*dpray_o_7).primal_0 - (pos_12 + make_float3 (_S7162));
    float3  k_22 = s_primal_ctx_abs_1(m_8) * make_float3 (_S7162);
    float3  _S7166 = - (m_8 * _S7165);
    float3  ta_6 = _S7166 - k_22;
    float3  tb_6 = _S7166 + k_22;
    float _S7167 = ta_6.x;
    float _S7168 = ta_6.y;
    float _S7169 = s_primal_ctx_max_0(_S7167, _S7168);
    float _S7170 = ta_6.z;
    float _S7171 = s_primal_ctx_max_0(_S7170, 0.0f);
    float _S7172 = s_primal_ctx_max_0(_S7169, _S7171);
    float _S7173 = tb_6.x;
    float _S7174 = tb_6.y;
    float _S7175 = s_primal_ctx_min_0(_S7173, _S7174);
    float _S7176 = tb_6.z;
    float _S7177 = s_primal_ctx_min_0(_S7175, _S7176);
    float _S7178 = _s_diff_ctx_17->_S7145 / _s_diff_ctx_17->_S7146;
    float _S7179 = _s_diff_ctx_17->_S7146 * _s_diff_ctx_17->_S7146;
    float3  _S7180 = make_float3 (0.0f);
    FixedArray<float, 8>  _S7181 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S7182;
    (&_S7182)->primal_0 = s_primal_ctx_max_0(_S7178, 0.0f) + 9.99999997475242708e-07f;
    (&_S7182)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S7182, dpdepth_4);
    DiffPair_float_0 _S7183;
    (&_S7183)->primal_0 = _S7178;
    (&_S7183)->differential_0 = 0.0f;
    DiffPair_float_0 _S7184;
    (&_S7184)->primal_0 = 0.0f;
    (&_S7184)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7183, &_S7184, _S7182.differential_0);
    float _S7185 = _S7183.differential_0 / _S7179;
    float _S7186 = _s_diff_ctx_17->_S7145 * - _S7185;
    float _S7187 = _s_diff_ctx_17->_S7146 * _S7185;
    int _dc_1 = int(8);
    float _S7188 = _S7186;
    float _S7189 = _S7187;
    float _S7190 = 0.0f;
    float _S7191 = 0.0f;
    float3  _S7192 = _S7180;
    float3  _S7193 = _S7180;
    FixedArray<float, 8>  _S7194;
    _S7194[int(0)] = 0.0f;
    _S7194[int(1)] = 0.0f;
    _S7194[int(2)] = 0.0f;
    _S7194[int(3)] = 0.0f;
    _S7194[int(4)] = 0.0f;
    _S7194[int(5)] = 0.0f;
    _S7194[int(6)] = 0.0f;
    _S7194[int(7)] = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S7195 = _dc_1 < int(8);
        int _S7196;
        float _S7197;
        float _S7198;
        float _S7199;
        float _S7200;
        float3  _S7201;
        float3  _S7202;
        bool _S7203;
        if(_S7195)
        {
            float _S7204 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S7205 = s_primal_ctx_lerp_0(_S7172, _S7177, _S7204);
            float3  _S7206 = make_float3 (_S7205);
            float3  _S7207 = (_S7159.primal_0 + _S7160.primal_0 * make_float3 (_S7205) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S7208 = _S7158;
            float _S7209 = s_primal_ctx_interp_0(&_S7208, _S7207);
            bool _S7210 = _S7209 > 1.10000002384185791f;
            if(_S7210)
            {
                _S7197 = _S7209;
                _S7198 = 0.0f;
            }
            else
            {
                float _S7211 = 0.90909093618392944f * _S7209 - 0.90468984842300415f;
                _S7197 = s_primal_ctx_exp_1(_S7211);
                _S7198 = _S7211;
            }
            float _S7212 = _S7197;
            float _S7213 = _S7198;
            _S7196 = int(1);
            _S7197 = _S7205;
            _S7198 = _S7212;
            _S7203 = _S7210;
            _S7199 = _S7213;
            _S7201 = _S7207;
            _S7202 = _S7206;
            _S7200 = _S7204;
        }
        else
        {
            _S7196 = int(0);
            _S7197 = 0.0f;
            _S7198 = 0.0f;
            _S7203 = false;
            _S7199 = 0.0f;
            _S7201 = _S7180;
            _S7202 = _S7180;
            _S7200 = 0.0f;
        }
        float _S7214;
        float _S7215;
        float _S7216;
        float _S7217;
        if(!(_S7196 != int(1)))
        {
            _S7214 = _S7188;
            _S7215 = _S7189;
            _S7216 = 0.0f;
            _S7217 = 0.0f;
        }
        else
        {
            _S7214 = 0.0f;
            _S7215 = 0.0f;
            _S7216 = _S7189;
            _S7217 = _S7188;
        }
        if(_S7195)
        {
            float _S7218 = _S7198 * _S7215;
            float _S7219 = _S7215 + _S7216;
            float _S7220 = _S7197 * _S7215 + _S7214;
            float _S7221 = _S7214 + _S7217;
            float _S7222;
            if(_S7203)
            {
                _S7222 = _S7220;
            }
            else
            {
                DiffPair_float_0 _S7223;
                (&_S7223)->primal_0 = _S7199;
                (&_S7223)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S7223, _S7220);
                _S7222 = 0.90909093618392944f * _S7223.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S7224;
            (&_S7224)->primal_0 = _S7158;
            (&_S7224)->differential_0 = _S7181;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S7225;
            (&_S7225)->primal_0 = _S7201;
            (&_S7225)->differential_0 = _S7180;
            s_bwd_prop_interp_0(&_S7224, &_S7225, _S7222);
            float3  _S7226 = _S7225.differential_0 / _S7161;
            float3  _S7227 = _S7160.primal_0 * _S7226;
            float3  _S7228 = _S7202 * _S7226;
            float _S7229 = _S7227.x + _S7227.y + _S7227.z + _S7218;
            DiffPair_float_0 _S7230;
            (&_S7230)->primal_0 = _S7172;
            (&_S7230)->differential_0 = 0.0f;
            DiffPair_float_0 _S7231;
            (&_S7231)->primal_0 = _S7177;
            (&_S7231)->differential_0 = 0.0f;
            DiffPair_float_0 _S7232;
            (&_S7232)->primal_0 = _S7200;
            (&_S7232)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S7230, &_S7231, &_S7232, _S7229);
            float _S7233 = (&_S7224)->differential_0[int(0)] + _S7194[int(0)];
            float _S7234 = (&_S7224)->differential_0[int(1)] + _S7194[int(1)];
            float _S7235 = (&_S7224)->differential_0[int(2)] + _S7194[int(2)];
            float _S7236 = (&_S7224)->differential_0[int(3)] + _S7194[int(3)];
            float _S7237 = (&_S7224)->differential_0[int(4)] + _S7194[int(4)];
            float _S7238 = (&_S7224)->differential_0[int(5)] + _S7194[int(5)];
            float _S7239 = (&_S7224)->differential_0[int(6)] + _S7194[int(6)];
            float _S7240 = (&_S7224)->differential_0[int(7)] + _S7194[int(7)];
            float3  _S7241 = _S7226 + _S7193;
            float3  _S7242 = _S7228 + _S7192;
            float _S7243 = _S7231.differential_0 + _S7190;
            float _S7244 = _S7230.differential_0 + _S7191;
            _S7188 = _S7221;
            _S7189 = _S7219;
            _S7190 = _S7243;
            _S7191 = _S7244;
            _S7192 = _S7242;
            _S7193 = _S7241;
            _S7194[int(0)] = _S7233;
            _S7194[int(1)] = _S7234;
            _S7194[int(2)] = _S7235;
            _S7194[int(3)] = _S7236;
            _S7194[int(4)] = _S7237;
            _S7194[int(5)] = _S7238;
            _S7194[int(6)] = _S7239;
            _S7194[int(7)] = _S7240;
        }
        else
        {
            _S7188 = _S7217;
            _S7189 = _S7216;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S7245;
    (&_S7245)->primal_0 = _S7175;
    (&_S7245)->differential_0 = 0.0f;
    DiffPair_float_0 _S7246;
    (&_S7246)->primal_0 = _S7176;
    (&_S7246)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7245, &_S7246, _S7190);
    DiffPair_float_0 _S7247;
    (&_S7247)->primal_0 = _S7173;
    (&_S7247)->differential_0 = 0.0f;
    DiffPair_float_0 _S7248;
    (&_S7248)->primal_0 = _S7174;
    (&_S7248)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7247, &_S7248, _S7245.differential_0);
    DiffPair_float_0 _S7249;
    (&_S7249)->primal_0 = _S7169;
    (&_S7249)->differential_0 = 0.0f;
    DiffPair_float_0 _S7250;
    (&_S7250)->primal_0 = _S7171;
    (&_S7250)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7249, &_S7250, _S7191);
    DiffPair_float_0 _S7251;
    (&_S7251)->primal_0 = _S7170;
    (&_S7251)->differential_0 = 0.0f;
    DiffPair_float_0 _S7252;
    (&_S7252)->primal_0 = 0.0f;
    (&_S7252)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7251, &_S7252, _S7250.differential_0);
    DiffPair_float_0 _S7253;
    (&_S7253)->primal_0 = _S7167;
    (&_S7253)->differential_0 = 0.0f;
    DiffPair_float_0 _S7254;
    (&_S7254)->primal_0 = _S7168;
    (&_S7254)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7253, &_S7254, _S7249.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S7247.differential_0, _S7248.differential_0, _S7246.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S7253.differential_0, _S7254.differential_0, _S7251.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S7255 = _S7163 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7256;
    (&_S7256)->primal_0 = m_8;
    (&_S7256)->differential_0 = _S7180;
    s_bwd_prop_abs_1(&_S7256, _S7255);
    float3  _S7257 = m_8 * s_diff_n_T_1;
    float3  _S7258 = - ((_S7256.differential_0 + _S7165 * s_diff_n_T_1) / _S7164) + _S7192;
    dpray_d_7->primal_0 = (*dpray_d_7).primal_0;
    dpray_d_7->differential_0 = _S7258;
    float3  _S7259 = _S7257 + _S7193;
    dpray_o_7->primal_0 = (*dpray_o_7).primal_0;
    dpray_o_7->differential_0 = _S7259;
    dprgb_2->primal_0 = (*dprgb_2).primal_0;
    dprgb_2->differential_0 = dpout_rgb_2;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S7194;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S7260, float _S7261, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S7262, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7263, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7264, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7265, float3  _S7266, float _S7267)
{
    FixedArray<float, 8>  _S7268 = _S7262->primal_0;
    float3  _S7269;
    float _S7270;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S7271;
    s_primal_ctx_evaluate_color_voxel_0(_S7260, _S7261, &_S7268, (*_S7263).primal_0, (*_S7264).primal_0, (*_S7265).primal_0, &_S7269, &_S7270, &_S7271);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S7272 = _S7271;
    s_bwd_prop_evaluate_color_voxel_0(_S7260, _S7261, _S7262, _S7263, _S7264, _S7265, _S7266, _S7267, &_S7272);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  * densities_10, float3  rgb_17, float3  ray_o_14, float3  ray_d_14, float3  v_out_rgb_1, float v_depth_12, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_10, float3  * v_ray_o_6, float3  * v_ray_d_6)
{
    FixedArray<float, 8>  _S7273 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = *densities_10;
    (&dp_densities_1)->differential_0 = _S7273;
    float3  _S7274 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_1;
    (&dp_rgb_1)->primal_0 = rgb_17;
    (&dp_rgb_1)->differential_0 = _S7274;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_5;
    (&dp_ray_o_5)->primal_0 = ray_o_14;
    (&dp_ray_o_5)->differential_0 = _S7274;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_5;
    (&dp_ray_d_5)->primal_0 = ray_d_14;
    (&dp_ray_d_5)->differential_0 = _S7274;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_1, &dp_ray_o_5, &dp_ray_d_5, v_out_rgb_1, v_depth_12);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_10 = dp_rgb_1.differential_0;
    *v_ray_o_6 = dp_ray_o_5.differential_0;
    *v_ray_d_6 = dp_ray_d_5.differential_0;
    return;
}

