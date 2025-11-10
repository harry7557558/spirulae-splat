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

struct s_bwd_prop_persp_proj_Intermediates_0
{
    bool _S174;
};

inline __device__ void persp_proj_3dgs(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  * dist_coeffs_8, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float2  _S175;
    bool _S176;
    float2  _S177;
    float _S178;
    bool _S179;
    for(;;)
    {
        float2  _S180 = float2 {mean3d_1.x, mean3d_1.y};
        _S177 = _S180;
        float _S181 = mean3d_1.z;
        _S178 = _S181;
        float2  uv_8 = _S180 / make_float2 (_S181);
        bool _S182 = _S181 < 0.0f;
        _S179 = _S182;
        if(_S182)
        {
            _S176 = true;
        }
        else
        {
            bool _S183 = is_valid_distortion(uv_8, dist_coeffs_8);
            _S176 = !_S183;
        }
        if(_S176)
        {
            _S175 = make_float2 (-10000.0f);
            break;
        }
        float u_4 = uv_8.x;
        float v_4 = uv_8.y;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float2  _S184 = uv_8 * make_float2 (1.0f + r2_4 * ((*dist_coeffs_8)[int(0)] + r2_4 * ((*dist_coeffs_8)[int(1)] + r2_4 * ((*dist_coeffs_8)[int(2)] + r2_4 * (*dist_coeffs_8)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_8)[int(4)] * u_4 * v_4 + (*dist_coeffs_8)[int(5)] * (r2_4 + 2.0f * u_4 * u_4) + (*dist_coeffs_8)[int(6)] * r2_4, 2.0f * (*dist_coeffs_8)[int(5)] * u_4 * v_4 + (*dist_coeffs_8)[int(4)] * (r2_4 + 2.0f * v_4 * v_4) + (*dist_coeffs_8)[int(7)] * r2_4);
        float2  _S185 = _S184 + make_float2 ((*dist_coeffs_8)[int(8)] * _S184.x + (*dist_coeffs_8)[int(9)] * _S184.y, 0.0f);
        _S175 = make_float2 (fx_1 * _S185.x + cx_1, fy_1 * _S185.y + cy_1);
        break;
    }
    *mean2d_1 = _S175;
    Matrix<float, 2, 3>  J_3;
    float2  _S186 = make_float2 (0.0f);
    float2  seed_4 = _S186;
    *&((&seed_4)->x) = 1.0f;
    float2  _S187 = seed_4;
    s_bwd_prop_persp_proj_Intermediates_0 _S188;
    (&_S188)->_S174 = false;
    (&_S188)->_S174 = false;
    float2  uv_9 = _S177 / make_float2 (_S178);
    if(_S179)
    {
    }
    else
    {
        bool _S189 = is_valid_distortion(uv_9, dist_coeffs_8);
        (&_S188)->_S174 = _S189;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S190 = _S188;
    float2  _S191 = make_float2 (0.0f);
    float2  _S192 = make_float2 (_S178);
    float2  uv_10 = _S177 / make_float2 (_S178);
    float2  _S193 = make_float2 (_S178 * _S178);
    if(_S179)
    {
        _S176 = true;
    }
    else
    {
        _S176 = !_S190._S174;
    }
    bool _S194 = !_S176;
    float _S195;
    float _S196;
    float _S197;
    float _S198;
    float _S199;
    float _S200;
    float _S201;
    float _S202;
    float _S203;
    float _S204;
    float _S205;
    float _S206;
    float _S207;
    float _S208;
    float _S209;
    float _S210;
    float _S211;
    float _S212;
    float _S213;
    float _S214;
    float _S215;
    if(_S194)
    {
        float u_5 = uv_10.x;
        float v_5 = uv_10.y;
        float r2_5 = u_5 * u_5 + v_5 * v_5;
        float _S216 = (*dist_coeffs_8)[int(2)] + r2_5 * (*dist_coeffs_8)[int(3)];
        float _S217 = (*dist_coeffs_8)[int(1)] + r2_5 * _S216;
        float _S218 = (*dist_coeffs_8)[int(0)] + r2_5 * _S217;
        float2  _S219 = make_float2 (1.0f + r2_5 * _S218);
        float _S220 = 2.0f * (*dist_coeffs_8)[int(4)];
        float _S221 = _S220 * u_5;
        float _S222 = 2.0f * u_5;
        float _S223 = 2.0f * (*dist_coeffs_8)[int(5)];
        float _S224 = _S223 * u_5;
        float _S225 = 2.0f * v_5;
        _S195 = fy_1;
        _S196 = fx_1;
        _S197 = (*dist_coeffs_8)[int(9)];
        _S198 = (*dist_coeffs_8)[int(8)];
        _S175 = _S219;
        _S199 = (*dist_coeffs_8)[int(7)];
        _S200 = (*dist_coeffs_8)[int(4)];
        _S201 = _S225;
        _S202 = v_5;
        _S203 = _S224;
        _S204 = _S223;
        _S205 = (*dist_coeffs_8)[int(6)];
        _S206 = (*dist_coeffs_8)[int(5)];
        _S207 = _S222;
        _S208 = u_5;
        _S209 = _S221;
        _S210 = _S220;
        _S211 = r2_5;
        _S212 = _S218;
        _S213 = _S217;
        _S214 = _S216;
        _S215 = (*dist_coeffs_8)[int(3)];
    }
    else
    {
        _S195 = 0.0f;
        _S196 = 0.0f;
        _S197 = 0.0f;
        _S198 = 0.0f;
        _S175 = _S191;
        _S199 = 0.0f;
        _S200 = 0.0f;
        _S201 = 0.0f;
        _S202 = 0.0f;
        _S203 = 0.0f;
        _S204 = 0.0f;
        _S205 = 0.0f;
        _S206 = 0.0f;
        _S207 = 0.0f;
        _S208 = 0.0f;
        _S209 = 0.0f;
        _S210 = 0.0f;
        _S211 = 0.0f;
        _S212 = 0.0f;
        _S213 = 0.0f;
        _S214 = 0.0f;
        _S215 = 0.0f;
    }
    if(_S194)
    {
        float _S226 = _S196 * _S187.x;
        float2  _S227 = make_float2 (_S226, _S195 * _S187.y) + make_float2 (_S198 * _S226, _S197 * _S226);
        float2  _S228 = uv_10 * _S227;
        float _S229 = _S200 * _S227.y;
        float _S230 = _S206 * _S227.x;
        float _S231 = _S228.x + _S228.y;
        float _S232 = _S211 * _S231;
        float _S233 = _S211 * _S232;
        float _S234 = _S199 * _S227.y + _S229 + _S205 * _S227.x + _S230 + _S212 * _S231 + _S213 * _S232 + _S214 * _S233 + _S215 * (_S211 * _S233);
        float _S235 = _S202 * _S234;
        float _S236 = _S208 * _S234;
        _S175 = _S175 * _S227 + make_float2 (_S204 * (_S202 * _S227.y) + _S207 * _S230 + 2.0f * (_S208 * _S230) + _S210 * (_S202 * _S227.x) + _S236 + _S236, _S201 * _S229 + 2.0f * (_S202 * _S229) + _S203 * _S227.y + _S209 * _S227.x + _S235 + _S235);
    }
    else
    {
        _S175 = _S191;
    }
    float2  _S237 = _S175 / _S193;
    float2  _S238 = _S177 * - _S237;
    float2  _S239 = _S192 * _S237;
    J_3[int(0)] = make_float3 (_S239.x, _S239.y, _S238.x + _S238.y);
    float2  seed_5 = _S186;
    *&((&seed_5)->y) = 1.0f;
    float2  _S240 = seed_5;
    s_bwd_prop_persp_proj_Intermediates_0 _S241;
    (&_S241)->_S174 = false;
    (&_S241)->_S174 = false;
    float2  uv_11 = _S177 / make_float2 (_S178);
    if(_S179)
    {
    }
    else
    {
        bool _S242 = is_valid_distortion(uv_11, dist_coeffs_8);
        (&_S241)->_S174 = _S242;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S243 = _S241;
    float2  uv_12 = _S177 / make_float2 (_S178);
    if(_S179)
    {
        _S176 = true;
    }
    else
    {
        _S176 = !_S243._S174;
    }
    bool _S244 = !_S176;
    if(_S244)
    {
        float u_6 = uv_12.x;
        float v_6 = uv_12.y;
        float r2_6 = u_6 * u_6 + v_6 * v_6;
        float _S245 = (*dist_coeffs_8)[int(2)] + r2_6 * (*dist_coeffs_8)[int(3)];
        float _S246 = (*dist_coeffs_8)[int(1)] + r2_6 * _S245;
        float _S247 = (*dist_coeffs_8)[int(0)] + r2_6 * _S246;
        float2  _S248 = make_float2 (1.0f + r2_6 * _S247);
        float _S249 = 2.0f * (*dist_coeffs_8)[int(4)];
        float _S250 = _S249 * u_6;
        float _S251 = 2.0f * u_6;
        float _S252 = 2.0f * (*dist_coeffs_8)[int(5)];
        float _S253 = _S252 * u_6;
        float _S254 = 2.0f * v_6;
        _S195 = fy_1;
        _S196 = fx_1;
        _S197 = (*dist_coeffs_8)[int(9)];
        _S198 = (*dist_coeffs_8)[int(8)];
        _S175 = _S248;
        _S199 = (*dist_coeffs_8)[int(7)];
        _S200 = (*dist_coeffs_8)[int(4)];
        _S201 = _S254;
        _S202 = v_6;
        _S203 = _S253;
        _S204 = _S252;
        _S205 = (*dist_coeffs_8)[int(6)];
        _S206 = (*dist_coeffs_8)[int(5)];
        _S207 = _S251;
        _S208 = u_6;
        _S209 = _S250;
        _S210 = _S249;
        _S211 = r2_6;
        _S212 = _S247;
        _S213 = _S246;
        _S214 = _S245;
        _S215 = (*dist_coeffs_8)[int(3)];
    }
    else
    {
        _S195 = 0.0f;
        _S196 = 0.0f;
        _S197 = 0.0f;
        _S198 = 0.0f;
        _S175 = _S191;
        _S199 = 0.0f;
        _S200 = 0.0f;
        _S201 = 0.0f;
        _S202 = 0.0f;
        _S203 = 0.0f;
        _S204 = 0.0f;
        _S205 = 0.0f;
        _S206 = 0.0f;
        _S207 = 0.0f;
        _S208 = 0.0f;
        _S209 = 0.0f;
        _S210 = 0.0f;
        _S211 = 0.0f;
        _S212 = 0.0f;
        _S213 = 0.0f;
        _S214 = 0.0f;
        _S215 = 0.0f;
    }
    if(_S244)
    {
        float _S255 = _S196 * _S240.x;
        float2  _S256 = make_float2 (_S255, _S195 * _S240.y) + make_float2 (_S198 * _S255, _S197 * _S255);
        float2  _S257 = uv_12 * _S256;
        float _S258 = _S200 * _S256.y;
        float _S259 = _S206 * _S256.x;
        float _S260 = _S257.x + _S257.y;
        float _S261 = _S211 * _S260;
        float _S262 = _S211 * _S261;
        float _S263 = _S199 * _S256.y + _S258 + _S205 * _S256.x + _S259 + _S212 * _S260 + _S213 * _S261 + _S214 * _S262 + _S215 * (_S211 * _S262);
        float _S264 = _S202 * _S263;
        float _S265 = _S208 * _S263;
        _S175 = _S175 * _S256 + make_float2 (_S204 * (_S202 * _S256.y) + _S207 * _S259 + 2.0f * (_S208 * _S259) + _S210 * (_S202 * _S256.x) + _S265 + _S265, _S201 * _S258 + 2.0f * (_S202 * _S258) + _S203 * _S256.y + _S209 * _S256.x + _S264 + _S264);
    }
    else
    {
        _S175 = _S191;
    }
    float2  _S266 = _S175 / _S193;
    float2  _S267 = _S177 * - _S266;
    float2  _S268 = _S192 * _S266;
    J_3[int(1)] = make_float3 (_S268.x, _S268.y, _S267.x + _S267.y);
    *cov2d_1 = mul_6(mul_5(J_3, cov3d_1), transpose_1(J_3));
    return;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_7, float dOut_10)
{
    DiffPair_float_0 _S269 = *dpx_7;
    float _S270 = - (*dpy_4).primal_0 / ((*dpx_7).primal_0 * (*dpx_7).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S270;
    float _S271 = _S269.primal_0 / (_S269.primal_0 * _S269.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S271;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S272, float _S273)
{
    return (F32_atan2((_S272), (_S273)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S274, DiffPair_float_0 * _S275, float _S276)
{
    _d_atan2_0(_S274, _S275, _S276);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8, float _s_dOut_2)
{
    float _S277 = (*dpx_8).primal_0.x;
    float _S278 = (*dpx_8).primal_0.y;
    DiffPair_float_0 _S279;
    (&_S279)->primal_0 = _S277 * _S277 + _S278 * _S278;
    (&_S279)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S279, _s_dOut_2);
    float _S280 = (*dpx_8).primal_0.y * _S279.differential_0;
    float _S281 = _S280 + _S280;
    float _S282 = (*dpx_8).primal_0.x * _S279.differential_0;
    float _S283 = _S282 + _S282;
    float2  _S284 = make_float2 (0.0f);
    *&((&_S284)->y) = _S281;
    *&((&_S284)->x) = _S283;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S284;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S285, float _S286)
{
    s_bwd_prop_length_impl_1(_S285, _S286);
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  * dist_coeffs_9, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    float2  _S287 = float2 {mean3d_2.x, mean3d_2.y};
    float r_8 = length_0(_S287);
    float _S288 = mean3d_2.z;
    float theta_3 = (F32_atan2((r_8), (_S288)));
    float k_0;
    if(theta_3 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S288;
    }
    else
    {
        k_0 = theta_3 / r_8;
    }
    float2  _S289 = _S287 * make_float2 (k_0);
    float u_7 = _S289.x;
    float v_7 = _S289.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S290 = 2.0f * (*dist_coeffs_9)[int(4)];
    float _S291 = 2.0f * (*dist_coeffs_9)[int(5)];
    float2  _S292 = _S289 * make_float2 (1.0f + r2_7 * ((*dist_coeffs_9)[int(0)] + r2_7 * ((*dist_coeffs_9)[int(1)] + r2_7 * ((*dist_coeffs_9)[int(2)] + r2_7 * (*dist_coeffs_9)[int(3)])))) + make_float2 (_S290 * u_7 * v_7 + (*dist_coeffs_9)[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + (*dist_coeffs_9)[int(6)] * r2_7, _S291 * u_7 * v_7 + (*dist_coeffs_9)[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + (*dist_coeffs_9)[int(7)] * r2_7);
    float2  _S293 = _S292 + make_float2 ((*dist_coeffs_9)[int(8)] * _S292.x + (*dist_coeffs_9)[int(9)] * _S292.y, 0.0f);
    *mean2d_2 = make_float2 (fx_2 * _S293.x + cx_2, fy_2 * _S293.y + cy_2);
    Matrix<float, 2, 3>  J_4;
    float2  _S294 = make_float2 (0.0f);
    float2  seed_6 = _S294;
    *&((&seed_6)->x) = 1.0f;
    float2  _S295 = seed_6;
    float _S296 = s_primal_ctx_atan2_0(r_8, _S288);
    bool _S297 = _S296 < 0.00100000004749745f;
    float _S298;
    float _S299;
    float _S300;
    if(_S297)
    {
        float _S301 = 1.0f - _S296 * _S296 / 3.0f;
        float _S302 = _S288 * _S288;
        k_0 = _S301 / _S288;
        _S298 = 0.0f;
        _S299 = _S302;
        _S300 = _S301;
    }
    else
    {
        float _S303 = r_8 * r_8;
        k_0 = _S296 / r_8;
        _S298 = _S303;
        _S299 = 0.0f;
        _S300 = 0.0f;
    }
    float2  _S304 = make_float2 (k_0);
    float2  _S305 = _S287 * make_float2 (k_0);
    float u_8 = _S305.x;
    float v_8 = _S305.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S306 = (*dist_coeffs_9)[int(2)] + r2_8 * (*dist_coeffs_9)[int(3)];
    float _S307 = (*dist_coeffs_9)[int(1)] + r2_8 * _S306;
    float _S308 = (*dist_coeffs_9)[int(0)] + r2_8 * _S307;
    float _S309 = fx_2 * _S295.x;
    float2  _S310 = make_float2 (_S309, fy_2 * _S295.y) + make_float2 ((*dist_coeffs_9)[int(8)] * _S309, (*dist_coeffs_9)[int(9)] * _S309);
    float2  _S311 = _S305 * _S310;
    float _S312 = (*dist_coeffs_9)[int(4)] * _S310.y;
    float _S313 = (*dist_coeffs_9)[int(5)] * _S310.x;
    float _S314 = _S311.x + _S311.y;
    float _S315 = r2_8 * _S314;
    float _S316 = r2_8 * _S315;
    float _S317 = (*dist_coeffs_9)[int(7)] * _S310.y + _S312 + (*dist_coeffs_9)[int(6)] * _S310.x + _S313 + _S308 * _S314 + _S307 * _S315 + _S306 * _S316 + (*dist_coeffs_9)[int(3)] * (r2_8 * _S316);
    float _S318 = v_8 * _S317;
    float _S319 = u_8 * _S317;
    float2  _S320 = make_float2 (1.0f + r2_8 * _S308) * _S310 + make_float2 (_S291 * (v_8 * _S310.y) + 2.0f * u_8 * _S313 + 2.0f * (u_8 * _S313) + _S290 * (v_8 * _S310.x) + _S319 + _S319, 2.0f * v_8 * _S312 + 2.0f * (v_8 * _S312) + _S291 * u_8 * _S310.y + _S290 * u_8 * _S310.x + _S318 + _S318);
    float2  _S321 = _S287 * _S320;
    float2  _S322 = _S304 * _S320;
    float _S323 = _S321.x + _S321.y;
    if(_S297)
    {
        float _S324 = _S323 / _S299;
        float _S325 = _S300 * - _S324;
        float _S326 = _S296 * (0.3333333432674408f * - (_S288 * _S324));
        k_0 = _S326 + _S326;
        _S298 = _S325;
        _S299 = 0.0f;
    }
    else
    {
        float _S327 = _S323 / _S298;
        float _S328 = _S296 * - _S327;
        k_0 = r_8 * _S327;
        _S298 = 0.0f;
        _S299 = _S328;
    }
    DiffPair_float_0 _S329;
    (&_S329)->primal_0 = r_8;
    (&_S329)->differential_0 = 0.0f;
    DiffPair_float_0 _S330;
    (&_S330)->primal_0 = _S288;
    (&_S330)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S329, &_S330, k_0);
    float _S331 = _S330.differential_0 + _S298;
    float _S332 = _S329.differential_0 + _S299;
    float2  _S333 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S334;
    (&_S334)->primal_0 = _S287;
    (&_S334)->differential_0 = _S333;
    s_bwd_length_impl_1(&_S334, _S332);
    float2  _S335 = _S334.differential_0 + _S322;
    float3  _S336 = make_float3 (_S335.x, _S335.y, _S331);
    J_4[int(0)] = _S336;
    float2  seed_7 = _S294;
    *&((&seed_7)->y) = 1.0f;
    float2  _S337 = seed_7;
    if(_S297)
    {
        float _S338 = 1.0f - _S296 * _S296 / 3.0f;
        float _S339 = _S288 * _S288;
        k_0 = _S338 / _S288;
        _S298 = 0.0f;
        _S299 = _S339;
        _S300 = _S338;
    }
    else
    {
        float _S340 = r_8 * r_8;
        k_0 = _S296 / r_8;
        _S298 = _S340;
        _S299 = 0.0f;
        _S300 = 0.0f;
    }
    float2  _S341 = make_float2 (k_0);
    float2  _S342 = _S287 * make_float2 (k_0);
    float u_9 = _S342.x;
    float v_9 = _S342.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S343 = (*dist_coeffs_9)[int(2)] + r2_9 * (*dist_coeffs_9)[int(3)];
    float _S344 = (*dist_coeffs_9)[int(1)] + r2_9 * _S343;
    float _S345 = (*dist_coeffs_9)[int(0)] + r2_9 * _S344;
    float _S346 = fx_2 * _S337.x;
    float2  _S347 = make_float2 (_S346, fy_2 * _S337.y) + make_float2 ((*dist_coeffs_9)[int(8)] * _S346, (*dist_coeffs_9)[int(9)] * _S346);
    float2  _S348 = _S342 * _S347;
    float _S349 = (*dist_coeffs_9)[int(4)] * _S347.y;
    float _S350 = (*dist_coeffs_9)[int(5)] * _S347.x;
    float _S351 = _S348.x + _S348.y;
    float _S352 = r2_9 * _S351;
    float _S353 = r2_9 * _S352;
    float _S354 = (*dist_coeffs_9)[int(7)] * _S347.y + _S349 + (*dist_coeffs_9)[int(6)] * _S347.x + _S350 + _S345 * _S351 + _S344 * _S352 + _S343 * _S353 + (*dist_coeffs_9)[int(3)] * (r2_9 * _S353);
    float _S355 = v_9 * _S354;
    float _S356 = u_9 * _S354;
    float2  _S357 = make_float2 (1.0f + r2_9 * _S345) * _S347 + make_float2 (_S291 * (v_9 * _S347.y) + 2.0f * u_9 * _S350 + 2.0f * (u_9 * _S350) + _S290 * (v_9 * _S347.x) + _S356 + _S356, 2.0f * v_9 * _S349 + 2.0f * (v_9 * _S349) + _S291 * u_9 * _S347.y + _S290 * u_9 * _S347.x + _S355 + _S355);
    float2  _S358 = _S287 * _S357;
    float2  _S359 = _S341 * _S357;
    float _S360 = _S358.x + _S358.y;
    if(_S297)
    {
        float _S361 = _S360 / _S299;
        float _S362 = _S300 * - _S361;
        float _S363 = _S296 * (0.3333333432674408f * - (_S288 * _S361));
        k_0 = _S363 + _S363;
        _S298 = _S362;
        _S299 = 0.0f;
    }
    else
    {
        float _S364 = _S360 / _S298;
        float _S365 = _S296 * - _S364;
        k_0 = r_8 * _S364;
        _S298 = 0.0f;
        _S299 = _S365;
    }
    DiffPair_float_0 _S366;
    (&_S366)->primal_0 = r_8;
    (&_S366)->differential_0 = 0.0f;
    DiffPair_float_0 _S367;
    (&_S367)->primal_0 = _S288;
    (&_S367)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S366, &_S367, k_0);
    float _S368 = _S367.differential_0 + _S298;
    float _S369 = _S366.differential_0 + _S299;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S370;
    (&_S370)->primal_0 = _S287;
    (&_S370)->differential_0 = _S333;
    s_bwd_length_impl_1(&_S370, _S369);
    float2  _S371 = _S370.differential_0 + _S359;
    float3  _S372 = make_float3 (_S371.x, _S371.y, _S368);
    J_4[int(1)] = _S372;
    *cov2d_2 = mul_6(mul_5(J_4, cov3d_2), transpose_1(J_4));
    return;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_3, float fy_3, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_3, 0.0f, 0.0f, 0.0f, fy_3, 0.0f);
    *cov2d_3 = mul_6(mul_5(J_5, cov3d_3), transpose_1(J_5));
    *mean2d_3 = make_float2 (fx_3 * mean3d_3.x + cx_3, fy_3 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    float _S373 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S373;
    float _S374 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S374;
    float det_blur_0 = _S373 * _S374 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_9, float dOut_11)
{
    float _S375 = (F32_exp(((*dpx_9).primal_0))) * dOut_11;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S375;
    return;
}

inline __device__ float3  exp_0(float3  x_12)
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
        *_slang_vector_get_element_ptr(&result_12, i_8) = (F32_exp((_slang_vector_get_element(x_12, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_12;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_12)
{
    float3  _S376 = exp_0((*dpx_10).primal_0) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S376;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_11, float dOut_13)
{
    float _S377 = 1.0f / (*dpx_11).primal_0 * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S377;
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
        *_slang_vector_get_element_ptr(&result_13, i_9) = (F32_max((_slang_vector_get_element(x_13, i_9)), (_slang_vector_get_element(y_2, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_13;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  * dist_coeffs_10, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S378 = mean_c_0.z;
        bool _S379;
        if(_S378 < near_plane_0)
        {
            _S379 = true;
        }
        else
        {
            _S379 = _S378 > far_plane_0;
        }
        if(_S379)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S380 = exp_0(scale_2);
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
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S380.x, 0.0f, 0.0f, 0.0f, _S380.y, 0.0f, 0.0f, 0.0f, _S380.z));
        Matrix<float, 3, 3>  _S381 = transpose_0(R_4);
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S381);
        Matrix<float, 2, 2>  covar2d_0;
        persp_proj_3dgs(mean_c_0, covar_c_0, fx_4, fy_4, cx_4, cy_4, dist_coeffs_10, &covar2d_0, mean2d_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S382 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S382;
        float _S383 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S383;
        float det_blur_1 = _S382 * _S383 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S384 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S379 = true;
        }
        else
        {
            _S379 = xmin_0 >= float(image_width_0);
        }
        if(_S379)
        {
            _S379 = true;
        }
        else
        {
            _S379 = ymax_0 <= 0.0f;
        }
        if(_S379)
        {
            _S379 = true;
        }
        else
        {
            _S379 = ymin_0 >= float(image_height_0);
        }
        if(_S379)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S384.rows[int(0)].x, _S384.rows[int(0)].y, _S384.rows[int(1)].y);
        float3  _S385 = mean_0 - - mul_0(_S381, t_3);
        float3  _S386 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S386;
        float _S387 = _S385.x;
        float _S388 = _S385.y;
        float _S389 = _S385.z;
        float norm_0 = (F32_sqrt((_S387 * _S387 + _S388 * _S388 + _S389 * _S389)));
        float x_15 = _S387 / norm_0;
        float y_3 = _S388 / norm_0;
        float z_0 = _S389 / norm_0;
        float3  _S390 = _S386 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_15) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S390;
        float z2_4 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_15 * x_15 - y_3 * y_3;
        float fS1_0 = 2.0f * x_15 * y_3;
        float3  _S391 = _S390 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_15) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S391;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S391 + (make_float3 (-0.59004360437393188f * (x_15 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_15) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  * dist_coeffs_11, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S392 = length_1(mean_c_1);
        bool _S393;
        if(_S392 < near_plane_1)
        {
            _S393 = true;
        }
        else
        {
            _S393 = _S392 > far_plane_1;
        }
        if(_S393)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S394 = exp_0(scale_3);
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
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S394.x, 0.0f, 0.0f, 0.0f, _S394.y, 0.0f, 0.0f, 0.0f, _S394.z));
        Matrix<float, 3, 3>  _S395 = transpose_0(R_5);
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S395);
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, dist_coeffs_11, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S396 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S396;
        float _S397 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S397;
        float det_blur_2 = _S396 * _S397 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S398 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S393 = true;
        }
        else
        {
            _S393 = xmin_1 >= float(image_width_1);
        }
        if(_S393)
        {
            _S393 = true;
        }
        else
        {
            _S393 = ymax_1 <= 0.0f;
        }
        if(_S393)
        {
            _S393 = true;
        }
        else
        {
            _S393 = ymin_1 >= float(image_height_1);
        }
        if(_S393)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S398.rows[int(0)].x, _S398.rows[int(0)].y, _S398.rows[int(1)].y);
        float3  _S399 = mean_1 - - mul_0(_S395, t_4);
        float3  _S400 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S400;
        float _S401 = _S399.x;
        float _S402 = _S399.y;
        float _S403 = _S399.z;
        float norm_1 = (F32_sqrt((_S401 * _S401 + _S402 * _S402 + _S403 * _S403)));
        float x_17 = _S401 / norm_1;
        float y_4 = _S402 / norm_1;
        float z_1 = _S403 / norm_1;
        float3  _S404 = _S400 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_17) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S404;
        float z2_6 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_17 * x_17 - y_4 * y_4;
        float fS1_1 = 2.0f * x_17 * y_4;
        float3  _S405 = _S404 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_17) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S405;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S405 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_17) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_6, float fy_6, float cx_6, float cy_6, FixedArray<float, 10>  * dist_coeffs_12, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S406 = mean_c_2.z;
        bool _S407;
        if(_S406 < near_plane_2)
        {
            _S407 = true;
        }
        else
        {
            _S407 = _S406 > far_plane_2;
        }
        if(_S407)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S408 = exp_0(scale_4);
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
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S408.x, 0.0f, 0.0f, 0.0f, _S408.y, 0.0f, 0.0f, 0.0f, _S408.z));
        Matrix<float, 3, 3>  _S409 = transpose_0(R_6);
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S409);
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S410 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S410;
        float _S411 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S411;
        float det_blur_3 = _S410 * _S411 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S412 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S407 = true;
        }
        else
        {
            _S407 = xmin_2 >= float(image_width_2);
        }
        if(_S407)
        {
            _S407 = true;
        }
        else
        {
            _S407 = ymax_2 <= 0.0f;
        }
        if(_S407)
        {
            _S407 = true;
        }
        else
        {
            _S407 = ymin_2 >= float(image_height_2);
        }
        if(_S407)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S412.rows[int(0)].x, _S412.rows[int(0)].y, _S412.rows[int(1)].y);
        float3  _S413 = mean_2 - - mul_0(_S409, t_5);
        float3  _S414 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S414;
        float _S415 = _S413.x;
        float _S416 = _S413.y;
        float _S417 = _S413.z;
        float norm_2 = (F32_sqrt((_S415 * _S415 + _S416 * _S416 + _S417 * _S417)));
        float x_19 = _S415 / norm_2;
        float y_5 = _S416 / norm_2;
        float z_2 = _S417 / norm_2;
        float3  _S418 = _S414 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_19) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S418;
        float z2_8 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_19 * x_19 - y_5 * y_5;
        float fS1_2 = 2.0f * x_19 * y_5;
        float3  _S419 = _S418 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_19) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S419;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S419 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_19) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_7, float fy_7, float cx_7, float cy_7, FixedArray<float, 10>  * dist_coeffs_13, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
        float _S420 = mean_c_3.z;
        bool _S421;
        if(_S420 < near_plane_3)
        {
            _S421 = true;
        }
        else
        {
            _S421 = _S420 > far_plane_3;
        }
        if(_S421)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S422 = exp_0(scale_5);
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
        Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S422.x, 0.0f, 0.0f, 0.0f, _S422.y, 0.0f, 0.0f, 0.0f, _S422.z));
        Matrix<float, 3, 3>  _S423 = transpose_0(R_7);
        Matrix<float, 3, 3>  covar_c_3 = mul_4(mul_4(R_7, mul_4(M_5, transpose_0(M_5))), _S423);
        Matrix<float, 2, 2>  covar2d_3;
        persp_proj_3dgs(mean_c_3, covar_c_3, fx_7, fy_7, cx_7, cy_7, dist_coeffs_13, &covar2d_3, mean2d_7);
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S424 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S424;
        float _S425 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S425;
        float det_blur_4 = _S424 * _S425 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
            _S421 = true;
        }
        else
        {
            _S421 = xmin_3 >= float(image_width_3);
        }
        if(_S421)
        {
            _S421 = true;
        }
        else
        {
            _S421 = ymax_3 <= 0.0f;
        }
        if(_S421)
        {
            _S421 = true;
        }
        else
        {
            _S421 = ymin_3 >= float(image_height_3);
        }
        if(_S421)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_5);
        float3  _S426 = mean_3 - - mul_0(_S423, t_6);
        float3  _S427 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S427;
        float _S428 = _S426.x;
        float _S429 = _S426.y;
        float _S430 = _S426.z;
        float norm_3 = (F32_sqrt((_S428 * _S428 + _S429 * _S429 + _S430 * _S430)));
        float x_21 = _S428 / norm_3;
        float y_6 = _S429 / norm_3;
        float z_3 = _S430 / norm_3;
        float3  _S431 = _S427 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_21) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S431;
        float z2_10 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_21 * x_21 - y_6 * y_6;
        float fS1_3 = 2.0f * x_21 * y_6;
        float3  _S432 = _S431 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_21) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S432;
        float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S432 + (make_float3 (-0.59004360437393188f * (x_21 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_21) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_21 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_8, float fy_8, float cx_8, float cy_8, FixedArray<float, 10>  * dist_coeffs_14, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
        float _S433 = length_1(mean_c_4);
        bool _S434;
        if(_S433 < near_plane_4)
        {
            _S434 = true;
        }
        else
        {
            _S434 = _S433 > far_plane_4;
        }
        if(_S434)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S435 = exp_0(scale_6);
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
        Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S435.x, 0.0f, 0.0f, 0.0f, _S435.y, 0.0f, 0.0f, 0.0f, _S435.z));
        Matrix<float, 3, 3>  _S436 = transpose_0(R_8);
        Matrix<float, 3, 3>  covar_c_4 = mul_4(mul_4(R_8, mul_4(M_6, transpose_0(M_6))), _S436);
        Matrix<float, 2, 2>  covar2d_4;
        fisheye_proj_3dgs(mean_c_4, covar_c_4, fx_8, fy_8, cx_8, cy_8, dist_coeffs_14, &covar2d_4, mean2d_8);
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S437 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S437;
        float _S438 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S438;
        float det_blur_5 = _S437 * _S438 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
            _S434 = true;
        }
        else
        {
            _S434 = xmin_4 >= float(image_width_4);
        }
        if(_S434)
        {
            _S434 = true;
        }
        else
        {
            _S434 = ymax_4 <= 0.0f;
        }
        if(_S434)
        {
            _S434 = true;
        }
        else
        {
            _S434 = ymin_4 >= float(image_height_4);
        }
        if(_S434)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_6);
        float3  _S439 = mean_4 - - mul_0(_S436, t_7);
        float3  _S440 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S440;
        float _S441 = _S439.x;
        float _S442 = _S439.y;
        float _S443 = _S439.z;
        float norm_4 = (F32_sqrt((_S441 * _S441 + _S442 * _S442 + _S443 * _S443)));
        float x_23 = _S441 / norm_4;
        float y_7 = _S442 / norm_4;
        float z_4 = _S443 / norm_4;
        float3  _S444 = _S440 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_23) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S444;
        float z2_12 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_23 * x_23 - y_7 * y_7;
        float fS1_4 = 2.0f * x_23 * y_7;
        float3  _S445 = _S444 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_23) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S445;
        float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S445 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_23) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_9, float fy_9, float cx_9, float cy_9, FixedArray<float, 10>  * dist_coeffs_15, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S446 = exp_0(scale_7);
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
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S446.x, 0.0f, 0.0f, 0.0f, _S446.y, 0.0f, 0.0f, 0.0f, _S446.z));
    Matrix<float, 3, 3>  _S447 = transpose_0(R_9);
    Matrix<float, 2, 2>  covar2d_5;
    persp_proj_3dgs(mean_c_5, mul_4(mul_4(R_9, mul_4(M_7, transpose_0(M_7))), _S447), fx_9, fy_9, cx_9, cy_9, dist_coeffs_15, &covar2d_5, mean2d_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S448 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S448;
    float _S449 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S449;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S448 * _S449 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S450 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
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
    *conic_5 = make_float3 (_S450.rows[int(0)].x, _S450.rows[int(0)].y, _S450.rows[int(1)].y);
    float3  _S451 = mean_5 - - mul_0(_S447, t_8);
    float3  _S452 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S452;
    float _S453 = _S451.x;
    float _S454 = _S451.y;
    float _S455 = _S451.z;
    float norm_5 = (F32_sqrt((_S453 * _S453 + _S454 * _S454 + _S455 * _S455)));
    float x_25 = _S453 / norm_5;
    float y_8 = _S454 / norm_5;
    float z_5 = _S455 / norm_5;
    float3  _S456 = _S452 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_25) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S456;
    float z2_14 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_25 * x_25 - y_8 * y_8;
    float fS1_5 = 2.0f * x_25 * y_8;
    float3  _S457 = _S456 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_25) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S457;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S457 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_25) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_10, float fy_10, float cx_10, float cy_10, FixedArray<float, 10>  * dist_coeffs_16, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_10, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_6) + t_9;
    float3  _S458 = exp_0(scale_8);
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
    Matrix<float, 3, 3>  M_8 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S458.x, 0.0f, 0.0f, 0.0f, _S458.y, 0.0f, 0.0f, 0.0f, _S458.z));
    Matrix<float, 3, 3>  _S459 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    fisheye_proj_3dgs(mean_c_6, mul_4(mul_4(R_10, mul_4(M_8, transpose_0(M_8))), _S459), fx_10, fy_10, cx_10, cy_10, dist_coeffs_16, &covar2d_6, mean2d_10);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S460 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S460;
    float _S461 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S461;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S460 * _S461 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S462 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
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
    *conic_6 = make_float3 (_S462.rows[int(0)].x, _S462.rows[int(0)].y, _S462.rows[int(1)].y);
    float3  _S463 = mean_6 - - mul_0(_S459, t_9);
    float3  _S464 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S464;
    float _S465 = _S463.x;
    float _S466 = _S463.y;
    float _S467 = _S463.z;
    float norm_6 = (F32_sqrt((_S465 * _S465 + _S466 * _S466 + _S467 * _S467)));
    float x_27 = _S465 / norm_6;
    float y_9 = _S466 / norm_6;
    float z_6 = _S467 / norm_6;
    float3  _S468 = _S464 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_27) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S468;
    float z2_16 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_27 * x_27 - y_9 * y_9;
    float fS1_6 = 2.0f * x_27 * y_9;
    float3  _S469 = _S468 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_16 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_27) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S469;
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S469 + (make_float3 (-0.59004360437393188f * (x_27 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_16 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_27) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_27 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_11, float fy_11, float cx_11, float cy_11, FixedArray<float, 10>  * dist_coeffs_17, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_11, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_7) + t_10;
    float3  _S470 = exp_0(scale_9);
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
    Matrix<float, 3, 3>  M_9 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S470.x, 0.0f, 0.0f, 0.0f, _S470.y, 0.0f, 0.0f, 0.0f, _S470.z));
    Matrix<float, 3, 3>  _S471 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_11, 0.0f, 0.0f, 0.0f, fy_11, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_7, mul_4(mul_4(R_11, mul_4(M_9, transpose_0(M_9))), _S471)), transpose_1(J_7));
    *mean2d_11 = make_float2 (fx_11 * mean_c_7.x + cx_11, fy_11 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S472 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S472;
    float _S473 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S473;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S472 * _S473 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S474 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
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
    *conic_7 = make_float3 (_S474.rows[int(0)].x, _S474.rows[int(0)].y, _S474.rows[int(1)].y);
    float3  _S475 = mean_7 - - mul_0(_S471, t_10);
    float3  _S476 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S476;
    float _S477 = _S475.x;
    float _S478 = _S475.y;
    float _S479 = _S475.z;
    float norm_7 = (F32_sqrt((_S477 * _S477 + _S478 * _S478 + _S479 * _S479)));
    float x_29 = _S477 / norm_7;
    float y_10 = _S478 / norm_7;
    float z_7 = _S479 / norm_7;
    float3  _S480 = _S476 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_29) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S480;
    float z2_18 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_29 * x_29 - y_10 * y_10;
    float fS1_7 = 2.0f * x_29 * y_10;
    float3  _S481 = _S480 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_18 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_29) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S481;
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S481 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_18 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_29) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_12, float fy_12, float cx_12, float cy_12, FixedArray<float, 10>  * dist_coeffs_18, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_12, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_8) + t_11;
    float3  _S482 = exp_0(scale_10);
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
    Matrix<float, 3, 3>  M_10 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S482.x, 0.0f, 0.0f, 0.0f, _S482.y, 0.0f, 0.0f, 0.0f, _S482.z));
    Matrix<float, 3, 3>  _S483 = transpose_0(R_12);
    Matrix<float, 2, 2>  covar2d_8;
    persp_proj_3dgs(mean_c_8, mul_4(mul_4(R_12, mul_4(M_10, transpose_0(M_10))), _S483), fx_12, fy_12, cx_12, cy_12, dist_coeffs_18, &covar2d_8, mean2d_12);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S484 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S484;
    float _S485 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S485;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S484 * _S485 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
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
    float3  _S486 = mean_8 - - mul_0(_S483, t_11);
    float3  _S487 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S487;
    float _S488 = _S486.x;
    float _S489 = _S486.y;
    float _S490 = _S486.z;
    float norm_8 = (F32_sqrt((_S488 * _S488 + _S489 * _S489 + _S490 * _S490)));
    float x_31 = _S488 / norm_8;
    float y_11 = _S489 / norm_8;
    float z_8 = _S490 / norm_8;
    float3  _S491 = _S487 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_31) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S491;
    float z2_20 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_31 * x_31 - y_11 * y_11;
    float fS1_8 = 2.0f * x_31 * y_11;
    float3  _S492 = _S491 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_20 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_31) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S492;
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S492 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_20 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_31) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_13, float fy_13, float cx_13, float cy_13, FixedArray<float, 10>  * dist_coeffs_19, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_13, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_9) + t_12;
    float3  _S493 = exp_0(scale_11);
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
    Matrix<float, 3, 3>  M_11 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))), makeMatrix<float, 3, 3> (_S493.x, 0.0f, 0.0f, 0.0f, _S493.y, 0.0f, 0.0f, 0.0f, _S493.z));
    Matrix<float, 3, 3>  _S494 = transpose_0(R_13);
    Matrix<float, 2, 2>  covar2d_9;
    fisheye_proj_3dgs(mean_c_9, mul_4(mul_4(R_13, mul_4(M_11, transpose_0(M_11))), _S494), fx_13, fy_13, cx_13, cy_13, dist_coeffs_19, &covar2d_9, mean2d_13);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S495 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S495;
    float _S496 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S496;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S495 * _S496 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
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
    float3  _S497 = mean_9 - - mul_0(_S494, t_12);
    float3  _S498 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S498;
    float _S499 = _S497.x;
    float _S500 = _S497.y;
    float _S501 = _S497.z;
    float norm_9 = (F32_sqrt((_S499 * _S499 + _S500 * _S500 + _S501 * _S501)));
    float x_33 = _S499 / norm_9;
    float y_12 = _S500 / norm_9;
    float z_9 = _S501 / norm_9;
    float3  _S502 = _S498 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_33) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S502;
    float z2_22 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_33 * x_33 - y_12 * y_12;
    float fS1_9 = 2.0f * x_33 * y_12;
    float3  _S503 = _S502 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_33) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S503;
    float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S503 + (make_float3 (-0.59004360437393188f * (x_33 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_33) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_33 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

struct s_bwd_prop_persp_proj_3dgs_Intermediates_0
{
    bool _S504;
    bool _S505;
    bool _S506;
};

struct s_bwd_prop_projection_3dgs_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S507;
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S508;
};

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S509, float3  _S510)
{
    return mul_0(_S509, _S510);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S511)
{
    return exp_0(_S511);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S512, Matrix<float, 3, 3>  _S513)
{
    return mul_4(_S512, _S513);
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S514, Matrix<float, 3, 3>  _S515)
{
    return mul_5(_S514, _S515);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S516, Matrix<float, 3, 2>  _S517)
{
    return mul_6(_S516, _S517);
}

inline __device__ void s_primal_ctx_persp_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, FixedArray<float, 10>  * dpdist_coeffs_0, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_persp_proj_3dgs_Intermediates_0 * _s_diff_ctx_2)
{
    _s_diff_ctx_2->_S504 = false;
    _s_diff_ctx_2->_S505 = false;
    _s_diff_ctx_2->_S506 = false;
    _s_diff_ctx_2->_S504 = false;
    _s_diff_ctx_2->_S505 = false;
    _s_diff_ctx_2->_S506 = false;
    float2  _S518 = make_float2 (0.0f);
    float2  _S519 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S520 = dpmean3d_0.z;
    float2  uv_13 = _S519 / make_float2 (_S520);
    bool _S521 = _S520 < 0.0f;
    bool _S522;
    if(_S521)
    {
        _S522 = true;
    }
    else
    {
        bool _S523 = is_valid_distortion(uv_13, dpdist_coeffs_0);
        _s_diff_ctx_2->_S504 = _S523;
        _S522 = !_S523;
    }
    float2  _S524;
    if(_S522)
    {
        _S524 = make_float2 (-10000.0f);
    }
    bool _S525 = !_S522;
    if(_S525)
    {
        float u_10 = uv_13.x;
        float v_10 = uv_13.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float2  _S526 = uv_13 * make_float2 (1.0f + r2_10 * ((*dpdist_coeffs_0)[int(0)] + r2_10 * ((*dpdist_coeffs_0)[int(1)] + r2_10 * ((*dpdist_coeffs_0)[int(2)] + r2_10 * (*dpdist_coeffs_0)[int(3)])))) + make_float2 (2.0f * (*dpdist_coeffs_0)[int(4)] * u_10 * v_10 + (*dpdist_coeffs_0)[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + (*dpdist_coeffs_0)[int(6)] * r2_10, 2.0f * (*dpdist_coeffs_0)[int(5)] * u_10 * v_10 + (*dpdist_coeffs_0)[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + (*dpdist_coeffs_0)[int(7)] * r2_10);
        float2  _S527 = _S526 + make_float2 ((*dpdist_coeffs_0)[int(8)] * _S526.x + (*dpdist_coeffs_0)[int(9)] * _S526.y, 0.0f);
        _S524 = make_float2 (dpfx_0 * _S527.x + dpcx_0, dpfy_0 * _S527.y + dpcy_0);
    }
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (0.0f);
    float2  uv_14 = _S519 / make_float2 (_S520);
    s_bwd_prop_persp_proj_Intermediates_0 _S528;
    if(_S521)
    {
        (&_S528)->_S174 = false;
    }
    else
    {
        bool _S529 = is_valid_distortion(uv_14, dpdist_coeffs_0);
        _s_diff_ctx_2->_S505 = _S529;
        (&_S528)->_S174 = _S529;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S530 = _S528;
    float2  _S531 = make_float2 (_S520);
    float2  uv_15 = _S519 / make_float2 (_S520);
    float2  _S532 = make_float2 (_S520 * _S520);
    if(_S521)
    {
        _S522 = true;
    }
    else
    {
        _S522 = !_S530._S174;
    }
    bool _S533 = !_S522;
    float2  _S534;
    float _S535;
    float _S536;
    float _S537;
    float _S538;
    float _S539;
    float _S540;
    float _S541;
    float _S542;
    float _S543;
    float _S544;
    float _S545;
    float _S546;
    float _S547;
    float _S548;
    float _S549;
    float _S550;
    float _S551;
    float _S552;
    float _S553;
    float _S554;
    if(_S533)
    {
        float u_11 = uv_15.x;
        float v_11 = uv_15.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float _S555 = (*dpdist_coeffs_0)[int(2)] + r2_11 * (*dpdist_coeffs_0)[int(3)];
        float _S556 = (*dpdist_coeffs_0)[int(1)] + r2_11 * _S555;
        float _S557 = (*dpdist_coeffs_0)[int(0)] + r2_11 * _S556;
        float2  _S558 = make_float2 (1.0f + r2_11 * _S557);
        float _S559 = 2.0f * (*dpdist_coeffs_0)[int(4)];
        float _S560 = _S559 * u_11;
        float _S561 = 2.0f * u_11;
        float _S562 = 2.0f * (*dpdist_coeffs_0)[int(5)];
        float _S563 = _S562 * u_11;
        float _S564 = 2.0f * v_11;
        _S535 = dpfx_0;
        _S536 = (*dpdist_coeffs_0)[int(9)];
        _S537 = (*dpdist_coeffs_0)[int(8)];
        _S534 = _S558;
        _S538 = (*dpdist_coeffs_0)[int(7)];
        _S539 = (*dpdist_coeffs_0)[int(4)];
        _S540 = _S564;
        _S541 = v_11;
        _S542 = _S563;
        _S543 = _S562;
        _S544 = (*dpdist_coeffs_0)[int(6)];
        _S545 = (*dpdist_coeffs_0)[int(5)];
        _S546 = _S561;
        _S547 = u_11;
        _S548 = _S560;
        _S549 = _S559;
        _S550 = r2_11;
        _S551 = _S557;
        _S552 = _S556;
        _S553 = _S555;
        _S554 = (*dpdist_coeffs_0)[int(3)];
    }
    else
    {
        _S535 = 0.0f;
        _S536 = 0.0f;
        _S537 = 0.0f;
        _S534 = _S518;
        _S538 = 0.0f;
        _S539 = 0.0f;
        _S540 = 0.0f;
        _S541 = 0.0f;
        _S542 = 0.0f;
        _S543 = 0.0f;
        _S544 = 0.0f;
        _S545 = 0.0f;
        _S546 = 0.0f;
        _S547 = 0.0f;
        _S548 = 0.0f;
        _S549 = 0.0f;
        _S550 = 0.0f;
        _S551 = 0.0f;
        _S552 = 0.0f;
        _S553 = 0.0f;
        _S554 = 0.0f;
    }
    if(_S533)
    {
        float2  _S565 = make_float2 (_S535, 0.0f) + make_float2 (_S537 * _S535, _S536 * _S535);
        float2  _S566 = uv_15 * _S565;
        float _S567 = _S539 * _S565.y;
        float _S568 = _S545 * _S565.x;
        float _S569 = _S566.x + _S566.y;
        float _S570 = _S550 * _S569;
        float _S571 = _S550 * _S570;
        float _S572 = _S538 * _S565.y + _S567 + _S544 * _S565.x + _S568 + _S551 * _S569 + _S552 * _S570 + _S553 * _S571 + _S554 * (_S550 * _S571);
        float _S573 = _S541 * _S572;
        float _S574 = _S547 * _S572;
        _S534 = _S534 * _S565 + make_float2 (_S543 * (_S541 * _S565.y) + _S546 * _S568 + 2.0f * (_S547 * _S568) + _S549 * (_S541 * _S565.x) + _S574 + _S574, _S540 * _S567 + 2.0f * (_S541 * _S567) + _S542 * _S565.y + _S548 * _S565.x + _S573 + _S573);
    }
    else
    {
        _S534 = _S518;
    }
    float2  _S575 = _S534 / _S532;
    float2  _S576 = _S519 * - _S575;
    float2  _S577 = _S531 * _S575;
    float3  _S578 = make_float3 (_S577.x, _S577.y, _S576.x + _S576.y);
    Matrix<float, 2, 3>  _S579 = J_8;
    _S579[int(0)] = _S578;
    float2  uv_16 = _S519 / make_float2 (_S520);
    if(_S521)
    {
        (&_S528)->_S174 = false;
    }
    else
    {
        bool _S580 = is_valid_distortion(uv_16, dpdist_coeffs_0);
        _s_diff_ctx_2->_S506 = _S580;
        (&_S528)->_S174 = _S580;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S581 = _S528;
    float2  uv_17 = _S519 / make_float2 (_S520);
    if(_S521)
    {
        _S522 = true;
    }
    else
    {
        _S522 = !_S581._S174;
    }
    bool _S582 = !_S522;
    if(_S582)
    {
        float u_12 = uv_17.x;
        float v_12 = uv_17.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float _S583 = (*dpdist_coeffs_0)[int(2)] + r2_12 * (*dpdist_coeffs_0)[int(3)];
        float _S584 = (*dpdist_coeffs_0)[int(1)] + r2_12 * _S583;
        float _S585 = (*dpdist_coeffs_0)[int(0)] + r2_12 * _S584;
        float2  _S586 = make_float2 (1.0f + r2_12 * _S585);
        float _S587 = 2.0f * (*dpdist_coeffs_0)[int(5)];
        float _S588 = _S587 * u_12;
        float _S589 = 2.0f * v_12;
        _S535 = dpfy_0;
        _S534 = _S586;
        _S536 = (*dpdist_coeffs_0)[int(7)];
        _S537 = (*dpdist_coeffs_0)[int(4)];
        _S538 = _S589;
        _S539 = v_12;
        _S540 = _S588;
        _S541 = _S587;
        _S542 = u_12;
        _S543 = r2_12;
        _S544 = _S585;
        _S545 = _S584;
        _S546 = _S583;
        _S547 = (*dpdist_coeffs_0)[int(3)];
    }
    else
    {
        _S535 = 0.0f;
        _S534 = _S518;
        _S536 = 0.0f;
        _S537 = 0.0f;
        _S538 = 0.0f;
        _S539 = 0.0f;
        _S540 = 0.0f;
        _S541 = 0.0f;
        _S542 = 0.0f;
        _S543 = 0.0f;
        _S544 = 0.0f;
        _S545 = 0.0f;
        _S546 = 0.0f;
        _S547 = 0.0f;
    }
    if(_S582)
    {
        float2  _S590 = make_float2 (0.0f, _S535);
        float2  _S591 = uv_17 * _S590;
        float _S592 = _S537 * _S535;
        float _S593 = _S591.x + _S591.y;
        float _S594 = _S543 * _S593;
        float _S595 = _S543 * _S594;
        float _S596 = _S536 * _S535 + _S592 + _S544 * _S593 + _S545 * _S594 + _S546 * _S595 + _S547 * (_S543 * _S595);
        float _S597 = _S539 * _S596;
        float _S598 = _S542 * _S596;
        _S534 = _S534 * _S590 + make_float2 (_S541 * (_S539 * _S535) + _S598 + _S598, _S538 * _S592 + 2.0f * (_S539 * _S592) + _S540 * _S535 + _S597 + _S597);
    }
    else
    {
        _S534 = _S518;
    }
    float2  _S599 = _S534 / _S532;
    float2  _S600 = _S519 * - _S599;
    float2  _S601 = _S531 * _S599;
    _S579[int(1)] = make_float3 (_S601.x, _S601.y, _S600.x + _S600.y);
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S579, dpcov3d_0), transpose_1(_S579));
    *dpmean2d_0 = _S524;
    return;
}

inline __device__ float s_primal_ctx_max_0(float _S602, float _S603)
{
    return (F32_max((_S602), (_S603)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S604)
{
    return (F32_sqrt((_S604)));
}

inline __device__ float s_primal_ctx_exp_1(float _S605)
{
    return (F32_exp((_S605)));
}

inline __device__ float s_primal_ctx_log_0(float _S606)
{
    return (F32_log((_S606)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S607, float3  _S608)
{
    return dot_0(_S607, _S608);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S609, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S610, float3  _S611)
{
    _d_max_vector_0(_S609, _S610, _S611);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S612, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S613, float3  _S614)
{
    _d_mul_0(_S612, _S613, _S614);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S615, float _S616)
{
    _d_log_0(_S615, _S616);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S617, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S618, float _S619)
{
    _d_dot_0(_S617, _S618, _S619);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S620, DiffPair_float_0 * _S621, float _S622)
{
    _d_min_0(_S620, _S621, _S622);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S623, float _S624)
{
    _d_exp_0(_S623, _S624);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S625, DiffPair_float_0 * _S626, float _S627)
{
    _d_max_0(_S625, _S626, _S627);
    return;
}

struct DiffPair_arrayx3Cfloatx2C10x3E_0
{
    FixedArray<float, 10>  primal_0;
    FixedArray<float, 10>  differential_0;
};

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S628, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S629, Matrix<float, 2, 2>  _S630)
{
    mul_3(_S628, _S629, _S630);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S631, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S632, Matrix<float, 2, 3>  _S633)
{
    mul_2(_S631, _S632, _S633);
    return;
}

inline __device__ void s_bwd_prop_persp_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_1, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_persp_proj_3dgs_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S634 = *dpcov3d_1;
    DiffPair_float_0 _S635 = *dpfx_1;
    DiffPair_float_0 _S636 = *dpfy_1;
    FixedArray<float, 10>  _S637 = dpdist_coeffs_1->primal_0;
    float2  _S638 = make_float2 (0.0f);
    float2  _S639 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S640 = (*dpmean3d_1).primal_0.z;
    float2  _S641 = make_float2 (_S640);
    float2  uv_18 = _S639 / make_float2 (_S640);
    float2  _S642 = make_float2 (_S640 * _S640);
    bool _S643 = _S640 < 0.0f;
    bool _S644;
    if(_S643)
    {
        _S644 = true;
    }
    else
    {
        _S644 = !_s_diff_ctx_3->_S504;
    }
    bool _S645 = !_S644;
    float _S646;
    float _S647;
    float _S648;
    float _S649;
    float _S650;
    float _S651;
    float _S652;
    float _S653;
    float _S654;
    float _S655;
    float _S656;
    float _S657;
    float _S658;
    float _S659;
    float _S660;
    float _S661;
    float _S662;
    float _S663;
    float _S664;
    float _S665;
    float _S666;
    float _S667;
    float _S668;
    float _S669;
    float _S670;
    float2  _S671;
    if(_S645)
    {
        float u_13 = uv_18.x;
        float v_13 = uv_18.y;
        float r2_13 = u_13 * u_13 + v_13 * v_13;
        float _S672 = _S637[int(2)] + r2_13 * _S637[int(3)];
        float _S673 = _S637[int(1)] + r2_13 * _S672;
        float _S674 = _S637[int(0)] + r2_13 * _S673;
        float radial_1 = 1.0f + r2_13 * _S674;
        float2  _S675 = make_float2 (radial_1);
        float _S676 = 2.0f * _S637[int(4)];
        float _S677 = _S676 * u_13;
        float _S678 = 2.0f * u_13;
        float _S679 = r2_13 + _S678 * u_13;
        float _S680 = 2.0f * _S637[int(5)];
        float _S681 = _S680 * u_13;
        float _S682 = 2.0f * v_13;
        float _S683 = r2_13 + _S682 * v_13;
        float2  _S684 = uv_18 * make_float2 (radial_1) + make_float2 (_S677 * v_13 + _S637[int(5)] * _S679 + _S637[int(6)] * r2_13, _S681 * v_13 + _S637[int(4)] * _S683 + _S637[int(7)] * r2_13);
        float _S685 = _S684.x;
        float _S686 = _S684.y;
        float2  _S687 = _S684 + make_float2 (_S637[int(8)] * _S685 + _S637[int(9)] * _S686, 0.0f);
        float _S688 = _S687.x;
        _S646 = _S687.y;
        _S647 = _S688;
        _S648 = _S637[int(9)];
        _S649 = _S686;
        _S650 = _S637[int(8)];
        _S651 = _S685;
        _S671 = _S675;
        _S652 = _S637[int(7)];
        _S653 = r2_13;
        _S654 = _S637[int(4)];
        _S655 = _S683;
        _S656 = _S682;
        _S657 = v_13;
        _S658 = _S681;
        _S659 = _S680;
        _S660 = u_13;
        _S661 = _S637[int(6)];
        _S662 = _S637[int(5)];
        _S663 = _S679;
        _S664 = _S678;
        _S665 = _S677;
        _S666 = _S676;
        _S667 = _S674;
        _S668 = _S673;
        _S669 = _S672;
        _S670 = _S637[int(3)];
    }
    else
    {
        _S646 = 0.0f;
        _S647 = 0.0f;
        _S648 = 0.0f;
        _S649 = 0.0f;
        _S650 = 0.0f;
        _S651 = 0.0f;
        _S671 = _S638;
        _S652 = 0.0f;
        _S653 = 0.0f;
        _S654 = 0.0f;
        _S655 = 0.0f;
        _S656 = 0.0f;
        _S657 = 0.0f;
        _S658 = 0.0f;
        _S659 = 0.0f;
        _S660 = 0.0f;
        _S661 = 0.0f;
        _S662 = 0.0f;
        _S663 = 0.0f;
        _S664 = 0.0f;
        _S665 = 0.0f;
        _S666 = 0.0f;
        _S667 = 0.0f;
        _S668 = 0.0f;
        _S669 = 0.0f;
        _S670 = 0.0f;
    }
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    s_bwd_prop_persp_proj_Intermediates_0 _S689;
    if(_S643)
    {
        (&_S689)->_S174 = false;
    }
    else
    {
        (&_S689)->_S174 = _s_diff_ctx_3->_S505;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S690 = _S689;
    float2  uv_19 = _S639 / make_float2 (_S640);
    if(_S643)
    {
        _S644 = true;
    }
    else
    {
        _S644 = !_S690._S174;
    }
    bool _S691 = !_S644;
    float _S692;
    float _S693;
    float _S694;
    float _S695;
    float _S696;
    float _S697;
    float _S698;
    float _S699;
    float _S700;
    float _S701;
    float _S702;
    float _S703;
    float _S704;
    float _S705;
    float _S706;
    float _S707;
    float _S708;
    float _S709;
    float _S710;
    float _S711;
    float2  _S712;
    if(_S691)
    {
        float u_14 = uv_19.x;
        float v_14 = uv_19.y;
        float r2_14 = u_14 * u_14 + v_14 * v_14;
        float _S713 = _S637[int(2)] + r2_14 * _S637[int(3)];
        float _S714 = _S637[int(1)] + r2_14 * _S713;
        float _S715 = _S637[int(0)] + r2_14 * _S714;
        float2  _S716 = make_float2 (1.0f + r2_14 * _S715);
        float _S717 = 2.0f * _S637[int(4)];
        float _S718 = _S717 * u_14;
        float _S719 = 2.0f * u_14;
        float _S720 = 2.0f * _S637[int(5)];
        float _S721 = _S720 * u_14;
        float _S722 = 2.0f * v_14;
        _S692 = _S635.primal_0;
        _S693 = _S637[int(9)];
        _S694 = _S637[int(8)];
        _S712 = _S716;
        _S695 = _S637[int(7)];
        _S696 = _S637[int(4)];
        _S697 = _S722;
        _S698 = v_14;
        _S699 = _S721;
        _S700 = _S720;
        _S701 = _S637[int(6)];
        _S702 = _S637[int(5)];
        _S703 = _S719;
        _S704 = u_14;
        _S705 = _S718;
        _S706 = _S717;
        _S707 = r2_14;
        _S708 = _S715;
        _S709 = _S714;
        _S710 = _S713;
        _S711 = _S637[int(3)];
    }
    else
    {
        _S692 = 0.0f;
        _S693 = 0.0f;
        _S694 = 0.0f;
        _S712 = _S638;
        _S695 = 0.0f;
        _S696 = 0.0f;
        _S697 = 0.0f;
        _S698 = 0.0f;
        _S699 = 0.0f;
        _S700 = 0.0f;
        _S701 = 0.0f;
        _S702 = 0.0f;
        _S703 = 0.0f;
        _S704 = 0.0f;
        _S705 = 0.0f;
        _S706 = 0.0f;
        _S707 = 0.0f;
        _S708 = 0.0f;
        _S709 = 0.0f;
        _S710 = 0.0f;
        _S711 = 0.0f;
    }
    float _S723;
    float _S724;
    float _S725;
    float _S726;
    float _S727;
    float _S728;
    float _S729;
    float _S730;
    float _S731;
    float _S732;
    float _S733;
    float2  _S734;
    float2  _S735;
    if(_S691)
    {
        float2  _S736 = make_float2 (_S692, 0.0f) + make_float2 (_S694 * _S692, _S693 * _S692);
        float2  _S737 = uv_19 * _S736;
        float _S738 = _S696 * _S736.y;
        float _S739 = _S698 * _S736.y;
        float _S740 = _S702 * _S736.x;
        float _S741 = _S698 * _S736.x;
        float _S742 = _S737.x + _S737.y;
        float _S743 = _S707 * _S742;
        float _S744 = _S707 * _S743;
        float _S745 = _S707 * _S744;
        float _S746 = _S695 * _S736.y + _S738 + _S701 * _S736.x + _S740 + _S708 * _S742 + _S709 * _S743 + _S710 * _S744 + _S711 * _S745;
        float _S747 = _S698 * _S746;
        float _S748 = _S704 * _S746;
        _S734 = _S712 * _S736 + make_float2 (_S700 * _S739 + _S703 * _S740 + 2.0f * (_S704 * _S740) + _S706 * _S741 + _S748 + _S748, _S697 * _S738 + 2.0f * (_S698 * _S738) + _S699 * _S736.y + _S705 * _S736.x + _S747 + _S747);
        _S723 = _S746;
        _S724 = _S745;
        _S725 = _S744;
        _S726 = _S743;
        _S727 = _S742;
        _S728 = _S741;
        _S729 = _S736.x;
        _S730 = _S740;
        _S731 = _S739;
        _S732 = _S736.y;
        _S733 = _S738;
        _S735 = _S736;
    }
    else
    {
        _S734 = _S638;
        _S723 = 0.0f;
        _S724 = 0.0f;
        _S725 = 0.0f;
        _S726 = 0.0f;
        _S727 = 0.0f;
        _S728 = 0.0f;
        _S729 = 0.0f;
        _S730 = 0.0f;
        _S731 = 0.0f;
        _S732 = 0.0f;
        _S733 = 0.0f;
        _S735 = _S638;
    }
    float2  _S749 = _S734 / _S642;
    float2  _S750 = _S642 * _S642;
    float2  _S751 = - _S749;
    float2  _S752 = _S639 * _S751;
    float2  _S753 = _S641 * _S749;
    float3  _S754 = make_float3 (_S753.x, _S753.y, _S752.x + _S752.y);
    Matrix<float, 2, 3>  _S755 = J_9;
    _S755[int(0)] = _S754;
    if(_S643)
    {
        (&_S689)->_S174 = false;
    }
    else
    {
        (&_S689)->_S174 = _s_diff_ctx_3->_S506;
    }
    s_bwd_prop_persp_proj_Intermediates_0 _S756 = _S689;
    float2  uv_20 = _S639 / make_float2 (_S640);
    if(_S643)
    {
        _S644 = true;
    }
    else
    {
        _S644 = !_S756._S174;
    }
    bool _S757 = !_S644;
    float _S758;
    float _S759;
    float _S760;
    float _S761;
    float _S762;
    float _S763;
    float _S764;
    float _S765;
    float _S766;
    float _S767;
    float _S768;
    float _S769;
    float _S770;
    float _S771;
    float _S772;
    float _S773;
    float _S774;
    float _S775;
    float2  _S776;
    if(_S757)
    {
        float u_15 = uv_20.x;
        float v_15 = uv_20.y;
        float r2_15 = u_15 * u_15 + v_15 * v_15;
        float _S777 = _S637[int(2)] + r2_15 * _S637[int(3)];
        float _S778 = _S637[int(1)] + r2_15 * _S777;
        float _S779 = _S637[int(0)] + r2_15 * _S778;
        float2  _S780 = make_float2 (1.0f + r2_15 * _S779);
        float _S781 = 2.0f * _S637[int(4)];
        float _S782 = _S781 * u_15;
        float _S783 = 2.0f * u_15;
        float _S784 = 2.0f * _S637[int(5)];
        float _S785 = _S784 * u_15;
        float _S786 = 2.0f * v_15;
        _S758 = _S636.primal_0;
        _S776 = _S780;
        _S759 = _S637[int(7)];
        _S760 = _S637[int(4)];
        _S761 = _S786;
        _S762 = v_15;
        _S763 = _S785;
        _S764 = _S784;
        _S765 = _S637[int(6)];
        _S766 = _S637[int(5)];
        _S767 = _S783;
        _S768 = u_15;
        _S769 = _S782;
        _S770 = _S781;
        _S771 = r2_15;
        _S772 = _S779;
        _S773 = _S778;
        _S774 = _S777;
        _S775 = _S637[int(3)];
    }
    else
    {
        _S758 = 0.0f;
        _S776 = _S638;
        _S759 = 0.0f;
        _S760 = 0.0f;
        _S761 = 0.0f;
        _S762 = 0.0f;
        _S763 = 0.0f;
        _S764 = 0.0f;
        _S765 = 0.0f;
        _S766 = 0.0f;
        _S767 = 0.0f;
        _S768 = 0.0f;
        _S769 = 0.0f;
        _S770 = 0.0f;
        _S771 = 0.0f;
        _S772 = 0.0f;
        _S773 = 0.0f;
        _S774 = 0.0f;
        _S775 = 0.0f;
    }
    float _S787;
    float _S788;
    float _S789;
    float _S790;
    float _S791;
    float _S792;
    float _S793;
    float2  _S794;
    float2  _S795;
    if(_S757)
    {
        float2  _S796 = make_float2 (0.0f, _S758);
        float2  _S797 = uv_20 * _S796;
        float _S798 = _S760 * _S758;
        float _S799 = _S762 * _S758;
        float _S800 = _S797.x + _S797.y;
        float _S801 = _S771 * _S800;
        float _S802 = _S771 * _S801;
        float _S803 = _S771 * _S802;
        float _S804 = _S759 * _S758 + _S798 + _S772 * _S800 + _S773 * _S801 + _S774 * _S802 + _S775 * _S803;
        float _S805 = _S762 * _S804;
        float _S806 = _S768 * _S804;
        float _S807 = _S758;
        _S794 = _S776 * _S796 + make_float2 (_S764 * _S799 + _S806 + _S806, _S761 * _S798 + 2.0f * (_S762 * _S798) + _S763 * _S758 + _S805 + _S805);
        _S758 = _S804;
        _S787 = _S803;
        _S788 = _S802;
        _S789 = _S801;
        _S790 = _S800;
        _S791 = _S799;
        _S792 = _S807;
        _S793 = _S798;
        _S795 = _S796;
    }
    else
    {
        _S794 = _S638;
        _S758 = 0.0f;
        _S787 = 0.0f;
        _S788 = 0.0f;
        _S789 = 0.0f;
        _S790 = 0.0f;
        _S791 = 0.0f;
        _S792 = 0.0f;
        _S793 = 0.0f;
        _S795 = _S638;
    }
    float2  _S808 = _S794 / _S642;
    float2  _S809 = - _S808;
    float2  _S810 = _S639 * _S809;
    float2  _S811 = _S641 * _S808;
    _S755[int(1)] = make_float3 (_S811.x, _S811.y, _S810.x + _S810.y);
    Matrix<float, 3, 2>  _S812 = transpose_1(_S755);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S813;
    (&_S813)->primal_0 = s_primal_ctx_mul_3(_S755, _S634.primal_0);
    (&_S813)->differential_0 = J_9;
    Matrix<float, 3, 2>  _S814 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S815;
    (&_S815)->primal_0 = _S812;
    (&_S815)->differential_0 = _S814;
    s_bwd_prop_mul_2(&_S813, &_S815, dpcov2d_1);
    Matrix<float, 2, 3>  _S816 = transpose_2(_S815.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S817;
    (&_S817)->primal_0 = _S755;
    (&_S817)->differential_0 = J_9;
    Matrix<float, 3, 3>  _S818 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S819;
    (&_S819)->primal_0 = _S634.primal_0;
    (&_S819)->differential_0 = _S818;
    s_bwd_prop_mul_3(&_S817, &_S819, _S813.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S820 = _S819;
    Matrix<float, 2, 3>  _S821 = _S816 + _S817.differential_0;
    float2  _S822 = _S638;
    *&((&_S822)->y) = _S821.rows[int(1)].y;
    *&((&_S822)->x) = _S821.rows[int(1)].x;
    float2  _S823 = _S641 * _S822;
    float2  _S824 = _S808 * _S822;
    float2  _S825 = _S638;
    *&((&_S825)->y) = _S821.rows[int(1)].z;
    *&((&_S825)->x) = _S821.rows[int(1)].z;
    float2  _S826 = _S809 * _S825;
    float2  _S827 = (_S823 + - (_S639 * _S825)) / _S750;
    float2  _S828 = _S794 * - _S827;
    float2  _S829 = _S642 * _S827;
    if(_S757)
    {
        float _S830 = _S829.x;
        float _S831 = _S830 + _S830;
        float _S832 = _S758 * _S831;
        float _S833 = _S829.y + _S829.y;
        float _S834 = _S758 * _S833;
        float _S835 = _S768 * _S831 + _S762 * _S833;
        float _S836 = _S775 * _S835;
        float _S837 = _S787 * _S835;
        float _S838 = _S788 * _S835;
        float _S839 = _S788 * _S836;
        float _S840 = _S789 * _S835;
        float _S841 = _S774 * _S835 + _S771 * _S836;
        float _S842 = _S789 * _S841;
        float _S843 = _S790 * _S835;
        float _S844 = _S773 * _S835 + _S771 * _S841;
        float _S845 = _S790 * _S844;
        float _S846 = _S772 * _S835 + _S771 * _S844;
        float _S847 = _S762 * (_S770 * _S829.x);
        float _S848 = _S769 * _S829.y;
        float _S849 = _S766 * (_S835 + _S768 * (2.0f * _S829.x) + _S767 * _S829.x);
        float _S850 = _S765 * _S835;
        float _S851 = _S764 * _S829.x;
        float _S852 = _S791 * _S829.x;
        float _S853 = _S762 * _S851;
        float _S854 = _S792 * _S851;
        float _S855 = _S763 * _S829.y;
        float _S856 = _S792 * _S829.y;
        float _S857 = 2.0f * _S829.y;
        float _S858 = _S793 * _S857;
        float _S859 = _S793 * _S829.y;
        float _S860 = _S835 + _S762 * _S857 + _S761 * _S829.y;
        float _S861 = _S760 * _S860;
        float _S862 = _S792 * _S860;
        float _S863 = _S759 * _S835;
        float _S864 = _S792 * _S835;
        float2  _S865 = _S776 * _S829;
        float2  _S866 = _S795 * _S829;
        float2  _S867 = _S638;
        *&((&_S867)->y) = _S846;
        *&((&_S867)->x) = _S846;
        float2  _S868 = _S795 * _S867;
        float _S869 = _S853 + _S855 + _S861 + _S863;
        float _S870 = _S847 + _S848 + _S849 + _S850;
        float2  _S871 = _S865 + uv_20 * _S867;
        float2  _S872 = _S638;
        *&((&_S872)->y) = _S869;
        *&((&_S872)->x) = _S870;
        float _S873 = _S834 + _S854 + _S858;
        float _S874 = _S839 + _S842 + _S845;
        _S758 = (_S871 + _S872).y;
        _S776 = _S866;
        _S759 = _S864;
        _S760 = _S862;
        _S761 = _S859;
        _S763 = _S873;
        _S765 = _S856;
        _S766 = _S852;
        _S767 = _S832;
        _S769 = _S874;
        _S770 = _S843;
        _S787 = _S840;
        _S788 = _S838;
        _S789 = _S837;
        _S794 = _S868;
    }
    else
    {
        _S758 = 0.0f;
        _S776 = _S638;
        _S759 = 0.0f;
        _S760 = 0.0f;
        _S761 = 0.0f;
        _S763 = 0.0f;
        _S765 = 0.0f;
        _S766 = 0.0f;
        _S767 = 0.0f;
        _S769 = 0.0f;
        _S770 = 0.0f;
        _S787 = 0.0f;
        _S788 = 0.0f;
        _S789 = 0.0f;
        _S794 = _S638;
    }
    FixedArray<float, 10>  _S875;
    if(_S757)
    {
        float _S876 = 2.0f * (_S766 + _S768 * _S765);
        float _S877 = _S776.x + _S776.y;
        float _S878 = _S770 + _S771 * _S877;
        float _S879 = _S787 + _S771 * _S878;
        float _S880 = _S788 + _S771 * _S879;
        float _S881 = _S769 + _S772 * _S877 + _S773 * _S878 + _S774 * _S879 + _S775 * _S880;
        float _S882 = _S762 * _S881;
        float _S883 = _S768 * _S881;
        float _S884 = _S789 + _S771 * _S880;
        float2  _S885 = _S794 + make_float2 (_S767 + _S764 * _S765 + _S883 + _S883, _S763 + 2.0f * _S761 + _S882 + _S882);
        FixedArray<float, 10>  _S886;
        _S886[int(0)] = 0.0f;
        _S886[int(1)] = 0.0f;
        _S886[int(2)] = 0.0f;
        _S886[int(3)] = 0.0f;
        _S886[int(4)] = 0.0f;
        _S886[int(5)] = 0.0f;
        _S886[int(6)] = 0.0f;
        _S886[int(7)] = 0.0f;
        _S886[int(8)] = 0.0f;
        _S886[int(9)] = 0.0f;
        _S886[int(9)] = 0.0f;
        _S886[int(8)] = 0.0f;
        _S886[int(7)] = _S759;
        _S886[int(6)] = 0.0f;
        _S886[int(5)] = _S876;
        _S886[int(4)] = _S760;
        _S886[int(3)] = _S884;
        _S886[int(2)] = _S880;
        _S886[int(1)] = _S879;
        _S886[int(0)] = _S878;
        _S776 = _S885;
        _S875[int(0)] = _S886[int(0)];
        _S875[int(1)] = _S886[int(1)];
        _S875[int(2)] = _S886[int(2)];
        _S875[int(3)] = _S886[int(3)];
        _S875[int(4)] = _S886[int(4)];
        _S875[int(5)] = _S886[int(5)];
        _S875[int(6)] = _S886[int(6)];
        _S875[int(7)] = _S886[int(7)];
        _S875[int(8)] = _S886[int(8)];
        _S875[int(9)] = _S886[int(9)];
    }
    else
    {
        _S776 = _S794;
        _S758 = 0.0f;
        _S875[int(0)] = 0.0f;
        _S875[int(1)] = 0.0f;
        _S875[int(2)] = 0.0f;
        _S875[int(3)] = 0.0f;
        _S875[int(4)] = 0.0f;
        _S875[int(5)] = 0.0f;
        _S875[int(6)] = 0.0f;
        _S875[int(7)] = 0.0f;
        _S875[int(8)] = 0.0f;
        _S875[int(9)] = 0.0f;
    }
    float _S887 = _S640 * (_S828.x + _S828.y);
    float2  _S888 = _S776 / _S642;
    float2  _S889 = _S639 * - _S888;
    float2  _S890 = _S641 * _S888 + _S826;
    float2  _S891 = _S639 * - _S638;
    float3  _S892 = make_float3 (0.0f, 0.0f, _S891.x + _S891.y);
    float3  _S893 = make_float3 (_S890.x, _S890.y, _S887 + _S887 + _S824.x + _S824.y + _S889.x + _S889.y) + _S892;
    float2  _S894 = _S638;
    *&((&_S894)->y) = _S821.rows[int(0)].y;
    *&((&_S894)->x) = _S821.rows[int(0)].x;
    float2  _S895 = _S641 * _S894;
    float2  _S896 = _S749 * _S894;
    float2  _S897 = _S638;
    *&((&_S897)->y) = _S821.rows[int(0)].z;
    *&((&_S897)->x) = _S821.rows[int(0)].z;
    float2  _S898 = _S751 * _S897;
    float2  _S899 = (_S895 + - (_S639 * _S897)) / _S750;
    float2  _S900 = _S734 * - _S899;
    float2  _S901 = _S642 * _S899;
    if(_S691)
    {
        float _S902 = _S901.x;
        float _S903 = _S902 + _S902;
        float _S904 = _S723 * _S903;
        float _S905 = _S901.y + _S901.y;
        float _S906 = _S723 * _S905;
        float _S907 = _S704 * _S903 + _S698 * _S905;
        float _S908 = _S711 * _S907;
        float _S909 = _S724 * _S907;
        float _S910 = _S725 * _S907;
        float _S911 = _S725 * _S908;
        float _S912 = _S726 * _S907;
        float _S913 = _S710 * _S907 + _S707 * _S908;
        float _S914 = _S726 * _S913;
        float _S915 = _S727 * _S907;
        float _S916 = _S709 * _S907 + _S707 * _S913;
        float _S917 = _S727 * _S916;
        float _S918 = _S708 * _S907 + _S707 * _S916;
        float _S919 = _S706 * _S901.x;
        float _S920 = _S728 * _S901.x;
        float _S921 = _S698 * _S919;
        float _S922 = _S729 * _S919;
        float _S923 = _S705 * _S901.y;
        float _S924 = _S729 * _S901.y;
        float _S925 = 2.0f * _S901.x;
        float _S926 = _S730 * _S925;
        float _S927 = _S730 * _S901.x;
        float _S928 = _S907 + _S704 * _S925 + _S703 * _S901.x;
        float _S929 = _S702 * _S928;
        float _S930 = _S729 * _S928;
        float _S931 = _S701 * _S907;
        float _S932 = _S729 * _S907;
        float _S933 = _S700 * _S901.x;
        float _S934 = _S731 * _S901.x;
        float _S935 = _S698 * _S933;
        float _S936 = _S732 * _S933;
        float _S937 = _S699 * _S901.y;
        float _S938 = _S732 * _S901.y;
        float _S939 = 2.0f * _S901.y;
        float _S940 = _S733 * _S939;
        float _S941 = _S733 * _S901.y;
        float _S942 = _S907 + _S698 * _S939 + _S697 * _S901.y;
        float _S943 = _S696 * _S942;
        float _S944 = _S732 * _S942;
        float _S945 = _S695 * _S907;
        float _S946 = _S732 * _S907;
        float2  _S947 = _S712 * _S901;
        float2  _S948 = _S735 * _S901;
        float2  _S949 = _S638;
        *&((&_S949)->y) = _S918;
        *&((&_S949)->x) = _S918;
        float2  _S950 = _S735 * _S949;
        float _S951 = _S935 + _S937 + _S943 + _S945;
        float _S952 = _S921 + _S923 + _S929 + _S931;
        float2  _S953 = _S947 + uv_19 * _S949;
        float2  _S954 = _S638;
        *&((&_S954)->y) = _S951;
        *&((&_S954)->x) = _S952;
        float2  _S955 = _S953 + _S954;
        float _S956 = _S692 * _S955.x;
        float _S957 = _S692 * _S955.y;
        float _S958 = _S904 + _S926;
        float _S959 = _S906 + _S922 + _S936 + _S940;
        float _S960 = _S911 + _S914 + _S917;
        _S692 = _S694 * _S955.x + _S693 * _S955.y + _S955.x;
        _S693 = _S957;
        _S694 = _S956;
        _S712 = _S948;
        _S695 = _S946;
        _S696 = _S944;
        _S697 = _S941;
        _S699 = _S959;
        _S701 = _S938;
        _S702 = _S934;
        _S703 = _S932;
        _S705 = _S930;
        _S723 = _S927;
        _S724 = _S958;
        _S725 = _S924;
        _S726 = _S920;
        _S727 = _S960;
        _S728 = _S915;
        _S729 = _S912;
        _S730 = _S910;
        _S731 = _S909;
        _S734 = _S950;
    }
    else
    {
        _S692 = 0.0f;
        _S693 = 0.0f;
        _S694 = 0.0f;
        _S712 = _S638;
        _S695 = 0.0f;
        _S696 = 0.0f;
        _S697 = 0.0f;
        _S699 = 0.0f;
        _S701 = 0.0f;
        _S702 = 0.0f;
        _S703 = 0.0f;
        _S705 = 0.0f;
        _S723 = 0.0f;
        _S724 = 0.0f;
        _S725 = 0.0f;
        _S726 = 0.0f;
        _S727 = 0.0f;
        _S728 = 0.0f;
        _S729 = 0.0f;
        _S730 = 0.0f;
        _S731 = 0.0f;
        _S734 = _S638;
    }
    if(_S691)
    {
        float _S961 = _S712.x + _S712.y;
        float _S962 = _S728 + _S707 * _S961;
        float _S963 = _S729 + _S707 * _S962;
        float _S964 = _S730 + _S707 * _S963;
        float _S965 = _S727 + _S708 * _S961 + _S709 * _S962 + _S710 * _S963 + _S711 * _S964;
        float _S966 = _S698 * _S965;
        float _S967 = _S704 * _S965;
        float _S968 = _S705 + 2.0f * (_S702 + _S704 * _S701);
        float _S969 = _S696 + 2.0f * (_S726 + _S704 * _S725);
        float _S970 = _S731 + _S707 * _S964;
        float2  _S971 = _S734 + make_float2 (_S724 + _S700 * _S701 + 2.0f * _S723 + _S706 * _S725 + _S967 + _S967, _S699 + 2.0f * _S697 + _S966 + _S966);
        FixedArray<float, 10>  _S972;
        _S972[int(0)] = 0.0f;
        _S972[int(1)] = 0.0f;
        _S972[int(2)] = 0.0f;
        _S972[int(3)] = 0.0f;
        _S972[int(4)] = 0.0f;
        _S972[int(5)] = 0.0f;
        _S972[int(6)] = 0.0f;
        _S972[int(7)] = 0.0f;
        _S972[int(8)] = 0.0f;
        _S972[int(9)] = 0.0f;
        _S972[int(9)] = _S693;
        _S972[int(8)] = _S694;
        _S972[int(7)] = _S695;
        _S972[int(6)] = _S703;
        _S972[int(5)] = _S968;
        _S972[int(4)] = _S969;
        _S972[int(3)] = _S970;
        _S972[int(2)] = _S964;
        _S972[int(1)] = _S963;
        _S972[int(0)] = _S962;
        float _S973 = _S875[int(0)] + _S972[int(0)];
        float _S974 = _S875[int(1)] + _S972[int(1)];
        float _S975 = _S875[int(2)] + _S972[int(2)];
        float _S976 = _S875[int(3)] + _S972[int(3)];
        float _S977 = _S875[int(4)] + _S972[int(4)];
        float _S978 = _S875[int(5)] + _S972[int(5)];
        float _S979 = _S875[int(6)] + _S972[int(6)];
        float _S980 = _S875[int(7)] + _S972[int(7)];
        float _S981 = _S875[int(8)] + _S972[int(8)];
        float _S982 = _S875[int(9)] + _S972[int(9)];
        _S712 = _S971;
        _S875[int(0)] = _S973;
        _S875[int(1)] = _S974;
        _S875[int(2)] = _S975;
        _S875[int(3)] = _S976;
        _S875[int(4)] = _S977;
        _S875[int(5)] = _S978;
        _S875[int(6)] = _S979;
        _S875[int(7)] = _S980;
        _S875[int(8)] = _S981;
        _S875[int(9)] = _S982;
    }
    else
    {
        _S712 = _S734;
        _S692 = 0.0f;
    }
    float _S983 = _S640 * (_S900.x + _S900.y);
    float2  _S984 = _S712 / _S642;
    float2  _S985 = _S639 * - _S984;
    float2  _S986 = _S641 * _S984 + _S898;
    float3  _S987 = make_float3 (_S986.x, _S986.y, _S983 + _S983 + _S896.x + _S896.y + _S985.x + _S985.y) + _S892 + _S893;
    if(_S645)
    {
        float _S988 = _S635.primal_0 * dpmean2d_1.x;
        float _S989 = _S649 * _S988;
        float _S990 = _S651 * _S988;
        float2  _S991 = make_float2 (_S988, _S636.primal_0 * dpmean2d_1.y) + make_float2 (_S650 * _S988, _S648 * _S988);
        float2  _S992 = uv_18 * _S991;
        float2  _S993 = _S671 * _S991;
        float _S994 = _S653 * _S991.y;
        float _S995 = _S654 * _S991.y;
        float _S996 = _S657 * _S991.y;
        float _S997 = _S653 * _S991.x;
        float _S998 = _S662 * _S991.x;
        float _S999 = _S657 * _S991.x;
        float _S1000 = _S992.x + _S992.y;
        float _S1001 = _S653 * _S1000;
        float _S1002 = _S653 * _S1001;
        float _S1003 = _S653 * _S1002;
        float _S1004 = _S653 * _S1003;
        float _S1005 = _S652 * _S991.y + _S995 + _S661 * _S991.x + _S998 + _S667 * _S1000 + _S668 * _S1001 + _S669 * _S1002 + _S670 * _S1003;
        float _S1006 = _S657 * _S1005;
        float _S1007 = _S660 * _S1005;
        float _S1008 = _S656 * _S995 + 2.0f * (_S657 * _S995) + _S658 * _S991.y + _S665 * _S991.x + _S1006 + _S1006;
        float _S1009 = _S659 * _S996 + _S664 * _S998 + 2.0f * (_S660 * _S998) + _S666 * _S999 + _S1007 + _S1007;
        float _S1010 = 2.0f * (_S660 * _S996) + _S663 * _S991.x;
        float _S1011 = _S655 * _S991.y + 2.0f * (_S660 * _S999);
        float _S1012 = _S647 * dpmean2d_1.x + _S692;
        float _S1013 = _S646 * dpmean2d_1.y + _S758;
        FixedArray<float, 10>  _S1014;
        _S1014[int(0)] = 0.0f;
        _S1014[int(1)] = 0.0f;
        _S1014[int(2)] = 0.0f;
        _S1014[int(3)] = 0.0f;
        _S1014[int(4)] = 0.0f;
        _S1014[int(5)] = 0.0f;
        _S1014[int(6)] = 0.0f;
        _S1014[int(7)] = 0.0f;
        _S1014[int(8)] = 0.0f;
        _S1014[int(9)] = 0.0f;
        _S1014[int(9)] = _S989;
        _S1014[int(8)] = _S990;
        _S1014[int(7)] = _S994;
        _S1014[int(6)] = _S997;
        _S1014[int(5)] = _S1010;
        _S1014[int(4)] = _S1011;
        _S1014[int(3)] = _S1004;
        _S1014[int(2)] = _S1003;
        _S1014[int(1)] = _S1002;
        _S1014[int(0)] = _S1001;
        float _S1015 = _S875[int(0)] + _S1014[int(0)];
        float _S1016 = _S875[int(1)] + _S1014[int(1)];
        float _S1017 = _S875[int(2)] + _S1014[int(2)];
        float _S1018 = _S875[int(3)] + _S1014[int(3)];
        float _S1019 = _S875[int(4)] + _S1014[int(4)];
        float _S1020 = _S875[int(5)] + _S1014[int(5)];
        float _S1021 = _S875[int(6)] + _S1014[int(6)];
        float _S1022 = _S875[int(7)] + _S1014[int(7)];
        float _S1023 = _S875[int(8)] + _S1014[int(8)];
        float _S1024 = _S875[int(9)] + _S1014[int(9)];
        _S671 = _S993 + make_float2 (_S1009, _S1008);
        _S875[int(0)] = _S1015;
        _S875[int(1)] = _S1016;
        _S875[int(2)] = _S1017;
        _S875[int(3)] = _S1018;
        _S875[int(4)] = _S1019;
        _S875[int(5)] = _S1020;
        _S875[int(6)] = _S1021;
        _S875[int(7)] = _S1022;
        _S875[int(8)] = _S1023;
        _S875[int(9)] = _S1024;
        _S646 = dpmean2d_1.y;
        _S647 = dpmean2d_1.x;
        _S648 = _S1013;
        _S649 = _S1012;
    }
    else
    {
        _S671 = _S638;
        _S646 = 0.0f;
        _S647 = 0.0f;
        _S648 = _S758;
        _S649 = _S692;
    }
    float2  _S1025 = _S671 / _S642;
    float2  _S1026 = _S639 * - _S1025;
    float2  _S1027 = _S641 * _S1025;
    float _S1028 = _S1026.x + _S1026.y;
    dpdist_coeffs_1->primal_0 = dpdist_coeffs_1->primal_0;
    dpdist_coeffs_1->differential_0 = _S875;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = _S646;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = _S647;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S648;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S649;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S820.differential_0;
    float3  _S1029 = _S987 + make_float3 (_S1027.x, _S1027.y, _S1028);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1029;
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1030, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1031, Matrix<float, 3, 3>  _S1032)
{
    mul_1(_S1030, _S1031, _S1032);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1033, float3  _S1034)
{
    _d_exp_vector_0(_S1033, _S1034);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_14, float fy_14, float cx_14, float cy_14, FixedArray<float, 10>  * dist_coeffs_20, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    Matrix<float, 2, 2>  _S1035 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S1036 = { false, false, false };
    s_bwd_prop_projection_3dgs_persp_differentiable_Intermediates_0 _S1037;
    (&_S1037)->_S507 = _S1035;
    (&_S1037)->_S508 = _S1036;
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S1038 = s_primal_ctx_exp_0(scale_12);
    float _S1039 = quat_13.y;
    float x2_13 = _S1039 * _S1039;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S1040 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S1038.x, 0.0f, 0.0f, 0.0f, _S1038.y, 0.0f, 0.0f, 0.0f, _S1038.z);
    Matrix<float, 3, 3>  _S1041 = s_primal_ctx_mul_2(_S1040, S_0);
    Matrix<float, 3, 3>  _S1042 = transpose_0(_S1041);
    Matrix<float, 3, 3>  _S1043 = s_primal_ctx_mul_2(_S1041, _S1042);
    Matrix<float, 3, 3>  _S1044 = s_primal_ctx_mul_2(R_14, _S1043);
    Matrix<float, 3, 3>  _S1045 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S1046 = s_primal_ctx_mul_2(_S1044, _S1045);
    Matrix<float, 2, 2>  _S1047 = _S1035;
    float2  _S1048 = make_float2 (0.0f);
    float2  _S1049 = _S1048;
    s_primal_ctx_persp_proj_3dgs_0(mean_c_10, _S1046, fx_14, fy_14, cx_14, cy_14, dist_coeffs_20, &_S1047, &_S1049, &(&_S1037)->_S508);
    (&_S1037)->_S507 = _S1047;
    s_bwd_prop_projection_3dgs_persp_differentiable_Intermediates_0 _S1050 = _S1037;
    float _S1051 = _S1037._S507.rows[int(0)].y * _S1037._S507.rows[int(1)].x;
    float det_orig_11 = _S1037._S507.rows[int(0)].x * _S1037._S507.rows[int(1)].y - _S1051;
    float _S1052 = _S1037._S507.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1053 = _S1037._S507;
    *&(((&_S1053)->rows + (int(0)))->x) = _S1052;
    float _S1054 = _S1037._S507.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1053)->rows + (int(1)))->y) = _S1054;
    Matrix<float, 2, 2>  _S1055 = _S1053;
    Matrix<float, 2, 2>  _S1056 = _S1053;
    float det_blur_6 = _S1052 * _S1054 - _S1051;
    float _S1057 = det_orig_11 / det_blur_6;
    float _S1058 = det_blur_6 * det_blur_6;
    float _S1059 = s_primal_ctx_max_0(0.0f, _S1057);
    float _S1060 = s_primal_ctx_sqrt_0(_S1059);
    float invdet_7 = 1.0f / det_blur_6;
    float _S1061 = - _S1037._S507.rows[int(0)].y;
    float _S1062 = - _S1037._S507.rows[int(1)].x;
    float _S1063 = - in_opacity_10;
    float _S1064 = 1.0f + s_primal_ctx_exp_1(_S1063);
    float _S1065 = 1.0f / _S1064;
    float _S1066 = _S1064 * _S1064;
    float _S1067;
    if(antialiased_10)
    {
        _S1067 = _S1065 * _S1060;
    }
    else
    {
        _S1067 = _S1065;
    }
    float _S1068 = _S1067 / 0.00392156885936856f;
    float _S1069 = 2.0f * s_primal_ctx_log_0(_S1068);
    float _S1070 = s_primal_ctx_sqrt_0(_S1069);
    float _S1071 = _S1055.rows[int(0)].x;
    float _S1072 = _S1056.rows[int(1)].y;
    float _S1073 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S1074 = mean_10 - - s_primal_ctx_mul_1(_S1045, t_13);
    float _S1075 = _S1074.x;
    float _S1076 = _S1074.y;
    float _S1077 = _S1074.z;
    float _S1078 = _S1075 * _S1075 + _S1076 * _S1076 + _S1077 * _S1077;
    float _S1079 = s_primal_ctx_sqrt_0(_S1078);
    float x_34 = _S1075 / _S1079;
    float3  _S1080 = make_float3 (x_34);
    float _S1081 = _S1079 * _S1079;
    float y_13 = _S1076 / _S1079;
    float z_10 = _S1077 / _S1079;
    float3  _S1082 = make_float3 (z_10);
    float _S1083 = - y_13;
    float3  _S1084 = make_float3 (_S1083);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_34 * x_34 - y_13 * y_13;
    float _S1085 = 2.0f * x_34;
    float fS1_10 = _S1085 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S1086 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_34;
    float3  _S1087 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S1088 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S1089 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S1090 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S1091 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S1091;
    float3  _S1092 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_34;
    float3  _S1093 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S1094 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S1095 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S1096 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_10 - y_13 * fS1_10);
    float3  _S1097 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_10 + y_13 * fC1_10);
    float3  _S1098 = make_float3 (pSH9_0);
    float3  _S1099 = make_float3 (0.0f);
    float3  _S1100 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1101;
    (&_S1101)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1083) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_34) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S1101)->differential_0 = _S1100;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1102;
    (&_S1102)->primal_0 = _S1099;
    (&_S1102)->differential_0 = _S1100;
    s_bwd_prop_max_0(&_S1101, &_S1102, v_rgb_0);
    float3  _S1103 = _S1097 * _S1101.differential_0;
    float3  _S1104 = (*sh_coeffs_10)[int(15)] * _S1101.differential_0;
    float3  _S1105 = _S1095 * _S1101.differential_0;
    float3  _S1106 = (*sh_coeffs_10)[int(14)] * _S1101.differential_0;
    float3  _S1107 = _S1093 * _S1101.differential_0;
    float3  _S1108 = (*sh_coeffs_10)[int(13)] * _S1101.differential_0;
    float3  _S1109 = _S1092 * _S1101.differential_0;
    float3  _S1110 = (*sh_coeffs_10)[int(12)] * _S1101.differential_0;
    float3  _S1111 = _S1094 * _S1101.differential_0;
    float3  _S1112 = (*sh_coeffs_10)[int(11)] * _S1101.differential_0;
    float3  _S1113 = _S1096 * _S1101.differential_0;
    float3  _S1114 = (*sh_coeffs_10)[int(10)] * _S1101.differential_0;
    float3  _S1115 = _S1098 * _S1101.differential_0;
    float3  _S1116 = (*sh_coeffs_10)[int(9)] * _S1101.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S1116.x + _S1116.y + _S1116.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S1104.x + _S1104.y + _S1104.z);
    float _S1117 = _S1114.x + _S1114.y + _S1114.z;
    float _S1118 = _S1106.x + _S1106.y + _S1106.z;
    float _S1119 = _S1112.x + _S1112.y + _S1112.z;
    float _S1120 = _S1108.x + _S1108.y + _S1108.z;
    float _S1121 = _S1110.x + _S1110.y + _S1110.z;
    float _S1122 = - s_diff_fC2_T_0;
    float3  _S1123 = _S1089 * _S1101.differential_0;
    float3  _S1124 = (*sh_coeffs_10)[int(8)] * _S1101.differential_0;
    float3  _S1125 = _S1087 * _S1101.differential_0;
    float3  _S1126 = (*sh_coeffs_10)[int(7)] * _S1101.differential_0;
    float3  _S1127 = _S1086 * _S1101.differential_0;
    float3  _S1128 = (*sh_coeffs_10)[int(6)] * _S1101.differential_0;
    float3  _S1129 = _S1088 * _S1101.differential_0;
    float3  _S1130 = (*sh_coeffs_10)[int(5)] * _S1101.differential_0;
    float3  _S1131 = _S1090 * _S1101.differential_0;
    float3  _S1132 = (*sh_coeffs_10)[int(4)] * _S1101.differential_0;
    float _S1133 = _S1130.x + _S1130.y + _S1130.z;
    float _S1134 = _S1126.x + _S1126.y + _S1126.z;
    float _S1135 = fTmp1B_10 * _S1117 + x_34 * s_diff_fS2_T_0 + y_13 * _S1122 + 0.54627424478530884f * (_S1132.x + _S1132.y + _S1132.z);
    float _S1136 = fTmp1B_10 * _S1118 + y_13 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S1124.x + _S1124.y + _S1124.z);
    float _S1137 = y_13 * - _S1136;
    float _S1138 = x_34 * _S1136;
    float _S1139 = z_10 * (1.86588168144226074f * (z_10 * _S1121) + -2.28522896766662598f * (y_13 * _S1119 + x_34 * _S1120) + 0.94617468118667603f * (_S1128.x + _S1128.y + _S1128.z));
    float3  _S1140 = make_float3 (0.48860251903533936f) * _S1101.differential_0;
    float3  _S1141 = - _S1140;
    float3  _S1142 = _S1080 * _S1141;
    float3  _S1143 = (*sh_coeffs_10)[int(3)] * _S1141;
    float3  _S1144 = _S1082 * _S1140;
    float3  _S1145 = (*sh_coeffs_10)[int(2)] * _S1140;
    float3  _S1146 = _S1084 * _S1140;
    float3  _S1147 = (*sh_coeffs_10)[int(1)] * _S1140;
    float _S1148 = (_S1091 * _S1121 + 1.44530570507049561f * (fS1_10 * _S1117 + fC1_10 * _S1118) + -1.09254848957061768f * (y_13 * _S1133 + x_34 * _S1134) + _S1139 + _S1139 + _S1145.x + _S1145.y + _S1145.z) / _S1081;
    float _S1149 = _S1079 * _S1148;
    float _S1150 = (fTmp0C_10 * _S1119 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S1122 + fTmp0B_10 * _S1133 + _S1085 * _S1135 + _S1137 + _S1137 + - (_S1147.x + _S1147.y + _S1147.z)) / _S1081;
    float _S1151 = _S1079 * _S1150;
    float _S1152 = (fTmp0C_10 * _S1120 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S1134 + 2.0f * (y_13 * _S1135) + _S1138 + _S1138 + _S1143.x + _S1143.y + _S1143.z) / _S1081;
    float _S1153 = _S1079 * _S1152;
    float _S1154 = _S1077 * - _S1148 + _S1076 * - _S1150 + _S1075 * - _S1152;
    DiffPair_float_0 _S1155;
    (&_S1155)->primal_0 = _S1078;
    (&_S1155)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1155, _S1154);
    float _S1156 = _S1077 * _S1155.differential_0;
    float _S1157 = _S1076 * _S1155.differential_0;
    float _S1158 = _S1075 * _S1155.differential_0;
    float3  _S1159 = make_float3 (0.282094806432724f) * _S1101.differential_0;
    float3  _S1160 = make_float3 (_S1153 + _S1158 + _S1158, _S1151 + _S1157 + _S1157, _S1149 + _S1156 + _S1156);
    float3  _S1161 = - - _S1160;
    Matrix<float, 3, 3>  _S1162 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1163;
    (&_S1163)->primal_0 = _S1045;
    (&_S1163)->differential_0 = _S1162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1164;
    (&_S1164)->primal_0 = t_13;
    (&_S1164)->differential_0 = _S1100;
    s_bwd_prop_mul_1(&_S1163, &_S1164, _S1161);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1165 = _S1163;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1166 = _S1164;
    float2  _S1167 = _S1048;
    *&((&_S1167)->y) = v_conic_0.z;
    float2  _S1168 = _S1048;
    *&((&_S1168)->y) = v_conic_0.y;
    *&((&_S1168)->x) = v_conic_0.x;
    float _S1169 = 0.5f * v_depth_0;
    DiffPair_float_0 _S1170;
    (&_S1170)->primal_0 = _S1073;
    (&_S1170)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1170, _S1169);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1171;
    (&_S1171)->primal_0 = mean_c_10;
    (&_S1171)->differential_0 = _S1100;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1172;
    (&_S1172)->primal_0 = mean_c_10;
    (&_S1172)->differential_0 = _S1100;
    s_bwd_prop_dot_0(&_S1171, &_S1172, _S1170.differential_0);
    DiffPair_float_0 _S1173;
    (&_S1173)->primal_0 = _S1072;
    (&_S1173)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1173, 0.0f);
    DiffPair_float_0 _S1174;
    (&_S1174)->primal_0 = _S1071;
    (&_S1174)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1174, 0.0f);
    DiffPair_float_0 _S1175;
    (&_S1175)->primal_0 = 3.32999992370605469f;
    (&_S1175)->differential_0 = 0.0f;
    DiffPair_float_0 _S1176;
    (&_S1176)->primal_0 = _S1070;
    (&_S1176)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1175, &_S1176, 0.0f);
    DiffPair_float_0 _S1177;
    (&_S1177)->primal_0 = _S1069;
    (&_S1177)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1177, _S1176.differential_0);
    float _S1178 = 2.0f * _S1177.differential_0;
    DiffPair_float_0 _S1179;
    (&_S1179)->primal_0 = _S1068;
    (&_S1179)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1179, _S1178);
    float2  _S1180 = make_float2 (_S1174.differential_0, 0.0f);
    float _S1181 = v_opacity_0 + 254.9999847412109375f * _S1179.differential_0;
    Matrix<float, 2, 2>  _S1182 = _S1035;
    _S1182[int(1)] = _S1167;
    _S1182[int(0)] = _S1168;
    Matrix<float, 2, 2>  _S1183 = _S1182;
    FixedArray<float3 , 16>  _S1184;
    _S1184[int(0)] = _S1100;
    _S1184[int(1)] = _S1100;
    _S1184[int(2)] = _S1100;
    _S1184[int(3)] = _S1100;
    _S1184[int(4)] = _S1100;
    _S1184[int(5)] = _S1100;
    _S1184[int(6)] = _S1100;
    _S1184[int(7)] = _S1100;
    _S1184[int(8)] = _S1100;
    _S1184[int(9)] = _S1100;
    _S1184[int(10)] = _S1100;
    _S1184[int(11)] = _S1100;
    _S1184[int(12)] = _S1100;
    _S1184[int(13)] = _S1100;
    _S1184[int(14)] = _S1100;
    _S1184[int(15)] = _S1100;
    _S1184[int(7)] = _S1125;
    _S1184[int(0)] = _S1159;
    _S1184[int(1)] = _S1146;
    _S1184[int(2)] = _S1144;
    _S1184[int(3)] = _S1142;
    _S1184[int(4)] = _S1131;
    _S1184[int(5)] = _S1129;
    _S1184[int(6)] = _S1127;
    _S1184[int(15)] = _S1103;
    _S1184[int(8)] = _S1123;
    _S1184[int(9)] = _S1115;
    _S1184[int(10)] = _S1113;
    _S1184[int(11)] = _S1111;
    _S1184[int(12)] = _S1109;
    _S1184[int(13)] = _S1107;
    _S1184[int(14)] = _S1105;
    float3  _S1185 = _S1184[int(0)];
    float3  _S1186 = _S1184[int(1)];
    float3  _S1187 = _S1184[int(2)];
    float3  _S1188 = _S1184[int(3)];
    float3  _S1189 = _S1184[int(4)];
    float3  _S1190 = _S1184[int(5)];
    float3  _S1191 = _S1184[int(6)];
    float3  _S1192 = _S1184[int(7)];
    float3  _S1193 = _S1184[int(8)];
    float3  _S1194 = _S1184[int(9)];
    float3  _S1195 = _S1184[int(10)];
    float3  _S1196 = _S1184[int(11)];
    float3  _S1197 = _S1184[int(12)];
    float3  _S1198 = _S1184[int(13)];
    float3  _S1199 = _S1184[int(14)];
    float3  _S1200 = _S1184[int(15)];
    float3  _S1201 = _S1172.differential_0 + _S1171.differential_0;
    float2  _S1202 = make_float2 (0.0f, _S1173.differential_0);
    float _S1203;
    if(antialiased_10)
    {
        float _S1204 = _S1065 * _S1181;
        _S1067 = _S1060 * _S1181;
        _S1203 = _S1204;
    }
    else
    {
        _S1067 = _S1181;
        _S1203 = 0.0f;
    }
    float _S1205 = - (_S1067 / _S1066);
    DiffPair_float_0 _S1206;
    (&_S1206)->primal_0 = _S1063;
    (&_S1206)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1206, _S1205);
    float _S1207 = - _S1206.differential_0;
    float _S1208 = invdet_7 * _S1183.rows[int(1)].y;
    float _S1209 = - (invdet_7 * _S1183.rows[int(1)].x);
    float _S1210 = - (invdet_7 * _S1183.rows[int(0)].y);
    float _S1211 = invdet_7 * _S1183.rows[int(0)].x;
    float _S1212 = - ((_S1052 * _S1183.rows[int(1)].y + _S1062 * _S1183.rows[int(1)].x + _S1061 * _S1183.rows[int(0)].y + _S1054 * _S1183.rows[int(0)].x) / _S1058);
    DiffPair_float_0 _S1213;
    (&_S1213)->primal_0 = _S1059;
    (&_S1213)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1213, _S1203);
    DiffPair_float_0 _S1214;
    (&_S1214)->primal_0 = 0.0f;
    (&_S1214)->differential_0 = 0.0f;
    DiffPair_float_0 _S1215;
    (&_S1215)->primal_0 = _S1057;
    (&_S1215)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1214, &_S1215, _S1213.differential_0);
    float _S1216 = _S1215.differential_0 / _S1058;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1216;
    float _S1217 = _S1212 + det_orig_11 * - _S1216;
    float _S1218 = - _S1217;
    float _S1219 = _S1052 * _S1217;
    float _S1220 = _S1054 * _S1217;
    Matrix<float, 2, 2>  _S1221 = _S1035;
    _S1221[int(1)] = _S1202;
    _S1221[int(0)] = _S1180;
    _S1053 = _S1221;
    *&(((&_S1053)->rows + (int(1)))->y) = 0.0f;
    float _S1222 = _S1211 + _S1219 + _S1221.rows[int(1)].y;
    *&(((&_S1053)->rows + (int(0)))->x) = 0.0f;
    float _S1223 = _S1208 + _S1220 + _S1221.rows[int(0)].x;
    float _S1224 = _S1218 + - s_diff_det_orig_T_0;
    float _S1225 = _S1209 + _S1050._S507.rows[int(0)].y * _S1224;
    float _S1226 = _S1210 + _S1050._S507.rows[int(1)].x * _S1224;
    float _S1227 = _S1050._S507.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1228 = _S1222 + _S1050._S507.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1229 = _S1048;
    *&((&_S1229)->x) = _S1225;
    *&((&_S1229)->y) = _S1228;
    float _S1230 = _S1223 + _S1227;
    float2  _S1231 = _S1048;
    *&((&_S1231)->y) = _S1226;
    *&((&_S1231)->x) = _S1230;
    Matrix<float, 2, 2>  _S1232 = _S1035;
    _S1232[int(1)] = _S1229;
    _S1232[int(0)] = _S1231;
    Matrix<float, 2, 2>  _S1233 = _S1053 + _S1232;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
    (&_S1234)->primal_0 = mean_c_10;
    (&_S1234)->differential_0 = _S1100;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1235;
    (&_S1235)->primal_0 = _S1046;
    (&_S1235)->differential_0 = _S1162;
    DiffPair_float_0 _S1236;
    (&_S1236)->primal_0 = fx_14;
    (&_S1236)->differential_0 = 0.0f;
    DiffPair_float_0 _S1237;
    (&_S1237)->primal_0 = fy_14;
    (&_S1237)->differential_0 = 0.0f;
    DiffPair_float_0 _S1238;
    (&_S1238)->primal_0 = cx_14;
    (&_S1238)->differential_0 = 0.0f;
    DiffPair_float_0 _S1239;
    (&_S1239)->primal_0 = cy_14;
    (&_S1239)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S1240 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S1241;
    (&_S1241)->primal_0 = *dist_coeffs_20;
    (&_S1241)->differential_0 = _S1240;
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S1242 = _S1050._S508;
    s_bwd_prop_persp_proj_3dgs_0(&_S1234, &_S1235, &_S1236, &_S1237, &_S1238, &_S1239, &_S1241, _S1233, v_mean2d_0, &_S1242);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1243;
    (&_S1243)->primal_0 = _S1044;
    (&_S1243)->differential_0 = _S1162;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1244;
    (&_S1244)->primal_0 = _S1045;
    (&_S1244)->differential_0 = _S1162;
    s_bwd_prop_mul_4(&_S1243, &_S1244, _S1235.differential_0);
    Matrix<float, 3, 3>  _S1245 = transpose_0(_S1244.differential_0 + _S1165.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1246;
    (&_S1246)->primal_0 = R_14;
    (&_S1246)->differential_0 = _S1162;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1247;
    (&_S1247)->primal_0 = _S1043;
    (&_S1247)->differential_0 = _S1162;
    s_bwd_prop_mul_4(&_S1246, &_S1247, _S1243.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1248;
    (&_S1248)->primal_0 = _S1041;
    (&_S1248)->differential_0 = _S1162;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1249;
    (&_S1249)->primal_0 = _S1042;
    (&_S1249)->differential_0 = _S1162;
    s_bwd_prop_mul_4(&_S1248, &_S1249, _S1247.differential_0);
    Matrix<float, 3, 3>  _S1250 = _S1248.differential_0 + transpose_0(_S1249.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1251;
    (&_S1251)->primal_0 = _S1040;
    (&_S1251)->differential_0 = _S1162;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1252;
    (&_S1252)->primal_0 = S_0;
    (&_S1252)->differential_0 = _S1162;
    s_bwd_prop_mul_4(&_S1251, &_S1252, _S1250);
    Matrix<float, 3, 3>  _S1253 = transpose_0(_S1251.differential_0);
    float _S1254 = 2.0f * - _S1253.rows[int(2)].z;
    float _S1255 = 2.0f * _S1253.rows[int(2)].y;
    float _S1256 = 2.0f * _S1253.rows[int(2)].x;
    float _S1257 = 2.0f * _S1253.rows[int(1)].z;
    float _S1258 = 2.0f * - _S1253.rows[int(1)].y;
    float _S1259 = 2.0f * _S1253.rows[int(1)].x;
    float _S1260 = 2.0f * _S1253.rows[int(0)].z;
    float _S1261 = 2.0f * _S1253.rows[int(0)].y;
    float _S1262 = 2.0f * - _S1253.rows[int(0)].x;
    float _S1263 = - _S1259 + _S1261;
    float _S1264 = _S1256 + - _S1260;
    float _S1265 = - _S1255 + _S1257;
    float _S1266 = _S1255 + _S1257;
    float _S1267 = _S1256 + _S1260;
    float _S1268 = _S1259 + _S1261;
    float _S1269 = quat_13.w * (_S1258 + _S1262);
    float _S1270 = quat_13.z * (_S1254 + _S1262);
    float _S1271 = quat_13.y * (_S1254 + _S1258);
    float _S1272 = quat_13.x * _S1263 + quat_13.z * _S1266 + quat_13.y * _S1267 + _S1269 + _S1269;
    float _S1273 = quat_13.x * _S1264 + quat_13.w * _S1266 + quat_13.y * _S1268 + _S1270 + _S1270;
    float _S1274 = quat_13.x * _S1265 + quat_13.w * _S1267 + quat_13.z * _S1268 + _S1271 + _S1271;
    float _S1275 = quat_13.w * _S1263 + quat_13.z * _S1264 + quat_13.y * _S1265;
    float3  _S1276 = _S1100;
    *&((&_S1276)->z) = _S1252.differential_0.rows[int(2)].z;
    *&((&_S1276)->y) = _S1252.differential_0.rows[int(1)].y;
    *&((&_S1276)->x) = _S1252.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1277;
    (&_S1277)->primal_0 = scale_12;
    (&_S1277)->differential_0 = _S1100;
    s_bwd_prop_exp_1(&_S1277, _S1276);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1278 = _S1277;
    float3  _S1279 = _S1234.differential_0 + _S1201;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1280;
    (&_S1280)->primal_0 = R_14;
    (&_S1280)->differential_0 = _S1162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1281;
    (&_S1281)->primal_0 = mean_10;
    (&_S1281)->differential_0 = _S1100;
    s_bwd_prop_mul_1(&_S1280, &_S1281, _S1279);
    float3  _S1282 = _S1279 + _S1166.differential_0;
    Matrix<float, 3, 3>  _S1283 = _S1245 + _S1246.differential_0 + _S1280.differential_0;
    float4  _S1284 = make_float4 (0.0f);
    *&((&_S1284)->w) = _S1272;
    *&((&_S1284)->z) = _S1273;
    *&((&_S1284)->y) = _S1274;
    *&((&_S1284)->x) = _S1275;
    float4  _S1285 = _S1284;
    float3  _S1286 = _S1281.differential_0 + _S1160;
    *v_mean_0 = _S1286;
    *v_quat_0 = _S1285;
    *v_scale_0 = _S1278.differential_0;
    *v_in_opacity_0 = _S1207;
    (*v_sh_coeffs_0)[int(0)] = _S1185;
    (*v_sh_coeffs_0)[int(1)] = _S1186;
    (*v_sh_coeffs_0)[int(2)] = _S1187;
    (*v_sh_coeffs_0)[int(3)] = _S1188;
    (*v_sh_coeffs_0)[int(4)] = _S1189;
    (*v_sh_coeffs_0)[int(5)] = _S1190;
    (*v_sh_coeffs_0)[int(6)] = _S1191;
    (*v_sh_coeffs_0)[int(7)] = _S1192;
    (*v_sh_coeffs_0)[int(8)] = _S1193;
    (*v_sh_coeffs_0)[int(9)] = _S1194;
    (*v_sh_coeffs_0)[int(10)] = _S1195;
    (*v_sh_coeffs_0)[int(11)] = _S1196;
    (*v_sh_coeffs_0)[int(12)] = _S1197;
    (*v_sh_coeffs_0)[int(13)] = _S1198;
    (*v_sh_coeffs_0)[int(14)] = _S1199;
    (*v_sh_coeffs_0)[int(15)] = _S1200;
    *v_R_1 = _S1283;
    *v_t_1 = _S1282;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S1287;
    DiffPair_float_0 _S1288;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S1289;
    DiffPair_float_0 _S1290;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1291;
    DiffPair_float_0 _S1292;
    DiffPair_float_0 _S1293;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1294;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S1295;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1296;
};

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S1297, float _S1298)
{
    return s_primal_ctx_atan2_0(_S1297, _S1298);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S1299;
    DiffPair_float_0 _S1300;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S1301 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S1299 = _S1301;
    _s_diff_ctx_4->_S1300 = _S1301;
    (&_s_diff_ctx_4->_S1299)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1299)->differential_0 = 0.0f;
    (&_s_diff_ctx_4->_S1300)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1300)->differential_0 = 0.0f;
    DiffPair_float_0 _S1302 = *dpdpy_0;
    _s_diff_ctx_4->_S1299 = *dpdpy_0;
    DiffPair_float_0 _S1303 = *dpdpx_0;
    _s_diff_ctx_4->_S1300 = *dpdpx_0;
    float _S1304 = _S1303.primal_0 * _S1303.primal_0 + _S1302.primal_0 * _S1302.primal_0;
    float _S1305 = - _S1302.primal_0 / _S1304 * dpdOut_0;
    float _S1306 = _S1303.primal_0 / _S1304 * dpdOut_0;
    dpdpy_0->primal_0 = _S1302.primal_0;
    dpdpy_0->differential_0 = _S1306;
    dpdpx_0->primal_0 = _S1303.primal_0;
    dpdpx_0->differential_0 = _S1305;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S1307, DiffPair_float_0 * _S1308, float _S1309, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S1310 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S1287 = _S1310;
    _s_diff_ctx_5->_S1288 = _S1310;
    (&_s_diff_ctx_5->_S1287)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S1287)->differential_0 = 0.0f;
    (&_s_diff_ctx_5->_S1288)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S1288)->differential_0 = 0.0f;
    DiffPair_float_0 _S1311 = *_S1307;
    _s_diff_ctx_5->_S1287 = *_S1307;
    DiffPair_float_0 _S1312 = *_S1308;
    _s_diff_ctx_5->_S1288 = *_S1308;
    DiffPair_float_0 _S1313 = _S1311;
    DiffPair_float_0 _S1314 = _S1312;
    s_bwd_prop_d_atan2_Intermediates_0 _S1315;
    (&_S1315)->_S1299 = _S1310;
    (&_S1315)->_S1300 = _S1310;
    s_primal_ctx_d_atan2_0(&_S1313, &_S1314, _S1309, &_S1315);
    *_S1307 = _S1313;
    *_S1308 = _S1314;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1316;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1317;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1318;
    DiffPair_float_0 _S1319;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1320;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1321;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_6)
{
    DiffPair_float_0 _S1322 = { 0.0f, 0.0f };
    _s_diff_ctx_6->_S1321 = _S1322;
    (&_s_diff_ctx_6->_S1321)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1321)->differential_0 = 0.0f;
    DiffPair_float_0 _S1323 = *dpdpx_1;
    _s_diff_ctx_6->_S1321 = *dpdpx_1;
    float _S1324 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S1323.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S1323.primal_0;
    dpdpx_1->differential_0 = _S1324;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1325, float _S1326, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_7)
{
    DiffPair_float_0 _S1327 = { 0.0f, 0.0f };
    _s_diff_ctx_7->_S1317 = _S1327;
    (&_s_diff_ctx_7->_S1317)->primal_0 = 0.0f;
    (&_s_diff_ctx_7->_S1317)->differential_0 = 0.0f;
    DiffPair_float_0 _S1328 = *_S1325;
    _s_diff_ctx_7->_S1317 = *_S1325;
    DiffPair_float_0 _S1329 = _S1328;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1330;
    (&_S1330)->_S1321 = _S1327;
    s_primal_ctx_d_sqrt_0(&_S1329, _S1326, &_S1330);
    *_S1325 = _S1329;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_8)
{
    float2  _S1331 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1332 = { _S1331, _S1331 };
    DiffPair_float_0 _S1333 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1334 = { _S1333 };
    _s_diff_ctx_8->_S1318 = _S1332;
    _s_diff_ctx_8->_S1319 = _S1333;
    _s_diff_ctx_8->_S1320 = _S1334;
    (&_s_diff_ctx_8->_S1318)->primal_0 = _S1331;
    (&_s_diff_ctx_8->_S1318)->differential_0 = _S1331;
    (&_s_diff_ctx_8->_S1319)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1319)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1335 = *dpdpx_2;
    _s_diff_ctx_8->_S1318 = *dpdpx_2;
    float _S1336 = _S1335.primal_0.x;
    float _S1337 = _S1335.primal_0.y;
    DiffPair_float_0 _S1338;
    (&_S1338)->primal_0 = _S1336 * _S1336 + _S1337 * _S1337;
    (&_S1338)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S1338, dp_s_dOut_0, &_s_diff_ctx_8->_S1320);
    _s_diff_ctx_8->_S1319 = _S1338;
    float _S1339 = _S1335.primal_0.y * _S1338.differential_0;
    float _S1340 = _S1339 + _S1339;
    float _S1341 = _S1335.primal_0.x * _S1338.differential_0;
    float _S1342 = _S1341 + _S1341;
    float2  _S1343 = _S1331;
    *&((&_S1343)->y) = _S1340;
    *&((&_S1343)->x) = _S1342;
    dpdpx_2->primal_0 = _S1335.primal_0;
    dpdpx_2->differential_0 = _S1343;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1344, float _S1345, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_9)
{
    float2  _S1346 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1347 = { _S1346, _S1346 };
    _s_diff_ctx_9->_S1316 = _S1347;
    (&_s_diff_ctx_9->_S1316)->primal_0 = _S1346;
    (&_s_diff_ctx_9->_S1316)->differential_0 = _S1346;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1348 = *_S1344;
    _s_diff_ctx_9->_S1316 = *_S1344;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1349 = _S1348;
    DiffPair_float_0 _S1350 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1351 = { _S1350 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1352;
    (&_S1352)->_S1318 = _S1347;
    (&_S1352)->_S1319 = _S1350;
    (&_S1352)->_S1320 = _S1351;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1349, _S1345, &_S1352);
    *_S1344 = _S1349;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_2, Matrix<float, 3, 3>  dpcov3d_2, float dpfx_2, float dpfy_2, float dpcx_2, float dpcy_2, FixedArray<float, 10>  * dpdist_coeffs_2, Matrix<float, 2, 2>  * dpcov2d_2, float2  * dpmean2d_2, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_float_0 _S1353 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1354 = { _S1353, _S1353 };
    _s_diff_ctx_10->_S1289 = _S1353;
    _s_diff_ctx_10->_S1290 = _S1353;
    _s_diff_ctx_10->_S1291 = _S1354;
    _s_diff_ctx_10->_S1292 = _S1353;
    _s_diff_ctx_10->_S1293 = _S1353;
    _s_diff_ctx_10->_S1294 = _S1354;
    (&_s_diff_ctx_10->_S1289)->primal_0 = 0.0f;
    (&_s_diff_ctx_10->_S1289)->differential_0 = 0.0f;
    (&_s_diff_ctx_10->_S1290)->primal_0 = 0.0f;
    (&_s_diff_ctx_10->_S1290)->differential_0 = 0.0f;
    (&_s_diff_ctx_10->_S1292)->primal_0 = 0.0f;
    (&_s_diff_ctx_10->_S1292)->differential_0 = 0.0f;
    (&_s_diff_ctx_10->_S1293)->primal_0 = 0.0f;
    (&_s_diff_ctx_10->_S1293)->differential_0 = 0.0f;
    float2  _S1355 = make_float2 (0.0f);
    float2  _S1356 = float2 {dpmean3d_2.x, dpmean3d_2.y};
    float _S1357 = length_0(_S1356);
    float _S1358 = dpmean3d_2.z;
    float _S1359 = s_primal_ctx_atan2_0(_S1357, _S1358);
    float k_1;
    if(_S1359 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S1359 * _S1359 / 3.0f) / _S1358;
    }
    else
    {
        k_1 = _S1359 / _S1357;
    }
    float2  _S1360 = _S1356 * make_float2 (k_1);
    float u_16 = _S1360.x;
    float v_16 = _S1360.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S1361 = 2.0f * (*dpdist_coeffs_2)[int(4)];
    float _S1362 = 2.0f * (*dpdist_coeffs_2)[int(5)];
    float2  _S1363 = _S1360 * make_float2 (1.0f + r2_16 * ((*dpdist_coeffs_2)[int(0)] + r2_16 * ((*dpdist_coeffs_2)[int(1)] + r2_16 * ((*dpdist_coeffs_2)[int(2)] + r2_16 * (*dpdist_coeffs_2)[int(3)])))) + make_float2 (_S1361 * u_16 * v_16 + (*dpdist_coeffs_2)[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + (*dpdist_coeffs_2)[int(6)] * r2_16, _S1362 * u_16 * v_16 + (*dpdist_coeffs_2)[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + (*dpdist_coeffs_2)[int(7)] * r2_16);
    float2  _S1364 = _S1363 + make_float2 ((*dpdist_coeffs_2)[int(8)] * _S1363.x + (*dpdist_coeffs_2)[int(9)] * _S1363.y, 0.0f);
    float2  _S1365 = make_float2 (dpfx_2 * _S1364.x + dpcx_2, dpfy_2 * _S1364.y + dpcy_2);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S1366 = s_primal_ctx_s_primal_ctx_atan2_0(_S1357, _S1358);
    bool _S1367 = _S1366 < 0.00100000004749745f;
    float _S1368;
    float _S1369;
    float _S1370;
    if(_S1367)
    {
        float _S1371 = 1.0f - _S1366 * _S1366 / 3.0f;
        float _S1372 = _S1358 * _S1358;
        k_1 = _S1371 / _S1358;
        _S1368 = 0.0f;
        _S1369 = _S1372;
        _S1370 = _S1371;
    }
    else
    {
        float _S1373 = _S1357 * _S1357;
        k_1 = _S1366 / _S1357;
        _S1368 = _S1373;
        _S1369 = 0.0f;
        _S1370 = 0.0f;
    }
    float2  _S1374 = make_float2 (k_1);
    float2  _S1375 = _S1356 * make_float2 (k_1);
    float u_17 = _S1375.x;
    float v_17 = _S1375.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S1376 = (*dpdist_coeffs_2)[int(2)] + r2_17 * (*dpdist_coeffs_2)[int(3)];
    float _S1377 = (*dpdist_coeffs_2)[int(1)] + r2_17 * _S1376;
    float _S1378 = (*dpdist_coeffs_2)[int(0)] + r2_17 * _S1377;
    float2  _S1379 = make_float2 (dpfx_2, 0.0f) + make_float2 ((*dpdist_coeffs_2)[int(8)] * dpfx_2, (*dpdist_coeffs_2)[int(9)] * dpfx_2);
    float2  _S1380 = _S1375 * _S1379;
    float _S1381 = (*dpdist_coeffs_2)[int(4)] * _S1379.y;
    float _S1382 = (*dpdist_coeffs_2)[int(5)] * _S1379.x;
    float _S1383 = _S1380.x + _S1380.y;
    float _S1384 = r2_17 * _S1383;
    float _S1385 = r2_17 * _S1384;
    float _S1386 = (*dpdist_coeffs_2)[int(7)] * _S1379.y + _S1381 + (*dpdist_coeffs_2)[int(6)] * _S1379.x + _S1382 + _S1378 * _S1383 + _S1377 * _S1384 + _S1376 * _S1385 + (*dpdist_coeffs_2)[int(3)] * (r2_17 * _S1385);
    float _S1387 = v_17 * _S1386;
    float _S1388 = u_17 * _S1386;
    float2  _S1389 = make_float2 (1.0f + r2_17 * _S1378) * _S1379 + make_float2 (_S1362 * (v_17 * _S1379.y) + 2.0f * u_17 * _S1382 + 2.0f * (u_17 * _S1382) + _S1361 * (v_17 * _S1379.x) + _S1388 + _S1388, 2.0f * v_17 * _S1381 + 2.0f * (v_17 * _S1381) + _S1362 * u_17 * _S1379.y + _S1361 * u_17 * _S1379.x + _S1387 + _S1387);
    float2  _S1390 = _S1356 * _S1389;
    float2  _S1391 = _S1374 * _S1389;
    float _S1392 = _S1390.x + _S1390.y;
    if(_S1367)
    {
        float _S1393 = _S1392 / _S1369;
        float _S1394 = _S1370 * - _S1393;
        float _S1395 = _S1366 * (0.3333333432674408f * - (_S1358 * _S1393));
        k_1 = _S1395 + _S1395;
        _S1368 = _S1394;
        _S1369 = 0.0f;
    }
    else
    {
        float _S1396 = _S1392 / _S1368;
        float _S1397 = _S1366 * - _S1396;
        k_1 = _S1357 * _S1396;
        _S1368 = 0.0f;
        _S1369 = _S1397;
    }
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1357;
    (&_S1398)->differential_0 = 0.0f;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1358;
    (&_S1399)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1398, &_S1399, k_1, &_s_diff_ctx_10->_S1291);
    _s_diff_ctx_10->_S1289 = _S1398;
    _s_diff_ctx_10->_S1290 = _S1399;
    float _S1400 = _S1399.differential_0 + _S1368;
    float _S1401 = _S1398.differential_0 + _S1369;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1402;
    (&_S1402)->primal_0 = _S1356;
    (&_S1402)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1403 = { _S1355, _S1355 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1404;
    (&_S1404)->_S1316 = _S1403;
    s_primal_ctx_s_bwd_length_impl_0(&_S1402, _S1401, &_S1404);
    float2  _S1405 = _S1402.differential_0 + _S1391;
    float3  _S1406 = make_float3 (_S1405.x, _S1405.y, _S1400);
    Matrix<float, 2, 3>  _S1407 = J_10;
    _S1407[int(0)] = _S1406;
    if(_S1367)
    {
        float _S1408 = 1.0f - _S1366 * _S1366 / 3.0f;
        float _S1409 = _S1358 * _S1358;
        k_1 = _S1408 / _S1358;
        _S1368 = 0.0f;
        _S1369 = _S1409;
        _S1370 = _S1408;
    }
    else
    {
        float _S1410 = _S1357 * _S1357;
        k_1 = _S1366 / _S1357;
        _S1368 = _S1410;
        _S1369 = 0.0f;
        _S1370 = 0.0f;
    }
    float2  _S1411 = make_float2 (k_1);
    float2  _S1412 = _S1356 * make_float2 (k_1);
    float u_18 = _S1412.x;
    float v_18 = _S1412.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S1413 = (*dpdist_coeffs_2)[int(2)] + r2_18 * (*dpdist_coeffs_2)[int(3)];
    float _S1414 = (*dpdist_coeffs_2)[int(1)] + r2_18 * _S1413;
    float _S1415 = (*dpdist_coeffs_2)[int(0)] + r2_18 * _S1414;
    float2  _S1416 = make_float2 (0.0f, dpfy_2);
    float2  _S1417 = _S1412 * _S1416;
    float _S1418 = (*dpdist_coeffs_2)[int(4)] * dpfy_2;
    float _S1419 = _S1417.x + _S1417.y;
    float _S1420 = r2_18 * _S1419;
    float _S1421 = r2_18 * _S1420;
    float _S1422 = (*dpdist_coeffs_2)[int(7)] * dpfy_2 + _S1418 + _S1415 * _S1419 + _S1414 * _S1420 + _S1413 * _S1421 + (*dpdist_coeffs_2)[int(3)] * (r2_18 * _S1421);
    float _S1423 = v_18 * _S1422;
    float _S1424 = u_18 * _S1422;
    float2  _S1425 = make_float2 (1.0f + r2_18 * _S1415) * _S1416 + make_float2 (_S1362 * (v_18 * dpfy_2) + _S1424 + _S1424, 2.0f * v_18 * _S1418 + 2.0f * (v_18 * _S1418) + _S1362 * u_18 * dpfy_2 + _S1423 + _S1423);
    float2  _S1426 = _S1356 * _S1425;
    float2  _S1427 = _S1411 * _S1425;
    float _S1428 = _S1426.x + _S1426.y;
    if(_S1367)
    {
        float _S1429 = _S1428 / _S1369;
        float _S1430 = _S1370 * - _S1429;
        float _S1431 = _S1366 * (0.3333333432674408f * - (_S1358 * _S1429));
        k_1 = _S1431 + _S1431;
        _S1368 = _S1430;
        _S1369 = 0.0f;
    }
    else
    {
        float _S1432 = _S1428 / _S1368;
        float _S1433 = _S1366 * - _S1432;
        k_1 = _S1357 * _S1432;
        _S1368 = 0.0f;
        _S1369 = _S1433;
    }
    DiffPair_float_0 _S1434;
    (&_S1434)->primal_0 = _S1357;
    (&_S1434)->differential_0 = 0.0f;
    DiffPair_float_0 _S1435;
    (&_S1435)->primal_0 = _S1358;
    (&_S1435)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1434, &_S1435, k_1, &_s_diff_ctx_10->_S1294);
    _s_diff_ctx_10->_S1292 = _S1434;
    _s_diff_ctx_10->_S1293 = _S1435;
    float _S1436 = _S1435.differential_0 + _S1368;
    float _S1437 = _S1434.differential_0 + _S1369;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1438;
    (&_S1438)->primal_0 = _S1356;
    (&_S1438)->differential_0 = _S1355;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1439;
    (&_S1439)->_S1316 = _S1403;
    s_primal_ctx_s_bwd_length_impl_0(&_S1438, _S1437, &_S1439);
    float2  _S1440 = _S1438.differential_0 + _S1427;
    float3  _S1441 = make_float3 (_S1440.x, _S1440.y, _S1436);
    _S1407[int(1)] = _S1441;
    *dpcov2d_2 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S1407, dpcov3d_2), transpose_1(_S1407));
    *dpmean2d_2 = _S1365;
    return;
}

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

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_1 * dpdpx_3, DiffPair_float_0 * dpdOut_2, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_1 _S1442 = *dpdpx_3;
    float _S1443 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_11->_S1321)->primal_0);
    float _S1444 = s_primal_ctx_sqrt_0(_S1443);
    float _S1445 = 0.5f / _S1444 * (*dpdpx_3).differential_0.differential_0;
    float _S1446 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S1444 * _S1444));
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1443;
    (&_S1447)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1447, _S1446);
    DiffPair_float_0 _S1448;
    (&_S1448)->primal_0 = 1.00000001168609742e-07f;
    (&_S1448)->differential_0 = 0.0f;
    DiffPair_float_0 _S1449;
    (&_S1449)->primal_0 = (&_s_diff_ctx_11->_S1321)->primal_0;
    (&_S1449)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1448, &_S1449, _S1447.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S1449.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S1445;
    dpdpx_3->primal_0 = _S1442.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S1450, DiffPair_float_0 * _S1451, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_1 _S1452 = *_S1450;
    DiffPair_float_0 _S1453 = _s_diff_ctx_12->_S1317;
    DiffPair_float_0 _S1454 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S1455;
    (&_S1455)->_S1321 = _S1454;
    s_primal_ctx_d_sqrt_0(&_S1453, (*_S1451).primal_0, &_S1455);
    DiffPair_float_0 _S1456 = { (*_S1450).differential_0.primal_0, (*_S1450).differential_0.differential_0 };
    DiffPair_1 _S1457;
    (&_S1457)->primal_0 = _s_diff_ctx_12->_S1317;
    (&_S1457)->differential_0 = _S1456;
    DiffPair_float_0 _S1458;
    (&_S1458)->primal_0 = (*_S1451).primal_0;
    (&_S1458)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1459 = _S1455;
    s_bwd_prop_d_sqrt_0(&_S1457, &_S1458, &_S1459);
    DiffPair_float_0 _S1460 = { _S1457.differential_0.primal_0, _S1457.differential_0.differential_0 };
    _S1451->primal_0 = (*_S1451).primal_0;
    _S1451->differential_0 = _S1458.differential_0;
    _S1450->primal_0 = _S1452.primal_0;
    _S1450->differential_0 = _S1460;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S1461, float _s_dOut_3)
{
    DiffPair_float_0 _S1462;
    (&_S1462)->primal_0 = (*_S1461).primal_0;
    (&_S1462)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1462, _s_dOut_3);
    _S1461->primal_0 = (*_S1461).primal_0;
    _S1461->differential_0 = _S1462.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_0 _S1463 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->x) * *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->x) + *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->y) * *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->y);
    DiffPair_float_0 _S1464 = { len_0, 0.0f };
    float2  _S1465 = make_float2 (0.0f);
    float _S1466 = (*dpdpx_5).differential_0.differential_0.x;
    float _S1467 = _S1466 + _S1466;
    float _S1468 = (&_s_diff_ctx_13->_S1319)->differential_0 * _S1467;
    float _S1469 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S1470 = (&_s_diff_ctx_13->_S1319)->differential_0 * _S1469;
    DiffPair_float_0 _S1471 = { 0.0f, *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->x) * _S1467 + *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->y) * _S1469 };
    DiffPair_1 _S1472;
    (&_S1472)->primal_0 = _S1464;
    (&_S1472)->differential_0 = _S1471;
    DiffPair_float_0 _S1473;
    (&_S1473)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S1473)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S1472, &_S1473, &_s_diff_ctx_13->_S1320);
    DiffPair_float_0 _S1474;
    (&_S1474)->primal_0 = len_0;
    (&_S1474)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1474, 0.0f);
    float _S1475 = _S1472.differential_0.primal_0 + _S1474.differential_0;
    float _S1476 = *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->y) * _S1475;
    float _S1477 = _S1470 + _S1476 + _S1476;
    float _S1478 = *&((&(&_s_diff_ctx_13->_S1318)->primal_0)->x) * _S1475;
    float _S1479 = _S1468 + _S1478 + _S1478;
    float2  _S1480 = _S1465;
    *&((&_S1480)->y) = _S1477;
    *&((&_S1480)->x) = _S1479;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S1463.differential_0.primal_0 + _S1480, _S1465 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S1473.differential_0;
    dpdpx_5->primal_0 = _S1463.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S1481 = (*dpdpx_7).primal_0.x;
    float _S1482 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S1483;
    (&_S1483)->primal_0 = _S1481 * _S1481 + _S1482 * _S1482;
    (&_S1483)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1483, _s_dOut_4);
    float _S1484 = (*dpdpx_7).primal_0.y * _S1483.differential_0;
    float _S1485 = _S1484 + _S1484;
    float _S1486 = (*dpdpx_7).primal_0.x * _S1483.differential_0;
    float _S1487 = _S1486 + _S1486;
    float2  _S1488 = make_float2 (0.0f);
    *&((&_S1488)->y) = _S1485;
    *&((&_S1488)->x) = _S1487;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S1488;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S1489, DiffPair_float_0 * _S1490, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_0 _S1491 = *_S1489;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1492 = _s_diff_ctx_14->_S1316;
    float2  _S1493 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1494 = { _S1493, _S1493 };
    DiffPair_float_0 _S1495 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1496 = { _S1495 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1497;
    (&_S1497)->_S1318 = _S1494;
    (&_S1497)->_S1319 = _S1495;
    (&_S1497)->_S1320 = _S1496;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1492, (*_S1490).primal_0, &_S1497);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1498 = { (*_S1489).differential_0.primal_0, (*_S1489).differential_0.differential_0 };
    DiffPair_0 _S1499;
    (&_S1499)->primal_0 = _s_diff_ctx_14->_S1316;
    (&_S1499)->differential_0 = _S1498;
    DiffPair_float_0 _S1500;
    (&_S1500)->primal_0 = (*_S1490).primal_0;
    (&_S1500)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1501 = _S1497;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S1499, &_S1500, &_S1501);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1502;
    (&_S1502)->primal_0 = (&_s_diff_ctx_14->_S1316)->primal_0;
    (&_S1502)->differential_0 = _S1493;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S1502, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1503 = { _S1499.differential_0.primal_0 + _S1502.differential_0, _S1499.differential_0.differential_0 };
    _S1490->primal_0 = (*_S1490).primal_0;
    _S1490->differential_0 = _S1500.differential_0;
    _S1489->primal_0 = _S1491.primal_0;
    _S1489->differential_0 = _S1503;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_1 _S1504 = *dpdpy_1;
    DiffPair_1 _S1505 = *dpdpx_8;
    float _S1506 = - (&_s_diff_ctx_15->_S1299)->primal_0;
    float _S1507 = (&_s_diff_ctx_15->_S1300)->primal_0 * (&_s_diff_ctx_15->_S1300)->primal_0 + (&_s_diff_ctx_15->_S1299)->primal_0 * (&_s_diff_ctx_15->_S1299)->primal_0;
    float _S1508 = _S1507 * _S1507;
    float _S1509 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S1508;
    float _S1510 = (&_s_diff_ctx_15->_S1300)->primal_0 * - _S1509;
    float _S1511 = (&_s_diff_ctx_15->_S1299)->primal_0 * _S1510;
    float _S1512 = (&_s_diff_ctx_15->_S1300)->primal_0 * _S1510;
    float _S1513 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S1508;
    float _S1514 = _S1506 * - _S1513;
    float _S1515 = (&_s_diff_ctx_15->_S1299)->primal_0 * _S1514;
    float _S1516 = (&_s_diff_ctx_15->_S1300)->primal_0 * _S1514;
    DiffPair_float_0 dpdpx_9 = { _S1516 + _S1516 + ((*dpdpx_8).differential_0.primal_0 + (_S1512 + _S1512 + _S1507 * _S1509)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S1511 + _S1511 + (*dpdpy_1).differential_0.primal_0 + _S1515 + _S1515 + - (_S1507 * _S1513), 0.0f };
    float _S1517 = (&_s_diff_ctx_15->_S1300)->primal_0 / _S1507 * (*dpdpy_1).differential_0.differential_0 + _S1506 / _S1507 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S1517;
    dpdpy_1->primal_0 = _S1504.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S1505.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S1518, DiffPair_1 * _S1519, DiffPair_float_0 * _S1520, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_16)
{
    DiffPair_1 _S1521 = *_S1518;
    DiffPair_1 _S1522 = *_S1519;
    DiffPair_float_0 _S1523 = _s_diff_ctx_16->_S1287;
    DiffPair_float_0 _S1524 = _s_diff_ctx_16->_S1288;
    DiffPair_float_0 _S1525 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S1526;
    (&_S1526)->_S1299 = _S1525;
    (&_S1526)->_S1300 = _S1525;
    s_primal_ctx_d_atan2_0(&_S1523, &_S1524, (*_S1520).primal_0, &_S1526);
    DiffPair_float_0 _S1527 = { (*_S1519).differential_0.primal_0, (*_S1519).differential_0.differential_0 };
    DiffPair_float_0 _S1528 = { (*_S1518).differential_0.primal_0, (*_S1518).differential_0.differential_0 };
    DiffPair_1 _S1529;
    (&_S1529)->primal_0 = _s_diff_ctx_16->_S1287;
    (&_S1529)->differential_0 = _S1528;
    DiffPair_1 _S1530;
    (&_S1530)->primal_0 = _s_diff_ctx_16->_S1288;
    (&_S1530)->differential_0 = _S1527;
    DiffPair_float_0 _S1531;
    (&_S1531)->primal_0 = (*_S1520).primal_0;
    (&_S1531)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S1532 = _S1526;
    s_bwd_prop_d_atan2_0(&_S1529, &_S1530, &_S1531, &_S1532);
    DiffPair_float_0 _S1533 = { _S1530.differential_0.primal_0, _S1530.differential_0.differential_0 };
    DiffPair_float_0 _S1534 = { _S1529.differential_0.primal_0, _S1529.differential_0.differential_0 };
    _S1520->primal_0 = (*_S1520).primal_0;
    _S1520->differential_0 = _S1531.differential_0;
    _S1518->primal_0 = _S1521.primal_0;
    _S1518->differential_0 = _S1534;
    _S1519->primal_0 = _S1522.primal_0;
    _S1519->differential_0 = _S1533;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S1535, DiffPair_float_0 * _S1536, float _s_dOut_5)
{
    DiffPair_float_0 _S1537;
    (&_S1537)->primal_0 = (*_S1535).primal_0;
    (&_S1537)->differential_0 = 0.0f;
    DiffPair_float_0 _S1538;
    (&_S1538)->primal_0 = (*_S1536).primal_0;
    (&_S1538)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1537, &_S1538, _s_dOut_5);
    _S1536->primal_0 = (*_S1536).primal_0;
    _S1536->differential_0 = _S1538.differential_0;
    _S1535->primal_0 = (*_S1535).primal_0;
    _S1535->differential_0 = _S1537.differential_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_3, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_3, DiffPair_float_0 * dpfx_3, DiffPair_float_0 * dpfy_3, DiffPair_float_0 * dpcx_3, DiffPair_float_0 * dpcy_3, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_3, Matrix<float, 2, 2>  dpcov2d_3, float2  dpmean2d_3, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_17)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1539 = *dpcov3d_3;
    DiffPair_float_0 _S1540 = *dpfx_3;
    DiffPair_float_0 _S1541 = *dpfy_3;
    FixedArray<float, 10>  _S1542 = dpdist_coeffs_3->primal_0;
    float2  _S1543 = make_float2 (0.0f);
    float2  _S1544 = float2 {(*dpmean3d_3).primal_0.x, (*dpmean3d_3).primal_0.y};
    float _S1545 = length_0(_S1544);
    float _S1546 = (*dpmean3d_3).primal_0.z;
    float _S1547 = s_primal_ctx_atan2_0(_S1545, _S1546);
    bool _S1548 = _S1547 < 0.00100000004749745f;
    float k_2;
    float _S1549;
    float _S1550;
    float _S1551;
    if(_S1548)
    {
        float _S1552 = 1.0f - _S1547 * _S1547 / 3.0f;
        float _S1553 = _S1546 * _S1546;
        k_2 = _S1552 / _S1546;
        _S1549 = 0.0f;
        _S1550 = _S1553;
        _S1551 = _S1552;
    }
    else
    {
        float _S1554 = _S1545 * _S1545;
        k_2 = _S1547 / _S1545;
        _S1549 = _S1554;
        _S1550 = 0.0f;
        _S1551 = 0.0f;
    }
    float2  _S1555 = make_float2 (k_2);
    float2  _S1556 = _S1544 * make_float2 (k_2);
    float u_19 = _S1556.x;
    float v_19 = _S1556.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S1557 = _S1542[int(2)] + r2_19 * _S1542[int(3)];
    float _S1558 = _S1542[int(1)] + r2_19 * _S1557;
    float _S1559 = _S1542[int(0)] + r2_19 * _S1558;
    float radial_2 = 1.0f + r2_19 * _S1559;
    float2  _S1560 = make_float2 (radial_2);
    float _S1561 = 2.0f * _S1542[int(4)];
    float _S1562 = _S1561 * u_19;
    float _S1563 = 2.0f * u_19;
    float _S1564 = r2_19 + _S1563 * u_19;
    float _S1565 = 2.0f * _S1542[int(5)];
    float _S1566 = _S1565 * u_19;
    float _S1567 = 2.0f * v_19;
    float _S1568 = r2_19 + _S1567 * v_19;
    float2  _S1569 = _S1556 * make_float2 (radial_2) + make_float2 (_S1562 * v_19 + _S1542[int(5)] * _S1564 + _S1542[int(6)] * r2_19, _S1566 * v_19 + _S1542[int(4)] * _S1568 + _S1542[int(7)] * r2_19);
    float _S1570 = _S1569.x;
    float _S1571 = _S1569.y;
    float2  _S1572 = _S1569 + make_float2 (_S1542[int(8)] * _S1570 + _S1542[int(9)] * _S1571, 0.0f);
    float _S1573 = _S1572.x;
    float _S1574 = _S1572.y;
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (0.0f);
    float _S1575 = s_primal_ctx_s_primal_ctx_atan2_0(_S1545, _S1546);
    bool _S1576 = _S1575 < 0.00100000004749745f;
    float _S1577;
    float _S1578;
    float _S1579;
    if(_S1576)
    {
        float _S1580 = 1.0f - _S1575 * _S1575 / 3.0f;
        float _S1581 = _S1546 * _S1546;
        k_2 = _S1580 / _S1546;
        _S1577 = 0.0f;
        _S1578 = _S1581;
        _S1579 = _S1580;
    }
    else
    {
        float _S1582 = _S1545 * _S1545;
        k_2 = _S1575 / _S1545;
        _S1577 = _S1582;
        _S1578 = 0.0f;
        _S1579 = 0.0f;
    }
    float2  _S1583 = make_float2 (k_2);
    float2  _S1584 = _S1544 * make_float2 (k_2);
    float u_20 = _S1584.x;
    float v_20 = _S1584.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S1585 = _S1542[int(2)] + r2_20 * _S1542[int(3)];
    float _S1586 = _S1542[int(1)] + r2_20 * _S1585;
    float _S1587 = _S1542[int(0)] + r2_20 * _S1586;
    float2  _S1588 = make_float2 (1.0f + r2_20 * _S1587);
    float _S1589 = _S1561 * u_20;
    float _S1590 = 2.0f * u_20;
    float _S1591 = _S1565 * u_20;
    float _S1592 = 2.0f * v_20;
    float2  _S1593 = make_float2 (_S1540.primal_0, 0.0f) + make_float2 (_S1542[int(8)] * _S1540.primal_0, _S1542[int(9)] * _S1540.primal_0);
    float2  _S1594 = _S1584 * _S1593;
    float _S1595 = _S1542[int(4)] * _S1593.y;
    float _S1596 = v_20 * _S1593.y;
    float _S1597 = _S1542[int(5)] * _S1593.x;
    float _S1598 = v_20 * _S1593.x;
    float _S1599 = _S1594.x + _S1594.y;
    float _S1600 = r2_20 * _S1599;
    float _S1601 = r2_20 * _S1600;
    float _S1602 = r2_20 * _S1601;
    float _S1603 = _S1542[int(7)] * _S1593.y + _S1595 + _S1542[int(6)] * _S1593.x + _S1597 + _S1587 * _S1599 + _S1586 * _S1600 + _S1585 * _S1601 + _S1542[int(3)] * _S1602;
    float _S1604 = v_20 * _S1603;
    float _S1605 = u_20 * _S1603;
    float2  _S1606 = _S1588 * _S1593 + make_float2 (_S1565 * _S1596 + _S1590 * _S1597 + 2.0f * (u_20 * _S1597) + _S1561 * _S1598 + _S1605 + _S1605, _S1592 * _S1595 + 2.0f * (v_20 * _S1595) + _S1591 * _S1593.y + _S1589 * _S1593.x + _S1604 + _S1604);
    float2  _S1607 = _S1544 * _S1606;
    float2  _S1608 = _S1583 * _S1606;
    float _S1609 = _S1607.x + _S1607.y;
    float k_3;
    float _S1610;
    float _S1611;
    float _S1612;
    float _S1613;
    float _S1614;
    float _S1615;
    float _S1616;
    float _S1617;
    if(_S1576)
    {
        float _S1618 = _S1609 / _S1578;
        float _S1619 = _S1578 * _S1578;
        float _S1620 = - _S1618;
        float _S1621 = _S1579 * _S1620;
        float _S1622 = 0.3333333432674408f * - (_S1546 * _S1618);
        float _S1623 = _S1575 * _S1622;
        k_2 = _S1623 + _S1623;
        k_3 = _S1621;
        _S1610 = 0.0f;
        _S1611 = 0.0f;
        _S1612 = 0.0f;
        _S1613 = 0.0f;
        _S1614 = _S1622;
        _S1615 = _S1618;
        _S1616 = _S1620;
        _S1617 = _S1619;
    }
    else
    {
        float _S1624 = _S1609 / _S1577;
        float _S1625 = _S1577 * _S1577;
        float _S1626 = - _S1624;
        float _S1627 = _S1575 * _S1626;
        k_2 = _S1545 * _S1624;
        k_3 = 0.0f;
        _S1610 = _S1627;
        _S1611 = _S1624;
        _S1612 = _S1626;
        _S1613 = _S1625;
        _S1614 = 0.0f;
        _S1615 = 0.0f;
        _S1616 = 0.0f;
        _S1617 = 0.0f;
    }
    DiffPair_float_0 _S1628 = { _S1545, 0.0f };
    DiffPair_float_0 _S1629 = { _S1546, 0.0f };
    float _S1630 = (&_s_diff_ctx_17->_S1290)->differential_0 + k_3;
    float _S1631 = (&_s_diff_ctx_17->_S1289)->differential_0 + _S1610;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1632 = { _S1544, _S1543 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1633;
    (&_S1633)->primal_0 = _S1544;
    (&_S1633)->differential_0 = _S1543;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1634 = { _S1543, _S1543 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1635;
    (&_S1635)->_S1316 = _S1634;
    s_primal_ctx_s_bwd_length_impl_0(&_S1633, _S1631, &_S1635);
    float2  _S1636 = _S1633.differential_0 + _S1608;
    float3  _S1637 = make_float3 (_S1636.x, _S1636.y, _S1630);
    Matrix<float, 2, 3>  _S1638 = J_11;
    _S1638[int(0)] = _S1637;
    float _S1639;
    float _S1640;
    if(_S1576)
    {
        float _S1641 = 1.0f - _S1575 * _S1575 / 3.0f;
        float _S1642 = _S1546 * _S1546;
        k_3 = _S1641 / _S1546;
        _S1610 = 0.0f;
        _S1639 = _S1642;
        _S1640 = _S1641;
    }
    else
    {
        float _S1643 = _S1545 * _S1545;
        k_3 = _S1575 / _S1545;
        _S1610 = _S1643;
        _S1639 = 0.0f;
        _S1640 = 0.0f;
    }
    float2  _S1644 = make_float2 (k_3);
    float2  _S1645 = _S1544 * make_float2 (k_3);
    float u_21 = _S1645.x;
    float v_21 = _S1645.y;
    float r2_21 = u_21 * u_21 + v_21 * v_21;
    float _S1646 = _S1542[int(2)] + r2_21 * _S1542[int(3)];
    float _S1647 = _S1542[int(1)] + r2_21 * _S1646;
    float _S1648 = _S1542[int(0)] + r2_21 * _S1647;
    float2  _S1649 = make_float2 (1.0f + r2_21 * _S1648);
    float _S1650 = _S1561 * u_21;
    float _S1651 = 2.0f * u_21;
    float _S1652 = _S1565 * u_21;
    float _S1653 = 2.0f * v_21;
    float2  _S1654 = make_float2 (0.0f, _S1541.primal_0);
    float2  _S1655 = _S1645 * _S1654;
    float _S1656 = _S1542[int(4)] * _S1541.primal_0;
    float _S1657 = v_21 * _S1541.primal_0;
    float _S1658 = _S1655.x + _S1655.y;
    float _S1659 = r2_21 * _S1658;
    float _S1660 = r2_21 * _S1659;
    float _S1661 = r2_21 * _S1660;
    float _S1662 = _S1542[int(7)] * _S1541.primal_0 + _S1656 + _S1648 * _S1658 + _S1647 * _S1659 + _S1646 * _S1660 + _S1542[int(3)] * _S1661;
    float _S1663 = v_21 * _S1662;
    float _S1664 = u_21 * _S1662;
    float2  _S1665 = _S1649 * _S1654 + make_float2 (_S1565 * _S1657 + _S1664 + _S1664, _S1653 * _S1656 + 2.0f * (v_21 * _S1656) + _S1652 * _S1541.primal_0 + _S1663 + _S1663);
    float2  _S1666 = _S1544 * _S1665;
    float2  _S1667 = _S1644 * _S1665;
    float _S1668 = _S1666.x + _S1666.y;
    float _S1669;
    float _S1670;
    float _S1671;
    float _S1672;
    float _S1673;
    float _S1674;
    float _S1675;
    float _S1676;
    float _S1677;
    if(_S1576)
    {
        float _S1678 = _S1668 / _S1639;
        float _S1679 = _S1639 * _S1639;
        float _S1680 = - _S1678;
        float _S1681 = _S1640 * _S1680;
        float _S1682 = 0.3333333432674408f * - (_S1546 * _S1678);
        float _S1683 = _S1575 * _S1682;
        k_3 = _S1683 + _S1683;
        _S1669 = _S1681;
        _S1670 = 0.0f;
        _S1671 = 0.0f;
        _S1672 = 0.0f;
        _S1673 = 0.0f;
        _S1674 = _S1682;
        _S1675 = _S1678;
        _S1676 = _S1680;
        _S1677 = _S1679;
    }
    else
    {
        float _S1684 = _S1668 / _S1610;
        float _S1685 = _S1610 * _S1610;
        float _S1686 = - _S1684;
        float _S1687 = _S1575 * _S1686;
        k_3 = _S1545 * _S1684;
        _S1669 = 0.0f;
        _S1670 = _S1687;
        _S1671 = _S1684;
        _S1672 = _S1686;
        _S1673 = _S1685;
        _S1674 = 0.0f;
        _S1675 = 0.0f;
        _S1676 = 0.0f;
        _S1677 = 0.0f;
    }
    float _S1688 = (&_s_diff_ctx_17->_S1293)->differential_0 + _S1669;
    float _S1689 = (&_s_diff_ctx_17->_S1292)->differential_0 + _S1670;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1690;
    (&_S1690)->primal_0 = _S1544;
    (&_S1690)->differential_0 = _S1543;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1691;
    (&_S1691)->_S1316 = _S1634;
    s_primal_ctx_s_bwd_length_impl_0(&_S1690, _S1689, &_S1691);
    float2  _S1692 = _S1690.differential_0 + _S1667;
    float3  _S1693 = make_float3 (_S1692.x, _S1692.y, _S1688);
    _S1638[int(1)] = _S1693;
    Matrix<float, 3, 2>  _S1694 = transpose_1(_S1638);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1695;
    (&_S1695)->primal_0 = s_primal_ctx_mul_3(_S1638, _S1539.primal_0);
    (&_S1695)->differential_0 = J_11;
    Matrix<float, 3, 2>  _S1696 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1697;
    (&_S1697)->primal_0 = _S1694;
    (&_S1697)->differential_0 = _S1696;
    s_bwd_prop_mul_2(&_S1695, &_S1697, dpcov2d_3);
    Matrix<float, 2, 3>  _S1698 = transpose_2(_S1697.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1699;
    (&_S1699)->primal_0 = _S1638;
    (&_S1699)->differential_0 = J_11;
    Matrix<float, 3, 3>  _S1700 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1701;
    (&_S1701)->primal_0 = _S1539.primal_0;
    (&_S1701)->differential_0 = _S1700;
    s_bwd_prop_mul_3(&_S1699, &_S1701, _S1695.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1702 = _S1701;
    Matrix<float, 2, 3>  _S1703 = _S1698 + _S1699.differential_0;
    float2  _S1704 = _S1543;
    *&((&_S1704)->y) = _S1703.rows[int(1)].y;
    *&((&_S1704)->x) = _S1703.rows[int(1)].x;
    float2  _S1705 = _S1704;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1706 = { _S1543, _S1704 };
    DiffPair_0 _S1707;
    (&_S1707)->primal_0 = _S1632;
    (&_S1707)->differential_0 = _S1706;
    DiffPair_float_0 _S1708;
    (&_S1708)->primal_0 = _S1689;
    (&_S1708)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1709 = _S1691;
    s_bwd_prop_s_bwd_length_impl_0(&_S1707, &_S1708, &_S1709);
    DiffPair_0 _S1710 = _S1707;
    DiffPair_float_0 _S1711 = _S1708;
    DiffPair_float_0 _S1712 = { 0.0f, _S1703.rows[int(1)].z };
    DiffPair_float_0 _S1713 = { 0.0f, _S1708.differential_0 };
    DiffPair_1 _S1714;
    (&_S1714)->primal_0 = _S1628;
    (&_S1714)->differential_0 = _S1713;
    DiffPair_1 _S1715;
    (&_S1715)->primal_0 = _S1629;
    (&_S1715)->differential_0 = _S1712;
    DiffPair_float_0 _S1716;
    (&_S1716)->primal_0 = k_3;
    (&_S1716)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1714, &_S1715, &_S1716, &_s_diff_ctx_17->_S1294);
    DiffPair_1 _S1717 = _S1714;
    DiffPair_1 _S1718 = _S1715;
    DiffPair_float_0 _S1719 = _S1716;
    if(_S1576)
    {
        float _S1720 = _S1719.differential_0 + _S1719.differential_0;
        float _S1721 = _S1674 * _S1720;
        float _S1722 = - (0.3333333432674408f * (_S1575 * _S1720));
        float _S1723 = _S1676 * _S1703.rows[int(1)].z;
        float _S1724 = (_S1546 * _S1722 + - (_S1640 * _S1703.rows[int(1)].z)) / _S1677;
        float _S1725 = _S1668 * - _S1724;
        float _S1726 = _S1675 * _S1722 + _S1718.differential_0.primal_0;
        k_3 = _S1639 * _S1724;
        _S1669 = 0.0f;
        _S1670 = _S1725;
        _S1671 = _S1723;
        _S1672 = _S1717.differential_0.primal_0;
        _S1673 = _S1721;
        _S1674 = _S1726;
    }
    else
    {
        float _S1727 = _S1672 * _S1711.differential_0;
        float _S1728 = (_S1545 * _S1719.differential_0 + - (_S1575 * _S1711.differential_0)) / _S1673;
        float _S1729 = _S1668 * - _S1728;
        float _S1730 = _S1671 * _S1719.differential_0 + _S1717.differential_0.primal_0;
        k_3 = _S1610 * _S1728;
        _S1669 = _S1729;
        _S1670 = 0.0f;
        _S1671 = 0.0f;
        _S1672 = _S1730;
        _S1673 = _S1727;
        _S1674 = _S1718.differential_0.primal_0;
    }
    float2  _S1731 = _S1644 * _S1705;
    float2  _S1732 = _S1665 * _S1705;
    float2  _S1733 = _S1543;
    *&((&_S1733)->y) = k_3;
    *&((&_S1733)->x) = k_3;
    float2  _S1734 = _S1665 * _S1733;
    float2  _S1735 = _S1731 + _S1544 * _S1733;
    float _S1736 = _S1735.x;
    float _S1737 = _S1736 + _S1736;
    float _S1738 = _S1662 * _S1737;
    float _S1739 = _S1735.y + _S1735.y;
    float _S1740 = _S1662 * _S1739;
    float _S1741 = u_21 * _S1737 + v_21 * _S1739;
    float _S1742 = _S1542[int(3)] * _S1741;
    float _S1743 = _S1661 * _S1741;
    float _S1744 = _S1660 * _S1741;
    float _S1745 = _S1660 * _S1742;
    float _S1746 = _S1659 * _S1741;
    float _S1747 = _S1646 * _S1741 + r2_21 * _S1742;
    float _S1748 = _S1659 * _S1747;
    float _S1749 = _S1658 * _S1741;
    float _S1750 = _S1647 * _S1741 + r2_21 * _S1747;
    float _S1751 = _S1658 * _S1750;
    float _S1752 = _S1648 * _S1741 + r2_21 * _S1750;
    float _S1753 = v_21 * (_S1561 * _S1735.x);
    float _S1754 = _S1650 * _S1735.y;
    float _S1755 = _S1542[int(5)] * (_S1741 + u_21 * (2.0f * _S1735.x) + _S1651 * _S1735.x);
    float _S1756 = _S1542[int(6)] * _S1741;
    float _S1757 = _S1565 * _S1735.x;
    float _S1758 = _S1657 * _S1735.x;
    float _S1759 = v_21 * _S1757;
    float _S1760 = _S1541.primal_0 * _S1757;
    float _S1761 = _S1652 * _S1735.y;
    float _S1762 = _S1541.primal_0 * _S1735.y;
    float _S1763 = 2.0f * _S1735.y;
    float _S1764 = _S1656 * _S1763;
    float _S1765 = _S1656 * _S1735.y;
    float _S1766 = _S1741 + v_21 * _S1763 + _S1653 * _S1735.y;
    float _S1767 = _S1542[int(4)] * _S1766;
    float _S1768 = _S1541.primal_0 * _S1766;
    float _S1769 = _S1542[int(7)] * _S1741;
    float _S1770 = _S1541.primal_0 * _S1741;
    float2  _S1771 = _S1649 * _S1735;
    float2  _S1772 = _S1654 * _S1735;
    float2  _S1773 = _S1543;
    *&((&_S1773)->y) = _S1752;
    *&((&_S1773)->x) = _S1752;
    float2  _S1774 = _S1654 * _S1773;
    float _S1775 = _S1759 + _S1761 + _S1767 + _S1769;
    float _S1776 = _S1753 + _S1754 + _S1755 + _S1756;
    float2  _S1777 = _S1771 + _S1645 * _S1773;
    float2  _S1778 = _S1543;
    *&((&_S1778)->y) = _S1775;
    *&((&_S1778)->x) = _S1776;
    float2  _S1779 = _S1777 + _S1778;
    float _S1780 = _S1772.x + _S1772.y;
    float _S1781 = _S1749 + r2_21 * _S1780;
    float _S1782 = _S1746 + r2_21 * _S1781;
    float _S1783 = _S1744 + r2_21 * _S1782;
    float _S1784 = _S1745 + _S1748 + _S1751 + _S1648 * _S1780 + _S1647 * _S1781 + _S1646 * _S1782 + _S1542[int(3)] * _S1783;
    float _S1785 = v_21 * _S1784;
    float _S1786 = u_21 * _S1784;
    float2  _S1787 = _S1734 + _S1710.differential_0.primal_0;
    float2  _S1788 = _S1774 + make_float2 (_S1738 + _S1565 * _S1762 + _S1786 + _S1786, _S1740 + _S1760 + _S1764 + 2.0f * _S1765 + _S1785 + _S1785);
    float _S1789 = _S1758 + u_21 * _S1762;
    float _S1790 = _S1743 + r2_21 * _S1783;
    float2  _S1791 = _S1544 * _S1788;
    float _S1792 = _S1732.x + _S1732.y + _S1791.x + _S1791.y;
    float2  _S1793 = _S1644 * _S1788 + _S1787;
    if(_S1576)
    {
        float _S1794 = _S1546 * _S1670;
        float _S1795 = _S1792 / _S1639;
        float _S1796 = _S1575 * (0.3333333432674408f * - (_S1671 + _S1546 * _S1795));
        float _S1797 = _S1794 + _S1794 + _S1640 * - _S1795 + _S1674;
        k_3 = _S1796 + _S1796 + _S1673;
        _S1610 = _S1797;
        _S1639 = _S1672;
    }
    else
    {
        float _S1798 = _S1545 * _S1669;
        float _S1799 = _S1792 / _S1610;
        float _S1800 = _S1798 + _S1798 + _S1575 * - _S1799 + _S1672;
        k_3 = _S1545 * _S1799 + _S1673;
        _S1610 = _S1674;
        _S1639 = _S1800;
    }
    DiffPair_float_0 _S1801;
    (&_S1801)->primal_0 = _S1545;
    (&_S1801)->differential_0 = 0.0f;
    DiffPair_float_0 _S1802;
    (&_S1802)->primal_0 = _S1546;
    (&_S1802)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1801, &_S1802, k_3);
    float _S1803 = _S1802.differential_0 + _S1610;
    float _S1804 = _S1801.differential_0 + _S1639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1805;
    (&_S1805)->primal_0 = _S1544;
    (&_S1805)->differential_0 = _S1543;
    s_bwd_length_impl_1(&_S1805, _S1804);
    float2  _S1806 = _S1805.differential_0 + _S1793;
    float3  _S1807 = make_float3 (_S1806.x, _S1806.y, _S1803);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1808;
    (&_S1808)->primal_0 = _S1544;
    (&_S1808)->differential_0 = _S1543;
    s_bwd_length_impl_1(&_S1808, 0.0f);
    float3  _S1809 = _S1807 + make_float3 (_S1808.differential_0.x, _S1808.differential_0.y, 0.0f);
    float2  _S1810 = _S1543;
    *&((&_S1810)->y) = _S1703.rows[int(0)].y;
    *&((&_S1810)->x) = _S1703.rows[int(0)].x;
    float2  _S1811 = _S1810;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1812 = { _S1543, _S1810 };
    DiffPair_0 _S1813;
    (&_S1813)->primal_0 = _S1632;
    (&_S1813)->differential_0 = _S1812;
    DiffPair_float_0 _S1814;
    (&_S1814)->primal_0 = _S1631;
    (&_S1814)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1815 = _S1635;
    s_bwd_prop_s_bwd_length_impl_0(&_S1813, &_S1814, &_S1815);
    DiffPair_0 _S1816 = _S1813;
    DiffPair_float_0 _S1817 = _S1814;
    DiffPair_float_0 _S1818 = { 0.0f, _S1703.rows[int(0)].z };
    DiffPair_float_0 _S1819 = { 0.0f, _S1814.differential_0 };
    DiffPair_1 _S1820;
    (&_S1820)->primal_0 = _S1628;
    (&_S1820)->differential_0 = _S1819;
    DiffPair_1 _S1821;
    (&_S1821)->primal_0 = _S1629;
    (&_S1821)->differential_0 = _S1818;
    DiffPair_float_0 _S1822;
    (&_S1822)->primal_0 = k_2;
    (&_S1822)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1820, &_S1821, &_S1822, &_s_diff_ctx_17->_S1291);
    DiffPair_1 _S1823 = _S1820;
    DiffPair_1 _S1824 = _S1821;
    DiffPair_float_0 _S1825 = _S1822;
    if(_S1576)
    {
        float _S1826 = _S1825.differential_0 + _S1825.differential_0;
        float _S1827 = _S1614 * _S1826;
        float _S1828 = - (0.3333333432674408f * (_S1575 * _S1826));
        float _S1829 = _S1616 * _S1703.rows[int(0)].z;
        float _S1830 = (_S1546 * _S1828 + - (_S1579 * _S1703.rows[int(0)].z)) / _S1617;
        float _S1831 = _S1609 * - _S1830;
        float _S1832 = _S1615 * _S1828 + _S1824.differential_0.primal_0;
        k_2 = _S1578 * _S1830;
        k_3 = 0.0f;
        _S1610 = _S1831;
        _S1611 = _S1829;
        _S1612 = _S1823.differential_0.primal_0;
        _S1613 = _S1827;
        _S1614 = _S1832;
    }
    else
    {
        float _S1833 = _S1612 * _S1817.differential_0;
        float _S1834 = (_S1545 * _S1825.differential_0 + - (_S1575 * _S1817.differential_0)) / _S1613;
        float _S1835 = _S1609 * - _S1834;
        float _S1836 = _S1611 * _S1825.differential_0 + _S1823.differential_0.primal_0;
        k_2 = _S1577 * _S1834;
        k_3 = _S1835;
        _S1610 = 0.0f;
        _S1611 = 0.0f;
        _S1612 = _S1836;
        _S1613 = _S1833;
        _S1614 = _S1824.differential_0.primal_0;
    }
    float2  _S1837 = _S1583 * _S1811;
    float2  _S1838 = _S1606 * _S1811;
    float2  _S1839 = _S1543;
    *&((&_S1839)->y) = k_2;
    *&((&_S1839)->x) = k_2;
    float2  _S1840 = _S1606 * _S1839;
    float2  _S1841 = _S1837 + _S1544 * _S1839;
    float _S1842 = _S1841.x;
    float _S1843 = _S1842 + _S1842;
    float _S1844 = _S1603 * _S1843;
    float _S1845 = _S1841.y + _S1841.y;
    float _S1846 = _S1603 * _S1845;
    float _S1847 = u_20 * _S1843 + v_20 * _S1845;
    float _S1848 = _S1542[int(3)] * _S1847;
    float _S1849 = _S1602 * _S1847;
    float _S1850 = _S1601 * _S1847;
    float _S1851 = _S1601 * _S1848;
    float _S1852 = _S1600 * _S1847;
    float _S1853 = _S1585 * _S1847 + r2_20 * _S1848;
    float _S1854 = _S1600 * _S1853;
    float _S1855 = _S1599 * _S1847;
    float _S1856 = _S1586 * _S1847 + r2_20 * _S1853;
    float _S1857 = _S1599 * _S1856;
    float _S1858 = _S1587 * _S1847 + r2_20 * _S1856;
    float _S1859 = _S1561 * _S1841.x;
    float _S1860 = _S1598 * _S1841.x;
    float _S1861 = v_20 * _S1859;
    float _S1862 = _S1593.x * _S1859;
    float _S1863 = _S1589 * _S1841.y;
    float _S1864 = _S1593.x * _S1841.y;
    float _S1865 = 2.0f * _S1841.x;
    float _S1866 = _S1597 * _S1865;
    float _S1867 = _S1597 * _S1841.x;
    float _S1868 = _S1847 + u_20 * _S1865 + _S1590 * _S1841.x;
    float _S1869 = _S1542[int(5)] * _S1868;
    float _S1870 = _S1593.x * _S1868;
    float _S1871 = _S1542[int(6)] * _S1847;
    float _S1872 = _S1593.x * _S1847;
    float _S1873 = _S1565 * _S1841.x;
    float _S1874 = _S1596 * _S1841.x;
    float _S1875 = v_20 * _S1873;
    float _S1876 = _S1593.y * _S1873;
    float _S1877 = _S1591 * _S1841.y;
    float _S1878 = _S1593.y * _S1841.y;
    float _S1879 = 2.0f * _S1841.y;
    float _S1880 = _S1595 * _S1879;
    float _S1881 = _S1595 * _S1841.y;
    float _S1882 = _S1847 + v_20 * _S1879 + _S1592 * _S1841.y;
    float _S1883 = _S1542[int(4)] * _S1882;
    float _S1884 = _S1593.y * _S1882;
    float _S1885 = _S1542[int(7)] * _S1847;
    float _S1886 = _S1593.y * _S1847;
    float2  _S1887 = _S1588 * _S1841;
    float2  _S1888 = _S1593 * _S1841;
    float2  _S1889 = _S1543;
    *&((&_S1889)->y) = _S1858;
    *&((&_S1889)->x) = _S1858;
    float2  _S1890 = _S1593 * _S1889;
    float _S1891 = _S1875 + _S1877 + _S1883 + _S1885;
    float _S1892 = _S1861 + _S1863 + _S1869 + _S1871;
    float2  _S1893 = _S1887 + _S1584 * _S1889;
    float2  _S1894 = _S1543;
    *&((&_S1894)->y) = _S1891;
    *&((&_S1894)->x) = _S1892;
    float2  _S1895 = _S1893 + _S1894;
    float _S1896 = _S1540.primal_0 * _S1895.x;
    float _S1897 = _S1540.primal_0 * _S1895.y;
    float _S1898 = _S1888.x + _S1888.y;
    float _S1899 = _S1855 + r2_20 * _S1898;
    float _S1900 = _S1852 + r2_20 * _S1899;
    float _S1901 = _S1850 + r2_20 * _S1900;
    float _S1902 = _S1851 + _S1854 + _S1857 + _S1587 * _S1898 + _S1586 * _S1899 + _S1585 * _S1900 + _S1542[int(3)] * _S1901;
    float _S1903 = v_20 * _S1902;
    float _S1904 = u_20 * _S1902;
    float2  _S1905 = _S1840 + _S1816.differential_0.primal_0;
    float _S1906 = _S1886 + _S1770;
    float _S1907 = _S1860 + u_20 * _S1864;
    float _S1908 = _S1542[int(8)] * _S1895.x + _S1542[int(9)] * _S1895.y + _S1895.x;
    float _S1909 = _S1900 + _S1782;
    float _S1910 = _S1899 + _S1781;
    float2  _S1911 = _S1890 + make_float2 (_S1844 + _S1866 + _S1565 * _S1878 + 2.0f * _S1867 + _S1561 * _S1864 + _S1904 + _S1904, _S1846 + _S1862 + _S1876 + _S1880 + 2.0f * _S1881 + _S1903 + _S1903);
    float _S1912 = _S1874 + u_20 * _S1878 + _S1789;
    float _S1913 = _S1849 + r2_20 * _S1901 + _S1790;
    float _S1914 = _S1901 + _S1783;
    float _S1915 = _S1884 + _S1768;
    float2  _S1916 = _S1544 * _S1911;
    float _S1917 = _S1838.x + _S1838.y + _S1916.x + _S1916.y;
    float2  _S1918 = _S1583 * _S1911 + _S1905;
    if(_S1576)
    {
        float _S1919 = _S1546 * _S1610;
        float _S1920 = _S1917 / _S1578;
        float _S1921 = _S1575 * (0.3333333432674408f * - (_S1611 + _S1546 * _S1920));
        float _S1922 = _S1919 + _S1919 + _S1579 * - _S1920 + _S1614;
        k_2 = _S1921 + _S1921 + _S1613;
        _S1577 = _S1922;
        _S1578 = _S1612;
    }
    else
    {
        float _S1923 = _S1545 * k_3;
        float _S1924 = _S1917 / _S1577;
        float _S1925 = _S1923 + _S1923 + _S1575 * - _S1924 + _S1612;
        k_2 = _S1545 * _S1924 + _S1613;
        _S1577 = _S1614;
        _S1578 = _S1925;
    }
    DiffPair_float_0 _S1926;
    (&_S1926)->primal_0 = _S1545;
    (&_S1926)->differential_0 = 0.0f;
    DiffPair_float_0 _S1927;
    (&_S1927)->primal_0 = _S1546;
    (&_S1927)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1926, &_S1927, k_2);
    float _S1928 = _S1927.differential_0 + _S1577;
    float _S1929 = _S1926.differential_0 + _S1578;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1930;
    (&_S1930)->primal_0 = _S1544;
    (&_S1930)->differential_0 = _S1543;
    s_bwd_length_impl_1(&_S1930, _S1929);
    float2  _S1931 = _S1930.differential_0 + _S1918;
    float3  _S1932 = make_float3 (_S1931.x, _S1931.y, _S1928);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1933;
    (&_S1933)->primal_0 = _S1544;
    (&_S1933)->differential_0 = _S1543;
    s_bwd_length_impl_1(&_S1933, 0.0f);
    float _S1934 = _S1540.primal_0 * dpmean2d_3.x;
    float2  _S1935 = make_float2 (_S1934, _S1541.primal_0 * dpmean2d_3.y) + make_float2 (_S1542[int(8)] * _S1934, _S1542[int(9)] * _S1934);
    float2  _S1936 = _S1556 * _S1935;
    float2  _S1937 = _S1560 * _S1935;
    float _S1938 = _S1542[int(4)] * _S1935.y;
    float _S1939 = v_19 * _S1935.y;
    float _S1940 = _S1542[int(5)] * _S1935.x;
    float _S1941 = v_19 * _S1935.x;
    float _S1942 = _S1936.x + _S1936.y;
    float _S1943 = r2_19 * _S1942;
    float _S1944 = r2_19 * _S1943;
    float _S1945 = r2_19 * _S1944;
    float _S1946 = _S1542[int(7)] * _S1935.y + _S1938 + _S1542[int(6)] * _S1935.x + _S1940 + _S1559 * _S1942 + _S1558 * _S1943 + _S1557 * _S1944 + _S1542[int(3)] * _S1945;
    float _S1947 = v_19 * _S1946;
    float _S1948 = u_19 * _S1946;
    float _S1949 = _S1567 * _S1938 + 2.0f * (v_19 * _S1938) + _S1566 * _S1935.y + _S1562 * _S1935.x + _S1947 + _S1947;
    float _S1950 = _S1565 * _S1939 + _S1563 * _S1940 + 2.0f * (u_19 * _S1940) + _S1561 * _S1941 + _S1948 + _S1948;
    float _S1951 = _S1571 * _S1934 + _S1897;
    float _S1952 = _S1570 * _S1934 + _S1896;
    float _S1953 = r2_19 * _S1935.y + _S1906;
    float _S1954 = r2_19 * _S1935.x + _S1872;
    float _S1955 = 2.0f * (u_19 * _S1939 + _S1912) + _S1564 * _S1935.x + _S1870;
    float _S1956 = _S1568 * _S1935.y + 2.0f * (u_19 * _S1941 + _S1907) + _S1915;
    float _S1957 = r2_19 * _S1945 + _S1913;
    float _S1958 = _S1945 + _S1914;
    float _S1959 = _S1944 + _S1909;
    float _S1960 = _S1943 + _S1910;
    float3  _S1961 = _S1932 + make_float3 (_S1933.differential_0.x, _S1933.differential_0.y, 0.0f) + _S1809;
    float _S1962 = _S1573 * dpmean2d_3.x + _S1908;
    float _S1963 = _S1574 * dpmean2d_3.y + _S1779.y;
    FixedArray<float, 10>  _S1964;
    _S1964[int(0)] = 0.0f;
    _S1964[int(1)] = 0.0f;
    _S1964[int(2)] = 0.0f;
    _S1964[int(3)] = 0.0f;
    _S1964[int(4)] = 0.0f;
    _S1964[int(5)] = 0.0f;
    _S1964[int(6)] = 0.0f;
    _S1964[int(7)] = 0.0f;
    _S1964[int(8)] = 0.0f;
    _S1964[int(9)] = 0.0f;
    _S1964[int(9)] = _S1951;
    _S1964[int(8)] = _S1952;
    _S1964[int(7)] = _S1953;
    _S1964[int(6)] = _S1954;
    _S1964[int(5)] = _S1955;
    _S1964[int(4)] = _S1956;
    _S1964[int(3)] = _S1957;
    _S1964[int(2)] = _S1958;
    _S1964[int(1)] = _S1959;
    _S1964[int(0)] = _S1960;
    FixedArray<float, 10>  _S1965 = {
        _S1964[int(0)], _S1964[int(1)], _S1964[int(2)], _S1964[int(3)], _S1964[int(4)], _S1964[int(5)], _S1964[int(6)], _S1964[int(7)], _S1964[int(8)], _S1964[int(9)]
    };
    float2  _S1966 = _S1937 + make_float2 (_S1950, _S1949);
    float2  _S1967 = _S1544 * _S1966;
    float2  _S1968 = _S1555 * _S1966;
    float _S1969 = _S1967.x + _S1967.y;
    if(_S1548)
    {
        float _S1970 = _S1969 / _S1550;
        float _S1971 = _S1551 * - _S1970;
        float _S1972 = _S1547 * (0.3333333432674408f * - (_S1546 * _S1970));
        k_2 = _S1972 + _S1972;
        _S1549 = _S1971;
        _S1550 = 0.0f;
    }
    else
    {
        float _S1973 = _S1969 / _S1549;
        float _S1974 = _S1547 * - _S1973;
        k_2 = _S1545 * _S1973;
        _S1549 = 0.0f;
        _S1550 = _S1974;
    }
    DiffPair_float_0 _S1975;
    (&_S1975)->primal_0 = _S1545;
    (&_S1975)->differential_0 = 0.0f;
    DiffPair_float_0 _S1976;
    (&_S1976)->primal_0 = _S1546;
    (&_S1976)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1975, &_S1976, k_2);
    float _S1977 = _S1976.differential_0 + _S1549;
    float _S1978 = _S1975.differential_0 + _S1550;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1979;
    (&_S1979)->primal_0 = _S1544;
    (&_S1979)->differential_0 = _S1543;
    s_bwd_length_impl_1(&_S1979, _S1978);
    float2  _S1980 = _S1979.differential_0 + _S1968;
    dpdist_coeffs_3->primal_0 = dpdist_coeffs_3->primal_0;
    dpdist_coeffs_3->differential_0 = _S1965;
    dpcy_3->primal_0 = (*dpcy_3).primal_0;
    dpcy_3->differential_0 = dpmean2d_3.y;
    dpcx_3->primal_0 = (*dpcx_3).primal_0;
    dpcx_3->differential_0 = dpmean2d_3.x;
    dpfy_3->primal_0 = (*dpfy_3).primal_0;
    dpfy_3->differential_0 = _S1963;
    dpfx_3->primal_0 = (*dpfx_3).primal_0;
    dpfx_3->differential_0 = _S1962;
    dpcov3d_3->primal_0 = (*dpcov3d_3).primal_0;
    dpcov3d_3->differential_0 = _S1702.differential_0;
    float3  _S1981 = _S1961 + make_float3 (_S1980.x, _S1980.y, _S1977);
    dpmean3d_3->primal_0 = (*dpmean3d_3).primal_0;
    dpmean3d_3->differential_0 = _S1981;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_15, float fy_15, float cx_15, float cy_15, FixedArray<float, 10>  * dist_coeffs_21, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1982 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1983 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1984 = { _S1983, _S1983 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1985 = { _S1983, _S1983, _S1984, _S1983, _S1983, _S1984 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1986;
    (&_S1986)->_S1295 = _S1982;
    (&_S1986)->_S1296 = _S1985;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float3  _S1987 = s_primal_ctx_exp_0(scale_13);
    float _S1988 = quat_14.y;
    float x2_14 = _S1988 * _S1988;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1989 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1987.x, 0.0f, 0.0f, 0.0f, _S1987.y, 0.0f, 0.0f, 0.0f, _S1987.z);
    Matrix<float, 3, 3>  _S1990 = s_primal_ctx_mul_2(_S1989, S_1);
    Matrix<float, 3, 3>  _S1991 = transpose_0(_S1990);
    Matrix<float, 3, 3>  _S1992 = s_primal_ctx_mul_2(_S1990, _S1991);
    Matrix<float, 3, 3>  _S1993 = s_primal_ctx_mul_2(R_15, _S1992);
    Matrix<float, 3, 3>  _S1994 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1995 = s_primal_ctx_mul_2(_S1993, _S1994);
    Matrix<float, 2, 2>  _S1996 = _S1982;
    float2  _S1997 = make_float2 (0.0f);
    float2  _S1998 = _S1997;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1995, fx_15, fy_15, cx_15, cy_15, dist_coeffs_21, &_S1996, &_S1998, &(&_S1986)->_S1296);
    (&_S1986)->_S1295 = _S1996;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1999 = _S1986;
    float _S2000 = _S1986._S1295.rows[int(0)].y * _S1986._S1295.rows[int(1)].x;
    float det_orig_12 = _S1986._S1295.rows[int(0)].x * _S1986._S1295.rows[int(1)].y - _S2000;
    float _S2001 = _S1986._S1295.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2002 = _S1986._S1295;
    *&(((&_S2002)->rows + (int(0)))->x) = _S2001;
    float _S2003 = _S1986._S1295.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2002)->rows + (int(1)))->y) = _S2003;
    Matrix<float, 2, 2>  _S2004 = _S2002;
    Matrix<float, 2, 2>  _S2005 = _S2002;
    float det_blur_7 = _S2001 * _S2003 - _S2000;
    float _S2006 = det_orig_12 / det_blur_7;
    float _S2007 = det_blur_7 * det_blur_7;
    float _S2008 = s_primal_ctx_max_0(0.0f, _S2006);
    float _S2009 = s_primal_ctx_sqrt_0(_S2008);
    float invdet_8 = 1.0f / det_blur_7;
    float _S2010 = - _S1986._S1295.rows[int(0)].y;
    float _S2011 = - _S1986._S1295.rows[int(1)].x;
    float _S2012 = - in_opacity_11;
    float _S2013 = 1.0f + s_primal_ctx_exp_1(_S2012);
    float _S2014 = 1.0f / _S2013;
    float _S2015 = _S2013 * _S2013;
    float _S2016;
    if(antialiased_11)
    {
        _S2016 = _S2014 * _S2009;
    }
    else
    {
        _S2016 = _S2014;
    }
    float _S2017 = _S2016 / 0.00392156885936856f;
    float _S2018 = 2.0f * s_primal_ctx_log_0(_S2017);
    float _S2019 = s_primal_ctx_sqrt_0(_S2018);
    float _S2020 = _S2004.rows[int(0)].x;
    float _S2021 = _S2005.rows[int(1)].y;
    float _S2022 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S2023 = mean_11 - - s_primal_ctx_mul_1(_S1994, t_14);
    float _S2024 = _S2023.x;
    float _S2025 = _S2023.y;
    float _S2026 = _S2023.z;
    float _S2027 = _S2024 * _S2024 + _S2025 * _S2025 + _S2026 * _S2026;
    float _S2028 = s_primal_ctx_sqrt_0(_S2027);
    float x_35 = _S2024 / _S2028;
    float3  _S2029 = make_float3 (x_35);
    float _S2030 = _S2028 * _S2028;
    float y_14 = _S2025 / _S2028;
    float z_11 = _S2026 / _S2028;
    float3  _S2031 = make_float3 (z_11);
    float _S2032 = - y_14;
    float3  _S2033 = make_float3 (_S2032);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_35 * x_35 - y_14 * y_14;
    float _S2034 = 2.0f * x_35;
    float fS1_11 = _S2034 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S2035 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_35;
    float3  _S2036 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S2037 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S2038 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S2039 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S2040 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S2040;
    float3  _S2041 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_35;
    float3  _S2042 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S2043 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S2044 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S2045 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_35 * fC1_11 - y_14 * fS1_11);
    float3  _S2046 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_35 * fS1_11 + y_14 * fC1_11);
    float3  _S2047 = make_float3 (pSH9_1);
    float3  _S2048 = make_float3 (0.0f);
    float3  _S2049 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2050;
    (&_S2050)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2032) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_35) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S2050)->differential_0 = _S2049;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2051;
    (&_S2051)->primal_0 = _S2048;
    (&_S2051)->differential_0 = _S2049;
    s_bwd_prop_max_0(&_S2050, &_S2051, v_rgb_1);
    float3  _S2052 = _S2046 * _S2050.differential_0;
    float3  _S2053 = (*sh_coeffs_11)[int(15)] * _S2050.differential_0;
    float3  _S2054 = _S2044 * _S2050.differential_0;
    float3  _S2055 = (*sh_coeffs_11)[int(14)] * _S2050.differential_0;
    float3  _S2056 = _S2042 * _S2050.differential_0;
    float3  _S2057 = (*sh_coeffs_11)[int(13)] * _S2050.differential_0;
    float3  _S2058 = _S2041 * _S2050.differential_0;
    float3  _S2059 = (*sh_coeffs_11)[int(12)] * _S2050.differential_0;
    float3  _S2060 = _S2043 * _S2050.differential_0;
    float3  _S2061 = (*sh_coeffs_11)[int(11)] * _S2050.differential_0;
    float3  _S2062 = _S2045 * _S2050.differential_0;
    float3  _S2063 = (*sh_coeffs_11)[int(10)] * _S2050.differential_0;
    float3  _S2064 = _S2047 * _S2050.differential_0;
    float3  _S2065 = (*sh_coeffs_11)[int(9)] * _S2050.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S2065.x + _S2065.y + _S2065.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S2053.x + _S2053.y + _S2053.z);
    float _S2066 = _S2063.x + _S2063.y + _S2063.z;
    float _S2067 = _S2055.x + _S2055.y + _S2055.z;
    float _S2068 = _S2061.x + _S2061.y + _S2061.z;
    float _S2069 = _S2057.x + _S2057.y + _S2057.z;
    float _S2070 = _S2059.x + _S2059.y + _S2059.z;
    float _S2071 = - s_diff_fC2_T_1;
    float3  _S2072 = _S2038 * _S2050.differential_0;
    float3  _S2073 = (*sh_coeffs_11)[int(8)] * _S2050.differential_0;
    float3  _S2074 = _S2036 * _S2050.differential_0;
    float3  _S2075 = (*sh_coeffs_11)[int(7)] * _S2050.differential_0;
    float3  _S2076 = _S2035 * _S2050.differential_0;
    float3  _S2077 = (*sh_coeffs_11)[int(6)] * _S2050.differential_0;
    float3  _S2078 = _S2037 * _S2050.differential_0;
    float3  _S2079 = (*sh_coeffs_11)[int(5)] * _S2050.differential_0;
    float3  _S2080 = _S2039 * _S2050.differential_0;
    float3  _S2081 = (*sh_coeffs_11)[int(4)] * _S2050.differential_0;
    float _S2082 = _S2079.x + _S2079.y + _S2079.z;
    float _S2083 = _S2075.x + _S2075.y + _S2075.z;
    float _S2084 = fTmp1B_11 * _S2066 + x_35 * s_diff_fS2_T_1 + y_14 * _S2071 + 0.54627424478530884f * (_S2081.x + _S2081.y + _S2081.z);
    float _S2085 = fTmp1B_11 * _S2067 + y_14 * s_diff_fS2_T_1 + x_35 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S2073.x + _S2073.y + _S2073.z);
    float _S2086 = y_14 * - _S2085;
    float _S2087 = x_35 * _S2085;
    float _S2088 = z_11 * (1.86588168144226074f * (z_11 * _S2070) + -2.28522896766662598f * (y_14 * _S2068 + x_35 * _S2069) + 0.94617468118667603f * (_S2077.x + _S2077.y + _S2077.z));
    float3  _S2089 = make_float3 (0.48860251903533936f) * _S2050.differential_0;
    float3  _S2090 = - _S2089;
    float3  _S2091 = _S2029 * _S2090;
    float3  _S2092 = (*sh_coeffs_11)[int(3)] * _S2090;
    float3  _S2093 = _S2031 * _S2089;
    float3  _S2094 = (*sh_coeffs_11)[int(2)] * _S2089;
    float3  _S2095 = _S2033 * _S2089;
    float3  _S2096 = (*sh_coeffs_11)[int(1)] * _S2089;
    float _S2097 = (_S2040 * _S2070 + 1.44530570507049561f * (fS1_11 * _S2066 + fC1_11 * _S2067) + -1.09254848957061768f * (y_14 * _S2082 + x_35 * _S2083) + _S2088 + _S2088 + _S2094.x + _S2094.y + _S2094.z) / _S2030;
    float _S2098 = _S2028 * _S2097;
    float _S2099 = (fTmp0C_11 * _S2068 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S2071 + fTmp0B_11 * _S2082 + _S2034 * _S2084 + _S2086 + _S2086 + - (_S2096.x + _S2096.y + _S2096.z)) / _S2030;
    float _S2100 = _S2028 * _S2099;
    float _S2101 = (fTmp0C_11 * _S2069 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S2083 + 2.0f * (y_14 * _S2084) + _S2087 + _S2087 + _S2092.x + _S2092.y + _S2092.z) / _S2030;
    float _S2102 = _S2028 * _S2101;
    float _S2103 = _S2026 * - _S2097 + _S2025 * - _S2099 + _S2024 * - _S2101;
    DiffPair_float_0 _S2104;
    (&_S2104)->primal_0 = _S2027;
    (&_S2104)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2104, _S2103);
    float _S2105 = _S2026 * _S2104.differential_0;
    float _S2106 = _S2025 * _S2104.differential_0;
    float _S2107 = _S2024 * _S2104.differential_0;
    float3  _S2108 = make_float3 (0.282094806432724f) * _S2050.differential_0;
    float3  _S2109 = make_float3 (_S2102 + _S2107 + _S2107, _S2100 + _S2106 + _S2106, _S2098 + _S2105 + _S2105);
    float3  _S2110 = - - _S2109;
    Matrix<float, 3, 3>  _S2111 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2112;
    (&_S2112)->primal_0 = _S1994;
    (&_S2112)->differential_0 = _S2111;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2113;
    (&_S2113)->primal_0 = t_14;
    (&_S2113)->differential_0 = _S2049;
    s_bwd_prop_mul_1(&_S2112, &_S2113, _S2110);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2114 = _S2112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2115 = _S2113;
    float2  _S2116 = _S1997;
    *&((&_S2116)->y) = v_conic_1.z;
    float2  _S2117 = _S1997;
    *&((&_S2117)->y) = v_conic_1.y;
    *&((&_S2117)->x) = v_conic_1.x;
    float _S2118 = 0.5f * v_depth_1;
    DiffPair_float_0 _S2119;
    (&_S2119)->primal_0 = _S2022;
    (&_S2119)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2119, _S2118);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2120;
    (&_S2120)->primal_0 = mean_c_11;
    (&_S2120)->differential_0 = _S2049;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2121;
    (&_S2121)->primal_0 = mean_c_11;
    (&_S2121)->differential_0 = _S2049;
    s_bwd_prop_dot_0(&_S2120, &_S2121, _S2119.differential_0);
    DiffPair_float_0 _S2122;
    (&_S2122)->primal_0 = _S2021;
    (&_S2122)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2122, 0.0f);
    DiffPair_float_0 _S2123;
    (&_S2123)->primal_0 = _S2020;
    (&_S2123)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2123, 0.0f);
    DiffPair_float_0 _S2124;
    (&_S2124)->primal_0 = 3.32999992370605469f;
    (&_S2124)->differential_0 = 0.0f;
    DiffPair_float_0 _S2125;
    (&_S2125)->primal_0 = _S2019;
    (&_S2125)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2124, &_S2125, 0.0f);
    DiffPair_float_0 _S2126;
    (&_S2126)->primal_0 = _S2018;
    (&_S2126)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2126, _S2125.differential_0);
    float _S2127 = 2.0f * _S2126.differential_0;
    DiffPair_float_0 _S2128;
    (&_S2128)->primal_0 = _S2017;
    (&_S2128)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2128, _S2127);
    float2  _S2129 = make_float2 (_S2123.differential_0, 0.0f);
    float _S2130 = v_opacity_1 + 254.9999847412109375f * _S2128.differential_0;
    Matrix<float, 2, 2>  _S2131 = _S1982;
    _S2131[int(1)] = _S2116;
    _S2131[int(0)] = _S2117;
    Matrix<float, 2, 2>  _S2132 = _S2131;
    FixedArray<float3 , 16>  _S2133;
    _S2133[int(0)] = _S2049;
    _S2133[int(1)] = _S2049;
    _S2133[int(2)] = _S2049;
    _S2133[int(3)] = _S2049;
    _S2133[int(4)] = _S2049;
    _S2133[int(5)] = _S2049;
    _S2133[int(6)] = _S2049;
    _S2133[int(7)] = _S2049;
    _S2133[int(8)] = _S2049;
    _S2133[int(9)] = _S2049;
    _S2133[int(10)] = _S2049;
    _S2133[int(11)] = _S2049;
    _S2133[int(12)] = _S2049;
    _S2133[int(13)] = _S2049;
    _S2133[int(14)] = _S2049;
    _S2133[int(15)] = _S2049;
    _S2133[int(7)] = _S2074;
    _S2133[int(0)] = _S2108;
    _S2133[int(1)] = _S2095;
    _S2133[int(2)] = _S2093;
    _S2133[int(3)] = _S2091;
    _S2133[int(4)] = _S2080;
    _S2133[int(5)] = _S2078;
    _S2133[int(6)] = _S2076;
    _S2133[int(15)] = _S2052;
    _S2133[int(8)] = _S2072;
    _S2133[int(9)] = _S2064;
    _S2133[int(10)] = _S2062;
    _S2133[int(11)] = _S2060;
    _S2133[int(12)] = _S2058;
    _S2133[int(13)] = _S2056;
    _S2133[int(14)] = _S2054;
    float3  _S2134 = _S2133[int(0)];
    float3  _S2135 = _S2133[int(1)];
    float3  _S2136 = _S2133[int(2)];
    float3  _S2137 = _S2133[int(3)];
    float3  _S2138 = _S2133[int(4)];
    float3  _S2139 = _S2133[int(5)];
    float3  _S2140 = _S2133[int(6)];
    float3  _S2141 = _S2133[int(7)];
    float3  _S2142 = _S2133[int(8)];
    float3  _S2143 = _S2133[int(9)];
    float3  _S2144 = _S2133[int(10)];
    float3  _S2145 = _S2133[int(11)];
    float3  _S2146 = _S2133[int(12)];
    float3  _S2147 = _S2133[int(13)];
    float3  _S2148 = _S2133[int(14)];
    float3  _S2149 = _S2133[int(15)];
    float3  _S2150 = _S2121.differential_0 + _S2120.differential_0;
    float2  _S2151 = make_float2 (0.0f, _S2122.differential_0);
    float _S2152;
    if(antialiased_11)
    {
        float _S2153 = _S2014 * _S2130;
        _S2016 = _S2009 * _S2130;
        _S2152 = _S2153;
    }
    else
    {
        _S2016 = _S2130;
        _S2152 = 0.0f;
    }
    float _S2154 = - (_S2016 / _S2015);
    DiffPair_float_0 _S2155;
    (&_S2155)->primal_0 = _S2012;
    (&_S2155)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2155, _S2154);
    float _S2156 = - _S2155.differential_0;
    float _S2157 = invdet_8 * _S2132.rows[int(1)].y;
    float _S2158 = - (invdet_8 * _S2132.rows[int(1)].x);
    float _S2159 = - (invdet_8 * _S2132.rows[int(0)].y);
    float _S2160 = invdet_8 * _S2132.rows[int(0)].x;
    float _S2161 = - ((_S2001 * _S2132.rows[int(1)].y + _S2011 * _S2132.rows[int(1)].x + _S2010 * _S2132.rows[int(0)].y + _S2003 * _S2132.rows[int(0)].x) / _S2007);
    DiffPair_float_0 _S2162;
    (&_S2162)->primal_0 = _S2008;
    (&_S2162)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2162, _S2152);
    DiffPair_float_0 _S2163;
    (&_S2163)->primal_0 = 0.0f;
    (&_S2163)->differential_0 = 0.0f;
    DiffPair_float_0 _S2164;
    (&_S2164)->primal_0 = _S2006;
    (&_S2164)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2163, &_S2164, _S2162.differential_0);
    float _S2165 = _S2164.differential_0 / _S2007;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2165;
    float _S2166 = _S2161 + det_orig_12 * - _S2165;
    float _S2167 = - _S2166;
    float _S2168 = _S2001 * _S2166;
    float _S2169 = _S2003 * _S2166;
    Matrix<float, 2, 2>  _S2170 = _S1982;
    _S2170[int(1)] = _S2151;
    _S2170[int(0)] = _S2129;
    _S2002 = _S2170;
    *&(((&_S2002)->rows + (int(1)))->y) = 0.0f;
    float _S2171 = _S2160 + _S2168 + _S2170.rows[int(1)].y;
    *&(((&_S2002)->rows + (int(0)))->x) = 0.0f;
    float _S2172 = _S2157 + _S2169 + _S2170.rows[int(0)].x;
    float _S2173 = _S2167 + - s_diff_det_orig_T_1;
    float _S2174 = _S2158 + _S1999._S1295.rows[int(0)].y * _S2173;
    float _S2175 = _S2159 + _S1999._S1295.rows[int(1)].x * _S2173;
    float _S2176 = _S1999._S1295.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2177 = _S2171 + _S1999._S1295.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2178 = _S1997;
    *&((&_S2178)->x) = _S2174;
    *&((&_S2178)->y) = _S2177;
    float _S2179 = _S2172 + _S2176;
    float2  _S2180 = _S1997;
    *&((&_S2180)->y) = _S2175;
    *&((&_S2180)->x) = _S2179;
    Matrix<float, 2, 2>  _S2181 = _S1982;
    _S2181[int(1)] = _S2178;
    _S2181[int(0)] = _S2180;
    Matrix<float, 2, 2>  _S2182 = _S2002 + _S2181;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2183;
    (&_S2183)->primal_0 = mean_c_11;
    (&_S2183)->differential_0 = _S2049;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2184;
    (&_S2184)->primal_0 = _S1995;
    (&_S2184)->differential_0 = _S2111;
    DiffPair_float_0 _S2185;
    (&_S2185)->primal_0 = fx_15;
    (&_S2185)->differential_0 = 0.0f;
    DiffPair_float_0 _S2186;
    (&_S2186)->primal_0 = fy_15;
    (&_S2186)->differential_0 = 0.0f;
    DiffPair_float_0 _S2187;
    (&_S2187)->primal_0 = cx_15;
    (&_S2187)->differential_0 = 0.0f;
    DiffPair_float_0 _S2188;
    (&_S2188)->primal_0 = cy_15;
    (&_S2188)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S2189 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2190;
    (&_S2190)->primal_0 = *dist_coeffs_21;
    (&_S2190)->differential_0 = _S2189;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2191 = _S1999._S1296;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2183, &_S2184, &_S2185, &_S2186, &_S2187, &_S2188, &_S2190, _S2182, v_mean2d_1, &_S2191);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2192;
    (&_S2192)->primal_0 = _S1993;
    (&_S2192)->differential_0 = _S2111;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2193;
    (&_S2193)->primal_0 = _S1994;
    (&_S2193)->differential_0 = _S2111;
    s_bwd_prop_mul_4(&_S2192, &_S2193, _S2184.differential_0);
    Matrix<float, 3, 3>  _S2194 = transpose_0(_S2193.differential_0 + _S2114.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2195;
    (&_S2195)->primal_0 = R_15;
    (&_S2195)->differential_0 = _S2111;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2196;
    (&_S2196)->primal_0 = _S1992;
    (&_S2196)->differential_0 = _S2111;
    s_bwd_prop_mul_4(&_S2195, &_S2196, _S2192.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2197;
    (&_S2197)->primal_0 = _S1990;
    (&_S2197)->differential_0 = _S2111;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2198;
    (&_S2198)->primal_0 = _S1991;
    (&_S2198)->differential_0 = _S2111;
    s_bwd_prop_mul_4(&_S2197, &_S2198, _S2196.differential_0);
    Matrix<float, 3, 3>  _S2199 = _S2197.differential_0 + transpose_0(_S2198.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2200;
    (&_S2200)->primal_0 = _S1989;
    (&_S2200)->differential_0 = _S2111;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2201;
    (&_S2201)->primal_0 = S_1;
    (&_S2201)->differential_0 = _S2111;
    s_bwd_prop_mul_4(&_S2200, &_S2201, _S2199);
    Matrix<float, 3, 3>  _S2202 = transpose_0(_S2200.differential_0);
    float _S2203 = 2.0f * - _S2202.rows[int(2)].z;
    float _S2204 = 2.0f * _S2202.rows[int(2)].y;
    float _S2205 = 2.0f * _S2202.rows[int(2)].x;
    float _S2206 = 2.0f * _S2202.rows[int(1)].z;
    float _S2207 = 2.0f * - _S2202.rows[int(1)].y;
    float _S2208 = 2.0f * _S2202.rows[int(1)].x;
    float _S2209 = 2.0f * _S2202.rows[int(0)].z;
    float _S2210 = 2.0f * _S2202.rows[int(0)].y;
    float _S2211 = 2.0f * - _S2202.rows[int(0)].x;
    float _S2212 = - _S2208 + _S2210;
    float _S2213 = _S2205 + - _S2209;
    float _S2214 = - _S2204 + _S2206;
    float _S2215 = _S2204 + _S2206;
    float _S2216 = _S2205 + _S2209;
    float _S2217 = _S2208 + _S2210;
    float _S2218 = quat_14.w * (_S2207 + _S2211);
    float _S2219 = quat_14.z * (_S2203 + _S2211);
    float _S2220 = quat_14.y * (_S2203 + _S2207);
    float _S2221 = quat_14.x * _S2212 + quat_14.z * _S2215 + quat_14.y * _S2216 + _S2218 + _S2218;
    float _S2222 = quat_14.x * _S2213 + quat_14.w * _S2215 + quat_14.y * _S2217 + _S2219 + _S2219;
    float _S2223 = quat_14.x * _S2214 + quat_14.w * _S2216 + quat_14.z * _S2217 + _S2220 + _S2220;
    float _S2224 = quat_14.w * _S2212 + quat_14.z * _S2213 + quat_14.y * _S2214;
    float3  _S2225 = _S2049;
    *&((&_S2225)->z) = _S2201.differential_0.rows[int(2)].z;
    *&((&_S2225)->y) = _S2201.differential_0.rows[int(1)].y;
    *&((&_S2225)->x) = _S2201.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2226;
    (&_S2226)->primal_0 = scale_13;
    (&_S2226)->differential_0 = _S2049;
    s_bwd_prop_exp_1(&_S2226, _S2225);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2227 = _S2226;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2228;
    (&_S2228)->primal_0 = mean_c_11;
    (&_S2228)->differential_0 = _S2049;
    s_bwd_length_impl_0(&_S2228, 0.0f);
    float3  _S2229 = _S2183.differential_0 + _S2228.differential_0 + _S2150;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2230;
    (&_S2230)->primal_0 = R_15;
    (&_S2230)->differential_0 = _S2111;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2231;
    (&_S2231)->primal_0 = mean_11;
    (&_S2231)->differential_0 = _S2049;
    s_bwd_prop_mul_1(&_S2230, &_S2231, _S2229);
    float3  _S2232 = _S2229 + _S2115.differential_0;
    Matrix<float, 3, 3>  _S2233 = _S2194 + _S2195.differential_0 + _S2230.differential_0;
    float4  _S2234 = make_float4 (0.0f);
    *&((&_S2234)->w) = _S2221;
    *&((&_S2234)->z) = _S2222;
    *&((&_S2234)->y) = _S2223;
    *&((&_S2234)->x) = _S2224;
    float4  _S2235 = _S2234;
    float3  _S2236 = _S2231.differential_0 + _S2109;
    *v_mean_1 = _S2236;
    *v_quat_1 = _S2235;
    *v_scale_1 = _S2227.differential_0;
    *v_in_opacity_1 = _S2156;
    (*v_sh_coeffs_1)[int(0)] = _S2134;
    (*v_sh_coeffs_1)[int(1)] = _S2135;
    (*v_sh_coeffs_1)[int(2)] = _S2136;
    (*v_sh_coeffs_1)[int(3)] = _S2137;
    (*v_sh_coeffs_1)[int(4)] = _S2138;
    (*v_sh_coeffs_1)[int(5)] = _S2139;
    (*v_sh_coeffs_1)[int(6)] = _S2140;
    (*v_sh_coeffs_1)[int(7)] = _S2141;
    (*v_sh_coeffs_1)[int(8)] = _S2142;
    (*v_sh_coeffs_1)[int(9)] = _S2143;
    (*v_sh_coeffs_1)[int(10)] = _S2144;
    (*v_sh_coeffs_1)[int(11)] = _S2145;
    (*v_sh_coeffs_1)[int(12)] = _S2146;
    (*v_sh_coeffs_1)[int(13)] = _S2147;
    (*v_sh_coeffs_1)[int(14)] = _S2148;
    (*v_sh_coeffs_1)[int(15)] = _S2149;
    *v_R_2 = _S2233;
    *v_t_2 = _S2232;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_16, float fy_16, float cx_16, float cy_16, FixedArray<float, 10>  * dist_coeffs_22, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S2237 = s_primal_ctx_exp_0(scale_14);
    float _S2238 = quat_15.y;
    float x2_15 = _S2238 * _S2238;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S2239 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2237.x, 0.0f, 0.0f, 0.0f, _S2237.y, 0.0f, 0.0f, 0.0f, _S2237.z);
    Matrix<float, 3, 3>  _S2240 = s_primal_ctx_mul_2(_S2239, S_2);
    Matrix<float, 3, 3>  _S2241 = transpose_0(_S2240);
    Matrix<float, 3, 3>  _S2242 = s_primal_ctx_mul_2(_S2240, _S2241);
    Matrix<float, 3, 3>  _S2243 = s_primal_ctx_mul_2(R_16, _S2242);
    Matrix<float, 3, 3>  _S2244 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S2245 = s_primal_ctx_mul_2(_S2243, _S2244);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (fx_16, 0.0f, 0.0f, 0.0f, fy_16, 0.0f);
    Matrix<float, 2, 3>  _S2246 = s_primal_ctx_mul_3(J_12, _S2245);
    Matrix<float, 3, 2>  _S2247 = transpose_1(J_12);
    Matrix<float, 2, 2>  _S2248 = s_primal_ctx_mul_4(_S2246, _S2247);
    float _S2249 = _S2248.rows[int(0)].y * _S2248.rows[int(1)].x;
    float det_orig_13 = _S2248.rows[int(0)].x * _S2248.rows[int(1)].y - _S2249;
    float _S2250 = _S2248.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2251 = _S2248;
    *&(((&_S2251)->rows + (int(0)))->x) = _S2250;
    float _S2252 = _S2248.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2251)->rows + (int(1)))->y) = _S2252;
    Matrix<float, 2, 2>  _S2253 = _S2251;
    Matrix<float, 2, 2>  _S2254 = _S2251;
    float det_blur_8 = _S2250 * _S2252 - _S2249;
    float _S2255 = det_orig_13 / det_blur_8;
    float _S2256 = det_blur_8 * det_blur_8;
    float _S2257 = s_primal_ctx_max_0(0.0f, _S2255);
    float _S2258 = s_primal_ctx_sqrt_0(_S2257);
    float invdet_9 = 1.0f / det_blur_8;
    float _S2259 = - _S2248.rows[int(0)].y;
    float _S2260 = - _S2248.rows[int(1)].x;
    float _S2261 = - in_opacity_12;
    float _S2262 = 1.0f + s_primal_ctx_exp_1(_S2261);
    float _S2263 = 1.0f / _S2262;
    float _S2264 = _S2262 * _S2262;
    float _S2265;
    if(antialiased_12)
    {
        _S2265 = _S2263 * _S2258;
    }
    else
    {
        _S2265 = _S2263;
    }
    float _S2266 = _S2265 / 0.00392156885936856f;
    float _S2267 = 2.0f * s_primal_ctx_log_0(_S2266);
    float _S2268 = s_primal_ctx_sqrt_0(_S2267);
    float _S2269 = _S2253.rows[int(0)].x;
    float _S2270 = _S2254.rows[int(1)].y;
    float _S2271 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S2272 = mean_12 - - s_primal_ctx_mul_1(_S2244, t_15);
    float _S2273 = _S2272.x;
    float _S2274 = _S2272.y;
    float _S2275 = _S2272.z;
    float _S2276 = _S2273 * _S2273 + _S2274 * _S2274 + _S2275 * _S2275;
    float _S2277 = s_primal_ctx_sqrt_0(_S2276);
    float x_36 = _S2273 / _S2277;
    float3  _S2278 = make_float3 (x_36);
    float _S2279 = _S2277 * _S2277;
    float y_15 = _S2274 / _S2277;
    float z_12 = _S2275 / _S2277;
    float3  _S2280 = make_float3 (z_12);
    float _S2281 = - y_15;
    float3  _S2282 = make_float3 (_S2281);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_36 * x_36 - y_15 * y_15;
    float _S2283 = 2.0f * x_36;
    float fS1_12 = _S2283 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S2284 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_36;
    float3  _S2285 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S2286 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S2287 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S2288 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S2289 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S2289;
    float3  _S2290 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_36;
    float3  _S2291 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S2292 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S2293 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S2294 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_36 * fC1_12 - y_15 * fS1_12);
    float3  _S2295 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_36 * fS1_12 + y_15 * fC1_12);
    float3  _S2296 = make_float3 (pSH9_2);
    float3  _S2297 = make_float3 (0.0f);
    float3  _S2298 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2299;
    (&_S2299)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2281) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_36) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S2299)->differential_0 = _S2298;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2300;
    (&_S2300)->primal_0 = _S2297;
    (&_S2300)->differential_0 = _S2298;
    s_bwd_prop_max_0(&_S2299, &_S2300, v_rgb_2);
    float3  _S2301 = _S2295 * _S2299.differential_0;
    float3  _S2302 = (*sh_coeffs_12)[int(15)] * _S2299.differential_0;
    float3  _S2303 = _S2293 * _S2299.differential_0;
    float3  _S2304 = (*sh_coeffs_12)[int(14)] * _S2299.differential_0;
    float3  _S2305 = _S2291 * _S2299.differential_0;
    float3  _S2306 = (*sh_coeffs_12)[int(13)] * _S2299.differential_0;
    float3  _S2307 = _S2290 * _S2299.differential_0;
    float3  _S2308 = (*sh_coeffs_12)[int(12)] * _S2299.differential_0;
    float3  _S2309 = _S2292 * _S2299.differential_0;
    float3  _S2310 = (*sh_coeffs_12)[int(11)] * _S2299.differential_0;
    float3  _S2311 = _S2294 * _S2299.differential_0;
    float3  _S2312 = (*sh_coeffs_12)[int(10)] * _S2299.differential_0;
    float3  _S2313 = _S2296 * _S2299.differential_0;
    float3  _S2314 = (*sh_coeffs_12)[int(9)] * _S2299.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S2314.x + _S2314.y + _S2314.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S2302.x + _S2302.y + _S2302.z);
    float _S2315 = _S2312.x + _S2312.y + _S2312.z;
    float _S2316 = _S2304.x + _S2304.y + _S2304.z;
    float _S2317 = _S2310.x + _S2310.y + _S2310.z;
    float _S2318 = _S2306.x + _S2306.y + _S2306.z;
    float _S2319 = _S2308.x + _S2308.y + _S2308.z;
    float _S2320 = - s_diff_fC2_T_2;
    float3  _S2321 = _S2287 * _S2299.differential_0;
    float3  _S2322 = (*sh_coeffs_12)[int(8)] * _S2299.differential_0;
    float3  _S2323 = _S2285 * _S2299.differential_0;
    float3  _S2324 = (*sh_coeffs_12)[int(7)] * _S2299.differential_0;
    float3  _S2325 = _S2284 * _S2299.differential_0;
    float3  _S2326 = (*sh_coeffs_12)[int(6)] * _S2299.differential_0;
    float3  _S2327 = _S2286 * _S2299.differential_0;
    float3  _S2328 = (*sh_coeffs_12)[int(5)] * _S2299.differential_0;
    float3  _S2329 = _S2288 * _S2299.differential_0;
    float3  _S2330 = (*sh_coeffs_12)[int(4)] * _S2299.differential_0;
    float _S2331 = _S2328.x + _S2328.y + _S2328.z;
    float _S2332 = _S2324.x + _S2324.y + _S2324.z;
    float _S2333 = fTmp1B_12 * _S2315 + x_36 * s_diff_fS2_T_2 + y_15 * _S2320 + 0.54627424478530884f * (_S2330.x + _S2330.y + _S2330.z);
    float _S2334 = fTmp1B_12 * _S2316 + y_15 * s_diff_fS2_T_2 + x_36 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S2322.x + _S2322.y + _S2322.z);
    float _S2335 = y_15 * - _S2334;
    float _S2336 = x_36 * _S2334;
    float _S2337 = z_12 * (1.86588168144226074f * (z_12 * _S2319) + -2.28522896766662598f * (y_15 * _S2317 + x_36 * _S2318) + 0.94617468118667603f * (_S2326.x + _S2326.y + _S2326.z));
    float3  _S2338 = make_float3 (0.48860251903533936f) * _S2299.differential_0;
    float3  _S2339 = - _S2338;
    float3  _S2340 = _S2278 * _S2339;
    float3  _S2341 = (*sh_coeffs_12)[int(3)] * _S2339;
    float3  _S2342 = _S2280 * _S2338;
    float3  _S2343 = (*sh_coeffs_12)[int(2)] * _S2338;
    float3  _S2344 = _S2282 * _S2338;
    float3  _S2345 = (*sh_coeffs_12)[int(1)] * _S2338;
    float _S2346 = (_S2289 * _S2319 + 1.44530570507049561f * (fS1_12 * _S2315 + fC1_12 * _S2316) + -1.09254848957061768f * (y_15 * _S2331 + x_36 * _S2332) + _S2337 + _S2337 + _S2343.x + _S2343.y + _S2343.z) / _S2279;
    float _S2347 = _S2277 * _S2346;
    float _S2348 = (fTmp0C_12 * _S2317 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S2320 + fTmp0B_12 * _S2331 + _S2283 * _S2333 + _S2335 + _S2335 + - (_S2345.x + _S2345.y + _S2345.z)) / _S2279;
    float _S2349 = _S2277 * _S2348;
    float _S2350 = (fTmp0C_12 * _S2318 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S2332 + 2.0f * (y_15 * _S2333) + _S2336 + _S2336 + _S2341.x + _S2341.y + _S2341.z) / _S2279;
    float _S2351 = _S2277 * _S2350;
    float _S2352 = _S2275 * - _S2346 + _S2274 * - _S2348 + _S2273 * - _S2350;
    DiffPair_float_0 _S2353;
    (&_S2353)->primal_0 = _S2276;
    (&_S2353)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2353, _S2352);
    float _S2354 = _S2275 * _S2353.differential_0;
    float _S2355 = _S2274 * _S2353.differential_0;
    float _S2356 = _S2273 * _S2353.differential_0;
    float3  _S2357 = make_float3 (0.282094806432724f) * _S2299.differential_0;
    float3  _S2358 = make_float3 (_S2351 + _S2356 + _S2356, _S2349 + _S2355 + _S2355, _S2347 + _S2354 + _S2354);
    float3  _S2359 = - - _S2358;
    Matrix<float, 3, 3>  _S2360 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2361;
    (&_S2361)->primal_0 = _S2244;
    (&_S2361)->differential_0 = _S2360;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2362;
    (&_S2362)->primal_0 = t_15;
    (&_S2362)->differential_0 = _S2298;
    s_bwd_prop_mul_1(&_S2361, &_S2362, _S2359);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2363 = _S2361;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2364 = _S2362;
    float2  _S2365 = make_float2 (0.0f);
    float2  _S2366 = _S2365;
    *&((&_S2366)->y) = v_conic_2.z;
    float2  _S2367 = _S2365;
    *&((&_S2367)->y) = v_conic_2.y;
    *&((&_S2367)->x) = v_conic_2.x;
    float _S2368 = 0.5f * v_depth_2;
    DiffPair_float_0 _S2369;
    (&_S2369)->primal_0 = _S2271;
    (&_S2369)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2369, _S2368);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2370;
    (&_S2370)->primal_0 = mean_c_12;
    (&_S2370)->differential_0 = _S2298;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2371;
    (&_S2371)->primal_0 = mean_c_12;
    (&_S2371)->differential_0 = _S2298;
    s_bwd_prop_dot_0(&_S2370, &_S2371, _S2369.differential_0);
    DiffPair_float_0 _S2372;
    (&_S2372)->primal_0 = _S2270;
    (&_S2372)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2372, 0.0f);
    DiffPair_float_0 _S2373;
    (&_S2373)->primal_0 = _S2269;
    (&_S2373)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2373, 0.0f);
    DiffPair_float_0 _S2374;
    (&_S2374)->primal_0 = 3.32999992370605469f;
    (&_S2374)->differential_0 = 0.0f;
    DiffPair_float_0 _S2375;
    (&_S2375)->primal_0 = _S2268;
    (&_S2375)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2374, &_S2375, 0.0f);
    DiffPair_float_0 _S2376;
    (&_S2376)->primal_0 = _S2267;
    (&_S2376)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2376, _S2375.differential_0);
    float _S2377 = 2.0f * _S2376.differential_0;
    DiffPair_float_0 _S2378;
    (&_S2378)->primal_0 = _S2266;
    (&_S2378)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2378, _S2377);
    float _S2379 = v_opacity_2 + 254.9999847412109375f * _S2378.differential_0;
    Matrix<float, 2, 2>  _S2380 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2381 = _S2380;
    _S2381[int(1)] = _S2366;
    _S2381[int(0)] = _S2367;
    Matrix<float, 2, 2>  _S2382 = _S2381;
    FixedArray<float3 , 16>  _S2383;
    _S2383[int(0)] = _S2298;
    _S2383[int(1)] = _S2298;
    _S2383[int(2)] = _S2298;
    _S2383[int(3)] = _S2298;
    _S2383[int(4)] = _S2298;
    _S2383[int(5)] = _S2298;
    _S2383[int(6)] = _S2298;
    _S2383[int(7)] = _S2298;
    _S2383[int(8)] = _S2298;
    _S2383[int(9)] = _S2298;
    _S2383[int(10)] = _S2298;
    _S2383[int(11)] = _S2298;
    _S2383[int(12)] = _S2298;
    _S2383[int(13)] = _S2298;
    _S2383[int(14)] = _S2298;
    _S2383[int(15)] = _S2298;
    _S2383[int(7)] = _S2323;
    _S2383[int(0)] = _S2357;
    _S2383[int(1)] = _S2344;
    _S2383[int(2)] = _S2342;
    _S2383[int(3)] = _S2340;
    _S2383[int(4)] = _S2329;
    _S2383[int(5)] = _S2327;
    _S2383[int(6)] = _S2325;
    _S2383[int(15)] = _S2301;
    _S2383[int(8)] = _S2321;
    _S2383[int(9)] = _S2313;
    _S2383[int(10)] = _S2311;
    _S2383[int(11)] = _S2309;
    _S2383[int(12)] = _S2307;
    _S2383[int(13)] = _S2305;
    _S2383[int(14)] = _S2303;
    float3  _S2384 = _S2383[int(0)];
    float3  _S2385 = _S2383[int(1)];
    float3  _S2386 = _S2383[int(2)];
    float3  _S2387 = _S2383[int(3)];
    float3  _S2388 = _S2383[int(4)];
    float3  _S2389 = _S2383[int(5)];
    float3  _S2390 = _S2383[int(6)];
    float3  _S2391 = _S2383[int(7)];
    float3  _S2392 = _S2383[int(8)];
    float3  _S2393 = _S2383[int(9)];
    float3  _S2394 = _S2383[int(10)];
    float3  _S2395 = _S2383[int(11)];
    float3  _S2396 = _S2383[int(12)];
    float3  _S2397 = _S2383[int(13)];
    float3  _S2398 = _S2383[int(14)];
    float3  _S2399 = _S2383[int(15)];
    float3  _S2400 = _S2371.differential_0 + _S2370.differential_0;
    float2  _S2401 = make_float2 (0.0f, _S2372.differential_0);
    float2  _S2402 = make_float2 (_S2373.differential_0, 0.0f);
    float _S2403;
    if(antialiased_12)
    {
        float _S2404 = _S2263 * _S2379;
        _S2265 = _S2258 * _S2379;
        _S2403 = _S2404;
    }
    else
    {
        _S2265 = _S2379;
        _S2403 = 0.0f;
    }
    float _S2405 = - (_S2265 / _S2264);
    DiffPair_float_0 _S2406;
    (&_S2406)->primal_0 = _S2261;
    (&_S2406)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2406, _S2405);
    float _S2407 = - _S2406.differential_0;
    float _S2408 = invdet_9 * _S2382.rows[int(1)].y;
    float _S2409 = - (invdet_9 * _S2382.rows[int(1)].x);
    float _S2410 = - (invdet_9 * _S2382.rows[int(0)].y);
    float _S2411 = invdet_9 * _S2382.rows[int(0)].x;
    float _S2412 = - ((_S2250 * _S2382.rows[int(1)].y + _S2260 * _S2382.rows[int(1)].x + _S2259 * _S2382.rows[int(0)].y + _S2252 * _S2382.rows[int(0)].x) / _S2256);
    DiffPair_float_0 _S2413;
    (&_S2413)->primal_0 = _S2257;
    (&_S2413)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2413, _S2403);
    DiffPair_float_0 _S2414;
    (&_S2414)->primal_0 = 0.0f;
    (&_S2414)->differential_0 = 0.0f;
    DiffPair_float_0 _S2415;
    (&_S2415)->primal_0 = _S2255;
    (&_S2415)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2414, &_S2415, _S2413.differential_0);
    float _S2416 = _S2415.differential_0 / _S2256;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2416;
    float _S2417 = _S2412 + det_orig_13 * - _S2416;
    float _S2418 = - _S2417;
    float _S2419 = _S2250 * _S2417;
    float _S2420 = _S2252 * _S2417;
    Matrix<float, 2, 2>  _S2421 = _S2380;
    _S2421[int(1)] = _S2401;
    _S2421[int(0)] = _S2402;
    _S2251 = _S2421;
    *&(((&_S2251)->rows + (int(1)))->y) = 0.0f;
    float _S2422 = _S2411 + _S2419 + _S2421.rows[int(1)].y;
    *&(((&_S2251)->rows + (int(0)))->x) = 0.0f;
    float _S2423 = _S2408 + _S2420 + _S2421.rows[int(0)].x;
    float _S2424 = _S2418 + - s_diff_det_orig_T_2;
    float _S2425 = _S2409 + _S2248.rows[int(0)].y * _S2424;
    float _S2426 = _S2410 + _S2248.rows[int(1)].x * _S2424;
    float _S2427 = _S2248.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2428 = _S2422 + _S2248.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2429 = _S2365;
    *&((&_S2429)->x) = _S2425;
    *&((&_S2429)->y) = _S2428;
    float _S2430 = _S2423 + _S2427;
    float2  _S2431 = _S2365;
    *&((&_S2431)->y) = _S2426;
    *&((&_S2431)->x) = _S2430;
    float _S2432 = fy_16 * v_mean2d_2.y;
    float _S2433 = fx_16 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S2434 = _S2380;
    _S2434[int(1)] = _S2429;
    _S2434[int(0)] = _S2431;
    Matrix<float, 2, 2>  _S2435 = _S2251 + _S2434;
    Matrix<float, 2, 3>  _S2436 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2437;
    (&_S2437)->primal_0 = _S2246;
    (&_S2437)->differential_0 = _S2436;
    Matrix<float, 3, 2>  _S2438 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2439;
    (&_S2439)->primal_0 = _S2247;
    (&_S2439)->differential_0 = _S2438;
    s_bwd_prop_mul_2(&_S2437, &_S2439, _S2435);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2440;
    (&_S2440)->primal_0 = J_12;
    (&_S2440)->differential_0 = _S2436;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2441;
    (&_S2441)->primal_0 = _S2245;
    (&_S2441)->differential_0 = _S2360;
    s_bwd_prop_mul_3(&_S2440, &_S2441, _S2437.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2442;
    (&_S2442)->primal_0 = _S2243;
    (&_S2442)->differential_0 = _S2360;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2443;
    (&_S2443)->primal_0 = _S2244;
    (&_S2443)->differential_0 = _S2360;
    s_bwd_prop_mul_4(&_S2442, &_S2443, _S2441.differential_0);
    Matrix<float, 3, 3>  _S2444 = transpose_0(_S2443.differential_0 + _S2363.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2445;
    (&_S2445)->primal_0 = R_16;
    (&_S2445)->differential_0 = _S2360;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2446;
    (&_S2446)->primal_0 = _S2242;
    (&_S2446)->differential_0 = _S2360;
    s_bwd_prop_mul_4(&_S2445, &_S2446, _S2442.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2447;
    (&_S2447)->primal_0 = _S2240;
    (&_S2447)->differential_0 = _S2360;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2448;
    (&_S2448)->primal_0 = _S2241;
    (&_S2448)->differential_0 = _S2360;
    s_bwd_prop_mul_4(&_S2447, &_S2448, _S2446.differential_0);
    Matrix<float, 3, 3>  _S2449 = _S2447.differential_0 + transpose_0(_S2448.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2450;
    (&_S2450)->primal_0 = _S2239;
    (&_S2450)->differential_0 = _S2360;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2451;
    (&_S2451)->primal_0 = S_2;
    (&_S2451)->differential_0 = _S2360;
    s_bwd_prop_mul_4(&_S2450, &_S2451, _S2449);
    Matrix<float, 3, 3>  _S2452 = transpose_0(_S2450.differential_0);
    float _S2453 = 2.0f * - _S2452.rows[int(2)].z;
    float _S2454 = 2.0f * _S2452.rows[int(2)].y;
    float _S2455 = 2.0f * _S2452.rows[int(2)].x;
    float _S2456 = 2.0f * _S2452.rows[int(1)].z;
    float _S2457 = 2.0f * - _S2452.rows[int(1)].y;
    float _S2458 = 2.0f * _S2452.rows[int(1)].x;
    float _S2459 = 2.0f * _S2452.rows[int(0)].z;
    float _S2460 = 2.0f * _S2452.rows[int(0)].y;
    float _S2461 = 2.0f * - _S2452.rows[int(0)].x;
    float _S2462 = - _S2458 + _S2460;
    float _S2463 = _S2455 + - _S2459;
    float _S2464 = - _S2454 + _S2456;
    float _S2465 = _S2454 + _S2456;
    float _S2466 = _S2455 + _S2459;
    float _S2467 = _S2458 + _S2460;
    float _S2468 = quat_15.w * (_S2457 + _S2461);
    float _S2469 = quat_15.z * (_S2453 + _S2461);
    float _S2470 = quat_15.y * (_S2453 + _S2457);
    float _S2471 = quat_15.x * _S2462 + quat_15.z * _S2465 + quat_15.y * _S2466 + _S2468 + _S2468;
    float _S2472 = quat_15.x * _S2463 + quat_15.w * _S2465 + quat_15.y * _S2467 + _S2469 + _S2469;
    float _S2473 = quat_15.x * _S2464 + quat_15.w * _S2466 + quat_15.z * _S2467 + _S2470 + _S2470;
    float _S2474 = quat_15.w * _S2462 + quat_15.z * _S2463 + quat_15.y * _S2464;
    float3  _S2475 = _S2298;
    *&((&_S2475)->z) = _S2451.differential_0.rows[int(2)].z;
    *&((&_S2475)->y) = _S2451.differential_0.rows[int(1)].y;
    *&((&_S2475)->x) = _S2451.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2476;
    (&_S2476)->primal_0 = scale_14;
    (&_S2476)->differential_0 = _S2298;
    s_bwd_prop_exp_1(&_S2476, _S2475);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2477 = _S2476;
    float3  _S2478 = _S2298;
    *&((&_S2478)->y) = _S2432;
    *&((&_S2478)->x) = _S2433;
    float3  _S2479 = _S2400 + _S2478;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2480;
    (&_S2480)->primal_0 = R_16;
    (&_S2480)->differential_0 = _S2360;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2481;
    (&_S2481)->primal_0 = mean_12;
    (&_S2481)->differential_0 = _S2298;
    s_bwd_prop_mul_1(&_S2480, &_S2481, _S2479);
    float3  _S2482 = _S2479 + _S2364.differential_0;
    Matrix<float, 3, 3>  _S2483 = _S2444 + _S2445.differential_0 + _S2480.differential_0;
    float4  _S2484 = make_float4 (0.0f);
    *&((&_S2484)->w) = _S2471;
    *&((&_S2484)->z) = _S2472;
    *&((&_S2484)->y) = _S2473;
    *&((&_S2484)->x) = _S2474;
    float4  _S2485 = _S2484;
    float3  _S2486 = _S2481.differential_0 + _S2358;
    *v_mean_2 = _S2486;
    *v_quat_2 = _S2485;
    *v_scale_2 = _S2477.differential_0;
    *v_in_opacity_2 = _S2407;
    (*v_sh_coeffs_2)[int(0)] = _S2384;
    (*v_sh_coeffs_2)[int(1)] = _S2385;
    (*v_sh_coeffs_2)[int(2)] = _S2386;
    (*v_sh_coeffs_2)[int(3)] = _S2387;
    (*v_sh_coeffs_2)[int(4)] = _S2388;
    (*v_sh_coeffs_2)[int(5)] = _S2389;
    (*v_sh_coeffs_2)[int(6)] = _S2390;
    (*v_sh_coeffs_2)[int(7)] = _S2391;
    (*v_sh_coeffs_2)[int(8)] = _S2392;
    (*v_sh_coeffs_2)[int(9)] = _S2393;
    (*v_sh_coeffs_2)[int(10)] = _S2394;
    (*v_sh_coeffs_2)[int(11)] = _S2395;
    (*v_sh_coeffs_2)[int(12)] = _S2396;
    (*v_sh_coeffs_2)[int(13)] = _S2397;
    (*v_sh_coeffs_2)[int(14)] = _S2398;
    (*v_sh_coeffs_2)[int(15)] = _S2399;
    *v_R_3 = _S2483;
    *v_t_3 = _S2482;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2487;
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S2488;
};

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_17, float fy_17, float cx_17, float cy_17, FixedArray<float, 10>  * dist_coeffs_23, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 2, 2>  _S2489 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S2490 = { false, false, false };
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2491;
    (&_S2491)->_S2487 = _S2489;
    (&_S2491)->_S2488 = _S2490;
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S2492 = s_primal_ctx_exp_0(scale_15);
    float _S2493 = quat_16.y;
    float x2_16 = _S2493 * _S2493;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S2494 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S2492.x, 0.0f, 0.0f, 0.0f, _S2492.y, 0.0f, 0.0f, 0.0f, _S2492.z);
    Matrix<float, 3, 3>  _S2495 = s_primal_ctx_mul_2(_S2494, S_3);
    Matrix<float, 3, 3>  _S2496 = transpose_0(_S2495);
    Matrix<float, 3, 3>  _S2497 = s_primal_ctx_mul_2(_S2495, _S2496);
    Matrix<float, 3, 3>  _S2498 = s_primal_ctx_mul_2(R_17, _S2497);
    Matrix<float, 3, 3>  _S2499 = transpose_0(R_17);
    Matrix<float, 3, 3>  _S2500 = s_primal_ctx_mul_2(_S2498, _S2499);
    Matrix<float, 2, 2>  _S2501 = _S2489;
    float2  _S2502 = make_float2 (0.0f);
    float2  _S2503 = _S2502;
    s_primal_ctx_persp_proj_3dgs_0(mean_c_13, _S2500, fx_17, fy_17, cx_17, cy_17, dist_coeffs_23, &_S2501, &_S2503, &(&_S2491)->_S2488);
    (&_S2491)->_S2487 = _S2501;
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2504 = _S2491;
    float _S2505 = _S2491._S2487.rows[int(0)].y * _S2491._S2487.rows[int(1)].x;
    float det_orig_14 = _S2491._S2487.rows[int(0)].x * _S2491._S2487.rows[int(1)].y - _S2505;
    float _S2506 = _S2491._S2487.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2507 = _S2491._S2487;
    *&(((&_S2507)->rows + (int(0)))->x) = _S2506;
    float _S2508 = _S2491._S2487.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2507)->rows + (int(1)))->y) = _S2508;
    Matrix<float, 2, 2>  _S2509 = _S2507;
    Matrix<float, 2, 2>  _S2510 = _S2507;
    float det_blur_9 = _S2506 * _S2508 - _S2505;
    float _S2511 = det_orig_14 / det_blur_9;
    float _S2512 = det_blur_9 * det_blur_9;
    float _S2513 = s_primal_ctx_max_0(0.0f, _S2511);
    float _S2514 = s_primal_ctx_sqrt_0(_S2513);
    float _S2515 = - in_opacity_13;
    float _S2516 = 1.0f + s_primal_ctx_exp_1(_S2515);
    float _S2517 = 1.0f / _S2516;
    float _S2518 = _S2516 * _S2516;
    float _S2519;
    if(antialiased_13)
    {
        _S2519 = _S2517 * _S2514;
    }
    else
    {
        _S2519 = _S2517;
    }
    float _S2520 = _S2519 / 0.00392156885936856f;
    float _S2521 = 2.0f * s_primal_ctx_log_0(_S2520);
    float _S2522 = s_primal_ctx_sqrt_0(_S2521);
    float _S2523 = _S2509.rows[int(0)].x;
    float _S2524 = _S2510.rows[int(1)].y;
    float _S2525 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S2526 = - scale_15;
    float3  _S2527 = mean_13 - - s_primal_ctx_mul_1(_S2499, t_16);
    float _S2528 = _S2527.x;
    float _S2529 = _S2527.y;
    float _S2530 = _S2527.z;
    float _S2531 = _S2528 * _S2528 + _S2529 * _S2529 + _S2530 * _S2530;
    float _S2532 = s_primal_ctx_sqrt_0(_S2531);
    float x_37 = _S2528 / _S2532;
    float3  _S2533 = make_float3 (x_37);
    float _S2534 = _S2532 * _S2532;
    float y_16 = _S2529 / _S2532;
    float z_13 = _S2530 / _S2532;
    float3  _S2535 = make_float3 (z_13);
    float _S2536 = - y_16;
    float3  _S2537 = make_float3 (_S2536);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_37 * x_37 - y_16 * y_16;
    float _S2538 = 2.0f * x_37;
    float fS1_13 = _S2538 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S2539 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_37;
    float3  _S2540 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S2541 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2542 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2543 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S2544 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S2544;
    float3  _S2545 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_37;
    float3  _S2546 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S2547 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2548 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2549 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_37 * fC1_13 - y_16 * fS1_13);
    float3  _S2550 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_37 * fS1_13 + y_16 * fC1_13);
    float3  _S2551 = make_float3 (pSH9_3);
    float3  _S2552 = make_float3 (0.0f);
    float3  _S2553 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2554;
    (&_S2554)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2536) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_37) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S2554)->differential_0 = _S2553;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2555;
    (&_S2555)->primal_0 = _S2552;
    (&_S2555)->differential_0 = _S2553;
    s_bwd_prop_max_0(&_S2554, &_S2555, v_rgb_3);
    float3  _S2556 = _S2550 * _S2554.differential_0;
    float3  _S2557 = (*sh_coeffs_13)[int(15)] * _S2554.differential_0;
    float3  _S2558 = _S2548 * _S2554.differential_0;
    float3  _S2559 = (*sh_coeffs_13)[int(14)] * _S2554.differential_0;
    float3  _S2560 = _S2546 * _S2554.differential_0;
    float3  _S2561 = (*sh_coeffs_13)[int(13)] * _S2554.differential_0;
    float3  _S2562 = _S2545 * _S2554.differential_0;
    float3  _S2563 = (*sh_coeffs_13)[int(12)] * _S2554.differential_0;
    float3  _S2564 = _S2547 * _S2554.differential_0;
    float3  _S2565 = (*sh_coeffs_13)[int(11)] * _S2554.differential_0;
    float3  _S2566 = _S2549 * _S2554.differential_0;
    float3  _S2567 = (*sh_coeffs_13)[int(10)] * _S2554.differential_0;
    float3  _S2568 = _S2551 * _S2554.differential_0;
    float3  _S2569 = (*sh_coeffs_13)[int(9)] * _S2554.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2569.x + _S2569.y + _S2569.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2557.x + _S2557.y + _S2557.z);
    float _S2570 = _S2567.x + _S2567.y + _S2567.z;
    float _S2571 = _S2559.x + _S2559.y + _S2559.z;
    float _S2572 = _S2565.x + _S2565.y + _S2565.z;
    float _S2573 = _S2561.x + _S2561.y + _S2561.z;
    float _S2574 = _S2563.x + _S2563.y + _S2563.z;
    float _S2575 = - s_diff_fC2_T_3;
    float3  _S2576 = _S2542 * _S2554.differential_0;
    float3  _S2577 = (*sh_coeffs_13)[int(8)] * _S2554.differential_0;
    float3  _S2578 = _S2540 * _S2554.differential_0;
    float3  _S2579 = (*sh_coeffs_13)[int(7)] * _S2554.differential_0;
    float3  _S2580 = _S2539 * _S2554.differential_0;
    float3  _S2581 = (*sh_coeffs_13)[int(6)] * _S2554.differential_0;
    float3  _S2582 = _S2541 * _S2554.differential_0;
    float3  _S2583 = (*sh_coeffs_13)[int(5)] * _S2554.differential_0;
    float3  _S2584 = _S2543 * _S2554.differential_0;
    float3  _S2585 = (*sh_coeffs_13)[int(4)] * _S2554.differential_0;
    float _S2586 = _S2583.x + _S2583.y + _S2583.z;
    float _S2587 = _S2579.x + _S2579.y + _S2579.z;
    float _S2588 = fTmp1B_13 * _S2570 + x_37 * s_diff_fS2_T_3 + y_16 * _S2575 + 0.54627424478530884f * (_S2585.x + _S2585.y + _S2585.z);
    float _S2589 = fTmp1B_13 * _S2571 + y_16 * s_diff_fS2_T_3 + x_37 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2577.x + _S2577.y + _S2577.z);
    float _S2590 = y_16 * - _S2589;
    float _S2591 = x_37 * _S2589;
    float _S2592 = z_13 * (1.86588168144226074f * (z_13 * _S2574) + -2.28522896766662598f * (y_16 * _S2572 + x_37 * _S2573) + 0.94617468118667603f * (_S2581.x + _S2581.y + _S2581.z));
    float3  _S2593 = make_float3 (0.48860251903533936f) * _S2554.differential_0;
    float3  _S2594 = - _S2593;
    float3  _S2595 = _S2533 * _S2594;
    float3  _S2596 = (*sh_coeffs_13)[int(3)] * _S2594;
    float3  _S2597 = _S2535 * _S2593;
    float3  _S2598 = (*sh_coeffs_13)[int(2)] * _S2593;
    float3  _S2599 = _S2537 * _S2593;
    float3  _S2600 = (*sh_coeffs_13)[int(1)] * _S2593;
    float _S2601 = (_S2544 * _S2574 + 1.44530570507049561f * (fS1_13 * _S2570 + fC1_13 * _S2571) + -1.09254848957061768f * (y_16 * _S2586 + x_37 * _S2587) + _S2592 + _S2592 + _S2598.x + _S2598.y + _S2598.z) / _S2534;
    float _S2602 = _S2532 * _S2601;
    float _S2603 = (fTmp0C_13 * _S2572 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2575 + fTmp0B_13 * _S2586 + _S2538 * _S2588 + _S2590 + _S2590 + - (_S2600.x + _S2600.y + _S2600.z)) / _S2534;
    float _S2604 = _S2532 * _S2603;
    float _S2605 = (fTmp0C_13 * _S2573 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2587 + 2.0f * (y_16 * _S2588) + _S2591 + _S2591 + _S2596.x + _S2596.y + _S2596.z) / _S2534;
    float _S2606 = _S2532 * _S2605;
    float _S2607 = _S2530 * - _S2601 + _S2529 * - _S2603 + _S2528 * - _S2605;
    DiffPair_float_0 _S2608;
    (&_S2608)->primal_0 = _S2531;
    (&_S2608)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2608, _S2607);
    float _S2609 = _S2530 * _S2608.differential_0;
    float _S2610 = _S2529 * _S2608.differential_0;
    float _S2611 = _S2528 * _S2608.differential_0;
    float3  _S2612 = make_float3 (0.282094806432724f) * _S2554.differential_0;
    float3  _S2613 = make_float3 (_S2606 + _S2611 + _S2611, _S2604 + _S2610 + _S2610, _S2602 + _S2609 + _S2609);
    float3  _S2614 = - - _S2613;
    Matrix<float, 3, 3>  _S2615 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2616;
    (&_S2616)->primal_0 = _S2499;
    (&_S2616)->differential_0 = _S2615;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2617;
    (&_S2617)->primal_0 = t_16;
    (&_S2617)->differential_0 = _S2553;
    s_bwd_prop_mul_1(&_S2616, &_S2617, _S2614);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2618 = _S2616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2619 = _S2617;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2620;
    (&_S2620)->primal_0 = _S2526;
    (&_S2620)->differential_0 = _S2553;
    s_bwd_prop_exp_1(&_S2620, v_conic_3);
    float3  _S2621 = - _S2620.differential_0;
    float _S2622 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2623;
    (&_S2623)->primal_0 = _S2525;
    (&_S2623)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2623, _S2622);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2624;
    (&_S2624)->primal_0 = mean_c_13;
    (&_S2624)->differential_0 = _S2553;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2625;
    (&_S2625)->primal_0 = mean_c_13;
    (&_S2625)->differential_0 = _S2553;
    s_bwd_prop_dot_0(&_S2624, &_S2625, _S2623.differential_0);
    DiffPair_float_0 _S2626;
    (&_S2626)->primal_0 = _S2524;
    (&_S2626)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2626, 0.0f);
    DiffPair_float_0 _S2627;
    (&_S2627)->primal_0 = _S2523;
    (&_S2627)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2627, 0.0f);
    DiffPair_float_0 _S2628;
    (&_S2628)->primal_0 = 3.32999992370605469f;
    (&_S2628)->differential_0 = 0.0f;
    DiffPair_float_0 _S2629;
    (&_S2629)->primal_0 = _S2522;
    (&_S2629)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2628, &_S2629, 0.0f);
    DiffPair_float_0 _S2630;
    (&_S2630)->primal_0 = _S2521;
    (&_S2630)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2630, _S2629.differential_0);
    float _S2631 = 2.0f * _S2630.differential_0;
    DiffPair_float_0 _S2632;
    (&_S2632)->primal_0 = _S2520;
    (&_S2632)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2632, _S2631);
    float2  _S2633 = make_float2 (_S2627.differential_0, 0.0f);
    float _S2634 = v_opacity_3 + 254.9999847412109375f * _S2632.differential_0;
    FixedArray<float3 , 16>  _S2635;
    _S2635[int(0)] = _S2553;
    _S2635[int(1)] = _S2553;
    _S2635[int(2)] = _S2553;
    _S2635[int(3)] = _S2553;
    _S2635[int(4)] = _S2553;
    _S2635[int(5)] = _S2553;
    _S2635[int(6)] = _S2553;
    _S2635[int(7)] = _S2553;
    _S2635[int(8)] = _S2553;
    _S2635[int(9)] = _S2553;
    _S2635[int(10)] = _S2553;
    _S2635[int(11)] = _S2553;
    _S2635[int(12)] = _S2553;
    _S2635[int(13)] = _S2553;
    _S2635[int(14)] = _S2553;
    _S2635[int(15)] = _S2553;
    _S2635[int(7)] = _S2578;
    _S2635[int(0)] = _S2612;
    _S2635[int(1)] = _S2599;
    _S2635[int(2)] = _S2597;
    _S2635[int(3)] = _S2595;
    _S2635[int(4)] = _S2584;
    _S2635[int(5)] = _S2582;
    _S2635[int(6)] = _S2580;
    _S2635[int(15)] = _S2556;
    _S2635[int(8)] = _S2576;
    _S2635[int(9)] = _S2568;
    _S2635[int(10)] = _S2566;
    _S2635[int(11)] = _S2564;
    _S2635[int(12)] = _S2562;
    _S2635[int(13)] = _S2560;
    _S2635[int(14)] = _S2558;
    float3  _S2636 = _S2635[int(0)];
    float3  _S2637 = _S2635[int(1)];
    float3  _S2638 = _S2635[int(2)];
    float3  _S2639 = _S2635[int(3)];
    float3  _S2640 = _S2635[int(4)];
    float3  _S2641 = _S2635[int(5)];
    float3  _S2642 = _S2635[int(6)];
    float3  _S2643 = _S2635[int(7)];
    float3  _S2644 = _S2635[int(8)];
    float3  _S2645 = _S2635[int(9)];
    float3  _S2646 = _S2635[int(10)];
    float3  _S2647 = _S2635[int(11)];
    float3  _S2648 = _S2635[int(12)];
    float3  _S2649 = _S2635[int(13)];
    float3  _S2650 = _S2635[int(14)];
    float3  _S2651 = _S2635[int(15)];
    float3  _S2652 = _S2625.differential_0 + _S2624.differential_0;
    float2  _S2653 = make_float2 (0.0f, _S2626.differential_0);
    float _S2654;
    if(antialiased_13)
    {
        float _S2655 = _S2517 * _S2634;
        _S2519 = _S2514 * _S2634;
        _S2654 = _S2655;
    }
    else
    {
        _S2519 = _S2634;
        _S2654 = 0.0f;
    }
    float _S2656 = - (_S2519 / _S2518);
    DiffPair_float_0 _S2657;
    (&_S2657)->primal_0 = _S2515;
    (&_S2657)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2657, _S2656);
    float _S2658 = - _S2657.differential_0;
    DiffPair_float_0 _S2659;
    (&_S2659)->primal_0 = _S2513;
    (&_S2659)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2659, _S2654);
    DiffPair_float_0 _S2660;
    (&_S2660)->primal_0 = 0.0f;
    (&_S2660)->differential_0 = 0.0f;
    DiffPair_float_0 _S2661;
    (&_S2661)->primal_0 = _S2511;
    (&_S2661)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2660, &_S2661, _S2659.differential_0);
    float _S2662 = _S2661.differential_0 / _S2512;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2662;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2662;
    float _S2663 = - s_diff_det_blur_T_0;
    float _S2664 = _S2506 * s_diff_det_blur_T_0;
    float _S2665 = _S2508 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2666 = _S2489;
    _S2666[int(1)] = _S2653;
    _S2666[int(0)] = _S2633;
    _S2507 = _S2666;
    *&(((&_S2507)->rows + (int(1)))->y) = 0.0f;
    float _S2667 = _S2664 + _S2666.rows[int(1)].y;
    *&(((&_S2507)->rows + (int(0)))->x) = 0.0f;
    float _S2668 = _S2665 + _S2666.rows[int(0)].x;
    float _S2669 = _S2663 + - s_diff_det_orig_T_3;
    float _S2670 = _S2504._S2487.rows[int(0)].y * _S2669;
    float _S2671 = _S2504._S2487.rows[int(1)].x * _S2669;
    float _S2672 = _S2504._S2487.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2673 = _S2667 + _S2504._S2487.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2674 = _S2502;
    *&((&_S2674)->x) = _S2670;
    *&((&_S2674)->y) = _S2673;
    float _S2675 = _S2668 + _S2672;
    float2  _S2676 = _S2502;
    *&((&_S2676)->y) = _S2671;
    *&((&_S2676)->x) = _S2675;
    Matrix<float, 2, 2>  _S2677 = _S2489;
    _S2677[int(1)] = _S2674;
    _S2677[int(0)] = _S2676;
    Matrix<float, 2, 2>  _S2678 = _S2507 + _S2677;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2679;
    (&_S2679)->primal_0 = mean_c_13;
    (&_S2679)->differential_0 = _S2553;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2680;
    (&_S2680)->primal_0 = _S2500;
    (&_S2680)->differential_0 = _S2615;
    DiffPair_float_0 _S2681;
    (&_S2681)->primal_0 = fx_17;
    (&_S2681)->differential_0 = 0.0f;
    DiffPair_float_0 _S2682;
    (&_S2682)->primal_0 = fy_17;
    (&_S2682)->differential_0 = 0.0f;
    DiffPair_float_0 _S2683;
    (&_S2683)->primal_0 = cx_17;
    (&_S2683)->differential_0 = 0.0f;
    DiffPair_float_0 _S2684;
    (&_S2684)->primal_0 = cy_17;
    (&_S2684)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S2685 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2686;
    (&_S2686)->primal_0 = *dist_coeffs_23;
    (&_S2686)->differential_0 = _S2685;
    s_bwd_prop_persp_proj_3dgs_Intermediates_0 _S2687 = _S2504._S2488;
    s_bwd_prop_persp_proj_3dgs_0(&_S2679, &_S2680, &_S2681, &_S2682, &_S2683, &_S2684, &_S2686, _S2678, v_mean2d_3, &_S2687);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2688;
    (&_S2688)->primal_0 = _S2498;
    (&_S2688)->differential_0 = _S2615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2689;
    (&_S2689)->primal_0 = _S2499;
    (&_S2689)->differential_0 = _S2615;
    s_bwd_prop_mul_4(&_S2688, &_S2689, _S2680.differential_0);
    Matrix<float, 3, 3>  _S2690 = transpose_0(_S2689.differential_0 + _S2618.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2691;
    (&_S2691)->primal_0 = R_17;
    (&_S2691)->differential_0 = _S2615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2692;
    (&_S2692)->primal_0 = _S2497;
    (&_S2692)->differential_0 = _S2615;
    s_bwd_prop_mul_4(&_S2691, &_S2692, _S2688.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2693;
    (&_S2693)->primal_0 = _S2495;
    (&_S2693)->differential_0 = _S2615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2694;
    (&_S2694)->primal_0 = _S2496;
    (&_S2694)->differential_0 = _S2615;
    s_bwd_prop_mul_4(&_S2693, &_S2694, _S2692.differential_0);
    Matrix<float, 3, 3>  _S2695 = _S2693.differential_0 + transpose_0(_S2694.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2696;
    (&_S2696)->primal_0 = _S2494;
    (&_S2696)->differential_0 = _S2615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2697;
    (&_S2697)->primal_0 = S_3;
    (&_S2697)->differential_0 = _S2615;
    s_bwd_prop_mul_4(&_S2696, &_S2697, _S2695);
    Matrix<float, 3, 3>  _S2698 = transpose_0(_S2696.differential_0);
    float _S2699 = 2.0f * - _S2698.rows[int(2)].z;
    float _S2700 = 2.0f * _S2698.rows[int(2)].y;
    float _S2701 = 2.0f * _S2698.rows[int(2)].x;
    float _S2702 = 2.0f * _S2698.rows[int(1)].z;
    float _S2703 = 2.0f * - _S2698.rows[int(1)].y;
    float _S2704 = 2.0f * _S2698.rows[int(1)].x;
    float _S2705 = 2.0f * _S2698.rows[int(0)].z;
    float _S2706 = 2.0f * _S2698.rows[int(0)].y;
    float _S2707 = 2.0f * - _S2698.rows[int(0)].x;
    float _S2708 = - _S2704 + _S2706;
    float _S2709 = _S2701 + - _S2705;
    float _S2710 = - _S2700 + _S2702;
    float _S2711 = _S2700 + _S2702;
    float _S2712 = _S2701 + _S2705;
    float _S2713 = _S2704 + _S2706;
    float _S2714 = quat_16.w * (_S2703 + _S2707);
    float _S2715 = quat_16.z * (_S2699 + _S2707);
    float _S2716 = quat_16.y * (_S2699 + _S2703);
    float _S2717 = quat_16.x * _S2708 + quat_16.z * _S2711 + quat_16.y * _S2712 + _S2714 + _S2714;
    float _S2718 = quat_16.x * _S2709 + quat_16.w * _S2711 + quat_16.y * _S2713 + _S2715 + _S2715;
    float _S2719 = quat_16.x * _S2710 + quat_16.w * _S2712 + quat_16.z * _S2713 + _S2716 + _S2716;
    float _S2720 = quat_16.w * _S2708 + quat_16.z * _S2709 + quat_16.y * _S2710;
    float3  _S2721 = _S2553;
    *&((&_S2721)->z) = _S2697.differential_0.rows[int(2)].z;
    *&((&_S2721)->y) = _S2697.differential_0.rows[int(1)].y;
    *&((&_S2721)->x) = _S2697.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2722;
    (&_S2722)->primal_0 = scale_15;
    (&_S2722)->differential_0 = _S2553;
    s_bwd_prop_exp_1(&_S2722, _S2721);
    float3  _S2723 = _S2679.differential_0 + _S2652;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2724;
    (&_S2724)->primal_0 = R_17;
    (&_S2724)->differential_0 = _S2615;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2725;
    (&_S2725)->primal_0 = mean_13;
    (&_S2725)->differential_0 = _S2553;
    s_bwd_prop_mul_1(&_S2724, &_S2725, _S2723);
    float3  _S2726 = _S2723 + _S2619.differential_0;
    Matrix<float, 3, 3>  _S2727 = _S2690 + _S2691.differential_0 + _S2724.differential_0;
    float3  _S2728 = _S2722.differential_0 + _S2621;
    float4  _S2729 = make_float4 (0.0f);
    *&((&_S2729)->w) = _S2717;
    *&((&_S2729)->z) = _S2718;
    *&((&_S2729)->y) = _S2719;
    *&((&_S2729)->x) = _S2720;
    float4  _S2730 = _S2729;
    float3  _S2731 = _S2725.differential_0 + _S2613;
    *v_mean_3 = _S2731;
    *v_quat_3 = _S2730;
    *v_scale_3 = _S2728;
    *v_in_opacity_3 = _S2658;
    (*v_sh_coeffs_3)[int(0)] = _S2636;
    (*v_sh_coeffs_3)[int(1)] = _S2637;
    (*v_sh_coeffs_3)[int(2)] = _S2638;
    (*v_sh_coeffs_3)[int(3)] = _S2639;
    (*v_sh_coeffs_3)[int(4)] = _S2640;
    (*v_sh_coeffs_3)[int(5)] = _S2641;
    (*v_sh_coeffs_3)[int(6)] = _S2642;
    (*v_sh_coeffs_3)[int(7)] = _S2643;
    (*v_sh_coeffs_3)[int(8)] = _S2644;
    (*v_sh_coeffs_3)[int(9)] = _S2645;
    (*v_sh_coeffs_3)[int(10)] = _S2646;
    (*v_sh_coeffs_3)[int(11)] = _S2647;
    (*v_sh_coeffs_3)[int(12)] = _S2648;
    (*v_sh_coeffs_3)[int(13)] = _S2649;
    (*v_sh_coeffs_3)[int(14)] = _S2650;
    (*v_sh_coeffs_3)[int(15)] = _S2651;
    *v_R_4 = _S2727;
    *v_t_4 = _S2726;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2732;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2733;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_18, float fy_18, float cx_18, float cy_18, FixedArray<float, 10>  * dist_coeffs_24, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2734 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S2735 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S2736 = { _S2735, _S2735 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2737 = { _S2735, _S2735, _S2736, _S2735, _S2735, _S2736 };
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2738;
    (&_S2738)->_S2732 = _S2734;
    (&_S2738)->_S2733 = _S2737;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2739 = s_primal_ctx_exp_0(scale_16);
    float _S2740 = quat_17.y;
    float x2_17 = _S2740 * _S2740;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2741 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    Matrix<float, 3, 3>  S_4 = makeMatrix<float, 3, 3> (_S2739.x, 0.0f, 0.0f, 0.0f, _S2739.y, 0.0f, 0.0f, 0.0f, _S2739.z);
    Matrix<float, 3, 3>  _S2742 = s_primal_ctx_mul_2(_S2741, S_4);
    Matrix<float, 3, 3>  _S2743 = transpose_0(_S2742);
    Matrix<float, 3, 3>  _S2744 = s_primal_ctx_mul_2(_S2742, _S2743);
    Matrix<float, 3, 3>  _S2745 = s_primal_ctx_mul_2(R_18, _S2744);
    Matrix<float, 3, 3>  _S2746 = transpose_0(R_18);
    Matrix<float, 3, 3>  _S2747 = s_primal_ctx_mul_2(_S2745, _S2746);
    Matrix<float, 2, 2>  _S2748 = _S2734;
    float2  _S2749 = make_float2 (0.0f);
    float2  _S2750 = _S2749;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_14, _S2747, fx_18, fy_18, cx_18, cy_18, dist_coeffs_24, &_S2748, &_S2750, &(&_S2738)->_S2733);
    (&_S2738)->_S2732 = _S2748;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2751 = _S2738;
    float _S2752 = _S2738._S2732.rows[int(0)].y * _S2738._S2732.rows[int(1)].x;
    float det_orig_15 = _S2738._S2732.rows[int(0)].x * _S2738._S2732.rows[int(1)].y - _S2752;
    float _S2753 = _S2738._S2732.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2754 = _S2738._S2732;
    *&(((&_S2754)->rows + (int(0)))->x) = _S2753;
    float _S2755 = _S2738._S2732.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2754)->rows + (int(1)))->y) = _S2755;
    Matrix<float, 2, 2>  _S2756 = _S2754;
    Matrix<float, 2, 2>  _S2757 = _S2754;
    float det_blur_10 = _S2753 * _S2755 - _S2752;
    float _S2758 = det_orig_15 / det_blur_10;
    float _S2759 = det_blur_10 * det_blur_10;
    float _S2760 = s_primal_ctx_max_0(0.0f, _S2758);
    float _S2761 = s_primal_ctx_sqrt_0(_S2760);
    float _S2762 = - in_opacity_14;
    float _S2763 = 1.0f + s_primal_ctx_exp_1(_S2762);
    float _S2764 = 1.0f / _S2763;
    float _S2765 = _S2763 * _S2763;
    float _S2766;
    if(antialiased_14)
    {
        _S2766 = _S2764 * _S2761;
    }
    else
    {
        _S2766 = _S2764;
    }
    float _S2767 = _S2766 / 0.00392156885936856f;
    float _S2768 = 2.0f * s_primal_ctx_log_0(_S2767);
    float _S2769 = s_primal_ctx_sqrt_0(_S2768);
    float _S2770 = _S2756.rows[int(0)].x;
    float _S2771 = _S2757.rows[int(1)].y;
    float _S2772 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2773 = - scale_16;
    float3  _S2774 = mean_14 - - s_primal_ctx_mul_1(_S2746, t_17);
    float _S2775 = _S2774.x;
    float _S2776 = _S2774.y;
    float _S2777 = _S2774.z;
    float _S2778 = _S2775 * _S2775 + _S2776 * _S2776 + _S2777 * _S2777;
    float _S2779 = s_primal_ctx_sqrt_0(_S2778);
    float x_38 = _S2775 / _S2779;
    float3  _S2780 = make_float3 (x_38);
    float _S2781 = _S2779 * _S2779;
    float y_17 = _S2776 / _S2779;
    float z_14 = _S2777 / _S2779;
    float3  _S2782 = make_float3 (z_14);
    float _S2783 = - y_17;
    float3  _S2784 = make_float3 (_S2783);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_38 * x_38 - y_17 * y_17;
    float _S2785 = 2.0f * x_38;
    float fS1_14 = _S2785 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2786 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_38;
    float3  _S2787 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2788 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2789 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2790 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2791 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2791;
    float3  _S2792 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_38;
    float3  _S2793 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2794 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2795 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2796 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_38 * fC1_14 - y_17 * fS1_14);
    float3  _S2797 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_38 * fS1_14 + y_17 * fC1_14);
    float3  _S2798 = make_float3 (pSH9_4);
    float3  _S2799 = make_float3 (0.0f);
    float3  _S2800 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2801;
    (&_S2801)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2783) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_38) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2801)->differential_0 = _S2800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2802;
    (&_S2802)->primal_0 = _S2799;
    (&_S2802)->differential_0 = _S2800;
    s_bwd_prop_max_0(&_S2801, &_S2802, v_rgb_4);
    float3  _S2803 = _S2797 * _S2801.differential_0;
    float3  _S2804 = (*sh_coeffs_14)[int(15)] * _S2801.differential_0;
    float3  _S2805 = _S2795 * _S2801.differential_0;
    float3  _S2806 = (*sh_coeffs_14)[int(14)] * _S2801.differential_0;
    float3  _S2807 = _S2793 * _S2801.differential_0;
    float3  _S2808 = (*sh_coeffs_14)[int(13)] * _S2801.differential_0;
    float3  _S2809 = _S2792 * _S2801.differential_0;
    float3  _S2810 = (*sh_coeffs_14)[int(12)] * _S2801.differential_0;
    float3  _S2811 = _S2794 * _S2801.differential_0;
    float3  _S2812 = (*sh_coeffs_14)[int(11)] * _S2801.differential_0;
    float3  _S2813 = _S2796 * _S2801.differential_0;
    float3  _S2814 = (*sh_coeffs_14)[int(10)] * _S2801.differential_0;
    float3  _S2815 = _S2798 * _S2801.differential_0;
    float3  _S2816 = (*sh_coeffs_14)[int(9)] * _S2801.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2816.x + _S2816.y + _S2816.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2804.x + _S2804.y + _S2804.z);
    float _S2817 = _S2814.x + _S2814.y + _S2814.z;
    float _S2818 = _S2806.x + _S2806.y + _S2806.z;
    float _S2819 = _S2812.x + _S2812.y + _S2812.z;
    float _S2820 = _S2808.x + _S2808.y + _S2808.z;
    float _S2821 = _S2810.x + _S2810.y + _S2810.z;
    float _S2822 = - s_diff_fC2_T_4;
    float3  _S2823 = _S2789 * _S2801.differential_0;
    float3  _S2824 = (*sh_coeffs_14)[int(8)] * _S2801.differential_0;
    float3  _S2825 = _S2787 * _S2801.differential_0;
    float3  _S2826 = (*sh_coeffs_14)[int(7)] * _S2801.differential_0;
    float3  _S2827 = _S2786 * _S2801.differential_0;
    float3  _S2828 = (*sh_coeffs_14)[int(6)] * _S2801.differential_0;
    float3  _S2829 = _S2788 * _S2801.differential_0;
    float3  _S2830 = (*sh_coeffs_14)[int(5)] * _S2801.differential_0;
    float3  _S2831 = _S2790 * _S2801.differential_0;
    float3  _S2832 = (*sh_coeffs_14)[int(4)] * _S2801.differential_0;
    float _S2833 = _S2830.x + _S2830.y + _S2830.z;
    float _S2834 = _S2826.x + _S2826.y + _S2826.z;
    float _S2835 = fTmp1B_14 * _S2817 + x_38 * s_diff_fS2_T_4 + y_17 * _S2822 + 0.54627424478530884f * (_S2832.x + _S2832.y + _S2832.z);
    float _S2836 = fTmp1B_14 * _S2818 + y_17 * s_diff_fS2_T_4 + x_38 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2824.x + _S2824.y + _S2824.z);
    float _S2837 = y_17 * - _S2836;
    float _S2838 = x_38 * _S2836;
    float _S2839 = z_14 * (1.86588168144226074f * (z_14 * _S2821) + -2.28522896766662598f * (y_17 * _S2819 + x_38 * _S2820) + 0.94617468118667603f * (_S2828.x + _S2828.y + _S2828.z));
    float3  _S2840 = make_float3 (0.48860251903533936f) * _S2801.differential_0;
    float3  _S2841 = - _S2840;
    float3  _S2842 = _S2780 * _S2841;
    float3  _S2843 = (*sh_coeffs_14)[int(3)] * _S2841;
    float3  _S2844 = _S2782 * _S2840;
    float3  _S2845 = (*sh_coeffs_14)[int(2)] * _S2840;
    float3  _S2846 = _S2784 * _S2840;
    float3  _S2847 = (*sh_coeffs_14)[int(1)] * _S2840;
    float _S2848 = (_S2791 * _S2821 + 1.44530570507049561f * (fS1_14 * _S2817 + fC1_14 * _S2818) + -1.09254848957061768f * (y_17 * _S2833 + x_38 * _S2834) + _S2839 + _S2839 + _S2845.x + _S2845.y + _S2845.z) / _S2781;
    float _S2849 = _S2779 * _S2848;
    float _S2850 = (fTmp0C_14 * _S2819 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2822 + fTmp0B_14 * _S2833 + _S2785 * _S2835 + _S2837 + _S2837 + - (_S2847.x + _S2847.y + _S2847.z)) / _S2781;
    float _S2851 = _S2779 * _S2850;
    float _S2852 = (fTmp0C_14 * _S2820 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2834 + 2.0f * (y_17 * _S2835) + _S2838 + _S2838 + _S2843.x + _S2843.y + _S2843.z) / _S2781;
    float _S2853 = _S2779 * _S2852;
    float _S2854 = _S2777 * - _S2848 + _S2776 * - _S2850 + _S2775 * - _S2852;
    DiffPair_float_0 _S2855;
    (&_S2855)->primal_0 = _S2778;
    (&_S2855)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2855, _S2854);
    float _S2856 = _S2777 * _S2855.differential_0;
    float _S2857 = _S2776 * _S2855.differential_0;
    float _S2858 = _S2775 * _S2855.differential_0;
    float3  _S2859 = make_float3 (0.282094806432724f) * _S2801.differential_0;
    float3  _S2860 = make_float3 (_S2853 + _S2858 + _S2858, _S2851 + _S2857 + _S2857, _S2849 + _S2856 + _S2856);
    float3  _S2861 = - - _S2860;
    Matrix<float, 3, 3>  _S2862 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2863;
    (&_S2863)->primal_0 = _S2746;
    (&_S2863)->differential_0 = _S2862;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2864;
    (&_S2864)->primal_0 = t_17;
    (&_S2864)->differential_0 = _S2800;
    s_bwd_prop_mul_1(&_S2863, &_S2864, _S2861);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2865 = _S2863;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2866 = _S2864;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2867;
    (&_S2867)->primal_0 = _S2773;
    (&_S2867)->differential_0 = _S2800;
    s_bwd_prop_exp_1(&_S2867, v_conic_4);
    float3  _S2868 = - _S2867.differential_0;
    float _S2869 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2870;
    (&_S2870)->primal_0 = _S2772;
    (&_S2870)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2870, _S2869);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2871;
    (&_S2871)->primal_0 = mean_c_14;
    (&_S2871)->differential_0 = _S2800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2872;
    (&_S2872)->primal_0 = mean_c_14;
    (&_S2872)->differential_0 = _S2800;
    s_bwd_prop_dot_0(&_S2871, &_S2872, _S2870.differential_0);
    DiffPair_float_0 _S2873;
    (&_S2873)->primal_0 = _S2771;
    (&_S2873)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2873, 0.0f);
    DiffPair_float_0 _S2874;
    (&_S2874)->primal_0 = _S2770;
    (&_S2874)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2874, 0.0f);
    DiffPair_float_0 _S2875;
    (&_S2875)->primal_0 = 3.32999992370605469f;
    (&_S2875)->differential_0 = 0.0f;
    DiffPair_float_0 _S2876;
    (&_S2876)->primal_0 = _S2769;
    (&_S2876)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2875, &_S2876, 0.0f);
    DiffPair_float_0 _S2877;
    (&_S2877)->primal_0 = _S2768;
    (&_S2877)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2877, _S2876.differential_0);
    float _S2878 = 2.0f * _S2877.differential_0;
    DiffPair_float_0 _S2879;
    (&_S2879)->primal_0 = _S2767;
    (&_S2879)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2879, _S2878);
    float2  _S2880 = make_float2 (_S2874.differential_0, 0.0f);
    float _S2881 = v_opacity_4 + 254.9999847412109375f * _S2879.differential_0;
    FixedArray<float3 , 16>  _S2882;
    _S2882[int(0)] = _S2800;
    _S2882[int(1)] = _S2800;
    _S2882[int(2)] = _S2800;
    _S2882[int(3)] = _S2800;
    _S2882[int(4)] = _S2800;
    _S2882[int(5)] = _S2800;
    _S2882[int(6)] = _S2800;
    _S2882[int(7)] = _S2800;
    _S2882[int(8)] = _S2800;
    _S2882[int(9)] = _S2800;
    _S2882[int(10)] = _S2800;
    _S2882[int(11)] = _S2800;
    _S2882[int(12)] = _S2800;
    _S2882[int(13)] = _S2800;
    _S2882[int(14)] = _S2800;
    _S2882[int(15)] = _S2800;
    _S2882[int(7)] = _S2825;
    _S2882[int(0)] = _S2859;
    _S2882[int(1)] = _S2846;
    _S2882[int(2)] = _S2844;
    _S2882[int(3)] = _S2842;
    _S2882[int(4)] = _S2831;
    _S2882[int(5)] = _S2829;
    _S2882[int(6)] = _S2827;
    _S2882[int(15)] = _S2803;
    _S2882[int(8)] = _S2823;
    _S2882[int(9)] = _S2815;
    _S2882[int(10)] = _S2813;
    _S2882[int(11)] = _S2811;
    _S2882[int(12)] = _S2809;
    _S2882[int(13)] = _S2807;
    _S2882[int(14)] = _S2805;
    float3  _S2883 = _S2882[int(0)];
    float3  _S2884 = _S2882[int(1)];
    float3  _S2885 = _S2882[int(2)];
    float3  _S2886 = _S2882[int(3)];
    float3  _S2887 = _S2882[int(4)];
    float3  _S2888 = _S2882[int(5)];
    float3  _S2889 = _S2882[int(6)];
    float3  _S2890 = _S2882[int(7)];
    float3  _S2891 = _S2882[int(8)];
    float3  _S2892 = _S2882[int(9)];
    float3  _S2893 = _S2882[int(10)];
    float3  _S2894 = _S2882[int(11)];
    float3  _S2895 = _S2882[int(12)];
    float3  _S2896 = _S2882[int(13)];
    float3  _S2897 = _S2882[int(14)];
    float3  _S2898 = _S2882[int(15)];
    float3  _S2899 = _S2872.differential_0 + _S2871.differential_0;
    float2  _S2900 = make_float2 (0.0f, _S2873.differential_0);
    float _S2901;
    if(antialiased_14)
    {
        float _S2902 = _S2764 * _S2881;
        _S2766 = _S2761 * _S2881;
        _S2901 = _S2902;
    }
    else
    {
        _S2766 = _S2881;
        _S2901 = 0.0f;
    }
    float _S2903 = - (_S2766 / _S2765);
    DiffPair_float_0 _S2904;
    (&_S2904)->primal_0 = _S2762;
    (&_S2904)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2904, _S2903);
    float _S2905 = - _S2904.differential_0;
    DiffPair_float_0 _S2906;
    (&_S2906)->primal_0 = _S2760;
    (&_S2906)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2906, _S2901);
    DiffPair_float_0 _S2907;
    (&_S2907)->primal_0 = 0.0f;
    (&_S2907)->differential_0 = 0.0f;
    DiffPair_float_0 _S2908;
    (&_S2908)->primal_0 = _S2758;
    (&_S2908)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2907, &_S2908, _S2906.differential_0);
    float _S2909 = _S2908.differential_0 / _S2759;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2909;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2909;
    float _S2910 = - s_diff_det_blur_T_1;
    float _S2911 = _S2753 * s_diff_det_blur_T_1;
    float _S2912 = _S2755 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2913 = _S2734;
    _S2913[int(1)] = _S2900;
    _S2913[int(0)] = _S2880;
    _S2754 = _S2913;
    *&(((&_S2754)->rows + (int(1)))->y) = 0.0f;
    float _S2914 = _S2911 + _S2913.rows[int(1)].y;
    *&(((&_S2754)->rows + (int(0)))->x) = 0.0f;
    float _S2915 = _S2912 + _S2913.rows[int(0)].x;
    float _S2916 = _S2910 + - s_diff_det_orig_T_4;
    float _S2917 = _S2751._S2732.rows[int(0)].y * _S2916;
    float _S2918 = _S2751._S2732.rows[int(1)].x * _S2916;
    float _S2919 = _S2751._S2732.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2920 = _S2914 + _S2751._S2732.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2921 = _S2749;
    *&((&_S2921)->x) = _S2917;
    *&((&_S2921)->y) = _S2920;
    float _S2922 = _S2915 + _S2919;
    float2  _S2923 = _S2749;
    *&((&_S2923)->y) = _S2918;
    *&((&_S2923)->x) = _S2922;
    Matrix<float, 2, 2>  _S2924 = _S2734;
    _S2924[int(1)] = _S2921;
    _S2924[int(0)] = _S2923;
    Matrix<float, 2, 2>  _S2925 = _S2754 + _S2924;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2926;
    (&_S2926)->primal_0 = mean_c_14;
    (&_S2926)->differential_0 = _S2800;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2927;
    (&_S2927)->primal_0 = _S2747;
    (&_S2927)->differential_0 = _S2862;
    DiffPair_float_0 _S2928;
    (&_S2928)->primal_0 = fx_18;
    (&_S2928)->differential_0 = 0.0f;
    DiffPair_float_0 _S2929;
    (&_S2929)->primal_0 = fy_18;
    (&_S2929)->differential_0 = 0.0f;
    DiffPair_float_0 _S2930;
    (&_S2930)->primal_0 = cx_18;
    (&_S2930)->differential_0 = 0.0f;
    DiffPair_float_0 _S2931;
    (&_S2931)->primal_0 = cy_18;
    (&_S2931)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S2932 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2933;
    (&_S2933)->primal_0 = *dist_coeffs_24;
    (&_S2933)->differential_0 = _S2932;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2934 = _S2751._S2733;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2926, &_S2927, &_S2928, &_S2929, &_S2930, &_S2931, &_S2933, _S2925, v_mean2d_4, &_S2934);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2935;
    (&_S2935)->primal_0 = _S2745;
    (&_S2935)->differential_0 = _S2862;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2936;
    (&_S2936)->primal_0 = _S2746;
    (&_S2936)->differential_0 = _S2862;
    s_bwd_prop_mul_4(&_S2935, &_S2936, _S2927.differential_0);
    Matrix<float, 3, 3>  _S2937 = transpose_0(_S2936.differential_0 + _S2865.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2938;
    (&_S2938)->primal_0 = R_18;
    (&_S2938)->differential_0 = _S2862;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2939;
    (&_S2939)->primal_0 = _S2744;
    (&_S2939)->differential_0 = _S2862;
    s_bwd_prop_mul_4(&_S2938, &_S2939, _S2935.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2940;
    (&_S2940)->primal_0 = _S2742;
    (&_S2940)->differential_0 = _S2862;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2941;
    (&_S2941)->primal_0 = _S2743;
    (&_S2941)->differential_0 = _S2862;
    s_bwd_prop_mul_4(&_S2940, &_S2941, _S2939.differential_0);
    Matrix<float, 3, 3>  _S2942 = _S2940.differential_0 + transpose_0(_S2941.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2943;
    (&_S2943)->primal_0 = _S2741;
    (&_S2943)->differential_0 = _S2862;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2944;
    (&_S2944)->primal_0 = S_4;
    (&_S2944)->differential_0 = _S2862;
    s_bwd_prop_mul_4(&_S2943, &_S2944, _S2942);
    Matrix<float, 3, 3>  _S2945 = transpose_0(_S2943.differential_0);
    float _S2946 = 2.0f * - _S2945.rows[int(2)].z;
    float _S2947 = 2.0f * _S2945.rows[int(2)].y;
    float _S2948 = 2.0f * _S2945.rows[int(2)].x;
    float _S2949 = 2.0f * _S2945.rows[int(1)].z;
    float _S2950 = 2.0f * - _S2945.rows[int(1)].y;
    float _S2951 = 2.0f * _S2945.rows[int(1)].x;
    float _S2952 = 2.0f * _S2945.rows[int(0)].z;
    float _S2953 = 2.0f * _S2945.rows[int(0)].y;
    float _S2954 = 2.0f * - _S2945.rows[int(0)].x;
    float _S2955 = - _S2951 + _S2953;
    float _S2956 = _S2948 + - _S2952;
    float _S2957 = - _S2947 + _S2949;
    float _S2958 = _S2947 + _S2949;
    float _S2959 = _S2948 + _S2952;
    float _S2960 = _S2951 + _S2953;
    float _S2961 = quat_17.w * (_S2950 + _S2954);
    float _S2962 = quat_17.z * (_S2946 + _S2954);
    float _S2963 = quat_17.y * (_S2946 + _S2950);
    float _S2964 = quat_17.x * _S2955 + quat_17.z * _S2958 + quat_17.y * _S2959 + _S2961 + _S2961;
    float _S2965 = quat_17.x * _S2956 + quat_17.w * _S2958 + quat_17.y * _S2960 + _S2962 + _S2962;
    float _S2966 = quat_17.x * _S2957 + quat_17.w * _S2959 + quat_17.z * _S2960 + _S2963 + _S2963;
    float _S2967 = quat_17.w * _S2955 + quat_17.z * _S2956 + quat_17.y * _S2957;
    float3  _S2968 = _S2800;
    *&((&_S2968)->z) = _S2944.differential_0.rows[int(2)].z;
    *&((&_S2968)->y) = _S2944.differential_0.rows[int(1)].y;
    *&((&_S2968)->x) = _S2944.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2969;
    (&_S2969)->primal_0 = scale_16;
    (&_S2969)->differential_0 = _S2800;
    s_bwd_prop_exp_1(&_S2969, _S2968);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2970;
    (&_S2970)->primal_0 = mean_c_14;
    (&_S2970)->differential_0 = _S2800;
    s_bwd_length_impl_0(&_S2970, 0.0f);
    float3  _S2971 = _S2926.differential_0 + _S2970.differential_0 + _S2899;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2972;
    (&_S2972)->primal_0 = R_18;
    (&_S2972)->differential_0 = _S2862;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2973;
    (&_S2973)->primal_0 = mean_14;
    (&_S2973)->differential_0 = _S2800;
    s_bwd_prop_mul_1(&_S2972, &_S2973, _S2971);
    float3  _S2974 = _S2971 + _S2866.differential_0;
    Matrix<float, 3, 3>  _S2975 = _S2937 + _S2938.differential_0 + _S2972.differential_0;
    float3  _S2976 = _S2969.differential_0 + _S2868;
    float4  _S2977 = make_float4 (0.0f);
    *&((&_S2977)->w) = _S2964;
    *&((&_S2977)->z) = _S2965;
    *&((&_S2977)->y) = _S2966;
    *&((&_S2977)->x) = _S2967;
    float4  _S2978 = _S2977;
    float3  _S2979 = _S2973.differential_0 + _S2860;
    *v_mean_4 = _S2979;
    *v_quat_4 = _S2978;
    *v_scale_4 = _S2976;
    *v_in_opacity_4 = _S2905;
    (*v_sh_coeffs_4)[int(0)] = _S2883;
    (*v_sh_coeffs_4)[int(1)] = _S2884;
    (*v_sh_coeffs_4)[int(2)] = _S2885;
    (*v_sh_coeffs_4)[int(3)] = _S2886;
    (*v_sh_coeffs_4)[int(4)] = _S2887;
    (*v_sh_coeffs_4)[int(5)] = _S2888;
    (*v_sh_coeffs_4)[int(6)] = _S2889;
    (*v_sh_coeffs_4)[int(7)] = _S2890;
    (*v_sh_coeffs_4)[int(8)] = _S2891;
    (*v_sh_coeffs_4)[int(9)] = _S2892;
    (*v_sh_coeffs_4)[int(10)] = _S2893;
    (*v_sh_coeffs_4)[int(11)] = _S2894;
    (*v_sh_coeffs_4)[int(12)] = _S2895;
    (*v_sh_coeffs_4)[int(13)] = _S2896;
    (*v_sh_coeffs_4)[int(14)] = _S2897;
    (*v_sh_coeffs_4)[int(15)] = _S2898;
    *v_R_5 = _S2975;
    *v_t_5 = _S2974;
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

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_6)
{
    float _S2980 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2980 * _S2980;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2981 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2982 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2983;
    (&_S2983)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2983)->differential_0 = _S2982;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2984;
    (&_S2984)->primal_0 = _S2981;
    (&_S2984)->differential_0 = _S2982;
    s_bwd_prop_mul_4(&_S2983, &_S2984, _s_dOut_6);
    Matrix<float, 3, 3>  _S2985 = transpose_0(transpose_0(_S2984.differential_0));
    float _S2986 = 2.0f * - _S2985.rows[int(2)].z;
    float _S2987 = 2.0f * _S2985.rows[int(2)].y;
    float _S2988 = 2.0f * _S2985.rows[int(2)].x;
    float _S2989 = 2.0f * _S2985.rows[int(1)].z;
    float _S2990 = 2.0f * - _S2985.rows[int(1)].y;
    float _S2991 = 2.0f * _S2985.rows[int(1)].x;
    float _S2992 = 2.0f * _S2985.rows[int(0)].z;
    float _S2993 = 2.0f * _S2985.rows[int(0)].y;
    float _S2994 = 2.0f * - _S2985.rows[int(0)].x;
    float _S2995 = - _S2991 + _S2993;
    float _S2996 = _S2988 + - _S2992;
    float _S2997 = - _S2987 + _S2989;
    float _S2998 = _S2987 + _S2989;
    float _S2999 = _S2988 + _S2992;
    float _S3000 = _S2991 + _S2993;
    float _S3001 = (*dpquat_0).primal_0.w * (_S2990 + _S2994);
    float _S3002 = (*dpquat_0).primal_0.z * (_S2986 + _S2994);
    float _S3003 = (*dpquat_0).primal_0.y * (_S2986 + _S2990);
    float _S3004 = (*dpquat_0).primal_0.x * _S2995 + (*dpquat_0).primal_0.z * _S2998 + (*dpquat_0).primal_0.y * _S2999 + _S3001 + _S3001;
    float _S3005 = (*dpquat_0).primal_0.x * _S2996 + (*dpquat_0).primal_0.w * _S2998 + (*dpquat_0).primal_0.y * _S3000 + _S3002 + _S3002;
    float _S3006 = (*dpquat_0).primal_0.x * _S2997 + (*dpquat_0).primal_0.w * _S2999 + (*dpquat_0).primal_0.z * _S3000 + _S3003 + _S3003;
    float _S3007 = (*dpquat_0).primal_0.w * _S2995 + (*dpquat_0).primal_0.z * _S2996 + (*dpquat_0).primal_0.y * _S2997;
    float3  _S3008 = make_float3 (_S2983.differential_0.rows[int(0)].x, _S2983.differential_0.rows[int(1)].y, _S2983.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S3008;
    float4  _S3009 = make_float4 (0.0f);
    *&((&_S3009)->w) = _S3004;
    *&((&_S3009)->z) = _S3005;
    *&((&_S3009)->y) = _S3006;
    *&((&_S3009)->x) = _S3007;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S3009;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S3010, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3011, Matrix<float, 3, 3>  _S3012)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S3010, _S3011, _S3012);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S3013 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S3013;
    float3  _S3014 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S3014;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S3015 = dOut_15.y;
    float _S3016 = dOut_15.z;
    float _S3017 = dOut_15.x;
    float _S3018 = (*a_0).primal_0.z * _S3015 + - (*a_0).primal_0.y * _S3016;
    float _S3019 = - (*a_0).primal_0.z * _S3017 + (*a_0).primal_0.x * _S3016;
    float _S3020 = (*a_0).primal_0.y * _S3017 + - (*a_0).primal_0.x * _S3015;
    float3  _S3021 = make_float3 (- (*b_0).primal_0.z * _S3015 + (*b_0).primal_0.y * _S3016, (*b_0).primal_0.z * _S3017 + - (*b_0).primal_0.x * _S3016, - (*b_0).primal_0.y * _S3017 + (*b_0).primal_0.x * _S3015);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S3021;
    float3  _S3022 = make_float3 (_S3018, _S3019, _S3020);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S3022;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S3023 = left_10.y;
    float _S3024 = right_10.z;
    float _S3025 = left_10.z;
    float _S3026 = right_10.y;
    float _S3027 = right_10.x;
    float _S3028 = left_10.x;
    return make_float3 (_S3023 * _S3024 - _S3025 * _S3026, _S3025 * _S3027 - _S3028 * _S3024, _S3028 * _S3026 - _S3023 * _S3027);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S3029, float3  _S3030)
{
    return cross_0(_S3029, _S3030);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3031, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3032, float3  _S3033)
{
    _d_cross_0(_S3031, _S3032, _S3033);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_7)
{
    float3  _S3034 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S3035 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S3034);
    float3  _S3036 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S3037 = s_primal_ctx_cross_0(_S3036, _S3035);
    float _S3038 = -0.5f * s_primal_ctx_dot_0(_S3037, _S3037);
    float _S3039 = s_primal_ctx_dot_0(_S3036, _S3036);
    float _S3040 = _S3038 / _S3039;
    float _S3041 = _S3039 * _S3039;
    float _S3042 = (*dpopacity_0).primal_0 * _s_dOut_7;
    float _S3043 = s_primal_ctx_exp_1(_S3040) * _s_dOut_7;
    DiffPair_float_0 _S3044;
    (&_S3044)->primal_0 = _S3040;
    (&_S3044)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3044, _S3042);
    float _S3045 = _S3044.differential_0 / _S3041;
    float _S3046 = _S3038 * - _S3045;
    float _S3047 = _S3039 * _S3045;
    float3  _S3048 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3049;
    (&_S3049)->primal_0 = _S3036;
    (&_S3049)->differential_0 = _S3048;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3050;
    (&_S3050)->primal_0 = _S3036;
    (&_S3050)->differential_0 = _S3048;
    s_bwd_prop_dot_0(&_S3049, &_S3050, _S3046);
    float _S3051 = -0.5f * _S3047;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3052;
    (&_S3052)->primal_0 = _S3037;
    (&_S3052)->differential_0 = _S3048;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3053;
    (&_S3053)->primal_0 = _S3037;
    (&_S3053)->differential_0 = _S3048;
    s_bwd_prop_dot_0(&_S3052, &_S3053, _S3051);
    float3  _S3054 = _S3053.differential_0 + _S3052.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3055;
    (&_S3055)->primal_0 = _S3036;
    (&_S3055)->differential_0 = _S3048;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3056;
    (&_S3056)->primal_0 = _S3035;
    (&_S3056)->differential_0 = _S3048;
    s_bwd_prop_cross_0(&_S3055, &_S3056, _S3054);
    float3  _S3057 = _S3050.differential_0 + _S3049.differential_0 + _S3055.differential_0;
    Matrix<float, 3, 3>  _S3058 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3059;
    (&_S3059)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3059)->differential_0 = _S3058;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3060;
    (&_S3060)->primal_0 = (*dpray_d_2).primal_0;
    (&_S3060)->differential_0 = _S3048;
    s_bwd_prop_mul_1(&_S3059, &_S3060, _S3057);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3061;
    (&_S3061)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3061)->differential_0 = _S3058;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3062;
    (&_S3062)->primal_0 = _S3034;
    (&_S3062)->differential_0 = _S3048;
    s_bwd_prop_mul_1(&_S3061, &_S3062, _S3056.differential_0);
    float3  _S3063 = - _S3062.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S3060.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S3062.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S3043;
    Matrix<float, 3, 3>  _S3064 = _S3059.differential_0 + _S3061.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S3064;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S3063;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3065, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3066, DiffPair_float_0 * _S3067, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3068, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3069, float _S3070)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S3065, _S3066, _S3067, _S3068, _S3069, _S3070);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S3071 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S3071;
    Matrix<float, 3, 3>  _S3072 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S3072;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S3071;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S3071;
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
    float3  _S3073 = mean_17 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S3073, _S3073) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3074 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S3075 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S3076;
    (&_S3076)->primal_0 = s_primal_ctx_dot_0(_S3074, _S3074) + 9.99999997475242708e-07f;
    (&_S3076)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3076, _S3075);
    float3  _S3077 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3078;
    (&_S3078)->primal_0 = _S3074;
    (&_S3078)->differential_0 = _S3077;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3079;
    (&_S3079)->primal_0 = _S3074;
    (&_S3079)->differential_0 = _S3077;
    s_bwd_prop_dot_0(&_S3078, &_S3079, _S3076.differential_0);
    float3  _S3080 = _S3079.differential_0 + _S3078.differential_0;
    float3  _S3081 = - _S3080;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S3077;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S3081;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3082 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3082;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3080;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3083, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3084, DiffPair_float_0 * _S3085, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3086, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3087, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3088, float3  _S3089, float _S3090)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3083, _S3084, _S3085, _S3086, _S3087, _S3088, _S3089, _S3090);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3091 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S3091;
    Matrix<float, 3, 3>  _S3092 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3092;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3091;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S3091;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S3091;
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
    float _S3093 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_16;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S3093;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S3094 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3094;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  dOut_18)
{
    float3  _S3095 = make_float3 (1.0f) / (*dpx_15).primal_0 * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3095;
    return;
}

inline __device__ float3  log_0(float3  x_40)
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
        *_slang_vector_get_element_ptr(&result_14, i_10) = (F32_log((_slang_vector_get_element(x_40, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_14;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_19, float fy_19, float cx_19, float cy_19, FixedArray<float, 10>  * dist_coeffs_25, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S3096 = mean_c_15.z;
        bool _S3097;
        if(_S3096 < near_plane_10)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3096 > far_plane_10;
        }
        if(_S3097)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3098 = scale_19.x;
        float sx_0 = (F32_exp((_S3098)));
        float _S3099 = scale_19.y;
        float sy_0 = (F32_exp((_S3099)));
        float sz_0 = scale_19.z - 0.5f * (_S3098 + _S3099);
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
        Matrix<float, 3, 3>  _S3100 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S3100, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S3100, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S3100, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S3101 = vert0_c_0.z;
        float _S3102 = vert1_c_0.z;
        float _S3103 = vert2_c_0.z;
        if(_S3101 < near_plane_10)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3101 > far_plane_10;
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3102 < near_plane_10;
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3102 > far_plane_10;
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3103 < near_plane_10;
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = _S3103 > far_plane_10;
        }
        if(_S3097)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S3104;
        for(;;)
        {
            float2  uv_21 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3101);
            if(_S3101 < 0.0f)
            {
                _S3097 = true;
            }
            else
            {
                bool _S3105 = is_valid_distortion(uv_21, dist_coeffs_25);
                _S3097 = !_S3105;
            }
            if(_S3097)
            {
                _S3104 = make_float2 (-10000.0f);
                break;
            }
            float u_22 = uv_21.x;
            float v_22 = uv_21.y;
            float r2_22 = u_22 * u_22 + v_22 * v_22;
            float2  _S3106 = uv_21 * make_float2 (1.0f + r2_22 * ((*dist_coeffs_25)[int(0)] + r2_22 * ((*dist_coeffs_25)[int(1)] + r2_22 * ((*dist_coeffs_25)[int(2)] + r2_22 * (*dist_coeffs_25)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_25)[int(4)] * u_22 * v_22 + (*dist_coeffs_25)[int(5)] * (r2_22 + 2.0f * u_22 * u_22) + (*dist_coeffs_25)[int(6)] * r2_22, 2.0f * (*dist_coeffs_25)[int(5)] * u_22 * v_22 + (*dist_coeffs_25)[int(4)] * (r2_22 + 2.0f * v_22 * v_22) + (*dist_coeffs_25)[int(7)] * r2_22);
            float2  _S3107 = _S3106 + make_float2 ((*dist_coeffs_25)[int(8)] * _S3106.x + (*dist_coeffs_25)[int(9)] * _S3106.y, 0.0f);
            _S3104 = make_float2 (fx_19 * _S3107.x + cx_19, fy_19 * _S3107.y + cy_19);
            break;
        }
        *uv0_0 = _S3104;
        for(;;)
        {
            float2  uv_22 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3102);
            if(_S3102 < 0.0f)
            {
                _S3097 = true;
            }
            else
            {
                bool _S3108 = is_valid_distortion(uv_22, dist_coeffs_25);
                _S3097 = !_S3108;
            }
            if(_S3097)
            {
                _S3104 = make_float2 (-10000.0f);
                break;
            }
            float u_23 = uv_22.x;
            float v_23 = uv_22.y;
            float r2_23 = u_23 * u_23 + v_23 * v_23;
            float2  _S3109 = uv_22 * make_float2 (1.0f + r2_23 * ((*dist_coeffs_25)[int(0)] + r2_23 * ((*dist_coeffs_25)[int(1)] + r2_23 * ((*dist_coeffs_25)[int(2)] + r2_23 * (*dist_coeffs_25)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_25)[int(4)] * u_23 * v_23 + (*dist_coeffs_25)[int(5)] * (r2_23 + 2.0f * u_23 * u_23) + (*dist_coeffs_25)[int(6)] * r2_23, 2.0f * (*dist_coeffs_25)[int(5)] * u_23 * v_23 + (*dist_coeffs_25)[int(4)] * (r2_23 + 2.0f * v_23 * v_23) + (*dist_coeffs_25)[int(7)] * r2_23);
            float2  _S3110 = _S3109 + make_float2 ((*dist_coeffs_25)[int(8)] * _S3109.x + (*dist_coeffs_25)[int(9)] * _S3109.y, 0.0f);
            _S3104 = make_float2 (fx_19 * _S3110.x + cx_19, fy_19 * _S3110.y + cy_19);
            break;
        }
        *uv1_0 = _S3104;
        for(;;)
        {
            float2  uv_23 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3103);
            if(_S3103 < 0.0f)
            {
                _S3097 = true;
            }
            else
            {
                bool _S3111 = is_valid_distortion(uv_23, dist_coeffs_25);
                _S3097 = !_S3111;
            }
            if(_S3097)
            {
                _S3104 = make_float2 (-10000.0f);
                break;
            }
            float u_24 = uv_23.x;
            float v_24 = uv_23.y;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float2  _S3112 = uv_23 * make_float2 (1.0f + r2_24 * ((*dist_coeffs_25)[int(0)] + r2_24 * ((*dist_coeffs_25)[int(1)] + r2_24 * ((*dist_coeffs_25)[int(2)] + r2_24 * (*dist_coeffs_25)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_25)[int(4)] * u_24 * v_24 + (*dist_coeffs_25)[int(5)] * (r2_24 + 2.0f * u_24 * u_24) + (*dist_coeffs_25)[int(6)] * r2_24, 2.0f * (*dist_coeffs_25)[int(5)] * u_24 * v_24 + (*dist_coeffs_25)[int(4)] * (r2_24 + 2.0f * v_24 * v_24) + (*dist_coeffs_25)[int(7)] * r2_24);
            float2  _S3113 = _S3112 + make_float2 ((*dist_coeffs_25)[int(8)] * _S3112.x + (*dist_coeffs_25)[int(9)] * _S3112.y, 0.0f);
            _S3104 = make_float2 (fx_19 * _S3113.x + cx_19, fy_19 * _S3113.y + cy_19);
            break;
        }
        *uv2_0 = _S3104;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S3104 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S3104)));
        float _S3114 = _S3104.x;
        float xmax_5 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S3114))) + offset_0;
        float xmin_5 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S3114))) - offset_0;
        float _S3115 = _S3104.y;
        float ymax_5 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S3115))) + offset_0;
        float ymin_5 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S3115))) - offset_0;
        if(xmax_5 <= 0.0f)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = xmin_5 >= float(image_width_15);
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = ymax_5 <= 0.0f;
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            _S3097 = ymin_5 >= float(image_height_15);
        }
        if(_S3097)
        {
            _S3097 = true;
        }
        else
        {
            if(_S3096 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3097 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3097 = false;
                }
                if(_S3097)
                {
                    _S3097 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3097 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3097 = false;
                    }
                }
            }
            else
            {
                _S3097 = false;
            }
        }
        if(_S3097)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3116 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S3117 = _S3116.x;
        float _S3118 = _S3116.y;
        float _S3119 = _S3116.z;
        float norm_10 = (F32_sqrt((_S3117 * _S3117 + _S3118 * _S3118 + _S3119 * _S3119)));
        float x_42 = _S3117 / norm_10;
        float y_18 = _S3118 / norm_10;
        float z_15 = _S3119 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_42 * x_42 - y_18 * y_18;
        float fS1_15 = 2.0f * x_42 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_42) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_42) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_42 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_42) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_42 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3120 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3120);
        float3  _S3121 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3122 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3121 + _S3122 + make_float3 (0.5f), _S3120);
        (*rgb_12)[int(2)] = max_0(_S3121 - _S3122 + make_float3 (0.5f), _S3120);
        float3  _S3123 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3123 * make_float3 (float(- (F32_sign((dot_0(_S3123, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_20, float fy_20, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_26, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S3124 = length_1(mean_c_16);
        bool _S3125;
        if(_S3124 < near_plane_11)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3124 > far_plane_11;
        }
        if(_S3125)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3126 = scale_20.x;
        float sx_1 = (F32_exp((_S3126)));
        float _S3127 = scale_20.y;
        float sy_1 = (F32_exp((_S3127)));
        float sz_1 = scale_20.z - 0.5f * (_S3126 + _S3127);
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
        Matrix<float, 3, 3>  _S3128 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S3128, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S3128, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S3128, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S3129 = length_1(vert0_c_1);
        float _S3130 = length_1(vert1_c_1);
        float _S3131 = length_1(vert2_c_1);
        if(_S3129 < near_plane_11)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3129 > far_plane_11;
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3130 < near_plane_11;
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3130 > far_plane_11;
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3131 < near_plane_11;
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = _S3131 > far_plane_11;
        }
        if(_S3125)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S3132 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_9 = length_0(_S3132);
        float _S3133 = vert0_c_1.z;
        float theta_4 = (F32_atan2((r_9), (_S3133)));
        float k_4;
        if(theta_4 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S3133;
        }
        else
        {
            k_4 = theta_4 / r_9;
        }
        float2  _S3134 = _S3132 * make_float2 (k_4);
        float u_25 = _S3134.x;
        float v_25 = _S3134.y;
        float r2_25 = u_25 * u_25 + v_25 * v_25;
        float _S3135 = 2.0f * (*dist_coeffs_26)[int(4)];
        float _S3136 = 2.0f * (*dist_coeffs_26)[int(5)];
        float2  _S3137 = _S3134 * make_float2 (1.0f + r2_25 * ((*dist_coeffs_26)[int(0)] + r2_25 * ((*dist_coeffs_26)[int(1)] + r2_25 * ((*dist_coeffs_26)[int(2)] + r2_25 * (*dist_coeffs_26)[int(3)])))) + make_float2 (_S3135 * u_25 * v_25 + (*dist_coeffs_26)[int(5)] * (r2_25 + 2.0f * u_25 * u_25) + (*dist_coeffs_26)[int(6)] * r2_25, _S3136 * u_25 * v_25 + (*dist_coeffs_26)[int(4)] * (r2_25 + 2.0f * v_25 * v_25) + (*dist_coeffs_26)[int(7)] * r2_25);
        float2  _S3138 = _S3137 + make_float2 ((*dist_coeffs_26)[int(8)] * _S3137.x + (*dist_coeffs_26)[int(9)] * _S3137.y, 0.0f);
        *uv0_1 = make_float2 (fx_20 * _S3138.x + cx_20, fy_20 * _S3138.y + cy_20);
        float2  _S3139 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_10 = length_0(_S3139);
        float _S3140 = vert1_c_1.z;
        float theta_5 = (F32_atan2((r_10), (_S3140)));
        if(theta_5 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_5 * theta_5 / 3.0f) / _S3140;
        }
        else
        {
            k_4 = theta_5 / r_10;
        }
        float2  _S3141 = _S3139 * make_float2 (k_4);
        float u_26 = _S3141.x;
        float v_26 = _S3141.y;
        float r2_26 = u_26 * u_26 + v_26 * v_26;
        float2  _S3142 = _S3141 * make_float2 (1.0f + r2_26 * ((*dist_coeffs_26)[int(0)] + r2_26 * ((*dist_coeffs_26)[int(1)] + r2_26 * ((*dist_coeffs_26)[int(2)] + r2_26 * (*dist_coeffs_26)[int(3)])))) + make_float2 (_S3135 * u_26 * v_26 + (*dist_coeffs_26)[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + (*dist_coeffs_26)[int(6)] * r2_26, _S3136 * u_26 * v_26 + (*dist_coeffs_26)[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + (*dist_coeffs_26)[int(7)] * r2_26);
        float2  _S3143 = _S3142 + make_float2 ((*dist_coeffs_26)[int(8)] * _S3142.x + (*dist_coeffs_26)[int(9)] * _S3142.y, 0.0f);
        *uv1_1 = make_float2 (fx_20 * _S3143.x + cx_20, fy_20 * _S3143.y + cy_20);
        float2  _S3144 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_11 = length_0(_S3144);
        float _S3145 = vert2_c_1.z;
        float theta_6 = (F32_atan2((r_11), (_S3145)));
        if(theta_6 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_6 * theta_6 / 3.0f) / _S3145;
        }
        else
        {
            k_4 = theta_6 / r_11;
        }
        float2  _S3146 = _S3144 * make_float2 (k_4);
        float u_27 = _S3146.x;
        float v_27 = _S3146.y;
        float r2_27 = u_27 * u_27 + v_27 * v_27;
        float2  _S3147 = _S3146 * make_float2 (1.0f + r2_27 * ((*dist_coeffs_26)[int(0)] + r2_27 * ((*dist_coeffs_26)[int(1)] + r2_27 * ((*dist_coeffs_26)[int(2)] + r2_27 * (*dist_coeffs_26)[int(3)])))) + make_float2 (_S3135 * u_27 * v_27 + (*dist_coeffs_26)[int(5)] * (r2_27 + 2.0f * u_27 * u_27) + (*dist_coeffs_26)[int(6)] * r2_27, _S3136 * u_27 * v_27 + (*dist_coeffs_26)[int(4)] * (r2_27 + 2.0f * v_27 * v_27) + (*dist_coeffs_26)[int(7)] * r2_27);
        float2  _S3148 = _S3147 + make_float2 ((*dist_coeffs_26)[int(8)] * _S3147.x + (*dist_coeffs_26)[int(9)] * _S3147.y, 0.0f);
        float _S3149 = fx_20 * _S3148.x + cx_20;
        float _S3150 = fy_20 * _S3148.y + cy_20;
        float2  _S3151 = make_float2 (_S3149, _S3150);
        *uv2_1 = _S3151;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S3151 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S3151)));
        float xmax_6 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S3149))) + offset_1;
        float xmin_6 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S3149))) - offset_1;
        float ymax_6 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S3150))) + offset_1;
        float ymin_6 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S3150))) - offset_1;
        if(xmax_6 <= 0.0f)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = xmin_6 >= float(image_width_16);
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = ymax_6 <= 0.0f;
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            _S3125 = ymin_6 >= float(image_height_16);
        }
        if(_S3125)
        {
            _S3125 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3125 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3125 = false;
                }
                if(_S3125)
                {
                    _S3125 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3125 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3125 = false;
                    }
                }
            }
            else
            {
                _S3125 = false;
            }
        }
        if(_S3125)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3129, _S3130, _S3131) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3152 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S3153 = _S3152.x;
        float _S3154 = _S3152.y;
        float _S3155 = _S3152.z;
        float norm_11 = (F32_sqrt((_S3153 * _S3153 + _S3154 * _S3154 + _S3155 * _S3155)));
        float x_44 = _S3153 / norm_11;
        float y_19 = _S3154 / norm_11;
        float z_16 = _S3155 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_44 * x_44 - y_19 * y_19;
        float fS1_16 = 2.0f * x_44 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_44) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_44) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_44) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3156 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3156);
        float3  _S3157 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3158 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3157 + _S3158 + make_float3 (0.5f), _S3156);
        (*rgb_13)[int(2)] = max_0(_S3157 - _S3158 + make_float3 (0.5f), _S3156);
        float3  _S3159 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3159 * make_float3 (float(- (F32_sign((dot_0(_S3159, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_21, float fy_21, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_27, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    for(;;)
    {
        float2  _S3160;
        bool _S3161;
        float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
        float _S3162 = scale_21.x;
        float sx_2 = (F32_exp((_S3162)));
        float _S3163 = scale_21.y;
        float sy_2 = (F32_exp((_S3163)));
        float sz_2 = scale_21.z - 0.5f * (_S3162 + _S3163);
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
        Matrix<float, 3, 3>  _S3164 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
        float3  vert0_c_2 = mul_0(R_21, mul_0(_S3164, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
        float3  vert1_c_2 = mul_0(R_21, mul_0(_S3164, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
        float3  vert2_c_2 = mul_0(R_21, mul_0(_S3164, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
        for(;;)
        {
            float _S3165 = vert0_c_2.z;
            float2  uv_24 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (_S3165);
            if(_S3165 < 0.0f)
            {
                _S3161 = true;
            }
            else
            {
                bool _S3166 = is_valid_distortion(uv_24, dist_coeffs_27);
                _S3161 = !_S3166;
            }
            if(_S3161)
            {
                _S3160 = make_float2 (-10000.0f);
                break;
            }
            float u_28 = uv_24.x;
            float v_28 = uv_24.y;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float2  _S3167 = uv_24 * make_float2 (1.0f + r2_28 * ((*dist_coeffs_27)[int(0)] + r2_28 * ((*dist_coeffs_27)[int(1)] + r2_28 * ((*dist_coeffs_27)[int(2)] + r2_28 * (*dist_coeffs_27)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_27)[int(4)] * u_28 * v_28 + (*dist_coeffs_27)[int(5)] * (r2_28 + 2.0f * u_28 * u_28) + (*dist_coeffs_27)[int(6)] * r2_28, 2.0f * (*dist_coeffs_27)[int(5)] * u_28 * v_28 + (*dist_coeffs_27)[int(4)] * (r2_28 + 2.0f * v_28 * v_28) + (*dist_coeffs_27)[int(7)] * r2_28);
            float2  _S3168 = _S3167 + make_float2 ((*dist_coeffs_27)[int(8)] * _S3167.x + (*dist_coeffs_27)[int(9)] * _S3167.y, 0.0f);
            _S3160 = make_float2 (fx_21 * _S3168.x + cx_21, fy_21 * _S3168.y + cy_21);
            break;
        }
        *uv0_2 = _S3160;
        for(;;)
        {
            float _S3169 = vert1_c_2.z;
            float2  uv_25 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (_S3169);
            if(_S3169 < 0.0f)
            {
                _S3161 = true;
            }
            else
            {
                bool _S3170 = is_valid_distortion(uv_25, dist_coeffs_27);
                _S3161 = !_S3170;
            }
            if(_S3161)
            {
                _S3160 = make_float2 (-10000.0f);
                break;
            }
            float u_29 = uv_25.x;
            float v_29 = uv_25.y;
            float r2_29 = u_29 * u_29 + v_29 * v_29;
            float2  _S3171 = uv_25 * make_float2 (1.0f + r2_29 * ((*dist_coeffs_27)[int(0)] + r2_29 * ((*dist_coeffs_27)[int(1)] + r2_29 * ((*dist_coeffs_27)[int(2)] + r2_29 * (*dist_coeffs_27)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_27)[int(4)] * u_29 * v_29 + (*dist_coeffs_27)[int(5)] * (r2_29 + 2.0f * u_29 * u_29) + (*dist_coeffs_27)[int(6)] * r2_29, 2.0f * (*dist_coeffs_27)[int(5)] * u_29 * v_29 + (*dist_coeffs_27)[int(4)] * (r2_29 + 2.0f * v_29 * v_29) + (*dist_coeffs_27)[int(7)] * r2_29);
            float2  _S3172 = _S3171 + make_float2 ((*dist_coeffs_27)[int(8)] * _S3171.x + (*dist_coeffs_27)[int(9)] * _S3171.y, 0.0f);
            _S3160 = make_float2 (fx_21 * _S3172.x + cx_21, fy_21 * _S3172.y + cy_21);
            break;
        }
        *uv1_2 = _S3160;
        for(;;)
        {
            float _S3173 = vert2_c_2.z;
            float2  uv_26 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (_S3173);
            if(_S3173 < 0.0f)
            {
                _S3161 = true;
            }
            else
            {
                bool _S3174 = is_valid_distortion(uv_26, dist_coeffs_27);
                _S3161 = !_S3174;
            }
            if(_S3161)
            {
                _S3160 = make_float2 (-10000.0f);
                break;
            }
            float u_30 = uv_26.x;
            float v_30 = uv_26.y;
            float r2_30 = u_30 * u_30 + v_30 * v_30;
            float2  _S3175 = uv_26 * make_float2 (1.0f + r2_30 * ((*dist_coeffs_27)[int(0)] + r2_30 * ((*dist_coeffs_27)[int(1)] + r2_30 * ((*dist_coeffs_27)[int(2)] + r2_30 * (*dist_coeffs_27)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_27)[int(4)] * u_30 * v_30 + (*dist_coeffs_27)[int(5)] * (r2_30 + 2.0f * u_30 * u_30) + (*dist_coeffs_27)[int(6)] * r2_30, 2.0f * (*dist_coeffs_27)[int(5)] * u_30 * v_30 + (*dist_coeffs_27)[int(4)] * (r2_30 + 2.0f * v_30 * v_30) + (*dist_coeffs_27)[int(7)] * r2_30);
            float2  _S3176 = _S3175 + make_float2 ((*dist_coeffs_27)[int(8)] * _S3175.x + (*dist_coeffs_27)[int(9)] * _S3175.y, 0.0f);
            _S3160 = make_float2 (fx_21 * _S3176.x + cx_21, fy_21 * _S3176.y + cy_21);
            break;
        }
        *uv2_2 = _S3160;
        float2  e0_2 = *uv1_2 - *uv0_2;
        float2  e1_2 = _S3160 - *uv1_2;
        float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3160)));
        float _S3177 = _S3160.x;
        float _S3178 = _S3160.y;
        *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3177))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3178))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3177))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3178))) + offset_2)))));
        *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_2 = hardness_2;
        float3  _S3179 = mean_21 - - mul_0(transpose_0(R_21), t_20);
        float _S3180 = _S3179.x;
        float _S3181 = _S3179.y;
        float _S3182 = _S3179.z;
        float norm_12 = (F32_sqrt((_S3180 * _S3180 + _S3181 * _S3181 + _S3182 * _S3182)));
        float x_46 = _S3180 / norm_12;
        float y_20 = _S3181 / norm_12;
        float z_17 = _S3182 / norm_12;
        float z2_40 = z_17 * z_17;
        float fTmp0B_17 = -1.09254848957061768f * z_17;
        float fC1_17 = x_46 * x_46 - y_20 * y_20;
        float fS1_17 = 2.0f * x_46 * y_20;
        float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
        float fTmp1B_17 = 1.44530570507049561f * z_17;
        float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_46) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_46) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_46) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
        float3  _S3183 = make_float3 (0.0f);
        (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3183);
        float3  _S3184 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
        float3  _S3185 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_14)[int(1)] = max_0(_S3184 + _S3185 + make_float3 (0.5f), _S3183);
        (*rgb_14)[int(2)] = max_0(_S3184 - _S3185 + make_float3 (0.5f), _S3183);
        float3  _S3186 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
        *normal_2 = _S3186 * make_float3 (float(- (F32_sign((dot_0(_S3186, mean_c_17))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_22, float fy_22, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_28, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S3187 = scale_22.x;
    float sx_3 = (F32_exp((_S3187)));
    float _S3188 = scale_22.y;
    float sy_3 = (F32_exp((_S3188)));
    float sz_3 = scale_22.z - 0.5f * (_S3187 + _S3188);
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
    Matrix<float, 3, 3>  _S3189 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S3189, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S3189, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S3189, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    float2  _S3190 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_12 = length_0(_S3190);
    float _S3191 = vert0_c_3.z;
    float theta_7 = (F32_atan2((r_12), (_S3191)));
    float k_5;
    if(theta_7 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_7 * theta_7 / 3.0f) / _S3191;
    }
    else
    {
        k_5 = theta_7 / r_12;
    }
    float2  _S3192 = _S3190 * make_float2 (k_5);
    float u_31 = _S3192.x;
    float v_31 = _S3192.y;
    float r2_31 = u_31 * u_31 + v_31 * v_31;
    float _S3193 = 2.0f * (*dist_coeffs_28)[int(4)];
    float _S3194 = 2.0f * (*dist_coeffs_28)[int(5)];
    float2  _S3195 = _S3192 * make_float2 (1.0f + r2_31 * ((*dist_coeffs_28)[int(0)] + r2_31 * ((*dist_coeffs_28)[int(1)] + r2_31 * ((*dist_coeffs_28)[int(2)] + r2_31 * (*dist_coeffs_28)[int(3)])))) + make_float2 (_S3193 * u_31 * v_31 + (*dist_coeffs_28)[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + (*dist_coeffs_28)[int(6)] * r2_31, _S3194 * u_31 * v_31 + (*dist_coeffs_28)[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + (*dist_coeffs_28)[int(7)] * r2_31);
    float2  _S3196 = _S3195 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3195.x + (*dist_coeffs_28)[int(9)] * _S3195.y, 0.0f);
    *uv0_3 = make_float2 (fx_22 * _S3196.x + cx_22, fy_22 * _S3196.y + cy_22);
    float2  _S3197 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_13 = length_0(_S3197);
    float _S3198 = vert1_c_3.z;
    float theta_8 = (F32_atan2((r_13), (_S3198)));
    if(theta_8 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_8 * theta_8 / 3.0f) / _S3198;
    }
    else
    {
        k_5 = theta_8 / r_13;
    }
    float2  _S3199 = _S3197 * make_float2 (k_5);
    float u_32 = _S3199.x;
    float v_32 = _S3199.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float2  _S3200 = _S3199 * make_float2 (1.0f + r2_32 * ((*dist_coeffs_28)[int(0)] + r2_32 * ((*dist_coeffs_28)[int(1)] + r2_32 * ((*dist_coeffs_28)[int(2)] + r2_32 * (*dist_coeffs_28)[int(3)])))) + make_float2 (_S3193 * u_32 * v_32 + (*dist_coeffs_28)[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + (*dist_coeffs_28)[int(6)] * r2_32, _S3194 * u_32 * v_32 + (*dist_coeffs_28)[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + (*dist_coeffs_28)[int(7)] * r2_32);
    float2  _S3201 = _S3200 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3200.x + (*dist_coeffs_28)[int(9)] * _S3200.y, 0.0f);
    *uv1_3 = make_float2 (fx_22 * _S3201.x + cx_22, fy_22 * _S3201.y + cy_22);
    float2  _S3202 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_14 = length_0(_S3202);
    float _S3203 = vert2_c_3.z;
    float theta_9 = (F32_atan2((r_14), (_S3203)));
    if(theta_9 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_9 * theta_9 / 3.0f) / _S3203;
    }
    else
    {
        k_5 = theta_9 / r_14;
    }
    float2  _S3204 = _S3202 * make_float2 (k_5);
    float u_33 = _S3204.x;
    float v_33 = _S3204.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S3205 = _S3204 * make_float2 (1.0f + r2_33 * ((*dist_coeffs_28)[int(0)] + r2_33 * ((*dist_coeffs_28)[int(1)] + r2_33 * ((*dist_coeffs_28)[int(2)] + r2_33 * (*dist_coeffs_28)[int(3)])))) + make_float2 (_S3193 * u_33 * v_33 + (*dist_coeffs_28)[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + (*dist_coeffs_28)[int(6)] * r2_33, _S3194 * u_33 * v_33 + (*dist_coeffs_28)[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + (*dist_coeffs_28)[int(7)] * r2_33);
    float2  _S3206 = _S3205 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3205.x + (*dist_coeffs_28)[int(9)] * _S3205.y, 0.0f);
    float _S3207 = fx_22 * _S3206.x + cx_22;
    float _S3208 = fy_22 * _S3206.y + cy_22;
    float2  _S3209 = make_float2 (_S3207, _S3208);
    *uv2_3 = _S3209;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3209 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3209)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3207))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3208))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3207))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3208))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3210 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S3211 = _S3210.x;
    float _S3212 = _S3210.y;
    float _S3213 = _S3210.z;
    float norm_13 = (F32_sqrt((_S3211 * _S3211 + _S3212 * _S3212 + _S3213 * _S3213)));
    float x_48 = _S3211 / norm_13;
    float y_21 = _S3212 / norm_13;
    float z_18 = _S3213 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_48 * x_48 - y_21 * y_21;
    float fS1_18 = 2.0f * x_48 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_48) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_48) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_48) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3214 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3214);
    float3  _S3215 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3216 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3215 + _S3216 + make_float3 (0.5f), _S3214);
    (*rgb_15)[int(2)] = max_0(_S3215 - _S3216 + make_float3 (0.5f), _S3214);
    float3  _S3217 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3217 * make_float3 (float(- (F32_sign((dot_0(_S3217, mean_c_18))))));
    return;
}

struct s_bwd_prop_projection_opaque_triangle_persp_differentiable_Intermediates_0
{
    bool _S3218;
    bool _S3219;
    bool _S3220;
};

inline __device__ float s_primal_ctx_min_0(float _S3221, float _S3222)
{
    return (F32_min((_S3221), (_S3222)));
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3223, float3  _S3224)
{
    _d_log_vector_0(_S3223, _S3224);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3225, float _S3226)
{
    _d_exp2_0(_S3225, _S3226);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3227, float _S3228)
{
    _d_abs_0(_S3227, _S3228);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_23, float fy_23, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    s_bwd_prop_projection_opaque_triangle_persp_differentiable_Intermediates_0 _S3229;
    (&_S3229)->_S3218 = false;
    (&_S3229)->_S3219 = false;
    (&_S3229)->_S3220 = false;
    (&_S3229)->_S3218 = false;
    (&_S3229)->_S3219 = false;
    (&_S3229)->_S3220 = false;
    float _S3230 = scale_23.x;
    float _S3231 = s_primal_ctx_exp_1(_S3230);
    float _S3232 = scale_23.y;
    float _S3233 = s_primal_ctx_exp_1(_S3232);
    float sz_4 = scale_23.z - 0.5f * (_S3230 + _S3232);
    float _S3234 = quat_24.y;
    float x2_24 = _S3234 * _S3234;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S3235 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S3236 = make_float3 (_S3231, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S3235, _S3236) + mean_23;
    float _S3237 = -0.5f + sz_4;
    float3  _S3238 = make_float3 (_S3231 * _S3237, _S3233, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S3235, _S3238) + mean_23;
    float _S3239 = -0.5f - sz_4;
    float3  _S3240 = make_float3 (_S3231 * _S3239, - _S3233, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S3235, _S3240) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S3241 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3242 = vert0_c_4.z;
    float2  uv_27 = _S3241 / make_float2 (_S3242);
    bool _S3243 = _S3242 < 0.0f;
    if(_S3243)
    {
    }
    else
    {
        bool _S3244 = is_valid_distortion(uv_27, dist_coeffs_29);
        (&_S3229)->_S3218 = _S3244;
    }
    float2  _S3245 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3246 = vert1_c_4.z;
    float2  uv_28 = _S3245 / make_float2 (_S3246);
    bool _S3247 = _S3246 < 0.0f;
    if(_S3247)
    {
    }
    else
    {
        bool _S3248 = is_valid_distortion(uv_28, dist_coeffs_29);
        (&_S3229)->_S3219 = _S3248;
    }
    float2  _S3249 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3250 = vert2_c_4.z;
    float2  uv_29 = _S3249 / make_float2 (_S3250);
    bool _S3251 = _S3250 < 0.0f;
    if(_S3251)
    {
    }
    else
    {
        bool _S3252 = is_valid_distortion(uv_29, dist_coeffs_29);
        (&_S3229)->_S3220 = _S3252;
    }
    s_bwd_prop_projection_opaque_triangle_persp_differentiable_Intermediates_0 _S3253 = _S3229;
    float2  _S3254 = make_float2 (0.0f);
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float2  _S3255 = make_float2 (_S3242);
    float2  uv_30 = _S3241 / make_float2 (_S3242);
    float2  _S3256 = make_float2 (_S3242 * _S3242);
    bool _S3257;
    if(_S3243)
    {
        _S3257 = true;
    }
    else
    {
        _S3257 = !_S3253._S3218;
    }
    float2  _S3258;
    if(_S3257)
    {
        _S3258 = make_float2 (-10000.0f);
    }
    bool _S3259 = !_S3257;
    float2  _S3260;
    float _S3261;
    float _S3262;
    float _S3263;
    float _S3264;
    float _S3265;
    float _S3266;
    float _S3267;
    float _S3268;
    float _S3269;
    float _S3270;
    float _S3271;
    float _S3272;
    float _S3273;
    float _S3274;
    float _S3275;
    float _S3276;
    float _S3277;
    float _S3278;
    float _S3279;
    if(_S3259)
    {
        float u_34 = uv_30.x;
        float v_34 = uv_30.y;
        float r2_34 = u_34 * u_34 + v_34 * v_34;
        float _S3280 = (*dist_coeffs_29)[int(2)] + r2_34 * (*dist_coeffs_29)[int(3)];
        float _S3281 = (*dist_coeffs_29)[int(1)] + r2_34 * _S3280;
        float _S3282 = (*dist_coeffs_29)[int(0)] + r2_34 * _S3281;
        float radial_3 = 1.0f + r2_34 * _S3282;
        float2  _S3283 = make_float2 (radial_3);
        float _S3284 = 2.0f * (*dist_coeffs_29)[int(4)];
        float _S3285 = _S3284 * u_34;
        float _S3286 = 2.0f * u_34;
        float _S3287 = 2.0f * (*dist_coeffs_29)[int(5)];
        float _S3288 = _S3287 * u_34;
        float _S3289 = 2.0f * v_34;
        float2  _S3290 = uv_30 * make_float2 (radial_3) + make_float2 (_S3285 * v_34 + (*dist_coeffs_29)[int(5)] * (r2_34 + _S3286 * u_34) + (*dist_coeffs_29)[int(6)] * r2_34, _S3288 * v_34 + (*dist_coeffs_29)[int(4)] * (r2_34 + _S3289 * v_34) + (*dist_coeffs_29)[int(7)] * r2_34);
        float2  _S3291 = _S3290 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3290.x + (*dist_coeffs_29)[int(9)] * _S3290.y, 0.0f);
        _S3258 = make_float2 (fx_23 * _S3291.x + cx_23, fy_23 * _S3291.y + cy_23);
        _S3261 = (*dist_coeffs_29)[int(9)];
        _S3262 = (*dist_coeffs_29)[int(8)];
        _S3260 = _S3283;
        _S3263 = (*dist_coeffs_29)[int(7)];
        _S3264 = (*dist_coeffs_29)[int(4)];
        _S3265 = _S3289;
        _S3266 = v_34;
        _S3267 = _S3288;
        _S3268 = _S3287;
        _S3269 = (*dist_coeffs_29)[int(6)];
        _S3270 = (*dist_coeffs_29)[int(5)];
        _S3271 = _S3286;
        _S3272 = u_34;
        _S3273 = _S3285;
        _S3274 = _S3284;
        _S3275 = r2_34;
        _S3276 = _S3282;
        _S3277 = _S3281;
        _S3278 = _S3280;
        _S3279 = (*dist_coeffs_29)[int(3)];
    }
    else
    {
        _S3261 = 0.0f;
        _S3262 = 0.0f;
        _S3260 = _S3254;
        _S3263 = 0.0f;
        _S3264 = 0.0f;
        _S3265 = 0.0f;
        _S3266 = 0.0f;
        _S3267 = 0.0f;
        _S3268 = 0.0f;
        _S3269 = 0.0f;
        _S3270 = 0.0f;
        _S3271 = 0.0f;
        _S3272 = 0.0f;
        _S3273 = 0.0f;
        _S3274 = 0.0f;
        _S3275 = 0.0f;
        _S3276 = 0.0f;
        _S3277 = 0.0f;
        _S3278 = 0.0f;
        _S3279 = 0.0f;
    }
    float2  _S3292 = make_float2 (_S3246);
    float2  uv_31 = _S3245 / make_float2 (_S3246);
    float2  _S3293 = make_float2 (_S3246 * _S3246);
    if(_S3247)
    {
        _S3257 = true;
    }
    else
    {
        _S3257 = !_S3253._S3219;
    }
    float2  _S3294;
    if(_S3257)
    {
        _S3294 = make_float2 (-10000.0f);
    }
    bool _S3295 = !_S3257;
    float2  _S3296;
    float _S3297;
    float _S3298;
    float _S3299;
    float _S3300;
    float _S3301;
    float _S3302;
    float _S3303;
    float _S3304;
    float _S3305;
    float _S3306;
    float _S3307;
    float _S3308;
    float _S3309;
    float _S3310;
    float _S3311;
    float _S3312;
    float _S3313;
    float _S3314;
    float _S3315;
    if(_S3295)
    {
        float u_35 = uv_31.x;
        float v_35 = uv_31.y;
        float r2_35 = u_35 * u_35 + v_35 * v_35;
        float _S3316 = (*dist_coeffs_29)[int(2)] + r2_35 * (*dist_coeffs_29)[int(3)];
        float _S3317 = (*dist_coeffs_29)[int(1)] + r2_35 * _S3316;
        float _S3318 = (*dist_coeffs_29)[int(0)] + r2_35 * _S3317;
        float radial_4 = 1.0f + r2_35 * _S3318;
        float2  _S3319 = make_float2 (radial_4);
        float _S3320 = 2.0f * (*dist_coeffs_29)[int(4)];
        float _S3321 = _S3320 * u_35;
        float _S3322 = 2.0f * u_35;
        float _S3323 = 2.0f * (*dist_coeffs_29)[int(5)];
        float _S3324 = _S3323 * u_35;
        float _S3325 = 2.0f * v_35;
        float2  _S3326 = uv_31 * make_float2 (radial_4) + make_float2 (_S3321 * v_35 + (*dist_coeffs_29)[int(5)] * (r2_35 + _S3322 * u_35) + (*dist_coeffs_29)[int(6)] * r2_35, _S3324 * v_35 + (*dist_coeffs_29)[int(4)] * (r2_35 + _S3325 * v_35) + (*dist_coeffs_29)[int(7)] * r2_35);
        float2  _S3327 = _S3326 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3326.x + (*dist_coeffs_29)[int(9)] * _S3326.y, 0.0f);
        _S3294 = make_float2 (fx_23 * _S3327.x + cx_23, fy_23 * _S3327.y + cy_23);
        _S3297 = (*dist_coeffs_29)[int(9)];
        _S3298 = (*dist_coeffs_29)[int(8)];
        _S3296 = _S3319;
        _S3299 = (*dist_coeffs_29)[int(7)];
        _S3300 = (*dist_coeffs_29)[int(4)];
        _S3301 = _S3325;
        _S3302 = v_35;
        _S3303 = _S3324;
        _S3304 = _S3323;
        _S3305 = (*dist_coeffs_29)[int(6)];
        _S3306 = (*dist_coeffs_29)[int(5)];
        _S3307 = _S3322;
        _S3308 = u_35;
        _S3309 = _S3321;
        _S3310 = _S3320;
        _S3311 = r2_35;
        _S3312 = _S3318;
        _S3313 = _S3317;
        _S3314 = _S3316;
        _S3315 = (*dist_coeffs_29)[int(3)];
    }
    else
    {
        _S3297 = 0.0f;
        _S3298 = 0.0f;
        _S3296 = _S3254;
        _S3299 = 0.0f;
        _S3300 = 0.0f;
        _S3301 = 0.0f;
        _S3302 = 0.0f;
        _S3303 = 0.0f;
        _S3304 = 0.0f;
        _S3305 = 0.0f;
        _S3306 = 0.0f;
        _S3307 = 0.0f;
        _S3308 = 0.0f;
        _S3309 = 0.0f;
        _S3310 = 0.0f;
        _S3311 = 0.0f;
        _S3312 = 0.0f;
        _S3313 = 0.0f;
        _S3314 = 0.0f;
        _S3315 = 0.0f;
    }
    float2  _S3328 = make_float2 (_S3250);
    float2  uv_32 = _S3249 / make_float2 (_S3250);
    float2  _S3329 = make_float2 (_S3250 * _S3250);
    if(_S3251)
    {
        _S3257 = true;
    }
    else
    {
        _S3257 = !_S3253._S3220;
    }
    float2  _S3330;
    if(_S3257)
    {
        _S3330 = make_float2 (-10000.0f);
    }
    bool _S3331 = !_S3257;
    float2  _S3332;
    float _S3333;
    float _S3334;
    float _S3335;
    float _S3336;
    float _S3337;
    float _S3338;
    float _S3339;
    float _S3340;
    float _S3341;
    float _S3342;
    float _S3343;
    float _S3344;
    float _S3345;
    float _S3346;
    float _S3347;
    float _S3348;
    float _S3349;
    float _S3350;
    float _S3351;
    if(_S3331)
    {
        float u_36 = uv_32.x;
        float v_36 = uv_32.y;
        float r2_36 = u_36 * u_36 + v_36 * v_36;
        float _S3352 = (*dist_coeffs_29)[int(2)] + r2_36 * (*dist_coeffs_29)[int(3)];
        float _S3353 = (*dist_coeffs_29)[int(1)] + r2_36 * _S3352;
        float _S3354 = (*dist_coeffs_29)[int(0)] + r2_36 * _S3353;
        float radial_5 = 1.0f + r2_36 * _S3354;
        float2  _S3355 = make_float2 (radial_5);
        float _S3356 = 2.0f * (*dist_coeffs_29)[int(4)];
        float _S3357 = _S3356 * u_36;
        float _S3358 = 2.0f * u_36;
        float _S3359 = 2.0f * (*dist_coeffs_29)[int(5)];
        float _S3360 = _S3359 * u_36;
        float _S3361 = 2.0f * v_36;
        float2  _S3362 = uv_32 * make_float2 (radial_5) + make_float2 (_S3357 * v_36 + (*dist_coeffs_29)[int(5)] * (r2_36 + _S3358 * u_36) + (*dist_coeffs_29)[int(6)] * r2_36, _S3360 * v_36 + (*dist_coeffs_29)[int(4)] * (r2_36 + _S3361 * v_36) + (*dist_coeffs_29)[int(7)] * r2_36);
        float2  _S3363 = _S3362 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3362.x + (*dist_coeffs_29)[int(9)] * _S3362.y, 0.0f);
        _S3330 = make_float2 (fx_23 * _S3363.x + cx_23, fy_23 * _S3363.y + cy_23);
        _S3333 = (*dist_coeffs_29)[int(9)];
        _S3334 = (*dist_coeffs_29)[int(8)];
        _S3332 = _S3355;
        _S3335 = (*dist_coeffs_29)[int(7)];
        _S3336 = (*dist_coeffs_29)[int(4)];
        _S3337 = _S3361;
        _S3338 = v_36;
        _S3339 = _S3360;
        _S3340 = _S3359;
        _S3341 = (*dist_coeffs_29)[int(6)];
        _S3342 = (*dist_coeffs_29)[int(5)];
        _S3343 = _S3358;
        _S3344 = u_36;
        _S3345 = _S3357;
        _S3346 = _S3356;
        _S3347 = r2_36;
        _S3348 = _S3354;
        _S3349 = _S3353;
        _S3350 = _S3352;
        _S3351 = (*dist_coeffs_29)[int(3)];
    }
    else
    {
        _S3333 = 0.0f;
        _S3334 = 0.0f;
        _S3332 = _S3254;
        _S3335 = 0.0f;
        _S3336 = 0.0f;
        _S3337 = 0.0f;
        _S3338 = 0.0f;
        _S3339 = 0.0f;
        _S3340 = 0.0f;
        _S3341 = 0.0f;
        _S3342 = 0.0f;
        _S3343 = 0.0f;
        _S3344 = 0.0f;
        _S3345 = 0.0f;
        _S3346 = 0.0f;
        _S3347 = 0.0f;
        _S3348 = 0.0f;
        _S3349 = 0.0f;
        _S3350 = 0.0f;
        _S3351 = 0.0f;
    }
    float2  e0_4 = _S3294 - _S3258;
    float2  e1_4 = _S3330 - _S3294;
    float2  e2_0 = _S3258 - _S3330;
    float _S3364 = e0_4.x;
    float _S3365 = e1_4.y;
    float _S3366 = e0_4.y;
    float _S3367 = e1_4.x;
    float _S3368 = _S3364 * _S3365 - _S3366 * _S3367;
    float _S3369 = 1.0f - hardness_4.y;
    float _S3370 = -1.0f / _S3369;
    float _S3371 = _S3369 * _S3369;
    float _S3372 = _S3258.x;
    float _S3373 = _S3294.x;
    float _S3374 = s_primal_ctx_max_0(_S3372, _S3373);
    float _S3375 = _S3330.x;
    float _S3376 = s_primal_ctx_min_0(_S3372, _S3373);
    float _S3377 = _S3258.y;
    float _S3378 = _S3294.y;
    float _S3379 = s_primal_ctx_max_0(_S3377, _S3378);
    float _S3380 = _S3330.y;
    float _S3381 = s_primal_ctx_min_0(_S3377, _S3378);
    float3  _S3382 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3383 = transpose_0(R_23);
    float3  _S3384 = mean_23 - - s_primal_ctx_mul_1(_S3383, t_22);
    float _S3385 = _S3384.x;
    float _S3386 = _S3384.y;
    float _S3387 = _S3384.z;
    float _S3388 = _S3385 * _S3385 + _S3386 * _S3386 + _S3387 * _S3387;
    float _S3389 = s_primal_ctx_sqrt_0(_S3388);
    float x_49 = _S3385 / _S3389;
    float3  _S3390 = make_float3 (x_49);
    float _S3391 = _S3389 * _S3389;
    float y_22 = _S3386 / _S3389;
    float z_19 = _S3387 / _S3389;
    float3  _S3392 = make_float3 (z_19);
    float _S3393 = - y_22;
    float3  _S3394 = make_float3 (_S3393);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_49 * x_49 - y_22 * y_22;
    float _S3395 = 2.0f * x_49;
    float fS1_19 = _S3395 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S3396 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_49;
    float3  _S3397 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3398 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3399 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3400 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3401 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3401;
    float3  _S3402 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_49;
    float3  _S3403 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3404 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3405 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3406 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_49 * fC1_19 - y_22 * fS1_19);
    float3  _S3407 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_49 * fS1_19 + y_22 * fC1_19);
    float3  _S3408 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3393) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_49) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3409 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3410 = make_float3 (0.0f);
    float3  _S3411 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3412 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3413 = make_float3 (_S3412);
    float3  _S3414 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3412);
    float3  _S3415 = _S3411 + _S3414 + make_float3 (0.5f);
    float3  _S3416 = _S3411 - _S3414 + make_float3 (0.5f);
    float3  _S3417 = vert1_c_4 - vert0_c_4;
    float3  _S3418 = vert2_c_4 - vert0_c_4;
    float3  _S3419 = s_primal_ctx_cross_0(_S3417, _S3418);
    float3  _S3420 = normalize_0(_S3419);
    float3  _S3421 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3420, mean_c_19)))))) * v_normal_0;
    float3  _S3422 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3423;
    (&_S3423)->primal_0 = _S3420;
    (&_S3423)->differential_0 = _S3422;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3424;
    (&_S3424)->primal_0 = mean_c_19;
    (&_S3424)->differential_0 = _S3422;
    s_bwd_prop_dot_0(&_S3423, &_S3424, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3425 = _S3424;
    float3  _S3426 = _S3421 + _S3423.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3427;
    (&_S3427)->primal_0 = _S3419;
    (&_S3427)->differential_0 = _S3422;
    s_bwd_normalize_impl_0(&_S3427, _S3426);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3428;
    (&_S3428)->primal_0 = _S3417;
    (&_S3428)->differential_0 = _S3422;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3429;
    (&_S3429)->primal_0 = _S3418;
    (&_S3429)->differential_0 = _S3422;
    s_bwd_prop_cross_0(&_S3428, &_S3429, _S3427.differential_0);
    float3  _S3430 = - _S3429.differential_0;
    float3  _S3431 = - _S3428.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3432;
    (&_S3432)->primal_0 = _S3416;
    (&_S3432)->differential_0 = _S3422;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3433;
    (&_S3433)->primal_0 = _S3410;
    (&_S3433)->differential_0 = _S3422;
    s_bwd_prop_max_0(&_S3432, &_S3433, (*v_rgb_6)[int(2)]);
    float3  _S3434 = - _S3432.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3435;
    (&_S3435)->primal_0 = _S3415;
    (&_S3435)->differential_0 = _S3422;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3436;
    (&_S3436)->primal_0 = _S3410;
    (&_S3436)->differential_0 = _S3422;
    s_bwd_prop_max_0(&_S3435, &_S3436, (*v_rgb_6)[int(1)]);
    float3  _S3437 = _S3413 * (_S3434 + _S3435.differential_0);
    float3  _S3438 = _S3432.differential_0 + _S3435.differential_0;
    float3  _S3439 = make_float3 (0.5f) * - _S3438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3440;
    (&_S3440)->primal_0 = _S3409;
    (&_S3440)->differential_0 = _S3422;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3441;
    (&_S3441)->primal_0 = _S3410;
    (&_S3441)->differential_0 = _S3422;
    s_bwd_prop_max_0(&_S3440, &_S3441, (*v_rgb_6)[int(0)]);
    float3  _S3442 = _S3439 + _S3440.differential_0;
    float3  _S3443 = _S3438 + _S3440.differential_0;
    float3  _S3444 = _S3407 * _S3443;
    float3  _S3445 = (*sh_coeffs_19)[int(15)] * _S3443;
    float3  _S3446 = _S3405 * _S3443;
    float3  _S3447 = (*sh_coeffs_19)[int(14)] * _S3443;
    float3  _S3448 = _S3403 * _S3443;
    float3  _S3449 = (*sh_coeffs_19)[int(13)] * _S3443;
    float3  _S3450 = _S3402 * _S3443;
    float3  _S3451 = (*sh_coeffs_19)[int(12)] * _S3443;
    float3  _S3452 = _S3404 * _S3443;
    float3  _S3453 = (*sh_coeffs_19)[int(11)] * _S3443;
    float3  _S3454 = _S3406 * _S3443;
    float3  _S3455 = (*sh_coeffs_19)[int(10)] * _S3443;
    float3  _S3456 = _S3408 * _S3443;
    float3  _S3457 = (*sh_coeffs_19)[int(9)] * _S3443;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3457.x + _S3457.y + _S3457.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3445.x + _S3445.y + _S3445.z);
    float _S3458 = _S3455.x + _S3455.y + _S3455.z;
    float _S3459 = _S3447.x + _S3447.y + _S3447.z;
    float _S3460 = _S3453.x + _S3453.y + _S3453.z;
    float _S3461 = _S3449.x + _S3449.y + _S3449.z;
    float _S3462 = _S3451.x + _S3451.y + _S3451.z;
    float _S3463 = - s_diff_fC2_T_5;
    float3  _S3464 = _S3399 * _S3443;
    float3  _S3465 = (*sh_coeffs_19)[int(8)] * _S3443;
    float3  _S3466 = _S3397 * _S3443;
    float3  _S3467 = (*sh_coeffs_19)[int(7)] * _S3443;
    float3  _S3468 = _S3396 * _S3443;
    float3  _S3469 = (*sh_coeffs_19)[int(6)] * _S3443;
    float3  _S3470 = _S3398 * _S3443;
    float3  _S3471 = (*sh_coeffs_19)[int(5)] * _S3443;
    float3  _S3472 = _S3400 * _S3443;
    float3  _S3473 = (*sh_coeffs_19)[int(4)] * _S3443;
    float _S3474 = _S3471.x + _S3471.y + _S3471.z;
    float _S3475 = _S3467.x + _S3467.y + _S3467.z;
    float _S3476 = fTmp1B_19 * _S3458 + x_49 * s_diff_fS2_T_5 + y_22 * _S3463 + 0.54627424478530884f * (_S3473.x + _S3473.y + _S3473.z);
    float _S3477 = fTmp1B_19 * _S3459 + y_22 * s_diff_fS2_T_5 + x_49 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3465.x + _S3465.y + _S3465.z);
    float _S3478 = y_22 * - _S3477;
    float _S3479 = x_49 * _S3477;
    float _S3480 = z_19 * (1.86588168144226074f * (z_19 * _S3462) + -2.28522896766662598f * (y_22 * _S3460 + x_49 * _S3461) + 0.94617468118667603f * (_S3469.x + _S3469.y + _S3469.z));
    float3  _S3481 = make_float3 (0.48860251903533936f) * _S3443;
    float3  _S3482 = - _S3481;
    float3  _S3483 = _S3390 * _S3482;
    float3  _S3484 = (*sh_coeffs_19)[int(3)] * _S3482;
    float3  _S3485 = _S3392 * _S3481;
    float3  _S3486 = (*sh_coeffs_19)[int(2)] * _S3481;
    float3  _S3487 = _S3394 * _S3481;
    float3  _S3488 = (*sh_coeffs_19)[int(1)] * _S3481;
    float _S3489 = (_S3401 * _S3462 + 1.44530570507049561f * (fS1_19 * _S3458 + fC1_19 * _S3459) + -1.09254848957061768f * (y_22 * _S3474 + x_49 * _S3475) + _S3480 + _S3480 + _S3486.x + _S3486.y + _S3486.z) / _S3391;
    float _S3490 = _S3389 * _S3489;
    float _S3491 = (fTmp0C_19 * _S3460 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3463 + fTmp0B_19 * _S3474 + _S3395 * _S3476 + _S3478 + _S3478 + - (_S3488.x + _S3488.y + _S3488.z)) / _S3391;
    float _S3492 = _S3389 * _S3491;
    float _S3493 = (fTmp0C_19 * _S3461 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3475 + 2.0f * (y_22 * _S3476) + _S3479 + _S3479 + _S3484.x + _S3484.y + _S3484.z) / _S3391;
    float _S3494 = _S3389 * _S3493;
    float _S3495 = _S3387 * - _S3489 + _S3386 * - _S3491 + _S3385 * - _S3493;
    DiffPair_float_0 _S3496;
    (&_S3496)->primal_0 = _S3388;
    (&_S3496)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3496, _S3495);
    float _S3497 = _S3387 * _S3496.differential_0;
    float _S3498 = _S3386 * _S3496.differential_0;
    float _S3499 = _S3385 * _S3496.differential_0;
    float3  _S3500 = make_float3 (0.282094806432724f) * _S3443;
    float3  _S3501 = make_float3 (_S3494 + _S3499 + _S3499, _S3492 + _S3498 + _S3498, _S3490 + _S3497 + _S3497);
    float3  _S3502 = - - _S3501;
    Matrix<float, 3, 3>  _S3503 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3504;
    (&_S3504)->primal_0 = _S3383;
    (&_S3504)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3505;
    (&_S3505)->primal_0 = t_22;
    (&_S3505)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3504, &_S3505, _S3502);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3506 = _S3505;
    Matrix<float, 3, 3>  _S3507 = transpose_0(_S3504.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3508;
    (&_S3508)->primal_0 = _S3382;
    (&_S3508)->differential_0 = _S3422;
    s_bwd_prop_log_1(&_S3508, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3509;
    (&_S3509)->primal_0 = vert2_c_4;
    (&_S3509)->differential_0 = _S3422;
    s_bwd_length_impl_0(&_S3509, _S3508.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3510;
    (&_S3510)->primal_0 = vert1_c_4;
    (&_S3510)->differential_0 = _S3422;
    s_bwd_length_impl_0(&_S3510, _S3508.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3511;
    (&_S3511)->primal_0 = vert0_c_4;
    (&_S3511)->differential_0 = _S3422;
    s_bwd_length_impl_0(&_S3511, _S3508.differential_0.x);
    DiffPair_float_0 _S3512;
    (&_S3512)->primal_0 = _S3381;
    (&_S3512)->differential_0 = 0.0f;
    DiffPair_float_0 _S3513;
    (&_S3513)->primal_0 = _S3380;
    (&_S3513)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3512, &_S3513, 0.0f);
    DiffPair_float_0 _S3514;
    (&_S3514)->primal_0 = _S3377;
    (&_S3514)->differential_0 = 0.0f;
    DiffPair_float_0 _S3515;
    (&_S3515)->primal_0 = _S3378;
    (&_S3515)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3514, &_S3515, _S3512.differential_0);
    DiffPair_float_0 _S3516;
    (&_S3516)->primal_0 = _S3379;
    (&_S3516)->differential_0 = 0.0f;
    DiffPair_float_0 _S3517;
    (&_S3517)->primal_0 = _S3380;
    (&_S3517)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3516, &_S3517, 0.0f);
    float _S3518 = _S3513.differential_0 + _S3517.differential_0;
    DiffPair_float_0 _S3519;
    (&_S3519)->primal_0 = _S3377;
    (&_S3519)->differential_0 = 0.0f;
    DiffPair_float_0 _S3520;
    (&_S3520)->primal_0 = _S3378;
    (&_S3520)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3519, &_S3520, _S3516.differential_0);
    float _S3521 = _S3515.differential_0 + _S3520.differential_0;
    float _S3522 = _S3514.differential_0 + _S3519.differential_0;
    DiffPair_float_0 _S3523;
    (&_S3523)->primal_0 = _S3376;
    (&_S3523)->differential_0 = 0.0f;
    DiffPair_float_0 _S3524;
    (&_S3524)->primal_0 = _S3375;
    (&_S3524)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3523, &_S3524, 0.0f);
    DiffPair_float_0 _S3525;
    (&_S3525)->primal_0 = _S3372;
    (&_S3525)->differential_0 = 0.0f;
    DiffPair_float_0 _S3526;
    (&_S3526)->primal_0 = _S3373;
    (&_S3526)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3525, &_S3526, _S3523.differential_0);
    DiffPair_float_0 _S3527;
    (&_S3527)->primal_0 = _S3374;
    (&_S3527)->differential_0 = 0.0f;
    DiffPair_float_0 _S3528;
    (&_S3528)->primal_0 = _S3375;
    (&_S3528)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3527, &_S3528, 0.0f);
    float _S3529 = _S3524.differential_0 + _S3528.differential_0;
    DiffPair_float_0 _S3530;
    (&_S3530)->primal_0 = _S3372;
    (&_S3530)->differential_0 = 0.0f;
    DiffPair_float_0 _S3531;
    (&_S3531)->primal_0 = _S3373;
    (&_S3531)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3530, &_S3531, _S3527.differential_0);
    float _S3532 = _S3526.differential_0 + _S3531.differential_0;
    float _S3533 = _S3525.differential_0 + _S3530.differential_0;
    DiffPair_float_0 _S3534;
    (&_S3534)->primal_0 = _S3370;
    (&_S3534)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3534, 0.0f);
    float _S3535 = - (-1.0f * - (_S3534.differential_0 / _S3371));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3536;
    (&_S3536)->primal_0 = e2_0;
    (&_S3536)->differential_0 = _S3254;
    s_bwd_length_impl_1(&_S3536, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3537;
    (&_S3537)->primal_0 = e1_4;
    (&_S3537)->differential_0 = _S3254;
    s_bwd_length_impl_1(&_S3537, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3538;
    (&_S3538)->primal_0 = e0_4;
    (&_S3538)->differential_0 = _S3254;
    s_bwd_length_impl_1(&_S3538, -0.0f);
    DiffPair_float_0 _S3539;
    (&_S3539)->primal_0 = _S3368;
    (&_S3539)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3539, 0.0f);
    float _S3540 = - _S3539.differential_0;
    float2  _S3541 = _S3537.differential_0 + make_float2 (_S3366 * _S3540, _S3364 * _S3539.differential_0);
    float2  _S3542 = _S3538.differential_0 + make_float2 (_S3365 * _S3539.differential_0, _S3367 * _S3540);
    float2  _S3543 = v_uv0_0 + _S3536.differential_0 + - _S3542 + make_float2 (_S3533, _S3522);
    float3  _S3544 = _S3428.differential_0 + _S3510.differential_0;
    float2  _S3545 = v_out_hardness_0 + make_float2 (0.0f, _S3535);
    float2  _S3546 = v_uv2_0 + - _S3536.differential_0 + _S3541 + make_float2 (_S3529, _S3518);
    float2  _S3547 = v_uv1_0 + - _S3541 + _S3542 + make_float2 (_S3532, _S3521);
    float3  _S3548 = _S3430 + _S3431 + _S3511.differential_0;
    float3  _S3549 = _S3429.differential_0 + _S3509.differential_0;
    FixedArray<float3 , 2>  _S3550;
    _S3550[int(0)] = _S3422;
    _S3550[int(1)] = _S3422;
    _S3550[int(1)] = _S3437;
    _S3550[int(0)] = _S3442;
    float3  _S3551 = _S3550[int(0)];
    float3  _S3552 = _S3550[int(1)];
    FixedArray<float3 , 16>  _S3553;
    _S3553[int(0)] = _S3422;
    _S3553[int(1)] = _S3422;
    _S3553[int(2)] = _S3422;
    _S3553[int(3)] = _S3422;
    _S3553[int(4)] = _S3422;
    _S3553[int(5)] = _S3422;
    _S3553[int(6)] = _S3422;
    _S3553[int(7)] = _S3422;
    _S3553[int(8)] = _S3422;
    _S3553[int(9)] = _S3422;
    _S3553[int(10)] = _S3422;
    _S3553[int(11)] = _S3422;
    _S3553[int(12)] = _S3422;
    _S3553[int(13)] = _S3422;
    _S3553[int(14)] = _S3422;
    _S3553[int(15)] = _S3422;
    _S3553[int(7)] = _S3466;
    _S3553[int(0)] = _S3500;
    _S3553[int(1)] = _S3487;
    _S3553[int(2)] = _S3485;
    _S3553[int(3)] = _S3483;
    _S3553[int(4)] = _S3472;
    _S3553[int(5)] = _S3470;
    _S3553[int(6)] = _S3468;
    _S3553[int(15)] = _S3444;
    _S3553[int(8)] = _S3464;
    _S3553[int(9)] = _S3456;
    _S3553[int(10)] = _S3454;
    _S3553[int(11)] = _S3452;
    _S3553[int(12)] = _S3450;
    _S3553[int(13)] = _S3448;
    _S3553[int(14)] = _S3446;
    float3  _S3554 = _S3553[int(0)];
    float3  _S3555 = _S3553[int(1)];
    float3  _S3556 = _S3553[int(2)];
    float3  _S3557 = _S3553[int(3)];
    float3  _S3558 = _S3553[int(4)];
    float3  _S3559 = _S3553[int(5)];
    float3  _S3560 = _S3553[int(6)];
    float3  _S3561 = _S3553[int(7)];
    float3  _S3562 = _S3553[int(8)];
    float3  _S3563 = _S3553[int(9)];
    float3  _S3564 = _S3553[int(10)];
    float3  _S3565 = _S3553[int(11)];
    float3  _S3566 = _S3553[int(12)];
    float3  _S3567 = _S3553[int(13)];
    float3  _S3568 = _S3553[int(14)];
    float3  _S3569 = _S3553[int(15)];
    if(_S3331)
    {
        float _S3570 = fx_23 * _S3546.x;
        float2  _S3571 = make_float2 (_S3570, fy_23 * _S3546.y) + make_float2 (_S3334 * _S3570, _S3333 * _S3570);
        float2  _S3572 = uv_32 * _S3571;
        float _S3573 = _S3336 * _S3571.y;
        float _S3574 = _S3342 * _S3571.x;
        float _S3575 = _S3572.x + _S3572.y;
        float _S3576 = _S3347 * _S3575;
        float _S3577 = _S3347 * _S3576;
        float _S3578 = _S3335 * _S3571.y + _S3573 + _S3341 * _S3571.x + _S3574 + _S3348 * _S3575 + _S3349 * _S3576 + _S3350 * _S3577 + _S3351 * (_S3347 * _S3577);
        float _S3579 = _S3338 * _S3578;
        float _S3580 = _S3344 * _S3578;
        _S3258 = _S3332 * _S3571 + make_float2 (_S3340 * (_S3338 * _S3571.y) + _S3343 * _S3574 + 2.0f * (_S3344 * _S3574) + _S3346 * (_S3338 * _S3571.x) + _S3580 + _S3580, _S3337 * _S3573 + 2.0f * (_S3338 * _S3573) + _S3339 * _S3571.y + _S3345 * _S3571.x + _S3579 + _S3579);
    }
    else
    {
        _S3258 = _S3254;
    }
    float2  _S3581 = _S3258 / _S3329;
    float2  _S3582 = _S3249 * - _S3581;
    float2  _S3583 = _S3328 * _S3581;
    float3  _S3584 = _S3549 + make_float3 (_S3583.x, _S3583.y, _S3582.x + _S3582.y);
    if(_S3295)
    {
        float _S3585 = fx_23 * _S3547.x;
        float2  _S3586 = make_float2 (_S3585, fy_23 * _S3547.y) + make_float2 (_S3298 * _S3585, _S3297 * _S3585);
        float2  _S3587 = uv_31 * _S3586;
        float _S3588 = _S3300 * _S3586.y;
        float _S3589 = _S3306 * _S3586.x;
        float _S3590 = _S3587.x + _S3587.y;
        float _S3591 = _S3311 * _S3590;
        float _S3592 = _S3311 * _S3591;
        float _S3593 = _S3299 * _S3586.y + _S3588 + _S3305 * _S3586.x + _S3589 + _S3312 * _S3590 + _S3313 * _S3591 + _S3314 * _S3592 + _S3315 * (_S3311 * _S3592);
        float _S3594 = _S3302 * _S3593;
        float _S3595 = _S3308 * _S3593;
        _S3258 = _S3296 * _S3586 + make_float2 (_S3304 * (_S3302 * _S3586.y) + _S3307 * _S3589 + 2.0f * (_S3308 * _S3589) + _S3310 * (_S3302 * _S3586.x) + _S3595 + _S3595, _S3301 * _S3588 + 2.0f * (_S3302 * _S3588) + _S3303 * _S3586.y + _S3309 * _S3586.x + _S3594 + _S3594);
    }
    else
    {
        _S3258 = _S3254;
    }
    float2  _S3596 = _S3258 / _S3293;
    float2  _S3597 = _S3245 * - _S3596;
    float2  _S3598 = _S3292 * _S3596;
    float3  _S3599 = _S3544 + make_float3 (_S3598.x, _S3598.y, _S3597.x + _S3597.y);
    if(_S3259)
    {
        float _S3600 = fx_23 * _S3543.x;
        float2  _S3601 = make_float2 (_S3600, fy_23 * _S3543.y) + make_float2 (_S3262 * _S3600, _S3261 * _S3600);
        float2  _S3602 = uv_30 * _S3601;
        float _S3603 = _S3264 * _S3601.y;
        float _S3604 = _S3270 * _S3601.x;
        float _S3605 = _S3602.x + _S3602.y;
        float _S3606 = _S3275 * _S3605;
        float _S3607 = _S3275 * _S3606;
        float _S3608 = _S3263 * _S3601.y + _S3603 + _S3269 * _S3601.x + _S3604 + _S3276 * _S3605 + _S3277 * _S3606 + _S3278 * _S3607 + _S3279 * (_S3275 * _S3607);
        float _S3609 = _S3266 * _S3608;
        float _S3610 = _S3272 * _S3608;
        _S3258 = _S3260 * _S3601 + make_float2 (_S3268 * (_S3266 * _S3601.y) + _S3271 * _S3604 + 2.0f * (_S3272 * _S3604) + _S3274 * (_S3266 * _S3601.x) + _S3610 + _S3610, _S3265 * _S3603 + 2.0f * (_S3266 * _S3603) + _S3267 * _S3601.y + _S3273 * _S3601.x + _S3609 + _S3609);
    }
    else
    {
        _S3258 = _S3254;
    }
    float2  _S3611 = _S3258 / _S3256;
    float2  _S3612 = _S3241 * - _S3611;
    float2  _S3613 = _S3255 * _S3611;
    float _S3614 = _S3612.x + _S3612.y;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3615;
    (&_S3615)->primal_0 = R_23;
    (&_S3615)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3616;
    (&_S3616)->primal_0 = vert2_0;
    (&_S3616)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3615, &_S3616, _S3584);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3617;
    (&_S3617)->primal_0 = R_23;
    (&_S3617)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3618;
    (&_S3618)->primal_0 = vert1_0;
    (&_S3618)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3617, &_S3618, _S3599);
    float3  _S3619 = _S3548 + make_float3 (_S3613.x, _S3613.y, _S3614);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3620;
    (&_S3620)->primal_0 = R_23;
    (&_S3620)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3621;
    (&_S3621)->primal_0 = vert0_0;
    (&_S3621)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3620, &_S3621, _S3619);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3622;
    (&_S3622)->primal_0 = _S3235;
    (&_S3622)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3623;
    (&_S3623)->primal_0 = _S3240;
    (&_S3623)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3622, &_S3623, _S3616.differential_0);
    float _S3624 = - _S3623.differential_0.y;
    float _S3625 = _S3239 * _S3623.differential_0.x;
    float _S3626 = - (_S3231 * _S3623.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3627;
    (&_S3627)->primal_0 = _S3235;
    (&_S3627)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3628;
    (&_S3628)->primal_0 = _S3238;
    (&_S3628)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3627, &_S3628, _S3618.differential_0);
    float _S3629 = _S3231 * _S3628.differential_0.x;
    float _S3630 = _S3237 * _S3628.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3631;
    (&_S3631)->primal_0 = _S3235;
    (&_S3631)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3632;
    (&_S3632)->primal_0 = _S3236;
    (&_S3632)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3631, &_S3632, _S3621.differential_0);
    Matrix<float, 3, 3>  _S3633 = transpose_0(_S3622.differential_0 + _S3627.differential_0 + _S3631.differential_0);
    float _S3634 = 2.0f * - _S3633.rows[int(2)].z;
    float _S3635 = 2.0f * _S3633.rows[int(2)].y;
    float _S3636 = 2.0f * _S3633.rows[int(2)].x;
    float _S3637 = 2.0f * _S3633.rows[int(1)].z;
    float _S3638 = 2.0f * - _S3633.rows[int(1)].y;
    float _S3639 = 2.0f * _S3633.rows[int(1)].x;
    float _S3640 = 2.0f * _S3633.rows[int(0)].z;
    float _S3641 = 2.0f * _S3633.rows[int(0)].y;
    float _S3642 = 2.0f * - _S3633.rows[int(0)].x;
    float _S3643 = - _S3639 + _S3641;
    float _S3644 = _S3636 + - _S3640;
    float _S3645 = - _S3635 + _S3637;
    float _S3646 = _S3635 + _S3637;
    float _S3647 = _S3636 + _S3640;
    float _S3648 = _S3639 + _S3641;
    float _S3649 = quat_24.w * (_S3638 + _S3642);
    float _S3650 = quat_24.z * (_S3634 + _S3642);
    float _S3651 = quat_24.y * (_S3634 + _S3638);
    float _S3652 = quat_24.x * _S3643 + quat_24.z * _S3646 + quat_24.y * _S3647 + _S3649 + _S3649;
    float _S3653 = quat_24.x * _S3644 + quat_24.w * _S3646 + quat_24.y * _S3648 + _S3650 + _S3650;
    float _S3654 = quat_24.x * _S3645 + quat_24.w * _S3647 + quat_24.z * _S3648 + _S3651 + _S3651;
    float _S3655 = quat_24.w * _S3643 + quat_24.z * _S3644 + quat_24.y * _S3645;
    float _S3656 = _S3626 + _S3629;
    float _S3657 = 0.5f * - _S3656;
    float _S3658 = _S3624 + _S3628.differential_0.y;
    DiffPair_float_0 _S3659;
    (&_S3659)->primal_0 = _S3232;
    (&_S3659)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3659, _S3658);
    float _S3660 = _S3657 + _S3659.differential_0;
    float _S3661 = _S3625 + _S3630 + _S3632.differential_0.x;
    DiffPair_float_0 _S3662;
    (&_S3662)->primal_0 = _S3230;
    (&_S3662)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3662, _S3661);
    float _S3663 = _S3657 + _S3662.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3664;
    (&_S3664)->primal_0 = R_23;
    (&_S3664)->differential_0 = _S3503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3665;
    (&_S3665)->primal_0 = mean_23;
    (&_S3665)->differential_0 = _S3422;
    s_bwd_prop_mul_1(&_S3664, &_S3665, _S3425.differential_0);
    float3  _S3666 = _S3584 + _S3599 + _S3619 + _S3425.differential_0 + _S3506.differential_0;
    Matrix<float, 3, 3>  _S3667 = _S3615.differential_0 + _S3617.differential_0 + _S3620.differential_0 + _S3664.differential_0 + _S3507;
    float3  _S3668 = make_float3 (_S3663, _S3660, _S3656);
    float4  _S3669 = make_float4 (0.0f);
    *&((&_S3669)->w) = _S3652;
    *&((&_S3669)->z) = _S3653;
    *&((&_S3669)->y) = _S3654;
    *&((&_S3669)->x) = _S3655;
    *v_mean_7 = _S3616.differential_0 + _S3618.differential_0 + _S3621.differential_0 + _S3665.differential_0 + _S3501;
    *v_quat_6 = _S3669;
    *v_scale_6 = _S3668;
    *v_hardness_0 = _S3545;
    (*v_sh_coeffs_5)[int(0)] = _S3554;
    (*v_sh_coeffs_5)[int(1)] = _S3555;
    (*v_sh_coeffs_5)[int(2)] = _S3556;
    (*v_sh_coeffs_5)[int(3)] = _S3557;
    (*v_sh_coeffs_5)[int(4)] = _S3558;
    (*v_sh_coeffs_5)[int(5)] = _S3559;
    (*v_sh_coeffs_5)[int(6)] = _S3560;
    (*v_sh_coeffs_5)[int(7)] = _S3561;
    (*v_sh_coeffs_5)[int(8)] = _S3562;
    (*v_sh_coeffs_5)[int(9)] = _S3563;
    (*v_sh_coeffs_5)[int(10)] = _S3564;
    (*v_sh_coeffs_5)[int(11)] = _S3565;
    (*v_sh_coeffs_5)[int(12)] = _S3566;
    (*v_sh_coeffs_5)[int(13)] = _S3567;
    (*v_sh_coeffs_5)[int(14)] = _S3568;
    (*v_sh_coeffs_5)[int(15)] = _S3569;
    (*v_ch_coeffs_0)[int(0)] = _S3551;
    (*v_ch_coeffs_0)[int(1)] = _S3552;
    *v_R_6 = _S3667;
    *v_t_6 = _S3666;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_24, float fy_24, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S3670 = scale_24.x;
    float _S3671 = s_primal_ctx_exp_1(_S3670);
    float _S3672 = scale_24.y;
    float _S3673 = s_primal_ctx_exp_1(_S3672);
    float sz_5 = scale_24.z - 0.5f * (_S3670 + _S3672);
    float _S3674 = quat_25.y;
    float x2_25 = _S3674 * _S3674;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3675 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3676 = make_float3 (_S3671, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S3675, _S3676) + mean_24;
    float _S3677 = -0.5f + sz_5;
    float3  _S3678 = make_float3 (_S3671 * _S3677, _S3673, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S3675, _S3678) + mean_24;
    float _S3679 = -0.5f - sz_5;
    float3  _S3680 = make_float3 (_S3671 * _S3679, - _S3673, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S3675, _S3680) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    float _S3681 = length_1(vert0_c_5);
    float _S3682 = length_1(vert1_c_5);
    float _S3683 = length_1(vert2_c_5);
    float2  _S3684 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3685 = length_0(_S3684);
    float _S3686 = vert0_c_5.z;
    float _S3687 = s_primal_ctx_atan2_0(_S3685, _S3686);
    bool _S3688 = _S3687 < 0.00100000004749745f;
    float k_6;
    float _S3689;
    float _S3690;
    float _S3691;
    if(_S3688)
    {
        float _S3692 = 1.0f - _S3687 * _S3687 / 3.0f;
        float _S3693 = _S3686 * _S3686;
        k_6 = _S3692 / _S3686;
        _S3689 = 0.0f;
        _S3690 = _S3693;
        _S3691 = _S3692;
    }
    else
    {
        float _S3694 = _S3685 * _S3685;
        k_6 = _S3687 / _S3685;
        _S3689 = _S3694;
        _S3690 = 0.0f;
        _S3691 = 0.0f;
    }
    float2  _S3695 = make_float2 (k_6);
    float2  _S3696 = _S3684 * make_float2 (k_6);
    float u_37 = _S3696.x;
    float v_37 = _S3696.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float _S3697 = (*dist_coeffs_30)[int(2)] + r2_37 * (*dist_coeffs_30)[int(3)];
    float _S3698 = (*dist_coeffs_30)[int(1)] + r2_37 * _S3697;
    float _S3699 = (*dist_coeffs_30)[int(0)] + r2_37 * _S3698;
    float radial_6 = 1.0f + r2_37 * _S3699;
    float2  _S3700 = make_float2 (radial_6);
    float _S3701 = 2.0f * (*dist_coeffs_30)[int(4)];
    float _S3702 = _S3701 * u_37;
    float _S3703 = 2.0f * u_37;
    float _S3704 = 2.0f * (*dist_coeffs_30)[int(5)];
    float _S3705 = _S3704 * u_37;
    float _S3706 = 2.0f * v_37;
    float2  _S3707 = _S3696 * make_float2 (radial_6) + make_float2 (_S3702 * v_37 + (*dist_coeffs_30)[int(5)] * (r2_37 + _S3703 * u_37) + (*dist_coeffs_30)[int(6)] * r2_37, _S3705 * v_37 + (*dist_coeffs_30)[int(4)] * (r2_37 + _S3706 * v_37) + (*dist_coeffs_30)[int(7)] * r2_37);
    float2  _S3708 = _S3707 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3707.x + (*dist_coeffs_30)[int(9)] * _S3707.y, 0.0f);
    float _S3709 = fx_24 * _S3708.x + cx_24;
    float _S3710 = fy_24 * _S3708.y + cy_24;
    float2  _S3711 = make_float2 (_S3709, _S3710);
    float2  _S3712 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3713 = length_0(_S3712);
    float _S3714 = vert1_c_5.z;
    float _S3715 = s_primal_ctx_atan2_0(_S3713, _S3714);
    bool _S3716 = _S3715 < 0.00100000004749745f;
    float _S3717;
    float _S3718;
    float _S3719;
    if(_S3716)
    {
        float _S3720 = 1.0f - _S3715 * _S3715 / 3.0f;
        float _S3721 = _S3714 * _S3714;
        k_6 = _S3720 / _S3714;
        _S3717 = 0.0f;
        _S3718 = _S3721;
        _S3719 = _S3720;
    }
    else
    {
        float _S3722 = _S3713 * _S3713;
        k_6 = _S3715 / _S3713;
        _S3717 = _S3722;
        _S3718 = 0.0f;
        _S3719 = 0.0f;
    }
    float2  _S3723 = make_float2 (k_6);
    float2  _S3724 = _S3712 * make_float2 (k_6);
    float u_38 = _S3724.x;
    float v_38 = _S3724.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float _S3725 = (*dist_coeffs_30)[int(2)] + r2_38 * (*dist_coeffs_30)[int(3)];
    float _S3726 = (*dist_coeffs_30)[int(1)] + r2_38 * _S3725;
    float _S3727 = (*dist_coeffs_30)[int(0)] + r2_38 * _S3726;
    float radial_7 = 1.0f + r2_38 * _S3727;
    float2  _S3728 = make_float2 (radial_7);
    float _S3729 = _S3701 * u_38;
    float _S3730 = 2.0f * u_38;
    float _S3731 = _S3704 * u_38;
    float _S3732 = 2.0f * v_38;
    float2  _S3733 = _S3724 * make_float2 (radial_7) + make_float2 (_S3729 * v_38 + (*dist_coeffs_30)[int(5)] * (r2_38 + _S3730 * u_38) + (*dist_coeffs_30)[int(6)] * r2_38, _S3731 * v_38 + (*dist_coeffs_30)[int(4)] * (r2_38 + _S3732 * v_38) + (*dist_coeffs_30)[int(7)] * r2_38);
    float2  _S3734 = _S3733 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3733.x + (*dist_coeffs_30)[int(9)] * _S3733.y, 0.0f);
    float _S3735 = fx_24 * _S3734.x + cx_24;
    float _S3736 = fy_24 * _S3734.y + cy_24;
    float2  _S3737 = make_float2 (_S3735, _S3736);
    float2  _S3738 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3739 = length_0(_S3738);
    float _S3740 = vert2_c_5.z;
    float _S3741 = s_primal_ctx_atan2_0(_S3739, _S3740);
    bool _S3742 = _S3741 < 0.00100000004749745f;
    float _S3743;
    float _S3744;
    float _S3745;
    if(_S3742)
    {
        float _S3746 = 1.0f - _S3741 * _S3741 / 3.0f;
        float _S3747 = _S3740 * _S3740;
        k_6 = _S3746 / _S3740;
        _S3743 = 0.0f;
        _S3744 = _S3747;
        _S3745 = _S3746;
    }
    else
    {
        float _S3748 = _S3739 * _S3739;
        k_6 = _S3741 / _S3739;
        _S3743 = _S3748;
        _S3744 = 0.0f;
        _S3745 = 0.0f;
    }
    float2  _S3749 = make_float2 (k_6);
    float2  _S3750 = _S3738 * make_float2 (k_6);
    float u_39 = _S3750.x;
    float v_39 = _S3750.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float _S3751 = (*dist_coeffs_30)[int(2)] + r2_39 * (*dist_coeffs_30)[int(3)];
    float _S3752 = (*dist_coeffs_30)[int(1)] + r2_39 * _S3751;
    float _S3753 = (*dist_coeffs_30)[int(0)] + r2_39 * _S3752;
    float radial_8 = 1.0f + r2_39 * _S3753;
    float2  _S3754 = make_float2 (radial_8);
    float _S3755 = _S3701 * u_39;
    float _S3756 = 2.0f * u_39;
    float _S3757 = _S3704 * u_39;
    float _S3758 = 2.0f * v_39;
    float2  _S3759 = _S3750 * make_float2 (radial_8) + make_float2 (_S3755 * v_39 + (*dist_coeffs_30)[int(5)] * (r2_39 + _S3756 * u_39) + (*dist_coeffs_30)[int(6)] * r2_39, _S3757 * v_39 + (*dist_coeffs_30)[int(4)] * (r2_39 + _S3758 * v_39) + (*dist_coeffs_30)[int(7)] * r2_39);
    float2  _S3760 = _S3759 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3759.x + (*dist_coeffs_30)[int(9)] * _S3759.y, 0.0f);
    float _S3761 = fx_24 * _S3760.x + cx_24;
    float _S3762 = fy_24 * _S3760.y + cy_24;
    float2  _S3763 = make_float2 (_S3761, _S3762);
    float2  e0_5 = _S3737 - _S3711;
    float2  e1_5 = _S3763 - _S3737;
    float2  e2_1 = _S3711 - _S3763;
    float _S3764 = e0_5.x;
    float _S3765 = e1_5.y;
    float _S3766 = e0_5.y;
    float _S3767 = e1_5.x;
    float _S3768 = _S3764 * _S3765 - _S3766 * _S3767;
    float _S3769 = 1.0f - hardness_5.y;
    float _S3770 = -1.0f / _S3769;
    float _S3771 = _S3769 * _S3769;
    float _S3772 = s_primal_ctx_max_0(_S3709, _S3735);
    float _S3773 = s_primal_ctx_min_0(_S3709, _S3735);
    float _S3774 = s_primal_ctx_max_0(_S3710, _S3736);
    float _S3775 = s_primal_ctx_min_0(_S3710, _S3736);
    float3  _S3776 = make_float3 (_S3681, _S3682, _S3683) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3777 = transpose_0(R_24);
    float3  _S3778 = mean_24 - - s_primal_ctx_mul_1(_S3777, t_23);
    float _S3779 = _S3778.x;
    float _S3780 = _S3778.y;
    float _S3781 = _S3778.z;
    float _S3782 = _S3779 * _S3779 + _S3780 * _S3780 + _S3781 * _S3781;
    float _S3783 = s_primal_ctx_sqrt_0(_S3782);
    float x_50 = _S3779 / _S3783;
    float3  _S3784 = make_float3 (x_50);
    float _S3785 = _S3783 * _S3783;
    float y_23 = _S3780 / _S3783;
    float z_20 = _S3781 / _S3783;
    float3  _S3786 = make_float3 (z_20);
    float _S3787 = - y_23;
    float3  _S3788 = make_float3 (_S3787);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_50 * x_50 - y_23 * y_23;
    float _S3789 = 2.0f * x_50;
    float fS1_20 = _S3789 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3790 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_50;
    float3  _S3791 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3792 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3793 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3794 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3795 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3795;
    float3  _S3796 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_50;
    float3  _S3797 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3798 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3799 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3800 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_50 * fC1_20 - y_23 * fS1_20);
    float3  _S3801 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_50 * fS1_20 + y_23 * fC1_20);
    float3  _S3802 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3787) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_50) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3803 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3804 = make_float3 (0.0f);
    float3  _S3805 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3806 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3807 = make_float3 (_S3806);
    float3  _S3808 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3806);
    float3  _S3809 = _S3805 + _S3808 + make_float3 (0.5f);
    float3  _S3810 = _S3805 - _S3808 + make_float3 (0.5f);
    float3  _S3811 = vert1_c_5 - vert0_c_5;
    float3  _S3812 = vert2_c_5 - vert0_c_5;
    float3  _S3813 = s_primal_ctx_cross_0(_S3811, _S3812);
    float3  _S3814 = normalize_0(_S3813);
    float3  _S3815 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3814, mean_c_20)))))) * v_normal_1;
    float3  _S3816 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3817;
    (&_S3817)->primal_0 = _S3814;
    (&_S3817)->differential_0 = _S3816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3818;
    (&_S3818)->primal_0 = mean_c_20;
    (&_S3818)->differential_0 = _S3816;
    s_bwd_prop_dot_0(&_S3817, &_S3818, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3819 = _S3818;
    float3  _S3820 = _S3815 + _S3817.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3821;
    (&_S3821)->primal_0 = _S3813;
    (&_S3821)->differential_0 = _S3816;
    s_bwd_normalize_impl_0(&_S3821, _S3820);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3822;
    (&_S3822)->primal_0 = _S3811;
    (&_S3822)->differential_0 = _S3816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3823;
    (&_S3823)->primal_0 = _S3812;
    (&_S3823)->differential_0 = _S3816;
    s_bwd_prop_cross_0(&_S3822, &_S3823, _S3821.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3824 = _S3822;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3825 = _S3823;
    float3  _S3826 = - _S3823.differential_0;
    float3  _S3827 = - _S3822.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3828;
    (&_S3828)->primal_0 = _S3810;
    (&_S3828)->differential_0 = _S3816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3829;
    (&_S3829)->primal_0 = _S3804;
    (&_S3829)->differential_0 = _S3816;
    s_bwd_prop_max_0(&_S3828, &_S3829, (*v_rgb_7)[int(2)]);
    float3  _S3830 = - _S3828.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3831;
    (&_S3831)->primal_0 = _S3809;
    (&_S3831)->differential_0 = _S3816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3832;
    (&_S3832)->primal_0 = _S3804;
    (&_S3832)->differential_0 = _S3816;
    s_bwd_prop_max_0(&_S3831, &_S3832, (*v_rgb_7)[int(1)]);
    float3  _S3833 = _S3807 * (_S3830 + _S3831.differential_0);
    float3  _S3834 = _S3828.differential_0 + _S3831.differential_0;
    float3  _S3835 = make_float3 (0.5f) * - _S3834;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3836;
    (&_S3836)->primal_0 = _S3803;
    (&_S3836)->differential_0 = _S3816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3837;
    (&_S3837)->primal_0 = _S3804;
    (&_S3837)->differential_0 = _S3816;
    s_bwd_prop_max_0(&_S3836, &_S3837, (*v_rgb_7)[int(0)]);
    float3  _S3838 = _S3835 + _S3836.differential_0;
    float3  _S3839 = _S3834 + _S3836.differential_0;
    float3  _S3840 = _S3801 * _S3839;
    float3  _S3841 = (*sh_coeffs_20)[int(15)] * _S3839;
    float3  _S3842 = _S3799 * _S3839;
    float3  _S3843 = (*sh_coeffs_20)[int(14)] * _S3839;
    float3  _S3844 = _S3797 * _S3839;
    float3  _S3845 = (*sh_coeffs_20)[int(13)] * _S3839;
    float3  _S3846 = _S3796 * _S3839;
    float3  _S3847 = (*sh_coeffs_20)[int(12)] * _S3839;
    float3  _S3848 = _S3798 * _S3839;
    float3  _S3849 = (*sh_coeffs_20)[int(11)] * _S3839;
    float3  _S3850 = _S3800 * _S3839;
    float3  _S3851 = (*sh_coeffs_20)[int(10)] * _S3839;
    float3  _S3852 = _S3802 * _S3839;
    float3  _S3853 = (*sh_coeffs_20)[int(9)] * _S3839;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3853.x + _S3853.y + _S3853.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3841.x + _S3841.y + _S3841.z);
    float _S3854 = _S3851.x + _S3851.y + _S3851.z;
    float _S3855 = _S3843.x + _S3843.y + _S3843.z;
    float _S3856 = _S3849.x + _S3849.y + _S3849.z;
    float _S3857 = _S3845.x + _S3845.y + _S3845.z;
    float _S3858 = _S3847.x + _S3847.y + _S3847.z;
    float _S3859 = - s_diff_fC2_T_6;
    float3  _S3860 = _S3793 * _S3839;
    float3  _S3861 = (*sh_coeffs_20)[int(8)] * _S3839;
    float3  _S3862 = _S3791 * _S3839;
    float3  _S3863 = (*sh_coeffs_20)[int(7)] * _S3839;
    float3  _S3864 = _S3790 * _S3839;
    float3  _S3865 = (*sh_coeffs_20)[int(6)] * _S3839;
    float3  _S3866 = _S3792 * _S3839;
    float3  _S3867 = (*sh_coeffs_20)[int(5)] * _S3839;
    float3  _S3868 = _S3794 * _S3839;
    float3  _S3869 = (*sh_coeffs_20)[int(4)] * _S3839;
    float _S3870 = _S3867.x + _S3867.y + _S3867.z;
    float _S3871 = _S3863.x + _S3863.y + _S3863.z;
    float _S3872 = fTmp1B_20 * _S3854 + x_50 * s_diff_fS2_T_6 + y_23 * _S3859 + 0.54627424478530884f * (_S3869.x + _S3869.y + _S3869.z);
    float _S3873 = fTmp1B_20 * _S3855 + y_23 * s_diff_fS2_T_6 + x_50 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3861.x + _S3861.y + _S3861.z);
    float _S3874 = y_23 * - _S3873;
    float _S3875 = x_50 * _S3873;
    float _S3876 = z_20 * (1.86588168144226074f * (z_20 * _S3858) + -2.28522896766662598f * (y_23 * _S3856 + x_50 * _S3857) + 0.94617468118667603f * (_S3865.x + _S3865.y + _S3865.z));
    float3  _S3877 = make_float3 (0.48860251903533936f) * _S3839;
    float3  _S3878 = - _S3877;
    float3  _S3879 = _S3784 * _S3878;
    float3  _S3880 = (*sh_coeffs_20)[int(3)] * _S3878;
    float3  _S3881 = _S3786 * _S3877;
    float3  _S3882 = (*sh_coeffs_20)[int(2)] * _S3877;
    float3  _S3883 = _S3788 * _S3877;
    float3  _S3884 = (*sh_coeffs_20)[int(1)] * _S3877;
    float _S3885 = (_S3795 * _S3858 + 1.44530570507049561f * (fS1_20 * _S3854 + fC1_20 * _S3855) + -1.09254848957061768f * (y_23 * _S3870 + x_50 * _S3871) + _S3876 + _S3876 + _S3882.x + _S3882.y + _S3882.z) / _S3785;
    float _S3886 = _S3783 * _S3885;
    float _S3887 = (fTmp0C_20 * _S3856 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3859 + fTmp0B_20 * _S3870 + _S3789 * _S3872 + _S3874 + _S3874 + - (_S3884.x + _S3884.y + _S3884.z)) / _S3785;
    float _S3888 = _S3783 * _S3887;
    float _S3889 = (fTmp0C_20 * _S3857 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3871 + 2.0f * (y_23 * _S3872) + _S3875 + _S3875 + _S3880.x + _S3880.y + _S3880.z) / _S3785;
    float _S3890 = _S3783 * _S3889;
    float _S3891 = _S3781 * - _S3885 + _S3780 * - _S3887 + _S3779 * - _S3889;
    DiffPair_float_0 _S3892;
    (&_S3892)->primal_0 = _S3782;
    (&_S3892)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3892, _S3891);
    float _S3893 = _S3781 * _S3892.differential_0;
    float _S3894 = _S3780 * _S3892.differential_0;
    float _S3895 = _S3779 * _S3892.differential_0;
    float3  _S3896 = make_float3 (0.282094806432724f) * _S3839;
    float3  _S3897 = make_float3 (_S3890 + _S3895 + _S3895, _S3888 + _S3894 + _S3894, _S3886 + _S3893 + _S3893);
    float3  _S3898 = - - _S3897;
    Matrix<float, 3, 3>  _S3899 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3900;
    (&_S3900)->primal_0 = _S3777;
    (&_S3900)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3901;
    (&_S3901)->primal_0 = t_23;
    (&_S3901)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S3900, &_S3901, _S3898);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3902 = _S3901;
    Matrix<float, 3, 3>  _S3903 = transpose_0(_S3900.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3904;
    (&_S3904)->primal_0 = _S3776;
    (&_S3904)->differential_0 = _S3816;
    s_bwd_prop_log_1(&_S3904, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3905 = _S3904;
    DiffPair_float_0 _S3906;
    (&_S3906)->primal_0 = _S3775;
    (&_S3906)->differential_0 = 0.0f;
    DiffPair_float_0 _S3907;
    (&_S3907)->primal_0 = _S3762;
    (&_S3907)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3906, &_S3907, 0.0f);
    DiffPair_float_0 _S3908;
    (&_S3908)->primal_0 = _S3710;
    (&_S3908)->differential_0 = 0.0f;
    DiffPair_float_0 _S3909;
    (&_S3909)->primal_0 = _S3736;
    (&_S3909)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3908, &_S3909, _S3906.differential_0);
    DiffPair_float_0 _S3910;
    (&_S3910)->primal_0 = _S3774;
    (&_S3910)->differential_0 = 0.0f;
    DiffPair_float_0 _S3911;
    (&_S3911)->primal_0 = _S3762;
    (&_S3911)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3910, &_S3911, 0.0f);
    DiffPair_float_0 _S3912;
    (&_S3912)->primal_0 = _S3710;
    (&_S3912)->differential_0 = 0.0f;
    DiffPair_float_0 _S3913;
    (&_S3913)->primal_0 = _S3736;
    (&_S3913)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3912, &_S3913, _S3910.differential_0);
    DiffPair_float_0 _S3914;
    (&_S3914)->primal_0 = _S3773;
    (&_S3914)->differential_0 = 0.0f;
    DiffPair_float_0 _S3915;
    (&_S3915)->primal_0 = _S3761;
    (&_S3915)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3914, &_S3915, 0.0f);
    DiffPair_float_0 _S3916;
    (&_S3916)->primal_0 = _S3709;
    (&_S3916)->differential_0 = 0.0f;
    DiffPair_float_0 _S3917;
    (&_S3917)->primal_0 = _S3735;
    (&_S3917)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3916, &_S3917, _S3914.differential_0);
    DiffPair_float_0 _S3918;
    (&_S3918)->primal_0 = _S3772;
    (&_S3918)->differential_0 = 0.0f;
    DiffPair_float_0 _S3919;
    (&_S3919)->primal_0 = _S3761;
    (&_S3919)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3918, &_S3919, 0.0f);
    DiffPair_float_0 _S3920;
    (&_S3920)->primal_0 = _S3709;
    (&_S3920)->differential_0 = 0.0f;
    DiffPair_float_0 _S3921;
    (&_S3921)->primal_0 = _S3735;
    (&_S3921)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3920, &_S3921, _S3918.differential_0);
    DiffPair_float_0 _S3922;
    (&_S3922)->primal_0 = _S3770;
    (&_S3922)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3922, 0.0f);
    float _S3923 = - (-1.0f * - (_S3922.differential_0 / _S3771));
    float2  _S3924 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3925;
    (&_S3925)->primal_0 = e2_1;
    (&_S3925)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S3925, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3926;
    (&_S3926)->primal_0 = e1_5;
    (&_S3926)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S3926, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3927;
    (&_S3927)->primal_0 = e0_5;
    (&_S3927)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S3927, -0.0f);
    DiffPair_float_0 _S3928;
    (&_S3928)->primal_0 = _S3768;
    (&_S3928)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3928, 0.0f);
    float _S3929 = - _S3928.differential_0;
    float2  _S3930 = _S3926.differential_0 + make_float2 (_S3766 * _S3929, _S3764 * _S3928.differential_0);
    float2  _S3931 = _S3927.differential_0 + make_float2 (_S3765 * _S3928.differential_0, _S3767 * _S3929);
    float2  _S3932 = v_uv2_1 + - _S3925.differential_0 + _S3930;
    float _S3933 = fx_24 * (_S3915.differential_0 + _S3919.differential_0 + _S3932.x);
    float2  _S3934 = make_float2 (_S3933, fy_24 * (_S3907.differential_0 + _S3911.differential_0 + _S3932.y)) + make_float2 ((*dist_coeffs_30)[int(8)] * _S3933, (*dist_coeffs_30)[int(9)] * _S3933);
    float2  _S3935 = _S3750 * _S3934;
    float2  _S3936 = _S3754 * _S3934;
    float _S3937 = (*dist_coeffs_30)[int(4)] * _S3934.y;
    float _S3938 = (*dist_coeffs_30)[int(5)] * _S3934.x;
    float _S3939 = _S3935.x + _S3935.y;
    float _S3940 = r2_39 * _S3939;
    float _S3941 = r2_39 * _S3940;
    float _S3942 = (*dist_coeffs_30)[int(7)] * _S3934.y + _S3937 + (*dist_coeffs_30)[int(6)] * _S3934.x + _S3938 + _S3753 * _S3939 + _S3752 * _S3940 + _S3751 * _S3941 + (*dist_coeffs_30)[int(3)] * (r2_39 * _S3941);
    float _S3943 = v_39 * _S3942;
    float _S3944 = u_39 * _S3942;
    float _S3945 = _S3758 * _S3937 + 2.0f * (v_39 * _S3937) + _S3757 * _S3934.y + _S3755 * _S3934.x + _S3943 + _S3943;
    float _S3946 = _S3704 * (v_39 * _S3934.y) + _S3756 * _S3938 + 2.0f * (u_39 * _S3938) + _S3701 * (v_39 * _S3934.x) + _S3944 + _S3944;
    float2  _S3947 = v_uv0_1 + _S3925.differential_0 + - _S3931;
    float2  _S3948 = v_out_hardness_1 + make_float2 (0.0f, _S3923);
    float _S3949 = _S3917.differential_0 + _S3921.differential_0;
    float2  _S3950 = v_uv1_1 + - _S3930 + _S3931;
    float3  _S3951 = _S3826 + _S3827;
    FixedArray<float3 , 2>  _S3952;
    _S3952[int(0)] = _S3816;
    _S3952[int(1)] = _S3816;
    _S3952[int(1)] = _S3833;
    _S3952[int(0)] = _S3838;
    float3  _S3953 = _S3952[int(0)];
    float3  _S3954 = _S3952[int(1)];
    FixedArray<float3 , 16>  _S3955;
    _S3955[int(0)] = _S3816;
    _S3955[int(1)] = _S3816;
    _S3955[int(2)] = _S3816;
    _S3955[int(3)] = _S3816;
    _S3955[int(4)] = _S3816;
    _S3955[int(5)] = _S3816;
    _S3955[int(6)] = _S3816;
    _S3955[int(7)] = _S3816;
    _S3955[int(8)] = _S3816;
    _S3955[int(9)] = _S3816;
    _S3955[int(10)] = _S3816;
    _S3955[int(11)] = _S3816;
    _S3955[int(12)] = _S3816;
    _S3955[int(13)] = _S3816;
    _S3955[int(14)] = _S3816;
    _S3955[int(15)] = _S3816;
    _S3955[int(7)] = _S3862;
    _S3955[int(0)] = _S3896;
    _S3955[int(1)] = _S3883;
    _S3955[int(2)] = _S3881;
    _S3955[int(3)] = _S3879;
    _S3955[int(4)] = _S3868;
    _S3955[int(5)] = _S3866;
    _S3955[int(6)] = _S3864;
    _S3955[int(15)] = _S3840;
    _S3955[int(8)] = _S3860;
    _S3955[int(9)] = _S3852;
    _S3955[int(10)] = _S3850;
    _S3955[int(11)] = _S3848;
    _S3955[int(12)] = _S3846;
    _S3955[int(13)] = _S3844;
    _S3955[int(14)] = _S3842;
    float3  _S3956 = _S3955[int(0)];
    float3  _S3957 = _S3955[int(1)];
    float3  _S3958 = _S3955[int(2)];
    float3  _S3959 = _S3955[int(3)];
    float3  _S3960 = _S3955[int(4)];
    float3  _S3961 = _S3955[int(5)];
    float3  _S3962 = _S3955[int(6)];
    float3  _S3963 = _S3955[int(7)];
    float3  _S3964 = _S3955[int(8)];
    float3  _S3965 = _S3955[int(9)];
    float3  _S3966 = _S3955[int(10)];
    float3  _S3967 = _S3955[int(11)];
    float3  _S3968 = _S3955[int(12)];
    float3  _S3969 = _S3955[int(13)];
    float3  _S3970 = _S3955[int(14)];
    float3  _S3971 = _S3955[int(15)];
    float _S3972 = _S3916.differential_0 + _S3920.differential_0;
    float _S3973 = _S3908.differential_0 + _S3912.differential_0;
    float _S3974 = _S3909.differential_0 + _S3913.differential_0;
    float2  _S3975 = _S3936 + make_float2 (_S3946, _S3945);
    float2  _S3976 = _S3738 * _S3975;
    float2  _S3977 = _S3749 * _S3975;
    float _S3978 = _S3976.x + _S3976.y;
    if(_S3742)
    {
        float _S3979 = _S3978 / _S3744;
        float _S3980 = _S3745 * - _S3979;
        float _S3981 = _S3741 * (0.3333333432674408f * - (_S3740 * _S3979));
        k_6 = _S3981 + _S3981;
        _S3743 = _S3980;
        _S3744 = 0.0f;
    }
    else
    {
        float _S3982 = _S3978 / _S3743;
        float _S3983 = _S3741 * - _S3982;
        k_6 = _S3739 * _S3982;
        _S3743 = 0.0f;
        _S3744 = _S3983;
    }
    DiffPair_float_0 _S3984;
    (&_S3984)->primal_0 = _S3739;
    (&_S3984)->differential_0 = 0.0f;
    DiffPair_float_0 _S3985;
    (&_S3985)->primal_0 = _S3740;
    (&_S3985)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3984, &_S3985, k_6);
    float _S3986 = _S3985.differential_0 + _S3743;
    float _S3987 = _S3984.differential_0 + _S3744;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3988;
    (&_S3988)->primal_0 = _S3738;
    (&_S3988)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S3988, _S3987);
    float2  _S3989 = _S3988.differential_0 + _S3977;
    float _S3990 = fx_24 * (_S3950.x + _S3949);
    float2  _S3991 = make_float2 (_S3990, fy_24 * (_S3950.y + _S3974)) + make_float2 ((*dist_coeffs_30)[int(8)] * _S3990, (*dist_coeffs_30)[int(9)] * _S3990);
    float2  _S3992 = _S3724 * _S3991;
    float _S3993 = (*dist_coeffs_30)[int(4)] * _S3991.y;
    float _S3994 = (*dist_coeffs_30)[int(5)] * _S3991.x;
    float _S3995 = _S3992.x + _S3992.y;
    float _S3996 = r2_38 * _S3995;
    float _S3997 = r2_38 * _S3996;
    float _S3998 = (*dist_coeffs_30)[int(7)] * _S3991.y + _S3993 + (*dist_coeffs_30)[int(6)] * _S3991.x + _S3994 + _S3727 * _S3995 + _S3726 * _S3996 + _S3725 * _S3997 + (*dist_coeffs_30)[int(3)] * (r2_38 * _S3997);
    float _S3999 = v_38 * _S3998;
    float _S4000 = u_38 * _S3998;
    float3  _S4001 = _S3825.differential_0 + make_float3 (_S3989.x, _S3989.y, _S3986);
    float2  _S4002 = _S3728 * _S3991 + make_float2 (_S3704 * (v_38 * _S3991.y) + _S3730 * _S3994 + 2.0f * (u_38 * _S3994) + _S3701 * (v_38 * _S3991.x) + _S4000 + _S4000, _S3732 * _S3993 + 2.0f * (v_38 * _S3993) + _S3731 * _S3991.y + _S3729 * _S3991.x + _S3999 + _S3999);
    float2  _S4003 = _S3712 * _S4002;
    float2  _S4004 = _S3723 * _S4002;
    float _S4005 = _S4003.x + _S4003.y;
    if(_S3716)
    {
        float _S4006 = _S4005 / _S3718;
        float _S4007 = _S3719 * - _S4006;
        float _S4008 = _S3715 * (0.3333333432674408f * - (_S3714 * _S4006));
        k_6 = _S4008 + _S4008;
        _S3717 = _S4007;
        _S3718 = 0.0f;
    }
    else
    {
        float _S4009 = _S4005 / _S3717;
        float _S4010 = _S3715 * - _S4009;
        k_6 = _S3713 * _S4009;
        _S3717 = 0.0f;
        _S3718 = _S4010;
    }
    DiffPair_float_0 _S4011;
    (&_S4011)->primal_0 = _S3713;
    (&_S4011)->differential_0 = 0.0f;
    DiffPair_float_0 _S4012;
    (&_S4012)->primal_0 = _S3714;
    (&_S4012)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4011, &_S4012, k_6);
    float _S4013 = _S4012.differential_0 + _S3717;
    float _S4014 = _S4011.differential_0 + _S3718;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4015;
    (&_S4015)->primal_0 = _S3712;
    (&_S4015)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S4015, _S4014);
    float2  _S4016 = _S4015.differential_0 + _S4004;
    float _S4017 = fx_24 * (_S3947.x + _S3972);
    float2  _S4018 = make_float2 (_S4017, fy_24 * (_S3947.y + _S3973)) + make_float2 ((*dist_coeffs_30)[int(8)] * _S4017, (*dist_coeffs_30)[int(9)] * _S4017);
    float2  _S4019 = _S3696 * _S4018;
    float _S4020 = (*dist_coeffs_30)[int(4)] * _S4018.y;
    float _S4021 = (*dist_coeffs_30)[int(5)] * _S4018.x;
    float _S4022 = _S4019.x + _S4019.y;
    float _S4023 = r2_37 * _S4022;
    float _S4024 = r2_37 * _S4023;
    float _S4025 = (*dist_coeffs_30)[int(7)] * _S4018.y + _S4020 + (*dist_coeffs_30)[int(6)] * _S4018.x + _S4021 + _S3699 * _S4022 + _S3698 * _S4023 + _S3697 * _S4024 + (*dist_coeffs_30)[int(3)] * (r2_37 * _S4024);
    float _S4026 = v_37 * _S4025;
    float _S4027 = u_37 * _S4025;
    float3  _S4028 = _S3824.differential_0 + make_float3 (_S4016.x, _S4016.y, _S4013);
    float2  _S4029 = _S3700 * _S4018 + make_float2 (_S3704 * (v_37 * _S4018.y) + _S3703 * _S4021 + 2.0f * (u_37 * _S4021) + _S3701 * (v_37 * _S4018.x) + _S4027 + _S4027, _S3706 * _S4020 + 2.0f * (v_37 * _S4020) + _S3705 * _S4018.y + _S3702 * _S4018.x + _S4026 + _S4026);
    float2  _S4030 = _S3684 * _S4029;
    float2  _S4031 = _S3695 * _S4029;
    float _S4032 = _S4030.x + _S4030.y;
    if(_S3688)
    {
        float _S4033 = _S4032 / _S3690;
        float _S4034 = _S3691 * - _S4033;
        float _S4035 = _S3687 * (0.3333333432674408f * - (_S3686 * _S4033));
        k_6 = _S4035 + _S4035;
        _S3689 = _S4034;
        _S3690 = 0.0f;
    }
    else
    {
        float _S4036 = _S4032 / _S3689;
        float _S4037 = _S3687 * - _S4036;
        k_6 = _S3685 * _S4036;
        _S3689 = 0.0f;
        _S3690 = _S4037;
    }
    DiffPair_float_0 _S4038;
    (&_S4038)->primal_0 = _S3685;
    (&_S4038)->differential_0 = 0.0f;
    DiffPair_float_0 _S4039;
    (&_S4039)->primal_0 = _S3686;
    (&_S4039)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4038, &_S4039, k_6);
    float _S4040 = _S4039.differential_0 + _S3689;
    float _S4041 = _S4038.differential_0 + _S3690;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4042;
    (&_S4042)->primal_0 = _S3684;
    (&_S4042)->differential_0 = _S3924;
    s_bwd_length_impl_1(&_S4042, _S4041);
    float2  _S4043 = _S4042.differential_0 + _S4031;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4044;
    (&_S4044)->primal_0 = vert2_c_5;
    (&_S4044)->differential_0 = _S3816;
    s_bwd_length_impl_0(&_S4044, _S3905.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4045;
    (&_S4045)->primal_0 = vert1_c_5;
    (&_S4045)->differential_0 = _S3816;
    s_bwd_length_impl_0(&_S4045, _S3905.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4046;
    (&_S4046)->primal_0 = vert0_c_5;
    (&_S4046)->differential_0 = _S3816;
    s_bwd_length_impl_0(&_S4046, _S3905.differential_0.x);
    float3  _S4047 = _S4044.differential_0 + _S4001;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4048;
    (&_S4048)->primal_0 = R_24;
    (&_S4048)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4049;
    (&_S4049)->primal_0 = vert2_1;
    (&_S4049)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4048, &_S4049, _S4047);
    float3  _S4050 = _S4045.differential_0 + _S4028;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4051;
    (&_S4051)->primal_0 = R_24;
    (&_S4051)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4052;
    (&_S4052)->primal_0 = vert1_1;
    (&_S4052)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4051, &_S4052, _S4050);
    float3  _S4053 = _S4046.differential_0 + _S3951 + make_float3 (_S4043.x, _S4043.y, _S4040);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4054;
    (&_S4054)->primal_0 = R_24;
    (&_S4054)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4055;
    (&_S4055)->primal_0 = vert0_1;
    (&_S4055)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4054, &_S4055, _S4053);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4056;
    (&_S4056)->primal_0 = _S3675;
    (&_S4056)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4057;
    (&_S4057)->primal_0 = _S3680;
    (&_S4057)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4056, &_S4057, _S4049.differential_0);
    float _S4058 = - _S4057.differential_0.y;
    float _S4059 = _S3679 * _S4057.differential_0.x;
    float _S4060 = - (_S3671 * _S4057.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4061;
    (&_S4061)->primal_0 = _S3675;
    (&_S4061)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4062;
    (&_S4062)->primal_0 = _S3678;
    (&_S4062)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4061, &_S4062, _S4052.differential_0);
    float _S4063 = _S3671 * _S4062.differential_0.x;
    float _S4064 = _S3677 * _S4062.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4065;
    (&_S4065)->primal_0 = _S3675;
    (&_S4065)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4066;
    (&_S4066)->primal_0 = _S3676;
    (&_S4066)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4065, &_S4066, _S4055.differential_0);
    Matrix<float, 3, 3>  _S4067 = transpose_0(_S4056.differential_0 + _S4061.differential_0 + _S4065.differential_0);
    float _S4068 = 2.0f * - _S4067.rows[int(2)].z;
    float _S4069 = 2.0f * _S4067.rows[int(2)].y;
    float _S4070 = 2.0f * _S4067.rows[int(2)].x;
    float _S4071 = 2.0f * _S4067.rows[int(1)].z;
    float _S4072 = 2.0f * - _S4067.rows[int(1)].y;
    float _S4073 = 2.0f * _S4067.rows[int(1)].x;
    float _S4074 = 2.0f * _S4067.rows[int(0)].z;
    float _S4075 = 2.0f * _S4067.rows[int(0)].y;
    float _S4076 = 2.0f * - _S4067.rows[int(0)].x;
    float _S4077 = - _S4073 + _S4075;
    float _S4078 = _S4070 + - _S4074;
    float _S4079 = - _S4069 + _S4071;
    float _S4080 = _S4069 + _S4071;
    float _S4081 = _S4070 + _S4074;
    float _S4082 = _S4073 + _S4075;
    float _S4083 = quat_25.w * (_S4072 + _S4076);
    float _S4084 = quat_25.z * (_S4068 + _S4076);
    float _S4085 = quat_25.y * (_S4068 + _S4072);
    float _S4086 = quat_25.x * _S4077 + quat_25.z * _S4080 + quat_25.y * _S4081 + _S4083 + _S4083;
    float _S4087 = quat_25.x * _S4078 + quat_25.w * _S4080 + quat_25.y * _S4082 + _S4084 + _S4084;
    float _S4088 = quat_25.x * _S4079 + quat_25.w * _S4081 + quat_25.z * _S4082 + _S4085 + _S4085;
    float _S4089 = quat_25.w * _S4077 + quat_25.z * _S4078 + quat_25.y * _S4079;
    float _S4090 = _S4060 + _S4063;
    float _S4091 = 0.5f * - _S4090;
    float _S4092 = _S4058 + _S4062.differential_0.y;
    DiffPair_float_0 _S4093;
    (&_S4093)->primal_0 = _S3672;
    (&_S4093)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4093, _S4092);
    float _S4094 = _S4091 + _S4093.differential_0;
    float _S4095 = _S4059 + _S4064 + _S4066.differential_0.x;
    DiffPair_float_0 _S4096;
    (&_S4096)->primal_0 = _S3670;
    (&_S4096)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4096, _S4095);
    float _S4097 = _S4091 + _S4096.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4098;
    (&_S4098)->primal_0 = mean_c_20;
    (&_S4098)->differential_0 = _S3816;
    s_bwd_length_impl_0(&_S4098, 0.0f);
    float3  _S4099 = _S4098.differential_0 + _S3819.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4100;
    (&_S4100)->primal_0 = R_24;
    (&_S4100)->differential_0 = _S3899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4101;
    (&_S4101)->primal_0 = mean_24;
    (&_S4101)->differential_0 = _S3816;
    s_bwd_prop_mul_1(&_S4100, &_S4101, _S4099);
    float3  _S4102 = _S4047 + _S4050 + _S4053 + _S4099 + _S3902.differential_0;
    Matrix<float, 3, 3>  _S4103 = _S4048.differential_0 + _S4051.differential_0 + _S4054.differential_0 + _S4100.differential_0 + _S3903;
    float3  _S4104 = make_float3 (_S4097, _S4094, _S4090);
    float4  _S4105 = make_float4 (0.0f);
    *&((&_S4105)->w) = _S4086;
    *&((&_S4105)->z) = _S4087;
    *&((&_S4105)->y) = _S4088;
    *&((&_S4105)->x) = _S4089;
    float4  _S4106 = _S4105;
    float3  _S4107 = _S4049.differential_0 + _S4052.differential_0 + _S4055.differential_0 + _S4101.differential_0 + _S3897;
    *v_mean_8 = _S4107;
    *v_quat_7 = _S4106;
    *v_scale_7 = _S4104;
    *v_hardness_1 = _S3948;
    (*v_sh_coeffs_6)[int(0)] = _S3956;
    (*v_sh_coeffs_6)[int(1)] = _S3957;
    (*v_sh_coeffs_6)[int(2)] = _S3958;
    (*v_sh_coeffs_6)[int(3)] = _S3959;
    (*v_sh_coeffs_6)[int(4)] = _S3960;
    (*v_sh_coeffs_6)[int(5)] = _S3961;
    (*v_sh_coeffs_6)[int(6)] = _S3962;
    (*v_sh_coeffs_6)[int(7)] = _S3963;
    (*v_sh_coeffs_6)[int(8)] = _S3964;
    (*v_sh_coeffs_6)[int(9)] = _S3965;
    (*v_sh_coeffs_6)[int(10)] = _S3966;
    (*v_sh_coeffs_6)[int(11)] = _S3967;
    (*v_sh_coeffs_6)[int(12)] = _S3968;
    (*v_sh_coeffs_6)[int(13)] = _S3969;
    (*v_sh_coeffs_6)[int(14)] = _S3970;
    (*v_sh_coeffs_6)[int(15)] = _S3971;
    (*v_ch_coeffs_1)[int(0)] = _S3953;
    (*v_ch_coeffs_1)[int(1)] = _S3954;
    *v_R_7 = _S4103;
    *v_t_7 = _S4102;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_19)
{
    DiffPair_float_0 _S4108 = *dpx_16;
    bool _S4109;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S4109 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S4109 = false;
    }
    float _S4110;
    if(_S4109)
    {
        _S4110 = dOut_19;
    }
    else
    {
        _S4110 = 0.0f;
    }
    dpx_16->primal_0 = _S4108.primal_0;
    dpx_16->differential_0 = _S4110;
    DiffPair_float_0 _S4111 = *dpMin_0;
    if((_S4108.primal_0) < ((*dpMin_0).primal_0))
    {
        _S4110 = dOut_19;
    }
    else
    {
        _S4110 = 0.0f;
    }
    dpMin_0->primal_0 = _S4111.primal_0;
    dpMin_0->differential_0 = _S4110;
    DiffPair_float_0 _S4112 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S4110 = dOut_19;
    }
    else
    {
        _S4110 = 0.0f;
    }
    dpMax_0->primal_0 = _S4112.primal_0;
    dpMax_0->differential_0 = _S4110;
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
        DiffPair_float_0 _S4113 = *dpx_17;
        float _S4114 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S4114;
        float _S4115 = val_0 * (F32_log((_S4113.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S4115;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S4116 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S4116))));
    float2  _S4117 = p_0 - v0_0;
    float2  _S4118 = normalize_1(e0_6);
    float2  _S4119 = p_0 - v1_0;
    float2  _S4120 = normalize_1(e1_6);
    float2  _S4121 = p_0 - v2_0;
    float2  _S4122 = normalize_1(e2_2);
    float _S4123 = hardness_6.x;
    float _S4124 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S4117.x * _S4118.y - _S4117.y * _S4118.x)), (se_0 * (_S4119.x * _S4120.y - _S4119.y * _S4120.x))))), (se_0 * (_S4121.x * _S4122.y - _S4121.y * _S4122.x)))) / ((F32_abs((_S4116))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S4124))));
    float _S4125;
    if(a_1 <= 0.0f)
    {
        _S4125 = 0.0f;
    }
    else
    {
        _S4125 = (F32_min(((F32_pow((a_1), (_S4124)))), (0.99900001287460327f)));
    }
    return _S4123 * _S4125;
}

inline __device__ float s_primal_ctx_abs_0(float _S4126)
{
    return (F32_abs((_S4126)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S4127, float _S4128, float _S4129)
{
    return clamp_0(_S4127, _S4128, _S4129);
}

inline __device__ float s_primal_ctx_exp2_0(float _S4130)
{
    return (F32_exp2((_S4130)));
}

inline __device__ float s_primal_ctx_pow_0(float _S4131, float _S4132)
{
    return (F32_pow((_S4131), (_S4132)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S4133, DiffPair_float_0 * _S4134, float _S4135)
{
    _d_pow_0(_S4133, _S4134, _S4135);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S4136, DiffPair_float_0 * _S4137, DiffPair_float_0 * _S4138, float _S4139)
{
    _d_clamp_0(_S4136, _S4137, _S4138, _S4139);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_8)
{
    float _S4140 = length_0((*dpx_18).primal_0);
    float2  _S4141 = (*dpx_18).primal_0 * _s_dOut_8;
    float2  _S4142 = make_float2 (1.0f / _S4140) * _s_dOut_8;
    float _S4143 = - ((_S4141.x + _S4141.y) / (_S4140 * _S4140));
    float2  _S4144 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4145;
    (&_S4145)->primal_0 = (*dpx_18).primal_0;
    (&_S4145)->differential_0 = _S4144;
    s_bwd_length_impl_1(&_S4145, _S4143);
    float2  _S4146 = _S4142 + _S4145.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S4146;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4147, float2  _S4148)
{
    s_bwd_prop_normalize_impl_1(_S4147, _S4148);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_9)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S4149 = e0_7.x;
    float _S4150 = e1_7.y;
    float _S4151 = e0_7.y;
    float _S4152 = e1_7.x;
    float _S4153 = _S4149 * _S4150 - _S4151 * _S4152;
    float se_1 = float((F32_sign((_S4153))));
    float2  _S4154 = p_1 - (*dpv0_0).primal_0;
    float2  _S4155 = normalize_1(e0_7);
    float _S4156 = _S4154.x;
    float _S4157 = _S4155.y;
    float _S4158 = _S4154.y;
    float _S4159 = _S4155.x;
    float de0_0 = se_1 * (_S4156 * _S4157 - _S4158 * _S4159);
    float2  _S4160 = p_1 - (*dpv1_0).primal_0;
    float2  _S4161 = normalize_1(e1_7);
    float _S4162 = _S4160.x;
    float _S4163 = _S4161.y;
    float _S4164 = _S4160.y;
    float _S4165 = _S4161.x;
    float de1_0 = se_1 * (_S4162 * _S4163 - _S4164 * _S4165);
    float2  _S4166 = p_1 - (*dpv2_0).primal_0;
    float2  _S4167 = normalize_1(e2_3);
    float _S4168 = _S4166.x;
    float _S4169 = _S4167.y;
    float _S4170 = _S4166.y;
    float _S4171 = _S4167.x;
    float de2_0 = se_1 * (_S4168 * _S4169 - _S4170 * _S4171);
    float _S4172 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4173 = s_primal_ctx_max_0(_S4172, de2_0);
    float _S4174 = s_primal_ctx_abs_0(_S4153);
    float _S4175 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4174 / _S4175;
    float _S4176 = _S4175 * _S4175;
    float _S4177 = (*dphardness_0).primal_0.x;
    float _S4178 = (*dphardness_0).primal_0.y;
    float _S4179 = dmax_0 * dmax_0;
    float _S4180 = 1.0f + _S4173 / dmax_0;
    float _S4181 = 1.0f - s_primal_ctx_clamp_0(_S4178, 0.00499999988824129f, 0.98000001907348633f);
    float _S4182 = -1.0f / _S4181;
    float _S4183 = _S4181 * _S4181;
    float _S4184 = 1.0f - s_primal_ctx_exp2_0(_S4182);
    float a_2 = 1.0f - _S4180 * _S4184;
    bool _S4185 = a_2 <= 0.0f;
    float _S4186;
    float _S4187;
    if(_S4185)
    {
        _S4186 = 0.0f;
        _S4187 = 0.0f;
    }
    else
    {
        float _S4188 = s_primal_ctx_pow_0(a_2, _S4181);
        _S4186 = s_primal_ctx_min_0(_S4188, 0.99900001287460327f);
        _S4187 = _S4188;
    }
    float _S4189 = _S4177 * _s_dOut_9;
    float _S4190 = _S4186 * _s_dOut_9;
    if(_S4185)
    {
        _S4186 = 0.0f;
        _S4187 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4191;
        (&_S4191)->primal_0 = _S4187;
        (&_S4191)->differential_0 = 0.0f;
        DiffPair_float_0 _S4192;
        (&_S4192)->primal_0 = 0.99900001287460327f;
        (&_S4192)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4191, &_S4192, _S4189);
        DiffPair_float_0 _S4193;
        (&_S4193)->primal_0 = a_2;
        (&_S4193)->differential_0 = 0.0f;
        DiffPair_float_0 _S4194;
        (&_S4194)->primal_0 = _S4181;
        (&_S4194)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4193, &_S4194, _S4191.differential_0);
        _S4186 = _S4193.differential_0;
        _S4187 = _S4194.differential_0;
    }
    float _S4195 = - _S4186;
    float _S4196 = _S4184 * _S4195;
    float _S4197 = - (_S4180 * _S4195);
    DiffPair_float_0 _S4198;
    (&_S4198)->primal_0 = _S4182;
    (&_S4198)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4198, _S4197);
    float _S4199 = - (-1.0f * - (_S4198.differential_0 / _S4183) + _S4187);
    float _S4200 = _S4196 / _S4179;
    float s_diff_dmax_T_0 = _S4173 * - _S4200;
    float _S4201 = dmax_0 * _S4200;
    DiffPair_float_0 _S4202;
    (&_S4202)->primal_0 = _S4178;
    (&_S4202)->differential_0 = 0.0f;
    DiffPair_float_0 _S4203;
    (&_S4203)->primal_0 = 0.00499999988824129f;
    (&_S4203)->differential_0 = 0.0f;
    DiffPair_float_0 _S4204;
    (&_S4204)->primal_0 = 0.98000001907348633f;
    (&_S4204)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4202, &_S4203, &_S4204, _S4199);
    float _S4205 = s_diff_dmax_T_0 / _S4176;
    float _S4206 = _S4174 * - _S4205;
    float _S4207 = _S4175 * _S4205;
    float2  _S4208 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4209;
    (&_S4209)->primal_0 = e2_3;
    (&_S4209)->differential_0 = _S4208;
    s_bwd_length_impl_1(&_S4209, _S4206);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4210;
    (&_S4210)->primal_0 = e1_7;
    (&_S4210)->differential_0 = _S4208;
    s_bwd_length_impl_1(&_S4210, _S4206);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4211;
    (&_S4211)->primal_0 = e0_7;
    (&_S4211)->differential_0 = _S4208;
    s_bwd_length_impl_1(&_S4211, _S4206);
    DiffPair_float_0 _S4212;
    (&_S4212)->primal_0 = _S4153;
    (&_S4212)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4212, _S4207);
    DiffPair_float_0 _S4213;
    (&_S4213)->primal_0 = _S4172;
    (&_S4213)->differential_0 = 0.0f;
    DiffPair_float_0 _S4214;
    (&_S4214)->primal_0 = de2_0;
    (&_S4214)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4213, &_S4214, _S4201);
    DiffPair_float_0 _S4215;
    (&_S4215)->primal_0 = de0_0;
    (&_S4215)->differential_0 = 0.0f;
    DiffPair_float_0 _S4216;
    (&_S4216)->primal_0 = de1_0;
    (&_S4216)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4215, &_S4216, _S4213.differential_0);
    float _S4217 = se_1 * _S4214.differential_0;
    float _S4218 = - _S4217;
    float _S4219 = _S4171 * _S4218;
    float _S4220 = _S4169 * _S4217;
    float2  _S4221 = make_float2 (_S4170 * _S4218, _S4168 * _S4217);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4222;
    (&_S4222)->primal_0 = e2_3;
    (&_S4222)->differential_0 = _S4208;
    s_bwd_normalize_impl_1(&_S4222, _S4221);
    float2  _S4223 = - make_float2 (_S4220, _S4219);
    float _S4224 = se_1 * _S4216.differential_0;
    float _S4225 = - _S4224;
    float _S4226 = _S4165 * _S4225;
    float _S4227 = _S4163 * _S4224;
    float2  _S4228 = make_float2 (_S4164 * _S4225, _S4162 * _S4224);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4229;
    (&_S4229)->primal_0 = e1_7;
    (&_S4229)->differential_0 = _S4208;
    s_bwd_normalize_impl_1(&_S4229, _S4228);
    float2  _S4230 = - make_float2 (_S4227, _S4226);
    float _S4231 = se_1 * _S4215.differential_0;
    float _S4232 = - _S4231;
    float _S4233 = _S4159 * _S4232;
    float _S4234 = _S4157 * _S4231;
    float2  _S4235 = make_float2 (_S4158 * _S4232, _S4156 * _S4231);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4236;
    (&_S4236)->primal_0 = e0_7;
    (&_S4236)->differential_0 = _S4208;
    s_bwd_normalize_impl_1(&_S4236, _S4235);
    float2  _S4237 = - make_float2 (_S4234, _S4233);
    float _S4238 = - _S4212.differential_0;
    float2  _S4239 = _S4209.differential_0 + _S4222.differential_0;
    float2  _S4240 = - _S4239;
    float2  _S4241 = _S4210.differential_0 + _S4229.differential_0 + make_float2 (_S4151 * _S4238, _S4149 * _S4212.differential_0);
    float2  _S4242 = - _S4241;
    float2  _S4243 = _S4211.differential_0 + _S4236.differential_0 + make_float2 (_S4150 * _S4212.differential_0, _S4152 * _S4238);
    float2  _S4244 = - _S4243;
    float2  _S4245 = make_float2 (_S4190, _S4202.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4245;
    float2  _S4246 = _S4223 + _S4240 + _S4241;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4246;
    float2  _S4247 = _S4230 + _S4242 + _S4243;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4247;
    float2  _S4248 = _S4237 + _S4239 + _S4244;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4248;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4249, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4250, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4251, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4252, float2  _S4253, float _S4254)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4249, _S4250, _S4251, _S4252, _S4253, _S4254);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4255 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4255;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4255;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4255;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4255;
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
    float2  _S4256 = p_3 - v0_2;
    float2  _S4257 = p_3 - v1_2;
    float2  _S4258 = p_3 - v2_2;
    float _S4259 = e0_8.x;
    float _S4260 = e1_8.y;
    float _S4261 = e0_8.y;
    float _S4262 = e1_8.x;
    float _S4263 = _S4259 * _S4260 - _S4261 * _S4262;
    float se_2 = float((F32_sign((_S4263))));
    float _S4264 = hardness_8.x;
    float _S4265 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4256.x * _S4261 - _S4256.y * _S4259)), (se_2 * (_S4257.x * _S4260 - _S4257.y * _S4262))))), (se_2 * (_S4258.x * e2_4.y - _S4258.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4256 - e0_8 * make_float2 (clamp_0(dot_1(_S4256, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4257 - e1_8 * make_float2 (clamp_0(dot_1(_S4257, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4258 - e2_4 * make_float2 (clamp_0(dot_1(_S4258, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4263))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4265))));
    float _S4266;
    if(a_3 <= 0.0f)
    {
        _S4266 = 0.0f;
    }
    else
    {
        _S4266 = (F32_min(((F32_pow((a_3), (_S4265)))), (0.99900001287460327f)));
    }
    return _S4264 * _S4266;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4267, float2  _S4268)
{
    return dot_1(_S4267, _S4268);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4269, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4270, float _S4271)
{
    _d_dot_1(_S4269, _S4270, _S4271);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_10)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4272 = p_4 - (*dpv0_1).primal_0;
    float _S4273 = s_primal_ctx_dot_1(_S4272, e0_9);
    float _S4274 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4275 = _S4273 / _S4274;
    float _S4276 = _S4274 * _S4274;
    float _S4277 = s_primal_ctx_clamp_0(_S4275, 0.0f, 1.0f);
    float2  _S4278 = make_float2 (_S4277);
    float2  _S4279 = _S4272 - e0_9 * make_float2 (_S4277);
    float _S4280 = length_0(_S4279);
    float2  _S4281 = p_4 - (*dpv1_1).primal_0;
    float _S4282 = s_primal_ctx_dot_1(_S4281, e1_9);
    float _S4283 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4284 = _S4282 / _S4283;
    float _S4285 = _S4283 * _S4283;
    float _S4286 = s_primal_ctx_clamp_0(_S4284, 0.0f, 1.0f);
    float2  _S4287 = make_float2 (_S4286);
    float2  _S4288 = _S4281 - e1_9 * make_float2 (_S4286);
    float _S4289 = length_0(_S4288);
    float2  _S4290 = p_4 - (*dpv2_1).primal_0;
    float _S4291 = s_primal_ctx_dot_1(_S4290, e2_5);
    float _S4292 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4293 = _S4291 / _S4292;
    float _S4294 = _S4292 * _S4292;
    float _S4295 = s_primal_ctx_clamp_0(_S4293, 0.0f, 1.0f);
    float2  _S4296 = make_float2 (_S4295);
    float2  _S4297 = _S4290 - e2_5 * make_float2 (_S4295);
    float _S4298 = length_0(_S4297);
    float _S4299 = e0_9.x;
    float _S4300 = e1_9.y;
    float _S4301 = e0_9.y;
    float _S4302 = e1_9.x;
    float _S4303 = _S4299 * _S4300 - _S4301 * _S4302;
    float se_3 = float((F32_sign((_S4303))));
    float _S4304 = _S4272.x;
    float _S4305 = _S4272.y;
    float s0_0 = se_3 * (_S4304 * _S4301 - _S4305 * _S4299);
    float _S4306 = _S4281.x;
    float _S4307 = _S4281.y;
    float s1_0 = se_3 * (_S4306 * _S4300 - _S4307 * _S4302);
    float _S4308 = _S4290.x;
    float _S4309 = e2_5.y;
    float _S4310 = _S4290.y;
    float _S4311 = e2_5.x;
    float s2_0 = se_3 * (_S4308 * _S4309 - _S4310 * _S4311);
    float _S4312 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4312, s2_0)))));
    float _S4313 = s_primal_ctx_min_0(_S4280, _S4289);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4313, _S4298);
    float _S4314 = s_primal_ctx_abs_0(_S4303);
    float _S4315 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4314 / _S4315;
    float _S4316 = _S4315 * _S4315;
    float _S4317 = (*dphardness_1).primal_0.x;
    float _S4318 = (*dphardness_1).primal_0.y;
    float _S4319 = dmax_1 * dmax_1;
    float _S4320 = 1.0f + dv_0 / dmax_1;
    float _S4321 = 1.0f - s_primal_ctx_clamp_0(_S4318, 0.00499999988824129f, 0.98000001907348633f);
    float _S4322 = -1.0f / _S4321;
    float _S4323 = _S4321 * _S4321;
    float _S4324 = 1.0f - s_primal_ctx_exp2_0(_S4322);
    float a_4 = 1.0f - _S4320 * _S4324;
    bool _S4325 = a_4 <= 0.0f;
    float _S4326;
    float _S4327;
    if(_S4325)
    {
        _S4326 = 0.0f;
        _S4327 = 0.0f;
    }
    else
    {
        float _S4328 = s_primal_ctx_pow_0(a_4, _S4321);
        _S4326 = s_primal_ctx_min_0(_S4328, 0.99900001287460327f);
        _S4327 = _S4328;
    }
    float _S4329 = _S4317 * _s_dOut_10;
    float _S4330 = _S4326 * _s_dOut_10;
    if(_S4325)
    {
        _S4326 = 0.0f;
        _S4327 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4331;
        (&_S4331)->primal_0 = _S4327;
        (&_S4331)->differential_0 = 0.0f;
        DiffPair_float_0 _S4332;
        (&_S4332)->primal_0 = 0.99900001287460327f;
        (&_S4332)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4331, &_S4332, _S4329);
        DiffPair_float_0 _S4333;
        (&_S4333)->primal_0 = a_4;
        (&_S4333)->differential_0 = 0.0f;
        DiffPair_float_0 _S4334;
        (&_S4334)->primal_0 = _S4321;
        (&_S4334)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4333, &_S4334, _S4331.differential_0);
        _S4326 = _S4333.differential_0;
        _S4327 = _S4334.differential_0;
    }
    float _S4335 = - _S4326;
    float _S4336 = _S4324 * _S4335;
    float _S4337 = - (_S4320 * _S4335);
    DiffPair_float_0 _S4338;
    (&_S4338)->primal_0 = _S4322;
    (&_S4338)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4338, _S4337);
    float _S4339 = - (-1.0f * - (_S4338.differential_0 / _S4323) + _S4327);
    float _S4340 = _S4336 / _S4319;
    float s_diff_dmax_T_1 = dv_0 * - _S4340;
    float s_diff_dv_T_0 = dmax_1 * _S4340;
    DiffPair_float_0 _S4341;
    (&_S4341)->primal_0 = _S4318;
    (&_S4341)->differential_0 = 0.0f;
    DiffPair_float_0 _S4342;
    (&_S4342)->primal_0 = 0.00499999988824129f;
    (&_S4342)->differential_0 = 0.0f;
    DiffPair_float_0 _S4343;
    (&_S4343)->primal_0 = 0.98000001907348633f;
    (&_S4343)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4341, &_S4342, &_S4343, _S4339);
    float _S4344 = s_diff_dmax_T_1 / _S4316;
    float _S4345 = _S4314 * - _S4344;
    float _S4346 = _S4315 * _S4344;
    float2  _S4347 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4348;
    (&_S4348)->primal_0 = e2_5;
    (&_S4348)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4348, _S4345);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4349;
    (&_S4349)->primal_0 = e1_9;
    (&_S4349)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4349, _S4345);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4350;
    (&_S4350)->primal_0 = e0_9;
    (&_S4350)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4350, _S4345);
    DiffPair_float_0 _S4351;
    (&_S4351)->primal_0 = _S4303;
    (&_S4351)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4351, _S4346);
    float _S4352 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4353;
    (&_S4353)->primal_0 = _S4313;
    (&_S4353)->differential_0 = 0.0f;
    DiffPair_float_0 _S4354;
    (&_S4354)->primal_0 = _S4298;
    (&_S4354)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4353, &_S4354, _S4352);
    DiffPair_float_0 _S4355;
    (&_S4355)->primal_0 = _S4280;
    (&_S4355)->differential_0 = 0.0f;
    DiffPair_float_0 _S4356;
    (&_S4356)->primal_0 = _S4289;
    (&_S4356)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4355, &_S4356, _S4353.differential_0);
    DiffPair_float_0 _S4357;
    (&_S4357)->primal_0 = _S4312;
    (&_S4357)->differential_0 = 0.0f;
    DiffPair_float_0 _S4358;
    (&_S4358)->primal_0 = s2_0;
    (&_S4358)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4357, &_S4358, 0.0f);
    DiffPair_float_0 _S4359;
    (&_S4359)->primal_0 = s0_0;
    (&_S4359)->differential_0 = 0.0f;
    DiffPair_float_0 _S4360;
    (&_S4360)->primal_0 = s1_0;
    (&_S4360)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4359, &_S4360, _S4357.differential_0);
    float _S4361 = se_3 * _S4358.differential_0;
    float _S4362 = - _S4361;
    float _S4363 = _S4310 * _S4362;
    float _S4364 = _S4311 * _S4362;
    float _S4365 = _S4308 * _S4361;
    float _S4366 = _S4309 * _S4361;
    float _S4367 = se_3 * _S4360.differential_0;
    float _S4368 = - _S4367;
    float _S4369 = _S4302 * _S4368;
    float _S4370 = _S4300 * _S4367;
    float _S4371 = se_3 * _S4359.differential_0;
    float _S4372 = - _S4371;
    float _S4373 = _S4299 * _S4372;
    float _S4374 = _S4301 * _S4371;
    float _S4375 = - _S4351.differential_0;
    float _S4376 = _S4307 * _S4368 + _S4301 * _S4375;
    float _S4377 = _S4304 * _S4371 + _S4302 * _S4375;
    float _S4378 = _S4306 * _S4367 + _S4299 * _S4351.differential_0;
    float _S4379 = _S4305 * _S4372 + _S4300 * _S4351.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4380;
    (&_S4380)->primal_0 = _S4297;
    (&_S4380)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4380, _S4354.differential_0);
    float2  _S4381 = - _S4380.differential_0;
    float2  _S4382 = e2_5 * _S4381;
    float2  _S4383 = _S4296 * _S4381;
    float _S4384 = _S4382.x + _S4382.y;
    DiffPair_float_0 _S4385;
    (&_S4385)->primal_0 = _S4293;
    (&_S4385)->differential_0 = 0.0f;
    DiffPair_float_0 _S4386;
    (&_S4386)->primal_0 = 0.0f;
    (&_S4386)->differential_0 = 0.0f;
    DiffPair_float_0 _S4387;
    (&_S4387)->primal_0 = 1.0f;
    (&_S4387)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4385, &_S4386, &_S4387, _S4384);
    float _S4388 = _S4385.differential_0 / _S4294;
    float _S4389 = _S4291 * - _S4388;
    float _S4390 = _S4292 * _S4388;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4391;
    (&_S4391)->primal_0 = e2_5;
    (&_S4391)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4392;
    (&_S4392)->primal_0 = e2_5;
    (&_S4392)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4391, &_S4392, _S4389);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4393;
    (&_S4393)->primal_0 = _S4290;
    (&_S4393)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4394;
    (&_S4394)->primal_0 = e2_5;
    (&_S4394)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4393, &_S4394, _S4390);
    float2  _S4395 = - (_S4380.differential_0 + _S4393.differential_0 + make_float2 (_S4366, _S4364));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4396;
    (&_S4396)->primal_0 = _S4288;
    (&_S4396)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4396, _S4356.differential_0);
    float2  _S4397 = - _S4396.differential_0;
    float2  _S4398 = e1_9 * _S4397;
    float2  _S4399 = _S4287 * _S4397;
    float _S4400 = _S4398.x + _S4398.y;
    DiffPair_float_0 _S4401;
    (&_S4401)->primal_0 = _S4284;
    (&_S4401)->differential_0 = 0.0f;
    DiffPair_float_0 _S4402;
    (&_S4402)->primal_0 = 0.0f;
    (&_S4402)->differential_0 = 0.0f;
    DiffPair_float_0 _S4403;
    (&_S4403)->primal_0 = 1.0f;
    (&_S4403)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4401, &_S4402, &_S4403, _S4400);
    float _S4404 = _S4401.differential_0 / _S4285;
    float _S4405 = _S4282 * - _S4404;
    float _S4406 = _S4283 * _S4404;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4407;
    (&_S4407)->primal_0 = e1_9;
    (&_S4407)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4408;
    (&_S4408)->primal_0 = e1_9;
    (&_S4408)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4407, &_S4408, _S4405);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4409;
    (&_S4409)->primal_0 = _S4281;
    (&_S4409)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4410;
    (&_S4410)->primal_0 = e1_9;
    (&_S4410)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4409, &_S4410, _S4406);
    float2  _S4411 = - (_S4396.differential_0 + _S4409.differential_0 + make_float2 (_S4370, _S4369));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4412;
    (&_S4412)->primal_0 = _S4279;
    (&_S4412)->differential_0 = _S4347;
    s_bwd_length_impl_1(&_S4412, _S4355.differential_0);
    float2  _S4413 = - _S4412.differential_0;
    float2  _S4414 = e0_9 * _S4413;
    float2  _S4415 = _S4278 * _S4413;
    float _S4416 = _S4414.x + _S4414.y;
    DiffPair_float_0 _S4417;
    (&_S4417)->primal_0 = _S4275;
    (&_S4417)->differential_0 = 0.0f;
    DiffPair_float_0 _S4418;
    (&_S4418)->primal_0 = 0.0f;
    (&_S4418)->differential_0 = 0.0f;
    DiffPair_float_0 _S4419;
    (&_S4419)->primal_0 = 1.0f;
    (&_S4419)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4417, &_S4418, &_S4419, _S4416);
    float _S4420 = _S4417.differential_0 / _S4276;
    float _S4421 = _S4273 * - _S4420;
    float _S4422 = _S4274 * _S4420;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4423;
    (&_S4423)->primal_0 = e0_9;
    (&_S4423)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4424;
    (&_S4424)->primal_0 = e0_9;
    (&_S4424)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4423, &_S4424, _S4421);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4425;
    (&_S4425)->primal_0 = _S4272;
    (&_S4425)->differential_0 = _S4347;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4426;
    (&_S4426)->primal_0 = e0_9;
    (&_S4426)->differential_0 = _S4347;
    s_bwd_prop_dot_1(&_S4425, &_S4426, _S4422);
    float2  _S4427 = - (_S4412.differential_0 + _S4425.differential_0 + make_float2 (_S4374, _S4373));
    float2  _S4428 = _S4348.differential_0 + _S4383 + _S4392.differential_0 + _S4391.differential_0 + _S4394.differential_0 + make_float2 (_S4363, _S4365);
    float2  _S4429 = - _S4428;
    float2  _S4430 = _S4349.differential_0 + _S4399 + _S4408.differential_0 + _S4407.differential_0 + _S4410.differential_0 + make_float2 (_S4376, _S4378);
    float2  _S4431 = - _S4430;
    float2  _S4432 = _S4350.differential_0 + _S4415 + _S4424.differential_0 + _S4423.differential_0 + _S4426.differential_0 + make_float2 (_S4379, _S4377);
    float2  _S4433 = - _S4432;
    float2  _S4434 = make_float2 (_S4330, _S4341.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4434;
    float2  _S4435 = _S4395 + _S4429 + _S4430;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4435;
    float2  _S4436 = _S4411 + _S4431 + _S4432;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4436;
    float2  _S4437 = _S4427 + _S4428 + _S4433;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4437;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4438, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4439, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4440, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4441, float2  _S4442, float _S4443)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4438, _S4439, _S4440, _S4441, _S4442, _S4443);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4444 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4444;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4444;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4444;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4444;
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
    float _S4445 = 0.3333333432674408f * dpdepth_1;
    float3  _S4446 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4447 = make_float3 (0.0f);
    float3  _S4448 = _S4447;
    *&((&_S4448)->z) = _S4445;
    *&((&_S4448)->y) = _S4445;
    *&((&_S4448)->x) = _S4445;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4448;
    FixedArray<float3 , 3>  _S4449;
    _S4449[int(0)] = _S4447;
    _S4449[int(1)] = _S4447;
    _S4449[int(2)] = _S4447;
    _S4449[int(2)] = _S4446;
    _S4449[int(1)] = _S4446;
    _S4449[int(0)] = _S4446;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4449;
    float2  _S4450 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4450;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4450;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4450;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4451, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4452, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4453, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4454, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4455, float2  _S4456, float3  _S4457, float _S4458)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4451, _S4452, _S4453, _S4454, _S4455, _S4456, _S4457, _S4458);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4459 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4459;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4459;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4459;
    float3  _S4460 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4461 = { _S4460, _S4460, _S4460 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4461;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4460;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_8, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_25, float fy_25, float cx_25, float cy_25, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_25) + t_24;
        float _S4462 = mean_c_21.z;
        bool _S4463;
        if(_S4462 < near_plane_14)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4462 > far_plane_14;
        }
        if(_S4463)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4464 = scale_25.x;
        float sx_4 = (F32_exp((_S4464)));
        float _S4465 = scale_25.y;
        float sy_4 = (F32_exp((_S4465)));
        float sz_6 = scale_25.z - 0.5f * (_S4464 + _S4465);
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
        Matrix<float, 3, 3>  _S4466 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
        float3  vert0_2 = mul_0(_S4466, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
        float3  vert1_2 = mul_0(_S4466, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
        float3  vert2_2 = mul_0(_S4466, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
        float3  vert0_c_6 = mul_0(R_25, vert0_2) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_2) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_2) + t_24;
        float _S4467 = vert0_c_6.z;
        float _S4468 = vert1_c_6.z;
        float _S4469 = vert2_c_6.z;
        if(_S4467 < near_plane_14)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4467 > far_plane_14;
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4468 < near_plane_14;
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4468 > far_plane_14;
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4469 < near_plane_14;
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = _S4469 > far_plane_14;
        }
        if(_S4463)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S4470;
        for(;;)
        {
            float2  uv_33 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4467);
            if(_S4467 < 0.0f)
            {
                _S4463 = true;
            }
            else
            {
                bool _S4471 = is_valid_distortion(uv_33, dist_coeffs_31);
                _S4463 = !_S4471;
            }
            if(_S4463)
            {
                _S4470 = make_float2 (-10000.0f);
                break;
            }
            float u_40 = uv_33.x;
            float v_40 = uv_33.y;
            float r2_40 = u_40 * u_40 + v_40 * v_40;
            float2  _S4472 = uv_33 * make_float2 (1.0f + r2_40 * ((*dist_coeffs_31)[int(0)] + r2_40 * ((*dist_coeffs_31)[int(1)] + r2_40 * ((*dist_coeffs_31)[int(2)] + r2_40 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_40 * v_40 + (*dist_coeffs_31)[int(5)] * (r2_40 + 2.0f * u_40 * u_40) + (*dist_coeffs_31)[int(6)] * r2_40, 2.0f * (*dist_coeffs_31)[int(5)] * u_40 * v_40 + (*dist_coeffs_31)[int(4)] * (r2_40 + 2.0f * v_40 * v_40) + (*dist_coeffs_31)[int(7)] * r2_40);
            float2  _S4473 = _S4472 + make_float2 ((*dist_coeffs_31)[int(8)] * _S4472.x + (*dist_coeffs_31)[int(9)] * _S4472.y, 0.0f);
            _S4470 = make_float2 (fx_25 * _S4473.x + cx_25, fy_25 * _S4473.y + cy_25);
            break;
        }
        float2  _S4474;
        for(;;)
        {
            float2  uv_34 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4468);
            if(_S4468 < 0.0f)
            {
                _S4463 = true;
            }
            else
            {
                bool _S4475 = is_valid_distortion(uv_34, dist_coeffs_31);
                _S4463 = !_S4475;
            }
            if(_S4463)
            {
                _S4474 = make_float2 (-10000.0f);
                break;
            }
            float u_41 = uv_34.x;
            float v_41 = uv_34.y;
            float r2_41 = u_41 * u_41 + v_41 * v_41;
            float2  _S4476 = uv_34 * make_float2 (1.0f + r2_41 * ((*dist_coeffs_31)[int(0)] + r2_41 * ((*dist_coeffs_31)[int(1)] + r2_41 * ((*dist_coeffs_31)[int(2)] + r2_41 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_41 * v_41 + (*dist_coeffs_31)[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + (*dist_coeffs_31)[int(6)] * r2_41, 2.0f * (*dist_coeffs_31)[int(5)] * u_41 * v_41 + (*dist_coeffs_31)[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + (*dist_coeffs_31)[int(7)] * r2_41);
            float2  _S4477 = _S4476 + make_float2 ((*dist_coeffs_31)[int(8)] * _S4476.x + (*dist_coeffs_31)[int(9)] * _S4476.y, 0.0f);
            _S4474 = make_float2 (fx_25 * _S4477.x + cx_25, fy_25 * _S4477.y + cy_25);
            break;
        }
        float2  _S4478;
        for(;;)
        {
            float2  uv_35 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4469);
            if(_S4469 < 0.0f)
            {
                _S4463 = true;
            }
            else
            {
                bool _S4479 = is_valid_distortion(uv_35, dist_coeffs_31);
                _S4463 = !_S4479;
            }
            if(_S4463)
            {
                _S4478 = make_float2 (-10000.0f);
                break;
            }
            float u_42 = uv_35.x;
            float v_42 = uv_35.y;
            float r2_42 = u_42 * u_42 + v_42 * v_42;
            float2  _S4480 = uv_35 * make_float2 (1.0f + r2_42 * ((*dist_coeffs_31)[int(0)] + r2_42 * ((*dist_coeffs_31)[int(1)] + r2_42 * ((*dist_coeffs_31)[int(2)] + r2_42 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_42 * v_42 + (*dist_coeffs_31)[int(5)] * (r2_42 + 2.0f * u_42 * u_42) + (*dist_coeffs_31)[int(6)] * r2_42, 2.0f * (*dist_coeffs_31)[int(5)] * u_42 * v_42 + (*dist_coeffs_31)[int(4)] * (r2_42 + 2.0f * v_42 * v_42) + (*dist_coeffs_31)[int(7)] * r2_42);
            float2  _S4481 = _S4480 + make_float2 ((*dist_coeffs_31)[int(8)] * _S4480.x + (*dist_coeffs_31)[int(9)] * _S4480.y, 0.0f);
            _S4478 = make_float2 (fx_25 * _S4481.x + cx_25, fy_25 * _S4481.y + cy_25);
            break;
        }
        float2  e0_10 = _S4474 - _S4470;
        float2  e1_10 = _S4478 - _S4474;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(_S4470 - _S4478)));
        float _S4482 = _S4470.x;
        float _S4483 = _S4474.x;
        float _S4484 = _S4478.x;
        float xmax_7 = (F32_max(((F32_max((_S4482), (_S4483)))), (_S4484))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4482), (_S4483)))), (_S4484))) - offset_4;
        float _S4485 = _S4470.y;
        float _S4486 = _S4474.y;
        float _S4487 = _S4478.y;
        float ymax_7 = (F32_max(((F32_max((_S4485), (_S4486)))), (_S4487))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4485), (_S4486)))), (_S4487))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = xmin_7 >= float(image_width_21);
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = ymax_7 <= 0.0f;
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            _S4463 = ymin_7 >= float(image_height_21);
        }
        if(_S4463)
        {
            _S4463 = true;
        }
        else
        {
            if(_S4462 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4463 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4463 = false;
                }
                if(_S4463)
                {
                    _S4463 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4463 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4463 = false;
                    }
                }
            }
            else
            {
                _S4463 = false;
            }
        }
        if(_S4463)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4488 = mean_25 - - mul_0(transpose_0(R_25), t_24);
        float _S4489 = _S4488.x;
        float _S4490 = _S4488.y;
        float _S4491 = _S4488.z;
        float norm_14 = (F32_sqrt((_S4489 * _S4489 + _S4490 * _S4490 + _S4491 * _S4491)));
        float x_53 = _S4489 / norm_14;
        float y_24 = _S4490 / norm_14;
        float z_21 = _S4491 / norm_14;
        float z2_48 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_53 * x_53 - y_24 * y_24;
        float fS1_21 = 2.0f * x_53 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_48 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_53) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_48 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_53) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_53 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_48 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_53) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_53 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4492 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4492);
        float3  _S4493 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4494 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4493 + _S4494 + make_float3 (0.5f), _S4492);
        (*rgbs_0)[int(2)] = max_0(_S4493 - _S4494 + make_float3 (0.5f), _S4492);
        (*verts_0)[int(0)] = vert0_2;
        (*verts_0)[int(1)] = vert1_2;
        (*verts_0)[int(2)] = vert2_2;
        float3  _S4495 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4495 * make_float3 (float(- (F32_sign((dot_0(_S4495, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_26, float fy_26, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_26) + t_25;
        float _S4496 = length_1(mean_c_22);
        bool _S4497;
        if(_S4496 < near_plane_15)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4496 > far_plane_15;
        }
        if(_S4497)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4498 = scale_26.x;
        float sx_5 = (F32_exp((_S4498)));
        float _S4499 = scale_26.y;
        float sy_5 = (F32_exp((_S4499)));
        float sz_7 = scale_26.z - 0.5f * (_S4498 + _S4499);
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
        Matrix<float, 3, 3>  _S4500 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_49), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_49), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4500, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4500, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4500, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_7 = mul_0(R_26, vert0_3) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_3) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_3) + t_25;
        float _S4501 = length_1(vert0_c_7);
        float _S4502 = length_1(vert1_c_7);
        float _S4503 = length_1(vert2_c_7);
        if(_S4501 < near_plane_15)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4501 > far_plane_15;
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4502 < near_plane_15;
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4502 > far_plane_15;
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4503 < near_plane_15;
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = _S4503 > far_plane_15;
        }
        if(_S4497)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S4504 = float2 {vert0_c_7.x, vert0_c_7.y};
        float r_15 = length_0(_S4504);
        float _S4505 = vert0_c_7.z;
        float theta_10 = (F32_atan2((r_15), (_S4505)));
        float k_7;
        if(theta_10 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_10 * theta_10 / 3.0f) / _S4505;
        }
        else
        {
            k_7 = theta_10 / r_15;
        }
        float2  _S4506 = _S4504 * make_float2 (k_7);
        float u_43 = _S4506.x;
        float v_43 = _S4506.y;
        float r2_43 = u_43 * u_43 + v_43 * v_43;
        float _S4507 = 2.0f * (*dist_coeffs_32)[int(4)];
        float _S4508 = 2.0f * (*dist_coeffs_32)[int(5)];
        float2  _S4509 = _S4506 * make_float2 (1.0f + r2_43 * ((*dist_coeffs_32)[int(0)] + r2_43 * ((*dist_coeffs_32)[int(1)] + r2_43 * ((*dist_coeffs_32)[int(2)] + r2_43 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S4507 * u_43 * v_43 + (*dist_coeffs_32)[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + (*dist_coeffs_32)[int(6)] * r2_43, _S4508 * u_43 * v_43 + (*dist_coeffs_32)[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + (*dist_coeffs_32)[int(7)] * r2_43);
        float2  _S4510 = _S4509 + make_float2 ((*dist_coeffs_32)[int(8)] * _S4509.x + (*dist_coeffs_32)[int(9)] * _S4509.y, 0.0f);
        float _S4511 = fx_26 * _S4510.x + cx_26;
        float _S4512 = fy_26 * _S4510.y + cy_26;
        float2  _S4513 = make_float2 (_S4511, _S4512);
        float2  _S4514 = float2 {vert1_c_7.x, vert1_c_7.y};
        float r_16 = length_0(_S4514);
        float _S4515 = vert1_c_7.z;
        float theta_11 = (F32_atan2((r_16), (_S4515)));
        if(theta_11 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_11 * theta_11 / 3.0f) / _S4515;
        }
        else
        {
            k_7 = theta_11 / r_16;
        }
        float2  _S4516 = _S4514 * make_float2 (k_7);
        float u_44 = _S4516.x;
        float v_44 = _S4516.y;
        float r2_44 = u_44 * u_44 + v_44 * v_44;
        float2  _S4517 = _S4516 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_32)[int(0)] + r2_44 * ((*dist_coeffs_32)[int(1)] + r2_44 * ((*dist_coeffs_32)[int(2)] + r2_44 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S4507 * u_44 * v_44 + (*dist_coeffs_32)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_32)[int(6)] * r2_44, _S4508 * u_44 * v_44 + (*dist_coeffs_32)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_32)[int(7)] * r2_44);
        float2  _S4518 = _S4517 + make_float2 ((*dist_coeffs_32)[int(8)] * _S4517.x + (*dist_coeffs_32)[int(9)] * _S4517.y, 0.0f);
        float _S4519 = fx_26 * _S4518.x + cx_26;
        float _S4520 = fy_26 * _S4518.y + cy_26;
        float2  _S4521 = make_float2 (_S4519, _S4520);
        float2  _S4522 = float2 {vert2_c_7.x, vert2_c_7.y};
        float r_17 = length_0(_S4522);
        float _S4523 = vert2_c_7.z;
        float theta_12 = (F32_atan2((r_17), (_S4523)));
        if(theta_12 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_12 * theta_12 / 3.0f) / _S4523;
        }
        else
        {
            k_7 = theta_12 / r_17;
        }
        float2  _S4524 = _S4522 * make_float2 (k_7);
        float u_45 = _S4524.x;
        float v_45 = _S4524.y;
        float r2_45 = u_45 * u_45 + v_45 * v_45;
        float2  _S4525 = _S4524 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_32)[int(0)] + r2_45 * ((*dist_coeffs_32)[int(1)] + r2_45 * ((*dist_coeffs_32)[int(2)] + r2_45 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S4507 * u_45 * v_45 + (*dist_coeffs_32)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_32)[int(6)] * r2_45, _S4508 * u_45 * v_45 + (*dist_coeffs_32)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_32)[int(7)] * r2_45);
        float2  _S4526 = _S4525 + make_float2 ((*dist_coeffs_32)[int(8)] * _S4525.x + (*dist_coeffs_32)[int(9)] * _S4525.y, 0.0f);
        float _S4527 = fx_26 * _S4526.x + cx_26;
        float _S4528 = fy_26 * _S4526.y + cy_26;
        float2  _S4529 = make_float2 (_S4527, _S4528);
        float2  e0_11 = _S4521 - _S4513;
        float2  e1_11 = _S4529 - _S4521;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(_S4513 - _S4529)));
        float xmax_8 = (F32_max(((F32_max((_S4511), (_S4519)))), (_S4527))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4511), (_S4519)))), (_S4527))) - offset_5;
        float ymax_8 = (F32_max(((F32_max((_S4512), (_S4520)))), (_S4528))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4512), (_S4520)))), (_S4528))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = xmin_8 >= float(image_width_22);
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = ymax_8 <= 0.0f;
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            _S4497 = ymin_8 >= float(image_height_22);
        }
        if(_S4497)
        {
            _S4497 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4497 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4497 = false;
                }
                if(_S4497)
                {
                    _S4497 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4497 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4497 = false;
                    }
                }
            }
            else
            {
                _S4497 = false;
            }
        }
        if(_S4497)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4530 = mean_26 - - mul_0(transpose_0(R_26), t_25);
        float _S4531 = _S4530.x;
        float _S4532 = _S4530.y;
        float _S4533 = _S4530.z;
        float norm_15 = (F32_sqrt((_S4531 * _S4531 + _S4532 * _S4532 + _S4533 * _S4533)));
        float x_55 = _S4531 / norm_15;
        float y_25 = _S4532 / norm_15;
        float z_22 = _S4533 / norm_15;
        float z2_50 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_55 * x_55 - y_25 * y_25;
        float fS1_22 = 2.0f * x_55 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_50 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_55) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_50 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_55) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_55 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_50 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_55) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_55 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4534 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4534);
        float3  _S4535 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4536 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4535 + _S4536 + make_float3 (0.5f), _S4534);
        (*rgbs_1)[int(2)] = max_0(_S4535 - _S4536 + make_float3 (0.5f), _S4534);
        (*verts_1)[int(0)] = vert0_3;
        (*verts_1)[int(1)] = vert1_3;
        (*verts_1)[int(2)] = vert2_3;
        float3  _S4537 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4537 * make_float3 (float(- (F32_sign((dot_0(_S4537, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_27, float fy_27, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    for(;;)
    {
        float2  _S4538;
        bool _S4539;
        float3  mean_c_23 = mul_0(R_27, mean_27) + t_26;
        float _S4540 = scale_27.x;
        float sx_6 = (F32_exp((_S4540)));
        float _S4541 = scale_27.y;
        float sy_6 = (F32_exp((_S4541)));
        float sz_8 = scale_27.z - 0.5f * (_S4540 + _S4541);
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
        Matrix<float, 3, 3>  _S4542 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_51), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_51), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
        float3  vert0_4 = mul_0(_S4542, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
        float3  vert1_4 = mul_0(_S4542, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
        float3  vert2_4 = mul_0(_S4542, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
        float3  vert0_c_8 = mul_0(R_27, vert0_4) + t_26;
        float3  vert1_c_8 = mul_0(R_27, vert1_4) + t_26;
        float3  vert2_c_8 = mul_0(R_27, vert2_4) + t_26;
        for(;;)
        {
            float _S4543 = vert0_c_8.z;
            float2  uv_36 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (_S4543);
            if(_S4543 < 0.0f)
            {
                _S4539 = true;
            }
            else
            {
                bool _S4544 = is_valid_distortion(uv_36, dist_coeffs_33);
                _S4539 = !_S4544;
            }
            if(_S4539)
            {
                _S4538 = make_float2 (-10000.0f);
                break;
            }
            float u_46 = uv_36.x;
            float v_46 = uv_36.y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S4545 = uv_36 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_33)[int(0)] + r2_46 * ((*dist_coeffs_33)[int(1)] + r2_46 * ((*dist_coeffs_33)[int(2)] + r2_46 * (*dist_coeffs_33)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_33)[int(4)] * u_46 * v_46 + (*dist_coeffs_33)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_33)[int(6)] * r2_46, 2.0f * (*dist_coeffs_33)[int(5)] * u_46 * v_46 + (*dist_coeffs_33)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_33)[int(7)] * r2_46);
            float2  _S4546 = _S4545 + make_float2 ((*dist_coeffs_33)[int(8)] * _S4545.x + (*dist_coeffs_33)[int(9)] * _S4545.y, 0.0f);
            _S4538 = make_float2 (fx_27 * _S4546.x + cx_27, fy_27 * _S4546.y + cy_27);
            break;
        }
        float2  _S4547;
        for(;;)
        {
            float _S4548 = vert1_c_8.z;
            float2  uv_37 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (_S4548);
            if(_S4548 < 0.0f)
            {
                _S4539 = true;
            }
            else
            {
                bool _S4549 = is_valid_distortion(uv_37, dist_coeffs_33);
                _S4539 = !_S4549;
            }
            if(_S4539)
            {
                _S4547 = make_float2 (-10000.0f);
                break;
            }
            float u_47 = uv_37.x;
            float v_47 = uv_37.y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S4550 = uv_37 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_33)[int(0)] + r2_47 * ((*dist_coeffs_33)[int(1)] + r2_47 * ((*dist_coeffs_33)[int(2)] + r2_47 * (*dist_coeffs_33)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_33)[int(4)] * u_47 * v_47 + (*dist_coeffs_33)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_33)[int(6)] * r2_47, 2.0f * (*dist_coeffs_33)[int(5)] * u_47 * v_47 + (*dist_coeffs_33)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_33)[int(7)] * r2_47);
            float2  _S4551 = _S4550 + make_float2 ((*dist_coeffs_33)[int(8)] * _S4550.x + (*dist_coeffs_33)[int(9)] * _S4550.y, 0.0f);
            _S4547 = make_float2 (fx_27 * _S4551.x + cx_27, fy_27 * _S4551.y + cy_27);
            break;
        }
        float2  _S4552;
        for(;;)
        {
            float _S4553 = vert2_c_8.z;
            float2  uv_38 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (_S4553);
            if(_S4553 < 0.0f)
            {
                _S4539 = true;
            }
            else
            {
                bool _S4554 = is_valid_distortion(uv_38, dist_coeffs_33);
                _S4539 = !_S4554;
            }
            if(_S4539)
            {
                _S4552 = make_float2 (-10000.0f);
                break;
            }
            float u_48 = uv_38.x;
            float v_48 = uv_38.y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S4555 = uv_38 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_33)[int(0)] + r2_48 * ((*dist_coeffs_33)[int(1)] + r2_48 * ((*dist_coeffs_33)[int(2)] + r2_48 * (*dist_coeffs_33)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_33)[int(4)] * u_48 * v_48 + (*dist_coeffs_33)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_33)[int(6)] * r2_48, 2.0f * (*dist_coeffs_33)[int(5)] * u_48 * v_48 + (*dist_coeffs_33)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_33)[int(7)] * r2_48);
            float2  _S4556 = _S4555 + make_float2 ((*dist_coeffs_33)[int(8)] * _S4555.x + (*dist_coeffs_33)[int(9)] * _S4555.y, 0.0f);
            _S4552 = make_float2 (fx_27 * _S4556.x + cx_27, fy_27 * _S4556.y + cy_27);
            break;
        }
        float2  e0_12 = _S4547 - _S4538;
        float2  e1_12 = _S4552 - _S4547;
        float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(_S4538 - _S4552)));
        float _S4557 = _S4538.x;
        float _S4558 = _S4547.x;
        float _S4559 = _S4552.x;
        float _S4560 = _S4538.y;
        float _S4561 = _S4547.y;
        float _S4562 = _S4552.y;
        *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4557), (_S4558)))), (_S4559))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4560), (_S4561)))), (_S4562))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4557), (_S4558)))), (_S4559))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4560), (_S4561)))), (_S4562))) + offset_6)))));
        *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4563 = mean_27 - - mul_0(transpose_0(R_27), t_26);
        float _S4564 = _S4563.x;
        float _S4565 = _S4563.y;
        float _S4566 = _S4563.z;
        float norm_16 = (F32_sqrt((_S4564 * _S4564 + _S4565 * _S4565 + _S4566 * _S4566)));
        float x_57 = _S4564 / norm_16;
        float y_26 = _S4565 / norm_16;
        float z_23 = _S4566 / norm_16;
        float z2_52 = z_23 * z_23;
        float fTmp0B_23 = -1.09254848957061768f * z_23;
        float fC1_23 = x_57 * x_57 - y_26 * y_26;
        float fS1_23 = 2.0f * x_57 * y_26;
        float fTmp0C_23 = -2.28522896766662598f * z2_52 + 0.4570457935333252f;
        float fTmp1B_23 = 1.44530570507049561f * z_23;
        float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_57) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_52 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_57) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_57 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_52 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_57) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_57 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
        float3  _S4567 = make_float3 (0.0f);
        (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4567);
        float3  _S4568 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
        float3  _S4569 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_2)[int(1)] = max_0(_S4568 + _S4569 + make_float3 (0.5f), _S4567);
        (*rgbs_2)[int(2)] = max_0(_S4568 - _S4569 + make_float3 (0.5f), _S4567);
        (*verts_2)[int(0)] = vert0_4;
        (*verts_2)[int(1)] = vert1_4;
        (*verts_2)[int(2)] = vert2_4;
        float3  _S4570 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
        *normal_6 = _S4570 * make_float3 (float(- (F32_sign((dot_0(_S4570, mean_c_23))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_28, float fy_28, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_34, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_28) + t_27;
    float _S4571 = scale_28.x;
    float sx_7 = (F32_exp((_S4571)));
    float _S4572 = scale_28.y;
    float sy_7 = (F32_exp((_S4572)));
    float sz_9 = scale_28.z - 0.5f * (_S4571 + _S4572);
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
    Matrix<float, 3, 3>  _S4573 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_53), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_53), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4573, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4573, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4573, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_9 = mul_0(R_28, vert0_5) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_5) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_5) + t_27;
    float2  _S4574 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_18 = length_0(_S4574);
    float _S4575 = vert0_c_9.z;
    float theta_13 = (F32_atan2((r_18), (_S4575)));
    float k_8;
    if(theta_13 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_13 * theta_13 / 3.0f) / _S4575;
    }
    else
    {
        k_8 = theta_13 / r_18;
    }
    float2  _S4576 = _S4574 * make_float2 (k_8);
    float u_49 = _S4576.x;
    float v_49 = _S4576.y;
    float r2_49 = u_49 * u_49 + v_49 * v_49;
    float _S4577 = 2.0f * (*dist_coeffs_34)[int(4)];
    float _S4578 = 2.0f * (*dist_coeffs_34)[int(5)];
    float2  _S4579 = _S4576 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_34)[int(0)] + r2_49 * ((*dist_coeffs_34)[int(1)] + r2_49 * ((*dist_coeffs_34)[int(2)] + r2_49 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S4577 * u_49 * v_49 + (*dist_coeffs_34)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_34)[int(6)] * r2_49, _S4578 * u_49 * v_49 + (*dist_coeffs_34)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_34)[int(7)] * r2_49);
    float2  _S4580 = _S4579 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4579.x + (*dist_coeffs_34)[int(9)] * _S4579.y, 0.0f);
    float _S4581 = fx_28 * _S4580.x + cx_28;
    float _S4582 = fy_28 * _S4580.y + cy_28;
    float2  _S4583 = make_float2 (_S4581, _S4582);
    float2  _S4584 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_19 = length_0(_S4584);
    float _S4585 = vert1_c_9.z;
    float theta_14 = (F32_atan2((r_19), (_S4585)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_14 * theta_14 / 3.0f) / _S4585;
    }
    else
    {
        k_8 = theta_14 / r_19;
    }
    float2  _S4586 = _S4584 * make_float2 (k_8);
    float u_50 = _S4586.x;
    float v_50 = _S4586.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float2  _S4587 = _S4586 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_34)[int(0)] + r2_50 * ((*dist_coeffs_34)[int(1)] + r2_50 * ((*dist_coeffs_34)[int(2)] + r2_50 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S4577 * u_50 * v_50 + (*dist_coeffs_34)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_34)[int(6)] * r2_50, _S4578 * u_50 * v_50 + (*dist_coeffs_34)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_34)[int(7)] * r2_50);
    float2  _S4588 = _S4587 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4587.x + (*dist_coeffs_34)[int(9)] * _S4587.y, 0.0f);
    float _S4589 = fx_28 * _S4588.x + cx_28;
    float _S4590 = fy_28 * _S4588.y + cy_28;
    float2  _S4591 = make_float2 (_S4589, _S4590);
    float2  _S4592 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_20 = length_0(_S4592);
    float _S4593 = vert2_c_9.z;
    float theta_15 = (F32_atan2((r_20), (_S4593)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_15 * theta_15 / 3.0f) / _S4593;
    }
    else
    {
        k_8 = theta_15 / r_20;
    }
    float2  _S4594 = _S4592 * make_float2 (k_8);
    float u_51 = _S4594.x;
    float v_51 = _S4594.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S4595 = _S4594 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_34)[int(0)] + r2_51 * ((*dist_coeffs_34)[int(1)] + r2_51 * ((*dist_coeffs_34)[int(2)] + r2_51 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S4577 * u_51 * v_51 + (*dist_coeffs_34)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_34)[int(6)] * r2_51, _S4578 * u_51 * v_51 + (*dist_coeffs_34)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_34)[int(7)] * r2_51);
    float2  _S4596 = _S4595 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4595.x + (*dist_coeffs_34)[int(9)] * _S4595.y, 0.0f);
    float _S4597 = fx_28 * _S4596.x + cx_28;
    float _S4598 = fy_28 * _S4596.y + cy_28;
    float2  _S4599 = make_float2 (_S4597, _S4598);
    float2  e0_13 = _S4591 - _S4583;
    float2  e1_13 = _S4599 - _S4591;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(_S4583 - _S4599)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4581), (_S4589)))), (_S4597))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4582), (_S4590)))), (_S4598))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4581), (_S4589)))), (_S4597))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4582), (_S4590)))), (_S4598))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4600 = mean_28 - - mul_0(transpose_0(R_28), t_27);
    float _S4601 = _S4600.x;
    float _S4602 = _S4600.y;
    float _S4603 = _S4600.z;
    float norm_17 = (F32_sqrt((_S4601 * _S4601 + _S4602 * _S4602 + _S4603 * _S4603)));
    float x_59 = _S4601 / norm_17;
    float y_27 = _S4602 / norm_17;
    float z_24 = _S4603 / norm_17;
    float z2_54 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_59 * x_59 - y_27 * y_27;
    float fS1_24 = 2.0f * x_59 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_54 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_59) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_54 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_59) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_59 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_54 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_59) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_59 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4604 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4604);
    float3  _S4605 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4606 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4605 + _S4606 + make_float3 (0.5f), _S4604);
    (*rgbs_3)[int(2)] = max_0(_S4605 - _S4606 + make_float3 (0.5f), _S4604);
    (*verts_3)[int(0)] = vert0_5;
    (*verts_3)[int(1)] = vert1_5;
    (*verts_3)[int(2)] = vert2_5;
    float3  _S4607 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4607 * make_float3 (float(- (F32_sign((dot_0(_S4607, mean_c_24))))));
    return;
}

struct s_bwd_prop_projection_opaque_triangle_eval3d_persp_differentiable_Intermediates_0
{
    bool _S4608;
    bool _S4609;
    bool _S4610;
};

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_29, float fy_29, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    s_bwd_prop_projection_opaque_triangle_eval3d_persp_differentiable_Intermediates_0 _S4611;
    (&_S4611)->_S4608 = false;
    (&_S4611)->_S4609 = false;
    (&_S4611)->_S4610 = false;
    (&_S4611)->_S4608 = false;
    (&_S4611)->_S4609 = false;
    (&_S4611)->_S4610 = false;
    float _S4612 = scale_29.x;
    float _S4613 = s_primal_ctx_exp_1(_S4612);
    float _S4614 = scale_29.y;
    float _S4615 = s_primal_ctx_exp_1(_S4614);
    float sz_10 = scale_29.z - 0.5f * (_S4612 + _S4614);
    float _S4616 = quat_30.y;
    float x2_30 = _S4616 * _S4616;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_55 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S4617 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_55), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_55), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  _S4618 = make_float3 (_S4613, 0.0f, 0.0f);
    float3  vert0_6 = s_primal_ctx_mul_1(_S4617, _S4618) + mean_29;
    float _S4619 = -0.5f + sz_10;
    float3  _S4620 = make_float3 (_S4613 * _S4619, _S4615, 0.0f);
    float3  vert1_6 = s_primal_ctx_mul_1(_S4617, _S4620) + mean_29;
    float _S4621 = -0.5f - sz_10;
    float3  _S4622 = make_float3 (_S4613 * _S4621, - _S4615, 0.0f);
    float3  vert2_6 = s_primal_ctx_mul_1(_S4617, _S4622) + mean_29;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_6) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_6) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_6) + t_28;
    float2  _S4623 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4624 = vert0_c_10.z;
    float2  uv_39 = _S4623 / make_float2 (_S4624);
    bool _S4625 = _S4624 < 0.0f;
    if(_S4625)
    {
    }
    else
    {
        bool _S4626 = is_valid_distortion(uv_39, dist_coeffs_35);
        (&_S4611)->_S4608 = _S4626;
    }
    float2  _S4627 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4628 = vert1_c_10.z;
    float2  uv_40 = _S4627 / make_float2 (_S4628);
    bool _S4629 = _S4628 < 0.0f;
    if(_S4629)
    {
    }
    else
    {
        bool _S4630 = is_valid_distortion(uv_40, dist_coeffs_35);
        (&_S4611)->_S4609 = _S4630;
    }
    float2  _S4631 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4632 = vert2_c_10.z;
    float2  uv_41 = _S4631 / make_float2 (_S4632);
    bool _S4633 = _S4632 < 0.0f;
    if(_S4633)
    {
    }
    else
    {
        bool _S4634 = is_valid_distortion(uv_41, dist_coeffs_35);
        (&_S4611)->_S4610 = _S4634;
    }
    s_bwd_prop_projection_opaque_triangle_eval3d_persp_differentiable_Intermediates_0 _S4635 = _S4611;
    float2  _S4636 = make_float2 (0.0f);
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_29) + t_28;
    float2  _S4637 = make_float2 (_S4624);
    float2  uv_42 = _S4623 / make_float2 (_S4624);
    float2  _S4638 = make_float2 (_S4624 * _S4624);
    bool _S4639;
    if(_S4625)
    {
        _S4639 = true;
    }
    else
    {
        _S4639 = !_S4635._S4608;
    }
    float2  _S4640;
    if(_S4639)
    {
        _S4640 = make_float2 (-10000.0f);
    }
    bool _S4641 = !_S4639;
    float2  _S4642;
    float _S4643;
    float _S4644;
    float _S4645;
    float _S4646;
    float _S4647;
    float _S4648;
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
    float _S4660;
    float _S4661;
    if(_S4641)
    {
        float u_52 = uv_42.x;
        float v_52 = uv_42.y;
        float r2_52 = u_52 * u_52 + v_52 * v_52;
        float _S4662 = (*dist_coeffs_35)[int(2)] + r2_52 * (*dist_coeffs_35)[int(3)];
        float _S4663 = (*dist_coeffs_35)[int(1)] + r2_52 * _S4662;
        float _S4664 = (*dist_coeffs_35)[int(0)] + r2_52 * _S4663;
        float radial_9 = 1.0f + r2_52 * _S4664;
        float2  _S4665 = make_float2 (radial_9);
        float _S4666 = 2.0f * (*dist_coeffs_35)[int(4)];
        float _S4667 = _S4666 * u_52;
        float _S4668 = 2.0f * u_52;
        float _S4669 = 2.0f * (*dist_coeffs_35)[int(5)];
        float _S4670 = _S4669 * u_52;
        float _S4671 = 2.0f * v_52;
        float2  _S4672 = uv_42 * make_float2 (radial_9) + make_float2 (_S4667 * v_52 + (*dist_coeffs_35)[int(5)] * (r2_52 + _S4668 * u_52) + (*dist_coeffs_35)[int(6)] * r2_52, _S4670 * v_52 + (*dist_coeffs_35)[int(4)] * (r2_52 + _S4671 * v_52) + (*dist_coeffs_35)[int(7)] * r2_52);
        float2  _S4673 = _S4672 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4672.x + (*dist_coeffs_35)[int(9)] * _S4672.y, 0.0f);
        _S4640 = make_float2 (fx_29 * _S4673.x + cx_29, fy_29 * _S4673.y + cy_29);
        _S4643 = (*dist_coeffs_35)[int(9)];
        _S4644 = (*dist_coeffs_35)[int(8)];
        _S4642 = _S4665;
        _S4645 = (*dist_coeffs_35)[int(7)];
        _S4646 = (*dist_coeffs_35)[int(4)];
        _S4647 = _S4671;
        _S4648 = v_52;
        _S4649 = _S4670;
        _S4650 = _S4669;
        _S4651 = (*dist_coeffs_35)[int(6)];
        _S4652 = (*dist_coeffs_35)[int(5)];
        _S4653 = _S4668;
        _S4654 = u_52;
        _S4655 = _S4667;
        _S4656 = _S4666;
        _S4657 = r2_52;
        _S4658 = _S4664;
        _S4659 = _S4663;
        _S4660 = _S4662;
        _S4661 = (*dist_coeffs_35)[int(3)];
    }
    else
    {
        _S4643 = 0.0f;
        _S4644 = 0.0f;
        _S4642 = _S4636;
        _S4645 = 0.0f;
        _S4646 = 0.0f;
        _S4647 = 0.0f;
        _S4648 = 0.0f;
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
        _S4660 = 0.0f;
        _S4661 = 0.0f;
    }
    float2  _S4674 = make_float2 (_S4628);
    float2  uv_43 = _S4627 / make_float2 (_S4628);
    float2  _S4675 = make_float2 (_S4628 * _S4628);
    if(_S4629)
    {
        _S4639 = true;
    }
    else
    {
        _S4639 = !_S4635._S4609;
    }
    float2  _S4676;
    if(_S4639)
    {
        _S4676 = make_float2 (-10000.0f);
    }
    bool _S4677 = !_S4639;
    float2  _S4678;
    float _S4679;
    float _S4680;
    float _S4681;
    float _S4682;
    float _S4683;
    float _S4684;
    float _S4685;
    float _S4686;
    float _S4687;
    float _S4688;
    float _S4689;
    float _S4690;
    float _S4691;
    float _S4692;
    float _S4693;
    float _S4694;
    float _S4695;
    float _S4696;
    float _S4697;
    if(_S4677)
    {
        float u_53 = uv_43.x;
        float v_53 = uv_43.y;
        float r2_53 = u_53 * u_53 + v_53 * v_53;
        float _S4698 = (*dist_coeffs_35)[int(2)] + r2_53 * (*dist_coeffs_35)[int(3)];
        float _S4699 = (*dist_coeffs_35)[int(1)] + r2_53 * _S4698;
        float _S4700 = (*dist_coeffs_35)[int(0)] + r2_53 * _S4699;
        float radial_10 = 1.0f + r2_53 * _S4700;
        float2  _S4701 = make_float2 (radial_10);
        float _S4702 = 2.0f * (*dist_coeffs_35)[int(4)];
        float _S4703 = _S4702 * u_53;
        float _S4704 = 2.0f * u_53;
        float _S4705 = 2.0f * (*dist_coeffs_35)[int(5)];
        float _S4706 = _S4705 * u_53;
        float _S4707 = 2.0f * v_53;
        float2  _S4708 = uv_43 * make_float2 (radial_10) + make_float2 (_S4703 * v_53 + (*dist_coeffs_35)[int(5)] * (r2_53 + _S4704 * u_53) + (*dist_coeffs_35)[int(6)] * r2_53, _S4706 * v_53 + (*dist_coeffs_35)[int(4)] * (r2_53 + _S4707 * v_53) + (*dist_coeffs_35)[int(7)] * r2_53);
        float2  _S4709 = _S4708 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4708.x + (*dist_coeffs_35)[int(9)] * _S4708.y, 0.0f);
        _S4676 = make_float2 (fx_29 * _S4709.x + cx_29, fy_29 * _S4709.y + cy_29);
        _S4679 = (*dist_coeffs_35)[int(9)];
        _S4680 = (*dist_coeffs_35)[int(8)];
        _S4678 = _S4701;
        _S4681 = (*dist_coeffs_35)[int(7)];
        _S4682 = (*dist_coeffs_35)[int(4)];
        _S4683 = _S4707;
        _S4684 = v_53;
        _S4685 = _S4706;
        _S4686 = _S4705;
        _S4687 = (*dist_coeffs_35)[int(6)];
        _S4688 = (*dist_coeffs_35)[int(5)];
        _S4689 = _S4704;
        _S4690 = u_53;
        _S4691 = _S4703;
        _S4692 = _S4702;
        _S4693 = r2_53;
        _S4694 = _S4700;
        _S4695 = _S4699;
        _S4696 = _S4698;
        _S4697 = (*dist_coeffs_35)[int(3)];
    }
    else
    {
        _S4679 = 0.0f;
        _S4680 = 0.0f;
        _S4678 = _S4636;
        _S4681 = 0.0f;
        _S4682 = 0.0f;
        _S4683 = 0.0f;
        _S4684 = 0.0f;
        _S4685 = 0.0f;
        _S4686 = 0.0f;
        _S4687 = 0.0f;
        _S4688 = 0.0f;
        _S4689 = 0.0f;
        _S4690 = 0.0f;
        _S4691 = 0.0f;
        _S4692 = 0.0f;
        _S4693 = 0.0f;
        _S4694 = 0.0f;
        _S4695 = 0.0f;
        _S4696 = 0.0f;
        _S4697 = 0.0f;
    }
    float2  _S4710 = make_float2 (_S4632);
    float2  uv_44 = _S4631 / make_float2 (_S4632);
    float2  _S4711 = make_float2 (_S4632 * _S4632);
    if(_S4633)
    {
        _S4639 = true;
    }
    else
    {
        _S4639 = !_S4635._S4610;
    }
    float2  _S4712;
    if(_S4639)
    {
        _S4712 = make_float2 (-10000.0f);
    }
    bool _S4713 = !_S4639;
    float2  _S4714;
    float _S4715;
    float _S4716;
    float _S4717;
    float _S4718;
    float _S4719;
    float _S4720;
    float _S4721;
    float _S4722;
    float _S4723;
    float _S4724;
    float _S4725;
    float _S4726;
    float _S4727;
    float _S4728;
    float _S4729;
    float _S4730;
    float _S4731;
    float _S4732;
    float _S4733;
    if(_S4713)
    {
        float u_54 = uv_44.x;
        float v_54 = uv_44.y;
        float r2_54 = u_54 * u_54 + v_54 * v_54;
        float _S4734 = (*dist_coeffs_35)[int(2)] + r2_54 * (*dist_coeffs_35)[int(3)];
        float _S4735 = (*dist_coeffs_35)[int(1)] + r2_54 * _S4734;
        float _S4736 = (*dist_coeffs_35)[int(0)] + r2_54 * _S4735;
        float radial_11 = 1.0f + r2_54 * _S4736;
        float2  _S4737 = make_float2 (radial_11);
        float _S4738 = 2.0f * (*dist_coeffs_35)[int(4)];
        float _S4739 = _S4738 * u_54;
        float _S4740 = 2.0f * u_54;
        float _S4741 = 2.0f * (*dist_coeffs_35)[int(5)];
        float _S4742 = _S4741 * u_54;
        float _S4743 = 2.0f * v_54;
        float2  _S4744 = uv_44 * make_float2 (radial_11) + make_float2 (_S4739 * v_54 + (*dist_coeffs_35)[int(5)] * (r2_54 + _S4740 * u_54) + (*dist_coeffs_35)[int(6)] * r2_54, _S4742 * v_54 + (*dist_coeffs_35)[int(4)] * (r2_54 + _S4743 * v_54) + (*dist_coeffs_35)[int(7)] * r2_54);
        float2  _S4745 = _S4744 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4744.x + (*dist_coeffs_35)[int(9)] * _S4744.y, 0.0f);
        _S4712 = make_float2 (fx_29 * _S4745.x + cx_29, fy_29 * _S4745.y + cy_29);
        _S4715 = (*dist_coeffs_35)[int(9)];
        _S4716 = (*dist_coeffs_35)[int(8)];
        _S4714 = _S4737;
        _S4717 = (*dist_coeffs_35)[int(7)];
        _S4718 = (*dist_coeffs_35)[int(4)];
        _S4719 = _S4743;
        _S4720 = v_54;
        _S4721 = _S4742;
        _S4722 = _S4741;
        _S4723 = (*dist_coeffs_35)[int(6)];
        _S4724 = (*dist_coeffs_35)[int(5)];
        _S4725 = _S4740;
        _S4726 = u_54;
        _S4727 = _S4739;
        _S4728 = _S4738;
        _S4729 = r2_54;
        _S4730 = _S4736;
        _S4731 = _S4735;
        _S4732 = _S4734;
        _S4733 = (*dist_coeffs_35)[int(3)];
    }
    else
    {
        _S4715 = 0.0f;
        _S4716 = 0.0f;
        _S4714 = _S4636;
        _S4717 = 0.0f;
        _S4718 = 0.0f;
        _S4719 = 0.0f;
        _S4720 = 0.0f;
        _S4721 = 0.0f;
        _S4722 = 0.0f;
        _S4723 = 0.0f;
        _S4724 = 0.0f;
        _S4725 = 0.0f;
        _S4726 = 0.0f;
        _S4727 = 0.0f;
        _S4728 = 0.0f;
        _S4729 = 0.0f;
        _S4730 = 0.0f;
        _S4731 = 0.0f;
        _S4732 = 0.0f;
        _S4733 = 0.0f;
    }
    float2  e0_14 = _S4676 - _S4640;
    float2  e1_14 = _S4712 - _S4676;
    float2  e2_6 = _S4640 - _S4712;
    float _S4746 = e0_14.x;
    float _S4747 = e1_14.y;
    float _S4748 = e0_14.y;
    float _S4749 = e1_14.x;
    float _S4750 = _S4746 * _S4747 - _S4748 * _S4749;
    float _S4751 = 1.0f - hardness_14.y;
    float _S4752 = -1.0f / _S4751;
    float _S4753 = _S4751 * _S4751;
    float _S4754 = _S4640.x;
    float _S4755 = _S4676.x;
    float _S4756 = s_primal_ctx_max_0(_S4754, _S4755);
    float _S4757 = _S4712.x;
    float _S4758 = s_primal_ctx_min_0(_S4754, _S4755);
    float _S4759 = _S4640.y;
    float _S4760 = _S4676.y;
    float _S4761 = s_primal_ctx_max_0(_S4759, _S4760);
    float _S4762 = _S4712.y;
    float _S4763 = s_primal_ctx_min_0(_S4759, _S4760);
    float3  _S4764 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4765 = length_1(_S4764) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4766 = transpose_0(R_29);
    float3  _S4767 = mean_29 - - s_primal_ctx_mul_1(_S4766, t_28);
    float _S4768 = _S4767.x;
    float _S4769 = _S4767.y;
    float _S4770 = _S4767.z;
    float _S4771 = _S4768 * _S4768 + _S4769 * _S4769 + _S4770 * _S4770;
    float _S4772 = s_primal_ctx_sqrt_0(_S4771);
    float x_60 = _S4768 / _S4772;
    float3  _S4773 = make_float3 (x_60);
    float _S4774 = _S4772 * _S4772;
    float y_28 = _S4769 / _S4772;
    float z_25 = _S4770 / _S4772;
    float3  _S4775 = make_float3 (z_25);
    float _S4776 = - y_28;
    float3  _S4777 = make_float3 (_S4776);
    float z2_56 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_60 * x_60 - y_28 * y_28;
    float _S4778 = 2.0f * x_60;
    float fS1_25 = _S4778 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_56 - 0.31539157032966614f;
    float3  _S4779 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_60;
    float3  _S4780 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4781 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4782 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4783 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_56 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4784 = 1.86588168144226074f * z2_56 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4784;
    float3  _S4785 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_60;
    float3  _S4786 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4787 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4788 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4789 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_60 * fC1_25 - y_28 * fS1_25);
    float3  _S4790 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_60 * fS1_25 + y_28 * fC1_25);
    float3  _S4791 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4776) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_60) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4792 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4793 = make_float3 (0.0f);
    float3  _S4794 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4795 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4796 = make_float3 (_S4795);
    float3  _S4797 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4795);
    float3  _S4798 = _S4794 + _S4797 + make_float3 (0.5f);
    float3  _S4799 = _S4794 - _S4797 + make_float3 (0.5f);
    float3  _S4800 = vert1_c_10 - vert0_c_10;
    float3  _S4801 = vert2_c_10 - vert0_c_10;
    float3  _S4802 = s_primal_ctx_cross_0(_S4800, _S4801);
    float3  _S4803 = normalize_0(_S4802);
    float3  _S4804 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4803, mean_c_25)))))) * v_normal_2;
    float3  _S4805 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4806;
    (&_S4806)->primal_0 = _S4803;
    (&_S4806)->differential_0 = _S4805;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4807;
    (&_S4807)->primal_0 = mean_c_25;
    (&_S4807)->differential_0 = _S4805;
    s_bwd_prop_dot_0(&_S4806, &_S4807, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4808 = _S4807;
    float3  _S4809 = _S4804 + _S4806.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4810;
    (&_S4810)->primal_0 = _S4802;
    (&_S4810)->differential_0 = _S4805;
    s_bwd_normalize_impl_0(&_S4810, _S4809);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4811;
    (&_S4811)->primal_0 = _S4800;
    (&_S4811)->differential_0 = _S4805;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4812;
    (&_S4812)->primal_0 = _S4801;
    (&_S4812)->differential_0 = _S4805;
    s_bwd_prop_cross_0(&_S4811, &_S4812, _S4810.differential_0);
    float3  _S4813 = - _S4812.differential_0;
    float3  _S4814 = - _S4811.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4815;
    (&_S4815)->primal_0 = _S4799;
    (&_S4815)->differential_0 = _S4805;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4816;
    (&_S4816)->primal_0 = _S4793;
    (&_S4816)->differential_0 = _S4805;
    s_bwd_prop_max_0(&_S4815, &_S4816, (*v_rgbs_0)[int(2)]);
    float3  _S4817 = - _S4815.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4818;
    (&_S4818)->primal_0 = _S4798;
    (&_S4818)->differential_0 = _S4805;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4819;
    (&_S4819)->primal_0 = _S4793;
    (&_S4819)->differential_0 = _S4805;
    s_bwd_prop_max_0(&_S4818, &_S4819, (*v_rgbs_0)[int(1)]);
    float3  _S4820 = _S4796 * (_S4817 + _S4818.differential_0);
    float3  _S4821 = _S4815.differential_0 + _S4818.differential_0;
    float3  _S4822 = make_float3 (0.5f) * - _S4821;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4823;
    (&_S4823)->primal_0 = _S4792;
    (&_S4823)->differential_0 = _S4805;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4824;
    (&_S4824)->primal_0 = _S4793;
    (&_S4824)->differential_0 = _S4805;
    s_bwd_prop_max_0(&_S4823, &_S4824, (*v_rgbs_0)[int(0)]);
    float3  _S4825 = _S4822 + _S4823.differential_0;
    float3  _S4826 = _S4821 + _S4823.differential_0;
    float3  _S4827 = _S4790 * _S4826;
    float3  _S4828 = (*sh_coeffs_25)[int(15)] * _S4826;
    float3  _S4829 = _S4788 * _S4826;
    float3  _S4830 = (*sh_coeffs_25)[int(14)] * _S4826;
    float3  _S4831 = _S4786 * _S4826;
    float3  _S4832 = (*sh_coeffs_25)[int(13)] * _S4826;
    float3  _S4833 = _S4785 * _S4826;
    float3  _S4834 = (*sh_coeffs_25)[int(12)] * _S4826;
    float3  _S4835 = _S4787 * _S4826;
    float3  _S4836 = (*sh_coeffs_25)[int(11)] * _S4826;
    float3  _S4837 = _S4789 * _S4826;
    float3  _S4838 = (*sh_coeffs_25)[int(10)] * _S4826;
    float3  _S4839 = _S4791 * _S4826;
    float3  _S4840 = (*sh_coeffs_25)[int(9)] * _S4826;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4840.x + _S4840.y + _S4840.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4828.x + _S4828.y + _S4828.z);
    float _S4841 = _S4838.x + _S4838.y + _S4838.z;
    float _S4842 = _S4830.x + _S4830.y + _S4830.z;
    float _S4843 = _S4836.x + _S4836.y + _S4836.z;
    float _S4844 = _S4832.x + _S4832.y + _S4832.z;
    float _S4845 = _S4834.x + _S4834.y + _S4834.z;
    float _S4846 = - s_diff_fC2_T_7;
    float3  _S4847 = _S4782 * _S4826;
    float3  _S4848 = (*sh_coeffs_25)[int(8)] * _S4826;
    float3  _S4849 = _S4780 * _S4826;
    float3  _S4850 = (*sh_coeffs_25)[int(7)] * _S4826;
    float3  _S4851 = _S4779 * _S4826;
    float3  _S4852 = (*sh_coeffs_25)[int(6)] * _S4826;
    float3  _S4853 = _S4781 * _S4826;
    float3  _S4854 = (*sh_coeffs_25)[int(5)] * _S4826;
    float3  _S4855 = _S4783 * _S4826;
    float3  _S4856 = (*sh_coeffs_25)[int(4)] * _S4826;
    float _S4857 = _S4854.x + _S4854.y + _S4854.z;
    float _S4858 = _S4850.x + _S4850.y + _S4850.z;
    float _S4859 = fTmp1B_25 * _S4841 + x_60 * s_diff_fS2_T_7 + y_28 * _S4846 + 0.54627424478530884f * (_S4856.x + _S4856.y + _S4856.z);
    float _S4860 = fTmp1B_25 * _S4842 + y_28 * s_diff_fS2_T_7 + x_60 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4848.x + _S4848.y + _S4848.z);
    float _S4861 = y_28 * - _S4860;
    float _S4862 = x_60 * _S4860;
    float _S4863 = z_25 * (1.86588168144226074f * (z_25 * _S4845) + -2.28522896766662598f * (y_28 * _S4843 + x_60 * _S4844) + 0.94617468118667603f * (_S4852.x + _S4852.y + _S4852.z));
    float3  _S4864 = make_float3 (0.48860251903533936f) * _S4826;
    float3  _S4865 = - _S4864;
    float3  _S4866 = _S4773 * _S4865;
    float3  _S4867 = (*sh_coeffs_25)[int(3)] * _S4865;
    float3  _S4868 = _S4775 * _S4864;
    float3  _S4869 = (*sh_coeffs_25)[int(2)] * _S4864;
    float3  _S4870 = _S4777 * _S4864;
    float3  _S4871 = (*sh_coeffs_25)[int(1)] * _S4864;
    float _S4872 = (_S4784 * _S4845 + 1.44530570507049561f * (fS1_25 * _S4841 + fC1_25 * _S4842) + -1.09254848957061768f * (y_28 * _S4857 + x_60 * _S4858) + _S4863 + _S4863 + _S4869.x + _S4869.y + _S4869.z) / _S4774;
    float _S4873 = _S4772 * _S4872;
    float _S4874 = (fTmp0C_25 * _S4843 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4846 + fTmp0B_25 * _S4857 + _S4778 * _S4859 + _S4861 + _S4861 + - (_S4871.x + _S4871.y + _S4871.z)) / _S4774;
    float _S4875 = _S4772 * _S4874;
    float _S4876 = (fTmp0C_25 * _S4844 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4858 + 2.0f * (y_28 * _S4859) + _S4862 + _S4862 + _S4867.x + _S4867.y + _S4867.z) / _S4774;
    float _S4877 = _S4772 * _S4876;
    float _S4878 = _S4770 * - _S4872 + _S4769 * - _S4874 + _S4768 * - _S4876;
    DiffPair_float_0 _S4879;
    (&_S4879)->primal_0 = _S4771;
    (&_S4879)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4879, _S4878);
    float _S4880 = _S4770 * _S4879.differential_0;
    float _S4881 = _S4769 * _S4879.differential_0;
    float _S4882 = _S4768 * _S4879.differential_0;
    float3  _S4883 = make_float3 (0.282094806432724f) * _S4826;
    float3  _S4884 = make_float3 (_S4877 + _S4882 + _S4882, _S4875 + _S4881 + _S4881, _S4873 + _S4880 + _S4880);
    float3  _S4885 = - - _S4884;
    Matrix<float, 3, 3>  _S4886 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4887;
    (&_S4887)->primal_0 = _S4766;
    (&_S4887)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4888;
    (&_S4888)->primal_0 = t_28;
    (&_S4888)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S4887, &_S4888, _S4885);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4889 = _S4888;
    Matrix<float, 3, 3>  _S4890 = transpose_0(_S4887.differential_0);
    DiffPair_float_0 _S4891;
    (&_S4891)->primal_0 = _S4765;
    (&_S4891)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4891, v_depth_9);
    float _S4892 = 0.3333333432674408f * _S4891.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4893;
    (&_S4893)->primal_0 = _S4764;
    (&_S4893)->differential_0 = _S4805;
    s_bwd_length_impl_0(&_S4893, _S4892);
    DiffPair_float_0 _S4894;
    (&_S4894)->primal_0 = _S4763;
    (&_S4894)->differential_0 = 0.0f;
    DiffPair_float_0 _S4895;
    (&_S4895)->primal_0 = _S4762;
    (&_S4895)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4894, &_S4895, 0.0f);
    DiffPair_float_0 _S4896;
    (&_S4896)->primal_0 = _S4759;
    (&_S4896)->differential_0 = 0.0f;
    DiffPair_float_0 _S4897;
    (&_S4897)->primal_0 = _S4760;
    (&_S4897)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4896, &_S4897, _S4894.differential_0);
    DiffPair_float_0 _S4898;
    (&_S4898)->primal_0 = _S4761;
    (&_S4898)->differential_0 = 0.0f;
    DiffPair_float_0 _S4899;
    (&_S4899)->primal_0 = _S4762;
    (&_S4899)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4898, &_S4899, 0.0f);
    float _S4900 = _S4895.differential_0 + _S4899.differential_0;
    DiffPair_float_0 _S4901;
    (&_S4901)->primal_0 = _S4759;
    (&_S4901)->differential_0 = 0.0f;
    DiffPair_float_0 _S4902;
    (&_S4902)->primal_0 = _S4760;
    (&_S4902)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4901, &_S4902, _S4898.differential_0);
    float _S4903 = _S4897.differential_0 + _S4902.differential_0;
    float _S4904 = _S4896.differential_0 + _S4901.differential_0;
    DiffPair_float_0 _S4905;
    (&_S4905)->primal_0 = _S4758;
    (&_S4905)->differential_0 = 0.0f;
    DiffPair_float_0 _S4906;
    (&_S4906)->primal_0 = _S4757;
    (&_S4906)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4905, &_S4906, 0.0f);
    DiffPair_float_0 _S4907;
    (&_S4907)->primal_0 = _S4754;
    (&_S4907)->differential_0 = 0.0f;
    DiffPair_float_0 _S4908;
    (&_S4908)->primal_0 = _S4755;
    (&_S4908)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4907, &_S4908, _S4905.differential_0);
    DiffPair_float_0 _S4909;
    (&_S4909)->primal_0 = _S4756;
    (&_S4909)->differential_0 = 0.0f;
    DiffPair_float_0 _S4910;
    (&_S4910)->primal_0 = _S4757;
    (&_S4910)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4909, &_S4910, 0.0f);
    float _S4911 = _S4906.differential_0 + _S4910.differential_0;
    DiffPair_float_0 _S4912;
    (&_S4912)->primal_0 = _S4754;
    (&_S4912)->differential_0 = 0.0f;
    DiffPair_float_0 _S4913;
    (&_S4913)->primal_0 = _S4755;
    (&_S4913)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4912, &_S4913, _S4909.differential_0);
    float _S4914 = _S4908.differential_0 + _S4913.differential_0;
    float _S4915 = _S4907.differential_0 + _S4912.differential_0;
    DiffPair_float_0 _S4916;
    (&_S4916)->primal_0 = _S4752;
    (&_S4916)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4916, 0.0f);
    float _S4917 = - (-1.0f * - (_S4916.differential_0 / _S4753));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4918;
    (&_S4918)->primal_0 = e2_6;
    (&_S4918)->differential_0 = _S4636;
    s_bwd_length_impl_1(&_S4918, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4919;
    (&_S4919)->primal_0 = e1_14;
    (&_S4919)->differential_0 = _S4636;
    s_bwd_length_impl_1(&_S4919, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4920;
    (&_S4920)->primal_0 = e0_14;
    (&_S4920)->differential_0 = _S4636;
    s_bwd_length_impl_1(&_S4920, -0.0f);
    DiffPair_float_0 _S4921;
    (&_S4921)->primal_0 = _S4750;
    (&_S4921)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4921, 0.0f);
    float _S4922 = - _S4921.differential_0;
    float2  _S4923 = - _S4918.differential_0;
    float2  _S4924 = _S4919.differential_0 + make_float2 (_S4748 * _S4922, _S4746 * _S4921.differential_0);
    float2  _S4925 = - _S4924;
    float2  _S4926 = _S4920.differential_0 + make_float2 (_S4747 * _S4921.differential_0, _S4749 * _S4922);
    float2  _S4927 = - _S4926;
    float3  _S4928 = _S4811.differential_0 + _S4893.differential_0;
    float3  _S4929 = _S4813 + _S4814 + _S4893.differential_0;
    float3  _S4930 = _S4812.differential_0 + _S4893.differential_0;
    FixedArray<float3 , 2>  _S4931;
    _S4931[int(0)] = _S4805;
    _S4931[int(1)] = _S4805;
    _S4931[int(1)] = _S4820;
    _S4931[int(0)] = _S4825;
    float3  _S4932 = _S4931[int(0)];
    float3  _S4933 = _S4931[int(1)];
    FixedArray<float3 , 16>  _S4934;
    _S4934[int(0)] = _S4805;
    _S4934[int(1)] = _S4805;
    _S4934[int(2)] = _S4805;
    _S4934[int(3)] = _S4805;
    _S4934[int(4)] = _S4805;
    _S4934[int(5)] = _S4805;
    _S4934[int(6)] = _S4805;
    _S4934[int(7)] = _S4805;
    _S4934[int(8)] = _S4805;
    _S4934[int(9)] = _S4805;
    _S4934[int(10)] = _S4805;
    _S4934[int(11)] = _S4805;
    _S4934[int(12)] = _S4805;
    _S4934[int(13)] = _S4805;
    _S4934[int(14)] = _S4805;
    _S4934[int(15)] = _S4805;
    _S4934[int(7)] = _S4849;
    _S4934[int(0)] = _S4883;
    _S4934[int(1)] = _S4870;
    _S4934[int(2)] = _S4868;
    _S4934[int(3)] = _S4866;
    _S4934[int(4)] = _S4855;
    _S4934[int(5)] = _S4853;
    _S4934[int(6)] = _S4851;
    _S4934[int(15)] = _S4827;
    _S4934[int(8)] = _S4847;
    _S4934[int(9)] = _S4839;
    _S4934[int(10)] = _S4837;
    _S4934[int(11)] = _S4835;
    _S4934[int(12)] = _S4833;
    _S4934[int(13)] = _S4831;
    _S4934[int(14)] = _S4829;
    float3  _S4935 = _S4934[int(0)];
    float3  _S4936 = _S4934[int(1)];
    float3  _S4937 = _S4934[int(2)];
    float3  _S4938 = _S4934[int(3)];
    float3  _S4939 = _S4934[int(4)];
    float3  _S4940 = _S4934[int(5)];
    float3  _S4941 = _S4934[int(6)];
    float3  _S4942 = _S4934[int(7)];
    float3  _S4943 = _S4934[int(8)];
    float3  _S4944 = _S4934[int(9)];
    float3  _S4945 = _S4934[int(10)];
    float3  _S4946 = _S4934[int(11)];
    float3  _S4947 = _S4934[int(12)];
    float3  _S4948 = _S4934[int(13)];
    float3  _S4949 = _S4934[int(14)];
    float3  _S4950 = _S4934[int(15)];
    float2  _S4951 = _S4918.differential_0 + _S4927 + make_float2 (_S4915, _S4904);
    float2  _S4952 = _S4923 + _S4924 + make_float2 (_S4911, _S4900);
    float2  _S4953 = _S4925 + _S4926 + make_float2 (_S4914, _S4903);
    float2  _S4954 = make_float2 (0.0f, _S4917);
    if(_S4713)
    {
        float _S4955 = fx_29 * _S4952.x;
        float2  _S4956 = make_float2 (_S4955, fy_29 * _S4952.y) + make_float2 (_S4716 * _S4955, _S4715 * _S4955);
        float2  _S4957 = uv_44 * _S4956;
        float _S4958 = _S4718 * _S4956.y;
        float _S4959 = _S4724 * _S4956.x;
        float _S4960 = _S4957.x + _S4957.y;
        float _S4961 = _S4729 * _S4960;
        float _S4962 = _S4729 * _S4961;
        float _S4963 = _S4717 * _S4956.y + _S4958 + _S4723 * _S4956.x + _S4959 + _S4730 * _S4960 + _S4731 * _S4961 + _S4732 * _S4962 + _S4733 * (_S4729 * _S4962);
        float _S4964 = _S4720 * _S4963;
        float _S4965 = _S4726 * _S4963;
        _S4640 = _S4714 * _S4956 + make_float2 (_S4722 * (_S4720 * _S4956.y) + _S4725 * _S4959 + 2.0f * (_S4726 * _S4959) + _S4728 * (_S4720 * _S4956.x) + _S4965 + _S4965, _S4719 * _S4958 + 2.0f * (_S4720 * _S4958) + _S4721 * _S4956.y + _S4727 * _S4956.x + _S4964 + _S4964);
    }
    else
    {
        _S4640 = _S4636;
    }
    float2  _S4966 = _S4640 / _S4711;
    float2  _S4967 = _S4631 * - _S4966;
    float2  _S4968 = _S4710 * _S4966;
    float3  _S4969 = _S4930 + make_float3 (_S4968.x, _S4968.y, _S4967.x + _S4967.y);
    if(_S4677)
    {
        float _S4970 = fx_29 * _S4953.x;
        float2  _S4971 = make_float2 (_S4970, fy_29 * _S4953.y) + make_float2 (_S4680 * _S4970, _S4679 * _S4970);
        float2  _S4972 = uv_43 * _S4971;
        float _S4973 = _S4682 * _S4971.y;
        float _S4974 = _S4688 * _S4971.x;
        float _S4975 = _S4972.x + _S4972.y;
        float _S4976 = _S4693 * _S4975;
        float _S4977 = _S4693 * _S4976;
        float _S4978 = _S4681 * _S4971.y + _S4973 + _S4687 * _S4971.x + _S4974 + _S4694 * _S4975 + _S4695 * _S4976 + _S4696 * _S4977 + _S4697 * (_S4693 * _S4977);
        float _S4979 = _S4684 * _S4978;
        float _S4980 = _S4690 * _S4978;
        _S4640 = _S4678 * _S4971 + make_float2 (_S4686 * (_S4684 * _S4971.y) + _S4689 * _S4974 + 2.0f * (_S4690 * _S4974) + _S4692 * (_S4684 * _S4971.x) + _S4980 + _S4980, _S4683 * _S4973 + 2.0f * (_S4684 * _S4973) + _S4685 * _S4971.y + _S4691 * _S4971.x + _S4979 + _S4979);
    }
    else
    {
        _S4640 = _S4636;
    }
    float2  _S4981 = _S4640 / _S4675;
    float2  _S4982 = _S4627 * - _S4981;
    float2  _S4983 = _S4674 * _S4981;
    float3  _S4984 = _S4928 + make_float3 (_S4983.x, _S4983.y, _S4982.x + _S4982.y);
    if(_S4641)
    {
        float _S4985 = fx_29 * _S4951.x;
        float2  _S4986 = make_float2 (_S4985, fy_29 * _S4951.y) + make_float2 (_S4644 * _S4985, _S4643 * _S4985);
        float2  _S4987 = uv_42 * _S4986;
        float _S4988 = _S4646 * _S4986.y;
        float _S4989 = _S4652 * _S4986.x;
        float _S4990 = _S4987.x + _S4987.y;
        float _S4991 = _S4657 * _S4990;
        float _S4992 = _S4657 * _S4991;
        float _S4993 = _S4645 * _S4986.y + _S4988 + _S4651 * _S4986.x + _S4989 + _S4658 * _S4990 + _S4659 * _S4991 + _S4660 * _S4992 + _S4661 * (_S4657 * _S4992);
        float _S4994 = _S4648 * _S4993;
        float _S4995 = _S4654 * _S4993;
        _S4640 = _S4642 * _S4986 + make_float2 (_S4650 * (_S4648 * _S4986.y) + _S4653 * _S4989 + 2.0f * (_S4654 * _S4989) + _S4656 * (_S4648 * _S4986.x) + _S4995 + _S4995, _S4647 * _S4988 + 2.0f * (_S4648 * _S4988) + _S4649 * _S4986.y + _S4655 * _S4986.x + _S4994 + _S4994);
    }
    else
    {
        _S4640 = _S4636;
    }
    float2  _S4996 = _S4640 / _S4638;
    float2  _S4997 = _S4623 * - _S4996;
    float2  _S4998 = _S4637 * _S4996;
    float _S4999 = _S4997.x + _S4997.y;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5000;
    (&_S5000)->primal_0 = R_29;
    (&_S5000)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5001;
    (&_S5001)->primal_0 = vert2_6;
    (&_S5001)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5000, &_S5001, _S4969);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5002;
    (&_S5002)->primal_0 = R_29;
    (&_S5002)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5003;
    (&_S5003)->primal_0 = vert1_6;
    (&_S5003)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5002, &_S5003, _S4984);
    float3  _S5004 = _S4929 + make_float3 (_S4998.x, _S4998.y, _S4999);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5005;
    (&_S5005)->primal_0 = R_29;
    (&_S5005)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5006;
    (&_S5006)->primal_0 = vert0_6;
    (&_S5006)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5005, &_S5006, _S5004);
    float3  _S5007 = _S5001.differential_0 + (*v_verts_0)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5008;
    (&_S5008)->primal_0 = _S4617;
    (&_S5008)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5009;
    (&_S5009)->primal_0 = _S4622;
    (&_S5009)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5008, &_S5009, _S5007);
    float _S5010 = - _S5009.differential_0.y;
    float _S5011 = _S4621 * _S5009.differential_0.x;
    float _S5012 = - (_S4613 * _S5009.differential_0.x);
    float3  _S5013 = _S5003.differential_0 + (*v_verts_0)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5014;
    (&_S5014)->primal_0 = _S4617;
    (&_S5014)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5015;
    (&_S5015)->primal_0 = _S4620;
    (&_S5015)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5014, &_S5015, _S5013);
    float _S5016 = _S4613 * _S5015.differential_0.x;
    float _S5017 = _S4619 * _S5015.differential_0.x;
    float3  _S5018 = _S5006.differential_0 + (*v_verts_0)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5019;
    (&_S5019)->primal_0 = _S4617;
    (&_S5019)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5020;
    (&_S5020)->primal_0 = _S4618;
    (&_S5020)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5019, &_S5020, _S5018);
    Matrix<float, 3, 3>  _S5021 = transpose_0(_S5008.differential_0 + _S5014.differential_0 + _S5019.differential_0);
    float _S5022 = 2.0f * - _S5021.rows[int(2)].z;
    float _S5023 = 2.0f * _S5021.rows[int(2)].y;
    float _S5024 = 2.0f * _S5021.rows[int(2)].x;
    float _S5025 = 2.0f * _S5021.rows[int(1)].z;
    float _S5026 = 2.0f * - _S5021.rows[int(1)].y;
    float _S5027 = 2.0f * _S5021.rows[int(1)].x;
    float _S5028 = 2.0f * _S5021.rows[int(0)].z;
    float _S5029 = 2.0f * _S5021.rows[int(0)].y;
    float _S5030 = 2.0f * - _S5021.rows[int(0)].x;
    float _S5031 = - _S5027 + _S5029;
    float _S5032 = _S5024 + - _S5028;
    float _S5033 = - _S5023 + _S5025;
    float _S5034 = _S5023 + _S5025;
    float _S5035 = _S5024 + _S5028;
    float _S5036 = _S5027 + _S5029;
    float _S5037 = quat_30.w * (_S5026 + _S5030);
    float _S5038 = quat_30.z * (_S5022 + _S5030);
    float _S5039 = quat_30.y * (_S5022 + _S5026);
    float _S5040 = quat_30.x * _S5031 + quat_30.z * _S5034 + quat_30.y * _S5035 + _S5037 + _S5037;
    float _S5041 = quat_30.x * _S5032 + quat_30.w * _S5034 + quat_30.y * _S5036 + _S5038 + _S5038;
    float _S5042 = quat_30.x * _S5033 + quat_30.w * _S5035 + quat_30.z * _S5036 + _S5039 + _S5039;
    float _S5043 = quat_30.w * _S5031 + quat_30.z * _S5032 + quat_30.y * _S5033;
    float _S5044 = _S5012 + _S5016;
    float _S5045 = 0.5f * - _S5044;
    float _S5046 = _S5010 + _S5015.differential_0.y;
    DiffPair_float_0 _S5047;
    (&_S5047)->primal_0 = _S4614;
    (&_S5047)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5047, _S5046);
    float _S5048 = _S5045 + _S5047.differential_0;
    float _S5049 = _S5011 + _S5017 + _S5020.differential_0.x;
    DiffPair_float_0 _S5050;
    (&_S5050)->primal_0 = _S4612;
    (&_S5050)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5050, _S5049);
    float _S5051 = _S5045 + _S5050.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5052;
    (&_S5052)->primal_0 = R_29;
    (&_S5052)->differential_0 = _S4886;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5053;
    (&_S5053)->primal_0 = mean_29;
    (&_S5053)->differential_0 = _S4805;
    s_bwd_prop_mul_1(&_S5052, &_S5053, _S4808.differential_0);
    float3  _S5054 = _S4969 + _S4984 + _S5004 + _S4808.differential_0 + _S4889.differential_0;
    Matrix<float, 3, 3>  _S5055 = _S5000.differential_0 + _S5002.differential_0 + _S5005.differential_0 + _S5052.differential_0 + _S4890;
    float3  _S5056 = make_float3 (_S5051, _S5048, _S5044);
    float4  _S5057 = make_float4 (0.0f);
    *&((&_S5057)->w) = _S5040;
    *&((&_S5057)->z) = _S5041;
    *&((&_S5057)->y) = _S5042;
    *&((&_S5057)->x) = _S5043;
    *v_mean_9 = _S5007 + _S5013 + _S5018 + _S5053.differential_0 + _S4884;
    *v_quat_8 = _S5057;
    *v_scale_8 = _S5056;
    *v_hardness_4 = _S4954;
    (*v_sh_coeffs_7)[int(0)] = _S4935;
    (*v_sh_coeffs_7)[int(1)] = _S4936;
    (*v_sh_coeffs_7)[int(2)] = _S4937;
    (*v_sh_coeffs_7)[int(3)] = _S4938;
    (*v_sh_coeffs_7)[int(4)] = _S4939;
    (*v_sh_coeffs_7)[int(5)] = _S4940;
    (*v_sh_coeffs_7)[int(6)] = _S4941;
    (*v_sh_coeffs_7)[int(7)] = _S4942;
    (*v_sh_coeffs_7)[int(8)] = _S4943;
    (*v_sh_coeffs_7)[int(9)] = _S4944;
    (*v_sh_coeffs_7)[int(10)] = _S4945;
    (*v_sh_coeffs_7)[int(11)] = _S4946;
    (*v_sh_coeffs_7)[int(12)] = _S4947;
    (*v_sh_coeffs_7)[int(13)] = _S4948;
    (*v_sh_coeffs_7)[int(14)] = _S4949;
    (*v_sh_coeffs_7)[int(15)] = _S4950;
    (*v_ch_coeffs_2)[int(0)] = _S4932;
    (*v_ch_coeffs_2)[int(1)] = _S4933;
    *v_R_8 = _S5055;
    *v_t_8 = _S5054;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_30, float fy_30, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_30) + t_29;
    float _S5058 = scale_30.x;
    float _S5059 = s_primal_ctx_exp_1(_S5058);
    float _S5060 = scale_30.y;
    float _S5061 = s_primal_ctx_exp_1(_S5060);
    float sz_11 = scale_30.z - 0.5f * (_S5058 + _S5060);
    float _S5062 = quat_31.y;
    float x2_31 = _S5062 * _S5062;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_57 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S5063 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_57), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_57), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S5064 = make_float3 (_S5059, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S5063, _S5064) + mean_30;
    float _S5065 = -0.5f + sz_11;
    float3  _S5066 = make_float3 (_S5059 * _S5065, _S5061, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S5063, _S5066) + mean_30;
    float _S5067 = -0.5f - sz_11;
    float3  _S5068 = make_float3 (_S5059 * _S5067, - _S5061, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S5063, _S5068) + mean_30;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_7) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_7) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_7) + t_29;
    float2  _S5069 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S5070 = length_0(_S5069);
    float _S5071 = vert0_c_11.z;
    float _S5072 = s_primal_ctx_atan2_0(_S5070, _S5071);
    bool _S5073 = _S5072 < 0.00100000004749745f;
    float k_9;
    float _S5074;
    float _S5075;
    float _S5076;
    if(_S5073)
    {
        float _S5077 = 1.0f - _S5072 * _S5072 / 3.0f;
        float _S5078 = _S5071 * _S5071;
        k_9 = _S5077 / _S5071;
        _S5074 = 0.0f;
        _S5075 = _S5078;
        _S5076 = _S5077;
    }
    else
    {
        float _S5079 = _S5070 * _S5070;
        k_9 = _S5072 / _S5070;
        _S5074 = _S5079;
        _S5075 = 0.0f;
        _S5076 = 0.0f;
    }
    float2  _S5080 = make_float2 (k_9);
    float2  _S5081 = _S5069 * make_float2 (k_9);
    float u_55 = _S5081.x;
    float v_55 = _S5081.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float _S5082 = (*dist_coeffs_36)[int(2)] + r2_55 * (*dist_coeffs_36)[int(3)];
    float _S5083 = (*dist_coeffs_36)[int(1)] + r2_55 * _S5082;
    float _S5084 = (*dist_coeffs_36)[int(0)] + r2_55 * _S5083;
    float radial_12 = 1.0f + r2_55 * _S5084;
    float2  _S5085 = make_float2 (radial_12);
    float _S5086 = 2.0f * (*dist_coeffs_36)[int(4)];
    float _S5087 = _S5086 * u_55;
    float _S5088 = 2.0f * u_55;
    float _S5089 = 2.0f * (*dist_coeffs_36)[int(5)];
    float _S5090 = _S5089 * u_55;
    float _S5091 = 2.0f * v_55;
    float2  _S5092 = _S5081 * make_float2 (radial_12) + make_float2 (_S5087 * v_55 + (*dist_coeffs_36)[int(5)] * (r2_55 + _S5088 * u_55) + (*dist_coeffs_36)[int(6)] * r2_55, _S5090 * v_55 + (*dist_coeffs_36)[int(4)] * (r2_55 + _S5091 * v_55) + (*dist_coeffs_36)[int(7)] * r2_55);
    float2  _S5093 = _S5092 + make_float2 ((*dist_coeffs_36)[int(8)] * _S5092.x + (*dist_coeffs_36)[int(9)] * _S5092.y, 0.0f);
    float _S5094 = fx_30 * _S5093.x + cx_30;
    float _S5095 = fy_30 * _S5093.y + cy_30;
    float2  _S5096 = make_float2 (_S5094, _S5095);
    float2  _S5097 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S5098 = length_0(_S5097);
    float _S5099 = vert1_c_11.z;
    float _S5100 = s_primal_ctx_atan2_0(_S5098, _S5099);
    bool _S5101 = _S5100 < 0.00100000004749745f;
    float _S5102;
    float _S5103;
    float _S5104;
    if(_S5101)
    {
        float _S5105 = 1.0f - _S5100 * _S5100 / 3.0f;
        float _S5106 = _S5099 * _S5099;
        k_9 = _S5105 / _S5099;
        _S5102 = 0.0f;
        _S5103 = _S5106;
        _S5104 = _S5105;
    }
    else
    {
        float _S5107 = _S5098 * _S5098;
        k_9 = _S5100 / _S5098;
        _S5102 = _S5107;
        _S5103 = 0.0f;
        _S5104 = 0.0f;
    }
    float2  _S5108 = make_float2 (k_9);
    float2  _S5109 = _S5097 * make_float2 (k_9);
    float u_56 = _S5109.x;
    float v_56 = _S5109.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S5110 = (*dist_coeffs_36)[int(2)] + r2_56 * (*dist_coeffs_36)[int(3)];
    float _S5111 = (*dist_coeffs_36)[int(1)] + r2_56 * _S5110;
    float _S5112 = (*dist_coeffs_36)[int(0)] + r2_56 * _S5111;
    float radial_13 = 1.0f + r2_56 * _S5112;
    float2  _S5113 = make_float2 (radial_13);
    float _S5114 = _S5086 * u_56;
    float _S5115 = 2.0f * u_56;
    float _S5116 = _S5089 * u_56;
    float _S5117 = 2.0f * v_56;
    float2  _S5118 = _S5109 * make_float2 (radial_13) + make_float2 (_S5114 * v_56 + (*dist_coeffs_36)[int(5)] * (r2_56 + _S5115 * u_56) + (*dist_coeffs_36)[int(6)] * r2_56, _S5116 * v_56 + (*dist_coeffs_36)[int(4)] * (r2_56 + _S5117 * v_56) + (*dist_coeffs_36)[int(7)] * r2_56);
    float2  _S5119 = _S5118 + make_float2 ((*dist_coeffs_36)[int(8)] * _S5118.x + (*dist_coeffs_36)[int(9)] * _S5118.y, 0.0f);
    float _S5120 = fx_30 * _S5119.x + cx_30;
    float _S5121 = fy_30 * _S5119.y + cy_30;
    float2  _S5122 = make_float2 (_S5120, _S5121);
    float2  _S5123 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S5124 = length_0(_S5123);
    float _S5125 = vert2_c_11.z;
    float _S5126 = s_primal_ctx_atan2_0(_S5124, _S5125);
    bool _S5127 = _S5126 < 0.00100000004749745f;
    float _S5128;
    float _S5129;
    float _S5130;
    if(_S5127)
    {
        float _S5131 = 1.0f - _S5126 * _S5126 / 3.0f;
        float _S5132 = _S5125 * _S5125;
        k_9 = _S5131 / _S5125;
        _S5128 = 0.0f;
        _S5129 = _S5132;
        _S5130 = _S5131;
    }
    else
    {
        float _S5133 = _S5124 * _S5124;
        k_9 = _S5126 / _S5124;
        _S5128 = _S5133;
        _S5129 = 0.0f;
        _S5130 = 0.0f;
    }
    float2  _S5134 = make_float2 (k_9);
    float2  _S5135 = _S5123 * make_float2 (k_9);
    float u_57 = _S5135.x;
    float v_57 = _S5135.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S5136 = (*dist_coeffs_36)[int(2)] + r2_57 * (*dist_coeffs_36)[int(3)];
    float _S5137 = (*dist_coeffs_36)[int(1)] + r2_57 * _S5136;
    float _S5138 = (*dist_coeffs_36)[int(0)] + r2_57 * _S5137;
    float radial_14 = 1.0f + r2_57 * _S5138;
    float2  _S5139 = make_float2 (radial_14);
    float _S5140 = _S5086 * u_57;
    float _S5141 = 2.0f * u_57;
    float _S5142 = _S5089 * u_57;
    float _S5143 = 2.0f * v_57;
    float2  _S5144 = _S5135 * make_float2 (radial_14) + make_float2 (_S5140 * v_57 + (*dist_coeffs_36)[int(5)] * (r2_57 + _S5141 * u_57) + (*dist_coeffs_36)[int(6)] * r2_57, _S5142 * v_57 + (*dist_coeffs_36)[int(4)] * (r2_57 + _S5143 * v_57) + (*dist_coeffs_36)[int(7)] * r2_57);
    float2  _S5145 = _S5144 + make_float2 ((*dist_coeffs_36)[int(8)] * _S5144.x + (*dist_coeffs_36)[int(9)] * _S5144.y, 0.0f);
    float _S5146 = fx_30 * _S5145.x + cx_30;
    float _S5147 = fy_30 * _S5145.y + cy_30;
    float2  _S5148 = make_float2 (_S5146, _S5147);
    float2  e0_15 = _S5122 - _S5096;
    float2  e1_15 = _S5148 - _S5122;
    float2  e2_7 = _S5096 - _S5148;
    float _S5149 = e0_15.x;
    float _S5150 = e1_15.y;
    float _S5151 = e0_15.y;
    float _S5152 = e1_15.x;
    float _S5153 = _S5149 * _S5150 - _S5151 * _S5152;
    float _S5154 = 1.0f - hardness_15.y;
    float _S5155 = -1.0f / _S5154;
    float _S5156 = _S5154 * _S5154;
    float _S5157 = s_primal_ctx_max_0(_S5094, _S5120);
    float _S5158 = s_primal_ctx_min_0(_S5094, _S5120);
    float _S5159 = s_primal_ctx_max_0(_S5095, _S5121);
    float _S5160 = s_primal_ctx_min_0(_S5095, _S5121);
    float3  _S5161 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S5162 = length_1(_S5161) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S5163 = transpose_0(R_30);
    float3  _S5164 = mean_30 - - s_primal_ctx_mul_1(_S5163, t_29);
    float _S5165 = _S5164.x;
    float _S5166 = _S5164.y;
    float _S5167 = _S5164.z;
    float _S5168 = _S5165 * _S5165 + _S5166 * _S5166 + _S5167 * _S5167;
    float _S5169 = s_primal_ctx_sqrt_0(_S5168);
    float x_61 = _S5165 / _S5169;
    float3  _S5170 = make_float3 (x_61);
    float _S5171 = _S5169 * _S5169;
    float y_29 = _S5166 / _S5169;
    float z_26 = _S5167 / _S5169;
    float3  _S5172 = make_float3 (z_26);
    float _S5173 = - y_29;
    float3  _S5174 = make_float3 (_S5173);
    float z2_58 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_61 * x_61 - y_29 * y_29;
    float _S5175 = 2.0f * x_61;
    float fS1_26 = _S5175 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_58 - 0.31539157032966614f;
    float3  _S5176 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_61;
    float3  _S5177 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S5178 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S5179 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S5180 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_58 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S5181 = 1.86588168144226074f * z2_58 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S5181;
    float3  _S5182 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_61;
    float3  _S5183 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S5184 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S5185 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S5186 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_61 * fC1_26 - y_29 * fS1_26);
    float3  _S5187 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_61 * fS1_26 + y_29 * fC1_26);
    float3  _S5188 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S5173) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_61) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S5189 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S5190 = make_float3 (0.0f);
    float3  _S5191 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S5192 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S5193 = make_float3 (_S5192);
    float3  _S5194 = (*ch_coeffs_11)[int(1)] * make_float3 (_S5192);
    float3  _S5195 = _S5191 + _S5194 + make_float3 (0.5f);
    float3  _S5196 = _S5191 - _S5194 + make_float3 (0.5f);
    float3  _S5197 = vert1_c_11 - vert0_c_11;
    float3  _S5198 = vert2_c_11 - vert0_c_11;
    float3  _S5199 = s_primal_ctx_cross_0(_S5197, _S5198);
    float3  _S5200 = normalize_0(_S5199);
    float3  _S5201 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S5200, mean_c_26)))))) * v_normal_3;
    float3  _S5202 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5203;
    (&_S5203)->primal_0 = _S5200;
    (&_S5203)->differential_0 = _S5202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5204;
    (&_S5204)->primal_0 = mean_c_26;
    (&_S5204)->differential_0 = _S5202;
    s_bwd_prop_dot_0(&_S5203, &_S5204, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5205 = _S5204;
    float3  _S5206 = _S5201 + _S5203.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5207;
    (&_S5207)->primal_0 = _S5199;
    (&_S5207)->differential_0 = _S5202;
    s_bwd_normalize_impl_0(&_S5207, _S5206);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5208;
    (&_S5208)->primal_0 = _S5197;
    (&_S5208)->differential_0 = _S5202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5209;
    (&_S5209)->primal_0 = _S5198;
    (&_S5209)->differential_0 = _S5202;
    s_bwd_prop_cross_0(&_S5208, &_S5209, _S5207.differential_0);
    float3  _S5210 = - _S5209.differential_0;
    float3  _S5211 = - _S5208.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5212;
    (&_S5212)->primal_0 = _S5196;
    (&_S5212)->differential_0 = _S5202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5213;
    (&_S5213)->primal_0 = _S5190;
    (&_S5213)->differential_0 = _S5202;
    s_bwd_prop_max_0(&_S5212, &_S5213, (*v_rgbs_1)[int(2)]);
    float3  _S5214 = - _S5212.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5215;
    (&_S5215)->primal_0 = _S5195;
    (&_S5215)->differential_0 = _S5202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5216;
    (&_S5216)->primal_0 = _S5190;
    (&_S5216)->differential_0 = _S5202;
    s_bwd_prop_max_0(&_S5215, &_S5216, (*v_rgbs_1)[int(1)]);
    float3  _S5217 = _S5193 * (_S5214 + _S5215.differential_0);
    float3  _S5218 = _S5212.differential_0 + _S5215.differential_0;
    float3  _S5219 = make_float3 (0.5f) * - _S5218;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5220;
    (&_S5220)->primal_0 = _S5189;
    (&_S5220)->differential_0 = _S5202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5221;
    (&_S5221)->primal_0 = _S5190;
    (&_S5221)->differential_0 = _S5202;
    s_bwd_prop_max_0(&_S5220, &_S5221, (*v_rgbs_1)[int(0)]);
    float3  _S5222 = _S5219 + _S5220.differential_0;
    float3  _S5223 = _S5218 + _S5220.differential_0;
    float3  _S5224 = _S5187 * _S5223;
    float3  _S5225 = (*sh_coeffs_26)[int(15)] * _S5223;
    float3  _S5226 = _S5185 * _S5223;
    float3  _S5227 = (*sh_coeffs_26)[int(14)] * _S5223;
    float3  _S5228 = _S5183 * _S5223;
    float3  _S5229 = (*sh_coeffs_26)[int(13)] * _S5223;
    float3  _S5230 = _S5182 * _S5223;
    float3  _S5231 = (*sh_coeffs_26)[int(12)] * _S5223;
    float3  _S5232 = _S5184 * _S5223;
    float3  _S5233 = (*sh_coeffs_26)[int(11)] * _S5223;
    float3  _S5234 = _S5186 * _S5223;
    float3  _S5235 = (*sh_coeffs_26)[int(10)] * _S5223;
    float3  _S5236 = _S5188 * _S5223;
    float3  _S5237 = (*sh_coeffs_26)[int(9)] * _S5223;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S5237.x + _S5237.y + _S5237.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S5225.x + _S5225.y + _S5225.z);
    float _S5238 = _S5235.x + _S5235.y + _S5235.z;
    float _S5239 = _S5227.x + _S5227.y + _S5227.z;
    float _S5240 = _S5233.x + _S5233.y + _S5233.z;
    float _S5241 = _S5229.x + _S5229.y + _S5229.z;
    float _S5242 = _S5231.x + _S5231.y + _S5231.z;
    float _S5243 = - s_diff_fC2_T_8;
    float3  _S5244 = _S5179 * _S5223;
    float3  _S5245 = (*sh_coeffs_26)[int(8)] * _S5223;
    float3  _S5246 = _S5177 * _S5223;
    float3  _S5247 = (*sh_coeffs_26)[int(7)] * _S5223;
    float3  _S5248 = _S5176 * _S5223;
    float3  _S5249 = (*sh_coeffs_26)[int(6)] * _S5223;
    float3  _S5250 = _S5178 * _S5223;
    float3  _S5251 = (*sh_coeffs_26)[int(5)] * _S5223;
    float3  _S5252 = _S5180 * _S5223;
    float3  _S5253 = (*sh_coeffs_26)[int(4)] * _S5223;
    float _S5254 = _S5251.x + _S5251.y + _S5251.z;
    float _S5255 = _S5247.x + _S5247.y + _S5247.z;
    float _S5256 = fTmp1B_26 * _S5238 + x_61 * s_diff_fS2_T_8 + y_29 * _S5243 + 0.54627424478530884f * (_S5253.x + _S5253.y + _S5253.z);
    float _S5257 = fTmp1B_26 * _S5239 + y_29 * s_diff_fS2_T_8 + x_61 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S5245.x + _S5245.y + _S5245.z);
    float _S5258 = y_29 * - _S5257;
    float _S5259 = x_61 * _S5257;
    float _S5260 = z_26 * (1.86588168144226074f * (z_26 * _S5242) + -2.28522896766662598f * (y_29 * _S5240 + x_61 * _S5241) + 0.94617468118667603f * (_S5249.x + _S5249.y + _S5249.z));
    float3  _S5261 = make_float3 (0.48860251903533936f) * _S5223;
    float3  _S5262 = - _S5261;
    float3  _S5263 = _S5170 * _S5262;
    float3  _S5264 = (*sh_coeffs_26)[int(3)] * _S5262;
    float3  _S5265 = _S5172 * _S5261;
    float3  _S5266 = (*sh_coeffs_26)[int(2)] * _S5261;
    float3  _S5267 = _S5174 * _S5261;
    float3  _S5268 = (*sh_coeffs_26)[int(1)] * _S5261;
    float _S5269 = (_S5181 * _S5242 + 1.44530570507049561f * (fS1_26 * _S5238 + fC1_26 * _S5239) + -1.09254848957061768f * (y_29 * _S5254 + x_61 * _S5255) + _S5260 + _S5260 + _S5266.x + _S5266.y + _S5266.z) / _S5171;
    float _S5270 = _S5169 * _S5269;
    float _S5271 = (fTmp0C_26 * _S5240 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S5243 + fTmp0B_26 * _S5254 + _S5175 * _S5256 + _S5258 + _S5258 + - (_S5268.x + _S5268.y + _S5268.z)) / _S5171;
    float _S5272 = _S5169 * _S5271;
    float _S5273 = (fTmp0C_26 * _S5241 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S5255 + 2.0f * (y_29 * _S5256) + _S5259 + _S5259 + _S5264.x + _S5264.y + _S5264.z) / _S5171;
    float _S5274 = _S5169 * _S5273;
    float _S5275 = _S5167 * - _S5269 + _S5166 * - _S5271 + _S5165 * - _S5273;
    DiffPair_float_0 _S5276;
    (&_S5276)->primal_0 = _S5168;
    (&_S5276)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5276, _S5275);
    float _S5277 = _S5167 * _S5276.differential_0;
    float _S5278 = _S5166 * _S5276.differential_0;
    float _S5279 = _S5165 * _S5276.differential_0;
    float3  _S5280 = make_float3 (0.282094806432724f) * _S5223;
    float3  _S5281 = make_float3 (_S5274 + _S5279 + _S5279, _S5272 + _S5278 + _S5278, _S5270 + _S5277 + _S5277);
    float3  _S5282 = - - _S5281;
    Matrix<float, 3, 3>  _S5283 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5284;
    (&_S5284)->primal_0 = _S5163;
    (&_S5284)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5285;
    (&_S5285)->primal_0 = t_29;
    (&_S5285)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5284, &_S5285, _S5282);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5286 = _S5285;
    Matrix<float, 3, 3>  _S5287 = transpose_0(_S5284.differential_0);
    DiffPair_float_0 _S5288;
    (&_S5288)->primal_0 = _S5162;
    (&_S5288)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5288, v_depth_10);
    float _S5289 = 0.3333333432674408f * _S5288.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5290;
    (&_S5290)->primal_0 = _S5161;
    (&_S5290)->differential_0 = _S5202;
    s_bwd_length_impl_0(&_S5290, _S5289);
    DiffPair_float_0 _S5291;
    (&_S5291)->primal_0 = _S5160;
    (&_S5291)->differential_0 = 0.0f;
    DiffPair_float_0 _S5292;
    (&_S5292)->primal_0 = _S5147;
    (&_S5292)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5291, &_S5292, 0.0f);
    DiffPair_float_0 _S5293;
    (&_S5293)->primal_0 = _S5095;
    (&_S5293)->differential_0 = 0.0f;
    DiffPair_float_0 _S5294;
    (&_S5294)->primal_0 = _S5121;
    (&_S5294)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5293, &_S5294, _S5291.differential_0);
    DiffPair_float_0 _S5295;
    (&_S5295)->primal_0 = _S5159;
    (&_S5295)->differential_0 = 0.0f;
    DiffPair_float_0 _S5296;
    (&_S5296)->primal_0 = _S5147;
    (&_S5296)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5295, &_S5296, 0.0f);
    DiffPair_float_0 _S5297;
    (&_S5297)->primal_0 = _S5095;
    (&_S5297)->differential_0 = 0.0f;
    DiffPair_float_0 _S5298;
    (&_S5298)->primal_0 = _S5121;
    (&_S5298)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5297, &_S5298, _S5295.differential_0);
    DiffPair_float_0 _S5299;
    (&_S5299)->primal_0 = _S5158;
    (&_S5299)->differential_0 = 0.0f;
    DiffPair_float_0 _S5300;
    (&_S5300)->primal_0 = _S5146;
    (&_S5300)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5299, &_S5300, 0.0f);
    DiffPair_float_0 _S5301;
    (&_S5301)->primal_0 = _S5094;
    (&_S5301)->differential_0 = 0.0f;
    DiffPair_float_0 _S5302;
    (&_S5302)->primal_0 = _S5120;
    (&_S5302)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5301, &_S5302, _S5299.differential_0);
    DiffPair_float_0 _S5303;
    (&_S5303)->primal_0 = _S5157;
    (&_S5303)->differential_0 = 0.0f;
    DiffPair_float_0 _S5304;
    (&_S5304)->primal_0 = _S5146;
    (&_S5304)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5303, &_S5304, 0.0f);
    DiffPair_float_0 _S5305;
    (&_S5305)->primal_0 = _S5094;
    (&_S5305)->differential_0 = 0.0f;
    DiffPair_float_0 _S5306;
    (&_S5306)->primal_0 = _S5120;
    (&_S5306)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5305, &_S5306, _S5303.differential_0);
    DiffPair_float_0 _S5307;
    (&_S5307)->primal_0 = _S5155;
    (&_S5307)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5307, 0.0f);
    float _S5308 = - (-1.0f * - (_S5307.differential_0 / _S5156));
    float2  _S5309 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5310;
    (&_S5310)->primal_0 = e2_7;
    (&_S5310)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5310, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5311;
    (&_S5311)->primal_0 = e1_15;
    (&_S5311)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5311, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5312;
    (&_S5312)->primal_0 = e0_15;
    (&_S5312)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5312, -0.0f);
    DiffPair_float_0 _S5313;
    (&_S5313)->primal_0 = _S5153;
    (&_S5313)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5313, 0.0f);
    float _S5314 = - _S5313.differential_0;
    float2  _S5315 = _S5311.differential_0 + make_float2 (_S5151 * _S5314, _S5149 * _S5313.differential_0);
    float2  _S5316 = - _S5315;
    float2  _S5317 = _S5312.differential_0 + make_float2 (_S5150 * _S5313.differential_0, _S5152 * _S5314);
    float2  _S5318 = - _S5317;
    float2  _S5319 = - _S5310.differential_0 + _S5315;
    float _S5320 = fx_30 * (_S5300.differential_0 + _S5304.differential_0 + _S5319.x);
    float2  _S5321 = make_float2 (_S5320, fy_30 * (_S5292.differential_0 + _S5296.differential_0 + _S5319.y)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S5320, (*dist_coeffs_36)[int(9)] * _S5320);
    float2  _S5322 = _S5135 * _S5321;
    float2  _S5323 = _S5139 * _S5321;
    float _S5324 = (*dist_coeffs_36)[int(4)] * _S5321.y;
    float _S5325 = (*dist_coeffs_36)[int(5)] * _S5321.x;
    float _S5326 = _S5322.x + _S5322.y;
    float _S5327 = r2_57 * _S5326;
    float _S5328 = r2_57 * _S5327;
    float _S5329 = (*dist_coeffs_36)[int(7)] * _S5321.y + _S5324 + (*dist_coeffs_36)[int(6)] * _S5321.x + _S5325 + _S5138 * _S5326 + _S5137 * _S5327 + _S5136 * _S5328 + (*dist_coeffs_36)[int(3)] * (r2_57 * _S5328);
    float _S5330 = v_57 * _S5329;
    float _S5331 = u_57 * _S5329;
    float _S5332 = _S5143 * _S5324 + 2.0f * (v_57 * _S5324) + _S5142 * _S5321.y + _S5140 * _S5321.x + _S5330 + _S5330;
    float _S5333 = _S5089 * (v_57 * _S5321.y) + _S5141 * _S5325 + 2.0f * (u_57 * _S5325) + _S5086 * (v_57 * _S5321.x) + _S5331 + _S5331;
    float3  _S5334 = _S5208.differential_0 + _S5290.differential_0;
    float3  _S5335 = _S5210 + _S5211 + _S5290.differential_0;
    float3  _S5336 = _S5209.differential_0 + _S5290.differential_0;
    FixedArray<float3 , 2>  _S5337;
    _S5337[int(0)] = _S5202;
    _S5337[int(1)] = _S5202;
    _S5337[int(1)] = _S5217;
    _S5337[int(0)] = _S5222;
    float3  _S5338 = _S5337[int(0)];
    float3  _S5339 = _S5337[int(1)];
    FixedArray<float3 , 16>  _S5340;
    _S5340[int(0)] = _S5202;
    _S5340[int(1)] = _S5202;
    _S5340[int(2)] = _S5202;
    _S5340[int(3)] = _S5202;
    _S5340[int(4)] = _S5202;
    _S5340[int(5)] = _S5202;
    _S5340[int(6)] = _S5202;
    _S5340[int(7)] = _S5202;
    _S5340[int(8)] = _S5202;
    _S5340[int(9)] = _S5202;
    _S5340[int(10)] = _S5202;
    _S5340[int(11)] = _S5202;
    _S5340[int(12)] = _S5202;
    _S5340[int(13)] = _S5202;
    _S5340[int(14)] = _S5202;
    _S5340[int(15)] = _S5202;
    _S5340[int(7)] = _S5246;
    _S5340[int(0)] = _S5280;
    _S5340[int(1)] = _S5267;
    _S5340[int(2)] = _S5265;
    _S5340[int(3)] = _S5263;
    _S5340[int(4)] = _S5252;
    _S5340[int(5)] = _S5250;
    _S5340[int(6)] = _S5248;
    _S5340[int(15)] = _S5224;
    _S5340[int(8)] = _S5244;
    _S5340[int(9)] = _S5236;
    _S5340[int(10)] = _S5234;
    _S5340[int(11)] = _S5232;
    _S5340[int(12)] = _S5230;
    _S5340[int(13)] = _S5228;
    _S5340[int(14)] = _S5226;
    float3  _S5341 = _S5340[int(0)];
    float3  _S5342 = _S5340[int(1)];
    float3  _S5343 = _S5340[int(2)];
    float3  _S5344 = _S5340[int(3)];
    float3  _S5345 = _S5340[int(4)];
    float3  _S5346 = _S5340[int(5)];
    float3  _S5347 = _S5340[int(6)];
    float3  _S5348 = _S5340[int(7)];
    float3  _S5349 = _S5340[int(8)];
    float3  _S5350 = _S5340[int(9)];
    float3  _S5351 = _S5340[int(10)];
    float3  _S5352 = _S5340[int(11)];
    float3  _S5353 = _S5340[int(12)];
    float3  _S5354 = _S5340[int(13)];
    float3  _S5355 = _S5340[int(14)];
    float3  _S5356 = _S5340[int(15)];
    float _S5357 = _S5301.differential_0 + _S5305.differential_0;
    float2  _S5358 = _S5310.differential_0 + _S5318;
    float _S5359 = _S5294.differential_0 + _S5298.differential_0;
    float _S5360 = _S5293.differential_0 + _S5297.differential_0;
    float2  _S5361 = _S5316 + _S5317;
    float _S5362 = _S5302.differential_0 + _S5306.differential_0;
    float2  _S5363 = make_float2 (0.0f, _S5308);
    float2  _S5364 = _S5323 + make_float2 (_S5333, _S5332);
    float2  _S5365 = _S5123 * _S5364;
    float2  _S5366 = _S5134 * _S5364;
    float _S5367 = _S5365.x + _S5365.y;
    if(_S5127)
    {
        float _S5368 = _S5367 / _S5129;
        float _S5369 = _S5130 * - _S5368;
        float _S5370 = _S5126 * (0.3333333432674408f * - (_S5125 * _S5368));
        k_9 = _S5370 + _S5370;
        _S5128 = _S5369;
        _S5129 = 0.0f;
    }
    else
    {
        float _S5371 = _S5367 / _S5128;
        float _S5372 = _S5126 * - _S5371;
        k_9 = _S5124 * _S5371;
        _S5128 = 0.0f;
        _S5129 = _S5372;
    }
    DiffPair_float_0 _S5373;
    (&_S5373)->primal_0 = _S5124;
    (&_S5373)->differential_0 = 0.0f;
    DiffPair_float_0 _S5374;
    (&_S5374)->primal_0 = _S5125;
    (&_S5374)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5373, &_S5374, k_9);
    float _S5375 = _S5374.differential_0 + _S5128;
    float _S5376 = _S5373.differential_0 + _S5129;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5377;
    (&_S5377)->primal_0 = _S5123;
    (&_S5377)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5377, _S5376);
    float2  _S5378 = _S5377.differential_0 + _S5366;
    float _S5379 = fx_30 * (_S5361.x + _S5362);
    float2  _S5380 = make_float2 (_S5379, fy_30 * (_S5361.y + _S5359)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S5379, (*dist_coeffs_36)[int(9)] * _S5379);
    float2  _S5381 = _S5109 * _S5380;
    float _S5382 = (*dist_coeffs_36)[int(4)] * _S5380.y;
    float _S5383 = (*dist_coeffs_36)[int(5)] * _S5380.x;
    float _S5384 = _S5381.x + _S5381.y;
    float _S5385 = r2_56 * _S5384;
    float _S5386 = r2_56 * _S5385;
    float _S5387 = (*dist_coeffs_36)[int(7)] * _S5380.y + _S5382 + (*dist_coeffs_36)[int(6)] * _S5380.x + _S5383 + _S5112 * _S5384 + _S5111 * _S5385 + _S5110 * _S5386 + (*dist_coeffs_36)[int(3)] * (r2_56 * _S5386);
    float _S5388 = v_56 * _S5387;
    float _S5389 = u_56 * _S5387;
    float3  _S5390 = _S5336 + make_float3 (_S5378.x, _S5378.y, _S5375);
    float2  _S5391 = _S5113 * _S5380 + make_float2 (_S5089 * (v_56 * _S5380.y) + _S5115 * _S5383 + 2.0f * (u_56 * _S5383) + _S5086 * (v_56 * _S5380.x) + _S5389 + _S5389, _S5117 * _S5382 + 2.0f * (v_56 * _S5382) + _S5116 * _S5380.y + _S5114 * _S5380.x + _S5388 + _S5388);
    float2  _S5392 = _S5097 * _S5391;
    float2  _S5393 = _S5108 * _S5391;
    float _S5394 = _S5392.x + _S5392.y;
    if(_S5101)
    {
        float _S5395 = _S5394 / _S5103;
        float _S5396 = _S5104 * - _S5395;
        float _S5397 = _S5100 * (0.3333333432674408f * - (_S5099 * _S5395));
        k_9 = _S5397 + _S5397;
        _S5102 = _S5396;
        _S5103 = 0.0f;
    }
    else
    {
        float _S5398 = _S5394 / _S5102;
        float _S5399 = _S5100 * - _S5398;
        k_9 = _S5098 * _S5398;
        _S5102 = 0.0f;
        _S5103 = _S5399;
    }
    DiffPair_float_0 _S5400;
    (&_S5400)->primal_0 = _S5098;
    (&_S5400)->differential_0 = 0.0f;
    DiffPair_float_0 _S5401;
    (&_S5401)->primal_0 = _S5099;
    (&_S5401)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5400, &_S5401, k_9);
    float _S5402 = _S5401.differential_0 + _S5102;
    float _S5403 = _S5400.differential_0 + _S5103;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5404;
    (&_S5404)->primal_0 = _S5097;
    (&_S5404)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5404, _S5403);
    float2  _S5405 = _S5404.differential_0 + _S5393;
    float _S5406 = fx_30 * (_S5358.x + _S5357);
    float2  _S5407 = make_float2 (_S5406, fy_30 * (_S5358.y + _S5360)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S5406, (*dist_coeffs_36)[int(9)] * _S5406);
    float2  _S5408 = _S5081 * _S5407;
    float _S5409 = (*dist_coeffs_36)[int(4)] * _S5407.y;
    float _S5410 = (*dist_coeffs_36)[int(5)] * _S5407.x;
    float _S5411 = _S5408.x + _S5408.y;
    float _S5412 = r2_55 * _S5411;
    float _S5413 = r2_55 * _S5412;
    float _S5414 = (*dist_coeffs_36)[int(7)] * _S5407.y + _S5409 + (*dist_coeffs_36)[int(6)] * _S5407.x + _S5410 + _S5084 * _S5411 + _S5083 * _S5412 + _S5082 * _S5413 + (*dist_coeffs_36)[int(3)] * (r2_55 * _S5413);
    float _S5415 = v_55 * _S5414;
    float _S5416 = u_55 * _S5414;
    float3  _S5417 = _S5334 + make_float3 (_S5405.x, _S5405.y, _S5402);
    float2  _S5418 = _S5085 * _S5407 + make_float2 (_S5089 * (v_55 * _S5407.y) + _S5088 * _S5410 + 2.0f * (u_55 * _S5410) + _S5086 * (v_55 * _S5407.x) + _S5416 + _S5416, _S5091 * _S5409 + 2.0f * (v_55 * _S5409) + _S5090 * _S5407.y + _S5087 * _S5407.x + _S5415 + _S5415);
    float2  _S5419 = _S5069 * _S5418;
    float2  _S5420 = _S5080 * _S5418;
    float _S5421 = _S5419.x + _S5419.y;
    if(_S5073)
    {
        float _S5422 = _S5421 / _S5075;
        float _S5423 = _S5076 * - _S5422;
        float _S5424 = _S5072 * (0.3333333432674408f * - (_S5071 * _S5422));
        k_9 = _S5424 + _S5424;
        _S5074 = _S5423;
        _S5075 = 0.0f;
    }
    else
    {
        float _S5425 = _S5421 / _S5074;
        float _S5426 = _S5072 * - _S5425;
        k_9 = _S5070 * _S5425;
        _S5074 = 0.0f;
        _S5075 = _S5426;
    }
    DiffPair_float_0 _S5427;
    (&_S5427)->primal_0 = _S5070;
    (&_S5427)->differential_0 = 0.0f;
    DiffPair_float_0 _S5428;
    (&_S5428)->primal_0 = _S5071;
    (&_S5428)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5427, &_S5428, k_9);
    float _S5429 = _S5428.differential_0 + _S5074;
    float _S5430 = _S5427.differential_0 + _S5075;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5431;
    (&_S5431)->primal_0 = _S5069;
    (&_S5431)->differential_0 = _S5309;
    s_bwd_length_impl_1(&_S5431, _S5430);
    float2  _S5432 = _S5431.differential_0 + _S5420;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5433;
    (&_S5433)->primal_0 = vert2_c_11;
    (&_S5433)->differential_0 = _S5202;
    s_bwd_length_impl_0(&_S5433, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5434;
    (&_S5434)->primal_0 = vert1_c_11;
    (&_S5434)->differential_0 = _S5202;
    s_bwd_length_impl_0(&_S5434, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5435;
    (&_S5435)->primal_0 = vert0_c_11;
    (&_S5435)->differential_0 = _S5202;
    s_bwd_length_impl_0(&_S5435, 0.0f);
    float3  _S5436 = _S5433.differential_0 + _S5390;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5437;
    (&_S5437)->primal_0 = R_30;
    (&_S5437)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5438;
    (&_S5438)->primal_0 = vert2_7;
    (&_S5438)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5437, &_S5438, _S5436);
    float3  _S5439 = _S5434.differential_0 + _S5417;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5440;
    (&_S5440)->primal_0 = R_30;
    (&_S5440)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5441;
    (&_S5441)->primal_0 = vert1_7;
    (&_S5441)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5440, &_S5441, _S5439);
    float3  _S5442 = _S5435.differential_0 + _S5335 + make_float3 (_S5432.x, _S5432.y, _S5429);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5443;
    (&_S5443)->primal_0 = R_30;
    (&_S5443)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5444;
    (&_S5444)->primal_0 = vert0_7;
    (&_S5444)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5443, &_S5444, _S5442);
    float3  _S5445 = _S5438.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5446;
    (&_S5446)->primal_0 = _S5063;
    (&_S5446)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5447;
    (&_S5447)->primal_0 = _S5068;
    (&_S5447)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5446, &_S5447, _S5445);
    float _S5448 = - _S5447.differential_0.y;
    float _S5449 = _S5067 * _S5447.differential_0.x;
    float _S5450 = - (_S5059 * _S5447.differential_0.x);
    float3  _S5451 = _S5441.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5452;
    (&_S5452)->primal_0 = _S5063;
    (&_S5452)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5453;
    (&_S5453)->primal_0 = _S5066;
    (&_S5453)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5452, &_S5453, _S5451);
    float _S5454 = _S5059 * _S5453.differential_0.x;
    float _S5455 = _S5065 * _S5453.differential_0.x;
    float3  _S5456 = _S5444.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5457;
    (&_S5457)->primal_0 = _S5063;
    (&_S5457)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5458;
    (&_S5458)->primal_0 = _S5064;
    (&_S5458)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5457, &_S5458, _S5456);
    Matrix<float, 3, 3>  _S5459 = transpose_0(_S5446.differential_0 + _S5452.differential_0 + _S5457.differential_0);
    float _S5460 = 2.0f * - _S5459.rows[int(2)].z;
    float _S5461 = 2.0f * _S5459.rows[int(2)].y;
    float _S5462 = 2.0f * _S5459.rows[int(2)].x;
    float _S5463 = 2.0f * _S5459.rows[int(1)].z;
    float _S5464 = 2.0f * - _S5459.rows[int(1)].y;
    float _S5465 = 2.0f * _S5459.rows[int(1)].x;
    float _S5466 = 2.0f * _S5459.rows[int(0)].z;
    float _S5467 = 2.0f * _S5459.rows[int(0)].y;
    float _S5468 = 2.0f * - _S5459.rows[int(0)].x;
    float _S5469 = - _S5465 + _S5467;
    float _S5470 = _S5462 + - _S5466;
    float _S5471 = - _S5461 + _S5463;
    float _S5472 = _S5461 + _S5463;
    float _S5473 = _S5462 + _S5466;
    float _S5474 = _S5465 + _S5467;
    float _S5475 = quat_31.w * (_S5464 + _S5468);
    float _S5476 = quat_31.z * (_S5460 + _S5468);
    float _S5477 = quat_31.y * (_S5460 + _S5464);
    float _S5478 = quat_31.x * _S5469 + quat_31.z * _S5472 + quat_31.y * _S5473 + _S5475 + _S5475;
    float _S5479 = quat_31.x * _S5470 + quat_31.w * _S5472 + quat_31.y * _S5474 + _S5476 + _S5476;
    float _S5480 = quat_31.x * _S5471 + quat_31.w * _S5473 + quat_31.z * _S5474 + _S5477 + _S5477;
    float _S5481 = quat_31.w * _S5469 + quat_31.z * _S5470 + quat_31.y * _S5471;
    float _S5482 = _S5450 + _S5454;
    float _S5483 = 0.5f * - _S5482;
    float _S5484 = _S5448 + _S5453.differential_0.y;
    DiffPair_float_0 _S5485;
    (&_S5485)->primal_0 = _S5060;
    (&_S5485)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5485, _S5484);
    float _S5486 = _S5483 + _S5485.differential_0;
    float _S5487 = _S5449 + _S5455 + _S5458.differential_0.x;
    DiffPair_float_0 _S5488;
    (&_S5488)->primal_0 = _S5058;
    (&_S5488)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5488, _S5487);
    float _S5489 = _S5483 + _S5488.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5490;
    (&_S5490)->primal_0 = mean_c_26;
    (&_S5490)->differential_0 = _S5202;
    s_bwd_length_impl_0(&_S5490, 0.0f);
    float3  _S5491 = _S5490.differential_0 + _S5205.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5492;
    (&_S5492)->primal_0 = R_30;
    (&_S5492)->differential_0 = _S5283;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5493;
    (&_S5493)->primal_0 = mean_30;
    (&_S5493)->differential_0 = _S5202;
    s_bwd_prop_mul_1(&_S5492, &_S5493, _S5491);
    float3  _S5494 = _S5436 + _S5439 + _S5442 + _S5491 + _S5286.differential_0;
    Matrix<float, 3, 3>  _S5495 = _S5437.differential_0 + _S5440.differential_0 + _S5443.differential_0 + _S5492.differential_0 + _S5287;
    float3  _S5496 = make_float3 (_S5489, _S5486, _S5482);
    float4  _S5497 = make_float4 (0.0f);
    *&((&_S5497)->w) = _S5478;
    *&((&_S5497)->z) = _S5479;
    *&((&_S5497)->y) = _S5480;
    *&((&_S5497)->x) = _S5481;
    float4  _S5498 = _S5497;
    float3  _S5499 = _S5445 + _S5451 + _S5456 + _S5493.differential_0 + _S5281;
    *v_mean_10 = _S5499;
    *v_quat_9 = _S5498;
    *v_scale_9 = _S5496;
    *v_hardness_5 = _S5363;
    (*v_sh_coeffs_8)[int(0)] = _S5341;
    (*v_sh_coeffs_8)[int(1)] = _S5342;
    (*v_sh_coeffs_8)[int(2)] = _S5343;
    (*v_sh_coeffs_8)[int(3)] = _S5344;
    (*v_sh_coeffs_8)[int(4)] = _S5345;
    (*v_sh_coeffs_8)[int(5)] = _S5346;
    (*v_sh_coeffs_8)[int(6)] = _S5347;
    (*v_sh_coeffs_8)[int(7)] = _S5348;
    (*v_sh_coeffs_8)[int(8)] = _S5349;
    (*v_sh_coeffs_8)[int(9)] = _S5350;
    (*v_sh_coeffs_8)[int(10)] = _S5351;
    (*v_sh_coeffs_8)[int(11)] = _S5352;
    (*v_sh_coeffs_8)[int(12)] = _S5353;
    (*v_sh_coeffs_8)[int(13)] = _S5354;
    (*v_sh_coeffs_8)[int(14)] = _S5355;
    (*v_sh_coeffs_8)[int(15)] = _S5356;
    (*v_ch_coeffs_3)[int(0)] = _S5338;
    (*v_ch_coeffs_3)[int(1)] = _S5339;
    *v_R_9 = _S5495;
    *v_t_9 = _S5494;
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
    float u_58 = d_0 * dot_0(- q_2, v2v0_0);
    float v_58 = d_0 * dot_0(q_2, v1v0_0);
    float t_30 = d_0 * dot_0(- n_0, rov0_0);
    bool _S5500;
    if(u_58 >= 0.0f)
    {
        _S5500 = v_58 >= 0.0f;
    }
    else
    {
        _S5500 = false;
    }
    if(_S5500)
    {
        _S5500 = (u_58 + v_58) <= 1.0f;
    }
    else
    {
        _S5500 = false;
    }
    if(_S5500)
    {
        _S5500 = t_30 >= 0.0f;
    }
    else
    {
        _S5500 = false;
    }
    if(!_S5500)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_58), (v_58)))), ((F32_sqrt((0.5f))) * (1.0f - u_58 - v_58)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5501;
    if(opac_0 < 0.0f)
    {
        _S5501 = 0.0f;
    }
    else
    {
        _S5501 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5501;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5502 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5503 = *dpray_d_4;
    float3  v1v0_1 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_1 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_1 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5504 = s_primal_ctx_cross_0(v1v0_1, v2v0_1);
    float3  _S5505 = s_primal_ctx_cross_0(rov0_1, (*dpray_d_4).primal_0);
    float _S5506 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S5504);
    float d_1 = 1.0f / _S5506;
    float _S5507 = _S5506 * _S5506;
    float3  _S5508 = - _S5505;
    float _S5509 = s_primal_ctx_dot_0(_S5508, v2v0_1);
    float u_59 = d_1 * _S5509;
    float _S5510 = s_primal_ctx_dot_0(_S5505, v1v0_1);
    float v_59 = d_1 * _S5510;
    float3  _S5511 = - _S5504;
    float t_31 = d_1 * s_primal_ctx_dot_0(_S5511, rov0_1);
    bool _S5512;
    if(u_59 >= 0.0f)
    {
        _S5512 = v_59 >= 0.0f;
    }
    else
    {
        _S5512 = false;
    }
    if(_S5512)
    {
        _S5512 = (u_59 + v_59) <= 1.0f;
    }
    else
    {
        _S5512 = false;
    }
    if(_S5512)
    {
        _S5512 = t_31 >= 0.0f;
    }
    else
    {
        _S5512 = false;
    }
    bool _S5513 = !!_S5512;
    float _S5514;
    float _S5515;
    float _S5516;
    float _S5517;
    float _S5518;
    float _S5519;
    float _S5520;
    float _S5521;
    float _S5522;
    float _S5523;
    float _S5524;
    if(_S5513)
    {
        float _S5525 = s_primal_ctx_min_0(u_59, v_59);
        float _S5526 = s_primal_ctx_sqrt_0(0.5f);
        float _S5527 = _S5526 * (1.0f - u_59 - v_59);
        float _S5528 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5525, _S5527) * _S5528;
        float _S5529 = _S5502.primal_0.y;
        float _S5530 = 1.0f - opac_1;
        float _S5531 = 1.0f - s_primal_ctx_clamp_0(_S5529, 0.0f, 0.99989998340606689f);
        float _S5532 = 1.0f / _S5531;
        float _S5533 = _S5531 * _S5531;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5530, _S5532);
        float o_1 = _S5502.primal_0.x;
        bool _S5534 = opac_1 < 0.0f;
        if(_S5534)
        {
            _S5514 = 0.0f;
        }
        else
        {
            _S5514 = o_1 * w_1;
        }
        _S5512 = _S5534;
        _S5515 = o_1;
        _S5516 = w_1;
        _S5517 = _S5530;
        _S5518 = _S5532;
        _S5519 = _S5533;
        _S5520 = _S5529;
        _S5521 = _S5528;
        _S5522 = _S5525;
        _S5523 = _S5527;
        _S5524 = _S5526;
    }
    else
    {
        _S5512 = false;
        _S5514 = 0.0f;
        _S5515 = 0.0f;
        _S5516 = 0.0f;
        _S5517 = 0.0f;
        _S5518 = 0.0f;
        _S5519 = 0.0f;
        _S5520 = 0.0f;
        _S5521 = 0.0f;
        _S5522 = 0.0f;
        _S5523 = 0.0f;
        _S5524 = 0.0f;
    }
    float2  _S5535 = make_float2 (0.0f);
    float2  _S5536;
    if(_S5513)
    {
        if(_S5512)
        {
            _S5514 = 0.0f;
            _S5515 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5537;
            (&_S5537)->primal_0 = _S5514;
            (&_S5537)->differential_0 = 0.0f;
            DiffPair_float_0 _S5538;
            (&_S5538)->primal_0 = 0.99500000476837158f;
            (&_S5538)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5537, &_S5538, _s_dOut_11);
            float _S5539 = _S5515 * _S5537.differential_0;
            _S5514 = _S5516 * _S5537.differential_0;
            _S5515 = _S5539;
        }
        float _S5540 = - _S5515;
        DiffPair_float_0 _S5541;
        (&_S5541)->primal_0 = _S5517;
        (&_S5541)->differential_0 = 0.0f;
        DiffPair_float_0 _S5542;
        (&_S5542)->primal_0 = _S5518;
        (&_S5542)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5541, &_S5542, _S5540);
        float _S5543 = - - (_S5542.differential_0 / _S5519);
        float s_diff_opac_T_0 = - _S5541.differential_0;
        DiffPair_float_0 _S5544;
        (&_S5544)->primal_0 = _S5520;
        (&_S5544)->differential_0 = 0.0f;
        DiffPair_float_0 _S5545;
        (&_S5545)->primal_0 = 0.0f;
        (&_S5545)->differential_0 = 0.0f;
        DiffPair_float_0 _S5546;
        (&_S5546)->primal_0 = 0.99989998340606689f;
        (&_S5546)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5544, &_S5545, &_S5546, _S5543);
        float _S5547 = _S5521 * s_diff_opac_T_0;
        DiffPair_float_0 _S5548;
        (&_S5548)->primal_0 = _S5522;
        (&_S5548)->differential_0 = 0.0f;
        DiffPair_float_0 _S5549;
        (&_S5549)->primal_0 = _S5523;
        (&_S5549)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5548, &_S5549, _S5547);
        float _S5550 = - (_S5524 * _S5549.differential_0);
        DiffPair_float_0 _S5551;
        (&_S5551)->primal_0 = u_59;
        (&_S5551)->differential_0 = 0.0f;
        DiffPair_float_0 _S5552;
        (&_S5552)->primal_0 = v_59;
        (&_S5552)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5551, &_S5552, _S5548.differential_0);
        float2  _S5553 = make_float2 (_S5514, _S5544.differential_0);
        float _S5554 = _S5550 + _S5552.differential_0;
        _S5514 = _S5550 + _S5551.differential_0;
        _S5515 = _S5554;
        _S5536 = _S5553;
    }
    else
    {
        _S5514 = 0.0f;
        _S5515 = 0.0f;
        _S5536 = _S5535;
    }
    float3  _S5555 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5556;
    (&_S5556)->primal_0 = _S5511;
    (&_S5556)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5557;
    (&_S5557)->primal_0 = rov0_1;
    (&_S5557)->differential_0 = _S5555;
    s_bwd_prop_dot_0(&_S5556, &_S5557, 0.0f);
    float3  _S5558 = - _S5556.differential_0;
    float _S5559 = d_1 * _S5515;
    float _S5560 = _S5510 * _S5515;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5561;
    (&_S5561)->primal_0 = _S5505;
    (&_S5561)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5562;
    (&_S5562)->primal_0 = v1v0_1;
    (&_S5562)->differential_0 = _S5555;
    s_bwd_prop_dot_0(&_S5561, &_S5562, _S5559);
    float _S5563 = d_1 * _S5514;
    float _S5564 = _S5509 * _S5514;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5565;
    (&_S5565)->primal_0 = _S5508;
    (&_S5565)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5566;
    (&_S5566)->primal_0 = v2v0_1;
    (&_S5566)->differential_0 = _S5555;
    s_bwd_prop_dot_0(&_S5565, &_S5566, _S5563);
    float3  _S5567 = - _S5565.differential_0;
    float _S5568 = - ((_S5560 + _S5564) / _S5507);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5569;
    (&_S5569)->primal_0 = _S5503.primal_0;
    (&_S5569)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5570;
    (&_S5570)->primal_0 = _S5504;
    (&_S5570)->differential_0 = _S5555;
    s_bwd_prop_dot_0(&_S5569, &_S5570, _S5568);
    float3  _S5571 = _S5561.differential_0 + _S5567;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5572;
    (&_S5572)->primal_0 = rov0_1;
    (&_S5572)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5573;
    (&_S5573)->primal_0 = _S5503.primal_0;
    (&_S5573)->differential_0 = _S5555;
    s_bwd_prop_cross_0(&_S5572, &_S5573, _S5571);
    float3  _S5574 = _S5558 + _S5570.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5575;
    (&_S5575)->primal_0 = v1v0_1;
    (&_S5575)->differential_0 = _S5555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5576;
    (&_S5576)->primal_0 = v2v0_1;
    (&_S5576)->differential_0 = _S5555;
    s_bwd_prop_cross_0(&_S5575, &_S5576, _S5574);
    float3  _S5577 = _S5557.differential_0 + _S5572.differential_0;
    float3  _S5578 = _S5566.differential_0 + _S5576.differential_0;
    float3  _S5579 = _S5562.differential_0 + _S5575.differential_0;
    float3  _S5580 = - _S5577 + - _S5578 + - _S5579;
    float3  _S5581 = _S5569.differential_0 + _S5573.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S5581;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S5577;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5536;
    FixedArray<float3 , 3>  _S5582;
    _S5582[int(0)] = _S5555;
    _S5582[int(1)] = _S5555;
    _S5582[int(2)] = _S5555;
    _S5582[int(2)] = _S5578;
    _S5582[int(0)] = _S5580;
    _S5582[int(1)] = _S5579;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5582;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5583, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5584, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5585, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5586, float _S5587)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5583, _S5584, _S5585, _S5586, _S5587);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_5, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5588 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5589 = { _S5588, _S5588, _S5588 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_5;
    (&dp_verts_0)->differential_0 = _S5589;
    float2  _S5590 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5590;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S5588;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S5588;
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
    float u_60 = d_2 * dot_0(- q_3, v2v0_2);
    float v_60 = d_2 * dot_0(q_3, v1v0_2);
    *depth_20 = d_2 * dot_0(- n_2, rov0_2);
    *color_13 = (*rgbs_5)[int(0)] * make_float3 (1.0f - u_60 - v_60) + (*rgbs_5)[int(1)] * make_float3 (u_60) + (*rgbs_5)[int(2)] * make_float3 (v_60);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_3 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_3 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_3 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S5591 = s_primal_ctx_cross_0(v1v0_3, v2v0_3);
    float3  _S5592 = s_primal_ctx_cross_0(rov0_3, (*dpray_d_5).primal_0);
    float _S5593 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S5591);
    float d_3 = 1.0f / _S5593;
    float _S5594 = _S5593 * _S5593;
    float3  _S5595 = - _S5592;
    float _S5596 = s_primal_ctx_dot_0(_S5595, v2v0_3);
    float u_61 = d_3 * _S5596;
    float3  _S5597 = make_float3 (u_61);
    float _S5598 = s_primal_ctx_dot_0(_S5592, v1v0_3);
    float v_61 = d_3 * _S5598;
    float3  _S5599 = make_float3 (v_61);
    float3  _S5600 = - _S5591;
    float _S5601 = s_primal_ctx_dot_0(_S5600, rov0_3);
    float _S5602 = d_3 * _S5601;
    float3  _S5603 = make_float3 (1.0f - u_61 - v_61);
    DiffPair_float_0 _S5604;
    (&_S5604)->primal_0 = s_primal_ctx_max_0(_S5602, 9.999999960041972e-13f);
    (&_S5604)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5604, dpdepth_2);
    DiffPair_float_0 _S5605;
    (&_S5605)->primal_0 = _S5602;
    (&_S5605)->differential_0 = 0.0f;
    DiffPair_float_0 _S5606;
    (&_S5606)->primal_0 = 9.999999960041972e-13f;
    (&_S5606)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5605, &_S5606, _S5604.differential_0);
    float3  _S5607 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5608 = _S5599 * dpcolor_1;
    float3  _S5609 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5610 = _S5597 * dpcolor_1;
    float3  _S5611 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5612 = _S5603 * dpcolor_1;
    float _S5613 = - (_S5611.x + _S5611.y + _S5611.z);
    float _S5614 = d_3 * _S5605.differential_0;
    float _S5615 = _S5601 * _S5605.differential_0;
    float3  _S5616 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5617;
    (&_S5617)->primal_0 = _S5600;
    (&_S5617)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5618;
    (&_S5618)->primal_0 = rov0_3;
    (&_S5618)->differential_0 = _S5616;
    s_bwd_prop_dot_0(&_S5617, &_S5618, _S5614);
    float3  _S5619 = - _S5617.differential_0;
    float _S5620 = _S5613 + _S5607.x + _S5607.y + _S5607.z;
    float _S5621 = d_3 * _S5620;
    float _S5622 = _S5598 * _S5620;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5623;
    (&_S5623)->primal_0 = _S5592;
    (&_S5623)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5624;
    (&_S5624)->primal_0 = v1v0_3;
    (&_S5624)->differential_0 = _S5616;
    s_bwd_prop_dot_0(&_S5623, &_S5624, _S5621);
    float _S5625 = _S5613 + _S5609.x + _S5609.y + _S5609.z;
    float _S5626 = d_3 * _S5625;
    float _S5627 = _S5596 * _S5625;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5628;
    (&_S5628)->primal_0 = _S5595;
    (&_S5628)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5629;
    (&_S5629)->primal_0 = v2v0_3;
    (&_S5629)->differential_0 = _S5616;
    s_bwd_prop_dot_0(&_S5628, &_S5629, _S5626);
    float3  _S5630 = - _S5628.differential_0;
    float _S5631 = - ((_S5615 + _S5622 + _S5627) / _S5594);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5632;
    (&_S5632)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5632)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5633;
    (&_S5633)->primal_0 = _S5591;
    (&_S5633)->differential_0 = _S5616;
    s_bwd_prop_dot_0(&_S5632, &_S5633, _S5631);
    float3  _S5634 = _S5623.differential_0 + _S5630;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5635;
    (&_S5635)->primal_0 = rov0_3;
    (&_S5635)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5636;
    (&_S5636)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5636)->differential_0 = _S5616;
    s_bwd_prop_cross_0(&_S5635, &_S5636, _S5634);
    float3  _S5637 = _S5619 + _S5633.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5638;
    (&_S5638)->primal_0 = v1v0_3;
    (&_S5638)->differential_0 = _S5616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5639;
    (&_S5639)->primal_0 = v2v0_3;
    (&_S5639)->differential_0 = _S5616;
    s_bwd_prop_cross_0(&_S5638, &_S5639, _S5637);
    float3  _S5640 = _S5618.differential_0 + _S5635.differential_0;
    float3  _S5641 = _S5629.differential_0 + _S5639.differential_0;
    float3  _S5642 = _S5624.differential_0 + _S5638.differential_0;
    float3  _S5643 = - _S5640 + - _S5641 + - _S5642;
    float3  _S5644 = _S5632.differential_0 + _S5636.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S5644;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S5640;
    FixedArray<float3 , 3>  _S5645;
    _S5645[int(0)] = _S5616;
    _S5645[int(1)] = _S5616;
    _S5645[int(2)] = _S5616;
    _S5645[int(2)] = _S5608;
    _S5645[int(1)] = _S5610;
    _S5645[int(0)] = _S5612;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5645;
    FixedArray<float3 , 3>  _S5646;
    _S5646[int(0)] = _S5616;
    _S5646[int(1)] = _S5616;
    _S5646[int(2)] = _S5616;
    _S5646[int(2)] = _S5641;
    _S5646[int(0)] = _S5643;
    _S5646[int(1)] = _S5642;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5646;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5647, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5648, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5649, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5650, float3  _S5651, float _S5652)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5647, _S5648, _S5649, _S5650, _S5651, _S5652);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_8, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_9, float3  ray_d_9, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5653 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5654 = { _S5653, _S5653, _S5653 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_8;
    (&dp_verts_1)->differential_0 = _S5654;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5654;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_9;
    (&dp_ray_o_3)->differential_0 = _S5653;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_9;
    (&dp_ray_d_3)->differential_0 = _S5653;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

