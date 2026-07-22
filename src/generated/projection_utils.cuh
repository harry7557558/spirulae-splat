#pragma once

#include "generated/slang.cuh"

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

inline __device__ Matrix<float, 2, 2>  transpose_1(Matrix<float, 2, 2>  x_1)
{
    Matrix<float, 2, 2>  result_1;
    int r_1 = int(0);
    for(;;)
    {
        if(r_1 < int(2))
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

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_0)
{
    float x_2 = quat_0.y;
    float x2_0 = x_2 * x_2;
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

inline __device__ void mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_0, Matrix<float, 3, 3>  dOut_0)
{
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_0;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_0.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(0)].x * dOut_0.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_0.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(0)].y * dOut_0.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_0.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(0)].z * dOut_0.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_0.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(0)].x * dOut_0.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_0.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(0)].y * dOut_0.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_0.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(0)].z * dOut_0.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_0.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(0)].x * dOut_0.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_0.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(0)].y * dOut_0.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_0.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(0)].z * dOut_0.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_0.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(1)].x * dOut_0.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_0.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(1)].y * dOut_0.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_0.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(1)].z * dOut_0.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_0.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(1)].x * dOut_0.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_0.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(1)].y * dOut_0.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_0.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(1)].z * dOut_0.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_0.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(1)].x * dOut_0.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_0.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(1)].y * dOut_0.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_0.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(1)].z * dOut_0.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_0.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(2)].x * dOut_0.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_0.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(2)].y * dOut_0.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_0.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(2)].z * dOut_0.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_0.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(2)].x * dOut_0.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_0.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(2)].y * dOut_0.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_0.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(2)].z * dOut_0.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_0.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(2)].x * dOut_0.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_0.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(2)].y * dOut_0.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_0.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(2)].z * dOut_0.rows[int(2)].z;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
{
    Matrix<float, 3, 3>  result_2;
    int r_2 = int(0);
    for(;;)
    {
        if(r_2 < int(3))
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
            int i_0 = int(0);
            float sum_0 = 0.0f;
            for(;;)
            {
                if(i_0 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_1 = sum_0 + _slang_vector_get_element(left_1.rows[r_2], i_0) * _slang_vector_get_element(right_1.rows[i_0], c_2);
                i_0 = i_0 + int(1);
                sum_0 = sum_1;
            }
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_2)), c_2) = sum_0;
            c_2 = c_2 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_2;
}

inline __device__ void quat_scale_to_covar(float4  quat_1, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_3 = quat_1.y;
    float x2_1 = x_3 * x_3;
    float y2_1 = quat_1.z * quat_1.z;
    float z2_1 = quat_1.w * quat_1.w;
    float xy_1 = quat_1.y * quat_1.z;
    float xz_1 = quat_1.y * quat_1.w;
    float yz_1 = quat_1.z * quat_1.w;
    float wx_1 = quat_1.x * quat_1.y;
    float wy_1 = quat_1.x * quat_1.z;
    float wz_1 = quat_1.x * quat_1.w;
    Matrix<float, 3, 3>  M_0 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_1(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_2, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_4 = quat_2.y;
    float x2_2 = x_4 * x_4;
    float y2_2 = quat_2.z * quat_2.z;
    float z2_2 = quat_2.w * quat_2.w;
    float xy_2 = quat_2.y * quat_2.z;
    float xz_2 = quat_2.y * quat_2.w;
    float yz_2 = quat_2.z * quat_2.w;
    float wx_2 = quat_2.x * quat_2.y;
    float wy_2 = quat_2.x * quat_2.z;
    float wz_2 = quat_2.x * quat_2.w;
    *M_1 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_2, float3  dOut_1)
{
    float _S1 = (*left_2).primal_0.rows[int(0)].x * dOut_1.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_2).primal_0.x * dOut_1.x;
    float sum_2 = _S1 + (*left_2).primal_0.rows[int(1)].x * dOut_1.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_2).primal_0.x * dOut_1.y;
    float sum_3 = sum_2 + (*left_2).primal_0.rows[int(2)].x * dOut_1.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_2).primal_0.x * dOut_1.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_3;
    float _S2 = (*left_2).primal_0.rows[int(0)].y * dOut_1.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_2).primal_0.y * dOut_1.x;
    float sum_4 = _S2 + (*left_2).primal_0.rows[int(1)].y * dOut_1.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_2).primal_0.y * dOut_1.y;
    float sum_5 = sum_4 + (*left_2).primal_0.rows[int(2)].y * dOut_1.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_2).primal_0.y * dOut_1.z;
    *&((&right_d_result_1)->y) = sum_5;
    float _S3 = (*left_2).primal_0.rows[int(0)].z * dOut_1.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_2).primal_0.z * dOut_1.x;
    float sum_6 = _S3 + (*left_2).primal_0.rows[int(1)].z * dOut_1.y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = (*right_2).primal_0.z * dOut_1.y;
    float sum_7 = sum_6 + (*left_2).primal_0.rows[int(2)].z * dOut_1.z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = (*right_2).primal_0.z * dOut_1.z;
    *&((&right_d_result_1)->z) = sum_7;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_2(Matrix<float, 3, 3>  left_3, float3  right_3)
{
    float3  result_3;
    int i_1 = int(0);
    for(;;)
    {
        if(i_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_3.rows[i_1], j_0) * _slang_vector_get_element(right_3, j_0);
            j_0 = j_0 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_3, i_1) = sum_8;
        i_1 = i_1 + int(1);
    }
    return result_3;
}

inline __device__ float3  apply_sqrt_covar_to_vec(float4  quat_3, float3  scale_2, float3  vec_0)
{
    float x_5 = quat_3.y;
    float x2_3 = x_5 * x_5;
    float y2_3 = quat_3.z * quat_3.z;
    float z2_3 = quat_3.w * quat_3.w;
    float xy_3 = quat_3.y * quat_3.z;
    float xz_3 = quat_3.y * quat_3.w;
    float yz_3 = quat_3.z * quat_3.w;
    float wx_3 = quat_3.x * quat_3.y;
    float wy_3 = quat_3.x * quat_3.z;
    float wz_3 = quat_3.x * quat_3.w;
    return mul_2(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), scale_2 * vec_0);
}

inline __device__ float3  apply_covar_to_vec(float4  quat_4, float3  scale_3, float3  vec_1)
{
    float x_6 = quat_4.y;
    float x2_4 = x_6 * x_6;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_4 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  M_2 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (scale_3.x, 0.0f, 0.0f, 0.0f, scale_3.y, 0.0f, 0.0f, 0.0f, scale_3.z));
    return mul_2(mul_1(M_2, transpose_0(M_2)), vec_1);
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_0, FixedArray<float, 10>  dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float _S4 = 0.0f * v_0;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float s_diff_r2_0 = u_0 + u_0 + (_S4 + _S4);
    float _S5 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
    float _S6 = dist_coeffs_0[int(1)] + r2_0 * _S5;
    float _S7 = dist_coeffs_0[int(0)] + r2_0 * _S6;
    float radial_0 = 1.0f + r2_0 * _S7;
    float _S8 = 2.0f * dist_coeffs_0[int(4)];
    float _S9 = _S8 * u_0;
    float _S10 = 2.0f * u_0;
    float _S11 = 2.0f * dist_coeffs_0[int(5)];
    float _S12 = _S11 * u_0;
    float _S13 = 2.0f * v_0;
    float2  _S14 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S7 + (s_diff_r2_0 * _S6 + (s_diff_r2_0 * _S5 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (_S8 * v_0 + 0.0f * _S9 + (s_diff_r2_0 + (_S10 + _S10)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S11 * v_0 + 0.0f * _S12 + (s_diff_r2_0 + (_S4 + 0.0f * _S13)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
    float _S15 = 0.0f * u_0;
    float s_diff_r2_1 = _S15 + _S15 + (v_0 + v_0);
    float2  _S16 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S7 + (s_diff_r2_1 * _S6 + (s_diff_r2_1 * _S5 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (0.0f * _S8 * v_0 + _S9 + (s_diff_r2_1 + (_S15 + 0.0f * _S10)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S11 * v_0 + _S12 + (s_diff_r2_1 + (_S13 + _S13)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
    Matrix<float, 2, 2>  _S17 = transpose_1(makeMatrix<float, 2, 2> (_S14 + make_float2 (_S14.x * dist_coeffs_0[int(8)] + _S14.y * dist_coeffs_0[int(9)], 0.0f), _S16 + make_float2 (_S16.x * dist_coeffs_0[int(8)] + _S16.y * dist_coeffs_0[int(9)], 0.0f)));
    return (F32_min((determinant_0(_S17)), ((F32_min((_S17.rows[int(0)].x), (_S17.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ bool persp_proj_nav(float3  p_view_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, float2  * uv_1)
{
    bool _S18;
    for(;;)
    {
        float _S19 = p_view_0.z;
        *uv_1 = float2 {p_view_0.x, p_view_0.y} / make_float2 (_S19);
        if(_S19 < 0.0f)
        {
            _S18 = true;
        }
        else
        {
            float u_1 = (*uv_1).x;
            float v_1 = (*uv_1).y;
            float _S20 = 0.0f * v_1;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float s_diff_r2_2 = u_1 + u_1 + (_S20 + _S20);
            float _S21 = dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)];
            float _S22 = dist_coeffs_1[int(1)] + r2_1 * _S21;
            float _S23 = dist_coeffs_1[int(0)] + r2_1 * _S22;
            float radial_1 = 1.0f + r2_1 * _S23;
            float _S24 = 2.0f * dist_coeffs_1[int(4)];
            float _S25 = _S24 * u_1;
            float _S26 = 2.0f * u_1;
            float _S27 = 2.0f * dist_coeffs_1[int(5)];
            float _S28 = _S27 * u_1;
            float _S29 = 2.0f * v_1;
            float2  _S30 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S23 + (s_diff_r2_2 * _S22 + (s_diff_r2_2 * _S21 + s_diff_r2_2 * dist_coeffs_1[int(3)] * r2_1) * r2_1) * r2_1) * *uv_1 + make_float2 (_S24 * v_1 + 0.0f * _S25 + (s_diff_r2_2 + (_S26 + _S26)) * dist_coeffs_1[int(5)] + s_diff_r2_2 * dist_coeffs_1[int(6)], _S27 * v_1 + 0.0f * _S28 + (s_diff_r2_2 + (_S20 + 0.0f * _S29)) * dist_coeffs_1[int(4)] + s_diff_r2_2 * dist_coeffs_1[int(7)]);
            float _S31 = 0.0f * u_1;
            float s_diff_r2_3 = _S31 + _S31 + (v_1 + v_1);
            float2  _S32 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S23 + (s_diff_r2_3 * _S22 + (s_diff_r2_3 * _S21 + s_diff_r2_3 * dist_coeffs_1[int(3)] * r2_1) * r2_1) * r2_1) * *uv_1 + make_float2 (0.0f * _S24 * v_1 + _S25 + (s_diff_r2_3 + (_S31 + 0.0f * _S26)) * dist_coeffs_1[int(5)] + s_diff_r2_3 * dist_coeffs_1[int(6)], 0.0f * _S27 * v_1 + _S28 + (s_diff_r2_3 + (_S29 + _S29)) * dist_coeffs_1[int(4)] + s_diff_r2_3 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S33 = transpose_1(makeMatrix<float, 2, 2> (_S30 + make_float2 (_S30.x * dist_coeffs_1[int(8)] + _S30.y * dist_coeffs_1[int(9)], 0.0f), _S32 + make_float2 (_S32.x * dist_coeffs_1[int(8)] + _S32.y * dist_coeffs_1[int(9)], 0.0f)));
            _S18 = !((F32_min((determinant_0(_S33)), ((F32_min((_S33.rows[int(0)].x), (_S33.rows[int(1)].y)))))) > 0.0f);
        }
        if(_S18)
        {
            break;
        }
        float u_2 = (*uv_1).x;
        float v_2 = (*uv_1).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S34 = *uv_1 * make_float2 (1.0f + r2_2 * (dist_coeffs_1[int(0)] + r2_2 * (dist_coeffs_1[int(1)] + r2_2 * (dist_coeffs_1[int(2)] + r2_2 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_2 * v_2 + dist_coeffs_1[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + dist_coeffs_1[int(6)] * r2_2, 2.0f * dist_coeffs_1[int(5)] * u_2 * v_2 + dist_coeffs_1[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + dist_coeffs_1[int(7)] * r2_2);
        float2  _S35 = _S34 + make_float2 (dist_coeffs_1[int(8)] * _S34.x + dist_coeffs_1[int(9)] * _S34.y, 0.0f);
        *uv_1 = make_float2 (intrins_0.x * _S35.x + intrins_0.z, intrins_0.y * _S35.y + intrins_0.w);
        break;
    }
    return !_S18;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac(float3  p_view_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2)
{
    float2  _S36 = float2 {p_view_1.x, p_view_1.y};
    float _S37 = p_view_1.z;
    float2  _S38 = _S36 / make_float2 (_S37);
    float2  _S39 = _S36 * make_float2 (0.0f);
    float _S40 = _S37 * _S37;
    float2  _S41 = (make_float2 (1.0f, 0.0f) * make_float2 (_S37) - _S39) / make_float2 (_S40);
    float u_3 = _S38.x;
    float s_diff_u_0 = _S41.x;
    float v_3 = _S38.y;
    float s_diff_v_0 = _S41.y;
    float _S42 = s_diff_u_0 * u_3;
    float _S43 = s_diff_v_0 * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = _S42 + _S42 + (_S43 + _S43);
    float _S44 = dist_coeffs_2[int(2)] + r2_3 * dist_coeffs_2[int(3)];
    float _S45 = dist_coeffs_2[int(1)] + r2_3 * _S44;
    float _S46 = dist_coeffs_2[int(0)] + r2_3 * _S45;
    float _S47 = 2.0f * dist_coeffs_2[int(4)];
    float _S48 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S49 = _S41 * make_float2 (1.0f + r2_3 * _S46) + make_float2 (s_diff_r2_4 * _S46 + (s_diff_r2_4 * _S45 + (s_diff_r2_4 * _S44 + s_diff_r2_4 * dist_coeffs_2[int(3)] * r2_3) * r2_3) * r2_3) * _S38 + make_float2 (s_diff_u_0 * _S47 * v_3 + s_diff_v_0 * (_S47 * u_3) + (s_diff_r2_4 + (s_diff_u_0 * 2.0f * u_3 + s_diff_u_0 * (2.0f * u_3))) * dist_coeffs_2[int(5)] + s_diff_r2_4 * dist_coeffs_2[int(6)], s_diff_u_0 * _S48 * v_3 + s_diff_v_0 * (_S48 * u_3) + (s_diff_r2_4 + (s_diff_v_0 * 2.0f * v_3 + s_diff_v_0 * (2.0f * v_3))) * dist_coeffs_2[int(4)] + s_diff_r2_4 * dist_coeffs_2[int(7)]);
    float2  _S50 = _S49 + make_float2 (_S49.x * dist_coeffs_2[int(8)] + _S49.y * dist_coeffs_2[int(9)], 0.0f);
    float fx_0 = intrins_1.x;
    float fy_0 = intrins_1.y;
    float _S51 = _S50.y * fy_0;
    Matrix<float, 2, 3>  J_0;
    *&(((&J_0)->rows + (int(0)))->x) = _S50.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->x) = _S51;
    float2  _S52 = _S36 / make_float2 (_S37);
    float2  _S53 = (make_float2 (0.0f, 1.0f) * make_float2 (_S37) - _S39) / make_float2 (_S40);
    float u_4 = _S52.x;
    float s_diff_u_1 = _S53.x;
    float v_4 = _S52.y;
    float s_diff_v_1 = _S53.y;
    float _S54 = s_diff_u_1 * u_4;
    float _S55 = s_diff_v_1 * v_4;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float s_diff_r2_5 = _S54 + _S54 + (_S55 + _S55);
    float _S56 = dist_coeffs_2[int(2)] + r2_4 * dist_coeffs_2[int(3)];
    float _S57 = dist_coeffs_2[int(1)] + r2_4 * _S56;
    float _S58 = dist_coeffs_2[int(0)] + r2_4 * _S57;
    float2  _S59 = _S53 * make_float2 (1.0f + r2_4 * _S58) + make_float2 (s_diff_r2_5 * _S58 + (s_diff_r2_5 * _S57 + (s_diff_r2_5 * _S56 + s_diff_r2_5 * dist_coeffs_2[int(3)] * r2_4) * r2_4) * r2_4) * _S52 + make_float2 (s_diff_u_1 * _S47 * v_4 + s_diff_v_1 * (_S47 * u_4) + (s_diff_r2_5 + (s_diff_u_1 * 2.0f * u_4 + s_diff_u_1 * (2.0f * u_4))) * dist_coeffs_2[int(5)] + s_diff_r2_5 * dist_coeffs_2[int(6)], s_diff_u_1 * _S48 * v_4 + s_diff_v_1 * (_S48 * u_4) + (s_diff_r2_5 + (s_diff_v_1 * 2.0f * v_4 + s_diff_v_1 * (2.0f * v_4))) * dist_coeffs_2[int(4)] + s_diff_r2_5 * dist_coeffs_2[int(7)]);
    float2  _S60 = _S59 + make_float2 (_S59.x * dist_coeffs_2[int(8)] + _S59.y * dist_coeffs_2[int(9)], 0.0f);
    float _S61 = _S60.y * fy_0;
    *&(((&J_0)->rows + (int(0)))->y) = _S60.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->y) = _S61;
    float2  _S62 = _S36 / make_float2 (_S37);
    float2  _S63 = (make_float2 (0.0f, 0.0f) * make_float2 (_S37) - _S36) / make_float2 (_S40);
    float u_5 = _S62.x;
    float s_diff_u_2 = _S63.x;
    float v_5 = _S62.y;
    float s_diff_v_2 = _S63.y;
    float _S64 = s_diff_u_2 * u_5;
    float _S65 = s_diff_v_2 * v_5;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float s_diff_r2_6 = _S64 + _S64 + (_S65 + _S65);
    float _S66 = dist_coeffs_2[int(2)] + r2_5 * dist_coeffs_2[int(3)];
    float _S67 = dist_coeffs_2[int(1)] + r2_5 * _S66;
    float _S68 = dist_coeffs_2[int(0)] + r2_5 * _S67;
    float2  _S69 = _S63 * make_float2 (1.0f + r2_5 * _S68) + make_float2 (s_diff_r2_6 * _S68 + (s_diff_r2_6 * _S67 + (s_diff_r2_6 * _S66 + s_diff_r2_6 * dist_coeffs_2[int(3)] * r2_5) * r2_5) * r2_5) * _S62 + make_float2 (s_diff_u_2 * _S47 * v_5 + s_diff_v_2 * (_S47 * u_5) + (s_diff_r2_6 + (s_diff_u_2 * 2.0f * u_5 + s_diff_u_2 * (2.0f * u_5))) * dist_coeffs_2[int(5)] + s_diff_r2_6 * dist_coeffs_2[int(6)], s_diff_u_2 * _S48 * v_5 + s_diff_v_2 * (_S48 * u_5) + (s_diff_r2_6 + (s_diff_v_2 * 2.0f * v_5 + s_diff_v_2 * (2.0f * v_5))) * dist_coeffs_2[int(4)] + s_diff_r2_6 * dist_coeffs_2[int(7)]);
    float2  _S70 = _S69 + make_float2 (_S69.x * dist_coeffs_2[int(8)] + _S69.y * dist_coeffs_2[int(9)], 0.0f);
    float _S71 = _S70.y * fy_0;
    *&(((&J_0)->rows + (int(0)))->z) = _S70.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->z) = _S71;
    return J_0;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ DiffPair_float_0 _d_sqrt_0(DiffPair_float_0 * dpx_0)
{
    DiffPair_float_0 _S72 = { (F32_sqrt((dpx_0->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_0->primal_0)))))) * dpx_0->differential_0 };
    return _S72;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, float dOut_2)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_0).primal_0.x * dOut_2;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_1).primal_0.x * dOut_2;
    *&((&x_d_result_0)->y) = (*dpy_0).primal_0.y * dOut_2;
    *&((&y_d_result_0)->y) = (*dpx_1).primal_0.y * dOut_2;
    *&((&x_d_result_0)->z) = (*dpy_0).primal_0.z * dOut_2;
    *&((&y_d_result_0)->z) = (*dpx_1).primal_0.z * dOut_2;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = x_d_result_0;
    dpy_0->primal_0 = (*dpy_0).primal_0;
    dpy_0->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_7, float3  y_0)
{
    int i_2 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_2 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_5 = result_4 + _slang_vector_get_element(x_7, i_2) * _slang_vector_get_element(y_0, i_2);
        i_2 = i_2 + int(1);
        result_4 = result_5;
    }
    return result_4;
}

inline __device__ float dot_1(float2  x_8, float2  y_1)
{
    int i_3 = int(0);
    float result_6 = 0.0f;
    for(;;)
    {
        if(i_3 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_7 = result_6 + _slang_vector_get_element(x_8, i_3) * _slang_vector_get_element(y_1, i_3);
        i_3 = i_3 + int(1);
        result_6 = result_7;
    }
    return result_6;
}

inline __device__ float length_0(float2  x_9)
{
    return (F32_sqrt((dot_1(x_9, x_9))));
}

inline __device__ float length_1(float3  x_10)
{
    return (F32_sqrt((dot_0(x_10, x_10))));
}

inline __device__ DiffPair_float_0 _d_atan2_0(DiffPair_float_0 * dpy_1, DiffPair_float_0 * dpx_2)
{
    float _S73 = dpx_2->primal_0 * dpx_2->primal_0 + dpy_1->primal_0 * dpy_1->primal_0;
    DiffPair_float_0 _S74 = { (F32_atan2((dpy_1->primal_0), (dpx_2->primal_0))), - dpy_1->primal_0 / _S73 * dpx_2->differential_0 + dpx_2->primal_0 / _S73 * dpy_1->differential_0 };
    return _S74;
}

inline __device__ bool fisheye_proj_nav(float3  p_view_2, float4  intrins_2, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_2)
{
    bool _S75;
    for(;;)
    {
        float2  _S76 = float2 {p_view_2.x, p_view_2.y};
        float r_3 = length_0(_S76);
        float _S77 = p_view_2.z;
        float theta_0 = (F32_atan2((r_3), (_S77)));
        float k_0;
        if(theta_0 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S77;
        }
        else
        {
            k_0 = theta_0 / r_3;
        }
        float2  _S78 = _S76 * make_float2 (k_0);
        *uv_2 = _S78;
        float u_6 = _S78.x;
        float v_6 = _S78.y;
        float _S79 = 0.0f * v_6;
        float r2_6 = u_6 * u_6 + v_6 * v_6;
        float s_diff_r2_7 = u_6 + u_6 + (_S79 + _S79);
        float _S80 = dist_coeffs_3[int(2)] + r2_6 * dist_coeffs_3[int(3)];
        float _S81 = dist_coeffs_3[int(1)] + r2_6 * _S80;
        float _S82 = dist_coeffs_3[int(0)] + r2_6 * _S81;
        float radial_2 = 1.0f + r2_6 * _S82;
        float _S83 = 2.0f * dist_coeffs_3[int(4)];
        float _S84 = _S83 * u_6;
        float _S85 = 2.0f * u_6;
        float _S86 = 2.0f * dist_coeffs_3[int(5)];
        float _S87 = _S86 * u_6;
        float _S88 = 2.0f * v_6;
        float2  _S89 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_7 * _S82 + (s_diff_r2_7 * _S81 + (s_diff_r2_7 * _S80 + s_diff_r2_7 * dist_coeffs_3[int(3)] * r2_6) * r2_6) * r2_6) * _S78 + make_float2 (_S83 * v_6 + 0.0f * _S84 + (s_diff_r2_7 + (_S85 + _S85)) * dist_coeffs_3[int(5)] + s_diff_r2_7 * dist_coeffs_3[int(6)], _S86 * v_6 + 0.0f * _S87 + (s_diff_r2_7 + (_S79 + 0.0f * _S88)) * dist_coeffs_3[int(4)] + s_diff_r2_7 * dist_coeffs_3[int(7)]);
        float _S90 = 0.0f * u_6;
        float s_diff_r2_8 = _S90 + _S90 + (v_6 + v_6);
        float2  _S91 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_8 * _S82 + (s_diff_r2_8 * _S81 + (s_diff_r2_8 * _S80 + s_diff_r2_8 * dist_coeffs_3[int(3)] * r2_6) * r2_6) * r2_6) * _S78 + make_float2 (0.0f * _S83 * v_6 + _S84 + (s_diff_r2_8 + (_S90 + 0.0f * _S85)) * dist_coeffs_3[int(5)] + s_diff_r2_8 * dist_coeffs_3[int(6)], 0.0f * _S86 * v_6 + _S87 + (s_diff_r2_8 + (_S88 + _S88)) * dist_coeffs_3[int(4)] + s_diff_r2_8 * dist_coeffs_3[int(7)]);
        Matrix<float, 2, 2>  _S92 = transpose_1(makeMatrix<float, 2, 2> (_S89 + make_float2 (_S89.x * dist_coeffs_3[int(8)] + _S89.y * dist_coeffs_3[int(9)], 0.0f), _S91 + make_float2 (_S91.x * dist_coeffs_3[int(8)] + _S91.y * dist_coeffs_3[int(9)], 0.0f)));
        bool _S93 = !((F32_min((determinant_0(_S92)), ((F32_min((_S92.rows[int(0)].x), (_S92.rows[int(1)].y)))))) > 0.0f);
        _S75 = _S93;
        if(_S93)
        {
            break;
        }
        float u_7 = (*uv_2).x;
        float v_7 = (*uv_2).y;
        float r2_7 = u_7 * u_7 + v_7 * v_7;
        float2  _S94 = *uv_2 * make_float2 (1.0f + r2_7 * (dist_coeffs_3[int(0)] + r2_7 * (dist_coeffs_3[int(1)] + r2_7 * (dist_coeffs_3[int(2)] + r2_7 * dist_coeffs_3[int(3)])))) + make_float2 (_S83 * u_7 * v_7 + dist_coeffs_3[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + dist_coeffs_3[int(6)] * r2_7, _S86 * u_7 * v_7 + dist_coeffs_3[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + dist_coeffs_3[int(7)] * r2_7);
        float2  _S95 = _S94 + make_float2 (dist_coeffs_3[int(8)] * _S94.x + dist_coeffs_3[int(9)] * _S94.y, 0.0f);
        *uv_2 = make_float2 (intrins_2.x * _S95.x + intrins_2.z, intrins_2.y * _S95.y + intrins_2.w);
        break;
    }
    return !_S75;
}

inline __device__ DiffPair_float_0 _d_sin_0(DiffPair_float_0 * dpx_3)
{
    DiffPair_float_0 _S96 = { (F32_sin((dpx_3->primal_0))), (F32_cos((dpx_3->primal_0))) * dpx_3->differential_0 };
    return _S96;
}

inline __device__ bool equisolid_proj_nav(float3  p_view_3, float4  intrins_3, FixedArray<float, 10>  dist_coeffs_4, float2  * uv_3)
{
    bool _S97;
    for(;;)
    {
        float2  _S98 = float2 {p_view_3.x, p_view_3.y};
        float r_4 = length_0(_S98);
        float _S99 = p_view_3.z;
        float theta_1 = (F32_atan2((r_4), (_S99)));
        float k_1;
        if(r_4 < 9.99999997475242708e-07f)
        {
            k_1 = (1.0f - theta_1 * theta_1 / 24.0f) / _S99;
        }
        else
        {
            k_1 = 2.0f * (F32_sin((0.5f * theta_1))) / r_4;
        }
        float2  _S100 = _S98 * make_float2 (k_1);
        *uv_3 = _S100;
        float u_8 = _S100.x;
        float v_8 = _S100.y;
        float _S101 = 0.0f * v_8;
        float r2_8 = u_8 * u_8 + v_8 * v_8;
        float s_diff_r2_9 = u_8 + u_8 + (_S101 + _S101);
        float _S102 = dist_coeffs_4[int(2)] + r2_8 * dist_coeffs_4[int(3)];
        float _S103 = dist_coeffs_4[int(1)] + r2_8 * _S102;
        float _S104 = dist_coeffs_4[int(0)] + r2_8 * _S103;
        float radial_3 = 1.0f + r2_8 * _S104;
        float _S105 = 2.0f * dist_coeffs_4[int(4)];
        float _S106 = _S105 * u_8;
        float _S107 = 2.0f * u_8;
        float _S108 = 2.0f * dist_coeffs_4[int(5)];
        float _S109 = _S108 * u_8;
        float _S110 = 2.0f * v_8;
        float2  _S111 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_9 * _S104 + (s_diff_r2_9 * _S103 + (s_diff_r2_9 * _S102 + s_diff_r2_9 * dist_coeffs_4[int(3)] * r2_8) * r2_8) * r2_8) * _S100 + make_float2 (_S105 * v_8 + 0.0f * _S106 + (s_diff_r2_9 + (_S107 + _S107)) * dist_coeffs_4[int(5)] + s_diff_r2_9 * dist_coeffs_4[int(6)], _S108 * v_8 + 0.0f * _S109 + (s_diff_r2_9 + (_S101 + 0.0f * _S110)) * dist_coeffs_4[int(4)] + s_diff_r2_9 * dist_coeffs_4[int(7)]);
        float _S112 = 0.0f * u_8;
        float s_diff_r2_10 = _S112 + _S112 + (v_8 + v_8);
        float2  _S113 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_10 * _S104 + (s_diff_r2_10 * _S103 + (s_diff_r2_10 * _S102 + s_diff_r2_10 * dist_coeffs_4[int(3)] * r2_8) * r2_8) * r2_8) * _S100 + make_float2 (0.0f * _S105 * v_8 + _S106 + (s_diff_r2_10 + (_S112 + 0.0f * _S107)) * dist_coeffs_4[int(5)] + s_diff_r2_10 * dist_coeffs_4[int(6)], 0.0f * _S108 * v_8 + _S109 + (s_diff_r2_10 + (_S110 + _S110)) * dist_coeffs_4[int(4)] + s_diff_r2_10 * dist_coeffs_4[int(7)]);
        Matrix<float, 2, 2>  _S114 = transpose_1(makeMatrix<float, 2, 2> (_S111 + make_float2 (_S111.x * dist_coeffs_4[int(8)] + _S111.y * dist_coeffs_4[int(9)], 0.0f), _S113 + make_float2 (_S113.x * dist_coeffs_4[int(8)] + _S113.y * dist_coeffs_4[int(9)], 0.0f)));
        bool _S115 = !((F32_min((determinant_0(_S114)), ((F32_min((_S114.rows[int(0)].x), (_S114.rows[int(1)].y)))))) > 0.0f);
        _S97 = _S115;
        if(_S115)
        {
            break;
        }
        float u_9 = (*uv_3).x;
        float v_9 = (*uv_3).y;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float2  _S116 = *uv_3 * make_float2 (1.0f + r2_9 * (dist_coeffs_4[int(0)] + r2_9 * (dist_coeffs_4[int(1)] + r2_9 * (dist_coeffs_4[int(2)] + r2_9 * dist_coeffs_4[int(3)])))) + make_float2 (_S105 * u_9 * v_9 + dist_coeffs_4[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_4[int(6)] * r2_9, _S108 * u_9 * v_9 + dist_coeffs_4[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_4[int(7)] * r2_9);
        float2  _S117 = _S116 + make_float2 (dist_coeffs_4[int(8)] * _S116.x + dist_coeffs_4[int(9)] * _S116.y, 0.0f);
        *uv_3 = make_float2 (intrins_3.x * _S117.x + intrins_3.z, intrins_3.y * _S117.y + intrins_3.w);
        break;
    }
    return !_S97;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_4)
{
    float _S118 = *&((&dpx_4->differential_0)->x) * *&((&dpx_4->primal_0)->x);
    float _S119 = *&((&dpx_4->differential_0)->y) * *&((&dpx_4->primal_0)->y);
    float s_diff_len_0 = _S118 + _S118 + (_S119 + _S119);
    DiffPair_float_0 _S120;
    (&_S120)->primal_0 = *&((&dpx_4->primal_0)->x) * *&((&dpx_4->primal_0)->x) + *&((&dpx_4->primal_0)->y) * *&((&dpx_4->primal_0)->y);
    (&_S120)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S121 = _d_sqrt_0(&_S120);
    DiffPair_float_0 _S122 = { _S121.primal_0, _S121.differential_0 };
    return _S122;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac(float3  p_view_4, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5)
{
    Matrix<float, 2, 3>  J_1;
    float2  _S123 = float2 {p_view_4.x, p_view_4.y};
    float2  _S124 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S125;
    (&_S125)->primal_0 = _S123;
    (&_S125)->differential_0 = _S124;
    DiffPair_float_0 _S126 = s_fwd_length_impl_0(&_S125);
    float _S127 = p_view_4.z;
    DiffPair_float_0 _S128;
    (&_S128)->primal_0 = _S126.primal_0;
    (&_S128)->differential_0 = _S126.differential_0;
    DiffPair_float_0 _S129;
    (&_S129)->primal_0 = _S127;
    (&_S129)->differential_0 = 0.0f;
    DiffPair_float_0 _S130 = _d_atan2_0(&_S128, &_S129);
    float k_2;
    float s_diff_k_0;
    if((_S130.primal_0) < 0.00100000004749745f)
    {
        float _S131 = _S130.differential_0 * _S130.primal_0;
        float _S132 = 1.0f - _S130.primal_0 * _S130.primal_0 / 3.0f;
        float _S133 = ((0.0f - (_S131 + _S131) * 0.3333333432674408f) * _S127 - _S132 * 0.0f) / (_S127 * _S127);
        k_2 = _S132 / _S127;
        s_diff_k_0 = _S133;
    }
    else
    {
        float _S134 = (_S130.differential_0 * _S126.primal_0 - _S130.primal_0 * _S126.differential_0) / (_S126.primal_0 * _S126.primal_0);
        k_2 = _S130.primal_0 / _S126.primal_0;
        s_diff_k_0 = _S134;
    }
    float2  _S135 = _S123 * make_float2 (k_2);
    float2  _S136 = _S124 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S123;
    float u_10 = _S135.x;
    float s_diff_u_3 = _S136.x;
    float v_10 = _S135.y;
    float s_diff_v_3 = _S136.y;
    float _S137 = s_diff_u_3 * u_10;
    float _S138 = s_diff_v_3 * v_10;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float s_diff_r2_11 = _S137 + _S137 + (_S138 + _S138);
    float _S139 = dist_coeffs_5[int(2)] + r2_10 * dist_coeffs_5[int(3)];
    float _S140 = dist_coeffs_5[int(1)] + r2_10 * _S139;
    float _S141 = dist_coeffs_5[int(0)] + r2_10 * _S140;
    float _S142 = 2.0f * dist_coeffs_5[int(4)];
    float _S143 = 2.0f * dist_coeffs_5[int(5)];
    float2  _S144 = _S136 * make_float2 (1.0f + r2_10 * _S141) + make_float2 (s_diff_r2_11 * _S141 + (s_diff_r2_11 * _S140 + (s_diff_r2_11 * _S139 + s_diff_r2_11 * dist_coeffs_5[int(3)] * r2_10) * r2_10) * r2_10) * _S135 + make_float2 (s_diff_u_3 * _S142 * v_10 + s_diff_v_3 * (_S142 * u_10) + (s_diff_r2_11 + (s_diff_u_3 * 2.0f * u_10 + s_diff_u_3 * (2.0f * u_10))) * dist_coeffs_5[int(5)] + s_diff_r2_11 * dist_coeffs_5[int(6)], s_diff_u_3 * _S143 * v_10 + s_diff_v_3 * (_S143 * u_10) + (s_diff_r2_11 + (s_diff_v_3 * 2.0f * v_10 + s_diff_v_3 * (2.0f * v_10))) * dist_coeffs_5[int(4)] + s_diff_r2_11 * dist_coeffs_5[int(7)]);
    float2  _S145 = _S144 + make_float2 (_S144.x * dist_coeffs_5[int(8)] + _S144.y * dist_coeffs_5[int(9)], 0.0f);
    float fx_1 = intrins_4.x;
    float fy_1 = intrins_4.y;
    float _S146 = _S145.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->x) = _S145.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->x) = _S146;
    float2  _S147 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S148;
    (&_S148)->primal_0 = _S123;
    (&_S148)->differential_0 = _S147;
    DiffPair_float_0 _S149 = s_fwd_length_impl_0(&_S148);
    DiffPair_float_0 _S150;
    (&_S150)->primal_0 = _S149.primal_0;
    (&_S150)->differential_0 = _S149.differential_0;
    DiffPair_float_0 _S151;
    (&_S151)->primal_0 = _S127;
    (&_S151)->differential_0 = 0.0f;
    DiffPair_float_0 _S152 = _d_atan2_0(&_S150, &_S151);
    if((_S152.primal_0) < 0.00100000004749745f)
    {
        float _S153 = _S152.differential_0 * _S152.primal_0;
        float _S154 = 1.0f - _S152.primal_0 * _S152.primal_0 / 3.0f;
        float _S155 = ((0.0f - (_S153 + _S153) * 0.3333333432674408f) * _S127 - _S154 * 0.0f) / (_S127 * _S127);
        k_2 = _S154 / _S127;
        s_diff_k_0 = _S155;
    }
    else
    {
        float _S156 = (_S152.differential_0 * _S149.primal_0 - _S152.primal_0 * _S149.differential_0) / (_S149.primal_0 * _S149.primal_0);
        k_2 = _S152.primal_0 / _S149.primal_0;
        s_diff_k_0 = _S156;
    }
    float2  _S157 = _S123 * make_float2 (k_2);
    float2  _S158 = _S147 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S123;
    float u_11 = _S157.x;
    float s_diff_u_4 = _S158.x;
    float v_11 = _S157.y;
    float s_diff_v_4 = _S158.y;
    float _S159 = s_diff_u_4 * u_11;
    float _S160 = s_diff_v_4 * v_11;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float s_diff_r2_12 = _S159 + _S159 + (_S160 + _S160);
    float _S161 = dist_coeffs_5[int(2)] + r2_11 * dist_coeffs_5[int(3)];
    float _S162 = dist_coeffs_5[int(1)] + r2_11 * _S161;
    float _S163 = dist_coeffs_5[int(0)] + r2_11 * _S162;
    float2  _S164 = _S158 * make_float2 (1.0f + r2_11 * _S163) + make_float2 (s_diff_r2_12 * _S163 + (s_diff_r2_12 * _S162 + (s_diff_r2_12 * _S161 + s_diff_r2_12 * dist_coeffs_5[int(3)] * r2_11) * r2_11) * r2_11) * _S157 + make_float2 (s_diff_u_4 * _S142 * v_11 + s_diff_v_4 * (_S142 * u_11) + (s_diff_r2_12 + (s_diff_u_4 * 2.0f * u_11 + s_diff_u_4 * (2.0f * u_11))) * dist_coeffs_5[int(5)] + s_diff_r2_12 * dist_coeffs_5[int(6)], s_diff_u_4 * _S143 * v_11 + s_diff_v_4 * (_S143 * u_11) + (s_diff_r2_12 + (s_diff_v_4 * 2.0f * v_11 + s_diff_v_4 * (2.0f * v_11))) * dist_coeffs_5[int(4)] + s_diff_r2_12 * dist_coeffs_5[int(7)]);
    float2  _S165 = _S164 + make_float2 (_S164.x * dist_coeffs_5[int(8)] + _S164.y * dist_coeffs_5[int(9)], 0.0f);
    float _S166 = _S165.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->y) = _S165.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->y) = _S166;
    float2  _S167 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S168;
    (&_S168)->primal_0 = _S123;
    (&_S168)->differential_0 = _S167;
    DiffPair_float_0 _S169 = s_fwd_length_impl_0(&_S168);
    DiffPair_float_0 _S170;
    (&_S170)->primal_0 = _S169.primal_0;
    (&_S170)->differential_0 = _S169.differential_0;
    DiffPair_float_0 _S171;
    (&_S171)->primal_0 = _S127;
    (&_S171)->differential_0 = 1.0f;
    DiffPair_float_0 _S172 = _d_atan2_0(&_S170, &_S171);
    if((_S172.primal_0) < 0.00100000004749745f)
    {
        float _S173 = _S172.differential_0 * _S172.primal_0;
        float _S174 = 1.0f - _S172.primal_0 * _S172.primal_0 / 3.0f;
        float _S175 = ((0.0f - (_S173 + _S173) * 0.3333333432674408f) * _S127 - _S174) / (_S127 * _S127);
        k_2 = _S174 / _S127;
        s_diff_k_0 = _S175;
    }
    else
    {
        float _S176 = (_S172.differential_0 * _S169.primal_0 - _S172.primal_0 * _S169.differential_0) / (_S169.primal_0 * _S169.primal_0);
        k_2 = _S172.primal_0 / _S169.primal_0;
        s_diff_k_0 = _S176;
    }
    float2  _S177 = _S123 * make_float2 (k_2);
    float2  _S178 = _S167 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S123;
    float u_12 = _S177.x;
    float s_diff_u_5 = _S178.x;
    float v_12 = _S177.y;
    float s_diff_v_5 = _S178.y;
    float _S179 = s_diff_u_5 * u_12;
    float _S180 = s_diff_v_5 * v_12;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float s_diff_r2_13 = _S179 + _S179 + (_S180 + _S180);
    float _S181 = dist_coeffs_5[int(2)] + r2_12 * dist_coeffs_5[int(3)];
    float _S182 = dist_coeffs_5[int(1)] + r2_12 * _S181;
    float _S183 = dist_coeffs_5[int(0)] + r2_12 * _S182;
    float2  _S184 = _S178 * make_float2 (1.0f + r2_12 * _S183) + make_float2 (s_diff_r2_13 * _S183 + (s_diff_r2_13 * _S182 + (s_diff_r2_13 * _S181 + s_diff_r2_13 * dist_coeffs_5[int(3)] * r2_12) * r2_12) * r2_12) * _S177 + make_float2 (s_diff_u_5 * _S142 * v_12 + s_diff_v_5 * (_S142 * u_12) + (s_diff_r2_13 + (s_diff_u_5 * 2.0f * u_12 + s_diff_u_5 * (2.0f * u_12))) * dist_coeffs_5[int(5)] + s_diff_r2_13 * dist_coeffs_5[int(6)], s_diff_u_5 * _S143 * v_12 + s_diff_v_5 * (_S143 * u_12) + (s_diff_r2_13 + (s_diff_v_5 * 2.0f * v_12 + s_diff_v_5 * (2.0f * v_12))) * dist_coeffs_5[int(4)] + s_diff_r2_13 * dist_coeffs_5[int(7)]);
    float2  _S185 = _S184 + make_float2 (_S184.x * dist_coeffs_5[int(8)] + _S184.y * dist_coeffs_5[int(9)], 0.0f);
    float _S186 = _S185.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->z) = _S185.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->z) = _S186;
    return J_1;
}

inline __device__ Matrix<float, 2, 3>  equisolid_proj_jac(float3  p_view_5, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6)
{
    Matrix<float, 2, 3>  J_2;
    float2  _S187 = float2 {p_view_5.x, p_view_5.y};
    float2  _S188 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S189;
    (&_S189)->primal_0 = _S187;
    (&_S189)->differential_0 = _S188;
    DiffPair_float_0 _S190 = s_fwd_length_impl_0(&_S189);
    float _S191 = p_view_5.z;
    DiffPair_float_0 _S192;
    (&_S192)->primal_0 = _S190.primal_0;
    (&_S192)->differential_0 = _S190.differential_0;
    DiffPair_float_0 _S193;
    (&_S193)->primal_0 = _S191;
    (&_S193)->differential_0 = 0.0f;
    DiffPair_float_0 _S194 = _d_atan2_0(&_S192, &_S193);
    float k_3;
    float s_diff_k_1;
    if((_S190.primal_0) < 9.99999997475242708e-07f)
    {
        float _S195 = _S194.differential_0 * _S194.primal_0;
        float _S196 = 1.0f - _S194.primal_0 * _S194.primal_0 / 24.0f;
        float _S197 = ((0.0f - (_S195 + _S195) * 0.0416666679084301f) * _S191 - _S196 * 0.0f) / (_S191 * _S191);
        k_3 = _S196 / _S191;
        s_diff_k_1 = _S197;
    }
    else
    {
        float _S198 = _S194.differential_0 * 0.5f;
        DiffPair_float_0 _S199;
        (&_S199)->primal_0 = 0.5f * _S194.primal_0;
        (&_S199)->differential_0 = _S198;
        DiffPair_float_0 _S200 = _d_sin_0(&_S199);
        float _S201 = 2.0f * _S200.primal_0;
        float _S202 = (_S200.differential_0 * 2.0f * _S190.primal_0 - _S201 * _S190.differential_0) / (_S190.primal_0 * _S190.primal_0);
        k_3 = _S201 / _S190.primal_0;
        s_diff_k_1 = _S202;
    }
    float2  _S203 = _S187 * make_float2 (k_3);
    float2  _S204 = _S188 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S187;
    float u_13 = _S203.x;
    float s_diff_u_6 = _S204.x;
    float v_13 = _S203.y;
    float s_diff_v_6 = _S204.y;
    float _S205 = s_diff_u_6 * u_13;
    float _S206 = s_diff_v_6 * v_13;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float s_diff_r2_14 = _S205 + _S205 + (_S206 + _S206);
    float _S207 = dist_coeffs_6[int(2)] + r2_13 * dist_coeffs_6[int(3)];
    float _S208 = dist_coeffs_6[int(1)] + r2_13 * _S207;
    float _S209 = dist_coeffs_6[int(0)] + r2_13 * _S208;
    float _S210 = 2.0f * dist_coeffs_6[int(4)];
    float _S211 = 2.0f * dist_coeffs_6[int(5)];
    float2  _S212 = _S204 * make_float2 (1.0f + r2_13 * _S209) + make_float2 (s_diff_r2_14 * _S209 + (s_diff_r2_14 * _S208 + (s_diff_r2_14 * _S207 + s_diff_r2_14 * dist_coeffs_6[int(3)] * r2_13) * r2_13) * r2_13) * _S203 + make_float2 (s_diff_u_6 * _S210 * v_13 + s_diff_v_6 * (_S210 * u_13) + (s_diff_r2_14 + (s_diff_u_6 * 2.0f * u_13 + s_diff_u_6 * (2.0f * u_13))) * dist_coeffs_6[int(5)] + s_diff_r2_14 * dist_coeffs_6[int(6)], s_diff_u_6 * _S211 * v_13 + s_diff_v_6 * (_S211 * u_13) + (s_diff_r2_14 + (s_diff_v_6 * 2.0f * v_13 + s_diff_v_6 * (2.0f * v_13))) * dist_coeffs_6[int(4)] + s_diff_r2_14 * dist_coeffs_6[int(7)]);
    float2  _S213 = _S212 + make_float2 (_S212.x * dist_coeffs_6[int(8)] + _S212.y * dist_coeffs_6[int(9)], 0.0f);
    float fx_2 = intrins_5.x;
    float fy_2 = intrins_5.y;
    float _S214 = _S213.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->x) = _S213.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->x) = _S214;
    float2  _S215 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S216;
    (&_S216)->primal_0 = _S187;
    (&_S216)->differential_0 = _S215;
    DiffPair_float_0 _S217 = s_fwd_length_impl_0(&_S216);
    DiffPair_float_0 _S218;
    (&_S218)->primal_0 = _S217.primal_0;
    (&_S218)->differential_0 = _S217.differential_0;
    DiffPair_float_0 _S219;
    (&_S219)->primal_0 = _S191;
    (&_S219)->differential_0 = 0.0f;
    DiffPair_float_0 _S220 = _d_atan2_0(&_S218, &_S219);
    if((_S217.primal_0) < 9.99999997475242708e-07f)
    {
        float _S221 = _S220.differential_0 * _S220.primal_0;
        float _S222 = 1.0f - _S220.primal_0 * _S220.primal_0 / 24.0f;
        float _S223 = ((0.0f - (_S221 + _S221) * 0.0416666679084301f) * _S191 - _S222 * 0.0f) / (_S191 * _S191);
        k_3 = _S222 / _S191;
        s_diff_k_1 = _S223;
    }
    else
    {
        float _S224 = _S220.differential_0 * 0.5f;
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = 0.5f * _S220.primal_0;
        (&_S225)->differential_0 = _S224;
        DiffPair_float_0 _S226 = _d_sin_0(&_S225);
        float _S227 = 2.0f * _S226.primal_0;
        float _S228 = (_S226.differential_0 * 2.0f * _S217.primal_0 - _S227 * _S217.differential_0) / (_S217.primal_0 * _S217.primal_0);
        k_3 = _S227 / _S217.primal_0;
        s_diff_k_1 = _S228;
    }
    float2  _S229 = _S187 * make_float2 (k_3);
    float2  _S230 = _S215 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S187;
    float u_14 = _S229.x;
    float s_diff_u_7 = _S230.x;
    float v_14 = _S229.y;
    float s_diff_v_7 = _S230.y;
    float _S231 = s_diff_u_7 * u_14;
    float _S232 = s_diff_v_7 * v_14;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float s_diff_r2_15 = _S231 + _S231 + (_S232 + _S232);
    float _S233 = dist_coeffs_6[int(2)] + r2_14 * dist_coeffs_6[int(3)];
    float _S234 = dist_coeffs_6[int(1)] + r2_14 * _S233;
    float _S235 = dist_coeffs_6[int(0)] + r2_14 * _S234;
    float2  _S236 = _S230 * make_float2 (1.0f + r2_14 * _S235) + make_float2 (s_diff_r2_15 * _S235 + (s_diff_r2_15 * _S234 + (s_diff_r2_15 * _S233 + s_diff_r2_15 * dist_coeffs_6[int(3)] * r2_14) * r2_14) * r2_14) * _S229 + make_float2 (s_diff_u_7 * _S210 * v_14 + s_diff_v_7 * (_S210 * u_14) + (s_diff_r2_15 + (s_diff_u_7 * 2.0f * u_14 + s_diff_u_7 * (2.0f * u_14))) * dist_coeffs_6[int(5)] + s_diff_r2_15 * dist_coeffs_6[int(6)], s_diff_u_7 * _S211 * v_14 + s_diff_v_7 * (_S211 * u_14) + (s_diff_r2_15 + (s_diff_v_7 * 2.0f * v_14 + s_diff_v_7 * (2.0f * v_14))) * dist_coeffs_6[int(4)] + s_diff_r2_15 * dist_coeffs_6[int(7)]);
    float2  _S237 = _S236 + make_float2 (_S236.x * dist_coeffs_6[int(8)] + _S236.y * dist_coeffs_6[int(9)], 0.0f);
    float _S238 = _S237.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->y) = _S237.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->y) = _S238;
    float2  _S239 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S240;
    (&_S240)->primal_0 = _S187;
    (&_S240)->differential_0 = _S239;
    DiffPair_float_0 _S241 = s_fwd_length_impl_0(&_S240);
    DiffPair_float_0 _S242;
    (&_S242)->primal_0 = _S241.primal_0;
    (&_S242)->differential_0 = _S241.differential_0;
    DiffPair_float_0 _S243;
    (&_S243)->primal_0 = _S191;
    (&_S243)->differential_0 = 1.0f;
    DiffPair_float_0 _S244 = _d_atan2_0(&_S242, &_S243);
    if((_S241.primal_0) < 9.99999997475242708e-07f)
    {
        float _S245 = _S244.differential_0 * _S244.primal_0;
        float _S246 = 1.0f - _S244.primal_0 * _S244.primal_0 / 24.0f;
        float _S247 = ((0.0f - (_S245 + _S245) * 0.0416666679084301f) * _S191 - _S246) / (_S191 * _S191);
        k_3 = _S246 / _S191;
        s_diff_k_1 = _S247;
    }
    else
    {
        float _S248 = _S244.differential_0 * 0.5f;
        DiffPair_float_0 _S249;
        (&_S249)->primal_0 = 0.5f * _S244.primal_0;
        (&_S249)->differential_0 = _S248;
        DiffPair_float_0 _S250 = _d_sin_0(&_S249);
        float _S251 = 2.0f * _S250.primal_0;
        float _S252 = (_S250.differential_0 * 2.0f * _S241.primal_0 - _S251 * _S241.differential_0) / (_S241.primal_0 * _S241.primal_0);
        k_3 = _S251 / _S241.primal_0;
        s_diff_k_1 = _S252;
    }
    float2  _S253 = _S187 * make_float2 (k_3);
    float2  _S254 = _S239 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S187;
    float u_15 = _S253.x;
    float s_diff_u_8 = _S254.x;
    float v_15 = _S253.y;
    float s_diff_v_8 = _S254.y;
    float _S255 = s_diff_u_8 * u_15;
    float _S256 = s_diff_v_8 * v_15;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float s_diff_r2_16 = _S255 + _S255 + (_S256 + _S256);
    float _S257 = dist_coeffs_6[int(2)] + r2_15 * dist_coeffs_6[int(3)];
    float _S258 = dist_coeffs_6[int(1)] + r2_15 * _S257;
    float _S259 = dist_coeffs_6[int(0)] + r2_15 * _S258;
    float2  _S260 = _S254 * make_float2 (1.0f + r2_15 * _S259) + make_float2 (s_diff_r2_16 * _S259 + (s_diff_r2_16 * _S258 + (s_diff_r2_16 * _S257 + s_diff_r2_16 * dist_coeffs_6[int(3)] * r2_15) * r2_15) * r2_15) * _S253 + make_float2 (s_diff_u_8 * _S210 * v_15 + s_diff_v_8 * (_S210 * u_15) + (s_diff_r2_16 + (s_diff_u_8 * 2.0f * u_15 + s_diff_u_8 * (2.0f * u_15))) * dist_coeffs_6[int(5)] + s_diff_r2_16 * dist_coeffs_6[int(6)], s_diff_u_8 * _S211 * v_15 + s_diff_v_8 * (_S211 * u_15) + (s_diff_r2_16 + (s_diff_v_8 * 2.0f * v_15 + s_diff_v_8 * (2.0f * v_15))) * dist_coeffs_6[int(4)] + s_diff_r2_16 * dist_coeffs_6[int(7)]);
    float2  _S261 = _S260 + make_float2 (_S260.x * dist_coeffs_6[int(8)] + _S260.y * dist_coeffs_6[int(9)], 0.0f);
    float _S262 = _S261.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->z) = _S261.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->z) = _S262;
    return J_2;
}

inline __device__ bool equirect_proj_nav(float3  p_view_6, float4  intrins_6, FixedArray<float, 10>  dist_coeffs_7, float2  * uv_4)
{
    *uv_4 = make_float2 (intrins_6.x * (F32_atan2((p_view_6.x), (p_view_6.z))) + intrins_6.z, intrins_6.y * (F32_atan2((p_view_6.y), (length_0(float2 {p_view_6.x, p_view_6.z})))) + intrins_6.w);
    return true;
}

inline __device__ Matrix<float, 2, 3>  equirect_proj_jac(float3  p_view_7, float4  intrins_7, FixedArray<float, 10>  dist_coeffs_8)
{
    float _S263 = p_view_7.x;
    float _S264 = p_view_7.z;
    DiffPair_float_0 _S265;
    (&_S265)->primal_0 = _S263;
    (&_S265)->differential_0 = 1.0f;
    DiffPair_float_0 _S266;
    (&_S266)->primal_0 = _S264;
    (&_S266)->differential_0 = 0.0f;
    DiffPair_float_0 _S267 = _d_atan2_0(&_S265, &_S266);
    float _S268 = p_view_7.y;
    float2  _S269 = float2 {p_view_7.x, p_view_7.z};
    float2  _S270 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S271;
    (&_S271)->primal_0 = _S269;
    (&_S271)->differential_0 = _S270;
    DiffPair_float_0 _S272 = s_fwd_length_impl_0(&_S271);
    DiffPair_float_0 _S273;
    (&_S273)->primal_0 = _S268;
    (&_S273)->differential_0 = 0.0f;
    DiffPair_float_0 _S274;
    (&_S274)->primal_0 = _S272.primal_0;
    (&_S274)->differential_0 = _S272.differential_0;
    DiffPair_float_0 _S275 = _d_atan2_0(&_S273, &_S274);
    float fx_3 = intrins_7.x;
    float fy_3 = intrins_7.y;
    float _S276 = _S275.differential_0 * fy_3;
    Matrix<float, 2, 3>  J_3;
    *&(((&J_3)->rows + (int(0)))->x) = _S267.differential_0 * fx_3;
    *&(((&J_3)->rows + (int(1)))->x) = _S276;
    DiffPair_float_0 _S277;
    (&_S277)->primal_0 = _S263;
    (&_S277)->differential_0 = 0.0f;
    DiffPair_float_0 _S278;
    (&_S278)->primal_0 = _S264;
    (&_S278)->differential_0 = 0.0f;
    DiffPair_float_0 _S279 = _d_atan2_0(&_S277, &_S278);
    float2  _S280 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S281;
    (&_S281)->primal_0 = _S269;
    (&_S281)->differential_0 = _S280;
    DiffPair_float_0 _S282 = s_fwd_length_impl_0(&_S281);
    DiffPair_float_0 _S283;
    (&_S283)->primal_0 = _S268;
    (&_S283)->differential_0 = 1.0f;
    DiffPair_float_0 _S284;
    (&_S284)->primal_0 = _S282.primal_0;
    (&_S284)->differential_0 = _S282.differential_0;
    DiffPair_float_0 _S285 = _d_atan2_0(&_S283, &_S284);
    float _S286 = _S285.differential_0 * fy_3;
    *&(((&J_3)->rows + (int(0)))->y) = _S279.differential_0 * fx_3;
    *&(((&J_3)->rows + (int(1)))->y) = _S286;
    DiffPair_float_0 _S287;
    (&_S287)->primal_0 = _S263;
    (&_S287)->differential_0 = 0.0f;
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = _S264;
    (&_S288)->differential_0 = 1.0f;
    DiffPair_float_0 _S289 = _d_atan2_0(&_S287, &_S288);
    float2  _S290 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S291;
    (&_S291)->primal_0 = _S269;
    (&_S291)->differential_0 = _S290;
    DiffPair_float_0 _S292 = s_fwd_length_impl_0(&_S291);
    DiffPair_float_0 _S293;
    (&_S293)->primal_0 = _S268;
    (&_S293)->differential_0 = 0.0f;
    DiffPair_float_0 _S294;
    (&_S294)->primal_0 = _S292.primal_0;
    (&_S294)->differential_0 = _S292.differential_0;
    DiffPair_float_0 _S295 = _d_atan2_0(&_S293, &_S294);
    float _S296 = _S295.differential_0 * fy_3;
    *&(((&J_3)->rows + (int(0)))->z) = _S289.differential_0 * fx_3;
    *&(((&J_3)->rows + (int(1)))->z) = _S296;
    return J_3;
}

inline __device__ float2  distort_point(float2  uv_5, int camera_model_0, FixedArray<float, 10>  dist_coeffs_9)
{
    if(camera_model_0 == int(3))
    {
        return uv_5;
    }
    float k_4;
    float2  _S297;
    if(camera_model_0 == int(1))
    {
        float r_5 = length_0(uv_5);
        float theta_2 = (F32_atan((r_5)));
        if(r_5 < 0.00100000004749745f)
        {
            k_4 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            k_4 = theta_2 / r_5;
        }
        _S297 = uv_5 * make_float2 (k_4);
    }
    else
    {
        if(camera_model_0 == int(2))
        {
            float r_6 = length_0(uv_5);
            float theta_3 = (F32_atan((r_6)));
            if(r_6 < 0.00100000004749745f)
            {
                k_4 = 1.0f - theta_3 * theta_3 / 24.0f;
            }
            else
            {
                k_4 = 2.0f * (F32_sin((0.5f * theta_3))) / r_6;
            }
            _S297 = uv_5 * make_float2 (k_4);
        }
        else
        {
            _S297 = uv_5;
        }
    }
    float u_16 = _S297.x;
    float v_16 = _S297.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S298 = _S297 * make_float2 (1.0f + r2_16 * (dist_coeffs_9[int(0)] + r2_16 * (dist_coeffs_9[int(1)] + r2_16 * (dist_coeffs_9[int(2)] + r2_16 * dist_coeffs_9[int(3)])))) + make_float2 (2.0f * dist_coeffs_9[int(4)] * u_16 * v_16 + dist_coeffs_9[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + dist_coeffs_9[int(6)] * r2_16, 2.0f * dist_coeffs_9[int(5)] * u_16 * v_16 + dist_coeffs_9[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + dist_coeffs_9[int(7)] * r2_16);
    return _S298 + make_float2 (dist_coeffs_9[int(8)] * _S298.x + dist_coeffs_9[int(9)] * _S298.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_6, FixedArray<float, 10>  * dist_coeffs_10, int maxiter_0, float2  * uv_undist_0)
{
    int i_4 = int(0);
    float2  q_0 = uv_6;
    for(;;)
    {
        if(i_4 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S299 = (*dist_coeffs_10)[int(3)];
        float _S300 = (*dist_coeffs_10)[int(4)];
        float _S301 = (*dist_coeffs_10)[int(5)];
        float _S302 = (*dist_coeffs_10)[int(6)];
        float _S303 = (*dist_coeffs_10)[int(7)];
        float _S304 = (*dist_coeffs_10)[int(8)];
        float _S305 = (*dist_coeffs_10)[int(9)];
        float u_17 = q_0.x;
        float v_17 = q_0.y;
        float r2_17 = u_17 * u_17 + v_17 * v_17;
        float _S306 = (*dist_coeffs_10)[int(2)] + r2_17 * (*dist_coeffs_10)[int(3)];
        float _S307 = (*dist_coeffs_10)[int(1)] + r2_17 * _S306;
        float _S308 = (*dist_coeffs_10)[int(0)] + r2_17 * _S307;
        float radial_4 = 1.0f + r2_17 * _S308;
        float _S309 = 2.0f * (*dist_coeffs_10)[int(4)];
        float _S310 = _S309 * u_17;
        float _S311 = 2.0f * u_17;
        float _S312 = 2.0f * (*dist_coeffs_10)[int(5)];
        float _S313 = _S312 * u_17;
        float _S314 = 2.0f * v_17;
        float2  _S315 = q_0 * make_float2 (radial_4) + make_float2 (_S310 * v_17 + (*dist_coeffs_10)[int(5)] * (r2_17 + _S311 * u_17) + (*dist_coeffs_10)[int(6)] * r2_17, _S313 * v_17 + (*dist_coeffs_10)[int(4)] * (r2_17 + _S314 * v_17) + (*dist_coeffs_10)[int(7)] * r2_17);
        float2  r_7 = _S315 + make_float2 ((*dist_coeffs_10)[int(8)] * _S315.x + (*dist_coeffs_10)[int(9)] * _S315.y, 0.0f) - uv_6;
        float _S316 = 0.0f * v_17;
        float s_diff_r2_17 = u_17 + u_17 + (_S316 + _S316);
        float2  _S317 = make_float2 (1.0f, 0.0f) * make_float2 (radial_4) + make_float2 (s_diff_r2_17 * _S308 + (s_diff_r2_17 * _S307 + (s_diff_r2_17 * _S306 + s_diff_r2_17 * _S299 * r2_17) * r2_17) * r2_17) * q_0 + make_float2 (_S309 * v_17 + 0.0f * _S310 + (s_diff_r2_17 + (_S311 + _S311)) * _S301 + s_diff_r2_17 * _S302, _S312 * v_17 + 0.0f * _S313 + (s_diff_r2_17 + (_S316 + 0.0f * _S314)) * _S300 + s_diff_r2_17 * _S303);
        float _S318 = 0.0f * u_17;
        float s_diff_r2_18 = _S318 + _S318 + (v_17 + v_17);
        float2  _S319 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (s_diff_r2_18 * _S308 + (s_diff_r2_18 * _S307 + (s_diff_r2_18 * _S306 + s_diff_r2_18 * _S299 * r2_17) * r2_17) * r2_17) * q_0 + make_float2 (0.0f * _S309 * v_17 + _S310 + (s_diff_r2_18 + (_S318 + 0.0f * _S311)) * _S301 + s_diff_r2_18 * _S302, 0.0f * _S312 * v_17 + _S313 + (s_diff_r2_18 + (_S314 + _S314)) * _S300 + s_diff_r2_18 * _S303);
        Matrix<float, 2, 2>  _S320 = transpose_1(makeMatrix<float, 2, 2> (_S317 + make_float2 (_S317.x * _S304 + _S317.y * _S305, 0.0f), _S319 + make_float2 (_S319.x * _S304 + _S319.y * _S305, 0.0f)));
        float inv_det_0 = 1.0f / (_S320.rows[int(0)].x * _S320.rows[int(1)].y - _S320.rows[int(0)].y * _S320.rows[int(1)].x);
        float _S321 = r_7.x;
        float _S322 = r_7.y;
        float2  q_1 = q_0 - make_float2 ((_S321 * _S320.rows[int(1)].y - _S322 * _S320.rows[int(0)].y) * inv_det_0, (- _S321 * _S320.rows[int(1)].x + _S322 * _S320.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S323 = (*dist_coeffs_10)[int(0)];
    float _S324 = (*dist_coeffs_10)[int(1)];
    float _S325 = (*dist_coeffs_10)[int(2)];
    float _S326 = (*dist_coeffs_10)[int(3)];
    float _S327 = (*dist_coeffs_10)[int(4)];
    float _S328 = (*dist_coeffs_10)[int(5)];
    float _S329 = (*dist_coeffs_10)[int(6)];
    float _S330 = (*dist_coeffs_10)[int(7)];
    float _S331 = (*dist_coeffs_10)[int(8)];
    float _S332 = (*dist_coeffs_10)[int(9)];
    float u_18 = q_0.x;
    float v_18 = q_0.y;
    float _S333 = 0.0f * v_18;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float s_diff_r2_19 = u_18 + u_18 + (_S333 + _S333);
    float _S334 = (*dist_coeffs_10)[int(2)] + r2_18 * (*dist_coeffs_10)[int(3)];
    float _S335 = (*dist_coeffs_10)[int(1)] + r2_18 * _S334;
    float _S336 = (*dist_coeffs_10)[int(0)] + r2_18 * _S335;
    float radial_5 = 1.0f + r2_18 * _S336;
    float _S337 = 2.0f * (*dist_coeffs_10)[int(4)];
    float _S338 = _S337 * u_18;
    float _S339 = 2.0f * u_18;
    float _S340 = 2.0f * (*dist_coeffs_10)[int(5)];
    float _S341 = _S340 * u_18;
    float _S342 = 2.0f * v_18;
    float2  _S343 = make_float2 (1.0f, 0.0f) * make_float2 (radial_5) + make_float2 (s_diff_r2_19 * _S336 + (s_diff_r2_19 * _S335 + (s_diff_r2_19 * _S334 + s_diff_r2_19 * (*dist_coeffs_10)[int(3)] * r2_18) * r2_18) * r2_18) * q_0 + make_float2 (_S337 * v_18 + 0.0f * _S338 + (s_diff_r2_19 + (_S339 + _S339)) * (*dist_coeffs_10)[int(5)] + s_diff_r2_19 * (*dist_coeffs_10)[int(6)], _S340 * v_18 + 0.0f * _S341 + (s_diff_r2_19 + (_S333 + 0.0f * _S342)) * (*dist_coeffs_10)[int(4)] + s_diff_r2_19 * (*dist_coeffs_10)[int(7)]);
    float _S344 = 0.0f * u_18;
    float s_diff_r2_20 = _S344 + _S344 + (v_18 + v_18);
    float2  _S345 = make_float2 (0.0f, 1.0f) * make_float2 (radial_5) + make_float2 (s_diff_r2_20 * _S336 + (s_diff_r2_20 * _S335 + (s_diff_r2_20 * _S334 + s_diff_r2_20 * (*dist_coeffs_10)[int(3)] * r2_18) * r2_18) * r2_18) * q_0 + make_float2 (0.0f * _S337 * v_18 + _S338 + (s_diff_r2_20 + (_S344 + 0.0f * _S339)) * (*dist_coeffs_10)[int(5)] + s_diff_r2_20 * (*dist_coeffs_10)[int(6)], 0.0f * _S340 * v_18 + _S341 + (s_diff_r2_20 + (_S342 + _S342)) * (*dist_coeffs_10)[int(4)] + s_diff_r2_20 * (*dist_coeffs_10)[int(7)]);
    Matrix<float, 2, 2>  _S346 = transpose_1(makeMatrix<float, 2, 2> (_S343 + make_float2 (_S343.x * (*dist_coeffs_10)[int(8)] + _S343.y * (*dist_coeffs_10)[int(9)], 0.0f), _S345 + make_float2 (_S345.x * (*dist_coeffs_10)[int(8)] + _S345.y * (*dist_coeffs_10)[int(9)], 0.0f)));
    bool _S347;
    if((F32_min((determinant_0(_S346)), ((F32_min((_S346.rows[int(0)].x), (_S346.rows[int(1)].y)))))) > 0.0f)
    {
        float u_19 = (*uv_undist_0).x;
        float v_19 = (*uv_undist_0).y;
        float r2_19 = u_19 * u_19 + v_19 * v_19;
        float2  _S348 = *uv_undist_0 * make_float2 (1.0f + r2_19 * (_S323 + r2_19 * (_S324 + r2_19 * (_S325 + r2_19 * _S326)))) + make_float2 (_S337 * u_19 * v_19 + _S328 * (r2_19 + 2.0f * u_19 * u_19) + _S329 * r2_19, _S340 * u_19 * v_19 + _S327 * (r2_19 + 2.0f * v_19 * v_19) + _S330 * r2_19);
        _S347 = (length_0(_S348 + make_float2 (_S331 * _S348.x + _S332 * _S348.y, 0.0f) - uv_6)) < 0.00999999977648258f;
    }
    else
    {
        _S347 = false;
    }
    return _S347;
}

inline __device__ bool undistort_point(float2  uv_7, int camera_model_1, FixedArray<float, 10>  dist_coeffs_11, float2  * uv_undist_1)
{
    float2  _S349 = uv_7;
    if(camera_model_1 == int(3))
    {
        float lon_0 = _S349.x;
        float lat_0 = _S349.y;
        float cl_0 = (F32_cos((lat_0)));
        *uv_undist_1 = make_float2 (cl_0 * (F32_sin((lon_0))), (F32_sin((lat_0)))) / make_float2 ((F32_max((cl_0 * (F32_cos((lon_0)))), (9.999999960041972e-13f))));
        return true;
    }
    FixedArray<float, 10>  _S350 = dist_coeffs_11;
    bool _S351 = undistort_point_0(_S349, &_S350, int(8), &_S349);
    if(!_S351)
    {
        return false;
    }
    float3  raydir_0;
    if(camera_model_1 == int(1))
    {
        float r_8 = length_0(_S349);
        float s_0;
        if(r_8 < 0.00100000004749745f)
        {
            s_0 = 1.0f - r_8 * r_8 / 6.0f;
        }
        else
        {
            s_0 = (F32_sin((r_8))) / r_8;
        }
        raydir_0 = make_float3 ((_S349 * make_float2 (s_0)).x, (_S349 * make_float2 (s_0)).y, (F32_cos((r_8))));
    }
    else
    {
        if(camera_model_1 == int(2))
        {
            float r_9 = length_0(_S349);
            raydir_0 = make_float3 ((_S349 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_9 * r_9)))))))).x, (_S349 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_9 * r_9)))))))).y, 1.0f - 0.5f * r_9 * r_9);
        }
        else
        {
            raydir_0 = make_float3 (_S349.x, _S349.y, 1.0f);
        }
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_8, int camera_model_2, FixedArray<float, 10>  dist_coeffs_12, float3  * raydir_1)
{
    float2  _S352 = uv_8;
    int3  _S353 = make_int3 (int(0));
    float3  _S354 = make_float3 ((float)_S353.x, (float)_S353.y, (float)_S353.z);
    *raydir_1 = _S354;
    if(camera_model_2 == int(3))
    {
        float lon_1 = _S352.x;
        float lat_1 = _S352.y;
        float cl_1 = (F32_cos((lat_1)));
        *raydir_1 = make_float3 (cl_1 * (F32_sin((lon_1))), (F32_sin((lat_1))), cl_1 * (F32_cos((lon_1))));
        return true;
    }
    FixedArray<float, 10>  _S355 = dist_coeffs_12;
    bool _S356 = undistort_point_0(_S352, &_S355, int(8), &_S352);
    if(!_S356)
    {
        return false;
    }
    if(camera_model_2 == int(1))
    {
        float r_10 = length_0(_S352);
        float s_1;
        if(r_10 < 0.00100000004749745f)
        {
            s_1 = 1.0f - r_10 * r_10 / 6.0f;
        }
        else
        {
            s_1 = (F32_sin((r_10))) / r_10;
        }
        *raydir_1 = make_float3 ((_S352 * make_float2 (s_1)).x, (_S352 * make_float2 (s_1)).y, (F32_cos((r_10))));
    }
    else
    {
        if(camera_model_2 == int(2))
        {
            float r_11 = length_0(_S352);
            *raydir_1 = make_float3 ((_S352 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_11 * r_11)))))))).x, (_S352 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_11 * r_11)))))))).y, 1.0f - 0.5f * r_11 * r_11);
        }
        else
        {
            *raydir_1 = make_float3 (_S352.x, _S352.y, 1.0f);
        }
    }
    return true;
}

inline __device__ float3  normalize_0(float3  x_11)
{
    return x_11 / make_float3 (length_1(x_11));
}

inline __device__ bool generate_ray(float2  uv_9, int camera_model_3, FixedArray<float, 10>  dist_coeffs_13, float3  * raydir_2)
{
    float2  _S357 = uv_9;
    if(camera_model_3 == int(3))
    {
        bool _S358;
        if((F32_abs((_S357.x))) > 3.14159274101257324f)
        {
            _S358 = true;
        }
        else
        {
            _S358 = (F32_abs((_S357.y))) > 1.57079637050628662f;
        }
        if(_S358)
        {
            int3  _S359 = make_int3 (int(0));
            float3  _S360 = make_float3 ((float)_S359.x, (float)_S359.y, (float)_S359.z);
            *raydir_2 = _S360;
            return false;
        }
        float lon_2 = _S357.x;
        float lat_2 = _S357.y;
        float cl_2 = (F32_cos((lat_2)));
        *raydir_2 = make_float3 (cl_2 * (F32_sin((lon_2))), (F32_sin((lat_2))), cl_2 * (F32_cos((lon_2))));
        return true;
    }
    FixedArray<float, 10>  _S361 = dist_coeffs_13;
    bool _S362 = undistort_point_0(_S357, &_S361, int(8), &_S357);
    if(!_S362)
    {
        int3  _S363 = make_int3 (int(0));
        float3  _S364 = make_float3 ((float)_S363.x, (float)_S363.y, (float)_S363.z);
        *raydir_2 = _S364;
        return false;
    }
    if(camera_model_3 == int(1))
    {
        float r_12 = length_0(_S357);
        if(r_12 >= 3.14159274101257324f)
        {
            int3  _S365 = make_int3 (int(0));
            float3  _S366 = make_float3 ((float)_S365.x, (float)_S365.y, (float)_S365.z);
            *raydir_2 = _S366;
            return false;
        }
        float s_2;
        if(r_12 < 0.00100000004749745f)
        {
            s_2 = 1.0f - r_12 * r_12 / 6.0f;
        }
        else
        {
            s_2 = (F32_sin((r_12))) / r_12;
        }
        *raydir_2 = make_float3 ((_S357 * make_float2 (s_2)).x, (_S357 * make_float2 (s_2)).y, (F32_cos((r_12))));
    }
    else
    {
        if(camera_model_3 == int(2))
        {
            float r_13 = length_0(_S357);
            if(r_13 >= 2.0f)
            {
                int3  _S367 = make_int3 (int(0));
                float3  _S368 = make_float3 ((float)_S367.x, (float)_S367.y, (float)_S367.z);
                *raydir_2 = _S368;
                return false;
            }
            *raydir_2 = make_float3 ((_S357 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_13 * r_13)))))))).x, (_S357 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_13 * r_13)))))))).y, 1.0f - 0.5f * r_13 * r_13);
        }
        else
        {
            *raydir_2 = make_float3 (_S357.x, _S357.y, 1.0f);
        }
    }
    *raydir_2 = normalize_0(*raydir_2);
    return true;
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_4, float3  dOut_3)
{
    float _S369 = (*right_4).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  right_d_result_2;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = (*left_4).primal_0.x * dOut_3.x;
    float sum_10 = _S369 + (*right_4).primal_0.rows[int(0)].y * dOut_3.y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = (*left_4).primal_0.x * dOut_3.y;
    float sum_11 = sum_10 + (*right_4).primal_0.rows[int(0)].z * dOut_3.z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = (*left_4).primal_0.x * dOut_3.z;
    float3  left_d_result_2;
    *&((&left_d_result_2)->x) = sum_11;
    float _S370 = (*right_4).primal_0.rows[int(1)].x * dOut_3.x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = (*left_4).primal_0.y * dOut_3.x;
    float sum_12 = _S370 + (*right_4).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = (*left_4).primal_0.y * dOut_3.y;
    float sum_13 = sum_12 + (*right_4).primal_0.rows[int(1)].z * dOut_3.z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = (*left_4).primal_0.y * dOut_3.z;
    *&((&left_d_result_2)->y) = sum_13;
    float _S371 = (*right_4).primal_0.rows[int(2)].x * dOut_3.x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = (*left_4).primal_0.z * dOut_3.x;
    float sum_14 = _S371 + (*right_4).primal_0.rows[int(2)].y * dOut_3.y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = (*left_4).primal_0.z * dOut_3.y;
    float sum_15 = sum_14 + (*right_4).primal_0.rows[int(2)].z * dOut_3.z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = (*left_4).primal_0.z * dOut_3.z;
    *&((&left_d_result_2)->z) = sum_15;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_2;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_2;
    return;
}

inline __device__ float3  mul_3(float3  left_5, Matrix<float, 3, 3>  right_5)
{
    float3  result_8;
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
        int i_5 = int(0);
        float sum_16 = 0.0f;
        for(;;)
        {
            if(i_5 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_17 = sum_16 + _slang_vector_get_element(left_5, i_5) * _slang_vector_get_element(right_5.rows[i_5], j_1);
            i_5 = i_5 + int(1);
            sum_16 = sum_17;
        }
        *_slang_vector_get_element_ptr(&result_8, j_1) = sum_16;
        j_1 = j_1 + int(1);
    }
    return result_8;
}

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_0, float3  t_0)
{
    return - mul_3(t_0, R_0);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_1, float3  raydir_3)
{
    return mul_3(raydir_3, R_1);
}

inline __device__ float3  undo_transform_ray_d(Matrix<float, 3, 3>  R_2, float3  raydir_4)
{
    return mul_3(raydir_4, transpose_0(R_2));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S372, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S373, float3  _S374)
{
    _d_mul_1(_S372, _S373, _S374);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S375 = - _s_dOut_0;
    float3  _S376 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S377;
    (&_S377)->primal_0 = (*dpt_0).primal_0;
    (&_S377)->differential_0 = _S376;
    Matrix<float, 3, 3>  _S378 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S379;
    (&_S379)->primal_0 = (*dpR_0).primal_0;
    (&_S379)->differential_0 = _S378;
    s_bwd_prop_mul_0(&_S377, &_S379, _S375);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S377.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S379.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S380, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S381, float3  _S382)
{
    s_bwd_prop_transform_ray_o_0(_S380, _S381, _S382);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S383 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S383;
    float3  _S384 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S384;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S385 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S386;
    (&_S386)->primal_0 = (*dpraydir_0).primal_0;
    (&_S386)->differential_0 = _S385;
    Matrix<float, 3, 3>  _S387 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S388;
    (&_S388)->primal_0 = (*dpR_1).primal_0;
    (&_S388)->differential_0 = _S387;
    s_bwd_prop_mul_0(&_S386, &_S388, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S386.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S388.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S389, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S390, float3  _S391)
{
    s_bwd_prop_transform_ray_d_0(_S389, _S390, _S391);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S392 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S392;
    float3  _S393 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S393;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_4)
{
    float _S394 = (F32_exp(((*dpx_5).primal_0))) * dOut_4;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S394;
    return;
}

inline __device__ float3  exp_0(float3  x_12)
{
    float3  result_9;
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
        *_slang_vector_get_element_ptr(&result_9, i_6) = (F32_exp((_slang_vector_get_element(x_12, i_6))));
        i_6 = i_6 + int(1);
    }
    return result_9;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  dOut_5)
{
    float3  _S395 = exp_0((*dpx_6).primal_0) * dOut_5;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S395;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_5, float3  scale_4)
{
    float x_13 = quat_5.y;
    float x2_5 = x_13 * x_13;
    float y2_5 = quat_5.z * quat_5.z;
    float z2_5 = quat_5.w * quat_5.w;
    float xy_5 = quat_5.y * quat_5.z;
    float xz_5 = quat_5.y * quat_5.w;
    float yz_5 = quat_5.z * quat_5.w;
    float wx_5 = quat_5.x * quat_5.y;
    float wy_5 = quat_5.x * quat_5.z;
    float wz_5 = quat_5.x * quat_5.w;
    float3  _S396 = exp_0(- scale_4);
    return mul_1(makeMatrix<float, 3, 3> (_S396.x, 0.0f, 0.0f, 0.0f, _S396.y, 0.0f, 0.0f, 0.0f, _S396.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)))));
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ float3  s_primal_ctx_exp_0(float3  _S397)
{
    return exp_0(_S397);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S398, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S399, Matrix<float, 3, 3>  _S400)
{
    mul_0(_S398, _S399, _S400);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S401, float3  _S402)
{
    _d_exp_vector_0(_S401, _S402);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_2)
{
    float _S403 = (*dpquat_0).primal_0.y;
    float x2_6 = _S403 * _S403;
    float y2_6 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_6 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_6 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_6 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_6 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    float3  _S404 = - (*dpscale_0).primal_0;
    float3  _S405 = s_primal_ctx_exp_0(_S404);
    Matrix<float, 3, 3>  _S406 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))));
    Matrix<float, 3, 3>  _S407 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S408;
    (&_S408)->primal_0 = makeMatrix<float, 3, 3> (_S405.x, 0.0f, 0.0f, 0.0f, _S405.y, 0.0f, 0.0f, 0.0f, _S405.z);
    (&_S408)->differential_0 = _S407;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S409;
    (&_S409)->primal_0 = _S406;
    (&_S409)->differential_0 = _S407;
    s_bwd_prop_mul_1(&_S408, &_S409, _s_dOut_2);
    Matrix<float, 3, 3>  _S410 = transpose_0(_S409.differential_0);
    float3  _S411 = make_float3 (_S408.differential_0.rows[int(0)].x, _S408.differential_0.rows[int(1)].y, _S408.differential_0.rows[int(2)].z);
    float3  _S412 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S413;
    (&_S413)->primal_0 = _S404;
    (&_S413)->differential_0 = _S412;
    s_bwd_prop_exp_0(&_S413, _S411);
    float3  _S414 = - _S413.differential_0;
    Matrix<float, 3, 3>  _S415 = transpose_0(_S410);
    float _S416 = 2.0f * - _S415.rows[int(2)].z;
    float _S417 = 2.0f * _S415.rows[int(2)].y;
    float _S418 = 2.0f * _S415.rows[int(2)].x;
    float _S419 = 2.0f * _S415.rows[int(1)].z;
    float _S420 = 2.0f * - _S415.rows[int(1)].y;
    float _S421 = 2.0f * _S415.rows[int(1)].x;
    float _S422 = 2.0f * _S415.rows[int(0)].z;
    float _S423 = 2.0f * _S415.rows[int(0)].y;
    float _S424 = 2.0f * - _S415.rows[int(0)].x;
    float _S425 = - _S421 + _S423;
    float _S426 = _S418 + - _S422;
    float _S427 = - _S417 + _S419;
    float _S428 = _S417 + _S419;
    float _S429 = _S418 + _S422;
    float _S430 = _S421 + _S423;
    float _S431 = (*dpquat_0).primal_0.w * (_S420 + _S424);
    float _S432 = (*dpquat_0).primal_0.z * (_S416 + _S424);
    float _S433 = (*dpquat_0).primal_0.y * (_S416 + _S420);
    float _S434 = (*dpquat_0).primal_0.x * _S425 + (*dpquat_0).primal_0.z * _S428 + (*dpquat_0).primal_0.y * _S429 + _S431 + _S431;
    float _S435 = (*dpquat_0).primal_0.x * _S426 + (*dpquat_0).primal_0.w * _S428 + (*dpquat_0).primal_0.y * _S430 + _S432 + _S432;
    float _S436 = (*dpquat_0).primal_0.x * _S427 + (*dpquat_0).primal_0.w * _S429 + (*dpquat_0).primal_0.z * _S430 + _S433 + _S433;
    float _S437 = (*dpquat_0).primal_0.w * _S425 + (*dpquat_0).primal_0.z * _S426 + (*dpquat_0).primal_0.y * _S427;
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S414;
    float4  _S438 = make_float4 (0.0f);
    *&((&_S438)->w) = _S434;
    *&((&_S438)->z) = _S435;
    *&((&_S438)->y) = _S436;
    *&((&_S438)->x) = _S437;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S438;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S439, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S440, Matrix<float, 3, 3>  _S441)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S439, _S440, _S441);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_6, float3  scale_5, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_0, float3  * v_scale_0)
{
    float4  _S442 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_6;
    (&dp_quat_0)->differential_0 = _S442;
    float3  _S443 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_5;
    (&dp_scale_0)->differential_0 = _S443;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_0 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_6)
{
    float _S444 = dOut_6.y;
    float _S445 = dOut_6.z;
    float _S446 = dOut_6.x;
    float _S447 = (*a_0).primal_0.z * _S444 + - (*a_0).primal_0.y * _S445;
    float _S448 = - (*a_0).primal_0.z * _S446 + (*a_0).primal_0.x * _S445;
    float _S449 = (*a_0).primal_0.y * _S446 + - (*a_0).primal_0.x * _S444;
    float3  _S450 = make_float3 (- (*b_0).primal_0.z * _S444 + (*b_0).primal_0.y * _S445, (*b_0).primal_0.z * _S446 + - (*b_0).primal_0.x * _S445, - (*b_0).primal_0.y * _S446 + (*b_0).primal_0.x * _S444);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S450;
    float3  _S451 = make_float3 (_S447, _S448, _S449);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S451;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S452 = left_6.y;
    float _S453 = right_6.z;
    float _S454 = left_6.z;
    float _S455 = right_6.y;
    float _S456 = right_6.x;
    float _S457 = left_6.x;
    return make_float3 (_S452 * _S453 - _S454 * _S455, _S454 * _S456 - _S457 * _S453, _S457 * _S455 - _S452 * _S456);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_0, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_2(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_2(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_0 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S458, float3  _S459)
{
    return mul_2(_S458, _S459);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S460, float3  _S461)
{
    return cross_0(_S460, _S461);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S462, float3  _S463)
{
    return dot_0(_S462, _S463);
}

inline __device__ float s_primal_ctx_exp_1(float _S464)
{
    return (F32_exp((_S464)));
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_float_0 * _S465, float _S466)
{
    _d_exp_0(_S465, _S466);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S467, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S468, float _S469)
{
    _d_dot_0(_S467, _S468, _S469);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S470, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S471, float3  _S472)
{
    _d_cross_0(_S470, _S471, _S472);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S473, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S474, float3  _S475)
{
    _d_mul_0(_S473, _S474, _S475);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    float3  _S476 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S477 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S476);
    float3  _S478 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S479 = s_primal_ctx_cross_0(_S478, _S477);
    float _S480 = -0.5f * s_primal_ctx_dot_0(_S479, _S479);
    float _S481 = s_primal_ctx_dot_0(_S478, _S478);
    float _S482 = _S480 / _S481;
    float _S483 = _S481 * _S481;
    float _S484 = (*dpopacity_0).primal_0 * _s_dOut_3;
    float _S485 = s_primal_ctx_exp_1(_S482) * _s_dOut_3;
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = _S482;
    (&_S486)->differential_0 = 0.0f;
    s_bwd_prop_exp_1(&_S486, _S484);
    float _S487 = _S486.differential_0 / _S483;
    float _S488 = _S480 * - _S487;
    float _S489 = _S481 * _S487;
    float3  _S490 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S491;
    (&_S491)->primal_0 = _S478;
    (&_S491)->differential_0 = _S490;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S492;
    (&_S492)->primal_0 = _S478;
    (&_S492)->differential_0 = _S490;
    s_bwd_prop_dot_0(&_S491, &_S492, _S488);
    float _S493 = -0.5f * _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S494;
    (&_S494)->primal_0 = _S479;
    (&_S494)->differential_0 = _S490;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S495;
    (&_S495)->primal_0 = _S479;
    (&_S495)->differential_0 = _S490;
    s_bwd_prop_dot_0(&_S494, &_S495, _S493);
    float3  _S496 = _S495.differential_0 + _S494.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S497;
    (&_S497)->primal_0 = _S478;
    (&_S497)->differential_0 = _S490;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S498;
    (&_S498)->primal_0 = _S477;
    (&_S498)->differential_0 = _S490;
    s_bwd_prop_cross_0(&_S497, &_S498, _S496);
    float3  _S499 = _S492.differential_0 + _S491.differential_0 + _S497.differential_0;
    Matrix<float, 3, 3>  _S500 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S501;
    (&_S501)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S501)->differential_0 = _S500;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S502;
    (&_S502)->primal_0 = (*dpray_d_0).primal_0;
    (&_S502)->differential_0 = _S490;
    s_bwd_prop_mul_2(&_S501, &_S502, _S499);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S503;
    (&_S503)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S503)->differential_0 = _S500;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S504;
    (&_S504)->primal_0 = _S476;
    (&_S504)->differential_0 = _S490;
    s_bwd_prop_mul_2(&_S503, &_S504, _S498.differential_0);
    float3  _S505 = - _S504.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S502.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S504.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S485;
    Matrix<float, 3, 3>  _S506 = _S501.differential_0 + _S503.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S506;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S505;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S507, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S508, DiffPair_float_0 * _S509, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S510, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S511, float _S512)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S507, _S508, _S509, _S510, _S511, _S512);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_1, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_0, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S513 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S513;
    Matrix<float, 3, 3>  _S514 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S514;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_1;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S513;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S513;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_iscl_rot_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_0 = dp_mean_0.differential_0;
    *v_iscl_rot_1 = dp_iscl_rot_0.differential_0;
    *v_opacity_0 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_2, Matrix<float, 3, 3>  iscl_rot_2, float opacity_2, float3  rgb_0, float3  ray_o_2, float3  ray_d_2, float3  * out_rgb_0, float * depth_0)
{
    *out_rgb_0 = rgb_0;
    float3  grd_1 = mul_2(iscl_rot_2, ray_d_2);
    *depth_0 = - dot_0(mul_2(iscl_rot_2, ray_o_2 - mean_2), grd_1) / dot_0(grd_1, grd_1);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S515 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S516 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S515);
    float3  _S517 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S518 = s_primal_ctx_dot_0(_S517, _S517);
    float _S519 = dpdepth_0 / (_S518 * _S518);
    float _S520 = - s_primal_ctx_dot_0(_S516, _S517) * - _S519;
    float _S521 = _S518 * _S519;
    float3  _S522 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S523;
    (&_S523)->primal_0 = _S517;
    (&_S523)->differential_0 = _S522;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S524;
    (&_S524)->primal_0 = _S517;
    (&_S524)->differential_0 = _S522;
    s_bwd_prop_dot_0(&_S523, &_S524, _S520);
    float _S525 = - _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S526;
    (&_S526)->primal_0 = _S516;
    (&_S526)->differential_0 = _S522;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S527;
    (&_S527)->primal_0 = _S517;
    (&_S527)->differential_0 = _S522;
    s_bwd_prop_dot_0(&_S526, &_S527, _S525);
    float3  _S528 = _S524.differential_0 + _S523.differential_0 + _S527.differential_0;
    Matrix<float, 3, 3>  _S529 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S530;
    (&_S530)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S530)->differential_0 = _S529;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S531;
    (&_S531)->primal_0 = (*dpray_d_1).primal_0;
    (&_S531)->differential_0 = _S522;
    s_bwd_prop_mul_2(&_S530, &_S531, _S528);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S532;
    (&_S532)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S532)->differential_0 = _S529;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S533;
    (&_S533)->primal_0 = _S515;
    (&_S533)->differential_0 = _S522;
    s_bwd_prop_mul_2(&_S532, &_S533, _S526.differential_0);
    float3  _S534 = - _S533.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S531.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S533.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S535 = _S530.differential_0 + _S532.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S535;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S534;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S536, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S537, DiffPair_float_0 * _S538, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S539, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S540, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S541, float3  _S542, float _S543)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S536, _S537, _S538, _S539, _S540, _S541, _S542, _S543);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_3, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_1, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S544 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S544;
    Matrix<float, 3, 3>  _S545 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S545;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_3;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S544;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S544;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S544;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_mean_1 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_1 = dp_opacity_1.differential_0;
    *v_rgb_0 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ float view_radius_3dgs(float3  mean_4, float3  log_scale_0, float logit_opacity_0, float3  campos_0)
{
    float radius_0 = (F32_exp(((F32_max((log_scale_0.x), ((F32_max((log_scale_0.y), (log_scale_0.z))))))))) * (F32_sqrt((2.0f * (F32_log(((F32_max((255.0f / (1.0f + (F32_exp((- logit_opacity_0))))), (1.0f)))))))));
    float dist_0 = length_1(mean_4 - campos_0);
    return radius_0 / ((F32_max((dist_0), (radius_0))) + (F32_sqrt(((F32_max((dist_0 * dist_0 - radius_0 * radius_0), (0.0f)))))));
}

