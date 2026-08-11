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

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ DiffPair_float_0 _d_atan2_0(DiffPair_float_0 * dpy_0, DiffPair_float_0 * dpx_0)
{
    float _S4 = dpx_0->primal_0 * dpx_0->primal_0 + dpy_0->primal_0 * dpy_0->primal_0;
    DiffPair_float_0 _S5 = { (F32_atan2((dpy_0->primal_0), (dpx_0->primal_0))), - dpy_0->primal_0 / _S4 * dpx_0->differential_0 + dpx_0->primal_0 / _S4 * dpy_0->differential_0 };
    return _S5;
}

inline __device__ DiffPair_float_0 _d_sqrt_0(DiffPair_float_0 * dpx_1)
{
    DiffPair_float_0 _S6 = { (F32_sqrt((dpx_1->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_1->primal_0)))))) * dpx_1->differential_0 };
    return _S6;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_2)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_2;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_2).primal_0.x * dOut_2;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_2;
    *&((&y_d_result_0)->y) = (*dpx_2).primal_0.y * dOut_2;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_2;
    *&((&y_d_result_0)->z) = (*dpx_2).primal_0.z * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = x_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float2  x_7, float2  y_0)
{
    int i_2 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_2 < int(2))
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

inline __device__ float dot_1(float3  x_8, float3  y_1)
{
    int i_3 = int(0);
    float result_6 = 0.0f;
    for(;;)
    {
        if(i_3 < int(3))
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
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ float length_1(float3  x_10)
{
    return (F32_sqrt((dot_1(x_10, x_10))));
}

inline __device__ bool equirect_proj_nav(float3  p_view_0, float4  intrins_0, float2  * uv_0)
{
    *uv_0 = make_float2 (intrins_0.x * (F32_atan2((p_view_0.x), (p_view_0.z))) + intrins_0.z, intrins_0.y * (F32_atan2((p_view_0.y), (length_0(float2 {p_view_0.x, p_view_0.z})))) + intrins_0.w);
    return true;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_3)
{
    float _S7 = *&((&dpx_3->differential_0)->x) * *&((&dpx_3->primal_0)->x);
    float _S8 = *&((&dpx_3->differential_0)->y) * *&((&dpx_3->primal_0)->y);
    float s_diff_len_0 = _S7 + _S7 + (_S8 + _S8);
    DiffPair_float_0 _S9;
    (&_S9)->primal_0 = *&((&dpx_3->primal_0)->x) * *&((&dpx_3->primal_0)->x) + *&((&dpx_3->primal_0)->y) * *&((&dpx_3->primal_0)->y);
    (&_S9)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S10 = _d_sqrt_0(&_S9);
    DiffPair_float_0 _S11 = { _S10.primal_0, _S10.differential_0 };
    return _S11;
}

inline __device__ Matrix<float, 2, 3>  equirect_proj_jac(float3  p_view_1, float4  intrins_1)
{
    float _S12 = p_view_1.x;
    float _S13 = p_view_1.z;
    DiffPair_float_0 _S14;
    (&_S14)->primal_0 = _S12;
    (&_S14)->differential_0 = 1.0f;
    DiffPair_float_0 _S15;
    (&_S15)->primal_0 = _S13;
    (&_S15)->differential_0 = 0.0f;
    DiffPair_float_0 _S16 = _d_atan2_0(&_S14, &_S15);
    float _S17 = p_view_1.y;
    float2  _S18 = float2 {p_view_1.x, p_view_1.z};
    float2  _S19 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S20;
    (&_S20)->primal_0 = _S18;
    (&_S20)->differential_0 = _S19;
    DiffPair_float_0 _S21 = s_fwd_length_impl_0(&_S20);
    DiffPair_float_0 _S22;
    (&_S22)->primal_0 = _S17;
    (&_S22)->differential_0 = 0.0f;
    DiffPair_float_0 _S23;
    (&_S23)->primal_0 = _S21.primal_0;
    (&_S23)->differential_0 = _S21.differential_0;
    DiffPair_float_0 _S24 = _d_atan2_0(&_S22, &_S23);
    float fx_0 = intrins_1.x;
    float fy_0 = intrins_1.y;
    float _S25 = _S24.differential_0 * fy_0;
    Matrix<float, 2, 3>  J_0;
    *&(((&J_0)->rows + (int(0)))->x) = _S16.differential_0 * fx_0;
    *&(((&J_0)->rows + (int(1)))->x) = _S25;
    DiffPair_float_0 _S26;
    (&_S26)->primal_0 = _S12;
    (&_S26)->differential_0 = 0.0f;
    DiffPair_float_0 _S27;
    (&_S27)->primal_0 = _S13;
    (&_S27)->differential_0 = 0.0f;
    DiffPair_float_0 _S28 = _d_atan2_0(&_S26, &_S27);
    float2  _S29 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S30;
    (&_S30)->primal_0 = _S18;
    (&_S30)->differential_0 = _S29;
    DiffPair_float_0 _S31 = s_fwd_length_impl_0(&_S30);
    DiffPair_float_0 _S32;
    (&_S32)->primal_0 = _S17;
    (&_S32)->differential_0 = 1.0f;
    DiffPair_float_0 _S33;
    (&_S33)->primal_0 = _S31.primal_0;
    (&_S33)->differential_0 = _S31.differential_0;
    DiffPair_float_0 _S34 = _d_atan2_0(&_S32, &_S33);
    float _S35 = _S34.differential_0 * fy_0;
    *&(((&J_0)->rows + (int(0)))->y) = _S28.differential_0 * fx_0;
    *&(((&J_0)->rows + (int(1)))->y) = _S35;
    DiffPair_float_0 _S36;
    (&_S36)->primal_0 = _S12;
    (&_S36)->differential_0 = 0.0f;
    DiffPair_float_0 _S37;
    (&_S37)->primal_0 = _S13;
    (&_S37)->differential_0 = 1.0f;
    DiffPair_float_0 _S38 = _d_atan2_0(&_S36, &_S37);
    float2  _S39 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S40;
    (&_S40)->primal_0 = _S18;
    (&_S40)->differential_0 = _S39;
    DiffPair_float_0 _S41 = s_fwd_length_impl_0(&_S40);
    DiffPair_float_0 _S42;
    (&_S42)->primal_0 = _S17;
    (&_S42)->differential_0 = 0.0f;
    DiffPair_float_0 _S43;
    (&_S43)->primal_0 = _S41.primal_0;
    (&_S43)->differential_0 = _S41.differential_0;
    DiffPair_float_0 _S44 = _d_atan2_0(&_S42, &_S43);
    float _S45 = _S44.differential_0 * fy_0;
    *&(((&J_0)->rows + (int(0)))->z) = _S38.differential_0 * fx_0;
    *&(((&J_0)->rows + (int(1)))->z) = _S45;
    return J_0;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion_none(float2  uv_1, FixedArray<float, 1>  dist_coeffs_0)
{
    return true;
}

inline __device__ float2  DistNone_distort_0(float2  uv_2, FixedArray<float, 1>  * coeffs_0)
{
    return uv_2;
}

inline __device__ bool persp_proj_nav_none(float3  p_view_2, float4  intrins_2, FixedArray<float, 1>  dist_coeffs_1, float2  * uv_3)
{
    bool _S46;
    for(;;)
    {
        float2  _S47 = float2 {p_view_2.x, p_view_2.y};
        float _S48 = p_view_2.z;
        float2  uv0_0 = _S47 / make_float2 (_S48);
        bool _S49 = _S48 < 0.0f;
        if(_S49)
        {
            *uv_3 = uv0_0;
            _S46 = false;
            break;
        }
        float2  uv_4 = _S47 / make_float2 (_S48);
        FixedArray<float, 1>  _S50 = dist_coeffs_1;
        float2  _S51 = DistNone_distort_0(uv_4, &_S50);
        *uv_3 = make_float2 (intrins_2.x * _S51.x + intrins_2.z, intrins_2.y * _S51.y + intrins_2.w);
        _S46 = true;
        break;
    }
    return _S46;
}

inline __device__ bool fisheye_proj_nav_none(float3  p_view_3, float4  intrins_3, FixedArray<float, 1>  dist_coeffs_2, float2  * uv_5)
{
    float2  _S52 = float2 {p_view_3.x, p_view_3.y};
    float r_3 = length_0(_S52);
    float _S53 = p_view_3.z;
    float theta_0 = (F32_atan2((r_3), (_S53)));
    float k_0;
    if(theta_0 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S53;
    }
    else
    {
        k_0 = theta_0 / r_3;
    }
    float2  _S54 = _S52 * make_float2 (k_0);
    FixedArray<float, 1>  _S55 = dist_coeffs_2;
    float2  _S56 = DistNone_distort_0(_S54, &_S55);
    *uv_5 = make_float2 (intrins_3.x * _S56.x + intrins_3.z, intrins_3.y * _S56.y + intrins_3.w);
    return true;
}

inline __device__ DiffPair_float_0 _d_sin_0(DiffPair_float_0 * dpx_4)
{
    DiffPair_float_0 _S57 = { (F32_sin((dpx_4->primal_0))), (F32_cos((dpx_4->primal_0))) * dpx_4->differential_0 };
    return _S57;
}

inline __device__ bool equisolid_proj_nav_none(float3  p_view_4, float4  intrins_4, FixedArray<float, 1>  dist_coeffs_3, float2  * uv_6)
{
    float2  _S58 = float2 {p_view_4.x, p_view_4.y};
    float r_4 = length_0(_S58);
    float _S59 = p_view_4.z;
    float theta_1 = (F32_atan2((r_4), (_S59)));
    float k_1;
    if(r_4 < 9.99999997475242708e-07f)
    {
        k_1 = (1.0f - theta_1 * theta_1 / 24.0f) / _S59;
    }
    else
    {
        k_1 = 2.0f * (F32_sin((0.5f * theta_1))) / r_4;
    }
    float2  _S60 = _S58 * make_float2 (k_1);
    FixedArray<float, 1>  _S61 = dist_coeffs_3;
    float2  _S62 = DistNone_distort_0(_S60, &_S61);
    *uv_6 = make_float2 (intrins_4.x * _S62.x + intrins_4.z, intrins_4.y * _S62.y + intrins_4.w);
    return true;
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistNone_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_0, FixedArray<float, 1>  * coeffs_1)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S63 = { dpuv_0->primal_0, dpuv_0->differential_0 };
    return _S63;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac_none(float3  p_view_5, float4  intrins_5, FixedArray<float, 1>  dist_coeffs_4)
{
    float2  _S64 = float2 {p_view_5.x, p_view_5.y};
    float _S65 = p_view_5.z;
    float2  _S66 = _S64 * make_float2 (0.0f);
    float _S67 = _S65 * _S65;
    float2  s_diff_uv_0 = (make_float2 (1.0f, 0.0f) * make_float2 (_S65) - _S66) / make_float2 (_S67);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S68;
    (&_S68)->primal_0 = _S64 / make_float2 (_S65);
    (&_S68)->differential_0 = s_diff_uv_0;
    FixedArray<float, 1>  _S69 = dist_coeffs_4;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S70 = s_fwd_DistNone_distort_0(&_S68, &_S69);
    float fx_1 = intrins_5.x;
    float fy_1 = intrins_5.y;
    float _S71 = _S70.differential_0.y * fy_1;
    Matrix<float, 2, 3>  J_1;
    *&(((&J_1)->rows + (int(0)))->x) = _S70.differential_0.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->x) = _S71;
    float2  s_diff_uv_1 = (make_float2 (0.0f, 1.0f) * make_float2 (_S65) - _S66) / make_float2 (_S67);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S72;
    (&_S72)->primal_0 = _S64 / make_float2 (_S65);
    (&_S72)->differential_0 = s_diff_uv_1;
    FixedArray<float, 1>  _S73 = dist_coeffs_4;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S74 = s_fwd_DistNone_distort_0(&_S72, &_S73);
    float _S75 = _S74.differential_0.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->y) = _S74.differential_0.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->y) = _S75;
    float2  s_diff_uv_2 = (make_float2 (0.0f, 0.0f) * make_float2 (_S65) - _S64) / make_float2 (_S67);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S76;
    (&_S76)->primal_0 = _S64 / make_float2 (_S65);
    (&_S76)->differential_0 = s_diff_uv_2;
    FixedArray<float, 1>  _S77 = dist_coeffs_4;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S78 = s_fwd_DistNone_distort_0(&_S76, &_S77);
    float _S79 = _S78.differential_0.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->z) = _S78.differential_0.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->z) = _S79;
    return J_1;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac_none(float3  p_view_6, float4  intrins_6, FixedArray<float, 1>  dist_coeffs_5)
{
    Matrix<float, 2, 3>  J_2;
    float2  _S80 = float2 {p_view_6.x, p_view_6.y};
    float2  _S81 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S82;
    (&_S82)->primal_0 = _S80;
    (&_S82)->differential_0 = _S81;
    DiffPair_float_0 _S83 = s_fwd_length_impl_0(&_S82);
    float _S84 = p_view_6.z;
    DiffPair_float_0 _S85;
    (&_S85)->primal_0 = _S83.primal_0;
    (&_S85)->differential_0 = _S83.differential_0;
    DiffPair_float_0 _S86;
    (&_S86)->primal_0 = _S84;
    (&_S86)->differential_0 = 0.0f;
    DiffPair_float_0 _S87 = _d_atan2_0(&_S85, &_S86);
    float k_2;
    float s_diff_k_0;
    if((_S87.primal_0) < 0.00100000004749745f)
    {
        float _S88 = _S87.differential_0 * _S87.primal_0;
        float _S89 = 1.0f - _S87.primal_0 * _S87.primal_0 / 3.0f;
        float _S90 = ((0.0f - (_S88 + _S88) * 0.3333333432674408f) * _S84 - _S89 * 0.0f) / (_S84 * _S84);
        k_2 = _S89 / _S84;
        s_diff_k_0 = _S90;
    }
    else
    {
        float _S91 = (_S87.differential_0 * _S83.primal_0 - _S87.primal_0 * _S83.differential_0) / (_S83.primal_0 * _S83.primal_0);
        k_2 = _S87.primal_0 / _S83.primal_0;
        s_diff_k_0 = _S91;
    }
    float2  _S92 = _S80 * make_float2 (k_2);
    float2  _S93 = _S81 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S80;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S94;
    (&_S94)->primal_0 = _S92;
    (&_S94)->differential_0 = _S93;
    FixedArray<float, 1>  _S95 = dist_coeffs_5;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S96 = s_fwd_DistNone_distort_0(&_S94, &_S95);
    float fx_2 = intrins_6.x;
    float fy_2 = intrins_6.y;
    float _S97 = _S96.differential_0.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->x) = _S96.differential_0.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->x) = _S97;
    float2  _S98 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S99;
    (&_S99)->primal_0 = _S80;
    (&_S99)->differential_0 = _S98;
    DiffPair_float_0 _S100 = s_fwd_length_impl_0(&_S99);
    DiffPair_float_0 _S101;
    (&_S101)->primal_0 = _S100.primal_0;
    (&_S101)->differential_0 = _S100.differential_0;
    DiffPair_float_0 _S102;
    (&_S102)->primal_0 = _S84;
    (&_S102)->differential_0 = 0.0f;
    DiffPair_float_0 _S103 = _d_atan2_0(&_S101, &_S102);
    if((_S103.primal_0) < 0.00100000004749745f)
    {
        float _S104 = _S103.differential_0 * _S103.primal_0;
        float _S105 = 1.0f - _S103.primal_0 * _S103.primal_0 / 3.0f;
        float _S106 = ((0.0f - (_S104 + _S104) * 0.3333333432674408f) * _S84 - _S105 * 0.0f) / (_S84 * _S84);
        k_2 = _S105 / _S84;
        s_diff_k_0 = _S106;
    }
    else
    {
        float _S107 = (_S103.differential_0 * _S100.primal_0 - _S103.primal_0 * _S100.differential_0) / (_S100.primal_0 * _S100.primal_0);
        k_2 = _S103.primal_0 / _S100.primal_0;
        s_diff_k_0 = _S107;
    }
    float2  _S108 = _S80 * make_float2 (k_2);
    float2  _S109 = _S98 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S80;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S110;
    (&_S110)->primal_0 = _S108;
    (&_S110)->differential_0 = _S109;
    FixedArray<float, 1>  _S111 = dist_coeffs_5;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S112 = s_fwd_DistNone_distort_0(&_S110, &_S111);
    float _S113 = _S112.differential_0.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->y) = _S112.differential_0.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->y) = _S113;
    float2  _S114 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S115;
    (&_S115)->primal_0 = _S80;
    (&_S115)->differential_0 = _S114;
    DiffPair_float_0 _S116 = s_fwd_length_impl_0(&_S115);
    DiffPair_float_0 _S117;
    (&_S117)->primal_0 = _S116.primal_0;
    (&_S117)->differential_0 = _S116.differential_0;
    DiffPair_float_0 _S118;
    (&_S118)->primal_0 = _S84;
    (&_S118)->differential_0 = 1.0f;
    DiffPair_float_0 _S119 = _d_atan2_0(&_S117, &_S118);
    if((_S119.primal_0) < 0.00100000004749745f)
    {
        float _S120 = _S119.differential_0 * _S119.primal_0;
        float _S121 = 1.0f - _S119.primal_0 * _S119.primal_0 / 3.0f;
        float _S122 = ((0.0f - (_S120 + _S120) * 0.3333333432674408f) * _S84 - _S121) / (_S84 * _S84);
        k_2 = _S121 / _S84;
        s_diff_k_0 = _S122;
    }
    else
    {
        float _S123 = (_S119.differential_0 * _S116.primal_0 - _S119.primal_0 * _S116.differential_0) / (_S116.primal_0 * _S116.primal_0);
        k_2 = _S119.primal_0 / _S116.primal_0;
        s_diff_k_0 = _S123;
    }
    float2  _S124 = _S80 * make_float2 (k_2);
    float2  _S125 = _S114 * make_float2 (k_2) + make_float2 (s_diff_k_0) * _S80;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S126;
    (&_S126)->primal_0 = _S124;
    (&_S126)->differential_0 = _S125;
    FixedArray<float, 1>  _S127 = dist_coeffs_5;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S128 = s_fwd_DistNone_distort_0(&_S126, &_S127);
    float _S129 = _S128.differential_0.y * fy_2;
    *&(((&J_2)->rows + (int(0)))->z) = _S128.differential_0.x * fx_2;
    *&(((&J_2)->rows + (int(1)))->z) = _S129;
    return J_2;
}

inline __device__ Matrix<float, 2, 3>  equisolid_proj_jac_none(float3  p_view_7, float4  intrins_7, FixedArray<float, 1>  dist_coeffs_6)
{
    Matrix<float, 2, 3>  J_3;
    float2  _S130 = float2 {p_view_7.x, p_view_7.y};
    float2  _S131 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S132;
    (&_S132)->primal_0 = _S130;
    (&_S132)->differential_0 = _S131;
    DiffPair_float_0 _S133 = s_fwd_length_impl_0(&_S132);
    float _S134 = p_view_7.z;
    DiffPair_float_0 _S135;
    (&_S135)->primal_0 = _S133.primal_0;
    (&_S135)->differential_0 = _S133.differential_0;
    DiffPair_float_0 _S136;
    (&_S136)->primal_0 = _S134;
    (&_S136)->differential_0 = 0.0f;
    DiffPair_float_0 _S137 = _d_atan2_0(&_S135, &_S136);
    float k_3;
    float s_diff_k_1;
    if((_S133.primal_0) < 9.99999997475242708e-07f)
    {
        float _S138 = _S137.differential_0 * _S137.primal_0;
        float _S139 = 1.0f - _S137.primal_0 * _S137.primal_0 / 24.0f;
        float _S140 = ((0.0f - (_S138 + _S138) * 0.0416666679084301f) * _S134 - _S139 * 0.0f) / (_S134 * _S134);
        k_3 = _S139 / _S134;
        s_diff_k_1 = _S140;
    }
    else
    {
        float _S141 = _S137.differential_0 * 0.5f;
        DiffPair_float_0 _S142;
        (&_S142)->primal_0 = 0.5f * _S137.primal_0;
        (&_S142)->differential_0 = _S141;
        DiffPair_float_0 _S143 = _d_sin_0(&_S142);
        float _S144 = 2.0f * _S143.primal_0;
        float _S145 = (_S143.differential_0 * 2.0f * _S133.primal_0 - _S144 * _S133.differential_0) / (_S133.primal_0 * _S133.primal_0);
        k_3 = _S144 / _S133.primal_0;
        s_diff_k_1 = _S145;
    }
    float2  _S146 = _S130 * make_float2 (k_3);
    float2  _S147 = _S131 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S130;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S148;
    (&_S148)->primal_0 = _S146;
    (&_S148)->differential_0 = _S147;
    FixedArray<float, 1>  _S149 = dist_coeffs_6;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S150 = s_fwd_DistNone_distort_0(&_S148, &_S149);
    float fx_3 = intrins_7.x;
    float fy_3 = intrins_7.y;
    float _S151 = _S150.differential_0.y * fy_3;
    *&(((&J_3)->rows + (int(0)))->x) = _S150.differential_0.x * fx_3;
    *&(((&J_3)->rows + (int(1)))->x) = _S151;
    float2  _S152 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S153;
    (&_S153)->primal_0 = _S130;
    (&_S153)->differential_0 = _S152;
    DiffPair_float_0 _S154 = s_fwd_length_impl_0(&_S153);
    DiffPair_float_0 _S155;
    (&_S155)->primal_0 = _S154.primal_0;
    (&_S155)->differential_0 = _S154.differential_0;
    DiffPair_float_0 _S156;
    (&_S156)->primal_0 = _S134;
    (&_S156)->differential_0 = 0.0f;
    DiffPair_float_0 _S157 = _d_atan2_0(&_S155, &_S156);
    if((_S154.primal_0) < 9.99999997475242708e-07f)
    {
        float _S158 = _S157.differential_0 * _S157.primal_0;
        float _S159 = 1.0f - _S157.primal_0 * _S157.primal_0 / 24.0f;
        float _S160 = ((0.0f - (_S158 + _S158) * 0.0416666679084301f) * _S134 - _S159 * 0.0f) / (_S134 * _S134);
        k_3 = _S159 / _S134;
        s_diff_k_1 = _S160;
    }
    else
    {
        float _S161 = _S157.differential_0 * 0.5f;
        DiffPair_float_0 _S162;
        (&_S162)->primal_0 = 0.5f * _S157.primal_0;
        (&_S162)->differential_0 = _S161;
        DiffPair_float_0 _S163 = _d_sin_0(&_S162);
        float _S164 = 2.0f * _S163.primal_0;
        float _S165 = (_S163.differential_0 * 2.0f * _S154.primal_0 - _S164 * _S154.differential_0) / (_S154.primal_0 * _S154.primal_0);
        k_3 = _S164 / _S154.primal_0;
        s_diff_k_1 = _S165;
    }
    float2  _S166 = _S130 * make_float2 (k_3);
    float2  _S167 = _S152 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S130;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S168;
    (&_S168)->primal_0 = _S166;
    (&_S168)->differential_0 = _S167;
    FixedArray<float, 1>  _S169 = dist_coeffs_6;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S170 = s_fwd_DistNone_distort_0(&_S168, &_S169);
    float _S171 = _S170.differential_0.y * fy_3;
    *&(((&J_3)->rows + (int(0)))->y) = _S170.differential_0.x * fx_3;
    *&(((&J_3)->rows + (int(1)))->y) = _S171;
    float2  _S172 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S173;
    (&_S173)->primal_0 = _S130;
    (&_S173)->differential_0 = _S172;
    DiffPair_float_0 _S174 = s_fwd_length_impl_0(&_S173);
    DiffPair_float_0 _S175;
    (&_S175)->primal_0 = _S174.primal_0;
    (&_S175)->differential_0 = _S174.differential_0;
    DiffPair_float_0 _S176;
    (&_S176)->primal_0 = _S134;
    (&_S176)->differential_0 = 1.0f;
    DiffPair_float_0 _S177 = _d_atan2_0(&_S175, &_S176);
    if((_S174.primal_0) < 9.99999997475242708e-07f)
    {
        float _S178 = _S177.differential_0 * _S177.primal_0;
        float _S179 = 1.0f - _S177.primal_0 * _S177.primal_0 / 24.0f;
        float _S180 = ((0.0f - (_S178 + _S178) * 0.0416666679084301f) * _S134 - _S179) / (_S134 * _S134);
        k_3 = _S179 / _S134;
        s_diff_k_1 = _S180;
    }
    else
    {
        float _S181 = _S177.differential_0 * 0.5f;
        DiffPair_float_0 _S182;
        (&_S182)->primal_0 = 0.5f * _S177.primal_0;
        (&_S182)->differential_0 = _S181;
        DiffPair_float_0 _S183 = _d_sin_0(&_S182);
        float _S184 = 2.0f * _S183.primal_0;
        float _S185 = (_S183.differential_0 * 2.0f * _S174.primal_0 - _S184 * _S174.differential_0) / (_S174.primal_0 * _S174.primal_0);
        k_3 = _S184 / _S174.primal_0;
        s_diff_k_1 = _S185;
    }
    float2  _S186 = _S130 * make_float2 (k_3);
    float2  _S187 = _S172 * make_float2 (k_3) + make_float2 (s_diff_k_1) * _S130;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S188;
    (&_S188)->primal_0 = _S186;
    (&_S188)->differential_0 = _S187;
    FixedArray<float, 1>  _S189 = dist_coeffs_6;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S190 = s_fwd_DistNone_distort_0(&_S188, &_S189);
    float _S191 = _S190.differential_0.y * fy_3;
    *&(((&J_3)->rows + (int(0)))->z) = _S190.differential_0.x * fx_3;
    *&(((&J_3)->rows + (int(1)))->z) = _S191;
    return J_3;
}

inline __device__ float2  distort_point_none(float2  uv_7, int camera_model_0, FixedArray<float, 1>  dist_coeffs_7)
{
    float2  _S192;
    for(;;)
    {
        if(camera_model_0 == int(3))
        {
            _S192 = uv_7;
            break;
        }
        float k_4;
        if(camera_model_0 == int(1))
        {
            float r_5 = length_0(uv_7);
            float theta_2 = (F32_atan((r_5)));
            if(r_5 < 0.00100000004749745f)
            {
                k_4 = 1.0f - theta_2 * theta_2 / 6.0f;
            }
            else
            {
                k_4 = theta_2 / r_5;
            }
            _S192 = uv_7 * make_float2 (k_4);
        }
        else
        {
            if(camera_model_0 == int(2))
            {
                float r_6 = length_0(uv_7);
                float theta_3 = (F32_atan((r_6)));
                if(r_6 < 0.00100000004749745f)
                {
                    k_4 = 1.0f - theta_3 * theta_3 / 24.0f;
                }
                else
                {
                    k_4 = 2.0f * (F32_sin((0.5f * theta_3))) / r_6;
                }
                _S192 = uv_7 * make_float2 (k_4);
            }
            else
            {
                _S192 = uv_7;
            }
        }
        FixedArray<float, 1>  _S193 = dist_coeffs_7;
        float2  _S194 = DistNone_distort_0(_S192, &_S193);
        _S192 = _S194;
        break;
    }
    return _S192;
}

inline __device__ bool undistort_point_0(float2  uv_8, FixedArray<float, 1>  * dist_coeffs_8, int maxiter_0, float2  * uv_undist_0)
{
    *uv_undist_0 = uv_8;
    return true;
}

inline __device__ float2  DistOpenCV_distort_0(float2  uv_9, FixedArray<float, 4>  * coeffs_2)
{
    float u_0 = uv_9.x;
    float v_0 = uv_9.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    return uv_9 * make_float2 (1.0f + r2_0 * ((*coeffs_2)[int(0)] + r2_0 * (*coeffs_2)[int(1)])) + make_float2 (2.0f * (*coeffs_2)[int(2)] * u_0 * v_0 + (*coeffs_2)[int(3)] * (r2_0 + 2.0f * u_0 * u_0), 2.0f * (*coeffs_2)[int(3)] * u_0 * v_0 + (*coeffs_2)[int(2)] * (r2_0 + 2.0f * v_0 * v_0));
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistOpenCV_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_1, FixedArray<float, 4>  * coeffs_3)
{
    float u_1 = dpuv_1->primal_0.x;
    float s_diff_u_0 = dpuv_1->differential_0.x;
    float v_1 = dpuv_1->primal_0.y;
    float s_diff_v_0 = dpuv_1->differential_0.y;
    float _S195 = s_diff_u_0 * u_1;
    float _S196 = s_diff_v_0 * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_0 = _S195 + _S195 + (_S196 + _S196);
    float _S197 = (*coeffs_3)[int(0)] + r2_1 * (*coeffs_3)[int(1)];
    float radial_0 = 1.0f + r2_1 * _S197;
    float _S198 = 2.0f * (*coeffs_3)[int(2)];
    float _S199 = _S198 * u_1;
    float _S200 = 2.0f * u_1;
    float _S201 = 2.0f * (*coeffs_3)[int(3)];
    float _S202 = _S201 * u_1;
    float _S203 = 2.0f * v_1;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S204 = { dpuv_1->primal_0 * make_float2 (radial_0) + make_float2 (_S199 * v_1 + (*coeffs_3)[int(3)] * (r2_1 + _S200 * u_1), _S202 * v_1 + (*coeffs_3)[int(2)] * (r2_1 + _S203 * v_1)), dpuv_1->differential_0 * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S197 + s_diff_r2_0 * (*coeffs_3)[int(1)] * r2_1) * dpuv_1->primal_0 + make_float2 (s_diff_u_0 * _S198 * v_1 + s_diff_v_0 * _S199 + (s_diff_r2_0 + (s_diff_u_0 * 2.0f * u_1 + s_diff_u_0 * _S200)) * (*coeffs_3)[int(3)], s_diff_u_0 * _S201 * v_1 + s_diff_v_0 * _S202 + (s_diff_r2_0 + (s_diff_v_0 * 2.0f * v_1 + s_diff_v_0 * _S203)) * (*coeffs_3)[int(2)]) };
    return _S204;
}

inline __device__ bool undistort_point_1(float2  uv_10, FixedArray<float, 4>  * dist_coeffs_9, int maxiter_1, float2  * uv_undist_1)
{
    int i_4 = int(0);
    float2  q_0 = uv_10;
    for(;;)
    {
        if(i_4 < maxiter_1)
        {
        }
        else
        {
            break;
        }
        float2  _S205 = DistOpenCV_distort_0(q_0, dist_coeffs_9);
        float2  r_7 = _S205 - uv_10;
        float2  _S206 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S207;
        (&_S207)->primal_0 = q_0;
        (&_S207)->differential_0 = _S206;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S208 = s_fwd_DistOpenCV_distort_0(&_S207, dist_coeffs_9);
        float2  _S209 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S210;
        (&_S210)->primal_0 = q_0;
        (&_S210)->differential_0 = _S209;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S211 = s_fwd_DistOpenCV_distort_0(&_S210, dist_coeffs_9);
        Matrix<float, 2, 2>  _S212 = transpose_1(makeMatrix<float, 2, 2> (_S208.differential_0, _S211.differential_0));
        float inv_det_0 = 1.0f / (_S212.rows[int(0)].x * _S212.rows[int(1)].y - _S212.rows[int(0)].y * _S212.rows[int(1)].x);
        float _S213 = r_7.x;
        float _S214 = r_7.y;
        float2  q_1 = q_0 - make_float2 ((_S213 * _S212.rows[int(1)].y - _S214 * _S212.rows[int(0)].y) * inv_det_0, (- _S213 * _S212.rows[int(1)].x + _S214 * _S212.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_1 = q_0;
    float2  _S215 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S216;
    (&_S216)->primal_0 = q_0;
    (&_S216)->differential_0 = _S215;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S217 = s_fwd_DistOpenCV_distort_0(&_S216, dist_coeffs_9);
    float2  _S218 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S219;
    (&_S219)->primal_0 = q_0;
    (&_S219)->differential_0 = _S218;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S220 = s_fwd_DistOpenCV_distort_0(&_S219, dist_coeffs_9);
    Matrix<float, 2, 2>  _S221 = transpose_1(makeMatrix<float, 2, 2> (_S217.differential_0, _S220.differential_0));
    float _S222 = (F32_min((determinant_0(_S221)), ((F32_min((_S221.rows[int(0)].x), (_S221.rows[int(1)].y))))));
    bool _S223;
    if(_S222 > 0.25f)
    {
        _S223 = _S222 < 4.0f;
    }
    else
    {
        _S223 = false;
    }
    if(_S223)
    {
        float2  _S224 = DistOpenCV_distort_0(q_0, dist_coeffs_9);
        _S223 = (dot_0(q_0, _S224)) >= 0.0f;
    }
    else
    {
        _S223 = false;
    }
    if(_S223)
    {
        float2  _S225 = DistOpenCV_distort_0(*uv_undist_1, dist_coeffs_9);
        _S223 = (length_0(_S225 - uv_10)) < 0.00999999977648258f;
    }
    else
    {
        _S223 = false;
    }
    return _S223;
}

inline __device__ float2  DistThinPrism_distort_0(float2  uv_11, FixedArray<float, 8>  * coeffs_4)
{
    float u_2 = uv_11.x;
    float v_2 = uv_11.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    return uv_11 * make_float2 (1.0f + r2_2 * ((*coeffs_4)[int(0)] + r2_2 * ((*coeffs_4)[int(1)] + r2_2 * ((*coeffs_4)[int(2)] + r2_2 * (*coeffs_4)[int(3)])))) + make_float2 (2.0f * (*coeffs_4)[int(4)] * u_2 * v_2 + (*coeffs_4)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*coeffs_4)[int(6)] * r2_2, 2.0f * (*coeffs_4)[int(5)] * u_2 * v_2 + (*coeffs_4)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*coeffs_4)[int(7)] * r2_2);
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistThinPrism_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_2, FixedArray<float, 8>  * coeffs_5)
{
    float u_3 = dpuv_2->primal_0.x;
    float s_diff_u_1 = dpuv_2->differential_0.x;
    float v_3 = dpuv_2->primal_0.y;
    float s_diff_v_1 = dpuv_2->differential_0.y;
    float _S226 = s_diff_u_1 * u_3;
    float _S227 = s_diff_v_1 * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_1 = _S226 + _S226 + (_S227 + _S227);
    float _S228 = (*coeffs_5)[int(2)] + r2_3 * (*coeffs_5)[int(3)];
    float _S229 = (*coeffs_5)[int(1)] + r2_3 * _S228;
    float _S230 = (*coeffs_5)[int(0)] + r2_3 * _S229;
    float radial_1 = 1.0f + r2_3 * _S230;
    float _S231 = 2.0f * (*coeffs_5)[int(4)];
    float _S232 = _S231 * u_3;
    float _S233 = 2.0f * u_3;
    float _S234 = 2.0f * (*coeffs_5)[int(5)];
    float _S235 = _S234 * u_3;
    float _S236 = 2.0f * v_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S237 = { dpuv_2->primal_0 * make_float2 (radial_1) + make_float2 (_S232 * v_3 + (*coeffs_5)[int(5)] * (r2_3 + _S233 * u_3) + (*coeffs_5)[int(6)] * r2_3, _S235 * v_3 + (*coeffs_5)[int(4)] * (r2_3 + _S236 * v_3) + (*coeffs_5)[int(7)] * r2_3), dpuv_2->differential_0 * make_float2 (radial_1) + make_float2 (s_diff_r2_1 * _S230 + (s_diff_r2_1 * _S229 + (s_diff_r2_1 * _S228 + s_diff_r2_1 * (*coeffs_5)[int(3)] * r2_3) * r2_3) * r2_3) * dpuv_2->primal_0 + make_float2 (s_diff_u_1 * _S231 * v_3 + s_diff_v_1 * _S232 + (s_diff_r2_1 + (s_diff_u_1 * 2.0f * u_3 + s_diff_u_1 * _S233)) * (*coeffs_5)[int(5)] + s_diff_r2_1 * (*coeffs_5)[int(6)], s_diff_u_1 * _S234 * v_3 + s_diff_v_1 * _S235 + (s_diff_r2_1 + (s_diff_v_1 * 2.0f * v_3 + s_diff_v_1 * _S236)) * (*coeffs_5)[int(4)] + s_diff_r2_1 * (*coeffs_5)[int(7)]) };
    return _S237;
}

inline __device__ bool undistort_point_2(float2  uv_12, FixedArray<float, 8>  * dist_coeffs_10, int maxiter_2, float2  * uv_undist_2)
{
    int i_5 = int(0);
    float2  q_2 = uv_12;
    for(;;)
    {
        if(i_5 < maxiter_2)
        {
        }
        else
        {
            break;
        }
        float2  _S238 = DistThinPrism_distort_0(q_2, dist_coeffs_10);
        float2  r_8 = _S238 - uv_12;
        float2  _S239 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S240;
        (&_S240)->primal_0 = q_2;
        (&_S240)->differential_0 = _S239;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S241 = s_fwd_DistThinPrism_distort_0(&_S240, dist_coeffs_10);
        float2  _S242 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S243;
        (&_S243)->primal_0 = q_2;
        (&_S243)->differential_0 = _S242;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S244 = s_fwd_DistThinPrism_distort_0(&_S243, dist_coeffs_10);
        Matrix<float, 2, 2>  _S245 = transpose_1(makeMatrix<float, 2, 2> (_S241.differential_0, _S244.differential_0));
        float inv_det_1 = 1.0f / (_S245.rows[int(0)].x * _S245.rows[int(1)].y - _S245.rows[int(0)].y * _S245.rows[int(1)].x);
        float _S246 = r_8.x;
        float _S247 = r_8.y;
        float2  q_3 = q_2 - make_float2 ((_S246 * _S245.rows[int(1)].y - _S247 * _S245.rows[int(0)].y) * inv_det_1, (- _S246 * _S245.rows[int(1)].x + _S247 * _S245.rows[int(0)].x) * inv_det_1);
        i_5 = i_5 + int(1);
        q_2 = q_3;
    }
    *uv_undist_2 = q_2;
    float2  _S248 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S249;
    (&_S249)->primal_0 = q_2;
    (&_S249)->differential_0 = _S248;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S250 = s_fwd_DistThinPrism_distort_0(&_S249, dist_coeffs_10);
    float2  _S251 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S252;
    (&_S252)->primal_0 = q_2;
    (&_S252)->differential_0 = _S251;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S253 = s_fwd_DistThinPrism_distort_0(&_S252, dist_coeffs_10);
    Matrix<float, 2, 2>  _S254 = transpose_1(makeMatrix<float, 2, 2> (_S250.differential_0, _S253.differential_0));
    float _S255 = (F32_min((determinant_0(_S254)), ((F32_min((_S254.rows[int(0)].x), (_S254.rows[int(1)].y))))));
    bool _S256;
    if(_S255 > 0.25f)
    {
        _S256 = _S255 < 4.0f;
    }
    else
    {
        _S256 = false;
    }
    if(_S256)
    {
        float2  _S257 = DistThinPrism_distort_0(q_2, dist_coeffs_10);
        _S256 = (dot_0(q_2, _S257)) >= 0.0f;
    }
    else
    {
        _S256 = false;
    }
    if(_S256)
    {
        float2  _S258 = DistThinPrism_distort_0(*uv_undist_2, dist_coeffs_10);
        _S256 = (length_0(_S258 - uv_12)) < 0.00999999977648258f;
    }
    else
    {
        _S256 = false;
    }
    return _S256;
}

inline __device__ float2  DistRational_distort_0(float2  uv_13, FixedArray<float, 8>  * coeffs_6)
{
    float u_4 = uv_13.x;
    float v_4 = uv_13.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    return uv_13 * make_float2 ((1.0f + r2_4 * ((*coeffs_6)[int(0)] + r2_4 * ((*coeffs_6)[int(1)] + r2_4 * (*coeffs_6)[int(2)]))) / (1.0f + r2_4 * ((*coeffs_6)[int(3)] + r2_4 * ((*coeffs_6)[int(4)] + r2_4 * (*coeffs_6)[int(5)])))) + make_float2 (2.0f * (*coeffs_6)[int(6)] * u_4 * v_4 + (*coeffs_6)[int(7)] * (r2_4 + 2.0f * u_4 * u_4), 2.0f * (*coeffs_6)[int(7)] * u_4 * v_4 + (*coeffs_6)[int(6)] * (r2_4 + 2.0f * v_4 * v_4));
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistRational_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_3, FixedArray<float, 8>  * coeffs_7)
{
    float u_5 = dpuv_3->primal_0.x;
    float s_diff_u_2 = dpuv_3->differential_0.x;
    float v_5 = dpuv_3->primal_0.y;
    float s_diff_v_2 = dpuv_3->differential_0.y;
    float _S259 = s_diff_u_2 * u_5;
    float _S260 = s_diff_v_2 * v_5;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float s_diff_r2_2 = _S259 + _S259 + (_S260 + _S260);
    float _S261 = (*coeffs_7)[int(1)] + r2_5 * (*coeffs_7)[int(2)];
    float _S262 = (*coeffs_7)[int(0)] + r2_5 * _S261;
    float _S263 = 1.0f + r2_5 * _S262;
    float _S264 = (*coeffs_7)[int(4)] + r2_5 * (*coeffs_7)[int(5)];
    float _S265 = (*coeffs_7)[int(3)] + r2_5 * _S264;
    float _S266 = 1.0f + r2_5 * _S265;
    float radial_2 = _S263 / _S266;
    float _S267 = 2.0f * (*coeffs_7)[int(6)];
    float _S268 = _S267 * u_5;
    float _S269 = 2.0f * u_5;
    float _S270 = 2.0f * (*coeffs_7)[int(7)];
    float _S271 = _S270 * u_5;
    float _S272 = 2.0f * v_5;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S273 = { dpuv_3->primal_0 * make_float2 (radial_2) + make_float2 (_S268 * v_5 + (*coeffs_7)[int(7)] * (r2_5 + _S269 * u_5), _S271 * v_5 + (*coeffs_7)[int(6)] * (r2_5 + _S272 * v_5)), dpuv_3->differential_0 * make_float2 (radial_2) + make_float2 (((s_diff_r2_2 * _S262 + (s_diff_r2_2 * _S261 + s_diff_r2_2 * (*coeffs_7)[int(2)] * r2_5) * r2_5) * _S266 - _S263 * (s_diff_r2_2 * _S265 + (s_diff_r2_2 * _S264 + s_diff_r2_2 * (*coeffs_7)[int(5)] * r2_5) * r2_5)) / (_S266 * _S266)) * dpuv_3->primal_0 + make_float2 (s_diff_u_2 * _S267 * v_5 + s_diff_v_2 * _S268 + (s_diff_r2_2 + (s_diff_u_2 * 2.0f * u_5 + s_diff_u_2 * _S269)) * (*coeffs_7)[int(7)], s_diff_u_2 * _S270 * v_5 + s_diff_v_2 * _S271 + (s_diff_r2_2 + (s_diff_v_2 * 2.0f * v_5 + s_diff_v_2 * _S272)) * (*coeffs_7)[int(6)]) };
    return _S273;
}

inline __device__ bool undistort_point_3(float2  uv_14, FixedArray<float, 8>  * dist_coeffs_11, int maxiter_3, float2  * uv_undist_3)
{
    int i_6 = int(0);
    float2  q_4 = uv_14;
    for(;;)
    {
        if(i_6 < maxiter_3)
        {
        }
        else
        {
            break;
        }
        float2  _S274 = DistRational_distort_0(q_4, dist_coeffs_11);
        float2  r_9 = _S274 - uv_14;
        float2  _S275 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S276;
        (&_S276)->primal_0 = q_4;
        (&_S276)->differential_0 = _S275;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S277 = s_fwd_DistRational_distort_0(&_S276, dist_coeffs_11);
        float2  _S278 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S279;
        (&_S279)->primal_0 = q_4;
        (&_S279)->differential_0 = _S278;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S280 = s_fwd_DistRational_distort_0(&_S279, dist_coeffs_11);
        Matrix<float, 2, 2>  _S281 = transpose_1(makeMatrix<float, 2, 2> (_S277.differential_0, _S280.differential_0));
        float inv_det_2 = 1.0f / (_S281.rows[int(0)].x * _S281.rows[int(1)].y - _S281.rows[int(0)].y * _S281.rows[int(1)].x);
        float _S282 = r_9.x;
        float _S283 = r_9.y;
        float2  q_5 = q_4 - make_float2 ((_S282 * _S281.rows[int(1)].y - _S283 * _S281.rows[int(0)].y) * inv_det_2, (- _S282 * _S281.rows[int(1)].x + _S283 * _S281.rows[int(0)].x) * inv_det_2);
        i_6 = i_6 + int(1);
        q_4 = q_5;
    }
    *uv_undist_3 = q_4;
    float2  _S284 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S285;
    (&_S285)->primal_0 = q_4;
    (&_S285)->differential_0 = _S284;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S286 = s_fwd_DistRational_distort_0(&_S285, dist_coeffs_11);
    float2  _S287 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S288;
    (&_S288)->primal_0 = q_4;
    (&_S288)->differential_0 = _S287;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S289 = s_fwd_DistRational_distort_0(&_S288, dist_coeffs_11);
    Matrix<float, 2, 2>  _S290 = transpose_1(makeMatrix<float, 2, 2> (_S286.differential_0, _S289.differential_0));
    float _S291 = (F32_min((determinant_0(_S290)), ((F32_min((_S290.rows[int(0)].x), (_S290.rows[int(1)].y))))));
    bool _S292;
    if(_S291 > 0.25f)
    {
        _S292 = _S291 < 4.0f;
    }
    else
    {
        _S292 = false;
    }
    if(_S292)
    {
        float2  _S293 = DistRational_distort_0(q_4, dist_coeffs_11);
        _S292 = (dot_0(q_4, _S293)) >= 0.0f;
    }
    else
    {
        _S292 = false;
    }
    if(_S292)
    {
        float2  _S294 = DistRational_distort_0(*uv_undist_3, dist_coeffs_11);
        _S292 = (length_0(_S294 - uv_14)) < 0.00999999977648258f;
    }
    else
    {
        _S292 = false;
    }
    return _S292;
}

inline __device__ bool undistort_point_none(float2  uv_15, int camera_model_1, FixedArray<float, 1>  dist_coeffs_12, float2  * uv_undist_4)
{
    bool _S295;
    for(;;)
    {
        *uv_undist_4 = make_float2 (0.0f);
        if(camera_model_1 == int(3))
        {
            float lon_0 = uv_15.x;
            float lat_0 = uv_15.y;
            float cl_0 = (F32_cos((lat_0)));
            *uv_undist_4 = make_float2 (cl_0 * (F32_sin((lon_0))), (F32_sin((lat_0)))) / make_float2 ((F32_max((cl_0 * (F32_cos((lon_0)))), (9.999999960041972e-13f))));
            _S295 = true;
            break;
        }
        FixedArray<float, 1>  _S296 = dist_coeffs_12;
        float2  uv_u_0;
        bool _S297 = undistort_point_0(uv_15, &_S296, int(8), &uv_u_0);
        if(!_S297)
        {
            _S295 = false;
            break;
        }
        float2  _S298 = uv_u_0;
        float3  raydir_0;
        if(camera_model_1 == int(1))
        {
            float r_10 = length_0(_S298);
            float s_0;
            if(r_10 < 0.00100000004749745f)
            {
                s_0 = 1.0f - r_10 * r_10 / 6.0f;
            }
            else
            {
                s_0 = (F32_sin((r_10))) / r_10;
            }
            raydir_0 = make_float3 ((_S298 * make_float2 (s_0)).x, (_S298 * make_float2 (s_0)).y, (F32_cos((r_10))));
        }
        else
        {
            if(camera_model_1 == int(2))
            {
                float r_11 = length_0(_S298);
                raydir_0 = make_float3 ((_S298 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_11 * r_11)))))))).x, (_S298 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_11 * r_11)))))))).y, 1.0f - 0.5f * r_11 * r_11);
            }
            else
            {
                raydir_0 = make_float3 (_S298.x, _S298.y, 1.0f);
            }
        }
        *uv_undist_4 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
        _S295 = true;
        break;
    }
    return _S295;
}

inline __device__ bool unproject_point_none(float2  uv_16, int camera_model_2, FixedArray<float, 1>  dist_coeffs_13, float3  * raydir_1)
{
    bool _S299;
    for(;;)
    {
        int3  _S300 = make_int3 (int(0));
        float3  _S301 = make_float3 ((float)_S300.x, (float)_S300.y, (float)_S300.z);
        *raydir_1 = _S301;
        if(camera_model_2 == int(3))
        {
            float lon_1 = uv_16.x;
            float lat_1 = uv_16.y;
            float cl_1 = (F32_cos((lat_1)));
            *raydir_1 = make_float3 (cl_1 * (F32_sin((lon_1))), (F32_sin((lat_1))), cl_1 * (F32_cos((lon_1))));
            _S299 = true;
            break;
        }
        FixedArray<float, 1>  _S302 = dist_coeffs_13;
        float2  uv_u_1;
        bool _S303 = undistort_point_0(uv_16, &_S302, int(8), &uv_u_1);
        if(!_S303)
        {
            _S299 = false;
            break;
        }
        float2  _S304 = uv_u_1;
        if(camera_model_2 == int(1))
        {
            float r_12 = length_0(_S304);
            float s_1;
            if(r_12 < 0.00100000004749745f)
            {
                s_1 = 1.0f - r_12 * r_12 / 6.0f;
            }
            else
            {
                s_1 = (F32_sin((r_12))) / r_12;
            }
            *raydir_1 = make_float3 ((_S304 * make_float2 (s_1)).x, (_S304 * make_float2 (s_1)).y, (F32_cos((r_12))));
        }
        else
        {
            if(camera_model_2 == int(2))
            {
                float r_13 = length_0(_S304);
                *raydir_1 = make_float3 ((_S304 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_13 * r_13)))))))).x, (_S304 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_13 * r_13)))))))).y, 1.0f - 0.5f * r_13 * r_13);
            }
            else
            {
                *raydir_1 = make_float3 (_S304.x, _S304.y, 1.0f);
            }
        }
        _S299 = true;
        break;
    }
    return _S299;
}

inline __device__ float3  normalize_0(float3  x_11)
{
    return x_11 / make_float3 (length_1(x_11));
}

inline __device__ bool generate_ray_none(float2  uv_17, int camera_model_3, FixedArray<float, 1>  dist_coeffs_14, float3  * raydir_2)
{
    bool _S305;
    for(;;)
    {
        if(camera_model_3 == int(3))
        {
            float _S306 = uv_17.x;
            if((F32_abs((_S306))) > 3.14159274101257324f)
            {
                _S305 = true;
            }
            else
            {
                _S305 = (F32_abs((uv_17.y))) > 1.57079637050628662f;
            }
            if(_S305)
            {
                int3  _S307 = make_int3 (int(0));
                float3  _S308 = make_float3 ((float)_S307.x, (float)_S307.y, (float)_S307.z);
                *raydir_2 = _S308;
                _S305 = false;
                break;
            }
            float lat_2 = uv_17.y;
            float cl_2 = (F32_cos((lat_2)));
            *raydir_2 = make_float3 (cl_2 * (F32_sin((_S306))), (F32_sin((lat_2))), cl_2 * (F32_cos((_S306))));
            _S305 = true;
            break;
        }
        FixedArray<float, 1>  _S309 = dist_coeffs_14;
        float2  uv_u_2;
        bool _S310 = undistort_point_0(uv_17, &_S309, int(8), &uv_u_2);
        if(!_S310)
        {
            int3  _S311 = make_int3 (int(0));
            float3  _S312 = make_float3 ((float)_S311.x, (float)_S311.y, (float)_S311.z);
            *raydir_2 = _S312;
            _S305 = false;
            break;
        }
        float2  _S313 = uv_u_2;
        if(camera_model_3 == int(1))
        {
            float r_14 = length_0(_S313);
            if(r_14 >= 3.14159274101257324f)
            {
                int3  _S314 = make_int3 (int(0));
                float3  _S315 = make_float3 ((float)_S314.x, (float)_S314.y, (float)_S314.z);
                *raydir_2 = _S315;
                _S305 = false;
                break;
            }
            float s_2;
            if(r_14 < 0.00100000004749745f)
            {
                s_2 = 1.0f - r_14 * r_14 / 6.0f;
            }
            else
            {
                s_2 = (F32_sin((r_14))) / r_14;
            }
            *raydir_2 = make_float3 ((_S313 * make_float2 (s_2)).x, (_S313 * make_float2 (s_2)).y, (F32_cos((r_14))));
        }
        else
        {
            if(camera_model_3 == int(2))
            {
                float r_15 = length_0(_S313);
                if(r_15 >= 2.0f)
                {
                    int3  _S316 = make_int3 (int(0));
                    float3  _S317 = make_float3 ((float)_S316.x, (float)_S316.y, (float)_S316.z);
                    *raydir_2 = _S317;
                    _S305 = false;
                    break;
                }
                *raydir_2 = make_float3 ((_S313 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_15 * r_15)))))))).x, (_S313 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_15 * r_15)))))))).y, 1.0f - 0.5f * r_15 * r_15);
            }
            else
            {
                *raydir_2 = make_float3 (_S313.x, _S313.y, 1.0f);
            }
        }
        *raydir_2 = normalize_0(*raydir_2);
        _S305 = true;
        break;
    }
    return _S305;
}

inline __device__ bool is_valid_distortion_opencv(float2  uv_18, FixedArray<float, 4>  dist_coeffs_15)
{
    float2  _S318 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S319;
    (&_S319)->primal_0 = uv_18;
    (&_S319)->differential_0 = _S318;
    FixedArray<float, 4>  _S320 = dist_coeffs_15;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S321 = s_fwd_DistOpenCV_distort_0(&_S319, &_S320);
    float2  _S322 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S323;
    (&_S323)->primal_0 = uv_18;
    (&_S323)->differential_0 = _S322;
    FixedArray<float, 4>  _S324 = dist_coeffs_15;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S325 = s_fwd_DistOpenCV_distort_0(&_S323, &_S324);
    Matrix<float, 2, 2>  _S326 = transpose_1(makeMatrix<float, 2, 2> (_S321.differential_0, _S325.differential_0));
    float _S327 = (F32_min((determinant_0(_S326)), ((F32_min((_S326.rows[int(0)].x), (_S326.rows[int(1)].y))))));
    bool _S328;
    if(_S327 > 0.25f)
    {
        _S328 = _S327 < 4.0f;
    }
    else
    {
        _S328 = false;
    }
    if(_S328)
    {
        FixedArray<float, 4>  _S329 = dist_coeffs_15;
        float2  _S330 = DistOpenCV_distort_0(uv_18, &_S329);
        _S328 = (dot_0(uv_18, _S330)) >= 0.0f;
    }
    else
    {
        _S328 = false;
    }
    return _S328;
}

inline __device__ bool persp_proj_nav_opencv(float3  p_view_8, float4  intrins_8, FixedArray<float, 4>  dist_coeffs_16, float2  * uv_19)
{
    bool _S331;
    for(;;)
    {
        float2  _S332 = float2 {p_view_8.x, p_view_8.y};
        float _S333 = p_view_8.z;
        float2  uv0_1 = _S332 / make_float2 (_S333);
        if(_S333 < 0.0f)
        {
            _S331 = true;
        }
        else
        {
            float2  _S334 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S335;
            (&_S335)->primal_0 = uv0_1;
            (&_S335)->differential_0 = _S334;
            FixedArray<float, 4>  _S336 = dist_coeffs_16;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S337 = s_fwd_DistOpenCV_distort_0(&_S335, &_S336);
            float2  _S338 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S339;
            (&_S339)->primal_0 = uv0_1;
            (&_S339)->differential_0 = _S338;
            FixedArray<float, 4>  _S340 = dist_coeffs_16;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S341 = s_fwd_DistOpenCV_distort_0(&_S339, &_S340);
            Matrix<float, 2, 2>  _S342 = transpose_1(makeMatrix<float, 2, 2> (_S337.differential_0, _S341.differential_0));
            float _S343 = (F32_min((determinant_0(_S342)), ((F32_min((_S342.rows[int(0)].x), (_S342.rows[int(1)].y))))));
            if(_S343 > 0.25f)
            {
                _S331 = _S343 < 4.0f;
            }
            else
            {
                _S331 = false;
            }
            if(_S331)
            {
                FixedArray<float, 4>  _S344 = dist_coeffs_16;
                float2  _S345 = DistOpenCV_distort_0(uv0_1, &_S344);
                _S331 = (dot_0(uv0_1, _S345)) >= 0.0f;
            }
            else
            {
                _S331 = false;
            }
            _S331 = !_S331;
        }
        if(_S331)
        {
            *uv_19 = uv0_1;
            _S331 = false;
            break;
        }
        float2  uv_20 = _S332 / make_float2 (_S333);
        FixedArray<float, 4>  _S346 = dist_coeffs_16;
        float2  _S347 = DistOpenCV_distort_0(uv_20, &_S346);
        *uv_19 = make_float2 (intrins_8.x * _S347.x + intrins_8.z, intrins_8.y * _S347.y + intrins_8.w);
        _S331 = true;
        break;
    }
    return _S331;
}

inline __device__ bool fisheye_proj_nav_opencv(float3  p_view_9, float4  intrins_9, FixedArray<float, 4>  dist_coeffs_17, float2  * uv_21)
{
    bool _S348;
    for(;;)
    {
        float2  _S349 = float2 {p_view_9.x, p_view_9.y};
        float r_16 = length_0(_S349);
        float _S350 = p_view_9.z;
        float theta_4 = (F32_atan2((r_16), (_S350)));
        bool _S351 = theta_4 < 0.00100000004749745f;
        float k_5;
        if(_S351)
        {
            k_5 = (1.0f - theta_4 * theta_4 / 3.0f) / _S350;
        }
        else
        {
            k_5 = theta_4 / r_16;
        }
        float2  _S352 = _S349 * make_float2 (k_5);
        float2  _S353 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S354;
        (&_S354)->primal_0 = _S352;
        (&_S354)->differential_0 = _S353;
        FixedArray<float, 4>  _S355 = dist_coeffs_17;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S356 = s_fwd_DistOpenCV_distort_0(&_S354, &_S355);
        float2  _S357 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S358;
        (&_S358)->primal_0 = _S352;
        (&_S358)->differential_0 = _S357;
        FixedArray<float, 4>  _S359 = dist_coeffs_17;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S360 = s_fwd_DistOpenCV_distort_0(&_S358, &_S359);
        Matrix<float, 2, 2>  _S361 = transpose_1(makeMatrix<float, 2, 2> (_S356.differential_0, _S360.differential_0));
        float _S362 = (F32_min((determinant_0(_S361)), ((F32_min((_S361.rows[int(0)].x), (_S361.rows[int(1)].y))))));
        if(_S362 > 0.25f)
        {
            _S348 = _S362 < 4.0f;
        }
        else
        {
            _S348 = false;
        }
        if(_S348)
        {
            FixedArray<float, 4>  _S363 = dist_coeffs_17;
            float2  _S364 = DistOpenCV_distort_0(_S352, &_S363);
            _S348 = (dot_0(_S352, _S364)) >= 0.0f;
        }
        else
        {
            _S348 = false;
        }
        if(!_S348)
        {
            *uv_21 = _S352;
            _S348 = false;
            break;
        }
        if(_S351)
        {
            k_5 = (1.0f - theta_4 * theta_4 / 3.0f) / _S350;
        }
        else
        {
            k_5 = theta_4 / r_16;
        }
        float2  _S365 = _S349 * make_float2 (k_5);
        FixedArray<float, 4>  _S366 = dist_coeffs_17;
        float2  _S367 = DistOpenCV_distort_0(_S365, &_S366);
        *uv_21 = make_float2 (intrins_9.x * _S367.x + intrins_9.z, intrins_9.y * _S367.y + intrins_9.w);
        _S348 = true;
        break;
    }
    return _S348;
}

inline __device__ bool equisolid_proj_nav_opencv(float3  p_view_10, float4  intrins_10, FixedArray<float, 4>  dist_coeffs_18, float2  * uv_22)
{
    bool _S368;
    for(;;)
    {
        float2  _S369 = float2 {p_view_10.x, p_view_10.y};
        float r_17 = length_0(_S369);
        float _S370 = p_view_10.z;
        float theta_5 = (F32_atan2((r_17), (_S370)));
        bool _S371 = r_17 < 9.99999997475242708e-07f;
        float k_6;
        if(_S371)
        {
            k_6 = (1.0f - theta_5 * theta_5 / 24.0f) / _S370;
        }
        else
        {
            k_6 = 2.0f * (F32_sin((0.5f * theta_5))) / r_17;
        }
        float2  _S372 = _S369 * make_float2 (k_6);
        float2  _S373 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S374;
        (&_S374)->primal_0 = _S372;
        (&_S374)->differential_0 = _S373;
        FixedArray<float, 4>  _S375 = dist_coeffs_18;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S376 = s_fwd_DistOpenCV_distort_0(&_S374, &_S375);
        float2  _S377 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S378;
        (&_S378)->primal_0 = _S372;
        (&_S378)->differential_0 = _S377;
        FixedArray<float, 4>  _S379 = dist_coeffs_18;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S380 = s_fwd_DistOpenCV_distort_0(&_S378, &_S379);
        Matrix<float, 2, 2>  _S381 = transpose_1(makeMatrix<float, 2, 2> (_S376.differential_0, _S380.differential_0));
        float _S382 = (F32_min((determinant_0(_S381)), ((F32_min((_S381.rows[int(0)].x), (_S381.rows[int(1)].y))))));
        if(_S382 > 0.25f)
        {
            _S368 = _S382 < 4.0f;
        }
        else
        {
            _S368 = false;
        }
        if(_S368)
        {
            FixedArray<float, 4>  _S383 = dist_coeffs_18;
            float2  _S384 = DistOpenCV_distort_0(_S372, &_S383);
            _S368 = (dot_0(_S372, _S384)) >= 0.0f;
        }
        else
        {
            _S368 = false;
        }
        if(!_S368)
        {
            *uv_22 = _S372;
            _S368 = false;
            break;
        }
        if(_S371)
        {
            k_6 = (1.0f - theta_5 * theta_5 / 24.0f) / _S370;
        }
        else
        {
            k_6 = 2.0f * (F32_sin((0.5f * theta_5))) / r_17;
        }
        float2  _S385 = _S369 * make_float2 (k_6);
        FixedArray<float, 4>  _S386 = dist_coeffs_18;
        float2  _S387 = DistOpenCV_distort_0(_S385, &_S386);
        *uv_22 = make_float2 (intrins_10.x * _S387.x + intrins_10.z, intrins_10.y * _S387.y + intrins_10.w);
        _S368 = true;
        break;
    }
    return _S368;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac_opencv(float3  p_view_11, float4  intrins_11, FixedArray<float, 4>  dist_coeffs_19)
{
    float2  _S388 = float2 {p_view_11.x, p_view_11.y};
    float _S389 = p_view_11.z;
    float2  _S390 = _S388 * make_float2 (0.0f);
    float _S391 = _S389 * _S389;
    float2  s_diff_uv_3 = (make_float2 (1.0f, 0.0f) * make_float2 (_S389) - _S390) / make_float2 (_S391);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S392;
    (&_S392)->primal_0 = _S388 / make_float2 (_S389);
    (&_S392)->differential_0 = s_diff_uv_3;
    FixedArray<float, 4>  _S393 = dist_coeffs_19;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S394 = s_fwd_DistOpenCV_distort_0(&_S392, &_S393);
    float fx_4 = intrins_11.x;
    float fy_4 = intrins_11.y;
    float _S395 = _S394.differential_0.y * fy_4;
    Matrix<float, 2, 3>  J_4;
    *&(((&J_4)->rows + (int(0)))->x) = _S394.differential_0.x * fx_4;
    *&(((&J_4)->rows + (int(1)))->x) = _S395;
    float2  s_diff_uv_4 = (make_float2 (0.0f, 1.0f) * make_float2 (_S389) - _S390) / make_float2 (_S391);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S396;
    (&_S396)->primal_0 = _S388 / make_float2 (_S389);
    (&_S396)->differential_0 = s_diff_uv_4;
    FixedArray<float, 4>  _S397 = dist_coeffs_19;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S398 = s_fwd_DistOpenCV_distort_0(&_S396, &_S397);
    float _S399 = _S398.differential_0.y * fy_4;
    *&(((&J_4)->rows + (int(0)))->y) = _S398.differential_0.x * fx_4;
    *&(((&J_4)->rows + (int(1)))->y) = _S399;
    float2  s_diff_uv_5 = (make_float2 (0.0f, 0.0f) * make_float2 (_S389) - _S388) / make_float2 (_S391);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S400;
    (&_S400)->primal_0 = _S388 / make_float2 (_S389);
    (&_S400)->differential_0 = s_diff_uv_5;
    FixedArray<float, 4>  _S401 = dist_coeffs_19;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S402 = s_fwd_DistOpenCV_distort_0(&_S400, &_S401);
    float _S403 = _S402.differential_0.y * fy_4;
    *&(((&J_4)->rows + (int(0)))->z) = _S402.differential_0.x * fx_4;
    *&(((&J_4)->rows + (int(1)))->z) = _S403;
    return J_4;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac_opencv(float3  p_view_12, float4  intrins_12, FixedArray<float, 4>  dist_coeffs_20)
{
    Matrix<float, 2, 3>  J_5;
    float2  _S404 = float2 {p_view_12.x, p_view_12.y};
    float2  _S405 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S406;
    (&_S406)->primal_0 = _S404;
    (&_S406)->differential_0 = _S405;
    DiffPair_float_0 _S407 = s_fwd_length_impl_0(&_S406);
    float _S408 = p_view_12.z;
    DiffPair_float_0 _S409;
    (&_S409)->primal_0 = _S407.primal_0;
    (&_S409)->differential_0 = _S407.differential_0;
    DiffPair_float_0 _S410;
    (&_S410)->primal_0 = _S408;
    (&_S410)->differential_0 = 0.0f;
    DiffPair_float_0 _S411 = _d_atan2_0(&_S409, &_S410);
    float k_7;
    float s_diff_k_2;
    if((_S411.primal_0) < 0.00100000004749745f)
    {
        float _S412 = _S411.differential_0 * _S411.primal_0;
        float _S413 = 1.0f - _S411.primal_0 * _S411.primal_0 / 3.0f;
        float _S414 = ((0.0f - (_S412 + _S412) * 0.3333333432674408f) * _S408 - _S413 * 0.0f) / (_S408 * _S408);
        k_7 = _S413 / _S408;
        s_diff_k_2 = _S414;
    }
    else
    {
        float _S415 = (_S411.differential_0 * _S407.primal_0 - _S411.primal_0 * _S407.differential_0) / (_S407.primal_0 * _S407.primal_0);
        k_7 = _S411.primal_0 / _S407.primal_0;
        s_diff_k_2 = _S415;
    }
    float2  _S416 = _S404 * make_float2 (k_7);
    float2  _S417 = _S405 * make_float2 (k_7) + make_float2 (s_diff_k_2) * _S404;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S418;
    (&_S418)->primal_0 = _S416;
    (&_S418)->differential_0 = _S417;
    FixedArray<float, 4>  _S419 = dist_coeffs_20;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S420 = s_fwd_DistOpenCV_distort_0(&_S418, &_S419);
    float fx_5 = intrins_12.x;
    float fy_5 = intrins_12.y;
    float _S421 = _S420.differential_0.y * fy_5;
    *&(((&J_5)->rows + (int(0)))->x) = _S420.differential_0.x * fx_5;
    *&(((&J_5)->rows + (int(1)))->x) = _S421;
    float2  _S422 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S423;
    (&_S423)->primal_0 = _S404;
    (&_S423)->differential_0 = _S422;
    DiffPair_float_0 _S424 = s_fwd_length_impl_0(&_S423);
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = _S424.primal_0;
    (&_S425)->differential_0 = _S424.differential_0;
    DiffPair_float_0 _S426;
    (&_S426)->primal_0 = _S408;
    (&_S426)->differential_0 = 0.0f;
    DiffPair_float_0 _S427 = _d_atan2_0(&_S425, &_S426);
    if((_S427.primal_0) < 0.00100000004749745f)
    {
        float _S428 = _S427.differential_0 * _S427.primal_0;
        float _S429 = 1.0f - _S427.primal_0 * _S427.primal_0 / 3.0f;
        float _S430 = ((0.0f - (_S428 + _S428) * 0.3333333432674408f) * _S408 - _S429 * 0.0f) / (_S408 * _S408);
        k_7 = _S429 / _S408;
        s_diff_k_2 = _S430;
    }
    else
    {
        float _S431 = (_S427.differential_0 * _S424.primal_0 - _S427.primal_0 * _S424.differential_0) / (_S424.primal_0 * _S424.primal_0);
        k_7 = _S427.primal_0 / _S424.primal_0;
        s_diff_k_2 = _S431;
    }
    float2  _S432 = _S404 * make_float2 (k_7);
    float2  _S433 = _S422 * make_float2 (k_7) + make_float2 (s_diff_k_2) * _S404;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S434;
    (&_S434)->primal_0 = _S432;
    (&_S434)->differential_0 = _S433;
    FixedArray<float, 4>  _S435 = dist_coeffs_20;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S436 = s_fwd_DistOpenCV_distort_0(&_S434, &_S435);
    float _S437 = _S436.differential_0.y * fy_5;
    *&(((&J_5)->rows + (int(0)))->y) = _S436.differential_0.x * fx_5;
    *&(((&J_5)->rows + (int(1)))->y) = _S437;
    float2  _S438 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S439;
    (&_S439)->primal_0 = _S404;
    (&_S439)->differential_0 = _S438;
    DiffPair_float_0 _S440 = s_fwd_length_impl_0(&_S439);
    DiffPair_float_0 _S441;
    (&_S441)->primal_0 = _S440.primal_0;
    (&_S441)->differential_0 = _S440.differential_0;
    DiffPair_float_0 _S442;
    (&_S442)->primal_0 = _S408;
    (&_S442)->differential_0 = 1.0f;
    DiffPair_float_0 _S443 = _d_atan2_0(&_S441, &_S442);
    if((_S443.primal_0) < 0.00100000004749745f)
    {
        float _S444 = _S443.differential_0 * _S443.primal_0;
        float _S445 = 1.0f - _S443.primal_0 * _S443.primal_0 / 3.0f;
        float _S446 = ((0.0f - (_S444 + _S444) * 0.3333333432674408f) * _S408 - _S445) / (_S408 * _S408);
        k_7 = _S445 / _S408;
        s_diff_k_2 = _S446;
    }
    else
    {
        float _S447 = (_S443.differential_0 * _S440.primal_0 - _S443.primal_0 * _S440.differential_0) / (_S440.primal_0 * _S440.primal_0);
        k_7 = _S443.primal_0 / _S440.primal_0;
        s_diff_k_2 = _S447;
    }
    float2  _S448 = _S404 * make_float2 (k_7);
    float2  _S449 = _S438 * make_float2 (k_7) + make_float2 (s_diff_k_2) * _S404;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S450;
    (&_S450)->primal_0 = _S448;
    (&_S450)->differential_0 = _S449;
    FixedArray<float, 4>  _S451 = dist_coeffs_20;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S452 = s_fwd_DistOpenCV_distort_0(&_S450, &_S451);
    float _S453 = _S452.differential_0.y * fy_5;
    *&(((&J_5)->rows + (int(0)))->z) = _S452.differential_0.x * fx_5;
    *&(((&J_5)->rows + (int(1)))->z) = _S453;
    return J_5;
}

inline __device__ Matrix<float, 2, 3>  equisolid_proj_jac_opencv(float3  p_view_13, float4  intrins_13, FixedArray<float, 4>  dist_coeffs_21)
{
    Matrix<float, 2, 3>  J_6;
    float2  _S454 = float2 {p_view_13.x, p_view_13.y};
    float2  _S455 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S456;
    (&_S456)->primal_0 = _S454;
    (&_S456)->differential_0 = _S455;
    DiffPair_float_0 _S457 = s_fwd_length_impl_0(&_S456);
    float _S458 = p_view_13.z;
    DiffPair_float_0 _S459;
    (&_S459)->primal_0 = _S457.primal_0;
    (&_S459)->differential_0 = _S457.differential_0;
    DiffPair_float_0 _S460;
    (&_S460)->primal_0 = _S458;
    (&_S460)->differential_0 = 0.0f;
    DiffPair_float_0 _S461 = _d_atan2_0(&_S459, &_S460);
    float k_8;
    float s_diff_k_3;
    if((_S457.primal_0) < 9.99999997475242708e-07f)
    {
        float _S462 = _S461.differential_0 * _S461.primal_0;
        float _S463 = 1.0f - _S461.primal_0 * _S461.primal_0 / 24.0f;
        float _S464 = ((0.0f - (_S462 + _S462) * 0.0416666679084301f) * _S458 - _S463 * 0.0f) / (_S458 * _S458);
        k_8 = _S463 / _S458;
        s_diff_k_3 = _S464;
    }
    else
    {
        float _S465 = _S461.differential_0 * 0.5f;
        DiffPair_float_0 _S466;
        (&_S466)->primal_0 = 0.5f * _S461.primal_0;
        (&_S466)->differential_0 = _S465;
        DiffPair_float_0 _S467 = _d_sin_0(&_S466);
        float _S468 = 2.0f * _S467.primal_0;
        float _S469 = (_S467.differential_0 * 2.0f * _S457.primal_0 - _S468 * _S457.differential_0) / (_S457.primal_0 * _S457.primal_0);
        k_8 = _S468 / _S457.primal_0;
        s_diff_k_3 = _S469;
    }
    float2  _S470 = _S454 * make_float2 (k_8);
    float2  _S471 = _S455 * make_float2 (k_8) + make_float2 (s_diff_k_3) * _S454;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S472;
    (&_S472)->primal_0 = _S470;
    (&_S472)->differential_0 = _S471;
    FixedArray<float, 4>  _S473 = dist_coeffs_21;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S474 = s_fwd_DistOpenCV_distort_0(&_S472, &_S473);
    float fx_6 = intrins_13.x;
    float fy_6 = intrins_13.y;
    float _S475 = _S474.differential_0.y * fy_6;
    *&(((&J_6)->rows + (int(0)))->x) = _S474.differential_0.x * fx_6;
    *&(((&J_6)->rows + (int(1)))->x) = _S475;
    float2  _S476 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S477;
    (&_S477)->primal_0 = _S454;
    (&_S477)->differential_0 = _S476;
    DiffPair_float_0 _S478 = s_fwd_length_impl_0(&_S477);
    DiffPair_float_0 _S479;
    (&_S479)->primal_0 = _S478.primal_0;
    (&_S479)->differential_0 = _S478.differential_0;
    DiffPair_float_0 _S480;
    (&_S480)->primal_0 = _S458;
    (&_S480)->differential_0 = 0.0f;
    DiffPair_float_0 _S481 = _d_atan2_0(&_S479, &_S480);
    if((_S478.primal_0) < 9.99999997475242708e-07f)
    {
        float _S482 = _S481.differential_0 * _S481.primal_0;
        float _S483 = 1.0f - _S481.primal_0 * _S481.primal_0 / 24.0f;
        float _S484 = ((0.0f - (_S482 + _S482) * 0.0416666679084301f) * _S458 - _S483 * 0.0f) / (_S458 * _S458);
        k_8 = _S483 / _S458;
        s_diff_k_3 = _S484;
    }
    else
    {
        float _S485 = _S481.differential_0 * 0.5f;
        DiffPair_float_0 _S486;
        (&_S486)->primal_0 = 0.5f * _S481.primal_0;
        (&_S486)->differential_0 = _S485;
        DiffPair_float_0 _S487 = _d_sin_0(&_S486);
        float _S488 = 2.0f * _S487.primal_0;
        float _S489 = (_S487.differential_0 * 2.0f * _S478.primal_0 - _S488 * _S478.differential_0) / (_S478.primal_0 * _S478.primal_0);
        k_8 = _S488 / _S478.primal_0;
        s_diff_k_3 = _S489;
    }
    float2  _S490 = _S454 * make_float2 (k_8);
    float2  _S491 = _S476 * make_float2 (k_8) + make_float2 (s_diff_k_3) * _S454;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S492;
    (&_S492)->primal_0 = _S490;
    (&_S492)->differential_0 = _S491;
    FixedArray<float, 4>  _S493 = dist_coeffs_21;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S494 = s_fwd_DistOpenCV_distort_0(&_S492, &_S493);
    float _S495 = _S494.differential_0.y * fy_6;
    *&(((&J_6)->rows + (int(0)))->y) = _S494.differential_0.x * fx_6;
    *&(((&J_6)->rows + (int(1)))->y) = _S495;
    float2  _S496 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S497;
    (&_S497)->primal_0 = _S454;
    (&_S497)->differential_0 = _S496;
    DiffPair_float_0 _S498 = s_fwd_length_impl_0(&_S497);
    DiffPair_float_0 _S499;
    (&_S499)->primal_0 = _S498.primal_0;
    (&_S499)->differential_0 = _S498.differential_0;
    DiffPair_float_0 _S500;
    (&_S500)->primal_0 = _S458;
    (&_S500)->differential_0 = 1.0f;
    DiffPair_float_0 _S501 = _d_atan2_0(&_S499, &_S500);
    if((_S498.primal_0) < 9.99999997475242708e-07f)
    {
        float _S502 = _S501.differential_0 * _S501.primal_0;
        float _S503 = 1.0f - _S501.primal_0 * _S501.primal_0 / 24.0f;
        float _S504 = ((0.0f - (_S502 + _S502) * 0.0416666679084301f) * _S458 - _S503) / (_S458 * _S458);
        k_8 = _S503 / _S458;
        s_diff_k_3 = _S504;
    }
    else
    {
        float _S505 = _S501.differential_0 * 0.5f;
        DiffPair_float_0 _S506;
        (&_S506)->primal_0 = 0.5f * _S501.primal_0;
        (&_S506)->differential_0 = _S505;
        DiffPair_float_0 _S507 = _d_sin_0(&_S506);
        float _S508 = 2.0f * _S507.primal_0;
        float _S509 = (_S507.differential_0 * 2.0f * _S498.primal_0 - _S508 * _S498.differential_0) / (_S498.primal_0 * _S498.primal_0);
        k_8 = _S508 / _S498.primal_0;
        s_diff_k_3 = _S509;
    }
    float2  _S510 = _S454 * make_float2 (k_8);
    float2  _S511 = _S496 * make_float2 (k_8) + make_float2 (s_diff_k_3) * _S454;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S512;
    (&_S512)->primal_0 = _S510;
    (&_S512)->differential_0 = _S511;
    FixedArray<float, 4>  _S513 = dist_coeffs_21;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S514 = s_fwd_DistOpenCV_distort_0(&_S512, &_S513);
    float _S515 = _S514.differential_0.y * fy_6;
    *&(((&J_6)->rows + (int(0)))->z) = _S514.differential_0.x * fx_6;
    *&(((&J_6)->rows + (int(1)))->z) = _S515;
    return J_6;
}

inline __device__ float2  distort_point_opencv(float2  uv_23, int camera_model_4, FixedArray<float, 4>  dist_coeffs_22)
{
    float2  _S516;
    for(;;)
    {
        if(camera_model_4 == int(3))
        {
            _S516 = uv_23;
            break;
        }
        float k_9;
        if(camera_model_4 == int(1))
        {
            float r_18 = length_0(uv_23);
            float theta_6 = (F32_atan((r_18)));
            if(r_18 < 0.00100000004749745f)
            {
                k_9 = 1.0f - theta_6 * theta_6 / 6.0f;
            }
            else
            {
                k_9 = theta_6 / r_18;
            }
            _S516 = uv_23 * make_float2 (k_9);
        }
        else
        {
            if(camera_model_4 == int(2))
            {
                float r_19 = length_0(uv_23);
                float theta_7 = (F32_atan((r_19)));
                if(r_19 < 0.00100000004749745f)
                {
                    k_9 = 1.0f - theta_7 * theta_7 / 24.0f;
                }
                else
                {
                    k_9 = 2.0f * (F32_sin((0.5f * theta_7))) / r_19;
                }
                _S516 = uv_23 * make_float2 (k_9);
            }
            else
            {
                _S516 = uv_23;
            }
        }
        FixedArray<float, 4>  _S517 = dist_coeffs_22;
        float2  _S518 = DistOpenCV_distort_0(_S516, &_S517);
        _S516 = _S518;
        break;
    }
    return _S516;
}

inline __device__ bool undistort_point_opencv(float2  uv_24, int camera_model_5, FixedArray<float, 4>  dist_coeffs_23, float2  * uv_undist_5)
{
    bool _S519;
    for(;;)
    {
        *uv_undist_5 = make_float2 (0.0f);
        if(camera_model_5 == int(3))
        {
            float lon_2 = uv_24.x;
            float lat_3 = uv_24.y;
            float cl_3 = (F32_cos((lat_3)));
            *uv_undist_5 = make_float2 (cl_3 * (F32_sin((lon_2))), (F32_sin((lat_3)))) / make_float2 ((F32_max((cl_3 * (F32_cos((lon_2)))), (9.999999960041972e-13f))));
            _S519 = true;
            break;
        }
        FixedArray<float, 4>  _S520 = dist_coeffs_23;
        float2  uv_u_3;
        bool _S521 = undistort_point_1(uv_24, &_S520, int(8), &uv_u_3);
        if(!_S521)
        {
            _S519 = false;
            break;
        }
        float2  _S522 = uv_u_3;
        float3  raydir_3;
        if(camera_model_5 == int(1))
        {
            float r_20 = length_0(_S522);
            float s_3;
            if(r_20 < 0.00100000004749745f)
            {
                s_3 = 1.0f - r_20 * r_20 / 6.0f;
            }
            else
            {
                s_3 = (F32_sin((r_20))) / r_20;
            }
            raydir_3 = make_float3 ((_S522 * make_float2 (s_3)).x, (_S522 * make_float2 (s_3)).y, (F32_cos((r_20))));
        }
        else
        {
            if(camera_model_5 == int(2))
            {
                float r_21 = length_0(_S522);
                raydir_3 = make_float3 ((_S522 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_21 * r_21)))))))).x, (_S522 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_21 * r_21)))))))).y, 1.0f - 0.5f * r_21 * r_21);
            }
            else
            {
                raydir_3 = make_float3 (_S522.x, _S522.y, 1.0f);
            }
        }
        *uv_undist_5 = float2 {raydir_3.x, raydir_3.y} / make_float2 ((F32_max((raydir_3.z), (9.999999960041972e-13f))));
        _S519 = true;
        break;
    }
    return _S519;
}

inline __device__ bool unproject_point_opencv(float2  uv_25, int camera_model_6, FixedArray<float, 4>  dist_coeffs_24, float3  * raydir_4)
{
    bool _S523;
    for(;;)
    {
        int3  _S524 = make_int3 (int(0));
        float3  _S525 = make_float3 ((float)_S524.x, (float)_S524.y, (float)_S524.z);
        *raydir_4 = _S525;
        if(camera_model_6 == int(3))
        {
            float lon_3 = uv_25.x;
            float lat_4 = uv_25.y;
            float cl_4 = (F32_cos((lat_4)));
            *raydir_4 = make_float3 (cl_4 * (F32_sin((lon_3))), (F32_sin((lat_4))), cl_4 * (F32_cos((lon_3))));
            _S523 = true;
            break;
        }
        FixedArray<float, 4>  _S526 = dist_coeffs_24;
        float2  uv_u_4;
        bool _S527 = undistort_point_1(uv_25, &_S526, int(8), &uv_u_4);
        if(!_S527)
        {
            _S523 = false;
            break;
        }
        float2  _S528 = uv_u_4;
        if(camera_model_6 == int(1))
        {
            float r_22 = length_0(_S528);
            float s_4;
            if(r_22 < 0.00100000004749745f)
            {
                s_4 = 1.0f - r_22 * r_22 / 6.0f;
            }
            else
            {
                s_4 = (F32_sin((r_22))) / r_22;
            }
            *raydir_4 = make_float3 ((_S528 * make_float2 (s_4)).x, (_S528 * make_float2 (s_4)).y, (F32_cos((r_22))));
        }
        else
        {
            if(camera_model_6 == int(2))
            {
                float r_23 = length_0(_S528);
                *raydir_4 = make_float3 ((_S528 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_23 * r_23)))))))).x, (_S528 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_23 * r_23)))))))).y, 1.0f - 0.5f * r_23 * r_23);
            }
            else
            {
                *raydir_4 = make_float3 (_S528.x, _S528.y, 1.0f);
            }
        }
        _S523 = true;
        break;
    }
    return _S523;
}

inline __device__ bool generate_ray_opencv(float2  uv_26, int camera_model_7, FixedArray<float, 4>  dist_coeffs_25, float3  * raydir_5)
{
    bool _S529;
    for(;;)
    {
        if(camera_model_7 == int(3))
        {
            float _S530 = uv_26.x;
            if((F32_abs((_S530))) > 3.14159274101257324f)
            {
                _S529 = true;
            }
            else
            {
                _S529 = (F32_abs((uv_26.y))) > 1.57079637050628662f;
            }
            if(_S529)
            {
                int3  _S531 = make_int3 (int(0));
                float3  _S532 = make_float3 ((float)_S531.x, (float)_S531.y, (float)_S531.z);
                *raydir_5 = _S532;
                _S529 = false;
                break;
            }
            float lat_5 = uv_26.y;
            float cl_5 = (F32_cos((lat_5)));
            *raydir_5 = make_float3 (cl_5 * (F32_sin((_S530))), (F32_sin((lat_5))), cl_5 * (F32_cos((_S530))));
            _S529 = true;
            break;
        }
        FixedArray<float, 4>  _S533 = dist_coeffs_25;
        float2  uv_u_5;
        bool _S534 = undistort_point_1(uv_26, &_S533, int(8), &uv_u_5);
        if(!_S534)
        {
            int3  _S535 = make_int3 (int(0));
            float3  _S536 = make_float3 ((float)_S535.x, (float)_S535.y, (float)_S535.z);
            *raydir_5 = _S536;
            _S529 = false;
            break;
        }
        float2  _S537 = uv_u_5;
        if(camera_model_7 == int(1))
        {
            float r_24 = length_0(_S537);
            if(r_24 >= 3.14159274101257324f)
            {
                int3  _S538 = make_int3 (int(0));
                float3  _S539 = make_float3 ((float)_S538.x, (float)_S538.y, (float)_S538.z);
                *raydir_5 = _S539;
                _S529 = false;
                break;
            }
            float s_5;
            if(r_24 < 0.00100000004749745f)
            {
                s_5 = 1.0f - r_24 * r_24 / 6.0f;
            }
            else
            {
                s_5 = (F32_sin((r_24))) / r_24;
            }
            *raydir_5 = make_float3 ((_S537 * make_float2 (s_5)).x, (_S537 * make_float2 (s_5)).y, (F32_cos((r_24))));
        }
        else
        {
            if(camera_model_7 == int(2))
            {
                float r_25 = length_0(_S537);
                if(r_25 >= 2.0f)
                {
                    int3  _S540 = make_int3 (int(0));
                    float3  _S541 = make_float3 ((float)_S540.x, (float)_S540.y, (float)_S540.z);
                    *raydir_5 = _S541;
                    _S529 = false;
                    break;
                }
                *raydir_5 = make_float3 ((_S537 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_25 * r_25)))))))).x, (_S537 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_25 * r_25)))))))).y, 1.0f - 0.5f * r_25 * r_25);
            }
            else
            {
                *raydir_5 = make_float3 (_S537.x, _S537.y, 1.0f);
            }
        }
        *raydir_5 = normalize_0(*raydir_5);
        _S529 = true;
        break;
    }
    return _S529;
}

inline __device__ bool is_valid_distortion_prism(float2  uv_27, FixedArray<float, 8>  dist_coeffs_26)
{
    float2  _S542 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S543;
    (&_S543)->primal_0 = uv_27;
    (&_S543)->differential_0 = _S542;
    FixedArray<float, 8>  _S544 = dist_coeffs_26;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S545 = s_fwd_DistThinPrism_distort_0(&_S543, &_S544);
    float2  _S546 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S547;
    (&_S547)->primal_0 = uv_27;
    (&_S547)->differential_0 = _S546;
    FixedArray<float, 8>  _S548 = dist_coeffs_26;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S549 = s_fwd_DistThinPrism_distort_0(&_S547, &_S548);
    Matrix<float, 2, 2>  _S550 = transpose_1(makeMatrix<float, 2, 2> (_S545.differential_0, _S549.differential_0));
    float _S551 = (F32_min((determinant_0(_S550)), ((F32_min((_S550.rows[int(0)].x), (_S550.rows[int(1)].y))))));
    bool _S552;
    if(_S551 > 0.25f)
    {
        _S552 = _S551 < 4.0f;
    }
    else
    {
        _S552 = false;
    }
    if(_S552)
    {
        FixedArray<float, 8>  _S553 = dist_coeffs_26;
        float2  _S554 = DistThinPrism_distort_0(uv_27, &_S553);
        _S552 = (dot_0(uv_27, _S554)) >= 0.0f;
    }
    else
    {
        _S552 = false;
    }
    return _S552;
}

inline __device__ bool persp_proj_nav_prism(float3  p_view_14, float4  intrins_14, FixedArray<float, 8>  dist_coeffs_27, float2  * uv_28)
{
    bool _S555;
    for(;;)
    {
        float2  _S556 = float2 {p_view_14.x, p_view_14.y};
        float _S557 = p_view_14.z;
        float2  uv0_2 = _S556 / make_float2 (_S557);
        if(_S557 < 0.0f)
        {
            _S555 = true;
        }
        else
        {
            float2  _S558 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S559;
            (&_S559)->primal_0 = uv0_2;
            (&_S559)->differential_0 = _S558;
            FixedArray<float, 8>  _S560 = dist_coeffs_27;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S561 = s_fwd_DistThinPrism_distort_0(&_S559, &_S560);
            float2  _S562 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S563;
            (&_S563)->primal_0 = uv0_2;
            (&_S563)->differential_0 = _S562;
            FixedArray<float, 8>  _S564 = dist_coeffs_27;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S565 = s_fwd_DistThinPrism_distort_0(&_S563, &_S564);
            Matrix<float, 2, 2>  _S566 = transpose_1(makeMatrix<float, 2, 2> (_S561.differential_0, _S565.differential_0));
            float _S567 = (F32_min((determinant_0(_S566)), ((F32_min((_S566.rows[int(0)].x), (_S566.rows[int(1)].y))))));
            if(_S567 > 0.25f)
            {
                _S555 = _S567 < 4.0f;
            }
            else
            {
                _S555 = false;
            }
            if(_S555)
            {
                FixedArray<float, 8>  _S568 = dist_coeffs_27;
                float2  _S569 = DistThinPrism_distort_0(uv0_2, &_S568);
                _S555 = (dot_0(uv0_2, _S569)) >= 0.0f;
            }
            else
            {
                _S555 = false;
            }
            _S555 = !_S555;
        }
        if(_S555)
        {
            *uv_28 = uv0_2;
            _S555 = false;
            break;
        }
        float2  uv_29 = _S556 / make_float2 (_S557);
        FixedArray<float, 8>  _S570 = dist_coeffs_27;
        float2  _S571 = DistThinPrism_distort_0(uv_29, &_S570);
        *uv_28 = make_float2 (intrins_14.x * _S571.x + intrins_14.z, intrins_14.y * _S571.y + intrins_14.w);
        _S555 = true;
        break;
    }
    return _S555;
}

inline __device__ bool fisheye_proj_nav_prism(float3  p_view_15, float4  intrins_15, FixedArray<float, 8>  dist_coeffs_28, float2  * uv_30)
{
    bool _S572;
    for(;;)
    {
        float2  _S573 = float2 {p_view_15.x, p_view_15.y};
        float r_26 = length_0(_S573);
        float _S574 = p_view_15.z;
        float theta_8 = (F32_atan2((r_26), (_S574)));
        bool _S575 = theta_8 < 0.00100000004749745f;
        float k_10;
        if(_S575)
        {
            k_10 = (1.0f - theta_8 * theta_8 / 3.0f) / _S574;
        }
        else
        {
            k_10 = theta_8 / r_26;
        }
        float2  _S576 = _S573 * make_float2 (k_10);
        float2  _S577 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S578;
        (&_S578)->primal_0 = _S576;
        (&_S578)->differential_0 = _S577;
        FixedArray<float, 8>  _S579 = dist_coeffs_28;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S580 = s_fwd_DistThinPrism_distort_0(&_S578, &_S579);
        float2  _S581 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S582;
        (&_S582)->primal_0 = _S576;
        (&_S582)->differential_0 = _S581;
        FixedArray<float, 8>  _S583 = dist_coeffs_28;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S584 = s_fwd_DistThinPrism_distort_0(&_S582, &_S583);
        Matrix<float, 2, 2>  _S585 = transpose_1(makeMatrix<float, 2, 2> (_S580.differential_0, _S584.differential_0));
        float _S586 = (F32_min((determinant_0(_S585)), ((F32_min((_S585.rows[int(0)].x), (_S585.rows[int(1)].y))))));
        if(_S586 > 0.25f)
        {
            _S572 = _S586 < 4.0f;
        }
        else
        {
            _S572 = false;
        }
        if(_S572)
        {
            FixedArray<float, 8>  _S587 = dist_coeffs_28;
            float2  _S588 = DistThinPrism_distort_0(_S576, &_S587);
            _S572 = (dot_0(_S576, _S588)) >= 0.0f;
        }
        else
        {
            _S572 = false;
        }
        if(!_S572)
        {
            *uv_30 = _S576;
            _S572 = false;
            break;
        }
        if(_S575)
        {
            k_10 = (1.0f - theta_8 * theta_8 / 3.0f) / _S574;
        }
        else
        {
            k_10 = theta_8 / r_26;
        }
        float2  _S589 = _S573 * make_float2 (k_10);
        FixedArray<float, 8>  _S590 = dist_coeffs_28;
        float2  _S591 = DistThinPrism_distort_0(_S589, &_S590);
        *uv_30 = make_float2 (intrins_15.x * _S591.x + intrins_15.z, intrins_15.y * _S591.y + intrins_15.w);
        _S572 = true;
        break;
    }
    return _S572;
}

inline __device__ bool equisolid_proj_nav_prism(float3  p_view_16, float4  intrins_16, FixedArray<float, 8>  dist_coeffs_29, float2  * uv_31)
{
    bool _S592;
    for(;;)
    {
        float2  _S593 = float2 {p_view_16.x, p_view_16.y};
        float r_27 = length_0(_S593);
        float _S594 = p_view_16.z;
        float theta_9 = (F32_atan2((r_27), (_S594)));
        bool _S595 = r_27 < 9.99999997475242708e-07f;
        float k_11;
        if(_S595)
        {
            k_11 = (1.0f - theta_9 * theta_9 / 24.0f) / _S594;
        }
        else
        {
            k_11 = 2.0f * (F32_sin((0.5f * theta_9))) / r_27;
        }
        float2  _S596 = _S593 * make_float2 (k_11);
        float2  _S597 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S598;
        (&_S598)->primal_0 = _S596;
        (&_S598)->differential_0 = _S597;
        FixedArray<float, 8>  _S599 = dist_coeffs_29;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S600 = s_fwd_DistThinPrism_distort_0(&_S598, &_S599);
        float2  _S601 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S602;
        (&_S602)->primal_0 = _S596;
        (&_S602)->differential_0 = _S601;
        FixedArray<float, 8>  _S603 = dist_coeffs_29;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S604 = s_fwd_DistThinPrism_distort_0(&_S602, &_S603);
        Matrix<float, 2, 2>  _S605 = transpose_1(makeMatrix<float, 2, 2> (_S600.differential_0, _S604.differential_0));
        float _S606 = (F32_min((determinant_0(_S605)), ((F32_min((_S605.rows[int(0)].x), (_S605.rows[int(1)].y))))));
        if(_S606 > 0.25f)
        {
            _S592 = _S606 < 4.0f;
        }
        else
        {
            _S592 = false;
        }
        if(_S592)
        {
            FixedArray<float, 8>  _S607 = dist_coeffs_29;
            float2  _S608 = DistThinPrism_distort_0(_S596, &_S607);
            _S592 = (dot_0(_S596, _S608)) >= 0.0f;
        }
        else
        {
            _S592 = false;
        }
        if(!_S592)
        {
            *uv_31 = _S596;
            _S592 = false;
            break;
        }
        if(_S595)
        {
            k_11 = (1.0f - theta_9 * theta_9 / 24.0f) / _S594;
        }
        else
        {
            k_11 = 2.0f * (F32_sin((0.5f * theta_9))) / r_27;
        }
        float2  _S609 = _S593 * make_float2 (k_11);
        FixedArray<float, 8>  _S610 = dist_coeffs_29;
        float2  _S611 = DistThinPrism_distort_0(_S609, &_S610);
        *uv_31 = make_float2 (intrins_16.x * _S611.x + intrins_16.z, intrins_16.y * _S611.y + intrins_16.w);
        _S592 = true;
        break;
    }
    return _S592;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac_prism(float3  p_view_17, float4  intrins_17, FixedArray<float, 8>  dist_coeffs_30)
{
    float2  _S612 = float2 {p_view_17.x, p_view_17.y};
    float _S613 = p_view_17.z;
    float2  _S614 = _S612 * make_float2 (0.0f);
    float _S615 = _S613 * _S613;
    float2  s_diff_uv_6 = (make_float2 (1.0f, 0.0f) * make_float2 (_S613) - _S614) / make_float2 (_S615);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S616;
    (&_S616)->primal_0 = _S612 / make_float2 (_S613);
    (&_S616)->differential_0 = s_diff_uv_6;
    FixedArray<float, 8>  _S617 = dist_coeffs_30;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S618 = s_fwd_DistThinPrism_distort_0(&_S616, &_S617);
    float fx_7 = intrins_17.x;
    float fy_7 = intrins_17.y;
    float _S619 = _S618.differential_0.y * fy_7;
    Matrix<float, 2, 3>  J_7;
    *&(((&J_7)->rows + (int(0)))->x) = _S618.differential_0.x * fx_7;
    *&(((&J_7)->rows + (int(1)))->x) = _S619;
    float2  s_diff_uv_7 = (make_float2 (0.0f, 1.0f) * make_float2 (_S613) - _S614) / make_float2 (_S615);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S620;
    (&_S620)->primal_0 = _S612 / make_float2 (_S613);
    (&_S620)->differential_0 = s_diff_uv_7;
    FixedArray<float, 8>  _S621 = dist_coeffs_30;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S622 = s_fwd_DistThinPrism_distort_0(&_S620, &_S621);
    float _S623 = _S622.differential_0.y * fy_7;
    *&(((&J_7)->rows + (int(0)))->y) = _S622.differential_0.x * fx_7;
    *&(((&J_7)->rows + (int(1)))->y) = _S623;
    float2  s_diff_uv_8 = (make_float2 (0.0f, 0.0f) * make_float2 (_S613) - _S612) / make_float2 (_S615);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S624;
    (&_S624)->primal_0 = _S612 / make_float2 (_S613);
    (&_S624)->differential_0 = s_diff_uv_8;
    FixedArray<float, 8>  _S625 = dist_coeffs_30;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S626 = s_fwd_DistThinPrism_distort_0(&_S624, &_S625);
    float _S627 = _S626.differential_0.y * fy_7;
    *&(((&J_7)->rows + (int(0)))->z) = _S626.differential_0.x * fx_7;
    *&(((&J_7)->rows + (int(1)))->z) = _S627;
    return J_7;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac_prism(float3  p_view_18, float4  intrins_18, FixedArray<float, 8>  dist_coeffs_31)
{
    Matrix<float, 2, 3>  J_8;
    float2  _S628 = float2 {p_view_18.x, p_view_18.y};
    float2  _S629 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S630;
    (&_S630)->primal_0 = _S628;
    (&_S630)->differential_0 = _S629;
    DiffPair_float_0 _S631 = s_fwd_length_impl_0(&_S630);
    float _S632 = p_view_18.z;
    DiffPair_float_0 _S633;
    (&_S633)->primal_0 = _S631.primal_0;
    (&_S633)->differential_0 = _S631.differential_0;
    DiffPair_float_0 _S634;
    (&_S634)->primal_0 = _S632;
    (&_S634)->differential_0 = 0.0f;
    DiffPair_float_0 _S635 = _d_atan2_0(&_S633, &_S634);
    float k_12;
    float s_diff_k_4;
    if((_S635.primal_0) < 0.00100000004749745f)
    {
        float _S636 = _S635.differential_0 * _S635.primal_0;
        float _S637 = 1.0f - _S635.primal_0 * _S635.primal_0 / 3.0f;
        float _S638 = ((0.0f - (_S636 + _S636) * 0.3333333432674408f) * _S632 - _S637 * 0.0f) / (_S632 * _S632);
        k_12 = _S637 / _S632;
        s_diff_k_4 = _S638;
    }
    else
    {
        float _S639 = (_S635.differential_0 * _S631.primal_0 - _S635.primal_0 * _S631.differential_0) / (_S631.primal_0 * _S631.primal_0);
        k_12 = _S635.primal_0 / _S631.primal_0;
        s_diff_k_4 = _S639;
    }
    float2  _S640 = _S628 * make_float2 (k_12);
    float2  _S641 = _S629 * make_float2 (k_12) + make_float2 (s_diff_k_4) * _S628;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S642;
    (&_S642)->primal_0 = _S640;
    (&_S642)->differential_0 = _S641;
    FixedArray<float, 8>  _S643 = dist_coeffs_31;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S644 = s_fwd_DistThinPrism_distort_0(&_S642, &_S643);
    float fx_8 = intrins_18.x;
    float fy_8 = intrins_18.y;
    float _S645 = _S644.differential_0.y * fy_8;
    *&(((&J_8)->rows + (int(0)))->x) = _S644.differential_0.x * fx_8;
    *&(((&J_8)->rows + (int(1)))->x) = _S645;
    float2  _S646 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S647;
    (&_S647)->primal_0 = _S628;
    (&_S647)->differential_0 = _S646;
    DiffPair_float_0 _S648 = s_fwd_length_impl_0(&_S647);
    DiffPair_float_0 _S649;
    (&_S649)->primal_0 = _S648.primal_0;
    (&_S649)->differential_0 = _S648.differential_0;
    DiffPair_float_0 _S650;
    (&_S650)->primal_0 = _S632;
    (&_S650)->differential_0 = 0.0f;
    DiffPair_float_0 _S651 = _d_atan2_0(&_S649, &_S650);
    if((_S651.primal_0) < 0.00100000004749745f)
    {
        float _S652 = _S651.differential_0 * _S651.primal_0;
        float _S653 = 1.0f - _S651.primal_0 * _S651.primal_0 / 3.0f;
        float _S654 = ((0.0f - (_S652 + _S652) * 0.3333333432674408f) * _S632 - _S653 * 0.0f) / (_S632 * _S632);
        k_12 = _S653 / _S632;
        s_diff_k_4 = _S654;
    }
    else
    {
        float _S655 = (_S651.differential_0 * _S648.primal_0 - _S651.primal_0 * _S648.differential_0) / (_S648.primal_0 * _S648.primal_0);
        k_12 = _S651.primal_0 / _S648.primal_0;
        s_diff_k_4 = _S655;
    }
    float2  _S656 = _S628 * make_float2 (k_12);
    float2  _S657 = _S646 * make_float2 (k_12) + make_float2 (s_diff_k_4) * _S628;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S658;
    (&_S658)->primal_0 = _S656;
    (&_S658)->differential_0 = _S657;
    FixedArray<float, 8>  _S659 = dist_coeffs_31;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S660 = s_fwd_DistThinPrism_distort_0(&_S658, &_S659);
    float _S661 = _S660.differential_0.y * fy_8;
    *&(((&J_8)->rows + (int(0)))->y) = _S660.differential_0.x * fx_8;
    *&(((&J_8)->rows + (int(1)))->y) = _S661;
    float2  _S662 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S663;
    (&_S663)->primal_0 = _S628;
    (&_S663)->differential_0 = _S662;
    DiffPair_float_0 _S664 = s_fwd_length_impl_0(&_S663);
    DiffPair_float_0 _S665;
    (&_S665)->primal_0 = _S664.primal_0;
    (&_S665)->differential_0 = _S664.differential_0;
    DiffPair_float_0 _S666;
    (&_S666)->primal_0 = _S632;
    (&_S666)->differential_0 = 1.0f;
    DiffPair_float_0 _S667 = _d_atan2_0(&_S665, &_S666);
    if((_S667.primal_0) < 0.00100000004749745f)
    {
        float _S668 = _S667.differential_0 * _S667.primal_0;
        float _S669 = 1.0f - _S667.primal_0 * _S667.primal_0 / 3.0f;
        float _S670 = ((0.0f - (_S668 + _S668) * 0.3333333432674408f) * _S632 - _S669) / (_S632 * _S632);
        k_12 = _S669 / _S632;
        s_diff_k_4 = _S670;
    }
    else
    {
        float _S671 = (_S667.differential_0 * _S664.primal_0 - _S667.primal_0 * _S664.differential_0) / (_S664.primal_0 * _S664.primal_0);
        k_12 = _S667.primal_0 / _S664.primal_0;
        s_diff_k_4 = _S671;
    }
    float2  _S672 = _S628 * make_float2 (k_12);
    float2  _S673 = _S662 * make_float2 (k_12) + make_float2 (s_diff_k_4) * _S628;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S674;
    (&_S674)->primal_0 = _S672;
    (&_S674)->differential_0 = _S673;
    FixedArray<float, 8>  _S675 = dist_coeffs_31;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S676 = s_fwd_DistThinPrism_distort_0(&_S674, &_S675);
    float _S677 = _S676.differential_0.y * fy_8;
    *&(((&J_8)->rows + (int(0)))->z) = _S676.differential_0.x * fx_8;
    *&(((&J_8)->rows + (int(1)))->z) = _S677;
    return J_8;
}

inline __device__ Matrix<float, 2, 3>  equisolid_proj_jac_prism(float3  p_view_19, float4  intrins_19, FixedArray<float, 8>  dist_coeffs_32)
{
    Matrix<float, 2, 3>  J_9;
    float2  _S678 = float2 {p_view_19.x, p_view_19.y};
    float2  _S679 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S680;
    (&_S680)->primal_0 = _S678;
    (&_S680)->differential_0 = _S679;
    DiffPair_float_0 _S681 = s_fwd_length_impl_0(&_S680);
    float _S682 = p_view_19.z;
    DiffPair_float_0 _S683;
    (&_S683)->primal_0 = _S681.primal_0;
    (&_S683)->differential_0 = _S681.differential_0;
    DiffPair_float_0 _S684;
    (&_S684)->primal_0 = _S682;
    (&_S684)->differential_0 = 0.0f;
    DiffPair_float_0 _S685 = _d_atan2_0(&_S683, &_S684);
    float k_13;
    float s_diff_k_5;
    if((_S681.primal_0) < 9.99999997475242708e-07f)
    {
        float _S686 = _S685.differential_0 * _S685.primal_0;
        float _S687 = 1.0f - _S685.primal_0 * _S685.primal_0 / 24.0f;
        float _S688 = ((0.0f - (_S686 + _S686) * 0.0416666679084301f) * _S682 - _S687 * 0.0f) / (_S682 * _S682);
        k_13 = _S687 / _S682;
        s_diff_k_5 = _S688;
    }
    else
    {
        float _S689 = _S685.differential_0 * 0.5f;
        DiffPair_float_0 _S690;
        (&_S690)->primal_0 = 0.5f * _S685.primal_0;
        (&_S690)->differential_0 = _S689;
        DiffPair_float_0 _S691 = _d_sin_0(&_S690);
        float _S692 = 2.0f * _S691.primal_0;
        float _S693 = (_S691.differential_0 * 2.0f * _S681.primal_0 - _S692 * _S681.differential_0) / (_S681.primal_0 * _S681.primal_0);
        k_13 = _S692 / _S681.primal_0;
        s_diff_k_5 = _S693;
    }
    float2  _S694 = _S678 * make_float2 (k_13);
    float2  _S695 = _S679 * make_float2 (k_13) + make_float2 (s_diff_k_5) * _S678;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S696;
    (&_S696)->primal_0 = _S694;
    (&_S696)->differential_0 = _S695;
    FixedArray<float, 8>  _S697 = dist_coeffs_32;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S698 = s_fwd_DistThinPrism_distort_0(&_S696, &_S697);
    float fx_9 = intrins_19.x;
    float fy_9 = intrins_19.y;
    float _S699 = _S698.differential_0.y * fy_9;
    *&(((&J_9)->rows + (int(0)))->x) = _S698.differential_0.x * fx_9;
    *&(((&J_9)->rows + (int(1)))->x) = _S699;
    float2  _S700 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S701;
    (&_S701)->primal_0 = _S678;
    (&_S701)->differential_0 = _S700;
    DiffPair_float_0 _S702 = s_fwd_length_impl_0(&_S701);
    DiffPair_float_0 _S703;
    (&_S703)->primal_0 = _S702.primal_0;
    (&_S703)->differential_0 = _S702.differential_0;
    DiffPair_float_0 _S704;
    (&_S704)->primal_0 = _S682;
    (&_S704)->differential_0 = 0.0f;
    DiffPair_float_0 _S705 = _d_atan2_0(&_S703, &_S704);
    if((_S702.primal_0) < 9.99999997475242708e-07f)
    {
        float _S706 = _S705.differential_0 * _S705.primal_0;
        float _S707 = 1.0f - _S705.primal_0 * _S705.primal_0 / 24.0f;
        float _S708 = ((0.0f - (_S706 + _S706) * 0.0416666679084301f) * _S682 - _S707 * 0.0f) / (_S682 * _S682);
        k_13 = _S707 / _S682;
        s_diff_k_5 = _S708;
    }
    else
    {
        float _S709 = _S705.differential_0 * 0.5f;
        DiffPair_float_0 _S710;
        (&_S710)->primal_0 = 0.5f * _S705.primal_0;
        (&_S710)->differential_0 = _S709;
        DiffPair_float_0 _S711 = _d_sin_0(&_S710);
        float _S712 = 2.0f * _S711.primal_0;
        float _S713 = (_S711.differential_0 * 2.0f * _S702.primal_0 - _S712 * _S702.differential_0) / (_S702.primal_0 * _S702.primal_0);
        k_13 = _S712 / _S702.primal_0;
        s_diff_k_5 = _S713;
    }
    float2  _S714 = _S678 * make_float2 (k_13);
    float2  _S715 = _S700 * make_float2 (k_13) + make_float2 (s_diff_k_5) * _S678;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S716;
    (&_S716)->primal_0 = _S714;
    (&_S716)->differential_0 = _S715;
    FixedArray<float, 8>  _S717 = dist_coeffs_32;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S718 = s_fwd_DistThinPrism_distort_0(&_S716, &_S717);
    float _S719 = _S718.differential_0.y * fy_9;
    *&(((&J_9)->rows + (int(0)))->y) = _S718.differential_0.x * fx_9;
    *&(((&J_9)->rows + (int(1)))->y) = _S719;
    float2  _S720 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S721;
    (&_S721)->primal_0 = _S678;
    (&_S721)->differential_0 = _S720;
    DiffPair_float_0 _S722 = s_fwd_length_impl_0(&_S721);
    DiffPair_float_0 _S723;
    (&_S723)->primal_0 = _S722.primal_0;
    (&_S723)->differential_0 = _S722.differential_0;
    DiffPair_float_0 _S724;
    (&_S724)->primal_0 = _S682;
    (&_S724)->differential_0 = 1.0f;
    DiffPair_float_0 _S725 = _d_atan2_0(&_S723, &_S724);
    if((_S722.primal_0) < 9.99999997475242708e-07f)
    {
        float _S726 = _S725.differential_0 * _S725.primal_0;
        float _S727 = 1.0f - _S725.primal_0 * _S725.primal_0 / 24.0f;
        float _S728 = ((0.0f - (_S726 + _S726) * 0.0416666679084301f) * _S682 - _S727) / (_S682 * _S682);
        k_13 = _S727 / _S682;
        s_diff_k_5 = _S728;
    }
    else
    {
        float _S729 = _S725.differential_0 * 0.5f;
        DiffPair_float_0 _S730;
        (&_S730)->primal_0 = 0.5f * _S725.primal_0;
        (&_S730)->differential_0 = _S729;
        DiffPair_float_0 _S731 = _d_sin_0(&_S730);
        float _S732 = 2.0f * _S731.primal_0;
        float _S733 = (_S731.differential_0 * 2.0f * _S722.primal_0 - _S732 * _S722.differential_0) / (_S722.primal_0 * _S722.primal_0);
        k_13 = _S732 / _S722.primal_0;
        s_diff_k_5 = _S733;
    }
    float2  _S734 = _S678 * make_float2 (k_13);
    float2  _S735 = _S720 * make_float2 (k_13) + make_float2 (s_diff_k_5) * _S678;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S736;
    (&_S736)->primal_0 = _S734;
    (&_S736)->differential_0 = _S735;
    FixedArray<float, 8>  _S737 = dist_coeffs_32;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S738 = s_fwd_DistThinPrism_distort_0(&_S736, &_S737);
    float _S739 = _S738.differential_0.y * fy_9;
    *&(((&J_9)->rows + (int(0)))->z) = _S738.differential_0.x * fx_9;
    *&(((&J_9)->rows + (int(1)))->z) = _S739;
    return J_9;
}

inline __device__ float2  distort_point_prism(float2  uv_32, int camera_model_8, FixedArray<float, 8>  dist_coeffs_33)
{
    float2  _S740;
    for(;;)
    {
        if(camera_model_8 == int(3))
        {
            _S740 = uv_32;
            break;
        }
        float k_14;
        if(camera_model_8 == int(1))
        {
            float r_28 = length_0(uv_32);
            float theta_10 = (F32_atan((r_28)));
            if(r_28 < 0.00100000004749745f)
            {
                k_14 = 1.0f - theta_10 * theta_10 / 6.0f;
            }
            else
            {
                k_14 = theta_10 / r_28;
            }
            _S740 = uv_32 * make_float2 (k_14);
        }
        else
        {
            if(camera_model_8 == int(2))
            {
                float r_29 = length_0(uv_32);
                float theta_11 = (F32_atan((r_29)));
                if(r_29 < 0.00100000004749745f)
                {
                    k_14 = 1.0f - theta_11 * theta_11 / 24.0f;
                }
                else
                {
                    k_14 = 2.0f * (F32_sin((0.5f * theta_11))) / r_29;
                }
                _S740 = uv_32 * make_float2 (k_14);
            }
            else
            {
                _S740 = uv_32;
            }
        }
        FixedArray<float, 8>  _S741 = dist_coeffs_33;
        float2  _S742 = DistThinPrism_distort_0(_S740, &_S741);
        _S740 = _S742;
        break;
    }
    return _S740;
}

inline __device__ bool undistort_point_prism(float2  uv_33, int camera_model_9, FixedArray<float, 8>  dist_coeffs_34, float2  * uv_undist_6)
{
    bool _S743;
    for(;;)
    {
        *uv_undist_6 = make_float2 (0.0f);
        if(camera_model_9 == int(3))
        {
            float lon_4 = uv_33.x;
            float lat_6 = uv_33.y;
            float cl_6 = (F32_cos((lat_6)));
            *uv_undist_6 = make_float2 (cl_6 * (F32_sin((lon_4))), (F32_sin((lat_6)))) / make_float2 ((F32_max((cl_6 * (F32_cos((lon_4)))), (9.999999960041972e-13f))));
            _S743 = true;
            break;
        }
        FixedArray<float, 8>  _S744 = dist_coeffs_34;
        float2  uv_u_6;
        bool _S745 = undistort_point_2(uv_33, &_S744, int(8), &uv_u_6);
        if(!_S745)
        {
            _S743 = false;
            break;
        }
        float2  _S746 = uv_u_6;
        float3  raydir_6;
        if(camera_model_9 == int(1))
        {
            float r_30 = length_0(_S746);
            float s_6;
            if(r_30 < 0.00100000004749745f)
            {
                s_6 = 1.0f - r_30 * r_30 / 6.0f;
            }
            else
            {
                s_6 = (F32_sin((r_30))) / r_30;
            }
            raydir_6 = make_float3 ((_S746 * make_float2 (s_6)).x, (_S746 * make_float2 (s_6)).y, (F32_cos((r_30))));
        }
        else
        {
            if(camera_model_9 == int(2))
            {
                float r_31 = length_0(_S746);
                raydir_6 = make_float3 ((_S746 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_31 * r_31)))))))).x, (_S746 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_31 * r_31)))))))).y, 1.0f - 0.5f * r_31 * r_31);
            }
            else
            {
                raydir_6 = make_float3 (_S746.x, _S746.y, 1.0f);
            }
        }
        *uv_undist_6 = float2 {raydir_6.x, raydir_6.y} / make_float2 ((F32_max((raydir_6.z), (9.999999960041972e-13f))));
        _S743 = true;
        break;
    }
    return _S743;
}

inline __device__ bool unproject_point_prism(float2  uv_34, int camera_model_10, FixedArray<float, 8>  dist_coeffs_35, float3  * raydir_7)
{
    bool _S747;
    for(;;)
    {
        int3  _S748 = make_int3 (int(0));
        float3  _S749 = make_float3 ((float)_S748.x, (float)_S748.y, (float)_S748.z);
        *raydir_7 = _S749;
        if(camera_model_10 == int(3))
        {
            float lon_5 = uv_34.x;
            float lat_7 = uv_34.y;
            float cl_7 = (F32_cos((lat_7)));
            *raydir_7 = make_float3 (cl_7 * (F32_sin((lon_5))), (F32_sin((lat_7))), cl_7 * (F32_cos((lon_5))));
            _S747 = true;
            break;
        }
        FixedArray<float, 8>  _S750 = dist_coeffs_35;
        float2  uv_u_7;
        bool _S751 = undistort_point_2(uv_34, &_S750, int(8), &uv_u_7);
        if(!_S751)
        {
            _S747 = false;
            break;
        }
        float2  _S752 = uv_u_7;
        if(camera_model_10 == int(1))
        {
            float r_32 = length_0(_S752);
            float s_7;
            if(r_32 < 0.00100000004749745f)
            {
                s_7 = 1.0f - r_32 * r_32 / 6.0f;
            }
            else
            {
                s_7 = (F32_sin((r_32))) / r_32;
            }
            *raydir_7 = make_float3 ((_S752 * make_float2 (s_7)).x, (_S752 * make_float2 (s_7)).y, (F32_cos((r_32))));
        }
        else
        {
            if(camera_model_10 == int(2))
            {
                float r_33 = length_0(_S752);
                *raydir_7 = make_float3 ((_S752 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_33 * r_33)))))))).x, (_S752 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_33 * r_33)))))))).y, 1.0f - 0.5f * r_33 * r_33);
            }
            else
            {
                *raydir_7 = make_float3 (_S752.x, _S752.y, 1.0f);
            }
        }
        _S747 = true;
        break;
    }
    return _S747;
}

inline __device__ bool generate_ray_prism(float2  uv_35, int camera_model_11, FixedArray<float, 8>  dist_coeffs_36, float3  * raydir_8)
{
    bool _S753;
    for(;;)
    {
        if(camera_model_11 == int(3))
        {
            float _S754 = uv_35.x;
            if((F32_abs((_S754))) > 3.14159274101257324f)
            {
                _S753 = true;
            }
            else
            {
                _S753 = (F32_abs((uv_35.y))) > 1.57079637050628662f;
            }
            if(_S753)
            {
                int3  _S755 = make_int3 (int(0));
                float3  _S756 = make_float3 ((float)_S755.x, (float)_S755.y, (float)_S755.z);
                *raydir_8 = _S756;
                _S753 = false;
                break;
            }
            float lat_8 = uv_35.y;
            float cl_8 = (F32_cos((lat_8)));
            *raydir_8 = make_float3 (cl_8 * (F32_sin((_S754))), (F32_sin((lat_8))), cl_8 * (F32_cos((_S754))));
            _S753 = true;
            break;
        }
        FixedArray<float, 8>  _S757 = dist_coeffs_36;
        float2  uv_u_8;
        bool _S758 = undistort_point_2(uv_35, &_S757, int(8), &uv_u_8);
        if(!_S758)
        {
            int3  _S759 = make_int3 (int(0));
            float3  _S760 = make_float3 ((float)_S759.x, (float)_S759.y, (float)_S759.z);
            *raydir_8 = _S760;
            _S753 = false;
            break;
        }
        float2  _S761 = uv_u_8;
        if(camera_model_11 == int(1))
        {
            float r_34 = length_0(_S761);
            if(r_34 >= 3.14159274101257324f)
            {
                int3  _S762 = make_int3 (int(0));
                float3  _S763 = make_float3 ((float)_S762.x, (float)_S762.y, (float)_S762.z);
                *raydir_8 = _S763;
                _S753 = false;
                break;
            }
            float s_8;
            if(r_34 < 0.00100000004749745f)
            {
                s_8 = 1.0f - r_34 * r_34 / 6.0f;
            }
            else
            {
                s_8 = (F32_sin((r_34))) / r_34;
            }
            *raydir_8 = make_float3 ((_S761 * make_float2 (s_8)).x, (_S761 * make_float2 (s_8)).y, (F32_cos((r_34))));
        }
        else
        {
            if(camera_model_11 == int(2))
            {
                float r_35 = length_0(_S761);
                if(r_35 >= 2.0f)
                {
                    int3  _S764 = make_int3 (int(0));
                    float3  _S765 = make_float3 ((float)_S764.x, (float)_S764.y, (float)_S764.z);
                    *raydir_8 = _S765;
                    _S753 = false;
                    break;
                }
                *raydir_8 = make_float3 ((_S761 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_35 * r_35)))))))).x, (_S761 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_35 * r_35)))))))).y, 1.0f - 0.5f * r_35 * r_35);
            }
            else
            {
                *raydir_8 = make_float3 (_S761.x, _S761.y, 1.0f);
            }
        }
        *raydir_8 = normalize_0(*raydir_8);
        _S753 = true;
        break;
    }
    return _S753;
}

inline __device__ bool is_valid_distortion_rational(float2  uv_36, FixedArray<float, 8>  dist_coeffs_37)
{
    float2  _S766 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S767;
    (&_S767)->primal_0 = uv_36;
    (&_S767)->differential_0 = _S766;
    FixedArray<float, 8>  _S768 = dist_coeffs_37;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S769 = s_fwd_DistRational_distort_0(&_S767, &_S768);
    float2  _S770 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S771;
    (&_S771)->primal_0 = uv_36;
    (&_S771)->differential_0 = _S770;
    FixedArray<float, 8>  _S772 = dist_coeffs_37;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S773 = s_fwd_DistRational_distort_0(&_S771, &_S772);
    Matrix<float, 2, 2>  _S774 = transpose_1(makeMatrix<float, 2, 2> (_S769.differential_0, _S773.differential_0));
    float _S775 = (F32_min((determinant_0(_S774)), ((F32_min((_S774.rows[int(0)].x), (_S774.rows[int(1)].y))))));
    bool _S776;
    if(_S775 > 0.25f)
    {
        _S776 = _S775 < 4.0f;
    }
    else
    {
        _S776 = false;
    }
    if(_S776)
    {
        FixedArray<float, 8>  _S777 = dist_coeffs_37;
        float2  _S778 = DistRational_distort_0(uv_36, &_S777);
        _S776 = (dot_0(uv_36, _S778)) >= 0.0f;
    }
    else
    {
        _S776 = false;
    }
    return _S776;
}

inline __device__ bool persp_proj_nav_rational(float3  p_view_20, float4  intrins_20, FixedArray<float, 8>  dist_coeffs_38, float2  * uv_37)
{
    bool _S779;
    for(;;)
    {
        float2  _S780 = float2 {p_view_20.x, p_view_20.y};
        float _S781 = p_view_20.z;
        float2  uv0_3 = _S780 / make_float2 (_S781);
        if(_S781 < 0.0f)
        {
            _S779 = true;
        }
        else
        {
            float2  _S782 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S783;
            (&_S783)->primal_0 = uv0_3;
            (&_S783)->differential_0 = _S782;
            FixedArray<float, 8>  _S784 = dist_coeffs_38;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S785 = s_fwd_DistRational_distort_0(&_S783, &_S784);
            float2  _S786 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S787;
            (&_S787)->primal_0 = uv0_3;
            (&_S787)->differential_0 = _S786;
            FixedArray<float, 8>  _S788 = dist_coeffs_38;
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S789 = s_fwd_DistRational_distort_0(&_S787, &_S788);
            Matrix<float, 2, 2>  _S790 = transpose_1(makeMatrix<float, 2, 2> (_S785.differential_0, _S789.differential_0));
            float _S791 = (F32_min((determinant_0(_S790)), ((F32_min((_S790.rows[int(0)].x), (_S790.rows[int(1)].y))))));
            if(_S791 > 0.25f)
            {
                _S779 = _S791 < 4.0f;
            }
            else
            {
                _S779 = false;
            }
            if(_S779)
            {
                FixedArray<float, 8>  _S792 = dist_coeffs_38;
                float2  _S793 = DistRational_distort_0(uv0_3, &_S792);
                _S779 = (dot_0(uv0_3, _S793)) >= 0.0f;
            }
            else
            {
                _S779 = false;
            }
            _S779 = !_S779;
        }
        if(_S779)
        {
            *uv_37 = uv0_3;
            _S779 = false;
            break;
        }
        float2  uv_38 = _S780 / make_float2 (_S781);
        FixedArray<float, 8>  _S794 = dist_coeffs_38;
        float2  _S795 = DistRational_distort_0(uv_38, &_S794);
        *uv_37 = make_float2 (intrins_20.x * _S795.x + intrins_20.z, intrins_20.y * _S795.y + intrins_20.w);
        _S779 = true;
        break;
    }
    return _S779;
}

inline __device__ bool fisheye_proj_nav_rational(float3  p_view_21, float4  intrins_21, FixedArray<float, 8>  dist_coeffs_39, float2  * uv_39)
{
    bool _S796;
    for(;;)
    {
        float2  _S797 = float2 {p_view_21.x, p_view_21.y};
        float r_36 = length_0(_S797);
        float _S798 = p_view_21.z;
        float theta_12 = (F32_atan2((r_36), (_S798)));
        bool _S799 = theta_12 < 0.00100000004749745f;
        float k_15;
        if(_S799)
        {
            k_15 = (1.0f - theta_12 * theta_12 / 3.0f) / _S798;
        }
        else
        {
            k_15 = theta_12 / r_36;
        }
        float2  _S800 = _S797 * make_float2 (k_15);
        float2  _S801 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S802;
        (&_S802)->primal_0 = _S800;
        (&_S802)->differential_0 = _S801;
        FixedArray<float, 8>  _S803 = dist_coeffs_39;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S804 = s_fwd_DistRational_distort_0(&_S802, &_S803);
        float2  _S805 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S806;
        (&_S806)->primal_0 = _S800;
        (&_S806)->differential_0 = _S805;
        FixedArray<float, 8>  _S807 = dist_coeffs_39;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S808 = s_fwd_DistRational_distort_0(&_S806, &_S807);
        Matrix<float, 2, 2>  _S809 = transpose_1(makeMatrix<float, 2, 2> (_S804.differential_0, _S808.differential_0));
        float _S810 = (F32_min((determinant_0(_S809)), ((F32_min((_S809.rows[int(0)].x), (_S809.rows[int(1)].y))))));
        if(_S810 > 0.25f)
        {
            _S796 = _S810 < 4.0f;
        }
        else
        {
            _S796 = false;
        }
        if(_S796)
        {
            FixedArray<float, 8>  _S811 = dist_coeffs_39;
            float2  _S812 = DistRational_distort_0(_S800, &_S811);
            _S796 = (dot_0(_S800, _S812)) >= 0.0f;
        }
        else
        {
            _S796 = false;
        }
        if(!_S796)
        {
            *uv_39 = _S800;
            _S796 = false;
            break;
        }
        if(_S799)
        {
            k_15 = (1.0f - theta_12 * theta_12 / 3.0f) / _S798;
        }
        else
        {
            k_15 = theta_12 / r_36;
        }
        float2  _S813 = _S797 * make_float2 (k_15);
        FixedArray<float, 8>  _S814 = dist_coeffs_39;
        float2  _S815 = DistRational_distort_0(_S813, &_S814);
        *uv_39 = make_float2 (intrins_21.x * _S815.x + intrins_21.z, intrins_21.y * _S815.y + intrins_21.w);
        _S796 = true;
        break;
    }
    return _S796;
}

inline __device__ bool equisolid_proj_nav_rational(float3  p_view_22, float4  intrins_22, FixedArray<float, 8>  dist_coeffs_40, float2  * uv_40)
{
    bool _S816;
    for(;;)
    {
        float2  _S817 = float2 {p_view_22.x, p_view_22.y};
        float r_37 = length_0(_S817);
        float _S818 = p_view_22.z;
        float theta_13 = (F32_atan2((r_37), (_S818)));
        bool _S819 = r_37 < 9.99999997475242708e-07f;
        float k_16;
        if(_S819)
        {
            k_16 = (1.0f - theta_13 * theta_13 / 24.0f) / _S818;
        }
        else
        {
            k_16 = 2.0f * (F32_sin((0.5f * theta_13))) / r_37;
        }
        float2  _S820 = _S817 * make_float2 (k_16);
        float2  _S821 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S822;
        (&_S822)->primal_0 = _S820;
        (&_S822)->differential_0 = _S821;
        FixedArray<float, 8>  _S823 = dist_coeffs_40;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S824 = s_fwd_DistRational_distort_0(&_S822, &_S823);
        float2  _S825 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S826;
        (&_S826)->primal_0 = _S820;
        (&_S826)->differential_0 = _S825;
        FixedArray<float, 8>  _S827 = dist_coeffs_40;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S828 = s_fwd_DistRational_distort_0(&_S826, &_S827);
        Matrix<float, 2, 2>  _S829 = transpose_1(makeMatrix<float, 2, 2> (_S824.differential_0, _S828.differential_0));
        float _S830 = (F32_min((determinant_0(_S829)), ((F32_min((_S829.rows[int(0)].x), (_S829.rows[int(1)].y))))));
        if(_S830 > 0.25f)
        {
            _S816 = _S830 < 4.0f;
        }
        else
        {
            _S816 = false;
        }
        if(_S816)
        {
            FixedArray<float, 8>  _S831 = dist_coeffs_40;
            float2  _S832 = DistRational_distort_0(_S820, &_S831);
            _S816 = (dot_0(_S820, _S832)) >= 0.0f;
        }
        else
        {
            _S816 = false;
        }
        if(!_S816)
        {
            *uv_40 = _S820;
            _S816 = false;
            break;
        }
        if(_S819)
        {
            k_16 = (1.0f - theta_13 * theta_13 / 24.0f) / _S818;
        }
        else
        {
            k_16 = 2.0f * (F32_sin((0.5f * theta_13))) / r_37;
        }
        float2  _S833 = _S817 * make_float2 (k_16);
        FixedArray<float, 8>  _S834 = dist_coeffs_40;
        float2  _S835 = DistRational_distort_0(_S833, &_S834);
        *uv_40 = make_float2 (intrins_22.x * _S835.x + intrins_22.z, intrins_22.y * _S835.y + intrins_22.w);
        _S816 = true;
        break;
    }
    return _S816;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac_rational(float3  p_view_23, float4  intrins_23, FixedArray<float, 8>  dist_coeffs_41)
{
    float2  _S836 = float2 {p_view_23.x, p_view_23.y};
    float _S837 = p_view_23.z;
    float2  _S838 = _S836 * make_float2 (0.0f);
    float _S839 = _S837 * _S837;
    float2  s_diff_uv_9 = (make_float2 (1.0f, 0.0f) * make_float2 (_S837) - _S838) / make_float2 (_S839);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S840;
    (&_S840)->primal_0 = _S836 / make_float2 (_S837);
    (&_S840)->differential_0 = s_diff_uv_9;
    FixedArray<float, 8>  _S841 = dist_coeffs_41;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S842 = s_fwd_DistRational_distort_0(&_S840, &_S841);
    float fx_10 = intrins_23.x;
    float fy_10 = intrins_23.y;
    float _S843 = _S842.differential_0.y * fy_10;
    Matrix<float, 2, 3>  J_10;
    *&(((&J_10)->rows + (int(0)))->x) = _S842.differential_0.x * fx_10;
    *&(((&J_10)->rows + (int(1)))->x) = _S843;
    float2  s_diff_uv_10 = (make_float2 (0.0f, 1.0f) * make_float2 (_S837) - _S838) / make_float2 (_S839);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S844;
    (&_S844)->primal_0 = _S836 / make_float2 (_S837);
    (&_S844)->differential_0 = s_diff_uv_10;
    FixedArray<float, 8>  _S845 = dist_coeffs_41;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S846 = s_fwd_DistRational_distort_0(&_S844, &_S845);
    float _S847 = _S846.differential_0.y * fy_10;
    *&(((&J_10)->rows + (int(0)))->y) = _S846.differential_0.x * fx_10;
    *&(((&J_10)->rows + (int(1)))->y) = _S847;
    float2  s_diff_uv_11 = (make_float2 (0.0f, 0.0f) * make_float2 (_S837) - _S836) / make_float2 (_S839);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S848;
    (&_S848)->primal_0 = _S836 / make_float2 (_S837);
    (&_S848)->differential_0 = s_diff_uv_11;
    FixedArray<float, 8>  _S849 = dist_coeffs_41;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S850 = s_fwd_DistRational_distort_0(&_S848, &_S849);
    float _S851 = _S850.differential_0.y * fy_10;
    *&(((&J_10)->rows + (int(0)))->z) = _S850.differential_0.x * fx_10;
    *&(((&J_10)->rows + (int(1)))->z) = _S851;
    return J_10;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac_rational(float3  p_view_24, float4  intrins_24, FixedArray<float, 8>  dist_coeffs_42)
{
    Matrix<float, 2, 3>  J_11;
    float2  _S852 = float2 {p_view_24.x, p_view_24.y};
    float2  _S853 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S854;
    (&_S854)->primal_0 = _S852;
    (&_S854)->differential_0 = _S853;
    DiffPair_float_0 _S855 = s_fwd_length_impl_0(&_S854);
    float _S856 = p_view_24.z;
    DiffPair_float_0 _S857;
    (&_S857)->primal_0 = _S855.primal_0;
    (&_S857)->differential_0 = _S855.differential_0;
    DiffPair_float_0 _S858;
    (&_S858)->primal_0 = _S856;
    (&_S858)->differential_0 = 0.0f;
    DiffPair_float_0 _S859 = _d_atan2_0(&_S857, &_S858);
    float k_17;
    float s_diff_k_6;
    if((_S859.primal_0) < 0.00100000004749745f)
    {
        float _S860 = _S859.differential_0 * _S859.primal_0;
        float _S861 = 1.0f - _S859.primal_0 * _S859.primal_0 / 3.0f;
        float _S862 = ((0.0f - (_S860 + _S860) * 0.3333333432674408f) * _S856 - _S861 * 0.0f) / (_S856 * _S856);
        k_17 = _S861 / _S856;
        s_diff_k_6 = _S862;
    }
    else
    {
        float _S863 = (_S859.differential_0 * _S855.primal_0 - _S859.primal_0 * _S855.differential_0) / (_S855.primal_0 * _S855.primal_0);
        k_17 = _S859.primal_0 / _S855.primal_0;
        s_diff_k_6 = _S863;
    }
    float2  _S864 = _S852 * make_float2 (k_17);
    float2  _S865 = _S853 * make_float2 (k_17) + make_float2 (s_diff_k_6) * _S852;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S866;
    (&_S866)->primal_0 = _S864;
    (&_S866)->differential_0 = _S865;
    FixedArray<float, 8>  _S867 = dist_coeffs_42;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S868 = s_fwd_DistRational_distort_0(&_S866, &_S867);
    float fx_11 = intrins_24.x;
    float fy_11 = intrins_24.y;
    float _S869 = _S868.differential_0.y * fy_11;
    *&(((&J_11)->rows + (int(0)))->x) = _S868.differential_0.x * fx_11;
    *&(((&J_11)->rows + (int(1)))->x) = _S869;
    float2  _S870 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S871;
    (&_S871)->primal_0 = _S852;
    (&_S871)->differential_0 = _S870;
    DiffPair_float_0 _S872 = s_fwd_length_impl_0(&_S871);
    DiffPair_float_0 _S873;
    (&_S873)->primal_0 = _S872.primal_0;
    (&_S873)->differential_0 = _S872.differential_0;
    DiffPair_float_0 _S874;
    (&_S874)->primal_0 = _S856;
    (&_S874)->differential_0 = 0.0f;
    DiffPair_float_0 _S875 = _d_atan2_0(&_S873, &_S874);
    if((_S875.primal_0) < 0.00100000004749745f)
    {
        float _S876 = _S875.differential_0 * _S875.primal_0;
        float _S877 = 1.0f - _S875.primal_0 * _S875.primal_0 / 3.0f;
        float _S878 = ((0.0f - (_S876 + _S876) * 0.3333333432674408f) * _S856 - _S877 * 0.0f) / (_S856 * _S856);
        k_17 = _S877 / _S856;
        s_diff_k_6 = _S878;
    }
    else
    {
        float _S879 = (_S875.differential_0 * _S872.primal_0 - _S875.primal_0 * _S872.differential_0) / (_S872.primal_0 * _S872.primal_0);
        k_17 = _S875.primal_0 / _S872.primal_0;
        s_diff_k_6 = _S879;
    }
    float2  _S880 = _S852 * make_float2 (k_17);
    float2  _S881 = _S870 * make_float2 (k_17) + make_float2 (s_diff_k_6) * _S852;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S882;
    (&_S882)->primal_0 = _S880;
    (&_S882)->differential_0 = _S881;
    FixedArray<float, 8>  _S883 = dist_coeffs_42;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S884 = s_fwd_DistRational_distort_0(&_S882, &_S883);
    float _S885 = _S884.differential_0.y * fy_11;
    *&(((&J_11)->rows + (int(0)))->y) = _S884.differential_0.x * fx_11;
    *&(((&J_11)->rows + (int(1)))->y) = _S885;
    float2  _S886 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S887;
    (&_S887)->primal_0 = _S852;
    (&_S887)->differential_0 = _S886;
    DiffPair_float_0 _S888 = s_fwd_length_impl_0(&_S887);
    DiffPair_float_0 _S889;
    (&_S889)->primal_0 = _S888.primal_0;
    (&_S889)->differential_0 = _S888.differential_0;
    DiffPair_float_0 _S890;
    (&_S890)->primal_0 = _S856;
    (&_S890)->differential_0 = 1.0f;
    DiffPair_float_0 _S891 = _d_atan2_0(&_S889, &_S890);
    if((_S891.primal_0) < 0.00100000004749745f)
    {
        float _S892 = _S891.differential_0 * _S891.primal_0;
        float _S893 = 1.0f - _S891.primal_0 * _S891.primal_0 / 3.0f;
        float _S894 = ((0.0f - (_S892 + _S892) * 0.3333333432674408f) * _S856 - _S893) / (_S856 * _S856);
        k_17 = _S893 / _S856;
        s_diff_k_6 = _S894;
    }
    else
    {
        float _S895 = (_S891.differential_0 * _S888.primal_0 - _S891.primal_0 * _S888.differential_0) / (_S888.primal_0 * _S888.primal_0);
        k_17 = _S891.primal_0 / _S888.primal_0;
        s_diff_k_6 = _S895;
    }
    float2  _S896 = _S852 * make_float2 (k_17);
    float2  _S897 = _S886 * make_float2 (k_17) + make_float2 (s_diff_k_6) * _S852;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S898;
    (&_S898)->primal_0 = _S896;
    (&_S898)->differential_0 = _S897;
    FixedArray<float, 8>  _S899 = dist_coeffs_42;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S900 = s_fwd_DistRational_distort_0(&_S898, &_S899);
    float _S901 = _S900.differential_0.y * fy_11;
    *&(((&J_11)->rows + (int(0)))->z) = _S900.differential_0.x * fx_11;
    *&(((&J_11)->rows + (int(1)))->z) = _S901;
    return J_11;
}

inline __device__ Matrix<float, 2, 3>  equisolid_proj_jac_rational(float3  p_view_25, float4  intrins_25, FixedArray<float, 8>  dist_coeffs_43)
{
    Matrix<float, 2, 3>  J_12;
    float2  _S902 = float2 {p_view_25.x, p_view_25.y};
    float2  _S903 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S904;
    (&_S904)->primal_0 = _S902;
    (&_S904)->differential_0 = _S903;
    DiffPair_float_0 _S905 = s_fwd_length_impl_0(&_S904);
    float _S906 = p_view_25.z;
    DiffPair_float_0 _S907;
    (&_S907)->primal_0 = _S905.primal_0;
    (&_S907)->differential_0 = _S905.differential_0;
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S906;
    (&_S908)->differential_0 = 0.0f;
    DiffPair_float_0 _S909 = _d_atan2_0(&_S907, &_S908);
    float k_18;
    float s_diff_k_7;
    if((_S905.primal_0) < 9.99999997475242708e-07f)
    {
        float _S910 = _S909.differential_0 * _S909.primal_0;
        float _S911 = 1.0f - _S909.primal_0 * _S909.primal_0 / 24.0f;
        float _S912 = ((0.0f - (_S910 + _S910) * 0.0416666679084301f) * _S906 - _S911 * 0.0f) / (_S906 * _S906);
        k_18 = _S911 / _S906;
        s_diff_k_7 = _S912;
    }
    else
    {
        float _S913 = _S909.differential_0 * 0.5f;
        DiffPair_float_0 _S914;
        (&_S914)->primal_0 = 0.5f * _S909.primal_0;
        (&_S914)->differential_0 = _S913;
        DiffPair_float_0 _S915 = _d_sin_0(&_S914);
        float _S916 = 2.0f * _S915.primal_0;
        float _S917 = (_S915.differential_0 * 2.0f * _S905.primal_0 - _S916 * _S905.differential_0) / (_S905.primal_0 * _S905.primal_0);
        k_18 = _S916 / _S905.primal_0;
        s_diff_k_7 = _S917;
    }
    float2  _S918 = _S902 * make_float2 (k_18);
    float2  _S919 = _S903 * make_float2 (k_18) + make_float2 (s_diff_k_7) * _S902;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S920;
    (&_S920)->primal_0 = _S918;
    (&_S920)->differential_0 = _S919;
    FixedArray<float, 8>  _S921 = dist_coeffs_43;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S922 = s_fwd_DistRational_distort_0(&_S920, &_S921);
    float fx_12 = intrins_25.x;
    float fy_12 = intrins_25.y;
    float _S923 = _S922.differential_0.y * fy_12;
    *&(((&J_12)->rows + (int(0)))->x) = _S922.differential_0.x * fx_12;
    *&(((&J_12)->rows + (int(1)))->x) = _S923;
    float2  _S924 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S925;
    (&_S925)->primal_0 = _S902;
    (&_S925)->differential_0 = _S924;
    DiffPair_float_0 _S926 = s_fwd_length_impl_0(&_S925);
    DiffPair_float_0 _S927;
    (&_S927)->primal_0 = _S926.primal_0;
    (&_S927)->differential_0 = _S926.differential_0;
    DiffPair_float_0 _S928;
    (&_S928)->primal_0 = _S906;
    (&_S928)->differential_0 = 0.0f;
    DiffPair_float_0 _S929 = _d_atan2_0(&_S927, &_S928);
    if((_S926.primal_0) < 9.99999997475242708e-07f)
    {
        float _S930 = _S929.differential_0 * _S929.primal_0;
        float _S931 = 1.0f - _S929.primal_0 * _S929.primal_0 / 24.0f;
        float _S932 = ((0.0f - (_S930 + _S930) * 0.0416666679084301f) * _S906 - _S931 * 0.0f) / (_S906 * _S906);
        k_18 = _S931 / _S906;
        s_diff_k_7 = _S932;
    }
    else
    {
        float _S933 = _S929.differential_0 * 0.5f;
        DiffPair_float_0 _S934;
        (&_S934)->primal_0 = 0.5f * _S929.primal_0;
        (&_S934)->differential_0 = _S933;
        DiffPair_float_0 _S935 = _d_sin_0(&_S934);
        float _S936 = 2.0f * _S935.primal_0;
        float _S937 = (_S935.differential_0 * 2.0f * _S926.primal_0 - _S936 * _S926.differential_0) / (_S926.primal_0 * _S926.primal_0);
        k_18 = _S936 / _S926.primal_0;
        s_diff_k_7 = _S937;
    }
    float2  _S938 = _S902 * make_float2 (k_18);
    float2  _S939 = _S924 * make_float2 (k_18) + make_float2 (s_diff_k_7) * _S902;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S940;
    (&_S940)->primal_0 = _S938;
    (&_S940)->differential_0 = _S939;
    FixedArray<float, 8>  _S941 = dist_coeffs_43;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S942 = s_fwd_DistRational_distort_0(&_S940, &_S941);
    float _S943 = _S942.differential_0.y * fy_12;
    *&(((&J_12)->rows + (int(0)))->y) = _S942.differential_0.x * fx_12;
    *&(((&J_12)->rows + (int(1)))->y) = _S943;
    float2  _S944 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S945;
    (&_S945)->primal_0 = _S902;
    (&_S945)->differential_0 = _S944;
    DiffPair_float_0 _S946 = s_fwd_length_impl_0(&_S945);
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S946.primal_0;
    (&_S947)->differential_0 = _S946.differential_0;
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = _S906;
    (&_S948)->differential_0 = 1.0f;
    DiffPair_float_0 _S949 = _d_atan2_0(&_S947, &_S948);
    if((_S946.primal_0) < 9.99999997475242708e-07f)
    {
        float _S950 = _S949.differential_0 * _S949.primal_0;
        float _S951 = 1.0f - _S949.primal_0 * _S949.primal_0 / 24.0f;
        float _S952 = ((0.0f - (_S950 + _S950) * 0.0416666679084301f) * _S906 - _S951) / (_S906 * _S906);
        k_18 = _S951 / _S906;
        s_diff_k_7 = _S952;
    }
    else
    {
        float _S953 = _S949.differential_0 * 0.5f;
        DiffPair_float_0 _S954;
        (&_S954)->primal_0 = 0.5f * _S949.primal_0;
        (&_S954)->differential_0 = _S953;
        DiffPair_float_0 _S955 = _d_sin_0(&_S954);
        float _S956 = 2.0f * _S955.primal_0;
        float _S957 = (_S955.differential_0 * 2.0f * _S946.primal_0 - _S956 * _S946.differential_0) / (_S946.primal_0 * _S946.primal_0);
        k_18 = _S956 / _S946.primal_0;
        s_diff_k_7 = _S957;
    }
    float2  _S958 = _S902 * make_float2 (k_18);
    float2  _S959 = _S944 * make_float2 (k_18) + make_float2 (s_diff_k_7) * _S902;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S960;
    (&_S960)->primal_0 = _S958;
    (&_S960)->differential_0 = _S959;
    FixedArray<float, 8>  _S961 = dist_coeffs_43;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S962 = s_fwd_DistRational_distort_0(&_S960, &_S961);
    float _S963 = _S962.differential_0.y * fy_12;
    *&(((&J_12)->rows + (int(0)))->z) = _S962.differential_0.x * fx_12;
    *&(((&J_12)->rows + (int(1)))->z) = _S963;
    return J_12;
}

inline __device__ float2  distort_point_rational(float2  uv_41, int camera_model_12, FixedArray<float, 8>  dist_coeffs_44)
{
    float2  _S964;
    for(;;)
    {
        if(camera_model_12 == int(3))
        {
            _S964 = uv_41;
            break;
        }
        float k_19;
        if(camera_model_12 == int(1))
        {
            float r_38 = length_0(uv_41);
            float theta_14 = (F32_atan((r_38)));
            if(r_38 < 0.00100000004749745f)
            {
                k_19 = 1.0f - theta_14 * theta_14 / 6.0f;
            }
            else
            {
                k_19 = theta_14 / r_38;
            }
            _S964 = uv_41 * make_float2 (k_19);
        }
        else
        {
            if(camera_model_12 == int(2))
            {
                float r_39 = length_0(uv_41);
                float theta_15 = (F32_atan((r_39)));
                if(r_39 < 0.00100000004749745f)
                {
                    k_19 = 1.0f - theta_15 * theta_15 / 24.0f;
                }
                else
                {
                    k_19 = 2.0f * (F32_sin((0.5f * theta_15))) / r_39;
                }
                _S964 = uv_41 * make_float2 (k_19);
            }
            else
            {
                _S964 = uv_41;
            }
        }
        FixedArray<float, 8>  _S965 = dist_coeffs_44;
        float2  _S966 = DistRational_distort_0(_S964, &_S965);
        _S964 = _S966;
        break;
    }
    return _S964;
}

inline __device__ bool undistort_point_rational(float2  uv_42, int camera_model_13, FixedArray<float, 8>  dist_coeffs_45, float2  * uv_undist_7)
{
    bool _S967;
    for(;;)
    {
        *uv_undist_7 = make_float2 (0.0f);
        if(camera_model_13 == int(3))
        {
            float lon_6 = uv_42.x;
            float lat_9 = uv_42.y;
            float cl_9 = (F32_cos((lat_9)));
            *uv_undist_7 = make_float2 (cl_9 * (F32_sin((lon_6))), (F32_sin((lat_9)))) / make_float2 ((F32_max((cl_9 * (F32_cos((lon_6)))), (9.999999960041972e-13f))));
            _S967 = true;
            break;
        }
        FixedArray<float, 8>  _S968 = dist_coeffs_45;
        float2  uv_u_9;
        bool _S969 = undistort_point_3(uv_42, &_S968, int(8), &uv_u_9);
        if(!_S969)
        {
            _S967 = false;
            break;
        }
        float2  _S970 = uv_u_9;
        float3  raydir_9;
        if(camera_model_13 == int(1))
        {
            float r_40 = length_0(_S970);
            float s_9;
            if(r_40 < 0.00100000004749745f)
            {
                s_9 = 1.0f - r_40 * r_40 / 6.0f;
            }
            else
            {
                s_9 = (F32_sin((r_40))) / r_40;
            }
            raydir_9 = make_float3 ((_S970 * make_float2 (s_9)).x, (_S970 * make_float2 (s_9)).y, (F32_cos((r_40))));
        }
        else
        {
            if(camera_model_13 == int(2))
            {
                float r_41 = length_0(_S970);
                raydir_9 = make_float3 ((_S970 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_41 * r_41)))))))).x, (_S970 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_41 * r_41)))))))).y, 1.0f - 0.5f * r_41 * r_41);
            }
            else
            {
                raydir_9 = make_float3 (_S970.x, _S970.y, 1.0f);
            }
        }
        *uv_undist_7 = float2 {raydir_9.x, raydir_9.y} / make_float2 ((F32_max((raydir_9.z), (9.999999960041972e-13f))));
        _S967 = true;
        break;
    }
    return _S967;
}

inline __device__ bool unproject_point_rational(float2  uv_43, int camera_model_14, FixedArray<float, 8>  dist_coeffs_46, float3  * raydir_10)
{
    bool _S971;
    for(;;)
    {
        int3  _S972 = make_int3 (int(0));
        float3  _S973 = make_float3 ((float)_S972.x, (float)_S972.y, (float)_S972.z);
        *raydir_10 = _S973;
        if(camera_model_14 == int(3))
        {
            float lon_7 = uv_43.x;
            float lat_10 = uv_43.y;
            float cl_10 = (F32_cos((lat_10)));
            *raydir_10 = make_float3 (cl_10 * (F32_sin((lon_7))), (F32_sin((lat_10))), cl_10 * (F32_cos((lon_7))));
            _S971 = true;
            break;
        }
        FixedArray<float, 8>  _S974 = dist_coeffs_46;
        float2  uv_u_10;
        bool _S975 = undistort_point_3(uv_43, &_S974, int(8), &uv_u_10);
        if(!_S975)
        {
            _S971 = false;
            break;
        }
        float2  _S976 = uv_u_10;
        if(camera_model_14 == int(1))
        {
            float r_42 = length_0(_S976);
            float s_10;
            if(r_42 < 0.00100000004749745f)
            {
                s_10 = 1.0f - r_42 * r_42 / 6.0f;
            }
            else
            {
                s_10 = (F32_sin((r_42))) / r_42;
            }
            *raydir_10 = make_float3 ((_S976 * make_float2 (s_10)).x, (_S976 * make_float2 (s_10)).y, (F32_cos((r_42))));
        }
        else
        {
            if(camera_model_14 == int(2))
            {
                float r_43 = length_0(_S976);
                *raydir_10 = make_float3 ((_S976 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_43 * r_43)))))))).x, (_S976 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_43 * r_43)))))))).y, 1.0f - 0.5f * r_43 * r_43);
            }
            else
            {
                *raydir_10 = make_float3 (_S976.x, _S976.y, 1.0f);
            }
        }
        _S971 = true;
        break;
    }
    return _S971;
}

inline __device__ bool generate_ray_rational(float2  uv_44, int camera_model_15, FixedArray<float, 8>  dist_coeffs_47, float3  * raydir_11)
{
    bool _S977;
    for(;;)
    {
        if(camera_model_15 == int(3))
        {
            float _S978 = uv_44.x;
            if((F32_abs((_S978))) > 3.14159274101257324f)
            {
                _S977 = true;
            }
            else
            {
                _S977 = (F32_abs((uv_44.y))) > 1.57079637050628662f;
            }
            if(_S977)
            {
                int3  _S979 = make_int3 (int(0));
                float3  _S980 = make_float3 ((float)_S979.x, (float)_S979.y, (float)_S979.z);
                *raydir_11 = _S980;
                _S977 = false;
                break;
            }
            float lat_11 = uv_44.y;
            float cl_11 = (F32_cos((lat_11)));
            *raydir_11 = make_float3 (cl_11 * (F32_sin((_S978))), (F32_sin((lat_11))), cl_11 * (F32_cos((_S978))));
            _S977 = true;
            break;
        }
        FixedArray<float, 8>  _S981 = dist_coeffs_47;
        float2  uv_u_11;
        bool _S982 = undistort_point_3(uv_44, &_S981, int(8), &uv_u_11);
        if(!_S982)
        {
            int3  _S983 = make_int3 (int(0));
            float3  _S984 = make_float3 ((float)_S983.x, (float)_S983.y, (float)_S983.z);
            *raydir_11 = _S984;
            _S977 = false;
            break;
        }
        float2  _S985 = uv_u_11;
        if(camera_model_15 == int(1))
        {
            float r_44 = length_0(_S985);
            if(r_44 >= 3.14159274101257324f)
            {
                int3  _S986 = make_int3 (int(0));
                float3  _S987 = make_float3 ((float)_S986.x, (float)_S986.y, (float)_S986.z);
                *raydir_11 = _S987;
                _S977 = false;
                break;
            }
            float s_11;
            if(r_44 < 0.00100000004749745f)
            {
                s_11 = 1.0f - r_44 * r_44 / 6.0f;
            }
            else
            {
                s_11 = (F32_sin((r_44))) / r_44;
            }
            *raydir_11 = make_float3 ((_S985 * make_float2 (s_11)).x, (_S985 * make_float2 (s_11)).y, (F32_cos((r_44))));
        }
        else
        {
            if(camera_model_15 == int(2))
            {
                float r_45 = length_0(_S985);
                if(r_45 >= 2.0f)
                {
                    int3  _S988 = make_int3 (int(0));
                    float3  _S989 = make_float3 ((float)_S988.x, (float)_S988.y, (float)_S988.z);
                    *raydir_11 = _S989;
                    _S977 = false;
                    break;
                }
                *raydir_11 = make_float3 ((_S985 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_45 * r_45)))))))).x, (_S985 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_45 * r_45)))))))).y, 1.0f - 0.5f * r_45 * r_45);
            }
            else
            {
                *raydir_11 = make_float3 (_S985.x, _S985.y, 1.0f);
            }
        }
        *raydir_11 = normalize_0(*raydir_11);
        _S977 = true;
        break;
    }
    return _S977;
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_4, float3  dOut_3)
{
    float _S990 = (*right_4).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  right_d_result_2;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = (*left_4).primal_0.x * dOut_3.x;
    float sum_10 = _S990 + (*right_4).primal_0.rows[int(0)].y * dOut_3.y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = (*left_4).primal_0.x * dOut_3.y;
    float sum_11 = sum_10 + (*right_4).primal_0.rows[int(0)].z * dOut_3.z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = (*left_4).primal_0.x * dOut_3.z;
    float3  left_d_result_2;
    *&((&left_d_result_2)->x) = sum_11;
    float _S991 = (*right_4).primal_0.rows[int(1)].x * dOut_3.x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = (*left_4).primal_0.y * dOut_3.x;
    float sum_12 = _S991 + (*right_4).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = (*left_4).primal_0.y * dOut_3.y;
    float sum_13 = sum_12 + (*right_4).primal_0.rows[int(1)].z * dOut_3.z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = (*left_4).primal_0.y * dOut_3.z;
    *&((&left_d_result_2)->y) = sum_13;
    float _S992 = (*right_4).primal_0.rows[int(2)].x * dOut_3.x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = (*left_4).primal_0.z * dOut_3.x;
    float sum_14 = _S992 + (*right_4).primal_0.rows[int(2)].y * dOut_3.y;
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
        int i_7 = int(0);
        float sum_16 = 0.0f;
        for(;;)
        {
            if(i_7 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_17 = sum_16 + _slang_vector_get_element(left_5, i_7) * _slang_vector_get_element(right_5.rows[i_7], j_1);
            i_7 = i_7 + int(1);
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

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_1, float3  raydir_12)
{
    return mul_3(raydir_12, R_1);
}

inline __device__ float3  undo_transform_ray_d(Matrix<float, 3, 3>  R_2, float3  raydir_13)
{
    return mul_3(raydir_13, transpose_0(R_2));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S993, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S994, float3  _S995)
{
    _d_mul_1(_S993, _S994, _S995);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S996 = - _s_dOut_0;
    float3  _S997 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S998;
    (&_S998)->primal_0 = (*dpt_0).primal_0;
    (&_S998)->differential_0 = _S997;
    Matrix<float, 3, 3>  _S999 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1000;
    (&_S1000)->primal_0 = (*dpR_0).primal_0;
    (&_S1000)->differential_0 = _S999;
    s_bwd_prop_mul_0(&_S998, &_S1000, _S996);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S998.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S1000.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1001, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1002, float3  _S1003)
{
    s_bwd_prop_transform_ray_o_0(_S1001, _S1002, _S1003);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S1004 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S1004;
    float3  _S1005 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S1005;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S1006 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1007;
    (&_S1007)->primal_0 = (*dpraydir_0).primal_0;
    (&_S1007)->differential_0 = _S1006;
    Matrix<float, 3, 3>  _S1008 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1009;
    (&_S1009)->primal_0 = (*dpR_1).primal_0;
    (&_S1009)->differential_0 = _S1008;
    s_bwd_prop_mul_0(&_S1007, &_S1009, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S1007.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S1009.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1010, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1011, float3  _S1012)
{
    s_bwd_prop_transform_ray_d_0(_S1010, _S1011, _S1012);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_14, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S1013 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S1013;
    float3  _S1014 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_14;
    (&dp_raydir_0)->differential_0 = _S1014;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_4)
{
    float _S1015 = (F32_exp(((*dpx_5).primal_0))) * dOut_4;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S1015;
    return;
}

inline __device__ float3  exp_0(float3  x_12)
{
    float3  result_9;
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
        *_slang_vector_get_element_ptr(&result_9, i_8) = (F32_exp((_slang_vector_get_element(x_12, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_9;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  dOut_5)
{
    float3  _S1016 = exp_0((*dpx_6).primal_0) * dOut_5;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S1016;
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
    float3  _S1017 = exp_0(- scale_4);
    return mul_1(makeMatrix<float, 3, 3> (_S1017.x, 0.0f, 0.0f, 0.0f, _S1017.y, 0.0f, 0.0f, 0.0f, _S1017.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)))));
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ float3  s_primal_ctx_exp_0(float3  _S1018)
{
    return exp_0(_S1018);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1019, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1020, Matrix<float, 3, 3>  _S1021)
{
    mul_0(_S1019, _S1020, _S1021);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1022, float3  _S1023)
{
    _d_exp_vector_0(_S1022, _S1023);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_2)
{
    float _S1024 = (*dpquat_0).primal_0.y;
    float x2_6 = _S1024 * _S1024;
    float y2_6 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_6 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_6 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_6 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_6 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_6 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    float3  _S1025 = - (*dpscale_0).primal_0;
    float3  _S1026 = s_primal_ctx_exp_0(_S1025);
    Matrix<float, 3, 3>  _S1027 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))));
    Matrix<float, 3, 3>  _S1028 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1029;
    (&_S1029)->primal_0 = makeMatrix<float, 3, 3> (_S1026.x, 0.0f, 0.0f, 0.0f, _S1026.y, 0.0f, 0.0f, 0.0f, _S1026.z);
    (&_S1029)->differential_0 = _S1028;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1030;
    (&_S1030)->primal_0 = _S1027;
    (&_S1030)->differential_0 = _S1028;
    s_bwd_prop_mul_1(&_S1029, &_S1030, _s_dOut_2);
    Matrix<float, 3, 3>  _S1031 = transpose_0(_S1030.differential_0);
    float3  _S1032 = make_float3 (_S1029.differential_0.rows[int(0)].x, _S1029.differential_0.rows[int(1)].y, _S1029.differential_0.rows[int(2)].z);
    float3  _S1033 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1034;
    (&_S1034)->primal_0 = _S1025;
    (&_S1034)->differential_0 = _S1033;
    s_bwd_prop_exp_0(&_S1034, _S1032);
    float3  _S1035 = - _S1034.differential_0;
    Matrix<float, 3, 3>  _S1036 = transpose_0(_S1031);
    float _S1037 = 2.0f * - _S1036.rows[int(2)].z;
    float _S1038 = 2.0f * _S1036.rows[int(2)].y;
    float _S1039 = 2.0f * _S1036.rows[int(2)].x;
    float _S1040 = 2.0f * _S1036.rows[int(1)].z;
    float _S1041 = 2.0f * - _S1036.rows[int(1)].y;
    float _S1042 = 2.0f * _S1036.rows[int(1)].x;
    float _S1043 = 2.0f * _S1036.rows[int(0)].z;
    float _S1044 = 2.0f * _S1036.rows[int(0)].y;
    float _S1045 = 2.0f * - _S1036.rows[int(0)].x;
    float _S1046 = - _S1042 + _S1044;
    float _S1047 = _S1039 + - _S1043;
    float _S1048 = - _S1038 + _S1040;
    float _S1049 = _S1038 + _S1040;
    float _S1050 = _S1039 + _S1043;
    float _S1051 = _S1042 + _S1044;
    float _S1052 = (*dpquat_0).primal_0.w * (_S1041 + _S1045);
    float _S1053 = (*dpquat_0).primal_0.z * (_S1037 + _S1045);
    float _S1054 = (*dpquat_0).primal_0.y * (_S1037 + _S1041);
    float _S1055 = (*dpquat_0).primal_0.x * _S1046 + (*dpquat_0).primal_0.z * _S1049 + (*dpquat_0).primal_0.y * _S1050 + _S1052 + _S1052;
    float _S1056 = (*dpquat_0).primal_0.x * _S1047 + (*dpquat_0).primal_0.w * _S1049 + (*dpquat_0).primal_0.y * _S1051 + _S1053 + _S1053;
    float _S1057 = (*dpquat_0).primal_0.x * _S1048 + (*dpquat_0).primal_0.w * _S1050 + (*dpquat_0).primal_0.z * _S1051 + _S1054 + _S1054;
    float _S1058 = (*dpquat_0).primal_0.w * _S1046 + (*dpquat_0).primal_0.z * _S1047 + (*dpquat_0).primal_0.y * _S1048;
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S1035;
    float4  _S1059 = make_float4 (0.0f);
    *&((&_S1059)->w) = _S1055;
    *&((&_S1059)->z) = _S1056;
    *&((&_S1059)->y) = _S1057;
    *&((&_S1059)->x) = _S1058;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S1059;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1060, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1061, Matrix<float, 3, 3>  _S1062)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S1060, _S1061, _S1062);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_6, float3  scale_5, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_0, float3  * v_scale_0)
{
    float4  _S1063 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_6;
    (&dp_quat_0)->differential_0 = _S1063;
    float3  _S1064 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_5;
    (&dp_scale_0)->differential_0 = _S1064;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_0 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_6)
{
    float _S1065 = dOut_6.y;
    float _S1066 = dOut_6.z;
    float _S1067 = dOut_6.x;
    float _S1068 = (*a_0).primal_0.z * _S1065 + - (*a_0).primal_0.y * _S1066;
    float _S1069 = - (*a_0).primal_0.z * _S1067 + (*a_0).primal_0.x * _S1066;
    float _S1070 = (*a_0).primal_0.y * _S1067 + - (*a_0).primal_0.x * _S1065;
    float3  _S1071 = make_float3 (- (*b_0).primal_0.z * _S1065 + (*b_0).primal_0.y * _S1066, (*b_0).primal_0.z * _S1067 + - (*b_0).primal_0.x * _S1066, - (*b_0).primal_0.y * _S1067 + (*b_0).primal_0.x * _S1065);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S1071;
    float3  _S1072 = make_float3 (_S1068, _S1069, _S1070);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S1072;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S1073 = left_6.y;
    float _S1074 = right_6.z;
    float _S1075 = left_6.z;
    float _S1076 = right_6.y;
    float _S1077 = right_6.x;
    float _S1078 = left_6.x;
    return make_float3 (_S1073 * _S1074 - _S1075 * _S1076, _S1075 * _S1077 - _S1078 * _S1074, _S1078 * _S1076 - _S1073 * _S1077);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_0, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_2(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_2(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_0 * (F32_exp((-0.5f * dot_1(gcrod_0, gcrod_0) / dot_1(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S1079, float3  _S1080)
{
    return mul_2(_S1079, _S1080);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1081, float3  _S1082)
{
    return cross_0(_S1081, _S1082);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S1083, float3  _S1084)
{
    return dot_1(_S1083, _S1084);
}

inline __device__ float s_primal_ctx_exp_1(float _S1085)
{
    return (F32_exp((_S1085)));
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_float_0 * _S1086, float _S1087)
{
    _d_exp_0(_S1086, _S1087);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1088, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1089, float _S1090)
{
    _d_dot_0(_S1088, _S1089, _S1090);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1091, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1092, float3  _S1093)
{
    _d_cross_0(_S1091, _S1092, _S1093);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1094, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1095, float3  _S1096)
{
    _d_mul_0(_S1094, _S1095, _S1096);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    float3  _S1097 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S1098 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S1097);
    float3  _S1099 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S1100 = s_primal_ctx_cross_0(_S1099, _S1098);
    float _S1101 = -0.5f * s_primal_ctx_dot_0(_S1100, _S1100);
    float _S1102 = s_primal_ctx_dot_0(_S1099, _S1099);
    float _S1103 = _S1101 / _S1102;
    float _S1104 = _S1102 * _S1102;
    float _S1105 = (*dpopacity_0).primal_0 * _s_dOut_3;
    float _S1106 = s_primal_ctx_exp_1(_S1103) * _s_dOut_3;
    DiffPair_float_0 _S1107;
    (&_S1107)->primal_0 = _S1103;
    (&_S1107)->differential_0 = 0.0f;
    s_bwd_prop_exp_1(&_S1107, _S1105);
    float _S1108 = _S1107.differential_0 / _S1104;
    float _S1109 = _S1101 * - _S1108;
    float _S1110 = _S1102 * _S1108;
    float3  _S1111 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1112;
    (&_S1112)->primal_0 = _S1099;
    (&_S1112)->differential_0 = _S1111;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1113;
    (&_S1113)->primal_0 = _S1099;
    (&_S1113)->differential_0 = _S1111;
    s_bwd_prop_dot_0(&_S1112, &_S1113, _S1109);
    float _S1114 = -0.5f * _S1110;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1115;
    (&_S1115)->primal_0 = _S1100;
    (&_S1115)->differential_0 = _S1111;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1116;
    (&_S1116)->primal_0 = _S1100;
    (&_S1116)->differential_0 = _S1111;
    s_bwd_prop_dot_0(&_S1115, &_S1116, _S1114);
    float3  _S1117 = _S1116.differential_0 + _S1115.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1118;
    (&_S1118)->primal_0 = _S1099;
    (&_S1118)->differential_0 = _S1111;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1119;
    (&_S1119)->primal_0 = _S1098;
    (&_S1119)->differential_0 = _S1111;
    s_bwd_prop_cross_0(&_S1118, &_S1119, _S1117);
    float3  _S1120 = _S1113.differential_0 + _S1112.differential_0 + _S1118.differential_0;
    Matrix<float, 3, 3>  _S1121 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1122;
    (&_S1122)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S1122)->differential_0 = _S1121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1123;
    (&_S1123)->primal_0 = (*dpray_d_0).primal_0;
    (&_S1123)->differential_0 = _S1111;
    s_bwd_prop_mul_2(&_S1122, &_S1123, _S1120);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1124;
    (&_S1124)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S1124)->differential_0 = _S1121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1125;
    (&_S1125)->primal_0 = _S1097;
    (&_S1125)->differential_0 = _S1111;
    s_bwd_prop_mul_2(&_S1124, &_S1125, _S1119.differential_0);
    float3  _S1126 = - _S1125.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1123.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1125.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S1106;
    Matrix<float, 3, 3>  _S1127 = _S1122.differential_0 + _S1124.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S1127;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S1126;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1128, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1129, DiffPair_float_0 * _S1130, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1131, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1132, float _S1133)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S1128, _S1129, _S1130, _S1131, _S1132, _S1133);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_1, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_0, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1134 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S1134;
    Matrix<float, 3, 3>  _S1135 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S1135;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_1;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S1134;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S1134;
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
    *depth_0 = - dot_1(mul_2(iscl_rot_2, ray_o_2 - mean_2), grd_1) / dot_1(grd_1, grd_1);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S1136 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S1137 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S1136);
    float3  _S1138 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S1139 = s_primal_ctx_dot_0(_S1138, _S1138);
    float _S1140 = dpdepth_0 / (_S1139 * _S1139);
    float _S1141 = - s_primal_ctx_dot_0(_S1137, _S1138) * - _S1140;
    float _S1142 = _S1139 * _S1140;
    float3  _S1143 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1144;
    (&_S1144)->primal_0 = _S1138;
    (&_S1144)->differential_0 = _S1143;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1145;
    (&_S1145)->primal_0 = _S1138;
    (&_S1145)->differential_0 = _S1143;
    s_bwd_prop_dot_0(&_S1144, &_S1145, _S1141);
    float _S1146 = - _S1142;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1147;
    (&_S1147)->primal_0 = _S1137;
    (&_S1147)->differential_0 = _S1143;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1148;
    (&_S1148)->primal_0 = _S1138;
    (&_S1148)->differential_0 = _S1143;
    s_bwd_prop_dot_0(&_S1147, &_S1148, _S1146);
    float3  _S1149 = _S1145.differential_0 + _S1144.differential_0 + _S1148.differential_0;
    Matrix<float, 3, 3>  _S1150 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1151;
    (&_S1151)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S1151)->differential_0 = _S1150;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1152;
    (&_S1152)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1152)->differential_0 = _S1143;
    s_bwd_prop_mul_2(&_S1151, &_S1152, _S1149);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1153;
    (&_S1153)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S1153)->differential_0 = _S1150;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1154;
    (&_S1154)->primal_0 = _S1136;
    (&_S1154)->differential_0 = _S1143;
    s_bwd_prop_mul_2(&_S1153, &_S1154, _S1147.differential_0);
    float3  _S1155 = - _S1154.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1152.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1154.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S1156 = _S1151.differential_0 + _S1153.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S1156;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S1155;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1157, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1158, DiffPair_float_0 * _S1159, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1160, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1161, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1162, float3  _S1163, float _S1164)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S1157, _S1158, _S1159, _S1160, _S1161, _S1162, _S1163, _S1164);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_3, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_1, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S1165 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S1165;
    Matrix<float, 3, 3>  _S1166 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S1166;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_3;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S1165;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S1165;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S1165;
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

