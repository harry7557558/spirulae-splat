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

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_0, FixedArray<float, 10>  dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float _S1 = 0.0f * v_0;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float s_diff_r2_0 = u_0 + u_0 + (_S1 + _S1);
    float _S2 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
    float _S3 = dist_coeffs_0[int(1)] + r2_0 * _S2;
    float _S4 = dist_coeffs_0[int(0)] + r2_0 * _S3;
    float radial_0 = 1.0f + r2_0 * _S4;
    float _S5 = 2.0f * dist_coeffs_0[int(4)];
    float _S6 = _S5 * u_0;
    float _S7 = 2.0f * u_0;
    float _S8 = 2.0f * dist_coeffs_0[int(5)];
    float _S9 = _S8 * u_0;
    float _S10 = 2.0f * v_0;
    float2  _S11 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S4 + (s_diff_r2_0 * _S3 + (s_diff_r2_0 * _S2 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (_S5 * v_0 + 0.0f * _S6 + (s_diff_r2_0 + (_S7 + _S7)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S8 * v_0 + 0.0f * _S9 + (s_diff_r2_0 + (_S1 + 0.0f * _S10)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
    float _S12 = 0.0f * u_0;
    float s_diff_r2_1 = _S12 + _S12 + (v_0 + v_0);
    float2  _S13 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S4 + (s_diff_r2_1 * _S3 + (s_diff_r2_1 * _S2 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (0.0f * _S5 * v_0 + _S6 + (s_diff_r2_1 + (_S12 + 0.0f * _S7)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S8 * v_0 + _S9 + (s_diff_r2_1 + (_S10 + _S10)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
    Matrix<float, 2, 2>  _S14 = transpose_1(makeMatrix<float, 2, 2> (_S11 + make_float2 (_S11.x * dist_coeffs_0[int(8)] + _S11.y * dist_coeffs_0[int(9)], 0.0f), _S13 + make_float2 (_S13.x * dist_coeffs_0[int(8)] + _S13.y * dist_coeffs_0[int(9)], 0.0f)));
    return (F32_min((determinant_0(_S14)), ((F32_min((_S14.rows[int(0)].x), (_S14.rows[int(1)].y)))))) > 0.0f;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_0, float dOut_1)
{
    float _S15 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_0).primal_0)))))) * dOut_1;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S15;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

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

inline __device__ float dot_0(float3  x_4, float3  y_0)
{
    int i_1 = int(0);
    float result_3 = 0.0f;
    for(;;)
    {
        if(i_1 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_4 = result_3 + _slang_vector_get_element(x_4, i_1) * _slang_vector_get_element(y_0, i_1);
        i_1 = i_1 + int(1);
        result_3 = result_4;
    }
    return result_3;
}

inline __device__ float dot_1(float2  x_5, float2  y_1)
{
    int i_2 = int(0);
    float result_5 = 0.0f;
    for(;;)
    {
        if(i_2 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_6 = result_5 + _slang_vector_get_element(x_5, i_2) * _slang_vector_get_element(y_1, i_2);
        i_2 = i_2 + int(1);
        result_5 = result_6;
    }
    return result_5;
}

inline __device__ float dot_2(float4  x_6, float4  y_2)
{
    int i_3 = int(0);
    float result_7 = 0.0f;
    for(;;)
    {
        if(i_3 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_8 = result_7 + _slang_vector_get_element(x_6, i_3) * _slang_vector_get_element(y_2, i_3);
        i_3 = i_3 + int(1);
        result_7 = result_8;
    }
    return result_7;
}

inline __device__ float length_0(float2  x_7)
{
    return (F32_sqrt((dot_1(x_7, x_7))));
}

inline __device__ float length_1(float3  x_8)
{
    return (F32_sqrt((dot_0(x_8, x_8))));
}

inline __device__ float length_2(float4  x_9)
{
    return (F32_sqrt((dot_2(x_9, x_9))));
}

inline __device__ float2  distort_point(float2  uv_1, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_1)
{
    float2  _S16;
    if(is_fisheye_0)
    {
        float r_3 = length_0(uv_1);
        float theta_0 = (F32_atan((r_3)));
        float _S17;
        if(r_3 < 0.00100000004749745f)
        {
            _S17 = 1.0f - r_3 * r_3 / 3.0f;
        }
        else
        {
            _S17 = theta_0 / r_3;
        }
        _S16 = uv_1 * make_float2 (_S17);
    }
    else
    {
        _S16 = uv_1;
    }
    float u_1 = _S16.x;
    float v_1 = _S16.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S18 = _S16 * make_float2 (1.0f + r2_1 * (dist_coeffs_1[int(0)] + r2_1 * (dist_coeffs_1[int(1)] + r2_1 * (dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_1 * v_1 + dist_coeffs_1[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_1[int(6)] * r2_1, 2.0f * dist_coeffs_1[int(5)] * u_1 * v_1 + dist_coeffs_1[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_1[int(7)] * r2_1);
    return _S18 + make_float2 (dist_coeffs_1[int(8)] * _S18.x + dist_coeffs_1[int(9)] * _S18.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_2, FixedArray<float, 10>  * dist_coeffs_2, int maxiter_0, float2  * uv_undist_0)
{
    int i_4 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_4 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S19 = (*dist_coeffs_2)[int(3)];
        float _S20 = (*dist_coeffs_2)[int(4)];
        float _S21 = (*dist_coeffs_2)[int(5)];
        float _S22 = (*dist_coeffs_2)[int(6)];
        float _S23 = (*dist_coeffs_2)[int(7)];
        float _S24 = (*dist_coeffs_2)[int(8)];
        float _S25 = (*dist_coeffs_2)[int(9)];
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float _S26 = (*dist_coeffs_2)[int(2)] + r2_2 * (*dist_coeffs_2)[int(3)];
        float _S27 = (*dist_coeffs_2)[int(1)] + r2_2 * _S26;
        float _S28 = (*dist_coeffs_2)[int(0)] + r2_2 * _S27;
        float radial_1 = 1.0f + r2_2 * _S28;
        float _S29 = 2.0f * (*dist_coeffs_2)[int(4)];
        float _S30 = _S29 * u_2;
        float _S31 = 2.0f * u_2;
        float _S32 = 2.0f * (*dist_coeffs_2)[int(5)];
        float _S33 = _S32 * u_2;
        float _S34 = 2.0f * v_2;
        float2  _S35 = q_0 * make_float2 (radial_1) + make_float2 (_S30 * v_2 + (*dist_coeffs_2)[int(5)] * (r2_2 + _S31 * u_2) + (*dist_coeffs_2)[int(6)] * r2_2, _S33 * v_2 + (*dist_coeffs_2)[int(4)] * (r2_2 + _S34 * v_2) + (*dist_coeffs_2)[int(7)] * r2_2);
        float2  r_4 = _S35 + make_float2 ((*dist_coeffs_2)[int(8)] * _S35.x + (*dist_coeffs_2)[int(9)] * _S35.y, 0.0f) - uv_2;
        float _S36 = 0.0f * v_2;
        float s_diff_r2_2 = u_2 + u_2 + (_S36 + _S36);
        float2  _S37 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S28 + (s_diff_r2_2 * _S27 + (s_diff_r2_2 * _S26 + s_diff_r2_2 * _S19 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (_S29 * v_2 + 0.0f * _S30 + (s_diff_r2_2 + (_S31 + _S31)) * _S21 + s_diff_r2_2 * _S22, _S32 * v_2 + 0.0f * _S33 + (s_diff_r2_2 + (_S36 + 0.0f * _S34)) * _S20 + s_diff_r2_2 * _S23);
        float _S38 = 0.0f * u_2;
        float s_diff_r2_3 = _S38 + _S38 + (v_2 + v_2);
        float2  _S39 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S28 + (s_diff_r2_3 * _S27 + (s_diff_r2_3 * _S26 + s_diff_r2_3 * _S19 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (0.0f * _S29 * v_2 + _S30 + (s_diff_r2_3 + (_S38 + 0.0f * _S31)) * _S21 + s_diff_r2_3 * _S22, 0.0f * _S32 * v_2 + _S33 + (s_diff_r2_3 + (_S34 + _S34)) * _S20 + s_diff_r2_3 * _S23);
        Matrix<float, 2, 2>  _S40 = transpose_1(makeMatrix<float, 2, 2> (_S37 + make_float2 (_S37.x * _S24 + _S37.y * _S25, 0.0f), _S39 + make_float2 (_S39.x * _S24 + _S39.y * _S25, 0.0f)));
        float inv_det_0 = 1.0f / (_S40.rows[int(0)].x * _S40.rows[int(1)].y - _S40.rows[int(0)].y * _S40.rows[int(1)].x);
        float _S41 = r_4.x;
        float _S42 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S41 * _S40.rows[int(1)].y - _S42 * _S40.rows[int(0)].y) * inv_det_0, (- _S41 * _S40.rows[int(1)].x + _S42 * _S40.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S43 = (*dist_coeffs_2)[int(0)];
    float _S44 = (*dist_coeffs_2)[int(1)];
    float _S45 = (*dist_coeffs_2)[int(2)];
    float _S46 = (*dist_coeffs_2)[int(3)];
    float _S47 = (*dist_coeffs_2)[int(4)];
    float _S48 = (*dist_coeffs_2)[int(5)];
    float _S49 = (*dist_coeffs_2)[int(6)];
    float _S50 = (*dist_coeffs_2)[int(7)];
    float _S51 = (*dist_coeffs_2)[int(8)];
    float _S52 = (*dist_coeffs_2)[int(9)];
    float u_3 = q_0.x;
    float v_3 = q_0.y;
    float _S53 = 0.0f * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = u_3 + u_3 + (_S53 + _S53);
    float _S54 = (*dist_coeffs_2)[int(2)] + r2_3 * (*dist_coeffs_2)[int(3)];
    float _S55 = (*dist_coeffs_2)[int(1)] + r2_3 * _S54;
    float _S56 = (*dist_coeffs_2)[int(0)] + r2_3 * _S55;
    float radial_2 = 1.0f + r2_3 * _S56;
    float _S57 = 2.0f * (*dist_coeffs_2)[int(4)];
    float _S58 = _S57 * u_3;
    float _S59 = 2.0f * u_3;
    float _S60 = 2.0f * (*dist_coeffs_2)[int(5)];
    float _S61 = _S60 * u_3;
    float _S62 = 2.0f * v_3;
    float2  _S63 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S56 + (s_diff_r2_4 * _S55 + (s_diff_r2_4 * _S54 + s_diff_r2_4 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (_S57 * v_3 + 0.0f * _S58 + (s_diff_r2_4 + (_S59 + _S59)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_4 * (*dist_coeffs_2)[int(6)], _S60 * v_3 + 0.0f * _S61 + (s_diff_r2_4 + (_S53 + 0.0f * _S62)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_4 * (*dist_coeffs_2)[int(7)]);
    float _S64 = 0.0f * u_3;
    float s_diff_r2_5 = _S64 + _S64 + (v_3 + v_3);
    float2  _S65 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S56 + (s_diff_r2_5 * _S55 + (s_diff_r2_5 * _S54 + s_diff_r2_5 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (0.0f * _S57 * v_3 + _S58 + (s_diff_r2_5 + (_S64 + 0.0f * _S59)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_5 * (*dist_coeffs_2)[int(6)], 0.0f * _S60 * v_3 + _S61 + (s_diff_r2_5 + (_S62 + _S62)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_5 * (*dist_coeffs_2)[int(7)]);
    Matrix<float, 2, 2>  _S66 = transpose_1(makeMatrix<float, 2, 2> (_S63 + make_float2 (_S63.x * (*dist_coeffs_2)[int(8)] + _S63.y * (*dist_coeffs_2)[int(9)], 0.0f), _S65 + make_float2 (_S65.x * (*dist_coeffs_2)[int(8)] + _S65.y * (*dist_coeffs_2)[int(9)], 0.0f)));
    bool _S67;
    if((F32_min((determinant_0(_S66)), ((F32_min((_S66.rows[int(0)].x), (_S66.rows[int(1)].y)))))) > 0.0f)
    {
        float u_4 = (*uv_undist_0).x;
        float v_4 = (*uv_undist_0).y;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float2  _S68 = *uv_undist_0 * make_float2 (1.0f + r2_4 * (_S43 + r2_4 * (_S44 + r2_4 * (_S45 + r2_4 * _S46)))) + make_float2 (_S57 * u_4 * v_4 + _S48 * (r2_4 + 2.0f * u_4 * u_4) + _S49 * r2_4, _S60 * u_4 * v_4 + _S47 * (r2_4 + 2.0f * v_4 * v_4) + _S50 * r2_4);
        _S67 = (length_0(_S68 + make_float2 (_S51 * _S68.x + _S52 * _S68.y, 0.0f) - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S67 = false;
    }
    return _S67;
}

inline __device__ bool undistort_point(float2  uv_3, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_undist_1)
{
    float2  _S69 = uv_3;
    FixedArray<float, 10>  _S70 = dist_coeffs_3;
    bool _S71 = undistort_point_0(uv_3, &_S70, int(8), &_S69);
    if(!_S71)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S72 = _S69;
        float theta_1 = length_0(_S69);
        float _S73;
        if(theta_1 < 0.00100000004749745f)
        {
            _S73 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S73 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S74 = make_float3 ((_S72 * make_float2 (_S73)).x, (_S72 * make_float2 (_S73)).y, (F32_cos((theta_1))));
        raydir_0 = _S74;
    }
    else
    {
        raydir_0 = make_float3 (_S69.x, _S69.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_4, float3  * raydir_1)
{
    float2  _S75 = uv_4;
    int3  _S76 = make_int3 (int(0));
    float3  _S77 = make_float3 ((float)_S76.x, (float)_S76.y, (float)_S76.z);
    *raydir_1 = _S77;
    FixedArray<float, 10>  _S78 = dist_coeffs_4;
    bool _S79 = undistort_point_0(uv_4, &_S78, int(8), &_S75);
    if(!_S79)
    {
        return false;
    }
    float3  _S80;
    if(is_fisheye_2)
    {
        float2  _S81 = _S75;
        float theta_2 = length_0(_S75);
        float _S82;
        if(theta_2 < 0.00100000004749745f)
        {
            _S82 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S82 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S83 = make_float3 ((_S81 * make_float2 (_S82)).x, (_S81 * make_float2 (_S82)).y, (F32_cos((theta_2))));
        _S80 = _S83;
    }
    else
    {
        _S80 = make_float3 (_S75.x, _S75.y, 1.0f);
    }
    *raydir_1 = _S80;
    return true;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float4  normalize_1(float4  x_11)
{
    return x_11 / make_float4 (length_2(x_11));
}

inline __device__ bool generate_ray(float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_5, float3  * raydir_2)
{
    float2  _S84 = uv_5;
    FixedArray<float, 10>  _S85 = dist_coeffs_5;
    bool _S86 = undistort_point_0(uv_5, &_S85, int(8), &_S84);
    if(!_S86)
    {
        int3  _S87 = make_int3 (int(0));
        float3  _S88 = make_float3 ((float)_S87.x, (float)_S87.y, (float)_S87.z);
        *raydir_2 = _S88;
        return false;
    }
    float3  _S89;
    if(is_fisheye_3)
    {
        float2  _S90 = _S84;
        float theta_3 = length_0(_S84);
        float _S91;
        if(theta_3 < 0.00100000004749745f)
        {
            _S91 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S91 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S92 = make_float3 ((_S90 * make_float2 (_S91)).x, (_S90 * make_float2 (_S91)).y, (F32_cos((theta_3))));
        _S89 = _S92;
    }
    else
    {
        _S89 = make_float3 (_S84.x, _S84.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S89);
    return true;
}

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_3)
{
    float _S93 = (*right_2).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_3.x;
    float sum_2 = _S93 + (*right_2).primal_0.rows[int(0)].y * dOut_3.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_3.y;
    float sum_3 = sum_2 + (*right_2).primal_0.rows[int(0)].z * dOut_3.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_3.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_3;
    float _S94 = (*right_2).primal_0.rows[int(1)].x * dOut_3.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_3.x;
    float sum_4 = _S94 + (*right_2).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_3.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(1)].z * dOut_3.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_3.z;
    *&((&left_d_result_1)->y) = sum_5;
    float _S95 = (*right_2).primal_0.rows[int(2)].x * dOut_3.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_3.x;
    float sum_6 = _S95 + (*right_2).primal_0.rows[int(2)].y * dOut_3.y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_3.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(2)].z * dOut_3.z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_3.z;
    *&((&left_d_result_1)->z) = sum_7;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_2(float3  left_3, Matrix<float, 3, 3>  right_3)
{
    float3  result_9;
    int j_0 = int(0);
    for(;;)
    {
        if(j_0 < int(3))
        {
        }
        else
        {
            break;
        }
        int i_5 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(i_5 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_3, i_5) * _slang_vector_get_element(right_3.rows[i_5], j_0);
            i_5 = i_5 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_9, j_0) = sum_8;
        j_0 = j_0 + int(1);
    }
    return result_9;
}

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_0, float3  t_0)
{
    return - mul_2(t_0, R_0);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_1, float3  raydir_3)
{
    return mul_2(raydir_3, R_1);
}

inline __device__ float3  undo_transform_ray_d(Matrix<float, 3, 3>  R_2, float3  raydir_4)
{
    return mul_2(raydir_4, transpose_0(R_2));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S96, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S97, float3  _S98)
{
    _d_mul_0(_S96, _S97, _S98);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S99 = - _s_dOut_0;
    float3  _S100 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S101;
    (&_S101)->primal_0 = (*dpt_0).primal_0;
    (&_S101)->differential_0 = _S100;
    Matrix<float, 3, 3>  _S102 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S103;
    (&_S103)->primal_0 = (*dpR_0).primal_0;
    (&_S103)->differential_0 = _S102;
    s_bwd_prop_mul_0(&_S101, &_S103, _S99);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S101.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S103.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S104, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S105, float3  _S106)
{
    s_bwd_prop_transform_ray_o_0(_S104, _S105, _S106);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S107 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S107;
    float3  _S108 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S108;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S109 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S110;
    (&_S110)->primal_0 = (*dpraydir_0).primal_0;
    (&_S110)->differential_0 = _S109;
    Matrix<float, 3, 3>  _S111 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S112;
    (&_S112)->primal_0 = (*dpR_1).primal_0;
    (&_S112)->differential_0 = _S111;
    s_bwd_prop_mul_0(&_S110, &_S112, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S110.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S112.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S113, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S114, float3  _S115)
{
    s_bwd_prop_transform_ray_d_0(_S113, _S114, _S115);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S116 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S116;
    float3  _S117 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S117;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_2, float3  scale_1)
{
    float4  _S118 = normalize_1(quat_2);
    float x_12 = _S118.y;
    float x2_2 = x_12 * x_12;
    float y2_2 = _S118.z * _S118.z;
    float z2_2 = _S118.w * _S118.w;
    float xy_2 = _S118.y * _S118.z;
    float xz_2 = _S118.y * _S118.w;
    float yz_2 = _S118.z * _S118.w;
    float wx_2 = _S118.x * _S118.y;
    float wy_2 = _S118.x * _S118.z;
    float wz_2 = _S118.x * _S118.w;
    return mul_1(makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)))));
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S119, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S120, Matrix<float, 3, 3>  _S121)
{
    mul_0(_S119, _S120, _S121);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S122, float _S123)
{
    _d_sqrt_0(_S122, _S123);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_2, float _s_dOut_2)
{
    float _S124 = (*dpx_2).primal_0.x;
    float _S125 = (*dpx_2).primal_0.y;
    float _S126 = (*dpx_2).primal_0.z;
    float _S127 = (*dpx_2).primal_0.w;
    DiffPair_float_0 _S128;
    (&_S128)->primal_0 = _S124 * _S124 + _S125 * _S125 + _S126 * _S126 + _S127 * _S127;
    (&_S128)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S128, _s_dOut_2);
    float _S129 = (*dpx_2).primal_0.w * _S128.differential_0;
    float _S130 = _S129 + _S129;
    float _S131 = (*dpx_2).primal_0.z * _S128.differential_0;
    float _S132 = _S131 + _S131;
    float _S133 = (*dpx_2).primal_0.y * _S128.differential_0;
    float _S134 = _S133 + _S133;
    float _S135 = (*dpx_2).primal_0.x * _S128.differential_0;
    float _S136 = _S135 + _S135;
    float4  _S137 = make_float4 (0.0f);
    *&((&_S137)->w) = _S130;
    *&((&_S137)->z) = _S132;
    *&((&_S137)->y) = _S134;
    *&((&_S137)->x) = _S136;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S137;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S138, float _S139)
{
    s_bwd_prop_length_impl_0(_S138, _S139);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_3, float4  _s_dOut_3)
{
    float _S140 = length_2((*dpx_3).primal_0);
    float4  _S141 = (*dpx_3).primal_0 * _s_dOut_3;
    float4  _S142 = make_float4 (1.0f / _S140) * _s_dOut_3;
    float _S143 = - ((_S141.x + _S141.y + _S141.z + _S141.w) / (_S140 * _S140));
    float4  _S144 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S145;
    (&_S145)->primal_0 = (*dpx_3).primal_0;
    (&_S145)->differential_0 = _S144;
    s_bwd_length_impl_0(&_S145, _S143);
    float4  _S146 = _S142 + _S145.differential_0;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S146;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S147, float4  _S148)
{
    s_bwd_prop_normalize_impl_0(_S147, _S148);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_4)
{
    float4  _S149 = normalize_1((*dpquat_0).primal_0);
    float _S150 = _S149.y;
    float x2_3 = _S150 * _S150;
    float y2_3 = _S149.z * _S149.z;
    float z2_3 = _S149.w * _S149.w;
    float xy_3 = _S149.y * _S149.z;
    float xz_3 = _S149.y * _S149.w;
    float yz_3 = _S149.z * _S149.w;
    float wx_3 = _S149.x * _S149.y;
    float wy_3 = _S149.x * _S149.z;
    float wz_3 = _S149.x * _S149.w;
    Matrix<float, 3, 3>  _S151 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))));
    Matrix<float, 3, 3>  _S152 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S153;
    (&_S153)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S153)->differential_0 = _S152;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S154;
    (&_S154)->primal_0 = _S151;
    (&_S154)->differential_0 = _S152;
    s_bwd_prop_mul_1(&_S153, &_S154, _s_dOut_4);
    Matrix<float, 3, 3>  _S155 = transpose_0(transpose_0(_S154.differential_0));
    float _S156 = 2.0f * - _S155.rows[int(2)].z;
    float _S157 = 2.0f * _S155.rows[int(2)].y;
    float _S158 = 2.0f * _S155.rows[int(2)].x;
    float _S159 = 2.0f * _S155.rows[int(1)].z;
    float _S160 = 2.0f * - _S155.rows[int(1)].y;
    float _S161 = 2.0f * _S155.rows[int(1)].x;
    float _S162 = 2.0f * _S155.rows[int(0)].z;
    float _S163 = 2.0f * _S155.rows[int(0)].y;
    float _S164 = 2.0f * - _S155.rows[int(0)].x;
    float _S165 = - _S161 + _S163;
    float _S166 = _S158 + - _S162;
    float _S167 = - _S157 + _S159;
    float _S168 = _S157 + _S159;
    float _S169 = _S158 + _S162;
    float _S170 = _S161 + _S163;
    float _S171 = _S149.w * (_S160 + _S164);
    float _S172 = _S149.z * (_S156 + _S164);
    float _S173 = _S149.y * (_S156 + _S160);
    float _S174 = _S149.x * _S165 + _S149.z * _S168 + _S149.y * _S169 + _S171 + _S171;
    float _S175 = _S149.x * _S166 + _S149.w * _S168 + _S149.y * _S170 + _S172 + _S172;
    float _S176 = _S149.x * _S167 + _S149.w * _S169 + _S149.z * _S170 + _S173 + _S173;
    float _S177 = _S149.w * _S165 + _S149.z * _S166 + _S149.y * _S167;
    float4  _S178 = make_float4 (0.0f);
    float4  _S179 = _S178;
    *&((&_S179)->w) = _S174;
    *&((&_S179)->z) = _S175;
    *&((&_S179)->y) = _S176;
    *&((&_S179)->x) = _S177;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S180;
    (&_S180)->primal_0 = (*dpquat_0).primal_0;
    (&_S180)->differential_0 = _S178;
    s_bwd_normalize_impl_0(&_S180, _S179);
    float3  _S181 = make_float3 (_S153.differential_0.rows[int(0)].x, _S153.differential_0.rows[int(1)].y, _S153.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S181;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S180.differential_0;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S182, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S183, Matrix<float, 3, 3>  _S184)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S182, _S183, _S184);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_3, float3  scale_2, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_0, float3  * v_scale_0)
{
    float4  _S185 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_3;
    (&dp_quat_0)->differential_0 = _S185;
    float3  _S186 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_2;
    (&dp_scale_0)->differential_0 = _S186;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_0 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_4, float3  dOut_4)
{
    float _S187 = (*left_4).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = (*right_4).primal_0.x * dOut_4.x;
    float sum_10 = _S187 + (*left_4).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = (*right_4).primal_0.x * dOut_4.y;
    float sum_11 = sum_10 + (*left_4).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = (*right_4).primal_0.x * dOut_4.z;
    float3  right_d_result_2;
    *&((&right_d_result_2)->x) = sum_11;
    float _S188 = (*left_4).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = (*right_4).primal_0.y * dOut_4.x;
    float sum_12 = _S188 + (*left_4).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = (*right_4).primal_0.y * dOut_4.y;
    float sum_13 = sum_12 + (*left_4).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = (*right_4).primal_0.y * dOut_4.z;
    *&((&right_d_result_2)->y) = sum_13;
    float _S189 = (*left_4).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = (*right_4).primal_0.z * dOut_4.x;
    float sum_14 = _S189 + (*left_4).primal_0.rows[int(1)].z * dOut_4.y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = (*right_4).primal_0.z * dOut_4.y;
    float sum_15 = sum_14 + (*left_4).primal_0.rows[int(2)].z * dOut_4.z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = (*right_4).primal_0.z * dOut_4.z;
    *&((&right_d_result_2)->z) = sum_15;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_2;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_2;
    return;
}

inline __device__ float3  mul_3(Matrix<float, 3, 3>  left_5, float3  right_5)
{
    float3  result_10;
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
        int j_1 = int(0);
        float sum_16 = 0.0f;
        for(;;)
        {
            if(j_1 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_17 = sum_16 + _slang_vector_get_element(left_5.rows[i_6], j_1) * _slang_vector_get_element(right_5, j_1);
            j_1 = j_1 + int(1);
            sum_16 = sum_17;
        }
        *_slang_vector_get_element_ptr(&result_10, i_6) = sum_16;
        i_6 = i_6 + int(1);
    }
    return result_10;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_5)
{
    float _S190 = dOut_5.y;
    float _S191 = dOut_5.z;
    float _S192 = dOut_5.x;
    float _S193 = (*a_0).primal_0.z * _S190 + - (*a_0).primal_0.y * _S191;
    float _S194 = - (*a_0).primal_0.z * _S192 + (*a_0).primal_0.x * _S191;
    float _S195 = (*a_0).primal_0.y * _S192 + - (*a_0).primal_0.x * _S190;
    float3  _S196 = make_float3 (- (*b_0).primal_0.z * _S190 + (*b_0).primal_0.y * _S191, (*b_0).primal_0.z * _S192 + - (*b_0).primal_0.x * _S191, - (*b_0).primal_0.y * _S192 + (*b_0).primal_0.x * _S190);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S196;
    float3  _S197 = make_float3 (_S193, _S194, _S195);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S197;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S198 = left_6.y;
    float _S199 = right_6.z;
    float _S200 = left_6.z;
    float _S201 = right_6.y;
    float _S202 = right_6.x;
    float _S203 = left_6.x;
    return make_float3 (_S198 * _S199 - _S200 * _S201, _S200 * _S202 - _S203 * _S199, _S203 * _S201 - _S198 * _S202);
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_4, float dOut_6)
{
    float _S204 = (F32_exp(((*dpx_4).primal_0))) * dOut_6;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S204;
    return;
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_0, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_3(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_3(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_0 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S205, float3  _S206)
{
    return mul_3(_S205, _S206);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S207, float3  _S208)
{
    return cross_0(_S207, _S208);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S209, float3  _S210)
{
    return dot_0(_S209, _S210);
}

inline __device__ float s_primal_ctx_exp_0(float _S211)
{
    return (F32_exp((_S211)));
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S212, float _S213)
{
    _d_exp_0(_S212, _S213);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S214, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S215, float _S216)
{
    _d_dot_0(_S214, _S215, _S216);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S217, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S218, float3  _S219)
{
    _d_cross_0(_S217, _S218, _S219);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S220, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S221, float3  _S222)
{
    _d_mul_1(_S220, _S221, _S222);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_5)
{
    float3  _S223 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S224 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S223);
    float3  _S225 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S226 = s_primal_ctx_cross_0(_S225, _S224);
    float _S227 = -0.5f * s_primal_ctx_dot_0(_S226, _S226);
    float _S228 = s_primal_ctx_dot_0(_S225, _S225);
    float _S229 = _S227 / _S228;
    float _S230 = _S228 * _S228;
    float _S231 = (*dpopacity_0).primal_0 * _s_dOut_5;
    float _S232 = s_primal_ctx_exp_0(_S229) * _s_dOut_5;
    DiffPair_float_0 _S233;
    (&_S233)->primal_0 = _S229;
    (&_S233)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S233, _S231);
    float _S234 = _S233.differential_0 / _S230;
    float _S235 = _S227 * - _S234;
    float _S236 = _S228 * _S234;
    float3  _S237 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S238;
    (&_S238)->primal_0 = _S225;
    (&_S238)->differential_0 = _S237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S239;
    (&_S239)->primal_0 = _S225;
    (&_S239)->differential_0 = _S237;
    s_bwd_prop_dot_0(&_S238, &_S239, _S235);
    float _S240 = -0.5f * _S236;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S241;
    (&_S241)->primal_0 = _S226;
    (&_S241)->differential_0 = _S237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S226;
    (&_S242)->differential_0 = _S237;
    s_bwd_prop_dot_0(&_S241, &_S242, _S240);
    float3  _S243 = _S242.differential_0 + _S241.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S244;
    (&_S244)->primal_0 = _S225;
    (&_S244)->differential_0 = _S237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S245;
    (&_S245)->primal_0 = _S224;
    (&_S245)->differential_0 = _S237;
    s_bwd_prop_cross_0(&_S244, &_S245, _S243);
    float3  _S246 = _S239.differential_0 + _S238.differential_0 + _S244.differential_0;
    Matrix<float, 3, 3>  _S247 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S248;
    (&_S248)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S248)->differential_0 = _S247;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S249;
    (&_S249)->primal_0 = (*dpray_d_0).primal_0;
    (&_S249)->differential_0 = _S237;
    s_bwd_prop_mul_2(&_S248, &_S249, _S246);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S250;
    (&_S250)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S250)->differential_0 = _S247;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S251;
    (&_S251)->primal_0 = _S223;
    (&_S251)->differential_0 = _S237;
    s_bwd_prop_mul_2(&_S250, &_S251, _S245.differential_0);
    float3  _S252 = - _S251.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S249.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S251.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S232;
    Matrix<float, 3, 3>  _S253 = _S248.differential_0 + _S250.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S253;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S252;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S254, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S255, DiffPair_float_0 * _S256, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S257, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S258, float _S259)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S254, _S255, _S256, _S257, _S258, _S259);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_1, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_0, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S260 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S260;
    Matrix<float, 3, 3>  _S261 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S261;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_1;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S260;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S260;
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
    float3  grd_1 = mul_3(iscl_rot_2, ray_d_2);
    *depth_0 = - dot_0(mul_3(iscl_rot_2, ray_o_2 - mean_2), grd_1) / dot_0(grd_1, grd_1);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S262 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S263 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S262);
    float3  _S264 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S265 = s_primal_ctx_dot_0(_S264, _S264);
    float _S266 = dpdepth_0 / (_S265 * _S265);
    float _S267 = - s_primal_ctx_dot_0(_S263, _S264) * - _S266;
    float _S268 = _S265 * _S266;
    float3  _S269 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S270;
    (&_S270)->primal_0 = _S264;
    (&_S270)->differential_0 = _S269;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S271;
    (&_S271)->primal_0 = _S264;
    (&_S271)->differential_0 = _S269;
    s_bwd_prop_dot_0(&_S270, &_S271, _S267);
    float _S272 = - _S268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S273;
    (&_S273)->primal_0 = _S263;
    (&_S273)->differential_0 = _S269;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S274;
    (&_S274)->primal_0 = _S264;
    (&_S274)->differential_0 = _S269;
    s_bwd_prop_dot_0(&_S273, &_S274, _S272);
    float3  _S275 = _S271.differential_0 + _S270.differential_0 + _S274.differential_0;
    Matrix<float, 3, 3>  _S276 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S277;
    (&_S277)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S277)->differential_0 = _S276;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S278;
    (&_S278)->primal_0 = (*dpray_d_1).primal_0;
    (&_S278)->differential_0 = _S269;
    s_bwd_prop_mul_2(&_S277, &_S278, _S275);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S279;
    (&_S279)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S279)->differential_0 = _S276;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S280;
    (&_S280)->primal_0 = _S262;
    (&_S280)->differential_0 = _S269;
    s_bwd_prop_mul_2(&_S279, &_S280, _S273.differential_0);
    float3  _S281 = - _S280.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S278.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S280.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S282 = _S277.differential_0 + _S279.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S282;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S281;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S283, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S284, DiffPair_float_0 * _S285, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S286, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S287, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S288, float3  _S289, float _S290)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S283, _S284, _S285, _S286, _S287, _S288, _S289, _S290);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_3, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_1, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S291 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S291;
    Matrix<float, 3, 3>  _S292 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S292;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_3;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S291;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S291;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S291;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_mean_1 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_1 = dp_opacity_1.differential_0;
    *v_rgb_0 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_4, float4  quat_4, float3  scale_3, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S293 = scale_3.x;
    float sx_0 = (F32_exp((_S293)));
    float _S294 = scale_3.y;
    float sy_0 = (F32_exp((_S294)));
    float sz_0 = scale_3.z - 0.5f * (_S293 + _S294);
    float x_13 = quat_4.y;
    float x2_4 = x_13 * x_13;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_4 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  _S295 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    *vert0_0 = mul_3(_S295, make_float3 (sx_0, 0.0f, 0.0f)) + mean_4;
    *vert1_0 = mul_3(_S295, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_4;
    *vert2_0 = mul_3(_S295, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_4;
    return;
}

