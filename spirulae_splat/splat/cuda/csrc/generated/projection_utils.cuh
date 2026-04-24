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

inline __device__ bool persp_proj_nav(float3  p_view_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, float2  * uv_1)
{
    bool _S15;
    for(;;)
    {
        float _S16 = p_view_0.z;
        *uv_1 = float2 {p_view_0.x, p_view_0.y} / make_float2 (_S16);
        if(_S16 < 0.0f)
        {
            _S15 = true;
        }
        else
        {
            float u_1 = (*uv_1).x;
            float v_1 = (*uv_1).y;
            float _S17 = 0.0f * v_1;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float s_diff_r2_2 = u_1 + u_1 + (_S17 + _S17);
            float _S18 = dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)];
            float _S19 = dist_coeffs_1[int(1)] + r2_1 * _S18;
            float _S20 = dist_coeffs_1[int(0)] + r2_1 * _S19;
            float radial_1 = 1.0f + r2_1 * _S20;
            float _S21 = 2.0f * dist_coeffs_1[int(4)];
            float _S22 = _S21 * u_1;
            float _S23 = 2.0f * u_1;
            float _S24 = 2.0f * dist_coeffs_1[int(5)];
            float _S25 = _S24 * u_1;
            float _S26 = 2.0f * v_1;
            float2  _S27 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S20 + (s_diff_r2_2 * _S19 + (s_diff_r2_2 * _S18 + s_diff_r2_2 * dist_coeffs_1[int(3)] * r2_1) * r2_1) * r2_1) * *uv_1 + make_float2 (_S21 * v_1 + 0.0f * _S22 + (s_diff_r2_2 + (_S23 + _S23)) * dist_coeffs_1[int(5)] + s_diff_r2_2 * dist_coeffs_1[int(6)], _S24 * v_1 + 0.0f * _S25 + (s_diff_r2_2 + (_S17 + 0.0f * _S26)) * dist_coeffs_1[int(4)] + s_diff_r2_2 * dist_coeffs_1[int(7)]);
            float _S28 = 0.0f * u_1;
            float s_diff_r2_3 = _S28 + _S28 + (v_1 + v_1);
            float2  _S29 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S20 + (s_diff_r2_3 * _S19 + (s_diff_r2_3 * _S18 + s_diff_r2_3 * dist_coeffs_1[int(3)] * r2_1) * r2_1) * r2_1) * *uv_1 + make_float2 (0.0f * _S21 * v_1 + _S22 + (s_diff_r2_3 + (_S28 + 0.0f * _S23)) * dist_coeffs_1[int(5)] + s_diff_r2_3 * dist_coeffs_1[int(6)], 0.0f * _S24 * v_1 + _S25 + (s_diff_r2_3 + (_S26 + _S26)) * dist_coeffs_1[int(4)] + s_diff_r2_3 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S30 = transpose_1(makeMatrix<float, 2, 2> (_S27 + make_float2 (_S27.x * dist_coeffs_1[int(8)] + _S27.y * dist_coeffs_1[int(9)], 0.0f), _S29 + make_float2 (_S29.x * dist_coeffs_1[int(8)] + _S29.y * dist_coeffs_1[int(9)], 0.0f)));
            _S15 = !((F32_min((determinant_0(_S30)), ((F32_min((_S30.rows[int(0)].x), (_S30.rows[int(1)].y)))))) > 0.0f);
        }
        if(_S15)
        {
            break;
        }
        float u_2 = (*uv_1).x;
        float v_2 = (*uv_1).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S31 = *uv_1 * make_float2 (1.0f + r2_2 * (dist_coeffs_1[int(0)] + r2_2 * (dist_coeffs_1[int(1)] + r2_2 * (dist_coeffs_1[int(2)] + r2_2 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_2 * v_2 + dist_coeffs_1[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + dist_coeffs_1[int(6)] * r2_2, 2.0f * dist_coeffs_1[int(5)] * u_2 * v_2 + dist_coeffs_1[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + dist_coeffs_1[int(7)] * r2_2);
        float2  _S32 = _S31 + make_float2 (dist_coeffs_1[int(8)] * _S31.x + dist_coeffs_1[int(9)] * _S31.y, 0.0f);
        *uv_1 = make_float2 (intrins_0.x * _S32.x + intrins_0.z, intrins_0.y * _S32.y + intrins_0.w);
        break;
    }
    return !_S15;
}

inline __device__ Matrix<float, 2, 3>  persp_proj_jac(float3  p_view_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2)
{
    float2  _S33 = float2 {p_view_1.x, p_view_1.y};
    float _S34 = p_view_1.z;
    float2  _S35 = _S33 / make_float2 (_S34);
    float2  _S36 = _S33 * make_float2 (0.0f);
    float _S37 = _S34 * _S34;
    float2  _S38 = (make_float2 (1.0f, 0.0f) * make_float2 (_S34) - _S36) / make_float2 (_S37);
    float u_3 = _S35.x;
    float s_diff_u_0 = _S38.x;
    float v_3 = _S35.y;
    float s_diff_v_0 = _S38.y;
    float _S39 = s_diff_u_0 * u_3;
    float _S40 = s_diff_v_0 * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = _S39 + _S39 + (_S40 + _S40);
    float _S41 = dist_coeffs_2[int(2)] + r2_3 * dist_coeffs_2[int(3)];
    float _S42 = dist_coeffs_2[int(1)] + r2_3 * _S41;
    float _S43 = dist_coeffs_2[int(0)] + r2_3 * _S42;
    float _S44 = 2.0f * dist_coeffs_2[int(4)];
    float _S45 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S46 = _S38 * make_float2 (1.0f + r2_3 * _S43) + make_float2 (s_diff_r2_4 * _S43 + (s_diff_r2_4 * _S42 + (s_diff_r2_4 * _S41 + s_diff_r2_4 * dist_coeffs_2[int(3)] * r2_3) * r2_3) * r2_3) * _S35 + make_float2 (s_diff_u_0 * _S44 * v_3 + s_diff_v_0 * (_S44 * u_3) + (s_diff_r2_4 + (s_diff_u_0 * 2.0f * u_3 + s_diff_u_0 * (2.0f * u_3))) * dist_coeffs_2[int(5)] + s_diff_r2_4 * dist_coeffs_2[int(6)], s_diff_u_0 * _S45 * v_3 + s_diff_v_0 * (_S45 * u_3) + (s_diff_r2_4 + (s_diff_v_0 * 2.0f * v_3 + s_diff_v_0 * (2.0f * v_3))) * dist_coeffs_2[int(4)] + s_diff_r2_4 * dist_coeffs_2[int(7)]);
    float2  _S47 = _S46 + make_float2 (_S46.x * dist_coeffs_2[int(8)] + _S46.y * dist_coeffs_2[int(9)], 0.0f);
    float fx_0 = intrins_1.x;
    float fy_0 = intrins_1.y;
    float _S48 = _S47.y * fy_0;
    Matrix<float, 2, 3>  J_0;
    *&(((&J_0)->rows + (int(0)))->x) = _S47.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->x) = _S48;
    float2  _S49 = _S33 / make_float2 (_S34);
    float2  _S50 = (make_float2 (0.0f, 1.0f) * make_float2 (_S34) - _S36) / make_float2 (_S37);
    float u_4 = _S49.x;
    float s_diff_u_1 = _S50.x;
    float v_4 = _S49.y;
    float s_diff_v_1 = _S50.y;
    float _S51 = s_diff_u_1 * u_4;
    float _S52 = s_diff_v_1 * v_4;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float s_diff_r2_5 = _S51 + _S51 + (_S52 + _S52);
    float _S53 = dist_coeffs_2[int(2)] + r2_4 * dist_coeffs_2[int(3)];
    float _S54 = dist_coeffs_2[int(1)] + r2_4 * _S53;
    float _S55 = dist_coeffs_2[int(0)] + r2_4 * _S54;
    float2  _S56 = _S50 * make_float2 (1.0f + r2_4 * _S55) + make_float2 (s_diff_r2_5 * _S55 + (s_diff_r2_5 * _S54 + (s_diff_r2_5 * _S53 + s_diff_r2_5 * dist_coeffs_2[int(3)] * r2_4) * r2_4) * r2_4) * _S49 + make_float2 (s_diff_u_1 * _S44 * v_4 + s_diff_v_1 * (_S44 * u_4) + (s_diff_r2_5 + (s_diff_u_1 * 2.0f * u_4 + s_diff_u_1 * (2.0f * u_4))) * dist_coeffs_2[int(5)] + s_diff_r2_5 * dist_coeffs_2[int(6)], s_diff_u_1 * _S45 * v_4 + s_diff_v_1 * (_S45 * u_4) + (s_diff_r2_5 + (s_diff_v_1 * 2.0f * v_4 + s_diff_v_1 * (2.0f * v_4))) * dist_coeffs_2[int(4)] + s_diff_r2_5 * dist_coeffs_2[int(7)]);
    float2  _S57 = _S56 + make_float2 (_S56.x * dist_coeffs_2[int(8)] + _S56.y * dist_coeffs_2[int(9)], 0.0f);
    float _S58 = _S57.y * fy_0;
    *&(((&J_0)->rows + (int(0)))->y) = _S57.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->y) = _S58;
    float2  _S59 = _S33 / make_float2 (_S34);
    float2  _S60 = (make_float2 (0.0f, 0.0f) * make_float2 (_S34) - _S33) / make_float2 (_S37);
    float u_5 = _S59.x;
    float s_diff_u_2 = _S60.x;
    float v_5 = _S59.y;
    float s_diff_v_2 = _S60.y;
    float _S61 = s_diff_u_2 * u_5;
    float _S62 = s_diff_v_2 * v_5;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float s_diff_r2_6 = _S61 + _S61 + (_S62 + _S62);
    float _S63 = dist_coeffs_2[int(2)] + r2_5 * dist_coeffs_2[int(3)];
    float _S64 = dist_coeffs_2[int(1)] + r2_5 * _S63;
    float _S65 = dist_coeffs_2[int(0)] + r2_5 * _S64;
    float2  _S66 = _S60 * make_float2 (1.0f + r2_5 * _S65) + make_float2 (s_diff_r2_6 * _S65 + (s_diff_r2_6 * _S64 + (s_diff_r2_6 * _S63 + s_diff_r2_6 * dist_coeffs_2[int(3)] * r2_5) * r2_5) * r2_5) * _S59 + make_float2 (s_diff_u_2 * _S44 * v_5 + s_diff_v_2 * (_S44 * u_5) + (s_diff_r2_6 + (s_diff_u_2 * 2.0f * u_5 + s_diff_u_2 * (2.0f * u_5))) * dist_coeffs_2[int(5)] + s_diff_r2_6 * dist_coeffs_2[int(6)], s_diff_u_2 * _S45 * v_5 + s_diff_v_2 * (_S45 * u_5) + (s_diff_r2_6 + (s_diff_v_2 * 2.0f * v_5 + s_diff_v_2 * (2.0f * v_5))) * dist_coeffs_2[int(4)] + s_diff_r2_6 * dist_coeffs_2[int(7)]);
    float2  _S67 = _S66 + make_float2 (_S66.x * dist_coeffs_2[int(8)] + _S66.y * dist_coeffs_2[int(9)], 0.0f);
    float _S68 = _S67.y * fy_0;
    *&(((&J_0)->rows + (int(0)))->z) = _S67.x * fx_0;
    *&(((&J_0)->rows + (int(1)))->z) = _S68;
    return J_0;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ DiffPair_float_0 _d_sqrt_0(DiffPair_float_0 * dpx_0)
{
    DiffPair_float_0 _S69 = { (F32_sqrt((dpx_0->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_0->primal_0)))))) * dpx_0->differential_0 };
    return _S69;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, float dOut_1)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_0).primal_0.x * dOut_1;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_1).primal_0.x * dOut_1;
    *&((&x_d_result_0)->y) = (*dpy_0).primal_0.y * dOut_1;
    *&((&y_d_result_0)->y) = (*dpx_1).primal_0.y * dOut_1;
    *&((&x_d_result_0)->z) = (*dpy_0).primal_0.z * dOut_1;
    *&((&y_d_result_0)->z) = (*dpx_1).primal_0.z * dOut_1;
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

inline __device__ DiffPair_float_0 _d_atan2_0(DiffPair_float_0 * dpy_1, DiffPair_float_0 * dpx_2)
{
    float _S70 = dpx_2->primal_0 * dpx_2->primal_0 + dpy_1->primal_0 * dpy_1->primal_0;
    DiffPair_float_0 _S71 = { (F32_atan2((dpy_1->primal_0), (dpx_2->primal_0))), - dpy_1->primal_0 / _S70 * dpx_2->differential_0 + dpx_2->primal_0 / _S70 * dpy_1->differential_0 };
    return _S71;
}

inline __device__ bool fisheye_proj_nav(float3  p_view_2, float4  intrins_2, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_2)
{
    bool _S72;
    for(;;)
    {
        float2  _S73 = float2 {p_view_2.x, p_view_2.y};
        float r_3 = length_0(_S73);
        float _S74 = p_view_2.z;
        float theta_0 = (F32_atan2((r_3), (_S74)));
        float k_0;
        if(theta_0 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S74;
        }
        else
        {
            k_0 = theta_0 / r_3;
        }
        float2  _S75 = _S73 * make_float2 (k_0);
        *uv_2 = _S75;
        float u_6 = _S75.x;
        float v_6 = _S75.y;
        float _S76 = 0.0f * v_6;
        float r2_6 = u_6 * u_6 + v_6 * v_6;
        float s_diff_r2_7 = u_6 + u_6 + (_S76 + _S76);
        float _S77 = dist_coeffs_3[int(2)] + r2_6 * dist_coeffs_3[int(3)];
        float _S78 = dist_coeffs_3[int(1)] + r2_6 * _S77;
        float _S79 = dist_coeffs_3[int(0)] + r2_6 * _S78;
        float radial_2 = 1.0f + r2_6 * _S79;
        float _S80 = 2.0f * dist_coeffs_3[int(4)];
        float _S81 = _S80 * u_6;
        float _S82 = 2.0f * u_6;
        float _S83 = 2.0f * dist_coeffs_3[int(5)];
        float _S84 = _S83 * u_6;
        float _S85 = 2.0f * v_6;
        float2  _S86 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_7 * _S79 + (s_diff_r2_7 * _S78 + (s_diff_r2_7 * _S77 + s_diff_r2_7 * dist_coeffs_3[int(3)] * r2_6) * r2_6) * r2_6) * _S75 + make_float2 (_S80 * v_6 + 0.0f * _S81 + (s_diff_r2_7 + (_S82 + _S82)) * dist_coeffs_3[int(5)] + s_diff_r2_7 * dist_coeffs_3[int(6)], _S83 * v_6 + 0.0f * _S84 + (s_diff_r2_7 + (_S76 + 0.0f * _S85)) * dist_coeffs_3[int(4)] + s_diff_r2_7 * dist_coeffs_3[int(7)]);
        float _S87 = 0.0f * u_6;
        float s_diff_r2_8 = _S87 + _S87 + (v_6 + v_6);
        float2  _S88 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_8 * _S79 + (s_diff_r2_8 * _S78 + (s_diff_r2_8 * _S77 + s_diff_r2_8 * dist_coeffs_3[int(3)] * r2_6) * r2_6) * r2_6) * _S75 + make_float2 (0.0f * _S80 * v_6 + _S81 + (s_diff_r2_8 + (_S87 + 0.0f * _S82)) * dist_coeffs_3[int(5)] + s_diff_r2_8 * dist_coeffs_3[int(6)], 0.0f * _S83 * v_6 + _S84 + (s_diff_r2_8 + (_S85 + _S85)) * dist_coeffs_3[int(4)] + s_diff_r2_8 * dist_coeffs_3[int(7)]);
        Matrix<float, 2, 2>  _S89 = transpose_1(makeMatrix<float, 2, 2> (_S86 + make_float2 (_S86.x * dist_coeffs_3[int(8)] + _S86.y * dist_coeffs_3[int(9)], 0.0f), _S88 + make_float2 (_S88.x * dist_coeffs_3[int(8)] + _S88.y * dist_coeffs_3[int(9)], 0.0f)));
        bool _S90 = !((F32_min((determinant_0(_S89)), ((F32_min((_S89.rows[int(0)].x), (_S89.rows[int(1)].y)))))) > 0.0f);
        _S72 = _S90;
        if(_S90)
        {
            break;
        }
        float u_7 = (*uv_2).x;
        float v_7 = (*uv_2).y;
        float r2_7 = u_7 * u_7 + v_7 * v_7;
        float2  _S91 = *uv_2 * make_float2 (1.0f + r2_7 * (dist_coeffs_3[int(0)] + r2_7 * (dist_coeffs_3[int(1)] + r2_7 * (dist_coeffs_3[int(2)] + r2_7 * dist_coeffs_3[int(3)])))) + make_float2 (_S80 * u_7 * v_7 + dist_coeffs_3[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + dist_coeffs_3[int(6)] * r2_7, _S83 * u_7 * v_7 + dist_coeffs_3[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + dist_coeffs_3[int(7)] * r2_7);
        float2  _S92 = _S91 + make_float2 (dist_coeffs_3[int(8)] * _S91.x + dist_coeffs_3[int(9)] * _S91.y, 0.0f);
        *uv_2 = make_float2 (intrins_2.x * _S92.x + intrins_2.z, intrins_2.y * _S92.y + intrins_2.w);
        break;
    }
    return !_S72;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_3)
{
    float _S93 = *&((&dpx_3->differential_0)->x) * *&((&dpx_3->primal_0)->x);
    float _S94 = *&((&dpx_3->differential_0)->y) * *&((&dpx_3->primal_0)->y);
    float s_diff_len_0 = _S93 + _S93 + (_S94 + _S94);
    DiffPair_float_0 _S95;
    (&_S95)->primal_0 = *&((&dpx_3->primal_0)->x) * *&((&dpx_3->primal_0)->x) + *&((&dpx_3->primal_0)->y) * *&((&dpx_3->primal_0)->y);
    (&_S95)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S96 = _d_sqrt_0(&_S95);
    DiffPair_float_0 _S97 = { _S96.primal_0, _S96.differential_0 };
    return _S97;
}

inline __device__ Matrix<float, 2, 3>  fisheye_proj_jac(float3  p_view_3, float4  intrins_3, FixedArray<float, 10>  dist_coeffs_4)
{
    Matrix<float, 2, 3>  J_1;
    float2  _S98 = float2 {p_view_3.x, p_view_3.y};
    float2  _S99 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S100;
    (&_S100)->primal_0 = _S98;
    (&_S100)->differential_0 = _S99;
    DiffPair_float_0 _S101 = s_fwd_length_impl_0(&_S100);
    float _S102 = p_view_3.z;
    DiffPair_float_0 _S103;
    (&_S103)->primal_0 = _S101.primal_0;
    (&_S103)->differential_0 = _S101.differential_0;
    DiffPair_float_0 _S104;
    (&_S104)->primal_0 = _S102;
    (&_S104)->differential_0 = 0.0f;
    DiffPair_float_0 _S105 = _d_atan2_0(&_S103, &_S104);
    float k_1;
    float s_diff_k_0;
    if((_S105.primal_0) < 0.00100000004749745f)
    {
        float _S106 = _S105.differential_0 * _S105.primal_0;
        float _S107 = 1.0f - _S105.primal_0 * _S105.primal_0 / 3.0f;
        float _S108 = ((0.0f - (_S106 + _S106) * 0.3333333432674408f) * _S102 - _S107 * 0.0f) / (_S102 * _S102);
        k_1 = _S107 / _S102;
        s_diff_k_0 = _S108;
    }
    else
    {
        float _S109 = (_S105.differential_0 * _S101.primal_0 - _S105.primal_0 * _S101.differential_0) / (_S101.primal_0 * _S101.primal_0);
        k_1 = _S105.primal_0 / _S101.primal_0;
        s_diff_k_0 = _S109;
    }
    float2  _S110 = _S98 * make_float2 (k_1);
    float2  _S111 = _S99 * make_float2 (k_1) + make_float2 (s_diff_k_0) * _S98;
    float u_8 = _S110.x;
    float s_diff_u_3 = _S111.x;
    float v_8 = _S110.y;
    float s_diff_v_3 = _S111.y;
    float _S112 = s_diff_u_3 * u_8;
    float _S113 = s_diff_v_3 * v_8;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float s_diff_r2_9 = _S112 + _S112 + (_S113 + _S113);
    float _S114 = dist_coeffs_4[int(2)] + r2_8 * dist_coeffs_4[int(3)];
    float _S115 = dist_coeffs_4[int(1)] + r2_8 * _S114;
    float _S116 = dist_coeffs_4[int(0)] + r2_8 * _S115;
    float _S117 = 2.0f * dist_coeffs_4[int(4)];
    float _S118 = 2.0f * dist_coeffs_4[int(5)];
    float2  _S119 = _S111 * make_float2 (1.0f + r2_8 * _S116) + make_float2 (s_diff_r2_9 * _S116 + (s_diff_r2_9 * _S115 + (s_diff_r2_9 * _S114 + s_diff_r2_9 * dist_coeffs_4[int(3)] * r2_8) * r2_8) * r2_8) * _S110 + make_float2 (s_diff_u_3 * _S117 * v_8 + s_diff_v_3 * (_S117 * u_8) + (s_diff_r2_9 + (s_diff_u_3 * 2.0f * u_8 + s_diff_u_3 * (2.0f * u_8))) * dist_coeffs_4[int(5)] + s_diff_r2_9 * dist_coeffs_4[int(6)], s_diff_u_3 * _S118 * v_8 + s_diff_v_3 * (_S118 * u_8) + (s_diff_r2_9 + (s_diff_v_3 * 2.0f * v_8 + s_diff_v_3 * (2.0f * v_8))) * dist_coeffs_4[int(4)] + s_diff_r2_9 * dist_coeffs_4[int(7)]);
    float2  _S120 = _S119 + make_float2 (_S119.x * dist_coeffs_4[int(8)] + _S119.y * dist_coeffs_4[int(9)], 0.0f);
    float fx_1 = intrins_3.x;
    float fy_1 = intrins_3.y;
    float _S121 = _S120.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->x) = _S120.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->x) = _S121;
    float2  _S122 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S123;
    (&_S123)->primal_0 = _S98;
    (&_S123)->differential_0 = _S122;
    DiffPair_float_0 _S124 = s_fwd_length_impl_0(&_S123);
    DiffPair_float_0 _S125;
    (&_S125)->primal_0 = _S124.primal_0;
    (&_S125)->differential_0 = _S124.differential_0;
    DiffPair_float_0 _S126;
    (&_S126)->primal_0 = _S102;
    (&_S126)->differential_0 = 0.0f;
    DiffPair_float_0 _S127 = _d_atan2_0(&_S125, &_S126);
    if((_S127.primal_0) < 0.00100000004749745f)
    {
        float _S128 = _S127.differential_0 * _S127.primal_0;
        float _S129 = 1.0f - _S127.primal_0 * _S127.primal_0 / 3.0f;
        float _S130 = ((0.0f - (_S128 + _S128) * 0.3333333432674408f) * _S102 - _S129 * 0.0f) / (_S102 * _S102);
        k_1 = _S129 / _S102;
        s_diff_k_0 = _S130;
    }
    else
    {
        float _S131 = (_S127.differential_0 * _S124.primal_0 - _S127.primal_0 * _S124.differential_0) / (_S124.primal_0 * _S124.primal_0);
        k_1 = _S127.primal_0 / _S124.primal_0;
        s_diff_k_0 = _S131;
    }
    float2  _S132 = _S98 * make_float2 (k_1);
    float2  _S133 = _S122 * make_float2 (k_1) + make_float2 (s_diff_k_0) * _S98;
    float u_9 = _S132.x;
    float s_diff_u_4 = _S133.x;
    float v_9 = _S132.y;
    float s_diff_v_4 = _S133.y;
    float _S134 = s_diff_u_4 * u_9;
    float _S135 = s_diff_v_4 * v_9;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float s_diff_r2_10 = _S134 + _S134 + (_S135 + _S135);
    float _S136 = dist_coeffs_4[int(2)] + r2_9 * dist_coeffs_4[int(3)];
    float _S137 = dist_coeffs_4[int(1)] + r2_9 * _S136;
    float _S138 = dist_coeffs_4[int(0)] + r2_9 * _S137;
    float2  _S139 = _S133 * make_float2 (1.0f + r2_9 * _S138) + make_float2 (s_diff_r2_10 * _S138 + (s_diff_r2_10 * _S137 + (s_diff_r2_10 * _S136 + s_diff_r2_10 * dist_coeffs_4[int(3)] * r2_9) * r2_9) * r2_9) * _S132 + make_float2 (s_diff_u_4 * _S117 * v_9 + s_diff_v_4 * (_S117 * u_9) + (s_diff_r2_10 + (s_diff_u_4 * 2.0f * u_9 + s_diff_u_4 * (2.0f * u_9))) * dist_coeffs_4[int(5)] + s_diff_r2_10 * dist_coeffs_4[int(6)], s_diff_u_4 * _S118 * v_9 + s_diff_v_4 * (_S118 * u_9) + (s_diff_r2_10 + (s_diff_v_4 * 2.0f * v_9 + s_diff_v_4 * (2.0f * v_9))) * dist_coeffs_4[int(4)] + s_diff_r2_10 * dist_coeffs_4[int(7)]);
    float2  _S140 = _S139 + make_float2 (_S139.x * dist_coeffs_4[int(8)] + _S139.y * dist_coeffs_4[int(9)], 0.0f);
    float _S141 = _S140.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->y) = _S140.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->y) = _S141;
    float2  _S142 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S143;
    (&_S143)->primal_0 = _S98;
    (&_S143)->differential_0 = _S142;
    DiffPair_float_0 _S144 = s_fwd_length_impl_0(&_S143);
    DiffPair_float_0 _S145;
    (&_S145)->primal_0 = _S144.primal_0;
    (&_S145)->differential_0 = _S144.differential_0;
    DiffPair_float_0 _S146;
    (&_S146)->primal_0 = _S102;
    (&_S146)->differential_0 = 1.0f;
    DiffPair_float_0 _S147 = _d_atan2_0(&_S145, &_S146);
    if((_S147.primal_0) < 0.00100000004749745f)
    {
        float _S148 = _S147.differential_0 * _S147.primal_0;
        float _S149 = 1.0f - _S147.primal_0 * _S147.primal_0 / 3.0f;
        float _S150 = ((0.0f - (_S148 + _S148) * 0.3333333432674408f) * _S102 - _S149) / (_S102 * _S102);
        k_1 = _S149 / _S102;
        s_diff_k_0 = _S150;
    }
    else
    {
        float _S151 = (_S147.differential_0 * _S144.primal_0 - _S147.primal_0 * _S144.differential_0) / (_S144.primal_0 * _S144.primal_0);
        k_1 = _S147.primal_0 / _S144.primal_0;
        s_diff_k_0 = _S151;
    }
    float2  _S152 = _S98 * make_float2 (k_1);
    float2  _S153 = _S142 * make_float2 (k_1) + make_float2 (s_diff_k_0) * _S98;
    float u_10 = _S152.x;
    float s_diff_u_5 = _S153.x;
    float v_10 = _S152.y;
    float s_diff_v_5 = _S153.y;
    float _S154 = s_diff_u_5 * u_10;
    float _S155 = s_diff_v_5 * v_10;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float s_diff_r2_11 = _S154 + _S154 + (_S155 + _S155);
    float _S156 = dist_coeffs_4[int(2)] + r2_10 * dist_coeffs_4[int(3)];
    float _S157 = dist_coeffs_4[int(1)] + r2_10 * _S156;
    float _S158 = dist_coeffs_4[int(0)] + r2_10 * _S157;
    float2  _S159 = _S153 * make_float2 (1.0f + r2_10 * _S158) + make_float2 (s_diff_r2_11 * _S158 + (s_diff_r2_11 * _S157 + (s_diff_r2_11 * _S156 + s_diff_r2_11 * dist_coeffs_4[int(3)] * r2_10) * r2_10) * r2_10) * _S152 + make_float2 (s_diff_u_5 * _S117 * v_10 + s_diff_v_5 * (_S117 * u_10) + (s_diff_r2_11 + (s_diff_u_5 * 2.0f * u_10 + s_diff_u_5 * (2.0f * u_10))) * dist_coeffs_4[int(5)] + s_diff_r2_11 * dist_coeffs_4[int(6)], s_diff_u_5 * _S118 * v_10 + s_diff_v_5 * (_S118 * u_10) + (s_diff_r2_11 + (s_diff_v_5 * 2.0f * v_10 + s_diff_v_5 * (2.0f * v_10))) * dist_coeffs_4[int(4)] + s_diff_r2_11 * dist_coeffs_4[int(7)]);
    float2  _S160 = _S159 + make_float2 (_S159.x * dist_coeffs_4[int(8)] + _S159.y * dist_coeffs_4[int(9)], 0.0f);
    float _S161 = _S160.y * fy_1;
    *&(((&J_1)->rows + (int(0)))->z) = _S160.x * fx_1;
    *&(((&J_1)->rows + (int(1)))->z) = _S161;
    return J_1;
}

inline __device__ float2  distort_point(float2  uv_3, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_5)
{
    float2  _S162;
    if(is_fisheye_0)
    {
        float r_4 = length_0(uv_3);
        float theta_1 = (F32_atan((r_4)));
        float _S163;
        if(r_4 < 0.00100000004749745f)
        {
            _S163 = 1.0f - r_4 * r_4 / 3.0f;
        }
        else
        {
            _S163 = theta_1 / r_4;
        }
        _S162 = uv_3 * make_float2 (_S163);
    }
    else
    {
        _S162 = uv_3;
    }
    float u_11 = _S162.x;
    float v_11 = _S162.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float2  _S164 = _S162 * make_float2 (1.0f + r2_11 * (dist_coeffs_5[int(0)] + r2_11 * (dist_coeffs_5[int(1)] + r2_11 * (dist_coeffs_5[int(2)] + r2_11 * dist_coeffs_5[int(3)])))) + make_float2 (2.0f * dist_coeffs_5[int(4)] * u_11 * v_11 + dist_coeffs_5[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_5[int(6)] * r2_11, 2.0f * dist_coeffs_5[int(5)] * u_11 * v_11 + dist_coeffs_5[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_5[int(7)] * r2_11);
    return _S164 + make_float2 (dist_coeffs_5[int(8)] * _S164.x + dist_coeffs_5[int(9)] * _S164.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_4, FixedArray<float, 10>  * dist_coeffs_6, int maxiter_0, float2  * uv_undist_0)
{
    int i_4 = int(0);
    float2  q_0 = uv_4;
    for(;;)
    {
        if(i_4 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S165 = (*dist_coeffs_6)[int(3)];
        float _S166 = (*dist_coeffs_6)[int(4)];
        float _S167 = (*dist_coeffs_6)[int(5)];
        float _S168 = (*dist_coeffs_6)[int(6)];
        float _S169 = (*dist_coeffs_6)[int(7)];
        float _S170 = (*dist_coeffs_6)[int(8)];
        float _S171 = (*dist_coeffs_6)[int(9)];
        float u_12 = q_0.x;
        float v_12 = q_0.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float _S172 = (*dist_coeffs_6)[int(2)] + r2_12 * (*dist_coeffs_6)[int(3)];
        float _S173 = (*dist_coeffs_6)[int(1)] + r2_12 * _S172;
        float _S174 = (*dist_coeffs_6)[int(0)] + r2_12 * _S173;
        float radial_3 = 1.0f + r2_12 * _S174;
        float _S175 = 2.0f * (*dist_coeffs_6)[int(4)];
        float _S176 = _S175 * u_12;
        float _S177 = 2.0f * u_12;
        float _S178 = 2.0f * (*dist_coeffs_6)[int(5)];
        float _S179 = _S178 * u_12;
        float _S180 = 2.0f * v_12;
        float2  _S181 = q_0 * make_float2 (radial_3) + make_float2 (_S176 * v_12 + (*dist_coeffs_6)[int(5)] * (r2_12 + _S177 * u_12) + (*dist_coeffs_6)[int(6)] * r2_12, _S179 * v_12 + (*dist_coeffs_6)[int(4)] * (r2_12 + _S180 * v_12) + (*dist_coeffs_6)[int(7)] * r2_12);
        float2  r_5 = _S181 + make_float2 ((*dist_coeffs_6)[int(8)] * _S181.x + (*dist_coeffs_6)[int(9)] * _S181.y, 0.0f) - uv_4;
        float _S182 = 0.0f * v_12;
        float s_diff_r2_12 = u_12 + u_12 + (_S182 + _S182);
        float2  _S183 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_12 * _S174 + (s_diff_r2_12 * _S173 + (s_diff_r2_12 * _S172 + s_diff_r2_12 * _S165 * r2_12) * r2_12) * r2_12) * q_0 + make_float2 (_S175 * v_12 + 0.0f * _S176 + (s_diff_r2_12 + (_S177 + _S177)) * _S167 + s_diff_r2_12 * _S168, _S178 * v_12 + 0.0f * _S179 + (s_diff_r2_12 + (_S182 + 0.0f * _S180)) * _S166 + s_diff_r2_12 * _S169);
        float _S184 = 0.0f * u_12;
        float s_diff_r2_13 = _S184 + _S184 + (v_12 + v_12);
        float2  _S185 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_13 * _S174 + (s_diff_r2_13 * _S173 + (s_diff_r2_13 * _S172 + s_diff_r2_13 * _S165 * r2_12) * r2_12) * r2_12) * q_0 + make_float2 (0.0f * _S175 * v_12 + _S176 + (s_diff_r2_13 + (_S184 + 0.0f * _S177)) * _S167 + s_diff_r2_13 * _S168, 0.0f * _S178 * v_12 + _S179 + (s_diff_r2_13 + (_S180 + _S180)) * _S166 + s_diff_r2_13 * _S169);
        Matrix<float, 2, 2>  _S186 = transpose_1(makeMatrix<float, 2, 2> (_S183 + make_float2 (_S183.x * _S170 + _S183.y * _S171, 0.0f), _S185 + make_float2 (_S185.x * _S170 + _S185.y * _S171, 0.0f)));
        float inv_det_0 = 1.0f / (_S186.rows[int(0)].x * _S186.rows[int(1)].y - _S186.rows[int(0)].y * _S186.rows[int(1)].x);
        float _S187 = r_5.x;
        float _S188 = r_5.y;
        float2  q_1 = q_0 - make_float2 ((_S187 * _S186.rows[int(1)].y - _S188 * _S186.rows[int(0)].y) * inv_det_0, (- _S187 * _S186.rows[int(1)].x + _S188 * _S186.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S189 = (*dist_coeffs_6)[int(0)];
    float _S190 = (*dist_coeffs_6)[int(1)];
    float _S191 = (*dist_coeffs_6)[int(2)];
    float _S192 = (*dist_coeffs_6)[int(3)];
    float _S193 = (*dist_coeffs_6)[int(4)];
    float _S194 = (*dist_coeffs_6)[int(5)];
    float _S195 = (*dist_coeffs_6)[int(6)];
    float _S196 = (*dist_coeffs_6)[int(7)];
    float _S197 = (*dist_coeffs_6)[int(8)];
    float _S198 = (*dist_coeffs_6)[int(9)];
    float u_13 = q_0.x;
    float v_13 = q_0.y;
    float _S199 = 0.0f * v_13;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float s_diff_r2_14 = u_13 + u_13 + (_S199 + _S199);
    float _S200 = (*dist_coeffs_6)[int(2)] + r2_13 * (*dist_coeffs_6)[int(3)];
    float _S201 = (*dist_coeffs_6)[int(1)] + r2_13 * _S200;
    float _S202 = (*dist_coeffs_6)[int(0)] + r2_13 * _S201;
    float radial_4 = 1.0f + r2_13 * _S202;
    float _S203 = 2.0f * (*dist_coeffs_6)[int(4)];
    float _S204 = _S203 * u_13;
    float _S205 = 2.0f * u_13;
    float _S206 = 2.0f * (*dist_coeffs_6)[int(5)];
    float _S207 = _S206 * u_13;
    float _S208 = 2.0f * v_13;
    float2  _S209 = make_float2 (1.0f, 0.0f) * make_float2 (radial_4) + make_float2 (s_diff_r2_14 * _S202 + (s_diff_r2_14 * _S201 + (s_diff_r2_14 * _S200 + s_diff_r2_14 * (*dist_coeffs_6)[int(3)] * r2_13) * r2_13) * r2_13) * q_0 + make_float2 (_S203 * v_13 + 0.0f * _S204 + (s_diff_r2_14 + (_S205 + _S205)) * (*dist_coeffs_6)[int(5)] + s_diff_r2_14 * (*dist_coeffs_6)[int(6)], _S206 * v_13 + 0.0f * _S207 + (s_diff_r2_14 + (_S199 + 0.0f * _S208)) * (*dist_coeffs_6)[int(4)] + s_diff_r2_14 * (*dist_coeffs_6)[int(7)]);
    float _S210 = 0.0f * u_13;
    float s_diff_r2_15 = _S210 + _S210 + (v_13 + v_13);
    float2  _S211 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (s_diff_r2_15 * _S202 + (s_diff_r2_15 * _S201 + (s_diff_r2_15 * _S200 + s_diff_r2_15 * (*dist_coeffs_6)[int(3)] * r2_13) * r2_13) * r2_13) * q_0 + make_float2 (0.0f * _S203 * v_13 + _S204 + (s_diff_r2_15 + (_S210 + 0.0f * _S205)) * (*dist_coeffs_6)[int(5)] + s_diff_r2_15 * (*dist_coeffs_6)[int(6)], 0.0f * _S206 * v_13 + _S207 + (s_diff_r2_15 + (_S208 + _S208)) * (*dist_coeffs_6)[int(4)] + s_diff_r2_15 * (*dist_coeffs_6)[int(7)]);
    Matrix<float, 2, 2>  _S212 = transpose_1(makeMatrix<float, 2, 2> (_S209 + make_float2 (_S209.x * (*dist_coeffs_6)[int(8)] + _S209.y * (*dist_coeffs_6)[int(9)], 0.0f), _S211 + make_float2 (_S211.x * (*dist_coeffs_6)[int(8)] + _S211.y * (*dist_coeffs_6)[int(9)], 0.0f)));
    bool _S213;
    if((F32_min((determinant_0(_S212)), ((F32_min((_S212.rows[int(0)].x), (_S212.rows[int(1)].y)))))) > 0.0f)
    {
        float u_14 = (*uv_undist_0).x;
        float v_14 = (*uv_undist_0).y;
        float r2_14 = u_14 * u_14 + v_14 * v_14;
        float2  _S214 = *uv_undist_0 * make_float2 (1.0f + r2_14 * (_S189 + r2_14 * (_S190 + r2_14 * (_S191 + r2_14 * _S192)))) + make_float2 (_S203 * u_14 * v_14 + _S194 * (r2_14 + 2.0f * u_14 * u_14) + _S195 * r2_14, _S206 * u_14 * v_14 + _S193 * (r2_14 + 2.0f * v_14 * v_14) + _S196 * r2_14);
        _S213 = (length_0(_S214 + make_float2 (_S197 * _S214.x + _S198 * _S214.y, 0.0f) - uv_4)) < 0.00999999977648258f;
    }
    else
    {
        _S213 = false;
    }
    return _S213;
}

inline __device__ bool undistort_point(float2  uv_5, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_7, float2  * uv_undist_1)
{
    float2  _S215 = uv_5;
    FixedArray<float, 10>  _S216 = dist_coeffs_7;
    bool _S217 = undistort_point_0(uv_5, &_S216, int(8), &_S215);
    if(!_S217)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S218 = _S215;
        float theta_2 = length_0(_S215);
        float _S219;
        if(theta_2 < 0.00100000004749745f)
        {
            _S219 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S219 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S220 = make_float3 ((_S218 * make_float2 (_S219)).x, (_S218 * make_float2 (_S219)).y, (F32_cos((theta_2))));
        raydir_0 = _S220;
    }
    else
    {
        raydir_0 = make_float3 (_S215.x, _S215.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_6, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_8, float3  * raydir_1)
{
    float2  _S221 = uv_6;
    int3  _S222 = make_int3 (int(0));
    float3  _S223 = make_float3 ((float)_S222.x, (float)_S222.y, (float)_S222.z);
    *raydir_1 = _S223;
    FixedArray<float, 10>  _S224 = dist_coeffs_8;
    bool _S225 = undistort_point_0(uv_6, &_S224, int(8), &_S221);
    if(!_S225)
    {
        return false;
    }
    float3  _S226;
    if(is_fisheye_2)
    {
        float2  _S227 = _S221;
        float theta_3 = length_0(_S221);
        float _S228;
        if(theta_3 < 0.00100000004749745f)
        {
            _S228 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S228 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S229 = make_float3 ((_S227 * make_float2 (_S228)).x, (_S227 * make_float2 (_S228)).y, (F32_cos((theta_3))));
        _S226 = _S229;
    }
    else
    {
        _S226 = make_float3 (_S221.x, _S221.y, 1.0f);
    }
    *raydir_1 = _S226;
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

inline __device__ bool generate_ray(float2  uv_7, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_9, float3  * raydir_2)
{
    float2  _S230 = uv_7;
    FixedArray<float, 10>  _S231 = dist_coeffs_9;
    bool _S232 = undistort_point_0(uv_7, &_S231, int(8), &_S230);
    if(!_S232)
    {
        int3  _S233 = make_int3 (int(0));
        float3  _S234 = make_float3 ((float)_S233.x, (float)_S233.y, (float)_S233.z);
        *raydir_2 = _S234;
        return false;
    }
    float3  _S235;
    if(is_fisheye_3)
    {
        float2  _S236 = _S230;
        float theta_4 = length_0(_S230);
        float _S237;
        if(theta_4 < 0.00100000004749745f)
        {
            _S237 = 1.0f - theta_4 * theta_4 / 6.0f;
        }
        else
        {
            _S237 = (F32_sin((theta_4))) / theta_4;
        }
        float3  _S238 = make_float3 ((_S236 * make_float2 (_S237)).x, (_S236 * make_float2 (_S237)).y, (F32_cos((theta_4))));
        _S235 = _S238;
    }
    else
    {
        _S235 = make_float3 (_S230.x, _S230.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S235);
    return true;
}

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_2)
{
    float _S239 = (*right_2).primal_0.rows[int(0)].x * dOut_2.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_2.x;
    float sum_2 = _S239 + (*right_2).primal_0.rows[int(0)].y * dOut_2.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_2.y;
    float sum_3 = sum_2 + (*right_2).primal_0.rows[int(0)].z * dOut_2.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_2.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_3;
    float _S240 = (*right_2).primal_0.rows[int(1)].x * dOut_2.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_2.x;
    float sum_4 = _S240 + (*right_2).primal_0.rows[int(1)].y * dOut_2.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_2.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(1)].z * dOut_2.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_2.z;
    *&((&left_d_result_1)->y) = sum_5;
    float _S241 = (*right_2).primal_0.rows[int(2)].x * dOut_2.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_2.x;
    float sum_6 = _S241 + (*right_2).primal_0.rows[int(2)].y * dOut_2.y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_2.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(2)].z * dOut_2.z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_2.z;
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

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S242, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S243, float3  _S244)
{
    _d_mul_0(_S242, _S243, _S244);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S245 = - _s_dOut_0;
    float3  _S246 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S247;
    (&_S247)->primal_0 = (*dpt_0).primal_0;
    (&_S247)->differential_0 = _S246;
    Matrix<float, 3, 3>  _S248 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S249;
    (&_S249)->primal_0 = (*dpR_0).primal_0;
    (&_S249)->differential_0 = _S248;
    s_bwd_prop_mul_0(&_S247, &_S249, _S245);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S247.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S249.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S250, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S251, float3  _S252)
{
    s_bwd_prop_transform_ray_o_0(_S250, _S251, _S252);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S253 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S253;
    float3  _S254 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S254;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S255 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S256;
    (&_S256)->primal_0 = (*dpraydir_0).primal_0;
    (&_S256)->differential_0 = _S255;
    Matrix<float, 3, 3>  _S257 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S258;
    (&_S258)->primal_0 = (*dpR_1).primal_0;
    (&_S258)->differential_0 = _S257;
    s_bwd_prop_mul_0(&_S256, &_S258, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S256.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S258.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S259, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S260, float3  _S261)
{
    s_bwd_prop_transform_ray_d_0(_S259, _S260, _S261);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S262 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S262;
    float3  _S263 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S263;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_2, float3  scale_1)
{
    float x_12 = quat_2.y;
    float x2_2 = x_12 * x_12;
    float y2_2 = quat_2.z * quat_2.z;
    float z2_2 = quat_2.w * quat_2.w;
    float xy_2 = quat_2.y * quat_2.z;
    float xz_2 = quat_2.y * quat_2.w;
    float yz_2 = quat_2.z * quat_2.w;
    float wx_2 = quat_2.x * quat_2.y;
    float wy_2 = quat_2.x * quat_2.z;
    float wz_2 = quat_2.x * quat_2.w;
    return mul_1(makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)))));
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S264, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S265, Matrix<float, 3, 3>  _S266)
{
    mul_0(_S264, _S265, _S266);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_2)
{
    float _S267 = (*dpquat_0).primal_0.y;
    float x2_3 = _S267 * _S267;
    float y2_3 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_3 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_3 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_3 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_3 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S268 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))));
    Matrix<float, 3, 3>  _S269 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S270;
    (&_S270)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S270)->differential_0 = _S269;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S271;
    (&_S271)->primal_0 = _S268;
    (&_S271)->differential_0 = _S269;
    s_bwd_prop_mul_1(&_S270, &_S271, _s_dOut_2);
    Matrix<float, 3, 3>  _S272 = transpose_0(transpose_0(_S271.differential_0));
    float _S273 = 2.0f * - _S272.rows[int(2)].z;
    float _S274 = 2.0f * _S272.rows[int(2)].y;
    float _S275 = 2.0f * _S272.rows[int(2)].x;
    float _S276 = 2.0f * _S272.rows[int(1)].z;
    float _S277 = 2.0f * - _S272.rows[int(1)].y;
    float _S278 = 2.0f * _S272.rows[int(1)].x;
    float _S279 = 2.0f * _S272.rows[int(0)].z;
    float _S280 = 2.0f * _S272.rows[int(0)].y;
    float _S281 = 2.0f * - _S272.rows[int(0)].x;
    float _S282 = - _S278 + _S280;
    float _S283 = _S275 + - _S279;
    float _S284 = - _S274 + _S276;
    float _S285 = _S274 + _S276;
    float _S286 = _S275 + _S279;
    float _S287 = _S278 + _S280;
    float _S288 = (*dpquat_0).primal_0.w * (_S277 + _S281);
    float _S289 = (*dpquat_0).primal_0.z * (_S273 + _S281);
    float _S290 = (*dpquat_0).primal_0.y * (_S273 + _S277);
    float _S291 = (*dpquat_0).primal_0.x * _S282 + (*dpquat_0).primal_0.z * _S285 + (*dpquat_0).primal_0.y * _S286 + _S288 + _S288;
    float _S292 = (*dpquat_0).primal_0.x * _S283 + (*dpquat_0).primal_0.w * _S285 + (*dpquat_0).primal_0.y * _S287 + _S289 + _S289;
    float _S293 = (*dpquat_0).primal_0.x * _S284 + (*dpquat_0).primal_0.w * _S286 + (*dpquat_0).primal_0.z * _S287 + _S290 + _S290;
    float _S294 = (*dpquat_0).primal_0.w * _S282 + (*dpquat_0).primal_0.z * _S283 + (*dpquat_0).primal_0.y * _S284;
    float3  _S295 = make_float3 (_S270.differential_0.rows[int(0)].x, _S270.differential_0.rows[int(1)].y, _S270.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S295;
    float4  _S296 = make_float4 (0.0f);
    *&((&_S296)->w) = _S291;
    *&((&_S296)->z) = _S292;
    *&((&_S296)->y) = _S293;
    *&((&_S296)->x) = _S294;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S296;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S297, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S298, Matrix<float, 3, 3>  _S299)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S297, _S298, _S299);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_3, float3  scale_2, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_0, float3  * v_scale_0)
{
    float4  _S300 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_3;
    (&dp_quat_0)->differential_0 = _S300;
    float3  _S301 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_2;
    (&dp_scale_0)->differential_0 = _S301;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_0 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_4, float3  dOut_3)
{
    float _S302 = (*left_4).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = (*right_4).primal_0.x * dOut_3.x;
    float sum_10 = _S302 + (*left_4).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = (*right_4).primal_0.x * dOut_3.y;
    float sum_11 = sum_10 + (*left_4).primal_0.rows[int(2)].x * dOut_3.z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = (*right_4).primal_0.x * dOut_3.z;
    float3  right_d_result_2;
    *&((&right_d_result_2)->x) = sum_11;
    float _S303 = (*left_4).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = (*right_4).primal_0.y * dOut_3.x;
    float sum_12 = _S303 + (*left_4).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = (*right_4).primal_0.y * dOut_3.y;
    float sum_13 = sum_12 + (*left_4).primal_0.rows[int(2)].y * dOut_3.z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = (*right_4).primal_0.y * dOut_3.z;
    *&((&right_d_result_2)->y) = sum_13;
    float _S304 = (*left_4).primal_0.rows[int(0)].z * dOut_3.x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = (*right_4).primal_0.z * dOut_3.x;
    float sum_14 = _S304 + (*left_4).primal_0.rows[int(1)].z * dOut_3.y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = (*right_4).primal_0.z * dOut_3.y;
    float sum_15 = sum_14 + (*left_4).primal_0.rows[int(2)].z * dOut_3.z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = (*right_4).primal_0.z * dOut_3.z;
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

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_4)
{
    float _S305 = dOut_4.y;
    float _S306 = dOut_4.z;
    float _S307 = dOut_4.x;
    float _S308 = (*a_0).primal_0.z * _S305 + - (*a_0).primal_0.y * _S306;
    float _S309 = - (*a_0).primal_0.z * _S307 + (*a_0).primal_0.x * _S306;
    float _S310 = (*a_0).primal_0.y * _S307 + - (*a_0).primal_0.x * _S305;
    float3  _S311 = make_float3 (- (*b_0).primal_0.z * _S305 + (*b_0).primal_0.y * _S306, (*b_0).primal_0.z * _S307 + - (*b_0).primal_0.x * _S306, - (*b_0).primal_0.y * _S307 + (*b_0).primal_0.x * _S305);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S311;
    float3  _S312 = make_float3 (_S308, _S309, _S310);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S312;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S313 = left_6.y;
    float _S314 = right_6.z;
    float _S315 = left_6.z;
    float _S316 = right_6.y;
    float _S317 = right_6.x;
    float _S318 = left_6.x;
    return make_float3 (_S313 * _S314 - _S315 * _S316, _S315 * _S317 - _S318 * _S314, _S318 * _S316 - _S313 * _S317);
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_4, float dOut_5)
{
    float _S319 = (F32_exp(((*dpx_4).primal_0))) * dOut_5;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S319;
    return;
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_0, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_3(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_3(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_0 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S320, float3  _S321)
{
    return mul_3(_S320, _S321);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S322, float3  _S323)
{
    return cross_0(_S322, _S323);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S324, float3  _S325)
{
    return dot_0(_S324, _S325);
}

inline __device__ float s_primal_ctx_exp_0(float _S326)
{
    return (F32_exp((_S326)));
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S327, float _S328)
{
    _d_exp_0(_S327, _S328);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S329, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S330, float _S331)
{
    _d_dot_0(_S329, _S330, _S331);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S332, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S333, float3  _S334)
{
    _d_cross_0(_S332, _S333, _S334);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S335, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S336, float3  _S337)
{
    _d_mul_1(_S335, _S336, _S337);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    float3  _S338 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S339 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S338);
    float3  _S340 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S341 = s_primal_ctx_cross_0(_S340, _S339);
    float _S342 = -0.5f * s_primal_ctx_dot_0(_S341, _S341);
    float _S343 = s_primal_ctx_dot_0(_S340, _S340);
    float _S344 = _S342 / _S343;
    float _S345 = _S343 * _S343;
    float _S346 = (*dpopacity_0).primal_0 * _s_dOut_3;
    float _S347 = s_primal_ctx_exp_0(_S344) * _s_dOut_3;
    DiffPair_float_0 _S348;
    (&_S348)->primal_0 = _S344;
    (&_S348)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S348, _S346);
    float _S349 = _S348.differential_0 / _S345;
    float _S350 = _S342 * - _S349;
    float _S351 = _S343 * _S349;
    float3  _S352 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S353;
    (&_S353)->primal_0 = _S340;
    (&_S353)->differential_0 = _S352;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S354;
    (&_S354)->primal_0 = _S340;
    (&_S354)->differential_0 = _S352;
    s_bwd_prop_dot_0(&_S353, &_S354, _S350);
    float _S355 = -0.5f * _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S356;
    (&_S356)->primal_0 = _S341;
    (&_S356)->differential_0 = _S352;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S357;
    (&_S357)->primal_0 = _S341;
    (&_S357)->differential_0 = _S352;
    s_bwd_prop_dot_0(&_S356, &_S357, _S355);
    float3  _S358 = _S357.differential_0 + _S356.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S359;
    (&_S359)->primal_0 = _S340;
    (&_S359)->differential_0 = _S352;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S360;
    (&_S360)->primal_0 = _S339;
    (&_S360)->differential_0 = _S352;
    s_bwd_prop_cross_0(&_S359, &_S360, _S358);
    float3  _S361 = _S354.differential_0 + _S353.differential_0 + _S359.differential_0;
    Matrix<float, 3, 3>  _S362 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S363;
    (&_S363)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S363)->differential_0 = _S362;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S364;
    (&_S364)->primal_0 = (*dpray_d_0).primal_0;
    (&_S364)->differential_0 = _S352;
    s_bwd_prop_mul_2(&_S363, &_S364, _S361);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S365;
    (&_S365)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S365)->differential_0 = _S362;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S366;
    (&_S366)->primal_0 = _S338;
    (&_S366)->differential_0 = _S352;
    s_bwd_prop_mul_2(&_S365, &_S366, _S360.differential_0);
    float3  _S367 = - _S366.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S364.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S366.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S347;
    Matrix<float, 3, 3>  _S368 = _S363.differential_0 + _S365.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S368;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S367;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S369, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S370, DiffPair_float_0 * _S371, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S372, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S373, float _S374)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S369, _S370, _S371, _S372, _S373, _S374);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_1, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_0, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S375 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S375;
    Matrix<float, 3, 3>  _S376 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S376;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_1;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S375;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S375;
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
    float3  _S377 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S378 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S377);
    float3  _S379 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S380 = s_primal_ctx_dot_0(_S379, _S379);
    float _S381 = dpdepth_0 / (_S380 * _S380);
    float _S382 = - s_primal_ctx_dot_0(_S378, _S379) * - _S381;
    float _S383 = _S380 * _S381;
    float3  _S384 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S385;
    (&_S385)->primal_0 = _S379;
    (&_S385)->differential_0 = _S384;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S386;
    (&_S386)->primal_0 = _S379;
    (&_S386)->differential_0 = _S384;
    s_bwd_prop_dot_0(&_S385, &_S386, _S382);
    float _S387 = - _S383;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S388;
    (&_S388)->primal_0 = _S378;
    (&_S388)->differential_0 = _S384;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S389;
    (&_S389)->primal_0 = _S379;
    (&_S389)->differential_0 = _S384;
    s_bwd_prop_dot_0(&_S388, &_S389, _S387);
    float3  _S390 = _S386.differential_0 + _S385.differential_0 + _S389.differential_0;
    Matrix<float, 3, 3>  _S391 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S392;
    (&_S392)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S392)->differential_0 = _S391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S393;
    (&_S393)->primal_0 = (*dpray_d_1).primal_0;
    (&_S393)->differential_0 = _S384;
    s_bwd_prop_mul_2(&_S392, &_S393, _S390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S394;
    (&_S394)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S394)->differential_0 = _S391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S395;
    (&_S395)->primal_0 = _S377;
    (&_S395)->differential_0 = _S384;
    s_bwd_prop_mul_2(&_S394, &_S395, _S388.differential_0);
    float3  _S396 = - _S395.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S393.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S395.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S397 = _S392.differential_0 + _S394.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S397;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S396;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S398, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S399, DiffPair_float_0 * _S400, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S401, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S402, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S403, float3  _S404, float _S405)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S398, _S399, _S400, _S401, _S402, _S403, _S404, _S405);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_3, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_1, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S406 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S406;
    Matrix<float, 3, 3>  _S407 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S407;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_3;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S406;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S406;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S406;
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
    float _S408 = scale_3.x;
    float sx_0 = (F32_exp((_S408)));
    float _S409 = scale_3.y;
    float sy_0 = (F32_exp((_S409)));
    float sz_0 = scale_3.z - 0.5f * (_S408 + _S409);
    float4  _S410 = normalize_1(quat_4);
    float x_13 = _S410.y;
    float x2_4 = x_13 * x_13;
    float y2_4 = _S410.z * _S410.z;
    float z2_4 = _S410.w * _S410.w;
    float xy_4 = _S410.y * _S410.z;
    float xz_4 = _S410.y * _S410.w;
    float yz_4 = _S410.z * _S410.w;
    float wx_4 = _S410.x * _S410.y;
    float wy_4 = _S410.x * _S410.z;
    float wz_4 = _S410.x * _S410.w;
    Matrix<float, 3, 3>  _S411 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    *vert0_0 = mul_3(_S411, make_float3 (sx_0, 0.0f, 0.0f)) + mean_4;
    *vert1_0 = mul_3(_S411, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_4;
    *vert2_0 = mul_3(_S411, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_4;
    return;
}

