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

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, float dOut_1)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_0).primal_0.x * dOut_1;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_0).primal_0.x * dOut_1;
    *&((&x_d_result_0)->y) = (*dpy_0).primal_0.y * dOut_1;
    *&((&y_d_result_0)->y) = (*dpx_0).primal_0.y * dOut_1;
    *&((&x_d_result_0)->z) = (*dpy_0).primal_0.z * dOut_1;
    *&((&y_d_result_0)->z) = (*dpx_0).primal_0.z * dOut_1;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = x_d_result_0;
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

inline __device__ float length_0(float2  x_6)
{
    return (F32_sqrt((dot_1(x_6, x_6))));
}

inline __device__ float length_1(float3  x_7)
{
    return (F32_sqrt((dot_0(x_7, x_7))));
}

inline __device__ float2  distort_point(float2  uv_1, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_1)
{
    float2  _S15;
    if(is_fisheye_0)
    {
        float r_3 = length_0(uv_1);
        float theta_0 = (F32_atan((r_3)));
        float _S16;
        if(r_3 < 0.00100000004749745f)
        {
            _S16 = 1.0f - r_3 * r_3 / 3.0f;
        }
        else
        {
            _S16 = theta_0 / r_3;
        }
        _S15 = uv_1 * make_float2 (_S16);
    }
    else
    {
        _S15 = uv_1;
    }
    float u_1 = _S15.x;
    float v_1 = _S15.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S17 = _S15 * make_float2 (1.0f + r2_1 * (dist_coeffs_1[int(0)] + r2_1 * (dist_coeffs_1[int(1)] + r2_1 * (dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_1 * v_1 + dist_coeffs_1[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_1[int(6)] * r2_1, 2.0f * dist_coeffs_1[int(5)] * u_1 * v_1 + dist_coeffs_1[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_1[int(7)] * r2_1);
    return _S17 + make_float2 (dist_coeffs_1[int(8)] * _S17.x + dist_coeffs_1[int(9)] * _S17.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_2, FixedArray<float, 10>  * dist_coeffs_2, int maxiter_0, float2  * uv_undist_0)
{
    int i_3 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_3 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S18 = (*dist_coeffs_2)[int(3)];
        float _S19 = (*dist_coeffs_2)[int(4)];
        float _S20 = (*dist_coeffs_2)[int(5)];
        float _S21 = (*dist_coeffs_2)[int(6)];
        float _S22 = (*dist_coeffs_2)[int(7)];
        float _S23 = (*dist_coeffs_2)[int(8)];
        float _S24 = (*dist_coeffs_2)[int(9)];
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float _S25 = (*dist_coeffs_2)[int(2)] + r2_2 * (*dist_coeffs_2)[int(3)];
        float _S26 = (*dist_coeffs_2)[int(1)] + r2_2 * _S25;
        float _S27 = (*dist_coeffs_2)[int(0)] + r2_2 * _S26;
        float radial_1 = 1.0f + r2_2 * _S27;
        float _S28 = 2.0f * (*dist_coeffs_2)[int(4)];
        float _S29 = _S28 * u_2;
        float _S30 = 2.0f * u_2;
        float _S31 = 2.0f * (*dist_coeffs_2)[int(5)];
        float _S32 = _S31 * u_2;
        float _S33 = 2.0f * v_2;
        float2  _S34 = q_0 * make_float2 (radial_1) + make_float2 (_S29 * v_2 + (*dist_coeffs_2)[int(5)] * (r2_2 + _S30 * u_2) + (*dist_coeffs_2)[int(6)] * r2_2, _S32 * v_2 + (*dist_coeffs_2)[int(4)] * (r2_2 + _S33 * v_2) + (*dist_coeffs_2)[int(7)] * r2_2);
        float2  r_4 = _S34 + make_float2 ((*dist_coeffs_2)[int(8)] * _S34.x + (*dist_coeffs_2)[int(9)] * _S34.y, 0.0f) - uv_2;
        float _S35 = 0.0f * v_2;
        float s_diff_r2_2 = u_2 + u_2 + (_S35 + _S35);
        float2  _S36 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S27 + (s_diff_r2_2 * _S26 + (s_diff_r2_2 * _S25 + s_diff_r2_2 * _S18 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (_S28 * v_2 + 0.0f * _S29 + (s_diff_r2_2 + (_S30 + _S30)) * _S20 + s_diff_r2_2 * _S21, _S31 * v_2 + 0.0f * _S32 + (s_diff_r2_2 + (_S35 + 0.0f * _S33)) * _S19 + s_diff_r2_2 * _S22);
        float _S37 = 0.0f * u_2;
        float s_diff_r2_3 = _S37 + _S37 + (v_2 + v_2);
        float2  _S38 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S27 + (s_diff_r2_3 * _S26 + (s_diff_r2_3 * _S25 + s_diff_r2_3 * _S18 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (0.0f * _S28 * v_2 + _S29 + (s_diff_r2_3 + (_S37 + 0.0f * _S30)) * _S20 + s_diff_r2_3 * _S21, 0.0f * _S31 * v_2 + _S32 + (s_diff_r2_3 + (_S33 + _S33)) * _S19 + s_diff_r2_3 * _S22);
        Matrix<float, 2, 2>  _S39 = transpose_1(makeMatrix<float, 2, 2> (_S36 + make_float2 (_S36.x * _S23 + _S36.y * _S24, 0.0f), _S38 + make_float2 (_S38.x * _S23 + _S38.y * _S24, 0.0f)));
        float inv_det_0 = 1.0f / (_S39.rows[int(0)].x * _S39.rows[int(1)].y - _S39.rows[int(0)].y * _S39.rows[int(1)].x);
        float _S40 = r_4.x;
        float _S41 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S40 * _S39.rows[int(1)].y - _S41 * _S39.rows[int(0)].y) * inv_det_0, (- _S40 * _S39.rows[int(1)].x + _S41 * _S39.rows[int(0)].x) * inv_det_0);
        i_3 = i_3 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S42 = (*dist_coeffs_2)[int(0)];
    float _S43 = (*dist_coeffs_2)[int(1)];
    float _S44 = (*dist_coeffs_2)[int(2)];
    float _S45 = (*dist_coeffs_2)[int(3)];
    float _S46 = (*dist_coeffs_2)[int(4)];
    float _S47 = (*dist_coeffs_2)[int(5)];
    float _S48 = (*dist_coeffs_2)[int(6)];
    float _S49 = (*dist_coeffs_2)[int(7)];
    float _S50 = (*dist_coeffs_2)[int(8)];
    float _S51 = (*dist_coeffs_2)[int(9)];
    float u_3 = q_0.x;
    float v_3 = q_0.y;
    float _S52 = 0.0f * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = u_3 + u_3 + (_S52 + _S52);
    float _S53 = (*dist_coeffs_2)[int(2)] + r2_3 * (*dist_coeffs_2)[int(3)];
    float _S54 = (*dist_coeffs_2)[int(1)] + r2_3 * _S53;
    float _S55 = (*dist_coeffs_2)[int(0)] + r2_3 * _S54;
    float radial_2 = 1.0f + r2_3 * _S55;
    float _S56 = 2.0f * (*dist_coeffs_2)[int(4)];
    float _S57 = _S56 * u_3;
    float _S58 = 2.0f * u_3;
    float _S59 = 2.0f * (*dist_coeffs_2)[int(5)];
    float _S60 = _S59 * u_3;
    float _S61 = 2.0f * v_3;
    float2  _S62 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S55 + (s_diff_r2_4 * _S54 + (s_diff_r2_4 * _S53 + s_diff_r2_4 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (_S56 * v_3 + 0.0f * _S57 + (s_diff_r2_4 + (_S58 + _S58)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_4 * (*dist_coeffs_2)[int(6)], _S59 * v_3 + 0.0f * _S60 + (s_diff_r2_4 + (_S52 + 0.0f * _S61)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_4 * (*dist_coeffs_2)[int(7)]);
    float _S63 = 0.0f * u_3;
    float s_diff_r2_5 = _S63 + _S63 + (v_3 + v_3);
    float2  _S64 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S55 + (s_diff_r2_5 * _S54 + (s_diff_r2_5 * _S53 + s_diff_r2_5 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (0.0f * _S56 * v_3 + _S57 + (s_diff_r2_5 + (_S63 + 0.0f * _S58)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_5 * (*dist_coeffs_2)[int(6)], 0.0f * _S59 * v_3 + _S60 + (s_diff_r2_5 + (_S61 + _S61)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_5 * (*dist_coeffs_2)[int(7)]);
    Matrix<float, 2, 2>  _S65 = transpose_1(makeMatrix<float, 2, 2> (_S62 + make_float2 (_S62.x * (*dist_coeffs_2)[int(8)] + _S62.y * (*dist_coeffs_2)[int(9)], 0.0f), _S64 + make_float2 (_S64.x * (*dist_coeffs_2)[int(8)] + _S64.y * (*dist_coeffs_2)[int(9)], 0.0f)));
    bool _S66;
    if((F32_min((determinant_0(_S65)), ((F32_min((_S65.rows[int(0)].x), (_S65.rows[int(1)].y)))))) > 0.0f)
    {
        float u_4 = (*uv_undist_0).x;
        float v_4 = (*uv_undist_0).y;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float2  _S67 = *uv_undist_0 * make_float2 (1.0f + r2_4 * (_S42 + r2_4 * (_S43 + r2_4 * (_S44 + r2_4 * _S45)))) + make_float2 (_S56 * u_4 * v_4 + _S47 * (r2_4 + 2.0f * u_4 * u_4) + _S48 * r2_4, _S59 * u_4 * v_4 + _S46 * (r2_4 + 2.0f * v_4 * v_4) + _S49 * r2_4);
        _S66 = (length_0(_S67 + make_float2 (_S50 * _S67.x + _S51 * _S67.y, 0.0f) - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S66 = false;
    }
    return _S66;
}

inline __device__ bool undistort_point(float2  uv_3, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_undist_1)
{
    float2  _S68 = uv_3;
    FixedArray<float, 10>  _S69 = dist_coeffs_3;
    bool _S70 = undistort_point_0(uv_3, &_S69, int(8), &_S68);
    if(!_S70)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S71 = _S68;
        float theta_1 = length_0(_S68);
        float _S72;
        if(theta_1 < 0.00100000004749745f)
        {
            _S72 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S72 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S73 = make_float3 ((_S71 * make_float2 (_S72)).x, (_S71 * make_float2 (_S72)).y, (F32_cos((theta_1))));
        raydir_0 = _S73;
    }
    else
    {
        raydir_0 = make_float3 (_S68.x, _S68.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_4, float3  * raydir_1)
{
    float2  _S74 = uv_4;
    int3  _S75 = make_int3 (int(0));
    float3  _S76 = make_float3 ((float)_S75.x, (float)_S75.y, (float)_S75.z);
    *raydir_1 = _S76;
    FixedArray<float, 10>  _S77 = dist_coeffs_4;
    bool _S78 = undistort_point_0(uv_4, &_S77, int(8), &_S74);
    if(!_S78)
    {
        return false;
    }
    float3  _S79;
    if(is_fisheye_2)
    {
        float2  _S80 = _S74;
        float theta_2 = length_0(_S74);
        float _S81;
        if(theta_2 < 0.00100000004749745f)
        {
            _S81 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S81 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S82 = make_float3 ((_S80 * make_float2 (_S81)).x, (_S80 * make_float2 (_S81)).y, (F32_cos((theta_2))));
        _S79 = _S82;
    }
    else
    {
        _S79 = make_float3 (_S74.x, _S74.y, 1.0f);
    }
    *raydir_1 = _S79;
    return true;
}

inline __device__ float3  normalize_0(float3  x_8)
{
    return x_8 / make_float3 (length_1(x_8));
}

inline __device__ bool generate_ray(float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_5, float3  * raydir_2)
{
    float2  _S83 = uv_5;
    FixedArray<float, 10>  _S84 = dist_coeffs_5;
    bool _S85 = undistort_point_0(uv_5, &_S84, int(8), &_S83);
    if(!_S85)
    {
        int3  _S86 = make_int3 (int(0));
        float3  _S87 = make_float3 ((float)_S86.x, (float)_S86.y, (float)_S86.z);
        *raydir_2 = _S87;
        return false;
    }
    float3  _S88;
    if(is_fisheye_3)
    {
        float2  _S89 = _S83;
        float theta_3 = length_0(_S83);
        float _S90;
        if(theta_3 < 0.00100000004749745f)
        {
            _S90 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S90 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S91 = make_float3 ((_S89 * make_float2 (_S90)).x, (_S89 * make_float2 (_S90)).y, (F32_cos((theta_3))));
        _S88 = _S91;
    }
    else
    {
        _S88 = make_float3 (_S83.x, _S83.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S88);
    return true;
}

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_2)
{
    float _S92 = (*right_2).primal_0.rows[int(0)].x * dOut_2.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_2.x;
    float sum_2 = _S92 + (*right_2).primal_0.rows[int(0)].y * dOut_2.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_2.y;
    float sum_3 = sum_2 + (*right_2).primal_0.rows[int(0)].z * dOut_2.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_2.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_3;
    float _S93 = (*right_2).primal_0.rows[int(1)].x * dOut_2.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_2.x;
    float sum_4 = _S93 + (*right_2).primal_0.rows[int(1)].y * dOut_2.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_2.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(1)].z * dOut_2.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_2.z;
    *&((&left_d_result_1)->y) = sum_5;
    float _S94 = (*right_2).primal_0.rows[int(2)].x * dOut_2.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_2.x;
    float sum_6 = _S94 + (*right_2).primal_0.rows[int(2)].y * dOut_2.y;
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
    float3  result_7;
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
        int i_4 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(i_4 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_3, i_4) * _slang_vector_get_element(right_3.rows[i_4], j_0);
            i_4 = i_4 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_7, j_0) = sum_8;
        j_0 = j_0 + int(1);
    }
    return result_7;
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

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S95, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S96, float3  _S97)
{
    _d_mul_0(_S95, _S96, _S97);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S98 = - _s_dOut_0;
    float3  _S99 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S100;
    (&_S100)->primal_0 = (*dpt_0).primal_0;
    (&_S100)->differential_0 = _S99;
    Matrix<float, 3, 3>  _S101 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S102;
    (&_S102)->primal_0 = (*dpR_0).primal_0;
    (&_S102)->differential_0 = _S101;
    s_bwd_prop_mul_0(&_S100, &_S102, _S98);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S100.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S102.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S103, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S104, float3  _S105)
{
    s_bwd_prop_transform_ray_o_0(_S103, _S104, _S105);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S106 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S106;
    float3  _S107 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S107;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S108 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S109;
    (&_S109)->primal_0 = (*dpraydir_0).primal_0;
    (&_S109)->differential_0 = _S108;
    Matrix<float, 3, 3>  _S110 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S111;
    (&_S111)->primal_0 = (*dpR_1).primal_0;
    (&_S111)->differential_0 = _S110;
    s_bwd_prop_mul_0(&_S109, &_S111, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S109.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S111.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S112, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S113, float3  _S114)
{
    s_bwd_prop_transform_ray_d_0(_S112, _S113, _S114);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S115 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S115;
    float3  _S116 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S116;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_2, float3  scale_1)
{
    float x_9 = quat_2.y;
    float x2_2 = x_9 * x_9;
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

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S117, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S118, Matrix<float, 3, 3>  _S119)
{
    mul_0(_S117, _S118, _S119);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_2)
{
    float _S120 = (*dpquat_0).primal_0.y;
    float x2_3 = _S120 * _S120;
    float y2_3 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_3 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_3 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_3 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_3 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_3 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S121 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))));
    Matrix<float, 3, 3>  _S122 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S123;
    (&_S123)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S123)->differential_0 = _S122;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S124;
    (&_S124)->primal_0 = _S121;
    (&_S124)->differential_0 = _S122;
    s_bwd_prop_mul_1(&_S123, &_S124, _s_dOut_2);
    Matrix<float, 3, 3>  _S125 = transpose_0(transpose_0(_S124.differential_0));
    float _S126 = 2.0f * - _S125.rows[int(2)].z;
    float _S127 = 2.0f * _S125.rows[int(2)].y;
    float _S128 = 2.0f * _S125.rows[int(2)].x;
    float _S129 = 2.0f * _S125.rows[int(1)].z;
    float _S130 = 2.0f * - _S125.rows[int(1)].y;
    float _S131 = 2.0f * _S125.rows[int(1)].x;
    float _S132 = 2.0f * _S125.rows[int(0)].z;
    float _S133 = 2.0f * _S125.rows[int(0)].y;
    float _S134 = 2.0f * - _S125.rows[int(0)].x;
    float _S135 = - _S131 + _S133;
    float _S136 = _S128 + - _S132;
    float _S137 = - _S127 + _S129;
    float _S138 = _S127 + _S129;
    float _S139 = _S128 + _S132;
    float _S140 = _S131 + _S133;
    float _S141 = (*dpquat_0).primal_0.w * (_S130 + _S134);
    float _S142 = (*dpquat_0).primal_0.z * (_S126 + _S134);
    float _S143 = (*dpquat_0).primal_0.y * (_S126 + _S130);
    float _S144 = (*dpquat_0).primal_0.x * _S135 + (*dpquat_0).primal_0.z * _S138 + (*dpquat_0).primal_0.y * _S139 + _S141 + _S141;
    float _S145 = (*dpquat_0).primal_0.x * _S136 + (*dpquat_0).primal_0.w * _S138 + (*dpquat_0).primal_0.y * _S140 + _S142 + _S142;
    float _S146 = (*dpquat_0).primal_0.x * _S137 + (*dpquat_0).primal_0.w * _S139 + (*dpquat_0).primal_0.z * _S140 + _S143 + _S143;
    float _S147 = (*dpquat_0).primal_0.w * _S135 + (*dpquat_0).primal_0.z * _S136 + (*dpquat_0).primal_0.y * _S137;
    float3  _S148 = make_float3 (_S123.differential_0.rows[int(0)].x, _S123.differential_0.rows[int(1)].y, _S123.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S148;
    float4  _S149 = make_float4 (0.0f);
    *&((&_S149)->w) = _S144;
    *&((&_S149)->z) = _S145;
    *&((&_S149)->y) = _S146;
    *&((&_S149)->x) = _S147;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S149;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S150, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S151, Matrix<float, 3, 3>  _S152)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S150, _S151, _S152);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_3, float3  scale_2, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_0, float3  * v_scale_0)
{
    float4  _S153 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_3;
    (&dp_quat_0)->differential_0 = _S153;
    float3  _S154 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_2;
    (&dp_scale_0)->differential_0 = _S154;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_0 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_4, float3  dOut_3)
{
    float _S155 = (*left_4).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = (*right_4).primal_0.x * dOut_3.x;
    float sum_10 = _S155 + (*left_4).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = (*right_4).primal_0.x * dOut_3.y;
    float sum_11 = sum_10 + (*left_4).primal_0.rows[int(2)].x * dOut_3.z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = (*right_4).primal_0.x * dOut_3.z;
    float3  right_d_result_2;
    *&((&right_d_result_2)->x) = sum_11;
    float _S156 = (*left_4).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = (*right_4).primal_0.y * dOut_3.x;
    float sum_12 = _S156 + (*left_4).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = (*right_4).primal_0.y * dOut_3.y;
    float sum_13 = sum_12 + (*left_4).primal_0.rows[int(2)].y * dOut_3.z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = (*right_4).primal_0.y * dOut_3.z;
    *&((&right_d_result_2)->y) = sum_13;
    float _S157 = (*left_4).primal_0.rows[int(0)].z * dOut_3.x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = (*right_4).primal_0.z * dOut_3.x;
    float sum_14 = _S157 + (*left_4).primal_0.rows[int(1)].z * dOut_3.y;
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
    float3  result_8;
    int i_5 = int(0);
    for(;;)
    {
        if(i_5 < int(3))
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
            float sum_17 = sum_16 + _slang_vector_get_element(left_5.rows[i_5], j_1) * _slang_vector_get_element(right_5, j_1);
            j_1 = j_1 + int(1);
            sum_16 = sum_17;
        }
        *_slang_vector_get_element_ptr(&result_8, i_5) = sum_16;
        i_5 = i_5 + int(1);
    }
    return result_8;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_4)
{
    float _S158 = dOut_4.y;
    float _S159 = dOut_4.z;
    float _S160 = dOut_4.x;
    float _S161 = (*a_0).primal_0.z * _S158 + - (*a_0).primal_0.y * _S159;
    float _S162 = - (*a_0).primal_0.z * _S160 + (*a_0).primal_0.x * _S159;
    float _S163 = (*a_0).primal_0.y * _S160 + - (*a_0).primal_0.x * _S158;
    float3  _S164 = make_float3 (- (*b_0).primal_0.z * _S158 + (*b_0).primal_0.y * _S159, (*b_0).primal_0.z * _S160 + - (*b_0).primal_0.x * _S159, - (*b_0).primal_0.y * _S160 + (*b_0).primal_0.x * _S158);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S164;
    float3  _S165 = make_float3 (_S161, _S162, _S163);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S165;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S166 = left_6.y;
    float _S167 = right_6.z;
    float _S168 = left_6.z;
    float _S169 = right_6.y;
    float _S170 = right_6.x;
    float _S171 = left_6.x;
    return make_float3 (_S166 * _S167 - _S168 * _S169, _S168 * _S170 - _S171 * _S167, _S171 * _S169 - _S166 * _S170);
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_1, float dOut_5)
{
    float _S172 = (F32_exp(((*dpx_1).primal_0))) * dOut_5;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S172;
    return;
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_0, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_3(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_3(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_0 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S173, float3  _S174)
{
    return mul_3(_S173, _S174);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S175, float3  _S176)
{
    return cross_0(_S175, _S176);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S177, float3  _S178)
{
    return dot_0(_S177, _S178);
}

inline __device__ float s_primal_ctx_exp_0(float _S179)
{
    return (F32_exp((_S179)));
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S180, float _S181)
{
    _d_exp_0(_S180, _S181);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S182, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S183, float _S184)
{
    _d_dot_0(_S182, _S183, _S184);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S185, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S186, float3  _S187)
{
    _d_cross_0(_S185, _S186, _S187);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S188, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S189, float3  _S190)
{
    _d_mul_1(_S188, _S189, _S190);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    float3  _S191 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S192 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S191);
    float3  _S193 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S194 = s_primal_ctx_cross_0(_S193, _S192);
    float _S195 = -0.5f * s_primal_ctx_dot_0(_S194, _S194);
    float _S196 = s_primal_ctx_dot_0(_S193, _S193);
    float _S197 = _S195 / _S196;
    float _S198 = _S196 * _S196;
    float _S199 = (*dpopacity_0).primal_0 * _s_dOut_3;
    float _S200 = s_primal_ctx_exp_0(_S197) * _s_dOut_3;
    DiffPair_float_0 _S201;
    (&_S201)->primal_0 = _S197;
    (&_S201)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S201, _S199);
    float _S202 = _S201.differential_0 / _S198;
    float _S203 = _S195 * - _S202;
    float _S204 = _S196 * _S202;
    float3  _S205 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S206;
    (&_S206)->primal_0 = _S193;
    (&_S206)->differential_0 = _S205;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S207;
    (&_S207)->primal_0 = _S193;
    (&_S207)->differential_0 = _S205;
    s_bwd_prop_dot_0(&_S206, &_S207, _S203);
    float _S208 = -0.5f * _S204;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S209;
    (&_S209)->primal_0 = _S194;
    (&_S209)->differential_0 = _S205;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S210;
    (&_S210)->primal_0 = _S194;
    (&_S210)->differential_0 = _S205;
    s_bwd_prop_dot_0(&_S209, &_S210, _S208);
    float3  _S211 = _S210.differential_0 + _S209.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S212;
    (&_S212)->primal_0 = _S193;
    (&_S212)->differential_0 = _S205;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S213;
    (&_S213)->primal_0 = _S192;
    (&_S213)->differential_0 = _S205;
    s_bwd_prop_cross_0(&_S212, &_S213, _S211);
    float3  _S214 = _S207.differential_0 + _S206.differential_0 + _S212.differential_0;
    Matrix<float, 3, 3>  _S215 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S216;
    (&_S216)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S216)->differential_0 = _S215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S217;
    (&_S217)->primal_0 = (*dpray_d_0).primal_0;
    (&_S217)->differential_0 = _S205;
    s_bwd_prop_mul_2(&_S216, &_S217, _S214);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S218;
    (&_S218)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S218)->differential_0 = _S215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S219;
    (&_S219)->primal_0 = _S191;
    (&_S219)->differential_0 = _S205;
    s_bwd_prop_mul_2(&_S218, &_S219, _S213.differential_0);
    float3  _S220 = - _S219.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S217.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S219.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S200;
    Matrix<float, 3, 3>  _S221 = _S216.differential_0 + _S218.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S221;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S220;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S222, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S223, DiffPair_float_0 * _S224, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S225, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S226, float _S227)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S222, _S223, _S224, _S225, _S226, _S227);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_1, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_0, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S228 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S228;
    Matrix<float, 3, 3>  _S229 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S229;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_1;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S228;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S228;
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
    float3  _S230 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S231 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S230);
    float3  _S232 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S233 = s_primal_ctx_dot_0(_S232, _S232);
    float _S234 = dpdepth_0 / (_S233 * _S233);
    float _S235 = - s_primal_ctx_dot_0(_S231, _S232) * - _S234;
    float _S236 = _S233 * _S234;
    float3  _S237 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S238;
    (&_S238)->primal_0 = _S232;
    (&_S238)->differential_0 = _S237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S239;
    (&_S239)->primal_0 = _S232;
    (&_S239)->differential_0 = _S237;
    s_bwd_prop_dot_0(&_S238, &_S239, _S235);
    float _S240 = - _S236;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S241;
    (&_S241)->primal_0 = _S231;
    (&_S241)->differential_0 = _S237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S232;
    (&_S242)->differential_0 = _S237;
    s_bwd_prop_dot_0(&_S241, &_S242, _S240);
    float3  _S243 = _S239.differential_0 + _S238.differential_0 + _S242.differential_0;
    Matrix<float, 3, 3>  _S244 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S245;
    (&_S245)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S245)->differential_0 = _S244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S246;
    (&_S246)->primal_0 = (*dpray_d_1).primal_0;
    (&_S246)->differential_0 = _S237;
    s_bwd_prop_mul_2(&_S245, &_S246, _S243);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S247;
    (&_S247)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S247)->differential_0 = _S244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S248;
    (&_S248)->primal_0 = _S230;
    (&_S248)->differential_0 = _S237;
    s_bwd_prop_mul_2(&_S247, &_S248, _S241.differential_0);
    float3  _S249 = - _S248.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S246.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S248.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S250 = _S245.differential_0 + _S247.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S250;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S249;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S251, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S252, DiffPair_float_0 * _S253, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S254, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S255, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S256, float3  _S257, float _S258)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S251, _S252, _S253, _S254, _S255, _S256, _S257, _S258);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_3, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_1, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S259 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S259;
    Matrix<float, 3, 3>  _S260 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S260;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_3;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S259;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S259;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S259;
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
    float _S261 = scale_3.x;
    float sx_0 = (F32_exp((_S261)));
    float _S262 = scale_3.y;
    float sy_0 = (F32_exp((_S262)));
    float sz_0 = scale_3.z - 0.5f * (_S261 + _S262);
    float x_10 = quat_4.y;
    float x2_4 = x_10 * x_10;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_4 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  _S263 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    *vert0_0 = mul_3(_S263, make_float3 (sx_0, 0.0f, 0.0f)) + mean_4;
    *vert1_0 = mul_3(_S263, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_4;
    *vert2_0 = mul_3(_S263, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_4;
    return;
}

