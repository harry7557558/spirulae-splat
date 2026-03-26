#pragma once

#include "slang.cuh"

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_0)
{
    DiffPair_float_0 _S1 = *dpx_0;
    float _S2;
    if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
        {
            _S2 = 0.0f;
        }
        else
        {
            _S2 = 0.5f * dOut_0;
        }
    }
    dpx_0->primal_0 = _S1.primal_0;
    dpx_0->differential_0 = _S2;
    DiffPair_float_0 _S3 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S1.primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_0).primal_0))
        {
            _S2 = 0.0f;
        }
        else
        {
            _S2 = 0.5f * dOut_0;
        }
    }
    dpy_0->primal_0 = _S3.primal_0;
    dpy_0->differential_0 = _S2;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_1, float dOut_1)
{
    float _S4 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_1).primal_0)))))) * dOut_1;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S4;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

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

inline __device__ float dot_0(float3  x_0, float3  y_0)
{
    int i_0 = int(0);
    float result_0 = 0.0f;
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_1 = result_0 + _slang_vector_get_element(x_0, i_0) * _slang_vector_get_element(y_0, i_0);
        i_0 = i_0 + int(1);
        result_0 = result_1;
    }
    return result_0;
}

inline __device__ float dot_1(float2  x_1, float2  y_1)
{
    int i_1 = int(0);
    float result_2 = 0.0f;
    for(;;)
    {
        if(i_1 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_3 = result_2 + _slang_vector_get_element(x_1, i_1) * _slang_vector_get_element(y_1, i_1);
        i_1 = i_1 + int(1);
        result_2 = result_3;
    }
    return result_2;
}

inline __device__ float length_0(float2  x_2)
{
    return (F32_sqrt((dot_1(x_2, x_2))));
}

inline __device__ float length_1(float3  x_3)
{
    return (F32_sqrt((dot_0(x_3, x_3))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_2, DiffPair_float_0 * dpx_3, float dOut_3)
{
    DiffPair_float_0 _S5 = *dpx_3;
    float _S6 = - (*dpy_2).primal_0 / ((*dpx_3).primal_0 * (*dpx_3).primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_3;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S6;
    float _S7 = _S5.primal_0 / (_S5.primal_0 * _S5.primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_3;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = _S7;
    return;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_4)
{
    Matrix<float, 2, 2>  result_4;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(2))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_4)->rows + (r_0)), c_0) = _slang_vector_get_element(x_4.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_4;
}

inline __device__ Matrix<float, 3, 3>  transpose_1(Matrix<float, 3, 3>  x_5)
{
    Matrix<float, 3, 3>  result_5;
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
            if(c_1 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_5)->rows + (r_1)), c_1) = _slang_vector_get_element(x_5.rows[c_1], r_1);
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_5;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_4, DiffPair_float_0 * dpy_3, float dOut_4)
{
    DiffPair_float_0 _S8 = *dpx_4;
    float _S9;
    if(((*dpx_4).primal_0) < ((*dpy_3).primal_0))
    {
        _S9 = dOut_4;
    }
    else
    {
        if(((*dpx_4).primal_0) > ((*dpy_3).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_4;
        }
    }
    dpx_4->primal_0 = _S8.primal_0;
    dpx_4->differential_0 = _S9;
    DiffPair_float_0 _S10 = *dpy_3;
    if(((*dpy_3).primal_0) < (_S8.primal_0))
    {
        _S9 = dOut_4;
    }
    else
    {
        if(((*dpy_3).primal_0) > ((*dpx_4).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_4;
        }
    }
    dpy_3->primal_0 = _S10.primal_0;
    dpy_3->differential_0 = _S9;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_5)
{
    float _S11 = (*left_0).primal_0.rows[int(0)].x * dOut_5.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_5.x;
    float sum_0 = _S11 + (*left_0).primal_0.rows[int(1)].x * dOut_5.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_5.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_5.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_5.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S12 = (*left_0).primal_0.rows[int(0)].y * dOut_5.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_5.x;
    float sum_2 = _S12 + (*left_0).primal_0.rows[int(1)].y * dOut_5.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_5.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_5.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_5.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S13 = (*left_0).primal_0.rows[int(0)].z * dOut_5.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_5.x;
    float sum_4 = _S13 + (*left_0).primal_0.rows[int(1)].z * dOut_5.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_5.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_5.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_5.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_6;
    int i_2 = int(0);
    for(;;)
    {
        if(i_2 < int(3))
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
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_2], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_6, i_2) = sum_6;
        i_2 = i_2 + int(1);
    }
    return result_6;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_6)
{
    float _S14 = (F32_exp(((*dpx_5).primal_0))) * dOut_6;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S14;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_6, float dOut_7)
{
    float _S15 = _slang_select(((*dpx_6).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_6).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_7;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S15;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_7, float dOut_8)
{
    float _S16 = (F32_exp2(((*dpx_7).primal_0))) * 0.69314718246459961f * dOut_8;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S16;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, float3  dOut_9)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_8).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_4).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_9.x);
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = left_dp_0.differential_0;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_8).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_4).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_9.y);
    *&((&left_d_result_1)->y) = left_dp_1.differential_0;
    *&((&right_d_result_1)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_8).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_4).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_9.z);
    *&((&left_d_result_1)->z) = left_dp_2.differential_0;
    *&((&right_d_result_1)->z) = right_dp_2.differential_0;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = left_d_result_1;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  max_0(float3  x_6, float3  y_2)
{
    float3  result_7;
    int i_3 = int(0);
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_7, i_3) = (F32_max((_slang_vector_get_element(x_6, i_3)), (_slang_vector_get_element(y_2, i_3))));
        i_3 = i_3 + int(1);
    }
    return result_7;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_10)
{
    float _S17 = dOut_10.y;
    float _S18 = dOut_10.z;
    float _S19 = dOut_10.x;
    float _S20 = (*a_0).primal_0.z * _S17 + - (*a_0).primal_0.y * _S18;
    float _S21 = - (*a_0).primal_0.z * _S19 + (*a_0).primal_0.x * _S18;
    float _S22 = (*a_0).primal_0.y * _S19 + - (*a_0).primal_0.x * _S17;
    float3  _S23 = make_float3 (- (*b_0).primal_0.z * _S17 + (*b_0).primal_0.y * _S18, (*b_0).primal_0.z * _S19 + - (*b_0).primal_0.x * _S18, - (*b_0).primal_0.y * _S19 + (*b_0).primal_0.x * _S17);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S23;
    float3  _S24 = make_float3 (_S20, _S21, _S22);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S24;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S25 = left_2.y;
    float _S26 = right_2.z;
    float _S27 = left_2.z;
    float _S28 = right_2.y;
    float _S29 = right_2.x;
    float _S30 = left_2.x;
    return make_float3 (_S25 * _S26 - _S27 * _S28, _S27 * _S29 - _S30 * _S26, _S30 * _S28 - _S25 * _S29);
}

inline __device__ float3  normalize_0(float3  x_7)
{
    return x_7 / make_float3 (length_1(x_7));
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_0, float4  quat_0, float3  scale_0, float2  hardness_0, FixedArray<float3 , 16>  sh_coeffs_0, FixedArray<float3 , 2>  ch_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, float4  * aabb_xyxy_0, float * depth_0, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_0, mean_0) + t_0;
        float _S31 = mean_c_0.z;
        bool _S32;
        if(_S31 < near_plane_0)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S31 > far_plane_0;
        }
        if(_S32)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float _S33 = scale_0.x;
        float sx_0 = (F32_exp((_S33)));
        float _S34 = scale_0.y;
        float sy_0 = (F32_exp((_S34)));
        float sz_0 = scale_0.z - 0.5f * (_S33 + _S34);
        float x_8 = quat_0.y;
        float x2_0 = x_8 * x_8;
        float y2_0 = quat_0.z * quat_0.z;
        float z2_0 = quat_0.w * quat_0.w;
        float xy_0 = quat_0.y * quat_0.z;
        float xz_0 = quat_0.y * quat_0.w;
        float yz_0 = quat_0.z * quat_0.w;
        float wx_0 = quat_0.x * quat_0.y;
        float wy_0 = quat_0.x * quat_0.z;
        float wz_0 = quat_0.x * quat_0.w;
        Matrix<float, 3, 3>  _S35 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0)));
        float3  vert0_0 = mul_0(_S35, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
        float3  vert1_0 = mul_0(_S35, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
        float3  vert2_0 = mul_0(_S35, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
        float3  vert0_c_0 = mul_0(R_0, vert0_0) + t_0;
        float3  vert1_c_0 = mul_0(R_0, vert1_0) + t_0;
        float3  vert2_c_0 = mul_0(R_0, vert2_0) + t_0;
        float _S36 = vert0_c_0.z;
        float _S37 = vert1_c_0.z;
        float _S38 = vert2_c_0.z;
        if(_S36 < near_plane_0)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S36 > far_plane_0;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S37 < near_plane_0;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S37 > far_plane_0;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S38 < near_plane_0;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S38 > far_plane_0;
        }
        if(_S32)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float2  uv0_0;
        for(;;)
        {
            float2  uv0_1 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S36);
            if(_S36 < 0.0f)
            {
                _S32 = true;
            }
            else
            {
                float u_0 = uv0_1.x;
                float v_0 = uv0_1.y;
                float _S39 = 0.0f * v_0;
                float r2_0 = u_0 * u_0 + v_0 * v_0;
                float s_diff_r2_0 = u_0 + u_0 + (_S39 + _S39);
                float _S40 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                float _S41 = dist_coeffs_0[int(1)] + r2_0 * _S40;
                float _S42 = dist_coeffs_0[int(0)] + r2_0 * _S41;
                float radial_0 = 1.0f + r2_0 * _S42;
                float _S43 = 2.0f * dist_coeffs_0[int(4)];
                float _S44 = _S43 * u_0;
                float _S45 = 2.0f * u_0;
                float _S46 = 2.0f * dist_coeffs_0[int(5)];
                float _S47 = _S46 * u_0;
                float _S48 = 2.0f * v_0;
                float2  _S49 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S42 + (s_diff_r2_0 * _S41 + (s_diff_r2_0 * _S40 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv0_1 + make_float2 (_S43 * v_0 + 0.0f * _S44 + (s_diff_r2_0 + (_S45 + _S45)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S46 * v_0 + 0.0f * _S47 + (s_diff_r2_0 + (_S39 + 0.0f * _S48)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
                float _S50 = 0.0f * u_0;
                float s_diff_r2_1 = _S50 + _S50 + (v_0 + v_0);
                float2  _S51 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S42 + (s_diff_r2_1 * _S41 + (s_diff_r2_1 * _S40 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv0_1 + make_float2 (0.0f * _S43 * v_0 + _S44 + (s_diff_r2_1 + (_S50 + 0.0f * _S45)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S46 * v_0 + _S47 + (s_diff_r2_1 + (_S48 + _S48)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S52 = transpose_0(makeMatrix<float, 2, 2> (_S49 + make_float2 (_S49.x * dist_coeffs_0[int(8)] + _S49.y * dist_coeffs_0[int(9)], 0.0f), _S51 + make_float2 (_S51.x * dist_coeffs_0[int(8)] + _S51.y * dist_coeffs_0[int(9)], 0.0f)));
                _S32 = !((F32_min((determinant_0(_S52)), ((F32_min((_S52.rows[int(0)].x), (_S52.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
            {
                uv0_0 = uv0_1;
                break;
            }
            float u_1 = uv0_1.x;
            float v_1 = uv0_1.y;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float2  _S53 = uv0_1 * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
            float2  _S54 = _S53 + make_float2 (dist_coeffs_0[int(8)] * _S53.x + dist_coeffs_0[int(9)] * _S53.y, 0.0f);
            uv0_0 = make_float2 (fx_0 * _S54.x + cx_0, fy_0 * _S54.y + cy_0);
            break;
        }
        float2  uv1_0;
        bool all_valid_0 = true & (!_S32);
        for(;;)
        {
            float2  uv1_1 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S37);
            if(_S37 < 0.0f)
            {
                _S32 = true;
            }
            else
            {
                float u_2 = uv1_1.x;
                float v_2 = uv1_1.y;
                float _S55 = 0.0f * v_2;
                float r2_2 = u_2 * u_2 + v_2 * v_2;
                float s_diff_r2_2 = u_2 + u_2 + (_S55 + _S55);
                float _S56 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
                float _S57 = dist_coeffs_0[int(1)] + r2_2 * _S56;
                float _S58 = dist_coeffs_0[int(0)] + r2_2 * _S57;
                float radial_1 = 1.0f + r2_2 * _S58;
                float _S59 = 2.0f * dist_coeffs_0[int(4)];
                float _S60 = _S59 * u_2;
                float _S61 = 2.0f * u_2;
                float _S62 = 2.0f * dist_coeffs_0[int(5)];
                float _S63 = _S62 * u_2;
                float _S64 = 2.0f * v_2;
                float2  _S65 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S58 + (s_diff_r2_2 * _S57 + (s_diff_r2_2 * _S56 + s_diff_r2_2 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv1_1 + make_float2 (_S59 * v_2 + 0.0f * _S60 + (s_diff_r2_2 + (_S61 + _S61)) * dist_coeffs_0[int(5)] + s_diff_r2_2 * dist_coeffs_0[int(6)], _S62 * v_2 + 0.0f * _S63 + (s_diff_r2_2 + (_S55 + 0.0f * _S64)) * dist_coeffs_0[int(4)] + s_diff_r2_2 * dist_coeffs_0[int(7)]);
                float _S66 = 0.0f * u_2;
                float s_diff_r2_3 = _S66 + _S66 + (v_2 + v_2);
                float2  _S67 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S58 + (s_diff_r2_3 * _S57 + (s_diff_r2_3 * _S56 + s_diff_r2_3 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv1_1 + make_float2 (0.0f * _S59 * v_2 + _S60 + (s_diff_r2_3 + (_S66 + 0.0f * _S61)) * dist_coeffs_0[int(5)] + s_diff_r2_3 * dist_coeffs_0[int(6)], 0.0f * _S62 * v_2 + _S63 + (s_diff_r2_3 + (_S64 + _S64)) * dist_coeffs_0[int(4)] + s_diff_r2_3 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S68 = transpose_0(makeMatrix<float, 2, 2> (_S65 + make_float2 (_S65.x * dist_coeffs_0[int(8)] + _S65.y * dist_coeffs_0[int(9)], 0.0f), _S67 + make_float2 (_S67.x * dist_coeffs_0[int(8)] + _S67.y * dist_coeffs_0[int(9)], 0.0f)));
                _S32 = !((F32_min((determinant_0(_S68)), ((F32_min((_S68.rows[int(0)].x), (_S68.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
            {
                uv1_0 = uv1_1;
                break;
            }
            float u_3 = uv1_1.x;
            float v_3 = uv1_1.y;
            float r2_3 = u_3 * u_3 + v_3 * v_3;
            float2  _S69 = uv1_1 * make_float2 (1.0f + r2_3 * (dist_coeffs_0[int(0)] + r2_3 * (dist_coeffs_0[int(1)] + r2_3 * (dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_3 * v_3 + dist_coeffs_0[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + dist_coeffs_0[int(6)] * r2_3, 2.0f * dist_coeffs_0[int(5)] * u_3 * v_3 + dist_coeffs_0[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + dist_coeffs_0[int(7)] * r2_3);
            float2  _S70 = _S69 + make_float2 (dist_coeffs_0[int(8)] * _S69.x + dist_coeffs_0[int(9)] * _S69.y, 0.0f);
            uv1_0 = make_float2 (fx_0 * _S70.x + cx_0, fy_0 * _S70.y + cy_0);
            break;
        }
        float2  uv2_0;
        bool all_valid_1 = all_valid_0 & (!_S32);
        for(;;)
        {
            float2  uv2_1 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S38);
            if(_S38 < 0.0f)
            {
                _S32 = true;
            }
            else
            {
                float u_4 = uv2_1.x;
                float v_4 = uv2_1.y;
                float _S71 = 0.0f * v_4;
                float r2_4 = u_4 * u_4 + v_4 * v_4;
                float s_diff_r2_4 = u_4 + u_4 + (_S71 + _S71);
                float _S72 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
                float _S73 = dist_coeffs_0[int(1)] + r2_4 * _S72;
                float _S74 = dist_coeffs_0[int(0)] + r2_4 * _S73;
                float radial_2 = 1.0f + r2_4 * _S74;
                float _S75 = 2.0f * dist_coeffs_0[int(4)];
                float _S76 = _S75 * u_4;
                float _S77 = 2.0f * u_4;
                float _S78 = 2.0f * dist_coeffs_0[int(5)];
                float _S79 = _S78 * u_4;
                float _S80 = 2.0f * v_4;
                float2  _S81 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S74 + (s_diff_r2_4 * _S73 + (s_diff_r2_4 * _S72 + s_diff_r2_4 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv2_1 + make_float2 (_S75 * v_4 + 0.0f * _S76 + (s_diff_r2_4 + (_S77 + _S77)) * dist_coeffs_0[int(5)] + s_diff_r2_4 * dist_coeffs_0[int(6)], _S78 * v_4 + 0.0f * _S79 + (s_diff_r2_4 + (_S71 + 0.0f * _S80)) * dist_coeffs_0[int(4)] + s_diff_r2_4 * dist_coeffs_0[int(7)]);
                float _S82 = 0.0f * u_4;
                float s_diff_r2_5 = _S82 + _S82 + (v_4 + v_4);
                float2  _S83 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S74 + (s_diff_r2_5 * _S73 + (s_diff_r2_5 * _S72 + s_diff_r2_5 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv2_1 + make_float2 (0.0f * _S75 * v_4 + _S76 + (s_diff_r2_5 + (_S82 + 0.0f * _S77)) * dist_coeffs_0[int(5)] + s_diff_r2_5 * dist_coeffs_0[int(6)], 0.0f * _S78 * v_4 + _S79 + (s_diff_r2_5 + (_S80 + _S80)) * dist_coeffs_0[int(4)] + s_diff_r2_5 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S84 = transpose_0(makeMatrix<float, 2, 2> (_S81 + make_float2 (_S81.x * dist_coeffs_0[int(8)] + _S81.y * dist_coeffs_0[int(9)], 0.0f), _S83 + make_float2 (_S83.x * dist_coeffs_0[int(8)] + _S83.y * dist_coeffs_0[int(9)], 0.0f)));
                _S32 = !((F32_min((determinant_0(_S84)), ((F32_min((_S84.rows[int(0)].x), (_S84.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
            {
                uv2_0 = uv2_1;
                break;
            }
            float u_5 = uv2_1.x;
            float v_5 = uv2_1.y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S85 = uv2_1 * make_float2 (1.0f + r2_5 * (dist_coeffs_0[int(0)] + r2_5 * (dist_coeffs_0[int(1)] + r2_5 * (dist_coeffs_0[int(2)] + r2_5 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_5 * v_5 + dist_coeffs_0[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + dist_coeffs_0[int(6)] * r2_5, 2.0f * dist_coeffs_0[int(5)] * u_5 * v_5 + dist_coeffs_0[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + dist_coeffs_0[int(7)] * r2_5);
            float2  _S86 = _S85 + make_float2 (dist_coeffs_0[int(8)] * _S85.x + dist_coeffs_0[int(9)] * _S85.y, 0.0f);
            uv2_0 = make_float2 (fx_0 * _S86.x + cx_0, fy_0 * _S86.y + cy_0);
            break;
        }
        if(!(all_valid_1 & (!_S32)))
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float2  e0_0 = uv1_0 - uv0_0;
        float2  e1_0 = uv2_0 - uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(uv0_0 - uv2_0)));
        float _S87 = uv0_0.x;
        float _S88 = uv1_0.x;
        float _S89 = uv2_0.x;
        float xmax_0 = (F32_max(((F32_max((_S87), (_S88)))), (_S89))) + offset_0;
        float xmin_0 = (F32_min(((F32_min((_S87), (_S88)))), (_S89))) - offset_0;
        float _S90 = uv0_0.y;
        float _S91 = uv1_0.y;
        float _S92 = uv2_0.y;
        float ymax_0 = (F32_max(((F32_max((_S90), (_S91)))), (_S92))) + offset_0;
        float ymin_0 = (F32_min(((F32_min((_S90), (_S91)))), (_S92))) - offset_0;
        if(xmax_0 <= 0.0f)
        {
            _S32 = true;
        }
        else
        {
            _S32 = xmin_0 >= float(image_width_0);
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = ymax_0 <= 0.0f;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = ymin_0 >= float(image_height_0);
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            if(_S31 <= 0.0f)
            {
                if(xmin_0 <= 0.0f)
                {
                    _S32 = xmax_0 >= float(image_width_0);
                }
                else
                {
                    _S32 = false;
                }
                if(_S32)
                {
                    _S32 = true;
                }
                else
                {
                    if(ymin_0 <= 0.0f)
                    {
                        _S32 = ymax_0 >= float(image_width_0);
                    }
                    else
                    {
                        _S32 = false;
                    }
                }
            }
            else
            {
                _S32 = false;
            }
        }
        if(_S32)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (xmin_0, ymin_0, xmax_0, ymax_0);
        float3  _S93 = (vert0_c_0 + vert1_c_0 + vert2_c_0) / make_float3 (3.0f);
        *depth_0 = _S93.z;
        float3  _S94 = mean_0 - - mul_0(transpose_1(R_0), t_0);
        float _S95 = _S94.x;
        float _S96 = _S94.y;
        float _S97 = _S94.z;
        float norm_0 = (F32_sqrt((_S95 * _S95 + _S96 * _S96 + _S97 * _S97)));
        float x_9 = _S95 / norm_0;
        float y_3 = _S96 / norm_0;
        float z_0 = _S97 / norm_0;
        float z2_1 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_9 * x_9 - y_3 * y_3;
        float fS1_0 = 2.0f * x_9 * y_3;
        float fTmp0C_0 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        float3  color_0 = make_float3 (0.282094806432724f) * sh_coeffs_0[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * sh_coeffs_0[int(1)] + make_float3 (z_0) * sh_coeffs_0[int(2)] - make_float3 (x_9) * sh_coeffs_0[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * sh_coeffs_0[int(4)] + make_float3 (fTmp0B_0 * y_3) * sh_coeffs_0[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * sh_coeffs_0[int(6)] + make_float3 (fTmp0B_0 * x_9) * sh_coeffs_0[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * sh_coeffs_0[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_0 + y_3 * fC1_0)) * sh_coeffs_0[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * sh_coeffs_0[int(10)] + make_float3 (fTmp0C_0 * y_3) * sh_coeffs_0[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * sh_coeffs_0[int(12)] + make_float3 (fTmp0C_0 * x_9) * sh_coeffs_0[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * sh_coeffs_0[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_0 - y_3 * fS1_0)) * sh_coeffs_0[int(15)]);
        float3  _S98 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_0 + ch_coeffs_0[int(0)] + make_float3 (0.5f), _S98);
        float3  _S99 = color_0 - ch_coeffs_0[int(0)] * make_float3 (0.5f);
        float3  _S100 = ch_coeffs_0[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S99 + _S100 + make_float3 (0.5f), _S98);
        (*rgbs_0)[int(2)] = max_0(_S99 - _S100 + make_float3 (0.5f), _S98);
        (*verts_0)[int(0)] = vert0_0;
        (*verts_0)[int(1)] = vert1_0;
        (*verts_0)[int(2)] = vert2_0;
        float3  _S101 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S101 * make_float3 (float(- (F32_sign((dot_0(_S101, mean_c_0))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_1, float4  quat_1, float3  scale_1, float2  hardness_1, FixedArray<float3 , 16>  sh_coeffs_1, FixedArray<float3 , 2>  ch_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, float4  * aabb_xyxy_1, float * depth_1, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_1)
{
    float2  _S102;
    float _S103;
    float _S104;
    float _S105;
    float _S106;
    float _S107;
    float _S108;
    float _S109;
    float _S110;
    float _S111;
    float _S112;
    float _S113;
    float _S114;
    float2  _S115;
    float _S116;
    float _S117;
    bool _S118;
    bool _S119;
    bool _S120;
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_1, mean_1) + t_1;
        float _S121 = length_1(mean_c_1);
        bool _S122;
        if(_S121 < near_plane_1)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S121 > far_plane_1;
        }
        if(_S122)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S123 = scale_1.x;
        float sx_1 = (F32_exp((_S123)));
        float _S124 = scale_1.y;
        float sy_1 = (F32_exp((_S124)));
        float sz_1 = scale_1.z - 0.5f * (_S123 + _S124);
        float x_10 = quat_1.y;
        float x2_1 = x_10 * x_10;
        float y2_1 = quat_1.z * quat_1.z;
        float z2_2 = quat_1.w * quat_1.w;
        float xy_1 = quat_1.y * quat_1.z;
        float xz_1 = quat_1.y * quat_1.w;
        float yz_1 = quat_1.z * quat_1.w;
        float wx_1 = quat_1.x * quat_1.y;
        float wy_1 = quat_1.x * quat_1.z;
        float wz_1 = quat_1.x * quat_1.w;
        Matrix<float, 3, 3>  _S125 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_2), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_2), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1)));
        float3  vert0_1 = mul_0(_S125, make_float3 (sx_1, 0.0f, 0.0f)) + mean_1;
        float3  vert1_1 = mul_0(_S125, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_1;
        float3  vert2_1 = mul_0(_S125, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_1;
        float3  vert0_c_1 = mul_0(R_1, vert0_1) + t_1;
        float3  vert1_c_1 = mul_0(R_1, vert1_1) + t_1;
        float3  vert2_c_1 = mul_0(R_1, vert2_1) + t_1;
        float _S126 = length_1(vert0_c_1);
        float _S127 = length_1(vert1_c_1);
        float _S128 = length_1(vert2_c_1);
        if(_S126 < near_plane_1)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S126 > far_plane_1;
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S127 < near_plane_1;
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S127 > far_plane_1;
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S128 < near_plane_1;
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = _S128 > far_plane_1;
        }
        if(_S122)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float2  uv0_2;
        float k_0;
        for(;;)
        {
            float2  _S129 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_2 = length_0(_S129);
            float _S130 = vert0_c_1.z;
            float theta_0 = (F32_atan2((r_2), (_S130)));
            if(theta_0 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S130;
            }
            else
            {
                k_0 = theta_0 / r_2;
            }
            float2  uv0_3 = _S129 * make_float2 (k_0);
            float2  _S131 = make_float2 (1.0f, 0.0f);
            _S102 = _S131;
            _S103 = dist_coeffs_1[int(0)];
            _S104 = dist_coeffs_1[int(1)];
            _S105 = dist_coeffs_1[int(2)];
            _S106 = dist_coeffs_1[int(3)];
            _S107 = dist_coeffs_1[int(4)];
            _S108 = dist_coeffs_1[int(5)];
            _S109 = dist_coeffs_1[int(6)];
            _S110 = dist_coeffs_1[int(7)];
            _S111 = dist_coeffs_1[int(8)];
            _S112 = dist_coeffs_1[int(9)];
            float u_6 = uv0_3.x;
            float v_6 = uv0_3.y;
            float _S132 = 0.0f * v_6;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float s_diff_r2_6 = u_6 + u_6 + (_S132 + _S132);
            float _S133 = dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)];
            float _S134 = dist_coeffs_1[int(1)] + r2_6 * _S133;
            float _S135 = dist_coeffs_1[int(0)] + r2_6 * _S134;
            float _S136 = s_diff_r2_6 * _S135 + (s_diff_r2_6 * _S134 + (s_diff_r2_6 * _S133 + s_diff_r2_6 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float radial_3 = 1.0f + r2_6 * _S135;
            float _S137 = 2.0f * dist_coeffs_1[int(4)];
            _S113 = _S137;
            float _S138 = _S137 * u_6;
            float _S139 = 2.0f * u_6;
            float s_diff_du_0 = _S137 * v_6 + 0.0f * _S138 + (s_diff_r2_6 + (_S139 + _S139)) * dist_coeffs_1[int(5)] + s_diff_r2_6 * dist_coeffs_1[int(6)];
            float _S140 = 2.0f * dist_coeffs_1[int(5)];
            _S114 = _S140;
            float _S141 = _S140 * u_6;
            float _S142 = 2.0f * v_6;
            float2  _S143 = _S131 * make_float2 (radial_3) + make_float2 (_S136) * uv0_3 + make_float2 (s_diff_du_0, _S140 * v_6 + 0.0f * _S141 + (s_diff_r2_6 + (_S132 + 0.0f * _S142)) * dist_coeffs_1[int(4)] + s_diff_r2_6 * dist_coeffs_1[int(7)]);
            float2  _S144 = _S143 + make_float2 (_S143.x * dist_coeffs_1[int(8)] + _S143.y * dist_coeffs_1[int(9)], 0.0f);
            float2  _S145 = make_float2 (0.0f, 1.0f);
            _S115 = _S145;
            float _S146 = 0.0f * u_6;
            float s_diff_r2_7 = _S146 + _S146 + (v_6 + v_6);
            float _S147 = s_diff_r2_7 * _S135 + (s_diff_r2_7 * _S134 + (s_diff_r2_7 * _S133 + s_diff_r2_7 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float _S148 = 0.0f * _S137;
            _S116 = _S148;
            float s_diff_du_1 = _S148 * v_6 + _S138 + (s_diff_r2_7 + (_S146 + 0.0f * _S139)) * dist_coeffs_1[int(5)] + s_diff_r2_7 * dist_coeffs_1[int(6)];
            float _S149 = 0.0f * _S140;
            _S117 = _S149;
            float2  _S150 = _S145 * make_float2 (radial_3) + make_float2 (_S147) * uv0_3 + make_float2 (s_diff_du_1, _S149 * v_6 + _S141 + (s_diff_r2_7 + (_S142 + _S142)) * dist_coeffs_1[int(4)] + s_diff_r2_7 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S151 = transpose_0(makeMatrix<float, 2, 2> (_S144, _S150 + make_float2 (_S150.x * dist_coeffs_1[int(8)] + _S150.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S152 = !((F32_min((determinant_0(_S151)), ((F32_min((_S151.rows[int(0)].x), (_S151.rows[int(1)].y)))))) > 0.0f);
            _S118 = _S152;
            if(_S152)
            {
                uv0_2 = uv0_3;
                break;
            }
            float2  _S153 = uv0_3 * make_float2 (radial_3) + make_float2 (_S138 * v_6 + dist_coeffs_1[int(5)] * (r2_6 + _S139 * u_6) + dist_coeffs_1[int(6)] * r2_6, _S141 * v_6 + dist_coeffs_1[int(4)] * (r2_6 + _S142 * v_6) + dist_coeffs_1[int(7)] * r2_6);
            float2  _S154 = _S153 + make_float2 (dist_coeffs_1[int(8)] * _S153.x + dist_coeffs_1[int(9)] * _S153.y, 0.0f);
            uv0_2 = make_float2 (fx_1 * _S154.x + cx_1, fy_1 * _S154.y + cy_1);
            break;
        }
        float2  uv1_2;
        bool all_valid_2 = true & (!_S118);
        for(;;)
        {
            float2  _S155 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_3 = length_0(_S155);
            float _S156 = vert1_c_1.z;
            float theta_1 = (F32_atan2((r_3), (_S156)));
            if(theta_1 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S156;
            }
            else
            {
                k_0 = theta_1 / r_3;
            }
            float2  uv1_3 = _S155 * make_float2 (k_0);
            float u_7 = uv1_3.x;
            float v_7 = uv1_3.y;
            float _S157 = 0.0f * v_7;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float s_diff_r2_8 = u_7 + u_7 + (_S157 + _S157);
            float _S158 = _S105 + r2_7 * _S106;
            float _S159 = _S104 + r2_7 * _S158;
            float _S160 = _S103 + r2_7 * _S159;
            float radial_4 = 1.0f + r2_7 * _S160;
            float _S161 = _S113 * u_7;
            float _S162 = 2.0f * u_7;
            float _S163 = _S114 * u_7;
            float _S164 = 2.0f * v_7;
            float2  _S165 = _S102 * make_float2 (radial_4) + make_float2 (s_diff_r2_8 * _S160 + (s_diff_r2_8 * _S159 + (s_diff_r2_8 * _S158 + s_diff_r2_8 * _S106 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S113 * v_7 + 0.0f * _S161 + (s_diff_r2_8 + (_S162 + _S162)) * _S108 + s_diff_r2_8 * _S109, _S114 * v_7 + 0.0f * _S163 + (s_diff_r2_8 + (_S157 + 0.0f * _S164)) * _S107 + s_diff_r2_8 * _S110);
            float _S166 = 0.0f * u_7;
            float s_diff_r2_9 = _S166 + _S166 + (v_7 + v_7);
            float2  _S167 = _S115 * make_float2 (radial_4) + make_float2 (s_diff_r2_9 * _S160 + (s_diff_r2_9 * _S159 + (s_diff_r2_9 * _S158 + s_diff_r2_9 * _S106 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S116 * v_7 + _S161 + (s_diff_r2_9 + (_S166 + 0.0f * _S162)) * _S108 + s_diff_r2_9 * _S109, _S117 * v_7 + _S163 + (s_diff_r2_9 + (_S164 + _S164)) * _S107 + s_diff_r2_9 * _S110);
            Matrix<float, 2, 2>  _S168 = transpose_0(makeMatrix<float, 2, 2> (_S165 + make_float2 (_S165.x * _S111 + _S165.y * _S112, 0.0f), _S167 + make_float2 (_S167.x * _S111 + _S167.y * _S112, 0.0f)));
            bool _S169 = !((F32_min((determinant_0(_S168)), ((F32_min((_S168.rows[int(0)].x), (_S168.rows[int(1)].y)))))) > 0.0f);
            _S119 = _S169;
            if(_S169)
            {
                uv1_2 = uv1_3;
                break;
            }
            float2  _S170 = uv1_3 * make_float2 (radial_4) + make_float2 (_S161 * v_7 + _S108 * (r2_7 + _S162 * u_7) + _S109 * r2_7, _S163 * v_7 + _S107 * (r2_7 + _S164 * v_7) + _S110 * r2_7);
            float2  _S171 = _S170 + make_float2 (_S111 * _S170.x + _S112 * _S170.y, 0.0f);
            uv1_2 = make_float2 (fx_1 * _S171.x + cx_1, fy_1 * _S171.y + cy_1);
            break;
        }
        float2  uv2_2;
        bool all_valid_3 = all_valid_2 & (!_S119);
        for(;;)
        {
            float2  _S172 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_4 = length_0(_S172);
            float _S173 = vert2_c_1.z;
            float theta_2 = (F32_atan2((r_4), (_S173)));
            if(theta_2 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S173;
            }
            else
            {
                k_0 = theta_2 / r_4;
            }
            float2  uv2_3 = _S172 * make_float2 (k_0);
            float u_8 = uv2_3.x;
            float v_8 = uv2_3.y;
            float _S174 = 0.0f * v_8;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float s_diff_r2_10 = u_8 + u_8 + (_S174 + _S174);
            float _S175 = _S105 + r2_8 * _S106;
            float _S176 = _S104 + r2_8 * _S175;
            float _S177 = _S103 + r2_8 * _S176;
            float radial_5 = 1.0f + r2_8 * _S177;
            float _S178 = _S113 * u_8;
            float _S179 = 2.0f * u_8;
            float _S180 = _S114 * u_8;
            float _S181 = 2.0f * v_8;
            float2  _S182 = _S102 * make_float2 (radial_5) + make_float2 (s_diff_r2_10 * _S177 + (s_diff_r2_10 * _S176 + (s_diff_r2_10 * _S175 + s_diff_r2_10 * _S106 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S113 * v_8 + 0.0f * _S178 + (s_diff_r2_10 + (_S179 + _S179)) * _S108 + s_diff_r2_10 * _S109, _S114 * v_8 + 0.0f * _S180 + (s_diff_r2_10 + (_S174 + 0.0f * _S181)) * _S107 + s_diff_r2_10 * _S110);
            float _S183 = 0.0f * u_8;
            float s_diff_r2_11 = _S183 + _S183 + (v_8 + v_8);
            float2  _S184 = _S115 * make_float2 (radial_5) + make_float2 (s_diff_r2_11 * _S177 + (s_diff_r2_11 * _S176 + (s_diff_r2_11 * _S175 + s_diff_r2_11 * _S106 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S116 * v_8 + _S178 + (s_diff_r2_11 + (_S183 + 0.0f * _S179)) * _S108 + s_diff_r2_11 * _S109, _S117 * v_8 + _S180 + (s_diff_r2_11 + (_S181 + _S181)) * _S107 + s_diff_r2_11 * _S110);
            Matrix<float, 2, 2>  _S185 = transpose_0(makeMatrix<float, 2, 2> (_S182 + make_float2 (_S182.x * _S111 + _S182.y * _S112, 0.0f), _S184 + make_float2 (_S184.x * _S111 + _S184.y * _S112, 0.0f)));
            bool _S186 = !((F32_min((determinant_0(_S185)), ((F32_min((_S185.rows[int(0)].x), (_S185.rows[int(1)].y)))))) > 0.0f);
            _S120 = _S186;
            if(_S186)
            {
                uv2_2 = uv2_3;
                break;
            }
            float2  _S187 = uv2_3 * make_float2 (radial_5) + make_float2 (_S178 * v_8 + _S108 * (r2_8 + _S179 * u_8) + _S109 * r2_8, _S180 * v_8 + _S107 * (r2_8 + _S181 * v_8) + _S110 * r2_8);
            float2  _S188 = _S187 + make_float2 (_S111 * _S187.x + _S112 * _S187.y, 0.0f);
            uv2_2 = make_float2 (fx_1 * _S188.x + cx_1, fy_1 * _S188.y + cy_1);
            break;
        }
        if(!(all_valid_3 & (!_S120)))
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float2  e0_1 = uv1_2 - uv0_2;
        float2  e1_1 = uv2_2 - uv1_2;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(uv0_2 - uv2_2)));
        float _S189 = uv0_2.x;
        float _S190 = uv1_2.x;
        float _S191 = uv2_2.x;
        float xmax_1 = (F32_max(((F32_max((_S189), (_S190)))), (_S191))) + offset_1;
        float xmin_1 = (F32_min(((F32_min((_S189), (_S190)))), (_S191))) - offset_1;
        float _S192 = uv0_2.y;
        float _S193 = uv1_2.y;
        float _S194 = uv2_2.y;
        float ymax_1 = (F32_max(((F32_max((_S192), (_S193)))), (_S194))) + offset_1;
        float ymin_1 = (F32_min(((F32_min((_S192), (_S193)))), (_S194))) - offset_1;
        if(xmax_1 <= 0.0f)
        {
            _S122 = true;
        }
        else
        {
            _S122 = xmin_1 >= float(image_width_1);
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = ymax_1 <= 0.0f;
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            _S122 = ymin_1 >= float(image_height_1);
        }
        if(_S122)
        {
            _S122 = true;
        }
        else
        {
            if((mean_c_1.z) <= 0.0f)
            {
                if(xmin_1 <= 0.0f)
                {
                    _S122 = xmax_1 >= float(image_width_1);
                }
                else
                {
                    _S122 = false;
                }
                if(_S122)
                {
                    _S122 = true;
                }
                else
                {
                    if(ymin_1 <= 0.0f)
                    {
                        _S122 = ymax_1 >= float(image_width_1);
                    }
                    else
                    {
                        _S122 = false;
                    }
                }
            }
            else
            {
                _S122 = false;
            }
        }
        if(_S122)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (xmin_1, ymin_1, xmax_1, ymax_1);
        float3  _S195 = (vert0_c_1 + vert1_c_1 + vert2_c_1) / make_float3 (3.0f);
        *depth_1 = _S195.z;
        float3  _S196 = mean_1 - - mul_0(transpose_1(R_1), t_1);
        float _S197 = _S196.x;
        float _S198 = _S196.y;
        float _S199 = _S196.z;
        float norm_1 = (F32_sqrt((_S197 * _S197 + _S198 * _S198 + _S199 * _S199)));
        float x_11 = _S197 / norm_1;
        float y_4 = _S198 / norm_1;
        float z_1 = _S199 / norm_1;
        float z2_3 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_11 * x_11 - y_4 * y_4;
        float fS1_1 = 2.0f * x_11 * y_4;
        float fTmp0C_1 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        float3  color_1 = make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * sh_coeffs_1[int(1)] + make_float3 (z_1) * sh_coeffs_1[int(2)] - make_float3 (x_11) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_4) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_11) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_11 * fS1_1 + y_4 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_4) * sh_coeffs_1[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_11) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_11 * fC1_1 - y_4 * fS1_1)) * sh_coeffs_1[int(15)]);
        float3  _S200 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_1 + ch_coeffs_1[int(0)] + make_float3 (0.5f), _S200);
        float3  _S201 = color_1 - ch_coeffs_1[int(0)] * make_float3 (0.5f);
        float3  _S202 = ch_coeffs_1[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S201 + _S202 + make_float3 (0.5f), _S200);
        (*rgbs_1)[int(2)] = max_0(_S201 - _S202 + make_float3 (0.5f), _S200);
        (*verts_1)[int(0)] = vert0_1;
        (*verts_1)[int(1)] = vert1_1;
        (*verts_1)[int(2)] = vert2_1;
        float3  _S203 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S203 * make_float3 (float(- (F32_sign((dot_0(_S203, mean_c_1))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, float4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
{
    float3  mean_c_2 = mul_0(R_2, mean_2) + t_2;
    float _S204 = scale_2.x;
    float sx_2 = (F32_exp((_S204)));
    float _S205 = scale_2.y;
    float sy_2 = (F32_exp((_S205)));
    float sz_2 = scale_2.z - 0.5f * (_S204 + _S205);
    float x_12 = quat_2.y;
    float x2_2 = x_12 * x_12;
    float y2_2 = quat_2.z * quat_2.z;
    float z2_4 = quat_2.w * quat_2.w;
    float xy_2 = quat_2.y * quat_2.z;
    float xz_2 = quat_2.y * quat_2.w;
    float yz_2 = quat_2.z * quat_2.w;
    float wx_2 = quat_2.x * quat_2.y;
    float wy_2 = quat_2.x * quat_2.z;
    float wz_2 = quat_2.x * quat_2.w;
    Matrix<float, 3, 3>  _S206 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_4), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_4), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)));
    float3  vert0_2 = mul_0(_S206, make_float3 (sx_2, 0.0f, 0.0f)) + mean_2;
    float3  vert1_2 = mul_0(_S206, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_2;
    float3  vert2_2 = mul_0(_S206, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_2;
    float3  vert0_c_2 = mul_0(R_2, vert0_2) + t_2;
    float3  vert1_c_2 = mul_0(R_2, vert1_2) + t_2;
    float3  vert2_c_2 = mul_0(R_2, vert2_2) + t_2;
    float2  _S207 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_9 = _S207.x;
    float v_9 = _S207.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S208 = 2.0f * dist_coeffs_2[int(4)];
    float _S209 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S210 = _S207 * make_float2 (1.0f + r2_9 * (dist_coeffs_2[int(0)] + r2_9 * (dist_coeffs_2[int(1)] + r2_9 * (dist_coeffs_2[int(2)] + r2_9 * dist_coeffs_2[int(3)])))) + make_float2 (_S208 * u_9 * v_9 + dist_coeffs_2[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_2[int(6)] * r2_9, _S209 * u_9 * v_9 + dist_coeffs_2[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_2[int(7)] * r2_9);
    float2  _S211 = _S210 + make_float2 (dist_coeffs_2[int(8)] * _S210.x + dist_coeffs_2[int(9)] * _S210.y, 0.0f);
    float _S212 = fx_2 * _S211.x + cx_2;
    float _S213 = fy_2 * _S211.y + cy_2;
    float2  uv0_4 = make_float2 (_S212, _S213);
    float2  _S214 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_10 = _S214.x;
    float v_10 = _S214.y;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float2  _S215 = _S214 * make_float2 (1.0f + r2_10 * (dist_coeffs_2[int(0)] + r2_10 * (dist_coeffs_2[int(1)] + r2_10 * (dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)])))) + make_float2 (_S208 * u_10 * v_10 + dist_coeffs_2[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + dist_coeffs_2[int(6)] * r2_10, _S209 * u_10 * v_10 + dist_coeffs_2[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + dist_coeffs_2[int(7)] * r2_10);
    float2  _S216 = _S215 + make_float2 (dist_coeffs_2[int(8)] * _S215.x + dist_coeffs_2[int(9)] * _S215.y, 0.0f);
    float _S217 = fx_2 * _S216.x + cx_2;
    float _S218 = fy_2 * _S216.y + cy_2;
    float2  uv1_4 = make_float2 (_S217, _S218);
    float2  _S219 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_11 = _S219.x;
    float v_11 = _S219.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float2  _S220 = _S219 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S208 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S209 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
    float2  _S221 = _S220 + make_float2 (dist_coeffs_2[int(8)] * _S220.x + dist_coeffs_2[int(9)] * _S220.y, 0.0f);
    float _S222 = fx_2 * _S221.x + cx_2;
    float _S223 = fy_2 * _S221.y + cy_2;
    float2  uv2_4 = make_float2 (_S222, _S223);
    float2  e0_2 = uv1_4 - uv0_4;
    float2  e1_2 = uv2_4 - uv1_4;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(uv0_4 - uv2_4)));
    *aabb_xyxy_2 = make_float4 ((F32_min(((F32_min((_S212), (_S217)))), (_S222))) - offset_2, (F32_min(((F32_min((_S213), (_S218)))), (_S223))) - offset_2, (F32_max(((F32_max((_S212), (_S217)))), (_S222))) + offset_2, (F32_max(((F32_max((_S213), (_S218)))), (_S223))) + offset_2);
    *depth_2 = ((vert0_c_2 + vert1_c_2 + vert2_c_2) / make_float3 (3.0f)).z;
    float3  _S224 = mean_2 - - mul_0(transpose_1(R_2), t_2);
    float _S225 = _S224.x;
    float _S226 = _S224.y;
    float _S227 = _S224.z;
    float norm_2 = (F32_sqrt((_S225 * _S225 + _S226 * _S226 + _S227 * _S227)));
    float x_13 = _S225 / norm_2;
    float y_5 = _S226 / norm_2;
    float z_2 = _S227 / norm_2;
    float z2_5 = z_2 * z_2;
    float fTmp0B_2 = -1.09254848957061768f * z_2;
    float fC1_2 = x_13 * x_13 - y_5 * y_5;
    float fS1_2 = 2.0f * x_13 * y_5;
    float fTmp0C_2 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_2;
    float3  color_2 = make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_2[int(1)] + make_float3 (z_2) * sh_coeffs_2[int(2)] - make_float3 (x_13) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_5) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_13) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_13 * fS1_2 + y_5 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_5) * sh_coeffs_2[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_13) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_13 * fC1_2 - y_5 * fS1_2)) * sh_coeffs_2[int(15)]);
    float3  _S228 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_2 + ch_coeffs_2[int(0)] + make_float3 (0.5f), _S228);
    float3  _S229 = color_2 - ch_coeffs_2[int(0)] * make_float3 (0.5f);
    float3  _S230 = ch_coeffs_2[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S229 + _S230 + make_float3 (0.5f), _S228);
    (*rgbs_2)[int(2)] = max_0(_S229 - _S230 + make_float3 (0.5f), _S228);
    (*verts_2)[int(0)] = vert0_2;
    (*verts_2)[int(1)] = vert1_2;
    (*verts_2)[int(2)] = vert2_2;
    float3  _S231 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S231 * make_float3 (float(- (F32_sign((dot_0(_S231, mean_c_2))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, float4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
{
    float3  mean_c_3 = mul_0(R_3, mean_3) + t_3;
    float _S232 = scale_3.x;
    float sx_3 = (F32_exp((_S232)));
    float _S233 = scale_3.y;
    float sy_3 = (F32_exp((_S233)));
    float sz_3 = scale_3.z - 0.5f * (_S232 + _S233);
    float x_14 = quat_3.y;
    float x2_3 = x_14 * x_14;
    float y2_3 = quat_3.z * quat_3.z;
    float z2_6 = quat_3.w * quat_3.w;
    float xy_3 = quat_3.y * quat_3.z;
    float xz_3 = quat_3.y * quat_3.w;
    float yz_3 = quat_3.z * quat_3.w;
    float wx_3 = quat_3.x * quat_3.y;
    float wy_3 = quat_3.x * quat_3.z;
    float wz_3 = quat_3.x * quat_3.w;
    Matrix<float, 3, 3>  _S234 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_6), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_6), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    float3  vert0_3 = mul_0(_S234, make_float3 (sx_3, 0.0f, 0.0f)) + mean_3;
    float3  vert1_3 = mul_0(_S234, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_3;
    float3  vert2_3 = mul_0(_S234, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_3;
    float3  vert0_c_3 = mul_0(R_3, vert0_3) + t_3;
    float3  vert1_c_3 = mul_0(R_3, vert1_3) + t_3;
    float3  vert2_c_3 = mul_0(R_3, vert2_3) + t_3;
    float2  _S235 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_5 = length_0(_S235);
    float _S236 = vert0_c_3.z;
    float theta_3 = (F32_atan2((r_5), (_S236)));
    float k_1;
    if(theta_3 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_3 * theta_3 / 3.0f) / _S236;
    }
    else
    {
        k_1 = theta_3 / r_5;
    }
    float2  _S237 = _S235 * make_float2 (k_1);
    float u_12 = _S237.x;
    float v_12 = _S237.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S238 = 2.0f * dist_coeffs_3[int(4)];
    float _S239 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S240 = _S237 * make_float2 (1.0f + r2_12 * (dist_coeffs_3[int(0)] + r2_12 * (dist_coeffs_3[int(1)] + r2_12 * (dist_coeffs_3[int(2)] + r2_12 * dist_coeffs_3[int(3)])))) + make_float2 (_S238 * u_12 * v_12 + dist_coeffs_3[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + dist_coeffs_3[int(6)] * r2_12, _S239 * u_12 * v_12 + dist_coeffs_3[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + dist_coeffs_3[int(7)] * r2_12);
    float2  _S241 = _S240 + make_float2 (dist_coeffs_3[int(8)] * _S240.x + dist_coeffs_3[int(9)] * _S240.y, 0.0f);
    float _S242 = fx_3 * _S241.x + cx_3;
    float _S243 = fy_3 * _S241.y + cy_3;
    float2  uv0_5 = make_float2 (_S242, _S243);
    float2  _S244 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_6 = length_0(_S244);
    float _S245 = vert1_c_3.z;
    float theta_4 = (F32_atan2((r_6), (_S245)));
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S245;
    }
    else
    {
        k_1 = theta_4 / r_6;
    }
    float2  _S246 = _S244 * make_float2 (k_1);
    float u_13 = _S246.x;
    float v_13 = _S246.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S247 = _S246 * make_float2 (1.0f + r2_13 * (dist_coeffs_3[int(0)] + r2_13 * (dist_coeffs_3[int(1)] + r2_13 * (dist_coeffs_3[int(2)] + r2_13 * dist_coeffs_3[int(3)])))) + make_float2 (_S238 * u_13 * v_13 + dist_coeffs_3[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_3[int(6)] * r2_13, _S239 * u_13 * v_13 + dist_coeffs_3[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_3[int(7)] * r2_13);
    float2  _S248 = _S247 + make_float2 (dist_coeffs_3[int(8)] * _S247.x + dist_coeffs_3[int(9)] * _S247.y, 0.0f);
    float _S249 = fx_3 * _S248.x + cx_3;
    float _S250 = fy_3 * _S248.y + cy_3;
    float2  uv1_5 = make_float2 (_S249, _S250);
    float2  _S251 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_7 = length_0(_S251);
    float _S252 = vert2_c_3.z;
    float theta_5 = (F32_atan2((r_7), (_S252)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S252;
    }
    else
    {
        k_1 = theta_5 / r_7;
    }
    float2  _S253 = _S251 * make_float2 (k_1);
    float u_14 = _S253.x;
    float v_14 = _S253.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S254 = _S253 * make_float2 (1.0f + r2_14 * (dist_coeffs_3[int(0)] + r2_14 * (dist_coeffs_3[int(1)] + r2_14 * (dist_coeffs_3[int(2)] + r2_14 * dist_coeffs_3[int(3)])))) + make_float2 (_S238 * u_14 * v_14 + dist_coeffs_3[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + dist_coeffs_3[int(6)] * r2_14, _S239 * u_14 * v_14 + dist_coeffs_3[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + dist_coeffs_3[int(7)] * r2_14);
    float2  _S255 = _S254 + make_float2 (dist_coeffs_3[int(8)] * _S254.x + dist_coeffs_3[int(9)] * _S254.y, 0.0f);
    float _S256 = fx_3 * _S255.x + cx_3;
    float _S257 = fy_3 * _S255.y + cy_3;
    float2  uv2_5 = make_float2 (_S256, _S257);
    float2  e0_3 = uv1_5 - uv0_5;
    float2  e1_3 = uv2_5 - uv1_5;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(uv0_5 - uv2_5)));
    *aabb_xyxy_3 = make_float4 ((F32_min(((F32_min((_S242), (_S249)))), (_S256))) - offset_3, (F32_min(((F32_min((_S243), (_S250)))), (_S257))) - offset_3, (F32_max(((F32_max((_S242), (_S249)))), (_S256))) + offset_3, (F32_max(((F32_max((_S243), (_S250)))), (_S257))) + offset_3);
    *depth_3 = ((vert0_c_3 + vert1_c_3 + vert2_c_3) / make_float3 (3.0f)).z;
    float3  _S258 = mean_3 - - mul_0(transpose_1(R_3), t_3);
    float _S259 = _S258.x;
    float _S260 = _S258.y;
    float _S261 = _S258.z;
    float norm_3 = (F32_sqrt((_S259 * _S259 + _S260 * _S260 + _S261 * _S261)));
    float x_15 = _S259 / norm_3;
    float y_6 = _S260 / norm_3;
    float z_3 = _S261 / norm_3;
    float z2_7 = z_3 * z_3;
    float fTmp0B_3 = -1.09254848957061768f * z_3;
    float fC1_3 = x_15 * x_15 - y_6 * y_6;
    float fS1_3 = 2.0f * x_15 * y_6;
    float fTmp0C_3 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_3;
    float3  color_3 = make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_3[int(1)] + make_float3 (z_3) * sh_coeffs_3[int(2)] - make_float3 (x_15) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_6) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_15) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_3 + y_6 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_6) * sh_coeffs_3[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_15) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_3 - y_6 * fS1_3)) * sh_coeffs_3[int(15)]);
    float3  _S262 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_3 + ch_coeffs_3[int(0)] + make_float3 (0.5f), _S262);
    float3  _S263 = color_3 - ch_coeffs_3[int(0)] * make_float3 (0.5f);
    float3  _S264 = ch_coeffs_3[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S263 + _S264 + make_float3 (0.5f), _S262);
    (*rgbs_3)[int(2)] = max_0(_S263 - _S264 + make_float3 (0.5f), _S262);
    (*verts_3)[int(0)] = vert0_3;
    (*verts_3)[int(1)] = vert1_3;
    (*verts_3)[int(2)] = vert2_3;
    float3  _S265 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S265 * make_float3 (float(- (F32_sign((dot_0(_S265, mean_c_3))))));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S266, float3  _S267)
{
    return mul_0(_S266, _S267);
}

inline __device__ float s_primal_ctx_exp_0(float _S268)
{
    return (F32_exp((_S268)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S269)
{
    return (F32_sqrt((_S269)));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S270, float3  _S271)
{
    return cross_0(_S270, _S271);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S272, float3  _S273)
{
    return dot_0(_S272, _S273);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S274, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S275, float _S276)
{
    _d_dot_0(_S274, _S275, _S276);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S277, float _S278)
{
    _d_sqrt_0(_S277, _S278);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S279 = (*dpx_9).primal_0.x;
    float _S280 = (*dpx_9).primal_0.y;
    float _S281 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S282;
    (&_S282)->primal_0 = _S279 * _S279 + _S280 * _S280 + _S281 * _S281;
    (&_S282)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S282, _s_dOut_0);
    float _S283 = (*dpx_9).primal_0.z * _S282.differential_0;
    float _S284 = _S283 + _S283;
    float _S285 = (*dpx_9).primal_0.y * _S282.differential_0;
    float _S286 = _S285 + _S285;
    float _S287 = (*dpx_9).primal_0.x * _S282.differential_0;
    float _S288 = _S287 + _S287;
    float3  _S289 = make_float3 (0.0f);
    *&((&_S289)->z) = _S284;
    *&((&_S289)->y) = _S286;
    *&((&_S289)->x) = _S288;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S289;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S290, float _S291)
{
    s_bwd_prop_length_impl_0(_S290, _S291);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_1)
{
    float _S292 = length_1((*dpx_10).primal_0);
    float3  _S293 = (*dpx_10).primal_0 * _s_dOut_1;
    float3  _S294 = make_float3 (1.0f / _S292) * _s_dOut_1;
    float _S295 = - ((_S293.x + _S293.y + _S293.z) / (_S292 * _S292));
    float3  _S296 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S297;
    (&_S297)->primal_0 = (*dpx_10).primal_0;
    (&_S297)->differential_0 = _S296;
    s_bwd_length_impl_0(&_S297, _S295);
    float3  _S298 = _S294 + _S297.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S298;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S299, float3  _S300)
{
    s_bwd_prop_normalize_impl_0(_S299, _S300);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S301, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S302, float3  _S303)
{
    _d_cross_0(_S301, _S302, _S303);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S304, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S305, float3  _S306)
{
    _d_max_vector_0(_S304, _S305, _S306);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S307, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S308, float3  _S309)
{
    _d_mul_0(_S307, _S308, _S309);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S310, float _S311)
{
    _d_exp2_0(_S310, _S311);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float _s_dOut_2)
{
    float _S312 = (*dpx_11).primal_0.x;
    float _S313 = (*dpx_11).primal_0.y;
    DiffPair_float_0 _S314;
    (&_S314)->primal_0 = _S312 * _S312 + _S313 * _S313;
    (&_S314)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S314, _s_dOut_2);
    float _S315 = (*dpx_11).primal_0.y * _S314.differential_0;
    float _S316 = _S315 + _S315;
    float _S317 = (*dpx_11).primal_0.x * _S314.differential_0;
    float _S318 = _S317 + _S317;
    float2  _S319 = make_float2 (0.0f);
    *&((&_S319)->y) = _S316;
    *&((&_S319)->x) = _S318;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S319;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S320, float _S321)
{
    s_bwd_prop_length_impl_1(_S320, _S321);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S322, float _S323)
{
    _d_abs_0(_S322, _S323);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S324, float _S325)
{
    _d_exp_0(_S324, _S325);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S326 = scale_4.x;
    float _S327 = s_primal_ctx_exp_0(_S326);
    float _S328 = scale_4.y;
    float _S329 = s_primal_ctx_exp_0(_S328);
    float sz_4 = scale_4.z - 0.5f * (_S326 + _S328);
    float _S330 = quat_4.y;
    float x2_4 = _S330 * _S330;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_8 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  _S331 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S332 = make_float3 (_S327, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S331, _S332) + mean_4;
    float _S333 = -0.5f + sz_4;
    float3  _S334 = make_float3 (_S327 * _S333, _S329, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S331, _S334) + mean_4;
    float _S335 = -0.5f - sz_4;
    float3  _S336 = make_float3 (_S327 * _S335, - _S329, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S331, _S336) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S337 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S338 = vert0_c_4.z;
    float2  _S339 = make_float2 (_S338);
    float2  _S340 = _S337 / make_float2 (_S338);
    float2  _S341 = make_float2 (_S338 * _S338);
    float u_15 = _S340.x;
    float v_15 = _S340.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S342 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S343 = dist_coeffs_4[int(1)] + r2_15 * _S342;
    float _S344 = dist_coeffs_4[int(0)] + r2_15 * _S343;
    float radial_6 = 1.0f + r2_15 * _S344;
    float _S345 = 2.0f * dist_coeffs_4[int(4)];
    float _S346 = _S345 * u_15;
    float _S347 = 2.0f * u_15;
    float _S348 = 2.0f * dist_coeffs_4[int(5)];
    float _S349 = _S348 * u_15;
    float _S350 = 2.0f * v_15;
    float2  _S351 = _S340 * make_float2 (radial_6) + make_float2 (_S346 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S347 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S349 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S350 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S352 = _S351 + make_float2 (dist_coeffs_4[int(8)] * _S351.x + dist_coeffs_4[int(9)] * _S351.y, 0.0f);
    float _S353 = fx_4 * _S352.x + cx_4;
    float _S354 = fy_4 * _S352.y + cy_4;
    float2  uv0_6 = make_float2 (_S353, _S354);
    float2  _S355 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S356 = vert1_c_4.z;
    float2  _S357 = make_float2 (_S356);
    float2  _S358 = _S355 / make_float2 (_S356);
    float2  _S359 = make_float2 (_S356 * _S356);
    float u_16 = _S358.x;
    float v_16 = _S358.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S360 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S361 = dist_coeffs_4[int(1)] + r2_16 * _S360;
    float _S362 = dist_coeffs_4[int(0)] + r2_16 * _S361;
    float radial_7 = 1.0f + r2_16 * _S362;
    float _S363 = _S345 * u_16;
    float _S364 = 2.0f * u_16;
    float _S365 = _S348 * u_16;
    float _S366 = 2.0f * v_16;
    float2  _S367 = _S358 * make_float2 (radial_7) + make_float2 (_S363 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S364 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S365 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S366 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S368 = _S367 + make_float2 (dist_coeffs_4[int(8)] * _S367.x + dist_coeffs_4[int(9)] * _S367.y, 0.0f);
    float _S369 = fx_4 * _S368.x + cx_4;
    float _S370 = fy_4 * _S368.y + cy_4;
    float2  uv1_6 = make_float2 (_S369, _S370);
    float2  _S371 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S372 = vert2_c_4.z;
    float2  _S373 = make_float2 (_S372);
    float2  _S374 = _S371 / make_float2 (_S372);
    float2  _S375 = make_float2 (_S372 * _S372);
    float u_17 = _S374.x;
    float v_17 = _S374.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S376 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S377 = dist_coeffs_4[int(1)] + r2_17 * _S376;
    float _S378 = dist_coeffs_4[int(0)] + r2_17 * _S377;
    float radial_8 = 1.0f + r2_17 * _S378;
    float _S379 = _S345 * u_17;
    float _S380 = 2.0f * u_17;
    float _S381 = _S348 * u_17;
    float _S382 = 2.0f * v_17;
    float2  _S383 = _S374 * make_float2 (radial_8) + make_float2 (_S379 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S380 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S381 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S382 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S384 = _S383 + make_float2 (dist_coeffs_4[int(8)] * _S383.x + dist_coeffs_4[int(9)] * _S383.y, 0.0f);
    float _S385 = fx_4 * _S384.x + cx_4;
    float _S386 = fy_4 * _S384.y + cy_4;
    float2  uv2_6 = make_float2 (_S385, _S386);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S387 = e0_4.x;
    float _S388 = e1_4.y;
    float _S389 = e0_4.y;
    float _S390 = e1_4.x;
    float _S391 = _S387 * _S388 - _S389 * _S390;
    float _S392 = 1.0f - hardness_4.y;
    float _S393 = -1.0f / _S392;
    float _S394 = _S392 * _S392;
    float _S395 = (F32_max((_S353), (_S369)));
    float _S396 = (F32_min((_S353), (_S369)));
    float _S397 = (F32_max((_S354), (_S370)));
    float _S398 = (F32_min((_S354), (_S370)));
    Matrix<float, 3, 3>  _S399 = transpose_1(R_4);
    float3  _S400 = mean_4 - - s_primal_ctx_mul_0(_S399, t_4);
    float _S401 = _S400.x;
    float _S402 = _S400.y;
    float _S403 = _S400.z;
    float _S404 = _S401 * _S401 + _S402 * _S402 + _S403 * _S403;
    float _S405 = s_primal_ctx_sqrt_0(_S404);
    float x_16 = _S401 / _S405;
    float3  _S406 = make_float3 (x_16);
    float _S407 = _S405 * _S405;
    float y_7 = _S402 / _S405;
    float z_4 = _S403 / _S405;
    float3  _S408 = make_float3 (z_4);
    float _S409 = - y_7;
    float3  _S410 = make_float3 (_S409);
    float z2_9 = z_4 * z_4;
    float fTmp0B_4 = -1.09254848957061768f * z_4;
    float fC1_4 = x_16 * x_16 - y_7 * y_7;
    float _S411 = 2.0f * x_16;
    float fS1_4 = _S411 * y_7;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S412 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_16;
    float3  _S413 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_7;
    float3  _S414 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S415 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S416 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_4;
    float _S417 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_4 * _S417;
    float3  _S418 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_16;
    float3  _S419 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_7;
    float3  _S420 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S421 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S422 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_16 * fC1_4 - y_7 * fS1_4);
    float3  _S423 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_16 * fS1_4 + y_7 * fC1_4);
    float3  _S424 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S409) * sh_coeffs_4[int(1)] + make_float3 (z_4) * sh_coeffs_4[int(2)] - make_float3 (x_16) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]);
    float3  _S425 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S426 = make_float3 (0.0f);
    float3  _S427 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S428 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S429 = make_float3 (_S428);
    float3  _S430 = ch_coeffs_4[int(1)] * make_float3 (_S428);
    float3  _S431 = _S427 + _S430 + make_float3 (0.5f);
    float3  _S432 = _S427 - _S430 + make_float3 (0.5f);
    float3  _S433 = vert1_c_4 - vert0_c_4;
    float3  _S434 = vert2_c_4 - vert0_c_4;
    float3  _S435 = s_primal_ctx_cross_0(_S433, _S434);
    float3  _S436 = normalize_0(_S435);
    float3  _S437 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S436, mean_c_4)))))) * v_normal_0;
    float3  _S438 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S439;
    (&_S439)->primal_0 = _S436;
    (&_S439)->differential_0 = _S438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S440;
    (&_S440)->primal_0 = mean_c_4;
    (&_S440)->differential_0 = _S438;
    s_bwd_prop_dot_0(&_S439, &_S440, 0.0f);
    float3  _S441 = _S437 + _S439.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S442;
    (&_S442)->primal_0 = _S435;
    (&_S442)->differential_0 = _S438;
    s_bwd_normalize_impl_0(&_S442, _S441);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S443;
    (&_S443)->primal_0 = _S433;
    (&_S443)->differential_0 = _S438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S444;
    (&_S444)->primal_0 = _S434;
    (&_S444)->differential_0 = _S438;
    s_bwd_prop_cross_0(&_S443, &_S444, _S442.differential_0);
    float3  _S445 = - _S444.differential_0;
    float3  _S446 = - _S443.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S447;
    (&_S447)->primal_0 = _S432;
    (&_S447)->differential_0 = _S438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S448;
    (&_S448)->primal_0 = _S426;
    (&_S448)->differential_0 = _S438;
    s_bwd_prop_max_0(&_S447, &_S448, v_rgbs_0[int(2)]);
    float3  _S449 = - _S447.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S450;
    (&_S450)->primal_0 = _S431;
    (&_S450)->differential_0 = _S438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S451;
    (&_S451)->primal_0 = _S426;
    (&_S451)->differential_0 = _S438;
    s_bwd_prop_max_0(&_S450, &_S451, v_rgbs_0[int(1)]);
    float3  _S452 = _S429 * (_S449 + _S450.differential_0);
    float3  _S453 = _S447.differential_0 + _S450.differential_0;
    float3  _S454 = make_float3 (0.5f) * - _S453;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S455;
    (&_S455)->primal_0 = _S425;
    (&_S455)->differential_0 = _S438;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S456;
    (&_S456)->primal_0 = _S426;
    (&_S456)->differential_0 = _S438;
    s_bwd_prop_max_0(&_S455, &_S456, v_rgbs_0[int(0)]);
    float3  _S457 = _S454 + _S455.differential_0;
    float3  _S458 = _S453 + _S455.differential_0;
    float3  _S459 = _S423 * _S458;
    float3  _S460 = sh_coeffs_4[int(15)] * _S458;
    float3  _S461 = _S421 * _S458;
    float3  _S462 = sh_coeffs_4[int(14)] * _S458;
    float3  _S463 = _S419 * _S458;
    float3  _S464 = sh_coeffs_4[int(13)] * _S458;
    float3  _S465 = _S418 * _S458;
    float3  _S466 = sh_coeffs_4[int(12)] * _S458;
    float3  _S467 = _S420 * _S458;
    float3  _S468 = sh_coeffs_4[int(11)] * _S458;
    float3  _S469 = _S422 * _S458;
    float3  _S470 = sh_coeffs_4[int(10)] * _S458;
    float3  _S471 = _S424 * _S458;
    float3  _S472 = sh_coeffs_4[int(9)] * _S458;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S472.x + _S472.y + _S472.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S460.x + _S460.y + _S460.z);
    float _S473 = _S470.x + _S470.y + _S470.z;
    float _S474 = _S462.x + _S462.y + _S462.z;
    float _S475 = _S468.x + _S468.y + _S468.z;
    float _S476 = _S464.x + _S464.y + _S464.z;
    float _S477 = _S466.x + _S466.y + _S466.z;
    float _S478 = - s_diff_fC2_T_0;
    float3  _S479 = _S415 * _S458;
    float3  _S480 = sh_coeffs_4[int(8)] * _S458;
    float3  _S481 = _S413 * _S458;
    float3  _S482 = sh_coeffs_4[int(7)] * _S458;
    float3  _S483 = _S412 * _S458;
    float3  _S484 = sh_coeffs_4[int(6)] * _S458;
    float3  _S485 = _S414 * _S458;
    float3  _S486 = sh_coeffs_4[int(5)] * _S458;
    float3  _S487 = _S416 * _S458;
    float3  _S488 = sh_coeffs_4[int(4)] * _S458;
    float _S489 = _S486.x + _S486.y + _S486.z;
    float _S490 = _S482.x + _S482.y + _S482.z;
    float _S491 = fTmp1B_4 * _S473 + x_16 * s_diff_fS2_T_0 + y_7 * _S478 + 0.54627424478530884f * (_S488.x + _S488.y + _S488.z);
    float _S492 = fTmp1B_4 * _S474 + y_7 * s_diff_fS2_T_0 + x_16 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S480.x + _S480.y + _S480.z);
    float _S493 = y_7 * - _S492;
    float _S494 = x_16 * _S492;
    float _S495 = z_4 * (1.86588168144226074f * (z_4 * _S477) + -2.28522896766662598f * (y_7 * _S475 + x_16 * _S476) + 0.94617468118667603f * (_S484.x + _S484.y + _S484.z));
    float3  _S496 = make_float3 (0.48860251903533936f) * _S458;
    float3  _S497 = - _S496;
    float3  _S498 = _S406 * _S497;
    float3  _S499 = sh_coeffs_4[int(3)] * _S497;
    float3  _S500 = _S408 * _S496;
    float3  _S501 = sh_coeffs_4[int(2)] * _S496;
    float3  _S502 = _S410 * _S496;
    float3  _S503 = sh_coeffs_4[int(1)] * _S496;
    float _S504 = (_S417 * _S477 + 1.44530570507049561f * (fS1_4 * _S473 + fC1_4 * _S474) + -1.09254848957061768f * (y_7 * _S489 + x_16 * _S490) + _S495 + _S495 + _S501.x + _S501.y + _S501.z) / _S407;
    float _S505 = _S405 * _S504;
    float _S506 = (fTmp0C_4 * _S475 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S478 + fTmp0B_4 * _S489 + _S411 * _S491 + _S493 + _S493 + - (_S503.x + _S503.y + _S503.z)) / _S407;
    float _S507 = _S405 * _S506;
    float _S508 = (fTmp0C_4 * _S476 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S490 + 2.0f * (y_7 * _S491) + _S494 + _S494 + _S499.x + _S499.y + _S499.z) / _S407;
    float _S509 = _S405 * _S508;
    float _S510 = _S403 * - _S504 + _S402 * - _S506 + _S401 * - _S508;
    DiffPair_float_0 _S511;
    (&_S511)->primal_0 = _S404;
    (&_S511)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S511, _S510);
    float _S512 = _S403 * _S511.differential_0;
    float _S513 = _S402 * _S511.differential_0;
    float _S514 = _S401 * _S511.differential_0;
    float3  _S515 = make_float3 (0.282094806432724f) * _S458;
    float3  _S516 = make_float3 (_S509 + _S514 + _S514, _S507 + _S513 + _S513, _S505 + _S512 + _S512);
    float3  _S517 = - - _S516;
    Matrix<float, 3, 3>  _S518 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S519;
    (&_S519)->primal_0 = _S399;
    (&_S519)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S520;
    (&_S520)->primal_0 = t_4;
    (&_S520)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S519, &_S520, _S517);
    Matrix<float, 3, 3>  _S521 = transpose_1(_S519.differential_0);
    float3  _S522 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_0);
    DiffPair_float_0 _S523;
    (&_S523)->primal_0 = _S398;
    (&_S523)->differential_0 = 0.0f;
    DiffPair_float_0 _S524;
    (&_S524)->primal_0 = _S386;
    (&_S524)->differential_0 = 0.0f;
    _d_min_0(&_S523, &_S524, 0.0f);
    DiffPair_float_0 _S525;
    (&_S525)->primal_0 = _S354;
    (&_S525)->differential_0 = 0.0f;
    DiffPair_float_0 _S526;
    (&_S526)->primal_0 = _S370;
    (&_S526)->differential_0 = 0.0f;
    _d_min_0(&_S525, &_S526, _S523.differential_0);
    DiffPair_float_0 _S527;
    (&_S527)->primal_0 = _S397;
    (&_S527)->differential_0 = 0.0f;
    DiffPair_float_0 _S528;
    (&_S528)->primal_0 = _S386;
    (&_S528)->differential_0 = 0.0f;
    _d_max_0(&_S527, &_S528, 0.0f);
    DiffPair_float_0 _S529;
    (&_S529)->primal_0 = _S354;
    (&_S529)->differential_0 = 0.0f;
    DiffPair_float_0 _S530;
    (&_S530)->primal_0 = _S370;
    (&_S530)->differential_0 = 0.0f;
    _d_max_0(&_S529, &_S530, _S527.differential_0);
    DiffPair_float_0 _S531;
    (&_S531)->primal_0 = _S396;
    (&_S531)->differential_0 = 0.0f;
    DiffPair_float_0 _S532;
    (&_S532)->primal_0 = _S385;
    (&_S532)->differential_0 = 0.0f;
    _d_min_0(&_S531, &_S532, 0.0f);
    DiffPair_float_0 _S533;
    (&_S533)->primal_0 = _S353;
    (&_S533)->differential_0 = 0.0f;
    DiffPair_float_0 _S534;
    (&_S534)->primal_0 = _S369;
    (&_S534)->differential_0 = 0.0f;
    _d_min_0(&_S533, &_S534, _S531.differential_0);
    DiffPair_float_0 _S535;
    (&_S535)->primal_0 = _S395;
    (&_S535)->differential_0 = 0.0f;
    DiffPair_float_0 _S536;
    (&_S536)->primal_0 = _S385;
    (&_S536)->differential_0 = 0.0f;
    _d_max_0(&_S535, &_S536, 0.0f);
    DiffPair_float_0 _S537;
    (&_S537)->primal_0 = _S353;
    (&_S537)->differential_0 = 0.0f;
    DiffPair_float_0 _S538;
    (&_S538)->primal_0 = _S369;
    (&_S538)->differential_0 = 0.0f;
    _d_max_0(&_S537, &_S538, _S535.differential_0);
    DiffPair_float_0 _S539;
    (&_S539)->primal_0 = _S393;
    (&_S539)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S539, -0.0f);
    float _S540 = - (-1.0f * - (_S539.differential_0 / _S394));
    float2  _S541 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S542;
    (&_S542)->primal_0 = e2_0;
    (&_S542)->differential_0 = _S541;
    s_bwd_length_impl_1(&_S542, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S543;
    (&_S543)->primal_0 = e1_4;
    (&_S543)->differential_0 = _S541;
    s_bwd_length_impl_1(&_S543, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S544;
    (&_S544)->primal_0 = e0_4;
    (&_S544)->differential_0 = _S541;
    s_bwd_length_impl_1(&_S544, 0.0f);
    DiffPair_float_0 _S545;
    (&_S545)->primal_0 = _S391;
    (&_S545)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S545, -0.0f);
    float _S546 = - _S545.differential_0;
    float2  _S547 = _S543.differential_0 + make_float2 (_S389 * _S546, _S387 * _S545.differential_0);
    float2  _S548 = _S544.differential_0 + make_float2 (_S388 * _S545.differential_0, _S390 * _S546);
    float2  _S549 = - _S542.differential_0 + _S547;
    float _S550 = fx_4 * (_S532.differential_0 + _S536.differential_0 + _S549.x);
    float2  _S551 = make_float2 (_S550, fy_4 * (_S524.differential_0 + _S528.differential_0 + _S549.y)) + make_float2 (dist_coeffs_4[int(8)] * _S550, dist_coeffs_4[int(9)] * _S550);
    float2  _S552 = _S374 * _S551;
    float _S553 = dist_coeffs_4[int(4)] * _S551.y;
    float _S554 = dist_coeffs_4[int(5)] * _S551.x;
    float _S555 = _S552.x + _S552.y;
    float _S556 = r2_17 * _S555;
    float _S557 = r2_17 * _S556;
    float _S558 = dist_coeffs_4[int(7)] * _S551.y + _S553 + dist_coeffs_4[int(6)] * _S551.x + _S554 + _S378 * _S555 + _S377 * _S556 + _S376 * _S557 + dist_coeffs_4[int(3)] * (r2_17 * _S557);
    float _S559 = v_17 * _S558;
    float _S560 = u_17 * _S558;
    float2  _S561 = (make_float2 (radial_8) * _S551 + make_float2 (_S348 * (v_17 * _S551.y) + _S380 * _S554 + 2.0f * (u_17 * _S554) + _S345 * (v_17 * _S551.x) + _S560 + _S560, _S382 * _S553 + 2.0f * (v_17 * _S553) + _S381 * _S551.y + _S379 * _S551.x + _S559 + _S559)) / _S375;
    float2  _S562 = _S371 * - _S561;
    float2  _S563 = _S373 * _S561;
    float2  _S564 = - _S547 + _S548;
    float _S565 = fx_4 * (_S534.differential_0 + _S538.differential_0 + _S564.x);
    float2  _S566 = make_float2 (_S565, fy_4 * (_S526.differential_0 + _S530.differential_0 + _S564.y)) + make_float2 (dist_coeffs_4[int(8)] * _S565, dist_coeffs_4[int(9)] * _S565);
    float2  _S567 = _S358 * _S566;
    float _S568 = dist_coeffs_4[int(4)] * _S566.y;
    float _S569 = dist_coeffs_4[int(5)] * _S566.x;
    float _S570 = _S567.x + _S567.y;
    float _S571 = r2_16 * _S570;
    float _S572 = r2_16 * _S571;
    float _S573 = dist_coeffs_4[int(7)] * _S566.y + _S568 + dist_coeffs_4[int(6)] * _S566.x + _S569 + _S362 * _S570 + _S361 * _S571 + _S360 * _S572 + dist_coeffs_4[int(3)] * (r2_16 * _S572);
    float _S574 = v_16 * _S573;
    float _S575 = u_16 * _S573;
    float2  _S576 = (make_float2 (radial_7) * _S566 + make_float2 (_S348 * (v_16 * _S566.y) + _S364 * _S569 + 2.0f * (u_16 * _S569) + _S345 * (v_16 * _S566.x) + _S575 + _S575, _S366 * _S568 + 2.0f * (v_16 * _S568) + _S365 * _S566.y + _S363 * _S566.x + _S574 + _S574)) / _S359;
    float2  _S577 = _S355 * - _S576;
    float2  _S578 = _S357 * _S576;
    float _S579 = _S577.x + _S577.y;
    float2  _S580 = _S542.differential_0 + - _S548;
    float _S581 = fx_4 * (_S533.differential_0 + _S537.differential_0 + _S580.x);
    float2  _S582 = make_float2 (_S581, fy_4 * (_S525.differential_0 + _S529.differential_0 + _S580.y)) + make_float2 (dist_coeffs_4[int(8)] * _S581, dist_coeffs_4[int(9)] * _S581);
    float2  _S583 = _S340 * _S582;
    float _S584 = dist_coeffs_4[int(4)] * _S582.y;
    float _S585 = dist_coeffs_4[int(5)] * _S582.x;
    float _S586 = _S583.x + _S583.y;
    float _S587 = r2_15 * _S586;
    float _S588 = r2_15 * _S587;
    float _S589 = dist_coeffs_4[int(7)] * _S582.y + _S584 + dist_coeffs_4[int(6)] * _S582.x + _S585 + _S344 * _S586 + _S343 * _S587 + _S342 * _S588 + dist_coeffs_4[int(3)] * (r2_15 * _S588);
    float _S590 = v_15 * _S589;
    float _S591 = u_15 * _S589;
    float2  _S592 = (make_float2 (radial_6) * _S582 + make_float2 (_S348 * (v_15 * _S582.y) + _S347 * _S585 + 2.0f * (u_15 * _S585) + _S345 * (v_15 * _S582.x) + _S591 + _S591, _S350 * _S584 + 2.0f * (v_15 * _S584) + _S349 * _S582.y + _S346 * _S582.x + _S590 + _S590)) / _S341;
    float2  _S593 = _S337 * - _S592;
    float2  _S594 = _S339 * _S592;
    float _S595 = _S593.x + _S593.y;
    float3  _S596 = _S444.differential_0 + _S522 + make_float3 (_S563.x, _S563.y, _S562.x + _S562.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S597;
    (&_S597)->primal_0 = R_4;
    (&_S597)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S598;
    (&_S598)->primal_0 = vert2_4;
    (&_S598)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S597, &_S598, _S596);
    float3  _S599 = _S443.differential_0 + _S522 + make_float3 (_S578.x, _S578.y, _S579);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S600;
    (&_S600)->primal_0 = R_4;
    (&_S600)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S601;
    (&_S601)->primal_0 = vert1_4;
    (&_S601)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S600, &_S601, _S599);
    float3  _S602 = _S445 + _S446 + _S522 + make_float3 (_S594.x, _S594.y, _S595);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S603;
    (&_S603)->primal_0 = R_4;
    (&_S603)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S604;
    (&_S604)->primal_0 = vert0_4;
    (&_S604)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S603, &_S604, _S602);
    float3  _S605 = v_verts_0[int(2)] + _S598.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S606;
    (&_S606)->primal_0 = _S331;
    (&_S606)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S607;
    (&_S607)->primal_0 = _S336;
    (&_S607)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S606, &_S607, _S605);
    float _S608 = - _S607.differential_0.y;
    float _S609 = _S335 * _S607.differential_0.x;
    float _S610 = - (_S327 * _S607.differential_0.x);
    float3  _S611 = v_verts_0[int(1)] + _S601.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S612;
    (&_S612)->primal_0 = _S331;
    (&_S612)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S613;
    (&_S613)->primal_0 = _S334;
    (&_S613)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S612, &_S613, _S611);
    float _S614 = _S327 * _S613.differential_0.x;
    float _S615 = _S333 * _S613.differential_0.x;
    float3  _S616 = v_verts_0[int(0)] + _S604.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S617;
    (&_S617)->primal_0 = _S331;
    (&_S617)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S618;
    (&_S618)->primal_0 = _S332;
    (&_S618)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S617, &_S618, _S616);
    Matrix<float, 3, 3>  _S619 = transpose_1(_S606.differential_0 + _S612.differential_0 + _S617.differential_0);
    float _S620 = 2.0f * - _S619.rows[int(2)].z;
    float _S621 = 2.0f * _S619.rows[int(2)].y;
    float _S622 = 2.0f * _S619.rows[int(2)].x;
    float _S623 = 2.0f * _S619.rows[int(1)].z;
    float _S624 = 2.0f * - _S619.rows[int(1)].y;
    float _S625 = 2.0f * _S619.rows[int(1)].x;
    float _S626 = 2.0f * _S619.rows[int(0)].z;
    float _S627 = 2.0f * _S619.rows[int(0)].y;
    float _S628 = 2.0f * - _S619.rows[int(0)].x;
    float _S629 = - _S625 + _S627;
    float _S630 = _S622 + - _S626;
    float _S631 = - _S621 + _S623;
    float _S632 = _S621 + _S623;
    float _S633 = _S622 + _S626;
    float _S634 = _S625 + _S627;
    float _S635 = quat_4.w * (_S624 + _S628);
    float _S636 = quat_4.z * (_S620 + _S628);
    float _S637 = quat_4.y * (_S620 + _S624);
    float _S638 = quat_4.x * _S629 + quat_4.z * _S632 + quat_4.y * _S633 + _S635 + _S635;
    float _S639 = quat_4.x * _S630 + quat_4.w * _S632 + quat_4.y * _S634 + _S636 + _S636;
    float _S640 = quat_4.x * _S631 + quat_4.w * _S633 + quat_4.z * _S634 + _S637 + _S637;
    float _S641 = quat_4.w * _S629 + quat_4.z * _S630 + quat_4.y * _S631;
    float _S642 = _S610 + _S614;
    float _S643 = 0.5f * - _S642;
    float _S644 = _S608 + _S613.differential_0.y;
    DiffPair_float_0 _S645;
    (&_S645)->primal_0 = _S328;
    (&_S645)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S645, _S644);
    float _S646 = _S643 + _S645.differential_0;
    float _S647 = _S609 + _S615 + _S618.differential_0.x;
    DiffPair_float_0 _S648;
    (&_S648)->primal_0 = _S326;
    (&_S648)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S648, _S647);
    float _S649 = _S643 + _S648.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S650;
    (&_S650)->primal_0 = R_4;
    (&_S650)->differential_0 = _S518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S651;
    (&_S651)->primal_0 = mean_4;
    (&_S651)->differential_0 = _S438;
    s_bwd_prop_mul_0(&_S650, &_S651, _S440.differential_0);
    float3  _S652 = _S520.differential_0 + _S596 + _S599 + _S602 + _S440.differential_0;
    Matrix<float, 3, 3>  _S653 = _S521 + _S597.differential_0 + _S600.differential_0 + _S603.differential_0 + _S650.differential_0;
    FixedArray<float3 , 2>  _S654;
    _S654[int(0)] = _S438;
    _S654[int(1)] = _S438;
    _S654[int(1)] = _S452;
    _S654[int(0)] = _S457;
    FixedArray<float3 , 16>  _S655;
    _S655[int(0)] = _S438;
    _S655[int(1)] = _S438;
    _S655[int(2)] = _S438;
    _S655[int(3)] = _S438;
    _S655[int(4)] = _S438;
    _S655[int(5)] = _S438;
    _S655[int(6)] = _S438;
    _S655[int(7)] = _S438;
    _S655[int(8)] = _S438;
    _S655[int(9)] = _S438;
    _S655[int(10)] = _S438;
    _S655[int(11)] = _S438;
    _S655[int(12)] = _S438;
    _S655[int(13)] = _S438;
    _S655[int(14)] = _S438;
    _S655[int(15)] = _S438;
    _S655[int(15)] = _S459;
    _S655[int(14)] = _S461;
    _S655[int(13)] = _S463;
    _S655[int(12)] = _S465;
    _S655[int(11)] = _S467;
    _S655[int(10)] = _S469;
    _S655[int(9)] = _S471;
    _S655[int(8)] = _S479;
    _S655[int(7)] = _S481;
    _S655[int(6)] = _S483;
    _S655[int(5)] = _S485;
    _S655[int(4)] = _S487;
    _S655[int(3)] = _S498;
    _S655[int(2)] = _S500;
    _S655[int(1)] = _S502;
    _S655[int(0)] = _S515;
    float2  _S656 = make_float2 (0.0f, _S540);
    float3  _S657 = make_float3 (_S649, _S646, _S642);
    float4  _S658 = make_float4 (0.0f);
    *&((&_S658)->w) = _S638;
    *&((&_S658)->z) = _S639;
    *&((&_S658)->y) = _S640;
    *&((&_S658)->x) = _S641;
    *v_mean_0 = _S516 + _S605 + _S611 + _S616 + _S651.differential_0;
    *v_quat_0 = _S658;
    *v_scale_0 = _S657;
    *v_hardness_0 = _S656;
    *v_sh_coeffs_0 = _S655;
    *v_ch_coeffs_0 = _S654;
    *v_R_0 = _S653;
    *v_t_0 = _S652;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S659, float _S660)
{
    return (F32_atan2((_S659), (_S660)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S661, DiffPair_float_0 * _S662, float _S663)
{
    _d_atan2_0(_S661, _S662, _S663);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S664 = scale_5.x;
    float _S665 = s_primal_ctx_exp_0(_S664);
    float _S666 = scale_5.y;
    float _S667 = s_primal_ctx_exp_0(_S666);
    float sz_5 = scale_5.z - 0.5f * (_S664 + _S666);
    float _S668 = quat_5.y;
    float x2_5 = _S668 * _S668;
    float y2_5 = quat_5.z * quat_5.z;
    float z2_10 = quat_5.w * quat_5.w;
    float xy_5 = quat_5.y * quat_5.z;
    float xz_5 = quat_5.y * quat_5.w;
    float yz_5 = quat_5.z * quat_5.w;
    float wx_5 = quat_5.x * quat_5.y;
    float wy_5 = quat_5.x * quat_5.z;
    float wz_5 = quat_5.x * quat_5.w;
    Matrix<float, 3, 3>  _S669 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S670 = make_float3 (_S665, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S669, _S670) + mean_5;
    float _S671 = -0.5f + sz_5;
    float3  _S672 = make_float3 (_S665 * _S671, _S667, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S669, _S672) + mean_5;
    float _S673 = -0.5f - sz_5;
    float3  _S674 = make_float3 (_S665 * _S673, - _S667, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S669, _S674) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S675 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S676 = length_0(_S675);
    float _S677 = vert0_c_5.z;
    float _S678 = s_primal_ctx_atan2_0(_S676, _S677);
    bool _S679 = _S678 < 0.00100000004749745f;
    float k_2;
    float _S680;
    float _S681;
    float _S682;
    if(_S679)
    {
        float _S683 = 1.0f - _S678 * _S678 / 3.0f;
        float _S684 = _S677 * _S677;
        k_2 = _S683 / _S677;
        _S680 = _S684;
        _S681 = _S683;
        _S682 = 0.0f;
    }
    else
    {
        float _S685 = _S676 * _S676;
        k_2 = _S678 / _S676;
        _S680 = 0.0f;
        _S681 = 0.0f;
        _S682 = _S685;
    }
    float2  _S686 = make_float2 (k_2);
    float2  _S687 = _S675 * make_float2 (k_2);
    float u_18 = _S687.x;
    float v_18 = _S687.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S688 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S689 = dist_coeffs_5[int(1)] + r2_18 * _S688;
    float _S690 = dist_coeffs_5[int(0)] + r2_18 * _S689;
    float radial_9 = 1.0f + r2_18 * _S690;
    float _S691 = 2.0f * dist_coeffs_5[int(4)];
    float _S692 = _S691 * u_18;
    float _S693 = 2.0f * u_18;
    float _S694 = 2.0f * dist_coeffs_5[int(5)];
    float _S695 = _S694 * u_18;
    float _S696 = 2.0f * v_18;
    float2  _S697 = _S687 * make_float2 (radial_9) + make_float2 (_S692 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S693 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S695 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S696 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S698 = _S697 + make_float2 (dist_coeffs_5[int(8)] * _S697.x + dist_coeffs_5[int(9)] * _S697.y, 0.0f);
    float _S699 = fx_5 * _S698.x + cx_5;
    float _S700 = fy_5 * _S698.y + cy_5;
    float2  uv0_7 = make_float2 (_S699, _S700);
    float2  _S701 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S702 = length_0(_S701);
    float _S703 = vert1_c_5.z;
    float _S704 = s_primal_ctx_atan2_0(_S702, _S703);
    bool _S705 = _S704 < 0.00100000004749745f;
    float _S706;
    float _S707;
    float _S708;
    if(_S705)
    {
        float _S709 = 1.0f - _S704 * _S704 / 3.0f;
        float _S710 = _S703 * _S703;
        k_2 = _S709 / _S703;
        _S706 = _S710;
        _S707 = _S709;
        _S708 = 0.0f;
    }
    else
    {
        float _S711 = _S702 * _S702;
        k_2 = _S704 / _S702;
        _S706 = 0.0f;
        _S707 = 0.0f;
        _S708 = _S711;
    }
    float2  _S712 = make_float2 (k_2);
    float2  _S713 = _S701 * make_float2 (k_2);
    float u_19 = _S713.x;
    float v_19 = _S713.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S714 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S715 = dist_coeffs_5[int(1)] + r2_19 * _S714;
    float _S716 = dist_coeffs_5[int(0)] + r2_19 * _S715;
    float radial_10 = 1.0f + r2_19 * _S716;
    float _S717 = _S691 * u_19;
    float _S718 = 2.0f * u_19;
    float _S719 = _S694 * u_19;
    float _S720 = 2.0f * v_19;
    float2  _S721 = _S713 * make_float2 (radial_10) + make_float2 (_S717 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S718 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S719 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S720 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S722 = _S721 + make_float2 (dist_coeffs_5[int(8)] * _S721.x + dist_coeffs_5[int(9)] * _S721.y, 0.0f);
    float _S723 = fx_5 * _S722.x + cx_5;
    float _S724 = fy_5 * _S722.y + cy_5;
    float2  uv1_7 = make_float2 (_S723, _S724);
    float2  _S725 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S726 = length_0(_S725);
    float _S727 = vert2_c_5.z;
    float _S728 = s_primal_ctx_atan2_0(_S726, _S727);
    bool _S729 = _S728 < 0.00100000004749745f;
    float _S730;
    float _S731;
    float _S732;
    if(_S729)
    {
        float _S733 = 1.0f - _S728 * _S728 / 3.0f;
        float _S734 = _S727 * _S727;
        k_2 = _S733 / _S727;
        _S730 = _S734;
        _S731 = _S733;
        _S732 = 0.0f;
    }
    else
    {
        float _S735 = _S726 * _S726;
        k_2 = _S728 / _S726;
        _S730 = 0.0f;
        _S731 = 0.0f;
        _S732 = _S735;
    }
    float2  _S736 = make_float2 (k_2);
    float2  _S737 = _S725 * make_float2 (k_2);
    float u_20 = _S737.x;
    float v_20 = _S737.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S738 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S739 = dist_coeffs_5[int(1)] + r2_20 * _S738;
    float _S740 = dist_coeffs_5[int(0)] + r2_20 * _S739;
    float radial_11 = 1.0f + r2_20 * _S740;
    float _S741 = _S691 * u_20;
    float _S742 = 2.0f * u_20;
    float _S743 = _S694 * u_20;
    float _S744 = 2.0f * v_20;
    float2  _S745 = _S737 * make_float2 (radial_11) + make_float2 (_S741 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S742 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S743 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S744 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S746 = _S745 + make_float2 (dist_coeffs_5[int(8)] * _S745.x + dist_coeffs_5[int(9)] * _S745.y, 0.0f);
    float _S747 = fx_5 * _S746.x + cx_5;
    float _S748 = fy_5 * _S746.y + cy_5;
    float2  uv2_7 = make_float2 (_S747, _S748);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S749 = e0_5.x;
    float _S750 = e1_5.y;
    float _S751 = e0_5.y;
    float _S752 = e1_5.x;
    float _S753 = _S749 * _S750 - _S751 * _S752;
    float _S754 = 1.0f - hardness_5.y;
    float _S755 = -1.0f / _S754;
    float _S756 = _S754 * _S754;
    float _S757 = (F32_max((_S699), (_S723)));
    float _S758 = (F32_min((_S699), (_S723)));
    float _S759 = (F32_max((_S700), (_S724)));
    float _S760 = (F32_min((_S700), (_S724)));
    Matrix<float, 3, 3>  _S761 = transpose_1(R_5);
    float3  _S762 = mean_5 - - s_primal_ctx_mul_0(_S761, t_5);
    float _S763 = _S762.x;
    float _S764 = _S762.y;
    float _S765 = _S762.z;
    float _S766 = _S763 * _S763 + _S764 * _S764 + _S765 * _S765;
    float _S767 = s_primal_ctx_sqrt_0(_S766);
    float x_17 = _S763 / _S767;
    float3  _S768 = make_float3 (x_17);
    float _S769 = _S767 * _S767;
    float y_8 = _S764 / _S767;
    float z_5 = _S765 / _S767;
    float3  _S770 = make_float3 (z_5);
    float _S771 = - y_8;
    float3  _S772 = make_float3 (_S771);
    float z2_11 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_17 * x_17 - y_8 * y_8;
    float _S773 = 2.0f * x_17;
    float fS1_5 = _S773 * y_8;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S774 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_17;
    float3  _S775 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_8;
    float3  _S776 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S777 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S778 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    float _S779 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_5 * _S779;
    float3  _S780 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_17;
    float3  _S781 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_8;
    float3  _S782 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S783 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S784 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_17 * fC1_5 - y_8 * fS1_5);
    float3  _S785 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_17 * fS1_5 + y_8 * fC1_5);
    float3  _S786 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S771) * sh_coeffs_5[int(1)] + make_float3 (z_5) * sh_coeffs_5[int(2)] - make_float3 (x_17) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]);
    float3  _S787 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S788 = make_float3 (0.0f);
    float3  _S789 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S790 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S791 = make_float3 (_S790);
    float3  _S792 = ch_coeffs_5[int(1)] * make_float3 (_S790);
    float3  _S793 = _S789 + _S792 + make_float3 (0.5f);
    float3  _S794 = _S789 - _S792 + make_float3 (0.5f);
    float3  _S795 = vert1_c_5 - vert0_c_5;
    float3  _S796 = vert2_c_5 - vert0_c_5;
    float3  _S797 = s_primal_ctx_cross_0(_S795, _S796);
    float3  _S798 = normalize_0(_S797);
    float3  _S799 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S798, mean_c_5)))))) * v_normal_1;
    float3  _S800 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S801;
    (&_S801)->primal_0 = _S798;
    (&_S801)->differential_0 = _S800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S802;
    (&_S802)->primal_0 = mean_c_5;
    (&_S802)->differential_0 = _S800;
    s_bwd_prop_dot_0(&_S801, &_S802, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S803 = _S802;
    float3  _S804 = _S799 + _S801.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S805;
    (&_S805)->primal_0 = _S797;
    (&_S805)->differential_0 = _S800;
    s_bwd_normalize_impl_0(&_S805, _S804);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S806;
    (&_S806)->primal_0 = _S795;
    (&_S806)->differential_0 = _S800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S807;
    (&_S807)->primal_0 = _S796;
    (&_S807)->differential_0 = _S800;
    s_bwd_prop_cross_0(&_S806, &_S807, _S805.differential_0);
    float3  _S808 = - _S807.differential_0;
    float3  _S809 = - _S806.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S810;
    (&_S810)->primal_0 = _S794;
    (&_S810)->differential_0 = _S800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S811;
    (&_S811)->primal_0 = _S788;
    (&_S811)->differential_0 = _S800;
    s_bwd_prop_max_0(&_S810, &_S811, v_rgbs_1[int(2)]);
    float3  _S812 = - _S810.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S813;
    (&_S813)->primal_0 = _S793;
    (&_S813)->differential_0 = _S800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S814;
    (&_S814)->primal_0 = _S788;
    (&_S814)->differential_0 = _S800;
    s_bwd_prop_max_0(&_S813, &_S814, v_rgbs_1[int(1)]);
    float3  _S815 = _S791 * (_S812 + _S813.differential_0);
    float3  _S816 = _S810.differential_0 + _S813.differential_0;
    float3  _S817 = make_float3 (0.5f) * - _S816;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S818;
    (&_S818)->primal_0 = _S787;
    (&_S818)->differential_0 = _S800;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S819;
    (&_S819)->primal_0 = _S788;
    (&_S819)->differential_0 = _S800;
    s_bwd_prop_max_0(&_S818, &_S819, v_rgbs_1[int(0)]);
    float3  _S820 = _S817 + _S818.differential_0;
    float3  _S821 = _S816 + _S818.differential_0;
    float3  _S822 = _S785 * _S821;
    float3  _S823 = sh_coeffs_5[int(15)] * _S821;
    float3  _S824 = _S783 * _S821;
    float3  _S825 = sh_coeffs_5[int(14)] * _S821;
    float3  _S826 = _S781 * _S821;
    float3  _S827 = sh_coeffs_5[int(13)] * _S821;
    float3  _S828 = _S780 * _S821;
    float3  _S829 = sh_coeffs_5[int(12)] * _S821;
    float3  _S830 = _S782 * _S821;
    float3  _S831 = sh_coeffs_5[int(11)] * _S821;
    float3  _S832 = _S784 * _S821;
    float3  _S833 = sh_coeffs_5[int(10)] * _S821;
    float3  _S834 = _S786 * _S821;
    float3  _S835 = sh_coeffs_5[int(9)] * _S821;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S835.x + _S835.y + _S835.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S823.x + _S823.y + _S823.z);
    float _S836 = _S833.x + _S833.y + _S833.z;
    float _S837 = _S825.x + _S825.y + _S825.z;
    float _S838 = _S831.x + _S831.y + _S831.z;
    float _S839 = _S827.x + _S827.y + _S827.z;
    float _S840 = _S829.x + _S829.y + _S829.z;
    float _S841 = - s_diff_fC2_T_1;
    float3  _S842 = _S777 * _S821;
    float3  _S843 = sh_coeffs_5[int(8)] * _S821;
    float3  _S844 = _S775 * _S821;
    float3  _S845 = sh_coeffs_5[int(7)] * _S821;
    float3  _S846 = _S774 * _S821;
    float3  _S847 = sh_coeffs_5[int(6)] * _S821;
    float3  _S848 = _S776 * _S821;
    float3  _S849 = sh_coeffs_5[int(5)] * _S821;
    float3  _S850 = _S778 * _S821;
    float3  _S851 = sh_coeffs_5[int(4)] * _S821;
    float _S852 = _S849.x + _S849.y + _S849.z;
    float _S853 = _S845.x + _S845.y + _S845.z;
    float _S854 = fTmp1B_5 * _S836 + x_17 * s_diff_fS2_T_1 + y_8 * _S841 + 0.54627424478530884f * (_S851.x + _S851.y + _S851.z);
    float _S855 = fTmp1B_5 * _S837 + y_8 * s_diff_fS2_T_1 + x_17 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S843.x + _S843.y + _S843.z);
    float _S856 = y_8 * - _S855;
    float _S857 = x_17 * _S855;
    float _S858 = z_5 * (1.86588168144226074f * (z_5 * _S840) + -2.28522896766662598f * (y_8 * _S838 + x_17 * _S839) + 0.94617468118667603f * (_S847.x + _S847.y + _S847.z));
    float3  _S859 = make_float3 (0.48860251903533936f) * _S821;
    float3  _S860 = - _S859;
    float3  _S861 = _S768 * _S860;
    float3  _S862 = sh_coeffs_5[int(3)] * _S860;
    float3  _S863 = _S770 * _S859;
    float3  _S864 = sh_coeffs_5[int(2)] * _S859;
    float3  _S865 = _S772 * _S859;
    float3  _S866 = sh_coeffs_5[int(1)] * _S859;
    float _S867 = (_S779 * _S840 + 1.44530570507049561f * (fS1_5 * _S836 + fC1_5 * _S837) + -1.09254848957061768f * (y_8 * _S852 + x_17 * _S853) + _S858 + _S858 + _S864.x + _S864.y + _S864.z) / _S769;
    float _S868 = _S767 * _S867;
    float _S869 = (fTmp0C_5 * _S838 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S841 + fTmp0B_5 * _S852 + _S773 * _S854 + _S856 + _S856 + - (_S866.x + _S866.y + _S866.z)) / _S769;
    float _S870 = _S767 * _S869;
    float _S871 = (fTmp0C_5 * _S839 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S853 + 2.0f * (y_8 * _S854) + _S857 + _S857 + _S862.x + _S862.y + _S862.z) / _S769;
    float _S872 = _S767 * _S871;
    float _S873 = _S765 * - _S867 + _S764 * - _S869 + _S763 * - _S871;
    DiffPair_float_0 _S874;
    (&_S874)->primal_0 = _S766;
    (&_S874)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S874, _S873);
    float _S875 = _S765 * _S874.differential_0;
    float _S876 = _S764 * _S874.differential_0;
    float _S877 = _S763 * _S874.differential_0;
    float3  _S878 = make_float3 (0.282094806432724f) * _S821;
    float3  _S879 = make_float3 (_S872 + _S877 + _S877, _S870 + _S876 + _S876, _S868 + _S875 + _S875);
    float3  _S880 = - - _S879;
    Matrix<float, 3, 3>  _S881 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S882;
    (&_S882)->primal_0 = _S761;
    (&_S882)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S883;
    (&_S883)->primal_0 = t_5;
    (&_S883)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S882, &_S883, _S880);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S884 = _S883;
    Matrix<float, 3, 3>  _S885 = transpose_1(_S882.differential_0);
    float3  _S886 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_1);
    DiffPair_float_0 _S887;
    (&_S887)->primal_0 = _S760;
    (&_S887)->differential_0 = 0.0f;
    DiffPair_float_0 _S888;
    (&_S888)->primal_0 = _S748;
    (&_S888)->differential_0 = 0.0f;
    _d_min_0(&_S887, &_S888, 0.0f);
    DiffPair_float_0 _S889;
    (&_S889)->primal_0 = _S700;
    (&_S889)->differential_0 = 0.0f;
    DiffPair_float_0 _S890;
    (&_S890)->primal_0 = _S724;
    (&_S890)->differential_0 = 0.0f;
    _d_min_0(&_S889, &_S890, _S887.differential_0);
    DiffPair_float_0 _S891;
    (&_S891)->primal_0 = _S759;
    (&_S891)->differential_0 = 0.0f;
    DiffPair_float_0 _S892;
    (&_S892)->primal_0 = _S748;
    (&_S892)->differential_0 = 0.0f;
    _d_max_0(&_S891, &_S892, 0.0f);
    DiffPair_float_0 _S893;
    (&_S893)->primal_0 = _S700;
    (&_S893)->differential_0 = 0.0f;
    DiffPair_float_0 _S894;
    (&_S894)->primal_0 = _S724;
    (&_S894)->differential_0 = 0.0f;
    _d_max_0(&_S893, &_S894, _S891.differential_0);
    DiffPair_float_0 _S895;
    (&_S895)->primal_0 = _S758;
    (&_S895)->differential_0 = 0.0f;
    DiffPair_float_0 _S896;
    (&_S896)->primal_0 = _S747;
    (&_S896)->differential_0 = 0.0f;
    _d_min_0(&_S895, &_S896, 0.0f);
    DiffPair_float_0 _S897;
    (&_S897)->primal_0 = _S699;
    (&_S897)->differential_0 = 0.0f;
    DiffPair_float_0 _S898;
    (&_S898)->primal_0 = _S723;
    (&_S898)->differential_0 = 0.0f;
    _d_min_0(&_S897, &_S898, _S895.differential_0);
    DiffPair_float_0 _S899;
    (&_S899)->primal_0 = _S757;
    (&_S899)->differential_0 = 0.0f;
    DiffPair_float_0 _S900;
    (&_S900)->primal_0 = _S747;
    (&_S900)->differential_0 = 0.0f;
    _d_max_0(&_S899, &_S900, 0.0f);
    DiffPair_float_0 _S901;
    (&_S901)->primal_0 = _S699;
    (&_S901)->differential_0 = 0.0f;
    DiffPair_float_0 _S902;
    (&_S902)->primal_0 = _S723;
    (&_S902)->differential_0 = 0.0f;
    _d_max_0(&_S901, &_S902, _S899.differential_0);
    DiffPair_float_0 _S903;
    (&_S903)->primal_0 = _S755;
    (&_S903)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S903, -0.0f);
    float _S904 = - (-1.0f * - (_S903.differential_0 / _S756));
    float2  _S905 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S906;
    (&_S906)->primal_0 = e2_1;
    (&_S906)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S906, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S907;
    (&_S907)->primal_0 = e1_5;
    (&_S907)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S907, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S908;
    (&_S908)->primal_0 = e0_5;
    (&_S908)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S908, 0.0f);
    DiffPair_float_0 _S909;
    (&_S909)->primal_0 = _S753;
    (&_S909)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S909, -0.0f);
    float _S910 = - _S909.differential_0;
    float2  _S911 = _S907.differential_0 + make_float2 (_S751 * _S910, _S749 * _S909.differential_0);
    float2  _S912 = - _S911;
    float2  _S913 = _S908.differential_0 + make_float2 (_S750 * _S909.differential_0, _S752 * _S910);
    float2  _S914 = - _S913;
    float2  _S915 = - _S906.differential_0 + _S911;
    float _S916 = fx_5 * (_S896.differential_0 + _S900.differential_0 + _S915.x);
    float2  _S917 = make_float2 (_S916, fy_5 * (_S888.differential_0 + _S892.differential_0 + _S915.y)) + make_float2 (dist_coeffs_5[int(8)] * _S916, dist_coeffs_5[int(9)] * _S916);
    float2  _S918 = _S737 * _S917;
    float _S919 = dist_coeffs_5[int(4)] * _S917.y;
    float _S920 = dist_coeffs_5[int(5)] * _S917.x;
    float _S921 = _S918.x + _S918.y;
    float _S922 = r2_20 * _S921;
    float _S923 = r2_20 * _S922;
    float _S924 = dist_coeffs_5[int(7)] * _S917.y + _S919 + dist_coeffs_5[int(6)] * _S917.x + _S920 + _S740 * _S921 + _S739 * _S922 + _S738 * _S923 + dist_coeffs_5[int(3)] * (r2_20 * _S923);
    float _S925 = v_20 * _S924;
    float _S926 = u_20 * _S924;
    float2  _S927 = make_float2 (radial_11) * _S917 + make_float2 (_S694 * (v_20 * _S917.y) + _S742 * _S920 + 2.0f * (u_20 * _S920) + _S691 * (v_20 * _S917.x) + _S926 + _S926, _S744 * _S919 + 2.0f * (v_20 * _S919) + _S743 * _S917.y + _S741 * _S917.x + _S925 + _S925);
    float3  _S928 = _S806.differential_0 + _S886;
    float3  _S929 = _S808 + _S809 + _S886;
    float3  _S930 = _S807.differential_0 + _S886;
    FixedArray<float3 , 2>  _S931;
    _S931[int(0)] = _S800;
    _S931[int(1)] = _S800;
    _S931[int(1)] = _S815;
    _S931[int(0)] = _S820;
    float3  _S932 = _S931[int(0)];
    float3  _S933 = _S931[int(1)];
    FixedArray<float3 , 16>  _S934;
    _S934[int(0)] = _S800;
    _S934[int(1)] = _S800;
    _S934[int(2)] = _S800;
    _S934[int(3)] = _S800;
    _S934[int(4)] = _S800;
    _S934[int(5)] = _S800;
    _S934[int(6)] = _S800;
    _S934[int(7)] = _S800;
    _S934[int(8)] = _S800;
    _S934[int(9)] = _S800;
    _S934[int(10)] = _S800;
    _S934[int(11)] = _S800;
    _S934[int(12)] = _S800;
    _S934[int(13)] = _S800;
    _S934[int(14)] = _S800;
    _S934[int(15)] = _S800;
    _S934[int(7)] = _S844;
    _S934[int(0)] = _S878;
    _S934[int(1)] = _S865;
    _S934[int(2)] = _S863;
    _S934[int(3)] = _S861;
    _S934[int(4)] = _S850;
    _S934[int(5)] = _S848;
    _S934[int(6)] = _S846;
    _S934[int(15)] = _S822;
    _S934[int(8)] = _S842;
    _S934[int(9)] = _S834;
    _S934[int(10)] = _S832;
    _S934[int(11)] = _S830;
    _S934[int(12)] = _S828;
    _S934[int(13)] = _S826;
    _S934[int(14)] = _S824;
    float3  _S935 = _S934[int(0)];
    float3  _S936 = _S934[int(1)];
    float3  _S937 = _S934[int(2)];
    float3  _S938 = _S934[int(3)];
    float3  _S939 = _S934[int(4)];
    float3  _S940 = _S934[int(5)];
    float3  _S941 = _S934[int(6)];
    float3  _S942 = _S934[int(7)];
    float3  _S943 = _S934[int(8)];
    float3  _S944 = _S934[int(9)];
    float3  _S945 = _S934[int(10)];
    float3  _S946 = _S934[int(11)];
    float3  _S947 = _S934[int(12)];
    float3  _S948 = _S934[int(13)];
    float3  _S949 = _S934[int(14)];
    float3  _S950 = _S934[int(15)];
    float _S951 = _S890.differential_0 + _S894.differential_0;
    float _S952 = _S889.differential_0 + _S893.differential_0;
    float _S953 = _S897.differential_0 + _S901.differential_0;
    float _S954 = _S898.differential_0 + _S902.differential_0;
    float2  _S955 = _S906.differential_0 + _S914;
    float2  _S956 = _S912 + _S913;
    float2  _S957 = make_float2 (0.0f, _S904);
    float2  _S958 = _S725 * _S927;
    float2  _S959 = _S736 * _S927;
    float _S960 = _S958.x + _S958.y;
    if(_S729)
    {
        float _S961 = _S960 / _S730;
        float _S962 = _S731 * - _S961;
        float _S963 = _S728 * (0.3333333432674408f * - (_S727 * _S961));
        k_2 = _S963 + _S963;
        _S730 = _S962;
        _S731 = 0.0f;
    }
    else
    {
        float _S964 = _S960 / _S732;
        float _S965 = _S728 * - _S964;
        k_2 = _S726 * _S964;
        _S730 = 0.0f;
        _S731 = _S965;
    }
    DiffPair_float_0 _S966;
    (&_S966)->primal_0 = _S726;
    (&_S966)->differential_0 = 0.0f;
    DiffPair_float_0 _S967;
    (&_S967)->primal_0 = _S727;
    (&_S967)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S966, &_S967, k_2);
    float _S968 = _S967.differential_0 + _S730;
    float _S969 = _S966.differential_0 + _S731;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S970;
    (&_S970)->primal_0 = _S725;
    (&_S970)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S970, _S969);
    float2  _S971 = _S970.differential_0 + _S959;
    float _S972 = fx_5 * (_S956.x + _S954);
    float2  _S973 = make_float2 (_S972, fy_5 * (_S956.y + _S951)) + make_float2 (dist_coeffs_5[int(8)] * _S972, dist_coeffs_5[int(9)] * _S972);
    float2  _S974 = _S713 * _S973;
    float _S975 = dist_coeffs_5[int(4)] * _S973.y;
    float _S976 = dist_coeffs_5[int(5)] * _S973.x;
    float _S977 = _S974.x + _S974.y;
    float _S978 = r2_19 * _S977;
    float _S979 = r2_19 * _S978;
    float _S980 = dist_coeffs_5[int(7)] * _S973.y + _S975 + dist_coeffs_5[int(6)] * _S973.x + _S976 + _S716 * _S977 + _S715 * _S978 + _S714 * _S979 + dist_coeffs_5[int(3)] * (r2_19 * _S979);
    float _S981 = v_19 * _S980;
    float _S982 = u_19 * _S980;
    float2  _S983 = make_float2 (radial_10) * _S973 + make_float2 (_S694 * (v_19 * _S973.y) + _S718 * _S976 + 2.0f * (u_19 * _S976) + _S691 * (v_19 * _S973.x) + _S982 + _S982, _S720 * _S975 + 2.0f * (v_19 * _S975) + _S719 * _S973.y + _S717 * _S973.x + _S981 + _S981);
    float3  _S984 = _S930 + make_float3 (_S971.x, _S971.y, _S968);
    float2  _S985 = _S701 * _S983;
    float2  _S986 = _S712 * _S983;
    float _S987 = _S985.x + _S985.y;
    if(_S705)
    {
        float _S988 = _S987 / _S706;
        float _S989 = _S707 * - _S988;
        float _S990 = _S704 * (0.3333333432674408f * - (_S703 * _S988));
        k_2 = _S990 + _S990;
        _S706 = _S989;
        _S707 = 0.0f;
    }
    else
    {
        float _S991 = _S987 / _S708;
        float _S992 = _S704 * - _S991;
        k_2 = _S702 * _S991;
        _S706 = 0.0f;
        _S707 = _S992;
    }
    DiffPair_float_0 _S993;
    (&_S993)->primal_0 = _S702;
    (&_S993)->differential_0 = 0.0f;
    DiffPair_float_0 _S994;
    (&_S994)->primal_0 = _S703;
    (&_S994)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S993, &_S994, k_2);
    float _S995 = _S994.differential_0 + _S706;
    float _S996 = _S993.differential_0 + _S707;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S997;
    (&_S997)->primal_0 = _S701;
    (&_S997)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S997, _S996);
    float2  _S998 = _S997.differential_0 + _S986;
    float _S999 = fx_5 * (_S955.x + _S953);
    float2  _S1000 = make_float2 (_S999, fy_5 * (_S955.y + _S952)) + make_float2 (dist_coeffs_5[int(8)] * _S999, dist_coeffs_5[int(9)] * _S999);
    float2  _S1001 = _S687 * _S1000;
    float _S1002 = dist_coeffs_5[int(4)] * _S1000.y;
    float _S1003 = dist_coeffs_5[int(5)] * _S1000.x;
    float _S1004 = _S1001.x + _S1001.y;
    float _S1005 = r2_18 * _S1004;
    float _S1006 = r2_18 * _S1005;
    float _S1007 = dist_coeffs_5[int(7)] * _S1000.y + _S1002 + dist_coeffs_5[int(6)] * _S1000.x + _S1003 + _S690 * _S1004 + _S689 * _S1005 + _S688 * _S1006 + dist_coeffs_5[int(3)] * (r2_18 * _S1006);
    float _S1008 = v_18 * _S1007;
    float _S1009 = u_18 * _S1007;
    float2  _S1010 = make_float2 (radial_9) * _S1000 + make_float2 (_S694 * (v_18 * _S1000.y) + _S693 * _S1003 + 2.0f * (u_18 * _S1003) + _S691 * (v_18 * _S1000.x) + _S1009 + _S1009, _S696 * _S1002 + 2.0f * (v_18 * _S1002) + _S695 * _S1000.y + _S692 * _S1000.x + _S1008 + _S1008);
    float3  _S1011 = _S928 + make_float3 (_S998.x, _S998.y, _S995);
    float2  _S1012 = _S675 * _S1010;
    float2  _S1013 = _S686 * _S1010;
    float _S1014 = _S1012.x + _S1012.y;
    if(_S679)
    {
        float _S1015 = _S1014 / _S680;
        float _S1016 = _S681 * - _S1015;
        float _S1017 = _S678 * (0.3333333432674408f * - (_S677 * _S1015));
        k_2 = _S1017 + _S1017;
        _S680 = _S1016;
        _S681 = 0.0f;
    }
    else
    {
        float _S1018 = _S1014 / _S682;
        float _S1019 = _S678 * - _S1018;
        k_2 = _S676 * _S1018;
        _S680 = 0.0f;
        _S681 = _S1019;
    }
    DiffPair_float_0 _S1020;
    (&_S1020)->primal_0 = _S676;
    (&_S1020)->differential_0 = 0.0f;
    DiffPair_float_0 _S1021;
    (&_S1021)->primal_0 = _S677;
    (&_S1021)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1020, &_S1021, k_2);
    float _S1022 = _S1021.differential_0 + _S680;
    float _S1023 = _S1020.differential_0 + _S681;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1024;
    (&_S1024)->primal_0 = _S675;
    (&_S1024)->differential_0 = _S905;
    s_bwd_length_impl_1(&_S1024, _S1023);
    float2  _S1025 = _S1024.differential_0 + _S1013;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1026;
    (&_S1026)->primal_0 = vert2_c_5;
    (&_S1026)->differential_0 = _S800;
    s_bwd_length_impl_0(&_S1026, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1027;
    (&_S1027)->primal_0 = vert1_c_5;
    (&_S1027)->differential_0 = _S800;
    s_bwd_length_impl_0(&_S1027, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1028;
    (&_S1028)->primal_0 = vert0_c_5;
    (&_S1028)->differential_0 = _S800;
    s_bwd_length_impl_0(&_S1028, 0.0f);
    float3  _S1029 = _S1026.differential_0 + _S984;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1030;
    (&_S1030)->primal_0 = R_5;
    (&_S1030)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1031;
    (&_S1031)->primal_0 = vert2_5;
    (&_S1031)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1030, &_S1031, _S1029);
    float3  _S1032 = _S1027.differential_0 + _S1011;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1033;
    (&_S1033)->primal_0 = R_5;
    (&_S1033)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1034;
    (&_S1034)->primal_0 = vert1_5;
    (&_S1034)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1033, &_S1034, _S1032);
    float3  _S1035 = _S1028.differential_0 + _S929 + make_float3 (_S1025.x, _S1025.y, _S1022);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1036;
    (&_S1036)->primal_0 = R_5;
    (&_S1036)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1037;
    (&_S1037)->primal_0 = vert0_5;
    (&_S1037)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1036, &_S1037, _S1035);
    float3  _S1038 = _S1031.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1039;
    (&_S1039)->primal_0 = _S669;
    (&_S1039)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1040;
    (&_S1040)->primal_0 = _S674;
    (&_S1040)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1039, &_S1040, _S1038);
    float _S1041 = - _S1040.differential_0.y;
    float _S1042 = _S673 * _S1040.differential_0.x;
    float _S1043 = - (_S665 * _S1040.differential_0.x);
    float3  _S1044 = _S1034.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1045;
    (&_S1045)->primal_0 = _S669;
    (&_S1045)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1046;
    (&_S1046)->primal_0 = _S672;
    (&_S1046)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1045, &_S1046, _S1044);
    float _S1047 = _S665 * _S1046.differential_0.x;
    float _S1048 = _S671 * _S1046.differential_0.x;
    float3  _S1049 = _S1037.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1050;
    (&_S1050)->primal_0 = _S669;
    (&_S1050)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1051;
    (&_S1051)->primal_0 = _S670;
    (&_S1051)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1050, &_S1051, _S1049);
    Matrix<float, 3, 3>  _S1052 = transpose_1(_S1039.differential_0 + _S1045.differential_0 + _S1050.differential_0);
    float _S1053 = 2.0f * - _S1052.rows[int(2)].z;
    float _S1054 = 2.0f * _S1052.rows[int(2)].y;
    float _S1055 = 2.0f * _S1052.rows[int(2)].x;
    float _S1056 = 2.0f * _S1052.rows[int(1)].z;
    float _S1057 = 2.0f * - _S1052.rows[int(1)].y;
    float _S1058 = 2.0f * _S1052.rows[int(1)].x;
    float _S1059 = 2.0f * _S1052.rows[int(0)].z;
    float _S1060 = 2.0f * _S1052.rows[int(0)].y;
    float _S1061 = 2.0f * - _S1052.rows[int(0)].x;
    float _S1062 = - _S1058 + _S1060;
    float _S1063 = _S1055 + - _S1059;
    float _S1064 = - _S1054 + _S1056;
    float _S1065 = _S1054 + _S1056;
    float _S1066 = _S1055 + _S1059;
    float _S1067 = _S1058 + _S1060;
    float _S1068 = quat_5.w * (_S1057 + _S1061);
    float _S1069 = quat_5.z * (_S1053 + _S1061);
    float _S1070 = quat_5.y * (_S1053 + _S1057);
    float _S1071 = quat_5.x * _S1062 + quat_5.z * _S1065 + quat_5.y * _S1066 + _S1068 + _S1068;
    float _S1072 = quat_5.x * _S1063 + quat_5.w * _S1065 + quat_5.y * _S1067 + _S1069 + _S1069;
    float _S1073 = quat_5.x * _S1064 + quat_5.w * _S1066 + quat_5.z * _S1067 + _S1070 + _S1070;
    float _S1074 = quat_5.w * _S1062 + quat_5.z * _S1063 + quat_5.y * _S1064;
    float _S1075 = _S1043 + _S1047;
    float _S1076 = 0.5f * - _S1075;
    float _S1077 = _S1041 + _S1046.differential_0.y;
    DiffPair_float_0 _S1078;
    (&_S1078)->primal_0 = _S666;
    (&_S1078)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1078, _S1077);
    float _S1079 = _S1076 + _S1078.differential_0;
    float _S1080 = _S1042 + _S1048 + _S1051.differential_0.x;
    DiffPair_float_0 _S1081;
    (&_S1081)->primal_0 = _S664;
    (&_S1081)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1081, _S1080);
    float _S1082 = _S1076 + _S1081.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1083;
    (&_S1083)->primal_0 = mean_c_5;
    (&_S1083)->differential_0 = _S800;
    s_bwd_length_impl_0(&_S1083, 0.0f);
    float3  _S1084 = _S1083.differential_0 + _S803.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1085;
    (&_S1085)->primal_0 = R_5;
    (&_S1085)->differential_0 = _S881;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1086;
    (&_S1086)->primal_0 = mean_5;
    (&_S1086)->differential_0 = _S800;
    s_bwd_prop_mul_0(&_S1085, &_S1086, _S1084);
    float3  _S1087 = _S1029 + _S1032 + _S1035 + _S1084 + _S884.differential_0;
    Matrix<float, 3, 3>  _S1088 = _S1030.differential_0 + _S1033.differential_0 + _S1036.differential_0 + _S1085.differential_0 + _S885;
    float3  _S1089 = make_float3 (_S1082, _S1079, _S1075);
    float4  _S1090 = make_float4 (0.0f);
    *&((&_S1090)->w) = _S1071;
    *&((&_S1090)->z) = _S1072;
    *&((&_S1090)->y) = _S1073;
    *&((&_S1090)->x) = _S1074;
    float4  _S1091 = _S1090;
    float3  _S1092 = _S1038 + _S1044 + _S1049 + _S1086.differential_0 + _S879;
    *v_mean_1 = _S1092;
    *v_quat_1 = _S1091;
    *v_scale_1 = _S1089;
    *v_hardness_1 = _S957;
    (*v_sh_coeffs_1)[int(0)] = _S935;
    (*v_sh_coeffs_1)[int(1)] = _S936;
    (*v_sh_coeffs_1)[int(2)] = _S937;
    (*v_sh_coeffs_1)[int(3)] = _S938;
    (*v_sh_coeffs_1)[int(4)] = _S939;
    (*v_sh_coeffs_1)[int(5)] = _S940;
    (*v_sh_coeffs_1)[int(6)] = _S941;
    (*v_sh_coeffs_1)[int(7)] = _S942;
    (*v_sh_coeffs_1)[int(8)] = _S943;
    (*v_sh_coeffs_1)[int(9)] = _S944;
    (*v_sh_coeffs_1)[int(10)] = _S945;
    (*v_sh_coeffs_1)[int(11)] = _S946;
    (*v_sh_coeffs_1)[int(12)] = _S947;
    (*v_sh_coeffs_1)[int(13)] = _S948;
    (*v_sh_coeffs_1)[int(14)] = _S949;
    (*v_sh_coeffs_1)[int(15)] = _S950;
    (*v_ch_coeffs_1)[int(0)] = _S932;
    (*v_ch_coeffs_1)[int(1)] = _S933;
    *v_R_1 = _S1088;
    *v_t_1 = _S1087;
    return;
}

inline __device__ bool ray_triangle_intersection_uvt(float3  ray_o_0, float3  ray_d_0, FixedArray<float3 , 3>  verts_4, float * u_21, float * v_21, float * t_6)
{
    float3  v1v0_0 = verts_4[int(1)] - verts_4[int(0)];
    float3  v2v0_0 = verts_4[int(2)] - verts_4[int(0)];
    float3  rov0_0 = ray_o_0 - verts_4[int(0)];
    float3  n_0 = cross_0(v1v0_0, v2v0_0);
    float3  q_0 = cross_0(rov0_0, ray_d_0);
    float d_0 = 1.0f / dot_0(ray_d_0, n_0);
    *u_21 = d_0 * dot_0(- q_0, v2v0_0);
    *v_21 = d_0 * dot_0(q_0, v1v0_0);
    *t_6 = d_0 * dot_0(- n_0, rov0_0);
    bool _S1093;
    if((*u_21) >= 0.0f)
    {
        _S1093 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1093 = false;
    }
    if(_S1093)
    {
        _S1093 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1093 = false;
    }
    if(_S1093)
    {
        _S1093 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1093 = false;
    }
    return _S1093;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1094 = *dpx_12;
    bool _S1095;
    if(((*dpx_12).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1095 = ((*dpx_12).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1095 = false;
    }
    float _S1096;
    if(_S1095)
    {
        _S1096 = dOut_11;
    }
    else
    {
        _S1096 = 0.0f;
    }
    dpx_12->primal_0 = _S1094.primal_0;
    dpx_12->differential_0 = _S1096;
    DiffPair_float_0 _S1097 = *dpMin_0;
    if((_S1094.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1096 = dOut_11;
    }
    else
    {
        _S1096 = 0.0f;
    }
    dpMin_0->primal_0 = _S1097.primal_0;
    dpMin_0->differential_0 = _S1096;
    DiffPair_float_0 _S1098 = *dpMax_0;
    if(((*dpx_12).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1096 = dOut_11;
    }
    else
    {
        _S1096 = 0.0f;
    }
    dpMax_0->primal_0 = _S1098.primal_0;
    dpMax_0->differential_0 = _S1096;
    return;
}

inline __device__ float clamp_0(float x_18, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_18), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_5, float dOut_12)
{
    if(((*dpx_13).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_13->primal_0 = (*dpx_13).primal_0;
        dpx_13->differential_0 = 0.0f;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_13).primal_0), ((*dpy_5).primal_0)));
        DiffPair_float_0 _S1099 = *dpx_13;
        float _S1100 = val_0 * (*dpy_5).primal_0 / (*dpx_13).primal_0 * dOut_12;
        dpx_13->primal_0 = (*dpx_13).primal_0;
        dpx_13->differential_0 = _S1100;
        float _S1101 = val_0 * (F32_log((_S1099.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1101;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle(FixedArray<float3 , 3>  verts_5, float2  hardness_6, float3  ray_o_1, float3  ray_d_1)
{
    float3  v1v0_1 = verts_5[int(1)] - verts_5[int(0)];
    float3  v2v0_1 = verts_5[int(2)] - verts_5[int(0)];
    float3  rov0_1 = ray_o_1 - verts_5[int(0)];
    float3  n_1 = cross_0(v1v0_1, v2v0_1);
    float3  q_1 = cross_0(rov0_1, ray_d_1);
    float d_1 = 1.0f / dot_0(ray_d_1, n_1);
    float u_22 = d_1 * dot_0(- q_1, v2v0_1);
    float v_22 = d_1 * dot_0(q_1, v1v0_1);
    float t_7 = d_1 * dot_0(- n_1, rov0_1);
    bool _S1102;
    if(u_22 >= 0.0f)
    {
        _S1102 = v_22 >= 0.0f;
    }
    else
    {
        _S1102 = false;
    }
    if(_S1102)
    {
        _S1102 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1102 = false;
    }
    if(_S1102)
    {
        _S1102 = t_7 >= 0.0f;
    }
    else
    {
        _S1102 = false;
    }
    if(!_S1102)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1103;
    if(opac_0 < 0.0f)
    {
        _S1103 = 0.0f;
    }
    else
    {
        _S1103 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1103;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1104, float _S1105, float _S1106)
{
    return clamp_0(_S1104, _S1105, _S1106);
}

inline __device__ float s_primal_ctx_pow_0(float _S1107, float _S1108)
{
    return (F32_pow((_S1107), (_S1108)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1109, DiffPair_float_0 * _S1110, float _S1111)
{
    _d_pow_0(_S1109, _S1110, _S1111);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1112, DiffPair_float_0 * _S1113, DiffPair_float_0 * _S1114, float _S1115)
{
    _d_clamp_0(_S1112, _S1113, _S1114, _S1115);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1116 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1117 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1118 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1119 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1120 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1118);
    float d_2 = 1.0f / _S1120;
    float _S1121 = _S1120 * _S1120;
    float3  _S1122 = - _S1119;
    float _S1123 = s_primal_ctx_dot_0(_S1122, v2v0_2);
    float u_23 = d_2 * _S1123;
    float _S1124 = s_primal_ctx_dot_0(_S1119, v1v0_2);
    float v_23 = d_2 * _S1124;
    float3  _S1125 = - _S1118;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1125, rov0_2);
    bool _S1126;
    if(u_23 >= 0.0f)
    {
        _S1126 = v_23 >= 0.0f;
    }
    else
    {
        _S1126 = false;
    }
    if(_S1126)
    {
        _S1126 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1126 = false;
    }
    if(_S1126)
    {
        _S1126 = t_8 >= 0.0f;
    }
    else
    {
        _S1126 = false;
    }
    bool _S1127 = !!_S1126;
    float _S1128;
    float _S1129;
    float _S1130;
    float _S1131;
    float _S1132;
    float _S1133;
    float _S1134;
    float _S1135;
    float _S1136;
    float _S1137;
    float _S1138;
    if(_S1127)
    {
        float _S1139 = (F32_min((u_23), (v_23)));
        float _S1140 = s_primal_ctx_sqrt_0(0.5f);
        float _S1141 = _S1140 * (1.0f - u_23 - v_23);
        float _S1142 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1139), (_S1141))) * _S1142;
        float _S1143 = _S1116.primal_0.y;
        float _S1144 = 1.0f - opac_1;
        float _S1145 = 1.0f - s_primal_ctx_clamp_0(_S1143, 0.0f, 0.99989998340606689f);
        float _S1146 = 1.0f / _S1145;
        float _S1147 = _S1145 * _S1145;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1144, _S1146);
        float o_1 = _S1116.primal_0.x;
        bool _S1148 = opac_1 < 0.0f;
        if(_S1148)
        {
            _S1128 = 0.0f;
        }
        else
        {
            _S1128 = o_1 * w_1;
        }
        _S1126 = _S1148;
        _S1129 = o_1;
        _S1130 = w_1;
        _S1131 = _S1144;
        _S1132 = _S1146;
        _S1133 = _S1147;
        _S1134 = _S1143;
        _S1135 = _S1142;
        _S1136 = _S1139;
        _S1137 = _S1141;
        _S1138 = _S1140;
    }
    else
    {
        _S1126 = false;
        _S1128 = 0.0f;
        _S1129 = 0.0f;
        _S1130 = 0.0f;
        _S1131 = 0.0f;
        _S1132 = 0.0f;
        _S1133 = 0.0f;
        _S1134 = 0.0f;
        _S1135 = 0.0f;
        _S1136 = 0.0f;
        _S1137 = 0.0f;
        _S1138 = 0.0f;
    }
    float2  _S1149 = make_float2 (0.0f);
    float2  _S1150;
    if(_S1127)
    {
        if(_S1126)
        {
            _S1128 = 0.0f;
            _S1129 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1151;
            (&_S1151)->primal_0 = _S1128;
            (&_S1151)->differential_0 = 0.0f;
            DiffPair_float_0 _S1152;
            (&_S1152)->primal_0 = 0.99500000476837158f;
            (&_S1152)->differential_0 = 0.0f;
            _d_min_0(&_S1151, &_S1152, _s_dOut_3);
            float _S1153 = _S1129 * _S1151.differential_0;
            _S1128 = _S1130 * _S1151.differential_0;
            _S1129 = _S1153;
        }
        float _S1154 = - _S1129;
        DiffPair_float_0 _S1155;
        (&_S1155)->primal_0 = _S1131;
        (&_S1155)->differential_0 = 0.0f;
        DiffPair_float_0 _S1156;
        (&_S1156)->primal_0 = _S1132;
        (&_S1156)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1155, &_S1156, _S1154);
        float _S1157 = - - (_S1156.differential_0 / _S1133);
        float s_diff_opac_T_0 = - _S1155.differential_0;
        DiffPair_float_0 _S1158;
        (&_S1158)->primal_0 = _S1134;
        (&_S1158)->differential_0 = 0.0f;
        DiffPair_float_0 _S1159;
        (&_S1159)->primal_0 = 0.0f;
        (&_S1159)->differential_0 = 0.0f;
        DiffPair_float_0 _S1160;
        (&_S1160)->primal_0 = 0.99989998340606689f;
        (&_S1160)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1158, &_S1159, &_S1160, _S1157);
        float _S1161 = _S1135 * s_diff_opac_T_0;
        DiffPair_float_0 _S1162;
        (&_S1162)->primal_0 = _S1136;
        (&_S1162)->differential_0 = 0.0f;
        DiffPair_float_0 _S1163;
        (&_S1163)->primal_0 = _S1137;
        (&_S1163)->differential_0 = 0.0f;
        _d_min_0(&_S1162, &_S1163, _S1161);
        float _S1164 = - (_S1138 * _S1163.differential_0);
        DiffPair_float_0 _S1165;
        (&_S1165)->primal_0 = u_23;
        (&_S1165)->differential_0 = 0.0f;
        DiffPair_float_0 _S1166;
        (&_S1166)->primal_0 = v_23;
        (&_S1166)->differential_0 = 0.0f;
        _d_min_0(&_S1165, &_S1166, _S1162.differential_0);
        float2  _S1167 = make_float2 (_S1128, _S1158.differential_0);
        float _S1168 = _S1164 + _S1166.differential_0;
        _S1128 = _S1164 + _S1165.differential_0;
        _S1129 = _S1168;
        _S1150 = _S1167;
    }
    else
    {
        _S1128 = 0.0f;
        _S1129 = 0.0f;
        _S1150 = _S1149;
    }
    float3  _S1169 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1170;
    (&_S1170)->primal_0 = _S1125;
    (&_S1170)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1171;
    (&_S1171)->primal_0 = rov0_2;
    (&_S1171)->differential_0 = _S1169;
    s_bwd_prop_dot_0(&_S1170, &_S1171, 0.0f);
    float3  _S1172 = - _S1170.differential_0;
    float _S1173 = d_2 * _S1129;
    float _S1174 = _S1124 * _S1129;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1175;
    (&_S1175)->primal_0 = _S1119;
    (&_S1175)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1176;
    (&_S1176)->primal_0 = v1v0_2;
    (&_S1176)->differential_0 = _S1169;
    s_bwd_prop_dot_0(&_S1175, &_S1176, _S1173);
    float _S1177 = d_2 * _S1128;
    float _S1178 = _S1123 * _S1128;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1179;
    (&_S1179)->primal_0 = _S1122;
    (&_S1179)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1180;
    (&_S1180)->primal_0 = v2v0_2;
    (&_S1180)->differential_0 = _S1169;
    s_bwd_prop_dot_0(&_S1179, &_S1180, _S1177);
    float3  _S1181 = - _S1179.differential_0;
    float _S1182 = - ((_S1174 + _S1178) / _S1121);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1183;
    (&_S1183)->primal_0 = _S1117.primal_0;
    (&_S1183)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1184;
    (&_S1184)->primal_0 = _S1118;
    (&_S1184)->differential_0 = _S1169;
    s_bwd_prop_dot_0(&_S1183, &_S1184, _S1182);
    float3  _S1185 = _S1175.differential_0 + _S1181;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1186;
    (&_S1186)->primal_0 = rov0_2;
    (&_S1186)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1187;
    (&_S1187)->primal_0 = _S1117.primal_0;
    (&_S1187)->differential_0 = _S1169;
    s_bwd_prop_cross_0(&_S1186, &_S1187, _S1185);
    float3  _S1188 = _S1172 + _S1184.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1189;
    (&_S1189)->primal_0 = v1v0_2;
    (&_S1189)->differential_0 = _S1169;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1190;
    (&_S1190)->primal_0 = v2v0_2;
    (&_S1190)->differential_0 = _S1169;
    s_bwd_prop_cross_0(&_S1189, &_S1190, _S1188);
    float3  _S1191 = _S1171.differential_0 + _S1186.differential_0;
    float3  _S1192 = _S1180.differential_0 + _S1190.differential_0;
    float3  _S1193 = _S1176.differential_0 + _S1189.differential_0;
    float3  _S1194 = - _S1191 + - _S1192 + - _S1193;
    float3  _S1195 = _S1183.differential_0 + _S1187.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1195;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1191;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1150;
    FixedArray<float3 , 3>  _S1196;
    _S1196[int(0)] = _S1169;
    _S1196[int(1)] = _S1169;
    _S1196[int(2)] = _S1169;
    _S1196[int(2)] = _S1192;
    _S1196[int(0)] = _S1194;
    _S1196[int(1)] = _S1193;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1196;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1197, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1198, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1199, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1200, float _S1201)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1197, _S1198, _S1199, _S1200, _S1201);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1202 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1203 = { _S1202, _S1202, _S1202 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1203;
    float2  _S1204 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1204;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1202;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1202;
    s_bwd_evaluate_alpha_opaque_triangle_0(&dp_verts_0, &dp_hardness_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_verts_2 = (&dp_verts_0)->differential_0;
    *v_hardness_2 = dp_hardness_0.differential_0;
    *v_ray_o_0 = dp_ray_o_0.differential_0;
    *v_ray_d_0 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ float evaluate_sorting_depth_opaque_triangle(FixedArray<float3 , 3>  verts_7, FixedArray<float3 , 3>  rgbs_4, float3  ray_o_3, float3  ray_d_3)
{
    float3  n_2 = cross_0(verts_7[int(1)] - verts_7[int(0)], verts_7[int(2)] - verts_7[int(0)]);
    return 1.0f / dot_0(ray_d_3, n_2) * dot_0(- n_2, ray_o_3 - verts_7[int(0)]);
}

inline __device__ void evaluate_color_opaque_triangle(FixedArray<float3 , 3>  verts_8, FixedArray<float3 , 3>  rgbs_5, float3  ray_o_4, float3  ray_d_4, float3  * color_6, float * depth_4)
{
    float3  v1v0_3 = verts_8[int(1)] - verts_8[int(0)];
    float3  v2v0_3 = verts_8[int(2)] - verts_8[int(0)];
    float3  rov0_3 = ray_o_4 - verts_8[int(0)];
    float3  n_3 = cross_0(v1v0_3, v2v0_3);
    float3  q_2 = cross_0(rov0_3, ray_d_4);
    float d_3 = 1.0f / dot_0(ray_d_4, n_3);
    float u_24 = d_3 * dot_0(- q_2, v2v0_3);
    float v_24 = d_3 * dot_0(q_2, v1v0_3);
    *depth_4 = d_3 * dot_0(- n_3, rov0_3);
    *color_6 = rgbs_5[int(0)] * make_float3 (1.0f - u_24 - v_24) + rgbs_5[int(1)] * make_float3 (u_24) + rgbs_5[int(2)] * make_float3 (v_24);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpcolor_0, float dpdepth_0)
{
    float3  v1v0_4 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_4 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_4 = (*dpray_o_1).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S1205 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1206 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1207 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1205);
    float d_4 = 1.0f / _S1207;
    float _S1208 = _S1207 * _S1207;
    float3  _S1209 = - _S1206;
    float _S1210 = s_primal_ctx_dot_0(_S1209, v2v0_4);
    float u_25 = d_4 * _S1210;
    float _S1211 = s_primal_ctx_dot_0(_S1206, v1v0_4);
    float v_25 = d_4 * _S1211;
    float3  _S1212 = - _S1205;
    float3  _S1213 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1214 = make_float3 (v_25) * dpcolor_0;
    float3  _S1215 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1216 = make_float3 (u_25) * dpcolor_0;
    float3  _S1217 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1218 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1219 = - (_S1217.x + _S1217.y + _S1217.z);
    float _S1220 = d_4 * dpdepth_0;
    float _S1221 = s_primal_ctx_dot_0(_S1212, rov0_4) * dpdepth_0;
    float3  _S1222 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1223;
    (&_S1223)->primal_0 = _S1212;
    (&_S1223)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1224;
    (&_S1224)->primal_0 = rov0_4;
    (&_S1224)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1223, &_S1224, _S1220);
    float3  _S1225 = - _S1223.differential_0;
    float _S1226 = _S1219 + _S1213.x + _S1213.y + _S1213.z;
    float _S1227 = d_4 * _S1226;
    float _S1228 = _S1211 * _S1226;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1229;
    (&_S1229)->primal_0 = _S1206;
    (&_S1229)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1230;
    (&_S1230)->primal_0 = v1v0_4;
    (&_S1230)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1229, &_S1230, _S1227);
    float _S1231 = _S1219 + _S1215.x + _S1215.y + _S1215.z;
    float _S1232 = d_4 * _S1231;
    float _S1233 = _S1210 * _S1231;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
    (&_S1234)->primal_0 = _S1209;
    (&_S1234)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1235;
    (&_S1235)->primal_0 = v2v0_4;
    (&_S1235)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1234, &_S1235, _S1232);
    float3  _S1236 = - _S1234.differential_0;
    float _S1237 = - ((_S1221 + _S1228 + _S1233) / _S1208);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1238;
    (&_S1238)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1238)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239;
    (&_S1239)->primal_0 = _S1205;
    (&_S1239)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1238, &_S1239, _S1237);
    float3  _S1240 = _S1229.differential_0 + _S1236;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1241;
    (&_S1241)->primal_0 = rov0_4;
    (&_S1241)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1242;
    (&_S1242)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1242)->differential_0 = _S1222;
    s_bwd_prop_cross_0(&_S1241, &_S1242, _S1240);
    float3  _S1243 = _S1225 + _S1239.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1244;
    (&_S1244)->primal_0 = v1v0_4;
    (&_S1244)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1245;
    (&_S1245)->primal_0 = v2v0_4;
    (&_S1245)->differential_0 = _S1222;
    s_bwd_prop_cross_0(&_S1244, &_S1245, _S1243);
    float3  _S1246 = _S1224.differential_0 + _S1241.differential_0;
    float3  _S1247 = _S1235.differential_0 + _S1245.differential_0;
    float3  _S1248 = _S1230.differential_0 + _S1244.differential_0;
    float3  _S1249 = - _S1246 + - _S1247 + - _S1248;
    float3  _S1250 = _S1238.differential_0 + _S1242.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1250;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1246;
    FixedArray<float3 , 3>  _S1251;
    _S1251[int(0)] = _S1222;
    _S1251[int(1)] = _S1222;
    _S1251[int(2)] = _S1222;
    _S1251[int(2)] = _S1214;
    _S1251[int(1)] = _S1216;
    _S1251[int(0)] = _S1218;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1251;
    FixedArray<float3 , 3>  _S1252;
    _S1252[int(0)] = _S1222;
    _S1252[int(1)] = _S1222;
    _S1252[int(2)] = _S1222;
    _S1252[int(2)] = _S1247;
    _S1252[int(0)] = _S1249;
    _S1252[int(1)] = _S1248;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1252;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1253, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1254, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1255, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1256, float3  _S1257, float _S1258)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1253, _S1254, _S1255, _S1256, _S1257, _S1258);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1259 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1260 = { _S1259, _S1259, _S1259 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1260;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1260;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1259;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1259;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

