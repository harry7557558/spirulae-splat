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
        float _S196 = _S195.z;
        *depth_1 = (F32_max((_S196), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S196 + length_1(_S195)))));
        float3  _S197 = mean_1 - - mul_0(transpose_1(R_1), t_1);
        float _S198 = _S197.x;
        float _S199 = _S197.y;
        float _S200 = _S197.z;
        float norm_1 = (F32_sqrt((_S198 * _S198 + _S199 * _S199 + _S200 * _S200)));
        float x_11 = _S198 / norm_1;
        float y_4 = _S199 / norm_1;
        float z_1 = _S200 / norm_1;
        float z2_3 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_11 * x_11 - y_4 * y_4;
        float fS1_1 = 2.0f * x_11 * y_4;
        float fTmp0C_1 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        float3  color_1 = make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * sh_coeffs_1[int(1)] + make_float3 (z_1) * sh_coeffs_1[int(2)] - make_float3 (x_11) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_4) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_11) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_11 * fS1_1 + y_4 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_4) * sh_coeffs_1[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_11) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_11 * fC1_1 - y_4 * fS1_1)) * sh_coeffs_1[int(15)]);
        float3  _S201 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_1 + ch_coeffs_1[int(0)] + make_float3 (0.5f), _S201);
        float3  _S202 = color_1 - ch_coeffs_1[int(0)] * make_float3 (0.5f);
        float3  _S203 = ch_coeffs_1[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S202 + _S203 + make_float3 (0.5f), _S201);
        (*rgbs_1)[int(2)] = max_0(_S202 - _S203 + make_float3 (0.5f), _S201);
        (*verts_1)[int(0)] = vert0_1;
        (*verts_1)[int(1)] = vert1_1;
        (*verts_1)[int(2)] = vert2_1;
        float3  _S204 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S204 * make_float3 (float(- (F32_sign((dot_0(_S204, mean_c_1))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, float4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
{
    float3  mean_c_2 = mul_0(R_2, mean_2) + t_2;
    float _S205 = scale_2.x;
    float sx_2 = (F32_exp((_S205)));
    float _S206 = scale_2.y;
    float sy_2 = (F32_exp((_S206)));
    float sz_2 = scale_2.z - 0.5f * (_S205 + _S206);
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
    Matrix<float, 3, 3>  _S207 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_4), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_4), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)));
    float3  vert0_2 = mul_0(_S207, make_float3 (sx_2, 0.0f, 0.0f)) + mean_2;
    float3  vert1_2 = mul_0(_S207, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_2;
    float3  vert2_2 = mul_0(_S207, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_2;
    float3  vert0_c_2 = mul_0(R_2, vert0_2) + t_2;
    float3  vert1_c_2 = mul_0(R_2, vert1_2) + t_2;
    float3  vert2_c_2 = mul_0(R_2, vert2_2) + t_2;
    float2  _S208 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_9 = _S208.x;
    float v_9 = _S208.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S209 = 2.0f * dist_coeffs_2[int(4)];
    float _S210 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S211 = _S208 * make_float2 (1.0f + r2_9 * (dist_coeffs_2[int(0)] + r2_9 * (dist_coeffs_2[int(1)] + r2_9 * (dist_coeffs_2[int(2)] + r2_9 * dist_coeffs_2[int(3)])))) + make_float2 (_S209 * u_9 * v_9 + dist_coeffs_2[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_2[int(6)] * r2_9, _S210 * u_9 * v_9 + dist_coeffs_2[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_2[int(7)] * r2_9);
    float2  _S212 = _S211 + make_float2 (dist_coeffs_2[int(8)] * _S211.x + dist_coeffs_2[int(9)] * _S211.y, 0.0f);
    float _S213 = fx_2 * _S212.x + cx_2;
    float _S214 = fy_2 * _S212.y + cy_2;
    float2  uv0_4 = make_float2 (_S213, _S214);
    float2  _S215 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_10 = _S215.x;
    float v_10 = _S215.y;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float2  _S216 = _S215 * make_float2 (1.0f + r2_10 * (dist_coeffs_2[int(0)] + r2_10 * (dist_coeffs_2[int(1)] + r2_10 * (dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)])))) + make_float2 (_S209 * u_10 * v_10 + dist_coeffs_2[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + dist_coeffs_2[int(6)] * r2_10, _S210 * u_10 * v_10 + dist_coeffs_2[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + dist_coeffs_2[int(7)] * r2_10);
    float2  _S217 = _S216 + make_float2 (dist_coeffs_2[int(8)] * _S216.x + dist_coeffs_2[int(9)] * _S216.y, 0.0f);
    float _S218 = fx_2 * _S217.x + cx_2;
    float _S219 = fy_2 * _S217.y + cy_2;
    float2  uv1_4 = make_float2 (_S218, _S219);
    float2  _S220 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_11 = _S220.x;
    float v_11 = _S220.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float2  _S221 = _S220 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S209 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S210 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
    float2  _S222 = _S221 + make_float2 (dist_coeffs_2[int(8)] * _S221.x + dist_coeffs_2[int(9)] * _S221.y, 0.0f);
    float _S223 = fx_2 * _S222.x + cx_2;
    float _S224 = fy_2 * _S222.y + cy_2;
    float2  uv2_4 = make_float2 (_S223, _S224);
    float2  e0_2 = uv1_4 - uv0_4;
    float2  e1_2 = uv2_4 - uv1_4;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(uv0_4 - uv2_4)));
    *aabb_xyxy_2 = make_float4 ((F32_min(((F32_min((_S213), (_S218)))), (_S223))) - offset_2, (F32_min(((F32_min((_S214), (_S219)))), (_S224))) - offset_2, (F32_max(((F32_max((_S213), (_S218)))), (_S223))) + offset_2, (F32_max(((F32_max((_S214), (_S219)))), (_S224))) + offset_2);
    *depth_2 = ((vert0_c_2 + vert1_c_2 + vert2_c_2) / make_float3 (3.0f)).z;
    float3  _S225 = mean_2 - - mul_0(transpose_1(R_2), t_2);
    float _S226 = _S225.x;
    float _S227 = _S225.y;
    float _S228 = _S225.z;
    float norm_2 = (F32_sqrt((_S226 * _S226 + _S227 * _S227 + _S228 * _S228)));
    float x_13 = _S226 / norm_2;
    float y_5 = _S227 / norm_2;
    float z_2 = _S228 / norm_2;
    float z2_5 = z_2 * z_2;
    float fTmp0B_2 = -1.09254848957061768f * z_2;
    float fC1_2 = x_13 * x_13 - y_5 * y_5;
    float fS1_2 = 2.0f * x_13 * y_5;
    float fTmp0C_2 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_2;
    float3  color_2 = make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_2[int(1)] + make_float3 (z_2) * sh_coeffs_2[int(2)] - make_float3 (x_13) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_5) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_13) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_13 * fS1_2 + y_5 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_5) * sh_coeffs_2[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_13) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_13 * fC1_2 - y_5 * fS1_2)) * sh_coeffs_2[int(15)]);
    float3  _S229 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_2 + ch_coeffs_2[int(0)] + make_float3 (0.5f), _S229);
    float3  _S230 = color_2 - ch_coeffs_2[int(0)] * make_float3 (0.5f);
    float3  _S231 = ch_coeffs_2[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S230 + _S231 + make_float3 (0.5f), _S229);
    (*rgbs_2)[int(2)] = max_0(_S230 - _S231 + make_float3 (0.5f), _S229);
    (*verts_2)[int(0)] = vert0_2;
    (*verts_2)[int(1)] = vert1_2;
    (*verts_2)[int(2)] = vert2_2;
    float3  _S232 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S232 * make_float3 (float(- (F32_sign((dot_0(_S232, mean_c_2))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, float4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
{
    float3  mean_c_3 = mul_0(R_3, mean_3) + t_3;
    float _S233 = scale_3.x;
    float sx_3 = (F32_exp((_S233)));
    float _S234 = scale_3.y;
    float sy_3 = (F32_exp((_S234)));
    float sz_3 = scale_3.z - 0.5f * (_S233 + _S234);
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
    Matrix<float, 3, 3>  _S235 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_6), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_6), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    float3  vert0_3 = mul_0(_S235, make_float3 (sx_3, 0.0f, 0.0f)) + mean_3;
    float3  vert1_3 = mul_0(_S235, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_3;
    float3  vert2_3 = mul_0(_S235, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_3;
    float3  vert0_c_3 = mul_0(R_3, vert0_3) + t_3;
    float3  vert1_c_3 = mul_0(R_3, vert1_3) + t_3;
    float3  vert2_c_3 = mul_0(R_3, vert2_3) + t_3;
    float2  _S236 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_5 = length_0(_S236);
    float _S237 = vert0_c_3.z;
    float theta_3 = (F32_atan2((r_5), (_S237)));
    float k_1;
    if(theta_3 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_3 * theta_3 / 3.0f) / _S237;
    }
    else
    {
        k_1 = theta_3 / r_5;
    }
    float2  _S238 = _S236 * make_float2 (k_1);
    float u_12 = _S238.x;
    float v_12 = _S238.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S239 = 2.0f * dist_coeffs_3[int(4)];
    float _S240 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S241 = _S238 * make_float2 (1.0f + r2_12 * (dist_coeffs_3[int(0)] + r2_12 * (dist_coeffs_3[int(1)] + r2_12 * (dist_coeffs_3[int(2)] + r2_12 * dist_coeffs_3[int(3)])))) + make_float2 (_S239 * u_12 * v_12 + dist_coeffs_3[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + dist_coeffs_3[int(6)] * r2_12, _S240 * u_12 * v_12 + dist_coeffs_3[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + dist_coeffs_3[int(7)] * r2_12);
    float2  _S242 = _S241 + make_float2 (dist_coeffs_3[int(8)] * _S241.x + dist_coeffs_3[int(9)] * _S241.y, 0.0f);
    float _S243 = fx_3 * _S242.x + cx_3;
    float _S244 = fy_3 * _S242.y + cy_3;
    float2  uv0_5 = make_float2 (_S243, _S244);
    float2  _S245 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_6 = length_0(_S245);
    float _S246 = vert1_c_3.z;
    float theta_4 = (F32_atan2((r_6), (_S246)));
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S246;
    }
    else
    {
        k_1 = theta_4 / r_6;
    }
    float2  _S247 = _S245 * make_float2 (k_1);
    float u_13 = _S247.x;
    float v_13 = _S247.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S248 = _S247 * make_float2 (1.0f + r2_13 * (dist_coeffs_3[int(0)] + r2_13 * (dist_coeffs_3[int(1)] + r2_13 * (dist_coeffs_3[int(2)] + r2_13 * dist_coeffs_3[int(3)])))) + make_float2 (_S239 * u_13 * v_13 + dist_coeffs_3[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_3[int(6)] * r2_13, _S240 * u_13 * v_13 + dist_coeffs_3[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_3[int(7)] * r2_13);
    float2  _S249 = _S248 + make_float2 (dist_coeffs_3[int(8)] * _S248.x + dist_coeffs_3[int(9)] * _S248.y, 0.0f);
    float _S250 = fx_3 * _S249.x + cx_3;
    float _S251 = fy_3 * _S249.y + cy_3;
    float2  uv1_5 = make_float2 (_S250, _S251);
    float2  _S252 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_7 = length_0(_S252);
    float _S253 = vert2_c_3.z;
    float theta_5 = (F32_atan2((r_7), (_S253)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S253;
    }
    else
    {
        k_1 = theta_5 / r_7;
    }
    float2  _S254 = _S252 * make_float2 (k_1);
    float u_14 = _S254.x;
    float v_14 = _S254.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S255 = _S254 * make_float2 (1.0f + r2_14 * (dist_coeffs_3[int(0)] + r2_14 * (dist_coeffs_3[int(1)] + r2_14 * (dist_coeffs_3[int(2)] + r2_14 * dist_coeffs_3[int(3)])))) + make_float2 (_S239 * u_14 * v_14 + dist_coeffs_3[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + dist_coeffs_3[int(6)] * r2_14, _S240 * u_14 * v_14 + dist_coeffs_3[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + dist_coeffs_3[int(7)] * r2_14);
    float2  _S256 = _S255 + make_float2 (dist_coeffs_3[int(8)] * _S255.x + dist_coeffs_3[int(9)] * _S255.y, 0.0f);
    float _S257 = fx_3 * _S256.x + cx_3;
    float _S258 = fy_3 * _S256.y + cy_3;
    float2  uv2_5 = make_float2 (_S257, _S258);
    float2  e0_3 = uv1_5 - uv0_5;
    float2  e1_3 = uv2_5 - uv1_5;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(uv0_5 - uv2_5)));
    *aabb_xyxy_3 = make_float4 ((F32_min(((F32_min((_S243), (_S250)))), (_S257))) - offset_3, (F32_min(((F32_min((_S244), (_S251)))), (_S258))) - offset_3, (F32_max(((F32_max((_S243), (_S250)))), (_S257))) + offset_3, (F32_max(((F32_max((_S244), (_S251)))), (_S258))) + offset_3);
    float3  _S259 = (vert0_c_3 + vert1_c_3 + vert2_c_3) / make_float3 (3.0f);
    float _S260 = _S259.z;
    *depth_3 = (F32_max((_S260), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S260 + length_1(_S259)))));
    float3  _S261 = mean_3 - - mul_0(transpose_1(R_3), t_3);
    float _S262 = _S261.x;
    float _S263 = _S261.y;
    float _S264 = _S261.z;
    float norm_3 = (F32_sqrt((_S262 * _S262 + _S263 * _S263 + _S264 * _S264)));
    float x_15 = _S262 / norm_3;
    float y_6 = _S263 / norm_3;
    float z_3 = _S264 / norm_3;
    float z2_7 = z_3 * z_3;
    float fTmp0B_3 = -1.09254848957061768f * z_3;
    float fC1_3 = x_15 * x_15 - y_6 * y_6;
    float fS1_3 = 2.0f * x_15 * y_6;
    float fTmp0C_3 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_3;
    float3  color_3 = make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_3[int(1)] + make_float3 (z_3) * sh_coeffs_3[int(2)] - make_float3 (x_15) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_6) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_15) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_3 + y_6 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_6) * sh_coeffs_3[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_15) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_3 - y_6 * fS1_3)) * sh_coeffs_3[int(15)]);
    float3  _S265 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_3 + ch_coeffs_3[int(0)] + make_float3 (0.5f), _S265);
    float3  _S266 = color_3 - ch_coeffs_3[int(0)] * make_float3 (0.5f);
    float3  _S267 = ch_coeffs_3[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S266 + _S267 + make_float3 (0.5f), _S265);
    (*rgbs_3)[int(2)] = max_0(_S266 - _S267 + make_float3 (0.5f), _S265);
    (*verts_3)[int(0)] = vert0_3;
    (*verts_3)[int(1)] = vert1_3;
    (*verts_3)[int(2)] = vert2_3;
    float3  _S268 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S268 * make_float3 (float(- (F32_sign((dot_0(_S268, mean_c_3))))));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S269, float3  _S270)
{
    return mul_0(_S269, _S270);
}

inline __device__ float s_primal_ctx_exp_0(float _S271)
{
    return (F32_exp((_S271)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S272)
{
    return (F32_sqrt((_S272)));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S273, float3  _S274)
{
    return cross_0(_S273, _S274);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S275, float3  _S276)
{
    return dot_0(_S275, _S276);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S277, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S278, float _S279)
{
    _d_dot_0(_S277, _S278, _S279);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S280, float _S281)
{
    _d_sqrt_0(_S280, _S281);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S282 = (*dpx_9).primal_0.x;
    float _S283 = (*dpx_9).primal_0.y;
    float _S284 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S285;
    (&_S285)->primal_0 = _S282 * _S282 + _S283 * _S283 + _S284 * _S284;
    (&_S285)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S285, _s_dOut_0);
    float _S286 = (*dpx_9).primal_0.z * _S285.differential_0;
    float _S287 = _S286 + _S286;
    float _S288 = (*dpx_9).primal_0.y * _S285.differential_0;
    float _S289 = _S288 + _S288;
    float _S290 = (*dpx_9).primal_0.x * _S285.differential_0;
    float _S291 = _S290 + _S290;
    float3  _S292 = make_float3 (0.0f);
    *&((&_S292)->z) = _S287;
    *&((&_S292)->y) = _S289;
    *&((&_S292)->x) = _S291;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S292;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S293, float _S294)
{
    s_bwd_prop_length_impl_0(_S293, _S294);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_1)
{
    float _S295 = length_1((*dpx_10).primal_0);
    float3  _S296 = (*dpx_10).primal_0 * _s_dOut_1;
    float3  _S297 = make_float3 (1.0f / _S295) * _s_dOut_1;
    float _S298 = - ((_S296.x + _S296.y + _S296.z) / (_S295 * _S295));
    float3  _S299 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S300;
    (&_S300)->primal_0 = (*dpx_10).primal_0;
    (&_S300)->differential_0 = _S299;
    s_bwd_length_impl_0(&_S300, _S298);
    float3  _S301 = _S297 + _S300.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S301;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S302, float3  _S303)
{
    s_bwd_prop_normalize_impl_0(_S302, _S303);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S304, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S305, float3  _S306)
{
    _d_cross_0(_S304, _S305, _S306);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S307, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S308, float3  _S309)
{
    _d_max_vector_0(_S307, _S308, _S309);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S310, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S311, float3  _S312)
{
    _d_mul_0(_S310, _S311, _S312);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S313, float _S314)
{
    _d_exp2_0(_S313, _S314);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float _s_dOut_2)
{
    float _S315 = (*dpx_11).primal_0.x;
    float _S316 = (*dpx_11).primal_0.y;
    DiffPair_float_0 _S317;
    (&_S317)->primal_0 = _S315 * _S315 + _S316 * _S316;
    (&_S317)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S317, _s_dOut_2);
    float _S318 = (*dpx_11).primal_0.y * _S317.differential_0;
    float _S319 = _S318 + _S318;
    float _S320 = (*dpx_11).primal_0.x * _S317.differential_0;
    float _S321 = _S320 + _S320;
    float2  _S322 = make_float2 (0.0f);
    *&((&_S322)->y) = _S319;
    *&((&_S322)->x) = _S321;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S322;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S323, float _S324)
{
    s_bwd_prop_length_impl_1(_S323, _S324);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S325, float _S326)
{
    _d_abs_0(_S325, _S326);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S327, float _S328)
{
    _d_exp_0(_S327, _S328);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S329 = scale_4.x;
    float _S330 = s_primal_ctx_exp_0(_S329);
    float _S331 = scale_4.y;
    float _S332 = s_primal_ctx_exp_0(_S331);
    float sz_4 = scale_4.z - 0.5f * (_S329 + _S331);
    float _S333 = quat_4.y;
    float x2_4 = _S333 * _S333;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_8 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  _S334 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S335 = make_float3 (_S330, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S334, _S335) + mean_4;
    float _S336 = -0.5f + sz_4;
    float3  _S337 = make_float3 (_S330 * _S336, _S332, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S334, _S337) + mean_4;
    float _S338 = -0.5f - sz_4;
    float3  _S339 = make_float3 (_S330 * _S338, - _S332, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S334, _S339) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S340 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S341 = vert0_c_4.z;
    float2  _S342 = make_float2 (_S341);
    float2  _S343 = _S340 / make_float2 (_S341);
    float2  _S344 = make_float2 (_S341 * _S341);
    float u_15 = _S343.x;
    float v_15 = _S343.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S345 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S346 = dist_coeffs_4[int(1)] + r2_15 * _S345;
    float _S347 = dist_coeffs_4[int(0)] + r2_15 * _S346;
    float radial_6 = 1.0f + r2_15 * _S347;
    float _S348 = 2.0f * dist_coeffs_4[int(4)];
    float _S349 = _S348 * u_15;
    float _S350 = 2.0f * u_15;
    float _S351 = 2.0f * dist_coeffs_4[int(5)];
    float _S352 = _S351 * u_15;
    float _S353 = 2.0f * v_15;
    float2  _S354 = _S343 * make_float2 (radial_6) + make_float2 (_S349 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S350 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S352 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S353 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S355 = _S354 + make_float2 (dist_coeffs_4[int(8)] * _S354.x + dist_coeffs_4[int(9)] * _S354.y, 0.0f);
    float _S356 = fx_4 * _S355.x + cx_4;
    float _S357 = fy_4 * _S355.y + cy_4;
    float2  uv0_6 = make_float2 (_S356, _S357);
    float2  _S358 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S359 = vert1_c_4.z;
    float2  _S360 = make_float2 (_S359);
    float2  _S361 = _S358 / make_float2 (_S359);
    float2  _S362 = make_float2 (_S359 * _S359);
    float u_16 = _S361.x;
    float v_16 = _S361.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S363 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S364 = dist_coeffs_4[int(1)] + r2_16 * _S363;
    float _S365 = dist_coeffs_4[int(0)] + r2_16 * _S364;
    float radial_7 = 1.0f + r2_16 * _S365;
    float _S366 = _S348 * u_16;
    float _S367 = 2.0f * u_16;
    float _S368 = _S351 * u_16;
    float _S369 = 2.0f * v_16;
    float2  _S370 = _S361 * make_float2 (radial_7) + make_float2 (_S366 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S367 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S368 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S369 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S371 = _S370 + make_float2 (dist_coeffs_4[int(8)] * _S370.x + dist_coeffs_4[int(9)] * _S370.y, 0.0f);
    float _S372 = fx_4 * _S371.x + cx_4;
    float _S373 = fy_4 * _S371.y + cy_4;
    float2  uv1_6 = make_float2 (_S372, _S373);
    float2  _S374 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S375 = vert2_c_4.z;
    float2  _S376 = make_float2 (_S375);
    float2  _S377 = _S374 / make_float2 (_S375);
    float2  _S378 = make_float2 (_S375 * _S375);
    float u_17 = _S377.x;
    float v_17 = _S377.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S379 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S380 = dist_coeffs_4[int(1)] + r2_17 * _S379;
    float _S381 = dist_coeffs_4[int(0)] + r2_17 * _S380;
    float radial_8 = 1.0f + r2_17 * _S381;
    float _S382 = _S348 * u_17;
    float _S383 = 2.0f * u_17;
    float _S384 = _S351 * u_17;
    float _S385 = 2.0f * v_17;
    float2  _S386 = _S377 * make_float2 (radial_8) + make_float2 (_S382 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S383 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S384 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S385 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S387 = _S386 + make_float2 (dist_coeffs_4[int(8)] * _S386.x + dist_coeffs_4[int(9)] * _S386.y, 0.0f);
    float _S388 = fx_4 * _S387.x + cx_4;
    float _S389 = fy_4 * _S387.y + cy_4;
    float2  uv2_6 = make_float2 (_S388, _S389);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S390 = e0_4.x;
    float _S391 = e1_4.y;
    float _S392 = e0_4.y;
    float _S393 = e1_4.x;
    float _S394 = _S390 * _S391 - _S392 * _S393;
    float _S395 = 1.0f - hardness_4.y;
    float _S396 = -1.0f / _S395;
    float _S397 = _S395 * _S395;
    float _S398 = (F32_max((_S356), (_S372)));
    float _S399 = (F32_min((_S356), (_S372)));
    float _S400 = (F32_max((_S357), (_S373)));
    float _S401 = (F32_min((_S357), (_S373)));
    Matrix<float, 3, 3>  _S402 = transpose_1(R_4);
    float3  _S403 = mean_4 - - s_primal_ctx_mul_0(_S402, t_4);
    float _S404 = _S403.x;
    float _S405 = _S403.y;
    float _S406 = _S403.z;
    float _S407 = _S404 * _S404 + _S405 * _S405 + _S406 * _S406;
    float _S408 = s_primal_ctx_sqrt_0(_S407);
    float x_16 = _S404 / _S408;
    float3  _S409 = make_float3 (x_16);
    float _S410 = _S408 * _S408;
    float y_7 = _S405 / _S408;
    float z_4 = _S406 / _S408;
    float3  _S411 = make_float3 (z_4);
    float _S412 = - y_7;
    float3  _S413 = make_float3 (_S412);
    float z2_9 = z_4 * z_4;
    float fTmp0B_4 = -1.09254848957061768f * z_4;
    float fC1_4 = x_16 * x_16 - y_7 * y_7;
    float _S414 = 2.0f * x_16;
    float fS1_4 = _S414 * y_7;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S415 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_16;
    float3  _S416 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_7;
    float3  _S417 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S418 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S419 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_4;
    float _S420 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_4 * _S420;
    float3  _S421 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_16;
    float3  _S422 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_7;
    float3  _S423 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S424 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S425 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_16 * fC1_4 - y_7 * fS1_4);
    float3  _S426 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_16 * fS1_4 + y_7 * fC1_4);
    float3  _S427 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S412) * sh_coeffs_4[int(1)] + make_float3 (z_4) * sh_coeffs_4[int(2)] - make_float3 (x_16) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]);
    float3  _S428 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S429 = make_float3 (0.0f);
    float3  _S430 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S431 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S432 = make_float3 (_S431);
    float3  _S433 = ch_coeffs_4[int(1)] * make_float3 (_S431);
    float3  _S434 = _S430 + _S433 + make_float3 (0.5f);
    float3  _S435 = _S430 - _S433 + make_float3 (0.5f);
    float3  _S436 = vert1_c_4 - vert0_c_4;
    float3  _S437 = vert2_c_4 - vert0_c_4;
    float3  _S438 = s_primal_ctx_cross_0(_S436, _S437);
    float3  _S439 = normalize_0(_S438);
    float3  _S440 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S439, mean_c_4)))))) * v_normal_0;
    float3  _S441 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S442;
    (&_S442)->primal_0 = _S439;
    (&_S442)->differential_0 = _S441;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S443;
    (&_S443)->primal_0 = mean_c_4;
    (&_S443)->differential_0 = _S441;
    s_bwd_prop_dot_0(&_S442, &_S443, 0.0f);
    float3  _S444 = _S440 + _S442.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S445;
    (&_S445)->primal_0 = _S438;
    (&_S445)->differential_0 = _S441;
    s_bwd_normalize_impl_0(&_S445, _S444);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S446;
    (&_S446)->primal_0 = _S436;
    (&_S446)->differential_0 = _S441;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S447;
    (&_S447)->primal_0 = _S437;
    (&_S447)->differential_0 = _S441;
    s_bwd_prop_cross_0(&_S446, &_S447, _S445.differential_0);
    float3  _S448 = - _S447.differential_0;
    float3  _S449 = - _S446.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S450;
    (&_S450)->primal_0 = _S435;
    (&_S450)->differential_0 = _S441;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S451;
    (&_S451)->primal_0 = _S429;
    (&_S451)->differential_0 = _S441;
    s_bwd_prop_max_0(&_S450, &_S451, v_rgbs_0[int(2)]);
    float3  _S452 = - _S450.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S453;
    (&_S453)->primal_0 = _S434;
    (&_S453)->differential_0 = _S441;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S454;
    (&_S454)->primal_0 = _S429;
    (&_S454)->differential_0 = _S441;
    s_bwd_prop_max_0(&_S453, &_S454, v_rgbs_0[int(1)]);
    float3  _S455 = _S432 * (_S452 + _S453.differential_0);
    float3  _S456 = _S450.differential_0 + _S453.differential_0;
    float3  _S457 = make_float3 (0.5f) * - _S456;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S458;
    (&_S458)->primal_0 = _S428;
    (&_S458)->differential_0 = _S441;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S459;
    (&_S459)->primal_0 = _S429;
    (&_S459)->differential_0 = _S441;
    s_bwd_prop_max_0(&_S458, &_S459, v_rgbs_0[int(0)]);
    float3  _S460 = _S457 + _S458.differential_0;
    float3  _S461 = _S456 + _S458.differential_0;
    float3  _S462 = _S426 * _S461;
    float3  _S463 = sh_coeffs_4[int(15)] * _S461;
    float3  _S464 = _S424 * _S461;
    float3  _S465 = sh_coeffs_4[int(14)] * _S461;
    float3  _S466 = _S422 * _S461;
    float3  _S467 = sh_coeffs_4[int(13)] * _S461;
    float3  _S468 = _S421 * _S461;
    float3  _S469 = sh_coeffs_4[int(12)] * _S461;
    float3  _S470 = _S423 * _S461;
    float3  _S471 = sh_coeffs_4[int(11)] * _S461;
    float3  _S472 = _S425 * _S461;
    float3  _S473 = sh_coeffs_4[int(10)] * _S461;
    float3  _S474 = _S427 * _S461;
    float3  _S475 = sh_coeffs_4[int(9)] * _S461;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S475.x + _S475.y + _S475.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S463.x + _S463.y + _S463.z);
    float _S476 = _S473.x + _S473.y + _S473.z;
    float _S477 = _S465.x + _S465.y + _S465.z;
    float _S478 = _S471.x + _S471.y + _S471.z;
    float _S479 = _S467.x + _S467.y + _S467.z;
    float _S480 = _S469.x + _S469.y + _S469.z;
    float _S481 = - s_diff_fC2_T_0;
    float3  _S482 = _S418 * _S461;
    float3  _S483 = sh_coeffs_4[int(8)] * _S461;
    float3  _S484 = _S416 * _S461;
    float3  _S485 = sh_coeffs_4[int(7)] * _S461;
    float3  _S486 = _S415 * _S461;
    float3  _S487 = sh_coeffs_4[int(6)] * _S461;
    float3  _S488 = _S417 * _S461;
    float3  _S489 = sh_coeffs_4[int(5)] * _S461;
    float3  _S490 = _S419 * _S461;
    float3  _S491 = sh_coeffs_4[int(4)] * _S461;
    float _S492 = _S489.x + _S489.y + _S489.z;
    float _S493 = _S485.x + _S485.y + _S485.z;
    float _S494 = fTmp1B_4 * _S476 + x_16 * s_diff_fS2_T_0 + y_7 * _S481 + 0.54627424478530884f * (_S491.x + _S491.y + _S491.z);
    float _S495 = fTmp1B_4 * _S477 + y_7 * s_diff_fS2_T_0 + x_16 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S483.x + _S483.y + _S483.z);
    float _S496 = y_7 * - _S495;
    float _S497 = x_16 * _S495;
    float _S498 = z_4 * (1.86588168144226074f * (z_4 * _S480) + -2.28522896766662598f * (y_7 * _S478 + x_16 * _S479) + 0.94617468118667603f * (_S487.x + _S487.y + _S487.z));
    float3  _S499 = make_float3 (0.48860251903533936f) * _S461;
    float3  _S500 = - _S499;
    float3  _S501 = _S409 * _S500;
    float3  _S502 = sh_coeffs_4[int(3)] * _S500;
    float3  _S503 = _S411 * _S499;
    float3  _S504 = sh_coeffs_4[int(2)] * _S499;
    float3  _S505 = _S413 * _S499;
    float3  _S506 = sh_coeffs_4[int(1)] * _S499;
    float _S507 = (_S420 * _S480 + 1.44530570507049561f * (fS1_4 * _S476 + fC1_4 * _S477) + -1.09254848957061768f * (y_7 * _S492 + x_16 * _S493) + _S498 + _S498 + _S504.x + _S504.y + _S504.z) / _S410;
    float _S508 = _S408 * _S507;
    float _S509 = (fTmp0C_4 * _S478 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S481 + fTmp0B_4 * _S492 + _S414 * _S494 + _S496 + _S496 + - (_S506.x + _S506.y + _S506.z)) / _S410;
    float _S510 = _S408 * _S509;
    float _S511 = (fTmp0C_4 * _S479 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S493 + 2.0f * (y_7 * _S494) + _S497 + _S497 + _S502.x + _S502.y + _S502.z) / _S410;
    float _S512 = _S408 * _S511;
    float _S513 = _S406 * - _S507 + _S405 * - _S509 + _S404 * - _S511;
    DiffPair_float_0 _S514;
    (&_S514)->primal_0 = _S407;
    (&_S514)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S514, _S513);
    float _S515 = _S406 * _S514.differential_0;
    float _S516 = _S405 * _S514.differential_0;
    float _S517 = _S404 * _S514.differential_0;
    float3  _S518 = make_float3 (0.282094806432724f) * _S461;
    float3  _S519 = make_float3 (_S512 + _S517 + _S517, _S510 + _S516 + _S516, _S508 + _S515 + _S515);
    float3  _S520 = - - _S519;
    Matrix<float, 3, 3>  _S521 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S522;
    (&_S522)->primal_0 = _S402;
    (&_S522)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S523;
    (&_S523)->primal_0 = t_4;
    (&_S523)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S522, &_S523, _S520);
    Matrix<float, 3, 3>  _S524 = transpose_1(_S522.differential_0);
    float3  _S525 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_0);
    DiffPair_float_0 _S526;
    (&_S526)->primal_0 = _S401;
    (&_S526)->differential_0 = 0.0f;
    DiffPair_float_0 _S527;
    (&_S527)->primal_0 = _S389;
    (&_S527)->differential_0 = 0.0f;
    _d_min_0(&_S526, &_S527, 0.0f);
    DiffPair_float_0 _S528;
    (&_S528)->primal_0 = _S357;
    (&_S528)->differential_0 = 0.0f;
    DiffPair_float_0 _S529;
    (&_S529)->primal_0 = _S373;
    (&_S529)->differential_0 = 0.0f;
    _d_min_0(&_S528, &_S529, _S526.differential_0);
    DiffPair_float_0 _S530;
    (&_S530)->primal_0 = _S400;
    (&_S530)->differential_0 = 0.0f;
    DiffPair_float_0 _S531;
    (&_S531)->primal_0 = _S389;
    (&_S531)->differential_0 = 0.0f;
    _d_max_0(&_S530, &_S531, 0.0f);
    DiffPair_float_0 _S532;
    (&_S532)->primal_0 = _S357;
    (&_S532)->differential_0 = 0.0f;
    DiffPair_float_0 _S533;
    (&_S533)->primal_0 = _S373;
    (&_S533)->differential_0 = 0.0f;
    _d_max_0(&_S532, &_S533, _S530.differential_0);
    DiffPair_float_0 _S534;
    (&_S534)->primal_0 = _S399;
    (&_S534)->differential_0 = 0.0f;
    DiffPair_float_0 _S535;
    (&_S535)->primal_0 = _S388;
    (&_S535)->differential_0 = 0.0f;
    _d_min_0(&_S534, &_S535, 0.0f);
    DiffPair_float_0 _S536;
    (&_S536)->primal_0 = _S356;
    (&_S536)->differential_0 = 0.0f;
    DiffPair_float_0 _S537;
    (&_S537)->primal_0 = _S372;
    (&_S537)->differential_0 = 0.0f;
    _d_min_0(&_S536, &_S537, _S534.differential_0);
    DiffPair_float_0 _S538;
    (&_S538)->primal_0 = _S398;
    (&_S538)->differential_0 = 0.0f;
    DiffPair_float_0 _S539;
    (&_S539)->primal_0 = _S388;
    (&_S539)->differential_0 = 0.0f;
    _d_max_0(&_S538, &_S539, 0.0f);
    DiffPair_float_0 _S540;
    (&_S540)->primal_0 = _S356;
    (&_S540)->differential_0 = 0.0f;
    DiffPair_float_0 _S541;
    (&_S541)->primal_0 = _S372;
    (&_S541)->differential_0 = 0.0f;
    _d_max_0(&_S540, &_S541, _S538.differential_0);
    DiffPair_float_0 _S542;
    (&_S542)->primal_0 = _S396;
    (&_S542)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S542, -0.0f);
    float _S543 = - (-1.0f * - (_S542.differential_0 / _S397));
    float2  _S544 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S545;
    (&_S545)->primal_0 = e2_0;
    (&_S545)->differential_0 = _S544;
    s_bwd_length_impl_1(&_S545, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S546;
    (&_S546)->primal_0 = e1_4;
    (&_S546)->differential_0 = _S544;
    s_bwd_length_impl_1(&_S546, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S547;
    (&_S547)->primal_0 = e0_4;
    (&_S547)->differential_0 = _S544;
    s_bwd_length_impl_1(&_S547, 0.0f);
    DiffPair_float_0 _S548;
    (&_S548)->primal_0 = _S394;
    (&_S548)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S548, -0.0f);
    float _S549 = - _S548.differential_0;
    float2  _S550 = _S546.differential_0 + make_float2 (_S392 * _S549, _S390 * _S548.differential_0);
    float2  _S551 = _S547.differential_0 + make_float2 (_S391 * _S548.differential_0, _S393 * _S549);
    float2  _S552 = - _S545.differential_0 + _S550;
    float _S553 = fx_4 * (_S535.differential_0 + _S539.differential_0 + _S552.x);
    float2  _S554 = make_float2 (_S553, fy_4 * (_S527.differential_0 + _S531.differential_0 + _S552.y)) + make_float2 (dist_coeffs_4[int(8)] * _S553, dist_coeffs_4[int(9)] * _S553);
    float2  _S555 = _S377 * _S554;
    float _S556 = dist_coeffs_4[int(4)] * _S554.y;
    float _S557 = dist_coeffs_4[int(5)] * _S554.x;
    float _S558 = _S555.x + _S555.y;
    float _S559 = r2_17 * _S558;
    float _S560 = r2_17 * _S559;
    float _S561 = dist_coeffs_4[int(7)] * _S554.y + _S556 + dist_coeffs_4[int(6)] * _S554.x + _S557 + _S381 * _S558 + _S380 * _S559 + _S379 * _S560 + dist_coeffs_4[int(3)] * (r2_17 * _S560);
    float _S562 = v_17 * _S561;
    float _S563 = u_17 * _S561;
    float2  _S564 = (make_float2 (radial_8) * _S554 + make_float2 (_S351 * (v_17 * _S554.y) + _S383 * _S557 + 2.0f * (u_17 * _S557) + _S348 * (v_17 * _S554.x) + _S563 + _S563, _S385 * _S556 + 2.0f * (v_17 * _S556) + _S384 * _S554.y + _S382 * _S554.x + _S562 + _S562)) / _S378;
    float2  _S565 = _S374 * - _S564;
    float2  _S566 = _S376 * _S564;
    float2  _S567 = - _S550 + _S551;
    float _S568 = fx_4 * (_S537.differential_0 + _S541.differential_0 + _S567.x);
    float2  _S569 = make_float2 (_S568, fy_4 * (_S529.differential_0 + _S533.differential_0 + _S567.y)) + make_float2 (dist_coeffs_4[int(8)] * _S568, dist_coeffs_4[int(9)] * _S568);
    float2  _S570 = _S361 * _S569;
    float _S571 = dist_coeffs_4[int(4)] * _S569.y;
    float _S572 = dist_coeffs_4[int(5)] * _S569.x;
    float _S573 = _S570.x + _S570.y;
    float _S574 = r2_16 * _S573;
    float _S575 = r2_16 * _S574;
    float _S576 = dist_coeffs_4[int(7)] * _S569.y + _S571 + dist_coeffs_4[int(6)] * _S569.x + _S572 + _S365 * _S573 + _S364 * _S574 + _S363 * _S575 + dist_coeffs_4[int(3)] * (r2_16 * _S575);
    float _S577 = v_16 * _S576;
    float _S578 = u_16 * _S576;
    float2  _S579 = (make_float2 (radial_7) * _S569 + make_float2 (_S351 * (v_16 * _S569.y) + _S367 * _S572 + 2.0f * (u_16 * _S572) + _S348 * (v_16 * _S569.x) + _S578 + _S578, _S369 * _S571 + 2.0f * (v_16 * _S571) + _S368 * _S569.y + _S366 * _S569.x + _S577 + _S577)) / _S362;
    float2  _S580 = _S358 * - _S579;
    float2  _S581 = _S360 * _S579;
    float _S582 = _S580.x + _S580.y;
    float2  _S583 = _S545.differential_0 + - _S551;
    float _S584 = fx_4 * (_S536.differential_0 + _S540.differential_0 + _S583.x);
    float2  _S585 = make_float2 (_S584, fy_4 * (_S528.differential_0 + _S532.differential_0 + _S583.y)) + make_float2 (dist_coeffs_4[int(8)] * _S584, dist_coeffs_4[int(9)] * _S584);
    float2  _S586 = _S343 * _S585;
    float _S587 = dist_coeffs_4[int(4)] * _S585.y;
    float _S588 = dist_coeffs_4[int(5)] * _S585.x;
    float _S589 = _S586.x + _S586.y;
    float _S590 = r2_15 * _S589;
    float _S591 = r2_15 * _S590;
    float _S592 = dist_coeffs_4[int(7)] * _S585.y + _S587 + dist_coeffs_4[int(6)] * _S585.x + _S588 + _S347 * _S589 + _S346 * _S590 + _S345 * _S591 + dist_coeffs_4[int(3)] * (r2_15 * _S591);
    float _S593 = v_15 * _S592;
    float _S594 = u_15 * _S592;
    float2  _S595 = (make_float2 (radial_6) * _S585 + make_float2 (_S351 * (v_15 * _S585.y) + _S350 * _S588 + 2.0f * (u_15 * _S588) + _S348 * (v_15 * _S585.x) + _S594 + _S594, _S353 * _S587 + 2.0f * (v_15 * _S587) + _S352 * _S585.y + _S349 * _S585.x + _S593 + _S593)) / _S344;
    float2  _S596 = _S340 * - _S595;
    float2  _S597 = _S342 * _S595;
    float _S598 = _S596.x + _S596.y;
    float3  _S599 = _S447.differential_0 + _S525 + make_float3 (_S566.x, _S566.y, _S565.x + _S565.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S600;
    (&_S600)->primal_0 = R_4;
    (&_S600)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S601;
    (&_S601)->primal_0 = vert2_4;
    (&_S601)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S600, &_S601, _S599);
    float3  _S602 = _S446.differential_0 + _S525 + make_float3 (_S581.x, _S581.y, _S582);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S603;
    (&_S603)->primal_0 = R_4;
    (&_S603)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S604;
    (&_S604)->primal_0 = vert1_4;
    (&_S604)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S603, &_S604, _S602);
    float3  _S605 = _S448 + _S449 + _S525 + make_float3 (_S597.x, _S597.y, _S598);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S606;
    (&_S606)->primal_0 = R_4;
    (&_S606)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S607;
    (&_S607)->primal_0 = vert0_4;
    (&_S607)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S606, &_S607, _S605);
    float3  _S608 = v_verts_0[int(2)] + _S601.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S609;
    (&_S609)->primal_0 = _S334;
    (&_S609)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S610;
    (&_S610)->primal_0 = _S339;
    (&_S610)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S609, &_S610, _S608);
    float _S611 = - _S610.differential_0.y;
    float _S612 = _S338 * _S610.differential_0.x;
    float _S613 = - (_S330 * _S610.differential_0.x);
    float3  _S614 = v_verts_0[int(1)] + _S604.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S615;
    (&_S615)->primal_0 = _S334;
    (&_S615)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S616;
    (&_S616)->primal_0 = _S337;
    (&_S616)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S615, &_S616, _S614);
    float _S617 = _S330 * _S616.differential_0.x;
    float _S618 = _S336 * _S616.differential_0.x;
    float3  _S619 = v_verts_0[int(0)] + _S607.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S620;
    (&_S620)->primal_0 = _S334;
    (&_S620)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S621;
    (&_S621)->primal_0 = _S335;
    (&_S621)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S620, &_S621, _S619);
    Matrix<float, 3, 3>  _S622 = transpose_1(_S609.differential_0 + _S615.differential_0 + _S620.differential_0);
    float _S623 = 2.0f * - _S622.rows[int(2)].z;
    float _S624 = 2.0f * _S622.rows[int(2)].y;
    float _S625 = 2.0f * _S622.rows[int(2)].x;
    float _S626 = 2.0f * _S622.rows[int(1)].z;
    float _S627 = 2.0f * - _S622.rows[int(1)].y;
    float _S628 = 2.0f * _S622.rows[int(1)].x;
    float _S629 = 2.0f * _S622.rows[int(0)].z;
    float _S630 = 2.0f * _S622.rows[int(0)].y;
    float _S631 = 2.0f * - _S622.rows[int(0)].x;
    float _S632 = - _S628 + _S630;
    float _S633 = _S625 + - _S629;
    float _S634 = - _S624 + _S626;
    float _S635 = _S624 + _S626;
    float _S636 = _S625 + _S629;
    float _S637 = _S628 + _S630;
    float _S638 = quat_4.w * (_S627 + _S631);
    float _S639 = quat_4.z * (_S623 + _S631);
    float _S640 = quat_4.y * (_S623 + _S627);
    float _S641 = quat_4.x * _S632 + quat_4.z * _S635 + quat_4.y * _S636 + _S638 + _S638;
    float _S642 = quat_4.x * _S633 + quat_4.w * _S635 + quat_4.y * _S637 + _S639 + _S639;
    float _S643 = quat_4.x * _S634 + quat_4.w * _S636 + quat_4.z * _S637 + _S640 + _S640;
    float _S644 = quat_4.w * _S632 + quat_4.z * _S633 + quat_4.y * _S634;
    float _S645 = _S613 + _S617;
    float _S646 = 0.5f * - _S645;
    float _S647 = _S611 + _S616.differential_0.y;
    DiffPair_float_0 _S648;
    (&_S648)->primal_0 = _S331;
    (&_S648)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S648, _S647);
    float _S649 = _S646 + _S648.differential_0;
    float _S650 = _S612 + _S618 + _S621.differential_0.x;
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = _S329;
    (&_S651)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S651, _S650);
    float _S652 = _S646 + _S651.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S653;
    (&_S653)->primal_0 = R_4;
    (&_S653)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S654;
    (&_S654)->primal_0 = mean_4;
    (&_S654)->differential_0 = _S441;
    s_bwd_prop_mul_0(&_S653, &_S654, _S443.differential_0);
    float3  _S655 = _S523.differential_0 + _S599 + _S602 + _S605 + _S443.differential_0;
    Matrix<float, 3, 3>  _S656 = _S524 + _S600.differential_0 + _S603.differential_0 + _S606.differential_0 + _S653.differential_0;
    FixedArray<float3 , 2>  _S657;
    _S657[int(0)] = _S441;
    _S657[int(1)] = _S441;
    _S657[int(1)] = _S455;
    _S657[int(0)] = _S460;
    FixedArray<float3 , 16>  _S658;
    _S658[int(0)] = _S441;
    _S658[int(1)] = _S441;
    _S658[int(2)] = _S441;
    _S658[int(3)] = _S441;
    _S658[int(4)] = _S441;
    _S658[int(5)] = _S441;
    _S658[int(6)] = _S441;
    _S658[int(7)] = _S441;
    _S658[int(8)] = _S441;
    _S658[int(9)] = _S441;
    _S658[int(10)] = _S441;
    _S658[int(11)] = _S441;
    _S658[int(12)] = _S441;
    _S658[int(13)] = _S441;
    _S658[int(14)] = _S441;
    _S658[int(15)] = _S441;
    _S658[int(15)] = _S462;
    _S658[int(14)] = _S464;
    _S658[int(13)] = _S466;
    _S658[int(12)] = _S468;
    _S658[int(11)] = _S470;
    _S658[int(10)] = _S472;
    _S658[int(9)] = _S474;
    _S658[int(8)] = _S482;
    _S658[int(7)] = _S484;
    _S658[int(6)] = _S486;
    _S658[int(5)] = _S488;
    _S658[int(4)] = _S490;
    _S658[int(3)] = _S501;
    _S658[int(2)] = _S503;
    _S658[int(1)] = _S505;
    _S658[int(0)] = _S518;
    float2  _S659 = make_float2 (0.0f, _S543);
    float3  _S660 = make_float3 (_S652, _S649, _S645);
    float4  _S661 = make_float4 (0.0f);
    *&((&_S661)->w) = _S641;
    *&((&_S661)->z) = _S642;
    *&((&_S661)->y) = _S643;
    *&((&_S661)->x) = _S644;
    *v_mean_0 = _S519 + _S608 + _S614 + _S619 + _S654.differential_0;
    *v_quat_0 = _S661;
    *v_scale_0 = _S660;
    *v_hardness_0 = _S659;
    *v_sh_coeffs_0 = _S658;
    *v_ch_coeffs_0 = _S657;
    *v_R_0 = _S656;
    *v_t_0 = _S655;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S662, float _S663)
{
    return (F32_atan2((_S662), (_S663)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S664, DiffPair_float_0 * _S665, float _S666)
{
    _d_atan2_0(_S664, _S665, _S666);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S667 = scale_5.x;
    float _S668 = s_primal_ctx_exp_0(_S667);
    float _S669 = scale_5.y;
    float _S670 = s_primal_ctx_exp_0(_S669);
    float sz_5 = scale_5.z - 0.5f * (_S667 + _S669);
    float _S671 = quat_5.y;
    float x2_5 = _S671 * _S671;
    float y2_5 = quat_5.z * quat_5.z;
    float z2_10 = quat_5.w * quat_5.w;
    float xy_5 = quat_5.y * quat_5.z;
    float xz_5 = quat_5.y * quat_5.w;
    float yz_5 = quat_5.z * quat_5.w;
    float wx_5 = quat_5.x * quat_5.y;
    float wy_5 = quat_5.x * quat_5.z;
    float wz_5 = quat_5.x * quat_5.w;
    Matrix<float, 3, 3>  _S672 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S673 = make_float3 (_S668, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S672, _S673) + mean_5;
    float _S674 = -0.5f + sz_5;
    float3  _S675 = make_float3 (_S668 * _S674, _S670, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S672, _S675) + mean_5;
    float _S676 = -0.5f - sz_5;
    float3  _S677 = make_float3 (_S668 * _S676, - _S670, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S672, _S677) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S678 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S679 = length_0(_S678);
    float _S680 = vert0_c_5.z;
    float _S681 = s_primal_ctx_atan2_0(_S679, _S680);
    bool _S682 = _S681 < 0.00100000004749745f;
    float k_2;
    float _S683;
    float _S684;
    float _S685;
    if(_S682)
    {
        float _S686 = 1.0f - _S681 * _S681 / 3.0f;
        float _S687 = _S680 * _S680;
        k_2 = _S686 / _S680;
        _S683 = _S687;
        _S684 = _S686;
        _S685 = 0.0f;
    }
    else
    {
        float _S688 = _S679 * _S679;
        k_2 = _S681 / _S679;
        _S683 = 0.0f;
        _S684 = 0.0f;
        _S685 = _S688;
    }
    float2  _S689 = make_float2 (k_2);
    float2  _S690 = _S678 * make_float2 (k_2);
    float u_18 = _S690.x;
    float v_18 = _S690.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S691 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S692 = dist_coeffs_5[int(1)] + r2_18 * _S691;
    float _S693 = dist_coeffs_5[int(0)] + r2_18 * _S692;
    float radial_9 = 1.0f + r2_18 * _S693;
    float _S694 = 2.0f * dist_coeffs_5[int(4)];
    float _S695 = _S694 * u_18;
    float _S696 = 2.0f * u_18;
    float _S697 = 2.0f * dist_coeffs_5[int(5)];
    float _S698 = _S697 * u_18;
    float _S699 = 2.0f * v_18;
    float2  _S700 = _S690 * make_float2 (radial_9) + make_float2 (_S695 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S696 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S698 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S699 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S701 = _S700 + make_float2 (dist_coeffs_5[int(8)] * _S700.x + dist_coeffs_5[int(9)] * _S700.y, 0.0f);
    float _S702 = fx_5 * _S701.x + cx_5;
    float _S703 = fy_5 * _S701.y + cy_5;
    float2  uv0_7 = make_float2 (_S702, _S703);
    float2  _S704 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S705 = length_0(_S704);
    float _S706 = vert1_c_5.z;
    float _S707 = s_primal_ctx_atan2_0(_S705, _S706);
    bool _S708 = _S707 < 0.00100000004749745f;
    float _S709;
    float _S710;
    float _S711;
    if(_S708)
    {
        float _S712 = 1.0f - _S707 * _S707 / 3.0f;
        float _S713 = _S706 * _S706;
        k_2 = _S712 / _S706;
        _S709 = _S713;
        _S710 = _S712;
        _S711 = 0.0f;
    }
    else
    {
        float _S714 = _S705 * _S705;
        k_2 = _S707 / _S705;
        _S709 = 0.0f;
        _S710 = 0.0f;
        _S711 = _S714;
    }
    float2  _S715 = make_float2 (k_2);
    float2  _S716 = _S704 * make_float2 (k_2);
    float u_19 = _S716.x;
    float v_19 = _S716.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S717 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S718 = dist_coeffs_5[int(1)] + r2_19 * _S717;
    float _S719 = dist_coeffs_5[int(0)] + r2_19 * _S718;
    float radial_10 = 1.0f + r2_19 * _S719;
    float _S720 = _S694 * u_19;
    float _S721 = 2.0f * u_19;
    float _S722 = _S697 * u_19;
    float _S723 = 2.0f * v_19;
    float2  _S724 = _S716 * make_float2 (radial_10) + make_float2 (_S720 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S721 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S722 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S723 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S725 = _S724 + make_float2 (dist_coeffs_5[int(8)] * _S724.x + dist_coeffs_5[int(9)] * _S724.y, 0.0f);
    float _S726 = fx_5 * _S725.x + cx_5;
    float _S727 = fy_5 * _S725.y + cy_5;
    float2  uv1_7 = make_float2 (_S726, _S727);
    float2  _S728 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S729 = length_0(_S728);
    float _S730 = vert2_c_5.z;
    float _S731 = s_primal_ctx_atan2_0(_S729, _S730);
    bool _S732 = _S731 < 0.00100000004749745f;
    float _S733;
    float _S734;
    float _S735;
    if(_S732)
    {
        float _S736 = 1.0f - _S731 * _S731 / 3.0f;
        float _S737 = _S730 * _S730;
        k_2 = _S736 / _S730;
        _S733 = _S737;
        _S734 = _S736;
        _S735 = 0.0f;
    }
    else
    {
        float _S738 = _S729 * _S729;
        k_2 = _S731 / _S729;
        _S733 = 0.0f;
        _S734 = 0.0f;
        _S735 = _S738;
    }
    float2  _S739 = make_float2 (k_2);
    float2  _S740 = _S728 * make_float2 (k_2);
    float u_20 = _S740.x;
    float v_20 = _S740.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S741 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S742 = dist_coeffs_5[int(1)] + r2_20 * _S741;
    float _S743 = dist_coeffs_5[int(0)] + r2_20 * _S742;
    float radial_11 = 1.0f + r2_20 * _S743;
    float _S744 = _S694 * u_20;
    float _S745 = 2.0f * u_20;
    float _S746 = _S697 * u_20;
    float _S747 = 2.0f * v_20;
    float2  _S748 = _S740 * make_float2 (radial_11) + make_float2 (_S744 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S745 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S746 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S747 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S749 = _S748 + make_float2 (dist_coeffs_5[int(8)] * _S748.x + dist_coeffs_5[int(9)] * _S748.y, 0.0f);
    float _S750 = fx_5 * _S749.x + cx_5;
    float _S751 = fy_5 * _S749.y + cy_5;
    float2  uv2_7 = make_float2 (_S750, _S751);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S752 = e0_5.x;
    float _S753 = e1_5.y;
    float _S754 = e0_5.y;
    float _S755 = e1_5.x;
    float _S756 = _S752 * _S753 - _S754 * _S755;
    float _S757 = 1.0f - hardness_5.y;
    float _S758 = -1.0f / _S757;
    float _S759 = _S757 * _S757;
    float _S760 = (F32_max((_S702), (_S726)));
    float _S761 = (F32_min((_S702), (_S726)));
    float _S762 = (F32_max((_S703), (_S727)));
    float _S763 = (F32_min((_S703), (_S727)));
    float3  _S764 = (vert0_c_5 + vert1_c_5 + vert2_c_5) / make_float3 (3.0f);
    float _S765 = _S764.z;
    float _S766 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S767 = _S766 * (_S765 + length_1(_S764));
    Matrix<float, 3, 3>  _S768 = transpose_1(R_5);
    float3  _S769 = mean_5 - - s_primal_ctx_mul_0(_S768, t_5);
    float _S770 = _S769.x;
    float _S771 = _S769.y;
    float _S772 = _S769.z;
    float _S773 = _S770 * _S770 + _S771 * _S771 + _S772 * _S772;
    float _S774 = s_primal_ctx_sqrt_0(_S773);
    float x_17 = _S770 / _S774;
    float3  _S775 = make_float3 (x_17);
    float _S776 = _S774 * _S774;
    float y_8 = _S771 / _S774;
    float z_5 = _S772 / _S774;
    float3  _S777 = make_float3 (z_5);
    float _S778 = - y_8;
    float3  _S779 = make_float3 (_S778);
    float z2_11 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_17 * x_17 - y_8 * y_8;
    float _S780 = 2.0f * x_17;
    float fS1_5 = _S780 * y_8;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S781 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_17;
    float3  _S782 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_8;
    float3  _S783 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S784 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S785 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    float _S786 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_5 * _S786;
    float3  _S787 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_17;
    float3  _S788 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_8;
    float3  _S789 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S790 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S791 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_17 * fC1_5 - y_8 * fS1_5);
    float3  _S792 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_17 * fS1_5 + y_8 * fC1_5);
    float3  _S793 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S778) * sh_coeffs_5[int(1)] + make_float3 (z_5) * sh_coeffs_5[int(2)] - make_float3 (x_17) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]);
    float3  _S794 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S795 = make_float3 (0.0f);
    float3  _S796 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S797 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S798 = make_float3 (_S797);
    float3  _S799 = ch_coeffs_5[int(1)] * make_float3 (_S797);
    float3  _S800 = _S796 + _S799 + make_float3 (0.5f);
    float3  _S801 = _S796 - _S799 + make_float3 (0.5f);
    float3  _S802 = vert1_c_5 - vert0_c_5;
    float3  _S803 = vert2_c_5 - vert0_c_5;
    float3  _S804 = s_primal_ctx_cross_0(_S802, _S803);
    float3  _S805 = normalize_0(_S804);
    float3  _S806 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S805, mean_c_5)))))) * v_normal_1;
    float3  _S807 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S808;
    (&_S808)->primal_0 = _S805;
    (&_S808)->differential_0 = _S807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S809;
    (&_S809)->primal_0 = mean_c_5;
    (&_S809)->differential_0 = _S807;
    s_bwd_prop_dot_0(&_S808, &_S809, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S810 = _S809;
    float3  _S811 = _S806 + _S808.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S812;
    (&_S812)->primal_0 = _S804;
    (&_S812)->differential_0 = _S807;
    s_bwd_normalize_impl_0(&_S812, _S811);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S813;
    (&_S813)->primal_0 = _S802;
    (&_S813)->differential_0 = _S807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S814;
    (&_S814)->primal_0 = _S803;
    (&_S814)->differential_0 = _S807;
    s_bwd_prop_cross_0(&_S813, &_S814, _S812.differential_0);
    float3  _S815 = - _S814.differential_0;
    float3  _S816 = - _S813.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S817;
    (&_S817)->primal_0 = _S801;
    (&_S817)->differential_0 = _S807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S818;
    (&_S818)->primal_0 = _S795;
    (&_S818)->differential_0 = _S807;
    s_bwd_prop_max_0(&_S817, &_S818, v_rgbs_1[int(2)]);
    float3  _S819 = - _S817.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S820;
    (&_S820)->primal_0 = _S800;
    (&_S820)->differential_0 = _S807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S821;
    (&_S821)->primal_0 = _S795;
    (&_S821)->differential_0 = _S807;
    s_bwd_prop_max_0(&_S820, &_S821, v_rgbs_1[int(1)]);
    float3  _S822 = _S798 * (_S819 + _S820.differential_0);
    float3  _S823 = _S817.differential_0 + _S820.differential_0;
    float3  _S824 = make_float3 (0.5f) * - _S823;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S825;
    (&_S825)->primal_0 = _S794;
    (&_S825)->differential_0 = _S807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S826;
    (&_S826)->primal_0 = _S795;
    (&_S826)->differential_0 = _S807;
    s_bwd_prop_max_0(&_S825, &_S826, v_rgbs_1[int(0)]);
    float3  _S827 = _S824 + _S825.differential_0;
    float3  _S828 = _S823 + _S825.differential_0;
    float3  _S829 = _S792 * _S828;
    float3  _S830 = sh_coeffs_5[int(15)] * _S828;
    float3  _S831 = _S790 * _S828;
    float3  _S832 = sh_coeffs_5[int(14)] * _S828;
    float3  _S833 = _S788 * _S828;
    float3  _S834 = sh_coeffs_5[int(13)] * _S828;
    float3  _S835 = _S787 * _S828;
    float3  _S836 = sh_coeffs_5[int(12)] * _S828;
    float3  _S837 = _S789 * _S828;
    float3  _S838 = sh_coeffs_5[int(11)] * _S828;
    float3  _S839 = _S791 * _S828;
    float3  _S840 = sh_coeffs_5[int(10)] * _S828;
    float3  _S841 = _S793 * _S828;
    float3  _S842 = sh_coeffs_5[int(9)] * _S828;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S842.x + _S842.y + _S842.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S830.x + _S830.y + _S830.z);
    float _S843 = _S840.x + _S840.y + _S840.z;
    float _S844 = _S832.x + _S832.y + _S832.z;
    float _S845 = _S838.x + _S838.y + _S838.z;
    float _S846 = _S834.x + _S834.y + _S834.z;
    float _S847 = _S836.x + _S836.y + _S836.z;
    float _S848 = - s_diff_fC2_T_1;
    float3  _S849 = _S784 * _S828;
    float3  _S850 = sh_coeffs_5[int(8)] * _S828;
    float3  _S851 = _S782 * _S828;
    float3  _S852 = sh_coeffs_5[int(7)] * _S828;
    float3  _S853 = _S781 * _S828;
    float3  _S854 = sh_coeffs_5[int(6)] * _S828;
    float3  _S855 = _S783 * _S828;
    float3  _S856 = sh_coeffs_5[int(5)] * _S828;
    float3  _S857 = _S785 * _S828;
    float3  _S858 = sh_coeffs_5[int(4)] * _S828;
    float _S859 = _S856.x + _S856.y + _S856.z;
    float _S860 = _S852.x + _S852.y + _S852.z;
    float _S861 = fTmp1B_5 * _S843 + x_17 * s_diff_fS2_T_1 + y_8 * _S848 + 0.54627424478530884f * (_S858.x + _S858.y + _S858.z);
    float _S862 = fTmp1B_5 * _S844 + y_8 * s_diff_fS2_T_1 + x_17 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S850.x + _S850.y + _S850.z);
    float _S863 = y_8 * - _S862;
    float _S864 = x_17 * _S862;
    float _S865 = z_5 * (1.86588168144226074f * (z_5 * _S847) + -2.28522896766662598f * (y_8 * _S845 + x_17 * _S846) + 0.94617468118667603f * (_S854.x + _S854.y + _S854.z));
    float3  _S866 = make_float3 (0.48860251903533936f) * _S828;
    float3  _S867 = - _S866;
    float3  _S868 = _S775 * _S867;
    float3  _S869 = sh_coeffs_5[int(3)] * _S867;
    float3  _S870 = _S777 * _S866;
    float3  _S871 = sh_coeffs_5[int(2)] * _S866;
    float3  _S872 = _S779 * _S866;
    float3  _S873 = sh_coeffs_5[int(1)] * _S866;
    float _S874 = (_S786 * _S847 + 1.44530570507049561f * (fS1_5 * _S843 + fC1_5 * _S844) + -1.09254848957061768f * (y_8 * _S859 + x_17 * _S860) + _S865 + _S865 + _S871.x + _S871.y + _S871.z) / _S776;
    float _S875 = _S774 * _S874;
    float _S876 = (fTmp0C_5 * _S845 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S848 + fTmp0B_5 * _S859 + _S780 * _S861 + _S863 + _S863 + - (_S873.x + _S873.y + _S873.z)) / _S776;
    float _S877 = _S774 * _S876;
    float _S878 = (fTmp0C_5 * _S846 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S860 + 2.0f * (y_8 * _S861) + _S864 + _S864 + _S869.x + _S869.y + _S869.z) / _S776;
    float _S879 = _S774 * _S878;
    float _S880 = _S772 * - _S874 + _S771 * - _S876 + _S770 * - _S878;
    DiffPair_float_0 _S881;
    (&_S881)->primal_0 = _S773;
    (&_S881)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S881, _S880);
    float _S882 = _S772 * _S881.differential_0;
    float _S883 = _S771 * _S881.differential_0;
    float _S884 = _S770 * _S881.differential_0;
    float3  _S885 = make_float3 (0.282094806432724f) * _S828;
    float3  _S886 = make_float3 (_S879 + _S884 + _S884, _S877 + _S883 + _S883, _S875 + _S882 + _S882);
    float3  _S887 = - - _S886;
    Matrix<float, 3, 3>  _S888 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S889;
    (&_S889)->primal_0 = _S768;
    (&_S889)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S890;
    (&_S890)->primal_0 = t_5;
    (&_S890)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S889, &_S890, _S887);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S891 = _S890;
    Matrix<float, 3, 3>  _S892 = transpose_1(_S889.differential_0);
    DiffPair_float_0 _S893;
    (&_S893)->primal_0 = _S765;
    (&_S893)->differential_0 = 0.0f;
    DiffPair_float_0 _S894;
    (&_S894)->primal_0 = _S767;
    (&_S894)->differential_0 = 0.0f;
    _d_max_0(&_S893, &_S894, v_depth_1);
    float _S895 = _S766 * _S894.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S896;
    (&_S896)->primal_0 = _S764;
    (&_S896)->differential_0 = _S807;
    s_bwd_length_impl_0(&_S896, _S895);
    float3  _S897 = make_float3 (0.3333333432674408f) * (_S896.differential_0 + make_float3 (0.0f, 0.0f, _S893.differential_0 + _S895));
    DiffPair_float_0 _S898;
    (&_S898)->primal_0 = _S763;
    (&_S898)->differential_0 = 0.0f;
    DiffPair_float_0 _S899;
    (&_S899)->primal_0 = _S751;
    (&_S899)->differential_0 = 0.0f;
    _d_min_0(&_S898, &_S899, 0.0f);
    DiffPair_float_0 _S900;
    (&_S900)->primal_0 = _S703;
    (&_S900)->differential_0 = 0.0f;
    DiffPair_float_0 _S901;
    (&_S901)->primal_0 = _S727;
    (&_S901)->differential_0 = 0.0f;
    _d_min_0(&_S900, &_S901, _S898.differential_0);
    DiffPair_float_0 _S902;
    (&_S902)->primal_0 = _S762;
    (&_S902)->differential_0 = 0.0f;
    DiffPair_float_0 _S903;
    (&_S903)->primal_0 = _S751;
    (&_S903)->differential_0 = 0.0f;
    _d_max_0(&_S902, &_S903, 0.0f);
    DiffPair_float_0 _S904;
    (&_S904)->primal_0 = _S703;
    (&_S904)->differential_0 = 0.0f;
    DiffPair_float_0 _S905;
    (&_S905)->primal_0 = _S727;
    (&_S905)->differential_0 = 0.0f;
    _d_max_0(&_S904, &_S905, _S902.differential_0);
    DiffPair_float_0 _S906;
    (&_S906)->primal_0 = _S761;
    (&_S906)->differential_0 = 0.0f;
    DiffPair_float_0 _S907;
    (&_S907)->primal_0 = _S750;
    (&_S907)->differential_0 = 0.0f;
    _d_min_0(&_S906, &_S907, 0.0f);
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S702;
    (&_S908)->differential_0 = 0.0f;
    DiffPair_float_0 _S909;
    (&_S909)->primal_0 = _S726;
    (&_S909)->differential_0 = 0.0f;
    _d_min_0(&_S908, &_S909, _S906.differential_0);
    DiffPair_float_0 _S910;
    (&_S910)->primal_0 = _S760;
    (&_S910)->differential_0 = 0.0f;
    DiffPair_float_0 _S911;
    (&_S911)->primal_0 = _S750;
    (&_S911)->differential_0 = 0.0f;
    _d_max_0(&_S910, &_S911, 0.0f);
    DiffPair_float_0 _S912;
    (&_S912)->primal_0 = _S702;
    (&_S912)->differential_0 = 0.0f;
    DiffPair_float_0 _S913;
    (&_S913)->primal_0 = _S726;
    (&_S913)->differential_0 = 0.0f;
    _d_max_0(&_S912, &_S913, _S910.differential_0);
    DiffPair_float_0 _S914;
    (&_S914)->primal_0 = _S758;
    (&_S914)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S914, -0.0f);
    float _S915 = - (-1.0f * - (_S914.differential_0 / _S759));
    float2  _S916 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S917;
    (&_S917)->primal_0 = e2_1;
    (&_S917)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S917, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S918;
    (&_S918)->primal_0 = e1_5;
    (&_S918)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S918, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S919;
    (&_S919)->primal_0 = e0_5;
    (&_S919)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S919, 0.0f);
    DiffPair_float_0 _S920;
    (&_S920)->primal_0 = _S756;
    (&_S920)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S920, -0.0f);
    float _S921 = - _S920.differential_0;
    float2  _S922 = _S918.differential_0 + make_float2 (_S754 * _S921, _S752 * _S920.differential_0);
    float2  _S923 = - _S922;
    float2  _S924 = _S919.differential_0 + make_float2 (_S753 * _S920.differential_0, _S755 * _S921);
    float2  _S925 = - _S924;
    float2  _S926 = - _S917.differential_0 + _S922;
    float _S927 = fx_5 * (_S907.differential_0 + _S911.differential_0 + _S926.x);
    float2  _S928 = make_float2 (_S927, fy_5 * (_S899.differential_0 + _S903.differential_0 + _S926.y)) + make_float2 (dist_coeffs_5[int(8)] * _S927, dist_coeffs_5[int(9)] * _S927);
    float2  _S929 = _S740 * _S928;
    float _S930 = dist_coeffs_5[int(4)] * _S928.y;
    float _S931 = dist_coeffs_5[int(5)] * _S928.x;
    float _S932 = _S929.x + _S929.y;
    float _S933 = r2_20 * _S932;
    float _S934 = r2_20 * _S933;
    float _S935 = dist_coeffs_5[int(7)] * _S928.y + _S930 + dist_coeffs_5[int(6)] * _S928.x + _S931 + _S743 * _S932 + _S742 * _S933 + _S741 * _S934 + dist_coeffs_5[int(3)] * (r2_20 * _S934);
    float _S936 = v_20 * _S935;
    float _S937 = u_20 * _S935;
    float2  _S938 = make_float2 (radial_11) * _S928 + make_float2 (_S697 * (v_20 * _S928.y) + _S745 * _S931 + 2.0f * (u_20 * _S931) + _S694 * (v_20 * _S928.x) + _S937 + _S937, _S747 * _S930 + 2.0f * (v_20 * _S930) + _S746 * _S928.y + _S744 * _S928.x + _S936 + _S936);
    float3  _S939 = _S813.differential_0 + _S897;
    float3  _S940 = _S815 + _S816 + _S897;
    float3  _S941 = _S814.differential_0 + _S897;
    FixedArray<float3 , 2>  _S942;
    _S942[int(0)] = _S807;
    _S942[int(1)] = _S807;
    _S942[int(1)] = _S822;
    _S942[int(0)] = _S827;
    float3  _S943 = _S942[int(0)];
    float3  _S944 = _S942[int(1)];
    FixedArray<float3 , 16>  _S945;
    _S945[int(0)] = _S807;
    _S945[int(1)] = _S807;
    _S945[int(2)] = _S807;
    _S945[int(3)] = _S807;
    _S945[int(4)] = _S807;
    _S945[int(5)] = _S807;
    _S945[int(6)] = _S807;
    _S945[int(7)] = _S807;
    _S945[int(8)] = _S807;
    _S945[int(9)] = _S807;
    _S945[int(10)] = _S807;
    _S945[int(11)] = _S807;
    _S945[int(12)] = _S807;
    _S945[int(13)] = _S807;
    _S945[int(14)] = _S807;
    _S945[int(15)] = _S807;
    _S945[int(7)] = _S851;
    _S945[int(0)] = _S885;
    _S945[int(1)] = _S872;
    _S945[int(2)] = _S870;
    _S945[int(3)] = _S868;
    _S945[int(4)] = _S857;
    _S945[int(5)] = _S855;
    _S945[int(6)] = _S853;
    _S945[int(15)] = _S829;
    _S945[int(8)] = _S849;
    _S945[int(9)] = _S841;
    _S945[int(10)] = _S839;
    _S945[int(11)] = _S837;
    _S945[int(12)] = _S835;
    _S945[int(13)] = _S833;
    _S945[int(14)] = _S831;
    float3  _S946 = _S945[int(0)];
    float3  _S947 = _S945[int(1)];
    float3  _S948 = _S945[int(2)];
    float3  _S949 = _S945[int(3)];
    float3  _S950 = _S945[int(4)];
    float3  _S951 = _S945[int(5)];
    float3  _S952 = _S945[int(6)];
    float3  _S953 = _S945[int(7)];
    float3  _S954 = _S945[int(8)];
    float3  _S955 = _S945[int(9)];
    float3  _S956 = _S945[int(10)];
    float3  _S957 = _S945[int(11)];
    float3  _S958 = _S945[int(12)];
    float3  _S959 = _S945[int(13)];
    float3  _S960 = _S945[int(14)];
    float3  _S961 = _S945[int(15)];
    float _S962 = _S901.differential_0 + _S905.differential_0;
    float _S963 = _S900.differential_0 + _S904.differential_0;
    float _S964 = _S908.differential_0 + _S912.differential_0;
    float _S965 = _S909.differential_0 + _S913.differential_0;
    float2  _S966 = _S917.differential_0 + _S925;
    float2  _S967 = _S923 + _S924;
    float2  _S968 = make_float2 (0.0f, _S915);
    float2  _S969 = _S728 * _S938;
    float2  _S970 = _S739 * _S938;
    float _S971 = _S969.x + _S969.y;
    if(_S732)
    {
        float _S972 = _S971 / _S733;
        float _S973 = _S734 * - _S972;
        float _S974 = _S731 * (0.3333333432674408f * - (_S730 * _S972));
        k_2 = _S974 + _S974;
        _S733 = _S973;
        _S734 = 0.0f;
    }
    else
    {
        float _S975 = _S971 / _S735;
        float _S976 = _S731 * - _S975;
        k_2 = _S729 * _S975;
        _S733 = 0.0f;
        _S734 = _S976;
    }
    DiffPair_float_0 _S977;
    (&_S977)->primal_0 = _S729;
    (&_S977)->differential_0 = 0.0f;
    DiffPair_float_0 _S978;
    (&_S978)->primal_0 = _S730;
    (&_S978)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S977, &_S978, k_2);
    float _S979 = _S978.differential_0 + _S733;
    float _S980 = _S977.differential_0 + _S734;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S981;
    (&_S981)->primal_0 = _S728;
    (&_S981)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S981, _S980);
    float2  _S982 = _S981.differential_0 + _S970;
    float _S983 = fx_5 * (_S967.x + _S965);
    float2  _S984 = make_float2 (_S983, fy_5 * (_S967.y + _S962)) + make_float2 (dist_coeffs_5[int(8)] * _S983, dist_coeffs_5[int(9)] * _S983);
    float2  _S985 = _S716 * _S984;
    float _S986 = dist_coeffs_5[int(4)] * _S984.y;
    float _S987 = dist_coeffs_5[int(5)] * _S984.x;
    float _S988 = _S985.x + _S985.y;
    float _S989 = r2_19 * _S988;
    float _S990 = r2_19 * _S989;
    float _S991 = dist_coeffs_5[int(7)] * _S984.y + _S986 + dist_coeffs_5[int(6)] * _S984.x + _S987 + _S719 * _S988 + _S718 * _S989 + _S717 * _S990 + dist_coeffs_5[int(3)] * (r2_19 * _S990);
    float _S992 = v_19 * _S991;
    float _S993 = u_19 * _S991;
    float2  _S994 = make_float2 (radial_10) * _S984 + make_float2 (_S697 * (v_19 * _S984.y) + _S721 * _S987 + 2.0f * (u_19 * _S987) + _S694 * (v_19 * _S984.x) + _S993 + _S993, _S723 * _S986 + 2.0f * (v_19 * _S986) + _S722 * _S984.y + _S720 * _S984.x + _S992 + _S992);
    float3  _S995 = _S941 + make_float3 (_S982.x, _S982.y, _S979);
    float2  _S996 = _S704 * _S994;
    float2  _S997 = _S715 * _S994;
    float _S998 = _S996.x + _S996.y;
    if(_S708)
    {
        float _S999 = _S998 / _S709;
        float _S1000 = _S710 * - _S999;
        float _S1001 = _S707 * (0.3333333432674408f * - (_S706 * _S999));
        k_2 = _S1001 + _S1001;
        _S709 = _S1000;
        _S710 = 0.0f;
    }
    else
    {
        float _S1002 = _S998 / _S711;
        float _S1003 = _S707 * - _S1002;
        k_2 = _S705 * _S1002;
        _S709 = 0.0f;
        _S710 = _S1003;
    }
    DiffPair_float_0 _S1004;
    (&_S1004)->primal_0 = _S705;
    (&_S1004)->differential_0 = 0.0f;
    DiffPair_float_0 _S1005;
    (&_S1005)->primal_0 = _S706;
    (&_S1005)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1004, &_S1005, k_2);
    float _S1006 = _S1005.differential_0 + _S709;
    float _S1007 = _S1004.differential_0 + _S710;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1008;
    (&_S1008)->primal_0 = _S704;
    (&_S1008)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S1008, _S1007);
    float2  _S1009 = _S1008.differential_0 + _S997;
    float _S1010 = fx_5 * (_S966.x + _S964);
    float2  _S1011 = make_float2 (_S1010, fy_5 * (_S966.y + _S963)) + make_float2 (dist_coeffs_5[int(8)] * _S1010, dist_coeffs_5[int(9)] * _S1010);
    float2  _S1012 = _S690 * _S1011;
    float _S1013 = dist_coeffs_5[int(4)] * _S1011.y;
    float _S1014 = dist_coeffs_5[int(5)] * _S1011.x;
    float _S1015 = _S1012.x + _S1012.y;
    float _S1016 = r2_18 * _S1015;
    float _S1017 = r2_18 * _S1016;
    float _S1018 = dist_coeffs_5[int(7)] * _S1011.y + _S1013 + dist_coeffs_5[int(6)] * _S1011.x + _S1014 + _S693 * _S1015 + _S692 * _S1016 + _S691 * _S1017 + dist_coeffs_5[int(3)] * (r2_18 * _S1017);
    float _S1019 = v_18 * _S1018;
    float _S1020 = u_18 * _S1018;
    float2  _S1021 = make_float2 (radial_9) * _S1011 + make_float2 (_S697 * (v_18 * _S1011.y) + _S696 * _S1014 + 2.0f * (u_18 * _S1014) + _S694 * (v_18 * _S1011.x) + _S1020 + _S1020, _S699 * _S1013 + 2.0f * (v_18 * _S1013) + _S698 * _S1011.y + _S695 * _S1011.x + _S1019 + _S1019);
    float3  _S1022 = _S939 + make_float3 (_S1009.x, _S1009.y, _S1006);
    float2  _S1023 = _S678 * _S1021;
    float2  _S1024 = _S689 * _S1021;
    float _S1025 = _S1023.x + _S1023.y;
    if(_S682)
    {
        float _S1026 = _S1025 / _S683;
        float _S1027 = _S684 * - _S1026;
        float _S1028 = _S681 * (0.3333333432674408f * - (_S680 * _S1026));
        k_2 = _S1028 + _S1028;
        _S683 = _S1027;
        _S684 = 0.0f;
    }
    else
    {
        float _S1029 = _S1025 / _S685;
        float _S1030 = _S681 * - _S1029;
        k_2 = _S679 * _S1029;
        _S683 = 0.0f;
        _S684 = _S1030;
    }
    DiffPair_float_0 _S1031;
    (&_S1031)->primal_0 = _S679;
    (&_S1031)->differential_0 = 0.0f;
    DiffPair_float_0 _S1032;
    (&_S1032)->primal_0 = _S680;
    (&_S1032)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1031, &_S1032, k_2);
    float _S1033 = _S1032.differential_0 + _S683;
    float _S1034 = _S1031.differential_0 + _S684;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1035;
    (&_S1035)->primal_0 = _S678;
    (&_S1035)->differential_0 = _S916;
    s_bwd_length_impl_1(&_S1035, _S1034);
    float2  _S1036 = _S1035.differential_0 + _S1024;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1037;
    (&_S1037)->primal_0 = vert2_c_5;
    (&_S1037)->differential_0 = _S807;
    s_bwd_length_impl_0(&_S1037, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1038;
    (&_S1038)->primal_0 = vert1_c_5;
    (&_S1038)->differential_0 = _S807;
    s_bwd_length_impl_0(&_S1038, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1039;
    (&_S1039)->primal_0 = vert0_c_5;
    (&_S1039)->differential_0 = _S807;
    s_bwd_length_impl_0(&_S1039, 0.0f);
    float3  _S1040 = _S1037.differential_0 + _S995;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1041;
    (&_S1041)->primal_0 = R_5;
    (&_S1041)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1042;
    (&_S1042)->primal_0 = vert2_5;
    (&_S1042)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1041, &_S1042, _S1040);
    float3  _S1043 = _S1038.differential_0 + _S1022;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1044;
    (&_S1044)->primal_0 = R_5;
    (&_S1044)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1045;
    (&_S1045)->primal_0 = vert1_5;
    (&_S1045)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1044, &_S1045, _S1043);
    float3  _S1046 = _S1039.differential_0 + _S940 + make_float3 (_S1036.x, _S1036.y, _S1033);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1047;
    (&_S1047)->primal_0 = R_5;
    (&_S1047)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1048;
    (&_S1048)->primal_0 = vert0_5;
    (&_S1048)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1047, &_S1048, _S1046);
    float3  _S1049 = _S1042.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1050;
    (&_S1050)->primal_0 = _S672;
    (&_S1050)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1051;
    (&_S1051)->primal_0 = _S677;
    (&_S1051)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1050, &_S1051, _S1049);
    float _S1052 = - _S1051.differential_0.y;
    float _S1053 = _S676 * _S1051.differential_0.x;
    float _S1054 = - (_S668 * _S1051.differential_0.x);
    float3  _S1055 = _S1045.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1056;
    (&_S1056)->primal_0 = _S672;
    (&_S1056)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1057;
    (&_S1057)->primal_0 = _S675;
    (&_S1057)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1056, &_S1057, _S1055);
    float _S1058 = _S668 * _S1057.differential_0.x;
    float _S1059 = _S674 * _S1057.differential_0.x;
    float3  _S1060 = _S1048.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1061;
    (&_S1061)->primal_0 = _S672;
    (&_S1061)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1062;
    (&_S1062)->primal_0 = _S673;
    (&_S1062)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1061, &_S1062, _S1060);
    Matrix<float, 3, 3>  _S1063 = transpose_1(_S1050.differential_0 + _S1056.differential_0 + _S1061.differential_0);
    float _S1064 = 2.0f * - _S1063.rows[int(2)].z;
    float _S1065 = 2.0f * _S1063.rows[int(2)].y;
    float _S1066 = 2.0f * _S1063.rows[int(2)].x;
    float _S1067 = 2.0f * _S1063.rows[int(1)].z;
    float _S1068 = 2.0f * - _S1063.rows[int(1)].y;
    float _S1069 = 2.0f * _S1063.rows[int(1)].x;
    float _S1070 = 2.0f * _S1063.rows[int(0)].z;
    float _S1071 = 2.0f * _S1063.rows[int(0)].y;
    float _S1072 = 2.0f * - _S1063.rows[int(0)].x;
    float _S1073 = - _S1069 + _S1071;
    float _S1074 = _S1066 + - _S1070;
    float _S1075 = - _S1065 + _S1067;
    float _S1076 = _S1065 + _S1067;
    float _S1077 = _S1066 + _S1070;
    float _S1078 = _S1069 + _S1071;
    float _S1079 = quat_5.w * (_S1068 + _S1072);
    float _S1080 = quat_5.z * (_S1064 + _S1072);
    float _S1081 = quat_5.y * (_S1064 + _S1068);
    float _S1082 = quat_5.x * _S1073 + quat_5.z * _S1076 + quat_5.y * _S1077 + _S1079 + _S1079;
    float _S1083 = quat_5.x * _S1074 + quat_5.w * _S1076 + quat_5.y * _S1078 + _S1080 + _S1080;
    float _S1084 = quat_5.x * _S1075 + quat_5.w * _S1077 + quat_5.z * _S1078 + _S1081 + _S1081;
    float _S1085 = quat_5.w * _S1073 + quat_5.z * _S1074 + quat_5.y * _S1075;
    float _S1086 = _S1054 + _S1058;
    float _S1087 = 0.5f * - _S1086;
    float _S1088 = _S1052 + _S1057.differential_0.y;
    DiffPair_float_0 _S1089;
    (&_S1089)->primal_0 = _S669;
    (&_S1089)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1089, _S1088);
    float _S1090 = _S1087 + _S1089.differential_0;
    float _S1091 = _S1053 + _S1059 + _S1062.differential_0.x;
    DiffPair_float_0 _S1092;
    (&_S1092)->primal_0 = _S667;
    (&_S1092)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1092, _S1091);
    float _S1093 = _S1087 + _S1092.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1094;
    (&_S1094)->primal_0 = mean_c_5;
    (&_S1094)->differential_0 = _S807;
    s_bwd_length_impl_0(&_S1094, 0.0f);
    float3  _S1095 = _S1094.differential_0 + _S810.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1096;
    (&_S1096)->primal_0 = R_5;
    (&_S1096)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1097;
    (&_S1097)->primal_0 = mean_5;
    (&_S1097)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S1096, &_S1097, _S1095);
    float3  _S1098 = _S1040 + _S1043 + _S1046 + _S1095 + _S891.differential_0;
    Matrix<float, 3, 3>  _S1099 = _S1041.differential_0 + _S1044.differential_0 + _S1047.differential_0 + _S1096.differential_0 + _S892;
    float3  _S1100 = make_float3 (_S1093, _S1090, _S1086);
    float4  _S1101 = make_float4 (0.0f);
    *&((&_S1101)->w) = _S1082;
    *&((&_S1101)->z) = _S1083;
    *&((&_S1101)->y) = _S1084;
    *&((&_S1101)->x) = _S1085;
    float4  _S1102 = _S1101;
    float3  _S1103 = _S1049 + _S1055 + _S1060 + _S1097.differential_0 + _S886;
    *v_mean_1 = _S1103;
    *v_quat_1 = _S1102;
    *v_scale_1 = _S1100;
    *v_hardness_1 = _S968;
    (*v_sh_coeffs_1)[int(0)] = _S946;
    (*v_sh_coeffs_1)[int(1)] = _S947;
    (*v_sh_coeffs_1)[int(2)] = _S948;
    (*v_sh_coeffs_1)[int(3)] = _S949;
    (*v_sh_coeffs_1)[int(4)] = _S950;
    (*v_sh_coeffs_1)[int(5)] = _S951;
    (*v_sh_coeffs_1)[int(6)] = _S952;
    (*v_sh_coeffs_1)[int(7)] = _S953;
    (*v_sh_coeffs_1)[int(8)] = _S954;
    (*v_sh_coeffs_1)[int(9)] = _S955;
    (*v_sh_coeffs_1)[int(10)] = _S956;
    (*v_sh_coeffs_1)[int(11)] = _S957;
    (*v_sh_coeffs_1)[int(12)] = _S958;
    (*v_sh_coeffs_1)[int(13)] = _S959;
    (*v_sh_coeffs_1)[int(14)] = _S960;
    (*v_sh_coeffs_1)[int(15)] = _S961;
    (*v_ch_coeffs_1)[int(0)] = _S943;
    (*v_ch_coeffs_1)[int(1)] = _S944;
    *v_R_1 = _S1099;
    *v_t_1 = _S1098;
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
    bool _S1104;
    if((*u_21) >= 0.0f)
    {
        _S1104 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1104 = false;
    }
    if(_S1104)
    {
        _S1104 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1104 = false;
    }
    if(_S1104)
    {
        _S1104 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1104 = false;
    }
    return _S1104;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1105 = *dpx_12;
    bool _S1106;
    if(((*dpx_12).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1106 = ((*dpx_12).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1106 = false;
    }
    float _S1107;
    if(_S1106)
    {
        _S1107 = dOut_11;
    }
    else
    {
        _S1107 = 0.0f;
    }
    dpx_12->primal_0 = _S1105.primal_0;
    dpx_12->differential_0 = _S1107;
    DiffPair_float_0 _S1108 = *dpMin_0;
    if((_S1105.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1107 = dOut_11;
    }
    else
    {
        _S1107 = 0.0f;
    }
    dpMin_0->primal_0 = _S1108.primal_0;
    dpMin_0->differential_0 = _S1107;
    DiffPair_float_0 _S1109 = *dpMax_0;
    if(((*dpx_12).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1107 = dOut_11;
    }
    else
    {
        _S1107 = 0.0f;
    }
    dpMax_0->primal_0 = _S1109.primal_0;
    dpMax_0->differential_0 = _S1107;
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
        DiffPair_float_0 _S1110 = *dpx_13;
        float _S1111 = val_0 * (*dpy_5).primal_0 / (*dpx_13).primal_0 * dOut_12;
        dpx_13->primal_0 = (*dpx_13).primal_0;
        dpx_13->differential_0 = _S1111;
        float _S1112 = val_0 * (F32_log((_S1110.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1112;
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
    bool _S1113;
    if(u_22 >= 0.0f)
    {
        _S1113 = v_22 >= 0.0f;
    }
    else
    {
        _S1113 = false;
    }
    if(_S1113)
    {
        _S1113 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1113 = false;
    }
    if(_S1113)
    {
        _S1113 = t_7 >= 0.0f;
    }
    else
    {
        _S1113 = false;
    }
    if(!_S1113)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1114;
    if(opac_0 < 0.0f)
    {
        _S1114 = 0.0f;
    }
    else
    {
        _S1114 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1114;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1115, float _S1116, float _S1117)
{
    return clamp_0(_S1115, _S1116, _S1117);
}

inline __device__ float s_primal_ctx_pow_0(float _S1118, float _S1119)
{
    return (F32_pow((_S1118), (_S1119)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1120, DiffPair_float_0 * _S1121, float _S1122)
{
    _d_pow_0(_S1120, _S1121, _S1122);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1123, DiffPair_float_0 * _S1124, DiffPair_float_0 * _S1125, float _S1126)
{
    _d_clamp_0(_S1123, _S1124, _S1125, _S1126);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1127 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1128 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1129 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1130 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1131 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1129);
    float d_2 = 1.0f / _S1131;
    float _S1132 = _S1131 * _S1131;
    float3  _S1133 = - _S1130;
    float _S1134 = s_primal_ctx_dot_0(_S1133, v2v0_2);
    float u_23 = d_2 * _S1134;
    float _S1135 = s_primal_ctx_dot_0(_S1130, v1v0_2);
    float v_23 = d_2 * _S1135;
    float3  _S1136 = - _S1129;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1136, rov0_2);
    bool _S1137;
    if(u_23 >= 0.0f)
    {
        _S1137 = v_23 >= 0.0f;
    }
    else
    {
        _S1137 = false;
    }
    if(_S1137)
    {
        _S1137 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1137 = false;
    }
    if(_S1137)
    {
        _S1137 = t_8 >= 0.0f;
    }
    else
    {
        _S1137 = false;
    }
    bool _S1138 = !!_S1137;
    float _S1139;
    float _S1140;
    float _S1141;
    float _S1142;
    float _S1143;
    float _S1144;
    float _S1145;
    float _S1146;
    float _S1147;
    float _S1148;
    float _S1149;
    if(_S1138)
    {
        float _S1150 = (F32_min((u_23), (v_23)));
        float _S1151 = s_primal_ctx_sqrt_0(0.5f);
        float _S1152 = _S1151 * (1.0f - u_23 - v_23);
        float _S1153 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1150), (_S1152))) * _S1153;
        float _S1154 = _S1127.primal_0.y;
        float _S1155 = 1.0f - opac_1;
        float _S1156 = 1.0f - s_primal_ctx_clamp_0(_S1154, 0.0f, 0.99989998340606689f);
        float _S1157 = 1.0f / _S1156;
        float _S1158 = _S1156 * _S1156;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1155, _S1157);
        float o_1 = _S1127.primal_0.x;
        bool _S1159 = opac_1 < 0.0f;
        if(_S1159)
        {
            _S1139 = 0.0f;
        }
        else
        {
            _S1139 = o_1 * w_1;
        }
        _S1137 = _S1159;
        _S1140 = o_1;
        _S1141 = w_1;
        _S1142 = _S1155;
        _S1143 = _S1157;
        _S1144 = _S1158;
        _S1145 = _S1154;
        _S1146 = _S1153;
        _S1147 = _S1150;
        _S1148 = _S1152;
        _S1149 = _S1151;
    }
    else
    {
        _S1137 = false;
        _S1139 = 0.0f;
        _S1140 = 0.0f;
        _S1141 = 0.0f;
        _S1142 = 0.0f;
        _S1143 = 0.0f;
        _S1144 = 0.0f;
        _S1145 = 0.0f;
        _S1146 = 0.0f;
        _S1147 = 0.0f;
        _S1148 = 0.0f;
        _S1149 = 0.0f;
    }
    float2  _S1160 = make_float2 (0.0f);
    float2  _S1161;
    if(_S1138)
    {
        if(_S1137)
        {
            _S1139 = 0.0f;
            _S1140 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1162;
            (&_S1162)->primal_0 = _S1139;
            (&_S1162)->differential_0 = 0.0f;
            DiffPair_float_0 _S1163;
            (&_S1163)->primal_0 = 0.99500000476837158f;
            (&_S1163)->differential_0 = 0.0f;
            _d_min_0(&_S1162, &_S1163, _s_dOut_3);
            float _S1164 = _S1140 * _S1162.differential_0;
            _S1139 = _S1141 * _S1162.differential_0;
            _S1140 = _S1164;
        }
        float _S1165 = - _S1140;
        DiffPair_float_0 _S1166;
        (&_S1166)->primal_0 = _S1142;
        (&_S1166)->differential_0 = 0.0f;
        DiffPair_float_0 _S1167;
        (&_S1167)->primal_0 = _S1143;
        (&_S1167)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1166, &_S1167, _S1165);
        float _S1168 = - - (_S1167.differential_0 / _S1144);
        float s_diff_opac_T_0 = - _S1166.differential_0;
        DiffPair_float_0 _S1169;
        (&_S1169)->primal_0 = _S1145;
        (&_S1169)->differential_0 = 0.0f;
        DiffPair_float_0 _S1170;
        (&_S1170)->primal_0 = 0.0f;
        (&_S1170)->differential_0 = 0.0f;
        DiffPair_float_0 _S1171;
        (&_S1171)->primal_0 = 0.99989998340606689f;
        (&_S1171)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1169, &_S1170, &_S1171, _S1168);
        float _S1172 = _S1146 * s_diff_opac_T_0;
        DiffPair_float_0 _S1173;
        (&_S1173)->primal_0 = _S1147;
        (&_S1173)->differential_0 = 0.0f;
        DiffPair_float_0 _S1174;
        (&_S1174)->primal_0 = _S1148;
        (&_S1174)->differential_0 = 0.0f;
        _d_min_0(&_S1173, &_S1174, _S1172);
        float _S1175 = - (_S1149 * _S1174.differential_0);
        DiffPair_float_0 _S1176;
        (&_S1176)->primal_0 = u_23;
        (&_S1176)->differential_0 = 0.0f;
        DiffPair_float_0 _S1177;
        (&_S1177)->primal_0 = v_23;
        (&_S1177)->differential_0 = 0.0f;
        _d_min_0(&_S1176, &_S1177, _S1173.differential_0);
        float2  _S1178 = make_float2 (_S1139, _S1169.differential_0);
        float _S1179 = _S1175 + _S1177.differential_0;
        _S1139 = _S1175 + _S1176.differential_0;
        _S1140 = _S1179;
        _S1161 = _S1178;
    }
    else
    {
        _S1139 = 0.0f;
        _S1140 = 0.0f;
        _S1161 = _S1160;
    }
    float3  _S1180 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1181;
    (&_S1181)->primal_0 = _S1136;
    (&_S1181)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1182;
    (&_S1182)->primal_0 = rov0_2;
    (&_S1182)->differential_0 = _S1180;
    s_bwd_prop_dot_0(&_S1181, &_S1182, 0.0f);
    float3  _S1183 = - _S1181.differential_0;
    float _S1184 = d_2 * _S1140;
    float _S1185 = _S1135 * _S1140;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1186;
    (&_S1186)->primal_0 = _S1130;
    (&_S1186)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1187;
    (&_S1187)->primal_0 = v1v0_2;
    (&_S1187)->differential_0 = _S1180;
    s_bwd_prop_dot_0(&_S1186, &_S1187, _S1184);
    float _S1188 = d_2 * _S1139;
    float _S1189 = _S1134 * _S1139;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1190;
    (&_S1190)->primal_0 = _S1133;
    (&_S1190)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1191;
    (&_S1191)->primal_0 = v2v0_2;
    (&_S1191)->differential_0 = _S1180;
    s_bwd_prop_dot_0(&_S1190, &_S1191, _S1188);
    float3  _S1192 = - _S1190.differential_0;
    float _S1193 = - ((_S1185 + _S1189) / _S1132);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1194;
    (&_S1194)->primal_0 = _S1128.primal_0;
    (&_S1194)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1195;
    (&_S1195)->primal_0 = _S1129;
    (&_S1195)->differential_0 = _S1180;
    s_bwd_prop_dot_0(&_S1194, &_S1195, _S1193);
    float3  _S1196 = _S1186.differential_0 + _S1192;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1197;
    (&_S1197)->primal_0 = rov0_2;
    (&_S1197)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1198;
    (&_S1198)->primal_0 = _S1128.primal_0;
    (&_S1198)->differential_0 = _S1180;
    s_bwd_prop_cross_0(&_S1197, &_S1198, _S1196);
    float3  _S1199 = _S1183 + _S1195.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1200;
    (&_S1200)->primal_0 = v1v0_2;
    (&_S1200)->differential_0 = _S1180;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1201;
    (&_S1201)->primal_0 = v2v0_2;
    (&_S1201)->differential_0 = _S1180;
    s_bwd_prop_cross_0(&_S1200, &_S1201, _S1199);
    float3  _S1202 = _S1182.differential_0 + _S1197.differential_0;
    float3  _S1203 = _S1191.differential_0 + _S1201.differential_0;
    float3  _S1204 = _S1187.differential_0 + _S1200.differential_0;
    float3  _S1205 = - _S1202 + - _S1203 + - _S1204;
    float3  _S1206 = _S1194.differential_0 + _S1198.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1206;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1202;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1161;
    FixedArray<float3 , 3>  _S1207;
    _S1207[int(0)] = _S1180;
    _S1207[int(1)] = _S1180;
    _S1207[int(2)] = _S1180;
    _S1207[int(2)] = _S1203;
    _S1207[int(0)] = _S1205;
    _S1207[int(1)] = _S1204;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1207;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1208, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1209, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1210, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1211, float _S1212)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1208, _S1209, _S1210, _S1211, _S1212);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1213 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1214 = { _S1213, _S1213, _S1213 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1214;
    float2  _S1215 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1213;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1213;
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
    float3  _S1216 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1217 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1218 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1216);
    float d_4 = 1.0f / _S1218;
    float _S1219 = _S1218 * _S1218;
    float3  _S1220 = - _S1217;
    float _S1221 = s_primal_ctx_dot_0(_S1220, v2v0_4);
    float u_25 = d_4 * _S1221;
    float _S1222 = s_primal_ctx_dot_0(_S1217, v1v0_4);
    float v_25 = d_4 * _S1222;
    float3  _S1223 = - _S1216;
    float3  _S1224 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1225 = make_float3 (v_25) * dpcolor_0;
    float3  _S1226 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1227 = make_float3 (u_25) * dpcolor_0;
    float3  _S1228 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1229 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1230 = - (_S1228.x + _S1228.y + _S1228.z);
    float _S1231 = d_4 * dpdepth_0;
    float _S1232 = s_primal_ctx_dot_0(_S1223, rov0_4) * dpdepth_0;
    float3  _S1233 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
    (&_S1234)->primal_0 = _S1223;
    (&_S1234)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1235;
    (&_S1235)->primal_0 = rov0_4;
    (&_S1235)->differential_0 = _S1233;
    s_bwd_prop_dot_0(&_S1234, &_S1235, _S1231);
    float3  _S1236 = - _S1234.differential_0;
    float _S1237 = _S1230 + _S1224.x + _S1224.y + _S1224.z;
    float _S1238 = d_4 * _S1237;
    float _S1239 = _S1222 * _S1237;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1240;
    (&_S1240)->primal_0 = _S1217;
    (&_S1240)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1241;
    (&_S1241)->primal_0 = v1v0_4;
    (&_S1241)->differential_0 = _S1233;
    s_bwd_prop_dot_0(&_S1240, &_S1241, _S1238);
    float _S1242 = _S1230 + _S1226.x + _S1226.y + _S1226.z;
    float _S1243 = d_4 * _S1242;
    float _S1244 = _S1221 * _S1242;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1245;
    (&_S1245)->primal_0 = _S1220;
    (&_S1245)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1246;
    (&_S1246)->primal_0 = v2v0_4;
    (&_S1246)->differential_0 = _S1233;
    s_bwd_prop_dot_0(&_S1245, &_S1246, _S1243);
    float3  _S1247 = - _S1245.differential_0;
    float _S1248 = - ((_S1232 + _S1239 + _S1244) / _S1219);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1249;
    (&_S1249)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1249)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1250;
    (&_S1250)->primal_0 = _S1216;
    (&_S1250)->differential_0 = _S1233;
    s_bwd_prop_dot_0(&_S1249, &_S1250, _S1248);
    float3  _S1251 = _S1240.differential_0 + _S1247;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1252;
    (&_S1252)->primal_0 = rov0_4;
    (&_S1252)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1253;
    (&_S1253)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1253)->differential_0 = _S1233;
    s_bwd_prop_cross_0(&_S1252, &_S1253, _S1251);
    float3  _S1254 = _S1236 + _S1250.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1255;
    (&_S1255)->primal_0 = v1v0_4;
    (&_S1255)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1256;
    (&_S1256)->primal_0 = v2v0_4;
    (&_S1256)->differential_0 = _S1233;
    s_bwd_prop_cross_0(&_S1255, &_S1256, _S1254);
    float3  _S1257 = _S1235.differential_0 + _S1252.differential_0;
    float3  _S1258 = _S1246.differential_0 + _S1256.differential_0;
    float3  _S1259 = _S1241.differential_0 + _S1255.differential_0;
    float3  _S1260 = - _S1257 + - _S1258 + - _S1259;
    float3  _S1261 = _S1249.differential_0 + _S1253.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1261;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1257;
    FixedArray<float3 , 3>  _S1262;
    _S1262[int(0)] = _S1233;
    _S1262[int(1)] = _S1233;
    _S1262[int(2)] = _S1233;
    _S1262[int(2)] = _S1225;
    _S1262[int(1)] = _S1227;
    _S1262[int(0)] = _S1229;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1262;
    FixedArray<float3 , 3>  _S1263;
    _S1263[int(0)] = _S1233;
    _S1263[int(1)] = _S1233;
    _S1263[int(2)] = _S1233;
    _S1263[int(2)] = _S1258;
    _S1263[int(0)] = _S1260;
    _S1263[int(1)] = _S1259;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1263;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1264, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1265, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1266, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1267, float3  _S1268, float _S1269)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1264, _S1265, _S1266, _S1267, _S1268, _S1269);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1270 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1271 = { _S1270, _S1270, _S1270 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1271;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1271;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1270;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

