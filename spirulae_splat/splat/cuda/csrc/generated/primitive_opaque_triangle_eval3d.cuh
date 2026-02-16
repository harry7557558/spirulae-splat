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

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_0, float4  quat_0, float3  scale_0, float2  hardness_0, FixedArray<float3 , 16>  sh_coeffs_0, FixedArray<float3 , 2>  ch_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float * depth_0, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_0)
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int((F32_floor((xmin_0)))), int((F32_floor((ymin_0)))), int((F32_ceil((xmax_0)))), int((F32_ceil((ymax_0)))));
        float3  _S93 = (vert0_c_0 + vert1_c_0 + vert2_c_0) / make_float3 (3.0f);
        float _S94 = _S93.z;
        *depth_0 = (F32_max((_S94), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S94 + length_1(_S93)))));
        float3  _S95 = mean_0 - - mul_0(transpose_1(R_0), t_0);
        float _S96 = _S95.x;
        float _S97 = _S95.y;
        float _S98 = _S95.z;
        float norm_0 = (F32_sqrt((_S96 * _S96 + _S97 * _S97 + _S98 * _S98)));
        float x_9 = _S96 / norm_0;
        float y_3 = _S97 / norm_0;
        float z_0 = _S98 / norm_0;
        float z2_1 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_9 * x_9 - y_3 * y_3;
        float fS1_0 = 2.0f * x_9 * y_3;
        float fTmp0C_0 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        float3  color_0 = make_float3 (0.282094806432724f) * sh_coeffs_0[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * sh_coeffs_0[int(1)] + make_float3 (z_0) * sh_coeffs_0[int(2)] - make_float3 (x_9) * sh_coeffs_0[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * sh_coeffs_0[int(4)] + make_float3 (fTmp0B_0 * y_3) * sh_coeffs_0[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * sh_coeffs_0[int(6)] + make_float3 (fTmp0B_0 * x_9) * sh_coeffs_0[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * sh_coeffs_0[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_0 + y_3 * fC1_0)) * sh_coeffs_0[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * sh_coeffs_0[int(10)] + make_float3 (fTmp0C_0 * y_3) * sh_coeffs_0[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * sh_coeffs_0[int(12)] + make_float3 (fTmp0C_0 * x_9) * sh_coeffs_0[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * sh_coeffs_0[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_0 - y_3 * fS1_0)) * sh_coeffs_0[int(15)]);
        float3  _S99 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_0 + ch_coeffs_0[int(0)] + make_float3 (0.5f), _S99);
        float3  _S100 = color_0 - ch_coeffs_0[int(0)] * make_float3 (0.5f);
        float3  _S101 = ch_coeffs_0[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S100 + _S101 + make_float3 (0.5f), _S99);
        (*rgbs_0)[int(2)] = max_0(_S100 - _S101 + make_float3 (0.5f), _S99);
        (*verts_0)[int(0)] = vert0_0;
        (*verts_0)[int(1)] = vert1_0;
        (*verts_0)[int(2)] = vert2_0;
        float3  _S102 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S102 * make_float3 (float(- (F32_sign((dot_0(_S102, mean_c_0))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_1, float4  quat_1, float3  scale_1, float2  hardness_1, FixedArray<float3 , 16>  sh_coeffs_1, FixedArray<float3 , 2>  ch_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float * depth_1, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_1)
{
    float2  _S103;
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
    float _S115;
    float2  _S116;
    float _S117;
    float _S118;
    bool _S119;
    bool _S120;
    bool _S121;
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_1, mean_1) + t_1;
        float _S122 = length_1(mean_c_1);
        bool _S123;
        if(_S122 < near_plane_1)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S122 > far_plane_1;
        }
        if(_S123)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S124 = scale_1.x;
        float sx_1 = (F32_exp((_S124)));
        float _S125 = scale_1.y;
        float sy_1 = (F32_exp((_S125)));
        float sz_1 = scale_1.z - 0.5f * (_S124 + _S125);
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
        Matrix<float, 3, 3>  _S126 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_2), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_2), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1)));
        float3  vert0_1 = mul_0(_S126, make_float3 (sx_1, 0.0f, 0.0f)) + mean_1;
        float3  vert1_1 = mul_0(_S126, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_1;
        float3  vert2_1 = mul_0(_S126, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_1;
        float3  vert0_c_1 = mul_0(R_1, vert0_1) + t_1;
        float3  vert1_c_1 = mul_0(R_1, vert1_1) + t_1;
        float3  vert2_c_1 = mul_0(R_1, vert2_1) + t_1;
        float _S127 = length_1(vert0_c_1);
        float _S128 = length_1(vert1_c_1);
        float _S129 = length_1(vert2_c_1);
        if(_S127 < near_plane_1)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S127 > far_plane_1;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S128 < near_plane_1;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S128 > far_plane_1;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S129 < near_plane_1;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S129 > far_plane_1;
        }
        if(_S123)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_2;
        float k_0;
        for(;;)
        {
            float2  _S130 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_2 = length_0(_S130);
            float _S131 = vert0_c_1.z;
            float theta_0 = (F32_atan2((r_2), (_S131)));
            if(theta_0 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S131;
            }
            else
            {
                k_0 = theta_0 / r_2;
            }
            float2  uv0_3 = _S130 * make_float2 (k_0);
            float2  _S132 = make_float2 (1.0f, 0.0f);
            _S103 = _S132;
            _S104 = dist_coeffs_1[int(0)];
            _S105 = dist_coeffs_1[int(1)];
            _S106 = dist_coeffs_1[int(2)];
            _S107 = dist_coeffs_1[int(3)];
            _S108 = dist_coeffs_1[int(4)];
            _S109 = dist_coeffs_1[int(5)];
            _S110 = dist_coeffs_1[int(6)];
            _S111 = dist_coeffs_1[int(7)];
            _S112 = dist_coeffs_1[int(8)];
            _S113 = dist_coeffs_1[int(9)];
            float u_6 = uv0_3.x;
            float v_6 = uv0_3.y;
            float _S133 = 0.0f * v_6;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float s_diff_r2_6 = u_6 + u_6 + (_S133 + _S133);
            float _S134 = dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)];
            float _S135 = dist_coeffs_1[int(1)] + r2_6 * _S134;
            float _S136 = dist_coeffs_1[int(0)] + r2_6 * _S135;
            float _S137 = s_diff_r2_6 * _S136 + (s_diff_r2_6 * _S135 + (s_diff_r2_6 * _S134 + s_diff_r2_6 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float radial_3 = 1.0f + r2_6 * _S136;
            float _S138 = 2.0f * dist_coeffs_1[int(4)];
            _S114 = _S138;
            float _S139 = _S138 * u_6;
            float _S140 = 2.0f * u_6;
            float s_diff_du_0 = _S138 * v_6 + 0.0f * _S139 + (s_diff_r2_6 + (_S140 + _S140)) * dist_coeffs_1[int(5)] + s_diff_r2_6 * dist_coeffs_1[int(6)];
            float _S141 = 2.0f * dist_coeffs_1[int(5)];
            _S115 = _S141;
            float _S142 = _S141 * u_6;
            float _S143 = 2.0f * v_6;
            float2  _S144 = _S132 * make_float2 (radial_3) + make_float2 (_S137) * uv0_3 + make_float2 (s_diff_du_0, _S141 * v_6 + 0.0f * _S142 + (s_diff_r2_6 + (_S133 + 0.0f * _S143)) * dist_coeffs_1[int(4)] + s_diff_r2_6 * dist_coeffs_1[int(7)]);
            float2  _S145 = _S144 + make_float2 (_S144.x * dist_coeffs_1[int(8)] + _S144.y * dist_coeffs_1[int(9)], 0.0f);
            float2  _S146 = make_float2 (0.0f, 1.0f);
            _S116 = _S146;
            float _S147 = 0.0f * u_6;
            float s_diff_r2_7 = _S147 + _S147 + (v_6 + v_6);
            float _S148 = s_diff_r2_7 * _S136 + (s_diff_r2_7 * _S135 + (s_diff_r2_7 * _S134 + s_diff_r2_7 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float _S149 = 0.0f * _S138;
            _S117 = _S149;
            float s_diff_du_1 = _S149 * v_6 + _S139 + (s_diff_r2_7 + (_S147 + 0.0f * _S140)) * dist_coeffs_1[int(5)] + s_diff_r2_7 * dist_coeffs_1[int(6)];
            float _S150 = 0.0f * _S141;
            _S118 = _S150;
            float2  _S151 = _S146 * make_float2 (radial_3) + make_float2 (_S148) * uv0_3 + make_float2 (s_diff_du_1, _S150 * v_6 + _S142 + (s_diff_r2_7 + (_S143 + _S143)) * dist_coeffs_1[int(4)] + s_diff_r2_7 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S152 = transpose_0(makeMatrix<float, 2, 2> (_S145, _S151 + make_float2 (_S151.x * dist_coeffs_1[int(8)] + _S151.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S153 = !((F32_min((determinant_0(_S152)), ((F32_min((_S152.rows[int(0)].x), (_S152.rows[int(1)].y)))))) > 0.0f);
            _S119 = _S153;
            if(_S153)
            {
                uv0_2 = uv0_3;
                break;
            }
            float2  _S154 = uv0_3 * make_float2 (radial_3) + make_float2 (_S139 * v_6 + dist_coeffs_1[int(5)] * (r2_6 + _S140 * u_6) + dist_coeffs_1[int(6)] * r2_6, _S142 * v_6 + dist_coeffs_1[int(4)] * (r2_6 + _S143 * v_6) + dist_coeffs_1[int(7)] * r2_6);
            float2  _S155 = _S154 + make_float2 (dist_coeffs_1[int(8)] * _S154.x + dist_coeffs_1[int(9)] * _S154.y, 0.0f);
            uv0_2 = make_float2 (fx_1 * _S155.x + cx_1, fy_1 * _S155.y + cy_1);
            break;
        }
        float2  uv1_2;
        bool all_valid_2 = true & (!_S119);
        for(;;)
        {
            float2  _S156 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_3 = length_0(_S156);
            float _S157 = vert1_c_1.z;
            float theta_1 = (F32_atan2((r_3), (_S157)));
            if(theta_1 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S157;
            }
            else
            {
                k_0 = theta_1 / r_3;
            }
            float2  uv1_3 = _S156 * make_float2 (k_0);
            float u_7 = uv1_3.x;
            float v_7 = uv1_3.y;
            float _S158 = 0.0f * v_7;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float s_diff_r2_8 = u_7 + u_7 + (_S158 + _S158);
            float _S159 = _S106 + r2_7 * _S107;
            float _S160 = _S105 + r2_7 * _S159;
            float _S161 = _S104 + r2_7 * _S160;
            float radial_4 = 1.0f + r2_7 * _S161;
            float _S162 = _S114 * u_7;
            float _S163 = 2.0f * u_7;
            float _S164 = _S115 * u_7;
            float _S165 = 2.0f * v_7;
            float2  _S166 = _S103 * make_float2 (radial_4) + make_float2 (s_diff_r2_8 * _S161 + (s_diff_r2_8 * _S160 + (s_diff_r2_8 * _S159 + s_diff_r2_8 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S114 * v_7 + 0.0f * _S162 + (s_diff_r2_8 + (_S163 + _S163)) * _S109 + s_diff_r2_8 * _S110, _S115 * v_7 + 0.0f * _S164 + (s_diff_r2_8 + (_S158 + 0.0f * _S165)) * _S108 + s_diff_r2_8 * _S111);
            float _S167 = 0.0f * u_7;
            float s_diff_r2_9 = _S167 + _S167 + (v_7 + v_7);
            float2  _S168 = _S116 * make_float2 (radial_4) + make_float2 (s_diff_r2_9 * _S161 + (s_diff_r2_9 * _S160 + (s_diff_r2_9 * _S159 + s_diff_r2_9 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S117 * v_7 + _S162 + (s_diff_r2_9 + (_S167 + 0.0f * _S163)) * _S109 + s_diff_r2_9 * _S110, _S118 * v_7 + _S164 + (s_diff_r2_9 + (_S165 + _S165)) * _S108 + s_diff_r2_9 * _S111);
            Matrix<float, 2, 2>  _S169 = transpose_0(makeMatrix<float, 2, 2> (_S166 + make_float2 (_S166.x * _S112 + _S166.y * _S113, 0.0f), _S168 + make_float2 (_S168.x * _S112 + _S168.y * _S113, 0.0f)));
            bool _S170 = !((F32_min((determinant_0(_S169)), ((F32_min((_S169.rows[int(0)].x), (_S169.rows[int(1)].y)))))) > 0.0f);
            _S120 = _S170;
            if(_S170)
            {
                uv1_2 = uv1_3;
                break;
            }
            float2  _S171 = uv1_3 * make_float2 (radial_4) + make_float2 (_S162 * v_7 + _S109 * (r2_7 + _S163 * u_7) + _S110 * r2_7, _S164 * v_7 + _S108 * (r2_7 + _S165 * v_7) + _S111 * r2_7);
            float2  _S172 = _S171 + make_float2 (_S112 * _S171.x + _S113 * _S171.y, 0.0f);
            uv1_2 = make_float2 (fx_1 * _S172.x + cx_1, fy_1 * _S172.y + cy_1);
            break;
        }
        float2  uv2_2;
        bool all_valid_3 = all_valid_2 & (!_S120);
        for(;;)
        {
            float2  _S173 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_4 = length_0(_S173);
            float _S174 = vert2_c_1.z;
            float theta_2 = (F32_atan2((r_4), (_S174)));
            if(theta_2 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S174;
            }
            else
            {
                k_0 = theta_2 / r_4;
            }
            float2  uv2_3 = _S173 * make_float2 (k_0);
            float u_8 = uv2_3.x;
            float v_8 = uv2_3.y;
            float _S175 = 0.0f * v_8;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float s_diff_r2_10 = u_8 + u_8 + (_S175 + _S175);
            float _S176 = _S106 + r2_8 * _S107;
            float _S177 = _S105 + r2_8 * _S176;
            float _S178 = _S104 + r2_8 * _S177;
            float radial_5 = 1.0f + r2_8 * _S178;
            float _S179 = _S114 * u_8;
            float _S180 = 2.0f * u_8;
            float _S181 = _S115 * u_8;
            float _S182 = 2.0f * v_8;
            float2  _S183 = _S103 * make_float2 (radial_5) + make_float2 (s_diff_r2_10 * _S178 + (s_diff_r2_10 * _S177 + (s_diff_r2_10 * _S176 + s_diff_r2_10 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S114 * v_8 + 0.0f * _S179 + (s_diff_r2_10 + (_S180 + _S180)) * _S109 + s_diff_r2_10 * _S110, _S115 * v_8 + 0.0f * _S181 + (s_diff_r2_10 + (_S175 + 0.0f * _S182)) * _S108 + s_diff_r2_10 * _S111);
            float _S184 = 0.0f * u_8;
            float s_diff_r2_11 = _S184 + _S184 + (v_8 + v_8);
            float2  _S185 = _S116 * make_float2 (radial_5) + make_float2 (s_diff_r2_11 * _S178 + (s_diff_r2_11 * _S177 + (s_diff_r2_11 * _S176 + s_diff_r2_11 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S117 * v_8 + _S179 + (s_diff_r2_11 + (_S184 + 0.0f * _S180)) * _S109 + s_diff_r2_11 * _S110, _S118 * v_8 + _S181 + (s_diff_r2_11 + (_S182 + _S182)) * _S108 + s_diff_r2_11 * _S111);
            Matrix<float, 2, 2>  _S186 = transpose_0(makeMatrix<float, 2, 2> (_S183 + make_float2 (_S183.x * _S112 + _S183.y * _S113, 0.0f), _S185 + make_float2 (_S185.x * _S112 + _S185.y * _S113, 0.0f)));
            bool _S187 = !((F32_min((determinant_0(_S186)), ((F32_min((_S186.rows[int(0)].x), (_S186.rows[int(1)].y)))))) > 0.0f);
            _S121 = _S187;
            if(_S187)
            {
                uv2_2 = uv2_3;
                break;
            }
            float2  _S188 = uv2_3 * make_float2 (radial_5) + make_float2 (_S179 * v_8 + _S109 * (r2_8 + _S180 * u_8) + _S110 * r2_8, _S181 * v_8 + _S108 * (r2_8 + _S182 * v_8) + _S111 * r2_8);
            float2  _S189 = _S188 + make_float2 (_S112 * _S188.x + _S113 * _S188.y, 0.0f);
            uv2_2 = make_float2 (fx_1 * _S189.x + cx_1, fy_1 * _S189.y + cy_1);
            break;
        }
        if(!(all_valid_3 & (!_S121)))
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_1 = uv1_2 - uv0_2;
        float2  e1_1 = uv2_2 - uv1_2;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(uv0_2 - uv2_2)));
        float _S190 = uv0_2.x;
        float _S191 = uv1_2.x;
        float _S192 = uv2_2.x;
        float xmax_1 = (F32_max(((F32_max((_S190), (_S191)))), (_S192))) + offset_1;
        float xmin_1 = (F32_min(((F32_min((_S190), (_S191)))), (_S192))) - offset_1;
        float _S193 = uv0_2.y;
        float _S194 = uv1_2.y;
        float _S195 = uv2_2.y;
        float ymax_1 = (F32_max(((F32_max((_S193), (_S194)))), (_S195))) + offset_1;
        float ymin_1 = (F32_min(((F32_min((_S193), (_S194)))), (_S195))) - offset_1;
        if(xmax_1 <= 0.0f)
        {
            _S123 = true;
        }
        else
        {
            _S123 = xmin_1 >= float(image_width_1);
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = ymax_1 <= 0.0f;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = ymin_1 >= float(image_height_1);
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            if((mean_c_1.z) <= 0.0f)
            {
                if(xmin_1 <= 0.0f)
                {
                    _S123 = xmax_1 >= float(image_width_1);
                }
                else
                {
                    _S123 = false;
                }
                if(_S123)
                {
                    _S123 = true;
                }
                else
                {
                    if(ymin_1 <= 0.0f)
                    {
                        _S123 = ymax_1 >= float(image_width_1);
                    }
                    else
                    {
                        _S123 = false;
                    }
                }
            }
            else
            {
                _S123 = false;
            }
        }
        if(_S123)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int((F32_floor((xmin_1)))), int((F32_floor((ymin_1)))), int((F32_ceil((xmax_1)))), int((F32_ceil((ymax_1)))));
        float3  _S196 = (vert0_c_1 + vert1_c_1 + vert2_c_1) / make_float3 (3.0f);
        float _S197 = _S196.z;
        *depth_1 = (F32_max((_S197), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S197 + length_1(_S196)))));
        float3  _S198 = mean_1 - - mul_0(transpose_1(R_1), t_1);
        float _S199 = _S198.x;
        float _S200 = _S198.y;
        float _S201 = _S198.z;
        float norm_1 = (F32_sqrt((_S199 * _S199 + _S200 * _S200 + _S201 * _S201)));
        float x_11 = _S199 / norm_1;
        float y_4 = _S200 / norm_1;
        float z_1 = _S201 / norm_1;
        float z2_3 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_11 * x_11 - y_4 * y_4;
        float fS1_1 = 2.0f * x_11 * y_4;
        float fTmp0C_1 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        float3  color_1 = make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * sh_coeffs_1[int(1)] + make_float3 (z_1) * sh_coeffs_1[int(2)] - make_float3 (x_11) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_4) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_11) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_11 * fS1_1 + y_4 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_4) * sh_coeffs_1[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_11) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_11 * fC1_1 - y_4 * fS1_1)) * sh_coeffs_1[int(15)]);
        float3  _S202 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_1 + ch_coeffs_1[int(0)] + make_float3 (0.5f), _S202);
        float3  _S203 = color_1 - ch_coeffs_1[int(0)] * make_float3 (0.5f);
        float3  _S204 = ch_coeffs_1[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S203 + _S204 + make_float3 (0.5f), _S202);
        (*rgbs_1)[int(2)] = max_0(_S203 - _S204 + make_float3 (0.5f), _S202);
        (*verts_1)[int(0)] = vert0_1;
        (*verts_1)[int(1)] = vert1_1;
        (*verts_1)[int(2)] = vert2_1;
        float3  _S205 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S205 * make_float3 (float(- (F32_sign((dot_0(_S205, mean_c_1))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
{
    float3  mean_c_2 = mul_0(R_2, mean_2) + t_2;
    float _S206 = scale_2.x;
    float sx_2 = (F32_exp((_S206)));
    float _S207 = scale_2.y;
    float sy_2 = (F32_exp((_S207)));
    float sz_2 = scale_2.z - 0.5f * (_S206 + _S207);
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
    Matrix<float, 3, 3>  _S208 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_4), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_4), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)));
    float3  vert0_2 = mul_0(_S208, make_float3 (sx_2, 0.0f, 0.0f)) + mean_2;
    float3  vert1_2 = mul_0(_S208, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_2;
    float3  vert2_2 = mul_0(_S208, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_2;
    float3  vert0_c_2 = mul_0(R_2, vert0_2) + t_2;
    float3  vert1_c_2 = mul_0(R_2, vert1_2) + t_2;
    float3  vert2_c_2 = mul_0(R_2, vert2_2) + t_2;
    float2  _S209 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_9 = _S209.x;
    float v_9 = _S209.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S210 = 2.0f * dist_coeffs_2[int(4)];
    float _S211 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S212 = _S209 * make_float2 (1.0f + r2_9 * (dist_coeffs_2[int(0)] + r2_9 * (dist_coeffs_2[int(1)] + r2_9 * (dist_coeffs_2[int(2)] + r2_9 * dist_coeffs_2[int(3)])))) + make_float2 (_S210 * u_9 * v_9 + dist_coeffs_2[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_2[int(6)] * r2_9, _S211 * u_9 * v_9 + dist_coeffs_2[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_2[int(7)] * r2_9);
    float2  _S213 = _S212 + make_float2 (dist_coeffs_2[int(8)] * _S212.x + dist_coeffs_2[int(9)] * _S212.y, 0.0f);
    float _S214 = fx_2 * _S213.x + cx_2;
    float _S215 = fy_2 * _S213.y + cy_2;
    float2  uv0_4 = make_float2 (_S214, _S215);
    float2  _S216 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_10 = _S216.x;
    float v_10 = _S216.y;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float2  _S217 = _S216 * make_float2 (1.0f + r2_10 * (dist_coeffs_2[int(0)] + r2_10 * (dist_coeffs_2[int(1)] + r2_10 * (dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)])))) + make_float2 (_S210 * u_10 * v_10 + dist_coeffs_2[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + dist_coeffs_2[int(6)] * r2_10, _S211 * u_10 * v_10 + dist_coeffs_2[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + dist_coeffs_2[int(7)] * r2_10);
    float2  _S218 = _S217 + make_float2 (dist_coeffs_2[int(8)] * _S217.x + dist_coeffs_2[int(9)] * _S217.y, 0.0f);
    float _S219 = fx_2 * _S218.x + cx_2;
    float _S220 = fy_2 * _S218.y + cy_2;
    float2  uv1_4 = make_float2 (_S219, _S220);
    float2  _S221 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_11 = _S221.x;
    float v_11 = _S221.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float2  _S222 = _S221 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S210 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S211 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
    float2  _S223 = _S222 + make_float2 (dist_coeffs_2[int(8)] * _S222.x + dist_coeffs_2[int(9)] * _S222.y, 0.0f);
    float _S224 = fx_2 * _S223.x + cx_2;
    float _S225 = fy_2 * _S223.y + cy_2;
    float2  uv2_4 = make_float2 (_S224, _S225);
    float2  e0_2 = uv1_4 - uv0_4;
    float2  e1_2 = uv2_4 - uv1_4;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(uv0_4 - uv2_4)));
    *aabb_xyxy_2 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S214), (_S219)))), (_S224))) - offset_2)))), int((F32_floor(((F32_min(((F32_min((_S215), (_S220)))), (_S225))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max((_S214), (_S219)))), (_S224))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max((_S215), (_S220)))), (_S225))) + offset_2)))));
    float3  _S226 = (vert0_c_2 + vert1_c_2 + vert2_c_2) / make_float3 (3.0f);
    float _S227 = _S226.z;
    *depth_2 = (F32_max((_S227), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S227 + length_1(_S226)))));
    float3  _S228 = mean_2 - - mul_0(transpose_1(R_2), t_2);
    float _S229 = _S228.x;
    float _S230 = _S228.y;
    float _S231 = _S228.z;
    float norm_2 = (F32_sqrt((_S229 * _S229 + _S230 * _S230 + _S231 * _S231)));
    float x_13 = _S229 / norm_2;
    float y_5 = _S230 / norm_2;
    float z_2 = _S231 / norm_2;
    float z2_5 = z_2 * z_2;
    float fTmp0B_2 = -1.09254848957061768f * z_2;
    float fC1_2 = x_13 * x_13 - y_5 * y_5;
    float fS1_2 = 2.0f * x_13 * y_5;
    float fTmp0C_2 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_2;
    float3  color_2 = make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_2[int(1)] + make_float3 (z_2) * sh_coeffs_2[int(2)] - make_float3 (x_13) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_5) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_13) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_13 * fS1_2 + y_5 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_5) * sh_coeffs_2[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_13) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_13 * fC1_2 - y_5 * fS1_2)) * sh_coeffs_2[int(15)]);
    float3  _S232 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_2 + ch_coeffs_2[int(0)] + make_float3 (0.5f), _S232);
    float3  _S233 = color_2 - ch_coeffs_2[int(0)] * make_float3 (0.5f);
    float3  _S234 = ch_coeffs_2[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S233 + _S234 + make_float3 (0.5f), _S232);
    (*rgbs_2)[int(2)] = max_0(_S233 - _S234 + make_float3 (0.5f), _S232);
    (*verts_2)[int(0)] = vert0_2;
    (*verts_2)[int(1)] = vert1_2;
    (*verts_2)[int(2)] = vert2_2;
    float3  _S235 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S235 * make_float3 (float(- (F32_sign((dot_0(_S235, mean_c_2))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
{
    float3  mean_c_3 = mul_0(R_3, mean_3) + t_3;
    float _S236 = scale_3.x;
    float sx_3 = (F32_exp((_S236)));
    float _S237 = scale_3.y;
    float sy_3 = (F32_exp((_S237)));
    float sz_3 = scale_3.z - 0.5f * (_S236 + _S237);
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
    Matrix<float, 3, 3>  _S238 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_6), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_6), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    float3  vert0_3 = mul_0(_S238, make_float3 (sx_3, 0.0f, 0.0f)) + mean_3;
    float3  vert1_3 = mul_0(_S238, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_3;
    float3  vert2_3 = mul_0(_S238, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_3;
    float3  vert0_c_3 = mul_0(R_3, vert0_3) + t_3;
    float3  vert1_c_3 = mul_0(R_3, vert1_3) + t_3;
    float3  vert2_c_3 = mul_0(R_3, vert2_3) + t_3;
    float2  _S239 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_5 = length_0(_S239);
    float _S240 = vert0_c_3.z;
    float theta_3 = (F32_atan2((r_5), (_S240)));
    float k_1;
    if(theta_3 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_3 * theta_3 / 3.0f) / _S240;
    }
    else
    {
        k_1 = theta_3 / r_5;
    }
    float2  _S241 = _S239 * make_float2 (k_1);
    float u_12 = _S241.x;
    float v_12 = _S241.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S242 = 2.0f * dist_coeffs_3[int(4)];
    float _S243 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S244 = _S241 * make_float2 (1.0f + r2_12 * (dist_coeffs_3[int(0)] + r2_12 * (dist_coeffs_3[int(1)] + r2_12 * (dist_coeffs_3[int(2)] + r2_12 * dist_coeffs_3[int(3)])))) + make_float2 (_S242 * u_12 * v_12 + dist_coeffs_3[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + dist_coeffs_3[int(6)] * r2_12, _S243 * u_12 * v_12 + dist_coeffs_3[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + dist_coeffs_3[int(7)] * r2_12);
    float2  _S245 = _S244 + make_float2 (dist_coeffs_3[int(8)] * _S244.x + dist_coeffs_3[int(9)] * _S244.y, 0.0f);
    float _S246 = fx_3 * _S245.x + cx_3;
    float _S247 = fy_3 * _S245.y + cy_3;
    float2  uv0_5 = make_float2 (_S246, _S247);
    float2  _S248 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_6 = length_0(_S248);
    float _S249 = vert1_c_3.z;
    float theta_4 = (F32_atan2((r_6), (_S249)));
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S249;
    }
    else
    {
        k_1 = theta_4 / r_6;
    }
    float2  _S250 = _S248 * make_float2 (k_1);
    float u_13 = _S250.x;
    float v_13 = _S250.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S251 = _S250 * make_float2 (1.0f + r2_13 * (dist_coeffs_3[int(0)] + r2_13 * (dist_coeffs_3[int(1)] + r2_13 * (dist_coeffs_3[int(2)] + r2_13 * dist_coeffs_3[int(3)])))) + make_float2 (_S242 * u_13 * v_13 + dist_coeffs_3[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_3[int(6)] * r2_13, _S243 * u_13 * v_13 + dist_coeffs_3[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_3[int(7)] * r2_13);
    float2  _S252 = _S251 + make_float2 (dist_coeffs_3[int(8)] * _S251.x + dist_coeffs_3[int(9)] * _S251.y, 0.0f);
    float _S253 = fx_3 * _S252.x + cx_3;
    float _S254 = fy_3 * _S252.y + cy_3;
    float2  uv1_5 = make_float2 (_S253, _S254);
    float2  _S255 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_7 = length_0(_S255);
    float _S256 = vert2_c_3.z;
    float theta_5 = (F32_atan2((r_7), (_S256)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S256;
    }
    else
    {
        k_1 = theta_5 / r_7;
    }
    float2  _S257 = _S255 * make_float2 (k_1);
    float u_14 = _S257.x;
    float v_14 = _S257.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S258 = _S257 * make_float2 (1.0f + r2_14 * (dist_coeffs_3[int(0)] + r2_14 * (dist_coeffs_3[int(1)] + r2_14 * (dist_coeffs_3[int(2)] + r2_14 * dist_coeffs_3[int(3)])))) + make_float2 (_S242 * u_14 * v_14 + dist_coeffs_3[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + dist_coeffs_3[int(6)] * r2_14, _S243 * u_14 * v_14 + dist_coeffs_3[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + dist_coeffs_3[int(7)] * r2_14);
    float2  _S259 = _S258 + make_float2 (dist_coeffs_3[int(8)] * _S258.x + dist_coeffs_3[int(9)] * _S258.y, 0.0f);
    float _S260 = fx_3 * _S259.x + cx_3;
    float _S261 = fy_3 * _S259.y + cy_3;
    float2  uv2_5 = make_float2 (_S260, _S261);
    float2  e0_3 = uv1_5 - uv0_5;
    float2  e1_3 = uv2_5 - uv1_5;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(uv0_5 - uv2_5)));
    *aabb_xyxy_3 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S246), (_S253)))), (_S260))) - offset_3)))), int((F32_floor(((F32_min(((F32_min((_S247), (_S254)))), (_S261))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max((_S246), (_S253)))), (_S260))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max((_S247), (_S254)))), (_S261))) + offset_3)))));
    float3  _S262 = (vert0_c_3 + vert1_c_3 + vert2_c_3) / make_float3 (3.0f);
    float _S263 = _S262.z;
    *depth_3 = (F32_max((_S263), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S263 + length_1(_S262)))));
    float3  _S264 = mean_3 - - mul_0(transpose_1(R_3), t_3);
    float _S265 = _S264.x;
    float _S266 = _S264.y;
    float _S267 = _S264.z;
    float norm_3 = (F32_sqrt((_S265 * _S265 + _S266 * _S266 + _S267 * _S267)));
    float x_15 = _S265 / norm_3;
    float y_6 = _S266 / norm_3;
    float z_3 = _S267 / norm_3;
    float z2_7 = z_3 * z_3;
    float fTmp0B_3 = -1.09254848957061768f * z_3;
    float fC1_3 = x_15 * x_15 - y_6 * y_6;
    float fS1_3 = 2.0f * x_15 * y_6;
    float fTmp0C_3 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_3;
    float3  color_3 = make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_3[int(1)] + make_float3 (z_3) * sh_coeffs_3[int(2)] - make_float3 (x_15) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_6) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_15) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_3 + y_6 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_6) * sh_coeffs_3[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_15) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_3 - y_6 * fS1_3)) * sh_coeffs_3[int(15)]);
    float3  _S268 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_3 + ch_coeffs_3[int(0)] + make_float3 (0.5f), _S268);
    float3  _S269 = color_3 - ch_coeffs_3[int(0)] * make_float3 (0.5f);
    float3  _S270 = ch_coeffs_3[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S269 + _S270 + make_float3 (0.5f), _S268);
    (*rgbs_3)[int(2)] = max_0(_S269 - _S270 + make_float3 (0.5f), _S268);
    (*verts_3)[int(0)] = vert0_3;
    (*verts_3)[int(1)] = vert1_3;
    (*verts_3)[int(2)] = vert2_3;
    float3  _S271 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S271 * make_float3 (float(- (F32_sign((dot_0(_S271, mean_c_3))))));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S272, float3  _S273)
{
    return mul_0(_S272, _S273);
}

inline __device__ float s_primal_ctx_exp_0(float _S274)
{
    return (F32_exp((_S274)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S275)
{
    return (F32_sqrt((_S275)));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S276, float3  _S277)
{
    return cross_0(_S276, _S277);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S278, float3  _S279)
{
    return dot_0(_S278, _S279);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S280, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S281, float _S282)
{
    _d_dot_0(_S280, _S281, _S282);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S283, float _S284)
{
    _d_sqrt_0(_S283, _S284);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S285 = (*dpx_9).primal_0.x;
    float _S286 = (*dpx_9).primal_0.y;
    float _S287 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = _S285 * _S285 + _S286 * _S286 + _S287 * _S287;
    (&_S288)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S288, _s_dOut_0);
    float _S289 = (*dpx_9).primal_0.z * _S288.differential_0;
    float _S290 = _S289 + _S289;
    float _S291 = (*dpx_9).primal_0.y * _S288.differential_0;
    float _S292 = _S291 + _S291;
    float _S293 = (*dpx_9).primal_0.x * _S288.differential_0;
    float _S294 = _S293 + _S293;
    float3  _S295 = make_float3 (0.0f);
    *&((&_S295)->z) = _S290;
    *&((&_S295)->y) = _S292;
    *&((&_S295)->x) = _S294;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S295;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S296, float _S297)
{
    s_bwd_prop_length_impl_0(_S296, _S297);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_1)
{
    float _S298 = length_1((*dpx_10).primal_0);
    float3  _S299 = (*dpx_10).primal_0 * _s_dOut_1;
    float3  _S300 = make_float3 (1.0f / _S298) * _s_dOut_1;
    float _S301 = - ((_S299.x + _S299.y + _S299.z) / (_S298 * _S298));
    float3  _S302 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S303;
    (&_S303)->primal_0 = (*dpx_10).primal_0;
    (&_S303)->differential_0 = _S302;
    s_bwd_length_impl_0(&_S303, _S301);
    float3  _S304 = _S300 + _S303.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S304;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S305, float3  _S306)
{
    s_bwd_prop_normalize_impl_0(_S305, _S306);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S307, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S308, float3  _S309)
{
    _d_cross_0(_S307, _S308, _S309);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S310, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S311, float3  _S312)
{
    _d_max_vector_0(_S310, _S311, _S312);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S313, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S314, float3  _S315)
{
    _d_mul_0(_S313, _S314, _S315);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S316, float _S317)
{
    _d_exp2_0(_S316, _S317);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float _s_dOut_2)
{
    float _S318 = (*dpx_11).primal_0.x;
    float _S319 = (*dpx_11).primal_0.y;
    DiffPair_float_0 _S320;
    (&_S320)->primal_0 = _S318 * _S318 + _S319 * _S319;
    (&_S320)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S320, _s_dOut_2);
    float _S321 = (*dpx_11).primal_0.y * _S320.differential_0;
    float _S322 = _S321 + _S321;
    float _S323 = (*dpx_11).primal_0.x * _S320.differential_0;
    float _S324 = _S323 + _S323;
    float2  _S325 = make_float2 (0.0f);
    *&((&_S325)->y) = _S322;
    *&((&_S325)->x) = _S324;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S325;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S326, float _S327)
{
    s_bwd_prop_length_impl_1(_S326, _S327);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S328, float _S329)
{
    _d_abs_0(_S328, _S329);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S330, float _S331)
{
    _d_exp_0(_S330, _S331);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S332 = scale_4.x;
    float _S333 = s_primal_ctx_exp_0(_S332);
    float _S334 = scale_4.y;
    float _S335 = s_primal_ctx_exp_0(_S334);
    float sz_4 = scale_4.z - 0.5f * (_S332 + _S334);
    float _S336 = quat_4.y;
    float x2_4 = _S336 * _S336;
    float y2_4 = quat_4.z * quat_4.z;
    float z2_8 = quat_4.w * quat_4.w;
    float xy_4 = quat_4.y * quat_4.z;
    float xz_4 = quat_4.y * quat_4.w;
    float yz_4 = quat_4.z * quat_4.w;
    float wx_4 = quat_4.x * quat_4.y;
    float wy_4 = quat_4.x * quat_4.z;
    float wz_4 = quat_4.x * quat_4.w;
    Matrix<float, 3, 3>  _S337 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S338 = make_float3 (_S333, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S337, _S338) + mean_4;
    float _S339 = -0.5f + sz_4;
    float3  _S340 = make_float3 (_S333 * _S339, _S335, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S337, _S340) + mean_4;
    float _S341 = -0.5f - sz_4;
    float3  _S342 = make_float3 (_S333 * _S341, - _S335, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S337, _S342) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S343 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S344 = vert0_c_4.z;
    float2  _S345 = make_float2 (_S344);
    float2  _S346 = _S343 / make_float2 (_S344);
    float2  _S347 = make_float2 (_S344 * _S344);
    float u_15 = _S346.x;
    float v_15 = _S346.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S348 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S349 = dist_coeffs_4[int(1)] + r2_15 * _S348;
    float _S350 = dist_coeffs_4[int(0)] + r2_15 * _S349;
    float radial_6 = 1.0f + r2_15 * _S350;
    float _S351 = 2.0f * dist_coeffs_4[int(4)];
    float _S352 = _S351 * u_15;
    float _S353 = 2.0f * u_15;
    float _S354 = 2.0f * dist_coeffs_4[int(5)];
    float _S355 = _S354 * u_15;
    float _S356 = 2.0f * v_15;
    float2  _S357 = _S346 * make_float2 (radial_6) + make_float2 (_S352 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S353 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S355 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S356 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S358 = _S357 + make_float2 (dist_coeffs_4[int(8)] * _S357.x + dist_coeffs_4[int(9)] * _S357.y, 0.0f);
    float _S359 = fx_4 * _S358.x + cx_4;
    float _S360 = fy_4 * _S358.y + cy_4;
    float2  uv0_6 = make_float2 (_S359, _S360);
    float2  _S361 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S362 = vert1_c_4.z;
    float2  _S363 = make_float2 (_S362);
    float2  _S364 = _S361 / make_float2 (_S362);
    float2  _S365 = make_float2 (_S362 * _S362);
    float u_16 = _S364.x;
    float v_16 = _S364.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S366 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S367 = dist_coeffs_4[int(1)] + r2_16 * _S366;
    float _S368 = dist_coeffs_4[int(0)] + r2_16 * _S367;
    float radial_7 = 1.0f + r2_16 * _S368;
    float _S369 = _S351 * u_16;
    float _S370 = 2.0f * u_16;
    float _S371 = _S354 * u_16;
    float _S372 = 2.0f * v_16;
    float2  _S373 = _S364 * make_float2 (radial_7) + make_float2 (_S369 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S370 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S371 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S372 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S374 = _S373 + make_float2 (dist_coeffs_4[int(8)] * _S373.x + dist_coeffs_4[int(9)] * _S373.y, 0.0f);
    float _S375 = fx_4 * _S374.x + cx_4;
    float _S376 = fy_4 * _S374.y + cy_4;
    float2  uv1_6 = make_float2 (_S375, _S376);
    float2  _S377 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S378 = vert2_c_4.z;
    float2  _S379 = make_float2 (_S378);
    float2  _S380 = _S377 / make_float2 (_S378);
    float2  _S381 = make_float2 (_S378 * _S378);
    float u_17 = _S380.x;
    float v_17 = _S380.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S382 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S383 = dist_coeffs_4[int(1)] + r2_17 * _S382;
    float _S384 = dist_coeffs_4[int(0)] + r2_17 * _S383;
    float radial_8 = 1.0f + r2_17 * _S384;
    float _S385 = _S351 * u_17;
    float _S386 = 2.0f * u_17;
    float _S387 = _S354 * u_17;
    float _S388 = 2.0f * v_17;
    float2  _S389 = _S380 * make_float2 (radial_8) + make_float2 (_S385 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S386 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S387 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S388 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S390 = _S389 + make_float2 (dist_coeffs_4[int(8)] * _S389.x + dist_coeffs_4[int(9)] * _S389.y, 0.0f);
    float _S391 = fx_4 * _S390.x + cx_4;
    float _S392 = fy_4 * _S390.y + cy_4;
    float2  uv2_6 = make_float2 (_S391, _S392);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S393 = e0_4.x;
    float _S394 = e1_4.y;
    float _S395 = e0_4.y;
    float _S396 = e1_4.x;
    float _S397 = _S393 * _S394 - _S395 * _S396;
    float _S398 = 1.0f - hardness_4.y;
    float _S399 = -1.0f / _S398;
    float _S400 = _S398 * _S398;
    float _S401 = (F32_max((_S359), (_S375)));
    float _S402 = (F32_min((_S359), (_S375)));
    float _S403 = (F32_max((_S360), (_S376)));
    float _S404 = (F32_min((_S360), (_S376)));
    float3  _S405 = (vert0_c_4 + vert1_c_4 + vert2_c_4) / make_float3 (3.0f);
    float _S406 = _S405.z;
    float _S407 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S408 = _S407 * (_S406 + length_1(_S405));
    Matrix<float, 3, 3>  _S409 = transpose_1(R_4);
    float3  _S410 = mean_4 - - s_primal_ctx_mul_0(_S409, t_4);
    float _S411 = _S410.x;
    float _S412 = _S410.y;
    float _S413 = _S410.z;
    float _S414 = _S411 * _S411 + _S412 * _S412 + _S413 * _S413;
    float _S415 = s_primal_ctx_sqrt_0(_S414);
    float x_16 = _S411 / _S415;
    float3  _S416 = make_float3 (x_16);
    float _S417 = _S415 * _S415;
    float y_7 = _S412 / _S415;
    float z_4 = _S413 / _S415;
    float3  _S418 = make_float3 (z_4);
    float _S419 = - y_7;
    float3  _S420 = make_float3 (_S419);
    float z2_9 = z_4 * z_4;
    float fTmp0B_4 = -1.09254848957061768f * z_4;
    float fC1_4 = x_16 * x_16 - y_7 * y_7;
    float _S421 = 2.0f * x_16;
    float fS1_4 = _S421 * y_7;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S422 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_16;
    float3  _S423 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_7;
    float3  _S424 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S425 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S426 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_4;
    float _S427 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_4 * _S427;
    float3  _S428 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_16;
    float3  _S429 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_7;
    float3  _S430 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S431 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S432 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_16 * fC1_4 - y_7 * fS1_4);
    float3  _S433 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_16 * fS1_4 + y_7 * fC1_4);
    float3  _S434 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S419) * sh_coeffs_4[int(1)] + make_float3 (z_4) * sh_coeffs_4[int(2)] - make_float3 (x_16) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]);
    float3  _S435 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S436 = make_float3 (0.0f);
    float3  _S437 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S438 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S439 = make_float3 (_S438);
    float3  _S440 = ch_coeffs_4[int(1)] * make_float3 (_S438);
    float3  _S441 = _S437 + _S440 + make_float3 (0.5f);
    float3  _S442 = _S437 - _S440 + make_float3 (0.5f);
    float3  _S443 = vert1_c_4 - vert0_c_4;
    float3  _S444 = vert2_c_4 - vert0_c_4;
    float3  _S445 = s_primal_ctx_cross_0(_S443, _S444);
    float3  _S446 = normalize_0(_S445);
    float3  _S447 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S446, mean_c_4)))))) * v_normal_0;
    float3  _S448 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S449;
    (&_S449)->primal_0 = _S446;
    (&_S449)->differential_0 = _S448;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S450;
    (&_S450)->primal_0 = mean_c_4;
    (&_S450)->differential_0 = _S448;
    s_bwd_prop_dot_0(&_S449, &_S450, 0.0f);
    float3  _S451 = _S447 + _S449.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S452;
    (&_S452)->primal_0 = _S445;
    (&_S452)->differential_0 = _S448;
    s_bwd_normalize_impl_0(&_S452, _S451);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S453;
    (&_S453)->primal_0 = _S443;
    (&_S453)->differential_0 = _S448;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S454;
    (&_S454)->primal_0 = _S444;
    (&_S454)->differential_0 = _S448;
    s_bwd_prop_cross_0(&_S453, &_S454, _S452.differential_0);
    float3  _S455 = - _S454.differential_0;
    float3  _S456 = - _S453.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S457;
    (&_S457)->primal_0 = _S442;
    (&_S457)->differential_0 = _S448;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S458;
    (&_S458)->primal_0 = _S436;
    (&_S458)->differential_0 = _S448;
    s_bwd_prop_max_0(&_S457, &_S458, v_rgbs_0[int(2)]);
    float3  _S459 = - _S457.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S460;
    (&_S460)->primal_0 = _S441;
    (&_S460)->differential_0 = _S448;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S461;
    (&_S461)->primal_0 = _S436;
    (&_S461)->differential_0 = _S448;
    s_bwd_prop_max_0(&_S460, &_S461, v_rgbs_0[int(1)]);
    float3  _S462 = _S439 * (_S459 + _S460.differential_0);
    float3  _S463 = _S457.differential_0 + _S460.differential_0;
    float3  _S464 = make_float3 (0.5f) * - _S463;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S465;
    (&_S465)->primal_0 = _S435;
    (&_S465)->differential_0 = _S448;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S466;
    (&_S466)->primal_0 = _S436;
    (&_S466)->differential_0 = _S448;
    s_bwd_prop_max_0(&_S465, &_S466, v_rgbs_0[int(0)]);
    float3  _S467 = _S464 + _S465.differential_0;
    float3  _S468 = _S463 + _S465.differential_0;
    float3  _S469 = _S433 * _S468;
    float3  _S470 = sh_coeffs_4[int(15)] * _S468;
    float3  _S471 = _S431 * _S468;
    float3  _S472 = sh_coeffs_4[int(14)] * _S468;
    float3  _S473 = _S429 * _S468;
    float3  _S474 = sh_coeffs_4[int(13)] * _S468;
    float3  _S475 = _S428 * _S468;
    float3  _S476 = sh_coeffs_4[int(12)] * _S468;
    float3  _S477 = _S430 * _S468;
    float3  _S478 = sh_coeffs_4[int(11)] * _S468;
    float3  _S479 = _S432 * _S468;
    float3  _S480 = sh_coeffs_4[int(10)] * _S468;
    float3  _S481 = _S434 * _S468;
    float3  _S482 = sh_coeffs_4[int(9)] * _S468;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S482.x + _S482.y + _S482.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S470.x + _S470.y + _S470.z);
    float _S483 = _S480.x + _S480.y + _S480.z;
    float _S484 = _S472.x + _S472.y + _S472.z;
    float _S485 = _S478.x + _S478.y + _S478.z;
    float _S486 = _S474.x + _S474.y + _S474.z;
    float _S487 = _S476.x + _S476.y + _S476.z;
    float _S488 = - s_diff_fC2_T_0;
    float3  _S489 = _S425 * _S468;
    float3  _S490 = sh_coeffs_4[int(8)] * _S468;
    float3  _S491 = _S423 * _S468;
    float3  _S492 = sh_coeffs_4[int(7)] * _S468;
    float3  _S493 = _S422 * _S468;
    float3  _S494 = sh_coeffs_4[int(6)] * _S468;
    float3  _S495 = _S424 * _S468;
    float3  _S496 = sh_coeffs_4[int(5)] * _S468;
    float3  _S497 = _S426 * _S468;
    float3  _S498 = sh_coeffs_4[int(4)] * _S468;
    float _S499 = _S496.x + _S496.y + _S496.z;
    float _S500 = _S492.x + _S492.y + _S492.z;
    float _S501 = fTmp1B_4 * _S483 + x_16 * s_diff_fS2_T_0 + y_7 * _S488 + 0.54627424478530884f * (_S498.x + _S498.y + _S498.z);
    float _S502 = fTmp1B_4 * _S484 + y_7 * s_diff_fS2_T_0 + x_16 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S490.x + _S490.y + _S490.z);
    float _S503 = y_7 * - _S502;
    float _S504 = x_16 * _S502;
    float _S505 = z_4 * (1.86588168144226074f * (z_4 * _S487) + -2.28522896766662598f * (y_7 * _S485 + x_16 * _S486) + 0.94617468118667603f * (_S494.x + _S494.y + _S494.z));
    float3  _S506 = make_float3 (0.48860251903533936f) * _S468;
    float3  _S507 = - _S506;
    float3  _S508 = _S416 * _S507;
    float3  _S509 = sh_coeffs_4[int(3)] * _S507;
    float3  _S510 = _S418 * _S506;
    float3  _S511 = sh_coeffs_4[int(2)] * _S506;
    float3  _S512 = _S420 * _S506;
    float3  _S513 = sh_coeffs_4[int(1)] * _S506;
    float _S514 = (_S427 * _S487 + 1.44530570507049561f * (fS1_4 * _S483 + fC1_4 * _S484) + -1.09254848957061768f * (y_7 * _S499 + x_16 * _S500) + _S505 + _S505 + _S511.x + _S511.y + _S511.z) / _S417;
    float _S515 = _S415 * _S514;
    float _S516 = (fTmp0C_4 * _S485 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S488 + fTmp0B_4 * _S499 + _S421 * _S501 + _S503 + _S503 + - (_S513.x + _S513.y + _S513.z)) / _S417;
    float _S517 = _S415 * _S516;
    float _S518 = (fTmp0C_4 * _S486 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S500 + 2.0f * (y_7 * _S501) + _S504 + _S504 + _S509.x + _S509.y + _S509.z) / _S417;
    float _S519 = _S415 * _S518;
    float _S520 = _S413 * - _S514 + _S412 * - _S516 + _S411 * - _S518;
    DiffPair_float_0 _S521;
    (&_S521)->primal_0 = _S414;
    (&_S521)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S521, _S520);
    float _S522 = _S413 * _S521.differential_0;
    float _S523 = _S412 * _S521.differential_0;
    float _S524 = _S411 * _S521.differential_0;
    float3  _S525 = make_float3 (0.282094806432724f) * _S468;
    float3  _S526 = make_float3 (_S519 + _S524 + _S524, _S517 + _S523 + _S523, _S515 + _S522 + _S522);
    float3  _S527 = - - _S526;
    Matrix<float, 3, 3>  _S528 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S529;
    (&_S529)->primal_0 = _S409;
    (&_S529)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S530;
    (&_S530)->primal_0 = t_4;
    (&_S530)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S529, &_S530, _S527);
    Matrix<float, 3, 3>  _S531 = transpose_1(_S529.differential_0);
    DiffPair_float_0 _S532;
    (&_S532)->primal_0 = _S406;
    (&_S532)->differential_0 = 0.0f;
    DiffPair_float_0 _S533;
    (&_S533)->primal_0 = _S408;
    (&_S533)->differential_0 = 0.0f;
    _d_max_0(&_S532, &_S533, v_depth_0);
    float _S534 = _S407 * _S533.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S535;
    (&_S535)->primal_0 = _S405;
    (&_S535)->differential_0 = _S448;
    s_bwd_length_impl_0(&_S535, _S534);
    float3  _S536 = make_float3 (0.3333333432674408f) * (_S535.differential_0 + make_float3 (0.0f, 0.0f, _S532.differential_0 + _S534));
    DiffPair_float_0 _S537;
    (&_S537)->primal_0 = _S404;
    (&_S537)->differential_0 = 0.0f;
    DiffPair_float_0 _S538;
    (&_S538)->primal_0 = _S392;
    (&_S538)->differential_0 = 0.0f;
    _d_min_0(&_S537, &_S538, 0.0f);
    DiffPair_float_0 _S539;
    (&_S539)->primal_0 = _S360;
    (&_S539)->differential_0 = 0.0f;
    DiffPair_float_0 _S540;
    (&_S540)->primal_0 = _S376;
    (&_S540)->differential_0 = 0.0f;
    _d_min_0(&_S539, &_S540, _S537.differential_0);
    DiffPair_float_0 _S541;
    (&_S541)->primal_0 = _S403;
    (&_S541)->differential_0 = 0.0f;
    DiffPair_float_0 _S542;
    (&_S542)->primal_0 = _S392;
    (&_S542)->differential_0 = 0.0f;
    _d_max_0(&_S541, &_S542, 0.0f);
    DiffPair_float_0 _S543;
    (&_S543)->primal_0 = _S360;
    (&_S543)->differential_0 = 0.0f;
    DiffPair_float_0 _S544;
    (&_S544)->primal_0 = _S376;
    (&_S544)->differential_0 = 0.0f;
    _d_max_0(&_S543, &_S544, _S541.differential_0);
    DiffPair_float_0 _S545;
    (&_S545)->primal_0 = _S402;
    (&_S545)->differential_0 = 0.0f;
    DiffPair_float_0 _S546;
    (&_S546)->primal_0 = _S391;
    (&_S546)->differential_0 = 0.0f;
    _d_min_0(&_S545, &_S546, 0.0f);
    DiffPair_float_0 _S547;
    (&_S547)->primal_0 = _S359;
    (&_S547)->differential_0 = 0.0f;
    DiffPair_float_0 _S548;
    (&_S548)->primal_0 = _S375;
    (&_S548)->differential_0 = 0.0f;
    _d_min_0(&_S547, &_S548, _S545.differential_0);
    DiffPair_float_0 _S549;
    (&_S549)->primal_0 = _S401;
    (&_S549)->differential_0 = 0.0f;
    DiffPair_float_0 _S550;
    (&_S550)->primal_0 = _S391;
    (&_S550)->differential_0 = 0.0f;
    _d_max_0(&_S549, &_S550, 0.0f);
    DiffPair_float_0 _S551;
    (&_S551)->primal_0 = _S359;
    (&_S551)->differential_0 = 0.0f;
    DiffPair_float_0 _S552;
    (&_S552)->primal_0 = _S375;
    (&_S552)->differential_0 = 0.0f;
    _d_max_0(&_S551, &_S552, _S549.differential_0);
    DiffPair_float_0 _S553;
    (&_S553)->primal_0 = _S399;
    (&_S553)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S553, 0.0f);
    float _S554 = - (-1.0f * - (_S553.differential_0 / _S400));
    float2  _S555 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S556;
    (&_S556)->primal_0 = e2_0;
    (&_S556)->differential_0 = _S555;
    s_bwd_length_impl_1(&_S556, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S557;
    (&_S557)->primal_0 = e1_4;
    (&_S557)->differential_0 = _S555;
    s_bwd_length_impl_1(&_S557, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S558;
    (&_S558)->primal_0 = e0_4;
    (&_S558)->differential_0 = _S555;
    s_bwd_length_impl_1(&_S558, -0.0f);
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = _S397;
    (&_S559)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S559, 0.0f);
    float _S560 = - _S559.differential_0;
    float2  _S561 = _S557.differential_0 + make_float2 (_S395 * _S560, _S393 * _S559.differential_0);
    float2  _S562 = _S558.differential_0 + make_float2 (_S394 * _S559.differential_0, _S396 * _S560);
    float2  _S563 = - _S556.differential_0 + _S561;
    float _S564 = fx_4 * (_S546.differential_0 + _S550.differential_0 + _S563.x);
    float2  _S565 = make_float2 (_S564, fy_4 * (_S538.differential_0 + _S542.differential_0 + _S563.y)) + make_float2 (dist_coeffs_4[int(8)] * _S564, dist_coeffs_4[int(9)] * _S564);
    float2  _S566 = _S380 * _S565;
    float _S567 = dist_coeffs_4[int(4)] * _S565.y;
    float _S568 = dist_coeffs_4[int(5)] * _S565.x;
    float _S569 = _S566.x + _S566.y;
    float _S570 = r2_17 * _S569;
    float _S571 = r2_17 * _S570;
    float _S572 = dist_coeffs_4[int(7)] * _S565.y + _S567 + dist_coeffs_4[int(6)] * _S565.x + _S568 + _S384 * _S569 + _S383 * _S570 + _S382 * _S571 + dist_coeffs_4[int(3)] * (r2_17 * _S571);
    float _S573 = v_17 * _S572;
    float _S574 = u_17 * _S572;
    float2  _S575 = (make_float2 (radial_8) * _S565 + make_float2 (_S354 * (v_17 * _S565.y) + _S386 * _S568 + 2.0f * (u_17 * _S568) + _S351 * (v_17 * _S565.x) + _S574 + _S574, _S388 * _S567 + 2.0f * (v_17 * _S567) + _S387 * _S565.y + _S385 * _S565.x + _S573 + _S573)) / _S381;
    float2  _S576 = _S377 * - _S575;
    float2  _S577 = _S379 * _S575;
    float2  _S578 = - _S561 + _S562;
    float _S579 = fx_4 * (_S548.differential_0 + _S552.differential_0 + _S578.x);
    float2  _S580 = make_float2 (_S579, fy_4 * (_S540.differential_0 + _S544.differential_0 + _S578.y)) + make_float2 (dist_coeffs_4[int(8)] * _S579, dist_coeffs_4[int(9)] * _S579);
    float2  _S581 = _S364 * _S580;
    float _S582 = dist_coeffs_4[int(4)] * _S580.y;
    float _S583 = dist_coeffs_4[int(5)] * _S580.x;
    float _S584 = _S581.x + _S581.y;
    float _S585 = r2_16 * _S584;
    float _S586 = r2_16 * _S585;
    float _S587 = dist_coeffs_4[int(7)] * _S580.y + _S582 + dist_coeffs_4[int(6)] * _S580.x + _S583 + _S368 * _S584 + _S367 * _S585 + _S366 * _S586 + dist_coeffs_4[int(3)] * (r2_16 * _S586);
    float _S588 = v_16 * _S587;
    float _S589 = u_16 * _S587;
    float2  _S590 = (make_float2 (radial_7) * _S580 + make_float2 (_S354 * (v_16 * _S580.y) + _S370 * _S583 + 2.0f * (u_16 * _S583) + _S351 * (v_16 * _S580.x) + _S589 + _S589, _S372 * _S582 + 2.0f * (v_16 * _S582) + _S371 * _S580.y + _S369 * _S580.x + _S588 + _S588)) / _S365;
    float2  _S591 = _S361 * - _S590;
    float2  _S592 = _S363 * _S590;
    float _S593 = _S591.x + _S591.y;
    float2  _S594 = _S556.differential_0 + - _S562;
    float _S595 = fx_4 * (_S547.differential_0 + _S551.differential_0 + _S594.x);
    float2  _S596 = make_float2 (_S595, fy_4 * (_S539.differential_0 + _S543.differential_0 + _S594.y)) + make_float2 (dist_coeffs_4[int(8)] * _S595, dist_coeffs_4[int(9)] * _S595);
    float2  _S597 = _S346 * _S596;
    float _S598 = dist_coeffs_4[int(4)] * _S596.y;
    float _S599 = dist_coeffs_4[int(5)] * _S596.x;
    float _S600 = _S597.x + _S597.y;
    float _S601 = r2_15 * _S600;
    float _S602 = r2_15 * _S601;
    float _S603 = dist_coeffs_4[int(7)] * _S596.y + _S598 + dist_coeffs_4[int(6)] * _S596.x + _S599 + _S350 * _S600 + _S349 * _S601 + _S348 * _S602 + dist_coeffs_4[int(3)] * (r2_15 * _S602);
    float _S604 = v_15 * _S603;
    float _S605 = u_15 * _S603;
    float2  _S606 = (make_float2 (radial_6) * _S596 + make_float2 (_S354 * (v_15 * _S596.y) + _S353 * _S599 + 2.0f * (u_15 * _S599) + _S351 * (v_15 * _S596.x) + _S605 + _S605, _S356 * _S598 + 2.0f * (v_15 * _S598) + _S355 * _S596.y + _S352 * _S596.x + _S604 + _S604)) / _S347;
    float2  _S607 = _S343 * - _S606;
    float2  _S608 = _S345 * _S606;
    float _S609 = _S607.x + _S607.y;
    float3  _S610 = _S454.differential_0 + _S536 + make_float3 (_S577.x, _S577.y, _S576.x + _S576.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S611;
    (&_S611)->primal_0 = R_4;
    (&_S611)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S612;
    (&_S612)->primal_0 = vert2_4;
    (&_S612)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S611, &_S612, _S610);
    float3  _S613 = _S453.differential_0 + _S536 + make_float3 (_S592.x, _S592.y, _S593);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S614;
    (&_S614)->primal_0 = R_4;
    (&_S614)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S615;
    (&_S615)->primal_0 = vert1_4;
    (&_S615)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S614, &_S615, _S613);
    float3  _S616 = _S455 + _S456 + _S536 + make_float3 (_S608.x, _S608.y, _S609);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S617;
    (&_S617)->primal_0 = R_4;
    (&_S617)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S618;
    (&_S618)->primal_0 = vert0_4;
    (&_S618)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S617, &_S618, _S616);
    float3  _S619 = v_verts_0[int(2)] + _S612.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S620;
    (&_S620)->primal_0 = _S337;
    (&_S620)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S621;
    (&_S621)->primal_0 = _S342;
    (&_S621)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S620, &_S621, _S619);
    float _S622 = - _S621.differential_0.y;
    float _S623 = _S341 * _S621.differential_0.x;
    float _S624 = - (_S333 * _S621.differential_0.x);
    float3  _S625 = v_verts_0[int(1)] + _S615.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S626;
    (&_S626)->primal_0 = _S337;
    (&_S626)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S627;
    (&_S627)->primal_0 = _S340;
    (&_S627)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S626, &_S627, _S625);
    float _S628 = _S333 * _S627.differential_0.x;
    float _S629 = _S339 * _S627.differential_0.x;
    float3  _S630 = v_verts_0[int(0)] + _S618.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S631;
    (&_S631)->primal_0 = _S337;
    (&_S631)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
    (&_S632)->primal_0 = _S338;
    (&_S632)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S631, &_S632, _S630);
    Matrix<float, 3, 3>  _S633 = transpose_1(_S620.differential_0 + _S626.differential_0 + _S631.differential_0);
    float _S634 = 2.0f * - _S633.rows[int(2)].z;
    float _S635 = 2.0f * _S633.rows[int(2)].y;
    float _S636 = 2.0f * _S633.rows[int(2)].x;
    float _S637 = 2.0f * _S633.rows[int(1)].z;
    float _S638 = 2.0f * - _S633.rows[int(1)].y;
    float _S639 = 2.0f * _S633.rows[int(1)].x;
    float _S640 = 2.0f * _S633.rows[int(0)].z;
    float _S641 = 2.0f * _S633.rows[int(0)].y;
    float _S642 = 2.0f * - _S633.rows[int(0)].x;
    float _S643 = - _S639 + _S641;
    float _S644 = _S636 + - _S640;
    float _S645 = - _S635 + _S637;
    float _S646 = _S635 + _S637;
    float _S647 = _S636 + _S640;
    float _S648 = _S639 + _S641;
    float _S649 = quat_4.w * (_S638 + _S642);
    float _S650 = quat_4.z * (_S634 + _S642);
    float _S651 = quat_4.y * (_S634 + _S638);
    float _S652 = quat_4.x * _S643 + quat_4.z * _S646 + quat_4.y * _S647 + _S649 + _S649;
    float _S653 = quat_4.x * _S644 + quat_4.w * _S646 + quat_4.y * _S648 + _S650 + _S650;
    float _S654 = quat_4.x * _S645 + quat_4.w * _S647 + quat_4.z * _S648 + _S651 + _S651;
    float _S655 = quat_4.w * _S643 + quat_4.z * _S644 + quat_4.y * _S645;
    float _S656 = _S624 + _S628;
    float _S657 = 0.5f * - _S656;
    float _S658 = _S622 + _S627.differential_0.y;
    DiffPair_float_0 _S659;
    (&_S659)->primal_0 = _S334;
    (&_S659)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S659, _S658);
    float _S660 = _S657 + _S659.differential_0;
    float _S661 = _S623 + _S629 + _S632.differential_0.x;
    DiffPair_float_0 _S662;
    (&_S662)->primal_0 = _S332;
    (&_S662)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S662, _S661);
    float _S663 = _S657 + _S662.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S664;
    (&_S664)->primal_0 = R_4;
    (&_S664)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S665;
    (&_S665)->primal_0 = mean_4;
    (&_S665)->differential_0 = _S448;
    s_bwd_prop_mul_0(&_S664, &_S665, _S450.differential_0);
    float3  _S666 = _S530.differential_0 + _S610 + _S613 + _S616 + _S450.differential_0;
    Matrix<float, 3, 3>  _S667 = _S531 + _S611.differential_0 + _S614.differential_0 + _S617.differential_0 + _S664.differential_0;
    FixedArray<float3 , 2>  _S668;
    _S668[int(0)] = _S448;
    _S668[int(1)] = _S448;
    _S668[int(1)] = _S462;
    _S668[int(0)] = _S467;
    FixedArray<float3 , 16>  _S669;
    _S669[int(0)] = _S448;
    _S669[int(1)] = _S448;
    _S669[int(2)] = _S448;
    _S669[int(3)] = _S448;
    _S669[int(4)] = _S448;
    _S669[int(5)] = _S448;
    _S669[int(6)] = _S448;
    _S669[int(7)] = _S448;
    _S669[int(8)] = _S448;
    _S669[int(9)] = _S448;
    _S669[int(10)] = _S448;
    _S669[int(11)] = _S448;
    _S669[int(12)] = _S448;
    _S669[int(13)] = _S448;
    _S669[int(14)] = _S448;
    _S669[int(15)] = _S448;
    _S669[int(15)] = _S469;
    _S669[int(14)] = _S471;
    _S669[int(13)] = _S473;
    _S669[int(12)] = _S475;
    _S669[int(11)] = _S477;
    _S669[int(10)] = _S479;
    _S669[int(9)] = _S481;
    _S669[int(8)] = _S489;
    _S669[int(7)] = _S491;
    _S669[int(6)] = _S493;
    _S669[int(5)] = _S495;
    _S669[int(4)] = _S497;
    _S669[int(3)] = _S508;
    _S669[int(2)] = _S510;
    _S669[int(1)] = _S512;
    _S669[int(0)] = _S525;
    float2  _S670 = make_float2 (0.0f, _S554);
    float3  _S671 = make_float3 (_S663, _S660, _S656);
    float4  _S672 = make_float4 (0.0f);
    *&((&_S672)->w) = _S652;
    *&((&_S672)->z) = _S653;
    *&((&_S672)->y) = _S654;
    *&((&_S672)->x) = _S655;
    *v_mean_0 = _S526 + _S619 + _S625 + _S630 + _S665.differential_0;
    *v_quat_0 = _S672;
    *v_scale_0 = _S671;
    *v_hardness_0 = _S670;
    *v_sh_coeffs_0 = _S669;
    *v_ch_coeffs_0 = _S668;
    *v_R_0 = _S667;
    *v_t_0 = _S666;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S673, float _S674)
{
    return (F32_atan2((_S673), (_S674)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S675, DiffPair_float_0 * _S676, float _S677)
{
    _d_atan2_0(_S675, _S676, _S677);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S678 = scale_5.x;
    float _S679 = s_primal_ctx_exp_0(_S678);
    float _S680 = scale_5.y;
    float _S681 = s_primal_ctx_exp_0(_S680);
    float sz_5 = scale_5.z - 0.5f * (_S678 + _S680);
    float _S682 = quat_5.y;
    float x2_5 = _S682 * _S682;
    float y2_5 = quat_5.z * quat_5.z;
    float z2_10 = quat_5.w * quat_5.w;
    float xy_5 = quat_5.y * quat_5.z;
    float xz_5 = quat_5.y * quat_5.w;
    float yz_5 = quat_5.z * quat_5.w;
    float wx_5 = quat_5.x * quat_5.y;
    float wy_5 = quat_5.x * quat_5.z;
    float wz_5 = quat_5.x * quat_5.w;
    Matrix<float, 3, 3>  _S683 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S684 = make_float3 (_S679, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S683, _S684) + mean_5;
    float _S685 = -0.5f + sz_5;
    float3  _S686 = make_float3 (_S679 * _S685, _S681, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S683, _S686) + mean_5;
    float _S687 = -0.5f - sz_5;
    float3  _S688 = make_float3 (_S679 * _S687, - _S681, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S683, _S688) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S689 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S690 = length_0(_S689);
    float _S691 = vert0_c_5.z;
    float _S692 = s_primal_ctx_atan2_0(_S690, _S691);
    bool _S693 = _S692 < 0.00100000004749745f;
    float k_2;
    float _S694;
    float _S695;
    float _S696;
    if(_S693)
    {
        float _S697 = 1.0f - _S692 * _S692 / 3.0f;
        float _S698 = _S691 * _S691;
        k_2 = _S697 / _S691;
        _S694 = _S698;
        _S695 = _S697;
        _S696 = 0.0f;
    }
    else
    {
        float _S699 = _S690 * _S690;
        k_2 = _S692 / _S690;
        _S694 = 0.0f;
        _S695 = 0.0f;
        _S696 = _S699;
    }
    float2  _S700 = make_float2 (k_2);
    float2  _S701 = _S689 * make_float2 (k_2);
    float u_18 = _S701.x;
    float v_18 = _S701.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S702 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S703 = dist_coeffs_5[int(1)] + r2_18 * _S702;
    float _S704 = dist_coeffs_5[int(0)] + r2_18 * _S703;
    float radial_9 = 1.0f + r2_18 * _S704;
    float _S705 = 2.0f * dist_coeffs_5[int(4)];
    float _S706 = _S705 * u_18;
    float _S707 = 2.0f * u_18;
    float _S708 = 2.0f * dist_coeffs_5[int(5)];
    float _S709 = _S708 * u_18;
    float _S710 = 2.0f * v_18;
    float2  _S711 = _S701 * make_float2 (radial_9) + make_float2 (_S706 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S707 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S709 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S710 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S712 = _S711 + make_float2 (dist_coeffs_5[int(8)] * _S711.x + dist_coeffs_5[int(9)] * _S711.y, 0.0f);
    float _S713 = fx_5 * _S712.x + cx_5;
    float _S714 = fy_5 * _S712.y + cy_5;
    float2  uv0_7 = make_float2 (_S713, _S714);
    float2  _S715 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S716 = length_0(_S715);
    float _S717 = vert1_c_5.z;
    float _S718 = s_primal_ctx_atan2_0(_S716, _S717);
    bool _S719 = _S718 < 0.00100000004749745f;
    float _S720;
    float _S721;
    float _S722;
    if(_S719)
    {
        float _S723 = 1.0f - _S718 * _S718 / 3.0f;
        float _S724 = _S717 * _S717;
        k_2 = _S723 / _S717;
        _S720 = _S724;
        _S721 = _S723;
        _S722 = 0.0f;
    }
    else
    {
        float _S725 = _S716 * _S716;
        k_2 = _S718 / _S716;
        _S720 = 0.0f;
        _S721 = 0.0f;
        _S722 = _S725;
    }
    float2  _S726 = make_float2 (k_2);
    float2  _S727 = _S715 * make_float2 (k_2);
    float u_19 = _S727.x;
    float v_19 = _S727.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S728 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S729 = dist_coeffs_5[int(1)] + r2_19 * _S728;
    float _S730 = dist_coeffs_5[int(0)] + r2_19 * _S729;
    float radial_10 = 1.0f + r2_19 * _S730;
    float _S731 = _S705 * u_19;
    float _S732 = 2.0f * u_19;
    float _S733 = _S708 * u_19;
    float _S734 = 2.0f * v_19;
    float2  _S735 = _S727 * make_float2 (radial_10) + make_float2 (_S731 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S732 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S733 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S734 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S736 = _S735 + make_float2 (dist_coeffs_5[int(8)] * _S735.x + dist_coeffs_5[int(9)] * _S735.y, 0.0f);
    float _S737 = fx_5 * _S736.x + cx_5;
    float _S738 = fy_5 * _S736.y + cy_5;
    float2  uv1_7 = make_float2 (_S737, _S738);
    float2  _S739 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S740 = length_0(_S739);
    float _S741 = vert2_c_5.z;
    float _S742 = s_primal_ctx_atan2_0(_S740, _S741);
    bool _S743 = _S742 < 0.00100000004749745f;
    float _S744;
    float _S745;
    float _S746;
    if(_S743)
    {
        float _S747 = 1.0f - _S742 * _S742 / 3.0f;
        float _S748 = _S741 * _S741;
        k_2 = _S747 / _S741;
        _S744 = _S748;
        _S745 = _S747;
        _S746 = 0.0f;
    }
    else
    {
        float _S749 = _S740 * _S740;
        k_2 = _S742 / _S740;
        _S744 = 0.0f;
        _S745 = 0.0f;
        _S746 = _S749;
    }
    float2  _S750 = make_float2 (k_2);
    float2  _S751 = _S739 * make_float2 (k_2);
    float u_20 = _S751.x;
    float v_20 = _S751.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S752 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S753 = dist_coeffs_5[int(1)] + r2_20 * _S752;
    float _S754 = dist_coeffs_5[int(0)] + r2_20 * _S753;
    float radial_11 = 1.0f + r2_20 * _S754;
    float _S755 = _S705 * u_20;
    float _S756 = 2.0f * u_20;
    float _S757 = _S708 * u_20;
    float _S758 = 2.0f * v_20;
    float2  _S759 = _S751 * make_float2 (radial_11) + make_float2 (_S755 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S756 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S757 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S758 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S760 = _S759 + make_float2 (dist_coeffs_5[int(8)] * _S759.x + dist_coeffs_5[int(9)] * _S759.y, 0.0f);
    float _S761 = fx_5 * _S760.x + cx_5;
    float _S762 = fy_5 * _S760.y + cy_5;
    float2  uv2_7 = make_float2 (_S761, _S762);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S763 = e0_5.x;
    float _S764 = e1_5.y;
    float _S765 = e0_5.y;
    float _S766 = e1_5.x;
    float _S767 = _S763 * _S764 - _S765 * _S766;
    float _S768 = 1.0f - hardness_5.y;
    float _S769 = -1.0f / _S768;
    float _S770 = _S768 * _S768;
    float _S771 = (F32_max((_S713), (_S737)));
    float _S772 = (F32_min((_S713), (_S737)));
    float _S773 = (F32_max((_S714), (_S738)));
    float _S774 = (F32_min((_S714), (_S738)));
    float3  _S775 = (vert0_c_5 + vert1_c_5 + vert2_c_5) / make_float3 (3.0f);
    float _S776 = _S775.z;
    float _S777 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S778 = _S777 * (_S776 + length_1(_S775));
    Matrix<float, 3, 3>  _S779 = transpose_1(R_5);
    float3  _S780 = mean_5 - - s_primal_ctx_mul_0(_S779, t_5);
    float _S781 = _S780.x;
    float _S782 = _S780.y;
    float _S783 = _S780.z;
    float _S784 = _S781 * _S781 + _S782 * _S782 + _S783 * _S783;
    float _S785 = s_primal_ctx_sqrt_0(_S784);
    float x_17 = _S781 / _S785;
    float3  _S786 = make_float3 (x_17);
    float _S787 = _S785 * _S785;
    float y_8 = _S782 / _S785;
    float z_5 = _S783 / _S785;
    float3  _S788 = make_float3 (z_5);
    float _S789 = - y_8;
    float3  _S790 = make_float3 (_S789);
    float z2_11 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_17 * x_17 - y_8 * y_8;
    float _S791 = 2.0f * x_17;
    float fS1_5 = _S791 * y_8;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S792 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_17;
    float3  _S793 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_8;
    float3  _S794 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S795 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S796 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    float _S797 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_5 * _S797;
    float3  _S798 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_17;
    float3  _S799 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_8;
    float3  _S800 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S801 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S802 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_17 * fC1_5 - y_8 * fS1_5);
    float3  _S803 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_17 * fS1_5 + y_8 * fC1_5);
    float3  _S804 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S789) * sh_coeffs_5[int(1)] + make_float3 (z_5) * sh_coeffs_5[int(2)] - make_float3 (x_17) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]);
    float3  _S805 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S806 = make_float3 (0.0f);
    float3  _S807 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S808 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S809 = make_float3 (_S808);
    float3  _S810 = ch_coeffs_5[int(1)] * make_float3 (_S808);
    float3  _S811 = _S807 + _S810 + make_float3 (0.5f);
    float3  _S812 = _S807 - _S810 + make_float3 (0.5f);
    float3  _S813 = vert1_c_5 - vert0_c_5;
    float3  _S814 = vert2_c_5 - vert0_c_5;
    float3  _S815 = s_primal_ctx_cross_0(_S813, _S814);
    float3  _S816 = normalize_0(_S815);
    float3  _S817 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S816, mean_c_5)))))) * v_normal_1;
    float3  _S818 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S819;
    (&_S819)->primal_0 = _S816;
    (&_S819)->differential_0 = _S818;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S820;
    (&_S820)->primal_0 = mean_c_5;
    (&_S820)->differential_0 = _S818;
    s_bwd_prop_dot_0(&_S819, &_S820, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S821 = _S820;
    float3  _S822 = _S817 + _S819.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S823;
    (&_S823)->primal_0 = _S815;
    (&_S823)->differential_0 = _S818;
    s_bwd_normalize_impl_0(&_S823, _S822);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S824;
    (&_S824)->primal_0 = _S813;
    (&_S824)->differential_0 = _S818;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S825;
    (&_S825)->primal_0 = _S814;
    (&_S825)->differential_0 = _S818;
    s_bwd_prop_cross_0(&_S824, &_S825, _S823.differential_0);
    float3  _S826 = - _S825.differential_0;
    float3  _S827 = - _S824.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S828;
    (&_S828)->primal_0 = _S812;
    (&_S828)->differential_0 = _S818;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S829;
    (&_S829)->primal_0 = _S806;
    (&_S829)->differential_0 = _S818;
    s_bwd_prop_max_0(&_S828, &_S829, v_rgbs_1[int(2)]);
    float3  _S830 = - _S828.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S831;
    (&_S831)->primal_0 = _S811;
    (&_S831)->differential_0 = _S818;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S832;
    (&_S832)->primal_0 = _S806;
    (&_S832)->differential_0 = _S818;
    s_bwd_prop_max_0(&_S831, &_S832, v_rgbs_1[int(1)]);
    float3  _S833 = _S809 * (_S830 + _S831.differential_0);
    float3  _S834 = _S828.differential_0 + _S831.differential_0;
    float3  _S835 = make_float3 (0.5f) * - _S834;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S836;
    (&_S836)->primal_0 = _S805;
    (&_S836)->differential_0 = _S818;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S837;
    (&_S837)->primal_0 = _S806;
    (&_S837)->differential_0 = _S818;
    s_bwd_prop_max_0(&_S836, &_S837, v_rgbs_1[int(0)]);
    float3  _S838 = _S835 + _S836.differential_0;
    float3  _S839 = _S834 + _S836.differential_0;
    float3  _S840 = _S803 * _S839;
    float3  _S841 = sh_coeffs_5[int(15)] * _S839;
    float3  _S842 = _S801 * _S839;
    float3  _S843 = sh_coeffs_5[int(14)] * _S839;
    float3  _S844 = _S799 * _S839;
    float3  _S845 = sh_coeffs_5[int(13)] * _S839;
    float3  _S846 = _S798 * _S839;
    float3  _S847 = sh_coeffs_5[int(12)] * _S839;
    float3  _S848 = _S800 * _S839;
    float3  _S849 = sh_coeffs_5[int(11)] * _S839;
    float3  _S850 = _S802 * _S839;
    float3  _S851 = sh_coeffs_5[int(10)] * _S839;
    float3  _S852 = _S804 * _S839;
    float3  _S853 = sh_coeffs_5[int(9)] * _S839;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S853.x + _S853.y + _S853.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S841.x + _S841.y + _S841.z);
    float _S854 = _S851.x + _S851.y + _S851.z;
    float _S855 = _S843.x + _S843.y + _S843.z;
    float _S856 = _S849.x + _S849.y + _S849.z;
    float _S857 = _S845.x + _S845.y + _S845.z;
    float _S858 = _S847.x + _S847.y + _S847.z;
    float _S859 = - s_diff_fC2_T_1;
    float3  _S860 = _S795 * _S839;
    float3  _S861 = sh_coeffs_5[int(8)] * _S839;
    float3  _S862 = _S793 * _S839;
    float3  _S863 = sh_coeffs_5[int(7)] * _S839;
    float3  _S864 = _S792 * _S839;
    float3  _S865 = sh_coeffs_5[int(6)] * _S839;
    float3  _S866 = _S794 * _S839;
    float3  _S867 = sh_coeffs_5[int(5)] * _S839;
    float3  _S868 = _S796 * _S839;
    float3  _S869 = sh_coeffs_5[int(4)] * _S839;
    float _S870 = _S867.x + _S867.y + _S867.z;
    float _S871 = _S863.x + _S863.y + _S863.z;
    float _S872 = fTmp1B_5 * _S854 + x_17 * s_diff_fS2_T_1 + y_8 * _S859 + 0.54627424478530884f * (_S869.x + _S869.y + _S869.z);
    float _S873 = fTmp1B_5 * _S855 + y_8 * s_diff_fS2_T_1 + x_17 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S861.x + _S861.y + _S861.z);
    float _S874 = y_8 * - _S873;
    float _S875 = x_17 * _S873;
    float _S876 = z_5 * (1.86588168144226074f * (z_5 * _S858) + -2.28522896766662598f * (y_8 * _S856 + x_17 * _S857) + 0.94617468118667603f * (_S865.x + _S865.y + _S865.z));
    float3  _S877 = make_float3 (0.48860251903533936f) * _S839;
    float3  _S878 = - _S877;
    float3  _S879 = _S786 * _S878;
    float3  _S880 = sh_coeffs_5[int(3)] * _S878;
    float3  _S881 = _S788 * _S877;
    float3  _S882 = sh_coeffs_5[int(2)] * _S877;
    float3  _S883 = _S790 * _S877;
    float3  _S884 = sh_coeffs_5[int(1)] * _S877;
    float _S885 = (_S797 * _S858 + 1.44530570507049561f * (fS1_5 * _S854 + fC1_5 * _S855) + -1.09254848957061768f * (y_8 * _S870 + x_17 * _S871) + _S876 + _S876 + _S882.x + _S882.y + _S882.z) / _S787;
    float _S886 = _S785 * _S885;
    float _S887 = (fTmp0C_5 * _S856 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S859 + fTmp0B_5 * _S870 + _S791 * _S872 + _S874 + _S874 + - (_S884.x + _S884.y + _S884.z)) / _S787;
    float _S888 = _S785 * _S887;
    float _S889 = (fTmp0C_5 * _S857 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S871 + 2.0f * (y_8 * _S872) + _S875 + _S875 + _S880.x + _S880.y + _S880.z) / _S787;
    float _S890 = _S785 * _S889;
    float _S891 = _S783 * - _S885 + _S782 * - _S887 + _S781 * - _S889;
    DiffPair_float_0 _S892;
    (&_S892)->primal_0 = _S784;
    (&_S892)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S892, _S891);
    float _S893 = _S783 * _S892.differential_0;
    float _S894 = _S782 * _S892.differential_0;
    float _S895 = _S781 * _S892.differential_0;
    float3  _S896 = make_float3 (0.282094806432724f) * _S839;
    float3  _S897 = make_float3 (_S890 + _S895 + _S895, _S888 + _S894 + _S894, _S886 + _S893 + _S893);
    float3  _S898 = - - _S897;
    Matrix<float, 3, 3>  _S899 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S900;
    (&_S900)->primal_0 = _S779;
    (&_S900)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S901;
    (&_S901)->primal_0 = t_5;
    (&_S901)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S900, &_S901, _S898);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S902 = _S901;
    Matrix<float, 3, 3>  _S903 = transpose_1(_S900.differential_0);
    DiffPair_float_0 _S904;
    (&_S904)->primal_0 = _S776;
    (&_S904)->differential_0 = 0.0f;
    DiffPair_float_0 _S905;
    (&_S905)->primal_0 = _S778;
    (&_S905)->differential_0 = 0.0f;
    _d_max_0(&_S904, &_S905, v_depth_1);
    float _S906 = _S777 * _S905.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S907;
    (&_S907)->primal_0 = _S775;
    (&_S907)->differential_0 = _S818;
    s_bwd_length_impl_0(&_S907, _S906);
    float3  _S908 = make_float3 (0.3333333432674408f) * (_S907.differential_0 + make_float3 (0.0f, 0.0f, _S904.differential_0 + _S906));
    DiffPair_float_0 _S909;
    (&_S909)->primal_0 = _S774;
    (&_S909)->differential_0 = 0.0f;
    DiffPair_float_0 _S910;
    (&_S910)->primal_0 = _S762;
    (&_S910)->differential_0 = 0.0f;
    _d_min_0(&_S909, &_S910, 0.0f);
    DiffPair_float_0 _S911;
    (&_S911)->primal_0 = _S714;
    (&_S911)->differential_0 = 0.0f;
    DiffPair_float_0 _S912;
    (&_S912)->primal_0 = _S738;
    (&_S912)->differential_0 = 0.0f;
    _d_min_0(&_S911, &_S912, _S909.differential_0);
    DiffPair_float_0 _S913;
    (&_S913)->primal_0 = _S773;
    (&_S913)->differential_0 = 0.0f;
    DiffPair_float_0 _S914;
    (&_S914)->primal_0 = _S762;
    (&_S914)->differential_0 = 0.0f;
    _d_max_0(&_S913, &_S914, 0.0f);
    DiffPair_float_0 _S915;
    (&_S915)->primal_0 = _S714;
    (&_S915)->differential_0 = 0.0f;
    DiffPair_float_0 _S916;
    (&_S916)->primal_0 = _S738;
    (&_S916)->differential_0 = 0.0f;
    _d_max_0(&_S915, &_S916, _S913.differential_0);
    DiffPair_float_0 _S917;
    (&_S917)->primal_0 = _S772;
    (&_S917)->differential_0 = 0.0f;
    DiffPair_float_0 _S918;
    (&_S918)->primal_0 = _S761;
    (&_S918)->differential_0 = 0.0f;
    _d_min_0(&_S917, &_S918, 0.0f);
    DiffPair_float_0 _S919;
    (&_S919)->primal_0 = _S713;
    (&_S919)->differential_0 = 0.0f;
    DiffPair_float_0 _S920;
    (&_S920)->primal_0 = _S737;
    (&_S920)->differential_0 = 0.0f;
    _d_min_0(&_S919, &_S920, _S917.differential_0);
    DiffPair_float_0 _S921;
    (&_S921)->primal_0 = _S771;
    (&_S921)->differential_0 = 0.0f;
    DiffPair_float_0 _S922;
    (&_S922)->primal_0 = _S761;
    (&_S922)->differential_0 = 0.0f;
    _d_max_0(&_S921, &_S922, 0.0f);
    DiffPair_float_0 _S923;
    (&_S923)->primal_0 = _S713;
    (&_S923)->differential_0 = 0.0f;
    DiffPair_float_0 _S924;
    (&_S924)->primal_0 = _S737;
    (&_S924)->differential_0 = 0.0f;
    _d_max_0(&_S923, &_S924, _S921.differential_0);
    DiffPair_float_0 _S925;
    (&_S925)->primal_0 = _S769;
    (&_S925)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S925, 0.0f);
    float _S926 = - (-1.0f * - (_S925.differential_0 / _S770));
    float2  _S927 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S928;
    (&_S928)->primal_0 = e2_1;
    (&_S928)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S928, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S929;
    (&_S929)->primal_0 = e1_5;
    (&_S929)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S929, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S930;
    (&_S930)->primal_0 = e0_5;
    (&_S930)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S930, -0.0f);
    DiffPair_float_0 _S931;
    (&_S931)->primal_0 = _S767;
    (&_S931)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S931, 0.0f);
    float _S932 = - _S931.differential_0;
    float2  _S933 = _S929.differential_0 + make_float2 (_S765 * _S932, _S763 * _S931.differential_0);
    float2  _S934 = - _S933;
    float2  _S935 = _S930.differential_0 + make_float2 (_S764 * _S931.differential_0, _S766 * _S932);
    float2  _S936 = - _S935;
    float2  _S937 = - _S928.differential_0 + _S933;
    float _S938 = fx_5 * (_S918.differential_0 + _S922.differential_0 + _S937.x);
    float2  _S939 = make_float2 (_S938, fy_5 * (_S910.differential_0 + _S914.differential_0 + _S937.y)) + make_float2 (dist_coeffs_5[int(8)] * _S938, dist_coeffs_5[int(9)] * _S938);
    float2  _S940 = _S751 * _S939;
    float _S941 = dist_coeffs_5[int(4)] * _S939.y;
    float _S942 = dist_coeffs_5[int(5)] * _S939.x;
    float _S943 = _S940.x + _S940.y;
    float _S944 = r2_20 * _S943;
    float _S945 = r2_20 * _S944;
    float _S946 = dist_coeffs_5[int(7)] * _S939.y + _S941 + dist_coeffs_5[int(6)] * _S939.x + _S942 + _S754 * _S943 + _S753 * _S944 + _S752 * _S945 + dist_coeffs_5[int(3)] * (r2_20 * _S945);
    float _S947 = v_20 * _S946;
    float _S948 = u_20 * _S946;
    float2  _S949 = make_float2 (radial_11) * _S939 + make_float2 (_S708 * (v_20 * _S939.y) + _S756 * _S942 + 2.0f * (u_20 * _S942) + _S705 * (v_20 * _S939.x) + _S948 + _S948, _S758 * _S941 + 2.0f * (v_20 * _S941) + _S757 * _S939.y + _S755 * _S939.x + _S947 + _S947);
    float3  _S950 = _S824.differential_0 + _S908;
    float3  _S951 = _S826 + _S827 + _S908;
    float3  _S952 = _S825.differential_0 + _S908;
    FixedArray<float3 , 2>  _S953;
    _S953[int(0)] = _S818;
    _S953[int(1)] = _S818;
    _S953[int(1)] = _S833;
    _S953[int(0)] = _S838;
    float3  _S954 = _S953[int(0)];
    float3  _S955 = _S953[int(1)];
    FixedArray<float3 , 16>  _S956;
    _S956[int(0)] = _S818;
    _S956[int(1)] = _S818;
    _S956[int(2)] = _S818;
    _S956[int(3)] = _S818;
    _S956[int(4)] = _S818;
    _S956[int(5)] = _S818;
    _S956[int(6)] = _S818;
    _S956[int(7)] = _S818;
    _S956[int(8)] = _S818;
    _S956[int(9)] = _S818;
    _S956[int(10)] = _S818;
    _S956[int(11)] = _S818;
    _S956[int(12)] = _S818;
    _S956[int(13)] = _S818;
    _S956[int(14)] = _S818;
    _S956[int(15)] = _S818;
    _S956[int(7)] = _S862;
    _S956[int(0)] = _S896;
    _S956[int(1)] = _S883;
    _S956[int(2)] = _S881;
    _S956[int(3)] = _S879;
    _S956[int(4)] = _S868;
    _S956[int(5)] = _S866;
    _S956[int(6)] = _S864;
    _S956[int(15)] = _S840;
    _S956[int(8)] = _S860;
    _S956[int(9)] = _S852;
    _S956[int(10)] = _S850;
    _S956[int(11)] = _S848;
    _S956[int(12)] = _S846;
    _S956[int(13)] = _S844;
    _S956[int(14)] = _S842;
    float3  _S957 = _S956[int(0)];
    float3  _S958 = _S956[int(1)];
    float3  _S959 = _S956[int(2)];
    float3  _S960 = _S956[int(3)];
    float3  _S961 = _S956[int(4)];
    float3  _S962 = _S956[int(5)];
    float3  _S963 = _S956[int(6)];
    float3  _S964 = _S956[int(7)];
    float3  _S965 = _S956[int(8)];
    float3  _S966 = _S956[int(9)];
    float3  _S967 = _S956[int(10)];
    float3  _S968 = _S956[int(11)];
    float3  _S969 = _S956[int(12)];
    float3  _S970 = _S956[int(13)];
    float3  _S971 = _S956[int(14)];
    float3  _S972 = _S956[int(15)];
    float _S973 = _S919.differential_0 + _S923.differential_0;
    float2  _S974 = _S928.differential_0 + _S936;
    float _S975 = _S912.differential_0 + _S916.differential_0;
    float _S976 = _S911.differential_0 + _S915.differential_0;
    float2  _S977 = _S934 + _S935;
    float _S978 = _S920.differential_0 + _S924.differential_0;
    float2  _S979 = make_float2 (0.0f, _S926);
    float2  _S980 = _S739 * _S949;
    float2  _S981 = _S750 * _S949;
    float _S982 = _S980.x + _S980.y;
    if(_S743)
    {
        float _S983 = _S982 / _S744;
        float _S984 = _S745 * - _S983;
        float _S985 = _S742 * (0.3333333432674408f * - (_S741 * _S983));
        k_2 = _S985 + _S985;
        _S744 = _S984;
        _S745 = 0.0f;
    }
    else
    {
        float _S986 = _S982 / _S746;
        float _S987 = _S742 * - _S986;
        k_2 = _S740 * _S986;
        _S744 = 0.0f;
        _S745 = _S987;
    }
    DiffPair_float_0 _S988;
    (&_S988)->primal_0 = _S740;
    (&_S988)->differential_0 = 0.0f;
    DiffPair_float_0 _S989;
    (&_S989)->primal_0 = _S741;
    (&_S989)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S988, &_S989, k_2);
    float _S990 = _S989.differential_0 + _S744;
    float _S991 = _S988.differential_0 + _S745;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S992;
    (&_S992)->primal_0 = _S739;
    (&_S992)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S992, _S991);
    float2  _S993 = _S992.differential_0 + _S981;
    float _S994 = fx_5 * (_S977.x + _S978);
    float2  _S995 = make_float2 (_S994, fy_5 * (_S977.y + _S975)) + make_float2 (dist_coeffs_5[int(8)] * _S994, dist_coeffs_5[int(9)] * _S994);
    float2  _S996 = _S727 * _S995;
    float _S997 = dist_coeffs_5[int(4)] * _S995.y;
    float _S998 = dist_coeffs_5[int(5)] * _S995.x;
    float _S999 = _S996.x + _S996.y;
    float _S1000 = r2_19 * _S999;
    float _S1001 = r2_19 * _S1000;
    float _S1002 = dist_coeffs_5[int(7)] * _S995.y + _S997 + dist_coeffs_5[int(6)] * _S995.x + _S998 + _S730 * _S999 + _S729 * _S1000 + _S728 * _S1001 + dist_coeffs_5[int(3)] * (r2_19 * _S1001);
    float _S1003 = v_19 * _S1002;
    float _S1004 = u_19 * _S1002;
    float2  _S1005 = make_float2 (radial_10) * _S995 + make_float2 (_S708 * (v_19 * _S995.y) + _S732 * _S998 + 2.0f * (u_19 * _S998) + _S705 * (v_19 * _S995.x) + _S1004 + _S1004, _S734 * _S997 + 2.0f * (v_19 * _S997) + _S733 * _S995.y + _S731 * _S995.x + _S1003 + _S1003);
    float3  _S1006 = _S952 + make_float3 (_S993.x, _S993.y, _S990);
    float2  _S1007 = _S715 * _S1005;
    float2  _S1008 = _S726 * _S1005;
    float _S1009 = _S1007.x + _S1007.y;
    if(_S719)
    {
        float _S1010 = _S1009 / _S720;
        float _S1011 = _S721 * - _S1010;
        float _S1012 = _S718 * (0.3333333432674408f * - (_S717 * _S1010));
        k_2 = _S1012 + _S1012;
        _S720 = _S1011;
        _S721 = 0.0f;
    }
    else
    {
        float _S1013 = _S1009 / _S722;
        float _S1014 = _S718 * - _S1013;
        k_2 = _S716 * _S1013;
        _S720 = 0.0f;
        _S721 = _S1014;
    }
    DiffPair_float_0 _S1015;
    (&_S1015)->primal_0 = _S716;
    (&_S1015)->differential_0 = 0.0f;
    DiffPair_float_0 _S1016;
    (&_S1016)->primal_0 = _S717;
    (&_S1016)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1015, &_S1016, k_2);
    float _S1017 = _S1016.differential_0 + _S720;
    float _S1018 = _S1015.differential_0 + _S721;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1019;
    (&_S1019)->primal_0 = _S715;
    (&_S1019)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S1019, _S1018);
    float2  _S1020 = _S1019.differential_0 + _S1008;
    float _S1021 = fx_5 * (_S974.x + _S973);
    float2  _S1022 = make_float2 (_S1021, fy_5 * (_S974.y + _S976)) + make_float2 (dist_coeffs_5[int(8)] * _S1021, dist_coeffs_5[int(9)] * _S1021);
    float2  _S1023 = _S701 * _S1022;
    float _S1024 = dist_coeffs_5[int(4)] * _S1022.y;
    float _S1025 = dist_coeffs_5[int(5)] * _S1022.x;
    float _S1026 = _S1023.x + _S1023.y;
    float _S1027 = r2_18 * _S1026;
    float _S1028 = r2_18 * _S1027;
    float _S1029 = dist_coeffs_5[int(7)] * _S1022.y + _S1024 + dist_coeffs_5[int(6)] * _S1022.x + _S1025 + _S704 * _S1026 + _S703 * _S1027 + _S702 * _S1028 + dist_coeffs_5[int(3)] * (r2_18 * _S1028);
    float _S1030 = v_18 * _S1029;
    float _S1031 = u_18 * _S1029;
    float2  _S1032 = make_float2 (radial_9) * _S1022 + make_float2 (_S708 * (v_18 * _S1022.y) + _S707 * _S1025 + 2.0f * (u_18 * _S1025) + _S705 * (v_18 * _S1022.x) + _S1031 + _S1031, _S710 * _S1024 + 2.0f * (v_18 * _S1024) + _S709 * _S1022.y + _S706 * _S1022.x + _S1030 + _S1030);
    float3  _S1033 = _S950 + make_float3 (_S1020.x, _S1020.y, _S1017);
    float2  _S1034 = _S689 * _S1032;
    float2  _S1035 = _S700 * _S1032;
    float _S1036 = _S1034.x + _S1034.y;
    if(_S693)
    {
        float _S1037 = _S1036 / _S694;
        float _S1038 = _S695 * - _S1037;
        float _S1039 = _S692 * (0.3333333432674408f * - (_S691 * _S1037));
        k_2 = _S1039 + _S1039;
        _S694 = _S1038;
        _S695 = 0.0f;
    }
    else
    {
        float _S1040 = _S1036 / _S696;
        float _S1041 = _S692 * - _S1040;
        k_2 = _S690 * _S1040;
        _S694 = 0.0f;
        _S695 = _S1041;
    }
    DiffPair_float_0 _S1042;
    (&_S1042)->primal_0 = _S690;
    (&_S1042)->differential_0 = 0.0f;
    DiffPair_float_0 _S1043;
    (&_S1043)->primal_0 = _S691;
    (&_S1043)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1042, &_S1043, k_2);
    float _S1044 = _S1043.differential_0 + _S694;
    float _S1045 = _S1042.differential_0 + _S695;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1046;
    (&_S1046)->primal_0 = _S689;
    (&_S1046)->differential_0 = _S927;
    s_bwd_length_impl_1(&_S1046, _S1045);
    float2  _S1047 = _S1046.differential_0 + _S1035;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1048;
    (&_S1048)->primal_0 = vert2_c_5;
    (&_S1048)->differential_0 = _S818;
    s_bwd_length_impl_0(&_S1048, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1049;
    (&_S1049)->primal_0 = vert1_c_5;
    (&_S1049)->differential_0 = _S818;
    s_bwd_length_impl_0(&_S1049, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1050;
    (&_S1050)->primal_0 = vert0_c_5;
    (&_S1050)->differential_0 = _S818;
    s_bwd_length_impl_0(&_S1050, 0.0f);
    float3  _S1051 = _S1048.differential_0 + _S1006;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1052;
    (&_S1052)->primal_0 = R_5;
    (&_S1052)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1053;
    (&_S1053)->primal_0 = vert2_5;
    (&_S1053)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1052, &_S1053, _S1051);
    float3  _S1054 = _S1049.differential_0 + _S1033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1055;
    (&_S1055)->primal_0 = R_5;
    (&_S1055)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1056;
    (&_S1056)->primal_0 = vert1_5;
    (&_S1056)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1055, &_S1056, _S1054);
    float3  _S1057 = _S1050.differential_0 + _S951 + make_float3 (_S1047.x, _S1047.y, _S1044);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1058;
    (&_S1058)->primal_0 = R_5;
    (&_S1058)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1059;
    (&_S1059)->primal_0 = vert0_5;
    (&_S1059)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1058, &_S1059, _S1057);
    float3  _S1060 = _S1053.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1061;
    (&_S1061)->primal_0 = _S683;
    (&_S1061)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1062;
    (&_S1062)->primal_0 = _S688;
    (&_S1062)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1061, &_S1062, _S1060);
    float _S1063 = - _S1062.differential_0.y;
    float _S1064 = _S687 * _S1062.differential_0.x;
    float _S1065 = - (_S679 * _S1062.differential_0.x);
    float3  _S1066 = _S1056.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1067;
    (&_S1067)->primal_0 = _S683;
    (&_S1067)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1068;
    (&_S1068)->primal_0 = _S686;
    (&_S1068)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1067, &_S1068, _S1066);
    float _S1069 = _S679 * _S1068.differential_0.x;
    float _S1070 = _S685 * _S1068.differential_0.x;
    float3  _S1071 = _S1059.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1072;
    (&_S1072)->primal_0 = _S683;
    (&_S1072)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1073;
    (&_S1073)->primal_0 = _S684;
    (&_S1073)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1072, &_S1073, _S1071);
    Matrix<float, 3, 3>  _S1074 = transpose_1(_S1061.differential_0 + _S1067.differential_0 + _S1072.differential_0);
    float _S1075 = 2.0f * - _S1074.rows[int(2)].z;
    float _S1076 = 2.0f * _S1074.rows[int(2)].y;
    float _S1077 = 2.0f * _S1074.rows[int(2)].x;
    float _S1078 = 2.0f * _S1074.rows[int(1)].z;
    float _S1079 = 2.0f * - _S1074.rows[int(1)].y;
    float _S1080 = 2.0f * _S1074.rows[int(1)].x;
    float _S1081 = 2.0f * _S1074.rows[int(0)].z;
    float _S1082 = 2.0f * _S1074.rows[int(0)].y;
    float _S1083 = 2.0f * - _S1074.rows[int(0)].x;
    float _S1084 = - _S1080 + _S1082;
    float _S1085 = _S1077 + - _S1081;
    float _S1086 = - _S1076 + _S1078;
    float _S1087 = _S1076 + _S1078;
    float _S1088 = _S1077 + _S1081;
    float _S1089 = _S1080 + _S1082;
    float _S1090 = quat_5.w * (_S1079 + _S1083);
    float _S1091 = quat_5.z * (_S1075 + _S1083);
    float _S1092 = quat_5.y * (_S1075 + _S1079);
    float _S1093 = quat_5.x * _S1084 + quat_5.z * _S1087 + quat_5.y * _S1088 + _S1090 + _S1090;
    float _S1094 = quat_5.x * _S1085 + quat_5.w * _S1087 + quat_5.y * _S1089 + _S1091 + _S1091;
    float _S1095 = quat_5.x * _S1086 + quat_5.w * _S1088 + quat_5.z * _S1089 + _S1092 + _S1092;
    float _S1096 = quat_5.w * _S1084 + quat_5.z * _S1085 + quat_5.y * _S1086;
    float _S1097 = _S1065 + _S1069;
    float _S1098 = 0.5f * - _S1097;
    float _S1099 = _S1063 + _S1068.differential_0.y;
    DiffPair_float_0 _S1100;
    (&_S1100)->primal_0 = _S680;
    (&_S1100)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1100, _S1099);
    float _S1101 = _S1098 + _S1100.differential_0;
    float _S1102 = _S1064 + _S1070 + _S1073.differential_0.x;
    DiffPair_float_0 _S1103;
    (&_S1103)->primal_0 = _S678;
    (&_S1103)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1103, _S1102);
    float _S1104 = _S1098 + _S1103.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1105;
    (&_S1105)->primal_0 = mean_c_5;
    (&_S1105)->differential_0 = _S818;
    s_bwd_length_impl_0(&_S1105, 0.0f);
    float3  _S1106 = _S1105.differential_0 + _S821.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1107;
    (&_S1107)->primal_0 = R_5;
    (&_S1107)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1108;
    (&_S1108)->primal_0 = mean_5;
    (&_S1108)->differential_0 = _S818;
    s_bwd_prop_mul_0(&_S1107, &_S1108, _S1106);
    float3  _S1109 = _S1051 + _S1054 + _S1057 + _S1106 + _S902.differential_0;
    Matrix<float, 3, 3>  _S1110 = _S1052.differential_0 + _S1055.differential_0 + _S1058.differential_0 + _S1107.differential_0 + _S903;
    float3  _S1111 = make_float3 (_S1104, _S1101, _S1097);
    float4  _S1112 = make_float4 (0.0f);
    *&((&_S1112)->w) = _S1093;
    *&((&_S1112)->z) = _S1094;
    *&((&_S1112)->y) = _S1095;
    *&((&_S1112)->x) = _S1096;
    float4  _S1113 = _S1112;
    float3  _S1114 = _S1060 + _S1066 + _S1071 + _S1108.differential_0 + _S897;
    *v_mean_1 = _S1114;
    *v_quat_1 = _S1113;
    *v_scale_1 = _S1111;
    *v_hardness_1 = _S979;
    (*v_sh_coeffs_1)[int(0)] = _S957;
    (*v_sh_coeffs_1)[int(1)] = _S958;
    (*v_sh_coeffs_1)[int(2)] = _S959;
    (*v_sh_coeffs_1)[int(3)] = _S960;
    (*v_sh_coeffs_1)[int(4)] = _S961;
    (*v_sh_coeffs_1)[int(5)] = _S962;
    (*v_sh_coeffs_1)[int(6)] = _S963;
    (*v_sh_coeffs_1)[int(7)] = _S964;
    (*v_sh_coeffs_1)[int(8)] = _S965;
    (*v_sh_coeffs_1)[int(9)] = _S966;
    (*v_sh_coeffs_1)[int(10)] = _S967;
    (*v_sh_coeffs_1)[int(11)] = _S968;
    (*v_sh_coeffs_1)[int(12)] = _S969;
    (*v_sh_coeffs_1)[int(13)] = _S970;
    (*v_sh_coeffs_1)[int(14)] = _S971;
    (*v_sh_coeffs_1)[int(15)] = _S972;
    (*v_ch_coeffs_1)[int(0)] = _S954;
    (*v_ch_coeffs_1)[int(1)] = _S955;
    *v_R_1 = _S1110;
    *v_t_1 = _S1109;
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
    bool _S1115;
    if((*u_21) >= 0.0f)
    {
        _S1115 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1115 = false;
    }
    if(_S1115)
    {
        _S1115 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1115 = false;
    }
    if(_S1115)
    {
        _S1115 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1115 = false;
    }
    return _S1115;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1116 = *dpx_12;
    bool _S1117;
    if(((*dpx_12).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1117 = ((*dpx_12).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1117 = false;
    }
    float _S1118;
    if(_S1117)
    {
        _S1118 = dOut_11;
    }
    else
    {
        _S1118 = 0.0f;
    }
    dpx_12->primal_0 = _S1116.primal_0;
    dpx_12->differential_0 = _S1118;
    DiffPair_float_0 _S1119 = *dpMin_0;
    if((_S1116.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1118 = dOut_11;
    }
    else
    {
        _S1118 = 0.0f;
    }
    dpMin_0->primal_0 = _S1119.primal_0;
    dpMin_0->differential_0 = _S1118;
    DiffPair_float_0 _S1120 = *dpMax_0;
    if(((*dpx_12).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1118 = dOut_11;
    }
    else
    {
        _S1118 = 0.0f;
    }
    dpMax_0->primal_0 = _S1120.primal_0;
    dpMax_0->differential_0 = _S1118;
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
        DiffPair_float_0 _S1121 = *dpx_13;
        float _S1122 = val_0 * (*dpy_5).primal_0 / (*dpx_13).primal_0 * dOut_12;
        dpx_13->primal_0 = (*dpx_13).primal_0;
        dpx_13->differential_0 = _S1122;
        float _S1123 = val_0 * (F32_log((_S1121.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1123;
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
    bool _S1124;
    if(u_22 >= 0.0f)
    {
        _S1124 = v_22 >= 0.0f;
    }
    else
    {
        _S1124 = false;
    }
    if(_S1124)
    {
        _S1124 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1124 = false;
    }
    if(_S1124)
    {
        _S1124 = t_7 >= 0.0f;
    }
    else
    {
        _S1124 = false;
    }
    if(!_S1124)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1125;
    if(opac_0 < 0.0f)
    {
        _S1125 = 0.0f;
    }
    else
    {
        _S1125 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1125;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1126, float _S1127, float _S1128)
{
    return clamp_0(_S1126, _S1127, _S1128);
}

inline __device__ float s_primal_ctx_pow_0(float _S1129, float _S1130)
{
    return (F32_pow((_S1129), (_S1130)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1131, DiffPair_float_0 * _S1132, float _S1133)
{
    _d_pow_0(_S1131, _S1132, _S1133);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1134, DiffPair_float_0 * _S1135, DiffPair_float_0 * _S1136, float _S1137)
{
    _d_clamp_0(_S1134, _S1135, _S1136, _S1137);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_3)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1138 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1139 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1140 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1141 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1142 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1140);
    float d_2 = 1.0f / _S1142;
    float _S1143 = _S1142 * _S1142;
    float3  _S1144 = - _S1141;
    float _S1145 = s_primal_ctx_dot_0(_S1144, v2v0_2);
    float u_23 = d_2 * _S1145;
    float _S1146 = s_primal_ctx_dot_0(_S1141, v1v0_2);
    float v_23 = d_2 * _S1146;
    float3  _S1147 = - _S1140;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1147, rov0_2);
    bool _S1148;
    if(u_23 >= 0.0f)
    {
        _S1148 = v_23 >= 0.0f;
    }
    else
    {
        _S1148 = false;
    }
    if(_S1148)
    {
        _S1148 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1148 = false;
    }
    if(_S1148)
    {
        _S1148 = t_8 >= 0.0f;
    }
    else
    {
        _S1148 = false;
    }
    bool _S1149 = !!_S1148;
    float _S1150;
    float _S1151;
    float _S1152;
    float _S1153;
    float _S1154;
    float _S1155;
    float _S1156;
    float _S1157;
    float _S1158;
    float _S1159;
    float _S1160;
    if(_S1149)
    {
        float _S1161 = (F32_min((u_23), (v_23)));
        float _S1162 = s_primal_ctx_sqrt_0(0.5f);
        float _S1163 = _S1162 * (1.0f - u_23 - v_23);
        float _S1164 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1161), (_S1163))) * _S1164;
        float _S1165 = _S1138.primal_0.y;
        float _S1166 = 1.0f - opac_1;
        float _S1167 = 1.0f - s_primal_ctx_clamp_0(_S1165, 0.0f, 0.99989998340606689f);
        float _S1168 = 1.0f / _S1167;
        float _S1169 = _S1167 * _S1167;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1166, _S1168);
        float o_1 = _S1138.primal_0.x;
        bool _S1170 = opac_1 < 0.0f;
        if(_S1170)
        {
            _S1150 = 0.0f;
        }
        else
        {
            _S1150 = o_1 * w_1;
        }
        _S1148 = _S1170;
        _S1151 = o_1;
        _S1152 = w_1;
        _S1153 = _S1166;
        _S1154 = _S1168;
        _S1155 = _S1169;
        _S1156 = _S1165;
        _S1157 = _S1164;
        _S1158 = _S1161;
        _S1159 = _S1163;
        _S1160 = _S1162;
    }
    else
    {
        _S1148 = false;
        _S1150 = 0.0f;
        _S1151 = 0.0f;
        _S1152 = 0.0f;
        _S1153 = 0.0f;
        _S1154 = 0.0f;
        _S1155 = 0.0f;
        _S1156 = 0.0f;
        _S1157 = 0.0f;
        _S1158 = 0.0f;
        _S1159 = 0.0f;
        _S1160 = 0.0f;
    }
    float2  _S1171 = make_float2 (0.0f);
    float2  _S1172;
    if(_S1149)
    {
        if(_S1148)
        {
            _S1150 = 0.0f;
            _S1151 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1173;
            (&_S1173)->primal_0 = _S1150;
            (&_S1173)->differential_0 = 0.0f;
            DiffPair_float_0 _S1174;
            (&_S1174)->primal_0 = 0.99500000476837158f;
            (&_S1174)->differential_0 = 0.0f;
            _d_min_0(&_S1173, &_S1174, _s_dOut_3);
            float _S1175 = _S1151 * _S1173.differential_0;
            _S1150 = _S1152 * _S1173.differential_0;
            _S1151 = _S1175;
        }
        float _S1176 = - _S1151;
        DiffPair_float_0 _S1177;
        (&_S1177)->primal_0 = _S1153;
        (&_S1177)->differential_0 = 0.0f;
        DiffPair_float_0 _S1178;
        (&_S1178)->primal_0 = _S1154;
        (&_S1178)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1177, &_S1178, _S1176);
        float _S1179 = - - (_S1178.differential_0 / _S1155);
        float s_diff_opac_T_0 = - _S1177.differential_0;
        DiffPair_float_0 _S1180;
        (&_S1180)->primal_0 = _S1156;
        (&_S1180)->differential_0 = 0.0f;
        DiffPair_float_0 _S1181;
        (&_S1181)->primal_0 = 0.0f;
        (&_S1181)->differential_0 = 0.0f;
        DiffPair_float_0 _S1182;
        (&_S1182)->primal_0 = 0.99989998340606689f;
        (&_S1182)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1180, &_S1181, &_S1182, _S1179);
        float _S1183 = _S1157 * s_diff_opac_T_0;
        DiffPair_float_0 _S1184;
        (&_S1184)->primal_0 = _S1158;
        (&_S1184)->differential_0 = 0.0f;
        DiffPair_float_0 _S1185;
        (&_S1185)->primal_0 = _S1159;
        (&_S1185)->differential_0 = 0.0f;
        _d_min_0(&_S1184, &_S1185, _S1183);
        float _S1186 = - (_S1160 * _S1185.differential_0);
        DiffPair_float_0 _S1187;
        (&_S1187)->primal_0 = u_23;
        (&_S1187)->differential_0 = 0.0f;
        DiffPair_float_0 _S1188;
        (&_S1188)->primal_0 = v_23;
        (&_S1188)->differential_0 = 0.0f;
        _d_min_0(&_S1187, &_S1188, _S1184.differential_0);
        float2  _S1189 = make_float2 (_S1150, _S1180.differential_0);
        float _S1190 = _S1186 + _S1188.differential_0;
        _S1150 = _S1186 + _S1187.differential_0;
        _S1151 = _S1190;
        _S1172 = _S1189;
    }
    else
    {
        _S1150 = 0.0f;
        _S1151 = 0.0f;
        _S1172 = _S1171;
    }
    float3  _S1191 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1192;
    (&_S1192)->primal_0 = _S1147;
    (&_S1192)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1193;
    (&_S1193)->primal_0 = rov0_2;
    (&_S1193)->differential_0 = _S1191;
    s_bwd_prop_dot_0(&_S1192, &_S1193, 0.0f);
    float3  _S1194 = - _S1192.differential_0;
    float _S1195 = d_2 * _S1151;
    float _S1196 = _S1146 * _S1151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1197;
    (&_S1197)->primal_0 = _S1141;
    (&_S1197)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1198;
    (&_S1198)->primal_0 = v1v0_2;
    (&_S1198)->differential_0 = _S1191;
    s_bwd_prop_dot_0(&_S1197, &_S1198, _S1195);
    float _S1199 = d_2 * _S1150;
    float _S1200 = _S1145 * _S1150;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1201;
    (&_S1201)->primal_0 = _S1144;
    (&_S1201)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1202;
    (&_S1202)->primal_0 = v2v0_2;
    (&_S1202)->differential_0 = _S1191;
    s_bwd_prop_dot_0(&_S1201, &_S1202, _S1199);
    float3  _S1203 = - _S1201.differential_0;
    float _S1204 = - ((_S1196 + _S1200) / _S1143);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1205;
    (&_S1205)->primal_0 = _S1139.primal_0;
    (&_S1205)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1206;
    (&_S1206)->primal_0 = _S1140;
    (&_S1206)->differential_0 = _S1191;
    s_bwd_prop_dot_0(&_S1205, &_S1206, _S1204);
    float3  _S1207 = _S1197.differential_0 + _S1203;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1208;
    (&_S1208)->primal_0 = rov0_2;
    (&_S1208)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1209;
    (&_S1209)->primal_0 = _S1139.primal_0;
    (&_S1209)->differential_0 = _S1191;
    s_bwd_prop_cross_0(&_S1208, &_S1209, _S1207);
    float3  _S1210 = _S1194 + _S1206.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1211;
    (&_S1211)->primal_0 = v1v0_2;
    (&_S1211)->differential_0 = _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1212;
    (&_S1212)->primal_0 = v2v0_2;
    (&_S1212)->differential_0 = _S1191;
    s_bwd_prop_cross_0(&_S1211, &_S1212, _S1210);
    float3  _S1213 = _S1193.differential_0 + _S1208.differential_0;
    float3  _S1214 = _S1202.differential_0 + _S1212.differential_0;
    float3  _S1215 = _S1198.differential_0 + _S1211.differential_0;
    float3  _S1216 = - _S1213 + - _S1214 + - _S1215;
    float3  _S1217 = _S1205.differential_0 + _S1209.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1217;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1213;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1172;
    FixedArray<float3 , 3>  _S1218;
    _S1218[int(0)] = _S1191;
    _S1218[int(1)] = _S1191;
    _S1218[int(2)] = _S1191;
    _S1218[int(2)] = _S1214;
    _S1218[int(0)] = _S1216;
    _S1218[int(1)] = _S1215;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1218;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1219, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1220, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1221, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1222, float _S1223)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1219, _S1220, _S1221, _S1222, _S1223);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1224 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1225 = { _S1224, _S1224, _S1224 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1225;
    float2  _S1226 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1226;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1224;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1224;
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
    float3  _S1227 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1228 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1229 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1227);
    float d_4 = 1.0f / _S1229;
    float _S1230 = _S1229 * _S1229;
    float3  _S1231 = - _S1228;
    float _S1232 = s_primal_ctx_dot_0(_S1231, v2v0_4);
    float u_25 = d_4 * _S1232;
    float _S1233 = s_primal_ctx_dot_0(_S1228, v1v0_4);
    float v_25 = d_4 * _S1233;
    float3  _S1234 = - _S1227;
    float3  _S1235 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1236 = make_float3 (v_25) * dpcolor_0;
    float3  _S1237 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1238 = make_float3 (u_25) * dpcolor_0;
    float3  _S1239 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1240 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1241 = - (_S1239.x + _S1239.y + _S1239.z);
    float _S1242 = d_4 * dpdepth_0;
    float _S1243 = s_primal_ctx_dot_0(_S1234, rov0_4) * dpdepth_0;
    float3  _S1244 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1245;
    (&_S1245)->primal_0 = _S1234;
    (&_S1245)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1246;
    (&_S1246)->primal_0 = rov0_4;
    (&_S1246)->differential_0 = _S1244;
    s_bwd_prop_dot_0(&_S1245, &_S1246, _S1242);
    float3  _S1247 = - _S1245.differential_0;
    float _S1248 = _S1241 + _S1235.x + _S1235.y + _S1235.z;
    float _S1249 = d_4 * _S1248;
    float _S1250 = _S1233 * _S1248;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1251;
    (&_S1251)->primal_0 = _S1228;
    (&_S1251)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1252;
    (&_S1252)->primal_0 = v1v0_4;
    (&_S1252)->differential_0 = _S1244;
    s_bwd_prop_dot_0(&_S1251, &_S1252, _S1249);
    float _S1253 = _S1241 + _S1237.x + _S1237.y + _S1237.z;
    float _S1254 = d_4 * _S1253;
    float _S1255 = _S1232 * _S1253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1256;
    (&_S1256)->primal_0 = _S1231;
    (&_S1256)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1257;
    (&_S1257)->primal_0 = v2v0_4;
    (&_S1257)->differential_0 = _S1244;
    s_bwd_prop_dot_0(&_S1256, &_S1257, _S1254);
    float3  _S1258 = - _S1256.differential_0;
    float _S1259 = - ((_S1243 + _S1250 + _S1255) / _S1230);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1260;
    (&_S1260)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1260)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1261;
    (&_S1261)->primal_0 = _S1227;
    (&_S1261)->differential_0 = _S1244;
    s_bwd_prop_dot_0(&_S1260, &_S1261, _S1259);
    float3  _S1262 = _S1251.differential_0 + _S1258;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1263;
    (&_S1263)->primal_0 = rov0_4;
    (&_S1263)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1264;
    (&_S1264)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1264)->differential_0 = _S1244;
    s_bwd_prop_cross_0(&_S1263, &_S1264, _S1262);
    float3  _S1265 = _S1247 + _S1261.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1266;
    (&_S1266)->primal_0 = v1v0_4;
    (&_S1266)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1267;
    (&_S1267)->primal_0 = v2v0_4;
    (&_S1267)->differential_0 = _S1244;
    s_bwd_prop_cross_0(&_S1266, &_S1267, _S1265);
    float3  _S1268 = _S1246.differential_0 + _S1263.differential_0;
    float3  _S1269 = _S1257.differential_0 + _S1267.differential_0;
    float3  _S1270 = _S1252.differential_0 + _S1266.differential_0;
    float3  _S1271 = - _S1268 + - _S1269 + - _S1270;
    float3  _S1272 = _S1260.differential_0 + _S1264.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1272;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1268;
    FixedArray<float3 , 3>  _S1273;
    _S1273[int(0)] = _S1244;
    _S1273[int(1)] = _S1244;
    _S1273[int(2)] = _S1244;
    _S1273[int(2)] = _S1236;
    _S1273[int(1)] = _S1238;
    _S1273[int(0)] = _S1240;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1273;
    FixedArray<float3 , 3>  _S1274;
    _S1274[int(0)] = _S1244;
    _S1274[int(1)] = _S1244;
    _S1274[int(2)] = _S1244;
    _S1274[int(2)] = _S1269;
    _S1274[int(0)] = _S1271;
    _S1274[int(1)] = _S1270;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1274;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1275, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1276, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1277, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1278, float3  _S1279, float _S1280)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1275, _S1276, _S1277, _S1278, _S1279, _S1280);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1281 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1282 = { _S1281, _S1281, _S1281 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1282;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1282;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1281;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

