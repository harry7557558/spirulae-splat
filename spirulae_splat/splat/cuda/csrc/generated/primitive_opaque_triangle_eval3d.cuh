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

inline __device__ float dot_2(float4  x_2, float4  y_2)
{
    int i_2 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_2 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_5 = result_4 + _slang_vector_get_element(x_2, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_4 = result_5;
    }
    return result_4;
}

inline __device__ float length_0(float2  x_3)
{
    return (F32_sqrt((dot_1(x_3, x_3))));
}

inline __device__ float length_1(float3  x_4)
{
    return (F32_sqrt((dot_0(x_4, x_4))));
}

inline __device__ float length_2(float4  x_5)
{
    return (F32_sqrt((dot_2(x_5, x_5))));
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

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_6)
{
    Matrix<float, 2, 2>  result_6;
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
            *_slang_vector_get_element_ptr(((&result_6)->rows + (r_0)), c_0) = _slang_vector_get_element(x_6.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_6;
}

inline __device__ Matrix<float, 3, 3>  transpose_1(Matrix<float, 3, 3>  x_7)
{
    Matrix<float, 3, 3>  result_7;
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
            *_slang_vector_get_element_ptr(((&result_7)->rows + (r_1)), c_1) = _slang_vector_get_element(x_7.rows[c_1], r_1);
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_7;
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
    float3  result_8;
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
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_3], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_8, i_3) = sum_6;
        i_3 = i_3 + int(1);
    }
    return result_8;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_6)
{
    float _S14 = (F32_exp(((*dpx_5).primal_0))) * dOut_6;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S14;
    return;
}

inline __device__ float4  normalize_0(float4  x_8)
{
    return x_8 / make_float4 (length_2(x_8));
}

inline __device__ float3  normalize_1(float3  x_9)
{
    return x_9 / make_float3 (length_1(x_9));
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

inline __device__ float3  max_0(float3  x_10, float3  y_3)
{
    float3  result_9;
    int i_4 = int(0);
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_9, i_4) = (F32_max((_slang_vector_get_element(x_10, i_4)), (_slang_vector_get_element(y_3, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_9;
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

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_0, float4  quat_0, float3  scale_0, float2  hardness_0, FixedArray<float3 , 16>  sh_coeffs_0, FixedArray<float3 , 2>  ch_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * depth_0, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_0, mean_0) + t_0;
        bool _S31 = (mean_c_0.z) <= 0.0f;
        if(_S31)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float _S32 = scale_0.x;
        float sx_0 = (F32_exp((_S32)));
        float _S33 = scale_0.y;
        float sy_0 = (F32_exp((_S33)));
        float sz_0 = scale_0.z - 0.5f * (_S32 + _S33);
        float4  _S34 = normalize_0(quat_0);
        float x_11 = _S34.y;
        float x2_0 = x_11 * x_11;
        float y2_0 = _S34.z * _S34.z;
        float z2_0 = _S34.w * _S34.w;
        float xy_0 = _S34.y * _S34.z;
        float xz_0 = _S34.y * _S34.w;
        float yz_0 = _S34.z * _S34.w;
        float wx_0 = _S34.x * _S34.y;
        float wy_0 = _S34.x * _S34.z;
        float wz_0 = _S34.x * _S34.w;
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
        bool _S39;
        if(_S36 <= 0.0f)
        {
            _S39 = true;
        }
        else
        {
            _S39 = _S37 <= 0.0f;
        }
        if(_S39)
        {
            _S39 = true;
        }
        else
        {
            _S39 = _S38 <= 0.0f;
        }
        if(_S39)
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
                _S39 = true;
            }
            else
            {
                float u_0 = uv0_1.x;
                float v_0 = uv0_1.y;
                float _S40 = 0.0f * v_0;
                float r2_0 = u_0 * u_0 + v_0 * v_0;
                float s_diff_r2_0 = u_0 + u_0 + (_S40 + _S40);
                float _S41 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                float _S42 = dist_coeffs_0[int(1)] + r2_0 * _S41;
                float _S43 = dist_coeffs_0[int(0)] + r2_0 * _S42;
                float radial_0 = 1.0f + r2_0 * _S43;
                float _S44 = 2.0f * dist_coeffs_0[int(4)];
                float _S45 = _S44 * u_0;
                float _S46 = 2.0f * u_0;
                float _S47 = 2.0f * dist_coeffs_0[int(5)];
                float _S48 = _S47 * u_0;
                float _S49 = 2.0f * v_0;
                float2  _S50 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S43 + (s_diff_r2_0 * _S42 + (s_diff_r2_0 * _S41 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv0_1 + make_float2 (_S44 * v_0 + 0.0f * _S45 + (s_diff_r2_0 + (_S46 + _S46)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S47 * v_0 + 0.0f * _S48 + (s_diff_r2_0 + (_S40 + 0.0f * _S49)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
                float _S51 = 0.0f * u_0;
                float s_diff_r2_1 = _S51 + _S51 + (v_0 + v_0);
                float2  _S52 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S43 + (s_diff_r2_1 * _S42 + (s_diff_r2_1 * _S41 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv0_1 + make_float2 (0.0f * _S44 * v_0 + _S45 + (s_diff_r2_1 + (_S51 + 0.0f * _S46)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S47 * v_0 + _S48 + (s_diff_r2_1 + (_S49 + _S49)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S53 = transpose_0(makeMatrix<float, 2, 2> (_S50 + make_float2 (_S50.x * dist_coeffs_0[int(8)] + _S50.y * dist_coeffs_0[int(9)], 0.0f), _S52 + make_float2 (_S52.x * dist_coeffs_0[int(8)] + _S52.y * dist_coeffs_0[int(9)], 0.0f)));
                _S39 = !((F32_min((determinant_0(_S53)), ((F32_min((_S53.rows[int(0)].x), (_S53.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S39)
            {
                uv0_0 = uv0_1;
                break;
            }
            float u_1 = uv0_1.x;
            float v_1 = uv0_1.y;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float2  _S54 = uv0_1 * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
            float2  _S55 = _S54 + make_float2 (dist_coeffs_0[int(8)] * _S54.x + dist_coeffs_0[int(9)] * _S54.y, 0.0f);
            uv0_0 = make_float2 (fx_0 * _S55.x + cx_0, fy_0 * _S55.y + cy_0);
            break;
        }
        float2  uv1_0;
        bool all_valid_0 = true & (!_S39);
        for(;;)
        {
            float2  uv1_1 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S37);
            if(_S37 < 0.0f)
            {
                _S39 = true;
            }
            else
            {
                float u_2 = uv1_1.x;
                float v_2 = uv1_1.y;
                float _S56 = 0.0f * v_2;
                float r2_2 = u_2 * u_2 + v_2 * v_2;
                float s_diff_r2_2 = u_2 + u_2 + (_S56 + _S56);
                float _S57 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
                float _S58 = dist_coeffs_0[int(1)] + r2_2 * _S57;
                float _S59 = dist_coeffs_0[int(0)] + r2_2 * _S58;
                float radial_1 = 1.0f + r2_2 * _S59;
                float _S60 = 2.0f * dist_coeffs_0[int(4)];
                float _S61 = _S60 * u_2;
                float _S62 = 2.0f * u_2;
                float _S63 = 2.0f * dist_coeffs_0[int(5)];
                float _S64 = _S63 * u_2;
                float _S65 = 2.0f * v_2;
                float2  _S66 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S59 + (s_diff_r2_2 * _S58 + (s_diff_r2_2 * _S57 + s_diff_r2_2 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv1_1 + make_float2 (_S60 * v_2 + 0.0f * _S61 + (s_diff_r2_2 + (_S62 + _S62)) * dist_coeffs_0[int(5)] + s_diff_r2_2 * dist_coeffs_0[int(6)], _S63 * v_2 + 0.0f * _S64 + (s_diff_r2_2 + (_S56 + 0.0f * _S65)) * dist_coeffs_0[int(4)] + s_diff_r2_2 * dist_coeffs_0[int(7)]);
                float _S67 = 0.0f * u_2;
                float s_diff_r2_3 = _S67 + _S67 + (v_2 + v_2);
                float2  _S68 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S59 + (s_diff_r2_3 * _S58 + (s_diff_r2_3 * _S57 + s_diff_r2_3 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv1_1 + make_float2 (0.0f * _S60 * v_2 + _S61 + (s_diff_r2_3 + (_S67 + 0.0f * _S62)) * dist_coeffs_0[int(5)] + s_diff_r2_3 * dist_coeffs_0[int(6)], 0.0f * _S63 * v_2 + _S64 + (s_diff_r2_3 + (_S65 + _S65)) * dist_coeffs_0[int(4)] + s_diff_r2_3 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S69 = transpose_0(makeMatrix<float, 2, 2> (_S66 + make_float2 (_S66.x * dist_coeffs_0[int(8)] + _S66.y * dist_coeffs_0[int(9)], 0.0f), _S68 + make_float2 (_S68.x * dist_coeffs_0[int(8)] + _S68.y * dist_coeffs_0[int(9)], 0.0f)));
                _S39 = !((F32_min((determinant_0(_S69)), ((F32_min((_S69.rows[int(0)].x), (_S69.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S39)
            {
                uv1_0 = uv1_1;
                break;
            }
            float u_3 = uv1_1.x;
            float v_3 = uv1_1.y;
            float r2_3 = u_3 * u_3 + v_3 * v_3;
            float2  _S70 = uv1_1 * make_float2 (1.0f + r2_3 * (dist_coeffs_0[int(0)] + r2_3 * (dist_coeffs_0[int(1)] + r2_3 * (dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_3 * v_3 + dist_coeffs_0[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + dist_coeffs_0[int(6)] * r2_3, 2.0f * dist_coeffs_0[int(5)] * u_3 * v_3 + dist_coeffs_0[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + dist_coeffs_0[int(7)] * r2_3);
            float2  _S71 = _S70 + make_float2 (dist_coeffs_0[int(8)] * _S70.x + dist_coeffs_0[int(9)] * _S70.y, 0.0f);
            uv1_0 = make_float2 (fx_0 * _S71.x + cx_0, fy_0 * _S71.y + cy_0);
            break;
        }
        float2  uv2_0;
        bool all_valid_1 = all_valid_0 & (!_S39);
        for(;;)
        {
            float2  uv2_1 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S38);
            if(_S38 < 0.0f)
            {
                _S39 = true;
            }
            else
            {
                float u_4 = uv2_1.x;
                float v_4 = uv2_1.y;
                float _S72 = 0.0f * v_4;
                float r2_4 = u_4 * u_4 + v_4 * v_4;
                float s_diff_r2_4 = u_4 + u_4 + (_S72 + _S72);
                float _S73 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
                float _S74 = dist_coeffs_0[int(1)] + r2_4 * _S73;
                float _S75 = dist_coeffs_0[int(0)] + r2_4 * _S74;
                float radial_2 = 1.0f + r2_4 * _S75;
                float _S76 = 2.0f * dist_coeffs_0[int(4)];
                float _S77 = _S76 * u_4;
                float _S78 = 2.0f * u_4;
                float _S79 = 2.0f * dist_coeffs_0[int(5)];
                float _S80 = _S79 * u_4;
                float _S81 = 2.0f * v_4;
                float2  _S82 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S75 + (s_diff_r2_4 * _S74 + (s_diff_r2_4 * _S73 + s_diff_r2_4 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv2_1 + make_float2 (_S76 * v_4 + 0.0f * _S77 + (s_diff_r2_4 + (_S78 + _S78)) * dist_coeffs_0[int(5)] + s_diff_r2_4 * dist_coeffs_0[int(6)], _S79 * v_4 + 0.0f * _S80 + (s_diff_r2_4 + (_S72 + 0.0f * _S81)) * dist_coeffs_0[int(4)] + s_diff_r2_4 * dist_coeffs_0[int(7)]);
                float _S83 = 0.0f * u_4;
                float s_diff_r2_5 = _S83 + _S83 + (v_4 + v_4);
                float2  _S84 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S75 + (s_diff_r2_5 * _S74 + (s_diff_r2_5 * _S73 + s_diff_r2_5 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv2_1 + make_float2 (0.0f * _S76 * v_4 + _S77 + (s_diff_r2_5 + (_S83 + 0.0f * _S78)) * dist_coeffs_0[int(5)] + s_diff_r2_5 * dist_coeffs_0[int(6)], 0.0f * _S79 * v_4 + _S80 + (s_diff_r2_5 + (_S81 + _S81)) * dist_coeffs_0[int(4)] + s_diff_r2_5 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S85 = transpose_0(makeMatrix<float, 2, 2> (_S82 + make_float2 (_S82.x * dist_coeffs_0[int(8)] + _S82.y * dist_coeffs_0[int(9)], 0.0f), _S84 + make_float2 (_S84.x * dist_coeffs_0[int(8)] + _S84.y * dist_coeffs_0[int(9)], 0.0f)));
                _S39 = !((F32_min((determinant_0(_S85)), ((F32_min((_S85.rows[int(0)].x), (_S85.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S39)
            {
                uv2_0 = uv2_1;
                break;
            }
            float u_5 = uv2_1.x;
            float v_5 = uv2_1.y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S86 = uv2_1 * make_float2 (1.0f + r2_5 * (dist_coeffs_0[int(0)] + r2_5 * (dist_coeffs_0[int(1)] + r2_5 * (dist_coeffs_0[int(2)] + r2_5 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_5 * v_5 + dist_coeffs_0[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + dist_coeffs_0[int(6)] * r2_5, 2.0f * dist_coeffs_0[int(5)] * u_5 * v_5 + dist_coeffs_0[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + dist_coeffs_0[int(7)] * r2_5);
            float2  _S87 = _S86 + make_float2 (dist_coeffs_0[int(8)] * _S86.x + dist_coeffs_0[int(9)] * _S86.y, 0.0f);
            uv2_0 = make_float2 (fx_0 * _S87.x + cx_0, fy_0 * _S87.y + cy_0);
            break;
        }
        if(!(all_valid_1 & (!_S39)))
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float2  e0_0 = uv1_0 - uv0_0;
        float2  e1_0 = uv2_0 - uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(uv0_0 - uv2_0)));
        float _S88 = uv0_0.x;
        float _S89 = uv1_0.x;
        float _S90 = uv2_0.x;
        float xmax_0 = (F32_max(((F32_max((_S88), (_S89)))), (_S90))) + offset_0;
        float xmin_0 = (F32_min(((F32_min((_S88), (_S89)))), (_S90))) - offset_0;
        float _S91 = uv0_0.y;
        float _S92 = uv1_0.y;
        float _S93 = uv2_0.y;
        float ymax_0 = (F32_max(((F32_max((_S91), (_S92)))), (_S93))) + offset_0;
        float ymin_0 = (F32_min(((F32_min((_S91), (_S92)))), (_S93))) - offset_0;
        if(xmax_0 <= 0.0f)
        {
            _S39 = true;
        }
        else
        {
            _S39 = xmin_0 >= float(image_width_0);
        }
        if(_S39)
        {
            _S39 = true;
        }
        else
        {
            _S39 = ymax_0 <= 0.0f;
        }
        if(_S39)
        {
            _S39 = true;
        }
        else
        {
            _S39 = ymin_0 >= float(image_height_0);
        }
        if(_S39)
        {
            _S39 = true;
        }
        else
        {
            if(_S31)
            {
                if(xmin_0 <= 0.0f)
                {
                    _S39 = xmax_0 >= float(image_width_0);
                }
                else
                {
                    _S39 = false;
                }
                if(_S39)
                {
                    _S39 = true;
                }
                else
                {
                    if(ymin_0 <= 0.0f)
                    {
                        _S39 = ymax_0 >= float(image_width_0);
                    }
                    else
                    {
                        _S39 = false;
                    }
                }
            }
            else
            {
                _S39 = false;
            }
        }
        if(_S39)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (xmin_0, ymin_0, xmax_0, ymax_0);
        float3  _S94 = (vert0_c_0 + vert1_c_0 + vert2_c_0) / make_float3 (3.0f);
        *depth_0 = _S94.z;
        float3  _S95 = mean_0 - - mul_0(transpose_1(R_0), t_0);
        float _S96 = _S95.x;
        float _S97 = _S95.y;
        float _S98 = _S95.z;
        float norm_0 = (F32_sqrt((_S96 * _S96 + _S97 * _S97 + _S98 * _S98)));
        float x_12 = _S96 / norm_0;
        float y_4 = _S97 / norm_0;
        float z_0 = _S98 / norm_0;
        float z2_1 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_12 * x_12 - y_4 * y_4;
        float fS1_0 = 2.0f * x_12 * y_4;
        float fTmp0C_0 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        float3  color_0 = make_float3 (0.282094806432724f) * sh_coeffs_0[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * sh_coeffs_0[int(1)] + make_float3 (z_0) * sh_coeffs_0[int(2)] - make_float3 (x_12) * sh_coeffs_0[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * sh_coeffs_0[int(4)] + make_float3 (fTmp0B_0 * y_4) * sh_coeffs_0[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * sh_coeffs_0[int(6)] + make_float3 (fTmp0B_0 * x_12) * sh_coeffs_0[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * sh_coeffs_0[int(8)]) + (make_float3 (-0.59004360437393188f * (x_12 * fS1_0 + y_4 * fC1_0)) * sh_coeffs_0[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * sh_coeffs_0[int(10)] + make_float3 (fTmp0C_0 * y_4) * sh_coeffs_0[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * sh_coeffs_0[int(12)] + make_float3 (fTmp0C_0 * x_12) * sh_coeffs_0[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * sh_coeffs_0[int(14)] + make_float3 (-0.59004360437393188f * (x_12 * fC1_0 - y_4 * fS1_0)) * sh_coeffs_0[int(15)]);
        float3  _S99 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_0 + ch_coeffs_0[int(0)] + make_float3 (0.5f), _S99);
        float3  _S100 = color_0 - ch_coeffs_0[int(0)] * make_float3 (0.5f);
        float3  _S101 = ch_coeffs_0[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S100 + _S101 + make_float3 (0.5f), _S99);
        (*rgbs_0)[int(2)] = max_0(_S100 - _S101 + make_float3 (0.5f), _S99);
        (*verts_0)[int(0)] = vert0_0;
        (*verts_0)[int(1)] = vert1_0;
        (*verts_0)[int(2)] = vert2_0;
        float3  _S102 = normalize_1(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S102 * make_float3 (float(- (F32_sign((dot_0(_S102, mean_c_0))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_1, float4  quat_1, float3  scale_1, float2  hardness_1, FixedArray<float3 , 16>  sh_coeffs_1, FixedArray<float3 , 2>  ch_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * depth_1, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_1)
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
        if((length_1(mean_c_1)) <= 0.0f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S122 = scale_1.x;
        float sx_1 = (F32_exp((_S122)));
        float _S123 = scale_1.y;
        float sy_1 = (F32_exp((_S123)));
        float sz_1 = scale_1.z - 0.5f * (_S122 + _S123);
        float4  _S124 = normalize_0(quat_1);
        float x_13 = _S124.y;
        float x2_1 = x_13 * x_13;
        float y2_1 = _S124.z * _S124.z;
        float z2_2 = _S124.w * _S124.w;
        float xy_1 = _S124.y * _S124.z;
        float xz_1 = _S124.y * _S124.w;
        float yz_1 = _S124.z * _S124.w;
        float wx_1 = _S124.x * _S124.y;
        float wy_1 = _S124.x * _S124.z;
        float wz_1 = _S124.x * _S124.w;
        Matrix<float, 3, 3>  _S125 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_2), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_2), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1)));
        float3  vert0_1 = mul_0(_S125, make_float3 (sx_1, 0.0f, 0.0f)) + mean_1;
        float3  vert1_1 = mul_0(_S125, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_1;
        float3  vert2_1 = mul_0(_S125, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_1;
        float3  vert0_c_1 = mul_0(R_1, vert0_1) + t_1;
        float3  vert1_c_1 = mul_0(R_1, vert1_1) + t_1;
        float3  vert2_c_1 = mul_0(R_1, vert2_1) + t_1;
        float _S126 = length_1(vert1_c_1);
        float _S127 = length_1(vert2_c_1);
        bool _S128;
        if((length_1(vert0_c_1)) <= 0.0f)
        {
            _S128 = true;
        }
        else
        {
            _S128 = _S126 <= 0.0f;
        }
        if(_S128)
        {
            _S128 = true;
        }
        else
        {
            _S128 = _S127 <= 0.0f;
        }
        if(_S128)
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
            _S103 = _S131;
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
            float _S132 = 0.0f * v_6;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float s_diff_r2_6 = u_6 + u_6 + (_S132 + _S132);
            float _S133 = dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)];
            float _S134 = dist_coeffs_1[int(1)] + r2_6 * _S133;
            float _S135 = dist_coeffs_1[int(0)] + r2_6 * _S134;
            float _S136 = s_diff_r2_6 * _S135 + (s_diff_r2_6 * _S134 + (s_diff_r2_6 * _S133 + s_diff_r2_6 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float radial_3 = 1.0f + r2_6 * _S135;
            float _S137 = 2.0f * dist_coeffs_1[int(4)];
            _S114 = _S137;
            float _S138 = _S137 * u_6;
            float _S139 = 2.0f * u_6;
            float s_diff_du_0 = _S137 * v_6 + 0.0f * _S138 + (s_diff_r2_6 + (_S139 + _S139)) * dist_coeffs_1[int(5)] + s_diff_r2_6 * dist_coeffs_1[int(6)];
            float _S140 = 2.0f * dist_coeffs_1[int(5)];
            _S115 = _S140;
            float _S141 = _S140 * u_6;
            float _S142 = 2.0f * v_6;
            float2  _S143 = _S131 * make_float2 (radial_3) + make_float2 (_S136) * uv0_3 + make_float2 (s_diff_du_0, _S140 * v_6 + 0.0f * _S141 + (s_diff_r2_6 + (_S132 + 0.0f * _S142)) * dist_coeffs_1[int(4)] + s_diff_r2_6 * dist_coeffs_1[int(7)]);
            float2  _S144 = _S143 + make_float2 (_S143.x * dist_coeffs_1[int(8)] + _S143.y * dist_coeffs_1[int(9)], 0.0f);
            float2  _S145 = make_float2 (0.0f, 1.0f);
            _S116 = _S145;
            float _S146 = 0.0f * u_6;
            float s_diff_r2_7 = _S146 + _S146 + (v_6 + v_6);
            float _S147 = s_diff_r2_7 * _S135 + (s_diff_r2_7 * _S134 + (s_diff_r2_7 * _S133 + s_diff_r2_7 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float _S148 = 0.0f * _S137;
            _S117 = _S148;
            float s_diff_du_1 = _S148 * v_6 + _S138 + (s_diff_r2_7 + (_S146 + 0.0f * _S139)) * dist_coeffs_1[int(5)] + s_diff_r2_7 * dist_coeffs_1[int(6)];
            float _S149 = 0.0f * _S140;
            _S118 = _S149;
            float2  _S150 = _S145 * make_float2 (radial_3) + make_float2 (_S147) * uv0_3 + make_float2 (s_diff_du_1, _S149 * v_6 + _S141 + (s_diff_r2_7 + (_S142 + _S142)) * dist_coeffs_1[int(4)] + s_diff_r2_7 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S151 = transpose_0(makeMatrix<float, 2, 2> (_S144, _S150 + make_float2 (_S150.x * dist_coeffs_1[int(8)] + _S150.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S152 = !((F32_min((determinant_0(_S151)), ((F32_min((_S151.rows[int(0)].x), (_S151.rows[int(1)].y)))))) > 0.0f);
            _S119 = _S152;
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
        bool all_valid_2 = true & (!_S119);
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
            float _S158 = _S106 + r2_7 * _S107;
            float _S159 = _S105 + r2_7 * _S158;
            float _S160 = _S104 + r2_7 * _S159;
            float radial_4 = 1.0f + r2_7 * _S160;
            float _S161 = _S114 * u_7;
            float _S162 = 2.0f * u_7;
            float _S163 = _S115 * u_7;
            float _S164 = 2.0f * v_7;
            float2  _S165 = _S103 * make_float2 (radial_4) + make_float2 (s_diff_r2_8 * _S160 + (s_diff_r2_8 * _S159 + (s_diff_r2_8 * _S158 + s_diff_r2_8 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S114 * v_7 + 0.0f * _S161 + (s_diff_r2_8 + (_S162 + _S162)) * _S109 + s_diff_r2_8 * _S110, _S115 * v_7 + 0.0f * _S163 + (s_diff_r2_8 + (_S157 + 0.0f * _S164)) * _S108 + s_diff_r2_8 * _S111);
            float _S166 = 0.0f * u_7;
            float s_diff_r2_9 = _S166 + _S166 + (v_7 + v_7);
            float2  _S167 = _S116 * make_float2 (radial_4) + make_float2 (s_diff_r2_9 * _S160 + (s_diff_r2_9 * _S159 + (s_diff_r2_9 * _S158 + s_diff_r2_9 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S117 * v_7 + _S161 + (s_diff_r2_9 + (_S166 + 0.0f * _S162)) * _S109 + s_diff_r2_9 * _S110, _S118 * v_7 + _S163 + (s_diff_r2_9 + (_S164 + _S164)) * _S108 + s_diff_r2_9 * _S111);
            Matrix<float, 2, 2>  _S168 = transpose_0(makeMatrix<float, 2, 2> (_S165 + make_float2 (_S165.x * _S112 + _S165.y * _S113, 0.0f), _S167 + make_float2 (_S167.x * _S112 + _S167.y * _S113, 0.0f)));
            bool _S169 = !((F32_min((determinant_0(_S168)), ((F32_min((_S168.rows[int(0)].x), (_S168.rows[int(1)].y)))))) > 0.0f);
            _S120 = _S169;
            if(_S169)
            {
                uv1_2 = uv1_3;
                break;
            }
            float2  _S170 = uv1_3 * make_float2 (radial_4) + make_float2 (_S161 * v_7 + _S109 * (r2_7 + _S162 * u_7) + _S110 * r2_7, _S163 * v_7 + _S108 * (r2_7 + _S164 * v_7) + _S111 * r2_7);
            float2  _S171 = _S170 + make_float2 (_S112 * _S170.x + _S113 * _S170.y, 0.0f);
            uv1_2 = make_float2 (fx_1 * _S171.x + cx_1, fy_1 * _S171.y + cy_1);
            break;
        }
        float2  uv2_2;
        bool all_valid_3 = all_valid_2 & (!_S120);
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
            float _S175 = _S106 + r2_8 * _S107;
            float _S176 = _S105 + r2_8 * _S175;
            float _S177 = _S104 + r2_8 * _S176;
            float radial_5 = 1.0f + r2_8 * _S177;
            float _S178 = _S114 * u_8;
            float _S179 = 2.0f * u_8;
            float _S180 = _S115 * u_8;
            float _S181 = 2.0f * v_8;
            float2  _S182 = _S103 * make_float2 (radial_5) + make_float2 (s_diff_r2_10 * _S177 + (s_diff_r2_10 * _S176 + (s_diff_r2_10 * _S175 + s_diff_r2_10 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S114 * v_8 + 0.0f * _S178 + (s_diff_r2_10 + (_S179 + _S179)) * _S109 + s_diff_r2_10 * _S110, _S115 * v_8 + 0.0f * _S180 + (s_diff_r2_10 + (_S174 + 0.0f * _S181)) * _S108 + s_diff_r2_10 * _S111);
            float _S183 = 0.0f * u_8;
            float s_diff_r2_11 = _S183 + _S183 + (v_8 + v_8);
            float2  _S184 = _S116 * make_float2 (radial_5) + make_float2 (s_diff_r2_11 * _S177 + (s_diff_r2_11 * _S176 + (s_diff_r2_11 * _S175 + s_diff_r2_11 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S117 * v_8 + _S178 + (s_diff_r2_11 + (_S183 + 0.0f * _S179)) * _S109 + s_diff_r2_11 * _S110, _S118 * v_8 + _S180 + (s_diff_r2_11 + (_S181 + _S181)) * _S108 + s_diff_r2_11 * _S111);
            Matrix<float, 2, 2>  _S185 = transpose_0(makeMatrix<float, 2, 2> (_S182 + make_float2 (_S182.x * _S112 + _S182.y * _S113, 0.0f), _S184 + make_float2 (_S184.x * _S112 + _S184.y * _S113, 0.0f)));
            bool _S186 = !((F32_min((determinant_0(_S185)), ((F32_min((_S185.rows[int(0)].x), (_S185.rows[int(1)].y)))))) > 0.0f);
            _S121 = _S186;
            if(_S186)
            {
                uv2_2 = uv2_3;
                break;
            }
            float2  _S187 = uv2_3 * make_float2 (radial_5) + make_float2 (_S178 * v_8 + _S109 * (r2_8 + _S179 * u_8) + _S110 * r2_8, _S180 * v_8 + _S108 * (r2_8 + _S181 * v_8) + _S111 * r2_8);
            float2  _S188 = _S187 + make_float2 (_S112 * _S187.x + _S113 * _S187.y, 0.0f);
            uv2_2 = make_float2 (fx_1 * _S188.x + cx_1, fy_1 * _S188.y + cy_1);
            break;
        }
        if(!(all_valid_3 & (!_S121)))
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
            _S128 = true;
        }
        else
        {
            _S128 = xmin_1 >= float(image_width_1);
        }
        if(_S128)
        {
            _S128 = true;
        }
        else
        {
            _S128 = ymax_1 <= 0.0f;
        }
        if(_S128)
        {
            _S128 = true;
        }
        else
        {
            _S128 = ymin_1 >= float(image_height_1);
        }
        if(_S128)
        {
            _S128 = true;
        }
        else
        {
            if((mean_c_1.z) <= 0.0f)
            {
                if(xmin_1 <= 0.0f)
                {
                    _S128 = xmax_1 >= float(image_width_1);
                }
                else
                {
                    _S128 = false;
                }
                if(_S128)
                {
                    _S128 = true;
                }
                else
                {
                    if(ymin_1 <= 0.0f)
                    {
                        _S128 = ymax_1 >= float(image_width_1);
                    }
                    else
                    {
                        _S128 = false;
                    }
                }
            }
            else
            {
                _S128 = false;
            }
        }
        if(_S128)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (xmin_1, ymin_1, xmax_1, ymax_1);
        float3  _S195 = (vert0_c_1 + vert1_c_1 + vert2_c_1) / make_float3 (3.0f);
        float x_14 = _S195.x;
        float y_5 = _S195.y;
        float z_1 = _S195.z;
        float _S196 = x_14 * x_14 + y_5 * y_5;
        *depth_1 = z_1 * z_1 * z_1 * z_1 + 0.001953125f * _S196 * _S196;
        float3  _S197 = mean_1 - - mul_0(transpose_1(R_1), t_1);
        float _S198 = _S197.x;
        float _S199 = _S197.y;
        float _S200 = _S197.z;
        float norm_1 = (F32_sqrt((_S198 * _S198 + _S199 * _S199 + _S200 * _S200)));
        float x_15 = _S198 / norm_1;
        float y_6 = _S199 / norm_1;
        float z_2 = _S200 / norm_1;
        float z2_3 = z_2 * z_2;
        float fTmp0B_1 = -1.09254848957061768f * z_2;
        float fC1_1 = x_15 * x_15 - y_6 * y_6;
        float fS1_1 = 2.0f * x_15 * y_6;
        float fTmp0C_1 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_2;
        float3  color_1 = make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_1[int(1)] + make_float3 (z_2) * sh_coeffs_1[int(2)] - make_float3 (x_15) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_6) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_15) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_1 + y_6 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_6) * sh_coeffs_1[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_15) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_1 - y_6 * fS1_1)) * sh_coeffs_1[int(15)]);
        float3  _S201 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_1 + ch_coeffs_1[int(0)] + make_float3 (0.5f), _S201);
        float3  _S202 = color_1 - ch_coeffs_1[int(0)] * make_float3 (0.5f);
        float3  _S203 = ch_coeffs_1[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S202 + _S203 + make_float3 (0.5f), _S201);
        (*rgbs_1)[int(2)] = max_0(_S202 - _S203 + make_float3 (0.5f), _S201);
        (*verts_1)[int(0)] = vert0_1;
        (*verts_1)[int(1)] = vert1_1;
        (*verts_1)[int(2)] = vert2_1;
        float3  _S204 = normalize_1(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S204 * make_float3 (float(- (F32_sign((dot_0(_S204, mean_c_1))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
{
    float3  mean_c_2 = mul_0(R_2, mean_2) + t_2;
    float _S205 = scale_2.x;
    float sx_2 = (F32_exp((_S205)));
    float _S206 = scale_2.y;
    float sy_2 = (F32_exp((_S206)));
    float sz_2 = scale_2.z - 0.5f * (_S205 + _S206);
    float4  _S207 = normalize_0(quat_2);
    float x_16 = _S207.y;
    float x2_2 = x_16 * x_16;
    float y2_2 = _S207.z * _S207.z;
    float z2_4 = _S207.w * _S207.w;
    float xy_2 = _S207.y * _S207.z;
    float xz_2 = _S207.y * _S207.w;
    float yz_2 = _S207.z * _S207.w;
    float wx_2 = _S207.x * _S207.y;
    float wy_2 = _S207.x * _S207.z;
    float wz_2 = _S207.x * _S207.w;
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
    *aabb_xyxy_2 = make_float4 ((F32_min(((F32_min((_S214), (_S219)))), (_S224))) - offset_2, (F32_min(((F32_min((_S215), (_S220)))), (_S225))) - offset_2, (F32_max(((F32_max((_S214), (_S219)))), (_S224))) + offset_2, (F32_max(((F32_max((_S215), (_S220)))), (_S225))) + offset_2);
    *depth_2 = ((vert0_c_2 + vert1_c_2 + vert2_c_2) / make_float3 (3.0f)).z;
    float3  _S226 = mean_2 - - mul_0(transpose_1(R_2), t_2);
    float _S227 = _S226.x;
    float _S228 = _S226.y;
    float _S229 = _S226.z;
    float norm_2 = (F32_sqrt((_S227 * _S227 + _S228 * _S228 + _S229 * _S229)));
    float x_17 = _S227 / norm_2;
    float y_7 = _S228 / norm_2;
    float z_3 = _S229 / norm_2;
    float z2_5 = z_3 * z_3;
    float fTmp0B_2 = -1.09254848957061768f * z_3;
    float fC1_2 = x_17 * x_17 - y_7 * y_7;
    float fS1_2 = 2.0f * x_17 * y_7;
    float fTmp0C_2 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_3;
    float3  color_2 = make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * sh_coeffs_2[int(1)] + make_float3 (z_3) * sh_coeffs_2[int(2)] - make_float3 (x_17) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_7) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_17) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_17 * fS1_2 + y_7 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_7) * sh_coeffs_2[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_17) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_2 - y_7 * fS1_2)) * sh_coeffs_2[int(15)]);
    float3  _S230 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_2 + ch_coeffs_2[int(0)] + make_float3 (0.5f), _S230);
    float3  _S231 = color_2 - ch_coeffs_2[int(0)] * make_float3 (0.5f);
    float3  _S232 = ch_coeffs_2[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S231 + _S232 + make_float3 (0.5f), _S230);
    (*rgbs_2)[int(2)] = max_0(_S231 - _S232 + make_float3 (0.5f), _S230);
    (*verts_2)[int(0)] = vert0_2;
    (*verts_2)[int(1)] = vert1_2;
    (*verts_2)[int(2)] = vert2_2;
    float3  _S233 = normalize_1(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S233 * make_float3 (float(- (F32_sign((dot_0(_S233, mean_c_2))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
{
    float3  mean_c_3 = mul_0(R_3, mean_3) + t_3;
    float _S234 = scale_3.x;
    float sx_3 = (F32_exp((_S234)));
    float _S235 = scale_3.y;
    float sy_3 = (F32_exp((_S235)));
    float sz_3 = scale_3.z - 0.5f * (_S234 + _S235);
    float4  _S236 = normalize_0(quat_3);
    float x_18 = _S236.y;
    float x2_3 = x_18 * x_18;
    float y2_3 = _S236.z * _S236.z;
    float z2_6 = _S236.w * _S236.w;
    float xy_3 = _S236.y * _S236.z;
    float xz_3 = _S236.y * _S236.w;
    float yz_3 = _S236.z * _S236.w;
    float wx_3 = _S236.x * _S236.y;
    float wy_3 = _S236.x * _S236.z;
    float wz_3 = _S236.x * _S236.w;
    Matrix<float, 3, 3>  _S237 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_6), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_6), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    float3  vert0_3 = mul_0(_S237, make_float3 (sx_3, 0.0f, 0.0f)) + mean_3;
    float3  vert1_3 = mul_0(_S237, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_3;
    float3  vert2_3 = mul_0(_S237, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_3;
    float3  vert0_c_3 = mul_0(R_3, vert0_3) + t_3;
    float3  vert1_c_3 = mul_0(R_3, vert1_3) + t_3;
    float3  vert2_c_3 = mul_0(R_3, vert2_3) + t_3;
    float2  _S238 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_5 = length_0(_S238);
    float _S239 = vert0_c_3.z;
    float theta_3 = (F32_atan2((r_5), (_S239)));
    float k_1;
    if(theta_3 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_3 * theta_3 / 3.0f) / _S239;
    }
    else
    {
        k_1 = theta_3 / r_5;
    }
    float2  _S240 = _S238 * make_float2 (k_1);
    float u_12 = _S240.x;
    float v_12 = _S240.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S241 = 2.0f * dist_coeffs_3[int(4)];
    float _S242 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S243 = _S240 * make_float2 (1.0f + r2_12 * (dist_coeffs_3[int(0)] + r2_12 * (dist_coeffs_3[int(1)] + r2_12 * (dist_coeffs_3[int(2)] + r2_12 * dist_coeffs_3[int(3)])))) + make_float2 (_S241 * u_12 * v_12 + dist_coeffs_3[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + dist_coeffs_3[int(6)] * r2_12, _S242 * u_12 * v_12 + dist_coeffs_3[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + dist_coeffs_3[int(7)] * r2_12);
    float2  _S244 = _S243 + make_float2 (dist_coeffs_3[int(8)] * _S243.x + dist_coeffs_3[int(9)] * _S243.y, 0.0f);
    float _S245 = fx_3 * _S244.x + cx_3;
    float _S246 = fy_3 * _S244.y + cy_3;
    float2  uv0_5 = make_float2 (_S245, _S246);
    float2  _S247 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_6 = length_0(_S247);
    float _S248 = vert1_c_3.z;
    float theta_4 = (F32_atan2((r_6), (_S248)));
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S248;
    }
    else
    {
        k_1 = theta_4 / r_6;
    }
    float2  _S249 = _S247 * make_float2 (k_1);
    float u_13 = _S249.x;
    float v_13 = _S249.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S250 = _S249 * make_float2 (1.0f + r2_13 * (dist_coeffs_3[int(0)] + r2_13 * (dist_coeffs_3[int(1)] + r2_13 * (dist_coeffs_3[int(2)] + r2_13 * dist_coeffs_3[int(3)])))) + make_float2 (_S241 * u_13 * v_13 + dist_coeffs_3[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_3[int(6)] * r2_13, _S242 * u_13 * v_13 + dist_coeffs_3[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_3[int(7)] * r2_13);
    float2  _S251 = _S250 + make_float2 (dist_coeffs_3[int(8)] * _S250.x + dist_coeffs_3[int(9)] * _S250.y, 0.0f);
    float _S252 = fx_3 * _S251.x + cx_3;
    float _S253 = fy_3 * _S251.y + cy_3;
    float2  uv1_5 = make_float2 (_S252, _S253);
    float2  _S254 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_7 = length_0(_S254);
    float _S255 = vert2_c_3.z;
    float theta_5 = (F32_atan2((r_7), (_S255)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S255;
    }
    else
    {
        k_1 = theta_5 / r_7;
    }
    float2  _S256 = _S254 * make_float2 (k_1);
    float u_14 = _S256.x;
    float v_14 = _S256.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S257 = _S256 * make_float2 (1.0f + r2_14 * (dist_coeffs_3[int(0)] + r2_14 * (dist_coeffs_3[int(1)] + r2_14 * (dist_coeffs_3[int(2)] + r2_14 * dist_coeffs_3[int(3)])))) + make_float2 (_S241 * u_14 * v_14 + dist_coeffs_3[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + dist_coeffs_3[int(6)] * r2_14, _S242 * u_14 * v_14 + dist_coeffs_3[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + dist_coeffs_3[int(7)] * r2_14);
    float2  _S258 = _S257 + make_float2 (dist_coeffs_3[int(8)] * _S257.x + dist_coeffs_3[int(9)] * _S257.y, 0.0f);
    float _S259 = fx_3 * _S258.x + cx_3;
    float _S260 = fy_3 * _S258.y + cy_3;
    float2  uv2_5 = make_float2 (_S259, _S260);
    float2  e0_3 = uv1_5 - uv0_5;
    float2  e1_3 = uv2_5 - uv1_5;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(uv0_5 - uv2_5)));
    *aabb_xyxy_3 = make_float4 ((F32_min(((F32_min((_S245), (_S252)))), (_S259))) - offset_3, (F32_min(((F32_min((_S246), (_S253)))), (_S260))) - offset_3, (F32_max(((F32_max((_S245), (_S252)))), (_S259))) + offset_3, (F32_max(((F32_max((_S246), (_S253)))), (_S260))) + offset_3);
    float3  _S261 = (vert0_c_3 + vert1_c_3 + vert2_c_3) / make_float3 (3.0f);
    float x_19 = _S261.x;
    float y_8 = _S261.y;
    float z_4 = _S261.z;
    float _S262 = x_19 * x_19 + y_8 * y_8;
    *depth_3 = z_4 * z_4 * z_4 * z_4 + 0.001953125f * _S262 * _S262;
    float3  _S263 = mean_3 - - mul_0(transpose_1(R_3), t_3);
    float _S264 = _S263.x;
    float _S265 = _S263.y;
    float _S266 = _S263.z;
    float norm_3 = (F32_sqrt((_S264 * _S264 + _S265 * _S265 + _S266 * _S266)));
    float x_20 = _S264 / norm_3;
    float y_9 = _S265 / norm_3;
    float z_5 = _S266 / norm_3;
    float z2_7 = z_5 * z_5;
    float fTmp0B_3 = -1.09254848957061768f * z_5;
    float fC1_3 = x_20 * x_20 - y_9 * y_9;
    float fS1_3 = 2.0f * x_20 * y_9;
    float fTmp0C_3 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_5;
    float3  color_3 = make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * sh_coeffs_3[int(1)] + make_float3 (z_5) * sh_coeffs_3[int(2)] - make_float3 (x_20) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_9) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_20) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_20 * fS1_3 + y_9 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_9) * sh_coeffs_3[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_20) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_3 - y_9 * fS1_3)) * sh_coeffs_3[int(15)]);
    float3  _S267 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_3 + ch_coeffs_3[int(0)] + make_float3 (0.5f), _S267);
    float3  _S268 = color_3 - ch_coeffs_3[int(0)] * make_float3 (0.5f);
    float3  _S269 = ch_coeffs_3[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S268 + _S269 + make_float3 (0.5f), _S267);
    (*rgbs_3)[int(2)] = max_0(_S268 - _S269 + make_float3 (0.5f), _S267);
    (*verts_3)[int(0)] = vert0_3;
    (*verts_3)[int(1)] = vert1_3;
    (*verts_3)[int(2)] = vert2_3;
    float3  _S270 = normalize_1(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S270 * make_float3 (float(- (F32_sign((dot_0(_S270, mean_c_3))))));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S271, float3  _S272)
{
    return mul_0(_S271, _S272);
}

inline __device__ float s_primal_ctx_exp_0(float _S273)
{
    return (F32_exp((_S273)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S274)
{
    return (F32_sqrt((_S274)));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S275, float3  _S276)
{
    return cross_0(_S275, _S276);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S277, float3  _S278)
{
    return dot_0(_S277, _S278);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S279, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S280, float _S281)
{
    _d_dot_0(_S279, _S280, _S281);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S282, float _S283)
{
    _d_sqrt_0(_S282, _S283);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S284 = (*dpx_9).primal_0.x;
    float _S285 = (*dpx_9).primal_0.y;
    float _S286 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S287;
    (&_S287)->primal_0 = _S284 * _S284 + _S285 * _S285 + _S286 * _S286;
    (&_S287)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S287, _s_dOut_0);
    float _S288 = (*dpx_9).primal_0.z * _S287.differential_0;
    float _S289 = _S288 + _S288;
    float _S290 = (*dpx_9).primal_0.y * _S287.differential_0;
    float _S291 = _S290 + _S290;
    float _S292 = (*dpx_9).primal_0.x * _S287.differential_0;
    float _S293 = _S292 + _S292;
    float3  _S294 = make_float3 (0.0f);
    *&((&_S294)->z) = _S289;
    *&((&_S294)->y) = _S291;
    *&((&_S294)->x) = _S293;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S294;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S295, float _S296)
{
    s_bwd_prop_length_impl_0(_S295, _S296);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_1)
{
    float _S297 = length_1((*dpx_10).primal_0);
    float3  _S298 = (*dpx_10).primal_0 * _s_dOut_1;
    float3  _S299 = make_float3 (1.0f / _S297) * _s_dOut_1;
    float _S300 = - ((_S298.x + _S298.y + _S298.z) / (_S297 * _S297));
    float3  _S301 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S302;
    (&_S302)->primal_0 = (*dpx_10).primal_0;
    (&_S302)->differential_0 = _S301;
    s_bwd_length_impl_0(&_S302, _S300);
    float3  _S303 = _S299 + _S302.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S303;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S304, float3  _S305)
{
    s_bwd_prop_normalize_impl_0(_S304, _S305);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S306, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S307, float3  _S308)
{
    _d_cross_0(_S306, _S307, _S308);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S309, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S310, float3  _S311)
{
    _d_max_vector_0(_S309, _S310, _S311);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S312, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S313, float3  _S314)
{
    _d_mul_0(_S312, _S313, _S314);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S315, float _S316)
{
    _d_exp2_0(_S315, _S316);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float _s_dOut_2)
{
    float _S317 = (*dpx_11).primal_0.x;
    float _S318 = (*dpx_11).primal_0.y;
    DiffPair_float_0 _S319;
    (&_S319)->primal_0 = _S317 * _S317 + _S318 * _S318;
    (&_S319)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S319, _s_dOut_2);
    float _S320 = (*dpx_11).primal_0.y * _S319.differential_0;
    float _S321 = _S320 + _S320;
    float _S322 = (*dpx_11).primal_0.x * _S319.differential_0;
    float _S323 = _S322 + _S322;
    float2  _S324 = make_float2 (0.0f);
    *&((&_S324)->y) = _S321;
    *&((&_S324)->x) = _S323;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S324;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S325, float _S326)
{
    s_bwd_prop_length_impl_1(_S325, _S326);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S327, float _S328)
{
    _d_abs_0(_S327, _S328);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_12, float _s_dOut_3)
{
    float _S329 = (*dpx_12).primal_0.x;
    float _S330 = (*dpx_12).primal_0.y;
    float _S331 = (*dpx_12).primal_0.z;
    float _S332 = (*dpx_12).primal_0.w;
    DiffPair_float_0 _S333;
    (&_S333)->primal_0 = _S329 * _S329 + _S330 * _S330 + _S331 * _S331 + _S332 * _S332;
    (&_S333)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S333, _s_dOut_3);
    float _S334 = (*dpx_12).primal_0.w * _S333.differential_0;
    float _S335 = _S334 + _S334;
    float _S336 = (*dpx_12).primal_0.z * _S333.differential_0;
    float _S337 = _S336 + _S336;
    float _S338 = (*dpx_12).primal_0.y * _S333.differential_0;
    float _S339 = _S338 + _S338;
    float _S340 = (*dpx_12).primal_0.x * _S333.differential_0;
    float _S341 = _S340 + _S340;
    float4  _S342 = make_float4 (0.0f);
    *&((&_S342)->w) = _S335;
    *&((&_S342)->z) = _S337;
    *&((&_S342)->y) = _S339;
    *&((&_S342)->x) = _S341;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S342;
    return;
}

inline __device__ void s_bwd_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S343, float _S344)
{
    s_bwd_prop_length_impl_2(_S343, _S344);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_13, float4  _s_dOut_4)
{
    float _S345 = length_2((*dpx_13).primal_0);
    float4  _S346 = (*dpx_13).primal_0 * _s_dOut_4;
    float4  _S347 = make_float4 (1.0f / _S345) * _s_dOut_4;
    float _S348 = - ((_S346.x + _S346.y + _S346.z + _S346.w) / (_S345 * _S345));
    float4  _S349 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S350;
    (&_S350)->primal_0 = (*dpx_13).primal_0;
    (&_S350)->differential_0 = _S349;
    s_bwd_length_impl_2(&_S350, _S348);
    float4  _S351 = _S347 + _S350.differential_0;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S351;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S352, float4  _S353)
{
    s_bwd_prop_normalize_impl_1(_S352, _S353);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S354, float _S355)
{
    _d_exp_0(_S354, _S355);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S356 = scale_4.x;
    float _S357 = s_primal_ctx_exp_0(_S356);
    float _S358 = scale_4.y;
    float _S359 = s_primal_ctx_exp_0(_S358);
    float sz_4 = scale_4.z - 0.5f * (_S356 + _S358);
    float4  _S360 = normalize_0(quat_4);
    float _S361 = _S360.y;
    float x2_4 = _S361 * _S361;
    float y2_4 = _S360.z * _S360.z;
    float z2_8 = _S360.w * _S360.w;
    float xy_4 = _S360.y * _S360.z;
    float xz_4 = _S360.y * _S360.w;
    float yz_4 = _S360.z * _S360.w;
    float wx_4 = _S360.x * _S360.y;
    float wy_4 = _S360.x * _S360.z;
    float wz_4 = _S360.x * _S360.w;
    Matrix<float, 3, 3>  _S362 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S363 = make_float3 (_S357, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S362, _S363) + mean_4;
    float _S364 = -0.5f + sz_4;
    float3  _S365 = make_float3 (_S357 * _S364, _S359, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S362, _S365) + mean_4;
    float _S366 = -0.5f - sz_4;
    float3  _S367 = make_float3 (_S357 * _S366, - _S359, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S362, _S367) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S368 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S369 = vert0_c_4.z;
    float2  _S370 = make_float2 (_S369);
    float2  _S371 = _S368 / make_float2 (_S369);
    float2  _S372 = make_float2 (_S369 * _S369);
    float u_15 = _S371.x;
    float v_15 = _S371.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S373 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S374 = dist_coeffs_4[int(1)] + r2_15 * _S373;
    float _S375 = dist_coeffs_4[int(0)] + r2_15 * _S374;
    float radial_6 = 1.0f + r2_15 * _S375;
    float _S376 = 2.0f * dist_coeffs_4[int(4)];
    float _S377 = _S376 * u_15;
    float _S378 = 2.0f * u_15;
    float _S379 = 2.0f * dist_coeffs_4[int(5)];
    float _S380 = _S379 * u_15;
    float _S381 = 2.0f * v_15;
    float2  _S382 = _S371 * make_float2 (radial_6) + make_float2 (_S377 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S378 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S380 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S381 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S383 = _S382 + make_float2 (dist_coeffs_4[int(8)] * _S382.x + dist_coeffs_4[int(9)] * _S382.y, 0.0f);
    float _S384 = fx_4 * _S383.x + cx_4;
    float _S385 = fy_4 * _S383.y + cy_4;
    float2  uv0_6 = make_float2 (_S384, _S385);
    float2  _S386 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S387 = vert1_c_4.z;
    float2  _S388 = make_float2 (_S387);
    float2  _S389 = _S386 / make_float2 (_S387);
    float2  _S390 = make_float2 (_S387 * _S387);
    float u_16 = _S389.x;
    float v_16 = _S389.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S391 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S392 = dist_coeffs_4[int(1)] + r2_16 * _S391;
    float _S393 = dist_coeffs_4[int(0)] + r2_16 * _S392;
    float radial_7 = 1.0f + r2_16 * _S393;
    float _S394 = _S376 * u_16;
    float _S395 = 2.0f * u_16;
    float _S396 = _S379 * u_16;
    float _S397 = 2.0f * v_16;
    float2  _S398 = _S389 * make_float2 (radial_7) + make_float2 (_S394 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S395 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S396 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S397 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S399 = _S398 + make_float2 (dist_coeffs_4[int(8)] * _S398.x + dist_coeffs_4[int(9)] * _S398.y, 0.0f);
    float _S400 = fx_4 * _S399.x + cx_4;
    float _S401 = fy_4 * _S399.y + cy_4;
    float2  uv1_6 = make_float2 (_S400, _S401);
    float2  _S402 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S403 = vert2_c_4.z;
    float2  _S404 = make_float2 (_S403);
    float2  _S405 = _S402 / make_float2 (_S403);
    float2  _S406 = make_float2 (_S403 * _S403);
    float u_17 = _S405.x;
    float v_17 = _S405.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S407 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S408 = dist_coeffs_4[int(1)] + r2_17 * _S407;
    float _S409 = dist_coeffs_4[int(0)] + r2_17 * _S408;
    float radial_8 = 1.0f + r2_17 * _S409;
    float _S410 = _S376 * u_17;
    float _S411 = 2.0f * u_17;
    float _S412 = _S379 * u_17;
    float _S413 = 2.0f * v_17;
    float2  _S414 = _S405 * make_float2 (radial_8) + make_float2 (_S410 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S411 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S412 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S413 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S415 = _S414 + make_float2 (dist_coeffs_4[int(8)] * _S414.x + dist_coeffs_4[int(9)] * _S414.y, 0.0f);
    float _S416 = fx_4 * _S415.x + cx_4;
    float _S417 = fy_4 * _S415.y + cy_4;
    float2  uv2_6 = make_float2 (_S416, _S417);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S418 = e0_4.x;
    float _S419 = e1_4.y;
    float _S420 = e0_4.y;
    float _S421 = e1_4.x;
    float _S422 = _S418 * _S419 - _S420 * _S421;
    float _S423 = 1.0f - hardness_4.y;
    float _S424 = -1.0f / _S423;
    float _S425 = _S423 * _S423;
    float _S426 = (F32_max((_S384), (_S400)));
    float _S427 = (F32_min((_S384), (_S400)));
    float _S428 = (F32_max((_S385), (_S401)));
    float _S429 = (F32_min((_S385), (_S401)));
    Matrix<float, 3, 3>  _S430 = transpose_1(R_4);
    float3  _S431 = mean_4 - - s_primal_ctx_mul_0(_S430, t_4);
    float _S432 = _S431.x;
    float _S433 = _S431.y;
    float _S434 = _S431.z;
    float _S435 = _S432 * _S432 + _S433 * _S433 + _S434 * _S434;
    float _S436 = s_primal_ctx_sqrt_0(_S435);
    float x_21 = _S432 / _S436;
    float3  _S437 = make_float3 (x_21);
    float _S438 = _S436 * _S436;
    float y_10 = _S433 / _S436;
    float z_6 = _S434 / _S436;
    float3  _S439 = make_float3 (z_6);
    float _S440 = - y_10;
    float3  _S441 = make_float3 (_S440);
    float z2_9 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_4 = x_21 * x_21 - y_10 * y_10;
    float _S442 = 2.0f * x_21;
    float fS1_4 = _S442 * y_10;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S443 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_21;
    float3  _S444 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_10;
    float3  _S445 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S446 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S447 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_6;
    float _S448 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_6 * _S448;
    float3  _S449 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_21;
    float3  _S450 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_10;
    float3  _S451 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S452 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S453 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_21 * fC1_4 - y_10 * fS1_4);
    float3  _S454 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_21 * fS1_4 + y_10 * fC1_4);
    float3  _S455 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S440) * sh_coeffs_4[int(1)] + make_float3 (z_6) * sh_coeffs_4[int(2)] - make_float3 (x_21) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]);
    float3  _S456 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S457 = make_float3 (0.0f);
    float3  _S458 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S459 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S460 = make_float3 (_S459);
    float3  _S461 = ch_coeffs_4[int(1)] * make_float3 (_S459);
    float3  _S462 = _S458 + _S461 + make_float3 (0.5f);
    float3  _S463 = _S458 - _S461 + make_float3 (0.5f);
    float3  _S464 = vert1_c_4 - vert0_c_4;
    float3  _S465 = vert2_c_4 - vert0_c_4;
    float3  _S466 = s_primal_ctx_cross_0(_S464, _S465);
    float3  _S467 = normalize_1(_S466);
    float3  _S468 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S467, mean_c_4)))))) * v_normal_0;
    float3  _S469 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S470;
    (&_S470)->primal_0 = _S467;
    (&_S470)->differential_0 = _S469;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S471;
    (&_S471)->primal_0 = mean_c_4;
    (&_S471)->differential_0 = _S469;
    s_bwd_prop_dot_0(&_S470, &_S471, 0.0f);
    float3  _S472 = _S468 + _S470.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S473;
    (&_S473)->primal_0 = _S466;
    (&_S473)->differential_0 = _S469;
    s_bwd_normalize_impl_0(&_S473, _S472);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S474;
    (&_S474)->primal_0 = _S464;
    (&_S474)->differential_0 = _S469;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S475;
    (&_S475)->primal_0 = _S465;
    (&_S475)->differential_0 = _S469;
    s_bwd_prop_cross_0(&_S474, &_S475, _S473.differential_0);
    float3  _S476 = - _S475.differential_0;
    float3  _S477 = - _S474.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S478;
    (&_S478)->primal_0 = _S463;
    (&_S478)->differential_0 = _S469;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S479;
    (&_S479)->primal_0 = _S457;
    (&_S479)->differential_0 = _S469;
    s_bwd_prop_max_0(&_S478, &_S479, v_rgbs_0[int(2)]);
    float3  _S480 = - _S478.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481;
    (&_S481)->primal_0 = _S462;
    (&_S481)->differential_0 = _S469;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S482;
    (&_S482)->primal_0 = _S457;
    (&_S482)->differential_0 = _S469;
    s_bwd_prop_max_0(&_S481, &_S482, v_rgbs_0[int(1)]);
    float3  _S483 = _S460 * (_S480 + _S481.differential_0);
    float3  _S484 = _S478.differential_0 + _S481.differential_0;
    float3  _S485 = make_float3 (0.5f) * - _S484;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S486;
    (&_S486)->primal_0 = _S456;
    (&_S486)->differential_0 = _S469;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S487;
    (&_S487)->primal_0 = _S457;
    (&_S487)->differential_0 = _S469;
    s_bwd_prop_max_0(&_S486, &_S487, v_rgbs_0[int(0)]);
    float3  _S488 = _S485 + _S486.differential_0;
    float3  _S489 = _S484 + _S486.differential_0;
    float3  _S490 = _S454 * _S489;
    float3  _S491 = sh_coeffs_4[int(15)] * _S489;
    float3  _S492 = _S452 * _S489;
    float3  _S493 = sh_coeffs_4[int(14)] * _S489;
    float3  _S494 = _S450 * _S489;
    float3  _S495 = sh_coeffs_4[int(13)] * _S489;
    float3  _S496 = _S449 * _S489;
    float3  _S497 = sh_coeffs_4[int(12)] * _S489;
    float3  _S498 = _S451 * _S489;
    float3  _S499 = sh_coeffs_4[int(11)] * _S489;
    float3  _S500 = _S453 * _S489;
    float3  _S501 = sh_coeffs_4[int(10)] * _S489;
    float3  _S502 = _S455 * _S489;
    float3  _S503 = sh_coeffs_4[int(9)] * _S489;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S503.x + _S503.y + _S503.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S491.x + _S491.y + _S491.z);
    float _S504 = _S501.x + _S501.y + _S501.z;
    float _S505 = _S493.x + _S493.y + _S493.z;
    float _S506 = _S499.x + _S499.y + _S499.z;
    float _S507 = _S495.x + _S495.y + _S495.z;
    float _S508 = _S497.x + _S497.y + _S497.z;
    float _S509 = - s_diff_fC2_T_0;
    float3  _S510 = _S446 * _S489;
    float3  _S511 = sh_coeffs_4[int(8)] * _S489;
    float3  _S512 = _S444 * _S489;
    float3  _S513 = sh_coeffs_4[int(7)] * _S489;
    float3  _S514 = _S443 * _S489;
    float3  _S515 = sh_coeffs_4[int(6)] * _S489;
    float3  _S516 = _S445 * _S489;
    float3  _S517 = sh_coeffs_4[int(5)] * _S489;
    float3  _S518 = _S447 * _S489;
    float3  _S519 = sh_coeffs_4[int(4)] * _S489;
    float _S520 = _S517.x + _S517.y + _S517.z;
    float _S521 = _S513.x + _S513.y + _S513.z;
    float _S522 = fTmp1B_4 * _S504 + x_21 * s_diff_fS2_T_0 + y_10 * _S509 + 0.54627424478530884f * (_S519.x + _S519.y + _S519.z);
    float _S523 = fTmp1B_4 * _S505 + y_10 * s_diff_fS2_T_0 + x_21 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S511.x + _S511.y + _S511.z);
    float _S524 = y_10 * - _S523;
    float _S525 = x_21 * _S523;
    float _S526 = z_6 * (1.86588168144226074f * (z_6 * _S508) + -2.28522896766662598f * (y_10 * _S506 + x_21 * _S507) + 0.94617468118667603f * (_S515.x + _S515.y + _S515.z));
    float3  _S527 = make_float3 (0.48860251903533936f) * _S489;
    float3  _S528 = - _S527;
    float3  _S529 = _S437 * _S528;
    float3  _S530 = sh_coeffs_4[int(3)] * _S528;
    float3  _S531 = _S439 * _S527;
    float3  _S532 = sh_coeffs_4[int(2)] * _S527;
    float3  _S533 = _S441 * _S527;
    float3  _S534 = sh_coeffs_4[int(1)] * _S527;
    float _S535 = (_S448 * _S508 + 1.44530570507049561f * (fS1_4 * _S504 + fC1_4 * _S505) + -1.09254848957061768f * (y_10 * _S520 + x_21 * _S521) + _S526 + _S526 + _S532.x + _S532.y + _S532.z) / _S438;
    float _S536 = _S436 * _S535;
    float _S537 = (fTmp0C_4 * _S506 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S509 + fTmp0B_4 * _S520 + _S442 * _S522 + _S524 + _S524 + - (_S534.x + _S534.y + _S534.z)) / _S438;
    float _S538 = _S436 * _S537;
    float _S539 = (fTmp0C_4 * _S507 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S521 + 2.0f * (y_10 * _S522) + _S525 + _S525 + _S530.x + _S530.y + _S530.z) / _S438;
    float _S540 = _S436 * _S539;
    float _S541 = _S434 * - _S535 + _S433 * - _S537 + _S432 * - _S539;
    DiffPair_float_0 _S542;
    (&_S542)->primal_0 = _S435;
    (&_S542)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S542, _S541);
    float _S543 = _S434 * _S542.differential_0;
    float _S544 = _S433 * _S542.differential_0;
    float _S545 = _S432 * _S542.differential_0;
    float3  _S546 = make_float3 (0.282094806432724f) * _S489;
    float3  _S547 = make_float3 (_S540 + _S545 + _S545, _S538 + _S544 + _S544, _S536 + _S543 + _S543);
    float3  _S548 = - - _S547;
    Matrix<float, 3, 3>  _S549 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S550;
    (&_S550)->primal_0 = _S430;
    (&_S550)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S551;
    (&_S551)->primal_0 = t_4;
    (&_S551)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S550, &_S551, _S548);
    Matrix<float, 3, 3>  _S552 = transpose_1(_S550.differential_0);
    float3  _S553 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_0);
    DiffPair_float_0 _S554;
    (&_S554)->primal_0 = _S429;
    (&_S554)->differential_0 = 0.0f;
    DiffPair_float_0 _S555;
    (&_S555)->primal_0 = _S417;
    (&_S555)->differential_0 = 0.0f;
    _d_min_0(&_S554, &_S555, 0.0f);
    DiffPair_float_0 _S556;
    (&_S556)->primal_0 = _S385;
    (&_S556)->differential_0 = 0.0f;
    DiffPair_float_0 _S557;
    (&_S557)->primal_0 = _S401;
    (&_S557)->differential_0 = 0.0f;
    _d_min_0(&_S556, &_S557, _S554.differential_0);
    DiffPair_float_0 _S558;
    (&_S558)->primal_0 = _S428;
    (&_S558)->differential_0 = 0.0f;
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = _S417;
    (&_S559)->differential_0 = 0.0f;
    _d_max_0(&_S558, &_S559, 0.0f);
    DiffPair_float_0 _S560;
    (&_S560)->primal_0 = _S385;
    (&_S560)->differential_0 = 0.0f;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S401;
    (&_S561)->differential_0 = 0.0f;
    _d_max_0(&_S560, &_S561, _S558.differential_0);
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S427;
    (&_S562)->differential_0 = 0.0f;
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = _S416;
    (&_S563)->differential_0 = 0.0f;
    _d_min_0(&_S562, &_S563, 0.0f);
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = _S384;
    (&_S564)->differential_0 = 0.0f;
    DiffPair_float_0 _S565;
    (&_S565)->primal_0 = _S400;
    (&_S565)->differential_0 = 0.0f;
    _d_min_0(&_S564, &_S565, _S562.differential_0);
    DiffPair_float_0 _S566;
    (&_S566)->primal_0 = _S426;
    (&_S566)->differential_0 = 0.0f;
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S416;
    (&_S567)->differential_0 = 0.0f;
    _d_max_0(&_S566, &_S567, 0.0f);
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S384;
    (&_S568)->differential_0 = 0.0f;
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S400;
    (&_S569)->differential_0 = 0.0f;
    _d_max_0(&_S568, &_S569, _S566.differential_0);
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S424;
    (&_S570)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S570, -0.0f);
    float _S571 = - (-1.0f * - (_S570.differential_0 / _S425));
    float2  _S572 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S573;
    (&_S573)->primal_0 = e2_0;
    (&_S573)->differential_0 = _S572;
    s_bwd_length_impl_1(&_S573, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S574;
    (&_S574)->primal_0 = e1_4;
    (&_S574)->differential_0 = _S572;
    s_bwd_length_impl_1(&_S574, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S575;
    (&_S575)->primal_0 = e0_4;
    (&_S575)->differential_0 = _S572;
    s_bwd_length_impl_1(&_S575, 0.0f);
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S422;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S576, -0.0f);
    float _S577 = - _S576.differential_0;
    float2  _S578 = _S574.differential_0 + make_float2 (_S420 * _S577, _S418 * _S576.differential_0);
    float2  _S579 = _S575.differential_0 + make_float2 (_S419 * _S576.differential_0, _S421 * _S577);
    float2  _S580 = - _S573.differential_0 + _S578;
    float _S581 = fx_4 * (_S563.differential_0 + _S567.differential_0 + _S580.x);
    float2  _S582 = make_float2 (_S581, fy_4 * (_S555.differential_0 + _S559.differential_0 + _S580.y)) + make_float2 (dist_coeffs_4[int(8)] * _S581, dist_coeffs_4[int(9)] * _S581);
    float2  _S583 = _S405 * _S582;
    float _S584 = dist_coeffs_4[int(4)] * _S582.y;
    float _S585 = dist_coeffs_4[int(5)] * _S582.x;
    float _S586 = _S583.x + _S583.y;
    float _S587 = r2_17 * _S586;
    float _S588 = r2_17 * _S587;
    float _S589 = dist_coeffs_4[int(7)] * _S582.y + _S584 + dist_coeffs_4[int(6)] * _S582.x + _S585 + _S409 * _S586 + _S408 * _S587 + _S407 * _S588 + dist_coeffs_4[int(3)] * (r2_17 * _S588);
    float _S590 = v_17 * _S589;
    float _S591 = u_17 * _S589;
    float2  _S592 = (make_float2 (radial_8) * _S582 + make_float2 (_S379 * (v_17 * _S582.y) + _S411 * _S585 + 2.0f * (u_17 * _S585) + _S376 * (v_17 * _S582.x) + _S591 + _S591, _S413 * _S584 + 2.0f * (v_17 * _S584) + _S412 * _S582.y + _S410 * _S582.x + _S590 + _S590)) / _S406;
    float2  _S593 = _S402 * - _S592;
    float2  _S594 = _S404 * _S592;
    float2  _S595 = - _S578 + _S579;
    float _S596 = fx_4 * (_S565.differential_0 + _S569.differential_0 + _S595.x);
    float2  _S597 = make_float2 (_S596, fy_4 * (_S557.differential_0 + _S561.differential_0 + _S595.y)) + make_float2 (dist_coeffs_4[int(8)] * _S596, dist_coeffs_4[int(9)] * _S596);
    float2  _S598 = _S389 * _S597;
    float _S599 = dist_coeffs_4[int(4)] * _S597.y;
    float _S600 = dist_coeffs_4[int(5)] * _S597.x;
    float _S601 = _S598.x + _S598.y;
    float _S602 = r2_16 * _S601;
    float _S603 = r2_16 * _S602;
    float _S604 = dist_coeffs_4[int(7)] * _S597.y + _S599 + dist_coeffs_4[int(6)] * _S597.x + _S600 + _S393 * _S601 + _S392 * _S602 + _S391 * _S603 + dist_coeffs_4[int(3)] * (r2_16 * _S603);
    float _S605 = v_16 * _S604;
    float _S606 = u_16 * _S604;
    float2  _S607 = (make_float2 (radial_7) * _S597 + make_float2 (_S379 * (v_16 * _S597.y) + _S395 * _S600 + 2.0f * (u_16 * _S600) + _S376 * (v_16 * _S597.x) + _S606 + _S606, _S397 * _S599 + 2.0f * (v_16 * _S599) + _S396 * _S597.y + _S394 * _S597.x + _S605 + _S605)) / _S390;
    float2  _S608 = _S386 * - _S607;
    float2  _S609 = _S388 * _S607;
    float _S610 = _S608.x + _S608.y;
    float2  _S611 = _S573.differential_0 + - _S579;
    float _S612 = fx_4 * (_S564.differential_0 + _S568.differential_0 + _S611.x);
    float2  _S613 = make_float2 (_S612, fy_4 * (_S556.differential_0 + _S560.differential_0 + _S611.y)) + make_float2 (dist_coeffs_4[int(8)] * _S612, dist_coeffs_4[int(9)] * _S612);
    float2  _S614 = _S371 * _S613;
    float _S615 = dist_coeffs_4[int(4)] * _S613.y;
    float _S616 = dist_coeffs_4[int(5)] * _S613.x;
    float _S617 = _S614.x + _S614.y;
    float _S618 = r2_15 * _S617;
    float _S619 = r2_15 * _S618;
    float _S620 = dist_coeffs_4[int(7)] * _S613.y + _S615 + dist_coeffs_4[int(6)] * _S613.x + _S616 + _S375 * _S617 + _S374 * _S618 + _S373 * _S619 + dist_coeffs_4[int(3)] * (r2_15 * _S619);
    float _S621 = v_15 * _S620;
    float _S622 = u_15 * _S620;
    float2  _S623 = (make_float2 (radial_6) * _S613 + make_float2 (_S379 * (v_15 * _S613.y) + _S378 * _S616 + 2.0f * (u_15 * _S616) + _S376 * (v_15 * _S613.x) + _S622 + _S622, _S381 * _S615 + 2.0f * (v_15 * _S615) + _S380 * _S613.y + _S377 * _S613.x + _S621 + _S621)) / _S372;
    float2  _S624 = _S368 * - _S623;
    float2  _S625 = _S370 * _S623;
    float _S626 = _S624.x + _S624.y;
    float3  _S627 = _S475.differential_0 + _S553 + make_float3 (_S594.x, _S594.y, _S593.x + _S593.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S628;
    (&_S628)->primal_0 = R_4;
    (&_S628)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S629;
    (&_S629)->primal_0 = vert2_4;
    (&_S629)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S628, &_S629, _S627);
    float3  _S630 = _S474.differential_0 + _S553 + make_float3 (_S609.x, _S609.y, _S610);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S631;
    (&_S631)->primal_0 = R_4;
    (&_S631)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
    (&_S632)->primal_0 = vert1_4;
    (&_S632)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S631, &_S632, _S630);
    float3  _S633 = _S476 + _S477 + _S553 + make_float3 (_S625.x, _S625.y, _S626);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S634;
    (&_S634)->primal_0 = R_4;
    (&_S634)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S635;
    (&_S635)->primal_0 = vert0_4;
    (&_S635)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S634, &_S635, _S633);
    float3  _S636 = v_verts_0[int(2)] + _S629.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S637;
    (&_S637)->primal_0 = _S362;
    (&_S637)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S367;
    (&_S638)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S637, &_S638, _S636);
    float _S639 = - _S638.differential_0.y;
    float _S640 = _S366 * _S638.differential_0.x;
    float _S641 = - (_S357 * _S638.differential_0.x);
    float3  _S642 = v_verts_0[int(1)] + _S632.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S643;
    (&_S643)->primal_0 = _S362;
    (&_S643)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S644;
    (&_S644)->primal_0 = _S365;
    (&_S644)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S643, &_S644, _S642);
    float _S645 = _S357 * _S644.differential_0.x;
    float _S646 = _S364 * _S644.differential_0.x;
    float3  _S647 = v_verts_0[int(0)] + _S635.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S648;
    (&_S648)->primal_0 = _S362;
    (&_S648)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S649;
    (&_S649)->primal_0 = _S363;
    (&_S649)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S648, &_S649, _S647);
    Matrix<float, 3, 3>  _S650 = transpose_1(_S637.differential_0 + _S643.differential_0 + _S648.differential_0);
    float _S651 = 2.0f * - _S650.rows[int(2)].z;
    float _S652 = 2.0f * _S650.rows[int(2)].y;
    float _S653 = 2.0f * _S650.rows[int(2)].x;
    float _S654 = 2.0f * _S650.rows[int(1)].z;
    float _S655 = 2.0f * - _S650.rows[int(1)].y;
    float _S656 = 2.0f * _S650.rows[int(1)].x;
    float _S657 = 2.0f * _S650.rows[int(0)].z;
    float _S658 = 2.0f * _S650.rows[int(0)].y;
    float _S659 = 2.0f * - _S650.rows[int(0)].x;
    float _S660 = - _S656 + _S658;
    float _S661 = _S653 + - _S657;
    float _S662 = - _S652 + _S654;
    float _S663 = _S652 + _S654;
    float _S664 = _S653 + _S657;
    float _S665 = _S656 + _S658;
    float _S666 = _S360.w * (_S655 + _S659);
    float _S667 = _S360.z * (_S651 + _S659);
    float _S668 = _S360.y * (_S651 + _S655);
    float _S669 = _S360.x * _S660 + _S360.z * _S663 + _S360.y * _S664 + _S666 + _S666;
    float _S670 = _S360.x * _S661 + _S360.w * _S663 + _S360.y * _S665 + _S667 + _S667;
    float _S671 = _S360.x * _S662 + _S360.w * _S664 + _S360.z * _S665 + _S668 + _S668;
    float _S672 = _S360.w * _S660 + _S360.z * _S661 + _S360.y * _S662;
    float4  _S673 = make_float4 (0.0f);
    float4  _S674 = _S673;
    *&((&_S674)->w) = _S669;
    *&((&_S674)->z) = _S670;
    *&((&_S674)->y) = _S671;
    *&((&_S674)->x) = _S672;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S675;
    (&_S675)->primal_0 = quat_4;
    (&_S675)->differential_0 = _S673;
    s_bwd_normalize_impl_1(&_S675, _S674);
    float _S676 = _S641 + _S645;
    float _S677 = 0.5f * - _S676;
    float _S678 = _S639 + _S644.differential_0.y;
    DiffPair_float_0 _S679;
    (&_S679)->primal_0 = _S358;
    (&_S679)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S679, _S678);
    float _S680 = _S677 + _S679.differential_0;
    float _S681 = _S640 + _S646 + _S649.differential_0.x;
    DiffPair_float_0 _S682;
    (&_S682)->primal_0 = _S356;
    (&_S682)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S682, _S681);
    float _S683 = _S677 + _S682.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S684;
    (&_S684)->primal_0 = R_4;
    (&_S684)->differential_0 = _S549;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S685;
    (&_S685)->primal_0 = mean_4;
    (&_S685)->differential_0 = _S469;
    s_bwd_prop_mul_0(&_S684, &_S685, _S471.differential_0);
    float3  _S686 = _S551.differential_0 + _S627 + _S630 + _S633 + _S471.differential_0;
    Matrix<float, 3, 3>  _S687 = _S552 + _S628.differential_0 + _S631.differential_0 + _S634.differential_0 + _S684.differential_0;
    FixedArray<float3 , 2>  _S688;
    _S688[int(0)] = _S469;
    _S688[int(1)] = _S469;
    _S688[int(1)] = _S483;
    _S688[int(0)] = _S488;
    FixedArray<float3 , 16>  _S689;
    _S689[int(0)] = _S469;
    _S689[int(1)] = _S469;
    _S689[int(2)] = _S469;
    _S689[int(3)] = _S469;
    _S689[int(4)] = _S469;
    _S689[int(5)] = _S469;
    _S689[int(6)] = _S469;
    _S689[int(7)] = _S469;
    _S689[int(8)] = _S469;
    _S689[int(9)] = _S469;
    _S689[int(10)] = _S469;
    _S689[int(11)] = _S469;
    _S689[int(12)] = _S469;
    _S689[int(13)] = _S469;
    _S689[int(14)] = _S469;
    _S689[int(15)] = _S469;
    _S689[int(15)] = _S490;
    _S689[int(14)] = _S492;
    _S689[int(13)] = _S494;
    _S689[int(12)] = _S496;
    _S689[int(11)] = _S498;
    _S689[int(10)] = _S500;
    _S689[int(9)] = _S502;
    _S689[int(8)] = _S510;
    _S689[int(7)] = _S512;
    _S689[int(6)] = _S514;
    _S689[int(5)] = _S516;
    _S689[int(4)] = _S518;
    _S689[int(3)] = _S529;
    _S689[int(2)] = _S531;
    _S689[int(1)] = _S533;
    _S689[int(0)] = _S546;
    float2  _S690 = make_float2 (0.0f, _S571);
    float3  _S691 = make_float3 (_S683, _S680, _S676);
    *v_mean_0 = _S547 + _S636 + _S642 + _S647 + _S685.differential_0;
    *v_quat_0 = _S675.differential_0;
    *v_scale_0 = _S691;
    *v_hardness_0 = _S690;
    *v_sh_coeffs_0 = _S689;
    *v_ch_coeffs_0 = _S688;
    *v_R_0 = _S687;
    *v_t_0 = _S686;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S692, float _S693)
{
    return (F32_atan2((_S692), (_S693)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S694, DiffPair_float_0 * _S695, float _S696)
{
    _d_atan2_0(_S694, _S695, _S696);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S697 = scale_5.x;
    float _S698 = s_primal_ctx_exp_0(_S697);
    float _S699 = scale_5.y;
    float _S700 = s_primal_ctx_exp_0(_S699);
    float sz_5 = scale_5.z - 0.5f * (_S697 + _S699);
    float4  _S701 = normalize_0(quat_5);
    float _S702 = _S701.y;
    float x2_5 = _S702 * _S702;
    float y2_5 = _S701.z * _S701.z;
    float z2_10 = _S701.w * _S701.w;
    float xy_5 = _S701.y * _S701.z;
    float xz_5 = _S701.y * _S701.w;
    float yz_5 = _S701.z * _S701.w;
    float wx_5 = _S701.x * _S701.y;
    float wy_5 = _S701.x * _S701.z;
    float wz_5 = _S701.x * _S701.w;
    Matrix<float, 3, 3>  _S703 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S704 = make_float3 (_S698, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S703, _S704) + mean_5;
    float _S705 = -0.5f + sz_5;
    float3  _S706 = make_float3 (_S698 * _S705, _S700, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S703, _S706) + mean_5;
    float _S707 = -0.5f - sz_5;
    float3  _S708 = make_float3 (_S698 * _S707, - _S700, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S703, _S708) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S709 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S710 = length_0(_S709);
    float _S711 = vert0_c_5.z;
    float _S712 = s_primal_ctx_atan2_0(_S710, _S711);
    bool _S713 = _S712 < 0.00100000004749745f;
    float k_2;
    float _S714;
    float _S715;
    float _S716;
    if(_S713)
    {
        float _S717 = 1.0f - _S712 * _S712 / 3.0f;
        float _S718 = _S711 * _S711;
        k_2 = _S717 / _S711;
        _S714 = _S718;
        _S715 = _S717;
        _S716 = 0.0f;
    }
    else
    {
        float _S719 = _S710 * _S710;
        k_2 = _S712 / _S710;
        _S714 = 0.0f;
        _S715 = 0.0f;
        _S716 = _S719;
    }
    float2  _S720 = make_float2 (k_2);
    float2  _S721 = _S709 * make_float2 (k_2);
    float u_18 = _S721.x;
    float v_18 = _S721.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S722 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S723 = dist_coeffs_5[int(1)] + r2_18 * _S722;
    float _S724 = dist_coeffs_5[int(0)] + r2_18 * _S723;
    float radial_9 = 1.0f + r2_18 * _S724;
    float _S725 = 2.0f * dist_coeffs_5[int(4)];
    float _S726 = _S725 * u_18;
    float _S727 = 2.0f * u_18;
    float _S728 = 2.0f * dist_coeffs_5[int(5)];
    float _S729 = _S728 * u_18;
    float _S730 = 2.0f * v_18;
    float2  _S731 = _S721 * make_float2 (radial_9) + make_float2 (_S726 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S727 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S729 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S730 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S732 = _S731 + make_float2 (dist_coeffs_5[int(8)] * _S731.x + dist_coeffs_5[int(9)] * _S731.y, 0.0f);
    float _S733 = fx_5 * _S732.x + cx_5;
    float _S734 = fy_5 * _S732.y + cy_5;
    float2  uv0_7 = make_float2 (_S733, _S734);
    float2  _S735 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S736 = length_0(_S735);
    float _S737 = vert1_c_5.z;
    float _S738 = s_primal_ctx_atan2_0(_S736, _S737);
    bool _S739 = _S738 < 0.00100000004749745f;
    float _S740;
    float _S741;
    float _S742;
    if(_S739)
    {
        float _S743 = 1.0f - _S738 * _S738 / 3.0f;
        float _S744 = _S737 * _S737;
        k_2 = _S743 / _S737;
        _S740 = _S744;
        _S741 = _S743;
        _S742 = 0.0f;
    }
    else
    {
        float _S745 = _S736 * _S736;
        k_2 = _S738 / _S736;
        _S740 = 0.0f;
        _S741 = 0.0f;
        _S742 = _S745;
    }
    float2  _S746 = make_float2 (k_2);
    float2  _S747 = _S735 * make_float2 (k_2);
    float u_19 = _S747.x;
    float v_19 = _S747.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S748 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S749 = dist_coeffs_5[int(1)] + r2_19 * _S748;
    float _S750 = dist_coeffs_5[int(0)] + r2_19 * _S749;
    float radial_10 = 1.0f + r2_19 * _S750;
    float _S751 = _S725 * u_19;
    float _S752 = 2.0f * u_19;
    float _S753 = _S728 * u_19;
    float _S754 = 2.0f * v_19;
    float2  _S755 = _S747 * make_float2 (radial_10) + make_float2 (_S751 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S752 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S753 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S754 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S756 = _S755 + make_float2 (dist_coeffs_5[int(8)] * _S755.x + dist_coeffs_5[int(9)] * _S755.y, 0.0f);
    float _S757 = fx_5 * _S756.x + cx_5;
    float _S758 = fy_5 * _S756.y + cy_5;
    float2  uv1_7 = make_float2 (_S757, _S758);
    float2  _S759 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S760 = length_0(_S759);
    float _S761 = vert2_c_5.z;
    float _S762 = s_primal_ctx_atan2_0(_S760, _S761);
    bool _S763 = _S762 < 0.00100000004749745f;
    float _S764;
    float _S765;
    float _S766;
    if(_S763)
    {
        float _S767 = 1.0f - _S762 * _S762 / 3.0f;
        float _S768 = _S761 * _S761;
        k_2 = _S767 / _S761;
        _S764 = _S768;
        _S765 = _S767;
        _S766 = 0.0f;
    }
    else
    {
        float _S769 = _S760 * _S760;
        k_2 = _S762 / _S760;
        _S764 = 0.0f;
        _S765 = 0.0f;
        _S766 = _S769;
    }
    float2  _S770 = make_float2 (k_2);
    float2  _S771 = _S759 * make_float2 (k_2);
    float u_20 = _S771.x;
    float v_20 = _S771.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S772 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S773 = dist_coeffs_5[int(1)] + r2_20 * _S772;
    float _S774 = dist_coeffs_5[int(0)] + r2_20 * _S773;
    float radial_11 = 1.0f + r2_20 * _S774;
    float _S775 = _S725 * u_20;
    float _S776 = 2.0f * u_20;
    float _S777 = _S728 * u_20;
    float _S778 = 2.0f * v_20;
    float2  _S779 = _S771 * make_float2 (radial_11) + make_float2 (_S775 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S776 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S777 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S778 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S780 = _S779 + make_float2 (dist_coeffs_5[int(8)] * _S779.x + dist_coeffs_5[int(9)] * _S779.y, 0.0f);
    float _S781 = fx_5 * _S780.x + cx_5;
    float _S782 = fy_5 * _S780.y + cy_5;
    float2  uv2_7 = make_float2 (_S781, _S782);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S783 = e0_5.x;
    float _S784 = e1_5.y;
    float _S785 = e0_5.y;
    float _S786 = e1_5.x;
    float _S787 = _S783 * _S784 - _S785 * _S786;
    float _S788 = 1.0f - hardness_5.y;
    float _S789 = -1.0f / _S788;
    float _S790 = _S788 * _S788;
    float _S791 = (F32_max((_S733), (_S757)));
    float _S792 = (F32_min((_S733), (_S757)));
    float _S793 = (F32_max((_S734), (_S758)));
    float _S794 = (F32_min((_S734), (_S758)));
    float3  _S795 = (vert0_c_5 + vert1_c_5 + vert2_c_5) / make_float3 (3.0f);
    float x_22 = _S795.x;
    float y_11 = _S795.y;
    float z_7 = _S795.z;
    float _S796 = z_7 * z_7;
    float _S797 = _S796 * z_7;
    float _S798 = x_22 * x_22 + y_11 * y_11;
    float _S799 = 0.001953125f * _S798;
    Matrix<float, 3, 3>  _S800 = transpose_1(R_5);
    float3  _S801 = mean_5 - - s_primal_ctx_mul_0(_S800, t_5);
    float _S802 = _S801.x;
    float _S803 = _S801.y;
    float _S804 = _S801.z;
    float _S805 = _S802 * _S802 + _S803 * _S803 + _S804 * _S804;
    float _S806 = s_primal_ctx_sqrt_0(_S805);
    float x_23 = _S802 / _S806;
    float3  _S807 = make_float3 (x_23);
    float _S808 = _S806 * _S806;
    float y_12 = _S803 / _S806;
    float z_8 = _S804 / _S806;
    float3  _S809 = make_float3 (z_8);
    float _S810 = - y_12;
    float3  _S811 = make_float3 (_S810);
    float z2_11 = z_8 * z_8;
    float fTmp0B_5 = -1.09254848957061768f * z_8;
    float fC1_5 = x_23 * x_23 - y_12 * y_12;
    float _S812 = 2.0f * x_23;
    float fS1_5 = _S812 * y_12;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S813 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_23;
    float3  _S814 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_12;
    float3  _S815 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S816 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S817 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_8;
    float _S818 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_8 * _S818;
    float3  _S819 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_23;
    float3  _S820 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_12;
    float3  _S821 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S822 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S823 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_23 * fC1_5 - y_12 * fS1_5);
    float3  _S824 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_23 * fS1_5 + y_12 * fC1_5);
    float3  _S825 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S810) * sh_coeffs_5[int(1)] + make_float3 (z_8) * sh_coeffs_5[int(2)] - make_float3 (x_23) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]);
    float3  _S826 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S827 = make_float3 (0.0f);
    float3  _S828 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S829 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S830 = make_float3 (_S829);
    float3  _S831 = ch_coeffs_5[int(1)] * make_float3 (_S829);
    float3  _S832 = _S828 + _S831 + make_float3 (0.5f);
    float3  _S833 = _S828 - _S831 + make_float3 (0.5f);
    float3  _S834 = vert1_c_5 - vert0_c_5;
    float3  _S835 = vert2_c_5 - vert0_c_5;
    float3  _S836 = s_primal_ctx_cross_0(_S834, _S835);
    float3  _S837 = normalize_1(_S836);
    float3  _S838 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S837, mean_c_5)))))) * v_normal_1;
    float3  _S839 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S840;
    (&_S840)->primal_0 = _S837;
    (&_S840)->differential_0 = _S839;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S841;
    (&_S841)->primal_0 = mean_c_5;
    (&_S841)->differential_0 = _S839;
    s_bwd_prop_dot_0(&_S840, &_S841, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S842 = _S841;
    float3  _S843 = _S838 + _S840.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S844;
    (&_S844)->primal_0 = _S836;
    (&_S844)->differential_0 = _S839;
    s_bwd_normalize_impl_0(&_S844, _S843);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S845;
    (&_S845)->primal_0 = _S834;
    (&_S845)->differential_0 = _S839;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S846;
    (&_S846)->primal_0 = _S835;
    (&_S846)->differential_0 = _S839;
    s_bwd_prop_cross_0(&_S845, &_S846, _S844.differential_0);
    float3  _S847 = - _S846.differential_0;
    float3  _S848 = - _S845.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S849;
    (&_S849)->primal_0 = _S833;
    (&_S849)->differential_0 = _S839;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S850;
    (&_S850)->primal_0 = _S827;
    (&_S850)->differential_0 = _S839;
    s_bwd_prop_max_0(&_S849, &_S850, v_rgbs_1[int(2)]);
    float3  _S851 = - _S849.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S852;
    (&_S852)->primal_0 = _S832;
    (&_S852)->differential_0 = _S839;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S853;
    (&_S853)->primal_0 = _S827;
    (&_S853)->differential_0 = _S839;
    s_bwd_prop_max_0(&_S852, &_S853, v_rgbs_1[int(1)]);
    float3  _S854 = _S830 * (_S851 + _S852.differential_0);
    float3  _S855 = _S849.differential_0 + _S852.differential_0;
    float3  _S856 = make_float3 (0.5f) * - _S855;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S857;
    (&_S857)->primal_0 = _S826;
    (&_S857)->differential_0 = _S839;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S858;
    (&_S858)->primal_0 = _S827;
    (&_S858)->differential_0 = _S839;
    s_bwd_prop_max_0(&_S857, &_S858, v_rgbs_1[int(0)]);
    float3  _S859 = _S856 + _S857.differential_0;
    float3  _S860 = _S855 + _S857.differential_0;
    float3  _S861 = _S824 * _S860;
    float3  _S862 = sh_coeffs_5[int(15)] * _S860;
    float3  _S863 = _S822 * _S860;
    float3  _S864 = sh_coeffs_5[int(14)] * _S860;
    float3  _S865 = _S820 * _S860;
    float3  _S866 = sh_coeffs_5[int(13)] * _S860;
    float3  _S867 = _S819 * _S860;
    float3  _S868 = sh_coeffs_5[int(12)] * _S860;
    float3  _S869 = _S821 * _S860;
    float3  _S870 = sh_coeffs_5[int(11)] * _S860;
    float3  _S871 = _S823 * _S860;
    float3  _S872 = sh_coeffs_5[int(10)] * _S860;
    float3  _S873 = _S825 * _S860;
    float3  _S874 = sh_coeffs_5[int(9)] * _S860;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S874.x + _S874.y + _S874.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S862.x + _S862.y + _S862.z);
    float _S875 = _S872.x + _S872.y + _S872.z;
    float _S876 = _S864.x + _S864.y + _S864.z;
    float _S877 = _S870.x + _S870.y + _S870.z;
    float _S878 = _S866.x + _S866.y + _S866.z;
    float _S879 = _S868.x + _S868.y + _S868.z;
    float _S880 = - s_diff_fC2_T_1;
    float3  _S881 = _S816 * _S860;
    float3  _S882 = sh_coeffs_5[int(8)] * _S860;
    float3  _S883 = _S814 * _S860;
    float3  _S884 = sh_coeffs_5[int(7)] * _S860;
    float3  _S885 = _S813 * _S860;
    float3  _S886 = sh_coeffs_5[int(6)] * _S860;
    float3  _S887 = _S815 * _S860;
    float3  _S888 = sh_coeffs_5[int(5)] * _S860;
    float3  _S889 = _S817 * _S860;
    float3  _S890 = sh_coeffs_5[int(4)] * _S860;
    float _S891 = _S888.x + _S888.y + _S888.z;
    float _S892 = _S884.x + _S884.y + _S884.z;
    float _S893 = fTmp1B_5 * _S875 + x_23 * s_diff_fS2_T_1 + y_12 * _S880 + 0.54627424478530884f * (_S890.x + _S890.y + _S890.z);
    float _S894 = fTmp1B_5 * _S876 + y_12 * s_diff_fS2_T_1 + x_23 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S882.x + _S882.y + _S882.z);
    float _S895 = y_12 * - _S894;
    float _S896 = x_23 * _S894;
    float _S897 = z_8 * (1.86588168144226074f * (z_8 * _S879) + -2.28522896766662598f * (y_12 * _S877 + x_23 * _S878) + 0.94617468118667603f * (_S886.x + _S886.y + _S886.z));
    float3  _S898 = make_float3 (0.48860251903533936f) * _S860;
    float3  _S899 = - _S898;
    float3  _S900 = _S807 * _S899;
    float3  _S901 = sh_coeffs_5[int(3)] * _S899;
    float3  _S902 = _S809 * _S898;
    float3  _S903 = sh_coeffs_5[int(2)] * _S898;
    float3  _S904 = _S811 * _S898;
    float3  _S905 = sh_coeffs_5[int(1)] * _S898;
    float _S906 = (_S818 * _S879 + 1.44530570507049561f * (fS1_5 * _S875 + fC1_5 * _S876) + -1.09254848957061768f * (y_12 * _S891 + x_23 * _S892) + _S897 + _S897 + _S903.x + _S903.y + _S903.z) / _S808;
    float _S907 = _S806 * _S906;
    float _S908 = (fTmp0C_5 * _S877 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S880 + fTmp0B_5 * _S891 + _S812 * _S893 + _S895 + _S895 + - (_S905.x + _S905.y + _S905.z)) / _S808;
    float _S909 = _S806 * _S908;
    float _S910 = (fTmp0C_5 * _S878 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S892 + 2.0f * (y_12 * _S893) + _S896 + _S896 + _S901.x + _S901.y + _S901.z) / _S808;
    float _S911 = _S806 * _S910;
    float _S912 = _S804 * - _S906 + _S803 * - _S908 + _S802 * - _S910;
    DiffPair_float_0 _S913;
    (&_S913)->primal_0 = _S805;
    (&_S913)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S913, _S912);
    float _S914 = _S804 * _S913.differential_0;
    float _S915 = _S803 * _S913.differential_0;
    float _S916 = _S802 * _S913.differential_0;
    float3  _S917 = make_float3 (0.282094806432724f) * _S860;
    float3  _S918 = make_float3 (_S911 + _S916 + _S916, _S909 + _S915 + _S915, _S907 + _S914 + _S914);
    float3  _S919 = - - _S918;
    Matrix<float, 3, 3>  _S920 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S921;
    (&_S921)->primal_0 = _S800;
    (&_S921)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S922;
    (&_S922)->primal_0 = t_5;
    (&_S922)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S921, &_S922, _S919);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S923 = _S922;
    Matrix<float, 3, 3>  _S924 = transpose_1(_S921.differential_0);
    float _S925 = _S799 * v_depth_1 + 0.001953125f * (_S798 * v_depth_1);
    float _S926 = y_11 * _S925;
    float _S927 = x_22 * _S925;
    float _S928 = z_7 * v_depth_1;
    float _S929 = z_7 * (z_7 * _S928);
    float3  _S930 = make_float3 (0.3333333432674408f) * make_float3 (_S927 + _S927, _S926 + _S926, _S797 * v_depth_1 + _S796 * _S928 + _S929 + _S929);
    DiffPair_float_0 _S931;
    (&_S931)->primal_0 = _S794;
    (&_S931)->differential_0 = 0.0f;
    DiffPair_float_0 _S932;
    (&_S932)->primal_0 = _S782;
    (&_S932)->differential_0 = 0.0f;
    _d_min_0(&_S931, &_S932, 0.0f);
    DiffPair_float_0 _S933;
    (&_S933)->primal_0 = _S734;
    (&_S933)->differential_0 = 0.0f;
    DiffPair_float_0 _S934;
    (&_S934)->primal_0 = _S758;
    (&_S934)->differential_0 = 0.0f;
    _d_min_0(&_S933, &_S934, _S931.differential_0);
    DiffPair_float_0 _S935;
    (&_S935)->primal_0 = _S793;
    (&_S935)->differential_0 = 0.0f;
    DiffPair_float_0 _S936;
    (&_S936)->primal_0 = _S782;
    (&_S936)->differential_0 = 0.0f;
    _d_max_0(&_S935, &_S936, 0.0f);
    DiffPair_float_0 _S937;
    (&_S937)->primal_0 = _S734;
    (&_S937)->differential_0 = 0.0f;
    DiffPair_float_0 _S938;
    (&_S938)->primal_0 = _S758;
    (&_S938)->differential_0 = 0.0f;
    _d_max_0(&_S937, &_S938, _S935.differential_0);
    DiffPair_float_0 _S939;
    (&_S939)->primal_0 = _S792;
    (&_S939)->differential_0 = 0.0f;
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S781;
    (&_S940)->differential_0 = 0.0f;
    _d_min_0(&_S939, &_S940, 0.0f);
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = _S733;
    (&_S941)->differential_0 = 0.0f;
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = _S757;
    (&_S942)->differential_0 = 0.0f;
    _d_min_0(&_S941, &_S942, _S939.differential_0);
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S791;
    (&_S943)->differential_0 = 0.0f;
    DiffPair_float_0 _S944;
    (&_S944)->primal_0 = _S781;
    (&_S944)->differential_0 = 0.0f;
    _d_max_0(&_S943, &_S944, 0.0f);
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S733;
    (&_S945)->differential_0 = 0.0f;
    DiffPair_float_0 _S946;
    (&_S946)->primal_0 = _S757;
    (&_S946)->differential_0 = 0.0f;
    _d_max_0(&_S945, &_S946, _S943.differential_0);
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S789;
    (&_S947)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S947, -0.0f);
    float _S948 = - (-1.0f * - (_S947.differential_0 / _S790));
    float2  _S949 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S950;
    (&_S950)->primal_0 = e2_1;
    (&_S950)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S950, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S951;
    (&_S951)->primal_0 = e1_5;
    (&_S951)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S951, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S952;
    (&_S952)->primal_0 = e0_5;
    (&_S952)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S952, 0.0f);
    DiffPair_float_0 _S953;
    (&_S953)->primal_0 = _S787;
    (&_S953)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S953, -0.0f);
    float _S954 = - _S953.differential_0;
    float2  _S955 = _S951.differential_0 + make_float2 (_S785 * _S954, _S783 * _S953.differential_0);
    float2  _S956 = - _S955;
    float2  _S957 = _S952.differential_0 + make_float2 (_S784 * _S953.differential_0, _S786 * _S954);
    float2  _S958 = - _S957;
    float2  _S959 = - _S950.differential_0 + _S955;
    float _S960 = fx_5 * (_S940.differential_0 + _S944.differential_0 + _S959.x);
    float2  _S961 = make_float2 (_S960, fy_5 * (_S932.differential_0 + _S936.differential_0 + _S959.y)) + make_float2 (dist_coeffs_5[int(8)] * _S960, dist_coeffs_5[int(9)] * _S960);
    float2  _S962 = _S771 * _S961;
    float _S963 = dist_coeffs_5[int(4)] * _S961.y;
    float _S964 = dist_coeffs_5[int(5)] * _S961.x;
    float _S965 = _S962.x + _S962.y;
    float _S966 = r2_20 * _S965;
    float _S967 = r2_20 * _S966;
    float _S968 = dist_coeffs_5[int(7)] * _S961.y + _S963 + dist_coeffs_5[int(6)] * _S961.x + _S964 + _S774 * _S965 + _S773 * _S966 + _S772 * _S967 + dist_coeffs_5[int(3)] * (r2_20 * _S967);
    float _S969 = v_20 * _S968;
    float _S970 = u_20 * _S968;
    float2  _S971 = make_float2 (radial_11) * _S961 + make_float2 (_S728 * (v_20 * _S961.y) + _S776 * _S964 + 2.0f * (u_20 * _S964) + _S725 * (v_20 * _S961.x) + _S970 + _S970, _S778 * _S963 + 2.0f * (v_20 * _S963) + _S777 * _S961.y + _S775 * _S961.x + _S969 + _S969);
    float3  _S972 = _S845.differential_0 + _S930;
    float3  _S973 = _S847 + _S848 + _S930;
    float3  _S974 = _S846.differential_0 + _S930;
    FixedArray<float3 , 2>  _S975;
    _S975[int(0)] = _S839;
    _S975[int(1)] = _S839;
    _S975[int(1)] = _S854;
    _S975[int(0)] = _S859;
    float3  _S976 = _S975[int(0)];
    float3  _S977 = _S975[int(1)];
    FixedArray<float3 , 16>  _S978;
    _S978[int(0)] = _S839;
    _S978[int(1)] = _S839;
    _S978[int(2)] = _S839;
    _S978[int(3)] = _S839;
    _S978[int(4)] = _S839;
    _S978[int(5)] = _S839;
    _S978[int(6)] = _S839;
    _S978[int(7)] = _S839;
    _S978[int(8)] = _S839;
    _S978[int(9)] = _S839;
    _S978[int(10)] = _S839;
    _S978[int(11)] = _S839;
    _S978[int(12)] = _S839;
    _S978[int(13)] = _S839;
    _S978[int(14)] = _S839;
    _S978[int(15)] = _S839;
    _S978[int(7)] = _S883;
    _S978[int(0)] = _S917;
    _S978[int(1)] = _S904;
    _S978[int(2)] = _S902;
    _S978[int(3)] = _S900;
    _S978[int(4)] = _S889;
    _S978[int(5)] = _S887;
    _S978[int(6)] = _S885;
    _S978[int(15)] = _S861;
    _S978[int(8)] = _S881;
    _S978[int(9)] = _S873;
    _S978[int(10)] = _S871;
    _S978[int(11)] = _S869;
    _S978[int(12)] = _S867;
    _S978[int(13)] = _S865;
    _S978[int(14)] = _S863;
    float3  _S979 = _S978[int(0)];
    float3  _S980 = _S978[int(1)];
    float3  _S981 = _S978[int(2)];
    float3  _S982 = _S978[int(3)];
    float3  _S983 = _S978[int(4)];
    float3  _S984 = _S978[int(5)];
    float3  _S985 = _S978[int(6)];
    float3  _S986 = _S978[int(7)];
    float3  _S987 = _S978[int(8)];
    float3  _S988 = _S978[int(9)];
    float3  _S989 = _S978[int(10)];
    float3  _S990 = _S978[int(11)];
    float3  _S991 = _S978[int(12)];
    float3  _S992 = _S978[int(13)];
    float3  _S993 = _S978[int(14)];
    float3  _S994 = _S978[int(15)];
    float _S995 = _S934.differential_0 + _S938.differential_0;
    float _S996 = _S933.differential_0 + _S937.differential_0;
    float _S997 = _S941.differential_0 + _S945.differential_0;
    float _S998 = _S942.differential_0 + _S946.differential_0;
    float2  _S999 = _S950.differential_0 + _S958;
    float2  _S1000 = _S956 + _S957;
    float2  _S1001 = make_float2 (0.0f, _S948);
    float2  _S1002 = _S759 * _S971;
    float2  _S1003 = _S770 * _S971;
    float _S1004 = _S1002.x + _S1002.y;
    if(_S763)
    {
        float _S1005 = _S1004 / _S764;
        float _S1006 = _S765 * - _S1005;
        float _S1007 = _S762 * (0.3333333432674408f * - (_S761 * _S1005));
        k_2 = _S1007 + _S1007;
        _S764 = _S1006;
        _S765 = 0.0f;
    }
    else
    {
        float _S1008 = _S1004 / _S766;
        float _S1009 = _S762 * - _S1008;
        k_2 = _S760 * _S1008;
        _S764 = 0.0f;
        _S765 = _S1009;
    }
    DiffPair_float_0 _S1010;
    (&_S1010)->primal_0 = _S760;
    (&_S1010)->differential_0 = 0.0f;
    DiffPair_float_0 _S1011;
    (&_S1011)->primal_0 = _S761;
    (&_S1011)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1010, &_S1011, k_2);
    float _S1012 = _S1011.differential_0 + _S764;
    float _S1013 = _S1010.differential_0 + _S765;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1014;
    (&_S1014)->primal_0 = _S759;
    (&_S1014)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S1014, _S1013);
    float2  _S1015 = _S1014.differential_0 + _S1003;
    float _S1016 = fx_5 * (_S1000.x + _S998);
    float2  _S1017 = make_float2 (_S1016, fy_5 * (_S1000.y + _S995)) + make_float2 (dist_coeffs_5[int(8)] * _S1016, dist_coeffs_5[int(9)] * _S1016);
    float2  _S1018 = _S747 * _S1017;
    float _S1019 = dist_coeffs_5[int(4)] * _S1017.y;
    float _S1020 = dist_coeffs_5[int(5)] * _S1017.x;
    float _S1021 = _S1018.x + _S1018.y;
    float _S1022 = r2_19 * _S1021;
    float _S1023 = r2_19 * _S1022;
    float _S1024 = dist_coeffs_5[int(7)] * _S1017.y + _S1019 + dist_coeffs_5[int(6)] * _S1017.x + _S1020 + _S750 * _S1021 + _S749 * _S1022 + _S748 * _S1023 + dist_coeffs_5[int(3)] * (r2_19 * _S1023);
    float _S1025 = v_19 * _S1024;
    float _S1026 = u_19 * _S1024;
    float2  _S1027 = make_float2 (radial_10) * _S1017 + make_float2 (_S728 * (v_19 * _S1017.y) + _S752 * _S1020 + 2.0f * (u_19 * _S1020) + _S725 * (v_19 * _S1017.x) + _S1026 + _S1026, _S754 * _S1019 + 2.0f * (v_19 * _S1019) + _S753 * _S1017.y + _S751 * _S1017.x + _S1025 + _S1025);
    float3  _S1028 = _S974 + make_float3 (_S1015.x, _S1015.y, _S1012);
    float2  _S1029 = _S735 * _S1027;
    float2  _S1030 = _S746 * _S1027;
    float _S1031 = _S1029.x + _S1029.y;
    if(_S739)
    {
        float _S1032 = _S1031 / _S740;
        float _S1033 = _S741 * - _S1032;
        float _S1034 = _S738 * (0.3333333432674408f * - (_S737 * _S1032));
        k_2 = _S1034 + _S1034;
        _S740 = _S1033;
        _S741 = 0.0f;
    }
    else
    {
        float _S1035 = _S1031 / _S742;
        float _S1036 = _S738 * - _S1035;
        k_2 = _S736 * _S1035;
        _S740 = 0.0f;
        _S741 = _S1036;
    }
    DiffPair_float_0 _S1037;
    (&_S1037)->primal_0 = _S736;
    (&_S1037)->differential_0 = 0.0f;
    DiffPair_float_0 _S1038;
    (&_S1038)->primal_0 = _S737;
    (&_S1038)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1037, &_S1038, k_2);
    float _S1039 = _S1038.differential_0 + _S740;
    float _S1040 = _S1037.differential_0 + _S741;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1041;
    (&_S1041)->primal_0 = _S735;
    (&_S1041)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S1041, _S1040);
    float2  _S1042 = _S1041.differential_0 + _S1030;
    float _S1043 = fx_5 * (_S999.x + _S997);
    float2  _S1044 = make_float2 (_S1043, fy_5 * (_S999.y + _S996)) + make_float2 (dist_coeffs_5[int(8)] * _S1043, dist_coeffs_5[int(9)] * _S1043);
    float2  _S1045 = _S721 * _S1044;
    float _S1046 = dist_coeffs_5[int(4)] * _S1044.y;
    float _S1047 = dist_coeffs_5[int(5)] * _S1044.x;
    float _S1048 = _S1045.x + _S1045.y;
    float _S1049 = r2_18 * _S1048;
    float _S1050 = r2_18 * _S1049;
    float _S1051 = dist_coeffs_5[int(7)] * _S1044.y + _S1046 + dist_coeffs_5[int(6)] * _S1044.x + _S1047 + _S724 * _S1048 + _S723 * _S1049 + _S722 * _S1050 + dist_coeffs_5[int(3)] * (r2_18 * _S1050);
    float _S1052 = v_18 * _S1051;
    float _S1053 = u_18 * _S1051;
    float2  _S1054 = make_float2 (radial_9) * _S1044 + make_float2 (_S728 * (v_18 * _S1044.y) + _S727 * _S1047 + 2.0f * (u_18 * _S1047) + _S725 * (v_18 * _S1044.x) + _S1053 + _S1053, _S730 * _S1046 + 2.0f * (v_18 * _S1046) + _S729 * _S1044.y + _S726 * _S1044.x + _S1052 + _S1052);
    float3  _S1055 = _S972 + make_float3 (_S1042.x, _S1042.y, _S1039);
    float2  _S1056 = _S709 * _S1054;
    float2  _S1057 = _S720 * _S1054;
    float _S1058 = _S1056.x + _S1056.y;
    if(_S713)
    {
        float _S1059 = _S1058 / _S714;
        float _S1060 = _S715 * - _S1059;
        float _S1061 = _S712 * (0.3333333432674408f * - (_S711 * _S1059));
        k_2 = _S1061 + _S1061;
        _S714 = _S1060;
        _S715 = 0.0f;
    }
    else
    {
        float _S1062 = _S1058 / _S716;
        float _S1063 = _S712 * - _S1062;
        k_2 = _S710 * _S1062;
        _S714 = 0.0f;
        _S715 = _S1063;
    }
    DiffPair_float_0 _S1064;
    (&_S1064)->primal_0 = _S710;
    (&_S1064)->differential_0 = 0.0f;
    DiffPair_float_0 _S1065;
    (&_S1065)->primal_0 = _S711;
    (&_S1065)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1064, &_S1065, k_2);
    float _S1066 = _S1065.differential_0 + _S714;
    float _S1067 = _S1064.differential_0 + _S715;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1068;
    (&_S1068)->primal_0 = _S709;
    (&_S1068)->differential_0 = _S949;
    s_bwd_length_impl_1(&_S1068, _S1067);
    float2  _S1069 = _S1068.differential_0 + _S1057;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1070;
    (&_S1070)->primal_0 = vert2_c_5;
    (&_S1070)->differential_0 = _S839;
    s_bwd_length_impl_0(&_S1070, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1071;
    (&_S1071)->primal_0 = vert1_c_5;
    (&_S1071)->differential_0 = _S839;
    s_bwd_length_impl_0(&_S1071, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1072;
    (&_S1072)->primal_0 = vert0_c_5;
    (&_S1072)->differential_0 = _S839;
    s_bwd_length_impl_0(&_S1072, 0.0f);
    float3  _S1073 = _S1070.differential_0 + _S1028;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1074;
    (&_S1074)->primal_0 = R_5;
    (&_S1074)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1075;
    (&_S1075)->primal_0 = vert2_5;
    (&_S1075)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1074, &_S1075, _S1073);
    float3  _S1076 = _S1071.differential_0 + _S1055;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1077;
    (&_S1077)->primal_0 = R_5;
    (&_S1077)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1078;
    (&_S1078)->primal_0 = vert1_5;
    (&_S1078)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1077, &_S1078, _S1076);
    float3  _S1079 = _S1072.differential_0 + _S973 + make_float3 (_S1069.x, _S1069.y, _S1066);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1080;
    (&_S1080)->primal_0 = R_5;
    (&_S1080)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1081;
    (&_S1081)->primal_0 = vert0_5;
    (&_S1081)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1080, &_S1081, _S1079);
    float3  _S1082 = _S1075.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1083;
    (&_S1083)->primal_0 = _S703;
    (&_S1083)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1084;
    (&_S1084)->primal_0 = _S708;
    (&_S1084)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1083, &_S1084, _S1082);
    float _S1085 = - _S1084.differential_0.y;
    float _S1086 = _S707 * _S1084.differential_0.x;
    float _S1087 = - (_S698 * _S1084.differential_0.x);
    float3  _S1088 = _S1078.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1089;
    (&_S1089)->primal_0 = _S703;
    (&_S1089)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1090;
    (&_S1090)->primal_0 = _S706;
    (&_S1090)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1089, &_S1090, _S1088);
    float _S1091 = _S698 * _S1090.differential_0.x;
    float _S1092 = _S705 * _S1090.differential_0.x;
    float3  _S1093 = _S1081.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1094;
    (&_S1094)->primal_0 = _S703;
    (&_S1094)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1095;
    (&_S1095)->primal_0 = _S704;
    (&_S1095)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1094, &_S1095, _S1093);
    Matrix<float, 3, 3>  _S1096 = transpose_1(_S1083.differential_0 + _S1089.differential_0 + _S1094.differential_0);
    float _S1097 = 2.0f * - _S1096.rows[int(2)].z;
    float _S1098 = 2.0f * _S1096.rows[int(2)].y;
    float _S1099 = 2.0f * _S1096.rows[int(2)].x;
    float _S1100 = 2.0f * _S1096.rows[int(1)].z;
    float _S1101 = 2.0f * - _S1096.rows[int(1)].y;
    float _S1102 = 2.0f * _S1096.rows[int(1)].x;
    float _S1103 = 2.0f * _S1096.rows[int(0)].z;
    float _S1104 = 2.0f * _S1096.rows[int(0)].y;
    float _S1105 = 2.0f * - _S1096.rows[int(0)].x;
    float _S1106 = - _S1102 + _S1104;
    float _S1107 = _S1099 + - _S1103;
    float _S1108 = - _S1098 + _S1100;
    float _S1109 = _S1098 + _S1100;
    float _S1110 = _S1099 + _S1103;
    float _S1111 = _S1102 + _S1104;
    float _S1112 = _S701.w * (_S1101 + _S1105);
    float _S1113 = _S701.z * (_S1097 + _S1105);
    float _S1114 = _S701.y * (_S1097 + _S1101);
    float _S1115 = _S701.x * _S1106 + _S701.z * _S1109 + _S701.y * _S1110 + _S1112 + _S1112;
    float _S1116 = _S701.x * _S1107 + _S701.w * _S1109 + _S701.y * _S1111 + _S1113 + _S1113;
    float _S1117 = _S701.x * _S1108 + _S701.w * _S1110 + _S701.z * _S1111 + _S1114 + _S1114;
    float _S1118 = _S701.w * _S1106 + _S701.z * _S1107 + _S701.y * _S1108;
    float4  _S1119 = make_float4 (0.0f);
    float4  _S1120 = _S1119;
    *&((&_S1120)->w) = _S1115;
    *&((&_S1120)->z) = _S1116;
    *&((&_S1120)->y) = _S1117;
    *&((&_S1120)->x) = _S1118;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1121;
    (&_S1121)->primal_0 = quat_5;
    (&_S1121)->differential_0 = _S1119;
    s_bwd_normalize_impl_1(&_S1121, _S1120);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1122 = _S1121;
    float _S1123 = _S1087 + _S1091;
    float _S1124 = 0.5f * - _S1123;
    float _S1125 = _S1085 + _S1090.differential_0.y;
    DiffPair_float_0 _S1126;
    (&_S1126)->primal_0 = _S699;
    (&_S1126)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1126, _S1125);
    float _S1127 = _S1124 + _S1126.differential_0;
    float _S1128 = _S1086 + _S1092 + _S1095.differential_0.x;
    DiffPair_float_0 _S1129;
    (&_S1129)->primal_0 = _S697;
    (&_S1129)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1129, _S1128);
    float _S1130 = _S1124 + _S1129.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1131;
    (&_S1131)->primal_0 = mean_c_5;
    (&_S1131)->differential_0 = _S839;
    s_bwd_length_impl_0(&_S1131, 0.0f);
    float3  _S1132 = _S1131.differential_0 + _S842.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1133;
    (&_S1133)->primal_0 = R_5;
    (&_S1133)->differential_0 = _S920;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1134;
    (&_S1134)->primal_0 = mean_5;
    (&_S1134)->differential_0 = _S839;
    s_bwd_prop_mul_0(&_S1133, &_S1134, _S1132);
    float3  _S1135 = _S1073 + _S1076 + _S1079 + _S1132 + _S923.differential_0;
    Matrix<float, 3, 3>  _S1136 = _S1074.differential_0 + _S1077.differential_0 + _S1080.differential_0 + _S1133.differential_0 + _S924;
    float3  _S1137 = make_float3 (_S1130, _S1127, _S1123);
    float3  _S1138 = _S1082 + _S1088 + _S1093 + _S1134.differential_0 + _S918;
    *v_mean_1 = _S1138;
    *v_quat_1 = _S1122.differential_0;
    *v_scale_1 = _S1137;
    *v_hardness_1 = _S1001;
    (*v_sh_coeffs_1)[int(0)] = _S979;
    (*v_sh_coeffs_1)[int(1)] = _S980;
    (*v_sh_coeffs_1)[int(2)] = _S981;
    (*v_sh_coeffs_1)[int(3)] = _S982;
    (*v_sh_coeffs_1)[int(4)] = _S983;
    (*v_sh_coeffs_1)[int(5)] = _S984;
    (*v_sh_coeffs_1)[int(6)] = _S985;
    (*v_sh_coeffs_1)[int(7)] = _S986;
    (*v_sh_coeffs_1)[int(8)] = _S987;
    (*v_sh_coeffs_1)[int(9)] = _S988;
    (*v_sh_coeffs_1)[int(10)] = _S989;
    (*v_sh_coeffs_1)[int(11)] = _S990;
    (*v_sh_coeffs_1)[int(12)] = _S991;
    (*v_sh_coeffs_1)[int(13)] = _S992;
    (*v_sh_coeffs_1)[int(14)] = _S993;
    (*v_sh_coeffs_1)[int(15)] = _S994;
    (*v_ch_coeffs_1)[int(0)] = _S976;
    (*v_ch_coeffs_1)[int(1)] = _S977;
    *v_R_1 = _S1136;
    *v_t_1 = _S1135;
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
    bool _S1139;
    if((*u_21) >= 0.0f)
    {
        _S1139 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1139 = false;
    }
    if(_S1139)
    {
        _S1139 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1139 = false;
    }
    if(_S1139)
    {
        _S1139 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1139 = false;
    }
    return _S1139;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_14, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1140 = *dpx_14;
    bool _S1141;
    if(((*dpx_14).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1141 = ((*dpx_14).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1141 = false;
    }
    float _S1142;
    if(_S1141)
    {
        _S1142 = dOut_11;
    }
    else
    {
        _S1142 = 0.0f;
    }
    dpx_14->primal_0 = _S1140.primal_0;
    dpx_14->differential_0 = _S1142;
    DiffPair_float_0 _S1143 = *dpMin_0;
    if((_S1140.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1142 = dOut_11;
    }
    else
    {
        _S1142 = 0.0f;
    }
    dpMin_0->primal_0 = _S1143.primal_0;
    dpMin_0->differential_0 = _S1142;
    DiffPair_float_0 _S1144 = *dpMax_0;
    if(((*dpx_14).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1142 = dOut_11;
    }
    else
    {
        _S1142 = 0.0f;
    }
    dpMax_0->primal_0 = _S1144.primal_0;
    dpMax_0->differential_0 = _S1142;
    return;
}

inline __device__ float clamp_0(float x_24, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_24), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_15, DiffPair_float_0 * dpy_5, float dOut_12)
{
    if(((*dpx_15).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_15->primal_0 = (*dpx_15).primal_0;
        dpx_15->differential_0 = 0.0f;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_15).primal_0), ((*dpy_5).primal_0)));
        DiffPair_float_0 _S1145 = *dpx_15;
        float _S1146 = val_0 * (*dpy_5).primal_0 / (*dpx_15).primal_0 * dOut_12;
        dpx_15->primal_0 = (*dpx_15).primal_0;
        dpx_15->differential_0 = _S1146;
        float _S1147 = val_0 * (F32_log((_S1145.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1147;
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
    bool _S1148;
    if(u_22 >= 0.0f)
    {
        _S1148 = v_22 >= 0.0f;
    }
    else
    {
        _S1148 = false;
    }
    if(_S1148)
    {
        _S1148 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1148 = false;
    }
    if(_S1148)
    {
        _S1148 = t_7 >= 0.0f;
    }
    else
    {
        _S1148 = false;
    }
    if(!_S1148)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1149;
    if(opac_0 < 0.0f)
    {
        _S1149 = 0.0f;
    }
    else
    {
        _S1149 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1149;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1150, float _S1151, float _S1152)
{
    return clamp_0(_S1150, _S1151, _S1152);
}

inline __device__ float s_primal_ctx_pow_0(float _S1153, float _S1154)
{
    return (F32_pow((_S1153), (_S1154)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1155, DiffPair_float_0 * _S1156, float _S1157)
{
    _d_pow_0(_S1155, _S1156, _S1157);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1158, DiffPair_float_0 * _S1159, DiffPair_float_0 * _S1160, float _S1161)
{
    _d_clamp_0(_S1158, _S1159, _S1160, _S1161);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_5)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1162 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1163 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1164 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1165 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1166 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1164);
    float d_2 = 1.0f / _S1166;
    float _S1167 = _S1166 * _S1166;
    float3  _S1168 = - _S1165;
    float _S1169 = s_primal_ctx_dot_0(_S1168, v2v0_2);
    float u_23 = d_2 * _S1169;
    float _S1170 = s_primal_ctx_dot_0(_S1165, v1v0_2);
    float v_23 = d_2 * _S1170;
    float3  _S1171 = - _S1164;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1171, rov0_2);
    bool _S1172;
    if(u_23 >= 0.0f)
    {
        _S1172 = v_23 >= 0.0f;
    }
    else
    {
        _S1172 = false;
    }
    if(_S1172)
    {
        _S1172 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1172 = false;
    }
    if(_S1172)
    {
        _S1172 = t_8 >= 0.0f;
    }
    else
    {
        _S1172 = false;
    }
    bool _S1173 = !!_S1172;
    float _S1174;
    float _S1175;
    float _S1176;
    float _S1177;
    float _S1178;
    float _S1179;
    float _S1180;
    float _S1181;
    float _S1182;
    float _S1183;
    float _S1184;
    if(_S1173)
    {
        float _S1185 = (F32_min((u_23), (v_23)));
        float _S1186 = s_primal_ctx_sqrt_0(0.5f);
        float _S1187 = _S1186 * (1.0f - u_23 - v_23);
        float _S1188 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1185), (_S1187))) * _S1188;
        float _S1189 = _S1162.primal_0.y;
        float _S1190 = 1.0f - opac_1;
        float _S1191 = 1.0f - s_primal_ctx_clamp_0(_S1189, 0.0f, 0.99989998340606689f);
        float _S1192 = 1.0f / _S1191;
        float _S1193 = _S1191 * _S1191;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1190, _S1192);
        float o_1 = _S1162.primal_0.x;
        bool _S1194 = opac_1 < 0.0f;
        if(_S1194)
        {
            _S1174 = 0.0f;
        }
        else
        {
            _S1174 = o_1 * w_1;
        }
        _S1172 = _S1194;
        _S1175 = o_1;
        _S1176 = w_1;
        _S1177 = _S1190;
        _S1178 = _S1192;
        _S1179 = _S1193;
        _S1180 = _S1189;
        _S1181 = _S1188;
        _S1182 = _S1185;
        _S1183 = _S1187;
        _S1184 = _S1186;
    }
    else
    {
        _S1172 = false;
        _S1174 = 0.0f;
        _S1175 = 0.0f;
        _S1176 = 0.0f;
        _S1177 = 0.0f;
        _S1178 = 0.0f;
        _S1179 = 0.0f;
        _S1180 = 0.0f;
        _S1181 = 0.0f;
        _S1182 = 0.0f;
        _S1183 = 0.0f;
        _S1184 = 0.0f;
    }
    float2  _S1195 = make_float2 (0.0f);
    float2  _S1196;
    if(_S1173)
    {
        if(_S1172)
        {
            _S1174 = 0.0f;
            _S1175 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1197;
            (&_S1197)->primal_0 = _S1174;
            (&_S1197)->differential_0 = 0.0f;
            DiffPair_float_0 _S1198;
            (&_S1198)->primal_0 = 0.99500000476837158f;
            (&_S1198)->differential_0 = 0.0f;
            _d_min_0(&_S1197, &_S1198, _s_dOut_5);
            float _S1199 = _S1175 * _S1197.differential_0;
            _S1174 = _S1176 * _S1197.differential_0;
            _S1175 = _S1199;
        }
        float _S1200 = - _S1175;
        DiffPair_float_0 _S1201;
        (&_S1201)->primal_0 = _S1177;
        (&_S1201)->differential_0 = 0.0f;
        DiffPair_float_0 _S1202;
        (&_S1202)->primal_0 = _S1178;
        (&_S1202)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1201, &_S1202, _S1200);
        float _S1203 = - - (_S1202.differential_0 / _S1179);
        float s_diff_opac_T_0 = - _S1201.differential_0;
        DiffPair_float_0 _S1204;
        (&_S1204)->primal_0 = _S1180;
        (&_S1204)->differential_0 = 0.0f;
        DiffPair_float_0 _S1205;
        (&_S1205)->primal_0 = 0.0f;
        (&_S1205)->differential_0 = 0.0f;
        DiffPair_float_0 _S1206;
        (&_S1206)->primal_0 = 0.99989998340606689f;
        (&_S1206)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1204, &_S1205, &_S1206, _S1203);
        float _S1207 = _S1181 * s_diff_opac_T_0;
        DiffPair_float_0 _S1208;
        (&_S1208)->primal_0 = _S1182;
        (&_S1208)->differential_0 = 0.0f;
        DiffPair_float_0 _S1209;
        (&_S1209)->primal_0 = _S1183;
        (&_S1209)->differential_0 = 0.0f;
        _d_min_0(&_S1208, &_S1209, _S1207);
        float _S1210 = - (_S1184 * _S1209.differential_0);
        DiffPair_float_0 _S1211;
        (&_S1211)->primal_0 = u_23;
        (&_S1211)->differential_0 = 0.0f;
        DiffPair_float_0 _S1212;
        (&_S1212)->primal_0 = v_23;
        (&_S1212)->differential_0 = 0.0f;
        _d_min_0(&_S1211, &_S1212, _S1208.differential_0);
        float2  _S1213 = make_float2 (_S1174, _S1204.differential_0);
        float _S1214 = _S1210 + _S1212.differential_0;
        _S1174 = _S1210 + _S1211.differential_0;
        _S1175 = _S1214;
        _S1196 = _S1213;
    }
    else
    {
        _S1174 = 0.0f;
        _S1175 = 0.0f;
        _S1196 = _S1195;
    }
    float3  _S1215 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1216;
    (&_S1216)->primal_0 = _S1171;
    (&_S1216)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1217;
    (&_S1217)->primal_0 = rov0_2;
    (&_S1217)->differential_0 = _S1215;
    s_bwd_prop_dot_0(&_S1216, &_S1217, 0.0f);
    float3  _S1218 = - _S1216.differential_0;
    float _S1219 = d_2 * _S1175;
    float _S1220 = _S1170 * _S1175;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1221;
    (&_S1221)->primal_0 = _S1165;
    (&_S1221)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1222;
    (&_S1222)->primal_0 = v1v0_2;
    (&_S1222)->differential_0 = _S1215;
    s_bwd_prop_dot_0(&_S1221, &_S1222, _S1219);
    float _S1223 = d_2 * _S1174;
    float _S1224 = _S1169 * _S1174;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1225;
    (&_S1225)->primal_0 = _S1168;
    (&_S1225)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1226;
    (&_S1226)->primal_0 = v2v0_2;
    (&_S1226)->differential_0 = _S1215;
    s_bwd_prop_dot_0(&_S1225, &_S1226, _S1223);
    float3  _S1227 = - _S1225.differential_0;
    float _S1228 = - ((_S1220 + _S1224) / _S1167);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1229;
    (&_S1229)->primal_0 = _S1163.primal_0;
    (&_S1229)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1230;
    (&_S1230)->primal_0 = _S1164;
    (&_S1230)->differential_0 = _S1215;
    s_bwd_prop_dot_0(&_S1229, &_S1230, _S1228);
    float3  _S1231 = _S1221.differential_0 + _S1227;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1232;
    (&_S1232)->primal_0 = rov0_2;
    (&_S1232)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1233;
    (&_S1233)->primal_0 = _S1163.primal_0;
    (&_S1233)->differential_0 = _S1215;
    s_bwd_prop_cross_0(&_S1232, &_S1233, _S1231);
    float3  _S1234 = _S1218 + _S1230.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1235;
    (&_S1235)->primal_0 = v1v0_2;
    (&_S1235)->differential_0 = _S1215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1236;
    (&_S1236)->primal_0 = v2v0_2;
    (&_S1236)->differential_0 = _S1215;
    s_bwd_prop_cross_0(&_S1235, &_S1236, _S1234);
    float3  _S1237 = _S1217.differential_0 + _S1232.differential_0;
    float3  _S1238 = _S1226.differential_0 + _S1236.differential_0;
    float3  _S1239 = _S1222.differential_0 + _S1235.differential_0;
    float3  _S1240 = - _S1237 + - _S1238 + - _S1239;
    float3  _S1241 = _S1229.differential_0 + _S1233.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1241;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1237;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1196;
    FixedArray<float3 , 3>  _S1242;
    _S1242[int(0)] = _S1215;
    _S1242[int(1)] = _S1215;
    _S1242[int(2)] = _S1215;
    _S1242[int(2)] = _S1238;
    _S1242[int(0)] = _S1240;
    _S1242[int(1)] = _S1239;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1242;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1243, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1244, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1245, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1246, float _S1247)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1243, _S1244, _S1245, _S1246, _S1247);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1248 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1249 = { _S1248, _S1248, _S1248 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1249;
    float2  _S1250 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1248;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1248;
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
    float3  _S1251 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1252 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1253 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1251);
    float d_4 = 1.0f / _S1253;
    float _S1254 = _S1253 * _S1253;
    float3  _S1255 = - _S1252;
    float _S1256 = s_primal_ctx_dot_0(_S1255, v2v0_4);
    float u_25 = d_4 * _S1256;
    float _S1257 = s_primal_ctx_dot_0(_S1252, v1v0_4);
    float v_25 = d_4 * _S1257;
    float3  _S1258 = - _S1251;
    float3  _S1259 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1260 = make_float3 (v_25) * dpcolor_0;
    float3  _S1261 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1262 = make_float3 (u_25) * dpcolor_0;
    float3  _S1263 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1264 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1265 = - (_S1263.x + _S1263.y + _S1263.z);
    float _S1266 = d_4 * dpdepth_0;
    float _S1267 = s_primal_ctx_dot_0(_S1258, rov0_4) * dpdepth_0;
    float3  _S1268 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1269;
    (&_S1269)->primal_0 = _S1258;
    (&_S1269)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1270;
    (&_S1270)->primal_0 = rov0_4;
    (&_S1270)->differential_0 = _S1268;
    s_bwd_prop_dot_0(&_S1269, &_S1270, _S1266);
    float3  _S1271 = - _S1269.differential_0;
    float _S1272 = _S1265 + _S1259.x + _S1259.y + _S1259.z;
    float _S1273 = d_4 * _S1272;
    float _S1274 = _S1257 * _S1272;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1275;
    (&_S1275)->primal_0 = _S1252;
    (&_S1275)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1276;
    (&_S1276)->primal_0 = v1v0_4;
    (&_S1276)->differential_0 = _S1268;
    s_bwd_prop_dot_0(&_S1275, &_S1276, _S1273);
    float _S1277 = _S1265 + _S1261.x + _S1261.y + _S1261.z;
    float _S1278 = d_4 * _S1277;
    float _S1279 = _S1256 * _S1277;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1280;
    (&_S1280)->primal_0 = _S1255;
    (&_S1280)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1281;
    (&_S1281)->primal_0 = v2v0_4;
    (&_S1281)->differential_0 = _S1268;
    s_bwd_prop_dot_0(&_S1280, &_S1281, _S1278);
    float3  _S1282 = - _S1280.differential_0;
    float _S1283 = - ((_S1267 + _S1274 + _S1279) / _S1254);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1284;
    (&_S1284)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1284)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1285;
    (&_S1285)->primal_0 = _S1251;
    (&_S1285)->differential_0 = _S1268;
    s_bwd_prop_dot_0(&_S1284, &_S1285, _S1283);
    float3  _S1286 = _S1275.differential_0 + _S1282;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1287;
    (&_S1287)->primal_0 = rov0_4;
    (&_S1287)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1288;
    (&_S1288)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1288)->differential_0 = _S1268;
    s_bwd_prop_cross_0(&_S1287, &_S1288, _S1286);
    float3  _S1289 = _S1271 + _S1285.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1290;
    (&_S1290)->primal_0 = v1v0_4;
    (&_S1290)->differential_0 = _S1268;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1291;
    (&_S1291)->primal_0 = v2v0_4;
    (&_S1291)->differential_0 = _S1268;
    s_bwd_prop_cross_0(&_S1290, &_S1291, _S1289);
    float3  _S1292 = _S1270.differential_0 + _S1287.differential_0;
    float3  _S1293 = _S1281.differential_0 + _S1291.differential_0;
    float3  _S1294 = _S1276.differential_0 + _S1290.differential_0;
    float3  _S1295 = - _S1292 + - _S1293 + - _S1294;
    float3  _S1296 = _S1284.differential_0 + _S1288.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1296;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1292;
    FixedArray<float3 , 3>  _S1297;
    _S1297[int(0)] = _S1268;
    _S1297[int(1)] = _S1268;
    _S1297[int(2)] = _S1268;
    _S1297[int(2)] = _S1260;
    _S1297[int(1)] = _S1262;
    _S1297[int(0)] = _S1264;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1297;
    FixedArray<float3 , 3>  _S1298;
    _S1298[int(0)] = _S1268;
    _S1298[int(1)] = _S1268;
    _S1298[int(2)] = _S1268;
    _S1298[int(2)] = _S1293;
    _S1298[int(0)] = _S1295;
    _S1298[int(1)] = _S1294;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1298;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1299, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1300, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1301, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1302, float3  _S1303, float _S1304)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1299, _S1300, _S1301, _S1302, _S1303, _S1304);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1305 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1306 = { _S1305, _S1305, _S1305 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1306;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1306;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1305;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1305;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

