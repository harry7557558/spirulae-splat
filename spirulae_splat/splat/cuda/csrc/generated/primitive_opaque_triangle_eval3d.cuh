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

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_0, float4  quat_0, float3  scale_0, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_0, FixedArray<float3 , 2>  ch_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * depth_0, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_0)
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
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_12) * (*sh_coeffs_0)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_4) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_12) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_12 * fS1_0 + y_4 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_4) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_12) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_12 * fC1_0 - y_4 * fS1_0)) * (*sh_coeffs_0)[int(15)]);
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

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_1, float4  quat_1, float3  scale_1, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_1, FixedArray<float3 , 2>  ch_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * depth_1, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_1)
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
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_1)[int(1)] + make_float3 (z_2) * (*sh_coeffs_1)[int(2)] - make_float3 (x_15) * (*sh_coeffs_1)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_6) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_15) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_1 + y_6 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_6) * (*sh_coeffs_1)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_15) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_1 - y_6 * fS1_1)) * (*sh_coeffs_1)[int(15)]);
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

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
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
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_2)[int(1)] + make_float3 (z_3) * (*sh_coeffs_2)[int(2)] - make_float3 (x_17) * (*sh_coeffs_2)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_7) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_17) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_17 * fS1_2 + y_7 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_7) * (*sh_coeffs_2)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_17) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_2 - y_7 * fS1_2)) * (*sh_coeffs_2)[int(15)]);
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

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
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
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_3)[int(1)] + make_float3 (z_5) * (*sh_coeffs_3)[int(2)] - make_float3 (x_20) * (*sh_coeffs_3)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_9) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_20) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_20 * fS1_3 + y_9 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_9) * (*sh_coeffs_3)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_20) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_3 - y_9 * fS1_3)) * (*sh_coeffs_3)[int(15)]);
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

struct s_bwd_prop_projection_opaque_triangle_eval3d_persp_differentiable_Intermediates_0
{
    FixedArray<float3 , 16>  _S271;
};

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

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_12, float _s_dOut_3)
{
    float _S330 = (*dpx_12).primal_0.x;
    float _S331 = (*dpx_12).primal_0.y;
    float _S332 = (*dpx_12).primal_0.z;
    float _S333 = (*dpx_12).primal_0.w;
    DiffPair_float_0 _S334;
    (&_S334)->primal_0 = _S330 * _S330 + _S331 * _S331 + _S332 * _S332 + _S333 * _S333;
    (&_S334)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S334, _s_dOut_3);
    float _S335 = (*dpx_12).primal_0.w * _S334.differential_0;
    float _S336 = _S335 + _S335;
    float _S337 = (*dpx_12).primal_0.z * _S334.differential_0;
    float _S338 = _S337 + _S337;
    float _S339 = (*dpx_12).primal_0.y * _S334.differential_0;
    float _S340 = _S339 + _S339;
    float _S341 = (*dpx_12).primal_0.x * _S334.differential_0;
    float _S342 = _S341 + _S341;
    float4  _S343 = make_float4 (0.0f);
    *&((&_S343)->w) = _S336;
    *&((&_S343)->z) = _S338;
    *&((&_S343)->y) = _S340;
    *&((&_S343)->x) = _S342;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S343;
    return;
}

inline __device__ void s_bwd_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S344, float _S345)
{
    s_bwd_prop_length_impl_2(_S344, _S345);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_13, float4  _s_dOut_4)
{
    float _S346 = length_2((*dpx_13).primal_0);
    float4  _S347 = (*dpx_13).primal_0 * _s_dOut_4;
    float4  _S348 = make_float4 (1.0f / _S346) * _s_dOut_4;
    float _S349 = - ((_S347.x + _S347.y + _S347.z + _S347.w) / (_S346 * _S346));
    float4  _S350 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S351;
    (&_S351)->primal_0 = (*dpx_13).primal_0;
    (&_S351)->differential_0 = _S350;
    s_bwd_length_impl_2(&_S351, _S349);
    float4  _S352 = _S348 + _S351.differential_0;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S352;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S353, float4  _S354)
{
    s_bwd_prop_normalize_impl_1(_S353, _S354);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S355, float _S356)
{
    _d_exp_0(_S355, _S356);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  _S357 = make_float3 (0.0f);
    FixedArray<float3 , 16>  _S358 = {
        _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357, _S357
    };
    s_bwd_prop_projection_opaque_triangle_eval3d_persp_differentiable_Intermediates_0 _S359;
    (&_S359)->_S271 = _S358;
    (&_S359)->_S271 = *sh_coeffs_4;
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S360 = scale_4.x;
    float _S361 = s_primal_ctx_exp_0(_S360);
    float _S362 = scale_4.y;
    float _S363 = s_primal_ctx_exp_0(_S362);
    float sz_4 = scale_4.z - 0.5f * (_S360 + _S362);
    float4  _S364 = normalize_0(quat_4);
    float _S365 = _S364.y;
    float x2_4 = _S365 * _S365;
    float y2_4 = _S364.z * _S364.z;
    float z2_8 = _S364.w * _S364.w;
    float xy_4 = _S364.y * _S364.z;
    float xz_4 = _S364.y * _S364.w;
    float yz_4 = _S364.z * _S364.w;
    float wx_4 = _S364.x * _S364.y;
    float wy_4 = _S364.x * _S364.z;
    float wz_4 = _S364.x * _S364.w;
    Matrix<float, 3, 3>  _S366 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S367 = make_float3 (_S361, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S366, _S367) + mean_4;
    float _S368 = -0.5f + sz_4;
    float3  _S369 = make_float3 (_S361 * _S368, _S363, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S366, _S369) + mean_4;
    float _S370 = -0.5f - sz_4;
    float3  _S371 = make_float3 (_S361 * _S370, - _S363, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S366, _S371) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S372 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S373 = vert0_c_4.z;
    float2  _S374 = make_float2 (_S373);
    float2  _S375 = _S372 / make_float2 (_S373);
    float2  _S376 = make_float2 (_S373 * _S373);
    float u_15 = _S375.x;
    float v_15 = _S375.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S377 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S378 = dist_coeffs_4[int(1)] + r2_15 * _S377;
    float _S379 = dist_coeffs_4[int(0)] + r2_15 * _S378;
    float radial_6 = 1.0f + r2_15 * _S379;
    float _S380 = 2.0f * dist_coeffs_4[int(4)];
    float _S381 = _S380 * u_15;
    float _S382 = 2.0f * u_15;
    float _S383 = 2.0f * dist_coeffs_4[int(5)];
    float _S384 = _S383 * u_15;
    float _S385 = 2.0f * v_15;
    float2  _S386 = _S375 * make_float2 (radial_6) + make_float2 (_S381 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S382 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S384 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S385 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S387 = _S386 + make_float2 (dist_coeffs_4[int(8)] * _S386.x + dist_coeffs_4[int(9)] * _S386.y, 0.0f);
    float _S388 = fx_4 * _S387.x + cx_4;
    float _S389 = fy_4 * _S387.y + cy_4;
    float2  uv0_6 = make_float2 (_S388, _S389);
    float2  _S390 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S391 = vert1_c_4.z;
    float2  _S392 = make_float2 (_S391);
    float2  _S393 = _S390 / make_float2 (_S391);
    float2  _S394 = make_float2 (_S391 * _S391);
    float u_16 = _S393.x;
    float v_16 = _S393.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S395 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S396 = dist_coeffs_4[int(1)] + r2_16 * _S395;
    float _S397 = dist_coeffs_4[int(0)] + r2_16 * _S396;
    float radial_7 = 1.0f + r2_16 * _S397;
    float _S398 = _S380 * u_16;
    float _S399 = 2.0f * u_16;
    float _S400 = _S383 * u_16;
    float _S401 = 2.0f * v_16;
    float2  _S402 = _S393 * make_float2 (radial_7) + make_float2 (_S398 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S399 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S400 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S401 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S403 = _S402 + make_float2 (dist_coeffs_4[int(8)] * _S402.x + dist_coeffs_4[int(9)] * _S402.y, 0.0f);
    float _S404 = fx_4 * _S403.x + cx_4;
    float _S405 = fy_4 * _S403.y + cy_4;
    float2  uv1_6 = make_float2 (_S404, _S405);
    float2  _S406 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S407 = vert2_c_4.z;
    float2  _S408 = make_float2 (_S407);
    float2  _S409 = _S406 / make_float2 (_S407);
    float2  _S410 = make_float2 (_S407 * _S407);
    float u_17 = _S409.x;
    float v_17 = _S409.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S411 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S412 = dist_coeffs_4[int(1)] + r2_17 * _S411;
    float _S413 = dist_coeffs_4[int(0)] + r2_17 * _S412;
    float radial_8 = 1.0f + r2_17 * _S413;
    float _S414 = _S380 * u_17;
    float _S415 = 2.0f * u_17;
    float _S416 = _S383 * u_17;
    float _S417 = 2.0f * v_17;
    float2  _S418 = _S409 * make_float2 (radial_8) + make_float2 (_S414 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S415 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S416 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S417 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S419 = _S418 + make_float2 (dist_coeffs_4[int(8)] * _S418.x + dist_coeffs_4[int(9)] * _S418.y, 0.0f);
    float _S420 = fx_4 * _S419.x + cx_4;
    float _S421 = fy_4 * _S419.y + cy_4;
    float2  uv2_6 = make_float2 (_S420, _S421);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S422 = e0_4.x;
    float _S423 = e1_4.y;
    float _S424 = e0_4.y;
    float _S425 = e1_4.x;
    float _S426 = _S422 * _S423 - _S424 * _S425;
    float _S427 = 1.0f - hardness_4.y;
    float _S428 = -1.0f / _S427;
    float _S429 = _S427 * _S427;
    float _S430 = (F32_max((_S388), (_S404)));
    float _S431 = (F32_min((_S388), (_S404)));
    float _S432 = (F32_max((_S389), (_S405)));
    float _S433 = (F32_min((_S389), (_S405)));
    Matrix<float, 3, 3>  _S434 = transpose_1(R_4);
    float3  _S435 = mean_4 - - s_primal_ctx_mul_0(_S434, t_4);
    float _S436 = _S435.x;
    float _S437 = _S435.y;
    float _S438 = _S435.z;
    float _S439 = _S436 * _S436 + _S437 * _S437 + _S438 * _S438;
    float _S440 = s_primal_ctx_sqrt_0(_S439);
    float x_21 = _S436 / _S440;
    float3  _S441 = make_float3 (x_21);
    float _S442 = _S440 * _S440;
    float y_10 = _S437 / _S440;
    float z_6 = _S438 / _S440;
    float3  _S443 = make_float3 (z_6);
    float _S444 = - y_10;
    float3  _S445 = make_float3 (_S444);
    float z2_9 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_4 = x_21 * x_21 - y_10 * y_10;
    float _S446 = 2.0f * x_21;
    float fS1_4 = _S446 * y_10;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S447 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_21;
    float3  _S448 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_10;
    float3  _S449 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S450 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S451 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_6;
    float _S452 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_6 * _S452;
    float3  _S453 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_21;
    float3  _S454 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_10;
    float3  _S455 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S456 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S457 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_21 * fC1_4 - y_10 * fS1_4);
    float3  _S458 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_21 * fS1_4 + y_10 * fC1_4);
    float3  _S459 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * (&_S359)->_S271[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S444) * (&_S359)->_S271[int(1)] + make_float3 (z_6) * (&_S359)->_S271[int(2)] - make_float3 (x_21) * (&_S359)->_S271[int(3)]) + (make_float3 (pSH4_0) * (&_S359)->_S271[int(4)] + make_float3 (pSH5_0) * (&_S359)->_S271[int(5)] + make_float3 (pSH6_0) * (&_S359)->_S271[int(6)] + make_float3 (pSH7_0) * (&_S359)->_S271[int(7)] + make_float3 (pSH8_0) * (&_S359)->_S271[int(8)]) + (make_float3 (pSH9_0) * (&_S359)->_S271[int(9)] + make_float3 (pSH10_0) * (&_S359)->_S271[int(10)] + make_float3 (pSH11_0) * (&_S359)->_S271[int(11)] + make_float3 (pSH12_0) * (&_S359)->_S271[int(12)] + make_float3 (pSH13_0) * (&_S359)->_S271[int(13)] + make_float3 (pSH14_0) * (&_S359)->_S271[int(14)] + make_float3 (pSH15_0) * (&_S359)->_S271[int(15)]);
    float3  _S460 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S461 = make_float3 (0.0f);
    float3  _S462 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S463 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S464 = make_float3 (_S463);
    float3  _S465 = ch_coeffs_4[int(1)] * make_float3 (_S463);
    float3  _S466 = _S462 + _S465 + make_float3 (0.5f);
    float3  _S467 = _S462 - _S465 + make_float3 (0.5f);
    float3  _S468 = vert1_c_4 - vert0_c_4;
    float3  _S469 = vert2_c_4 - vert0_c_4;
    float3  _S470 = s_primal_ctx_cross_0(_S468, _S469);
    float3  _S471 = normalize_1(_S470);
    float3  _S472 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S471, mean_c_4)))))) * v_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S473;
    (&_S473)->primal_0 = _S471;
    (&_S473)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S474;
    (&_S474)->primal_0 = mean_c_4;
    (&_S474)->differential_0 = _S357;
    s_bwd_prop_dot_0(&_S473, &_S474, 0.0f);
    float3  _S475 = _S472 + _S473.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S476;
    (&_S476)->primal_0 = _S470;
    (&_S476)->differential_0 = _S357;
    s_bwd_normalize_impl_0(&_S476, _S475);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S477;
    (&_S477)->primal_0 = _S468;
    (&_S477)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S478;
    (&_S478)->primal_0 = _S469;
    (&_S478)->differential_0 = _S357;
    s_bwd_prop_cross_0(&_S477, &_S478, _S476.differential_0);
    float3  _S479 = - _S478.differential_0;
    float3  _S480 = - _S477.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481;
    (&_S481)->primal_0 = _S467;
    (&_S481)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S482;
    (&_S482)->primal_0 = _S461;
    (&_S482)->differential_0 = _S357;
    s_bwd_prop_max_0(&_S481, &_S482, v_rgbs_0[int(2)]);
    float3  _S483 = - _S481.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S466;
    (&_S484)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S485;
    (&_S485)->primal_0 = _S461;
    (&_S485)->differential_0 = _S357;
    s_bwd_prop_max_0(&_S484, &_S485, v_rgbs_0[int(1)]);
    float3  _S486 = _S464 * (_S483 + _S484.differential_0);
    float3  _S487 = _S481.differential_0 + _S484.differential_0;
    float3  _S488 = make_float3 (0.5f) * - _S487;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S489;
    (&_S489)->primal_0 = _S460;
    (&_S489)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S490;
    (&_S490)->primal_0 = _S461;
    (&_S490)->differential_0 = _S357;
    s_bwd_prop_max_0(&_S489, &_S490, v_rgbs_0[int(0)]);
    float3  _S491 = _S488 + _S489.differential_0;
    float3  _S492 = _S487 + _S489.differential_0;
    float3  _S493 = _S458 * _S492;
    float3  _S494 = (&_S359)->_S271[int(15)] * _S492;
    float3  _S495 = _S456 * _S492;
    float3  _S496 = (&_S359)->_S271[int(14)] * _S492;
    float3  _S497 = _S454 * _S492;
    float3  _S498 = (&_S359)->_S271[int(13)] * _S492;
    float3  _S499 = _S453 * _S492;
    float3  _S500 = (&_S359)->_S271[int(12)] * _S492;
    float3  _S501 = _S455 * _S492;
    float3  _S502 = (&_S359)->_S271[int(11)] * _S492;
    float3  _S503 = _S457 * _S492;
    float3  _S504 = (&_S359)->_S271[int(10)] * _S492;
    float3  _S505 = _S459 * _S492;
    float3  _S506 = (&_S359)->_S271[int(9)] * _S492;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S506.x + _S506.y + _S506.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S494.x + _S494.y + _S494.z);
    float _S507 = _S504.x + _S504.y + _S504.z;
    float _S508 = _S496.x + _S496.y + _S496.z;
    float _S509 = _S502.x + _S502.y + _S502.z;
    float _S510 = _S498.x + _S498.y + _S498.z;
    float _S511 = _S500.x + _S500.y + _S500.z;
    float _S512 = - s_diff_fC2_T_0;
    float3  _S513 = _S450 * _S492;
    float3  _S514 = (&_S359)->_S271[int(8)] * _S492;
    float3  _S515 = _S448 * _S492;
    float3  _S516 = (&_S359)->_S271[int(7)] * _S492;
    float3  _S517 = _S447 * _S492;
    float3  _S518 = (&_S359)->_S271[int(6)] * _S492;
    float3  _S519 = _S449 * _S492;
    float3  _S520 = (&_S359)->_S271[int(5)] * _S492;
    float3  _S521 = _S451 * _S492;
    float3  _S522 = (&_S359)->_S271[int(4)] * _S492;
    float _S523 = _S520.x + _S520.y + _S520.z;
    float _S524 = _S516.x + _S516.y + _S516.z;
    float _S525 = fTmp1B_4 * _S507 + x_21 * s_diff_fS2_T_0 + y_10 * _S512 + 0.54627424478530884f * (_S522.x + _S522.y + _S522.z);
    float _S526 = fTmp1B_4 * _S508 + y_10 * s_diff_fS2_T_0 + x_21 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S514.x + _S514.y + _S514.z);
    float _S527 = y_10 * - _S526;
    float _S528 = x_21 * _S526;
    float _S529 = z_6 * (1.86588168144226074f * (z_6 * _S511) + -2.28522896766662598f * (y_10 * _S509 + x_21 * _S510) + 0.94617468118667603f * (_S518.x + _S518.y + _S518.z));
    float3  _S530 = make_float3 (0.48860251903533936f) * _S492;
    float3  _S531 = - _S530;
    float3  _S532 = _S441 * _S531;
    float3  _S533 = (&_S359)->_S271[int(3)] * _S531;
    float3  _S534 = _S443 * _S530;
    float3  _S535 = (&_S359)->_S271[int(2)] * _S530;
    float3  _S536 = _S445 * _S530;
    float3  _S537 = (&_S359)->_S271[int(1)] * _S530;
    float _S538 = (_S452 * _S511 + 1.44530570507049561f * (fS1_4 * _S507 + fC1_4 * _S508) + -1.09254848957061768f * (y_10 * _S523 + x_21 * _S524) + _S529 + _S529 + _S535.x + _S535.y + _S535.z) / _S442;
    float _S539 = _S440 * _S538;
    float _S540 = (fTmp0C_4 * _S509 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S512 + fTmp0B_4 * _S523 + _S446 * _S525 + _S527 + _S527 + - (_S537.x + _S537.y + _S537.z)) / _S442;
    float _S541 = _S440 * _S540;
    float _S542 = (fTmp0C_4 * _S510 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S524 + 2.0f * (y_10 * _S525) + _S528 + _S528 + _S533.x + _S533.y + _S533.z) / _S442;
    float _S543 = _S440 * _S542;
    float _S544 = _S438 * - _S538 + _S437 * - _S540 + _S436 * - _S542;
    DiffPair_float_0 _S545;
    (&_S545)->primal_0 = _S439;
    (&_S545)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S545, _S544);
    float _S546 = _S438 * _S545.differential_0;
    float _S547 = _S437 * _S545.differential_0;
    float _S548 = _S436 * _S545.differential_0;
    float3  _S549 = make_float3 (0.282094806432724f) * _S492;
    float3  _S550 = make_float3 (_S543 + _S548 + _S548, _S541 + _S547 + _S547, _S539 + _S546 + _S546);
    float3  _S551 = - - _S550;
    Matrix<float, 3, 3>  _S552 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S553;
    (&_S553)->primal_0 = _S434;
    (&_S553)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S554;
    (&_S554)->primal_0 = t_4;
    (&_S554)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S553, &_S554, _S551);
    Matrix<float, 3, 3>  _S555 = transpose_1(_S553.differential_0);
    float3  _S556 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_0);
    DiffPair_float_0 _S557;
    (&_S557)->primal_0 = _S433;
    (&_S557)->differential_0 = 0.0f;
    DiffPair_float_0 _S558;
    (&_S558)->primal_0 = _S421;
    (&_S558)->differential_0 = 0.0f;
    _d_min_0(&_S557, &_S558, 0.0f);
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = _S389;
    (&_S559)->differential_0 = 0.0f;
    DiffPair_float_0 _S560;
    (&_S560)->primal_0 = _S405;
    (&_S560)->differential_0 = 0.0f;
    _d_min_0(&_S559, &_S560, _S557.differential_0);
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S432;
    (&_S561)->differential_0 = 0.0f;
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S421;
    (&_S562)->differential_0 = 0.0f;
    _d_max_0(&_S561, &_S562, 0.0f);
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = _S389;
    (&_S563)->differential_0 = 0.0f;
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = _S405;
    (&_S564)->differential_0 = 0.0f;
    _d_max_0(&_S563, &_S564, _S561.differential_0);
    DiffPair_float_0 _S565;
    (&_S565)->primal_0 = _S431;
    (&_S565)->differential_0 = 0.0f;
    DiffPair_float_0 _S566;
    (&_S566)->primal_0 = _S420;
    (&_S566)->differential_0 = 0.0f;
    _d_min_0(&_S565, &_S566, 0.0f);
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S388;
    (&_S567)->differential_0 = 0.0f;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S404;
    (&_S568)->differential_0 = 0.0f;
    _d_min_0(&_S567, &_S568, _S565.differential_0);
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S430;
    (&_S569)->differential_0 = 0.0f;
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S420;
    (&_S570)->differential_0 = 0.0f;
    _d_max_0(&_S569, &_S570, 0.0f);
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S388;
    (&_S571)->differential_0 = 0.0f;
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S404;
    (&_S572)->differential_0 = 0.0f;
    _d_max_0(&_S571, &_S572, _S569.differential_0);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S428;
    (&_S573)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S573, -0.0f);
    float _S574 = - (-1.0f * - (_S573.differential_0 / _S429));
    float2  _S575 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S576;
    (&_S576)->primal_0 = e2_0;
    (&_S576)->differential_0 = _S575;
    s_bwd_length_impl_1(&_S576, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S577;
    (&_S577)->primal_0 = e1_4;
    (&_S577)->differential_0 = _S575;
    s_bwd_length_impl_1(&_S577, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S578;
    (&_S578)->primal_0 = e0_4;
    (&_S578)->differential_0 = _S575;
    s_bwd_length_impl_1(&_S578, 0.0f);
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S426;
    (&_S579)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S579, -0.0f);
    float _S580 = - _S579.differential_0;
    float2  _S581 = _S577.differential_0 + make_float2 (_S424 * _S580, _S422 * _S579.differential_0);
    float2  _S582 = _S578.differential_0 + make_float2 (_S423 * _S579.differential_0, _S425 * _S580);
    float2  _S583 = - _S576.differential_0 + _S581;
    float _S584 = fx_4 * (_S566.differential_0 + _S570.differential_0 + _S583.x);
    float2  _S585 = make_float2 (_S584, fy_4 * (_S558.differential_0 + _S562.differential_0 + _S583.y)) + make_float2 (dist_coeffs_4[int(8)] * _S584, dist_coeffs_4[int(9)] * _S584);
    float2  _S586 = _S409 * _S585;
    float _S587 = dist_coeffs_4[int(4)] * _S585.y;
    float _S588 = dist_coeffs_4[int(5)] * _S585.x;
    float _S589 = _S586.x + _S586.y;
    float _S590 = r2_17 * _S589;
    float _S591 = r2_17 * _S590;
    float _S592 = dist_coeffs_4[int(7)] * _S585.y + _S587 + dist_coeffs_4[int(6)] * _S585.x + _S588 + _S413 * _S589 + _S412 * _S590 + _S411 * _S591 + dist_coeffs_4[int(3)] * (r2_17 * _S591);
    float _S593 = v_17 * _S592;
    float _S594 = u_17 * _S592;
    float2  _S595 = (make_float2 (radial_8) * _S585 + make_float2 (_S383 * (v_17 * _S585.y) + _S415 * _S588 + 2.0f * (u_17 * _S588) + _S380 * (v_17 * _S585.x) + _S594 + _S594, _S417 * _S587 + 2.0f * (v_17 * _S587) + _S416 * _S585.y + _S414 * _S585.x + _S593 + _S593)) / _S410;
    float2  _S596 = _S406 * - _S595;
    float2  _S597 = _S408 * _S595;
    float2  _S598 = - _S581 + _S582;
    float _S599 = fx_4 * (_S568.differential_0 + _S572.differential_0 + _S598.x);
    float2  _S600 = make_float2 (_S599, fy_4 * (_S560.differential_0 + _S564.differential_0 + _S598.y)) + make_float2 (dist_coeffs_4[int(8)] * _S599, dist_coeffs_4[int(9)] * _S599);
    float2  _S601 = _S393 * _S600;
    float _S602 = dist_coeffs_4[int(4)] * _S600.y;
    float _S603 = dist_coeffs_4[int(5)] * _S600.x;
    float _S604 = _S601.x + _S601.y;
    float _S605 = r2_16 * _S604;
    float _S606 = r2_16 * _S605;
    float _S607 = dist_coeffs_4[int(7)] * _S600.y + _S602 + dist_coeffs_4[int(6)] * _S600.x + _S603 + _S397 * _S604 + _S396 * _S605 + _S395 * _S606 + dist_coeffs_4[int(3)] * (r2_16 * _S606);
    float _S608 = v_16 * _S607;
    float _S609 = u_16 * _S607;
    float2  _S610 = (make_float2 (radial_7) * _S600 + make_float2 (_S383 * (v_16 * _S600.y) + _S399 * _S603 + 2.0f * (u_16 * _S603) + _S380 * (v_16 * _S600.x) + _S609 + _S609, _S401 * _S602 + 2.0f * (v_16 * _S602) + _S400 * _S600.y + _S398 * _S600.x + _S608 + _S608)) / _S394;
    float2  _S611 = _S390 * - _S610;
    float2  _S612 = _S392 * _S610;
    float _S613 = _S611.x + _S611.y;
    float2  _S614 = _S576.differential_0 + - _S582;
    float _S615 = fx_4 * (_S567.differential_0 + _S571.differential_0 + _S614.x);
    float2  _S616 = make_float2 (_S615, fy_4 * (_S559.differential_0 + _S563.differential_0 + _S614.y)) + make_float2 (dist_coeffs_4[int(8)] * _S615, dist_coeffs_4[int(9)] * _S615);
    float2  _S617 = _S375 * _S616;
    float _S618 = dist_coeffs_4[int(4)] * _S616.y;
    float _S619 = dist_coeffs_4[int(5)] * _S616.x;
    float _S620 = _S617.x + _S617.y;
    float _S621 = r2_15 * _S620;
    float _S622 = r2_15 * _S621;
    float _S623 = dist_coeffs_4[int(7)] * _S616.y + _S618 + dist_coeffs_4[int(6)] * _S616.x + _S619 + _S379 * _S620 + _S378 * _S621 + _S377 * _S622 + dist_coeffs_4[int(3)] * (r2_15 * _S622);
    float _S624 = v_15 * _S623;
    float _S625 = u_15 * _S623;
    float2  _S626 = (make_float2 (radial_6) * _S616 + make_float2 (_S383 * (v_15 * _S616.y) + _S382 * _S619 + 2.0f * (u_15 * _S619) + _S380 * (v_15 * _S616.x) + _S625 + _S625, _S385 * _S618 + 2.0f * (v_15 * _S618) + _S384 * _S616.y + _S381 * _S616.x + _S624 + _S624)) / _S376;
    float2  _S627 = _S372 * - _S626;
    float2  _S628 = _S374 * _S626;
    float _S629 = _S627.x + _S627.y;
    float3  _S630 = _S478.differential_0 + _S556 + make_float3 (_S597.x, _S597.y, _S596.x + _S596.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S631;
    (&_S631)->primal_0 = R_4;
    (&_S631)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
    (&_S632)->primal_0 = vert2_4;
    (&_S632)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S631, &_S632, _S630);
    float3  _S633 = _S477.differential_0 + _S556 + make_float3 (_S612.x, _S612.y, _S613);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S634;
    (&_S634)->primal_0 = R_4;
    (&_S634)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S635;
    (&_S635)->primal_0 = vert1_4;
    (&_S635)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S634, &_S635, _S633);
    float3  _S636 = _S479 + _S480 + _S556 + make_float3 (_S628.x, _S628.y, _S629);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S637;
    (&_S637)->primal_0 = R_4;
    (&_S637)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
    (&_S638)->primal_0 = vert0_4;
    (&_S638)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S637, &_S638, _S636);
    float3  _S639 = v_verts_0[int(2)] + _S632.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S640;
    (&_S640)->primal_0 = _S366;
    (&_S640)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S641;
    (&_S641)->primal_0 = _S371;
    (&_S641)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S640, &_S641, _S639);
    float _S642 = - _S641.differential_0.y;
    float _S643 = _S370 * _S641.differential_0.x;
    float _S644 = - (_S361 * _S641.differential_0.x);
    float3  _S645 = v_verts_0[int(1)] + _S635.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S646;
    (&_S646)->primal_0 = _S366;
    (&_S646)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S647;
    (&_S647)->primal_0 = _S369;
    (&_S647)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S646, &_S647, _S645);
    float _S648 = _S361 * _S647.differential_0.x;
    float _S649 = _S368 * _S647.differential_0.x;
    float3  _S650 = v_verts_0[int(0)] + _S638.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S651;
    (&_S651)->primal_0 = _S366;
    (&_S651)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S652;
    (&_S652)->primal_0 = _S367;
    (&_S652)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S651, &_S652, _S650);
    Matrix<float, 3, 3>  _S653 = transpose_1(_S640.differential_0 + _S646.differential_0 + _S651.differential_0);
    float _S654 = 2.0f * - _S653.rows[int(2)].z;
    float _S655 = 2.0f * _S653.rows[int(2)].y;
    float _S656 = 2.0f * _S653.rows[int(2)].x;
    float _S657 = 2.0f * _S653.rows[int(1)].z;
    float _S658 = 2.0f * - _S653.rows[int(1)].y;
    float _S659 = 2.0f * _S653.rows[int(1)].x;
    float _S660 = 2.0f * _S653.rows[int(0)].z;
    float _S661 = 2.0f * _S653.rows[int(0)].y;
    float _S662 = 2.0f * - _S653.rows[int(0)].x;
    float _S663 = - _S659 + _S661;
    float _S664 = _S656 + - _S660;
    float _S665 = - _S655 + _S657;
    float _S666 = _S655 + _S657;
    float _S667 = _S656 + _S660;
    float _S668 = _S659 + _S661;
    float _S669 = _S364.w * (_S658 + _S662);
    float _S670 = _S364.z * (_S654 + _S662);
    float _S671 = _S364.y * (_S654 + _S658);
    float _S672 = _S364.x * _S663 + _S364.z * _S666 + _S364.y * _S667 + _S669 + _S669;
    float _S673 = _S364.x * _S664 + _S364.w * _S666 + _S364.y * _S668 + _S670 + _S670;
    float _S674 = _S364.x * _S665 + _S364.w * _S667 + _S364.z * _S668 + _S671 + _S671;
    float _S675 = _S364.w * _S663 + _S364.z * _S664 + _S364.y * _S665;
    float4  _S676 = make_float4 (0.0f);
    float4  _S677 = _S676;
    *&((&_S677)->w) = _S672;
    *&((&_S677)->z) = _S673;
    *&((&_S677)->y) = _S674;
    *&((&_S677)->x) = _S675;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S678;
    (&_S678)->primal_0 = quat_4;
    (&_S678)->differential_0 = _S676;
    s_bwd_normalize_impl_1(&_S678, _S677);
    float _S679 = _S644 + _S648;
    float _S680 = 0.5f * - _S679;
    float _S681 = _S642 + _S647.differential_0.y;
    DiffPair_float_0 _S682;
    (&_S682)->primal_0 = _S362;
    (&_S682)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S682, _S681);
    float _S683 = _S680 + _S682.differential_0;
    float _S684 = _S643 + _S649 + _S652.differential_0.x;
    DiffPair_float_0 _S685;
    (&_S685)->primal_0 = _S360;
    (&_S685)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S685, _S684);
    float _S686 = _S680 + _S685.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S687;
    (&_S687)->primal_0 = R_4;
    (&_S687)->differential_0 = _S552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S688;
    (&_S688)->primal_0 = mean_4;
    (&_S688)->differential_0 = _S357;
    s_bwd_prop_mul_0(&_S687, &_S688, _S474.differential_0);
    FixedArray<float3 , 16>  _S689;
    _S689[int(0)] = _S357;
    _S689[int(1)] = _S357;
    _S689[int(2)] = _S357;
    _S689[int(3)] = _S357;
    _S689[int(4)] = _S357;
    _S689[int(5)] = _S357;
    _S689[int(6)] = _S357;
    _S689[int(7)] = _S357;
    _S689[int(8)] = _S357;
    _S689[int(9)] = _S357;
    _S689[int(10)] = _S357;
    _S689[int(11)] = _S357;
    _S689[int(12)] = _S357;
    _S689[int(13)] = _S357;
    _S689[int(14)] = _S357;
    _S689[int(15)] = _S357;
    _S689[int(7)] = _S515;
    _S689[int(0)] = _S549;
    _S689[int(1)] = _S536;
    _S689[int(2)] = _S534;
    _S689[int(3)] = _S532;
    _S689[int(4)] = _S521;
    _S689[int(5)] = _S519;
    _S689[int(6)] = _S517;
    _S689[int(8)] = _S513;
    _S689[int(9)] = _S505;
    _S689[int(10)] = _S503;
    _S689[int(11)] = _S501;
    _S689[int(12)] = _S499;
    _S689[int(13)] = _S497;
    _S689[int(14)] = _S495;
    _S689[int(15)] = _S493;
    float3  _S690 = _S554.differential_0 + _S630 + _S633 + _S636 + _S474.differential_0;
    Matrix<float, 3, 3>  _S691 = _S555 + _S631.differential_0 + _S634.differential_0 + _S637.differential_0 + _S687.differential_0;
    FixedArray<float3 , 2>  _S692;
    _S692[int(0)] = _S357;
    _S692[int(1)] = _S357;
    _S692[int(1)] = _S486;
    _S692[int(0)] = _S491;
    float2  _S693 = make_float2 (0.0f, _S574);
    float3  _S694 = make_float3 (_S686, _S683, _S679);
    *v_mean_0 = _S550 + _S639 + _S645 + _S650 + _S688.differential_0;
    *v_quat_0 = _S678.differential_0;
    *v_scale_0 = _S694;
    *v_hardness_0 = _S693;
    (*v_sh_coeffs_0)[int(0)] = _S689[int(0)];
    (*v_sh_coeffs_0)[int(1)] = _S689[int(1)];
    (*v_sh_coeffs_0)[int(2)] = _S689[int(2)];
    (*v_sh_coeffs_0)[int(3)] = _S689[int(3)];
    (*v_sh_coeffs_0)[int(4)] = _S689[int(4)];
    (*v_sh_coeffs_0)[int(5)] = _S689[int(5)];
    (*v_sh_coeffs_0)[int(6)] = _S689[int(6)];
    (*v_sh_coeffs_0)[int(7)] = _S689[int(7)];
    (*v_sh_coeffs_0)[int(8)] = _S689[int(8)];
    (*v_sh_coeffs_0)[int(9)] = _S689[int(9)];
    (*v_sh_coeffs_0)[int(10)] = _S689[int(10)];
    (*v_sh_coeffs_0)[int(11)] = _S689[int(11)];
    (*v_sh_coeffs_0)[int(12)] = _S689[int(12)];
    (*v_sh_coeffs_0)[int(13)] = _S689[int(13)];
    (*v_sh_coeffs_0)[int(14)] = _S689[int(14)];
    (*v_sh_coeffs_0)[int(15)] = _S689[int(15)];
    *v_ch_coeffs_0 = _S692;
    *v_R_0 = _S691;
    *v_t_0 = _S690;
    return;
}

struct s_bwd_prop_projection_opaque_triangle_eval3d_fisheye_differentiable_Intermediates_0
{
    FixedArray<float3 , 16>  _S695;
};

inline __device__ float s_primal_ctx_atan2_0(float _S696, float _S697)
{
    return (F32_atan2((_S696), (_S697)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S698, DiffPair_float_0 * _S699, float _S700)
{
    _d_atan2_0(_S698, _S699, _S700);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  _S701 = make_float3 (0.0f);
    FixedArray<float3 , 16>  _S702 = {
        _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701, _S701
    };
    s_bwd_prop_projection_opaque_triangle_eval3d_fisheye_differentiable_Intermediates_0 _S703;
    (&_S703)->_S695 = _S702;
    (&_S703)->_S695 = *sh_coeffs_5;
    s_bwd_prop_projection_opaque_triangle_eval3d_fisheye_differentiable_Intermediates_0 _S704 = _S703;
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S705 = scale_5.x;
    float _S706 = s_primal_ctx_exp_0(_S705);
    float _S707 = scale_5.y;
    float _S708 = s_primal_ctx_exp_0(_S707);
    float sz_5 = scale_5.z - 0.5f * (_S705 + _S707);
    float4  _S709 = normalize_0(quat_5);
    float _S710 = _S709.y;
    float x2_5 = _S710 * _S710;
    float y2_5 = _S709.z * _S709.z;
    float z2_10 = _S709.w * _S709.w;
    float xy_5 = _S709.y * _S709.z;
    float xz_5 = _S709.y * _S709.w;
    float yz_5 = _S709.z * _S709.w;
    float wx_5 = _S709.x * _S709.y;
    float wy_5 = _S709.x * _S709.z;
    float wz_5 = _S709.x * _S709.w;
    Matrix<float, 3, 3>  _S711 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S712 = make_float3 (_S706, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S711, _S712) + mean_5;
    float _S713 = -0.5f + sz_5;
    float3  _S714 = make_float3 (_S706 * _S713, _S708, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S711, _S714) + mean_5;
    float _S715 = -0.5f - sz_5;
    float3  _S716 = make_float3 (_S706 * _S715, - _S708, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S711, _S716) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S717 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S718 = length_0(_S717);
    float _S719 = vert0_c_5.z;
    float _S720 = s_primal_ctx_atan2_0(_S718, _S719);
    bool _S721 = _S720 < 0.00100000004749745f;
    float k_2;
    float _S722;
    float _S723;
    float _S724;
    if(_S721)
    {
        float _S725 = 1.0f - _S720 * _S720 / 3.0f;
        float _S726 = _S719 * _S719;
        k_2 = _S725 / _S719;
        _S722 = _S726;
        _S723 = _S725;
        _S724 = 0.0f;
    }
    else
    {
        float _S727 = _S718 * _S718;
        k_2 = _S720 / _S718;
        _S722 = 0.0f;
        _S723 = 0.0f;
        _S724 = _S727;
    }
    float2  _S728 = make_float2 (k_2);
    float2  _S729 = _S717 * make_float2 (k_2);
    float u_18 = _S729.x;
    float v_18 = _S729.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S730 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S731 = dist_coeffs_5[int(1)] + r2_18 * _S730;
    float _S732 = dist_coeffs_5[int(0)] + r2_18 * _S731;
    float radial_9 = 1.0f + r2_18 * _S732;
    float _S733 = 2.0f * dist_coeffs_5[int(4)];
    float _S734 = _S733 * u_18;
    float _S735 = 2.0f * u_18;
    float _S736 = 2.0f * dist_coeffs_5[int(5)];
    float _S737 = _S736 * u_18;
    float _S738 = 2.0f * v_18;
    float2  _S739 = _S729 * make_float2 (radial_9) + make_float2 (_S734 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S735 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S737 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S738 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S740 = _S739 + make_float2 (dist_coeffs_5[int(8)] * _S739.x + dist_coeffs_5[int(9)] * _S739.y, 0.0f);
    float _S741 = fx_5 * _S740.x + cx_5;
    float _S742 = fy_5 * _S740.y + cy_5;
    float2  uv0_7 = make_float2 (_S741, _S742);
    float2  _S743 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S744 = length_0(_S743);
    float _S745 = vert1_c_5.z;
    float _S746 = s_primal_ctx_atan2_0(_S744, _S745);
    bool _S747 = _S746 < 0.00100000004749745f;
    float _S748;
    float _S749;
    float _S750;
    if(_S747)
    {
        float _S751 = 1.0f - _S746 * _S746 / 3.0f;
        float _S752 = _S745 * _S745;
        k_2 = _S751 / _S745;
        _S748 = _S752;
        _S749 = _S751;
        _S750 = 0.0f;
    }
    else
    {
        float _S753 = _S744 * _S744;
        k_2 = _S746 / _S744;
        _S748 = 0.0f;
        _S749 = 0.0f;
        _S750 = _S753;
    }
    float2  _S754 = make_float2 (k_2);
    float2  _S755 = _S743 * make_float2 (k_2);
    float u_19 = _S755.x;
    float v_19 = _S755.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S756 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S757 = dist_coeffs_5[int(1)] + r2_19 * _S756;
    float _S758 = dist_coeffs_5[int(0)] + r2_19 * _S757;
    float radial_10 = 1.0f + r2_19 * _S758;
    float _S759 = _S733 * u_19;
    float _S760 = 2.0f * u_19;
    float _S761 = _S736 * u_19;
    float _S762 = 2.0f * v_19;
    float2  _S763 = _S755 * make_float2 (radial_10) + make_float2 (_S759 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S760 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S761 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S762 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S764 = _S763 + make_float2 (dist_coeffs_5[int(8)] * _S763.x + dist_coeffs_5[int(9)] * _S763.y, 0.0f);
    float _S765 = fx_5 * _S764.x + cx_5;
    float _S766 = fy_5 * _S764.y + cy_5;
    float2  uv1_7 = make_float2 (_S765, _S766);
    float2  _S767 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S768 = length_0(_S767);
    float _S769 = vert2_c_5.z;
    float _S770 = s_primal_ctx_atan2_0(_S768, _S769);
    bool _S771 = _S770 < 0.00100000004749745f;
    float _S772;
    float _S773;
    float _S774;
    if(_S771)
    {
        float _S775 = 1.0f - _S770 * _S770 / 3.0f;
        float _S776 = _S769 * _S769;
        k_2 = _S775 / _S769;
        _S772 = _S776;
        _S773 = _S775;
        _S774 = 0.0f;
    }
    else
    {
        float _S777 = _S768 * _S768;
        k_2 = _S770 / _S768;
        _S772 = 0.0f;
        _S773 = 0.0f;
        _S774 = _S777;
    }
    float2  _S778 = make_float2 (k_2);
    float2  _S779 = _S767 * make_float2 (k_2);
    float u_20 = _S779.x;
    float v_20 = _S779.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S780 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S781 = dist_coeffs_5[int(1)] + r2_20 * _S780;
    float _S782 = dist_coeffs_5[int(0)] + r2_20 * _S781;
    float radial_11 = 1.0f + r2_20 * _S782;
    float _S783 = _S733 * u_20;
    float _S784 = 2.0f * u_20;
    float _S785 = _S736 * u_20;
    float _S786 = 2.0f * v_20;
    float2  _S787 = _S779 * make_float2 (radial_11) + make_float2 (_S783 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S784 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S785 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S786 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S788 = _S787 + make_float2 (dist_coeffs_5[int(8)] * _S787.x + dist_coeffs_5[int(9)] * _S787.y, 0.0f);
    float _S789 = fx_5 * _S788.x + cx_5;
    float _S790 = fy_5 * _S788.y + cy_5;
    float2  uv2_7 = make_float2 (_S789, _S790);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S791 = e0_5.x;
    float _S792 = e1_5.y;
    float _S793 = e0_5.y;
    float _S794 = e1_5.x;
    float _S795 = _S791 * _S792 - _S793 * _S794;
    float _S796 = 1.0f - hardness_5.y;
    float _S797 = -1.0f / _S796;
    float _S798 = _S796 * _S796;
    float _S799 = (F32_max((_S741), (_S765)));
    float _S800 = (F32_min((_S741), (_S765)));
    float _S801 = (F32_max((_S742), (_S766)));
    float _S802 = (F32_min((_S742), (_S766)));
    float3  _S803 = (vert0_c_5 + vert1_c_5 + vert2_c_5) / make_float3 (3.0f);
    float x_22 = _S803.x;
    float y_11 = _S803.y;
    float z_7 = _S803.z;
    float _S804 = z_7 * z_7;
    float _S805 = _S804 * z_7;
    float _S806 = x_22 * x_22 + y_11 * y_11;
    float _S807 = 0.001953125f * _S806;
    Matrix<float, 3, 3>  _S808 = transpose_1(R_5);
    float3  _S809 = mean_5 - - s_primal_ctx_mul_0(_S808, t_5);
    float _S810 = _S809.x;
    float _S811 = _S809.y;
    float _S812 = _S809.z;
    float _S813 = _S810 * _S810 + _S811 * _S811 + _S812 * _S812;
    float _S814 = s_primal_ctx_sqrt_0(_S813);
    float x_23 = _S810 / _S814;
    float3  _S815 = make_float3 (x_23);
    float _S816 = _S814 * _S814;
    float y_12 = _S811 / _S814;
    float z_8 = _S812 / _S814;
    float3  _S817 = make_float3 (z_8);
    float _S818 = - y_12;
    float3  _S819 = make_float3 (_S818);
    float z2_11 = z_8 * z_8;
    float fTmp0B_5 = -1.09254848957061768f * z_8;
    float fC1_5 = x_23 * x_23 - y_12 * y_12;
    float _S820 = 2.0f * x_23;
    float fS1_5 = _S820 * y_12;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S821 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_23;
    float3  _S822 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_12;
    float3  _S823 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S824 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S825 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_8;
    float _S826 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_8 * _S826;
    float3  _S827 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_23;
    float3  _S828 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_12;
    float3  _S829 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S830 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S831 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_23 * fC1_5 - y_12 * fS1_5);
    float3  _S832 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_23 * fS1_5 + y_12 * fC1_5);
    float3  _S833 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * _S704._S695[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S818) * _S704._S695[int(1)] + make_float3 (z_8) * _S704._S695[int(2)] - make_float3 (x_23) * _S704._S695[int(3)]) + (make_float3 (pSH4_1) * _S704._S695[int(4)] + make_float3 (pSH5_1) * _S704._S695[int(5)] + make_float3 (pSH6_1) * _S704._S695[int(6)] + make_float3 (pSH7_1) * _S704._S695[int(7)] + make_float3 (pSH8_1) * _S704._S695[int(8)]) + (make_float3 (pSH9_1) * _S704._S695[int(9)] + make_float3 (pSH10_1) * _S704._S695[int(10)] + make_float3 (pSH11_1) * _S704._S695[int(11)] + make_float3 (pSH12_1) * _S704._S695[int(12)] + make_float3 (pSH13_1) * _S704._S695[int(13)] + make_float3 (pSH14_1) * _S704._S695[int(14)] + make_float3 (pSH15_1) * _S704._S695[int(15)]);
    float3  _S834 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S835 = make_float3 (0.0f);
    float3  _S836 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S837 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S838 = make_float3 (_S837);
    float3  _S839 = ch_coeffs_5[int(1)] * make_float3 (_S837);
    float3  _S840 = _S836 + _S839 + make_float3 (0.5f);
    float3  _S841 = _S836 - _S839 + make_float3 (0.5f);
    float3  _S842 = vert1_c_5 - vert0_c_5;
    float3  _S843 = vert2_c_5 - vert0_c_5;
    float3  _S844 = s_primal_ctx_cross_0(_S842, _S843);
    float3  _S845 = normalize_1(_S844);
    float3  _S846 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S845, mean_c_5)))))) * v_normal_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S847;
    (&_S847)->primal_0 = _S845;
    (&_S847)->differential_0 = _S701;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S848;
    (&_S848)->primal_0 = mean_c_5;
    (&_S848)->differential_0 = _S701;
    s_bwd_prop_dot_0(&_S847, &_S848, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S849 = _S848;
    float3  _S850 = _S846 + _S847.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S851;
    (&_S851)->primal_0 = _S844;
    (&_S851)->differential_0 = _S701;
    s_bwd_normalize_impl_0(&_S851, _S850);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S852;
    (&_S852)->primal_0 = _S842;
    (&_S852)->differential_0 = _S701;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S853;
    (&_S853)->primal_0 = _S843;
    (&_S853)->differential_0 = _S701;
    s_bwd_prop_cross_0(&_S852, &_S853, _S851.differential_0);
    float3  _S854 = - _S853.differential_0;
    float3  _S855 = - _S852.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S856;
    (&_S856)->primal_0 = _S841;
    (&_S856)->differential_0 = _S701;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S857;
    (&_S857)->primal_0 = _S835;
    (&_S857)->differential_0 = _S701;
    s_bwd_prop_max_0(&_S856, &_S857, v_rgbs_1[int(2)]);
    float3  _S858 = - _S856.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S859;
    (&_S859)->primal_0 = _S840;
    (&_S859)->differential_0 = _S701;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S860;
    (&_S860)->primal_0 = _S835;
    (&_S860)->differential_0 = _S701;
    s_bwd_prop_max_0(&_S859, &_S860, v_rgbs_1[int(1)]);
    float3  _S861 = _S838 * (_S858 + _S859.differential_0);
    float3  _S862 = _S856.differential_0 + _S859.differential_0;
    float3  _S863 = make_float3 (0.5f) * - _S862;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S864;
    (&_S864)->primal_0 = _S834;
    (&_S864)->differential_0 = _S701;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S865;
    (&_S865)->primal_0 = _S835;
    (&_S865)->differential_0 = _S701;
    s_bwd_prop_max_0(&_S864, &_S865, v_rgbs_1[int(0)]);
    float3  _S866 = _S863 + _S864.differential_0;
    float3  _S867 = _S862 + _S864.differential_0;
    float3  _S868 = _S832 * _S867;
    float3  _S869 = _S704._S695[int(15)] * _S867;
    float3  _S870 = _S830 * _S867;
    float3  _S871 = _S704._S695[int(14)] * _S867;
    float3  _S872 = _S828 * _S867;
    float3  _S873 = _S704._S695[int(13)] * _S867;
    float3  _S874 = _S827 * _S867;
    float3  _S875 = _S704._S695[int(12)] * _S867;
    float3  _S876 = _S829 * _S867;
    float3  _S877 = _S704._S695[int(11)] * _S867;
    float3  _S878 = _S831 * _S867;
    float3  _S879 = _S704._S695[int(10)] * _S867;
    float3  _S880 = _S833 * _S867;
    float3  _S881 = _S704._S695[int(9)] * _S867;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S881.x + _S881.y + _S881.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S869.x + _S869.y + _S869.z);
    float _S882 = _S879.x + _S879.y + _S879.z;
    float _S883 = _S871.x + _S871.y + _S871.z;
    float _S884 = _S877.x + _S877.y + _S877.z;
    float _S885 = _S873.x + _S873.y + _S873.z;
    float _S886 = _S875.x + _S875.y + _S875.z;
    float _S887 = - s_diff_fC2_T_1;
    float3  _S888 = _S824 * _S867;
    float3  _S889 = _S704._S695[int(8)] * _S867;
    float3  _S890 = _S822 * _S867;
    float3  _S891 = _S704._S695[int(7)] * _S867;
    float3  _S892 = _S821 * _S867;
    float3  _S893 = _S704._S695[int(6)] * _S867;
    float3  _S894 = _S823 * _S867;
    float3  _S895 = _S704._S695[int(5)] * _S867;
    float3  _S896 = _S825 * _S867;
    float3  _S897 = _S704._S695[int(4)] * _S867;
    float _S898 = _S895.x + _S895.y + _S895.z;
    float _S899 = _S891.x + _S891.y + _S891.z;
    float _S900 = fTmp1B_5 * _S882 + x_23 * s_diff_fS2_T_1 + y_12 * _S887 + 0.54627424478530884f * (_S897.x + _S897.y + _S897.z);
    float _S901 = fTmp1B_5 * _S883 + y_12 * s_diff_fS2_T_1 + x_23 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S889.x + _S889.y + _S889.z);
    float _S902 = y_12 * - _S901;
    float _S903 = x_23 * _S901;
    float _S904 = z_8 * (1.86588168144226074f * (z_8 * _S886) + -2.28522896766662598f * (y_12 * _S884 + x_23 * _S885) + 0.94617468118667603f * (_S893.x + _S893.y + _S893.z));
    float3  _S905 = make_float3 (0.48860251903533936f) * _S867;
    float3  _S906 = - _S905;
    float3  _S907 = _S815 * _S906;
    float3  _S908 = _S704._S695[int(3)] * _S906;
    float3  _S909 = _S817 * _S905;
    float3  _S910 = _S704._S695[int(2)] * _S905;
    float3  _S911 = _S819 * _S905;
    float3  _S912 = _S704._S695[int(1)] * _S905;
    float _S913 = (_S826 * _S886 + 1.44530570507049561f * (fS1_5 * _S882 + fC1_5 * _S883) + -1.09254848957061768f * (y_12 * _S898 + x_23 * _S899) + _S904 + _S904 + _S910.x + _S910.y + _S910.z) / _S816;
    float _S914 = _S814 * _S913;
    float _S915 = (fTmp0C_5 * _S884 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S887 + fTmp0B_5 * _S898 + _S820 * _S900 + _S902 + _S902 + - (_S912.x + _S912.y + _S912.z)) / _S816;
    float _S916 = _S814 * _S915;
    float _S917 = (fTmp0C_5 * _S885 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S899 + 2.0f * (y_12 * _S900) + _S903 + _S903 + _S908.x + _S908.y + _S908.z) / _S816;
    float _S918 = _S814 * _S917;
    float _S919 = _S812 * - _S913 + _S811 * - _S915 + _S810 * - _S917;
    DiffPair_float_0 _S920;
    (&_S920)->primal_0 = _S813;
    (&_S920)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S920, _S919);
    float _S921 = _S812 * _S920.differential_0;
    float _S922 = _S811 * _S920.differential_0;
    float _S923 = _S810 * _S920.differential_0;
    float3  _S924 = make_float3 (0.282094806432724f) * _S867;
    float3  _S925 = make_float3 (_S918 + _S923 + _S923, _S916 + _S922 + _S922, _S914 + _S921 + _S921);
    float3  _S926 = - - _S925;
    Matrix<float, 3, 3>  _S927 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S928;
    (&_S928)->primal_0 = _S808;
    (&_S928)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S929;
    (&_S929)->primal_0 = t_5;
    (&_S929)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S928, &_S929, _S926);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S930 = _S929;
    Matrix<float, 3, 3>  _S931 = transpose_1(_S928.differential_0);
    float _S932 = _S807 * v_depth_1 + 0.001953125f * (_S806 * v_depth_1);
    float _S933 = y_11 * _S932;
    float _S934 = x_22 * _S932;
    float _S935 = z_7 * v_depth_1;
    float _S936 = z_7 * (z_7 * _S935);
    float3  _S937 = make_float3 (0.3333333432674408f) * make_float3 (_S934 + _S934, _S933 + _S933, _S805 * v_depth_1 + _S804 * _S935 + _S936 + _S936);
    DiffPair_float_0 _S938;
    (&_S938)->primal_0 = _S802;
    (&_S938)->differential_0 = 0.0f;
    DiffPair_float_0 _S939;
    (&_S939)->primal_0 = _S790;
    (&_S939)->differential_0 = 0.0f;
    _d_min_0(&_S938, &_S939, 0.0f);
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S742;
    (&_S940)->differential_0 = 0.0f;
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = _S766;
    (&_S941)->differential_0 = 0.0f;
    _d_min_0(&_S940, &_S941, _S938.differential_0);
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = _S801;
    (&_S942)->differential_0 = 0.0f;
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S790;
    (&_S943)->differential_0 = 0.0f;
    _d_max_0(&_S942, &_S943, 0.0f);
    DiffPair_float_0 _S944;
    (&_S944)->primal_0 = _S742;
    (&_S944)->differential_0 = 0.0f;
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S766;
    (&_S945)->differential_0 = 0.0f;
    _d_max_0(&_S944, &_S945, _S942.differential_0);
    DiffPair_float_0 _S946;
    (&_S946)->primal_0 = _S800;
    (&_S946)->differential_0 = 0.0f;
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S789;
    (&_S947)->differential_0 = 0.0f;
    _d_min_0(&_S946, &_S947, 0.0f);
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = _S741;
    (&_S948)->differential_0 = 0.0f;
    DiffPair_float_0 _S949;
    (&_S949)->primal_0 = _S765;
    (&_S949)->differential_0 = 0.0f;
    _d_min_0(&_S948, &_S949, _S946.differential_0);
    DiffPair_float_0 _S950;
    (&_S950)->primal_0 = _S799;
    (&_S950)->differential_0 = 0.0f;
    DiffPair_float_0 _S951;
    (&_S951)->primal_0 = _S789;
    (&_S951)->differential_0 = 0.0f;
    _d_max_0(&_S950, &_S951, 0.0f);
    DiffPair_float_0 _S952;
    (&_S952)->primal_0 = _S741;
    (&_S952)->differential_0 = 0.0f;
    DiffPair_float_0 _S953;
    (&_S953)->primal_0 = _S765;
    (&_S953)->differential_0 = 0.0f;
    _d_max_0(&_S952, &_S953, _S950.differential_0);
    DiffPair_float_0 _S954;
    (&_S954)->primal_0 = _S797;
    (&_S954)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S954, -0.0f);
    float _S955 = - (-1.0f * - (_S954.differential_0 / _S798));
    float2  _S956 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S957;
    (&_S957)->primal_0 = e2_1;
    (&_S957)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S957, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S958;
    (&_S958)->primal_0 = e1_5;
    (&_S958)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S958, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S959;
    (&_S959)->primal_0 = e0_5;
    (&_S959)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S959, 0.0f);
    DiffPair_float_0 _S960;
    (&_S960)->primal_0 = _S795;
    (&_S960)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S960, -0.0f);
    float _S961 = - _S960.differential_0;
    float2  _S962 = _S958.differential_0 + make_float2 (_S793 * _S961, _S791 * _S960.differential_0);
    float2  _S963 = - _S962;
    float2  _S964 = _S959.differential_0 + make_float2 (_S792 * _S960.differential_0, _S794 * _S961);
    float2  _S965 = - _S964;
    float2  _S966 = - _S957.differential_0 + _S962;
    float _S967 = fx_5 * (_S947.differential_0 + _S951.differential_0 + _S966.x);
    float2  _S968 = make_float2 (_S967, fy_5 * (_S939.differential_0 + _S943.differential_0 + _S966.y)) + make_float2 (dist_coeffs_5[int(8)] * _S967, dist_coeffs_5[int(9)] * _S967);
    float2  _S969 = _S779 * _S968;
    float _S970 = dist_coeffs_5[int(4)] * _S968.y;
    float _S971 = dist_coeffs_5[int(5)] * _S968.x;
    float _S972 = _S969.x + _S969.y;
    float _S973 = r2_20 * _S972;
    float _S974 = r2_20 * _S973;
    float _S975 = dist_coeffs_5[int(7)] * _S968.y + _S970 + dist_coeffs_5[int(6)] * _S968.x + _S971 + _S782 * _S972 + _S781 * _S973 + _S780 * _S974 + dist_coeffs_5[int(3)] * (r2_20 * _S974);
    float _S976 = v_20 * _S975;
    float _S977 = u_20 * _S975;
    float2  _S978 = make_float2 (radial_11) * _S968 + make_float2 (_S736 * (v_20 * _S968.y) + _S784 * _S971 + 2.0f * (u_20 * _S971) + _S733 * (v_20 * _S968.x) + _S977 + _S977, _S786 * _S970 + 2.0f * (v_20 * _S970) + _S785 * _S968.y + _S783 * _S968.x + _S976 + _S976);
    FixedArray<float3 , 16>  _S979;
    _S979[int(0)] = _S701;
    _S979[int(1)] = _S701;
    _S979[int(2)] = _S701;
    _S979[int(3)] = _S701;
    _S979[int(4)] = _S701;
    _S979[int(5)] = _S701;
    _S979[int(6)] = _S701;
    _S979[int(7)] = _S701;
    _S979[int(8)] = _S701;
    _S979[int(9)] = _S701;
    _S979[int(10)] = _S701;
    _S979[int(11)] = _S701;
    _S979[int(12)] = _S701;
    _S979[int(13)] = _S701;
    _S979[int(14)] = _S701;
    _S979[int(15)] = _S701;
    _S979[int(0)] = _S924;
    _S979[int(1)] = _S911;
    _S979[int(2)] = _S909;
    _S979[int(3)] = _S907;
    _S979[int(4)] = _S896;
    _S979[int(5)] = _S894;
    _S979[int(6)] = _S892;
    _S979[int(7)] = _S890;
    _S979[int(8)] = _S888;
    _S979[int(9)] = _S880;
    _S979[int(10)] = _S878;
    _S979[int(11)] = _S876;
    _S979[int(12)] = _S874;
    _S979[int(13)] = _S872;
    _S979[int(14)] = _S870;
    _S979[int(15)] = _S868;
    float3  _S980 = _S979[int(0)];
    float3  _S981 = _S979[int(1)];
    float3  _S982 = _S979[int(2)];
    float3  _S983 = _S979[int(3)];
    float3  _S984 = _S979[int(4)];
    float3  _S985 = _S979[int(5)];
    float3  _S986 = _S979[int(6)];
    float3  _S987 = _S979[int(7)];
    float3  _S988 = _S979[int(8)];
    float3  _S989 = _S979[int(9)];
    float3  _S990 = _S979[int(10)];
    float3  _S991 = _S979[int(11)];
    float3  _S992 = _S979[int(12)];
    float3  _S993 = _S979[int(13)];
    float3  _S994 = _S979[int(14)];
    float3  _S995 = _S979[int(15)];
    float3  _S996 = _S852.differential_0 + _S937;
    float3  _S997 = _S854 + _S855 + _S937;
    float3  _S998 = _S853.differential_0 + _S937;
    FixedArray<float3 , 2>  _S999;
    _S999[int(0)] = _S701;
    _S999[int(1)] = _S701;
    _S999[int(1)] = _S861;
    _S999[int(0)] = _S866;
    float3  _S1000 = _S999[int(0)];
    float3  _S1001 = _S999[int(1)];
    float _S1002 = _S941.differential_0 + _S945.differential_0;
    float _S1003 = _S940.differential_0 + _S944.differential_0;
    float _S1004 = _S948.differential_0 + _S952.differential_0;
    float _S1005 = _S949.differential_0 + _S953.differential_0;
    float2  _S1006 = _S957.differential_0 + _S965;
    float2  _S1007 = _S963 + _S964;
    float2  _S1008 = make_float2 (0.0f, _S955);
    float2  _S1009 = _S767 * _S978;
    float2  _S1010 = _S778 * _S978;
    float _S1011 = _S1009.x + _S1009.y;
    if(_S771)
    {
        float _S1012 = _S1011 / _S772;
        float _S1013 = _S773 * - _S1012;
        float _S1014 = _S770 * (0.3333333432674408f * - (_S769 * _S1012));
        k_2 = _S1014 + _S1014;
        _S772 = _S1013;
        _S773 = 0.0f;
    }
    else
    {
        float _S1015 = _S1011 / _S774;
        float _S1016 = _S770 * - _S1015;
        k_2 = _S768 * _S1015;
        _S772 = 0.0f;
        _S773 = _S1016;
    }
    DiffPair_float_0 _S1017;
    (&_S1017)->primal_0 = _S768;
    (&_S1017)->differential_0 = 0.0f;
    DiffPair_float_0 _S1018;
    (&_S1018)->primal_0 = _S769;
    (&_S1018)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1017, &_S1018, k_2);
    float _S1019 = _S1018.differential_0 + _S772;
    float _S1020 = _S1017.differential_0 + _S773;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1021;
    (&_S1021)->primal_0 = _S767;
    (&_S1021)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S1021, _S1020);
    float2  _S1022 = _S1021.differential_0 + _S1010;
    float _S1023 = fx_5 * (_S1007.x + _S1005);
    float2  _S1024 = make_float2 (_S1023, fy_5 * (_S1007.y + _S1002)) + make_float2 (dist_coeffs_5[int(8)] * _S1023, dist_coeffs_5[int(9)] * _S1023);
    float2  _S1025 = _S755 * _S1024;
    float _S1026 = dist_coeffs_5[int(4)] * _S1024.y;
    float _S1027 = dist_coeffs_5[int(5)] * _S1024.x;
    float _S1028 = _S1025.x + _S1025.y;
    float _S1029 = r2_19 * _S1028;
    float _S1030 = r2_19 * _S1029;
    float _S1031 = dist_coeffs_5[int(7)] * _S1024.y + _S1026 + dist_coeffs_5[int(6)] * _S1024.x + _S1027 + _S758 * _S1028 + _S757 * _S1029 + _S756 * _S1030 + dist_coeffs_5[int(3)] * (r2_19 * _S1030);
    float _S1032 = v_19 * _S1031;
    float _S1033 = u_19 * _S1031;
    float2  _S1034 = make_float2 (radial_10) * _S1024 + make_float2 (_S736 * (v_19 * _S1024.y) + _S760 * _S1027 + 2.0f * (u_19 * _S1027) + _S733 * (v_19 * _S1024.x) + _S1033 + _S1033, _S762 * _S1026 + 2.0f * (v_19 * _S1026) + _S761 * _S1024.y + _S759 * _S1024.x + _S1032 + _S1032);
    float3  _S1035 = _S998 + make_float3 (_S1022.x, _S1022.y, _S1019);
    float2  _S1036 = _S743 * _S1034;
    float2  _S1037 = _S754 * _S1034;
    float _S1038 = _S1036.x + _S1036.y;
    if(_S747)
    {
        float _S1039 = _S1038 / _S748;
        float _S1040 = _S749 * - _S1039;
        float _S1041 = _S746 * (0.3333333432674408f * - (_S745 * _S1039));
        k_2 = _S1041 + _S1041;
        _S748 = _S1040;
        _S749 = 0.0f;
    }
    else
    {
        float _S1042 = _S1038 / _S750;
        float _S1043 = _S746 * - _S1042;
        k_2 = _S744 * _S1042;
        _S748 = 0.0f;
        _S749 = _S1043;
    }
    DiffPair_float_0 _S1044;
    (&_S1044)->primal_0 = _S744;
    (&_S1044)->differential_0 = 0.0f;
    DiffPair_float_0 _S1045;
    (&_S1045)->primal_0 = _S745;
    (&_S1045)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1044, &_S1045, k_2);
    float _S1046 = _S1045.differential_0 + _S748;
    float _S1047 = _S1044.differential_0 + _S749;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1048;
    (&_S1048)->primal_0 = _S743;
    (&_S1048)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S1048, _S1047);
    float2  _S1049 = _S1048.differential_0 + _S1037;
    float _S1050 = fx_5 * (_S1006.x + _S1004);
    float2  _S1051 = make_float2 (_S1050, fy_5 * (_S1006.y + _S1003)) + make_float2 (dist_coeffs_5[int(8)] * _S1050, dist_coeffs_5[int(9)] * _S1050);
    float2  _S1052 = _S729 * _S1051;
    float _S1053 = dist_coeffs_5[int(4)] * _S1051.y;
    float _S1054 = dist_coeffs_5[int(5)] * _S1051.x;
    float _S1055 = _S1052.x + _S1052.y;
    float _S1056 = r2_18 * _S1055;
    float _S1057 = r2_18 * _S1056;
    float _S1058 = dist_coeffs_5[int(7)] * _S1051.y + _S1053 + dist_coeffs_5[int(6)] * _S1051.x + _S1054 + _S732 * _S1055 + _S731 * _S1056 + _S730 * _S1057 + dist_coeffs_5[int(3)] * (r2_18 * _S1057);
    float _S1059 = v_18 * _S1058;
    float _S1060 = u_18 * _S1058;
    float2  _S1061 = make_float2 (radial_9) * _S1051 + make_float2 (_S736 * (v_18 * _S1051.y) + _S735 * _S1054 + 2.0f * (u_18 * _S1054) + _S733 * (v_18 * _S1051.x) + _S1060 + _S1060, _S738 * _S1053 + 2.0f * (v_18 * _S1053) + _S737 * _S1051.y + _S734 * _S1051.x + _S1059 + _S1059);
    float3  _S1062 = _S996 + make_float3 (_S1049.x, _S1049.y, _S1046);
    float2  _S1063 = _S717 * _S1061;
    float2  _S1064 = _S728 * _S1061;
    float _S1065 = _S1063.x + _S1063.y;
    if(_S721)
    {
        float _S1066 = _S1065 / _S722;
        float _S1067 = _S723 * - _S1066;
        float _S1068 = _S720 * (0.3333333432674408f * - (_S719 * _S1066));
        k_2 = _S1068 + _S1068;
        _S722 = _S1067;
        _S723 = 0.0f;
    }
    else
    {
        float _S1069 = _S1065 / _S724;
        float _S1070 = _S720 * - _S1069;
        k_2 = _S718 * _S1069;
        _S722 = 0.0f;
        _S723 = _S1070;
    }
    DiffPair_float_0 _S1071;
    (&_S1071)->primal_0 = _S718;
    (&_S1071)->differential_0 = 0.0f;
    DiffPair_float_0 _S1072;
    (&_S1072)->primal_0 = _S719;
    (&_S1072)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1071, &_S1072, k_2);
    float _S1073 = _S1072.differential_0 + _S722;
    float _S1074 = _S1071.differential_0 + _S723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1075;
    (&_S1075)->primal_0 = _S717;
    (&_S1075)->differential_0 = _S956;
    s_bwd_length_impl_1(&_S1075, _S1074);
    float2  _S1076 = _S1075.differential_0 + _S1064;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1077;
    (&_S1077)->primal_0 = vert2_c_5;
    (&_S1077)->differential_0 = _S701;
    s_bwd_length_impl_0(&_S1077, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1078;
    (&_S1078)->primal_0 = vert1_c_5;
    (&_S1078)->differential_0 = _S701;
    s_bwd_length_impl_0(&_S1078, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1079;
    (&_S1079)->primal_0 = vert0_c_5;
    (&_S1079)->differential_0 = _S701;
    s_bwd_length_impl_0(&_S1079, 0.0f);
    float3  _S1080 = _S1077.differential_0 + _S1035;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1081;
    (&_S1081)->primal_0 = R_5;
    (&_S1081)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1082;
    (&_S1082)->primal_0 = vert2_5;
    (&_S1082)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1081, &_S1082, _S1080);
    float3  _S1083 = _S1078.differential_0 + _S1062;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1084;
    (&_S1084)->primal_0 = R_5;
    (&_S1084)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1085;
    (&_S1085)->primal_0 = vert1_5;
    (&_S1085)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1084, &_S1085, _S1083);
    float3  _S1086 = _S1079.differential_0 + _S997 + make_float3 (_S1076.x, _S1076.y, _S1073);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1087;
    (&_S1087)->primal_0 = R_5;
    (&_S1087)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1088;
    (&_S1088)->primal_0 = vert0_5;
    (&_S1088)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1087, &_S1088, _S1086);
    float3  _S1089 = _S1082.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1090;
    (&_S1090)->primal_0 = _S711;
    (&_S1090)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1091;
    (&_S1091)->primal_0 = _S716;
    (&_S1091)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1090, &_S1091, _S1089);
    float _S1092 = - _S1091.differential_0.y;
    float _S1093 = _S715 * _S1091.differential_0.x;
    float _S1094 = - (_S706 * _S1091.differential_0.x);
    float3  _S1095 = _S1085.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1096;
    (&_S1096)->primal_0 = _S711;
    (&_S1096)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1097;
    (&_S1097)->primal_0 = _S714;
    (&_S1097)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1096, &_S1097, _S1095);
    float _S1098 = _S706 * _S1097.differential_0.x;
    float _S1099 = _S713 * _S1097.differential_0.x;
    float3  _S1100 = _S1088.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1101;
    (&_S1101)->primal_0 = _S711;
    (&_S1101)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1102;
    (&_S1102)->primal_0 = _S712;
    (&_S1102)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1101, &_S1102, _S1100);
    Matrix<float, 3, 3>  _S1103 = transpose_1(_S1090.differential_0 + _S1096.differential_0 + _S1101.differential_0);
    float _S1104 = 2.0f * - _S1103.rows[int(2)].z;
    float _S1105 = 2.0f * _S1103.rows[int(2)].y;
    float _S1106 = 2.0f * _S1103.rows[int(2)].x;
    float _S1107 = 2.0f * _S1103.rows[int(1)].z;
    float _S1108 = 2.0f * - _S1103.rows[int(1)].y;
    float _S1109 = 2.0f * _S1103.rows[int(1)].x;
    float _S1110 = 2.0f * _S1103.rows[int(0)].z;
    float _S1111 = 2.0f * _S1103.rows[int(0)].y;
    float _S1112 = 2.0f * - _S1103.rows[int(0)].x;
    float _S1113 = - _S1109 + _S1111;
    float _S1114 = _S1106 + - _S1110;
    float _S1115 = - _S1105 + _S1107;
    float _S1116 = _S1105 + _S1107;
    float _S1117 = _S1106 + _S1110;
    float _S1118 = _S1109 + _S1111;
    float _S1119 = _S709.w * (_S1108 + _S1112);
    float _S1120 = _S709.z * (_S1104 + _S1112);
    float _S1121 = _S709.y * (_S1104 + _S1108);
    float _S1122 = _S709.x * _S1113 + _S709.z * _S1116 + _S709.y * _S1117 + _S1119 + _S1119;
    float _S1123 = _S709.x * _S1114 + _S709.w * _S1116 + _S709.y * _S1118 + _S1120 + _S1120;
    float _S1124 = _S709.x * _S1115 + _S709.w * _S1117 + _S709.z * _S1118 + _S1121 + _S1121;
    float _S1125 = _S709.w * _S1113 + _S709.z * _S1114 + _S709.y * _S1115;
    float4  _S1126 = make_float4 (0.0f);
    float4  _S1127 = _S1126;
    *&((&_S1127)->w) = _S1122;
    *&((&_S1127)->z) = _S1123;
    *&((&_S1127)->y) = _S1124;
    *&((&_S1127)->x) = _S1125;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1128;
    (&_S1128)->primal_0 = quat_5;
    (&_S1128)->differential_0 = _S1126;
    s_bwd_normalize_impl_1(&_S1128, _S1127);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1129 = _S1128;
    float _S1130 = _S1094 + _S1098;
    float _S1131 = 0.5f * - _S1130;
    float _S1132 = _S1092 + _S1097.differential_0.y;
    DiffPair_float_0 _S1133;
    (&_S1133)->primal_0 = _S707;
    (&_S1133)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1133, _S1132);
    float _S1134 = _S1131 + _S1133.differential_0;
    float _S1135 = _S1093 + _S1099 + _S1102.differential_0.x;
    DiffPair_float_0 _S1136;
    (&_S1136)->primal_0 = _S705;
    (&_S1136)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1136, _S1135);
    float _S1137 = _S1131 + _S1136.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1138;
    (&_S1138)->primal_0 = mean_c_5;
    (&_S1138)->differential_0 = _S701;
    s_bwd_length_impl_0(&_S1138, 0.0f);
    float3  _S1139 = _S1138.differential_0 + _S849.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1140;
    (&_S1140)->primal_0 = R_5;
    (&_S1140)->differential_0 = _S927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1141;
    (&_S1141)->primal_0 = mean_5;
    (&_S1141)->differential_0 = _S701;
    s_bwd_prop_mul_0(&_S1140, &_S1141, _S1139);
    float3  _S1142 = _S1080 + _S1083 + _S1086 + _S1139 + _S930.differential_0;
    Matrix<float, 3, 3>  _S1143 = _S1081.differential_0 + _S1084.differential_0 + _S1087.differential_0 + _S1140.differential_0 + _S931;
    float3  _S1144 = make_float3 (_S1137, _S1134, _S1130);
    float3  _S1145 = _S1089 + _S1095 + _S1100 + _S1141.differential_0 + _S925;
    *v_mean_1 = _S1145;
    *v_quat_1 = _S1129.differential_0;
    *v_scale_1 = _S1144;
    *v_hardness_1 = _S1008;
    (*v_sh_coeffs_1)[int(0)] = _S980;
    (*v_sh_coeffs_1)[int(1)] = _S981;
    (*v_sh_coeffs_1)[int(2)] = _S982;
    (*v_sh_coeffs_1)[int(3)] = _S983;
    (*v_sh_coeffs_1)[int(4)] = _S984;
    (*v_sh_coeffs_1)[int(5)] = _S985;
    (*v_sh_coeffs_1)[int(6)] = _S986;
    (*v_sh_coeffs_1)[int(7)] = _S987;
    (*v_sh_coeffs_1)[int(8)] = _S988;
    (*v_sh_coeffs_1)[int(9)] = _S989;
    (*v_sh_coeffs_1)[int(10)] = _S990;
    (*v_sh_coeffs_1)[int(11)] = _S991;
    (*v_sh_coeffs_1)[int(12)] = _S992;
    (*v_sh_coeffs_1)[int(13)] = _S993;
    (*v_sh_coeffs_1)[int(14)] = _S994;
    (*v_sh_coeffs_1)[int(15)] = _S995;
    (*v_ch_coeffs_1)[int(0)] = _S1000;
    (*v_ch_coeffs_1)[int(1)] = _S1001;
    *v_R_1 = _S1143;
    *v_t_1 = _S1142;
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
    bool _S1146;
    if((*u_21) >= 0.0f)
    {
        _S1146 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1146 = false;
    }
    if(_S1146)
    {
        _S1146 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1146 = false;
    }
    if(_S1146)
    {
        _S1146 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1146 = false;
    }
    return _S1146;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_14, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1147 = *dpx_14;
    bool _S1148;
    if(((*dpx_14).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1148 = ((*dpx_14).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1148 = false;
    }
    float _S1149;
    if(_S1148)
    {
        _S1149 = dOut_11;
    }
    else
    {
        _S1149 = 0.0f;
    }
    dpx_14->primal_0 = _S1147.primal_0;
    dpx_14->differential_0 = _S1149;
    DiffPair_float_0 _S1150 = *dpMin_0;
    if((_S1147.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1149 = dOut_11;
    }
    else
    {
        _S1149 = 0.0f;
    }
    dpMin_0->primal_0 = _S1150.primal_0;
    dpMin_0->differential_0 = _S1149;
    DiffPair_float_0 _S1151 = *dpMax_0;
    if(((*dpx_14).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1149 = dOut_11;
    }
    else
    {
        _S1149 = 0.0f;
    }
    dpMax_0->primal_0 = _S1151.primal_0;
    dpMax_0->differential_0 = _S1149;
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
        DiffPair_float_0 _S1152 = *dpx_15;
        float _S1153 = val_0 * (*dpy_5).primal_0 / (*dpx_15).primal_0 * dOut_12;
        dpx_15->primal_0 = (*dpx_15).primal_0;
        dpx_15->differential_0 = _S1153;
        float _S1154 = val_0 * (F32_log((_S1152.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1154;
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
    bool _S1155;
    if(u_22 >= 0.0f)
    {
        _S1155 = v_22 >= 0.0f;
    }
    else
    {
        _S1155 = false;
    }
    if(_S1155)
    {
        _S1155 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1155 = false;
    }
    if(_S1155)
    {
        _S1155 = t_7 >= 0.0f;
    }
    else
    {
        _S1155 = false;
    }
    if(!_S1155)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1156;
    if(opac_0 < 0.0f)
    {
        _S1156 = 0.0f;
    }
    else
    {
        _S1156 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1156;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1157, float _S1158, float _S1159)
{
    return clamp_0(_S1157, _S1158, _S1159);
}

inline __device__ float s_primal_ctx_pow_0(float _S1160, float _S1161)
{
    return (F32_pow((_S1160), (_S1161)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1162, DiffPair_float_0 * _S1163, float _S1164)
{
    _d_pow_0(_S1162, _S1163, _S1164);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1165, DiffPair_float_0 * _S1166, DiffPair_float_0 * _S1167, float _S1168)
{
    _d_clamp_0(_S1165, _S1166, _S1167, _S1168);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_5)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1169 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1170 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1171 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1172 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1173 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1171);
    float d_2 = 1.0f / _S1173;
    float _S1174 = _S1173 * _S1173;
    float3  _S1175 = - _S1172;
    float _S1176 = s_primal_ctx_dot_0(_S1175, v2v0_2);
    float u_23 = d_2 * _S1176;
    float _S1177 = s_primal_ctx_dot_0(_S1172, v1v0_2);
    float v_23 = d_2 * _S1177;
    float3  _S1178 = - _S1171;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1178, rov0_2);
    bool _S1179;
    if(u_23 >= 0.0f)
    {
        _S1179 = v_23 >= 0.0f;
    }
    else
    {
        _S1179 = false;
    }
    if(_S1179)
    {
        _S1179 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1179 = false;
    }
    if(_S1179)
    {
        _S1179 = t_8 >= 0.0f;
    }
    else
    {
        _S1179 = false;
    }
    bool _S1180 = !!_S1179;
    float _S1181;
    float _S1182;
    float _S1183;
    float _S1184;
    float _S1185;
    float _S1186;
    float _S1187;
    float _S1188;
    float _S1189;
    float _S1190;
    float _S1191;
    if(_S1180)
    {
        float _S1192 = (F32_min((u_23), (v_23)));
        float _S1193 = s_primal_ctx_sqrt_0(0.5f);
        float _S1194 = _S1193 * (1.0f - u_23 - v_23);
        float _S1195 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1192), (_S1194))) * _S1195;
        float _S1196 = _S1169.primal_0.y;
        float _S1197 = 1.0f - opac_1;
        float _S1198 = 1.0f - s_primal_ctx_clamp_0(_S1196, 0.0f, 0.99989998340606689f);
        float _S1199 = 1.0f / _S1198;
        float _S1200 = _S1198 * _S1198;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1197, _S1199);
        float o_1 = _S1169.primal_0.x;
        bool _S1201 = opac_1 < 0.0f;
        if(_S1201)
        {
            _S1181 = 0.0f;
        }
        else
        {
            _S1181 = o_1 * w_1;
        }
        _S1179 = _S1201;
        _S1182 = o_1;
        _S1183 = w_1;
        _S1184 = _S1197;
        _S1185 = _S1199;
        _S1186 = _S1200;
        _S1187 = _S1196;
        _S1188 = _S1195;
        _S1189 = _S1192;
        _S1190 = _S1194;
        _S1191 = _S1193;
    }
    else
    {
        _S1179 = false;
        _S1181 = 0.0f;
        _S1182 = 0.0f;
        _S1183 = 0.0f;
        _S1184 = 0.0f;
        _S1185 = 0.0f;
        _S1186 = 0.0f;
        _S1187 = 0.0f;
        _S1188 = 0.0f;
        _S1189 = 0.0f;
        _S1190 = 0.0f;
        _S1191 = 0.0f;
    }
    float2  _S1202 = make_float2 (0.0f);
    float2  _S1203;
    if(_S1180)
    {
        if(_S1179)
        {
            _S1181 = 0.0f;
            _S1182 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1204;
            (&_S1204)->primal_0 = _S1181;
            (&_S1204)->differential_0 = 0.0f;
            DiffPair_float_0 _S1205;
            (&_S1205)->primal_0 = 0.99500000476837158f;
            (&_S1205)->differential_0 = 0.0f;
            _d_min_0(&_S1204, &_S1205, _s_dOut_5);
            float _S1206 = _S1182 * _S1204.differential_0;
            _S1181 = _S1183 * _S1204.differential_0;
            _S1182 = _S1206;
        }
        float _S1207 = - _S1182;
        DiffPair_float_0 _S1208;
        (&_S1208)->primal_0 = _S1184;
        (&_S1208)->differential_0 = 0.0f;
        DiffPair_float_0 _S1209;
        (&_S1209)->primal_0 = _S1185;
        (&_S1209)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1208, &_S1209, _S1207);
        float _S1210 = - - (_S1209.differential_0 / _S1186);
        float s_diff_opac_T_0 = - _S1208.differential_0;
        DiffPair_float_0 _S1211;
        (&_S1211)->primal_0 = _S1187;
        (&_S1211)->differential_0 = 0.0f;
        DiffPair_float_0 _S1212;
        (&_S1212)->primal_0 = 0.0f;
        (&_S1212)->differential_0 = 0.0f;
        DiffPair_float_0 _S1213;
        (&_S1213)->primal_0 = 0.99989998340606689f;
        (&_S1213)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1211, &_S1212, &_S1213, _S1210);
        float _S1214 = _S1188 * s_diff_opac_T_0;
        DiffPair_float_0 _S1215;
        (&_S1215)->primal_0 = _S1189;
        (&_S1215)->differential_0 = 0.0f;
        DiffPair_float_0 _S1216;
        (&_S1216)->primal_0 = _S1190;
        (&_S1216)->differential_0 = 0.0f;
        _d_min_0(&_S1215, &_S1216, _S1214);
        float _S1217 = - (_S1191 * _S1216.differential_0);
        DiffPair_float_0 _S1218;
        (&_S1218)->primal_0 = u_23;
        (&_S1218)->differential_0 = 0.0f;
        DiffPair_float_0 _S1219;
        (&_S1219)->primal_0 = v_23;
        (&_S1219)->differential_0 = 0.0f;
        _d_min_0(&_S1218, &_S1219, _S1215.differential_0);
        float2  _S1220 = make_float2 (_S1181, _S1211.differential_0);
        float _S1221 = _S1217 + _S1219.differential_0;
        _S1181 = _S1217 + _S1218.differential_0;
        _S1182 = _S1221;
        _S1203 = _S1220;
    }
    else
    {
        _S1181 = 0.0f;
        _S1182 = 0.0f;
        _S1203 = _S1202;
    }
    float3  _S1222 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1223;
    (&_S1223)->primal_0 = _S1178;
    (&_S1223)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1224;
    (&_S1224)->primal_0 = rov0_2;
    (&_S1224)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1223, &_S1224, 0.0f);
    float3  _S1225 = - _S1223.differential_0;
    float _S1226 = d_2 * _S1182;
    float _S1227 = _S1177 * _S1182;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1228;
    (&_S1228)->primal_0 = _S1172;
    (&_S1228)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1229;
    (&_S1229)->primal_0 = v1v0_2;
    (&_S1229)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1228, &_S1229, _S1226);
    float _S1230 = d_2 * _S1181;
    float _S1231 = _S1176 * _S1181;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1232;
    (&_S1232)->primal_0 = _S1175;
    (&_S1232)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1233;
    (&_S1233)->primal_0 = v2v0_2;
    (&_S1233)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1232, &_S1233, _S1230);
    float3  _S1234 = - _S1232.differential_0;
    float _S1235 = - ((_S1227 + _S1231) / _S1174);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1236;
    (&_S1236)->primal_0 = _S1170.primal_0;
    (&_S1236)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1237;
    (&_S1237)->primal_0 = _S1171;
    (&_S1237)->differential_0 = _S1222;
    s_bwd_prop_dot_0(&_S1236, &_S1237, _S1235);
    float3  _S1238 = _S1228.differential_0 + _S1234;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239;
    (&_S1239)->primal_0 = rov0_2;
    (&_S1239)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1240;
    (&_S1240)->primal_0 = _S1170.primal_0;
    (&_S1240)->differential_0 = _S1222;
    s_bwd_prop_cross_0(&_S1239, &_S1240, _S1238);
    float3  _S1241 = _S1225 + _S1237.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1242;
    (&_S1242)->primal_0 = v1v0_2;
    (&_S1242)->differential_0 = _S1222;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1243;
    (&_S1243)->primal_0 = v2v0_2;
    (&_S1243)->differential_0 = _S1222;
    s_bwd_prop_cross_0(&_S1242, &_S1243, _S1241);
    float3  _S1244 = _S1224.differential_0 + _S1239.differential_0;
    float3  _S1245 = _S1233.differential_0 + _S1243.differential_0;
    float3  _S1246 = _S1229.differential_0 + _S1242.differential_0;
    float3  _S1247 = - _S1244 + - _S1245 + - _S1246;
    float3  _S1248 = _S1236.differential_0 + _S1240.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1248;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1244;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1203;
    FixedArray<float3 , 3>  _S1249;
    _S1249[int(0)] = _S1222;
    _S1249[int(1)] = _S1222;
    _S1249[int(2)] = _S1222;
    _S1249[int(2)] = _S1245;
    _S1249[int(0)] = _S1247;
    _S1249[int(1)] = _S1246;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1249;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1250, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1251, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1252, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1253, float _S1254)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1250, _S1251, _S1252, _S1253, _S1254);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1255 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1256 = { _S1255, _S1255, _S1255 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1256;
    float2  _S1257 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1257;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1255;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1255;
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
    float3  _S1258 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1259 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1260 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1258);
    float d_4 = 1.0f / _S1260;
    float _S1261 = _S1260 * _S1260;
    float3  _S1262 = - _S1259;
    float _S1263 = s_primal_ctx_dot_0(_S1262, v2v0_4);
    float u_25 = d_4 * _S1263;
    float _S1264 = s_primal_ctx_dot_0(_S1259, v1v0_4);
    float v_25 = d_4 * _S1264;
    float3  _S1265 = - _S1258;
    float3  _S1266 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1267 = make_float3 (v_25) * dpcolor_0;
    float3  _S1268 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1269 = make_float3 (u_25) * dpcolor_0;
    float3  _S1270 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1271 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1272 = - (_S1270.x + _S1270.y + _S1270.z);
    float _S1273 = d_4 * dpdepth_0;
    float _S1274 = s_primal_ctx_dot_0(_S1265, rov0_4) * dpdepth_0;
    float3  _S1275 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1276;
    (&_S1276)->primal_0 = _S1265;
    (&_S1276)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1277;
    (&_S1277)->primal_0 = rov0_4;
    (&_S1277)->differential_0 = _S1275;
    s_bwd_prop_dot_0(&_S1276, &_S1277, _S1273);
    float3  _S1278 = - _S1276.differential_0;
    float _S1279 = _S1272 + _S1266.x + _S1266.y + _S1266.z;
    float _S1280 = d_4 * _S1279;
    float _S1281 = _S1264 * _S1279;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1282;
    (&_S1282)->primal_0 = _S1259;
    (&_S1282)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1283;
    (&_S1283)->primal_0 = v1v0_4;
    (&_S1283)->differential_0 = _S1275;
    s_bwd_prop_dot_0(&_S1282, &_S1283, _S1280);
    float _S1284 = _S1272 + _S1268.x + _S1268.y + _S1268.z;
    float _S1285 = d_4 * _S1284;
    float _S1286 = _S1263 * _S1284;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1287;
    (&_S1287)->primal_0 = _S1262;
    (&_S1287)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1288;
    (&_S1288)->primal_0 = v2v0_4;
    (&_S1288)->differential_0 = _S1275;
    s_bwd_prop_dot_0(&_S1287, &_S1288, _S1285);
    float3  _S1289 = - _S1287.differential_0;
    float _S1290 = - ((_S1274 + _S1281 + _S1286) / _S1261);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1291;
    (&_S1291)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1291)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1292;
    (&_S1292)->primal_0 = _S1258;
    (&_S1292)->differential_0 = _S1275;
    s_bwd_prop_dot_0(&_S1291, &_S1292, _S1290);
    float3  _S1293 = _S1282.differential_0 + _S1289;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1294;
    (&_S1294)->primal_0 = rov0_4;
    (&_S1294)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1295;
    (&_S1295)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1295)->differential_0 = _S1275;
    s_bwd_prop_cross_0(&_S1294, &_S1295, _S1293);
    float3  _S1296 = _S1278 + _S1292.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1297;
    (&_S1297)->primal_0 = v1v0_4;
    (&_S1297)->differential_0 = _S1275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1298;
    (&_S1298)->primal_0 = v2v0_4;
    (&_S1298)->differential_0 = _S1275;
    s_bwd_prop_cross_0(&_S1297, &_S1298, _S1296);
    float3  _S1299 = _S1277.differential_0 + _S1294.differential_0;
    float3  _S1300 = _S1288.differential_0 + _S1298.differential_0;
    float3  _S1301 = _S1283.differential_0 + _S1297.differential_0;
    float3  _S1302 = - _S1299 + - _S1300 + - _S1301;
    float3  _S1303 = _S1291.differential_0 + _S1295.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1303;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1299;
    FixedArray<float3 , 3>  _S1304;
    _S1304[int(0)] = _S1275;
    _S1304[int(1)] = _S1275;
    _S1304[int(2)] = _S1275;
    _S1304[int(2)] = _S1267;
    _S1304[int(1)] = _S1269;
    _S1304[int(0)] = _S1271;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1304;
    FixedArray<float3 , 3>  _S1305;
    _S1305[int(0)] = _S1275;
    _S1305[int(1)] = _S1275;
    _S1305[int(2)] = _S1275;
    _S1305[int(2)] = _S1300;
    _S1305[int(0)] = _S1302;
    _S1305[int(1)] = _S1301;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1305;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1306, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1307, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1308, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1309, float3  _S1310, float _S1311)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1306, _S1307, _S1308, _S1309, _S1310, _S1311);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1312 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1313 = { _S1312, _S1312, _S1312 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1313;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1313;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1312;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1312;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

