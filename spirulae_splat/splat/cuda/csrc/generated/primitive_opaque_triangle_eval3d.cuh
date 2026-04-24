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
        float4  _S35 = normalize_0(quat_0);
        float x_11 = _S35.y;
        float x2_0 = x_11 * x_11;
        float y2_0 = _S35.z * _S35.z;
        float z2_0 = _S35.w * _S35.w;
        float xy_0 = _S35.y * _S35.z;
        float xz_0 = _S35.y * _S35.w;
        float yz_0 = _S35.z * _S35.w;
        float wx_0 = _S35.x * _S35.y;
        float wy_0 = _S35.x * _S35.z;
        float wz_0 = _S35.x * _S35.w;
        Matrix<float, 3, 3>  _S36 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0)));
        float3  vert0_0 = mul_0(_S36, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
        float3  vert1_0 = mul_0(_S36, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
        float3  vert2_0 = mul_0(_S36, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
        float3  vert0_c_0 = mul_0(R_0, vert0_0) + t_0;
        float3  vert1_c_0 = mul_0(R_0, vert1_0) + t_0;
        float3  vert2_c_0 = mul_0(R_0, vert2_0) + t_0;
        float _S37 = vert0_c_0.z;
        float _S38 = vert1_c_0.z;
        float _S39 = vert2_c_0.z;
        if(_S37 < near_plane_0)
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
            _S32 = true;
        }
        else
        {
            _S32 = _S39 < near_plane_0;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S39 > far_plane_0;
        }
        if(_S32)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float2  uv0_0;
        for(;;)
        {
            float2  uv0_1 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S37);
            if(_S37 < 0.0f)
            {
                _S32 = true;
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
                _S32 = !((F32_min((determinant_0(_S53)), ((F32_min((_S53.rows[int(0)].x), (_S53.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
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
        bool all_valid_0 = true & (!_S32);
        for(;;)
        {
            float2  uv1_1 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S38);
            if(_S38 < 0.0f)
            {
                _S32 = true;
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
                _S32 = !((F32_min((determinant_0(_S69)), ((F32_min((_S69.rows[int(0)].x), (_S69.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
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
        bool all_valid_1 = all_valid_0 & (!_S32);
        for(;;)
        {
            float2  uv2_1 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S39);
            if(_S39 < 0.0f)
            {
                _S32 = true;
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
                _S32 = !((F32_min((determinant_0(_S85)), ((F32_min((_S85.rows[int(0)].x), (_S85.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S32)
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
        if(!(all_valid_1 & (!_S32)))
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

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_1, float4  quat_1, float3  scale_1, float2  hardness_1, FixedArray<float3 , 16>  sh_coeffs_1, FixedArray<float3 , 2>  ch_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, float4  * aabb_xyxy_1, float * depth_1, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_1)
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
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S124 = scale_1.x;
        float sx_1 = (F32_exp((_S124)));
        float _S125 = scale_1.y;
        float sy_1 = (F32_exp((_S125)));
        float sz_1 = scale_1.z - 0.5f * (_S124 + _S125);
        float4  _S126 = normalize_0(quat_1);
        float x_13 = _S126.y;
        float x2_1 = x_13 * x_13;
        float y2_1 = _S126.z * _S126.z;
        float z2_2 = _S126.w * _S126.w;
        float xy_1 = _S126.y * _S126.z;
        float xz_1 = _S126.y * _S126.w;
        float yz_1 = _S126.z * _S126.w;
        float wx_1 = _S126.x * _S126.y;
        float wy_1 = _S126.x * _S126.z;
        float wz_1 = _S126.x * _S126.w;
        Matrix<float, 3, 3>  _S127 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_2), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_2), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1)));
        float3  vert0_1 = mul_0(_S127, make_float3 (sx_1, 0.0f, 0.0f)) + mean_1;
        float3  vert1_1 = mul_0(_S127, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_1;
        float3  vert2_1 = mul_0(_S127, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_1;
        float3  vert0_c_1 = mul_0(R_1, vert0_1) + t_1;
        float3  vert1_c_1 = mul_0(R_1, vert1_1) + t_1;
        float3  vert2_c_1 = mul_0(R_1, vert2_1) + t_1;
        float _S128 = length_1(vert0_c_1);
        float _S129 = length_1(vert1_c_1);
        float _S130 = length_1(vert2_c_1);
        if(_S128 < near_plane_1)
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
            _S123 = true;
        }
        else
        {
            _S123 = _S130 < near_plane_1;
        }
        if(_S123)
        {
            _S123 = true;
        }
        else
        {
            _S123 = _S130 > far_plane_1;
        }
        if(_S123)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float2  uv0_2;
        float k_0;
        for(;;)
        {
            float2  _S131 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_2 = length_0(_S131);
            float _S132 = vert0_c_1.z;
            float theta_0 = (F32_atan2((r_2), (_S132)));
            if(theta_0 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S132;
            }
            else
            {
                k_0 = theta_0 / r_2;
            }
            float2  uv0_3 = _S131 * make_float2 (k_0);
            float2  _S133 = make_float2 (1.0f, 0.0f);
            _S103 = _S133;
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
            float _S134 = 0.0f * v_6;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float s_diff_r2_6 = u_6 + u_6 + (_S134 + _S134);
            float _S135 = dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)];
            float _S136 = dist_coeffs_1[int(1)] + r2_6 * _S135;
            float _S137 = dist_coeffs_1[int(0)] + r2_6 * _S136;
            float _S138 = s_diff_r2_6 * _S137 + (s_diff_r2_6 * _S136 + (s_diff_r2_6 * _S135 + s_diff_r2_6 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float radial_3 = 1.0f + r2_6 * _S137;
            float _S139 = 2.0f * dist_coeffs_1[int(4)];
            _S114 = _S139;
            float _S140 = _S139 * u_6;
            float _S141 = 2.0f * u_6;
            float s_diff_du_0 = _S139 * v_6 + 0.0f * _S140 + (s_diff_r2_6 + (_S141 + _S141)) * dist_coeffs_1[int(5)] + s_diff_r2_6 * dist_coeffs_1[int(6)];
            float _S142 = 2.0f * dist_coeffs_1[int(5)];
            _S115 = _S142;
            float _S143 = _S142 * u_6;
            float _S144 = 2.0f * v_6;
            float2  _S145 = _S133 * make_float2 (radial_3) + make_float2 (_S138) * uv0_3 + make_float2 (s_diff_du_0, _S142 * v_6 + 0.0f * _S143 + (s_diff_r2_6 + (_S134 + 0.0f * _S144)) * dist_coeffs_1[int(4)] + s_diff_r2_6 * dist_coeffs_1[int(7)]);
            float2  _S146 = _S145 + make_float2 (_S145.x * dist_coeffs_1[int(8)] + _S145.y * dist_coeffs_1[int(9)], 0.0f);
            float2  _S147 = make_float2 (0.0f, 1.0f);
            _S116 = _S147;
            float _S148 = 0.0f * u_6;
            float s_diff_r2_7 = _S148 + _S148 + (v_6 + v_6);
            float _S149 = s_diff_r2_7 * _S137 + (s_diff_r2_7 * _S136 + (s_diff_r2_7 * _S135 + s_diff_r2_7 * dist_coeffs_1[int(3)] * r2_6) * r2_6) * r2_6;
            float _S150 = 0.0f * _S139;
            _S117 = _S150;
            float s_diff_du_1 = _S150 * v_6 + _S140 + (s_diff_r2_7 + (_S148 + 0.0f * _S141)) * dist_coeffs_1[int(5)] + s_diff_r2_7 * dist_coeffs_1[int(6)];
            float _S151 = 0.0f * _S142;
            _S118 = _S151;
            float2  _S152 = _S147 * make_float2 (radial_3) + make_float2 (_S149) * uv0_3 + make_float2 (s_diff_du_1, _S151 * v_6 + _S143 + (s_diff_r2_7 + (_S144 + _S144)) * dist_coeffs_1[int(4)] + s_diff_r2_7 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S153 = transpose_0(makeMatrix<float, 2, 2> (_S146, _S152 + make_float2 (_S152.x * dist_coeffs_1[int(8)] + _S152.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S154 = !((F32_min((determinant_0(_S153)), ((F32_min((_S153.rows[int(0)].x), (_S153.rows[int(1)].y)))))) > 0.0f);
            _S119 = _S154;
            if(_S154)
            {
                uv0_2 = uv0_3;
                break;
            }
            float2  _S155 = uv0_3 * make_float2 (radial_3) + make_float2 (_S140 * v_6 + dist_coeffs_1[int(5)] * (r2_6 + _S141 * u_6) + dist_coeffs_1[int(6)] * r2_6, _S143 * v_6 + dist_coeffs_1[int(4)] * (r2_6 + _S144 * v_6) + dist_coeffs_1[int(7)] * r2_6);
            float2  _S156 = _S155 + make_float2 (dist_coeffs_1[int(8)] * _S155.x + dist_coeffs_1[int(9)] * _S155.y, 0.0f);
            uv0_2 = make_float2 (fx_1 * _S156.x + cx_1, fy_1 * _S156.y + cy_1);
            break;
        }
        float2  uv1_2;
        bool all_valid_2 = true & (!_S119);
        for(;;)
        {
            float2  _S157 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_3 = length_0(_S157);
            float _S158 = vert1_c_1.z;
            float theta_1 = (F32_atan2((r_3), (_S158)));
            if(theta_1 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S158;
            }
            else
            {
                k_0 = theta_1 / r_3;
            }
            float2  uv1_3 = _S157 * make_float2 (k_0);
            float u_7 = uv1_3.x;
            float v_7 = uv1_3.y;
            float _S159 = 0.0f * v_7;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float s_diff_r2_8 = u_7 + u_7 + (_S159 + _S159);
            float _S160 = _S106 + r2_7 * _S107;
            float _S161 = _S105 + r2_7 * _S160;
            float _S162 = _S104 + r2_7 * _S161;
            float radial_4 = 1.0f + r2_7 * _S162;
            float _S163 = _S114 * u_7;
            float _S164 = 2.0f * u_7;
            float _S165 = _S115 * u_7;
            float _S166 = 2.0f * v_7;
            float2  _S167 = _S103 * make_float2 (radial_4) + make_float2 (s_diff_r2_8 * _S162 + (s_diff_r2_8 * _S161 + (s_diff_r2_8 * _S160 + s_diff_r2_8 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S114 * v_7 + 0.0f * _S163 + (s_diff_r2_8 + (_S164 + _S164)) * _S109 + s_diff_r2_8 * _S110, _S115 * v_7 + 0.0f * _S165 + (s_diff_r2_8 + (_S159 + 0.0f * _S166)) * _S108 + s_diff_r2_8 * _S111);
            float _S168 = 0.0f * u_7;
            float s_diff_r2_9 = _S168 + _S168 + (v_7 + v_7);
            float2  _S169 = _S116 * make_float2 (radial_4) + make_float2 (s_diff_r2_9 * _S162 + (s_diff_r2_9 * _S161 + (s_diff_r2_9 * _S160 + s_diff_r2_9 * _S107 * r2_7) * r2_7) * r2_7) * uv1_3 + make_float2 (_S117 * v_7 + _S163 + (s_diff_r2_9 + (_S168 + 0.0f * _S164)) * _S109 + s_diff_r2_9 * _S110, _S118 * v_7 + _S165 + (s_diff_r2_9 + (_S166 + _S166)) * _S108 + s_diff_r2_9 * _S111);
            Matrix<float, 2, 2>  _S170 = transpose_0(makeMatrix<float, 2, 2> (_S167 + make_float2 (_S167.x * _S112 + _S167.y * _S113, 0.0f), _S169 + make_float2 (_S169.x * _S112 + _S169.y * _S113, 0.0f)));
            bool _S171 = !((F32_min((determinant_0(_S170)), ((F32_min((_S170.rows[int(0)].x), (_S170.rows[int(1)].y)))))) > 0.0f);
            _S120 = _S171;
            if(_S171)
            {
                uv1_2 = uv1_3;
                break;
            }
            float2  _S172 = uv1_3 * make_float2 (radial_4) + make_float2 (_S163 * v_7 + _S109 * (r2_7 + _S164 * u_7) + _S110 * r2_7, _S165 * v_7 + _S108 * (r2_7 + _S166 * v_7) + _S111 * r2_7);
            float2  _S173 = _S172 + make_float2 (_S112 * _S172.x + _S113 * _S172.y, 0.0f);
            uv1_2 = make_float2 (fx_1 * _S173.x + cx_1, fy_1 * _S173.y + cy_1);
            break;
        }
        float2  uv2_2;
        bool all_valid_3 = all_valid_2 & (!_S120);
        for(;;)
        {
            float2  _S174 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_4 = length_0(_S174);
            float _S175 = vert2_c_1.z;
            float theta_2 = (F32_atan2((r_4), (_S175)));
            if(theta_2 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S175;
            }
            else
            {
                k_0 = theta_2 / r_4;
            }
            float2  uv2_3 = _S174 * make_float2 (k_0);
            float u_8 = uv2_3.x;
            float v_8 = uv2_3.y;
            float _S176 = 0.0f * v_8;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float s_diff_r2_10 = u_8 + u_8 + (_S176 + _S176);
            float _S177 = _S106 + r2_8 * _S107;
            float _S178 = _S105 + r2_8 * _S177;
            float _S179 = _S104 + r2_8 * _S178;
            float radial_5 = 1.0f + r2_8 * _S179;
            float _S180 = _S114 * u_8;
            float _S181 = 2.0f * u_8;
            float _S182 = _S115 * u_8;
            float _S183 = 2.0f * v_8;
            float2  _S184 = _S103 * make_float2 (radial_5) + make_float2 (s_diff_r2_10 * _S179 + (s_diff_r2_10 * _S178 + (s_diff_r2_10 * _S177 + s_diff_r2_10 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S114 * v_8 + 0.0f * _S180 + (s_diff_r2_10 + (_S181 + _S181)) * _S109 + s_diff_r2_10 * _S110, _S115 * v_8 + 0.0f * _S182 + (s_diff_r2_10 + (_S176 + 0.0f * _S183)) * _S108 + s_diff_r2_10 * _S111);
            float _S185 = 0.0f * u_8;
            float s_diff_r2_11 = _S185 + _S185 + (v_8 + v_8);
            float2  _S186 = _S116 * make_float2 (radial_5) + make_float2 (s_diff_r2_11 * _S179 + (s_diff_r2_11 * _S178 + (s_diff_r2_11 * _S177 + s_diff_r2_11 * _S107 * r2_8) * r2_8) * r2_8) * uv2_3 + make_float2 (_S117 * v_8 + _S180 + (s_diff_r2_11 + (_S185 + 0.0f * _S181)) * _S109 + s_diff_r2_11 * _S110, _S118 * v_8 + _S182 + (s_diff_r2_11 + (_S183 + _S183)) * _S108 + s_diff_r2_11 * _S111);
            Matrix<float, 2, 2>  _S187 = transpose_0(makeMatrix<float, 2, 2> (_S184 + make_float2 (_S184.x * _S112 + _S184.y * _S113, 0.0f), _S186 + make_float2 (_S186.x * _S112 + _S186.y * _S113, 0.0f)));
            bool _S188 = !((F32_min((determinant_0(_S187)), ((F32_min((_S187.rows[int(0)].x), (_S187.rows[int(1)].y)))))) > 0.0f);
            _S121 = _S188;
            if(_S188)
            {
                uv2_2 = uv2_3;
                break;
            }
            float2  _S189 = uv2_3 * make_float2 (radial_5) + make_float2 (_S180 * v_8 + _S109 * (r2_8 + _S181 * u_8) + _S110 * r2_8, _S182 * v_8 + _S108 * (r2_8 + _S183 * v_8) + _S111 * r2_8);
            float2  _S190 = _S189 + make_float2 (_S112 * _S189.x + _S113 * _S189.y, 0.0f);
            uv2_2 = make_float2 (fx_1 * _S190.x + cx_1, fy_1 * _S190.y + cy_1);
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
        float _S191 = uv0_2.x;
        float _S192 = uv1_2.x;
        float _S193 = uv2_2.x;
        float xmax_1 = (F32_max(((F32_max((_S191), (_S192)))), (_S193))) + offset_1;
        float xmin_1 = (F32_min(((F32_min((_S191), (_S192)))), (_S193))) - offset_1;
        float _S194 = uv0_2.y;
        float _S195 = uv1_2.y;
        float _S196 = uv2_2.y;
        float ymax_1 = (F32_max(((F32_max((_S194), (_S195)))), (_S196))) + offset_1;
        float ymin_1 = (F32_min(((F32_min((_S194), (_S195)))), (_S196))) - offset_1;
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
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (xmin_1, ymin_1, xmax_1, ymax_1);
        float3  _S197 = (vert0_c_1 + vert1_c_1 + vert2_c_1) / make_float3 (3.0f);
        float x_14 = _S197.x;
        float y_5 = _S197.y;
        float z_1 = _S197.z;
        float _S198 = x_14 * x_14 + y_5 * y_5;
        *depth_1 = z_1 * z_1 * z_1 * z_1 + 0.001953125f * _S198 * _S198;
        float3  _S199 = mean_1 - - mul_0(transpose_1(R_1), t_1);
        float _S200 = _S199.x;
        float _S201 = _S199.y;
        float _S202 = _S199.z;
        float norm_1 = (F32_sqrt((_S200 * _S200 + _S201 * _S201 + _S202 * _S202)));
        float x_15 = _S200 / norm_1;
        float y_6 = _S201 / norm_1;
        float z_2 = _S202 / norm_1;
        float z2_3 = z_2 * z_2;
        float fTmp0B_1 = -1.09254848957061768f * z_2;
        float fC1_1 = x_15 * x_15 - y_6 * y_6;
        float fS1_1 = 2.0f * x_15 * y_6;
        float fTmp0C_1 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_2;
        float3  color_1 = make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_1[int(1)] + make_float3 (z_2) * sh_coeffs_1[int(2)] - make_float3 (x_15) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_6) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_15) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_1 + y_6 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_6) * sh_coeffs_1[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_15) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_1 - y_6 * fS1_1)) * sh_coeffs_1[int(15)]);
        float3  _S203 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_1 + ch_coeffs_1[int(0)] + make_float3 (0.5f), _S203);
        float3  _S204 = color_1 - ch_coeffs_1[int(0)] * make_float3 (0.5f);
        float3  _S205 = ch_coeffs_1[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S204 + _S205 + make_float3 (0.5f), _S203);
        (*rgbs_1)[int(2)] = max_0(_S204 - _S205 + make_float3 (0.5f), _S203);
        (*verts_1)[int(0)] = vert0_1;
        (*verts_1)[int(1)] = vert1_1;
        (*verts_1)[int(2)] = vert2_1;
        float3  _S206 = normalize_1(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S206 * make_float3 (float(- (F32_sign((dot_0(_S206, mean_c_1))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_2, float4  quat_2, float3  scale_2, float2  hardness_2, FixedArray<float3 , 16>  sh_coeffs_2, FixedArray<float3 , 2>  ch_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, float4  * aabb_xyxy_2, float * depth_2, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_2)
{
    float3  mean_c_2 = mul_0(R_2, mean_2) + t_2;
    float _S207 = scale_2.x;
    float sx_2 = (F32_exp((_S207)));
    float _S208 = scale_2.y;
    float sy_2 = (F32_exp((_S208)));
    float sz_2 = scale_2.z - 0.5f * (_S207 + _S208);
    float4  _S209 = normalize_0(quat_2);
    float x_16 = _S209.y;
    float x2_2 = x_16 * x_16;
    float y2_2 = _S209.z * _S209.z;
    float z2_4 = _S209.w * _S209.w;
    float xy_2 = _S209.y * _S209.z;
    float xz_2 = _S209.y * _S209.w;
    float yz_2 = _S209.z * _S209.w;
    float wx_2 = _S209.x * _S209.y;
    float wy_2 = _S209.x * _S209.z;
    float wz_2 = _S209.x * _S209.w;
    Matrix<float, 3, 3>  _S210 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_4), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_4), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)));
    float3  vert0_2 = mul_0(_S210, make_float3 (sx_2, 0.0f, 0.0f)) + mean_2;
    float3  vert1_2 = mul_0(_S210, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_2;
    float3  vert2_2 = mul_0(_S210, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_2;
    float3  vert0_c_2 = mul_0(R_2, vert0_2) + t_2;
    float3  vert1_c_2 = mul_0(R_2, vert1_2) + t_2;
    float3  vert2_c_2 = mul_0(R_2, vert2_2) + t_2;
    float2  _S211 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_9 = _S211.x;
    float v_9 = _S211.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S212 = 2.0f * dist_coeffs_2[int(4)];
    float _S213 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S214 = _S211 * make_float2 (1.0f + r2_9 * (dist_coeffs_2[int(0)] + r2_9 * (dist_coeffs_2[int(1)] + r2_9 * (dist_coeffs_2[int(2)] + r2_9 * dist_coeffs_2[int(3)])))) + make_float2 (_S212 * u_9 * v_9 + dist_coeffs_2[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_2[int(6)] * r2_9, _S213 * u_9 * v_9 + dist_coeffs_2[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_2[int(7)] * r2_9);
    float2  _S215 = _S214 + make_float2 (dist_coeffs_2[int(8)] * _S214.x + dist_coeffs_2[int(9)] * _S214.y, 0.0f);
    float _S216 = fx_2 * _S215.x + cx_2;
    float _S217 = fy_2 * _S215.y + cy_2;
    float2  uv0_4 = make_float2 (_S216, _S217);
    float2  _S218 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_10 = _S218.x;
    float v_10 = _S218.y;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float2  _S219 = _S218 * make_float2 (1.0f + r2_10 * (dist_coeffs_2[int(0)] + r2_10 * (dist_coeffs_2[int(1)] + r2_10 * (dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)])))) + make_float2 (_S212 * u_10 * v_10 + dist_coeffs_2[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + dist_coeffs_2[int(6)] * r2_10, _S213 * u_10 * v_10 + dist_coeffs_2[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + dist_coeffs_2[int(7)] * r2_10);
    float2  _S220 = _S219 + make_float2 (dist_coeffs_2[int(8)] * _S219.x + dist_coeffs_2[int(9)] * _S219.y, 0.0f);
    float _S221 = fx_2 * _S220.x + cx_2;
    float _S222 = fy_2 * _S220.y + cy_2;
    float2  uv1_4 = make_float2 (_S221, _S222);
    float2  _S223 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_11 = _S223.x;
    float v_11 = _S223.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float2  _S224 = _S223 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S212 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S213 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
    float2  _S225 = _S224 + make_float2 (dist_coeffs_2[int(8)] * _S224.x + dist_coeffs_2[int(9)] * _S224.y, 0.0f);
    float _S226 = fx_2 * _S225.x + cx_2;
    float _S227 = fy_2 * _S225.y + cy_2;
    float2  uv2_4 = make_float2 (_S226, _S227);
    float2  e0_2 = uv1_4 - uv0_4;
    float2  e1_2 = uv2_4 - uv1_4;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(uv0_4 - uv2_4)));
    *aabb_xyxy_2 = make_float4 ((F32_min(((F32_min((_S216), (_S221)))), (_S226))) - offset_2, (F32_min(((F32_min((_S217), (_S222)))), (_S227))) - offset_2, (F32_max(((F32_max((_S216), (_S221)))), (_S226))) + offset_2, (F32_max(((F32_max((_S217), (_S222)))), (_S227))) + offset_2);
    *depth_2 = ((vert0_c_2 + vert1_c_2 + vert2_c_2) / make_float3 (3.0f)).z;
    float3  _S228 = mean_2 - - mul_0(transpose_1(R_2), t_2);
    float _S229 = _S228.x;
    float _S230 = _S228.y;
    float _S231 = _S228.z;
    float norm_2 = (F32_sqrt((_S229 * _S229 + _S230 * _S230 + _S231 * _S231)));
    float x_17 = _S229 / norm_2;
    float y_7 = _S230 / norm_2;
    float z_3 = _S231 / norm_2;
    float z2_5 = z_3 * z_3;
    float fTmp0B_2 = -1.09254848957061768f * z_3;
    float fC1_2 = x_17 * x_17 - y_7 * y_7;
    float fS1_2 = 2.0f * x_17 * y_7;
    float fTmp0C_2 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_3;
    float3  color_2 = make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * sh_coeffs_2[int(1)] + make_float3 (z_3) * sh_coeffs_2[int(2)] - make_float3 (x_17) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_7) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_17) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_17 * fS1_2 + y_7 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_7) * sh_coeffs_2[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_17) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_2 - y_7 * fS1_2)) * sh_coeffs_2[int(15)]);
    float3  _S232 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_2 + ch_coeffs_2[int(0)] + make_float3 (0.5f), _S232);
    float3  _S233 = color_2 - ch_coeffs_2[int(0)] * make_float3 (0.5f);
    float3  _S234 = ch_coeffs_2[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S233 + _S234 + make_float3 (0.5f), _S232);
    (*rgbs_2)[int(2)] = max_0(_S233 - _S234 + make_float3 (0.5f), _S232);
    (*verts_2)[int(0)] = vert0_2;
    (*verts_2)[int(1)] = vert1_2;
    (*verts_2)[int(2)] = vert2_2;
    float3  _S235 = normalize_1(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S235 * make_float3 (float(- (F32_sign((dot_0(_S235, mean_c_2))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_3, float4  quat_3, float3  scale_3, float2  hardness_3, FixedArray<float3 , 16>  sh_coeffs_3, FixedArray<float3 , 2>  ch_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, float4  * aabb_xyxy_3, float * depth_3, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_3)
{
    float3  mean_c_3 = mul_0(R_3, mean_3) + t_3;
    float _S236 = scale_3.x;
    float sx_3 = (F32_exp((_S236)));
    float _S237 = scale_3.y;
    float sy_3 = (F32_exp((_S237)));
    float sz_3 = scale_3.z - 0.5f * (_S236 + _S237);
    float4  _S238 = normalize_0(quat_3);
    float x_18 = _S238.y;
    float x2_3 = x_18 * x_18;
    float y2_3 = _S238.z * _S238.z;
    float z2_6 = _S238.w * _S238.w;
    float xy_3 = _S238.y * _S238.z;
    float xz_3 = _S238.y * _S238.w;
    float yz_3 = _S238.z * _S238.w;
    float wx_3 = _S238.x * _S238.y;
    float wy_3 = _S238.x * _S238.z;
    float wz_3 = _S238.x * _S238.w;
    Matrix<float, 3, 3>  _S239 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_6), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_6), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    float3  vert0_3 = mul_0(_S239, make_float3 (sx_3, 0.0f, 0.0f)) + mean_3;
    float3  vert1_3 = mul_0(_S239, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_3;
    float3  vert2_3 = mul_0(_S239, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_3;
    float3  vert0_c_3 = mul_0(R_3, vert0_3) + t_3;
    float3  vert1_c_3 = mul_0(R_3, vert1_3) + t_3;
    float3  vert2_c_3 = mul_0(R_3, vert2_3) + t_3;
    float2  _S240 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_5 = length_0(_S240);
    float _S241 = vert0_c_3.z;
    float theta_3 = (F32_atan2((r_5), (_S241)));
    float k_1;
    if(theta_3 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_3 * theta_3 / 3.0f) / _S241;
    }
    else
    {
        k_1 = theta_3 / r_5;
    }
    float2  _S242 = _S240 * make_float2 (k_1);
    float u_12 = _S242.x;
    float v_12 = _S242.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S243 = 2.0f * dist_coeffs_3[int(4)];
    float _S244 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S245 = _S242 * make_float2 (1.0f + r2_12 * (dist_coeffs_3[int(0)] + r2_12 * (dist_coeffs_3[int(1)] + r2_12 * (dist_coeffs_3[int(2)] + r2_12 * dist_coeffs_3[int(3)])))) + make_float2 (_S243 * u_12 * v_12 + dist_coeffs_3[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + dist_coeffs_3[int(6)] * r2_12, _S244 * u_12 * v_12 + dist_coeffs_3[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + dist_coeffs_3[int(7)] * r2_12);
    float2  _S246 = _S245 + make_float2 (dist_coeffs_3[int(8)] * _S245.x + dist_coeffs_3[int(9)] * _S245.y, 0.0f);
    float _S247 = fx_3 * _S246.x + cx_3;
    float _S248 = fy_3 * _S246.y + cy_3;
    float2  uv0_5 = make_float2 (_S247, _S248);
    float2  _S249 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_6 = length_0(_S249);
    float _S250 = vert1_c_3.z;
    float theta_4 = (F32_atan2((r_6), (_S250)));
    if(theta_4 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_4 * theta_4 / 3.0f) / _S250;
    }
    else
    {
        k_1 = theta_4 / r_6;
    }
    float2  _S251 = _S249 * make_float2 (k_1);
    float u_13 = _S251.x;
    float v_13 = _S251.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S252 = _S251 * make_float2 (1.0f + r2_13 * (dist_coeffs_3[int(0)] + r2_13 * (dist_coeffs_3[int(1)] + r2_13 * (dist_coeffs_3[int(2)] + r2_13 * dist_coeffs_3[int(3)])))) + make_float2 (_S243 * u_13 * v_13 + dist_coeffs_3[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_3[int(6)] * r2_13, _S244 * u_13 * v_13 + dist_coeffs_3[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_3[int(7)] * r2_13);
    float2  _S253 = _S252 + make_float2 (dist_coeffs_3[int(8)] * _S252.x + dist_coeffs_3[int(9)] * _S252.y, 0.0f);
    float _S254 = fx_3 * _S253.x + cx_3;
    float _S255 = fy_3 * _S253.y + cy_3;
    float2  uv1_5 = make_float2 (_S254, _S255);
    float2  _S256 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_7 = length_0(_S256);
    float _S257 = vert2_c_3.z;
    float theta_5 = (F32_atan2((r_7), (_S257)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S257;
    }
    else
    {
        k_1 = theta_5 / r_7;
    }
    float2  _S258 = _S256 * make_float2 (k_1);
    float u_14 = _S258.x;
    float v_14 = _S258.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S259 = _S258 * make_float2 (1.0f + r2_14 * (dist_coeffs_3[int(0)] + r2_14 * (dist_coeffs_3[int(1)] + r2_14 * (dist_coeffs_3[int(2)] + r2_14 * dist_coeffs_3[int(3)])))) + make_float2 (_S243 * u_14 * v_14 + dist_coeffs_3[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + dist_coeffs_3[int(6)] * r2_14, _S244 * u_14 * v_14 + dist_coeffs_3[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + dist_coeffs_3[int(7)] * r2_14);
    float2  _S260 = _S259 + make_float2 (dist_coeffs_3[int(8)] * _S259.x + dist_coeffs_3[int(9)] * _S259.y, 0.0f);
    float _S261 = fx_3 * _S260.x + cx_3;
    float _S262 = fy_3 * _S260.y + cy_3;
    float2  uv2_5 = make_float2 (_S261, _S262);
    float2  e0_3 = uv1_5 - uv0_5;
    float2  e1_3 = uv2_5 - uv1_5;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(uv0_5 - uv2_5)));
    *aabb_xyxy_3 = make_float4 ((F32_min(((F32_min((_S247), (_S254)))), (_S261))) - offset_3, (F32_min(((F32_min((_S248), (_S255)))), (_S262))) - offset_3, (F32_max(((F32_max((_S247), (_S254)))), (_S261))) + offset_3, (F32_max(((F32_max((_S248), (_S255)))), (_S262))) + offset_3);
    float3  _S263 = (vert0_c_3 + vert1_c_3 + vert2_c_3) / make_float3 (3.0f);
    float x_19 = _S263.x;
    float y_8 = _S263.y;
    float z_4 = _S263.z;
    float _S264 = x_19 * x_19 + y_8 * y_8;
    *depth_3 = z_4 * z_4 * z_4 * z_4 + 0.001953125f * _S264 * _S264;
    float3  _S265 = mean_3 - - mul_0(transpose_1(R_3), t_3);
    float _S266 = _S265.x;
    float _S267 = _S265.y;
    float _S268 = _S265.z;
    float norm_3 = (F32_sqrt((_S266 * _S266 + _S267 * _S267 + _S268 * _S268)));
    float x_20 = _S266 / norm_3;
    float y_9 = _S267 / norm_3;
    float z_5 = _S268 / norm_3;
    float z2_7 = z_5 * z_5;
    float fTmp0B_3 = -1.09254848957061768f * z_5;
    float fC1_3 = x_20 * x_20 - y_9 * y_9;
    float fS1_3 = 2.0f * x_20 * y_9;
    float fTmp0C_3 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_5;
    float3  color_3 = make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * sh_coeffs_3[int(1)] + make_float3 (z_5) * sh_coeffs_3[int(2)] - make_float3 (x_20) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_9) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_20) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_20 * fS1_3 + y_9 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_9) * sh_coeffs_3[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_20) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_3 - y_9 * fS1_3)) * sh_coeffs_3[int(15)]);
    float3  _S269 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_3 + ch_coeffs_3[int(0)] + make_float3 (0.5f), _S269);
    float3  _S270 = color_3 - ch_coeffs_3[int(0)] * make_float3 (0.5f);
    float3  _S271 = ch_coeffs_3[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S270 + _S271 + make_float3 (0.5f), _S269);
    (*rgbs_3)[int(2)] = max_0(_S270 - _S271 + make_float3 (0.5f), _S269);
    (*verts_3)[int(0)] = vert0_3;
    (*verts_3)[int(1)] = vert1_3;
    (*verts_3)[int(2)] = vert2_3;
    float3  _S272 = normalize_1(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S272 * make_float3 (float(- (F32_sign((dot_0(_S272, mean_c_3))))));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S273, float3  _S274)
{
    return mul_0(_S273, _S274);
}

inline __device__ float s_primal_ctx_exp_0(float _S275)
{
    return (F32_exp((_S275)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S276)
{
    return (F32_sqrt((_S276)));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S277, float3  _S278)
{
    return cross_0(_S277, _S278);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S279, float3  _S280)
{
    return dot_0(_S279, _S280);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S281, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S282, float _S283)
{
    _d_dot_0(_S281, _S282, _S283);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S284, float _S285)
{
    _d_sqrt_0(_S284, _S285);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S286 = (*dpx_9).primal_0.x;
    float _S287 = (*dpx_9).primal_0.y;
    float _S288 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S289;
    (&_S289)->primal_0 = _S286 * _S286 + _S287 * _S287 + _S288 * _S288;
    (&_S289)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S289, _s_dOut_0);
    float _S290 = (*dpx_9).primal_0.z * _S289.differential_0;
    float _S291 = _S290 + _S290;
    float _S292 = (*dpx_9).primal_0.y * _S289.differential_0;
    float _S293 = _S292 + _S292;
    float _S294 = (*dpx_9).primal_0.x * _S289.differential_0;
    float _S295 = _S294 + _S294;
    float3  _S296 = make_float3 (0.0f);
    *&((&_S296)->z) = _S291;
    *&((&_S296)->y) = _S293;
    *&((&_S296)->x) = _S295;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S296;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S297, float _S298)
{
    s_bwd_prop_length_impl_0(_S297, _S298);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_1)
{
    float _S299 = length_1((*dpx_10).primal_0);
    float3  _S300 = (*dpx_10).primal_0 * _s_dOut_1;
    float3  _S301 = make_float3 (1.0f / _S299) * _s_dOut_1;
    float _S302 = - ((_S300.x + _S300.y + _S300.z) / (_S299 * _S299));
    float3  _S303 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S304;
    (&_S304)->primal_0 = (*dpx_10).primal_0;
    (&_S304)->differential_0 = _S303;
    s_bwd_length_impl_0(&_S304, _S302);
    float3  _S305 = _S301 + _S304.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S305;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S306, float3  _S307)
{
    s_bwd_prop_normalize_impl_0(_S306, _S307);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S308, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S309, float3  _S310)
{
    _d_cross_0(_S308, _S309, _S310);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S311, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S312, float3  _S313)
{
    _d_max_vector_0(_S311, _S312, _S313);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S314, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S315, float3  _S316)
{
    _d_mul_0(_S314, _S315, _S316);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S317, float _S318)
{
    _d_exp2_0(_S317, _S318);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float _s_dOut_2)
{
    float _S319 = (*dpx_11).primal_0.x;
    float _S320 = (*dpx_11).primal_0.y;
    DiffPair_float_0 _S321;
    (&_S321)->primal_0 = _S319 * _S319 + _S320 * _S320;
    (&_S321)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S321, _s_dOut_2);
    float _S322 = (*dpx_11).primal_0.y * _S321.differential_0;
    float _S323 = _S322 + _S322;
    float _S324 = (*dpx_11).primal_0.x * _S321.differential_0;
    float _S325 = _S324 + _S324;
    float2  _S326 = make_float2 (0.0f);
    *&((&_S326)->y) = _S323;
    *&((&_S326)->x) = _S325;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S326;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S327, float _S328)
{
    s_bwd_prop_length_impl_1(_S327, _S328);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S329, float _S330)
{
    _d_abs_0(_S329, _S330);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_12, float _s_dOut_3)
{
    float _S331 = (*dpx_12).primal_0.x;
    float _S332 = (*dpx_12).primal_0.y;
    float _S333 = (*dpx_12).primal_0.z;
    float _S334 = (*dpx_12).primal_0.w;
    DiffPair_float_0 _S335;
    (&_S335)->primal_0 = _S331 * _S331 + _S332 * _S332 + _S333 * _S333 + _S334 * _S334;
    (&_S335)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S335, _s_dOut_3);
    float _S336 = (*dpx_12).primal_0.w * _S335.differential_0;
    float _S337 = _S336 + _S336;
    float _S338 = (*dpx_12).primal_0.z * _S335.differential_0;
    float _S339 = _S338 + _S338;
    float _S340 = (*dpx_12).primal_0.y * _S335.differential_0;
    float _S341 = _S340 + _S340;
    float _S342 = (*dpx_12).primal_0.x * _S335.differential_0;
    float _S343 = _S342 + _S342;
    float4  _S344 = make_float4 (0.0f);
    *&((&_S344)->w) = _S337;
    *&((&_S344)->z) = _S339;
    *&((&_S344)->y) = _S341;
    *&((&_S344)->x) = _S343;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S344;
    return;
}

inline __device__ void s_bwd_length_impl_2(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S345, float _S346)
{
    s_bwd_prop_length_impl_2(_S345, _S346);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_13, float4  _s_dOut_4)
{
    float _S347 = length_2((*dpx_13).primal_0);
    float4  _S348 = (*dpx_13).primal_0 * _s_dOut_4;
    float4  _S349 = make_float4 (1.0f / _S347) * _s_dOut_4;
    float _S350 = - ((_S348.x + _S348.y + _S348.z + _S348.w) / (_S347 * _S347));
    float4  _S351 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S352;
    (&_S352)->primal_0 = (*dpx_13).primal_0;
    (&_S352)->differential_0 = _S351;
    s_bwd_length_impl_2(&_S352, _S350);
    float4  _S353 = _S349 + _S352.differential_0;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S353;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S354, float4  _S355)
{
    s_bwd_prop_normalize_impl_1(_S354, _S355);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S356, float _S357)
{
    _d_exp_0(_S356, _S357);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_4, float4  quat_4, float3  scale_4, float2  hardness_4, FixedArray<float3 , 16>  sh_coeffs_4, FixedArray<float3 , 2>  ch_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float v_depth_0, FixedArray<float3 , 3>  v_verts_0, FixedArray<float3 , 3>  v_rgbs_0, float3  v_normal_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, mean_4) + t_4;
    float _S358 = scale_4.x;
    float _S359 = s_primal_ctx_exp_0(_S358);
    float _S360 = scale_4.y;
    float _S361 = s_primal_ctx_exp_0(_S360);
    float sz_4 = scale_4.z - 0.5f * (_S358 + _S360);
    float4  _S362 = normalize_0(quat_4);
    float _S363 = _S362.y;
    float x2_4 = _S363 * _S363;
    float y2_4 = _S362.z * _S362.z;
    float z2_8 = _S362.w * _S362.w;
    float xy_4 = _S362.y * _S362.z;
    float xz_4 = _S362.y * _S362.w;
    float yz_4 = _S362.z * _S362.w;
    float wx_4 = _S362.x * _S362.y;
    float wy_4 = _S362.x * _S362.z;
    float wz_4 = _S362.x * _S362.w;
    Matrix<float, 3, 3>  _S364 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_8), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_8), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    float3  _S365 = make_float3 (_S359, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S364, _S365) + mean_4;
    float _S366 = -0.5f + sz_4;
    float3  _S367 = make_float3 (_S359 * _S366, _S361, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S364, _S367) + mean_4;
    float _S368 = -0.5f - sz_4;
    float3  _S369 = make_float3 (_S359 * _S368, - _S361, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S364, _S369) + mean_4;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_4, vert0_4) + t_4;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_4, vert1_4) + t_4;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_4, vert2_4) + t_4;
    float2  _S370 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S371 = vert0_c_4.z;
    float2  _S372 = make_float2 (_S371);
    float2  _S373 = _S370 / make_float2 (_S371);
    float2  _S374 = make_float2 (_S371 * _S371);
    float u_15 = _S373.x;
    float v_15 = _S373.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S375 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
    float _S376 = dist_coeffs_4[int(1)] + r2_15 * _S375;
    float _S377 = dist_coeffs_4[int(0)] + r2_15 * _S376;
    float radial_6 = 1.0f + r2_15 * _S377;
    float _S378 = 2.0f * dist_coeffs_4[int(4)];
    float _S379 = _S378 * u_15;
    float _S380 = 2.0f * u_15;
    float _S381 = 2.0f * dist_coeffs_4[int(5)];
    float _S382 = _S381 * u_15;
    float _S383 = 2.0f * v_15;
    float2  _S384 = _S373 * make_float2 (radial_6) + make_float2 (_S379 * v_15 + dist_coeffs_4[int(5)] * (r2_15 + _S380 * u_15) + dist_coeffs_4[int(6)] * r2_15, _S382 * v_15 + dist_coeffs_4[int(4)] * (r2_15 + _S383 * v_15) + dist_coeffs_4[int(7)] * r2_15);
    float2  _S385 = _S384 + make_float2 (dist_coeffs_4[int(8)] * _S384.x + dist_coeffs_4[int(9)] * _S384.y, 0.0f);
    float _S386 = fx_4 * _S385.x + cx_4;
    float _S387 = fy_4 * _S385.y + cy_4;
    float2  uv0_6 = make_float2 (_S386, _S387);
    float2  _S388 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S389 = vert1_c_4.z;
    float2  _S390 = make_float2 (_S389);
    float2  _S391 = _S388 / make_float2 (_S389);
    float2  _S392 = make_float2 (_S389 * _S389);
    float u_16 = _S391.x;
    float v_16 = _S391.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S393 = dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)];
    float _S394 = dist_coeffs_4[int(1)] + r2_16 * _S393;
    float _S395 = dist_coeffs_4[int(0)] + r2_16 * _S394;
    float radial_7 = 1.0f + r2_16 * _S395;
    float _S396 = _S378 * u_16;
    float _S397 = 2.0f * u_16;
    float _S398 = _S381 * u_16;
    float _S399 = 2.0f * v_16;
    float2  _S400 = _S391 * make_float2 (radial_7) + make_float2 (_S396 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + _S397 * u_16) + dist_coeffs_4[int(6)] * r2_16, _S398 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + _S399 * v_16) + dist_coeffs_4[int(7)] * r2_16);
    float2  _S401 = _S400 + make_float2 (dist_coeffs_4[int(8)] * _S400.x + dist_coeffs_4[int(9)] * _S400.y, 0.0f);
    float _S402 = fx_4 * _S401.x + cx_4;
    float _S403 = fy_4 * _S401.y + cy_4;
    float2  uv1_6 = make_float2 (_S402, _S403);
    float2  _S404 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S405 = vert2_c_4.z;
    float2  _S406 = make_float2 (_S405);
    float2  _S407 = _S404 / make_float2 (_S405);
    float2  _S408 = make_float2 (_S405 * _S405);
    float u_17 = _S407.x;
    float v_17 = _S407.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S409 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
    float _S410 = dist_coeffs_4[int(1)] + r2_17 * _S409;
    float _S411 = dist_coeffs_4[int(0)] + r2_17 * _S410;
    float radial_8 = 1.0f + r2_17 * _S411;
    float _S412 = _S378 * u_17;
    float _S413 = 2.0f * u_17;
    float _S414 = _S381 * u_17;
    float _S415 = 2.0f * v_17;
    float2  _S416 = _S407 * make_float2 (radial_8) + make_float2 (_S412 * v_17 + dist_coeffs_4[int(5)] * (r2_17 + _S413 * u_17) + dist_coeffs_4[int(6)] * r2_17, _S414 * v_17 + dist_coeffs_4[int(4)] * (r2_17 + _S415 * v_17) + dist_coeffs_4[int(7)] * r2_17);
    float2  _S417 = _S416 + make_float2 (dist_coeffs_4[int(8)] * _S416.x + dist_coeffs_4[int(9)] * _S416.y, 0.0f);
    float _S418 = fx_4 * _S417.x + cx_4;
    float _S419 = fy_4 * _S417.y + cy_4;
    float2  uv2_6 = make_float2 (_S418, _S419);
    float2  e0_4 = uv1_6 - uv0_6;
    float2  e1_4 = uv2_6 - uv1_6;
    float2  e2_0 = uv0_6 - uv2_6;
    float _S420 = e0_4.x;
    float _S421 = e1_4.y;
    float _S422 = e0_4.y;
    float _S423 = e1_4.x;
    float _S424 = _S420 * _S421 - _S422 * _S423;
    float _S425 = 1.0f - hardness_4.y;
    float _S426 = -1.0f / _S425;
    float _S427 = _S425 * _S425;
    float _S428 = (F32_max((_S386), (_S402)));
    float _S429 = (F32_min((_S386), (_S402)));
    float _S430 = (F32_max((_S387), (_S403)));
    float _S431 = (F32_min((_S387), (_S403)));
    Matrix<float, 3, 3>  _S432 = transpose_1(R_4);
    float3  _S433 = mean_4 - - s_primal_ctx_mul_0(_S432, t_4);
    float _S434 = _S433.x;
    float _S435 = _S433.y;
    float _S436 = _S433.z;
    float _S437 = _S434 * _S434 + _S435 * _S435 + _S436 * _S436;
    float _S438 = s_primal_ctx_sqrt_0(_S437);
    float x_21 = _S434 / _S438;
    float3  _S439 = make_float3 (x_21);
    float _S440 = _S438 * _S438;
    float y_10 = _S435 / _S438;
    float z_6 = _S436 / _S438;
    float3  _S441 = make_float3 (z_6);
    float _S442 = - y_10;
    float3  _S443 = make_float3 (_S442);
    float z2_9 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_4 = x_21 * x_21 - y_10 * y_10;
    float _S444 = 2.0f * x_21;
    float fS1_4 = _S444 * y_10;
    float pSH6_0 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float3  _S445 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_21;
    float3  _S446 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_10;
    float3  _S447 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S448 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S449 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_6;
    float _S450 = 1.86588168144226074f * z2_9 - 1.11952900886535645f;
    float pSH12_0 = z_6 * _S450;
    float3  _S451 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_21;
    float3  _S452 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_10;
    float3  _S453 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S454 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S455 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_21 * fC1_4 - y_10 * fS1_4);
    float3  _S456 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_21 * fS1_4 + y_10 * fC1_4);
    float3  _S457 = make_float3 (pSH9_0);
    float3  color_4 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S442) * sh_coeffs_4[int(1)] + make_float3 (z_6) * sh_coeffs_4[int(2)] - make_float3 (x_21) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]);
    float3  _S458 = color_4 + ch_coeffs_4[int(0)] + make_float3 (0.5f);
    float3  _S459 = make_float3 (0.0f);
    float3  _S460 = color_4 - ch_coeffs_4[int(0)] * make_float3 (0.5f);
    float _S461 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S462 = make_float3 (_S461);
    float3  _S463 = ch_coeffs_4[int(1)] * make_float3 (_S461);
    float3  _S464 = _S460 + _S463 + make_float3 (0.5f);
    float3  _S465 = _S460 - _S463 + make_float3 (0.5f);
    float3  _S466 = vert1_c_4 - vert0_c_4;
    float3  _S467 = vert2_c_4 - vert0_c_4;
    float3  _S468 = s_primal_ctx_cross_0(_S466, _S467);
    float3  _S469 = normalize_1(_S468);
    float3  _S470 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S469, mean_c_4)))))) * v_normal_0;
    float3  _S471 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S472;
    (&_S472)->primal_0 = _S469;
    (&_S472)->differential_0 = _S471;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S473;
    (&_S473)->primal_0 = mean_c_4;
    (&_S473)->differential_0 = _S471;
    s_bwd_prop_dot_0(&_S472, &_S473, 0.0f);
    float3  _S474 = _S470 + _S472.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S475;
    (&_S475)->primal_0 = _S468;
    (&_S475)->differential_0 = _S471;
    s_bwd_normalize_impl_0(&_S475, _S474);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S476;
    (&_S476)->primal_0 = _S466;
    (&_S476)->differential_0 = _S471;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S477;
    (&_S477)->primal_0 = _S467;
    (&_S477)->differential_0 = _S471;
    s_bwd_prop_cross_0(&_S476, &_S477, _S475.differential_0);
    float3  _S478 = - _S477.differential_0;
    float3  _S479 = - _S476.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S480;
    (&_S480)->primal_0 = _S465;
    (&_S480)->differential_0 = _S471;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481;
    (&_S481)->primal_0 = _S459;
    (&_S481)->differential_0 = _S471;
    s_bwd_prop_max_0(&_S480, &_S481, v_rgbs_0[int(2)]);
    float3  _S482 = - _S480.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S483;
    (&_S483)->primal_0 = _S464;
    (&_S483)->differential_0 = _S471;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S459;
    (&_S484)->differential_0 = _S471;
    s_bwd_prop_max_0(&_S483, &_S484, v_rgbs_0[int(1)]);
    float3  _S485 = _S462 * (_S482 + _S483.differential_0);
    float3  _S486 = _S480.differential_0 + _S483.differential_0;
    float3  _S487 = make_float3 (0.5f) * - _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S488;
    (&_S488)->primal_0 = _S458;
    (&_S488)->differential_0 = _S471;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S489;
    (&_S489)->primal_0 = _S459;
    (&_S489)->differential_0 = _S471;
    s_bwd_prop_max_0(&_S488, &_S489, v_rgbs_0[int(0)]);
    float3  _S490 = _S487 + _S488.differential_0;
    float3  _S491 = _S486 + _S488.differential_0;
    float3  _S492 = _S456 * _S491;
    float3  _S493 = sh_coeffs_4[int(15)] * _S491;
    float3  _S494 = _S454 * _S491;
    float3  _S495 = sh_coeffs_4[int(14)] * _S491;
    float3  _S496 = _S452 * _S491;
    float3  _S497 = sh_coeffs_4[int(13)] * _S491;
    float3  _S498 = _S451 * _S491;
    float3  _S499 = sh_coeffs_4[int(12)] * _S491;
    float3  _S500 = _S453 * _S491;
    float3  _S501 = sh_coeffs_4[int(11)] * _S491;
    float3  _S502 = _S455 * _S491;
    float3  _S503 = sh_coeffs_4[int(10)] * _S491;
    float3  _S504 = _S457 * _S491;
    float3  _S505 = sh_coeffs_4[int(9)] * _S491;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S505.x + _S505.y + _S505.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S493.x + _S493.y + _S493.z);
    float _S506 = _S503.x + _S503.y + _S503.z;
    float _S507 = _S495.x + _S495.y + _S495.z;
    float _S508 = _S501.x + _S501.y + _S501.z;
    float _S509 = _S497.x + _S497.y + _S497.z;
    float _S510 = _S499.x + _S499.y + _S499.z;
    float _S511 = - s_diff_fC2_T_0;
    float3  _S512 = _S448 * _S491;
    float3  _S513 = sh_coeffs_4[int(8)] * _S491;
    float3  _S514 = _S446 * _S491;
    float3  _S515 = sh_coeffs_4[int(7)] * _S491;
    float3  _S516 = _S445 * _S491;
    float3  _S517 = sh_coeffs_4[int(6)] * _S491;
    float3  _S518 = _S447 * _S491;
    float3  _S519 = sh_coeffs_4[int(5)] * _S491;
    float3  _S520 = _S449 * _S491;
    float3  _S521 = sh_coeffs_4[int(4)] * _S491;
    float _S522 = _S519.x + _S519.y + _S519.z;
    float _S523 = _S515.x + _S515.y + _S515.z;
    float _S524 = fTmp1B_4 * _S506 + x_21 * s_diff_fS2_T_0 + y_10 * _S511 + 0.54627424478530884f * (_S521.x + _S521.y + _S521.z);
    float _S525 = fTmp1B_4 * _S507 + y_10 * s_diff_fS2_T_0 + x_21 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S513.x + _S513.y + _S513.z);
    float _S526 = y_10 * - _S525;
    float _S527 = x_21 * _S525;
    float _S528 = z_6 * (1.86588168144226074f * (z_6 * _S510) + -2.28522896766662598f * (y_10 * _S508 + x_21 * _S509) + 0.94617468118667603f * (_S517.x + _S517.y + _S517.z));
    float3  _S529 = make_float3 (0.48860251903533936f) * _S491;
    float3  _S530 = - _S529;
    float3  _S531 = _S439 * _S530;
    float3  _S532 = sh_coeffs_4[int(3)] * _S530;
    float3  _S533 = _S441 * _S529;
    float3  _S534 = sh_coeffs_4[int(2)] * _S529;
    float3  _S535 = _S443 * _S529;
    float3  _S536 = sh_coeffs_4[int(1)] * _S529;
    float _S537 = (_S450 * _S510 + 1.44530570507049561f * (fS1_4 * _S506 + fC1_4 * _S507) + -1.09254848957061768f * (y_10 * _S522 + x_21 * _S523) + _S528 + _S528 + _S534.x + _S534.y + _S534.z) / _S440;
    float _S538 = _S438 * _S537;
    float _S539 = (fTmp0C_4 * _S508 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S511 + fTmp0B_4 * _S522 + _S444 * _S524 + _S526 + _S526 + - (_S536.x + _S536.y + _S536.z)) / _S440;
    float _S540 = _S438 * _S539;
    float _S541 = (fTmp0C_4 * _S509 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S523 + 2.0f * (y_10 * _S524) + _S527 + _S527 + _S532.x + _S532.y + _S532.z) / _S440;
    float _S542 = _S438 * _S541;
    float _S543 = _S436 * - _S537 + _S435 * - _S539 + _S434 * - _S541;
    DiffPair_float_0 _S544;
    (&_S544)->primal_0 = _S437;
    (&_S544)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S544, _S543);
    float _S545 = _S436 * _S544.differential_0;
    float _S546 = _S435 * _S544.differential_0;
    float _S547 = _S434 * _S544.differential_0;
    float3  _S548 = make_float3 (0.282094806432724f) * _S491;
    float3  _S549 = make_float3 (_S542 + _S547 + _S547, _S540 + _S546 + _S546, _S538 + _S545 + _S545);
    float3  _S550 = - - _S549;
    Matrix<float, 3, 3>  _S551 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S552;
    (&_S552)->primal_0 = _S432;
    (&_S552)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S553;
    (&_S553)->primal_0 = t_4;
    (&_S553)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S552, &_S553, _S550);
    Matrix<float, 3, 3>  _S554 = transpose_1(_S552.differential_0);
    float3  _S555 = make_float3 (0.3333333432674408f) * make_float3 (0.0f, 0.0f, v_depth_0);
    DiffPair_float_0 _S556;
    (&_S556)->primal_0 = _S431;
    (&_S556)->differential_0 = 0.0f;
    DiffPair_float_0 _S557;
    (&_S557)->primal_0 = _S419;
    (&_S557)->differential_0 = 0.0f;
    _d_min_0(&_S556, &_S557, 0.0f);
    DiffPair_float_0 _S558;
    (&_S558)->primal_0 = _S387;
    (&_S558)->differential_0 = 0.0f;
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = _S403;
    (&_S559)->differential_0 = 0.0f;
    _d_min_0(&_S558, &_S559, _S556.differential_0);
    DiffPair_float_0 _S560;
    (&_S560)->primal_0 = _S430;
    (&_S560)->differential_0 = 0.0f;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S419;
    (&_S561)->differential_0 = 0.0f;
    _d_max_0(&_S560, &_S561, 0.0f);
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S387;
    (&_S562)->differential_0 = 0.0f;
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = _S403;
    (&_S563)->differential_0 = 0.0f;
    _d_max_0(&_S562, &_S563, _S560.differential_0);
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = _S429;
    (&_S564)->differential_0 = 0.0f;
    DiffPair_float_0 _S565;
    (&_S565)->primal_0 = _S418;
    (&_S565)->differential_0 = 0.0f;
    _d_min_0(&_S564, &_S565, 0.0f);
    DiffPair_float_0 _S566;
    (&_S566)->primal_0 = _S386;
    (&_S566)->differential_0 = 0.0f;
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S402;
    (&_S567)->differential_0 = 0.0f;
    _d_min_0(&_S566, &_S567, _S564.differential_0);
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S428;
    (&_S568)->differential_0 = 0.0f;
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S418;
    (&_S569)->differential_0 = 0.0f;
    _d_max_0(&_S568, &_S569, 0.0f);
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S386;
    (&_S570)->differential_0 = 0.0f;
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S402;
    (&_S571)->differential_0 = 0.0f;
    _d_max_0(&_S570, &_S571, _S568.differential_0);
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S426;
    (&_S572)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S572, -0.0f);
    float _S573 = - (-1.0f * - (_S572.differential_0 / _S427));
    float2  _S574 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S575;
    (&_S575)->primal_0 = e2_0;
    (&_S575)->differential_0 = _S574;
    s_bwd_length_impl_1(&_S575, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S576;
    (&_S576)->primal_0 = e1_4;
    (&_S576)->differential_0 = _S574;
    s_bwd_length_impl_1(&_S576, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S577;
    (&_S577)->primal_0 = e0_4;
    (&_S577)->differential_0 = _S574;
    s_bwd_length_impl_1(&_S577, 0.0f);
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = _S424;
    (&_S578)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S578, -0.0f);
    float _S579 = - _S578.differential_0;
    float2  _S580 = _S576.differential_0 + make_float2 (_S422 * _S579, _S420 * _S578.differential_0);
    float2  _S581 = _S577.differential_0 + make_float2 (_S421 * _S578.differential_0, _S423 * _S579);
    float2  _S582 = - _S575.differential_0 + _S580;
    float _S583 = fx_4 * (_S565.differential_0 + _S569.differential_0 + _S582.x);
    float2  _S584 = make_float2 (_S583, fy_4 * (_S557.differential_0 + _S561.differential_0 + _S582.y)) + make_float2 (dist_coeffs_4[int(8)] * _S583, dist_coeffs_4[int(9)] * _S583);
    float2  _S585 = _S407 * _S584;
    float _S586 = dist_coeffs_4[int(4)] * _S584.y;
    float _S587 = dist_coeffs_4[int(5)] * _S584.x;
    float _S588 = _S585.x + _S585.y;
    float _S589 = r2_17 * _S588;
    float _S590 = r2_17 * _S589;
    float _S591 = dist_coeffs_4[int(7)] * _S584.y + _S586 + dist_coeffs_4[int(6)] * _S584.x + _S587 + _S411 * _S588 + _S410 * _S589 + _S409 * _S590 + dist_coeffs_4[int(3)] * (r2_17 * _S590);
    float _S592 = v_17 * _S591;
    float _S593 = u_17 * _S591;
    float2  _S594 = (make_float2 (radial_8) * _S584 + make_float2 (_S381 * (v_17 * _S584.y) + _S413 * _S587 + 2.0f * (u_17 * _S587) + _S378 * (v_17 * _S584.x) + _S593 + _S593, _S415 * _S586 + 2.0f * (v_17 * _S586) + _S414 * _S584.y + _S412 * _S584.x + _S592 + _S592)) / _S408;
    float2  _S595 = _S404 * - _S594;
    float2  _S596 = _S406 * _S594;
    float2  _S597 = - _S580 + _S581;
    float _S598 = fx_4 * (_S567.differential_0 + _S571.differential_0 + _S597.x);
    float2  _S599 = make_float2 (_S598, fy_4 * (_S559.differential_0 + _S563.differential_0 + _S597.y)) + make_float2 (dist_coeffs_4[int(8)] * _S598, dist_coeffs_4[int(9)] * _S598);
    float2  _S600 = _S391 * _S599;
    float _S601 = dist_coeffs_4[int(4)] * _S599.y;
    float _S602 = dist_coeffs_4[int(5)] * _S599.x;
    float _S603 = _S600.x + _S600.y;
    float _S604 = r2_16 * _S603;
    float _S605 = r2_16 * _S604;
    float _S606 = dist_coeffs_4[int(7)] * _S599.y + _S601 + dist_coeffs_4[int(6)] * _S599.x + _S602 + _S395 * _S603 + _S394 * _S604 + _S393 * _S605 + dist_coeffs_4[int(3)] * (r2_16 * _S605);
    float _S607 = v_16 * _S606;
    float _S608 = u_16 * _S606;
    float2  _S609 = (make_float2 (radial_7) * _S599 + make_float2 (_S381 * (v_16 * _S599.y) + _S397 * _S602 + 2.0f * (u_16 * _S602) + _S378 * (v_16 * _S599.x) + _S608 + _S608, _S399 * _S601 + 2.0f * (v_16 * _S601) + _S398 * _S599.y + _S396 * _S599.x + _S607 + _S607)) / _S392;
    float2  _S610 = _S388 * - _S609;
    float2  _S611 = _S390 * _S609;
    float _S612 = _S610.x + _S610.y;
    float2  _S613 = _S575.differential_0 + - _S581;
    float _S614 = fx_4 * (_S566.differential_0 + _S570.differential_0 + _S613.x);
    float2  _S615 = make_float2 (_S614, fy_4 * (_S558.differential_0 + _S562.differential_0 + _S613.y)) + make_float2 (dist_coeffs_4[int(8)] * _S614, dist_coeffs_4[int(9)] * _S614);
    float2  _S616 = _S373 * _S615;
    float _S617 = dist_coeffs_4[int(4)] * _S615.y;
    float _S618 = dist_coeffs_4[int(5)] * _S615.x;
    float _S619 = _S616.x + _S616.y;
    float _S620 = r2_15 * _S619;
    float _S621 = r2_15 * _S620;
    float _S622 = dist_coeffs_4[int(7)] * _S615.y + _S617 + dist_coeffs_4[int(6)] * _S615.x + _S618 + _S377 * _S619 + _S376 * _S620 + _S375 * _S621 + dist_coeffs_4[int(3)] * (r2_15 * _S621);
    float _S623 = v_15 * _S622;
    float _S624 = u_15 * _S622;
    float2  _S625 = (make_float2 (radial_6) * _S615 + make_float2 (_S381 * (v_15 * _S615.y) + _S380 * _S618 + 2.0f * (u_15 * _S618) + _S378 * (v_15 * _S615.x) + _S624 + _S624, _S383 * _S617 + 2.0f * (v_15 * _S617) + _S382 * _S615.y + _S379 * _S615.x + _S623 + _S623)) / _S374;
    float2  _S626 = _S370 * - _S625;
    float2  _S627 = _S372 * _S625;
    float _S628 = _S626.x + _S626.y;
    float3  _S629 = _S477.differential_0 + _S555 + make_float3 (_S596.x, _S596.y, _S595.x + _S595.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S630;
    (&_S630)->primal_0 = R_4;
    (&_S630)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S631;
    (&_S631)->primal_0 = vert2_4;
    (&_S631)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S630, &_S631, _S629);
    float3  _S632 = _S476.differential_0 + _S555 + make_float3 (_S611.x, _S611.y, _S612);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S633;
    (&_S633)->primal_0 = R_4;
    (&_S633)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S634;
    (&_S634)->primal_0 = vert1_4;
    (&_S634)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S633, &_S634, _S632);
    float3  _S635 = _S478 + _S479 + _S555 + make_float3 (_S627.x, _S627.y, _S628);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S636;
    (&_S636)->primal_0 = R_4;
    (&_S636)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S637;
    (&_S637)->primal_0 = vert0_4;
    (&_S637)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S636, &_S637, _S635);
    float3  _S638 = v_verts_0[int(2)] + _S631.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S639;
    (&_S639)->primal_0 = _S364;
    (&_S639)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S640;
    (&_S640)->primal_0 = _S369;
    (&_S640)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S639, &_S640, _S638);
    float _S641 = - _S640.differential_0.y;
    float _S642 = _S368 * _S640.differential_0.x;
    float _S643 = - (_S359 * _S640.differential_0.x);
    float3  _S644 = v_verts_0[int(1)] + _S634.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S645;
    (&_S645)->primal_0 = _S364;
    (&_S645)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S646;
    (&_S646)->primal_0 = _S367;
    (&_S646)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S645, &_S646, _S644);
    float _S647 = _S359 * _S646.differential_0.x;
    float _S648 = _S366 * _S646.differential_0.x;
    float3  _S649 = v_verts_0[int(0)] + _S637.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S650;
    (&_S650)->primal_0 = _S364;
    (&_S650)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S651;
    (&_S651)->primal_0 = _S365;
    (&_S651)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S650, &_S651, _S649);
    Matrix<float, 3, 3>  _S652 = transpose_1(_S639.differential_0 + _S645.differential_0 + _S650.differential_0);
    float _S653 = 2.0f * - _S652.rows[int(2)].z;
    float _S654 = 2.0f * _S652.rows[int(2)].y;
    float _S655 = 2.0f * _S652.rows[int(2)].x;
    float _S656 = 2.0f * _S652.rows[int(1)].z;
    float _S657 = 2.0f * - _S652.rows[int(1)].y;
    float _S658 = 2.0f * _S652.rows[int(1)].x;
    float _S659 = 2.0f * _S652.rows[int(0)].z;
    float _S660 = 2.0f * _S652.rows[int(0)].y;
    float _S661 = 2.0f * - _S652.rows[int(0)].x;
    float _S662 = - _S658 + _S660;
    float _S663 = _S655 + - _S659;
    float _S664 = - _S654 + _S656;
    float _S665 = _S654 + _S656;
    float _S666 = _S655 + _S659;
    float _S667 = _S658 + _S660;
    float _S668 = _S362.w * (_S657 + _S661);
    float _S669 = _S362.z * (_S653 + _S661);
    float _S670 = _S362.y * (_S653 + _S657);
    float _S671 = _S362.x * _S662 + _S362.z * _S665 + _S362.y * _S666 + _S668 + _S668;
    float _S672 = _S362.x * _S663 + _S362.w * _S665 + _S362.y * _S667 + _S669 + _S669;
    float _S673 = _S362.x * _S664 + _S362.w * _S666 + _S362.z * _S667 + _S670 + _S670;
    float _S674 = _S362.w * _S662 + _S362.z * _S663 + _S362.y * _S664;
    float4  _S675 = make_float4 (0.0f);
    float4  _S676 = _S675;
    *&((&_S676)->w) = _S671;
    *&((&_S676)->z) = _S672;
    *&((&_S676)->y) = _S673;
    *&((&_S676)->x) = _S674;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S677;
    (&_S677)->primal_0 = quat_4;
    (&_S677)->differential_0 = _S675;
    s_bwd_normalize_impl_1(&_S677, _S676);
    float _S678 = _S643 + _S647;
    float _S679 = 0.5f * - _S678;
    float _S680 = _S641 + _S646.differential_0.y;
    DiffPair_float_0 _S681;
    (&_S681)->primal_0 = _S360;
    (&_S681)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S681, _S680);
    float _S682 = _S679 + _S681.differential_0;
    float _S683 = _S642 + _S648 + _S651.differential_0.x;
    DiffPair_float_0 _S684;
    (&_S684)->primal_0 = _S358;
    (&_S684)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S684, _S683);
    float _S685 = _S679 + _S684.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S686;
    (&_S686)->primal_0 = R_4;
    (&_S686)->differential_0 = _S551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S687;
    (&_S687)->primal_0 = mean_4;
    (&_S687)->differential_0 = _S471;
    s_bwd_prop_mul_0(&_S686, &_S687, _S473.differential_0);
    float3  _S688 = _S553.differential_0 + _S629 + _S632 + _S635 + _S473.differential_0;
    Matrix<float, 3, 3>  _S689 = _S554 + _S630.differential_0 + _S633.differential_0 + _S636.differential_0 + _S686.differential_0;
    FixedArray<float3 , 2>  _S690;
    _S690[int(0)] = _S471;
    _S690[int(1)] = _S471;
    _S690[int(1)] = _S485;
    _S690[int(0)] = _S490;
    FixedArray<float3 , 16>  _S691;
    _S691[int(0)] = _S471;
    _S691[int(1)] = _S471;
    _S691[int(2)] = _S471;
    _S691[int(3)] = _S471;
    _S691[int(4)] = _S471;
    _S691[int(5)] = _S471;
    _S691[int(6)] = _S471;
    _S691[int(7)] = _S471;
    _S691[int(8)] = _S471;
    _S691[int(9)] = _S471;
    _S691[int(10)] = _S471;
    _S691[int(11)] = _S471;
    _S691[int(12)] = _S471;
    _S691[int(13)] = _S471;
    _S691[int(14)] = _S471;
    _S691[int(15)] = _S471;
    _S691[int(15)] = _S492;
    _S691[int(14)] = _S494;
    _S691[int(13)] = _S496;
    _S691[int(12)] = _S498;
    _S691[int(11)] = _S500;
    _S691[int(10)] = _S502;
    _S691[int(9)] = _S504;
    _S691[int(8)] = _S512;
    _S691[int(7)] = _S514;
    _S691[int(6)] = _S516;
    _S691[int(5)] = _S518;
    _S691[int(4)] = _S520;
    _S691[int(3)] = _S531;
    _S691[int(2)] = _S533;
    _S691[int(1)] = _S535;
    _S691[int(0)] = _S548;
    float2  _S692 = make_float2 (0.0f, _S573);
    float3  _S693 = make_float3 (_S685, _S682, _S678);
    *v_mean_0 = _S549 + _S638 + _S644 + _S649 + _S687.differential_0;
    *v_quat_0 = _S677.differential_0;
    *v_scale_0 = _S693;
    *v_hardness_0 = _S692;
    *v_sh_coeffs_0 = _S691;
    *v_ch_coeffs_0 = _S690;
    *v_R_0 = _S689;
    *v_t_0 = _S688;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S694, float _S695)
{
    return (F32_atan2((_S694), (_S695)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S696, DiffPair_float_0 * _S697, float _S698)
{
    _d_atan2_0(_S696, _S697, _S698);
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_5, float4  quat_5, float3  scale_5, float2  hardness_5, FixedArray<float3 , 16>  sh_coeffs_5, FixedArray<float3 , 2>  ch_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float v_depth_1, FixedArray<float3 , 3>  v_verts_1, FixedArray<float3 , 3>  v_rgbs_1, float3  v_normal_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, mean_5) + t_5;
    float _S699 = scale_5.x;
    float _S700 = s_primal_ctx_exp_0(_S699);
    float _S701 = scale_5.y;
    float _S702 = s_primal_ctx_exp_0(_S701);
    float sz_5 = scale_5.z - 0.5f * (_S699 + _S701);
    float4  _S703 = normalize_0(quat_5);
    float _S704 = _S703.y;
    float x2_5 = _S704 * _S704;
    float y2_5 = _S703.z * _S703.z;
    float z2_10 = _S703.w * _S703.w;
    float xy_5 = _S703.y * _S703.z;
    float xz_5 = _S703.y * _S703.w;
    float yz_5 = _S703.z * _S703.w;
    float wx_5 = _S703.x * _S703.y;
    float wy_5 = _S703.x * _S703.z;
    float wz_5 = _S703.x * _S703.w;
    Matrix<float, 3, 3>  _S705 = transpose_1(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_10), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_10), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  _S706 = make_float3 (_S700, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S705, _S706) + mean_5;
    float _S707 = -0.5f + sz_5;
    float3  _S708 = make_float3 (_S700 * _S707, _S702, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S705, _S708) + mean_5;
    float _S709 = -0.5f - sz_5;
    float3  _S710 = make_float3 (_S700 * _S709, - _S702, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S705, _S710) + mean_5;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_5, vert0_5) + t_5;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_5, vert1_5) + t_5;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_5, vert2_5) + t_5;
    float2  _S711 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S712 = length_0(_S711);
    float _S713 = vert0_c_5.z;
    float _S714 = s_primal_ctx_atan2_0(_S712, _S713);
    bool _S715 = _S714 < 0.00100000004749745f;
    float k_2;
    float _S716;
    float _S717;
    float _S718;
    if(_S715)
    {
        float _S719 = 1.0f - _S714 * _S714 / 3.0f;
        float _S720 = _S713 * _S713;
        k_2 = _S719 / _S713;
        _S716 = _S720;
        _S717 = _S719;
        _S718 = 0.0f;
    }
    else
    {
        float _S721 = _S712 * _S712;
        k_2 = _S714 / _S712;
        _S716 = 0.0f;
        _S717 = 0.0f;
        _S718 = _S721;
    }
    float2  _S722 = make_float2 (k_2);
    float2  _S723 = _S711 * make_float2 (k_2);
    float u_18 = _S723.x;
    float v_18 = _S723.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S724 = dist_coeffs_5[int(2)] + r2_18 * dist_coeffs_5[int(3)];
    float _S725 = dist_coeffs_5[int(1)] + r2_18 * _S724;
    float _S726 = dist_coeffs_5[int(0)] + r2_18 * _S725;
    float radial_9 = 1.0f + r2_18 * _S726;
    float _S727 = 2.0f * dist_coeffs_5[int(4)];
    float _S728 = _S727 * u_18;
    float _S729 = 2.0f * u_18;
    float _S730 = 2.0f * dist_coeffs_5[int(5)];
    float _S731 = _S730 * u_18;
    float _S732 = 2.0f * v_18;
    float2  _S733 = _S723 * make_float2 (radial_9) + make_float2 (_S728 * v_18 + dist_coeffs_5[int(5)] * (r2_18 + _S729 * u_18) + dist_coeffs_5[int(6)] * r2_18, _S731 * v_18 + dist_coeffs_5[int(4)] * (r2_18 + _S732 * v_18) + dist_coeffs_5[int(7)] * r2_18);
    float2  _S734 = _S733 + make_float2 (dist_coeffs_5[int(8)] * _S733.x + dist_coeffs_5[int(9)] * _S733.y, 0.0f);
    float _S735 = fx_5 * _S734.x + cx_5;
    float _S736 = fy_5 * _S734.y + cy_5;
    float2  uv0_7 = make_float2 (_S735, _S736);
    float2  _S737 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S738 = length_0(_S737);
    float _S739 = vert1_c_5.z;
    float _S740 = s_primal_ctx_atan2_0(_S738, _S739);
    bool _S741 = _S740 < 0.00100000004749745f;
    float _S742;
    float _S743;
    float _S744;
    if(_S741)
    {
        float _S745 = 1.0f - _S740 * _S740 / 3.0f;
        float _S746 = _S739 * _S739;
        k_2 = _S745 / _S739;
        _S742 = _S746;
        _S743 = _S745;
        _S744 = 0.0f;
    }
    else
    {
        float _S747 = _S738 * _S738;
        k_2 = _S740 / _S738;
        _S742 = 0.0f;
        _S743 = 0.0f;
        _S744 = _S747;
    }
    float2  _S748 = make_float2 (k_2);
    float2  _S749 = _S737 * make_float2 (k_2);
    float u_19 = _S749.x;
    float v_19 = _S749.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S750 = dist_coeffs_5[int(2)] + r2_19 * dist_coeffs_5[int(3)];
    float _S751 = dist_coeffs_5[int(1)] + r2_19 * _S750;
    float _S752 = dist_coeffs_5[int(0)] + r2_19 * _S751;
    float radial_10 = 1.0f + r2_19 * _S752;
    float _S753 = _S727 * u_19;
    float _S754 = 2.0f * u_19;
    float _S755 = _S730 * u_19;
    float _S756 = 2.0f * v_19;
    float2  _S757 = _S749 * make_float2 (radial_10) + make_float2 (_S753 * v_19 + dist_coeffs_5[int(5)] * (r2_19 + _S754 * u_19) + dist_coeffs_5[int(6)] * r2_19, _S755 * v_19 + dist_coeffs_5[int(4)] * (r2_19 + _S756 * v_19) + dist_coeffs_5[int(7)] * r2_19);
    float2  _S758 = _S757 + make_float2 (dist_coeffs_5[int(8)] * _S757.x + dist_coeffs_5[int(9)] * _S757.y, 0.0f);
    float _S759 = fx_5 * _S758.x + cx_5;
    float _S760 = fy_5 * _S758.y + cy_5;
    float2  uv1_7 = make_float2 (_S759, _S760);
    float2  _S761 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S762 = length_0(_S761);
    float _S763 = vert2_c_5.z;
    float _S764 = s_primal_ctx_atan2_0(_S762, _S763);
    bool _S765 = _S764 < 0.00100000004749745f;
    float _S766;
    float _S767;
    float _S768;
    if(_S765)
    {
        float _S769 = 1.0f - _S764 * _S764 / 3.0f;
        float _S770 = _S763 * _S763;
        k_2 = _S769 / _S763;
        _S766 = _S770;
        _S767 = _S769;
        _S768 = 0.0f;
    }
    else
    {
        float _S771 = _S762 * _S762;
        k_2 = _S764 / _S762;
        _S766 = 0.0f;
        _S767 = 0.0f;
        _S768 = _S771;
    }
    float2  _S772 = make_float2 (k_2);
    float2  _S773 = _S761 * make_float2 (k_2);
    float u_20 = _S773.x;
    float v_20 = _S773.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S774 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
    float _S775 = dist_coeffs_5[int(1)] + r2_20 * _S774;
    float _S776 = dist_coeffs_5[int(0)] + r2_20 * _S775;
    float radial_11 = 1.0f + r2_20 * _S776;
    float _S777 = _S727 * u_20;
    float _S778 = 2.0f * u_20;
    float _S779 = _S730 * u_20;
    float _S780 = 2.0f * v_20;
    float2  _S781 = _S773 * make_float2 (radial_11) + make_float2 (_S777 * v_20 + dist_coeffs_5[int(5)] * (r2_20 + _S778 * u_20) + dist_coeffs_5[int(6)] * r2_20, _S779 * v_20 + dist_coeffs_5[int(4)] * (r2_20 + _S780 * v_20) + dist_coeffs_5[int(7)] * r2_20);
    float2  _S782 = _S781 + make_float2 (dist_coeffs_5[int(8)] * _S781.x + dist_coeffs_5[int(9)] * _S781.y, 0.0f);
    float _S783 = fx_5 * _S782.x + cx_5;
    float _S784 = fy_5 * _S782.y + cy_5;
    float2  uv2_7 = make_float2 (_S783, _S784);
    float2  e0_5 = uv1_7 - uv0_7;
    float2  e1_5 = uv2_7 - uv1_7;
    float2  e2_1 = uv0_7 - uv2_7;
    float _S785 = e0_5.x;
    float _S786 = e1_5.y;
    float _S787 = e0_5.y;
    float _S788 = e1_5.x;
    float _S789 = _S785 * _S786 - _S787 * _S788;
    float _S790 = 1.0f - hardness_5.y;
    float _S791 = -1.0f / _S790;
    float _S792 = _S790 * _S790;
    float _S793 = (F32_max((_S735), (_S759)));
    float _S794 = (F32_min((_S735), (_S759)));
    float _S795 = (F32_max((_S736), (_S760)));
    float _S796 = (F32_min((_S736), (_S760)));
    float3  _S797 = (vert0_c_5 + vert1_c_5 + vert2_c_5) / make_float3 (3.0f);
    float x_22 = _S797.x;
    float y_11 = _S797.y;
    float z_7 = _S797.z;
    float _S798 = z_7 * z_7;
    float _S799 = _S798 * z_7;
    float _S800 = x_22 * x_22 + y_11 * y_11;
    float _S801 = 0.001953125f * _S800;
    Matrix<float, 3, 3>  _S802 = transpose_1(R_5);
    float3  _S803 = mean_5 - - s_primal_ctx_mul_0(_S802, t_5);
    float _S804 = _S803.x;
    float _S805 = _S803.y;
    float _S806 = _S803.z;
    float _S807 = _S804 * _S804 + _S805 * _S805 + _S806 * _S806;
    float _S808 = s_primal_ctx_sqrt_0(_S807);
    float x_23 = _S804 / _S808;
    float3  _S809 = make_float3 (x_23);
    float _S810 = _S808 * _S808;
    float y_12 = _S805 / _S808;
    float z_8 = _S806 / _S808;
    float3  _S811 = make_float3 (z_8);
    float _S812 = - y_12;
    float3  _S813 = make_float3 (_S812);
    float z2_11 = z_8 * z_8;
    float fTmp0B_5 = -1.09254848957061768f * z_8;
    float fC1_5 = x_23 * x_23 - y_12 * y_12;
    float _S814 = 2.0f * x_23;
    float fS1_5 = _S814 * y_12;
    float pSH6_1 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float3  _S815 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_23;
    float3  _S816 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_12;
    float3  _S817 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S818 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S819 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_8;
    float _S820 = 1.86588168144226074f * z2_11 - 1.11952900886535645f;
    float pSH12_1 = z_8 * _S820;
    float3  _S821 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_23;
    float3  _S822 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_12;
    float3  _S823 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S824 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S825 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_23 * fC1_5 - y_12 * fS1_5);
    float3  _S826 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_23 * fS1_5 + y_12 * fC1_5);
    float3  _S827 = make_float3 (pSH9_1);
    float3  color_5 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S812) * sh_coeffs_5[int(1)] + make_float3 (z_8) * sh_coeffs_5[int(2)] - make_float3 (x_23) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]);
    float3  _S828 = color_5 + ch_coeffs_5[int(0)] + make_float3 (0.5f);
    float3  _S829 = make_float3 (0.0f);
    float3  _S830 = color_5 - ch_coeffs_5[int(0)] * make_float3 (0.5f);
    float _S831 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S832 = make_float3 (_S831);
    float3  _S833 = ch_coeffs_5[int(1)] * make_float3 (_S831);
    float3  _S834 = _S830 + _S833 + make_float3 (0.5f);
    float3  _S835 = _S830 - _S833 + make_float3 (0.5f);
    float3  _S836 = vert1_c_5 - vert0_c_5;
    float3  _S837 = vert2_c_5 - vert0_c_5;
    float3  _S838 = s_primal_ctx_cross_0(_S836, _S837);
    float3  _S839 = normalize_1(_S838);
    float3  _S840 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S839, mean_c_5)))))) * v_normal_1;
    float3  _S841 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S842;
    (&_S842)->primal_0 = _S839;
    (&_S842)->differential_0 = _S841;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S843;
    (&_S843)->primal_0 = mean_c_5;
    (&_S843)->differential_0 = _S841;
    s_bwd_prop_dot_0(&_S842, &_S843, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S844 = _S843;
    float3  _S845 = _S840 + _S842.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S846;
    (&_S846)->primal_0 = _S838;
    (&_S846)->differential_0 = _S841;
    s_bwd_normalize_impl_0(&_S846, _S845);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S847;
    (&_S847)->primal_0 = _S836;
    (&_S847)->differential_0 = _S841;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S848;
    (&_S848)->primal_0 = _S837;
    (&_S848)->differential_0 = _S841;
    s_bwd_prop_cross_0(&_S847, &_S848, _S846.differential_0);
    float3  _S849 = - _S848.differential_0;
    float3  _S850 = - _S847.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S851;
    (&_S851)->primal_0 = _S835;
    (&_S851)->differential_0 = _S841;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S852;
    (&_S852)->primal_0 = _S829;
    (&_S852)->differential_0 = _S841;
    s_bwd_prop_max_0(&_S851, &_S852, v_rgbs_1[int(2)]);
    float3  _S853 = - _S851.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S854;
    (&_S854)->primal_0 = _S834;
    (&_S854)->differential_0 = _S841;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S855;
    (&_S855)->primal_0 = _S829;
    (&_S855)->differential_0 = _S841;
    s_bwd_prop_max_0(&_S854, &_S855, v_rgbs_1[int(1)]);
    float3  _S856 = _S832 * (_S853 + _S854.differential_0);
    float3  _S857 = _S851.differential_0 + _S854.differential_0;
    float3  _S858 = make_float3 (0.5f) * - _S857;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S859;
    (&_S859)->primal_0 = _S828;
    (&_S859)->differential_0 = _S841;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S860;
    (&_S860)->primal_0 = _S829;
    (&_S860)->differential_0 = _S841;
    s_bwd_prop_max_0(&_S859, &_S860, v_rgbs_1[int(0)]);
    float3  _S861 = _S858 + _S859.differential_0;
    float3  _S862 = _S857 + _S859.differential_0;
    float3  _S863 = _S826 * _S862;
    float3  _S864 = sh_coeffs_5[int(15)] * _S862;
    float3  _S865 = _S824 * _S862;
    float3  _S866 = sh_coeffs_5[int(14)] * _S862;
    float3  _S867 = _S822 * _S862;
    float3  _S868 = sh_coeffs_5[int(13)] * _S862;
    float3  _S869 = _S821 * _S862;
    float3  _S870 = sh_coeffs_5[int(12)] * _S862;
    float3  _S871 = _S823 * _S862;
    float3  _S872 = sh_coeffs_5[int(11)] * _S862;
    float3  _S873 = _S825 * _S862;
    float3  _S874 = sh_coeffs_5[int(10)] * _S862;
    float3  _S875 = _S827 * _S862;
    float3  _S876 = sh_coeffs_5[int(9)] * _S862;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S876.x + _S876.y + _S876.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S864.x + _S864.y + _S864.z);
    float _S877 = _S874.x + _S874.y + _S874.z;
    float _S878 = _S866.x + _S866.y + _S866.z;
    float _S879 = _S872.x + _S872.y + _S872.z;
    float _S880 = _S868.x + _S868.y + _S868.z;
    float _S881 = _S870.x + _S870.y + _S870.z;
    float _S882 = - s_diff_fC2_T_1;
    float3  _S883 = _S818 * _S862;
    float3  _S884 = sh_coeffs_5[int(8)] * _S862;
    float3  _S885 = _S816 * _S862;
    float3  _S886 = sh_coeffs_5[int(7)] * _S862;
    float3  _S887 = _S815 * _S862;
    float3  _S888 = sh_coeffs_5[int(6)] * _S862;
    float3  _S889 = _S817 * _S862;
    float3  _S890 = sh_coeffs_5[int(5)] * _S862;
    float3  _S891 = _S819 * _S862;
    float3  _S892 = sh_coeffs_5[int(4)] * _S862;
    float _S893 = _S890.x + _S890.y + _S890.z;
    float _S894 = _S886.x + _S886.y + _S886.z;
    float _S895 = fTmp1B_5 * _S877 + x_23 * s_diff_fS2_T_1 + y_12 * _S882 + 0.54627424478530884f * (_S892.x + _S892.y + _S892.z);
    float _S896 = fTmp1B_5 * _S878 + y_12 * s_diff_fS2_T_1 + x_23 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S884.x + _S884.y + _S884.z);
    float _S897 = y_12 * - _S896;
    float _S898 = x_23 * _S896;
    float _S899 = z_8 * (1.86588168144226074f * (z_8 * _S881) + -2.28522896766662598f * (y_12 * _S879 + x_23 * _S880) + 0.94617468118667603f * (_S888.x + _S888.y + _S888.z));
    float3  _S900 = make_float3 (0.48860251903533936f) * _S862;
    float3  _S901 = - _S900;
    float3  _S902 = _S809 * _S901;
    float3  _S903 = sh_coeffs_5[int(3)] * _S901;
    float3  _S904 = _S811 * _S900;
    float3  _S905 = sh_coeffs_5[int(2)] * _S900;
    float3  _S906 = _S813 * _S900;
    float3  _S907 = sh_coeffs_5[int(1)] * _S900;
    float _S908 = (_S820 * _S881 + 1.44530570507049561f * (fS1_5 * _S877 + fC1_5 * _S878) + -1.09254848957061768f * (y_12 * _S893 + x_23 * _S894) + _S899 + _S899 + _S905.x + _S905.y + _S905.z) / _S810;
    float _S909 = _S808 * _S908;
    float _S910 = (fTmp0C_5 * _S879 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S882 + fTmp0B_5 * _S893 + _S814 * _S895 + _S897 + _S897 + - (_S907.x + _S907.y + _S907.z)) / _S810;
    float _S911 = _S808 * _S910;
    float _S912 = (fTmp0C_5 * _S880 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S894 + 2.0f * (y_12 * _S895) + _S898 + _S898 + _S903.x + _S903.y + _S903.z) / _S810;
    float _S913 = _S808 * _S912;
    float _S914 = _S806 * - _S908 + _S805 * - _S910 + _S804 * - _S912;
    DiffPair_float_0 _S915;
    (&_S915)->primal_0 = _S807;
    (&_S915)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S915, _S914);
    float _S916 = _S806 * _S915.differential_0;
    float _S917 = _S805 * _S915.differential_0;
    float _S918 = _S804 * _S915.differential_0;
    float3  _S919 = make_float3 (0.282094806432724f) * _S862;
    float3  _S920 = make_float3 (_S913 + _S918 + _S918, _S911 + _S917 + _S917, _S909 + _S916 + _S916);
    float3  _S921 = - - _S920;
    Matrix<float, 3, 3>  _S922 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S923;
    (&_S923)->primal_0 = _S802;
    (&_S923)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S924;
    (&_S924)->primal_0 = t_5;
    (&_S924)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S923, &_S924, _S921);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S925 = _S924;
    Matrix<float, 3, 3>  _S926 = transpose_1(_S923.differential_0);
    float _S927 = _S801 * v_depth_1 + 0.001953125f * (_S800 * v_depth_1);
    float _S928 = y_11 * _S927;
    float _S929 = x_22 * _S927;
    float _S930 = z_7 * v_depth_1;
    float _S931 = z_7 * (z_7 * _S930);
    float3  _S932 = make_float3 (0.3333333432674408f) * make_float3 (_S929 + _S929, _S928 + _S928, _S799 * v_depth_1 + _S798 * _S930 + _S931 + _S931);
    DiffPair_float_0 _S933;
    (&_S933)->primal_0 = _S796;
    (&_S933)->differential_0 = 0.0f;
    DiffPair_float_0 _S934;
    (&_S934)->primal_0 = _S784;
    (&_S934)->differential_0 = 0.0f;
    _d_min_0(&_S933, &_S934, 0.0f);
    DiffPair_float_0 _S935;
    (&_S935)->primal_0 = _S736;
    (&_S935)->differential_0 = 0.0f;
    DiffPair_float_0 _S936;
    (&_S936)->primal_0 = _S760;
    (&_S936)->differential_0 = 0.0f;
    _d_min_0(&_S935, &_S936, _S933.differential_0);
    DiffPair_float_0 _S937;
    (&_S937)->primal_0 = _S795;
    (&_S937)->differential_0 = 0.0f;
    DiffPair_float_0 _S938;
    (&_S938)->primal_0 = _S784;
    (&_S938)->differential_0 = 0.0f;
    _d_max_0(&_S937, &_S938, 0.0f);
    DiffPair_float_0 _S939;
    (&_S939)->primal_0 = _S736;
    (&_S939)->differential_0 = 0.0f;
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S760;
    (&_S940)->differential_0 = 0.0f;
    _d_max_0(&_S939, &_S940, _S937.differential_0);
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = _S794;
    (&_S941)->differential_0 = 0.0f;
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = _S783;
    (&_S942)->differential_0 = 0.0f;
    _d_min_0(&_S941, &_S942, 0.0f);
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S735;
    (&_S943)->differential_0 = 0.0f;
    DiffPair_float_0 _S944;
    (&_S944)->primal_0 = _S759;
    (&_S944)->differential_0 = 0.0f;
    _d_min_0(&_S943, &_S944, _S941.differential_0);
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S793;
    (&_S945)->differential_0 = 0.0f;
    DiffPair_float_0 _S946;
    (&_S946)->primal_0 = _S783;
    (&_S946)->differential_0 = 0.0f;
    _d_max_0(&_S945, &_S946, 0.0f);
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S735;
    (&_S947)->differential_0 = 0.0f;
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = _S759;
    (&_S948)->differential_0 = 0.0f;
    _d_max_0(&_S947, &_S948, _S945.differential_0);
    DiffPair_float_0 _S949;
    (&_S949)->primal_0 = _S791;
    (&_S949)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S949, -0.0f);
    float _S950 = - (-1.0f * - (_S949.differential_0 / _S792));
    float2  _S951 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S952;
    (&_S952)->primal_0 = e2_1;
    (&_S952)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S952, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S953;
    (&_S953)->primal_0 = e1_5;
    (&_S953)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S953, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S954;
    (&_S954)->primal_0 = e0_5;
    (&_S954)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S954, 0.0f);
    DiffPair_float_0 _S955;
    (&_S955)->primal_0 = _S789;
    (&_S955)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S955, -0.0f);
    float _S956 = - _S955.differential_0;
    float2  _S957 = _S953.differential_0 + make_float2 (_S787 * _S956, _S785 * _S955.differential_0);
    float2  _S958 = - _S957;
    float2  _S959 = _S954.differential_0 + make_float2 (_S786 * _S955.differential_0, _S788 * _S956);
    float2  _S960 = - _S959;
    float2  _S961 = - _S952.differential_0 + _S957;
    float _S962 = fx_5 * (_S942.differential_0 + _S946.differential_0 + _S961.x);
    float2  _S963 = make_float2 (_S962, fy_5 * (_S934.differential_0 + _S938.differential_0 + _S961.y)) + make_float2 (dist_coeffs_5[int(8)] * _S962, dist_coeffs_5[int(9)] * _S962);
    float2  _S964 = _S773 * _S963;
    float _S965 = dist_coeffs_5[int(4)] * _S963.y;
    float _S966 = dist_coeffs_5[int(5)] * _S963.x;
    float _S967 = _S964.x + _S964.y;
    float _S968 = r2_20 * _S967;
    float _S969 = r2_20 * _S968;
    float _S970 = dist_coeffs_5[int(7)] * _S963.y + _S965 + dist_coeffs_5[int(6)] * _S963.x + _S966 + _S776 * _S967 + _S775 * _S968 + _S774 * _S969 + dist_coeffs_5[int(3)] * (r2_20 * _S969);
    float _S971 = v_20 * _S970;
    float _S972 = u_20 * _S970;
    float2  _S973 = make_float2 (radial_11) * _S963 + make_float2 (_S730 * (v_20 * _S963.y) + _S778 * _S966 + 2.0f * (u_20 * _S966) + _S727 * (v_20 * _S963.x) + _S972 + _S972, _S780 * _S965 + 2.0f * (v_20 * _S965) + _S779 * _S963.y + _S777 * _S963.x + _S971 + _S971);
    float3  _S974 = _S847.differential_0 + _S932;
    float3  _S975 = _S849 + _S850 + _S932;
    float3  _S976 = _S848.differential_0 + _S932;
    FixedArray<float3 , 2>  _S977;
    _S977[int(0)] = _S841;
    _S977[int(1)] = _S841;
    _S977[int(1)] = _S856;
    _S977[int(0)] = _S861;
    float3  _S978 = _S977[int(0)];
    float3  _S979 = _S977[int(1)];
    FixedArray<float3 , 16>  _S980;
    _S980[int(0)] = _S841;
    _S980[int(1)] = _S841;
    _S980[int(2)] = _S841;
    _S980[int(3)] = _S841;
    _S980[int(4)] = _S841;
    _S980[int(5)] = _S841;
    _S980[int(6)] = _S841;
    _S980[int(7)] = _S841;
    _S980[int(8)] = _S841;
    _S980[int(9)] = _S841;
    _S980[int(10)] = _S841;
    _S980[int(11)] = _S841;
    _S980[int(12)] = _S841;
    _S980[int(13)] = _S841;
    _S980[int(14)] = _S841;
    _S980[int(15)] = _S841;
    _S980[int(7)] = _S885;
    _S980[int(0)] = _S919;
    _S980[int(1)] = _S906;
    _S980[int(2)] = _S904;
    _S980[int(3)] = _S902;
    _S980[int(4)] = _S891;
    _S980[int(5)] = _S889;
    _S980[int(6)] = _S887;
    _S980[int(15)] = _S863;
    _S980[int(8)] = _S883;
    _S980[int(9)] = _S875;
    _S980[int(10)] = _S873;
    _S980[int(11)] = _S871;
    _S980[int(12)] = _S869;
    _S980[int(13)] = _S867;
    _S980[int(14)] = _S865;
    float3  _S981 = _S980[int(0)];
    float3  _S982 = _S980[int(1)];
    float3  _S983 = _S980[int(2)];
    float3  _S984 = _S980[int(3)];
    float3  _S985 = _S980[int(4)];
    float3  _S986 = _S980[int(5)];
    float3  _S987 = _S980[int(6)];
    float3  _S988 = _S980[int(7)];
    float3  _S989 = _S980[int(8)];
    float3  _S990 = _S980[int(9)];
    float3  _S991 = _S980[int(10)];
    float3  _S992 = _S980[int(11)];
    float3  _S993 = _S980[int(12)];
    float3  _S994 = _S980[int(13)];
    float3  _S995 = _S980[int(14)];
    float3  _S996 = _S980[int(15)];
    float _S997 = _S936.differential_0 + _S940.differential_0;
    float _S998 = _S935.differential_0 + _S939.differential_0;
    float _S999 = _S943.differential_0 + _S947.differential_0;
    float _S1000 = _S944.differential_0 + _S948.differential_0;
    float2  _S1001 = _S952.differential_0 + _S960;
    float2  _S1002 = _S958 + _S959;
    float2  _S1003 = make_float2 (0.0f, _S950);
    float2  _S1004 = _S761 * _S973;
    float2  _S1005 = _S772 * _S973;
    float _S1006 = _S1004.x + _S1004.y;
    if(_S765)
    {
        float _S1007 = _S1006 / _S766;
        float _S1008 = _S767 * - _S1007;
        float _S1009 = _S764 * (0.3333333432674408f * - (_S763 * _S1007));
        k_2 = _S1009 + _S1009;
        _S766 = _S1008;
        _S767 = 0.0f;
    }
    else
    {
        float _S1010 = _S1006 / _S768;
        float _S1011 = _S764 * - _S1010;
        k_2 = _S762 * _S1010;
        _S766 = 0.0f;
        _S767 = _S1011;
    }
    DiffPair_float_0 _S1012;
    (&_S1012)->primal_0 = _S762;
    (&_S1012)->differential_0 = 0.0f;
    DiffPair_float_0 _S1013;
    (&_S1013)->primal_0 = _S763;
    (&_S1013)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1012, &_S1013, k_2);
    float _S1014 = _S1013.differential_0 + _S766;
    float _S1015 = _S1012.differential_0 + _S767;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1016;
    (&_S1016)->primal_0 = _S761;
    (&_S1016)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S1016, _S1015);
    float2  _S1017 = _S1016.differential_0 + _S1005;
    float _S1018 = fx_5 * (_S1002.x + _S1000);
    float2  _S1019 = make_float2 (_S1018, fy_5 * (_S1002.y + _S997)) + make_float2 (dist_coeffs_5[int(8)] * _S1018, dist_coeffs_5[int(9)] * _S1018);
    float2  _S1020 = _S749 * _S1019;
    float _S1021 = dist_coeffs_5[int(4)] * _S1019.y;
    float _S1022 = dist_coeffs_5[int(5)] * _S1019.x;
    float _S1023 = _S1020.x + _S1020.y;
    float _S1024 = r2_19 * _S1023;
    float _S1025 = r2_19 * _S1024;
    float _S1026 = dist_coeffs_5[int(7)] * _S1019.y + _S1021 + dist_coeffs_5[int(6)] * _S1019.x + _S1022 + _S752 * _S1023 + _S751 * _S1024 + _S750 * _S1025 + dist_coeffs_5[int(3)] * (r2_19 * _S1025);
    float _S1027 = v_19 * _S1026;
    float _S1028 = u_19 * _S1026;
    float2  _S1029 = make_float2 (radial_10) * _S1019 + make_float2 (_S730 * (v_19 * _S1019.y) + _S754 * _S1022 + 2.0f * (u_19 * _S1022) + _S727 * (v_19 * _S1019.x) + _S1028 + _S1028, _S756 * _S1021 + 2.0f * (v_19 * _S1021) + _S755 * _S1019.y + _S753 * _S1019.x + _S1027 + _S1027);
    float3  _S1030 = _S976 + make_float3 (_S1017.x, _S1017.y, _S1014);
    float2  _S1031 = _S737 * _S1029;
    float2  _S1032 = _S748 * _S1029;
    float _S1033 = _S1031.x + _S1031.y;
    if(_S741)
    {
        float _S1034 = _S1033 / _S742;
        float _S1035 = _S743 * - _S1034;
        float _S1036 = _S740 * (0.3333333432674408f * - (_S739 * _S1034));
        k_2 = _S1036 + _S1036;
        _S742 = _S1035;
        _S743 = 0.0f;
    }
    else
    {
        float _S1037 = _S1033 / _S744;
        float _S1038 = _S740 * - _S1037;
        k_2 = _S738 * _S1037;
        _S742 = 0.0f;
        _S743 = _S1038;
    }
    DiffPair_float_0 _S1039;
    (&_S1039)->primal_0 = _S738;
    (&_S1039)->differential_0 = 0.0f;
    DiffPair_float_0 _S1040;
    (&_S1040)->primal_0 = _S739;
    (&_S1040)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1039, &_S1040, k_2);
    float _S1041 = _S1040.differential_0 + _S742;
    float _S1042 = _S1039.differential_0 + _S743;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1043;
    (&_S1043)->primal_0 = _S737;
    (&_S1043)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S1043, _S1042);
    float2  _S1044 = _S1043.differential_0 + _S1032;
    float _S1045 = fx_5 * (_S1001.x + _S999);
    float2  _S1046 = make_float2 (_S1045, fy_5 * (_S1001.y + _S998)) + make_float2 (dist_coeffs_5[int(8)] * _S1045, dist_coeffs_5[int(9)] * _S1045);
    float2  _S1047 = _S723 * _S1046;
    float _S1048 = dist_coeffs_5[int(4)] * _S1046.y;
    float _S1049 = dist_coeffs_5[int(5)] * _S1046.x;
    float _S1050 = _S1047.x + _S1047.y;
    float _S1051 = r2_18 * _S1050;
    float _S1052 = r2_18 * _S1051;
    float _S1053 = dist_coeffs_5[int(7)] * _S1046.y + _S1048 + dist_coeffs_5[int(6)] * _S1046.x + _S1049 + _S726 * _S1050 + _S725 * _S1051 + _S724 * _S1052 + dist_coeffs_5[int(3)] * (r2_18 * _S1052);
    float _S1054 = v_18 * _S1053;
    float _S1055 = u_18 * _S1053;
    float2  _S1056 = make_float2 (radial_9) * _S1046 + make_float2 (_S730 * (v_18 * _S1046.y) + _S729 * _S1049 + 2.0f * (u_18 * _S1049) + _S727 * (v_18 * _S1046.x) + _S1055 + _S1055, _S732 * _S1048 + 2.0f * (v_18 * _S1048) + _S731 * _S1046.y + _S728 * _S1046.x + _S1054 + _S1054);
    float3  _S1057 = _S974 + make_float3 (_S1044.x, _S1044.y, _S1041);
    float2  _S1058 = _S711 * _S1056;
    float2  _S1059 = _S722 * _S1056;
    float _S1060 = _S1058.x + _S1058.y;
    if(_S715)
    {
        float _S1061 = _S1060 / _S716;
        float _S1062 = _S717 * - _S1061;
        float _S1063 = _S714 * (0.3333333432674408f * - (_S713 * _S1061));
        k_2 = _S1063 + _S1063;
        _S716 = _S1062;
        _S717 = 0.0f;
    }
    else
    {
        float _S1064 = _S1060 / _S718;
        float _S1065 = _S714 * - _S1064;
        k_2 = _S712 * _S1064;
        _S716 = 0.0f;
        _S717 = _S1065;
    }
    DiffPair_float_0 _S1066;
    (&_S1066)->primal_0 = _S712;
    (&_S1066)->differential_0 = 0.0f;
    DiffPair_float_0 _S1067;
    (&_S1067)->primal_0 = _S713;
    (&_S1067)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1066, &_S1067, k_2);
    float _S1068 = _S1067.differential_0 + _S716;
    float _S1069 = _S1066.differential_0 + _S717;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1070;
    (&_S1070)->primal_0 = _S711;
    (&_S1070)->differential_0 = _S951;
    s_bwd_length_impl_1(&_S1070, _S1069);
    float2  _S1071 = _S1070.differential_0 + _S1059;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1072;
    (&_S1072)->primal_0 = vert2_c_5;
    (&_S1072)->differential_0 = _S841;
    s_bwd_length_impl_0(&_S1072, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1073;
    (&_S1073)->primal_0 = vert1_c_5;
    (&_S1073)->differential_0 = _S841;
    s_bwd_length_impl_0(&_S1073, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1074;
    (&_S1074)->primal_0 = vert0_c_5;
    (&_S1074)->differential_0 = _S841;
    s_bwd_length_impl_0(&_S1074, 0.0f);
    float3  _S1075 = _S1072.differential_0 + _S1030;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1076;
    (&_S1076)->primal_0 = R_5;
    (&_S1076)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1077;
    (&_S1077)->primal_0 = vert2_5;
    (&_S1077)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1076, &_S1077, _S1075);
    float3  _S1078 = _S1073.differential_0 + _S1057;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1079;
    (&_S1079)->primal_0 = R_5;
    (&_S1079)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1080;
    (&_S1080)->primal_0 = vert1_5;
    (&_S1080)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1079, &_S1080, _S1078);
    float3  _S1081 = _S1074.differential_0 + _S975 + make_float3 (_S1071.x, _S1071.y, _S1068);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1082;
    (&_S1082)->primal_0 = R_5;
    (&_S1082)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1083;
    (&_S1083)->primal_0 = vert0_5;
    (&_S1083)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1082, &_S1083, _S1081);
    float3  _S1084 = _S1077.differential_0 + v_verts_1[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1085;
    (&_S1085)->primal_0 = _S705;
    (&_S1085)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1086;
    (&_S1086)->primal_0 = _S710;
    (&_S1086)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1085, &_S1086, _S1084);
    float _S1087 = - _S1086.differential_0.y;
    float _S1088 = _S709 * _S1086.differential_0.x;
    float _S1089 = - (_S700 * _S1086.differential_0.x);
    float3  _S1090 = _S1080.differential_0 + v_verts_1[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1091;
    (&_S1091)->primal_0 = _S705;
    (&_S1091)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1092;
    (&_S1092)->primal_0 = _S708;
    (&_S1092)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1091, &_S1092, _S1090);
    float _S1093 = _S700 * _S1092.differential_0.x;
    float _S1094 = _S707 * _S1092.differential_0.x;
    float3  _S1095 = _S1083.differential_0 + v_verts_1[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1096;
    (&_S1096)->primal_0 = _S705;
    (&_S1096)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1097;
    (&_S1097)->primal_0 = _S706;
    (&_S1097)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1096, &_S1097, _S1095);
    Matrix<float, 3, 3>  _S1098 = transpose_1(_S1085.differential_0 + _S1091.differential_0 + _S1096.differential_0);
    float _S1099 = 2.0f * - _S1098.rows[int(2)].z;
    float _S1100 = 2.0f * _S1098.rows[int(2)].y;
    float _S1101 = 2.0f * _S1098.rows[int(2)].x;
    float _S1102 = 2.0f * _S1098.rows[int(1)].z;
    float _S1103 = 2.0f * - _S1098.rows[int(1)].y;
    float _S1104 = 2.0f * _S1098.rows[int(1)].x;
    float _S1105 = 2.0f * _S1098.rows[int(0)].z;
    float _S1106 = 2.0f * _S1098.rows[int(0)].y;
    float _S1107 = 2.0f * - _S1098.rows[int(0)].x;
    float _S1108 = - _S1104 + _S1106;
    float _S1109 = _S1101 + - _S1105;
    float _S1110 = - _S1100 + _S1102;
    float _S1111 = _S1100 + _S1102;
    float _S1112 = _S1101 + _S1105;
    float _S1113 = _S1104 + _S1106;
    float _S1114 = _S703.w * (_S1103 + _S1107);
    float _S1115 = _S703.z * (_S1099 + _S1107);
    float _S1116 = _S703.y * (_S1099 + _S1103);
    float _S1117 = _S703.x * _S1108 + _S703.z * _S1111 + _S703.y * _S1112 + _S1114 + _S1114;
    float _S1118 = _S703.x * _S1109 + _S703.w * _S1111 + _S703.y * _S1113 + _S1115 + _S1115;
    float _S1119 = _S703.x * _S1110 + _S703.w * _S1112 + _S703.z * _S1113 + _S1116 + _S1116;
    float _S1120 = _S703.w * _S1108 + _S703.z * _S1109 + _S703.y * _S1110;
    float4  _S1121 = make_float4 (0.0f);
    float4  _S1122 = _S1121;
    *&((&_S1122)->w) = _S1117;
    *&((&_S1122)->z) = _S1118;
    *&((&_S1122)->y) = _S1119;
    *&((&_S1122)->x) = _S1120;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1123;
    (&_S1123)->primal_0 = quat_5;
    (&_S1123)->differential_0 = _S1121;
    s_bwd_normalize_impl_1(&_S1123, _S1122);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1124 = _S1123;
    float _S1125 = _S1089 + _S1093;
    float _S1126 = 0.5f * - _S1125;
    float _S1127 = _S1087 + _S1092.differential_0.y;
    DiffPair_float_0 _S1128;
    (&_S1128)->primal_0 = _S701;
    (&_S1128)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1128, _S1127);
    float _S1129 = _S1126 + _S1128.differential_0;
    float _S1130 = _S1088 + _S1094 + _S1097.differential_0.x;
    DiffPair_float_0 _S1131;
    (&_S1131)->primal_0 = _S699;
    (&_S1131)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1131, _S1130);
    float _S1132 = _S1126 + _S1131.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1133;
    (&_S1133)->primal_0 = mean_c_5;
    (&_S1133)->differential_0 = _S841;
    s_bwd_length_impl_0(&_S1133, 0.0f);
    float3  _S1134 = _S1133.differential_0 + _S844.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1135;
    (&_S1135)->primal_0 = R_5;
    (&_S1135)->differential_0 = _S922;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1136;
    (&_S1136)->primal_0 = mean_5;
    (&_S1136)->differential_0 = _S841;
    s_bwd_prop_mul_0(&_S1135, &_S1136, _S1134);
    float3  _S1137 = _S1075 + _S1078 + _S1081 + _S1134 + _S925.differential_0;
    Matrix<float, 3, 3>  _S1138 = _S1076.differential_0 + _S1079.differential_0 + _S1082.differential_0 + _S1135.differential_0 + _S926;
    float3  _S1139 = make_float3 (_S1132, _S1129, _S1125);
    float3  _S1140 = _S1084 + _S1090 + _S1095 + _S1136.differential_0 + _S920;
    *v_mean_1 = _S1140;
    *v_quat_1 = _S1124.differential_0;
    *v_scale_1 = _S1139;
    *v_hardness_1 = _S1003;
    (*v_sh_coeffs_1)[int(0)] = _S981;
    (*v_sh_coeffs_1)[int(1)] = _S982;
    (*v_sh_coeffs_1)[int(2)] = _S983;
    (*v_sh_coeffs_1)[int(3)] = _S984;
    (*v_sh_coeffs_1)[int(4)] = _S985;
    (*v_sh_coeffs_1)[int(5)] = _S986;
    (*v_sh_coeffs_1)[int(6)] = _S987;
    (*v_sh_coeffs_1)[int(7)] = _S988;
    (*v_sh_coeffs_1)[int(8)] = _S989;
    (*v_sh_coeffs_1)[int(9)] = _S990;
    (*v_sh_coeffs_1)[int(10)] = _S991;
    (*v_sh_coeffs_1)[int(11)] = _S992;
    (*v_sh_coeffs_1)[int(12)] = _S993;
    (*v_sh_coeffs_1)[int(13)] = _S994;
    (*v_sh_coeffs_1)[int(14)] = _S995;
    (*v_sh_coeffs_1)[int(15)] = _S996;
    (*v_ch_coeffs_1)[int(0)] = _S978;
    (*v_ch_coeffs_1)[int(1)] = _S979;
    *v_R_1 = _S1138;
    *v_t_1 = _S1137;
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
    bool _S1141;
    if((*u_21) >= 0.0f)
    {
        _S1141 = (*v_21) >= 0.0f;
    }
    else
    {
        _S1141 = false;
    }
    if(_S1141)
    {
        _S1141 = (*u_21 + *v_21) <= 1.0f;
    }
    else
    {
        _S1141 = false;
    }
    if(_S1141)
    {
        _S1141 = (*t_6) >= 0.0f;
    }
    else
    {
        _S1141 = false;
    }
    return _S1141;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_14, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S1142 = *dpx_14;
    bool _S1143;
    if(((*dpx_14).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S1143 = ((*dpx_14).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S1143 = false;
    }
    float _S1144;
    if(_S1143)
    {
        _S1144 = dOut_11;
    }
    else
    {
        _S1144 = 0.0f;
    }
    dpx_14->primal_0 = _S1142.primal_0;
    dpx_14->differential_0 = _S1144;
    DiffPair_float_0 _S1145 = *dpMin_0;
    if((_S1142.primal_0) < ((*dpMin_0).primal_0))
    {
        _S1144 = dOut_11;
    }
    else
    {
        _S1144 = 0.0f;
    }
    dpMin_0->primal_0 = _S1145.primal_0;
    dpMin_0->differential_0 = _S1144;
    DiffPair_float_0 _S1146 = *dpMax_0;
    if(((*dpx_14).primal_0) > ((*dpMax_0).primal_0))
    {
        _S1144 = dOut_11;
    }
    else
    {
        _S1144 = 0.0f;
    }
    dpMax_0->primal_0 = _S1146.primal_0;
    dpMax_0->differential_0 = _S1144;
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
        DiffPair_float_0 _S1147 = *dpx_15;
        float _S1148 = val_0 * (*dpy_5).primal_0 / (*dpx_15).primal_0 * dOut_12;
        dpx_15->primal_0 = (*dpx_15).primal_0;
        dpx_15->differential_0 = _S1148;
        float _S1149 = val_0 * (F32_log((_S1147.primal_0))) * dOut_12;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1149;
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
    bool _S1150;
    if(u_22 >= 0.0f)
    {
        _S1150 = v_22 >= 0.0f;
    }
    else
    {
        _S1150 = false;
    }
    if(_S1150)
    {
        _S1150 = (u_22 + v_22) <= 1.0f;
    }
    else
    {
        _S1150 = false;
    }
    if(_S1150)
    {
        _S1150 = t_7 >= 0.0f;
    }
    else
    {
        _S1150 = false;
    }
    if(!_S1150)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_22), (v_22)))), ((F32_sqrt((0.5f))) * (1.0f - u_22 - v_22)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_6.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_6.x;
    float _S1151;
    if(opac_0 < 0.0f)
    {
        _S1151 = 0.0f;
    }
    else
    {
        _S1151 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S1151;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ float s_primal_ctx_clamp_0(float _S1152, float _S1153, float _S1154)
{
    return clamp_0(_S1152, _S1153, _S1154);
}

inline __device__ float s_primal_ctx_pow_0(float _S1155, float _S1156)
{
    return (F32_pow((_S1155), (_S1156)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1157, DiffPair_float_0 * _S1158, float _S1159)
{
    _d_pow_0(_S1157, _S1158, _S1159);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1160, DiffPair_float_0 * _S1161, DiffPair_float_0 * _S1162, float _S1163)
{
    _d_clamp_0(_S1160, _S1161, _S1162, _S1163);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_5)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1164 = *dphardness_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1165 = *dpray_d_0;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_0).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S1166 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S1167 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_0).primal_0);
    float _S1168 = s_primal_ctx_dot_0((*dpray_d_0).primal_0, _S1166);
    float d_2 = 1.0f / _S1168;
    float _S1169 = _S1168 * _S1168;
    float3  _S1170 = - _S1167;
    float _S1171 = s_primal_ctx_dot_0(_S1170, v2v0_2);
    float u_23 = d_2 * _S1171;
    float _S1172 = s_primal_ctx_dot_0(_S1167, v1v0_2);
    float v_23 = d_2 * _S1172;
    float3  _S1173 = - _S1166;
    float t_8 = d_2 * s_primal_ctx_dot_0(_S1173, rov0_2);
    bool _S1174;
    if(u_23 >= 0.0f)
    {
        _S1174 = v_23 >= 0.0f;
    }
    else
    {
        _S1174 = false;
    }
    if(_S1174)
    {
        _S1174 = (u_23 + v_23) <= 1.0f;
    }
    else
    {
        _S1174 = false;
    }
    if(_S1174)
    {
        _S1174 = t_8 >= 0.0f;
    }
    else
    {
        _S1174 = false;
    }
    bool _S1175 = !!_S1174;
    float _S1176;
    float _S1177;
    float _S1178;
    float _S1179;
    float _S1180;
    float _S1181;
    float _S1182;
    float _S1183;
    float _S1184;
    float _S1185;
    float _S1186;
    if(_S1175)
    {
        float _S1187 = (F32_min((u_23), (v_23)));
        float _S1188 = s_primal_ctx_sqrt_0(0.5f);
        float _S1189 = _S1188 * (1.0f - u_23 - v_23);
        float _S1190 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = (F32_min((_S1187), (_S1189))) * _S1190;
        float _S1191 = _S1164.primal_0.y;
        float _S1192 = 1.0f - opac_1;
        float _S1193 = 1.0f - s_primal_ctx_clamp_0(_S1191, 0.0f, 0.99989998340606689f);
        float _S1194 = 1.0f / _S1193;
        float _S1195 = _S1193 * _S1193;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S1192, _S1194);
        float o_1 = _S1164.primal_0.x;
        bool _S1196 = opac_1 < 0.0f;
        if(_S1196)
        {
            _S1176 = 0.0f;
        }
        else
        {
            _S1176 = o_1 * w_1;
        }
        _S1174 = _S1196;
        _S1177 = o_1;
        _S1178 = w_1;
        _S1179 = _S1192;
        _S1180 = _S1194;
        _S1181 = _S1195;
        _S1182 = _S1191;
        _S1183 = _S1190;
        _S1184 = _S1187;
        _S1185 = _S1189;
        _S1186 = _S1188;
    }
    else
    {
        _S1174 = false;
        _S1176 = 0.0f;
        _S1177 = 0.0f;
        _S1178 = 0.0f;
        _S1179 = 0.0f;
        _S1180 = 0.0f;
        _S1181 = 0.0f;
        _S1182 = 0.0f;
        _S1183 = 0.0f;
        _S1184 = 0.0f;
        _S1185 = 0.0f;
        _S1186 = 0.0f;
    }
    float2  _S1197 = make_float2 (0.0f);
    float2  _S1198;
    if(_S1175)
    {
        if(_S1174)
        {
            _S1176 = 0.0f;
            _S1177 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S1199;
            (&_S1199)->primal_0 = _S1176;
            (&_S1199)->differential_0 = 0.0f;
            DiffPair_float_0 _S1200;
            (&_S1200)->primal_0 = 0.99500000476837158f;
            (&_S1200)->differential_0 = 0.0f;
            _d_min_0(&_S1199, &_S1200, _s_dOut_5);
            float _S1201 = _S1177 * _S1199.differential_0;
            _S1176 = _S1178 * _S1199.differential_0;
            _S1177 = _S1201;
        }
        float _S1202 = - _S1177;
        DiffPair_float_0 _S1203;
        (&_S1203)->primal_0 = _S1179;
        (&_S1203)->differential_0 = 0.0f;
        DiffPair_float_0 _S1204;
        (&_S1204)->primal_0 = _S1180;
        (&_S1204)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1203, &_S1204, _S1202);
        float _S1205 = - - (_S1204.differential_0 / _S1181);
        float s_diff_opac_T_0 = - _S1203.differential_0;
        DiffPair_float_0 _S1206;
        (&_S1206)->primal_0 = _S1182;
        (&_S1206)->differential_0 = 0.0f;
        DiffPair_float_0 _S1207;
        (&_S1207)->primal_0 = 0.0f;
        (&_S1207)->differential_0 = 0.0f;
        DiffPair_float_0 _S1208;
        (&_S1208)->primal_0 = 0.99989998340606689f;
        (&_S1208)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1206, &_S1207, &_S1208, _S1205);
        float _S1209 = _S1183 * s_diff_opac_T_0;
        DiffPair_float_0 _S1210;
        (&_S1210)->primal_0 = _S1184;
        (&_S1210)->differential_0 = 0.0f;
        DiffPair_float_0 _S1211;
        (&_S1211)->primal_0 = _S1185;
        (&_S1211)->differential_0 = 0.0f;
        _d_min_0(&_S1210, &_S1211, _S1209);
        float _S1212 = - (_S1186 * _S1211.differential_0);
        DiffPair_float_0 _S1213;
        (&_S1213)->primal_0 = u_23;
        (&_S1213)->differential_0 = 0.0f;
        DiffPair_float_0 _S1214;
        (&_S1214)->primal_0 = v_23;
        (&_S1214)->differential_0 = 0.0f;
        _d_min_0(&_S1213, &_S1214, _S1210.differential_0);
        float2  _S1215 = make_float2 (_S1176, _S1206.differential_0);
        float _S1216 = _S1212 + _S1214.differential_0;
        _S1176 = _S1212 + _S1213.differential_0;
        _S1177 = _S1216;
        _S1198 = _S1215;
    }
    else
    {
        _S1176 = 0.0f;
        _S1177 = 0.0f;
        _S1198 = _S1197;
    }
    float3  _S1217 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1218;
    (&_S1218)->primal_0 = _S1173;
    (&_S1218)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1219;
    (&_S1219)->primal_0 = rov0_2;
    (&_S1219)->differential_0 = _S1217;
    s_bwd_prop_dot_0(&_S1218, &_S1219, 0.0f);
    float3  _S1220 = - _S1218.differential_0;
    float _S1221 = d_2 * _S1177;
    float _S1222 = _S1172 * _S1177;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1223;
    (&_S1223)->primal_0 = _S1167;
    (&_S1223)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1224;
    (&_S1224)->primal_0 = v1v0_2;
    (&_S1224)->differential_0 = _S1217;
    s_bwd_prop_dot_0(&_S1223, &_S1224, _S1221);
    float _S1225 = d_2 * _S1176;
    float _S1226 = _S1171 * _S1176;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1227;
    (&_S1227)->primal_0 = _S1170;
    (&_S1227)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1228;
    (&_S1228)->primal_0 = v2v0_2;
    (&_S1228)->differential_0 = _S1217;
    s_bwd_prop_dot_0(&_S1227, &_S1228, _S1225);
    float3  _S1229 = - _S1227.differential_0;
    float _S1230 = - ((_S1222 + _S1226) / _S1169);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1231;
    (&_S1231)->primal_0 = _S1165.primal_0;
    (&_S1231)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1232;
    (&_S1232)->primal_0 = _S1166;
    (&_S1232)->differential_0 = _S1217;
    s_bwd_prop_dot_0(&_S1231, &_S1232, _S1230);
    float3  _S1233 = _S1223.differential_0 + _S1229;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
    (&_S1234)->primal_0 = rov0_2;
    (&_S1234)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1235;
    (&_S1235)->primal_0 = _S1165.primal_0;
    (&_S1235)->differential_0 = _S1217;
    s_bwd_prop_cross_0(&_S1234, &_S1235, _S1233);
    float3  _S1236 = _S1220 + _S1232.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1237;
    (&_S1237)->primal_0 = v1v0_2;
    (&_S1237)->differential_0 = _S1217;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1238;
    (&_S1238)->primal_0 = v2v0_2;
    (&_S1238)->differential_0 = _S1217;
    s_bwd_prop_cross_0(&_S1237, &_S1238, _S1236);
    float3  _S1239 = _S1219.differential_0 + _S1234.differential_0;
    float3  _S1240 = _S1228.differential_0 + _S1238.differential_0;
    float3  _S1241 = _S1224.differential_0 + _S1237.differential_0;
    float3  _S1242 = - _S1239 + - _S1240 + - _S1241;
    float3  _S1243 = _S1231.differential_0 + _S1235.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S1243;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S1239;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S1198;
    FixedArray<float3 , 3>  _S1244;
    _S1244[int(0)] = _S1217;
    _S1244[int(1)] = _S1217;
    _S1244[int(2)] = _S1217;
    _S1244[int(2)] = _S1240;
    _S1244[int(0)] = _S1242;
    _S1244[int(1)] = _S1241;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S1244;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1245, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1246, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1247, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1248, float _S1249)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S1245, _S1246, _S1247, _S1248, _S1249);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_6, float2  hardness_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    float3  _S1250 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1251 = { _S1250, _S1250, _S1250 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = verts_6;
    (&dp_verts_0)->differential_0 = _S1251;
    float2  _S1252 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S1252;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1250;
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
    float3  _S1253 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S1254 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_1).primal_0);
    float _S1255 = s_primal_ctx_dot_0((*dpray_d_1).primal_0, _S1253);
    float d_4 = 1.0f / _S1255;
    float _S1256 = _S1255 * _S1255;
    float3  _S1257 = - _S1254;
    float _S1258 = s_primal_ctx_dot_0(_S1257, v2v0_4);
    float u_25 = d_4 * _S1258;
    float _S1259 = s_primal_ctx_dot_0(_S1254, v1v0_4);
    float v_25 = d_4 * _S1259;
    float3  _S1260 = - _S1253;
    float3  _S1261 = dprgbs_0->primal_0[int(2)] * dpcolor_0;
    float3  _S1262 = make_float3 (v_25) * dpcolor_0;
    float3  _S1263 = dprgbs_0->primal_0[int(1)] * dpcolor_0;
    float3  _S1264 = make_float3 (u_25) * dpcolor_0;
    float3  _S1265 = dprgbs_0->primal_0[int(0)] * dpcolor_0;
    float3  _S1266 = make_float3 (1.0f - u_25 - v_25) * dpcolor_0;
    float _S1267 = - (_S1265.x + _S1265.y + _S1265.z);
    float _S1268 = d_4 * dpdepth_0;
    float _S1269 = s_primal_ctx_dot_0(_S1260, rov0_4) * dpdepth_0;
    float3  _S1270 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1271;
    (&_S1271)->primal_0 = _S1260;
    (&_S1271)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1272;
    (&_S1272)->primal_0 = rov0_4;
    (&_S1272)->differential_0 = _S1270;
    s_bwd_prop_dot_0(&_S1271, &_S1272, _S1268);
    float3  _S1273 = - _S1271.differential_0;
    float _S1274 = _S1267 + _S1261.x + _S1261.y + _S1261.z;
    float _S1275 = d_4 * _S1274;
    float _S1276 = _S1259 * _S1274;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1277;
    (&_S1277)->primal_0 = _S1254;
    (&_S1277)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1278;
    (&_S1278)->primal_0 = v1v0_4;
    (&_S1278)->differential_0 = _S1270;
    s_bwd_prop_dot_0(&_S1277, &_S1278, _S1275);
    float _S1279 = _S1267 + _S1263.x + _S1263.y + _S1263.z;
    float _S1280 = d_4 * _S1279;
    float _S1281 = _S1258 * _S1279;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1282;
    (&_S1282)->primal_0 = _S1257;
    (&_S1282)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1283;
    (&_S1283)->primal_0 = v2v0_4;
    (&_S1283)->differential_0 = _S1270;
    s_bwd_prop_dot_0(&_S1282, &_S1283, _S1280);
    float3  _S1284 = - _S1282.differential_0;
    float _S1285 = - ((_S1269 + _S1276 + _S1281) / _S1256);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1286;
    (&_S1286)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1286)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1287;
    (&_S1287)->primal_0 = _S1253;
    (&_S1287)->differential_0 = _S1270;
    s_bwd_prop_dot_0(&_S1286, &_S1287, _S1285);
    float3  _S1288 = _S1277.differential_0 + _S1284;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1289;
    (&_S1289)->primal_0 = rov0_4;
    (&_S1289)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1290;
    (&_S1290)->primal_0 = (*dpray_d_1).primal_0;
    (&_S1290)->differential_0 = _S1270;
    s_bwd_prop_cross_0(&_S1289, &_S1290, _S1288);
    float3  _S1291 = _S1273 + _S1287.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1292;
    (&_S1292)->primal_0 = v1v0_4;
    (&_S1292)->differential_0 = _S1270;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1293;
    (&_S1293)->primal_0 = v2v0_4;
    (&_S1293)->differential_0 = _S1270;
    s_bwd_prop_cross_0(&_S1292, &_S1293, _S1291);
    float3  _S1294 = _S1272.differential_0 + _S1289.differential_0;
    float3  _S1295 = _S1283.differential_0 + _S1293.differential_0;
    float3  _S1296 = _S1278.differential_0 + _S1292.differential_0;
    float3  _S1297 = - _S1294 + - _S1295 + - _S1296;
    float3  _S1298 = _S1286.differential_0 + _S1290.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1298;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1294;
    FixedArray<float3 , 3>  _S1299;
    _S1299[int(0)] = _S1270;
    _S1299[int(1)] = _S1270;
    _S1299[int(2)] = _S1270;
    _S1299[int(2)] = _S1262;
    _S1299[int(1)] = _S1264;
    _S1299[int(0)] = _S1266;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S1299;
    FixedArray<float3 , 3>  _S1300;
    _S1300[int(0)] = _S1270;
    _S1300[int(1)] = _S1270;
    _S1300[int(2)] = _S1270;
    _S1300[int(2)] = _S1295;
    _S1300[int(0)] = _S1297;
    _S1300[int(1)] = _S1296;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S1300;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1301, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S1302, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1303, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1304, float3  _S1305, float _S1306)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S1301, _S1302, _S1303, _S1304, _S1305, _S1306);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  verts_9, FixedArray<float3 , 3>  rgbs_6, float3  ray_o_5, float3  ray_d_5, float3  v_color_0, float v_depth_2, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1307 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S1308 = { _S1307, _S1307, _S1307 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = verts_9;
    (&dp_verts_1)->differential_0 = _S1308;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S1308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_5;
    (&dp_ray_o_1)->differential_0 = _S1307;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_5;
    (&dp_ray_d_1)->differential_0 = _S1307;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_1, &dp_ray_d_1, v_color_0, v_depth_2);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

