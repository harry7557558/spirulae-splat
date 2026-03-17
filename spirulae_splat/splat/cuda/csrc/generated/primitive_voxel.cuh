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

inline __device__ float dot_0(float2  x_0, float2  y_0)
{
    int i_0 = int(0);
    float result_0 = 0.0f;
    for(;;)
    {
        if(i_0 < int(2))
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

inline __device__ float dot_1(float3  x_1, float3  y_1)
{
    int i_1 = int(0);
    float result_2 = 0.0f;
    for(;;)
    {
        if(i_1 < int(3))
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
    return (F32_sqrt((dot_0(x_2, x_2))));
}

inline __device__ float length_1(float3  x_3)
{
    return (F32_sqrt((dot_1(x_3, x_3))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_1, DiffPair_float_0 * dpx_2, float dOut_2)
{
    DiffPair_float_0 _S5 = *dpx_2;
    float _S6 = - (*dpy_1).primal_0 / ((*dpx_2).primal_0 * (*dpx_2).primal_0 + (*dpy_1).primal_0 * (*dpy_1).primal_0) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S6;
    float _S7 = _S5.primal_0 / (_S5.primal_0 * _S5.primal_0 + (*dpy_1).primal_0 * (*dpy_1).primal_0) * dOut_2;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = _S7;
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

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_2, float dOut_3)
{
    DiffPair_float_0 _S8 = *dpx_3;
    float _S9;
    if(((*dpx_3).primal_0) < ((*dpy_2).primal_0))
    {
        _S9 = dOut_3;
    }
    else
    {
        if(((*dpx_3).primal_0) > ((*dpy_2).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_3;
        }
    }
    dpx_3->primal_0 = _S8.primal_0;
    dpx_3->differential_0 = _S9;
    DiffPair_float_0 _S10 = *dpy_2;
    if(((*dpy_2).primal_0) < (_S8.primal_0))
    {
        _S9 = dOut_3;
    }
    else
    {
        if(((*dpy_2).primal_0) > ((*dpx_3).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_3;
        }
    }
    dpy_2->primal_0 = _S10.primal_0;
    dpy_2->differential_0 = _S9;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_4)
{
    float _S11 = (*left_0).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_4.x;
    float sum_0 = _S11 + (*left_0).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_4.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_4.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S12 = (*left_0).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_4.x;
    float sum_2 = _S12 + (*left_0).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_4.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_4.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S13 = (*left_0).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_4.x;
    float sum_4 = _S13 + (*left_0).primal_0.rows[int(1)].z * dOut_4.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_4.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_4.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_4.z;
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

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_4, float dOut_5)
{
    float _S14 = (F32_exp(((*dpx_4).primal_0))) * dOut_5;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S14;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_3, float3  dOut_6)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_5).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_3).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_6.x);
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = left_dp_0.differential_0;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_5).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_3).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_6.y);
    *&((&left_d_result_1)->y) = left_dp_1.differential_0;
    *&((&right_d_result_1)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_5).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_3).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_6.z);
    *&((&left_d_result_1)->z) = left_dp_2.differential_0;
    *&((&right_d_result_1)->z) = right_dp_2.differential_0;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = left_d_result_1;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = right_d_result_1;
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

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  densities_0, FixedArray<float3 , 16>  sh_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, float4  * aabb_xyxy_0, float * depth_0, float3  * rgbs_0)
{
    float2  * _S15;
    float2  * _S16;
    float2  * _S17;
    float2  * _S18;
    float2  * _S19;
    float2  * _S20;
    float2  * _S21;
    float2  * _S22;
    bool _S23;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_0;
        float3  _S24 = mul_0(R_0, pos_0) + t_0;
        pos_c_0[int(0)] = _S24;
        float _S25 = _S24.z;
        float _S26 = (F32_min((far_plane_0), (_S25)));
        float _S27 = (F32_max((near_plane_0), (_S25)));
        float3  _S28 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 0.0f)) + t_0;
        pos_c_0[int(1)] = _S28;
        float _S29 = _S28.z;
        float _S30 = (F32_min((_S26), (_S29)));
        float _S31 = (F32_max((_S27), (_S29)));
        float3  _S32 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 0.0f)) + t_0;
        pos_c_0[int(2)] = _S32;
        float _S33 = _S32.z;
        float _S34 = (F32_min((_S30), (_S33)));
        float _S35 = (F32_max((_S31), (_S33)));
        float3  _S36 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 1.0f, 0.0f)) + t_0;
        pos_c_0[int(3)] = _S36;
        float _S37 = _S36.z;
        float _S38 = (F32_min((_S34), (_S37)));
        float _S39 = (F32_max((_S35), (_S37)));
        float3  _S40 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 0.0f, 1.0f)) + t_0;
        pos_c_0[int(4)] = _S40;
        float _S41 = _S40.z;
        float _S42 = (F32_min((_S38), (_S41)));
        float _S43 = (F32_max((_S39), (_S41)));
        float3  _S44 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 1.0f)) + t_0;
        pos_c_0[int(5)] = _S44;
        float _S45 = _S44.z;
        float _S46 = (F32_min((_S42), (_S45)));
        float _S47 = (F32_max((_S43), (_S45)));
        float3  _S48 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 1.0f)) + t_0;
        pos_c_0[int(6)] = _S48;
        float _S49 = _S48.z;
        float _S50 = (F32_min((_S46), (_S49)));
        float _S51 = (F32_max((_S47), (_S49)));
        float3  _S52 = mul_0(R_0, pos_0 + make_float3 (size_0)) + t_0;
        pos_c_0[int(7)] = _S52;
        float _S53 = _S52.z;
        float _S54 = (F32_min((_S50), (_S53)));
        float _S55 = (F32_max((_S51), (_S53)));
        bool _S56;
        if(_S54 < near_plane_0)
        {
            _S56 = true;
        }
        else
        {
            _S56 = _S55 > far_plane_0;
        }
        if(_S56)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float3  mean_c_0 = mul_0(R_0, pos_0 + make_float3 (0.5f * size_0)) + t_0;
        FixedArray<float2 , 8>  uv_0;
        for(;;)
        {
            float3  _S57 = pos_c_0[int(0)];
            _S15 = &uv_0[int(0)];
            for(;;)
            {
                float _S58 = _S57.z;
                uv_0[int(0)] = float2 {_S57.x, _S57.y} / make_float2 (_S58);
                if(_S58 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_0 = uv_0[int(0)].x;
                    float v_0 = uv_0[int(0)].y;
                    float _S59 = u_0 + u_0;
                    float r2_0 = u_0 * u_0 + v_0 * v_0;
                    float _S60 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                    float _S61 = dist_coeffs_0[int(1)] + r2_0 * _S60;
                    float _S62 = dist_coeffs_0[int(0)] + r2_0 * _S61;
                    float radial_0 = 1.0f + r2_0 * _S62;
                    float _S63 = 2.0f * dist_coeffs_0[int(4)];
                    float _S64 = 2.0f * u_0;
                    float _S65 = 2.0f * dist_coeffs_0[int(5)];
                    float _S66 = 2.0f * v_0;
                    float2  _S67 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (_S59 * _S62 + (_S59 * _S61 + (_S59 * _S60 + _S59 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0[int(0)] + make_float2 (_S63 * v_0 + (_S59 + (_S64 + _S64)) * dist_coeffs_0[int(5)] + _S59 * dist_coeffs_0[int(6)], _S65 * v_0 + _S59 * dist_coeffs_0[int(4)] + _S59 * dist_coeffs_0[int(7)]);
                    float _S68 = v_0 + v_0;
                    float2  _S69 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (_S68 * _S62 + (_S68 * _S61 + (_S68 * _S60 + _S68 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0[int(0)] + make_float2 (_S63 * u_0 + _S68 * dist_coeffs_0[int(5)] + _S68 * dist_coeffs_0[int(6)], _S65 * u_0 + (_S68 + (_S66 + _S66)) * dist_coeffs_0[int(4)] + _S68 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S70 = transpose_0(makeMatrix<float, 2, 2> (_S67 + make_float2 (_S67.x * dist_coeffs_0[int(8)] + _S67.y * dist_coeffs_0[int(9)], 0.0f), _S69 + make_float2 (_S69.x * dist_coeffs_0[int(8)] + _S69.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S70)), ((F32_min((_S70.rows[int(0)].x), (_S70.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_1 = uv_0[int(0)].x;
                float v_1 = uv_0[int(0)].y;
                float r2_1 = u_1 * u_1 + v_1 * v_1;
                float2  _S71 = uv_0[int(0)] * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
                float2  _S72 = _S71 + make_float2 (dist_coeffs_0[int(8)] * _S71.x + dist_coeffs_0[int(9)] * _S71.y, 0.0f);
                uv_0[int(0)] = make_float2 (fx_0 * _S72.x + cx_0, fy_0 * _S72.y + cy_0);
                break;
            }
            bool all_valid_0 = true & (!_S56);
            float3  _S73 = pos_c_0[int(1)];
            _S16 = &uv_0[int(1)];
            for(;;)
            {
                float _S74 = _S73.z;
                uv_0[int(1)] = float2 {_S73.x, _S73.y} / make_float2 (_S74);
                if(_S74 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_2 = uv_0[int(1)].x;
                    float v_2 = uv_0[int(1)].y;
                    float _S75 = u_2 + u_2;
                    float r2_2 = u_2 * u_2 + v_2 * v_2;
                    float _S76 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
                    float _S77 = dist_coeffs_0[int(1)] + r2_2 * _S76;
                    float _S78 = dist_coeffs_0[int(0)] + r2_2 * _S77;
                    float radial_1 = 1.0f + r2_2 * _S78;
                    float _S79 = 2.0f * dist_coeffs_0[int(4)];
                    float _S80 = 2.0f * u_2;
                    float _S81 = 2.0f * dist_coeffs_0[int(5)];
                    float _S82 = 2.0f * v_2;
                    float2  _S83 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (_S75 * _S78 + (_S75 * _S77 + (_S75 * _S76 + _S75 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv_0[int(1)] + make_float2 (_S79 * v_2 + (_S75 + (_S80 + _S80)) * dist_coeffs_0[int(5)] + _S75 * dist_coeffs_0[int(6)], _S81 * v_2 + _S75 * dist_coeffs_0[int(4)] + _S75 * dist_coeffs_0[int(7)]);
                    float _S84 = v_2 + v_2;
                    float2  _S85 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (_S84 * _S78 + (_S84 * _S77 + (_S84 * _S76 + _S84 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv_0[int(1)] + make_float2 (_S79 * u_2 + _S84 * dist_coeffs_0[int(5)] + _S84 * dist_coeffs_0[int(6)], _S81 * u_2 + (_S84 + (_S82 + _S82)) * dist_coeffs_0[int(4)] + _S84 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S86 = transpose_0(makeMatrix<float, 2, 2> (_S83 + make_float2 (_S83.x * dist_coeffs_0[int(8)] + _S83.y * dist_coeffs_0[int(9)], 0.0f), _S85 + make_float2 (_S85.x * dist_coeffs_0[int(8)] + _S85.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S86)), ((F32_min((_S86.rows[int(0)].x), (_S86.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_3 = uv_0[int(1)].x;
                float v_3 = uv_0[int(1)].y;
                float r2_3 = u_3 * u_3 + v_3 * v_3;
                float2  _S87 = uv_0[int(1)] * make_float2 (1.0f + r2_3 * (dist_coeffs_0[int(0)] + r2_3 * (dist_coeffs_0[int(1)] + r2_3 * (dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_3 * v_3 + dist_coeffs_0[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + dist_coeffs_0[int(6)] * r2_3, 2.0f * dist_coeffs_0[int(5)] * u_3 * v_3 + dist_coeffs_0[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + dist_coeffs_0[int(7)] * r2_3);
                float2  _S88 = _S87 + make_float2 (dist_coeffs_0[int(8)] * _S87.x + dist_coeffs_0[int(9)] * _S87.y, 0.0f);
                uv_0[int(1)] = make_float2 (fx_0 * _S88.x + cx_0, fy_0 * _S88.y + cy_0);
                break;
            }
            bool all_valid_1 = all_valid_0 & (!_S56);
            float3  _S89 = pos_c_0[int(2)];
            _S17 = &uv_0[int(2)];
            for(;;)
            {
                float _S90 = _S89.z;
                uv_0[int(2)] = float2 {_S89.x, _S89.y} / make_float2 (_S90);
                if(_S90 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_4 = uv_0[int(2)].x;
                    float v_4 = uv_0[int(2)].y;
                    float _S91 = u_4 + u_4;
                    float r2_4 = u_4 * u_4 + v_4 * v_4;
                    float _S92 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
                    float _S93 = dist_coeffs_0[int(1)] + r2_4 * _S92;
                    float _S94 = dist_coeffs_0[int(0)] + r2_4 * _S93;
                    float radial_2 = 1.0f + r2_4 * _S94;
                    float _S95 = 2.0f * dist_coeffs_0[int(4)];
                    float _S96 = 2.0f * u_4;
                    float _S97 = 2.0f * dist_coeffs_0[int(5)];
                    float _S98 = 2.0f * v_4;
                    float2  _S99 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (_S91 * _S94 + (_S91 * _S93 + (_S91 * _S92 + _S91 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv_0[int(2)] + make_float2 (_S95 * v_4 + (_S91 + (_S96 + _S96)) * dist_coeffs_0[int(5)] + _S91 * dist_coeffs_0[int(6)], _S97 * v_4 + _S91 * dist_coeffs_0[int(4)] + _S91 * dist_coeffs_0[int(7)]);
                    float _S100 = v_4 + v_4;
                    float2  _S101 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (_S100 * _S94 + (_S100 * _S93 + (_S100 * _S92 + _S100 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv_0[int(2)] + make_float2 (_S95 * u_4 + _S100 * dist_coeffs_0[int(5)] + _S100 * dist_coeffs_0[int(6)], _S97 * u_4 + (_S100 + (_S98 + _S98)) * dist_coeffs_0[int(4)] + _S100 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S102 = transpose_0(makeMatrix<float, 2, 2> (_S99 + make_float2 (_S99.x * dist_coeffs_0[int(8)] + _S99.y * dist_coeffs_0[int(9)], 0.0f), _S101 + make_float2 (_S101.x * dist_coeffs_0[int(8)] + _S101.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S102)), ((F32_min((_S102.rows[int(0)].x), (_S102.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_5 = uv_0[int(2)].x;
                float v_5 = uv_0[int(2)].y;
                float r2_5 = u_5 * u_5 + v_5 * v_5;
                float2  _S103 = uv_0[int(2)] * make_float2 (1.0f + r2_5 * (dist_coeffs_0[int(0)] + r2_5 * (dist_coeffs_0[int(1)] + r2_5 * (dist_coeffs_0[int(2)] + r2_5 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_5 * v_5 + dist_coeffs_0[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + dist_coeffs_0[int(6)] * r2_5, 2.0f * dist_coeffs_0[int(5)] * u_5 * v_5 + dist_coeffs_0[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + dist_coeffs_0[int(7)] * r2_5);
                float2  _S104 = _S103 + make_float2 (dist_coeffs_0[int(8)] * _S103.x + dist_coeffs_0[int(9)] * _S103.y, 0.0f);
                uv_0[int(2)] = make_float2 (fx_0 * _S104.x + cx_0, fy_0 * _S104.y + cy_0);
                break;
            }
            bool all_valid_2 = all_valid_1 & (!_S56);
            float3  _S105 = pos_c_0[int(3)];
            _S18 = &uv_0[int(3)];
            for(;;)
            {
                float _S106 = _S105.z;
                uv_0[int(3)] = float2 {_S105.x, _S105.y} / make_float2 (_S106);
                if(_S106 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_6 = uv_0[int(3)].x;
                    float v_6 = uv_0[int(3)].y;
                    float _S107 = u_6 + u_6;
                    float r2_6 = u_6 * u_6 + v_6 * v_6;
                    float _S108 = dist_coeffs_0[int(2)] + r2_6 * dist_coeffs_0[int(3)];
                    float _S109 = dist_coeffs_0[int(1)] + r2_6 * _S108;
                    float _S110 = dist_coeffs_0[int(0)] + r2_6 * _S109;
                    float radial_3 = 1.0f + r2_6 * _S110;
                    float _S111 = 2.0f * dist_coeffs_0[int(4)];
                    float _S112 = 2.0f * u_6;
                    float _S113 = 2.0f * dist_coeffs_0[int(5)];
                    float _S114 = 2.0f * v_6;
                    float2  _S115 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (_S107 * _S110 + (_S107 * _S109 + (_S107 * _S108 + _S107 * dist_coeffs_0[int(3)] * r2_6) * r2_6) * r2_6) * uv_0[int(3)] + make_float2 (_S111 * v_6 + (_S107 + (_S112 + _S112)) * dist_coeffs_0[int(5)] + _S107 * dist_coeffs_0[int(6)], _S113 * v_6 + _S107 * dist_coeffs_0[int(4)] + _S107 * dist_coeffs_0[int(7)]);
                    float _S116 = v_6 + v_6;
                    float2  _S117 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (_S116 * _S110 + (_S116 * _S109 + (_S116 * _S108 + _S116 * dist_coeffs_0[int(3)] * r2_6) * r2_6) * r2_6) * uv_0[int(3)] + make_float2 (_S111 * u_6 + _S116 * dist_coeffs_0[int(5)] + _S116 * dist_coeffs_0[int(6)], _S113 * u_6 + (_S116 + (_S114 + _S114)) * dist_coeffs_0[int(4)] + _S116 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S118 = transpose_0(makeMatrix<float, 2, 2> (_S115 + make_float2 (_S115.x * dist_coeffs_0[int(8)] + _S115.y * dist_coeffs_0[int(9)], 0.0f), _S117 + make_float2 (_S117.x * dist_coeffs_0[int(8)] + _S117.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S118)), ((F32_min((_S118.rows[int(0)].x), (_S118.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_7 = uv_0[int(3)].x;
                float v_7 = uv_0[int(3)].y;
                float r2_7 = u_7 * u_7 + v_7 * v_7;
                float2  _S119 = uv_0[int(3)] * make_float2 (1.0f + r2_7 * (dist_coeffs_0[int(0)] + r2_7 * (dist_coeffs_0[int(1)] + r2_7 * (dist_coeffs_0[int(2)] + r2_7 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_7 * v_7 + dist_coeffs_0[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + dist_coeffs_0[int(6)] * r2_7, 2.0f * dist_coeffs_0[int(5)] * u_7 * v_7 + dist_coeffs_0[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + dist_coeffs_0[int(7)] * r2_7);
                float2  _S120 = _S119 + make_float2 (dist_coeffs_0[int(8)] * _S119.x + dist_coeffs_0[int(9)] * _S119.y, 0.0f);
                uv_0[int(3)] = make_float2 (fx_0 * _S120.x + cx_0, fy_0 * _S120.y + cy_0);
                break;
            }
            bool all_valid_3 = all_valid_2 & (!_S56);
            float3  _S121 = pos_c_0[int(4)];
            _S19 = &uv_0[int(4)];
            for(;;)
            {
                float _S122 = _S121.z;
                uv_0[int(4)] = float2 {_S121.x, _S121.y} / make_float2 (_S122);
                if(_S122 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_8 = uv_0[int(4)].x;
                    float v_8 = uv_0[int(4)].y;
                    float _S123 = u_8 + u_8;
                    float r2_8 = u_8 * u_8 + v_8 * v_8;
                    float _S124 = dist_coeffs_0[int(2)] + r2_8 * dist_coeffs_0[int(3)];
                    float _S125 = dist_coeffs_0[int(1)] + r2_8 * _S124;
                    float _S126 = dist_coeffs_0[int(0)] + r2_8 * _S125;
                    float radial_4 = 1.0f + r2_8 * _S126;
                    float _S127 = 2.0f * dist_coeffs_0[int(4)];
                    float _S128 = 2.0f * u_8;
                    float _S129 = 2.0f * dist_coeffs_0[int(5)];
                    float _S130 = 2.0f * v_8;
                    float2  _S131 = make_float2 (1.0f, 0.0f) * make_float2 (radial_4) + make_float2 (_S123 * _S126 + (_S123 * _S125 + (_S123 * _S124 + _S123 * dist_coeffs_0[int(3)] * r2_8) * r2_8) * r2_8) * uv_0[int(4)] + make_float2 (_S127 * v_8 + (_S123 + (_S128 + _S128)) * dist_coeffs_0[int(5)] + _S123 * dist_coeffs_0[int(6)], _S129 * v_8 + _S123 * dist_coeffs_0[int(4)] + _S123 * dist_coeffs_0[int(7)]);
                    float _S132 = v_8 + v_8;
                    float2  _S133 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (_S132 * _S126 + (_S132 * _S125 + (_S132 * _S124 + _S132 * dist_coeffs_0[int(3)] * r2_8) * r2_8) * r2_8) * uv_0[int(4)] + make_float2 (_S127 * u_8 + _S132 * dist_coeffs_0[int(5)] + _S132 * dist_coeffs_0[int(6)], _S129 * u_8 + (_S132 + (_S130 + _S130)) * dist_coeffs_0[int(4)] + _S132 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S134 = transpose_0(makeMatrix<float, 2, 2> (_S131 + make_float2 (_S131.x * dist_coeffs_0[int(8)] + _S131.y * dist_coeffs_0[int(9)], 0.0f), _S133 + make_float2 (_S133.x * dist_coeffs_0[int(8)] + _S133.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S134)), ((F32_min((_S134.rows[int(0)].x), (_S134.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_9 = uv_0[int(4)].x;
                float v_9 = uv_0[int(4)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S135 = uv_0[int(4)] * make_float2 (1.0f + r2_9 * (dist_coeffs_0[int(0)] + r2_9 * (dist_coeffs_0[int(1)] + r2_9 * (dist_coeffs_0[int(2)] + r2_9 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_9 * v_9 + dist_coeffs_0[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_0[int(6)] * r2_9, 2.0f * dist_coeffs_0[int(5)] * u_9 * v_9 + dist_coeffs_0[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_0[int(7)] * r2_9);
                float2  _S136 = _S135 + make_float2 (dist_coeffs_0[int(8)] * _S135.x + dist_coeffs_0[int(9)] * _S135.y, 0.0f);
                uv_0[int(4)] = make_float2 (fx_0 * _S136.x + cx_0, fy_0 * _S136.y + cy_0);
                break;
            }
            bool all_valid_4 = all_valid_3 & (!_S56);
            float3  _S137 = pos_c_0[int(5)];
            _S20 = &uv_0[int(5)];
            for(;;)
            {
                float _S138 = _S137.z;
                uv_0[int(5)] = float2 {_S137.x, _S137.y} / make_float2 (_S138);
                if(_S138 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_10 = uv_0[int(5)].x;
                    float v_10 = uv_0[int(5)].y;
                    float _S139 = u_10 + u_10;
                    float r2_10 = u_10 * u_10 + v_10 * v_10;
                    float _S140 = dist_coeffs_0[int(2)] + r2_10 * dist_coeffs_0[int(3)];
                    float _S141 = dist_coeffs_0[int(1)] + r2_10 * _S140;
                    float _S142 = dist_coeffs_0[int(0)] + r2_10 * _S141;
                    float radial_5 = 1.0f + r2_10 * _S142;
                    float _S143 = 2.0f * dist_coeffs_0[int(4)];
                    float _S144 = 2.0f * u_10;
                    float _S145 = 2.0f * dist_coeffs_0[int(5)];
                    float _S146 = 2.0f * v_10;
                    float2  _S147 = make_float2 (1.0f, 0.0f) * make_float2 (radial_5) + make_float2 (_S139 * _S142 + (_S139 * _S141 + (_S139 * _S140 + _S139 * dist_coeffs_0[int(3)] * r2_10) * r2_10) * r2_10) * uv_0[int(5)] + make_float2 (_S143 * v_10 + (_S139 + (_S144 + _S144)) * dist_coeffs_0[int(5)] + _S139 * dist_coeffs_0[int(6)], _S145 * v_10 + _S139 * dist_coeffs_0[int(4)] + _S139 * dist_coeffs_0[int(7)]);
                    float _S148 = v_10 + v_10;
                    float2  _S149 = make_float2 (0.0f, 1.0f) * make_float2 (radial_5) + make_float2 (_S148 * _S142 + (_S148 * _S141 + (_S148 * _S140 + _S148 * dist_coeffs_0[int(3)] * r2_10) * r2_10) * r2_10) * uv_0[int(5)] + make_float2 (_S143 * u_10 + _S148 * dist_coeffs_0[int(5)] + _S148 * dist_coeffs_0[int(6)], _S145 * u_10 + (_S148 + (_S146 + _S146)) * dist_coeffs_0[int(4)] + _S148 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S150 = transpose_0(makeMatrix<float, 2, 2> (_S147 + make_float2 (_S147.x * dist_coeffs_0[int(8)] + _S147.y * dist_coeffs_0[int(9)], 0.0f), _S149 + make_float2 (_S149.x * dist_coeffs_0[int(8)] + _S149.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S150)), ((F32_min((_S150.rows[int(0)].x), (_S150.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_11 = uv_0[int(5)].x;
                float v_11 = uv_0[int(5)].y;
                float r2_11 = u_11 * u_11 + v_11 * v_11;
                float2  _S151 = uv_0[int(5)] * make_float2 (1.0f + r2_11 * (dist_coeffs_0[int(0)] + r2_11 * (dist_coeffs_0[int(1)] + r2_11 * (dist_coeffs_0[int(2)] + r2_11 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_11 * v_11 + dist_coeffs_0[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_0[int(6)] * r2_11, 2.0f * dist_coeffs_0[int(5)] * u_11 * v_11 + dist_coeffs_0[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_0[int(7)] * r2_11);
                float2  _S152 = _S151 + make_float2 (dist_coeffs_0[int(8)] * _S151.x + dist_coeffs_0[int(9)] * _S151.y, 0.0f);
                uv_0[int(5)] = make_float2 (fx_0 * _S152.x + cx_0, fy_0 * _S152.y + cy_0);
                break;
            }
            bool all_valid_5 = all_valid_4 & (!_S56);
            float3  _S153 = pos_c_0[int(6)];
            _S21 = &uv_0[int(6)];
            for(;;)
            {
                float _S154 = _S153.z;
                uv_0[int(6)] = float2 {_S153.x, _S153.y} / make_float2 (_S154);
                if(_S154 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_12 = uv_0[int(6)].x;
                    float v_12 = uv_0[int(6)].y;
                    float _S155 = u_12 + u_12;
                    float r2_12 = u_12 * u_12 + v_12 * v_12;
                    float _S156 = dist_coeffs_0[int(2)] + r2_12 * dist_coeffs_0[int(3)];
                    float _S157 = dist_coeffs_0[int(1)] + r2_12 * _S156;
                    float _S158 = dist_coeffs_0[int(0)] + r2_12 * _S157;
                    float radial_6 = 1.0f + r2_12 * _S158;
                    float _S159 = 2.0f * dist_coeffs_0[int(4)];
                    float _S160 = 2.0f * u_12;
                    float _S161 = 2.0f * dist_coeffs_0[int(5)];
                    float _S162 = 2.0f * v_12;
                    float2  _S163 = make_float2 (1.0f, 0.0f) * make_float2 (radial_6) + make_float2 (_S155 * _S158 + (_S155 * _S157 + (_S155 * _S156 + _S155 * dist_coeffs_0[int(3)] * r2_12) * r2_12) * r2_12) * uv_0[int(6)] + make_float2 (_S159 * v_12 + (_S155 + (_S160 + _S160)) * dist_coeffs_0[int(5)] + _S155 * dist_coeffs_0[int(6)], _S161 * v_12 + _S155 * dist_coeffs_0[int(4)] + _S155 * dist_coeffs_0[int(7)]);
                    float _S164 = v_12 + v_12;
                    float2  _S165 = make_float2 (0.0f, 1.0f) * make_float2 (radial_6) + make_float2 (_S164 * _S158 + (_S164 * _S157 + (_S164 * _S156 + _S164 * dist_coeffs_0[int(3)] * r2_12) * r2_12) * r2_12) * uv_0[int(6)] + make_float2 (_S159 * u_12 + _S164 * dist_coeffs_0[int(5)] + _S164 * dist_coeffs_0[int(6)], _S161 * u_12 + (_S164 + (_S162 + _S162)) * dist_coeffs_0[int(4)] + _S164 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S166 = transpose_0(makeMatrix<float, 2, 2> (_S163 + make_float2 (_S163.x * dist_coeffs_0[int(8)] + _S163.y * dist_coeffs_0[int(9)], 0.0f), _S165 + make_float2 (_S165.x * dist_coeffs_0[int(8)] + _S165.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S166)), ((F32_min((_S166.rows[int(0)].x), (_S166.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_13 = uv_0[int(6)].x;
                float v_13 = uv_0[int(6)].y;
                float r2_13 = u_13 * u_13 + v_13 * v_13;
                float2  _S167 = uv_0[int(6)] * make_float2 (1.0f + r2_13 * (dist_coeffs_0[int(0)] + r2_13 * (dist_coeffs_0[int(1)] + r2_13 * (dist_coeffs_0[int(2)] + r2_13 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_13 * v_13 + dist_coeffs_0[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_0[int(6)] * r2_13, 2.0f * dist_coeffs_0[int(5)] * u_13 * v_13 + dist_coeffs_0[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_0[int(7)] * r2_13);
                float2  _S168 = _S167 + make_float2 (dist_coeffs_0[int(8)] * _S167.x + dist_coeffs_0[int(9)] * _S167.y, 0.0f);
                uv_0[int(6)] = make_float2 (fx_0 * _S168.x + cx_0, fy_0 * _S168.y + cy_0);
                break;
            }
            bool all_valid_6 = all_valid_5 & (!_S56);
            float3  _S169 = pos_c_0[int(7)];
            _S22 = &uv_0[int(7)];
            for(;;)
            {
                float _S170 = _S169.z;
                uv_0[int(7)] = float2 {_S169.x, _S169.y} / make_float2 (_S170);
                if(_S170 < 0.0f)
                {
                    _S56 = true;
                }
                else
                {
                    float u_14 = uv_0[int(7)].x;
                    float v_14 = uv_0[int(7)].y;
                    float _S171 = u_14 + u_14;
                    float r2_14 = u_14 * u_14 + v_14 * v_14;
                    float _S172 = dist_coeffs_0[int(2)] + r2_14 * dist_coeffs_0[int(3)];
                    float _S173 = dist_coeffs_0[int(1)] + r2_14 * _S172;
                    float _S174 = dist_coeffs_0[int(0)] + r2_14 * _S173;
                    float radial_7 = 1.0f + r2_14 * _S174;
                    float _S175 = 2.0f * dist_coeffs_0[int(4)];
                    float _S176 = 2.0f * u_14;
                    float _S177 = 2.0f * dist_coeffs_0[int(5)];
                    float _S178 = 2.0f * v_14;
                    float2  _S179 = make_float2 (1.0f, 0.0f) * make_float2 (radial_7) + make_float2 (_S171 * _S174 + (_S171 * _S173 + (_S171 * _S172 + _S171 * dist_coeffs_0[int(3)] * r2_14) * r2_14) * r2_14) * uv_0[int(7)] + make_float2 (_S175 * v_14 + (_S171 + (_S176 + _S176)) * dist_coeffs_0[int(5)] + _S171 * dist_coeffs_0[int(6)], _S177 * v_14 + _S171 * dist_coeffs_0[int(4)] + _S171 * dist_coeffs_0[int(7)]);
                    float _S180 = v_14 + v_14;
                    float2  _S181 = make_float2 (0.0f, 1.0f) * make_float2 (radial_7) + make_float2 (_S180 * _S174 + (_S180 * _S173 + (_S180 * _S172 + _S180 * dist_coeffs_0[int(3)] * r2_14) * r2_14) * r2_14) * uv_0[int(7)] + make_float2 (_S175 * u_14 + _S180 * dist_coeffs_0[int(5)] + _S180 * dist_coeffs_0[int(6)], _S177 * u_14 + (_S180 + (_S178 + _S178)) * dist_coeffs_0[int(4)] + _S180 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S182 = transpose_0(makeMatrix<float, 2, 2> (_S179 + make_float2 (_S179.x * dist_coeffs_0[int(8)] + _S179.y * dist_coeffs_0[int(9)], 0.0f), _S181 + make_float2 (_S181.x * dist_coeffs_0[int(8)] + _S181.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S56 = !((F32_min((determinant_0(_S182)), ((F32_min((_S182.rows[int(0)].x), (_S182.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S56)
                {
                    break;
                }
                float u_15 = uv_0[int(7)].x;
                float v_15 = uv_0[int(7)].y;
                float r2_15 = u_15 * u_15 + v_15 * v_15;
                float2  _S183 = uv_0[int(7)] * make_float2 (1.0f + r2_15 * (dist_coeffs_0[int(0)] + r2_15 * (dist_coeffs_0[int(1)] + r2_15 * (dist_coeffs_0[int(2)] + r2_15 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_15 * v_15 + dist_coeffs_0[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + dist_coeffs_0[int(6)] * r2_15, 2.0f * dist_coeffs_0[int(5)] * u_15 * v_15 + dist_coeffs_0[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + dist_coeffs_0[int(7)] * r2_15);
                float2  _S184 = _S183 + make_float2 (dist_coeffs_0[int(8)] * _S183.x + dist_coeffs_0[int(9)] * _S183.y, 0.0f);
                uv_0[int(7)] = make_float2 (fx_0 * _S184.x + cx_0, fy_0 * _S184.y + cy_0);
                break;
            }
            _S23 = all_valid_6 & (!_S56);
            break;
        }
        if(!_S23)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((densities_0[int(0)]), (densities_0[int(1)])))), (densities_0[int(2)])))), (densities_0[int(3)])))), (densities_0[int(4)])))), (densities_0[int(5)])))), (densities_0[int(6)])))), (densities_0[int(7)]))) * size_0 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float _S185 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S15).x), ((*_S16).x)))), ((*_S17).x)))), ((*_S18).x)))), ((*_S19).x)))), ((*_S20).x)))), ((*_S21).x)))), ((*_S22).x)));
        float _S186 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S15).x), ((*_S16).x)))), ((*_S17).x)))), ((*_S18).x)))), ((*_S19).x)))), ((*_S20).x)))), ((*_S21).x)))), ((*_S22).x)));
        float _S187 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S15).y), ((*_S16).y)))), ((*_S17).y)))), ((*_S18).y)))), ((*_S19).y)))), ((*_S20).y)))), ((*_S21).y)))), ((*_S22).y)));
        float _S188 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S15).y), ((*_S16).y)))), ((*_S17).y)))), ((*_S18).y)))), ((*_S19).y)))), ((*_S20).y)))), ((*_S21).y)))), ((*_S22).y)));
        if(_S185 <= 0.0f)
        {
            _S56 = true;
        }
        else
        {
            _S56 = _S186 >= float(image_width_0);
        }
        if(_S56)
        {
            _S56 = true;
        }
        else
        {
            _S56 = _S187 <= 0.0f;
        }
        if(_S56)
        {
            _S56 = true;
        }
        else
        {
            _S56 = _S188 >= float(image_height_0);
        }
        if(_S56)
        {
            _S56 = true;
        }
        else
        {
            if(_S54 <= 0.0f)
            {
                if(_S186 <= 0.0f)
                {
                    _S56 = _S185 >= float(image_width_0);
                }
                else
                {
                    _S56 = false;
                }
                if(_S56)
                {
                    _S56 = true;
                }
                else
                {
                    if(_S188 <= 0.0f)
                    {
                        _S56 = _S187 >= float(image_width_0);
                    }
                    else
                    {
                        _S56 = false;
                    }
                }
            }
            else
            {
                _S56 = false;
            }
        }
        if(_S56)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (float(int((F32_floor((_S186))))), float(int((F32_floor((_S188))))), float(int((F32_ceil((_S185))))), float(int((F32_ceil((_S187))))));
        *depth_0 = mean_c_0.z;
        float3  _S189 = mean_c_0 - - mul_0(transpose_1(R_0), t_0);
        float _S190 = _S189.x;
        float _S191 = _S189.y;
        float _S192 = _S189.z;
        float norm_0 = (F32_sqrt((_S190 * _S190 + _S191 * _S191 + _S192 * _S192)));
        float x_7 = _S190 / norm_0;
        float y_3 = _S191 / norm_0;
        float z_0 = _S192 / norm_0;
        float z2_0 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_7 * x_7 - y_3 * y_3;
        float fS1_0 = 2.0f * x_7 * y_3;
        float fTmp0C_0 = -2.28522896766662598f * z2_0 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgbs_0 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_0[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * sh_coeffs_0[int(1)] + make_float3 (z_0) * sh_coeffs_0[int(2)] - make_float3 (x_7) * sh_coeffs_0[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * sh_coeffs_0[int(4)] + make_float3 (fTmp0B_0 * y_3) * sh_coeffs_0[int(5)] + make_float3 (0.94617468118667603f * z2_0 - 0.31539157032966614f) * sh_coeffs_0[int(6)] + make_float3 (fTmp0B_0 * x_7) * sh_coeffs_0[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * sh_coeffs_0[int(8)]) + (make_float3 (-0.59004360437393188f * (x_7 * fS1_0 + y_3 * fC1_0)) * sh_coeffs_0[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * sh_coeffs_0[int(10)] + make_float3 (fTmp0C_0 * y_3) * sh_coeffs_0[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_0 - 1.11952900886535645f)) * sh_coeffs_0[int(12)] + make_float3 (fTmp0C_0 * x_7) * sh_coeffs_0[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * sh_coeffs_0[int(14)] + make_float3 (-0.59004360437393188f * (x_7 * fC1_0 - y_3 * fS1_0)) * sh_coeffs_0[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  densities_1, FixedArray<float3 , 16>  sh_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, float4  * aabb_xyxy_1, float * depth_1, float3  * rgbs_1)
{
    float2  * _S193;
    float2  _S194;
    float _S195;
    float _S196;
    float _S197;
    float _S198;
    float _S199;
    float _S200;
    float _S201;
    float _S202;
    float _S203;
    float _S204;
    float _S205;
    float _S206;
    float2  _S207;
    bool _S208;
    float2  * _S209;
    bool _S210;
    float2  * _S211;
    bool _S212;
    float2  * _S213;
    bool _S214;
    float2  * _S215;
    bool _S216;
    float2  * _S217;
    bool _S218;
    float2  * _S219;
    bool _S220;
    float2  * _S221;
    bool _S222;
    bool _S223;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S224 = mul_0(R_1, pos_1) + t_1;
        pos_c_1[int(0)] = _S224;
        float _S225 = length_1(_S224);
        float _S226 = (F32_min((far_plane_1), (_S225)));
        float _S227 = (F32_max((near_plane_1), (_S225)));
        float3  _S228 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_1;
        pos_c_1[int(1)] = _S228;
        float _S229 = length_1(_S228);
        float _S230 = (F32_min((_S226), (_S229)));
        float _S231 = (F32_max((_S227), (_S229)));
        float3  _S232 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(2)] = _S232;
        float _S233 = length_1(_S232);
        float _S234 = (F32_min((_S230), (_S233)));
        float _S235 = (F32_max((_S231), (_S233)));
        float3  _S236 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(3)] = _S236;
        float _S237 = length_1(_S236);
        float _S238 = (F32_min((_S234), (_S237)));
        float _S239 = (F32_max((_S235), (_S237)));
        float3  _S240 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(4)] = _S240;
        float _S241 = length_1(_S240);
        float _S242 = (F32_min((_S238), (_S241)));
        float _S243 = (F32_max((_S239), (_S241)));
        float3  _S244 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(5)] = _S244;
        float _S245 = length_1(_S244);
        float _S246 = (F32_min((_S242), (_S245)));
        float _S247 = (F32_max((_S243), (_S245)));
        float3  _S248 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_1;
        pos_c_1[int(6)] = _S248;
        float _S249 = length_1(_S248);
        float _S250 = (F32_min((_S246), (_S249)));
        float _S251 = (F32_max((_S247), (_S249)));
        float3  _S252 = mul_0(R_1, pos_1 + make_float3 (size_1)) + t_1;
        pos_c_1[int(7)] = _S252;
        float _S253 = length_1(_S252);
        float _S254 = (F32_min((_S250), (_S253)));
        float _S255 = (F32_max((_S251), (_S253)));
        bool _S256;
        if(_S254 < near_plane_1)
        {
            _S256 = true;
        }
        else
        {
            _S256 = _S255 > far_plane_1;
        }
        if(_S256)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float3  mean_c_1 = mul_0(R_1, pos_1 + make_float3 (0.5f * size_1)) + t_1;
        FixedArray<float2 , 8>  uv_1;
        for(;;)
        {
            float k_0;
            float3  _S257 = pos_c_1[int(0)];
            _S193 = &uv_1[int(0)];
            for(;;)
            {
                float2  _S258 = float2 {_S257.x, _S257.y};
                float r_2 = length_0(_S258);
                float _S259 = _S257.z;
                float theta_0 = (F32_atan2((r_2), (_S259)));
                if(theta_0 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S259;
                }
                else
                {
                    k_0 = theta_0 / r_2;
                }
                float2  _S260 = _S258 * make_float2 (k_0);
                uv_1[int(0)] = _S260;
                float2  _S261 = make_float2 (1.0f, 0.0f);
                _S194 = _S261;
                _S195 = dist_coeffs_1[int(0)];
                _S196 = dist_coeffs_1[int(1)];
                _S197 = dist_coeffs_1[int(2)];
                _S198 = dist_coeffs_1[int(3)];
                _S199 = dist_coeffs_1[int(4)];
                _S200 = dist_coeffs_1[int(5)];
                _S201 = dist_coeffs_1[int(6)];
                _S202 = dist_coeffs_1[int(7)];
                _S203 = dist_coeffs_1[int(8)];
                _S204 = dist_coeffs_1[int(9)];
                float u_16 = _S260.x;
                float v_16 = _S260.y;
                float _S262 = u_16 + u_16;
                float r2_16 = u_16 * u_16 + v_16 * v_16;
                float _S263 = dist_coeffs_1[int(2)] + r2_16 * dist_coeffs_1[int(3)];
                float _S264 = dist_coeffs_1[int(1)] + r2_16 * _S263;
                float _S265 = dist_coeffs_1[int(0)] + r2_16 * _S264;
                float _S266 = _S262 * _S265 + (_S262 * _S264 + (_S262 * _S263 + _S262 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16;
                float radial_8 = 1.0f + r2_16 * _S265;
                float _S267 = 2.0f * dist_coeffs_1[int(4)];
                _S205 = _S267;
                float _S268 = _S267 * u_16;
                float _S269 = 2.0f * u_16;
                float s_diff_du_0 = _S267 * v_16 + (_S262 + (_S269 + _S269)) * dist_coeffs_1[int(5)] + _S262 * dist_coeffs_1[int(6)];
                float _S270 = 2.0f * dist_coeffs_1[int(5)];
                _S206 = _S270;
                float _S271 = _S270 * u_16;
                float _S272 = 2.0f * v_16;
                float2  _S273 = _S261 * make_float2 (radial_8) + make_float2 (_S266) * _S260 + make_float2 (s_diff_du_0, _S270 * v_16 + _S262 * dist_coeffs_1[int(4)] + _S262 * dist_coeffs_1[int(7)]);
                float2  _S274 = _S273 + make_float2 (_S273.x * dist_coeffs_1[int(8)] + _S273.y * dist_coeffs_1[int(9)], 0.0f);
                float2  _S275 = make_float2 (0.0f, 1.0f);
                _S207 = _S275;
                float _S276 = v_16 + v_16;
                float2  _S277 = _S275 * make_float2 (radial_8) + make_float2 (_S276 * _S265 + (_S276 * _S264 + (_S276 * _S263 + _S276 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16) * _S260 + make_float2 (_S268 + _S276 * dist_coeffs_1[int(5)] + _S276 * dist_coeffs_1[int(6)], _S271 + (_S276 + (_S272 + _S272)) * dist_coeffs_1[int(4)] + _S276 * dist_coeffs_1[int(7)]);
                Matrix<float, 2, 2>  _S278 = transpose_0(makeMatrix<float, 2, 2> (_S274, _S277 + make_float2 (_S277.x * dist_coeffs_1[int(8)] + _S277.y * dist_coeffs_1[int(9)], 0.0f)));
                bool _S279 = !((F32_min((determinant_0(_S278)), ((F32_min((_S278.rows[int(0)].x), (_S278.rows[int(1)].y)))))) > 0.0f);
                _S208 = _S279;
                if(_S279)
                {
                    break;
                }
                float u_17 = uv_1[int(0)].x;
                float v_17 = uv_1[int(0)].y;
                float r2_17 = u_17 * u_17 + v_17 * v_17;
                float2  _S280 = uv_1[int(0)] * make_float2 (1.0f + r2_17 * (dist_coeffs_1[int(0)] + r2_17 * (dist_coeffs_1[int(1)] + r2_17 * (dist_coeffs_1[int(2)] + r2_17 * dist_coeffs_1[int(3)])))) + make_float2 (_S267 * u_17 * v_17 + dist_coeffs_1[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + dist_coeffs_1[int(6)] * r2_17, _S270 * u_17 * v_17 + dist_coeffs_1[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + dist_coeffs_1[int(7)] * r2_17);
                float2  _S281 = _S280 + make_float2 (dist_coeffs_1[int(8)] * _S280.x + dist_coeffs_1[int(9)] * _S280.y, 0.0f);
                uv_1[int(0)] = make_float2 (fx_1 * _S281.x + cx_1, fy_1 * _S281.y + cy_1);
                break;
            }
            bool all_valid_7 = true & (!_S208);
            float3  _S282 = pos_c_1[int(1)];
            _S209 = &uv_1[int(1)];
            for(;;)
            {
                float2  _S283 = float2 {_S282.x, _S282.y};
                float r_3 = length_0(_S283);
                float _S284 = _S282.z;
                float theta_1 = (F32_atan2((r_3), (_S284)));
                if(theta_1 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S284;
                }
                else
                {
                    k_0 = theta_1 / r_3;
                }
                float2  _S285 = _S283 * make_float2 (k_0);
                uv_1[int(1)] = _S285;
                float u_18 = _S285.x;
                float v_18 = _S285.y;
                float _S286 = u_18 + u_18;
                float r2_18 = u_18 * u_18 + v_18 * v_18;
                float _S287 = _S197 + r2_18 * _S198;
                float _S288 = _S196 + r2_18 * _S287;
                float _S289 = _S195 + r2_18 * _S288;
                float radial_9 = 1.0f + r2_18 * _S289;
                float _S290 = 2.0f * u_18;
                float _S291 = 2.0f * v_18;
                float2  _S292 = _S194 * make_float2 (radial_9) + make_float2 (_S286 * _S289 + (_S286 * _S288 + (_S286 * _S287 + _S286 * _S198 * r2_18) * r2_18) * r2_18) * _S285 + make_float2 (_S205 * v_18 + (_S286 + (_S290 + _S290)) * _S200 + _S286 * _S201, _S206 * v_18 + _S286 * _S199 + _S286 * _S202);
                float _S293 = v_18 + v_18;
                float2  _S294 = _S207 * make_float2 (radial_9) + make_float2 (_S293 * _S289 + (_S293 * _S288 + (_S293 * _S287 + _S293 * _S198 * r2_18) * r2_18) * r2_18) * _S285 + make_float2 (_S205 * u_18 + _S293 * _S200 + _S293 * _S201, _S206 * u_18 + (_S293 + (_S291 + _S291)) * _S199 + _S293 * _S202);
                Matrix<float, 2, 2>  _S295 = transpose_0(makeMatrix<float, 2, 2> (_S292 + make_float2 (_S292.x * _S203 + _S292.y * _S204, 0.0f), _S294 + make_float2 (_S294.x * _S203 + _S294.y * _S204, 0.0f)));
                bool _S296 = !((F32_min((determinant_0(_S295)), ((F32_min((_S295.rows[int(0)].x), (_S295.rows[int(1)].y)))))) > 0.0f);
                _S210 = _S296;
                if(_S296)
                {
                    break;
                }
                float u_19 = uv_1[int(1)].x;
                float v_19 = uv_1[int(1)].y;
                float r2_19 = u_19 * u_19 + v_19 * v_19;
                float2  _S297 = uv_1[int(1)] * make_float2 (1.0f + r2_19 * (_S195 + r2_19 * (_S196 + r2_19 * (_S197 + r2_19 * _S198)))) + make_float2 (_S205 * u_19 * v_19 + _S200 * (r2_19 + 2.0f * u_19 * u_19) + _S201 * r2_19, _S206 * u_19 * v_19 + _S199 * (r2_19 + 2.0f * v_19 * v_19) + _S202 * r2_19);
                float2  _S298 = _S297 + make_float2 (_S203 * _S297.x + _S204 * _S297.y, 0.0f);
                uv_1[int(1)] = make_float2 (fx_1 * _S298.x + cx_1, fy_1 * _S298.y + cy_1);
                break;
            }
            bool all_valid_8 = all_valid_7 & (!_S210);
            float3  _S299 = pos_c_1[int(2)];
            _S211 = &uv_1[int(2)];
            for(;;)
            {
                float2  _S300 = float2 {_S299.x, _S299.y};
                float r_4 = length_0(_S300);
                float _S301 = _S299.z;
                float theta_2 = (F32_atan2((r_4), (_S301)));
                if(theta_2 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S301;
                }
                else
                {
                    k_0 = theta_2 / r_4;
                }
                float2  _S302 = _S300 * make_float2 (k_0);
                uv_1[int(2)] = _S302;
                float u_20 = _S302.x;
                float v_20 = _S302.y;
                float _S303 = u_20 + u_20;
                float r2_20 = u_20 * u_20 + v_20 * v_20;
                float _S304 = _S197 + r2_20 * _S198;
                float _S305 = _S196 + r2_20 * _S304;
                float _S306 = _S195 + r2_20 * _S305;
                float radial_10 = 1.0f + r2_20 * _S306;
                float _S307 = 2.0f * u_20;
                float _S308 = 2.0f * v_20;
                float2  _S309 = _S194 * make_float2 (radial_10) + make_float2 (_S303 * _S306 + (_S303 * _S305 + (_S303 * _S304 + _S303 * _S198 * r2_20) * r2_20) * r2_20) * _S302 + make_float2 (_S205 * v_20 + (_S303 + (_S307 + _S307)) * _S200 + _S303 * _S201, _S206 * v_20 + _S303 * _S199 + _S303 * _S202);
                float _S310 = v_20 + v_20;
                float2  _S311 = _S207 * make_float2 (radial_10) + make_float2 (_S310 * _S306 + (_S310 * _S305 + (_S310 * _S304 + _S310 * _S198 * r2_20) * r2_20) * r2_20) * _S302 + make_float2 (_S205 * u_20 + _S310 * _S200 + _S310 * _S201, _S206 * u_20 + (_S310 + (_S308 + _S308)) * _S199 + _S310 * _S202);
                Matrix<float, 2, 2>  _S312 = transpose_0(makeMatrix<float, 2, 2> (_S309 + make_float2 (_S309.x * _S203 + _S309.y * _S204, 0.0f), _S311 + make_float2 (_S311.x * _S203 + _S311.y * _S204, 0.0f)));
                bool _S313 = !((F32_min((determinant_0(_S312)), ((F32_min((_S312.rows[int(0)].x), (_S312.rows[int(1)].y)))))) > 0.0f);
                _S212 = _S313;
                if(_S313)
                {
                    break;
                }
                float u_21 = uv_1[int(2)].x;
                float v_21 = uv_1[int(2)].y;
                float r2_21 = u_21 * u_21 + v_21 * v_21;
                float2  _S314 = uv_1[int(2)] * make_float2 (1.0f + r2_21 * (_S195 + r2_21 * (_S196 + r2_21 * (_S197 + r2_21 * _S198)))) + make_float2 (_S205 * u_21 * v_21 + _S200 * (r2_21 + 2.0f * u_21 * u_21) + _S201 * r2_21, _S206 * u_21 * v_21 + _S199 * (r2_21 + 2.0f * v_21 * v_21) + _S202 * r2_21);
                float2  _S315 = _S314 + make_float2 (_S203 * _S314.x + _S204 * _S314.y, 0.0f);
                uv_1[int(2)] = make_float2 (fx_1 * _S315.x + cx_1, fy_1 * _S315.y + cy_1);
                break;
            }
            bool all_valid_9 = all_valid_8 & (!_S212);
            float3  _S316 = pos_c_1[int(3)];
            _S213 = &uv_1[int(3)];
            for(;;)
            {
                float2  _S317 = float2 {_S316.x, _S316.y};
                float r_5 = length_0(_S317);
                float _S318 = _S316.z;
                float theta_3 = (F32_atan2((r_5), (_S318)));
                if(theta_3 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S318;
                }
                else
                {
                    k_0 = theta_3 / r_5;
                }
                float2  _S319 = _S317 * make_float2 (k_0);
                uv_1[int(3)] = _S319;
                float u_22 = _S319.x;
                float v_22 = _S319.y;
                float _S320 = u_22 + u_22;
                float r2_22 = u_22 * u_22 + v_22 * v_22;
                float _S321 = _S197 + r2_22 * _S198;
                float _S322 = _S196 + r2_22 * _S321;
                float _S323 = _S195 + r2_22 * _S322;
                float radial_11 = 1.0f + r2_22 * _S323;
                float _S324 = 2.0f * u_22;
                float _S325 = 2.0f * v_22;
                float2  _S326 = _S194 * make_float2 (radial_11) + make_float2 (_S320 * _S323 + (_S320 * _S322 + (_S320 * _S321 + _S320 * _S198 * r2_22) * r2_22) * r2_22) * _S319 + make_float2 (_S205 * v_22 + (_S320 + (_S324 + _S324)) * _S200 + _S320 * _S201, _S206 * v_22 + _S320 * _S199 + _S320 * _S202);
                float _S327 = v_22 + v_22;
                float2  _S328 = _S207 * make_float2 (radial_11) + make_float2 (_S327 * _S323 + (_S327 * _S322 + (_S327 * _S321 + _S327 * _S198 * r2_22) * r2_22) * r2_22) * _S319 + make_float2 (_S205 * u_22 + _S327 * _S200 + _S327 * _S201, _S206 * u_22 + (_S327 + (_S325 + _S325)) * _S199 + _S327 * _S202);
                Matrix<float, 2, 2>  _S329 = transpose_0(makeMatrix<float, 2, 2> (_S326 + make_float2 (_S326.x * _S203 + _S326.y * _S204, 0.0f), _S328 + make_float2 (_S328.x * _S203 + _S328.y * _S204, 0.0f)));
                bool _S330 = !((F32_min((determinant_0(_S329)), ((F32_min((_S329.rows[int(0)].x), (_S329.rows[int(1)].y)))))) > 0.0f);
                _S214 = _S330;
                if(_S330)
                {
                    break;
                }
                float u_23 = uv_1[int(3)].x;
                float v_23 = uv_1[int(3)].y;
                float r2_23 = u_23 * u_23 + v_23 * v_23;
                float2  _S331 = uv_1[int(3)] * make_float2 (1.0f + r2_23 * (_S195 + r2_23 * (_S196 + r2_23 * (_S197 + r2_23 * _S198)))) + make_float2 (_S205 * u_23 * v_23 + _S200 * (r2_23 + 2.0f * u_23 * u_23) + _S201 * r2_23, _S206 * u_23 * v_23 + _S199 * (r2_23 + 2.0f * v_23 * v_23) + _S202 * r2_23);
                float2  _S332 = _S331 + make_float2 (_S203 * _S331.x + _S204 * _S331.y, 0.0f);
                uv_1[int(3)] = make_float2 (fx_1 * _S332.x + cx_1, fy_1 * _S332.y + cy_1);
                break;
            }
            bool all_valid_10 = all_valid_9 & (!_S214);
            float3  _S333 = pos_c_1[int(4)];
            _S215 = &uv_1[int(4)];
            for(;;)
            {
                float2  _S334 = float2 {_S333.x, _S333.y};
                float r_6 = length_0(_S334);
                float _S335 = _S333.z;
                float theta_4 = (F32_atan2((r_6), (_S335)));
                if(theta_4 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S335;
                }
                else
                {
                    k_0 = theta_4 / r_6;
                }
                float2  _S336 = _S334 * make_float2 (k_0);
                uv_1[int(4)] = _S336;
                float u_24 = _S336.x;
                float v_24 = _S336.y;
                float _S337 = u_24 + u_24;
                float r2_24 = u_24 * u_24 + v_24 * v_24;
                float _S338 = _S197 + r2_24 * _S198;
                float _S339 = _S196 + r2_24 * _S338;
                float _S340 = _S195 + r2_24 * _S339;
                float radial_12 = 1.0f + r2_24 * _S340;
                float _S341 = 2.0f * u_24;
                float _S342 = 2.0f * v_24;
                float2  _S343 = _S194 * make_float2 (radial_12) + make_float2 (_S337 * _S340 + (_S337 * _S339 + (_S337 * _S338 + _S337 * _S198 * r2_24) * r2_24) * r2_24) * _S336 + make_float2 (_S205 * v_24 + (_S337 + (_S341 + _S341)) * _S200 + _S337 * _S201, _S206 * v_24 + _S337 * _S199 + _S337 * _S202);
                float _S344 = v_24 + v_24;
                float2  _S345 = _S207 * make_float2 (radial_12) + make_float2 (_S344 * _S340 + (_S344 * _S339 + (_S344 * _S338 + _S344 * _S198 * r2_24) * r2_24) * r2_24) * _S336 + make_float2 (_S205 * u_24 + _S344 * _S200 + _S344 * _S201, _S206 * u_24 + (_S344 + (_S342 + _S342)) * _S199 + _S344 * _S202);
                Matrix<float, 2, 2>  _S346 = transpose_0(makeMatrix<float, 2, 2> (_S343 + make_float2 (_S343.x * _S203 + _S343.y * _S204, 0.0f), _S345 + make_float2 (_S345.x * _S203 + _S345.y * _S204, 0.0f)));
                bool _S347 = !((F32_min((determinant_0(_S346)), ((F32_min((_S346.rows[int(0)].x), (_S346.rows[int(1)].y)))))) > 0.0f);
                _S216 = _S347;
                if(_S347)
                {
                    break;
                }
                float u_25 = uv_1[int(4)].x;
                float v_25 = uv_1[int(4)].y;
                float r2_25 = u_25 * u_25 + v_25 * v_25;
                float2  _S348 = uv_1[int(4)] * make_float2 (1.0f + r2_25 * (_S195 + r2_25 * (_S196 + r2_25 * (_S197 + r2_25 * _S198)))) + make_float2 (_S205 * u_25 * v_25 + _S200 * (r2_25 + 2.0f * u_25 * u_25) + _S201 * r2_25, _S206 * u_25 * v_25 + _S199 * (r2_25 + 2.0f * v_25 * v_25) + _S202 * r2_25);
                float2  _S349 = _S348 + make_float2 (_S203 * _S348.x + _S204 * _S348.y, 0.0f);
                uv_1[int(4)] = make_float2 (fx_1 * _S349.x + cx_1, fy_1 * _S349.y + cy_1);
                break;
            }
            bool all_valid_11 = all_valid_10 & (!_S216);
            float3  _S350 = pos_c_1[int(5)];
            _S217 = &uv_1[int(5)];
            for(;;)
            {
                float2  _S351 = float2 {_S350.x, _S350.y};
                float r_7 = length_0(_S351);
                float _S352 = _S350.z;
                float theta_5 = (F32_atan2((r_7), (_S352)));
                if(theta_5 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_5 * theta_5 / 3.0f) / _S352;
                }
                else
                {
                    k_0 = theta_5 / r_7;
                }
                float2  _S353 = _S351 * make_float2 (k_0);
                uv_1[int(5)] = _S353;
                float u_26 = _S353.x;
                float v_26 = _S353.y;
                float _S354 = u_26 + u_26;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float _S355 = _S197 + r2_26 * _S198;
                float _S356 = _S196 + r2_26 * _S355;
                float _S357 = _S195 + r2_26 * _S356;
                float radial_13 = 1.0f + r2_26 * _S357;
                float _S358 = 2.0f * u_26;
                float _S359 = 2.0f * v_26;
                float2  _S360 = _S194 * make_float2 (radial_13) + make_float2 (_S354 * _S357 + (_S354 * _S356 + (_S354 * _S355 + _S354 * _S198 * r2_26) * r2_26) * r2_26) * _S353 + make_float2 (_S205 * v_26 + (_S354 + (_S358 + _S358)) * _S200 + _S354 * _S201, _S206 * v_26 + _S354 * _S199 + _S354 * _S202);
                float _S361 = v_26 + v_26;
                float2  _S362 = _S207 * make_float2 (radial_13) + make_float2 (_S361 * _S357 + (_S361 * _S356 + (_S361 * _S355 + _S361 * _S198 * r2_26) * r2_26) * r2_26) * _S353 + make_float2 (_S205 * u_26 + _S361 * _S200 + _S361 * _S201, _S206 * u_26 + (_S361 + (_S359 + _S359)) * _S199 + _S361 * _S202);
                Matrix<float, 2, 2>  _S363 = transpose_0(makeMatrix<float, 2, 2> (_S360 + make_float2 (_S360.x * _S203 + _S360.y * _S204, 0.0f), _S362 + make_float2 (_S362.x * _S203 + _S362.y * _S204, 0.0f)));
                bool _S364 = !((F32_min((determinant_0(_S363)), ((F32_min((_S363.rows[int(0)].x), (_S363.rows[int(1)].y)))))) > 0.0f);
                _S218 = _S364;
                if(_S364)
                {
                    break;
                }
                float u_27 = uv_1[int(5)].x;
                float v_27 = uv_1[int(5)].y;
                float r2_27 = u_27 * u_27 + v_27 * v_27;
                float2  _S365 = uv_1[int(5)] * make_float2 (1.0f + r2_27 * (_S195 + r2_27 * (_S196 + r2_27 * (_S197 + r2_27 * _S198)))) + make_float2 (_S205 * u_27 * v_27 + _S200 * (r2_27 + 2.0f * u_27 * u_27) + _S201 * r2_27, _S206 * u_27 * v_27 + _S199 * (r2_27 + 2.0f * v_27 * v_27) + _S202 * r2_27);
                float2  _S366 = _S365 + make_float2 (_S203 * _S365.x + _S204 * _S365.y, 0.0f);
                uv_1[int(5)] = make_float2 (fx_1 * _S366.x + cx_1, fy_1 * _S366.y + cy_1);
                break;
            }
            bool all_valid_12 = all_valid_11 & (!_S218);
            float3  _S367 = pos_c_1[int(6)];
            _S219 = &uv_1[int(6)];
            for(;;)
            {
                float2  _S368 = float2 {_S367.x, _S367.y};
                float r_8 = length_0(_S368);
                float _S369 = _S367.z;
                float theta_6 = (F32_atan2((r_8), (_S369)));
                if(theta_6 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_6 * theta_6 / 3.0f) / _S369;
                }
                else
                {
                    k_0 = theta_6 / r_8;
                }
                float2  _S370 = _S368 * make_float2 (k_0);
                uv_1[int(6)] = _S370;
                float u_28 = _S370.x;
                float v_28 = _S370.y;
                float _S371 = u_28 + u_28;
                float r2_28 = u_28 * u_28 + v_28 * v_28;
                float _S372 = _S197 + r2_28 * _S198;
                float _S373 = _S196 + r2_28 * _S372;
                float _S374 = _S195 + r2_28 * _S373;
                float radial_14 = 1.0f + r2_28 * _S374;
                float _S375 = 2.0f * u_28;
                float _S376 = 2.0f * v_28;
                float2  _S377 = _S194 * make_float2 (radial_14) + make_float2 (_S371 * _S374 + (_S371 * _S373 + (_S371 * _S372 + _S371 * _S198 * r2_28) * r2_28) * r2_28) * _S370 + make_float2 (_S205 * v_28 + (_S371 + (_S375 + _S375)) * _S200 + _S371 * _S201, _S206 * v_28 + _S371 * _S199 + _S371 * _S202);
                float _S378 = v_28 + v_28;
                float2  _S379 = _S207 * make_float2 (radial_14) + make_float2 (_S378 * _S374 + (_S378 * _S373 + (_S378 * _S372 + _S378 * _S198 * r2_28) * r2_28) * r2_28) * _S370 + make_float2 (_S205 * u_28 + _S378 * _S200 + _S378 * _S201, _S206 * u_28 + (_S378 + (_S376 + _S376)) * _S199 + _S378 * _S202);
                Matrix<float, 2, 2>  _S380 = transpose_0(makeMatrix<float, 2, 2> (_S377 + make_float2 (_S377.x * _S203 + _S377.y * _S204, 0.0f), _S379 + make_float2 (_S379.x * _S203 + _S379.y * _S204, 0.0f)));
                bool _S381 = !((F32_min((determinant_0(_S380)), ((F32_min((_S380.rows[int(0)].x), (_S380.rows[int(1)].y)))))) > 0.0f);
                _S220 = _S381;
                if(_S381)
                {
                    break;
                }
                float u_29 = uv_1[int(6)].x;
                float v_29 = uv_1[int(6)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S382 = uv_1[int(6)] * make_float2 (1.0f + r2_29 * (_S195 + r2_29 * (_S196 + r2_29 * (_S197 + r2_29 * _S198)))) + make_float2 (_S205 * u_29 * v_29 + _S200 * (r2_29 + 2.0f * u_29 * u_29) + _S201 * r2_29, _S206 * u_29 * v_29 + _S199 * (r2_29 + 2.0f * v_29 * v_29) + _S202 * r2_29);
                float2  _S383 = _S382 + make_float2 (_S203 * _S382.x + _S204 * _S382.y, 0.0f);
                uv_1[int(6)] = make_float2 (fx_1 * _S383.x + cx_1, fy_1 * _S383.y + cy_1);
                break;
            }
            bool all_valid_13 = all_valid_12 & (!_S220);
            float3  _S384 = pos_c_1[int(7)];
            _S221 = &uv_1[int(7)];
            for(;;)
            {
                float2  _S385 = float2 {_S384.x, _S384.y};
                float r_9 = length_0(_S385);
                float _S386 = _S384.z;
                float theta_7 = (F32_atan2((r_9), (_S386)));
                if(theta_7 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_7 * theta_7 / 3.0f) / _S386;
                }
                else
                {
                    k_0 = theta_7 / r_9;
                }
                float2  _S387 = _S385 * make_float2 (k_0);
                uv_1[int(7)] = _S387;
                float u_30 = _S387.x;
                float v_30 = _S387.y;
                float _S388 = u_30 + u_30;
                float r2_30 = u_30 * u_30 + v_30 * v_30;
                float _S389 = _S197 + r2_30 * _S198;
                float _S390 = _S196 + r2_30 * _S389;
                float _S391 = _S195 + r2_30 * _S390;
                float radial_15 = 1.0f + r2_30 * _S391;
                float _S392 = 2.0f * u_30;
                float _S393 = 2.0f * v_30;
                float2  _S394 = _S194 * make_float2 (radial_15) + make_float2 (_S388 * _S391 + (_S388 * _S390 + (_S388 * _S389 + _S388 * _S198 * r2_30) * r2_30) * r2_30) * _S387 + make_float2 (_S205 * v_30 + (_S388 + (_S392 + _S392)) * _S200 + _S388 * _S201, _S206 * v_30 + _S388 * _S199 + _S388 * _S202);
                float _S395 = v_30 + v_30;
                float2  _S396 = _S207 * make_float2 (radial_15) + make_float2 (_S395 * _S391 + (_S395 * _S390 + (_S395 * _S389 + _S395 * _S198 * r2_30) * r2_30) * r2_30) * _S387 + make_float2 (_S205 * u_30 + _S395 * _S200 + _S395 * _S201, _S206 * u_30 + (_S395 + (_S393 + _S393)) * _S199 + _S395 * _S202);
                Matrix<float, 2, 2>  _S397 = transpose_0(makeMatrix<float, 2, 2> (_S394 + make_float2 (_S394.x * _S203 + _S394.y * _S204, 0.0f), _S396 + make_float2 (_S396.x * _S203 + _S396.y * _S204, 0.0f)));
                bool _S398 = !((F32_min((determinant_0(_S397)), ((F32_min((_S397.rows[int(0)].x), (_S397.rows[int(1)].y)))))) > 0.0f);
                _S222 = _S398;
                if(_S398)
                {
                    break;
                }
                float u_31 = uv_1[int(7)].x;
                float v_31 = uv_1[int(7)].y;
                float r2_31 = u_31 * u_31 + v_31 * v_31;
                float2  _S399 = uv_1[int(7)] * make_float2 (1.0f + r2_31 * (_S195 + r2_31 * (_S196 + r2_31 * (_S197 + r2_31 * _S198)))) + make_float2 (_S205 * u_31 * v_31 + _S200 * (r2_31 + 2.0f * u_31 * u_31) + _S201 * r2_31, _S206 * u_31 * v_31 + _S199 * (r2_31 + 2.0f * v_31 * v_31) + _S202 * r2_31);
                float2  _S400 = _S399 + make_float2 (_S203 * _S399.x + _S204 * _S399.y, 0.0f);
                uv_1[int(7)] = make_float2 (fx_1 * _S400.x + cx_1, fy_1 * _S400.y + cy_1);
                break;
            }
            _S223 = all_valid_13 & (!_S222);
            break;
        }
        if(!_S223)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((densities_1[int(0)]), (densities_1[int(1)])))), (densities_1[int(2)])))), (densities_1[int(3)])))), (densities_1[int(4)])))), (densities_1[int(5)])))), (densities_1[int(6)])))), (densities_1[int(7)]))) * size_1 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S401 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S193).x), ((*_S209).x)))), ((*_S211).x)))), ((*_S213).x)))), ((*_S215).x)))), ((*_S217).x)))), ((*_S219).x)))), ((*_S221).x)));
        float _S402 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S193).x), ((*_S209).x)))), ((*_S211).x)))), ((*_S213).x)))), ((*_S215).x)))), ((*_S217).x)))), ((*_S219).x)))), ((*_S221).x)));
        float _S403 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S193).y), ((*_S209).y)))), ((*_S211).y)))), ((*_S213).y)))), ((*_S215).y)))), ((*_S217).y)))), ((*_S219).y)))), ((*_S221).y)));
        float _S404 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S193).y), ((*_S209).y)))), ((*_S211).y)))), ((*_S213).y)))), ((*_S215).y)))), ((*_S217).y)))), ((*_S219).y)))), ((*_S221).y)));
        if(_S401 <= 0.0f)
        {
            _S256 = true;
        }
        else
        {
            _S256 = _S402 >= float(image_width_1);
        }
        if(_S256)
        {
            _S256 = true;
        }
        else
        {
            _S256 = _S403 <= 0.0f;
        }
        if(_S256)
        {
            _S256 = true;
        }
        else
        {
            _S256 = _S404 >= float(image_height_1);
        }
        if(_S256)
        {
            _S256 = true;
        }
        else
        {
            if(_S254 <= 0.0f)
            {
                if(_S402 <= 0.0f)
                {
                    _S256 = _S401 >= float(image_width_1);
                }
                else
                {
                    _S256 = false;
                }
                if(_S256)
                {
                    _S256 = true;
                }
                else
                {
                    if(_S404 <= 0.0f)
                    {
                        _S256 = _S403 >= float(image_width_1);
                    }
                    else
                    {
                        _S256 = false;
                    }
                }
            }
            else
            {
                _S256 = false;
            }
        }
        if(_S256)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (float(int((F32_floor((_S402))))), float(int((F32_floor((_S404))))), float(int((F32_ceil((_S401))))), float(int((F32_ceil((_S403))))));
        float _S405 = mean_c_1.z;
        *depth_1 = (F32_max((_S405), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S405 + length_1(mean_c_1)))));
        float3  _S406 = mean_c_1 - - mul_0(transpose_1(R_1), t_1);
        float _S407 = _S406.x;
        float _S408 = _S406.y;
        float _S409 = _S406.z;
        float norm_1 = (F32_sqrt((_S407 * _S407 + _S408 * _S408 + _S409 * _S409)));
        float x_8 = _S407 / norm_1;
        float y_4 = _S408 / norm_1;
        float z_1 = _S409 / norm_1;
        float z2_1 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_8 * x_8 - y_4 * y_4;
        float fS1_1 = 2.0f * x_8 * y_4;
        float fTmp0C_1 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgbs_1 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * sh_coeffs_1[int(1)] + make_float3 (z_1) * sh_coeffs_1[int(2)] - make_float3 (x_8) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_4) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_8) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_8 * fS1_1 + y_4 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_4) * sh_coeffs_1[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_8) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_8 * fC1_1 - y_4 * fS1_1)) * sh_coeffs_1[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  densities_2, FixedArray<float3 , 16>  sh_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, float4  * aabb_xyxy_2, float * depth_2, float3  * rgbs_2)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S410 = mul_0(R_2, pos_2) + t_2;
    pos_c_2[int(0)] = _S410;
    float3  _S411 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_2;
    pos_c_2[int(1)] = _S411;
    float3  _S412 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(2)] = _S412;
    float3  _S413 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(3)] = _S413;
    float3  _S414 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(4)] = _S414;
    float3  _S415 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(5)] = _S415;
    float3  _S416 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_2;
    pos_c_2[int(6)] = _S416;
    float3  _S417 = mul_0(R_2, pos_2 + make_float3 (size_2)) + t_2;
    pos_c_2[int(7)] = _S417;
    float3  mean_c_2 = mul_0(R_2, pos_2 + make_float3 (0.5f * size_2)) + t_2;
    FixedArray<float2 , 8>  uv_2;
    float2  _S418 = float2 {_S410.x, _S410.y} / make_float2 (_S410.z);
    float u_32 = _S418.x;
    float v_32 = _S418.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float _S419 = 2.0f * dist_coeffs_2[int(4)];
    float _S420 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S421 = _S418 * make_float2 (1.0f + r2_32 * (dist_coeffs_2[int(0)] + r2_32 * (dist_coeffs_2[int(1)] + r2_32 * (dist_coeffs_2[int(2)] + r2_32 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_32 * v_32 + dist_coeffs_2[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + dist_coeffs_2[int(6)] * r2_32, _S420 * u_32 * v_32 + dist_coeffs_2[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + dist_coeffs_2[int(7)] * r2_32);
    float2  _S422 = _S421 + make_float2 (dist_coeffs_2[int(8)] * _S421.x + dist_coeffs_2[int(9)] * _S421.y, 0.0f);
    float _S423 = fx_2 * _S422.x + cx_2;
    float _S424 = fy_2 * _S422.y + cy_2;
    uv_2[int(0)] = make_float2 (_S423, _S424);
    float2  _S425 = float2 {_S411.x, _S411.y} / make_float2 (_S411.z);
    float u_33 = _S425.x;
    float v_33 = _S425.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S426 = _S425 * make_float2 (1.0f + r2_33 * (dist_coeffs_2[int(0)] + r2_33 * (dist_coeffs_2[int(1)] + r2_33 * (dist_coeffs_2[int(2)] + r2_33 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_33 * v_33 + dist_coeffs_2[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + dist_coeffs_2[int(6)] * r2_33, _S420 * u_33 * v_33 + dist_coeffs_2[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + dist_coeffs_2[int(7)] * r2_33);
    float2  _S427 = _S426 + make_float2 (dist_coeffs_2[int(8)] * _S426.x + dist_coeffs_2[int(9)] * _S426.y, 0.0f);
    float _S428 = fx_2 * _S427.x + cx_2;
    float _S429 = fy_2 * _S427.y + cy_2;
    uv_2[int(1)] = make_float2 (_S428, _S429);
    float2  _S430 = float2 {_S412.x, _S412.y} / make_float2 (_S412.z);
    float u_34 = _S430.x;
    float v_34 = _S430.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S431 = _S430 * make_float2 (1.0f + r2_34 * (dist_coeffs_2[int(0)] + r2_34 * (dist_coeffs_2[int(1)] + r2_34 * (dist_coeffs_2[int(2)] + r2_34 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_34 * v_34 + dist_coeffs_2[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + dist_coeffs_2[int(6)] * r2_34, _S420 * u_34 * v_34 + dist_coeffs_2[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + dist_coeffs_2[int(7)] * r2_34);
    float2  _S432 = _S431 + make_float2 (dist_coeffs_2[int(8)] * _S431.x + dist_coeffs_2[int(9)] * _S431.y, 0.0f);
    float _S433 = fx_2 * _S432.x + cx_2;
    float _S434 = fy_2 * _S432.y + cy_2;
    uv_2[int(2)] = make_float2 (_S433, _S434);
    float2  _S435 = float2 {_S413.x, _S413.y} / make_float2 (_S413.z);
    float u_35 = _S435.x;
    float v_35 = _S435.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S436 = _S435 * make_float2 (1.0f + r2_35 * (dist_coeffs_2[int(0)] + r2_35 * (dist_coeffs_2[int(1)] + r2_35 * (dist_coeffs_2[int(2)] + r2_35 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_35 * v_35 + dist_coeffs_2[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + dist_coeffs_2[int(6)] * r2_35, _S420 * u_35 * v_35 + dist_coeffs_2[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + dist_coeffs_2[int(7)] * r2_35);
    float2  _S437 = _S436 + make_float2 (dist_coeffs_2[int(8)] * _S436.x + dist_coeffs_2[int(9)] * _S436.y, 0.0f);
    float _S438 = fx_2 * _S437.x + cx_2;
    float _S439 = fy_2 * _S437.y + cy_2;
    uv_2[int(3)] = make_float2 (_S438, _S439);
    float2  _S440 = float2 {_S414.x, _S414.y} / make_float2 (_S414.z);
    float u_36 = _S440.x;
    float v_36 = _S440.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S441 = _S440 * make_float2 (1.0f + r2_36 * (dist_coeffs_2[int(0)] + r2_36 * (dist_coeffs_2[int(1)] + r2_36 * (dist_coeffs_2[int(2)] + r2_36 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_36 * v_36 + dist_coeffs_2[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + dist_coeffs_2[int(6)] * r2_36, _S420 * u_36 * v_36 + dist_coeffs_2[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + dist_coeffs_2[int(7)] * r2_36);
    float2  _S442 = _S441 + make_float2 (dist_coeffs_2[int(8)] * _S441.x + dist_coeffs_2[int(9)] * _S441.y, 0.0f);
    float _S443 = fx_2 * _S442.x + cx_2;
    float _S444 = fy_2 * _S442.y + cy_2;
    uv_2[int(4)] = make_float2 (_S443, _S444);
    float2  _S445 = float2 {_S415.x, _S415.y} / make_float2 (_S415.z);
    float u_37 = _S445.x;
    float v_37 = _S445.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S446 = _S445 * make_float2 (1.0f + r2_37 * (dist_coeffs_2[int(0)] + r2_37 * (dist_coeffs_2[int(1)] + r2_37 * (dist_coeffs_2[int(2)] + r2_37 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_37 * v_37 + dist_coeffs_2[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + dist_coeffs_2[int(6)] * r2_37, _S420 * u_37 * v_37 + dist_coeffs_2[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + dist_coeffs_2[int(7)] * r2_37);
    float2  _S447 = _S446 + make_float2 (dist_coeffs_2[int(8)] * _S446.x + dist_coeffs_2[int(9)] * _S446.y, 0.0f);
    float _S448 = fx_2 * _S447.x + cx_2;
    float _S449 = fy_2 * _S447.y + cy_2;
    uv_2[int(5)] = make_float2 (_S448, _S449);
    float2  _S450 = float2 {_S416.x, _S416.y} / make_float2 (_S416.z);
    float u_38 = _S450.x;
    float v_38 = _S450.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float2  _S451 = _S450 * make_float2 (1.0f + r2_38 * (dist_coeffs_2[int(0)] + r2_38 * (dist_coeffs_2[int(1)] + r2_38 * (dist_coeffs_2[int(2)] + r2_38 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_38 * v_38 + dist_coeffs_2[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + dist_coeffs_2[int(6)] * r2_38, _S420 * u_38 * v_38 + dist_coeffs_2[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + dist_coeffs_2[int(7)] * r2_38);
    float2  _S452 = _S451 + make_float2 (dist_coeffs_2[int(8)] * _S451.x + dist_coeffs_2[int(9)] * _S451.y, 0.0f);
    float _S453 = fx_2 * _S452.x + cx_2;
    float _S454 = fy_2 * _S452.y + cy_2;
    uv_2[int(6)] = make_float2 (_S453, _S454);
    float2  _S455 = float2 {_S417.x, _S417.y} / make_float2 (_S417.z);
    float u_39 = _S455.x;
    float v_39 = _S455.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float2  _S456 = _S455 * make_float2 (1.0f + r2_39 * (dist_coeffs_2[int(0)] + r2_39 * (dist_coeffs_2[int(1)] + r2_39 * (dist_coeffs_2[int(2)] + r2_39 * dist_coeffs_2[int(3)])))) + make_float2 (_S419 * u_39 * v_39 + dist_coeffs_2[int(5)] * (r2_39 + 2.0f * u_39 * u_39) + dist_coeffs_2[int(6)] * r2_39, _S420 * u_39 * v_39 + dist_coeffs_2[int(4)] * (r2_39 + 2.0f * v_39 * v_39) + dist_coeffs_2[int(7)] * r2_39);
    float2  _S457 = _S456 + make_float2 (dist_coeffs_2[int(8)] * _S456.x + dist_coeffs_2[int(9)] * _S456.y, 0.0f);
    float _S458 = fx_2 * _S457.x + cx_2;
    float _S459 = fy_2 * _S457.y + cy_2;
    uv_2[int(7)] = make_float2 (_S458, _S459);
    *aabb_xyxy_2 = make_float4 (float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S423), (_S428)))), (_S433)))), (_S438)))), (_S443)))), (_S448)))), (_S453)))), (_S458)))))))), float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S424), (_S429)))), (_S434)))), (_S439)))), (_S444)))), (_S449)))), (_S454)))), (_S459)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S423), (_S428)))), (_S433)))), (_S438)))), (_S443)))), (_S448)))), (_S453)))), (_S458)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S424), (_S429)))), (_S434)))), (_S439)))), (_S444)))), (_S449)))), (_S454)))), (_S459)))))))));
    *depth_2 = mean_c_2.z;
    float3  _S460 = mean_c_2 - - mul_0(transpose_1(R_2), t_2);
    float _S461 = _S460.x;
    float _S462 = _S460.y;
    float _S463 = _S460.z;
    float norm_2 = (F32_sqrt((_S461 * _S461 + _S462 * _S462 + _S463 * _S463)));
    float x_9 = _S461 / norm_2;
    float y_5 = _S462 / norm_2;
    float z_2 = _S463 / norm_2;
    float z2_2 = z_2 * z_2;
    float fTmp0B_2 = -1.09254848957061768f * z_2;
    float fC1_2 = x_9 * x_9 - y_5 * y_5;
    float fS1_2 = 2.0f * x_9 * y_5;
    float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_2;
    *rgbs_2 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_2[int(1)] + make_float3 (z_2) * sh_coeffs_2[int(2)] - make_float3 (x_9) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_5) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_2 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_9) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_2 + y_5 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_5) * sh_coeffs_2[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_2 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_9) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_2 - y_5 * fS1_2)) * sh_coeffs_2[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  densities_3, FixedArray<float3 , 16>  sh_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, float4  * aabb_xyxy_3, float * depth_3, float3  * rgbs_3)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S464 = mul_0(R_3, pos_3) + t_3;
    pos_c_3[int(0)] = _S464;
    pos_c_3[int(1)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_3;
    pos_c_3[int(2)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(3)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(4)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(5)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(6)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_3;
    pos_c_3[int(7)] = mul_0(R_3, pos_3 + make_float3 (size_3)) + t_3;
    float3  mean_c_3 = mul_0(R_3, pos_3 + make_float3 (0.5f * size_3)) + t_3;
    FixedArray<float2 , 8>  uv_3;
    float2  _S465 = float2 {_S464.x, _S464.y};
    float r_10 = length_0(_S465);
    float _S466 = _S464.z;
    float theta_8 = (F32_atan2((r_10), (_S466)));
    float k_1;
    if(theta_8 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_8 * theta_8 / 3.0f) / _S466;
    }
    else
    {
        k_1 = theta_8 / r_10;
    }
    float2  _S467 = _S465 * make_float2 (k_1);
    float u_40 = _S467.x;
    float v_40 = _S467.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S468 = 2.0f * dist_coeffs_3[int(4)];
    float _S469 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S470 = _S467 * make_float2 (1.0f + r2_40 * (dist_coeffs_3[int(0)] + r2_40 * (dist_coeffs_3[int(1)] + r2_40 * (dist_coeffs_3[int(2)] + r2_40 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_40 * v_40 + dist_coeffs_3[int(5)] * (r2_40 + 2.0f * u_40 * u_40) + dist_coeffs_3[int(6)] * r2_40, _S469 * u_40 * v_40 + dist_coeffs_3[int(4)] * (r2_40 + 2.0f * v_40 * v_40) + dist_coeffs_3[int(7)] * r2_40);
    float2  _S471 = _S470 + make_float2 (dist_coeffs_3[int(8)] * _S470.x + dist_coeffs_3[int(9)] * _S470.y, 0.0f);
    uv_3[int(0)] = make_float2 (fx_3 * _S471.x + cx_3, fy_3 * _S471.y + cy_3);
    float2  _S472 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_11 = length_0(_S472);
    float _S473 = pos_c_3[int(1)].z;
    float theta_9 = (F32_atan2((r_11), (_S473)));
    if(theta_9 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_9 * theta_9 / 3.0f) / _S473;
    }
    else
    {
        k_1 = theta_9 / r_11;
    }
    float2  _S474 = _S472 * make_float2 (k_1);
    float u_41 = _S474.x;
    float v_41 = _S474.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float2  _S475 = _S474 * make_float2 (1.0f + r2_41 * (dist_coeffs_3[int(0)] + r2_41 * (dist_coeffs_3[int(1)] + r2_41 * (dist_coeffs_3[int(2)] + r2_41 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_41 * v_41 + dist_coeffs_3[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + dist_coeffs_3[int(6)] * r2_41, _S469 * u_41 * v_41 + dist_coeffs_3[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + dist_coeffs_3[int(7)] * r2_41);
    float2  _S476 = _S475 + make_float2 (dist_coeffs_3[int(8)] * _S475.x + dist_coeffs_3[int(9)] * _S475.y, 0.0f);
    uv_3[int(1)] = make_float2 (fx_3 * _S476.x + cx_3, fy_3 * _S476.y + cy_3);
    float2  _S477 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_12 = length_0(_S477);
    float _S478 = pos_c_3[int(2)].z;
    float theta_10 = (F32_atan2((r_12), (_S478)));
    if(theta_10 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_10 * theta_10 / 3.0f) / _S478;
    }
    else
    {
        k_1 = theta_10 / r_12;
    }
    float2  _S479 = _S477 * make_float2 (k_1);
    float u_42 = _S479.x;
    float v_42 = _S479.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float2  _S480 = _S479 * make_float2 (1.0f + r2_42 * (dist_coeffs_3[int(0)] + r2_42 * (dist_coeffs_3[int(1)] + r2_42 * (dist_coeffs_3[int(2)] + r2_42 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_42 * v_42 + dist_coeffs_3[int(5)] * (r2_42 + 2.0f * u_42 * u_42) + dist_coeffs_3[int(6)] * r2_42, _S469 * u_42 * v_42 + dist_coeffs_3[int(4)] * (r2_42 + 2.0f * v_42 * v_42) + dist_coeffs_3[int(7)] * r2_42);
    float2  _S481 = _S480 + make_float2 (dist_coeffs_3[int(8)] * _S480.x + dist_coeffs_3[int(9)] * _S480.y, 0.0f);
    uv_3[int(2)] = make_float2 (fx_3 * _S481.x + cx_3, fy_3 * _S481.y + cy_3);
    float2  _S482 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_13 = length_0(_S482);
    float _S483 = pos_c_3[int(3)].z;
    float theta_11 = (F32_atan2((r_13), (_S483)));
    if(theta_11 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_11 * theta_11 / 3.0f) / _S483;
    }
    else
    {
        k_1 = theta_11 / r_13;
    }
    float2  _S484 = _S482 * make_float2 (k_1);
    float u_43 = _S484.x;
    float v_43 = _S484.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float2  _S485 = _S484 * make_float2 (1.0f + r2_43 * (dist_coeffs_3[int(0)] + r2_43 * (dist_coeffs_3[int(1)] + r2_43 * (dist_coeffs_3[int(2)] + r2_43 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_43 * v_43 + dist_coeffs_3[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + dist_coeffs_3[int(6)] * r2_43, _S469 * u_43 * v_43 + dist_coeffs_3[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + dist_coeffs_3[int(7)] * r2_43);
    float2  _S486 = _S485 + make_float2 (dist_coeffs_3[int(8)] * _S485.x + dist_coeffs_3[int(9)] * _S485.y, 0.0f);
    uv_3[int(3)] = make_float2 (fx_3 * _S486.x + cx_3, fy_3 * _S486.y + cy_3);
    float2  _S487 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_14 = length_0(_S487);
    float _S488 = pos_c_3[int(4)].z;
    float theta_12 = (F32_atan2((r_14), (_S488)));
    if(theta_12 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_12 * theta_12 / 3.0f) / _S488;
    }
    else
    {
        k_1 = theta_12 / r_14;
    }
    float2  _S489 = _S487 * make_float2 (k_1);
    float u_44 = _S489.x;
    float v_44 = _S489.y;
    float r2_44 = u_44 * u_44 + v_44 * v_44;
    float2  _S490 = _S489 * make_float2 (1.0f + r2_44 * (dist_coeffs_3[int(0)] + r2_44 * (dist_coeffs_3[int(1)] + r2_44 * (dist_coeffs_3[int(2)] + r2_44 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_44 * v_44 + dist_coeffs_3[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + dist_coeffs_3[int(6)] * r2_44, _S469 * u_44 * v_44 + dist_coeffs_3[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + dist_coeffs_3[int(7)] * r2_44);
    float2  _S491 = _S490 + make_float2 (dist_coeffs_3[int(8)] * _S490.x + dist_coeffs_3[int(9)] * _S490.y, 0.0f);
    uv_3[int(4)] = make_float2 (fx_3 * _S491.x + cx_3, fy_3 * _S491.y + cy_3);
    float2  _S492 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_15 = length_0(_S492);
    float _S493 = pos_c_3[int(5)].z;
    float theta_13 = (F32_atan2((r_15), (_S493)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_13 * theta_13 / 3.0f) / _S493;
    }
    else
    {
        k_1 = theta_13 / r_15;
    }
    float2  _S494 = _S492 * make_float2 (k_1);
    float u_45 = _S494.x;
    float v_45 = _S494.y;
    float r2_45 = u_45 * u_45 + v_45 * v_45;
    float2  _S495 = _S494 * make_float2 (1.0f + r2_45 * (dist_coeffs_3[int(0)] + r2_45 * (dist_coeffs_3[int(1)] + r2_45 * (dist_coeffs_3[int(2)] + r2_45 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_45 * v_45 + dist_coeffs_3[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + dist_coeffs_3[int(6)] * r2_45, _S469 * u_45 * v_45 + dist_coeffs_3[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + dist_coeffs_3[int(7)] * r2_45);
    float2  _S496 = _S495 + make_float2 (dist_coeffs_3[int(8)] * _S495.x + dist_coeffs_3[int(9)] * _S495.y, 0.0f);
    uv_3[int(5)] = make_float2 (fx_3 * _S496.x + cx_3, fy_3 * _S496.y + cy_3);
    float2  _S497 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_16 = length_0(_S497);
    float _S498 = pos_c_3[int(6)].z;
    float theta_14 = (F32_atan2((r_16), (_S498)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_14 * theta_14 / 3.0f) / _S498;
    }
    else
    {
        k_1 = theta_14 / r_16;
    }
    float2  _S499 = _S497 * make_float2 (k_1);
    float u_46 = _S499.x;
    float v_46 = _S499.y;
    float r2_46 = u_46 * u_46 + v_46 * v_46;
    float2  _S500 = _S499 * make_float2 (1.0f + r2_46 * (dist_coeffs_3[int(0)] + r2_46 * (dist_coeffs_3[int(1)] + r2_46 * (dist_coeffs_3[int(2)] + r2_46 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_46 * v_46 + dist_coeffs_3[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + dist_coeffs_3[int(6)] * r2_46, _S469 * u_46 * v_46 + dist_coeffs_3[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + dist_coeffs_3[int(7)] * r2_46);
    float2  _S501 = _S500 + make_float2 (dist_coeffs_3[int(8)] * _S500.x + dist_coeffs_3[int(9)] * _S500.y, 0.0f);
    uv_3[int(6)] = make_float2 (fx_3 * _S501.x + cx_3, fy_3 * _S501.y + cy_3);
    float2  _S502 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_17 = length_0(_S502);
    float _S503 = pos_c_3[int(7)].z;
    float theta_15 = (F32_atan2((r_17), (_S503)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_15 * theta_15 / 3.0f) / _S503;
    }
    else
    {
        k_1 = theta_15 / r_17;
    }
    float2  _S504 = _S502 * make_float2 (k_1);
    float u_47 = _S504.x;
    float v_47 = _S504.y;
    float r2_47 = u_47 * u_47 + v_47 * v_47;
    float2  _S505 = _S504 * make_float2 (1.0f + r2_47 * (dist_coeffs_3[int(0)] + r2_47 * (dist_coeffs_3[int(1)] + r2_47 * (dist_coeffs_3[int(2)] + r2_47 * dist_coeffs_3[int(3)])))) + make_float2 (_S468 * u_47 * v_47 + dist_coeffs_3[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + dist_coeffs_3[int(6)] * r2_47, _S469 * u_47 * v_47 + dist_coeffs_3[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + dist_coeffs_3[int(7)] * r2_47);
    float2  _S506 = _S505 + make_float2 (dist_coeffs_3[int(8)] * _S505.x + dist_coeffs_3[int(9)] * _S505.y, 0.0f);
    float _S507 = fx_3 * _S506.x + cx_3;
    float _S508 = fy_3 * _S506.y + cy_3;
    uv_3[int(7)] = make_float2 (_S507, _S508);
    *aabb_xyxy_3 = make_float4 (float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S507)))))))), float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S508)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S507)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S508)))))))));
    float _S509 = mean_c_3.z;
    *depth_3 = (F32_max((_S509), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S509 + length_1(mean_c_3)))));
    float3  _S510 = mean_c_3 - - mul_0(transpose_1(R_3), t_3);
    float _S511 = _S510.x;
    float _S512 = _S510.y;
    float _S513 = _S510.z;
    float norm_3 = (F32_sqrt((_S511 * _S511 + _S512 * _S512 + _S513 * _S513)));
    float x_10 = _S511 / norm_3;
    float y_6 = _S512 / norm_3;
    float z_3 = _S513 / norm_3;
    float z2_3 = z_3 * z_3;
    float fTmp0B_3 = -1.09254848957061768f * z_3;
    float fC1_3 = x_10 * x_10 - y_6 * y_6;
    float fS1_3 = 2.0f * x_10 * y_6;
    float fTmp0C_3 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_3;
    *rgbs_3 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_3[int(1)] + make_float3 (z_3) * sh_coeffs_3[int(2)] - make_float3 (x_10) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_6) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_10) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_10 * fS1_3 + y_6 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_6) * sh_coeffs_3[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_10) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_10 * fC1_3 - y_6 * fS1_3)) * sh_coeffs_3[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S514, float3  _S515)
{
    return mul_0(_S514, _S515);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S516)
{
    return (F32_sqrt((_S516)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S517, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S518, float3  _S519)
{
    _d_max_vector_0(_S517, _S518, _S519);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S520, float _S521)
{
    _d_sqrt_0(_S520, _S521);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S522, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S523, float3  _S524)
{
    _d_mul_0(_S522, _S523, _S524);
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  densities_4, FixedArray<float3 , 16>  sh_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float3  v_rgb_0, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  _S525 = s_primal_ctx_mul_0(R_4, pos_4) + t_4;
    float _S526 = _S525.z;
    float2  _S527 = make_float2 (_S526);
    float _S528 = (F32_min((1.00000001504746622e+30f), (_S526)));
    float _S529 = (F32_max((0.0f), (_S526)));
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S530 = s_primal_ctx_mul_0(R_4, pos_i_0) + t_4;
    float _S531 = _S530.z;
    float2  _S532 = make_float2 (_S531);
    float _S533 = (F32_min((_S528), (_S531)));
    float _S534 = (F32_max((_S529), (_S531)));
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S535 = s_primal_ctx_mul_0(R_4, pos_i_1) + t_4;
    float _S536 = _S535.z;
    float2  _S537 = make_float2 (_S536);
    float _S538 = (F32_min((_S533), (_S536)));
    float _S539 = (F32_max((_S534), (_S536)));
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S540 = s_primal_ctx_mul_0(R_4, pos_i_2) + t_4;
    float _S541 = _S540.z;
    float2  _S542 = make_float2 (_S541);
    float _S543 = (F32_min((_S538), (_S541)));
    float _S544 = (F32_max((_S539), (_S541)));
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S545 = s_primal_ctx_mul_0(R_4, pos_i_3) + t_4;
    float _S546 = _S545.z;
    float2  _S547 = make_float2 (_S546);
    float _S548 = (F32_min((_S543), (_S546)));
    float _S549 = (F32_max((_S544), (_S546)));
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S550 = s_primal_ctx_mul_0(R_4, pos_i_4) + t_4;
    float _S551 = _S550.z;
    float2  _S552 = make_float2 (_S551);
    float _S553 = (F32_min((_S548), (_S551)));
    float _S554 = (F32_max((_S549), (_S551)));
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S555 = s_primal_ctx_mul_0(R_4, pos_i_5) + t_4;
    float _S556 = _S555.z;
    float2  _S557 = make_float2 (_S556);
    float _S558 = (F32_min((_S553), (_S556)));
    float _S559 = (F32_max((_S554), (_S556)));
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S560 = s_primal_ctx_mul_0(R_4, pos_i_6) + t_4;
    float _S561 = _S560.z;
    float2  _S562 = make_float2 (_S561);
    float3  _S563 = pos_4 + make_float3 (0.5f * size_4);
    float2  _S564 = float2 {_S525.x, _S525.y};
    float2  _S565 = _S564 / make_float2 (_S526);
    float2  _S566 = make_float2 (_S526 * _S526);
    float u_48 = _S565.x;
    float v_48 = _S565.y;
    float r2_48 = u_48 * u_48 + v_48 * v_48;
    float _S567 = dist_coeffs_4[int(2)] + r2_48 * dist_coeffs_4[int(3)];
    float _S568 = dist_coeffs_4[int(1)] + r2_48 * _S567;
    float _S569 = dist_coeffs_4[int(0)] + r2_48 * _S568;
    float radial_16 = 1.0f + r2_48 * _S569;
    float _S570 = 2.0f * dist_coeffs_4[int(4)];
    float _S571 = _S570 * u_48;
    float _S572 = 2.0f * u_48;
    float _S573 = 2.0f * dist_coeffs_4[int(5)];
    float _S574 = _S573 * u_48;
    float _S575 = 2.0f * v_48;
    float2  _S576 = _S565 * make_float2 (radial_16) + make_float2 (_S571 * v_48 + dist_coeffs_4[int(5)] * (r2_48 + _S572 * u_48) + dist_coeffs_4[int(6)] * r2_48, _S574 * v_48 + dist_coeffs_4[int(4)] * (r2_48 + _S575 * v_48) + dist_coeffs_4[int(7)] * r2_48);
    float2  _S577 = _S576 + make_float2 (dist_coeffs_4[int(8)] * _S576.x + dist_coeffs_4[int(9)] * _S576.y, 0.0f);
    float _S578 = fx_4 * _S577.x + cx_4;
    float _S579 = fy_4 * _S577.y + cy_4;
    float2  _S580 = float2 {_S530.x, _S530.y};
    float2  _S581 = _S580 / make_float2 (_S531);
    float2  _S582 = make_float2 (_S531 * _S531);
    float u_49 = _S581.x;
    float v_49 = _S581.y;
    float r2_49 = u_49 * u_49 + v_49 * v_49;
    float _S583 = dist_coeffs_4[int(2)] + r2_49 * dist_coeffs_4[int(3)];
    float _S584 = dist_coeffs_4[int(1)] + r2_49 * _S583;
    float _S585 = dist_coeffs_4[int(0)] + r2_49 * _S584;
    float radial_17 = 1.0f + r2_49 * _S585;
    float _S586 = _S570 * u_49;
    float _S587 = 2.0f * u_49;
    float _S588 = _S573 * u_49;
    float _S589 = 2.0f * v_49;
    float2  _S590 = _S581 * make_float2 (radial_17) + make_float2 (_S586 * v_49 + dist_coeffs_4[int(5)] * (r2_49 + _S587 * u_49) + dist_coeffs_4[int(6)] * r2_49, _S588 * v_49 + dist_coeffs_4[int(4)] * (r2_49 + _S589 * v_49) + dist_coeffs_4[int(7)] * r2_49);
    float2  _S591 = _S590 + make_float2 (dist_coeffs_4[int(8)] * _S590.x + dist_coeffs_4[int(9)] * _S590.y, 0.0f);
    float _S592 = fx_4 * _S591.x + cx_4;
    float _S593 = fy_4 * _S591.y + cy_4;
    float2  _S594 = float2 {_S535.x, _S535.y};
    float2  _S595 = _S594 / make_float2 (_S536);
    float2  _S596 = make_float2 (_S536 * _S536);
    float u_50 = _S595.x;
    float v_50 = _S595.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S597 = dist_coeffs_4[int(2)] + r2_50 * dist_coeffs_4[int(3)];
    float _S598 = dist_coeffs_4[int(1)] + r2_50 * _S597;
    float _S599 = dist_coeffs_4[int(0)] + r2_50 * _S598;
    float radial_18 = 1.0f + r2_50 * _S599;
    float _S600 = _S570 * u_50;
    float _S601 = 2.0f * u_50;
    float _S602 = _S573 * u_50;
    float _S603 = 2.0f * v_50;
    float2  _S604 = _S595 * make_float2 (radial_18) + make_float2 (_S600 * v_50 + dist_coeffs_4[int(5)] * (r2_50 + _S601 * u_50) + dist_coeffs_4[int(6)] * r2_50, _S602 * v_50 + dist_coeffs_4[int(4)] * (r2_50 + _S603 * v_50) + dist_coeffs_4[int(7)] * r2_50);
    float2  _S605 = _S604 + make_float2 (dist_coeffs_4[int(8)] * _S604.x + dist_coeffs_4[int(9)] * _S604.y, 0.0f);
    float _S606 = fx_4 * _S605.x + cx_4;
    float _S607 = fy_4 * _S605.y + cy_4;
    float2  _S608 = float2 {_S540.x, _S540.y};
    float2  _S609 = _S608 / make_float2 (_S541);
    float2  _S610 = make_float2 (_S541 * _S541);
    float u_51 = _S609.x;
    float v_51 = _S609.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float _S611 = dist_coeffs_4[int(2)] + r2_51 * dist_coeffs_4[int(3)];
    float _S612 = dist_coeffs_4[int(1)] + r2_51 * _S611;
    float _S613 = dist_coeffs_4[int(0)] + r2_51 * _S612;
    float radial_19 = 1.0f + r2_51 * _S613;
    float _S614 = _S570 * u_51;
    float _S615 = 2.0f * u_51;
    float _S616 = _S573 * u_51;
    float _S617 = 2.0f * v_51;
    float2  _S618 = _S609 * make_float2 (radial_19) + make_float2 (_S614 * v_51 + dist_coeffs_4[int(5)] * (r2_51 + _S615 * u_51) + dist_coeffs_4[int(6)] * r2_51, _S616 * v_51 + dist_coeffs_4[int(4)] * (r2_51 + _S617 * v_51) + dist_coeffs_4[int(7)] * r2_51);
    float2  _S619 = _S618 + make_float2 (dist_coeffs_4[int(8)] * _S618.x + dist_coeffs_4[int(9)] * _S618.y, 0.0f);
    float _S620 = fx_4 * _S619.x + cx_4;
    float _S621 = fy_4 * _S619.y + cy_4;
    float2  _S622 = float2 {_S545.x, _S545.y};
    float2  _S623 = _S622 / make_float2 (_S546);
    float2  _S624 = make_float2 (_S546 * _S546);
    float u_52 = _S623.x;
    float v_52 = _S623.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float _S625 = dist_coeffs_4[int(2)] + r2_52 * dist_coeffs_4[int(3)];
    float _S626 = dist_coeffs_4[int(1)] + r2_52 * _S625;
    float _S627 = dist_coeffs_4[int(0)] + r2_52 * _S626;
    float radial_20 = 1.0f + r2_52 * _S627;
    float _S628 = _S570 * u_52;
    float _S629 = 2.0f * u_52;
    float _S630 = _S573 * u_52;
    float _S631 = 2.0f * v_52;
    float2  _S632 = _S623 * make_float2 (radial_20) + make_float2 (_S628 * v_52 + dist_coeffs_4[int(5)] * (r2_52 + _S629 * u_52) + dist_coeffs_4[int(6)] * r2_52, _S630 * v_52 + dist_coeffs_4[int(4)] * (r2_52 + _S631 * v_52) + dist_coeffs_4[int(7)] * r2_52);
    float2  _S633 = _S632 + make_float2 (dist_coeffs_4[int(8)] * _S632.x + dist_coeffs_4[int(9)] * _S632.y, 0.0f);
    float _S634 = fx_4 * _S633.x + cx_4;
    float _S635 = fy_4 * _S633.y + cy_4;
    float2  _S636 = float2 {_S550.x, _S550.y};
    float2  _S637 = _S636 / make_float2 (_S551);
    float2  _S638 = make_float2 (_S551 * _S551);
    float u_53 = _S637.x;
    float v_53 = _S637.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S639 = dist_coeffs_4[int(2)] + r2_53 * dist_coeffs_4[int(3)];
    float _S640 = dist_coeffs_4[int(1)] + r2_53 * _S639;
    float _S641 = dist_coeffs_4[int(0)] + r2_53 * _S640;
    float radial_21 = 1.0f + r2_53 * _S641;
    float _S642 = _S570 * u_53;
    float _S643 = 2.0f * u_53;
    float _S644 = _S573 * u_53;
    float _S645 = 2.0f * v_53;
    float2  _S646 = _S637 * make_float2 (radial_21) + make_float2 (_S642 * v_53 + dist_coeffs_4[int(5)] * (r2_53 + _S643 * u_53) + dist_coeffs_4[int(6)] * r2_53, _S644 * v_53 + dist_coeffs_4[int(4)] * (r2_53 + _S645 * v_53) + dist_coeffs_4[int(7)] * r2_53);
    float2  _S647 = _S646 + make_float2 (dist_coeffs_4[int(8)] * _S646.x + dist_coeffs_4[int(9)] * _S646.y, 0.0f);
    float _S648 = fx_4 * _S647.x + cx_4;
    float _S649 = fy_4 * _S647.y + cy_4;
    float2  _S650 = float2 {_S555.x, _S555.y};
    float2  _S651 = _S650 / make_float2 (_S556);
    float2  _S652 = make_float2 (_S556 * _S556);
    float u_54 = _S651.x;
    float v_54 = _S651.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float _S653 = dist_coeffs_4[int(2)] + r2_54 * dist_coeffs_4[int(3)];
    float _S654 = dist_coeffs_4[int(1)] + r2_54 * _S653;
    float _S655 = dist_coeffs_4[int(0)] + r2_54 * _S654;
    float radial_22 = 1.0f + r2_54 * _S655;
    float _S656 = _S570 * u_54;
    float _S657 = 2.0f * u_54;
    float _S658 = _S573 * u_54;
    float _S659 = 2.0f * v_54;
    float2  _S660 = _S651 * make_float2 (radial_22) + make_float2 (_S656 * v_54 + dist_coeffs_4[int(5)] * (r2_54 + _S657 * u_54) + dist_coeffs_4[int(6)] * r2_54, _S658 * v_54 + dist_coeffs_4[int(4)] * (r2_54 + _S659 * v_54) + dist_coeffs_4[int(7)] * r2_54);
    float2  _S661 = _S660 + make_float2 (dist_coeffs_4[int(8)] * _S660.x + dist_coeffs_4[int(9)] * _S660.y, 0.0f);
    float _S662 = fx_4 * _S661.x + cx_4;
    float _S663 = fy_4 * _S661.y + cy_4;
    float2  _S664 = float2 {_S560.x, _S560.y};
    float2  _S665 = _S664 / make_float2 (_S561);
    float2  _S666 = make_float2 (_S561 * _S561);
    float u_55 = _S665.x;
    float v_55 = _S665.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float _S667 = dist_coeffs_4[int(2)] + r2_55 * dist_coeffs_4[int(3)];
    float _S668 = dist_coeffs_4[int(1)] + r2_55 * _S667;
    float _S669 = dist_coeffs_4[int(0)] + r2_55 * _S668;
    float radial_23 = 1.0f + r2_55 * _S669;
    float _S670 = _S570 * u_55;
    float _S671 = 2.0f * u_55;
    float _S672 = _S573 * u_55;
    float _S673 = 2.0f * v_55;
    float2  _S674 = _S665 * make_float2 (radial_23) + make_float2 (_S670 * v_55 + dist_coeffs_4[int(5)] * (r2_55 + _S671 * u_55) + dist_coeffs_4[int(6)] * r2_55, _S672 * v_55 + dist_coeffs_4[int(4)] * (r2_55 + _S673 * v_55) + dist_coeffs_4[int(7)] * r2_55);
    float2  _S675 = _S674 + make_float2 (dist_coeffs_4[int(8)] * _S674.x + dist_coeffs_4[int(9)] * _S674.y, 0.0f);
    float _S676 = fx_4 * _S675.x + cx_4;
    float _S677 = fy_4 * _S675.y + cy_4;
    float _S678 = (F32_max((_S578), (_S592)));
    float _S679 = (F32_min((_S578), (_S592)));
    float _S680 = (F32_max((_S579), (_S593)));
    float _S681 = (F32_min((_S579), (_S593)));
    float _S682 = (F32_max((_S678), (_S606)));
    float _S683 = (F32_min((_S679), (_S606)));
    float _S684 = (F32_max((_S680), (_S607)));
    float _S685 = (F32_min((_S681), (_S607)));
    float _S686 = (F32_max((_S682), (_S620)));
    float _S687 = (F32_min((_S683), (_S620)));
    float _S688 = (F32_max((_S684), (_S621)));
    float _S689 = (F32_min((_S685), (_S621)));
    float _S690 = (F32_max((_S686), (_S634)));
    float _S691 = (F32_min((_S687), (_S634)));
    float _S692 = (F32_max((_S688), (_S635)));
    float _S693 = (F32_min((_S689), (_S635)));
    float _S694 = (F32_max((_S690), (_S648)));
    float _S695 = (F32_min((_S691), (_S648)));
    float _S696 = (F32_max((_S692), (_S649)));
    float _S697 = (F32_min((_S693), (_S649)));
    float _S698 = (F32_max((_S694), (_S662)));
    float _S699 = (F32_min((_S695), (_S662)));
    float _S700 = (F32_max((_S696), (_S663)));
    float _S701 = (F32_min((_S697), (_S663)));
    Matrix<float, 3, 3>  _S702 = transpose_1(R_4);
    float3  _S703 = s_primal_ctx_mul_0(R_4, _S563) + t_4 - - s_primal_ctx_mul_0(_S702, t_4);
    float _S704 = _S703.x;
    float _S705 = _S703.y;
    float _S706 = _S703.z;
    float _S707 = _S704 * _S704 + _S705 * _S705 + _S706 * _S706;
    float _S708 = s_primal_ctx_sqrt_0(_S707);
    float x_11 = _S704 / _S708;
    float3  _S709 = make_float3 (x_11);
    float _S710 = _S708 * _S708;
    float y_7 = _S705 / _S708;
    float z_4 = _S706 / _S708;
    float3  _S711 = make_float3 (z_4);
    float _S712 = - y_7;
    float3  _S713 = make_float3 (_S712);
    float z2_4 = z_4 * z_4;
    float fTmp0B_4 = -1.09254848957061768f * z_4;
    float fC1_4 = x_11 * x_11 - y_7 * y_7;
    float _S714 = 2.0f * x_11;
    float fS1_4 = _S714 * y_7;
    float pSH6_0 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float3  _S715 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_11;
    float3  _S716 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_7;
    float3  _S717 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S718 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S719 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_4;
    float _S720 = 1.86588168144226074f * z2_4 - 1.11952900886535645f;
    float pSH12_0 = z_4 * _S720;
    float3  _S721 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_11;
    float3  _S722 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_7;
    float3  _S723 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S724 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S725 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_11 * fC1_4 - y_7 * fS1_4);
    float3  _S726 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_11 * fS1_4 + y_7 * fC1_4);
    float3  _S727 = make_float3 (pSH9_0);
    float3  _S728 = make_float3 (0.0f);
    float3  _S729 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S730;
    (&_S730)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S712) * sh_coeffs_4[int(1)] + make_float3 (z_4) * sh_coeffs_4[int(2)] - make_float3 (x_11) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]) + make_float3 (0.5f);
    (&_S730)->differential_0 = _S729;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S731;
    (&_S731)->primal_0 = _S728;
    (&_S731)->differential_0 = _S729;
    s_bwd_prop_max_0(&_S730, &_S731, v_rgb_0);
    float3  _S732 = _S726 * _S730.differential_0;
    float3  _S733 = sh_coeffs_4[int(15)] * _S730.differential_0;
    float3  _S734 = _S724 * _S730.differential_0;
    float3  _S735 = sh_coeffs_4[int(14)] * _S730.differential_0;
    float3  _S736 = _S722 * _S730.differential_0;
    float3  _S737 = sh_coeffs_4[int(13)] * _S730.differential_0;
    float3  _S738 = _S721 * _S730.differential_0;
    float3  _S739 = sh_coeffs_4[int(12)] * _S730.differential_0;
    float3  _S740 = _S723 * _S730.differential_0;
    float3  _S741 = sh_coeffs_4[int(11)] * _S730.differential_0;
    float3  _S742 = _S725 * _S730.differential_0;
    float3  _S743 = sh_coeffs_4[int(10)] * _S730.differential_0;
    float3  _S744 = _S727 * _S730.differential_0;
    float3  _S745 = sh_coeffs_4[int(9)] * _S730.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S745.x + _S745.y + _S745.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S733.x + _S733.y + _S733.z);
    float _S746 = _S743.x + _S743.y + _S743.z;
    float _S747 = _S735.x + _S735.y + _S735.z;
    float _S748 = _S741.x + _S741.y + _S741.z;
    float _S749 = _S737.x + _S737.y + _S737.z;
    float _S750 = _S739.x + _S739.y + _S739.z;
    float _S751 = - s_diff_fC2_T_0;
    float3  _S752 = _S718 * _S730.differential_0;
    float3  _S753 = sh_coeffs_4[int(8)] * _S730.differential_0;
    float3  _S754 = _S716 * _S730.differential_0;
    float3  _S755 = sh_coeffs_4[int(7)] * _S730.differential_0;
    float3  _S756 = _S715 * _S730.differential_0;
    float3  _S757 = sh_coeffs_4[int(6)] * _S730.differential_0;
    float3  _S758 = _S717 * _S730.differential_0;
    float3  _S759 = sh_coeffs_4[int(5)] * _S730.differential_0;
    float3  _S760 = _S719 * _S730.differential_0;
    float3  _S761 = sh_coeffs_4[int(4)] * _S730.differential_0;
    float _S762 = _S759.x + _S759.y + _S759.z;
    float _S763 = _S755.x + _S755.y + _S755.z;
    float _S764 = fTmp1B_4 * _S746 + x_11 * s_diff_fS2_T_0 + y_7 * _S751 + 0.54627424478530884f * (_S761.x + _S761.y + _S761.z);
    float _S765 = fTmp1B_4 * _S747 + y_7 * s_diff_fS2_T_0 + x_11 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S753.x + _S753.y + _S753.z);
    float _S766 = y_7 * - _S765;
    float _S767 = x_11 * _S765;
    float _S768 = z_4 * (1.86588168144226074f * (z_4 * _S750) + -2.28522896766662598f * (y_7 * _S748 + x_11 * _S749) + 0.94617468118667603f * (_S757.x + _S757.y + _S757.z));
    float3  _S769 = make_float3 (0.48860251903533936f) * _S730.differential_0;
    float3  _S770 = - _S769;
    float3  _S771 = _S709 * _S770;
    float3  _S772 = sh_coeffs_4[int(3)] * _S770;
    float3  _S773 = _S711 * _S769;
    float3  _S774 = sh_coeffs_4[int(2)] * _S769;
    float3  _S775 = _S713 * _S769;
    float3  _S776 = sh_coeffs_4[int(1)] * _S769;
    float _S777 = (_S720 * _S750 + 1.44530570507049561f * (fS1_4 * _S746 + fC1_4 * _S747) + -1.09254848957061768f * (y_7 * _S762 + x_11 * _S763) + _S768 + _S768 + _S774.x + _S774.y + _S774.z) / _S710;
    float _S778 = _S708 * _S777;
    float _S779 = (fTmp0C_4 * _S748 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S751 + fTmp0B_4 * _S762 + _S714 * _S764 + _S766 + _S766 + - (_S776.x + _S776.y + _S776.z)) / _S710;
    float _S780 = _S708 * _S779;
    float _S781 = (fTmp0C_4 * _S749 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S763 + 2.0f * (y_7 * _S764) + _S767 + _S767 + _S772.x + _S772.y + _S772.z) / _S710;
    float _S782 = _S708 * _S781;
    float _S783 = _S706 * - _S777 + _S705 * - _S779 + _S704 * - _S781;
    DiffPair_float_0 _S784;
    (&_S784)->primal_0 = _S707;
    (&_S784)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S784, _S783);
    float _S785 = _S706 * _S784.differential_0;
    float _S786 = _S705 * _S784.differential_0;
    float _S787 = _S704 * _S784.differential_0;
    float3  _S788 = make_float3 (0.282094806432724f) * _S730.differential_0;
    float3  _S789 = make_float3 (_S782 + _S787 + _S787, _S780 + _S786 + _S786, _S778 + _S785 + _S785);
    float3  _S790 = - - _S789;
    Matrix<float, 3, 3>  _S791 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S792;
    (&_S792)->primal_0 = _S702;
    (&_S792)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S793;
    (&_S793)->primal_0 = t_4;
    (&_S793)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S792, &_S793, _S790);
    Matrix<float, 3, 3>  _S794 = transpose_1(_S792.differential_0);
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = _S701;
    (&_S795)->differential_0 = 0.0f;
    DiffPair_float_0 _S796;
    (&_S796)->primal_0 = _S677;
    (&_S796)->differential_0 = 0.0f;
    _d_min_0(&_S795, &_S796, 0.0f);
    DiffPair_float_0 _S797;
    (&_S797)->primal_0 = _S700;
    (&_S797)->differential_0 = 0.0f;
    DiffPair_float_0 _S798;
    (&_S798)->primal_0 = _S677;
    (&_S798)->differential_0 = 0.0f;
    _d_max_0(&_S797, &_S798, 0.0f);
    DiffPair_float_0 _S799;
    (&_S799)->primal_0 = _S699;
    (&_S799)->differential_0 = 0.0f;
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = _S676;
    (&_S800)->differential_0 = 0.0f;
    _d_min_0(&_S799, &_S800, 0.0f);
    DiffPair_float_0 _S801;
    (&_S801)->primal_0 = _S698;
    (&_S801)->differential_0 = 0.0f;
    DiffPair_float_0 _S802;
    (&_S802)->primal_0 = _S676;
    (&_S802)->differential_0 = 0.0f;
    _d_max_0(&_S801, &_S802, 0.0f);
    DiffPair_float_0 _S803;
    (&_S803)->primal_0 = _S697;
    (&_S803)->differential_0 = 0.0f;
    DiffPair_float_0 _S804;
    (&_S804)->primal_0 = _S663;
    (&_S804)->differential_0 = 0.0f;
    _d_min_0(&_S803, &_S804, _S795.differential_0);
    DiffPair_float_0 _S805;
    (&_S805)->primal_0 = _S696;
    (&_S805)->differential_0 = 0.0f;
    DiffPair_float_0 _S806;
    (&_S806)->primal_0 = _S663;
    (&_S806)->differential_0 = 0.0f;
    _d_max_0(&_S805, &_S806, _S797.differential_0);
    DiffPair_float_0 _S807;
    (&_S807)->primal_0 = _S695;
    (&_S807)->differential_0 = 0.0f;
    DiffPair_float_0 _S808;
    (&_S808)->primal_0 = _S662;
    (&_S808)->differential_0 = 0.0f;
    _d_min_0(&_S807, &_S808, _S799.differential_0);
    DiffPair_float_0 _S809;
    (&_S809)->primal_0 = _S694;
    (&_S809)->differential_0 = 0.0f;
    DiffPair_float_0 _S810;
    (&_S810)->primal_0 = _S662;
    (&_S810)->differential_0 = 0.0f;
    _d_max_0(&_S809, &_S810, _S801.differential_0);
    DiffPair_float_0 _S811;
    (&_S811)->primal_0 = _S693;
    (&_S811)->differential_0 = 0.0f;
    DiffPair_float_0 _S812;
    (&_S812)->primal_0 = _S649;
    (&_S812)->differential_0 = 0.0f;
    _d_min_0(&_S811, &_S812, _S803.differential_0);
    DiffPair_float_0 _S813;
    (&_S813)->primal_0 = _S692;
    (&_S813)->differential_0 = 0.0f;
    DiffPair_float_0 _S814;
    (&_S814)->primal_0 = _S649;
    (&_S814)->differential_0 = 0.0f;
    _d_max_0(&_S813, &_S814, _S805.differential_0);
    DiffPair_float_0 _S815;
    (&_S815)->primal_0 = _S691;
    (&_S815)->differential_0 = 0.0f;
    DiffPair_float_0 _S816;
    (&_S816)->primal_0 = _S648;
    (&_S816)->differential_0 = 0.0f;
    _d_min_0(&_S815, &_S816, _S807.differential_0);
    DiffPair_float_0 _S817;
    (&_S817)->primal_0 = _S690;
    (&_S817)->differential_0 = 0.0f;
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = _S648;
    (&_S818)->differential_0 = 0.0f;
    _d_max_0(&_S817, &_S818, _S809.differential_0);
    DiffPair_float_0 _S819;
    (&_S819)->primal_0 = _S689;
    (&_S819)->differential_0 = 0.0f;
    DiffPair_float_0 _S820;
    (&_S820)->primal_0 = _S635;
    (&_S820)->differential_0 = 0.0f;
    _d_min_0(&_S819, &_S820, _S811.differential_0);
    DiffPair_float_0 _S821;
    (&_S821)->primal_0 = _S688;
    (&_S821)->differential_0 = 0.0f;
    DiffPair_float_0 _S822;
    (&_S822)->primal_0 = _S635;
    (&_S822)->differential_0 = 0.0f;
    _d_max_0(&_S821, &_S822, _S813.differential_0);
    DiffPair_float_0 _S823;
    (&_S823)->primal_0 = _S687;
    (&_S823)->differential_0 = 0.0f;
    DiffPair_float_0 _S824;
    (&_S824)->primal_0 = _S634;
    (&_S824)->differential_0 = 0.0f;
    _d_min_0(&_S823, &_S824, _S815.differential_0);
    DiffPair_float_0 _S825;
    (&_S825)->primal_0 = _S686;
    (&_S825)->differential_0 = 0.0f;
    DiffPair_float_0 _S826;
    (&_S826)->primal_0 = _S634;
    (&_S826)->differential_0 = 0.0f;
    _d_max_0(&_S825, &_S826, _S817.differential_0);
    DiffPair_float_0 _S827;
    (&_S827)->primal_0 = _S685;
    (&_S827)->differential_0 = 0.0f;
    DiffPair_float_0 _S828;
    (&_S828)->primal_0 = _S621;
    (&_S828)->differential_0 = 0.0f;
    _d_min_0(&_S827, &_S828, _S819.differential_0);
    DiffPair_float_0 _S829;
    (&_S829)->primal_0 = _S684;
    (&_S829)->differential_0 = 0.0f;
    DiffPair_float_0 _S830;
    (&_S830)->primal_0 = _S621;
    (&_S830)->differential_0 = 0.0f;
    _d_max_0(&_S829, &_S830, _S821.differential_0);
    DiffPair_float_0 _S831;
    (&_S831)->primal_0 = _S683;
    (&_S831)->differential_0 = 0.0f;
    DiffPair_float_0 _S832;
    (&_S832)->primal_0 = _S620;
    (&_S832)->differential_0 = 0.0f;
    _d_min_0(&_S831, &_S832, _S823.differential_0);
    DiffPair_float_0 _S833;
    (&_S833)->primal_0 = _S682;
    (&_S833)->differential_0 = 0.0f;
    DiffPair_float_0 _S834;
    (&_S834)->primal_0 = _S620;
    (&_S834)->differential_0 = 0.0f;
    _d_max_0(&_S833, &_S834, _S825.differential_0);
    DiffPair_float_0 _S835;
    (&_S835)->primal_0 = _S681;
    (&_S835)->differential_0 = 0.0f;
    DiffPair_float_0 _S836;
    (&_S836)->primal_0 = _S607;
    (&_S836)->differential_0 = 0.0f;
    _d_min_0(&_S835, &_S836, _S827.differential_0);
    DiffPair_float_0 _S837;
    (&_S837)->primal_0 = _S680;
    (&_S837)->differential_0 = 0.0f;
    DiffPair_float_0 _S838;
    (&_S838)->primal_0 = _S607;
    (&_S838)->differential_0 = 0.0f;
    _d_max_0(&_S837, &_S838, _S829.differential_0);
    DiffPair_float_0 _S839;
    (&_S839)->primal_0 = _S679;
    (&_S839)->differential_0 = 0.0f;
    DiffPair_float_0 _S840;
    (&_S840)->primal_0 = _S606;
    (&_S840)->differential_0 = 0.0f;
    _d_min_0(&_S839, &_S840, _S831.differential_0);
    DiffPair_float_0 _S841;
    (&_S841)->primal_0 = _S678;
    (&_S841)->differential_0 = 0.0f;
    DiffPair_float_0 _S842;
    (&_S842)->primal_0 = _S606;
    (&_S842)->differential_0 = 0.0f;
    _d_max_0(&_S841, &_S842, _S833.differential_0);
    DiffPair_float_0 _S843;
    (&_S843)->primal_0 = _S579;
    (&_S843)->differential_0 = 0.0f;
    DiffPair_float_0 _S844;
    (&_S844)->primal_0 = _S593;
    (&_S844)->differential_0 = 0.0f;
    _d_min_0(&_S843, &_S844, _S835.differential_0);
    DiffPair_float_0 _S845;
    (&_S845)->primal_0 = _S579;
    (&_S845)->differential_0 = 0.0f;
    DiffPair_float_0 _S846;
    (&_S846)->primal_0 = _S593;
    (&_S846)->differential_0 = 0.0f;
    _d_max_0(&_S845, &_S846, _S837.differential_0);
    DiffPair_float_0 _S847;
    (&_S847)->primal_0 = _S578;
    (&_S847)->differential_0 = 0.0f;
    DiffPair_float_0 _S848;
    (&_S848)->primal_0 = _S592;
    (&_S848)->differential_0 = 0.0f;
    _d_min_0(&_S847, &_S848, _S839.differential_0);
    DiffPair_float_0 _S849;
    (&_S849)->primal_0 = _S578;
    (&_S849)->differential_0 = 0.0f;
    DiffPair_float_0 _S850;
    (&_S850)->primal_0 = _S592;
    (&_S850)->differential_0 = 0.0f;
    _d_max_0(&_S849, &_S850, _S841.differential_0);
    float _S851 = fx_4 * (_S800.differential_0 + _S802.differential_0);
    float2  _S852 = make_float2 (_S851, fy_4 * (_S796.differential_0 + _S798.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S851, dist_coeffs_4[int(9)] * _S851);
    float2  _S853 = _S665 * _S852;
    float _S854 = dist_coeffs_4[int(4)] * _S852.y;
    float _S855 = dist_coeffs_4[int(5)] * _S852.x;
    float _S856 = _S853.x + _S853.y;
    float _S857 = r2_55 * _S856;
    float _S858 = r2_55 * _S857;
    float _S859 = dist_coeffs_4[int(7)] * _S852.y + _S854 + dist_coeffs_4[int(6)] * _S852.x + _S855 + _S669 * _S856 + _S668 * _S857 + _S667 * _S858 + dist_coeffs_4[int(3)] * (r2_55 * _S858);
    float _S860 = v_55 * _S859;
    float _S861 = u_55 * _S859;
    float2  _S862 = (make_float2 (radial_23) * _S852 + make_float2 (_S573 * (v_55 * _S852.y) + _S671 * _S855 + 2.0f * (u_55 * _S855) + _S570 * (v_55 * _S852.x) + _S861 + _S861, _S673 * _S854 + 2.0f * (v_55 * _S854) + _S672 * _S852.y + _S670 * _S852.x + _S860 + _S860)) / _S666;
    float2  _S863 = _S664 * - _S862;
    float2  _S864 = _S562 * _S862;
    float _S865 = fx_4 * (_S808.differential_0 + _S810.differential_0);
    float2  _S866 = make_float2 (_S865, fy_4 * (_S804.differential_0 + _S806.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S865, dist_coeffs_4[int(9)] * _S865);
    float2  _S867 = _S651 * _S866;
    float _S868 = dist_coeffs_4[int(4)] * _S866.y;
    float _S869 = dist_coeffs_4[int(5)] * _S866.x;
    float _S870 = _S867.x + _S867.y;
    float _S871 = r2_54 * _S870;
    float _S872 = r2_54 * _S871;
    float _S873 = dist_coeffs_4[int(7)] * _S866.y + _S868 + dist_coeffs_4[int(6)] * _S866.x + _S869 + _S655 * _S870 + _S654 * _S871 + _S653 * _S872 + dist_coeffs_4[int(3)] * (r2_54 * _S872);
    float _S874 = v_54 * _S873;
    float _S875 = u_54 * _S873;
    float2  _S876 = (make_float2 (radial_22) * _S866 + make_float2 (_S573 * (v_54 * _S866.y) + _S657 * _S869 + 2.0f * (u_54 * _S869) + _S570 * (v_54 * _S866.x) + _S875 + _S875, _S659 * _S868 + 2.0f * (v_54 * _S868) + _S658 * _S866.y + _S656 * _S866.x + _S874 + _S874)) / _S652;
    float2  _S877 = _S650 * - _S876;
    float2  _S878 = _S557 * _S876;
    float _S879 = fx_4 * (_S816.differential_0 + _S818.differential_0);
    float2  _S880 = make_float2 (_S879, fy_4 * (_S812.differential_0 + _S814.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S879, dist_coeffs_4[int(9)] * _S879);
    float2  _S881 = _S637 * _S880;
    float _S882 = dist_coeffs_4[int(4)] * _S880.y;
    float _S883 = dist_coeffs_4[int(5)] * _S880.x;
    float _S884 = _S881.x + _S881.y;
    float _S885 = r2_53 * _S884;
    float _S886 = r2_53 * _S885;
    float _S887 = dist_coeffs_4[int(7)] * _S880.y + _S882 + dist_coeffs_4[int(6)] * _S880.x + _S883 + _S641 * _S884 + _S640 * _S885 + _S639 * _S886 + dist_coeffs_4[int(3)] * (r2_53 * _S886);
    float _S888 = v_53 * _S887;
    float _S889 = u_53 * _S887;
    float2  _S890 = (make_float2 (radial_21) * _S880 + make_float2 (_S573 * (v_53 * _S880.y) + _S643 * _S883 + 2.0f * (u_53 * _S883) + _S570 * (v_53 * _S880.x) + _S889 + _S889, _S645 * _S882 + 2.0f * (v_53 * _S882) + _S644 * _S880.y + _S642 * _S880.x + _S888 + _S888)) / _S638;
    float2  _S891 = _S636 * - _S890;
    float2  _S892 = _S552 * _S890;
    float _S893 = fx_4 * (_S824.differential_0 + _S826.differential_0);
    float2  _S894 = make_float2 (_S893, fy_4 * (_S820.differential_0 + _S822.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S893, dist_coeffs_4[int(9)] * _S893);
    float2  _S895 = _S623 * _S894;
    float _S896 = dist_coeffs_4[int(4)] * _S894.y;
    float _S897 = dist_coeffs_4[int(5)] * _S894.x;
    float _S898 = _S895.x + _S895.y;
    float _S899 = r2_52 * _S898;
    float _S900 = r2_52 * _S899;
    float _S901 = dist_coeffs_4[int(7)] * _S894.y + _S896 + dist_coeffs_4[int(6)] * _S894.x + _S897 + _S627 * _S898 + _S626 * _S899 + _S625 * _S900 + dist_coeffs_4[int(3)] * (r2_52 * _S900);
    float _S902 = v_52 * _S901;
    float _S903 = u_52 * _S901;
    float2  _S904 = (make_float2 (radial_20) * _S894 + make_float2 (_S573 * (v_52 * _S894.y) + _S629 * _S897 + 2.0f * (u_52 * _S897) + _S570 * (v_52 * _S894.x) + _S903 + _S903, _S631 * _S896 + 2.0f * (v_52 * _S896) + _S630 * _S894.y + _S628 * _S894.x + _S902 + _S902)) / _S624;
    float2  _S905 = _S622 * - _S904;
    float2  _S906 = _S547 * _S904;
    float _S907 = fx_4 * (_S832.differential_0 + _S834.differential_0);
    float2  _S908 = make_float2 (_S907, fy_4 * (_S828.differential_0 + _S830.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S907, dist_coeffs_4[int(9)] * _S907);
    float2  _S909 = _S609 * _S908;
    float _S910 = dist_coeffs_4[int(4)] * _S908.y;
    float _S911 = dist_coeffs_4[int(5)] * _S908.x;
    float _S912 = _S909.x + _S909.y;
    float _S913 = r2_51 * _S912;
    float _S914 = r2_51 * _S913;
    float _S915 = dist_coeffs_4[int(7)] * _S908.y + _S910 + dist_coeffs_4[int(6)] * _S908.x + _S911 + _S613 * _S912 + _S612 * _S913 + _S611 * _S914 + dist_coeffs_4[int(3)] * (r2_51 * _S914);
    float _S916 = v_51 * _S915;
    float _S917 = u_51 * _S915;
    float2  _S918 = (make_float2 (radial_19) * _S908 + make_float2 (_S573 * (v_51 * _S908.y) + _S615 * _S911 + 2.0f * (u_51 * _S911) + _S570 * (v_51 * _S908.x) + _S917 + _S917, _S617 * _S910 + 2.0f * (v_51 * _S910) + _S616 * _S908.y + _S614 * _S908.x + _S916 + _S916)) / _S610;
    float2  _S919 = _S608 * - _S918;
    float2  _S920 = _S542 * _S918;
    float _S921 = fx_4 * (_S840.differential_0 + _S842.differential_0);
    float2  _S922 = make_float2 (_S921, fy_4 * (_S836.differential_0 + _S838.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S921, dist_coeffs_4[int(9)] * _S921);
    float2  _S923 = _S595 * _S922;
    float _S924 = dist_coeffs_4[int(4)] * _S922.y;
    float _S925 = dist_coeffs_4[int(5)] * _S922.x;
    float _S926 = _S923.x + _S923.y;
    float _S927 = r2_50 * _S926;
    float _S928 = r2_50 * _S927;
    float _S929 = dist_coeffs_4[int(7)] * _S922.y + _S924 + dist_coeffs_4[int(6)] * _S922.x + _S925 + _S599 * _S926 + _S598 * _S927 + _S597 * _S928 + dist_coeffs_4[int(3)] * (r2_50 * _S928);
    float _S930 = v_50 * _S929;
    float _S931 = u_50 * _S929;
    float2  _S932 = (make_float2 (radial_18) * _S922 + make_float2 (_S573 * (v_50 * _S922.y) + _S601 * _S925 + 2.0f * (u_50 * _S925) + _S570 * (v_50 * _S922.x) + _S931 + _S931, _S603 * _S924 + 2.0f * (v_50 * _S924) + _S602 * _S922.y + _S600 * _S922.x + _S930 + _S930)) / _S596;
    float2  _S933 = _S594 * - _S932;
    float2  _S934 = _S537 * _S932;
    float _S935 = fx_4 * (_S848.differential_0 + _S850.differential_0);
    float2  _S936 = make_float2 (_S935, fy_4 * (_S844.differential_0 + _S846.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S935, dist_coeffs_4[int(9)] * _S935);
    float2  _S937 = _S581 * _S936;
    float _S938 = dist_coeffs_4[int(4)] * _S936.y;
    float _S939 = dist_coeffs_4[int(5)] * _S936.x;
    float _S940 = _S937.x + _S937.y;
    float _S941 = r2_49 * _S940;
    float _S942 = r2_49 * _S941;
    float _S943 = dist_coeffs_4[int(7)] * _S936.y + _S938 + dist_coeffs_4[int(6)] * _S936.x + _S939 + _S585 * _S940 + _S584 * _S941 + _S583 * _S942 + dist_coeffs_4[int(3)] * (r2_49 * _S942);
    float _S944 = v_49 * _S943;
    float _S945 = u_49 * _S943;
    float2  _S946 = (make_float2 (radial_17) * _S936 + make_float2 (_S573 * (v_49 * _S936.y) + _S587 * _S939 + 2.0f * (u_49 * _S939) + _S570 * (v_49 * _S936.x) + _S945 + _S945, _S589 * _S938 + 2.0f * (v_49 * _S938) + _S588 * _S936.y + _S586 * _S936.x + _S944 + _S944)) / _S582;
    float2  _S947 = _S580 * - _S946;
    float2  _S948 = _S532 * _S946;
    float _S949 = fx_4 * (_S847.differential_0 + _S849.differential_0);
    float2  _S950 = make_float2 (_S949, fy_4 * (_S843.differential_0 + _S845.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S949, dist_coeffs_4[int(9)] * _S949);
    float2  _S951 = _S565 * _S950;
    float _S952 = dist_coeffs_4[int(4)] * _S950.y;
    float _S953 = dist_coeffs_4[int(5)] * _S950.x;
    float _S954 = _S951.x + _S951.y;
    float _S955 = r2_48 * _S954;
    float _S956 = r2_48 * _S955;
    float _S957 = dist_coeffs_4[int(7)] * _S950.y + _S952 + dist_coeffs_4[int(6)] * _S950.x + _S953 + _S569 * _S954 + _S568 * _S955 + _S567 * _S956 + dist_coeffs_4[int(3)] * (r2_48 * _S956);
    float _S958 = v_48 * _S957;
    float _S959 = u_48 * _S957;
    float2  _S960 = (make_float2 (radial_16) * _S950 + make_float2 (_S573 * (v_48 * _S950.y) + _S572 * _S953 + 2.0f * (u_48 * _S953) + _S570 * (v_48 * _S950.x) + _S959 + _S959, _S575 * _S952 + 2.0f * (v_48 * _S952) + _S574 * _S950.y + _S571 * _S950.x + _S958 + _S958)) / _S566;
    float2  _S961 = _S564 * - _S960;
    float2  _S962 = _S527 * _S960;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S963;
    (&_S963)->primal_0 = R_4;
    (&_S963)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S964;
    (&_S964)->primal_0 = _S563;
    (&_S964)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S963, &_S964, _S789);
    DiffPair_float_0 _S965;
    (&_S965)->primal_0 = _S559;
    (&_S965)->differential_0 = 0.0f;
    DiffPair_float_0 _S966;
    (&_S966)->primal_0 = _S561;
    (&_S966)->differential_0 = 0.0f;
    _d_max_0(&_S965, &_S966, 0.0f);
    DiffPair_float_0 _S967;
    (&_S967)->primal_0 = _S558;
    (&_S967)->differential_0 = 0.0f;
    DiffPair_float_0 _S968;
    (&_S968)->primal_0 = _S561;
    (&_S968)->differential_0 = 0.0f;
    _d_min_0(&_S967, &_S968, 0.0f);
    float3  _S969 = make_float3 (_S864.x, _S864.y, _S966.differential_0 + _S968.differential_0 + _S863.x + _S863.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S970;
    (&_S970)->primal_0 = R_4;
    (&_S970)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S971;
    (&_S971)->primal_0 = pos_i_6;
    (&_S971)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S970, &_S971, _S969);
    DiffPair_float_0 _S972;
    (&_S972)->primal_0 = _S554;
    (&_S972)->differential_0 = 0.0f;
    DiffPair_float_0 _S973;
    (&_S973)->primal_0 = _S556;
    (&_S973)->differential_0 = 0.0f;
    _d_max_0(&_S972, &_S973, _S965.differential_0);
    DiffPair_float_0 _S974;
    (&_S974)->primal_0 = _S553;
    (&_S974)->differential_0 = 0.0f;
    DiffPair_float_0 _S975;
    (&_S975)->primal_0 = _S556;
    (&_S975)->differential_0 = 0.0f;
    _d_min_0(&_S974, &_S975, _S967.differential_0);
    float3  _S976 = make_float3 (_S878.x, _S878.y, _S973.differential_0 + _S975.differential_0 + _S877.x + _S877.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S977;
    (&_S977)->primal_0 = R_4;
    (&_S977)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S978;
    (&_S978)->primal_0 = pos_i_5;
    (&_S978)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S977, &_S978, _S976);
    DiffPair_float_0 _S979;
    (&_S979)->primal_0 = _S549;
    (&_S979)->differential_0 = 0.0f;
    DiffPair_float_0 _S980;
    (&_S980)->primal_0 = _S551;
    (&_S980)->differential_0 = 0.0f;
    _d_max_0(&_S979, &_S980, _S972.differential_0);
    DiffPair_float_0 _S981;
    (&_S981)->primal_0 = _S548;
    (&_S981)->differential_0 = 0.0f;
    DiffPair_float_0 _S982;
    (&_S982)->primal_0 = _S551;
    (&_S982)->differential_0 = 0.0f;
    _d_min_0(&_S981, &_S982, _S974.differential_0);
    float3  _S983 = make_float3 (_S892.x, _S892.y, _S980.differential_0 + _S982.differential_0 + _S891.x + _S891.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S984;
    (&_S984)->primal_0 = R_4;
    (&_S984)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S985;
    (&_S985)->primal_0 = pos_i_4;
    (&_S985)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S984, &_S985, _S983);
    DiffPair_float_0 _S986;
    (&_S986)->primal_0 = _S544;
    (&_S986)->differential_0 = 0.0f;
    DiffPair_float_0 _S987;
    (&_S987)->primal_0 = _S546;
    (&_S987)->differential_0 = 0.0f;
    _d_max_0(&_S986, &_S987, _S979.differential_0);
    DiffPair_float_0 _S988;
    (&_S988)->primal_0 = _S543;
    (&_S988)->differential_0 = 0.0f;
    DiffPair_float_0 _S989;
    (&_S989)->primal_0 = _S546;
    (&_S989)->differential_0 = 0.0f;
    _d_min_0(&_S988, &_S989, _S981.differential_0);
    float3  _S990 = make_float3 (_S906.x, _S906.y, _S987.differential_0 + _S989.differential_0 + _S905.x + _S905.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S991;
    (&_S991)->primal_0 = R_4;
    (&_S991)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S992;
    (&_S992)->primal_0 = pos_i_3;
    (&_S992)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S991, &_S992, _S990);
    DiffPair_float_0 _S993;
    (&_S993)->primal_0 = _S539;
    (&_S993)->differential_0 = 0.0f;
    DiffPair_float_0 _S994;
    (&_S994)->primal_0 = _S541;
    (&_S994)->differential_0 = 0.0f;
    _d_max_0(&_S993, &_S994, _S986.differential_0);
    DiffPair_float_0 _S995;
    (&_S995)->primal_0 = _S538;
    (&_S995)->differential_0 = 0.0f;
    DiffPair_float_0 _S996;
    (&_S996)->primal_0 = _S541;
    (&_S996)->differential_0 = 0.0f;
    _d_min_0(&_S995, &_S996, _S988.differential_0);
    float3  _S997 = make_float3 (_S920.x, _S920.y, _S994.differential_0 + _S996.differential_0 + _S919.x + _S919.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S998;
    (&_S998)->primal_0 = R_4;
    (&_S998)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S999;
    (&_S999)->primal_0 = pos_i_2;
    (&_S999)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S998, &_S999, _S997);
    DiffPair_float_0 _S1000;
    (&_S1000)->primal_0 = _S534;
    (&_S1000)->differential_0 = 0.0f;
    DiffPair_float_0 _S1001;
    (&_S1001)->primal_0 = _S536;
    (&_S1001)->differential_0 = 0.0f;
    _d_max_0(&_S1000, &_S1001, _S993.differential_0);
    DiffPair_float_0 _S1002;
    (&_S1002)->primal_0 = _S533;
    (&_S1002)->differential_0 = 0.0f;
    DiffPair_float_0 _S1003;
    (&_S1003)->primal_0 = _S536;
    (&_S1003)->differential_0 = 0.0f;
    _d_min_0(&_S1002, &_S1003, _S995.differential_0);
    float3  _S1004 = make_float3 (_S934.x, _S934.y, _S1001.differential_0 + _S1003.differential_0 + _S933.x + _S933.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1005;
    (&_S1005)->primal_0 = R_4;
    (&_S1005)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1006;
    (&_S1006)->primal_0 = pos_i_1;
    (&_S1006)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S1005, &_S1006, _S1004);
    DiffPair_float_0 _S1007;
    (&_S1007)->primal_0 = _S529;
    (&_S1007)->differential_0 = 0.0f;
    DiffPair_float_0 _S1008;
    (&_S1008)->primal_0 = _S531;
    (&_S1008)->differential_0 = 0.0f;
    _d_max_0(&_S1007, &_S1008, _S1000.differential_0);
    DiffPair_float_0 _S1009;
    (&_S1009)->primal_0 = _S528;
    (&_S1009)->differential_0 = 0.0f;
    DiffPair_float_0 _S1010;
    (&_S1010)->primal_0 = _S531;
    (&_S1010)->differential_0 = 0.0f;
    _d_min_0(&_S1009, &_S1010, _S1002.differential_0);
    float3  _S1011 = make_float3 (_S948.x, _S948.y, _S1008.differential_0 + _S1010.differential_0 + _S947.x + _S947.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1012;
    (&_S1012)->primal_0 = R_4;
    (&_S1012)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1013;
    (&_S1013)->primal_0 = pos_i_0;
    (&_S1013)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S1012, &_S1013, _S1011);
    DiffPair_float_0 _S1014;
    (&_S1014)->primal_0 = 0.0f;
    (&_S1014)->differential_0 = 0.0f;
    DiffPair_float_0 _S1015;
    (&_S1015)->primal_0 = _S526;
    (&_S1015)->differential_0 = 0.0f;
    _d_max_0(&_S1014, &_S1015, _S1007.differential_0);
    DiffPair_float_0 _S1016;
    (&_S1016)->primal_0 = 1.00000001504746622e+30f;
    (&_S1016)->differential_0 = 0.0f;
    DiffPair_float_0 _S1017;
    (&_S1017)->primal_0 = _S526;
    (&_S1017)->differential_0 = 0.0f;
    _d_min_0(&_S1016, &_S1017, _S1009.differential_0);
    float3  _S1018 = make_float3 (_S962.x, _S962.y, _S1015.differential_0 + _S1017.differential_0 + _S961.x + _S961.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1019;
    (&_S1019)->primal_0 = R_4;
    (&_S1019)->differential_0 = _S791;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1020;
    (&_S1020)->primal_0 = pos_4;
    (&_S1020)->differential_0 = _S729;
    s_bwd_prop_mul_0(&_S1019, &_S1020, _S1018);
    float3  _S1021 = _S793.differential_0 + _S789 + _S969 + _S976 + _S983 + _S990 + _S997 + _S1004 + _S1011 + _S1018;
    Matrix<float, 3, 3>  _S1022 = _S794 + _S963.differential_0 + _S970.differential_0 + _S977.differential_0 + _S984.differential_0 + _S991.differential_0 + _S998.differential_0 + _S1005.differential_0 + _S1012.differential_0 + _S1019.differential_0;
    FixedArray<float3 , 16>  _S1023;
    _S1023[int(0)] = _S729;
    _S1023[int(1)] = _S729;
    _S1023[int(2)] = _S729;
    _S1023[int(3)] = _S729;
    _S1023[int(4)] = _S729;
    _S1023[int(5)] = _S729;
    _S1023[int(6)] = _S729;
    _S1023[int(7)] = _S729;
    _S1023[int(8)] = _S729;
    _S1023[int(9)] = _S729;
    _S1023[int(10)] = _S729;
    _S1023[int(11)] = _S729;
    _S1023[int(12)] = _S729;
    _S1023[int(13)] = _S729;
    _S1023[int(14)] = _S729;
    _S1023[int(15)] = _S729;
    _S1023[int(15)] = _S732;
    _S1023[int(14)] = _S734;
    _S1023[int(13)] = _S736;
    _S1023[int(12)] = _S738;
    _S1023[int(11)] = _S740;
    _S1023[int(10)] = _S742;
    _S1023[int(9)] = _S744;
    _S1023[int(8)] = _S752;
    _S1023[int(7)] = _S754;
    _S1023[int(6)] = _S756;
    _S1023[int(5)] = _S758;
    _S1023[int(4)] = _S760;
    _S1023[int(3)] = _S771;
    _S1023[int(2)] = _S773;
    _S1023[int(1)] = _S775;
    _S1023[int(0)] = _S788;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_0 = _S1023;
    *v_R_0 = _S1022;
    *v_t_0 = _S1021;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S1024, float _S1025)
{
    return (F32_atan2((_S1024), (_S1025)));
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float _s_dOut_0)
{
    float _S1026 = (*dpx_6).primal_0.x;
    float _S1027 = (*dpx_6).primal_0.y;
    float _S1028 = (*dpx_6).primal_0.z;
    DiffPair_float_0 _S1029;
    (&_S1029)->primal_0 = _S1026 * _S1026 + _S1027 * _S1027 + _S1028 * _S1028;
    (&_S1029)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1029, _s_dOut_0);
    float _S1030 = (*dpx_6).primal_0.z * _S1029.differential_0;
    float _S1031 = _S1030 + _S1030;
    float _S1032 = (*dpx_6).primal_0.y * _S1029.differential_0;
    float _S1033 = _S1032 + _S1032;
    float _S1034 = (*dpx_6).primal_0.x * _S1029.differential_0;
    float _S1035 = _S1034 + _S1034;
    float3  _S1036 = make_float3 (0.0f);
    *&((&_S1036)->z) = _S1031;
    *&((&_S1036)->y) = _S1033;
    *&((&_S1036)->x) = _S1035;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S1036;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1037, float _S1038)
{
    s_bwd_prop_length_impl_0(_S1037, _S1038);
    return;
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S1039, DiffPair_float_0 * _S1040, float _S1041)
{
    _d_atan2_0(_S1039, _S1040, _S1041);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_7, float _s_dOut_1)
{
    float _S1042 = (*dpx_7).primal_0.x;
    float _S1043 = (*dpx_7).primal_0.y;
    DiffPair_float_0 _S1044;
    (&_S1044)->primal_0 = _S1042 * _S1042 + _S1043 * _S1043;
    (&_S1044)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1044, _s_dOut_1);
    float _S1045 = (*dpx_7).primal_0.y * _S1044.differential_0;
    float _S1046 = _S1045 + _S1045;
    float _S1047 = (*dpx_7).primal_0.x * _S1044.differential_0;
    float _S1048 = _S1047 + _S1047;
    float2  _S1049 = make_float2 (0.0f);
    *&((&_S1049)->y) = _S1046;
    *&((&_S1049)->x) = _S1048;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S1049;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1050, float _S1051)
{
    s_bwd_prop_length_impl_1(_S1050, _S1051);
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  densities_5, FixedArray<float3 , 16>  sh_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float3  v_rgb_1, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  _S1052 = s_primal_ctx_mul_0(R_5, pos_5) + t_5;
    float _S1053 = length_1(_S1052);
    float _S1054 = (F32_min((1.00000001504746622e+30f), (_S1053)));
    float _S1055 = (F32_max((0.0f), (_S1053)));
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S1056 = s_primal_ctx_mul_0(R_5, pos_i_7) + t_5;
    float _S1057 = length_1(_S1056);
    float _S1058 = (F32_min((_S1054), (_S1057)));
    float _S1059 = (F32_max((_S1055), (_S1057)));
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S1060 = s_primal_ctx_mul_0(R_5, pos_i_8) + t_5;
    float _S1061 = length_1(_S1060);
    float _S1062 = (F32_min((_S1058), (_S1061)));
    float _S1063 = (F32_max((_S1059), (_S1061)));
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S1064 = s_primal_ctx_mul_0(R_5, pos_i_9) + t_5;
    float _S1065 = length_1(_S1064);
    float _S1066 = (F32_min((_S1062), (_S1065)));
    float _S1067 = (F32_max((_S1063), (_S1065)));
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S1068 = s_primal_ctx_mul_0(R_5, pos_i_10) + t_5;
    float _S1069 = length_1(_S1068);
    float _S1070 = (F32_min((_S1066), (_S1069)));
    float _S1071 = (F32_max((_S1067), (_S1069)));
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S1072 = s_primal_ctx_mul_0(R_5, pos_i_11) + t_5;
    float _S1073 = length_1(_S1072);
    float _S1074 = (F32_min((_S1070), (_S1073)));
    float _S1075 = (F32_max((_S1071), (_S1073)));
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S1076 = s_primal_ctx_mul_0(R_5, pos_i_12) + t_5;
    float _S1077 = length_1(_S1076);
    float _S1078 = (F32_min((_S1074), (_S1077)));
    float _S1079 = (F32_max((_S1075), (_S1077)));
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S1080 = s_primal_ctx_mul_0(R_5, pos_i_13) + t_5;
    float _S1081 = length_1(_S1080);
    float3  _S1082 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_4 = s_primal_ctx_mul_0(R_5, _S1082) + t_5;
    float2  _S1083 = float2 {_S1052.x, _S1052.y};
    float _S1084 = length_0(_S1083);
    float _S1085 = _S1052.z;
    float _S1086 = s_primal_ctx_atan2_0(_S1084, _S1085);
    bool _S1087 = _S1086 < 0.00100000004749745f;
    float k_2;
    float _S1088;
    float _S1089;
    float _S1090;
    if(_S1087)
    {
        float _S1091 = 1.0f - _S1086 * _S1086 / 3.0f;
        float _S1092 = _S1085 * _S1085;
        k_2 = _S1091 / _S1085;
        _S1088 = _S1092;
        _S1089 = _S1091;
        _S1090 = 0.0f;
    }
    else
    {
        float _S1093 = _S1084 * _S1084;
        k_2 = _S1086 / _S1084;
        _S1088 = 0.0f;
        _S1089 = 0.0f;
        _S1090 = _S1093;
    }
    float2  _S1094 = make_float2 (k_2);
    float2  _S1095 = _S1083 * make_float2 (k_2);
    float u_56 = _S1095.x;
    float v_56 = _S1095.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S1096 = dist_coeffs_5[int(2)] + r2_56 * dist_coeffs_5[int(3)];
    float _S1097 = dist_coeffs_5[int(1)] + r2_56 * _S1096;
    float _S1098 = dist_coeffs_5[int(0)] + r2_56 * _S1097;
    float radial_24 = 1.0f + r2_56 * _S1098;
    float _S1099 = 2.0f * dist_coeffs_5[int(4)];
    float _S1100 = _S1099 * u_56;
    float _S1101 = 2.0f * u_56;
    float _S1102 = 2.0f * dist_coeffs_5[int(5)];
    float _S1103 = _S1102 * u_56;
    float _S1104 = 2.0f * v_56;
    float2  _S1105 = _S1095 * make_float2 (radial_24) + make_float2 (_S1100 * v_56 + dist_coeffs_5[int(5)] * (r2_56 + _S1101 * u_56) + dist_coeffs_5[int(6)] * r2_56, _S1103 * v_56 + dist_coeffs_5[int(4)] * (r2_56 + _S1104 * v_56) + dist_coeffs_5[int(7)] * r2_56);
    float2  _S1106 = _S1105 + make_float2 (dist_coeffs_5[int(8)] * _S1105.x + dist_coeffs_5[int(9)] * _S1105.y, 0.0f);
    float _S1107 = fx_5 * _S1106.x + cx_5;
    float _S1108 = fy_5 * _S1106.y + cy_5;
    float2  _S1109 = float2 {_S1056.x, _S1056.y};
    float _S1110 = length_0(_S1109);
    float _S1111 = _S1056.z;
    float _S1112 = s_primal_ctx_atan2_0(_S1110, _S1111);
    bool _S1113 = _S1112 < 0.00100000004749745f;
    float _S1114;
    float _S1115;
    float _S1116;
    if(_S1113)
    {
        float _S1117 = 1.0f - _S1112 * _S1112 / 3.0f;
        float _S1118 = _S1111 * _S1111;
        k_2 = _S1117 / _S1111;
        _S1114 = _S1118;
        _S1115 = _S1117;
        _S1116 = 0.0f;
    }
    else
    {
        float _S1119 = _S1110 * _S1110;
        k_2 = _S1112 / _S1110;
        _S1114 = 0.0f;
        _S1115 = 0.0f;
        _S1116 = _S1119;
    }
    float2  _S1120 = make_float2 (k_2);
    float2  _S1121 = _S1109 * make_float2 (k_2);
    float u_57 = _S1121.x;
    float v_57 = _S1121.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S1122 = dist_coeffs_5[int(2)] + r2_57 * dist_coeffs_5[int(3)];
    float _S1123 = dist_coeffs_5[int(1)] + r2_57 * _S1122;
    float _S1124 = dist_coeffs_5[int(0)] + r2_57 * _S1123;
    float radial_25 = 1.0f + r2_57 * _S1124;
    float _S1125 = _S1099 * u_57;
    float _S1126 = 2.0f * u_57;
    float _S1127 = _S1102 * u_57;
    float _S1128 = 2.0f * v_57;
    float2  _S1129 = _S1121 * make_float2 (radial_25) + make_float2 (_S1125 * v_57 + dist_coeffs_5[int(5)] * (r2_57 + _S1126 * u_57) + dist_coeffs_5[int(6)] * r2_57, _S1127 * v_57 + dist_coeffs_5[int(4)] * (r2_57 + _S1128 * v_57) + dist_coeffs_5[int(7)] * r2_57);
    float2  _S1130 = _S1129 + make_float2 (dist_coeffs_5[int(8)] * _S1129.x + dist_coeffs_5[int(9)] * _S1129.y, 0.0f);
    float _S1131 = fx_5 * _S1130.x + cx_5;
    float _S1132 = fy_5 * _S1130.y + cy_5;
    float2  _S1133 = float2 {_S1060.x, _S1060.y};
    float _S1134 = length_0(_S1133);
    float _S1135 = _S1060.z;
    float _S1136 = s_primal_ctx_atan2_0(_S1134, _S1135);
    bool _S1137 = _S1136 < 0.00100000004749745f;
    float _S1138;
    float _S1139;
    float _S1140;
    if(_S1137)
    {
        float _S1141 = 1.0f - _S1136 * _S1136 / 3.0f;
        float _S1142 = _S1135 * _S1135;
        k_2 = _S1141 / _S1135;
        _S1138 = _S1142;
        _S1139 = _S1141;
        _S1140 = 0.0f;
    }
    else
    {
        float _S1143 = _S1134 * _S1134;
        k_2 = _S1136 / _S1134;
        _S1138 = 0.0f;
        _S1139 = 0.0f;
        _S1140 = _S1143;
    }
    float2  _S1144 = make_float2 (k_2);
    float2  _S1145 = _S1133 * make_float2 (k_2);
    float u_58 = _S1145.x;
    float v_58 = _S1145.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S1146 = dist_coeffs_5[int(2)] + r2_58 * dist_coeffs_5[int(3)];
    float _S1147 = dist_coeffs_5[int(1)] + r2_58 * _S1146;
    float _S1148 = dist_coeffs_5[int(0)] + r2_58 * _S1147;
    float radial_26 = 1.0f + r2_58 * _S1148;
    float _S1149 = _S1099 * u_58;
    float _S1150 = 2.0f * u_58;
    float _S1151 = _S1102 * u_58;
    float _S1152 = 2.0f * v_58;
    float2  _S1153 = _S1145 * make_float2 (radial_26) + make_float2 (_S1149 * v_58 + dist_coeffs_5[int(5)] * (r2_58 + _S1150 * u_58) + dist_coeffs_5[int(6)] * r2_58, _S1151 * v_58 + dist_coeffs_5[int(4)] * (r2_58 + _S1152 * v_58) + dist_coeffs_5[int(7)] * r2_58);
    float2  _S1154 = _S1153 + make_float2 (dist_coeffs_5[int(8)] * _S1153.x + dist_coeffs_5[int(9)] * _S1153.y, 0.0f);
    float _S1155 = fx_5 * _S1154.x + cx_5;
    float _S1156 = fy_5 * _S1154.y + cy_5;
    float2  _S1157 = float2 {_S1064.x, _S1064.y};
    float _S1158 = length_0(_S1157);
    float _S1159 = _S1064.z;
    float _S1160 = s_primal_ctx_atan2_0(_S1158, _S1159);
    bool _S1161 = _S1160 < 0.00100000004749745f;
    float _S1162;
    float _S1163;
    float _S1164;
    if(_S1161)
    {
        float _S1165 = 1.0f - _S1160 * _S1160 / 3.0f;
        float _S1166 = _S1159 * _S1159;
        k_2 = _S1165 / _S1159;
        _S1162 = _S1166;
        _S1163 = _S1165;
        _S1164 = 0.0f;
    }
    else
    {
        float _S1167 = _S1158 * _S1158;
        k_2 = _S1160 / _S1158;
        _S1162 = 0.0f;
        _S1163 = 0.0f;
        _S1164 = _S1167;
    }
    float2  _S1168 = make_float2 (k_2);
    float2  _S1169 = _S1157 * make_float2 (k_2);
    float u_59 = _S1169.x;
    float v_59 = _S1169.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S1170 = dist_coeffs_5[int(2)] + r2_59 * dist_coeffs_5[int(3)];
    float _S1171 = dist_coeffs_5[int(1)] + r2_59 * _S1170;
    float _S1172 = dist_coeffs_5[int(0)] + r2_59 * _S1171;
    float radial_27 = 1.0f + r2_59 * _S1172;
    float _S1173 = _S1099 * u_59;
    float _S1174 = 2.0f * u_59;
    float _S1175 = _S1102 * u_59;
    float _S1176 = 2.0f * v_59;
    float2  _S1177 = _S1169 * make_float2 (radial_27) + make_float2 (_S1173 * v_59 + dist_coeffs_5[int(5)] * (r2_59 + _S1174 * u_59) + dist_coeffs_5[int(6)] * r2_59, _S1175 * v_59 + dist_coeffs_5[int(4)] * (r2_59 + _S1176 * v_59) + dist_coeffs_5[int(7)] * r2_59);
    float2  _S1178 = _S1177 + make_float2 (dist_coeffs_5[int(8)] * _S1177.x + dist_coeffs_5[int(9)] * _S1177.y, 0.0f);
    float _S1179 = fx_5 * _S1178.x + cx_5;
    float _S1180 = fy_5 * _S1178.y + cy_5;
    float2  _S1181 = float2 {_S1068.x, _S1068.y};
    float _S1182 = length_0(_S1181);
    float _S1183 = _S1068.z;
    float _S1184 = s_primal_ctx_atan2_0(_S1182, _S1183);
    bool _S1185 = _S1184 < 0.00100000004749745f;
    float _S1186;
    float _S1187;
    float _S1188;
    if(_S1185)
    {
        float _S1189 = 1.0f - _S1184 * _S1184 / 3.0f;
        float _S1190 = _S1183 * _S1183;
        k_2 = _S1189 / _S1183;
        _S1186 = _S1190;
        _S1187 = _S1189;
        _S1188 = 0.0f;
    }
    else
    {
        float _S1191 = _S1182 * _S1182;
        k_2 = _S1184 / _S1182;
        _S1186 = 0.0f;
        _S1187 = 0.0f;
        _S1188 = _S1191;
    }
    float2  _S1192 = make_float2 (k_2);
    float2  _S1193 = _S1181 * make_float2 (k_2);
    float u_60 = _S1193.x;
    float v_60 = _S1193.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S1194 = dist_coeffs_5[int(2)] + r2_60 * dist_coeffs_5[int(3)];
    float _S1195 = dist_coeffs_5[int(1)] + r2_60 * _S1194;
    float _S1196 = dist_coeffs_5[int(0)] + r2_60 * _S1195;
    float radial_28 = 1.0f + r2_60 * _S1196;
    float _S1197 = _S1099 * u_60;
    float _S1198 = 2.0f * u_60;
    float _S1199 = _S1102 * u_60;
    float _S1200 = 2.0f * v_60;
    float2  _S1201 = _S1193 * make_float2 (radial_28) + make_float2 (_S1197 * v_60 + dist_coeffs_5[int(5)] * (r2_60 + _S1198 * u_60) + dist_coeffs_5[int(6)] * r2_60, _S1199 * v_60 + dist_coeffs_5[int(4)] * (r2_60 + _S1200 * v_60) + dist_coeffs_5[int(7)] * r2_60);
    float2  _S1202 = _S1201 + make_float2 (dist_coeffs_5[int(8)] * _S1201.x + dist_coeffs_5[int(9)] * _S1201.y, 0.0f);
    float _S1203 = fx_5 * _S1202.x + cx_5;
    float _S1204 = fy_5 * _S1202.y + cy_5;
    float2  _S1205 = float2 {_S1072.x, _S1072.y};
    float _S1206 = length_0(_S1205);
    float _S1207 = _S1072.z;
    float _S1208 = s_primal_ctx_atan2_0(_S1206, _S1207);
    bool _S1209 = _S1208 < 0.00100000004749745f;
    float _S1210;
    float _S1211;
    float _S1212;
    if(_S1209)
    {
        float _S1213 = 1.0f - _S1208 * _S1208 / 3.0f;
        float _S1214 = _S1207 * _S1207;
        k_2 = _S1213 / _S1207;
        _S1210 = _S1214;
        _S1211 = _S1213;
        _S1212 = 0.0f;
    }
    else
    {
        float _S1215 = _S1206 * _S1206;
        k_2 = _S1208 / _S1206;
        _S1210 = 0.0f;
        _S1211 = 0.0f;
        _S1212 = _S1215;
    }
    float2  _S1216 = make_float2 (k_2);
    float2  _S1217 = _S1205 * make_float2 (k_2);
    float u_61 = _S1217.x;
    float v_61 = _S1217.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S1218 = dist_coeffs_5[int(2)] + r2_61 * dist_coeffs_5[int(3)];
    float _S1219 = dist_coeffs_5[int(1)] + r2_61 * _S1218;
    float _S1220 = dist_coeffs_5[int(0)] + r2_61 * _S1219;
    float radial_29 = 1.0f + r2_61 * _S1220;
    float _S1221 = _S1099 * u_61;
    float _S1222 = 2.0f * u_61;
    float _S1223 = _S1102 * u_61;
    float _S1224 = 2.0f * v_61;
    float2  _S1225 = _S1217 * make_float2 (radial_29) + make_float2 (_S1221 * v_61 + dist_coeffs_5[int(5)] * (r2_61 + _S1222 * u_61) + dist_coeffs_5[int(6)] * r2_61, _S1223 * v_61 + dist_coeffs_5[int(4)] * (r2_61 + _S1224 * v_61) + dist_coeffs_5[int(7)] * r2_61);
    float2  _S1226 = _S1225 + make_float2 (dist_coeffs_5[int(8)] * _S1225.x + dist_coeffs_5[int(9)] * _S1225.y, 0.0f);
    float _S1227 = fx_5 * _S1226.x + cx_5;
    float _S1228 = fy_5 * _S1226.y + cy_5;
    float2  _S1229 = float2 {_S1076.x, _S1076.y};
    float _S1230 = length_0(_S1229);
    float _S1231 = _S1076.z;
    float _S1232 = s_primal_ctx_atan2_0(_S1230, _S1231);
    bool _S1233 = _S1232 < 0.00100000004749745f;
    float _S1234;
    float _S1235;
    float _S1236;
    if(_S1233)
    {
        float _S1237 = 1.0f - _S1232 * _S1232 / 3.0f;
        float _S1238 = _S1231 * _S1231;
        k_2 = _S1237 / _S1231;
        _S1234 = _S1238;
        _S1235 = _S1237;
        _S1236 = 0.0f;
    }
    else
    {
        float _S1239 = _S1230 * _S1230;
        k_2 = _S1232 / _S1230;
        _S1234 = 0.0f;
        _S1235 = 0.0f;
        _S1236 = _S1239;
    }
    float2  _S1240 = make_float2 (k_2);
    float2  _S1241 = _S1229 * make_float2 (k_2);
    float u_62 = _S1241.x;
    float v_62 = _S1241.y;
    float r2_62 = u_62 * u_62 + v_62 * v_62;
    float _S1242 = dist_coeffs_5[int(2)] + r2_62 * dist_coeffs_5[int(3)];
    float _S1243 = dist_coeffs_5[int(1)] + r2_62 * _S1242;
    float _S1244 = dist_coeffs_5[int(0)] + r2_62 * _S1243;
    float radial_30 = 1.0f + r2_62 * _S1244;
    float _S1245 = _S1099 * u_62;
    float _S1246 = 2.0f * u_62;
    float _S1247 = _S1102 * u_62;
    float _S1248 = 2.0f * v_62;
    float2  _S1249 = _S1241 * make_float2 (radial_30) + make_float2 (_S1245 * v_62 + dist_coeffs_5[int(5)] * (r2_62 + _S1246 * u_62) + dist_coeffs_5[int(6)] * r2_62, _S1247 * v_62 + dist_coeffs_5[int(4)] * (r2_62 + _S1248 * v_62) + dist_coeffs_5[int(7)] * r2_62);
    float2  _S1250 = _S1249 + make_float2 (dist_coeffs_5[int(8)] * _S1249.x + dist_coeffs_5[int(9)] * _S1249.y, 0.0f);
    float _S1251 = fx_5 * _S1250.x + cx_5;
    float _S1252 = fy_5 * _S1250.y + cy_5;
    float2  _S1253 = float2 {_S1080.x, _S1080.y};
    float _S1254 = length_0(_S1253);
    float _S1255 = _S1080.z;
    float _S1256 = s_primal_ctx_atan2_0(_S1254, _S1255);
    bool _S1257 = _S1256 < 0.00100000004749745f;
    float _S1258;
    float _S1259;
    float _S1260;
    if(_S1257)
    {
        float _S1261 = 1.0f - _S1256 * _S1256 / 3.0f;
        float _S1262 = _S1255 * _S1255;
        k_2 = _S1261 / _S1255;
        _S1258 = _S1262;
        _S1259 = _S1261;
        _S1260 = 0.0f;
    }
    else
    {
        float _S1263 = _S1254 * _S1254;
        k_2 = _S1256 / _S1254;
        _S1258 = 0.0f;
        _S1259 = 0.0f;
        _S1260 = _S1263;
    }
    float2  _S1264 = make_float2 (k_2);
    float2  _S1265 = _S1253 * make_float2 (k_2);
    float u_63 = _S1265.x;
    float v_63 = _S1265.y;
    float r2_63 = u_63 * u_63 + v_63 * v_63;
    float _S1266 = dist_coeffs_5[int(2)] + r2_63 * dist_coeffs_5[int(3)];
    float _S1267 = dist_coeffs_5[int(1)] + r2_63 * _S1266;
    float _S1268 = dist_coeffs_5[int(0)] + r2_63 * _S1267;
    float radial_31 = 1.0f + r2_63 * _S1268;
    float _S1269 = _S1099 * u_63;
    float _S1270 = 2.0f * u_63;
    float _S1271 = _S1102 * u_63;
    float _S1272 = 2.0f * v_63;
    float2  _S1273 = _S1265 * make_float2 (radial_31) + make_float2 (_S1269 * v_63 + dist_coeffs_5[int(5)] * (r2_63 + _S1270 * u_63) + dist_coeffs_5[int(6)] * r2_63, _S1271 * v_63 + dist_coeffs_5[int(4)] * (r2_63 + _S1272 * v_63) + dist_coeffs_5[int(7)] * r2_63);
    float2  _S1274 = _S1273 + make_float2 (dist_coeffs_5[int(8)] * _S1273.x + dist_coeffs_5[int(9)] * _S1273.y, 0.0f);
    float _S1275 = fx_5 * _S1274.x + cx_5;
    float _S1276 = fy_5 * _S1274.y + cy_5;
    float _S1277 = (F32_max((_S1107), (_S1131)));
    float _S1278 = (F32_min((_S1107), (_S1131)));
    float _S1279 = (F32_max((_S1108), (_S1132)));
    float _S1280 = (F32_min((_S1108), (_S1132)));
    float _S1281 = (F32_max((_S1277), (_S1155)));
    float _S1282 = (F32_min((_S1278), (_S1155)));
    float _S1283 = (F32_max((_S1279), (_S1156)));
    float _S1284 = (F32_min((_S1280), (_S1156)));
    float _S1285 = (F32_max((_S1281), (_S1179)));
    float _S1286 = (F32_min((_S1282), (_S1179)));
    float _S1287 = (F32_max((_S1283), (_S1180)));
    float _S1288 = (F32_min((_S1284), (_S1180)));
    float _S1289 = (F32_max((_S1285), (_S1203)));
    float _S1290 = (F32_min((_S1286), (_S1203)));
    float _S1291 = (F32_max((_S1287), (_S1204)));
    float _S1292 = (F32_min((_S1288), (_S1204)));
    float _S1293 = (F32_max((_S1289), (_S1227)));
    float _S1294 = (F32_min((_S1290), (_S1227)));
    float _S1295 = (F32_max((_S1291), (_S1228)));
    float _S1296 = (F32_min((_S1292), (_S1228)));
    float _S1297 = (F32_max((_S1293), (_S1251)));
    float _S1298 = (F32_min((_S1294), (_S1251)));
    float _S1299 = (F32_max((_S1295), (_S1252)));
    float _S1300 = (F32_min((_S1296), (_S1252)));
    float _S1301 = mean_c_4.z;
    float _S1302 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S1303 = _S1302 * (_S1301 + length_1(mean_c_4));
    Matrix<float, 3, 3>  _S1304 = transpose_1(R_5);
    float3  _S1305 = mean_c_4 - - s_primal_ctx_mul_0(_S1304, t_5);
    float _S1306 = _S1305.x;
    float _S1307 = _S1305.y;
    float _S1308 = _S1305.z;
    float _S1309 = _S1306 * _S1306 + _S1307 * _S1307 + _S1308 * _S1308;
    float _S1310 = s_primal_ctx_sqrt_0(_S1309);
    float x_12 = _S1306 / _S1310;
    float3  _S1311 = make_float3 (x_12);
    float _S1312 = _S1310 * _S1310;
    float y_8 = _S1307 / _S1310;
    float z_5 = _S1308 / _S1310;
    float3  _S1313 = make_float3 (z_5);
    float _S1314 = - y_8;
    float3  _S1315 = make_float3 (_S1314);
    float z2_5 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_12 * x_12 - y_8 * y_8;
    float _S1316 = 2.0f * x_12;
    float fS1_5 = _S1316 * y_8;
    float pSH6_1 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float3  _S1317 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_12;
    float3  _S1318 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_8;
    float3  _S1319 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S1320 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S1321 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    float _S1322 = 1.86588168144226074f * z2_5 - 1.11952900886535645f;
    float pSH12_1 = z_5 * _S1322;
    float3  _S1323 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_12;
    float3  _S1324 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_8;
    float3  _S1325 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S1326 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S1327 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_12 * fC1_5 - y_8 * fS1_5);
    float3  _S1328 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_12 * fS1_5 + y_8 * fC1_5);
    float3  _S1329 = make_float3 (pSH9_1);
    float3  _S1330 = make_float3 (0.0f);
    float3  _S1331 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1332;
    (&_S1332)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1314) * sh_coeffs_5[int(1)] + make_float3 (z_5) * sh_coeffs_5[int(2)] - make_float3 (x_12) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]) + make_float3 (0.5f);
    (&_S1332)->differential_0 = _S1331;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1333;
    (&_S1333)->primal_0 = _S1330;
    (&_S1333)->differential_0 = _S1331;
    s_bwd_prop_max_0(&_S1332, &_S1333, v_rgb_1);
    float3  _S1334 = _S1328 * _S1332.differential_0;
    float3  _S1335 = sh_coeffs_5[int(15)] * _S1332.differential_0;
    float3  _S1336 = _S1326 * _S1332.differential_0;
    float3  _S1337 = sh_coeffs_5[int(14)] * _S1332.differential_0;
    float3  _S1338 = _S1324 * _S1332.differential_0;
    float3  _S1339 = sh_coeffs_5[int(13)] * _S1332.differential_0;
    float3  _S1340 = _S1323 * _S1332.differential_0;
    float3  _S1341 = sh_coeffs_5[int(12)] * _S1332.differential_0;
    float3  _S1342 = _S1325 * _S1332.differential_0;
    float3  _S1343 = sh_coeffs_5[int(11)] * _S1332.differential_0;
    float3  _S1344 = _S1327 * _S1332.differential_0;
    float3  _S1345 = sh_coeffs_5[int(10)] * _S1332.differential_0;
    float3  _S1346 = _S1329 * _S1332.differential_0;
    float3  _S1347 = sh_coeffs_5[int(9)] * _S1332.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1347.x + _S1347.y + _S1347.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1335.x + _S1335.y + _S1335.z);
    float _S1348 = _S1345.x + _S1345.y + _S1345.z;
    float _S1349 = _S1337.x + _S1337.y + _S1337.z;
    float _S1350 = _S1343.x + _S1343.y + _S1343.z;
    float _S1351 = _S1339.x + _S1339.y + _S1339.z;
    float _S1352 = _S1341.x + _S1341.y + _S1341.z;
    float _S1353 = - s_diff_fC2_T_1;
    float3  _S1354 = _S1320 * _S1332.differential_0;
    float3  _S1355 = sh_coeffs_5[int(8)] * _S1332.differential_0;
    float3  _S1356 = _S1318 * _S1332.differential_0;
    float3  _S1357 = sh_coeffs_5[int(7)] * _S1332.differential_0;
    float3  _S1358 = _S1317 * _S1332.differential_0;
    float3  _S1359 = sh_coeffs_5[int(6)] * _S1332.differential_0;
    float3  _S1360 = _S1319 * _S1332.differential_0;
    float3  _S1361 = sh_coeffs_5[int(5)] * _S1332.differential_0;
    float3  _S1362 = _S1321 * _S1332.differential_0;
    float3  _S1363 = sh_coeffs_5[int(4)] * _S1332.differential_0;
    float _S1364 = _S1361.x + _S1361.y + _S1361.z;
    float _S1365 = _S1357.x + _S1357.y + _S1357.z;
    float _S1366 = fTmp1B_5 * _S1348 + x_12 * s_diff_fS2_T_1 + y_8 * _S1353 + 0.54627424478530884f * (_S1363.x + _S1363.y + _S1363.z);
    float _S1367 = fTmp1B_5 * _S1349 + y_8 * s_diff_fS2_T_1 + x_12 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1355.x + _S1355.y + _S1355.z);
    float _S1368 = y_8 * - _S1367;
    float _S1369 = x_12 * _S1367;
    float _S1370 = z_5 * (1.86588168144226074f * (z_5 * _S1352) + -2.28522896766662598f * (y_8 * _S1350 + x_12 * _S1351) + 0.94617468118667603f * (_S1359.x + _S1359.y + _S1359.z));
    float3  _S1371 = make_float3 (0.48860251903533936f) * _S1332.differential_0;
    float3  _S1372 = - _S1371;
    float3  _S1373 = _S1311 * _S1372;
    float3  _S1374 = sh_coeffs_5[int(3)] * _S1372;
    float3  _S1375 = _S1313 * _S1371;
    float3  _S1376 = sh_coeffs_5[int(2)] * _S1371;
    float3  _S1377 = _S1315 * _S1371;
    float3  _S1378 = sh_coeffs_5[int(1)] * _S1371;
    float _S1379 = (_S1322 * _S1352 + 1.44530570507049561f * (fS1_5 * _S1348 + fC1_5 * _S1349) + -1.09254848957061768f * (y_8 * _S1364 + x_12 * _S1365) + _S1370 + _S1370 + _S1376.x + _S1376.y + _S1376.z) / _S1312;
    float _S1380 = _S1310 * _S1379;
    float _S1381 = (fTmp0C_5 * _S1350 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S1353 + fTmp0B_5 * _S1364 + _S1316 * _S1366 + _S1368 + _S1368 + - (_S1378.x + _S1378.y + _S1378.z)) / _S1312;
    float _S1382 = _S1310 * _S1381;
    float _S1383 = (fTmp0C_5 * _S1351 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S1365 + 2.0f * (y_8 * _S1366) + _S1369 + _S1369 + _S1374.x + _S1374.y + _S1374.z) / _S1312;
    float _S1384 = _S1310 * _S1383;
    float _S1385 = _S1308 * - _S1379 + _S1307 * - _S1381 + _S1306 * - _S1383;
    DiffPair_float_0 _S1386;
    (&_S1386)->primal_0 = _S1309;
    (&_S1386)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1386, _S1385);
    float _S1387 = _S1308 * _S1386.differential_0;
    float _S1388 = _S1307 * _S1386.differential_0;
    float _S1389 = _S1306 * _S1386.differential_0;
    float3  _S1390 = make_float3 (0.282094806432724f) * _S1332.differential_0;
    float3  _S1391 = make_float3 (_S1384 + _S1389 + _S1389, _S1382 + _S1388 + _S1388, _S1380 + _S1387 + _S1387);
    float3  _S1392 = - - _S1391;
    Matrix<float, 3, 3>  _S1393 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1394;
    (&_S1394)->primal_0 = _S1304;
    (&_S1394)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1395;
    (&_S1395)->primal_0 = t_5;
    (&_S1395)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1394, &_S1395, _S1392);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1396 = _S1395;
    Matrix<float, 3, 3>  _S1397 = transpose_1(_S1394.differential_0);
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1301;
    (&_S1398)->differential_0 = 0.0f;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1303;
    (&_S1399)->differential_0 = 0.0f;
    _d_max_0(&_S1398, &_S1399, 0.0f);
    float _S1400 = _S1302 * _S1399.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1401;
    (&_S1401)->primal_0 = mean_c_4;
    (&_S1401)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1401, _S1400);
    float _S1402 = _S1398.differential_0 + _S1400;
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = _S1300;
    (&_S1403)->differential_0 = 0.0f;
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = _S1276;
    (&_S1404)->differential_0 = 0.0f;
    _d_min_0(&_S1403, &_S1404, 0.0f);
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = _S1299;
    (&_S1405)->differential_0 = 0.0f;
    DiffPair_float_0 _S1406;
    (&_S1406)->primal_0 = _S1276;
    (&_S1406)->differential_0 = 0.0f;
    _d_max_0(&_S1405, &_S1406, 0.0f);
    DiffPair_float_0 _S1407;
    (&_S1407)->primal_0 = _S1298;
    (&_S1407)->differential_0 = 0.0f;
    DiffPair_float_0 _S1408;
    (&_S1408)->primal_0 = _S1275;
    (&_S1408)->differential_0 = 0.0f;
    _d_min_0(&_S1407, &_S1408, 0.0f);
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1297;
    (&_S1409)->differential_0 = 0.0f;
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = _S1275;
    (&_S1410)->differential_0 = 0.0f;
    _d_max_0(&_S1409, &_S1410, 0.0f);
    DiffPair_float_0 _S1411;
    (&_S1411)->primal_0 = _S1296;
    (&_S1411)->differential_0 = 0.0f;
    DiffPair_float_0 _S1412;
    (&_S1412)->primal_0 = _S1252;
    (&_S1412)->differential_0 = 0.0f;
    _d_min_0(&_S1411, &_S1412, _S1403.differential_0);
    DiffPair_float_0 _S1413;
    (&_S1413)->primal_0 = _S1295;
    (&_S1413)->differential_0 = 0.0f;
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = _S1252;
    (&_S1414)->differential_0 = 0.0f;
    _d_max_0(&_S1413, &_S1414, _S1405.differential_0);
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1294;
    (&_S1415)->differential_0 = 0.0f;
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1251;
    (&_S1416)->differential_0 = 0.0f;
    _d_min_0(&_S1415, &_S1416, _S1407.differential_0);
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1293;
    (&_S1417)->differential_0 = 0.0f;
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = _S1251;
    (&_S1418)->differential_0 = 0.0f;
    _d_max_0(&_S1417, &_S1418, _S1409.differential_0);
    DiffPair_float_0 _S1419;
    (&_S1419)->primal_0 = _S1292;
    (&_S1419)->differential_0 = 0.0f;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = _S1228;
    (&_S1420)->differential_0 = 0.0f;
    _d_min_0(&_S1419, &_S1420, _S1411.differential_0);
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1291;
    (&_S1421)->differential_0 = 0.0f;
    DiffPair_float_0 _S1422;
    (&_S1422)->primal_0 = _S1228;
    (&_S1422)->differential_0 = 0.0f;
    _d_max_0(&_S1421, &_S1422, _S1413.differential_0);
    DiffPair_float_0 _S1423;
    (&_S1423)->primal_0 = _S1290;
    (&_S1423)->differential_0 = 0.0f;
    DiffPair_float_0 _S1424;
    (&_S1424)->primal_0 = _S1227;
    (&_S1424)->differential_0 = 0.0f;
    _d_min_0(&_S1423, &_S1424, _S1415.differential_0);
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = _S1289;
    (&_S1425)->differential_0 = 0.0f;
    DiffPair_float_0 _S1426;
    (&_S1426)->primal_0 = _S1227;
    (&_S1426)->differential_0 = 0.0f;
    _d_max_0(&_S1425, &_S1426, _S1417.differential_0);
    DiffPair_float_0 _S1427;
    (&_S1427)->primal_0 = _S1288;
    (&_S1427)->differential_0 = 0.0f;
    DiffPair_float_0 _S1428;
    (&_S1428)->primal_0 = _S1204;
    (&_S1428)->differential_0 = 0.0f;
    _d_min_0(&_S1427, &_S1428, _S1419.differential_0);
    DiffPair_float_0 _S1429;
    (&_S1429)->primal_0 = _S1287;
    (&_S1429)->differential_0 = 0.0f;
    DiffPair_float_0 _S1430;
    (&_S1430)->primal_0 = _S1204;
    (&_S1430)->differential_0 = 0.0f;
    _d_max_0(&_S1429, &_S1430, _S1421.differential_0);
    DiffPair_float_0 _S1431;
    (&_S1431)->primal_0 = _S1286;
    (&_S1431)->differential_0 = 0.0f;
    DiffPair_float_0 _S1432;
    (&_S1432)->primal_0 = _S1203;
    (&_S1432)->differential_0 = 0.0f;
    _d_min_0(&_S1431, &_S1432, _S1423.differential_0);
    DiffPair_float_0 _S1433;
    (&_S1433)->primal_0 = _S1285;
    (&_S1433)->differential_0 = 0.0f;
    DiffPair_float_0 _S1434;
    (&_S1434)->primal_0 = _S1203;
    (&_S1434)->differential_0 = 0.0f;
    _d_max_0(&_S1433, &_S1434, _S1425.differential_0);
    DiffPair_float_0 _S1435;
    (&_S1435)->primal_0 = _S1284;
    (&_S1435)->differential_0 = 0.0f;
    DiffPair_float_0 _S1436;
    (&_S1436)->primal_0 = _S1180;
    (&_S1436)->differential_0 = 0.0f;
    _d_min_0(&_S1435, &_S1436, _S1427.differential_0);
    DiffPair_float_0 _S1437;
    (&_S1437)->primal_0 = _S1283;
    (&_S1437)->differential_0 = 0.0f;
    DiffPair_float_0 _S1438;
    (&_S1438)->primal_0 = _S1180;
    (&_S1438)->differential_0 = 0.0f;
    _d_max_0(&_S1437, &_S1438, _S1429.differential_0);
    DiffPair_float_0 _S1439;
    (&_S1439)->primal_0 = _S1282;
    (&_S1439)->differential_0 = 0.0f;
    DiffPair_float_0 _S1440;
    (&_S1440)->primal_0 = _S1179;
    (&_S1440)->differential_0 = 0.0f;
    _d_min_0(&_S1439, &_S1440, _S1431.differential_0);
    DiffPair_float_0 _S1441;
    (&_S1441)->primal_0 = _S1281;
    (&_S1441)->differential_0 = 0.0f;
    DiffPair_float_0 _S1442;
    (&_S1442)->primal_0 = _S1179;
    (&_S1442)->differential_0 = 0.0f;
    _d_max_0(&_S1441, &_S1442, _S1433.differential_0);
    DiffPair_float_0 _S1443;
    (&_S1443)->primal_0 = _S1280;
    (&_S1443)->differential_0 = 0.0f;
    DiffPair_float_0 _S1444;
    (&_S1444)->primal_0 = _S1156;
    (&_S1444)->differential_0 = 0.0f;
    _d_min_0(&_S1443, &_S1444, _S1435.differential_0);
    DiffPair_float_0 _S1445;
    (&_S1445)->primal_0 = _S1279;
    (&_S1445)->differential_0 = 0.0f;
    DiffPair_float_0 _S1446;
    (&_S1446)->primal_0 = _S1156;
    (&_S1446)->differential_0 = 0.0f;
    _d_max_0(&_S1445, &_S1446, _S1437.differential_0);
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1278;
    (&_S1447)->differential_0 = 0.0f;
    DiffPair_float_0 _S1448;
    (&_S1448)->primal_0 = _S1155;
    (&_S1448)->differential_0 = 0.0f;
    _d_min_0(&_S1447, &_S1448, _S1439.differential_0);
    DiffPair_float_0 _S1449;
    (&_S1449)->primal_0 = _S1277;
    (&_S1449)->differential_0 = 0.0f;
    DiffPair_float_0 _S1450;
    (&_S1450)->primal_0 = _S1155;
    (&_S1450)->differential_0 = 0.0f;
    _d_max_0(&_S1449, &_S1450, _S1441.differential_0);
    DiffPair_float_0 _S1451;
    (&_S1451)->primal_0 = _S1108;
    (&_S1451)->differential_0 = 0.0f;
    DiffPair_float_0 _S1452;
    (&_S1452)->primal_0 = _S1132;
    (&_S1452)->differential_0 = 0.0f;
    _d_min_0(&_S1451, &_S1452, _S1443.differential_0);
    DiffPair_float_0 _S1453;
    (&_S1453)->primal_0 = _S1108;
    (&_S1453)->differential_0 = 0.0f;
    DiffPair_float_0 _S1454;
    (&_S1454)->primal_0 = _S1132;
    (&_S1454)->differential_0 = 0.0f;
    _d_max_0(&_S1453, &_S1454, _S1445.differential_0);
    DiffPair_float_0 _S1455;
    (&_S1455)->primal_0 = _S1107;
    (&_S1455)->differential_0 = 0.0f;
    DiffPair_float_0 _S1456;
    (&_S1456)->primal_0 = _S1131;
    (&_S1456)->differential_0 = 0.0f;
    _d_min_0(&_S1455, &_S1456, _S1447.differential_0);
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1107;
    (&_S1457)->differential_0 = 0.0f;
    DiffPair_float_0 _S1458;
    (&_S1458)->primal_0 = _S1131;
    (&_S1458)->differential_0 = 0.0f;
    _d_max_0(&_S1457, &_S1458, _S1449.differential_0);
    float _S1459 = fx_5 * (_S1408.differential_0 + _S1410.differential_0);
    float2  _S1460 = make_float2 (_S1459, fy_5 * (_S1404.differential_0 + _S1406.differential_0)) + make_float2 (dist_coeffs_5[int(8)] * _S1459, dist_coeffs_5[int(9)] * _S1459);
    float2  _S1461 = _S1265 * _S1460;
    float _S1462 = dist_coeffs_5[int(4)] * _S1460.y;
    float _S1463 = dist_coeffs_5[int(5)] * _S1460.x;
    float _S1464 = _S1461.x + _S1461.y;
    float _S1465 = r2_63 * _S1464;
    float _S1466 = r2_63 * _S1465;
    float _S1467 = dist_coeffs_5[int(7)] * _S1460.y + _S1462 + dist_coeffs_5[int(6)] * _S1460.x + _S1463 + _S1268 * _S1464 + _S1267 * _S1465 + _S1266 * _S1466 + dist_coeffs_5[int(3)] * (r2_63 * _S1466);
    float _S1468 = v_63 * _S1467;
    float _S1469 = u_63 * _S1467;
    float2  _S1470 = make_float2 (radial_31) * _S1460 + make_float2 (_S1102 * (v_63 * _S1460.y) + _S1270 * _S1463 + 2.0f * (u_63 * _S1463) + _S1099 * (v_63 * _S1460.x) + _S1469 + _S1469, _S1272 * _S1462 + 2.0f * (v_63 * _S1462) + _S1271 * _S1460.y + _S1269 * _S1460.x + _S1468 + _S1468);
    FixedArray<float3 , 16>  _S1471;
    _S1471[int(0)] = _S1331;
    _S1471[int(1)] = _S1331;
    _S1471[int(2)] = _S1331;
    _S1471[int(3)] = _S1331;
    _S1471[int(4)] = _S1331;
    _S1471[int(5)] = _S1331;
    _S1471[int(6)] = _S1331;
    _S1471[int(7)] = _S1331;
    _S1471[int(8)] = _S1331;
    _S1471[int(9)] = _S1331;
    _S1471[int(10)] = _S1331;
    _S1471[int(11)] = _S1331;
    _S1471[int(12)] = _S1331;
    _S1471[int(13)] = _S1331;
    _S1471[int(14)] = _S1331;
    _S1471[int(15)] = _S1331;
    _S1471[int(7)] = _S1356;
    _S1471[int(0)] = _S1390;
    _S1471[int(1)] = _S1377;
    _S1471[int(2)] = _S1375;
    _S1471[int(3)] = _S1373;
    _S1471[int(4)] = _S1362;
    _S1471[int(5)] = _S1360;
    _S1471[int(6)] = _S1358;
    _S1471[int(15)] = _S1334;
    _S1471[int(8)] = _S1354;
    _S1471[int(9)] = _S1346;
    _S1471[int(10)] = _S1344;
    _S1471[int(11)] = _S1342;
    _S1471[int(12)] = _S1340;
    _S1471[int(13)] = _S1338;
    _S1471[int(14)] = _S1336;
    float3  _S1472 = _S1471[int(0)];
    float3  _S1473 = _S1471[int(1)];
    float3  _S1474 = _S1471[int(2)];
    float3  _S1475 = _S1471[int(3)];
    float3  _S1476 = _S1471[int(4)];
    float3  _S1477 = _S1471[int(5)];
    float3  _S1478 = _S1471[int(6)];
    float3  _S1479 = _S1471[int(7)];
    float3  _S1480 = _S1471[int(8)];
    float3  _S1481 = _S1471[int(9)];
    float3  _S1482 = _S1471[int(10)];
    float3  _S1483 = _S1471[int(11)];
    float3  _S1484 = _S1471[int(12)];
    float3  _S1485 = _S1471[int(13)];
    float3  _S1486 = _S1471[int(14)];
    float3  _S1487 = _S1471[int(15)];
    float3  _S1488 = _S1391 + _S1401.differential_0 + make_float3 (0.0f, 0.0f, _S1402);
    float _S1489 = _S1456.differential_0 + _S1458.differential_0;
    float _S1490 = _S1451.differential_0 + _S1453.differential_0;
    float _S1491 = _S1452.differential_0 + _S1454.differential_0;
    float _S1492 = _S1448.differential_0 + _S1450.differential_0;
    float _S1493 = _S1455.differential_0 + _S1457.differential_0;
    float _S1494 = _S1412.differential_0 + _S1414.differential_0;
    float _S1495 = _S1416.differential_0 + _S1418.differential_0;
    float _S1496 = _S1420.differential_0 + _S1422.differential_0;
    float _S1497 = _S1424.differential_0 + _S1426.differential_0;
    float _S1498 = _S1428.differential_0 + _S1430.differential_0;
    float _S1499 = _S1432.differential_0 + _S1434.differential_0;
    float _S1500 = _S1436.differential_0 + _S1438.differential_0;
    float _S1501 = _S1440.differential_0 + _S1442.differential_0;
    float _S1502 = _S1444.differential_0 + _S1446.differential_0;
    float2  _S1503 = _S1253 * _S1470;
    float2  _S1504 = _S1264 * _S1470;
    float _S1505 = _S1503.x + _S1503.y;
    if(_S1257)
    {
        float _S1506 = _S1505 / _S1258;
        float _S1507 = _S1259 * - _S1506;
        float _S1508 = _S1256 * (0.3333333432674408f * - (_S1255 * _S1506));
        k_2 = _S1508 + _S1508;
        _S1258 = _S1507;
        _S1259 = 0.0f;
    }
    else
    {
        float _S1509 = _S1505 / _S1260;
        float _S1510 = _S1256 * - _S1509;
        k_2 = _S1254 * _S1509;
        _S1258 = 0.0f;
        _S1259 = _S1510;
    }
    DiffPair_float_0 _S1511;
    (&_S1511)->primal_0 = _S1254;
    (&_S1511)->differential_0 = 0.0f;
    DiffPair_float_0 _S1512;
    (&_S1512)->primal_0 = _S1255;
    (&_S1512)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1511, &_S1512, k_2);
    float _S1513 = _S1512.differential_0 + _S1258;
    float _S1514 = _S1511.differential_0 + _S1259;
    float2  _S1515 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1516;
    (&_S1516)->primal_0 = _S1253;
    (&_S1516)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1516, _S1514);
    float2  _S1517 = _S1516.differential_0 + _S1504;
    float _S1518 = fx_5 * _S1495;
    float2  _S1519 = make_float2 (_S1518, fy_5 * _S1494) + make_float2 (dist_coeffs_5[int(8)] * _S1518, dist_coeffs_5[int(9)] * _S1518);
    float2  _S1520 = _S1241 * _S1519;
    float _S1521 = dist_coeffs_5[int(4)] * _S1519.y;
    float _S1522 = dist_coeffs_5[int(5)] * _S1519.x;
    float _S1523 = _S1520.x + _S1520.y;
    float _S1524 = r2_62 * _S1523;
    float _S1525 = r2_62 * _S1524;
    float _S1526 = dist_coeffs_5[int(7)] * _S1519.y + _S1521 + dist_coeffs_5[int(6)] * _S1519.x + _S1522 + _S1244 * _S1523 + _S1243 * _S1524 + _S1242 * _S1525 + dist_coeffs_5[int(3)] * (r2_62 * _S1525);
    float _S1527 = v_62 * _S1526;
    float _S1528 = u_62 * _S1526;
    float2  _S1529 = make_float2 (radial_30) * _S1519 + make_float2 (_S1102 * (v_62 * _S1519.y) + _S1246 * _S1522 + 2.0f * (u_62 * _S1522) + _S1099 * (v_62 * _S1519.x) + _S1528 + _S1528, _S1248 * _S1521 + 2.0f * (v_62 * _S1521) + _S1247 * _S1519.y + _S1245 * _S1519.x + _S1527 + _S1527);
    float3  _S1530 = make_float3 (_S1517.x, _S1517.y, _S1513);
    float2  _S1531 = _S1229 * _S1529;
    float2  _S1532 = _S1240 * _S1529;
    float _S1533 = _S1531.x + _S1531.y;
    if(_S1233)
    {
        float _S1534 = _S1533 / _S1234;
        float _S1535 = _S1235 * - _S1534;
        float _S1536 = _S1232 * (0.3333333432674408f * - (_S1231 * _S1534));
        k_2 = _S1536 + _S1536;
        _S1234 = _S1535;
        _S1235 = 0.0f;
    }
    else
    {
        float _S1537 = _S1533 / _S1236;
        float _S1538 = _S1232 * - _S1537;
        k_2 = _S1230 * _S1537;
        _S1234 = 0.0f;
        _S1235 = _S1538;
    }
    DiffPair_float_0 _S1539;
    (&_S1539)->primal_0 = _S1230;
    (&_S1539)->differential_0 = 0.0f;
    DiffPair_float_0 _S1540;
    (&_S1540)->primal_0 = _S1231;
    (&_S1540)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1539, &_S1540, k_2);
    float _S1541 = _S1540.differential_0 + _S1234;
    float _S1542 = _S1539.differential_0 + _S1235;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1543;
    (&_S1543)->primal_0 = _S1229;
    (&_S1543)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1543, _S1542);
    float2  _S1544 = _S1543.differential_0 + _S1532;
    float _S1545 = fx_5 * _S1497;
    float2  _S1546 = make_float2 (_S1545, fy_5 * _S1496) + make_float2 (dist_coeffs_5[int(8)] * _S1545, dist_coeffs_5[int(9)] * _S1545);
    float2  _S1547 = _S1217 * _S1546;
    float _S1548 = dist_coeffs_5[int(4)] * _S1546.y;
    float _S1549 = dist_coeffs_5[int(5)] * _S1546.x;
    float _S1550 = _S1547.x + _S1547.y;
    float _S1551 = r2_61 * _S1550;
    float _S1552 = r2_61 * _S1551;
    float _S1553 = dist_coeffs_5[int(7)] * _S1546.y + _S1548 + dist_coeffs_5[int(6)] * _S1546.x + _S1549 + _S1220 * _S1550 + _S1219 * _S1551 + _S1218 * _S1552 + dist_coeffs_5[int(3)] * (r2_61 * _S1552);
    float _S1554 = v_61 * _S1553;
    float _S1555 = u_61 * _S1553;
    float2  _S1556 = make_float2 (radial_29) * _S1546 + make_float2 (_S1102 * (v_61 * _S1546.y) + _S1222 * _S1549 + 2.0f * (u_61 * _S1549) + _S1099 * (v_61 * _S1546.x) + _S1555 + _S1555, _S1224 * _S1548 + 2.0f * (v_61 * _S1548) + _S1223 * _S1546.y + _S1221 * _S1546.x + _S1554 + _S1554);
    float3  _S1557 = make_float3 (_S1544.x, _S1544.y, _S1541);
    float2  _S1558 = _S1205 * _S1556;
    float2  _S1559 = _S1216 * _S1556;
    float _S1560 = _S1558.x + _S1558.y;
    if(_S1209)
    {
        float _S1561 = _S1560 / _S1210;
        float _S1562 = _S1211 * - _S1561;
        float _S1563 = _S1208 * (0.3333333432674408f * - (_S1207 * _S1561));
        k_2 = _S1563 + _S1563;
        _S1210 = _S1562;
        _S1211 = 0.0f;
    }
    else
    {
        float _S1564 = _S1560 / _S1212;
        float _S1565 = _S1208 * - _S1564;
        k_2 = _S1206 * _S1564;
        _S1210 = 0.0f;
        _S1211 = _S1565;
    }
    DiffPair_float_0 _S1566;
    (&_S1566)->primal_0 = _S1206;
    (&_S1566)->differential_0 = 0.0f;
    DiffPair_float_0 _S1567;
    (&_S1567)->primal_0 = _S1207;
    (&_S1567)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1566, &_S1567, k_2);
    float _S1568 = _S1567.differential_0 + _S1210;
    float _S1569 = _S1566.differential_0 + _S1211;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1570;
    (&_S1570)->primal_0 = _S1205;
    (&_S1570)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1570, _S1569);
    float2  _S1571 = _S1570.differential_0 + _S1559;
    float _S1572 = fx_5 * _S1499;
    float2  _S1573 = make_float2 (_S1572, fy_5 * _S1498) + make_float2 (dist_coeffs_5[int(8)] * _S1572, dist_coeffs_5[int(9)] * _S1572);
    float2  _S1574 = _S1193 * _S1573;
    float _S1575 = dist_coeffs_5[int(4)] * _S1573.y;
    float _S1576 = dist_coeffs_5[int(5)] * _S1573.x;
    float _S1577 = _S1574.x + _S1574.y;
    float _S1578 = r2_60 * _S1577;
    float _S1579 = r2_60 * _S1578;
    float _S1580 = dist_coeffs_5[int(7)] * _S1573.y + _S1575 + dist_coeffs_5[int(6)] * _S1573.x + _S1576 + _S1196 * _S1577 + _S1195 * _S1578 + _S1194 * _S1579 + dist_coeffs_5[int(3)] * (r2_60 * _S1579);
    float _S1581 = v_60 * _S1580;
    float _S1582 = u_60 * _S1580;
    float2  _S1583 = make_float2 (radial_28) * _S1573 + make_float2 (_S1102 * (v_60 * _S1573.y) + _S1198 * _S1576 + 2.0f * (u_60 * _S1576) + _S1099 * (v_60 * _S1573.x) + _S1582 + _S1582, _S1200 * _S1575 + 2.0f * (v_60 * _S1575) + _S1199 * _S1573.y + _S1197 * _S1573.x + _S1581 + _S1581);
    float3  _S1584 = make_float3 (_S1571.x, _S1571.y, _S1568);
    float2  _S1585 = _S1181 * _S1583;
    float2  _S1586 = _S1192 * _S1583;
    float _S1587 = _S1585.x + _S1585.y;
    if(_S1185)
    {
        float _S1588 = _S1587 / _S1186;
        float _S1589 = _S1187 * - _S1588;
        float _S1590 = _S1184 * (0.3333333432674408f * - (_S1183 * _S1588));
        k_2 = _S1590 + _S1590;
        _S1186 = _S1589;
        _S1187 = 0.0f;
    }
    else
    {
        float _S1591 = _S1587 / _S1188;
        float _S1592 = _S1184 * - _S1591;
        k_2 = _S1182 * _S1591;
        _S1186 = 0.0f;
        _S1187 = _S1592;
    }
    DiffPair_float_0 _S1593;
    (&_S1593)->primal_0 = _S1182;
    (&_S1593)->differential_0 = 0.0f;
    DiffPair_float_0 _S1594;
    (&_S1594)->primal_0 = _S1183;
    (&_S1594)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1593, &_S1594, k_2);
    float _S1595 = _S1594.differential_0 + _S1186;
    float _S1596 = _S1593.differential_0 + _S1187;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1597;
    (&_S1597)->primal_0 = _S1181;
    (&_S1597)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1597, _S1596);
    float2  _S1598 = _S1597.differential_0 + _S1586;
    float _S1599 = fx_5 * _S1501;
    float2  _S1600 = make_float2 (_S1599, fy_5 * _S1500) + make_float2 (dist_coeffs_5[int(8)] * _S1599, dist_coeffs_5[int(9)] * _S1599);
    float2  _S1601 = _S1169 * _S1600;
    float _S1602 = dist_coeffs_5[int(4)] * _S1600.y;
    float _S1603 = dist_coeffs_5[int(5)] * _S1600.x;
    float _S1604 = _S1601.x + _S1601.y;
    float _S1605 = r2_59 * _S1604;
    float _S1606 = r2_59 * _S1605;
    float _S1607 = dist_coeffs_5[int(7)] * _S1600.y + _S1602 + dist_coeffs_5[int(6)] * _S1600.x + _S1603 + _S1172 * _S1604 + _S1171 * _S1605 + _S1170 * _S1606 + dist_coeffs_5[int(3)] * (r2_59 * _S1606);
    float _S1608 = v_59 * _S1607;
    float _S1609 = u_59 * _S1607;
    float2  _S1610 = make_float2 (radial_27) * _S1600 + make_float2 (_S1102 * (v_59 * _S1600.y) + _S1174 * _S1603 + 2.0f * (u_59 * _S1603) + _S1099 * (v_59 * _S1600.x) + _S1609 + _S1609, _S1176 * _S1602 + 2.0f * (v_59 * _S1602) + _S1175 * _S1600.y + _S1173 * _S1600.x + _S1608 + _S1608);
    float3  _S1611 = make_float3 (_S1598.x, _S1598.y, _S1595);
    float2  _S1612 = _S1157 * _S1610;
    float2  _S1613 = _S1168 * _S1610;
    float _S1614 = _S1612.x + _S1612.y;
    if(_S1161)
    {
        float _S1615 = _S1614 / _S1162;
        float _S1616 = _S1163 * - _S1615;
        float _S1617 = _S1160 * (0.3333333432674408f * - (_S1159 * _S1615));
        k_2 = _S1617 + _S1617;
        _S1162 = _S1616;
        _S1163 = 0.0f;
    }
    else
    {
        float _S1618 = _S1614 / _S1164;
        float _S1619 = _S1160 * - _S1618;
        k_2 = _S1158 * _S1618;
        _S1162 = 0.0f;
        _S1163 = _S1619;
    }
    DiffPair_float_0 _S1620;
    (&_S1620)->primal_0 = _S1158;
    (&_S1620)->differential_0 = 0.0f;
    DiffPair_float_0 _S1621;
    (&_S1621)->primal_0 = _S1159;
    (&_S1621)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1620, &_S1621, k_2);
    float _S1622 = _S1621.differential_0 + _S1162;
    float _S1623 = _S1620.differential_0 + _S1163;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1624;
    (&_S1624)->primal_0 = _S1157;
    (&_S1624)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1624, _S1623);
    float2  _S1625 = _S1624.differential_0 + _S1613;
    float _S1626 = fx_5 * _S1492;
    float2  _S1627 = make_float2 (_S1626, fy_5 * _S1502) + make_float2 (dist_coeffs_5[int(8)] * _S1626, dist_coeffs_5[int(9)] * _S1626);
    float2  _S1628 = _S1145 * _S1627;
    float _S1629 = dist_coeffs_5[int(4)] * _S1627.y;
    float _S1630 = dist_coeffs_5[int(5)] * _S1627.x;
    float _S1631 = _S1628.x + _S1628.y;
    float _S1632 = r2_58 * _S1631;
    float _S1633 = r2_58 * _S1632;
    float _S1634 = dist_coeffs_5[int(7)] * _S1627.y + _S1629 + dist_coeffs_5[int(6)] * _S1627.x + _S1630 + _S1148 * _S1631 + _S1147 * _S1632 + _S1146 * _S1633 + dist_coeffs_5[int(3)] * (r2_58 * _S1633);
    float _S1635 = v_58 * _S1634;
    float _S1636 = u_58 * _S1634;
    float2  _S1637 = make_float2 (radial_26) * _S1627 + make_float2 (_S1102 * (v_58 * _S1627.y) + _S1150 * _S1630 + 2.0f * (u_58 * _S1630) + _S1099 * (v_58 * _S1627.x) + _S1636 + _S1636, _S1152 * _S1629 + 2.0f * (v_58 * _S1629) + _S1151 * _S1627.y + _S1149 * _S1627.x + _S1635 + _S1635);
    float3  _S1638 = make_float3 (_S1625.x, _S1625.y, _S1622);
    float2  _S1639 = _S1133 * _S1637;
    float2  _S1640 = _S1144 * _S1637;
    float _S1641 = _S1639.x + _S1639.y;
    if(_S1137)
    {
        float _S1642 = _S1641 / _S1138;
        float _S1643 = _S1139 * - _S1642;
        float _S1644 = _S1136 * (0.3333333432674408f * - (_S1135 * _S1642));
        k_2 = _S1644 + _S1644;
        _S1138 = _S1643;
        _S1139 = 0.0f;
    }
    else
    {
        float _S1645 = _S1641 / _S1140;
        float _S1646 = _S1136 * - _S1645;
        k_2 = _S1134 * _S1645;
        _S1138 = 0.0f;
        _S1139 = _S1646;
    }
    DiffPair_float_0 _S1647;
    (&_S1647)->primal_0 = _S1134;
    (&_S1647)->differential_0 = 0.0f;
    DiffPair_float_0 _S1648;
    (&_S1648)->primal_0 = _S1135;
    (&_S1648)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1647, &_S1648, k_2);
    float _S1649 = _S1648.differential_0 + _S1138;
    float _S1650 = _S1647.differential_0 + _S1139;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1651;
    (&_S1651)->primal_0 = _S1133;
    (&_S1651)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1651, _S1650);
    float2  _S1652 = _S1651.differential_0 + _S1640;
    float _S1653 = fx_5 * _S1489;
    float2  _S1654 = make_float2 (_S1653, fy_5 * _S1491) + make_float2 (dist_coeffs_5[int(8)] * _S1653, dist_coeffs_5[int(9)] * _S1653);
    float2  _S1655 = _S1121 * _S1654;
    float _S1656 = dist_coeffs_5[int(4)] * _S1654.y;
    float _S1657 = dist_coeffs_5[int(5)] * _S1654.x;
    float _S1658 = _S1655.x + _S1655.y;
    float _S1659 = r2_57 * _S1658;
    float _S1660 = r2_57 * _S1659;
    float _S1661 = dist_coeffs_5[int(7)] * _S1654.y + _S1656 + dist_coeffs_5[int(6)] * _S1654.x + _S1657 + _S1124 * _S1658 + _S1123 * _S1659 + _S1122 * _S1660 + dist_coeffs_5[int(3)] * (r2_57 * _S1660);
    float _S1662 = v_57 * _S1661;
    float _S1663 = u_57 * _S1661;
    float2  _S1664 = make_float2 (radial_25) * _S1654 + make_float2 (_S1102 * (v_57 * _S1654.y) + _S1126 * _S1657 + 2.0f * (u_57 * _S1657) + _S1099 * (v_57 * _S1654.x) + _S1663 + _S1663, _S1128 * _S1656 + 2.0f * (v_57 * _S1656) + _S1127 * _S1654.y + _S1125 * _S1654.x + _S1662 + _S1662);
    float3  _S1665 = make_float3 (_S1652.x, _S1652.y, _S1649);
    float2  _S1666 = _S1109 * _S1664;
    float2  _S1667 = _S1120 * _S1664;
    float _S1668 = _S1666.x + _S1666.y;
    if(_S1113)
    {
        float _S1669 = _S1668 / _S1114;
        float _S1670 = _S1115 * - _S1669;
        float _S1671 = _S1112 * (0.3333333432674408f * - (_S1111 * _S1669));
        k_2 = _S1671 + _S1671;
        _S1114 = _S1670;
        _S1115 = 0.0f;
    }
    else
    {
        float _S1672 = _S1668 / _S1116;
        float _S1673 = _S1112 * - _S1672;
        k_2 = _S1110 * _S1672;
        _S1114 = 0.0f;
        _S1115 = _S1673;
    }
    DiffPair_float_0 _S1674;
    (&_S1674)->primal_0 = _S1110;
    (&_S1674)->differential_0 = 0.0f;
    DiffPair_float_0 _S1675;
    (&_S1675)->primal_0 = _S1111;
    (&_S1675)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1674, &_S1675, k_2);
    float _S1676 = _S1675.differential_0 + _S1114;
    float _S1677 = _S1674.differential_0 + _S1115;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1678;
    (&_S1678)->primal_0 = _S1109;
    (&_S1678)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1678, _S1677);
    float2  _S1679 = _S1678.differential_0 + _S1667;
    float _S1680 = fx_5 * _S1493;
    float2  _S1681 = make_float2 (_S1680, fy_5 * _S1490) + make_float2 (dist_coeffs_5[int(8)] * _S1680, dist_coeffs_5[int(9)] * _S1680);
    float2  _S1682 = _S1095 * _S1681;
    float _S1683 = dist_coeffs_5[int(4)] * _S1681.y;
    float _S1684 = dist_coeffs_5[int(5)] * _S1681.x;
    float _S1685 = _S1682.x + _S1682.y;
    float _S1686 = r2_56 * _S1685;
    float _S1687 = r2_56 * _S1686;
    float _S1688 = dist_coeffs_5[int(7)] * _S1681.y + _S1683 + dist_coeffs_5[int(6)] * _S1681.x + _S1684 + _S1098 * _S1685 + _S1097 * _S1686 + _S1096 * _S1687 + dist_coeffs_5[int(3)] * (r2_56 * _S1687);
    float _S1689 = v_56 * _S1688;
    float _S1690 = u_56 * _S1688;
    float2  _S1691 = make_float2 (radial_24) * _S1681 + make_float2 (_S1102 * (v_56 * _S1681.y) + _S1101 * _S1684 + 2.0f * (u_56 * _S1684) + _S1099 * (v_56 * _S1681.x) + _S1690 + _S1690, _S1104 * _S1683 + 2.0f * (v_56 * _S1683) + _S1103 * _S1681.y + _S1100 * _S1681.x + _S1689 + _S1689);
    float3  _S1692 = make_float3 (_S1679.x, _S1679.y, _S1676);
    float2  _S1693 = _S1083 * _S1691;
    float2  _S1694 = _S1094 * _S1691;
    float _S1695 = _S1693.x + _S1693.y;
    if(_S1087)
    {
        float _S1696 = _S1695 / _S1088;
        float _S1697 = _S1089 * - _S1696;
        float _S1698 = _S1086 * (0.3333333432674408f * - (_S1085 * _S1696));
        k_2 = _S1698 + _S1698;
        _S1088 = _S1697;
        _S1089 = 0.0f;
    }
    else
    {
        float _S1699 = _S1695 / _S1090;
        float _S1700 = _S1086 * - _S1699;
        k_2 = _S1084 * _S1699;
        _S1088 = 0.0f;
        _S1089 = _S1700;
    }
    DiffPair_float_0 _S1701;
    (&_S1701)->primal_0 = _S1084;
    (&_S1701)->differential_0 = 0.0f;
    DiffPair_float_0 _S1702;
    (&_S1702)->primal_0 = _S1085;
    (&_S1702)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1701, &_S1702, k_2);
    float _S1703 = _S1702.differential_0 + _S1088;
    float _S1704 = _S1701.differential_0 + _S1089;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1705;
    (&_S1705)->primal_0 = _S1083;
    (&_S1705)->differential_0 = _S1515;
    s_bwd_length_impl_1(&_S1705, _S1704);
    float2  _S1706 = _S1705.differential_0 + _S1694;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1707;
    (&_S1707)->primal_0 = R_5;
    (&_S1707)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1708;
    (&_S1708)->primal_0 = _S1082;
    (&_S1708)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1707, &_S1708, _S1488);
    DiffPair_float_0 _S1709;
    (&_S1709)->primal_0 = _S1079;
    (&_S1709)->differential_0 = 0.0f;
    DiffPair_float_0 _S1710;
    (&_S1710)->primal_0 = _S1081;
    (&_S1710)->differential_0 = 0.0f;
    _d_max_0(&_S1709, &_S1710, 0.0f);
    DiffPair_float_0 _S1711;
    (&_S1711)->primal_0 = _S1078;
    (&_S1711)->differential_0 = 0.0f;
    DiffPair_float_0 _S1712;
    (&_S1712)->primal_0 = _S1081;
    (&_S1712)->differential_0 = 0.0f;
    _d_min_0(&_S1711, &_S1712, 0.0f);
    float _S1713 = _S1710.differential_0 + _S1712.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1714;
    (&_S1714)->primal_0 = _S1080;
    (&_S1714)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1714, _S1713);
    float3  _S1715 = _S1714.differential_0 + _S1530;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1716;
    (&_S1716)->primal_0 = R_5;
    (&_S1716)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1717;
    (&_S1717)->primal_0 = pos_i_13;
    (&_S1717)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1716, &_S1717, _S1715);
    DiffPair_float_0 _S1718;
    (&_S1718)->primal_0 = _S1075;
    (&_S1718)->differential_0 = 0.0f;
    DiffPair_float_0 _S1719;
    (&_S1719)->primal_0 = _S1077;
    (&_S1719)->differential_0 = 0.0f;
    _d_max_0(&_S1718, &_S1719, _S1709.differential_0);
    DiffPair_float_0 _S1720;
    (&_S1720)->primal_0 = _S1074;
    (&_S1720)->differential_0 = 0.0f;
    DiffPair_float_0 _S1721;
    (&_S1721)->primal_0 = _S1077;
    (&_S1721)->differential_0 = 0.0f;
    _d_min_0(&_S1720, &_S1721, _S1711.differential_0);
    float _S1722 = _S1719.differential_0 + _S1721.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1723;
    (&_S1723)->primal_0 = _S1076;
    (&_S1723)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1723, _S1722);
    float3  _S1724 = _S1723.differential_0 + _S1557;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1725;
    (&_S1725)->primal_0 = R_5;
    (&_S1725)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1726;
    (&_S1726)->primal_0 = pos_i_12;
    (&_S1726)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1725, &_S1726, _S1724);
    DiffPair_float_0 _S1727;
    (&_S1727)->primal_0 = _S1071;
    (&_S1727)->differential_0 = 0.0f;
    DiffPair_float_0 _S1728;
    (&_S1728)->primal_0 = _S1073;
    (&_S1728)->differential_0 = 0.0f;
    _d_max_0(&_S1727, &_S1728, _S1718.differential_0);
    DiffPair_float_0 _S1729;
    (&_S1729)->primal_0 = _S1070;
    (&_S1729)->differential_0 = 0.0f;
    DiffPair_float_0 _S1730;
    (&_S1730)->primal_0 = _S1073;
    (&_S1730)->differential_0 = 0.0f;
    _d_min_0(&_S1729, &_S1730, _S1720.differential_0);
    float _S1731 = _S1728.differential_0 + _S1730.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1732;
    (&_S1732)->primal_0 = _S1072;
    (&_S1732)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1732, _S1731);
    float3  _S1733 = _S1732.differential_0 + _S1584;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1734;
    (&_S1734)->primal_0 = R_5;
    (&_S1734)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1735;
    (&_S1735)->primal_0 = pos_i_11;
    (&_S1735)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1734, &_S1735, _S1733);
    DiffPair_float_0 _S1736;
    (&_S1736)->primal_0 = _S1067;
    (&_S1736)->differential_0 = 0.0f;
    DiffPair_float_0 _S1737;
    (&_S1737)->primal_0 = _S1069;
    (&_S1737)->differential_0 = 0.0f;
    _d_max_0(&_S1736, &_S1737, _S1727.differential_0);
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = _S1066;
    (&_S1738)->differential_0 = 0.0f;
    DiffPair_float_0 _S1739;
    (&_S1739)->primal_0 = _S1069;
    (&_S1739)->differential_0 = 0.0f;
    _d_min_0(&_S1738, &_S1739, _S1729.differential_0);
    float _S1740 = _S1737.differential_0 + _S1739.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1741;
    (&_S1741)->primal_0 = _S1068;
    (&_S1741)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1741, _S1740);
    float3  _S1742 = _S1741.differential_0 + _S1611;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1743;
    (&_S1743)->primal_0 = R_5;
    (&_S1743)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1744;
    (&_S1744)->primal_0 = pos_i_10;
    (&_S1744)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1743, &_S1744, _S1742);
    DiffPair_float_0 _S1745;
    (&_S1745)->primal_0 = _S1063;
    (&_S1745)->differential_0 = 0.0f;
    DiffPair_float_0 _S1746;
    (&_S1746)->primal_0 = _S1065;
    (&_S1746)->differential_0 = 0.0f;
    _d_max_0(&_S1745, &_S1746, _S1736.differential_0);
    DiffPair_float_0 _S1747;
    (&_S1747)->primal_0 = _S1062;
    (&_S1747)->differential_0 = 0.0f;
    DiffPair_float_0 _S1748;
    (&_S1748)->primal_0 = _S1065;
    (&_S1748)->differential_0 = 0.0f;
    _d_min_0(&_S1747, &_S1748, _S1738.differential_0);
    float _S1749 = _S1746.differential_0 + _S1748.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1750;
    (&_S1750)->primal_0 = _S1064;
    (&_S1750)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1750, _S1749);
    float3  _S1751 = _S1750.differential_0 + _S1638;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1752;
    (&_S1752)->primal_0 = R_5;
    (&_S1752)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1753;
    (&_S1753)->primal_0 = pos_i_9;
    (&_S1753)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1752, &_S1753, _S1751);
    DiffPair_float_0 _S1754;
    (&_S1754)->primal_0 = _S1059;
    (&_S1754)->differential_0 = 0.0f;
    DiffPair_float_0 _S1755;
    (&_S1755)->primal_0 = _S1061;
    (&_S1755)->differential_0 = 0.0f;
    _d_max_0(&_S1754, &_S1755, _S1745.differential_0);
    DiffPair_float_0 _S1756;
    (&_S1756)->primal_0 = _S1058;
    (&_S1756)->differential_0 = 0.0f;
    DiffPair_float_0 _S1757;
    (&_S1757)->primal_0 = _S1061;
    (&_S1757)->differential_0 = 0.0f;
    _d_min_0(&_S1756, &_S1757, _S1747.differential_0);
    float _S1758 = _S1755.differential_0 + _S1757.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1759;
    (&_S1759)->primal_0 = _S1060;
    (&_S1759)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1759, _S1758);
    float3  _S1760 = _S1759.differential_0 + _S1665;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1761;
    (&_S1761)->primal_0 = R_5;
    (&_S1761)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1762;
    (&_S1762)->primal_0 = pos_i_8;
    (&_S1762)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1761, &_S1762, _S1760);
    DiffPair_float_0 _S1763;
    (&_S1763)->primal_0 = _S1055;
    (&_S1763)->differential_0 = 0.0f;
    DiffPair_float_0 _S1764;
    (&_S1764)->primal_0 = _S1057;
    (&_S1764)->differential_0 = 0.0f;
    _d_max_0(&_S1763, &_S1764, _S1754.differential_0);
    DiffPair_float_0 _S1765;
    (&_S1765)->primal_0 = _S1054;
    (&_S1765)->differential_0 = 0.0f;
    DiffPair_float_0 _S1766;
    (&_S1766)->primal_0 = _S1057;
    (&_S1766)->differential_0 = 0.0f;
    _d_min_0(&_S1765, &_S1766, _S1756.differential_0);
    float _S1767 = _S1764.differential_0 + _S1766.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1768;
    (&_S1768)->primal_0 = _S1056;
    (&_S1768)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1768, _S1767);
    float3  _S1769 = _S1768.differential_0 + _S1692;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1770;
    (&_S1770)->primal_0 = R_5;
    (&_S1770)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1771;
    (&_S1771)->primal_0 = pos_i_7;
    (&_S1771)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1770, &_S1771, _S1769);
    DiffPair_float_0 _S1772;
    (&_S1772)->primal_0 = 0.0f;
    (&_S1772)->differential_0 = 0.0f;
    DiffPair_float_0 _S1773;
    (&_S1773)->primal_0 = _S1053;
    (&_S1773)->differential_0 = 0.0f;
    _d_max_0(&_S1772, &_S1773, _S1763.differential_0);
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = 1.00000001504746622e+30f;
    (&_S1774)->differential_0 = 0.0f;
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = _S1053;
    (&_S1775)->differential_0 = 0.0f;
    _d_min_0(&_S1774, &_S1775, _S1765.differential_0);
    float _S1776 = _S1773.differential_0 + _S1775.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1777;
    (&_S1777)->primal_0 = _S1052;
    (&_S1777)->differential_0 = _S1331;
    s_bwd_length_impl_0(&_S1777, _S1776);
    float3  _S1778 = _S1777.differential_0 + make_float3 (_S1706.x, _S1706.y, _S1703);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1779;
    (&_S1779)->primal_0 = R_5;
    (&_S1779)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1780;
    (&_S1780)->primal_0 = pos_5;
    (&_S1780)->differential_0 = _S1331;
    s_bwd_prop_mul_0(&_S1779, &_S1780, _S1778);
    float3  _S1781 = _S1488 + _S1715 + _S1724 + _S1733 + _S1742 + _S1751 + _S1760 + _S1769 + _S1778 + _S1396.differential_0;
    Matrix<float, 3, 3>  _S1782 = _S1707.differential_0 + _S1716.differential_0 + _S1725.differential_0 + _S1734.differential_0 + _S1743.differential_0 + _S1752.differential_0 + _S1761.differential_0 + _S1770.differential_0 + _S1779.differential_0 + _S1397;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_1)[int(0)] = _S1472;
    (*v_sh_coeffs_1)[int(1)] = _S1473;
    (*v_sh_coeffs_1)[int(2)] = _S1474;
    (*v_sh_coeffs_1)[int(3)] = _S1475;
    (*v_sh_coeffs_1)[int(4)] = _S1476;
    (*v_sh_coeffs_1)[int(5)] = _S1477;
    (*v_sh_coeffs_1)[int(6)] = _S1478;
    (*v_sh_coeffs_1)[int(7)] = _S1479;
    (*v_sh_coeffs_1)[int(8)] = _S1480;
    (*v_sh_coeffs_1)[int(9)] = _S1481;
    (*v_sh_coeffs_1)[int(10)] = _S1482;
    (*v_sh_coeffs_1)[int(11)] = _S1483;
    (*v_sh_coeffs_1)[int(12)] = _S1484;
    (*v_sh_coeffs_1)[int(13)] = _S1485;
    (*v_sh_coeffs_1)[int(14)] = _S1486;
    (*v_sh_coeffs_1)[int(15)] = _S1487;
    *v_R_1 = _S1782;
    *v_t_1 = _S1781;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  dOut_7)
{
    float3  _S1783 = _slang_select(((*dpx_8).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_8).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_7;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S1783;
    return;
}

inline __device__ float3  abs_0(float3  x_13)
{
    float3  result_8;
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
        *_slang_vector_get_element_ptr(&result_8, i_4) = (F32_abs((_slang_vector_get_element(x_13, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_8;
}

inline __device__ bool ray_aabb_intersection(float3  ray_o_0, float3  ray_d_0, float3  center_0, float radius_0, float * t0_0, float * t1_0)
{
    float3  m_1 = make_float3 (1.0f) / ray_d_0;
    float3  k_3 = abs_0(m_1) * make_float3 (radius_0);
    float3  _S1784 = - (m_1 * (ray_o_0 - center_0));
    float3  ta_0 = _S1784 - k_3;
    float3  tb_0 = _S1784 + k_3;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S1785 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S1785;
    return (*t0_0) < _S1785;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_9, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_8)
{
    float _S1786 = (1.0f - (*dps_0).primal_0) * dOut_8;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S1786;
    DiffPair_float_0 _S1787 = *dpy_4;
    float _S1788 = (*dps_0).primal_0 * dOut_8;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S1788;
    float _S1789 = (_S1787.primal_0 - (*dpx_9).primal_0) * dOut_8;
    dps_0->primal_0 = _S1787.primal_0;
    dps_0->differential_0 = _S1789;
    return;
}

inline __device__ float lerp_0(float x_14, float y_9, float s_0)
{
    return x_14 + (y_9 - x_14) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_0)
{
    float _S1790 = w_0.z;
    float _S1791 = 1.0f - _S1790;
    float _S1792 = w_0.y;
    float _S1793 = 1.0f - _S1792;
    float _S1794 = _S1791 * _S1793;
    float _S1795 = w_0.x;
    float _S1796 = 1.0f - _S1795;
    float _S1797 = _S1791 * _S1792;
    float _S1798 = _S1790 * _S1793;
    float _S1799 = _S1790 * _S1792;
    return _S1794 * _S1796 * (*densities_6)[int(0)] + _S1794 * _S1795 * (*densities_6)[int(1)] + _S1797 * _S1796 * (*densities_6)[int(2)] + _S1797 * _S1795 * (*densities_6)[int(3)] + _S1798 * _S1796 * (*densities_6)[int(4)] + _S1798 * _S1795 * (*densities_6)[int(5)] + _S1799 * _S1796 * (*densities_6)[int(6)] + _S1799 * _S1795 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  densities_7, float3  ray_o_1, float3  ray_d_1)
{
    float _S1800 = 0.5f * size_6;
    float3  m_2 = make_float3 (1.0f) / ray_d_1;
    float3  k_4 = abs_0(m_2) * make_float3 (_S1800);
    float3  _S1801 = - (m_2 * (ray_o_1 - (pos_6 + make_float3 (_S1800))));
    float3  ta_1 = _S1801 - k_4;
    float3  tb_1 = _S1801 + k_4;
    float t0_1 = (F32_max(((F32_max((ta_1.x), (ta_1.y)))), ((F32_max((ta_1.z), (0.0f))))));
    float t1_1 = (F32_min(((F32_min((tb_1.x), (tb_1.y)))), (tb_1.z)));
    if(!(t0_1 < t1_1))
    {
        return 0.0f;
    }
    int i_5 = int(0);
    float accum_0 = 0.0f;
    for(;;)
    {
        if(i_5 < int(8))
        {
        }
        else
        {
            break;
        }
        float3  _S1802 = (ray_o_1 + ray_d_1 * make_float3 (lerp_0(t0_1, t1_1, (float(i_5) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6);
        FixedArray<float, 8>  _S1803 = densities_7;
        float _S1804 = interp_0(&_S1803, _S1802);
        float _S1805;
        if(_S1804 > 1.10000002384185791f)
        {
            _S1805 = _S1804;
        }
        else
        {
            _S1805 = (F32_exp((0.90909093618392944f * _S1804 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S1805;
        i_5 = i_5 + int(1);
        accum_0 = accum_1;
    }
    return (F32_min((1.0f - (F32_exp((- (t1_1 - t0_1) / 8.0f * accum_0)))), (0.99900001287460327f)));
}

struct DiffPair_arrayx3Cfloatx2C8x3E_0
{
    FixedArray<float, 8>  primal_0;
    FixedArray<float, 8>  differential_0;
};

struct s_bwd_prop_evaluate_alpha_voxel_Intermediates_0
{
    float _S1806;
};

inline __device__ float3  s_primal_ctx_abs_0(float3  _S1807)
{
    return abs_0(_S1807);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1808, float _S1809, float _S1810)
{
    return lerp_0(_S1808, _S1809, _S1810);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S1811 = dpw_0.z;
    float _S1812 = 1.0f - _S1811;
    float _S1813 = dpw_0.y;
    float _S1814 = 1.0f - _S1813;
    float _S1815 = _S1812 * _S1814;
    float _S1816 = dpw_0.x;
    float _S1817 = 1.0f - _S1816;
    float _S1818 = _S1812 * _S1813;
    float _S1819 = _S1811 * _S1814;
    float _S1820 = _S1811 * _S1813;
    return _S1815 * _S1817 * (*dpdensities_0)[int(0)] + _S1815 * _S1816 * (*dpdensities_0)[int(1)] + _S1818 * _S1817 * (*dpdensities_0)[int(2)] + _S1818 * _S1816 * (*dpdensities_0)[int(3)] + _S1819 * _S1817 * (*dpdensities_0)[int(4)] + _S1819 * _S1816 * (*dpdensities_0)[int(5)] + _S1820 * _S1817 * (*dpdensities_0)[int(6)] + _S1820 * _S1816 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_exp_0(float _S1821)
{
    return (F32_exp((_S1821)));
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_0, float3  dpray_d_0, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S1806 = 0.0f;
    _s_diff_ctx_0->_S1806 = 0.0f;
    float _S1822 = 0.5f * size_7;
    float3  m_3 = make_float3 (1.0f) / dpray_d_0;
    float3  k_5 = s_primal_ctx_abs_0(m_3) * make_float3 (_S1822);
    float3  _S1823 = - (m_3 * (dpray_o_0 - (pos_7 + make_float3 (_S1822))));
    float3  ta_2 = _S1823 - k_5;
    float3  tb_2 = _S1823 + k_5;
    float _S1824 = (F32_max(((F32_max((ta_2.x), (ta_2.y)))), ((F32_max((ta_2.z), (0.0f))))));
    float _S1825 = (F32_min(((F32_min((tb_2.x), (tb_2.y)))), (tb_2.z)));
    float accum_2;
    if(!!(_S1824 < _S1825))
    {
        float _S1826 = - (_S1825 - _S1824);
        bool _runFlag_0 = true;
        int i_6 = int(0);
        accum_2 = 0.0f;
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_0->_S1806 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            int _S1827;
            float _S1828;
            if(i_6 < int(8))
            {
                float _S1829 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_0 + dpray_d_0 * make_float3 (s_primal_ctx_lerp_0(_S1824, _S1825, (float(i_6) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S1829 > 1.10000002384185791f)
                {
                    _S1828 = _S1829;
                }
                else
                {
                    _S1828 = s_primal_ctx_exp_0(0.90909093618392944f * _S1829 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S1828;
                _S1827 = int(2);
                _S1828 = accum_3;
            }
            else
            {
                _S1827 = int(1);
                _S1828 = 0.0f;
            }
            if(_S1827 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                i_6 = i_6 + int(1);
                accum_2 = _S1828;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = (F32_min((1.0f - s_primal_ctx_exp_0(_S1826 / 8.0f * accum_2)), (0.99900001287460327f)));
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1830, float _S1831)
{
    _d_exp_0(_S1830, _S1831);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_2)
{
    float _S1832 = (*dpw_1).primal_0.z;
    float _S1833 = 1.0f - _S1832;
    float _S1834 = (*dpw_1).primal_0.y;
    float _S1835 = 1.0f - _S1834;
    float _S1836 = _S1833 * _S1835;
    float _S1837 = (*dpw_1).primal_0.x;
    float _S1838 = 1.0f - _S1837;
    float _S1839 = _S1833 * _S1834;
    float _S1840 = _S1832 * _S1835;
    float _S1841 = _S1832 * _S1834;
    float _S1842 = _S1841 * _S1837 * _s_dOut_2;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_2;
    float _S1843 = _S1841 * _S1838 * _s_dOut_2;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_2;
    float _S1844 = _S1840 * _S1837 * _s_dOut_2;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_2;
    float _S1845 = _S1840 * _S1838 * _s_dOut_2;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_2;
    float _S1846 = _S1839 * _S1837 * _s_dOut_2;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_2;
    float _S1847 = _S1839 * _S1838 * _s_dOut_2;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_2;
    float _S1848 = _S1836 * _S1837 * _s_dOut_2;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_2;
    float _S1849 = _S1836 * _S1838 * _s_dOut_2;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_2;
    float _S1850 = _S1837 * s_diff_w7_T_0 + _S1838 * s_diff_w6_T_0;
    float _S1851 = _S1837 * s_diff_w5_T_0 + _S1838 * s_diff_w4_T_0;
    float _S1852 = _S1837 * s_diff_w3_T_0 + _S1838 * s_diff_w2_T_0;
    float _S1853 = _S1837 * s_diff_w1_T_0 + _S1838 * s_diff_w0_T_0;
    float3  _S1854 = make_float3 (_S1841 * s_diff_w7_T_0 + _S1840 * s_diff_w5_T_0 + _S1839 * s_diff_w3_T_0 + _S1836 * s_diff_w1_T_0 + - (_S1841 * s_diff_w6_T_0 + _S1840 * s_diff_w4_T_0 + _S1839 * s_diff_w2_T_0 + _S1836 * s_diff_w0_T_0), _S1832 * _S1850 + _S1833 * _S1852 + - (_S1832 * _S1851 + _S1833 * _S1853), _S1834 * _S1850 + _S1835 * _S1851 + - (_S1834 * _S1852 + _S1835 * _S1853));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S1854;
    FixedArray<float, 8>  _S1855;
    _S1855[int(0)] = 0.0f;
    _S1855[int(1)] = 0.0f;
    _S1855[int(2)] = 0.0f;
    _S1855[int(3)] = 0.0f;
    _S1855[int(4)] = 0.0f;
    _S1855[int(5)] = 0.0f;
    _S1855[int(6)] = 0.0f;
    _S1855[int(7)] = 0.0f;
    _S1855[int(7)] = _S1842;
    _S1855[int(6)] = _S1843;
    _S1855[int(5)] = _S1844;
    _S1855[int(4)] = _S1845;
    _S1855[int(3)] = _S1846;
    _S1855[int(2)] = _S1847;
    _S1855[int(1)] = _S1848;
    _S1855[int(0)] = _S1849;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S1855;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S1856, DiffPair_float_0 * _S1857, DiffPair_float_0 * _S1858, float _S1859)
{
    _d_lerp_0(_S1856, _S1857, _S1858, _S1859);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1860, float3  _S1861)
{
    _d_abs_vector_0(_S1860, _S1861);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float _s_dOut_3, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_1)
{
    FixedArray<float, 8>  _S1862 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1863 = *dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1864 = *dpray_d_1;
    float _S1865 = 0.5f * size_8;
    float3  _S1866 = make_float3 (_S1865);
    float3  _S1867 = make_float3 (size_8);
    float3  m_4 = make_float3 (1.0f) / (*dpray_d_1).primal_0;
    float3  _S1868 = (*dpray_d_1).primal_0 * (*dpray_d_1).primal_0;
    float3  _S1869 = (*dpray_o_1).primal_0 - (pos_8 + make_float3 (_S1865));
    float3  k_6 = s_primal_ctx_abs_0(m_4) * make_float3 (_S1865);
    float3  _S1870 = - (m_4 * _S1869);
    float3  ta_3 = _S1870 - k_6;
    float3  tb_3 = _S1870 + k_6;
    float _S1871 = ta_3.x;
    float _S1872 = ta_3.y;
    float _S1873 = (F32_max((_S1871), (_S1872)));
    float _S1874 = ta_3.z;
    float _S1875 = (F32_max((_S1874), (0.0f)));
    float _S1876 = (F32_max((_S1873), (_S1875)));
    float _S1877 = tb_3.x;
    float _S1878 = tb_3.y;
    float _S1879 = (F32_min((_S1877), (_S1878)));
    float _S1880 = tb_3.z;
    float _S1881 = (F32_min((_S1879), (_S1880)));
    bool _S1882 = !!(_S1876 < _S1881);
    float _S1883;
    float _S1884;
    float _S1885;
    if(_S1882)
    {
        float _S1886 = - (_S1881 - _S1876) / 8.0f;
        float _S1887 = _S1886 * _s_diff_ctx_1->_S1806;
        _S1883 = 1.0f - s_primal_ctx_exp_0(_S1887);
        _S1884 = _S1887;
        _S1885 = _S1886;
    }
    else
    {
        _S1883 = 0.0f;
        _S1884 = 0.0f;
        _S1885 = 0.0f;
    }
    float3  _S1888 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1889 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<float, 8>  _S1890;
    float3  _S1891;
    float3  _S1892;
    if(_S1882)
    {
        DiffPair_float_0 _S1893;
        (&_S1893)->primal_0 = _S1883;
        (&_S1893)->differential_0 = 0.0f;
        DiffPair_float_0 _S1894;
        (&_S1894)->primal_0 = 0.99900001287460327f;
        (&_S1894)->differential_0 = 0.0f;
        _d_min_0(&_S1893, &_S1894, _s_dOut_3);
        float _S1895 = - _S1893.differential_0;
        DiffPair_float_0 _S1896;
        (&_S1896)->primal_0 = _S1884;
        (&_S1896)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S1896, _S1895);
        float _S1897 = _S1885 * _S1896.differential_0;
        float _S1898 = 0.125f * (_s_diff_ctx_1->_S1806 * _S1896.differential_0);
        int _dc_0 = int(8);
        _S1883 = _S1897;
        _S1890[int(0)] = 0.0f;
        _S1890[int(1)] = 0.0f;
        _S1890[int(2)] = 0.0f;
        _S1890[int(3)] = 0.0f;
        _S1890[int(4)] = 0.0f;
        _S1890[int(5)] = 0.0f;
        _S1890[int(6)] = 0.0f;
        _S1890[int(7)] = 0.0f;
        _S1891 = _S1888;
        _S1892 = _S1888;
        _S1884 = 0.0f;
        _S1885 = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S1899 = _dc_0 < int(8);
            float _S1900;
            float _S1901;
            int _S1902;
            float3  _S1903;
            float3  _S1904;
            bool _S1905;
            if(_S1899)
            {
                float _S1906 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S1907 = s_primal_ctx_lerp_0(_S1876, _S1881, _S1906);
                float3  _S1908 = make_float3 (_S1907);
                float3  _S1909 = (_S1863.primal_0 + _S1864.primal_0 * make_float3 (_S1907) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S1910 = _S1862;
                float _S1911 = s_primal_ctx_interp_0(&_S1910, _S1909);
                bool _S1912 = _S1911 > 1.10000002384185791f;
                if(_S1912)
                {
                    _S1900 = 0.0f;
                }
                else
                {
                    _S1900 = 0.90909093618392944f * _S1911 - 0.90468984842300415f;
                }
                _S1902 = int(2);
                _S1905 = _S1912;
                _S1903 = _S1909;
                _S1904 = _S1908;
                _S1901 = _S1906;
            }
            else
            {
                _S1902 = int(1);
                _S1905 = false;
                _S1900 = 0.0f;
                _S1903 = _S1888;
                _S1904 = _S1888;
                _S1901 = 0.0f;
            }
            float _S1913;
            float _S1914;
            if(!(_S1902 != int(2)))
            {
                _S1913 = _S1883;
                _S1914 = 0.0f;
            }
            else
            {
                _S1913 = 0.0f;
                _S1914 = _S1883;
            }
            if(_S1899)
            {
                float _S1915 = _S1913 + _S1914;
                float _S1916;
                if(_S1905)
                {
                    _S1916 = _S1913;
                }
                else
                {
                    DiffPair_float_0 _S1917;
                    (&_S1917)->primal_0 = _S1900;
                    (&_S1917)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S1917, _S1913);
                    _S1916 = 0.90909093618392944f * _S1917.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S1918;
                (&_S1918)->primal_0 = _S1862;
                (&_S1918)->differential_0 = _S1889;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1919;
                (&_S1919)->primal_0 = _S1903;
                (&_S1919)->differential_0 = _S1888;
                s_bwd_prop_interp_0(&_S1918, &_S1919, _S1916);
                float3  _S1920 = _S1919.differential_0 / _S1867;
                float3  _S1921 = _S1864.primal_0 * _S1920;
                float3  _S1922 = _S1904 * _S1920;
                float _S1923 = _S1921.x + _S1921.y + _S1921.z;
                DiffPair_float_0 _S1924;
                (&_S1924)->primal_0 = _S1876;
                (&_S1924)->differential_0 = 0.0f;
                DiffPair_float_0 _S1925;
                (&_S1925)->primal_0 = _S1881;
                (&_S1925)->differential_0 = 0.0f;
                DiffPair_float_0 _S1926;
                (&_S1926)->primal_0 = _S1901;
                (&_S1926)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S1924, &_S1925, &_S1926, _S1923);
                float _S1927 = (&_S1918)->differential_0[int(0)] + _S1890[int(0)];
                float _S1928 = (&_S1918)->differential_0[int(1)] + _S1890[int(1)];
                float _S1929 = (&_S1918)->differential_0[int(2)] + _S1890[int(2)];
                float _S1930 = (&_S1918)->differential_0[int(3)] + _S1890[int(3)];
                float _S1931 = (&_S1918)->differential_0[int(4)] + _S1890[int(4)];
                float _S1932 = (&_S1918)->differential_0[int(5)] + _S1890[int(5)];
                float _S1933 = (&_S1918)->differential_0[int(6)] + _S1890[int(6)];
                float _S1934 = (&_S1918)->differential_0[int(7)] + _S1890[int(7)];
                float3  _S1935 = _S1920 + _S1891;
                float3  _S1936 = _S1922 + _S1892;
                float _S1937 = _S1925.differential_0 + _S1884;
                float _S1938 = _S1924.differential_0 + _S1885;
                _S1883 = _S1915;
                _S1890[int(0)] = _S1927;
                _S1890[int(1)] = _S1928;
                _S1890[int(2)] = _S1929;
                _S1890[int(3)] = _S1930;
                _S1890[int(4)] = _S1931;
                _S1890[int(5)] = _S1932;
                _S1890[int(6)] = _S1933;
                _S1890[int(7)] = _S1934;
                _S1891 = _S1935;
                _S1892 = _S1936;
                _S1884 = _S1937;
                _S1885 = _S1938;
            }
            else
            {
                _S1883 = _S1914;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S1939 = - _S1898;
        float _S1940 = - _S1939 + _S1885;
        float3  _S1941 = _S1891;
        _S1883 = _S1939 + _S1884;
        _S1884 = _S1940;
        _S1891 = _S1892;
        _S1892 = _S1941;
    }
    else
    {
        _S1883 = 0.0f;
        _S1884 = 0.0f;
        _S1891 = _S1888;
        _S1892 = _S1888;
        _S1890[int(0)] = 0.0f;
        _S1890[int(1)] = 0.0f;
        _S1890[int(2)] = 0.0f;
        _S1890[int(3)] = 0.0f;
        _S1890[int(4)] = 0.0f;
        _S1890[int(5)] = 0.0f;
        _S1890[int(6)] = 0.0f;
        _S1890[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S1942;
    (&_S1942)->primal_0 = _S1879;
    (&_S1942)->differential_0 = 0.0f;
    DiffPair_float_0 _S1943;
    (&_S1943)->primal_0 = _S1880;
    (&_S1943)->differential_0 = 0.0f;
    _d_min_0(&_S1942, &_S1943, _S1883);
    DiffPair_float_0 _S1944;
    (&_S1944)->primal_0 = _S1877;
    (&_S1944)->differential_0 = 0.0f;
    DiffPair_float_0 _S1945;
    (&_S1945)->primal_0 = _S1878;
    (&_S1945)->differential_0 = 0.0f;
    _d_min_0(&_S1944, &_S1945, _S1942.differential_0);
    DiffPair_float_0 _S1946;
    (&_S1946)->primal_0 = _S1873;
    (&_S1946)->differential_0 = 0.0f;
    DiffPair_float_0 _S1947;
    (&_S1947)->primal_0 = _S1875;
    (&_S1947)->differential_0 = 0.0f;
    _d_max_0(&_S1946, &_S1947, _S1884);
    DiffPair_float_0 _S1948;
    (&_S1948)->primal_0 = _S1874;
    (&_S1948)->differential_0 = 0.0f;
    DiffPair_float_0 _S1949;
    (&_S1949)->primal_0 = 0.0f;
    (&_S1949)->differential_0 = 0.0f;
    _d_max_0(&_S1948, &_S1949, _S1947.differential_0);
    DiffPair_float_0 _S1950;
    (&_S1950)->primal_0 = _S1871;
    (&_S1950)->differential_0 = 0.0f;
    DiffPair_float_0 _S1951;
    (&_S1951)->primal_0 = _S1872;
    (&_S1951)->differential_0 = 0.0f;
    _d_max_0(&_S1950, &_S1951, _S1946.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S1944.differential_0, _S1945.differential_0, _S1943.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S1950.differential_0, _S1951.differential_0, _S1948.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S1952 = _S1866 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1953;
    (&_S1953)->primal_0 = m_4;
    (&_S1953)->differential_0 = _S1888;
    s_bwd_prop_abs_0(&_S1953, _S1952);
    float3  _S1954 = m_4 * s_diff_n_T_0;
    float3  _S1955 = - ((_S1953.differential_0 + _S1869 * s_diff_n_T_0) / _S1868) + _S1891;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1955;
    float3  _S1956 = _S1954 + _S1892;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1956;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S1890;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S1957, float _S1958, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S1959, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1960, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1961, float _S1962)
{
    FixedArray<float, 8>  _S1963 = _S1959->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1964;
    float _S1965 = s_primal_ctx_evaluate_alpha_voxel_0(_S1957, _S1958, &_S1963, (*_S1960).primal_0, (*_S1961).primal_0, &_S1964);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1966 = _S1964;
    s_bwd_prop_evaluate_alpha_voxel_0(_S1957, _S1958, _S1959, _S1960, _S1961, _S1962, &_S1966);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  densities_8, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    FixedArray<float, 8>  _S1967 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = densities_8;
    (&dp_densities_0)->differential_0 = _S1967;
    float3  _S1968 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1968;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1968;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_0 = dp_ray_o_0.differential_0;
    *v_ray_d_0 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  densities_9, float3  rgb_0, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_4)
{
    *out_rgb_0 = rgb_0;
    float _S1969 = 0.5f * size_10;
    float3  m_5 = make_float3 (1.0f) / ray_d_3;
    float3  k_7 = abs_0(m_5) * make_float3 (_S1969);
    float3  _S1970 = - (m_5 * (ray_o_3 - (pos_10 + make_float3 (_S1969))));
    float3  ta_4 = _S1970 - k_7;
    float3  tb_4 = _S1970 + k_7;
    float _S1971 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S1972 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
    int i_7 = int(0);
    float accum_4 = 0.0f;
    float depth_accum_0 = 0.0f;
    for(;;)
    {
        if(i_7 < int(8))
        {
        }
        else
        {
            break;
        }
        float t_6 = lerp_0(_S1971, _S1972, (float(i_7) + 0.5f) / 8.0f);
        float3  _S1973 = (ray_o_3 + ray_d_3 * make_float3 (t_6) - pos_10) / make_float3 (size_10);
        FixedArray<float, 8>  _S1974 = densities_9;
        float _S1975 = interp_0(&_S1974, _S1973);
        float _S1976;
        if(_S1975 > 1.10000002384185791f)
        {
            _S1976 = _S1975;
        }
        else
        {
            _S1976 = (F32_exp((0.90909093618392944f * _S1975 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S1976;
        float depth_accum_1 = depth_accum_0 + t_6 * _S1976;
        i_7 = i_7 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_4 = (F32_max((depth_accum_0 / accum_4), (0.0f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S1977;
    float _S1978;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_0, float3  dpray_o_2, float3  dpray_d_2, float3  * dpout_rgb_0, float * dpdepth_0, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_2)
{
    _s_diff_ctx_2->_S1977 = 0.0f;
    _s_diff_ctx_2->_S1978 = 0.0f;
    float _S1979 = 0.5f * size_11;
    float3  m_6 = make_float3 (1.0f) / dpray_d_2;
    float3  k_8 = s_primal_ctx_abs_0(m_6) * make_float3 (_S1979);
    float3  _S1980 = - (m_6 * (dpray_o_2 - (pos_11 + make_float3 (_S1979))));
    float3  ta_5 = _S1980 - k_8;
    float3  tb_5 = _S1980 + k_8;
    float _S1981 = (F32_max(((F32_max((ta_5.x), (ta_5.y)))), ((F32_max((ta_5.z), (0.0f))))));
    float _S1982 = (F32_min(((F32_min((tb_5.x), (tb_5.y)))), (tb_5.z)));
    bool _runFlag_1 = true;
    int i_8 = int(0);
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_2->_S1977 = depth_accum_2;
        _s_diff_ctx_2->_S1978 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        int _S1983;
        float _S1984;
        float _S1985;
        if(i_8 < int(8))
        {
            float _S1986 = s_primal_ctx_lerp_0(_S1981, _S1982, (float(i_8) + 0.5f) / 8.0f);
            float _S1987 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_2 + dpray_d_2 * make_float3 (_S1986) - pos_11) / make_float3 (size_11));
            if(_S1987 > 1.10000002384185791f)
            {
                _S1984 = _S1987;
            }
            else
            {
                _S1984 = s_primal_ctx_exp_0(0.90909093618392944f * _S1987 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S1984;
            float depth_accum_3 = depth_accum_2 + _S1986 * _S1984;
            _S1983 = int(1);
            _S1984 = accum_7;
            _S1985 = depth_accum_3;
        }
        else
        {
            _S1983 = int(0);
            _S1984 = 0.0f;
            _S1985 = 0.0f;
        }
        if(_S1983 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            i_8 = i_8 + int(1);
            accum_6 = _S1984;
            depth_accum_2 = _S1985;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S1988 = (F32_max((depth_accum_2 / accum_6), (0.0f)));
    *dpout_rgb_0 = dprgb_0;
    *dpdepth_0 = _S1988;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_1, float dpdepth_1, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_3)
{
    FixedArray<float, 8>  _S1989 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1990 = *dpray_o_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1991 = *dpray_d_3;
    float3  _S1992 = make_float3 (size_12);
    float _S1993 = 0.5f * size_12;
    float3  _S1994 = make_float3 (_S1993);
    float3  m_7 = make_float3 (1.0f) / (*dpray_d_3).primal_0;
    float3  _S1995 = (*dpray_d_3).primal_0 * (*dpray_d_3).primal_0;
    float3  _S1996 = (*dpray_o_3).primal_0 - (pos_12 + make_float3 (_S1993));
    float3  k_9 = s_primal_ctx_abs_0(m_7) * make_float3 (_S1993);
    float3  _S1997 = - (m_7 * _S1996);
    float3  ta_6 = _S1997 - k_9;
    float3  tb_6 = _S1997 + k_9;
    float _S1998 = ta_6.x;
    float _S1999 = ta_6.y;
    float _S2000 = (F32_max((_S1998), (_S1999)));
    float _S2001 = ta_6.z;
    float _S2002 = (F32_max((_S2001), (0.0f)));
    float _S2003 = (F32_max((_S2000), (_S2002)));
    float _S2004 = tb_6.x;
    float _S2005 = tb_6.y;
    float _S2006 = (F32_min((_S2004), (_S2005)));
    float _S2007 = tb_6.z;
    float _S2008 = (F32_min((_S2006), (_S2007)));
    float _S2009 = _s_diff_ctx_3->_S1978 * _s_diff_ctx_3->_S1978;
    float3  _S2010 = make_float3 (0.0f);
    FixedArray<float, 8>  _S2011 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S2012;
    (&_S2012)->primal_0 = _s_diff_ctx_3->_S1977 / _s_diff_ctx_3->_S1978;
    (&_S2012)->differential_0 = 0.0f;
    DiffPair_float_0 _S2013;
    (&_S2013)->primal_0 = 0.0f;
    (&_S2013)->differential_0 = 0.0f;
    _d_max_0(&_S2012, &_S2013, dpdepth_1);
    float _S2014 = _S2012.differential_0 / _S2009;
    float _S2015 = _s_diff_ctx_3->_S1977 * - _S2014;
    float _S2016 = _s_diff_ctx_3->_S1978 * _S2014;
    int _dc_1 = int(8);
    float _S2017 = _S2015;
    float _S2018 = _S2016;
    FixedArray<float, 8>  _S2019;
    _S2019[int(0)] = 0.0f;
    _S2019[int(1)] = 0.0f;
    _S2019[int(2)] = 0.0f;
    _S2019[int(3)] = 0.0f;
    _S2019[int(4)] = 0.0f;
    _S2019[int(5)] = 0.0f;
    _S2019[int(6)] = 0.0f;
    _S2019[int(7)] = 0.0f;
    float3  _S2020 = _S2010;
    float3  _S2021 = _S2010;
    float _S2022 = 0.0f;
    float _S2023 = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S2024 = _dc_1 < int(8);
        int _S2025;
        float _S2026;
        float _S2027;
        float _S2028;
        float _S2029;
        float3  _S2030;
        float3  _S2031;
        bool _S2032;
        if(_S2024)
        {
            float _S2033 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S2034 = s_primal_ctx_lerp_0(_S2003, _S2008, _S2033);
            float3  _S2035 = make_float3 (_S2034);
            float3  _S2036 = (_S1990.primal_0 + _S1991.primal_0 * make_float3 (_S2034) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S2037 = _S1989;
            float _S2038 = s_primal_ctx_interp_0(&_S2037, _S2036);
            bool _S2039 = _S2038 > 1.10000002384185791f;
            if(_S2039)
            {
                _S2026 = _S2038;
                _S2027 = 0.0f;
            }
            else
            {
                float _S2040 = 0.90909093618392944f * _S2038 - 0.90468984842300415f;
                _S2026 = s_primal_ctx_exp_0(_S2040);
                _S2027 = _S2040;
            }
            float _S2041 = _S2026;
            float _S2042 = _S2027;
            _S2025 = int(1);
            _S2026 = _S2034;
            _S2027 = _S2041;
            _S2032 = _S2039;
            _S2028 = _S2042;
            _S2030 = _S2036;
            _S2031 = _S2035;
            _S2029 = _S2033;
        }
        else
        {
            _S2025 = int(0);
            _S2026 = 0.0f;
            _S2027 = 0.0f;
            _S2032 = false;
            _S2028 = 0.0f;
            _S2030 = _S2010;
            _S2031 = _S2010;
            _S2029 = 0.0f;
        }
        float _S2043;
        float _S2044;
        float _S2045;
        float _S2046;
        if(!(_S2025 != int(1)))
        {
            _S2043 = _S2017;
            _S2044 = _S2018;
            _S2045 = 0.0f;
            _S2046 = 0.0f;
        }
        else
        {
            _S2043 = 0.0f;
            _S2044 = 0.0f;
            _S2045 = _S2018;
            _S2046 = _S2017;
        }
        if(_S2024)
        {
            float _S2047 = _S2027 * _S2044;
            float _S2048 = _S2044 + _S2045;
            float _S2049 = _S2026 * _S2044 + _S2043;
            float _S2050 = _S2043 + _S2046;
            float _S2051;
            if(_S2032)
            {
                _S2051 = _S2049;
            }
            else
            {
                DiffPair_float_0 _S2052;
                (&_S2052)->primal_0 = _S2028;
                (&_S2052)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S2052, _S2049);
                _S2051 = 0.90909093618392944f * _S2052.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S2053;
            (&_S2053)->primal_0 = _S1989;
            (&_S2053)->differential_0 = _S2011;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2054;
            (&_S2054)->primal_0 = _S2030;
            (&_S2054)->differential_0 = _S2010;
            s_bwd_prop_interp_0(&_S2053, &_S2054, _S2051);
            float3  _S2055 = _S2054.differential_0 / _S1992;
            float3  _S2056 = _S1991.primal_0 * _S2055;
            float3  _S2057 = _S2031 * _S2055;
            float _S2058 = _S2056.x + _S2056.y + _S2056.z + _S2047;
            DiffPair_float_0 _S2059;
            (&_S2059)->primal_0 = _S2003;
            (&_S2059)->differential_0 = 0.0f;
            DiffPair_float_0 _S2060;
            (&_S2060)->primal_0 = _S2008;
            (&_S2060)->differential_0 = 0.0f;
            DiffPair_float_0 _S2061;
            (&_S2061)->primal_0 = _S2029;
            (&_S2061)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S2059, &_S2060, &_S2061, _S2058);
            float _S2062 = (&_S2053)->differential_0[int(0)] + _S2019[int(0)];
            float _S2063 = (&_S2053)->differential_0[int(1)] + _S2019[int(1)];
            float _S2064 = (&_S2053)->differential_0[int(2)] + _S2019[int(2)];
            float _S2065 = (&_S2053)->differential_0[int(3)] + _S2019[int(3)];
            float _S2066 = (&_S2053)->differential_0[int(4)] + _S2019[int(4)];
            float _S2067 = (&_S2053)->differential_0[int(5)] + _S2019[int(5)];
            float _S2068 = (&_S2053)->differential_0[int(6)] + _S2019[int(6)];
            float _S2069 = (&_S2053)->differential_0[int(7)] + _S2019[int(7)];
            float3  _S2070 = _S2055 + _S2020;
            float3  _S2071 = _S2057 + _S2021;
            float _S2072 = _S2060.differential_0 + _S2022;
            float _S2073 = _S2059.differential_0 + _S2023;
            _S2017 = _S2050;
            _S2018 = _S2048;
            _S2019[int(0)] = _S2062;
            _S2019[int(1)] = _S2063;
            _S2019[int(2)] = _S2064;
            _S2019[int(3)] = _S2065;
            _S2019[int(4)] = _S2066;
            _S2019[int(5)] = _S2067;
            _S2019[int(6)] = _S2068;
            _S2019[int(7)] = _S2069;
            _S2020 = _S2070;
            _S2021 = _S2071;
            _S2022 = _S2072;
            _S2023 = _S2073;
        }
        else
        {
            _S2017 = _S2046;
            _S2018 = _S2045;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S2074;
    (&_S2074)->primal_0 = _S2006;
    (&_S2074)->differential_0 = 0.0f;
    DiffPair_float_0 _S2075;
    (&_S2075)->primal_0 = _S2007;
    (&_S2075)->differential_0 = 0.0f;
    _d_min_0(&_S2074, &_S2075, _S2022);
    DiffPair_float_0 _S2076;
    (&_S2076)->primal_0 = _S2004;
    (&_S2076)->differential_0 = 0.0f;
    DiffPair_float_0 _S2077;
    (&_S2077)->primal_0 = _S2005;
    (&_S2077)->differential_0 = 0.0f;
    _d_min_0(&_S2076, &_S2077, _S2074.differential_0);
    DiffPair_float_0 _S2078;
    (&_S2078)->primal_0 = _S2000;
    (&_S2078)->differential_0 = 0.0f;
    DiffPair_float_0 _S2079;
    (&_S2079)->primal_0 = _S2002;
    (&_S2079)->differential_0 = 0.0f;
    _d_max_0(&_S2078, &_S2079, _S2023);
    DiffPair_float_0 _S2080;
    (&_S2080)->primal_0 = _S2001;
    (&_S2080)->differential_0 = 0.0f;
    DiffPair_float_0 _S2081;
    (&_S2081)->primal_0 = 0.0f;
    (&_S2081)->differential_0 = 0.0f;
    _d_max_0(&_S2080, &_S2081, _S2079.differential_0);
    DiffPair_float_0 _S2082;
    (&_S2082)->primal_0 = _S1998;
    (&_S2082)->differential_0 = 0.0f;
    DiffPair_float_0 _S2083;
    (&_S2083)->primal_0 = _S1999;
    (&_S2083)->differential_0 = 0.0f;
    _d_max_0(&_S2082, &_S2083, _S2078.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S2076.differential_0, _S2077.differential_0, _S2075.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S2082.differential_0, _S2083.differential_0, _S2080.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S2084 = _S1994 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2085;
    (&_S2085)->primal_0 = m_7;
    (&_S2085)->differential_0 = _S2010;
    s_bwd_prop_abs_0(&_S2085, _S2084);
    float3  _S2086 = m_7 * s_diff_n_T_1;
    float3  _S2087 = - ((_S2085.differential_0 + _S1996 * s_diff_n_T_1) / _S1995) + _S2021;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2087;
    float3  _S2088 = _S2086 + _S2020;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2088;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = dpout_rgb_1;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S2019;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S2089, float _S2090, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S2091, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2092, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2093, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2094, float3  _S2095, float _S2096)
{
    FixedArray<float, 8>  _S2097 = _S2091->primal_0;
    float3  _S2098;
    float _S2099;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2100;
    s_primal_ctx_evaluate_color_voxel_0(_S2089, _S2090, &_S2097, (*_S2092).primal_0, (*_S2093).primal_0, (*_S2094).primal_0, &_S2098, &_S2099, &_S2100);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2101 = _S2100;
    s_bwd_prop_evaluate_color_voxel_0(_S2089, _S2090, _S2091, _S2092, _S2093, _S2094, _S2095, _S2096, &_S2101);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  densities_10, float3  rgb_1, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_0, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    FixedArray<float, 8>  _S2102 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = densities_10;
    (&dp_densities_1)->differential_0 = _S2102;
    float3  _S2103 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S2103;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2103;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2103;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_2 = dp_rgb_0.differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

