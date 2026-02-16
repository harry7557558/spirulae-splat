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

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  densities_0, FixedArray<float3 , 16>  sh_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float * depth_0, float3  * rgbs_0)
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((densities_0[int(0)]), (densities_0[int(1)])))), (densities_0[int(2)])))), (densities_0[int(3)])))), (densities_0[int(4)])))), (densities_0[int(5)])))), (densities_0[int(6)])))), (densities_0[int(7)]))) * size_0 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
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
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int((F32_floor((_S186)))), int((F32_floor((_S188)))), int((F32_ceil((_S185)))), int((F32_ceil((_S187)))));
        float _S189 = mean_c_0.z;
        *depth_0 = (F32_max((_S189), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S189 + length_1(mean_c_0)))));
        float3  _S190 = mean_c_0 - - mul_0(transpose_1(R_0), t_0);
        float _S191 = _S190.x;
        float _S192 = _S190.y;
        float _S193 = _S190.z;
        float norm_0 = (F32_sqrt((_S191 * _S191 + _S192 * _S192 + _S193 * _S193)));
        float x_7 = _S191 / norm_0;
        float y_3 = _S192 / norm_0;
        float z_0 = _S193 / norm_0;
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

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  densities_1, FixedArray<float3 , 16>  sh_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float * depth_1, float3  * rgbs_1)
{
    float2  * _S194;
    float2  _S195;
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
    float _S207;
    float2  _S208;
    bool _S209;
    float2  * _S210;
    bool _S211;
    float2  * _S212;
    bool _S213;
    float2  * _S214;
    bool _S215;
    float2  * _S216;
    bool _S217;
    float2  * _S218;
    bool _S219;
    float2  * _S220;
    bool _S221;
    float2  * _S222;
    bool _S223;
    bool _S224;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S225 = mul_0(R_1, pos_1) + t_1;
        pos_c_1[int(0)] = _S225;
        float _S226 = length_1(_S225);
        float _S227 = (F32_min((far_plane_1), (_S226)));
        float _S228 = (F32_max((near_plane_1), (_S226)));
        float3  _S229 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_1;
        pos_c_1[int(1)] = _S229;
        float _S230 = length_1(_S229);
        float _S231 = (F32_min((_S227), (_S230)));
        float _S232 = (F32_max((_S228), (_S230)));
        float3  _S233 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(2)] = _S233;
        float _S234 = length_1(_S233);
        float _S235 = (F32_min((_S231), (_S234)));
        float _S236 = (F32_max((_S232), (_S234)));
        float3  _S237 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(3)] = _S237;
        float _S238 = length_1(_S237);
        float _S239 = (F32_min((_S235), (_S238)));
        float _S240 = (F32_max((_S236), (_S238)));
        float3  _S241 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(4)] = _S241;
        float _S242 = length_1(_S241);
        float _S243 = (F32_min((_S239), (_S242)));
        float _S244 = (F32_max((_S240), (_S242)));
        float3  _S245 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(5)] = _S245;
        float _S246 = length_1(_S245);
        float _S247 = (F32_min((_S243), (_S246)));
        float _S248 = (F32_max((_S244), (_S246)));
        float3  _S249 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_1;
        pos_c_1[int(6)] = _S249;
        float _S250 = length_1(_S249);
        float _S251 = (F32_min((_S247), (_S250)));
        float _S252 = (F32_max((_S248), (_S250)));
        float3  _S253 = mul_0(R_1, pos_1 + make_float3 (size_1)) + t_1;
        pos_c_1[int(7)] = _S253;
        float _S254 = length_1(_S253);
        float _S255 = (F32_min((_S251), (_S254)));
        float _S256 = (F32_max((_S252), (_S254)));
        bool _S257;
        if(_S255 < near_plane_1)
        {
            _S257 = true;
        }
        else
        {
            _S257 = _S256 > far_plane_1;
        }
        if(_S257)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_1 = mul_0(R_1, pos_1 + make_float3 (0.5f * size_1)) + t_1;
        FixedArray<float2 , 8>  uv_1;
        for(;;)
        {
            float k_0;
            float3  _S258 = pos_c_1[int(0)];
            _S194 = &uv_1[int(0)];
            for(;;)
            {
                float2  _S259 = float2 {_S258.x, _S258.y};
                float r_2 = length_0(_S259);
                float _S260 = _S258.z;
                float theta_0 = (F32_atan2((r_2), (_S260)));
                if(theta_0 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S260;
                }
                else
                {
                    k_0 = theta_0 / r_2;
                }
                float2  _S261 = _S259 * make_float2 (k_0);
                uv_1[int(0)] = _S261;
                float2  _S262 = make_float2 (1.0f, 0.0f);
                _S195 = _S262;
                _S196 = dist_coeffs_1[int(0)];
                _S197 = dist_coeffs_1[int(1)];
                _S198 = dist_coeffs_1[int(2)];
                _S199 = dist_coeffs_1[int(3)];
                _S200 = dist_coeffs_1[int(4)];
                _S201 = dist_coeffs_1[int(5)];
                _S202 = dist_coeffs_1[int(6)];
                _S203 = dist_coeffs_1[int(7)];
                _S204 = dist_coeffs_1[int(8)];
                _S205 = dist_coeffs_1[int(9)];
                float u_16 = _S261.x;
                float v_16 = _S261.y;
                float _S263 = u_16 + u_16;
                float r2_16 = u_16 * u_16 + v_16 * v_16;
                float _S264 = dist_coeffs_1[int(2)] + r2_16 * dist_coeffs_1[int(3)];
                float _S265 = dist_coeffs_1[int(1)] + r2_16 * _S264;
                float _S266 = dist_coeffs_1[int(0)] + r2_16 * _S265;
                float _S267 = _S263 * _S266 + (_S263 * _S265 + (_S263 * _S264 + _S263 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16;
                float radial_8 = 1.0f + r2_16 * _S266;
                float _S268 = 2.0f * dist_coeffs_1[int(4)];
                _S206 = _S268;
                float _S269 = _S268 * u_16;
                float _S270 = 2.0f * u_16;
                float s_diff_du_0 = _S268 * v_16 + (_S263 + (_S270 + _S270)) * dist_coeffs_1[int(5)] + _S263 * dist_coeffs_1[int(6)];
                float _S271 = 2.0f * dist_coeffs_1[int(5)];
                _S207 = _S271;
                float _S272 = _S271 * u_16;
                float _S273 = 2.0f * v_16;
                float2  _S274 = _S262 * make_float2 (radial_8) + make_float2 (_S267) * _S261 + make_float2 (s_diff_du_0, _S271 * v_16 + _S263 * dist_coeffs_1[int(4)] + _S263 * dist_coeffs_1[int(7)]);
                float2  _S275 = _S274 + make_float2 (_S274.x * dist_coeffs_1[int(8)] + _S274.y * dist_coeffs_1[int(9)], 0.0f);
                float2  _S276 = make_float2 (0.0f, 1.0f);
                _S208 = _S276;
                float _S277 = v_16 + v_16;
                float2  _S278 = _S276 * make_float2 (radial_8) + make_float2 (_S277 * _S266 + (_S277 * _S265 + (_S277 * _S264 + _S277 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16) * _S261 + make_float2 (_S269 + _S277 * dist_coeffs_1[int(5)] + _S277 * dist_coeffs_1[int(6)], _S272 + (_S277 + (_S273 + _S273)) * dist_coeffs_1[int(4)] + _S277 * dist_coeffs_1[int(7)]);
                Matrix<float, 2, 2>  _S279 = transpose_0(makeMatrix<float, 2, 2> (_S275, _S278 + make_float2 (_S278.x * dist_coeffs_1[int(8)] + _S278.y * dist_coeffs_1[int(9)], 0.0f)));
                bool _S280 = !((F32_min((determinant_0(_S279)), ((F32_min((_S279.rows[int(0)].x), (_S279.rows[int(1)].y)))))) > 0.0f);
                _S209 = _S280;
                if(_S280)
                {
                    break;
                }
                float u_17 = uv_1[int(0)].x;
                float v_17 = uv_1[int(0)].y;
                float r2_17 = u_17 * u_17 + v_17 * v_17;
                float2  _S281 = uv_1[int(0)] * make_float2 (1.0f + r2_17 * (dist_coeffs_1[int(0)] + r2_17 * (dist_coeffs_1[int(1)] + r2_17 * (dist_coeffs_1[int(2)] + r2_17 * dist_coeffs_1[int(3)])))) + make_float2 (_S268 * u_17 * v_17 + dist_coeffs_1[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + dist_coeffs_1[int(6)] * r2_17, _S271 * u_17 * v_17 + dist_coeffs_1[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + dist_coeffs_1[int(7)] * r2_17);
                float2  _S282 = _S281 + make_float2 (dist_coeffs_1[int(8)] * _S281.x + dist_coeffs_1[int(9)] * _S281.y, 0.0f);
                uv_1[int(0)] = make_float2 (fx_1 * _S282.x + cx_1, fy_1 * _S282.y + cy_1);
                break;
            }
            bool all_valid_7 = true & (!_S209);
            float3  _S283 = pos_c_1[int(1)];
            _S210 = &uv_1[int(1)];
            for(;;)
            {
                float2  _S284 = float2 {_S283.x, _S283.y};
                float r_3 = length_0(_S284);
                float _S285 = _S283.z;
                float theta_1 = (F32_atan2((r_3), (_S285)));
                if(theta_1 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S285;
                }
                else
                {
                    k_0 = theta_1 / r_3;
                }
                float2  _S286 = _S284 * make_float2 (k_0);
                uv_1[int(1)] = _S286;
                float u_18 = _S286.x;
                float v_18 = _S286.y;
                float _S287 = u_18 + u_18;
                float r2_18 = u_18 * u_18 + v_18 * v_18;
                float _S288 = _S198 + r2_18 * _S199;
                float _S289 = _S197 + r2_18 * _S288;
                float _S290 = _S196 + r2_18 * _S289;
                float radial_9 = 1.0f + r2_18 * _S290;
                float _S291 = 2.0f * u_18;
                float _S292 = 2.0f * v_18;
                float2  _S293 = _S195 * make_float2 (radial_9) + make_float2 (_S287 * _S290 + (_S287 * _S289 + (_S287 * _S288 + _S287 * _S199 * r2_18) * r2_18) * r2_18) * _S286 + make_float2 (_S206 * v_18 + (_S287 + (_S291 + _S291)) * _S201 + _S287 * _S202, _S207 * v_18 + _S287 * _S200 + _S287 * _S203);
                float _S294 = v_18 + v_18;
                float2  _S295 = _S208 * make_float2 (radial_9) + make_float2 (_S294 * _S290 + (_S294 * _S289 + (_S294 * _S288 + _S294 * _S199 * r2_18) * r2_18) * r2_18) * _S286 + make_float2 (_S206 * u_18 + _S294 * _S201 + _S294 * _S202, _S207 * u_18 + (_S294 + (_S292 + _S292)) * _S200 + _S294 * _S203);
                Matrix<float, 2, 2>  _S296 = transpose_0(makeMatrix<float, 2, 2> (_S293 + make_float2 (_S293.x * _S204 + _S293.y * _S205, 0.0f), _S295 + make_float2 (_S295.x * _S204 + _S295.y * _S205, 0.0f)));
                bool _S297 = !((F32_min((determinant_0(_S296)), ((F32_min((_S296.rows[int(0)].x), (_S296.rows[int(1)].y)))))) > 0.0f);
                _S211 = _S297;
                if(_S297)
                {
                    break;
                }
                float u_19 = uv_1[int(1)].x;
                float v_19 = uv_1[int(1)].y;
                float r2_19 = u_19 * u_19 + v_19 * v_19;
                float2  _S298 = uv_1[int(1)] * make_float2 (1.0f + r2_19 * (_S196 + r2_19 * (_S197 + r2_19 * (_S198 + r2_19 * _S199)))) + make_float2 (_S206 * u_19 * v_19 + _S201 * (r2_19 + 2.0f * u_19 * u_19) + _S202 * r2_19, _S207 * u_19 * v_19 + _S200 * (r2_19 + 2.0f * v_19 * v_19) + _S203 * r2_19);
                float2  _S299 = _S298 + make_float2 (_S204 * _S298.x + _S205 * _S298.y, 0.0f);
                uv_1[int(1)] = make_float2 (fx_1 * _S299.x + cx_1, fy_1 * _S299.y + cy_1);
                break;
            }
            bool all_valid_8 = all_valid_7 & (!_S211);
            float3  _S300 = pos_c_1[int(2)];
            _S212 = &uv_1[int(2)];
            for(;;)
            {
                float2  _S301 = float2 {_S300.x, _S300.y};
                float r_4 = length_0(_S301);
                float _S302 = _S300.z;
                float theta_2 = (F32_atan2((r_4), (_S302)));
                if(theta_2 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S302;
                }
                else
                {
                    k_0 = theta_2 / r_4;
                }
                float2  _S303 = _S301 * make_float2 (k_0);
                uv_1[int(2)] = _S303;
                float u_20 = _S303.x;
                float v_20 = _S303.y;
                float _S304 = u_20 + u_20;
                float r2_20 = u_20 * u_20 + v_20 * v_20;
                float _S305 = _S198 + r2_20 * _S199;
                float _S306 = _S197 + r2_20 * _S305;
                float _S307 = _S196 + r2_20 * _S306;
                float radial_10 = 1.0f + r2_20 * _S307;
                float _S308 = 2.0f * u_20;
                float _S309 = 2.0f * v_20;
                float2  _S310 = _S195 * make_float2 (radial_10) + make_float2 (_S304 * _S307 + (_S304 * _S306 + (_S304 * _S305 + _S304 * _S199 * r2_20) * r2_20) * r2_20) * _S303 + make_float2 (_S206 * v_20 + (_S304 + (_S308 + _S308)) * _S201 + _S304 * _S202, _S207 * v_20 + _S304 * _S200 + _S304 * _S203);
                float _S311 = v_20 + v_20;
                float2  _S312 = _S208 * make_float2 (radial_10) + make_float2 (_S311 * _S307 + (_S311 * _S306 + (_S311 * _S305 + _S311 * _S199 * r2_20) * r2_20) * r2_20) * _S303 + make_float2 (_S206 * u_20 + _S311 * _S201 + _S311 * _S202, _S207 * u_20 + (_S311 + (_S309 + _S309)) * _S200 + _S311 * _S203);
                Matrix<float, 2, 2>  _S313 = transpose_0(makeMatrix<float, 2, 2> (_S310 + make_float2 (_S310.x * _S204 + _S310.y * _S205, 0.0f), _S312 + make_float2 (_S312.x * _S204 + _S312.y * _S205, 0.0f)));
                bool _S314 = !((F32_min((determinant_0(_S313)), ((F32_min((_S313.rows[int(0)].x), (_S313.rows[int(1)].y)))))) > 0.0f);
                _S213 = _S314;
                if(_S314)
                {
                    break;
                }
                float u_21 = uv_1[int(2)].x;
                float v_21 = uv_1[int(2)].y;
                float r2_21 = u_21 * u_21 + v_21 * v_21;
                float2  _S315 = uv_1[int(2)] * make_float2 (1.0f + r2_21 * (_S196 + r2_21 * (_S197 + r2_21 * (_S198 + r2_21 * _S199)))) + make_float2 (_S206 * u_21 * v_21 + _S201 * (r2_21 + 2.0f * u_21 * u_21) + _S202 * r2_21, _S207 * u_21 * v_21 + _S200 * (r2_21 + 2.0f * v_21 * v_21) + _S203 * r2_21);
                float2  _S316 = _S315 + make_float2 (_S204 * _S315.x + _S205 * _S315.y, 0.0f);
                uv_1[int(2)] = make_float2 (fx_1 * _S316.x + cx_1, fy_1 * _S316.y + cy_1);
                break;
            }
            bool all_valid_9 = all_valid_8 & (!_S213);
            float3  _S317 = pos_c_1[int(3)];
            _S214 = &uv_1[int(3)];
            for(;;)
            {
                float2  _S318 = float2 {_S317.x, _S317.y};
                float r_5 = length_0(_S318);
                float _S319 = _S317.z;
                float theta_3 = (F32_atan2((r_5), (_S319)));
                if(theta_3 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S319;
                }
                else
                {
                    k_0 = theta_3 / r_5;
                }
                float2  _S320 = _S318 * make_float2 (k_0);
                uv_1[int(3)] = _S320;
                float u_22 = _S320.x;
                float v_22 = _S320.y;
                float _S321 = u_22 + u_22;
                float r2_22 = u_22 * u_22 + v_22 * v_22;
                float _S322 = _S198 + r2_22 * _S199;
                float _S323 = _S197 + r2_22 * _S322;
                float _S324 = _S196 + r2_22 * _S323;
                float radial_11 = 1.0f + r2_22 * _S324;
                float _S325 = 2.0f * u_22;
                float _S326 = 2.0f * v_22;
                float2  _S327 = _S195 * make_float2 (radial_11) + make_float2 (_S321 * _S324 + (_S321 * _S323 + (_S321 * _S322 + _S321 * _S199 * r2_22) * r2_22) * r2_22) * _S320 + make_float2 (_S206 * v_22 + (_S321 + (_S325 + _S325)) * _S201 + _S321 * _S202, _S207 * v_22 + _S321 * _S200 + _S321 * _S203);
                float _S328 = v_22 + v_22;
                float2  _S329 = _S208 * make_float2 (radial_11) + make_float2 (_S328 * _S324 + (_S328 * _S323 + (_S328 * _S322 + _S328 * _S199 * r2_22) * r2_22) * r2_22) * _S320 + make_float2 (_S206 * u_22 + _S328 * _S201 + _S328 * _S202, _S207 * u_22 + (_S328 + (_S326 + _S326)) * _S200 + _S328 * _S203);
                Matrix<float, 2, 2>  _S330 = transpose_0(makeMatrix<float, 2, 2> (_S327 + make_float2 (_S327.x * _S204 + _S327.y * _S205, 0.0f), _S329 + make_float2 (_S329.x * _S204 + _S329.y * _S205, 0.0f)));
                bool _S331 = !((F32_min((determinant_0(_S330)), ((F32_min((_S330.rows[int(0)].x), (_S330.rows[int(1)].y)))))) > 0.0f);
                _S215 = _S331;
                if(_S331)
                {
                    break;
                }
                float u_23 = uv_1[int(3)].x;
                float v_23 = uv_1[int(3)].y;
                float r2_23 = u_23 * u_23 + v_23 * v_23;
                float2  _S332 = uv_1[int(3)] * make_float2 (1.0f + r2_23 * (_S196 + r2_23 * (_S197 + r2_23 * (_S198 + r2_23 * _S199)))) + make_float2 (_S206 * u_23 * v_23 + _S201 * (r2_23 + 2.0f * u_23 * u_23) + _S202 * r2_23, _S207 * u_23 * v_23 + _S200 * (r2_23 + 2.0f * v_23 * v_23) + _S203 * r2_23);
                float2  _S333 = _S332 + make_float2 (_S204 * _S332.x + _S205 * _S332.y, 0.0f);
                uv_1[int(3)] = make_float2 (fx_1 * _S333.x + cx_1, fy_1 * _S333.y + cy_1);
                break;
            }
            bool all_valid_10 = all_valid_9 & (!_S215);
            float3  _S334 = pos_c_1[int(4)];
            _S216 = &uv_1[int(4)];
            for(;;)
            {
                float2  _S335 = float2 {_S334.x, _S334.y};
                float r_6 = length_0(_S335);
                float _S336 = _S334.z;
                float theta_4 = (F32_atan2((r_6), (_S336)));
                if(theta_4 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S336;
                }
                else
                {
                    k_0 = theta_4 / r_6;
                }
                float2  _S337 = _S335 * make_float2 (k_0);
                uv_1[int(4)] = _S337;
                float u_24 = _S337.x;
                float v_24 = _S337.y;
                float _S338 = u_24 + u_24;
                float r2_24 = u_24 * u_24 + v_24 * v_24;
                float _S339 = _S198 + r2_24 * _S199;
                float _S340 = _S197 + r2_24 * _S339;
                float _S341 = _S196 + r2_24 * _S340;
                float radial_12 = 1.0f + r2_24 * _S341;
                float _S342 = 2.0f * u_24;
                float _S343 = 2.0f * v_24;
                float2  _S344 = _S195 * make_float2 (radial_12) + make_float2 (_S338 * _S341 + (_S338 * _S340 + (_S338 * _S339 + _S338 * _S199 * r2_24) * r2_24) * r2_24) * _S337 + make_float2 (_S206 * v_24 + (_S338 + (_S342 + _S342)) * _S201 + _S338 * _S202, _S207 * v_24 + _S338 * _S200 + _S338 * _S203);
                float _S345 = v_24 + v_24;
                float2  _S346 = _S208 * make_float2 (radial_12) + make_float2 (_S345 * _S341 + (_S345 * _S340 + (_S345 * _S339 + _S345 * _S199 * r2_24) * r2_24) * r2_24) * _S337 + make_float2 (_S206 * u_24 + _S345 * _S201 + _S345 * _S202, _S207 * u_24 + (_S345 + (_S343 + _S343)) * _S200 + _S345 * _S203);
                Matrix<float, 2, 2>  _S347 = transpose_0(makeMatrix<float, 2, 2> (_S344 + make_float2 (_S344.x * _S204 + _S344.y * _S205, 0.0f), _S346 + make_float2 (_S346.x * _S204 + _S346.y * _S205, 0.0f)));
                bool _S348 = !((F32_min((determinant_0(_S347)), ((F32_min((_S347.rows[int(0)].x), (_S347.rows[int(1)].y)))))) > 0.0f);
                _S217 = _S348;
                if(_S348)
                {
                    break;
                }
                float u_25 = uv_1[int(4)].x;
                float v_25 = uv_1[int(4)].y;
                float r2_25 = u_25 * u_25 + v_25 * v_25;
                float2  _S349 = uv_1[int(4)] * make_float2 (1.0f + r2_25 * (_S196 + r2_25 * (_S197 + r2_25 * (_S198 + r2_25 * _S199)))) + make_float2 (_S206 * u_25 * v_25 + _S201 * (r2_25 + 2.0f * u_25 * u_25) + _S202 * r2_25, _S207 * u_25 * v_25 + _S200 * (r2_25 + 2.0f * v_25 * v_25) + _S203 * r2_25);
                float2  _S350 = _S349 + make_float2 (_S204 * _S349.x + _S205 * _S349.y, 0.0f);
                uv_1[int(4)] = make_float2 (fx_1 * _S350.x + cx_1, fy_1 * _S350.y + cy_1);
                break;
            }
            bool all_valid_11 = all_valid_10 & (!_S217);
            float3  _S351 = pos_c_1[int(5)];
            _S218 = &uv_1[int(5)];
            for(;;)
            {
                float2  _S352 = float2 {_S351.x, _S351.y};
                float r_7 = length_0(_S352);
                float _S353 = _S351.z;
                float theta_5 = (F32_atan2((r_7), (_S353)));
                if(theta_5 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_5 * theta_5 / 3.0f) / _S353;
                }
                else
                {
                    k_0 = theta_5 / r_7;
                }
                float2  _S354 = _S352 * make_float2 (k_0);
                uv_1[int(5)] = _S354;
                float u_26 = _S354.x;
                float v_26 = _S354.y;
                float _S355 = u_26 + u_26;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float _S356 = _S198 + r2_26 * _S199;
                float _S357 = _S197 + r2_26 * _S356;
                float _S358 = _S196 + r2_26 * _S357;
                float radial_13 = 1.0f + r2_26 * _S358;
                float _S359 = 2.0f * u_26;
                float _S360 = 2.0f * v_26;
                float2  _S361 = _S195 * make_float2 (radial_13) + make_float2 (_S355 * _S358 + (_S355 * _S357 + (_S355 * _S356 + _S355 * _S199 * r2_26) * r2_26) * r2_26) * _S354 + make_float2 (_S206 * v_26 + (_S355 + (_S359 + _S359)) * _S201 + _S355 * _S202, _S207 * v_26 + _S355 * _S200 + _S355 * _S203);
                float _S362 = v_26 + v_26;
                float2  _S363 = _S208 * make_float2 (radial_13) + make_float2 (_S362 * _S358 + (_S362 * _S357 + (_S362 * _S356 + _S362 * _S199 * r2_26) * r2_26) * r2_26) * _S354 + make_float2 (_S206 * u_26 + _S362 * _S201 + _S362 * _S202, _S207 * u_26 + (_S362 + (_S360 + _S360)) * _S200 + _S362 * _S203);
                Matrix<float, 2, 2>  _S364 = transpose_0(makeMatrix<float, 2, 2> (_S361 + make_float2 (_S361.x * _S204 + _S361.y * _S205, 0.0f), _S363 + make_float2 (_S363.x * _S204 + _S363.y * _S205, 0.0f)));
                bool _S365 = !((F32_min((determinant_0(_S364)), ((F32_min((_S364.rows[int(0)].x), (_S364.rows[int(1)].y)))))) > 0.0f);
                _S219 = _S365;
                if(_S365)
                {
                    break;
                }
                float u_27 = uv_1[int(5)].x;
                float v_27 = uv_1[int(5)].y;
                float r2_27 = u_27 * u_27 + v_27 * v_27;
                float2  _S366 = uv_1[int(5)] * make_float2 (1.0f + r2_27 * (_S196 + r2_27 * (_S197 + r2_27 * (_S198 + r2_27 * _S199)))) + make_float2 (_S206 * u_27 * v_27 + _S201 * (r2_27 + 2.0f * u_27 * u_27) + _S202 * r2_27, _S207 * u_27 * v_27 + _S200 * (r2_27 + 2.0f * v_27 * v_27) + _S203 * r2_27);
                float2  _S367 = _S366 + make_float2 (_S204 * _S366.x + _S205 * _S366.y, 0.0f);
                uv_1[int(5)] = make_float2 (fx_1 * _S367.x + cx_1, fy_1 * _S367.y + cy_1);
                break;
            }
            bool all_valid_12 = all_valid_11 & (!_S219);
            float3  _S368 = pos_c_1[int(6)];
            _S220 = &uv_1[int(6)];
            for(;;)
            {
                float2  _S369 = float2 {_S368.x, _S368.y};
                float r_8 = length_0(_S369);
                float _S370 = _S368.z;
                float theta_6 = (F32_atan2((r_8), (_S370)));
                if(theta_6 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_6 * theta_6 / 3.0f) / _S370;
                }
                else
                {
                    k_0 = theta_6 / r_8;
                }
                float2  _S371 = _S369 * make_float2 (k_0);
                uv_1[int(6)] = _S371;
                float u_28 = _S371.x;
                float v_28 = _S371.y;
                float _S372 = u_28 + u_28;
                float r2_28 = u_28 * u_28 + v_28 * v_28;
                float _S373 = _S198 + r2_28 * _S199;
                float _S374 = _S197 + r2_28 * _S373;
                float _S375 = _S196 + r2_28 * _S374;
                float radial_14 = 1.0f + r2_28 * _S375;
                float _S376 = 2.0f * u_28;
                float _S377 = 2.0f * v_28;
                float2  _S378 = _S195 * make_float2 (radial_14) + make_float2 (_S372 * _S375 + (_S372 * _S374 + (_S372 * _S373 + _S372 * _S199 * r2_28) * r2_28) * r2_28) * _S371 + make_float2 (_S206 * v_28 + (_S372 + (_S376 + _S376)) * _S201 + _S372 * _S202, _S207 * v_28 + _S372 * _S200 + _S372 * _S203);
                float _S379 = v_28 + v_28;
                float2  _S380 = _S208 * make_float2 (radial_14) + make_float2 (_S379 * _S375 + (_S379 * _S374 + (_S379 * _S373 + _S379 * _S199 * r2_28) * r2_28) * r2_28) * _S371 + make_float2 (_S206 * u_28 + _S379 * _S201 + _S379 * _S202, _S207 * u_28 + (_S379 + (_S377 + _S377)) * _S200 + _S379 * _S203);
                Matrix<float, 2, 2>  _S381 = transpose_0(makeMatrix<float, 2, 2> (_S378 + make_float2 (_S378.x * _S204 + _S378.y * _S205, 0.0f), _S380 + make_float2 (_S380.x * _S204 + _S380.y * _S205, 0.0f)));
                bool _S382 = !((F32_min((determinant_0(_S381)), ((F32_min((_S381.rows[int(0)].x), (_S381.rows[int(1)].y)))))) > 0.0f);
                _S221 = _S382;
                if(_S382)
                {
                    break;
                }
                float u_29 = uv_1[int(6)].x;
                float v_29 = uv_1[int(6)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S383 = uv_1[int(6)] * make_float2 (1.0f + r2_29 * (_S196 + r2_29 * (_S197 + r2_29 * (_S198 + r2_29 * _S199)))) + make_float2 (_S206 * u_29 * v_29 + _S201 * (r2_29 + 2.0f * u_29 * u_29) + _S202 * r2_29, _S207 * u_29 * v_29 + _S200 * (r2_29 + 2.0f * v_29 * v_29) + _S203 * r2_29);
                float2  _S384 = _S383 + make_float2 (_S204 * _S383.x + _S205 * _S383.y, 0.0f);
                uv_1[int(6)] = make_float2 (fx_1 * _S384.x + cx_1, fy_1 * _S384.y + cy_1);
                break;
            }
            bool all_valid_13 = all_valid_12 & (!_S221);
            float3  _S385 = pos_c_1[int(7)];
            _S222 = &uv_1[int(7)];
            for(;;)
            {
                float2  _S386 = float2 {_S385.x, _S385.y};
                float r_9 = length_0(_S386);
                float _S387 = _S385.z;
                float theta_7 = (F32_atan2((r_9), (_S387)));
                if(theta_7 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_7 * theta_7 / 3.0f) / _S387;
                }
                else
                {
                    k_0 = theta_7 / r_9;
                }
                float2  _S388 = _S386 * make_float2 (k_0);
                uv_1[int(7)] = _S388;
                float u_30 = _S388.x;
                float v_30 = _S388.y;
                float _S389 = u_30 + u_30;
                float r2_30 = u_30 * u_30 + v_30 * v_30;
                float _S390 = _S198 + r2_30 * _S199;
                float _S391 = _S197 + r2_30 * _S390;
                float _S392 = _S196 + r2_30 * _S391;
                float radial_15 = 1.0f + r2_30 * _S392;
                float _S393 = 2.0f * u_30;
                float _S394 = 2.0f * v_30;
                float2  _S395 = _S195 * make_float2 (radial_15) + make_float2 (_S389 * _S392 + (_S389 * _S391 + (_S389 * _S390 + _S389 * _S199 * r2_30) * r2_30) * r2_30) * _S388 + make_float2 (_S206 * v_30 + (_S389 + (_S393 + _S393)) * _S201 + _S389 * _S202, _S207 * v_30 + _S389 * _S200 + _S389 * _S203);
                float _S396 = v_30 + v_30;
                float2  _S397 = _S208 * make_float2 (radial_15) + make_float2 (_S396 * _S392 + (_S396 * _S391 + (_S396 * _S390 + _S396 * _S199 * r2_30) * r2_30) * r2_30) * _S388 + make_float2 (_S206 * u_30 + _S396 * _S201 + _S396 * _S202, _S207 * u_30 + (_S396 + (_S394 + _S394)) * _S200 + _S396 * _S203);
                Matrix<float, 2, 2>  _S398 = transpose_0(makeMatrix<float, 2, 2> (_S395 + make_float2 (_S395.x * _S204 + _S395.y * _S205, 0.0f), _S397 + make_float2 (_S397.x * _S204 + _S397.y * _S205, 0.0f)));
                bool _S399 = !((F32_min((determinant_0(_S398)), ((F32_min((_S398.rows[int(0)].x), (_S398.rows[int(1)].y)))))) > 0.0f);
                _S223 = _S399;
                if(_S399)
                {
                    break;
                }
                float u_31 = uv_1[int(7)].x;
                float v_31 = uv_1[int(7)].y;
                float r2_31 = u_31 * u_31 + v_31 * v_31;
                float2  _S400 = uv_1[int(7)] * make_float2 (1.0f + r2_31 * (_S196 + r2_31 * (_S197 + r2_31 * (_S198 + r2_31 * _S199)))) + make_float2 (_S206 * u_31 * v_31 + _S201 * (r2_31 + 2.0f * u_31 * u_31) + _S202 * r2_31, _S207 * u_31 * v_31 + _S200 * (r2_31 + 2.0f * v_31 * v_31) + _S203 * r2_31);
                float2  _S401 = _S400 + make_float2 (_S204 * _S400.x + _S205 * _S400.y, 0.0f);
                uv_1[int(7)] = make_float2 (fx_1 * _S401.x + cx_1, fy_1 * _S401.y + cy_1);
                break;
            }
            _S224 = all_valid_13 & (!_S223);
            break;
        }
        if(!_S224)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((densities_1[int(0)]), (densities_1[int(1)])))), (densities_1[int(2)])))), (densities_1[int(3)])))), (densities_1[int(4)])))), (densities_1[int(5)])))), (densities_1[int(6)])))), (densities_1[int(7)]))) * size_1 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S402 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S194).x), ((*_S210).x)))), ((*_S212).x)))), ((*_S214).x)))), ((*_S216).x)))), ((*_S218).x)))), ((*_S220).x)))), ((*_S222).x)));
        float _S403 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S194).x), ((*_S210).x)))), ((*_S212).x)))), ((*_S214).x)))), ((*_S216).x)))), ((*_S218).x)))), ((*_S220).x)))), ((*_S222).x)));
        float _S404 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S194).y), ((*_S210).y)))), ((*_S212).y)))), ((*_S214).y)))), ((*_S216).y)))), ((*_S218).y)))), ((*_S220).y)))), ((*_S222).y)));
        float _S405 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S194).y), ((*_S210).y)))), ((*_S212).y)))), ((*_S214).y)))), ((*_S216).y)))), ((*_S218).y)))), ((*_S220).y)))), ((*_S222).y)));
        if(_S402 <= 0.0f)
        {
            _S257 = true;
        }
        else
        {
            _S257 = _S403 >= float(image_width_1);
        }
        if(_S257)
        {
            _S257 = true;
        }
        else
        {
            _S257 = _S404 <= 0.0f;
        }
        if(_S257)
        {
            _S257 = true;
        }
        else
        {
            _S257 = _S405 >= float(image_height_1);
        }
        if(_S257)
        {
            _S257 = true;
        }
        else
        {
            if(_S255 <= 0.0f)
            {
                if(_S403 <= 0.0f)
                {
                    _S257 = _S402 >= float(image_width_1);
                }
                else
                {
                    _S257 = false;
                }
                if(_S257)
                {
                    _S257 = true;
                }
                else
                {
                    if(_S405 <= 0.0f)
                    {
                        _S257 = _S404 >= float(image_width_1);
                    }
                    else
                    {
                        _S257 = false;
                    }
                }
            }
            else
            {
                _S257 = false;
            }
        }
        if(_S257)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int((F32_floor((_S403)))), int((F32_floor((_S405)))), int((F32_ceil((_S402)))), int((F32_ceil((_S404)))));
        float _S406 = mean_c_1.z;
        *depth_1 = (F32_max((_S406), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S406 + length_1(mean_c_1)))));
        float3  _S407 = mean_c_1 - - mul_0(transpose_1(R_1), t_1);
        float _S408 = _S407.x;
        float _S409 = _S407.y;
        float _S410 = _S407.z;
        float norm_1 = (F32_sqrt((_S408 * _S408 + _S409 * _S409 + _S410 * _S410)));
        float x_8 = _S408 / norm_1;
        float y_4 = _S409 / norm_1;
        float z_1 = _S410 / norm_1;
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

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  densities_2, FixedArray<float3 , 16>  sh_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float * depth_2, float3  * rgbs_2)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S411 = mul_0(R_2, pos_2) + t_2;
    pos_c_2[int(0)] = _S411;
    float3  _S412 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_2;
    pos_c_2[int(1)] = _S412;
    float3  _S413 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(2)] = _S413;
    float3  _S414 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(3)] = _S414;
    float3  _S415 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(4)] = _S415;
    float3  _S416 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(5)] = _S416;
    float3  _S417 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_2;
    pos_c_2[int(6)] = _S417;
    float3  _S418 = mul_0(R_2, pos_2 + make_float3 (size_2)) + t_2;
    pos_c_2[int(7)] = _S418;
    float3  mean_c_2 = mul_0(R_2, pos_2 + make_float3 (0.5f * size_2)) + t_2;
    FixedArray<float2 , 8>  uv_2;
    float2  _S419 = float2 {_S411.x, _S411.y} / make_float2 (_S411.z);
    float u_32 = _S419.x;
    float v_32 = _S419.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float _S420 = 2.0f * dist_coeffs_2[int(4)];
    float _S421 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S422 = _S419 * make_float2 (1.0f + r2_32 * (dist_coeffs_2[int(0)] + r2_32 * (dist_coeffs_2[int(1)] + r2_32 * (dist_coeffs_2[int(2)] + r2_32 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_32 * v_32 + dist_coeffs_2[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + dist_coeffs_2[int(6)] * r2_32, _S421 * u_32 * v_32 + dist_coeffs_2[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + dist_coeffs_2[int(7)] * r2_32);
    float2  _S423 = _S422 + make_float2 (dist_coeffs_2[int(8)] * _S422.x + dist_coeffs_2[int(9)] * _S422.y, 0.0f);
    float _S424 = fx_2 * _S423.x + cx_2;
    float _S425 = fy_2 * _S423.y + cy_2;
    uv_2[int(0)] = make_float2 (_S424, _S425);
    float2  _S426 = float2 {_S412.x, _S412.y} / make_float2 (_S412.z);
    float u_33 = _S426.x;
    float v_33 = _S426.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S427 = _S426 * make_float2 (1.0f + r2_33 * (dist_coeffs_2[int(0)] + r2_33 * (dist_coeffs_2[int(1)] + r2_33 * (dist_coeffs_2[int(2)] + r2_33 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_33 * v_33 + dist_coeffs_2[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + dist_coeffs_2[int(6)] * r2_33, _S421 * u_33 * v_33 + dist_coeffs_2[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + dist_coeffs_2[int(7)] * r2_33);
    float2  _S428 = _S427 + make_float2 (dist_coeffs_2[int(8)] * _S427.x + dist_coeffs_2[int(9)] * _S427.y, 0.0f);
    float _S429 = fx_2 * _S428.x + cx_2;
    float _S430 = fy_2 * _S428.y + cy_2;
    uv_2[int(1)] = make_float2 (_S429, _S430);
    float2  _S431 = float2 {_S413.x, _S413.y} / make_float2 (_S413.z);
    float u_34 = _S431.x;
    float v_34 = _S431.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S432 = _S431 * make_float2 (1.0f + r2_34 * (dist_coeffs_2[int(0)] + r2_34 * (dist_coeffs_2[int(1)] + r2_34 * (dist_coeffs_2[int(2)] + r2_34 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_34 * v_34 + dist_coeffs_2[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + dist_coeffs_2[int(6)] * r2_34, _S421 * u_34 * v_34 + dist_coeffs_2[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + dist_coeffs_2[int(7)] * r2_34);
    float2  _S433 = _S432 + make_float2 (dist_coeffs_2[int(8)] * _S432.x + dist_coeffs_2[int(9)] * _S432.y, 0.0f);
    float _S434 = fx_2 * _S433.x + cx_2;
    float _S435 = fy_2 * _S433.y + cy_2;
    uv_2[int(2)] = make_float2 (_S434, _S435);
    float2  _S436 = float2 {_S414.x, _S414.y} / make_float2 (_S414.z);
    float u_35 = _S436.x;
    float v_35 = _S436.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S437 = _S436 * make_float2 (1.0f + r2_35 * (dist_coeffs_2[int(0)] + r2_35 * (dist_coeffs_2[int(1)] + r2_35 * (dist_coeffs_2[int(2)] + r2_35 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_35 * v_35 + dist_coeffs_2[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + dist_coeffs_2[int(6)] * r2_35, _S421 * u_35 * v_35 + dist_coeffs_2[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + dist_coeffs_2[int(7)] * r2_35);
    float2  _S438 = _S437 + make_float2 (dist_coeffs_2[int(8)] * _S437.x + dist_coeffs_2[int(9)] * _S437.y, 0.0f);
    float _S439 = fx_2 * _S438.x + cx_2;
    float _S440 = fy_2 * _S438.y + cy_2;
    uv_2[int(3)] = make_float2 (_S439, _S440);
    float2  _S441 = float2 {_S415.x, _S415.y} / make_float2 (_S415.z);
    float u_36 = _S441.x;
    float v_36 = _S441.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S442 = _S441 * make_float2 (1.0f + r2_36 * (dist_coeffs_2[int(0)] + r2_36 * (dist_coeffs_2[int(1)] + r2_36 * (dist_coeffs_2[int(2)] + r2_36 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_36 * v_36 + dist_coeffs_2[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + dist_coeffs_2[int(6)] * r2_36, _S421 * u_36 * v_36 + dist_coeffs_2[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + dist_coeffs_2[int(7)] * r2_36);
    float2  _S443 = _S442 + make_float2 (dist_coeffs_2[int(8)] * _S442.x + dist_coeffs_2[int(9)] * _S442.y, 0.0f);
    float _S444 = fx_2 * _S443.x + cx_2;
    float _S445 = fy_2 * _S443.y + cy_2;
    uv_2[int(4)] = make_float2 (_S444, _S445);
    float2  _S446 = float2 {_S416.x, _S416.y} / make_float2 (_S416.z);
    float u_37 = _S446.x;
    float v_37 = _S446.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S447 = _S446 * make_float2 (1.0f + r2_37 * (dist_coeffs_2[int(0)] + r2_37 * (dist_coeffs_2[int(1)] + r2_37 * (dist_coeffs_2[int(2)] + r2_37 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_37 * v_37 + dist_coeffs_2[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + dist_coeffs_2[int(6)] * r2_37, _S421 * u_37 * v_37 + dist_coeffs_2[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + dist_coeffs_2[int(7)] * r2_37);
    float2  _S448 = _S447 + make_float2 (dist_coeffs_2[int(8)] * _S447.x + dist_coeffs_2[int(9)] * _S447.y, 0.0f);
    float _S449 = fx_2 * _S448.x + cx_2;
    float _S450 = fy_2 * _S448.y + cy_2;
    uv_2[int(5)] = make_float2 (_S449, _S450);
    float2  _S451 = float2 {_S417.x, _S417.y} / make_float2 (_S417.z);
    float u_38 = _S451.x;
    float v_38 = _S451.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float2  _S452 = _S451 * make_float2 (1.0f + r2_38 * (dist_coeffs_2[int(0)] + r2_38 * (dist_coeffs_2[int(1)] + r2_38 * (dist_coeffs_2[int(2)] + r2_38 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_38 * v_38 + dist_coeffs_2[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + dist_coeffs_2[int(6)] * r2_38, _S421 * u_38 * v_38 + dist_coeffs_2[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + dist_coeffs_2[int(7)] * r2_38);
    float2  _S453 = _S452 + make_float2 (dist_coeffs_2[int(8)] * _S452.x + dist_coeffs_2[int(9)] * _S452.y, 0.0f);
    float _S454 = fx_2 * _S453.x + cx_2;
    float _S455 = fy_2 * _S453.y + cy_2;
    uv_2[int(6)] = make_float2 (_S454, _S455);
    float2  _S456 = float2 {_S418.x, _S418.y} / make_float2 (_S418.z);
    float u_39 = _S456.x;
    float v_39 = _S456.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float2  _S457 = _S456 * make_float2 (1.0f + r2_39 * (dist_coeffs_2[int(0)] + r2_39 * (dist_coeffs_2[int(1)] + r2_39 * (dist_coeffs_2[int(2)] + r2_39 * dist_coeffs_2[int(3)])))) + make_float2 (_S420 * u_39 * v_39 + dist_coeffs_2[int(5)] * (r2_39 + 2.0f * u_39 * u_39) + dist_coeffs_2[int(6)] * r2_39, _S421 * u_39 * v_39 + dist_coeffs_2[int(4)] * (r2_39 + 2.0f * v_39 * v_39) + dist_coeffs_2[int(7)] * r2_39);
    float2  _S458 = _S457 + make_float2 (dist_coeffs_2[int(8)] * _S457.x + dist_coeffs_2[int(9)] * _S457.y, 0.0f);
    float _S459 = fx_2 * _S458.x + cx_2;
    float _S460 = fy_2 * _S458.y + cy_2;
    uv_2[int(7)] = make_float2 (_S459, _S460);
    *aabb_xyxy_2 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S424), (_S429)))), (_S434)))), (_S439)))), (_S444)))), (_S449)))), (_S454)))), (_S459))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S425), (_S430)))), (_S435)))), (_S440)))), (_S445)))), (_S450)))), (_S455)))), (_S460))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S424), (_S429)))), (_S434)))), (_S439)))), (_S444)))), (_S449)))), (_S454)))), (_S459))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S425), (_S430)))), (_S435)))), (_S440)))), (_S445)))), (_S450)))), (_S455)))), (_S460))))))));
    float _S461 = mean_c_2.z;
    *depth_2 = (F32_max((_S461), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S461 + length_1(mean_c_2)))));
    float3  _S462 = mean_c_2 - - mul_0(transpose_1(R_2), t_2);
    float _S463 = _S462.x;
    float _S464 = _S462.y;
    float _S465 = _S462.z;
    float norm_2 = (F32_sqrt((_S463 * _S463 + _S464 * _S464 + _S465 * _S465)));
    float x_9 = _S463 / norm_2;
    float y_5 = _S464 / norm_2;
    float z_2 = _S465 / norm_2;
    float z2_2 = z_2 * z_2;
    float fTmp0B_2 = -1.09254848957061768f * z_2;
    float fC1_2 = x_9 * x_9 - y_5 * y_5;
    float fS1_2 = 2.0f * x_9 * y_5;
    float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_2;
    *rgbs_2 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_2[int(1)] + make_float3 (z_2) * sh_coeffs_2[int(2)] - make_float3 (x_9) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_5) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_2 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_9) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_2 + y_5 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_5) * sh_coeffs_2[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_2 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_9) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_2 - y_5 * fS1_2)) * sh_coeffs_2[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  densities_3, FixedArray<float3 , 16>  sh_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float * depth_3, float3  * rgbs_3)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S466 = mul_0(R_3, pos_3) + t_3;
    pos_c_3[int(0)] = _S466;
    pos_c_3[int(1)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_3;
    pos_c_3[int(2)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(3)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(4)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(5)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(6)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_3;
    pos_c_3[int(7)] = mul_0(R_3, pos_3 + make_float3 (size_3)) + t_3;
    float3  mean_c_3 = mul_0(R_3, pos_3 + make_float3 (0.5f * size_3)) + t_3;
    FixedArray<float2 , 8>  uv_3;
    float2  _S467 = float2 {_S466.x, _S466.y};
    float r_10 = length_0(_S467);
    float _S468 = _S466.z;
    float theta_8 = (F32_atan2((r_10), (_S468)));
    float k_1;
    if(theta_8 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_8 * theta_8 / 3.0f) / _S468;
    }
    else
    {
        k_1 = theta_8 / r_10;
    }
    float2  _S469 = _S467 * make_float2 (k_1);
    float u_40 = _S469.x;
    float v_40 = _S469.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S470 = 2.0f * dist_coeffs_3[int(4)];
    float _S471 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S472 = _S469 * make_float2 (1.0f + r2_40 * (dist_coeffs_3[int(0)] + r2_40 * (dist_coeffs_3[int(1)] + r2_40 * (dist_coeffs_3[int(2)] + r2_40 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_40 * v_40 + dist_coeffs_3[int(5)] * (r2_40 + 2.0f * u_40 * u_40) + dist_coeffs_3[int(6)] * r2_40, _S471 * u_40 * v_40 + dist_coeffs_3[int(4)] * (r2_40 + 2.0f * v_40 * v_40) + dist_coeffs_3[int(7)] * r2_40);
    float2  _S473 = _S472 + make_float2 (dist_coeffs_3[int(8)] * _S472.x + dist_coeffs_3[int(9)] * _S472.y, 0.0f);
    uv_3[int(0)] = make_float2 (fx_3 * _S473.x + cx_3, fy_3 * _S473.y + cy_3);
    float2  _S474 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_11 = length_0(_S474);
    float _S475 = pos_c_3[int(1)].z;
    float theta_9 = (F32_atan2((r_11), (_S475)));
    if(theta_9 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_9 * theta_9 / 3.0f) / _S475;
    }
    else
    {
        k_1 = theta_9 / r_11;
    }
    float2  _S476 = _S474 * make_float2 (k_1);
    float u_41 = _S476.x;
    float v_41 = _S476.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float2  _S477 = _S476 * make_float2 (1.0f + r2_41 * (dist_coeffs_3[int(0)] + r2_41 * (dist_coeffs_3[int(1)] + r2_41 * (dist_coeffs_3[int(2)] + r2_41 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_41 * v_41 + dist_coeffs_3[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + dist_coeffs_3[int(6)] * r2_41, _S471 * u_41 * v_41 + dist_coeffs_3[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + dist_coeffs_3[int(7)] * r2_41);
    float2  _S478 = _S477 + make_float2 (dist_coeffs_3[int(8)] * _S477.x + dist_coeffs_3[int(9)] * _S477.y, 0.0f);
    uv_3[int(1)] = make_float2 (fx_3 * _S478.x + cx_3, fy_3 * _S478.y + cy_3);
    float2  _S479 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_12 = length_0(_S479);
    float _S480 = pos_c_3[int(2)].z;
    float theta_10 = (F32_atan2((r_12), (_S480)));
    if(theta_10 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_10 * theta_10 / 3.0f) / _S480;
    }
    else
    {
        k_1 = theta_10 / r_12;
    }
    float2  _S481 = _S479 * make_float2 (k_1);
    float u_42 = _S481.x;
    float v_42 = _S481.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float2  _S482 = _S481 * make_float2 (1.0f + r2_42 * (dist_coeffs_3[int(0)] + r2_42 * (dist_coeffs_3[int(1)] + r2_42 * (dist_coeffs_3[int(2)] + r2_42 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_42 * v_42 + dist_coeffs_3[int(5)] * (r2_42 + 2.0f * u_42 * u_42) + dist_coeffs_3[int(6)] * r2_42, _S471 * u_42 * v_42 + dist_coeffs_3[int(4)] * (r2_42 + 2.0f * v_42 * v_42) + dist_coeffs_3[int(7)] * r2_42);
    float2  _S483 = _S482 + make_float2 (dist_coeffs_3[int(8)] * _S482.x + dist_coeffs_3[int(9)] * _S482.y, 0.0f);
    uv_3[int(2)] = make_float2 (fx_3 * _S483.x + cx_3, fy_3 * _S483.y + cy_3);
    float2  _S484 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_13 = length_0(_S484);
    float _S485 = pos_c_3[int(3)].z;
    float theta_11 = (F32_atan2((r_13), (_S485)));
    if(theta_11 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_11 * theta_11 / 3.0f) / _S485;
    }
    else
    {
        k_1 = theta_11 / r_13;
    }
    float2  _S486 = _S484 * make_float2 (k_1);
    float u_43 = _S486.x;
    float v_43 = _S486.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float2  _S487 = _S486 * make_float2 (1.0f + r2_43 * (dist_coeffs_3[int(0)] + r2_43 * (dist_coeffs_3[int(1)] + r2_43 * (dist_coeffs_3[int(2)] + r2_43 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_43 * v_43 + dist_coeffs_3[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + dist_coeffs_3[int(6)] * r2_43, _S471 * u_43 * v_43 + dist_coeffs_3[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + dist_coeffs_3[int(7)] * r2_43);
    float2  _S488 = _S487 + make_float2 (dist_coeffs_3[int(8)] * _S487.x + dist_coeffs_3[int(9)] * _S487.y, 0.0f);
    uv_3[int(3)] = make_float2 (fx_3 * _S488.x + cx_3, fy_3 * _S488.y + cy_3);
    float2  _S489 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_14 = length_0(_S489);
    float _S490 = pos_c_3[int(4)].z;
    float theta_12 = (F32_atan2((r_14), (_S490)));
    if(theta_12 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_12 * theta_12 / 3.0f) / _S490;
    }
    else
    {
        k_1 = theta_12 / r_14;
    }
    float2  _S491 = _S489 * make_float2 (k_1);
    float u_44 = _S491.x;
    float v_44 = _S491.y;
    float r2_44 = u_44 * u_44 + v_44 * v_44;
    float2  _S492 = _S491 * make_float2 (1.0f + r2_44 * (dist_coeffs_3[int(0)] + r2_44 * (dist_coeffs_3[int(1)] + r2_44 * (dist_coeffs_3[int(2)] + r2_44 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_44 * v_44 + dist_coeffs_3[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + dist_coeffs_3[int(6)] * r2_44, _S471 * u_44 * v_44 + dist_coeffs_3[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + dist_coeffs_3[int(7)] * r2_44);
    float2  _S493 = _S492 + make_float2 (dist_coeffs_3[int(8)] * _S492.x + dist_coeffs_3[int(9)] * _S492.y, 0.0f);
    uv_3[int(4)] = make_float2 (fx_3 * _S493.x + cx_3, fy_3 * _S493.y + cy_3);
    float2  _S494 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_15 = length_0(_S494);
    float _S495 = pos_c_3[int(5)].z;
    float theta_13 = (F32_atan2((r_15), (_S495)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_13 * theta_13 / 3.0f) / _S495;
    }
    else
    {
        k_1 = theta_13 / r_15;
    }
    float2  _S496 = _S494 * make_float2 (k_1);
    float u_45 = _S496.x;
    float v_45 = _S496.y;
    float r2_45 = u_45 * u_45 + v_45 * v_45;
    float2  _S497 = _S496 * make_float2 (1.0f + r2_45 * (dist_coeffs_3[int(0)] + r2_45 * (dist_coeffs_3[int(1)] + r2_45 * (dist_coeffs_3[int(2)] + r2_45 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_45 * v_45 + dist_coeffs_3[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + dist_coeffs_3[int(6)] * r2_45, _S471 * u_45 * v_45 + dist_coeffs_3[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + dist_coeffs_3[int(7)] * r2_45);
    float2  _S498 = _S497 + make_float2 (dist_coeffs_3[int(8)] * _S497.x + dist_coeffs_3[int(9)] * _S497.y, 0.0f);
    uv_3[int(5)] = make_float2 (fx_3 * _S498.x + cx_3, fy_3 * _S498.y + cy_3);
    float2  _S499 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_16 = length_0(_S499);
    float _S500 = pos_c_3[int(6)].z;
    float theta_14 = (F32_atan2((r_16), (_S500)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_14 * theta_14 / 3.0f) / _S500;
    }
    else
    {
        k_1 = theta_14 / r_16;
    }
    float2  _S501 = _S499 * make_float2 (k_1);
    float u_46 = _S501.x;
    float v_46 = _S501.y;
    float r2_46 = u_46 * u_46 + v_46 * v_46;
    float2  _S502 = _S501 * make_float2 (1.0f + r2_46 * (dist_coeffs_3[int(0)] + r2_46 * (dist_coeffs_3[int(1)] + r2_46 * (dist_coeffs_3[int(2)] + r2_46 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_46 * v_46 + dist_coeffs_3[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + dist_coeffs_3[int(6)] * r2_46, _S471 * u_46 * v_46 + dist_coeffs_3[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + dist_coeffs_3[int(7)] * r2_46);
    float2  _S503 = _S502 + make_float2 (dist_coeffs_3[int(8)] * _S502.x + dist_coeffs_3[int(9)] * _S502.y, 0.0f);
    uv_3[int(6)] = make_float2 (fx_3 * _S503.x + cx_3, fy_3 * _S503.y + cy_3);
    float2  _S504 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_17 = length_0(_S504);
    float _S505 = pos_c_3[int(7)].z;
    float theta_15 = (F32_atan2((r_17), (_S505)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_15 * theta_15 / 3.0f) / _S505;
    }
    else
    {
        k_1 = theta_15 / r_17;
    }
    float2  _S506 = _S504 * make_float2 (k_1);
    float u_47 = _S506.x;
    float v_47 = _S506.y;
    float r2_47 = u_47 * u_47 + v_47 * v_47;
    float2  _S507 = _S506 * make_float2 (1.0f + r2_47 * (dist_coeffs_3[int(0)] + r2_47 * (dist_coeffs_3[int(1)] + r2_47 * (dist_coeffs_3[int(2)] + r2_47 * dist_coeffs_3[int(3)])))) + make_float2 (_S470 * u_47 * v_47 + dist_coeffs_3[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + dist_coeffs_3[int(6)] * r2_47, _S471 * u_47 * v_47 + dist_coeffs_3[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + dist_coeffs_3[int(7)] * r2_47);
    float2  _S508 = _S507 + make_float2 (dist_coeffs_3[int(8)] * _S507.x + dist_coeffs_3[int(9)] * _S507.y, 0.0f);
    float _S509 = fx_3 * _S508.x + cx_3;
    float _S510 = fy_3 * _S508.y + cy_3;
    uv_3[int(7)] = make_float2 (_S509, _S510);
    *aabb_xyxy_3 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S509))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S510))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S509))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S510))))))));
    float _S511 = mean_c_3.z;
    *depth_3 = (F32_max((_S511), (1.0f / (1.0f + (F32_sqrt((2.0f)))) * (_S511 + length_1(mean_c_3)))));
    float3  _S512 = mean_c_3 - - mul_0(transpose_1(R_3), t_3);
    float _S513 = _S512.x;
    float _S514 = _S512.y;
    float _S515 = _S512.z;
    float norm_3 = (F32_sqrt((_S513 * _S513 + _S514 * _S514 + _S515 * _S515)));
    float x_10 = _S513 / norm_3;
    float y_6 = _S514 / norm_3;
    float z_3 = _S515 / norm_3;
    float z2_3 = z_3 * z_3;
    float fTmp0B_3 = -1.09254848957061768f * z_3;
    float fC1_3 = x_10 * x_10 - y_6 * y_6;
    float fS1_3 = 2.0f * x_10 * y_6;
    float fTmp0C_3 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_3;
    *rgbs_3 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_3[int(1)] + make_float3 (z_3) * sh_coeffs_3[int(2)] - make_float3 (x_10) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_6) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_10) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_10 * fS1_3 + y_6 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_6) * sh_coeffs_3[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_10) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_10 * fC1_3 - y_6 * fS1_3)) * sh_coeffs_3[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S516, float3  _S517)
{
    return mul_0(_S516, _S517);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S518)
{
    return (F32_sqrt((_S518)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S519, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S520, float3  _S521)
{
    _d_max_vector_0(_S519, _S520, _S521);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S522, float _S523)
{
    _d_sqrt_0(_S522, _S523);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S524, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S525, float3  _S526)
{
    _d_mul_0(_S524, _S525, _S526);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float _s_dOut_0)
{
    float _S527 = (*dpx_6).primal_0.x;
    float _S528 = (*dpx_6).primal_0.y;
    float _S529 = (*dpx_6).primal_0.z;
    DiffPair_float_0 _S530;
    (&_S530)->primal_0 = _S527 * _S527 + _S528 * _S528 + _S529 * _S529;
    (&_S530)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S530, _s_dOut_0);
    float _S531 = (*dpx_6).primal_0.z * _S530.differential_0;
    float _S532 = _S531 + _S531;
    float _S533 = (*dpx_6).primal_0.y * _S530.differential_0;
    float _S534 = _S533 + _S533;
    float _S535 = (*dpx_6).primal_0.x * _S530.differential_0;
    float _S536 = _S535 + _S535;
    float3  _S537 = make_float3 (0.0f);
    *&((&_S537)->z) = _S532;
    *&((&_S537)->y) = _S534;
    *&((&_S537)->x) = _S536;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S537;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S538, float _S539)
{
    s_bwd_prop_length_impl_0(_S538, _S539);
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  densities_4, FixedArray<float3 , 16>  sh_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float3  v_rgb_0, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  _S540 = s_primal_ctx_mul_0(R_4, pos_4) + t_4;
    float _S541 = _S540.z;
    float2  _S542 = make_float2 (_S541);
    float _S543 = (F32_min((1.00000001504746622e+30f), (_S541)));
    float _S544 = (F32_max((0.0f), (_S541)));
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S545 = s_primal_ctx_mul_0(R_4, pos_i_0) + t_4;
    float _S546 = _S545.z;
    float2  _S547 = make_float2 (_S546);
    float _S548 = (F32_min((_S543), (_S546)));
    float _S549 = (F32_max((_S544), (_S546)));
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S550 = s_primal_ctx_mul_0(R_4, pos_i_1) + t_4;
    float _S551 = _S550.z;
    float2  _S552 = make_float2 (_S551);
    float _S553 = (F32_min((_S548), (_S551)));
    float _S554 = (F32_max((_S549), (_S551)));
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S555 = s_primal_ctx_mul_0(R_4, pos_i_2) + t_4;
    float _S556 = _S555.z;
    float2  _S557 = make_float2 (_S556);
    float _S558 = (F32_min((_S553), (_S556)));
    float _S559 = (F32_max((_S554), (_S556)));
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S560 = s_primal_ctx_mul_0(R_4, pos_i_3) + t_4;
    float _S561 = _S560.z;
    float2  _S562 = make_float2 (_S561);
    float _S563 = (F32_min((_S558), (_S561)));
    float _S564 = (F32_max((_S559), (_S561)));
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S565 = s_primal_ctx_mul_0(R_4, pos_i_4) + t_4;
    float _S566 = _S565.z;
    float2  _S567 = make_float2 (_S566);
    float _S568 = (F32_min((_S563), (_S566)));
    float _S569 = (F32_max((_S564), (_S566)));
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S570 = s_primal_ctx_mul_0(R_4, pos_i_5) + t_4;
    float _S571 = _S570.z;
    float2  _S572 = make_float2 (_S571);
    float _S573 = (F32_min((_S568), (_S571)));
    float _S574 = (F32_max((_S569), (_S571)));
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S575 = s_primal_ctx_mul_0(R_4, pos_i_6) + t_4;
    float _S576 = _S575.z;
    float2  _S577 = make_float2 (_S576);
    float3  _S578 = pos_4 + make_float3 (0.5f * size_4);
    float3  mean_c_4 = s_primal_ctx_mul_0(R_4, _S578) + t_4;
    float2  _S579 = float2 {_S540.x, _S540.y};
    float2  _S580 = _S579 / make_float2 (_S541);
    float2  _S581 = make_float2 (_S541 * _S541);
    float u_48 = _S580.x;
    float v_48 = _S580.y;
    float r2_48 = u_48 * u_48 + v_48 * v_48;
    float _S582 = dist_coeffs_4[int(2)] + r2_48 * dist_coeffs_4[int(3)];
    float _S583 = dist_coeffs_4[int(1)] + r2_48 * _S582;
    float _S584 = dist_coeffs_4[int(0)] + r2_48 * _S583;
    float radial_16 = 1.0f + r2_48 * _S584;
    float _S585 = 2.0f * dist_coeffs_4[int(4)];
    float _S586 = _S585 * u_48;
    float _S587 = 2.0f * u_48;
    float _S588 = 2.0f * dist_coeffs_4[int(5)];
    float _S589 = _S588 * u_48;
    float _S590 = 2.0f * v_48;
    float2  _S591 = _S580 * make_float2 (radial_16) + make_float2 (_S586 * v_48 + dist_coeffs_4[int(5)] * (r2_48 + _S587 * u_48) + dist_coeffs_4[int(6)] * r2_48, _S589 * v_48 + dist_coeffs_4[int(4)] * (r2_48 + _S590 * v_48) + dist_coeffs_4[int(7)] * r2_48);
    float2  _S592 = _S591 + make_float2 (dist_coeffs_4[int(8)] * _S591.x + dist_coeffs_4[int(9)] * _S591.y, 0.0f);
    float _S593 = fx_4 * _S592.x + cx_4;
    float _S594 = fy_4 * _S592.y + cy_4;
    float2  _S595 = float2 {_S545.x, _S545.y};
    float2  _S596 = _S595 / make_float2 (_S546);
    float2  _S597 = make_float2 (_S546 * _S546);
    float u_49 = _S596.x;
    float v_49 = _S596.y;
    float r2_49 = u_49 * u_49 + v_49 * v_49;
    float _S598 = dist_coeffs_4[int(2)] + r2_49 * dist_coeffs_4[int(3)];
    float _S599 = dist_coeffs_4[int(1)] + r2_49 * _S598;
    float _S600 = dist_coeffs_4[int(0)] + r2_49 * _S599;
    float radial_17 = 1.0f + r2_49 * _S600;
    float _S601 = _S585 * u_49;
    float _S602 = 2.0f * u_49;
    float _S603 = _S588 * u_49;
    float _S604 = 2.0f * v_49;
    float2  _S605 = _S596 * make_float2 (radial_17) + make_float2 (_S601 * v_49 + dist_coeffs_4[int(5)] * (r2_49 + _S602 * u_49) + dist_coeffs_4[int(6)] * r2_49, _S603 * v_49 + dist_coeffs_4[int(4)] * (r2_49 + _S604 * v_49) + dist_coeffs_4[int(7)] * r2_49);
    float2  _S606 = _S605 + make_float2 (dist_coeffs_4[int(8)] * _S605.x + dist_coeffs_4[int(9)] * _S605.y, 0.0f);
    float _S607 = fx_4 * _S606.x + cx_4;
    float _S608 = fy_4 * _S606.y + cy_4;
    float2  _S609 = float2 {_S550.x, _S550.y};
    float2  _S610 = _S609 / make_float2 (_S551);
    float2  _S611 = make_float2 (_S551 * _S551);
    float u_50 = _S610.x;
    float v_50 = _S610.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S612 = dist_coeffs_4[int(2)] + r2_50 * dist_coeffs_4[int(3)];
    float _S613 = dist_coeffs_4[int(1)] + r2_50 * _S612;
    float _S614 = dist_coeffs_4[int(0)] + r2_50 * _S613;
    float radial_18 = 1.0f + r2_50 * _S614;
    float _S615 = _S585 * u_50;
    float _S616 = 2.0f * u_50;
    float _S617 = _S588 * u_50;
    float _S618 = 2.0f * v_50;
    float2  _S619 = _S610 * make_float2 (radial_18) + make_float2 (_S615 * v_50 + dist_coeffs_4[int(5)] * (r2_50 + _S616 * u_50) + dist_coeffs_4[int(6)] * r2_50, _S617 * v_50 + dist_coeffs_4[int(4)] * (r2_50 + _S618 * v_50) + dist_coeffs_4[int(7)] * r2_50);
    float2  _S620 = _S619 + make_float2 (dist_coeffs_4[int(8)] * _S619.x + dist_coeffs_4[int(9)] * _S619.y, 0.0f);
    float _S621 = fx_4 * _S620.x + cx_4;
    float _S622 = fy_4 * _S620.y + cy_4;
    float2  _S623 = float2 {_S555.x, _S555.y};
    float2  _S624 = _S623 / make_float2 (_S556);
    float2  _S625 = make_float2 (_S556 * _S556);
    float u_51 = _S624.x;
    float v_51 = _S624.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float _S626 = dist_coeffs_4[int(2)] + r2_51 * dist_coeffs_4[int(3)];
    float _S627 = dist_coeffs_4[int(1)] + r2_51 * _S626;
    float _S628 = dist_coeffs_4[int(0)] + r2_51 * _S627;
    float radial_19 = 1.0f + r2_51 * _S628;
    float _S629 = _S585 * u_51;
    float _S630 = 2.0f * u_51;
    float _S631 = _S588 * u_51;
    float _S632 = 2.0f * v_51;
    float2  _S633 = _S624 * make_float2 (radial_19) + make_float2 (_S629 * v_51 + dist_coeffs_4[int(5)] * (r2_51 + _S630 * u_51) + dist_coeffs_4[int(6)] * r2_51, _S631 * v_51 + dist_coeffs_4[int(4)] * (r2_51 + _S632 * v_51) + dist_coeffs_4[int(7)] * r2_51);
    float2  _S634 = _S633 + make_float2 (dist_coeffs_4[int(8)] * _S633.x + dist_coeffs_4[int(9)] * _S633.y, 0.0f);
    float _S635 = fx_4 * _S634.x + cx_4;
    float _S636 = fy_4 * _S634.y + cy_4;
    float2  _S637 = float2 {_S560.x, _S560.y};
    float2  _S638 = _S637 / make_float2 (_S561);
    float2  _S639 = make_float2 (_S561 * _S561);
    float u_52 = _S638.x;
    float v_52 = _S638.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float _S640 = dist_coeffs_4[int(2)] + r2_52 * dist_coeffs_4[int(3)];
    float _S641 = dist_coeffs_4[int(1)] + r2_52 * _S640;
    float _S642 = dist_coeffs_4[int(0)] + r2_52 * _S641;
    float radial_20 = 1.0f + r2_52 * _S642;
    float _S643 = _S585 * u_52;
    float _S644 = 2.0f * u_52;
    float _S645 = _S588 * u_52;
    float _S646 = 2.0f * v_52;
    float2  _S647 = _S638 * make_float2 (radial_20) + make_float2 (_S643 * v_52 + dist_coeffs_4[int(5)] * (r2_52 + _S644 * u_52) + dist_coeffs_4[int(6)] * r2_52, _S645 * v_52 + dist_coeffs_4[int(4)] * (r2_52 + _S646 * v_52) + dist_coeffs_4[int(7)] * r2_52);
    float2  _S648 = _S647 + make_float2 (dist_coeffs_4[int(8)] * _S647.x + dist_coeffs_4[int(9)] * _S647.y, 0.0f);
    float _S649 = fx_4 * _S648.x + cx_4;
    float _S650 = fy_4 * _S648.y + cy_4;
    float2  _S651 = float2 {_S565.x, _S565.y};
    float2  _S652 = _S651 / make_float2 (_S566);
    float2  _S653 = make_float2 (_S566 * _S566);
    float u_53 = _S652.x;
    float v_53 = _S652.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S654 = dist_coeffs_4[int(2)] + r2_53 * dist_coeffs_4[int(3)];
    float _S655 = dist_coeffs_4[int(1)] + r2_53 * _S654;
    float _S656 = dist_coeffs_4[int(0)] + r2_53 * _S655;
    float radial_21 = 1.0f + r2_53 * _S656;
    float _S657 = _S585 * u_53;
    float _S658 = 2.0f * u_53;
    float _S659 = _S588 * u_53;
    float _S660 = 2.0f * v_53;
    float2  _S661 = _S652 * make_float2 (radial_21) + make_float2 (_S657 * v_53 + dist_coeffs_4[int(5)] * (r2_53 + _S658 * u_53) + dist_coeffs_4[int(6)] * r2_53, _S659 * v_53 + dist_coeffs_4[int(4)] * (r2_53 + _S660 * v_53) + dist_coeffs_4[int(7)] * r2_53);
    float2  _S662 = _S661 + make_float2 (dist_coeffs_4[int(8)] * _S661.x + dist_coeffs_4[int(9)] * _S661.y, 0.0f);
    float _S663 = fx_4 * _S662.x + cx_4;
    float _S664 = fy_4 * _S662.y + cy_4;
    float2  _S665 = float2 {_S570.x, _S570.y};
    float2  _S666 = _S665 / make_float2 (_S571);
    float2  _S667 = make_float2 (_S571 * _S571);
    float u_54 = _S666.x;
    float v_54 = _S666.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float _S668 = dist_coeffs_4[int(2)] + r2_54 * dist_coeffs_4[int(3)];
    float _S669 = dist_coeffs_4[int(1)] + r2_54 * _S668;
    float _S670 = dist_coeffs_4[int(0)] + r2_54 * _S669;
    float radial_22 = 1.0f + r2_54 * _S670;
    float _S671 = _S585 * u_54;
    float _S672 = 2.0f * u_54;
    float _S673 = _S588 * u_54;
    float _S674 = 2.0f * v_54;
    float2  _S675 = _S666 * make_float2 (radial_22) + make_float2 (_S671 * v_54 + dist_coeffs_4[int(5)] * (r2_54 + _S672 * u_54) + dist_coeffs_4[int(6)] * r2_54, _S673 * v_54 + dist_coeffs_4[int(4)] * (r2_54 + _S674 * v_54) + dist_coeffs_4[int(7)] * r2_54);
    float2  _S676 = _S675 + make_float2 (dist_coeffs_4[int(8)] * _S675.x + dist_coeffs_4[int(9)] * _S675.y, 0.0f);
    float _S677 = fx_4 * _S676.x + cx_4;
    float _S678 = fy_4 * _S676.y + cy_4;
    float2  _S679 = float2 {_S575.x, _S575.y};
    float2  _S680 = _S679 / make_float2 (_S576);
    float2  _S681 = make_float2 (_S576 * _S576);
    float u_55 = _S680.x;
    float v_55 = _S680.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float _S682 = dist_coeffs_4[int(2)] + r2_55 * dist_coeffs_4[int(3)];
    float _S683 = dist_coeffs_4[int(1)] + r2_55 * _S682;
    float _S684 = dist_coeffs_4[int(0)] + r2_55 * _S683;
    float radial_23 = 1.0f + r2_55 * _S684;
    float _S685 = _S585 * u_55;
    float _S686 = 2.0f * u_55;
    float _S687 = _S588 * u_55;
    float _S688 = 2.0f * v_55;
    float2  _S689 = _S680 * make_float2 (radial_23) + make_float2 (_S685 * v_55 + dist_coeffs_4[int(5)] * (r2_55 + _S686 * u_55) + dist_coeffs_4[int(6)] * r2_55, _S687 * v_55 + dist_coeffs_4[int(4)] * (r2_55 + _S688 * v_55) + dist_coeffs_4[int(7)] * r2_55);
    float2  _S690 = _S689 + make_float2 (dist_coeffs_4[int(8)] * _S689.x + dist_coeffs_4[int(9)] * _S689.y, 0.0f);
    float _S691 = fx_4 * _S690.x + cx_4;
    float _S692 = fy_4 * _S690.y + cy_4;
    float _S693 = (F32_max((_S593), (_S607)));
    float _S694 = (F32_min((_S593), (_S607)));
    float _S695 = (F32_max((_S594), (_S608)));
    float _S696 = (F32_min((_S594), (_S608)));
    float _S697 = (F32_max((_S693), (_S621)));
    float _S698 = (F32_min((_S694), (_S621)));
    float _S699 = (F32_max((_S695), (_S622)));
    float _S700 = (F32_min((_S696), (_S622)));
    float _S701 = (F32_max((_S697), (_S635)));
    float _S702 = (F32_min((_S698), (_S635)));
    float _S703 = (F32_max((_S699), (_S636)));
    float _S704 = (F32_min((_S700), (_S636)));
    float _S705 = (F32_max((_S701), (_S649)));
    float _S706 = (F32_min((_S702), (_S649)));
    float _S707 = (F32_max((_S703), (_S650)));
    float _S708 = (F32_min((_S704), (_S650)));
    float _S709 = (F32_max((_S705), (_S663)));
    float _S710 = (F32_min((_S706), (_S663)));
    float _S711 = (F32_max((_S707), (_S664)));
    float _S712 = (F32_min((_S708), (_S664)));
    float _S713 = (F32_max((_S709), (_S677)));
    float _S714 = (F32_min((_S710), (_S677)));
    float _S715 = (F32_max((_S711), (_S678)));
    float _S716 = (F32_min((_S712), (_S678)));
    float _S717 = mean_c_4.z;
    float _S718 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S719 = _S718 * (_S717 + length_1(mean_c_4));
    Matrix<float, 3, 3>  _S720 = transpose_1(R_4);
    float3  _S721 = mean_c_4 - - s_primal_ctx_mul_0(_S720, t_4);
    float _S722 = _S721.x;
    float _S723 = _S721.y;
    float _S724 = _S721.z;
    float _S725 = _S722 * _S722 + _S723 * _S723 + _S724 * _S724;
    float _S726 = s_primal_ctx_sqrt_0(_S725);
    float x_11 = _S722 / _S726;
    float3  _S727 = make_float3 (x_11);
    float _S728 = _S726 * _S726;
    float y_7 = _S723 / _S726;
    float z_4 = _S724 / _S726;
    float3  _S729 = make_float3 (z_4);
    float _S730 = - y_7;
    float3  _S731 = make_float3 (_S730);
    float z2_4 = z_4 * z_4;
    float fTmp0B_4 = -1.09254848957061768f * z_4;
    float fC1_4 = x_11 * x_11 - y_7 * y_7;
    float _S732 = 2.0f * x_11;
    float fS1_4 = _S732 * y_7;
    float pSH6_0 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float3  _S733 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_11;
    float3  _S734 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_7;
    float3  _S735 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S736 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S737 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_4;
    float _S738 = 1.86588168144226074f * z2_4 - 1.11952900886535645f;
    float pSH12_0 = z_4 * _S738;
    float3  _S739 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_11;
    float3  _S740 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_7;
    float3  _S741 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S742 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S743 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_11 * fC1_4 - y_7 * fS1_4);
    float3  _S744 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_11 * fS1_4 + y_7 * fC1_4);
    float3  _S745 = make_float3 (pSH9_0);
    float3  _S746 = make_float3 (0.0f);
    float3  _S747 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S748;
    (&_S748)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S730) * sh_coeffs_4[int(1)] + make_float3 (z_4) * sh_coeffs_4[int(2)] - make_float3 (x_11) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]) + make_float3 (0.5f);
    (&_S748)->differential_0 = _S747;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S749;
    (&_S749)->primal_0 = _S746;
    (&_S749)->differential_0 = _S747;
    s_bwd_prop_max_0(&_S748, &_S749, v_rgb_0);
    float3  _S750 = _S744 * _S748.differential_0;
    float3  _S751 = sh_coeffs_4[int(15)] * _S748.differential_0;
    float3  _S752 = _S742 * _S748.differential_0;
    float3  _S753 = sh_coeffs_4[int(14)] * _S748.differential_0;
    float3  _S754 = _S740 * _S748.differential_0;
    float3  _S755 = sh_coeffs_4[int(13)] * _S748.differential_0;
    float3  _S756 = _S739 * _S748.differential_0;
    float3  _S757 = sh_coeffs_4[int(12)] * _S748.differential_0;
    float3  _S758 = _S741 * _S748.differential_0;
    float3  _S759 = sh_coeffs_4[int(11)] * _S748.differential_0;
    float3  _S760 = _S743 * _S748.differential_0;
    float3  _S761 = sh_coeffs_4[int(10)] * _S748.differential_0;
    float3  _S762 = _S745 * _S748.differential_0;
    float3  _S763 = sh_coeffs_4[int(9)] * _S748.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S763.x + _S763.y + _S763.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S751.x + _S751.y + _S751.z);
    float _S764 = _S761.x + _S761.y + _S761.z;
    float _S765 = _S753.x + _S753.y + _S753.z;
    float _S766 = _S759.x + _S759.y + _S759.z;
    float _S767 = _S755.x + _S755.y + _S755.z;
    float _S768 = _S757.x + _S757.y + _S757.z;
    float _S769 = - s_diff_fC2_T_0;
    float3  _S770 = _S736 * _S748.differential_0;
    float3  _S771 = sh_coeffs_4[int(8)] * _S748.differential_0;
    float3  _S772 = _S734 * _S748.differential_0;
    float3  _S773 = sh_coeffs_4[int(7)] * _S748.differential_0;
    float3  _S774 = _S733 * _S748.differential_0;
    float3  _S775 = sh_coeffs_4[int(6)] * _S748.differential_0;
    float3  _S776 = _S735 * _S748.differential_0;
    float3  _S777 = sh_coeffs_4[int(5)] * _S748.differential_0;
    float3  _S778 = _S737 * _S748.differential_0;
    float3  _S779 = sh_coeffs_4[int(4)] * _S748.differential_0;
    float _S780 = _S777.x + _S777.y + _S777.z;
    float _S781 = _S773.x + _S773.y + _S773.z;
    float _S782 = fTmp1B_4 * _S764 + x_11 * s_diff_fS2_T_0 + y_7 * _S769 + 0.54627424478530884f * (_S779.x + _S779.y + _S779.z);
    float _S783 = fTmp1B_4 * _S765 + y_7 * s_diff_fS2_T_0 + x_11 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S771.x + _S771.y + _S771.z);
    float _S784 = y_7 * - _S783;
    float _S785 = x_11 * _S783;
    float _S786 = z_4 * (1.86588168144226074f * (z_4 * _S768) + -2.28522896766662598f * (y_7 * _S766 + x_11 * _S767) + 0.94617468118667603f * (_S775.x + _S775.y + _S775.z));
    float3  _S787 = make_float3 (0.48860251903533936f) * _S748.differential_0;
    float3  _S788 = - _S787;
    float3  _S789 = _S727 * _S788;
    float3  _S790 = sh_coeffs_4[int(3)] * _S788;
    float3  _S791 = _S729 * _S787;
    float3  _S792 = sh_coeffs_4[int(2)] * _S787;
    float3  _S793 = _S731 * _S787;
    float3  _S794 = sh_coeffs_4[int(1)] * _S787;
    float _S795 = (_S738 * _S768 + 1.44530570507049561f * (fS1_4 * _S764 + fC1_4 * _S765) + -1.09254848957061768f * (y_7 * _S780 + x_11 * _S781) + _S786 + _S786 + _S792.x + _S792.y + _S792.z) / _S728;
    float _S796 = _S726 * _S795;
    float _S797 = (fTmp0C_4 * _S766 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S769 + fTmp0B_4 * _S780 + _S732 * _S782 + _S784 + _S784 + - (_S794.x + _S794.y + _S794.z)) / _S728;
    float _S798 = _S726 * _S797;
    float _S799 = (fTmp0C_4 * _S767 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S781 + 2.0f * (y_7 * _S782) + _S785 + _S785 + _S790.x + _S790.y + _S790.z) / _S728;
    float _S800 = _S726 * _S799;
    float _S801 = _S724 * - _S795 + _S723 * - _S797 + _S722 * - _S799;
    DiffPair_float_0 _S802;
    (&_S802)->primal_0 = _S725;
    (&_S802)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S802, _S801);
    float _S803 = _S724 * _S802.differential_0;
    float _S804 = _S723 * _S802.differential_0;
    float _S805 = _S722 * _S802.differential_0;
    float3  _S806 = make_float3 (0.282094806432724f) * _S748.differential_0;
    float3  _S807 = make_float3 (_S800 + _S805 + _S805, _S798 + _S804 + _S804, _S796 + _S803 + _S803);
    float3  _S808 = - - _S807;
    Matrix<float, 3, 3>  _S809 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S810;
    (&_S810)->primal_0 = _S720;
    (&_S810)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S811;
    (&_S811)->primal_0 = t_4;
    (&_S811)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S810, &_S811, _S808);
    Matrix<float, 3, 3>  _S812 = transpose_1(_S810.differential_0);
    DiffPair_float_0 _S813;
    (&_S813)->primal_0 = _S717;
    (&_S813)->differential_0 = 0.0f;
    DiffPair_float_0 _S814;
    (&_S814)->primal_0 = _S719;
    (&_S814)->differential_0 = 0.0f;
    _d_max_0(&_S813, &_S814, 0.0f);
    float _S815 = _S718 * _S814.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S816;
    (&_S816)->primal_0 = mean_c_4;
    (&_S816)->differential_0 = _S747;
    s_bwd_length_impl_0(&_S816, _S815);
    float _S817 = _S813.differential_0 + _S815;
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = _S716;
    (&_S818)->differential_0 = 0.0f;
    DiffPair_float_0 _S819;
    (&_S819)->primal_0 = _S692;
    (&_S819)->differential_0 = 0.0f;
    _d_min_0(&_S818, &_S819, 0.0f);
    DiffPair_float_0 _S820;
    (&_S820)->primal_0 = _S715;
    (&_S820)->differential_0 = 0.0f;
    DiffPair_float_0 _S821;
    (&_S821)->primal_0 = _S692;
    (&_S821)->differential_0 = 0.0f;
    _d_max_0(&_S820, &_S821, 0.0f);
    DiffPair_float_0 _S822;
    (&_S822)->primal_0 = _S714;
    (&_S822)->differential_0 = 0.0f;
    DiffPair_float_0 _S823;
    (&_S823)->primal_0 = _S691;
    (&_S823)->differential_0 = 0.0f;
    _d_min_0(&_S822, &_S823, 0.0f);
    DiffPair_float_0 _S824;
    (&_S824)->primal_0 = _S713;
    (&_S824)->differential_0 = 0.0f;
    DiffPair_float_0 _S825;
    (&_S825)->primal_0 = _S691;
    (&_S825)->differential_0 = 0.0f;
    _d_max_0(&_S824, &_S825, 0.0f);
    DiffPair_float_0 _S826;
    (&_S826)->primal_0 = _S712;
    (&_S826)->differential_0 = 0.0f;
    DiffPair_float_0 _S827;
    (&_S827)->primal_0 = _S678;
    (&_S827)->differential_0 = 0.0f;
    _d_min_0(&_S826, &_S827, _S818.differential_0);
    DiffPair_float_0 _S828;
    (&_S828)->primal_0 = _S711;
    (&_S828)->differential_0 = 0.0f;
    DiffPair_float_0 _S829;
    (&_S829)->primal_0 = _S678;
    (&_S829)->differential_0 = 0.0f;
    _d_max_0(&_S828, &_S829, _S820.differential_0);
    DiffPair_float_0 _S830;
    (&_S830)->primal_0 = _S710;
    (&_S830)->differential_0 = 0.0f;
    DiffPair_float_0 _S831;
    (&_S831)->primal_0 = _S677;
    (&_S831)->differential_0 = 0.0f;
    _d_min_0(&_S830, &_S831, _S822.differential_0);
    DiffPair_float_0 _S832;
    (&_S832)->primal_0 = _S709;
    (&_S832)->differential_0 = 0.0f;
    DiffPair_float_0 _S833;
    (&_S833)->primal_0 = _S677;
    (&_S833)->differential_0 = 0.0f;
    _d_max_0(&_S832, &_S833, _S824.differential_0);
    DiffPair_float_0 _S834;
    (&_S834)->primal_0 = _S708;
    (&_S834)->differential_0 = 0.0f;
    DiffPair_float_0 _S835;
    (&_S835)->primal_0 = _S664;
    (&_S835)->differential_0 = 0.0f;
    _d_min_0(&_S834, &_S835, _S826.differential_0);
    DiffPair_float_0 _S836;
    (&_S836)->primal_0 = _S707;
    (&_S836)->differential_0 = 0.0f;
    DiffPair_float_0 _S837;
    (&_S837)->primal_0 = _S664;
    (&_S837)->differential_0 = 0.0f;
    _d_max_0(&_S836, &_S837, _S828.differential_0);
    DiffPair_float_0 _S838;
    (&_S838)->primal_0 = _S706;
    (&_S838)->differential_0 = 0.0f;
    DiffPair_float_0 _S839;
    (&_S839)->primal_0 = _S663;
    (&_S839)->differential_0 = 0.0f;
    _d_min_0(&_S838, &_S839, _S830.differential_0);
    DiffPair_float_0 _S840;
    (&_S840)->primal_0 = _S705;
    (&_S840)->differential_0 = 0.0f;
    DiffPair_float_0 _S841;
    (&_S841)->primal_0 = _S663;
    (&_S841)->differential_0 = 0.0f;
    _d_max_0(&_S840, &_S841, _S832.differential_0);
    DiffPair_float_0 _S842;
    (&_S842)->primal_0 = _S704;
    (&_S842)->differential_0 = 0.0f;
    DiffPair_float_0 _S843;
    (&_S843)->primal_0 = _S650;
    (&_S843)->differential_0 = 0.0f;
    _d_min_0(&_S842, &_S843, _S834.differential_0);
    DiffPair_float_0 _S844;
    (&_S844)->primal_0 = _S703;
    (&_S844)->differential_0 = 0.0f;
    DiffPair_float_0 _S845;
    (&_S845)->primal_0 = _S650;
    (&_S845)->differential_0 = 0.0f;
    _d_max_0(&_S844, &_S845, _S836.differential_0);
    DiffPair_float_0 _S846;
    (&_S846)->primal_0 = _S702;
    (&_S846)->differential_0 = 0.0f;
    DiffPair_float_0 _S847;
    (&_S847)->primal_0 = _S649;
    (&_S847)->differential_0 = 0.0f;
    _d_min_0(&_S846, &_S847, _S838.differential_0);
    DiffPair_float_0 _S848;
    (&_S848)->primal_0 = _S701;
    (&_S848)->differential_0 = 0.0f;
    DiffPair_float_0 _S849;
    (&_S849)->primal_0 = _S649;
    (&_S849)->differential_0 = 0.0f;
    _d_max_0(&_S848, &_S849, _S840.differential_0);
    DiffPair_float_0 _S850;
    (&_S850)->primal_0 = _S700;
    (&_S850)->differential_0 = 0.0f;
    DiffPair_float_0 _S851;
    (&_S851)->primal_0 = _S636;
    (&_S851)->differential_0 = 0.0f;
    _d_min_0(&_S850, &_S851, _S842.differential_0);
    DiffPair_float_0 _S852;
    (&_S852)->primal_0 = _S699;
    (&_S852)->differential_0 = 0.0f;
    DiffPair_float_0 _S853;
    (&_S853)->primal_0 = _S636;
    (&_S853)->differential_0 = 0.0f;
    _d_max_0(&_S852, &_S853, _S844.differential_0);
    DiffPair_float_0 _S854;
    (&_S854)->primal_0 = _S698;
    (&_S854)->differential_0 = 0.0f;
    DiffPair_float_0 _S855;
    (&_S855)->primal_0 = _S635;
    (&_S855)->differential_0 = 0.0f;
    _d_min_0(&_S854, &_S855, _S846.differential_0);
    DiffPair_float_0 _S856;
    (&_S856)->primal_0 = _S697;
    (&_S856)->differential_0 = 0.0f;
    DiffPair_float_0 _S857;
    (&_S857)->primal_0 = _S635;
    (&_S857)->differential_0 = 0.0f;
    _d_max_0(&_S856, &_S857, _S848.differential_0);
    DiffPair_float_0 _S858;
    (&_S858)->primal_0 = _S696;
    (&_S858)->differential_0 = 0.0f;
    DiffPair_float_0 _S859;
    (&_S859)->primal_0 = _S622;
    (&_S859)->differential_0 = 0.0f;
    _d_min_0(&_S858, &_S859, _S850.differential_0);
    DiffPair_float_0 _S860;
    (&_S860)->primal_0 = _S695;
    (&_S860)->differential_0 = 0.0f;
    DiffPair_float_0 _S861;
    (&_S861)->primal_0 = _S622;
    (&_S861)->differential_0 = 0.0f;
    _d_max_0(&_S860, &_S861, _S852.differential_0);
    DiffPair_float_0 _S862;
    (&_S862)->primal_0 = _S694;
    (&_S862)->differential_0 = 0.0f;
    DiffPair_float_0 _S863;
    (&_S863)->primal_0 = _S621;
    (&_S863)->differential_0 = 0.0f;
    _d_min_0(&_S862, &_S863, _S854.differential_0);
    DiffPair_float_0 _S864;
    (&_S864)->primal_0 = _S693;
    (&_S864)->differential_0 = 0.0f;
    DiffPair_float_0 _S865;
    (&_S865)->primal_0 = _S621;
    (&_S865)->differential_0 = 0.0f;
    _d_max_0(&_S864, &_S865, _S856.differential_0);
    DiffPair_float_0 _S866;
    (&_S866)->primal_0 = _S594;
    (&_S866)->differential_0 = 0.0f;
    DiffPair_float_0 _S867;
    (&_S867)->primal_0 = _S608;
    (&_S867)->differential_0 = 0.0f;
    _d_min_0(&_S866, &_S867, _S858.differential_0);
    DiffPair_float_0 _S868;
    (&_S868)->primal_0 = _S594;
    (&_S868)->differential_0 = 0.0f;
    DiffPair_float_0 _S869;
    (&_S869)->primal_0 = _S608;
    (&_S869)->differential_0 = 0.0f;
    _d_max_0(&_S868, &_S869, _S860.differential_0);
    DiffPair_float_0 _S870;
    (&_S870)->primal_0 = _S593;
    (&_S870)->differential_0 = 0.0f;
    DiffPair_float_0 _S871;
    (&_S871)->primal_0 = _S607;
    (&_S871)->differential_0 = 0.0f;
    _d_min_0(&_S870, &_S871, _S862.differential_0);
    DiffPair_float_0 _S872;
    (&_S872)->primal_0 = _S593;
    (&_S872)->differential_0 = 0.0f;
    DiffPair_float_0 _S873;
    (&_S873)->primal_0 = _S607;
    (&_S873)->differential_0 = 0.0f;
    _d_max_0(&_S872, &_S873, _S864.differential_0);
    float _S874 = fx_4 * (_S823.differential_0 + _S825.differential_0);
    float2  _S875 = make_float2 (_S874, fy_4 * (_S819.differential_0 + _S821.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S874, dist_coeffs_4[int(9)] * _S874);
    float2  _S876 = _S680 * _S875;
    float _S877 = dist_coeffs_4[int(4)] * _S875.y;
    float _S878 = dist_coeffs_4[int(5)] * _S875.x;
    float _S879 = _S876.x + _S876.y;
    float _S880 = r2_55 * _S879;
    float _S881 = r2_55 * _S880;
    float _S882 = dist_coeffs_4[int(7)] * _S875.y + _S877 + dist_coeffs_4[int(6)] * _S875.x + _S878 + _S684 * _S879 + _S683 * _S880 + _S682 * _S881 + dist_coeffs_4[int(3)] * (r2_55 * _S881);
    float _S883 = v_55 * _S882;
    float _S884 = u_55 * _S882;
    float2  _S885 = (make_float2 (radial_23) * _S875 + make_float2 (_S588 * (v_55 * _S875.y) + _S686 * _S878 + 2.0f * (u_55 * _S878) + _S585 * (v_55 * _S875.x) + _S884 + _S884, _S688 * _S877 + 2.0f * (v_55 * _S877) + _S687 * _S875.y + _S685 * _S875.x + _S883 + _S883)) / _S681;
    float2  _S886 = _S679 * - _S885;
    float2  _S887 = _S577 * _S885;
    float _S888 = fx_4 * (_S831.differential_0 + _S833.differential_0);
    float2  _S889 = make_float2 (_S888, fy_4 * (_S827.differential_0 + _S829.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S888, dist_coeffs_4[int(9)] * _S888);
    float2  _S890 = _S666 * _S889;
    float _S891 = dist_coeffs_4[int(4)] * _S889.y;
    float _S892 = dist_coeffs_4[int(5)] * _S889.x;
    float _S893 = _S890.x + _S890.y;
    float _S894 = r2_54 * _S893;
    float _S895 = r2_54 * _S894;
    float _S896 = dist_coeffs_4[int(7)] * _S889.y + _S891 + dist_coeffs_4[int(6)] * _S889.x + _S892 + _S670 * _S893 + _S669 * _S894 + _S668 * _S895 + dist_coeffs_4[int(3)] * (r2_54 * _S895);
    float _S897 = v_54 * _S896;
    float _S898 = u_54 * _S896;
    float2  _S899 = (make_float2 (radial_22) * _S889 + make_float2 (_S588 * (v_54 * _S889.y) + _S672 * _S892 + 2.0f * (u_54 * _S892) + _S585 * (v_54 * _S889.x) + _S898 + _S898, _S674 * _S891 + 2.0f * (v_54 * _S891) + _S673 * _S889.y + _S671 * _S889.x + _S897 + _S897)) / _S667;
    float2  _S900 = _S665 * - _S899;
    float2  _S901 = _S572 * _S899;
    float _S902 = fx_4 * (_S839.differential_0 + _S841.differential_0);
    float2  _S903 = make_float2 (_S902, fy_4 * (_S835.differential_0 + _S837.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S902, dist_coeffs_4[int(9)] * _S902);
    float2  _S904 = _S652 * _S903;
    float _S905 = dist_coeffs_4[int(4)] * _S903.y;
    float _S906 = dist_coeffs_4[int(5)] * _S903.x;
    float _S907 = _S904.x + _S904.y;
    float _S908 = r2_53 * _S907;
    float _S909 = r2_53 * _S908;
    float _S910 = dist_coeffs_4[int(7)] * _S903.y + _S905 + dist_coeffs_4[int(6)] * _S903.x + _S906 + _S656 * _S907 + _S655 * _S908 + _S654 * _S909 + dist_coeffs_4[int(3)] * (r2_53 * _S909);
    float _S911 = v_53 * _S910;
    float _S912 = u_53 * _S910;
    float2  _S913 = (make_float2 (radial_21) * _S903 + make_float2 (_S588 * (v_53 * _S903.y) + _S658 * _S906 + 2.0f * (u_53 * _S906) + _S585 * (v_53 * _S903.x) + _S912 + _S912, _S660 * _S905 + 2.0f * (v_53 * _S905) + _S659 * _S903.y + _S657 * _S903.x + _S911 + _S911)) / _S653;
    float2  _S914 = _S651 * - _S913;
    float2  _S915 = _S567 * _S913;
    float _S916 = fx_4 * (_S847.differential_0 + _S849.differential_0);
    float2  _S917 = make_float2 (_S916, fy_4 * (_S843.differential_0 + _S845.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S916, dist_coeffs_4[int(9)] * _S916);
    float2  _S918 = _S638 * _S917;
    float _S919 = dist_coeffs_4[int(4)] * _S917.y;
    float _S920 = dist_coeffs_4[int(5)] * _S917.x;
    float _S921 = _S918.x + _S918.y;
    float _S922 = r2_52 * _S921;
    float _S923 = r2_52 * _S922;
    float _S924 = dist_coeffs_4[int(7)] * _S917.y + _S919 + dist_coeffs_4[int(6)] * _S917.x + _S920 + _S642 * _S921 + _S641 * _S922 + _S640 * _S923 + dist_coeffs_4[int(3)] * (r2_52 * _S923);
    float _S925 = v_52 * _S924;
    float _S926 = u_52 * _S924;
    float2  _S927 = (make_float2 (radial_20) * _S917 + make_float2 (_S588 * (v_52 * _S917.y) + _S644 * _S920 + 2.0f * (u_52 * _S920) + _S585 * (v_52 * _S917.x) + _S926 + _S926, _S646 * _S919 + 2.0f * (v_52 * _S919) + _S645 * _S917.y + _S643 * _S917.x + _S925 + _S925)) / _S639;
    float2  _S928 = _S637 * - _S927;
    float2  _S929 = _S562 * _S927;
    float _S930 = fx_4 * (_S855.differential_0 + _S857.differential_0);
    float2  _S931 = make_float2 (_S930, fy_4 * (_S851.differential_0 + _S853.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S930, dist_coeffs_4[int(9)] * _S930);
    float2  _S932 = _S624 * _S931;
    float _S933 = dist_coeffs_4[int(4)] * _S931.y;
    float _S934 = dist_coeffs_4[int(5)] * _S931.x;
    float _S935 = _S932.x + _S932.y;
    float _S936 = r2_51 * _S935;
    float _S937 = r2_51 * _S936;
    float _S938 = dist_coeffs_4[int(7)] * _S931.y + _S933 + dist_coeffs_4[int(6)] * _S931.x + _S934 + _S628 * _S935 + _S627 * _S936 + _S626 * _S937 + dist_coeffs_4[int(3)] * (r2_51 * _S937);
    float _S939 = v_51 * _S938;
    float _S940 = u_51 * _S938;
    float2  _S941 = (make_float2 (radial_19) * _S931 + make_float2 (_S588 * (v_51 * _S931.y) + _S630 * _S934 + 2.0f * (u_51 * _S934) + _S585 * (v_51 * _S931.x) + _S940 + _S940, _S632 * _S933 + 2.0f * (v_51 * _S933) + _S631 * _S931.y + _S629 * _S931.x + _S939 + _S939)) / _S625;
    float2  _S942 = _S623 * - _S941;
    float2  _S943 = _S557 * _S941;
    float _S944 = fx_4 * (_S863.differential_0 + _S865.differential_0);
    float2  _S945 = make_float2 (_S944, fy_4 * (_S859.differential_0 + _S861.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S944, dist_coeffs_4[int(9)] * _S944);
    float2  _S946 = _S610 * _S945;
    float _S947 = dist_coeffs_4[int(4)] * _S945.y;
    float _S948 = dist_coeffs_4[int(5)] * _S945.x;
    float _S949 = _S946.x + _S946.y;
    float _S950 = r2_50 * _S949;
    float _S951 = r2_50 * _S950;
    float _S952 = dist_coeffs_4[int(7)] * _S945.y + _S947 + dist_coeffs_4[int(6)] * _S945.x + _S948 + _S614 * _S949 + _S613 * _S950 + _S612 * _S951 + dist_coeffs_4[int(3)] * (r2_50 * _S951);
    float _S953 = v_50 * _S952;
    float _S954 = u_50 * _S952;
    float2  _S955 = (make_float2 (radial_18) * _S945 + make_float2 (_S588 * (v_50 * _S945.y) + _S616 * _S948 + 2.0f * (u_50 * _S948) + _S585 * (v_50 * _S945.x) + _S954 + _S954, _S618 * _S947 + 2.0f * (v_50 * _S947) + _S617 * _S945.y + _S615 * _S945.x + _S953 + _S953)) / _S611;
    float2  _S956 = _S609 * - _S955;
    float2  _S957 = _S552 * _S955;
    float _S958 = fx_4 * (_S871.differential_0 + _S873.differential_0);
    float2  _S959 = make_float2 (_S958, fy_4 * (_S867.differential_0 + _S869.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S958, dist_coeffs_4[int(9)] * _S958);
    float2  _S960 = _S596 * _S959;
    float _S961 = dist_coeffs_4[int(4)] * _S959.y;
    float _S962 = dist_coeffs_4[int(5)] * _S959.x;
    float _S963 = _S960.x + _S960.y;
    float _S964 = r2_49 * _S963;
    float _S965 = r2_49 * _S964;
    float _S966 = dist_coeffs_4[int(7)] * _S959.y + _S961 + dist_coeffs_4[int(6)] * _S959.x + _S962 + _S600 * _S963 + _S599 * _S964 + _S598 * _S965 + dist_coeffs_4[int(3)] * (r2_49 * _S965);
    float _S967 = v_49 * _S966;
    float _S968 = u_49 * _S966;
    float2  _S969 = (make_float2 (radial_17) * _S959 + make_float2 (_S588 * (v_49 * _S959.y) + _S602 * _S962 + 2.0f * (u_49 * _S962) + _S585 * (v_49 * _S959.x) + _S968 + _S968, _S604 * _S961 + 2.0f * (v_49 * _S961) + _S603 * _S959.y + _S601 * _S959.x + _S967 + _S967)) / _S597;
    float2  _S970 = _S595 * - _S969;
    float2  _S971 = _S547 * _S969;
    float _S972 = fx_4 * (_S870.differential_0 + _S872.differential_0);
    float2  _S973 = make_float2 (_S972, fy_4 * (_S866.differential_0 + _S868.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S972, dist_coeffs_4[int(9)] * _S972);
    float2  _S974 = _S580 * _S973;
    float _S975 = dist_coeffs_4[int(4)] * _S973.y;
    float _S976 = dist_coeffs_4[int(5)] * _S973.x;
    float _S977 = _S974.x + _S974.y;
    float _S978 = r2_48 * _S977;
    float _S979 = r2_48 * _S978;
    float _S980 = dist_coeffs_4[int(7)] * _S973.y + _S975 + dist_coeffs_4[int(6)] * _S973.x + _S976 + _S584 * _S977 + _S583 * _S978 + _S582 * _S979 + dist_coeffs_4[int(3)] * (r2_48 * _S979);
    float _S981 = v_48 * _S980;
    float _S982 = u_48 * _S980;
    float2  _S983 = (make_float2 (radial_16) * _S973 + make_float2 (_S588 * (v_48 * _S973.y) + _S587 * _S976 + 2.0f * (u_48 * _S976) + _S585 * (v_48 * _S973.x) + _S982 + _S982, _S590 * _S975 + 2.0f * (v_48 * _S975) + _S589 * _S973.y + _S586 * _S973.x + _S981 + _S981)) / _S581;
    float2  _S984 = _S579 * - _S983;
    float2  _S985 = _S542 * _S983;
    float3  _S986 = _S807 + _S816.differential_0 + make_float3 (0.0f, 0.0f, _S817);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S987;
    (&_S987)->primal_0 = R_4;
    (&_S987)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S988;
    (&_S988)->primal_0 = _S578;
    (&_S988)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S987, &_S988, _S986);
    DiffPair_float_0 _S989;
    (&_S989)->primal_0 = _S574;
    (&_S989)->differential_0 = 0.0f;
    DiffPair_float_0 _S990;
    (&_S990)->primal_0 = _S576;
    (&_S990)->differential_0 = 0.0f;
    _d_max_0(&_S989, &_S990, 0.0f);
    DiffPair_float_0 _S991;
    (&_S991)->primal_0 = _S573;
    (&_S991)->differential_0 = 0.0f;
    DiffPair_float_0 _S992;
    (&_S992)->primal_0 = _S576;
    (&_S992)->differential_0 = 0.0f;
    _d_min_0(&_S991, &_S992, 0.0f);
    float3  _S993 = make_float3 (_S887.x, _S887.y, _S990.differential_0 + _S992.differential_0 + _S886.x + _S886.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S994;
    (&_S994)->primal_0 = R_4;
    (&_S994)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S995;
    (&_S995)->primal_0 = pos_i_6;
    (&_S995)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S994, &_S995, _S993);
    DiffPair_float_0 _S996;
    (&_S996)->primal_0 = _S569;
    (&_S996)->differential_0 = 0.0f;
    DiffPair_float_0 _S997;
    (&_S997)->primal_0 = _S571;
    (&_S997)->differential_0 = 0.0f;
    _d_max_0(&_S996, &_S997, _S989.differential_0);
    DiffPair_float_0 _S998;
    (&_S998)->primal_0 = _S568;
    (&_S998)->differential_0 = 0.0f;
    DiffPair_float_0 _S999;
    (&_S999)->primal_0 = _S571;
    (&_S999)->differential_0 = 0.0f;
    _d_min_0(&_S998, &_S999, _S991.differential_0);
    float3  _S1000 = make_float3 (_S901.x, _S901.y, _S997.differential_0 + _S999.differential_0 + _S900.x + _S900.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1001;
    (&_S1001)->primal_0 = R_4;
    (&_S1001)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1002;
    (&_S1002)->primal_0 = pos_i_5;
    (&_S1002)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1001, &_S1002, _S1000);
    DiffPair_float_0 _S1003;
    (&_S1003)->primal_0 = _S564;
    (&_S1003)->differential_0 = 0.0f;
    DiffPair_float_0 _S1004;
    (&_S1004)->primal_0 = _S566;
    (&_S1004)->differential_0 = 0.0f;
    _d_max_0(&_S1003, &_S1004, _S996.differential_0);
    DiffPair_float_0 _S1005;
    (&_S1005)->primal_0 = _S563;
    (&_S1005)->differential_0 = 0.0f;
    DiffPair_float_0 _S1006;
    (&_S1006)->primal_0 = _S566;
    (&_S1006)->differential_0 = 0.0f;
    _d_min_0(&_S1005, &_S1006, _S998.differential_0);
    float3  _S1007 = make_float3 (_S915.x, _S915.y, _S1004.differential_0 + _S1006.differential_0 + _S914.x + _S914.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1008;
    (&_S1008)->primal_0 = R_4;
    (&_S1008)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1009;
    (&_S1009)->primal_0 = pos_i_4;
    (&_S1009)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1008, &_S1009, _S1007);
    DiffPair_float_0 _S1010;
    (&_S1010)->primal_0 = _S559;
    (&_S1010)->differential_0 = 0.0f;
    DiffPair_float_0 _S1011;
    (&_S1011)->primal_0 = _S561;
    (&_S1011)->differential_0 = 0.0f;
    _d_max_0(&_S1010, &_S1011, _S1003.differential_0);
    DiffPair_float_0 _S1012;
    (&_S1012)->primal_0 = _S558;
    (&_S1012)->differential_0 = 0.0f;
    DiffPair_float_0 _S1013;
    (&_S1013)->primal_0 = _S561;
    (&_S1013)->differential_0 = 0.0f;
    _d_min_0(&_S1012, &_S1013, _S1005.differential_0);
    float3  _S1014 = make_float3 (_S929.x, _S929.y, _S1011.differential_0 + _S1013.differential_0 + _S928.x + _S928.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1015;
    (&_S1015)->primal_0 = R_4;
    (&_S1015)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1016;
    (&_S1016)->primal_0 = pos_i_3;
    (&_S1016)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1015, &_S1016, _S1014);
    DiffPair_float_0 _S1017;
    (&_S1017)->primal_0 = _S554;
    (&_S1017)->differential_0 = 0.0f;
    DiffPair_float_0 _S1018;
    (&_S1018)->primal_0 = _S556;
    (&_S1018)->differential_0 = 0.0f;
    _d_max_0(&_S1017, &_S1018, _S1010.differential_0);
    DiffPair_float_0 _S1019;
    (&_S1019)->primal_0 = _S553;
    (&_S1019)->differential_0 = 0.0f;
    DiffPair_float_0 _S1020;
    (&_S1020)->primal_0 = _S556;
    (&_S1020)->differential_0 = 0.0f;
    _d_min_0(&_S1019, &_S1020, _S1012.differential_0);
    float3  _S1021 = make_float3 (_S943.x, _S943.y, _S1018.differential_0 + _S1020.differential_0 + _S942.x + _S942.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1022;
    (&_S1022)->primal_0 = R_4;
    (&_S1022)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1023;
    (&_S1023)->primal_0 = pos_i_2;
    (&_S1023)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1022, &_S1023, _S1021);
    DiffPair_float_0 _S1024;
    (&_S1024)->primal_0 = _S549;
    (&_S1024)->differential_0 = 0.0f;
    DiffPair_float_0 _S1025;
    (&_S1025)->primal_0 = _S551;
    (&_S1025)->differential_0 = 0.0f;
    _d_max_0(&_S1024, &_S1025, _S1017.differential_0);
    DiffPair_float_0 _S1026;
    (&_S1026)->primal_0 = _S548;
    (&_S1026)->differential_0 = 0.0f;
    DiffPair_float_0 _S1027;
    (&_S1027)->primal_0 = _S551;
    (&_S1027)->differential_0 = 0.0f;
    _d_min_0(&_S1026, &_S1027, _S1019.differential_0);
    float3  _S1028 = make_float3 (_S957.x, _S957.y, _S1025.differential_0 + _S1027.differential_0 + _S956.x + _S956.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1029;
    (&_S1029)->primal_0 = R_4;
    (&_S1029)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1030;
    (&_S1030)->primal_0 = pos_i_1;
    (&_S1030)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1029, &_S1030, _S1028);
    DiffPair_float_0 _S1031;
    (&_S1031)->primal_0 = _S544;
    (&_S1031)->differential_0 = 0.0f;
    DiffPair_float_0 _S1032;
    (&_S1032)->primal_0 = _S546;
    (&_S1032)->differential_0 = 0.0f;
    _d_max_0(&_S1031, &_S1032, _S1024.differential_0);
    DiffPair_float_0 _S1033;
    (&_S1033)->primal_0 = _S543;
    (&_S1033)->differential_0 = 0.0f;
    DiffPair_float_0 _S1034;
    (&_S1034)->primal_0 = _S546;
    (&_S1034)->differential_0 = 0.0f;
    _d_min_0(&_S1033, &_S1034, _S1026.differential_0);
    float3  _S1035 = make_float3 (_S971.x, _S971.y, _S1032.differential_0 + _S1034.differential_0 + _S970.x + _S970.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1036;
    (&_S1036)->primal_0 = R_4;
    (&_S1036)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1037;
    (&_S1037)->primal_0 = pos_i_0;
    (&_S1037)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1036, &_S1037, _S1035);
    DiffPair_float_0 _S1038;
    (&_S1038)->primal_0 = 0.0f;
    (&_S1038)->differential_0 = 0.0f;
    DiffPair_float_0 _S1039;
    (&_S1039)->primal_0 = _S541;
    (&_S1039)->differential_0 = 0.0f;
    _d_max_0(&_S1038, &_S1039, _S1031.differential_0);
    DiffPair_float_0 _S1040;
    (&_S1040)->primal_0 = 1.00000001504746622e+30f;
    (&_S1040)->differential_0 = 0.0f;
    DiffPair_float_0 _S1041;
    (&_S1041)->primal_0 = _S541;
    (&_S1041)->differential_0 = 0.0f;
    _d_min_0(&_S1040, &_S1041, _S1033.differential_0);
    float3  _S1042 = make_float3 (_S985.x, _S985.y, _S1039.differential_0 + _S1041.differential_0 + _S984.x + _S984.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1043;
    (&_S1043)->primal_0 = R_4;
    (&_S1043)->differential_0 = _S809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1044;
    (&_S1044)->primal_0 = pos_4;
    (&_S1044)->differential_0 = _S747;
    s_bwd_prop_mul_0(&_S1043, &_S1044, _S1042);
    float3  _S1045 = _S811.differential_0 + _S986 + _S993 + _S1000 + _S1007 + _S1014 + _S1021 + _S1028 + _S1035 + _S1042;
    Matrix<float, 3, 3>  _S1046 = _S812 + _S987.differential_0 + _S994.differential_0 + _S1001.differential_0 + _S1008.differential_0 + _S1015.differential_0 + _S1022.differential_0 + _S1029.differential_0 + _S1036.differential_0 + _S1043.differential_0;
    FixedArray<float3 , 16>  _S1047;
    _S1047[int(0)] = _S747;
    _S1047[int(1)] = _S747;
    _S1047[int(2)] = _S747;
    _S1047[int(3)] = _S747;
    _S1047[int(4)] = _S747;
    _S1047[int(5)] = _S747;
    _S1047[int(6)] = _S747;
    _S1047[int(7)] = _S747;
    _S1047[int(8)] = _S747;
    _S1047[int(9)] = _S747;
    _S1047[int(10)] = _S747;
    _S1047[int(11)] = _S747;
    _S1047[int(12)] = _S747;
    _S1047[int(13)] = _S747;
    _S1047[int(14)] = _S747;
    _S1047[int(15)] = _S747;
    _S1047[int(15)] = _S750;
    _S1047[int(14)] = _S752;
    _S1047[int(13)] = _S754;
    _S1047[int(12)] = _S756;
    _S1047[int(11)] = _S758;
    _S1047[int(10)] = _S760;
    _S1047[int(9)] = _S762;
    _S1047[int(8)] = _S770;
    _S1047[int(7)] = _S772;
    _S1047[int(6)] = _S774;
    _S1047[int(5)] = _S776;
    _S1047[int(4)] = _S778;
    _S1047[int(3)] = _S789;
    _S1047[int(2)] = _S791;
    _S1047[int(1)] = _S793;
    _S1047[int(0)] = _S806;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_0 = _S1047;
    *v_R_0 = _S1046;
    *v_t_0 = _S1045;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S1048, float _S1049)
{
    return (F32_atan2((_S1048), (_S1049)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S1050, DiffPair_float_0 * _S1051, float _S1052)
{
    _d_atan2_0(_S1050, _S1051, _S1052);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_7, float _s_dOut_1)
{
    float _S1053 = (*dpx_7).primal_0.x;
    float _S1054 = (*dpx_7).primal_0.y;
    DiffPair_float_0 _S1055;
    (&_S1055)->primal_0 = _S1053 * _S1053 + _S1054 * _S1054;
    (&_S1055)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1055, _s_dOut_1);
    float _S1056 = (*dpx_7).primal_0.y * _S1055.differential_0;
    float _S1057 = _S1056 + _S1056;
    float _S1058 = (*dpx_7).primal_0.x * _S1055.differential_0;
    float _S1059 = _S1058 + _S1058;
    float2  _S1060 = make_float2 (0.0f);
    *&((&_S1060)->y) = _S1057;
    *&((&_S1060)->x) = _S1059;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S1060;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1061, float _S1062)
{
    s_bwd_prop_length_impl_1(_S1061, _S1062);
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  densities_5, FixedArray<float3 , 16>  sh_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float3  v_rgb_1, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  _S1063 = s_primal_ctx_mul_0(R_5, pos_5) + t_5;
    float _S1064 = length_1(_S1063);
    float _S1065 = (F32_min((1.00000001504746622e+30f), (_S1064)));
    float _S1066 = (F32_max((0.0f), (_S1064)));
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S1067 = s_primal_ctx_mul_0(R_5, pos_i_7) + t_5;
    float _S1068 = length_1(_S1067);
    float _S1069 = (F32_min((_S1065), (_S1068)));
    float _S1070 = (F32_max((_S1066), (_S1068)));
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S1071 = s_primal_ctx_mul_0(R_5, pos_i_8) + t_5;
    float _S1072 = length_1(_S1071);
    float _S1073 = (F32_min((_S1069), (_S1072)));
    float _S1074 = (F32_max((_S1070), (_S1072)));
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S1075 = s_primal_ctx_mul_0(R_5, pos_i_9) + t_5;
    float _S1076 = length_1(_S1075);
    float _S1077 = (F32_min((_S1073), (_S1076)));
    float _S1078 = (F32_max((_S1074), (_S1076)));
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S1079 = s_primal_ctx_mul_0(R_5, pos_i_10) + t_5;
    float _S1080 = length_1(_S1079);
    float _S1081 = (F32_min((_S1077), (_S1080)));
    float _S1082 = (F32_max((_S1078), (_S1080)));
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S1083 = s_primal_ctx_mul_0(R_5, pos_i_11) + t_5;
    float _S1084 = length_1(_S1083);
    float _S1085 = (F32_min((_S1081), (_S1084)));
    float _S1086 = (F32_max((_S1082), (_S1084)));
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S1087 = s_primal_ctx_mul_0(R_5, pos_i_12) + t_5;
    float _S1088 = length_1(_S1087);
    float _S1089 = (F32_min((_S1085), (_S1088)));
    float _S1090 = (F32_max((_S1086), (_S1088)));
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S1091 = s_primal_ctx_mul_0(R_5, pos_i_13) + t_5;
    float _S1092 = length_1(_S1091);
    float3  _S1093 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_5 = s_primal_ctx_mul_0(R_5, _S1093) + t_5;
    float2  _S1094 = float2 {_S1063.x, _S1063.y};
    float _S1095 = length_0(_S1094);
    float _S1096 = _S1063.z;
    float _S1097 = s_primal_ctx_atan2_0(_S1095, _S1096);
    bool _S1098 = _S1097 < 0.00100000004749745f;
    float k_2;
    float _S1099;
    float _S1100;
    float _S1101;
    if(_S1098)
    {
        float _S1102 = 1.0f - _S1097 * _S1097 / 3.0f;
        float _S1103 = _S1096 * _S1096;
        k_2 = _S1102 / _S1096;
        _S1099 = _S1103;
        _S1100 = _S1102;
        _S1101 = 0.0f;
    }
    else
    {
        float _S1104 = _S1095 * _S1095;
        k_2 = _S1097 / _S1095;
        _S1099 = 0.0f;
        _S1100 = 0.0f;
        _S1101 = _S1104;
    }
    float2  _S1105 = make_float2 (k_2);
    float2  _S1106 = _S1094 * make_float2 (k_2);
    float u_56 = _S1106.x;
    float v_56 = _S1106.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S1107 = dist_coeffs_5[int(2)] + r2_56 * dist_coeffs_5[int(3)];
    float _S1108 = dist_coeffs_5[int(1)] + r2_56 * _S1107;
    float _S1109 = dist_coeffs_5[int(0)] + r2_56 * _S1108;
    float radial_24 = 1.0f + r2_56 * _S1109;
    float _S1110 = 2.0f * dist_coeffs_5[int(4)];
    float _S1111 = _S1110 * u_56;
    float _S1112 = 2.0f * u_56;
    float _S1113 = 2.0f * dist_coeffs_5[int(5)];
    float _S1114 = _S1113 * u_56;
    float _S1115 = 2.0f * v_56;
    float2  _S1116 = _S1106 * make_float2 (radial_24) + make_float2 (_S1111 * v_56 + dist_coeffs_5[int(5)] * (r2_56 + _S1112 * u_56) + dist_coeffs_5[int(6)] * r2_56, _S1114 * v_56 + dist_coeffs_5[int(4)] * (r2_56 + _S1115 * v_56) + dist_coeffs_5[int(7)] * r2_56);
    float2  _S1117 = _S1116 + make_float2 (dist_coeffs_5[int(8)] * _S1116.x + dist_coeffs_5[int(9)] * _S1116.y, 0.0f);
    float _S1118 = fx_5 * _S1117.x + cx_5;
    float _S1119 = fy_5 * _S1117.y + cy_5;
    float2  _S1120 = float2 {_S1067.x, _S1067.y};
    float _S1121 = length_0(_S1120);
    float _S1122 = _S1067.z;
    float _S1123 = s_primal_ctx_atan2_0(_S1121, _S1122);
    bool _S1124 = _S1123 < 0.00100000004749745f;
    float _S1125;
    float _S1126;
    float _S1127;
    if(_S1124)
    {
        float _S1128 = 1.0f - _S1123 * _S1123 / 3.0f;
        float _S1129 = _S1122 * _S1122;
        k_2 = _S1128 / _S1122;
        _S1125 = _S1129;
        _S1126 = _S1128;
        _S1127 = 0.0f;
    }
    else
    {
        float _S1130 = _S1121 * _S1121;
        k_2 = _S1123 / _S1121;
        _S1125 = 0.0f;
        _S1126 = 0.0f;
        _S1127 = _S1130;
    }
    float2  _S1131 = make_float2 (k_2);
    float2  _S1132 = _S1120 * make_float2 (k_2);
    float u_57 = _S1132.x;
    float v_57 = _S1132.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S1133 = dist_coeffs_5[int(2)] + r2_57 * dist_coeffs_5[int(3)];
    float _S1134 = dist_coeffs_5[int(1)] + r2_57 * _S1133;
    float _S1135 = dist_coeffs_5[int(0)] + r2_57 * _S1134;
    float radial_25 = 1.0f + r2_57 * _S1135;
    float _S1136 = _S1110 * u_57;
    float _S1137 = 2.0f * u_57;
    float _S1138 = _S1113 * u_57;
    float _S1139 = 2.0f * v_57;
    float2  _S1140 = _S1132 * make_float2 (radial_25) + make_float2 (_S1136 * v_57 + dist_coeffs_5[int(5)] * (r2_57 + _S1137 * u_57) + dist_coeffs_5[int(6)] * r2_57, _S1138 * v_57 + dist_coeffs_5[int(4)] * (r2_57 + _S1139 * v_57) + dist_coeffs_5[int(7)] * r2_57);
    float2  _S1141 = _S1140 + make_float2 (dist_coeffs_5[int(8)] * _S1140.x + dist_coeffs_5[int(9)] * _S1140.y, 0.0f);
    float _S1142 = fx_5 * _S1141.x + cx_5;
    float _S1143 = fy_5 * _S1141.y + cy_5;
    float2  _S1144 = float2 {_S1071.x, _S1071.y};
    float _S1145 = length_0(_S1144);
    float _S1146 = _S1071.z;
    float _S1147 = s_primal_ctx_atan2_0(_S1145, _S1146);
    bool _S1148 = _S1147 < 0.00100000004749745f;
    float _S1149;
    float _S1150;
    float _S1151;
    if(_S1148)
    {
        float _S1152 = 1.0f - _S1147 * _S1147 / 3.0f;
        float _S1153 = _S1146 * _S1146;
        k_2 = _S1152 / _S1146;
        _S1149 = _S1153;
        _S1150 = _S1152;
        _S1151 = 0.0f;
    }
    else
    {
        float _S1154 = _S1145 * _S1145;
        k_2 = _S1147 / _S1145;
        _S1149 = 0.0f;
        _S1150 = 0.0f;
        _S1151 = _S1154;
    }
    float2  _S1155 = make_float2 (k_2);
    float2  _S1156 = _S1144 * make_float2 (k_2);
    float u_58 = _S1156.x;
    float v_58 = _S1156.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S1157 = dist_coeffs_5[int(2)] + r2_58 * dist_coeffs_5[int(3)];
    float _S1158 = dist_coeffs_5[int(1)] + r2_58 * _S1157;
    float _S1159 = dist_coeffs_5[int(0)] + r2_58 * _S1158;
    float radial_26 = 1.0f + r2_58 * _S1159;
    float _S1160 = _S1110 * u_58;
    float _S1161 = 2.0f * u_58;
    float _S1162 = _S1113 * u_58;
    float _S1163 = 2.0f * v_58;
    float2  _S1164 = _S1156 * make_float2 (radial_26) + make_float2 (_S1160 * v_58 + dist_coeffs_5[int(5)] * (r2_58 + _S1161 * u_58) + dist_coeffs_5[int(6)] * r2_58, _S1162 * v_58 + dist_coeffs_5[int(4)] * (r2_58 + _S1163 * v_58) + dist_coeffs_5[int(7)] * r2_58);
    float2  _S1165 = _S1164 + make_float2 (dist_coeffs_5[int(8)] * _S1164.x + dist_coeffs_5[int(9)] * _S1164.y, 0.0f);
    float _S1166 = fx_5 * _S1165.x + cx_5;
    float _S1167 = fy_5 * _S1165.y + cy_5;
    float2  _S1168 = float2 {_S1075.x, _S1075.y};
    float _S1169 = length_0(_S1168);
    float _S1170 = _S1075.z;
    float _S1171 = s_primal_ctx_atan2_0(_S1169, _S1170);
    bool _S1172 = _S1171 < 0.00100000004749745f;
    float _S1173;
    float _S1174;
    float _S1175;
    if(_S1172)
    {
        float _S1176 = 1.0f - _S1171 * _S1171 / 3.0f;
        float _S1177 = _S1170 * _S1170;
        k_2 = _S1176 / _S1170;
        _S1173 = _S1177;
        _S1174 = _S1176;
        _S1175 = 0.0f;
    }
    else
    {
        float _S1178 = _S1169 * _S1169;
        k_2 = _S1171 / _S1169;
        _S1173 = 0.0f;
        _S1174 = 0.0f;
        _S1175 = _S1178;
    }
    float2  _S1179 = make_float2 (k_2);
    float2  _S1180 = _S1168 * make_float2 (k_2);
    float u_59 = _S1180.x;
    float v_59 = _S1180.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S1181 = dist_coeffs_5[int(2)] + r2_59 * dist_coeffs_5[int(3)];
    float _S1182 = dist_coeffs_5[int(1)] + r2_59 * _S1181;
    float _S1183 = dist_coeffs_5[int(0)] + r2_59 * _S1182;
    float radial_27 = 1.0f + r2_59 * _S1183;
    float _S1184 = _S1110 * u_59;
    float _S1185 = 2.0f * u_59;
    float _S1186 = _S1113 * u_59;
    float _S1187 = 2.0f * v_59;
    float2  _S1188 = _S1180 * make_float2 (radial_27) + make_float2 (_S1184 * v_59 + dist_coeffs_5[int(5)] * (r2_59 + _S1185 * u_59) + dist_coeffs_5[int(6)] * r2_59, _S1186 * v_59 + dist_coeffs_5[int(4)] * (r2_59 + _S1187 * v_59) + dist_coeffs_5[int(7)] * r2_59);
    float2  _S1189 = _S1188 + make_float2 (dist_coeffs_5[int(8)] * _S1188.x + dist_coeffs_5[int(9)] * _S1188.y, 0.0f);
    float _S1190 = fx_5 * _S1189.x + cx_5;
    float _S1191 = fy_5 * _S1189.y + cy_5;
    float2  _S1192 = float2 {_S1079.x, _S1079.y};
    float _S1193 = length_0(_S1192);
    float _S1194 = _S1079.z;
    float _S1195 = s_primal_ctx_atan2_0(_S1193, _S1194);
    bool _S1196 = _S1195 < 0.00100000004749745f;
    float _S1197;
    float _S1198;
    float _S1199;
    if(_S1196)
    {
        float _S1200 = 1.0f - _S1195 * _S1195 / 3.0f;
        float _S1201 = _S1194 * _S1194;
        k_2 = _S1200 / _S1194;
        _S1197 = _S1201;
        _S1198 = _S1200;
        _S1199 = 0.0f;
    }
    else
    {
        float _S1202 = _S1193 * _S1193;
        k_2 = _S1195 / _S1193;
        _S1197 = 0.0f;
        _S1198 = 0.0f;
        _S1199 = _S1202;
    }
    float2  _S1203 = make_float2 (k_2);
    float2  _S1204 = _S1192 * make_float2 (k_2);
    float u_60 = _S1204.x;
    float v_60 = _S1204.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S1205 = dist_coeffs_5[int(2)] + r2_60 * dist_coeffs_5[int(3)];
    float _S1206 = dist_coeffs_5[int(1)] + r2_60 * _S1205;
    float _S1207 = dist_coeffs_5[int(0)] + r2_60 * _S1206;
    float radial_28 = 1.0f + r2_60 * _S1207;
    float _S1208 = _S1110 * u_60;
    float _S1209 = 2.0f * u_60;
    float _S1210 = _S1113 * u_60;
    float _S1211 = 2.0f * v_60;
    float2  _S1212 = _S1204 * make_float2 (radial_28) + make_float2 (_S1208 * v_60 + dist_coeffs_5[int(5)] * (r2_60 + _S1209 * u_60) + dist_coeffs_5[int(6)] * r2_60, _S1210 * v_60 + dist_coeffs_5[int(4)] * (r2_60 + _S1211 * v_60) + dist_coeffs_5[int(7)] * r2_60);
    float2  _S1213 = _S1212 + make_float2 (dist_coeffs_5[int(8)] * _S1212.x + dist_coeffs_5[int(9)] * _S1212.y, 0.0f);
    float _S1214 = fx_5 * _S1213.x + cx_5;
    float _S1215 = fy_5 * _S1213.y + cy_5;
    float2  _S1216 = float2 {_S1083.x, _S1083.y};
    float _S1217 = length_0(_S1216);
    float _S1218 = _S1083.z;
    float _S1219 = s_primal_ctx_atan2_0(_S1217, _S1218);
    bool _S1220 = _S1219 < 0.00100000004749745f;
    float _S1221;
    float _S1222;
    float _S1223;
    if(_S1220)
    {
        float _S1224 = 1.0f - _S1219 * _S1219 / 3.0f;
        float _S1225 = _S1218 * _S1218;
        k_2 = _S1224 / _S1218;
        _S1221 = _S1225;
        _S1222 = _S1224;
        _S1223 = 0.0f;
    }
    else
    {
        float _S1226 = _S1217 * _S1217;
        k_2 = _S1219 / _S1217;
        _S1221 = 0.0f;
        _S1222 = 0.0f;
        _S1223 = _S1226;
    }
    float2  _S1227 = make_float2 (k_2);
    float2  _S1228 = _S1216 * make_float2 (k_2);
    float u_61 = _S1228.x;
    float v_61 = _S1228.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S1229 = dist_coeffs_5[int(2)] + r2_61 * dist_coeffs_5[int(3)];
    float _S1230 = dist_coeffs_5[int(1)] + r2_61 * _S1229;
    float _S1231 = dist_coeffs_5[int(0)] + r2_61 * _S1230;
    float radial_29 = 1.0f + r2_61 * _S1231;
    float _S1232 = _S1110 * u_61;
    float _S1233 = 2.0f * u_61;
    float _S1234 = _S1113 * u_61;
    float _S1235 = 2.0f * v_61;
    float2  _S1236 = _S1228 * make_float2 (radial_29) + make_float2 (_S1232 * v_61 + dist_coeffs_5[int(5)] * (r2_61 + _S1233 * u_61) + dist_coeffs_5[int(6)] * r2_61, _S1234 * v_61 + dist_coeffs_5[int(4)] * (r2_61 + _S1235 * v_61) + dist_coeffs_5[int(7)] * r2_61);
    float2  _S1237 = _S1236 + make_float2 (dist_coeffs_5[int(8)] * _S1236.x + dist_coeffs_5[int(9)] * _S1236.y, 0.0f);
    float _S1238 = fx_5 * _S1237.x + cx_5;
    float _S1239 = fy_5 * _S1237.y + cy_5;
    float2  _S1240 = float2 {_S1087.x, _S1087.y};
    float _S1241 = length_0(_S1240);
    float _S1242 = _S1087.z;
    float _S1243 = s_primal_ctx_atan2_0(_S1241, _S1242);
    bool _S1244 = _S1243 < 0.00100000004749745f;
    float _S1245;
    float _S1246;
    float _S1247;
    if(_S1244)
    {
        float _S1248 = 1.0f - _S1243 * _S1243 / 3.0f;
        float _S1249 = _S1242 * _S1242;
        k_2 = _S1248 / _S1242;
        _S1245 = _S1249;
        _S1246 = _S1248;
        _S1247 = 0.0f;
    }
    else
    {
        float _S1250 = _S1241 * _S1241;
        k_2 = _S1243 / _S1241;
        _S1245 = 0.0f;
        _S1246 = 0.0f;
        _S1247 = _S1250;
    }
    float2  _S1251 = make_float2 (k_2);
    float2  _S1252 = _S1240 * make_float2 (k_2);
    float u_62 = _S1252.x;
    float v_62 = _S1252.y;
    float r2_62 = u_62 * u_62 + v_62 * v_62;
    float _S1253 = dist_coeffs_5[int(2)] + r2_62 * dist_coeffs_5[int(3)];
    float _S1254 = dist_coeffs_5[int(1)] + r2_62 * _S1253;
    float _S1255 = dist_coeffs_5[int(0)] + r2_62 * _S1254;
    float radial_30 = 1.0f + r2_62 * _S1255;
    float _S1256 = _S1110 * u_62;
    float _S1257 = 2.0f * u_62;
    float _S1258 = _S1113 * u_62;
    float _S1259 = 2.0f * v_62;
    float2  _S1260 = _S1252 * make_float2 (radial_30) + make_float2 (_S1256 * v_62 + dist_coeffs_5[int(5)] * (r2_62 + _S1257 * u_62) + dist_coeffs_5[int(6)] * r2_62, _S1258 * v_62 + dist_coeffs_5[int(4)] * (r2_62 + _S1259 * v_62) + dist_coeffs_5[int(7)] * r2_62);
    float2  _S1261 = _S1260 + make_float2 (dist_coeffs_5[int(8)] * _S1260.x + dist_coeffs_5[int(9)] * _S1260.y, 0.0f);
    float _S1262 = fx_5 * _S1261.x + cx_5;
    float _S1263 = fy_5 * _S1261.y + cy_5;
    float2  _S1264 = float2 {_S1091.x, _S1091.y};
    float _S1265 = length_0(_S1264);
    float _S1266 = _S1091.z;
    float _S1267 = s_primal_ctx_atan2_0(_S1265, _S1266);
    bool _S1268 = _S1267 < 0.00100000004749745f;
    float _S1269;
    float _S1270;
    float _S1271;
    if(_S1268)
    {
        float _S1272 = 1.0f - _S1267 * _S1267 / 3.0f;
        float _S1273 = _S1266 * _S1266;
        k_2 = _S1272 / _S1266;
        _S1269 = _S1273;
        _S1270 = _S1272;
        _S1271 = 0.0f;
    }
    else
    {
        float _S1274 = _S1265 * _S1265;
        k_2 = _S1267 / _S1265;
        _S1269 = 0.0f;
        _S1270 = 0.0f;
        _S1271 = _S1274;
    }
    float2  _S1275 = make_float2 (k_2);
    float2  _S1276 = _S1264 * make_float2 (k_2);
    float u_63 = _S1276.x;
    float v_63 = _S1276.y;
    float r2_63 = u_63 * u_63 + v_63 * v_63;
    float _S1277 = dist_coeffs_5[int(2)] + r2_63 * dist_coeffs_5[int(3)];
    float _S1278 = dist_coeffs_5[int(1)] + r2_63 * _S1277;
    float _S1279 = dist_coeffs_5[int(0)] + r2_63 * _S1278;
    float radial_31 = 1.0f + r2_63 * _S1279;
    float _S1280 = _S1110 * u_63;
    float _S1281 = 2.0f * u_63;
    float _S1282 = _S1113 * u_63;
    float _S1283 = 2.0f * v_63;
    float2  _S1284 = _S1276 * make_float2 (radial_31) + make_float2 (_S1280 * v_63 + dist_coeffs_5[int(5)] * (r2_63 + _S1281 * u_63) + dist_coeffs_5[int(6)] * r2_63, _S1282 * v_63 + dist_coeffs_5[int(4)] * (r2_63 + _S1283 * v_63) + dist_coeffs_5[int(7)] * r2_63);
    float2  _S1285 = _S1284 + make_float2 (dist_coeffs_5[int(8)] * _S1284.x + dist_coeffs_5[int(9)] * _S1284.y, 0.0f);
    float _S1286 = fx_5 * _S1285.x + cx_5;
    float _S1287 = fy_5 * _S1285.y + cy_5;
    float _S1288 = (F32_max((_S1118), (_S1142)));
    float _S1289 = (F32_min((_S1118), (_S1142)));
    float _S1290 = (F32_max((_S1119), (_S1143)));
    float _S1291 = (F32_min((_S1119), (_S1143)));
    float _S1292 = (F32_max((_S1288), (_S1166)));
    float _S1293 = (F32_min((_S1289), (_S1166)));
    float _S1294 = (F32_max((_S1290), (_S1167)));
    float _S1295 = (F32_min((_S1291), (_S1167)));
    float _S1296 = (F32_max((_S1292), (_S1190)));
    float _S1297 = (F32_min((_S1293), (_S1190)));
    float _S1298 = (F32_max((_S1294), (_S1191)));
    float _S1299 = (F32_min((_S1295), (_S1191)));
    float _S1300 = (F32_max((_S1296), (_S1214)));
    float _S1301 = (F32_min((_S1297), (_S1214)));
    float _S1302 = (F32_max((_S1298), (_S1215)));
    float _S1303 = (F32_min((_S1299), (_S1215)));
    float _S1304 = (F32_max((_S1300), (_S1238)));
    float _S1305 = (F32_min((_S1301), (_S1238)));
    float _S1306 = (F32_max((_S1302), (_S1239)));
    float _S1307 = (F32_min((_S1303), (_S1239)));
    float _S1308 = (F32_max((_S1304), (_S1262)));
    float _S1309 = (F32_min((_S1305), (_S1262)));
    float _S1310 = (F32_max((_S1306), (_S1263)));
    float _S1311 = (F32_min((_S1307), (_S1263)));
    float _S1312 = mean_c_5.z;
    float _S1313 = 1.0f / (1.0f + s_primal_ctx_sqrt_0(2.0f));
    float _S1314 = _S1313 * (_S1312 + length_1(mean_c_5));
    Matrix<float, 3, 3>  _S1315 = transpose_1(R_5);
    float3  _S1316 = mean_c_5 - - s_primal_ctx_mul_0(_S1315, t_5);
    float _S1317 = _S1316.x;
    float _S1318 = _S1316.y;
    float _S1319 = _S1316.z;
    float _S1320 = _S1317 * _S1317 + _S1318 * _S1318 + _S1319 * _S1319;
    float _S1321 = s_primal_ctx_sqrt_0(_S1320);
    float x_12 = _S1317 / _S1321;
    float3  _S1322 = make_float3 (x_12);
    float _S1323 = _S1321 * _S1321;
    float y_8 = _S1318 / _S1321;
    float z_5 = _S1319 / _S1321;
    float3  _S1324 = make_float3 (z_5);
    float _S1325 = - y_8;
    float3  _S1326 = make_float3 (_S1325);
    float z2_5 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_12 * x_12 - y_8 * y_8;
    float _S1327 = 2.0f * x_12;
    float fS1_5 = _S1327 * y_8;
    float pSH6_1 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float3  _S1328 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_12;
    float3  _S1329 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_8;
    float3  _S1330 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S1331 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S1332 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    float _S1333 = 1.86588168144226074f * z2_5 - 1.11952900886535645f;
    float pSH12_1 = z_5 * _S1333;
    float3  _S1334 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_12;
    float3  _S1335 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_8;
    float3  _S1336 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S1337 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S1338 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_12 * fC1_5 - y_8 * fS1_5);
    float3  _S1339 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_12 * fS1_5 + y_8 * fC1_5);
    float3  _S1340 = make_float3 (pSH9_1);
    float3  _S1341 = make_float3 (0.0f);
    float3  _S1342 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1343;
    (&_S1343)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1325) * sh_coeffs_5[int(1)] + make_float3 (z_5) * sh_coeffs_5[int(2)] - make_float3 (x_12) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]) + make_float3 (0.5f);
    (&_S1343)->differential_0 = _S1342;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1344;
    (&_S1344)->primal_0 = _S1341;
    (&_S1344)->differential_0 = _S1342;
    s_bwd_prop_max_0(&_S1343, &_S1344, v_rgb_1);
    float3  _S1345 = _S1339 * _S1343.differential_0;
    float3  _S1346 = sh_coeffs_5[int(15)] * _S1343.differential_0;
    float3  _S1347 = _S1337 * _S1343.differential_0;
    float3  _S1348 = sh_coeffs_5[int(14)] * _S1343.differential_0;
    float3  _S1349 = _S1335 * _S1343.differential_0;
    float3  _S1350 = sh_coeffs_5[int(13)] * _S1343.differential_0;
    float3  _S1351 = _S1334 * _S1343.differential_0;
    float3  _S1352 = sh_coeffs_5[int(12)] * _S1343.differential_0;
    float3  _S1353 = _S1336 * _S1343.differential_0;
    float3  _S1354 = sh_coeffs_5[int(11)] * _S1343.differential_0;
    float3  _S1355 = _S1338 * _S1343.differential_0;
    float3  _S1356 = sh_coeffs_5[int(10)] * _S1343.differential_0;
    float3  _S1357 = _S1340 * _S1343.differential_0;
    float3  _S1358 = sh_coeffs_5[int(9)] * _S1343.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1358.x + _S1358.y + _S1358.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1346.x + _S1346.y + _S1346.z);
    float _S1359 = _S1356.x + _S1356.y + _S1356.z;
    float _S1360 = _S1348.x + _S1348.y + _S1348.z;
    float _S1361 = _S1354.x + _S1354.y + _S1354.z;
    float _S1362 = _S1350.x + _S1350.y + _S1350.z;
    float _S1363 = _S1352.x + _S1352.y + _S1352.z;
    float _S1364 = - s_diff_fC2_T_1;
    float3  _S1365 = _S1331 * _S1343.differential_0;
    float3  _S1366 = sh_coeffs_5[int(8)] * _S1343.differential_0;
    float3  _S1367 = _S1329 * _S1343.differential_0;
    float3  _S1368 = sh_coeffs_5[int(7)] * _S1343.differential_0;
    float3  _S1369 = _S1328 * _S1343.differential_0;
    float3  _S1370 = sh_coeffs_5[int(6)] * _S1343.differential_0;
    float3  _S1371 = _S1330 * _S1343.differential_0;
    float3  _S1372 = sh_coeffs_5[int(5)] * _S1343.differential_0;
    float3  _S1373 = _S1332 * _S1343.differential_0;
    float3  _S1374 = sh_coeffs_5[int(4)] * _S1343.differential_0;
    float _S1375 = _S1372.x + _S1372.y + _S1372.z;
    float _S1376 = _S1368.x + _S1368.y + _S1368.z;
    float _S1377 = fTmp1B_5 * _S1359 + x_12 * s_diff_fS2_T_1 + y_8 * _S1364 + 0.54627424478530884f * (_S1374.x + _S1374.y + _S1374.z);
    float _S1378 = fTmp1B_5 * _S1360 + y_8 * s_diff_fS2_T_1 + x_12 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1366.x + _S1366.y + _S1366.z);
    float _S1379 = y_8 * - _S1378;
    float _S1380 = x_12 * _S1378;
    float _S1381 = z_5 * (1.86588168144226074f * (z_5 * _S1363) + -2.28522896766662598f * (y_8 * _S1361 + x_12 * _S1362) + 0.94617468118667603f * (_S1370.x + _S1370.y + _S1370.z));
    float3  _S1382 = make_float3 (0.48860251903533936f) * _S1343.differential_0;
    float3  _S1383 = - _S1382;
    float3  _S1384 = _S1322 * _S1383;
    float3  _S1385 = sh_coeffs_5[int(3)] * _S1383;
    float3  _S1386 = _S1324 * _S1382;
    float3  _S1387 = sh_coeffs_5[int(2)] * _S1382;
    float3  _S1388 = _S1326 * _S1382;
    float3  _S1389 = sh_coeffs_5[int(1)] * _S1382;
    float _S1390 = (_S1333 * _S1363 + 1.44530570507049561f * (fS1_5 * _S1359 + fC1_5 * _S1360) + -1.09254848957061768f * (y_8 * _S1375 + x_12 * _S1376) + _S1381 + _S1381 + _S1387.x + _S1387.y + _S1387.z) / _S1323;
    float _S1391 = _S1321 * _S1390;
    float _S1392 = (fTmp0C_5 * _S1361 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S1364 + fTmp0B_5 * _S1375 + _S1327 * _S1377 + _S1379 + _S1379 + - (_S1389.x + _S1389.y + _S1389.z)) / _S1323;
    float _S1393 = _S1321 * _S1392;
    float _S1394 = (fTmp0C_5 * _S1362 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S1376 + 2.0f * (y_8 * _S1377) + _S1380 + _S1380 + _S1385.x + _S1385.y + _S1385.z) / _S1323;
    float _S1395 = _S1321 * _S1394;
    float _S1396 = _S1319 * - _S1390 + _S1318 * - _S1392 + _S1317 * - _S1394;
    DiffPair_float_0 _S1397;
    (&_S1397)->primal_0 = _S1320;
    (&_S1397)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1397, _S1396);
    float _S1398 = _S1319 * _S1397.differential_0;
    float _S1399 = _S1318 * _S1397.differential_0;
    float _S1400 = _S1317 * _S1397.differential_0;
    float3  _S1401 = make_float3 (0.282094806432724f) * _S1343.differential_0;
    float3  _S1402 = make_float3 (_S1395 + _S1400 + _S1400, _S1393 + _S1399 + _S1399, _S1391 + _S1398 + _S1398);
    float3  _S1403 = - - _S1402;
    Matrix<float, 3, 3>  _S1404 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1405;
    (&_S1405)->primal_0 = _S1315;
    (&_S1405)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1406;
    (&_S1406)->primal_0 = t_5;
    (&_S1406)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1405, &_S1406, _S1403);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1407 = _S1406;
    Matrix<float, 3, 3>  _S1408 = transpose_1(_S1405.differential_0);
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1312;
    (&_S1409)->differential_0 = 0.0f;
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = _S1314;
    (&_S1410)->differential_0 = 0.0f;
    _d_max_0(&_S1409, &_S1410, 0.0f);
    float _S1411 = _S1313 * _S1410.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1412;
    (&_S1412)->primal_0 = mean_c_5;
    (&_S1412)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1412, _S1411);
    float _S1413 = _S1409.differential_0 + _S1411;
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = _S1311;
    (&_S1414)->differential_0 = 0.0f;
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1287;
    (&_S1415)->differential_0 = 0.0f;
    _d_min_0(&_S1414, &_S1415, 0.0f);
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1310;
    (&_S1416)->differential_0 = 0.0f;
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1287;
    (&_S1417)->differential_0 = 0.0f;
    _d_max_0(&_S1416, &_S1417, 0.0f);
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = _S1309;
    (&_S1418)->differential_0 = 0.0f;
    DiffPair_float_0 _S1419;
    (&_S1419)->primal_0 = _S1286;
    (&_S1419)->differential_0 = 0.0f;
    _d_min_0(&_S1418, &_S1419, 0.0f);
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = _S1308;
    (&_S1420)->differential_0 = 0.0f;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1286;
    (&_S1421)->differential_0 = 0.0f;
    _d_max_0(&_S1420, &_S1421, 0.0f);
    DiffPair_float_0 _S1422;
    (&_S1422)->primal_0 = _S1307;
    (&_S1422)->differential_0 = 0.0f;
    DiffPair_float_0 _S1423;
    (&_S1423)->primal_0 = _S1263;
    (&_S1423)->differential_0 = 0.0f;
    _d_min_0(&_S1422, &_S1423, _S1414.differential_0);
    DiffPair_float_0 _S1424;
    (&_S1424)->primal_0 = _S1306;
    (&_S1424)->differential_0 = 0.0f;
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = _S1263;
    (&_S1425)->differential_0 = 0.0f;
    _d_max_0(&_S1424, &_S1425, _S1416.differential_0);
    DiffPair_float_0 _S1426;
    (&_S1426)->primal_0 = _S1305;
    (&_S1426)->differential_0 = 0.0f;
    DiffPair_float_0 _S1427;
    (&_S1427)->primal_0 = _S1262;
    (&_S1427)->differential_0 = 0.0f;
    _d_min_0(&_S1426, &_S1427, _S1418.differential_0);
    DiffPair_float_0 _S1428;
    (&_S1428)->primal_0 = _S1304;
    (&_S1428)->differential_0 = 0.0f;
    DiffPair_float_0 _S1429;
    (&_S1429)->primal_0 = _S1262;
    (&_S1429)->differential_0 = 0.0f;
    _d_max_0(&_S1428, &_S1429, _S1420.differential_0);
    DiffPair_float_0 _S1430;
    (&_S1430)->primal_0 = _S1303;
    (&_S1430)->differential_0 = 0.0f;
    DiffPair_float_0 _S1431;
    (&_S1431)->primal_0 = _S1239;
    (&_S1431)->differential_0 = 0.0f;
    _d_min_0(&_S1430, &_S1431, _S1422.differential_0);
    DiffPair_float_0 _S1432;
    (&_S1432)->primal_0 = _S1302;
    (&_S1432)->differential_0 = 0.0f;
    DiffPair_float_0 _S1433;
    (&_S1433)->primal_0 = _S1239;
    (&_S1433)->differential_0 = 0.0f;
    _d_max_0(&_S1432, &_S1433, _S1424.differential_0);
    DiffPair_float_0 _S1434;
    (&_S1434)->primal_0 = _S1301;
    (&_S1434)->differential_0 = 0.0f;
    DiffPair_float_0 _S1435;
    (&_S1435)->primal_0 = _S1238;
    (&_S1435)->differential_0 = 0.0f;
    _d_min_0(&_S1434, &_S1435, _S1426.differential_0);
    DiffPair_float_0 _S1436;
    (&_S1436)->primal_0 = _S1300;
    (&_S1436)->differential_0 = 0.0f;
    DiffPair_float_0 _S1437;
    (&_S1437)->primal_0 = _S1238;
    (&_S1437)->differential_0 = 0.0f;
    _d_max_0(&_S1436, &_S1437, _S1428.differential_0);
    DiffPair_float_0 _S1438;
    (&_S1438)->primal_0 = _S1299;
    (&_S1438)->differential_0 = 0.0f;
    DiffPair_float_0 _S1439;
    (&_S1439)->primal_0 = _S1215;
    (&_S1439)->differential_0 = 0.0f;
    _d_min_0(&_S1438, &_S1439, _S1430.differential_0);
    DiffPair_float_0 _S1440;
    (&_S1440)->primal_0 = _S1298;
    (&_S1440)->differential_0 = 0.0f;
    DiffPair_float_0 _S1441;
    (&_S1441)->primal_0 = _S1215;
    (&_S1441)->differential_0 = 0.0f;
    _d_max_0(&_S1440, &_S1441, _S1432.differential_0);
    DiffPair_float_0 _S1442;
    (&_S1442)->primal_0 = _S1297;
    (&_S1442)->differential_0 = 0.0f;
    DiffPair_float_0 _S1443;
    (&_S1443)->primal_0 = _S1214;
    (&_S1443)->differential_0 = 0.0f;
    _d_min_0(&_S1442, &_S1443, _S1434.differential_0);
    DiffPair_float_0 _S1444;
    (&_S1444)->primal_0 = _S1296;
    (&_S1444)->differential_0 = 0.0f;
    DiffPair_float_0 _S1445;
    (&_S1445)->primal_0 = _S1214;
    (&_S1445)->differential_0 = 0.0f;
    _d_max_0(&_S1444, &_S1445, _S1436.differential_0);
    DiffPair_float_0 _S1446;
    (&_S1446)->primal_0 = _S1295;
    (&_S1446)->differential_0 = 0.0f;
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1191;
    (&_S1447)->differential_0 = 0.0f;
    _d_min_0(&_S1446, &_S1447, _S1438.differential_0);
    DiffPair_float_0 _S1448;
    (&_S1448)->primal_0 = _S1294;
    (&_S1448)->differential_0 = 0.0f;
    DiffPair_float_0 _S1449;
    (&_S1449)->primal_0 = _S1191;
    (&_S1449)->differential_0 = 0.0f;
    _d_max_0(&_S1448, &_S1449, _S1440.differential_0);
    DiffPair_float_0 _S1450;
    (&_S1450)->primal_0 = _S1293;
    (&_S1450)->differential_0 = 0.0f;
    DiffPair_float_0 _S1451;
    (&_S1451)->primal_0 = _S1190;
    (&_S1451)->differential_0 = 0.0f;
    _d_min_0(&_S1450, &_S1451, _S1442.differential_0);
    DiffPair_float_0 _S1452;
    (&_S1452)->primal_0 = _S1292;
    (&_S1452)->differential_0 = 0.0f;
    DiffPair_float_0 _S1453;
    (&_S1453)->primal_0 = _S1190;
    (&_S1453)->differential_0 = 0.0f;
    _d_max_0(&_S1452, &_S1453, _S1444.differential_0);
    DiffPair_float_0 _S1454;
    (&_S1454)->primal_0 = _S1291;
    (&_S1454)->differential_0 = 0.0f;
    DiffPair_float_0 _S1455;
    (&_S1455)->primal_0 = _S1167;
    (&_S1455)->differential_0 = 0.0f;
    _d_min_0(&_S1454, &_S1455, _S1446.differential_0);
    DiffPair_float_0 _S1456;
    (&_S1456)->primal_0 = _S1290;
    (&_S1456)->differential_0 = 0.0f;
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1167;
    (&_S1457)->differential_0 = 0.0f;
    _d_max_0(&_S1456, &_S1457, _S1448.differential_0);
    DiffPair_float_0 _S1458;
    (&_S1458)->primal_0 = _S1289;
    (&_S1458)->differential_0 = 0.0f;
    DiffPair_float_0 _S1459;
    (&_S1459)->primal_0 = _S1166;
    (&_S1459)->differential_0 = 0.0f;
    _d_min_0(&_S1458, &_S1459, _S1450.differential_0);
    DiffPair_float_0 _S1460;
    (&_S1460)->primal_0 = _S1288;
    (&_S1460)->differential_0 = 0.0f;
    DiffPair_float_0 _S1461;
    (&_S1461)->primal_0 = _S1166;
    (&_S1461)->differential_0 = 0.0f;
    _d_max_0(&_S1460, &_S1461, _S1452.differential_0);
    DiffPair_float_0 _S1462;
    (&_S1462)->primal_0 = _S1119;
    (&_S1462)->differential_0 = 0.0f;
    DiffPair_float_0 _S1463;
    (&_S1463)->primal_0 = _S1143;
    (&_S1463)->differential_0 = 0.0f;
    _d_min_0(&_S1462, &_S1463, _S1454.differential_0);
    DiffPair_float_0 _S1464;
    (&_S1464)->primal_0 = _S1119;
    (&_S1464)->differential_0 = 0.0f;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = _S1143;
    (&_S1465)->differential_0 = 0.0f;
    _d_max_0(&_S1464, &_S1465, _S1456.differential_0);
    DiffPair_float_0 _S1466;
    (&_S1466)->primal_0 = _S1118;
    (&_S1466)->differential_0 = 0.0f;
    DiffPair_float_0 _S1467;
    (&_S1467)->primal_0 = _S1142;
    (&_S1467)->differential_0 = 0.0f;
    _d_min_0(&_S1466, &_S1467, _S1458.differential_0);
    DiffPair_float_0 _S1468;
    (&_S1468)->primal_0 = _S1118;
    (&_S1468)->differential_0 = 0.0f;
    DiffPair_float_0 _S1469;
    (&_S1469)->primal_0 = _S1142;
    (&_S1469)->differential_0 = 0.0f;
    _d_max_0(&_S1468, &_S1469, _S1460.differential_0);
    float _S1470 = fx_5 * (_S1419.differential_0 + _S1421.differential_0);
    float2  _S1471 = make_float2 (_S1470, fy_5 * (_S1415.differential_0 + _S1417.differential_0)) + make_float2 (dist_coeffs_5[int(8)] * _S1470, dist_coeffs_5[int(9)] * _S1470);
    float2  _S1472 = _S1276 * _S1471;
    float _S1473 = dist_coeffs_5[int(4)] * _S1471.y;
    float _S1474 = dist_coeffs_5[int(5)] * _S1471.x;
    float _S1475 = _S1472.x + _S1472.y;
    float _S1476 = r2_63 * _S1475;
    float _S1477 = r2_63 * _S1476;
    float _S1478 = dist_coeffs_5[int(7)] * _S1471.y + _S1473 + dist_coeffs_5[int(6)] * _S1471.x + _S1474 + _S1279 * _S1475 + _S1278 * _S1476 + _S1277 * _S1477 + dist_coeffs_5[int(3)] * (r2_63 * _S1477);
    float _S1479 = v_63 * _S1478;
    float _S1480 = u_63 * _S1478;
    float2  _S1481 = make_float2 (radial_31) * _S1471 + make_float2 (_S1113 * (v_63 * _S1471.y) + _S1281 * _S1474 + 2.0f * (u_63 * _S1474) + _S1110 * (v_63 * _S1471.x) + _S1480 + _S1480, _S1283 * _S1473 + 2.0f * (v_63 * _S1473) + _S1282 * _S1471.y + _S1280 * _S1471.x + _S1479 + _S1479);
    FixedArray<float3 , 16>  _S1482;
    _S1482[int(0)] = _S1342;
    _S1482[int(1)] = _S1342;
    _S1482[int(2)] = _S1342;
    _S1482[int(3)] = _S1342;
    _S1482[int(4)] = _S1342;
    _S1482[int(5)] = _S1342;
    _S1482[int(6)] = _S1342;
    _S1482[int(7)] = _S1342;
    _S1482[int(8)] = _S1342;
    _S1482[int(9)] = _S1342;
    _S1482[int(10)] = _S1342;
    _S1482[int(11)] = _S1342;
    _S1482[int(12)] = _S1342;
    _S1482[int(13)] = _S1342;
    _S1482[int(14)] = _S1342;
    _S1482[int(15)] = _S1342;
    _S1482[int(7)] = _S1367;
    _S1482[int(0)] = _S1401;
    _S1482[int(1)] = _S1388;
    _S1482[int(2)] = _S1386;
    _S1482[int(3)] = _S1384;
    _S1482[int(4)] = _S1373;
    _S1482[int(5)] = _S1371;
    _S1482[int(6)] = _S1369;
    _S1482[int(15)] = _S1345;
    _S1482[int(8)] = _S1365;
    _S1482[int(9)] = _S1357;
    _S1482[int(10)] = _S1355;
    _S1482[int(11)] = _S1353;
    _S1482[int(12)] = _S1351;
    _S1482[int(13)] = _S1349;
    _S1482[int(14)] = _S1347;
    float3  _S1483 = _S1482[int(0)];
    float3  _S1484 = _S1482[int(1)];
    float3  _S1485 = _S1482[int(2)];
    float3  _S1486 = _S1482[int(3)];
    float3  _S1487 = _S1482[int(4)];
    float3  _S1488 = _S1482[int(5)];
    float3  _S1489 = _S1482[int(6)];
    float3  _S1490 = _S1482[int(7)];
    float3  _S1491 = _S1482[int(8)];
    float3  _S1492 = _S1482[int(9)];
    float3  _S1493 = _S1482[int(10)];
    float3  _S1494 = _S1482[int(11)];
    float3  _S1495 = _S1482[int(12)];
    float3  _S1496 = _S1482[int(13)];
    float3  _S1497 = _S1482[int(14)];
    float3  _S1498 = _S1482[int(15)];
    float3  _S1499 = _S1402 + _S1412.differential_0 + make_float3 (0.0f, 0.0f, _S1413);
    float _S1500 = _S1467.differential_0 + _S1469.differential_0;
    float _S1501 = _S1462.differential_0 + _S1464.differential_0;
    float _S1502 = _S1463.differential_0 + _S1465.differential_0;
    float _S1503 = _S1459.differential_0 + _S1461.differential_0;
    float _S1504 = _S1466.differential_0 + _S1468.differential_0;
    float _S1505 = _S1423.differential_0 + _S1425.differential_0;
    float _S1506 = _S1427.differential_0 + _S1429.differential_0;
    float _S1507 = _S1431.differential_0 + _S1433.differential_0;
    float _S1508 = _S1435.differential_0 + _S1437.differential_0;
    float _S1509 = _S1439.differential_0 + _S1441.differential_0;
    float _S1510 = _S1443.differential_0 + _S1445.differential_0;
    float _S1511 = _S1447.differential_0 + _S1449.differential_0;
    float _S1512 = _S1451.differential_0 + _S1453.differential_0;
    float _S1513 = _S1455.differential_0 + _S1457.differential_0;
    float2  _S1514 = _S1264 * _S1481;
    float2  _S1515 = _S1275 * _S1481;
    float _S1516 = _S1514.x + _S1514.y;
    if(_S1268)
    {
        float _S1517 = _S1516 / _S1269;
        float _S1518 = _S1270 * - _S1517;
        float _S1519 = _S1267 * (0.3333333432674408f * - (_S1266 * _S1517));
        k_2 = _S1519 + _S1519;
        _S1269 = _S1518;
        _S1270 = 0.0f;
    }
    else
    {
        float _S1520 = _S1516 / _S1271;
        float _S1521 = _S1267 * - _S1520;
        k_2 = _S1265 * _S1520;
        _S1269 = 0.0f;
        _S1270 = _S1521;
    }
    DiffPair_float_0 _S1522;
    (&_S1522)->primal_0 = _S1265;
    (&_S1522)->differential_0 = 0.0f;
    DiffPair_float_0 _S1523;
    (&_S1523)->primal_0 = _S1266;
    (&_S1523)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1522, &_S1523, k_2);
    float _S1524 = _S1523.differential_0 + _S1269;
    float _S1525 = _S1522.differential_0 + _S1270;
    float2  _S1526 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1527;
    (&_S1527)->primal_0 = _S1264;
    (&_S1527)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1527, _S1525);
    float2  _S1528 = _S1527.differential_0 + _S1515;
    float _S1529 = fx_5 * _S1506;
    float2  _S1530 = make_float2 (_S1529, fy_5 * _S1505) + make_float2 (dist_coeffs_5[int(8)] * _S1529, dist_coeffs_5[int(9)] * _S1529);
    float2  _S1531 = _S1252 * _S1530;
    float _S1532 = dist_coeffs_5[int(4)] * _S1530.y;
    float _S1533 = dist_coeffs_5[int(5)] * _S1530.x;
    float _S1534 = _S1531.x + _S1531.y;
    float _S1535 = r2_62 * _S1534;
    float _S1536 = r2_62 * _S1535;
    float _S1537 = dist_coeffs_5[int(7)] * _S1530.y + _S1532 + dist_coeffs_5[int(6)] * _S1530.x + _S1533 + _S1255 * _S1534 + _S1254 * _S1535 + _S1253 * _S1536 + dist_coeffs_5[int(3)] * (r2_62 * _S1536);
    float _S1538 = v_62 * _S1537;
    float _S1539 = u_62 * _S1537;
    float2  _S1540 = make_float2 (radial_30) * _S1530 + make_float2 (_S1113 * (v_62 * _S1530.y) + _S1257 * _S1533 + 2.0f * (u_62 * _S1533) + _S1110 * (v_62 * _S1530.x) + _S1539 + _S1539, _S1259 * _S1532 + 2.0f * (v_62 * _S1532) + _S1258 * _S1530.y + _S1256 * _S1530.x + _S1538 + _S1538);
    float3  _S1541 = make_float3 (_S1528.x, _S1528.y, _S1524);
    float2  _S1542 = _S1240 * _S1540;
    float2  _S1543 = _S1251 * _S1540;
    float _S1544 = _S1542.x + _S1542.y;
    if(_S1244)
    {
        float _S1545 = _S1544 / _S1245;
        float _S1546 = _S1246 * - _S1545;
        float _S1547 = _S1243 * (0.3333333432674408f * - (_S1242 * _S1545));
        k_2 = _S1547 + _S1547;
        _S1245 = _S1546;
        _S1246 = 0.0f;
    }
    else
    {
        float _S1548 = _S1544 / _S1247;
        float _S1549 = _S1243 * - _S1548;
        k_2 = _S1241 * _S1548;
        _S1245 = 0.0f;
        _S1246 = _S1549;
    }
    DiffPair_float_0 _S1550;
    (&_S1550)->primal_0 = _S1241;
    (&_S1550)->differential_0 = 0.0f;
    DiffPair_float_0 _S1551;
    (&_S1551)->primal_0 = _S1242;
    (&_S1551)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1550, &_S1551, k_2);
    float _S1552 = _S1551.differential_0 + _S1245;
    float _S1553 = _S1550.differential_0 + _S1246;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1554;
    (&_S1554)->primal_0 = _S1240;
    (&_S1554)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1554, _S1553);
    float2  _S1555 = _S1554.differential_0 + _S1543;
    float _S1556 = fx_5 * _S1508;
    float2  _S1557 = make_float2 (_S1556, fy_5 * _S1507) + make_float2 (dist_coeffs_5[int(8)] * _S1556, dist_coeffs_5[int(9)] * _S1556);
    float2  _S1558 = _S1228 * _S1557;
    float _S1559 = dist_coeffs_5[int(4)] * _S1557.y;
    float _S1560 = dist_coeffs_5[int(5)] * _S1557.x;
    float _S1561 = _S1558.x + _S1558.y;
    float _S1562 = r2_61 * _S1561;
    float _S1563 = r2_61 * _S1562;
    float _S1564 = dist_coeffs_5[int(7)] * _S1557.y + _S1559 + dist_coeffs_5[int(6)] * _S1557.x + _S1560 + _S1231 * _S1561 + _S1230 * _S1562 + _S1229 * _S1563 + dist_coeffs_5[int(3)] * (r2_61 * _S1563);
    float _S1565 = v_61 * _S1564;
    float _S1566 = u_61 * _S1564;
    float2  _S1567 = make_float2 (radial_29) * _S1557 + make_float2 (_S1113 * (v_61 * _S1557.y) + _S1233 * _S1560 + 2.0f * (u_61 * _S1560) + _S1110 * (v_61 * _S1557.x) + _S1566 + _S1566, _S1235 * _S1559 + 2.0f * (v_61 * _S1559) + _S1234 * _S1557.y + _S1232 * _S1557.x + _S1565 + _S1565);
    float3  _S1568 = make_float3 (_S1555.x, _S1555.y, _S1552);
    float2  _S1569 = _S1216 * _S1567;
    float2  _S1570 = _S1227 * _S1567;
    float _S1571 = _S1569.x + _S1569.y;
    if(_S1220)
    {
        float _S1572 = _S1571 / _S1221;
        float _S1573 = _S1222 * - _S1572;
        float _S1574 = _S1219 * (0.3333333432674408f * - (_S1218 * _S1572));
        k_2 = _S1574 + _S1574;
        _S1221 = _S1573;
        _S1222 = 0.0f;
    }
    else
    {
        float _S1575 = _S1571 / _S1223;
        float _S1576 = _S1219 * - _S1575;
        k_2 = _S1217 * _S1575;
        _S1221 = 0.0f;
        _S1222 = _S1576;
    }
    DiffPair_float_0 _S1577;
    (&_S1577)->primal_0 = _S1217;
    (&_S1577)->differential_0 = 0.0f;
    DiffPair_float_0 _S1578;
    (&_S1578)->primal_0 = _S1218;
    (&_S1578)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1577, &_S1578, k_2);
    float _S1579 = _S1578.differential_0 + _S1221;
    float _S1580 = _S1577.differential_0 + _S1222;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1581;
    (&_S1581)->primal_0 = _S1216;
    (&_S1581)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1581, _S1580);
    float2  _S1582 = _S1581.differential_0 + _S1570;
    float _S1583 = fx_5 * _S1510;
    float2  _S1584 = make_float2 (_S1583, fy_5 * _S1509) + make_float2 (dist_coeffs_5[int(8)] * _S1583, dist_coeffs_5[int(9)] * _S1583);
    float2  _S1585 = _S1204 * _S1584;
    float _S1586 = dist_coeffs_5[int(4)] * _S1584.y;
    float _S1587 = dist_coeffs_5[int(5)] * _S1584.x;
    float _S1588 = _S1585.x + _S1585.y;
    float _S1589 = r2_60 * _S1588;
    float _S1590 = r2_60 * _S1589;
    float _S1591 = dist_coeffs_5[int(7)] * _S1584.y + _S1586 + dist_coeffs_5[int(6)] * _S1584.x + _S1587 + _S1207 * _S1588 + _S1206 * _S1589 + _S1205 * _S1590 + dist_coeffs_5[int(3)] * (r2_60 * _S1590);
    float _S1592 = v_60 * _S1591;
    float _S1593 = u_60 * _S1591;
    float2  _S1594 = make_float2 (radial_28) * _S1584 + make_float2 (_S1113 * (v_60 * _S1584.y) + _S1209 * _S1587 + 2.0f * (u_60 * _S1587) + _S1110 * (v_60 * _S1584.x) + _S1593 + _S1593, _S1211 * _S1586 + 2.0f * (v_60 * _S1586) + _S1210 * _S1584.y + _S1208 * _S1584.x + _S1592 + _S1592);
    float3  _S1595 = make_float3 (_S1582.x, _S1582.y, _S1579);
    float2  _S1596 = _S1192 * _S1594;
    float2  _S1597 = _S1203 * _S1594;
    float _S1598 = _S1596.x + _S1596.y;
    if(_S1196)
    {
        float _S1599 = _S1598 / _S1197;
        float _S1600 = _S1198 * - _S1599;
        float _S1601 = _S1195 * (0.3333333432674408f * - (_S1194 * _S1599));
        k_2 = _S1601 + _S1601;
        _S1197 = _S1600;
        _S1198 = 0.0f;
    }
    else
    {
        float _S1602 = _S1598 / _S1199;
        float _S1603 = _S1195 * - _S1602;
        k_2 = _S1193 * _S1602;
        _S1197 = 0.0f;
        _S1198 = _S1603;
    }
    DiffPair_float_0 _S1604;
    (&_S1604)->primal_0 = _S1193;
    (&_S1604)->differential_0 = 0.0f;
    DiffPair_float_0 _S1605;
    (&_S1605)->primal_0 = _S1194;
    (&_S1605)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1604, &_S1605, k_2);
    float _S1606 = _S1605.differential_0 + _S1197;
    float _S1607 = _S1604.differential_0 + _S1198;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1608;
    (&_S1608)->primal_0 = _S1192;
    (&_S1608)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1608, _S1607);
    float2  _S1609 = _S1608.differential_0 + _S1597;
    float _S1610 = fx_5 * _S1512;
    float2  _S1611 = make_float2 (_S1610, fy_5 * _S1511) + make_float2 (dist_coeffs_5[int(8)] * _S1610, dist_coeffs_5[int(9)] * _S1610);
    float2  _S1612 = _S1180 * _S1611;
    float _S1613 = dist_coeffs_5[int(4)] * _S1611.y;
    float _S1614 = dist_coeffs_5[int(5)] * _S1611.x;
    float _S1615 = _S1612.x + _S1612.y;
    float _S1616 = r2_59 * _S1615;
    float _S1617 = r2_59 * _S1616;
    float _S1618 = dist_coeffs_5[int(7)] * _S1611.y + _S1613 + dist_coeffs_5[int(6)] * _S1611.x + _S1614 + _S1183 * _S1615 + _S1182 * _S1616 + _S1181 * _S1617 + dist_coeffs_5[int(3)] * (r2_59 * _S1617);
    float _S1619 = v_59 * _S1618;
    float _S1620 = u_59 * _S1618;
    float2  _S1621 = make_float2 (radial_27) * _S1611 + make_float2 (_S1113 * (v_59 * _S1611.y) + _S1185 * _S1614 + 2.0f * (u_59 * _S1614) + _S1110 * (v_59 * _S1611.x) + _S1620 + _S1620, _S1187 * _S1613 + 2.0f * (v_59 * _S1613) + _S1186 * _S1611.y + _S1184 * _S1611.x + _S1619 + _S1619);
    float3  _S1622 = make_float3 (_S1609.x, _S1609.y, _S1606);
    float2  _S1623 = _S1168 * _S1621;
    float2  _S1624 = _S1179 * _S1621;
    float _S1625 = _S1623.x + _S1623.y;
    if(_S1172)
    {
        float _S1626 = _S1625 / _S1173;
        float _S1627 = _S1174 * - _S1626;
        float _S1628 = _S1171 * (0.3333333432674408f * - (_S1170 * _S1626));
        k_2 = _S1628 + _S1628;
        _S1173 = _S1627;
        _S1174 = 0.0f;
    }
    else
    {
        float _S1629 = _S1625 / _S1175;
        float _S1630 = _S1171 * - _S1629;
        k_2 = _S1169 * _S1629;
        _S1173 = 0.0f;
        _S1174 = _S1630;
    }
    DiffPair_float_0 _S1631;
    (&_S1631)->primal_0 = _S1169;
    (&_S1631)->differential_0 = 0.0f;
    DiffPair_float_0 _S1632;
    (&_S1632)->primal_0 = _S1170;
    (&_S1632)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1631, &_S1632, k_2);
    float _S1633 = _S1632.differential_0 + _S1173;
    float _S1634 = _S1631.differential_0 + _S1174;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1635;
    (&_S1635)->primal_0 = _S1168;
    (&_S1635)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1635, _S1634);
    float2  _S1636 = _S1635.differential_0 + _S1624;
    float _S1637 = fx_5 * _S1503;
    float2  _S1638 = make_float2 (_S1637, fy_5 * _S1513) + make_float2 (dist_coeffs_5[int(8)] * _S1637, dist_coeffs_5[int(9)] * _S1637);
    float2  _S1639 = _S1156 * _S1638;
    float _S1640 = dist_coeffs_5[int(4)] * _S1638.y;
    float _S1641 = dist_coeffs_5[int(5)] * _S1638.x;
    float _S1642 = _S1639.x + _S1639.y;
    float _S1643 = r2_58 * _S1642;
    float _S1644 = r2_58 * _S1643;
    float _S1645 = dist_coeffs_5[int(7)] * _S1638.y + _S1640 + dist_coeffs_5[int(6)] * _S1638.x + _S1641 + _S1159 * _S1642 + _S1158 * _S1643 + _S1157 * _S1644 + dist_coeffs_5[int(3)] * (r2_58 * _S1644);
    float _S1646 = v_58 * _S1645;
    float _S1647 = u_58 * _S1645;
    float2  _S1648 = make_float2 (radial_26) * _S1638 + make_float2 (_S1113 * (v_58 * _S1638.y) + _S1161 * _S1641 + 2.0f * (u_58 * _S1641) + _S1110 * (v_58 * _S1638.x) + _S1647 + _S1647, _S1163 * _S1640 + 2.0f * (v_58 * _S1640) + _S1162 * _S1638.y + _S1160 * _S1638.x + _S1646 + _S1646);
    float3  _S1649 = make_float3 (_S1636.x, _S1636.y, _S1633);
    float2  _S1650 = _S1144 * _S1648;
    float2  _S1651 = _S1155 * _S1648;
    float _S1652 = _S1650.x + _S1650.y;
    if(_S1148)
    {
        float _S1653 = _S1652 / _S1149;
        float _S1654 = _S1150 * - _S1653;
        float _S1655 = _S1147 * (0.3333333432674408f * - (_S1146 * _S1653));
        k_2 = _S1655 + _S1655;
        _S1149 = _S1654;
        _S1150 = 0.0f;
    }
    else
    {
        float _S1656 = _S1652 / _S1151;
        float _S1657 = _S1147 * - _S1656;
        k_2 = _S1145 * _S1656;
        _S1149 = 0.0f;
        _S1150 = _S1657;
    }
    DiffPair_float_0 _S1658;
    (&_S1658)->primal_0 = _S1145;
    (&_S1658)->differential_0 = 0.0f;
    DiffPair_float_0 _S1659;
    (&_S1659)->primal_0 = _S1146;
    (&_S1659)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1658, &_S1659, k_2);
    float _S1660 = _S1659.differential_0 + _S1149;
    float _S1661 = _S1658.differential_0 + _S1150;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1662;
    (&_S1662)->primal_0 = _S1144;
    (&_S1662)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1662, _S1661);
    float2  _S1663 = _S1662.differential_0 + _S1651;
    float _S1664 = fx_5 * _S1500;
    float2  _S1665 = make_float2 (_S1664, fy_5 * _S1502) + make_float2 (dist_coeffs_5[int(8)] * _S1664, dist_coeffs_5[int(9)] * _S1664);
    float2  _S1666 = _S1132 * _S1665;
    float _S1667 = dist_coeffs_5[int(4)] * _S1665.y;
    float _S1668 = dist_coeffs_5[int(5)] * _S1665.x;
    float _S1669 = _S1666.x + _S1666.y;
    float _S1670 = r2_57 * _S1669;
    float _S1671 = r2_57 * _S1670;
    float _S1672 = dist_coeffs_5[int(7)] * _S1665.y + _S1667 + dist_coeffs_5[int(6)] * _S1665.x + _S1668 + _S1135 * _S1669 + _S1134 * _S1670 + _S1133 * _S1671 + dist_coeffs_5[int(3)] * (r2_57 * _S1671);
    float _S1673 = v_57 * _S1672;
    float _S1674 = u_57 * _S1672;
    float2  _S1675 = make_float2 (radial_25) * _S1665 + make_float2 (_S1113 * (v_57 * _S1665.y) + _S1137 * _S1668 + 2.0f * (u_57 * _S1668) + _S1110 * (v_57 * _S1665.x) + _S1674 + _S1674, _S1139 * _S1667 + 2.0f * (v_57 * _S1667) + _S1138 * _S1665.y + _S1136 * _S1665.x + _S1673 + _S1673);
    float3  _S1676 = make_float3 (_S1663.x, _S1663.y, _S1660);
    float2  _S1677 = _S1120 * _S1675;
    float2  _S1678 = _S1131 * _S1675;
    float _S1679 = _S1677.x + _S1677.y;
    if(_S1124)
    {
        float _S1680 = _S1679 / _S1125;
        float _S1681 = _S1126 * - _S1680;
        float _S1682 = _S1123 * (0.3333333432674408f * - (_S1122 * _S1680));
        k_2 = _S1682 + _S1682;
        _S1125 = _S1681;
        _S1126 = 0.0f;
    }
    else
    {
        float _S1683 = _S1679 / _S1127;
        float _S1684 = _S1123 * - _S1683;
        k_2 = _S1121 * _S1683;
        _S1125 = 0.0f;
        _S1126 = _S1684;
    }
    DiffPair_float_0 _S1685;
    (&_S1685)->primal_0 = _S1121;
    (&_S1685)->differential_0 = 0.0f;
    DiffPair_float_0 _S1686;
    (&_S1686)->primal_0 = _S1122;
    (&_S1686)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1685, &_S1686, k_2);
    float _S1687 = _S1686.differential_0 + _S1125;
    float _S1688 = _S1685.differential_0 + _S1126;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1689;
    (&_S1689)->primal_0 = _S1120;
    (&_S1689)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1689, _S1688);
    float2  _S1690 = _S1689.differential_0 + _S1678;
    float _S1691 = fx_5 * _S1504;
    float2  _S1692 = make_float2 (_S1691, fy_5 * _S1501) + make_float2 (dist_coeffs_5[int(8)] * _S1691, dist_coeffs_5[int(9)] * _S1691);
    float2  _S1693 = _S1106 * _S1692;
    float _S1694 = dist_coeffs_5[int(4)] * _S1692.y;
    float _S1695 = dist_coeffs_5[int(5)] * _S1692.x;
    float _S1696 = _S1693.x + _S1693.y;
    float _S1697 = r2_56 * _S1696;
    float _S1698 = r2_56 * _S1697;
    float _S1699 = dist_coeffs_5[int(7)] * _S1692.y + _S1694 + dist_coeffs_5[int(6)] * _S1692.x + _S1695 + _S1109 * _S1696 + _S1108 * _S1697 + _S1107 * _S1698 + dist_coeffs_5[int(3)] * (r2_56 * _S1698);
    float _S1700 = v_56 * _S1699;
    float _S1701 = u_56 * _S1699;
    float2  _S1702 = make_float2 (radial_24) * _S1692 + make_float2 (_S1113 * (v_56 * _S1692.y) + _S1112 * _S1695 + 2.0f * (u_56 * _S1695) + _S1110 * (v_56 * _S1692.x) + _S1701 + _S1701, _S1115 * _S1694 + 2.0f * (v_56 * _S1694) + _S1114 * _S1692.y + _S1111 * _S1692.x + _S1700 + _S1700);
    float3  _S1703 = make_float3 (_S1690.x, _S1690.y, _S1687);
    float2  _S1704 = _S1094 * _S1702;
    float2  _S1705 = _S1105 * _S1702;
    float _S1706 = _S1704.x + _S1704.y;
    if(_S1098)
    {
        float _S1707 = _S1706 / _S1099;
        float _S1708 = _S1100 * - _S1707;
        float _S1709 = _S1097 * (0.3333333432674408f * - (_S1096 * _S1707));
        k_2 = _S1709 + _S1709;
        _S1099 = _S1708;
        _S1100 = 0.0f;
    }
    else
    {
        float _S1710 = _S1706 / _S1101;
        float _S1711 = _S1097 * - _S1710;
        k_2 = _S1095 * _S1710;
        _S1099 = 0.0f;
        _S1100 = _S1711;
    }
    DiffPair_float_0 _S1712;
    (&_S1712)->primal_0 = _S1095;
    (&_S1712)->differential_0 = 0.0f;
    DiffPair_float_0 _S1713;
    (&_S1713)->primal_0 = _S1096;
    (&_S1713)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1712, &_S1713, k_2);
    float _S1714 = _S1713.differential_0 + _S1099;
    float _S1715 = _S1712.differential_0 + _S1100;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1716;
    (&_S1716)->primal_0 = _S1094;
    (&_S1716)->differential_0 = _S1526;
    s_bwd_length_impl_1(&_S1716, _S1715);
    float2  _S1717 = _S1716.differential_0 + _S1705;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1718;
    (&_S1718)->primal_0 = R_5;
    (&_S1718)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1719;
    (&_S1719)->primal_0 = _S1093;
    (&_S1719)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1718, &_S1719, _S1499);
    DiffPair_float_0 _S1720;
    (&_S1720)->primal_0 = _S1090;
    (&_S1720)->differential_0 = 0.0f;
    DiffPair_float_0 _S1721;
    (&_S1721)->primal_0 = _S1092;
    (&_S1721)->differential_0 = 0.0f;
    _d_max_0(&_S1720, &_S1721, 0.0f);
    DiffPair_float_0 _S1722;
    (&_S1722)->primal_0 = _S1089;
    (&_S1722)->differential_0 = 0.0f;
    DiffPair_float_0 _S1723;
    (&_S1723)->primal_0 = _S1092;
    (&_S1723)->differential_0 = 0.0f;
    _d_min_0(&_S1722, &_S1723, 0.0f);
    float _S1724 = _S1721.differential_0 + _S1723.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1725;
    (&_S1725)->primal_0 = _S1091;
    (&_S1725)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1725, _S1724);
    float3  _S1726 = _S1725.differential_0 + _S1541;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1727;
    (&_S1727)->primal_0 = R_5;
    (&_S1727)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1728;
    (&_S1728)->primal_0 = pos_i_13;
    (&_S1728)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1727, &_S1728, _S1726);
    DiffPair_float_0 _S1729;
    (&_S1729)->primal_0 = _S1086;
    (&_S1729)->differential_0 = 0.0f;
    DiffPair_float_0 _S1730;
    (&_S1730)->primal_0 = _S1088;
    (&_S1730)->differential_0 = 0.0f;
    _d_max_0(&_S1729, &_S1730, _S1720.differential_0);
    DiffPair_float_0 _S1731;
    (&_S1731)->primal_0 = _S1085;
    (&_S1731)->differential_0 = 0.0f;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = _S1088;
    (&_S1732)->differential_0 = 0.0f;
    _d_min_0(&_S1731, &_S1732, _S1722.differential_0);
    float _S1733 = _S1730.differential_0 + _S1732.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1734;
    (&_S1734)->primal_0 = _S1087;
    (&_S1734)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1734, _S1733);
    float3  _S1735 = _S1734.differential_0 + _S1568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1736;
    (&_S1736)->primal_0 = R_5;
    (&_S1736)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1737;
    (&_S1737)->primal_0 = pos_i_12;
    (&_S1737)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1736, &_S1737, _S1735);
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = _S1082;
    (&_S1738)->differential_0 = 0.0f;
    DiffPair_float_0 _S1739;
    (&_S1739)->primal_0 = _S1084;
    (&_S1739)->differential_0 = 0.0f;
    _d_max_0(&_S1738, &_S1739, _S1729.differential_0);
    DiffPair_float_0 _S1740;
    (&_S1740)->primal_0 = _S1081;
    (&_S1740)->differential_0 = 0.0f;
    DiffPair_float_0 _S1741;
    (&_S1741)->primal_0 = _S1084;
    (&_S1741)->differential_0 = 0.0f;
    _d_min_0(&_S1740, &_S1741, _S1731.differential_0);
    float _S1742 = _S1739.differential_0 + _S1741.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1743;
    (&_S1743)->primal_0 = _S1083;
    (&_S1743)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1743, _S1742);
    float3  _S1744 = _S1743.differential_0 + _S1595;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1745;
    (&_S1745)->primal_0 = R_5;
    (&_S1745)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1746;
    (&_S1746)->primal_0 = pos_i_11;
    (&_S1746)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1745, &_S1746, _S1744);
    DiffPair_float_0 _S1747;
    (&_S1747)->primal_0 = _S1078;
    (&_S1747)->differential_0 = 0.0f;
    DiffPair_float_0 _S1748;
    (&_S1748)->primal_0 = _S1080;
    (&_S1748)->differential_0 = 0.0f;
    _d_max_0(&_S1747, &_S1748, _S1738.differential_0);
    DiffPair_float_0 _S1749;
    (&_S1749)->primal_0 = _S1077;
    (&_S1749)->differential_0 = 0.0f;
    DiffPair_float_0 _S1750;
    (&_S1750)->primal_0 = _S1080;
    (&_S1750)->differential_0 = 0.0f;
    _d_min_0(&_S1749, &_S1750, _S1740.differential_0);
    float _S1751 = _S1748.differential_0 + _S1750.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1752;
    (&_S1752)->primal_0 = _S1079;
    (&_S1752)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1752, _S1751);
    float3  _S1753 = _S1752.differential_0 + _S1622;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1754;
    (&_S1754)->primal_0 = R_5;
    (&_S1754)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1755;
    (&_S1755)->primal_0 = pos_i_10;
    (&_S1755)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1754, &_S1755, _S1753);
    DiffPair_float_0 _S1756;
    (&_S1756)->primal_0 = _S1074;
    (&_S1756)->differential_0 = 0.0f;
    DiffPair_float_0 _S1757;
    (&_S1757)->primal_0 = _S1076;
    (&_S1757)->differential_0 = 0.0f;
    _d_max_0(&_S1756, &_S1757, _S1747.differential_0);
    DiffPair_float_0 _S1758;
    (&_S1758)->primal_0 = _S1073;
    (&_S1758)->differential_0 = 0.0f;
    DiffPair_float_0 _S1759;
    (&_S1759)->primal_0 = _S1076;
    (&_S1759)->differential_0 = 0.0f;
    _d_min_0(&_S1758, &_S1759, _S1749.differential_0);
    float _S1760 = _S1757.differential_0 + _S1759.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1761;
    (&_S1761)->primal_0 = _S1075;
    (&_S1761)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1761, _S1760);
    float3  _S1762 = _S1761.differential_0 + _S1649;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1763;
    (&_S1763)->primal_0 = R_5;
    (&_S1763)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1764;
    (&_S1764)->primal_0 = pos_i_9;
    (&_S1764)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1763, &_S1764, _S1762);
    DiffPair_float_0 _S1765;
    (&_S1765)->primal_0 = _S1070;
    (&_S1765)->differential_0 = 0.0f;
    DiffPair_float_0 _S1766;
    (&_S1766)->primal_0 = _S1072;
    (&_S1766)->differential_0 = 0.0f;
    _d_max_0(&_S1765, &_S1766, _S1756.differential_0);
    DiffPair_float_0 _S1767;
    (&_S1767)->primal_0 = _S1069;
    (&_S1767)->differential_0 = 0.0f;
    DiffPair_float_0 _S1768;
    (&_S1768)->primal_0 = _S1072;
    (&_S1768)->differential_0 = 0.0f;
    _d_min_0(&_S1767, &_S1768, _S1758.differential_0);
    float _S1769 = _S1766.differential_0 + _S1768.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1770;
    (&_S1770)->primal_0 = _S1071;
    (&_S1770)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1770, _S1769);
    float3  _S1771 = _S1770.differential_0 + _S1676;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1772;
    (&_S1772)->primal_0 = R_5;
    (&_S1772)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1773;
    (&_S1773)->primal_0 = pos_i_8;
    (&_S1773)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1772, &_S1773, _S1771);
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = _S1066;
    (&_S1774)->differential_0 = 0.0f;
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = _S1068;
    (&_S1775)->differential_0 = 0.0f;
    _d_max_0(&_S1774, &_S1775, _S1765.differential_0);
    DiffPair_float_0 _S1776;
    (&_S1776)->primal_0 = _S1065;
    (&_S1776)->differential_0 = 0.0f;
    DiffPair_float_0 _S1777;
    (&_S1777)->primal_0 = _S1068;
    (&_S1777)->differential_0 = 0.0f;
    _d_min_0(&_S1776, &_S1777, _S1767.differential_0);
    float _S1778 = _S1775.differential_0 + _S1777.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1779;
    (&_S1779)->primal_0 = _S1067;
    (&_S1779)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1779, _S1778);
    float3  _S1780 = _S1779.differential_0 + _S1703;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1781;
    (&_S1781)->primal_0 = R_5;
    (&_S1781)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1782;
    (&_S1782)->primal_0 = pos_i_7;
    (&_S1782)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1781, &_S1782, _S1780);
    DiffPair_float_0 _S1783;
    (&_S1783)->primal_0 = 0.0f;
    (&_S1783)->differential_0 = 0.0f;
    DiffPair_float_0 _S1784;
    (&_S1784)->primal_0 = _S1064;
    (&_S1784)->differential_0 = 0.0f;
    _d_max_0(&_S1783, &_S1784, _S1774.differential_0);
    DiffPair_float_0 _S1785;
    (&_S1785)->primal_0 = 1.00000001504746622e+30f;
    (&_S1785)->differential_0 = 0.0f;
    DiffPair_float_0 _S1786;
    (&_S1786)->primal_0 = _S1064;
    (&_S1786)->differential_0 = 0.0f;
    _d_min_0(&_S1785, &_S1786, _S1776.differential_0);
    float _S1787 = _S1784.differential_0 + _S1786.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1788;
    (&_S1788)->primal_0 = _S1063;
    (&_S1788)->differential_0 = _S1342;
    s_bwd_length_impl_0(&_S1788, _S1787);
    float3  _S1789 = _S1788.differential_0 + make_float3 (_S1717.x, _S1717.y, _S1714);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1790;
    (&_S1790)->primal_0 = R_5;
    (&_S1790)->differential_0 = _S1404;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1791;
    (&_S1791)->primal_0 = pos_5;
    (&_S1791)->differential_0 = _S1342;
    s_bwd_prop_mul_0(&_S1790, &_S1791, _S1789);
    float3  _S1792 = _S1499 + _S1726 + _S1735 + _S1744 + _S1753 + _S1762 + _S1771 + _S1780 + _S1789 + _S1407.differential_0;
    Matrix<float, 3, 3>  _S1793 = _S1718.differential_0 + _S1727.differential_0 + _S1736.differential_0 + _S1745.differential_0 + _S1754.differential_0 + _S1763.differential_0 + _S1772.differential_0 + _S1781.differential_0 + _S1790.differential_0 + _S1408;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_1)[int(0)] = _S1483;
    (*v_sh_coeffs_1)[int(1)] = _S1484;
    (*v_sh_coeffs_1)[int(2)] = _S1485;
    (*v_sh_coeffs_1)[int(3)] = _S1486;
    (*v_sh_coeffs_1)[int(4)] = _S1487;
    (*v_sh_coeffs_1)[int(5)] = _S1488;
    (*v_sh_coeffs_1)[int(6)] = _S1489;
    (*v_sh_coeffs_1)[int(7)] = _S1490;
    (*v_sh_coeffs_1)[int(8)] = _S1491;
    (*v_sh_coeffs_1)[int(9)] = _S1492;
    (*v_sh_coeffs_1)[int(10)] = _S1493;
    (*v_sh_coeffs_1)[int(11)] = _S1494;
    (*v_sh_coeffs_1)[int(12)] = _S1495;
    (*v_sh_coeffs_1)[int(13)] = _S1496;
    (*v_sh_coeffs_1)[int(14)] = _S1497;
    (*v_sh_coeffs_1)[int(15)] = _S1498;
    *v_R_1 = _S1793;
    *v_t_1 = _S1792;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  dOut_7)
{
    float3  _S1794 = _slang_select(((*dpx_8).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_8).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_7;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S1794;
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
    float3  _S1795 = - (m_1 * (ray_o_0 - center_0));
    float3  ta_0 = _S1795 - k_3;
    float3  tb_0 = _S1795 + k_3;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S1796 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S1796;
    return (*t0_0) < _S1796;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_9, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_8)
{
    float _S1797 = (1.0f - (*dps_0).primal_0) * dOut_8;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S1797;
    DiffPair_float_0 _S1798 = *dpy_4;
    float _S1799 = (*dps_0).primal_0 * dOut_8;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S1799;
    float _S1800 = (_S1798.primal_0 - (*dpx_9).primal_0) * dOut_8;
    dps_0->primal_0 = _S1798.primal_0;
    dps_0->differential_0 = _S1800;
    return;
}

inline __device__ float lerp_0(float x_14, float y_9, float s_0)
{
    return x_14 + (y_9 - x_14) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_0)
{
    float _S1801 = w_0.z;
    float _S1802 = 1.0f - _S1801;
    float _S1803 = w_0.y;
    float _S1804 = 1.0f - _S1803;
    float _S1805 = _S1802 * _S1804;
    float _S1806 = w_0.x;
    float _S1807 = 1.0f - _S1806;
    float _S1808 = _S1802 * _S1803;
    float _S1809 = _S1801 * _S1804;
    float _S1810 = _S1801 * _S1803;
    return _S1805 * _S1807 * (*densities_6)[int(0)] + _S1805 * _S1806 * (*densities_6)[int(1)] + _S1808 * _S1807 * (*densities_6)[int(2)] + _S1808 * _S1806 * (*densities_6)[int(3)] + _S1809 * _S1807 * (*densities_6)[int(4)] + _S1809 * _S1806 * (*densities_6)[int(5)] + _S1810 * _S1807 * (*densities_6)[int(6)] + _S1810 * _S1806 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  densities_7, float3  ray_o_1, float3  ray_d_1)
{
    float _S1811 = 0.5f * size_6;
    float3  m_2 = make_float3 (1.0f) / ray_d_1;
    float3  k_4 = abs_0(m_2) * make_float3 (_S1811);
    float3  _S1812 = - (m_2 * (ray_o_1 - (pos_6 + make_float3 (_S1811))));
    float3  ta_1 = _S1812 - k_4;
    float3  tb_1 = _S1812 + k_4;
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
        float3  _S1813 = (ray_o_1 + ray_d_1 * make_float3 (lerp_0(t0_1, t1_1, (float(i_5) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6);
        FixedArray<float, 8>  _S1814 = densities_7;
        float _S1815 = interp_0(&_S1814, _S1813);
        float _S1816;
        if(_S1815 > 1.10000002384185791f)
        {
            _S1816 = _S1815;
        }
        else
        {
            _S1816 = (F32_exp((0.90909093618392944f * _S1815 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S1816;
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
    float _S1817;
};

inline __device__ float3  s_primal_ctx_abs_0(float3  _S1818)
{
    return abs_0(_S1818);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1819, float _S1820, float _S1821)
{
    return lerp_0(_S1819, _S1820, _S1821);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S1822 = dpw_0.z;
    float _S1823 = 1.0f - _S1822;
    float _S1824 = dpw_0.y;
    float _S1825 = 1.0f - _S1824;
    float _S1826 = _S1823 * _S1825;
    float _S1827 = dpw_0.x;
    float _S1828 = 1.0f - _S1827;
    float _S1829 = _S1823 * _S1824;
    float _S1830 = _S1822 * _S1825;
    float _S1831 = _S1822 * _S1824;
    return _S1826 * _S1828 * (*dpdensities_0)[int(0)] + _S1826 * _S1827 * (*dpdensities_0)[int(1)] + _S1829 * _S1828 * (*dpdensities_0)[int(2)] + _S1829 * _S1827 * (*dpdensities_0)[int(3)] + _S1830 * _S1828 * (*dpdensities_0)[int(4)] + _S1830 * _S1827 * (*dpdensities_0)[int(5)] + _S1831 * _S1828 * (*dpdensities_0)[int(6)] + _S1831 * _S1827 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_exp_0(float _S1832)
{
    return (F32_exp((_S1832)));
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_0, float3  dpray_d_0, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S1817 = 0.0f;
    _s_diff_ctx_0->_S1817 = 0.0f;
    float _S1833 = 0.5f * size_7;
    float3  m_3 = make_float3 (1.0f) / dpray_d_0;
    float3  k_5 = s_primal_ctx_abs_0(m_3) * make_float3 (_S1833);
    float3  _S1834 = - (m_3 * (dpray_o_0 - (pos_7 + make_float3 (_S1833))));
    float3  ta_2 = _S1834 - k_5;
    float3  tb_2 = _S1834 + k_5;
    float _S1835 = (F32_max(((F32_max((ta_2.x), (ta_2.y)))), ((F32_max((ta_2.z), (0.0f))))));
    float _S1836 = (F32_min(((F32_min((tb_2.x), (tb_2.y)))), (tb_2.z)));
    float accum_2;
    if(!!(_S1835 < _S1836))
    {
        float _S1837 = - (_S1836 - _S1835);
        bool _runFlag_0 = true;
        int i_6 = int(0);
        accum_2 = 0.0f;
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_0->_S1817 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            int _S1838;
            float _S1839;
            if(i_6 < int(8))
            {
                float _S1840 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_0 + dpray_d_0 * make_float3 (s_primal_ctx_lerp_0(_S1835, _S1836, (float(i_6) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S1840 > 1.10000002384185791f)
                {
                    _S1839 = _S1840;
                }
                else
                {
                    _S1839 = s_primal_ctx_exp_0(0.90909093618392944f * _S1840 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S1839;
                _S1838 = int(2);
                _S1839 = accum_3;
            }
            else
            {
                _S1838 = int(1);
                _S1839 = 0.0f;
            }
            if(_S1838 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                i_6 = i_6 + int(1);
                accum_2 = _S1839;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = (F32_min((1.0f - s_primal_ctx_exp_0(_S1837 / 8.0f * accum_2)), (0.99900001287460327f)));
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1841, float _S1842)
{
    _d_exp_0(_S1841, _S1842);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_2)
{
    float _S1843 = (*dpw_1).primal_0.z;
    float _S1844 = 1.0f - _S1843;
    float _S1845 = (*dpw_1).primal_0.y;
    float _S1846 = 1.0f - _S1845;
    float _S1847 = _S1844 * _S1846;
    float _S1848 = (*dpw_1).primal_0.x;
    float _S1849 = 1.0f - _S1848;
    float _S1850 = _S1844 * _S1845;
    float _S1851 = _S1843 * _S1846;
    float _S1852 = _S1843 * _S1845;
    float _S1853 = _S1852 * _S1848 * _s_dOut_2;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_2;
    float _S1854 = _S1852 * _S1849 * _s_dOut_2;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_2;
    float _S1855 = _S1851 * _S1848 * _s_dOut_2;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_2;
    float _S1856 = _S1851 * _S1849 * _s_dOut_2;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_2;
    float _S1857 = _S1850 * _S1848 * _s_dOut_2;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_2;
    float _S1858 = _S1850 * _S1849 * _s_dOut_2;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_2;
    float _S1859 = _S1847 * _S1848 * _s_dOut_2;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_2;
    float _S1860 = _S1847 * _S1849 * _s_dOut_2;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_2;
    float _S1861 = _S1848 * s_diff_w7_T_0 + _S1849 * s_diff_w6_T_0;
    float _S1862 = _S1848 * s_diff_w5_T_0 + _S1849 * s_diff_w4_T_0;
    float _S1863 = _S1848 * s_diff_w3_T_0 + _S1849 * s_diff_w2_T_0;
    float _S1864 = _S1848 * s_diff_w1_T_0 + _S1849 * s_diff_w0_T_0;
    float3  _S1865 = make_float3 (_S1852 * s_diff_w7_T_0 + _S1851 * s_diff_w5_T_0 + _S1850 * s_diff_w3_T_0 + _S1847 * s_diff_w1_T_0 + - (_S1852 * s_diff_w6_T_0 + _S1851 * s_diff_w4_T_0 + _S1850 * s_diff_w2_T_0 + _S1847 * s_diff_w0_T_0), _S1843 * _S1861 + _S1844 * _S1863 + - (_S1843 * _S1862 + _S1844 * _S1864), _S1845 * _S1861 + _S1846 * _S1862 + - (_S1845 * _S1863 + _S1846 * _S1864));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S1865;
    FixedArray<float, 8>  _S1866;
    _S1866[int(0)] = 0.0f;
    _S1866[int(1)] = 0.0f;
    _S1866[int(2)] = 0.0f;
    _S1866[int(3)] = 0.0f;
    _S1866[int(4)] = 0.0f;
    _S1866[int(5)] = 0.0f;
    _S1866[int(6)] = 0.0f;
    _S1866[int(7)] = 0.0f;
    _S1866[int(7)] = _S1853;
    _S1866[int(6)] = _S1854;
    _S1866[int(5)] = _S1855;
    _S1866[int(4)] = _S1856;
    _S1866[int(3)] = _S1857;
    _S1866[int(2)] = _S1858;
    _S1866[int(1)] = _S1859;
    _S1866[int(0)] = _S1860;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S1866;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S1867, DiffPair_float_0 * _S1868, DiffPair_float_0 * _S1869, float _S1870)
{
    _d_lerp_0(_S1867, _S1868, _S1869, _S1870);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1871, float3  _S1872)
{
    _d_abs_vector_0(_S1871, _S1872);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float _s_dOut_3, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_1)
{
    FixedArray<float, 8>  _S1873 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1874 = *dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1875 = *dpray_d_1;
    float _S1876 = 0.5f * size_8;
    float3  _S1877 = make_float3 (_S1876);
    float3  _S1878 = make_float3 (size_8);
    float3  m_4 = make_float3 (1.0f) / (*dpray_d_1).primal_0;
    float3  _S1879 = (*dpray_d_1).primal_0 * (*dpray_d_1).primal_0;
    float3  _S1880 = (*dpray_o_1).primal_0 - (pos_8 + make_float3 (_S1876));
    float3  k_6 = s_primal_ctx_abs_0(m_4) * make_float3 (_S1876);
    float3  _S1881 = - (m_4 * _S1880);
    float3  ta_3 = _S1881 - k_6;
    float3  tb_3 = _S1881 + k_6;
    float _S1882 = ta_3.x;
    float _S1883 = ta_3.y;
    float _S1884 = (F32_max((_S1882), (_S1883)));
    float _S1885 = ta_3.z;
    float _S1886 = (F32_max((_S1885), (0.0f)));
    float _S1887 = (F32_max((_S1884), (_S1886)));
    float _S1888 = tb_3.x;
    float _S1889 = tb_3.y;
    float _S1890 = (F32_min((_S1888), (_S1889)));
    float _S1891 = tb_3.z;
    float _S1892 = (F32_min((_S1890), (_S1891)));
    bool _S1893 = !!(_S1887 < _S1892);
    float _S1894;
    float _S1895;
    float _S1896;
    if(_S1893)
    {
        float _S1897 = - (_S1892 - _S1887) / 8.0f;
        float _S1898 = _S1897 * _s_diff_ctx_1->_S1817;
        _S1894 = 1.0f - s_primal_ctx_exp_0(_S1898);
        _S1895 = _S1898;
        _S1896 = _S1897;
    }
    else
    {
        _S1894 = 0.0f;
        _S1895 = 0.0f;
        _S1896 = 0.0f;
    }
    float3  _S1899 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1900 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<float, 8>  _S1901;
    float3  _S1902;
    float3  _S1903;
    if(_S1893)
    {
        DiffPair_float_0 _S1904;
        (&_S1904)->primal_0 = _S1894;
        (&_S1904)->differential_0 = 0.0f;
        DiffPair_float_0 _S1905;
        (&_S1905)->primal_0 = 0.99900001287460327f;
        (&_S1905)->differential_0 = 0.0f;
        _d_min_0(&_S1904, &_S1905, _s_dOut_3);
        float _S1906 = - _S1904.differential_0;
        DiffPair_float_0 _S1907;
        (&_S1907)->primal_0 = _S1895;
        (&_S1907)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S1907, _S1906);
        float _S1908 = _S1896 * _S1907.differential_0;
        float _S1909 = 0.125f * (_s_diff_ctx_1->_S1817 * _S1907.differential_0);
        int _dc_0 = int(8);
        _S1894 = _S1908;
        _S1901[int(0)] = 0.0f;
        _S1901[int(1)] = 0.0f;
        _S1901[int(2)] = 0.0f;
        _S1901[int(3)] = 0.0f;
        _S1901[int(4)] = 0.0f;
        _S1901[int(5)] = 0.0f;
        _S1901[int(6)] = 0.0f;
        _S1901[int(7)] = 0.0f;
        _S1902 = _S1899;
        _S1903 = _S1899;
        _S1895 = 0.0f;
        _S1896 = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S1910 = _dc_0 < int(8);
            float _S1911;
            float _S1912;
            int _S1913;
            float3  _S1914;
            float3  _S1915;
            bool _S1916;
            if(_S1910)
            {
                float _S1917 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S1918 = s_primal_ctx_lerp_0(_S1887, _S1892, _S1917);
                float3  _S1919 = make_float3 (_S1918);
                float3  _S1920 = (_S1874.primal_0 + _S1875.primal_0 * make_float3 (_S1918) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S1921 = _S1873;
                float _S1922 = s_primal_ctx_interp_0(&_S1921, _S1920);
                bool _S1923 = _S1922 > 1.10000002384185791f;
                if(_S1923)
                {
                    _S1911 = 0.0f;
                }
                else
                {
                    _S1911 = 0.90909093618392944f * _S1922 - 0.90468984842300415f;
                }
                _S1913 = int(2);
                _S1916 = _S1923;
                _S1914 = _S1920;
                _S1915 = _S1919;
                _S1912 = _S1917;
            }
            else
            {
                _S1913 = int(1);
                _S1916 = false;
                _S1911 = 0.0f;
                _S1914 = _S1899;
                _S1915 = _S1899;
                _S1912 = 0.0f;
            }
            float _S1924;
            float _S1925;
            if(!(_S1913 != int(2)))
            {
                _S1924 = _S1894;
                _S1925 = 0.0f;
            }
            else
            {
                _S1924 = 0.0f;
                _S1925 = _S1894;
            }
            if(_S1910)
            {
                float _S1926 = _S1924 + _S1925;
                float _S1927;
                if(_S1916)
                {
                    _S1927 = _S1924;
                }
                else
                {
                    DiffPair_float_0 _S1928;
                    (&_S1928)->primal_0 = _S1911;
                    (&_S1928)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S1928, _S1924);
                    _S1927 = 0.90909093618392944f * _S1928.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S1929;
                (&_S1929)->primal_0 = _S1873;
                (&_S1929)->differential_0 = _S1900;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1930;
                (&_S1930)->primal_0 = _S1914;
                (&_S1930)->differential_0 = _S1899;
                s_bwd_prop_interp_0(&_S1929, &_S1930, _S1927);
                float3  _S1931 = _S1930.differential_0 / _S1878;
                float3  _S1932 = _S1875.primal_0 * _S1931;
                float3  _S1933 = _S1915 * _S1931;
                float _S1934 = _S1932.x + _S1932.y + _S1932.z;
                DiffPair_float_0 _S1935;
                (&_S1935)->primal_0 = _S1887;
                (&_S1935)->differential_0 = 0.0f;
                DiffPair_float_0 _S1936;
                (&_S1936)->primal_0 = _S1892;
                (&_S1936)->differential_0 = 0.0f;
                DiffPair_float_0 _S1937;
                (&_S1937)->primal_0 = _S1912;
                (&_S1937)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S1935, &_S1936, &_S1937, _S1934);
                float _S1938 = (&_S1929)->differential_0[int(0)] + _S1901[int(0)];
                float _S1939 = (&_S1929)->differential_0[int(1)] + _S1901[int(1)];
                float _S1940 = (&_S1929)->differential_0[int(2)] + _S1901[int(2)];
                float _S1941 = (&_S1929)->differential_0[int(3)] + _S1901[int(3)];
                float _S1942 = (&_S1929)->differential_0[int(4)] + _S1901[int(4)];
                float _S1943 = (&_S1929)->differential_0[int(5)] + _S1901[int(5)];
                float _S1944 = (&_S1929)->differential_0[int(6)] + _S1901[int(6)];
                float _S1945 = (&_S1929)->differential_0[int(7)] + _S1901[int(7)];
                float3  _S1946 = _S1931 + _S1902;
                float3  _S1947 = _S1933 + _S1903;
                float _S1948 = _S1936.differential_0 + _S1895;
                float _S1949 = _S1935.differential_0 + _S1896;
                _S1894 = _S1926;
                _S1901[int(0)] = _S1938;
                _S1901[int(1)] = _S1939;
                _S1901[int(2)] = _S1940;
                _S1901[int(3)] = _S1941;
                _S1901[int(4)] = _S1942;
                _S1901[int(5)] = _S1943;
                _S1901[int(6)] = _S1944;
                _S1901[int(7)] = _S1945;
                _S1902 = _S1946;
                _S1903 = _S1947;
                _S1895 = _S1948;
                _S1896 = _S1949;
            }
            else
            {
                _S1894 = _S1925;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S1950 = - _S1909;
        float _S1951 = - _S1950 + _S1896;
        float3  _S1952 = _S1902;
        _S1894 = _S1950 + _S1895;
        _S1895 = _S1951;
        _S1902 = _S1903;
        _S1903 = _S1952;
    }
    else
    {
        _S1894 = 0.0f;
        _S1895 = 0.0f;
        _S1902 = _S1899;
        _S1903 = _S1899;
        _S1901[int(0)] = 0.0f;
        _S1901[int(1)] = 0.0f;
        _S1901[int(2)] = 0.0f;
        _S1901[int(3)] = 0.0f;
        _S1901[int(4)] = 0.0f;
        _S1901[int(5)] = 0.0f;
        _S1901[int(6)] = 0.0f;
        _S1901[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S1953;
    (&_S1953)->primal_0 = _S1890;
    (&_S1953)->differential_0 = 0.0f;
    DiffPair_float_0 _S1954;
    (&_S1954)->primal_0 = _S1891;
    (&_S1954)->differential_0 = 0.0f;
    _d_min_0(&_S1953, &_S1954, _S1894);
    DiffPair_float_0 _S1955;
    (&_S1955)->primal_0 = _S1888;
    (&_S1955)->differential_0 = 0.0f;
    DiffPair_float_0 _S1956;
    (&_S1956)->primal_0 = _S1889;
    (&_S1956)->differential_0 = 0.0f;
    _d_min_0(&_S1955, &_S1956, _S1953.differential_0);
    DiffPair_float_0 _S1957;
    (&_S1957)->primal_0 = _S1884;
    (&_S1957)->differential_0 = 0.0f;
    DiffPair_float_0 _S1958;
    (&_S1958)->primal_0 = _S1886;
    (&_S1958)->differential_0 = 0.0f;
    _d_max_0(&_S1957, &_S1958, _S1895);
    DiffPair_float_0 _S1959;
    (&_S1959)->primal_0 = _S1885;
    (&_S1959)->differential_0 = 0.0f;
    DiffPair_float_0 _S1960;
    (&_S1960)->primal_0 = 0.0f;
    (&_S1960)->differential_0 = 0.0f;
    _d_max_0(&_S1959, &_S1960, _S1958.differential_0);
    DiffPair_float_0 _S1961;
    (&_S1961)->primal_0 = _S1882;
    (&_S1961)->differential_0 = 0.0f;
    DiffPair_float_0 _S1962;
    (&_S1962)->primal_0 = _S1883;
    (&_S1962)->differential_0 = 0.0f;
    _d_max_0(&_S1961, &_S1962, _S1957.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S1955.differential_0, _S1956.differential_0, _S1954.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S1961.differential_0, _S1962.differential_0, _S1959.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S1963 = _S1877 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1964;
    (&_S1964)->primal_0 = m_4;
    (&_S1964)->differential_0 = _S1899;
    s_bwd_prop_abs_0(&_S1964, _S1963);
    float3  _S1965 = m_4 * s_diff_n_T_0;
    float3  _S1966 = - ((_S1964.differential_0 + _S1880 * s_diff_n_T_0) / _S1879) + _S1902;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1966;
    float3  _S1967 = _S1965 + _S1903;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1967;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S1901;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S1968, float _S1969, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S1970, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1971, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1972, float _S1973)
{
    FixedArray<float, 8>  _S1974 = _S1970->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1975;
    float _S1976 = s_primal_ctx_evaluate_alpha_voxel_0(_S1968, _S1969, &_S1974, (*_S1971).primal_0, (*_S1972).primal_0, &_S1975);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1977 = _S1975;
    s_bwd_prop_evaluate_alpha_voxel_0(_S1968, _S1969, _S1970, _S1971, _S1972, _S1973, &_S1977);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  densities_8, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    FixedArray<float, 8>  _S1978 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = densities_8;
    (&dp_densities_0)->differential_0 = _S1978;
    float3  _S1979 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1979;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1979;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_0 = dp_ray_o_0.differential_0;
    *v_ray_d_0 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  densities_9, float3  rgb_0, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_4)
{
    *out_rgb_0 = rgb_0;
    float _S1980 = 0.5f * size_10;
    float3  m_5 = make_float3 (1.0f) / ray_d_3;
    float3  k_7 = abs_0(m_5) * make_float3 (_S1980);
    float3  _S1981 = - (m_5 * (ray_o_3 - (pos_10 + make_float3 (_S1980))));
    float3  ta_4 = _S1981 - k_7;
    float3  tb_4 = _S1981 + k_7;
    float _S1982 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S1983 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
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
        float t_6 = lerp_0(_S1982, _S1983, (float(i_7) + 0.5f) / 8.0f);
        float3  _S1984 = (ray_o_3 + ray_d_3 * make_float3 (t_6) - pos_10) / make_float3 (size_10);
        FixedArray<float, 8>  _S1985 = densities_9;
        float _S1986 = interp_0(&_S1985, _S1984);
        float _S1987;
        if(_S1986 > 1.10000002384185791f)
        {
            _S1987 = _S1986;
        }
        else
        {
            _S1987 = (F32_exp((0.90909093618392944f * _S1986 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S1987;
        float depth_accum_1 = depth_accum_0 + t_6 * _S1987;
        i_7 = i_7 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_4 = (F32_max((depth_accum_0 / accum_4), (0.0f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S1988;
    float _S1989;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_0, float3  dpray_o_2, float3  dpray_d_2, float3  * dpout_rgb_0, float * dpdepth_0, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_2)
{
    _s_diff_ctx_2->_S1988 = 0.0f;
    _s_diff_ctx_2->_S1989 = 0.0f;
    float _S1990 = 0.5f * size_11;
    float3  m_6 = make_float3 (1.0f) / dpray_d_2;
    float3  k_8 = s_primal_ctx_abs_0(m_6) * make_float3 (_S1990);
    float3  _S1991 = - (m_6 * (dpray_o_2 - (pos_11 + make_float3 (_S1990))));
    float3  ta_5 = _S1991 - k_8;
    float3  tb_5 = _S1991 + k_8;
    float _S1992 = (F32_max(((F32_max((ta_5.x), (ta_5.y)))), ((F32_max((ta_5.z), (0.0f))))));
    float _S1993 = (F32_min(((F32_min((tb_5.x), (tb_5.y)))), (tb_5.z)));
    bool _runFlag_1 = true;
    int i_8 = int(0);
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_2->_S1988 = depth_accum_2;
        _s_diff_ctx_2->_S1989 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        int _S1994;
        float _S1995;
        float _S1996;
        if(i_8 < int(8))
        {
            float _S1997 = s_primal_ctx_lerp_0(_S1992, _S1993, (float(i_8) + 0.5f) / 8.0f);
            float _S1998 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_2 + dpray_d_2 * make_float3 (_S1997) - pos_11) / make_float3 (size_11));
            if(_S1998 > 1.10000002384185791f)
            {
                _S1995 = _S1998;
            }
            else
            {
                _S1995 = s_primal_ctx_exp_0(0.90909093618392944f * _S1998 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S1995;
            float depth_accum_3 = depth_accum_2 + _S1997 * _S1995;
            _S1994 = int(1);
            _S1995 = accum_7;
            _S1996 = depth_accum_3;
        }
        else
        {
            _S1994 = int(0);
            _S1995 = 0.0f;
            _S1996 = 0.0f;
        }
        if(_S1994 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            i_8 = i_8 + int(1);
            accum_6 = _S1995;
            depth_accum_2 = _S1996;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S1999 = (F32_max((depth_accum_2 / accum_6), (0.0f)));
    *dpout_rgb_0 = dprgb_0;
    *dpdepth_0 = _S1999;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_1, float dpdepth_1, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_3)
{
    FixedArray<float, 8>  _S2000 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2001 = *dpray_o_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2002 = *dpray_d_3;
    float3  _S2003 = make_float3 (size_12);
    float _S2004 = 0.5f * size_12;
    float3  _S2005 = make_float3 (_S2004);
    float3  m_7 = make_float3 (1.0f) / (*dpray_d_3).primal_0;
    float3  _S2006 = (*dpray_d_3).primal_0 * (*dpray_d_3).primal_0;
    float3  _S2007 = (*dpray_o_3).primal_0 - (pos_12 + make_float3 (_S2004));
    float3  k_9 = s_primal_ctx_abs_0(m_7) * make_float3 (_S2004);
    float3  _S2008 = - (m_7 * _S2007);
    float3  ta_6 = _S2008 - k_9;
    float3  tb_6 = _S2008 + k_9;
    float _S2009 = ta_6.x;
    float _S2010 = ta_6.y;
    float _S2011 = (F32_max((_S2009), (_S2010)));
    float _S2012 = ta_6.z;
    float _S2013 = (F32_max((_S2012), (0.0f)));
    float _S2014 = (F32_max((_S2011), (_S2013)));
    float _S2015 = tb_6.x;
    float _S2016 = tb_6.y;
    float _S2017 = (F32_min((_S2015), (_S2016)));
    float _S2018 = tb_6.z;
    float _S2019 = (F32_min((_S2017), (_S2018)));
    float _S2020 = _s_diff_ctx_3->_S1989 * _s_diff_ctx_3->_S1989;
    float3  _S2021 = make_float3 (0.0f);
    FixedArray<float, 8>  _S2022 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S2023;
    (&_S2023)->primal_0 = _s_diff_ctx_3->_S1988 / _s_diff_ctx_3->_S1989;
    (&_S2023)->differential_0 = 0.0f;
    DiffPair_float_0 _S2024;
    (&_S2024)->primal_0 = 0.0f;
    (&_S2024)->differential_0 = 0.0f;
    _d_max_0(&_S2023, &_S2024, dpdepth_1);
    float _S2025 = _S2023.differential_0 / _S2020;
    float _S2026 = _s_diff_ctx_3->_S1988 * - _S2025;
    float _S2027 = _s_diff_ctx_3->_S1989 * _S2025;
    int _dc_1 = int(8);
    float _S2028 = _S2026;
    float _S2029 = _S2027;
    FixedArray<float, 8>  _S2030;
    _S2030[int(0)] = 0.0f;
    _S2030[int(1)] = 0.0f;
    _S2030[int(2)] = 0.0f;
    _S2030[int(3)] = 0.0f;
    _S2030[int(4)] = 0.0f;
    _S2030[int(5)] = 0.0f;
    _S2030[int(6)] = 0.0f;
    _S2030[int(7)] = 0.0f;
    float3  _S2031 = _S2021;
    float3  _S2032 = _S2021;
    float _S2033 = 0.0f;
    float _S2034 = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S2035 = _dc_1 < int(8);
        int _S2036;
        float _S2037;
        float _S2038;
        float _S2039;
        float _S2040;
        float3  _S2041;
        float3  _S2042;
        bool _S2043;
        if(_S2035)
        {
            float _S2044 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S2045 = s_primal_ctx_lerp_0(_S2014, _S2019, _S2044);
            float3  _S2046 = make_float3 (_S2045);
            float3  _S2047 = (_S2001.primal_0 + _S2002.primal_0 * make_float3 (_S2045) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S2048 = _S2000;
            float _S2049 = s_primal_ctx_interp_0(&_S2048, _S2047);
            bool _S2050 = _S2049 > 1.10000002384185791f;
            if(_S2050)
            {
                _S2037 = _S2049;
                _S2038 = 0.0f;
            }
            else
            {
                float _S2051 = 0.90909093618392944f * _S2049 - 0.90468984842300415f;
                _S2037 = s_primal_ctx_exp_0(_S2051);
                _S2038 = _S2051;
            }
            float _S2052 = _S2037;
            float _S2053 = _S2038;
            _S2036 = int(1);
            _S2037 = _S2045;
            _S2038 = _S2052;
            _S2043 = _S2050;
            _S2039 = _S2053;
            _S2041 = _S2047;
            _S2042 = _S2046;
            _S2040 = _S2044;
        }
        else
        {
            _S2036 = int(0);
            _S2037 = 0.0f;
            _S2038 = 0.0f;
            _S2043 = false;
            _S2039 = 0.0f;
            _S2041 = _S2021;
            _S2042 = _S2021;
            _S2040 = 0.0f;
        }
        float _S2054;
        float _S2055;
        float _S2056;
        float _S2057;
        if(!(_S2036 != int(1)))
        {
            _S2054 = _S2028;
            _S2055 = _S2029;
            _S2056 = 0.0f;
            _S2057 = 0.0f;
        }
        else
        {
            _S2054 = 0.0f;
            _S2055 = 0.0f;
            _S2056 = _S2029;
            _S2057 = _S2028;
        }
        if(_S2035)
        {
            float _S2058 = _S2038 * _S2055;
            float _S2059 = _S2055 + _S2056;
            float _S2060 = _S2037 * _S2055 + _S2054;
            float _S2061 = _S2054 + _S2057;
            float _S2062;
            if(_S2043)
            {
                _S2062 = _S2060;
            }
            else
            {
                DiffPair_float_0 _S2063;
                (&_S2063)->primal_0 = _S2039;
                (&_S2063)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S2063, _S2060);
                _S2062 = 0.90909093618392944f * _S2063.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S2064;
            (&_S2064)->primal_0 = _S2000;
            (&_S2064)->differential_0 = _S2022;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2065;
            (&_S2065)->primal_0 = _S2041;
            (&_S2065)->differential_0 = _S2021;
            s_bwd_prop_interp_0(&_S2064, &_S2065, _S2062);
            float3  _S2066 = _S2065.differential_0 / _S2003;
            float3  _S2067 = _S2002.primal_0 * _S2066;
            float3  _S2068 = _S2042 * _S2066;
            float _S2069 = _S2067.x + _S2067.y + _S2067.z + _S2058;
            DiffPair_float_0 _S2070;
            (&_S2070)->primal_0 = _S2014;
            (&_S2070)->differential_0 = 0.0f;
            DiffPair_float_0 _S2071;
            (&_S2071)->primal_0 = _S2019;
            (&_S2071)->differential_0 = 0.0f;
            DiffPair_float_0 _S2072;
            (&_S2072)->primal_0 = _S2040;
            (&_S2072)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S2070, &_S2071, &_S2072, _S2069);
            float _S2073 = (&_S2064)->differential_0[int(0)] + _S2030[int(0)];
            float _S2074 = (&_S2064)->differential_0[int(1)] + _S2030[int(1)];
            float _S2075 = (&_S2064)->differential_0[int(2)] + _S2030[int(2)];
            float _S2076 = (&_S2064)->differential_0[int(3)] + _S2030[int(3)];
            float _S2077 = (&_S2064)->differential_0[int(4)] + _S2030[int(4)];
            float _S2078 = (&_S2064)->differential_0[int(5)] + _S2030[int(5)];
            float _S2079 = (&_S2064)->differential_0[int(6)] + _S2030[int(6)];
            float _S2080 = (&_S2064)->differential_0[int(7)] + _S2030[int(7)];
            float3  _S2081 = _S2066 + _S2031;
            float3  _S2082 = _S2068 + _S2032;
            float _S2083 = _S2071.differential_0 + _S2033;
            float _S2084 = _S2070.differential_0 + _S2034;
            _S2028 = _S2061;
            _S2029 = _S2059;
            _S2030[int(0)] = _S2073;
            _S2030[int(1)] = _S2074;
            _S2030[int(2)] = _S2075;
            _S2030[int(3)] = _S2076;
            _S2030[int(4)] = _S2077;
            _S2030[int(5)] = _S2078;
            _S2030[int(6)] = _S2079;
            _S2030[int(7)] = _S2080;
            _S2031 = _S2081;
            _S2032 = _S2082;
            _S2033 = _S2083;
            _S2034 = _S2084;
        }
        else
        {
            _S2028 = _S2057;
            _S2029 = _S2056;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S2085;
    (&_S2085)->primal_0 = _S2017;
    (&_S2085)->differential_0 = 0.0f;
    DiffPair_float_0 _S2086;
    (&_S2086)->primal_0 = _S2018;
    (&_S2086)->differential_0 = 0.0f;
    _d_min_0(&_S2085, &_S2086, _S2033);
    DiffPair_float_0 _S2087;
    (&_S2087)->primal_0 = _S2015;
    (&_S2087)->differential_0 = 0.0f;
    DiffPair_float_0 _S2088;
    (&_S2088)->primal_0 = _S2016;
    (&_S2088)->differential_0 = 0.0f;
    _d_min_0(&_S2087, &_S2088, _S2085.differential_0);
    DiffPair_float_0 _S2089;
    (&_S2089)->primal_0 = _S2011;
    (&_S2089)->differential_0 = 0.0f;
    DiffPair_float_0 _S2090;
    (&_S2090)->primal_0 = _S2013;
    (&_S2090)->differential_0 = 0.0f;
    _d_max_0(&_S2089, &_S2090, _S2034);
    DiffPair_float_0 _S2091;
    (&_S2091)->primal_0 = _S2012;
    (&_S2091)->differential_0 = 0.0f;
    DiffPair_float_0 _S2092;
    (&_S2092)->primal_0 = 0.0f;
    (&_S2092)->differential_0 = 0.0f;
    _d_max_0(&_S2091, &_S2092, _S2090.differential_0);
    DiffPair_float_0 _S2093;
    (&_S2093)->primal_0 = _S2009;
    (&_S2093)->differential_0 = 0.0f;
    DiffPair_float_0 _S2094;
    (&_S2094)->primal_0 = _S2010;
    (&_S2094)->differential_0 = 0.0f;
    _d_max_0(&_S2093, &_S2094, _S2089.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S2087.differential_0, _S2088.differential_0, _S2086.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S2093.differential_0, _S2094.differential_0, _S2091.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S2095 = _S2005 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2096;
    (&_S2096)->primal_0 = m_7;
    (&_S2096)->differential_0 = _S2021;
    s_bwd_prop_abs_0(&_S2096, _S2095);
    float3  _S2097 = m_7 * s_diff_n_T_1;
    float3  _S2098 = - ((_S2096.differential_0 + _S2007 * s_diff_n_T_1) / _S2006) + _S2032;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2098;
    float3  _S2099 = _S2097 + _S2031;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2099;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = dpout_rgb_1;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S2030;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S2100, float _S2101, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S2102, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2103, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2104, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2105, float3  _S2106, float _S2107)
{
    FixedArray<float, 8>  _S2108 = _S2102->primal_0;
    float3  _S2109;
    float _S2110;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2111;
    s_primal_ctx_evaluate_color_voxel_0(_S2100, _S2101, &_S2108, (*_S2103).primal_0, (*_S2104).primal_0, (*_S2105).primal_0, &_S2109, &_S2110, &_S2111);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2112 = _S2111;
    s_bwd_prop_evaluate_color_voxel_0(_S2100, _S2101, _S2102, _S2103, _S2104, _S2105, _S2106, _S2107, &_S2112);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  densities_10, float3  rgb_1, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_0, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    FixedArray<float, 8>  _S2113 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = densities_10;
    (&dp_densities_1)->differential_0 = _S2113;
    float3  _S2114 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S2114;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2114;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2114;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_2 = dp_rgb_0.differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

