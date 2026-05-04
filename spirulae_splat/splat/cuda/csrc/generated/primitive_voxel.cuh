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

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  densities_0, FixedArray<float3 , 16>  sh_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * depth_0, float3  * rgbs_0)
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
        float _S25 = (F32_min((1.00000001504746622e+30f), (_S24.z)));
        float3  _S26 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 0.0f)) + t_0;
        pos_c_0[int(1)] = _S26;
        float _S27 = (F32_min((_S25), (_S26.z)));
        float3  _S28 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 0.0f)) + t_0;
        pos_c_0[int(2)] = _S28;
        float _S29 = (F32_min((_S27), (_S28.z)));
        float3  _S30 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 1.0f, 0.0f)) + t_0;
        pos_c_0[int(3)] = _S30;
        float _S31 = (F32_min((_S29), (_S30.z)));
        float3  _S32 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 0.0f, 1.0f)) + t_0;
        pos_c_0[int(4)] = _S32;
        float _S33 = (F32_min((_S31), (_S32.z)));
        float3  _S34 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 1.0f)) + t_0;
        pos_c_0[int(5)] = _S34;
        float _S35 = (F32_min((_S33), (_S34.z)));
        float3  _S36 = mul_0(R_0, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 1.0f)) + t_0;
        pos_c_0[int(6)] = _S36;
        float _S37 = (F32_min((_S35), (_S36.z)));
        float3  _S38 = mul_0(R_0, pos_0 + make_float3 (size_0)) + t_0;
        pos_c_0[int(7)] = _S38;
        bool _S39 = (F32_min((_S37), (_S38.z))) <= 0.0f;
        if(_S39)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        bool _S40;
        float3  mean_c_0 = mul_0(R_0, pos_0 + make_float3 (0.5f * size_0)) + t_0;
        FixedArray<float2 , 8>  uv_0;
        for(;;)
        {
            float3  _S41 = pos_c_0[int(0)];
            _S15 = &uv_0[int(0)];
            for(;;)
            {
                float _S42 = _S41.z;
                uv_0[int(0)] = float2 {_S41.x, _S41.y} / make_float2 (_S42);
                if(_S42 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_0 = uv_0[int(0)].x;
                    float v_0 = uv_0[int(0)].y;
                    float _S43 = u_0 + u_0;
                    float r2_0 = u_0 * u_0 + v_0 * v_0;
                    float _S44 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                    float _S45 = dist_coeffs_0[int(1)] + r2_0 * _S44;
                    float _S46 = dist_coeffs_0[int(0)] + r2_0 * _S45;
                    float radial_0 = 1.0f + r2_0 * _S46;
                    float _S47 = 2.0f * dist_coeffs_0[int(4)];
                    float _S48 = 2.0f * u_0;
                    float _S49 = 2.0f * dist_coeffs_0[int(5)];
                    float _S50 = 2.0f * v_0;
                    float2  _S51 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (_S43 * _S46 + (_S43 * _S45 + (_S43 * _S44 + _S43 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0[int(0)] + make_float2 (_S47 * v_0 + (_S43 + (_S48 + _S48)) * dist_coeffs_0[int(5)] + _S43 * dist_coeffs_0[int(6)], _S49 * v_0 + _S43 * dist_coeffs_0[int(4)] + _S43 * dist_coeffs_0[int(7)]);
                    float _S52 = v_0 + v_0;
                    float2  _S53 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (_S52 * _S46 + (_S52 * _S45 + (_S52 * _S44 + _S52 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0[int(0)] + make_float2 (_S47 * u_0 + _S52 * dist_coeffs_0[int(5)] + _S52 * dist_coeffs_0[int(6)], _S49 * u_0 + (_S52 + (_S50 + _S50)) * dist_coeffs_0[int(4)] + _S52 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S54 = transpose_0(makeMatrix<float, 2, 2> (_S51 + make_float2 (_S51.x * dist_coeffs_0[int(8)] + _S51.y * dist_coeffs_0[int(9)], 0.0f), _S53 + make_float2 (_S53.x * dist_coeffs_0[int(8)] + _S53.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S54)), ((F32_min((_S54.rows[int(0)].x), (_S54.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_1 = uv_0[int(0)].x;
                float v_1 = uv_0[int(0)].y;
                float r2_1 = u_1 * u_1 + v_1 * v_1;
                float2  _S55 = uv_0[int(0)] * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
                float2  _S56 = _S55 + make_float2 (dist_coeffs_0[int(8)] * _S55.x + dist_coeffs_0[int(9)] * _S55.y, 0.0f);
                uv_0[int(0)] = make_float2 (fx_0 * _S56.x + cx_0, fy_0 * _S56.y + cy_0);
                break;
            }
            bool all_valid_0 = true & (!_S40);
            float3  _S57 = pos_c_0[int(1)];
            _S16 = &uv_0[int(1)];
            for(;;)
            {
                float _S58 = _S57.z;
                uv_0[int(1)] = float2 {_S57.x, _S57.y} / make_float2 (_S58);
                if(_S58 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_2 = uv_0[int(1)].x;
                    float v_2 = uv_0[int(1)].y;
                    float _S59 = u_2 + u_2;
                    float r2_2 = u_2 * u_2 + v_2 * v_2;
                    float _S60 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
                    float _S61 = dist_coeffs_0[int(1)] + r2_2 * _S60;
                    float _S62 = dist_coeffs_0[int(0)] + r2_2 * _S61;
                    float radial_1 = 1.0f + r2_2 * _S62;
                    float _S63 = 2.0f * dist_coeffs_0[int(4)];
                    float _S64 = 2.0f * u_2;
                    float _S65 = 2.0f * dist_coeffs_0[int(5)];
                    float _S66 = 2.0f * v_2;
                    float2  _S67 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (_S59 * _S62 + (_S59 * _S61 + (_S59 * _S60 + _S59 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv_0[int(1)] + make_float2 (_S63 * v_2 + (_S59 + (_S64 + _S64)) * dist_coeffs_0[int(5)] + _S59 * dist_coeffs_0[int(6)], _S65 * v_2 + _S59 * dist_coeffs_0[int(4)] + _S59 * dist_coeffs_0[int(7)]);
                    float _S68 = v_2 + v_2;
                    float2  _S69 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (_S68 * _S62 + (_S68 * _S61 + (_S68 * _S60 + _S68 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * uv_0[int(1)] + make_float2 (_S63 * u_2 + _S68 * dist_coeffs_0[int(5)] + _S68 * dist_coeffs_0[int(6)], _S65 * u_2 + (_S68 + (_S66 + _S66)) * dist_coeffs_0[int(4)] + _S68 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S70 = transpose_0(makeMatrix<float, 2, 2> (_S67 + make_float2 (_S67.x * dist_coeffs_0[int(8)] + _S67.y * dist_coeffs_0[int(9)], 0.0f), _S69 + make_float2 (_S69.x * dist_coeffs_0[int(8)] + _S69.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S70)), ((F32_min((_S70.rows[int(0)].x), (_S70.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_3 = uv_0[int(1)].x;
                float v_3 = uv_0[int(1)].y;
                float r2_3 = u_3 * u_3 + v_3 * v_3;
                float2  _S71 = uv_0[int(1)] * make_float2 (1.0f + r2_3 * (dist_coeffs_0[int(0)] + r2_3 * (dist_coeffs_0[int(1)] + r2_3 * (dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_3 * v_3 + dist_coeffs_0[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + dist_coeffs_0[int(6)] * r2_3, 2.0f * dist_coeffs_0[int(5)] * u_3 * v_3 + dist_coeffs_0[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + dist_coeffs_0[int(7)] * r2_3);
                float2  _S72 = _S71 + make_float2 (dist_coeffs_0[int(8)] * _S71.x + dist_coeffs_0[int(9)] * _S71.y, 0.0f);
                uv_0[int(1)] = make_float2 (fx_0 * _S72.x + cx_0, fy_0 * _S72.y + cy_0);
                break;
            }
            bool all_valid_1 = all_valid_0 & (!_S40);
            float3  _S73 = pos_c_0[int(2)];
            _S17 = &uv_0[int(2)];
            for(;;)
            {
                float _S74 = _S73.z;
                uv_0[int(2)] = float2 {_S73.x, _S73.y} / make_float2 (_S74);
                if(_S74 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_4 = uv_0[int(2)].x;
                    float v_4 = uv_0[int(2)].y;
                    float _S75 = u_4 + u_4;
                    float r2_4 = u_4 * u_4 + v_4 * v_4;
                    float _S76 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
                    float _S77 = dist_coeffs_0[int(1)] + r2_4 * _S76;
                    float _S78 = dist_coeffs_0[int(0)] + r2_4 * _S77;
                    float radial_2 = 1.0f + r2_4 * _S78;
                    float _S79 = 2.0f * dist_coeffs_0[int(4)];
                    float _S80 = 2.0f * u_4;
                    float _S81 = 2.0f * dist_coeffs_0[int(5)];
                    float _S82 = 2.0f * v_4;
                    float2  _S83 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (_S75 * _S78 + (_S75 * _S77 + (_S75 * _S76 + _S75 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv_0[int(2)] + make_float2 (_S79 * v_4 + (_S75 + (_S80 + _S80)) * dist_coeffs_0[int(5)] + _S75 * dist_coeffs_0[int(6)], _S81 * v_4 + _S75 * dist_coeffs_0[int(4)] + _S75 * dist_coeffs_0[int(7)]);
                    float _S84 = v_4 + v_4;
                    float2  _S85 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (_S84 * _S78 + (_S84 * _S77 + (_S84 * _S76 + _S84 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * uv_0[int(2)] + make_float2 (_S79 * u_4 + _S84 * dist_coeffs_0[int(5)] + _S84 * dist_coeffs_0[int(6)], _S81 * u_4 + (_S84 + (_S82 + _S82)) * dist_coeffs_0[int(4)] + _S84 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S86 = transpose_0(makeMatrix<float, 2, 2> (_S83 + make_float2 (_S83.x * dist_coeffs_0[int(8)] + _S83.y * dist_coeffs_0[int(9)], 0.0f), _S85 + make_float2 (_S85.x * dist_coeffs_0[int(8)] + _S85.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S86)), ((F32_min((_S86.rows[int(0)].x), (_S86.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_5 = uv_0[int(2)].x;
                float v_5 = uv_0[int(2)].y;
                float r2_5 = u_5 * u_5 + v_5 * v_5;
                float2  _S87 = uv_0[int(2)] * make_float2 (1.0f + r2_5 * (dist_coeffs_0[int(0)] + r2_5 * (dist_coeffs_0[int(1)] + r2_5 * (dist_coeffs_0[int(2)] + r2_5 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_5 * v_5 + dist_coeffs_0[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + dist_coeffs_0[int(6)] * r2_5, 2.0f * dist_coeffs_0[int(5)] * u_5 * v_5 + dist_coeffs_0[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + dist_coeffs_0[int(7)] * r2_5);
                float2  _S88 = _S87 + make_float2 (dist_coeffs_0[int(8)] * _S87.x + dist_coeffs_0[int(9)] * _S87.y, 0.0f);
                uv_0[int(2)] = make_float2 (fx_0 * _S88.x + cx_0, fy_0 * _S88.y + cy_0);
                break;
            }
            bool all_valid_2 = all_valid_1 & (!_S40);
            float3  _S89 = pos_c_0[int(3)];
            _S18 = &uv_0[int(3)];
            for(;;)
            {
                float _S90 = _S89.z;
                uv_0[int(3)] = float2 {_S89.x, _S89.y} / make_float2 (_S90);
                if(_S90 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_6 = uv_0[int(3)].x;
                    float v_6 = uv_0[int(3)].y;
                    float _S91 = u_6 + u_6;
                    float r2_6 = u_6 * u_6 + v_6 * v_6;
                    float _S92 = dist_coeffs_0[int(2)] + r2_6 * dist_coeffs_0[int(3)];
                    float _S93 = dist_coeffs_0[int(1)] + r2_6 * _S92;
                    float _S94 = dist_coeffs_0[int(0)] + r2_6 * _S93;
                    float radial_3 = 1.0f + r2_6 * _S94;
                    float _S95 = 2.0f * dist_coeffs_0[int(4)];
                    float _S96 = 2.0f * u_6;
                    float _S97 = 2.0f * dist_coeffs_0[int(5)];
                    float _S98 = 2.0f * v_6;
                    float2  _S99 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (_S91 * _S94 + (_S91 * _S93 + (_S91 * _S92 + _S91 * dist_coeffs_0[int(3)] * r2_6) * r2_6) * r2_6) * uv_0[int(3)] + make_float2 (_S95 * v_6 + (_S91 + (_S96 + _S96)) * dist_coeffs_0[int(5)] + _S91 * dist_coeffs_0[int(6)], _S97 * v_6 + _S91 * dist_coeffs_0[int(4)] + _S91 * dist_coeffs_0[int(7)]);
                    float _S100 = v_6 + v_6;
                    float2  _S101 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (_S100 * _S94 + (_S100 * _S93 + (_S100 * _S92 + _S100 * dist_coeffs_0[int(3)] * r2_6) * r2_6) * r2_6) * uv_0[int(3)] + make_float2 (_S95 * u_6 + _S100 * dist_coeffs_0[int(5)] + _S100 * dist_coeffs_0[int(6)], _S97 * u_6 + (_S100 + (_S98 + _S98)) * dist_coeffs_0[int(4)] + _S100 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S102 = transpose_0(makeMatrix<float, 2, 2> (_S99 + make_float2 (_S99.x * dist_coeffs_0[int(8)] + _S99.y * dist_coeffs_0[int(9)], 0.0f), _S101 + make_float2 (_S101.x * dist_coeffs_0[int(8)] + _S101.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S102)), ((F32_min((_S102.rows[int(0)].x), (_S102.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_7 = uv_0[int(3)].x;
                float v_7 = uv_0[int(3)].y;
                float r2_7 = u_7 * u_7 + v_7 * v_7;
                float2  _S103 = uv_0[int(3)] * make_float2 (1.0f + r2_7 * (dist_coeffs_0[int(0)] + r2_7 * (dist_coeffs_0[int(1)] + r2_7 * (dist_coeffs_0[int(2)] + r2_7 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_7 * v_7 + dist_coeffs_0[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + dist_coeffs_0[int(6)] * r2_7, 2.0f * dist_coeffs_0[int(5)] * u_7 * v_7 + dist_coeffs_0[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + dist_coeffs_0[int(7)] * r2_7);
                float2  _S104 = _S103 + make_float2 (dist_coeffs_0[int(8)] * _S103.x + dist_coeffs_0[int(9)] * _S103.y, 0.0f);
                uv_0[int(3)] = make_float2 (fx_0 * _S104.x + cx_0, fy_0 * _S104.y + cy_0);
                break;
            }
            bool all_valid_3 = all_valid_2 & (!_S40);
            float3  _S105 = pos_c_0[int(4)];
            _S19 = &uv_0[int(4)];
            for(;;)
            {
                float _S106 = _S105.z;
                uv_0[int(4)] = float2 {_S105.x, _S105.y} / make_float2 (_S106);
                if(_S106 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_8 = uv_0[int(4)].x;
                    float v_8 = uv_0[int(4)].y;
                    float _S107 = u_8 + u_8;
                    float r2_8 = u_8 * u_8 + v_8 * v_8;
                    float _S108 = dist_coeffs_0[int(2)] + r2_8 * dist_coeffs_0[int(3)];
                    float _S109 = dist_coeffs_0[int(1)] + r2_8 * _S108;
                    float _S110 = dist_coeffs_0[int(0)] + r2_8 * _S109;
                    float radial_4 = 1.0f + r2_8 * _S110;
                    float _S111 = 2.0f * dist_coeffs_0[int(4)];
                    float _S112 = 2.0f * u_8;
                    float _S113 = 2.0f * dist_coeffs_0[int(5)];
                    float _S114 = 2.0f * v_8;
                    float2  _S115 = make_float2 (1.0f, 0.0f) * make_float2 (radial_4) + make_float2 (_S107 * _S110 + (_S107 * _S109 + (_S107 * _S108 + _S107 * dist_coeffs_0[int(3)] * r2_8) * r2_8) * r2_8) * uv_0[int(4)] + make_float2 (_S111 * v_8 + (_S107 + (_S112 + _S112)) * dist_coeffs_0[int(5)] + _S107 * dist_coeffs_0[int(6)], _S113 * v_8 + _S107 * dist_coeffs_0[int(4)] + _S107 * dist_coeffs_0[int(7)]);
                    float _S116 = v_8 + v_8;
                    float2  _S117 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (_S116 * _S110 + (_S116 * _S109 + (_S116 * _S108 + _S116 * dist_coeffs_0[int(3)] * r2_8) * r2_8) * r2_8) * uv_0[int(4)] + make_float2 (_S111 * u_8 + _S116 * dist_coeffs_0[int(5)] + _S116 * dist_coeffs_0[int(6)], _S113 * u_8 + (_S116 + (_S114 + _S114)) * dist_coeffs_0[int(4)] + _S116 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S118 = transpose_0(makeMatrix<float, 2, 2> (_S115 + make_float2 (_S115.x * dist_coeffs_0[int(8)] + _S115.y * dist_coeffs_0[int(9)], 0.0f), _S117 + make_float2 (_S117.x * dist_coeffs_0[int(8)] + _S117.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S118)), ((F32_min((_S118.rows[int(0)].x), (_S118.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_9 = uv_0[int(4)].x;
                float v_9 = uv_0[int(4)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S119 = uv_0[int(4)] * make_float2 (1.0f + r2_9 * (dist_coeffs_0[int(0)] + r2_9 * (dist_coeffs_0[int(1)] + r2_9 * (dist_coeffs_0[int(2)] + r2_9 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_9 * v_9 + dist_coeffs_0[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + dist_coeffs_0[int(6)] * r2_9, 2.0f * dist_coeffs_0[int(5)] * u_9 * v_9 + dist_coeffs_0[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + dist_coeffs_0[int(7)] * r2_9);
                float2  _S120 = _S119 + make_float2 (dist_coeffs_0[int(8)] * _S119.x + dist_coeffs_0[int(9)] * _S119.y, 0.0f);
                uv_0[int(4)] = make_float2 (fx_0 * _S120.x + cx_0, fy_0 * _S120.y + cy_0);
                break;
            }
            bool all_valid_4 = all_valid_3 & (!_S40);
            float3  _S121 = pos_c_0[int(5)];
            _S20 = &uv_0[int(5)];
            for(;;)
            {
                float _S122 = _S121.z;
                uv_0[int(5)] = float2 {_S121.x, _S121.y} / make_float2 (_S122);
                if(_S122 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_10 = uv_0[int(5)].x;
                    float v_10 = uv_0[int(5)].y;
                    float _S123 = u_10 + u_10;
                    float r2_10 = u_10 * u_10 + v_10 * v_10;
                    float _S124 = dist_coeffs_0[int(2)] + r2_10 * dist_coeffs_0[int(3)];
                    float _S125 = dist_coeffs_0[int(1)] + r2_10 * _S124;
                    float _S126 = dist_coeffs_0[int(0)] + r2_10 * _S125;
                    float radial_5 = 1.0f + r2_10 * _S126;
                    float _S127 = 2.0f * dist_coeffs_0[int(4)];
                    float _S128 = 2.0f * u_10;
                    float _S129 = 2.0f * dist_coeffs_0[int(5)];
                    float _S130 = 2.0f * v_10;
                    float2  _S131 = make_float2 (1.0f, 0.0f) * make_float2 (radial_5) + make_float2 (_S123 * _S126 + (_S123 * _S125 + (_S123 * _S124 + _S123 * dist_coeffs_0[int(3)] * r2_10) * r2_10) * r2_10) * uv_0[int(5)] + make_float2 (_S127 * v_10 + (_S123 + (_S128 + _S128)) * dist_coeffs_0[int(5)] + _S123 * dist_coeffs_0[int(6)], _S129 * v_10 + _S123 * dist_coeffs_0[int(4)] + _S123 * dist_coeffs_0[int(7)]);
                    float _S132 = v_10 + v_10;
                    float2  _S133 = make_float2 (0.0f, 1.0f) * make_float2 (radial_5) + make_float2 (_S132 * _S126 + (_S132 * _S125 + (_S132 * _S124 + _S132 * dist_coeffs_0[int(3)] * r2_10) * r2_10) * r2_10) * uv_0[int(5)] + make_float2 (_S127 * u_10 + _S132 * dist_coeffs_0[int(5)] + _S132 * dist_coeffs_0[int(6)], _S129 * u_10 + (_S132 + (_S130 + _S130)) * dist_coeffs_0[int(4)] + _S132 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S134 = transpose_0(makeMatrix<float, 2, 2> (_S131 + make_float2 (_S131.x * dist_coeffs_0[int(8)] + _S131.y * dist_coeffs_0[int(9)], 0.0f), _S133 + make_float2 (_S133.x * dist_coeffs_0[int(8)] + _S133.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S134)), ((F32_min((_S134.rows[int(0)].x), (_S134.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_11 = uv_0[int(5)].x;
                float v_11 = uv_0[int(5)].y;
                float r2_11 = u_11 * u_11 + v_11 * v_11;
                float2  _S135 = uv_0[int(5)] * make_float2 (1.0f + r2_11 * (dist_coeffs_0[int(0)] + r2_11 * (dist_coeffs_0[int(1)] + r2_11 * (dist_coeffs_0[int(2)] + r2_11 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_11 * v_11 + dist_coeffs_0[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_0[int(6)] * r2_11, 2.0f * dist_coeffs_0[int(5)] * u_11 * v_11 + dist_coeffs_0[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_0[int(7)] * r2_11);
                float2  _S136 = _S135 + make_float2 (dist_coeffs_0[int(8)] * _S135.x + dist_coeffs_0[int(9)] * _S135.y, 0.0f);
                uv_0[int(5)] = make_float2 (fx_0 * _S136.x + cx_0, fy_0 * _S136.y + cy_0);
                break;
            }
            bool all_valid_5 = all_valid_4 & (!_S40);
            float3  _S137 = pos_c_0[int(6)];
            _S21 = &uv_0[int(6)];
            for(;;)
            {
                float _S138 = _S137.z;
                uv_0[int(6)] = float2 {_S137.x, _S137.y} / make_float2 (_S138);
                if(_S138 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_12 = uv_0[int(6)].x;
                    float v_12 = uv_0[int(6)].y;
                    float _S139 = u_12 + u_12;
                    float r2_12 = u_12 * u_12 + v_12 * v_12;
                    float _S140 = dist_coeffs_0[int(2)] + r2_12 * dist_coeffs_0[int(3)];
                    float _S141 = dist_coeffs_0[int(1)] + r2_12 * _S140;
                    float _S142 = dist_coeffs_0[int(0)] + r2_12 * _S141;
                    float radial_6 = 1.0f + r2_12 * _S142;
                    float _S143 = 2.0f * dist_coeffs_0[int(4)];
                    float _S144 = 2.0f * u_12;
                    float _S145 = 2.0f * dist_coeffs_0[int(5)];
                    float _S146 = 2.0f * v_12;
                    float2  _S147 = make_float2 (1.0f, 0.0f) * make_float2 (radial_6) + make_float2 (_S139 * _S142 + (_S139 * _S141 + (_S139 * _S140 + _S139 * dist_coeffs_0[int(3)] * r2_12) * r2_12) * r2_12) * uv_0[int(6)] + make_float2 (_S143 * v_12 + (_S139 + (_S144 + _S144)) * dist_coeffs_0[int(5)] + _S139 * dist_coeffs_0[int(6)], _S145 * v_12 + _S139 * dist_coeffs_0[int(4)] + _S139 * dist_coeffs_0[int(7)]);
                    float _S148 = v_12 + v_12;
                    float2  _S149 = make_float2 (0.0f, 1.0f) * make_float2 (radial_6) + make_float2 (_S148 * _S142 + (_S148 * _S141 + (_S148 * _S140 + _S148 * dist_coeffs_0[int(3)] * r2_12) * r2_12) * r2_12) * uv_0[int(6)] + make_float2 (_S143 * u_12 + _S148 * dist_coeffs_0[int(5)] + _S148 * dist_coeffs_0[int(6)], _S145 * u_12 + (_S148 + (_S146 + _S146)) * dist_coeffs_0[int(4)] + _S148 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S150 = transpose_0(makeMatrix<float, 2, 2> (_S147 + make_float2 (_S147.x * dist_coeffs_0[int(8)] + _S147.y * dist_coeffs_0[int(9)], 0.0f), _S149 + make_float2 (_S149.x * dist_coeffs_0[int(8)] + _S149.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S150)), ((F32_min((_S150.rows[int(0)].x), (_S150.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_13 = uv_0[int(6)].x;
                float v_13 = uv_0[int(6)].y;
                float r2_13 = u_13 * u_13 + v_13 * v_13;
                float2  _S151 = uv_0[int(6)] * make_float2 (1.0f + r2_13 * (dist_coeffs_0[int(0)] + r2_13 * (dist_coeffs_0[int(1)] + r2_13 * (dist_coeffs_0[int(2)] + r2_13 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_13 * v_13 + dist_coeffs_0[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + dist_coeffs_0[int(6)] * r2_13, 2.0f * dist_coeffs_0[int(5)] * u_13 * v_13 + dist_coeffs_0[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + dist_coeffs_0[int(7)] * r2_13);
                float2  _S152 = _S151 + make_float2 (dist_coeffs_0[int(8)] * _S151.x + dist_coeffs_0[int(9)] * _S151.y, 0.0f);
                uv_0[int(6)] = make_float2 (fx_0 * _S152.x + cx_0, fy_0 * _S152.y + cy_0);
                break;
            }
            bool all_valid_6 = all_valid_5 & (!_S40);
            float3  _S153 = pos_c_0[int(7)];
            _S22 = &uv_0[int(7)];
            for(;;)
            {
                float _S154 = _S153.z;
                uv_0[int(7)] = float2 {_S153.x, _S153.y} / make_float2 (_S154);
                if(_S154 < 0.0f)
                {
                    _S40 = true;
                }
                else
                {
                    float u_14 = uv_0[int(7)].x;
                    float v_14 = uv_0[int(7)].y;
                    float _S155 = u_14 + u_14;
                    float r2_14 = u_14 * u_14 + v_14 * v_14;
                    float _S156 = dist_coeffs_0[int(2)] + r2_14 * dist_coeffs_0[int(3)];
                    float _S157 = dist_coeffs_0[int(1)] + r2_14 * _S156;
                    float _S158 = dist_coeffs_0[int(0)] + r2_14 * _S157;
                    float radial_7 = 1.0f + r2_14 * _S158;
                    float _S159 = 2.0f * dist_coeffs_0[int(4)];
                    float _S160 = 2.0f * u_14;
                    float _S161 = 2.0f * dist_coeffs_0[int(5)];
                    float _S162 = 2.0f * v_14;
                    float2  _S163 = make_float2 (1.0f, 0.0f) * make_float2 (radial_7) + make_float2 (_S155 * _S158 + (_S155 * _S157 + (_S155 * _S156 + _S155 * dist_coeffs_0[int(3)] * r2_14) * r2_14) * r2_14) * uv_0[int(7)] + make_float2 (_S159 * v_14 + (_S155 + (_S160 + _S160)) * dist_coeffs_0[int(5)] + _S155 * dist_coeffs_0[int(6)], _S161 * v_14 + _S155 * dist_coeffs_0[int(4)] + _S155 * dist_coeffs_0[int(7)]);
                    float _S164 = v_14 + v_14;
                    float2  _S165 = make_float2 (0.0f, 1.0f) * make_float2 (radial_7) + make_float2 (_S164 * _S158 + (_S164 * _S157 + (_S164 * _S156 + _S164 * dist_coeffs_0[int(3)] * r2_14) * r2_14) * r2_14) * uv_0[int(7)] + make_float2 (_S159 * u_14 + _S164 * dist_coeffs_0[int(5)] + _S164 * dist_coeffs_0[int(6)], _S161 * u_14 + (_S164 + (_S162 + _S162)) * dist_coeffs_0[int(4)] + _S164 * dist_coeffs_0[int(7)]);
                    Matrix<float, 2, 2>  _S166 = transpose_0(makeMatrix<float, 2, 2> (_S163 + make_float2 (_S163.x * dist_coeffs_0[int(8)] + _S163.y * dist_coeffs_0[int(9)], 0.0f), _S165 + make_float2 (_S165.x * dist_coeffs_0[int(8)] + _S165.y * dist_coeffs_0[int(9)], 0.0f)));
                    _S40 = !((F32_min((determinant_0(_S166)), ((F32_min((_S166.rows[int(0)].x), (_S166.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S40)
                {
                    break;
                }
                float u_15 = uv_0[int(7)].x;
                float v_15 = uv_0[int(7)].y;
                float r2_15 = u_15 * u_15 + v_15 * v_15;
                float2  _S167 = uv_0[int(7)] * make_float2 (1.0f + r2_15 * (dist_coeffs_0[int(0)] + r2_15 * (dist_coeffs_0[int(1)] + r2_15 * (dist_coeffs_0[int(2)] + r2_15 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_15 * v_15 + dist_coeffs_0[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + dist_coeffs_0[int(6)] * r2_15, 2.0f * dist_coeffs_0[int(5)] * u_15 * v_15 + dist_coeffs_0[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + dist_coeffs_0[int(7)] * r2_15);
                float2  _S168 = _S167 + make_float2 (dist_coeffs_0[int(8)] * _S167.x + dist_coeffs_0[int(9)] * _S167.y, 0.0f);
                uv_0[int(7)] = make_float2 (fx_0 * _S168.x + cx_0, fy_0 * _S168.y + cy_0);
                break;
            }
            _S23 = all_valid_6 & (!_S40);
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
        float _S169 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S15).x), ((*_S16).x)))), ((*_S17).x)))), ((*_S18).x)))), ((*_S19).x)))), ((*_S20).x)))), ((*_S21).x)))), ((*_S22).x)));
        float _S170 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S15).x), ((*_S16).x)))), ((*_S17).x)))), ((*_S18).x)))), ((*_S19).x)))), ((*_S20).x)))), ((*_S21).x)))), ((*_S22).x)));
        float _S171 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S15).y), ((*_S16).y)))), ((*_S17).y)))), ((*_S18).y)))), ((*_S19).y)))), ((*_S20).y)))), ((*_S21).y)))), ((*_S22).y)));
        float _S172 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S15).y), ((*_S16).y)))), ((*_S17).y)))), ((*_S18).y)))), ((*_S19).y)))), ((*_S20).y)))), ((*_S21).y)))), ((*_S22).y)));
        if(_S169 <= 0.0f)
        {
            _S40 = true;
        }
        else
        {
            _S40 = _S170 >= float(image_width_0);
        }
        if(_S40)
        {
            _S40 = true;
        }
        else
        {
            _S40 = _S171 <= 0.0f;
        }
        if(_S40)
        {
            _S40 = true;
        }
        else
        {
            _S40 = _S172 >= float(image_height_0);
        }
        if(_S40)
        {
            _S40 = true;
        }
        else
        {
            if(_S39)
            {
                if(_S170 <= 0.0f)
                {
                    _S40 = _S169 >= float(image_width_0);
                }
                else
                {
                    _S40 = false;
                }
                if(_S40)
                {
                    _S40 = true;
                }
                else
                {
                    if(_S172 <= 0.0f)
                    {
                        _S40 = _S171 >= float(image_width_0);
                    }
                    else
                    {
                        _S40 = false;
                    }
                }
            }
            else
            {
                _S40 = false;
            }
        }
        if(_S40)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (float(int((F32_floor((_S170))))), float(int((F32_floor((_S172))))), float(int((F32_ceil((_S169))))), float(int((F32_ceil((_S171))))));
        *depth_0 = mean_c_0.z;
        float3  _S173 = mean_c_0 - - mul_0(transpose_1(R_0), t_0);
        float _S174 = _S173.x;
        float _S175 = _S173.y;
        float _S176 = _S173.z;
        float norm_0 = (F32_sqrt((_S174 * _S174 + _S175 * _S175 + _S176 * _S176)));
        float x_7 = _S174 / norm_0;
        float y_3 = _S175 / norm_0;
        float z_0 = _S176 / norm_0;
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

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  densities_1, FixedArray<float3 , 16>  sh_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * depth_1, float3  * rgbs_1)
{
    float2  * _S177;
    float2  _S178;
    float _S179;
    float _S180;
    float _S181;
    float _S182;
    float _S183;
    float _S184;
    float _S185;
    float _S186;
    float _S187;
    float _S188;
    float _S189;
    float _S190;
    float2  _S191;
    bool _S192;
    float2  * _S193;
    bool _S194;
    float2  * _S195;
    bool _S196;
    float2  * _S197;
    bool _S198;
    float2  * _S199;
    bool _S200;
    float2  * _S201;
    bool _S202;
    float2  * _S203;
    bool _S204;
    float2  * _S205;
    bool _S206;
    bool _S207;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S208 = mul_0(R_1, pos_1) + t_1;
        pos_c_1[int(0)] = _S208;
        float _S209 = (F32_min((1.00000001504746622e+30f), (length_1(_S208))));
        float3  _S210 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_1;
        pos_c_1[int(1)] = _S210;
        float _S211 = (F32_min((_S209), (length_1(_S210))));
        float3  _S212 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(2)] = _S212;
        float _S213 = (F32_min((_S211), (length_1(_S212))));
        float3  _S214 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_1;
        pos_c_1[int(3)] = _S214;
        float _S215 = (F32_min((_S213), (length_1(_S214))));
        float3  _S216 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(4)] = _S216;
        float _S217 = (F32_min((_S215), (length_1(_S216))));
        float3  _S218 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_1;
        pos_c_1[int(5)] = _S218;
        float _S219 = (F32_min((_S217), (length_1(_S218))));
        float3  _S220 = mul_0(R_1, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_1;
        pos_c_1[int(6)] = _S220;
        float _S221 = (F32_min((_S219), (length_1(_S220))));
        float3  _S222 = mul_0(R_1, pos_1 + make_float3 (size_1)) + t_1;
        pos_c_1[int(7)] = _S222;
        bool _S223 = (F32_min((_S221), (length_1(_S222)))) <= 0.0f;
        if(_S223)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float3  mean_c_1 = mul_0(R_1, pos_1 + make_float3 (0.5f * size_1)) + t_1;
        FixedArray<float2 , 8>  uv_1;
        for(;;)
        {
            float k_0;
            float3  _S224 = pos_c_1[int(0)];
            _S177 = &uv_1[int(0)];
            for(;;)
            {
                float2  _S225 = float2 {_S224.x, _S224.y};
                float r_2 = length_0(_S225);
                float _S226 = _S224.z;
                float theta_0 = (F32_atan2((r_2), (_S226)));
                if(theta_0 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S226;
                }
                else
                {
                    k_0 = theta_0 / r_2;
                }
                float2  _S227 = _S225 * make_float2 (k_0);
                uv_1[int(0)] = _S227;
                float2  _S228 = make_float2 (1.0f, 0.0f);
                _S178 = _S228;
                _S179 = dist_coeffs_1[int(0)];
                _S180 = dist_coeffs_1[int(1)];
                _S181 = dist_coeffs_1[int(2)];
                _S182 = dist_coeffs_1[int(3)];
                _S183 = dist_coeffs_1[int(4)];
                _S184 = dist_coeffs_1[int(5)];
                _S185 = dist_coeffs_1[int(6)];
                _S186 = dist_coeffs_1[int(7)];
                _S187 = dist_coeffs_1[int(8)];
                _S188 = dist_coeffs_1[int(9)];
                float u_16 = _S227.x;
                float v_16 = _S227.y;
                float _S229 = u_16 + u_16;
                float r2_16 = u_16 * u_16 + v_16 * v_16;
                float _S230 = dist_coeffs_1[int(2)] + r2_16 * dist_coeffs_1[int(3)];
                float _S231 = dist_coeffs_1[int(1)] + r2_16 * _S230;
                float _S232 = dist_coeffs_1[int(0)] + r2_16 * _S231;
                float _S233 = _S229 * _S232 + (_S229 * _S231 + (_S229 * _S230 + _S229 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16;
                float radial_8 = 1.0f + r2_16 * _S232;
                float _S234 = 2.0f * dist_coeffs_1[int(4)];
                _S189 = _S234;
                float _S235 = _S234 * u_16;
                float _S236 = 2.0f * u_16;
                float s_diff_du_0 = _S234 * v_16 + (_S229 + (_S236 + _S236)) * dist_coeffs_1[int(5)] + _S229 * dist_coeffs_1[int(6)];
                float _S237 = 2.0f * dist_coeffs_1[int(5)];
                _S190 = _S237;
                float _S238 = _S237 * u_16;
                float _S239 = 2.0f * v_16;
                float2  _S240 = _S228 * make_float2 (radial_8) + make_float2 (_S233) * _S227 + make_float2 (s_diff_du_0, _S237 * v_16 + _S229 * dist_coeffs_1[int(4)] + _S229 * dist_coeffs_1[int(7)]);
                float2  _S241 = _S240 + make_float2 (_S240.x * dist_coeffs_1[int(8)] + _S240.y * dist_coeffs_1[int(9)], 0.0f);
                float2  _S242 = make_float2 (0.0f, 1.0f);
                _S191 = _S242;
                float _S243 = v_16 + v_16;
                float2  _S244 = _S242 * make_float2 (radial_8) + make_float2 (_S243 * _S232 + (_S243 * _S231 + (_S243 * _S230 + _S243 * dist_coeffs_1[int(3)] * r2_16) * r2_16) * r2_16) * _S227 + make_float2 (_S235 + _S243 * dist_coeffs_1[int(5)] + _S243 * dist_coeffs_1[int(6)], _S238 + (_S243 + (_S239 + _S239)) * dist_coeffs_1[int(4)] + _S243 * dist_coeffs_1[int(7)]);
                Matrix<float, 2, 2>  _S245 = transpose_0(makeMatrix<float, 2, 2> (_S241, _S244 + make_float2 (_S244.x * dist_coeffs_1[int(8)] + _S244.y * dist_coeffs_1[int(9)], 0.0f)));
                bool _S246 = !((F32_min((determinant_0(_S245)), ((F32_min((_S245.rows[int(0)].x), (_S245.rows[int(1)].y)))))) > 0.0f);
                _S192 = _S246;
                if(_S246)
                {
                    break;
                }
                float u_17 = uv_1[int(0)].x;
                float v_17 = uv_1[int(0)].y;
                float r2_17 = u_17 * u_17 + v_17 * v_17;
                float2  _S247 = uv_1[int(0)] * make_float2 (1.0f + r2_17 * (dist_coeffs_1[int(0)] + r2_17 * (dist_coeffs_1[int(1)] + r2_17 * (dist_coeffs_1[int(2)] + r2_17 * dist_coeffs_1[int(3)])))) + make_float2 (_S234 * u_17 * v_17 + dist_coeffs_1[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + dist_coeffs_1[int(6)] * r2_17, _S237 * u_17 * v_17 + dist_coeffs_1[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + dist_coeffs_1[int(7)] * r2_17);
                float2  _S248 = _S247 + make_float2 (dist_coeffs_1[int(8)] * _S247.x + dist_coeffs_1[int(9)] * _S247.y, 0.0f);
                uv_1[int(0)] = make_float2 (fx_1 * _S248.x + cx_1, fy_1 * _S248.y + cy_1);
                break;
            }
            bool all_valid_7 = true & (!_S192);
            float3  _S249 = pos_c_1[int(1)];
            _S193 = &uv_1[int(1)];
            for(;;)
            {
                float2  _S250 = float2 {_S249.x, _S249.y};
                float r_3 = length_0(_S250);
                float _S251 = _S249.z;
                float theta_1 = (F32_atan2((r_3), (_S251)));
                if(theta_1 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S251;
                }
                else
                {
                    k_0 = theta_1 / r_3;
                }
                float2  _S252 = _S250 * make_float2 (k_0);
                uv_1[int(1)] = _S252;
                float u_18 = _S252.x;
                float v_18 = _S252.y;
                float _S253 = u_18 + u_18;
                float r2_18 = u_18 * u_18 + v_18 * v_18;
                float _S254 = _S181 + r2_18 * _S182;
                float _S255 = _S180 + r2_18 * _S254;
                float _S256 = _S179 + r2_18 * _S255;
                float radial_9 = 1.0f + r2_18 * _S256;
                float _S257 = 2.0f * u_18;
                float _S258 = 2.0f * v_18;
                float2  _S259 = _S178 * make_float2 (radial_9) + make_float2 (_S253 * _S256 + (_S253 * _S255 + (_S253 * _S254 + _S253 * _S182 * r2_18) * r2_18) * r2_18) * _S252 + make_float2 (_S189 * v_18 + (_S253 + (_S257 + _S257)) * _S184 + _S253 * _S185, _S190 * v_18 + _S253 * _S183 + _S253 * _S186);
                float _S260 = v_18 + v_18;
                float2  _S261 = _S191 * make_float2 (radial_9) + make_float2 (_S260 * _S256 + (_S260 * _S255 + (_S260 * _S254 + _S260 * _S182 * r2_18) * r2_18) * r2_18) * _S252 + make_float2 (_S189 * u_18 + _S260 * _S184 + _S260 * _S185, _S190 * u_18 + (_S260 + (_S258 + _S258)) * _S183 + _S260 * _S186);
                Matrix<float, 2, 2>  _S262 = transpose_0(makeMatrix<float, 2, 2> (_S259 + make_float2 (_S259.x * _S187 + _S259.y * _S188, 0.0f), _S261 + make_float2 (_S261.x * _S187 + _S261.y * _S188, 0.0f)));
                bool _S263 = !((F32_min((determinant_0(_S262)), ((F32_min((_S262.rows[int(0)].x), (_S262.rows[int(1)].y)))))) > 0.0f);
                _S194 = _S263;
                if(_S263)
                {
                    break;
                }
                float u_19 = uv_1[int(1)].x;
                float v_19 = uv_1[int(1)].y;
                float r2_19 = u_19 * u_19 + v_19 * v_19;
                float2  _S264 = uv_1[int(1)] * make_float2 (1.0f + r2_19 * (_S179 + r2_19 * (_S180 + r2_19 * (_S181 + r2_19 * _S182)))) + make_float2 (_S189 * u_19 * v_19 + _S184 * (r2_19 + 2.0f * u_19 * u_19) + _S185 * r2_19, _S190 * u_19 * v_19 + _S183 * (r2_19 + 2.0f * v_19 * v_19) + _S186 * r2_19);
                float2  _S265 = _S264 + make_float2 (_S187 * _S264.x + _S188 * _S264.y, 0.0f);
                uv_1[int(1)] = make_float2 (fx_1 * _S265.x + cx_1, fy_1 * _S265.y + cy_1);
                break;
            }
            bool all_valid_8 = all_valid_7 & (!_S194);
            float3  _S266 = pos_c_1[int(2)];
            _S195 = &uv_1[int(2)];
            for(;;)
            {
                float2  _S267 = float2 {_S266.x, _S266.y};
                float r_4 = length_0(_S267);
                float _S268 = _S266.z;
                float theta_2 = (F32_atan2((r_4), (_S268)));
                if(theta_2 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_2 * theta_2 / 3.0f) / _S268;
                }
                else
                {
                    k_0 = theta_2 / r_4;
                }
                float2  _S269 = _S267 * make_float2 (k_0);
                uv_1[int(2)] = _S269;
                float u_20 = _S269.x;
                float v_20 = _S269.y;
                float _S270 = u_20 + u_20;
                float r2_20 = u_20 * u_20 + v_20 * v_20;
                float _S271 = _S181 + r2_20 * _S182;
                float _S272 = _S180 + r2_20 * _S271;
                float _S273 = _S179 + r2_20 * _S272;
                float radial_10 = 1.0f + r2_20 * _S273;
                float _S274 = 2.0f * u_20;
                float _S275 = 2.0f * v_20;
                float2  _S276 = _S178 * make_float2 (radial_10) + make_float2 (_S270 * _S273 + (_S270 * _S272 + (_S270 * _S271 + _S270 * _S182 * r2_20) * r2_20) * r2_20) * _S269 + make_float2 (_S189 * v_20 + (_S270 + (_S274 + _S274)) * _S184 + _S270 * _S185, _S190 * v_20 + _S270 * _S183 + _S270 * _S186);
                float _S277 = v_20 + v_20;
                float2  _S278 = _S191 * make_float2 (radial_10) + make_float2 (_S277 * _S273 + (_S277 * _S272 + (_S277 * _S271 + _S277 * _S182 * r2_20) * r2_20) * r2_20) * _S269 + make_float2 (_S189 * u_20 + _S277 * _S184 + _S277 * _S185, _S190 * u_20 + (_S277 + (_S275 + _S275)) * _S183 + _S277 * _S186);
                Matrix<float, 2, 2>  _S279 = transpose_0(makeMatrix<float, 2, 2> (_S276 + make_float2 (_S276.x * _S187 + _S276.y * _S188, 0.0f), _S278 + make_float2 (_S278.x * _S187 + _S278.y * _S188, 0.0f)));
                bool _S280 = !((F32_min((determinant_0(_S279)), ((F32_min((_S279.rows[int(0)].x), (_S279.rows[int(1)].y)))))) > 0.0f);
                _S196 = _S280;
                if(_S280)
                {
                    break;
                }
                float u_21 = uv_1[int(2)].x;
                float v_21 = uv_1[int(2)].y;
                float r2_21 = u_21 * u_21 + v_21 * v_21;
                float2  _S281 = uv_1[int(2)] * make_float2 (1.0f + r2_21 * (_S179 + r2_21 * (_S180 + r2_21 * (_S181 + r2_21 * _S182)))) + make_float2 (_S189 * u_21 * v_21 + _S184 * (r2_21 + 2.0f * u_21 * u_21) + _S185 * r2_21, _S190 * u_21 * v_21 + _S183 * (r2_21 + 2.0f * v_21 * v_21) + _S186 * r2_21);
                float2  _S282 = _S281 + make_float2 (_S187 * _S281.x + _S188 * _S281.y, 0.0f);
                uv_1[int(2)] = make_float2 (fx_1 * _S282.x + cx_1, fy_1 * _S282.y + cy_1);
                break;
            }
            bool all_valid_9 = all_valid_8 & (!_S196);
            float3  _S283 = pos_c_1[int(3)];
            _S197 = &uv_1[int(3)];
            for(;;)
            {
                float2  _S284 = float2 {_S283.x, _S283.y};
                float r_5 = length_0(_S284);
                float _S285 = _S283.z;
                float theta_3 = (F32_atan2((r_5), (_S285)));
                if(theta_3 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S285;
                }
                else
                {
                    k_0 = theta_3 / r_5;
                }
                float2  _S286 = _S284 * make_float2 (k_0);
                uv_1[int(3)] = _S286;
                float u_22 = _S286.x;
                float v_22 = _S286.y;
                float _S287 = u_22 + u_22;
                float r2_22 = u_22 * u_22 + v_22 * v_22;
                float _S288 = _S181 + r2_22 * _S182;
                float _S289 = _S180 + r2_22 * _S288;
                float _S290 = _S179 + r2_22 * _S289;
                float radial_11 = 1.0f + r2_22 * _S290;
                float _S291 = 2.0f * u_22;
                float _S292 = 2.0f * v_22;
                float2  _S293 = _S178 * make_float2 (radial_11) + make_float2 (_S287 * _S290 + (_S287 * _S289 + (_S287 * _S288 + _S287 * _S182 * r2_22) * r2_22) * r2_22) * _S286 + make_float2 (_S189 * v_22 + (_S287 + (_S291 + _S291)) * _S184 + _S287 * _S185, _S190 * v_22 + _S287 * _S183 + _S287 * _S186);
                float _S294 = v_22 + v_22;
                float2  _S295 = _S191 * make_float2 (radial_11) + make_float2 (_S294 * _S290 + (_S294 * _S289 + (_S294 * _S288 + _S294 * _S182 * r2_22) * r2_22) * r2_22) * _S286 + make_float2 (_S189 * u_22 + _S294 * _S184 + _S294 * _S185, _S190 * u_22 + (_S294 + (_S292 + _S292)) * _S183 + _S294 * _S186);
                Matrix<float, 2, 2>  _S296 = transpose_0(makeMatrix<float, 2, 2> (_S293 + make_float2 (_S293.x * _S187 + _S293.y * _S188, 0.0f), _S295 + make_float2 (_S295.x * _S187 + _S295.y * _S188, 0.0f)));
                bool _S297 = !((F32_min((determinant_0(_S296)), ((F32_min((_S296.rows[int(0)].x), (_S296.rows[int(1)].y)))))) > 0.0f);
                _S198 = _S297;
                if(_S297)
                {
                    break;
                }
                float u_23 = uv_1[int(3)].x;
                float v_23 = uv_1[int(3)].y;
                float r2_23 = u_23 * u_23 + v_23 * v_23;
                float2  _S298 = uv_1[int(3)] * make_float2 (1.0f + r2_23 * (_S179 + r2_23 * (_S180 + r2_23 * (_S181 + r2_23 * _S182)))) + make_float2 (_S189 * u_23 * v_23 + _S184 * (r2_23 + 2.0f * u_23 * u_23) + _S185 * r2_23, _S190 * u_23 * v_23 + _S183 * (r2_23 + 2.0f * v_23 * v_23) + _S186 * r2_23);
                float2  _S299 = _S298 + make_float2 (_S187 * _S298.x + _S188 * _S298.y, 0.0f);
                uv_1[int(3)] = make_float2 (fx_1 * _S299.x + cx_1, fy_1 * _S299.y + cy_1);
                break;
            }
            bool all_valid_10 = all_valid_9 & (!_S198);
            float3  _S300 = pos_c_1[int(4)];
            _S199 = &uv_1[int(4)];
            for(;;)
            {
                float2  _S301 = float2 {_S300.x, _S300.y};
                float r_6 = length_0(_S301);
                float _S302 = _S300.z;
                float theta_4 = (F32_atan2((r_6), (_S302)));
                if(theta_4 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S302;
                }
                else
                {
                    k_0 = theta_4 / r_6;
                }
                float2  _S303 = _S301 * make_float2 (k_0);
                uv_1[int(4)] = _S303;
                float u_24 = _S303.x;
                float v_24 = _S303.y;
                float _S304 = u_24 + u_24;
                float r2_24 = u_24 * u_24 + v_24 * v_24;
                float _S305 = _S181 + r2_24 * _S182;
                float _S306 = _S180 + r2_24 * _S305;
                float _S307 = _S179 + r2_24 * _S306;
                float radial_12 = 1.0f + r2_24 * _S307;
                float _S308 = 2.0f * u_24;
                float _S309 = 2.0f * v_24;
                float2  _S310 = _S178 * make_float2 (radial_12) + make_float2 (_S304 * _S307 + (_S304 * _S306 + (_S304 * _S305 + _S304 * _S182 * r2_24) * r2_24) * r2_24) * _S303 + make_float2 (_S189 * v_24 + (_S304 + (_S308 + _S308)) * _S184 + _S304 * _S185, _S190 * v_24 + _S304 * _S183 + _S304 * _S186);
                float _S311 = v_24 + v_24;
                float2  _S312 = _S191 * make_float2 (radial_12) + make_float2 (_S311 * _S307 + (_S311 * _S306 + (_S311 * _S305 + _S311 * _S182 * r2_24) * r2_24) * r2_24) * _S303 + make_float2 (_S189 * u_24 + _S311 * _S184 + _S311 * _S185, _S190 * u_24 + (_S311 + (_S309 + _S309)) * _S183 + _S311 * _S186);
                Matrix<float, 2, 2>  _S313 = transpose_0(makeMatrix<float, 2, 2> (_S310 + make_float2 (_S310.x * _S187 + _S310.y * _S188, 0.0f), _S312 + make_float2 (_S312.x * _S187 + _S312.y * _S188, 0.0f)));
                bool _S314 = !((F32_min((determinant_0(_S313)), ((F32_min((_S313.rows[int(0)].x), (_S313.rows[int(1)].y)))))) > 0.0f);
                _S200 = _S314;
                if(_S314)
                {
                    break;
                }
                float u_25 = uv_1[int(4)].x;
                float v_25 = uv_1[int(4)].y;
                float r2_25 = u_25 * u_25 + v_25 * v_25;
                float2  _S315 = uv_1[int(4)] * make_float2 (1.0f + r2_25 * (_S179 + r2_25 * (_S180 + r2_25 * (_S181 + r2_25 * _S182)))) + make_float2 (_S189 * u_25 * v_25 + _S184 * (r2_25 + 2.0f * u_25 * u_25) + _S185 * r2_25, _S190 * u_25 * v_25 + _S183 * (r2_25 + 2.0f * v_25 * v_25) + _S186 * r2_25);
                float2  _S316 = _S315 + make_float2 (_S187 * _S315.x + _S188 * _S315.y, 0.0f);
                uv_1[int(4)] = make_float2 (fx_1 * _S316.x + cx_1, fy_1 * _S316.y + cy_1);
                break;
            }
            bool all_valid_11 = all_valid_10 & (!_S200);
            float3  _S317 = pos_c_1[int(5)];
            _S201 = &uv_1[int(5)];
            for(;;)
            {
                float2  _S318 = float2 {_S317.x, _S317.y};
                float r_7 = length_0(_S318);
                float _S319 = _S317.z;
                float theta_5 = (F32_atan2((r_7), (_S319)));
                if(theta_5 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_5 * theta_5 / 3.0f) / _S319;
                }
                else
                {
                    k_0 = theta_5 / r_7;
                }
                float2  _S320 = _S318 * make_float2 (k_0);
                uv_1[int(5)] = _S320;
                float u_26 = _S320.x;
                float v_26 = _S320.y;
                float _S321 = u_26 + u_26;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float _S322 = _S181 + r2_26 * _S182;
                float _S323 = _S180 + r2_26 * _S322;
                float _S324 = _S179 + r2_26 * _S323;
                float radial_13 = 1.0f + r2_26 * _S324;
                float _S325 = 2.0f * u_26;
                float _S326 = 2.0f * v_26;
                float2  _S327 = _S178 * make_float2 (radial_13) + make_float2 (_S321 * _S324 + (_S321 * _S323 + (_S321 * _S322 + _S321 * _S182 * r2_26) * r2_26) * r2_26) * _S320 + make_float2 (_S189 * v_26 + (_S321 + (_S325 + _S325)) * _S184 + _S321 * _S185, _S190 * v_26 + _S321 * _S183 + _S321 * _S186);
                float _S328 = v_26 + v_26;
                float2  _S329 = _S191 * make_float2 (radial_13) + make_float2 (_S328 * _S324 + (_S328 * _S323 + (_S328 * _S322 + _S328 * _S182 * r2_26) * r2_26) * r2_26) * _S320 + make_float2 (_S189 * u_26 + _S328 * _S184 + _S328 * _S185, _S190 * u_26 + (_S328 + (_S326 + _S326)) * _S183 + _S328 * _S186);
                Matrix<float, 2, 2>  _S330 = transpose_0(makeMatrix<float, 2, 2> (_S327 + make_float2 (_S327.x * _S187 + _S327.y * _S188, 0.0f), _S329 + make_float2 (_S329.x * _S187 + _S329.y * _S188, 0.0f)));
                bool _S331 = !((F32_min((determinant_0(_S330)), ((F32_min((_S330.rows[int(0)].x), (_S330.rows[int(1)].y)))))) > 0.0f);
                _S202 = _S331;
                if(_S331)
                {
                    break;
                }
                float u_27 = uv_1[int(5)].x;
                float v_27 = uv_1[int(5)].y;
                float r2_27 = u_27 * u_27 + v_27 * v_27;
                float2  _S332 = uv_1[int(5)] * make_float2 (1.0f + r2_27 * (_S179 + r2_27 * (_S180 + r2_27 * (_S181 + r2_27 * _S182)))) + make_float2 (_S189 * u_27 * v_27 + _S184 * (r2_27 + 2.0f * u_27 * u_27) + _S185 * r2_27, _S190 * u_27 * v_27 + _S183 * (r2_27 + 2.0f * v_27 * v_27) + _S186 * r2_27);
                float2  _S333 = _S332 + make_float2 (_S187 * _S332.x + _S188 * _S332.y, 0.0f);
                uv_1[int(5)] = make_float2 (fx_1 * _S333.x + cx_1, fy_1 * _S333.y + cy_1);
                break;
            }
            bool all_valid_12 = all_valid_11 & (!_S202);
            float3  _S334 = pos_c_1[int(6)];
            _S203 = &uv_1[int(6)];
            for(;;)
            {
                float2  _S335 = float2 {_S334.x, _S334.y};
                float r_8 = length_0(_S335);
                float _S336 = _S334.z;
                float theta_6 = (F32_atan2((r_8), (_S336)));
                if(theta_6 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_6 * theta_6 / 3.0f) / _S336;
                }
                else
                {
                    k_0 = theta_6 / r_8;
                }
                float2  _S337 = _S335 * make_float2 (k_0);
                uv_1[int(6)] = _S337;
                float u_28 = _S337.x;
                float v_28 = _S337.y;
                float _S338 = u_28 + u_28;
                float r2_28 = u_28 * u_28 + v_28 * v_28;
                float _S339 = _S181 + r2_28 * _S182;
                float _S340 = _S180 + r2_28 * _S339;
                float _S341 = _S179 + r2_28 * _S340;
                float radial_14 = 1.0f + r2_28 * _S341;
                float _S342 = 2.0f * u_28;
                float _S343 = 2.0f * v_28;
                float2  _S344 = _S178 * make_float2 (radial_14) + make_float2 (_S338 * _S341 + (_S338 * _S340 + (_S338 * _S339 + _S338 * _S182 * r2_28) * r2_28) * r2_28) * _S337 + make_float2 (_S189 * v_28 + (_S338 + (_S342 + _S342)) * _S184 + _S338 * _S185, _S190 * v_28 + _S338 * _S183 + _S338 * _S186);
                float _S345 = v_28 + v_28;
                float2  _S346 = _S191 * make_float2 (radial_14) + make_float2 (_S345 * _S341 + (_S345 * _S340 + (_S345 * _S339 + _S345 * _S182 * r2_28) * r2_28) * r2_28) * _S337 + make_float2 (_S189 * u_28 + _S345 * _S184 + _S345 * _S185, _S190 * u_28 + (_S345 + (_S343 + _S343)) * _S183 + _S345 * _S186);
                Matrix<float, 2, 2>  _S347 = transpose_0(makeMatrix<float, 2, 2> (_S344 + make_float2 (_S344.x * _S187 + _S344.y * _S188, 0.0f), _S346 + make_float2 (_S346.x * _S187 + _S346.y * _S188, 0.0f)));
                bool _S348 = !((F32_min((determinant_0(_S347)), ((F32_min((_S347.rows[int(0)].x), (_S347.rows[int(1)].y)))))) > 0.0f);
                _S204 = _S348;
                if(_S348)
                {
                    break;
                }
                float u_29 = uv_1[int(6)].x;
                float v_29 = uv_1[int(6)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S349 = uv_1[int(6)] * make_float2 (1.0f + r2_29 * (_S179 + r2_29 * (_S180 + r2_29 * (_S181 + r2_29 * _S182)))) + make_float2 (_S189 * u_29 * v_29 + _S184 * (r2_29 + 2.0f * u_29 * u_29) + _S185 * r2_29, _S190 * u_29 * v_29 + _S183 * (r2_29 + 2.0f * v_29 * v_29) + _S186 * r2_29);
                float2  _S350 = _S349 + make_float2 (_S187 * _S349.x + _S188 * _S349.y, 0.0f);
                uv_1[int(6)] = make_float2 (fx_1 * _S350.x + cx_1, fy_1 * _S350.y + cy_1);
                break;
            }
            bool all_valid_13 = all_valid_12 & (!_S204);
            float3  _S351 = pos_c_1[int(7)];
            _S205 = &uv_1[int(7)];
            for(;;)
            {
                float2  _S352 = float2 {_S351.x, _S351.y};
                float r_9 = length_0(_S352);
                float _S353 = _S351.z;
                float theta_7 = (F32_atan2((r_9), (_S353)));
                if(theta_7 < 0.00100000004749745f)
                {
                    k_0 = (1.0f - theta_7 * theta_7 / 3.0f) / _S353;
                }
                else
                {
                    k_0 = theta_7 / r_9;
                }
                float2  _S354 = _S352 * make_float2 (k_0);
                uv_1[int(7)] = _S354;
                float u_30 = _S354.x;
                float v_30 = _S354.y;
                float _S355 = u_30 + u_30;
                float r2_30 = u_30 * u_30 + v_30 * v_30;
                float _S356 = _S181 + r2_30 * _S182;
                float _S357 = _S180 + r2_30 * _S356;
                float _S358 = _S179 + r2_30 * _S357;
                float radial_15 = 1.0f + r2_30 * _S358;
                float _S359 = 2.0f * u_30;
                float _S360 = 2.0f * v_30;
                float2  _S361 = _S178 * make_float2 (radial_15) + make_float2 (_S355 * _S358 + (_S355 * _S357 + (_S355 * _S356 + _S355 * _S182 * r2_30) * r2_30) * r2_30) * _S354 + make_float2 (_S189 * v_30 + (_S355 + (_S359 + _S359)) * _S184 + _S355 * _S185, _S190 * v_30 + _S355 * _S183 + _S355 * _S186);
                float _S362 = v_30 + v_30;
                float2  _S363 = _S191 * make_float2 (radial_15) + make_float2 (_S362 * _S358 + (_S362 * _S357 + (_S362 * _S356 + _S362 * _S182 * r2_30) * r2_30) * r2_30) * _S354 + make_float2 (_S189 * u_30 + _S362 * _S184 + _S362 * _S185, _S190 * u_30 + (_S362 + (_S360 + _S360)) * _S183 + _S362 * _S186);
                Matrix<float, 2, 2>  _S364 = transpose_0(makeMatrix<float, 2, 2> (_S361 + make_float2 (_S361.x * _S187 + _S361.y * _S188, 0.0f), _S363 + make_float2 (_S363.x * _S187 + _S363.y * _S188, 0.0f)));
                bool _S365 = !((F32_min((determinant_0(_S364)), ((F32_min((_S364.rows[int(0)].x), (_S364.rows[int(1)].y)))))) > 0.0f);
                _S206 = _S365;
                if(_S365)
                {
                    break;
                }
                float u_31 = uv_1[int(7)].x;
                float v_31 = uv_1[int(7)].y;
                float r2_31 = u_31 * u_31 + v_31 * v_31;
                float2  _S366 = uv_1[int(7)] * make_float2 (1.0f + r2_31 * (_S179 + r2_31 * (_S180 + r2_31 * (_S181 + r2_31 * _S182)))) + make_float2 (_S189 * u_31 * v_31 + _S184 * (r2_31 + 2.0f * u_31 * u_31) + _S185 * r2_31, _S190 * u_31 * v_31 + _S183 * (r2_31 + 2.0f * v_31 * v_31) + _S186 * r2_31);
                float2  _S367 = _S366 + make_float2 (_S187 * _S366.x + _S188 * _S366.y, 0.0f);
                uv_1[int(7)] = make_float2 (fx_1 * _S367.x + cx_1, fy_1 * _S367.y + cy_1);
                break;
            }
            _S207 = all_valid_13 & (!_S206);
            break;
        }
        if(!_S207)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((densities_1[int(0)]), (densities_1[int(1)])))), (densities_1[int(2)])))), (densities_1[int(3)])))), (densities_1[int(4)])))), (densities_1[int(5)])))), (densities_1[int(6)])))), (densities_1[int(7)]))) * size_1 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S368 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S177).x), ((*_S193).x)))), ((*_S195).x)))), ((*_S197).x)))), ((*_S199).x)))), ((*_S201).x)))), ((*_S203).x)))), ((*_S205).x)));
        float _S369 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S177).x), ((*_S193).x)))), ((*_S195).x)))), ((*_S197).x)))), ((*_S199).x)))), ((*_S201).x)))), ((*_S203).x)))), ((*_S205).x)));
        float _S370 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S177).y), ((*_S193).y)))), ((*_S195).y)))), ((*_S197).y)))), ((*_S199).y)))), ((*_S201).y)))), ((*_S203).y)))), ((*_S205).y)));
        float _S371 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S177).y), ((*_S193).y)))), ((*_S195).y)))), ((*_S197).y)))), ((*_S199).y)))), ((*_S201).y)))), ((*_S203).y)))), ((*_S205).y)));
        bool _S372;
        if(_S368 <= 0.0f)
        {
            _S372 = true;
        }
        else
        {
            _S372 = _S369 >= float(image_width_1);
        }
        if(_S372)
        {
            _S372 = true;
        }
        else
        {
            _S372 = _S370 <= 0.0f;
        }
        if(_S372)
        {
            _S372 = true;
        }
        else
        {
            _S372 = _S371 >= float(image_height_1);
        }
        if(_S372)
        {
            _S372 = true;
        }
        else
        {
            if(_S223)
            {
                if(_S369 <= 0.0f)
                {
                    _S372 = _S368 >= float(image_width_1);
                }
                else
                {
                    _S372 = false;
                }
                if(_S372)
                {
                    _S372 = true;
                }
                else
                {
                    if(_S371 <= 0.0f)
                    {
                        _S372 = _S370 >= float(image_width_1);
                    }
                    else
                    {
                        _S372 = false;
                    }
                }
            }
            else
            {
                _S372 = false;
            }
        }
        if(_S372)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (float(int((F32_floor((_S369))))), float(int((F32_floor((_S371))))), float(int((F32_ceil((_S368))))), float(int((F32_ceil((_S370))))));
        float x_8 = mean_c_1.x;
        float y_4 = mean_c_1.y;
        float z_1 = mean_c_1.z;
        float _S373 = x_8 * x_8 + y_4 * y_4;
        *depth_1 = z_1 * z_1 * z_1 * z_1 + 0.001953125f * _S373 * _S373;
        float3  _S374 = mean_c_1 - - mul_0(transpose_1(R_1), t_1);
        float _S375 = _S374.x;
        float _S376 = _S374.y;
        float _S377 = _S374.z;
        float norm_1 = (F32_sqrt((_S375 * _S375 + _S376 * _S376 + _S377 * _S377)));
        float x_9 = _S375 / norm_1;
        float y_5 = _S376 / norm_1;
        float z_2 = _S377 / norm_1;
        float z2_1 = z_2 * z_2;
        float fTmp0B_1 = -1.09254848957061768f * z_2;
        float fC1_1 = x_9 * x_9 - y_5 * y_5;
        float fS1_1 = 2.0f * x_9 * y_5;
        float fTmp0C_1 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_2;
        *rgbs_1 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_1[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * sh_coeffs_1[int(1)] + make_float3 (z_2) * sh_coeffs_1[int(2)] - make_float3 (x_9) * sh_coeffs_1[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * sh_coeffs_1[int(4)] + make_float3 (fTmp0B_1 * y_5) * sh_coeffs_1[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * sh_coeffs_1[int(6)] + make_float3 (fTmp0B_1 * x_9) * sh_coeffs_1[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * sh_coeffs_1[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_1 + y_5 * fC1_1)) * sh_coeffs_1[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * sh_coeffs_1[int(10)] + make_float3 (fTmp0C_1 * y_5) * sh_coeffs_1[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * sh_coeffs_1[int(12)] + make_float3 (fTmp0C_1 * x_9) * sh_coeffs_1[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * sh_coeffs_1[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_1 - y_5 * fS1_1)) * sh_coeffs_1[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  densities_2, FixedArray<float3 , 16>  sh_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * depth_2, float3  * rgbs_2)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S378 = mul_0(R_2, pos_2) + t_2;
    pos_c_2[int(0)] = _S378;
    float3  _S379 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_2;
    pos_c_2[int(1)] = _S379;
    float3  _S380 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(2)] = _S380;
    float3  _S381 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_2;
    pos_c_2[int(3)] = _S381;
    float3  _S382 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(4)] = _S382;
    float3  _S383 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_2;
    pos_c_2[int(5)] = _S383;
    float3  _S384 = mul_0(R_2, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_2;
    pos_c_2[int(6)] = _S384;
    float3  _S385 = mul_0(R_2, pos_2 + make_float3 (size_2)) + t_2;
    pos_c_2[int(7)] = _S385;
    float3  mean_c_2 = mul_0(R_2, pos_2 + make_float3 (0.5f * size_2)) + t_2;
    FixedArray<float2 , 8>  uv_2;
    float2  _S386 = float2 {_S378.x, _S378.y} / make_float2 (_S378.z);
    float u_32 = _S386.x;
    float v_32 = _S386.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float _S387 = 2.0f * dist_coeffs_2[int(4)];
    float _S388 = 2.0f * dist_coeffs_2[int(5)];
    float2  _S389 = _S386 * make_float2 (1.0f + r2_32 * (dist_coeffs_2[int(0)] + r2_32 * (dist_coeffs_2[int(1)] + r2_32 * (dist_coeffs_2[int(2)] + r2_32 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_32 * v_32 + dist_coeffs_2[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + dist_coeffs_2[int(6)] * r2_32, _S388 * u_32 * v_32 + dist_coeffs_2[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + dist_coeffs_2[int(7)] * r2_32);
    float2  _S390 = _S389 + make_float2 (dist_coeffs_2[int(8)] * _S389.x + dist_coeffs_2[int(9)] * _S389.y, 0.0f);
    float _S391 = fx_2 * _S390.x + cx_2;
    float _S392 = fy_2 * _S390.y + cy_2;
    uv_2[int(0)] = make_float2 (_S391, _S392);
    float2  _S393 = float2 {_S379.x, _S379.y} / make_float2 (_S379.z);
    float u_33 = _S393.x;
    float v_33 = _S393.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S394 = _S393 * make_float2 (1.0f + r2_33 * (dist_coeffs_2[int(0)] + r2_33 * (dist_coeffs_2[int(1)] + r2_33 * (dist_coeffs_2[int(2)] + r2_33 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_33 * v_33 + dist_coeffs_2[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + dist_coeffs_2[int(6)] * r2_33, _S388 * u_33 * v_33 + dist_coeffs_2[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + dist_coeffs_2[int(7)] * r2_33);
    float2  _S395 = _S394 + make_float2 (dist_coeffs_2[int(8)] * _S394.x + dist_coeffs_2[int(9)] * _S394.y, 0.0f);
    float _S396 = fx_2 * _S395.x + cx_2;
    float _S397 = fy_2 * _S395.y + cy_2;
    uv_2[int(1)] = make_float2 (_S396, _S397);
    float2  _S398 = float2 {_S380.x, _S380.y} / make_float2 (_S380.z);
    float u_34 = _S398.x;
    float v_34 = _S398.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S399 = _S398 * make_float2 (1.0f + r2_34 * (dist_coeffs_2[int(0)] + r2_34 * (dist_coeffs_2[int(1)] + r2_34 * (dist_coeffs_2[int(2)] + r2_34 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_34 * v_34 + dist_coeffs_2[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + dist_coeffs_2[int(6)] * r2_34, _S388 * u_34 * v_34 + dist_coeffs_2[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + dist_coeffs_2[int(7)] * r2_34);
    float2  _S400 = _S399 + make_float2 (dist_coeffs_2[int(8)] * _S399.x + dist_coeffs_2[int(9)] * _S399.y, 0.0f);
    float _S401 = fx_2 * _S400.x + cx_2;
    float _S402 = fy_2 * _S400.y + cy_2;
    uv_2[int(2)] = make_float2 (_S401, _S402);
    float2  _S403 = float2 {_S381.x, _S381.y} / make_float2 (_S381.z);
    float u_35 = _S403.x;
    float v_35 = _S403.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S404 = _S403 * make_float2 (1.0f + r2_35 * (dist_coeffs_2[int(0)] + r2_35 * (dist_coeffs_2[int(1)] + r2_35 * (dist_coeffs_2[int(2)] + r2_35 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_35 * v_35 + dist_coeffs_2[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + dist_coeffs_2[int(6)] * r2_35, _S388 * u_35 * v_35 + dist_coeffs_2[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + dist_coeffs_2[int(7)] * r2_35);
    float2  _S405 = _S404 + make_float2 (dist_coeffs_2[int(8)] * _S404.x + dist_coeffs_2[int(9)] * _S404.y, 0.0f);
    float _S406 = fx_2 * _S405.x + cx_2;
    float _S407 = fy_2 * _S405.y + cy_2;
    uv_2[int(3)] = make_float2 (_S406, _S407);
    float2  _S408 = float2 {_S382.x, _S382.y} / make_float2 (_S382.z);
    float u_36 = _S408.x;
    float v_36 = _S408.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S409 = _S408 * make_float2 (1.0f + r2_36 * (dist_coeffs_2[int(0)] + r2_36 * (dist_coeffs_2[int(1)] + r2_36 * (dist_coeffs_2[int(2)] + r2_36 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_36 * v_36 + dist_coeffs_2[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + dist_coeffs_2[int(6)] * r2_36, _S388 * u_36 * v_36 + dist_coeffs_2[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + dist_coeffs_2[int(7)] * r2_36);
    float2  _S410 = _S409 + make_float2 (dist_coeffs_2[int(8)] * _S409.x + dist_coeffs_2[int(9)] * _S409.y, 0.0f);
    float _S411 = fx_2 * _S410.x + cx_2;
    float _S412 = fy_2 * _S410.y + cy_2;
    uv_2[int(4)] = make_float2 (_S411, _S412);
    float2  _S413 = float2 {_S383.x, _S383.y} / make_float2 (_S383.z);
    float u_37 = _S413.x;
    float v_37 = _S413.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S414 = _S413 * make_float2 (1.0f + r2_37 * (dist_coeffs_2[int(0)] + r2_37 * (dist_coeffs_2[int(1)] + r2_37 * (dist_coeffs_2[int(2)] + r2_37 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_37 * v_37 + dist_coeffs_2[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + dist_coeffs_2[int(6)] * r2_37, _S388 * u_37 * v_37 + dist_coeffs_2[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + dist_coeffs_2[int(7)] * r2_37);
    float2  _S415 = _S414 + make_float2 (dist_coeffs_2[int(8)] * _S414.x + dist_coeffs_2[int(9)] * _S414.y, 0.0f);
    float _S416 = fx_2 * _S415.x + cx_2;
    float _S417 = fy_2 * _S415.y + cy_2;
    uv_2[int(5)] = make_float2 (_S416, _S417);
    float2  _S418 = float2 {_S384.x, _S384.y} / make_float2 (_S384.z);
    float u_38 = _S418.x;
    float v_38 = _S418.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float2  _S419 = _S418 * make_float2 (1.0f + r2_38 * (dist_coeffs_2[int(0)] + r2_38 * (dist_coeffs_2[int(1)] + r2_38 * (dist_coeffs_2[int(2)] + r2_38 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_38 * v_38 + dist_coeffs_2[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + dist_coeffs_2[int(6)] * r2_38, _S388 * u_38 * v_38 + dist_coeffs_2[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + dist_coeffs_2[int(7)] * r2_38);
    float2  _S420 = _S419 + make_float2 (dist_coeffs_2[int(8)] * _S419.x + dist_coeffs_2[int(9)] * _S419.y, 0.0f);
    float _S421 = fx_2 * _S420.x + cx_2;
    float _S422 = fy_2 * _S420.y + cy_2;
    uv_2[int(6)] = make_float2 (_S421, _S422);
    float2  _S423 = float2 {_S385.x, _S385.y} / make_float2 (_S385.z);
    float u_39 = _S423.x;
    float v_39 = _S423.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float2  _S424 = _S423 * make_float2 (1.0f + r2_39 * (dist_coeffs_2[int(0)] + r2_39 * (dist_coeffs_2[int(1)] + r2_39 * (dist_coeffs_2[int(2)] + r2_39 * dist_coeffs_2[int(3)])))) + make_float2 (_S387 * u_39 * v_39 + dist_coeffs_2[int(5)] * (r2_39 + 2.0f * u_39 * u_39) + dist_coeffs_2[int(6)] * r2_39, _S388 * u_39 * v_39 + dist_coeffs_2[int(4)] * (r2_39 + 2.0f * v_39 * v_39) + dist_coeffs_2[int(7)] * r2_39);
    float2  _S425 = _S424 + make_float2 (dist_coeffs_2[int(8)] * _S424.x + dist_coeffs_2[int(9)] * _S424.y, 0.0f);
    float _S426 = fx_2 * _S425.x + cx_2;
    float _S427 = fy_2 * _S425.y + cy_2;
    uv_2[int(7)] = make_float2 (_S426, _S427);
    *aabb_xyxy_2 = make_float4 (float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S391), (_S396)))), (_S401)))), (_S406)))), (_S411)))), (_S416)))), (_S421)))), (_S426)))))))), float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S392), (_S397)))), (_S402)))), (_S407)))), (_S412)))), (_S417)))), (_S422)))), (_S427)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S391), (_S396)))), (_S401)))), (_S406)))), (_S411)))), (_S416)))), (_S421)))), (_S426)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S392), (_S397)))), (_S402)))), (_S407)))), (_S412)))), (_S417)))), (_S422)))), (_S427)))))))));
    *depth_2 = mean_c_2.z;
    float3  _S428 = mean_c_2 - - mul_0(transpose_1(R_2), t_2);
    float _S429 = _S428.x;
    float _S430 = _S428.y;
    float _S431 = _S428.z;
    float norm_2 = (F32_sqrt((_S429 * _S429 + _S430 * _S430 + _S431 * _S431)));
    float x_10 = _S429 / norm_2;
    float y_6 = _S430 / norm_2;
    float z_3 = _S431 / norm_2;
    float z2_2 = z_3 * z_3;
    float fTmp0B_2 = -1.09254848957061768f * z_3;
    float fC1_2 = x_10 * x_10 - y_6 * y_6;
    float fS1_2 = 2.0f * x_10 * y_6;
    float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_3;
    *rgbs_2 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_2[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * sh_coeffs_2[int(1)] + make_float3 (z_3) * sh_coeffs_2[int(2)] - make_float3 (x_10) * sh_coeffs_2[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * sh_coeffs_2[int(4)] + make_float3 (fTmp0B_2 * y_6) * sh_coeffs_2[int(5)] + make_float3 (0.94617468118667603f * z2_2 - 0.31539157032966614f) * sh_coeffs_2[int(6)] + make_float3 (fTmp0B_2 * x_10) * sh_coeffs_2[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * sh_coeffs_2[int(8)]) + (make_float3 (-0.59004360437393188f * (x_10 * fS1_2 + y_6 * fC1_2)) * sh_coeffs_2[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * sh_coeffs_2[int(10)] + make_float3 (fTmp0C_2 * y_6) * sh_coeffs_2[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_2 - 1.11952900886535645f)) * sh_coeffs_2[int(12)] + make_float3 (fTmp0C_2 * x_10) * sh_coeffs_2[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * sh_coeffs_2[int(14)] + make_float3 (-0.59004360437393188f * (x_10 * fC1_2 - y_6 * fS1_2)) * sh_coeffs_2[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  densities_3, FixedArray<float3 , 16>  sh_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * depth_3, float3  * rgbs_3)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S432 = mul_0(R_3, pos_3) + t_3;
    pos_c_3[int(0)] = _S432;
    pos_c_3[int(1)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_3;
    pos_c_3[int(2)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(3)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_3;
    pos_c_3[int(4)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(5)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_3;
    pos_c_3[int(6)] = mul_0(R_3, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_3;
    pos_c_3[int(7)] = mul_0(R_3, pos_3 + make_float3 (size_3)) + t_3;
    float3  mean_c_3 = mul_0(R_3, pos_3 + make_float3 (0.5f * size_3)) + t_3;
    FixedArray<float2 , 8>  uv_3;
    float2  _S433 = float2 {_S432.x, _S432.y};
    float r_10 = length_0(_S433);
    float _S434 = _S432.z;
    float theta_8 = (F32_atan2((r_10), (_S434)));
    float k_1;
    if(theta_8 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_8 * theta_8 / 3.0f) / _S434;
    }
    else
    {
        k_1 = theta_8 / r_10;
    }
    float2  _S435 = _S433 * make_float2 (k_1);
    float u_40 = _S435.x;
    float v_40 = _S435.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S436 = 2.0f * dist_coeffs_3[int(4)];
    float _S437 = 2.0f * dist_coeffs_3[int(5)];
    float2  _S438 = _S435 * make_float2 (1.0f + r2_40 * (dist_coeffs_3[int(0)] + r2_40 * (dist_coeffs_3[int(1)] + r2_40 * (dist_coeffs_3[int(2)] + r2_40 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_40 * v_40 + dist_coeffs_3[int(5)] * (r2_40 + 2.0f * u_40 * u_40) + dist_coeffs_3[int(6)] * r2_40, _S437 * u_40 * v_40 + dist_coeffs_3[int(4)] * (r2_40 + 2.0f * v_40 * v_40) + dist_coeffs_3[int(7)] * r2_40);
    float2  _S439 = _S438 + make_float2 (dist_coeffs_3[int(8)] * _S438.x + dist_coeffs_3[int(9)] * _S438.y, 0.0f);
    uv_3[int(0)] = make_float2 (fx_3 * _S439.x + cx_3, fy_3 * _S439.y + cy_3);
    float2  _S440 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_11 = length_0(_S440);
    float _S441 = pos_c_3[int(1)].z;
    float theta_9 = (F32_atan2((r_11), (_S441)));
    if(theta_9 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_9 * theta_9 / 3.0f) / _S441;
    }
    else
    {
        k_1 = theta_9 / r_11;
    }
    float2  _S442 = _S440 * make_float2 (k_1);
    float u_41 = _S442.x;
    float v_41 = _S442.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float2  _S443 = _S442 * make_float2 (1.0f + r2_41 * (dist_coeffs_3[int(0)] + r2_41 * (dist_coeffs_3[int(1)] + r2_41 * (dist_coeffs_3[int(2)] + r2_41 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_41 * v_41 + dist_coeffs_3[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + dist_coeffs_3[int(6)] * r2_41, _S437 * u_41 * v_41 + dist_coeffs_3[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + dist_coeffs_3[int(7)] * r2_41);
    float2  _S444 = _S443 + make_float2 (dist_coeffs_3[int(8)] * _S443.x + dist_coeffs_3[int(9)] * _S443.y, 0.0f);
    uv_3[int(1)] = make_float2 (fx_3 * _S444.x + cx_3, fy_3 * _S444.y + cy_3);
    float2  _S445 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_12 = length_0(_S445);
    float _S446 = pos_c_3[int(2)].z;
    float theta_10 = (F32_atan2((r_12), (_S446)));
    if(theta_10 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_10 * theta_10 / 3.0f) / _S446;
    }
    else
    {
        k_1 = theta_10 / r_12;
    }
    float2  _S447 = _S445 * make_float2 (k_1);
    float u_42 = _S447.x;
    float v_42 = _S447.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float2  _S448 = _S447 * make_float2 (1.0f + r2_42 * (dist_coeffs_3[int(0)] + r2_42 * (dist_coeffs_3[int(1)] + r2_42 * (dist_coeffs_3[int(2)] + r2_42 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_42 * v_42 + dist_coeffs_3[int(5)] * (r2_42 + 2.0f * u_42 * u_42) + dist_coeffs_3[int(6)] * r2_42, _S437 * u_42 * v_42 + dist_coeffs_3[int(4)] * (r2_42 + 2.0f * v_42 * v_42) + dist_coeffs_3[int(7)] * r2_42);
    float2  _S449 = _S448 + make_float2 (dist_coeffs_3[int(8)] * _S448.x + dist_coeffs_3[int(9)] * _S448.y, 0.0f);
    uv_3[int(2)] = make_float2 (fx_3 * _S449.x + cx_3, fy_3 * _S449.y + cy_3);
    float2  _S450 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_13 = length_0(_S450);
    float _S451 = pos_c_3[int(3)].z;
    float theta_11 = (F32_atan2((r_13), (_S451)));
    if(theta_11 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_11 * theta_11 / 3.0f) / _S451;
    }
    else
    {
        k_1 = theta_11 / r_13;
    }
    float2  _S452 = _S450 * make_float2 (k_1);
    float u_43 = _S452.x;
    float v_43 = _S452.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float2  _S453 = _S452 * make_float2 (1.0f + r2_43 * (dist_coeffs_3[int(0)] + r2_43 * (dist_coeffs_3[int(1)] + r2_43 * (dist_coeffs_3[int(2)] + r2_43 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_43 * v_43 + dist_coeffs_3[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + dist_coeffs_3[int(6)] * r2_43, _S437 * u_43 * v_43 + dist_coeffs_3[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + dist_coeffs_3[int(7)] * r2_43);
    float2  _S454 = _S453 + make_float2 (dist_coeffs_3[int(8)] * _S453.x + dist_coeffs_3[int(9)] * _S453.y, 0.0f);
    uv_3[int(3)] = make_float2 (fx_3 * _S454.x + cx_3, fy_3 * _S454.y + cy_3);
    float2  _S455 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_14 = length_0(_S455);
    float _S456 = pos_c_3[int(4)].z;
    float theta_12 = (F32_atan2((r_14), (_S456)));
    if(theta_12 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_12 * theta_12 / 3.0f) / _S456;
    }
    else
    {
        k_1 = theta_12 / r_14;
    }
    float2  _S457 = _S455 * make_float2 (k_1);
    float u_44 = _S457.x;
    float v_44 = _S457.y;
    float r2_44 = u_44 * u_44 + v_44 * v_44;
    float2  _S458 = _S457 * make_float2 (1.0f + r2_44 * (dist_coeffs_3[int(0)] + r2_44 * (dist_coeffs_3[int(1)] + r2_44 * (dist_coeffs_3[int(2)] + r2_44 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_44 * v_44 + dist_coeffs_3[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + dist_coeffs_3[int(6)] * r2_44, _S437 * u_44 * v_44 + dist_coeffs_3[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + dist_coeffs_3[int(7)] * r2_44);
    float2  _S459 = _S458 + make_float2 (dist_coeffs_3[int(8)] * _S458.x + dist_coeffs_3[int(9)] * _S458.y, 0.0f);
    uv_3[int(4)] = make_float2 (fx_3 * _S459.x + cx_3, fy_3 * _S459.y + cy_3);
    float2  _S460 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_15 = length_0(_S460);
    float _S461 = pos_c_3[int(5)].z;
    float theta_13 = (F32_atan2((r_15), (_S461)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_13 * theta_13 / 3.0f) / _S461;
    }
    else
    {
        k_1 = theta_13 / r_15;
    }
    float2  _S462 = _S460 * make_float2 (k_1);
    float u_45 = _S462.x;
    float v_45 = _S462.y;
    float r2_45 = u_45 * u_45 + v_45 * v_45;
    float2  _S463 = _S462 * make_float2 (1.0f + r2_45 * (dist_coeffs_3[int(0)] + r2_45 * (dist_coeffs_3[int(1)] + r2_45 * (dist_coeffs_3[int(2)] + r2_45 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_45 * v_45 + dist_coeffs_3[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + dist_coeffs_3[int(6)] * r2_45, _S437 * u_45 * v_45 + dist_coeffs_3[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + dist_coeffs_3[int(7)] * r2_45);
    float2  _S464 = _S463 + make_float2 (dist_coeffs_3[int(8)] * _S463.x + dist_coeffs_3[int(9)] * _S463.y, 0.0f);
    uv_3[int(5)] = make_float2 (fx_3 * _S464.x + cx_3, fy_3 * _S464.y + cy_3);
    float2  _S465 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_16 = length_0(_S465);
    float _S466 = pos_c_3[int(6)].z;
    float theta_14 = (F32_atan2((r_16), (_S466)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_14 * theta_14 / 3.0f) / _S466;
    }
    else
    {
        k_1 = theta_14 / r_16;
    }
    float2  _S467 = _S465 * make_float2 (k_1);
    float u_46 = _S467.x;
    float v_46 = _S467.y;
    float r2_46 = u_46 * u_46 + v_46 * v_46;
    float2  _S468 = _S467 * make_float2 (1.0f + r2_46 * (dist_coeffs_3[int(0)] + r2_46 * (dist_coeffs_3[int(1)] + r2_46 * (dist_coeffs_3[int(2)] + r2_46 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_46 * v_46 + dist_coeffs_3[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + dist_coeffs_3[int(6)] * r2_46, _S437 * u_46 * v_46 + dist_coeffs_3[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + dist_coeffs_3[int(7)] * r2_46);
    float2  _S469 = _S468 + make_float2 (dist_coeffs_3[int(8)] * _S468.x + dist_coeffs_3[int(9)] * _S468.y, 0.0f);
    uv_3[int(6)] = make_float2 (fx_3 * _S469.x + cx_3, fy_3 * _S469.y + cy_3);
    float2  _S470 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_17 = length_0(_S470);
    float _S471 = pos_c_3[int(7)].z;
    float theta_15 = (F32_atan2((r_17), (_S471)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_15 * theta_15 / 3.0f) / _S471;
    }
    else
    {
        k_1 = theta_15 / r_17;
    }
    float2  _S472 = _S470 * make_float2 (k_1);
    float u_47 = _S472.x;
    float v_47 = _S472.y;
    float r2_47 = u_47 * u_47 + v_47 * v_47;
    float2  _S473 = _S472 * make_float2 (1.0f + r2_47 * (dist_coeffs_3[int(0)] + r2_47 * (dist_coeffs_3[int(1)] + r2_47 * (dist_coeffs_3[int(2)] + r2_47 * dist_coeffs_3[int(3)])))) + make_float2 (_S436 * u_47 * v_47 + dist_coeffs_3[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + dist_coeffs_3[int(6)] * r2_47, _S437 * u_47 * v_47 + dist_coeffs_3[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + dist_coeffs_3[int(7)] * r2_47);
    float2  _S474 = _S473 + make_float2 (dist_coeffs_3[int(8)] * _S473.x + dist_coeffs_3[int(9)] * _S473.y, 0.0f);
    float _S475 = fx_3 * _S474.x + cx_3;
    float _S476 = fy_3 * _S474.y + cy_3;
    uv_3[int(7)] = make_float2 (_S475, _S476);
    *aabb_xyxy_3 = make_float4 (float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S475)))))))), float(int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S476)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S475)))))))), float(int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S476)))))))));
    float x_11 = mean_c_3.x;
    float y_7 = mean_c_3.y;
    float z_4 = mean_c_3.z;
    float _S477 = x_11 * x_11 + y_7 * y_7;
    *depth_3 = z_4 * z_4 * z_4 * z_4 + 0.001953125f * _S477 * _S477;
    float3  _S478 = mean_c_3 - - mul_0(transpose_1(R_3), t_3);
    float _S479 = _S478.x;
    float _S480 = _S478.y;
    float _S481 = _S478.z;
    float norm_3 = (F32_sqrt((_S479 * _S479 + _S480 * _S480 + _S481 * _S481)));
    float x_12 = _S479 / norm_3;
    float y_8 = _S480 / norm_3;
    float z_5 = _S481 / norm_3;
    float z2_3 = z_5 * z_5;
    float fTmp0B_3 = -1.09254848957061768f * z_5;
    float fC1_3 = x_12 * x_12 - y_8 * y_8;
    float fS1_3 = 2.0f * x_12 * y_8;
    float fTmp0C_3 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_5;
    *rgbs_3 = max_0(make_float3 (0.282094806432724f) * sh_coeffs_3[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * sh_coeffs_3[int(1)] + make_float3 (z_5) * sh_coeffs_3[int(2)] - make_float3 (x_12) * sh_coeffs_3[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * sh_coeffs_3[int(4)] + make_float3 (fTmp0B_3 * y_8) * sh_coeffs_3[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * sh_coeffs_3[int(6)] + make_float3 (fTmp0B_3 * x_12) * sh_coeffs_3[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * sh_coeffs_3[int(8)]) + (make_float3 (-0.59004360437393188f * (x_12 * fS1_3 + y_8 * fC1_3)) * sh_coeffs_3[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * sh_coeffs_3[int(10)] + make_float3 (fTmp0C_3 * y_8) * sh_coeffs_3[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * sh_coeffs_3[int(12)] + make_float3 (fTmp0C_3 * x_12) * sh_coeffs_3[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * sh_coeffs_3[int(14)] + make_float3 (-0.59004360437393188f * (x_12 * fC1_3 - y_8 * fS1_3)) * sh_coeffs_3[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S482, float3  _S483)
{
    return mul_0(_S482, _S483);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S484)
{
    return (F32_sqrt((_S484)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S485, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S486, float3  _S487)
{
    _d_max_vector_0(_S485, _S486, _S487);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S488, float _S489)
{
    _d_sqrt_0(_S488, _S489);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S490, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S491, float3  _S492)
{
    _d_mul_0(_S490, _S491, _S492);
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  densities_4, FixedArray<float3 , 16>  sh_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float3  v_rgb_0, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  _S493 = s_primal_ctx_mul_0(R_4, pos_4) + t_4;
    float _S494 = _S493.z;
    float2  _S495 = make_float2 (_S494);
    float _S496 = (F32_min((1.00000001504746622e+30f), (_S494)));
    float _S497 = (F32_max((0.0f), (_S494)));
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S498 = s_primal_ctx_mul_0(R_4, pos_i_0) + t_4;
    float _S499 = _S498.z;
    float2  _S500 = make_float2 (_S499);
    float _S501 = (F32_min((_S496), (_S499)));
    float _S502 = (F32_max((_S497), (_S499)));
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S503 = s_primal_ctx_mul_0(R_4, pos_i_1) + t_4;
    float _S504 = _S503.z;
    float2  _S505 = make_float2 (_S504);
    float _S506 = (F32_min((_S501), (_S504)));
    float _S507 = (F32_max((_S502), (_S504)));
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S508 = s_primal_ctx_mul_0(R_4, pos_i_2) + t_4;
    float _S509 = _S508.z;
    float2  _S510 = make_float2 (_S509);
    float _S511 = (F32_min((_S506), (_S509)));
    float _S512 = (F32_max((_S507), (_S509)));
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S513 = s_primal_ctx_mul_0(R_4, pos_i_3) + t_4;
    float _S514 = _S513.z;
    float2  _S515 = make_float2 (_S514);
    float _S516 = (F32_min((_S511), (_S514)));
    float _S517 = (F32_max((_S512), (_S514)));
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S518 = s_primal_ctx_mul_0(R_4, pos_i_4) + t_4;
    float _S519 = _S518.z;
    float2  _S520 = make_float2 (_S519);
    float _S521 = (F32_min((_S516), (_S519)));
    float _S522 = (F32_max((_S517), (_S519)));
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S523 = s_primal_ctx_mul_0(R_4, pos_i_5) + t_4;
    float _S524 = _S523.z;
    float2  _S525 = make_float2 (_S524);
    float _S526 = (F32_min((_S521), (_S524)));
    float _S527 = (F32_max((_S522), (_S524)));
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S528 = s_primal_ctx_mul_0(R_4, pos_i_6) + t_4;
    float _S529 = _S528.z;
    float2  _S530 = make_float2 (_S529);
    float3  _S531 = pos_4 + make_float3 (0.5f * size_4);
    float2  _S532 = float2 {_S493.x, _S493.y};
    float2  _S533 = _S532 / make_float2 (_S494);
    float2  _S534 = make_float2 (_S494 * _S494);
    float u_48 = _S533.x;
    float v_48 = _S533.y;
    float r2_48 = u_48 * u_48 + v_48 * v_48;
    float _S535 = dist_coeffs_4[int(2)] + r2_48 * dist_coeffs_4[int(3)];
    float _S536 = dist_coeffs_4[int(1)] + r2_48 * _S535;
    float _S537 = dist_coeffs_4[int(0)] + r2_48 * _S536;
    float radial_16 = 1.0f + r2_48 * _S537;
    float _S538 = 2.0f * dist_coeffs_4[int(4)];
    float _S539 = _S538 * u_48;
    float _S540 = 2.0f * u_48;
    float _S541 = 2.0f * dist_coeffs_4[int(5)];
    float _S542 = _S541 * u_48;
    float _S543 = 2.0f * v_48;
    float2  _S544 = _S533 * make_float2 (radial_16) + make_float2 (_S539 * v_48 + dist_coeffs_4[int(5)] * (r2_48 + _S540 * u_48) + dist_coeffs_4[int(6)] * r2_48, _S542 * v_48 + dist_coeffs_4[int(4)] * (r2_48 + _S543 * v_48) + dist_coeffs_4[int(7)] * r2_48);
    float2  _S545 = _S544 + make_float2 (dist_coeffs_4[int(8)] * _S544.x + dist_coeffs_4[int(9)] * _S544.y, 0.0f);
    float _S546 = fx_4 * _S545.x + cx_4;
    float _S547 = fy_4 * _S545.y + cy_4;
    float2  _S548 = float2 {_S498.x, _S498.y};
    float2  _S549 = _S548 / make_float2 (_S499);
    float2  _S550 = make_float2 (_S499 * _S499);
    float u_49 = _S549.x;
    float v_49 = _S549.y;
    float r2_49 = u_49 * u_49 + v_49 * v_49;
    float _S551 = dist_coeffs_4[int(2)] + r2_49 * dist_coeffs_4[int(3)];
    float _S552 = dist_coeffs_4[int(1)] + r2_49 * _S551;
    float _S553 = dist_coeffs_4[int(0)] + r2_49 * _S552;
    float radial_17 = 1.0f + r2_49 * _S553;
    float _S554 = _S538 * u_49;
    float _S555 = 2.0f * u_49;
    float _S556 = _S541 * u_49;
    float _S557 = 2.0f * v_49;
    float2  _S558 = _S549 * make_float2 (radial_17) + make_float2 (_S554 * v_49 + dist_coeffs_4[int(5)] * (r2_49 + _S555 * u_49) + dist_coeffs_4[int(6)] * r2_49, _S556 * v_49 + dist_coeffs_4[int(4)] * (r2_49 + _S557 * v_49) + dist_coeffs_4[int(7)] * r2_49);
    float2  _S559 = _S558 + make_float2 (dist_coeffs_4[int(8)] * _S558.x + dist_coeffs_4[int(9)] * _S558.y, 0.0f);
    float _S560 = fx_4 * _S559.x + cx_4;
    float _S561 = fy_4 * _S559.y + cy_4;
    float2  _S562 = float2 {_S503.x, _S503.y};
    float2  _S563 = _S562 / make_float2 (_S504);
    float2  _S564 = make_float2 (_S504 * _S504);
    float u_50 = _S563.x;
    float v_50 = _S563.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S565 = dist_coeffs_4[int(2)] + r2_50 * dist_coeffs_4[int(3)];
    float _S566 = dist_coeffs_4[int(1)] + r2_50 * _S565;
    float _S567 = dist_coeffs_4[int(0)] + r2_50 * _S566;
    float radial_18 = 1.0f + r2_50 * _S567;
    float _S568 = _S538 * u_50;
    float _S569 = 2.0f * u_50;
    float _S570 = _S541 * u_50;
    float _S571 = 2.0f * v_50;
    float2  _S572 = _S563 * make_float2 (radial_18) + make_float2 (_S568 * v_50 + dist_coeffs_4[int(5)] * (r2_50 + _S569 * u_50) + dist_coeffs_4[int(6)] * r2_50, _S570 * v_50 + dist_coeffs_4[int(4)] * (r2_50 + _S571 * v_50) + dist_coeffs_4[int(7)] * r2_50);
    float2  _S573 = _S572 + make_float2 (dist_coeffs_4[int(8)] * _S572.x + dist_coeffs_4[int(9)] * _S572.y, 0.0f);
    float _S574 = fx_4 * _S573.x + cx_4;
    float _S575 = fy_4 * _S573.y + cy_4;
    float2  _S576 = float2 {_S508.x, _S508.y};
    float2  _S577 = _S576 / make_float2 (_S509);
    float2  _S578 = make_float2 (_S509 * _S509);
    float u_51 = _S577.x;
    float v_51 = _S577.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float _S579 = dist_coeffs_4[int(2)] + r2_51 * dist_coeffs_4[int(3)];
    float _S580 = dist_coeffs_4[int(1)] + r2_51 * _S579;
    float _S581 = dist_coeffs_4[int(0)] + r2_51 * _S580;
    float radial_19 = 1.0f + r2_51 * _S581;
    float _S582 = _S538 * u_51;
    float _S583 = 2.0f * u_51;
    float _S584 = _S541 * u_51;
    float _S585 = 2.0f * v_51;
    float2  _S586 = _S577 * make_float2 (radial_19) + make_float2 (_S582 * v_51 + dist_coeffs_4[int(5)] * (r2_51 + _S583 * u_51) + dist_coeffs_4[int(6)] * r2_51, _S584 * v_51 + dist_coeffs_4[int(4)] * (r2_51 + _S585 * v_51) + dist_coeffs_4[int(7)] * r2_51);
    float2  _S587 = _S586 + make_float2 (dist_coeffs_4[int(8)] * _S586.x + dist_coeffs_4[int(9)] * _S586.y, 0.0f);
    float _S588 = fx_4 * _S587.x + cx_4;
    float _S589 = fy_4 * _S587.y + cy_4;
    float2  _S590 = float2 {_S513.x, _S513.y};
    float2  _S591 = _S590 / make_float2 (_S514);
    float2  _S592 = make_float2 (_S514 * _S514);
    float u_52 = _S591.x;
    float v_52 = _S591.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float _S593 = dist_coeffs_4[int(2)] + r2_52 * dist_coeffs_4[int(3)];
    float _S594 = dist_coeffs_4[int(1)] + r2_52 * _S593;
    float _S595 = dist_coeffs_4[int(0)] + r2_52 * _S594;
    float radial_20 = 1.0f + r2_52 * _S595;
    float _S596 = _S538 * u_52;
    float _S597 = 2.0f * u_52;
    float _S598 = _S541 * u_52;
    float _S599 = 2.0f * v_52;
    float2  _S600 = _S591 * make_float2 (radial_20) + make_float2 (_S596 * v_52 + dist_coeffs_4[int(5)] * (r2_52 + _S597 * u_52) + dist_coeffs_4[int(6)] * r2_52, _S598 * v_52 + dist_coeffs_4[int(4)] * (r2_52 + _S599 * v_52) + dist_coeffs_4[int(7)] * r2_52);
    float2  _S601 = _S600 + make_float2 (dist_coeffs_4[int(8)] * _S600.x + dist_coeffs_4[int(9)] * _S600.y, 0.0f);
    float _S602 = fx_4 * _S601.x + cx_4;
    float _S603 = fy_4 * _S601.y + cy_4;
    float2  _S604 = float2 {_S518.x, _S518.y};
    float2  _S605 = _S604 / make_float2 (_S519);
    float2  _S606 = make_float2 (_S519 * _S519);
    float u_53 = _S605.x;
    float v_53 = _S605.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S607 = dist_coeffs_4[int(2)] + r2_53 * dist_coeffs_4[int(3)];
    float _S608 = dist_coeffs_4[int(1)] + r2_53 * _S607;
    float _S609 = dist_coeffs_4[int(0)] + r2_53 * _S608;
    float radial_21 = 1.0f + r2_53 * _S609;
    float _S610 = _S538 * u_53;
    float _S611 = 2.0f * u_53;
    float _S612 = _S541 * u_53;
    float _S613 = 2.0f * v_53;
    float2  _S614 = _S605 * make_float2 (radial_21) + make_float2 (_S610 * v_53 + dist_coeffs_4[int(5)] * (r2_53 + _S611 * u_53) + dist_coeffs_4[int(6)] * r2_53, _S612 * v_53 + dist_coeffs_4[int(4)] * (r2_53 + _S613 * v_53) + dist_coeffs_4[int(7)] * r2_53);
    float2  _S615 = _S614 + make_float2 (dist_coeffs_4[int(8)] * _S614.x + dist_coeffs_4[int(9)] * _S614.y, 0.0f);
    float _S616 = fx_4 * _S615.x + cx_4;
    float _S617 = fy_4 * _S615.y + cy_4;
    float2  _S618 = float2 {_S523.x, _S523.y};
    float2  _S619 = _S618 / make_float2 (_S524);
    float2  _S620 = make_float2 (_S524 * _S524);
    float u_54 = _S619.x;
    float v_54 = _S619.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float _S621 = dist_coeffs_4[int(2)] + r2_54 * dist_coeffs_4[int(3)];
    float _S622 = dist_coeffs_4[int(1)] + r2_54 * _S621;
    float _S623 = dist_coeffs_4[int(0)] + r2_54 * _S622;
    float radial_22 = 1.0f + r2_54 * _S623;
    float _S624 = _S538 * u_54;
    float _S625 = 2.0f * u_54;
    float _S626 = _S541 * u_54;
    float _S627 = 2.0f * v_54;
    float2  _S628 = _S619 * make_float2 (radial_22) + make_float2 (_S624 * v_54 + dist_coeffs_4[int(5)] * (r2_54 + _S625 * u_54) + dist_coeffs_4[int(6)] * r2_54, _S626 * v_54 + dist_coeffs_4[int(4)] * (r2_54 + _S627 * v_54) + dist_coeffs_4[int(7)] * r2_54);
    float2  _S629 = _S628 + make_float2 (dist_coeffs_4[int(8)] * _S628.x + dist_coeffs_4[int(9)] * _S628.y, 0.0f);
    float _S630 = fx_4 * _S629.x + cx_4;
    float _S631 = fy_4 * _S629.y + cy_4;
    float2  _S632 = float2 {_S528.x, _S528.y};
    float2  _S633 = _S632 / make_float2 (_S529);
    float2  _S634 = make_float2 (_S529 * _S529);
    float u_55 = _S633.x;
    float v_55 = _S633.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float _S635 = dist_coeffs_4[int(2)] + r2_55 * dist_coeffs_4[int(3)];
    float _S636 = dist_coeffs_4[int(1)] + r2_55 * _S635;
    float _S637 = dist_coeffs_4[int(0)] + r2_55 * _S636;
    float radial_23 = 1.0f + r2_55 * _S637;
    float _S638 = _S538 * u_55;
    float _S639 = 2.0f * u_55;
    float _S640 = _S541 * u_55;
    float _S641 = 2.0f * v_55;
    float2  _S642 = _S633 * make_float2 (radial_23) + make_float2 (_S638 * v_55 + dist_coeffs_4[int(5)] * (r2_55 + _S639 * u_55) + dist_coeffs_4[int(6)] * r2_55, _S640 * v_55 + dist_coeffs_4[int(4)] * (r2_55 + _S641 * v_55) + dist_coeffs_4[int(7)] * r2_55);
    float2  _S643 = _S642 + make_float2 (dist_coeffs_4[int(8)] * _S642.x + dist_coeffs_4[int(9)] * _S642.y, 0.0f);
    float _S644 = fx_4 * _S643.x + cx_4;
    float _S645 = fy_4 * _S643.y + cy_4;
    float _S646 = (F32_max((_S546), (_S560)));
    float _S647 = (F32_min((_S546), (_S560)));
    float _S648 = (F32_max((_S547), (_S561)));
    float _S649 = (F32_min((_S547), (_S561)));
    float _S650 = (F32_max((_S646), (_S574)));
    float _S651 = (F32_min((_S647), (_S574)));
    float _S652 = (F32_max((_S648), (_S575)));
    float _S653 = (F32_min((_S649), (_S575)));
    float _S654 = (F32_max((_S650), (_S588)));
    float _S655 = (F32_min((_S651), (_S588)));
    float _S656 = (F32_max((_S652), (_S589)));
    float _S657 = (F32_min((_S653), (_S589)));
    float _S658 = (F32_max((_S654), (_S602)));
    float _S659 = (F32_min((_S655), (_S602)));
    float _S660 = (F32_max((_S656), (_S603)));
    float _S661 = (F32_min((_S657), (_S603)));
    float _S662 = (F32_max((_S658), (_S616)));
    float _S663 = (F32_min((_S659), (_S616)));
    float _S664 = (F32_max((_S660), (_S617)));
    float _S665 = (F32_min((_S661), (_S617)));
    float _S666 = (F32_max((_S662), (_S630)));
    float _S667 = (F32_min((_S663), (_S630)));
    float _S668 = (F32_max((_S664), (_S631)));
    float _S669 = (F32_min((_S665), (_S631)));
    Matrix<float, 3, 3>  _S670 = transpose_1(R_4);
    float3  _S671 = s_primal_ctx_mul_0(R_4, _S531) + t_4 - - s_primal_ctx_mul_0(_S670, t_4);
    float _S672 = _S671.x;
    float _S673 = _S671.y;
    float _S674 = _S671.z;
    float _S675 = _S672 * _S672 + _S673 * _S673 + _S674 * _S674;
    float _S676 = s_primal_ctx_sqrt_0(_S675);
    float x_13 = _S672 / _S676;
    float3  _S677 = make_float3 (x_13);
    float _S678 = _S676 * _S676;
    float y_9 = _S673 / _S676;
    float z_6 = _S674 / _S676;
    float3  _S679 = make_float3 (z_6);
    float _S680 = - y_9;
    float3  _S681 = make_float3 (_S680);
    float z2_4 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_4 = x_13 * x_13 - y_9 * y_9;
    float _S682 = 2.0f * x_13;
    float fS1_4 = _S682 * y_9;
    float pSH6_0 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float3  _S683 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_13;
    float3  _S684 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_9;
    float3  _S685 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S686 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S687 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_6;
    float _S688 = 1.86588168144226074f * z2_4 - 1.11952900886535645f;
    float pSH12_0 = z_6 * _S688;
    float3  _S689 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_13;
    float3  _S690 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_9;
    float3  _S691 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S692 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S693 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_13 * fC1_4 - y_9 * fS1_4);
    float3  _S694 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_13 * fS1_4 + y_9 * fC1_4);
    float3  _S695 = make_float3 (pSH9_0);
    float3  _S696 = make_float3 (0.0f);
    float3  _S697 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S698;
    (&_S698)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_4[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S680) * sh_coeffs_4[int(1)] + make_float3 (z_6) * sh_coeffs_4[int(2)] - make_float3 (x_13) * sh_coeffs_4[int(3)]) + (make_float3 (pSH4_0) * sh_coeffs_4[int(4)] + make_float3 (pSH5_0) * sh_coeffs_4[int(5)] + make_float3 (pSH6_0) * sh_coeffs_4[int(6)] + make_float3 (pSH7_0) * sh_coeffs_4[int(7)] + make_float3 (pSH8_0) * sh_coeffs_4[int(8)]) + (make_float3 (pSH9_0) * sh_coeffs_4[int(9)] + make_float3 (pSH10_0) * sh_coeffs_4[int(10)] + make_float3 (pSH11_0) * sh_coeffs_4[int(11)] + make_float3 (pSH12_0) * sh_coeffs_4[int(12)] + make_float3 (pSH13_0) * sh_coeffs_4[int(13)] + make_float3 (pSH14_0) * sh_coeffs_4[int(14)] + make_float3 (pSH15_0) * sh_coeffs_4[int(15)]) + make_float3 (0.5f);
    (&_S698)->differential_0 = _S697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S699;
    (&_S699)->primal_0 = _S696;
    (&_S699)->differential_0 = _S697;
    s_bwd_prop_max_0(&_S698, &_S699, v_rgb_0);
    float3  _S700 = _S694 * _S698.differential_0;
    float3  _S701 = sh_coeffs_4[int(15)] * _S698.differential_0;
    float3  _S702 = _S692 * _S698.differential_0;
    float3  _S703 = sh_coeffs_4[int(14)] * _S698.differential_0;
    float3  _S704 = _S690 * _S698.differential_0;
    float3  _S705 = sh_coeffs_4[int(13)] * _S698.differential_0;
    float3  _S706 = _S689 * _S698.differential_0;
    float3  _S707 = sh_coeffs_4[int(12)] * _S698.differential_0;
    float3  _S708 = _S691 * _S698.differential_0;
    float3  _S709 = sh_coeffs_4[int(11)] * _S698.differential_0;
    float3  _S710 = _S693 * _S698.differential_0;
    float3  _S711 = sh_coeffs_4[int(10)] * _S698.differential_0;
    float3  _S712 = _S695 * _S698.differential_0;
    float3  _S713 = sh_coeffs_4[int(9)] * _S698.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S713.x + _S713.y + _S713.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S701.x + _S701.y + _S701.z);
    float _S714 = _S711.x + _S711.y + _S711.z;
    float _S715 = _S703.x + _S703.y + _S703.z;
    float _S716 = _S709.x + _S709.y + _S709.z;
    float _S717 = _S705.x + _S705.y + _S705.z;
    float _S718 = _S707.x + _S707.y + _S707.z;
    float _S719 = - s_diff_fC2_T_0;
    float3  _S720 = _S686 * _S698.differential_0;
    float3  _S721 = sh_coeffs_4[int(8)] * _S698.differential_0;
    float3  _S722 = _S684 * _S698.differential_0;
    float3  _S723 = sh_coeffs_4[int(7)] * _S698.differential_0;
    float3  _S724 = _S683 * _S698.differential_0;
    float3  _S725 = sh_coeffs_4[int(6)] * _S698.differential_0;
    float3  _S726 = _S685 * _S698.differential_0;
    float3  _S727 = sh_coeffs_4[int(5)] * _S698.differential_0;
    float3  _S728 = _S687 * _S698.differential_0;
    float3  _S729 = sh_coeffs_4[int(4)] * _S698.differential_0;
    float _S730 = _S727.x + _S727.y + _S727.z;
    float _S731 = _S723.x + _S723.y + _S723.z;
    float _S732 = fTmp1B_4 * _S714 + x_13 * s_diff_fS2_T_0 + y_9 * _S719 + 0.54627424478530884f * (_S729.x + _S729.y + _S729.z);
    float _S733 = fTmp1B_4 * _S715 + y_9 * s_diff_fS2_T_0 + x_13 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S721.x + _S721.y + _S721.z);
    float _S734 = y_9 * - _S733;
    float _S735 = x_13 * _S733;
    float _S736 = z_6 * (1.86588168144226074f * (z_6 * _S718) + -2.28522896766662598f * (y_9 * _S716 + x_13 * _S717) + 0.94617468118667603f * (_S725.x + _S725.y + _S725.z));
    float3  _S737 = make_float3 (0.48860251903533936f) * _S698.differential_0;
    float3  _S738 = - _S737;
    float3  _S739 = _S677 * _S738;
    float3  _S740 = sh_coeffs_4[int(3)] * _S738;
    float3  _S741 = _S679 * _S737;
    float3  _S742 = sh_coeffs_4[int(2)] * _S737;
    float3  _S743 = _S681 * _S737;
    float3  _S744 = sh_coeffs_4[int(1)] * _S737;
    float _S745 = (_S688 * _S718 + 1.44530570507049561f * (fS1_4 * _S714 + fC1_4 * _S715) + -1.09254848957061768f * (y_9 * _S730 + x_13 * _S731) + _S736 + _S736 + _S742.x + _S742.y + _S742.z) / _S678;
    float _S746 = _S676 * _S745;
    float _S747 = (fTmp0C_4 * _S716 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S719 + fTmp0B_4 * _S730 + _S682 * _S732 + _S734 + _S734 + - (_S744.x + _S744.y + _S744.z)) / _S678;
    float _S748 = _S676 * _S747;
    float _S749 = (fTmp0C_4 * _S717 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S731 + 2.0f * (y_9 * _S732) + _S735 + _S735 + _S740.x + _S740.y + _S740.z) / _S678;
    float _S750 = _S676 * _S749;
    float _S751 = _S674 * - _S745 + _S673 * - _S747 + _S672 * - _S749;
    DiffPair_float_0 _S752;
    (&_S752)->primal_0 = _S675;
    (&_S752)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S752, _S751);
    float _S753 = _S674 * _S752.differential_0;
    float _S754 = _S673 * _S752.differential_0;
    float _S755 = _S672 * _S752.differential_0;
    float3  _S756 = make_float3 (0.282094806432724f) * _S698.differential_0;
    float3  _S757 = make_float3 (_S750 + _S755 + _S755, _S748 + _S754 + _S754, _S746 + _S753 + _S753);
    float3  _S758 = - - _S757;
    Matrix<float, 3, 3>  _S759 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S760;
    (&_S760)->primal_0 = _S670;
    (&_S760)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S761;
    (&_S761)->primal_0 = t_4;
    (&_S761)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S760, &_S761, _S758);
    Matrix<float, 3, 3>  _S762 = transpose_1(_S760.differential_0);
    DiffPair_float_0 _S763;
    (&_S763)->primal_0 = _S669;
    (&_S763)->differential_0 = 0.0f;
    DiffPair_float_0 _S764;
    (&_S764)->primal_0 = _S645;
    (&_S764)->differential_0 = 0.0f;
    _d_min_0(&_S763, &_S764, 0.0f);
    DiffPair_float_0 _S765;
    (&_S765)->primal_0 = _S668;
    (&_S765)->differential_0 = 0.0f;
    DiffPair_float_0 _S766;
    (&_S766)->primal_0 = _S645;
    (&_S766)->differential_0 = 0.0f;
    _d_max_0(&_S765, &_S766, 0.0f);
    DiffPair_float_0 _S767;
    (&_S767)->primal_0 = _S667;
    (&_S767)->differential_0 = 0.0f;
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = _S644;
    (&_S768)->differential_0 = 0.0f;
    _d_min_0(&_S767, &_S768, 0.0f);
    DiffPair_float_0 _S769;
    (&_S769)->primal_0 = _S666;
    (&_S769)->differential_0 = 0.0f;
    DiffPair_float_0 _S770;
    (&_S770)->primal_0 = _S644;
    (&_S770)->differential_0 = 0.0f;
    _d_max_0(&_S769, &_S770, 0.0f);
    DiffPair_float_0 _S771;
    (&_S771)->primal_0 = _S665;
    (&_S771)->differential_0 = 0.0f;
    DiffPair_float_0 _S772;
    (&_S772)->primal_0 = _S631;
    (&_S772)->differential_0 = 0.0f;
    _d_min_0(&_S771, &_S772, _S763.differential_0);
    DiffPair_float_0 _S773;
    (&_S773)->primal_0 = _S664;
    (&_S773)->differential_0 = 0.0f;
    DiffPair_float_0 _S774;
    (&_S774)->primal_0 = _S631;
    (&_S774)->differential_0 = 0.0f;
    _d_max_0(&_S773, &_S774, _S765.differential_0);
    DiffPair_float_0 _S775;
    (&_S775)->primal_0 = _S663;
    (&_S775)->differential_0 = 0.0f;
    DiffPair_float_0 _S776;
    (&_S776)->primal_0 = _S630;
    (&_S776)->differential_0 = 0.0f;
    _d_min_0(&_S775, &_S776, _S767.differential_0);
    DiffPair_float_0 _S777;
    (&_S777)->primal_0 = _S662;
    (&_S777)->differential_0 = 0.0f;
    DiffPair_float_0 _S778;
    (&_S778)->primal_0 = _S630;
    (&_S778)->differential_0 = 0.0f;
    _d_max_0(&_S777, &_S778, _S769.differential_0);
    DiffPair_float_0 _S779;
    (&_S779)->primal_0 = _S661;
    (&_S779)->differential_0 = 0.0f;
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = _S617;
    (&_S780)->differential_0 = 0.0f;
    _d_min_0(&_S779, &_S780, _S771.differential_0);
    DiffPair_float_0 _S781;
    (&_S781)->primal_0 = _S660;
    (&_S781)->differential_0 = 0.0f;
    DiffPair_float_0 _S782;
    (&_S782)->primal_0 = _S617;
    (&_S782)->differential_0 = 0.0f;
    _d_max_0(&_S781, &_S782, _S773.differential_0);
    DiffPair_float_0 _S783;
    (&_S783)->primal_0 = _S659;
    (&_S783)->differential_0 = 0.0f;
    DiffPair_float_0 _S784;
    (&_S784)->primal_0 = _S616;
    (&_S784)->differential_0 = 0.0f;
    _d_min_0(&_S783, &_S784, _S775.differential_0);
    DiffPair_float_0 _S785;
    (&_S785)->primal_0 = _S658;
    (&_S785)->differential_0 = 0.0f;
    DiffPair_float_0 _S786;
    (&_S786)->primal_0 = _S616;
    (&_S786)->differential_0 = 0.0f;
    _d_max_0(&_S785, &_S786, _S777.differential_0);
    DiffPair_float_0 _S787;
    (&_S787)->primal_0 = _S657;
    (&_S787)->differential_0 = 0.0f;
    DiffPair_float_0 _S788;
    (&_S788)->primal_0 = _S603;
    (&_S788)->differential_0 = 0.0f;
    _d_min_0(&_S787, &_S788, _S779.differential_0);
    DiffPair_float_0 _S789;
    (&_S789)->primal_0 = _S656;
    (&_S789)->differential_0 = 0.0f;
    DiffPair_float_0 _S790;
    (&_S790)->primal_0 = _S603;
    (&_S790)->differential_0 = 0.0f;
    _d_max_0(&_S789, &_S790, _S781.differential_0);
    DiffPair_float_0 _S791;
    (&_S791)->primal_0 = _S655;
    (&_S791)->differential_0 = 0.0f;
    DiffPair_float_0 _S792;
    (&_S792)->primal_0 = _S602;
    (&_S792)->differential_0 = 0.0f;
    _d_min_0(&_S791, &_S792, _S783.differential_0);
    DiffPair_float_0 _S793;
    (&_S793)->primal_0 = _S654;
    (&_S793)->differential_0 = 0.0f;
    DiffPair_float_0 _S794;
    (&_S794)->primal_0 = _S602;
    (&_S794)->differential_0 = 0.0f;
    _d_max_0(&_S793, &_S794, _S785.differential_0);
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = _S653;
    (&_S795)->differential_0 = 0.0f;
    DiffPair_float_0 _S796;
    (&_S796)->primal_0 = _S589;
    (&_S796)->differential_0 = 0.0f;
    _d_min_0(&_S795, &_S796, _S787.differential_0);
    DiffPair_float_0 _S797;
    (&_S797)->primal_0 = _S652;
    (&_S797)->differential_0 = 0.0f;
    DiffPair_float_0 _S798;
    (&_S798)->primal_0 = _S589;
    (&_S798)->differential_0 = 0.0f;
    _d_max_0(&_S797, &_S798, _S789.differential_0);
    DiffPair_float_0 _S799;
    (&_S799)->primal_0 = _S651;
    (&_S799)->differential_0 = 0.0f;
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = _S588;
    (&_S800)->differential_0 = 0.0f;
    _d_min_0(&_S799, &_S800, _S791.differential_0);
    DiffPair_float_0 _S801;
    (&_S801)->primal_0 = _S650;
    (&_S801)->differential_0 = 0.0f;
    DiffPair_float_0 _S802;
    (&_S802)->primal_0 = _S588;
    (&_S802)->differential_0 = 0.0f;
    _d_max_0(&_S801, &_S802, _S793.differential_0);
    DiffPair_float_0 _S803;
    (&_S803)->primal_0 = _S649;
    (&_S803)->differential_0 = 0.0f;
    DiffPair_float_0 _S804;
    (&_S804)->primal_0 = _S575;
    (&_S804)->differential_0 = 0.0f;
    _d_min_0(&_S803, &_S804, _S795.differential_0);
    DiffPair_float_0 _S805;
    (&_S805)->primal_0 = _S648;
    (&_S805)->differential_0 = 0.0f;
    DiffPair_float_0 _S806;
    (&_S806)->primal_0 = _S575;
    (&_S806)->differential_0 = 0.0f;
    _d_max_0(&_S805, &_S806, _S797.differential_0);
    DiffPair_float_0 _S807;
    (&_S807)->primal_0 = _S647;
    (&_S807)->differential_0 = 0.0f;
    DiffPair_float_0 _S808;
    (&_S808)->primal_0 = _S574;
    (&_S808)->differential_0 = 0.0f;
    _d_min_0(&_S807, &_S808, _S799.differential_0);
    DiffPair_float_0 _S809;
    (&_S809)->primal_0 = _S646;
    (&_S809)->differential_0 = 0.0f;
    DiffPair_float_0 _S810;
    (&_S810)->primal_0 = _S574;
    (&_S810)->differential_0 = 0.0f;
    _d_max_0(&_S809, &_S810, _S801.differential_0);
    DiffPair_float_0 _S811;
    (&_S811)->primal_0 = _S547;
    (&_S811)->differential_0 = 0.0f;
    DiffPair_float_0 _S812;
    (&_S812)->primal_0 = _S561;
    (&_S812)->differential_0 = 0.0f;
    _d_min_0(&_S811, &_S812, _S803.differential_0);
    DiffPair_float_0 _S813;
    (&_S813)->primal_0 = _S547;
    (&_S813)->differential_0 = 0.0f;
    DiffPair_float_0 _S814;
    (&_S814)->primal_0 = _S561;
    (&_S814)->differential_0 = 0.0f;
    _d_max_0(&_S813, &_S814, _S805.differential_0);
    DiffPair_float_0 _S815;
    (&_S815)->primal_0 = _S546;
    (&_S815)->differential_0 = 0.0f;
    DiffPair_float_0 _S816;
    (&_S816)->primal_0 = _S560;
    (&_S816)->differential_0 = 0.0f;
    _d_min_0(&_S815, &_S816, _S807.differential_0);
    DiffPair_float_0 _S817;
    (&_S817)->primal_0 = _S546;
    (&_S817)->differential_0 = 0.0f;
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = _S560;
    (&_S818)->differential_0 = 0.0f;
    _d_max_0(&_S817, &_S818, _S809.differential_0);
    float _S819 = fx_4 * (_S768.differential_0 + _S770.differential_0);
    float2  _S820 = make_float2 (_S819, fy_4 * (_S764.differential_0 + _S766.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S819, dist_coeffs_4[int(9)] * _S819);
    float2  _S821 = _S633 * _S820;
    float _S822 = dist_coeffs_4[int(4)] * _S820.y;
    float _S823 = dist_coeffs_4[int(5)] * _S820.x;
    float _S824 = _S821.x + _S821.y;
    float _S825 = r2_55 * _S824;
    float _S826 = r2_55 * _S825;
    float _S827 = dist_coeffs_4[int(7)] * _S820.y + _S822 + dist_coeffs_4[int(6)] * _S820.x + _S823 + _S637 * _S824 + _S636 * _S825 + _S635 * _S826 + dist_coeffs_4[int(3)] * (r2_55 * _S826);
    float _S828 = v_55 * _S827;
    float _S829 = u_55 * _S827;
    float2  _S830 = (make_float2 (radial_23) * _S820 + make_float2 (_S541 * (v_55 * _S820.y) + _S639 * _S823 + 2.0f * (u_55 * _S823) + _S538 * (v_55 * _S820.x) + _S829 + _S829, _S641 * _S822 + 2.0f * (v_55 * _S822) + _S640 * _S820.y + _S638 * _S820.x + _S828 + _S828)) / _S634;
    float2  _S831 = _S632 * - _S830;
    float2  _S832 = _S530 * _S830;
    float _S833 = fx_4 * (_S776.differential_0 + _S778.differential_0);
    float2  _S834 = make_float2 (_S833, fy_4 * (_S772.differential_0 + _S774.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S833, dist_coeffs_4[int(9)] * _S833);
    float2  _S835 = _S619 * _S834;
    float _S836 = dist_coeffs_4[int(4)] * _S834.y;
    float _S837 = dist_coeffs_4[int(5)] * _S834.x;
    float _S838 = _S835.x + _S835.y;
    float _S839 = r2_54 * _S838;
    float _S840 = r2_54 * _S839;
    float _S841 = dist_coeffs_4[int(7)] * _S834.y + _S836 + dist_coeffs_4[int(6)] * _S834.x + _S837 + _S623 * _S838 + _S622 * _S839 + _S621 * _S840 + dist_coeffs_4[int(3)] * (r2_54 * _S840);
    float _S842 = v_54 * _S841;
    float _S843 = u_54 * _S841;
    float2  _S844 = (make_float2 (radial_22) * _S834 + make_float2 (_S541 * (v_54 * _S834.y) + _S625 * _S837 + 2.0f * (u_54 * _S837) + _S538 * (v_54 * _S834.x) + _S843 + _S843, _S627 * _S836 + 2.0f * (v_54 * _S836) + _S626 * _S834.y + _S624 * _S834.x + _S842 + _S842)) / _S620;
    float2  _S845 = _S618 * - _S844;
    float2  _S846 = _S525 * _S844;
    float _S847 = fx_4 * (_S784.differential_0 + _S786.differential_0);
    float2  _S848 = make_float2 (_S847, fy_4 * (_S780.differential_0 + _S782.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S847, dist_coeffs_4[int(9)] * _S847);
    float2  _S849 = _S605 * _S848;
    float _S850 = dist_coeffs_4[int(4)] * _S848.y;
    float _S851 = dist_coeffs_4[int(5)] * _S848.x;
    float _S852 = _S849.x + _S849.y;
    float _S853 = r2_53 * _S852;
    float _S854 = r2_53 * _S853;
    float _S855 = dist_coeffs_4[int(7)] * _S848.y + _S850 + dist_coeffs_4[int(6)] * _S848.x + _S851 + _S609 * _S852 + _S608 * _S853 + _S607 * _S854 + dist_coeffs_4[int(3)] * (r2_53 * _S854);
    float _S856 = v_53 * _S855;
    float _S857 = u_53 * _S855;
    float2  _S858 = (make_float2 (radial_21) * _S848 + make_float2 (_S541 * (v_53 * _S848.y) + _S611 * _S851 + 2.0f * (u_53 * _S851) + _S538 * (v_53 * _S848.x) + _S857 + _S857, _S613 * _S850 + 2.0f * (v_53 * _S850) + _S612 * _S848.y + _S610 * _S848.x + _S856 + _S856)) / _S606;
    float2  _S859 = _S604 * - _S858;
    float2  _S860 = _S520 * _S858;
    float _S861 = fx_4 * (_S792.differential_0 + _S794.differential_0);
    float2  _S862 = make_float2 (_S861, fy_4 * (_S788.differential_0 + _S790.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S861, dist_coeffs_4[int(9)] * _S861);
    float2  _S863 = _S591 * _S862;
    float _S864 = dist_coeffs_4[int(4)] * _S862.y;
    float _S865 = dist_coeffs_4[int(5)] * _S862.x;
    float _S866 = _S863.x + _S863.y;
    float _S867 = r2_52 * _S866;
    float _S868 = r2_52 * _S867;
    float _S869 = dist_coeffs_4[int(7)] * _S862.y + _S864 + dist_coeffs_4[int(6)] * _S862.x + _S865 + _S595 * _S866 + _S594 * _S867 + _S593 * _S868 + dist_coeffs_4[int(3)] * (r2_52 * _S868);
    float _S870 = v_52 * _S869;
    float _S871 = u_52 * _S869;
    float2  _S872 = (make_float2 (radial_20) * _S862 + make_float2 (_S541 * (v_52 * _S862.y) + _S597 * _S865 + 2.0f * (u_52 * _S865) + _S538 * (v_52 * _S862.x) + _S871 + _S871, _S599 * _S864 + 2.0f * (v_52 * _S864) + _S598 * _S862.y + _S596 * _S862.x + _S870 + _S870)) / _S592;
    float2  _S873 = _S590 * - _S872;
    float2  _S874 = _S515 * _S872;
    float _S875 = fx_4 * (_S800.differential_0 + _S802.differential_0);
    float2  _S876 = make_float2 (_S875, fy_4 * (_S796.differential_0 + _S798.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S875, dist_coeffs_4[int(9)] * _S875);
    float2  _S877 = _S577 * _S876;
    float _S878 = dist_coeffs_4[int(4)] * _S876.y;
    float _S879 = dist_coeffs_4[int(5)] * _S876.x;
    float _S880 = _S877.x + _S877.y;
    float _S881 = r2_51 * _S880;
    float _S882 = r2_51 * _S881;
    float _S883 = dist_coeffs_4[int(7)] * _S876.y + _S878 + dist_coeffs_4[int(6)] * _S876.x + _S879 + _S581 * _S880 + _S580 * _S881 + _S579 * _S882 + dist_coeffs_4[int(3)] * (r2_51 * _S882);
    float _S884 = v_51 * _S883;
    float _S885 = u_51 * _S883;
    float2  _S886 = (make_float2 (radial_19) * _S876 + make_float2 (_S541 * (v_51 * _S876.y) + _S583 * _S879 + 2.0f * (u_51 * _S879) + _S538 * (v_51 * _S876.x) + _S885 + _S885, _S585 * _S878 + 2.0f * (v_51 * _S878) + _S584 * _S876.y + _S582 * _S876.x + _S884 + _S884)) / _S578;
    float2  _S887 = _S576 * - _S886;
    float2  _S888 = _S510 * _S886;
    float _S889 = fx_4 * (_S808.differential_0 + _S810.differential_0);
    float2  _S890 = make_float2 (_S889, fy_4 * (_S804.differential_0 + _S806.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S889, dist_coeffs_4[int(9)] * _S889);
    float2  _S891 = _S563 * _S890;
    float _S892 = dist_coeffs_4[int(4)] * _S890.y;
    float _S893 = dist_coeffs_4[int(5)] * _S890.x;
    float _S894 = _S891.x + _S891.y;
    float _S895 = r2_50 * _S894;
    float _S896 = r2_50 * _S895;
    float _S897 = dist_coeffs_4[int(7)] * _S890.y + _S892 + dist_coeffs_4[int(6)] * _S890.x + _S893 + _S567 * _S894 + _S566 * _S895 + _S565 * _S896 + dist_coeffs_4[int(3)] * (r2_50 * _S896);
    float _S898 = v_50 * _S897;
    float _S899 = u_50 * _S897;
    float2  _S900 = (make_float2 (radial_18) * _S890 + make_float2 (_S541 * (v_50 * _S890.y) + _S569 * _S893 + 2.0f * (u_50 * _S893) + _S538 * (v_50 * _S890.x) + _S899 + _S899, _S571 * _S892 + 2.0f * (v_50 * _S892) + _S570 * _S890.y + _S568 * _S890.x + _S898 + _S898)) / _S564;
    float2  _S901 = _S562 * - _S900;
    float2  _S902 = _S505 * _S900;
    float _S903 = fx_4 * (_S816.differential_0 + _S818.differential_0);
    float2  _S904 = make_float2 (_S903, fy_4 * (_S812.differential_0 + _S814.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S903, dist_coeffs_4[int(9)] * _S903);
    float2  _S905 = _S549 * _S904;
    float _S906 = dist_coeffs_4[int(4)] * _S904.y;
    float _S907 = dist_coeffs_4[int(5)] * _S904.x;
    float _S908 = _S905.x + _S905.y;
    float _S909 = r2_49 * _S908;
    float _S910 = r2_49 * _S909;
    float _S911 = dist_coeffs_4[int(7)] * _S904.y + _S906 + dist_coeffs_4[int(6)] * _S904.x + _S907 + _S553 * _S908 + _S552 * _S909 + _S551 * _S910 + dist_coeffs_4[int(3)] * (r2_49 * _S910);
    float _S912 = v_49 * _S911;
    float _S913 = u_49 * _S911;
    float2  _S914 = (make_float2 (radial_17) * _S904 + make_float2 (_S541 * (v_49 * _S904.y) + _S555 * _S907 + 2.0f * (u_49 * _S907) + _S538 * (v_49 * _S904.x) + _S913 + _S913, _S557 * _S906 + 2.0f * (v_49 * _S906) + _S556 * _S904.y + _S554 * _S904.x + _S912 + _S912)) / _S550;
    float2  _S915 = _S548 * - _S914;
    float2  _S916 = _S500 * _S914;
    float _S917 = fx_4 * (_S815.differential_0 + _S817.differential_0);
    float2  _S918 = make_float2 (_S917, fy_4 * (_S811.differential_0 + _S813.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S917, dist_coeffs_4[int(9)] * _S917);
    float2  _S919 = _S533 * _S918;
    float _S920 = dist_coeffs_4[int(4)] * _S918.y;
    float _S921 = dist_coeffs_4[int(5)] * _S918.x;
    float _S922 = _S919.x + _S919.y;
    float _S923 = r2_48 * _S922;
    float _S924 = r2_48 * _S923;
    float _S925 = dist_coeffs_4[int(7)] * _S918.y + _S920 + dist_coeffs_4[int(6)] * _S918.x + _S921 + _S537 * _S922 + _S536 * _S923 + _S535 * _S924 + dist_coeffs_4[int(3)] * (r2_48 * _S924);
    float _S926 = v_48 * _S925;
    float _S927 = u_48 * _S925;
    float2  _S928 = (make_float2 (radial_16) * _S918 + make_float2 (_S541 * (v_48 * _S918.y) + _S540 * _S921 + 2.0f * (u_48 * _S921) + _S538 * (v_48 * _S918.x) + _S927 + _S927, _S543 * _S920 + 2.0f * (v_48 * _S920) + _S542 * _S918.y + _S539 * _S918.x + _S926 + _S926)) / _S534;
    float2  _S929 = _S532 * - _S928;
    float2  _S930 = _S495 * _S928;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S931;
    (&_S931)->primal_0 = R_4;
    (&_S931)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S932;
    (&_S932)->primal_0 = _S531;
    (&_S932)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S931, &_S932, _S757);
    DiffPair_float_0 _S933;
    (&_S933)->primal_0 = _S527;
    (&_S933)->differential_0 = 0.0f;
    DiffPair_float_0 _S934;
    (&_S934)->primal_0 = _S529;
    (&_S934)->differential_0 = 0.0f;
    _d_max_0(&_S933, &_S934, 0.0f);
    DiffPair_float_0 _S935;
    (&_S935)->primal_0 = _S526;
    (&_S935)->differential_0 = 0.0f;
    DiffPair_float_0 _S936;
    (&_S936)->primal_0 = _S529;
    (&_S936)->differential_0 = 0.0f;
    _d_min_0(&_S935, &_S936, 0.0f);
    float3  _S937 = make_float3 (_S832.x, _S832.y, _S934.differential_0 + _S936.differential_0 + _S831.x + _S831.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S938;
    (&_S938)->primal_0 = R_4;
    (&_S938)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S939;
    (&_S939)->primal_0 = pos_i_6;
    (&_S939)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S938, &_S939, _S937);
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S522;
    (&_S940)->differential_0 = 0.0f;
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = _S524;
    (&_S941)->differential_0 = 0.0f;
    _d_max_0(&_S940, &_S941, _S933.differential_0);
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = _S521;
    (&_S942)->differential_0 = 0.0f;
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S524;
    (&_S943)->differential_0 = 0.0f;
    _d_min_0(&_S942, &_S943, _S935.differential_0);
    float3  _S944 = make_float3 (_S846.x, _S846.y, _S941.differential_0 + _S943.differential_0 + _S845.x + _S845.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S945;
    (&_S945)->primal_0 = R_4;
    (&_S945)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S946;
    (&_S946)->primal_0 = pos_i_5;
    (&_S946)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S945, &_S946, _S944);
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S517;
    (&_S947)->differential_0 = 0.0f;
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = _S519;
    (&_S948)->differential_0 = 0.0f;
    _d_max_0(&_S947, &_S948, _S940.differential_0);
    DiffPair_float_0 _S949;
    (&_S949)->primal_0 = _S516;
    (&_S949)->differential_0 = 0.0f;
    DiffPair_float_0 _S950;
    (&_S950)->primal_0 = _S519;
    (&_S950)->differential_0 = 0.0f;
    _d_min_0(&_S949, &_S950, _S942.differential_0);
    float3  _S951 = make_float3 (_S860.x, _S860.y, _S948.differential_0 + _S950.differential_0 + _S859.x + _S859.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S952;
    (&_S952)->primal_0 = R_4;
    (&_S952)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S953;
    (&_S953)->primal_0 = pos_i_4;
    (&_S953)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S952, &_S953, _S951);
    DiffPair_float_0 _S954;
    (&_S954)->primal_0 = _S512;
    (&_S954)->differential_0 = 0.0f;
    DiffPair_float_0 _S955;
    (&_S955)->primal_0 = _S514;
    (&_S955)->differential_0 = 0.0f;
    _d_max_0(&_S954, &_S955, _S947.differential_0);
    DiffPair_float_0 _S956;
    (&_S956)->primal_0 = _S511;
    (&_S956)->differential_0 = 0.0f;
    DiffPair_float_0 _S957;
    (&_S957)->primal_0 = _S514;
    (&_S957)->differential_0 = 0.0f;
    _d_min_0(&_S956, &_S957, _S949.differential_0);
    float3  _S958 = make_float3 (_S874.x, _S874.y, _S955.differential_0 + _S957.differential_0 + _S873.x + _S873.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S959;
    (&_S959)->primal_0 = R_4;
    (&_S959)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S960;
    (&_S960)->primal_0 = pos_i_3;
    (&_S960)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S959, &_S960, _S958);
    DiffPair_float_0 _S961;
    (&_S961)->primal_0 = _S507;
    (&_S961)->differential_0 = 0.0f;
    DiffPair_float_0 _S962;
    (&_S962)->primal_0 = _S509;
    (&_S962)->differential_0 = 0.0f;
    _d_max_0(&_S961, &_S962, _S954.differential_0);
    DiffPair_float_0 _S963;
    (&_S963)->primal_0 = _S506;
    (&_S963)->differential_0 = 0.0f;
    DiffPair_float_0 _S964;
    (&_S964)->primal_0 = _S509;
    (&_S964)->differential_0 = 0.0f;
    _d_min_0(&_S963, &_S964, _S956.differential_0);
    float3  _S965 = make_float3 (_S888.x, _S888.y, _S962.differential_0 + _S964.differential_0 + _S887.x + _S887.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S966;
    (&_S966)->primal_0 = R_4;
    (&_S966)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S967;
    (&_S967)->primal_0 = pos_i_2;
    (&_S967)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S966, &_S967, _S965);
    DiffPair_float_0 _S968;
    (&_S968)->primal_0 = _S502;
    (&_S968)->differential_0 = 0.0f;
    DiffPair_float_0 _S969;
    (&_S969)->primal_0 = _S504;
    (&_S969)->differential_0 = 0.0f;
    _d_max_0(&_S968, &_S969, _S961.differential_0);
    DiffPair_float_0 _S970;
    (&_S970)->primal_0 = _S501;
    (&_S970)->differential_0 = 0.0f;
    DiffPair_float_0 _S971;
    (&_S971)->primal_0 = _S504;
    (&_S971)->differential_0 = 0.0f;
    _d_min_0(&_S970, &_S971, _S963.differential_0);
    float3  _S972 = make_float3 (_S902.x, _S902.y, _S969.differential_0 + _S971.differential_0 + _S901.x + _S901.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S973;
    (&_S973)->primal_0 = R_4;
    (&_S973)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S974;
    (&_S974)->primal_0 = pos_i_1;
    (&_S974)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S973, &_S974, _S972);
    DiffPair_float_0 _S975;
    (&_S975)->primal_0 = _S497;
    (&_S975)->differential_0 = 0.0f;
    DiffPair_float_0 _S976;
    (&_S976)->primal_0 = _S499;
    (&_S976)->differential_0 = 0.0f;
    _d_max_0(&_S975, &_S976, _S968.differential_0);
    DiffPair_float_0 _S977;
    (&_S977)->primal_0 = _S496;
    (&_S977)->differential_0 = 0.0f;
    DiffPair_float_0 _S978;
    (&_S978)->primal_0 = _S499;
    (&_S978)->differential_0 = 0.0f;
    _d_min_0(&_S977, &_S978, _S970.differential_0);
    float3  _S979 = make_float3 (_S916.x, _S916.y, _S976.differential_0 + _S978.differential_0 + _S915.x + _S915.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S980;
    (&_S980)->primal_0 = R_4;
    (&_S980)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S981;
    (&_S981)->primal_0 = pos_i_0;
    (&_S981)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S980, &_S981, _S979);
    DiffPair_float_0 _S982;
    (&_S982)->primal_0 = 0.0f;
    (&_S982)->differential_0 = 0.0f;
    DiffPair_float_0 _S983;
    (&_S983)->primal_0 = _S494;
    (&_S983)->differential_0 = 0.0f;
    _d_max_0(&_S982, &_S983, _S975.differential_0);
    DiffPair_float_0 _S984;
    (&_S984)->primal_0 = 1.00000001504746622e+30f;
    (&_S984)->differential_0 = 0.0f;
    DiffPair_float_0 _S985;
    (&_S985)->primal_0 = _S494;
    (&_S985)->differential_0 = 0.0f;
    _d_min_0(&_S984, &_S985, _S977.differential_0);
    float3  _S986 = make_float3 (_S930.x, _S930.y, _S983.differential_0 + _S985.differential_0 + _S929.x + _S929.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S987;
    (&_S987)->primal_0 = R_4;
    (&_S987)->differential_0 = _S759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S988;
    (&_S988)->primal_0 = pos_4;
    (&_S988)->differential_0 = _S697;
    s_bwd_prop_mul_0(&_S987, &_S988, _S986);
    float3  _S989 = _S761.differential_0 + _S757 + _S937 + _S944 + _S951 + _S958 + _S965 + _S972 + _S979 + _S986;
    Matrix<float, 3, 3>  _S990 = _S762 + _S931.differential_0 + _S938.differential_0 + _S945.differential_0 + _S952.differential_0 + _S959.differential_0 + _S966.differential_0 + _S973.differential_0 + _S980.differential_0 + _S987.differential_0;
    FixedArray<float3 , 16>  _S991;
    _S991[int(0)] = _S697;
    _S991[int(1)] = _S697;
    _S991[int(2)] = _S697;
    _S991[int(3)] = _S697;
    _S991[int(4)] = _S697;
    _S991[int(5)] = _S697;
    _S991[int(6)] = _S697;
    _S991[int(7)] = _S697;
    _S991[int(8)] = _S697;
    _S991[int(9)] = _S697;
    _S991[int(10)] = _S697;
    _S991[int(11)] = _S697;
    _S991[int(12)] = _S697;
    _S991[int(13)] = _S697;
    _S991[int(14)] = _S697;
    _S991[int(15)] = _S697;
    _S991[int(15)] = _S700;
    _S991[int(14)] = _S702;
    _S991[int(13)] = _S704;
    _S991[int(12)] = _S706;
    _S991[int(11)] = _S708;
    _S991[int(10)] = _S710;
    _S991[int(9)] = _S712;
    _S991[int(8)] = _S720;
    _S991[int(7)] = _S722;
    _S991[int(6)] = _S724;
    _S991[int(5)] = _S726;
    _S991[int(4)] = _S728;
    _S991[int(3)] = _S739;
    _S991[int(2)] = _S741;
    _S991[int(1)] = _S743;
    _S991[int(0)] = _S756;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_0 = _S991;
    *v_R_0 = _S990;
    *v_t_0 = _S989;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S992, float _S993)
{
    return (F32_atan2((_S992), (_S993)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S994, DiffPair_float_0 * _S995, float _S996)
{
    _d_atan2_0(_S994, _S995, _S996);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_6, float _s_dOut_0)
{
    float _S997 = (*dpx_6).primal_0.x;
    float _S998 = (*dpx_6).primal_0.y;
    DiffPair_float_0 _S999;
    (&_S999)->primal_0 = _S997 * _S997 + _S998 * _S998;
    (&_S999)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S999, _s_dOut_0);
    float _S1000 = (*dpx_6).primal_0.y * _S999.differential_0;
    float _S1001 = _S1000 + _S1000;
    float _S1002 = (*dpx_6).primal_0.x * _S999.differential_0;
    float _S1003 = _S1002 + _S1002;
    float2  _S1004 = make_float2 (0.0f);
    *&((&_S1004)->y) = _S1001;
    *&((&_S1004)->x) = _S1003;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S1004;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1005, float _S1006)
{
    s_bwd_prop_length_impl_0(_S1005, _S1006);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_1)
{
    float _S1007 = (*dpx_7).primal_0.x;
    float _S1008 = (*dpx_7).primal_0.y;
    float _S1009 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S1010;
    (&_S1010)->primal_0 = _S1007 * _S1007 + _S1008 * _S1008 + _S1009 * _S1009;
    (&_S1010)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1010, _s_dOut_1);
    float _S1011 = (*dpx_7).primal_0.z * _S1010.differential_0;
    float _S1012 = _S1011 + _S1011;
    float _S1013 = (*dpx_7).primal_0.y * _S1010.differential_0;
    float _S1014 = _S1013 + _S1013;
    float _S1015 = (*dpx_7).primal_0.x * _S1010.differential_0;
    float _S1016 = _S1015 + _S1015;
    float3  _S1017 = make_float3 (0.0f);
    *&((&_S1017)->z) = _S1012;
    *&((&_S1017)->y) = _S1014;
    *&((&_S1017)->x) = _S1016;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S1017;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1018, float _S1019)
{
    s_bwd_prop_length_impl_1(_S1018, _S1019);
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  densities_5, FixedArray<float3 , 16>  sh_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float3  v_rgb_1, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  _S1020 = s_primal_ctx_mul_0(R_5, pos_5) + t_5;
    float _S1021 = length_1(_S1020);
    float _S1022 = (F32_min((1.00000001504746622e+30f), (_S1021)));
    float _S1023 = (F32_max((0.0f), (_S1021)));
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S1024 = s_primal_ctx_mul_0(R_5, pos_i_7) + t_5;
    float _S1025 = length_1(_S1024);
    float _S1026 = (F32_min((_S1022), (_S1025)));
    float _S1027 = (F32_max((_S1023), (_S1025)));
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S1028 = s_primal_ctx_mul_0(R_5, pos_i_8) + t_5;
    float _S1029 = length_1(_S1028);
    float _S1030 = (F32_min((_S1026), (_S1029)));
    float _S1031 = (F32_max((_S1027), (_S1029)));
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S1032 = s_primal_ctx_mul_0(R_5, pos_i_9) + t_5;
    float _S1033 = length_1(_S1032);
    float _S1034 = (F32_min((_S1030), (_S1033)));
    float _S1035 = (F32_max((_S1031), (_S1033)));
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S1036 = s_primal_ctx_mul_0(R_5, pos_i_10) + t_5;
    float _S1037 = length_1(_S1036);
    float _S1038 = (F32_min((_S1034), (_S1037)));
    float _S1039 = (F32_max((_S1035), (_S1037)));
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S1040 = s_primal_ctx_mul_0(R_5, pos_i_11) + t_5;
    float _S1041 = length_1(_S1040);
    float _S1042 = (F32_min((_S1038), (_S1041)));
    float _S1043 = (F32_max((_S1039), (_S1041)));
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S1044 = s_primal_ctx_mul_0(R_5, pos_i_12) + t_5;
    float _S1045 = length_1(_S1044);
    float _S1046 = (F32_min((_S1042), (_S1045)));
    float _S1047 = (F32_max((_S1043), (_S1045)));
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S1048 = s_primal_ctx_mul_0(R_5, pos_i_13) + t_5;
    float _S1049 = length_1(_S1048);
    float3  _S1050 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_4 = s_primal_ctx_mul_0(R_5, _S1050) + t_5;
    float2  _S1051 = float2 {_S1020.x, _S1020.y};
    float _S1052 = length_0(_S1051);
    float _S1053 = _S1020.z;
    float _S1054 = s_primal_ctx_atan2_0(_S1052, _S1053);
    bool _S1055 = _S1054 < 0.00100000004749745f;
    float k_2;
    float _S1056;
    float _S1057;
    float _S1058;
    if(_S1055)
    {
        float _S1059 = 1.0f - _S1054 * _S1054 / 3.0f;
        float _S1060 = _S1053 * _S1053;
        k_2 = _S1059 / _S1053;
        _S1056 = _S1060;
        _S1057 = _S1059;
        _S1058 = 0.0f;
    }
    else
    {
        float _S1061 = _S1052 * _S1052;
        k_2 = _S1054 / _S1052;
        _S1056 = 0.0f;
        _S1057 = 0.0f;
        _S1058 = _S1061;
    }
    float2  _S1062 = make_float2 (k_2);
    float2  _S1063 = _S1051 * make_float2 (k_2);
    float u_56 = _S1063.x;
    float v_56 = _S1063.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S1064 = dist_coeffs_5[int(2)] + r2_56 * dist_coeffs_5[int(3)];
    float _S1065 = dist_coeffs_5[int(1)] + r2_56 * _S1064;
    float _S1066 = dist_coeffs_5[int(0)] + r2_56 * _S1065;
    float radial_24 = 1.0f + r2_56 * _S1066;
    float _S1067 = 2.0f * dist_coeffs_5[int(4)];
    float _S1068 = _S1067 * u_56;
    float _S1069 = 2.0f * u_56;
    float _S1070 = 2.0f * dist_coeffs_5[int(5)];
    float _S1071 = _S1070 * u_56;
    float _S1072 = 2.0f * v_56;
    float2  _S1073 = _S1063 * make_float2 (radial_24) + make_float2 (_S1068 * v_56 + dist_coeffs_5[int(5)] * (r2_56 + _S1069 * u_56) + dist_coeffs_5[int(6)] * r2_56, _S1071 * v_56 + dist_coeffs_5[int(4)] * (r2_56 + _S1072 * v_56) + dist_coeffs_5[int(7)] * r2_56);
    float2  _S1074 = _S1073 + make_float2 (dist_coeffs_5[int(8)] * _S1073.x + dist_coeffs_5[int(9)] * _S1073.y, 0.0f);
    float _S1075 = fx_5 * _S1074.x + cx_5;
    float _S1076 = fy_5 * _S1074.y + cy_5;
    float2  _S1077 = float2 {_S1024.x, _S1024.y};
    float _S1078 = length_0(_S1077);
    float _S1079 = _S1024.z;
    float _S1080 = s_primal_ctx_atan2_0(_S1078, _S1079);
    bool _S1081 = _S1080 < 0.00100000004749745f;
    float _S1082;
    float _S1083;
    float _S1084;
    if(_S1081)
    {
        float _S1085 = 1.0f - _S1080 * _S1080 / 3.0f;
        float _S1086 = _S1079 * _S1079;
        k_2 = _S1085 / _S1079;
        _S1082 = _S1086;
        _S1083 = _S1085;
        _S1084 = 0.0f;
    }
    else
    {
        float _S1087 = _S1078 * _S1078;
        k_2 = _S1080 / _S1078;
        _S1082 = 0.0f;
        _S1083 = 0.0f;
        _S1084 = _S1087;
    }
    float2  _S1088 = make_float2 (k_2);
    float2  _S1089 = _S1077 * make_float2 (k_2);
    float u_57 = _S1089.x;
    float v_57 = _S1089.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S1090 = dist_coeffs_5[int(2)] + r2_57 * dist_coeffs_5[int(3)];
    float _S1091 = dist_coeffs_5[int(1)] + r2_57 * _S1090;
    float _S1092 = dist_coeffs_5[int(0)] + r2_57 * _S1091;
    float radial_25 = 1.0f + r2_57 * _S1092;
    float _S1093 = _S1067 * u_57;
    float _S1094 = 2.0f * u_57;
    float _S1095 = _S1070 * u_57;
    float _S1096 = 2.0f * v_57;
    float2  _S1097 = _S1089 * make_float2 (radial_25) + make_float2 (_S1093 * v_57 + dist_coeffs_5[int(5)] * (r2_57 + _S1094 * u_57) + dist_coeffs_5[int(6)] * r2_57, _S1095 * v_57 + dist_coeffs_5[int(4)] * (r2_57 + _S1096 * v_57) + dist_coeffs_5[int(7)] * r2_57);
    float2  _S1098 = _S1097 + make_float2 (dist_coeffs_5[int(8)] * _S1097.x + dist_coeffs_5[int(9)] * _S1097.y, 0.0f);
    float _S1099 = fx_5 * _S1098.x + cx_5;
    float _S1100 = fy_5 * _S1098.y + cy_5;
    float2  _S1101 = float2 {_S1028.x, _S1028.y};
    float _S1102 = length_0(_S1101);
    float _S1103 = _S1028.z;
    float _S1104 = s_primal_ctx_atan2_0(_S1102, _S1103);
    bool _S1105 = _S1104 < 0.00100000004749745f;
    float _S1106;
    float _S1107;
    float _S1108;
    if(_S1105)
    {
        float _S1109 = 1.0f - _S1104 * _S1104 / 3.0f;
        float _S1110 = _S1103 * _S1103;
        k_2 = _S1109 / _S1103;
        _S1106 = _S1110;
        _S1107 = _S1109;
        _S1108 = 0.0f;
    }
    else
    {
        float _S1111 = _S1102 * _S1102;
        k_2 = _S1104 / _S1102;
        _S1106 = 0.0f;
        _S1107 = 0.0f;
        _S1108 = _S1111;
    }
    float2  _S1112 = make_float2 (k_2);
    float2  _S1113 = _S1101 * make_float2 (k_2);
    float u_58 = _S1113.x;
    float v_58 = _S1113.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S1114 = dist_coeffs_5[int(2)] + r2_58 * dist_coeffs_5[int(3)];
    float _S1115 = dist_coeffs_5[int(1)] + r2_58 * _S1114;
    float _S1116 = dist_coeffs_5[int(0)] + r2_58 * _S1115;
    float radial_26 = 1.0f + r2_58 * _S1116;
    float _S1117 = _S1067 * u_58;
    float _S1118 = 2.0f * u_58;
    float _S1119 = _S1070 * u_58;
    float _S1120 = 2.0f * v_58;
    float2  _S1121 = _S1113 * make_float2 (radial_26) + make_float2 (_S1117 * v_58 + dist_coeffs_5[int(5)] * (r2_58 + _S1118 * u_58) + dist_coeffs_5[int(6)] * r2_58, _S1119 * v_58 + dist_coeffs_5[int(4)] * (r2_58 + _S1120 * v_58) + dist_coeffs_5[int(7)] * r2_58);
    float2  _S1122 = _S1121 + make_float2 (dist_coeffs_5[int(8)] * _S1121.x + dist_coeffs_5[int(9)] * _S1121.y, 0.0f);
    float _S1123 = fx_5 * _S1122.x + cx_5;
    float _S1124 = fy_5 * _S1122.y + cy_5;
    float2  _S1125 = float2 {_S1032.x, _S1032.y};
    float _S1126 = length_0(_S1125);
    float _S1127 = _S1032.z;
    float _S1128 = s_primal_ctx_atan2_0(_S1126, _S1127);
    bool _S1129 = _S1128 < 0.00100000004749745f;
    float _S1130;
    float _S1131;
    float _S1132;
    if(_S1129)
    {
        float _S1133 = 1.0f - _S1128 * _S1128 / 3.0f;
        float _S1134 = _S1127 * _S1127;
        k_2 = _S1133 / _S1127;
        _S1130 = _S1134;
        _S1131 = _S1133;
        _S1132 = 0.0f;
    }
    else
    {
        float _S1135 = _S1126 * _S1126;
        k_2 = _S1128 / _S1126;
        _S1130 = 0.0f;
        _S1131 = 0.0f;
        _S1132 = _S1135;
    }
    float2  _S1136 = make_float2 (k_2);
    float2  _S1137 = _S1125 * make_float2 (k_2);
    float u_59 = _S1137.x;
    float v_59 = _S1137.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S1138 = dist_coeffs_5[int(2)] + r2_59 * dist_coeffs_5[int(3)];
    float _S1139 = dist_coeffs_5[int(1)] + r2_59 * _S1138;
    float _S1140 = dist_coeffs_5[int(0)] + r2_59 * _S1139;
    float radial_27 = 1.0f + r2_59 * _S1140;
    float _S1141 = _S1067 * u_59;
    float _S1142 = 2.0f * u_59;
    float _S1143 = _S1070 * u_59;
    float _S1144 = 2.0f * v_59;
    float2  _S1145 = _S1137 * make_float2 (radial_27) + make_float2 (_S1141 * v_59 + dist_coeffs_5[int(5)] * (r2_59 + _S1142 * u_59) + dist_coeffs_5[int(6)] * r2_59, _S1143 * v_59 + dist_coeffs_5[int(4)] * (r2_59 + _S1144 * v_59) + dist_coeffs_5[int(7)] * r2_59);
    float2  _S1146 = _S1145 + make_float2 (dist_coeffs_5[int(8)] * _S1145.x + dist_coeffs_5[int(9)] * _S1145.y, 0.0f);
    float _S1147 = fx_5 * _S1146.x + cx_5;
    float _S1148 = fy_5 * _S1146.y + cy_5;
    float2  _S1149 = float2 {_S1036.x, _S1036.y};
    float _S1150 = length_0(_S1149);
    float _S1151 = _S1036.z;
    float _S1152 = s_primal_ctx_atan2_0(_S1150, _S1151);
    bool _S1153 = _S1152 < 0.00100000004749745f;
    float _S1154;
    float _S1155;
    float _S1156;
    if(_S1153)
    {
        float _S1157 = 1.0f - _S1152 * _S1152 / 3.0f;
        float _S1158 = _S1151 * _S1151;
        k_2 = _S1157 / _S1151;
        _S1154 = _S1158;
        _S1155 = _S1157;
        _S1156 = 0.0f;
    }
    else
    {
        float _S1159 = _S1150 * _S1150;
        k_2 = _S1152 / _S1150;
        _S1154 = 0.0f;
        _S1155 = 0.0f;
        _S1156 = _S1159;
    }
    float2  _S1160 = make_float2 (k_2);
    float2  _S1161 = _S1149 * make_float2 (k_2);
    float u_60 = _S1161.x;
    float v_60 = _S1161.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S1162 = dist_coeffs_5[int(2)] + r2_60 * dist_coeffs_5[int(3)];
    float _S1163 = dist_coeffs_5[int(1)] + r2_60 * _S1162;
    float _S1164 = dist_coeffs_5[int(0)] + r2_60 * _S1163;
    float radial_28 = 1.0f + r2_60 * _S1164;
    float _S1165 = _S1067 * u_60;
    float _S1166 = 2.0f * u_60;
    float _S1167 = _S1070 * u_60;
    float _S1168 = 2.0f * v_60;
    float2  _S1169 = _S1161 * make_float2 (radial_28) + make_float2 (_S1165 * v_60 + dist_coeffs_5[int(5)] * (r2_60 + _S1166 * u_60) + dist_coeffs_5[int(6)] * r2_60, _S1167 * v_60 + dist_coeffs_5[int(4)] * (r2_60 + _S1168 * v_60) + dist_coeffs_5[int(7)] * r2_60);
    float2  _S1170 = _S1169 + make_float2 (dist_coeffs_5[int(8)] * _S1169.x + dist_coeffs_5[int(9)] * _S1169.y, 0.0f);
    float _S1171 = fx_5 * _S1170.x + cx_5;
    float _S1172 = fy_5 * _S1170.y + cy_5;
    float2  _S1173 = float2 {_S1040.x, _S1040.y};
    float _S1174 = length_0(_S1173);
    float _S1175 = _S1040.z;
    float _S1176 = s_primal_ctx_atan2_0(_S1174, _S1175);
    bool _S1177 = _S1176 < 0.00100000004749745f;
    float _S1178;
    float _S1179;
    float _S1180;
    if(_S1177)
    {
        float _S1181 = 1.0f - _S1176 * _S1176 / 3.0f;
        float _S1182 = _S1175 * _S1175;
        k_2 = _S1181 / _S1175;
        _S1178 = _S1182;
        _S1179 = _S1181;
        _S1180 = 0.0f;
    }
    else
    {
        float _S1183 = _S1174 * _S1174;
        k_2 = _S1176 / _S1174;
        _S1178 = 0.0f;
        _S1179 = 0.0f;
        _S1180 = _S1183;
    }
    float2  _S1184 = make_float2 (k_2);
    float2  _S1185 = _S1173 * make_float2 (k_2);
    float u_61 = _S1185.x;
    float v_61 = _S1185.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S1186 = dist_coeffs_5[int(2)] + r2_61 * dist_coeffs_5[int(3)];
    float _S1187 = dist_coeffs_5[int(1)] + r2_61 * _S1186;
    float _S1188 = dist_coeffs_5[int(0)] + r2_61 * _S1187;
    float radial_29 = 1.0f + r2_61 * _S1188;
    float _S1189 = _S1067 * u_61;
    float _S1190 = 2.0f * u_61;
    float _S1191 = _S1070 * u_61;
    float _S1192 = 2.0f * v_61;
    float2  _S1193 = _S1185 * make_float2 (radial_29) + make_float2 (_S1189 * v_61 + dist_coeffs_5[int(5)] * (r2_61 + _S1190 * u_61) + dist_coeffs_5[int(6)] * r2_61, _S1191 * v_61 + dist_coeffs_5[int(4)] * (r2_61 + _S1192 * v_61) + dist_coeffs_5[int(7)] * r2_61);
    float2  _S1194 = _S1193 + make_float2 (dist_coeffs_5[int(8)] * _S1193.x + dist_coeffs_5[int(9)] * _S1193.y, 0.0f);
    float _S1195 = fx_5 * _S1194.x + cx_5;
    float _S1196 = fy_5 * _S1194.y + cy_5;
    float2  _S1197 = float2 {_S1044.x, _S1044.y};
    float _S1198 = length_0(_S1197);
    float _S1199 = _S1044.z;
    float _S1200 = s_primal_ctx_atan2_0(_S1198, _S1199);
    bool _S1201 = _S1200 < 0.00100000004749745f;
    float _S1202;
    float _S1203;
    float _S1204;
    if(_S1201)
    {
        float _S1205 = 1.0f - _S1200 * _S1200 / 3.0f;
        float _S1206 = _S1199 * _S1199;
        k_2 = _S1205 / _S1199;
        _S1202 = _S1206;
        _S1203 = _S1205;
        _S1204 = 0.0f;
    }
    else
    {
        float _S1207 = _S1198 * _S1198;
        k_2 = _S1200 / _S1198;
        _S1202 = 0.0f;
        _S1203 = 0.0f;
        _S1204 = _S1207;
    }
    float2  _S1208 = make_float2 (k_2);
    float2  _S1209 = _S1197 * make_float2 (k_2);
    float u_62 = _S1209.x;
    float v_62 = _S1209.y;
    float r2_62 = u_62 * u_62 + v_62 * v_62;
    float _S1210 = dist_coeffs_5[int(2)] + r2_62 * dist_coeffs_5[int(3)];
    float _S1211 = dist_coeffs_5[int(1)] + r2_62 * _S1210;
    float _S1212 = dist_coeffs_5[int(0)] + r2_62 * _S1211;
    float radial_30 = 1.0f + r2_62 * _S1212;
    float _S1213 = _S1067 * u_62;
    float _S1214 = 2.0f * u_62;
    float _S1215 = _S1070 * u_62;
    float _S1216 = 2.0f * v_62;
    float2  _S1217 = _S1209 * make_float2 (radial_30) + make_float2 (_S1213 * v_62 + dist_coeffs_5[int(5)] * (r2_62 + _S1214 * u_62) + dist_coeffs_5[int(6)] * r2_62, _S1215 * v_62 + dist_coeffs_5[int(4)] * (r2_62 + _S1216 * v_62) + dist_coeffs_5[int(7)] * r2_62);
    float2  _S1218 = _S1217 + make_float2 (dist_coeffs_5[int(8)] * _S1217.x + dist_coeffs_5[int(9)] * _S1217.y, 0.0f);
    float _S1219 = fx_5 * _S1218.x + cx_5;
    float _S1220 = fy_5 * _S1218.y + cy_5;
    float2  _S1221 = float2 {_S1048.x, _S1048.y};
    float _S1222 = length_0(_S1221);
    float _S1223 = _S1048.z;
    float _S1224 = s_primal_ctx_atan2_0(_S1222, _S1223);
    bool _S1225 = _S1224 < 0.00100000004749745f;
    float _S1226;
    float _S1227;
    float _S1228;
    if(_S1225)
    {
        float _S1229 = 1.0f - _S1224 * _S1224 / 3.0f;
        float _S1230 = _S1223 * _S1223;
        k_2 = _S1229 / _S1223;
        _S1226 = _S1230;
        _S1227 = _S1229;
        _S1228 = 0.0f;
    }
    else
    {
        float _S1231 = _S1222 * _S1222;
        k_2 = _S1224 / _S1222;
        _S1226 = 0.0f;
        _S1227 = 0.0f;
        _S1228 = _S1231;
    }
    float2  _S1232 = make_float2 (k_2);
    float2  _S1233 = _S1221 * make_float2 (k_2);
    float u_63 = _S1233.x;
    float v_63 = _S1233.y;
    float r2_63 = u_63 * u_63 + v_63 * v_63;
    float _S1234 = dist_coeffs_5[int(2)] + r2_63 * dist_coeffs_5[int(3)];
    float _S1235 = dist_coeffs_5[int(1)] + r2_63 * _S1234;
    float _S1236 = dist_coeffs_5[int(0)] + r2_63 * _S1235;
    float radial_31 = 1.0f + r2_63 * _S1236;
    float _S1237 = _S1067 * u_63;
    float _S1238 = 2.0f * u_63;
    float _S1239 = _S1070 * u_63;
    float _S1240 = 2.0f * v_63;
    float2  _S1241 = _S1233 * make_float2 (radial_31) + make_float2 (_S1237 * v_63 + dist_coeffs_5[int(5)] * (r2_63 + _S1238 * u_63) + dist_coeffs_5[int(6)] * r2_63, _S1239 * v_63 + dist_coeffs_5[int(4)] * (r2_63 + _S1240 * v_63) + dist_coeffs_5[int(7)] * r2_63);
    float2  _S1242 = _S1241 + make_float2 (dist_coeffs_5[int(8)] * _S1241.x + dist_coeffs_5[int(9)] * _S1241.y, 0.0f);
    float _S1243 = fx_5 * _S1242.x + cx_5;
    float _S1244 = fy_5 * _S1242.y + cy_5;
    float _S1245 = (F32_max((_S1075), (_S1099)));
    float _S1246 = (F32_min((_S1075), (_S1099)));
    float _S1247 = (F32_max((_S1076), (_S1100)));
    float _S1248 = (F32_min((_S1076), (_S1100)));
    float _S1249 = (F32_max((_S1245), (_S1123)));
    float _S1250 = (F32_min((_S1246), (_S1123)));
    float _S1251 = (F32_max((_S1247), (_S1124)));
    float _S1252 = (F32_min((_S1248), (_S1124)));
    float _S1253 = (F32_max((_S1249), (_S1147)));
    float _S1254 = (F32_min((_S1250), (_S1147)));
    float _S1255 = (F32_max((_S1251), (_S1148)));
    float _S1256 = (F32_min((_S1252), (_S1148)));
    float _S1257 = (F32_max((_S1253), (_S1171)));
    float _S1258 = (F32_min((_S1254), (_S1171)));
    float _S1259 = (F32_max((_S1255), (_S1172)));
    float _S1260 = (F32_min((_S1256), (_S1172)));
    float _S1261 = (F32_max((_S1257), (_S1195)));
    float _S1262 = (F32_min((_S1258), (_S1195)));
    float _S1263 = (F32_max((_S1259), (_S1196)));
    float _S1264 = (F32_min((_S1260), (_S1196)));
    float _S1265 = (F32_max((_S1261), (_S1219)));
    float _S1266 = (F32_min((_S1262), (_S1219)));
    float _S1267 = (F32_max((_S1263), (_S1220)));
    float _S1268 = (F32_min((_S1264), (_S1220)));
    Matrix<float, 3, 3>  _S1269 = transpose_1(R_5);
    float3  _S1270 = mean_c_4 - - s_primal_ctx_mul_0(_S1269, t_5);
    float _S1271 = _S1270.x;
    float _S1272 = _S1270.y;
    float _S1273 = _S1270.z;
    float _S1274 = _S1271 * _S1271 + _S1272 * _S1272 + _S1273 * _S1273;
    float _S1275 = s_primal_ctx_sqrt_0(_S1274);
    float x_14 = _S1271 / _S1275;
    float3  _S1276 = make_float3 (x_14);
    float _S1277 = _S1275 * _S1275;
    float y_10 = _S1272 / _S1275;
    float z_7 = _S1273 / _S1275;
    float3  _S1278 = make_float3 (z_7);
    float _S1279 = - y_10;
    float3  _S1280 = make_float3 (_S1279);
    float z2_5 = z_7 * z_7;
    float fTmp0B_5 = -1.09254848957061768f * z_7;
    float fC1_5 = x_14 * x_14 - y_10 * y_10;
    float _S1281 = 2.0f * x_14;
    float fS1_5 = _S1281 * y_10;
    float pSH6_1 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float3  _S1282 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_14;
    float3  _S1283 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_10;
    float3  _S1284 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S1285 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S1286 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_7;
    float _S1287 = 1.86588168144226074f * z2_5 - 1.11952900886535645f;
    float pSH12_1 = z_7 * _S1287;
    float3  _S1288 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_14;
    float3  _S1289 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_10;
    float3  _S1290 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S1291 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S1292 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_14 * fC1_5 - y_10 * fS1_5);
    float3  _S1293 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_14 * fS1_5 + y_10 * fC1_5);
    float3  _S1294 = make_float3 (pSH9_1);
    float3  _S1295 = make_float3 (0.0f);
    float3  _S1296 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1297;
    (&_S1297)->primal_0 = make_float3 (0.282094806432724f) * sh_coeffs_5[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1279) * sh_coeffs_5[int(1)] + make_float3 (z_7) * sh_coeffs_5[int(2)] - make_float3 (x_14) * sh_coeffs_5[int(3)]) + (make_float3 (pSH4_1) * sh_coeffs_5[int(4)] + make_float3 (pSH5_1) * sh_coeffs_5[int(5)] + make_float3 (pSH6_1) * sh_coeffs_5[int(6)] + make_float3 (pSH7_1) * sh_coeffs_5[int(7)] + make_float3 (pSH8_1) * sh_coeffs_5[int(8)]) + (make_float3 (pSH9_1) * sh_coeffs_5[int(9)] + make_float3 (pSH10_1) * sh_coeffs_5[int(10)] + make_float3 (pSH11_1) * sh_coeffs_5[int(11)] + make_float3 (pSH12_1) * sh_coeffs_5[int(12)] + make_float3 (pSH13_1) * sh_coeffs_5[int(13)] + make_float3 (pSH14_1) * sh_coeffs_5[int(14)] + make_float3 (pSH15_1) * sh_coeffs_5[int(15)]) + make_float3 (0.5f);
    (&_S1297)->differential_0 = _S1296;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1298;
    (&_S1298)->primal_0 = _S1295;
    (&_S1298)->differential_0 = _S1296;
    s_bwd_prop_max_0(&_S1297, &_S1298, v_rgb_1);
    float3  _S1299 = _S1293 * _S1297.differential_0;
    float3  _S1300 = sh_coeffs_5[int(15)] * _S1297.differential_0;
    float3  _S1301 = _S1291 * _S1297.differential_0;
    float3  _S1302 = sh_coeffs_5[int(14)] * _S1297.differential_0;
    float3  _S1303 = _S1289 * _S1297.differential_0;
    float3  _S1304 = sh_coeffs_5[int(13)] * _S1297.differential_0;
    float3  _S1305 = _S1288 * _S1297.differential_0;
    float3  _S1306 = sh_coeffs_5[int(12)] * _S1297.differential_0;
    float3  _S1307 = _S1290 * _S1297.differential_0;
    float3  _S1308 = sh_coeffs_5[int(11)] * _S1297.differential_0;
    float3  _S1309 = _S1292 * _S1297.differential_0;
    float3  _S1310 = sh_coeffs_5[int(10)] * _S1297.differential_0;
    float3  _S1311 = _S1294 * _S1297.differential_0;
    float3  _S1312 = sh_coeffs_5[int(9)] * _S1297.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1312.x + _S1312.y + _S1312.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1300.x + _S1300.y + _S1300.z);
    float _S1313 = _S1310.x + _S1310.y + _S1310.z;
    float _S1314 = _S1302.x + _S1302.y + _S1302.z;
    float _S1315 = _S1308.x + _S1308.y + _S1308.z;
    float _S1316 = _S1304.x + _S1304.y + _S1304.z;
    float _S1317 = _S1306.x + _S1306.y + _S1306.z;
    float _S1318 = - s_diff_fC2_T_1;
    float3  _S1319 = _S1285 * _S1297.differential_0;
    float3  _S1320 = sh_coeffs_5[int(8)] * _S1297.differential_0;
    float3  _S1321 = _S1283 * _S1297.differential_0;
    float3  _S1322 = sh_coeffs_5[int(7)] * _S1297.differential_0;
    float3  _S1323 = _S1282 * _S1297.differential_0;
    float3  _S1324 = sh_coeffs_5[int(6)] * _S1297.differential_0;
    float3  _S1325 = _S1284 * _S1297.differential_0;
    float3  _S1326 = sh_coeffs_5[int(5)] * _S1297.differential_0;
    float3  _S1327 = _S1286 * _S1297.differential_0;
    float3  _S1328 = sh_coeffs_5[int(4)] * _S1297.differential_0;
    float _S1329 = _S1326.x + _S1326.y + _S1326.z;
    float _S1330 = _S1322.x + _S1322.y + _S1322.z;
    float _S1331 = fTmp1B_5 * _S1313 + x_14 * s_diff_fS2_T_1 + y_10 * _S1318 + 0.54627424478530884f * (_S1328.x + _S1328.y + _S1328.z);
    float _S1332 = fTmp1B_5 * _S1314 + y_10 * s_diff_fS2_T_1 + x_14 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1320.x + _S1320.y + _S1320.z);
    float _S1333 = y_10 * - _S1332;
    float _S1334 = x_14 * _S1332;
    float _S1335 = z_7 * (1.86588168144226074f * (z_7 * _S1317) + -2.28522896766662598f * (y_10 * _S1315 + x_14 * _S1316) + 0.94617468118667603f * (_S1324.x + _S1324.y + _S1324.z));
    float3  _S1336 = make_float3 (0.48860251903533936f) * _S1297.differential_0;
    float3  _S1337 = - _S1336;
    float3  _S1338 = _S1276 * _S1337;
    float3  _S1339 = sh_coeffs_5[int(3)] * _S1337;
    float3  _S1340 = _S1278 * _S1336;
    float3  _S1341 = sh_coeffs_5[int(2)] * _S1336;
    float3  _S1342 = _S1280 * _S1336;
    float3  _S1343 = sh_coeffs_5[int(1)] * _S1336;
    float _S1344 = (_S1287 * _S1317 + 1.44530570507049561f * (fS1_5 * _S1313 + fC1_5 * _S1314) + -1.09254848957061768f * (y_10 * _S1329 + x_14 * _S1330) + _S1335 + _S1335 + _S1341.x + _S1341.y + _S1341.z) / _S1277;
    float _S1345 = _S1275 * _S1344;
    float _S1346 = (fTmp0C_5 * _S1315 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S1318 + fTmp0B_5 * _S1329 + _S1281 * _S1331 + _S1333 + _S1333 + - (_S1343.x + _S1343.y + _S1343.z)) / _S1277;
    float _S1347 = _S1275 * _S1346;
    float _S1348 = (fTmp0C_5 * _S1316 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S1330 + 2.0f * (y_10 * _S1331) + _S1334 + _S1334 + _S1339.x + _S1339.y + _S1339.z) / _S1277;
    float _S1349 = _S1275 * _S1348;
    float _S1350 = _S1273 * - _S1344 + _S1272 * - _S1346 + _S1271 * - _S1348;
    DiffPair_float_0 _S1351;
    (&_S1351)->primal_0 = _S1274;
    (&_S1351)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1351, _S1350);
    float _S1352 = _S1273 * _S1351.differential_0;
    float _S1353 = _S1272 * _S1351.differential_0;
    float _S1354 = _S1271 * _S1351.differential_0;
    float3  _S1355 = make_float3 (0.282094806432724f) * _S1297.differential_0;
    float3  _S1356 = make_float3 (_S1349 + _S1354 + _S1354, _S1347 + _S1353 + _S1353, _S1345 + _S1352 + _S1352);
    float3  _S1357 = - - _S1356;
    Matrix<float, 3, 3>  _S1358 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1359;
    (&_S1359)->primal_0 = _S1269;
    (&_S1359)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1360;
    (&_S1360)->primal_0 = t_5;
    (&_S1360)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1359, &_S1360, _S1357);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1361 = _S1360;
    Matrix<float, 3, 3>  _S1362 = transpose_1(_S1359.differential_0);
    DiffPair_float_0 _S1363;
    (&_S1363)->primal_0 = _S1268;
    (&_S1363)->differential_0 = 0.0f;
    DiffPair_float_0 _S1364;
    (&_S1364)->primal_0 = _S1244;
    (&_S1364)->differential_0 = 0.0f;
    _d_min_0(&_S1363, &_S1364, 0.0f);
    DiffPair_float_0 _S1365;
    (&_S1365)->primal_0 = _S1267;
    (&_S1365)->differential_0 = 0.0f;
    DiffPair_float_0 _S1366;
    (&_S1366)->primal_0 = _S1244;
    (&_S1366)->differential_0 = 0.0f;
    _d_max_0(&_S1365, &_S1366, 0.0f);
    DiffPair_float_0 _S1367;
    (&_S1367)->primal_0 = _S1266;
    (&_S1367)->differential_0 = 0.0f;
    DiffPair_float_0 _S1368;
    (&_S1368)->primal_0 = _S1243;
    (&_S1368)->differential_0 = 0.0f;
    _d_min_0(&_S1367, &_S1368, 0.0f);
    DiffPair_float_0 _S1369;
    (&_S1369)->primal_0 = _S1265;
    (&_S1369)->differential_0 = 0.0f;
    DiffPair_float_0 _S1370;
    (&_S1370)->primal_0 = _S1243;
    (&_S1370)->differential_0 = 0.0f;
    _d_max_0(&_S1369, &_S1370, 0.0f);
    DiffPair_float_0 _S1371;
    (&_S1371)->primal_0 = _S1264;
    (&_S1371)->differential_0 = 0.0f;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1220;
    (&_S1372)->differential_0 = 0.0f;
    _d_min_0(&_S1371, &_S1372, _S1363.differential_0);
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S1263;
    (&_S1373)->differential_0 = 0.0f;
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = _S1220;
    (&_S1374)->differential_0 = 0.0f;
    _d_max_0(&_S1373, &_S1374, _S1365.differential_0);
    DiffPair_float_0 _S1375;
    (&_S1375)->primal_0 = _S1262;
    (&_S1375)->differential_0 = 0.0f;
    DiffPair_float_0 _S1376;
    (&_S1376)->primal_0 = _S1219;
    (&_S1376)->differential_0 = 0.0f;
    _d_min_0(&_S1375, &_S1376, _S1367.differential_0);
    DiffPair_float_0 _S1377;
    (&_S1377)->primal_0 = _S1261;
    (&_S1377)->differential_0 = 0.0f;
    DiffPair_float_0 _S1378;
    (&_S1378)->primal_0 = _S1219;
    (&_S1378)->differential_0 = 0.0f;
    _d_max_0(&_S1377, &_S1378, _S1369.differential_0);
    DiffPair_float_0 _S1379;
    (&_S1379)->primal_0 = _S1260;
    (&_S1379)->differential_0 = 0.0f;
    DiffPair_float_0 _S1380;
    (&_S1380)->primal_0 = _S1196;
    (&_S1380)->differential_0 = 0.0f;
    _d_min_0(&_S1379, &_S1380, _S1371.differential_0);
    DiffPair_float_0 _S1381;
    (&_S1381)->primal_0 = _S1259;
    (&_S1381)->differential_0 = 0.0f;
    DiffPair_float_0 _S1382;
    (&_S1382)->primal_0 = _S1196;
    (&_S1382)->differential_0 = 0.0f;
    _d_max_0(&_S1381, &_S1382, _S1373.differential_0);
    DiffPair_float_0 _S1383;
    (&_S1383)->primal_0 = _S1258;
    (&_S1383)->differential_0 = 0.0f;
    DiffPair_float_0 _S1384;
    (&_S1384)->primal_0 = _S1195;
    (&_S1384)->differential_0 = 0.0f;
    _d_min_0(&_S1383, &_S1384, _S1375.differential_0);
    DiffPair_float_0 _S1385;
    (&_S1385)->primal_0 = _S1257;
    (&_S1385)->differential_0 = 0.0f;
    DiffPair_float_0 _S1386;
    (&_S1386)->primal_0 = _S1195;
    (&_S1386)->differential_0 = 0.0f;
    _d_max_0(&_S1385, &_S1386, _S1377.differential_0);
    DiffPair_float_0 _S1387;
    (&_S1387)->primal_0 = _S1256;
    (&_S1387)->differential_0 = 0.0f;
    DiffPair_float_0 _S1388;
    (&_S1388)->primal_0 = _S1172;
    (&_S1388)->differential_0 = 0.0f;
    _d_min_0(&_S1387, &_S1388, _S1379.differential_0);
    DiffPair_float_0 _S1389;
    (&_S1389)->primal_0 = _S1255;
    (&_S1389)->differential_0 = 0.0f;
    DiffPair_float_0 _S1390;
    (&_S1390)->primal_0 = _S1172;
    (&_S1390)->differential_0 = 0.0f;
    _d_max_0(&_S1389, &_S1390, _S1381.differential_0);
    DiffPair_float_0 _S1391;
    (&_S1391)->primal_0 = _S1254;
    (&_S1391)->differential_0 = 0.0f;
    DiffPair_float_0 _S1392;
    (&_S1392)->primal_0 = _S1171;
    (&_S1392)->differential_0 = 0.0f;
    _d_min_0(&_S1391, &_S1392, _S1383.differential_0);
    DiffPair_float_0 _S1393;
    (&_S1393)->primal_0 = _S1253;
    (&_S1393)->differential_0 = 0.0f;
    DiffPair_float_0 _S1394;
    (&_S1394)->primal_0 = _S1171;
    (&_S1394)->differential_0 = 0.0f;
    _d_max_0(&_S1393, &_S1394, _S1385.differential_0);
    DiffPair_float_0 _S1395;
    (&_S1395)->primal_0 = _S1252;
    (&_S1395)->differential_0 = 0.0f;
    DiffPair_float_0 _S1396;
    (&_S1396)->primal_0 = _S1148;
    (&_S1396)->differential_0 = 0.0f;
    _d_min_0(&_S1395, &_S1396, _S1387.differential_0);
    DiffPair_float_0 _S1397;
    (&_S1397)->primal_0 = _S1251;
    (&_S1397)->differential_0 = 0.0f;
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1148;
    (&_S1398)->differential_0 = 0.0f;
    _d_max_0(&_S1397, &_S1398, _S1389.differential_0);
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1250;
    (&_S1399)->differential_0 = 0.0f;
    DiffPair_float_0 _S1400;
    (&_S1400)->primal_0 = _S1147;
    (&_S1400)->differential_0 = 0.0f;
    _d_min_0(&_S1399, &_S1400, _S1391.differential_0);
    DiffPair_float_0 _S1401;
    (&_S1401)->primal_0 = _S1249;
    (&_S1401)->differential_0 = 0.0f;
    DiffPair_float_0 _S1402;
    (&_S1402)->primal_0 = _S1147;
    (&_S1402)->differential_0 = 0.0f;
    _d_max_0(&_S1401, &_S1402, _S1393.differential_0);
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = _S1248;
    (&_S1403)->differential_0 = 0.0f;
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = _S1124;
    (&_S1404)->differential_0 = 0.0f;
    _d_min_0(&_S1403, &_S1404, _S1395.differential_0);
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = _S1247;
    (&_S1405)->differential_0 = 0.0f;
    DiffPair_float_0 _S1406;
    (&_S1406)->primal_0 = _S1124;
    (&_S1406)->differential_0 = 0.0f;
    _d_max_0(&_S1405, &_S1406, _S1397.differential_0);
    DiffPair_float_0 _S1407;
    (&_S1407)->primal_0 = _S1246;
    (&_S1407)->differential_0 = 0.0f;
    DiffPair_float_0 _S1408;
    (&_S1408)->primal_0 = _S1123;
    (&_S1408)->differential_0 = 0.0f;
    _d_min_0(&_S1407, &_S1408, _S1399.differential_0);
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1245;
    (&_S1409)->differential_0 = 0.0f;
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = _S1123;
    (&_S1410)->differential_0 = 0.0f;
    _d_max_0(&_S1409, &_S1410, _S1401.differential_0);
    DiffPair_float_0 _S1411;
    (&_S1411)->primal_0 = _S1076;
    (&_S1411)->differential_0 = 0.0f;
    DiffPair_float_0 _S1412;
    (&_S1412)->primal_0 = _S1100;
    (&_S1412)->differential_0 = 0.0f;
    _d_min_0(&_S1411, &_S1412, _S1403.differential_0);
    DiffPair_float_0 _S1413;
    (&_S1413)->primal_0 = _S1076;
    (&_S1413)->differential_0 = 0.0f;
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = _S1100;
    (&_S1414)->differential_0 = 0.0f;
    _d_max_0(&_S1413, &_S1414, _S1405.differential_0);
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1075;
    (&_S1415)->differential_0 = 0.0f;
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1099;
    (&_S1416)->differential_0 = 0.0f;
    _d_min_0(&_S1415, &_S1416, _S1407.differential_0);
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1075;
    (&_S1417)->differential_0 = 0.0f;
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = _S1099;
    (&_S1418)->differential_0 = 0.0f;
    _d_max_0(&_S1417, &_S1418, _S1409.differential_0);
    float _S1419 = fx_5 * (_S1368.differential_0 + _S1370.differential_0);
    float2  _S1420 = make_float2 (_S1419, fy_5 * (_S1364.differential_0 + _S1366.differential_0)) + make_float2 (dist_coeffs_5[int(8)] * _S1419, dist_coeffs_5[int(9)] * _S1419);
    float2  _S1421 = _S1233 * _S1420;
    float _S1422 = dist_coeffs_5[int(4)] * _S1420.y;
    float _S1423 = dist_coeffs_5[int(5)] * _S1420.x;
    float _S1424 = _S1421.x + _S1421.y;
    float _S1425 = r2_63 * _S1424;
    float _S1426 = r2_63 * _S1425;
    float _S1427 = dist_coeffs_5[int(7)] * _S1420.y + _S1422 + dist_coeffs_5[int(6)] * _S1420.x + _S1423 + _S1236 * _S1424 + _S1235 * _S1425 + _S1234 * _S1426 + dist_coeffs_5[int(3)] * (r2_63 * _S1426);
    float _S1428 = v_63 * _S1427;
    float _S1429 = u_63 * _S1427;
    float2  _S1430 = make_float2 (radial_31) * _S1420 + make_float2 (_S1070 * (v_63 * _S1420.y) + _S1238 * _S1423 + 2.0f * (u_63 * _S1423) + _S1067 * (v_63 * _S1420.x) + _S1429 + _S1429, _S1240 * _S1422 + 2.0f * (v_63 * _S1422) + _S1239 * _S1420.y + _S1237 * _S1420.x + _S1428 + _S1428);
    FixedArray<float3 , 16>  _S1431;
    _S1431[int(0)] = _S1296;
    _S1431[int(1)] = _S1296;
    _S1431[int(2)] = _S1296;
    _S1431[int(3)] = _S1296;
    _S1431[int(4)] = _S1296;
    _S1431[int(5)] = _S1296;
    _S1431[int(6)] = _S1296;
    _S1431[int(7)] = _S1296;
    _S1431[int(8)] = _S1296;
    _S1431[int(9)] = _S1296;
    _S1431[int(10)] = _S1296;
    _S1431[int(11)] = _S1296;
    _S1431[int(12)] = _S1296;
    _S1431[int(13)] = _S1296;
    _S1431[int(14)] = _S1296;
    _S1431[int(15)] = _S1296;
    _S1431[int(7)] = _S1321;
    _S1431[int(0)] = _S1355;
    _S1431[int(1)] = _S1342;
    _S1431[int(2)] = _S1340;
    _S1431[int(3)] = _S1338;
    _S1431[int(4)] = _S1327;
    _S1431[int(5)] = _S1325;
    _S1431[int(6)] = _S1323;
    _S1431[int(15)] = _S1299;
    _S1431[int(8)] = _S1319;
    _S1431[int(9)] = _S1311;
    _S1431[int(10)] = _S1309;
    _S1431[int(11)] = _S1307;
    _S1431[int(12)] = _S1305;
    _S1431[int(13)] = _S1303;
    _S1431[int(14)] = _S1301;
    float3  _S1432 = _S1431[int(0)];
    float3  _S1433 = _S1431[int(1)];
    float3  _S1434 = _S1431[int(2)];
    float3  _S1435 = _S1431[int(3)];
    float3  _S1436 = _S1431[int(4)];
    float3  _S1437 = _S1431[int(5)];
    float3  _S1438 = _S1431[int(6)];
    float3  _S1439 = _S1431[int(7)];
    float3  _S1440 = _S1431[int(8)];
    float3  _S1441 = _S1431[int(9)];
    float3  _S1442 = _S1431[int(10)];
    float3  _S1443 = _S1431[int(11)];
    float3  _S1444 = _S1431[int(12)];
    float3  _S1445 = _S1431[int(13)];
    float3  _S1446 = _S1431[int(14)];
    float3  _S1447 = _S1431[int(15)];
    float _S1448 = _S1416.differential_0 + _S1418.differential_0;
    float _S1449 = _S1411.differential_0 + _S1413.differential_0;
    float _S1450 = _S1412.differential_0 + _S1414.differential_0;
    float _S1451 = _S1408.differential_0 + _S1410.differential_0;
    float _S1452 = _S1415.differential_0 + _S1417.differential_0;
    float _S1453 = _S1372.differential_0 + _S1374.differential_0;
    float _S1454 = _S1376.differential_0 + _S1378.differential_0;
    float _S1455 = _S1380.differential_0 + _S1382.differential_0;
    float _S1456 = _S1384.differential_0 + _S1386.differential_0;
    float _S1457 = _S1388.differential_0 + _S1390.differential_0;
    float _S1458 = _S1392.differential_0 + _S1394.differential_0;
    float _S1459 = _S1396.differential_0 + _S1398.differential_0;
    float _S1460 = _S1400.differential_0 + _S1402.differential_0;
    float _S1461 = _S1404.differential_0 + _S1406.differential_0;
    float2  _S1462 = _S1221 * _S1430;
    float2  _S1463 = _S1232 * _S1430;
    float _S1464 = _S1462.x + _S1462.y;
    if(_S1225)
    {
        float _S1465 = _S1464 / _S1226;
        float _S1466 = _S1227 * - _S1465;
        float _S1467 = _S1224 * (0.3333333432674408f * - (_S1223 * _S1465));
        k_2 = _S1467 + _S1467;
        _S1226 = _S1466;
        _S1227 = 0.0f;
    }
    else
    {
        float _S1468 = _S1464 / _S1228;
        float _S1469 = _S1224 * - _S1468;
        k_2 = _S1222 * _S1468;
        _S1226 = 0.0f;
        _S1227 = _S1469;
    }
    DiffPair_float_0 _S1470;
    (&_S1470)->primal_0 = _S1222;
    (&_S1470)->differential_0 = 0.0f;
    DiffPair_float_0 _S1471;
    (&_S1471)->primal_0 = _S1223;
    (&_S1471)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1470, &_S1471, k_2);
    float _S1472 = _S1471.differential_0 + _S1226;
    float _S1473 = _S1470.differential_0 + _S1227;
    float2  _S1474 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1475;
    (&_S1475)->primal_0 = _S1221;
    (&_S1475)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1475, _S1473);
    float2  _S1476 = _S1475.differential_0 + _S1463;
    float _S1477 = fx_5 * _S1454;
    float2  _S1478 = make_float2 (_S1477, fy_5 * _S1453) + make_float2 (dist_coeffs_5[int(8)] * _S1477, dist_coeffs_5[int(9)] * _S1477);
    float2  _S1479 = _S1209 * _S1478;
    float _S1480 = dist_coeffs_5[int(4)] * _S1478.y;
    float _S1481 = dist_coeffs_5[int(5)] * _S1478.x;
    float _S1482 = _S1479.x + _S1479.y;
    float _S1483 = r2_62 * _S1482;
    float _S1484 = r2_62 * _S1483;
    float _S1485 = dist_coeffs_5[int(7)] * _S1478.y + _S1480 + dist_coeffs_5[int(6)] * _S1478.x + _S1481 + _S1212 * _S1482 + _S1211 * _S1483 + _S1210 * _S1484 + dist_coeffs_5[int(3)] * (r2_62 * _S1484);
    float _S1486 = v_62 * _S1485;
    float _S1487 = u_62 * _S1485;
    float2  _S1488 = make_float2 (radial_30) * _S1478 + make_float2 (_S1070 * (v_62 * _S1478.y) + _S1214 * _S1481 + 2.0f * (u_62 * _S1481) + _S1067 * (v_62 * _S1478.x) + _S1487 + _S1487, _S1216 * _S1480 + 2.0f * (v_62 * _S1480) + _S1215 * _S1478.y + _S1213 * _S1478.x + _S1486 + _S1486);
    float3  _S1489 = make_float3 (_S1476.x, _S1476.y, _S1472);
    float2  _S1490 = _S1197 * _S1488;
    float2  _S1491 = _S1208 * _S1488;
    float _S1492 = _S1490.x + _S1490.y;
    if(_S1201)
    {
        float _S1493 = _S1492 / _S1202;
        float _S1494 = _S1203 * - _S1493;
        float _S1495 = _S1200 * (0.3333333432674408f * - (_S1199 * _S1493));
        k_2 = _S1495 + _S1495;
        _S1202 = _S1494;
        _S1203 = 0.0f;
    }
    else
    {
        float _S1496 = _S1492 / _S1204;
        float _S1497 = _S1200 * - _S1496;
        k_2 = _S1198 * _S1496;
        _S1202 = 0.0f;
        _S1203 = _S1497;
    }
    DiffPair_float_0 _S1498;
    (&_S1498)->primal_0 = _S1198;
    (&_S1498)->differential_0 = 0.0f;
    DiffPair_float_0 _S1499;
    (&_S1499)->primal_0 = _S1199;
    (&_S1499)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1498, &_S1499, k_2);
    float _S1500 = _S1499.differential_0 + _S1202;
    float _S1501 = _S1498.differential_0 + _S1203;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1502;
    (&_S1502)->primal_0 = _S1197;
    (&_S1502)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1502, _S1501);
    float2  _S1503 = _S1502.differential_0 + _S1491;
    float _S1504 = fx_5 * _S1456;
    float2  _S1505 = make_float2 (_S1504, fy_5 * _S1455) + make_float2 (dist_coeffs_5[int(8)] * _S1504, dist_coeffs_5[int(9)] * _S1504);
    float2  _S1506 = _S1185 * _S1505;
    float _S1507 = dist_coeffs_5[int(4)] * _S1505.y;
    float _S1508 = dist_coeffs_5[int(5)] * _S1505.x;
    float _S1509 = _S1506.x + _S1506.y;
    float _S1510 = r2_61 * _S1509;
    float _S1511 = r2_61 * _S1510;
    float _S1512 = dist_coeffs_5[int(7)] * _S1505.y + _S1507 + dist_coeffs_5[int(6)] * _S1505.x + _S1508 + _S1188 * _S1509 + _S1187 * _S1510 + _S1186 * _S1511 + dist_coeffs_5[int(3)] * (r2_61 * _S1511);
    float _S1513 = v_61 * _S1512;
    float _S1514 = u_61 * _S1512;
    float2  _S1515 = make_float2 (radial_29) * _S1505 + make_float2 (_S1070 * (v_61 * _S1505.y) + _S1190 * _S1508 + 2.0f * (u_61 * _S1508) + _S1067 * (v_61 * _S1505.x) + _S1514 + _S1514, _S1192 * _S1507 + 2.0f * (v_61 * _S1507) + _S1191 * _S1505.y + _S1189 * _S1505.x + _S1513 + _S1513);
    float3  _S1516 = make_float3 (_S1503.x, _S1503.y, _S1500);
    float2  _S1517 = _S1173 * _S1515;
    float2  _S1518 = _S1184 * _S1515;
    float _S1519 = _S1517.x + _S1517.y;
    if(_S1177)
    {
        float _S1520 = _S1519 / _S1178;
        float _S1521 = _S1179 * - _S1520;
        float _S1522 = _S1176 * (0.3333333432674408f * - (_S1175 * _S1520));
        k_2 = _S1522 + _S1522;
        _S1178 = _S1521;
        _S1179 = 0.0f;
    }
    else
    {
        float _S1523 = _S1519 / _S1180;
        float _S1524 = _S1176 * - _S1523;
        k_2 = _S1174 * _S1523;
        _S1178 = 0.0f;
        _S1179 = _S1524;
    }
    DiffPair_float_0 _S1525;
    (&_S1525)->primal_0 = _S1174;
    (&_S1525)->differential_0 = 0.0f;
    DiffPair_float_0 _S1526;
    (&_S1526)->primal_0 = _S1175;
    (&_S1526)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1525, &_S1526, k_2);
    float _S1527 = _S1526.differential_0 + _S1178;
    float _S1528 = _S1525.differential_0 + _S1179;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1529;
    (&_S1529)->primal_0 = _S1173;
    (&_S1529)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1529, _S1528);
    float2  _S1530 = _S1529.differential_0 + _S1518;
    float _S1531 = fx_5 * _S1458;
    float2  _S1532 = make_float2 (_S1531, fy_5 * _S1457) + make_float2 (dist_coeffs_5[int(8)] * _S1531, dist_coeffs_5[int(9)] * _S1531);
    float2  _S1533 = _S1161 * _S1532;
    float _S1534 = dist_coeffs_5[int(4)] * _S1532.y;
    float _S1535 = dist_coeffs_5[int(5)] * _S1532.x;
    float _S1536 = _S1533.x + _S1533.y;
    float _S1537 = r2_60 * _S1536;
    float _S1538 = r2_60 * _S1537;
    float _S1539 = dist_coeffs_5[int(7)] * _S1532.y + _S1534 + dist_coeffs_5[int(6)] * _S1532.x + _S1535 + _S1164 * _S1536 + _S1163 * _S1537 + _S1162 * _S1538 + dist_coeffs_5[int(3)] * (r2_60 * _S1538);
    float _S1540 = v_60 * _S1539;
    float _S1541 = u_60 * _S1539;
    float2  _S1542 = make_float2 (radial_28) * _S1532 + make_float2 (_S1070 * (v_60 * _S1532.y) + _S1166 * _S1535 + 2.0f * (u_60 * _S1535) + _S1067 * (v_60 * _S1532.x) + _S1541 + _S1541, _S1168 * _S1534 + 2.0f * (v_60 * _S1534) + _S1167 * _S1532.y + _S1165 * _S1532.x + _S1540 + _S1540);
    float3  _S1543 = make_float3 (_S1530.x, _S1530.y, _S1527);
    float2  _S1544 = _S1149 * _S1542;
    float2  _S1545 = _S1160 * _S1542;
    float _S1546 = _S1544.x + _S1544.y;
    if(_S1153)
    {
        float _S1547 = _S1546 / _S1154;
        float _S1548 = _S1155 * - _S1547;
        float _S1549 = _S1152 * (0.3333333432674408f * - (_S1151 * _S1547));
        k_2 = _S1549 + _S1549;
        _S1154 = _S1548;
        _S1155 = 0.0f;
    }
    else
    {
        float _S1550 = _S1546 / _S1156;
        float _S1551 = _S1152 * - _S1550;
        k_2 = _S1150 * _S1550;
        _S1154 = 0.0f;
        _S1155 = _S1551;
    }
    DiffPair_float_0 _S1552;
    (&_S1552)->primal_0 = _S1150;
    (&_S1552)->differential_0 = 0.0f;
    DiffPair_float_0 _S1553;
    (&_S1553)->primal_0 = _S1151;
    (&_S1553)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1552, &_S1553, k_2);
    float _S1554 = _S1553.differential_0 + _S1154;
    float _S1555 = _S1552.differential_0 + _S1155;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1556;
    (&_S1556)->primal_0 = _S1149;
    (&_S1556)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1556, _S1555);
    float2  _S1557 = _S1556.differential_0 + _S1545;
    float _S1558 = fx_5 * _S1460;
    float2  _S1559 = make_float2 (_S1558, fy_5 * _S1459) + make_float2 (dist_coeffs_5[int(8)] * _S1558, dist_coeffs_5[int(9)] * _S1558);
    float2  _S1560 = _S1137 * _S1559;
    float _S1561 = dist_coeffs_5[int(4)] * _S1559.y;
    float _S1562 = dist_coeffs_5[int(5)] * _S1559.x;
    float _S1563 = _S1560.x + _S1560.y;
    float _S1564 = r2_59 * _S1563;
    float _S1565 = r2_59 * _S1564;
    float _S1566 = dist_coeffs_5[int(7)] * _S1559.y + _S1561 + dist_coeffs_5[int(6)] * _S1559.x + _S1562 + _S1140 * _S1563 + _S1139 * _S1564 + _S1138 * _S1565 + dist_coeffs_5[int(3)] * (r2_59 * _S1565);
    float _S1567 = v_59 * _S1566;
    float _S1568 = u_59 * _S1566;
    float2  _S1569 = make_float2 (radial_27) * _S1559 + make_float2 (_S1070 * (v_59 * _S1559.y) + _S1142 * _S1562 + 2.0f * (u_59 * _S1562) + _S1067 * (v_59 * _S1559.x) + _S1568 + _S1568, _S1144 * _S1561 + 2.0f * (v_59 * _S1561) + _S1143 * _S1559.y + _S1141 * _S1559.x + _S1567 + _S1567);
    float3  _S1570 = make_float3 (_S1557.x, _S1557.y, _S1554);
    float2  _S1571 = _S1125 * _S1569;
    float2  _S1572 = _S1136 * _S1569;
    float _S1573 = _S1571.x + _S1571.y;
    if(_S1129)
    {
        float _S1574 = _S1573 / _S1130;
        float _S1575 = _S1131 * - _S1574;
        float _S1576 = _S1128 * (0.3333333432674408f * - (_S1127 * _S1574));
        k_2 = _S1576 + _S1576;
        _S1130 = _S1575;
        _S1131 = 0.0f;
    }
    else
    {
        float _S1577 = _S1573 / _S1132;
        float _S1578 = _S1128 * - _S1577;
        k_2 = _S1126 * _S1577;
        _S1130 = 0.0f;
        _S1131 = _S1578;
    }
    DiffPair_float_0 _S1579;
    (&_S1579)->primal_0 = _S1126;
    (&_S1579)->differential_0 = 0.0f;
    DiffPair_float_0 _S1580;
    (&_S1580)->primal_0 = _S1127;
    (&_S1580)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1579, &_S1580, k_2);
    float _S1581 = _S1580.differential_0 + _S1130;
    float _S1582 = _S1579.differential_0 + _S1131;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1583;
    (&_S1583)->primal_0 = _S1125;
    (&_S1583)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1583, _S1582);
    float2  _S1584 = _S1583.differential_0 + _S1572;
    float _S1585 = fx_5 * _S1451;
    float2  _S1586 = make_float2 (_S1585, fy_5 * _S1461) + make_float2 (dist_coeffs_5[int(8)] * _S1585, dist_coeffs_5[int(9)] * _S1585);
    float2  _S1587 = _S1113 * _S1586;
    float _S1588 = dist_coeffs_5[int(4)] * _S1586.y;
    float _S1589 = dist_coeffs_5[int(5)] * _S1586.x;
    float _S1590 = _S1587.x + _S1587.y;
    float _S1591 = r2_58 * _S1590;
    float _S1592 = r2_58 * _S1591;
    float _S1593 = dist_coeffs_5[int(7)] * _S1586.y + _S1588 + dist_coeffs_5[int(6)] * _S1586.x + _S1589 + _S1116 * _S1590 + _S1115 * _S1591 + _S1114 * _S1592 + dist_coeffs_5[int(3)] * (r2_58 * _S1592);
    float _S1594 = v_58 * _S1593;
    float _S1595 = u_58 * _S1593;
    float2  _S1596 = make_float2 (radial_26) * _S1586 + make_float2 (_S1070 * (v_58 * _S1586.y) + _S1118 * _S1589 + 2.0f * (u_58 * _S1589) + _S1067 * (v_58 * _S1586.x) + _S1595 + _S1595, _S1120 * _S1588 + 2.0f * (v_58 * _S1588) + _S1119 * _S1586.y + _S1117 * _S1586.x + _S1594 + _S1594);
    float3  _S1597 = make_float3 (_S1584.x, _S1584.y, _S1581);
    float2  _S1598 = _S1101 * _S1596;
    float2  _S1599 = _S1112 * _S1596;
    float _S1600 = _S1598.x + _S1598.y;
    if(_S1105)
    {
        float _S1601 = _S1600 / _S1106;
        float _S1602 = _S1107 * - _S1601;
        float _S1603 = _S1104 * (0.3333333432674408f * - (_S1103 * _S1601));
        k_2 = _S1603 + _S1603;
        _S1106 = _S1602;
        _S1107 = 0.0f;
    }
    else
    {
        float _S1604 = _S1600 / _S1108;
        float _S1605 = _S1104 * - _S1604;
        k_2 = _S1102 * _S1604;
        _S1106 = 0.0f;
        _S1107 = _S1605;
    }
    DiffPair_float_0 _S1606;
    (&_S1606)->primal_0 = _S1102;
    (&_S1606)->differential_0 = 0.0f;
    DiffPair_float_0 _S1607;
    (&_S1607)->primal_0 = _S1103;
    (&_S1607)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1606, &_S1607, k_2);
    float _S1608 = _S1607.differential_0 + _S1106;
    float _S1609 = _S1606.differential_0 + _S1107;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1610;
    (&_S1610)->primal_0 = _S1101;
    (&_S1610)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1610, _S1609);
    float2  _S1611 = _S1610.differential_0 + _S1599;
    float _S1612 = fx_5 * _S1448;
    float2  _S1613 = make_float2 (_S1612, fy_5 * _S1450) + make_float2 (dist_coeffs_5[int(8)] * _S1612, dist_coeffs_5[int(9)] * _S1612);
    float2  _S1614 = _S1089 * _S1613;
    float _S1615 = dist_coeffs_5[int(4)] * _S1613.y;
    float _S1616 = dist_coeffs_5[int(5)] * _S1613.x;
    float _S1617 = _S1614.x + _S1614.y;
    float _S1618 = r2_57 * _S1617;
    float _S1619 = r2_57 * _S1618;
    float _S1620 = dist_coeffs_5[int(7)] * _S1613.y + _S1615 + dist_coeffs_5[int(6)] * _S1613.x + _S1616 + _S1092 * _S1617 + _S1091 * _S1618 + _S1090 * _S1619 + dist_coeffs_5[int(3)] * (r2_57 * _S1619);
    float _S1621 = v_57 * _S1620;
    float _S1622 = u_57 * _S1620;
    float2  _S1623 = make_float2 (radial_25) * _S1613 + make_float2 (_S1070 * (v_57 * _S1613.y) + _S1094 * _S1616 + 2.0f * (u_57 * _S1616) + _S1067 * (v_57 * _S1613.x) + _S1622 + _S1622, _S1096 * _S1615 + 2.0f * (v_57 * _S1615) + _S1095 * _S1613.y + _S1093 * _S1613.x + _S1621 + _S1621);
    float3  _S1624 = make_float3 (_S1611.x, _S1611.y, _S1608);
    float2  _S1625 = _S1077 * _S1623;
    float2  _S1626 = _S1088 * _S1623;
    float _S1627 = _S1625.x + _S1625.y;
    if(_S1081)
    {
        float _S1628 = _S1627 / _S1082;
        float _S1629 = _S1083 * - _S1628;
        float _S1630 = _S1080 * (0.3333333432674408f * - (_S1079 * _S1628));
        k_2 = _S1630 + _S1630;
        _S1082 = _S1629;
        _S1083 = 0.0f;
    }
    else
    {
        float _S1631 = _S1627 / _S1084;
        float _S1632 = _S1080 * - _S1631;
        k_2 = _S1078 * _S1631;
        _S1082 = 0.0f;
        _S1083 = _S1632;
    }
    DiffPair_float_0 _S1633;
    (&_S1633)->primal_0 = _S1078;
    (&_S1633)->differential_0 = 0.0f;
    DiffPair_float_0 _S1634;
    (&_S1634)->primal_0 = _S1079;
    (&_S1634)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1633, &_S1634, k_2);
    float _S1635 = _S1634.differential_0 + _S1082;
    float _S1636 = _S1633.differential_0 + _S1083;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1637;
    (&_S1637)->primal_0 = _S1077;
    (&_S1637)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1637, _S1636);
    float2  _S1638 = _S1637.differential_0 + _S1626;
    float _S1639 = fx_5 * _S1452;
    float2  _S1640 = make_float2 (_S1639, fy_5 * _S1449) + make_float2 (dist_coeffs_5[int(8)] * _S1639, dist_coeffs_5[int(9)] * _S1639);
    float2  _S1641 = _S1063 * _S1640;
    float _S1642 = dist_coeffs_5[int(4)] * _S1640.y;
    float _S1643 = dist_coeffs_5[int(5)] * _S1640.x;
    float _S1644 = _S1641.x + _S1641.y;
    float _S1645 = r2_56 * _S1644;
    float _S1646 = r2_56 * _S1645;
    float _S1647 = dist_coeffs_5[int(7)] * _S1640.y + _S1642 + dist_coeffs_5[int(6)] * _S1640.x + _S1643 + _S1066 * _S1644 + _S1065 * _S1645 + _S1064 * _S1646 + dist_coeffs_5[int(3)] * (r2_56 * _S1646);
    float _S1648 = v_56 * _S1647;
    float _S1649 = u_56 * _S1647;
    float2  _S1650 = make_float2 (radial_24) * _S1640 + make_float2 (_S1070 * (v_56 * _S1640.y) + _S1069 * _S1643 + 2.0f * (u_56 * _S1643) + _S1067 * (v_56 * _S1640.x) + _S1649 + _S1649, _S1072 * _S1642 + 2.0f * (v_56 * _S1642) + _S1071 * _S1640.y + _S1068 * _S1640.x + _S1648 + _S1648);
    float3  _S1651 = make_float3 (_S1638.x, _S1638.y, _S1635);
    float2  _S1652 = _S1051 * _S1650;
    float2  _S1653 = _S1062 * _S1650;
    float _S1654 = _S1652.x + _S1652.y;
    if(_S1055)
    {
        float _S1655 = _S1654 / _S1056;
        float _S1656 = _S1057 * - _S1655;
        float _S1657 = _S1054 * (0.3333333432674408f * - (_S1053 * _S1655));
        k_2 = _S1657 + _S1657;
        _S1056 = _S1656;
        _S1057 = 0.0f;
    }
    else
    {
        float _S1658 = _S1654 / _S1058;
        float _S1659 = _S1054 * - _S1658;
        k_2 = _S1052 * _S1658;
        _S1056 = 0.0f;
        _S1057 = _S1659;
    }
    DiffPair_float_0 _S1660;
    (&_S1660)->primal_0 = _S1052;
    (&_S1660)->differential_0 = 0.0f;
    DiffPair_float_0 _S1661;
    (&_S1661)->primal_0 = _S1053;
    (&_S1661)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1660, &_S1661, k_2);
    float _S1662 = _S1661.differential_0 + _S1056;
    float _S1663 = _S1660.differential_0 + _S1057;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1664;
    (&_S1664)->primal_0 = _S1051;
    (&_S1664)->differential_0 = _S1474;
    s_bwd_length_impl_0(&_S1664, _S1663);
    float2  _S1665 = _S1664.differential_0 + _S1653;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1666;
    (&_S1666)->primal_0 = R_5;
    (&_S1666)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1667;
    (&_S1667)->primal_0 = _S1050;
    (&_S1667)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1666, &_S1667, _S1356);
    DiffPair_float_0 _S1668;
    (&_S1668)->primal_0 = _S1047;
    (&_S1668)->differential_0 = 0.0f;
    DiffPair_float_0 _S1669;
    (&_S1669)->primal_0 = _S1049;
    (&_S1669)->differential_0 = 0.0f;
    _d_max_0(&_S1668, &_S1669, 0.0f);
    DiffPair_float_0 _S1670;
    (&_S1670)->primal_0 = _S1046;
    (&_S1670)->differential_0 = 0.0f;
    DiffPair_float_0 _S1671;
    (&_S1671)->primal_0 = _S1049;
    (&_S1671)->differential_0 = 0.0f;
    _d_min_0(&_S1670, &_S1671, 0.0f);
    float _S1672 = _S1669.differential_0 + _S1671.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1673;
    (&_S1673)->primal_0 = _S1048;
    (&_S1673)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1673, _S1672);
    float3  _S1674 = _S1673.differential_0 + _S1489;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1675;
    (&_S1675)->primal_0 = R_5;
    (&_S1675)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1676;
    (&_S1676)->primal_0 = pos_i_13;
    (&_S1676)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1675, &_S1676, _S1674);
    DiffPair_float_0 _S1677;
    (&_S1677)->primal_0 = _S1043;
    (&_S1677)->differential_0 = 0.0f;
    DiffPair_float_0 _S1678;
    (&_S1678)->primal_0 = _S1045;
    (&_S1678)->differential_0 = 0.0f;
    _d_max_0(&_S1677, &_S1678, _S1668.differential_0);
    DiffPair_float_0 _S1679;
    (&_S1679)->primal_0 = _S1042;
    (&_S1679)->differential_0 = 0.0f;
    DiffPair_float_0 _S1680;
    (&_S1680)->primal_0 = _S1045;
    (&_S1680)->differential_0 = 0.0f;
    _d_min_0(&_S1679, &_S1680, _S1670.differential_0);
    float _S1681 = _S1678.differential_0 + _S1680.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1682;
    (&_S1682)->primal_0 = _S1044;
    (&_S1682)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1682, _S1681);
    float3  _S1683 = _S1682.differential_0 + _S1516;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1684;
    (&_S1684)->primal_0 = R_5;
    (&_S1684)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1685;
    (&_S1685)->primal_0 = pos_i_12;
    (&_S1685)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1684, &_S1685, _S1683);
    DiffPair_float_0 _S1686;
    (&_S1686)->primal_0 = _S1039;
    (&_S1686)->differential_0 = 0.0f;
    DiffPair_float_0 _S1687;
    (&_S1687)->primal_0 = _S1041;
    (&_S1687)->differential_0 = 0.0f;
    _d_max_0(&_S1686, &_S1687, _S1677.differential_0);
    DiffPair_float_0 _S1688;
    (&_S1688)->primal_0 = _S1038;
    (&_S1688)->differential_0 = 0.0f;
    DiffPair_float_0 _S1689;
    (&_S1689)->primal_0 = _S1041;
    (&_S1689)->differential_0 = 0.0f;
    _d_min_0(&_S1688, &_S1689, _S1679.differential_0);
    float _S1690 = _S1687.differential_0 + _S1689.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1691;
    (&_S1691)->primal_0 = _S1040;
    (&_S1691)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1691, _S1690);
    float3  _S1692 = _S1691.differential_0 + _S1543;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1693;
    (&_S1693)->primal_0 = R_5;
    (&_S1693)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1694;
    (&_S1694)->primal_0 = pos_i_11;
    (&_S1694)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1693, &_S1694, _S1692);
    DiffPair_float_0 _S1695;
    (&_S1695)->primal_0 = _S1035;
    (&_S1695)->differential_0 = 0.0f;
    DiffPair_float_0 _S1696;
    (&_S1696)->primal_0 = _S1037;
    (&_S1696)->differential_0 = 0.0f;
    _d_max_0(&_S1695, &_S1696, _S1686.differential_0);
    DiffPair_float_0 _S1697;
    (&_S1697)->primal_0 = _S1034;
    (&_S1697)->differential_0 = 0.0f;
    DiffPair_float_0 _S1698;
    (&_S1698)->primal_0 = _S1037;
    (&_S1698)->differential_0 = 0.0f;
    _d_min_0(&_S1697, &_S1698, _S1688.differential_0);
    float _S1699 = _S1696.differential_0 + _S1698.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1700;
    (&_S1700)->primal_0 = _S1036;
    (&_S1700)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1700, _S1699);
    float3  _S1701 = _S1700.differential_0 + _S1570;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1702;
    (&_S1702)->primal_0 = R_5;
    (&_S1702)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1703;
    (&_S1703)->primal_0 = pos_i_10;
    (&_S1703)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1702, &_S1703, _S1701);
    DiffPair_float_0 _S1704;
    (&_S1704)->primal_0 = _S1031;
    (&_S1704)->differential_0 = 0.0f;
    DiffPair_float_0 _S1705;
    (&_S1705)->primal_0 = _S1033;
    (&_S1705)->differential_0 = 0.0f;
    _d_max_0(&_S1704, &_S1705, _S1695.differential_0);
    DiffPair_float_0 _S1706;
    (&_S1706)->primal_0 = _S1030;
    (&_S1706)->differential_0 = 0.0f;
    DiffPair_float_0 _S1707;
    (&_S1707)->primal_0 = _S1033;
    (&_S1707)->differential_0 = 0.0f;
    _d_min_0(&_S1706, &_S1707, _S1697.differential_0);
    float _S1708 = _S1705.differential_0 + _S1707.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1709;
    (&_S1709)->primal_0 = _S1032;
    (&_S1709)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1709, _S1708);
    float3  _S1710 = _S1709.differential_0 + _S1597;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1711;
    (&_S1711)->primal_0 = R_5;
    (&_S1711)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1712;
    (&_S1712)->primal_0 = pos_i_9;
    (&_S1712)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1711, &_S1712, _S1710);
    DiffPair_float_0 _S1713;
    (&_S1713)->primal_0 = _S1027;
    (&_S1713)->differential_0 = 0.0f;
    DiffPair_float_0 _S1714;
    (&_S1714)->primal_0 = _S1029;
    (&_S1714)->differential_0 = 0.0f;
    _d_max_0(&_S1713, &_S1714, _S1704.differential_0);
    DiffPair_float_0 _S1715;
    (&_S1715)->primal_0 = _S1026;
    (&_S1715)->differential_0 = 0.0f;
    DiffPair_float_0 _S1716;
    (&_S1716)->primal_0 = _S1029;
    (&_S1716)->differential_0 = 0.0f;
    _d_min_0(&_S1715, &_S1716, _S1706.differential_0);
    float _S1717 = _S1714.differential_0 + _S1716.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1718;
    (&_S1718)->primal_0 = _S1028;
    (&_S1718)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1718, _S1717);
    float3  _S1719 = _S1718.differential_0 + _S1624;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1720;
    (&_S1720)->primal_0 = R_5;
    (&_S1720)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1721;
    (&_S1721)->primal_0 = pos_i_8;
    (&_S1721)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1720, &_S1721, _S1719);
    DiffPair_float_0 _S1722;
    (&_S1722)->primal_0 = _S1023;
    (&_S1722)->differential_0 = 0.0f;
    DiffPair_float_0 _S1723;
    (&_S1723)->primal_0 = _S1025;
    (&_S1723)->differential_0 = 0.0f;
    _d_max_0(&_S1722, &_S1723, _S1713.differential_0);
    DiffPair_float_0 _S1724;
    (&_S1724)->primal_0 = _S1022;
    (&_S1724)->differential_0 = 0.0f;
    DiffPair_float_0 _S1725;
    (&_S1725)->primal_0 = _S1025;
    (&_S1725)->differential_0 = 0.0f;
    _d_min_0(&_S1724, &_S1725, _S1715.differential_0);
    float _S1726 = _S1723.differential_0 + _S1725.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1727;
    (&_S1727)->primal_0 = _S1024;
    (&_S1727)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1727, _S1726);
    float3  _S1728 = _S1727.differential_0 + _S1651;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1729;
    (&_S1729)->primal_0 = R_5;
    (&_S1729)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1730;
    (&_S1730)->primal_0 = pos_i_7;
    (&_S1730)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1729, &_S1730, _S1728);
    DiffPair_float_0 _S1731;
    (&_S1731)->primal_0 = 0.0f;
    (&_S1731)->differential_0 = 0.0f;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = _S1021;
    (&_S1732)->differential_0 = 0.0f;
    _d_max_0(&_S1731, &_S1732, _S1722.differential_0);
    DiffPair_float_0 _S1733;
    (&_S1733)->primal_0 = 1.00000001504746622e+30f;
    (&_S1733)->differential_0 = 0.0f;
    DiffPair_float_0 _S1734;
    (&_S1734)->primal_0 = _S1021;
    (&_S1734)->differential_0 = 0.0f;
    _d_min_0(&_S1733, &_S1734, _S1724.differential_0);
    float _S1735 = _S1732.differential_0 + _S1734.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1736;
    (&_S1736)->primal_0 = _S1020;
    (&_S1736)->differential_0 = _S1296;
    s_bwd_length_impl_1(&_S1736, _S1735);
    float3  _S1737 = _S1736.differential_0 + make_float3 (_S1665.x, _S1665.y, _S1662);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1738;
    (&_S1738)->primal_0 = R_5;
    (&_S1738)->differential_0 = _S1358;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1739;
    (&_S1739)->primal_0 = pos_5;
    (&_S1739)->differential_0 = _S1296;
    s_bwd_prop_mul_0(&_S1738, &_S1739, _S1737);
    float3  _S1740 = _S1356 + _S1674 + _S1683 + _S1692 + _S1701 + _S1710 + _S1719 + _S1728 + _S1737 + _S1361.differential_0;
    Matrix<float, 3, 3>  _S1741 = _S1666.differential_0 + _S1675.differential_0 + _S1684.differential_0 + _S1693.differential_0 + _S1702.differential_0 + _S1711.differential_0 + _S1720.differential_0 + _S1729.differential_0 + _S1738.differential_0 + _S1362;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_1)[int(0)] = _S1432;
    (*v_sh_coeffs_1)[int(1)] = _S1433;
    (*v_sh_coeffs_1)[int(2)] = _S1434;
    (*v_sh_coeffs_1)[int(3)] = _S1435;
    (*v_sh_coeffs_1)[int(4)] = _S1436;
    (*v_sh_coeffs_1)[int(5)] = _S1437;
    (*v_sh_coeffs_1)[int(6)] = _S1438;
    (*v_sh_coeffs_1)[int(7)] = _S1439;
    (*v_sh_coeffs_1)[int(8)] = _S1440;
    (*v_sh_coeffs_1)[int(9)] = _S1441;
    (*v_sh_coeffs_1)[int(10)] = _S1442;
    (*v_sh_coeffs_1)[int(11)] = _S1443;
    (*v_sh_coeffs_1)[int(12)] = _S1444;
    (*v_sh_coeffs_1)[int(13)] = _S1445;
    (*v_sh_coeffs_1)[int(14)] = _S1446;
    (*v_sh_coeffs_1)[int(15)] = _S1447;
    *v_R_1 = _S1741;
    *v_t_1 = _S1740;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  dOut_7)
{
    float3  _S1742 = _slang_select(((*dpx_8).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_8).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_7;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S1742;
    return;
}

inline __device__ float3  abs_0(float3  x_15)
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
        *_slang_vector_get_element_ptr(&result_8, i_4) = (F32_abs((_slang_vector_get_element(x_15, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_8;
}

inline __device__ bool ray_aabb_intersection(float3  ray_o_0, float3  ray_d_0, float3  center_0, float radius_0, float * t0_0, float * t1_0)
{
    float3  m_1 = make_float3 (1.0f) / ray_d_0;
    float3  k_3 = abs_0(m_1) * make_float3 (radius_0);
    float3  _S1743 = - (m_1 * (ray_o_0 - center_0));
    float3  ta_0 = _S1743 - k_3;
    float3  tb_0 = _S1743 + k_3;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S1744 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S1744;
    return (*t0_0) < _S1744;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_9, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_8)
{
    float _S1745 = (1.0f - (*dps_0).primal_0) * dOut_8;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S1745;
    DiffPair_float_0 _S1746 = *dpy_4;
    float _S1747 = (*dps_0).primal_0 * dOut_8;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S1747;
    float _S1748 = (_S1746.primal_0 - (*dpx_9).primal_0) * dOut_8;
    dps_0->primal_0 = _S1746.primal_0;
    dps_0->differential_0 = _S1748;
    return;
}

inline __device__ float lerp_0(float x_16, float y_11, float s_0)
{
    return x_16 + (y_11 - x_16) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_0)
{
    float _S1749 = w_0.z;
    float _S1750 = 1.0f - _S1749;
    float _S1751 = w_0.y;
    float _S1752 = 1.0f - _S1751;
    float _S1753 = _S1750 * _S1752;
    float _S1754 = w_0.x;
    float _S1755 = 1.0f - _S1754;
    float _S1756 = _S1750 * _S1751;
    float _S1757 = _S1749 * _S1752;
    float _S1758 = _S1749 * _S1751;
    return _S1753 * _S1755 * (*densities_6)[int(0)] + _S1753 * _S1754 * (*densities_6)[int(1)] + _S1756 * _S1755 * (*densities_6)[int(2)] + _S1756 * _S1754 * (*densities_6)[int(3)] + _S1757 * _S1755 * (*densities_6)[int(4)] + _S1757 * _S1754 * (*densities_6)[int(5)] + _S1758 * _S1755 * (*densities_6)[int(6)] + _S1758 * _S1754 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  densities_7, float3  ray_o_1, float3  ray_d_1)
{
    float _S1759 = 0.5f * size_6;
    float3  m_2 = make_float3 (1.0f) / ray_d_1;
    float3  k_4 = abs_0(m_2) * make_float3 (_S1759);
    float3  _S1760 = - (m_2 * (ray_o_1 - (pos_6 + make_float3 (_S1759))));
    float3  ta_1 = _S1760 - k_4;
    float3  tb_1 = _S1760 + k_4;
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
        float3  _S1761 = (ray_o_1 + ray_d_1 * make_float3 (lerp_0(t0_1, t1_1, (float(i_5) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6);
        FixedArray<float, 8>  _S1762 = densities_7;
        float _S1763 = interp_0(&_S1762, _S1761);
        float _S1764;
        if(_S1763 > 1.10000002384185791f)
        {
            _S1764 = _S1763;
        }
        else
        {
            _S1764 = (F32_exp((0.90909093618392944f * _S1763 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S1764;
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
    float _S1765;
};

inline __device__ float3  s_primal_ctx_abs_0(float3  _S1766)
{
    return abs_0(_S1766);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1767, float _S1768, float _S1769)
{
    return lerp_0(_S1767, _S1768, _S1769);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S1770 = dpw_0.z;
    float _S1771 = 1.0f - _S1770;
    float _S1772 = dpw_0.y;
    float _S1773 = 1.0f - _S1772;
    float _S1774 = _S1771 * _S1773;
    float _S1775 = dpw_0.x;
    float _S1776 = 1.0f - _S1775;
    float _S1777 = _S1771 * _S1772;
    float _S1778 = _S1770 * _S1773;
    float _S1779 = _S1770 * _S1772;
    return _S1774 * _S1776 * (*dpdensities_0)[int(0)] + _S1774 * _S1775 * (*dpdensities_0)[int(1)] + _S1777 * _S1776 * (*dpdensities_0)[int(2)] + _S1777 * _S1775 * (*dpdensities_0)[int(3)] + _S1778 * _S1776 * (*dpdensities_0)[int(4)] + _S1778 * _S1775 * (*dpdensities_0)[int(5)] + _S1779 * _S1776 * (*dpdensities_0)[int(6)] + _S1779 * _S1775 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_exp_0(float _S1780)
{
    return (F32_exp((_S1780)));
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_0, float3  dpray_d_0, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S1765 = 0.0f;
    _s_diff_ctx_0->_S1765 = 0.0f;
    float _S1781 = 0.5f * size_7;
    float3  m_3 = make_float3 (1.0f) / dpray_d_0;
    float3  k_5 = s_primal_ctx_abs_0(m_3) * make_float3 (_S1781);
    float3  _S1782 = - (m_3 * (dpray_o_0 - (pos_7 + make_float3 (_S1781))));
    float3  ta_2 = _S1782 - k_5;
    float3  tb_2 = _S1782 + k_5;
    float _S1783 = (F32_max(((F32_max((ta_2.x), (ta_2.y)))), ((F32_max((ta_2.z), (0.0f))))));
    float _S1784 = (F32_min(((F32_min((tb_2.x), (tb_2.y)))), (tb_2.z)));
    float accum_2;
    if(!!(_S1783 < _S1784))
    {
        float _S1785 = - (_S1784 - _S1783);
        bool _runFlag_0 = true;
        int i_6 = int(0);
        accum_2 = 0.0f;
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_0->_S1765 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            int _S1786;
            float _S1787;
            if(i_6 < int(8))
            {
                float _S1788 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_0 + dpray_d_0 * make_float3 (s_primal_ctx_lerp_0(_S1783, _S1784, (float(i_6) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S1788 > 1.10000002384185791f)
                {
                    _S1787 = _S1788;
                }
                else
                {
                    _S1787 = s_primal_ctx_exp_0(0.90909093618392944f * _S1788 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S1787;
                _S1786 = int(2);
                _S1787 = accum_3;
            }
            else
            {
                _S1786 = int(1);
                _S1787 = 0.0f;
            }
            if(_S1786 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                i_6 = i_6 + int(1);
                accum_2 = _S1787;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = (F32_min((1.0f - s_primal_ctx_exp_0(_S1785 / 8.0f * accum_2)), (0.99900001287460327f)));
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1789, float _S1790)
{
    _d_exp_0(_S1789, _S1790);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_2)
{
    float _S1791 = (*dpw_1).primal_0.z;
    float _S1792 = 1.0f - _S1791;
    float _S1793 = (*dpw_1).primal_0.y;
    float _S1794 = 1.0f - _S1793;
    float _S1795 = _S1792 * _S1794;
    float _S1796 = (*dpw_1).primal_0.x;
    float _S1797 = 1.0f - _S1796;
    float _S1798 = _S1792 * _S1793;
    float _S1799 = _S1791 * _S1794;
    float _S1800 = _S1791 * _S1793;
    float _S1801 = _S1800 * _S1796 * _s_dOut_2;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_2;
    float _S1802 = _S1800 * _S1797 * _s_dOut_2;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_2;
    float _S1803 = _S1799 * _S1796 * _s_dOut_2;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_2;
    float _S1804 = _S1799 * _S1797 * _s_dOut_2;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_2;
    float _S1805 = _S1798 * _S1796 * _s_dOut_2;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_2;
    float _S1806 = _S1798 * _S1797 * _s_dOut_2;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_2;
    float _S1807 = _S1795 * _S1796 * _s_dOut_2;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_2;
    float _S1808 = _S1795 * _S1797 * _s_dOut_2;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_2;
    float _S1809 = _S1796 * s_diff_w7_T_0 + _S1797 * s_diff_w6_T_0;
    float _S1810 = _S1796 * s_diff_w5_T_0 + _S1797 * s_diff_w4_T_0;
    float _S1811 = _S1796 * s_diff_w3_T_0 + _S1797 * s_diff_w2_T_0;
    float _S1812 = _S1796 * s_diff_w1_T_0 + _S1797 * s_diff_w0_T_0;
    float3  _S1813 = make_float3 (_S1800 * s_diff_w7_T_0 + _S1799 * s_diff_w5_T_0 + _S1798 * s_diff_w3_T_0 + _S1795 * s_diff_w1_T_0 + - (_S1800 * s_diff_w6_T_0 + _S1799 * s_diff_w4_T_0 + _S1798 * s_diff_w2_T_0 + _S1795 * s_diff_w0_T_0), _S1791 * _S1809 + _S1792 * _S1811 + - (_S1791 * _S1810 + _S1792 * _S1812), _S1793 * _S1809 + _S1794 * _S1810 + - (_S1793 * _S1811 + _S1794 * _S1812));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S1813;
    FixedArray<float, 8>  _S1814;
    _S1814[int(0)] = 0.0f;
    _S1814[int(1)] = 0.0f;
    _S1814[int(2)] = 0.0f;
    _S1814[int(3)] = 0.0f;
    _S1814[int(4)] = 0.0f;
    _S1814[int(5)] = 0.0f;
    _S1814[int(6)] = 0.0f;
    _S1814[int(7)] = 0.0f;
    _S1814[int(7)] = _S1801;
    _S1814[int(6)] = _S1802;
    _S1814[int(5)] = _S1803;
    _S1814[int(4)] = _S1804;
    _S1814[int(3)] = _S1805;
    _S1814[int(2)] = _S1806;
    _S1814[int(1)] = _S1807;
    _S1814[int(0)] = _S1808;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S1814;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S1815, DiffPair_float_0 * _S1816, DiffPair_float_0 * _S1817, float _S1818)
{
    _d_lerp_0(_S1815, _S1816, _S1817, _S1818);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1819, float3  _S1820)
{
    _d_abs_vector_0(_S1819, _S1820);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float _s_dOut_3, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_1)
{
    FixedArray<float, 8>  _S1821 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1822 = *dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1823 = *dpray_d_1;
    float _S1824 = 0.5f * size_8;
    float3  _S1825 = make_float3 (_S1824);
    float3  _S1826 = make_float3 (size_8);
    float3  m_4 = make_float3 (1.0f) / (*dpray_d_1).primal_0;
    float3  _S1827 = (*dpray_d_1).primal_0 * (*dpray_d_1).primal_0;
    float3  _S1828 = (*dpray_o_1).primal_0 - (pos_8 + make_float3 (_S1824));
    float3  k_6 = s_primal_ctx_abs_0(m_4) * make_float3 (_S1824);
    float3  _S1829 = - (m_4 * _S1828);
    float3  ta_3 = _S1829 - k_6;
    float3  tb_3 = _S1829 + k_6;
    float _S1830 = ta_3.x;
    float _S1831 = ta_3.y;
    float _S1832 = (F32_max((_S1830), (_S1831)));
    float _S1833 = ta_3.z;
    float _S1834 = (F32_max((_S1833), (0.0f)));
    float _S1835 = (F32_max((_S1832), (_S1834)));
    float _S1836 = tb_3.x;
    float _S1837 = tb_3.y;
    float _S1838 = (F32_min((_S1836), (_S1837)));
    float _S1839 = tb_3.z;
    float _S1840 = (F32_min((_S1838), (_S1839)));
    bool _S1841 = !!(_S1835 < _S1840);
    float _S1842;
    float _S1843;
    float _S1844;
    if(_S1841)
    {
        float _S1845 = - (_S1840 - _S1835) / 8.0f;
        float _S1846 = _S1845 * _s_diff_ctx_1->_S1765;
        _S1842 = 1.0f - s_primal_ctx_exp_0(_S1846);
        _S1843 = _S1846;
        _S1844 = _S1845;
    }
    else
    {
        _S1842 = 0.0f;
        _S1843 = 0.0f;
        _S1844 = 0.0f;
    }
    float3  _S1847 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1848 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<float, 8>  _S1849;
    float3  _S1850;
    float3  _S1851;
    if(_S1841)
    {
        DiffPair_float_0 _S1852;
        (&_S1852)->primal_0 = _S1842;
        (&_S1852)->differential_0 = 0.0f;
        DiffPair_float_0 _S1853;
        (&_S1853)->primal_0 = 0.99900001287460327f;
        (&_S1853)->differential_0 = 0.0f;
        _d_min_0(&_S1852, &_S1853, _s_dOut_3);
        float _S1854 = - _S1852.differential_0;
        DiffPair_float_0 _S1855;
        (&_S1855)->primal_0 = _S1843;
        (&_S1855)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S1855, _S1854);
        float _S1856 = _S1844 * _S1855.differential_0;
        float _S1857 = 0.125f * (_s_diff_ctx_1->_S1765 * _S1855.differential_0);
        int _dc_0 = int(8);
        _S1842 = _S1856;
        _S1849[int(0)] = 0.0f;
        _S1849[int(1)] = 0.0f;
        _S1849[int(2)] = 0.0f;
        _S1849[int(3)] = 0.0f;
        _S1849[int(4)] = 0.0f;
        _S1849[int(5)] = 0.0f;
        _S1849[int(6)] = 0.0f;
        _S1849[int(7)] = 0.0f;
        _S1850 = _S1847;
        _S1851 = _S1847;
        _S1843 = 0.0f;
        _S1844 = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S1858 = _dc_0 < int(8);
            float _S1859;
            float _S1860;
            int _S1861;
            float3  _S1862;
            float3  _S1863;
            bool _S1864;
            if(_S1858)
            {
                float _S1865 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S1866 = s_primal_ctx_lerp_0(_S1835, _S1840, _S1865);
                float3  _S1867 = make_float3 (_S1866);
                float3  _S1868 = (_S1822.primal_0 + _S1823.primal_0 * make_float3 (_S1866) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S1869 = _S1821;
                float _S1870 = s_primal_ctx_interp_0(&_S1869, _S1868);
                bool _S1871 = _S1870 > 1.10000002384185791f;
                if(_S1871)
                {
                    _S1859 = 0.0f;
                }
                else
                {
                    _S1859 = 0.90909093618392944f * _S1870 - 0.90468984842300415f;
                }
                _S1861 = int(2);
                _S1864 = _S1871;
                _S1862 = _S1868;
                _S1863 = _S1867;
                _S1860 = _S1865;
            }
            else
            {
                _S1861 = int(1);
                _S1864 = false;
                _S1859 = 0.0f;
                _S1862 = _S1847;
                _S1863 = _S1847;
                _S1860 = 0.0f;
            }
            float _S1872;
            float _S1873;
            if(!(_S1861 != int(2)))
            {
                _S1872 = _S1842;
                _S1873 = 0.0f;
            }
            else
            {
                _S1872 = 0.0f;
                _S1873 = _S1842;
            }
            if(_S1858)
            {
                float _S1874 = _S1872 + _S1873;
                float _S1875;
                if(_S1864)
                {
                    _S1875 = _S1872;
                }
                else
                {
                    DiffPair_float_0 _S1876;
                    (&_S1876)->primal_0 = _S1859;
                    (&_S1876)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S1876, _S1872);
                    _S1875 = 0.90909093618392944f * _S1876.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S1877;
                (&_S1877)->primal_0 = _S1821;
                (&_S1877)->differential_0 = _S1848;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1878;
                (&_S1878)->primal_0 = _S1862;
                (&_S1878)->differential_0 = _S1847;
                s_bwd_prop_interp_0(&_S1877, &_S1878, _S1875);
                float3  _S1879 = _S1878.differential_0 / _S1826;
                float3  _S1880 = _S1823.primal_0 * _S1879;
                float3  _S1881 = _S1863 * _S1879;
                float _S1882 = _S1880.x + _S1880.y + _S1880.z;
                DiffPair_float_0 _S1883;
                (&_S1883)->primal_0 = _S1835;
                (&_S1883)->differential_0 = 0.0f;
                DiffPair_float_0 _S1884;
                (&_S1884)->primal_0 = _S1840;
                (&_S1884)->differential_0 = 0.0f;
                DiffPair_float_0 _S1885;
                (&_S1885)->primal_0 = _S1860;
                (&_S1885)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S1883, &_S1884, &_S1885, _S1882);
                float _S1886 = (&_S1877)->differential_0[int(0)] + _S1849[int(0)];
                float _S1887 = (&_S1877)->differential_0[int(1)] + _S1849[int(1)];
                float _S1888 = (&_S1877)->differential_0[int(2)] + _S1849[int(2)];
                float _S1889 = (&_S1877)->differential_0[int(3)] + _S1849[int(3)];
                float _S1890 = (&_S1877)->differential_0[int(4)] + _S1849[int(4)];
                float _S1891 = (&_S1877)->differential_0[int(5)] + _S1849[int(5)];
                float _S1892 = (&_S1877)->differential_0[int(6)] + _S1849[int(6)];
                float _S1893 = (&_S1877)->differential_0[int(7)] + _S1849[int(7)];
                float3  _S1894 = _S1879 + _S1850;
                float3  _S1895 = _S1881 + _S1851;
                float _S1896 = _S1884.differential_0 + _S1843;
                float _S1897 = _S1883.differential_0 + _S1844;
                _S1842 = _S1874;
                _S1849[int(0)] = _S1886;
                _S1849[int(1)] = _S1887;
                _S1849[int(2)] = _S1888;
                _S1849[int(3)] = _S1889;
                _S1849[int(4)] = _S1890;
                _S1849[int(5)] = _S1891;
                _S1849[int(6)] = _S1892;
                _S1849[int(7)] = _S1893;
                _S1850 = _S1894;
                _S1851 = _S1895;
                _S1843 = _S1896;
                _S1844 = _S1897;
            }
            else
            {
                _S1842 = _S1873;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S1898 = - _S1857;
        float _S1899 = - _S1898 + _S1844;
        float3  _S1900 = _S1850;
        _S1842 = _S1898 + _S1843;
        _S1843 = _S1899;
        _S1850 = _S1851;
        _S1851 = _S1900;
    }
    else
    {
        _S1842 = 0.0f;
        _S1843 = 0.0f;
        _S1850 = _S1847;
        _S1851 = _S1847;
        _S1849[int(0)] = 0.0f;
        _S1849[int(1)] = 0.0f;
        _S1849[int(2)] = 0.0f;
        _S1849[int(3)] = 0.0f;
        _S1849[int(4)] = 0.0f;
        _S1849[int(5)] = 0.0f;
        _S1849[int(6)] = 0.0f;
        _S1849[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S1901;
    (&_S1901)->primal_0 = _S1838;
    (&_S1901)->differential_0 = 0.0f;
    DiffPair_float_0 _S1902;
    (&_S1902)->primal_0 = _S1839;
    (&_S1902)->differential_0 = 0.0f;
    _d_min_0(&_S1901, &_S1902, _S1842);
    DiffPair_float_0 _S1903;
    (&_S1903)->primal_0 = _S1836;
    (&_S1903)->differential_0 = 0.0f;
    DiffPair_float_0 _S1904;
    (&_S1904)->primal_0 = _S1837;
    (&_S1904)->differential_0 = 0.0f;
    _d_min_0(&_S1903, &_S1904, _S1901.differential_0);
    DiffPair_float_0 _S1905;
    (&_S1905)->primal_0 = _S1832;
    (&_S1905)->differential_0 = 0.0f;
    DiffPair_float_0 _S1906;
    (&_S1906)->primal_0 = _S1834;
    (&_S1906)->differential_0 = 0.0f;
    _d_max_0(&_S1905, &_S1906, _S1843);
    DiffPair_float_0 _S1907;
    (&_S1907)->primal_0 = _S1833;
    (&_S1907)->differential_0 = 0.0f;
    DiffPair_float_0 _S1908;
    (&_S1908)->primal_0 = 0.0f;
    (&_S1908)->differential_0 = 0.0f;
    _d_max_0(&_S1907, &_S1908, _S1906.differential_0);
    DiffPair_float_0 _S1909;
    (&_S1909)->primal_0 = _S1830;
    (&_S1909)->differential_0 = 0.0f;
    DiffPair_float_0 _S1910;
    (&_S1910)->primal_0 = _S1831;
    (&_S1910)->differential_0 = 0.0f;
    _d_max_0(&_S1909, &_S1910, _S1905.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S1903.differential_0, _S1904.differential_0, _S1902.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S1909.differential_0, _S1910.differential_0, _S1907.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S1911 = _S1825 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1912;
    (&_S1912)->primal_0 = m_4;
    (&_S1912)->differential_0 = _S1847;
    s_bwd_prop_abs_0(&_S1912, _S1911);
    float3  _S1913 = m_4 * s_diff_n_T_0;
    float3  _S1914 = - ((_S1912.differential_0 + _S1828 * s_diff_n_T_0) / _S1827) + _S1850;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1914;
    float3  _S1915 = _S1913 + _S1851;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1915;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S1849;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S1916, float _S1917, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S1918, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1919, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1920, float _S1921)
{
    FixedArray<float, 8>  _S1922 = _S1918->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1923;
    float _S1924 = s_primal_ctx_evaluate_alpha_voxel_0(_S1916, _S1917, &_S1922, (*_S1919).primal_0, (*_S1920).primal_0, &_S1923);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1925 = _S1923;
    s_bwd_prop_evaluate_alpha_voxel_0(_S1916, _S1917, _S1918, _S1919, _S1920, _S1921, &_S1925);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  densities_8, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    FixedArray<float, 8>  _S1926 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = densities_8;
    (&dp_densities_0)->differential_0 = _S1926;
    float3  _S1927 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1927;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1927;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_0 = dp_ray_o_0.differential_0;
    *v_ray_d_0 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  densities_9, float3  rgb_0, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_4)
{
    *out_rgb_0 = rgb_0;
    float _S1928 = 0.5f * size_10;
    float3  m_5 = make_float3 (1.0f) / ray_d_3;
    float3  k_7 = abs_0(m_5) * make_float3 (_S1928);
    float3  _S1929 = - (m_5 * (ray_o_3 - (pos_10 + make_float3 (_S1928))));
    float3  ta_4 = _S1929 - k_7;
    float3  tb_4 = _S1929 + k_7;
    float _S1930 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S1931 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
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
        float t_6 = lerp_0(_S1930, _S1931, (float(i_7) + 0.5f) / 8.0f);
        float3  _S1932 = (ray_o_3 + ray_d_3 * make_float3 (t_6) - pos_10) / make_float3 (size_10);
        FixedArray<float, 8>  _S1933 = densities_9;
        float _S1934 = interp_0(&_S1933, _S1932);
        float _S1935;
        if(_S1934 > 1.10000002384185791f)
        {
            _S1935 = _S1934;
        }
        else
        {
            _S1935 = (F32_exp((0.90909093618392944f * _S1934 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S1935;
        float depth_accum_1 = depth_accum_0 + t_6 * _S1935;
        i_7 = i_7 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_4 = (F32_max((depth_accum_0 / accum_4), (0.0f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S1936;
    float _S1937;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_0, float3  dpray_o_2, float3  dpray_d_2, float3  * dpout_rgb_0, float * dpdepth_0, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_2)
{
    _s_diff_ctx_2->_S1936 = 0.0f;
    _s_diff_ctx_2->_S1937 = 0.0f;
    float _S1938 = 0.5f * size_11;
    float3  m_6 = make_float3 (1.0f) / dpray_d_2;
    float3  k_8 = s_primal_ctx_abs_0(m_6) * make_float3 (_S1938);
    float3  _S1939 = - (m_6 * (dpray_o_2 - (pos_11 + make_float3 (_S1938))));
    float3  ta_5 = _S1939 - k_8;
    float3  tb_5 = _S1939 + k_8;
    float _S1940 = (F32_max(((F32_max((ta_5.x), (ta_5.y)))), ((F32_max((ta_5.z), (0.0f))))));
    float _S1941 = (F32_min(((F32_min((tb_5.x), (tb_5.y)))), (tb_5.z)));
    bool _runFlag_1 = true;
    int i_8 = int(0);
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_2->_S1936 = depth_accum_2;
        _s_diff_ctx_2->_S1937 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        int _S1942;
        float _S1943;
        float _S1944;
        if(i_8 < int(8))
        {
            float _S1945 = s_primal_ctx_lerp_0(_S1940, _S1941, (float(i_8) + 0.5f) / 8.0f);
            float _S1946 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_2 + dpray_d_2 * make_float3 (_S1945) - pos_11) / make_float3 (size_11));
            if(_S1946 > 1.10000002384185791f)
            {
                _S1943 = _S1946;
            }
            else
            {
                _S1943 = s_primal_ctx_exp_0(0.90909093618392944f * _S1946 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S1943;
            float depth_accum_3 = depth_accum_2 + _S1945 * _S1943;
            _S1942 = int(1);
            _S1943 = accum_7;
            _S1944 = depth_accum_3;
        }
        else
        {
            _S1942 = int(0);
            _S1943 = 0.0f;
            _S1944 = 0.0f;
        }
        if(_S1942 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            i_8 = i_8 + int(1);
            accum_6 = _S1943;
            depth_accum_2 = _S1944;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S1947 = (F32_max((depth_accum_2 / accum_6), (0.0f)));
    *dpout_rgb_0 = dprgb_0;
    *dpdepth_0 = _S1947;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_1, float dpdepth_1, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_3)
{
    FixedArray<float, 8>  _S1948 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1949 = *dpray_o_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1950 = *dpray_d_3;
    float3  _S1951 = make_float3 (size_12);
    float _S1952 = 0.5f * size_12;
    float3  _S1953 = make_float3 (_S1952);
    float3  m_7 = make_float3 (1.0f) / (*dpray_d_3).primal_0;
    float3  _S1954 = (*dpray_d_3).primal_0 * (*dpray_d_3).primal_0;
    float3  _S1955 = (*dpray_o_3).primal_0 - (pos_12 + make_float3 (_S1952));
    float3  k_9 = s_primal_ctx_abs_0(m_7) * make_float3 (_S1952);
    float3  _S1956 = - (m_7 * _S1955);
    float3  ta_6 = _S1956 - k_9;
    float3  tb_6 = _S1956 + k_9;
    float _S1957 = ta_6.x;
    float _S1958 = ta_6.y;
    float _S1959 = (F32_max((_S1957), (_S1958)));
    float _S1960 = ta_6.z;
    float _S1961 = (F32_max((_S1960), (0.0f)));
    float _S1962 = (F32_max((_S1959), (_S1961)));
    float _S1963 = tb_6.x;
    float _S1964 = tb_6.y;
    float _S1965 = (F32_min((_S1963), (_S1964)));
    float _S1966 = tb_6.z;
    float _S1967 = (F32_min((_S1965), (_S1966)));
    float _S1968 = _s_diff_ctx_3->_S1937 * _s_diff_ctx_3->_S1937;
    float3  _S1969 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1970 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S1971;
    (&_S1971)->primal_0 = _s_diff_ctx_3->_S1936 / _s_diff_ctx_3->_S1937;
    (&_S1971)->differential_0 = 0.0f;
    DiffPair_float_0 _S1972;
    (&_S1972)->primal_0 = 0.0f;
    (&_S1972)->differential_0 = 0.0f;
    _d_max_0(&_S1971, &_S1972, dpdepth_1);
    float _S1973 = _S1971.differential_0 / _S1968;
    float _S1974 = _s_diff_ctx_3->_S1936 * - _S1973;
    float _S1975 = _s_diff_ctx_3->_S1937 * _S1973;
    int _dc_1 = int(8);
    float _S1976 = _S1974;
    float _S1977 = _S1975;
    FixedArray<float, 8>  _S1978;
    _S1978[int(0)] = 0.0f;
    _S1978[int(1)] = 0.0f;
    _S1978[int(2)] = 0.0f;
    _S1978[int(3)] = 0.0f;
    _S1978[int(4)] = 0.0f;
    _S1978[int(5)] = 0.0f;
    _S1978[int(6)] = 0.0f;
    _S1978[int(7)] = 0.0f;
    float3  _S1979 = _S1969;
    float3  _S1980 = _S1969;
    float _S1981 = 0.0f;
    float _S1982 = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S1983 = _dc_1 < int(8);
        int _S1984;
        float _S1985;
        float _S1986;
        float _S1987;
        float _S1988;
        float3  _S1989;
        float3  _S1990;
        bool _S1991;
        if(_S1983)
        {
            float _S1992 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S1993 = s_primal_ctx_lerp_0(_S1962, _S1967, _S1992);
            float3  _S1994 = make_float3 (_S1993);
            float3  _S1995 = (_S1949.primal_0 + _S1950.primal_0 * make_float3 (_S1993) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S1996 = _S1948;
            float _S1997 = s_primal_ctx_interp_0(&_S1996, _S1995);
            bool _S1998 = _S1997 > 1.10000002384185791f;
            if(_S1998)
            {
                _S1985 = _S1997;
                _S1986 = 0.0f;
            }
            else
            {
                float _S1999 = 0.90909093618392944f * _S1997 - 0.90468984842300415f;
                _S1985 = s_primal_ctx_exp_0(_S1999);
                _S1986 = _S1999;
            }
            float _S2000 = _S1985;
            float _S2001 = _S1986;
            _S1984 = int(1);
            _S1985 = _S1993;
            _S1986 = _S2000;
            _S1991 = _S1998;
            _S1987 = _S2001;
            _S1989 = _S1995;
            _S1990 = _S1994;
            _S1988 = _S1992;
        }
        else
        {
            _S1984 = int(0);
            _S1985 = 0.0f;
            _S1986 = 0.0f;
            _S1991 = false;
            _S1987 = 0.0f;
            _S1989 = _S1969;
            _S1990 = _S1969;
            _S1988 = 0.0f;
        }
        float _S2002;
        float _S2003;
        float _S2004;
        float _S2005;
        if(!(_S1984 != int(1)))
        {
            _S2002 = _S1976;
            _S2003 = _S1977;
            _S2004 = 0.0f;
            _S2005 = 0.0f;
        }
        else
        {
            _S2002 = 0.0f;
            _S2003 = 0.0f;
            _S2004 = _S1977;
            _S2005 = _S1976;
        }
        if(_S1983)
        {
            float _S2006 = _S1986 * _S2003;
            float _S2007 = _S2003 + _S2004;
            float _S2008 = _S1985 * _S2003 + _S2002;
            float _S2009 = _S2002 + _S2005;
            float _S2010;
            if(_S1991)
            {
                _S2010 = _S2008;
            }
            else
            {
                DiffPair_float_0 _S2011;
                (&_S2011)->primal_0 = _S1987;
                (&_S2011)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S2011, _S2008);
                _S2010 = 0.90909093618392944f * _S2011.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S2012;
            (&_S2012)->primal_0 = _S1948;
            (&_S2012)->differential_0 = _S1970;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2013;
            (&_S2013)->primal_0 = _S1989;
            (&_S2013)->differential_0 = _S1969;
            s_bwd_prop_interp_0(&_S2012, &_S2013, _S2010);
            float3  _S2014 = _S2013.differential_0 / _S1951;
            float3  _S2015 = _S1950.primal_0 * _S2014;
            float3  _S2016 = _S1990 * _S2014;
            float _S2017 = _S2015.x + _S2015.y + _S2015.z + _S2006;
            DiffPair_float_0 _S2018;
            (&_S2018)->primal_0 = _S1962;
            (&_S2018)->differential_0 = 0.0f;
            DiffPair_float_0 _S2019;
            (&_S2019)->primal_0 = _S1967;
            (&_S2019)->differential_0 = 0.0f;
            DiffPair_float_0 _S2020;
            (&_S2020)->primal_0 = _S1988;
            (&_S2020)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S2018, &_S2019, &_S2020, _S2017);
            float _S2021 = (&_S2012)->differential_0[int(0)] + _S1978[int(0)];
            float _S2022 = (&_S2012)->differential_0[int(1)] + _S1978[int(1)];
            float _S2023 = (&_S2012)->differential_0[int(2)] + _S1978[int(2)];
            float _S2024 = (&_S2012)->differential_0[int(3)] + _S1978[int(3)];
            float _S2025 = (&_S2012)->differential_0[int(4)] + _S1978[int(4)];
            float _S2026 = (&_S2012)->differential_0[int(5)] + _S1978[int(5)];
            float _S2027 = (&_S2012)->differential_0[int(6)] + _S1978[int(6)];
            float _S2028 = (&_S2012)->differential_0[int(7)] + _S1978[int(7)];
            float3  _S2029 = _S2014 + _S1979;
            float3  _S2030 = _S2016 + _S1980;
            float _S2031 = _S2019.differential_0 + _S1981;
            float _S2032 = _S2018.differential_0 + _S1982;
            _S1976 = _S2009;
            _S1977 = _S2007;
            _S1978[int(0)] = _S2021;
            _S1978[int(1)] = _S2022;
            _S1978[int(2)] = _S2023;
            _S1978[int(3)] = _S2024;
            _S1978[int(4)] = _S2025;
            _S1978[int(5)] = _S2026;
            _S1978[int(6)] = _S2027;
            _S1978[int(7)] = _S2028;
            _S1979 = _S2029;
            _S1980 = _S2030;
            _S1981 = _S2031;
            _S1982 = _S2032;
        }
        else
        {
            _S1976 = _S2005;
            _S1977 = _S2004;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S2033;
    (&_S2033)->primal_0 = _S1965;
    (&_S2033)->differential_0 = 0.0f;
    DiffPair_float_0 _S2034;
    (&_S2034)->primal_0 = _S1966;
    (&_S2034)->differential_0 = 0.0f;
    _d_min_0(&_S2033, &_S2034, _S1981);
    DiffPair_float_0 _S2035;
    (&_S2035)->primal_0 = _S1963;
    (&_S2035)->differential_0 = 0.0f;
    DiffPair_float_0 _S2036;
    (&_S2036)->primal_0 = _S1964;
    (&_S2036)->differential_0 = 0.0f;
    _d_min_0(&_S2035, &_S2036, _S2033.differential_0);
    DiffPair_float_0 _S2037;
    (&_S2037)->primal_0 = _S1959;
    (&_S2037)->differential_0 = 0.0f;
    DiffPair_float_0 _S2038;
    (&_S2038)->primal_0 = _S1961;
    (&_S2038)->differential_0 = 0.0f;
    _d_max_0(&_S2037, &_S2038, _S1982);
    DiffPair_float_0 _S2039;
    (&_S2039)->primal_0 = _S1960;
    (&_S2039)->differential_0 = 0.0f;
    DiffPair_float_0 _S2040;
    (&_S2040)->primal_0 = 0.0f;
    (&_S2040)->differential_0 = 0.0f;
    _d_max_0(&_S2039, &_S2040, _S2038.differential_0);
    DiffPair_float_0 _S2041;
    (&_S2041)->primal_0 = _S1957;
    (&_S2041)->differential_0 = 0.0f;
    DiffPair_float_0 _S2042;
    (&_S2042)->primal_0 = _S1958;
    (&_S2042)->differential_0 = 0.0f;
    _d_max_0(&_S2041, &_S2042, _S2037.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S2035.differential_0, _S2036.differential_0, _S2034.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S2041.differential_0, _S2042.differential_0, _S2039.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S2043 = _S1953 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2044;
    (&_S2044)->primal_0 = m_7;
    (&_S2044)->differential_0 = _S1969;
    s_bwd_prop_abs_0(&_S2044, _S2043);
    float3  _S2045 = m_7 * s_diff_n_T_1;
    float3  _S2046 = - ((_S2044.differential_0 + _S1955 * s_diff_n_T_1) / _S1954) + _S1980;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2046;
    float3  _S2047 = _S2045 + _S1979;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2047;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = dpout_rgb_1;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S1978;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S2048, float _S2049, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S2050, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2051, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2052, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2053, float3  _S2054, float _S2055)
{
    FixedArray<float, 8>  _S2056 = _S2050->primal_0;
    float3  _S2057;
    float _S2058;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2059;
    s_primal_ctx_evaluate_color_voxel_0(_S2048, _S2049, &_S2056, (*_S2051).primal_0, (*_S2052).primal_0, (*_S2053).primal_0, &_S2057, &_S2058, &_S2059);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2060 = _S2059;
    s_bwd_prop_evaluate_color_voxel_0(_S2048, _S2049, _S2050, _S2051, _S2052, _S2053, _S2054, _S2055, &_S2060);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  densities_10, float3  rgb_1, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_0, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    FixedArray<float, 8>  _S2061 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = densities_10;
    (&dp_densities_1)->differential_0 = _S2061;
    float3  _S2062 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S2062;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2062;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2062;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_2 = dp_rgb_0.differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

