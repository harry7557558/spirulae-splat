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

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  densities_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_0, float fy_0, float cx_0, float cy_0, FixedArray<float, 10>  dist_coeffs_0, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * depth_0, float3  * rgbs_0)
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
        *aabb_xyxy_0 = make_float4 (_S170, _S172, _S169, _S171);
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
        *rgbs_0 = max_0(make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_7) * (*sh_coeffs_0)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_0 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_7) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_7 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_0 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_7) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_7 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  densities_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_1, float fy_1, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_1, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * depth_1, float3  * rgbs_1)
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
        *aabb_xyxy_1 = make_float4 (_S369, _S371, _S368, _S370);
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
        *rgbs_1 = max_0(make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_1)[int(1)] + make_float3 (z_2) * (*sh_coeffs_1)[int(2)] - make_float3 (x_9) * (*sh_coeffs_1)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_5) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_1 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_9) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_9 * fS1_1 + y_5 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_5) * (*sh_coeffs_1)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_1 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_9) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_9 * fC1_1 - y_5 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  densities_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_2, float fy_2, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_2, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * depth_2, float3  * rgbs_2)
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
    *aabb_xyxy_2 = make_float4 ((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S391), (_S396)))), (_S401)))), (_S406)))), (_S411)))), (_S416)))), (_S421)))), (_S426))), (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S392), (_S397)))), (_S402)))), (_S407)))), (_S412)))), (_S417)))), (_S422)))), (_S427))), (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S391), (_S396)))), (_S401)))), (_S406)))), (_S411)))), (_S416)))), (_S421)))), (_S426))), (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S392), (_S397)))), (_S402)))), (_S407)))), (_S412)))), (_S417)))), (_S422)))), (_S427))));
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
    *rgbs_2 = max_0(make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_2)[int(1)] + make_float3 (z_3) * (*sh_coeffs_2)[int(2)] - make_float3 (x_10) * (*sh_coeffs_2)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_6) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_2 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_10) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_10 * fS1_2 + y_6 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_6) * (*sh_coeffs_2)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_2 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_10) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_10 * fC1_2 - y_6 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  densities_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_3, float fy_3, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_3, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * depth_3, float3  * rgbs_3)
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
    *aabb_xyxy_3 = make_float4 ((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S475))), (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S476))), (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].x), (uv_3[int(1)].x)))), (uv_3[int(2)].x)))), (uv_3[int(3)].x)))), (uv_3[int(4)].x)))), (uv_3[int(5)].x)))), (uv_3[int(6)].x)))), (_S475))), (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_3[int(0)].y), (uv_3[int(1)].y)))), (uv_3[int(2)].y)))), (uv_3[int(3)].y)))), (uv_3[int(4)].y)))), (uv_3[int(5)].y)))), (uv_3[int(6)].y)))), (_S476))));
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
    *rgbs_3 = max_0(make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_3)[int(1)] + make_float3 (z_5) * (*sh_coeffs_3)[int(2)] - make_float3 (x_12) * (*sh_coeffs_3)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_8) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_3 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_12) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_12 * fS1_3 + y_8 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_8) * (*sh_coeffs_3)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_3 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_12) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_12 * fC1_3 - y_8 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

struct s_bwd_prop_projection_voxel_eval3d_persp_differentiable_Intermediates_0
{
    FixedArray<float3 , 16>  _S482;
};

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S483, float3  _S484)
{
    return mul_0(_S483, _S484);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S485)
{
    return (F32_sqrt((_S485)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S486, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S487, float3  _S488)
{
    _d_max_vector_0(_S486, _S487, _S488);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S489, float _S490)
{
    _d_sqrt_0(_S489, _S490);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S491, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S492, float3  _S493)
{
    _d_mul_0(_S491, _S492, _S493);
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  * densities_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_4, float fy_4, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_4, uint image_width_4, uint image_height_4, float3  v_rgb_0, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  _S494 = make_float3 (0.0f);
    FixedArray<float3 , 16>  _S495 = {
        _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494, _S494
    };
    s_bwd_prop_projection_voxel_eval3d_persp_differentiable_Intermediates_0 _S496;
    (&_S496)->_S482 = _S495;
    (&_S496)->_S482 = *sh_coeffs_4;
    float3  _S497 = s_primal_ctx_mul_0(R_4, pos_4) + t_4;
    float _S498 = _S497.z;
    float2  _S499 = make_float2 (_S498);
    float _S500 = (F32_min((1.00000001504746622e+30f), (_S498)));
    float _S501 = (F32_max((0.0f), (_S498)));
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S502 = s_primal_ctx_mul_0(R_4, pos_i_0) + t_4;
    float _S503 = _S502.z;
    float2  _S504 = make_float2 (_S503);
    float _S505 = (F32_min((_S500), (_S503)));
    float _S506 = (F32_max((_S501), (_S503)));
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S507 = s_primal_ctx_mul_0(R_4, pos_i_1) + t_4;
    float _S508 = _S507.z;
    float2  _S509 = make_float2 (_S508);
    float _S510 = (F32_min((_S505), (_S508)));
    float _S511 = (F32_max((_S506), (_S508)));
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S512 = s_primal_ctx_mul_0(R_4, pos_i_2) + t_4;
    float _S513 = _S512.z;
    float2  _S514 = make_float2 (_S513);
    float _S515 = (F32_min((_S510), (_S513)));
    float _S516 = (F32_max((_S511), (_S513)));
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S517 = s_primal_ctx_mul_0(R_4, pos_i_3) + t_4;
    float _S518 = _S517.z;
    float2  _S519 = make_float2 (_S518);
    float _S520 = (F32_min((_S515), (_S518)));
    float _S521 = (F32_max((_S516), (_S518)));
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S522 = s_primal_ctx_mul_0(R_4, pos_i_4) + t_4;
    float _S523 = _S522.z;
    float2  _S524 = make_float2 (_S523);
    float _S525 = (F32_min((_S520), (_S523)));
    float _S526 = (F32_max((_S521), (_S523)));
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S527 = s_primal_ctx_mul_0(R_4, pos_i_5) + t_4;
    float _S528 = _S527.z;
    float2  _S529 = make_float2 (_S528);
    float _S530 = (F32_min((_S525), (_S528)));
    float _S531 = (F32_max((_S526), (_S528)));
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S532 = s_primal_ctx_mul_0(R_4, pos_i_6) + t_4;
    float _S533 = _S532.z;
    float2  _S534 = make_float2 (_S533);
    float3  _S535 = pos_4 + make_float3 (0.5f * size_4);
    float2  _S536 = float2 {_S497.x, _S497.y};
    float2  _S537 = _S536 / make_float2 (_S498);
    float2  _S538 = make_float2 (_S498 * _S498);
    float u_48 = _S537.x;
    float v_48 = _S537.y;
    float r2_48 = u_48 * u_48 + v_48 * v_48;
    float _S539 = dist_coeffs_4[int(2)] + r2_48 * dist_coeffs_4[int(3)];
    float _S540 = dist_coeffs_4[int(1)] + r2_48 * _S539;
    float _S541 = dist_coeffs_4[int(0)] + r2_48 * _S540;
    float radial_16 = 1.0f + r2_48 * _S541;
    float _S542 = 2.0f * dist_coeffs_4[int(4)];
    float _S543 = _S542 * u_48;
    float _S544 = 2.0f * u_48;
    float _S545 = 2.0f * dist_coeffs_4[int(5)];
    float _S546 = _S545 * u_48;
    float _S547 = 2.0f * v_48;
    float2  _S548 = _S537 * make_float2 (radial_16) + make_float2 (_S543 * v_48 + dist_coeffs_4[int(5)] * (r2_48 + _S544 * u_48) + dist_coeffs_4[int(6)] * r2_48, _S546 * v_48 + dist_coeffs_4[int(4)] * (r2_48 + _S547 * v_48) + dist_coeffs_4[int(7)] * r2_48);
    float2  _S549 = _S548 + make_float2 (dist_coeffs_4[int(8)] * _S548.x + dist_coeffs_4[int(9)] * _S548.y, 0.0f);
    float _S550 = fx_4 * _S549.x + cx_4;
    float _S551 = fy_4 * _S549.y + cy_4;
    float2  _S552 = float2 {_S502.x, _S502.y};
    float2  _S553 = _S552 / make_float2 (_S503);
    float2  _S554 = make_float2 (_S503 * _S503);
    float u_49 = _S553.x;
    float v_49 = _S553.y;
    float r2_49 = u_49 * u_49 + v_49 * v_49;
    float _S555 = dist_coeffs_4[int(2)] + r2_49 * dist_coeffs_4[int(3)];
    float _S556 = dist_coeffs_4[int(1)] + r2_49 * _S555;
    float _S557 = dist_coeffs_4[int(0)] + r2_49 * _S556;
    float radial_17 = 1.0f + r2_49 * _S557;
    float _S558 = _S542 * u_49;
    float _S559 = 2.0f * u_49;
    float _S560 = _S545 * u_49;
    float _S561 = 2.0f * v_49;
    float2  _S562 = _S553 * make_float2 (radial_17) + make_float2 (_S558 * v_49 + dist_coeffs_4[int(5)] * (r2_49 + _S559 * u_49) + dist_coeffs_4[int(6)] * r2_49, _S560 * v_49 + dist_coeffs_4[int(4)] * (r2_49 + _S561 * v_49) + dist_coeffs_4[int(7)] * r2_49);
    float2  _S563 = _S562 + make_float2 (dist_coeffs_4[int(8)] * _S562.x + dist_coeffs_4[int(9)] * _S562.y, 0.0f);
    float _S564 = fx_4 * _S563.x + cx_4;
    float _S565 = fy_4 * _S563.y + cy_4;
    float2  _S566 = float2 {_S507.x, _S507.y};
    float2  _S567 = _S566 / make_float2 (_S508);
    float2  _S568 = make_float2 (_S508 * _S508);
    float u_50 = _S567.x;
    float v_50 = _S567.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S569 = dist_coeffs_4[int(2)] + r2_50 * dist_coeffs_4[int(3)];
    float _S570 = dist_coeffs_4[int(1)] + r2_50 * _S569;
    float _S571 = dist_coeffs_4[int(0)] + r2_50 * _S570;
    float radial_18 = 1.0f + r2_50 * _S571;
    float _S572 = _S542 * u_50;
    float _S573 = 2.0f * u_50;
    float _S574 = _S545 * u_50;
    float _S575 = 2.0f * v_50;
    float2  _S576 = _S567 * make_float2 (radial_18) + make_float2 (_S572 * v_50 + dist_coeffs_4[int(5)] * (r2_50 + _S573 * u_50) + dist_coeffs_4[int(6)] * r2_50, _S574 * v_50 + dist_coeffs_4[int(4)] * (r2_50 + _S575 * v_50) + dist_coeffs_4[int(7)] * r2_50);
    float2  _S577 = _S576 + make_float2 (dist_coeffs_4[int(8)] * _S576.x + dist_coeffs_4[int(9)] * _S576.y, 0.0f);
    float _S578 = fx_4 * _S577.x + cx_4;
    float _S579 = fy_4 * _S577.y + cy_4;
    float2  _S580 = float2 {_S512.x, _S512.y};
    float2  _S581 = _S580 / make_float2 (_S513);
    float2  _S582 = make_float2 (_S513 * _S513);
    float u_51 = _S581.x;
    float v_51 = _S581.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float _S583 = dist_coeffs_4[int(2)] + r2_51 * dist_coeffs_4[int(3)];
    float _S584 = dist_coeffs_4[int(1)] + r2_51 * _S583;
    float _S585 = dist_coeffs_4[int(0)] + r2_51 * _S584;
    float radial_19 = 1.0f + r2_51 * _S585;
    float _S586 = _S542 * u_51;
    float _S587 = 2.0f * u_51;
    float _S588 = _S545 * u_51;
    float _S589 = 2.0f * v_51;
    float2  _S590 = _S581 * make_float2 (radial_19) + make_float2 (_S586 * v_51 + dist_coeffs_4[int(5)] * (r2_51 + _S587 * u_51) + dist_coeffs_4[int(6)] * r2_51, _S588 * v_51 + dist_coeffs_4[int(4)] * (r2_51 + _S589 * v_51) + dist_coeffs_4[int(7)] * r2_51);
    float2  _S591 = _S590 + make_float2 (dist_coeffs_4[int(8)] * _S590.x + dist_coeffs_4[int(9)] * _S590.y, 0.0f);
    float _S592 = fx_4 * _S591.x + cx_4;
    float _S593 = fy_4 * _S591.y + cy_4;
    float2  _S594 = float2 {_S517.x, _S517.y};
    float2  _S595 = _S594 / make_float2 (_S518);
    float2  _S596 = make_float2 (_S518 * _S518);
    float u_52 = _S595.x;
    float v_52 = _S595.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float _S597 = dist_coeffs_4[int(2)] + r2_52 * dist_coeffs_4[int(3)];
    float _S598 = dist_coeffs_4[int(1)] + r2_52 * _S597;
    float _S599 = dist_coeffs_4[int(0)] + r2_52 * _S598;
    float radial_20 = 1.0f + r2_52 * _S599;
    float _S600 = _S542 * u_52;
    float _S601 = 2.0f * u_52;
    float _S602 = _S545 * u_52;
    float _S603 = 2.0f * v_52;
    float2  _S604 = _S595 * make_float2 (radial_20) + make_float2 (_S600 * v_52 + dist_coeffs_4[int(5)] * (r2_52 + _S601 * u_52) + dist_coeffs_4[int(6)] * r2_52, _S602 * v_52 + dist_coeffs_4[int(4)] * (r2_52 + _S603 * v_52) + dist_coeffs_4[int(7)] * r2_52);
    float2  _S605 = _S604 + make_float2 (dist_coeffs_4[int(8)] * _S604.x + dist_coeffs_4[int(9)] * _S604.y, 0.0f);
    float _S606 = fx_4 * _S605.x + cx_4;
    float _S607 = fy_4 * _S605.y + cy_4;
    float2  _S608 = float2 {_S522.x, _S522.y};
    float2  _S609 = _S608 / make_float2 (_S523);
    float2  _S610 = make_float2 (_S523 * _S523);
    float u_53 = _S609.x;
    float v_53 = _S609.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S611 = dist_coeffs_4[int(2)] + r2_53 * dist_coeffs_4[int(3)];
    float _S612 = dist_coeffs_4[int(1)] + r2_53 * _S611;
    float _S613 = dist_coeffs_4[int(0)] + r2_53 * _S612;
    float radial_21 = 1.0f + r2_53 * _S613;
    float _S614 = _S542 * u_53;
    float _S615 = 2.0f * u_53;
    float _S616 = _S545 * u_53;
    float _S617 = 2.0f * v_53;
    float2  _S618 = _S609 * make_float2 (radial_21) + make_float2 (_S614 * v_53 + dist_coeffs_4[int(5)] * (r2_53 + _S615 * u_53) + dist_coeffs_4[int(6)] * r2_53, _S616 * v_53 + dist_coeffs_4[int(4)] * (r2_53 + _S617 * v_53) + dist_coeffs_4[int(7)] * r2_53);
    float2  _S619 = _S618 + make_float2 (dist_coeffs_4[int(8)] * _S618.x + dist_coeffs_4[int(9)] * _S618.y, 0.0f);
    float _S620 = fx_4 * _S619.x + cx_4;
    float _S621 = fy_4 * _S619.y + cy_4;
    float2  _S622 = float2 {_S527.x, _S527.y};
    float2  _S623 = _S622 / make_float2 (_S528);
    float2  _S624 = make_float2 (_S528 * _S528);
    float u_54 = _S623.x;
    float v_54 = _S623.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float _S625 = dist_coeffs_4[int(2)] + r2_54 * dist_coeffs_4[int(3)];
    float _S626 = dist_coeffs_4[int(1)] + r2_54 * _S625;
    float _S627 = dist_coeffs_4[int(0)] + r2_54 * _S626;
    float radial_22 = 1.0f + r2_54 * _S627;
    float _S628 = _S542 * u_54;
    float _S629 = 2.0f * u_54;
    float _S630 = _S545 * u_54;
    float _S631 = 2.0f * v_54;
    float2  _S632 = _S623 * make_float2 (radial_22) + make_float2 (_S628 * v_54 + dist_coeffs_4[int(5)] * (r2_54 + _S629 * u_54) + dist_coeffs_4[int(6)] * r2_54, _S630 * v_54 + dist_coeffs_4[int(4)] * (r2_54 + _S631 * v_54) + dist_coeffs_4[int(7)] * r2_54);
    float2  _S633 = _S632 + make_float2 (dist_coeffs_4[int(8)] * _S632.x + dist_coeffs_4[int(9)] * _S632.y, 0.0f);
    float _S634 = fx_4 * _S633.x + cx_4;
    float _S635 = fy_4 * _S633.y + cy_4;
    float2  _S636 = float2 {_S532.x, _S532.y};
    float2  _S637 = _S636 / make_float2 (_S533);
    float2  _S638 = make_float2 (_S533 * _S533);
    float u_55 = _S637.x;
    float v_55 = _S637.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float _S639 = dist_coeffs_4[int(2)] + r2_55 * dist_coeffs_4[int(3)];
    float _S640 = dist_coeffs_4[int(1)] + r2_55 * _S639;
    float _S641 = dist_coeffs_4[int(0)] + r2_55 * _S640;
    float radial_23 = 1.0f + r2_55 * _S641;
    float _S642 = _S542 * u_55;
    float _S643 = 2.0f * u_55;
    float _S644 = _S545 * u_55;
    float _S645 = 2.0f * v_55;
    float2  _S646 = _S637 * make_float2 (radial_23) + make_float2 (_S642 * v_55 + dist_coeffs_4[int(5)] * (r2_55 + _S643 * u_55) + dist_coeffs_4[int(6)] * r2_55, _S644 * v_55 + dist_coeffs_4[int(4)] * (r2_55 + _S645 * v_55) + dist_coeffs_4[int(7)] * r2_55);
    float2  _S647 = _S646 + make_float2 (dist_coeffs_4[int(8)] * _S646.x + dist_coeffs_4[int(9)] * _S646.y, 0.0f);
    float _S648 = fx_4 * _S647.x + cx_4;
    float _S649 = fy_4 * _S647.y + cy_4;
    float _S650 = (F32_max((_S550), (_S564)));
    float _S651 = (F32_min((_S550), (_S564)));
    float _S652 = (F32_max((_S551), (_S565)));
    float _S653 = (F32_min((_S551), (_S565)));
    float _S654 = (F32_max((_S650), (_S578)));
    float _S655 = (F32_min((_S651), (_S578)));
    float _S656 = (F32_max((_S652), (_S579)));
    float _S657 = (F32_min((_S653), (_S579)));
    float _S658 = (F32_max((_S654), (_S592)));
    float _S659 = (F32_min((_S655), (_S592)));
    float _S660 = (F32_max((_S656), (_S593)));
    float _S661 = (F32_min((_S657), (_S593)));
    float _S662 = (F32_max((_S658), (_S606)));
    float _S663 = (F32_min((_S659), (_S606)));
    float _S664 = (F32_max((_S660), (_S607)));
    float _S665 = (F32_min((_S661), (_S607)));
    float _S666 = (F32_max((_S662), (_S620)));
    float _S667 = (F32_min((_S663), (_S620)));
    float _S668 = (F32_max((_S664), (_S621)));
    float _S669 = (F32_min((_S665), (_S621)));
    float _S670 = (F32_max((_S666), (_S634)));
    float _S671 = (F32_min((_S667), (_S634)));
    float _S672 = (F32_max((_S668), (_S635)));
    float _S673 = (F32_min((_S669), (_S635)));
    Matrix<float, 3, 3>  _S674 = transpose_1(R_4);
    float3  _S675 = s_primal_ctx_mul_0(R_4, _S535) + t_4 - - s_primal_ctx_mul_0(_S674, t_4);
    float _S676 = _S675.x;
    float _S677 = _S675.y;
    float _S678 = _S675.z;
    float _S679 = _S676 * _S676 + _S677 * _S677 + _S678 * _S678;
    float _S680 = s_primal_ctx_sqrt_0(_S679);
    float x_13 = _S676 / _S680;
    float3  _S681 = make_float3 (x_13);
    float _S682 = _S680 * _S680;
    float y_9 = _S677 / _S680;
    float z_6 = _S678 / _S680;
    float3  _S683 = make_float3 (z_6);
    float _S684 = - y_9;
    float3  _S685 = make_float3 (_S684);
    float z2_4 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_4 = x_13 * x_13 - y_9 * y_9;
    float _S686 = 2.0f * x_13;
    float fS1_4 = _S686 * y_9;
    float pSH6_0 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float3  _S687 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_4 * x_13;
    float3  _S688 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_4 * y_9;
    float3  _S689 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_4;
    float3  _S690 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_4;
    float3  _S691 = make_float3 (pSH4_0);
    float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_6;
    float _S692 = 1.86588168144226074f * z2_4 - 1.11952900886535645f;
    float pSH12_0 = z_6 * _S692;
    float3  _S693 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_4 * x_13;
    float3  _S694 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_4 * y_9;
    float3  _S695 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_4 * fC1_4;
    float3  _S696 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_4 * fS1_4;
    float3  _S697 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_13 * fC1_4 - y_9 * fS1_4);
    float3  _S698 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_13 * fS1_4 + y_9 * fC1_4);
    float3  _S699 = make_float3 (pSH9_0);
    float3  _S700 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S701;
    (&_S701)->primal_0 = make_float3 (0.282094806432724f) * (&_S496)->_S482[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S684) * (&_S496)->_S482[int(1)] + make_float3 (z_6) * (&_S496)->_S482[int(2)] - make_float3 (x_13) * (&_S496)->_S482[int(3)]) + (make_float3 (pSH4_0) * (&_S496)->_S482[int(4)] + make_float3 (pSH5_0) * (&_S496)->_S482[int(5)] + make_float3 (pSH6_0) * (&_S496)->_S482[int(6)] + make_float3 (pSH7_0) * (&_S496)->_S482[int(7)] + make_float3 (pSH8_0) * (&_S496)->_S482[int(8)]) + (make_float3 (pSH9_0) * (&_S496)->_S482[int(9)] + make_float3 (pSH10_0) * (&_S496)->_S482[int(10)] + make_float3 (pSH11_0) * (&_S496)->_S482[int(11)] + make_float3 (pSH12_0) * (&_S496)->_S482[int(12)] + make_float3 (pSH13_0) * (&_S496)->_S482[int(13)] + make_float3 (pSH14_0) * (&_S496)->_S482[int(14)] + make_float3 (pSH15_0) * (&_S496)->_S482[int(15)]) + make_float3 (0.5f);
    (&_S701)->differential_0 = _S494;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S702;
    (&_S702)->primal_0 = _S700;
    (&_S702)->differential_0 = _S494;
    s_bwd_prop_max_0(&_S701, &_S702, v_rgb_0);
    float3  _S703 = _S698 * _S701.differential_0;
    float3  _S704 = (&_S496)->_S482[int(15)] * _S701.differential_0;
    float3  _S705 = _S696 * _S701.differential_0;
    float3  _S706 = (&_S496)->_S482[int(14)] * _S701.differential_0;
    float3  _S707 = _S694 * _S701.differential_0;
    float3  _S708 = (&_S496)->_S482[int(13)] * _S701.differential_0;
    float3  _S709 = _S693 * _S701.differential_0;
    float3  _S710 = (&_S496)->_S482[int(12)] * _S701.differential_0;
    float3  _S711 = _S695 * _S701.differential_0;
    float3  _S712 = (&_S496)->_S482[int(11)] * _S701.differential_0;
    float3  _S713 = _S697 * _S701.differential_0;
    float3  _S714 = (&_S496)->_S482[int(10)] * _S701.differential_0;
    float3  _S715 = _S699 * _S701.differential_0;
    float3  _S716 = (&_S496)->_S482[int(9)] * _S701.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S716.x + _S716.y + _S716.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S704.x + _S704.y + _S704.z);
    float _S717 = _S714.x + _S714.y + _S714.z;
    float _S718 = _S706.x + _S706.y + _S706.z;
    float _S719 = _S712.x + _S712.y + _S712.z;
    float _S720 = _S708.x + _S708.y + _S708.z;
    float _S721 = _S710.x + _S710.y + _S710.z;
    float _S722 = - s_diff_fC2_T_0;
    float3  _S723 = _S690 * _S701.differential_0;
    float3  _S724 = (&_S496)->_S482[int(8)] * _S701.differential_0;
    float3  _S725 = _S688 * _S701.differential_0;
    float3  _S726 = (&_S496)->_S482[int(7)] * _S701.differential_0;
    float3  _S727 = _S687 * _S701.differential_0;
    float3  _S728 = (&_S496)->_S482[int(6)] * _S701.differential_0;
    float3  _S729 = _S689 * _S701.differential_0;
    float3  _S730 = (&_S496)->_S482[int(5)] * _S701.differential_0;
    float3  _S731 = _S691 * _S701.differential_0;
    float3  _S732 = (&_S496)->_S482[int(4)] * _S701.differential_0;
    float _S733 = _S730.x + _S730.y + _S730.z;
    float _S734 = _S726.x + _S726.y + _S726.z;
    float _S735 = fTmp1B_4 * _S717 + x_13 * s_diff_fS2_T_0 + y_9 * _S722 + 0.54627424478530884f * (_S732.x + _S732.y + _S732.z);
    float _S736 = fTmp1B_4 * _S718 + y_9 * s_diff_fS2_T_0 + x_13 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S724.x + _S724.y + _S724.z);
    float _S737 = y_9 * - _S736;
    float _S738 = x_13 * _S736;
    float _S739 = z_6 * (1.86588168144226074f * (z_6 * _S721) + -2.28522896766662598f * (y_9 * _S719 + x_13 * _S720) + 0.94617468118667603f * (_S728.x + _S728.y + _S728.z));
    float3  _S740 = make_float3 (0.48860251903533936f) * _S701.differential_0;
    float3  _S741 = - _S740;
    float3  _S742 = _S681 * _S741;
    float3  _S743 = (&_S496)->_S482[int(3)] * _S741;
    float3  _S744 = _S683 * _S740;
    float3  _S745 = (&_S496)->_S482[int(2)] * _S740;
    float3  _S746 = _S685 * _S740;
    float3  _S747 = (&_S496)->_S482[int(1)] * _S740;
    float _S748 = (_S692 * _S721 + 1.44530570507049561f * (fS1_4 * _S717 + fC1_4 * _S718) + -1.09254848957061768f * (y_9 * _S733 + x_13 * _S734) + _S739 + _S739 + _S745.x + _S745.y + _S745.z) / _S682;
    float _S749 = _S680 * _S748;
    float _S750 = (fTmp0C_4 * _S719 + fC1_4 * s_diff_fS2_T_0 + fS1_4 * _S722 + fTmp0B_4 * _S733 + _S686 * _S735 + _S737 + _S737 + - (_S747.x + _S747.y + _S747.z)) / _S682;
    float _S751 = _S680 * _S750;
    float _S752 = (fTmp0C_4 * _S720 + fS1_4 * s_diff_fS2_T_0 + fC1_4 * s_diff_fC2_T_0 + fTmp0B_4 * _S734 + 2.0f * (y_9 * _S735) + _S738 + _S738 + _S743.x + _S743.y + _S743.z) / _S682;
    float _S753 = _S680 * _S752;
    float _S754 = _S678 * - _S748 + _S677 * - _S750 + _S676 * - _S752;
    DiffPair_float_0 _S755;
    (&_S755)->primal_0 = _S679;
    (&_S755)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S755, _S754);
    float _S756 = _S678 * _S755.differential_0;
    float _S757 = _S677 * _S755.differential_0;
    float _S758 = _S676 * _S755.differential_0;
    float3  _S759 = make_float3 (0.282094806432724f) * _S701.differential_0;
    float3  _S760 = make_float3 (_S753 + _S758 + _S758, _S751 + _S757 + _S757, _S749 + _S756 + _S756);
    float3  _S761 = - - _S760;
    Matrix<float, 3, 3>  _S762 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S763;
    (&_S763)->primal_0 = _S674;
    (&_S763)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S764;
    (&_S764)->primal_0 = t_4;
    (&_S764)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S763, &_S764, _S761);
    Matrix<float, 3, 3>  _S765 = transpose_1(_S763.differential_0);
    DiffPair_float_0 _S766;
    (&_S766)->primal_0 = _S673;
    (&_S766)->differential_0 = 0.0f;
    DiffPair_float_0 _S767;
    (&_S767)->primal_0 = _S649;
    (&_S767)->differential_0 = 0.0f;
    _d_min_0(&_S766, &_S767, 0.0f);
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = _S672;
    (&_S768)->differential_0 = 0.0f;
    DiffPair_float_0 _S769;
    (&_S769)->primal_0 = _S649;
    (&_S769)->differential_0 = 0.0f;
    _d_max_0(&_S768, &_S769, 0.0f);
    DiffPair_float_0 _S770;
    (&_S770)->primal_0 = _S671;
    (&_S770)->differential_0 = 0.0f;
    DiffPair_float_0 _S771;
    (&_S771)->primal_0 = _S648;
    (&_S771)->differential_0 = 0.0f;
    _d_min_0(&_S770, &_S771, 0.0f);
    DiffPair_float_0 _S772;
    (&_S772)->primal_0 = _S670;
    (&_S772)->differential_0 = 0.0f;
    DiffPair_float_0 _S773;
    (&_S773)->primal_0 = _S648;
    (&_S773)->differential_0 = 0.0f;
    _d_max_0(&_S772, &_S773, 0.0f);
    DiffPair_float_0 _S774;
    (&_S774)->primal_0 = _S669;
    (&_S774)->differential_0 = 0.0f;
    DiffPair_float_0 _S775;
    (&_S775)->primal_0 = _S635;
    (&_S775)->differential_0 = 0.0f;
    _d_min_0(&_S774, &_S775, _S766.differential_0);
    DiffPair_float_0 _S776;
    (&_S776)->primal_0 = _S668;
    (&_S776)->differential_0 = 0.0f;
    DiffPair_float_0 _S777;
    (&_S777)->primal_0 = _S635;
    (&_S777)->differential_0 = 0.0f;
    _d_max_0(&_S776, &_S777, _S768.differential_0);
    DiffPair_float_0 _S778;
    (&_S778)->primal_0 = _S667;
    (&_S778)->differential_0 = 0.0f;
    DiffPair_float_0 _S779;
    (&_S779)->primal_0 = _S634;
    (&_S779)->differential_0 = 0.0f;
    _d_min_0(&_S778, &_S779, _S770.differential_0);
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = _S666;
    (&_S780)->differential_0 = 0.0f;
    DiffPair_float_0 _S781;
    (&_S781)->primal_0 = _S634;
    (&_S781)->differential_0 = 0.0f;
    _d_max_0(&_S780, &_S781, _S772.differential_0);
    DiffPair_float_0 _S782;
    (&_S782)->primal_0 = _S665;
    (&_S782)->differential_0 = 0.0f;
    DiffPair_float_0 _S783;
    (&_S783)->primal_0 = _S621;
    (&_S783)->differential_0 = 0.0f;
    _d_min_0(&_S782, &_S783, _S774.differential_0);
    DiffPair_float_0 _S784;
    (&_S784)->primal_0 = _S664;
    (&_S784)->differential_0 = 0.0f;
    DiffPair_float_0 _S785;
    (&_S785)->primal_0 = _S621;
    (&_S785)->differential_0 = 0.0f;
    _d_max_0(&_S784, &_S785, _S776.differential_0);
    DiffPair_float_0 _S786;
    (&_S786)->primal_0 = _S663;
    (&_S786)->differential_0 = 0.0f;
    DiffPair_float_0 _S787;
    (&_S787)->primal_0 = _S620;
    (&_S787)->differential_0 = 0.0f;
    _d_min_0(&_S786, &_S787, _S778.differential_0);
    DiffPair_float_0 _S788;
    (&_S788)->primal_0 = _S662;
    (&_S788)->differential_0 = 0.0f;
    DiffPair_float_0 _S789;
    (&_S789)->primal_0 = _S620;
    (&_S789)->differential_0 = 0.0f;
    _d_max_0(&_S788, &_S789, _S780.differential_0);
    DiffPair_float_0 _S790;
    (&_S790)->primal_0 = _S661;
    (&_S790)->differential_0 = 0.0f;
    DiffPair_float_0 _S791;
    (&_S791)->primal_0 = _S607;
    (&_S791)->differential_0 = 0.0f;
    _d_min_0(&_S790, &_S791, _S782.differential_0);
    DiffPair_float_0 _S792;
    (&_S792)->primal_0 = _S660;
    (&_S792)->differential_0 = 0.0f;
    DiffPair_float_0 _S793;
    (&_S793)->primal_0 = _S607;
    (&_S793)->differential_0 = 0.0f;
    _d_max_0(&_S792, &_S793, _S784.differential_0);
    DiffPair_float_0 _S794;
    (&_S794)->primal_0 = _S659;
    (&_S794)->differential_0 = 0.0f;
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = _S606;
    (&_S795)->differential_0 = 0.0f;
    _d_min_0(&_S794, &_S795, _S786.differential_0);
    DiffPair_float_0 _S796;
    (&_S796)->primal_0 = _S658;
    (&_S796)->differential_0 = 0.0f;
    DiffPair_float_0 _S797;
    (&_S797)->primal_0 = _S606;
    (&_S797)->differential_0 = 0.0f;
    _d_max_0(&_S796, &_S797, _S788.differential_0);
    DiffPair_float_0 _S798;
    (&_S798)->primal_0 = _S657;
    (&_S798)->differential_0 = 0.0f;
    DiffPair_float_0 _S799;
    (&_S799)->primal_0 = _S593;
    (&_S799)->differential_0 = 0.0f;
    _d_min_0(&_S798, &_S799, _S790.differential_0);
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = _S656;
    (&_S800)->differential_0 = 0.0f;
    DiffPair_float_0 _S801;
    (&_S801)->primal_0 = _S593;
    (&_S801)->differential_0 = 0.0f;
    _d_max_0(&_S800, &_S801, _S792.differential_0);
    DiffPair_float_0 _S802;
    (&_S802)->primal_0 = _S655;
    (&_S802)->differential_0 = 0.0f;
    DiffPair_float_0 _S803;
    (&_S803)->primal_0 = _S592;
    (&_S803)->differential_0 = 0.0f;
    _d_min_0(&_S802, &_S803, _S794.differential_0);
    DiffPair_float_0 _S804;
    (&_S804)->primal_0 = _S654;
    (&_S804)->differential_0 = 0.0f;
    DiffPair_float_0 _S805;
    (&_S805)->primal_0 = _S592;
    (&_S805)->differential_0 = 0.0f;
    _d_max_0(&_S804, &_S805, _S796.differential_0);
    DiffPair_float_0 _S806;
    (&_S806)->primal_0 = _S653;
    (&_S806)->differential_0 = 0.0f;
    DiffPair_float_0 _S807;
    (&_S807)->primal_0 = _S579;
    (&_S807)->differential_0 = 0.0f;
    _d_min_0(&_S806, &_S807, _S798.differential_0);
    DiffPair_float_0 _S808;
    (&_S808)->primal_0 = _S652;
    (&_S808)->differential_0 = 0.0f;
    DiffPair_float_0 _S809;
    (&_S809)->primal_0 = _S579;
    (&_S809)->differential_0 = 0.0f;
    _d_max_0(&_S808, &_S809, _S800.differential_0);
    DiffPair_float_0 _S810;
    (&_S810)->primal_0 = _S651;
    (&_S810)->differential_0 = 0.0f;
    DiffPair_float_0 _S811;
    (&_S811)->primal_0 = _S578;
    (&_S811)->differential_0 = 0.0f;
    _d_min_0(&_S810, &_S811, _S802.differential_0);
    DiffPair_float_0 _S812;
    (&_S812)->primal_0 = _S650;
    (&_S812)->differential_0 = 0.0f;
    DiffPair_float_0 _S813;
    (&_S813)->primal_0 = _S578;
    (&_S813)->differential_0 = 0.0f;
    _d_max_0(&_S812, &_S813, _S804.differential_0);
    DiffPair_float_0 _S814;
    (&_S814)->primal_0 = _S551;
    (&_S814)->differential_0 = 0.0f;
    DiffPair_float_0 _S815;
    (&_S815)->primal_0 = _S565;
    (&_S815)->differential_0 = 0.0f;
    _d_min_0(&_S814, &_S815, _S806.differential_0);
    DiffPair_float_0 _S816;
    (&_S816)->primal_0 = _S551;
    (&_S816)->differential_0 = 0.0f;
    DiffPair_float_0 _S817;
    (&_S817)->primal_0 = _S565;
    (&_S817)->differential_0 = 0.0f;
    _d_max_0(&_S816, &_S817, _S808.differential_0);
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = _S550;
    (&_S818)->differential_0 = 0.0f;
    DiffPair_float_0 _S819;
    (&_S819)->primal_0 = _S564;
    (&_S819)->differential_0 = 0.0f;
    _d_min_0(&_S818, &_S819, _S810.differential_0);
    DiffPair_float_0 _S820;
    (&_S820)->primal_0 = _S550;
    (&_S820)->differential_0 = 0.0f;
    DiffPair_float_0 _S821;
    (&_S821)->primal_0 = _S564;
    (&_S821)->differential_0 = 0.0f;
    _d_max_0(&_S820, &_S821, _S812.differential_0);
    float _S822 = fx_4 * (_S771.differential_0 + _S773.differential_0);
    float2  _S823 = make_float2 (_S822, fy_4 * (_S767.differential_0 + _S769.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S822, dist_coeffs_4[int(9)] * _S822);
    float2  _S824 = _S637 * _S823;
    float _S825 = dist_coeffs_4[int(4)] * _S823.y;
    float _S826 = dist_coeffs_4[int(5)] * _S823.x;
    float _S827 = _S824.x + _S824.y;
    float _S828 = r2_55 * _S827;
    float _S829 = r2_55 * _S828;
    float _S830 = dist_coeffs_4[int(7)] * _S823.y + _S825 + dist_coeffs_4[int(6)] * _S823.x + _S826 + _S641 * _S827 + _S640 * _S828 + _S639 * _S829 + dist_coeffs_4[int(3)] * (r2_55 * _S829);
    float _S831 = v_55 * _S830;
    float _S832 = u_55 * _S830;
    float2  _S833 = (make_float2 (radial_23) * _S823 + make_float2 (_S545 * (v_55 * _S823.y) + _S643 * _S826 + 2.0f * (u_55 * _S826) + _S542 * (v_55 * _S823.x) + _S832 + _S832, _S645 * _S825 + 2.0f * (v_55 * _S825) + _S644 * _S823.y + _S642 * _S823.x + _S831 + _S831)) / _S638;
    float2  _S834 = _S636 * - _S833;
    float2  _S835 = _S534 * _S833;
    float _S836 = fx_4 * (_S779.differential_0 + _S781.differential_0);
    float2  _S837 = make_float2 (_S836, fy_4 * (_S775.differential_0 + _S777.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S836, dist_coeffs_4[int(9)] * _S836);
    float2  _S838 = _S623 * _S837;
    float _S839 = dist_coeffs_4[int(4)] * _S837.y;
    float _S840 = dist_coeffs_4[int(5)] * _S837.x;
    float _S841 = _S838.x + _S838.y;
    float _S842 = r2_54 * _S841;
    float _S843 = r2_54 * _S842;
    float _S844 = dist_coeffs_4[int(7)] * _S837.y + _S839 + dist_coeffs_4[int(6)] * _S837.x + _S840 + _S627 * _S841 + _S626 * _S842 + _S625 * _S843 + dist_coeffs_4[int(3)] * (r2_54 * _S843);
    float _S845 = v_54 * _S844;
    float _S846 = u_54 * _S844;
    float2  _S847 = (make_float2 (radial_22) * _S837 + make_float2 (_S545 * (v_54 * _S837.y) + _S629 * _S840 + 2.0f * (u_54 * _S840) + _S542 * (v_54 * _S837.x) + _S846 + _S846, _S631 * _S839 + 2.0f * (v_54 * _S839) + _S630 * _S837.y + _S628 * _S837.x + _S845 + _S845)) / _S624;
    float2  _S848 = _S622 * - _S847;
    float2  _S849 = _S529 * _S847;
    float _S850 = fx_4 * (_S787.differential_0 + _S789.differential_0);
    float2  _S851 = make_float2 (_S850, fy_4 * (_S783.differential_0 + _S785.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S850, dist_coeffs_4[int(9)] * _S850);
    float2  _S852 = _S609 * _S851;
    float _S853 = dist_coeffs_4[int(4)] * _S851.y;
    float _S854 = dist_coeffs_4[int(5)] * _S851.x;
    float _S855 = _S852.x + _S852.y;
    float _S856 = r2_53 * _S855;
    float _S857 = r2_53 * _S856;
    float _S858 = dist_coeffs_4[int(7)] * _S851.y + _S853 + dist_coeffs_4[int(6)] * _S851.x + _S854 + _S613 * _S855 + _S612 * _S856 + _S611 * _S857 + dist_coeffs_4[int(3)] * (r2_53 * _S857);
    float _S859 = v_53 * _S858;
    float _S860 = u_53 * _S858;
    float2  _S861 = (make_float2 (radial_21) * _S851 + make_float2 (_S545 * (v_53 * _S851.y) + _S615 * _S854 + 2.0f * (u_53 * _S854) + _S542 * (v_53 * _S851.x) + _S860 + _S860, _S617 * _S853 + 2.0f * (v_53 * _S853) + _S616 * _S851.y + _S614 * _S851.x + _S859 + _S859)) / _S610;
    float2  _S862 = _S608 * - _S861;
    float2  _S863 = _S524 * _S861;
    float _S864 = fx_4 * (_S795.differential_0 + _S797.differential_0);
    float2  _S865 = make_float2 (_S864, fy_4 * (_S791.differential_0 + _S793.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S864, dist_coeffs_4[int(9)] * _S864);
    float2  _S866 = _S595 * _S865;
    float _S867 = dist_coeffs_4[int(4)] * _S865.y;
    float _S868 = dist_coeffs_4[int(5)] * _S865.x;
    float _S869 = _S866.x + _S866.y;
    float _S870 = r2_52 * _S869;
    float _S871 = r2_52 * _S870;
    float _S872 = dist_coeffs_4[int(7)] * _S865.y + _S867 + dist_coeffs_4[int(6)] * _S865.x + _S868 + _S599 * _S869 + _S598 * _S870 + _S597 * _S871 + dist_coeffs_4[int(3)] * (r2_52 * _S871);
    float _S873 = v_52 * _S872;
    float _S874 = u_52 * _S872;
    float2  _S875 = (make_float2 (radial_20) * _S865 + make_float2 (_S545 * (v_52 * _S865.y) + _S601 * _S868 + 2.0f * (u_52 * _S868) + _S542 * (v_52 * _S865.x) + _S874 + _S874, _S603 * _S867 + 2.0f * (v_52 * _S867) + _S602 * _S865.y + _S600 * _S865.x + _S873 + _S873)) / _S596;
    float2  _S876 = _S594 * - _S875;
    float2  _S877 = _S519 * _S875;
    float _S878 = fx_4 * (_S803.differential_0 + _S805.differential_0);
    float2  _S879 = make_float2 (_S878, fy_4 * (_S799.differential_0 + _S801.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S878, dist_coeffs_4[int(9)] * _S878);
    float2  _S880 = _S581 * _S879;
    float _S881 = dist_coeffs_4[int(4)] * _S879.y;
    float _S882 = dist_coeffs_4[int(5)] * _S879.x;
    float _S883 = _S880.x + _S880.y;
    float _S884 = r2_51 * _S883;
    float _S885 = r2_51 * _S884;
    float _S886 = dist_coeffs_4[int(7)] * _S879.y + _S881 + dist_coeffs_4[int(6)] * _S879.x + _S882 + _S585 * _S883 + _S584 * _S884 + _S583 * _S885 + dist_coeffs_4[int(3)] * (r2_51 * _S885);
    float _S887 = v_51 * _S886;
    float _S888 = u_51 * _S886;
    float2  _S889 = (make_float2 (radial_19) * _S879 + make_float2 (_S545 * (v_51 * _S879.y) + _S587 * _S882 + 2.0f * (u_51 * _S882) + _S542 * (v_51 * _S879.x) + _S888 + _S888, _S589 * _S881 + 2.0f * (v_51 * _S881) + _S588 * _S879.y + _S586 * _S879.x + _S887 + _S887)) / _S582;
    float2  _S890 = _S580 * - _S889;
    float2  _S891 = _S514 * _S889;
    float _S892 = fx_4 * (_S811.differential_0 + _S813.differential_0);
    float2  _S893 = make_float2 (_S892, fy_4 * (_S807.differential_0 + _S809.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S892, dist_coeffs_4[int(9)] * _S892);
    float2  _S894 = _S567 * _S893;
    float _S895 = dist_coeffs_4[int(4)] * _S893.y;
    float _S896 = dist_coeffs_4[int(5)] * _S893.x;
    float _S897 = _S894.x + _S894.y;
    float _S898 = r2_50 * _S897;
    float _S899 = r2_50 * _S898;
    float _S900 = dist_coeffs_4[int(7)] * _S893.y + _S895 + dist_coeffs_4[int(6)] * _S893.x + _S896 + _S571 * _S897 + _S570 * _S898 + _S569 * _S899 + dist_coeffs_4[int(3)] * (r2_50 * _S899);
    float _S901 = v_50 * _S900;
    float _S902 = u_50 * _S900;
    float2  _S903 = (make_float2 (radial_18) * _S893 + make_float2 (_S545 * (v_50 * _S893.y) + _S573 * _S896 + 2.0f * (u_50 * _S896) + _S542 * (v_50 * _S893.x) + _S902 + _S902, _S575 * _S895 + 2.0f * (v_50 * _S895) + _S574 * _S893.y + _S572 * _S893.x + _S901 + _S901)) / _S568;
    float2  _S904 = _S566 * - _S903;
    float2  _S905 = _S509 * _S903;
    float _S906 = fx_4 * (_S819.differential_0 + _S821.differential_0);
    float2  _S907 = make_float2 (_S906, fy_4 * (_S815.differential_0 + _S817.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S906, dist_coeffs_4[int(9)] * _S906);
    float2  _S908 = _S553 * _S907;
    float _S909 = dist_coeffs_4[int(4)] * _S907.y;
    float _S910 = dist_coeffs_4[int(5)] * _S907.x;
    float _S911 = _S908.x + _S908.y;
    float _S912 = r2_49 * _S911;
    float _S913 = r2_49 * _S912;
    float _S914 = dist_coeffs_4[int(7)] * _S907.y + _S909 + dist_coeffs_4[int(6)] * _S907.x + _S910 + _S557 * _S911 + _S556 * _S912 + _S555 * _S913 + dist_coeffs_4[int(3)] * (r2_49 * _S913);
    float _S915 = v_49 * _S914;
    float _S916 = u_49 * _S914;
    float2  _S917 = (make_float2 (radial_17) * _S907 + make_float2 (_S545 * (v_49 * _S907.y) + _S559 * _S910 + 2.0f * (u_49 * _S910) + _S542 * (v_49 * _S907.x) + _S916 + _S916, _S561 * _S909 + 2.0f * (v_49 * _S909) + _S560 * _S907.y + _S558 * _S907.x + _S915 + _S915)) / _S554;
    float2  _S918 = _S552 * - _S917;
    float2  _S919 = _S504 * _S917;
    float _S920 = fx_4 * (_S818.differential_0 + _S820.differential_0);
    float2  _S921 = make_float2 (_S920, fy_4 * (_S814.differential_0 + _S816.differential_0)) + make_float2 (dist_coeffs_4[int(8)] * _S920, dist_coeffs_4[int(9)] * _S920);
    float2  _S922 = _S537 * _S921;
    float _S923 = dist_coeffs_4[int(4)] * _S921.y;
    float _S924 = dist_coeffs_4[int(5)] * _S921.x;
    float _S925 = _S922.x + _S922.y;
    float _S926 = r2_48 * _S925;
    float _S927 = r2_48 * _S926;
    float _S928 = dist_coeffs_4[int(7)] * _S921.y + _S923 + dist_coeffs_4[int(6)] * _S921.x + _S924 + _S541 * _S925 + _S540 * _S926 + _S539 * _S927 + dist_coeffs_4[int(3)] * (r2_48 * _S927);
    float _S929 = v_48 * _S928;
    float _S930 = u_48 * _S928;
    float2  _S931 = (make_float2 (radial_16) * _S921 + make_float2 (_S545 * (v_48 * _S921.y) + _S544 * _S924 + 2.0f * (u_48 * _S924) + _S542 * (v_48 * _S921.x) + _S930 + _S930, _S547 * _S923 + 2.0f * (v_48 * _S923) + _S546 * _S921.y + _S543 * _S921.x + _S929 + _S929)) / _S538;
    float2  _S932 = _S536 * - _S931;
    float2  _S933 = _S499 * _S931;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S934;
    (&_S934)->primal_0 = R_4;
    (&_S934)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S935;
    (&_S935)->primal_0 = _S535;
    (&_S935)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S934, &_S935, _S760);
    DiffPair_float_0 _S936;
    (&_S936)->primal_0 = _S531;
    (&_S936)->differential_0 = 0.0f;
    DiffPair_float_0 _S937;
    (&_S937)->primal_0 = _S533;
    (&_S937)->differential_0 = 0.0f;
    _d_max_0(&_S936, &_S937, 0.0f);
    DiffPair_float_0 _S938;
    (&_S938)->primal_0 = _S530;
    (&_S938)->differential_0 = 0.0f;
    DiffPair_float_0 _S939;
    (&_S939)->primal_0 = _S533;
    (&_S939)->differential_0 = 0.0f;
    _d_min_0(&_S938, &_S939, 0.0f);
    float3  _S940 = make_float3 (_S835.x, _S835.y, _S937.differential_0 + _S939.differential_0 + _S834.x + _S834.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S941;
    (&_S941)->primal_0 = R_4;
    (&_S941)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S942;
    (&_S942)->primal_0 = pos_i_6;
    (&_S942)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S941, &_S942, _S940);
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S526;
    (&_S943)->differential_0 = 0.0f;
    DiffPair_float_0 _S944;
    (&_S944)->primal_0 = _S528;
    (&_S944)->differential_0 = 0.0f;
    _d_max_0(&_S943, &_S944, _S936.differential_0);
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S525;
    (&_S945)->differential_0 = 0.0f;
    DiffPair_float_0 _S946;
    (&_S946)->primal_0 = _S528;
    (&_S946)->differential_0 = 0.0f;
    _d_min_0(&_S945, &_S946, _S938.differential_0);
    float3  _S947 = make_float3 (_S849.x, _S849.y, _S944.differential_0 + _S946.differential_0 + _S848.x + _S848.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S948;
    (&_S948)->primal_0 = R_4;
    (&_S948)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S949;
    (&_S949)->primal_0 = pos_i_5;
    (&_S949)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S948, &_S949, _S947);
    DiffPair_float_0 _S950;
    (&_S950)->primal_0 = _S521;
    (&_S950)->differential_0 = 0.0f;
    DiffPair_float_0 _S951;
    (&_S951)->primal_0 = _S523;
    (&_S951)->differential_0 = 0.0f;
    _d_max_0(&_S950, &_S951, _S943.differential_0);
    DiffPair_float_0 _S952;
    (&_S952)->primal_0 = _S520;
    (&_S952)->differential_0 = 0.0f;
    DiffPair_float_0 _S953;
    (&_S953)->primal_0 = _S523;
    (&_S953)->differential_0 = 0.0f;
    _d_min_0(&_S952, &_S953, _S945.differential_0);
    float3  _S954 = make_float3 (_S863.x, _S863.y, _S951.differential_0 + _S953.differential_0 + _S862.x + _S862.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S955;
    (&_S955)->primal_0 = R_4;
    (&_S955)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S956;
    (&_S956)->primal_0 = pos_i_4;
    (&_S956)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S955, &_S956, _S954);
    DiffPair_float_0 _S957;
    (&_S957)->primal_0 = _S516;
    (&_S957)->differential_0 = 0.0f;
    DiffPair_float_0 _S958;
    (&_S958)->primal_0 = _S518;
    (&_S958)->differential_0 = 0.0f;
    _d_max_0(&_S957, &_S958, _S950.differential_0);
    DiffPair_float_0 _S959;
    (&_S959)->primal_0 = _S515;
    (&_S959)->differential_0 = 0.0f;
    DiffPair_float_0 _S960;
    (&_S960)->primal_0 = _S518;
    (&_S960)->differential_0 = 0.0f;
    _d_min_0(&_S959, &_S960, _S952.differential_0);
    float3  _S961 = make_float3 (_S877.x, _S877.y, _S958.differential_0 + _S960.differential_0 + _S876.x + _S876.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S962;
    (&_S962)->primal_0 = R_4;
    (&_S962)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S963;
    (&_S963)->primal_0 = pos_i_3;
    (&_S963)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S962, &_S963, _S961);
    DiffPair_float_0 _S964;
    (&_S964)->primal_0 = _S511;
    (&_S964)->differential_0 = 0.0f;
    DiffPair_float_0 _S965;
    (&_S965)->primal_0 = _S513;
    (&_S965)->differential_0 = 0.0f;
    _d_max_0(&_S964, &_S965, _S957.differential_0);
    DiffPair_float_0 _S966;
    (&_S966)->primal_0 = _S510;
    (&_S966)->differential_0 = 0.0f;
    DiffPair_float_0 _S967;
    (&_S967)->primal_0 = _S513;
    (&_S967)->differential_0 = 0.0f;
    _d_min_0(&_S966, &_S967, _S959.differential_0);
    float3  _S968 = make_float3 (_S891.x, _S891.y, _S965.differential_0 + _S967.differential_0 + _S890.x + _S890.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S969;
    (&_S969)->primal_0 = R_4;
    (&_S969)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S970;
    (&_S970)->primal_0 = pos_i_2;
    (&_S970)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S969, &_S970, _S968);
    DiffPair_float_0 _S971;
    (&_S971)->primal_0 = _S506;
    (&_S971)->differential_0 = 0.0f;
    DiffPair_float_0 _S972;
    (&_S972)->primal_0 = _S508;
    (&_S972)->differential_0 = 0.0f;
    _d_max_0(&_S971, &_S972, _S964.differential_0);
    DiffPair_float_0 _S973;
    (&_S973)->primal_0 = _S505;
    (&_S973)->differential_0 = 0.0f;
    DiffPair_float_0 _S974;
    (&_S974)->primal_0 = _S508;
    (&_S974)->differential_0 = 0.0f;
    _d_min_0(&_S973, &_S974, _S966.differential_0);
    float3  _S975 = make_float3 (_S905.x, _S905.y, _S972.differential_0 + _S974.differential_0 + _S904.x + _S904.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S976;
    (&_S976)->primal_0 = R_4;
    (&_S976)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S977;
    (&_S977)->primal_0 = pos_i_1;
    (&_S977)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S976, &_S977, _S975);
    DiffPair_float_0 _S978;
    (&_S978)->primal_0 = _S501;
    (&_S978)->differential_0 = 0.0f;
    DiffPair_float_0 _S979;
    (&_S979)->primal_0 = _S503;
    (&_S979)->differential_0 = 0.0f;
    _d_max_0(&_S978, &_S979, _S971.differential_0);
    DiffPair_float_0 _S980;
    (&_S980)->primal_0 = _S500;
    (&_S980)->differential_0 = 0.0f;
    DiffPair_float_0 _S981;
    (&_S981)->primal_0 = _S503;
    (&_S981)->differential_0 = 0.0f;
    _d_min_0(&_S980, &_S981, _S973.differential_0);
    float3  _S982 = make_float3 (_S919.x, _S919.y, _S979.differential_0 + _S981.differential_0 + _S918.x + _S918.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S983;
    (&_S983)->primal_0 = R_4;
    (&_S983)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S984;
    (&_S984)->primal_0 = pos_i_0;
    (&_S984)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S983, &_S984, _S982);
    DiffPair_float_0 _S985;
    (&_S985)->primal_0 = 0.0f;
    (&_S985)->differential_0 = 0.0f;
    DiffPair_float_0 _S986;
    (&_S986)->primal_0 = _S498;
    (&_S986)->differential_0 = 0.0f;
    _d_max_0(&_S985, &_S986, _S978.differential_0);
    DiffPair_float_0 _S987;
    (&_S987)->primal_0 = 1.00000001504746622e+30f;
    (&_S987)->differential_0 = 0.0f;
    DiffPair_float_0 _S988;
    (&_S988)->primal_0 = _S498;
    (&_S988)->differential_0 = 0.0f;
    _d_min_0(&_S987, &_S988, _S980.differential_0);
    float3  _S989 = make_float3 (_S933.x, _S933.y, _S986.differential_0 + _S988.differential_0 + _S932.x + _S932.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S990;
    (&_S990)->primal_0 = R_4;
    (&_S990)->differential_0 = _S762;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S991;
    (&_S991)->primal_0 = pos_4;
    (&_S991)->differential_0 = _S494;
    s_bwd_prop_mul_0(&_S990, &_S991, _S989);
    FixedArray<float3 , 16>  _S992;
    _S992[int(0)] = _S494;
    _S992[int(1)] = _S494;
    _S992[int(2)] = _S494;
    _S992[int(3)] = _S494;
    _S992[int(4)] = _S494;
    _S992[int(5)] = _S494;
    _S992[int(6)] = _S494;
    _S992[int(7)] = _S494;
    _S992[int(8)] = _S494;
    _S992[int(9)] = _S494;
    _S992[int(10)] = _S494;
    _S992[int(11)] = _S494;
    _S992[int(12)] = _S494;
    _S992[int(13)] = _S494;
    _S992[int(14)] = _S494;
    _S992[int(15)] = _S494;
    _S992[int(7)] = _S725;
    _S992[int(0)] = _S759;
    _S992[int(1)] = _S746;
    _S992[int(2)] = _S744;
    _S992[int(3)] = _S742;
    _S992[int(4)] = _S731;
    _S992[int(5)] = _S729;
    _S992[int(6)] = _S727;
    _S992[int(8)] = _S723;
    _S992[int(9)] = _S715;
    _S992[int(10)] = _S713;
    _S992[int(11)] = _S711;
    _S992[int(12)] = _S709;
    _S992[int(13)] = _S707;
    _S992[int(14)] = _S705;
    _S992[int(15)] = _S703;
    float3  _S993 = _S764.differential_0 + _S760 + _S940 + _S947 + _S954 + _S961 + _S968 + _S975 + _S982 + _S989;
    Matrix<float, 3, 3>  _S994 = _S765 + _S934.differential_0 + _S941.differential_0 + _S948.differential_0 + _S955.differential_0 + _S962.differential_0 + _S969.differential_0 + _S976.differential_0 + _S983.differential_0 + _S990.differential_0;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    (*v_sh_coeffs_0)[int(0)] = _S992[int(0)];
    (*v_sh_coeffs_0)[int(1)] = _S992[int(1)];
    (*v_sh_coeffs_0)[int(2)] = _S992[int(2)];
    (*v_sh_coeffs_0)[int(3)] = _S992[int(3)];
    (*v_sh_coeffs_0)[int(4)] = _S992[int(4)];
    (*v_sh_coeffs_0)[int(5)] = _S992[int(5)];
    (*v_sh_coeffs_0)[int(6)] = _S992[int(6)];
    (*v_sh_coeffs_0)[int(7)] = _S992[int(7)];
    (*v_sh_coeffs_0)[int(8)] = _S992[int(8)];
    (*v_sh_coeffs_0)[int(9)] = _S992[int(9)];
    (*v_sh_coeffs_0)[int(10)] = _S992[int(10)];
    (*v_sh_coeffs_0)[int(11)] = _S992[int(11)];
    (*v_sh_coeffs_0)[int(12)] = _S992[int(12)];
    (*v_sh_coeffs_0)[int(13)] = _S992[int(13)];
    (*v_sh_coeffs_0)[int(14)] = _S992[int(14)];
    (*v_sh_coeffs_0)[int(15)] = _S992[int(15)];
    *v_R_0 = _S994;
    *v_t_0 = _S993;
    return;
}

struct s_bwd_prop_projection_voxel_eval3d_fisheye_differentiable_Intermediates_0
{
    FixedArray<float3 , 16>  _S995;
};

inline __device__ float s_primal_ctx_atan2_0(float _S996, float _S997)
{
    return (F32_atan2((_S996), (_S997)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S998, DiffPair_float_0 * _S999, float _S1000)
{
    _d_atan2_0(_S998, _S999, _S1000);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_6, float _s_dOut_0)
{
    float _S1001 = (*dpx_6).primal_0.x;
    float _S1002 = (*dpx_6).primal_0.y;
    DiffPair_float_0 _S1003;
    (&_S1003)->primal_0 = _S1001 * _S1001 + _S1002 * _S1002;
    (&_S1003)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1003, _s_dOut_0);
    float _S1004 = (*dpx_6).primal_0.y * _S1003.differential_0;
    float _S1005 = _S1004 + _S1004;
    float _S1006 = (*dpx_6).primal_0.x * _S1003.differential_0;
    float _S1007 = _S1006 + _S1006;
    float2  _S1008 = make_float2 (0.0f);
    *&((&_S1008)->y) = _S1005;
    *&((&_S1008)->x) = _S1007;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S1008;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1009, float _S1010)
{
    s_bwd_prop_length_impl_0(_S1009, _S1010);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_1)
{
    float _S1011 = (*dpx_7).primal_0.x;
    float _S1012 = (*dpx_7).primal_0.y;
    float _S1013 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S1014;
    (&_S1014)->primal_0 = _S1011 * _S1011 + _S1012 * _S1012 + _S1013 * _S1013;
    (&_S1014)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1014, _s_dOut_1);
    float _S1015 = (*dpx_7).primal_0.z * _S1014.differential_0;
    float _S1016 = _S1015 + _S1015;
    float _S1017 = (*dpx_7).primal_0.y * _S1014.differential_0;
    float _S1018 = _S1017 + _S1017;
    float _S1019 = (*dpx_7).primal_0.x * _S1014.differential_0;
    float _S1020 = _S1019 + _S1019;
    float3  _S1021 = make_float3 (0.0f);
    *&((&_S1021)->z) = _S1016;
    *&((&_S1021)->y) = _S1018;
    *&((&_S1021)->x) = _S1020;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S1021;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1022, float _S1023)
{
    s_bwd_prop_length_impl_1(_S1022, _S1023);
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  * densities_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_5, float fy_5, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_5, uint image_width_5, uint image_height_5, float3  v_rgb_1, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  _S1024 = make_float3 (0.0f);
    FixedArray<float3 , 16>  _S1025 = {
        _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024, _S1024
    };
    s_bwd_prop_projection_voxel_eval3d_fisheye_differentiable_Intermediates_0 _S1026;
    (&_S1026)->_S995 = _S1025;
    (&_S1026)->_S995 = *sh_coeffs_5;
    s_bwd_prop_projection_voxel_eval3d_fisheye_differentiable_Intermediates_0 _S1027 = _S1026;
    float3  _S1028 = s_primal_ctx_mul_0(R_5, pos_5) + t_5;
    float _S1029 = length_1(_S1028);
    float _S1030 = (F32_min((1.00000001504746622e+30f), (_S1029)));
    float _S1031 = (F32_max((0.0f), (_S1029)));
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S1032 = s_primal_ctx_mul_0(R_5, pos_i_7) + t_5;
    float _S1033 = length_1(_S1032);
    float _S1034 = (F32_min((_S1030), (_S1033)));
    float _S1035 = (F32_max((_S1031), (_S1033)));
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S1036 = s_primal_ctx_mul_0(R_5, pos_i_8) + t_5;
    float _S1037 = length_1(_S1036);
    float _S1038 = (F32_min((_S1034), (_S1037)));
    float _S1039 = (F32_max((_S1035), (_S1037)));
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S1040 = s_primal_ctx_mul_0(R_5, pos_i_9) + t_5;
    float _S1041 = length_1(_S1040);
    float _S1042 = (F32_min((_S1038), (_S1041)));
    float _S1043 = (F32_max((_S1039), (_S1041)));
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S1044 = s_primal_ctx_mul_0(R_5, pos_i_10) + t_5;
    float _S1045 = length_1(_S1044);
    float _S1046 = (F32_min((_S1042), (_S1045)));
    float _S1047 = (F32_max((_S1043), (_S1045)));
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S1048 = s_primal_ctx_mul_0(R_5, pos_i_11) + t_5;
    float _S1049 = length_1(_S1048);
    float _S1050 = (F32_min((_S1046), (_S1049)));
    float _S1051 = (F32_max((_S1047), (_S1049)));
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S1052 = s_primal_ctx_mul_0(R_5, pos_i_12) + t_5;
    float _S1053 = length_1(_S1052);
    float _S1054 = (F32_min((_S1050), (_S1053)));
    float _S1055 = (F32_max((_S1051), (_S1053)));
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S1056 = s_primal_ctx_mul_0(R_5, pos_i_13) + t_5;
    float _S1057 = length_1(_S1056);
    float3  _S1058 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_4 = s_primal_ctx_mul_0(R_5, _S1058) + t_5;
    float2  _S1059 = float2 {_S1028.x, _S1028.y};
    float _S1060 = length_0(_S1059);
    float _S1061 = _S1028.z;
    float _S1062 = s_primal_ctx_atan2_0(_S1060, _S1061);
    bool _S1063 = _S1062 < 0.00100000004749745f;
    float k_2;
    float _S1064;
    float _S1065;
    float _S1066;
    if(_S1063)
    {
        float _S1067 = 1.0f - _S1062 * _S1062 / 3.0f;
        float _S1068 = _S1061 * _S1061;
        k_2 = _S1067 / _S1061;
        _S1064 = _S1068;
        _S1065 = _S1067;
        _S1066 = 0.0f;
    }
    else
    {
        float _S1069 = _S1060 * _S1060;
        k_2 = _S1062 / _S1060;
        _S1064 = 0.0f;
        _S1065 = 0.0f;
        _S1066 = _S1069;
    }
    float2  _S1070 = make_float2 (k_2);
    float2  _S1071 = _S1059 * make_float2 (k_2);
    float u_56 = _S1071.x;
    float v_56 = _S1071.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S1072 = dist_coeffs_5[int(2)] + r2_56 * dist_coeffs_5[int(3)];
    float _S1073 = dist_coeffs_5[int(1)] + r2_56 * _S1072;
    float _S1074 = dist_coeffs_5[int(0)] + r2_56 * _S1073;
    float radial_24 = 1.0f + r2_56 * _S1074;
    float _S1075 = 2.0f * dist_coeffs_5[int(4)];
    float _S1076 = _S1075 * u_56;
    float _S1077 = 2.0f * u_56;
    float _S1078 = 2.0f * dist_coeffs_5[int(5)];
    float _S1079 = _S1078 * u_56;
    float _S1080 = 2.0f * v_56;
    float2  _S1081 = _S1071 * make_float2 (radial_24) + make_float2 (_S1076 * v_56 + dist_coeffs_5[int(5)] * (r2_56 + _S1077 * u_56) + dist_coeffs_5[int(6)] * r2_56, _S1079 * v_56 + dist_coeffs_5[int(4)] * (r2_56 + _S1080 * v_56) + dist_coeffs_5[int(7)] * r2_56);
    float2  _S1082 = _S1081 + make_float2 (dist_coeffs_5[int(8)] * _S1081.x + dist_coeffs_5[int(9)] * _S1081.y, 0.0f);
    float _S1083 = fx_5 * _S1082.x + cx_5;
    float _S1084 = fy_5 * _S1082.y + cy_5;
    float2  _S1085 = float2 {_S1032.x, _S1032.y};
    float _S1086 = length_0(_S1085);
    float _S1087 = _S1032.z;
    float _S1088 = s_primal_ctx_atan2_0(_S1086, _S1087);
    bool _S1089 = _S1088 < 0.00100000004749745f;
    float _S1090;
    float _S1091;
    float _S1092;
    if(_S1089)
    {
        float _S1093 = 1.0f - _S1088 * _S1088 / 3.0f;
        float _S1094 = _S1087 * _S1087;
        k_2 = _S1093 / _S1087;
        _S1090 = _S1094;
        _S1091 = _S1093;
        _S1092 = 0.0f;
    }
    else
    {
        float _S1095 = _S1086 * _S1086;
        k_2 = _S1088 / _S1086;
        _S1090 = 0.0f;
        _S1091 = 0.0f;
        _S1092 = _S1095;
    }
    float2  _S1096 = make_float2 (k_2);
    float2  _S1097 = _S1085 * make_float2 (k_2);
    float u_57 = _S1097.x;
    float v_57 = _S1097.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S1098 = dist_coeffs_5[int(2)] + r2_57 * dist_coeffs_5[int(3)];
    float _S1099 = dist_coeffs_5[int(1)] + r2_57 * _S1098;
    float _S1100 = dist_coeffs_5[int(0)] + r2_57 * _S1099;
    float radial_25 = 1.0f + r2_57 * _S1100;
    float _S1101 = _S1075 * u_57;
    float _S1102 = 2.0f * u_57;
    float _S1103 = _S1078 * u_57;
    float _S1104 = 2.0f * v_57;
    float2  _S1105 = _S1097 * make_float2 (radial_25) + make_float2 (_S1101 * v_57 + dist_coeffs_5[int(5)] * (r2_57 + _S1102 * u_57) + dist_coeffs_5[int(6)] * r2_57, _S1103 * v_57 + dist_coeffs_5[int(4)] * (r2_57 + _S1104 * v_57) + dist_coeffs_5[int(7)] * r2_57);
    float2  _S1106 = _S1105 + make_float2 (dist_coeffs_5[int(8)] * _S1105.x + dist_coeffs_5[int(9)] * _S1105.y, 0.0f);
    float _S1107 = fx_5 * _S1106.x + cx_5;
    float _S1108 = fy_5 * _S1106.y + cy_5;
    float2  _S1109 = float2 {_S1036.x, _S1036.y};
    float _S1110 = length_0(_S1109);
    float _S1111 = _S1036.z;
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
    float u_58 = _S1121.x;
    float v_58 = _S1121.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S1122 = dist_coeffs_5[int(2)] + r2_58 * dist_coeffs_5[int(3)];
    float _S1123 = dist_coeffs_5[int(1)] + r2_58 * _S1122;
    float _S1124 = dist_coeffs_5[int(0)] + r2_58 * _S1123;
    float radial_26 = 1.0f + r2_58 * _S1124;
    float _S1125 = _S1075 * u_58;
    float _S1126 = 2.0f * u_58;
    float _S1127 = _S1078 * u_58;
    float _S1128 = 2.0f * v_58;
    float2  _S1129 = _S1121 * make_float2 (radial_26) + make_float2 (_S1125 * v_58 + dist_coeffs_5[int(5)] * (r2_58 + _S1126 * u_58) + dist_coeffs_5[int(6)] * r2_58, _S1127 * v_58 + dist_coeffs_5[int(4)] * (r2_58 + _S1128 * v_58) + dist_coeffs_5[int(7)] * r2_58);
    float2  _S1130 = _S1129 + make_float2 (dist_coeffs_5[int(8)] * _S1129.x + dist_coeffs_5[int(9)] * _S1129.y, 0.0f);
    float _S1131 = fx_5 * _S1130.x + cx_5;
    float _S1132 = fy_5 * _S1130.y + cy_5;
    float2  _S1133 = float2 {_S1040.x, _S1040.y};
    float _S1134 = length_0(_S1133);
    float _S1135 = _S1040.z;
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
    float u_59 = _S1145.x;
    float v_59 = _S1145.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S1146 = dist_coeffs_5[int(2)] + r2_59 * dist_coeffs_5[int(3)];
    float _S1147 = dist_coeffs_5[int(1)] + r2_59 * _S1146;
    float _S1148 = dist_coeffs_5[int(0)] + r2_59 * _S1147;
    float radial_27 = 1.0f + r2_59 * _S1148;
    float _S1149 = _S1075 * u_59;
    float _S1150 = 2.0f * u_59;
    float _S1151 = _S1078 * u_59;
    float _S1152 = 2.0f * v_59;
    float2  _S1153 = _S1145 * make_float2 (radial_27) + make_float2 (_S1149 * v_59 + dist_coeffs_5[int(5)] * (r2_59 + _S1150 * u_59) + dist_coeffs_5[int(6)] * r2_59, _S1151 * v_59 + dist_coeffs_5[int(4)] * (r2_59 + _S1152 * v_59) + dist_coeffs_5[int(7)] * r2_59);
    float2  _S1154 = _S1153 + make_float2 (dist_coeffs_5[int(8)] * _S1153.x + dist_coeffs_5[int(9)] * _S1153.y, 0.0f);
    float _S1155 = fx_5 * _S1154.x + cx_5;
    float _S1156 = fy_5 * _S1154.y + cy_5;
    float2  _S1157 = float2 {_S1044.x, _S1044.y};
    float _S1158 = length_0(_S1157);
    float _S1159 = _S1044.z;
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
    float u_60 = _S1169.x;
    float v_60 = _S1169.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S1170 = dist_coeffs_5[int(2)] + r2_60 * dist_coeffs_5[int(3)];
    float _S1171 = dist_coeffs_5[int(1)] + r2_60 * _S1170;
    float _S1172 = dist_coeffs_5[int(0)] + r2_60 * _S1171;
    float radial_28 = 1.0f + r2_60 * _S1172;
    float _S1173 = _S1075 * u_60;
    float _S1174 = 2.0f * u_60;
    float _S1175 = _S1078 * u_60;
    float _S1176 = 2.0f * v_60;
    float2  _S1177 = _S1169 * make_float2 (radial_28) + make_float2 (_S1173 * v_60 + dist_coeffs_5[int(5)] * (r2_60 + _S1174 * u_60) + dist_coeffs_5[int(6)] * r2_60, _S1175 * v_60 + dist_coeffs_5[int(4)] * (r2_60 + _S1176 * v_60) + dist_coeffs_5[int(7)] * r2_60);
    float2  _S1178 = _S1177 + make_float2 (dist_coeffs_5[int(8)] * _S1177.x + dist_coeffs_5[int(9)] * _S1177.y, 0.0f);
    float _S1179 = fx_5 * _S1178.x + cx_5;
    float _S1180 = fy_5 * _S1178.y + cy_5;
    float2  _S1181 = float2 {_S1048.x, _S1048.y};
    float _S1182 = length_0(_S1181);
    float _S1183 = _S1048.z;
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
    float u_61 = _S1193.x;
    float v_61 = _S1193.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S1194 = dist_coeffs_5[int(2)] + r2_61 * dist_coeffs_5[int(3)];
    float _S1195 = dist_coeffs_5[int(1)] + r2_61 * _S1194;
    float _S1196 = dist_coeffs_5[int(0)] + r2_61 * _S1195;
    float radial_29 = 1.0f + r2_61 * _S1196;
    float _S1197 = _S1075 * u_61;
    float _S1198 = 2.0f * u_61;
    float _S1199 = _S1078 * u_61;
    float _S1200 = 2.0f * v_61;
    float2  _S1201 = _S1193 * make_float2 (radial_29) + make_float2 (_S1197 * v_61 + dist_coeffs_5[int(5)] * (r2_61 + _S1198 * u_61) + dist_coeffs_5[int(6)] * r2_61, _S1199 * v_61 + dist_coeffs_5[int(4)] * (r2_61 + _S1200 * v_61) + dist_coeffs_5[int(7)] * r2_61);
    float2  _S1202 = _S1201 + make_float2 (dist_coeffs_5[int(8)] * _S1201.x + dist_coeffs_5[int(9)] * _S1201.y, 0.0f);
    float _S1203 = fx_5 * _S1202.x + cx_5;
    float _S1204 = fy_5 * _S1202.y + cy_5;
    float2  _S1205 = float2 {_S1052.x, _S1052.y};
    float _S1206 = length_0(_S1205);
    float _S1207 = _S1052.z;
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
    float u_62 = _S1217.x;
    float v_62 = _S1217.y;
    float r2_62 = u_62 * u_62 + v_62 * v_62;
    float _S1218 = dist_coeffs_5[int(2)] + r2_62 * dist_coeffs_5[int(3)];
    float _S1219 = dist_coeffs_5[int(1)] + r2_62 * _S1218;
    float _S1220 = dist_coeffs_5[int(0)] + r2_62 * _S1219;
    float radial_30 = 1.0f + r2_62 * _S1220;
    float _S1221 = _S1075 * u_62;
    float _S1222 = 2.0f * u_62;
    float _S1223 = _S1078 * u_62;
    float _S1224 = 2.0f * v_62;
    float2  _S1225 = _S1217 * make_float2 (radial_30) + make_float2 (_S1221 * v_62 + dist_coeffs_5[int(5)] * (r2_62 + _S1222 * u_62) + dist_coeffs_5[int(6)] * r2_62, _S1223 * v_62 + dist_coeffs_5[int(4)] * (r2_62 + _S1224 * v_62) + dist_coeffs_5[int(7)] * r2_62);
    float2  _S1226 = _S1225 + make_float2 (dist_coeffs_5[int(8)] * _S1225.x + dist_coeffs_5[int(9)] * _S1225.y, 0.0f);
    float _S1227 = fx_5 * _S1226.x + cx_5;
    float _S1228 = fy_5 * _S1226.y + cy_5;
    float2  _S1229 = float2 {_S1056.x, _S1056.y};
    float _S1230 = length_0(_S1229);
    float _S1231 = _S1056.z;
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
    float u_63 = _S1241.x;
    float v_63 = _S1241.y;
    float r2_63 = u_63 * u_63 + v_63 * v_63;
    float _S1242 = dist_coeffs_5[int(2)] + r2_63 * dist_coeffs_5[int(3)];
    float _S1243 = dist_coeffs_5[int(1)] + r2_63 * _S1242;
    float _S1244 = dist_coeffs_5[int(0)] + r2_63 * _S1243;
    float radial_31 = 1.0f + r2_63 * _S1244;
    float _S1245 = _S1075 * u_63;
    float _S1246 = 2.0f * u_63;
    float _S1247 = _S1078 * u_63;
    float _S1248 = 2.0f * v_63;
    float2  _S1249 = _S1241 * make_float2 (radial_31) + make_float2 (_S1245 * v_63 + dist_coeffs_5[int(5)] * (r2_63 + _S1246 * u_63) + dist_coeffs_5[int(6)] * r2_63, _S1247 * v_63 + dist_coeffs_5[int(4)] * (r2_63 + _S1248 * v_63) + dist_coeffs_5[int(7)] * r2_63);
    float2  _S1250 = _S1249 + make_float2 (dist_coeffs_5[int(8)] * _S1249.x + dist_coeffs_5[int(9)] * _S1249.y, 0.0f);
    float _S1251 = fx_5 * _S1250.x + cx_5;
    float _S1252 = fy_5 * _S1250.y + cy_5;
    float _S1253 = (F32_max((_S1083), (_S1107)));
    float _S1254 = (F32_min((_S1083), (_S1107)));
    float _S1255 = (F32_max((_S1084), (_S1108)));
    float _S1256 = (F32_min((_S1084), (_S1108)));
    float _S1257 = (F32_max((_S1253), (_S1131)));
    float _S1258 = (F32_min((_S1254), (_S1131)));
    float _S1259 = (F32_max((_S1255), (_S1132)));
    float _S1260 = (F32_min((_S1256), (_S1132)));
    float _S1261 = (F32_max((_S1257), (_S1155)));
    float _S1262 = (F32_min((_S1258), (_S1155)));
    float _S1263 = (F32_max((_S1259), (_S1156)));
    float _S1264 = (F32_min((_S1260), (_S1156)));
    float _S1265 = (F32_max((_S1261), (_S1179)));
    float _S1266 = (F32_min((_S1262), (_S1179)));
    float _S1267 = (F32_max((_S1263), (_S1180)));
    float _S1268 = (F32_min((_S1264), (_S1180)));
    float _S1269 = (F32_max((_S1265), (_S1203)));
    float _S1270 = (F32_min((_S1266), (_S1203)));
    float _S1271 = (F32_max((_S1267), (_S1204)));
    float _S1272 = (F32_min((_S1268), (_S1204)));
    float _S1273 = (F32_max((_S1269), (_S1227)));
    float _S1274 = (F32_min((_S1270), (_S1227)));
    float _S1275 = (F32_max((_S1271), (_S1228)));
    float _S1276 = (F32_min((_S1272), (_S1228)));
    Matrix<float, 3, 3>  _S1277 = transpose_1(R_5);
    float3  _S1278 = mean_c_4 - - s_primal_ctx_mul_0(_S1277, t_5);
    float _S1279 = _S1278.x;
    float _S1280 = _S1278.y;
    float _S1281 = _S1278.z;
    float _S1282 = _S1279 * _S1279 + _S1280 * _S1280 + _S1281 * _S1281;
    float _S1283 = s_primal_ctx_sqrt_0(_S1282);
    float x_14 = _S1279 / _S1283;
    float3  _S1284 = make_float3 (x_14);
    float _S1285 = _S1283 * _S1283;
    float y_10 = _S1280 / _S1283;
    float z_7 = _S1281 / _S1283;
    float3  _S1286 = make_float3 (z_7);
    float _S1287 = - y_10;
    float3  _S1288 = make_float3 (_S1287);
    float z2_5 = z_7 * z_7;
    float fTmp0B_5 = -1.09254848957061768f * z_7;
    float fC1_5 = x_14 * x_14 - y_10 * y_10;
    float _S1289 = 2.0f * x_14;
    float fS1_5 = _S1289 * y_10;
    float pSH6_1 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float3  _S1290 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_5 * x_14;
    float3  _S1291 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_5 * y_10;
    float3  _S1292 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_5;
    float3  _S1293 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_5;
    float3  _S1294 = make_float3 (pSH4_1);
    float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_7;
    float _S1295 = 1.86588168144226074f * z2_5 - 1.11952900886535645f;
    float pSH12_1 = z_7 * _S1295;
    float3  _S1296 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_5 * x_14;
    float3  _S1297 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_5 * y_10;
    float3  _S1298 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_5 * fC1_5;
    float3  _S1299 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_5 * fS1_5;
    float3  _S1300 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_14 * fC1_5 - y_10 * fS1_5);
    float3  _S1301 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_14 * fS1_5 + y_10 * fC1_5);
    float3  _S1302 = make_float3 (pSH9_1);
    float3  _S1303 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1304;
    (&_S1304)->primal_0 = make_float3 (0.282094806432724f) * _S1027._S995[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1287) * _S1027._S995[int(1)] + make_float3 (z_7) * _S1027._S995[int(2)] - make_float3 (x_14) * _S1027._S995[int(3)]) + (make_float3 (pSH4_1) * _S1027._S995[int(4)] + make_float3 (pSH5_1) * _S1027._S995[int(5)] + make_float3 (pSH6_1) * _S1027._S995[int(6)] + make_float3 (pSH7_1) * _S1027._S995[int(7)] + make_float3 (pSH8_1) * _S1027._S995[int(8)]) + (make_float3 (pSH9_1) * _S1027._S995[int(9)] + make_float3 (pSH10_1) * _S1027._S995[int(10)] + make_float3 (pSH11_1) * _S1027._S995[int(11)] + make_float3 (pSH12_1) * _S1027._S995[int(12)] + make_float3 (pSH13_1) * _S1027._S995[int(13)] + make_float3 (pSH14_1) * _S1027._S995[int(14)] + make_float3 (pSH15_1) * _S1027._S995[int(15)]) + make_float3 (0.5f);
    (&_S1304)->differential_0 = _S1024;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1305;
    (&_S1305)->primal_0 = _S1303;
    (&_S1305)->differential_0 = _S1024;
    s_bwd_prop_max_0(&_S1304, &_S1305, v_rgb_1);
    float3  _S1306 = _S1301 * _S1304.differential_0;
    float3  _S1307 = _S1027._S995[int(15)] * _S1304.differential_0;
    float3  _S1308 = _S1299 * _S1304.differential_0;
    float3  _S1309 = _S1027._S995[int(14)] * _S1304.differential_0;
    float3  _S1310 = _S1297 * _S1304.differential_0;
    float3  _S1311 = _S1027._S995[int(13)] * _S1304.differential_0;
    float3  _S1312 = _S1296 * _S1304.differential_0;
    float3  _S1313 = _S1027._S995[int(12)] * _S1304.differential_0;
    float3  _S1314 = _S1298 * _S1304.differential_0;
    float3  _S1315 = _S1027._S995[int(11)] * _S1304.differential_0;
    float3  _S1316 = _S1300 * _S1304.differential_0;
    float3  _S1317 = _S1027._S995[int(10)] * _S1304.differential_0;
    float3  _S1318 = _S1302 * _S1304.differential_0;
    float3  _S1319 = _S1027._S995[int(9)] * _S1304.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1319.x + _S1319.y + _S1319.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1307.x + _S1307.y + _S1307.z);
    float _S1320 = _S1317.x + _S1317.y + _S1317.z;
    float _S1321 = _S1309.x + _S1309.y + _S1309.z;
    float _S1322 = _S1315.x + _S1315.y + _S1315.z;
    float _S1323 = _S1311.x + _S1311.y + _S1311.z;
    float _S1324 = _S1313.x + _S1313.y + _S1313.z;
    float _S1325 = - s_diff_fC2_T_1;
    float3  _S1326 = _S1293 * _S1304.differential_0;
    float3  _S1327 = _S1027._S995[int(8)] * _S1304.differential_0;
    float3  _S1328 = _S1291 * _S1304.differential_0;
    float3  _S1329 = _S1027._S995[int(7)] * _S1304.differential_0;
    float3  _S1330 = _S1290 * _S1304.differential_0;
    float3  _S1331 = _S1027._S995[int(6)] * _S1304.differential_0;
    float3  _S1332 = _S1292 * _S1304.differential_0;
    float3  _S1333 = _S1027._S995[int(5)] * _S1304.differential_0;
    float3  _S1334 = _S1294 * _S1304.differential_0;
    float3  _S1335 = _S1027._S995[int(4)] * _S1304.differential_0;
    float _S1336 = _S1333.x + _S1333.y + _S1333.z;
    float _S1337 = _S1329.x + _S1329.y + _S1329.z;
    float _S1338 = fTmp1B_5 * _S1320 + x_14 * s_diff_fS2_T_1 + y_10 * _S1325 + 0.54627424478530884f * (_S1335.x + _S1335.y + _S1335.z);
    float _S1339 = fTmp1B_5 * _S1321 + y_10 * s_diff_fS2_T_1 + x_14 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1327.x + _S1327.y + _S1327.z);
    float _S1340 = y_10 * - _S1339;
    float _S1341 = x_14 * _S1339;
    float _S1342 = z_7 * (1.86588168144226074f * (z_7 * _S1324) + -2.28522896766662598f * (y_10 * _S1322 + x_14 * _S1323) + 0.94617468118667603f * (_S1331.x + _S1331.y + _S1331.z));
    float3  _S1343 = make_float3 (0.48860251903533936f) * _S1304.differential_0;
    float3  _S1344 = - _S1343;
    float3  _S1345 = _S1284 * _S1344;
    float3  _S1346 = _S1027._S995[int(3)] * _S1344;
    float3  _S1347 = _S1286 * _S1343;
    float3  _S1348 = _S1027._S995[int(2)] * _S1343;
    float3  _S1349 = _S1288 * _S1343;
    float3  _S1350 = _S1027._S995[int(1)] * _S1343;
    float _S1351 = (_S1295 * _S1324 + 1.44530570507049561f * (fS1_5 * _S1320 + fC1_5 * _S1321) + -1.09254848957061768f * (y_10 * _S1336 + x_14 * _S1337) + _S1342 + _S1342 + _S1348.x + _S1348.y + _S1348.z) / _S1285;
    float _S1352 = _S1283 * _S1351;
    float _S1353 = (fTmp0C_5 * _S1322 + fC1_5 * s_diff_fS2_T_1 + fS1_5 * _S1325 + fTmp0B_5 * _S1336 + _S1289 * _S1338 + _S1340 + _S1340 + - (_S1350.x + _S1350.y + _S1350.z)) / _S1285;
    float _S1354 = _S1283 * _S1353;
    float _S1355 = (fTmp0C_5 * _S1323 + fS1_5 * s_diff_fS2_T_1 + fC1_5 * s_diff_fC2_T_1 + fTmp0B_5 * _S1337 + 2.0f * (y_10 * _S1338) + _S1341 + _S1341 + _S1346.x + _S1346.y + _S1346.z) / _S1285;
    float _S1356 = _S1283 * _S1355;
    float _S1357 = _S1281 * - _S1351 + _S1280 * - _S1353 + _S1279 * - _S1355;
    DiffPair_float_0 _S1358;
    (&_S1358)->primal_0 = _S1282;
    (&_S1358)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1358, _S1357);
    float _S1359 = _S1281 * _S1358.differential_0;
    float _S1360 = _S1280 * _S1358.differential_0;
    float _S1361 = _S1279 * _S1358.differential_0;
    float3  _S1362 = make_float3 (0.282094806432724f) * _S1304.differential_0;
    float3  _S1363 = make_float3 (_S1356 + _S1361 + _S1361, _S1354 + _S1360 + _S1360, _S1352 + _S1359 + _S1359);
    float3  _S1364 = - - _S1363;
    Matrix<float, 3, 3>  _S1365 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1366;
    (&_S1366)->primal_0 = _S1277;
    (&_S1366)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1367;
    (&_S1367)->primal_0 = t_5;
    (&_S1367)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1366, &_S1367, _S1364);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1368 = _S1367;
    Matrix<float, 3, 3>  _S1369 = transpose_1(_S1366.differential_0);
    DiffPair_float_0 _S1370;
    (&_S1370)->primal_0 = _S1276;
    (&_S1370)->differential_0 = 0.0f;
    DiffPair_float_0 _S1371;
    (&_S1371)->primal_0 = _S1252;
    (&_S1371)->differential_0 = 0.0f;
    _d_min_0(&_S1370, &_S1371, 0.0f);
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1275;
    (&_S1372)->differential_0 = 0.0f;
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S1252;
    (&_S1373)->differential_0 = 0.0f;
    _d_max_0(&_S1372, &_S1373, 0.0f);
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = _S1274;
    (&_S1374)->differential_0 = 0.0f;
    DiffPair_float_0 _S1375;
    (&_S1375)->primal_0 = _S1251;
    (&_S1375)->differential_0 = 0.0f;
    _d_min_0(&_S1374, &_S1375, 0.0f);
    DiffPair_float_0 _S1376;
    (&_S1376)->primal_0 = _S1273;
    (&_S1376)->differential_0 = 0.0f;
    DiffPair_float_0 _S1377;
    (&_S1377)->primal_0 = _S1251;
    (&_S1377)->differential_0 = 0.0f;
    _d_max_0(&_S1376, &_S1377, 0.0f);
    DiffPair_float_0 _S1378;
    (&_S1378)->primal_0 = _S1272;
    (&_S1378)->differential_0 = 0.0f;
    DiffPair_float_0 _S1379;
    (&_S1379)->primal_0 = _S1228;
    (&_S1379)->differential_0 = 0.0f;
    _d_min_0(&_S1378, &_S1379, _S1370.differential_0);
    DiffPair_float_0 _S1380;
    (&_S1380)->primal_0 = _S1271;
    (&_S1380)->differential_0 = 0.0f;
    DiffPair_float_0 _S1381;
    (&_S1381)->primal_0 = _S1228;
    (&_S1381)->differential_0 = 0.0f;
    _d_max_0(&_S1380, &_S1381, _S1372.differential_0);
    DiffPair_float_0 _S1382;
    (&_S1382)->primal_0 = _S1270;
    (&_S1382)->differential_0 = 0.0f;
    DiffPair_float_0 _S1383;
    (&_S1383)->primal_0 = _S1227;
    (&_S1383)->differential_0 = 0.0f;
    _d_min_0(&_S1382, &_S1383, _S1374.differential_0);
    DiffPair_float_0 _S1384;
    (&_S1384)->primal_0 = _S1269;
    (&_S1384)->differential_0 = 0.0f;
    DiffPair_float_0 _S1385;
    (&_S1385)->primal_0 = _S1227;
    (&_S1385)->differential_0 = 0.0f;
    _d_max_0(&_S1384, &_S1385, _S1376.differential_0);
    DiffPair_float_0 _S1386;
    (&_S1386)->primal_0 = _S1268;
    (&_S1386)->differential_0 = 0.0f;
    DiffPair_float_0 _S1387;
    (&_S1387)->primal_0 = _S1204;
    (&_S1387)->differential_0 = 0.0f;
    _d_min_0(&_S1386, &_S1387, _S1378.differential_0);
    DiffPair_float_0 _S1388;
    (&_S1388)->primal_0 = _S1267;
    (&_S1388)->differential_0 = 0.0f;
    DiffPair_float_0 _S1389;
    (&_S1389)->primal_0 = _S1204;
    (&_S1389)->differential_0 = 0.0f;
    _d_max_0(&_S1388, &_S1389, _S1380.differential_0);
    DiffPair_float_0 _S1390;
    (&_S1390)->primal_0 = _S1266;
    (&_S1390)->differential_0 = 0.0f;
    DiffPair_float_0 _S1391;
    (&_S1391)->primal_0 = _S1203;
    (&_S1391)->differential_0 = 0.0f;
    _d_min_0(&_S1390, &_S1391, _S1382.differential_0);
    DiffPair_float_0 _S1392;
    (&_S1392)->primal_0 = _S1265;
    (&_S1392)->differential_0 = 0.0f;
    DiffPair_float_0 _S1393;
    (&_S1393)->primal_0 = _S1203;
    (&_S1393)->differential_0 = 0.0f;
    _d_max_0(&_S1392, &_S1393, _S1384.differential_0);
    DiffPair_float_0 _S1394;
    (&_S1394)->primal_0 = _S1264;
    (&_S1394)->differential_0 = 0.0f;
    DiffPair_float_0 _S1395;
    (&_S1395)->primal_0 = _S1180;
    (&_S1395)->differential_0 = 0.0f;
    _d_min_0(&_S1394, &_S1395, _S1386.differential_0);
    DiffPair_float_0 _S1396;
    (&_S1396)->primal_0 = _S1263;
    (&_S1396)->differential_0 = 0.0f;
    DiffPair_float_0 _S1397;
    (&_S1397)->primal_0 = _S1180;
    (&_S1397)->differential_0 = 0.0f;
    _d_max_0(&_S1396, &_S1397, _S1388.differential_0);
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1262;
    (&_S1398)->differential_0 = 0.0f;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1179;
    (&_S1399)->differential_0 = 0.0f;
    _d_min_0(&_S1398, &_S1399, _S1390.differential_0);
    DiffPair_float_0 _S1400;
    (&_S1400)->primal_0 = _S1261;
    (&_S1400)->differential_0 = 0.0f;
    DiffPair_float_0 _S1401;
    (&_S1401)->primal_0 = _S1179;
    (&_S1401)->differential_0 = 0.0f;
    _d_max_0(&_S1400, &_S1401, _S1392.differential_0);
    DiffPair_float_0 _S1402;
    (&_S1402)->primal_0 = _S1260;
    (&_S1402)->differential_0 = 0.0f;
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = _S1156;
    (&_S1403)->differential_0 = 0.0f;
    _d_min_0(&_S1402, &_S1403, _S1394.differential_0);
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = _S1259;
    (&_S1404)->differential_0 = 0.0f;
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = _S1156;
    (&_S1405)->differential_0 = 0.0f;
    _d_max_0(&_S1404, &_S1405, _S1396.differential_0);
    DiffPair_float_0 _S1406;
    (&_S1406)->primal_0 = _S1258;
    (&_S1406)->differential_0 = 0.0f;
    DiffPair_float_0 _S1407;
    (&_S1407)->primal_0 = _S1155;
    (&_S1407)->differential_0 = 0.0f;
    _d_min_0(&_S1406, &_S1407, _S1398.differential_0);
    DiffPair_float_0 _S1408;
    (&_S1408)->primal_0 = _S1257;
    (&_S1408)->differential_0 = 0.0f;
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1155;
    (&_S1409)->differential_0 = 0.0f;
    _d_max_0(&_S1408, &_S1409, _S1400.differential_0);
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = _S1256;
    (&_S1410)->differential_0 = 0.0f;
    DiffPair_float_0 _S1411;
    (&_S1411)->primal_0 = _S1132;
    (&_S1411)->differential_0 = 0.0f;
    _d_min_0(&_S1410, &_S1411, _S1402.differential_0);
    DiffPair_float_0 _S1412;
    (&_S1412)->primal_0 = _S1255;
    (&_S1412)->differential_0 = 0.0f;
    DiffPair_float_0 _S1413;
    (&_S1413)->primal_0 = _S1132;
    (&_S1413)->differential_0 = 0.0f;
    _d_max_0(&_S1412, &_S1413, _S1404.differential_0);
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = _S1254;
    (&_S1414)->differential_0 = 0.0f;
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1131;
    (&_S1415)->differential_0 = 0.0f;
    _d_min_0(&_S1414, &_S1415, _S1406.differential_0);
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1253;
    (&_S1416)->differential_0 = 0.0f;
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1131;
    (&_S1417)->differential_0 = 0.0f;
    _d_max_0(&_S1416, &_S1417, _S1408.differential_0);
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = _S1084;
    (&_S1418)->differential_0 = 0.0f;
    DiffPair_float_0 _S1419;
    (&_S1419)->primal_0 = _S1108;
    (&_S1419)->differential_0 = 0.0f;
    _d_min_0(&_S1418, &_S1419, _S1410.differential_0);
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = _S1084;
    (&_S1420)->differential_0 = 0.0f;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1108;
    (&_S1421)->differential_0 = 0.0f;
    _d_max_0(&_S1420, &_S1421, _S1412.differential_0);
    DiffPair_float_0 _S1422;
    (&_S1422)->primal_0 = _S1083;
    (&_S1422)->differential_0 = 0.0f;
    DiffPair_float_0 _S1423;
    (&_S1423)->primal_0 = _S1107;
    (&_S1423)->differential_0 = 0.0f;
    _d_min_0(&_S1422, &_S1423, _S1414.differential_0);
    DiffPair_float_0 _S1424;
    (&_S1424)->primal_0 = _S1083;
    (&_S1424)->differential_0 = 0.0f;
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = _S1107;
    (&_S1425)->differential_0 = 0.0f;
    _d_max_0(&_S1424, &_S1425, _S1416.differential_0);
    float _S1426 = fx_5 * (_S1375.differential_0 + _S1377.differential_0);
    float2  _S1427 = make_float2 (_S1426, fy_5 * (_S1371.differential_0 + _S1373.differential_0)) + make_float2 (dist_coeffs_5[int(8)] * _S1426, dist_coeffs_5[int(9)] * _S1426);
    float2  _S1428 = _S1241 * _S1427;
    float _S1429 = dist_coeffs_5[int(4)] * _S1427.y;
    float _S1430 = dist_coeffs_5[int(5)] * _S1427.x;
    float _S1431 = _S1428.x + _S1428.y;
    float _S1432 = r2_63 * _S1431;
    float _S1433 = r2_63 * _S1432;
    float _S1434 = dist_coeffs_5[int(7)] * _S1427.y + _S1429 + dist_coeffs_5[int(6)] * _S1427.x + _S1430 + _S1244 * _S1431 + _S1243 * _S1432 + _S1242 * _S1433 + dist_coeffs_5[int(3)] * (r2_63 * _S1433);
    float _S1435 = v_63 * _S1434;
    float _S1436 = u_63 * _S1434;
    float2  _S1437 = make_float2 (radial_31) * _S1427 + make_float2 (_S1078 * (v_63 * _S1427.y) + _S1246 * _S1430 + 2.0f * (u_63 * _S1430) + _S1075 * (v_63 * _S1427.x) + _S1436 + _S1436, _S1248 * _S1429 + 2.0f * (v_63 * _S1429) + _S1247 * _S1427.y + _S1245 * _S1427.x + _S1435 + _S1435);
    FixedArray<float3 , 16>  _S1438;
    _S1438[int(0)] = _S1024;
    _S1438[int(1)] = _S1024;
    _S1438[int(2)] = _S1024;
    _S1438[int(3)] = _S1024;
    _S1438[int(4)] = _S1024;
    _S1438[int(5)] = _S1024;
    _S1438[int(6)] = _S1024;
    _S1438[int(7)] = _S1024;
    _S1438[int(8)] = _S1024;
    _S1438[int(9)] = _S1024;
    _S1438[int(10)] = _S1024;
    _S1438[int(11)] = _S1024;
    _S1438[int(12)] = _S1024;
    _S1438[int(13)] = _S1024;
    _S1438[int(14)] = _S1024;
    _S1438[int(15)] = _S1024;
    _S1438[int(0)] = _S1362;
    _S1438[int(1)] = _S1349;
    _S1438[int(2)] = _S1347;
    _S1438[int(3)] = _S1345;
    _S1438[int(4)] = _S1334;
    _S1438[int(5)] = _S1332;
    _S1438[int(6)] = _S1330;
    _S1438[int(7)] = _S1328;
    _S1438[int(8)] = _S1326;
    _S1438[int(9)] = _S1318;
    _S1438[int(10)] = _S1316;
    _S1438[int(11)] = _S1314;
    _S1438[int(12)] = _S1312;
    _S1438[int(13)] = _S1310;
    _S1438[int(14)] = _S1308;
    _S1438[int(15)] = _S1306;
    float3  _S1439 = _S1438[int(0)];
    float3  _S1440 = _S1438[int(1)];
    float3  _S1441 = _S1438[int(2)];
    float3  _S1442 = _S1438[int(3)];
    float3  _S1443 = _S1438[int(4)];
    float3  _S1444 = _S1438[int(5)];
    float3  _S1445 = _S1438[int(6)];
    float3  _S1446 = _S1438[int(7)];
    float3  _S1447 = _S1438[int(8)];
    float3  _S1448 = _S1438[int(9)];
    float3  _S1449 = _S1438[int(10)];
    float3  _S1450 = _S1438[int(11)];
    float3  _S1451 = _S1438[int(12)];
    float3  _S1452 = _S1438[int(13)];
    float3  _S1453 = _S1438[int(14)];
    float3  _S1454 = _S1438[int(15)];
    float _S1455 = _S1419.differential_0 + _S1421.differential_0;
    float _S1456 = _S1418.differential_0 + _S1420.differential_0;
    float _S1457 = _S1422.differential_0 + _S1424.differential_0;
    float _S1458 = _S1423.differential_0 + _S1425.differential_0;
    float _S1459 = _S1415.differential_0 + _S1417.differential_0;
    float _S1460 = _S1379.differential_0 + _S1381.differential_0;
    float _S1461 = _S1383.differential_0 + _S1385.differential_0;
    float _S1462 = _S1387.differential_0 + _S1389.differential_0;
    float _S1463 = _S1391.differential_0 + _S1393.differential_0;
    float _S1464 = _S1395.differential_0 + _S1397.differential_0;
    float _S1465 = _S1399.differential_0 + _S1401.differential_0;
    float _S1466 = _S1403.differential_0 + _S1405.differential_0;
    float _S1467 = _S1407.differential_0 + _S1409.differential_0;
    float _S1468 = _S1411.differential_0 + _S1413.differential_0;
    float2  _S1469 = _S1229 * _S1437;
    float2  _S1470 = _S1240 * _S1437;
    float _S1471 = _S1469.x + _S1469.y;
    if(_S1233)
    {
        float _S1472 = _S1471 / _S1234;
        float _S1473 = _S1235 * - _S1472;
        float _S1474 = _S1232 * (0.3333333432674408f * - (_S1231 * _S1472));
        k_2 = _S1474 + _S1474;
        _S1234 = _S1473;
        _S1235 = 0.0f;
    }
    else
    {
        float _S1475 = _S1471 / _S1236;
        float _S1476 = _S1232 * - _S1475;
        k_2 = _S1230 * _S1475;
        _S1234 = 0.0f;
        _S1235 = _S1476;
    }
    DiffPair_float_0 _S1477;
    (&_S1477)->primal_0 = _S1230;
    (&_S1477)->differential_0 = 0.0f;
    DiffPair_float_0 _S1478;
    (&_S1478)->primal_0 = _S1231;
    (&_S1478)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1477, &_S1478, k_2);
    float _S1479 = _S1478.differential_0 + _S1234;
    float _S1480 = _S1477.differential_0 + _S1235;
    float2  _S1481 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1482;
    (&_S1482)->primal_0 = _S1229;
    (&_S1482)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1482, _S1480);
    float2  _S1483 = _S1482.differential_0 + _S1470;
    float _S1484 = fx_5 * _S1461;
    float2  _S1485 = make_float2 (_S1484, fy_5 * _S1460) + make_float2 (dist_coeffs_5[int(8)] * _S1484, dist_coeffs_5[int(9)] * _S1484);
    float2  _S1486 = _S1217 * _S1485;
    float _S1487 = dist_coeffs_5[int(4)] * _S1485.y;
    float _S1488 = dist_coeffs_5[int(5)] * _S1485.x;
    float _S1489 = _S1486.x + _S1486.y;
    float _S1490 = r2_62 * _S1489;
    float _S1491 = r2_62 * _S1490;
    float _S1492 = dist_coeffs_5[int(7)] * _S1485.y + _S1487 + dist_coeffs_5[int(6)] * _S1485.x + _S1488 + _S1220 * _S1489 + _S1219 * _S1490 + _S1218 * _S1491 + dist_coeffs_5[int(3)] * (r2_62 * _S1491);
    float _S1493 = v_62 * _S1492;
    float _S1494 = u_62 * _S1492;
    float2  _S1495 = make_float2 (radial_30) * _S1485 + make_float2 (_S1078 * (v_62 * _S1485.y) + _S1222 * _S1488 + 2.0f * (u_62 * _S1488) + _S1075 * (v_62 * _S1485.x) + _S1494 + _S1494, _S1224 * _S1487 + 2.0f * (v_62 * _S1487) + _S1223 * _S1485.y + _S1221 * _S1485.x + _S1493 + _S1493);
    float3  _S1496 = make_float3 (_S1483.x, _S1483.y, _S1479);
    float2  _S1497 = _S1205 * _S1495;
    float2  _S1498 = _S1216 * _S1495;
    float _S1499 = _S1497.x + _S1497.y;
    if(_S1209)
    {
        float _S1500 = _S1499 / _S1210;
        float _S1501 = _S1211 * - _S1500;
        float _S1502 = _S1208 * (0.3333333432674408f * - (_S1207 * _S1500));
        k_2 = _S1502 + _S1502;
        _S1210 = _S1501;
        _S1211 = 0.0f;
    }
    else
    {
        float _S1503 = _S1499 / _S1212;
        float _S1504 = _S1208 * - _S1503;
        k_2 = _S1206 * _S1503;
        _S1210 = 0.0f;
        _S1211 = _S1504;
    }
    DiffPair_float_0 _S1505;
    (&_S1505)->primal_0 = _S1206;
    (&_S1505)->differential_0 = 0.0f;
    DiffPair_float_0 _S1506;
    (&_S1506)->primal_0 = _S1207;
    (&_S1506)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1505, &_S1506, k_2);
    float _S1507 = _S1506.differential_0 + _S1210;
    float _S1508 = _S1505.differential_0 + _S1211;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1509;
    (&_S1509)->primal_0 = _S1205;
    (&_S1509)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1509, _S1508);
    float2  _S1510 = _S1509.differential_0 + _S1498;
    float _S1511 = fx_5 * _S1463;
    float2  _S1512 = make_float2 (_S1511, fy_5 * _S1462) + make_float2 (dist_coeffs_5[int(8)] * _S1511, dist_coeffs_5[int(9)] * _S1511);
    float2  _S1513 = _S1193 * _S1512;
    float _S1514 = dist_coeffs_5[int(4)] * _S1512.y;
    float _S1515 = dist_coeffs_5[int(5)] * _S1512.x;
    float _S1516 = _S1513.x + _S1513.y;
    float _S1517 = r2_61 * _S1516;
    float _S1518 = r2_61 * _S1517;
    float _S1519 = dist_coeffs_5[int(7)] * _S1512.y + _S1514 + dist_coeffs_5[int(6)] * _S1512.x + _S1515 + _S1196 * _S1516 + _S1195 * _S1517 + _S1194 * _S1518 + dist_coeffs_5[int(3)] * (r2_61 * _S1518);
    float _S1520 = v_61 * _S1519;
    float _S1521 = u_61 * _S1519;
    float2  _S1522 = make_float2 (radial_29) * _S1512 + make_float2 (_S1078 * (v_61 * _S1512.y) + _S1198 * _S1515 + 2.0f * (u_61 * _S1515) + _S1075 * (v_61 * _S1512.x) + _S1521 + _S1521, _S1200 * _S1514 + 2.0f * (v_61 * _S1514) + _S1199 * _S1512.y + _S1197 * _S1512.x + _S1520 + _S1520);
    float3  _S1523 = make_float3 (_S1510.x, _S1510.y, _S1507);
    float2  _S1524 = _S1181 * _S1522;
    float2  _S1525 = _S1192 * _S1522;
    float _S1526 = _S1524.x + _S1524.y;
    if(_S1185)
    {
        float _S1527 = _S1526 / _S1186;
        float _S1528 = _S1187 * - _S1527;
        float _S1529 = _S1184 * (0.3333333432674408f * - (_S1183 * _S1527));
        k_2 = _S1529 + _S1529;
        _S1186 = _S1528;
        _S1187 = 0.0f;
    }
    else
    {
        float _S1530 = _S1526 / _S1188;
        float _S1531 = _S1184 * - _S1530;
        k_2 = _S1182 * _S1530;
        _S1186 = 0.0f;
        _S1187 = _S1531;
    }
    DiffPair_float_0 _S1532;
    (&_S1532)->primal_0 = _S1182;
    (&_S1532)->differential_0 = 0.0f;
    DiffPair_float_0 _S1533;
    (&_S1533)->primal_0 = _S1183;
    (&_S1533)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1532, &_S1533, k_2);
    float _S1534 = _S1533.differential_0 + _S1186;
    float _S1535 = _S1532.differential_0 + _S1187;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1536;
    (&_S1536)->primal_0 = _S1181;
    (&_S1536)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1536, _S1535);
    float2  _S1537 = _S1536.differential_0 + _S1525;
    float _S1538 = fx_5 * _S1465;
    float2  _S1539 = make_float2 (_S1538, fy_5 * _S1464) + make_float2 (dist_coeffs_5[int(8)] * _S1538, dist_coeffs_5[int(9)] * _S1538);
    float2  _S1540 = _S1169 * _S1539;
    float _S1541 = dist_coeffs_5[int(4)] * _S1539.y;
    float _S1542 = dist_coeffs_5[int(5)] * _S1539.x;
    float _S1543 = _S1540.x + _S1540.y;
    float _S1544 = r2_60 * _S1543;
    float _S1545 = r2_60 * _S1544;
    float _S1546 = dist_coeffs_5[int(7)] * _S1539.y + _S1541 + dist_coeffs_5[int(6)] * _S1539.x + _S1542 + _S1172 * _S1543 + _S1171 * _S1544 + _S1170 * _S1545 + dist_coeffs_5[int(3)] * (r2_60 * _S1545);
    float _S1547 = v_60 * _S1546;
    float _S1548 = u_60 * _S1546;
    float2  _S1549 = make_float2 (radial_28) * _S1539 + make_float2 (_S1078 * (v_60 * _S1539.y) + _S1174 * _S1542 + 2.0f * (u_60 * _S1542) + _S1075 * (v_60 * _S1539.x) + _S1548 + _S1548, _S1176 * _S1541 + 2.0f * (v_60 * _S1541) + _S1175 * _S1539.y + _S1173 * _S1539.x + _S1547 + _S1547);
    float3  _S1550 = make_float3 (_S1537.x, _S1537.y, _S1534);
    float2  _S1551 = _S1157 * _S1549;
    float2  _S1552 = _S1168 * _S1549;
    float _S1553 = _S1551.x + _S1551.y;
    if(_S1161)
    {
        float _S1554 = _S1553 / _S1162;
        float _S1555 = _S1163 * - _S1554;
        float _S1556 = _S1160 * (0.3333333432674408f * - (_S1159 * _S1554));
        k_2 = _S1556 + _S1556;
        _S1162 = _S1555;
        _S1163 = 0.0f;
    }
    else
    {
        float _S1557 = _S1553 / _S1164;
        float _S1558 = _S1160 * - _S1557;
        k_2 = _S1158 * _S1557;
        _S1162 = 0.0f;
        _S1163 = _S1558;
    }
    DiffPair_float_0 _S1559;
    (&_S1559)->primal_0 = _S1158;
    (&_S1559)->differential_0 = 0.0f;
    DiffPair_float_0 _S1560;
    (&_S1560)->primal_0 = _S1159;
    (&_S1560)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1559, &_S1560, k_2);
    float _S1561 = _S1560.differential_0 + _S1162;
    float _S1562 = _S1559.differential_0 + _S1163;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1563;
    (&_S1563)->primal_0 = _S1157;
    (&_S1563)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1563, _S1562);
    float2  _S1564 = _S1563.differential_0 + _S1552;
    float _S1565 = fx_5 * _S1467;
    float2  _S1566 = make_float2 (_S1565, fy_5 * _S1466) + make_float2 (dist_coeffs_5[int(8)] * _S1565, dist_coeffs_5[int(9)] * _S1565);
    float2  _S1567 = _S1145 * _S1566;
    float _S1568 = dist_coeffs_5[int(4)] * _S1566.y;
    float _S1569 = dist_coeffs_5[int(5)] * _S1566.x;
    float _S1570 = _S1567.x + _S1567.y;
    float _S1571 = r2_59 * _S1570;
    float _S1572 = r2_59 * _S1571;
    float _S1573 = dist_coeffs_5[int(7)] * _S1566.y + _S1568 + dist_coeffs_5[int(6)] * _S1566.x + _S1569 + _S1148 * _S1570 + _S1147 * _S1571 + _S1146 * _S1572 + dist_coeffs_5[int(3)] * (r2_59 * _S1572);
    float _S1574 = v_59 * _S1573;
    float _S1575 = u_59 * _S1573;
    float2  _S1576 = make_float2 (radial_27) * _S1566 + make_float2 (_S1078 * (v_59 * _S1566.y) + _S1150 * _S1569 + 2.0f * (u_59 * _S1569) + _S1075 * (v_59 * _S1566.x) + _S1575 + _S1575, _S1152 * _S1568 + 2.0f * (v_59 * _S1568) + _S1151 * _S1566.y + _S1149 * _S1566.x + _S1574 + _S1574);
    float3  _S1577 = make_float3 (_S1564.x, _S1564.y, _S1561);
    float2  _S1578 = _S1133 * _S1576;
    float2  _S1579 = _S1144 * _S1576;
    float _S1580 = _S1578.x + _S1578.y;
    if(_S1137)
    {
        float _S1581 = _S1580 / _S1138;
        float _S1582 = _S1139 * - _S1581;
        float _S1583 = _S1136 * (0.3333333432674408f * - (_S1135 * _S1581));
        k_2 = _S1583 + _S1583;
        _S1138 = _S1582;
        _S1139 = 0.0f;
    }
    else
    {
        float _S1584 = _S1580 / _S1140;
        float _S1585 = _S1136 * - _S1584;
        k_2 = _S1134 * _S1584;
        _S1138 = 0.0f;
        _S1139 = _S1585;
    }
    DiffPair_float_0 _S1586;
    (&_S1586)->primal_0 = _S1134;
    (&_S1586)->differential_0 = 0.0f;
    DiffPair_float_0 _S1587;
    (&_S1587)->primal_0 = _S1135;
    (&_S1587)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1586, &_S1587, k_2);
    float _S1588 = _S1587.differential_0 + _S1138;
    float _S1589 = _S1586.differential_0 + _S1139;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1590;
    (&_S1590)->primal_0 = _S1133;
    (&_S1590)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1590, _S1589);
    float2  _S1591 = _S1590.differential_0 + _S1579;
    float _S1592 = fx_5 * _S1459;
    float2  _S1593 = make_float2 (_S1592, fy_5 * _S1468) + make_float2 (dist_coeffs_5[int(8)] * _S1592, dist_coeffs_5[int(9)] * _S1592);
    float2  _S1594 = _S1121 * _S1593;
    float _S1595 = dist_coeffs_5[int(4)] * _S1593.y;
    float _S1596 = dist_coeffs_5[int(5)] * _S1593.x;
    float _S1597 = _S1594.x + _S1594.y;
    float _S1598 = r2_58 * _S1597;
    float _S1599 = r2_58 * _S1598;
    float _S1600 = dist_coeffs_5[int(7)] * _S1593.y + _S1595 + dist_coeffs_5[int(6)] * _S1593.x + _S1596 + _S1124 * _S1597 + _S1123 * _S1598 + _S1122 * _S1599 + dist_coeffs_5[int(3)] * (r2_58 * _S1599);
    float _S1601 = v_58 * _S1600;
    float _S1602 = u_58 * _S1600;
    float2  _S1603 = make_float2 (radial_26) * _S1593 + make_float2 (_S1078 * (v_58 * _S1593.y) + _S1126 * _S1596 + 2.0f * (u_58 * _S1596) + _S1075 * (v_58 * _S1593.x) + _S1602 + _S1602, _S1128 * _S1595 + 2.0f * (v_58 * _S1595) + _S1127 * _S1593.y + _S1125 * _S1593.x + _S1601 + _S1601);
    float3  _S1604 = make_float3 (_S1591.x, _S1591.y, _S1588);
    float2  _S1605 = _S1109 * _S1603;
    float2  _S1606 = _S1120 * _S1603;
    float _S1607 = _S1605.x + _S1605.y;
    if(_S1113)
    {
        float _S1608 = _S1607 / _S1114;
        float _S1609 = _S1115 * - _S1608;
        float _S1610 = _S1112 * (0.3333333432674408f * - (_S1111 * _S1608));
        k_2 = _S1610 + _S1610;
        _S1114 = _S1609;
        _S1115 = 0.0f;
    }
    else
    {
        float _S1611 = _S1607 / _S1116;
        float _S1612 = _S1112 * - _S1611;
        k_2 = _S1110 * _S1611;
        _S1114 = 0.0f;
        _S1115 = _S1612;
    }
    DiffPair_float_0 _S1613;
    (&_S1613)->primal_0 = _S1110;
    (&_S1613)->differential_0 = 0.0f;
    DiffPair_float_0 _S1614;
    (&_S1614)->primal_0 = _S1111;
    (&_S1614)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1613, &_S1614, k_2);
    float _S1615 = _S1614.differential_0 + _S1114;
    float _S1616 = _S1613.differential_0 + _S1115;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1617;
    (&_S1617)->primal_0 = _S1109;
    (&_S1617)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1617, _S1616);
    float2  _S1618 = _S1617.differential_0 + _S1606;
    float _S1619 = fx_5 * _S1458;
    float2  _S1620 = make_float2 (_S1619, fy_5 * _S1455) + make_float2 (dist_coeffs_5[int(8)] * _S1619, dist_coeffs_5[int(9)] * _S1619);
    float2  _S1621 = _S1097 * _S1620;
    float _S1622 = dist_coeffs_5[int(4)] * _S1620.y;
    float _S1623 = dist_coeffs_5[int(5)] * _S1620.x;
    float _S1624 = _S1621.x + _S1621.y;
    float _S1625 = r2_57 * _S1624;
    float _S1626 = r2_57 * _S1625;
    float _S1627 = dist_coeffs_5[int(7)] * _S1620.y + _S1622 + dist_coeffs_5[int(6)] * _S1620.x + _S1623 + _S1100 * _S1624 + _S1099 * _S1625 + _S1098 * _S1626 + dist_coeffs_5[int(3)] * (r2_57 * _S1626);
    float _S1628 = v_57 * _S1627;
    float _S1629 = u_57 * _S1627;
    float2  _S1630 = make_float2 (radial_25) * _S1620 + make_float2 (_S1078 * (v_57 * _S1620.y) + _S1102 * _S1623 + 2.0f * (u_57 * _S1623) + _S1075 * (v_57 * _S1620.x) + _S1629 + _S1629, _S1104 * _S1622 + 2.0f * (v_57 * _S1622) + _S1103 * _S1620.y + _S1101 * _S1620.x + _S1628 + _S1628);
    float3  _S1631 = make_float3 (_S1618.x, _S1618.y, _S1615);
    float2  _S1632 = _S1085 * _S1630;
    float2  _S1633 = _S1096 * _S1630;
    float _S1634 = _S1632.x + _S1632.y;
    if(_S1089)
    {
        float _S1635 = _S1634 / _S1090;
        float _S1636 = _S1091 * - _S1635;
        float _S1637 = _S1088 * (0.3333333432674408f * - (_S1087 * _S1635));
        k_2 = _S1637 + _S1637;
        _S1090 = _S1636;
        _S1091 = 0.0f;
    }
    else
    {
        float _S1638 = _S1634 / _S1092;
        float _S1639 = _S1088 * - _S1638;
        k_2 = _S1086 * _S1638;
        _S1090 = 0.0f;
        _S1091 = _S1639;
    }
    DiffPair_float_0 _S1640;
    (&_S1640)->primal_0 = _S1086;
    (&_S1640)->differential_0 = 0.0f;
    DiffPair_float_0 _S1641;
    (&_S1641)->primal_0 = _S1087;
    (&_S1641)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1640, &_S1641, k_2);
    float _S1642 = _S1641.differential_0 + _S1090;
    float _S1643 = _S1640.differential_0 + _S1091;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1644;
    (&_S1644)->primal_0 = _S1085;
    (&_S1644)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1644, _S1643);
    float2  _S1645 = _S1644.differential_0 + _S1633;
    float _S1646 = fx_5 * _S1457;
    float2  _S1647 = make_float2 (_S1646, fy_5 * _S1456) + make_float2 (dist_coeffs_5[int(8)] * _S1646, dist_coeffs_5[int(9)] * _S1646);
    float2  _S1648 = _S1071 * _S1647;
    float _S1649 = dist_coeffs_5[int(4)] * _S1647.y;
    float _S1650 = dist_coeffs_5[int(5)] * _S1647.x;
    float _S1651 = _S1648.x + _S1648.y;
    float _S1652 = r2_56 * _S1651;
    float _S1653 = r2_56 * _S1652;
    float _S1654 = dist_coeffs_5[int(7)] * _S1647.y + _S1649 + dist_coeffs_5[int(6)] * _S1647.x + _S1650 + _S1074 * _S1651 + _S1073 * _S1652 + _S1072 * _S1653 + dist_coeffs_5[int(3)] * (r2_56 * _S1653);
    float _S1655 = v_56 * _S1654;
    float _S1656 = u_56 * _S1654;
    float2  _S1657 = make_float2 (radial_24) * _S1647 + make_float2 (_S1078 * (v_56 * _S1647.y) + _S1077 * _S1650 + 2.0f * (u_56 * _S1650) + _S1075 * (v_56 * _S1647.x) + _S1656 + _S1656, _S1080 * _S1649 + 2.0f * (v_56 * _S1649) + _S1079 * _S1647.y + _S1076 * _S1647.x + _S1655 + _S1655);
    float3  _S1658 = make_float3 (_S1645.x, _S1645.y, _S1642);
    float2  _S1659 = _S1059 * _S1657;
    float2  _S1660 = _S1070 * _S1657;
    float _S1661 = _S1659.x + _S1659.y;
    if(_S1063)
    {
        float _S1662 = _S1661 / _S1064;
        float _S1663 = _S1065 * - _S1662;
        float _S1664 = _S1062 * (0.3333333432674408f * - (_S1061 * _S1662));
        k_2 = _S1664 + _S1664;
        _S1064 = _S1663;
        _S1065 = 0.0f;
    }
    else
    {
        float _S1665 = _S1661 / _S1066;
        float _S1666 = _S1062 * - _S1665;
        k_2 = _S1060 * _S1665;
        _S1064 = 0.0f;
        _S1065 = _S1666;
    }
    DiffPair_float_0 _S1667;
    (&_S1667)->primal_0 = _S1060;
    (&_S1667)->differential_0 = 0.0f;
    DiffPair_float_0 _S1668;
    (&_S1668)->primal_0 = _S1061;
    (&_S1668)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1667, &_S1668, k_2);
    float _S1669 = _S1668.differential_0 + _S1064;
    float _S1670 = _S1667.differential_0 + _S1065;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1671;
    (&_S1671)->primal_0 = _S1059;
    (&_S1671)->differential_0 = _S1481;
    s_bwd_length_impl_0(&_S1671, _S1670);
    float2  _S1672 = _S1671.differential_0 + _S1660;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1673;
    (&_S1673)->primal_0 = R_5;
    (&_S1673)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1674;
    (&_S1674)->primal_0 = _S1058;
    (&_S1674)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1673, &_S1674, _S1363);
    DiffPair_float_0 _S1675;
    (&_S1675)->primal_0 = _S1055;
    (&_S1675)->differential_0 = 0.0f;
    DiffPair_float_0 _S1676;
    (&_S1676)->primal_0 = _S1057;
    (&_S1676)->differential_0 = 0.0f;
    _d_max_0(&_S1675, &_S1676, 0.0f);
    DiffPair_float_0 _S1677;
    (&_S1677)->primal_0 = _S1054;
    (&_S1677)->differential_0 = 0.0f;
    DiffPair_float_0 _S1678;
    (&_S1678)->primal_0 = _S1057;
    (&_S1678)->differential_0 = 0.0f;
    _d_min_0(&_S1677, &_S1678, 0.0f);
    float _S1679 = _S1676.differential_0 + _S1678.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1680;
    (&_S1680)->primal_0 = _S1056;
    (&_S1680)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1680, _S1679);
    float3  _S1681 = _S1680.differential_0 + _S1496;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1682;
    (&_S1682)->primal_0 = R_5;
    (&_S1682)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1683;
    (&_S1683)->primal_0 = pos_i_13;
    (&_S1683)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1682, &_S1683, _S1681);
    DiffPair_float_0 _S1684;
    (&_S1684)->primal_0 = _S1051;
    (&_S1684)->differential_0 = 0.0f;
    DiffPair_float_0 _S1685;
    (&_S1685)->primal_0 = _S1053;
    (&_S1685)->differential_0 = 0.0f;
    _d_max_0(&_S1684, &_S1685, _S1675.differential_0);
    DiffPair_float_0 _S1686;
    (&_S1686)->primal_0 = _S1050;
    (&_S1686)->differential_0 = 0.0f;
    DiffPair_float_0 _S1687;
    (&_S1687)->primal_0 = _S1053;
    (&_S1687)->differential_0 = 0.0f;
    _d_min_0(&_S1686, &_S1687, _S1677.differential_0);
    float _S1688 = _S1685.differential_0 + _S1687.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1689;
    (&_S1689)->primal_0 = _S1052;
    (&_S1689)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1689, _S1688);
    float3  _S1690 = _S1689.differential_0 + _S1523;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1691;
    (&_S1691)->primal_0 = R_5;
    (&_S1691)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1692;
    (&_S1692)->primal_0 = pos_i_12;
    (&_S1692)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1691, &_S1692, _S1690);
    DiffPair_float_0 _S1693;
    (&_S1693)->primal_0 = _S1047;
    (&_S1693)->differential_0 = 0.0f;
    DiffPair_float_0 _S1694;
    (&_S1694)->primal_0 = _S1049;
    (&_S1694)->differential_0 = 0.0f;
    _d_max_0(&_S1693, &_S1694, _S1684.differential_0);
    DiffPair_float_0 _S1695;
    (&_S1695)->primal_0 = _S1046;
    (&_S1695)->differential_0 = 0.0f;
    DiffPair_float_0 _S1696;
    (&_S1696)->primal_0 = _S1049;
    (&_S1696)->differential_0 = 0.0f;
    _d_min_0(&_S1695, &_S1696, _S1686.differential_0);
    float _S1697 = _S1694.differential_0 + _S1696.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1698;
    (&_S1698)->primal_0 = _S1048;
    (&_S1698)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1698, _S1697);
    float3  _S1699 = _S1698.differential_0 + _S1550;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1700;
    (&_S1700)->primal_0 = R_5;
    (&_S1700)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1701;
    (&_S1701)->primal_0 = pos_i_11;
    (&_S1701)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1700, &_S1701, _S1699);
    DiffPair_float_0 _S1702;
    (&_S1702)->primal_0 = _S1043;
    (&_S1702)->differential_0 = 0.0f;
    DiffPair_float_0 _S1703;
    (&_S1703)->primal_0 = _S1045;
    (&_S1703)->differential_0 = 0.0f;
    _d_max_0(&_S1702, &_S1703, _S1693.differential_0);
    DiffPair_float_0 _S1704;
    (&_S1704)->primal_0 = _S1042;
    (&_S1704)->differential_0 = 0.0f;
    DiffPair_float_0 _S1705;
    (&_S1705)->primal_0 = _S1045;
    (&_S1705)->differential_0 = 0.0f;
    _d_min_0(&_S1704, &_S1705, _S1695.differential_0);
    float _S1706 = _S1703.differential_0 + _S1705.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1707;
    (&_S1707)->primal_0 = _S1044;
    (&_S1707)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1707, _S1706);
    float3  _S1708 = _S1707.differential_0 + _S1577;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1709;
    (&_S1709)->primal_0 = R_5;
    (&_S1709)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1710;
    (&_S1710)->primal_0 = pos_i_10;
    (&_S1710)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1709, &_S1710, _S1708);
    DiffPair_float_0 _S1711;
    (&_S1711)->primal_0 = _S1039;
    (&_S1711)->differential_0 = 0.0f;
    DiffPair_float_0 _S1712;
    (&_S1712)->primal_0 = _S1041;
    (&_S1712)->differential_0 = 0.0f;
    _d_max_0(&_S1711, &_S1712, _S1702.differential_0);
    DiffPair_float_0 _S1713;
    (&_S1713)->primal_0 = _S1038;
    (&_S1713)->differential_0 = 0.0f;
    DiffPair_float_0 _S1714;
    (&_S1714)->primal_0 = _S1041;
    (&_S1714)->differential_0 = 0.0f;
    _d_min_0(&_S1713, &_S1714, _S1704.differential_0);
    float _S1715 = _S1712.differential_0 + _S1714.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1716;
    (&_S1716)->primal_0 = _S1040;
    (&_S1716)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1716, _S1715);
    float3  _S1717 = _S1716.differential_0 + _S1604;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1718;
    (&_S1718)->primal_0 = R_5;
    (&_S1718)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1719;
    (&_S1719)->primal_0 = pos_i_9;
    (&_S1719)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1718, &_S1719, _S1717);
    DiffPair_float_0 _S1720;
    (&_S1720)->primal_0 = _S1035;
    (&_S1720)->differential_0 = 0.0f;
    DiffPair_float_0 _S1721;
    (&_S1721)->primal_0 = _S1037;
    (&_S1721)->differential_0 = 0.0f;
    _d_max_0(&_S1720, &_S1721, _S1711.differential_0);
    DiffPair_float_0 _S1722;
    (&_S1722)->primal_0 = _S1034;
    (&_S1722)->differential_0 = 0.0f;
    DiffPair_float_0 _S1723;
    (&_S1723)->primal_0 = _S1037;
    (&_S1723)->differential_0 = 0.0f;
    _d_min_0(&_S1722, &_S1723, _S1713.differential_0);
    float _S1724 = _S1721.differential_0 + _S1723.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1725;
    (&_S1725)->primal_0 = _S1036;
    (&_S1725)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1725, _S1724);
    float3  _S1726 = _S1725.differential_0 + _S1631;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1727;
    (&_S1727)->primal_0 = R_5;
    (&_S1727)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1728;
    (&_S1728)->primal_0 = pos_i_8;
    (&_S1728)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1727, &_S1728, _S1726);
    DiffPair_float_0 _S1729;
    (&_S1729)->primal_0 = _S1031;
    (&_S1729)->differential_0 = 0.0f;
    DiffPair_float_0 _S1730;
    (&_S1730)->primal_0 = _S1033;
    (&_S1730)->differential_0 = 0.0f;
    _d_max_0(&_S1729, &_S1730, _S1720.differential_0);
    DiffPair_float_0 _S1731;
    (&_S1731)->primal_0 = _S1030;
    (&_S1731)->differential_0 = 0.0f;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = _S1033;
    (&_S1732)->differential_0 = 0.0f;
    _d_min_0(&_S1731, &_S1732, _S1722.differential_0);
    float _S1733 = _S1730.differential_0 + _S1732.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1734;
    (&_S1734)->primal_0 = _S1032;
    (&_S1734)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1734, _S1733);
    float3  _S1735 = _S1734.differential_0 + _S1658;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1736;
    (&_S1736)->primal_0 = R_5;
    (&_S1736)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1737;
    (&_S1737)->primal_0 = pos_i_7;
    (&_S1737)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1736, &_S1737, _S1735);
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = 0.0f;
    (&_S1738)->differential_0 = 0.0f;
    DiffPair_float_0 _S1739;
    (&_S1739)->primal_0 = _S1029;
    (&_S1739)->differential_0 = 0.0f;
    _d_max_0(&_S1738, &_S1739, _S1729.differential_0);
    DiffPair_float_0 _S1740;
    (&_S1740)->primal_0 = 1.00000001504746622e+30f;
    (&_S1740)->differential_0 = 0.0f;
    DiffPair_float_0 _S1741;
    (&_S1741)->primal_0 = _S1029;
    (&_S1741)->differential_0 = 0.0f;
    _d_min_0(&_S1740, &_S1741, _S1731.differential_0);
    float _S1742 = _S1739.differential_0 + _S1741.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1743;
    (&_S1743)->primal_0 = _S1028;
    (&_S1743)->differential_0 = _S1024;
    s_bwd_length_impl_1(&_S1743, _S1742);
    float3  _S1744 = _S1743.differential_0 + make_float3 (_S1672.x, _S1672.y, _S1669);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1745;
    (&_S1745)->primal_0 = R_5;
    (&_S1745)->differential_0 = _S1365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1746;
    (&_S1746)->primal_0 = pos_5;
    (&_S1746)->differential_0 = _S1024;
    s_bwd_prop_mul_0(&_S1745, &_S1746, _S1744);
    float3  _S1747 = _S1363 + _S1681 + _S1690 + _S1699 + _S1708 + _S1717 + _S1726 + _S1735 + _S1744 + _S1368.differential_0;
    Matrix<float, 3, 3>  _S1748 = _S1673.differential_0 + _S1682.differential_0 + _S1691.differential_0 + _S1700.differential_0 + _S1709.differential_0 + _S1718.differential_0 + _S1727.differential_0 + _S1736.differential_0 + _S1745.differential_0 + _S1369;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_1)[int(0)] = _S1439;
    (*v_sh_coeffs_1)[int(1)] = _S1440;
    (*v_sh_coeffs_1)[int(2)] = _S1441;
    (*v_sh_coeffs_1)[int(3)] = _S1442;
    (*v_sh_coeffs_1)[int(4)] = _S1443;
    (*v_sh_coeffs_1)[int(5)] = _S1444;
    (*v_sh_coeffs_1)[int(6)] = _S1445;
    (*v_sh_coeffs_1)[int(7)] = _S1446;
    (*v_sh_coeffs_1)[int(8)] = _S1447;
    (*v_sh_coeffs_1)[int(9)] = _S1448;
    (*v_sh_coeffs_1)[int(10)] = _S1449;
    (*v_sh_coeffs_1)[int(11)] = _S1450;
    (*v_sh_coeffs_1)[int(12)] = _S1451;
    (*v_sh_coeffs_1)[int(13)] = _S1452;
    (*v_sh_coeffs_1)[int(14)] = _S1453;
    (*v_sh_coeffs_1)[int(15)] = _S1454;
    *v_R_1 = _S1748;
    *v_t_1 = _S1747;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  dOut_7)
{
    float3  _S1749 = _slang_select(((*dpx_8).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_8).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_7;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S1749;
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
    float3  _S1750 = - (m_1 * (ray_o_0 - center_0));
    float3  ta_0 = _S1750 - k_3;
    float3  tb_0 = _S1750 + k_3;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S1751 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S1751;
    return (*t0_0) < _S1751;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_9, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_8)
{
    float _S1752 = (1.0f - (*dps_0).primal_0) * dOut_8;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S1752;
    DiffPair_float_0 _S1753 = *dpy_4;
    float _S1754 = (*dps_0).primal_0 * dOut_8;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S1754;
    float _S1755 = (_S1753.primal_0 - (*dpx_9).primal_0) * dOut_8;
    dps_0->primal_0 = _S1753.primal_0;
    dps_0->differential_0 = _S1755;
    return;
}

inline __device__ float lerp_0(float x_16, float y_11, float s_0)
{
    return x_16 + (y_11 - x_16) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_0)
{
    float _S1756 = w_0.z;
    float _S1757 = 1.0f - _S1756;
    float _S1758 = w_0.y;
    float _S1759 = 1.0f - _S1758;
    float _S1760 = _S1757 * _S1759;
    float _S1761 = w_0.x;
    float _S1762 = 1.0f - _S1761;
    float _S1763 = _S1757 * _S1758;
    float _S1764 = _S1756 * _S1759;
    float _S1765 = _S1756 * _S1758;
    return _S1760 * _S1762 * (*densities_6)[int(0)] + _S1760 * _S1761 * (*densities_6)[int(1)] + _S1763 * _S1762 * (*densities_6)[int(2)] + _S1763 * _S1761 * (*densities_6)[int(3)] + _S1764 * _S1762 * (*densities_6)[int(4)] + _S1764 * _S1761 * (*densities_6)[int(5)] + _S1765 * _S1762 * (*densities_6)[int(6)] + _S1765 * _S1761 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  densities_7, float3  ray_o_1, float3  ray_d_1)
{
    float _S1766 = 0.5f * size_6;
    float3  m_2 = make_float3 (1.0f) / ray_d_1;
    float3  k_4 = abs_0(m_2) * make_float3 (_S1766);
    float3  _S1767 = - (m_2 * (ray_o_1 - (pos_6 + make_float3 (_S1766))));
    float3  ta_1 = _S1767 - k_4;
    float3  tb_1 = _S1767 + k_4;
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
        float3  _S1768 = (ray_o_1 + ray_d_1 * make_float3 (lerp_0(t0_1, t1_1, (float(i_5) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6);
        FixedArray<float, 8>  _S1769 = densities_7;
        float _S1770 = interp_0(&_S1769, _S1768);
        float _S1771;
        if(_S1770 > 1.10000002384185791f)
        {
            _S1771 = _S1770;
        }
        else
        {
            _S1771 = (F32_exp((0.90909093618392944f * _S1770 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S1771;
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
    float _S1772;
};

inline __device__ float3  s_primal_ctx_abs_0(float3  _S1773)
{
    return abs_0(_S1773);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1774, float _S1775, float _S1776)
{
    return lerp_0(_S1774, _S1775, _S1776);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S1777 = dpw_0.z;
    float _S1778 = 1.0f - _S1777;
    float _S1779 = dpw_0.y;
    float _S1780 = 1.0f - _S1779;
    float _S1781 = _S1778 * _S1780;
    float _S1782 = dpw_0.x;
    float _S1783 = 1.0f - _S1782;
    float _S1784 = _S1778 * _S1779;
    float _S1785 = _S1777 * _S1780;
    float _S1786 = _S1777 * _S1779;
    return _S1781 * _S1783 * (*dpdensities_0)[int(0)] + _S1781 * _S1782 * (*dpdensities_0)[int(1)] + _S1784 * _S1783 * (*dpdensities_0)[int(2)] + _S1784 * _S1782 * (*dpdensities_0)[int(3)] + _S1785 * _S1783 * (*dpdensities_0)[int(4)] + _S1785 * _S1782 * (*dpdensities_0)[int(5)] + _S1786 * _S1783 * (*dpdensities_0)[int(6)] + _S1786 * _S1782 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_exp_0(float _S1787)
{
    return (F32_exp((_S1787)));
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_0, float3  dpray_d_0, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S1772 = 0.0f;
    _s_diff_ctx_0->_S1772 = 0.0f;
    float _S1788 = 0.5f * size_7;
    float3  m_3 = make_float3 (1.0f) / dpray_d_0;
    float3  k_5 = s_primal_ctx_abs_0(m_3) * make_float3 (_S1788);
    float3  _S1789 = - (m_3 * (dpray_o_0 - (pos_7 + make_float3 (_S1788))));
    float3  ta_2 = _S1789 - k_5;
    float3  tb_2 = _S1789 + k_5;
    float _S1790 = (F32_max(((F32_max((ta_2.x), (ta_2.y)))), ((F32_max((ta_2.z), (0.0f))))));
    float _S1791 = (F32_min(((F32_min((tb_2.x), (tb_2.y)))), (tb_2.z)));
    float accum_2;
    if(!!(_S1790 < _S1791))
    {
        float _S1792 = - (_S1791 - _S1790);
        bool _runFlag_0 = true;
        int i_6 = int(0);
        accum_2 = 0.0f;
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_0->_S1772 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            int _S1793;
            float _S1794;
            if(i_6 < int(8))
            {
                float _S1795 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_0 + dpray_d_0 * make_float3 (s_primal_ctx_lerp_0(_S1790, _S1791, (float(i_6) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S1795 > 1.10000002384185791f)
                {
                    _S1794 = _S1795;
                }
                else
                {
                    _S1794 = s_primal_ctx_exp_0(0.90909093618392944f * _S1795 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S1794;
                _S1793 = int(2);
                _S1794 = accum_3;
            }
            else
            {
                _S1793 = int(1);
                _S1794 = 0.0f;
            }
            if(_S1793 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                i_6 = i_6 + int(1);
                accum_2 = _S1794;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = (F32_min((1.0f - s_primal_ctx_exp_0(_S1792 / 8.0f * accum_2)), (0.99900001287460327f)));
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1796, float _S1797)
{
    _d_exp_0(_S1796, _S1797);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_2)
{
    float _S1798 = (*dpw_1).primal_0.z;
    float _S1799 = 1.0f - _S1798;
    float _S1800 = (*dpw_1).primal_0.y;
    float _S1801 = 1.0f - _S1800;
    float _S1802 = _S1799 * _S1801;
    float _S1803 = (*dpw_1).primal_0.x;
    float _S1804 = 1.0f - _S1803;
    float _S1805 = _S1799 * _S1800;
    float _S1806 = _S1798 * _S1801;
    float _S1807 = _S1798 * _S1800;
    float _S1808 = _S1807 * _S1803 * _s_dOut_2;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_2;
    float _S1809 = _S1807 * _S1804 * _s_dOut_2;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_2;
    float _S1810 = _S1806 * _S1803 * _s_dOut_2;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_2;
    float _S1811 = _S1806 * _S1804 * _s_dOut_2;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_2;
    float _S1812 = _S1805 * _S1803 * _s_dOut_2;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_2;
    float _S1813 = _S1805 * _S1804 * _s_dOut_2;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_2;
    float _S1814 = _S1802 * _S1803 * _s_dOut_2;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_2;
    float _S1815 = _S1802 * _S1804 * _s_dOut_2;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_2;
    float _S1816 = _S1803 * s_diff_w7_T_0 + _S1804 * s_diff_w6_T_0;
    float _S1817 = _S1803 * s_diff_w5_T_0 + _S1804 * s_diff_w4_T_0;
    float _S1818 = _S1803 * s_diff_w3_T_0 + _S1804 * s_diff_w2_T_0;
    float _S1819 = _S1803 * s_diff_w1_T_0 + _S1804 * s_diff_w0_T_0;
    float3  _S1820 = make_float3 (_S1807 * s_diff_w7_T_0 + _S1806 * s_diff_w5_T_0 + _S1805 * s_diff_w3_T_0 + _S1802 * s_diff_w1_T_0 + - (_S1807 * s_diff_w6_T_0 + _S1806 * s_diff_w4_T_0 + _S1805 * s_diff_w2_T_0 + _S1802 * s_diff_w0_T_0), _S1798 * _S1816 + _S1799 * _S1818 + - (_S1798 * _S1817 + _S1799 * _S1819), _S1800 * _S1816 + _S1801 * _S1817 + - (_S1800 * _S1818 + _S1801 * _S1819));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S1820;
    FixedArray<float, 8>  _S1821;
    _S1821[int(0)] = 0.0f;
    _S1821[int(1)] = 0.0f;
    _S1821[int(2)] = 0.0f;
    _S1821[int(3)] = 0.0f;
    _S1821[int(4)] = 0.0f;
    _S1821[int(5)] = 0.0f;
    _S1821[int(6)] = 0.0f;
    _S1821[int(7)] = 0.0f;
    _S1821[int(7)] = _S1808;
    _S1821[int(6)] = _S1809;
    _S1821[int(5)] = _S1810;
    _S1821[int(4)] = _S1811;
    _S1821[int(3)] = _S1812;
    _S1821[int(2)] = _S1813;
    _S1821[int(1)] = _S1814;
    _S1821[int(0)] = _S1815;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S1821;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S1822, DiffPair_float_0 * _S1823, DiffPair_float_0 * _S1824, float _S1825)
{
    _d_lerp_0(_S1822, _S1823, _S1824, _S1825);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1826, float3  _S1827)
{
    _d_abs_vector_0(_S1826, _S1827);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float _s_dOut_3, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_1)
{
    FixedArray<float, 8>  _S1828 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1829 = *dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1830 = *dpray_d_1;
    float _S1831 = 0.5f * size_8;
    float3  _S1832 = make_float3 (_S1831);
    float3  _S1833 = make_float3 (size_8);
    float3  m_4 = make_float3 (1.0f) / (*dpray_d_1).primal_0;
    float3  _S1834 = (*dpray_d_1).primal_0 * (*dpray_d_1).primal_0;
    float3  _S1835 = (*dpray_o_1).primal_0 - (pos_8 + make_float3 (_S1831));
    float3  k_6 = s_primal_ctx_abs_0(m_4) * make_float3 (_S1831);
    float3  _S1836 = - (m_4 * _S1835);
    float3  ta_3 = _S1836 - k_6;
    float3  tb_3 = _S1836 + k_6;
    float _S1837 = ta_3.x;
    float _S1838 = ta_3.y;
    float _S1839 = (F32_max((_S1837), (_S1838)));
    float _S1840 = ta_3.z;
    float _S1841 = (F32_max((_S1840), (0.0f)));
    float _S1842 = (F32_max((_S1839), (_S1841)));
    float _S1843 = tb_3.x;
    float _S1844 = tb_3.y;
    float _S1845 = (F32_min((_S1843), (_S1844)));
    float _S1846 = tb_3.z;
    float _S1847 = (F32_min((_S1845), (_S1846)));
    bool _S1848 = !!(_S1842 < _S1847);
    float _S1849;
    float _S1850;
    float _S1851;
    if(_S1848)
    {
        float _S1852 = - (_S1847 - _S1842) / 8.0f;
        float _S1853 = _S1852 * _s_diff_ctx_1->_S1772;
        _S1849 = 1.0f - s_primal_ctx_exp_0(_S1853);
        _S1850 = _S1853;
        _S1851 = _S1852;
    }
    else
    {
        _S1849 = 0.0f;
        _S1850 = 0.0f;
        _S1851 = 0.0f;
    }
    float3  _S1854 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1855 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<float, 8>  _S1856;
    float3  _S1857;
    float3  _S1858;
    if(_S1848)
    {
        DiffPair_float_0 _S1859;
        (&_S1859)->primal_0 = _S1849;
        (&_S1859)->differential_0 = 0.0f;
        DiffPair_float_0 _S1860;
        (&_S1860)->primal_0 = 0.99900001287460327f;
        (&_S1860)->differential_0 = 0.0f;
        _d_min_0(&_S1859, &_S1860, _s_dOut_3);
        float _S1861 = - _S1859.differential_0;
        DiffPair_float_0 _S1862;
        (&_S1862)->primal_0 = _S1850;
        (&_S1862)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S1862, _S1861);
        float _S1863 = _S1851 * _S1862.differential_0;
        float _S1864 = 0.125f * (_s_diff_ctx_1->_S1772 * _S1862.differential_0);
        int _dc_0 = int(8);
        _S1849 = _S1863;
        _S1856[int(0)] = 0.0f;
        _S1856[int(1)] = 0.0f;
        _S1856[int(2)] = 0.0f;
        _S1856[int(3)] = 0.0f;
        _S1856[int(4)] = 0.0f;
        _S1856[int(5)] = 0.0f;
        _S1856[int(6)] = 0.0f;
        _S1856[int(7)] = 0.0f;
        _S1857 = _S1854;
        _S1858 = _S1854;
        _S1850 = 0.0f;
        _S1851 = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S1865 = _dc_0 < int(8);
            float _S1866;
            float _S1867;
            int _S1868;
            float3  _S1869;
            float3  _S1870;
            bool _S1871;
            if(_S1865)
            {
                float _S1872 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S1873 = s_primal_ctx_lerp_0(_S1842, _S1847, _S1872);
                float3  _S1874 = make_float3 (_S1873);
                float3  _S1875 = (_S1829.primal_0 + _S1830.primal_0 * make_float3 (_S1873) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S1876 = _S1828;
                float _S1877 = s_primal_ctx_interp_0(&_S1876, _S1875);
                bool _S1878 = _S1877 > 1.10000002384185791f;
                if(_S1878)
                {
                    _S1866 = 0.0f;
                }
                else
                {
                    _S1866 = 0.90909093618392944f * _S1877 - 0.90468984842300415f;
                }
                _S1868 = int(2);
                _S1871 = _S1878;
                _S1869 = _S1875;
                _S1870 = _S1874;
                _S1867 = _S1872;
            }
            else
            {
                _S1868 = int(1);
                _S1871 = false;
                _S1866 = 0.0f;
                _S1869 = _S1854;
                _S1870 = _S1854;
                _S1867 = 0.0f;
            }
            float _S1879;
            float _S1880;
            if(!(_S1868 != int(2)))
            {
                _S1879 = _S1849;
                _S1880 = 0.0f;
            }
            else
            {
                _S1879 = 0.0f;
                _S1880 = _S1849;
            }
            if(_S1865)
            {
                float _S1881 = _S1879 + _S1880;
                float _S1882;
                if(_S1871)
                {
                    _S1882 = _S1879;
                }
                else
                {
                    DiffPair_float_0 _S1883;
                    (&_S1883)->primal_0 = _S1866;
                    (&_S1883)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S1883, _S1879);
                    _S1882 = 0.90909093618392944f * _S1883.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S1884;
                (&_S1884)->primal_0 = _S1828;
                (&_S1884)->differential_0 = _S1855;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1885;
                (&_S1885)->primal_0 = _S1869;
                (&_S1885)->differential_0 = _S1854;
                s_bwd_prop_interp_0(&_S1884, &_S1885, _S1882);
                float3  _S1886 = _S1885.differential_0 / _S1833;
                float3  _S1887 = _S1830.primal_0 * _S1886;
                float3  _S1888 = _S1870 * _S1886;
                float _S1889 = _S1887.x + _S1887.y + _S1887.z;
                DiffPair_float_0 _S1890;
                (&_S1890)->primal_0 = _S1842;
                (&_S1890)->differential_0 = 0.0f;
                DiffPair_float_0 _S1891;
                (&_S1891)->primal_0 = _S1847;
                (&_S1891)->differential_0 = 0.0f;
                DiffPair_float_0 _S1892;
                (&_S1892)->primal_0 = _S1867;
                (&_S1892)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S1890, &_S1891, &_S1892, _S1889);
                float _S1893 = (&_S1884)->differential_0[int(0)] + _S1856[int(0)];
                float _S1894 = (&_S1884)->differential_0[int(1)] + _S1856[int(1)];
                float _S1895 = (&_S1884)->differential_0[int(2)] + _S1856[int(2)];
                float _S1896 = (&_S1884)->differential_0[int(3)] + _S1856[int(3)];
                float _S1897 = (&_S1884)->differential_0[int(4)] + _S1856[int(4)];
                float _S1898 = (&_S1884)->differential_0[int(5)] + _S1856[int(5)];
                float _S1899 = (&_S1884)->differential_0[int(6)] + _S1856[int(6)];
                float _S1900 = (&_S1884)->differential_0[int(7)] + _S1856[int(7)];
                float3  _S1901 = _S1886 + _S1857;
                float3  _S1902 = _S1888 + _S1858;
                float _S1903 = _S1891.differential_0 + _S1850;
                float _S1904 = _S1890.differential_0 + _S1851;
                _S1849 = _S1881;
                _S1856[int(0)] = _S1893;
                _S1856[int(1)] = _S1894;
                _S1856[int(2)] = _S1895;
                _S1856[int(3)] = _S1896;
                _S1856[int(4)] = _S1897;
                _S1856[int(5)] = _S1898;
                _S1856[int(6)] = _S1899;
                _S1856[int(7)] = _S1900;
                _S1857 = _S1901;
                _S1858 = _S1902;
                _S1850 = _S1903;
                _S1851 = _S1904;
            }
            else
            {
                _S1849 = _S1880;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S1905 = - _S1864;
        float _S1906 = - _S1905 + _S1851;
        float3  _S1907 = _S1857;
        _S1849 = _S1905 + _S1850;
        _S1850 = _S1906;
        _S1857 = _S1858;
        _S1858 = _S1907;
    }
    else
    {
        _S1849 = 0.0f;
        _S1850 = 0.0f;
        _S1857 = _S1854;
        _S1858 = _S1854;
        _S1856[int(0)] = 0.0f;
        _S1856[int(1)] = 0.0f;
        _S1856[int(2)] = 0.0f;
        _S1856[int(3)] = 0.0f;
        _S1856[int(4)] = 0.0f;
        _S1856[int(5)] = 0.0f;
        _S1856[int(6)] = 0.0f;
        _S1856[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S1908;
    (&_S1908)->primal_0 = _S1845;
    (&_S1908)->differential_0 = 0.0f;
    DiffPair_float_0 _S1909;
    (&_S1909)->primal_0 = _S1846;
    (&_S1909)->differential_0 = 0.0f;
    _d_min_0(&_S1908, &_S1909, _S1849);
    DiffPair_float_0 _S1910;
    (&_S1910)->primal_0 = _S1843;
    (&_S1910)->differential_0 = 0.0f;
    DiffPair_float_0 _S1911;
    (&_S1911)->primal_0 = _S1844;
    (&_S1911)->differential_0 = 0.0f;
    _d_min_0(&_S1910, &_S1911, _S1908.differential_0);
    DiffPair_float_0 _S1912;
    (&_S1912)->primal_0 = _S1839;
    (&_S1912)->differential_0 = 0.0f;
    DiffPair_float_0 _S1913;
    (&_S1913)->primal_0 = _S1841;
    (&_S1913)->differential_0 = 0.0f;
    _d_max_0(&_S1912, &_S1913, _S1850);
    DiffPair_float_0 _S1914;
    (&_S1914)->primal_0 = _S1840;
    (&_S1914)->differential_0 = 0.0f;
    DiffPair_float_0 _S1915;
    (&_S1915)->primal_0 = 0.0f;
    (&_S1915)->differential_0 = 0.0f;
    _d_max_0(&_S1914, &_S1915, _S1913.differential_0);
    DiffPair_float_0 _S1916;
    (&_S1916)->primal_0 = _S1837;
    (&_S1916)->differential_0 = 0.0f;
    DiffPair_float_0 _S1917;
    (&_S1917)->primal_0 = _S1838;
    (&_S1917)->differential_0 = 0.0f;
    _d_max_0(&_S1916, &_S1917, _S1912.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S1910.differential_0, _S1911.differential_0, _S1909.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S1916.differential_0, _S1917.differential_0, _S1914.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S1918 = _S1832 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1919;
    (&_S1919)->primal_0 = m_4;
    (&_S1919)->differential_0 = _S1854;
    s_bwd_prop_abs_0(&_S1919, _S1918);
    float3  _S1920 = m_4 * s_diff_n_T_0;
    float3  _S1921 = - ((_S1919.differential_0 + _S1835 * s_diff_n_T_0) / _S1834) + _S1857;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S1921;
    float3  _S1922 = _S1920 + _S1858;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S1922;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S1856;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S1923, float _S1924, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S1925, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1926, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1927, float _S1928)
{
    FixedArray<float, 8>  _S1929 = _S1925->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1930;
    float _S1931 = s_primal_ctx_evaluate_alpha_voxel_0(_S1923, _S1924, &_S1929, (*_S1926).primal_0, (*_S1927).primal_0, &_S1930);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S1932 = _S1930;
    s_bwd_prop_evaluate_alpha_voxel_0(_S1923, _S1924, _S1925, _S1926, _S1927, _S1928, &_S1932);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  densities_8, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_0, float3  * v_ray_d_0)
{
    FixedArray<float, 8>  _S1933 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = densities_8;
    (&dp_densities_0)->differential_0 = _S1933;
    float3  _S1934 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1934;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1934;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_0 = dp_ray_o_0.differential_0;
    *v_ray_d_0 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  densities_9, float3  rgb_0, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_4)
{
    *out_rgb_0 = rgb_0;
    float _S1935 = 0.5f * size_10;
    float3  m_5 = make_float3 (1.0f) / ray_d_3;
    float3  k_7 = abs_0(m_5) * make_float3 (_S1935);
    float3  _S1936 = - (m_5 * (ray_o_3 - (pos_10 + make_float3 (_S1935))));
    float3  ta_4 = _S1936 - k_7;
    float3  tb_4 = _S1936 + k_7;
    float _S1937 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S1938 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
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
        float t_6 = lerp_0(_S1937, _S1938, (float(i_7) + 0.5f) / 8.0f);
        float3  _S1939 = (ray_o_3 + ray_d_3 * make_float3 (t_6) - pos_10) / make_float3 (size_10);
        FixedArray<float, 8>  _S1940 = densities_9;
        float _S1941 = interp_0(&_S1940, _S1939);
        float _S1942;
        if(_S1941 > 1.10000002384185791f)
        {
            _S1942 = _S1941;
        }
        else
        {
            _S1942 = (F32_exp((0.90909093618392944f * _S1941 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S1942;
        float depth_accum_1 = depth_accum_0 + t_6 * _S1942;
        i_7 = i_7 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_4 = (F32_max((depth_accum_0 / accum_4), (0.0f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S1943;
    float _S1944;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_0, float3  dpray_o_2, float3  dpray_d_2, float3  * dpout_rgb_0, float * dpdepth_0, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_2)
{
    _s_diff_ctx_2->_S1943 = 0.0f;
    _s_diff_ctx_2->_S1944 = 0.0f;
    float _S1945 = 0.5f * size_11;
    float3  m_6 = make_float3 (1.0f) / dpray_d_2;
    float3  k_8 = s_primal_ctx_abs_0(m_6) * make_float3 (_S1945);
    float3  _S1946 = - (m_6 * (dpray_o_2 - (pos_11 + make_float3 (_S1945))));
    float3  ta_5 = _S1946 - k_8;
    float3  tb_5 = _S1946 + k_8;
    float _S1947 = (F32_max(((F32_max((ta_5.x), (ta_5.y)))), ((F32_max((ta_5.z), (0.0f))))));
    float _S1948 = (F32_min(((F32_min((tb_5.x), (tb_5.y)))), (tb_5.z)));
    bool _runFlag_1 = true;
    int i_8 = int(0);
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_2->_S1943 = depth_accum_2;
        _s_diff_ctx_2->_S1944 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        int _S1949;
        float _S1950;
        float _S1951;
        if(i_8 < int(8))
        {
            float _S1952 = s_primal_ctx_lerp_0(_S1947, _S1948, (float(i_8) + 0.5f) / 8.0f);
            float _S1953 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_2 + dpray_d_2 * make_float3 (_S1952) - pos_11) / make_float3 (size_11));
            if(_S1953 > 1.10000002384185791f)
            {
                _S1950 = _S1953;
            }
            else
            {
                _S1950 = s_primal_ctx_exp_0(0.90909093618392944f * _S1953 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S1950;
            float depth_accum_3 = depth_accum_2 + _S1952 * _S1950;
            _S1949 = int(1);
            _S1950 = accum_7;
            _S1951 = depth_accum_3;
        }
        else
        {
            _S1949 = int(0);
            _S1950 = 0.0f;
            _S1951 = 0.0f;
        }
        if(_S1949 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            i_8 = i_8 + int(1);
            accum_6 = _S1950;
            depth_accum_2 = _S1951;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S1954 = (F32_max((depth_accum_2 / accum_6), (0.0f)));
    *dpout_rgb_0 = dprgb_0;
    *dpdepth_0 = _S1954;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_1, float dpdepth_1, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_3)
{
    FixedArray<float, 8>  _S1955 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1956 = *dpray_o_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1957 = *dpray_d_3;
    float3  _S1958 = make_float3 (size_12);
    float _S1959 = 0.5f * size_12;
    float3  _S1960 = make_float3 (_S1959);
    float3  m_7 = make_float3 (1.0f) / (*dpray_d_3).primal_0;
    float3  _S1961 = (*dpray_d_3).primal_0 * (*dpray_d_3).primal_0;
    float3  _S1962 = (*dpray_o_3).primal_0 - (pos_12 + make_float3 (_S1959));
    float3  k_9 = s_primal_ctx_abs_0(m_7) * make_float3 (_S1959);
    float3  _S1963 = - (m_7 * _S1962);
    float3  ta_6 = _S1963 - k_9;
    float3  tb_6 = _S1963 + k_9;
    float _S1964 = ta_6.x;
    float _S1965 = ta_6.y;
    float _S1966 = (F32_max((_S1964), (_S1965)));
    float _S1967 = ta_6.z;
    float _S1968 = (F32_max((_S1967), (0.0f)));
    float _S1969 = (F32_max((_S1966), (_S1968)));
    float _S1970 = tb_6.x;
    float _S1971 = tb_6.y;
    float _S1972 = (F32_min((_S1970), (_S1971)));
    float _S1973 = tb_6.z;
    float _S1974 = (F32_min((_S1972), (_S1973)));
    float _S1975 = _s_diff_ctx_3->_S1944 * _s_diff_ctx_3->_S1944;
    float3  _S1976 = make_float3 (0.0f);
    FixedArray<float, 8>  _S1977 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S1978;
    (&_S1978)->primal_0 = _s_diff_ctx_3->_S1943 / _s_diff_ctx_3->_S1944;
    (&_S1978)->differential_0 = 0.0f;
    DiffPair_float_0 _S1979;
    (&_S1979)->primal_0 = 0.0f;
    (&_S1979)->differential_0 = 0.0f;
    _d_max_0(&_S1978, &_S1979, dpdepth_1);
    float _S1980 = _S1978.differential_0 / _S1975;
    float _S1981 = _s_diff_ctx_3->_S1943 * - _S1980;
    float _S1982 = _s_diff_ctx_3->_S1944 * _S1980;
    int _dc_1 = int(8);
    float _S1983 = _S1981;
    float _S1984 = _S1982;
    FixedArray<float, 8>  _S1985;
    _S1985[int(0)] = 0.0f;
    _S1985[int(1)] = 0.0f;
    _S1985[int(2)] = 0.0f;
    _S1985[int(3)] = 0.0f;
    _S1985[int(4)] = 0.0f;
    _S1985[int(5)] = 0.0f;
    _S1985[int(6)] = 0.0f;
    _S1985[int(7)] = 0.0f;
    float3  _S1986 = _S1976;
    float3  _S1987 = _S1976;
    float _S1988 = 0.0f;
    float _S1989 = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S1990 = _dc_1 < int(8);
        int _S1991;
        float _S1992;
        float _S1993;
        float _S1994;
        float _S1995;
        float3  _S1996;
        float3  _S1997;
        bool _S1998;
        if(_S1990)
        {
            float _S1999 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S2000 = s_primal_ctx_lerp_0(_S1969, _S1974, _S1999);
            float3  _S2001 = make_float3 (_S2000);
            float3  _S2002 = (_S1956.primal_0 + _S1957.primal_0 * make_float3 (_S2000) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S2003 = _S1955;
            float _S2004 = s_primal_ctx_interp_0(&_S2003, _S2002);
            bool _S2005 = _S2004 > 1.10000002384185791f;
            if(_S2005)
            {
                _S1992 = _S2004;
                _S1993 = 0.0f;
            }
            else
            {
                float _S2006 = 0.90909093618392944f * _S2004 - 0.90468984842300415f;
                _S1992 = s_primal_ctx_exp_0(_S2006);
                _S1993 = _S2006;
            }
            float _S2007 = _S1992;
            float _S2008 = _S1993;
            _S1991 = int(1);
            _S1992 = _S2000;
            _S1993 = _S2007;
            _S1998 = _S2005;
            _S1994 = _S2008;
            _S1996 = _S2002;
            _S1997 = _S2001;
            _S1995 = _S1999;
        }
        else
        {
            _S1991 = int(0);
            _S1992 = 0.0f;
            _S1993 = 0.0f;
            _S1998 = false;
            _S1994 = 0.0f;
            _S1996 = _S1976;
            _S1997 = _S1976;
            _S1995 = 0.0f;
        }
        float _S2009;
        float _S2010;
        float _S2011;
        float _S2012;
        if(!(_S1991 != int(1)))
        {
            _S2009 = _S1983;
            _S2010 = _S1984;
            _S2011 = 0.0f;
            _S2012 = 0.0f;
        }
        else
        {
            _S2009 = 0.0f;
            _S2010 = 0.0f;
            _S2011 = _S1984;
            _S2012 = _S1983;
        }
        if(_S1990)
        {
            float _S2013 = _S1993 * _S2010;
            float _S2014 = _S2010 + _S2011;
            float _S2015 = _S1992 * _S2010 + _S2009;
            float _S2016 = _S2009 + _S2012;
            float _S2017;
            if(_S1998)
            {
                _S2017 = _S2015;
            }
            else
            {
                DiffPair_float_0 _S2018;
                (&_S2018)->primal_0 = _S1994;
                (&_S2018)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S2018, _S2015);
                _S2017 = 0.90909093618392944f * _S2018.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S2019;
            (&_S2019)->primal_0 = _S1955;
            (&_S2019)->differential_0 = _S1977;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2020;
            (&_S2020)->primal_0 = _S1996;
            (&_S2020)->differential_0 = _S1976;
            s_bwd_prop_interp_0(&_S2019, &_S2020, _S2017);
            float3  _S2021 = _S2020.differential_0 / _S1958;
            float3  _S2022 = _S1957.primal_0 * _S2021;
            float3  _S2023 = _S1997 * _S2021;
            float _S2024 = _S2022.x + _S2022.y + _S2022.z + _S2013;
            DiffPair_float_0 _S2025;
            (&_S2025)->primal_0 = _S1969;
            (&_S2025)->differential_0 = 0.0f;
            DiffPair_float_0 _S2026;
            (&_S2026)->primal_0 = _S1974;
            (&_S2026)->differential_0 = 0.0f;
            DiffPair_float_0 _S2027;
            (&_S2027)->primal_0 = _S1995;
            (&_S2027)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S2025, &_S2026, &_S2027, _S2024);
            float _S2028 = (&_S2019)->differential_0[int(0)] + _S1985[int(0)];
            float _S2029 = (&_S2019)->differential_0[int(1)] + _S1985[int(1)];
            float _S2030 = (&_S2019)->differential_0[int(2)] + _S1985[int(2)];
            float _S2031 = (&_S2019)->differential_0[int(3)] + _S1985[int(3)];
            float _S2032 = (&_S2019)->differential_0[int(4)] + _S1985[int(4)];
            float _S2033 = (&_S2019)->differential_0[int(5)] + _S1985[int(5)];
            float _S2034 = (&_S2019)->differential_0[int(6)] + _S1985[int(6)];
            float _S2035 = (&_S2019)->differential_0[int(7)] + _S1985[int(7)];
            float3  _S2036 = _S2021 + _S1986;
            float3  _S2037 = _S2023 + _S1987;
            float _S2038 = _S2026.differential_0 + _S1988;
            float _S2039 = _S2025.differential_0 + _S1989;
            _S1983 = _S2016;
            _S1984 = _S2014;
            _S1985[int(0)] = _S2028;
            _S1985[int(1)] = _S2029;
            _S1985[int(2)] = _S2030;
            _S1985[int(3)] = _S2031;
            _S1985[int(4)] = _S2032;
            _S1985[int(5)] = _S2033;
            _S1985[int(6)] = _S2034;
            _S1985[int(7)] = _S2035;
            _S1986 = _S2036;
            _S1987 = _S2037;
            _S1988 = _S2038;
            _S1989 = _S2039;
        }
        else
        {
            _S1983 = _S2012;
            _S1984 = _S2011;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S2040;
    (&_S2040)->primal_0 = _S1972;
    (&_S2040)->differential_0 = 0.0f;
    DiffPair_float_0 _S2041;
    (&_S2041)->primal_0 = _S1973;
    (&_S2041)->differential_0 = 0.0f;
    _d_min_0(&_S2040, &_S2041, _S1988);
    DiffPair_float_0 _S2042;
    (&_S2042)->primal_0 = _S1970;
    (&_S2042)->differential_0 = 0.0f;
    DiffPair_float_0 _S2043;
    (&_S2043)->primal_0 = _S1971;
    (&_S2043)->differential_0 = 0.0f;
    _d_min_0(&_S2042, &_S2043, _S2040.differential_0);
    DiffPair_float_0 _S2044;
    (&_S2044)->primal_0 = _S1966;
    (&_S2044)->differential_0 = 0.0f;
    DiffPair_float_0 _S2045;
    (&_S2045)->primal_0 = _S1968;
    (&_S2045)->differential_0 = 0.0f;
    _d_max_0(&_S2044, &_S2045, _S1989);
    DiffPair_float_0 _S2046;
    (&_S2046)->primal_0 = _S1967;
    (&_S2046)->differential_0 = 0.0f;
    DiffPair_float_0 _S2047;
    (&_S2047)->primal_0 = 0.0f;
    (&_S2047)->differential_0 = 0.0f;
    _d_max_0(&_S2046, &_S2047, _S2045.differential_0);
    DiffPair_float_0 _S2048;
    (&_S2048)->primal_0 = _S1964;
    (&_S2048)->differential_0 = 0.0f;
    DiffPair_float_0 _S2049;
    (&_S2049)->primal_0 = _S1965;
    (&_S2049)->differential_0 = 0.0f;
    _d_max_0(&_S2048, &_S2049, _S2044.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S2042.differential_0, _S2043.differential_0, _S2041.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S2048.differential_0, _S2049.differential_0, _S2046.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S2050 = _S1960 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2051;
    (&_S2051)->primal_0 = m_7;
    (&_S2051)->differential_0 = _S1976;
    s_bwd_prop_abs_0(&_S2051, _S2050);
    float3  _S2052 = m_7 * s_diff_n_T_1;
    float3  _S2053 = - ((_S2051.differential_0 + _S1962 * s_diff_n_T_1) / _S1961) + _S1987;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2053;
    float3  _S2054 = _S2052 + _S1986;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2054;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = dpout_rgb_1;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S1985;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S2055, float _S2056, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S2057, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2058, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2059, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2060, float3  _S2061, float _S2062)
{
    FixedArray<float, 8>  _S2063 = _S2057->primal_0;
    float3  _S2064;
    float _S2065;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2066;
    s_primal_ctx_evaluate_color_voxel_0(_S2055, _S2056, &_S2063, (*_S2058).primal_0, (*_S2059).primal_0, (*_S2060).primal_0, &_S2064, &_S2065, &_S2066);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S2067 = _S2066;
    s_bwd_prop_evaluate_color_voxel_0(_S2055, _S2056, _S2057, _S2058, _S2059, _S2060, _S2061, _S2062, &_S2067);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  densities_10, float3  rgb_1, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_0, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_2, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    FixedArray<float, 8>  _S2068 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = densities_10;
    (&dp_densities_1)->differential_0 = _S2068;
    float3  _S2069 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S2069;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2069;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2069;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_2 = dp_rgb_0.differential_0;
    *v_ray_o_1 = dp_ray_o_1.differential_0;
    *v_ray_d_1 = dp_ray_d_1.differential_0;
    return;
}

