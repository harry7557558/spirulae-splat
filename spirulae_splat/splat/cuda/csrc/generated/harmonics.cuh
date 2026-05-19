#pragma once

#include "slang.cuh"

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_0)
{
    Matrix<float, 3, 3>  result_0;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_0)->rows + (r_0)), c_0) = _slang_vector_get_element(x_0.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_0;
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_0)
{
    float _S1 = (*left_0).primal_0.rows[int(0)].x * dOut_0.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_0.x;
    float sum_0 = _S1 + (*left_0).primal_0.rows[int(1)].x * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_0.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_0.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S2 = (*left_0).primal_0.rows[int(0)].y * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_0.x;
    float sum_2 = _S2 + (*left_0).primal_0.rows[int(1)].y * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_0.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_0.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S3 = (*left_0).primal_0.rows[int(0)].z * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_0.x;
    float sum_4 = _S3 + (*left_0).primal_0.rows[int(1)].z * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_0.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_0.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_1;
    int i_0 = int(0);
    for(;;)
    {
        if(i_0 < int(3))
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
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_0], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_1, i_0) = sum_6;
        i_0 = i_0 + int(1);
    }
    return result_1;
}

inline __device__ float3  max_0(float3  x_1, float3  y_0)
{
    float3  result_2;
    int i_1 = int(0);
    for(;;)
    {
        if(i_1 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_2, i_1) = (F32_max((_slang_vector_get_element(x_1, i_1)), (_slang_vector_get_element(y_0, i_1))));
        i_1 = i_1 + int(1);
    }
    return result_2;
}

inline __device__ float3  sh0_to_color(float3  mean_0, Matrix<float, 3, 3>  R_0, float3  t_0, float3  coeff_dc_0, float3  * coeffs_0)
{
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_0 + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ float dot_0(float3  x_2, float3  y_1)
{
    int i_2 = int(0);
    float result_3 = 0.0f;
    for(;;)
    {
        if(i_2 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_4 = result_3 + _slang_vector_get_element(x_2, i_2) * _slang_vector_get_element(y_1, i_2);
        i_2 = i_2 + int(1);
        result_3 = result_4;
    }
    return result_3;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S4, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5, float3  _S6)
{
    _d_mul_0(_S4, _S5, _S6);
    return;
}

inline __device__ void sh0_to_color_vjp_inplace(float3  mean_1, Matrix<float, 3, 3>  R_1, float3  t_1, float3  coeff_dc_1, float3  * coeffs_1, float3  v_colors_0, float3  * v_coeff_dc_0, float3  * v_coeffs_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  colors_0 = make_float3 (0.282094806432724f) * coeff_dc_1;
    *v_coeff_dc_0 = *v_coeff_dc_0 + make_float3 (0.282094806432724f) * (v_colors_0 * make_float3 (float((colors_0.x) >= -0.5f), float((colors_0.y) >= -0.5f), float((colors_0.z) >= -0.5f)));
    float3  v_viewdir_0 = {};
    Matrix<float, 3, 3>  _S7 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S8;
    (&_S8)->primal_0 = transpose_0(R_1);
    (&_S8)->differential_0 = _S7;
    float3  _S9 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S10;
    (&_S10)->primal_0 = t_1;
    (&_S10)->differential_0 = _S9;
    s_bwd_prop_mul_0(&_S8, &_S10, v_viewdir_0);
    Matrix<float, 3, 3>  _S11 = transpose_0(_S8.differential_0);
    *v_mean_0 = *v_mean_0 + v_viewdir_0;
    *v_R_0 = *v_R_0 + _S11;
    *v_t_0 = *v_t_0 + _S10.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_atomic(float3  mean_2, Matrix<float, 3, 3>  R_2, float3  t_2, float3  coeff_dc_2, float3  * coeffs_2, float3  v_colors_1, float3  * v_coeff_dc_1, float3  * v_coeffs_1, float3  * v_mean_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  colors_1 = make_float3 (0.282094806432724f) * coeff_dc_2;
    *v_coeff_dc_1 = *v_coeff_dc_1 + make_float3 (0.282094806432724f) * (v_colors_1 * make_float3 (float((colors_1.x) >= -0.5f), float((colors_1.y) >= -0.5f), float((colors_1.z) >= -0.5f)));
    float3  v_viewdir_1 = {};
    Matrix<float, 3, 3>  _S12 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S13;
    (&_S13)->primal_0 = transpose_0(R_2);
    (&_S13)->differential_0 = _S12;
    float3  _S14 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S15;
    (&_S15)->primal_0 = t_2;
    (&_S15)->differential_0 = _S14;
    s_bwd_prop_mul_0(&_S13, &_S15, v_viewdir_1);
    Matrix<float, 3, 3>  _S16 = transpose_0(_S13.differential_0);
    *v_mean_1 = *v_mean_1 + v_viewdir_1;
    *v_R_1 = *v_R_1 + _S16;
    *v_t_1 = *v_t_1 + _S15.differential_0;
    return;
}

inline __device__ float3  sh1_to_color(float3  mean_3, Matrix<float, 3, 3>  R_3, float3  t_3, float3  coeff_dc_3, float3  * coeffs_3)
{
    float3  _S17 = mean_3 + mul_0(transpose_0(R_3), t_3);
    float _S18 = _S17.x;
    float _S19 = _S17.y;
    float _S20 = _S17.z;
    float inv_norm_0 = (F32_rsqrt((_S18 * _S18 + _S19 * _S19 + _S20 * _S20)));
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_3 + make_float3 (0.48860251903533936f) * (make_float3 (- (_S19 * inv_norm_0)) * *(coeffs_3 + int(0)) + make_float3 (_S20 * inv_norm_0) * *(coeffs_3 + int(1)) - make_float3 (_S18 * inv_norm_0) * *(coeffs_3 + int(2))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh1_to_color_vjp_inplace(float3  mean_4, Matrix<float, 3, 3>  R_4, float3  t_4, float3  coeff_dc_4, float3  * coeffs_4, float3  v_colors_2, float3  * v_coeff_dc_2, float3  * v_coeffs_2, float3  * v_mean_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 3, 3>  _S21 = transpose_0(R_4);
    float3  _S22 = mean_4 + mul_0(_S21, t_4);
    float _S23 = _S22.x;
    float _S24 = _S22.y;
    float _S25 = _S22.z;
    float inv_norm_1 = (F32_rsqrt((_S23 * _S23 + _S24 * _S24 + _S25 * _S25)));
    float x_3 = _S23 * inv_norm_1;
    float y_2 = _S24 * inv_norm_1;
    float z_0 = _S25 * inv_norm_1;
    float3  * _S26 = coeffs_4 + int(0);
    float3  * _S27 = coeffs_4 + int(1);
    float3  * _S28 = coeffs_4 + int(2);
    float3  colors_2 = make_float3 (0.282094806432724f) * coeff_dc_4 + make_float3 (0.48860251903533936f) * (make_float3 (- y_2) * *_S26 + make_float3 (z_0) * *_S27 - make_float3 (x_3) * *_S28);
    float3  _S29 = v_colors_2 * make_float3 (float((colors_2.x) >= -0.5f), float((colors_2.y) >= -0.5f), float((colors_2.z) >= -0.5f));
    float3  v_viewdir_2 = {};
    *v_coeff_dc_2 = *v_coeff_dc_2 + make_float3 (0.282094806432724f) * _S29;
    float3  * _S30 = v_coeffs_2 + int(0);
    *_S30 = *_S30 + make_float3 (-0.48860251903533936f * y_2) * _S29;
    float3  * _S31 = v_coeffs_2 + int(1);
    *_S31 = *_S31 + make_float3 (0.48860251903533936f * z_0) * _S29;
    float3  * _S32 = v_coeffs_2 + int(2);
    *_S32 = *_S32 + make_float3 (-0.48860251903533936f * x_3) * _S29;
    float3  dir_n_0 = make_float3 (x_3, y_2, z_0);
    float3  v_dir_n_0 = make_float3 (-0.48860251903533936f * dot_0(*_S28, _S29), -0.48860251903533936f * dot_0(*_S26, _S29), 0.48860251903533936f * dot_0(*_S27, _S29));
    float3  v_viewdir_3 = v_viewdir_2 + (v_dir_n_0 - make_float3 (dot_0(v_dir_n_0, dir_n_0)) * dir_n_0) * make_float3 (inv_norm_1);
    Matrix<float, 3, 3>  _S33 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S34;
    (&_S34)->primal_0 = _S21;
    (&_S34)->differential_0 = _S33;
    float3  _S35 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S36;
    (&_S36)->primal_0 = t_4;
    (&_S36)->differential_0 = _S35;
    s_bwd_prop_mul_0(&_S34, &_S36, v_viewdir_3);
    Matrix<float, 3, 3>  _S37 = transpose_0(_S34.differential_0);
    *v_mean_2 = *v_mean_2 + v_viewdir_3;
    *v_R_2 = *v_R_2 + _S37;
    *v_t_2 = *v_t_2 + _S36.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_atomic(float3  mean_5, Matrix<float, 3, 3>  R_5, float3  t_5, float3  coeff_dc_5, float3  * coeffs_5, float3  v_colors_3, float3  * v_coeff_dc_3, float3  * v_coeffs_3, float3  * v_mean_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    Matrix<float, 3, 3>  _S38 = transpose_0(R_5);
    float3  _S39 = mean_5 + mul_0(_S38, t_5);
    float _S40 = _S39.x;
    float _S41 = _S39.y;
    float _S42 = _S39.z;
    float inv_norm_2 = (F32_rsqrt((_S40 * _S40 + _S41 * _S41 + _S42 * _S42)));
    float x_4 = _S40 * inv_norm_2;
    float y_3 = _S41 * inv_norm_2;
    float z_1 = _S42 * inv_norm_2;
    float3  * _S43 = coeffs_5 + int(0);
    float3  * _S44 = coeffs_5 + int(1);
    float3  * _S45 = coeffs_5 + int(2);
    float3  colors_3 = make_float3 (0.282094806432724f) * coeff_dc_5 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * *_S43 + make_float3 (z_1) * *_S44 - make_float3 (x_4) * *_S45);
    float3  _S46 = v_colors_3 * make_float3 (float((colors_3.x) >= -0.5f), float((colors_3.y) >= -0.5f), float((colors_3.z) >= -0.5f));
    float3  v_viewdir_4 = {};
    *v_coeff_dc_3 = *v_coeff_dc_3 + make_float3 (0.282094806432724f) * _S46;
    float3  temp_0 = make_float3 (-0.48860251903533936f * y_3) * _S46;
    float _S47 = dot_0(temp_0, temp_0);
    bool _S48;
    if((F32_isfinite((_S47))))
    {
        _S48 = _S47 != 0.0f;
    }
    else
    {
        _S48 = false;
    }
    if(_S48)
    {
        float3  * _S49 = v_coeffs_3 + int(0);
        float _S50 = atomicAdd(&(_S49->x), temp_0.x);
        float _S51 = atomicAdd(&(_S49->y), temp_0.y);
        float _S52 = atomicAdd(&(_S49->z), temp_0.z);
    }
    float3  temp_1 = make_float3 (0.48860251903533936f * z_1) * _S46;
    float _S53 = dot_0(temp_1, temp_1);
    if((F32_isfinite((_S53))))
    {
        _S48 = _S53 != 0.0f;
    }
    else
    {
        _S48 = false;
    }
    if(_S48)
    {
        float3  * _S54 = v_coeffs_3 + int(1);
        float _S55 = atomicAdd(&(_S54->x), temp_1.x);
        float _S56 = atomicAdd(&(_S54->y), temp_1.y);
        float _S57 = atomicAdd(&(_S54->z), temp_1.z);
    }
    float3  temp_2 = make_float3 (-0.48860251903533936f * x_4) * _S46;
    float _S58 = dot_0(temp_2, temp_2);
    if((F32_isfinite((_S58))))
    {
        _S48 = _S58 != 0.0f;
    }
    else
    {
        _S48 = false;
    }
    if(_S48)
    {
        float3  * _S59 = v_coeffs_3 + int(2);
        float _S60 = atomicAdd(&(_S59->x), temp_2.x);
        float _S61 = atomicAdd(&(_S59->y), temp_2.y);
        float _S62 = atomicAdd(&(_S59->z), temp_2.z);
    }
    float3  dir_n_1 = make_float3 (x_4, y_3, z_1);
    float3  v_dir_n_1 = make_float3 (-0.48860251903533936f * dot_0(*_S45, _S46), -0.48860251903533936f * dot_0(*_S43, _S46), 0.48860251903533936f * dot_0(*_S44, _S46));
    float3  v_viewdir_5 = v_viewdir_4 + (v_dir_n_1 - make_float3 (dot_0(v_dir_n_1, dir_n_1)) * dir_n_1) * make_float3 (inv_norm_2);
    Matrix<float, 3, 3>  _S63 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S64;
    (&_S64)->primal_0 = _S38;
    (&_S64)->differential_0 = _S63;
    float3  _S65 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S66;
    (&_S66)->primal_0 = t_5;
    (&_S66)->differential_0 = _S65;
    s_bwd_prop_mul_0(&_S64, &_S66, v_viewdir_5);
    Matrix<float, 3, 3>  _S67 = transpose_0(_S64.differential_0);
    *v_mean_3 = *v_mean_3 + v_viewdir_5;
    *v_R_3 = *v_R_3 + _S67;
    *v_t_3 = *v_t_3 + _S66.differential_0;
    return;
}

inline __device__ float3  sh2_to_color(float3  mean_6, Matrix<float, 3, 3>  R_6, float3  t_6, float3  coeff_dc_6, float3  * coeffs_6)
{
    float3  _S68 = mean_6 + mul_0(transpose_0(R_6), t_6);
    float _S69 = _S68.x;
    float _S70 = _S68.y;
    float _S71 = _S68.z;
    float inv_norm_3 = (F32_rsqrt((_S69 * _S69 + _S70 * _S70 + _S71 * _S71)));
    float x_5 = _S69 * inv_norm_3;
    float y_4 = _S70 * inv_norm_3;
    float z_2 = _S71 * inv_norm_3;
    float fTmp0B_0 = -1.09254848957061768f * z_2;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_6 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * *(coeffs_6 + int(0)) + make_float3 (z_2) * *(coeffs_6 + int(1)) - make_float3 (x_5) * *(coeffs_6 + int(2))) + (make_float3 (0.54627424478530884f * (2.0f * x_5 * y_4)) * *(coeffs_6 + int(3)) + make_float3 (fTmp0B_0 * y_4) * *(coeffs_6 + int(4)) + make_float3 (0.94617468118667603f * (z_2 * z_2) - 0.31539157032966614f) * *(coeffs_6 + int(5)) + make_float3 (fTmp0B_0 * x_5) * *(coeffs_6 + int(6)) + make_float3 (0.54627424478530884f * (x_5 * x_5 - y_4 * y_4)) * *(coeffs_6 + int(7))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh2_to_color_vjp_inplace(float3  mean_7, Matrix<float, 3, 3>  R_7, float3  t_7, float3  coeff_dc_7, float3  * coeffs_7, float3  v_colors_4, float3  * v_coeff_dc_4, float3  * v_coeffs_4, float3  * v_mean_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 3, 3>  _S72 = transpose_0(R_7);
    float3  _S73 = mean_7 + mul_0(_S72, t_7);
    float _S74 = _S73.x;
    float _S75 = _S73.y;
    float _S76 = _S73.z;
    float inv_norm_4 = (F32_rsqrt((_S74 * _S74 + _S75 * _S75 + _S76 * _S76)));
    float x_6 = _S74 * inv_norm_4;
    float y_5 = _S75 * inv_norm_4;
    float z_3 = _S76 * inv_norm_4;
    float3  * _S77 = coeffs_7 + int(0);
    float3  * _S78 = coeffs_7 + int(1);
    float3  * _S79 = coeffs_7 + int(2);
    float fTmp0B_1 = -1.09254848957061768f * z_3;
    float _S80 = 2.0f * x_6;
    float pSH6_0 = 0.94617468118667603f * (z_3 * z_3) - 0.31539157032966614f;
    float pSH7_0 = fTmp0B_1 * x_6;
    float pSH5_0 = fTmp0B_1 * y_5;
    float pSH8_0 = 0.54627424478530884f * (x_6 * x_6 - y_5 * y_5);
    float pSH4_0 = 0.54627424478530884f * (_S80 * y_5);
    float3  * _S81 = coeffs_7 + int(3);
    float3  * _S82 = coeffs_7 + int(4);
    float3  * _S83 = coeffs_7 + int(5);
    float3  * _S84 = coeffs_7 + int(6);
    float3  * _S85 = coeffs_7 + int(7);
    float3  colors_4 = make_float3 (0.282094806432724f) * coeff_dc_7 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * *_S77 + make_float3 (z_3) * *_S78 - make_float3 (x_6) * *_S79) + (make_float3 (pSH4_0) * *_S81 + make_float3 (pSH5_0) * *_S82 + make_float3 (pSH6_0) * *_S83 + make_float3 (pSH7_0) * *_S84 + make_float3 (pSH8_0) * *_S85);
    float3  _S86 = v_colors_4 * make_float3 (float((colors_4.x) >= -0.5f), float((colors_4.y) >= -0.5f), float((colors_4.z) >= -0.5f));
    *v_coeff_dc_4 = *v_coeff_dc_4 + make_float3 (0.282094806432724f) * _S86;
    float3  v_viewdir_6 = {};
    float3  * _S87 = v_coeffs_4 + int(0);
    *_S87 = *_S87 + make_float3 (-0.48860251903533936f * y_5) * _S86;
    float3  * _S88 = v_coeffs_4 + int(1);
    *_S88 = *_S88 + make_float3 (0.48860251903533936f * z_3) * _S86;
    float3  * _S89 = v_coeffs_4 + int(2);
    *_S89 = *_S89 + make_float3 (-0.48860251903533936f * x_6) * _S86;
    float _S90 = -0.48860251903533936f * dot_0(*_S79, _S86);
    float _S91 = -0.48860251903533936f * dot_0(*_S77, _S86);
    float _S92 = 0.48860251903533936f * dot_0(*_S78, _S86);
    float3  * _S93 = v_coeffs_4 + int(3);
    *_S93 = *_S93 + make_float3 (pSH4_0) * _S86;
    float3  * _S94 = v_coeffs_4 + int(4);
    *_S94 = *_S94 + make_float3 (pSH5_0) * _S86;
    float3  * _S95 = v_coeffs_4 + int(5);
    *_S95 = *_S95 + make_float3 (pSH6_0) * _S86;
    float3  * _S96 = v_coeffs_4 + int(6);
    *_S96 = *_S96 + make_float3 (pSH7_0) * _S86;
    float3  * _S97 = v_coeffs_4 + int(7);
    *_S97 = *_S97 + make_float3 (pSH8_0) * _S86;
    float pSH8_x_0 = 0.54627424478530884f * _S80;
    float3  dir_n_2 = make_float3 (x_6, y_5, z_3);
    float3  v_dir_n_2 = make_float3 (_S90 + dot_0(_S86, make_float3 (0.54627424478530884f * (2.0f * y_5)) * *_S81 + make_float3 (pSH8_x_0) * *_S85 + make_float3 (fTmp0B_1) * *_S84), _S91 + dot_0(_S86, make_float3 (pSH8_x_0) * *_S81 + make_float3 (0.54627424478530884f * (-2.0f * y_5)) * *_S85 + make_float3 (fTmp0B_1) * *_S82), _S92 + dot_0(_S86, make_float3 (1.89234936237335205f * z_3) * *_S83 + make_float3 (-1.09254848957061768f * x_6) * *_S84 + make_float3 (-1.09254848957061768f * y_5) * *_S82));
    float3  v_viewdir_7 = v_viewdir_6 + (v_dir_n_2 - make_float3 (dot_0(v_dir_n_2, dir_n_2)) * dir_n_2) * make_float3 (inv_norm_4);
    Matrix<float, 3, 3>  _S98 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S99;
    (&_S99)->primal_0 = _S72;
    (&_S99)->differential_0 = _S98;
    float3  _S100 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S101;
    (&_S101)->primal_0 = t_7;
    (&_S101)->differential_0 = _S100;
    s_bwd_prop_mul_0(&_S99, &_S101, v_viewdir_7);
    Matrix<float, 3, 3>  _S102 = transpose_0(_S99.differential_0);
    *v_mean_4 = *v_mean_4 + v_viewdir_7;
    *v_R_4 = *v_R_4 + _S102;
    *v_t_4 = *v_t_4 + _S101.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_atomic(float3  mean_8, Matrix<float, 3, 3>  R_8, float3  t_8, float3  coeff_dc_8, float3  * coeffs_8, float3  v_colors_5, float3  * v_coeff_dc_5, float3  * v_coeffs_5, float3  * v_mean_5, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 3, 3>  _S103 = transpose_0(R_8);
    float3  _S104 = mean_8 + mul_0(_S103, t_8);
    float _S105 = _S104.x;
    float _S106 = _S104.y;
    float _S107 = _S104.z;
    float inv_norm_5 = (F32_rsqrt((_S105 * _S105 + _S106 * _S106 + _S107 * _S107)));
    float x_7 = _S105 * inv_norm_5;
    float y_6 = _S106 * inv_norm_5;
    float z_4 = _S107 * inv_norm_5;
    float3  * _S108 = coeffs_8 + int(0);
    float3  * _S109 = coeffs_8 + int(1);
    float3  * _S110 = coeffs_8 + int(2);
    float fTmp0B_2 = -1.09254848957061768f * z_4;
    float _S111 = 2.0f * x_7;
    float pSH6_1 = 0.94617468118667603f * (z_4 * z_4) - 0.31539157032966614f;
    float pSH7_1 = fTmp0B_2 * x_7;
    float pSH5_1 = fTmp0B_2 * y_6;
    float pSH8_1 = 0.54627424478530884f * (x_7 * x_7 - y_6 * y_6);
    float pSH4_1 = 0.54627424478530884f * (_S111 * y_6);
    float3  * _S112 = coeffs_8 + int(3);
    float3  * _S113 = coeffs_8 + int(4);
    float3  * _S114 = coeffs_8 + int(5);
    float3  * _S115 = coeffs_8 + int(6);
    float3  * _S116 = coeffs_8 + int(7);
    float3  colors_5 = make_float3 (0.282094806432724f) * coeff_dc_8 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * *_S108 + make_float3 (z_4) * *_S109 - make_float3 (x_7) * *_S110) + (make_float3 (pSH4_1) * *_S112 + make_float3 (pSH5_1) * *_S113 + make_float3 (pSH6_1) * *_S114 + make_float3 (pSH7_1) * *_S115 + make_float3 (pSH8_1) * *_S116);
    float3  _S117 = v_colors_5 * make_float3 (float((colors_5.x) >= -0.5f), float((colors_5.y) >= -0.5f), float((colors_5.z) >= -0.5f));
    *v_coeff_dc_5 = *v_coeff_dc_5 + make_float3 (0.282094806432724f) * _S117;
    float3  v_viewdir_8 = {};
    float3  temp_3 = make_float3 (-0.48860251903533936f * y_6) * _S117;
    float _S118 = dot_0(temp_3, temp_3);
    bool _S119;
    if((F32_isfinite((_S118))))
    {
        _S119 = _S118 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S120 = v_coeffs_5 + int(0);
        float _S121 = atomicAdd(&(_S120->x), temp_3.x);
        float _S122 = atomicAdd(&(_S120->y), temp_3.y);
        float _S123 = atomicAdd(&(_S120->z), temp_3.z);
    }
    float3  temp_4 = make_float3 (0.48860251903533936f * z_4) * _S117;
    float _S124 = dot_0(temp_4, temp_4);
    if((F32_isfinite((_S124))))
    {
        _S119 = _S124 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S125 = v_coeffs_5 + int(1);
        float _S126 = atomicAdd(&(_S125->x), temp_4.x);
        float _S127 = atomicAdd(&(_S125->y), temp_4.y);
        float _S128 = atomicAdd(&(_S125->z), temp_4.z);
    }
    float3  temp_5 = make_float3 (-0.48860251903533936f * x_7) * _S117;
    float _S129 = dot_0(temp_5, temp_5);
    if((F32_isfinite((_S129))))
    {
        _S119 = _S129 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S130 = v_coeffs_5 + int(2);
        float _S131 = atomicAdd(&(_S130->x), temp_5.x);
        float _S132 = atomicAdd(&(_S130->y), temp_5.y);
        float _S133 = atomicAdd(&(_S130->z), temp_5.z);
    }
    float _S134 = -0.48860251903533936f * dot_0(*_S110, _S117);
    float _S135 = -0.48860251903533936f * dot_0(*_S108, _S117);
    float _S136 = 0.48860251903533936f * dot_0(*_S109, _S117);
    float3  temp_6 = make_float3 (pSH4_1) * _S117;
    float _S137 = dot_0(temp_6, temp_6);
    if((F32_isfinite((_S137))))
    {
        _S119 = _S137 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S138 = v_coeffs_5 + int(3);
        float _S139 = atomicAdd(&(_S138->x), temp_6.x);
        float _S140 = atomicAdd(&(_S138->y), temp_6.y);
        float _S141 = atomicAdd(&(_S138->z), temp_6.z);
    }
    float3  temp_7 = make_float3 (pSH5_1) * _S117;
    float _S142 = dot_0(temp_7, temp_7);
    if((F32_isfinite((_S142))))
    {
        _S119 = _S142 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S143 = v_coeffs_5 + int(4);
        float _S144 = atomicAdd(&(_S143->x), temp_7.x);
        float _S145 = atomicAdd(&(_S143->y), temp_7.y);
        float _S146 = atomicAdd(&(_S143->z), temp_7.z);
    }
    float3  temp_8 = make_float3 (pSH6_1) * _S117;
    float _S147 = dot_0(temp_8, temp_8);
    if((F32_isfinite((_S147))))
    {
        _S119 = _S147 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S148 = v_coeffs_5 + int(5);
        float _S149 = atomicAdd(&(_S148->x), temp_8.x);
        float _S150 = atomicAdd(&(_S148->y), temp_8.y);
        float _S151 = atomicAdd(&(_S148->z), temp_8.z);
    }
    float3  temp_9 = make_float3 (pSH7_1) * _S117;
    float _S152 = dot_0(temp_9, temp_9);
    if((F32_isfinite((_S152))))
    {
        _S119 = _S152 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S153 = v_coeffs_5 + int(6);
        float _S154 = atomicAdd(&(_S153->x), temp_9.x);
        float _S155 = atomicAdd(&(_S153->y), temp_9.y);
        float _S156 = atomicAdd(&(_S153->z), temp_9.z);
    }
    float3  temp_10 = make_float3 (pSH8_1) * _S117;
    float _S157 = dot_0(temp_10, temp_10);
    if((F32_isfinite((_S157))))
    {
        _S119 = _S157 != 0.0f;
    }
    else
    {
        _S119 = false;
    }
    if(_S119)
    {
        float3  * _S158 = v_coeffs_5 + int(7);
        float _S159 = atomicAdd(&(_S158->x), temp_10.x);
        float _S160 = atomicAdd(&(_S158->y), temp_10.y);
        float _S161 = atomicAdd(&(_S158->z), temp_10.z);
    }
    float pSH8_x_1 = 0.54627424478530884f * _S111;
    float3  dir_n_3 = make_float3 (x_7, y_6, z_4);
    float3  v_dir_n_3 = make_float3 (_S134 + dot_0(_S117, make_float3 (0.54627424478530884f * (2.0f * y_6)) * *_S112 + make_float3 (pSH8_x_1) * *_S116 + make_float3 (fTmp0B_2) * *_S115), _S135 + dot_0(_S117, make_float3 (pSH8_x_1) * *_S112 + make_float3 (0.54627424478530884f * (-2.0f * y_6)) * *_S116 + make_float3 (fTmp0B_2) * *_S113), _S136 + dot_0(_S117, make_float3 (1.89234936237335205f * z_4) * *_S114 + make_float3 (-1.09254848957061768f * x_7) * *_S115 + make_float3 (-1.09254848957061768f * y_6) * *_S113));
    float3  v_viewdir_9 = v_viewdir_8 + (v_dir_n_3 - make_float3 (dot_0(v_dir_n_3, dir_n_3)) * dir_n_3) * make_float3 (inv_norm_5);
    Matrix<float, 3, 3>  _S162 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S163;
    (&_S163)->primal_0 = _S103;
    (&_S163)->differential_0 = _S162;
    float3  _S164 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S165;
    (&_S165)->primal_0 = t_8;
    (&_S165)->differential_0 = _S164;
    s_bwd_prop_mul_0(&_S163, &_S165, v_viewdir_9);
    Matrix<float, 3, 3>  _S166 = transpose_0(_S163.differential_0);
    *v_mean_5 = *v_mean_5 + v_viewdir_9;
    *v_R_5 = *v_R_5 + _S166;
    *v_t_5 = *v_t_5 + _S165.differential_0;
    return;
}

inline __device__ float3  sh3_to_color(float3  mean_9, Matrix<float, 3, 3>  R_9, float3  t_9, float3  coeff_dc_9, float3  * coeffs_9)
{
    float3  _S167 = mean_9 + mul_0(transpose_0(R_9), t_9);
    float _S168 = _S167.x;
    float _S169 = _S167.y;
    float _S170 = _S167.z;
    float inv_norm_6 = (F32_rsqrt((_S168 * _S168 + _S169 * _S169 + _S170 * _S170)));
    float x_8 = _S168 * inv_norm_6;
    float y_7 = _S169 * inv_norm_6;
    float z_5 = _S170 * inv_norm_6;
    float z2_0 = z_5 * z_5;
    float fTmp0B_3 = -1.09254848957061768f * z_5;
    float fC1_0 = x_8 * x_8 - y_7 * y_7;
    float fS1_0 = 2.0f * x_8 * y_7;
    float fTmp0C_0 = -2.28522896766662598f * z2_0 + 0.4570457935333252f;
    float fTmp1B_0 = 1.44530570507049561f * z_5;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_9 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * *(coeffs_9 + int(0)) + make_float3 (z_5) * *(coeffs_9 + int(1)) - make_float3 (x_8) * *(coeffs_9 + int(2))) + (make_float3 (0.54627424478530884f * fS1_0) * *(coeffs_9 + int(3)) + make_float3 (fTmp0B_3 * y_7) * *(coeffs_9 + int(4)) + make_float3 (0.94617468118667603f * z2_0 - 0.31539157032966614f) * *(coeffs_9 + int(5)) + make_float3 (fTmp0B_3 * x_8) * *(coeffs_9 + int(6)) + make_float3 (0.54627424478530884f * fC1_0) * *(coeffs_9 + int(7))) + (make_float3 (-0.59004360437393188f * (x_8 * fS1_0 + y_7 * fC1_0)) * *(coeffs_9 + int(8)) + make_float3 (fTmp1B_0 * fS1_0) * *(coeffs_9 + int(9)) + make_float3 (fTmp0C_0 * y_7) * *(coeffs_9 + int(10)) + make_float3 (z_5 * (1.86588168144226074f * z2_0 - 1.11952900886535645f)) * *(coeffs_9 + int(11)) + make_float3 (fTmp0C_0 * x_8) * *(coeffs_9 + int(12)) + make_float3 (fTmp1B_0 * fC1_0) * *(coeffs_9 + int(13)) + make_float3 (-0.59004360437393188f * (x_8 * fC1_0 - y_7 * fS1_0)) * *(coeffs_9 + int(14))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh3_to_color_vjp_inplace(float3  mean_10, Matrix<float, 3, 3>  R_10, float3  t_10, float3  coeff_dc_10, float3  * coeffs_10, float3  v_colors_6, float3  * v_coeff_dc_6, float3  * v_coeffs_6, float3  * v_mean_6, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    Matrix<float, 3, 3>  _S171 = transpose_0(R_10);
    float3  _S172 = mean_10 + mul_0(_S171, t_10);
    float _S173 = _S172.x;
    float _S174 = _S172.y;
    float _S175 = _S172.z;
    float inv_norm_7 = (F32_rsqrt((_S173 * _S173 + _S174 * _S174 + _S175 * _S175)));
    float x_9 = _S173 * inv_norm_7;
    float y_8 = _S174 * inv_norm_7;
    float z_6 = _S175 * inv_norm_7;
    float3  * _S176 = coeffs_10 + int(0);
    float3  * _S177 = coeffs_10 + int(1);
    float3  * _S178 = coeffs_10 + int(2);
    float z2_1 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_1 = x_9 * x_9 - y_8 * y_8;
    float _S179 = 2.0f * x_9;
    float fS1_1 = _S179 * y_8;
    float pSH6_2 = 0.94617468118667603f * z2_1 - 0.31539157032966614f;
    float pSH7_2 = fTmp0B_4 * x_9;
    float pSH5_2 = fTmp0B_4 * y_8;
    float pSH8_2 = 0.54627424478530884f * fC1_1;
    float pSH4_2 = 0.54627424478530884f * fS1_1;
    float3  * _S180 = coeffs_10 + int(3);
    float3  * _S181 = coeffs_10 + int(4);
    float3  * _S182 = coeffs_10 + int(5);
    float3  * _S183 = coeffs_10 + int(6);
    float3  * _S184 = coeffs_10 + int(7);
    float fTmp0C_1 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
    float fTmp1B_1 = 1.44530570507049561f * z_6;
    float pSH12_0 = z_6 * (1.86588168144226074f * z2_1 - 1.11952900886535645f);
    float pSH13_0 = fTmp0C_1 * x_9;
    float pSH11_0 = fTmp0C_1 * y_8;
    float pSH14_0 = fTmp1B_1 * fC1_1;
    float pSH10_0 = fTmp1B_1 * fS1_1;
    float pSH15_0 = -0.59004360437393188f * (x_9 * fC1_1 - y_8 * fS1_1);
    float pSH9_0 = -0.59004360437393188f * (x_9 * fS1_1 + y_8 * fC1_1);
    float3  * _S185 = coeffs_10 + int(8);
    float3  * _S186 = coeffs_10 + int(9);
    float3  * _S187 = coeffs_10 + int(10);
    float3  * _S188 = coeffs_10 + int(11);
    float3  * _S189 = coeffs_10 + int(12);
    float3  * _S190 = coeffs_10 + int(13);
    float3  * _S191 = coeffs_10 + int(14);
    float3  colors_6 = make_float3 (0.282094806432724f) * coeff_dc_10 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * *_S176 + make_float3 (z_6) * *_S177 - make_float3 (x_9) * *_S178) + (make_float3 (pSH4_2) * *_S180 + make_float3 (pSH5_2) * *_S181 + make_float3 (pSH6_2) * *_S182 + make_float3 (pSH7_2) * *_S183 + make_float3 (pSH8_2) * *_S184) + (make_float3 (pSH9_0) * *_S185 + make_float3 (pSH10_0) * *_S186 + make_float3 (pSH11_0) * *_S187 + make_float3 (pSH12_0) * *_S188 + make_float3 (pSH13_0) * *_S189 + make_float3 (pSH14_0) * *_S190 + make_float3 (pSH15_0) * *_S191);
    float3  _S192 = v_colors_6 * make_float3 (float((colors_6.x) >= -0.5f), float((colors_6.y) >= -0.5f), float((colors_6.z) >= -0.5f));
    float3  v_viewdir_10 = {};
    *v_coeff_dc_6 = *v_coeff_dc_6 + make_float3 (0.282094806432724f) * _S192;
    float3  * _S193 = v_coeffs_6 + int(0);
    *_S193 = *_S193 + make_float3 (-0.48860251903533936f * y_8) * _S192;
    float3  * _S194 = v_coeffs_6 + int(1);
    *_S194 = *_S194 + make_float3 (0.48860251903533936f * z_6) * _S192;
    float3  * _S195 = v_coeffs_6 + int(2);
    *_S195 = *_S195 + make_float3 (-0.48860251903533936f * x_9) * _S192;
    float _S196 = -0.48860251903533936f * dot_0(*_S178, _S192);
    float _S197 = -0.48860251903533936f * dot_0(*_S176, _S192);
    float _S198 = 0.48860251903533936f * dot_0(*_S177, _S192);
    float3  * _S199 = v_coeffs_6 + int(3);
    *_S199 = *_S199 + make_float3 (pSH4_2) * _S192;
    float3  * _S200 = v_coeffs_6 + int(4);
    *_S200 = *_S200 + make_float3 (pSH5_2) * _S192;
    float3  * _S201 = v_coeffs_6 + int(5);
    *_S201 = *_S201 + make_float3 (pSH6_2) * _S192;
    float3  * _S202 = v_coeffs_6 + int(6);
    *_S202 = *_S202 + make_float3 (pSH7_2) * _S192;
    float3  * _S203 = v_coeffs_6 + int(7);
    *_S203 = *_S203 + make_float3 (pSH8_2) * _S192;
    float fC1_y_0 = -2.0f * y_8;
    float fS1_x_0 = 2.0f * y_8;
    float pSH8_x_2 = 0.54627424478530884f * _S179;
    float v_x_0 = _S196 + dot_0(_S192, make_float3 (0.54627424478530884f * fS1_x_0) * *_S180 + make_float3 (pSH8_x_2) * *_S184 + make_float3 (fTmp0B_4) * *_S183);
    float v_y_0 = _S197 + dot_0(_S192, make_float3 (pSH8_x_2) * *_S180 + make_float3 (0.54627424478530884f * fC1_y_0) * *_S184 + make_float3 (fTmp0B_4) * *_S181);
    float v_z_0 = _S198 + dot_0(_S192, make_float3 (1.89234936237335205f * z_6) * *_S182 + make_float3 (-1.09254848957061768f * x_9) * *_S183 + make_float3 (-1.09254848957061768f * y_8) * *_S181);
    float3  * _S204 = v_coeffs_6 + int(8);
    *_S204 = *_S204 + make_float3 (pSH9_0) * _S192;
    float3  * _S205 = v_coeffs_6 + int(9);
    *_S205 = *_S205 + make_float3 (pSH10_0) * _S192;
    float3  * _S206 = v_coeffs_6 + int(10);
    *_S206 = *_S206 + make_float3 (pSH11_0) * _S192;
    float3  * _S207 = v_coeffs_6 + int(11);
    *_S207 = *_S207 + make_float3 (pSH12_0) * _S192;
    float3  * _S208 = v_coeffs_6 + int(12);
    *_S208 = *_S208 + make_float3 (pSH13_0) * _S192;
    float3  * _S209 = v_coeffs_6 + int(13);
    *_S209 = *_S209 + make_float3 (pSH14_0) * _S192;
    float3  * _S210 = v_coeffs_6 + int(14);
    *_S210 = *_S210 + make_float3 (pSH15_0) * _S192;
    float fTmp0C_z_0 = -4.57045793533325195f * z_6;
    float _S211 = x_9 * _S179;
    float _S212 = y_8 * _S179;
    float pSH14_x_0 = fTmp1B_1 * _S179;
    float3  dir_n_4 = make_float3 (x_9, y_8, z_6);
    float3  v_dir_n_4 = make_float3 (v_x_0 + dot_0(_S192, make_float3 (-0.59004360437393188f * (fS1_1 + x_9 * fS1_x_0 + _S212)) * *_S185 + make_float3 (-0.59004360437393188f * (fC1_1 + _S211 - y_8 * fS1_x_0)) * *_S191 + make_float3 (fTmp1B_1 * fS1_x_0) * *_S186 + make_float3 (pSH14_x_0) * *_S190 + make_float3 (fTmp0C_1) * *_S189), v_y_0 + dot_0(_S192, make_float3 (-0.59004360437393188f * (_S211 + fC1_1 + y_8 * fC1_y_0)) * *_S185 + make_float3 (-0.59004360437393188f * (x_9 * fC1_y_0 - fS1_1 - _S212)) * *_S191 + make_float3 (pSH14_x_0) * *_S186 + make_float3 (fTmp1B_1 * fC1_y_0) * *_S190 + make_float3 (fTmp0C_1) * *_S187), v_z_0 + dot_0(_S192, make_float3 (5.59764480590820312f * z2_1 - 1.11952900886535645f) * *_S188 + make_float3 (fTmp0C_z_0 * x_9) * *_S189 + make_float3 (fTmp0C_z_0 * y_8) * *_S187 + make_float3 (1.44530570507049561f * fC1_1) * *_S190 + make_float3 (1.44530570507049561f * fS1_1) * *_S186));
    float3  v_viewdir_11 = v_viewdir_10 + (v_dir_n_4 - make_float3 (dot_0(v_dir_n_4, dir_n_4)) * dir_n_4) * make_float3 (inv_norm_7);
    Matrix<float, 3, 3>  _S213 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S214;
    (&_S214)->primal_0 = _S171;
    (&_S214)->differential_0 = _S213;
    float3  _S215 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S216;
    (&_S216)->primal_0 = t_10;
    (&_S216)->differential_0 = _S215;
    s_bwd_prop_mul_0(&_S214, &_S216, v_viewdir_11);
    Matrix<float, 3, 3>  _S217 = transpose_0(_S214.differential_0);
    *v_mean_6 = *v_mean_6 + v_viewdir_11;
    *v_R_6 = *v_R_6 + _S217;
    *v_t_6 = *v_t_6 + _S216.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_atomic(float3  mean_11, Matrix<float, 3, 3>  R_11, float3  t_11, float3  coeff_dc_11, float3  * coeffs_11, float3  v_colors_7, float3  * v_coeff_dc_7, float3  * v_coeffs_7, float3  * v_mean_7, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    Matrix<float, 3, 3>  _S218 = transpose_0(R_11);
    float3  _S219 = mean_11 + mul_0(_S218, t_11);
    float _S220 = _S219.x;
    float _S221 = _S219.y;
    float _S222 = _S219.z;
    float inv_norm_8 = (F32_rsqrt((_S220 * _S220 + _S221 * _S221 + _S222 * _S222)));
    float x_10 = _S220 * inv_norm_8;
    float y_9 = _S221 * inv_norm_8;
    float z_7 = _S222 * inv_norm_8;
    float3  * _S223 = coeffs_11 + int(0);
    float3  * _S224 = coeffs_11 + int(1);
    float3  * _S225 = coeffs_11 + int(2);
    float z2_2 = z_7 * z_7;
    float fTmp0B_5 = -1.09254848957061768f * z_7;
    float fC1_2 = x_10 * x_10 - y_9 * y_9;
    float _S226 = 2.0f * x_10;
    float fS1_2 = _S226 * y_9;
    float pSH6_3 = 0.94617468118667603f * z2_2 - 0.31539157032966614f;
    float pSH7_3 = fTmp0B_5 * x_10;
    float pSH5_3 = fTmp0B_5 * y_9;
    float pSH8_3 = 0.54627424478530884f * fC1_2;
    float pSH4_3 = 0.54627424478530884f * fS1_2;
    float3  * _S227 = coeffs_11 + int(3);
    float3  * _S228 = coeffs_11 + int(4);
    float3  * _S229 = coeffs_11 + int(5);
    float3  * _S230 = coeffs_11 + int(6);
    float3  * _S231 = coeffs_11 + int(7);
    float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_7;
    float pSH12_1 = z_7 * (1.86588168144226074f * z2_2 - 1.11952900886535645f);
    float pSH13_1 = fTmp0C_2 * x_10;
    float pSH11_1 = fTmp0C_2 * y_9;
    float pSH14_1 = fTmp1B_2 * fC1_2;
    float pSH10_1 = fTmp1B_2 * fS1_2;
    float pSH15_1 = -0.59004360437393188f * (x_10 * fC1_2 - y_9 * fS1_2);
    float pSH9_1 = -0.59004360437393188f * (x_10 * fS1_2 + y_9 * fC1_2);
    float3  * _S232 = coeffs_11 + int(8);
    float3  * _S233 = coeffs_11 + int(9);
    float3  * _S234 = coeffs_11 + int(10);
    float3  * _S235 = coeffs_11 + int(11);
    float3  * _S236 = coeffs_11 + int(12);
    float3  * _S237 = coeffs_11 + int(13);
    float3  * _S238 = coeffs_11 + int(14);
    float3  colors_7 = make_float3 (0.282094806432724f) * coeff_dc_11 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * *_S223 + make_float3 (z_7) * *_S224 - make_float3 (x_10) * *_S225) + (make_float3 (pSH4_3) * *_S227 + make_float3 (pSH5_3) * *_S228 + make_float3 (pSH6_3) * *_S229 + make_float3 (pSH7_3) * *_S230 + make_float3 (pSH8_3) * *_S231) + (make_float3 (pSH9_1) * *_S232 + make_float3 (pSH10_1) * *_S233 + make_float3 (pSH11_1) * *_S234 + make_float3 (pSH12_1) * *_S235 + make_float3 (pSH13_1) * *_S236 + make_float3 (pSH14_1) * *_S237 + make_float3 (pSH15_1) * *_S238);
    float3  _S239 = v_colors_7 * make_float3 (float((colors_7.x) >= -0.5f), float((colors_7.y) >= -0.5f), float((colors_7.z) >= -0.5f));
    float3  v_viewdir_12 = {};
    *v_coeff_dc_7 = *v_coeff_dc_7 + make_float3 (0.282094806432724f) * _S239;
    float3  temp_11 = make_float3 (-0.48860251903533936f * y_9) * _S239;
    float _S240 = dot_0(temp_11, temp_11);
    bool _S241;
    if((F32_isfinite((_S240))))
    {
        _S241 = _S240 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S242 = v_coeffs_7 + int(0);
        float _S243 = atomicAdd(&(_S242->x), temp_11.x);
        float _S244 = atomicAdd(&(_S242->y), temp_11.y);
        float _S245 = atomicAdd(&(_S242->z), temp_11.z);
    }
    float3  temp_12 = make_float3 (0.48860251903533936f * z_7) * _S239;
    float _S246 = dot_0(temp_12, temp_12);
    if((F32_isfinite((_S246))))
    {
        _S241 = _S246 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S247 = v_coeffs_7 + int(1);
        float _S248 = atomicAdd(&(_S247->x), temp_12.x);
        float _S249 = atomicAdd(&(_S247->y), temp_12.y);
        float _S250 = atomicAdd(&(_S247->z), temp_12.z);
    }
    float3  temp_13 = make_float3 (-0.48860251903533936f * x_10) * _S239;
    float _S251 = dot_0(temp_13, temp_13);
    if((F32_isfinite((_S251))))
    {
        _S241 = _S251 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S252 = v_coeffs_7 + int(2);
        float _S253 = atomicAdd(&(_S252->x), temp_13.x);
        float _S254 = atomicAdd(&(_S252->y), temp_13.y);
        float _S255 = atomicAdd(&(_S252->z), temp_13.z);
    }
    float _S256 = -0.48860251903533936f * dot_0(*_S225, _S239);
    float _S257 = -0.48860251903533936f * dot_0(*_S223, _S239);
    float _S258 = 0.48860251903533936f * dot_0(*_S224, _S239);
    float3  temp_14 = make_float3 (pSH4_3) * _S239;
    float _S259 = dot_0(temp_14, temp_14);
    if((F32_isfinite((_S259))))
    {
        _S241 = _S259 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S260 = v_coeffs_7 + int(3);
        float _S261 = atomicAdd(&(_S260->x), temp_14.x);
        float _S262 = atomicAdd(&(_S260->y), temp_14.y);
        float _S263 = atomicAdd(&(_S260->z), temp_14.z);
    }
    float3  temp_15 = make_float3 (pSH5_3) * _S239;
    float _S264 = dot_0(temp_15, temp_15);
    if((F32_isfinite((_S264))))
    {
        _S241 = _S264 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S265 = v_coeffs_7 + int(4);
        float _S266 = atomicAdd(&(_S265->x), temp_15.x);
        float _S267 = atomicAdd(&(_S265->y), temp_15.y);
        float _S268 = atomicAdd(&(_S265->z), temp_15.z);
    }
    float3  temp_16 = make_float3 (pSH6_3) * _S239;
    float _S269 = dot_0(temp_16, temp_16);
    if((F32_isfinite((_S269))))
    {
        _S241 = _S269 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S270 = v_coeffs_7 + int(5);
        float _S271 = atomicAdd(&(_S270->x), temp_16.x);
        float _S272 = atomicAdd(&(_S270->y), temp_16.y);
        float _S273 = atomicAdd(&(_S270->z), temp_16.z);
    }
    float3  temp_17 = make_float3 (pSH7_3) * _S239;
    float _S274 = dot_0(temp_17, temp_17);
    if((F32_isfinite((_S274))))
    {
        _S241 = _S274 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S275 = v_coeffs_7 + int(6);
        float _S276 = atomicAdd(&(_S275->x), temp_17.x);
        float _S277 = atomicAdd(&(_S275->y), temp_17.y);
        float _S278 = atomicAdd(&(_S275->z), temp_17.z);
    }
    float3  temp_18 = make_float3 (pSH8_3) * _S239;
    float _S279 = dot_0(temp_18, temp_18);
    if((F32_isfinite((_S279))))
    {
        _S241 = _S279 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S280 = v_coeffs_7 + int(7);
        float _S281 = atomicAdd(&(_S280->x), temp_18.x);
        float _S282 = atomicAdd(&(_S280->y), temp_18.y);
        float _S283 = atomicAdd(&(_S280->z), temp_18.z);
    }
    float fC1_y_1 = -2.0f * y_9;
    float fS1_x_1 = 2.0f * y_9;
    float pSH8_x_3 = 0.54627424478530884f * _S226;
    float v_x_1 = _S256 + dot_0(_S239, make_float3 (0.54627424478530884f * fS1_x_1) * *_S227 + make_float3 (pSH8_x_3) * *_S231 + make_float3 (fTmp0B_5) * *_S230);
    float v_y_1 = _S257 + dot_0(_S239, make_float3 (pSH8_x_3) * *_S227 + make_float3 (0.54627424478530884f * fC1_y_1) * *_S231 + make_float3 (fTmp0B_5) * *_S228);
    float v_z_1 = _S258 + dot_0(_S239, make_float3 (1.89234936237335205f * z_7) * *_S229 + make_float3 (-1.09254848957061768f * x_10) * *_S230 + make_float3 (-1.09254848957061768f * y_9) * *_S228);
    float3  temp_19 = make_float3 (pSH9_1) * _S239;
    float _S284 = dot_0(temp_19, temp_19);
    if((F32_isfinite((_S284))))
    {
        _S241 = _S284 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S285 = v_coeffs_7 + int(8);
        float _S286 = atomicAdd(&(_S285->x), temp_19.x);
        float _S287 = atomicAdd(&(_S285->y), temp_19.y);
        float _S288 = atomicAdd(&(_S285->z), temp_19.z);
    }
    float3  temp_20 = make_float3 (pSH10_1) * _S239;
    float _S289 = dot_0(temp_20, temp_20);
    if((F32_isfinite((_S289))))
    {
        _S241 = _S289 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S290 = v_coeffs_7 + int(9);
        float _S291 = atomicAdd(&(_S290->x), temp_20.x);
        float _S292 = atomicAdd(&(_S290->y), temp_20.y);
        float _S293 = atomicAdd(&(_S290->z), temp_20.z);
    }
    float3  temp_21 = make_float3 (pSH11_1) * _S239;
    float _S294 = dot_0(temp_21, temp_21);
    if((F32_isfinite((_S294))))
    {
        _S241 = _S294 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S295 = v_coeffs_7 + int(10);
        float _S296 = atomicAdd(&(_S295->x), temp_21.x);
        float _S297 = atomicAdd(&(_S295->y), temp_21.y);
        float _S298 = atomicAdd(&(_S295->z), temp_21.z);
    }
    float3  temp_22 = make_float3 (pSH12_1) * _S239;
    float _S299 = dot_0(temp_22, temp_22);
    if((F32_isfinite((_S299))))
    {
        _S241 = _S299 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S300 = v_coeffs_7 + int(11);
        float _S301 = atomicAdd(&(_S300->x), temp_22.x);
        float _S302 = atomicAdd(&(_S300->y), temp_22.y);
        float _S303 = atomicAdd(&(_S300->z), temp_22.z);
    }
    float3  temp_23 = make_float3 (pSH13_1) * _S239;
    float _S304 = dot_0(temp_23, temp_23);
    if((F32_isfinite((_S304))))
    {
        _S241 = _S304 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S305 = v_coeffs_7 + int(12);
        float _S306 = atomicAdd(&(_S305->x), temp_23.x);
        float _S307 = atomicAdd(&(_S305->y), temp_23.y);
        float _S308 = atomicAdd(&(_S305->z), temp_23.z);
    }
    float3  temp_24 = make_float3 (pSH14_1) * _S239;
    float _S309 = dot_0(temp_24, temp_24);
    if((F32_isfinite((_S309))))
    {
        _S241 = _S309 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S310 = v_coeffs_7 + int(13);
        float _S311 = atomicAdd(&(_S310->x), temp_24.x);
        float _S312 = atomicAdd(&(_S310->y), temp_24.y);
        float _S313 = atomicAdd(&(_S310->z), temp_24.z);
    }
    float3  temp_25 = make_float3 (pSH15_1) * _S239;
    float _S314 = dot_0(temp_25, temp_25);
    if((F32_isfinite((_S314))))
    {
        _S241 = _S314 != 0.0f;
    }
    else
    {
        _S241 = false;
    }
    if(_S241)
    {
        float3  * _S315 = v_coeffs_7 + int(14);
        float _S316 = atomicAdd(&(_S315->x), temp_25.x);
        float _S317 = atomicAdd(&(_S315->y), temp_25.y);
        float _S318 = atomicAdd(&(_S315->z), temp_25.z);
    }
    float fTmp0C_z_1 = -4.57045793533325195f * z_7;
    float _S319 = x_10 * _S226;
    float _S320 = y_9 * _S226;
    float pSH14_x_1 = fTmp1B_2 * _S226;
    float3  dir_n_5 = make_float3 (x_10, y_9, z_7);
    float3  v_dir_n_5 = make_float3 (v_x_1 + dot_0(_S239, make_float3 (-0.59004360437393188f * (fS1_2 + x_10 * fS1_x_1 + _S320)) * *_S232 + make_float3 (-0.59004360437393188f * (fC1_2 + _S319 - y_9 * fS1_x_1)) * *_S238 + make_float3 (fTmp1B_2 * fS1_x_1) * *_S233 + make_float3 (pSH14_x_1) * *_S237 + make_float3 (fTmp0C_2) * *_S236), v_y_1 + dot_0(_S239, make_float3 (-0.59004360437393188f * (_S319 + fC1_2 + y_9 * fC1_y_1)) * *_S232 + make_float3 (-0.59004360437393188f * (x_10 * fC1_y_1 - fS1_2 - _S320)) * *_S238 + make_float3 (pSH14_x_1) * *_S233 + make_float3 (fTmp1B_2 * fC1_y_1) * *_S237 + make_float3 (fTmp0C_2) * *_S234), v_z_1 + dot_0(_S239, make_float3 (5.59764480590820312f * z2_2 - 1.11952900886535645f) * *_S235 + make_float3 (fTmp0C_z_1 * x_10) * *_S236 + make_float3 (fTmp0C_z_1 * y_9) * *_S234 + make_float3 (1.44530570507049561f * fC1_2) * *_S237 + make_float3 (1.44530570507049561f * fS1_2) * *_S233));
    float3  v_viewdir_13 = v_viewdir_12 + (v_dir_n_5 - make_float3 (dot_0(v_dir_n_5, dir_n_5)) * dir_n_5) * make_float3 (inv_norm_8);
    Matrix<float, 3, 3>  _S321 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S322;
    (&_S322)->primal_0 = _S218;
    (&_S322)->differential_0 = _S321;
    float3  _S323 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S324;
    (&_S324)->primal_0 = t_11;
    (&_S324)->differential_0 = _S323;
    s_bwd_prop_mul_0(&_S322, &_S324, v_viewdir_13);
    Matrix<float, 3, 3>  _S325 = transpose_0(_S322.differential_0);
    *v_mean_7 = *v_mean_7 + v_viewdir_13;
    *v_R_7 = *v_R_7 + _S325;
    *v_t_7 = *v_t_7 + _S324.differential_0;
    return;
}

inline __device__ float3  sh4_to_color(float3  mean_12, Matrix<float, 3, 3>  R_12, float3  t_12, float3  coeff_dc_12, float3  * coeffs_12)
{
    float3  _S326 = mean_12 + mul_0(transpose_0(R_12), t_12);
    float _S327 = _S326.x;
    float _S328 = _S326.y;
    float _S329 = _S326.z;
    float inv_norm_9 = (F32_rsqrt((_S327 * _S327 + _S328 * _S328 + _S329 * _S329)));
    float x_11 = _S327 * inv_norm_9;
    float y_10 = _S328 * inv_norm_9;
    float z_8 = _S329 * inv_norm_9;
    float z2_3 = z_8 * z_8;
    float fTmp0B_6 = -1.09254848957061768f * z_8;
    float fC1_3 = x_11 * x_11 - y_10 * y_10;
    float fS1_3 = 2.0f * x_11 * y_10;
    float pSH6_4 = 0.94617468118667603f * z2_3 - 0.31539157032966614f;
    float fTmp0C_3 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_8;
    float fC2_0 = x_11 * fC1_3 - y_10 * fS1_3;
    float fS2_0 = x_11 * fS1_3 + y_10 * fC1_3;
    float pSH12_2 = z_8 * (1.86588168144226074f * z2_3 - 1.11952900886535645f);
    float fTmp0D_0 = z_8 * (-4.68332576751708984f * z2_3 + 2.00713968276977539f);
    float fTmp1C_0 = 3.31161141395568848f * z2_3 - 0.47308734059333801f;
    float fTmp2B_0 = -1.77013075351715088f * z_8;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_12 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * *(coeffs_12 + int(0)) + make_float3 (z_8) * *(coeffs_12 + int(1)) - make_float3 (x_11) * *(coeffs_12 + int(2))) + (make_float3 (0.54627424478530884f * fS1_3) * *(coeffs_12 + int(3)) + make_float3 (fTmp0B_6 * y_10) * *(coeffs_12 + int(4)) + make_float3 (pSH6_4) * *(coeffs_12 + int(5)) + make_float3 (fTmp0B_6 * x_11) * *(coeffs_12 + int(6)) + make_float3 (0.54627424478530884f * fC1_3) * *(coeffs_12 + int(7))) + (make_float3 (-0.59004360437393188f * fS2_0) * *(coeffs_12 + int(8)) + make_float3 (fTmp1B_3 * fS1_3) * *(coeffs_12 + int(9)) + make_float3 (fTmp0C_3 * y_10) * *(coeffs_12 + int(10)) + make_float3 (pSH12_2) * *(coeffs_12 + int(11)) + make_float3 (fTmp0C_3 * x_11) * *(coeffs_12 + int(12)) + make_float3 (fTmp1B_3 * fC1_3) * *(coeffs_12 + int(13)) + make_float3 (-0.59004360437393188f * fC2_0) * *(coeffs_12 + int(14))) + (make_float3 (0.62583571672439575f * (x_11 * fS2_0 + y_10 * fC2_0)) * *(coeffs_12 + int(15)) + make_float3 (fTmp2B_0 * fS2_0) * *(coeffs_12 + int(16)) + make_float3 (fTmp1C_0 * fS1_3) * *(coeffs_12 + int(17)) + make_float3 (fTmp0D_0 * y_10) * *(coeffs_12 + int(18)) + make_float3 (1.9843134880065918f * z_8 * pSH12_2 - 1.00623059272766113f * pSH6_4) * *(coeffs_12 + int(19)) + make_float3 (fTmp0D_0 * x_11) * *(coeffs_12 + int(20)) + make_float3 (fTmp1C_0 * fC1_3) * *(coeffs_12 + int(21)) + make_float3 (fTmp2B_0 * fC2_0) * *(coeffs_12 + int(22)) + make_float3 (0.62583571672439575f * (x_11 * fC2_0 - y_10 * fS2_0)) * *(coeffs_12 + int(23))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh4_to_color_vjp_inplace(float3  mean_13, Matrix<float, 3, 3>  R_13, float3  t_13, float3  coeff_dc_13, float3  * coeffs_13, float3  v_colors_8, float3  * v_coeff_dc_8, float3  * v_coeffs_8, float3  * v_mean_8, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    Matrix<float, 3, 3>  _S330 = transpose_0(R_13);
    float3  _S331 = mean_13 + mul_0(_S330, t_13);
    float _S332 = _S331.x;
    float _S333 = _S331.y;
    float _S334 = _S331.z;
    float inv_norm_10 = (F32_rsqrt((_S332 * _S332 + _S333 * _S333 + _S334 * _S334)));
    float x_12 = _S332 * inv_norm_10;
    float y_11 = _S333 * inv_norm_10;
    float z_9 = _S334 * inv_norm_10;
    float3  * _S335 = coeffs_13 + int(0);
    float3  * _S336 = coeffs_13 + int(1);
    float3  * _S337 = coeffs_13 + int(2);
    float z2_4 = z_9 * z_9;
    float fTmp0B_7 = -1.09254848957061768f * z_9;
    float fC1_4 = x_12 * x_12 - y_11 * y_11;
    float _S338 = 2.0f * x_12;
    float fS1_4 = _S338 * y_11;
    float pSH6_5 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float pSH7_4 = fTmp0B_7 * x_12;
    float pSH5_4 = fTmp0B_7 * y_11;
    float pSH8_4 = 0.54627424478530884f * fC1_4;
    float pSH4_4 = 0.54627424478530884f * fS1_4;
    float3  * _S339 = coeffs_13 + int(3);
    float3  * _S340 = coeffs_13 + int(4);
    float3  * _S341 = coeffs_13 + int(5);
    float3  * _S342 = coeffs_13 + int(6);
    float3  * _S343 = coeffs_13 + int(7);
    float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_9;
    float fC2_1 = x_12 * fC1_4 - y_11 * fS1_4;
    float fS2_1 = x_12 * fS1_4 + y_11 * fC1_4;
    float pSH12_3 = z_9 * (1.86588168144226074f * z2_4 - 1.11952900886535645f);
    float pSH13_2 = fTmp0C_4 * x_12;
    float pSH11_2 = fTmp0C_4 * y_11;
    float pSH14_2 = fTmp1B_4 * fC1_4;
    float pSH10_2 = fTmp1B_4 * fS1_4;
    float pSH15_2 = -0.59004360437393188f * fC2_1;
    float pSH9_2 = -0.59004360437393188f * fS2_1;
    float3  * _S344 = coeffs_13 + int(8);
    float3  * _S345 = coeffs_13 + int(9);
    float3  * _S346 = coeffs_13 + int(10);
    float3  * _S347 = coeffs_13 + int(11);
    float3  * _S348 = coeffs_13 + int(12);
    float3  * _S349 = coeffs_13 + int(13);
    float3  * _S350 = coeffs_13 + int(14);
    float fTmp0D_1 = z_9 * (-4.68332576751708984f * z2_4 + 2.00713968276977539f);
    float fTmp1C_1 = 3.31161141395568848f * z2_4 - 0.47308734059333801f;
    float fTmp2B_1 = -1.77013075351715088f * z_9;
    float _S351 = 1.9843134880065918f * z_9 * pSH12_3;
    float pSH21_0 = fTmp0D_1 * x_12;
    float pSH19_0 = fTmp0D_1 * y_11;
    float pSH22_0 = fTmp1C_1 * fC1_4;
    float pSH18_0 = fTmp1C_1 * fS1_4;
    float pSH23_0 = fTmp2B_1 * fC2_1;
    float pSH17_0 = fTmp2B_1 * fS2_1;
    float pSH24_0 = 0.62583571672439575f * (x_12 * fC2_1 - y_11 * fS2_1);
    float pSH16_0 = 0.62583571672439575f * (x_12 * fS2_1 + y_11 * fC2_1);
    float3  * _S352 = coeffs_13 + int(15);
    float3  * _S353 = coeffs_13 + int(16);
    float3  * _S354 = coeffs_13 + int(17);
    float3  * _S355 = coeffs_13 + int(18);
    float3  * _S356 = coeffs_13 + int(19);
    float3  * _S357 = coeffs_13 + int(20);
    float3  * _S358 = coeffs_13 + int(21);
    float3  * _S359 = coeffs_13 + int(22);
    float3  * _S360 = coeffs_13 + int(23);
    float3  colors_8 = make_float3 (0.282094806432724f) * coeff_dc_13 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * *_S335 + make_float3 (z_9) * *_S336 - make_float3 (x_12) * *_S337) + (make_float3 (pSH4_4) * *_S339 + make_float3 (pSH5_4) * *_S340 + make_float3 (pSH6_5) * *_S341 + make_float3 (pSH7_4) * *_S342 + make_float3 (pSH8_4) * *_S343) + (make_float3 (pSH9_2) * *_S344 + make_float3 (pSH10_2) * *_S345 + make_float3 (pSH11_2) * *_S346 + make_float3 (pSH12_3) * *_S347 + make_float3 (pSH13_2) * *_S348 + make_float3 (pSH14_2) * *_S349 + make_float3 (pSH15_2) * *_S350) + (make_float3 (pSH16_0) * *_S352 + make_float3 (pSH17_0) * *_S353 + make_float3 (pSH18_0) * *_S354 + make_float3 (pSH19_0) * *_S355 + make_float3 (_S351 - 1.00623059272766113f * pSH6_5) * *_S356 + make_float3 (pSH21_0) * *_S357 + make_float3 (pSH22_0) * *_S358 + make_float3 (pSH23_0) * *_S359 + make_float3 (pSH24_0) * *_S360);
    float3  _S361 = v_colors_8 * make_float3 (float((colors_8.x) >= -0.5f), float((colors_8.y) >= -0.5f), float((colors_8.z) >= -0.5f));
    float3  v_viewdir_14 = {};
    *v_coeff_dc_8 = *v_coeff_dc_8 + make_float3 (0.282094806432724f) * _S361;
    float3  * _S362 = v_coeffs_8 + int(0);
    *_S362 = *_S362 + make_float3 (-0.48860251903533936f * y_11) * _S361;
    float3  * _S363 = v_coeffs_8 + int(1);
    *_S363 = *_S363 + make_float3 (0.48860251903533936f * z_9) * _S361;
    float3  * _S364 = v_coeffs_8 + int(2);
    *_S364 = *_S364 + make_float3 (-0.48860251903533936f * x_12) * _S361;
    float _S365 = -0.48860251903533936f * dot_0(*_S337, _S361);
    float _S366 = -0.48860251903533936f * dot_0(*_S335, _S361);
    float _S367 = 0.48860251903533936f * dot_0(*_S336, _S361);
    float3  * _S368 = v_coeffs_8 + int(3);
    *_S368 = *_S368 + make_float3 (pSH4_4) * _S361;
    float3  * _S369 = v_coeffs_8 + int(4);
    *_S369 = *_S369 + make_float3 (pSH5_4) * _S361;
    float3  * _S370 = v_coeffs_8 + int(5);
    *_S370 = *_S370 + make_float3 (pSH6_5) * _S361;
    float3  * _S371 = v_coeffs_8 + int(6);
    *_S371 = *_S371 + make_float3 (pSH7_4) * _S361;
    float3  * _S372 = v_coeffs_8 + int(7);
    *_S372 = *_S372 + make_float3 (pSH8_4) * _S361;
    float fC1_y_2 = -2.0f * y_11;
    float fS1_x_2 = 2.0f * y_11;
    float pSH6_z_0 = 1.89234936237335205f * z_9;
    float pSH8_x_4 = 0.54627424478530884f * _S338;
    float v_x_2 = _S365 + dot_0(_S361, make_float3 (0.54627424478530884f * fS1_x_2) * *_S339 + make_float3 (pSH8_x_4) * *_S343 + make_float3 (fTmp0B_7) * *_S342);
    float v_y_2 = _S366 + dot_0(_S361, make_float3 (pSH8_x_4) * *_S339 + make_float3 (0.54627424478530884f * fC1_y_2) * *_S343 + make_float3 (fTmp0B_7) * *_S340);
    float v_z_2 = _S367 + dot_0(_S361, make_float3 (pSH6_z_0) * *_S341 + make_float3 (-1.09254848957061768f * x_12) * *_S342 + make_float3 (-1.09254848957061768f * y_11) * *_S340);
    float3  * _S373 = v_coeffs_8 + int(8);
    *_S373 = *_S373 + make_float3 (pSH9_2) * _S361;
    float3  * _S374 = v_coeffs_8 + int(9);
    *_S374 = *_S374 + make_float3 (pSH10_2) * _S361;
    float3  * _S375 = v_coeffs_8 + int(10);
    *_S375 = *_S375 + make_float3 (pSH11_2) * _S361;
    float3  * _S376 = v_coeffs_8 + int(11);
    *_S376 = *_S376 + make_float3 (pSH12_3) * _S361;
    float3  * _S377 = v_coeffs_8 + int(12);
    *_S377 = *_S377 + make_float3 (pSH13_2) * _S361;
    float3  * _S378 = v_coeffs_8 + int(13);
    *_S378 = *_S378 + make_float3 (pSH14_2) * _S361;
    float3  * _S379 = v_coeffs_8 + int(14);
    *_S379 = *_S379 + make_float3 (pSH15_2) * _S361;
    float fTmp0C_z_2 = -4.57045793533325195f * z_9;
    float _S380 = x_12 * _S338;
    float fC2_x_0 = fC1_4 + _S380 - y_11 * fS1_x_2;
    float _S381 = y_11 * _S338;
    float fC2_y_0 = x_12 * fC1_y_2 - fS1_4 - _S381;
    float fS2_x_0 = fS1_4 + x_12 * fS1_x_2 + _S381;
    float fS2_y_0 = _S380 + fC1_4 + y_11 * fC1_y_2;
    float pSH12_z_0 = 5.59764480590820312f * z2_4 - 1.11952900886535645f;
    float pSH14_x_2 = fTmp1B_4 * _S338;
    float v_x_3 = v_x_2 + dot_0(_S361, make_float3 (-0.59004360437393188f * fS2_x_0) * *_S344 + make_float3 (-0.59004360437393188f * fC2_x_0) * *_S350 + make_float3 (fTmp1B_4 * fS1_x_2) * *_S345 + make_float3 (pSH14_x_2) * *_S349 + make_float3 (fTmp0C_4) * *_S348);
    float v_y_3 = v_y_2 + dot_0(_S361, make_float3 (-0.59004360437393188f * fS2_y_0) * *_S344 + make_float3 (-0.59004360437393188f * fC2_y_0) * *_S350 + make_float3 (pSH14_x_2) * *_S345 + make_float3 (fTmp1B_4 * fC1_y_2) * *_S349 + make_float3 (fTmp0C_4) * *_S346);
    float v_z_3 = v_z_2 + dot_0(_S361, make_float3 (pSH12_z_0) * *_S347 + make_float3 (fTmp0C_z_2 * x_12) * *_S348 + make_float3 (fTmp0C_z_2 * y_11) * *_S346 + make_float3 (1.44530570507049561f * fC1_4) * *_S349 + make_float3 (1.44530570507049561f * fS1_4) * *_S345);
    float pSH20_0 = _S351 + -1.00623059272766113f * pSH6_5;
    float3  * _S382 = v_coeffs_8 + int(15);
    *_S382 = *_S382 + make_float3 (pSH16_0) * _S361;
    float3  * _S383 = v_coeffs_8 + int(16);
    *_S383 = *_S383 + make_float3 (pSH17_0) * _S361;
    float3  * _S384 = v_coeffs_8 + int(17);
    *_S384 = *_S384 + make_float3 (pSH18_0) * _S361;
    float3  * _S385 = v_coeffs_8 + int(18);
    *_S385 = *_S385 + make_float3 (pSH19_0) * _S361;
    float3  * _S386 = v_coeffs_8 + int(19);
    *_S386 = *_S386 + make_float3 (pSH20_0) * _S361;
    float3  * _S387 = v_coeffs_8 + int(20);
    *_S387 = *_S387 + make_float3 (pSH21_0) * _S361;
    float3  * _S388 = v_coeffs_8 + int(21);
    *_S388 = *_S388 + make_float3 (pSH22_0) * _S361;
    float3  * _S389 = v_coeffs_8 + int(22);
    *_S389 = *_S389 + make_float3 (pSH23_0) * _S361;
    float3  * _S390 = v_coeffs_8 + int(23);
    *_S390 = *_S390 + make_float3 (pSH24_0) * _S361;
    float fTmp0D_z_0 = -14.04997730255126953f * z2_4 + 2.00713968276977539f;
    float fTmp1C_z_0 = 6.62322282791137695f * z_9;
    float pSH22_x_0 = fTmp1C_1 * _S338;
    float3  dir_n_6 = make_float3 (x_12, y_11, z_9);
    float3  v_dir_n_6 = make_float3 (v_x_3 + dot_0(_S361, make_float3 (0.62583571672439575f * (fS2_1 + y_11 * fC2_x_0 + x_12 * fS2_x_0)) * *_S352 + make_float3 (0.62583571672439575f * (fC2_1 + x_12 * fC2_x_0 - y_11 * fS2_x_0)) * *_S360 + make_float3 (fTmp2B_1 * fS2_x_0) * *_S353 + make_float3 (fTmp2B_1 * fC2_x_0) * *_S359 + make_float3 (fTmp1C_1 * fS1_x_2) * *_S354 + make_float3 (pSH22_x_0) * *_S358 + make_float3 (fTmp0D_1) * *_S357), v_y_3 + dot_0(_S361, make_float3 (0.62583571672439575f * (x_12 * fS2_y_0 + fC2_1 + y_11 * fC2_y_0)) * *_S352 + make_float3 (0.62583571672439575f * (x_12 * fC2_y_0 - fS2_1 - y_11 * fS2_y_0)) * *_S360 + make_float3 (fTmp2B_1 * fS2_y_0) * *_S353 + make_float3 (fTmp2B_1 * fC2_y_0) * *_S359 + make_float3 (pSH22_x_0) * *_S354 + make_float3 (fTmp1C_1 * fC1_y_2) * *_S358 + make_float3 (fTmp0D_1) * *_S355), v_z_3 + dot_0(_S361, make_float3 (1.9843134880065918f * (pSH12_3 + z_9 * pSH12_z_0) + -1.00623059272766113f * pSH6_z_0) * *_S356 + make_float3 (fTmp0D_z_0 * x_12) * *_S357 + make_float3 (fTmp0D_z_0 * y_11) * *_S355 + make_float3 (fTmp1C_z_0 * fC1_4) * *_S358 + make_float3 (fTmp1C_z_0 * fS1_4) * *_S354 + make_float3 (-1.77013075351715088f * fC2_1) * *_S359 + make_float3 (-1.77013075351715088f * fS2_1) * *_S353));
    float3  v_viewdir_15 = v_viewdir_14 + (v_dir_n_6 - make_float3 (dot_0(v_dir_n_6, dir_n_6)) * dir_n_6) * make_float3 (inv_norm_10);
    Matrix<float, 3, 3>  _S391 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S392;
    (&_S392)->primal_0 = _S330;
    (&_S392)->differential_0 = _S391;
    float3  _S393 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S394;
    (&_S394)->primal_0 = t_13;
    (&_S394)->differential_0 = _S393;
    s_bwd_prop_mul_0(&_S392, &_S394, v_viewdir_15);
    Matrix<float, 3, 3>  _S395 = transpose_0(_S392.differential_0);
    *v_mean_8 = *v_mean_8 + v_viewdir_15;
    *v_R_8 = *v_R_8 + _S395;
    *v_t_8 = *v_t_8 + _S394.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_atomic(float3  mean_14, Matrix<float, 3, 3>  R_14, float3  t_14, float3  coeff_dc_14, float3  * coeffs_14, float3  v_colors_9, float3  * v_coeff_dc_9, float3  * v_coeffs_9, float3  * v_mean_9, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    Matrix<float, 3, 3>  _S396 = transpose_0(R_14);
    float3  _S397 = mean_14 + mul_0(_S396, t_14);
    float _S398 = _S397.x;
    float _S399 = _S397.y;
    float _S400 = _S397.z;
    float inv_norm_11 = (F32_rsqrt((_S398 * _S398 + _S399 * _S399 + _S400 * _S400)));
    float x_13 = _S398 * inv_norm_11;
    float y_12 = _S399 * inv_norm_11;
    float z_10 = _S400 * inv_norm_11;
    float3  * _S401 = coeffs_14 + int(0);
    float3  * _S402 = coeffs_14 + int(1);
    float3  * _S403 = coeffs_14 + int(2);
    float z2_5 = z_10 * z_10;
    float fTmp0B_8 = -1.09254848957061768f * z_10;
    float fC1_5 = x_13 * x_13 - y_12 * y_12;
    float _S404 = 2.0f * x_13;
    float fS1_5 = _S404 * y_12;
    float pSH6_6 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float pSH7_5 = fTmp0B_8 * x_13;
    float pSH5_5 = fTmp0B_8 * y_12;
    float pSH8_5 = 0.54627424478530884f * fC1_5;
    float pSH4_5 = 0.54627424478530884f * fS1_5;
    float3  * _S405 = coeffs_14 + int(3);
    float3  * _S406 = coeffs_14 + int(4);
    float3  * _S407 = coeffs_14 + int(5);
    float3  * _S408 = coeffs_14 + int(6);
    float3  * _S409 = coeffs_14 + int(7);
    float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_10;
    float fC2_2 = x_13 * fC1_5 - y_12 * fS1_5;
    float fS2_2 = x_13 * fS1_5 + y_12 * fC1_5;
    float pSH12_4 = z_10 * (1.86588168144226074f * z2_5 - 1.11952900886535645f);
    float pSH13_3 = fTmp0C_5 * x_13;
    float pSH11_3 = fTmp0C_5 * y_12;
    float pSH14_3 = fTmp1B_5 * fC1_5;
    float pSH10_3 = fTmp1B_5 * fS1_5;
    float pSH15_3 = -0.59004360437393188f * fC2_2;
    float pSH9_3 = -0.59004360437393188f * fS2_2;
    float3  * _S410 = coeffs_14 + int(8);
    float3  * _S411 = coeffs_14 + int(9);
    float3  * _S412 = coeffs_14 + int(10);
    float3  * _S413 = coeffs_14 + int(11);
    float3  * _S414 = coeffs_14 + int(12);
    float3  * _S415 = coeffs_14 + int(13);
    float3  * _S416 = coeffs_14 + int(14);
    float fTmp0D_2 = z_10 * (-4.68332576751708984f * z2_5 + 2.00713968276977539f);
    float fTmp1C_2 = 3.31161141395568848f * z2_5 - 0.47308734059333801f;
    float fTmp2B_2 = -1.77013075351715088f * z_10;
    float _S417 = 1.9843134880065918f * z_10 * pSH12_4;
    float pSH21_1 = fTmp0D_2 * x_13;
    float pSH19_1 = fTmp0D_2 * y_12;
    float pSH22_1 = fTmp1C_2 * fC1_5;
    float pSH18_1 = fTmp1C_2 * fS1_5;
    float pSH23_1 = fTmp2B_2 * fC2_2;
    float pSH17_1 = fTmp2B_2 * fS2_2;
    float pSH24_1 = 0.62583571672439575f * (x_13 * fC2_2 - y_12 * fS2_2);
    float pSH16_1 = 0.62583571672439575f * (x_13 * fS2_2 + y_12 * fC2_2);
    float3  * _S418 = coeffs_14 + int(15);
    float3  * _S419 = coeffs_14 + int(16);
    float3  * _S420 = coeffs_14 + int(17);
    float3  * _S421 = coeffs_14 + int(18);
    float3  * _S422 = coeffs_14 + int(19);
    float3  * _S423 = coeffs_14 + int(20);
    float3  * _S424 = coeffs_14 + int(21);
    float3  * _S425 = coeffs_14 + int(22);
    float3  * _S426 = coeffs_14 + int(23);
    float3  colors_9 = make_float3 (0.282094806432724f) * coeff_dc_14 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * *_S401 + make_float3 (z_10) * *_S402 - make_float3 (x_13) * *_S403) + (make_float3 (pSH4_5) * *_S405 + make_float3 (pSH5_5) * *_S406 + make_float3 (pSH6_6) * *_S407 + make_float3 (pSH7_5) * *_S408 + make_float3 (pSH8_5) * *_S409) + (make_float3 (pSH9_3) * *_S410 + make_float3 (pSH10_3) * *_S411 + make_float3 (pSH11_3) * *_S412 + make_float3 (pSH12_4) * *_S413 + make_float3 (pSH13_3) * *_S414 + make_float3 (pSH14_3) * *_S415 + make_float3 (pSH15_3) * *_S416) + (make_float3 (pSH16_1) * *_S418 + make_float3 (pSH17_1) * *_S419 + make_float3 (pSH18_1) * *_S420 + make_float3 (pSH19_1) * *_S421 + make_float3 (_S417 - 1.00623059272766113f * pSH6_6) * *_S422 + make_float3 (pSH21_1) * *_S423 + make_float3 (pSH22_1) * *_S424 + make_float3 (pSH23_1) * *_S425 + make_float3 (pSH24_1) * *_S426);
    float3  _S427 = v_colors_9 * make_float3 (float((colors_9.x) >= -0.5f), float((colors_9.y) >= -0.5f), float((colors_9.z) >= -0.5f));
    float3  v_viewdir_16 = {};
    *v_coeff_dc_9 = *v_coeff_dc_9 + make_float3 (0.282094806432724f) * _S427;
    float3  temp_26 = make_float3 (-0.48860251903533936f * y_12) * _S427;
    float _S428 = dot_0(temp_26, temp_26);
    bool _S429;
    if((F32_isfinite((_S428))))
    {
        _S429 = _S428 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S430 = v_coeffs_9 + int(0);
        float _S431 = atomicAdd(&(_S430->x), temp_26.x);
        float _S432 = atomicAdd(&(_S430->y), temp_26.y);
        float _S433 = atomicAdd(&(_S430->z), temp_26.z);
    }
    float3  temp_27 = make_float3 (0.48860251903533936f * z_10) * _S427;
    float _S434 = dot_0(temp_27, temp_27);
    if((F32_isfinite((_S434))))
    {
        _S429 = _S434 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S435 = v_coeffs_9 + int(1);
        float _S436 = atomicAdd(&(_S435->x), temp_27.x);
        float _S437 = atomicAdd(&(_S435->y), temp_27.y);
        float _S438 = atomicAdd(&(_S435->z), temp_27.z);
    }
    float3  temp_28 = make_float3 (-0.48860251903533936f * x_13) * _S427;
    float _S439 = dot_0(temp_28, temp_28);
    if((F32_isfinite((_S439))))
    {
        _S429 = _S439 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S440 = v_coeffs_9 + int(2);
        float _S441 = atomicAdd(&(_S440->x), temp_28.x);
        float _S442 = atomicAdd(&(_S440->y), temp_28.y);
        float _S443 = atomicAdd(&(_S440->z), temp_28.z);
    }
    float _S444 = -0.48860251903533936f * dot_0(*_S403, _S427);
    float _S445 = -0.48860251903533936f * dot_0(*_S401, _S427);
    float _S446 = 0.48860251903533936f * dot_0(*_S402, _S427);
    float3  temp_29 = make_float3 (pSH4_5) * _S427;
    float _S447 = dot_0(temp_29, temp_29);
    if((F32_isfinite((_S447))))
    {
        _S429 = _S447 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S448 = v_coeffs_9 + int(3);
        float _S449 = atomicAdd(&(_S448->x), temp_29.x);
        float _S450 = atomicAdd(&(_S448->y), temp_29.y);
        float _S451 = atomicAdd(&(_S448->z), temp_29.z);
    }
    float3  temp_30 = make_float3 (pSH5_5) * _S427;
    float _S452 = dot_0(temp_30, temp_30);
    if((F32_isfinite((_S452))))
    {
        _S429 = _S452 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S453 = v_coeffs_9 + int(4);
        float _S454 = atomicAdd(&(_S453->x), temp_30.x);
        float _S455 = atomicAdd(&(_S453->y), temp_30.y);
        float _S456 = atomicAdd(&(_S453->z), temp_30.z);
    }
    float3  temp_31 = make_float3 (pSH6_6) * _S427;
    float _S457 = dot_0(temp_31, temp_31);
    if((F32_isfinite((_S457))))
    {
        _S429 = _S457 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S458 = v_coeffs_9 + int(5);
        float _S459 = atomicAdd(&(_S458->x), temp_31.x);
        float _S460 = atomicAdd(&(_S458->y), temp_31.y);
        float _S461 = atomicAdd(&(_S458->z), temp_31.z);
    }
    float3  temp_32 = make_float3 (pSH7_5) * _S427;
    float _S462 = dot_0(temp_32, temp_32);
    if((F32_isfinite((_S462))))
    {
        _S429 = _S462 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S463 = v_coeffs_9 + int(6);
        float _S464 = atomicAdd(&(_S463->x), temp_32.x);
        float _S465 = atomicAdd(&(_S463->y), temp_32.y);
        float _S466 = atomicAdd(&(_S463->z), temp_32.z);
    }
    float3  temp_33 = make_float3 (pSH8_5) * _S427;
    float _S467 = dot_0(temp_33, temp_33);
    if((F32_isfinite((_S467))))
    {
        _S429 = _S467 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S468 = v_coeffs_9 + int(7);
        float _S469 = atomicAdd(&(_S468->x), temp_33.x);
        float _S470 = atomicAdd(&(_S468->y), temp_33.y);
        float _S471 = atomicAdd(&(_S468->z), temp_33.z);
    }
    float fC1_y_3 = -2.0f * y_12;
    float fS1_x_3 = 2.0f * y_12;
    float pSH6_z_1 = 1.89234936237335205f * z_10;
    float pSH8_x_5 = 0.54627424478530884f * _S404;
    float v_x_4 = _S444 + dot_0(_S427, make_float3 (0.54627424478530884f * fS1_x_3) * *_S405 + make_float3 (pSH8_x_5) * *_S409 + make_float3 (fTmp0B_8) * *_S408);
    float v_y_4 = _S445 + dot_0(_S427, make_float3 (pSH8_x_5) * *_S405 + make_float3 (0.54627424478530884f * fC1_y_3) * *_S409 + make_float3 (fTmp0B_8) * *_S406);
    float v_z_4 = _S446 + dot_0(_S427, make_float3 (pSH6_z_1) * *_S407 + make_float3 (-1.09254848957061768f * x_13) * *_S408 + make_float3 (-1.09254848957061768f * y_12) * *_S406);
    float3  temp_34 = make_float3 (pSH9_3) * _S427;
    float _S472 = dot_0(temp_34, temp_34);
    if((F32_isfinite((_S472))))
    {
        _S429 = _S472 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S473 = v_coeffs_9 + int(8);
        float _S474 = atomicAdd(&(_S473->x), temp_34.x);
        float _S475 = atomicAdd(&(_S473->y), temp_34.y);
        float _S476 = atomicAdd(&(_S473->z), temp_34.z);
    }
    float3  temp_35 = make_float3 (pSH10_3) * _S427;
    float _S477 = dot_0(temp_35, temp_35);
    if((F32_isfinite((_S477))))
    {
        _S429 = _S477 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S478 = v_coeffs_9 + int(9);
        float _S479 = atomicAdd(&(_S478->x), temp_35.x);
        float _S480 = atomicAdd(&(_S478->y), temp_35.y);
        float _S481 = atomicAdd(&(_S478->z), temp_35.z);
    }
    float3  temp_36 = make_float3 (pSH11_3) * _S427;
    float _S482 = dot_0(temp_36, temp_36);
    if((F32_isfinite((_S482))))
    {
        _S429 = _S482 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S483 = v_coeffs_9 + int(10);
        float _S484 = atomicAdd(&(_S483->x), temp_36.x);
        float _S485 = atomicAdd(&(_S483->y), temp_36.y);
        float _S486 = atomicAdd(&(_S483->z), temp_36.z);
    }
    float3  temp_37 = make_float3 (pSH12_4) * _S427;
    float _S487 = dot_0(temp_37, temp_37);
    if((F32_isfinite((_S487))))
    {
        _S429 = _S487 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S488 = v_coeffs_9 + int(11);
        float _S489 = atomicAdd(&(_S488->x), temp_37.x);
        float _S490 = atomicAdd(&(_S488->y), temp_37.y);
        float _S491 = atomicAdd(&(_S488->z), temp_37.z);
    }
    float3  temp_38 = make_float3 (pSH13_3) * _S427;
    float _S492 = dot_0(temp_38, temp_38);
    if((F32_isfinite((_S492))))
    {
        _S429 = _S492 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S493 = v_coeffs_9 + int(12);
        float _S494 = atomicAdd(&(_S493->x), temp_38.x);
        float _S495 = atomicAdd(&(_S493->y), temp_38.y);
        float _S496 = atomicAdd(&(_S493->z), temp_38.z);
    }
    float3  temp_39 = make_float3 (pSH14_3) * _S427;
    float _S497 = dot_0(temp_39, temp_39);
    if((F32_isfinite((_S497))))
    {
        _S429 = _S497 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S498 = v_coeffs_9 + int(13);
        float _S499 = atomicAdd(&(_S498->x), temp_39.x);
        float _S500 = atomicAdd(&(_S498->y), temp_39.y);
        float _S501 = atomicAdd(&(_S498->z), temp_39.z);
    }
    float3  temp_40 = make_float3 (pSH15_3) * _S427;
    float _S502 = dot_0(temp_40, temp_40);
    if((F32_isfinite((_S502))))
    {
        _S429 = _S502 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S503 = v_coeffs_9 + int(14);
        float _S504 = atomicAdd(&(_S503->x), temp_40.x);
        float _S505 = atomicAdd(&(_S503->y), temp_40.y);
        float _S506 = atomicAdd(&(_S503->z), temp_40.z);
    }
    float fTmp0C_z_3 = -4.57045793533325195f * z_10;
    float _S507 = x_13 * _S404;
    float fC2_x_1 = fC1_5 + _S507 - y_12 * fS1_x_3;
    float _S508 = y_12 * _S404;
    float fC2_y_1 = x_13 * fC1_y_3 - fS1_5 - _S508;
    float fS2_x_1 = fS1_5 + x_13 * fS1_x_3 + _S508;
    float fS2_y_1 = _S507 + fC1_5 + y_12 * fC1_y_3;
    float pSH12_z_1 = 5.59764480590820312f * z2_5 - 1.11952900886535645f;
    float pSH14_x_3 = fTmp1B_5 * _S404;
    float v_x_5 = v_x_4 + dot_0(_S427, make_float3 (-0.59004360437393188f * fS2_x_1) * *_S410 + make_float3 (-0.59004360437393188f * fC2_x_1) * *_S416 + make_float3 (fTmp1B_5 * fS1_x_3) * *_S411 + make_float3 (pSH14_x_3) * *_S415 + make_float3 (fTmp0C_5) * *_S414);
    float v_y_5 = v_y_4 + dot_0(_S427, make_float3 (-0.59004360437393188f * fS2_y_1) * *_S410 + make_float3 (-0.59004360437393188f * fC2_y_1) * *_S416 + make_float3 (pSH14_x_3) * *_S411 + make_float3 (fTmp1B_5 * fC1_y_3) * *_S415 + make_float3 (fTmp0C_5) * *_S412);
    float v_z_5 = v_z_4 + dot_0(_S427, make_float3 (pSH12_z_1) * *_S413 + make_float3 (fTmp0C_z_3 * x_13) * *_S414 + make_float3 (fTmp0C_z_3 * y_12) * *_S412 + make_float3 (1.44530570507049561f * fC1_5) * *_S415 + make_float3 (1.44530570507049561f * fS1_5) * *_S411);
    float pSH20_1 = _S417 + -1.00623059272766113f * pSH6_6;
    float3  temp_41 = make_float3 (pSH16_1) * _S427;
    float _S509 = dot_0(temp_41, temp_41);
    if((F32_isfinite((_S509))))
    {
        _S429 = _S509 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S510 = v_coeffs_9 + int(15);
        float _S511 = atomicAdd(&(_S510->x), temp_41.x);
        float _S512 = atomicAdd(&(_S510->y), temp_41.y);
        float _S513 = atomicAdd(&(_S510->z), temp_41.z);
    }
    float3  temp_42 = make_float3 (pSH17_1) * _S427;
    float _S514 = dot_0(temp_42, temp_42);
    if((F32_isfinite((_S514))))
    {
        _S429 = _S514 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S515 = v_coeffs_9 + int(16);
        float _S516 = atomicAdd(&(_S515->x), temp_42.x);
        float _S517 = atomicAdd(&(_S515->y), temp_42.y);
        float _S518 = atomicAdd(&(_S515->z), temp_42.z);
    }
    float3  temp_43 = make_float3 (pSH18_1) * _S427;
    float _S519 = dot_0(temp_43, temp_43);
    if((F32_isfinite((_S519))))
    {
        _S429 = _S519 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S520 = v_coeffs_9 + int(17);
        float _S521 = atomicAdd(&(_S520->x), temp_43.x);
        float _S522 = atomicAdd(&(_S520->y), temp_43.y);
        float _S523 = atomicAdd(&(_S520->z), temp_43.z);
    }
    float3  temp_44 = make_float3 (pSH19_1) * _S427;
    float _S524 = dot_0(temp_44, temp_44);
    if((F32_isfinite((_S524))))
    {
        _S429 = _S524 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S525 = v_coeffs_9 + int(18);
        float _S526 = atomicAdd(&(_S525->x), temp_44.x);
        float _S527 = atomicAdd(&(_S525->y), temp_44.y);
        float _S528 = atomicAdd(&(_S525->z), temp_44.z);
    }
    float3  temp_45 = make_float3 (pSH20_1) * _S427;
    float _S529 = dot_0(temp_45, temp_45);
    if((F32_isfinite((_S529))))
    {
        _S429 = _S529 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S530 = v_coeffs_9 + int(19);
        float _S531 = atomicAdd(&(_S530->x), temp_45.x);
        float _S532 = atomicAdd(&(_S530->y), temp_45.y);
        float _S533 = atomicAdd(&(_S530->z), temp_45.z);
    }
    float3  temp_46 = make_float3 (pSH21_1) * _S427;
    float _S534 = dot_0(temp_46, temp_46);
    if((F32_isfinite((_S534))))
    {
        _S429 = _S534 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S535 = v_coeffs_9 + int(20);
        float _S536 = atomicAdd(&(_S535->x), temp_46.x);
        float _S537 = atomicAdd(&(_S535->y), temp_46.y);
        float _S538 = atomicAdd(&(_S535->z), temp_46.z);
    }
    float3  temp_47 = make_float3 (pSH22_1) * _S427;
    float _S539 = dot_0(temp_47, temp_47);
    if((F32_isfinite((_S539))))
    {
        _S429 = _S539 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S540 = v_coeffs_9 + int(21);
        float _S541 = atomicAdd(&(_S540->x), temp_47.x);
        float _S542 = atomicAdd(&(_S540->y), temp_47.y);
        float _S543 = atomicAdd(&(_S540->z), temp_47.z);
    }
    float3  temp_48 = make_float3 (pSH23_1) * _S427;
    float _S544 = dot_0(temp_48, temp_48);
    if((F32_isfinite((_S544))))
    {
        _S429 = _S544 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S545 = v_coeffs_9 + int(22);
        float _S546 = atomicAdd(&(_S545->x), temp_48.x);
        float _S547 = atomicAdd(&(_S545->y), temp_48.y);
        float _S548 = atomicAdd(&(_S545->z), temp_48.z);
    }
    float3  temp_49 = make_float3 (pSH24_1) * _S427;
    float _S549 = dot_0(temp_49, temp_49);
    if((F32_isfinite((_S549))))
    {
        _S429 = _S549 != 0.0f;
    }
    else
    {
        _S429 = false;
    }
    if(_S429)
    {
        float3  * _S550 = v_coeffs_9 + int(23);
        float _S551 = atomicAdd(&(_S550->x), temp_49.x);
        float _S552 = atomicAdd(&(_S550->y), temp_49.y);
        float _S553 = atomicAdd(&(_S550->z), temp_49.z);
    }
    float fTmp0D_z_1 = -14.04997730255126953f * z2_5 + 2.00713968276977539f;
    float fTmp1C_z_1 = 6.62322282791137695f * z_10;
    float pSH22_x_1 = fTmp1C_2 * _S404;
    float3  dir_n_7 = make_float3 (x_13, y_12, z_10);
    float3  v_dir_n_7 = make_float3 (v_x_5 + dot_0(_S427, make_float3 (0.62583571672439575f * (fS2_2 + y_12 * fC2_x_1 + x_13 * fS2_x_1)) * *_S418 + make_float3 (0.62583571672439575f * (fC2_2 + x_13 * fC2_x_1 - y_12 * fS2_x_1)) * *_S426 + make_float3 (fTmp2B_2 * fS2_x_1) * *_S419 + make_float3 (fTmp2B_2 * fC2_x_1) * *_S425 + make_float3 (fTmp1C_2 * fS1_x_3) * *_S420 + make_float3 (pSH22_x_1) * *_S424 + make_float3 (fTmp0D_2) * *_S423), v_y_5 + dot_0(_S427, make_float3 (0.62583571672439575f * (x_13 * fS2_y_1 + fC2_2 + y_12 * fC2_y_1)) * *_S418 + make_float3 (0.62583571672439575f * (x_13 * fC2_y_1 - fS2_2 - y_12 * fS2_y_1)) * *_S426 + make_float3 (fTmp2B_2 * fS2_y_1) * *_S419 + make_float3 (fTmp2B_2 * fC2_y_1) * *_S425 + make_float3 (pSH22_x_1) * *_S420 + make_float3 (fTmp1C_2 * fC1_y_3) * *_S424 + make_float3 (fTmp0D_2) * *_S421), v_z_5 + dot_0(_S427, make_float3 (1.9843134880065918f * (pSH12_4 + z_10 * pSH12_z_1) + -1.00623059272766113f * pSH6_z_1) * *_S422 + make_float3 (fTmp0D_z_1 * x_13) * *_S423 + make_float3 (fTmp0D_z_1 * y_12) * *_S421 + make_float3 (fTmp1C_z_1 * fC1_5) * *_S424 + make_float3 (fTmp1C_z_1 * fS1_5) * *_S420 + make_float3 (-1.77013075351715088f * fC2_2) * *_S425 + make_float3 (-1.77013075351715088f * fS2_2) * *_S419));
    float3  v_viewdir_17 = v_viewdir_16 + (v_dir_n_7 - make_float3 (dot_0(v_dir_n_7, dir_n_7)) * dir_n_7) * make_float3 (inv_norm_11);
    Matrix<float, 3, 3>  _S554 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S555;
    (&_S555)->primal_0 = _S396;
    (&_S555)->differential_0 = _S554;
    float3  _S556 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S557;
    (&_S557)->primal_0 = t_14;
    (&_S557)->differential_0 = _S556;
    s_bwd_prop_mul_0(&_S555, &_S557, v_viewdir_17);
    Matrix<float, 3, 3>  _S558 = transpose_0(_S555.differential_0);
    *v_mean_9 = *v_mean_9 + v_viewdir_17;
    *v_R_9 = *v_R_9 + _S558;
    *v_t_9 = *v_t_9 + _S557.differential_0;
    return;
}

inline __device__ float3  sh0_to_color(float3  mean_15, Matrix<float, 3, 3>  R_15, float3  t_15, float3  coeff_dc_15, float * coeffs_15)
{
    return max_0(make_float3 (0.282094806432724f * coeff_dc_15.x, 0.282094806432724f * coeff_dc_15.y, 0.282094806432724f * coeff_dc_15.z) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh0_to_color_vjp_inplace(float3  mean_16, Matrix<float, 3, 3>  R_16, float3  t_16, float3  coeff_dc_16, float * coeffs_16, float3  v_colors_10, float3  * v_coeff_dc_10, float * v_coeffs_10, float3  * v_mean_10, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_10)
{
    float3  _S559 = v_colors_10 * make_float3 (float((0.282094806432724f * coeff_dc_16.x) >= -0.5f), float((0.282094806432724f * coeff_dc_16.y) >= -0.5f), float((0.282094806432724f * coeff_dc_16.z) >= -0.5f));
    float3  v_viewdir_18 = {};
    *&(v_coeff_dc_10->x) = *&(v_coeff_dc_10->x) + 0.282094806432724f * _S559.x;
    *&(v_coeff_dc_10->y) = *&(v_coeff_dc_10->y) + 0.282094806432724f * _S559.y;
    *&(v_coeff_dc_10->z) = *&(v_coeff_dc_10->z) + 0.282094806432724f * _S559.z;
    Matrix<float, 3, 3>  _S560 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S561;
    (&_S561)->primal_0 = transpose_0(R_16);
    (&_S561)->differential_0 = _S560;
    float3  _S562 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S563;
    (&_S563)->primal_0 = t_16;
    (&_S563)->differential_0 = _S562;
    s_bwd_prop_mul_0(&_S561, &_S563, v_viewdir_18);
    Matrix<float, 3, 3>  _S564 = transpose_0(_S561.differential_0);
    *v_mean_10 = *v_mean_10 + v_viewdir_18;
    *v_R_10 = *v_R_10 + _S564;
    *v_t_10 = *v_t_10 + _S563.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_atomic(float3  mean_17, Matrix<float, 3, 3>  R_17, float3  t_17, float3  coeff_dc_17, float * coeffs_17, float3  v_colors_11, float3  * v_coeff_dc_11, float * v_coeffs_11, float3  * v_mean_11, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_11)
{
    float3  _S565 = v_colors_11 * make_float3 (float((0.282094806432724f * coeff_dc_17.x) >= -0.5f), float((0.282094806432724f * coeff_dc_17.y) >= -0.5f), float((0.282094806432724f * coeff_dc_17.z) >= -0.5f));
    float3  v_viewdir_19 = {};
    *&(v_coeff_dc_11->x) = *&(v_coeff_dc_11->x) + 0.282094806432724f * _S565.x;
    *&(v_coeff_dc_11->y) = *&(v_coeff_dc_11->y) + 0.282094806432724f * _S565.y;
    *&(v_coeff_dc_11->z) = *&(v_coeff_dc_11->z) + 0.282094806432724f * _S565.z;
    Matrix<float, 3, 3>  _S566 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S567;
    (&_S567)->primal_0 = transpose_0(R_17);
    (&_S567)->differential_0 = _S566;
    float3  _S568 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S569;
    (&_S569)->primal_0 = t_17;
    (&_S569)->differential_0 = _S568;
    s_bwd_prop_mul_0(&_S567, &_S569, v_viewdir_19);
    Matrix<float, 3, 3>  _S570 = transpose_0(_S567.differential_0);
    *v_mean_11 = *v_mean_11 + v_viewdir_19;
    *v_R_11 = *v_R_11 + _S570;
    *v_t_11 = *v_t_11 + _S569.differential_0;
    return;
}

inline __device__ float3  sh1_to_color(float3  mean_18, Matrix<float, 3, 3>  R_18, float3  t_18, float3  coeff_dc_18, float * coeffs_18)
{
    float3  _S571 = mean_18 + mul_0(transpose_0(R_18), t_18);
    float _S572 = _S571.x;
    float _S573 = _S571.y;
    float _S574 = _S571.z;
    float inv_norm_12 = (F32_rsqrt((_S572 * _S572 + _S573 * _S573 + _S574 * _S574)));
    float x_14 = _S572 * inv_norm_12;
    float z_11 = _S574 * inv_norm_12;
    float _S575 = - (_S573 * inv_norm_12);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_18.x + 0.48860251903533936f * (_S575 * *(coeffs_18 + int(0)) + z_11 * *(coeffs_18 + int(3)) - x_14 * *(coeffs_18 + int(6))), 0.282094806432724f * coeff_dc_18.y + 0.48860251903533936f * (_S575 * *(coeffs_18 + int(1)) + z_11 * *(coeffs_18 + int(4)) - x_14 * *(coeffs_18 + int(7))), 0.282094806432724f * coeff_dc_18.z + 0.48860251903533936f * (_S575 * *(coeffs_18 + int(2)) + z_11 * *(coeffs_18 + int(5)) - x_14 * *(coeffs_18 + int(8)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh1_to_color_vjp_inplace(float3  mean_19, Matrix<float, 3, 3>  R_19, float3  t_19, float3  coeff_dc_19, float * coeffs_19, float3  v_colors_12, float3  * v_coeff_dc_12, float * v_coeffs_12, float3  * v_mean_12, Matrix<float, 3, 3>  * v_R_12, float3  * v_t_12)
{
    Matrix<float, 3, 3>  _S576 = transpose_0(R_19);
    float3  _S577 = mean_19 + mul_0(_S576, t_19);
    float _S578 = _S577.x;
    float _S579 = _S577.y;
    float _S580 = _S577.z;
    float inv_norm_13 = (F32_rsqrt((_S578 * _S578 + _S579 * _S579 + _S580 * _S580)));
    float x_15 = _S578 * inv_norm_13;
    float y_13 = _S579 * inv_norm_13;
    float z_12 = _S580 * inv_norm_13;
    float _S581 = - y_13;
    float * _S582 = coeffs_19 + int(0);
    float * _S583 = coeffs_19 + int(3);
    float * _S584 = coeffs_19 + int(6);
    float * _S585 = coeffs_19 + int(1);
    float * _S586 = coeffs_19 + int(4);
    float * _S587 = coeffs_19 + int(7);
    float * _S588 = coeffs_19 + int(2);
    float * _S589 = coeffs_19 + int(5);
    float * _S590 = coeffs_19 + int(8);
    float3  _S591 = v_colors_12 * make_float3 (float((0.282094806432724f * coeff_dc_19.x + 0.48860251903533936f * (_S581 * *_S582 + z_12 * *_S583 - x_15 * *_S584)) >= -0.5f), float((0.282094806432724f * coeff_dc_19.y + 0.48860251903533936f * (_S581 * *_S585 + z_12 * *_S586 - x_15 * *_S587)) >= -0.5f), float((0.282094806432724f * coeff_dc_19.z + 0.48860251903533936f * (_S581 * *_S588 + z_12 * *_S589 - x_15 * *_S590)) >= -0.5f));
    float3  v_viewdir_20 = {};
    float _S592 = _S591.x;
    *&(v_coeff_dc_12->x) = *&(v_coeff_dc_12->x) + 0.282094806432724f * _S592;
    float * _S593 = v_coeffs_12 + int(0);
    float _S594 = -0.48860251903533936f * y_13;
    *_S593 = *_S593 + _S594 * _S592;
    float * _S595 = v_coeffs_12 + int(3);
    float _S596 = 0.48860251903533936f * z_12;
    *_S595 = *_S595 + _S596 * _S592;
    float * _S597 = v_coeffs_12 + int(6);
    float _S598 = -0.48860251903533936f * x_15;
    *_S597 = *_S597 + _S598 * _S592;
    float3  dir_n_8 = make_float3 (x_15, y_13, z_12);
    float3  v_dir_n_8 = make_float3 (-0.48860251903533936f * *_S584 * _S592, -0.48860251903533936f * *_S582 * _S592, 0.48860251903533936f * *_S583 * _S592);
    float3  v_viewdir_21 = v_viewdir_20 + (v_dir_n_8 - make_float3 (dot_0(v_dir_n_8, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    float _S599 = _S591.y;
    *&(v_coeff_dc_12->y) = *&(v_coeff_dc_12->y) + 0.282094806432724f * _S599;
    float * _S600 = v_coeffs_12 + int(1);
    *_S600 = *_S600 + _S594 * _S599;
    float * _S601 = v_coeffs_12 + int(4);
    *_S601 = *_S601 + _S596 * _S599;
    float * _S602 = v_coeffs_12 + int(7);
    *_S602 = *_S602 + _S598 * _S599;
    float3  v_dir_n_9 = make_float3 (-0.48860251903533936f * *_S587 * _S599, -0.48860251903533936f * *_S585 * _S599, 0.48860251903533936f * *_S586 * _S599);
    float3  v_viewdir_22 = v_viewdir_21 + (v_dir_n_9 - make_float3 (dot_0(v_dir_n_9, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    float _S603 = _S591.z;
    *&(v_coeff_dc_12->z) = *&(v_coeff_dc_12->z) + 0.282094806432724f * _S603;
    float * _S604 = v_coeffs_12 + int(2);
    *_S604 = *_S604 + _S594 * _S603;
    float * _S605 = v_coeffs_12 + int(5);
    *_S605 = *_S605 + _S596 * _S603;
    float * _S606 = v_coeffs_12 + int(8);
    *_S606 = *_S606 + _S598 * _S603;
    float3  v_dir_n_10 = make_float3 (-0.48860251903533936f * *_S590 * _S603, -0.48860251903533936f * *_S588 * _S603, 0.48860251903533936f * *_S589 * _S603);
    float3  v_viewdir_23 = v_viewdir_22 + (v_dir_n_10 - make_float3 (dot_0(v_dir_n_10, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    Matrix<float, 3, 3>  _S607 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S608;
    (&_S608)->primal_0 = _S576;
    (&_S608)->differential_0 = _S607;
    float3  _S609 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S610;
    (&_S610)->primal_0 = t_19;
    (&_S610)->differential_0 = _S609;
    s_bwd_prop_mul_0(&_S608, &_S610, v_viewdir_23);
    Matrix<float, 3, 3>  _S611 = transpose_0(_S608.differential_0);
    *v_mean_12 = *v_mean_12 + v_viewdir_23;
    *v_R_12 = *v_R_12 + _S611;
    *v_t_12 = *v_t_12 + _S610.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_atomic(float3  mean_20, Matrix<float, 3, 3>  R_20, float3  t_20, float3  coeff_dc_20, float * coeffs_20, float3  v_colors_13, float3  * v_coeff_dc_13, float * v_coeffs_13, float3  * v_mean_13, Matrix<float, 3, 3>  * v_R_13, float3  * v_t_13)
{
    Matrix<float, 3, 3>  _S612 = transpose_0(R_20);
    float3  _S613 = mean_20 + mul_0(_S612, t_20);
    float _S614 = _S613.x;
    float _S615 = _S613.y;
    float _S616 = _S613.z;
    float inv_norm_14 = (F32_rsqrt((_S614 * _S614 + _S615 * _S615 + _S616 * _S616)));
    float x_16 = _S614 * inv_norm_14;
    float y_14 = _S615 * inv_norm_14;
    float z_13 = _S616 * inv_norm_14;
    float _S617 = - y_14;
    float * _S618 = coeffs_20 + int(0);
    float * _S619 = coeffs_20 + int(3);
    float * _S620 = coeffs_20 + int(6);
    float * _S621 = coeffs_20 + int(1);
    float * _S622 = coeffs_20 + int(4);
    float * _S623 = coeffs_20 + int(7);
    float * _S624 = coeffs_20 + int(2);
    float * _S625 = coeffs_20 + int(5);
    float * _S626 = coeffs_20 + int(8);
    float3  _S627 = v_colors_13 * make_float3 (float((0.282094806432724f * coeff_dc_20.x + 0.48860251903533936f * (_S617 * *_S618 + z_13 * *_S619 - x_16 * *_S620)) >= -0.5f), float((0.282094806432724f * coeff_dc_20.y + 0.48860251903533936f * (_S617 * *_S621 + z_13 * *_S622 - x_16 * *_S623)) >= -0.5f), float((0.282094806432724f * coeff_dc_20.z + 0.48860251903533936f * (_S617 * *_S624 + z_13 * *_S625 - x_16 * *_S626)) >= -0.5f));
    float3  v_viewdir_24 = {};
    float _S628 = _S627.x;
    *&(v_coeff_dc_13->x) = *&(v_coeff_dc_13->x) + 0.282094806432724f * _S628;
    float _S629 = -0.48860251903533936f * y_14;
    float temp_50 = _S629 * _S628;
    bool _S630;
    if((F32_isfinite((temp_50))))
    {
        _S630 = temp_50 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S631 = atomicAdd(v_coeffs_13 + int(0), temp_50);
    }
    float _S632 = 0.48860251903533936f * z_13;
    float temp_51 = _S632 * _S628;
    if((F32_isfinite((temp_51))))
    {
        _S630 = temp_51 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S633 = atomicAdd(v_coeffs_13 + int(3), temp_51);
    }
    float _S634 = -0.48860251903533936f * x_16;
    float temp_52 = _S634 * _S628;
    if((F32_isfinite((temp_52))))
    {
        _S630 = temp_52 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S635 = atomicAdd(v_coeffs_13 + int(6), temp_52);
    }
    float3  dir_n_9 = make_float3 (x_16, y_14, z_13);
    float3  v_dir_n_11 = make_float3 (-0.48860251903533936f * *_S620 * _S628, -0.48860251903533936f * *_S618 * _S628, 0.48860251903533936f * *_S619 * _S628);
    float3  v_viewdir_25 = v_viewdir_24 + (v_dir_n_11 - make_float3 (dot_0(v_dir_n_11, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    float _S636 = _S627.y;
    *&(v_coeff_dc_13->y) = *&(v_coeff_dc_13->y) + 0.282094806432724f * _S636;
    float temp_53 = _S629 * _S636;
    if((F32_isfinite((temp_53))))
    {
        _S630 = temp_53 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S637 = atomicAdd(v_coeffs_13 + int(1), temp_53);
    }
    float temp_54 = _S632 * _S636;
    if((F32_isfinite((temp_54))))
    {
        _S630 = temp_54 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S638 = atomicAdd(v_coeffs_13 + int(4), temp_54);
    }
    float temp_55 = _S634 * _S636;
    if((F32_isfinite((temp_55))))
    {
        _S630 = temp_55 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S639 = atomicAdd(v_coeffs_13 + int(7), temp_55);
    }
    float3  v_dir_n_12 = make_float3 (-0.48860251903533936f * *_S623 * _S636, -0.48860251903533936f * *_S621 * _S636, 0.48860251903533936f * *_S622 * _S636);
    float3  v_viewdir_26 = v_viewdir_25 + (v_dir_n_12 - make_float3 (dot_0(v_dir_n_12, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    float _S640 = _S627.z;
    *&(v_coeff_dc_13->z) = *&(v_coeff_dc_13->z) + 0.282094806432724f * _S640;
    float temp_56 = _S629 * _S640;
    if((F32_isfinite((temp_56))))
    {
        _S630 = temp_56 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S641 = atomicAdd(v_coeffs_13 + int(2), temp_56);
    }
    float temp_57 = _S632 * _S640;
    if((F32_isfinite((temp_57))))
    {
        _S630 = temp_57 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S642 = atomicAdd(v_coeffs_13 + int(5), temp_57);
    }
    float temp_58 = _S634 * _S640;
    if((F32_isfinite((temp_58))))
    {
        _S630 = temp_58 != 0.0f;
    }
    else
    {
        _S630 = false;
    }
    if(_S630)
    {
        float _S643 = atomicAdd(v_coeffs_13 + int(8), temp_58);
    }
    float3  v_dir_n_13 = make_float3 (-0.48860251903533936f * *_S626 * _S640, -0.48860251903533936f * *_S624 * _S640, 0.48860251903533936f * *_S625 * _S640);
    float3  v_viewdir_27 = v_viewdir_26 + (v_dir_n_13 - make_float3 (dot_0(v_dir_n_13, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    Matrix<float, 3, 3>  _S644 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S645;
    (&_S645)->primal_0 = _S612;
    (&_S645)->differential_0 = _S644;
    float3  _S646 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S647;
    (&_S647)->primal_0 = t_20;
    (&_S647)->differential_0 = _S646;
    s_bwd_prop_mul_0(&_S645, &_S647, v_viewdir_27);
    Matrix<float, 3, 3>  _S648 = transpose_0(_S645.differential_0);
    *v_mean_13 = *v_mean_13 + v_viewdir_27;
    *v_R_13 = *v_R_13 + _S648;
    *v_t_13 = *v_t_13 + _S647.differential_0;
    return;
}

inline __device__ float3  sh2_to_color(float3  mean_21, Matrix<float, 3, 3>  R_21, float3  t_21, float3  coeff_dc_21, float * coeffs_21)
{
    float3  _S649 = mean_21 + mul_0(transpose_0(R_21), t_21);
    float _S650 = _S649.x;
    float _S651 = _S649.y;
    float _S652 = _S649.z;
    float inv_norm_15 = (F32_rsqrt((_S650 * _S650 + _S651 * _S651 + _S652 * _S652)));
    float x_17 = _S650 * inv_norm_15;
    float y_15 = _S651 * inv_norm_15;
    float z_14 = _S652 * inv_norm_15;
    float _S653 = - y_15;
    float fTmp0B_9 = -1.09254848957061768f * z_14;
    float pSH6_7 = 0.94617468118667603f * (z_14 * z_14) - 0.31539157032966614f;
    float pSH7_6 = fTmp0B_9 * x_17;
    float pSH5_6 = fTmp0B_9 * y_15;
    float pSH8_6 = 0.54627424478530884f * (x_17 * x_17 - y_15 * y_15);
    float pSH4_6 = 0.54627424478530884f * (2.0f * x_17 * y_15);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_21.x + 0.48860251903533936f * (_S653 * *(coeffs_21 + int(0)) + z_14 * *(coeffs_21 + int(3)) - x_17 * *(coeffs_21 + int(6))) + (pSH4_6 * *(coeffs_21 + int(9)) + pSH5_6 * *(coeffs_21 + int(12)) + pSH6_7 * *(coeffs_21 + int(15)) + pSH7_6 * *(coeffs_21 + int(18)) + pSH8_6 * *(coeffs_21 + int(21))), 0.282094806432724f * coeff_dc_21.y + 0.48860251903533936f * (_S653 * *(coeffs_21 + int(1)) + z_14 * *(coeffs_21 + int(4)) - x_17 * *(coeffs_21 + int(7))) + (pSH4_6 * *(coeffs_21 + int(10)) + pSH5_6 * *(coeffs_21 + int(13)) + pSH6_7 * *(coeffs_21 + int(16)) + pSH7_6 * *(coeffs_21 + int(19)) + pSH8_6 * *(coeffs_21 + int(22))), 0.282094806432724f * coeff_dc_21.z + 0.48860251903533936f * (_S653 * *(coeffs_21 + int(2)) + z_14 * *(coeffs_21 + int(5)) - x_17 * *(coeffs_21 + int(8))) + (pSH4_6 * *(coeffs_21 + int(11)) + pSH5_6 * *(coeffs_21 + int(14)) + pSH6_7 * *(coeffs_21 + int(17)) + pSH7_6 * *(coeffs_21 + int(20)) + pSH8_6 * *(coeffs_21 + int(23)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh2_to_color_vjp_inplace(float3  mean_22, Matrix<float, 3, 3>  R_22, float3  t_22, float3  coeff_dc_22, float * coeffs_22, float3  v_colors_14, float3  * v_coeff_dc_14, float * v_coeffs_14, float3  * v_mean_14, Matrix<float, 3, 3>  * v_R_14, float3  * v_t_14)
{
    Matrix<float, 3, 3>  _S654 = transpose_0(R_22);
    float3  _S655 = mean_22 + mul_0(_S654, t_22);
    float _S656 = _S655.x;
    float _S657 = _S655.y;
    float _S658 = _S655.z;
    float inv_norm_16 = (F32_rsqrt((_S656 * _S656 + _S657 * _S657 + _S658 * _S658)));
    float x_18 = _S656 * inv_norm_16;
    float y_16 = _S657 * inv_norm_16;
    float z_15 = _S658 * inv_norm_16;
    float _S659 = - y_16;
    float * _S660 = coeffs_22 + int(0);
    float * _S661 = coeffs_22 + int(3);
    float * _S662 = coeffs_22 + int(6);
    float fTmp0B_10 = -1.09254848957061768f * z_15;
    float _S663 = 2.0f * x_18;
    float pSH6_8 = 0.94617468118667603f * (z_15 * z_15) - 0.31539157032966614f;
    float pSH7_7 = fTmp0B_10 * x_18;
    float pSH5_7 = fTmp0B_10 * y_16;
    float pSH8_7 = 0.54627424478530884f * (x_18 * x_18 - y_16 * y_16);
    float pSH4_7 = 0.54627424478530884f * (_S663 * y_16);
    float * _S664 = coeffs_22 + int(9);
    float * _S665 = coeffs_22 + int(12);
    float * _S666 = coeffs_22 + int(15);
    float * _S667 = coeffs_22 + int(18);
    float * _S668 = coeffs_22 + int(21);
    float * _S669 = coeffs_22 + int(1);
    float * _S670 = coeffs_22 + int(4);
    float * _S671 = coeffs_22 + int(7);
    float * _S672 = coeffs_22 + int(10);
    float * _S673 = coeffs_22 + int(13);
    float * _S674 = coeffs_22 + int(16);
    float * _S675 = coeffs_22 + int(19);
    float * _S676 = coeffs_22 + int(22);
    float * _S677 = coeffs_22 + int(2);
    float * _S678 = coeffs_22 + int(5);
    float * _S679 = coeffs_22 + int(8);
    float * _S680 = coeffs_22 + int(11);
    float * _S681 = coeffs_22 + int(14);
    float * _S682 = coeffs_22 + int(17);
    float * _S683 = coeffs_22 + int(20);
    float * _S684 = coeffs_22 + int(23);
    float3  _S685 = v_colors_14 * make_float3 (float((0.282094806432724f * coeff_dc_22.x + 0.48860251903533936f * (_S659 * *_S660 + z_15 * *_S661 - x_18 * *_S662) + (pSH4_7 * *_S664 + pSH5_7 * *_S665 + pSH6_8 * *_S666 + pSH7_7 * *_S667 + pSH8_7 * *_S668)) >= -0.5f), float((0.282094806432724f * coeff_dc_22.y + 0.48860251903533936f * (_S659 * *_S669 + z_15 * *_S670 - x_18 * *_S671) + (pSH4_7 * *_S672 + pSH5_7 * *_S673 + pSH6_8 * *_S674 + pSH7_7 * *_S675 + pSH8_7 * *_S676)) >= -0.5f), float((0.282094806432724f * coeff_dc_22.z + 0.48860251903533936f * (_S659 * *_S677 + z_15 * *_S678 - x_18 * *_S679) + (pSH4_7 * *_S680 + pSH5_7 * *_S681 + pSH6_8 * *_S682 + pSH7_7 * *_S683 + pSH8_7 * *_S684)) >= -0.5f));
    float3  v_viewdir_28 = {};
    float _S686 = _S685.x;
    *&(v_coeff_dc_14->x) = *&(v_coeff_dc_14->x) + 0.282094806432724f * _S686;
    float * _S687 = v_coeffs_14 + int(0);
    float _S688 = -0.48860251903533936f * y_16;
    *_S687 = *_S687 + _S688 * _S686;
    float * _S689 = v_coeffs_14 + int(3);
    float _S690 = 0.48860251903533936f * z_15;
    *_S689 = *_S689 + _S690 * _S686;
    float * _S691 = v_coeffs_14 + int(6);
    float _S692 = -0.48860251903533936f * x_18;
    *_S691 = *_S691 + _S692 * _S686;
    float _S693 = -0.48860251903533936f * *_S662 * _S686;
    float _S694 = -0.48860251903533936f * *_S660 * _S686;
    float _S695 = 0.48860251903533936f * *_S661 * _S686;
    float * _S696 = v_coeffs_14 + int(9);
    *_S696 = *_S696 + pSH4_7 * _S686;
    float * _S697 = v_coeffs_14 + int(12);
    *_S697 = *_S697 + pSH5_7 * _S686;
    float * _S698 = v_coeffs_14 + int(15);
    *_S698 = *_S698 + pSH6_8 * _S686;
    float * _S699 = v_coeffs_14 + int(18);
    *_S699 = *_S699 + pSH7_7 * _S686;
    float * _S700 = v_coeffs_14 + int(21);
    *_S700 = *_S700 + pSH8_7 * _S686;
    float pSH6_z_2 = 1.89234936237335205f * z_15;
    float pSH7_z_0 = -1.09254848957061768f * x_18;
    float pSH5_z_0 = -1.09254848957061768f * y_16;
    float pSH8_x_6 = 0.54627424478530884f * _S663;
    float pSH8_y_0 = 0.54627424478530884f * (-2.0f * y_16);
    float pSH4_x_0 = 0.54627424478530884f * (2.0f * y_16);
    float3  dir_n_10 = make_float3 (x_18, y_16, z_15);
    float3  v_dir_n_14 = make_float3 (_S693 + _S686 * (pSH4_x_0 * *_S664 + pSH8_x_6 * *_S668 + fTmp0B_10 * *_S667), _S694 + _S686 * (pSH8_x_6 * *_S664 + pSH8_y_0 * *_S668 + fTmp0B_10 * *_S665), _S695 + _S686 * (pSH6_z_2 * *_S666 + pSH7_z_0 * *_S667 + pSH5_z_0 * *_S665));
    float3  v_viewdir_29 = v_viewdir_28 + (v_dir_n_14 - make_float3 (dot_0(v_dir_n_14, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    float _S701 = _S685.y;
    *&(v_coeff_dc_14->y) = *&(v_coeff_dc_14->y) + 0.282094806432724f * _S701;
    float * _S702 = v_coeffs_14 + int(1);
    *_S702 = *_S702 + _S688 * _S701;
    float * _S703 = v_coeffs_14 + int(4);
    *_S703 = *_S703 + _S690 * _S701;
    float * _S704 = v_coeffs_14 + int(7);
    *_S704 = *_S704 + _S692 * _S701;
    float _S705 = -0.48860251903533936f * *_S671 * _S701;
    float _S706 = -0.48860251903533936f * *_S669 * _S701;
    float _S707 = 0.48860251903533936f * *_S670 * _S701;
    float * _S708 = v_coeffs_14 + int(10);
    *_S708 = *_S708 + pSH4_7 * _S701;
    float * _S709 = v_coeffs_14 + int(13);
    *_S709 = *_S709 + pSH5_7 * _S701;
    float * _S710 = v_coeffs_14 + int(16);
    *_S710 = *_S710 + pSH6_8 * _S701;
    float * _S711 = v_coeffs_14 + int(19);
    *_S711 = *_S711 + pSH7_7 * _S701;
    float * _S712 = v_coeffs_14 + int(22);
    *_S712 = *_S712 + pSH8_7 * _S701;
    float3  v_dir_n_15 = make_float3 (_S705 + _S701 * (pSH4_x_0 * *_S672 + pSH8_x_6 * *_S676 + fTmp0B_10 * *_S675), _S706 + _S701 * (pSH8_x_6 * *_S672 + pSH8_y_0 * *_S676 + fTmp0B_10 * *_S673), _S707 + _S701 * (pSH6_z_2 * *_S674 + pSH7_z_0 * *_S675 + pSH5_z_0 * *_S673));
    float3  v_viewdir_30 = v_viewdir_29 + (v_dir_n_15 - make_float3 (dot_0(v_dir_n_15, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    float _S713 = _S685.z;
    *&(v_coeff_dc_14->z) = *&(v_coeff_dc_14->z) + 0.282094806432724f * _S713;
    float * _S714 = v_coeffs_14 + int(2);
    *_S714 = *_S714 + _S688 * _S713;
    float * _S715 = v_coeffs_14 + int(5);
    *_S715 = *_S715 + _S690 * _S713;
    float * _S716 = v_coeffs_14 + int(8);
    *_S716 = *_S716 + _S692 * _S713;
    float _S717 = -0.48860251903533936f * *_S679 * _S713;
    float _S718 = -0.48860251903533936f * *_S677 * _S713;
    float _S719 = 0.48860251903533936f * *_S678 * _S713;
    float * _S720 = v_coeffs_14 + int(11);
    *_S720 = *_S720 + pSH4_7 * _S713;
    float * _S721 = v_coeffs_14 + int(14);
    *_S721 = *_S721 + pSH5_7 * _S713;
    float * _S722 = v_coeffs_14 + int(17);
    *_S722 = *_S722 + pSH6_8 * _S713;
    float * _S723 = v_coeffs_14 + int(20);
    *_S723 = *_S723 + pSH7_7 * _S713;
    float * _S724 = v_coeffs_14 + int(23);
    *_S724 = *_S724 + pSH8_7 * _S713;
    float3  v_dir_n_16 = make_float3 (_S717 + _S713 * (pSH4_x_0 * *_S680 + pSH8_x_6 * *_S684 + fTmp0B_10 * *_S683), _S718 + _S713 * (pSH8_x_6 * *_S680 + pSH8_y_0 * *_S684 + fTmp0B_10 * *_S681), _S719 + _S713 * (pSH6_z_2 * *_S682 + pSH7_z_0 * *_S683 + pSH5_z_0 * *_S681));
    float3  v_viewdir_31 = v_viewdir_30 + (v_dir_n_16 - make_float3 (dot_0(v_dir_n_16, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    Matrix<float, 3, 3>  _S725 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S726;
    (&_S726)->primal_0 = _S654;
    (&_S726)->differential_0 = _S725;
    float3  _S727 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S728;
    (&_S728)->primal_0 = t_22;
    (&_S728)->differential_0 = _S727;
    s_bwd_prop_mul_0(&_S726, &_S728, v_viewdir_31);
    Matrix<float, 3, 3>  _S729 = transpose_0(_S726.differential_0);
    *v_mean_14 = *v_mean_14 + v_viewdir_31;
    *v_R_14 = *v_R_14 + _S729;
    *v_t_14 = *v_t_14 + _S728.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_atomic(float3  mean_23, Matrix<float, 3, 3>  R_23, float3  t_23, float3  coeff_dc_23, float * coeffs_23, float3  v_colors_15, float3  * v_coeff_dc_15, float * v_coeffs_15, float3  * v_mean_15, Matrix<float, 3, 3>  * v_R_15, float3  * v_t_15)
{
    Matrix<float, 3, 3>  _S730 = transpose_0(R_23);
    float3  _S731 = mean_23 + mul_0(_S730, t_23);
    float _S732 = _S731.x;
    float _S733 = _S731.y;
    float _S734 = _S731.z;
    float inv_norm_17 = (F32_rsqrt((_S732 * _S732 + _S733 * _S733 + _S734 * _S734)));
    float x_19 = _S732 * inv_norm_17;
    float y_17 = _S733 * inv_norm_17;
    float z_16 = _S734 * inv_norm_17;
    float _S735 = - y_17;
    float * _S736 = coeffs_23 + int(0);
    float * _S737 = coeffs_23 + int(3);
    float * _S738 = coeffs_23 + int(6);
    float fTmp0B_11 = -1.09254848957061768f * z_16;
    float _S739 = 2.0f * x_19;
    float pSH6_9 = 0.94617468118667603f * (z_16 * z_16) - 0.31539157032966614f;
    float pSH7_8 = fTmp0B_11 * x_19;
    float pSH5_8 = fTmp0B_11 * y_17;
    float pSH8_8 = 0.54627424478530884f * (x_19 * x_19 - y_17 * y_17);
    float pSH4_8 = 0.54627424478530884f * (_S739 * y_17);
    float * _S740 = coeffs_23 + int(9);
    float * _S741 = coeffs_23 + int(12);
    float * _S742 = coeffs_23 + int(15);
    float * _S743 = coeffs_23 + int(18);
    float * _S744 = coeffs_23 + int(21);
    float * _S745 = coeffs_23 + int(1);
    float * _S746 = coeffs_23 + int(4);
    float * _S747 = coeffs_23 + int(7);
    float * _S748 = coeffs_23 + int(10);
    float * _S749 = coeffs_23 + int(13);
    float * _S750 = coeffs_23 + int(16);
    float * _S751 = coeffs_23 + int(19);
    float * _S752 = coeffs_23 + int(22);
    float * _S753 = coeffs_23 + int(2);
    float * _S754 = coeffs_23 + int(5);
    float * _S755 = coeffs_23 + int(8);
    float * _S756 = coeffs_23 + int(11);
    float * _S757 = coeffs_23 + int(14);
    float * _S758 = coeffs_23 + int(17);
    float * _S759 = coeffs_23 + int(20);
    float * _S760 = coeffs_23 + int(23);
    float3  _S761 = v_colors_15 * make_float3 (float((0.282094806432724f * coeff_dc_23.x + 0.48860251903533936f * (_S735 * *_S736 + z_16 * *_S737 - x_19 * *_S738) + (pSH4_8 * *_S740 + pSH5_8 * *_S741 + pSH6_9 * *_S742 + pSH7_8 * *_S743 + pSH8_8 * *_S744)) >= -0.5f), float((0.282094806432724f * coeff_dc_23.y + 0.48860251903533936f * (_S735 * *_S745 + z_16 * *_S746 - x_19 * *_S747) + (pSH4_8 * *_S748 + pSH5_8 * *_S749 + pSH6_9 * *_S750 + pSH7_8 * *_S751 + pSH8_8 * *_S752)) >= -0.5f), float((0.282094806432724f * coeff_dc_23.z + 0.48860251903533936f * (_S735 * *_S753 + z_16 * *_S754 - x_19 * *_S755) + (pSH4_8 * *_S756 + pSH5_8 * *_S757 + pSH6_9 * *_S758 + pSH7_8 * *_S759 + pSH8_8 * *_S760)) >= -0.5f));
    float3  v_viewdir_32 = {};
    float _S762 = _S761.x;
    *&(v_coeff_dc_15->x) = *&(v_coeff_dc_15->x) + 0.282094806432724f * _S762;
    float _S763 = -0.48860251903533936f * y_17;
    float temp_59 = _S763 * _S762;
    bool _S764;
    if((F32_isfinite((temp_59))))
    {
        _S764 = temp_59 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S765 = atomicAdd(v_coeffs_15 + int(0), temp_59);
    }
    float _S766 = 0.48860251903533936f * z_16;
    float temp_60 = _S766 * _S762;
    if((F32_isfinite((temp_60))))
    {
        _S764 = temp_60 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S767 = atomicAdd(v_coeffs_15 + int(3), temp_60);
    }
    float _S768 = -0.48860251903533936f * x_19;
    float temp_61 = _S768 * _S762;
    if((F32_isfinite((temp_61))))
    {
        _S764 = temp_61 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S769 = atomicAdd(v_coeffs_15 + int(6), temp_61);
    }
    float _S770 = -0.48860251903533936f * *_S738 * _S762;
    float _S771 = -0.48860251903533936f * *_S736 * _S762;
    float _S772 = 0.48860251903533936f * *_S737 * _S762;
    float temp_62 = pSH4_8 * _S762;
    if((F32_isfinite((temp_62))))
    {
        _S764 = temp_62 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S773 = atomicAdd(v_coeffs_15 + int(9), temp_62);
    }
    float temp_63 = pSH5_8 * _S762;
    if((F32_isfinite((temp_63))))
    {
        _S764 = temp_63 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S774 = atomicAdd(v_coeffs_15 + int(12), temp_63);
    }
    float temp_64 = pSH6_9 * _S762;
    if((F32_isfinite((temp_64))))
    {
        _S764 = temp_64 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S775 = atomicAdd(v_coeffs_15 + int(15), temp_64);
    }
    float temp_65 = pSH7_8 * _S762;
    if((F32_isfinite((temp_65))))
    {
        _S764 = temp_65 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S776 = atomicAdd(v_coeffs_15 + int(18), temp_65);
    }
    float temp_66 = pSH8_8 * _S762;
    if((F32_isfinite((temp_66))))
    {
        _S764 = temp_66 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S777 = atomicAdd(v_coeffs_15 + int(21), temp_66);
    }
    float pSH6_z_3 = 1.89234936237335205f * z_16;
    float pSH7_z_1 = -1.09254848957061768f * x_19;
    float pSH5_z_1 = -1.09254848957061768f * y_17;
    float pSH8_x_7 = 0.54627424478530884f * _S739;
    float pSH8_y_1 = 0.54627424478530884f * (-2.0f * y_17);
    float pSH4_x_1 = 0.54627424478530884f * (2.0f * y_17);
    float3  dir_n_11 = make_float3 (x_19, y_17, z_16);
    float3  v_dir_n_17 = make_float3 (_S770 + _S762 * (pSH4_x_1 * *_S740 + pSH8_x_7 * *_S744 + fTmp0B_11 * *_S743), _S771 + _S762 * (pSH8_x_7 * *_S740 + pSH8_y_1 * *_S744 + fTmp0B_11 * *_S741), _S772 + _S762 * (pSH6_z_3 * *_S742 + pSH7_z_1 * *_S743 + pSH5_z_1 * *_S741));
    float3  v_viewdir_33 = v_viewdir_32 + (v_dir_n_17 - make_float3 (dot_0(v_dir_n_17, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    float _S778 = _S761.y;
    *&(v_coeff_dc_15->y) = *&(v_coeff_dc_15->y) + 0.282094806432724f * _S778;
    float temp_67 = _S763 * _S778;
    if((F32_isfinite((temp_67))))
    {
        _S764 = temp_67 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S779 = atomicAdd(v_coeffs_15 + int(1), temp_67);
    }
    float temp_68 = _S766 * _S778;
    if((F32_isfinite((temp_68))))
    {
        _S764 = temp_68 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S780 = atomicAdd(v_coeffs_15 + int(4), temp_68);
    }
    float temp_69 = _S768 * _S778;
    if((F32_isfinite((temp_69))))
    {
        _S764 = temp_69 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S781 = atomicAdd(v_coeffs_15 + int(7), temp_69);
    }
    float _S782 = -0.48860251903533936f * *_S747 * _S778;
    float _S783 = -0.48860251903533936f * *_S745 * _S778;
    float _S784 = 0.48860251903533936f * *_S746 * _S778;
    float temp_70 = pSH4_8 * _S778;
    if((F32_isfinite((temp_70))))
    {
        _S764 = temp_70 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S785 = atomicAdd(v_coeffs_15 + int(10), temp_70);
    }
    float temp_71 = pSH5_8 * _S778;
    if((F32_isfinite((temp_71))))
    {
        _S764 = temp_71 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S786 = atomicAdd(v_coeffs_15 + int(13), temp_71);
    }
    float temp_72 = pSH6_9 * _S778;
    if((F32_isfinite((temp_72))))
    {
        _S764 = temp_72 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S787 = atomicAdd(v_coeffs_15 + int(16), temp_72);
    }
    float temp_73 = pSH7_8 * _S778;
    if((F32_isfinite((temp_73))))
    {
        _S764 = temp_73 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S788 = atomicAdd(v_coeffs_15 + int(19), temp_73);
    }
    float temp_74 = pSH8_8 * _S778;
    if((F32_isfinite((temp_74))))
    {
        _S764 = temp_74 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S789 = atomicAdd(v_coeffs_15 + int(22), temp_74);
    }
    float3  v_dir_n_18 = make_float3 (_S782 + _S778 * (pSH4_x_1 * *_S748 + pSH8_x_7 * *_S752 + fTmp0B_11 * *_S751), _S783 + _S778 * (pSH8_x_7 * *_S748 + pSH8_y_1 * *_S752 + fTmp0B_11 * *_S749), _S784 + _S778 * (pSH6_z_3 * *_S750 + pSH7_z_1 * *_S751 + pSH5_z_1 * *_S749));
    float3  v_viewdir_34 = v_viewdir_33 + (v_dir_n_18 - make_float3 (dot_0(v_dir_n_18, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    float _S790 = _S761.z;
    *&(v_coeff_dc_15->z) = *&(v_coeff_dc_15->z) + 0.282094806432724f * _S790;
    float temp_75 = _S763 * _S790;
    if((F32_isfinite((temp_75))))
    {
        _S764 = temp_75 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S791 = atomicAdd(v_coeffs_15 + int(2), temp_75);
    }
    float temp_76 = _S766 * _S790;
    if((F32_isfinite((temp_76))))
    {
        _S764 = temp_76 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S792 = atomicAdd(v_coeffs_15 + int(5), temp_76);
    }
    float temp_77 = _S768 * _S790;
    if((F32_isfinite((temp_77))))
    {
        _S764 = temp_77 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S793 = atomicAdd(v_coeffs_15 + int(8), temp_77);
    }
    float _S794 = -0.48860251903533936f * *_S755 * _S790;
    float _S795 = -0.48860251903533936f * *_S753 * _S790;
    float _S796 = 0.48860251903533936f * *_S754 * _S790;
    float temp_78 = pSH4_8 * _S790;
    if((F32_isfinite((temp_78))))
    {
        _S764 = temp_78 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S797 = atomicAdd(v_coeffs_15 + int(11), temp_78);
    }
    float temp_79 = pSH5_8 * _S790;
    if((F32_isfinite((temp_79))))
    {
        _S764 = temp_79 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S798 = atomicAdd(v_coeffs_15 + int(14), temp_79);
    }
    float temp_80 = pSH6_9 * _S790;
    if((F32_isfinite((temp_80))))
    {
        _S764 = temp_80 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S799 = atomicAdd(v_coeffs_15 + int(17), temp_80);
    }
    float temp_81 = pSH7_8 * _S790;
    if((F32_isfinite((temp_81))))
    {
        _S764 = temp_81 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S800 = atomicAdd(v_coeffs_15 + int(20), temp_81);
    }
    float temp_82 = pSH8_8 * _S790;
    if((F32_isfinite((temp_82))))
    {
        _S764 = temp_82 != 0.0f;
    }
    else
    {
        _S764 = false;
    }
    if(_S764)
    {
        float _S801 = atomicAdd(v_coeffs_15 + int(23), temp_82);
    }
    float3  v_dir_n_19 = make_float3 (_S794 + _S790 * (pSH4_x_1 * *_S756 + pSH8_x_7 * *_S760 + fTmp0B_11 * *_S759), _S795 + _S790 * (pSH8_x_7 * *_S756 + pSH8_y_1 * *_S760 + fTmp0B_11 * *_S757), _S796 + _S790 * (pSH6_z_3 * *_S758 + pSH7_z_1 * *_S759 + pSH5_z_1 * *_S757));
    float3  v_viewdir_35 = v_viewdir_34 + (v_dir_n_19 - make_float3 (dot_0(v_dir_n_19, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    Matrix<float, 3, 3>  _S802 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S803;
    (&_S803)->primal_0 = _S730;
    (&_S803)->differential_0 = _S802;
    float3  _S804 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S805;
    (&_S805)->primal_0 = t_23;
    (&_S805)->differential_0 = _S804;
    s_bwd_prop_mul_0(&_S803, &_S805, v_viewdir_35);
    Matrix<float, 3, 3>  _S806 = transpose_0(_S803.differential_0);
    *v_mean_15 = *v_mean_15 + v_viewdir_35;
    *v_R_15 = *v_R_15 + _S806;
    *v_t_15 = *v_t_15 + _S805.differential_0;
    return;
}

inline __device__ float3  sh3_to_color(float3  mean_24, Matrix<float, 3, 3>  R_24, float3  t_24, float3  coeff_dc_24, float * coeffs_24)
{
    float3  _S807 = mean_24 + mul_0(transpose_0(R_24), t_24);
    float _S808 = _S807.x;
    float _S809 = _S807.y;
    float _S810 = _S807.z;
    float inv_norm_18 = (F32_rsqrt((_S808 * _S808 + _S809 * _S809 + _S810 * _S810)));
    float x_20 = _S808 * inv_norm_18;
    float y_18 = _S809 * inv_norm_18;
    float z_17 = _S810 * inv_norm_18;
    float _S811 = - y_18;
    float z2_6 = z_17 * z_17;
    float fTmp0B_12 = -1.09254848957061768f * z_17;
    float fC1_6 = x_20 * x_20 - y_18 * y_18;
    float fS1_6 = 2.0f * x_20 * y_18;
    float pSH6_10 = 0.94617468118667603f * z2_6 - 0.31539157032966614f;
    float pSH7_9 = fTmp0B_12 * x_20;
    float pSH5_9 = fTmp0B_12 * y_18;
    float pSH8_9 = 0.54627424478530884f * fC1_6;
    float pSH4_9 = 0.54627424478530884f * fS1_6;
    float fTmp0C_6 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_17;
    float pSH12_5 = z_17 * (1.86588168144226074f * z2_6 - 1.11952900886535645f);
    float pSH13_4 = fTmp0C_6 * x_20;
    float pSH11_4 = fTmp0C_6 * y_18;
    float pSH14_4 = fTmp1B_6 * fC1_6;
    float pSH10_4 = fTmp1B_6 * fS1_6;
    float pSH15_4 = -0.59004360437393188f * (x_20 * fC1_6 - y_18 * fS1_6);
    float pSH9_4 = -0.59004360437393188f * (x_20 * fS1_6 + y_18 * fC1_6);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_24.x + 0.48860251903533936f * (_S811 * *(coeffs_24 + int(0)) + z_17 * *(coeffs_24 + int(3)) - x_20 * *(coeffs_24 + int(6))) + (pSH4_9 * *(coeffs_24 + int(9)) + pSH5_9 * *(coeffs_24 + int(12)) + pSH6_10 * *(coeffs_24 + int(15)) + pSH7_9 * *(coeffs_24 + int(18)) + pSH8_9 * *(coeffs_24 + int(21))) + (pSH9_4 * *(coeffs_24 + int(24)) + pSH10_4 * *(coeffs_24 + int(27)) + pSH11_4 * *(coeffs_24 + int(30)) + pSH12_5 * *(coeffs_24 + int(33)) + pSH13_4 * *(coeffs_24 + int(36)) + pSH14_4 * *(coeffs_24 + int(39)) + pSH15_4 * *(coeffs_24 + int(42))), 0.282094806432724f * coeff_dc_24.y + 0.48860251903533936f * (_S811 * *(coeffs_24 + int(1)) + z_17 * *(coeffs_24 + int(4)) - x_20 * *(coeffs_24 + int(7))) + (pSH4_9 * *(coeffs_24 + int(10)) + pSH5_9 * *(coeffs_24 + int(13)) + pSH6_10 * *(coeffs_24 + int(16)) + pSH7_9 * *(coeffs_24 + int(19)) + pSH8_9 * *(coeffs_24 + int(22))) + (pSH9_4 * *(coeffs_24 + int(25)) + pSH10_4 * *(coeffs_24 + int(28)) + pSH11_4 * *(coeffs_24 + int(31)) + pSH12_5 * *(coeffs_24 + int(34)) + pSH13_4 * *(coeffs_24 + int(37)) + pSH14_4 * *(coeffs_24 + int(40)) + pSH15_4 * *(coeffs_24 + int(43))), 0.282094806432724f * coeff_dc_24.z + 0.48860251903533936f * (_S811 * *(coeffs_24 + int(2)) + z_17 * *(coeffs_24 + int(5)) - x_20 * *(coeffs_24 + int(8))) + (pSH4_9 * *(coeffs_24 + int(11)) + pSH5_9 * *(coeffs_24 + int(14)) + pSH6_10 * *(coeffs_24 + int(17)) + pSH7_9 * *(coeffs_24 + int(20)) + pSH8_9 * *(coeffs_24 + int(23))) + (pSH9_4 * *(coeffs_24 + int(26)) + pSH10_4 * *(coeffs_24 + int(29)) + pSH11_4 * *(coeffs_24 + int(32)) + pSH12_5 * *(coeffs_24 + int(35)) + pSH13_4 * *(coeffs_24 + int(38)) + pSH14_4 * *(coeffs_24 + int(41)) + pSH15_4 * *(coeffs_24 + int(44)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh3_to_color_vjp_inplace(float3  mean_25, Matrix<float, 3, 3>  R_25, float3  t_25, float3  coeff_dc_25, float * coeffs_25, float3  v_colors_16, float3  * v_coeff_dc_16, float * v_coeffs_16, float3  * v_mean_16, Matrix<float, 3, 3>  * v_R_16, float3  * v_t_16)
{
    Matrix<float, 3, 3>  _S812 = transpose_0(R_25);
    float3  _S813 = mean_25 + mul_0(_S812, t_25);
    float _S814 = _S813.x;
    float _S815 = _S813.y;
    float _S816 = _S813.z;
    float inv_norm_19 = (F32_rsqrt((_S814 * _S814 + _S815 * _S815 + _S816 * _S816)));
    float x_21 = _S814 * inv_norm_19;
    float y_19 = _S815 * inv_norm_19;
    float z_18 = _S816 * inv_norm_19;
    float _S817 = - y_19;
    float * _S818 = coeffs_25 + int(0);
    float * _S819 = coeffs_25 + int(3);
    float * _S820 = coeffs_25 + int(6);
    float z2_7 = z_18 * z_18;
    float fTmp0B_13 = -1.09254848957061768f * z_18;
    float fC1_7 = x_21 * x_21 - y_19 * y_19;
    float _S821 = 2.0f * x_21;
    float fS1_7 = _S821 * y_19;
    float pSH6_11 = 0.94617468118667603f * z2_7 - 0.31539157032966614f;
    float pSH7_10 = fTmp0B_13 * x_21;
    float pSH5_10 = fTmp0B_13 * y_19;
    float pSH8_10 = 0.54627424478530884f * fC1_7;
    float pSH4_10 = 0.54627424478530884f * fS1_7;
    float * _S822 = coeffs_25 + int(9);
    float * _S823 = coeffs_25 + int(12);
    float * _S824 = coeffs_25 + int(15);
    float * _S825 = coeffs_25 + int(18);
    float * _S826 = coeffs_25 + int(21);
    float fTmp0C_7 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_18;
    float pSH12_6 = z_18 * (1.86588168144226074f * z2_7 - 1.11952900886535645f);
    float pSH13_5 = fTmp0C_7 * x_21;
    float pSH11_5 = fTmp0C_7 * y_19;
    float pSH14_5 = fTmp1B_7 * fC1_7;
    float pSH10_5 = fTmp1B_7 * fS1_7;
    float pSH15_5 = -0.59004360437393188f * (x_21 * fC1_7 - y_19 * fS1_7);
    float pSH9_5 = -0.59004360437393188f * (x_21 * fS1_7 + y_19 * fC1_7);
    float * _S827 = coeffs_25 + int(24);
    float * _S828 = coeffs_25 + int(27);
    float * _S829 = coeffs_25 + int(30);
    float * _S830 = coeffs_25 + int(33);
    float * _S831 = coeffs_25 + int(36);
    float * _S832 = coeffs_25 + int(39);
    float * _S833 = coeffs_25 + int(42);
    float * _S834 = coeffs_25 + int(1);
    float * _S835 = coeffs_25 + int(4);
    float * _S836 = coeffs_25 + int(7);
    float * _S837 = coeffs_25 + int(10);
    float * _S838 = coeffs_25 + int(13);
    float * _S839 = coeffs_25 + int(16);
    float * _S840 = coeffs_25 + int(19);
    float * _S841 = coeffs_25 + int(22);
    float * _S842 = coeffs_25 + int(25);
    float * _S843 = coeffs_25 + int(28);
    float * _S844 = coeffs_25 + int(31);
    float * _S845 = coeffs_25 + int(34);
    float * _S846 = coeffs_25 + int(37);
    float * _S847 = coeffs_25 + int(40);
    float * _S848 = coeffs_25 + int(43);
    float * _S849 = coeffs_25 + int(2);
    float * _S850 = coeffs_25 + int(5);
    float * _S851 = coeffs_25 + int(8);
    float * _S852 = coeffs_25 + int(11);
    float * _S853 = coeffs_25 + int(14);
    float * _S854 = coeffs_25 + int(17);
    float * _S855 = coeffs_25 + int(20);
    float * _S856 = coeffs_25 + int(23);
    float * _S857 = coeffs_25 + int(26);
    float * _S858 = coeffs_25 + int(29);
    float * _S859 = coeffs_25 + int(32);
    float * _S860 = coeffs_25 + int(35);
    float * _S861 = coeffs_25 + int(38);
    float * _S862 = coeffs_25 + int(41);
    float * _S863 = coeffs_25 + int(44);
    float3  _S864 = v_colors_16 * make_float3 (float((0.282094806432724f * coeff_dc_25.x + 0.48860251903533936f * (_S817 * *_S818 + z_18 * *_S819 - x_21 * *_S820) + (pSH4_10 * *_S822 + pSH5_10 * *_S823 + pSH6_11 * *_S824 + pSH7_10 * *_S825 + pSH8_10 * *_S826) + (pSH9_5 * *_S827 + pSH10_5 * *_S828 + pSH11_5 * *_S829 + pSH12_6 * *_S830 + pSH13_5 * *_S831 + pSH14_5 * *_S832 + pSH15_5 * *_S833)) >= -0.5f), float((0.282094806432724f * coeff_dc_25.y + 0.48860251903533936f * (_S817 * *_S834 + z_18 * *_S835 - x_21 * *_S836) + (pSH4_10 * *_S837 + pSH5_10 * *_S838 + pSH6_11 * *_S839 + pSH7_10 * *_S840 + pSH8_10 * *_S841) + (pSH9_5 * *_S842 + pSH10_5 * *_S843 + pSH11_5 * *_S844 + pSH12_6 * *_S845 + pSH13_5 * *_S846 + pSH14_5 * *_S847 + pSH15_5 * *_S848)) >= -0.5f), float((0.282094806432724f * coeff_dc_25.z + 0.48860251903533936f * (_S817 * *_S849 + z_18 * *_S850 - x_21 * *_S851) + (pSH4_10 * *_S852 + pSH5_10 * *_S853 + pSH6_11 * *_S854 + pSH7_10 * *_S855 + pSH8_10 * *_S856) + (pSH9_5 * *_S857 + pSH10_5 * *_S858 + pSH11_5 * *_S859 + pSH12_6 * *_S860 + pSH13_5 * *_S861 + pSH14_5 * *_S862 + pSH15_5 * *_S863)) >= -0.5f));
    float3  v_viewdir_36 = {};
    float _S865 = _S864.x;
    *&(v_coeff_dc_16->x) = *&(v_coeff_dc_16->x) + 0.282094806432724f * _S865;
    float * _S866 = v_coeffs_16 + int(0);
    float _S867 = -0.48860251903533936f * y_19;
    *_S866 = *_S866 + _S867 * _S865;
    float * _S868 = v_coeffs_16 + int(3);
    float _S869 = 0.48860251903533936f * z_18;
    *_S868 = *_S868 + _S869 * _S865;
    float * _S870 = v_coeffs_16 + int(6);
    float _S871 = -0.48860251903533936f * x_21;
    *_S870 = *_S870 + _S871 * _S865;
    float _S872 = -0.48860251903533936f * *_S820 * _S865;
    float _S873 = -0.48860251903533936f * *_S818 * _S865;
    float _S874 = 0.48860251903533936f * *_S819 * _S865;
    float * _S875 = v_coeffs_16 + int(9);
    *_S875 = *_S875 + pSH4_10 * _S865;
    float * _S876 = v_coeffs_16 + int(12);
    *_S876 = *_S876 + pSH5_10 * _S865;
    float * _S877 = v_coeffs_16 + int(15);
    *_S877 = *_S877 + pSH6_11 * _S865;
    float * _S878 = v_coeffs_16 + int(18);
    *_S878 = *_S878 + pSH7_10 * _S865;
    float * _S879 = v_coeffs_16 + int(21);
    *_S879 = *_S879 + pSH8_10 * _S865;
    float fC1_y_4 = -2.0f * y_19;
    float fS1_x_4 = 2.0f * y_19;
    float pSH6_z_4 = 1.89234936237335205f * z_18;
    float pSH7_z_2 = -1.09254848957061768f * x_21;
    float pSH5_z_2 = -1.09254848957061768f * y_19;
    float pSH8_x_8 = 0.54627424478530884f * _S821;
    float pSH8_y_2 = 0.54627424478530884f * fC1_y_4;
    float pSH4_x_2 = 0.54627424478530884f * fS1_x_4;
    float v_x_6 = _S872 + _S865 * (pSH4_x_2 * *_S822 + pSH8_x_8 * *_S826 + fTmp0B_13 * *_S825);
    float v_y_6 = _S873 + _S865 * (pSH8_x_8 * *_S822 + pSH8_y_2 * *_S826 + fTmp0B_13 * *_S823);
    float v_z_6 = _S874 + _S865 * (pSH6_z_4 * *_S824 + pSH7_z_2 * *_S825 + pSH5_z_2 * *_S823);
    float * _S880 = v_coeffs_16 + int(24);
    *_S880 = *_S880 + pSH9_5 * _S865;
    float * _S881 = v_coeffs_16 + int(27);
    *_S881 = *_S881 + pSH10_5 * _S865;
    float * _S882 = v_coeffs_16 + int(30);
    *_S882 = *_S882 + pSH11_5 * _S865;
    float * _S883 = v_coeffs_16 + int(33);
    *_S883 = *_S883 + pSH12_6 * _S865;
    float * _S884 = v_coeffs_16 + int(36);
    *_S884 = *_S884 + pSH13_5 * _S865;
    float * _S885 = v_coeffs_16 + int(39);
    *_S885 = *_S885 + pSH14_5 * _S865;
    float * _S886 = v_coeffs_16 + int(42);
    *_S886 = *_S886 + pSH15_5 * _S865;
    float fTmp0C_z_4 = -4.57045793533325195f * z_18;
    float _S887 = x_21 * _S821;
    float _S888 = y_19 * _S821;
    float pSH12_z_2 = 5.59764480590820312f * z2_7 - 1.11952900886535645f;
    float pSH13_z_0 = fTmp0C_z_4 * x_21;
    float pSH11_z_0 = fTmp0C_z_4 * y_19;
    float pSH14_x_4 = fTmp1B_7 * _S821;
    float pSH14_y_0 = fTmp1B_7 * fC1_y_4;
    float pSH14_z_0 = 1.44530570507049561f * fC1_7;
    float pSH10_x_0 = fTmp1B_7 * fS1_x_4;
    float pSH10_z_0 = 1.44530570507049561f * fS1_7;
    float pSH15_x_0 = -0.59004360437393188f * (fC1_7 + _S887 - y_19 * fS1_x_4);
    float pSH15_y_0 = -0.59004360437393188f * (x_21 * fC1_y_4 - fS1_7 - _S888);
    float pSH9_x_0 = -0.59004360437393188f * (fS1_7 + x_21 * fS1_x_4 + _S888);
    float pSH9_y_0 = -0.59004360437393188f * (_S887 + fC1_7 + y_19 * fC1_y_4);
    float3  dir_n_12 = make_float3 (x_21, y_19, z_18);
    float3  v_dir_n_20 = make_float3 (v_x_6 + _S865 * (pSH9_x_0 * *_S827 + pSH15_x_0 * *_S833 + pSH10_x_0 * *_S828 + pSH14_x_4 * *_S832 + fTmp0C_7 * *_S831), v_y_6 + _S865 * (pSH9_y_0 * *_S827 + pSH15_y_0 * *_S833 + pSH14_x_4 * *_S828 + pSH14_y_0 * *_S832 + fTmp0C_7 * *_S829), v_z_6 + _S865 * (pSH12_z_2 * *_S830 + pSH13_z_0 * *_S831 + pSH11_z_0 * *_S829 + pSH14_z_0 * *_S832 + pSH10_z_0 * *_S828));
    float3  v_viewdir_37 = v_viewdir_36 + (v_dir_n_20 - make_float3 (dot_0(v_dir_n_20, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    float _S889 = _S864.y;
    *&(v_coeff_dc_16->y) = *&(v_coeff_dc_16->y) + 0.282094806432724f * _S889;
    float * _S890 = v_coeffs_16 + int(1);
    *_S890 = *_S890 + _S867 * _S889;
    float * _S891 = v_coeffs_16 + int(4);
    *_S891 = *_S891 + _S869 * _S889;
    float * _S892 = v_coeffs_16 + int(7);
    *_S892 = *_S892 + _S871 * _S889;
    float _S893 = -0.48860251903533936f * *_S836 * _S889;
    float _S894 = -0.48860251903533936f * *_S834 * _S889;
    float _S895 = 0.48860251903533936f * *_S835 * _S889;
    float * _S896 = v_coeffs_16 + int(10);
    *_S896 = *_S896 + pSH4_10 * _S889;
    float * _S897 = v_coeffs_16 + int(13);
    *_S897 = *_S897 + pSH5_10 * _S889;
    float * _S898 = v_coeffs_16 + int(16);
    *_S898 = *_S898 + pSH6_11 * _S889;
    float * _S899 = v_coeffs_16 + int(19);
    *_S899 = *_S899 + pSH7_10 * _S889;
    float * _S900 = v_coeffs_16 + int(22);
    *_S900 = *_S900 + pSH8_10 * _S889;
    float v_x_7 = _S893 + _S889 * (pSH4_x_2 * *_S837 + pSH8_x_8 * *_S841 + fTmp0B_13 * *_S840);
    float v_y_7 = _S894 + _S889 * (pSH8_x_8 * *_S837 + pSH8_y_2 * *_S841 + fTmp0B_13 * *_S838);
    float v_z_7 = _S895 + _S889 * (pSH6_z_4 * *_S839 + pSH7_z_2 * *_S840 + pSH5_z_2 * *_S838);
    float * _S901 = v_coeffs_16 + int(25);
    *_S901 = *_S901 + pSH9_5 * _S889;
    float * _S902 = v_coeffs_16 + int(28);
    *_S902 = *_S902 + pSH10_5 * _S889;
    float * _S903 = v_coeffs_16 + int(31);
    *_S903 = *_S903 + pSH11_5 * _S889;
    float * _S904 = v_coeffs_16 + int(34);
    *_S904 = *_S904 + pSH12_6 * _S889;
    float * _S905 = v_coeffs_16 + int(37);
    *_S905 = *_S905 + pSH13_5 * _S889;
    float * _S906 = v_coeffs_16 + int(40);
    *_S906 = *_S906 + pSH14_5 * _S889;
    float * _S907 = v_coeffs_16 + int(43);
    *_S907 = *_S907 + pSH15_5 * _S889;
    float3  v_dir_n_21 = make_float3 (v_x_7 + _S889 * (pSH9_x_0 * *_S842 + pSH15_x_0 * *_S848 + pSH10_x_0 * *_S843 + pSH14_x_4 * *_S847 + fTmp0C_7 * *_S846), v_y_7 + _S889 * (pSH9_y_0 * *_S842 + pSH15_y_0 * *_S848 + pSH14_x_4 * *_S843 + pSH14_y_0 * *_S847 + fTmp0C_7 * *_S844), v_z_7 + _S889 * (pSH12_z_2 * *_S845 + pSH13_z_0 * *_S846 + pSH11_z_0 * *_S844 + pSH14_z_0 * *_S847 + pSH10_z_0 * *_S843));
    float3  v_viewdir_38 = v_viewdir_37 + (v_dir_n_21 - make_float3 (dot_0(v_dir_n_21, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    float _S908 = _S864.z;
    *&(v_coeff_dc_16->z) = *&(v_coeff_dc_16->z) + 0.282094806432724f * _S908;
    float * _S909 = v_coeffs_16 + int(2);
    *_S909 = *_S909 + _S867 * _S908;
    float * _S910 = v_coeffs_16 + int(5);
    *_S910 = *_S910 + _S869 * _S908;
    float * _S911 = v_coeffs_16 + int(8);
    *_S911 = *_S911 + _S871 * _S908;
    float _S912 = -0.48860251903533936f * *_S851 * _S908;
    float _S913 = -0.48860251903533936f * *_S849 * _S908;
    float _S914 = 0.48860251903533936f * *_S850 * _S908;
    float * _S915 = v_coeffs_16 + int(11);
    *_S915 = *_S915 + pSH4_10 * _S908;
    float * _S916 = v_coeffs_16 + int(14);
    *_S916 = *_S916 + pSH5_10 * _S908;
    float * _S917 = v_coeffs_16 + int(17);
    *_S917 = *_S917 + pSH6_11 * _S908;
    float * _S918 = v_coeffs_16 + int(20);
    *_S918 = *_S918 + pSH7_10 * _S908;
    float * _S919 = v_coeffs_16 + int(23);
    *_S919 = *_S919 + pSH8_10 * _S908;
    float v_x_8 = _S912 + _S908 * (pSH4_x_2 * *_S852 + pSH8_x_8 * *_S856 + fTmp0B_13 * *_S855);
    float v_y_8 = _S913 + _S908 * (pSH8_x_8 * *_S852 + pSH8_y_2 * *_S856 + fTmp0B_13 * *_S853);
    float v_z_8 = _S914 + _S908 * (pSH6_z_4 * *_S854 + pSH7_z_2 * *_S855 + pSH5_z_2 * *_S853);
    float * _S920 = v_coeffs_16 + int(26);
    *_S920 = *_S920 + pSH9_5 * _S908;
    float * _S921 = v_coeffs_16 + int(29);
    *_S921 = *_S921 + pSH10_5 * _S908;
    float * _S922 = v_coeffs_16 + int(32);
    *_S922 = *_S922 + pSH11_5 * _S908;
    float * _S923 = v_coeffs_16 + int(35);
    *_S923 = *_S923 + pSH12_6 * _S908;
    float * _S924 = v_coeffs_16 + int(38);
    *_S924 = *_S924 + pSH13_5 * _S908;
    float * _S925 = v_coeffs_16 + int(41);
    *_S925 = *_S925 + pSH14_5 * _S908;
    float * _S926 = v_coeffs_16 + int(44);
    *_S926 = *_S926 + pSH15_5 * _S908;
    float3  v_dir_n_22 = make_float3 (v_x_8 + _S908 * (pSH9_x_0 * *_S857 + pSH15_x_0 * *_S863 + pSH10_x_0 * *_S858 + pSH14_x_4 * *_S862 + fTmp0C_7 * *_S861), v_y_8 + _S908 * (pSH9_y_0 * *_S857 + pSH15_y_0 * *_S863 + pSH14_x_4 * *_S858 + pSH14_y_0 * *_S862 + fTmp0C_7 * *_S859), v_z_8 + _S908 * (pSH12_z_2 * *_S860 + pSH13_z_0 * *_S861 + pSH11_z_0 * *_S859 + pSH14_z_0 * *_S862 + pSH10_z_0 * *_S858));
    float3  v_viewdir_39 = v_viewdir_38 + (v_dir_n_22 - make_float3 (dot_0(v_dir_n_22, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    Matrix<float, 3, 3>  _S927 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S928;
    (&_S928)->primal_0 = _S812;
    (&_S928)->differential_0 = _S927;
    float3  _S929 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S930;
    (&_S930)->primal_0 = t_25;
    (&_S930)->differential_0 = _S929;
    s_bwd_prop_mul_0(&_S928, &_S930, v_viewdir_39);
    Matrix<float, 3, 3>  _S931 = transpose_0(_S928.differential_0);
    *v_mean_16 = *v_mean_16 + v_viewdir_39;
    *v_R_16 = *v_R_16 + _S931;
    *v_t_16 = *v_t_16 + _S930.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_atomic(float3  mean_26, Matrix<float, 3, 3>  R_26, float3  t_26, float3  coeff_dc_26, float * coeffs_26, float3  v_colors_17, float3  * v_coeff_dc_17, float * v_coeffs_17, float3  * v_mean_17, Matrix<float, 3, 3>  * v_R_17, float3  * v_t_17)
{
    Matrix<float, 3, 3>  _S932 = transpose_0(R_26);
    float3  _S933 = mean_26 + mul_0(_S932, t_26);
    float _S934 = _S933.x;
    float _S935 = _S933.y;
    float _S936 = _S933.z;
    float inv_norm_20 = (F32_rsqrt((_S934 * _S934 + _S935 * _S935 + _S936 * _S936)));
    float x_22 = _S934 * inv_norm_20;
    float y_20 = _S935 * inv_norm_20;
    float z_19 = _S936 * inv_norm_20;
    float _S937 = - y_20;
    float * _S938 = coeffs_26 + int(0);
    float * _S939 = coeffs_26 + int(3);
    float * _S940 = coeffs_26 + int(6);
    float z2_8 = z_19 * z_19;
    float fTmp0B_14 = -1.09254848957061768f * z_19;
    float fC1_8 = x_22 * x_22 - y_20 * y_20;
    float _S941 = 2.0f * x_22;
    float fS1_8 = _S941 * y_20;
    float pSH6_12 = 0.94617468118667603f * z2_8 - 0.31539157032966614f;
    float pSH7_11 = fTmp0B_14 * x_22;
    float pSH5_11 = fTmp0B_14 * y_20;
    float pSH8_11 = 0.54627424478530884f * fC1_8;
    float pSH4_11 = 0.54627424478530884f * fS1_8;
    float * _S942 = coeffs_26 + int(9);
    float * _S943 = coeffs_26 + int(12);
    float * _S944 = coeffs_26 + int(15);
    float * _S945 = coeffs_26 + int(18);
    float * _S946 = coeffs_26 + int(21);
    float fTmp0C_8 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_19;
    float pSH12_7 = z_19 * (1.86588168144226074f * z2_8 - 1.11952900886535645f);
    float pSH13_6 = fTmp0C_8 * x_22;
    float pSH11_6 = fTmp0C_8 * y_20;
    float pSH14_6 = fTmp1B_8 * fC1_8;
    float pSH10_6 = fTmp1B_8 * fS1_8;
    float pSH15_6 = -0.59004360437393188f * (x_22 * fC1_8 - y_20 * fS1_8);
    float pSH9_6 = -0.59004360437393188f * (x_22 * fS1_8 + y_20 * fC1_8);
    float * _S947 = coeffs_26 + int(24);
    float * _S948 = coeffs_26 + int(27);
    float * _S949 = coeffs_26 + int(30);
    float * _S950 = coeffs_26 + int(33);
    float * _S951 = coeffs_26 + int(36);
    float * _S952 = coeffs_26 + int(39);
    float * _S953 = coeffs_26 + int(42);
    float * _S954 = coeffs_26 + int(1);
    float * _S955 = coeffs_26 + int(4);
    float * _S956 = coeffs_26 + int(7);
    float * _S957 = coeffs_26 + int(10);
    float * _S958 = coeffs_26 + int(13);
    float * _S959 = coeffs_26 + int(16);
    float * _S960 = coeffs_26 + int(19);
    float * _S961 = coeffs_26 + int(22);
    float * _S962 = coeffs_26 + int(25);
    float * _S963 = coeffs_26 + int(28);
    float * _S964 = coeffs_26 + int(31);
    float * _S965 = coeffs_26 + int(34);
    float * _S966 = coeffs_26 + int(37);
    float * _S967 = coeffs_26 + int(40);
    float * _S968 = coeffs_26 + int(43);
    float * _S969 = coeffs_26 + int(2);
    float * _S970 = coeffs_26 + int(5);
    float * _S971 = coeffs_26 + int(8);
    float * _S972 = coeffs_26 + int(11);
    float * _S973 = coeffs_26 + int(14);
    float * _S974 = coeffs_26 + int(17);
    float * _S975 = coeffs_26 + int(20);
    float * _S976 = coeffs_26 + int(23);
    float * _S977 = coeffs_26 + int(26);
    float * _S978 = coeffs_26 + int(29);
    float * _S979 = coeffs_26 + int(32);
    float * _S980 = coeffs_26 + int(35);
    float * _S981 = coeffs_26 + int(38);
    float * _S982 = coeffs_26 + int(41);
    float * _S983 = coeffs_26 + int(44);
    float3  _S984 = v_colors_17 * make_float3 (float((0.282094806432724f * coeff_dc_26.x + 0.48860251903533936f * (_S937 * *_S938 + z_19 * *_S939 - x_22 * *_S940) + (pSH4_11 * *_S942 + pSH5_11 * *_S943 + pSH6_12 * *_S944 + pSH7_11 * *_S945 + pSH8_11 * *_S946) + (pSH9_6 * *_S947 + pSH10_6 * *_S948 + pSH11_6 * *_S949 + pSH12_7 * *_S950 + pSH13_6 * *_S951 + pSH14_6 * *_S952 + pSH15_6 * *_S953)) >= -0.5f), float((0.282094806432724f * coeff_dc_26.y + 0.48860251903533936f * (_S937 * *_S954 + z_19 * *_S955 - x_22 * *_S956) + (pSH4_11 * *_S957 + pSH5_11 * *_S958 + pSH6_12 * *_S959 + pSH7_11 * *_S960 + pSH8_11 * *_S961) + (pSH9_6 * *_S962 + pSH10_6 * *_S963 + pSH11_6 * *_S964 + pSH12_7 * *_S965 + pSH13_6 * *_S966 + pSH14_6 * *_S967 + pSH15_6 * *_S968)) >= -0.5f), float((0.282094806432724f * coeff_dc_26.z + 0.48860251903533936f * (_S937 * *_S969 + z_19 * *_S970 - x_22 * *_S971) + (pSH4_11 * *_S972 + pSH5_11 * *_S973 + pSH6_12 * *_S974 + pSH7_11 * *_S975 + pSH8_11 * *_S976) + (pSH9_6 * *_S977 + pSH10_6 * *_S978 + pSH11_6 * *_S979 + pSH12_7 * *_S980 + pSH13_6 * *_S981 + pSH14_6 * *_S982 + pSH15_6 * *_S983)) >= -0.5f));
    float3  v_viewdir_40 = {};
    float _S985 = _S984.x;
    *&(v_coeff_dc_17->x) = *&(v_coeff_dc_17->x) + 0.282094806432724f * _S985;
    float _S986 = -0.48860251903533936f * y_20;
    float temp_83 = _S986 * _S985;
    bool _S987;
    if((F32_isfinite((temp_83))))
    {
        _S987 = temp_83 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S988 = atomicAdd(v_coeffs_17 + int(0), temp_83);
    }
    float _S989 = 0.48860251903533936f * z_19;
    float temp_84 = _S989 * _S985;
    if((F32_isfinite((temp_84))))
    {
        _S987 = temp_84 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S990 = atomicAdd(v_coeffs_17 + int(3), temp_84);
    }
    float _S991 = -0.48860251903533936f * x_22;
    float temp_85 = _S991 * _S985;
    if((F32_isfinite((temp_85))))
    {
        _S987 = temp_85 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S992 = atomicAdd(v_coeffs_17 + int(6), temp_85);
    }
    float _S993 = -0.48860251903533936f * *_S940 * _S985;
    float _S994 = -0.48860251903533936f * *_S938 * _S985;
    float _S995 = 0.48860251903533936f * *_S939 * _S985;
    float temp_86 = pSH4_11 * _S985;
    if((F32_isfinite((temp_86))))
    {
        _S987 = temp_86 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S996 = atomicAdd(v_coeffs_17 + int(9), temp_86);
    }
    float temp_87 = pSH5_11 * _S985;
    if((F32_isfinite((temp_87))))
    {
        _S987 = temp_87 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S997 = atomicAdd(v_coeffs_17 + int(12), temp_87);
    }
    float temp_88 = pSH6_12 * _S985;
    if((F32_isfinite((temp_88))))
    {
        _S987 = temp_88 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S998 = atomicAdd(v_coeffs_17 + int(15), temp_88);
    }
    float temp_89 = pSH7_11 * _S985;
    if((F32_isfinite((temp_89))))
    {
        _S987 = temp_89 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S999 = atomicAdd(v_coeffs_17 + int(18), temp_89);
    }
    float temp_90 = pSH8_11 * _S985;
    if((F32_isfinite((temp_90))))
    {
        _S987 = temp_90 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1000 = atomicAdd(v_coeffs_17 + int(21), temp_90);
    }
    float fC1_y_5 = -2.0f * y_20;
    float fS1_x_5 = 2.0f * y_20;
    float pSH6_z_5 = 1.89234936237335205f * z_19;
    float pSH7_z_3 = -1.09254848957061768f * x_22;
    float pSH5_z_3 = -1.09254848957061768f * y_20;
    float pSH8_x_9 = 0.54627424478530884f * _S941;
    float pSH8_y_3 = 0.54627424478530884f * fC1_y_5;
    float pSH4_x_3 = 0.54627424478530884f * fS1_x_5;
    float v_x_9 = _S993 + _S985 * (pSH4_x_3 * *_S942 + pSH8_x_9 * *_S946 + fTmp0B_14 * *_S945);
    float v_y_9 = _S994 + _S985 * (pSH8_x_9 * *_S942 + pSH8_y_3 * *_S946 + fTmp0B_14 * *_S943);
    float v_z_9 = _S995 + _S985 * (pSH6_z_5 * *_S944 + pSH7_z_3 * *_S945 + pSH5_z_3 * *_S943);
    float temp_91 = pSH9_6 * _S985;
    if((F32_isfinite((temp_91))))
    {
        _S987 = temp_91 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1001 = atomicAdd(v_coeffs_17 + int(24), temp_91);
    }
    float temp_92 = pSH10_6 * _S985;
    if((F32_isfinite((temp_92))))
    {
        _S987 = temp_92 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1002 = atomicAdd(v_coeffs_17 + int(27), temp_92);
    }
    float temp_93 = pSH11_6 * _S985;
    if((F32_isfinite((temp_93))))
    {
        _S987 = temp_93 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1003 = atomicAdd(v_coeffs_17 + int(30), temp_93);
    }
    float temp_94 = pSH12_7 * _S985;
    if((F32_isfinite((temp_94))))
    {
        _S987 = temp_94 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1004 = atomicAdd(v_coeffs_17 + int(33), temp_94);
    }
    float temp_95 = pSH13_6 * _S985;
    if((F32_isfinite((temp_95))))
    {
        _S987 = temp_95 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1005 = atomicAdd(v_coeffs_17 + int(36), temp_95);
    }
    float temp_96 = pSH14_6 * _S985;
    if((F32_isfinite((temp_96))))
    {
        _S987 = temp_96 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1006 = atomicAdd(v_coeffs_17 + int(39), temp_96);
    }
    float temp_97 = pSH15_6 * _S985;
    if((F32_isfinite((temp_97))))
    {
        _S987 = temp_97 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1007 = atomicAdd(v_coeffs_17 + int(42), temp_97);
    }
    float fTmp0C_z_5 = -4.57045793533325195f * z_19;
    float _S1008 = x_22 * _S941;
    float _S1009 = y_20 * _S941;
    float pSH12_z_3 = 5.59764480590820312f * z2_8 - 1.11952900886535645f;
    float pSH13_z_1 = fTmp0C_z_5 * x_22;
    float pSH11_z_1 = fTmp0C_z_5 * y_20;
    float pSH14_x_5 = fTmp1B_8 * _S941;
    float pSH14_y_1 = fTmp1B_8 * fC1_y_5;
    float pSH14_z_1 = 1.44530570507049561f * fC1_8;
    float pSH10_x_1 = fTmp1B_8 * fS1_x_5;
    float pSH10_z_1 = 1.44530570507049561f * fS1_8;
    float pSH15_x_1 = -0.59004360437393188f * (fC1_8 + _S1008 - y_20 * fS1_x_5);
    float pSH15_y_1 = -0.59004360437393188f * (x_22 * fC1_y_5 - fS1_8 - _S1009);
    float pSH9_x_1 = -0.59004360437393188f * (fS1_8 + x_22 * fS1_x_5 + _S1009);
    float pSH9_y_1 = -0.59004360437393188f * (_S1008 + fC1_8 + y_20 * fC1_y_5);
    float3  dir_n_13 = make_float3 (x_22, y_20, z_19);
    float3  v_dir_n_23 = make_float3 (v_x_9 + _S985 * (pSH9_x_1 * *_S947 + pSH15_x_1 * *_S953 + pSH10_x_1 * *_S948 + pSH14_x_5 * *_S952 + fTmp0C_8 * *_S951), v_y_9 + _S985 * (pSH9_y_1 * *_S947 + pSH15_y_1 * *_S953 + pSH14_x_5 * *_S948 + pSH14_y_1 * *_S952 + fTmp0C_8 * *_S949), v_z_9 + _S985 * (pSH12_z_3 * *_S950 + pSH13_z_1 * *_S951 + pSH11_z_1 * *_S949 + pSH14_z_1 * *_S952 + pSH10_z_1 * *_S948));
    float3  v_viewdir_41 = v_viewdir_40 + (v_dir_n_23 - make_float3 (dot_0(v_dir_n_23, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    float _S1010 = _S984.y;
    *&(v_coeff_dc_17->y) = *&(v_coeff_dc_17->y) + 0.282094806432724f * _S1010;
    float temp_98 = _S986 * _S1010;
    if((F32_isfinite((temp_98))))
    {
        _S987 = temp_98 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1011 = atomicAdd(v_coeffs_17 + int(1), temp_98);
    }
    float temp_99 = _S989 * _S1010;
    if((F32_isfinite((temp_99))))
    {
        _S987 = temp_99 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1012 = atomicAdd(v_coeffs_17 + int(4), temp_99);
    }
    float temp_100 = _S991 * _S1010;
    if((F32_isfinite((temp_100))))
    {
        _S987 = temp_100 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1013 = atomicAdd(v_coeffs_17 + int(7), temp_100);
    }
    float _S1014 = -0.48860251903533936f * *_S956 * _S1010;
    float _S1015 = -0.48860251903533936f * *_S954 * _S1010;
    float _S1016 = 0.48860251903533936f * *_S955 * _S1010;
    float temp_101 = pSH4_11 * _S1010;
    if((F32_isfinite((temp_101))))
    {
        _S987 = temp_101 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1017 = atomicAdd(v_coeffs_17 + int(10), temp_101);
    }
    float temp_102 = pSH5_11 * _S1010;
    if((F32_isfinite((temp_102))))
    {
        _S987 = temp_102 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1018 = atomicAdd(v_coeffs_17 + int(13), temp_102);
    }
    float temp_103 = pSH6_12 * _S1010;
    if((F32_isfinite((temp_103))))
    {
        _S987 = temp_103 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1019 = atomicAdd(v_coeffs_17 + int(16), temp_103);
    }
    float temp_104 = pSH7_11 * _S1010;
    if((F32_isfinite((temp_104))))
    {
        _S987 = temp_104 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1020 = atomicAdd(v_coeffs_17 + int(19), temp_104);
    }
    float temp_105 = pSH8_11 * _S1010;
    if((F32_isfinite((temp_105))))
    {
        _S987 = temp_105 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1021 = atomicAdd(v_coeffs_17 + int(22), temp_105);
    }
    float v_x_10 = _S1014 + _S1010 * (pSH4_x_3 * *_S957 + pSH8_x_9 * *_S961 + fTmp0B_14 * *_S960);
    float v_y_10 = _S1015 + _S1010 * (pSH8_x_9 * *_S957 + pSH8_y_3 * *_S961 + fTmp0B_14 * *_S958);
    float v_z_10 = _S1016 + _S1010 * (pSH6_z_5 * *_S959 + pSH7_z_3 * *_S960 + pSH5_z_3 * *_S958);
    float temp_106 = pSH9_6 * _S1010;
    if((F32_isfinite((temp_106))))
    {
        _S987 = temp_106 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1022 = atomicAdd(v_coeffs_17 + int(25), temp_106);
    }
    float temp_107 = pSH10_6 * _S1010;
    if((F32_isfinite((temp_107))))
    {
        _S987 = temp_107 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1023 = atomicAdd(v_coeffs_17 + int(28), temp_107);
    }
    float temp_108 = pSH11_6 * _S1010;
    if((F32_isfinite((temp_108))))
    {
        _S987 = temp_108 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1024 = atomicAdd(v_coeffs_17 + int(31), temp_108);
    }
    float temp_109 = pSH12_7 * _S1010;
    if((F32_isfinite((temp_109))))
    {
        _S987 = temp_109 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1025 = atomicAdd(v_coeffs_17 + int(34), temp_109);
    }
    float temp_110 = pSH13_6 * _S1010;
    if((F32_isfinite((temp_110))))
    {
        _S987 = temp_110 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1026 = atomicAdd(v_coeffs_17 + int(37), temp_110);
    }
    float temp_111 = pSH14_6 * _S1010;
    if((F32_isfinite((temp_111))))
    {
        _S987 = temp_111 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1027 = atomicAdd(v_coeffs_17 + int(40), temp_111);
    }
    float temp_112 = pSH15_6 * _S1010;
    if((F32_isfinite((temp_112))))
    {
        _S987 = temp_112 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1028 = atomicAdd(v_coeffs_17 + int(43), temp_112);
    }
    float3  v_dir_n_24 = make_float3 (v_x_10 + _S1010 * (pSH9_x_1 * *_S962 + pSH15_x_1 * *_S968 + pSH10_x_1 * *_S963 + pSH14_x_5 * *_S967 + fTmp0C_8 * *_S966), v_y_10 + _S1010 * (pSH9_y_1 * *_S962 + pSH15_y_1 * *_S968 + pSH14_x_5 * *_S963 + pSH14_y_1 * *_S967 + fTmp0C_8 * *_S964), v_z_10 + _S1010 * (pSH12_z_3 * *_S965 + pSH13_z_1 * *_S966 + pSH11_z_1 * *_S964 + pSH14_z_1 * *_S967 + pSH10_z_1 * *_S963));
    float3  v_viewdir_42 = v_viewdir_41 + (v_dir_n_24 - make_float3 (dot_0(v_dir_n_24, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    float _S1029 = _S984.z;
    *&(v_coeff_dc_17->z) = *&(v_coeff_dc_17->z) + 0.282094806432724f * _S1029;
    float temp_113 = _S986 * _S1029;
    if((F32_isfinite((temp_113))))
    {
        _S987 = temp_113 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1030 = atomicAdd(v_coeffs_17 + int(2), temp_113);
    }
    float temp_114 = _S989 * _S1029;
    if((F32_isfinite((temp_114))))
    {
        _S987 = temp_114 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1031 = atomicAdd(v_coeffs_17 + int(5), temp_114);
    }
    float temp_115 = _S991 * _S1029;
    if((F32_isfinite((temp_115))))
    {
        _S987 = temp_115 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1032 = atomicAdd(v_coeffs_17 + int(8), temp_115);
    }
    float _S1033 = -0.48860251903533936f * *_S971 * _S1029;
    float _S1034 = -0.48860251903533936f * *_S969 * _S1029;
    float _S1035 = 0.48860251903533936f * *_S970 * _S1029;
    float temp_116 = pSH4_11 * _S1029;
    if((F32_isfinite((temp_116))))
    {
        _S987 = temp_116 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1036 = atomicAdd(v_coeffs_17 + int(11), temp_116);
    }
    float temp_117 = pSH5_11 * _S1029;
    if((F32_isfinite((temp_117))))
    {
        _S987 = temp_117 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1037 = atomicAdd(v_coeffs_17 + int(14), temp_117);
    }
    float temp_118 = pSH6_12 * _S1029;
    if((F32_isfinite((temp_118))))
    {
        _S987 = temp_118 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1038 = atomicAdd(v_coeffs_17 + int(17), temp_118);
    }
    float temp_119 = pSH7_11 * _S1029;
    if((F32_isfinite((temp_119))))
    {
        _S987 = temp_119 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1039 = atomicAdd(v_coeffs_17 + int(20), temp_119);
    }
    float temp_120 = pSH8_11 * _S1029;
    if((F32_isfinite((temp_120))))
    {
        _S987 = temp_120 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1040 = atomicAdd(v_coeffs_17 + int(23), temp_120);
    }
    float v_x_11 = _S1033 + _S1029 * (pSH4_x_3 * *_S972 + pSH8_x_9 * *_S976 + fTmp0B_14 * *_S975);
    float v_y_11 = _S1034 + _S1029 * (pSH8_x_9 * *_S972 + pSH8_y_3 * *_S976 + fTmp0B_14 * *_S973);
    float v_z_11 = _S1035 + _S1029 * (pSH6_z_5 * *_S974 + pSH7_z_3 * *_S975 + pSH5_z_3 * *_S973);
    float temp_121 = pSH9_6 * _S1029;
    if((F32_isfinite((temp_121))))
    {
        _S987 = temp_121 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1041 = atomicAdd(v_coeffs_17 + int(26), temp_121);
    }
    float temp_122 = pSH10_6 * _S1029;
    if((F32_isfinite((temp_122))))
    {
        _S987 = temp_122 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1042 = atomicAdd(v_coeffs_17 + int(29), temp_122);
    }
    float temp_123 = pSH11_6 * _S1029;
    if((F32_isfinite((temp_123))))
    {
        _S987 = temp_123 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1043 = atomicAdd(v_coeffs_17 + int(32), temp_123);
    }
    float temp_124 = pSH12_7 * _S1029;
    if((F32_isfinite((temp_124))))
    {
        _S987 = temp_124 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1044 = atomicAdd(v_coeffs_17 + int(35), temp_124);
    }
    float temp_125 = pSH13_6 * _S1029;
    if((F32_isfinite((temp_125))))
    {
        _S987 = temp_125 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1045 = atomicAdd(v_coeffs_17 + int(38), temp_125);
    }
    float temp_126 = pSH14_6 * _S1029;
    if((F32_isfinite((temp_126))))
    {
        _S987 = temp_126 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1046 = atomicAdd(v_coeffs_17 + int(41), temp_126);
    }
    float temp_127 = pSH15_6 * _S1029;
    if((F32_isfinite((temp_127))))
    {
        _S987 = temp_127 != 0.0f;
    }
    else
    {
        _S987 = false;
    }
    if(_S987)
    {
        float _S1047 = atomicAdd(v_coeffs_17 + int(44), temp_127);
    }
    float3  v_dir_n_25 = make_float3 (v_x_11 + _S1029 * (pSH9_x_1 * *_S977 + pSH15_x_1 * *_S983 + pSH10_x_1 * *_S978 + pSH14_x_5 * *_S982 + fTmp0C_8 * *_S981), v_y_11 + _S1029 * (pSH9_y_1 * *_S977 + pSH15_y_1 * *_S983 + pSH14_x_5 * *_S978 + pSH14_y_1 * *_S982 + fTmp0C_8 * *_S979), v_z_11 + _S1029 * (pSH12_z_3 * *_S980 + pSH13_z_1 * *_S981 + pSH11_z_1 * *_S979 + pSH14_z_1 * *_S982 + pSH10_z_1 * *_S978));
    float3  v_viewdir_43 = v_viewdir_42 + (v_dir_n_25 - make_float3 (dot_0(v_dir_n_25, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    Matrix<float, 3, 3>  _S1048 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1049;
    (&_S1049)->primal_0 = _S932;
    (&_S1049)->differential_0 = _S1048;
    float3  _S1050 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1051;
    (&_S1051)->primal_0 = t_26;
    (&_S1051)->differential_0 = _S1050;
    s_bwd_prop_mul_0(&_S1049, &_S1051, v_viewdir_43);
    Matrix<float, 3, 3>  _S1052 = transpose_0(_S1049.differential_0);
    *v_mean_17 = *v_mean_17 + v_viewdir_43;
    *v_R_17 = *v_R_17 + _S1052;
    *v_t_17 = *v_t_17 + _S1051.differential_0;
    return;
}

inline __device__ float3  sh4_to_color(float3  mean_27, Matrix<float, 3, 3>  R_27, float3  t_27, float3  coeff_dc_27, float * coeffs_27)
{
    float3  _S1053 = mean_27 + mul_0(transpose_0(R_27), t_27);
    float _S1054 = _S1053.x;
    float _S1055 = _S1053.y;
    float _S1056 = _S1053.z;
    float inv_norm_21 = (F32_rsqrt((_S1054 * _S1054 + _S1055 * _S1055 + _S1056 * _S1056)));
    float x_23 = _S1054 * inv_norm_21;
    float y_21 = _S1055 * inv_norm_21;
    float z_20 = _S1056 * inv_norm_21;
    float _S1057 = - y_21;
    float z2_9 = z_20 * z_20;
    float fTmp0B_15 = -1.09254848957061768f * z_20;
    float fC1_9 = x_23 * x_23 - y_21 * y_21;
    float fS1_9 = 2.0f * x_23 * y_21;
    float pSH6_13 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
    float pSH7_12 = fTmp0B_15 * x_23;
    float pSH5_12 = fTmp0B_15 * y_21;
    float pSH8_12 = 0.54627424478530884f * fC1_9;
    float pSH4_12 = 0.54627424478530884f * fS1_9;
    float fTmp0C_9 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_20;
    float fC2_3 = x_23 * fC1_9 - y_21 * fS1_9;
    float fS2_3 = x_23 * fS1_9 + y_21 * fC1_9;
    float pSH12_8 = z_20 * (1.86588168144226074f * z2_9 - 1.11952900886535645f);
    float pSH13_7 = fTmp0C_9 * x_23;
    float pSH11_7 = fTmp0C_9 * y_21;
    float pSH14_7 = fTmp1B_9 * fC1_9;
    float pSH10_7 = fTmp1B_9 * fS1_9;
    float pSH15_7 = -0.59004360437393188f * fC2_3;
    float pSH9_7 = -0.59004360437393188f * fS2_3;
    float fTmp0D_3 = z_20 * (-4.68332576751708984f * z2_9 + 2.00713968276977539f);
    float fTmp1C_3 = 3.31161141395568848f * z2_9 - 0.47308734059333801f;
    float fTmp2B_3 = -1.77013075351715088f * z_20;
    float pSH20_2 = 1.9843134880065918f * z_20 * pSH12_8 - 1.00623059272766113f * pSH6_13;
    float pSH21_2 = fTmp0D_3 * x_23;
    float pSH19_2 = fTmp0D_3 * y_21;
    float pSH22_2 = fTmp1C_3 * fC1_9;
    float pSH18_2 = fTmp1C_3 * fS1_9;
    float pSH23_2 = fTmp2B_3 * fC2_3;
    float pSH17_2 = fTmp2B_3 * fS2_3;
    float pSH24_2 = 0.62583571672439575f * (x_23 * fC2_3 - y_21 * fS2_3);
    float pSH16_2 = 0.62583571672439575f * (x_23 * fS2_3 + y_21 * fC2_3);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_27.x + 0.48860251903533936f * (_S1057 * *(coeffs_27 + int(0)) + z_20 * *(coeffs_27 + int(3)) - x_23 * *(coeffs_27 + int(6))) + (pSH4_12 * *(coeffs_27 + int(9)) + pSH5_12 * *(coeffs_27 + int(12)) + pSH6_13 * *(coeffs_27 + int(15)) + pSH7_12 * *(coeffs_27 + int(18)) + pSH8_12 * *(coeffs_27 + int(21))) + (pSH9_7 * *(coeffs_27 + int(24)) + pSH10_7 * *(coeffs_27 + int(27)) + pSH11_7 * *(coeffs_27 + int(30)) + pSH12_8 * *(coeffs_27 + int(33)) + pSH13_7 * *(coeffs_27 + int(36)) + pSH14_7 * *(coeffs_27 + int(39)) + pSH15_7 * *(coeffs_27 + int(42))) + (pSH16_2 * *(coeffs_27 + int(45)) + pSH17_2 * *(coeffs_27 + int(48)) + pSH18_2 * *(coeffs_27 + int(51)) + pSH19_2 * *(coeffs_27 + int(54)) + pSH20_2 * *(coeffs_27 + int(57)) + pSH21_2 * *(coeffs_27 + int(60)) + pSH22_2 * *(coeffs_27 + int(63)) + pSH23_2 * *(coeffs_27 + int(66)) + pSH24_2 * *(coeffs_27 + int(69))), 0.282094806432724f * coeff_dc_27.y + 0.48860251903533936f * (_S1057 * *(coeffs_27 + int(1)) + z_20 * *(coeffs_27 + int(4)) - x_23 * *(coeffs_27 + int(7))) + (pSH4_12 * *(coeffs_27 + int(10)) + pSH5_12 * *(coeffs_27 + int(13)) + pSH6_13 * *(coeffs_27 + int(16)) + pSH7_12 * *(coeffs_27 + int(19)) + pSH8_12 * *(coeffs_27 + int(22))) + (pSH9_7 * *(coeffs_27 + int(25)) + pSH10_7 * *(coeffs_27 + int(28)) + pSH11_7 * *(coeffs_27 + int(31)) + pSH12_8 * *(coeffs_27 + int(34)) + pSH13_7 * *(coeffs_27 + int(37)) + pSH14_7 * *(coeffs_27 + int(40)) + pSH15_7 * *(coeffs_27 + int(43))) + (pSH16_2 * *(coeffs_27 + int(46)) + pSH17_2 * *(coeffs_27 + int(49)) + pSH18_2 * *(coeffs_27 + int(52)) + pSH19_2 * *(coeffs_27 + int(55)) + pSH20_2 * *(coeffs_27 + int(58)) + pSH21_2 * *(coeffs_27 + int(61)) + pSH22_2 * *(coeffs_27 + int(64)) + pSH23_2 * *(coeffs_27 + int(67)) + pSH24_2 * *(coeffs_27 + int(70))), 0.282094806432724f * coeff_dc_27.z + 0.48860251903533936f * (_S1057 * *(coeffs_27 + int(2)) + z_20 * *(coeffs_27 + int(5)) - x_23 * *(coeffs_27 + int(8))) + (pSH4_12 * *(coeffs_27 + int(11)) + pSH5_12 * *(coeffs_27 + int(14)) + pSH6_13 * *(coeffs_27 + int(17)) + pSH7_12 * *(coeffs_27 + int(20)) + pSH8_12 * *(coeffs_27 + int(23))) + (pSH9_7 * *(coeffs_27 + int(26)) + pSH10_7 * *(coeffs_27 + int(29)) + pSH11_7 * *(coeffs_27 + int(32)) + pSH12_8 * *(coeffs_27 + int(35)) + pSH13_7 * *(coeffs_27 + int(38)) + pSH14_7 * *(coeffs_27 + int(41)) + pSH15_7 * *(coeffs_27 + int(44))) + (pSH16_2 * *(coeffs_27 + int(47)) + pSH17_2 * *(coeffs_27 + int(50)) + pSH18_2 * *(coeffs_27 + int(53)) + pSH19_2 * *(coeffs_27 + int(56)) + pSH20_2 * *(coeffs_27 + int(59)) + pSH21_2 * *(coeffs_27 + int(62)) + pSH22_2 * *(coeffs_27 + int(65)) + pSH23_2 * *(coeffs_27 + int(68)) + pSH24_2 * *(coeffs_27 + int(71)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh4_to_color_vjp_inplace(float3  mean_28, Matrix<float, 3, 3>  R_28, float3  t_28, float3  coeff_dc_28, float * coeffs_28, float3  v_colors_18, float3  * v_coeff_dc_18, float * v_coeffs_18, float3  * v_mean_18, Matrix<float, 3, 3>  * v_R_18, float3  * v_t_18)
{
    Matrix<float, 3, 3>  _S1058 = transpose_0(R_28);
    float3  _S1059 = mean_28 + mul_0(_S1058, t_28);
    float _S1060 = _S1059.x;
    float _S1061 = _S1059.y;
    float _S1062 = _S1059.z;
    float inv_norm_22 = (F32_rsqrt((_S1060 * _S1060 + _S1061 * _S1061 + _S1062 * _S1062)));
    float x_24 = _S1060 * inv_norm_22;
    float y_22 = _S1061 * inv_norm_22;
    float z_21 = _S1062 * inv_norm_22;
    float _S1063 = - y_22;
    float * _S1064 = coeffs_28 + int(0);
    float * _S1065 = coeffs_28 + int(3);
    float * _S1066 = coeffs_28 + int(6);
    float z2_10 = z_21 * z_21;
    float fTmp0B_16 = -1.09254848957061768f * z_21;
    float fC1_10 = x_24 * x_24 - y_22 * y_22;
    float _S1067 = 2.0f * x_24;
    float fS1_10 = _S1067 * y_22;
    float pSH6_14 = 0.94617468118667603f * z2_10 - 0.31539157032966614f;
    float pSH7_13 = fTmp0B_16 * x_24;
    float pSH5_13 = fTmp0B_16 * y_22;
    float pSH8_13 = 0.54627424478530884f * fC1_10;
    float pSH4_13 = 0.54627424478530884f * fS1_10;
    float * _S1068 = coeffs_28 + int(9);
    float * _S1069 = coeffs_28 + int(12);
    float * _S1070 = coeffs_28 + int(15);
    float * _S1071 = coeffs_28 + int(18);
    float * _S1072 = coeffs_28 + int(21);
    float fTmp0C_10 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_21;
    float fC2_4 = x_24 * fC1_10 - y_22 * fS1_10;
    float fS2_4 = x_24 * fS1_10 + y_22 * fC1_10;
    float pSH12_9 = z_21 * (1.86588168144226074f * z2_10 - 1.11952900886535645f);
    float pSH13_8 = fTmp0C_10 * x_24;
    float pSH11_8 = fTmp0C_10 * y_22;
    float pSH14_8 = fTmp1B_10 * fC1_10;
    float pSH10_8 = fTmp1B_10 * fS1_10;
    float pSH15_8 = -0.59004360437393188f * fC2_4;
    float pSH9_8 = -0.59004360437393188f * fS2_4;
    float * _S1073 = coeffs_28 + int(24);
    float * _S1074 = coeffs_28 + int(27);
    float * _S1075 = coeffs_28 + int(30);
    float * _S1076 = coeffs_28 + int(33);
    float * _S1077 = coeffs_28 + int(36);
    float * _S1078 = coeffs_28 + int(39);
    float * _S1079 = coeffs_28 + int(42);
    float fTmp0D_4 = z_21 * (-4.68332576751708984f * z2_10 + 2.00713968276977539f);
    float fTmp1C_4 = 3.31161141395568848f * z2_10 - 0.47308734059333801f;
    float fTmp2B_4 = -1.77013075351715088f * z_21;
    float _S1080 = 1.9843134880065918f * z_21 * pSH12_9;
    float pSH20_3 = _S1080 - 1.00623059272766113f * pSH6_14;
    float pSH21_3 = fTmp0D_4 * x_24;
    float pSH19_3 = fTmp0D_4 * y_22;
    float pSH22_3 = fTmp1C_4 * fC1_10;
    float pSH18_3 = fTmp1C_4 * fS1_10;
    float pSH23_3 = fTmp2B_4 * fC2_4;
    float pSH17_3 = fTmp2B_4 * fS2_4;
    float pSH24_3 = 0.62583571672439575f * (x_24 * fC2_4 - y_22 * fS2_4);
    float pSH16_3 = 0.62583571672439575f * (x_24 * fS2_4 + y_22 * fC2_4);
    float * _S1081 = coeffs_28 + int(45);
    float * _S1082 = coeffs_28 + int(48);
    float * _S1083 = coeffs_28 + int(51);
    float * _S1084 = coeffs_28 + int(54);
    float * _S1085 = coeffs_28 + int(57);
    float * _S1086 = coeffs_28 + int(60);
    float * _S1087 = coeffs_28 + int(63);
    float * _S1088 = coeffs_28 + int(66);
    float * _S1089 = coeffs_28 + int(69);
    float * _S1090 = coeffs_28 + int(1);
    float * _S1091 = coeffs_28 + int(4);
    float * _S1092 = coeffs_28 + int(7);
    float * _S1093 = coeffs_28 + int(10);
    float * _S1094 = coeffs_28 + int(13);
    float * _S1095 = coeffs_28 + int(16);
    float * _S1096 = coeffs_28 + int(19);
    float * _S1097 = coeffs_28 + int(22);
    float * _S1098 = coeffs_28 + int(25);
    float * _S1099 = coeffs_28 + int(28);
    float * _S1100 = coeffs_28 + int(31);
    float * _S1101 = coeffs_28 + int(34);
    float * _S1102 = coeffs_28 + int(37);
    float * _S1103 = coeffs_28 + int(40);
    float * _S1104 = coeffs_28 + int(43);
    float * _S1105 = coeffs_28 + int(46);
    float * _S1106 = coeffs_28 + int(49);
    float * _S1107 = coeffs_28 + int(52);
    float * _S1108 = coeffs_28 + int(55);
    float * _S1109 = coeffs_28 + int(58);
    float * _S1110 = coeffs_28 + int(61);
    float * _S1111 = coeffs_28 + int(64);
    float * _S1112 = coeffs_28 + int(67);
    float * _S1113 = coeffs_28 + int(70);
    float * _S1114 = coeffs_28 + int(2);
    float * _S1115 = coeffs_28 + int(5);
    float * _S1116 = coeffs_28 + int(8);
    float * _S1117 = coeffs_28 + int(11);
    float * _S1118 = coeffs_28 + int(14);
    float * _S1119 = coeffs_28 + int(17);
    float * _S1120 = coeffs_28 + int(20);
    float * _S1121 = coeffs_28 + int(23);
    float * _S1122 = coeffs_28 + int(26);
    float * _S1123 = coeffs_28 + int(29);
    float * _S1124 = coeffs_28 + int(32);
    float * _S1125 = coeffs_28 + int(35);
    float * _S1126 = coeffs_28 + int(38);
    float * _S1127 = coeffs_28 + int(41);
    float * _S1128 = coeffs_28 + int(44);
    float * _S1129 = coeffs_28 + int(47);
    float * _S1130 = coeffs_28 + int(50);
    float * _S1131 = coeffs_28 + int(53);
    float * _S1132 = coeffs_28 + int(56);
    float * _S1133 = coeffs_28 + int(59);
    float * _S1134 = coeffs_28 + int(62);
    float * _S1135 = coeffs_28 + int(65);
    float * _S1136 = coeffs_28 + int(68);
    float * _S1137 = coeffs_28 + int(71);
    float3  _S1138 = v_colors_18 * make_float3 (float((0.282094806432724f * coeff_dc_28.x + 0.48860251903533936f * (_S1063 * *_S1064 + z_21 * *_S1065 - x_24 * *_S1066) + (pSH4_13 * *_S1068 + pSH5_13 * *_S1069 + pSH6_14 * *_S1070 + pSH7_13 * *_S1071 + pSH8_13 * *_S1072) + (pSH9_8 * *_S1073 + pSH10_8 * *_S1074 + pSH11_8 * *_S1075 + pSH12_9 * *_S1076 + pSH13_8 * *_S1077 + pSH14_8 * *_S1078 + pSH15_8 * *_S1079) + (pSH16_3 * *_S1081 + pSH17_3 * *_S1082 + pSH18_3 * *_S1083 + pSH19_3 * *_S1084 + pSH20_3 * *_S1085 + pSH21_3 * *_S1086 + pSH22_3 * *_S1087 + pSH23_3 * *_S1088 + pSH24_3 * *_S1089)) >= -0.5f), float((0.282094806432724f * coeff_dc_28.y + 0.48860251903533936f * (_S1063 * *_S1090 + z_21 * *_S1091 - x_24 * *_S1092) + (pSH4_13 * *_S1093 + pSH5_13 * *_S1094 + pSH6_14 * *_S1095 + pSH7_13 * *_S1096 + pSH8_13 * *_S1097) + (pSH9_8 * *_S1098 + pSH10_8 * *_S1099 + pSH11_8 * *_S1100 + pSH12_9 * *_S1101 + pSH13_8 * *_S1102 + pSH14_8 * *_S1103 + pSH15_8 * *_S1104) + (pSH16_3 * *_S1105 + pSH17_3 * *_S1106 + pSH18_3 * *_S1107 + pSH19_3 * *_S1108 + pSH20_3 * *_S1109 + pSH21_3 * *_S1110 + pSH22_3 * *_S1111 + pSH23_3 * *_S1112 + pSH24_3 * *_S1113)) >= -0.5f), float((0.282094806432724f * coeff_dc_28.z + 0.48860251903533936f * (_S1063 * *_S1114 + z_21 * *_S1115 - x_24 * *_S1116) + (pSH4_13 * *_S1117 + pSH5_13 * *_S1118 + pSH6_14 * *_S1119 + pSH7_13 * *_S1120 + pSH8_13 * *_S1121) + (pSH9_8 * *_S1122 + pSH10_8 * *_S1123 + pSH11_8 * *_S1124 + pSH12_9 * *_S1125 + pSH13_8 * *_S1126 + pSH14_8 * *_S1127 + pSH15_8 * *_S1128) + (pSH16_3 * *_S1129 + pSH17_3 * *_S1130 + pSH18_3 * *_S1131 + pSH19_3 * *_S1132 + pSH20_3 * *_S1133 + pSH21_3 * *_S1134 + pSH22_3 * *_S1135 + pSH23_3 * *_S1136 + pSH24_3 * *_S1137)) >= -0.5f));
    float3  v_viewdir_44 = {};
    float _S1139 = _S1138.x;
    *&(v_coeff_dc_18->x) = *&(v_coeff_dc_18->x) + 0.282094806432724f * _S1139;
    float * _S1140 = v_coeffs_18 + int(0);
    float _S1141 = -0.48860251903533936f * y_22;
    *_S1140 = *_S1140 + _S1141 * _S1139;
    float * _S1142 = v_coeffs_18 + int(3);
    float _S1143 = 0.48860251903533936f * z_21;
    *_S1142 = *_S1142 + _S1143 * _S1139;
    float * _S1144 = v_coeffs_18 + int(6);
    float _S1145 = -0.48860251903533936f * x_24;
    *_S1144 = *_S1144 + _S1145 * _S1139;
    float _S1146 = -0.48860251903533936f * *_S1066 * _S1139;
    float _S1147 = -0.48860251903533936f * *_S1064 * _S1139;
    float _S1148 = 0.48860251903533936f * *_S1065 * _S1139;
    float * _S1149 = v_coeffs_18 + int(9);
    *_S1149 = *_S1149 + pSH4_13 * _S1139;
    float * _S1150 = v_coeffs_18 + int(12);
    *_S1150 = *_S1150 + pSH5_13 * _S1139;
    float * _S1151 = v_coeffs_18 + int(15);
    *_S1151 = *_S1151 + pSH6_14 * _S1139;
    float * _S1152 = v_coeffs_18 + int(18);
    *_S1152 = *_S1152 + pSH7_13 * _S1139;
    float * _S1153 = v_coeffs_18 + int(21);
    *_S1153 = *_S1153 + pSH8_13 * _S1139;
    float fC1_y_6 = -2.0f * y_22;
    float fS1_x_6 = 2.0f * y_22;
    float pSH6_z_6 = 1.89234936237335205f * z_21;
    float pSH7_z_4 = -1.09254848957061768f * x_24;
    float pSH5_z_4 = -1.09254848957061768f * y_22;
    float pSH8_x_10 = 0.54627424478530884f * _S1067;
    float pSH8_y_4 = 0.54627424478530884f * fC1_y_6;
    float pSH4_x_4 = 0.54627424478530884f * fS1_x_6;
    float v_x_12 = _S1146 + _S1139 * (pSH4_x_4 * *_S1068 + pSH8_x_10 * *_S1072 + fTmp0B_16 * *_S1071);
    float v_y_12 = _S1147 + _S1139 * (pSH8_x_10 * *_S1068 + pSH8_y_4 * *_S1072 + fTmp0B_16 * *_S1069);
    float v_z_12 = _S1148 + _S1139 * (pSH6_z_6 * *_S1070 + pSH7_z_4 * *_S1071 + pSH5_z_4 * *_S1069);
    float * _S1154 = v_coeffs_18 + int(24);
    *_S1154 = *_S1154 + pSH9_8 * _S1139;
    float * _S1155 = v_coeffs_18 + int(27);
    *_S1155 = *_S1155 + pSH10_8 * _S1139;
    float * _S1156 = v_coeffs_18 + int(30);
    *_S1156 = *_S1156 + pSH11_8 * _S1139;
    float * _S1157 = v_coeffs_18 + int(33);
    *_S1157 = *_S1157 + pSH12_9 * _S1139;
    float * _S1158 = v_coeffs_18 + int(36);
    *_S1158 = *_S1158 + pSH13_8 * _S1139;
    float * _S1159 = v_coeffs_18 + int(39);
    *_S1159 = *_S1159 + pSH14_8 * _S1139;
    float * _S1160 = v_coeffs_18 + int(42);
    *_S1160 = *_S1160 + pSH15_8 * _S1139;
    float fTmp0C_z_6 = -4.57045793533325195f * z_21;
    float _S1161 = x_24 * _S1067;
    float fC2_x_2 = fC1_10 + _S1161 - y_22 * fS1_x_6;
    float _S1162 = y_22 * _S1067;
    float fC2_y_2 = x_24 * fC1_y_6 - fS1_10 - _S1162;
    float fS2_x_2 = fS1_10 + x_24 * fS1_x_6 + _S1162;
    float fS2_y_2 = _S1161 + fC1_10 + y_22 * fC1_y_6;
    float pSH12_z_4 = 5.59764480590820312f * z2_10 - 1.11952900886535645f;
    float pSH13_z_2 = fTmp0C_z_6 * x_24;
    float pSH11_z_2 = fTmp0C_z_6 * y_22;
    float pSH14_x_6 = fTmp1B_10 * _S1067;
    float pSH14_y_2 = fTmp1B_10 * fC1_y_6;
    float pSH14_z_2 = 1.44530570507049561f * fC1_10;
    float pSH10_x_2 = fTmp1B_10 * fS1_x_6;
    float pSH10_z_2 = 1.44530570507049561f * fS1_10;
    float pSH15_x_2 = -0.59004360437393188f * fC2_x_2;
    float pSH15_y_2 = -0.59004360437393188f * fC2_y_2;
    float pSH9_x_2 = -0.59004360437393188f * fS2_x_2;
    float pSH9_y_2 = -0.59004360437393188f * fS2_y_2;
    float v_x_13 = v_x_12 + _S1139 * (pSH9_x_2 * *_S1073 + pSH15_x_2 * *_S1079 + pSH10_x_2 * *_S1074 + pSH14_x_6 * *_S1078 + fTmp0C_10 * *_S1077);
    float v_y_13 = v_y_12 + _S1139 * (pSH9_y_2 * *_S1073 + pSH15_y_2 * *_S1079 + pSH14_x_6 * *_S1074 + pSH14_y_2 * *_S1078 + fTmp0C_10 * *_S1075);
    float v_z_13 = v_z_12 + _S1139 * (pSH12_z_4 * *_S1076 + pSH13_z_2 * *_S1077 + pSH11_z_2 * *_S1075 + pSH14_z_2 * *_S1078 + pSH10_z_2 * *_S1074);
    float pSH20_4 = _S1080 + -1.00623059272766113f * pSH6_14;
    float * _S1163 = v_coeffs_18 + int(45);
    *_S1163 = *_S1163 + pSH16_3 * _S1139;
    float * _S1164 = v_coeffs_18 + int(48);
    *_S1164 = *_S1164 + pSH17_3 * _S1139;
    float * _S1165 = v_coeffs_18 + int(51);
    *_S1165 = *_S1165 + pSH18_3 * _S1139;
    float * _S1166 = v_coeffs_18 + int(54);
    *_S1166 = *_S1166 + pSH19_3 * _S1139;
    float * _S1167 = v_coeffs_18 + int(57);
    *_S1167 = *_S1167 + pSH20_4 * _S1139;
    float * _S1168 = v_coeffs_18 + int(60);
    *_S1168 = *_S1168 + pSH21_3 * _S1139;
    float * _S1169 = v_coeffs_18 + int(63);
    *_S1169 = *_S1169 + pSH22_3 * _S1139;
    float * _S1170 = v_coeffs_18 + int(66);
    *_S1170 = *_S1170 + pSH23_3 * _S1139;
    float * _S1171 = v_coeffs_18 + int(69);
    *_S1171 = *_S1171 + pSH24_3 * _S1139;
    float fTmp0D_z_2 = -14.04997730255126953f * z2_10 + 2.00713968276977539f;
    float fTmp1C_z_2 = 6.62322282791137695f * z_21;
    float pSH20_z_0 = 1.9843134880065918f * (pSH12_9 + z_21 * pSH12_z_4) + -1.00623059272766113f * pSH6_z_6;
    float pSH21_z_0 = fTmp0D_z_2 * x_24;
    float pSH19_z_0 = fTmp0D_z_2 * y_22;
    float pSH22_x_2 = fTmp1C_4 * _S1067;
    float pSH22_y_0 = fTmp1C_4 * fC1_y_6;
    float pSH22_z_0 = fTmp1C_z_2 * fC1_10;
    float pSH18_x_0 = fTmp1C_4 * fS1_x_6;
    float pSH18_z_0 = fTmp1C_z_2 * fS1_10;
    float pSH23_x_0 = fTmp2B_4 * fC2_x_2;
    float pSH23_y_0 = fTmp2B_4 * fC2_y_2;
    float pSH23_z_0 = -1.77013075351715088f * fC2_4;
    float pSH17_x_0 = fTmp2B_4 * fS2_x_2;
    float pSH17_y_0 = fTmp2B_4 * fS2_y_2;
    float pSH17_z_0 = -1.77013075351715088f * fS2_4;
    float pSH24_x_0 = 0.62583571672439575f * (fC2_4 + x_24 * fC2_x_2 - y_22 * fS2_x_2);
    float pSH24_y_0 = 0.62583571672439575f * (x_24 * fC2_y_2 - fS2_4 - y_22 * fS2_y_2);
    float pSH16_x_0 = 0.62583571672439575f * (fS2_4 + y_22 * fC2_x_2 + x_24 * fS2_x_2);
    float pSH16_y_0 = 0.62583571672439575f * (x_24 * fS2_y_2 + fC2_4 + y_22 * fC2_y_2);
    float3  dir_n_14 = make_float3 (x_24, y_22, z_21);
    float3  v_dir_n_26 = make_float3 (v_x_13 + _S1139 * (pSH16_x_0 * *_S1081 + pSH24_x_0 * *_S1089 + pSH17_x_0 * *_S1082 + pSH23_x_0 * *_S1088 + pSH18_x_0 * *_S1083 + pSH22_x_2 * *_S1087 + fTmp0D_4 * *_S1086), v_y_13 + _S1139 * (pSH16_y_0 * *_S1081 + pSH24_y_0 * *_S1089 + pSH17_y_0 * *_S1082 + pSH23_y_0 * *_S1088 + pSH22_x_2 * *_S1083 + pSH22_y_0 * *_S1087 + fTmp0D_4 * *_S1084), v_z_13 + _S1139 * (pSH20_z_0 * *_S1085 + pSH21_z_0 * *_S1086 + pSH19_z_0 * *_S1084 + pSH22_z_0 * *_S1087 + pSH18_z_0 * *_S1083 + pSH23_z_0 * *_S1088 + pSH17_z_0 * *_S1082));
    float3  v_viewdir_45 = v_viewdir_44 + (v_dir_n_26 - make_float3 (dot_0(v_dir_n_26, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    float _S1172 = _S1138.y;
    *&(v_coeff_dc_18->y) = *&(v_coeff_dc_18->y) + 0.282094806432724f * _S1172;
    float * _S1173 = v_coeffs_18 + int(1);
    *_S1173 = *_S1173 + _S1141 * _S1172;
    float * _S1174 = v_coeffs_18 + int(4);
    *_S1174 = *_S1174 + _S1143 * _S1172;
    float * _S1175 = v_coeffs_18 + int(7);
    *_S1175 = *_S1175 + _S1145 * _S1172;
    float _S1176 = -0.48860251903533936f * *_S1092 * _S1172;
    float _S1177 = -0.48860251903533936f * *_S1090 * _S1172;
    float _S1178 = 0.48860251903533936f * *_S1091 * _S1172;
    float * _S1179 = v_coeffs_18 + int(10);
    *_S1179 = *_S1179 + pSH4_13 * _S1172;
    float * _S1180 = v_coeffs_18 + int(13);
    *_S1180 = *_S1180 + pSH5_13 * _S1172;
    float * _S1181 = v_coeffs_18 + int(16);
    *_S1181 = *_S1181 + pSH6_14 * _S1172;
    float * _S1182 = v_coeffs_18 + int(19);
    *_S1182 = *_S1182 + pSH7_13 * _S1172;
    float * _S1183 = v_coeffs_18 + int(22);
    *_S1183 = *_S1183 + pSH8_13 * _S1172;
    float v_x_14 = _S1176 + _S1172 * (pSH4_x_4 * *_S1093 + pSH8_x_10 * *_S1097 + fTmp0B_16 * *_S1096);
    float v_y_14 = _S1177 + _S1172 * (pSH8_x_10 * *_S1093 + pSH8_y_4 * *_S1097 + fTmp0B_16 * *_S1094);
    float v_z_14 = _S1178 + _S1172 * (pSH6_z_6 * *_S1095 + pSH7_z_4 * *_S1096 + pSH5_z_4 * *_S1094);
    float * _S1184 = v_coeffs_18 + int(25);
    *_S1184 = *_S1184 + pSH9_8 * _S1172;
    float * _S1185 = v_coeffs_18 + int(28);
    *_S1185 = *_S1185 + pSH10_8 * _S1172;
    float * _S1186 = v_coeffs_18 + int(31);
    *_S1186 = *_S1186 + pSH11_8 * _S1172;
    float * _S1187 = v_coeffs_18 + int(34);
    *_S1187 = *_S1187 + pSH12_9 * _S1172;
    float * _S1188 = v_coeffs_18 + int(37);
    *_S1188 = *_S1188 + pSH13_8 * _S1172;
    float * _S1189 = v_coeffs_18 + int(40);
    *_S1189 = *_S1189 + pSH14_8 * _S1172;
    float * _S1190 = v_coeffs_18 + int(43);
    *_S1190 = *_S1190 + pSH15_8 * _S1172;
    float v_x_15 = v_x_14 + _S1172 * (pSH9_x_2 * *_S1098 + pSH15_x_2 * *_S1104 + pSH10_x_2 * *_S1099 + pSH14_x_6 * *_S1103 + fTmp0C_10 * *_S1102);
    float v_y_15 = v_y_14 + _S1172 * (pSH9_y_2 * *_S1098 + pSH15_y_2 * *_S1104 + pSH14_x_6 * *_S1099 + pSH14_y_2 * *_S1103 + fTmp0C_10 * *_S1100);
    float v_z_15 = v_z_14 + _S1172 * (pSH12_z_4 * *_S1101 + pSH13_z_2 * *_S1102 + pSH11_z_2 * *_S1100 + pSH14_z_2 * *_S1103 + pSH10_z_2 * *_S1099);
    float * _S1191 = v_coeffs_18 + int(46);
    *_S1191 = *_S1191 + pSH16_3 * _S1172;
    float * _S1192 = v_coeffs_18 + int(49);
    *_S1192 = *_S1192 + pSH17_3 * _S1172;
    float * _S1193 = v_coeffs_18 + int(52);
    *_S1193 = *_S1193 + pSH18_3 * _S1172;
    float * _S1194 = v_coeffs_18 + int(55);
    *_S1194 = *_S1194 + pSH19_3 * _S1172;
    float * _S1195 = v_coeffs_18 + int(58);
    *_S1195 = *_S1195 + pSH20_4 * _S1172;
    float * _S1196 = v_coeffs_18 + int(61);
    *_S1196 = *_S1196 + pSH21_3 * _S1172;
    float * _S1197 = v_coeffs_18 + int(64);
    *_S1197 = *_S1197 + pSH22_3 * _S1172;
    float * _S1198 = v_coeffs_18 + int(67);
    *_S1198 = *_S1198 + pSH23_3 * _S1172;
    float * _S1199 = v_coeffs_18 + int(70);
    *_S1199 = *_S1199 + pSH24_3 * _S1172;
    float3  v_dir_n_27 = make_float3 (v_x_15 + _S1172 * (pSH16_x_0 * *_S1105 + pSH24_x_0 * *_S1113 + pSH17_x_0 * *_S1106 + pSH23_x_0 * *_S1112 + pSH18_x_0 * *_S1107 + pSH22_x_2 * *_S1111 + fTmp0D_4 * *_S1110), v_y_15 + _S1172 * (pSH16_y_0 * *_S1105 + pSH24_y_0 * *_S1113 + pSH17_y_0 * *_S1106 + pSH23_y_0 * *_S1112 + pSH22_x_2 * *_S1107 + pSH22_y_0 * *_S1111 + fTmp0D_4 * *_S1108), v_z_15 + _S1172 * (pSH20_z_0 * *_S1109 + pSH21_z_0 * *_S1110 + pSH19_z_0 * *_S1108 + pSH22_z_0 * *_S1111 + pSH18_z_0 * *_S1107 + pSH23_z_0 * *_S1112 + pSH17_z_0 * *_S1106));
    float3  v_viewdir_46 = v_viewdir_45 + (v_dir_n_27 - make_float3 (dot_0(v_dir_n_27, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    float _S1200 = _S1138.z;
    *&(v_coeff_dc_18->z) = *&(v_coeff_dc_18->z) + 0.282094806432724f * _S1200;
    float * _S1201 = v_coeffs_18 + int(2);
    *_S1201 = *_S1201 + _S1141 * _S1200;
    float * _S1202 = v_coeffs_18 + int(5);
    *_S1202 = *_S1202 + _S1143 * _S1200;
    float * _S1203 = v_coeffs_18 + int(8);
    *_S1203 = *_S1203 + _S1145 * _S1200;
    float _S1204 = -0.48860251903533936f * *_S1116 * _S1200;
    float _S1205 = -0.48860251903533936f * *_S1114 * _S1200;
    float _S1206 = 0.48860251903533936f * *_S1115 * _S1200;
    float * _S1207 = v_coeffs_18 + int(11);
    *_S1207 = *_S1207 + pSH4_13 * _S1200;
    float * _S1208 = v_coeffs_18 + int(14);
    *_S1208 = *_S1208 + pSH5_13 * _S1200;
    float * _S1209 = v_coeffs_18 + int(17);
    *_S1209 = *_S1209 + pSH6_14 * _S1200;
    float * _S1210 = v_coeffs_18 + int(20);
    *_S1210 = *_S1210 + pSH7_13 * _S1200;
    float * _S1211 = v_coeffs_18 + int(23);
    *_S1211 = *_S1211 + pSH8_13 * _S1200;
    float v_x_16 = _S1204 + _S1200 * (pSH4_x_4 * *_S1117 + pSH8_x_10 * *_S1121 + fTmp0B_16 * *_S1120);
    float v_y_16 = _S1205 + _S1200 * (pSH8_x_10 * *_S1117 + pSH8_y_4 * *_S1121 + fTmp0B_16 * *_S1118);
    float v_z_16 = _S1206 + _S1200 * (pSH6_z_6 * *_S1119 + pSH7_z_4 * *_S1120 + pSH5_z_4 * *_S1118);
    float * _S1212 = v_coeffs_18 + int(26);
    *_S1212 = *_S1212 + pSH9_8 * _S1200;
    float * _S1213 = v_coeffs_18 + int(29);
    *_S1213 = *_S1213 + pSH10_8 * _S1200;
    float * _S1214 = v_coeffs_18 + int(32);
    *_S1214 = *_S1214 + pSH11_8 * _S1200;
    float * _S1215 = v_coeffs_18 + int(35);
    *_S1215 = *_S1215 + pSH12_9 * _S1200;
    float * _S1216 = v_coeffs_18 + int(38);
    *_S1216 = *_S1216 + pSH13_8 * _S1200;
    float * _S1217 = v_coeffs_18 + int(41);
    *_S1217 = *_S1217 + pSH14_8 * _S1200;
    float * _S1218 = v_coeffs_18 + int(44);
    *_S1218 = *_S1218 + pSH15_8 * _S1200;
    float v_x_17 = v_x_16 + _S1200 * (pSH9_x_2 * *_S1122 + pSH15_x_2 * *_S1128 + pSH10_x_2 * *_S1123 + pSH14_x_6 * *_S1127 + fTmp0C_10 * *_S1126);
    float v_y_17 = v_y_16 + _S1200 * (pSH9_y_2 * *_S1122 + pSH15_y_2 * *_S1128 + pSH14_x_6 * *_S1123 + pSH14_y_2 * *_S1127 + fTmp0C_10 * *_S1124);
    float v_z_17 = v_z_16 + _S1200 * (pSH12_z_4 * *_S1125 + pSH13_z_2 * *_S1126 + pSH11_z_2 * *_S1124 + pSH14_z_2 * *_S1127 + pSH10_z_2 * *_S1123);
    float * _S1219 = v_coeffs_18 + int(47);
    *_S1219 = *_S1219 + pSH16_3 * _S1200;
    float * _S1220 = v_coeffs_18 + int(50);
    *_S1220 = *_S1220 + pSH17_3 * _S1200;
    float * _S1221 = v_coeffs_18 + int(53);
    *_S1221 = *_S1221 + pSH18_3 * _S1200;
    float * _S1222 = v_coeffs_18 + int(56);
    *_S1222 = *_S1222 + pSH19_3 * _S1200;
    float * _S1223 = v_coeffs_18 + int(59);
    *_S1223 = *_S1223 + pSH20_4 * _S1200;
    float * _S1224 = v_coeffs_18 + int(62);
    *_S1224 = *_S1224 + pSH21_3 * _S1200;
    float * _S1225 = v_coeffs_18 + int(65);
    *_S1225 = *_S1225 + pSH22_3 * _S1200;
    float * _S1226 = v_coeffs_18 + int(68);
    *_S1226 = *_S1226 + pSH23_3 * _S1200;
    float * _S1227 = v_coeffs_18 + int(71);
    *_S1227 = *_S1227 + pSH24_3 * _S1200;
    float3  v_dir_n_28 = make_float3 (v_x_17 + _S1200 * (pSH16_x_0 * *_S1129 + pSH24_x_0 * *_S1137 + pSH17_x_0 * *_S1130 + pSH23_x_0 * *_S1136 + pSH18_x_0 * *_S1131 + pSH22_x_2 * *_S1135 + fTmp0D_4 * *_S1134), v_y_17 + _S1200 * (pSH16_y_0 * *_S1129 + pSH24_y_0 * *_S1137 + pSH17_y_0 * *_S1130 + pSH23_y_0 * *_S1136 + pSH22_x_2 * *_S1131 + pSH22_y_0 * *_S1135 + fTmp0D_4 * *_S1132), v_z_17 + _S1200 * (pSH20_z_0 * *_S1133 + pSH21_z_0 * *_S1134 + pSH19_z_0 * *_S1132 + pSH22_z_0 * *_S1135 + pSH18_z_0 * *_S1131 + pSH23_z_0 * *_S1136 + pSH17_z_0 * *_S1130));
    float3  v_viewdir_47 = v_viewdir_46 + (v_dir_n_28 - make_float3 (dot_0(v_dir_n_28, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    Matrix<float, 3, 3>  _S1228 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1229;
    (&_S1229)->primal_0 = _S1058;
    (&_S1229)->differential_0 = _S1228;
    float3  _S1230 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1231;
    (&_S1231)->primal_0 = t_28;
    (&_S1231)->differential_0 = _S1230;
    s_bwd_prop_mul_0(&_S1229, &_S1231, v_viewdir_47);
    Matrix<float, 3, 3>  _S1232 = transpose_0(_S1229.differential_0);
    *v_mean_18 = *v_mean_18 + v_viewdir_47;
    *v_R_18 = *v_R_18 + _S1232;
    *v_t_18 = *v_t_18 + _S1231.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_atomic(float3  mean_29, Matrix<float, 3, 3>  R_29, float3  t_29, float3  coeff_dc_29, float * coeffs_29, float3  v_colors_19, float3  * v_coeff_dc_19, float * v_coeffs_19, float3  * v_mean_19, Matrix<float, 3, 3>  * v_R_19, float3  * v_t_19)
{
    Matrix<float, 3, 3>  _S1233 = transpose_0(R_29);
    float3  _S1234 = mean_29 + mul_0(_S1233, t_29);
    float _S1235 = _S1234.x;
    float _S1236 = _S1234.y;
    float _S1237 = _S1234.z;
    float inv_norm_23 = (F32_rsqrt((_S1235 * _S1235 + _S1236 * _S1236 + _S1237 * _S1237)));
    float x_25 = _S1235 * inv_norm_23;
    float y_23 = _S1236 * inv_norm_23;
    float z_22 = _S1237 * inv_norm_23;
    float _S1238 = - y_23;
    float * _S1239 = coeffs_29 + int(0);
    float * _S1240 = coeffs_29 + int(3);
    float * _S1241 = coeffs_29 + int(6);
    float z2_11 = z_22 * z_22;
    float fTmp0B_17 = -1.09254848957061768f * z_22;
    float fC1_11 = x_25 * x_25 - y_23 * y_23;
    float _S1242 = 2.0f * x_25;
    float fS1_11 = _S1242 * y_23;
    float pSH6_15 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float pSH7_14 = fTmp0B_17 * x_25;
    float pSH5_14 = fTmp0B_17 * y_23;
    float pSH8_14 = 0.54627424478530884f * fC1_11;
    float pSH4_14 = 0.54627424478530884f * fS1_11;
    float * _S1243 = coeffs_29 + int(9);
    float * _S1244 = coeffs_29 + int(12);
    float * _S1245 = coeffs_29 + int(15);
    float * _S1246 = coeffs_29 + int(18);
    float * _S1247 = coeffs_29 + int(21);
    float fTmp0C_11 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_22;
    float fC2_5 = x_25 * fC1_11 - y_23 * fS1_11;
    float fS2_5 = x_25 * fS1_11 + y_23 * fC1_11;
    float pSH12_10 = z_22 * (1.86588168144226074f * z2_11 - 1.11952900886535645f);
    float pSH13_9 = fTmp0C_11 * x_25;
    float pSH11_9 = fTmp0C_11 * y_23;
    float pSH14_9 = fTmp1B_11 * fC1_11;
    float pSH10_9 = fTmp1B_11 * fS1_11;
    float pSH15_9 = -0.59004360437393188f * fC2_5;
    float pSH9_9 = -0.59004360437393188f * fS2_5;
    float * _S1248 = coeffs_29 + int(24);
    float * _S1249 = coeffs_29 + int(27);
    float * _S1250 = coeffs_29 + int(30);
    float * _S1251 = coeffs_29 + int(33);
    float * _S1252 = coeffs_29 + int(36);
    float * _S1253 = coeffs_29 + int(39);
    float * _S1254 = coeffs_29 + int(42);
    float fTmp0D_5 = z_22 * (-4.68332576751708984f * z2_11 + 2.00713968276977539f);
    float fTmp1C_5 = 3.31161141395568848f * z2_11 - 0.47308734059333801f;
    float fTmp2B_5 = -1.77013075351715088f * z_22;
    float _S1255 = 1.9843134880065918f * z_22 * pSH12_10;
    float pSH20_5 = _S1255 - 1.00623059272766113f * pSH6_15;
    float pSH21_4 = fTmp0D_5 * x_25;
    float pSH19_4 = fTmp0D_5 * y_23;
    float pSH22_4 = fTmp1C_5 * fC1_11;
    float pSH18_4 = fTmp1C_5 * fS1_11;
    float pSH23_4 = fTmp2B_5 * fC2_5;
    float pSH17_4 = fTmp2B_5 * fS2_5;
    float pSH24_4 = 0.62583571672439575f * (x_25 * fC2_5 - y_23 * fS2_5);
    float pSH16_4 = 0.62583571672439575f * (x_25 * fS2_5 + y_23 * fC2_5);
    float * _S1256 = coeffs_29 + int(45);
    float * _S1257 = coeffs_29 + int(48);
    float * _S1258 = coeffs_29 + int(51);
    float * _S1259 = coeffs_29 + int(54);
    float * _S1260 = coeffs_29 + int(57);
    float * _S1261 = coeffs_29 + int(60);
    float * _S1262 = coeffs_29 + int(63);
    float * _S1263 = coeffs_29 + int(66);
    float * _S1264 = coeffs_29 + int(69);
    float * _S1265 = coeffs_29 + int(1);
    float * _S1266 = coeffs_29 + int(4);
    float * _S1267 = coeffs_29 + int(7);
    float * _S1268 = coeffs_29 + int(10);
    float * _S1269 = coeffs_29 + int(13);
    float * _S1270 = coeffs_29 + int(16);
    float * _S1271 = coeffs_29 + int(19);
    float * _S1272 = coeffs_29 + int(22);
    float * _S1273 = coeffs_29 + int(25);
    float * _S1274 = coeffs_29 + int(28);
    float * _S1275 = coeffs_29 + int(31);
    float * _S1276 = coeffs_29 + int(34);
    float * _S1277 = coeffs_29 + int(37);
    float * _S1278 = coeffs_29 + int(40);
    float * _S1279 = coeffs_29 + int(43);
    float * _S1280 = coeffs_29 + int(46);
    float * _S1281 = coeffs_29 + int(49);
    float * _S1282 = coeffs_29 + int(52);
    float * _S1283 = coeffs_29 + int(55);
    float * _S1284 = coeffs_29 + int(58);
    float * _S1285 = coeffs_29 + int(61);
    float * _S1286 = coeffs_29 + int(64);
    float * _S1287 = coeffs_29 + int(67);
    float * _S1288 = coeffs_29 + int(70);
    float * _S1289 = coeffs_29 + int(2);
    float * _S1290 = coeffs_29 + int(5);
    float * _S1291 = coeffs_29 + int(8);
    float * _S1292 = coeffs_29 + int(11);
    float * _S1293 = coeffs_29 + int(14);
    float * _S1294 = coeffs_29 + int(17);
    float * _S1295 = coeffs_29 + int(20);
    float * _S1296 = coeffs_29 + int(23);
    float * _S1297 = coeffs_29 + int(26);
    float * _S1298 = coeffs_29 + int(29);
    float * _S1299 = coeffs_29 + int(32);
    float * _S1300 = coeffs_29 + int(35);
    float * _S1301 = coeffs_29 + int(38);
    float * _S1302 = coeffs_29 + int(41);
    float * _S1303 = coeffs_29 + int(44);
    float * _S1304 = coeffs_29 + int(47);
    float * _S1305 = coeffs_29 + int(50);
    float * _S1306 = coeffs_29 + int(53);
    float * _S1307 = coeffs_29 + int(56);
    float * _S1308 = coeffs_29 + int(59);
    float * _S1309 = coeffs_29 + int(62);
    float * _S1310 = coeffs_29 + int(65);
    float * _S1311 = coeffs_29 + int(68);
    float * _S1312 = coeffs_29 + int(71);
    float3  _S1313 = v_colors_19 * make_float3 (float((0.282094806432724f * coeff_dc_29.x + 0.48860251903533936f * (_S1238 * *_S1239 + z_22 * *_S1240 - x_25 * *_S1241) + (pSH4_14 * *_S1243 + pSH5_14 * *_S1244 + pSH6_15 * *_S1245 + pSH7_14 * *_S1246 + pSH8_14 * *_S1247) + (pSH9_9 * *_S1248 + pSH10_9 * *_S1249 + pSH11_9 * *_S1250 + pSH12_10 * *_S1251 + pSH13_9 * *_S1252 + pSH14_9 * *_S1253 + pSH15_9 * *_S1254) + (pSH16_4 * *_S1256 + pSH17_4 * *_S1257 + pSH18_4 * *_S1258 + pSH19_4 * *_S1259 + pSH20_5 * *_S1260 + pSH21_4 * *_S1261 + pSH22_4 * *_S1262 + pSH23_4 * *_S1263 + pSH24_4 * *_S1264)) >= -0.5f), float((0.282094806432724f * coeff_dc_29.y + 0.48860251903533936f * (_S1238 * *_S1265 + z_22 * *_S1266 - x_25 * *_S1267) + (pSH4_14 * *_S1268 + pSH5_14 * *_S1269 + pSH6_15 * *_S1270 + pSH7_14 * *_S1271 + pSH8_14 * *_S1272) + (pSH9_9 * *_S1273 + pSH10_9 * *_S1274 + pSH11_9 * *_S1275 + pSH12_10 * *_S1276 + pSH13_9 * *_S1277 + pSH14_9 * *_S1278 + pSH15_9 * *_S1279) + (pSH16_4 * *_S1280 + pSH17_4 * *_S1281 + pSH18_4 * *_S1282 + pSH19_4 * *_S1283 + pSH20_5 * *_S1284 + pSH21_4 * *_S1285 + pSH22_4 * *_S1286 + pSH23_4 * *_S1287 + pSH24_4 * *_S1288)) >= -0.5f), float((0.282094806432724f * coeff_dc_29.z + 0.48860251903533936f * (_S1238 * *_S1289 + z_22 * *_S1290 - x_25 * *_S1291) + (pSH4_14 * *_S1292 + pSH5_14 * *_S1293 + pSH6_15 * *_S1294 + pSH7_14 * *_S1295 + pSH8_14 * *_S1296) + (pSH9_9 * *_S1297 + pSH10_9 * *_S1298 + pSH11_9 * *_S1299 + pSH12_10 * *_S1300 + pSH13_9 * *_S1301 + pSH14_9 * *_S1302 + pSH15_9 * *_S1303) + (pSH16_4 * *_S1304 + pSH17_4 * *_S1305 + pSH18_4 * *_S1306 + pSH19_4 * *_S1307 + pSH20_5 * *_S1308 + pSH21_4 * *_S1309 + pSH22_4 * *_S1310 + pSH23_4 * *_S1311 + pSH24_4 * *_S1312)) >= -0.5f));
    float3  v_viewdir_48 = {};
    float _S1314 = _S1313.x;
    *&(v_coeff_dc_19->x) = *&(v_coeff_dc_19->x) + 0.282094806432724f * _S1314;
    float _S1315 = -0.48860251903533936f * y_23;
    float temp_128 = _S1315 * _S1314;
    bool _S1316;
    if((F32_isfinite((temp_128))))
    {
        _S1316 = temp_128 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1317 = atomicAdd(v_coeffs_19 + int(0), temp_128);
    }
    float _S1318 = 0.48860251903533936f * z_22;
    float temp_129 = _S1318 * _S1314;
    if((F32_isfinite((temp_129))))
    {
        _S1316 = temp_129 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1319 = atomicAdd(v_coeffs_19 + int(3), temp_129);
    }
    float _S1320 = -0.48860251903533936f * x_25;
    float temp_130 = _S1320 * _S1314;
    if((F32_isfinite((temp_130))))
    {
        _S1316 = temp_130 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1321 = atomicAdd(v_coeffs_19 + int(6), temp_130);
    }
    float _S1322 = -0.48860251903533936f * *_S1241 * _S1314;
    float _S1323 = -0.48860251903533936f * *_S1239 * _S1314;
    float _S1324 = 0.48860251903533936f * *_S1240 * _S1314;
    float temp_131 = pSH4_14 * _S1314;
    if((F32_isfinite((temp_131))))
    {
        _S1316 = temp_131 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1325 = atomicAdd(v_coeffs_19 + int(9), temp_131);
    }
    float temp_132 = pSH5_14 * _S1314;
    if((F32_isfinite((temp_132))))
    {
        _S1316 = temp_132 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1326 = atomicAdd(v_coeffs_19 + int(12), temp_132);
    }
    float temp_133 = pSH6_15 * _S1314;
    if((F32_isfinite((temp_133))))
    {
        _S1316 = temp_133 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1327 = atomicAdd(v_coeffs_19 + int(15), temp_133);
    }
    float temp_134 = pSH7_14 * _S1314;
    if((F32_isfinite((temp_134))))
    {
        _S1316 = temp_134 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1328 = atomicAdd(v_coeffs_19 + int(18), temp_134);
    }
    float temp_135 = pSH8_14 * _S1314;
    if((F32_isfinite((temp_135))))
    {
        _S1316 = temp_135 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1329 = atomicAdd(v_coeffs_19 + int(21), temp_135);
    }
    float fC1_y_7 = -2.0f * y_23;
    float fS1_x_7 = 2.0f * y_23;
    float pSH6_z_7 = 1.89234936237335205f * z_22;
    float pSH7_z_5 = -1.09254848957061768f * x_25;
    float pSH5_z_5 = -1.09254848957061768f * y_23;
    float pSH8_x_11 = 0.54627424478530884f * _S1242;
    float pSH8_y_5 = 0.54627424478530884f * fC1_y_7;
    float pSH4_x_5 = 0.54627424478530884f * fS1_x_7;
    float v_x_18 = _S1322 + _S1314 * (pSH4_x_5 * *_S1243 + pSH8_x_11 * *_S1247 + fTmp0B_17 * *_S1246);
    float v_y_18 = _S1323 + _S1314 * (pSH8_x_11 * *_S1243 + pSH8_y_5 * *_S1247 + fTmp0B_17 * *_S1244);
    float v_z_18 = _S1324 + _S1314 * (pSH6_z_7 * *_S1245 + pSH7_z_5 * *_S1246 + pSH5_z_5 * *_S1244);
    float temp_136 = pSH9_9 * _S1314;
    if((F32_isfinite((temp_136))))
    {
        _S1316 = temp_136 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1330 = atomicAdd(v_coeffs_19 + int(24), temp_136);
    }
    float temp_137 = pSH10_9 * _S1314;
    if((F32_isfinite((temp_137))))
    {
        _S1316 = temp_137 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1331 = atomicAdd(v_coeffs_19 + int(27), temp_137);
    }
    float temp_138 = pSH11_9 * _S1314;
    if((F32_isfinite((temp_138))))
    {
        _S1316 = temp_138 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1332 = atomicAdd(v_coeffs_19 + int(30), temp_138);
    }
    float temp_139 = pSH12_10 * _S1314;
    if((F32_isfinite((temp_139))))
    {
        _S1316 = temp_139 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1333 = atomicAdd(v_coeffs_19 + int(33), temp_139);
    }
    float temp_140 = pSH13_9 * _S1314;
    if((F32_isfinite((temp_140))))
    {
        _S1316 = temp_140 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1334 = atomicAdd(v_coeffs_19 + int(36), temp_140);
    }
    float temp_141 = pSH14_9 * _S1314;
    if((F32_isfinite((temp_141))))
    {
        _S1316 = temp_141 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1335 = atomicAdd(v_coeffs_19 + int(39), temp_141);
    }
    float temp_142 = pSH15_9 * _S1314;
    if((F32_isfinite((temp_142))))
    {
        _S1316 = temp_142 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1336 = atomicAdd(v_coeffs_19 + int(42), temp_142);
    }
    float fTmp0C_z_7 = -4.57045793533325195f * z_22;
    float _S1337 = x_25 * _S1242;
    float fC2_x_3 = fC1_11 + _S1337 - y_23 * fS1_x_7;
    float _S1338 = y_23 * _S1242;
    float fC2_y_3 = x_25 * fC1_y_7 - fS1_11 - _S1338;
    float fS2_x_3 = fS1_11 + x_25 * fS1_x_7 + _S1338;
    float fS2_y_3 = _S1337 + fC1_11 + y_23 * fC1_y_7;
    float pSH12_z_5 = 5.59764480590820312f * z2_11 - 1.11952900886535645f;
    float pSH13_z_3 = fTmp0C_z_7 * x_25;
    float pSH11_z_3 = fTmp0C_z_7 * y_23;
    float pSH14_x_7 = fTmp1B_11 * _S1242;
    float pSH14_y_3 = fTmp1B_11 * fC1_y_7;
    float pSH14_z_3 = 1.44530570507049561f * fC1_11;
    float pSH10_x_3 = fTmp1B_11 * fS1_x_7;
    float pSH10_z_3 = 1.44530570507049561f * fS1_11;
    float pSH15_x_3 = -0.59004360437393188f * fC2_x_3;
    float pSH15_y_3 = -0.59004360437393188f * fC2_y_3;
    float pSH9_x_3 = -0.59004360437393188f * fS2_x_3;
    float pSH9_y_3 = -0.59004360437393188f * fS2_y_3;
    float v_x_19 = v_x_18 + _S1314 * (pSH9_x_3 * *_S1248 + pSH15_x_3 * *_S1254 + pSH10_x_3 * *_S1249 + pSH14_x_7 * *_S1253 + fTmp0C_11 * *_S1252);
    float v_y_19 = v_y_18 + _S1314 * (pSH9_y_3 * *_S1248 + pSH15_y_3 * *_S1254 + pSH14_x_7 * *_S1249 + pSH14_y_3 * *_S1253 + fTmp0C_11 * *_S1250);
    float v_z_19 = v_z_18 + _S1314 * (pSH12_z_5 * *_S1251 + pSH13_z_3 * *_S1252 + pSH11_z_3 * *_S1250 + pSH14_z_3 * *_S1253 + pSH10_z_3 * *_S1249);
    float pSH20_6 = _S1255 + -1.00623059272766113f * pSH6_15;
    float temp_143 = pSH16_4 * _S1314;
    if((F32_isfinite((temp_143))))
    {
        _S1316 = temp_143 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1339 = atomicAdd(v_coeffs_19 + int(45), temp_143);
    }
    float temp_144 = pSH17_4 * _S1314;
    if((F32_isfinite((temp_144))))
    {
        _S1316 = temp_144 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1340 = atomicAdd(v_coeffs_19 + int(48), temp_144);
    }
    float temp_145 = pSH18_4 * _S1314;
    if((F32_isfinite((temp_145))))
    {
        _S1316 = temp_145 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1341 = atomicAdd(v_coeffs_19 + int(51), temp_145);
    }
    float temp_146 = pSH19_4 * _S1314;
    if((F32_isfinite((temp_146))))
    {
        _S1316 = temp_146 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1342 = atomicAdd(v_coeffs_19 + int(54), temp_146);
    }
    float temp_147 = pSH20_6 * _S1314;
    if((F32_isfinite((temp_147))))
    {
        _S1316 = temp_147 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1343 = atomicAdd(v_coeffs_19 + int(57), temp_147);
    }
    float temp_148 = pSH21_4 * _S1314;
    if((F32_isfinite((temp_148))))
    {
        _S1316 = temp_148 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1344 = atomicAdd(v_coeffs_19 + int(60), temp_148);
    }
    float temp_149 = pSH22_4 * _S1314;
    if((F32_isfinite((temp_149))))
    {
        _S1316 = temp_149 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1345 = atomicAdd(v_coeffs_19 + int(63), temp_149);
    }
    float temp_150 = pSH23_4 * _S1314;
    if((F32_isfinite((temp_150))))
    {
        _S1316 = temp_150 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1346 = atomicAdd(v_coeffs_19 + int(66), temp_150);
    }
    float temp_151 = pSH24_4 * _S1314;
    if((F32_isfinite((temp_151))))
    {
        _S1316 = temp_151 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1347 = atomicAdd(v_coeffs_19 + int(69), temp_151);
    }
    float fTmp0D_z_3 = -14.04997730255126953f * z2_11 + 2.00713968276977539f;
    float fTmp1C_z_3 = 6.62322282791137695f * z_22;
    float pSH20_z_1 = 1.9843134880065918f * (pSH12_10 + z_22 * pSH12_z_5) + -1.00623059272766113f * pSH6_z_7;
    float pSH21_z_1 = fTmp0D_z_3 * x_25;
    float pSH19_z_1 = fTmp0D_z_3 * y_23;
    float pSH22_x_3 = fTmp1C_5 * _S1242;
    float pSH22_y_1 = fTmp1C_5 * fC1_y_7;
    float pSH22_z_1 = fTmp1C_z_3 * fC1_11;
    float pSH18_x_1 = fTmp1C_5 * fS1_x_7;
    float pSH18_z_1 = fTmp1C_z_3 * fS1_11;
    float pSH23_x_1 = fTmp2B_5 * fC2_x_3;
    float pSH23_y_1 = fTmp2B_5 * fC2_y_3;
    float pSH23_z_1 = -1.77013075351715088f * fC2_5;
    float pSH17_x_1 = fTmp2B_5 * fS2_x_3;
    float pSH17_y_1 = fTmp2B_5 * fS2_y_3;
    float pSH17_z_1 = -1.77013075351715088f * fS2_5;
    float pSH24_x_1 = 0.62583571672439575f * (fC2_5 + x_25 * fC2_x_3 - y_23 * fS2_x_3);
    float pSH24_y_1 = 0.62583571672439575f * (x_25 * fC2_y_3 - fS2_5 - y_23 * fS2_y_3);
    float pSH16_x_1 = 0.62583571672439575f * (fS2_5 + y_23 * fC2_x_3 + x_25 * fS2_x_3);
    float pSH16_y_1 = 0.62583571672439575f * (x_25 * fS2_y_3 + fC2_5 + y_23 * fC2_y_3);
    float3  dir_n_15 = make_float3 (x_25, y_23, z_22);
    float3  v_dir_n_29 = make_float3 (v_x_19 + _S1314 * (pSH16_x_1 * *_S1256 + pSH24_x_1 * *_S1264 + pSH17_x_1 * *_S1257 + pSH23_x_1 * *_S1263 + pSH18_x_1 * *_S1258 + pSH22_x_3 * *_S1262 + fTmp0D_5 * *_S1261), v_y_19 + _S1314 * (pSH16_y_1 * *_S1256 + pSH24_y_1 * *_S1264 + pSH17_y_1 * *_S1257 + pSH23_y_1 * *_S1263 + pSH22_x_3 * *_S1258 + pSH22_y_1 * *_S1262 + fTmp0D_5 * *_S1259), v_z_19 + _S1314 * (pSH20_z_1 * *_S1260 + pSH21_z_1 * *_S1261 + pSH19_z_1 * *_S1259 + pSH22_z_1 * *_S1262 + pSH18_z_1 * *_S1258 + pSH23_z_1 * *_S1263 + pSH17_z_1 * *_S1257));
    float3  v_viewdir_49 = v_viewdir_48 + (v_dir_n_29 - make_float3 (dot_0(v_dir_n_29, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    float _S1348 = _S1313.y;
    *&(v_coeff_dc_19->y) = *&(v_coeff_dc_19->y) + 0.282094806432724f * _S1348;
    float temp_152 = _S1315 * _S1348;
    if((F32_isfinite((temp_152))))
    {
        _S1316 = temp_152 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1349 = atomicAdd(v_coeffs_19 + int(1), temp_152);
    }
    float temp_153 = _S1318 * _S1348;
    if((F32_isfinite((temp_153))))
    {
        _S1316 = temp_153 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1350 = atomicAdd(v_coeffs_19 + int(4), temp_153);
    }
    float temp_154 = _S1320 * _S1348;
    if((F32_isfinite((temp_154))))
    {
        _S1316 = temp_154 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1351 = atomicAdd(v_coeffs_19 + int(7), temp_154);
    }
    float _S1352 = -0.48860251903533936f * *_S1267 * _S1348;
    float _S1353 = -0.48860251903533936f * *_S1265 * _S1348;
    float _S1354 = 0.48860251903533936f * *_S1266 * _S1348;
    float temp_155 = pSH4_14 * _S1348;
    if((F32_isfinite((temp_155))))
    {
        _S1316 = temp_155 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1355 = atomicAdd(v_coeffs_19 + int(10), temp_155);
    }
    float temp_156 = pSH5_14 * _S1348;
    if((F32_isfinite((temp_156))))
    {
        _S1316 = temp_156 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1356 = atomicAdd(v_coeffs_19 + int(13), temp_156);
    }
    float temp_157 = pSH6_15 * _S1348;
    if((F32_isfinite((temp_157))))
    {
        _S1316 = temp_157 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1357 = atomicAdd(v_coeffs_19 + int(16), temp_157);
    }
    float temp_158 = pSH7_14 * _S1348;
    if((F32_isfinite((temp_158))))
    {
        _S1316 = temp_158 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1358 = atomicAdd(v_coeffs_19 + int(19), temp_158);
    }
    float temp_159 = pSH8_14 * _S1348;
    if((F32_isfinite((temp_159))))
    {
        _S1316 = temp_159 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1359 = atomicAdd(v_coeffs_19 + int(22), temp_159);
    }
    float v_x_20 = _S1352 + _S1348 * (pSH4_x_5 * *_S1268 + pSH8_x_11 * *_S1272 + fTmp0B_17 * *_S1271);
    float v_y_20 = _S1353 + _S1348 * (pSH8_x_11 * *_S1268 + pSH8_y_5 * *_S1272 + fTmp0B_17 * *_S1269);
    float v_z_20 = _S1354 + _S1348 * (pSH6_z_7 * *_S1270 + pSH7_z_5 * *_S1271 + pSH5_z_5 * *_S1269);
    float temp_160 = pSH9_9 * _S1348;
    if((F32_isfinite((temp_160))))
    {
        _S1316 = temp_160 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1360 = atomicAdd(v_coeffs_19 + int(25), temp_160);
    }
    float temp_161 = pSH10_9 * _S1348;
    if((F32_isfinite((temp_161))))
    {
        _S1316 = temp_161 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1361 = atomicAdd(v_coeffs_19 + int(28), temp_161);
    }
    float temp_162 = pSH11_9 * _S1348;
    if((F32_isfinite((temp_162))))
    {
        _S1316 = temp_162 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1362 = atomicAdd(v_coeffs_19 + int(31), temp_162);
    }
    float temp_163 = pSH12_10 * _S1348;
    if((F32_isfinite((temp_163))))
    {
        _S1316 = temp_163 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1363 = atomicAdd(v_coeffs_19 + int(34), temp_163);
    }
    float temp_164 = pSH13_9 * _S1348;
    if((F32_isfinite((temp_164))))
    {
        _S1316 = temp_164 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1364 = atomicAdd(v_coeffs_19 + int(37), temp_164);
    }
    float temp_165 = pSH14_9 * _S1348;
    if((F32_isfinite((temp_165))))
    {
        _S1316 = temp_165 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1365 = atomicAdd(v_coeffs_19 + int(40), temp_165);
    }
    float temp_166 = pSH15_9 * _S1348;
    if((F32_isfinite((temp_166))))
    {
        _S1316 = temp_166 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1366 = atomicAdd(v_coeffs_19 + int(43), temp_166);
    }
    float v_x_21 = v_x_20 + _S1348 * (pSH9_x_3 * *_S1273 + pSH15_x_3 * *_S1279 + pSH10_x_3 * *_S1274 + pSH14_x_7 * *_S1278 + fTmp0C_11 * *_S1277);
    float v_y_21 = v_y_20 + _S1348 * (pSH9_y_3 * *_S1273 + pSH15_y_3 * *_S1279 + pSH14_x_7 * *_S1274 + pSH14_y_3 * *_S1278 + fTmp0C_11 * *_S1275);
    float v_z_21 = v_z_20 + _S1348 * (pSH12_z_5 * *_S1276 + pSH13_z_3 * *_S1277 + pSH11_z_3 * *_S1275 + pSH14_z_3 * *_S1278 + pSH10_z_3 * *_S1274);
    float temp_167 = pSH16_4 * _S1348;
    if((F32_isfinite((temp_167))))
    {
        _S1316 = temp_167 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1367 = atomicAdd(v_coeffs_19 + int(46), temp_167);
    }
    float temp_168 = pSH17_4 * _S1348;
    if((F32_isfinite((temp_168))))
    {
        _S1316 = temp_168 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1368 = atomicAdd(v_coeffs_19 + int(49), temp_168);
    }
    float temp_169 = pSH18_4 * _S1348;
    if((F32_isfinite((temp_169))))
    {
        _S1316 = temp_169 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1369 = atomicAdd(v_coeffs_19 + int(52), temp_169);
    }
    float temp_170 = pSH19_4 * _S1348;
    if((F32_isfinite((temp_170))))
    {
        _S1316 = temp_170 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1370 = atomicAdd(v_coeffs_19 + int(55), temp_170);
    }
    float temp_171 = pSH20_6 * _S1348;
    if((F32_isfinite((temp_171))))
    {
        _S1316 = temp_171 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1371 = atomicAdd(v_coeffs_19 + int(58), temp_171);
    }
    float temp_172 = pSH21_4 * _S1348;
    if((F32_isfinite((temp_172))))
    {
        _S1316 = temp_172 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1372 = atomicAdd(v_coeffs_19 + int(61), temp_172);
    }
    float temp_173 = pSH22_4 * _S1348;
    if((F32_isfinite((temp_173))))
    {
        _S1316 = temp_173 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1373 = atomicAdd(v_coeffs_19 + int(64), temp_173);
    }
    float temp_174 = pSH23_4 * _S1348;
    if((F32_isfinite((temp_174))))
    {
        _S1316 = temp_174 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1374 = atomicAdd(v_coeffs_19 + int(67), temp_174);
    }
    float temp_175 = pSH24_4 * _S1348;
    if((F32_isfinite((temp_175))))
    {
        _S1316 = temp_175 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1375 = atomicAdd(v_coeffs_19 + int(70), temp_175);
    }
    float3  v_dir_n_30 = make_float3 (v_x_21 + _S1348 * (pSH16_x_1 * *_S1280 + pSH24_x_1 * *_S1288 + pSH17_x_1 * *_S1281 + pSH23_x_1 * *_S1287 + pSH18_x_1 * *_S1282 + pSH22_x_3 * *_S1286 + fTmp0D_5 * *_S1285), v_y_21 + _S1348 * (pSH16_y_1 * *_S1280 + pSH24_y_1 * *_S1288 + pSH17_y_1 * *_S1281 + pSH23_y_1 * *_S1287 + pSH22_x_3 * *_S1282 + pSH22_y_1 * *_S1286 + fTmp0D_5 * *_S1283), v_z_21 + _S1348 * (pSH20_z_1 * *_S1284 + pSH21_z_1 * *_S1285 + pSH19_z_1 * *_S1283 + pSH22_z_1 * *_S1286 + pSH18_z_1 * *_S1282 + pSH23_z_1 * *_S1287 + pSH17_z_1 * *_S1281));
    float3  v_viewdir_50 = v_viewdir_49 + (v_dir_n_30 - make_float3 (dot_0(v_dir_n_30, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    float _S1376 = _S1313.z;
    *&(v_coeff_dc_19->z) = *&(v_coeff_dc_19->z) + 0.282094806432724f * _S1376;
    float temp_176 = _S1315 * _S1376;
    if((F32_isfinite((temp_176))))
    {
        _S1316 = temp_176 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1377 = atomicAdd(v_coeffs_19 + int(2), temp_176);
    }
    float temp_177 = _S1318 * _S1376;
    if((F32_isfinite((temp_177))))
    {
        _S1316 = temp_177 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1378 = atomicAdd(v_coeffs_19 + int(5), temp_177);
    }
    float temp_178 = _S1320 * _S1376;
    if((F32_isfinite((temp_178))))
    {
        _S1316 = temp_178 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1379 = atomicAdd(v_coeffs_19 + int(8), temp_178);
    }
    float _S1380 = -0.48860251903533936f * *_S1291 * _S1376;
    float _S1381 = -0.48860251903533936f * *_S1289 * _S1376;
    float _S1382 = 0.48860251903533936f * *_S1290 * _S1376;
    float temp_179 = pSH4_14 * _S1376;
    if((F32_isfinite((temp_179))))
    {
        _S1316 = temp_179 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1383 = atomicAdd(v_coeffs_19 + int(11), temp_179);
    }
    float temp_180 = pSH5_14 * _S1376;
    if((F32_isfinite((temp_180))))
    {
        _S1316 = temp_180 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1384 = atomicAdd(v_coeffs_19 + int(14), temp_180);
    }
    float temp_181 = pSH6_15 * _S1376;
    if((F32_isfinite((temp_181))))
    {
        _S1316 = temp_181 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1385 = atomicAdd(v_coeffs_19 + int(17), temp_181);
    }
    float temp_182 = pSH7_14 * _S1376;
    if((F32_isfinite((temp_182))))
    {
        _S1316 = temp_182 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1386 = atomicAdd(v_coeffs_19 + int(20), temp_182);
    }
    float temp_183 = pSH8_14 * _S1376;
    if((F32_isfinite((temp_183))))
    {
        _S1316 = temp_183 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1387 = atomicAdd(v_coeffs_19 + int(23), temp_183);
    }
    float v_x_22 = _S1380 + _S1376 * (pSH4_x_5 * *_S1292 + pSH8_x_11 * *_S1296 + fTmp0B_17 * *_S1295);
    float v_y_22 = _S1381 + _S1376 * (pSH8_x_11 * *_S1292 + pSH8_y_5 * *_S1296 + fTmp0B_17 * *_S1293);
    float v_z_22 = _S1382 + _S1376 * (pSH6_z_7 * *_S1294 + pSH7_z_5 * *_S1295 + pSH5_z_5 * *_S1293);
    float temp_184 = pSH9_9 * _S1376;
    if((F32_isfinite((temp_184))))
    {
        _S1316 = temp_184 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1388 = atomicAdd(v_coeffs_19 + int(26), temp_184);
    }
    float temp_185 = pSH10_9 * _S1376;
    if((F32_isfinite((temp_185))))
    {
        _S1316 = temp_185 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1389 = atomicAdd(v_coeffs_19 + int(29), temp_185);
    }
    float temp_186 = pSH11_9 * _S1376;
    if((F32_isfinite((temp_186))))
    {
        _S1316 = temp_186 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1390 = atomicAdd(v_coeffs_19 + int(32), temp_186);
    }
    float temp_187 = pSH12_10 * _S1376;
    if((F32_isfinite((temp_187))))
    {
        _S1316 = temp_187 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1391 = atomicAdd(v_coeffs_19 + int(35), temp_187);
    }
    float temp_188 = pSH13_9 * _S1376;
    if((F32_isfinite((temp_188))))
    {
        _S1316 = temp_188 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1392 = atomicAdd(v_coeffs_19 + int(38), temp_188);
    }
    float temp_189 = pSH14_9 * _S1376;
    if((F32_isfinite((temp_189))))
    {
        _S1316 = temp_189 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1393 = atomicAdd(v_coeffs_19 + int(41), temp_189);
    }
    float temp_190 = pSH15_9 * _S1376;
    if((F32_isfinite((temp_190))))
    {
        _S1316 = temp_190 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1394 = atomicAdd(v_coeffs_19 + int(44), temp_190);
    }
    float v_x_23 = v_x_22 + _S1376 * (pSH9_x_3 * *_S1297 + pSH15_x_3 * *_S1303 + pSH10_x_3 * *_S1298 + pSH14_x_7 * *_S1302 + fTmp0C_11 * *_S1301);
    float v_y_23 = v_y_22 + _S1376 * (pSH9_y_3 * *_S1297 + pSH15_y_3 * *_S1303 + pSH14_x_7 * *_S1298 + pSH14_y_3 * *_S1302 + fTmp0C_11 * *_S1299);
    float v_z_23 = v_z_22 + _S1376 * (pSH12_z_5 * *_S1300 + pSH13_z_3 * *_S1301 + pSH11_z_3 * *_S1299 + pSH14_z_3 * *_S1302 + pSH10_z_3 * *_S1298);
    float temp_191 = pSH16_4 * _S1376;
    if((F32_isfinite((temp_191))))
    {
        _S1316 = temp_191 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1395 = atomicAdd(v_coeffs_19 + int(47), temp_191);
    }
    float temp_192 = pSH17_4 * _S1376;
    if((F32_isfinite((temp_192))))
    {
        _S1316 = temp_192 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1396 = atomicAdd(v_coeffs_19 + int(50), temp_192);
    }
    float temp_193 = pSH18_4 * _S1376;
    if((F32_isfinite((temp_193))))
    {
        _S1316 = temp_193 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1397 = atomicAdd(v_coeffs_19 + int(53), temp_193);
    }
    float temp_194 = pSH19_4 * _S1376;
    if((F32_isfinite((temp_194))))
    {
        _S1316 = temp_194 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1398 = atomicAdd(v_coeffs_19 + int(56), temp_194);
    }
    float temp_195 = pSH20_6 * _S1376;
    if((F32_isfinite((temp_195))))
    {
        _S1316 = temp_195 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1399 = atomicAdd(v_coeffs_19 + int(59), temp_195);
    }
    float temp_196 = pSH21_4 * _S1376;
    if((F32_isfinite((temp_196))))
    {
        _S1316 = temp_196 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1400 = atomicAdd(v_coeffs_19 + int(62), temp_196);
    }
    float temp_197 = pSH22_4 * _S1376;
    if((F32_isfinite((temp_197))))
    {
        _S1316 = temp_197 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1401 = atomicAdd(v_coeffs_19 + int(65), temp_197);
    }
    float temp_198 = pSH23_4 * _S1376;
    if((F32_isfinite((temp_198))))
    {
        _S1316 = temp_198 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1402 = atomicAdd(v_coeffs_19 + int(68), temp_198);
    }
    float temp_199 = pSH24_4 * _S1376;
    if((F32_isfinite((temp_199))))
    {
        _S1316 = temp_199 != 0.0f;
    }
    else
    {
        _S1316 = false;
    }
    if(_S1316)
    {
        float _S1403 = atomicAdd(v_coeffs_19 + int(71), temp_199);
    }
    float3  v_dir_n_31 = make_float3 (v_x_23 + _S1376 * (pSH16_x_1 * *_S1304 + pSH24_x_1 * *_S1312 + pSH17_x_1 * *_S1305 + pSH23_x_1 * *_S1311 + pSH18_x_1 * *_S1306 + pSH22_x_3 * *_S1310 + fTmp0D_5 * *_S1309), v_y_23 + _S1376 * (pSH16_y_1 * *_S1304 + pSH24_y_1 * *_S1312 + pSH17_y_1 * *_S1305 + pSH23_y_1 * *_S1311 + pSH22_x_3 * *_S1306 + pSH22_y_1 * *_S1310 + fTmp0D_5 * *_S1307), v_z_23 + _S1376 * (pSH20_z_1 * *_S1308 + pSH21_z_1 * *_S1309 + pSH19_z_1 * *_S1307 + pSH22_z_1 * *_S1310 + pSH18_z_1 * *_S1306 + pSH23_z_1 * *_S1311 + pSH17_z_1 * *_S1305));
    float3  v_viewdir_51 = v_viewdir_50 + (v_dir_n_31 - make_float3 (dot_0(v_dir_n_31, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    Matrix<float, 3, 3>  _S1404 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1405;
    (&_S1405)->primal_0 = _S1233;
    (&_S1405)->differential_0 = _S1404;
    float3  _S1406 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1407;
    (&_S1407)->primal_0 = t_29;
    (&_S1407)->differential_0 = _S1406;
    s_bwd_prop_mul_0(&_S1405, &_S1407, v_viewdir_51);
    Matrix<float, 3, 3>  _S1408 = transpose_0(_S1405.differential_0);
    *v_mean_19 = *v_mean_19 + v_viewdir_51;
    *v_R_19 = *v_R_19 + _S1408;
    *v_t_19 = *v_t_19 + _S1407.differential_0;
    return;
}

