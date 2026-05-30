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

inline __device__ uint WaveGetActiveMask_0(uint _S4)
{
    return _S4;
}

inline __device__ float WaveActiveSum_0(float expr_0, uint _S5)
{
    float _S6 = (_waveSum((make_uint4 (WaveGetActiveMask_0(_S5), 0U, 0U, 0U)).x, (expr_0)));
    return _S6;
}

inline __device__ __shared__ FixedArray<float, 48>  _sh_block_reduce_shared_0;

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S7, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S8, float3  _S9)
{
    _d_mul_0(_S7, _S8, _S9);
    return;
}

inline __device__ void sh0_to_color_vjp_inplace(float3  mean_1, Matrix<float, 3, 3>  R_1, float3  t_1, float3  coeff_dc_1, float3  * coeffs_1, float3  v_colors_0, float3  * v_coeff_dc_0, float3  * v_coeffs_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  colors_0 = make_float3 (0.282094806432724f) * coeff_dc_1;
    *v_coeff_dc_0 = *v_coeff_dc_0 + make_float3 (0.282094806432724f) * (v_colors_0 * make_float3 (float((colors_0.x) >= -0.5f), float((colors_0.y) >= -0.5f), float((colors_0.z) >= -0.5f)));
    float3  v_viewdir_0 = {};
    Matrix<float, 3, 3>  _S10 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S11;
    (&_S11)->primal_0 = transpose_0(R_1);
    (&_S11)->differential_0 = _S10;
    float3  _S12 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S13;
    (&_S13)->primal_0 = t_1;
    (&_S13)->differential_0 = _S12;
    s_bwd_prop_mul_0(&_S11, &_S13, v_viewdir_0);
    Matrix<float, 3, 3>  _S14 = transpose_0(_S11.differential_0);
    *v_mean_0 = *v_mean_0 + v_viewdir_0;
    *v_R_0 = *v_R_0 + _S14;
    *v_t_0 = *v_t_0 + _S13.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_atomic(float3  mean_2, Matrix<float, 3, 3>  R_2, float3  t_2, float3  coeff_dc_2, float3  * coeffs_2, float3  v_colors_1, float3  * v_coeff_dc_1, float3  * v_coeffs_1, float3  * v_mean_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  colors_1 = make_float3 (0.282094806432724f) * coeff_dc_2;
    *v_coeff_dc_1 = *v_coeff_dc_1 + make_float3 (0.282094806432724f) * (v_colors_1 * make_float3 (float((colors_1.x) >= -0.5f), float((colors_1.y) >= -0.5f), float((colors_1.z) >= -0.5f)));
    float3  v_viewdir_1 = {};
    Matrix<float, 3, 3>  _S15 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S16;
    (&_S16)->primal_0 = transpose_0(R_2);
    (&_S16)->differential_0 = _S15;
    float3  _S17 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S18;
    (&_S18)->primal_0 = t_2;
    (&_S18)->differential_0 = _S17;
    s_bwd_prop_mul_0(&_S16, &_S18, v_viewdir_1);
    Matrix<float, 3, 3>  _S19 = transpose_0(_S16.differential_0);
    *v_mean_1 = *v_mean_1 + v_viewdir_1;
    *v_R_1 = *v_R_1 + _S19;
    *v_t_1 = *v_t_1 + _S18.differential_0;
    return;
}

inline __device__ float3  sh1_to_color(float3  mean_3, Matrix<float, 3, 3>  R_3, float3  t_3, float3  coeff_dc_3, float3  * coeffs_3)
{
    float3  _S20 = mean_3 + mul_0(transpose_0(R_3), t_3);
    float _S21 = _S20.x;
    float _S22 = _S20.y;
    float _S23 = _S20.z;
    float inv_norm_0 = (F32_rsqrt((_S21 * _S21 + _S22 * _S22 + _S23 * _S23)));
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_3 + make_float3 (0.48860251903533936f) * (make_float3 (- (_S22 * inv_norm_0)) * *(coeffs_3 + int(0)) + make_float3 (_S23 * inv_norm_0) * *(coeffs_3 + int(1)) - make_float3 (_S21 * inv_norm_0) * *(coeffs_3 + int(2))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh1_to_color_vjp_inplace(float3  mean_4, Matrix<float, 3, 3>  R_4, float3  t_4, float3  coeff_dc_4, float3  * coeffs_4, float3  v_colors_2, float3  * v_coeff_dc_2, float3  * v_coeffs_2, float3  * v_mean_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 3, 3>  _S24 = transpose_0(R_4);
    float3  _S25 = mean_4 + mul_0(_S24, t_4);
    float _S26 = _S25.x;
    float _S27 = _S25.y;
    float _S28 = _S25.z;
    float inv_norm_1 = (F32_rsqrt((_S26 * _S26 + _S27 * _S27 + _S28 * _S28)));
    float x_3 = _S26 * inv_norm_1;
    float y_2 = _S27 * inv_norm_1;
    float z_0 = _S28 * inv_norm_1;
    float3  * _S29 = coeffs_4 + int(0);
    float3  * _S30 = coeffs_4 + int(1);
    float3  * _S31 = coeffs_4 + int(2);
    float3  colors_2 = make_float3 (0.282094806432724f) * coeff_dc_4 + make_float3 (0.48860251903533936f) * (make_float3 (- y_2) * *_S29 + make_float3 (z_0) * *_S30 - make_float3 (x_3) * *_S31);
    float3  _S32 = v_colors_2 * make_float3 (float((colors_2.x) >= -0.5f), float((colors_2.y) >= -0.5f), float((colors_2.z) >= -0.5f));
    float3  v_viewdir_2 = {};
    *v_coeff_dc_2 = *v_coeff_dc_2 + make_float3 (0.282094806432724f) * _S32;
    float3  * _S33 = v_coeffs_2 + int(0);
    *_S33 = *_S33 + make_float3 (-0.48860251903533936f * y_2) * _S32;
    float3  * _S34 = v_coeffs_2 + int(1);
    *_S34 = *_S34 + make_float3 (0.48860251903533936f * z_0) * _S32;
    float3  * _S35 = v_coeffs_2 + int(2);
    *_S35 = *_S35 + make_float3 (-0.48860251903533936f * x_3) * _S32;
    float3  dir_n_0 = make_float3 (x_3, y_2, z_0);
    float3  v_dir_n_0 = make_float3 (-0.48860251903533936f * dot_0(*_S31, _S32), -0.48860251903533936f * dot_0(*_S29, _S32), 0.48860251903533936f * dot_0(*_S30, _S32));
    float3  v_viewdir_3 = v_viewdir_2 + (v_dir_n_0 - make_float3 (dot_0(v_dir_n_0, dir_n_0)) * dir_n_0) * make_float3 (inv_norm_1);
    Matrix<float, 3, 3>  _S36 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S37;
    (&_S37)->primal_0 = _S24;
    (&_S37)->differential_0 = _S36;
    float3  _S38 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S39;
    (&_S39)->primal_0 = t_4;
    (&_S39)->differential_0 = _S38;
    s_bwd_prop_mul_0(&_S37, &_S39, v_viewdir_3);
    Matrix<float, 3, 3>  _S40 = transpose_0(_S37.differential_0);
    *v_mean_2 = *v_mean_2 + v_viewdir_3;
    *v_R_2 = *v_R_2 + _S40;
    *v_t_2 = *v_t_2 + _S39.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_atomic(float3  mean_5, Matrix<float, 3, 3>  R_5, float3  t_5, float3  coeff_dc_5, float3  * coeffs_5, float3  v_colors_3, float3  * v_coeff_dc_3, float3  * v_coeffs_3, float3  * v_mean_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    Matrix<float, 3, 3>  _S41 = transpose_0(R_5);
    float3  _S42 = mean_5 + mul_0(_S41, t_5);
    float _S43 = _S42.x;
    float _S44 = _S42.y;
    float _S45 = _S42.z;
    float inv_norm_2 = (F32_rsqrt((_S43 * _S43 + _S44 * _S44 + _S45 * _S45)));
    float x_4 = _S43 * inv_norm_2;
    float y_3 = _S44 * inv_norm_2;
    float z_1 = _S45 * inv_norm_2;
    float3  * _S46 = coeffs_5 + int(0);
    float3  * _S47 = coeffs_5 + int(1);
    float3  * _S48 = coeffs_5 + int(2);
    float3  colors_3 = make_float3 (0.282094806432724f) * coeff_dc_5 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * *_S46 + make_float3 (z_1) * *_S47 - make_float3 (x_4) * *_S48);
    float3  _S49 = v_colors_3 * make_float3 (float((colors_3.x) >= -0.5f), float((colors_3.y) >= -0.5f), float((colors_3.z) >= -0.5f));
    float3  v_viewdir_4 = {};
    *v_coeff_dc_3 = *v_coeff_dc_3 + make_float3 (0.282094806432724f) * _S49;
    float3  temp_0 = make_float3 (-0.48860251903533936f * y_3) * _S49;
    float _S50 = dot_0(temp_0, temp_0);
    bool _S51;
    if((F32_isfinite((_S50))))
    {
        _S51 = _S50 != 0.0f;
    }
    else
    {
        _S51 = false;
    }
    if(_S51)
    {
        float3  * _S52 = v_coeffs_3 + int(0);
        float _S53 = atomicAdd(&(_S52->x), temp_0.x);
        float _S54 = atomicAdd(&(_S52->y), temp_0.y);
        float _S55 = atomicAdd(&(_S52->z), temp_0.z);
    }
    float3  temp_1 = make_float3 (0.48860251903533936f * z_1) * _S49;
    float _S56 = dot_0(temp_1, temp_1);
    if((F32_isfinite((_S56))))
    {
        _S51 = _S56 != 0.0f;
    }
    else
    {
        _S51 = false;
    }
    if(_S51)
    {
        float3  * _S57 = v_coeffs_3 + int(1);
        float _S58 = atomicAdd(&(_S57->x), temp_1.x);
        float _S59 = atomicAdd(&(_S57->y), temp_1.y);
        float _S60 = atomicAdd(&(_S57->z), temp_1.z);
    }
    float3  temp_2 = make_float3 (-0.48860251903533936f * x_4) * _S49;
    float _S61 = dot_0(temp_2, temp_2);
    if((F32_isfinite((_S61))))
    {
        _S51 = _S61 != 0.0f;
    }
    else
    {
        _S51 = false;
    }
    if(_S51)
    {
        float3  * _S62 = v_coeffs_3 + int(2);
        float _S63 = atomicAdd(&(_S62->x), temp_2.x);
        float _S64 = atomicAdd(&(_S62->y), temp_2.y);
        float _S65 = atomicAdd(&(_S62->z), temp_2.z);
    }
    float3  dir_n_1 = make_float3 (x_4, y_3, z_1);
    float3  v_dir_n_1 = make_float3 (-0.48860251903533936f * dot_0(*_S48, _S49), -0.48860251903533936f * dot_0(*_S46, _S49), 0.48860251903533936f * dot_0(*_S47, _S49));
    float3  v_viewdir_5 = v_viewdir_4 + (v_dir_n_1 - make_float3 (dot_0(v_dir_n_1, dir_n_1)) * dir_n_1) * make_float3 (inv_norm_2);
    Matrix<float, 3, 3>  _S66 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S67;
    (&_S67)->primal_0 = _S41;
    (&_S67)->differential_0 = _S66;
    float3  _S68 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S69;
    (&_S69)->primal_0 = t_5;
    (&_S69)->differential_0 = _S68;
    s_bwd_prop_mul_0(&_S67, &_S69, v_viewdir_5);
    Matrix<float, 3, 3>  _S70 = transpose_0(_S67.differential_0);
    *v_mean_3 = *v_mean_3 + v_viewdir_5;
    *v_R_3 = *v_R_3 + _S70;
    *v_t_3 = *v_t_3 + _S69.differential_0;
    return;
}

inline __device__ float3  sh2_to_color(float3  mean_6, Matrix<float, 3, 3>  R_6, float3  t_6, float3  coeff_dc_6, float3  * coeffs_6)
{
    float3  _S71 = mean_6 + mul_0(transpose_0(R_6), t_6);
    float _S72 = _S71.x;
    float _S73 = _S71.y;
    float _S74 = _S71.z;
    float inv_norm_3 = (F32_rsqrt((_S72 * _S72 + _S73 * _S73 + _S74 * _S74)));
    float x_5 = _S72 * inv_norm_3;
    float y_4 = _S73 * inv_norm_3;
    float z_2 = _S74 * inv_norm_3;
    float fTmp0B_0 = -1.09254848957061768f * z_2;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_6 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * *(coeffs_6 + int(0)) + make_float3 (z_2) * *(coeffs_6 + int(1)) - make_float3 (x_5) * *(coeffs_6 + int(2))) + (make_float3 (0.54627424478530884f * (2.0f * x_5 * y_4)) * *(coeffs_6 + int(3)) + make_float3 (fTmp0B_0 * y_4) * *(coeffs_6 + int(4)) + make_float3 (0.94617468118667603f * (z_2 * z_2) - 0.31539157032966614f) * *(coeffs_6 + int(5)) + make_float3 (fTmp0B_0 * x_5) * *(coeffs_6 + int(6)) + make_float3 (0.54627424478530884f * (x_5 * x_5 - y_4 * y_4)) * *(coeffs_6 + int(7))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh2_to_color_vjp_inplace(float3  mean_7, Matrix<float, 3, 3>  R_7, float3  t_7, float3  coeff_dc_7, float3  * coeffs_7, float3  v_colors_4, float3  * v_coeff_dc_4, float3  * v_coeffs_4, float3  * v_mean_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 3, 3>  _S75 = transpose_0(R_7);
    float3  _S76 = mean_7 + mul_0(_S75, t_7);
    float _S77 = _S76.x;
    float _S78 = _S76.y;
    float _S79 = _S76.z;
    float inv_norm_4 = (F32_rsqrt((_S77 * _S77 + _S78 * _S78 + _S79 * _S79)));
    float x_6 = _S77 * inv_norm_4;
    float y_5 = _S78 * inv_norm_4;
    float z_3 = _S79 * inv_norm_4;
    float3  * _S80 = coeffs_7 + int(0);
    float3  * _S81 = coeffs_7 + int(1);
    float3  * _S82 = coeffs_7 + int(2);
    float fTmp0B_1 = -1.09254848957061768f * z_3;
    float _S83 = 2.0f * x_6;
    float pSH6_0 = 0.94617468118667603f * (z_3 * z_3) - 0.31539157032966614f;
    float pSH7_0 = fTmp0B_1 * x_6;
    float pSH5_0 = fTmp0B_1 * y_5;
    float pSH8_0 = 0.54627424478530884f * (x_6 * x_6 - y_5 * y_5);
    float pSH4_0 = 0.54627424478530884f * (_S83 * y_5);
    float3  * _S84 = coeffs_7 + int(3);
    float3  * _S85 = coeffs_7 + int(4);
    float3  * _S86 = coeffs_7 + int(5);
    float3  * _S87 = coeffs_7 + int(6);
    float3  * _S88 = coeffs_7 + int(7);
    float3  colors_4 = make_float3 (0.282094806432724f) * coeff_dc_7 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * *_S80 + make_float3 (z_3) * *_S81 - make_float3 (x_6) * *_S82) + (make_float3 (pSH4_0) * *_S84 + make_float3 (pSH5_0) * *_S85 + make_float3 (pSH6_0) * *_S86 + make_float3 (pSH7_0) * *_S87 + make_float3 (pSH8_0) * *_S88);
    float3  _S89 = v_colors_4 * make_float3 (float((colors_4.x) >= -0.5f), float((colors_4.y) >= -0.5f), float((colors_4.z) >= -0.5f));
    *v_coeff_dc_4 = *v_coeff_dc_4 + make_float3 (0.282094806432724f) * _S89;
    float3  v_viewdir_6 = {};
    float3  * _S90 = v_coeffs_4 + int(0);
    *_S90 = *_S90 + make_float3 (-0.48860251903533936f * y_5) * _S89;
    float3  * _S91 = v_coeffs_4 + int(1);
    *_S91 = *_S91 + make_float3 (0.48860251903533936f * z_3) * _S89;
    float3  * _S92 = v_coeffs_4 + int(2);
    *_S92 = *_S92 + make_float3 (-0.48860251903533936f * x_6) * _S89;
    float _S93 = -0.48860251903533936f * dot_0(*_S82, _S89);
    float _S94 = -0.48860251903533936f * dot_0(*_S80, _S89);
    float _S95 = 0.48860251903533936f * dot_0(*_S81, _S89);
    float3  * _S96 = v_coeffs_4 + int(3);
    *_S96 = *_S96 + make_float3 (pSH4_0) * _S89;
    float3  * _S97 = v_coeffs_4 + int(4);
    *_S97 = *_S97 + make_float3 (pSH5_0) * _S89;
    float3  * _S98 = v_coeffs_4 + int(5);
    *_S98 = *_S98 + make_float3 (pSH6_0) * _S89;
    float3  * _S99 = v_coeffs_4 + int(6);
    *_S99 = *_S99 + make_float3 (pSH7_0) * _S89;
    float3  * _S100 = v_coeffs_4 + int(7);
    *_S100 = *_S100 + make_float3 (pSH8_0) * _S89;
    float pSH8_x_0 = 0.54627424478530884f * _S83;
    float3  dir_n_2 = make_float3 (x_6, y_5, z_3);
    float3  v_dir_n_2 = make_float3 (_S93 + dot_0(_S89, make_float3 (0.54627424478530884f * (2.0f * y_5)) * *_S84 + make_float3 (pSH8_x_0) * *_S88 + make_float3 (fTmp0B_1) * *_S87), _S94 + dot_0(_S89, make_float3 (pSH8_x_0) * *_S84 + make_float3 (0.54627424478530884f * (-2.0f * y_5)) * *_S88 + make_float3 (fTmp0B_1) * *_S85), _S95 + dot_0(_S89, make_float3 (1.89234936237335205f * z_3) * *_S86 + make_float3 (-1.09254848957061768f * x_6) * *_S87 + make_float3 (-1.09254848957061768f * y_5) * *_S85));
    float3  v_viewdir_7 = v_viewdir_6 + (v_dir_n_2 - make_float3 (dot_0(v_dir_n_2, dir_n_2)) * dir_n_2) * make_float3 (inv_norm_4);
    Matrix<float, 3, 3>  _S101 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S102;
    (&_S102)->primal_0 = _S75;
    (&_S102)->differential_0 = _S101;
    float3  _S103 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S104;
    (&_S104)->primal_0 = t_7;
    (&_S104)->differential_0 = _S103;
    s_bwd_prop_mul_0(&_S102, &_S104, v_viewdir_7);
    Matrix<float, 3, 3>  _S105 = transpose_0(_S102.differential_0);
    *v_mean_4 = *v_mean_4 + v_viewdir_7;
    *v_R_4 = *v_R_4 + _S105;
    *v_t_4 = *v_t_4 + _S104.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_atomic(float3  mean_8, Matrix<float, 3, 3>  R_8, float3  t_8, float3  coeff_dc_8, float3  * coeffs_8, float3  v_colors_5, float3  * v_coeff_dc_5, float3  * v_coeffs_5, float3  * v_mean_5, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 3, 3>  _S106 = transpose_0(R_8);
    float3  _S107 = mean_8 + mul_0(_S106, t_8);
    float _S108 = _S107.x;
    float _S109 = _S107.y;
    float _S110 = _S107.z;
    float inv_norm_5 = (F32_rsqrt((_S108 * _S108 + _S109 * _S109 + _S110 * _S110)));
    float x_7 = _S108 * inv_norm_5;
    float y_6 = _S109 * inv_norm_5;
    float z_4 = _S110 * inv_norm_5;
    float3  * _S111 = coeffs_8 + int(0);
    float3  * _S112 = coeffs_8 + int(1);
    float3  * _S113 = coeffs_8 + int(2);
    float fTmp0B_2 = -1.09254848957061768f * z_4;
    float _S114 = 2.0f * x_7;
    float pSH6_1 = 0.94617468118667603f * (z_4 * z_4) - 0.31539157032966614f;
    float pSH7_1 = fTmp0B_2 * x_7;
    float pSH5_1 = fTmp0B_2 * y_6;
    float pSH8_1 = 0.54627424478530884f * (x_7 * x_7 - y_6 * y_6);
    float pSH4_1 = 0.54627424478530884f * (_S114 * y_6);
    float3  * _S115 = coeffs_8 + int(3);
    float3  * _S116 = coeffs_8 + int(4);
    float3  * _S117 = coeffs_8 + int(5);
    float3  * _S118 = coeffs_8 + int(6);
    float3  * _S119 = coeffs_8 + int(7);
    float3  colors_5 = make_float3 (0.282094806432724f) * coeff_dc_8 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * *_S111 + make_float3 (z_4) * *_S112 - make_float3 (x_7) * *_S113) + (make_float3 (pSH4_1) * *_S115 + make_float3 (pSH5_1) * *_S116 + make_float3 (pSH6_1) * *_S117 + make_float3 (pSH7_1) * *_S118 + make_float3 (pSH8_1) * *_S119);
    float3  _S120 = v_colors_5 * make_float3 (float((colors_5.x) >= -0.5f), float((colors_5.y) >= -0.5f), float((colors_5.z) >= -0.5f));
    *v_coeff_dc_5 = *v_coeff_dc_5 + make_float3 (0.282094806432724f) * _S120;
    float3  v_viewdir_8 = {};
    float3  temp_3 = make_float3 (-0.48860251903533936f * y_6) * _S120;
    float _S121 = dot_0(temp_3, temp_3);
    bool _S122;
    if((F32_isfinite((_S121))))
    {
        _S122 = _S121 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S123 = v_coeffs_5 + int(0);
        float _S124 = atomicAdd(&(_S123->x), temp_3.x);
        float _S125 = atomicAdd(&(_S123->y), temp_3.y);
        float _S126 = atomicAdd(&(_S123->z), temp_3.z);
    }
    float3  temp_4 = make_float3 (0.48860251903533936f * z_4) * _S120;
    float _S127 = dot_0(temp_4, temp_4);
    if((F32_isfinite((_S127))))
    {
        _S122 = _S127 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S128 = v_coeffs_5 + int(1);
        float _S129 = atomicAdd(&(_S128->x), temp_4.x);
        float _S130 = atomicAdd(&(_S128->y), temp_4.y);
        float _S131 = atomicAdd(&(_S128->z), temp_4.z);
    }
    float3  temp_5 = make_float3 (-0.48860251903533936f * x_7) * _S120;
    float _S132 = dot_0(temp_5, temp_5);
    if((F32_isfinite((_S132))))
    {
        _S122 = _S132 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S133 = v_coeffs_5 + int(2);
        float _S134 = atomicAdd(&(_S133->x), temp_5.x);
        float _S135 = atomicAdd(&(_S133->y), temp_5.y);
        float _S136 = atomicAdd(&(_S133->z), temp_5.z);
    }
    float _S137 = -0.48860251903533936f * dot_0(*_S113, _S120);
    float _S138 = -0.48860251903533936f * dot_0(*_S111, _S120);
    float _S139 = 0.48860251903533936f * dot_0(*_S112, _S120);
    float3  temp_6 = make_float3 (pSH4_1) * _S120;
    float _S140 = dot_0(temp_6, temp_6);
    if((F32_isfinite((_S140))))
    {
        _S122 = _S140 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S141 = v_coeffs_5 + int(3);
        float _S142 = atomicAdd(&(_S141->x), temp_6.x);
        float _S143 = atomicAdd(&(_S141->y), temp_6.y);
        float _S144 = atomicAdd(&(_S141->z), temp_6.z);
    }
    float3  temp_7 = make_float3 (pSH5_1) * _S120;
    float _S145 = dot_0(temp_7, temp_7);
    if((F32_isfinite((_S145))))
    {
        _S122 = _S145 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S146 = v_coeffs_5 + int(4);
        float _S147 = atomicAdd(&(_S146->x), temp_7.x);
        float _S148 = atomicAdd(&(_S146->y), temp_7.y);
        float _S149 = atomicAdd(&(_S146->z), temp_7.z);
    }
    float3  temp_8 = make_float3 (pSH6_1) * _S120;
    float _S150 = dot_0(temp_8, temp_8);
    if((F32_isfinite((_S150))))
    {
        _S122 = _S150 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S151 = v_coeffs_5 + int(5);
        float _S152 = atomicAdd(&(_S151->x), temp_8.x);
        float _S153 = atomicAdd(&(_S151->y), temp_8.y);
        float _S154 = atomicAdd(&(_S151->z), temp_8.z);
    }
    float3  temp_9 = make_float3 (pSH7_1) * _S120;
    float _S155 = dot_0(temp_9, temp_9);
    if((F32_isfinite((_S155))))
    {
        _S122 = _S155 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S156 = v_coeffs_5 + int(6);
        float _S157 = atomicAdd(&(_S156->x), temp_9.x);
        float _S158 = atomicAdd(&(_S156->y), temp_9.y);
        float _S159 = atomicAdd(&(_S156->z), temp_9.z);
    }
    float3  temp_10 = make_float3 (pSH8_1) * _S120;
    float _S160 = dot_0(temp_10, temp_10);
    if((F32_isfinite((_S160))))
    {
        _S122 = _S160 != 0.0f;
    }
    else
    {
        _S122 = false;
    }
    if(_S122)
    {
        float3  * _S161 = v_coeffs_5 + int(7);
        float _S162 = atomicAdd(&(_S161->x), temp_10.x);
        float _S163 = atomicAdd(&(_S161->y), temp_10.y);
        float _S164 = atomicAdd(&(_S161->z), temp_10.z);
    }
    float pSH8_x_1 = 0.54627424478530884f * _S114;
    float3  dir_n_3 = make_float3 (x_7, y_6, z_4);
    float3  v_dir_n_3 = make_float3 (_S137 + dot_0(_S120, make_float3 (0.54627424478530884f * (2.0f * y_6)) * *_S115 + make_float3 (pSH8_x_1) * *_S119 + make_float3 (fTmp0B_2) * *_S118), _S138 + dot_0(_S120, make_float3 (pSH8_x_1) * *_S115 + make_float3 (0.54627424478530884f * (-2.0f * y_6)) * *_S119 + make_float3 (fTmp0B_2) * *_S116), _S139 + dot_0(_S120, make_float3 (1.89234936237335205f * z_4) * *_S117 + make_float3 (-1.09254848957061768f * x_7) * *_S118 + make_float3 (-1.09254848957061768f * y_6) * *_S116));
    float3  v_viewdir_9 = v_viewdir_8 + (v_dir_n_3 - make_float3 (dot_0(v_dir_n_3, dir_n_3)) * dir_n_3) * make_float3 (inv_norm_5);
    Matrix<float, 3, 3>  _S165 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S166;
    (&_S166)->primal_0 = _S106;
    (&_S166)->differential_0 = _S165;
    float3  _S167 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S168;
    (&_S168)->primal_0 = t_8;
    (&_S168)->differential_0 = _S167;
    s_bwd_prop_mul_0(&_S166, &_S168, v_viewdir_9);
    Matrix<float, 3, 3>  _S169 = transpose_0(_S166.differential_0);
    *v_mean_5 = *v_mean_5 + v_viewdir_9;
    *v_R_5 = *v_R_5 + _S169;
    *v_t_5 = *v_t_5 + _S168.differential_0;
    return;
}

inline __device__ float3  sh3_to_color(float3  mean_9, Matrix<float, 3, 3>  R_9, float3  t_9, float3  coeff_dc_9, float3  * coeffs_9)
{
    float3  _S170 = mean_9 + mul_0(transpose_0(R_9), t_9);
    float _S171 = _S170.x;
    float _S172 = _S170.y;
    float _S173 = _S170.z;
    float inv_norm_6 = (F32_rsqrt((_S171 * _S171 + _S172 * _S172 + _S173 * _S173)));
    float x_8 = _S171 * inv_norm_6;
    float y_7 = _S172 * inv_norm_6;
    float z_5 = _S173 * inv_norm_6;
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
    Matrix<float, 3, 3>  _S174 = transpose_0(R_10);
    float3  _S175 = mean_10 + mul_0(_S174, t_10);
    float _S176 = _S175.x;
    float _S177 = _S175.y;
    float _S178 = _S175.z;
    float inv_norm_7 = (F32_rsqrt((_S176 * _S176 + _S177 * _S177 + _S178 * _S178)));
    float x_9 = _S176 * inv_norm_7;
    float y_8 = _S177 * inv_norm_7;
    float z_6 = _S178 * inv_norm_7;
    float3  * _S179 = coeffs_10 + int(0);
    float3  * _S180 = coeffs_10 + int(1);
    float3  * _S181 = coeffs_10 + int(2);
    float z2_1 = z_6 * z_6;
    float fTmp0B_4 = -1.09254848957061768f * z_6;
    float fC1_1 = x_9 * x_9 - y_8 * y_8;
    float _S182 = 2.0f * x_9;
    float fS1_1 = _S182 * y_8;
    float pSH6_2 = 0.94617468118667603f * z2_1 - 0.31539157032966614f;
    float pSH7_2 = fTmp0B_4 * x_9;
    float pSH5_2 = fTmp0B_4 * y_8;
    float pSH8_2 = 0.54627424478530884f * fC1_1;
    float pSH4_2 = 0.54627424478530884f * fS1_1;
    float3  * _S183 = coeffs_10 + int(3);
    float3  * _S184 = coeffs_10 + int(4);
    float3  * _S185 = coeffs_10 + int(5);
    float3  * _S186 = coeffs_10 + int(6);
    float3  * _S187 = coeffs_10 + int(7);
    float fTmp0C_1 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
    float fTmp1B_1 = 1.44530570507049561f * z_6;
    float pSH12_0 = z_6 * (1.86588168144226074f * z2_1 - 1.11952900886535645f);
    float pSH13_0 = fTmp0C_1 * x_9;
    float pSH11_0 = fTmp0C_1 * y_8;
    float pSH14_0 = fTmp1B_1 * fC1_1;
    float pSH10_0 = fTmp1B_1 * fS1_1;
    float pSH15_0 = -0.59004360437393188f * (x_9 * fC1_1 - y_8 * fS1_1);
    float pSH9_0 = -0.59004360437393188f * (x_9 * fS1_1 + y_8 * fC1_1);
    float3  * _S188 = coeffs_10 + int(8);
    float3  * _S189 = coeffs_10 + int(9);
    float3  * _S190 = coeffs_10 + int(10);
    float3  * _S191 = coeffs_10 + int(11);
    float3  * _S192 = coeffs_10 + int(12);
    float3  * _S193 = coeffs_10 + int(13);
    float3  * _S194 = coeffs_10 + int(14);
    float3  colors_6 = make_float3 (0.282094806432724f) * coeff_dc_10 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * *_S179 + make_float3 (z_6) * *_S180 - make_float3 (x_9) * *_S181) + (make_float3 (pSH4_2) * *_S183 + make_float3 (pSH5_2) * *_S184 + make_float3 (pSH6_2) * *_S185 + make_float3 (pSH7_2) * *_S186 + make_float3 (pSH8_2) * *_S187) + (make_float3 (pSH9_0) * *_S188 + make_float3 (pSH10_0) * *_S189 + make_float3 (pSH11_0) * *_S190 + make_float3 (pSH12_0) * *_S191 + make_float3 (pSH13_0) * *_S192 + make_float3 (pSH14_0) * *_S193 + make_float3 (pSH15_0) * *_S194);
    float3  _S195 = v_colors_6 * make_float3 (float((colors_6.x) >= -0.5f), float((colors_6.y) >= -0.5f), float((colors_6.z) >= -0.5f));
    float3  v_viewdir_10 = {};
    *v_coeff_dc_6 = *v_coeff_dc_6 + make_float3 (0.282094806432724f) * _S195;
    float3  * _S196 = v_coeffs_6 + int(0);
    *_S196 = *_S196 + make_float3 (-0.48860251903533936f * y_8) * _S195;
    float3  * _S197 = v_coeffs_6 + int(1);
    *_S197 = *_S197 + make_float3 (0.48860251903533936f * z_6) * _S195;
    float3  * _S198 = v_coeffs_6 + int(2);
    *_S198 = *_S198 + make_float3 (-0.48860251903533936f * x_9) * _S195;
    float _S199 = -0.48860251903533936f * dot_0(*_S181, _S195);
    float _S200 = -0.48860251903533936f * dot_0(*_S179, _S195);
    float _S201 = 0.48860251903533936f * dot_0(*_S180, _S195);
    float3  * _S202 = v_coeffs_6 + int(3);
    *_S202 = *_S202 + make_float3 (pSH4_2) * _S195;
    float3  * _S203 = v_coeffs_6 + int(4);
    *_S203 = *_S203 + make_float3 (pSH5_2) * _S195;
    float3  * _S204 = v_coeffs_6 + int(5);
    *_S204 = *_S204 + make_float3 (pSH6_2) * _S195;
    float3  * _S205 = v_coeffs_6 + int(6);
    *_S205 = *_S205 + make_float3 (pSH7_2) * _S195;
    float3  * _S206 = v_coeffs_6 + int(7);
    *_S206 = *_S206 + make_float3 (pSH8_2) * _S195;
    float fC1_y_0 = -2.0f * y_8;
    float fS1_x_0 = 2.0f * y_8;
    float pSH8_x_2 = 0.54627424478530884f * _S182;
    float v_x_0 = _S199 + dot_0(_S195, make_float3 (0.54627424478530884f * fS1_x_0) * *_S183 + make_float3 (pSH8_x_2) * *_S187 + make_float3 (fTmp0B_4) * *_S186);
    float v_y_0 = _S200 + dot_0(_S195, make_float3 (pSH8_x_2) * *_S183 + make_float3 (0.54627424478530884f * fC1_y_0) * *_S187 + make_float3 (fTmp0B_4) * *_S184);
    float v_z_0 = _S201 + dot_0(_S195, make_float3 (1.89234936237335205f * z_6) * *_S185 + make_float3 (-1.09254848957061768f * x_9) * *_S186 + make_float3 (-1.09254848957061768f * y_8) * *_S184);
    float3  * _S207 = v_coeffs_6 + int(8);
    *_S207 = *_S207 + make_float3 (pSH9_0) * _S195;
    float3  * _S208 = v_coeffs_6 + int(9);
    *_S208 = *_S208 + make_float3 (pSH10_0) * _S195;
    float3  * _S209 = v_coeffs_6 + int(10);
    *_S209 = *_S209 + make_float3 (pSH11_0) * _S195;
    float3  * _S210 = v_coeffs_6 + int(11);
    *_S210 = *_S210 + make_float3 (pSH12_0) * _S195;
    float3  * _S211 = v_coeffs_6 + int(12);
    *_S211 = *_S211 + make_float3 (pSH13_0) * _S195;
    float3  * _S212 = v_coeffs_6 + int(13);
    *_S212 = *_S212 + make_float3 (pSH14_0) * _S195;
    float3  * _S213 = v_coeffs_6 + int(14);
    *_S213 = *_S213 + make_float3 (pSH15_0) * _S195;
    float fTmp0C_z_0 = -4.57045793533325195f * z_6;
    float _S214 = x_9 * _S182;
    float _S215 = y_8 * _S182;
    float pSH14_x_0 = fTmp1B_1 * _S182;
    float3  dir_n_4 = make_float3 (x_9, y_8, z_6);
    float3  v_dir_n_4 = make_float3 (v_x_0 + dot_0(_S195, make_float3 (-0.59004360437393188f * (fS1_1 + x_9 * fS1_x_0 + _S215)) * *_S188 + make_float3 (-0.59004360437393188f * (fC1_1 + _S214 - y_8 * fS1_x_0)) * *_S194 + make_float3 (fTmp1B_1 * fS1_x_0) * *_S189 + make_float3 (pSH14_x_0) * *_S193 + make_float3 (fTmp0C_1) * *_S192), v_y_0 + dot_0(_S195, make_float3 (-0.59004360437393188f * (_S214 + fC1_1 + y_8 * fC1_y_0)) * *_S188 + make_float3 (-0.59004360437393188f * (x_9 * fC1_y_0 - fS1_1 - _S215)) * *_S194 + make_float3 (pSH14_x_0) * *_S189 + make_float3 (fTmp1B_1 * fC1_y_0) * *_S193 + make_float3 (fTmp0C_1) * *_S190), v_z_0 + dot_0(_S195, make_float3 (5.59764480590820312f * z2_1 - 1.11952900886535645f) * *_S191 + make_float3 (fTmp0C_z_0 * x_9) * *_S192 + make_float3 (fTmp0C_z_0 * y_8) * *_S190 + make_float3 (1.44530570507049561f * fC1_1) * *_S193 + make_float3 (1.44530570507049561f * fS1_1) * *_S189));
    float3  v_viewdir_11 = v_viewdir_10 + (v_dir_n_4 - make_float3 (dot_0(v_dir_n_4, dir_n_4)) * dir_n_4) * make_float3 (inv_norm_7);
    Matrix<float, 3, 3>  _S216 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S217;
    (&_S217)->primal_0 = _S174;
    (&_S217)->differential_0 = _S216;
    float3  _S218 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S219;
    (&_S219)->primal_0 = t_10;
    (&_S219)->differential_0 = _S218;
    s_bwd_prop_mul_0(&_S217, &_S219, v_viewdir_11);
    Matrix<float, 3, 3>  _S220 = transpose_0(_S217.differential_0);
    *v_mean_6 = *v_mean_6 + v_viewdir_11;
    *v_R_6 = *v_R_6 + _S220;
    *v_t_6 = *v_t_6 + _S219.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_atomic(float3  mean_11, Matrix<float, 3, 3>  R_11, float3  t_11, float3  coeff_dc_11, float3  * coeffs_11, float3  v_colors_7, float3  * v_coeff_dc_7, float3  * v_coeffs_7, float3  * v_mean_7, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    Matrix<float, 3, 3>  _S221 = transpose_0(R_11);
    float3  _S222 = mean_11 + mul_0(_S221, t_11);
    float _S223 = _S222.x;
    float _S224 = _S222.y;
    float _S225 = _S222.z;
    float inv_norm_8 = (F32_rsqrt((_S223 * _S223 + _S224 * _S224 + _S225 * _S225)));
    float x_10 = _S223 * inv_norm_8;
    float y_9 = _S224 * inv_norm_8;
    float z_7 = _S225 * inv_norm_8;
    float3  * _S226 = coeffs_11 + int(0);
    float3  * _S227 = coeffs_11 + int(1);
    float3  * _S228 = coeffs_11 + int(2);
    float z2_2 = z_7 * z_7;
    float fTmp0B_5 = -1.09254848957061768f * z_7;
    float fC1_2 = x_10 * x_10 - y_9 * y_9;
    float _S229 = 2.0f * x_10;
    float fS1_2 = _S229 * y_9;
    float pSH6_3 = 0.94617468118667603f * z2_2 - 0.31539157032966614f;
    float pSH7_3 = fTmp0B_5 * x_10;
    float pSH5_3 = fTmp0B_5 * y_9;
    float pSH8_3 = 0.54627424478530884f * fC1_2;
    float pSH4_3 = 0.54627424478530884f * fS1_2;
    float3  * _S230 = coeffs_11 + int(3);
    float3  * _S231 = coeffs_11 + int(4);
    float3  * _S232 = coeffs_11 + int(5);
    float3  * _S233 = coeffs_11 + int(6);
    float3  * _S234 = coeffs_11 + int(7);
    float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
    float fTmp1B_2 = 1.44530570507049561f * z_7;
    float pSH12_1 = z_7 * (1.86588168144226074f * z2_2 - 1.11952900886535645f);
    float pSH13_1 = fTmp0C_2 * x_10;
    float pSH11_1 = fTmp0C_2 * y_9;
    float pSH14_1 = fTmp1B_2 * fC1_2;
    float pSH10_1 = fTmp1B_2 * fS1_2;
    float pSH15_1 = -0.59004360437393188f * (x_10 * fC1_2 - y_9 * fS1_2);
    float pSH9_1 = -0.59004360437393188f * (x_10 * fS1_2 + y_9 * fC1_2);
    float3  * _S235 = coeffs_11 + int(8);
    float3  * _S236 = coeffs_11 + int(9);
    float3  * _S237 = coeffs_11 + int(10);
    float3  * _S238 = coeffs_11 + int(11);
    float3  * _S239 = coeffs_11 + int(12);
    float3  * _S240 = coeffs_11 + int(13);
    float3  * _S241 = coeffs_11 + int(14);
    float3  colors_7 = make_float3 (0.282094806432724f) * coeff_dc_11 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * *_S226 + make_float3 (z_7) * *_S227 - make_float3 (x_10) * *_S228) + (make_float3 (pSH4_3) * *_S230 + make_float3 (pSH5_3) * *_S231 + make_float3 (pSH6_3) * *_S232 + make_float3 (pSH7_3) * *_S233 + make_float3 (pSH8_3) * *_S234) + (make_float3 (pSH9_1) * *_S235 + make_float3 (pSH10_1) * *_S236 + make_float3 (pSH11_1) * *_S237 + make_float3 (pSH12_1) * *_S238 + make_float3 (pSH13_1) * *_S239 + make_float3 (pSH14_1) * *_S240 + make_float3 (pSH15_1) * *_S241);
    float3  _S242 = v_colors_7 * make_float3 (float((colors_7.x) >= -0.5f), float((colors_7.y) >= -0.5f), float((colors_7.z) >= -0.5f));
    float3  v_viewdir_12 = {};
    *v_coeff_dc_7 = *v_coeff_dc_7 + make_float3 (0.282094806432724f) * _S242;
    float3  temp_11 = make_float3 (-0.48860251903533936f * y_9) * _S242;
    float _S243 = dot_0(temp_11, temp_11);
    bool _S244;
    if((F32_isfinite((_S243))))
    {
        _S244 = _S243 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S245 = v_coeffs_7 + int(0);
        float _S246 = atomicAdd(&(_S245->x), temp_11.x);
        float _S247 = atomicAdd(&(_S245->y), temp_11.y);
        float _S248 = atomicAdd(&(_S245->z), temp_11.z);
    }
    float3  temp_12 = make_float3 (0.48860251903533936f * z_7) * _S242;
    float _S249 = dot_0(temp_12, temp_12);
    if((F32_isfinite((_S249))))
    {
        _S244 = _S249 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S250 = v_coeffs_7 + int(1);
        float _S251 = atomicAdd(&(_S250->x), temp_12.x);
        float _S252 = atomicAdd(&(_S250->y), temp_12.y);
        float _S253 = atomicAdd(&(_S250->z), temp_12.z);
    }
    float3  temp_13 = make_float3 (-0.48860251903533936f * x_10) * _S242;
    float _S254 = dot_0(temp_13, temp_13);
    if((F32_isfinite((_S254))))
    {
        _S244 = _S254 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S255 = v_coeffs_7 + int(2);
        float _S256 = atomicAdd(&(_S255->x), temp_13.x);
        float _S257 = atomicAdd(&(_S255->y), temp_13.y);
        float _S258 = atomicAdd(&(_S255->z), temp_13.z);
    }
    float _S259 = -0.48860251903533936f * dot_0(*_S228, _S242);
    float _S260 = -0.48860251903533936f * dot_0(*_S226, _S242);
    float _S261 = 0.48860251903533936f * dot_0(*_S227, _S242);
    float3  temp_14 = make_float3 (pSH4_3) * _S242;
    float _S262 = dot_0(temp_14, temp_14);
    if((F32_isfinite((_S262))))
    {
        _S244 = _S262 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S263 = v_coeffs_7 + int(3);
        float _S264 = atomicAdd(&(_S263->x), temp_14.x);
        float _S265 = atomicAdd(&(_S263->y), temp_14.y);
        float _S266 = atomicAdd(&(_S263->z), temp_14.z);
    }
    float3  temp_15 = make_float3 (pSH5_3) * _S242;
    float _S267 = dot_0(temp_15, temp_15);
    if((F32_isfinite((_S267))))
    {
        _S244 = _S267 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S268 = v_coeffs_7 + int(4);
        float _S269 = atomicAdd(&(_S268->x), temp_15.x);
        float _S270 = atomicAdd(&(_S268->y), temp_15.y);
        float _S271 = atomicAdd(&(_S268->z), temp_15.z);
    }
    float3  temp_16 = make_float3 (pSH6_3) * _S242;
    float _S272 = dot_0(temp_16, temp_16);
    if((F32_isfinite((_S272))))
    {
        _S244 = _S272 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S273 = v_coeffs_7 + int(5);
        float _S274 = atomicAdd(&(_S273->x), temp_16.x);
        float _S275 = atomicAdd(&(_S273->y), temp_16.y);
        float _S276 = atomicAdd(&(_S273->z), temp_16.z);
    }
    float3  temp_17 = make_float3 (pSH7_3) * _S242;
    float _S277 = dot_0(temp_17, temp_17);
    if((F32_isfinite((_S277))))
    {
        _S244 = _S277 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S278 = v_coeffs_7 + int(6);
        float _S279 = atomicAdd(&(_S278->x), temp_17.x);
        float _S280 = atomicAdd(&(_S278->y), temp_17.y);
        float _S281 = atomicAdd(&(_S278->z), temp_17.z);
    }
    float3  temp_18 = make_float3 (pSH8_3) * _S242;
    float _S282 = dot_0(temp_18, temp_18);
    if((F32_isfinite((_S282))))
    {
        _S244 = _S282 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S283 = v_coeffs_7 + int(7);
        float _S284 = atomicAdd(&(_S283->x), temp_18.x);
        float _S285 = atomicAdd(&(_S283->y), temp_18.y);
        float _S286 = atomicAdd(&(_S283->z), temp_18.z);
    }
    float fC1_y_1 = -2.0f * y_9;
    float fS1_x_1 = 2.0f * y_9;
    float pSH8_x_3 = 0.54627424478530884f * _S229;
    float v_x_1 = _S259 + dot_0(_S242, make_float3 (0.54627424478530884f * fS1_x_1) * *_S230 + make_float3 (pSH8_x_3) * *_S234 + make_float3 (fTmp0B_5) * *_S233);
    float v_y_1 = _S260 + dot_0(_S242, make_float3 (pSH8_x_3) * *_S230 + make_float3 (0.54627424478530884f * fC1_y_1) * *_S234 + make_float3 (fTmp0B_5) * *_S231);
    float v_z_1 = _S261 + dot_0(_S242, make_float3 (1.89234936237335205f * z_7) * *_S232 + make_float3 (-1.09254848957061768f * x_10) * *_S233 + make_float3 (-1.09254848957061768f * y_9) * *_S231);
    float3  temp_19 = make_float3 (pSH9_1) * _S242;
    float _S287 = dot_0(temp_19, temp_19);
    if((F32_isfinite((_S287))))
    {
        _S244 = _S287 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S288 = v_coeffs_7 + int(8);
        float _S289 = atomicAdd(&(_S288->x), temp_19.x);
        float _S290 = atomicAdd(&(_S288->y), temp_19.y);
        float _S291 = atomicAdd(&(_S288->z), temp_19.z);
    }
    float3  temp_20 = make_float3 (pSH10_1) * _S242;
    float _S292 = dot_0(temp_20, temp_20);
    if((F32_isfinite((_S292))))
    {
        _S244 = _S292 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S293 = v_coeffs_7 + int(9);
        float _S294 = atomicAdd(&(_S293->x), temp_20.x);
        float _S295 = atomicAdd(&(_S293->y), temp_20.y);
        float _S296 = atomicAdd(&(_S293->z), temp_20.z);
    }
    float3  temp_21 = make_float3 (pSH11_1) * _S242;
    float _S297 = dot_0(temp_21, temp_21);
    if((F32_isfinite((_S297))))
    {
        _S244 = _S297 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S298 = v_coeffs_7 + int(10);
        float _S299 = atomicAdd(&(_S298->x), temp_21.x);
        float _S300 = atomicAdd(&(_S298->y), temp_21.y);
        float _S301 = atomicAdd(&(_S298->z), temp_21.z);
    }
    float3  temp_22 = make_float3 (pSH12_1) * _S242;
    float _S302 = dot_0(temp_22, temp_22);
    if((F32_isfinite((_S302))))
    {
        _S244 = _S302 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S303 = v_coeffs_7 + int(11);
        float _S304 = atomicAdd(&(_S303->x), temp_22.x);
        float _S305 = atomicAdd(&(_S303->y), temp_22.y);
        float _S306 = atomicAdd(&(_S303->z), temp_22.z);
    }
    float3  temp_23 = make_float3 (pSH13_1) * _S242;
    float _S307 = dot_0(temp_23, temp_23);
    if((F32_isfinite((_S307))))
    {
        _S244 = _S307 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S308 = v_coeffs_7 + int(12);
        float _S309 = atomicAdd(&(_S308->x), temp_23.x);
        float _S310 = atomicAdd(&(_S308->y), temp_23.y);
        float _S311 = atomicAdd(&(_S308->z), temp_23.z);
    }
    float3  temp_24 = make_float3 (pSH14_1) * _S242;
    float _S312 = dot_0(temp_24, temp_24);
    if((F32_isfinite((_S312))))
    {
        _S244 = _S312 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S313 = v_coeffs_7 + int(13);
        float _S314 = atomicAdd(&(_S313->x), temp_24.x);
        float _S315 = atomicAdd(&(_S313->y), temp_24.y);
        float _S316 = atomicAdd(&(_S313->z), temp_24.z);
    }
    float3  temp_25 = make_float3 (pSH15_1) * _S242;
    float _S317 = dot_0(temp_25, temp_25);
    if((F32_isfinite((_S317))))
    {
        _S244 = _S317 != 0.0f;
    }
    else
    {
        _S244 = false;
    }
    if(_S244)
    {
        float3  * _S318 = v_coeffs_7 + int(14);
        float _S319 = atomicAdd(&(_S318->x), temp_25.x);
        float _S320 = atomicAdd(&(_S318->y), temp_25.y);
        float _S321 = atomicAdd(&(_S318->z), temp_25.z);
    }
    float fTmp0C_z_1 = -4.57045793533325195f * z_7;
    float _S322 = x_10 * _S229;
    float _S323 = y_9 * _S229;
    float pSH14_x_1 = fTmp1B_2 * _S229;
    float3  dir_n_5 = make_float3 (x_10, y_9, z_7);
    float3  v_dir_n_5 = make_float3 (v_x_1 + dot_0(_S242, make_float3 (-0.59004360437393188f * (fS1_2 + x_10 * fS1_x_1 + _S323)) * *_S235 + make_float3 (-0.59004360437393188f * (fC1_2 + _S322 - y_9 * fS1_x_1)) * *_S241 + make_float3 (fTmp1B_2 * fS1_x_1) * *_S236 + make_float3 (pSH14_x_1) * *_S240 + make_float3 (fTmp0C_2) * *_S239), v_y_1 + dot_0(_S242, make_float3 (-0.59004360437393188f * (_S322 + fC1_2 + y_9 * fC1_y_1)) * *_S235 + make_float3 (-0.59004360437393188f * (x_10 * fC1_y_1 - fS1_2 - _S323)) * *_S241 + make_float3 (pSH14_x_1) * *_S236 + make_float3 (fTmp1B_2 * fC1_y_1) * *_S240 + make_float3 (fTmp0C_2) * *_S237), v_z_1 + dot_0(_S242, make_float3 (5.59764480590820312f * z2_2 - 1.11952900886535645f) * *_S238 + make_float3 (fTmp0C_z_1 * x_10) * *_S239 + make_float3 (fTmp0C_z_1 * y_9) * *_S237 + make_float3 (1.44530570507049561f * fC1_2) * *_S240 + make_float3 (1.44530570507049561f * fS1_2) * *_S236));
    float3  v_viewdir_13 = v_viewdir_12 + (v_dir_n_5 - make_float3 (dot_0(v_dir_n_5, dir_n_5)) * dir_n_5) * make_float3 (inv_norm_8);
    Matrix<float, 3, 3>  _S324 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S325;
    (&_S325)->primal_0 = _S221;
    (&_S325)->differential_0 = _S324;
    float3  _S326 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S327;
    (&_S327)->primal_0 = t_11;
    (&_S327)->differential_0 = _S326;
    s_bwd_prop_mul_0(&_S325, &_S327, v_viewdir_13);
    Matrix<float, 3, 3>  _S328 = transpose_0(_S325.differential_0);
    *v_mean_7 = *v_mean_7 + v_viewdir_13;
    *v_R_7 = *v_R_7 + _S328;
    *v_t_7 = *v_t_7 + _S327.differential_0;
    return;
}

inline __device__ float3  sh4_to_color(float3  mean_12, Matrix<float, 3, 3>  R_12, float3  t_12, float3  coeff_dc_12, float3  * coeffs_12)
{
    float3  _S329 = mean_12 + mul_0(transpose_0(R_12), t_12);
    float _S330 = _S329.x;
    float _S331 = _S329.y;
    float _S332 = _S329.z;
    float inv_norm_9 = (F32_rsqrt((_S330 * _S330 + _S331 * _S331 + _S332 * _S332)));
    float x_11 = _S330 * inv_norm_9;
    float y_10 = _S331 * inv_norm_9;
    float z_8 = _S332 * inv_norm_9;
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
    Matrix<float, 3, 3>  _S333 = transpose_0(R_13);
    float3  _S334 = mean_13 + mul_0(_S333, t_13);
    float _S335 = _S334.x;
    float _S336 = _S334.y;
    float _S337 = _S334.z;
    float inv_norm_10 = (F32_rsqrt((_S335 * _S335 + _S336 * _S336 + _S337 * _S337)));
    float x_12 = _S335 * inv_norm_10;
    float y_11 = _S336 * inv_norm_10;
    float z_9 = _S337 * inv_norm_10;
    float3  * _S338 = coeffs_13 + int(0);
    float3  * _S339 = coeffs_13 + int(1);
    float3  * _S340 = coeffs_13 + int(2);
    float z2_4 = z_9 * z_9;
    float fTmp0B_7 = -1.09254848957061768f * z_9;
    float fC1_4 = x_12 * x_12 - y_11 * y_11;
    float _S341 = 2.0f * x_12;
    float fS1_4 = _S341 * y_11;
    float pSH6_5 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
    float pSH7_4 = fTmp0B_7 * x_12;
    float pSH5_4 = fTmp0B_7 * y_11;
    float pSH8_4 = 0.54627424478530884f * fC1_4;
    float pSH4_4 = 0.54627424478530884f * fS1_4;
    float3  * _S342 = coeffs_13 + int(3);
    float3  * _S343 = coeffs_13 + int(4);
    float3  * _S344 = coeffs_13 + int(5);
    float3  * _S345 = coeffs_13 + int(6);
    float3  * _S346 = coeffs_13 + int(7);
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
    float3  * _S347 = coeffs_13 + int(8);
    float3  * _S348 = coeffs_13 + int(9);
    float3  * _S349 = coeffs_13 + int(10);
    float3  * _S350 = coeffs_13 + int(11);
    float3  * _S351 = coeffs_13 + int(12);
    float3  * _S352 = coeffs_13 + int(13);
    float3  * _S353 = coeffs_13 + int(14);
    float fTmp0D_1 = z_9 * (-4.68332576751708984f * z2_4 + 2.00713968276977539f);
    float fTmp1C_1 = 3.31161141395568848f * z2_4 - 0.47308734059333801f;
    float fTmp2B_1 = -1.77013075351715088f * z_9;
    float _S354 = 1.9843134880065918f * z_9 * pSH12_3;
    float pSH21_0 = fTmp0D_1 * x_12;
    float pSH19_0 = fTmp0D_1 * y_11;
    float pSH22_0 = fTmp1C_1 * fC1_4;
    float pSH18_0 = fTmp1C_1 * fS1_4;
    float pSH23_0 = fTmp2B_1 * fC2_1;
    float pSH17_0 = fTmp2B_1 * fS2_1;
    float pSH24_0 = 0.62583571672439575f * (x_12 * fC2_1 - y_11 * fS2_1);
    float pSH16_0 = 0.62583571672439575f * (x_12 * fS2_1 + y_11 * fC2_1);
    float3  * _S355 = coeffs_13 + int(15);
    float3  * _S356 = coeffs_13 + int(16);
    float3  * _S357 = coeffs_13 + int(17);
    float3  * _S358 = coeffs_13 + int(18);
    float3  * _S359 = coeffs_13 + int(19);
    float3  * _S360 = coeffs_13 + int(20);
    float3  * _S361 = coeffs_13 + int(21);
    float3  * _S362 = coeffs_13 + int(22);
    float3  * _S363 = coeffs_13 + int(23);
    float3  colors_8 = make_float3 (0.282094806432724f) * coeff_dc_13 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * *_S338 + make_float3 (z_9) * *_S339 - make_float3 (x_12) * *_S340) + (make_float3 (pSH4_4) * *_S342 + make_float3 (pSH5_4) * *_S343 + make_float3 (pSH6_5) * *_S344 + make_float3 (pSH7_4) * *_S345 + make_float3 (pSH8_4) * *_S346) + (make_float3 (pSH9_2) * *_S347 + make_float3 (pSH10_2) * *_S348 + make_float3 (pSH11_2) * *_S349 + make_float3 (pSH12_3) * *_S350 + make_float3 (pSH13_2) * *_S351 + make_float3 (pSH14_2) * *_S352 + make_float3 (pSH15_2) * *_S353) + (make_float3 (pSH16_0) * *_S355 + make_float3 (pSH17_0) * *_S356 + make_float3 (pSH18_0) * *_S357 + make_float3 (pSH19_0) * *_S358 + make_float3 (_S354 - 1.00623059272766113f * pSH6_5) * *_S359 + make_float3 (pSH21_0) * *_S360 + make_float3 (pSH22_0) * *_S361 + make_float3 (pSH23_0) * *_S362 + make_float3 (pSH24_0) * *_S363);
    float3  _S364 = v_colors_8 * make_float3 (float((colors_8.x) >= -0.5f), float((colors_8.y) >= -0.5f), float((colors_8.z) >= -0.5f));
    float3  v_viewdir_14 = {};
    *v_coeff_dc_8 = *v_coeff_dc_8 + make_float3 (0.282094806432724f) * _S364;
    float3  * _S365 = v_coeffs_8 + int(0);
    *_S365 = *_S365 + make_float3 (-0.48860251903533936f * y_11) * _S364;
    float3  * _S366 = v_coeffs_8 + int(1);
    *_S366 = *_S366 + make_float3 (0.48860251903533936f * z_9) * _S364;
    float3  * _S367 = v_coeffs_8 + int(2);
    *_S367 = *_S367 + make_float3 (-0.48860251903533936f * x_12) * _S364;
    float _S368 = -0.48860251903533936f * dot_0(*_S340, _S364);
    float _S369 = -0.48860251903533936f * dot_0(*_S338, _S364);
    float _S370 = 0.48860251903533936f * dot_0(*_S339, _S364);
    float3  * _S371 = v_coeffs_8 + int(3);
    *_S371 = *_S371 + make_float3 (pSH4_4) * _S364;
    float3  * _S372 = v_coeffs_8 + int(4);
    *_S372 = *_S372 + make_float3 (pSH5_4) * _S364;
    float3  * _S373 = v_coeffs_8 + int(5);
    *_S373 = *_S373 + make_float3 (pSH6_5) * _S364;
    float3  * _S374 = v_coeffs_8 + int(6);
    *_S374 = *_S374 + make_float3 (pSH7_4) * _S364;
    float3  * _S375 = v_coeffs_8 + int(7);
    *_S375 = *_S375 + make_float3 (pSH8_4) * _S364;
    float fC1_y_2 = -2.0f * y_11;
    float fS1_x_2 = 2.0f * y_11;
    float pSH6_z_0 = 1.89234936237335205f * z_9;
    float pSH8_x_4 = 0.54627424478530884f * _S341;
    float v_x_2 = _S368 + dot_0(_S364, make_float3 (0.54627424478530884f * fS1_x_2) * *_S342 + make_float3 (pSH8_x_4) * *_S346 + make_float3 (fTmp0B_7) * *_S345);
    float v_y_2 = _S369 + dot_0(_S364, make_float3 (pSH8_x_4) * *_S342 + make_float3 (0.54627424478530884f * fC1_y_2) * *_S346 + make_float3 (fTmp0B_7) * *_S343);
    float v_z_2 = _S370 + dot_0(_S364, make_float3 (pSH6_z_0) * *_S344 + make_float3 (-1.09254848957061768f * x_12) * *_S345 + make_float3 (-1.09254848957061768f * y_11) * *_S343);
    float3  * _S376 = v_coeffs_8 + int(8);
    *_S376 = *_S376 + make_float3 (pSH9_2) * _S364;
    float3  * _S377 = v_coeffs_8 + int(9);
    *_S377 = *_S377 + make_float3 (pSH10_2) * _S364;
    float3  * _S378 = v_coeffs_8 + int(10);
    *_S378 = *_S378 + make_float3 (pSH11_2) * _S364;
    float3  * _S379 = v_coeffs_8 + int(11);
    *_S379 = *_S379 + make_float3 (pSH12_3) * _S364;
    float3  * _S380 = v_coeffs_8 + int(12);
    *_S380 = *_S380 + make_float3 (pSH13_2) * _S364;
    float3  * _S381 = v_coeffs_8 + int(13);
    *_S381 = *_S381 + make_float3 (pSH14_2) * _S364;
    float3  * _S382 = v_coeffs_8 + int(14);
    *_S382 = *_S382 + make_float3 (pSH15_2) * _S364;
    float fTmp0C_z_2 = -4.57045793533325195f * z_9;
    float _S383 = x_12 * _S341;
    float fC2_x_0 = fC1_4 + _S383 - y_11 * fS1_x_2;
    float _S384 = y_11 * _S341;
    float fC2_y_0 = x_12 * fC1_y_2 - fS1_4 - _S384;
    float fS2_x_0 = fS1_4 + x_12 * fS1_x_2 + _S384;
    float fS2_y_0 = _S383 + fC1_4 + y_11 * fC1_y_2;
    float pSH12_z_0 = 5.59764480590820312f * z2_4 - 1.11952900886535645f;
    float pSH14_x_2 = fTmp1B_4 * _S341;
    float v_x_3 = v_x_2 + dot_0(_S364, make_float3 (-0.59004360437393188f * fS2_x_0) * *_S347 + make_float3 (-0.59004360437393188f * fC2_x_0) * *_S353 + make_float3 (fTmp1B_4 * fS1_x_2) * *_S348 + make_float3 (pSH14_x_2) * *_S352 + make_float3 (fTmp0C_4) * *_S351);
    float v_y_3 = v_y_2 + dot_0(_S364, make_float3 (-0.59004360437393188f * fS2_y_0) * *_S347 + make_float3 (-0.59004360437393188f * fC2_y_0) * *_S353 + make_float3 (pSH14_x_2) * *_S348 + make_float3 (fTmp1B_4 * fC1_y_2) * *_S352 + make_float3 (fTmp0C_4) * *_S349);
    float v_z_3 = v_z_2 + dot_0(_S364, make_float3 (pSH12_z_0) * *_S350 + make_float3 (fTmp0C_z_2 * x_12) * *_S351 + make_float3 (fTmp0C_z_2 * y_11) * *_S349 + make_float3 (1.44530570507049561f * fC1_4) * *_S352 + make_float3 (1.44530570507049561f * fS1_4) * *_S348);
    float pSH20_0 = _S354 + -1.00623059272766113f * pSH6_5;
    float3  * _S385 = v_coeffs_8 + int(15);
    *_S385 = *_S385 + make_float3 (pSH16_0) * _S364;
    float3  * _S386 = v_coeffs_8 + int(16);
    *_S386 = *_S386 + make_float3 (pSH17_0) * _S364;
    float3  * _S387 = v_coeffs_8 + int(17);
    *_S387 = *_S387 + make_float3 (pSH18_0) * _S364;
    float3  * _S388 = v_coeffs_8 + int(18);
    *_S388 = *_S388 + make_float3 (pSH19_0) * _S364;
    float3  * _S389 = v_coeffs_8 + int(19);
    *_S389 = *_S389 + make_float3 (pSH20_0) * _S364;
    float3  * _S390 = v_coeffs_8 + int(20);
    *_S390 = *_S390 + make_float3 (pSH21_0) * _S364;
    float3  * _S391 = v_coeffs_8 + int(21);
    *_S391 = *_S391 + make_float3 (pSH22_0) * _S364;
    float3  * _S392 = v_coeffs_8 + int(22);
    *_S392 = *_S392 + make_float3 (pSH23_0) * _S364;
    float3  * _S393 = v_coeffs_8 + int(23);
    *_S393 = *_S393 + make_float3 (pSH24_0) * _S364;
    float fTmp0D_z_0 = -14.04997730255126953f * z2_4 + 2.00713968276977539f;
    float fTmp1C_z_0 = 6.62322282791137695f * z_9;
    float pSH22_x_0 = fTmp1C_1 * _S341;
    float3  dir_n_6 = make_float3 (x_12, y_11, z_9);
    float3  v_dir_n_6 = make_float3 (v_x_3 + dot_0(_S364, make_float3 (0.62583571672439575f * (fS2_1 + y_11 * fC2_x_0 + x_12 * fS2_x_0)) * *_S355 + make_float3 (0.62583571672439575f * (fC2_1 + x_12 * fC2_x_0 - y_11 * fS2_x_0)) * *_S363 + make_float3 (fTmp2B_1 * fS2_x_0) * *_S356 + make_float3 (fTmp2B_1 * fC2_x_0) * *_S362 + make_float3 (fTmp1C_1 * fS1_x_2) * *_S357 + make_float3 (pSH22_x_0) * *_S361 + make_float3 (fTmp0D_1) * *_S360), v_y_3 + dot_0(_S364, make_float3 (0.62583571672439575f * (x_12 * fS2_y_0 + fC2_1 + y_11 * fC2_y_0)) * *_S355 + make_float3 (0.62583571672439575f * (x_12 * fC2_y_0 - fS2_1 - y_11 * fS2_y_0)) * *_S363 + make_float3 (fTmp2B_1 * fS2_y_0) * *_S356 + make_float3 (fTmp2B_1 * fC2_y_0) * *_S362 + make_float3 (pSH22_x_0) * *_S357 + make_float3 (fTmp1C_1 * fC1_y_2) * *_S361 + make_float3 (fTmp0D_1) * *_S358), v_z_3 + dot_0(_S364, make_float3 (1.9843134880065918f * (pSH12_3 + z_9 * pSH12_z_0) + -1.00623059272766113f * pSH6_z_0) * *_S359 + make_float3 (fTmp0D_z_0 * x_12) * *_S360 + make_float3 (fTmp0D_z_0 * y_11) * *_S358 + make_float3 (fTmp1C_z_0 * fC1_4) * *_S361 + make_float3 (fTmp1C_z_0 * fS1_4) * *_S357 + make_float3 (-1.77013075351715088f * fC2_1) * *_S362 + make_float3 (-1.77013075351715088f * fS2_1) * *_S356));
    float3  v_viewdir_15 = v_viewdir_14 + (v_dir_n_6 - make_float3 (dot_0(v_dir_n_6, dir_n_6)) * dir_n_6) * make_float3 (inv_norm_10);
    Matrix<float, 3, 3>  _S394 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S395;
    (&_S395)->primal_0 = _S333;
    (&_S395)->differential_0 = _S394;
    float3  _S396 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S397;
    (&_S397)->primal_0 = t_13;
    (&_S397)->differential_0 = _S396;
    s_bwd_prop_mul_0(&_S395, &_S397, v_viewdir_15);
    Matrix<float, 3, 3>  _S398 = transpose_0(_S395.differential_0);
    *v_mean_8 = *v_mean_8 + v_viewdir_15;
    *v_R_8 = *v_R_8 + _S398;
    *v_t_8 = *v_t_8 + _S397.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_atomic(float3  mean_14, Matrix<float, 3, 3>  R_14, float3  t_14, float3  coeff_dc_14, float3  * coeffs_14, float3  v_colors_9, float3  * v_coeff_dc_9, float3  * v_coeffs_9, float3  * v_mean_9, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    Matrix<float, 3, 3>  _S399 = transpose_0(R_14);
    float3  _S400 = mean_14 + mul_0(_S399, t_14);
    float _S401 = _S400.x;
    float _S402 = _S400.y;
    float _S403 = _S400.z;
    float inv_norm_11 = (F32_rsqrt((_S401 * _S401 + _S402 * _S402 + _S403 * _S403)));
    float x_13 = _S401 * inv_norm_11;
    float y_12 = _S402 * inv_norm_11;
    float z_10 = _S403 * inv_norm_11;
    float3  * _S404 = coeffs_14 + int(0);
    float3  * _S405 = coeffs_14 + int(1);
    float3  * _S406 = coeffs_14 + int(2);
    float z2_5 = z_10 * z_10;
    float fTmp0B_8 = -1.09254848957061768f * z_10;
    float fC1_5 = x_13 * x_13 - y_12 * y_12;
    float _S407 = 2.0f * x_13;
    float fS1_5 = _S407 * y_12;
    float pSH6_6 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
    float pSH7_5 = fTmp0B_8 * x_13;
    float pSH5_5 = fTmp0B_8 * y_12;
    float pSH8_5 = 0.54627424478530884f * fC1_5;
    float pSH4_5 = 0.54627424478530884f * fS1_5;
    float3  * _S408 = coeffs_14 + int(3);
    float3  * _S409 = coeffs_14 + int(4);
    float3  * _S410 = coeffs_14 + int(5);
    float3  * _S411 = coeffs_14 + int(6);
    float3  * _S412 = coeffs_14 + int(7);
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
    float3  * _S413 = coeffs_14 + int(8);
    float3  * _S414 = coeffs_14 + int(9);
    float3  * _S415 = coeffs_14 + int(10);
    float3  * _S416 = coeffs_14 + int(11);
    float3  * _S417 = coeffs_14 + int(12);
    float3  * _S418 = coeffs_14 + int(13);
    float3  * _S419 = coeffs_14 + int(14);
    float fTmp0D_2 = z_10 * (-4.68332576751708984f * z2_5 + 2.00713968276977539f);
    float fTmp1C_2 = 3.31161141395568848f * z2_5 - 0.47308734059333801f;
    float fTmp2B_2 = -1.77013075351715088f * z_10;
    float _S420 = 1.9843134880065918f * z_10 * pSH12_4;
    float pSH21_1 = fTmp0D_2 * x_13;
    float pSH19_1 = fTmp0D_2 * y_12;
    float pSH22_1 = fTmp1C_2 * fC1_5;
    float pSH18_1 = fTmp1C_2 * fS1_5;
    float pSH23_1 = fTmp2B_2 * fC2_2;
    float pSH17_1 = fTmp2B_2 * fS2_2;
    float pSH24_1 = 0.62583571672439575f * (x_13 * fC2_2 - y_12 * fS2_2);
    float pSH16_1 = 0.62583571672439575f * (x_13 * fS2_2 + y_12 * fC2_2);
    float3  * _S421 = coeffs_14 + int(15);
    float3  * _S422 = coeffs_14 + int(16);
    float3  * _S423 = coeffs_14 + int(17);
    float3  * _S424 = coeffs_14 + int(18);
    float3  * _S425 = coeffs_14 + int(19);
    float3  * _S426 = coeffs_14 + int(20);
    float3  * _S427 = coeffs_14 + int(21);
    float3  * _S428 = coeffs_14 + int(22);
    float3  * _S429 = coeffs_14 + int(23);
    float3  colors_9 = make_float3 (0.282094806432724f) * coeff_dc_14 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * *_S404 + make_float3 (z_10) * *_S405 - make_float3 (x_13) * *_S406) + (make_float3 (pSH4_5) * *_S408 + make_float3 (pSH5_5) * *_S409 + make_float3 (pSH6_6) * *_S410 + make_float3 (pSH7_5) * *_S411 + make_float3 (pSH8_5) * *_S412) + (make_float3 (pSH9_3) * *_S413 + make_float3 (pSH10_3) * *_S414 + make_float3 (pSH11_3) * *_S415 + make_float3 (pSH12_4) * *_S416 + make_float3 (pSH13_3) * *_S417 + make_float3 (pSH14_3) * *_S418 + make_float3 (pSH15_3) * *_S419) + (make_float3 (pSH16_1) * *_S421 + make_float3 (pSH17_1) * *_S422 + make_float3 (pSH18_1) * *_S423 + make_float3 (pSH19_1) * *_S424 + make_float3 (_S420 - 1.00623059272766113f * pSH6_6) * *_S425 + make_float3 (pSH21_1) * *_S426 + make_float3 (pSH22_1) * *_S427 + make_float3 (pSH23_1) * *_S428 + make_float3 (pSH24_1) * *_S429);
    float3  _S430 = v_colors_9 * make_float3 (float((colors_9.x) >= -0.5f), float((colors_9.y) >= -0.5f), float((colors_9.z) >= -0.5f));
    float3  v_viewdir_16 = {};
    *v_coeff_dc_9 = *v_coeff_dc_9 + make_float3 (0.282094806432724f) * _S430;
    float3  temp_26 = make_float3 (-0.48860251903533936f * y_12) * _S430;
    float _S431 = dot_0(temp_26, temp_26);
    bool _S432;
    if((F32_isfinite((_S431))))
    {
        _S432 = _S431 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S433 = v_coeffs_9 + int(0);
        float _S434 = atomicAdd(&(_S433->x), temp_26.x);
        float _S435 = atomicAdd(&(_S433->y), temp_26.y);
        float _S436 = atomicAdd(&(_S433->z), temp_26.z);
    }
    float3  temp_27 = make_float3 (0.48860251903533936f * z_10) * _S430;
    float _S437 = dot_0(temp_27, temp_27);
    if((F32_isfinite((_S437))))
    {
        _S432 = _S437 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S438 = v_coeffs_9 + int(1);
        float _S439 = atomicAdd(&(_S438->x), temp_27.x);
        float _S440 = atomicAdd(&(_S438->y), temp_27.y);
        float _S441 = atomicAdd(&(_S438->z), temp_27.z);
    }
    float3  temp_28 = make_float3 (-0.48860251903533936f * x_13) * _S430;
    float _S442 = dot_0(temp_28, temp_28);
    if((F32_isfinite((_S442))))
    {
        _S432 = _S442 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S443 = v_coeffs_9 + int(2);
        float _S444 = atomicAdd(&(_S443->x), temp_28.x);
        float _S445 = atomicAdd(&(_S443->y), temp_28.y);
        float _S446 = atomicAdd(&(_S443->z), temp_28.z);
    }
    float _S447 = -0.48860251903533936f * dot_0(*_S406, _S430);
    float _S448 = -0.48860251903533936f * dot_0(*_S404, _S430);
    float _S449 = 0.48860251903533936f * dot_0(*_S405, _S430);
    float3  temp_29 = make_float3 (pSH4_5) * _S430;
    float _S450 = dot_0(temp_29, temp_29);
    if((F32_isfinite((_S450))))
    {
        _S432 = _S450 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S451 = v_coeffs_9 + int(3);
        float _S452 = atomicAdd(&(_S451->x), temp_29.x);
        float _S453 = atomicAdd(&(_S451->y), temp_29.y);
        float _S454 = atomicAdd(&(_S451->z), temp_29.z);
    }
    float3  temp_30 = make_float3 (pSH5_5) * _S430;
    float _S455 = dot_0(temp_30, temp_30);
    if((F32_isfinite((_S455))))
    {
        _S432 = _S455 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S456 = v_coeffs_9 + int(4);
        float _S457 = atomicAdd(&(_S456->x), temp_30.x);
        float _S458 = atomicAdd(&(_S456->y), temp_30.y);
        float _S459 = atomicAdd(&(_S456->z), temp_30.z);
    }
    float3  temp_31 = make_float3 (pSH6_6) * _S430;
    float _S460 = dot_0(temp_31, temp_31);
    if((F32_isfinite((_S460))))
    {
        _S432 = _S460 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S461 = v_coeffs_9 + int(5);
        float _S462 = atomicAdd(&(_S461->x), temp_31.x);
        float _S463 = atomicAdd(&(_S461->y), temp_31.y);
        float _S464 = atomicAdd(&(_S461->z), temp_31.z);
    }
    float3  temp_32 = make_float3 (pSH7_5) * _S430;
    float _S465 = dot_0(temp_32, temp_32);
    if((F32_isfinite((_S465))))
    {
        _S432 = _S465 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S466 = v_coeffs_9 + int(6);
        float _S467 = atomicAdd(&(_S466->x), temp_32.x);
        float _S468 = atomicAdd(&(_S466->y), temp_32.y);
        float _S469 = atomicAdd(&(_S466->z), temp_32.z);
    }
    float3  temp_33 = make_float3 (pSH8_5) * _S430;
    float _S470 = dot_0(temp_33, temp_33);
    if((F32_isfinite((_S470))))
    {
        _S432 = _S470 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S471 = v_coeffs_9 + int(7);
        float _S472 = atomicAdd(&(_S471->x), temp_33.x);
        float _S473 = atomicAdd(&(_S471->y), temp_33.y);
        float _S474 = atomicAdd(&(_S471->z), temp_33.z);
    }
    float fC1_y_3 = -2.0f * y_12;
    float fS1_x_3 = 2.0f * y_12;
    float pSH6_z_1 = 1.89234936237335205f * z_10;
    float pSH8_x_5 = 0.54627424478530884f * _S407;
    float v_x_4 = _S447 + dot_0(_S430, make_float3 (0.54627424478530884f * fS1_x_3) * *_S408 + make_float3 (pSH8_x_5) * *_S412 + make_float3 (fTmp0B_8) * *_S411);
    float v_y_4 = _S448 + dot_0(_S430, make_float3 (pSH8_x_5) * *_S408 + make_float3 (0.54627424478530884f * fC1_y_3) * *_S412 + make_float3 (fTmp0B_8) * *_S409);
    float v_z_4 = _S449 + dot_0(_S430, make_float3 (pSH6_z_1) * *_S410 + make_float3 (-1.09254848957061768f * x_13) * *_S411 + make_float3 (-1.09254848957061768f * y_12) * *_S409);
    float3  temp_34 = make_float3 (pSH9_3) * _S430;
    float _S475 = dot_0(temp_34, temp_34);
    if((F32_isfinite((_S475))))
    {
        _S432 = _S475 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S476 = v_coeffs_9 + int(8);
        float _S477 = atomicAdd(&(_S476->x), temp_34.x);
        float _S478 = atomicAdd(&(_S476->y), temp_34.y);
        float _S479 = atomicAdd(&(_S476->z), temp_34.z);
    }
    float3  temp_35 = make_float3 (pSH10_3) * _S430;
    float _S480 = dot_0(temp_35, temp_35);
    if((F32_isfinite((_S480))))
    {
        _S432 = _S480 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S481 = v_coeffs_9 + int(9);
        float _S482 = atomicAdd(&(_S481->x), temp_35.x);
        float _S483 = atomicAdd(&(_S481->y), temp_35.y);
        float _S484 = atomicAdd(&(_S481->z), temp_35.z);
    }
    float3  temp_36 = make_float3 (pSH11_3) * _S430;
    float _S485 = dot_0(temp_36, temp_36);
    if((F32_isfinite((_S485))))
    {
        _S432 = _S485 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S486 = v_coeffs_9 + int(10);
        float _S487 = atomicAdd(&(_S486->x), temp_36.x);
        float _S488 = atomicAdd(&(_S486->y), temp_36.y);
        float _S489 = atomicAdd(&(_S486->z), temp_36.z);
    }
    float3  temp_37 = make_float3 (pSH12_4) * _S430;
    float _S490 = dot_0(temp_37, temp_37);
    if((F32_isfinite((_S490))))
    {
        _S432 = _S490 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S491 = v_coeffs_9 + int(11);
        float _S492 = atomicAdd(&(_S491->x), temp_37.x);
        float _S493 = atomicAdd(&(_S491->y), temp_37.y);
        float _S494 = atomicAdd(&(_S491->z), temp_37.z);
    }
    float3  temp_38 = make_float3 (pSH13_3) * _S430;
    float _S495 = dot_0(temp_38, temp_38);
    if((F32_isfinite((_S495))))
    {
        _S432 = _S495 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S496 = v_coeffs_9 + int(12);
        float _S497 = atomicAdd(&(_S496->x), temp_38.x);
        float _S498 = atomicAdd(&(_S496->y), temp_38.y);
        float _S499 = atomicAdd(&(_S496->z), temp_38.z);
    }
    float3  temp_39 = make_float3 (pSH14_3) * _S430;
    float _S500 = dot_0(temp_39, temp_39);
    if((F32_isfinite((_S500))))
    {
        _S432 = _S500 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S501 = v_coeffs_9 + int(13);
        float _S502 = atomicAdd(&(_S501->x), temp_39.x);
        float _S503 = atomicAdd(&(_S501->y), temp_39.y);
        float _S504 = atomicAdd(&(_S501->z), temp_39.z);
    }
    float3  temp_40 = make_float3 (pSH15_3) * _S430;
    float _S505 = dot_0(temp_40, temp_40);
    if((F32_isfinite((_S505))))
    {
        _S432 = _S505 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S506 = v_coeffs_9 + int(14);
        float _S507 = atomicAdd(&(_S506->x), temp_40.x);
        float _S508 = atomicAdd(&(_S506->y), temp_40.y);
        float _S509 = atomicAdd(&(_S506->z), temp_40.z);
    }
    float fTmp0C_z_3 = -4.57045793533325195f * z_10;
    float _S510 = x_13 * _S407;
    float fC2_x_1 = fC1_5 + _S510 - y_12 * fS1_x_3;
    float _S511 = y_12 * _S407;
    float fC2_y_1 = x_13 * fC1_y_3 - fS1_5 - _S511;
    float fS2_x_1 = fS1_5 + x_13 * fS1_x_3 + _S511;
    float fS2_y_1 = _S510 + fC1_5 + y_12 * fC1_y_3;
    float pSH12_z_1 = 5.59764480590820312f * z2_5 - 1.11952900886535645f;
    float pSH14_x_3 = fTmp1B_5 * _S407;
    float v_x_5 = v_x_4 + dot_0(_S430, make_float3 (-0.59004360437393188f * fS2_x_1) * *_S413 + make_float3 (-0.59004360437393188f * fC2_x_1) * *_S419 + make_float3 (fTmp1B_5 * fS1_x_3) * *_S414 + make_float3 (pSH14_x_3) * *_S418 + make_float3 (fTmp0C_5) * *_S417);
    float v_y_5 = v_y_4 + dot_0(_S430, make_float3 (-0.59004360437393188f * fS2_y_1) * *_S413 + make_float3 (-0.59004360437393188f * fC2_y_1) * *_S419 + make_float3 (pSH14_x_3) * *_S414 + make_float3 (fTmp1B_5 * fC1_y_3) * *_S418 + make_float3 (fTmp0C_5) * *_S415);
    float v_z_5 = v_z_4 + dot_0(_S430, make_float3 (pSH12_z_1) * *_S416 + make_float3 (fTmp0C_z_3 * x_13) * *_S417 + make_float3 (fTmp0C_z_3 * y_12) * *_S415 + make_float3 (1.44530570507049561f * fC1_5) * *_S418 + make_float3 (1.44530570507049561f * fS1_5) * *_S414);
    float pSH20_1 = _S420 + -1.00623059272766113f * pSH6_6;
    float3  temp_41 = make_float3 (pSH16_1) * _S430;
    float _S512 = dot_0(temp_41, temp_41);
    if((F32_isfinite((_S512))))
    {
        _S432 = _S512 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S513 = v_coeffs_9 + int(15);
        float _S514 = atomicAdd(&(_S513->x), temp_41.x);
        float _S515 = atomicAdd(&(_S513->y), temp_41.y);
        float _S516 = atomicAdd(&(_S513->z), temp_41.z);
    }
    float3  temp_42 = make_float3 (pSH17_1) * _S430;
    float _S517 = dot_0(temp_42, temp_42);
    if((F32_isfinite((_S517))))
    {
        _S432 = _S517 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S518 = v_coeffs_9 + int(16);
        float _S519 = atomicAdd(&(_S518->x), temp_42.x);
        float _S520 = atomicAdd(&(_S518->y), temp_42.y);
        float _S521 = atomicAdd(&(_S518->z), temp_42.z);
    }
    float3  temp_43 = make_float3 (pSH18_1) * _S430;
    float _S522 = dot_0(temp_43, temp_43);
    if((F32_isfinite((_S522))))
    {
        _S432 = _S522 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S523 = v_coeffs_9 + int(17);
        float _S524 = atomicAdd(&(_S523->x), temp_43.x);
        float _S525 = atomicAdd(&(_S523->y), temp_43.y);
        float _S526 = atomicAdd(&(_S523->z), temp_43.z);
    }
    float3  temp_44 = make_float3 (pSH19_1) * _S430;
    float _S527 = dot_0(temp_44, temp_44);
    if((F32_isfinite((_S527))))
    {
        _S432 = _S527 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S528 = v_coeffs_9 + int(18);
        float _S529 = atomicAdd(&(_S528->x), temp_44.x);
        float _S530 = atomicAdd(&(_S528->y), temp_44.y);
        float _S531 = atomicAdd(&(_S528->z), temp_44.z);
    }
    float3  temp_45 = make_float3 (pSH20_1) * _S430;
    float _S532 = dot_0(temp_45, temp_45);
    if((F32_isfinite((_S532))))
    {
        _S432 = _S532 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S533 = v_coeffs_9 + int(19);
        float _S534 = atomicAdd(&(_S533->x), temp_45.x);
        float _S535 = atomicAdd(&(_S533->y), temp_45.y);
        float _S536 = atomicAdd(&(_S533->z), temp_45.z);
    }
    float3  temp_46 = make_float3 (pSH21_1) * _S430;
    float _S537 = dot_0(temp_46, temp_46);
    if((F32_isfinite((_S537))))
    {
        _S432 = _S537 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S538 = v_coeffs_9 + int(20);
        float _S539 = atomicAdd(&(_S538->x), temp_46.x);
        float _S540 = atomicAdd(&(_S538->y), temp_46.y);
        float _S541 = atomicAdd(&(_S538->z), temp_46.z);
    }
    float3  temp_47 = make_float3 (pSH22_1) * _S430;
    float _S542 = dot_0(temp_47, temp_47);
    if((F32_isfinite((_S542))))
    {
        _S432 = _S542 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S543 = v_coeffs_9 + int(21);
        float _S544 = atomicAdd(&(_S543->x), temp_47.x);
        float _S545 = atomicAdd(&(_S543->y), temp_47.y);
        float _S546 = atomicAdd(&(_S543->z), temp_47.z);
    }
    float3  temp_48 = make_float3 (pSH23_1) * _S430;
    float _S547 = dot_0(temp_48, temp_48);
    if((F32_isfinite((_S547))))
    {
        _S432 = _S547 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S548 = v_coeffs_9 + int(22);
        float _S549 = atomicAdd(&(_S548->x), temp_48.x);
        float _S550 = atomicAdd(&(_S548->y), temp_48.y);
        float _S551 = atomicAdd(&(_S548->z), temp_48.z);
    }
    float3  temp_49 = make_float3 (pSH24_1) * _S430;
    float _S552 = dot_0(temp_49, temp_49);
    if((F32_isfinite((_S552))))
    {
        _S432 = _S552 != 0.0f;
    }
    else
    {
        _S432 = false;
    }
    if(_S432)
    {
        float3  * _S553 = v_coeffs_9 + int(23);
        float _S554 = atomicAdd(&(_S553->x), temp_49.x);
        float _S555 = atomicAdd(&(_S553->y), temp_49.y);
        float _S556 = atomicAdd(&(_S553->z), temp_49.z);
    }
    float fTmp0D_z_1 = -14.04997730255126953f * z2_5 + 2.00713968276977539f;
    float fTmp1C_z_1 = 6.62322282791137695f * z_10;
    float pSH22_x_1 = fTmp1C_2 * _S407;
    float3  dir_n_7 = make_float3 (x_13, y_12, z_10);
    float3  v_dir_n_7 = make_float3 (v_x_5 + dot_0(_S430, make_float3 (0.62583571672439575f * (fS2_2 + y_12 * fC2_x_1 + x_13 * fS2_x_1)) * *_S421 + make_float3 (0.62583571672439575f * (fC2_2 + x_13 * fC2_x_1 - y_12 * fS2_x_1)) * *_S429 + make_float3 (fTmp2B_2 * fS2_x_1) * *_S422 + make_float3 (fTmp2B_2 * fC2_x_1) * *_S428 + make_float3 (fTmp1C_2 * fS1_x_3) * *_S423 + make_float3 (pSH22_x_1) * *_S427 + make_float3 (fTmp0D_2) * *_S426), v_y_5 + dot_0(_S430, make_float3 (0.62583571672439575f * (x_13 * fS2_y_1 + fC2_2 + y_12 * fC2_y_1)) * *_S421 + make_float3 (0.62583571672439575f * (x_13 * fC2_y_1 - fS2_2 - y_12 * fS2_y_1)) * *_S429 + make_float3 (fTmp2B_2 * fS2_y_1) * *_S422 + make_float3 (fTmp2B_2 * fC2_y_1) * *_S428 + make_float3 (pSH22_x_1) * *_S423 + make_float3 (fTmp1C_2 * fC1_y_3) * *_S427 + make_float3 (fTmp0D_2) * *_S424), v_z_5 + dot_0(_S430, make_float3 (1.9843134880065918f * (pSH12_4 + z_10 * pSH12_z_1) + -1.00623059272766113f * pSH6_z_1) * *_S425 + make_float3 (fTmp0D_z_1 * x_13) * *_S426 + make_float3 (fTmp0D_z_1 * y_12) * *_S424 + make_float3 (fTmp1C_z_1 * fC1_5) * *_S427 + make_float3 (fTmp1C_z_1 * fS1_5) * *_S423 + make_float3 (-1.77013075351715088f * fC2_2) * *_S428 + make_float3 (-1.77013075351715088f * fS2_2) * *_S422));
    float3  v_viewdir_17 = v_viewdir_16 + (v_dir_n_7 - make_float3 (dot_0(v_dir_n_7, dir_n_7)) * dir_n_7) * make_float3 (inv_norm_11);
    Matrix<float, 3, 3>  _S557 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S558;
    (&_S558)->primal_0 = _S399;
    (&_S558)->differential_0 = _S557;
    float3  _S559 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S560;
    (&_S560)->primal_0 = t_14;
    (&_S560)->differential_0 = _S559;
    s_bwd_prop_mul_0(&_S558, &_S560, v_viewdir_17);
    Matrix<float, 3, 3>  _S561 = transpose_0(_S558.differential_0);
    *v_mean_9 = *v_mean_9 + v_viewdir_17;
    *v_R_9 = *v_R_9 + _S561;
    *v_t_9 = *v_t_9 + _S560.differential_0;
    return;
}

inline __device__ float3  sh0_to_color(float3  mean_15, Matrix<float, 3, 3>  R_15, float3  t_15, float3  coeff_dc_15, float * coeffs_15)
{
    return max_0(make_float3 (0.282094806432724f * coeff_dc_15.x, 0.282094806432724f * coeff_dc_15.y, 0.282094806432724f * coeff_dc_15.z) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh0_to_color_vjp_inplace(float3  mean_16, Matrix<float, 3, 3>  R_16, float3  t_16, float3  coeff_dc_16, float * coeffs_16, float3  v_colors_10, float3  * v_coeff_dc_10, float * v_coeffs_10, float3  * v_mean_10, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_10)
{
    float3  _S562 = v_colors_10 * make_float3 (float((0.282094806432724f * coeff_dc_16.x) >= -0.5f), float((0.282094806432724f * coeff_dc_16.y) >= -0.5f), float((0.282094806432724f * coeff_dc_16.z) >= -0.5f));
    float3  v_viewdir_18 = {};
    *&(v_coeff_dc_10->x) = *&(v_coeff_dc_10->x) + 0.282094806432724f * _S562.x;
    *&(v_coeff_dc_10->y) = *&(v_coeff_dc_10->y) + 0.282094806432724f * _S562.y;
    *&(v_coeff_dc_10->z) = *&(v_coeff_dc_10->z) + 0.282094806432724f * _S562.z;
    Matrix<float, 3, 3>  _S563 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S564;
    (&_S564)->primal_0 = transpose_0(R_16);
    (&_S564)->differential_0 = _S563;
    float3  _S565 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S566;
    (&_S566)->primal_0 = t_16;
    (&_S566)->differential_0 = _S565;
    s_bwd_prop_mul_0(&_S564, &_S566, v_viewdir_18);
    Matrix<float, 3, 3>  _S567 = transpose_0(_S564.differential_0);
    *v_mean_10 = *v_mean_10 + v_viewdir_18;
    *v_R_10 = *v_R_10 + _S567;
    *v_t_10 = *v_t_10 + _S566.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_atomic(float3  mean_17, Matrix<float, 3, 3>  R_17, float3  t_17, float3  coeff_dc_17, float * coeffs_17, float3  v_colors_11, float3  * v_coeff_dc_11, float * v_coeffs_11, float3  * v_mean_11, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_11)
{
    float3  _S568 = v_colors_11 * make_float3 (float((0.282094806432724f * coeff_dc_17.x) >= -0.5f), float((0.282094806432724f * coeff_dc_17.y) >= -0.5f), float((0.282094806432724f * coeff_dc_17.z) >= -0.5f));
    float3  v_viewdir_19 = {};
    *&(v_coeff_dc_11->x) = *&(v_coeff_dc_11->x) + 0.282094806432724f * _S568.x;
    *&(v_coeff_dc_11->y) = *&(v_coeff_dc_11->y) + 0.282094806432724f * _S568.y;
    *&(v_coeff_dc_11->z) = *&(v_coeff_dc_11->z) + 0.282094806432724f * _S568.z;
    Matrix<float, 3, 3>  _S569 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S570;
    (&_S570)->primal_0 = transpose_0(R_17);
    (&_S570)->differential_0 = _S569;
    float3  _S571 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S572;
    (&_S572)->primal_0 = t_17;
    (&_S572)->differential_0 = _S571;
    s_bwd_prop_mul_0(&_S570, &_S572, v_viewdir_19);
    Matrix<float, 3, 3>  _S573 = transpose_0(_S570.differential_0);
    *v_mean_11 = *v_mean_11 + v_viewdir_19;
    *v_R_11 = *v_R_11 + _S573;
    *v_t_11 = *v_t_11 + _S572.differential_0;
    return;
}

inline __device__ float3  sh1_to_color(float3  mean_18, Matrix<float, 3, 3>  R_18, float3  t_18, float3  coeff_dc_18, float * coeffs_18)
{
    float3  _S574 = mean_18 + mul_0(transpose_0(R_18), t_18);
    float _S575 = _S574.x;
    float _S576 = _S574.y;
    float _S577 = _S574.z;
    float inv_norm_12 = (F32_rsqrt((_S575 * _S575 + _S576 * _S576 + _S577 * _S577)));
    float x_14 = _S575 * inv_norm_12;
    float z_11 = _S577 * inv_norm_12;
    float _S578 = - (_S576 * inv_norm_12);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_18.x + 0.48860251903533936f * (_S578 * *(coeffs_18 + int(0)) + z_11 * *(coeffs_18 + int(3)) - x_14 * *(coeffs_18 + int(6))), 0.282094806432724f * coeff_dc_18.y + 0.48860251903533936f * (_S578 * *(coeffs_18 + int(1)) + z_11 * *(coeffs_18 + int(4)) - x_14 * *(coeffs_18 + int(7))), 0.282094806432724f * coeff_dc_18.z + 0.48860251903533936f * (_S578 * *(coeffs_18 + int(2)) + z_11 * *(coeffs_18 + int(5)) - x_14 * *(coeffs_18 + int(8)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh1_to_color_vjp_inplace(float3  mean_19, Matrix<float, 3, 3>  R_19, float3  t_19, float3  coeff_dc_19, float * coeffs_19, float3  v_colors_12, float3  * v_coeff_dc_12, float * v_coeffs_12, float3  * v_mean_12, Matrix<float, 3, 3>  * v_R_12, float3  * v_t_12)
{
    Matrix<float, 3, 3>  _S579 = transpose_0(R_19);
    float3  _S580 = mean_19 + mul_0(_S579, t_19);
    float _S581 = _S580.x;
    float _S582 = _S580.y;
    float _S583 = _S580.z;
    float inv_norm_13 = (F32_rsqrt((_S581 * _S581 + _S582 * _S582 + _S583 * _S583)));
    float x_15 = _S581 * inv_norm_13;
    float y_13 = _S582 * inv_norm_13;
    float z_12 = _S583 * inv_norm_13;
    float _S584 = - y_13;
    float * _S585 = coeffs_19 + int(0);
    float * _S586 = coeffs_19 + int(3);
    float * _S587 = coeffs_19 + int(6);
    float * _S588 = coeffs_19 + int(1);
    float * _S589 = coeffs_19 + int(4);
    float * _S590 = coeffs_19 + int(7);
    float * _S591 = coeffs_19 + int(2);
    float * _S592 = coeffs_19 + int(5);
    float * _S593 = coeffs_19 + int(8);
    float3  _S594 = v_colors_12 * make_float3 (float((0.282094806432724f * coeff_dc_19.x + 0.48860251903533936f * (_S584 * *_S585 + z_12 * *_S586 - x_15 * *_S587)) >= -0.5f), float((0.282094806432724f * coeff_dc_19.y + 0.48860251903533936f * (_S584 * *_S588 + z_12 * *_S589 - x_15 * *_S590)) >= -0.5f), float((0.282094806432724f * coeff_dc_19.z + 0.48860251903533936f * (_S584 * *_S591 + z_12 * *_S592 - x_15 * *_S593)) >= -0.5f));
    float3  v_viewdir_20 = {};
    float _S595 = _S594.x;
    *&(v_coeff_dc_12->x) = *&(v_coeff_dc_12->x) + 0.282094806432724f * _S595;
    float * _S596 = v_coeffs_12 + int(0);
    float _S597 = -0.48860251903533936f * y_13;
    *_S596 = *_S596 + _S597 * _S595;
    float * _S598 = v_coeffs_12 + int(3);
    float _S599 = 0.48860251903533936f * z_12;
    *_S598 = *_S598 + _S599 * _S595;
    float * _S600 = v_coeffs_12 + int(6);
    float _S601 = -0.48860251903533936f * x_15;
    *_S600 = *_S600 + _S601 * _S595;
    float3  dir_n_8 = make_float3 (x_15, y_13, z_12);
    float3  v_dir_n_8 = make_float3 (-0.48860251903533936f * *_S587 * _S595, -0.48860251903533936f * *_S585 * _S595, 0.48860251903533936f * *_S586 * _S595);
    float3  v_viewdir_21 = v_viewdir_20 + (v_dir_n_8 - make_float3 (dot_0(v_dir_n_8, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    float _S602 = _S594.y;
    *&(v_coeff_dc_12->y) = *&(v_coeff_dc_12->y) + 0.282094806432724f * _S602;
    float * _S603 = v_coeffs_12 + int(1);
    *_S603 = *_S603 + _S597 * _S602;
    float * _S604 = v_coeffs_12 + int(4);
    *_S604 = *_S604 + _S599 * _S602;
    float * _S605 = v_coeffs_12 + int(7);
    *_S605 = *_S605 + _S601 * _S602;
    float3  v_dir_n_9 = make_float3 (-0.48860251903533936f * *_S590 * _S602, -0.48860251903533936f * *_S588 * _S602, 0.48860251903533936f * *_S589 * _S602);
    float3  v_viewdir_22 = v_viewdir_21 + (v_dir_n_9 - make_float3 (dot_0(v_dir_n_9, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    float _S606 = _S594.z;
    *&(v_coeff_dc_12->z) = *&(v_coeff_dc_12->z) + 0.282094806432724f * _S606;
    float * _S607 = v_coeffs_12 + int(2);
    *_S607 = *_S607 + _S597 * _S606;
    float * _S608 = v_coeffs_12 + int(5);
    *_S608 = *_S608 + _S599 * _S606;
    float * _S609 = v_coeffs_12 + int(8);
    *_S609 = *_S609 + _S601 * _S606;
    float3  v_dir_n_10 = make_float3 (-0.48860251903533936f * *_S593 * _S606, -0.48860251903533936f * *_S591 * _S606, 0.48860251903533936f * *_S592 * _S606);
    float3  v_viewdir_23 = v_viewdir_22 + (v_dir_n_10 - make_float3 (dot_0(v_dir_n_10, dir_n_8)) * dir_n_8) * make_float3 (inv_norm_13);
    Matrix<float, 3, 3>  _S610 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S611;
    (&_S611)->primal_0 = _S579;
    (&_S611)->differential_0 = _S610;
    float3  _S612 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S613;
    (&_S613)->primal_0 = t_19;
    (&_S613)->differential_0 = _S612;
    s_bwd_prop_mul_0(&_S611, &_S613, v_viewdir_23);
    Matrix<float, 3, 3>  _S614 = transpose_0(_S611.differential_0);
    *v_mean_12 = *v_mean_12 + v_viewdir_23;
    *v_R_12 = *v_R_12 + _S614;
    *v_t_12 = *v_t_12 + _S613.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_atomic(float3  mean_20, Matrix<float, 3, 3>  R_20, float3  t_20, float3  coeff_dc_20, float * coeffs_20, float3  v_colors_13, float3  * v_coeff_dc_13, float * v_coeffs_13, float3  * v_mean_13, Matrix<float, 3, 3>  * v_R_13, float3  * v_t_13)
{
    Matrix<float, 3, 3>  _S615 = transpose_0(R_20);
    float3  _S616 = mean_20 + mul_0(_S615, t_20);
    float _S617 = _S616.x;
    float _S618 = _S616.y;
    float _S619 = _S616.z;
    float inv_norm_14 = (F32_rsqrt((_S617 * _S617 + _S618 * _S618 + _S619 * _S619)));
    float x_16 = _S617 * inv_norm_14;
    float y_14 = _S618 * inv_norm_14;
    float z_13 = _S619 * inv_norm_14;
    float _S620 = - y_14;
    float * _S621 = coeffs_20 + int(0);
    float * _S622 = coeffs_20 + int(3);
    float * _S623 = coeffs_20 + int(6);
    float * _S624 = coeffs_20 + int(1);
    float * _S625 = coeffs_20 + int(4);
    float * _S626 = coeffs_20 + int(7);
    float * _S627 = coeffs_20 + int(2);
    float * _S628 = coeffs_20 + int(5);
    float * _S629 = coeffs_20 + int(8);
    float3  _S630 = v_colors_13 * make_float3 (float((0.282094806432724f * coeff_dc_20.x + 0.48860251903533936f * (_S620 * *_S621 + z_13 * *_S622 - x_16 * *_S623)) >= -0.5f), float((0.282094806432724f * coeff_dc_20.y + 0.48860251903533936f * (_S620 * *_S624 + z_13 * *_S625 - x_16 * *_S626)) >= -0.5f), float((0.282094806432724f * coeff_dc_20.z + 0.48860251903533936f * (_S620 * *_S627 + z_13 * *_S628 - x_16 * *_S629)) >= -0.5f));
    float3  v_viewdir_24 = {};
    float _S631 = _S630.x;
    *&(v_coeff_dc_13->x) = *&(v_coeff_dc_13->x) + 0.282094806432724f * _S631;
    float _S632 = -0.48860251903533936f * y_14;
    float temp_50 = _S632 * _S631;
    bool _S633;
    if((F32_isfinite((temp_50))))
    {
        _S633 = temp_50 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S634 = atomicAdd(v_coeffs_13 + int(0), temp_50);
    }
    float _S635 = 0.48860251903533936f * z_13;
    float temp_51 = _S635 * _S631;
    if((F32_isfinite((temp_51))))
    {
        _S633 = temp_51 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S636 = atomicAdd(v_coeffs_13 + int(3), temp_51);
    }
    float _S637 = -0.48860251903533936f * x_16;
    float temp_52 = _S637 * _S631;
    if((F32_isfinite((temp_52))))
    {
        _S633 = temp_52 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S638 = atomicAdd(v_coeffs_13 + int(6), temp_52);
    }
    float3  dir_n_9 = make_float3 (x_16, y_14, z_13);
    float3  v_dir_n_11 = make_float3 (-0.48860251903533936f * *_S623 * _S631, -0.48860251903533936f * *_S621 * _S631, 0.48860251903533936f * *_S622 * _S631);
    float3  v_viewdir_25 = v_viewdir_24 + (v_dir_n_11 - make_float3 (dot_0(v_dir_n_11, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    float _S639 = _S630.y;
    *&(v_coeff_dc_13->y) = *&(v_coeff_dc_13->y) + 0.282094806432724f * _S639;
    float temp_53 = _S632 * _S639;
    if((F32_isfinite((temp_53))))
    {
        _S633 = temp_53 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S640 = atomicAdd(v_coeffs_13 + int(1), temp_53);
    }
    float temp_54 = _S635 * _S639;
    if((F32_isfinite((temp_54))))
    {
        _S633 = temp_54 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S641 = atomicAdd(v_coeffs_13 + int(4), temp_54);
    }
    float temp_55 = _S637 * _S639;
    if((F32_isfinite((temp_55))))
    {
        _S633 = temp_55 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S642 = atomicAdd(v_coeffs_13 + int(7), temp_55);
    }
    float3  v_dir_n_12 = make_float3 (-0.48860251903533936f * *_S626 * _S639, -0.48860251903533936f * *_S624 * _S639, 0.48860251903533936f * *_S625 * _S639);
    float3  v_viewdir_26 = v_viewdir_25 + (v_dir_n_12 - make_float3 (dot_0(v_dir_n_12, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    float _S643 = _S630.z;
    *&(v_coeff_dc_13->z) = *&(v_coeff_dc_13->z) + 0.282094806432724f * _S643;
    float temp_56 = _S632 * _S643;
    if((F32_isfinite((temp_56))))
    {
        _S633 = temp_56 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S644 = atomicAdd(v_coeffs_13 + int(2), temp_56);
    }
    float temp_57 = _S635 * _S643;
    if((F32_isfinite((temp_57))))
    {
        _S633 = temp_57 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S645 = atomicAdd(v_coeffs_13 + int(5), temp_57);
    }
    float temp_58 = _S637 * _S643;
    if((F32_isfinite((temp_58))))
    {
        _S633 = temp_58 != 0.0f;
    }
    else
    {
        _S633 = false;
    }
    if(_S633)
    {
        float _S646 = atomicAdd(v_coeffs_13 + int(8), temp_58);
    }
    float3  v_dir_n_13 = make_float3 (-0.48860251903533936f * *_S629 * _S643, -0.48860251903533936f * *_S627 * _S643, 0.48860251903533936f * *_S628 * _S643);
    float3  v_viewdir_27 = v_viewdir_26 + (v_dir_n_13 - make_float3 (dot_0(v_dir_n_13, dir_n_9)) * dir_n_9) * make_float3 (inv_norm_14);
    Matrix<float, 3, 3>  _S647 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S648;
    (&_S648)->primal_0 = _S615;
    (&_S648)->differential_0 = _S647;
    float3  _S649 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S650;
    (&_S650)->primal_0 = t_20;
    (&_S650)->differential_0 = _S649;
    s_bwd_prop_mul_0(&_S648, &_S650, v_viewdir_27);
    Matrix<float, 3, 3>  _S651 = transpose_0(_S648.differential_0);
    *v_mean_13 = *v_mean_13 + v_viewdir_27;
    *v_R_13 = *v_R_13 + _S651;
    *v_t_13 = *v_t_13 + _S650.differential_0;
    return;
}

inline __device__ float3  sh2_to_color(float3  mean_21, Matrix<float, 3, 3>  R_21, float3  t_21, float3  coeff_dc_21, float * coeffs_21)
{
    float3  _S652 = mean_21 + mul_0(transpose_0(R_21), t_21);
    float _S653 = _S652.x;
    float _S654 = _S652.y;
    float _S655 = _S652.z;
    float inv_norm_15 = (F32_rsqrt((_S653 * _S653 + _S654 * _S654 + _S655 * _S655)));
    float x_17 = _S653 * inv_norm_15;
    float y_15 = _S654 * inv_norm_15;
    float z_14 = _S655 * inv_norm_15;
    float _S656 = - y_15;
    float fTmp0B_9 = -1.09254848957061768f * z_14;
    float pSH6_7 = 0.94617468118667603f * (z_14 * z_14) - 0.31539157032966614f;
    float pSH7_6 = fTmp0B_9 * x_17;
    float pSH5_6 = fTmp0B_9 * y_15;
    float pSH8_6 = 0.54627424478530884f * (x_17 * x_17 - y_15 * y_15);
    float pSH4_6 = 0.54627424478530884f * (2.0f * x_17 * y_15);
    return max_0(make_float3 (0.282094806432724f * coeff_dc_21.x + 0.48860251903533936f * (_S656 * *(coeffs_21 + int(0)) + z_14 * *(coeffs_21 + int(3)) - x_17 * *(coeffs_21 + int(6))) + (pSH4_6 * *(coeffs_21 + int(9)) + pSH5_6 * *(coeffs_21 + int(12)) + pSH6_7 * *(coeffs_21 + int(15)) + pSH7_6 * *(coeffs_21 + int(18)) + pSH8_6 * *(coeffs_21 + int(21))), 0.282094806432724f * coeff_dc_21.y + 0.48860251903533936f * (_S656 * *(coeffs_21 + int(1)) + z_14 * *(coeffs_21 + int(4)) - x_17 * *(coeffs_21 + int(7))) + (pSH4_6 * *(coeffs_21 + int(10)) + pSH5_6 * *(coeffs_21 + int(13)) + pSH6_7 * *(coeffs_21 + int(16)) + pSH7_6 * *(coeffs_21 + int(19)) + pSH8_6 * *(coeffs_21 + int(22))), 0.282094806432724f * coeff_dc_21.z + 0.48860251903533936f * (_S656 * *(coeffs_21 + int(2)) + z_14 * *(coeffs_21 + int(5)) - x_17 * *(coeffs_21 + int(8))) + (pSH4_6 * *(coeffs_21 + int(11)) + pSH5_6 * *(coeffs_21 + int(14)) + pSH6_7 * *(coeffs_21 + int(17)) + pSH7_6 * *(coeffs_21 + int(20)) + pSH8_6 * *(coeffs_21 + int(23)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh2_to_color_vjp_inplace(float3  mean_22, Matrix<float, 3, 3>  R_22, float3  t_22, float3  coeff_dc_22, float * coeffs_22, float3  v_colors_14, float3  * v_coeff_dc_14, float * v_coeffs_14, float3  * v_mean_14, Matrix<float, 3, 3>  * v_R_14, float3  * v_t_14)
{
    Matrix<float, 3, 3>  _S657 = transpose_0(R_22);
    float3  _S658 = mean_22 + mul_0(_S657, t_22);
    float _S659 = _S658.x;
    float _S660 = _S658.y;
    float _S661 = _S658.z;
    float inv_norm_16 = (F32_rsqrt((_S659 * _S659 + _S660 * _S660 + _S661 * _S661)));
    float x_18 = _S659 * inv_norm_16;
    float y_16 = _S660 * inv_norm_16;
    float z_15 = _S661 * inv_norm_16;
    float _S662 = - y_16;
    float * _S663 = coeffs_22 + int(0);
    float * _S664 = coeffs_22 + int(3);
    float * _S665 = coeffs_22 + int(6);
    float fTmp0B_10 = -1.09254848957061768f * z_15;
    float _S666 = 2.0f * x_18;
    float pSH6_8 = 0.94617468118667603f * (z_15 * z_15) - 0.31539157032966614f;
    float pSH7_7 = fTmp0B_10 * x_18;
    float pSH5_7 = fTmp0B_10 * y_16;
    float pSH8_7 = 0.54627424478530884f * (x_18 * x_18 - y_16 * y_16);
    float pSH4_7 = 0.54627424478530884f * (_S666 * y_16);
    float * _S667 = coeffs_22 + int(9);
    float * _S668 = coeffs_22 + int(12);
    float * _S669 = coeffs_22 + int(15);
    float * _S670 = coeffs_22 + int(18);
    float * _S671 = coeffs_22 + int(21);
    float * _S672 = coeffs_22 + int(1);
    float * _S673 = coeffs_22 + int(4);
    float * _S674 = coeffs_22 + int(7);
    float * _S675 = coeffs_22 + int(10);
    float * _S676 = coeffs_22 + int(13);
    float * _S677 = coeffs_22 + int(16);
    float * _S678 = coeffs_22 + int(19);
    float * _S679 = coeffs_22 + int(22);
    float * _S680 = coeffs_22 + int(2);
    float * _S681 = coeffs_22 + int(5);
    float * _S682 = coeffs_22 + int(8);
    float * _S683 = coeffs_22 + int(11);
    float * _S684 = coeffs_22 + int(14);
    float * _S685 = coeffs_22 + int(17);
    float * _S686 = coeffs_22 + int(20);
    float * _S687 = coeffs_22 + int(23);
    float3  _S688 = v_colors_14 * make_float3 (float((0.282094806432724f * coeff_dc_22.x + 0.48860251903533936f * (_S662 * *_S663 + z_15 * *_S664 - x_18 * *_S665) + (pSH4_7 * *_S667 + pSH5_7 * *_S668 + pSH6_8 * *_S669 + pSH7_7 * *_S670 + pSH8_7 * *_S671)) >= -0.5f), float((0.282094806432724f * coeff_dc_22.y + 0.48860251903533936f * (_S662 * *_S672 + z_15 * *_S673 - x_18 * *_S674) + (pSH4_7 * *_S675 + pSH5_7 * *_S676 + pSH6_8 * *_S677 + pSH7_7 * *_S678 + pSH8_7 * *_S679)) >= -0.5f), float((0.282094806432724f * coeff_dc_22.z + 0.48860251903533936f * (_S662 * *_S680 + z_15 * *_S681 - x_18 * *_S682) + (pSH4_7 * *_S683 + pSH5_7 * *_S684 + pSH6_8 * *_S685 + pSH7_7 * *_S686 + pSH8_7 * *_S687)) >= -0.5f));
    float3  v_viewdir_28 = {};
    float _S689 = _S688.x;
    *&(v_coeff_dc_14->x) = *&(v_coeff_dc_14->x) + 0.282094806432724f * _S689;
    float * _S690 = v_coeffs_14 + int(0);
    float _S691 = -0.48860251903533936f * y_16;
    *_S690 = *_S690 + _S691 * _S689;
    float * _S692 = v_coeffs_14 + int(3);
    float _S693 = 0.48860251903533936f * z_15;
    *_S692 = *_S692 + _S693 * _S689;
    float * _S694 = v_coeffs_14 + int(6);
    float _S695 = -0.48860251903533936f * x_18;
    *_S694 = *_S694 + _S695 * _S689;
    float _S696 = -0.48860251903533936f * *_S665 * _S689;
    float _S697 = -0.48860251903533936f * *_S663 * _S689;
    float _S698 = 0.48860251903533936f * *_S664 * _S689;
    float * _S699 = v_coeffs_14 + int(9);
    *_S699 = *_S699 + pSH4_7 * _S689;
    float * _S700 = v_coeffs_14 + int(12);
    *_S700 = *_S700 + pSH5_7 * _S689;
    float * _S701 = v_coeffs_14 + int(15);
    *_S701 = *_S701 + pSH6_8 * _S689;
    float * _S702 = v_coeffs_14 + int(18);
    *_S702 = *_S702 + pSH7_7 * _S689;
    float * _S703 = v_coeffs_14 + int(21);
    *_S703 = *_S703 + pSH8_7 * _S689;
    float pSH6_z_2 = 1.89234936237335205f * z_15;
    float pSH7_z_0 = -1.09254848957061768f * x_18;
    float pSH5_z_0 = -1.09254848957061768f * y_16;
    float pSH8_x_6 = 0.54627424478530884f * _S666;
    float pSH8_y_0 = 0.54627424478530884f * (-2.0f * y_16);
    float pSH4_x_0 = 0.54627424478530884f * (2.0f * y_16);
    float3  dir_n_10 = make_float3 (x_18, y_16, z_15);
    float3  v_dir_n_14 = make_float3 (_S696 + _S689 * (pSH4_x_0 * *_S667 + pSH8_x_6 * *_S671 + fTmp0B_10 * *_S670), _S697 + _S689 * (pSH8_x_6 * *_S667 + pSH8_y_0 * *_S671 + fTmp0B_10 * *_S668), _S698 + _S689 * (pSH6_z_2 * *_S669 + pSH7_z_0 * *_S670 + pSH5_z_0 * *_S668));
    float3  v_viewdir_29 = v_viewdir_28 + (v_dir_n_14 - make_float3 (dot_0(v_dir_n_14, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    float _S704 = _S688.y;
    *&(v_coeff_dc_14->y) = *&(v_coeff_dc_14->y) + 0.282094806432724f * _S704;
    float * _S705 = v_coeffs_14 + int(1);
    *_S705 = *_S705 + _S691 * _S704;
    float * _S706 = v_coeffs_14 + int(4);
    *_S706 = *_S706 + _S693 * _S704;
    float * _S707 = v_coeffs_14 + int(7);
    *_S707 = *_S707 + _S695 * _S704;
    float _S708 = -0.48860251903533936f * *_S674 * _S704;
    float _S709 = -0.48860251903533936f * *_S672 * _S704;
    float _S710 = 0.48860251903533936f * *_S673 * _S704;
    float * _S711 = v_coeffs_14 + int(10);
    *_S711 = *_S711 + pSH4_7 * _S704;
    float * _S712 = v_coeffs_14 + int(13);
    *_S712 = *_S712 + pSH5_7 * _S704;
    float * _S713 = v_coeffs_14 + int(16);
    *_S713 = *_S713 + pSH6_8 * _S704;
    float * _S714 = v_coeffs_14 + int(19);
    *_S714 = *_S714 + pSH7_7 * _S704;
    float * _S715 = v_coeffs_14 + int(22);
    *_S715 = *_S715 + pSH8_7 * _S704;
    float3  v_dir_n_15 = make_float3 (_S708 + _S704 * (pSH4_x_0 * *_S675 + pSH8_x_6 * *_S679 + fTmp0B_10 * *_S678), _S709 + _S704 * (pSH8_x_6 * *_S675 + pSH8_y_0 * *_S679 + fTmp0B_10 * *_S676), _S710 + _S704 * (pSH6_z_2 * *_S677 + pSH7_z_0 * *_S678 + pSH5_z_0 * *_S676));
    float3  v_viewdir_30 = v_viewdir_29 + (v_dir_n_15 - make_float3 (dot_0(v_dir_n_15, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    float _S716 = _S688.z;
    *&(v_coeff_dc_14->z) = *&(v_coeff_dc_14->z) + 0.282094806432724f * _S716;
    float * _S717 = v_coeffs_14 + int(2);
    *_S717 = *_S717 + _S691 * _S716;
    float * _S718 = v_coeffs_14 + int(5);
    *_S718 = *_S718 + _S693 * _S716;
    float * _S719 = v_coeffs_14 + int(8);
    *_S719 = *_S719 + _S695 * _S716;
    float _S720 = -0.48860251903533936f * *_S682 * _S716;
    float _S721 = -0.48860251903533936f * *_S680 * _S716;
    float _S722 = 0.48860251903533936f * *_S681 * _S716;
    float * _S723 = v_coeffs_14 + int(11);
    *_S723 = *_S723 + pSH4_7 * _S716;
    float * _S724 = v_coeffs_14 + int(14);
    *_S724 = *_S724 + pSH5_7 * _S716;
    float * _S725 = v_coeffs_14 + int(17);
    *_S725 = *_S725 + pSH6_8 * _S716;
    float * _S726 = v_coeffs_14 + int(20);
    *_S726 = *_S726 + pSH7_7 * _S716;
    float * _S727 = v_coeffs_14 + int(23);
    *_S727 = *_S727 + pSH8_7 * _S716;
    float3  v_dir_n_16 = make_float3 (_S720 + _S716 * (pSH4_x_0 * *_S683 + pSH8_x_6 * *_S687 + fTmp0B_10 * *_S686), _S721 + _S716 * (pSH8_x_6 * *_S683 + pSH8_y_0 * *_S687 + fTmp0B_10 * *_S684), _S722 + _S716 * (pSH6_z_2 * *_S685 + pSH7_z_0 * *_S686 + pSH5_z_0 * *_S684));
    float3  v_viewdir_31 = v_viewdir_30 + (v_dir_n_16 - make_float3 (dot_0(v_dir_n_16, dir_n_10)) * dir_n_10) * make_float3 (inv_norm_16);
    Matrix<float, 3, 3>  _S728 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S729;
    (&_S729)->primal_0 = _S657;
    (&_S729)->differential_0 = _S728;
    float3  _S730 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S731;
    (&_S731)->primal_0 = t_22;
    (&_S731)->differential_0 = _S730;
    s_bwd_prop_mul_0(&_S729, &_S731, v_viewdir_31);
    Matrix<float, 3, 3>  _S732 = transpose_0(_S729.differential_0);
    *v_mean_14 = *v_mean_14 + v_viewdir_31;
    *v_R_14 = *v_R_14 + _S732;
    *v_t_14 = *v_t_14 + _S731.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_atomic(float3  mean_23, Matrix<float, 3, 3>  R_23, float3  t_23, float3  coeff_dc_23, float * coeffs_23, float3  v_colors_15, float3  * v_coeff_dc_15, float * v_coeffs_15, float3  * v_mean_15, Matrix<float, 3, 3>  * v_R_15, float3  * v_t_15)
{
    Matrix<float, 3, 3>  _S733 = transpose_0(R_23);
    float3  _S734 = mean_23 + mul_0(_S733, t_23);
    float _S735 = _S734.x;
    float _S736 = _S734.y;
    float _S737 = _S734.z;
    float inv_norm_17 = (F32_rsqrt((_S735 * _S735 + _S736 * _S736 + _S737 * _S737)));
    float x_19 = _S735 * inv_norm_17;
    float y_17 = _S736 * inv_norm_17;
    float z_16 = _S737 * inv_norm_17;
    float _S738 = - y_17;
    float * _S739 = coeffs_23 + int(0);
    float * _S740 = coeffs_23 + int(3);
    float * _S741 = coeffs_23 + int(6);
    float fTmp0B_11 = -1.09254848957061768f * z_16;
    float _S742 = 2.0f * x_19;
    float pSH6_9 = 0.94617468118667603f * (z_16 * z_16) - 0.31539157032966614f;
    float pSH7_8 = fTmp0B_11 * x_19;
    float pSH5_8 = fTmp0B_11 * y_17;
    float pSH8_8 = 0.54627424478530884f * (x_19 * x_19 - y_17 * y_17);
    float pSH4_8 = 0.54627424478530884f * (_S742 * y_17);
    float * _S743 = coeffs_23 + int(9);
    float * _S744 = coeffs_23 + int(12);
    float * _S745 = coeffs_23 + int(15);
    float * _S746 = coeffs_23 + int(18);
    float * _S747 = coeffs_23 + int(21);
    float * _S748 = coeffs_23 + int(1);
    float * _S749 = coeffs_23 + int(4);
    float * _S750 = coeffs_23 + int(7);
    float * _S751 = coeffs_23 + int(10);
    float * _S752 = coeffs_23 + int(13);
    float * _S753 = coeffs_23 + int(16);
    float * _S754 = coeffs_23 + int(19);
    float * _S755 = coeffs_23 + int(22);
    float * _S756 = coeffs_23 + int(2);
    float * _S757 = coeffs_23 + int(5);
    float * _S758 = coeffs_23 + int(8);
    float * _S759 = coeffs_23 + int(11);
    float * _S760 = coeffs_23 + int(14);
    float * _S761 = coeffs_23 + int(17);
    float * _S762 = coeffs_23 + int(20);
    float * _S763 = coeffs_23 + int(23);
    float3  _S764 = v_colors_15 * make_float3 (float((0.282094806432724f * coeff_dc_23.x + 0.48860251903533936f * (_S738 * *_S739 + z_16 * *_S740 - x_19 * *_S741) + (pSH4_8 * *_S743 + pSH5_8 * *_S744 + pSH6_9 * *_S745 + pSH7_8 * *_S746 + pSH8_8 * *_S747)) >= -0.5f), float((0.282094806432724f * coeff_dc_23.y + 0.48860251903533936f * (_S738 * *_S748 + z_16 * *_S749 - x_19 * *_S750) + (pSH4_8 * *_S751 + pSH5_8 * *_S752 + pSH6_9 * *_S753 + pSH7_8 * *_S754 + pSH8_8 * *_S755)) >= -0.5f), float((0.282094806432724f * coeff_dc_23.z + 0.48860251903533936f * (_S738 * *_S756 + z_16 * *_S757 - x_19 * *_S758) + (pSH4_8 * *_S759 + pSH5_8 * *_S760 + pSH6_9 * *_S761 + pSH7_8 * *_S762 + pSH8_8 * *_S763)) >= -0.5f));
    float3  v_viewdir_32 = {};
    float _S765 = _S764.x;
    *&(v_coeff_dc_15->x) = *&(v_coeff_dc_15->x) + 0.282094806432724f * _S765;
    float _S766 = -0.48860251903533936f * y_17;
    float temp_59 = _S766 * _S765;
    bool _S767;
    if((F32_isfinite((temp_59))))
    {
        _S767 = temp_59 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S768 = atomicAdd(v_coeffs_15 + int(0), temp_59);
    }
    float _S769 = 0.48860251903533936f * z_16;
    float temp_60 = _S769 * _S765;
    if((F32_isfinite((temp_60))))
    {
        _S767 = temp_60 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S770 = atomicAdd(v_coeffs_15 + int(3), temp_60);
    }
    float _S771 = -0.48860251903533936f * x_19;
    float temp_61 = _S771 * _S765;
    if((F32_isfinite((temp_61))))
    {
        _S767 = temp_61 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S772 = atomicAdd(v_coeffs_15 + int(6), temp_61);
    }
    float _S773 = -0.48860251903533936f * *_S741 * _S765;
    float _S774 = -0.48860251903533936f * *_S739 * _S765;
    float _S775 = 0.48860251903533936f * *_S740 * _S765;
    float temp_62 = pSH4_8 * _S765;
    if((F32_isfinite((temp_62))))
    {
        _S767 = temp_62 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S776 = atomicAdd(v_coeffs_15 + int(9), temp_62);
    }
    float temp_63 = pSH5_8 * _S765;
    if((F32_isfinite((temp_63))))
    {
        _S767 = temp_63 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S777 = atomicAdd(v_coeffs_15 + int(12), temp_63);
    }
    float temp_64 = pSH6_9 * _S765;
    if((F32_isfinite((temp_64))))
    {
        _S767 = temp_64 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S778 = atomicAdd(v_coeffs_15 + int(15), temp_64);
    }
    float temp_65 = pSH7_8 * _S765;
    if((F32_isfinite((temp_65))))
    {
        _S767 = temp_65 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S779 = atomicAdd(v_coeffs_15 + int(18), temp_65);
    }
    float temp_66 = pSH8_8 * _S765;
    if((F32_isfinite((temp_66))))
    {
        _S767 = temp_66 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S780 = atomicAdd(v_coeffs_15 + int(21), temp_66);
    }
    float pSH6_z_3 = 1.89234936237335205f * z_16;
    float pSH7_z_1 = -1.09254848957061768f * x_19;
    float pSH5_z_1 = -1.09254848957061768f * y_17;
    float pSH8_x_7 = 0.54627424478530884f * _S742;
    float pSH8_y_1 = 0.54627424478530884f * (-2.0f * y_17);
    float pSH4_x_1 = 0.54627424478530884f * (2.0f * y_17);
    float3  dir_n_11 = make_float3 (x_19, y_17, z_16);
    float3  v_dir_n_17 = make_float3 (_S773 + _S765 * (pSH4_x_1 * *_S743 + pSH8_x_7 * *_S747 + fTmp0B_11 * *_S746), _S774 + _S765 * (pSH8_x_7 * *_S743 + pSH8_y_1 * *_S747 + fTmp0B_11 * *_S744), _S775 + _S765 * (pSH6_z_3 * *_S745 + pSH7_z_1 * *_S746 + pSH5_z_1 * *_S744));
    float3  v_viewdir_33 = v_viewdir_32 + (v_dir_n_17 - make_float3 (dot_0(v_dir_n_17, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    float _S781 = _S764.y;
    *&(v_coeff_dc_15->y) = *&(v_coeff_dc_15->y) + 0.282094806432724f * _S781;
    float temp_67 = _S766 * _S781;
    if((F32_isfinite((temp_67))))
    {
        _S767 = temp_67 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S782 = atomicAdd(v_coeffs_15 + int(1), temp_67);
    }
    float temp_68 = _S769 * _S781;
    if((F32_isfinite((temp_68))))
    {
        _S767 = temp_68 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S783 = atomicAdd(v_coeffs_15 + int(4), temp_68);
    }
    float temp_69 = _S771 * _S781;
    if((F32_isfinite((temp_69))))
    {
        _S767 = temp_69 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S784 = atomicAdd(v_coeffs_15 + int(7), temp_69);
    }
    float _S785 = -0.48860251903533936f * *_S750 * _S781;
    float _S786 = -0.48860251903533936f * *_S748 * _S781;
    float _S787 = 0.48860251903533936f * *_S749 * _S781;
    float temp_70 = pSH4_8 * _S781;
    if((F32_isfinite((temp_70))))
    {
        _S767 = temp_70 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S788 = atomicAdd(v_coeffs_15 + int(10), temp_70);
    }
    float temp_71 = pSH5_8 * _S781;
    if((F32_isfinite((temp_71))))
    {
        _S767 = temp_71 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S789 = atomicAdd(v_coeffs_15 + int(13), temp_71);
    }
    float temp_72 = pSH6_9 * _S781;
    if((F32_isfinite((temp_72))))
    {
        _S767 = temp_72 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S790 = atomicAdd(v_coeffs_15 + int(16), temp_72);
    }
    float temp_73 = pSH7_8 * _S781;
    if((F32_isfinite((temp_73))))
    {
        _S767 = temp_73 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S791 = atomicAdd(v_coeffs_15 + int(19), temp_73);
    }
    float temp_74 = pSH8_8 * _S781;
    if((F32_isfinite((temp_74))))
    {
        _S767 = temp_74 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S792 = atomicAdd(v_coeffs_15 + int(22), temp_74);
    }
    float3  v_dir_n_18 = make_float3 (_S785 + _S781 * (pSH4_x_1 * *_S751 + pSH8_x_7 * *_S755 + fTmp0B_11 * *_S754), _S786 + _S781 * (pSH8_x_7 * *_S751 + pSH8_y_1 * *_S755 + fTmp0B_11 * *_S752), _S787 + _S781 * (pSH6_z_3 * *_S753 + pSH7_z_1 * *_S754 + pSH5_z_1 * *_S752));
    float3  v_viewdir_34 = v_viewdir_33 + (v_dir_n_18 - make_float3 (dot_0(v_dir_n_18, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    float _S793 = _S764.z;
    *&(v_coeff_dc_15->z) = *&(v_coeff_dc_15->z) + 0.282094806432724f * _S793;
    float temp_75 = _S766 * _S793;
    if((F32_isfinite((temp_75))))
    {
        _S767 = temp_75 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S794 = atomicAdd(v_coeffs_15 + int(2), temp_75);
    }
    float temp_76 = _S769 * _S793;
    if((F32_isfinite((temp_76))))
    {
        _S767 = temp_76 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S795 = atomicAdd(v_coeffs_15 + int(5), temp_76);
    }
    float temp_77 = _S771 * _S793;
    if((F32_isfinite((temp_77))))
    {
        _S767 = temp_77 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S796 = atomicAdd(v_coeffs_15 + int(8), temp_77);
    }
    float _S797 = -0.48860251903533936f * *_S758 * _S793;
    float _S798 = -0.48860251903533936f * *_S756 * _S793;
    float _S799 = 0.48860251903533936f * *_S757 * _S793;
    float temp_78 = pSH4_8 * _S793;
    if((F32_isfinite((temp_78))))
    {
        _S767 = temp_78 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S800 = atomicAdd(v_coeffs_15 + int(11), temp_78);
    }
    float temp_79 = pSH5_8 * _S793;
    if((F32_isfinite((temp_79))))
    {
        _S767 = temp_79 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S801 = atomicAdd(v_coeffs_15 + int(14), temp_79);
    }
    float temp_80 = pSH6_9 * _S793;
    if((F32_isfinite((temp_80))))
    {
        _S767 = temp_80 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S802 = atomicAdd(v_coeffs_15 + int(17), temp_80);
    }
    float temp_81 = pSH7_8 * _S793;
    if((F32_isfinite((temp_81))))
    {
        _S767 = temp_81 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S803 = atomicAdd(v_coeffs_15 + int(20), temp_81);
    }
    float temp_82 = pSH8_8 * _S793;
    if((F32_isfinite((temp_82))))
    {
        _S767 = temp_82 != 0.0f;
    }
    else
    {
        _S767 = false;
    }
    if(_S767)
    {
        float _S804 = atomicAdd(v_coeffs_15 + int(23), temp_82);
    }
    float3  v_dir_n_19 = make_float3 (_S797 + _S793 * (pSH4_x_1 * *_S759 + pSH8_x_7 * *_S763 + fTmp0B_11 * *_S762), _S798 + _S793 * (pSH8_x_7 * *_S759 + pSH8_y_1 * *_S763 + fTmp0B_11 * *_S760), _S799 + _S793 * (pSH6_z_3 * *_S761 + pSH7_z_1 * *_S762 + pSH5_z_1 * *_S760));
    float3  v_viewdir_35 = v_viewdir_34 + (v_dir_n_19 - make_float3 (dot_0(v_dir_n_19, dir_n_11)) * dir_n_11) * make_float3 (inv_norm_17);
    Matrix<float, 3, 3>  _S805 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S806;
    (&_S806)->primal_0 = _S733;
    (&_S806)->differential_0 = _S805;
    float3  _S807 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S808;
    (&_S808)->primal_0 = t_23;
    (&_S808)->differential_0 = _S807;
    s_bwd_prop_mul_0(&_S806, &_S808, v_viewdir_35);
    Matrix<float, 3, 3>  _S809 = transpose_0(_S806.differential_0);
    *v_mean_15 = *v_mean_15 + v_viewdir_35;
    *v_R_15 = *v_R_15 + _S809;
    *v_t_15 = *v_t_15 + _S808.differential_0;
    return;
}

inline __device__ float3  sh3_to_color(float3  mean_24, Matrix<float, 3, 3>  R_24, float3  t_24, float3  coeff_dc_24, float * coeffs_24)
{
    float3  _S810 = mean_24 + mul_0(transpose_0(R_24), t_24);
    float _S811 = _S810.x;
    float _S812 = _S810.y;
    float _S813 = _S810.z;
    float inv_norm_18 = (F32_rsqrt((_S811 * _S811 + _S812 * _S812 + _S813 * _S813)));
    float x_20 = _S811 * inv_norm_18;
    float y_18 = _S812 * inv_norm_18;
    float z_17 = _S813 * inv_norm_18;
    float _S814 = - y_18;
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
    return max_0(make_float3 (0.282094806432724f * coeff_dc_24.x + 0.48860251903533936f * (_S814 * *(coeffs_24 + int(0)) + z_17 * *(coeffs_24 + int(3)) - x_20 * *(coeffs_24 + int(6))) + (pSH4_9 * *(coeffs_24 + int(9)) + pSH5_9 * *(coeffs_24 + int(12)) + pSH6_10 * *(coeffs_24 + int(15)) + pSH7_9 * *(coeffs_24 + int(18)) + pSH8_9 * *(coeffs_24 + int(21))) + (pSH9_4 * *(coeffs_24 + int(24)) + pSH10_4 * *(coeffs_24 + int(27)) + pSH11_4 * *(coeffs_24 + int(30)) + pSH12_5 * *(coeffs_24 + int(33)) + pSH13_4 * *(coeffs_24 + int(36)) + pSH14_4 * *(coeffs_24 + int(39)) + pSH15_4 * *(coeffs_24 + int(42))), 0.282094806432724f * coeff_dc_24.y + 0.48860251903533936f * (_S814 * *(coeffs_24 + int(1)) + z_17 * *(coeffs_24 + int(4)) - x_20 * *(coeffs_24 + int(7))) + (pSH4_9 * *(coeffs_24 + int(10)) + pSH5_9 * *(coeffs_24 + int(13)) + pSH6_10 * *(coeffs_24 + int(16)) + pSH7_9 * *(coeffs_24 + int(19)) + pSH8_9 * *(coeffs_24 + int(22))) + (pSH9_4 * *(coeffs_24 + int(25)) + pSH10_4 * *(coeffs_24 + int(28)) + pSH11_4 * *(coeffs_24 + int(31)) + pSH12_5 * *(coeffs_24 + int(34)) + pSH13_4 * *(coeffs_24 + int(37)) + pSH14_4 * *(coeffs_24 + int(40)) + pSH15_4 * *(coeffs_24 + int(43))), 0.282094806432724f * coeff_dc_24.z + 0.48860251903533936f * (_S814 * *(coeffs_24 + int(2)) + z_17 * *(coeffs_24 + int(5)) - x_20 * *(coeffs_24 + int(8))) + (pSH4_9 * *(coeffs_24 + int(11)) + pSH5_9 * *(coeffs_24 + int(14)) + pSH6_10 * *(coeffs_24 + int(17)) + pSH7_9 * *(coeffs_24 + int(20)) + pSH8_9 * *(coeffs_24 + int(23))) + (pSH9_4 * *(coeffs_24 + int(26)) + pSH10_4 * *(coeffs_24 + int(29)) + pSH11_4 * *(coeffs_24 + int(32)) + pSH12_5 * *(coeffs_24 + int(35)) + pSH13_4 * *(coeffs_24 + int(38)) + pSH14_4 * *(coeffs_24 + int(41)) + pSH15_4 * *(coeffs_24 + int(44)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh3_to_color_vjp_inplace(float3  mean_25, Matrix<float, 3, 3>  R_25, float3  t_25, float3  coeff_dc_25, float * coeffs_25, float3  v_colors_16, float3  * v_coeff_dc_16, float * v_coeffs_16, float3  * v_mean_16, Matrix<float, 3, 3>  * v_R_16, float3  * v_t_16)
{
    Matrix<float, 3, 3>  _S815 = transpose_0(R_25);
    float3  _S816 = mean_25 + mul_0(_S815, t_25);
    float _S817 = _S816.x;
    float _S818 = _S816.y;
    float _S819 = _S816.z;
    float inv_norm_19 = (F32_rsqrt((_S817 * _S817 + _S818 * _S818 + _S819 * _S819)));
    float x_21 = _S817 * inv_norm_19;
    float y_19 = _S818 * inv_norm_19;
    float z_18 = _S819 * inv_norm_19;
    float _S820 = - y_19;
    float * _S821 = coeffs_25 + int(0);
    float * _S822 = coeffs_25 + int(3);
    float * _S823 = coeffs_25 + int(6);
    float z2_7 = z_18 * z_18;
    float fTmp0B_13 = -1.09254848957061768f * z_18;
    float fC1_7 = x_21 * x_21 - y_19 * y_19;
    float _S824 = 2.0f * x_21;
    float fS1_7 = _S824 * y_19;
    float pSH6_11 = 0.94617468118667603f * z2_7 - 0.31539157032966614f;
    float pSH7_10 = fTmp0B_13 * x_21;
    float pSH5_10 = fTmp0B_13 * y_19;
    float pSH8_10 = 0.54627424478530884f * fC1_7;
    float pSH4_10 = 0.54627424478530884f * fS1_7;
    float * _S825 = coeffs_25 + int(9);
    float * _S826 = coeffs_25 + int(12);
    float * _S827 = coeffs_25 + int(15);
    float * _S828 = coeffs_25 + int(18);
    float * _S829 = coeffs_25 + int(21);
    float fTmp0C_7 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_18;
    float pSH12_6 = z_18 * (1.86588168144226074f * z2_7 - 1.11952900886535645f);
    float pSH13_5 = fTmp0C_7 * x_21;
    float pSH11_5 = fTmp0C_7 * y_19;
    float pSH14_5 = fTmp1B_7 * fC1_7;
    float pSH10_5 = fTmp1B_7 * fS1_7;
    float pSH15_5 = -0.59004360437393188f * (x_21 * fC1_7 - y_19 * fS1_7);
    float pSH9_5 = -0.59004360437393188f * (x_21 * fS1_7 + y_19 * fC1_7);
    float * _S830 = coeffs_25 + int(24);
    float * _S831 = coeffs_25 + int(27);
    float * _S832 = coeffs_25 + int(30);
    float * _S833 = coeffs_25 + int(33);
    float * _S834 = coeffs_25 + int(36);
    float * _S835 = coeffs_25 + int(39);
    float * _S836 = coeffs_25 + int(42);
    float * _S837 = coeffs_25 + int(1);
    float * _S838 = coeffs_25 + int(4);
    float * _S839 = coeffs_25 + int(7);
    float * _S840 = coeffs_25 + int(10);
    float * _S841 = coeffs_25 + int(13);
    float * _S842 = coeffs_25 + int(16);
    float * _S843 = coeffs_25 + int(19);
    float * _S844 = coeffs_25 + int(22);
    float * _S845 = coeffs_25 + int(25);
    float * _S846 = coeffs_25 + int(28);
    float * _S847 = coeffs_25 + int(31);
    float * _S848 = coeffs_25 + int(34);
    float * _S849 = coeffs_25 + int(37);
    float * _S850 = coeffs_25 + int(40);
    float * _S851 = coeffs_25 + int(43);
    float * _S852 = coeffs_25 + int(2);
    float * _S853 = coeffs_25 + int(5);
    float * _S854 = coeffs_25 + int(8);
    float * _S855 = coeffs_25 + int(11);
    float * _S856 = coeffs_25 + int(14);
    float * _S857 = coeffs_25 + int(17);
    float * _S858 = coeffs_25 + int(20);
    float * _S859 = coeffs_25 + int(23);
    float * _S860 = coeffs_25 + int(26);
    float * _S861 = coeffs_25 + int(29);
    float * _S862 = coeffs_25 + int(32);
    float * _S863 = coeffs_25 + int(35);
    float * _S864 = coeffs_25 + int(38);
    float * _S865 = coeffs_25 + int(41);
    float * _S866 = coeffs_25 + int(44);
    float3  _S867 = v_colors_16 * make_float3 (float((0.282094806432724f * coeff_dc_25.x + 0.48860251903533936f * (_S820 * *_S821 + z_18 * *_S822 - x_21 * *_S823) + (pSH4_10 * *_S825 + pSH5_10 * *_S826 + pSH6_11 * *_S827 + pSH7_10 * *_S828 + pSH8_10 * *_S829) + (pSH9_5 * *_S830 + pSH10_5 * *_S831 + pSH11_5 * *_S832 + pSH12_6 * *_S833 + pSH13_5 * *_S834 + pSH14_5 * *_S835 + pSH15_5 * *_S836)) >= -0.5f), float((0.282094806432724f * coeff_dc_25.y + 0.48860251903533936f * (_S820 * *_S837 + z_18 * *_S838 - x_21 * *_S839) + (pSH4_10 * *_S840 + pSH5_10 * *_S841 + pSH6_11 * *_S842 + pSH7_10 * *_S843 + pSH8_10 * *_S844) + (pSH9_5 * *_S845 + pSH10_5 * *_S846 + pSH11_5 * *_S847 + pSH12_6 * *_S848 + pSH13_5 * *_S849 + pSH14_5 * *_S850 + pSH15_5 * *_S851)) >= -0.5f), float((0.282094806432724f * coeff_dc_25.z + 0.48860251903533936f * (_S820 * *_S852 + z_18 * *_S853 - x_21 * *_S854) + (pSH4_10 * *_S855 + pSH5_10 * *_S856 + pSH6_11 * *_S857 + pSH7_10 * *_S858 + pSH8_10 * *_S859) + (pSH9_5 * *_S860 + pSH10_5 * *_S861 + pSH11_5 * *_S862 + pSH12_6 * *_S863 + pSH13_5 * *_S864 + pSH14_5 * *_S865 + pSH15_5 * *_S866)) >= -0.5f));
    float3  v_viewdir_36 = {};
    float _S868 = _S867.x;
    *&(v_coeff_dc_16->x) = *&(v_coeff_dc_16->x) + 0.282094806432724f * _S868;
    float * _S869 = v_coeffs_16 + int(0);
    float _S870 = -0.48860251903533936f * y_19;
    *_S869 = *_S869 + _S870 * _S868;
    float * _S871 = v_coeffs_16 + int(3);
    float _S872 = 0.48860251903533936f * z_18;
    *_S871 = *_S871 + _S872 * _S868;
    float * _S873 = v_coeffs_16 + int(6);
    float _S874 = -0.48860251903533936f * x_21;
    *_S873 = *_S873 + _S874 * _S868;
    float _S875 = -0.48860251903533936f * *_S823 * _S868;
    float _S876 = -0.48860251903533936f * *_S821 * _S868;
    float _S877 = 0.48860251903533936f * *_S822 * _S868;
    float * _S878 = v_coeffs_16 + int(9);
    *_S878 = *_S878 + pSH4_10 * _S868;
    float * _S879 = v_coeffs_16 + int(12);
    *_S879 = *_S879 + pSH5_10 * _S868;
    float * _S880 = v_coeffs_16 + int(15);
    *_S880 = *_S880 + pSH6_11 * _S868;
    float * _S881 = v_coeffs_16 + int(18);
    *_S881 = *_S881 + pSH7_10 * _S868;
    float * _S882 = v_coeffs_16 + int(21);
    *_S882 = *_S882 + pSH8_10 * _S868;
    float fC1_y_4 = -2.0f * y_19;
    float fS1_x_4 = 2.0f * y_19;
    float pSH6_z_4 = 1.89234936237335205f * z_18;
    float pSH7_z_2 = -1.09254848957061768f * x_21;
    float pSH5_z_2 = -1.09254848957061768f * y_19;
    float pSH8_x_8 = 0.54627424478530884f * _S824;
    float pSH8_y_2 = 0.54627424478530884f * fC1_y_4;
    float pSH4_x_2 = 0.54627424478530884f * fS1_x_4;
    float v_x_6 = _S875 + _S868 * (pSH4_x_2 * *_S825 + pSH8_x_8 * *_S829 + fTmp0B_13 * *_S828);
    float v_y_6 = _S876 + _S868 * (pSH8_x_8 * *_S825 + pSH8_y_2 * *_S829 + fTmp0B_13 * *_S826);
    float v_z_6 = _S877 + _S868 * (pSH6_z_4 * *_S827 + pSH7_z_2 * *_S828 + pSH5_z_2 * *_S826);
    float * _S883 = v_coeffs_16 + int(24);
    *_S883 = *_S883 + pSH9_5 * _S868;
    float * _S884 = v_coeffs_16 + int(27);
    *_S884 = *_S884 + pSH10_5 * _S868;
    float * _S885 = v_coeffs_16 + int(30);
    *_S885 = *_S885 + pSH11_5 * _S868;
    float * _S886 = v_coeffs_16 + int(33);
    *_S886 = *_S886 + pSH12_6 * _S868;
    float * _S887 = v_coeffs_16 + int(36);
    *_S887 = *_S887 + pSH13_5 * _S868;
    float * _S888 = v_coeffs_16 + int(39);
    *_S888 = *_S888 + pSH14_5 * _S868;
    float * _S889 = v_coeffs_16 + int(42);
    *_S889 = *_S889 + pSH15_5 * _S868;
    float fTmp0C_z_4 = -4.57045793533325195f * z_18;
    float _S890 = x_21 * _S824;
    float _S891 = y_19 * _S824;
    float pSH12_z_2 = 5.59764480590820312f * z2_7 - 1.11952900886535645f;
    float pSH13_z_0 = fTmp0C_z_4 * x_21;
    float pSH11_z_0 = fTmp0C_z_4 * y_19;
    float pSH14_x_4 = fTmp1B_7 * _S824;
    float pSH14_y_0 = fTmp1B_7 * fC1_y_4;
    float pSH14_z_0 = 1.44530570507049561f * fC1_7;
    float pSH10_x_0 = fTmp1B_7 * fS1_x_4;
    float pSH10_z_0 = 1.44530570507049561f * fS1_7;
    float pSH15_x_0 = -0.59004360437393188f * (fC1_7 + _S890 - y_19 * fS1_x_4);
    float pSH15_y_0 = -0.59004360437393188f * (x_21 * fC1_y_4 - fS1_7 - _S891);
    float pSH9_x_0 = -0.59004360437393188f * (fS1_7 + x_21 * fS1_x_4 + _S891);
    float pSH9_y_0 = -0.59004360437393188f * (_S890 + fC1_7 + y_19 * fC1_y_4);
    float3  dir_n_12 = make_float3 (x_21, y_19, z_18);
    float3  v_dir_n_20 = make_float3 (v_x_6 + _S868 * (pSH9_x_0 * *_S830 + pSH15_x_0 * *_S836 + pSH10_x_0 * *_S831 + pSH14_x_4 * *_S835 + fTmp0C_7 * *_S834), v_y_6 + _S868 * (pSH9_y_0 * *_S830 + pSH15_y_0 * *_S836 + pSH14_x_4 * *_S831 + pSH14_y_0 * *_S835 + fTmp0C_7 * *_S832), v_z_6 + _S868 * (pSH12_z_2 * *_S833 + pSH13_z_0 * *_S834 + pSH11_z_0 * *_S832 + pSH14_z_0 * *_S835 + pSH10_z_0 * *_S831));
    float3  v_viewdir_37 = v_viewdir_36 + (v_dir_n_20 - make_float3 (dot_0(v_dir_n_20, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    float _S892 = _S867.y;
    *&(v_coeff_dc_16->y) = *&(v_coeff_dc_16->y) + 0.282094806432724f * _S892;
    float * _S893 = v_coeffs_16 + int(1);
    *_S893 = *_S893 + _S870 * _S892;
    float * _S894 = v_coeffs_16 + int(4);
    *_S894 = *_S894 + _S872 * _S892;
    float * _S895 = v_coeffs_16 + int(7);
    *_S895 = *_S895 + _S874 * _S892;
    float _S896 = -0.48860251903533936f * *_S839 * _S892;
    float _S897 = -0.48860251903533936f * *_S837 * _S892;
    float _S898 = 0.48860251903533936f * *_S838 * _S892;
    float * _S899 = v_coeffs_16 + int(10);
    *_S899 = *_S899 + pSH4_10 * _S892;
    float * _S900 = v_coeffs_16 + int(13);
    *_S900 = *_S900 + pSH5_10 * _S892;
    float * _S901 = v_coeffs_16 + int(16);
    *_S901 = *_S901 + pSH6_11 * _S892;
    float * _S902 = v_coeffs_16 + int(19);
    *_S902 = *_S902 + pSH7_10 * _S892;
    float * _S903 = v_coeffs_16 + int(22);
    *_S903 = *_S903 + pSH8_10 * _S892;
    float v_x_7 = _S896 + _S892 * (pSH4_x_2 * *_S840 + pSH8_x_8 * *_S844 + fTmp0B_13 * *_S843);
    float v_y_7 = _S897 + _S892 * (pSH8_x_8 * *_S840 + pSH8_y_2 * *_S844 + fTmp0B_13 * *_S841);
    float v_z_7 = _S898 + _S892 * (pSH6_z_4 * *_S842 + pSH7_z_2 * *_S843 + pSH5_z_2 * *_S841);
    float * _S904 = v_coeffs_16 + int(25);
    *_S904 = *_S904 + pSH9_5 * _S892;
    float * _S905 = v_coeffs_16 + int(28);
    *_S905 = *_S905 + pSH10_5 * _S892;
    float * _S906 = v_coeffs_16 + int(31);
    *_S906 = *_S906 + pSH11_5 * _S892;
    float * _S907 = v_coeffs_16 + int(34);
    *_S907 = *_S907 + pSH12_6 * _S892;
    float * _S908 = v_coeffs_16 + int(37);
    *_S908 = *_S908 + pSH13_5 * _S892;
    float * _S909 = v_coeffs_16 + int(40);
    *_S909 = *_S909 + pSH14_5 * _S892;
    float * _S910 = v_coeffs_16 + int(43);
    *_S910 = *_S910 + pSH15_5 * _S892;
    float3  v_dir_n_21 = make_float3 (v_x_7 + _S892 * (pSH9_x_0 * *_S845 + pSH15_x_0 * *_S851 + pSH10_x_0 * *_S846 + pSH14_x_4 * *_S850 + fTmp0C_7 * *_S849), v_y_7 + _S892 * (pSH9_y_0 * *_S845 + pSH15_y_0 * *_S851 + pSH14_x_4 * *_S846 + pSH14_y_0 * *_S850 + fTmp0C_7 * *_S847), v_z_7 + _S892 * (pSH12_z_2 * *_S848 + pSH13_z_0 * *_S849 + pSH11_z_0 * *_S847 + pSH14_z_0 * *_S850 + pSH10_z_0 * *_S846));
    float3  v_viewdir_38 = v_viewdir_37 + (v_dir_n_21 - make_float3 (dot_0(v_dir_n_21, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    float _S911 = _S867.z;
    *&(v_coeff_dc_16->z) = *&(v_coeff_dc_16->z) + 0.282094806432724f * _S911;
    float * _S912 = v_coeffs_16 + int(2);
    *_S912 = *_S912 + _S870 * _S911;
    float * _S913 = v_coeffs_16 + int(5);
    *_S913 = *_S913 + _S872 * _S911;
    float * _S914 = v_coeffs_16 + int(8);
    *_S914 = *_S914 + _S874 * _S911;
    float _S915 = -0.48860251903533936f * *_S854 * _S911;
    float _S916 = -0.48860251903533936f * *_S852 * _S911;
    float _S917 = 0.48860251903533936f * *_S853 * _S911;
    float * _S918 = v_coeffs_16 + int(11);
    *_S918 = *_S918 + pSH4_10 * _S911;
    float * _S919 = v_coeffs_16 + int(14);
    *_S919 = *_S919 + pSH5_10 * _S911;
    float * _S920 = v_coeffs_16 + int(17);
    *_S920 = *_S920 + pSH6_11 * _S911;
    float * _S921 = v_coeffs_16 + int(20);
    *_S921 = *_S921 + pSH7_10 * _S911;
    float * _S922 = v_coeffs_16 + int(23);
    *_S922 = *_S922 + pSH8_10 * _S911;
    float v_x_8 = _S915 + _S911 * (pSH4_x_2 * *_S855 + pSH8_x_8 * *_S859 + fTmp0B_13 * *_S858);
    float v_y_8 = _S916 + _S911 * (pSH8_x_8 * *_S855 + pSH8_y_2 * *_S859 + fTmp0B_13 * *_S856);
    float v_z_8 = _S917 + _S911 * (pSH6_z_4 * *_S857 + pSH7_z_2 * *_S858 + pSH5_z_2 * *_S856);
    float * _S923 = v_coeffs_16 + int(26);
    *_S923 = *_S923 + pSH9_5 * _S911;
    float * _S924 = v_coeffs_16 + int(29);
    *_S924 = *_S924 + pSH10_5 * _S911;
    float * _S925 = v_coeffs_16 + int(32);
    *_S925 = *_S925 + pSH11_5 * _S911;
    float * _S926 = v_coeffs_16 + int(35);
    *_S926 = *_S926 + pSH12_6 * _S911;
    float * _S927 = v_coeffs_16 + int(38);
    *_S927 = *_S927 + pSH13_5 * _S911;
    float * _S928 = v_coeffs_16 + int(41);
    *_S928 = *_S928 + pSH14_5 * _S911;
    float * _S929 = v_coeffs_16 + int(44);
    *_S929 = *_S929 + pSH15_5 * _S911;
    float3  v_dir_n_22 = make_float3 (v_x_8 + _S911 * (pSH9_x_0 * *_S860 + pSH15_x_0 * *_S866 + pSH10_x_0 * *_S861 + pSH14_x_4 * *_S865 + fTmp0C_7 * *_S864), v_y_8 + _S911 * (pSH9_y_0 * *_S860 + pSH15_y_0 * *_S866 + pSH14_x_4 * *_S861 + pSH14_y_0 * *_S865 + fTmp0C_7 * *_S862), v_z_8 + _S911 * (pSH12_z_2 * *_S863 + pSH13_z_0 * *_S864 + pSH11_z_0 * *_S862 + pSH14_z_0 * *_S865 + pSH10_z_0 * *_S861));
    float3  v_viewdir_39 = v_viewdir_38 + (v_dir_n_22 - make_float3 (dot_0(v_dir_n_22, dir_n_12)) * dir_n_12) * make_float3 (inv_norm_19);
    Matrix<float, 3, 3>  _S930 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S931;
    (&_S931)->primal_0 = _S815;
    (&_S931)->differential_0 = _S930;
    float3  _S932 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S933;
    (&_S933)->primal_0 = t_25;
    (&_S933)->differential_0 = _S932;
    s_bwd_prop_mul_0(&_S931, &_S933, v_viewdir_39);
    Matrix<float, 3, 3>  _S934 = transpose_0(_S931.differential_0);
    *v_mean_16 = *v_mean_16 + v_viewdir_39;
    *v_R_16 = *v_R_16 + _S934;
    *v_t_16 = *v_t_16 + _S933.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_atomic(float3  mean_26, Matrix<float, 3, 3>  R_26, float3  t_26, float3  coeff_dc_26, float * coeffs_26, float3  v_colors_17, float3  * v_coeff_dc_17, float * v_coeffs_17, float3  * v_mean_17, Matrix<float, 3, 3>  * v_R_17, float3  * v_t_17)
{
    Matrix<float, 3, 3>  _S935 = transpose_0(R_26);
    float3  _S936 = mean_26 + mul_0(_S935, t_26);
    float _S937 = _S936.x;
    float _S938 = _S936.y;
    float _S939 = _S936.z;
    float inv_norm_20 = (F32_rsqrt((_S937 * _S937 + _S938 * _S938 + _S939 * _S939)));
    float x_22 = _S937 * inv_norm_20;
    float y_20 = _S938 * inv_norm_20;
    float z_19 = _S939 * inv_norm_20;
    float _S940 = - y_20;
    float * _S941 = coeffs_26 + int(0);
    float * _S942 = coeffs_26 + int(3);
    float * _S943 = coeffs_26 + int(6);
    float z2_8 = z_19 * z_19;
    float fTmp0B_14 = -1.09254848957061768f * z_19;
    float fC1_8 = x_22 * x_22 - y_20 * y_20;
    float _S944 = 2.0f * x_22;
    float fS1_8 = _S944 * y_20;
    float pSH6_12 = 0.94617468118667603f * z2_8 - 0.31539157032966614f;
    float pSH7_11 = fTmp0B_14 * x_22;
    float pSH5_11 = fTmp0B_14 * y_20;
    float pSH8_11 = 0.54627424478530884f * fC1_8;
    float pSH4_11 = 0.54627424478530884f * fS1_8;
    float * _S945 = coeffs_26 + int(9);
    float * _S946 = coeffs_26 + int(12);
    float * _S947 = coeffs_26 + int(15);
    float * _S948 = coeffs_26 + int(18);
    float * _S949 = coeffs_26 + int(21);
    float fTmp0C_8 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_19;
    float pSH12_7 = z_19 * (1.86588168144226074f * z2_8 - 1.11952900886535645f);
    float pSH13_6 = fTmp0C_8 * x_22;
    float pSH11_6 = fTmp0C_8 * y_20;
    float pSH14_6 = fTmp1B_8 * fC1_8;
    float pSH10_6 = fTmp1B_8 * fS1_8;
    float pSH15_6 = -0.59004360437393188f * (x_22 * fC1_8 - y_20 * fS1_8);
    float pSH9_6 = -0.59004360437393188f * (x_22 * fS1_8 + y_20 * fC1_8);
    float * _S950 = coeffs_26 + int(24);
    float * _S951 = coeffs_26 + int(27);
    float * _S952 = coeffs_26 + int(30);
    float * _S953 = coeffs_26 + int(33);
    float * _S954 = coeffs_26 + int(36);
    float * _S955 = coeffs_26 + int(39);
    float * _S956 = coeffs_26 + int(42);
    float * _S957 = coeffs_26 + int(1);
    float * _S958 = coeffs_26 + int(4);
    float * _S959 = coeffs_26 + int(7);
    float * _S960 = coeffs_26 + int(10);
    float * _S961 = coeffs_26 + int(13);
    float * _S962 = coeffs_26 + int(16);
    float * _S963 = coeffs_26 + int(19);
    float * _S964 = coeffs_26 + int(22);
    float * _S965 = coeffs_26 + int(25);
    float * _S966 = coeffs_26 + int(28);
    float * _S967 = coeffs_26 + int(31);
    float * _S968 = coeffs_26 + int(34);
    float * _S969 = coeffs_26 + int(37);
    float * _S970 = coeffs_26 + int(40);
    float * _S971 = coeffs_26 + int(43);
    float * _S972 = coeffs_26 + int(2);
    float * _S973 = coeffs_26 + int(5);
    float * _S974 = coeffs_26 + int(8);
    float * _S975 = coeffs_26 + int(11);
    float * _S976 = coeffs_26 + int(14);
    float * _S977 = coeffs_26 + int(17);
    float * _S978 = coeffs_26 + int(20);
    float * _S979 = coeffs_26 + int(23);
    float * _S980 = coeffs_26 + int(26);
    float * _S981 = coeffs_26 + int(29);
    float * _S982 = coeffs_26 + int(32);
    float * _S983 = coeffs_26 + int(35);
    float * _S984 = coeffs_26 + int(38);
    float * _S985 = coeffs_26 + int(41);
    float * _S986 = coeffs_26 + int(44);
    float3  _S987 = v_colors_17 * make_float3 (float((0.282094806432724f * coeff_dc_26.x + 0.48860251903533936f * (_S940 * *_S941 + z_19 * *_S942 - x_22 * *_S943) + (pSH4_11 * *_S945 + pSH5_11 * *_S946 + pSH6_12 * *_S947 + pSH7_11 * *_S948 + pSH8_11 * *_S949) + (pSH9_6 * *_S950 + pSH10_6 * *_S951 + pSH11_6 * *_S952 + pSH12_7 * *_S953 + pSH13_6 * *_S954 + pSH14_6 * *_S955 + pSH15_6 * *_S956)) >= -0.5f), float((0.282094806432724f * coeff_dc_26.y + 0.48860251903533936f * (_S940 * *_S957 + z_19 * *_S958 - x_22 * *_S959) + (pSH4_11 * *_S960 + pSH5_11 * *_S961 + pSH6_12 * *_S962 + pSH7_11 * *_S963 + pSH8_11 * *_S964) + (pSH9_6 * *_S965 + pSH10_6 * *_S966 + pSH11_6 * *_S967 + pSH12_7 * *_S968 + pSH13_6 * *_S969 + pSH14_6 * *_S970 + pSH15_6 * *_S971)) >= -0.5f), float((0.282094806432724f * coeff_dc_26.z + 0.48860251903533936f * (_S940 * *_S972 + z_19 * *_S973 - x_22 * *_S974) + (pSH4_11 * *_S975 + pSH5_11 * *_S976 + pSH6_12 * *_S977 + pSH7_11 * *_S978 + pSH8_11 * *_S979) + (pSH9_6 * *_S980 + pSH10_6 * *_S981 + pSH11_6 * *_S982 + pSH12_7 * *_S983 + pSH13_6 * *_S984 + pSH14_6 * *_S985 + pSH15_6 * *_S986)) >= -0.5f));
    float3  v_viewdir_40 = {};
    float _S988 = _S987.x;
    *&(v_coeff_dc_17->x) = *&(v_coeff_dc_17->x) + 0.282094806432724f * _S988;
    float _S989 = -0.48860251903533936f * y_20;
    float temp_83 = _S989 * _S988;
    bool _S990;
    if((F32_isfinite((temp_83))))
    {
        _S990 = temp_83 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S991 = atomicAdd(v_coeffs_17 + int(0), temp_83);
    }
    float _S992 = 0.48860251903533936f * z_19;
    float temp_84 = _S992 * _S988;
    if((F32_isfinite((temp_84))))
    {
        _S990 = temp_84 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S993 = atomicAdd(v_coeffs_17 + int(3), temp_84);
    }
    float _S994 = -0.48860251903533936f * x_22;
    float temp_85 = _S994 * _S988;
    if((F32_isfinite((temp_85))))
    {
        _S990 = temp_85 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S995 = atomicAdd(v_coeffs_17 + int(6), temp_85);
    }
    float _S996 = -0.48860251903533936f * *_S943 * _S988;
    float _S997 = -0.48860251903533936f * *_S941 * _S988;
    float _S998 = 0.48860251903533936f * *_S942 * _S988;
    float temp_86 = pSH4_11 * _S988;
    if((F32_isfinite((temp_86))))
    {
        _S990 = temp_86 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S999 = atomicAdd(v_coeffs_17 + int(9), temp_86);
    }
    float temp_87 = pSH5_11 * _S988;
    if((F32_isfinite((temp_87))))
    {
        _S990 = temp_87 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1000 = atomicAdd(v_coeffs_17 + int(12), temp_87);
    }
    float temp_88 = pSH6_12 * _S988;
    if((F32_isfinite((temp_88))))
    {
        _S990 = temp_88 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1001 = atomicAdd(v_coeffs_17 + int(15), temp_88);
    }
    float temp_89 = pSH7_11 * _S988;
    if((F32_isfinite((temp_89))))
    {
        _S990 = temp_89 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1002 = atomicAdd(v_coeffs_17 + int(18), temp_89);
    }
    float temp_90 = pSH8_11 * _S988;
    if((F32_isfinite((temp_90))))
    {
        _S990 = temp_90 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1003 = atomicAdd(v_coeffs_17 + int(21), temp_90);
    }
    float fC1_y_5 = -2.0f * y_20;
    float fS1_x_5 = 2.0f * y_20;
    float pSH6_z_5 = 1.89234936237335205f * z_19;
    float pSH7_z_3 = -1.09254848957061768f * x_22;
    float pSH5_z_3 = -1.09254848957061768f * y_20;
    float pSH8_x_9 = 0.54627424478530884f * _S944;
    float pSH8_y_3 = 0.54627424478530884f * fC1_y_5;
    float pSH4_x_3 = 0.54627424478530884f * fS1_x_5;
    float v_x_9 = _S996 + _S988 * (pSH4_x_3 * *_S945 + pSH8_x_9 * *_S949 + fTmp0B_14 * *_S948);
    float v_y_9 = _S997 + _S988 * (pSH8_x_9 * *_S945 + pSH8_y_3 * *_S949 + fTmp0B_14 * *_S946);
    float v_z_9 = _S998 + _S988 * (pSH6_z_5 * *_S947 + pSH7_z_3 * *_S948 + pSH5_z_3 * *_S946);
    float temp_91 = pSH9_6 * _S988;
    if((F32_isfinite((temp_91))))
    {
        _S990 = temp_91 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1004 = atomicAdd(v_coeffs_17 + int(24), temp_91);
    }
    float temp_92 = pSH10_6 * _S988;
    if((F32_isfinite((temp_92))))
    {
        _S990 = temp_92 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1005 = atomicAdd(v_coeffs_17 + int(27), temp_92);
    }
    float temp_93 = pSH11_6 * _S988;
    if((F32_isfinite((temp_93))))
    {
        _S990 = temp_93 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1006 = atomicAdd(v_coeffs_17 + int(30), temp_93);
    }
    float temp_94 = pSH12_7 * _S988;
    if((F32_isfinite((temp_94))))
    {
        _S990 = temp_94 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1007 = atomicAdd(v_coeffs_17 + int(33), temp_94);
    }
    float temp_95 = pSH13_6 * _S988;
    if((F32_isfinite((temp_95))))
    {
        _S990 = temp_95 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1008 = atomicAdd(v_coeffs_17 + int(36), temp_95);
    }
    float temp_96 = pSH14_6 * _S988;
    if((F32_isfinite((temp_96))))
    {
        _S990 = temp_96 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1009 = atomicAdd(v_coeffs_17 + int(39), temp_96);
    }
    float temp_97 = pSH15_6 * _S988;
    if((F32_isfinite((temp_97))))
    {
        _S990 = temp_97 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1010 = atomicAdd(v_coeffs_17 + int(42), temp_97);
    }
    float fTmp0C_z_5 = -4.57045793533325195f * z_19;
    float _S1011 = x_22 * _S944;
    float _S1012 = y_20 * _S944;
    float pSH12_z_3 = 5.59764480590820312f * z2_8 - 1.11952900886535645f;
    float pSH13_z_1 = fTmp0C_z_5 * x_22;
    float pSH11_z_1 = fTmp0C_z_5 * y_20;
    float pSH14_x_5 = fTmp1B_8 * _S944;
    float pSH14_y_1 = fTmp1B_8 * fC1_y_5;
    float pSH14_z_1 = 1.44530570507049561f * fC1_8;
    float pSH10_x_1 = fTmp1B_8 * fS1_x_5;
    float pSH10_z_1 = 1.44530570507049561f * fS1_8;
    float pSH15_x_1 = -0.59004360437393188f * (fC1_8 + _S1011 - y_20 * fS1_x_5);
    float pSH15_y_1 = -0.59004360437393188f * (x_22 * fC1_y_5 - fS1_8 - _S1012);
    float pSH9_x_1 = -0.59004360437393188f * (fS1_8 + x_22 * fS1_x_5 + _S1012);
    float pSH9_y_1 = -0.59004360437393188f * (_S1011 + fC1_8 + y_20 * fC1_y_5);
    float3  dir_n_13 = make_float3 (x_22, y_20, z_19);
    float3  v_dir_n_23 = make_float3 (v_x_9 + _S988 * (pSH9_x_1 * *_S950 + pSH15_x_1 * *_S956 + pSH10_x_1 * *_S951 + pSH14_x_5 * *_S955 + fTmp0C_8 * *_S954), v_y_9 + _S988 * (pSH9_y_1 * *_S950 + pSH15_y_1 * *_S956 + pSH14_x_5 * *_S951 + pSH14_y_1 * *_S955 + fTmp0C_8 * *_S952), v_z_9 + _S988 * (pSH12_z_3 * *_S953 + pSH13_z_1 * *_S954 + pSH11_z_1 * *_S952 + pSH14_z_1 * *_S955 + pSH10_z_1 * *_S951));
    float3  v_viewdir_41 = v_viewdir_40 + (v_dir_n_23 - make_float3 (dot_0(v_dir_n_23, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    float _S1013 = _S987.y;
    *&(v_coeff_dc_17->y) = *&(v_coeff_dc_17->y) + 0.282094806432724f * _S1013;
    float temp_98 = _S989 * _S1013;
    if((F32_isfinite((temp_98))))
    {
        _S990 = temp_98 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1014 = atomicAdd(v_coeffs_17 + int(1), temp_98);
    }
    float temp_99 = _S992 * _S1013;
    if((F32_isfinite((temp_99))))
    {
        _S990 = temp_99 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1015 = atomicAdd(v_coeffs_17 + int(4), temp_99);
    }
    float temp_100 = _S994 * _S1013;
    if((F32_isfinite((temp_100))))
    {
        _S990 = temp_100 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1016 = atomicAdd(v_coeffs_17 + int(7), temp_100);
    }
    float _S1017 = -0.48860251903533936f * *_S959 * _S1013;
    float _S1018 = -0.48860251903533936f * *_S957 * _S1013;
    float _S1019 = 0.48860251903533936f * *_S958 * _S1013;
    float temp_101 = pSH4_11 * _S1013;
    if((F32_isfinite((temp_101))))
    {
        _S990 = temp_101 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1020 = atomicAdd(v_coeffs_17 + int(10), temp_101);
    }
    float temp_102 = pSH5_11 * _S1013;
    if((F32_isfinite((temp_102))))
    {
        _S990 = temp_102 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1021 = atomicAdd(v_coeffs_17 + int(13), temp_102);
    }
    float temp_103 = pSH6_12 * _S1013;
    if((F32_isfinite((temp_103))))
    {
        _S990 = temp_103 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1022 = atomicAdd(v_coeffs_17 + int(16), temp_103);
    }
    float temp_104 = pSH7_11 * _S1013;
    if((F32_isfinite((temp_104))))
    {
        _S990 = temp_104 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1023 = atomicAdd(v_coeffs_17 + int(19), temp_104);
    }
    float temp_105 = pSH8_11 * _S1013;
    if((F32_isfinite((temp_105))))
    {
        _S990 = temp_105 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1024 = atomicAdd(v_coeffs_17 + int(22), temp_105);
    }
    float v_x_10 = _S1017 + _S1013 * (pSH4_x_3 * *_S960 + pSH8_x_9 * *_S964 + fTmp0B_14 * *_S963);
    float v_y_10 = _S1018 + _S1013 * (pSH8_x_9 * *_S960 + pSH8_y_3 * *_S964 + fTmp0B_14 * *_S961);
    float v_z_10 = _S1019 + _S1013 * (pSH6_z_5 * *_S962 + pSH7_z_3 * *_S963 + pSH5_z_3 * *_S961);
    float temp_106 = pSH9_6 * _S1013;
    if((F32_isfinite((temp_106))))
    {
        _S990 = temp_106 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1025 = atomicAdd(v_coeffs_17 + int(25), temp_106);
    }
    float temp_107 = pSH10_6 * _S1013;
    if((F32_isfinite((temp_107))))
    {
        _S990 = temp_107 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1026 = atomicAdd(v_coeffs_17 + int(28), temp_107);
    }
    float temp_108 = pSH11_6 * _S1013;
    if((F32_isfinite((temp_108))))
    {
        _S990 = temp_108 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1027 = atomicAdd(v_coeffs_17 + int(31), temp_108);
    }
    float temp_109 = pSH12_7 * _S1013;
    if((F32_isfinite((temp_109))))
    {
        _S990 = temp_109 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1028 = atomicAdd(v_coeffs_17 + int(34), temp_109);
    }
    float temp_110 = pSH13_6 * _S1013;
    if((F32_isfinite((temp_110))))
    {
        _S990 = temp_110 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1029 = atomicAdd(v_coeffs_17 + int(37), temp_110);
    }
    float temp_111 = pSH14_6 * _S1013;
    if((F32_isfinite((temp_111))))
    {
        _S990 = temp_111 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1030 = atomicAdd(v_coeffs_17 + int(40), temp_111);
    }
    float temp_112 = pSH15_6 * _S1013;
    if((F32_isfinite((temp_112))))
    {
        _S990 = temp_112 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1031 = atomicAdd(v_coeffs_17 + int(43), temp_112);
    }
    float3  v_dir_n_24 = make_float3 (v_x_10 + _S1013 * (pSH9_x_1 * *_S965 + pSH15_x_1 * *_S971 + pSH10_x_1 * *_S966 + pSH14_x_5 * *_S970 + fTmp0C_8 * *_S969), v_y_10 + _S1013 * (pSH9_y_1 * *_S965 + pSH15_y_1 * *_S971 + pSH14_x_5 * *_S966 + pSH14_y_1 * *_S970 + fTmp0C_8 * *_S967), v_z_10 + _S1013 * (pSH12_z_3 * *_S968 + pSH13_z_1 * *_S969 + pSH11_z_1 * *_S967 + pSH14_z_1 * *_S970 + pSH10_z_1 * *_S966));
    float3  v_viewdir_42 = v_viewdir_41 + (v_dir_n_24 - make_float3 (dot_0(v_dir_n_24, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    float _S1032 = _S987.z;
    *&(v_coeff_dc_17->z) = *&(v_coeff_dc_17->z) + 0.282094806432724f * _S1032;
    float temp_113 = _S989 * _S1032;
    if((F32_isfinite((temp_113))))
    {
        _S990 = temp_113 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1033 = atomicAdd(v_coeffs_17 + int(2), temp_113);
    }
    float temp_114 = _S992 * _S1032;
    if((F32_isfinite((temp_114))))
    {
        _S990 = temp_114 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1034 = atomicAdd(v_coeffs_17 + int(5), temp_114);
    }
    float temp_115 = _S994 * _S1032;
    if((F32_isfinite((temp_115))))
    {
        _S990 = temp_115 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1035 = atomicAdd(v_coeffs_17 + int(8), temp_115);
    }
    float _S1036 = -0.48860251903533936f * *_S974 * _S1032;
    float _S1037 = -0.48860251903533936f * *_S972 * _S1032;
    float _S1038 = 0.48860251903533936f * *_S973 * _S1032;
    float temp_116 = pSH4_11 * _S1032;
    if((F32_isfinite((temp_116))))
    {
        _S990 = temp_116 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1039 = atomicAdd(v_coeffs_17 + int(11), temp_116);
    }
    float temp_117 = pSH5_11 * _S1032;
    if((F32_isfinite((temp_117))))
    {
        _S990 = temp_117 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1040 = atomicAdd(v_coeffs_17 + int(14), temp_117);
    }
    float temp_118 = pSH6_12 * _S1032;
    if((F32_isfinite((temp_118))))
    {
        _S990 = temp_118 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1041 = atomicAdd(v_coeffs_17 + int(17), temp_118);
    }
    float temp_119 = pSH7_11 * _S1032;
    if((F32_isfinite((temp_119))))
    {
        _S990 = temp_119 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1042 = atomicAdd(v_coeffs_17 + int(20), temp_119);
    }
    float temp_120 = pSH8_11 * _S1032;
    if((F32_isfinite((temp_120))))
    {
        _S990 = temp_120 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1043 = atomicAdd(v_coeffs_17 + int(23), temp_120);
    }
    float v_x_11 = _S1036 + _S1032 * (pSH4_x_3 * *_S975 + pSH8_x_9 * *_S979 + fTmp0B_14 * *_S978);
    float v_y_11 = _S1037 + _S1032 * (pSH8_x_9 * *_S975 + pSH8_y_3 * *_S979 + fTmp0B_14 * *_S976);
    float v_z_11 = _S1038 + _S1032 * (pSH6_z_5 * *_S977 + pSH7_z_3 * *_S978 + pSH5_z_3 * *_S976);
    float temp_121 = pSH9_6 * _S1032;
    if((F32_isfinite((temp_121))))
    {
        _S990 = temp_121 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1044 = atomicAdd(v_coeffs_17 + int(26), temp_121);
    }
    float temp_122 = pSH10_6 * _S1032;
    if((F32_isfinite((temp_122))))
    {
        _S990 = temp_122 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1045 = atomicAdd(v_coeffs_17 + int(29), temp_122);
    }
    float temp_123 = pSH11_6 * _S1032;
    if((F32_isfinite((temp_123))))
    {
        _S990 = temp_123 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1046 = atomicAdd(v_coeffs_17 + int(32), temp_123);
    }
    float temp_124 = pSH12_7 * _S1032;
    if((F32_isfinite((temp_124))))
    {
        _S990 = temp_124 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1047 = atomicAdd(v_coeffs_17 + int(35), temp_124);
    }
    float temp_125 = pSH13_6 * _S1032;
    if((F32_isfinite((temp_125))))
    {
        _S990 = temp_125 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1048 = atomicAdd(v_coeffs_17 + int(38), temp_125);
    }
    float temp_126 = pSH14_6 * _S1032;
    if((F32_isfinite((temp_126))))
    {
        _S990 = temp_126 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1049 = atomicAdd(v_coeffs_17 + int(41), temp_126);
    }
    float temp_127 = pSH15_6 * _S1032;
    if((F32_isfinite((temp_127))))
    {
        _S990 = temp_127 != 0.0f;
    }
    else
    {
        _S990 = false;
    }
    if(_S990)
    {
        float _S1050 = atomicAdd(v_coeffs_17 + int(44), temp_127);
    }
    float3  v_dir_n_25 = make_float3 (v_x_11 + _S1032 * (pSH9_x_1 * *_S980 + pSH15_x_1 * *_S986 + pSH10_x_1 * *_S981 + pSH14_x_5 * *_S985 + fTmp0C_8 * *_S984), v_y_11 + _S1032 * (pSH9_y_1 * *_S980 + pSH15_y_1 * *_S986 + pSH14_x_5 * *_S981 + pSH14_y_1 * *_S985 + fTmp0C_8 * *_S982), v_z_11 + _S1032 * (pSH12_z_3 * *_S983 + pSH13_z_1 * *_S984 + pSH11_z_1 * *_S982 + pSH14_z_1 * *_S985 + pSH10_z_1 * *_S981));
    float3  v_viewdir_43 = v_viewdir_42 + (v_dir_n_25 - make_float3 (dot_0(v_dir_n_25, dir_n_13)) * dir_n_13) * make_float3 (inv_norm_20);
    Matrix<float, 3, 3>  _S1051 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1052;
    (&_S1052)->primal_0 = _S935;
    (&_S1052)->differential_0 = _S1051;
    float3  _S1053 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1054;
    (&_S1054)->primal_0 = t_26;
    (&_S1054)->differential_0 = _S1053;
    s_bwd_prop_mul_0(&_S1052, &_S1054, v_viewdir_43);
    Matrix<float, 3, 3>  _S1055 = transpose_0(_S1052.differential_0);
    *v_mean_17 = *v_mean_17 + v_viewdir_43;
    *v_R_17 = *v_R_17 + _S1055;
    *v_t_17 = *v_t_17 + _S1054.differential_0;
    return;
}

inline __device__ float3  sh4_to_color(float3  mean_27, Matrix<float, 3, 3>  R_27, float3  t_27, float3  coeff_dc_27, float * coeffs_27)
{
    float3  _S1056 = mean_27 + mul_0(transpose_0(R_27), t_27);
    float _S1057 = _S1056.x;
    float _S1058 = _S1056.y;
    float _S1059 = _S1056.z;
    float inv_norm_21 = (F32_rsqrt((_S1057 * _S1057 + _S1058 * _S1058 + _S1059 * _S1059)));
    float x_23 = _S1057 * inv_norm_21;
    float y_21 = _S1058 * inv_norm_21;
    float z_20 = _S1059 * inv_norm_21;
    float _S1060 = - y_21;
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
    return max_0(make_float3 (0.282094806432724f * coeff_dc_27.x + 0.48860251903533936f * (_S1060 * *(coeffs_27 + int(0)) + z_20 * *(coeffs_27 + int(3)) - x_23 * *(coeffs_27 + int(6))) + (pSH4_12 * *(coeffs_27 + int(9)) + pSH5_12 * *(coeffs_27 + int(12)) + pSH6_13 * *(coeffs_27 + int(15)) + pSH7_12 * *(coeffs_27 + int(18)) + pSH8_12 * *(coeffs_27 + int(21))) + (pSH9_7 * *(coeffs_27 + int(24)) + pSH10_7 * *(coeffs_27 + int(27)) + pSH11_7 * *(coeffs_27 + int(30)) + pSH12_8 * *(coeffs_27 + int(33)) + pSH13_7 * *(coeffs_27 + int(36)) + pSH14_7 * *(coeffs_27 + int(39)) + pSH15_7 * *(coeffs_27 + int(42))) + (pSH16_2 * *(coeffs_27 + int(45)) + pSH17_2 * *(coeffs_27 + int(48)) + pSH18_2 * *(coeffs_27 + int(51)) + pSH19_2 * *(coeffs_27 + int(54)) + pSH20_2 * *(coeffs_27 + int(57)) + pSH21_2 * *(coeffs_27 + int(60)) + pSH22_2 * *(coeffs_27 + int(63)) + pSH23_2 * *(coeffs_27 + int(66)) + pSH24_2 * *(coeffs_27 + int(69))), 0.282094806432724f * coeff_dc_27.y + 0.48860251903533936f * (_S1060 * *(coeffs_27 + int(1)) + z_20 * *(coeffs_27 + int(4)) - x_23 * *(coeffs_27 + int(7))) + (pSH4_12 * *(coeffs_27 + int(10)) + pSH5_12 * *(coeffs_27 + int(13)) + pSH6_13 * *(coeffs_27 + int(16)) + pSH7_12 * *(coeffs_27 + int(19)) + pSH8_12 * *(coeffs_27 + int(22))) + (pSH9_7 * *(coeffs_27 + int(25)) + pSH10_7 * *(coeffs_27 + int(28)) + pSH11_7 * *(coeffs_27 + int(31)) + pSH12_8 * *(coeffs_27 + int(34)) + pSH13_7 * *(coeffs_27 + int(37)) + pSH14_7 * *(coeffs_27 + int(40)) + pSH15_7 * *(coeffs_27 + int(43))) + (pSH16_2 * *(coeffs_27 + int(46)) + pSH17_2 * *(coeffs_27 + int(49)) + pSH18_2 * *(coeffs_27 + int(52)) + pSH19_2 * *(coeffs_27 + int(55)) + pSH20_2 * *(coeffs_27 + int(58)) + pSH21_2 * *(coeffs_27 + int(61)) + pSH22_2 * *(coeffs_27 + int(64)) + pSH23_2 * *(coeffs_27 + int(67)) + pSH24_2 * *(coeffs_27 + int(70))), 0.282094806432724f * coeff_dc_27.z + 0.48860251903533936f * (_S1060 * *(coeffs_27 + int(2)) + z_20 * *(coeffs_27 + int(5)) - x_23 * *(coeffs_27 + int(8))) + (pSH4_12 * *(coeffs_27 + int(11)) + pSH5_12 * *(coeffs_27 + int(14)) + pSH6_13 * *(coeffs_27 + int(17)) + pSH7_12 * *(coeffs_27 + int(20)) + pSH8_12 * *(coeffs_27 + int(23))) + (pSH9_7 * *(coeffs_27 + int(26)) + pSH10_7 * *(coeffs_27 + int(29)) + pSH11_7 * *(coeffs_27 + int(32)) + pSH12_8 * *(coeffs_27 + int(35)) + pSH13_7 * *(coeffs_27 + int(38)) + pSH14_7 * *(coeffs_27 + int(41)) + pSH15_7 * *(coeffs_27 + int(44))) + (pSH16_2 * *(coeffs_27 + int(47)) + pSH17_2 * *(coeffs_27 + int(50)) + pSH18_2 * *(coeffs_27 + int(53)) + pSH19_2 * *(coeffs_27 + int(56)) + pSH20_2 * *(coeffs_27 + int(59)) + pSH21_2 * *(coeffs_27 + int(62)) + pSH22_2 * *(coeffs_27 + int(65)) + pSH23_2 * *(coeffs_27 + int(68)) + pSH24_2 * *(coeffs_27 + int(71)))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh4_to_color_vjp_inplace(float3  mean_28, Matrix<float, 3, 3>  R_28, float3  t_28, float3  coeff_dc_28, float * coeffs_28, float3  v_colors_18, float3  * v_coeff_dc_18, float * v_coeffs_18, float3  * v_mean_18, Matrix<float, 3, 3>  * v_R_18, float3  * v_t_18)
{
    Matrix<float, 3, 3>  _S1061 = transpose_0(R_28);
    float3  _S1062 = mean_28 + mul_0(_S1061, t_28);
    float _S1063 = _S1062.x;
    float _S1064 = _S1062.y;
    float _S1065 = _S1062.z;
    float inv_norm_22 = (F32_rsqrt((_S1063 * _S1063 + _S1064 * _S1064 + _S1065 * _S1065)));
    float x_24 = _S1063 * inv_norm_22;
    float y_22 = _S1064 * inv_norm_22;
    float z_21 = _S1065 * inv_norm_22;
    float _S1066 = - y_22;
    float * _S1067 = coeffs_28 + int(0);
    float * _S1068 = coeffs_28 + int(3);
    float * _S1069 = coeffs_28 + int(6);
    float z2_10 = z_21 * z_21;
    float fTmp0B_16 = -1.09254848957061768f * z_21;
    float fC1_10 = x_24 * x_24 - y_22 * y_22;
    float _S1070 = 2.0f * x_24;
    float fS1_10 = _S1070 * y_22;
    float pSH6_14 = 0.94617468118667603f * z2_10 - 0.31539157032966614f;
    float pSH7_13 = fTmp0B_16 * x_24;
    float pSH5_13 = fTmp0B_16 * y_22;
    float pSH8_13 = 0.54627424478530884f * fC1_10;
    float pSH4_13 = 0.54627424478530884f * fS1_10;
    float * _S1071 = coeffs_28 + int(9);
    float * _S1072 = coeffs_28 + int(12);
    float * _S1073 = coeffs_28 + int(15);
    float * _S1074 = coeffs_28 + int(18);
    float * _S1075 = coeffs_28 + int(21);
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
    float * _S1076 = coeffs_28 + int(24);
    float * _S1077 = coeffs_28 + int(27);
    float * _S1078 = coeffs_28 + int(30);
    float * _S1079 = coeffs_28 + int(33);
    float * _S1080 = coeffs_28 + int(36);
    float * _S1081 = coeffs_28 + int(39);
    float * _S1082 = coeffs_28 + int(42);
    float fTmp0D_4 = z_21 * (-4.68332576751708984f * z2_10 + 2.00713968276977539f);
    float fTmp1C_4 = 3.31161141395568848f * z2_10 - 0.47308734059333801f;
    float fTmp2B_4 = -1.77013075351715088f * z_21;
    float _S1083 = 1.9843134880065918f * z_21 * pSH12_9;
    float pSH20_3 = _S1083 - 1.00623059272766113f * pSH6_14;
    float pSH21_3 = fTmp0D_4 * x_24;
    float pSH19_3 = fTmp0D_4 * y_22;
    float pSH22_3 = fTmp1C_4 * fC1_10;
    float pSH18_3 = fTmp1C_4 * fS1_10;
    float pSH23_3 = fTmp2B_4 * fC2_4;
    float pSH17_3 = fTmp2B_4 * fS2_4;
    float pSH24_3 = 0.62583571672439575f * (x_24 * fC2_4 - y_22 * fS2_4);
    float pSH16_3 = 0.62583571672439575f * (x_24 * fS2_4 + y_22 * fC2_4);
    float * _S1084 = coeffs_28 + int(45);
    float * _S1085 = coeffs_28 + int(48);
    float * _S1086 = coeffs_28 + int(51);
    float * _S1087 = coeffs_28 + int(54);
    float * _S1088 = coeffs_28 + int(57);
    float * _S1089 = coeffs_28 + int(60);
    float * _S1090 = coeffs_28 + int(63);
    float * _S1091 = coeffs_28 + int(66);
    float * _S1092 = coeffs_28 + int(69);
    float * _S1093 = coeffs_28 + int(1);
    float * _S1094 = coeffs_28 + int(4);
    float * _S1095 = coeffs_28 + int(7);
    float * _S1096 = coeffs_28 + int(10);
    float * _S1097 = coeffs_28 + int(13);
    float * _S1098 = coeffs_28 + int(16);
    float * _S1099 = coeffs_28 + int(19);
    float * _S1100 = coeffs_28 + int(22);
    float * _S1101 = coeffs_28 + int(25);
    float * _S1102 = coeffs_28 + int(28);
    float * _S1103 = coeffs_28 + int(31);
    float * _S1104 = coeffs_28 + int(34);
    float * _S1105 = coeffs_28 + int(37);
    float * _S1106 = coeffs_28 + int(40);
    float * _S1107 = coeffs_28 + int(43);
    float * _S1108 = coeffs_28 + int(46);
    float * _S1109 = coeffs_28 + int(49);
    float * _S1110 = coeffs_28 + int(52);
    float * _S1111 = coeffs_28 + int(55);
    float * _S1112 = coeffs_28 + int(58);
    float * _S1113 = coeffs_28 + int(61);
    float * _S1114 = coeffs_28 + int(64);
    float * _S1115 = coeffs_28 + int(67);
    float * _S1116 = coeffs_28 + int(70);
    float * _S1117 = coeffs_28 + int(2);
    float * _S1118 = coeffs_28 + int(5);
    float * _S1119 = coeffs_28 + int(8);
    float * _S1120 = coeffs_28 + int(11);
    float * _S1121 = coeffs_28 + int(14);
    float * _S1122 = coeffs_28 + int(17);
    float * _S1123 = coeffs_28 + int(20);
    float * _S1124 = coeffs_28 + int(23);
    float * _S1125 = coeffs_28 + int(26);
    float * _S1126 = coeffs_28 + int(29);
    float * _S1127 = coeffs_28 + int(32);
    float * _S1128 = coeffs_28 + int(35);
    float * _S1129 = coeffs_28 + int(38);
    float * _S1130 = coeffs_28 + int(41);
    float * _S1131 = coeffs_28 + int(44);
    float * _S1132 = coeffs_28 + int(47);
    float * _S1133 = coeffs_28 + int(50);
    float * _S1134 = coeffs_28 + int(53);
    float * _S1135 = coeffs_28 + int(56);
    float * _S1136 = coeffs_28 + int(59);
    float * _S1137 = coeffs_28 + int(62);
    float * _S1138 = coeffs_28 + int(65);
    float * _S1139 = coeffs_28 + int(68);
    float * _S1140 = coeffs_28 + int(71);
    float3  _S1141 = v_colors_18 * make_float3 (float((0.282094806432724f * coeff_dc_28.x + 0.48860251903533936f * (_S1066 * *_S1067 + z_21 * *_S1068 - x_24 * *_S1069) + (pSH4_13 * *_S1071 + pSH5_13 * *_S1072 + pSH6_14 * *_S1073 + pSH7_13 * *_S1074 + pSH8_13 * *_S1075) + (pSH9_8 * *_S1076 + pSH10_8 * *_S1077 + pSH11_8 * *_S1078 + pSH12_9 * *_S1079 + pSH13_8 * *_S1080 + pSH14_8 * *_S1081 + pSH15_8 * *_S1082) + (pSH16_3 * *_S1084 + pSH17_3 * *_S1085 + pSH18_3 * *_S1086 + pSH19_3 * *_S1087 + pSH20_3 * *_S1088 + pSH21_3 * *_S1089 + pSH22_3 * *_S1090 + pSH23_3 * *_S1091 + pSH24_3 * *_S1092)) >= -0.5f), float((0.282094806432724f * coeff_dc_28.y + 0.48860251903533936f * (_S1066 * *_S1093 + z_21 * *_S1094 - x_24 * *_S1095) + (pSH4_13 * *_S1096 + pSH5_13 * *_S1097 + pSH6_14 * *_S1098 + pSH7_13 * *_S1099 + pSH8_13 * *_S1100) + (pSH9_8 * *_S1101 + pSH10_8 * *_S1102 + pSH11_8 * *_S1103 + pSH12_9 * *_S1104 + pSH13_8 * *_S1105 + pSH14_8 * *_S1106 + pSH15_8 * *_S1107) + (pSH16_3 * *_S1108 + pSH17_3 * *_S1109 + pSH18_3 * *_S1110 + pSH19_3 * *_S1111 + pSH20_3 * *_S1112 + pSH21_3 * *_S1113 + pSH22_3 * *_S1114 + pSH23_3 * *_S1115 + pSH24_3 * *_S1116)) >= -0.5f), float((0.282094806432724f * coeff_dc_28.z + 0.48860251903533936f * (_S1066 * *_S1117 + z_21 * *_S1118 - x_24 * *_S1119) + (pSH4_13 * *_S1120 + pSH5_13 * *_S1121 + pSH6_14 * *_S1122 + pSH7_13 * *_S1123 + pSH8_13 * *_S1124) + (pSH9_8 * *_S1125 + pSH10_8 * *_S1126 + pSH11_8 * *_S1127 + pSH12_9 * *_S1128 + pSH13_8 * *_S1129 + pSH14_8 * *_S1130 + pSH15_8 * *_S1131) + (pSH16_3 * *_S1132 + pSH17_3 * *_S1133 + pSH18_3 * *_S1134 + pSH19_3 * *_S1135 + pSH20_3 * *_S1136 + pSH21_3 * *_S1137 + pSH22_3 * *_S1138 + pSH23_3 * *_S1139 + pSH24_3 * *_S1140)) >= -0.5f));
    float3  v_viewdir_44 = {};
    float _S1142 = _S1141.x;
    *&(v_coeff_dc_18->x) = *&(v_coeff_dc_18->x) + 0.282094806432724f * _S1142;
    float * _S1143 = v_coeffs_18 + int(0);
    float _S1144 = -0.48860251903533936f * y_22;
    *_S1143 = *_S1143 + _S1144 * _S1142;
    float * _S1145 = v_coeffs_18 + int(3);
    float _S1146 = 0.48860251903533936f * z_21;
    *_S1145 = *_S1145 + _S1146 * _S1142;
    float * _S1147 = v_coeffs_18 + int(6);
    float _S1148 = -0.48860251903533936f * x_24;
    *_S1147 = *_S1147 + _S1148 * _S1142;
    float _S1149 = -0.48860251903533936f * *_S1069 * _S1142;
    float _S1150 = -0.48860251903533936f * *_S1067 * _S1142;
    float _S1151 = 0.48860251903533936f * *_S1068 * _S1142;
    float * _S1152 = v_coeffs_18 + int(9);
    *_S1152 = *_S1152 + pSH4_13 * _S1142;
    float * _S1153 = v_coeffs_18 + int(12);
    *_S1153 = *_S1153 + pSH5_13 * _S1142;
    float * _S1154 = v_coeffs_18 + int(15);
    *_S1154 = *_S1154 + pSH6_14 * _S1142;
    float * _S1155 = v_coeffs_18 + int(18);
    *_S1155 = *_S1155 + pSH7_13 * _S1142;
    float * _S1156 = v_coeffs_18 + int(21);
    *_S1156 = *_S1156 + pSH8_13 * _S1142;
    float fC1_y_6 = -2.0f * y_22;
    float fS1_x_6 = 2.0f * y_22;
    float pSH6_z_6 = 1.89234936237335205f * z_21;
    float pSH7_z_4 = -1.09254848957061768f * x_24;
    float pSH5_z_4 = -1.09254848957061768f * y_22;
    float pSH8_x_10 = 0.54627424478530884f * _S1070;
    float pSH8_y_4 = 0.54627424478530884f * fC1_y_6;
    float pSH4_x_4 = 0.54627424478530884f * fS1_x_6;
    float v_x_12 = _S1149 + _S1142 * (pSH4_x_4 * *_S1071 + pSH8_x_10 * *_S1075 + fTmp0B_16 * *_S1074);
    float v_y_12 = _S1150 + _S1142 * (pSH8_x_10 * *_S1071 + pSH8_y_4 * *_S1075 + fTmp0B_16 * *_S1072);
    float v_z_12 = _S1151 + _S1142 * (pSH6_z_6 * *_S1073 + pSH7_z_4 * *_S1074 + pSH5_z_4 * *_S1072);
    float * _S1157 = v_coeffs_18 + int(24);
    *_S1157 = *_S1157 + pSH9_8 * _S1142;
    float * _S1158 = v_coeffs_18 + int(27);
    *_S1158 = *_S1158 + pSH10_8 * _S1142;
    float * _S1159 = v_coeffs_18 + int(30);
    *_S1159 = *_S1159 + pSH11_8 * _S1142;
    float * _S1160 = v_coeffs_18 + int(33);
    *_S1160 = *_S1160 + pSH12_9 * _S1142;
    float * _S1161 = v_coeffs_18 + int(36);
    *_S1161 = *_S1161 + pSH13_8 * _S1142;
    float * _S1162 = v_coeffs_18 + int(39);
    *_S1162 = *_S1162 + pSH14_8 * _S1142;
    float * _S1163 = v_coeffs_18 + int(42);
    *_S1163 = *_S1163 + pSH15_8 * _S1142;
    float fTmp0C_z_6 = -4.57045793533325195f * z_21;
    float _S1164 = x_24 * _S1070;
    float fC2_x_2 = fC1_10 + _S1164 - y_22 * fS1_x_6;
    float _S1165 = y_22 * _S1070;
    float fC2_y_2 = x_24 * fC1_y_6 - fS1_10 - _S1165;
    float fS2_x_2 = fS1_10 + x_24 * fS1_x_6 + _S1165;
    float fS2_y_2 = _S1164 + fC1_10 + y_22 * fC1_y_6;
    float pSH12_z_4 = 5.59764480590820312f * z2_10 - 1.11952900886535645f;
    float pSH13_z_2 = fTmp0C_z_6 * x_24;
    float pSH11_z_2 = fTmp0C_z_6 * y_22;
    float pSH14_x_6 = fTmp1B_10 * _S1070;
    float pSH14_y_2 = fTmp1B_10 * fC1_y_6;
    float pSH14_z_2 = 1.44530570507049561f * fC1_10;
    float pSH10_x_2 = fTmp1B_10 * fS1_x_6;
    float pSH10_z_2 = 1.44530570507049561f * fS1_10;
    float pSH15_x_2 = -0.59004360437393188f * fC2_x_2;
    float pSH15_y_2 = -0.59004360437393188f * fC2_y_2;
    float pSH9_x_2 = -0.59004360437393188f * fS2_x_2;
    float pSH9_y_2 = -0.59004360437393188f * fS2_y_2;
    float v_x_13 = v_x_12 + _S1142 * (pSH9_x_2 * *_S1076 + pSH15_x_2 * *_S1082 + pSH10_x_2 * *_S1077 + pSH14_x_6 * *_S1081 + fTmp0C_10 * *_S1080);
    float v_y_13 = v_y_12 + _S1142 * (pSH9_y_2 * *_S1076 + pSH15_y_2 * *_S1082 + pSH14_x_6 * *_S1077 + pSH14_y_2 * *_S1081 + fTmp0C_10 * *_S1078);
    float v_z_13 = v_z_12 + _S1142 * (pSH12_z_4 * *_S1079 + pSH13_z_2 * *_S1080 + pSH11_z_2 * *_S1078 + pSH14_z_2 * *_S1081 + pSH10_z_2 * *_S1077);
    float pSH20_4 = _S1083 + -1.00623059272766113f * pSH6_14;
    float * _S1166 = v_coeffs_18 + int(45);
    *_S1166 = *_S1166 + pSH16_3 * _S1142;
    float * _S1167 = v_coeffs_18 + int(48);
    *_S1167 = *_S1167 + pSH17_3 * _S1142;
    float * _S1168 = v_coeffs_18 + int(51);
    *_S1168 = *_S1168 + pSH18_3 * _S1142;
    float * _S1169 = v_coeffs_18 + int(54);
    *_S1169 = *_S1169 + pSH19_3 * _S1142;
    float * _S1170 = v_coeffs_18 + int(57);
    *_S1170 = *_S1170 + pSH20_4 * _S1142;
    float * _S1171 = v_coeffs_18 + int(60);
    *_S1171 = *_S1171 + pSH21_3 * _S1142;
    float * _S1172 = v_coeffs_18 + int(63);
    *_S1172 = *_S1172 + pSH22_3 * _S1142;
    float * _S1173 = v_coeffs_18 + int(66);
    *_S1173 = *_S1173 + pSH23_3 * _S1142;
    float * _S1174 = v_coeffs_18 + int(69);
    *_S1174 = *_S1174 + pSH24_3 * _S1142;
    float fTmp0D_z_2 = -14.04997730255126953f * z2_10 + 2.00713968276977539f;
    float fTmp1C_z_2 = 6.62322282791137695f * z_21;
    float pSH20_z_0 = 1.9843134880065918f * (pSH12_9 + z_21 * pSH12_z_4) + -1.00623059272766113f * pSH6_z_6;
    float pSH21_z_0 = fTmp0D_z_2 * x_24;
    float pSH19_z_0 = fTmp0D_z_2 * y_22;
    float pSH22_x_2 = fTmp1C_4 * _S1070;
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
    float3  v_dir_n_26 = make_float3 (v_x_13 + _S1142 * (pSH16_x_0 * *_S1084 + pSH24_x_0 * *_S1092 + pSH17_x_0 * *_S1085 + pSH23_x_0 * *_S1091 + pSH18_x_0 * *_S1086 + pSH22_x_2 * *_S1090 + fTmp0D_4 * *_S1089), v_y_13 + _S1142 * (pSH16_y_0 * *_S1084 + pSH24_y_0 * *_S1092 + pSH17_y_0 * *_S1085 + pSH23_y_0 * *_S1091 + pSH22_x_2 * *_S1086 + pSH22_y_0 * *_S1090 + fTmp0D_4 * *_S1087), v_z_13 + _S1142 * (pSH20_z_0 * *_S1088 + pSH21_z_0 * *_S1089 + pSH19_z_0 * *_S1087 + pSH22_z_0 * *_S1090 + pSH18_z_0 * *_S1086 + pSH23_z_0 * *_S1091 + pSH17_z_0 * *_S1085));
    float3  v_viewdir_45 = v_viewdir_44 + (v_dir_n_26 - make_float3 (dot_0(v_dir_n_26, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    float _S1175 = _S1141.y;
    *&(v_coeff_dc_18->y) = *&(v_coeff_dc_18->y) + 0.282094806432724f * _S1175;
    float * _S1176 = v_coeffs_18 + int(1);
    *_S1176 = *_S1176 + _S1144 * _S1175;
    float * _S1177 = v_coeffs_18 + int(4);
    *_S1177 = *_S1177 + _S1146 * _S1175;
    float * _S1178 = v_coeffs_18 + int(7);
    *_S1178 = *_S1178 + _S1148 * _S1175;
    float _S1179 = -0.48860251903533936f * *_S1095 * _S1175;
    float _S1180 = -0.48860251903533936f * *_S1093 * _S1175;
    float _S1181 = 0.48860251903533936f * *_S1094 * _S1175;
    float * _S1182 = v_coeffs_18 + int(10);
    *_S1182 = *_S1182 + pSH4_13 * _S1175;
    float * _S1183 = v_coeffs_18 + int(13);
    *_S1183 = *_S1183 + pSH5_13 * _S1175;
    float * _S1184 = v_coeffs_18 + int(16);
    *_S1184 = *_S1184 + pSH6_14 * _S1175;
    float * _S1185 = v_coeffs_18 + int(19);
    *_S1185 = *_S1185 + pSH7_13 * _S1175;
    float * _S1186 = v_coeffs_18 + int(22);
    *_S1186 = *_S1186 + pSH8_13 * _S1175;
    float v_x_14 = _S1179 + _S1175 * (pSH4_x_4 * *_S1096 + pSH8_x_10 * *_S1100 + fTmp0B_16 * *_S1099);
    float v_y_14 = _S1180 + _S1175 * (pSH8_x_10 * *_S1096 + pSH8_y_4 * *_S1100 + fTmp0B_16 * *_S1097);
    float v_z_14 = _S1181 + _S1175 * (pSH6_z_6 * *_S1098 + pSH7_z_4 * *_S1099 + pSH5_z_4 * *_S1097);
    float * _S1187 = v_coeffs_18 + int(25);
    *_S1187 = *_S1187 + pSH9_8 * _S1175;
    float * _S1188 = v_coeffs_18 + int(28);
    *_S1188 = *_S1188 + pSH10_8 * _S1175;
    float * _S1189 = v_coeffs_18 + int(31);
    *_S1189 = *_S1189 + pSH11_8 * _S1175;
    float * _S1190 = v_coeffs_18 + int(34);
    *_S1190 = *_S1190 + pSH12_9 * _S1175;
    float * _S1191 = v_coeffs_18 + int(37);
    *_S1191 = *_S1191 + pSH13_8 * _S1175;
    float * _S1192 = v_coeffs_18 + int(40);
    *_S1192 = *_S1192 + pSH14_8 * _S1175;
    float * _S1193 = v_coeffs_18 + int(43);
    *_S1193 = *_S1193 + pSH15_8 * _S1175;
    float v_x_15 = v_x_14 + _S1175 * (pSH9_x_2 * *_S1101 + pSH15_x_2 * *_S1107 + pSH10_x_2 * *_S1102 + pSH14_x_6 * *_S1106 + fTmp0C_10 * *_S1105);
    float v_y_15 = v_y_14 + _S1175 * (pSH9_y_2 * *_S1101 + pSH15_y_2 * *_S1107 + pSH14_x_6 * *_S1102 + pSH14_y_2 * *_S1106 + fTmp0C_10 * *_S1103);
    float v_z_15 = v_z_14 + _S1175 * (pSH12_z_4 * *_S1104 + pSH13_z_2 * *_S1105 + pSH11_z_2 * *_S1103 + pSH14_z_2 * *_S1106 + pSH10_z_2 * *_S1102);
    float * _S1194 = v_coeffs_18 + int(46);
    *_S1194 = *_S1194 + pSH16_3 * _S1175;
    float * _S1195 = v_coeffs_18 + int(49);
    *_S1195 = *_S1195 + pSH17_3 * _S1175;
    float * _S1196 = v_coeffs_18 + int(52);
    *_S1196 = *_S1196 + pSH18_3 * _S1175;
    float * _S1197 = v_coeffs_18 + int(55);
    *_S1197 = *_S1197 + pSH19_3 * _S1175;
    float * _S1198 = v_coeffs_18 + int(58);
    *_S1198 = *_S1198 + pSH20_4 * _S1175;
    float * _S1199 = v_coeffs_18 + int(61);
    *_S1199 = *_S1199 + pSH21_3 * _S1175;
    float * _S1200 = v_coeffs_18 + int(64);
    *_S1200 = *_S1200 + pSH22_3 * _S1175;
    float * _S1201 = v_coeffs_18 + int(67);
    *_S1201 = *_S1201 + pSH23_3 * _S1175;
    float * _S1202 = v_coeffs_18 + int(70);
    *_S1202 = *_S1202 + pSH24_3 * _S1175;
    float3  v_dir_n_27 = make_float3 (v_x_15 + _S1175 * (pSH16_x_0 * *_S1108 + pSH24_x_0 * *_S1116 + pSH17_x_0 * *_S1109 + pSH23_x_0 * *_S1115 + pSH18_x_0 * *_S1110 + pSH22_x_2 * *_S1114 + fTmp0D_4 * *_S1113), v_y_15 + _S1175 * (pSH16_y_0 * *_S1108 + pSH24_y_0 * *_S1116 + pSH17_y_0 * *_S1109 + pSH23_y_0 * *_S1115 + pSH22_x_2 * *_S1110 + pSH22_y_0 * *_S1114 + fTmp0D_4 * *_S1111), v_z_15 + _S1175 * (pSH20_z_0 * *_S1112 + pSH21_z_0 * *_S1113 + pSH19_z_0 * *_S1111 + pSH22_z_0 * *_S1114 + pSH18_z_0 * *_S1110 + pSH23_z_0 * *_S1115 + pSH17_z_0 * *_S1109));
    float3  v_viewdir_46 = v_viewdir_45 + (v_dir_n_27 - make_float3 (dot_0(v_dir_n_27, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    float _S1203 = _S1141.z;
    *&(v_coeff_dc_18->z) = *&(v_coeff_dc_18->z) + 0.282094806432724f * _S1203;
    float * _S1204 = v_coeffs_18 + int(2);
    *_S1204 = *_S1204 + _S1144 * _S1203;
    float * _S1205 = v_coeffs_18 + int(5);
    *_S1205 = *_S1205 + _S1146 * _S1203;
    float * _S1206 = v_coeffs_18 + int(8);
    *_S1206 = *_S1206 + _S1148 * _S1203;
    float _S1207 = -0.48860251903533936f * *_S1119 * _S1203;
    float _S1208 = -0.48860251903533936f * *_S1117 * _S1203;
    float _S1209 = 0.48860251903533936f * *_S1118 * _S1203;
    float * _S1210 = v_coeffs_18 + int(11);
    *_S1210 = *_S1210 + pSH4_13 * _S1203;
    float * _S1211 = v_coeffs_18 + int(14);
    *_S1211 = *_S1211 + pSH5_13 * _S1203;
    float * _S1212 = v_coeffs_18 + int(17);
    *_S1212 = *_S1212 + pSH6_14 * _S1203;
    float * _S1213 = v_coeffs_18 + int(20);
    *_S1213 = *_S1213 + pSH7_13 * _S1203;
    float * _S1214 = v_coeffs_18 + int(23);
    *_S1214 = *_S1214 + pSH8_13 * _S1203;
    float v_x_16 = _S1207 + _S1203 * (pSH4_x_4 * *_S1120 + pSH8_x_10 * *_S1124 + fTmp0B_16 * *_S1123);
    float v_y_16 = _S1208 + _S1203 * (pSH8_x_10 * *_S1120 + pSH8_y_4 * *_S1124 + fTmp0B_16 * *_S1121);
    float v_z_16 = _S1209 + _S1203 * (pSH6_z_6 * *_S1122 + pSH7_z_4 * *_S1123 + pSH5_z_4 * *_S1121);
    float * _S1215 = v_coeffs_18 + int(26);
    *_S1215 = *_S1215 + pSH9_8 * _S1203;
    float * _S1216 = v_coeffs_18 + int(29);
    *_S1216 = *_S1216 + pSH10_8 * _S1203;
    float * _S1217 = v_coeffs_18 + int(32);
    *_S1217 = *_S1217 + pSH11_8 * _S1203;
    float * _S1218 = v_coeffs_18 + int(35);
    *_S1218 = *_S1218 + pSH12_9 * _S1203;
    float * _S1219 = v_coeffs_18 + int(38);
    *_S1219 = *_S1219 + pSH13_8 * _S1203;
    float * _S1220 = v_coeffs_18 + int(41);
    *_S1220 = *_S1220 + pSH14_8 * _S1203;
    float * _S1221 = v_coeffs_18 + int(44);
    *_S1221 = *_S1221 + pSH15_8 * _S1203;
    float v_x_17 = v_x_16 + _S1203 * (pSH9_x_2 * *_S1125 + pSH15_x_2 * *_S1131 + pSH10_x_2 * *_S1126 + pSH14_x_6 * *_S1130 + fTmp0C_10 * *_S1129);
    float v_y_17 = v_y_16 + _S1203 * (pSH9_y_2 * *_S1125 + pSH15_y_2 * *_S1131 + pSH14_x_6 * *_S1126 + pSH14_y_2 * *_S1130 + fTmp0C_10 * *_S1127);
    float v_z_17 = v_z_16 + _S1203 * (pSH12_z_4 * *_S1128 + pSH13_z_2 * *_S1129 + pSH11_z_2 * *_S1127 + pSH14_z_2 * *_S1130 + pSH10_z_2 * *_S1126);
    float * _S1222 = v_coeffs_18 + int(47);
    *_S1222 = *_S1222 + pSH16_3 * _S1203;
    float * _S1223 = v_coeffs_18 + int(50);
    *_S1223 = *_S1223 + pSH17_3 * _S1203;
    float * _S1224 = v_coeffs_18 + int(53);
    *_S1224 = *_S1224 + pSH18_3 * _S1203;
    float * _S1225 = v_coeffs_18 + int(56);
    *_S1225 = *_S1225 + pSH19_3 * _S1203;
    float * _S1226 = v_coeffs_18 + int(59);
    *_S1226 = *_S1226 + pSH20_4 * _S1203;
    float * _S1227 = v_coeffs_18 + int(62);
    *_S1227 = *_S1227 + pSH21_3 * _S1203;
    float * _S1228 = v_coeffs_18 + int(65);
    *_S1228 = *_S1228 + pSH22_3 * _S1203;
    float * _S1229 = v_coeffs_18 + int(68);
    *_S1229 = *_S1229 + pSH23_3 * _S1203;
    float * _S1230 = v_coeffs_18 + int(71);
    *_S1230 = *_S1230 + pSH24_3 * _S1203;
    float3  v_dir_n_28 = make_float3 (v_x_17 + _S1203 * (pSH16_x_0 * *_S1132 + pSH24_x_0 * *_S1140 + pSH17_x_0 * *_S1133 + pSH23_x_0 * *_S1139 + pSH18_x_0 * *_S1134 + pSH22_x_2 * *_S1138 + fTmp0D_4 * *_S1137), v_y_17 + _S1203 * (pSH16_y_0 * *_S1132 + pSH24_y_0 * *_S1140 + pSH17_y_0 * *_S1133 + pSH23_y_0 * *_S1139 + pSH22_x_2 * *_S1134 + pSH22_y_0 * *_S1138 + fTmp0D_4 * *_S1135), v_z_17 + _S1203 * (pSH20_z_0 * *_S1136 + pSH21_z_0 * *_S1137 + pSH19_z_0 * *_S1135 + pSH22_z_0 * *_S1138 + pSH18_z_0 * *_S1134 + pSH23_z_0 * *_S1139 + pSH17_z_0 * *_S1133));
    float3  v_viewdir_47 = v_viewdir_46 + (v_dir_n_28 - make_float3 (dot_0(v_dir_n_28, dir_n_14)) * dir_n_14) * make_float3 (inv_norm_22);
    Matrix<float, 3, 3>  _S1231 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1232;
    (&_S1232)->primal_0 = _S1061;
    (&_S1232)->differential_0 = _S1231;
    float3  _S1233 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
    (&_S1234)->primal_0 = t_28;
    (&_S1234)->differential_0 = _S1233;
    s_bwd_prop_mul_0(&_S1232, &_S1234, v_viewdir_47);
    Matrix<float, 3, 3>  _S1235 = transpose_0(_S1232.differential_0);
    *v_mean_18 = *v_mean_18 + v_viewdir_47;
    *v_R_18 = *v_R_18 + _S1235;
    *v_t_18 = *v_t_18 + _S1234.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_atomic(float3  mean_29, Matrix<float, 3, 3>  R_29, float3  t_29, float3  coeff_dc_29, float * coeffs_29, float3  v_colors_19, float3  * v_coeff_dc_19, float * v_coeffs_19, float3  * v_mean_19, Matrix<float, 3, 3>  * v_R_19, float3  * v_t_19)
{
    Matrix<float, 3, 3>  _S1236 = transpose_0(R_29);
    float3  _S1237 = mean_29 + mul_0(_S1236, t_29);
    float _S1238 = _S1237.x;
    float _S1239 = _S1237.y;
    float _S1240 = _S1237.z;
    float inv_norm_23 = (F32_rsqrt((_S1238 * _S1238 + _S1239 * _S1239 + _S1240 * _S1240)));
    float x_25 = _S1238 * inv_norm_23;
    float y_23 = _S1239 * inv_norm_23;
    float z_22 = _S1240 * inv_norm_23;
    float _S1241 = - y_23;
    float * _S1242 = coeffs_29 + int(0);
    float * _S1243 = coeffs_29 + int(3);
    float * _S1244 = coeffs_29 + int(6);
    float z2_11 = z_22 * z_22;
    float fTmp0B_17 = -1.09254848957061768f * z_22;
    float fC1_11 = x_25 * x_25 - y_23 * y_23;
    float _S1245 = 2.0f * x_25;
    float fS1_11 = _S1245 * y_23;
    float pSH6_15 = 0.94617468118667603f * z2_11 - 0.31539157032966614f;
    float pSH7_14 = fTmp0B_17 * x_25;
    float pSH5_14 = fTmp0B_17 * y_23;
    float pSH8_14 = 0.54627424478530884f * fC1_11;
    float pSH4_14 = 0.54627424478530884f * fS1_11;
    float * _S1246 = coeffs_29 + int(9);
    float * _S1247 = coeffs_29 + int(12);
    float * _S1248 = coeffs_29 + int(15);
    float * _S1249 = coeffs_29 + int(18);
    float * _S1250 = coeffs_29 + int(21);
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
    float * _S1251 = coeffs_29 + int(24);
    float * _S1252 = coeffs_29 + int(27);
    float * _S1253 = coeffs_29 + int(30);
    float * _S1254 = coeffs_29 + int(33);
    float * _S1255 = coeffs_29 + int(36);
    float * _S1256 = coeffs_29 + int(39);
    float * _S1257 = coeffs_29 + int(42);
    float fTmp0D_5 = z_22 * (-4.68332576751708984f * z2_11 + 2.00713968276977539f);
    float fTmp1C_5 = 3.31161141395568848f * z2_11 - 0.47308734059333801f;
    float fTmp2B_5 = -1.77013075351715088f * z_22;
    float _S1258 = 1.9843134880065918f * z_22 * pSH12_10;
    float pSH20_5 = _S1258 - 1.00623059272766113f * pSH6_15;
    float pSH21_4 = fTmp0D_5 * x_25;
    float pSH19_4 = fTmp0D_5 * y_23;
    float pSH22_4 = fTmp1C_5 * fC1_11;
    float pSH18_4 = fTmp1C_5 * fS1_11;
    float pSH23_4 = fTmp2B_5 * fC2_5;
    float pSH17_4 = fTmp2B_5 * fS2_5;
    float pSH24_4 = 0.62583571672439575f * (x_25 * fC2_5 - y_23 * fS2_5);
    float pSH16_4 = 0.62583571672439575f * (x_25 * fS2_5 + y_23 * fC2_5);
    float * _S1259 = coeffs_29 + int(45);
    float * _S1260 = coeffs_29 + int(48);
    float * _S1261 = coeffs_29 + int(51);
    float * _S1262 = coeffs_29 + int(54);
    float * _S1263 = coeffs_29 + int(57);
    float * _S1264 = coeffs_29 + int(60);
    float * _S1265 = coeffs_29 + int(63);
    float * _S1266 = coeffs_29 + int(66);
    float * _S1267 = coeffs_29 + int(69);
    float * _S1268 = coeffs_29 + int(1);
    float * _S1269 = coeffs_29 + int(4);
    float * _S1270 = coeffs_29 + int(7);
    float * _S1271 = coeffs_29 + int(10);
    float * _S1272 = coeffs_29 + int(13);
    float * _S1273 = coeffs_29 + int(16);
    float * _S1274 = coeffs_29 + int(19);
    float * _S1275 = coeffs_29 + int(22);
    float * _S1276 = coeffs_29 + int(25);
    float * _S1277 = coeffs_29 + int(28);
    float * _S1278 = coeffs_29 + int(31);
    float * _S1279 = coeffs_29 + int(34);
    float * _S1280 = coeffs_29 + int(37);
    float * _S1281 = coeffs_29 + int(40);
    float * _S1282 = coeffs_29 + int(43);
    float * _S1283 = coeffs_29 + int(46);
    float * _S1284 = coeffs_29 + int(49);
    float * _S1285 = coeffs_29 + int(52);
    float * _S1286 = coeffs_29 + int(55);
    float * _S1287 = coeffs_29 + int(58);
    float * _S1288 = coeffs_29 + int(61);
    float * _S1289 = coeffs_29 + int(64);
    float * _S1290 = coeffs_29 + int(67);
    float * _S1291 = coeffs_29 + int(70);
    float * _S1292 = coeffs_29 + int(2);
    float * _S1293 = coeffs_29 + int(5);
    float * _S1294 = coeffs_29 + int(8);
    float * _S1295 = coeffs_29 + int(11);
    float * _S1296 = coeffs_29 + int(14);
    float * _S1297 = coeffs_29 + int(17);
    float * _S1298 = coeffs_29 + int(20);
    float * _S1299 = coeffs_29 + int(23);
    float * _S1300 = coeffs_29 + int(26);
    float * _S1301 = coeffs_29 + int(29);
    float * _S1302 = coeffs_29 + int(32);
    float * _S1303 = coeffs_29 + int(35);
    float * _S1304 = coeffs_29 + int(38);
    float * _S1305 = coeffs_29 + int(41);
    float * _S1306 = coeffs_29 + int(44);
    float * _S1307 = coeffs_29 + int(47);
    float * _S1308 = coeffs_29 + int(50);
    float * _S1309 = coeffs_29 + int(53);
    float * _S1310 = coeffs_29 + int(56);
    float * _S1311 = coeffs_29 + int(59);
    float * _S1312 = coeffs_29 + int(62);
    float * _S1313 = coeffs_29 + int(65);
    float * _S1314 = coeffs_29 + int(68);
    float * _S1315 = coeffs_29 + int(71);
    float3  _S1316 = v_colors_19 * make_float3 (float((0.282094806432724f * coeff_dc_29.x + 0.48860251903533936f * (_S1241 * *_S1242 + z_22 * *_S1243 - x_25 * *_S1244) + (pSH4_14 * *_S1246 + pSH5_14 * *_S1247 + pSH6_15 * *_S1248 + pSH7_14 * *_S1249 + pSH8_14 * *_S1250) + (pSH9_9 * *_S1251 + pSH10_9 * *_S1252 + pSH11_9 * *_S1253 + pSH12_10 * *_S1254 + pSH13_9 * *_S1255 + pSH14_9 * *_S1256 + pSH15_9 * *_S1257) + (pSH16_4 * *_S1259 + pSH17_4 * *_S1260 + pSH18_4 * *_S1261 + pSH19_4 * *_S1262 + pSH20_5 * *_S1263 + pSH21_4 * *_S1264 + pSH22_4 * *_S1265 + pSH23_4 * *_S1266 + pSH24_4 * *_S1267)) >= -0.5f), float((0.282094806432724f * coeff_dc_29.y + 0.48860251903533936f * (_S1241 * *_S1268 + z_22 * *_S1269 - x_25 * *_S1270) + (pSH4_14 * *_S1271 + pSH5_14 * *_S1272 + pSH6_15 * *_S1273 + pSH7_14 * *_S1274 + pSH8_14 * *_S1275) + (pSH9_9 * *_S1276 + pSH10_9 * *_S1277 + pSH11_9 * *_S1278 + pSH12_10 * *_S1279 + pSH13_9 * *_S1280 + pSH14_9 * *_S1281 + pSH15_9 * *_S1282) + (pSH16_4 * *_S1283 + pSH17_4 * *_S1284 + pSH18_4 * *_S1285 + pSH19_4 * *_S1286 + pSH20_5 * *_S1287 + pSH21_4 * *_S1288 + pSH22_4 * *_S1289 + pSH23_4 * *_S1290 + pSH24_4 * *_S1291)) >= -0.5f), float((0.282094806432724f * coeff_dc_29.z + 0.48860251903533936f * (_S1241 * *_S1292 + z_22 * *_S1293 - x_25 * *_S1294) + (pSH4_14 * *_S1295 + pSH5_14 * *_S1296 + pSH6_15 * *_S1297 + pSH7_14 * *_S1298 + pSH8_14 * *_S1299) + (pSH9_9 * *_S1300 + pSH10_9 * *_S1301 + pSH11_9 * *_S1302 + pSH12_10 * *_S1303 + pSH13_9 * *_S1304 + pSH14_9 * *_S1305 + pSH15_9 * *_S1306) + (pSH16_4 * *_S1307 + pSH17_4 * *_S1308 + pSH18_4 * *_S1309 + pSH19_4 * *_S1310 + pSH20_5 * *_S1311 + pSH21_4 * *_S1312 + pSH22_4 * *_S1313 + pSH23_4 * *_S1314 + pSH24_4 * *_S1315)) >= -0.5f));
    float3  v_viewdir_48 = {};
    float _S1317 = _S1316.x;
    *&(v_coeff_dc_19->x) = *&(v_coeff_dc_19->x) + 0.282094806432724f * _S1317;
    float _S1318 = -0.48860251903533936f * y_23;
    float temp_128 = _S1318 * _S1317;
    bool _S1319;
    if((F32_isfinite((temp_128))))
    {
        _S1319 = temp_128 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1320 = atomicAdd(v_coeffs_19 + int(0), temp_128);
    }
    float _S1321 = 0.48860251903533936f * z_22;
    float temp_129 = _S1321 * _S1317;
    if((F32_isfinite((temp_129))))
    {
        _S1319 = temp_129 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1322 = atomicAdd(v_coeffs_19 + int(3), temp_129);
    }
    float _S1323 = -0.48860251903533936f * x_25;
    float temp_130 = _S1323 * _S1317;
    if((F32_isfinite((temp_130))))
    {
        _S1319 = temp_130 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1324 = atomicAdd(v_coeffs_19 + int(6), temp_130);
    }
    float _S1325 = -0.48860251903533936f * *_S1244 * _S1317;
    float _S1326 = -0.48860251903533936f * *_S1242 * _S1317;
    float _S1327 = 0.48860251903533936f * *_S1243 * _S1317;
    float temp_131 = pSH4_14 * _S1317;
    if((F32_isfinite((temp_131))))
    {
        _S1319 = temp_131 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1328 = atomicAdd(v_coeffs_19 + int(9), temp_131);
    }
    float temp_132 = pSH5_14 * _S1317;
    if((F32_isfinite((temp_132))))
    {
        _S1319 = temp_132 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1329 = atomicAdd(v_coeffs_19 + int(12), temp_132);
    }
    float temp_133 = pSH6_15 * _S1317;
    if((F32_isfinite((temp_133))))
    {
        _S1319 = temp_133 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1330 = atomicAdd(v_coeffs_19 + int(15), temp_133);
    }
    float temp_134 = pSH7_14 * _S1317;
    if((F32_isfinite((temp_134))))
    {
        _S1319 = temp_134 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1331 = atomicAdd(v_coeffs_19 + int(18), temp_134);
    }
    float temp_135 = pSH8_14 * _S1317;
    if((F32_isfinite((temp_135))))
    {
        _S1319 = temp_135 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1332 = atomicAdd(v_coeffs_19 + int(21), temp_135);
    }
    float fC1_y_7 = -2.0f * y_23;
    float fS1_x_7 = 2.0f * y_23;
    float pSH6_z_7 = 1.89234936237335205f * z_22;
    float pSH7_z_5 = -1.09254848957061768f * x_25;
    float pSH5_z_5 = -1.09254848957061768f * y_23;
    float pSH8_x_11 = 0.54627424478530884f * _S1245;
    float pSH8_y_5 = 0.54627424478530884f * fC1_y_7;
    float pSH4_x_5 = 0.54627424478530884f * fS1_x_7;
    float v_x_18 = _S1325 + _S1317 * (pSH4_x_5 * *_S1246 + pSH8_x_11 * *_S1250 + fTmp0B_17 * *_S1249);
    float v_y_18 = _S1326 + _S1317 * (pSH8_x_11 * *_S1246 + pSH8_y_5 * *_S1250 + fTmp0B_17 * *_S1247);
    float v_z_18 = _S1327 + _S1317 * (pSH6_z_7 * *_S1248 + pSH7_z_5 * *_S1249 + pSH5_z_5 * *_S1247);
    float temp_136 = pSH9_9 * _S1317;
    if((F32_isfinite((temp_136))))
    {
        _S1319 = temp_136 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1333 = atomicAdd(v_coeffs_19 + int(24), temp_136);
    }
    float temp_137 = pSH10_9 * _S1317;
    if((F32_isfinite((temp_137))))
    {
        _S1319 = temp_137 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1334 = atomicAdd(v_coeffs_19 + int(27), temp_137);
    }
    float temp_138 = pSH11_9 * _S1317;
    if((F32_isfinite((temp_138))))
    {
        _S1319 = temp_138 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1335 = atomicAdd(v_coeffs_19 + int(30), temp_138);
    }
    float temp_139 = pSH12_10 * _S1317;
    if((F32_isfinite((temp_139))))
    {
        _S1319 = temp_139 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1336 = atomicAdd(v_coeffs_19 + int(33), temp_139);
    }
    float temp_140 = pSH13_9 * _S1317;
    if((F32_isfinite((temp_140))))
    {
        _S1319 = temp_140 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1337 = atomicAdd(v_coeffs_19 + int(36), temp_140);
    }
    float temp_141 = pSH14_9 * _S1317;
    if((F32_isfinite((temp_141))))
    {
        _S1319 = temp_141 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1338 = atomicAdd(v_coeffs_19 + int(39), temp_141);
    }
    float temp_142 = pSH15_9 * _S1317;
    if((F32_isfinite((temp_142))))
    {
        _S1319 = temp_142 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1339 = atomicAdd(v_coeffs_19 + int(42), temp_142);
    }
    float fTmp0C_z_7 = -4.57045793533325195f * z_22;
    float _S1340 = x_25 * _S1245;
    float fC2_x_3 = fC1_11 + _S1340 - y_23 * fS1_x_7;
    float _S1341 = y_23 * _S1245;
    float fC2_y_3 = x_25 * fC1_y_7 - fS1_11 - _S1341;
    float fS2_x_3 = fS1_11 + x_25 * fS1_x_7 + _S1341;
    float fS2_y_3 = _S1340 + fC1_11 + y_23 * fC1_y_7;
    float pSH12_z_5 = 5.59764480590820312f * z2_11 - 1.11952900886535645f;
    float pSH13_z_3 = fTmp0C_z_7 * x_25;
    float pSH11_z_3 = fTmp0C_z_7 * y_23;
    float pSH14_x_7 = fTmp1B_11 * _S1245;
    float pSH14_y_3 = fTmp1B_11 * fC1_y_7;
    float pSH14_z_3 = 1.44530570507049561f * fC1_11;
    float pSH10_x_3 = fTmp1B_11 * fS1_x_7;
    float pSH10_z_3 = 1.44530570507049561f * fS1_11;
    float pSH15_x_3 = -0.59004360437393188f * fC2_x_3;
    float pSH15_y_3 = -0.59004360437393188f * fC2_y_3;
    float pSH9_x_3 = -0.59004360437393188f * fS2_x_3;
    float pSH9_y_3 = -0.59004360437393188f * fS2_y_3;
    float v_x_19 = v_x_18 + _S1317 * (pSH9_x_3 * *_S1251 + pSH15_x_3 * *_S1257 + pSH10_x_3 * *_S1252 + pSH14_x_7 * *_S1256 + fTmp0C_11 * *_S1255);
    float v_y_19 = v_y_18 + _S1317 * (pSH9_y_3 * *_S1251 + pSH15_y_3 * *_S1257 + pSH14_x_7 * *_S1252 + pSH14_y_3 * *_S1256 + fTmp0C_11 * *_S1253);
    float v_z_19 = v_z_18 + _S1317 * (pSH12_z_5 * *_S1254 + pSH13_z_3 * *_S1255 + pSH11_z_3 * *_S1253 + pSH14_z_3 * *_S1256 + pSH10_z_3 * *_S1252);
    float pSH20_6 = _S1258 + -1.00623059272766113f * pSH6_15;
    float temp_143 = pSH16_4 * _S1317;
    if((F32_isfinite((temp_143))))
    {
        _S1319 = temp_143 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1342 = atomicAdd(v_coeffs_19 + int(45), temp_143);
    }
    float temp_144 = pSH17_4 * _S1317;
    if((F32_isfinite((temp_144))))
    {
        _S1319 = temp_144 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1343 = atomicAdd(v_coeffs_19 + int(48), temp_144);
    }
    float temp_145 = pSH18_4 * _S1317;
    if((F32_isfinite((temp_145))))
    {
        _S1319 = temp_145 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1344 = atomicAdd(v_coeffs_19 + int(51), temp_145);
    }
    float temp_146 = pSH19_4 * _S1317;
    if((F32_isfinite((temp_146))))
    {
        _S1319 = temp_146 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1345 = atomicAdd(v_coeffs_19 + int(54), temp_146);
    }
    float temp_147 = pSH20_6 * _S1317;
    if((F32_isfinite((temp_147))))
    {
        _S1319 = temp_147 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1346 = atomicAdd(v_coeffs_19 + int(57), temp_147);
    }
    float temp_148 = pSH21_4 * _S1317;
    if((F32_isfinite((temp_148))))
    {
        _S1319 = temp_148 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1347 = atomicAdd(v_coeffs_19 + int(60), temp_148);
    }
    float temp_149 = pSH22_4 * _S1317;
    if((F32_isfinite((temp_149))))
    {
        _S1319 = temp_149 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1348 = atomicAdd(v_coeffs_19 + int(63), temp_149);
    }
    float temp_150 = pSH23_4 * _S1317;
    if((F32_isfinite((temp_150))))
    {
        _S1319 = temp_150 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1349 = atomicAdd(v_coeffs_19 + int(66), temp_150);
    }
    float temp_151 = pSH24_4 * _S1317;
    if((F32_isfinite((temp_151))))
    {
        _S1319 = temp_151 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1350 = atomicAdd(v_coeffs_19 + int(69), temp_151);
    }
    float fTmp0D_z_3 = -14.04997730255126953f * z2_11 + 2.00713968276977539f;
    float fTmp1C_z_3 = 6.62322282791137695f * z_22;
    float pSH20_z_1 = 1.9843134880065918f * (pSH12_10 + z_22 * pSH12_z_5) + -1.00623059272766113f * pSH6_z_7;
    float pSH21_z_1 = fTmp0D_z_3 * x_25;
    float pSH19_z_1 = fTmp0D_z_3 * y_23;
    float pSH22_x_3 = fTmp1C_5 * _S1245;
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
    float3  v_dir_n_29 = make_float3 (v_x_19 + _S1317 * (pSH16_x_1 * *_S1259 + pSH24_x_1 * *_S1267 + pSH17_x_1 * *_S1260 + pSH23_x_1 * *_S1266 + pSH18_x_1 * *_S1261 + pSH22_x_3 * *_S1265 + fTmp0D_5 * *_S1264), v_y_19 + _S1317 * (pSH16_y_1 * *_S1259 + pSH24_y_1 * *_S1267 + pSH17_y_1 * *_S1260 + pSH23_y_1 * *_S1266 + pSH22_x_3 * *_S1261 + pSH22_y_1 * *_S1265 + fTmp0D_5 * *_S1262), v_z_19 + _S1317 * (pSH20_z_1 * *_S1263 + pSH21_z_1 * *_S1264 + pSH19_z_1 * *_S1262 + pSH22_z_1 * *_S1265 + pSH18_z_1 * *_S1261 + pSH23_z_1 * *_S1266 + pSH17_z_1 * *_S1260));
    float3  v_viewdir_49 = v_viewdir_48 + (v_dir_n_29 - make_float3 (dot_0(v_dir_n_29, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    float _S1351 = _S1316.y;
    *&(v_coeff_dc_19->y) = *&(v_coeff_dc_19->y) + 0.282094806432724f * _S1351;
    float temp_152 = _S1318 * _S1351;
    if((F32_isfinite((temp_152))))
    {
        _S1319 = temp_152 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1352 = atomicAdd(v_coeffs_19 + int(1), temp_152);
    }
    float temp_153 = _S1321 * _S1351;
    if((F32_isfinite((temp_153))))
    {
        _S1319 = temp_153 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1353 = atomicAdd(v_coeffs_19 + int(4), temp_153);
    }
    float temp_154 = _S1323 * _S1351;
    if((F32_isfinite((temp_154))))
    {
        _S1319 = temp_154 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1354 = atomicAdd(v_coeffs_19 + int(7), temp_154);
    }
    float _S1355 = -0.48860251903533936f * *_S1270 * _S1351;
    float _S1356 = -0.48860251903533936f * *_S1268 * _S1351;
    float _S1357 = 0.48860251903533936f * *_S1269 * _S1351;
    float temp_155 = pSH4_14 * _S1351;
    if((F32_isfinite((temp_155))))
    {
        _S1319 = temp_155 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1358 = atomicAdd(v_coeffs_19 + int(10), temp_155);
    }
    float temp_156 = pSH5_14 * _S1351;
    if((F32_isfinite((temp_156))))
    {
        _S1319 = temp_156 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1359 = atomicAdd(v_coeffs_19 + int(13), temp_156);
    }
    float temp_157 = pSH6_15 * _S1351;
    if((F32_isfinite((temp_157))))
    {
        _S1319 = temp_157 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1360 = atomicAdd(v_coeffs_19 + int(16), temp_157);
    }
    float temp_158 = pSH7_14 * _S1351;
    if((F32_isfinite((temp_158))))
    {
        _S1319 = temp_158 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1361 = atomicAdd(v_coeffs_19 + int(19), temp_158);
    }
    float temp_159 = pSH8_14 * _S1351;
    if((F32_isfinite((temp_159))))
    {
        _S1319 = temp_159 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1362 = atomicAdd(v_coeffs_19 + int(22), temp_159);
    }
    float v_x_20 = _S1355 + _S1351 * (pSH4_x_5 * *_S1271 + pSH8_x_11 * *_S1275 + fTmp0B_17 * *_S1274);
    float v_y_20 = _S1356 + _S1351 * (pSH8_x_11 * *_S1271 + pSH8_y_5 * *_S1275 + fTmp0B_17 * *_S1272);
    float v_z_20 = _S1357 + _S1351 * (pSH6_z_7 * *_S1273 + pSH7_z_5 * *_S1274 + pSH5_z_5 * *_S1272);
    float temp_160 = pSH9_9 * _S1351;
    if((F32_isfinite((temp_160))))
    {
        _S1319 = temp_160 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1363 = atomicAdd(v_coeffs_19 + int(25), temp_160);
    }
    float temp_161 = pSH10_9 * _S1351;
    if((F32_isfinite((temp_161))))
    {
        _S1319 = temp_161 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1364 = atomicAdd(v_coeffs_19 + int(28), temp_161);
    }
    float temp_162 = pSH11_9 * _S1351;
    if((F32_isfinite((temp_162))))
    {
        _S1319 = temp_162 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1365 = atomicAdd(v_coeffs_19 + int(31), temp_162);
    }
    float temp_163 = pSH12_10 * _S1351;
    if((F32_isfinite((temp_163))))
    {
        _S1319 = temp_163 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1366 = atomicAdd(v_coeffs_19 + int(34), temp_163);
    }
    float temp_164 = pSH13_9 * _S1351;
    if((F32_isfinite((temp_164))))
    {
        _S1319 = temp_164 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1367 = atomicAdd(v_coeffs_19 + int(37), temp_164);
    }
    float temp_165 = pSH14_9 * _S1351;
    if((F32_isfinite((temp_165))))
    {
        _S1319 = temp_165 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1368 = atomicAdd(v_coeffs_19 + int(40), temp_165);
    }
    float temp_166 = pSH15_9 * _S1351;
    if((F32_isfinite((temp_166))))
    {
        _S1319 = temp_166 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1369 = atomicAdd(v_coeffs_19 + int(43), temp_166);
    }
    float v_x_21 = v_x_20 + _S1351 * (pSH9_x_3 * *_S1276 + pSH15_x_3 * *_S1282 + pSH10_x_3 * *_S1277 + pSH14_x_7 * *_S1281 + fTmp0C_11 * *_S1280);
    float v_y_21 = v_y_20 + _S1351 * (pSH9_y_3 * *_S1276 + pSH15_y_3 * *_S1282 + pSH14_x_7 * *_S1277 + pSH14_y_3 * *_S1281 + fTmp0C_11 * *_S1278);
    float v_z_21 = v_z_20 + _S1351 * (pSH12_z_5 * *_S1279 + pSH13_z_3 * *_S1280 + pSH11_z_3 * *_S1278 + pSH14_z_3 * *_S1281 + pSH10_z_3 * *_S1277);
    float temp_167 = pSH16_4 * _S1351;
    if((F32_isfinite((temp_167))))
    {
        _S1319 = temp_167 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1370 = atomicAdd(v_coeffs_19 + int(46), temp_167);
    }
    float temp_168 = pSH17_4 * _S1351;
    if((F32_isfinite((temp_168))))
    {
        _S1319 = temp_168 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1371 = atomicAdd(v_coeffs_19 + int(49), temp_168);
    }
    float temp_169 = pSH18_4 * _S1351;
    if((F32_isfinite((temp_169))))
    {
        _S1319 = temp_169 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1372 = atomicAdd(v_coeffs_19 + int(52), temp_169);
    }
    float temp_170 = pSH19_4 * _S1351;
    if((F32_isfinite((temp_170))))
    {
        _S1319 = temp_170 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1373 = atomicAdd(v_coeffs_19 + int(55), temp_170);
    }
    float temp_171 = pSH20_6 * _S1351;
    if((F32_isfinite((temp_171))))
    {
        _S1319 = temp_171 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1374 = atomicAdd(v_coeffs_19 + int(58), temp_171);
    }
    float temp_172 = pSH21_4 * _S1351;
    if((F32_isfinite((temp_172))))
    {
        _S1319 = temp_172 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1375 = atomicAdd(v_coeffs_19 + int(61), temp_172);
    }
    float temp_173 = pSH22_4 * _S1351;
    if((F32_isfinite((temp_173))))
    {
        _S1319 = temp_173 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1376 = atomicAdd(v_coeffs_19 + int(64), temp_173);
    }
    float temp_174 = pSH23_4 * _S1351;
    if((F32_isfinite((temp_174))))
    {
        _S1319 = temp_174 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1377 = atomicAdd(v_coeffs_19 + int(67), temp_174);
    }
    float temp_175 = pSH24_4 * _S1351;
    if((F32_isfinite((temp_175))))
    {
        _S1319 = temp_175 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1378 = atomicAdd(v_coeffs_19 + int(70), temp_175);
    }
    float3  v_dir_n_30 = make_float3 (v_x_21 + _S1351 * (pSH16_x_1 * *_S1283 + pSH24_x_1 * *_S1291 + pSH17_x_1 * *_S1284 + pSH23_x_1 * *_S1290 + pSH18_x_1 * *_S1285 + pSH22_x_3 * *_S1289 + fTmp0D_5 * *_S1288), v_y_21 + _S1351 * (pSH16_y_1 * *_S1283 + pSH24_y_1 * *_S1291 + pSH17_y_1 * *_S1284 + pSH23_y_1 * *_S1290 + pSH22_x_3 * *_S1285 + pSH22_y_1 * *_S1289 + fTmp0D_5 * *_S1286), v_z_21 + _S1351 * (pSH20_z_1 * *_S1287 + pSH21_z_1 * *_S1288 + pSH19_z_1 * *_S1286 + pSH22_z_1 * *_S1289 + pSH18_z_1 * *_S1285 + pSH23_z_1 * *_S1290 + pSH17_z_1 * *_S1284));
    float3  v_viewdir_50 = v_viewdir_49 + (v_dir_n_30 - make_float3 (dot_0(v_dir_n_30, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    float _S1379 = _S1316.z;
    *&(v_coeff_dc_19->z) = *&(v_coeff_dc_19->z) + 0.282094806432724f * _S1379;
    float temp_176 = _S1318 * _S1379;
    if((F32_isfinite((temp_176))))
    {
        _S1319 = temp_176 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1380 = atomicAdd(v_coeffs_19 + int(2), temp_176);
    }
    float temp_177 = _S1321 * _S1379;
    if((F32_isfinite((temp_177))))
    {
        _S1319 = temp_177 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1381 = atomicAdd(v_coeffs_19 + int(5), temp_177);
    }
    float temp_178 = _S1323 * _S1379;
    if((F32_isfinite((temp_178))))
    {
        _S1319 = temp_178 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1382 = atomicAdd(v_coeffs_19 + int(8), temp_178);
    }
    float _S1383 = -0.48860251903533936f * *_S1294 * _S1379;
    float _S1384 = -0.48860251903533936f * *_S1292 * _S1379;
    float _S1385 = 0.48860251903533936f * *_S1293 * _S1379;
    float temp_179 = pSH4_14 * _S1379;
    if((F32_isfinite((temp_179))))
    {
        _S1319 = temp_179 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1386 = atomicAdd(v_coeffs_19 + int(11), temp_179);
    }
    float temp_180 = pSH5_14 * _S1379;
    if((F32_isfinite((temp_180))))
    {
        _S1319 = temp_180 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1387 = atomicAdd(v_coeffs_19 + int(14), temp_180);
    }
    float temp_181 = pSH6_15 * _S1379;
    if((F32_isfinite((temp_181))))
    {
        _S1319 = temp_181 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1388 = atomicAdd(v_coeffs_19 + int(17), temp_181);
    }
    float temp_182 = pSH7_14 * _S1379;
    if((F32_isfinite((temp_182))))
    {
        _S1319 = temp_182 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1389 = atomicAdd(v_coeffs_19 + int(20), temp_182);
    }
    float temp_183 = pSH8_14 * _S1379;
    if((F32_isfinite((temp_183))))
    {
        _S1319 = temp_183 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1390 = atomicAdd(v_coeffs_19 + int(23), temp_183);
    }
    float v_x_22 = _S1383 + _S1379 * (pSH4_x_5 * *_S1295 + pSH8_x_11 * *_S1299 + fTmp0B_17 * *_S1298);
    float v_y_22 = _S1384 + _S1379 * (pSH8_x_11 * *_S1295 + pSH8_y_5 * *_S1299 + fTmp0B_17 * *_S1296);
    float v_z_22 = _S1385 + _S1379 * (pSH6_z_7 * *_S1297 + pSH7_z_5 * *_S1298 + pSH5_z_5 * *_S1296);
    float temp_184 = pSH9_9 * _S1379;
    if((F32_isfinite((temp_184))))
    {
        _S1319 = temp_184 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1391 = atomicAdd(v_coeffs_19 + int(26), temp_184);
    }
    float temp_185 = pSH10_9 * _S1379;
    if((F32_isfinite((temp_185))))
    {
        _S1319 = temp_185 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1392 = atomicAdd(v_coeffs_19 + int(29), temp_185);
    }
    float temp_186 = pSH11_9 * _S1379;
    if((F32_isfinite((temp_186))))
    {
        _S1319 = temp_186 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1393 = atomicAdd(v_coeffs_19 + int(32), temp_186);
    }
    float temp_187 = pSH12_10 * _S1379;
    if((F32_isfinite((temp_187))))
    {
        _S1319 = temp_187 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1394 = atomicAdd(v_coeffs_19 + int(35), temp_187);
    }
    float temp_188 = pSH13_9 * _S1379;
    if((F32_isfinite((temp_188))))
    {
        _S1319 = temp_188 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1395 = atomicAdd(v_coeffs_19 + int(38), temp_188);
    }
    float temp_189 = pSH14_9 * _S1379;
    if((F32_isfinite((temp_189))))
    {
        _S1319 = temp_189 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1396 = atomicAdd(v_coeffs_19 + int(41), temp_189);
    }
    float temp_190 = pSH15_9 * _S1379;
    if((F32_isfinite((temp_190))))
    {
        _S1319 = temp_190 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1397 = atomicAdd(v_coeffs_19 + int(44), temp_190);
    }
    float v_x_23 = v_x_22 + _S1379 * (pSH9_x_3 * *_S1300 + pSH15_x_3 * *_S1306 + pSH10_x_3 * *_S1301 + pSH14_x_7 * *_S1305 + fTmp0C_11 * *_S1304);
    float v_y_23 = v_y_22 + _S1379 * (pSH9_y_3 * *_S1300 + pSH15_y_3 * *_S1306 + pSH14_x_7 * *_S1301 + pSH14_y_3 * *_S1305 + fTmp0C_11 * *_S1302);
    float v_z_23 = v_z_22 + _S1379 * (pSH12_z_5 * *_S1303 + pSH13_z_3 * *_S1304 + pSH11_z_3 * *_S1302 + pSH14_z_3 * *_S1305 + pSH10_z_3 * *_S1301);
    float temp_191 = pSH16_4 * _S1379;
    if((F32_isfinite((temp_191))))
    {
        _S1319 = temp_191 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1398 = atomicAdd(v_coeffs_19 + int(47), temp_191);
    }
    float temp_192 = pSH17_4 * _S1379;
    if((F32_isfinite((temp_192))))
    {
        _S1319 = temp_192 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1399 = atomicAdd(v_coeffs_19 + int(50), temp_192);
    }
    float temp_193 = pSH18_4 * _S1379;
    if((F32_isfinite((temp_193))))
    {
        _S1319 = temp_193 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1400 = atomicAdd(v_coeffs_19 + int(53), temp_193);
    }
    float temp_194 = pSH19_4 * _S1379;
    if((F32_isfinite((temp_194))))
    {
        _S1319 = temp_194 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1401 = atomicAdd(v_coeffs_19 + int(56), temp_194);
    }
    float temp_195 = pSH20_6 * _S1379;
    if((F32_isfinite((temp_195))))
    {
        _S1319 = temp_195 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1402 = atomicAdd(v_coeffs_19 + int(59), temp_195);
    }
    float temp_196 = pSH21_4 * _S1379;
    if((F32_isfinite((temp_196))))
    {
        _S1319 = temp_196 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1403 = atomicAdd(v_coeffs_19 + int(62), temp_196);
    }
    float temp_197 = pSH22_4 * _S1379;
    if((F32_isfinite((temp_197))))
    {
        _S1319 = temp_197 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1404 = atomicAdd(v_coeffs_19 + int(65), temp_197);
    }
    float temp_198 = pSH23_4 * _S1379;
    if((F32_isfinite((temp_198))))
    {
        _S1319 = temp_198 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1405 = atomicAdd(v_coeffs_19 + int(68), temp_198);
    }
    float temp_199 = pSH24_4 * _S1379;
    if((F32_isfinite((temp_199))))
    {
        _S1319 = temp_199 != 0.0f;
    }
    else
    {
        _S1319 = false;
    }
    if(_S1319)
    {
        float _S1406 = atomicAdd(v_coeffs_19 + int(71), temp_199);
    }
    float3  v_dir_n_31 = make_float3 (v_x_23 + _S1379 * (pSH16_x_1 * *_S1307 + pSH24_x_1 * *_S1315 + pSH17_x_1 * *_S1308 + pSH23_x_1 * *_S1314 + pSH18_x_1 * *_S1309 + pSH22_x_3 * *_S1313 + fTmp0D_5 * *_S1312), v_y_23 + _S1379 * (pSH16_y_1 * *_S1307 + pSH24_y_1 * *_S1315 + pSH17_y_1 * *_S1308 + pSH23_y_1 * *_S1314 + pSH22_x_3 * *_S1309 + pSH22_y_1 * *_S1313 + fTmp0D_5 * *_S1310), v_z_23 + _S1379 * (pSH20_z_1 * *_S1311 + pSH21_z_1 * *_S1312 + pSH19_z_1 * *_S1310 + pSH22_z_1 * *_S1313 + pSH18_z_1 * *_S1309 + pSH23_z_1 * *_S1314 + pSH17_z_1 * *_S1308));
    float3  v_viewdir_51 = v_viewdir_50 + (v_dir_n_31 - make_float3 (dot_0(v_dir_n_31, dir_n_15)) * dir_n_15) * make_float3 (inv_norm_23);
    Matrix<float, 3, 3>  _S1407 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1408;
    (&_S1408)->primal_0 = _S1236;
    (&_S1408)->differential_0 = _S1407;
    float3  _S1409 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1410;
    (&_S1410)->primal_0 = t_29;
    (&_S1410)->differential_0 = _S1409;
    s_bwd_prop_mul_0(&_S1408, &_S1410, v_viewdir_51);
    Matrix<float, 3, 3>  _S1411 = transpose_0(_S1408.differential_0);
    *v_mean_19 = *v_mean_19 + v_viewdir_51;
    *v_R_19 = *v_R_19 + _S1411;
    *v_t_19 = *v_t_19 + _S1410.differential_0;
    return;
}

inline __device__ float3  sh0_to_color_dir(float3  dir_0, float3  coeff_dc_30, float3  * coeffs_30)
{
    return make_float3 (0.282094806432724f) * coeff_dc_30;
}

inline __device__ void sh0_to_color_dir_vjp_inplace(float3  dir_1, float3  coeff_dc_31, float3  * coeffs_31, float3  v_colors_20, float3  * v_coeff_dc_20, float3  * v_coeffs_20, float3  * v_dir_0)
{
    *v_coeff_dc_20 = *v_coeff_dc_20 + make_float3 (0.282094806432724f) * v_colors_20;
    return;
}

inline __device__ void sh0_to_color_dir_vjp_atomic(float3  dir_2, float3  coeff_dc_32, float3  * coeffs_32, float3  v_colors_21, float3  * v_coeff_dc_21, float3  * v_coeffs_21, float3  * v_dir_1)
{
    *v_coeff_dc_21 = *v_coeff_dc_21 + make_float3 (0.282094806432724f) * v_colors_21;
    return;
}

inline __device__ void sh0_to_color_dir_vjp_block_atomic(float3  dir_3, float3  coeff_dc_33, float3  * coeffs_33, float3  v_colors_22, float3  * v_coeff_dc_22, float3  * v_coeffs_22, float3  * v_dir_2, uint thread_id_0)
{
    *v_coeff_dc_22 = *v_coeff_dc_22 + make_float3 (0.282094806432724f) * v_colors_22;
    return;
}

inline __device__ float3  sh1_to_color_dir(float3  dir_4, float3  coeff_dc_34, float3  * coeffs_34)
{
    float _S1412 = dir_4.x;
    float _S1413 = dir_4.y;
    float _S1414 = dir_4.z;
    float inv_norm_24 = (F32_rsqrt((_S1412 * _S1412 + _S1413 * _S1413 + _S1414 * _S1414)));
    return make_float3 (0.282094806432724f) * coeff_dc_34 + make_float3 (0.48860251903533936f) * (make_float3 (- (_S1413 * inv_norm_24)) * *(coeffs_34 + int(0)) + make_float3 (_S1414 * inv_norm_24) * *(coeffs_34 + int(1)) - make_float3 (_S1412 * inv_norm_24) * *(coeffs_34 + int(2)));
}

inline __device__ void sh1_to_color_dir_vjp_inplace(float3  dir_5, float3  coeff_dc_35, float3  * coeffs_35, float3  v_colors_23, float3  * v_coeff_dc_23, float3  * v_coeffs_23, float3  * v_dir_3)
{
    *v_coeff_dc_23 = *v_coeff_dc_23 + make_float3 (0.282094806432724f) * v_colors_23;
    float _S1415 = dir_5.x;
    float _S1416 = dir_5.y;
    float _S1417 = dir_5.z;
    float inorm_0 = (F32_rsqrt((_S1415 * _S1415 + _S1416 * _S1416 + _S1417 * _S1417)));
    float x_26 = _S1415 * inorm_0;
    float y_24 = _S1416 * inorm_0;
    float z_23 = _S1417 * inorm_0;
    float3  * _S1418 = v_coeffs_23 + int(0);
    *_S1418 = *_S1418 + make_float3 (-0.48860251903533936f * y_24) * v_colors_23;
    float3  * _S1419 = v_coeffs_23 + int(1);
    *_S1419 = *_S1419 + make_float3 (0.48860251903533936f * z_23) * v_colors_23;
    float3  * _S1420 = v_coeffs_23 + int(2);
    *_S1420 = *_S1420 + make_float3 (-0.48860251903533936f * x_26) * v_colors_23;
    float3  dir_n_16 = make_float3 (x_26, y_24, z_23);
    float3  v_dir_n_32 = make_float3 (-0.48860251903533936f * dot_0(*(coeffs_35 + int(2)), v_colors_23), -0.48860251903533936f * dot_0(*(coeffs_35 + int(0)), v_colors_23), 0.48860251903533936f * dot_0(*(coeffs_35 + int(1)), v_colors_23));
    *v_dir_3 = *v_dir_3 + (v_dir_n_32 - make_float3 (dot_0(v_dir_n_32, dir_n_16)) * dir_n_16) * make_float3 (inorm_0);
    return;
}

inline __device__ void sh1_to_color_dir_vjp_atomic(float3  dir_6, float3  coeff_dc_36, float3  * coeffs_36, float3  v_colors_24, float3  * v_coeff_dc_24, float3  * v_coeffs_24, float3  * v_dir_4)
{
    *v_coeff_dc_24 = *v_coeff_dc_24 + make_float3 (0.282094806432724f) * v_colors_24;
    float _S1421 = dir_6.x;
    float _S1422 = dir_6.y;
    float _S1423 = dir_6.z;
    float inorm_1 = (F32_rsqrt((_S1421 * _S1421 + _S1422 * _S1422 + _S1423 * _S1423)));
    float x_27 = _S1421 * inorm_1;
    float y_25 = _S1422 * inorm_1;
    float z_24 = _S1423 * inorm_1;
    float3  temp_200 = make_float3 (-0.48860251903533936f * y_25) * v_colors_24;
    float _S1424 = dot_0(temp_200, temp_200);
    bool _S1425;
    if((F32_isfinite((_S1424))))
    {
        _S1425 = _S1424 != 0.0f;
    }
    else
    {
        _S1425 = false;
    }
    if(_S1425)
    {
        float3  * _S1426 = v_coeffs_24 + int(0);
        float _S1427 = atomicAdd(&(_S1426->x), temp_200.x);
        float _S1428 = atomicAdd(&(_S1426->y), temp_200.y);
        float _S1429 = atomicAdd(&(_S1426->z), temp_200.z);
    }
    float3  temp_201 = make_float3 (0.48860251903533936f * z_24) * v_colors_24;
    float _S1430 = dot_0(temp_201, temp_201);
    if((F32_isfinite((_S1430))))
    {
        _S1425 = _S1430 != 0.0f;
    }
    else
    {
        _S1425 = false;
    }
    if(_S1425)
    {
        float3  * _S1431 = v_coeffs_24 + int(1);
        float _S1432 = atomicAdd(&(_S1431->x), temp_201.x);
        float _S1433 = atomicAdd(&(_S1431->y), temp_201.y);
        float _S1434 = atomicAdd(&(_S1431->z), temp_201.z);
    }
    float3  temp_202 = make_float3 (-0.48860251903533936f * x_27) * v_colors_24;
    float _S1435 = dot_0(temp_202, temp_202);
    if((F32_isfinite((_S1435))))
    {
        _S1425 = _S1435 != 0.0f;
    }
    else
    {
        _S1425 = false;
    }
    if(_S1425)
    {
        float3  * _S1436 = v_coeffs_24 + int(2);
        float _S1437 = atomicAdd(&(_S1436->x), temp_202.x);
        float _S1438 = atomicAdd(&(_S1436->y), temp_202.y);
        float _S1439 = atomicAdd(&(_S1436->z), temp_202.z);
    }
    float3  dir_n_17 = make_float3 (x_27, y_25, z_24);
    float3  v_dir_n_33 = make_float3 (-0.48860251903533936f * dot_0(*(coeffs_36 + int(2)), v_colors_24), -0.48860251903533936f * dot_0(*(coeffs_36 + int(0)), v_colors_24), 0.48860251903533936f * dot_0(*(coeffs_36 + int(1)), v_colors_24));
    *v_dir_4 = *v_dir_4 + (v_dir_n_33 - make_float3 (dot_0(v_dir_n_33, dir_n_17)) * dir_n_17) * make_float3 (inorm_1);
    return;
}

inline __device__ void sh1_to_color_dir_vjp_block_atomic(float3  dir_7, float3  coeff_dc_37, float3  * coeffs_37, float3  v_colors_25, float3  * v_coeff_dc_25, float3  * v_coeffs_25, float3  * v_dir_5, uint thread_id_1, uint _S1440)
{
    *v_coeff_dc_25 = *v_coeff_dc_25 + make_float3 (0.282094806432724f) * v_colors_25;
    float _S1441 = dir_7.x;
    float _S1442 = dir_7.y;
    float _S1443 = dir_7.z;
    float inorm_2 = (F32_rsqrt((_S1441 * _S1441 + _S1442 * _S1442 + _S1443 * _S1443)));
    float x_28 = _S1441 * inorm_2;
    float y_26 = _S1442 * inorm_2;
    float z_25 = _S1443 * inorm_2;
    float3  _S1444 = make_float3 (-0.48860251903533936f * y_26) * v_colors_25;
    float3  _S1445 = _S1444;
    bool _S1446 = (F32_isfinite((_S1444.x)));
    uint _S1447 = __ballot_sync(_S1440, _S1446);
    float v_0;
    uint _S1448;
    if(_S1446)
    {
        float _S1449 = _S1445.x;
        uint _S1450 = __ballot_sync(_S1440, true);
        v_0 = _S1449;
        _S1448 = _S1450;
    }
    else
    {
        uint _S1451 = __ballot_sync(_S1440, true);
        v_0 = 0.0f;
        _S1448 = _S1451;
    }
    *&((&_S1445)->x) = v_0;
    bool _S1452 = (F32_isfinite((_S1445.y)));
    uint _S1453 = __ballot_sync(_S1448, _S1452);
    if(_S1452)
    {
        float _S1454 = _S1445.y;
        uint _S1455 = __ballot_sync(_S1448, true);
        v_0 = _S1454;
        _S1448 = _S1455;
    }
    else
    {
        uint _S1456 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1456;
    }
    *&((&_S1445)->y) = v_0;
    bool _S1457 = (F32_isfinite((_S1445.z)));
    uint _S1458 = __ballot_sync(_S1448, _S1457);
    if(_S1457)
    {
        float _S1459 = _S1445.z;
        uint _S1460 = __ballot_sync(_S1448, true);
        v_0 = _S1459;
        _S1448 = _S1460;
    }
    else
    {
        uint _S1461 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1461;
    }
    *&((&_S1445)->z) = v_0;
    float _S1462 = WaveActiveSum_0(_S1445.x, _S1448);
    *&((&_S1445)->x) = _S1462;
    float _S1463 = WaveActiveSum_0(_S1445.y, _S1448);
    *&((&_S1445)->y) = _S1463;
    float _S1464 = WaveActiveSum_0(_S1445.z, _S1448);
    *&((&_S1445)->z) = _S1464;
    uint warp_id_0 = thread_id_1 / 32U;
    uint lane_id_0 = thread_id_1 % 32U;
    bool _S1465 = lane_id_0 == 0U;
    uint _S1466 = __ballot_sync(_S1448, _S1465);
    if(_S1465)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_0] = _S1445.x;
        (*&_sh_block_reduce_shared_0)[warp_id_0 + 16U] = _S1445.y;
        (*&_sh_block_reduce_shared_0)[warp_id_0 + 32U] = _S1445.z;
        uint _S1467 = __ballot_sync(_S1448, true);
        _S1448 = _S1467;
    }
    else
    {
        uint _S1468 = __ballot_sync(_S1448, true);
        _S1448 = _S1468;
    }
    __syncthreads();
    bool _S1469 = warp_id_0 < 3U;
    uint _S1470 = __ballot_sync(_S1448, _S1469);
    uint _S1471;
    bool _S1472;
    if(_S1469)
    {
        bool _S1473 = lane_id_0 < 16U;
        uint _S1474 = __ballot_sync(_S1470, _S1473);
        if(_S1473)
        {
            float _S1475 = (*&_sh_block_reduce_shared_0)[lane_id_0 + warp_id_0 * 16U];
            uint _S1476 = __ballot_sync(_S1470, true);
            v_0 = _S1475;
            _S1471 = _S1476;
        }
        else
        {
            uint _S1477 = __ballot_sync(_S1470, true);
            v_0 = 0.0f;
            _S1471 = _S1477;
        }
        float _S1478 = WaveActiveSum_0(v_0, _S1471);
        uint _S1479 = __ballot_sync(_S1471, _S1465);
        if(_S1465)
        {
            bool _S1480 = _S1478 != 0.0f;
            uint _S1481 = __ballot_sync(_S1471, true);
            _S1472 = _S1480;
            _S1471 = _S1481;
        }
        else
        {
            uint _S1482 = __ballot_sync(_S1471, true);
            _S1472 = false;
            _S1471 = _S1482;
        }
        uint _S1483 = __ballot_sync(_S1471, _S1472);
        if(_S1472)
        {
            bool _S1484 = warp_id_0 == 0U;
            uint _S1485 = __ballot_sync(_S1483, _S1484);
            if(_S1484)
            {
                float _S1486 = atomicAdd(&((v_coeffs_25 + 0U)->x), _S1478);
                uint _S1487 = __ballot_sync(_S1483, true);
            }
            else
            {
                uint _S1488 = _S1483 & (~_S1485);
                bool _S1489 = warp_id_0 == 1U;
                uint _S1490 = __ballot_sync(_S1488, _S1489);
                if(_S1489)
                {
                    float _S1491 = atomicAdd(&((v_coeffs_25 + 0U)->y), _S1478);
                    uint _S1492 = __ballot_sync(_S1488, true);
                }
                else
                {
                    float _S1493 = atomicAdd(&((v_coeffs_25 + 0U)->z), _S1478);
                    uint _S1494 = __ballot_sync(_S1488, true);
                }
                uint _S1495 = __ballot_sync(_S1483, true);
            }
            uint _S1496 = __ballot_sync(_S1471, true);
        }
        else
        {
            uint _S1497 = __ballot_sync(_S1471, true);
        }
        uint _S1498 = __ballot_sync(_S1448, true);
        _S1448 = _S1498;
    }
    else
    {
        uint _S1499 = __ballot_sync(_S1448, true);
        _S1448 = _S1499;
    }
    float3  _S1500 = make_float3 (0.48860251903533936f * z_25) * v_colors_25;
    float3  _S1501 = _S1500;
    bool _S1502 = (F32_isfinite((_S1500.x)));
    uint _S1503 = __ballot_sync(_S1448, _S1502);
    if(_S1502)
    {
        float _S1504 = _S1501.x;
        uint _S1505 = __ballot_sync(_S1448, true);
        v_0 = _S1504;
        _S1448 = _S1505;
    }
    else
    {
        uint _S1506 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1506;
    }
    *&((&_S1501)->x) = v_0;
    bool _S1507 = (F32_isfinite((_S1501.y)));
    uint _S1508 = __ballot_sync(_S1448, _S1507);
    if(_S1507)
    {
        float _S1509 = _S1501.y;
        uint _S1510 = __ballot_sync(_S1448, true);
        v_0 = _S1509;
        _S1448 = _S1510;
    }
    else
    {
        uint _S1511 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1511;
    }
    *&((&_S1501)->y) = v_0;
    bool _S1512 = (F32_isfinite((_S1501.z)));
    uint _S1513 = __ballot_sync(_S1448, _S1512);
    if(_S1512)
    {
        float _S1514 = _S1501.z;
        uint _S1515 = __ballot_sync(_S1448, true);
        v_0 = _S1514;
        _S1448 = _S1515;
    }
    else
    {
        uint _S1516 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1516;
    }
    *&((&_S1501)->z) = v_0;
    float _S1517 = WaveActiveSum_0(_S1501.x, _S1448);
    *&((&_S1501)->x) = _S1517;
    float _S1518 = WaveActiveSum_0(_S1501.y, _S1448);
    *&((&_S1501)->y) = _S1518;
    float _S1519 = WaveActiveSum_0(_S1501.z, _S1448);
    *&((&_S1501)->z) = _S1519;
    uint warp_id_1 = thread_id_1 / 32U;
    uint _S1520 = __ballot_sync(_S1448, _S1465);
    if(_S1465)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_1] = _S1501.x;
        (*&_sh_block_reduce_shared_0)[warp_id_1 + 16U] = _S1501.y;
        (*&_sh_block_reduce_shared_0)[warp_id_1 + 32U] = _S1501.z;
        uint _S1521 = __ballot_sync(_S1448, true);
        _S1448 = _S1521;
    }
    else
    {
        uint _S1522 = __ballot_sync(_S1448, true);
        _S1448 = _S1522;
    }
    __syncthreads();
    bool _S1523 = warp_id_1 < 3U;
    uint _S1524 = __ballot_sync(_S1448, _S1523);
    if(_S1523)
    {
        bool _S1525 = lane_id_0 < 16U;
        uint _S1526 = __ballot_sync(_S1524, _S1525);
        if(_S1525)
        {
            float _S1527 = (*&_sh_block_reduce_shared_0)[lane_id_0 + warp_id_1 * 16U];
            uint _S1528 = __ballot_sync(_S1524, true);
            v_0 = _S1527;
            _S1471 = _S1528;
        }
        else
        {
            uint _S1529 = __ballot_sync(_S1524, true);
            v_0 = 0.0f;
            _S1471 = _S1529;
        }
        float _S1530 = WaveActiveSum_0(v_0, _S1471);
        uint _S1531 = __ballot_sync(_S1471, _S1465);
        if(_S1465)
        {
            bool _S1532 = _S1530 != 0.0f;
            uint _S1533 = __ballot_sync(_S1471, true);
            _S1472 = _S1532;
            _S1471 = _S1533;
        }
        else
        {
            uint _S1534 = __ballot_sync(_S1471, true);
            _S1472 = false;
            _S1471 = _S1534;
        }
        uint _S1535 = __ballot_sync(_S1471, _S1472);
        if(_S1472)
        {
            bool _S1536 = warp_id_1 == 0U;
            uint _S1537 = __ballot_sync(_S1535, _S1536);
            if(_S1536)
            {
                float _S1538 = atomicAdd(&((v_coeffs_25 + 1U)->x), _S1530);
                uint _S1539 = __ballot_sync(_S1535, true);
            }
            else
            {
                uint _S1540 = _S1535 & (~_S1537);
                bool _S1541 = warp_id_1 == 1U;
                uint _S1542 = __ballot_sync(_S1540, _S1541);
                if(_S1541)
                {
                    float _S1543 = atomicAdd(&((v_coeffs_25 + 1U)->y), _S1530);
                    uint _S1544 = __ballot_sync(_S1540, true);
                }
                else
                {
                    float _S1545 = atomicAdd(&((v_coeffs_25 + 1U)->z), _S1530);
                    uint _S1546 = __ballot_sync(_S1540, true);
                }
                uint _S1547 = __ballot_sync(_S1535, true);
            }
            uint _S1548 = __ballot_sync(_S1471, true);
        }
        else
        {
            uint _S1549 = __ballot_sync(_S1471, true);
        }
        uint _S1550 = __ballot_sync(_S1448, true);
        _S1448 = _S1550;
    }
    else
    {
        uint _S1551 = __ballot_sync(_S1448, true);
        _S1448 = _S1551;
    }
    float3  _S1552 = make_float3 (-0.48860251903533936f * x_28) * v_colors_25;
    float3  _S1553 = _S1552;
    bool _S1554 = (F32_isfinite((_S1552.x)));
    uint _S1555 = __ballot_sync(_S1448, _S1554);
    if(_S1554)
    {
        float _S1556 = _S1553.x;
        uint _S1557 = __ballot_sync(_S1448, true);
        v_0 = _S1556;
        _S1448 = _S1557;
    }
    else
    {
        uint _S1558 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1558;
    }
    *&((&_S1553)->x) = v_0;
    bool _S1559 = (F32_isfinite((_S1553.y)));
    uint _S1560 = __ballot_sync(_S1448, _S1559);
    if(_S1559)
    {
        float _S1561 = _S1553.y;
        uint _S1562 = __ballot_sync(_S1448, true);
        v_0 = _S1561;
        _S1448 = _S1562;
    }
    else
    {
        uint _S1563 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1563;
    }
    *&((&_S1553)->y) = v_0;
    bool _S1564 = (F32_isfinite((_S1553.z)));
    uint _S1565 = __ballot_sync(_S1448, _S1564);
    if(_S1564)
    {
        float _S1566 = _S1553.z;
        uint _S1567 = __ballot_sync(_S1448, true);
        v_0 = _S1566;
        _S1448 = _S1567;
    }
    else
    {
        uint _S1568 = __ballot_sync(_S1448, true);
        v_0 = 0.0f;
        _S1448 = _S1568;
    }
    *&((&_S1553)->z) = v_0;
    float _S1569 = WaveActiveSum_0(_S1553.x, _S1448);
    *&((&_S1553)->x) = _S1569;
    float _S1570 = WaveActiveSum_0(_S1553.y, _S1448);
    *&((&_S1553)->y) = _S1570;
    float _S1571 = WaveActiveSum_0(_S1553.z, _S1448);
    *&((&_S1553)->z) = _S1571;
    uint warp_id_2 = thread_id_1 / 32U;
    uint _S1572 = __ballot_sync(_S1448, _S1465);
    if(_S1465)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_2] = _S1553.x;
        (*&_sh_block_reduce_shared_0)[warp_id_2 + 16U] = _S1553.y;
        (*&_sh_block_reduce_shared_0)[warp_id_2 + 32U] = _S1553.z;
        uint _S1573 = __ballot_sync(_S1448, true);
        _S1448 = _S1573;
    }
    else
    {
        uint _S1574 = __ballot_sync(_S1448, true);
        _S1448 = _S1574;
    }
    __syncthreads();
    bool _S1575 = warp_id_2 < 3U;
    uint _S1576 = __ballot_sync(_S1448, _S1575);
    if(_S1575)
    {
        bool _S1577 = lane_id_0 < 16U;
        uint _S1578 = __ballot_sync(_S1576, _S1577);
        if(_S1577)
        {
            float _S1579 = (*&_sh_block_reduce_shared_0)[lane_id_0 + warp_id_2 * 16U];
            uint _S1580 = __ballot_sync(_S1576, true);
            v_0 = _S1579;
            _S1448 = _S1580;
        }
        else
        {
            uint _S1581 = __ballot_sync(_S1576, true);
            v_0 = 0.0f;
            _S1448 = _S1581;
        }
        float _S1582 = WaveActiveSum_0(v_0, _S1448);
        uint _S1583 = __ballot_sync(_S1448, _S1465);
        if(_S1465)
        {
            _S1472 = _S1582 != 0.0f;
        }
        else
        {
            _S1472 = false;
        }
        if(_S1472)
        {
            if(warp_id_2 == 0U)
            {
                float _S1584 = atomicAdd(&((v_coeffs_25 + 2U)->x), _S1582);
            }
            else
            {
                if(warp_id_2 == 1U)
                {
                    float _S1585 = atomicAdd(&((v_coeffs_25 + 2U)->y), _S1582);
                }
                else
                {
                    float _S1586 = atomicAdd(&((v_coeffs_25 + 2U)->z), _S1582);
                }
            }
        }
    }
    float3  dir_n_18 = make_float3 (x_28, y_26, z_25);
    float3  v_dir_n_34 = make_float3 (-0.48860251903533936f * dot_0(*(coeffs_37 + int(2)), v_colors_25), -0.48860251903533936f * dot_0(*(coeffs_37 + int(0)), v_colors_25), 0.48860251903533936f * dot_0(*(coeffs_37 + int(1)), v_colors_25));
    *v_dir_5 = *v_dir_5 + (v_dir_n_34 - make_float3 (dot_0(v_dir_n_34, dir_n_18)) * dir_n_18) * make_float3 (inorm_2);
    return;
}

inline __device__ float3  sh2_to_color_dir(float3  dir_8, float3  coeff_dc_38, float3  * coeffs_38)
{
    float _S1587 = dir_8.x;
    float _S1588 = dir_8.y;
    float _S1589 = dir_8.z;
    float inv_norm_25 = (F32_rsqrt((_S1587 * _S1587 + _S1588 * _S1588 + _S1589 * _S1589)));
    float x_29 = _S1587 * inv_norm_25;
    float y_27 = _S1588 * inv_norm_25;
    float z_26 = _S1589 * inv_norm_25;
    float fTmp0B_18 = -1.09254848957061768f * z_26;
    return make_float3 (0.282094806432724f) * coeff_dc_38 + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * *(coeffs_38 + int(0)) + make_float3 (z_26) * *(coeffs_38 + int(1)) - make_float3 (x_29) * *(coeffs_38 + int(2))) + (make_float3 (0.54627424478530884f * (2.0f * x_29 * y_27)) * *(coeffs_38 + int(3)) + make_float3 (fTmp0B_18 * y_27) * *(coeffs_38 + int(4)) + make_float3 (0.94617468118667603f * (z_26 * z_26) - 0.31539157032966614f) * *(coeffs_38 + int(5)) + make_float3 (fTmp0B_18 * x_29) * *(coeffs_38 + int(6)) + make_float3 (0.54627424478530884f * (x_29 * x_29 - y_27 * y_27)) * *(coeffs_38 + int(7)));
}

inline __device__ void sh2_to_color_dir_vjp_inplace(float3  dir_9, float3  coeff_dc_39, float3  * coeffs_39, float3  v_colors_26, float3  * v_coeff_dc_26, float3  * v_coeffs_26, float3  * v_dir_6)
{
    *v_coeff_dc_26 = *v_coeff_dc_26 + make_float3 (0.282094806432724f) * v_colors_26;
    float _S1590 = dir_9.x;
    float _S1591 = dir_9.y;
    float _S1592 = dir_9.z;
    float inorm_3 = (F32_rsqrt((_S1590 * _S1590 + _S1591 * _S1591 + _S1592 * _S1592)));
    float x_30 = _S1590 * inorm_3;
    float y_28 = _S1591 * inorm_3;
    float z_27 = _S1592 * inorm_3;
    float3  * _S1593 = v_coeffs_26 + int(0);
    *_S1593 = *_S1593 + make_float3 (-0.48860251903533936f * y_28) * v_colors_26;
    float3  * _S1594 = v_coeffs_26 + int(1);
    *_S1594 = *_S1594 + make_float3 (0.48860251903533936f * z_27) * v_colors_26;
    float3  * _S1595 = v_coeffs_26 + int(2);
    *_S1595 = *_S1595 + make_float3 (-0.48860251903533936f * x_30) * v_colors_26;
    float _S1596 = -0.48860251903533936f * dot_0(*(coeffs_39 + int(2)), v_colors_26);
    float _S1597 = -0.48860251903533936f * dot_0(*(coeffs_39 + int(0)), v_colors_26);
    float _S1598 = 0.48860251903533936f * dot_0(*(coeffs_39 + int(1)), v_colors_26);
    float fTmp0B_19 = -1.09254848957061768f * z_27;
    float _S1599 = 2.0f * x_30;
    float pSH6_16 = 0.94617468118667603f * (z_27 * z_27) - 0.31539157032966614f;
    float pSH7_15 = fTmp0B_19 * x_30;
    float pSH5_15 = fTmp0B_19 * y_28;
    float pSH8_15 = 0.54627424478530884f * (x_30 * x_30 - y_28 * y_28);
    float3  * _S1600 = v_coeffs_26 + int(3);
    *_S1600 = *_S1600 + make_float3 (0.54627424478530884f * (_S1599 * y_28)) * v_colors_26;
    float3  * _S1601 = v_coeffs_26 + int(4);
    *_S1601 = *_S1601 + make_float3 (pSH5_15) * v_colors_26;
    float3  * _S1602 = v_coeffs_26 + int(5);
    *_S1602 = *_S1602 + make_float3 (pSH6_16) * v_colors_26;
    float3  * _S1603 = v_coeffs_26 + int(6);
    *_S1603 = *_S1603 + make_float3 (pSH7_15) * v_colors_26;
    float3  * _S1604 = v_coeffs_26 + int(7);
    *_S1604 = *_S1604 + make_float3 (pSH8_15) * v_colors_26;
    float pSH8_x_12 = 0.54627424478530884f * _S1599;
    float3  * _S1605 = coeffs_39 + int(3);
    float3  * _S1606 = coeffs_39 + int(7);
    float3  * _S1607 = coeffs_39 + int(6);
    float3  * _S1608 = coeffs_39 + int(4);
    float3  dir_n_19 = make_float3 (x_30, y_28, z_27);
    float3  v_dir_n_35 = make_float3 (_S1596 + dot_0(v_colors_26, make_float3 (0.54627424478530884f * (2.0f * y_28)) * *_S1605 + make_float3 (pSH8_x_12) * *_S1606 + make_float3 (fTmp0B_19) * *_S1607), _S1597 + dot_0(v_colors_26, make_float3 (pSH8_x_12) * *_S1605 + make_float3 (0.54627424478530884f * (-2.0f * y_28)) * *_S1606 + make_float3 (fTmp0B_19) * *_S1608), _S1598 + dot_0(v_colors_26, make_float3 (1.89234936237335205f * z_27) * *(coeffs_39 + int(5)) + make_float3 (-1.09254848957061768f * x_30) * *_S1607 + make_float3 (-1.09254848957061768f * y_28) * *_S1608));
    *v_dir_6 = *v_dir_6 + (v_dir_n_35 - make_float3 (dot_0(v_dir_n_35, dir_n_19)) * dir_n_19) * make_float3 (inorm_3);
    return;
}

inline __device__ void sh2_to_color_dir_vjp_atomic(float3  dir_10, float3  coeff_dc_40, float3  * coeffs_40, float3  v_colors_27, float3  * v_coeff_dc_27, float3  * v_coeffs_27, float3  * v_dir_7)
{
    *v_coeff_dc_27 = *v_coeff_dc_27 + make_float3 (0.282094806432724f) * v_colors_27;
    float _S1609 = dir_10.x;
    float _S1610 = dir_10.y;
    float _S1611 = dir_10.z;
    float inorm_4 = (F32_rsqrt((_S1609 * _S1609 + _S1610 * _S1610 + _S1611 * _S1611)));
    float x_31 = _S1609 * inorm_4;
    float y_29 = _S1610 * inorm_4;
    float z_28 = _S1611 * inorm_4;
    float3  temp_203 = make_float3 (-0.48860251903533936f * y_29) * v_colors_27;
    float _S1612 = dot_0(temp_203, temp_203);
    bool _S1613;
    if((F32_isfinite((_S1612))))
    {
        _S1613 = _S1612 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1614 = v_coeffs_27 + int(0);
        float _S1615 = atomicAdd(&(_S1614->x), temp_203.x);
        float _S1616 = atomicAdd(&(_S1614->y), temp_203.y);
        float _S1617 = atomicAdd(&(_S1614->z), temp_203.z);
    }
    float3  temp_204 = make_float3 (0.48860251903533936f * z_28) * v_colors_27;
    float _S1618 = dot_0(temp_204, temp_204);
    if((F32_isfinite((_S1618))))
    {
        _S1613 = _S1618 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1619 = v_coeffs_27 + int(1);
        float _S1620 = atomicAdd(&(_S1619->x), temp_204.x);
        float _S1621 = atomicAdd(&(_S1619->y), temp_204.y);
        float _S1622 = atomicAdd(&(_S1619->z), temp_204.z);
    }
    float3  temp_205 = make_float3 (-0.48860251903533936f * x_31) * v_colors_27;
    float _S1623 = dot_0(temp_205, temp_205);
    if((F32_isfinite((_S1623))))
    {
        _S1613 = _S1623 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1624 = v_coeffs_27 + int(2);
        float _S1625 = atomicAdd(&(_S1624->x), temp_205.x);
        float _S1626 = atomicAdd(&(_S1624->y), temp_205.y);
        float _S1627 = atomicAdd(&(_S1624->z), temp_205.z);
    }
    float _S1628 = -0.48860251903533936f * dot_0(*(coeffs_40 + int(2)), v_colors_27);
    float _S1629 = -0.48860251903533936f * dot_0(*(coeffs_40 + int(0)), v_colors_27);
    float _S1630 = 0.48860251903533936f * dot_0(*(coeffs_40 + int(1)), v_colors_27);
    float fTmp0B_20 = -1.09254848957061768f * z_28;
    float _S1631 = 2.0f * x_31;
    float pSH6_17 = 0.94617468118667603f * (z_28 * z_28) - 0.31539157032966614f;
    float pSH7_16 = fTmp0B_20 * x_31;
    float pSH5_16 = fTmp0B_20 * y_29;
    float pSH8_16 = 0.54627424478530884f * (x_31 * x_31 - y_29 * y_29);
    float3  temp_206 = make_float3 (0.54627424478530884f * (_S1631 * y_29)) * v_colors_27;
    float _S1632 = dot_0(temp_206, temp_206);
    if((F32_isfinite((_S1632))))
    {
        _S1613 = _S1632 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1633 = v_coeffs_27 + int(3);
        float _S1634 = atomicAdd(&(_S1633->x), temp_206.x);
        float _S1635 = atomicAdd(&(_S1633->y), temp_206.y);
        float _S1636 = atomicAdd(&(_S1633->z), temp_206.z);
    }
    float3  temp_207 = make_float3 (pSH5_16) * v_colors_27;
    float _S1637 = dot_0(temp_207, temp_207);
    if((F32_isfinite((_S1637))))
    {
        _S1613 = _S1637 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1638 = v_coeffs_27 + int(4);
        float _S1639 = atomicAdd(&(_S1638->x), temp_207.x);
        float _S1640 = atomicAdd(&(_S1638->y), temp_207.y);
        float _S1641 = atomicAdd(&(_S1638->z), temp_207.z);
    }
    float3  temp_208 = make_float3 (pSH6_17) * v_colors_27;
    float _S1642 = dot_0(temp_208, temp_208);
    if((F32_isfinite((_S1642))))
    {
        _S1613 = _S1642 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1643 = v_coeffs_27 + int(5);
        float _S1644 = atomicAdd(&(_S1643->x), temp_208.x);
        float _S1645 = atomicAdd(&(_S1643->y), temp_208.y);
        float _S1646 = atomicAdd(&(_S1643->z), temp_208.z);
    }
    float3  temp_209 = make_float3 (pSH7_16) * v_colors_27;
    float _S1647 = dot_0(temp_209, temp_209);
    if((F32_isfinite((_S1647))))
    {
        _S1613 = _S1647 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1648 = v_coeffs_27 + int(6);
        float _S1649 = atomicAdd(&(_S1648->x), temp_209.x);
        float _S1650 = atomicAdd(&(_S1648->y), temp_209.y);
        float _S1651 = atomicAdd(&(_S1648->z), temp_209.z);
    }
    float3  temp_210 = make_float3 (pSH8_16) * v_colors_27;
    float _S1652 = dot_0(temp_210, temp_210);
    if((F32_isfinite((_S1652))))
    {
        _S1613 = _S1652 != 0.0f;
    }
    else
    {
        _S1613 = false;
    }
    if(_S1613)
    {
        float3  * _S1653 = v_coeffs_27 + int(7);
        float _S1654 = atomicAdd(&(_S1653->x), temp_210.x);
        float _S1655 = atomicAdd(&(_S1653->y), temp_210.y);
        float _S1656 = atomicAdd(&(_S1653->z), temp_210.z);
    }
    float pSH8_x_13 = 0.54627424478530884f * _S1631;
    float3  * _S1657 = coeffs_40 + int(3);
    float3  * _S1658 = coeffs_40 + int(7);
    float3  * _S1659 = coeffs_40 + int(6);
    float3  * _S1660 = coeffs_40 + int(4);
    float3  dir_n_20 = make_float3 (x_31, y_29, z_28);
    float3  v_dir_n_36 = make_float3 (_S1628 + dot_0(v_colors_27, make_float3 (0.54627424478530884f * (2.0f * y_29)) * *_S1657 + make_float3 (pSH8_x_13) * *_S1658 + make_float3 (fTmp0B_20) * *_S1659), _S1629 + dot_0(v_colors_27, make_float3 (pSH8_x_13) * *_S1657 + make_float3 (0.54627424478530884f * (-2.0f * y_29)) * *_S1658 + make_float3 (fTmp0B_20) * *_S1660), _S1630 + dot_0(v_colors_27, make_float3 (1.89234936237335205f * z_28) * *(coeffs_40 + int(5)) + make_float3 (-1.09254848957061768f * x_31) * *_S1659 + make_float3 (-1.09254848957061768f * y_29) * *_S1660));
    *v_dir_7 = *v_dir_7 + (v_dir_n_36 - make_float3 (dot_0(v_dir_n_36, dir_n_20)) * dir_n_20) * make_float3 (inorm_4);
    return;
}

inline __device__ void sh2_to_color_dir_vjp_block_atomic(float3  dir_11, float3  coeff_dc_41, float3  * coeffs_41, float3  v_colors_28, float3  * v_coeff_dc_28, float3  * v_coeffs_28, float3  * v_dir_8, uint thread_id_2, uint _S1661)
{
    *v_coeff_dc_28 = *v_coeff_dc_28 + make_float3 (0.282094806432724f) * v_colors_28;
    float _S1662 = dir_11.x;
    float _S1663 = dir_11.y;
    float _S1664 = dir_11.z;
    float inorm_5 = (F32_rsqrt((_S1662 * _S1662 + _S1663 * _S1663 + _S1664 * _S1664)));
    float x_32 = _S1662 * inorm_5;
    float y_30 = _S1663 * inorm_5;
    float z_29 = _S1664 * inorm_5;
    float3  _S1665 = make_float3 (-0.48860251903533936f * y_30) * v_colors_28;
    float3  _S1666 = _S1665;
    bool _S1667 = (F32_isfinite((_S1665.x)));
    uint _S1668 = __ballot_sync(_S1661, _S1667);
    float v_1;
    uint _S1669;
    if(_S1667)
    {
        float _S1670 = _S1666.x;
        uint _S1671 = __ballot_sync(_S1661, true);
        v_1 = _S1670;
        _S1669 = _S1671;
    }
    else
    {
        uint _S1672 = __ballot_sync(_S1661, true);
        v_1 = 0.0f;
        _S1669 = _S1672;
    }
    *&((&_S1666)->x) = v_1;
    bool _S1673 = (F32_isfinite((_S1666.y)));
    uint _S1674 = __ballot_sync(_S1669, _S1673);
    if(_S1673)
    {
        float _S1675 = _S1666.y;
        uint _S1676 = __ballot_sync(_S1669, true);
        v_1 = _S1675;
        _S1669 = _S1676;
    }
    else
    {
        uint _S1677 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1677;
    }
    *&((&_S1666)->y) = v_1;
    bool _S1678 = (F32_isfinite((_S1666.z)));
    uint _S1679 = __ballot_sync(_S1669, _S1678);
    if(_S1678)
    {
        float _S1680 = _S1666.z;
        uint _S1681 = __ballot_sync(_S1669, true);
        v_1 = _S1680;
        _S1669 = _S1681;
    }
    else
    {
        uint _S1682 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1682;
    }
    *&((&_S1666)->z) = v_1;
    float _S1683 = WaveActiveSum_0(_S1666.x, _S1669);
    *&((&_S1666)->x) = _S1683;
    float _S1684 = WaveActiveSum_0(_S1666.y, _S1669);
    *&((&_S1666)->y) = _S1684;
    float _S1685 = WaveActiveSum_0(_S1666.z, _S1669);
    *&((&_S1666)->z) = _S1685;
    uint warp_id_3 = thread_id_2 / 32U;
    uint lane_id_1 = thread_id_2 % 32U;
    bool _S1686 = lane_id_1 == 0U;
    uint _S1687 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_3] = _S1666.x;
        (*&_sh_block_reduce_shared_0)[warp_id_3 + 16U] = _S1666.y;
        (*&_sh_block_reduce_shared_0)[warp_id_3 + 32U] = _S1666.z;
        uint _S1688 = __ballot_sync(_S1669, true);
        _S1669 = _S1688;
    }
    else
    {
        uint _S1689 = __ballot_sync(_S1669, true);
        _S1669 = _S1689;
    }
    __syncthreads();
    bool _S1690 = warp_id_3 < 3U;
    uint _S1691 = __ballot_sync(_S1669, _S1690);
    uint _S1692;
    bool _S1693;
    if(_S1690)
    {
        bool _S1694 = lane_id_1 < 16U;
        uint _S1695 = __ballot_sync(_S1691, _S1694);
        if(_S1694)
        {
            float _S1696 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_3 * 16U];
            uint _S1697 = __ballot_sync(_S1691, true);
            v_1 = _S1696;
            _S1692 = _S1697;
        }
        else
        {
            uint _S1698 = __ballot_sync(_S1691, true);
            v_1 = 0.0f;
            _S1692 = _S1698;
        }
        float _S1699 = WaveActiveSum_0(v_1, _S1692);
        uint _S1700 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1701 = _S1699 != 0.0f;
            uint _S1702 = __ballot_sync(_S1692, true);
            _S1693 = _S1701;
            _S1692 = _S1702;
        }
        else
        {
            uint _S1703 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1703;
        }
        uint _S1704 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1705 = warp_id_3 == 0U;
            uint _S1706 = __ballot_sync(_S1704, _S1705);
            if(_S1705)
            {
                float _S1707 = atomicAdd(&((v_coeffs_28 + 0U)->x), _S1699);
                uint _S1708 = __ballot_sync(_S1704, true);
            }
            else
            {
                uint _S1709 = _S1704 & (~_S1706);
                bool _S1710 = warp_id_3 == 1U;
                uint _S1711 = __ballot_sync(_S1709, _S1710);
                if(_S1710)
                {
                    float _S1712 = atomicAdd(&((v_coeffs_28 + 0U)->y), _S1699);
                    uint _S1713 = __ballot_sync(_S1709, true);
                }
                else
                {
                    float _S1714 = atomicAdd(&((v_coeffs_28 + 0U)->z), _S1699);
                    uint _S1715 = __ballot_sync(_S1709, true);
                }
                uint _S1716 = __ballot_sync(_S1704, true);
            }
            uint _S1717 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1718 = __ballot_sync(_S1692, true);
        }
        uint _S1719 = __ballot_sync(_S1669, true);
        _S1669 = _S1719;
    }
    else
    {
        uint _S1720 = __ballot_sync(_S1669, true);
        _S1669 = _S1720;
    }
    float3  _S1721 = make_float3 (0.48860251903533936f * z_29) * v_colors_28;
    float3  _S1722 = _S1721;
    bool _S1723 = (F32_isfinite((_S1721.x)));
    uint _S1724 = __ballot_sync(_S1669, _S1723);
    if(_S1723)
    {
        float _S1725 = _S1722.x;
        uint _S1726 = __ballot_sync(_S1669, true);
        v_1 = _S1725;
        _S1669 = _S1726;
    }
    else
    {
        uint _S1727 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1727;
    }
    *&((&_S1722)->x) = v_1;
    bool _S1728 = (F32_isfinite((_S1722.y)));
    uint _S1729 = __ballot_sync(_S1669, _S1728);
    if(_S1728)
    {
        float _S1730 = _S1722.y;
        uint _S1731 = __ballot_sync(_S1669, true);
        v_1 = _S1730;
        _S1669 = _S1731;
    }
    else
    {
        uint _S1732 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1732;
    }
    *&((&_S1722)->y) = v_1;
    bool _S1733 = (F32_isfinite((_S1722.z)));
    uint _S1734 = __ballot_sync(_S1669, _S1733);
    if(_S1733)
    {
        float _S1735 = _S1722.z;
        uint _S1736 = __ballot_sync(_S1669, true);
        v_1 = _S1735;
        _S1669 = _S1736;
    }
    else
    {
        uint _S1737 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1737;
    }
    *&((&_S1722)->z) = v_1;
    float _S1738 = WaveActiveSum_0(_S1722.x, _S1669);
    *&((&_S1722)->x) = _S1738;
    float _S1739 = WaveActiveSum_0(_S1722.y, _S1669);
    *&((&_S1722)->y) = _S1739;
    float _S1740 = WaveActiveSum_0(_S1722.z, _S1669);
    *&((&_S1722)->z) = _S1740;
    uint warp_id_4 = thread_id_2 / 32U;
    uint _S1741 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_4] = _S1722.x;
        (*&_sh_block_reduce_shared_0)[warp_id_4 + 16U] = _S1722.y;
        (*&_sh_block_reduce_shared_0)[warp_id_4 + 32U] = _S1722.z;
        uint _S1742 = __ballot_sync(_S1669, true);
        _S1669 = _S1742;
    }
    else
    {
        uint _S1743 = __ballot_sync(_S1669, true);
        _S1669 = _S1743;
    }
    __syncthreads();
    bool _S1744 = warp_id_4 < 3U;
    uint _S1745 = __ballot_sync(_S1669, _S1744);
    if(_S1744)
    {
        bool _S1746 = lane_id_1 < 16U;
        uint _S1747 = __ballot_sync(_S1745, _S1746);
        if(_S1746)
        {
            float _S1748 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_4 * 16U];
            uint _S1749 = __ballot_sync(_S1745, true);
            v_1 = _S1748;
            _S1692 = _S1749;
        }
        else
        {
            uint _S1750 = __ballot_sync(_S1745, true);
            v_1 = 0.0f;
            _S1692 = _S1750;
        }
        float _S1751 = WaveActiveSum_0(v_1, _S1692);
        uint _S1752 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1753 = _S1751 != 0.0f;
            uint _S1754 = __ballot_sync(_S1692, true);
            _S1693 = _S1753;
            _S1692 = _S1754;
        }
        else
        {
            uint _S1755 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1755;
        }
        uint _S1756 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1757 = warp_id_4 == 0U;
            uint _S1758 = __ballot_sync(_S1756, _S1757);
            if(_S1757)
            {
                float _S1759 = atomicAdd(&((v_coeffs_28 + 1U)->x), _S1751);
                uint _S1760 = __ballot_sync(_S1756, true);
            }
            else
            {
                uint _S1761 = _S1756 & (~_S1758);
                bool _S1762 = warp_id_4 == 1U;
                uint _S1763 = __ballot_sync(_S1761, _S1762);
                if(_S1762)
                {
                    float _S1764 = atomicAdd(&((v_coeffs_28 + 1U)->y), _S1751);
                    uint _S1765 = __ballot_sync(_S1761, true);
                }
                else
                {
                    float _S1766 = atomicAdd(&((v_coeffs_28 + 1U)->z), _S1751);
                    uint _S1767 = __ballot_sync(_S1761, true);
                }
                uint _S1768 = __ballot_sync(_S1756, true);
            }
            uint _S1769 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1770 = __ballot_sync(_S1692, true);
        }
        uint _S1771 = __ballot_sync(_S1669, true);
        _S1669 = _S1771;
    }
    else
    {
        uint _S1772 = __ballot_sync(_S1669, true);
        _S1669 = _S1772;
    }
    float3  _S1773 = make_float3 (-0.48860251903533936f * x_32) * v_colors_28;
    float3  _S1774 = _S1773;
    bool _S1775 = (F32_isfinite((_S1773.x)));
    uint _S1776 = __ballot_sync(_S1669, _S1775);
    if(_S1775)
    {
        float _S1777 = _S1774.x;
        uint _S1778 = __ballot_sync(_S1669, true);
        v_1 = _S1777;
        _S1669 = _S1778;
    }
    else
    {
        uint _S1779 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1779;
    }
    *&((&_S1774)->x) = v_1;
    bool _S1780 = (F32_isfinite((_S1774.y)));
    uint _S1781 = __ballot_sync(_S1669, _S1780);
    if(_S1780)
    {
        float _S1782 = _S1774.y;
        uint _S1783 = __ballot_sync(_S1669, true);
        v_1 = _S1782;
        _S1669 = _S1783;
    }
    else
    {
        uint _S1784 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1784;
    }
    *&((&_S1774)->y) = v_1;
    bool _S1785 = (F32_isfinite((_S1774.z)));
    uint _S1786 = __ballot_sync(_S1669, _S1785);
    if(_S1785)
    {
        float _S1787 = _S1774.z;
        uint _S1788 = __ballot_sync(_S1669, true);
        v_1 = _S1787;
        _S1669 = _S1788;
    }
    else
    {
        uint _S1789 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1789;
    }
    *&((&_S1774)->z) = v_1;
    float _S1790 = WaveActiveSum_0(_S1774.x, _S1669);
    *&((&_S1774)->x) = _S1790;
    float _S1791 = WaveActiveSum_0(_S1774.y, _S1669);
    *&((&_S1774)->y) = _S1791;
    float _S1792 = WaveActiveSum_0(_S1774.z, _S1669);
    *&((&_S1774)->z) = _S1792;
    uint warp_id_5 = thread_id_2 / 32U;
    uint _S1793 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_5] = _S1774.x;
        (*&_sh_block_reduce_shared_0)[warp_id_5 + 16U] = _S1774.y;
        (*&_sh_block_reduce_shared_0)[warp_id_5 + 32U] = _S1774.z;
        uint _S1794 = __ballot_sync(_S1669, true);
        _S1669 = _S1794;
    }
    else
    {
        uint _S1795 = __ballot_sync(_S1669, true);
        _S1669 = _S1795;
    }
    __syncthreads();
    bool _S1796 = warp_id_5 < 3U;
    uint _S1797 = __ballot_sync(_S1669, _S1796);
    if(_S1796)
    {
        bool _S1798 = lane_id_1 < 16U;
        uint _S1799 = __ballot_sync(_S1797, _S1798);
        if(_S1798)
        {
            float _S1800 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_5 * 16U];
            uint _S1801 = __ballot_sync(_S1797, true);
            v_1 = _S1800;
            _S1692 = _S1801;
        }
        else
        {
            uint _S1802 = __ballot_sync(_S1797, true);
            v_1 = 0.0f;
            _S1692 = _S1802;
        }
        float _S1803 = WaveActiveSum_0(v_1, _S1692);
        uint _S1804 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1805 = _S1803 != 0.0f;
            uint _S1806 = __ballot_sync(_S1692, true);
            _S1693 = _S1805;
            _S1692 = _S1806;
        }
        else
        {
            uint _S1807 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1807;
        }
        uint _S1808 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1809 = warp_id_5 == 0U;
            uint _S1810 = __ballot_sync(_S1808, _S1809);
            if(_S1809)
            {
                float _S1811 = atomicAdd(&((v_coeffs_28 + 2U)->x), _S1803);
                uint _S1812 = __ballot_sync(_S1808, true);
            }
            else
            {
                uint _S1813 = _S1808 & (~_S1810);
                bool _S1814 = warp_id_5 == 1U;
                uint _S1815 = __ballot_sync(_S1813, _S1814);
                if(_S1814)
                {
                    float _S1816 = atomicAdd(&((v_coeffs_28 + 2U)->y), _S1803);
                    uint _S1817 = __ballot_sync(_S1813, true);
                }
                else
                {
                    float _S1818 = atomicAdd(&((v_coeffs_28 + 2U)->z), _S1803);
                    uint _S1819 = __ballot_sync(_S1813, true);
                }
                uint _S1820 = __ballot_sync(_S1808, true);
            }
            uint _S1821 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1822 = __ballot_sync(_S1692, true);
        }
        uint _S1823 = __ballot_sync(_S1669, true);
        _S1669 = _S1823;
    }
    else
    {
        uint _S1824 = __ballot_sync(_S1669, true);
        _S1669 = _S1824;
    }
    float _S1825 = -0.48860251903533936f * dot_0(*(coeffs_41 + int(2)), v_colors_28);
    float _S1826 = -0.48860251903533936f * dot_0(*(coeffs_41 + int(0)), v_colors_28);
    float _S1827 = 0.48860251903533936f * dot_0(*(coeffs_41 + int(1)), v_colors_28);
    float fTmp0B_21 = -1.09254848957061768f * z_29;
    float _S1828 = 2.0f * x_32;
    float pSH6_18 = 0.94617468118667603f * (z_29 * z_29) - 0.31539157032966614f;
    float pSH7_17 = fTmp0B_21 * x_32;
    float pSH5_17 = fTmp0B_21 * y_30;
    float pSH8_17 = 0.54627424478530884f * (x_32 * x_32 - y_30 * y_30);
    float3  _S1829 = make_float3 (0.54627424478530884f * (_S1828 * y_30)) * v_colors_28;
    float3  _S1830 = _S1829;
    bool _S1831 = (F32_isfinite((_S1829.x)));
    uint _S1832 = __ballot_sync(_S1669, _S1831);
    if(_S1831)
    {
        float _S1833 = _S1830.x;
        uint _S1834 = __ballot_sync(_S1669, true);
        v_1 = _S1833;
        _S1669 = _S1834;
    }
    else
    {
        uint _S1835 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1835;
    }
    *&((&_S1830)->x) = v_1;
    bool _S1836 = (F32_isfinite((_S1830.y)));
    uint _S1837 = __ballot_sync(_S1669, _S1836);
    if(_S1836)
    {
        float _S1838 = _S1830.y;
        uint _S1839 = __ballot_sync(_S1669, true);
        v_1 = _S1838;
        _S1669 = _S1839;
    }
    else
    {
        uint _S1840 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1840;
    }
    *&((&_S1830)->y) = v_1;
    bool _S1841 = (F32_isfinite((_S1830.z)));
    uint _S1842 = __ballot_sync(_S1669, _S1841);
    if(_S1841)
    {
        float _S1843 = _S1830.z;
        uint _S1844 = __ballot_sync(_S1669, true);
        v_1 = _S1843;
        _S1669 = _S1844;
    }
    else
    {
        uint _S1845 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1845;
    }
    *&((&_S1830)->z) = v_1;
    float _S1846 = WaveActiveSum_0(_S1830.x, _S1669);
    *&((&_S1830)->x) = _S1846;
    float _S1847 = WaveActiveSum_0(_S1830.y, _S1669);
    *&((&_S1830)->y) = _S1847;
    float _S1848 = WaveActiveSum_0(_S1830.z, _S1669);
    *&((&_S1830)->z) = _S1848;
    uint warp_id_6 = thread_id_2 / 32U;
    uint _S1849 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_6] = _S1830.x;
        (*&_sh_block_reduce_shared_0)[warp_id_6 + 16U] = _S1830.y;
        (*&_sh_block_reduce_shared_0)[warp_id_6 + 32U] = _S1830.z;
        uint _S1850 = __ballot_sync(_S1669, true);
        _S1669 = _S1850;
    }
    else
    {
        uint _S1851 = __ballot_sync(_S1669, true);
        _S1669 = _S1851;
    }
    __syncthreads();
    bool _S1852 = warp_id_6 < 3U;
    uint _S1853 = __ballot_sync(_S1669, _S1852);
    if(_S1852)
    {
        bool _S1854 = lane_id_1 < 16U;
        uint _S1855 = __ballot_sync(_S1853, _S1854);
        if(_S1854)
        {
            float _S1856 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_6 * 16U];
            uint _S1857 = __ballot_sync(_S1853, true);
            v_1 = _S1856;
            _S1692 = _S1857;
        }
        else
        {
            uint _S1858 = __ballot_sync(_S1853, true);
            v_1 = 0.0f;
            _S1692 = _S1858;
        }
        float _S1859 = WaveActiveSum_0(v_1, _S1692);
        uint _S1860 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1861 = _S1859 != 0.0f;
            uint _S1862 = __ballot_sync(_S1692, true);
            _S1693 = _S1861;
            _S1692 = _S1862;
        }
        else
        {
            uint _S1863 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1863;
        }
        uint _S1864 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1865 = warp_id_6 == 0U;
            uint _S1866 = __ballot_sync(_S1864, _S1865);
            if(_S1865)
            {
                float _S1867 = atomicAdd(&((v_coeffs_28 + 3U)->x), _S1859);
                uint _S1868 = __ballot_sync(_S1864, true);
            }
            else
            {
                uint _S1869 = _S1864 & (~_S1866);
                bool _S1870 = warp_id_6 == 1U;
                uint _S1871 = __ballot_sync(_S1869, _S1870);
                if(_S1870)
                {
                    float _S1872 = atomicAdd(&((v_coeffs_28 + 3U)->y), _S1859);
                    uint _S1873 = __ballot_sync(_S1869, true);
                }
                else
                {
                    float _S1874 = atomicAdd(&((v_coeffs_28 + 3U)->z), _S1859);
                    uint _S1875 = __ballot_sync(_S1869, true);
                }
                uint _S1876 = __ballot_sync(_S1864, true);
            }
            uint _S1877 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1878 = __ballot_sync(_S1692, true);
        }
        uint _S1879 = __ballot_sync(_S1669, true);
        _S1669 = _S1879;
    }
    else
    {
        uint _S1880 = __ballot_sync(_S1669, true);
        _S1669 = _S1880;
    }
    float3  _S1881 = make_float3 (pSH5_17) * v_colors_28;
    float3  _S1882 = _S1881;
    bool _S1883 = (F32_isfinite((_S1881.x)));
    uint _S1884 = __ballot_sync(_S1669, _S1883);
    if(_S1883)
    {
        float _S1885 = _S1882.x;
        uint _S1886 = __ballot_sync(_S1669, true);
        v_1 = _S1885;
        _S1669 = _S1886;
    }
    else
    {
        uint _S1887 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1887;
    }
    *&((&_S1882)->x) = v_1;
    bool _S1888 = (F32_isfinite((_S1882.y)));
    uint _S1889 = __ballot_sync(_S1669, _S1888);
    if(_S1888)
    {
        float _S1890 = _S1882.y;
        uint _S1891 = __ballot_sync(_S1669, true);
        v_1 = _S1890;
        _S1669 = _S1891;
    }
    else
    {
        uint _S1892 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1892;
    }
    *&((&_S1882)->y) = v_1;
    bool _S1893 = (F32_isfinite((_S1882.z)));
    uint _S1894 = __ballot_sync(_S1669, _S1893);
    if(_S1893)
    {
        float _S1895 = _S1882.z;
        uint _S1896 = __ballot_sync(_S1669, true);
        v_1 = _S1895;
        _S1669 = _S1896;
    }
    else
    {
        uint _S1897 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1897;
    }
    *&((&_S1882)->z) = v_1;
    float _S1898 = WaveActiveSum_0(_S1882.x, _S1669);
    *&((&_S1882)->x) = _S1898;
    float _S1899 = WaveActiveSum_0(_S1882.y, _S1669);
    *&((&_S1882)->y) = _S1899;
    float _S1900 = WaveActiveSum_0(_S1882.z, _S1669);
    *&((&_S1882)->z) = _S1900;
    uint warp_id_7 = thread_id_2 / 32U;
    uint _S1901 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_7] = _S1882.x;
        (*&_sh_block_reduce_shared_0)[warp_id_7 + 16U] = _S1882.y;
        (*&_sh_block_reduce_shared_0)[warp_id_7 + 32U] = _S1882.z;
        uint _S1902 = __ballot_sync(_S1669, true);
        _S1669 = _S1902;
    }
    else
    {
        uint _S1903 = __ballot_sync(_S1669, true);
        _S1669 = _S1903;
    }
    __syncthreads();
    bool _S1904 = warp_id_7 < 3U;
    uint _S1905 = __ballot_sync(_S1669, _S1904);
    if(_S1904)
    {
        bool _S1906 = lane_id_1 < 16U;
        uint _S1907 = __ballot_sync(_S1905, _S1906);
        if(_S1906)
        {
            float _S1908 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_7 * 16U];
            uint _S1909 = __ballot_sync(_S1905, true);
            v_1 = _S1908;
            _S1692 = _S1909;
        }
        else
        {
            uint _S1910 = __ballot_sync(_S1905, true);
            v_1 = 0.0f;
            _S1692 = _S1910;
        }
        float _S1911 = WaveActiveSum_0(v_1, _S1692);
        uint _S1912 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1913 = _S1911 != 0.0f;
            uint _S1914 = __ballot_sync(_S1692, true);
            _S1693 = _S1913;
            _S1692 = _S1914;
        }
        else
        {
            uint _S1915 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1915;
        }
        uint _S1916 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1917 = warp_id_7 == 0U;
            uint _S1918 = __ballot_sync(_S1916, _S1917);
            if(_S1917)
            {
                float _S1919 = atomicAdd(&((v_coeffs_28 + 4U)->x), _S1911);
                uint _S1920 = __ballot_sync(_S1916, true);
            }
            else
            {
                uint _S1921 = _S1916 & (~_S1918);
                bool _S1922 = warp_id_7 == 1U;
                uint _S1923 = __ballot_sync(_S1921, _S1922);
                if(_S1922)
                {
                    float _S1924 = atomicAdd(&((v_coeffs_28 + 4U)->y), _S1911);
                    uint _S1925 = __ballot_sync(_S1921, true);
                }
                else
                {
                    float _S1926 = atomicAdd(&((v_coeffs_28 + 4U)->z), _S1911);
                    uint _S1927 = __ballot_sync(_S1921, true);
                }
                uint _S1928 = __ballot_sync(_S1916, true);
            }
            uint _S1929 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1930 = __ballot_sync(_S1692, true);
        }
        uint _S1931 = __ballot_sync(_S1669, true);
        _S1669 = _S1931;
    }
    else
    {
        uint _S1932 = __ballot_sync(_S1669, true);
        _S1669 = _S1932;
    }
    float3  _S1933 = make_float3 (pSH6_18) * v_colors_28;
    float3  _S1934 = _S1933;
    bool _S1935 = (F32_isfinite((_S1933.x)));
    uint _S1936 = __ballot_sync(_S1669, _S1935);
    if(_S1935)
    {
        float _S1937 = _S1934.x;
        uint _S1938 = __ballot_sync(_S1669, true);
        v_1 = _S1937;
        _S1669 = _S1938;
    }
    else
    {
        uint _S1939 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1939;
    }
    *&((&_S1934)->x) = v_1;
    bool _S1940 = (F32_isfinite((_S1934.y)));
    uint _S1941 = __ballot_sync(_S1669, _S1940);
    if(_S1940)
    {
        float _S1942 = _S1934.y;
        uint _S1943 = __ballot_sync(_S1669, true);
        v_1 = _S1942;
        _S1669 = _S1943;
    }
    else
    {
        uint _S1944 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1944;
    }
    *&((&_S1934)->y) = v_1;
    bool _S1945 = (F32_isfinite((_S1934.z)));
    uint _S1946 = __ballot_sync(_S1669, _S1945);
    if(_S1945)
    {
        float _S1947 = _S1934.z;
        uint _S1948 = __ballot_sync(_S1669, true);
        v_1 = _S1947;
        _S1669 = _S1948;
    }
    else
    {
        uint _S1949 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1949;
    }
    *&((&_S1934)->z) = v_1;
    float _S1950 = WaveActiveSum_0(_S1934.x, _S1669);
    *&((&_S1934)->x) = _S1950;
    float _S1951 = WaveActiveSum_0(_S1934.y, _S1669);
    *&((&_S1934)->y) = _S1951;
    float _S1952 = WaveActiveSum_0(_S1934.z, _S1669);
    *&((&_S1934)->z) = _S1952;
    uint warp_id_8 = thread_id_2 / 32U;
    uint _S1953 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_8] = _S1934.x;
        (*&_sh_block_reduce_shared_0)[warp_id_8 + 16U] = _S1934.y;
        (*&_sh_block_reduce_shared_0)[warp_id_8 + 32U] = _S1934.z;
        uint _S1954 = __ballot_sync(_S1669, true);
        _S1669 = _S1954;
    }
    else
    {
        uint _S1955 = __ballot_sync(_S1669, true);
        _S1669 = _S1955;
    }
    __syncthreads();
    bool _S1956 = warp_id_8 < 3U;
    uint _S1957 = __ballot_sync(_S1669, _S1956);
    if(_S1956)
    {
        bool _S1958 = lane_id_1 < 16U;
        uint _S1959 = __ballot_sync(_S1957, _S1958);
        if(_S1958)
        {
            float _S1960 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_8 * 16U];
            uint _S1961 = __ballot_sync(_S1957, true);
            v_1 = _S1960;
            _S1692 = _S1961;
        }
        else
        {
            uint _S1962 = __ballot_sync(_S1957, true);
            v_1 = 0.0f;
            _S1692 = _S1962;
        }
        float _S1963 = WaveActiveSum_0(v_1, _S1692);
        uint _S1964 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S1965 = _S1963 != 0.0f;
            uint _S1966 = __ballot_sync(_S1692, true);
            _S1693 = _S1965;
            _S1692 = _S1966;
        }
        else
        {
            uint _S1967 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S1967;
        }
        uint _S1968 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S1969 = warp_id_8 == 0U;
            uint _S1970 = __ballot_sync(_S1968, _S1969);
            if(_S1969)
            {
                float _S1971 = atomicAdd(&((v_coeffs_28 + 5U)->x), _S1963);
                uint _S1972 = __ballot_sync(_S1968, true);
            }
            else
            {
                uint _S1973 = _S1968 & (~_S1970);
                bool _S1974 = warp_id_8 == 1U;
                uint _S1975 = __ballot_sync(_S1973, _S1974);
                if(_S1974)
                {
                    float _S1976 = atomicAdd(&((v_coeffs_28 + 5U)->y), _S1963);
                    uint _S1977 = __ballot_sync(_S1973, true);
                }
                else
                {
                    float _S1978 = atomicAdd(&((v_coeffs_28 + 5U)->z), _S1963);
                    uint _S1979 = __ballot_sync(_S1973, true);
                }
                uint _S1980 = __ballot_sync(_S1968, true);
            }
            uint _S1981 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S1982 = __ballot_sync(_S1692, true);
        }
        uint _S1983 = __ballot_sync(_S1669, true);
        _S1669 = _S1983;
    }
    else
    {
        uint _S1984 = __ballot_sync(_S1669, true);
        _S1669 = _S1984;
    }
    float3  _S1985 = make_float3 (pSH7_17) * v_colors_28;
    float3  _S1986 = _S1985;
    bool _S1987 = (F32_isfinite((_S1985.x)));
    uint _S1988 = __ballot_sync(_S1669, _S1987);
    if(_S1987)
    {
        float _S1989 = _S1986.x;
        uint _S1990 = __ballot_sync(_S1669, true);
        v_1 = _S1989;
        _S1669 = _S1990;
    }
    else
    {
        uint _S1991 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1991;
    }
    *&((&_S1986)->x) = v_1;
    bool _S1992 = (F32_isfinite((_S1986.y)));
    uint _S1993 = __ballot_sync(_S1669, _S1992);
    if(_S1992)
    {
        float _S1994 = _S1986.y;
        uint _S1995 = __ballot_sync(_S1669, true);
        v_1 = _S1994;
        _S1669 = _S1995;
    }
    else
    {
        uint _S1996 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S1996;
    }
    *&((&_S1986)->y) = v_1;
    bool _S1997 = (F32_isfinite((_S1986.z)));
    uint _S1998 = __ballot_sync(_S1669, _S1997);
    if(_S1997)
    {
        float _S1999 = _S1986.z;
        uint _S2000 = __ballot_sync(_S1669, true);
        v_1 = _S1999;
        _S1669 = _S2000;
    }
    else
    {
        uint _S2001 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S2001;
    }
    *&((&_S1986)->z) = v_1;
    float _S2002 = WaveActiveSum_0(_S1986.x, _S1669);
    *&((&_S1986)->x) = _S2002;
    float _S2003 = WaveActiveSum_0(_S1986.y, _S1669);
    *&((&_S1986)->y) = _S2003;
    float _S2004 = WaveActiveSum_0(_S1986.z, _S1669);
    *&((&_S1986)->z) = _S2004;
    uint warp_id_9 = thread_id_2 / 32U;
    uint _S2005 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_9] = _S1986.x;
        (*&_sh_block_reduce_shared_0)[warp_id_9 + 16U] = _S1986.y;
        (*&_sh_block_reduce_shared_0)[warp_id_9 + 32U] = _S1986.z;
        uint _S2006 = __ballot_sync(_S1669, true);
        _S1669 = _S2006;
    }
    else
    {
        uint _S2007 = __ballot_sync(_S1669, true);
        _S1669 = _S2007;
    }
    __syncthreads();
    bool _S2008 = warp_id_9 < 3U;
    uint _S2009 = __ballot_sync(_S1669, _S2008);
    if(_S2008)
    {
        bool _S2010 = lane_id_1 < 16U;
        uint _S2011 = __ballot_sync(_S2009, _S2010);
        if(_S2010)
        {
            float _S2012 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_9 * 16U];
            uint _S2013 = __ballot_sync(_S2009, true);
            v_1 = _S2012;
            _S1692 = _S2013;
        }
        else
        {
            uint _S2014 = __ballot_sync(_S2009, true);
            v_1 = 0.0f;
            _S1692 = _S2014;
        }
        float _S2015 = WaveActiveSum_0(v_1, _S1692);
        uint _S2016 = __ballot_sync(_S1692, _S1686);
        if(_S1686)
        {
            bool _S2017 = _S2015 != 0.0f;
            uint _S2018 = __ballot_sync(_S1692, true);
            _S1693 = _S2017;
            _S1692 = _S2018;
        }
        else
        {
            uint _S2019 = __ballot_sync(_S1692, true);
            _S1693 = false;
            _S1692 = _S2019;
        }
        uint _S2020 = __ballot_sync(_S1692, _S1693);
        if(_S1693)
        {
            bool _S2021 = warp_id_9 == 0U;
            uint _S2022 = __ballot_sync(_S2020, _S2021);
            if(_S2021)
            {
                float _S2023 = atomicAdd(&((v_coeffs_28 + 6U)->x), _S2015);
                uint _S2024 = __ballot_sync(_S2020, true);
            }
            else
            {
                uint _S2025 = _S2020 & (~_S2022);
                bool _S2026 = warp_id_9 == 1U;
                uint _S2027 = __ballot_sync(_S2025, _S2026);
                if(_S2026)
                {
                    float _S2028 = atomicAdd(&((v_coeffs_28 + 6U)->y), _S2015);
                    uint _S2029 = __ballot_sync(_S2025, true);
                }
                else
                {
                    float _S2030 = atomicAdd(&((v_coeffs_28 + 6U)->z), _S2015);
                    uint _S2031 = __ballot_sync(_S2025, true);
                }
                uint _S2032 = __ballot_sync(_S2020, true);
            }
            uint _S2033 = __ballot_sync(_S1692, true);
        }
        else
        {
            uint _S2034 = __ballot_sync(_S1692, true);
        }
        uint _S2035 = __ballot_sync(_S1669, true);
        _S1669 = _S2035;
    }
    else
    {
        uint _S2036 = __ballot_sync(_S1669, true);
        _S1669 = _S2036;
    }
    float3  _S2037 = make_float3 (pSH8_17) * v_colors_28;
    float3  _S2038 = _S2037;
    bool _S2039 = (F32_isfinite((_S2037.x)));
    uint _S2040 = __ballot_sync(_S1669, _S2039);
    if(_S2039)
    {
        float _S2041 = _S2038.x;
        uint _S2042 = __ballot_sync(_S1669, true);
        v_1 = _S2041;
        _S1669 = _S2042;
    }
    else
    {
        uint _S2043 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S2043;
    }
    *&((&_S2038)->x) = v_1;
    bool _S2044 = (F32_isfinite((_S2038.y)));
    uint _S2045 = __ballot_sync(_S1669, _S2044);
    if(_S2044)
    {
        float _S2046 = _S2038.y;
        uint _S2047 = __ballot_sync(_S1669, true);
        v_1 = _S2046;
        _S1669 = _S2047;
    }
    else
    {
        uint _S2048 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S2048;
    }
    *&((&_S2038)->y) = v_1;
    bool _S2049 = (F32_isfinite((_S2038.z)));
    uint _S2050 = __ballot_sync(_S1669, _S2049);
    if(_S2049)
    {
        float _S2051 = _S2038.z;
        uint _S2052 = __ballot_sync(_S1669, true);
        v_1 = _S2051;
        _S1669 = _S2052;
    }
    else
    {
        uint _S2053 = __ballot_sync(_S1669, true);
        v_1 = 0.0f;
        _S1669 = _S2053;
    }
    *&((&_S2038)->z) = v_1;
    float _S2054 = WaveActiveSum_0(_S2038.x, _S1669);
    *&((&_S2038)->x) = _S2054;
    float _S2055 = WaveActiveSum_0(_S2038.y, _S1669);
    *&((&_S2038)->y) = _S2055;
    float _S2056 = WaveActiveSum_0(_S2038.z, _S1669);
    *&((&_S2038)->z) = _S2056;
    uint warp_id_10 = thread_id_2 / 32U;
    uint _S2057 = __ballot_sync(_S1669, _S1686);
    if(_S1686)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_10] = _S2038.x;
        (*&_sh_block_reduce_shared_0)[warp_id_10 + 16U] = _S2038.y;
        (*&_sh_block_reduce_shared_0)[warp_id_10 + 32U] = _S2038.z;
        uint _S2058 = __ballot_sync(_S1669, true);
        _S1669 = _S2058;
    }
    else
    {
        uint _S2059 = __ballot_sync(_S1669, true);
        _S1669 = _S2059;
    }
    __syncthreads();
    bool _S2060 = warp_id_10 < 3U;
    uint _S2061 = __ballot_sync(_S1669, _S2060);
    if(_S2060)
    {
        bool _S2062 = lane_id_1 < 16U;
        uint _S2063 = __ballot_sync(_S2061, _S2062);
        if(_S2062)
        {
            float _S2064 = (*&_sh_block_reduce_shared_0)[lane_id_1 + warp_id_10 * 16U];
            uint _S2065 = __ballot_sync(_S2061, true);
            v_1 = _S2064;
            _S1669 = _S2065;
        }
        else
        {
            uint _S2066 = __ballot_sync(_S2061, true);
            v_1 = 0.0f;
            _S1669 = _S2066;
        }
        float _S2067 = WaveActiveSum_0(v_1, _S1669);
        uint _S2068 = __ballot_sync(_S1669, _S1686);
        if(_S1686)
        {
            _S1693 = _S2067 != 0.0f;
        }
        else
        {
            _S1693 = false;
        }
        if(_S1693)
        {
            if(warp_id_10 == 0U)
            {
                float _S2069 = atomicAdd(&((v_coeffs_28 + 7U)->x), _S2067);
            }
            else
            {
                if(warp_id_10 == 1U)
                {
                    float _S2070 = atomicAdd(&((v_coeffs_28 + 7U)->y), _S2067);
                }
                else
                {
                    float _S2071 = atomicAdd(&((v_coeffs_28 + 7U)->z), _S2067);
                }
            }
        }
    }
    float pSH8_x_14 = 0.54627424478530884f * _S1828;
    float3  * _S2072 = coeffs_41 + int(3);
    float3  * _S2073 = coeffs_41 + int(7);
    float3  * _S2074 = coeffs_41 + int(6);
    float3  * _S2075 = coeffs_41 + int(4);
    float3  dir_n_21 = make_float3 (x_32, y_30, z_29);
    float3  v_dir_n_37 = make_float3 (_S1825 + dot_0(v_colors_28, make_float3 (0.54627424478530884f * (2.0f * y_30)) * *_S2072 + make_float3 (pSH8_x_14) * *_S2073 + make_float3 (fTmp0B_21) * *_S2074), _S1826 + dot_0(v_colors_28, make_float3 (pSH8_x_14) * *_S2072 + make_float3 (0.54627424478530884f * (-2.0f * y_30)) * *_S2073 + make_float3 (fTmp0B_21) * *_S2075), _S1827 + dot_0(v_colors_28, make_float3 (1.89234936237335205f * z_29) * *(coeffs_41 + int(5)) + make_float3 (-1.09254848957061768f * x_32) * *_S2074 + make_float3 (-1.09254848957061768f * y_30) * *_S2075));
    *v_dir_8 = *v_dir_8 + (v_dir_n_37 - make_float3 (dot_0(v_dir_n_37, dir_n_21)) * dir_n_21) * make_float3 (inorm_5);
    return;
}

inline __device__ float3  sh3_to_color_dir(float3  dir_12, float3  coeff_dc_42, float3  * coeffs_42)
{
    float _S2076 = dir_12.x;
    float _S2077 = dir_12.y;
    float _S2078 = dir_12.z;
    float inv_norm_26 = (F32_rsqrt((_S2076 * _S2076 + _S2077 * _S2077 + _S2078 * _S2078)));
    float x_33 = _S2076 * inv_norm_26;
    float y_31 = _S2077 * inv_norm_26;
    float z_30 = _S2078 * inv_norm_26;
    float z2_12 = z_30 * z_30;
    float fTmp0B_22 = -1.09254848957061768f * z_30;
    float fC1_12 = x_33 * x_33 - y_31 * y_31;
    float fS1_12 = 2.0f * x_33 * y_31;
    float fTmp0C_12 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_30;
    return make_float3 (0.282094806432724f) * coeff_dc_42 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * *(coeffs_42 + int(0)) + make_float3 (z_30) * *(coeffs_42 + int(1)) - make_float3 (x_33) * *(coeffs_42 + int(2))) + (make_float3 (0.54627424478530884f * fS1_12) * *(coeffs_42 + int(3)) + make_float3 (fTmp0B_22 * y_31) * *(coeffs_42 + int(4)) + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * *(coeffs_42 + int(5)) + make_float3 (fTmp0B_22 * x_33) * *(coeffs_42 + int(6)) + make_float3 (0.54627424478530884f * fC1_12) * *(coeffs_42 + int(7))) + (make_float3 (-0.59004360437393188f * (x_33 * fS1_12 + y_31 * fC1_12)) * *(coeffs_42 + int(8)) + make_float3 (fTmp1B_12 * fS1_12) * *(coeffs_42 + int(9)) + make_float3 (fTmp0C_12 * y_31) * *(coeffs_42 + int(10)) + make_float3 (z_30 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * *(coeffs_42 + int(11)) + make_float3 (fTmp0C_12 * x_33) * *(coeffs_42 + int(12)) + make_float3 (fTmp1B_12 * fC1_12) * *(coeffs_42 + int(13)) + make_float3 (-0.59004360437393188f * (x_33 * fC1_12 - y_31 * fS1_12)) * *(coeffs_42 + int(14)));
}

inline __device__ void sh3_to_color_dir_vjp_inplace(float3  dir_13, float3  coeff_dc_43, float3  * coeffs_43, float3  v_colors_29, float3  * v_coeff_dc_29, float3  * v_coeffs_29, float3  * v_dir_9)
{
    *v_coeff_dc_29 = *v_coeff_dc_29 + make_float3 (0.282094806432724f) * v_colors_29;
    float _S2079 = dir_13.x;
    float _S2080 = dir_13.y;
    float _S2081 = dir_13.z;
    float inorm_6 = (F32_rsqrt((_S2079 * _S2079 + _S2080 * _S2080 + _S2081 * _S2081)));
    float x_34 = _S2079 * inorm_6;
    float y_32 = _S2080 * inorm_6;
    float z_31 = _S2081 * inorm_6;
    float3  * _S2082 = v_coeffs_29 + int(0);
    *_S2082 = *_S2082 + make_float3 (-0.48860251903533936f * y_32) * v_colors_29;
    float3  * _S2083 = v_coeffs_29 + int(1);
    *_S2083 = *_S2083 + make_float3 (0.48860251903533936f * z_31) * v_colors_29;
    float3  * _S2084 = v_coeffs_29 + int(2);
    *_S2084 = *_S2084 + make_float3 (-0.48860251903533936f * x_34) * v_colors_29;
    float _S2085 = -0.48860251903533936f * dot_0(*(coeffs_43 + int(2)), v_colors_29);
    float _S2086 = -0.48860251903533936f * dot_0(*(coeffs_43 + int(0)), v_colors_29);
    float _S2087 = 0.48860251903533936f * dot_0(*(coeffs_43 + int(1)), v_colors_29);
    float z2_13 = z_31 * z_31;
    float fTmp0B_23 = -1.09254848957061768f * z_31;
    float fC1_13 = x_34 * x_34 - y_32 * y_32;
    float _S2088 = 2.0f * x_34;
    float fS1_13 = _S2088 * y_32;
    float pSH6_19 = 0.94617468118667603f * z2_13 - 0.31539157032966614f;
    float pSH7_18 = fTmp0B_23 * x_34;
    float pSH5_18 = fTmp0B_23 * y_32;
    float pSH8_18 = 0.54627424478530884f * fC1_13;
    float3  * _S2089 = v_coeffs_29 + int(3);
    *_S2089 = *_S2089 + make_float3 (0.54627424478530884f * fS1_13) * v_colors_29;
    float3  * _S2090 = v_coeffs_29 + int(4);
    *_S2090 = *_S2090 + make_float3 (pSH5_18) * v_colors_29;
    float3  * _S2091 = v_coeffs_29 + int(5);
    *_S2091 = *_S2091 + make_float3 (pSH6_19) * v_colors_29;
    float3  * _S2092 = v_coeffs_29 + int(6);
    *_S2092 = *_S2092 + make_float3 (pSH7_18) * v_colors_29;
    float3  * _S2093 = v_coeffs_29 + int(7);
    *_S2093 = *_S2093 + make_float3 (pSH8_18) * v_colors_29;
    float fC1_y_8 = -2.0f * y_32;
    float fS1_x_8 = 2.0f * y_32;
    float pSH8_x_15 = 0.54627424478530884f * _S2088;
    float3  * _S2094 = coeffs_43 + int(3);
    float3  * _S2095 = coeffs_43 + int(7);
    float3  * _S2096 = coeffs_43 + int(6);
    float v_x_24 = _S2085 + dot_0(v_colors_29, make_float3 (0.54627424478530884f * fS1_x_8) * *_S2094 + make_float3 (pSH8_x_15) * *_S2095 + make_float3 (fTmp0B_23) * *_S2096);
    float3  * _S2097 = coeffs_43 + int(4);
    float v_y_24 = _S2086 + dot_0(v_colors_29, make_float3 (pSH8_x_15) * *_S2094 + make_float3 (0.54627424478530884f * fC1_y_8) * *_S2095 + make_float3 (fTmp0B_23) * *_S2097);
    float v_z_24 = _S2087 + dot_0(v_colors_29, make_float3 (1.89234936237335205f * z_31) * *(coeffs_43 + int(5)) + make_float3 (-1.09254848957061768f * x_34) * *_S2096 + make_float3 (-1.09254848957061768f * y_32) * *_S2097);
    float fTmp0C_13 = -2.28522896766662598f * z2_13 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_31;
    float pSH12_11 = z_31 * (1.86588168144226074f * z2_13 - 1.11952900886535645f);
    float pSH13_10 = fTmp0C_13 * x_34;
    float pSH11_10 = fTmp0C_13 * y_32;
    float pSH14_10 = fTmp1B_13 * fC1_13;
    float pSH10_10 = fTmp1B_13 * fS1_13;
    float pSH15_10 = -0.59004360437393188f * (x_34 * fC1_13 - y_32 * fS1_13);
    float3  * _S2098 = v_coeffs_29 + int(8);
    *_S2098 = *_S2098 + make_float3 (-0.59004360437393188f * (x_34 * fS1_13 + y_32 * fC1_13)) * v_colors_29;
    float3  * _S2099 = v_coeffs_29 + int(9);
    *_S2099 = *_S2099 + make_float3 (pSH10_10) * v_colors_29;
    float3  * _S2100 = v_coeffs_29 + int(10);
    *_S2100 = *_S2100 + make_float3 (pSH11_10) * v_colors_29;
    float3  * _S2101 = v_coeffs_29 + int(11);
    *_S2101 = *_S2101 + make_float3 (pSH12_11) * v_colors_29;
    float3  * _S2102 = v_coeffs_29 + int(12);
    *_S2102 = *_S2102 + make_float3 (pSH13_10) * v_colors_29;
    float3  * _S2103 = v_coeffs_29 + int(13);
    *_S2103 = *_S2103 + make_float3 (pSH14_10) * v_colors_29;
    float3  * _S2104 = v_coeffs_29 + int(14);
    *_S2104 = *_S2104 + make_float3 (pSH15_10) * v_colors_29;
    float fTmp0C_z_8 = -4.57045793533325195f * z_31;
    float _S2105 = x_34 * _S2088;
    float _S2106 = y_32 * _S2088;
    float pSH14_x_8 = fTmp1B_13 * _S2088;
    float3  * _S2107 = coeffs_43 + int(8);
    float3  * _S2108 = coeffs_43 + int(14);
    float3  * _S2109 = coeffs_43 + int(9);
    float3  * _S2110 = coeffs_43 + int(13);
    float3  * _S2111 = coeffs_43 + int(12);
    float3  * _S2112 = coeffs_43 + int(10);
    float3  dir_n_22 = make_float3 (x_34, y_32, z_31);
    float3  v_dir_n_38 = make_float3 (v_x_24 + dot_0(v_colors_29, make_float3 (-0.59004360437393188f * (fS1_13 + x_34 * fS1_x_8 + _S2106)) * *_S2107 + make_float3 (-0.59004360437393188f * (fC1_13 + _S2105 - y_32 * fS1_x_8)) * *_S2108 + make_float3 (fTmp1B_13 * fS1_x_8) * *_S2109 + make_float3 (pSH14_x_8) * *_S2110 + make_float3 (fTmp0C_13) * *_S2111), v_y_24 + dot_0(v_colors_29, make_float3 (-0.59004360437393188f * (_S2105 + fC1_13 + y_32 * fC1_y_8)) * *_S2107 + make_float3 (-0.59004360437393188f * (x_34 * fC1_y_8 - fS1_13 - _S2106)) * *_S2108 + make_float3 (pSH14_x_8) * *_S2109 + make_float3 (fTmp1B_13 * fC1_y_8) * *_S2110 + make_float3 (fTmp0C_13) * *_S2112), v_z_24 + dot_0(v_colors_29, make_float3 (5.59764480590820312f * z2_13 - 1.11952900886535645f) * *(coeffs_43 + int(11)) + make_float3 (fTmp0C_z_8 * x_34) * *_S2111 + make_float3 (fTmp0C_z_8 * y_32) * *_S2112 + make_float3 (1.44530570507049561f * fC1_13) * *_S2110 + make_float3 (1.44530570507049561f * fS1_13) * *_S2109));
    *v_dir_9 = *v_dir_9 + (v_dir_n_38 - make_float3 (dot_0(v_dir_n_38, dir_n_22)) * dir_n_22) * make_float3 (inorm_6);
    return;
}

inline __device__ void sh3_to_color_dir_vjp_atomic(float3  dir_14, float3  coeff_dc_44, float3  * coeffs_44, float3  v_colors_30, float3  * v_coeff_dc_30, float3  * v_coeffs_30, float3  * v_dir_10)
{
    *v_coeff_dc_30 = *v_coeff_dc_30 + make_float3 (0.282094806432724f) * v_colors_30;
    float _S2113 = dir_14.x;
    float _S2114 = dir_14.y;
    float _S2115 = dir_14.z;
    float inorm_7 = (F32_rsqrt((_S2113 * _S2113 + _S2114 * _S2114 + _S2115 * _S2115)));
    float x_35 = _S2113 * inorm_7;
    float y_33 = _S2114 * inorm_7;
    float z_32 = _S2115 * inorm_7;
    float3  temp_211 = make_float3 (-0.48860251903533936f * y_33) * v_colors_30;
    float _S2116 = dot_0(temp_211, temp_211);
    bool _S2117;
    if((F32_isfinite((_S2116))))
    {
        _S2117 = _S2116 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2118 = v_coeffs_30 + int(0);
        float _S2119 = atomicAdd(&(_S2118->x), temp_211.x);
        float _S2120 = atomicAdd(&(_S2118->y), temp_211.y);
        float _S2121 = atomicAdd(&(_S2118->z), temp_211.z);
    }
    float3  temp_212 = make_float3 (0.48860251903533936f * z_32) * v_colors_30;
    float _S2122 = dot_0(temp_212, temp_212);
    if((F32_isfinite((_S2122))))
    {
        _S2117 = _S2122 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2123 = v_coeffs_30 + int(1);
        float _S2124 = atomicAdd(&(_S2123->x), temp_212.x);
        float _S2125 = atomicAdd(&(_S2123->y), temp_212.y);
        float _S2126 = atomicAdd(&(_S2123->z), temp_212.z);
    }
    float3  temp_213 = make_float3 (-0.48860251903533936f * x_35) * v_colors_30;
    float _S2127 = dot_0(temp_213, temp_213);
    if((F32_isfinite((_S2127))))
    {
        _S2117 = _S2127 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2128 = v_coeffs_30 + int(2);
        float _S2129 = atomicAdd(&(_S2128->x), temp_213.x);
        float _S2130 = atomicAdd(&(_S2128->y), temp_213.y);
        float _S2131 = atomicAdd(&(_S2128->z), temp_213.z);
    }
    float _S2132 = -0.48860251903533936f * dot_0(*(coeffs_44 + int(2)), v_colors_30);
    float _S2133 = -0.48860251903533936f * dot_0(*(coeffs_44 + int(0)), v_colors_30);
    float _S2134 = 0.48860251903533936f * dot_0(*(coeffs_44 + int(1)), v_colors_30);
    float z2_14 = z_32 * z_32;
    float fTmp0B_24 = -1.09254848957061768f * z_32;
    float fC1_14 = x_35 * x_35 - y_33 * y_33;
    float _S2135 = 2.0f * x_35;
    float fS1_14 = _S2135 * y_33;
    float pSH6_20 = 0.94617468118667603f * z2_14 - 0.31539157032966614f;
    float pSH7_19 = fTmp0B_24 * x_35;
    float pSH5_19 = fTmp0B_24 * y_33;
    float pSH8_19 = 0.54627424478530884f * fC1_14;
    float3  temp_214 = make_float3 (0.54627424478530884f * fS1_14) * v_colors_30;
    float _S2136 = dot_0(temp_214, temp_214);
    if((F32_isfinite((_S2136))))
    {
        _S2117 = _S2136 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2137 = v_coeffs_30 + int(3);
        float _S2138 = atomicAdd(&(_S2137->x), temp_214.x);
        float _S2139 = atomicAdd(&(_S2137->y), temp_214.y);
        float _S2140 = atomicAdd(&(_S2137->z), temp_214.z);
    }
    float3  temp_215 = make_float3 (pSH5_19) * v_colors_30;
    float _S2141 = dot_0(temp_215, temp_215);
    if((F32_isfinite((_S2141))))
    {
        _S2117 = _S2141 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2142 = v_coeffs_30 + int(4);
        float _S2143 = atomicAdd(&(_S2142->x), temp_215.x);
        float _S2144 = atomicAdd(&(_S2142->y), temp_215.y);
        float _S2145 = atomicAdd(&(_S2142->z), temp_215.z);
    }
    float3  temp_216 = make_float3 (pSH6_20) * v_colors_30;
    float _S2146 = dot_0(temp_216, temp_216);
    if((F32_isfinite((_S2146))))
    {
        _S2117 = _S2146 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2147 = v_coeffs_30 + int(5);
        float _S2148 = atomicAdd(&(_S2147->x), temp_216.x);
        float _S2149 = atomicAdd(&(_S2147->y), temp_216.y);
        float _S2150 = atomicAdd(&(_S2147->z), temp_216.z);
    }
    float3  temp_217 = make_float3 (pSH7_19) * v_colors_30;
    float _S2151 = dot_0(temp_217, temp_217);
    if((F32_isfinite((_S2151))))
    {
        _S2117 = _S2151 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2152 = v_coeffs_30 + int(6);
        float _S2153 = atomicAdd(&(_S2152->x), temp_217.x);
        float _S2154 = atomicAdd(&(_S2152->y), temp_217.y);
        float _S2155 = atomicAdd(&(_S2152->z), temp_217.z);
    }
    float3  temp_218 = make_float3 (pSH8_19) * v_colors_30;
    float _S2156 = dot_0(temp_218, temp_218);
    if((F32_isfinite((_S2156))))
    {
        _S2117 = _S2156 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2157 = v_coeffs_30 + int(7);
        float _S2158 = atomicAdd(&(_S2157->x), temp_218.x);
        float _S2159 = atomicAdd(&(_S2157->y), temp_218.y);
        float _S2160 = atomicAdd(&(_S2157->z), temp_218.z);
    }
    float fC1_y_9 = -2.0f * y_33;
    float fS1_x_9 = 2.0f * y_33;
    float pSH8_x_16 = 0.54627424478530884f * _S2135;
    float3  * _S2161 = coeffs_44 + int(3);
    float3  * _S2162 = coeffs_44 + int(7);
    float3  * _S2163 = coeffs_44 + int(6);
    float v_x_25 = _S2132 + dot_0(v_colors_30, make_float3 (0.54627424478530884f * fS1_x_9) * *_S2161 + make_float3 (pSH8_x_16) * *_S2162 + make_float3 (fTmp0B_24) * *_S2163);
    float3  * _S2164 = coeffs_44 + int(4);
    float v_y_25 = _S2133 + dot_0(v_colors_30, make_float3 (pSH8_x_16) * *_S2161 + make_float3 (0.54627424478530884f * fC1_y_9) * *_S2162 + make_float3 (fTmp0B_24) * *_S2164);
    float v_z_25 = _S2134 + dot_0(v_colors_30, make_float3 (1.89234936237335205f * z_32) * *(coeffs_44 + int(5)) + make_float3 (-1.09254848957061768f * x_35) * *_S2163 + make_float3 (-1.09254848957061768f * y_33) * *_S2164);
    float fTmp0C_14 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_32;
    float pSH12_12 = z_32 * (1.86588168144226074f * z2_14 - 1.11952900886535645f);
    float pSH13_11 = fTmp0C_14 * x_35;
    float pSH11_11 = fTmp0C_14 * y_33;
    float pSH14_11 = fTmp1B_14 * fC1_14;
    float pSH10_11 = fTmp1B_14 * fS1_14;
    float pSH15_11 = -0.59004360437393188f * (x_35 * fC1_14 - y_33 * fS1_14);
    float3  temp_219 = make_float3 (-0.59004360437393188f * (x_35 * fS1_14 + y_33 * fC1_14)) * v_colors_30;
    float _S2165 = dot_0(temp_219, temp_219);
    if((F32_isfinite((_S2165))))
    {
        _S2117 = _S2165 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2166 = v_coeffs_30 + int(8);
        float _S2167 = atomicAdd(&(_S2166->x), temp_219.x);
        float _S2168 = atomicAdd(&(_S2166->y), temp_219.y);
        float _S2169 = atomicAdd(&(_S2166->z), temp_219.z);
    }
    float3  temp_220 = make_float3 (pSH10_11) * v_colors_30;
    float _S2170 = dot_0(temp_220, temp_220);
    if((F32_isfinite((_S2170))))
    {
        _S2117 = _S2170 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2171 = v_coeffs_30 + int(9);
        float _S2172 = atomicAdd(&(_S2171->x), temp_220.x);
        float _S2173 = atomicAdd(&(_S2171->y), temp_220.y);
        float _S2174 = atomicAdd(&(_S2171->z), temp_220.z);
    }
    float3  temp_221 = make_float3 (pSH11_11) * v_colors_30;
    float _S2175 = dot_0(temp_221, temp_221);
    if((F32_isfinite((_S2175))))
    {
        _S2117 = _S2175 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2176 = v_coeffs_30 + int(10);
        float _S2177 = atomicAdd(&(_S2176->x), temp_221.x);
        float _S2178 = atomicAdd(&(_S2176->y), temp_221.y);
        float _S2179 = atomicAdd(&(_S2176->z), temp_221.z);
    }
    float3  temp_222 = make_float3 (pSH12_12) * v_colors_30;
    float _S2180 = dot_0(temp_222, temp_222);
    if((F32_isfinite((_S2180))))
    {
        _S2117 = _S2180 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2181 = v_coeffs_30 + int(11);
        float _S2182 = atomicAdd(&(_S2181->x), temp_222.x);
        float _S2183 = atomicAdd(&(_S2181->y), temp_222.y);
        float _S2184 = atomicAdd(&(_S2181->z), temp_222.z);
    }
    float3  temp_223 = make_float3 (pSH13_11) * v_colors_30;
    float _S2185 = dot_0(temp_223, temp_223);
    if((F32_isfinite((_S2185))))
    {
        _S2117 = _S2185 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2186 = v_coeffs_30 + int(12);
        float _S2187 = atomicAdd(&(_S2186->x), temp_223.x);
        float _S2188 = atomicAdd(&(_S2186->y), temp_223.y);
        float _S2189 = atomicAdd(&(_S2186->z), temp_223.z);
    }
    float3  temp_224 = make_float3 (pSH14_11) * v_colors_30;
    float _S2190 = dot_0(temp_224, temp_224);
    if((F32_isfinite((_S2190))))
    {
        _S2117 = _S2190 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2191 = v_coeffs_30 + int(13);
        float _S2192 = atomicAdd(&(_S2191->x), temp_224.x);
        float _S2193 = atomicAdd(&(_S2191->y), temp_224.y);
        float _S2194 = atomicAdd(&(_S2191->z), temp_224.z);
    }
    float3  temp_225 = make_float3 (pSH15_11) * v_colors_30;
    float _S2195 = dot_0(temp_225, temp_225);
    if((F32_isfinite((_S2195))))
    {
        _S2117 = _S2195 != 0.0f;
    }
    else
    {
        _S2117 = false;
    }
    if(_S2117)
    {
        float3  * _S2196 = v_coeffs_30 + int(14);
        float _S2197 = atomicAdd(&(_S2196->x), temp_225.x);
        float _S2198 = atomicAdd(&(_S2196->y), temp_225.y);
        float _S2199 = atomicAdd(&(_S2196->z), temp_225.z);
    }
    float fTmp0C_z_9 = -4.57045793533325195f * z_32;
    float _S2200 = x_35 * _S2135;
    float _S2201 = y_33 * _S2135;
    float pSH14_x_9 = fTmp1B_14 * _S2135;
    float3  * _S2202 = coeffs_44 + int(8);
    float3  * _S2203 = coeffs_44 + int(14);
    float3  * _S2204 = coeffs_44 + int(9);
    float3  * _S2205 = coeffs_44 + int(13);
    float3  * _S2206 = coeffs_44 + int(12);
    float3  * _S2207 = coeffs_44 + int(10);
    float3  dir_n_23 = make_float3 (x_35, y_33, z_32);
    float3  v_dir_n_39 = make_float3 (v_x_25 + dot_0(v_colors_30, make_float3 (-0.59004360437393188f * (fS1_14 + x_35 * fS1_x_9 + _S2201)) * *_S2202 + make_float3 (-0.59004360437393188f * (fC1_14 + _S2200 - y_33 * fS1_x_9)) * *_S2203 + make_float3 (fTmp1B_14 * fS1_x_9) * *_S2204 + make_float3 (pSH14_x_9) * *_S2205 + make_float3 (fTmp0C_14) * *_S2206), v_y_25 + dot_0(v_colors_30, make_float3 (-0.59004360437393188f * (_S2200 + fC1_14 + y_33 * fC1_y_9)) * *_S2202 + make_float3 (-0.59004360437393188f * (x_35 * fC1_y_9 - fS1_14 - _S2201)) * *_S2203 + make_float3 (pSH14_x_9) * *_S2204 + make_float3 (fTmp1B_14 * fC1_y_9) * *_S2205 + make_float3 (fTmp0C_14) * *_S2207), v_z_25 + dot_0(v_colors_30, make_float3 (5.59764480590820312f * z2_14 - 1.11952900886535645f) * *(coeffs_44 + int(11)) + make_float3 (fTmp0C_z_9 * x_35) * *_S2206 + make_float3 (fTmp0C_z_9 * y_33) * *_S2207 + make_float3 (1.44530570507049561f * fC1_14) * *_S2205 + make_float3 (1.44530570507049561f * fS1_14) * *_S2204));
    *v_dir_10 = *v_dir_10 + (v_dir_n_39 - make_float3 (dot_0(v_dir_n_39, dir_n_23)) * dir_n_23) * make_float3 (inorm_7);
    return;
}

inline __device__ void sh3_to_color_dir_vjp_block_atomic(float3  dir_15, float3  coeff_dc_45, float3  * coeffs_45, float3  v_colors_31, float3  * v_coeff_dc_31, float3  * v_coeffs_31, float3  * v_dir_11, uint thread_id_3, uint _S2208)
{
    *v_coeff_dc_31 = *v_coeff_dc_31 + make_float3 (0.282094806432724f) * v_colors_31;
    float _S2209 = dir_15.x;
    float _S2210 = dir_15.y;
    float _S2211 = dir_15.z;
    float inorm_8 = (F32_rsqrt((_S2209 * _S2209 + _S2210 * _S2210 + _S2211 * _S2211)));
    float x_36 = _S2209 * inorm_8;
    float y_34 = _S2210 * inorm_8;
    float z_33 = _S2211 * inorm_8;
    float3  _S2212 = make_float3 (-0.48860251903533936f * y_34) * v_colors_31;
    float3  _S2213 = _S2212;
    bool _S2214 = (F32_isfinite((_S2212.x)));
    uint _S2215 = __ballot_sync(_S2208, _S2214);
    float v_2;
    uint _S2216;
    if(_S2214)
    {
        float _S2217 = _S2213.x;
        uint _S2218 = __ballot_sync(_S2208, true);
        v_2 = _S2217;
        _S2216 = _S2218;
    }
    else
    {
        uint _S2219 = __ballot_sync(_S2208, true);
        v_2 = 0.0f;
        _S2216 = _S2219;
    }
    *&((&_S2213)->x) = v_2;
    bool _S2220 = (F32_isfinite((_S2213.y)));
    uint _S2221 = __ballot_sync(_S2216, _S2220);
    if(_S2220)
    {
        float _S2222 = _S2213.y;
        uint _S2223 = __ballot_sync(_S2216, true);
        v_2 = _S2222;
        _S2216 = _S2223;
    }
    else
    {
        uint _S2224 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2224;
    }
    *&((&_S2213)->y) = v_2;
    bool _S2225 = (F32_isfinite((_S2213.z)));
    uint _S2226 = __ballot_sync(_S2216, _S2225);
    if(_S2225)
    {
        float _S2227 = _S2213.z;
        uint _S2228 = __ballot_sync(_S2216, true);
        v_2 = _S2227;
        _S2216 = _S2228;
    }
    else
    {
        uint _S2229 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2229;
    }
    *&((&_S2213)->z) = v_2;
    float _S2230 = WaveActiveSum_0(_S2213.x, _S2216);
    *&((&_S2213)->x) = _S2230;
    float _S2231 = WaveActiveSum_0(_S2213.y, _S2216);
    *&((&_S2213)->y) = _S2231;
    float _S2232 = WaveActiveSum_0(_S2213.z, _S2216);
    *&((&_S2213)->z) = _S2232;
    uint warp_id_11 = thread_id_3 / 32U;
    uint lane_id_2 = thread_id_3 % 32U;
    bool _S2233 = lane_id_2 == 0U;
    uint _S2234 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_11] = _S2213.x;
        (*&_sh_block_reduce_shared_0)[warp_id_11 + 16U] = _S2213.y;
        (*&_sh_block_reduce_shared_0)[warp_id_11 + 32U] = _S2213.z;
        uint _S2235 = __ballot_sync(_S2216, true);
        _S2216 = _S2235;
    }
    else
    {
        uint _S2236 = __ballot_sync(_S2216, true);
        _S2216 = _S2236;
    }
    __syncthreads();
    bool _S2237 = warp_id_11 < 3U;
    uint _S2238 = __ballot_sync(_S2216, _S2237);
    uint _S2239;
    bool _S2240;
    if(_S2237)
    {
        bool _S2241 = lane_id_2 < 16U;
        uint _S2242 = __ballot_sync(_S2238, _S2241);
        if(_S2241)
        {
            float _S2243 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_11 * 16U];
            uint _S2244 = __ballot_sync(_S2238, true);
            v_2 = _S2243;
            _S2239 = _S2244;
        }
        else
        {
            uint _S2245 = __ballot_sync(_S2238, true);
            v_2 = 0.0f;
            _S2239 = _S2245;
        }
        float _S2246 = WaveActiveSum_0(v_2, _S2239);
        uint _S2247 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2248 = _S2246 != 0.0f;
            uint _S2249 = __ballot_sync(_S2239, true);
            _S2240 = _S2248;
            _S2239 = _S2249;
        }
        else
        {
            uint _S2250 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2250;
        }
        uint _S2251 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2252 = warp_id_11 == 0U;
            uint _S2253 = __ballot_sync(_S2251, _S2252);
            if(_S2252)
            {
                float _S2254 = atomicAdd(&((v_coeffs_31 + 0U)->x), _S2246);
                uint _S2255 = __ballot_sync(_S2251, true);
            }
            else
            {
                uint _S2256 = _S2251 & (~_S2253);
                bool _S2257 = warp_id_11 == 1U;
                uint _S2258 = __ballot_sync(_S2256, _S2257);
                if(_S2257)
                {
                    float _S2259 = atomicAdd(&((v_coeffs_31 + 0U)->y), _S2246);
                    uint _S2260 = __ballot_sync(_S2256, true);
                }
                else
                {
                    float _S2261 = atomicAdd(&((v_coeffs_31 + 0U)->z), _S2246);
                    uint _S2262 = __ballot_sync(_S2256, true);
                }
                uint _S2263 = __ballot_sync(_S2251, true);
            }
            uint _S2264 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2265 = __ballot_sync(_S2239, true);
        }
        uint _S2266 = __ballot_sync(_S2216, true);
        _S2216 = _S2266;
    }
    else
    {
        uint _S2267 = __ballot_sync(_S2216, true);
        _S2216 = _S2267;
    }
    float3  _S2268 = make_float3 (0.48860251903533936f * z_33) * v_colors_31;
    float3  _S2269 = _S2268;
    bool _S2270 = (F32_isfinite((_S2268.x)));
    uint _S2271 = __ballot_sync(_S2216, _S2270);
    if(_S2270)
    {
        float _S2272 = _S2269.x;
        uint _S2273 = __ballot_sync(_S2216, true);
        v_2 = _S2272;
        _S2216 = _S2273;
    }
    else
    {
        uint _S2274 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2274;
    }
    *&((&_S2269)->x) = v_2;
    bool _S2275 = (F32_isfinite((_S2269.y)));
    uint _S2276 = __ballot_sync(_S2216, _S2275);
    if(_S2275)
    {
        float _S2277 = _S2269.y;
        uint _S2278 = __ballot_sync(_S2216, true);
        v_2 = _S2277;
        _S2216 = _S2278;
    }
    else
    {
        uint _S2279 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2279;
    }
    *&((&_S2269)->y) = v_2;
    bool _S2280 = (F32_isfinite((_S2269.z)));
    uint _S2281 = __ballot_sync(_S2216, _S2280);
    if(_S2280)
    {
        float _S2282 = _S2269.z;
        uint _S2283 = __ballot_sync(_S2216, true);
        v_2 = _S2282;
        _S2216 = _S2283;
    }
    else
    {
        uint _S2284 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2284;
    }
    *&((&_S2269)->z) = v_2;
    float _S2285 = WaveActiveSum_0(_S2269.x, _S2216);
    *&((&_S2269)->x) = _S2285;
    float _S2286 = WaveActiveSum_0(_S2269.y, _S2216);
    *&((&_S2269)->y) = _S2286;
    float _S2287 = WaveActiveSum_0(_S2269.z, _S2216);
    *&((&_S2269)->z) = _S2287;
    uint warp_id_12 = thread_id_3 / 32U;
    uint _S2288 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_12] = _S2269.x;
        (*&_sh_block_reduce_shared_0)[warp_id_12 + 16U] = _S2269.y;
        (*&_sh_block_reduce_shared_0)[warp_id_12 + 32U] = _S2269.z;
        uint _S2289 = __ballot_sync(_S2216, true);
        _S2216 = _S2289;
    }
    else
    {
        uint _S2290 = __ballot_sync(_S2216, true);
        _S2216 = _S2290;
    }
    __syncthreads();
    bool _S2291 = warp_id_12 < 3U;
    uint _S2292 = __ballot_sync(_S2216, _S2291);
    if(_S2291)
    {
        bool _S2293 = lane_id_2 < 16U;
        uint _S2294 = __ballot_sync(_S2292, _S2293);
        if(_S2293)
        {
            float _S2295 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_12 * 16U];
            uint _S2296 = __ballot_sync(_S2292, true);
            v_2 = _S2295;
            _S2239 = _S2296;
        }
        else
        {
            uint _S2297 = __ballot_sync(_S2292, true);
            v_2 = 0.0f;
            _S2239 = _S2297;
        }
        float _S2298 = WaveActiveSum_0(v_2, _S2239);
        uint _S2299 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2300 = _S2298 != 0.0f;
            uint _S2301 = __ballot_sync(_S2239, true);
            _S2240 = _S2300;
            _S2239 = _S2301;
        }
        else
        {
            uint _S2302 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2302;
        }
        uint _S2303 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2304 = warp_id_12 == 0U;
            uint _S2305 = __ballot_sync(_S2303, _S2304);
            if(_S2304)
            {
                float _S2306 = atomicAdd(&((v_coeffs_31 + 1U)->x), _S2298);
                uint _S2307 = __ballot_sync(_S2303, true);
            }
            else
            {
                uint _S2308 = _S2303 & (~_S2305);
                bool _S2309 = warp_id_12 == 1U;
                uint _S2310 = __ballot_sync(_S2308, _S2309);
                if(_S2309)
                {
                    float _S2311 = atomicAdd(&((v_coeffs_31 + 1U)->y), _S2298);
                    uint _S2312 = __ballot_sync(_S2308, true);
                }
                else
                {
                    float _S2313 = atomicAdd(&((v_coeffs_31 + 1U)->z), _S2298);
                    uint _S2314 = __ballot_sync(_S2308, true);
                }
                uint _S2315 = __ballot_sync(_S2303, true);
            }
            uint _S2316 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2317 = __ballot_sync(_S2239, true);
        }
        uint _S2318 = __ballot_sync(_S2216, true);
        _S2216 = _S2318;
    }
    else
    {
        uint _S2319 = __ballot_sync(_S2216, true);
        _S2216 = _S2319;
    }
    float3  _S2320 = make_float3 (-0.48860251903533936f * x_36) * v_colors_31;
    float3  _S2321 = _S2320;
    bool _S2322 = (F32_isfinite((_S2320.x)));
    uint _S2323 = __ballot_sync(_S2216, _S2322);
    if(_S2322)
    {
        float _S2324 = _S2321.x;
        uint _S2325 = __ballot_sync(_S2216, true);
        v_2 = _S2324;
        _S2216 = _S2325;
    }
    else
    {
        uint _S2326 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2326;
    }
    *&((&_S2321)->x) = v_2;
    bool _S2327 = (F32_isfinite((_S2321.y)));
    uint _S2328 = __ballot_sync(_S2216, _S2327);
    if(_S2327)
    {
        float _S2329 = _S2321.y;
        uint _S2330 = __ballot_sync(_S2216, true);
        v_2 = _S2329;
        _S2216 = _S2330;
    }
    else
    {
        uint _S2331 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2331;
    }
    *&((&_S2321)->y) = v_2;
    bool _S2332 = (F32_isfinite((_S2321.z)));
    uint _S2333 = __ballot_sync(_S2216, _S2332);
    if(_S2332)
    {
        float _S2334 = _S2321.z;
        uint _S2335 = __ballot_sync(_S2216, true);
        v_2 = _S2334;
        _S2216 = _S2335;
    }
    else
    {
        uint _S2336 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2336;
    }
    *&((&_S2321)->z) = v_2;
    float _S2337 = WaveActiveSum_0(_S2321.x, _S2216);
    *&((&_S2321)->x) = _S2337;
    float _S2338 = WaveActiveSum_0(_S2321.y, _S2216);
    *&((&_S2321)->y) = _S2338;
    float _S2339 = WaveActiveSum_0(_S2321.z, _S2216);
    *&((&_S2321)->z) = _S2339;
    uint warp_id_13 = thread_id_3 / 32U;
    uint _S2340 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_13] = _S2321.x;
        (*&_sh_block_reduce_shared_0)[warp_id_13 + 16U] = _S2321.y;
        (*&_sh_block_reduce_shared_0)[warp_id_13 + 32U] = _S2321.z;
        uint _S2341 = __ballot_sync(_S2216, true);
        _S2216 = _S2341;
    }
    else
    {
        uint _S2342 = __ballot_sync(_S2216, true);
        _S2216 = _S2342;
    }
    __syncthreads();
    bool _S2343 = warp_id_13 < 3U;
    uint _S2344 = __ballot_sync(_S2216, _S2343);
    if(_S2343)
    {
        bool _S2345 = lane_id_2 < 16U;
        uint _S2346 = __ballot_sync(_S2344, _S2345);
        if(_S2345)
        {
            float _S2347 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_13 * 16U];
            uint _S2348 = __ballot_sync(_S2344, true);
            v_2 = _S2347;
            _S2239 = _S2348;
        }
        else
        {
            uint _S2349 = __ballot_sync(_S2344, true);
            v_2 = 0.0f;
            _S2239 = _S2349;
        }
        float _S2350 = WaveActiveSum_0(v_2, _S2239);
        uint _S2351 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2352 = _S2350 != 0.0f;
            uint _S2353 = __ballot_sync(_S2239, true);
            _S2240 = _S2352;
            _S2239 = _S2353;
        }
        else
        {
            uint _S2354 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2354;
        }
        uint _S2355 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2356 = warp_id_13 == 0U;
            uint _S2357 = __ballot_sync(_S2355, _S2356);
            if(_S2356)
            {
                float _S2358 = atomicAdd(&((v_coeffs_31 + 2U)->x), _S2350);
                uint _S2359 = __ballot_sync(_S2355, true);
            }
            else
            {
                uint _S2360 = _S2355 & (~_S2357);
                bool _S2361 = warp_id_13 == 1U;
                uint _S2362 = __ballot_sync(_S2360, _S2361);
                if(_S2361)
                {
                    float _S2363 = atomicAdd(&((v_coeffs_31 + 2U)->y), _S2350);
                    uint _S2364 = __ballot_sync(_S2360, true);
                }
                else
                {
                    float _S2365 = atomicAdd(&((v_coeffs_31 + 2U)->z), _S2350);
                    uint _S2366 = __ballot_sync(_S2360, true);
                }
                uint _S2367 = __ballot_sync(_S2355, true);
            }
            uint _S2368 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2369 = __ballot_sync(_S2239, true);
        }
        uint _S2370 = __ballot_sync(_S2216, true);
        _S2216 = _S2370;
    }
    else
    {
        uint _S2371 = __ballot_sync(_S2216, true);
        _S2216 = _S2371;
    }
    float _S2372 = -0.48860251903533936f * dot_0(*(coeffs_45 + int(2)), v_colors_31);
    float _S2373 = -0.48860251903533936f * dot_0(*(coeffs_45 + int(0)), v_colors_31);
    float _S2374 = 0.48860251903533936f * dot_0(*(coeffs_45 + int(1)), v_colors_31);
    float z2_15 = z_33 * z_33;
    float fTmp0B_25 = -1.09254848957061768f * z_33;
    float fC1_15 = x_36 * x_36 - y_34 * y_34;
    float _S2375 = 2.0f * x_36;
    float fS1_15 = _S2375 * y_34;
    float pSH6_21 = 0.94617468118667603f * z2_15 - 0.31539157032966614f;
    float pSH7_20 = fTmp0B_25 * x_36;
    float pSH5_20 = fTmp0B_25 * y_34;
    float pSH8_20 = 0.54627424478530884f * fC1_15;
    float3  _S2376 = make_float3 (0.54627424478530884f * fS1_15) * v_colors_31;
    float3  _S2377 = _S2376;
    bool _S2378 = (F32_isfinite((_S2376.x)));
    uint _S2379 = __ballot_sync(_S2216, _S2378);
    if(_S2378)
    {
        float _S2380 = _S2377.x;
        uint _S2381 = __ballot_sync(_S2216, true);
        v_2 = _S2380;
        _S2216 = _S2381;
    }
    else
    {
        uint _S2382 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2382;
    }
    *&((&_S2377)->x) = v_2;
    bool _S2383 = (F32_isfinite((_S2377.y)));
    uint _S2384 = __ballot_sync(_S2216, _S2383);
    if(_S2383)
    {
        float _S2385 = _S2377.y;
        uint _S2386 = __ballot_sync(_S2216, true);
        v_2 = _S2385;
        _S2216 = _S2386;
    }
    else
    {
        uint _S2387 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2387;
    }
    *&((&_S2377)->y) = v_2;
    bool _S2388 = (F32_isfinite((_S2377.z)));
    uint _S2389 = __ballot_sync(_S2216, _S2388);
    if(_S2388)
    {
        float _S2390 = _S2377.z;
        uint _S2391 = __ballot_sync(_S2216, true);
        v_2 = _S2390;
        _S2216 = _S2391;
    }
    else
    {
        uint _S2392 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2392;
    }
    *&((&_S2377)->z) = v_2;
    float _S2393 = WaveActiveSum_0(_S2377.x, _S2216);
    *&((&_S2377)->x) = _S2393;
    float _S2394 = WaveActiveSum_0(_S2377.y, _S2216);
    *&((&_S2377)->y) = _S2394;
    float _S2395 = WaveActiveSum_0(_S2377.z, _S2216);
    *&((&_S2377)->z) = _S2395;
    uint warp_id_14 = thread_id_3 / 32U;
    uint _S2396 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_14] = _S2377.x;
        (*&_sh_block_reduce_shared_0)[warp_id_14 + 16U] = _S2377.y;
        (*&_sh_block_reduce_shared_0)[warp_id_14 + 32U] = _S2377.z;
        uint _S2397 = __ballot_sync(_S2216, true);
        _S2216 = _S2397;
    }
    else
    {
        uint _S2398 = __ballot_sync(_S2216, true);
        _S2216 = _S2398;
    }
    __syncthreads();
    bool _S2399 = warp_id_14 < 3U;
    uint _S2400 = __ballot_sync(_S2216, _S2399);
    if(_S2399)
    {
        bool _S2401 = lane_id_2 < 16U;
        uint _S2402 = __ballot_sync(_S2400, _S2401);
        if(_S2401)
        {
            float _S2403 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_14 * 16U];
            uint _S2404 = __ballot_sync(_S2400, true);
            v_2 = _S2403;
            _S2239 = _S2404;
        }
        else
        {
            uint _S2405 = __ballot_sync(_S2400, true);
            v_2 = 0.0f;
            _S2239 = _S2405;
        }
        float _S2406 = WaveActiveSum_0(v_2, _S2239);
        uint _S2407 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2408 = _S2406 != 0.0f;
            uint _S2409 = __ballot_sync(_S2239, true);
            _S2240 = _S2408;
            _S2239 = _S2409;
        }
        else
        {
            uint _S2410 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2410;
        }
        uint _S2411 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2412 = warp_id_14 == 0U;
            uint _S2413 = __ballot_sync(_S2411, _S2412);
            if(_S2412)
            {
                float _S2414 = atomicAdd(&((v_coeffs_31 + 3U)->x), _S2406);
                uint _S2415 = __ballot_sync(_S2411, true);
            }
            else
            {
                uint _S2416 = _S2411 & (~_S2413);
                bool _S2417 = warp_id_14 == 1U;
                uint _S2418 = __ballot_sync(_S2416, _S2417);
                if(_S2417)
                {
                    float _S2419 = atomicAdd(&((v_coeffs_31 + 3U)->y), _S2406);
                    uint _S2420 = __ballot_sync(_S2416, true);
                }
                else
                {
                    float _S2421 = atomicAdd(&((v_coeffs_31 + 3U)->z), _S2406);
                    uint _S2422 = __ballot_sync(_S2416, true);
                }
                uint _S2423 = __ballot_sync(_S2411, true);
            }
            uint _S2424 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2425 = __ballot_sync(_S2239, true);
        }
        uint _S2426 = __ballot_sync(_S2216, true);
        _S2216 = _S2426;
    }
    else
    {
        uint _S2427 = __ballot_sync(_S2216, true);
        _S2216 = _S2427;
    }
    float3  _S2428 = make_float3 (pSH5_20) * v_colors_31;
    float3  _S2429 = _S2428;
    bool _S2430 = (F32_isfinite((_S2428.x)));
    uint _S2431 = __ballot_sync(_S2216, _S2430);
    if(_S2430)
    {
        float _S2432 = _S2429.x;
        uint _S2433 = __ballot_sync(_S2216, true);
        v_2 = _S2432;
        _S2216 = _S2433;
    }
    else
    {
        uint _S2434 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2434;
    }
    *&((&_S2429)->x) = v_2;
    bool _S2435 = (F32_isfinite((_S2429.y)));
    uint _S2436 = __ballot_sync(_S2216, _S2435);
    if(_S2435)
    {
        float _S2437 = _S2429.y;
        uint _S2438 = __ballot_sync(_S2216, true);
        v_2 = _S2437;
        _S2216 = _S2438;
    }
    else
    {
        uint _S2439 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2439;
    }
    *&((&_S2429)->y) = v_2;
    bool _S2440 = (F32_isfinite((_S2429.z)));
    uint _S2441 = __ballot_sync(_S2216, _S2440);
    if(_S2440)
    {
        float _S2442 = _S2429.z;
        uint _S2443 = __ballot_sync(_S2216, true);
        v_2 = _S2442;
        _S2216 = _S2443;
    }
    else
    {
        uint _S2444 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2444;
    }
    *&((&_S2429)->z) = v_2;
    float _S2445 = WaveActiveSum_0(_S2429.x, _S2216);
    *&((&_S2429)->x) = _S2445;
    float _S2446 = WaveActiveSum_0(_S2429.y, _S2216);
    *&((&_S2429)->y) = _S2446;
    float _S2447 = WaveActiveSum_0(_S2429.z, _S2216);
    *&((&_S2429)->z) = _S2447;
    uint warp_id_15 = thread_id_3 / 32U;
    uint _S2448 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_15] = _S2429.x;
        (*&_sh_block_reduce_shared_0)[warp_id_15 + 16U] = _S2429.y;
        (*&_sh_block_reduce_shared_0)[warp_id_15 + 32U] = _S2429.z;
        uint _S2449 = __ballot_sync(_S2216, true);
        _S2216 = _S2449;
    }
    else
    {
        uint _S2450 = __ballot_sync(_S2216, true);
        _S2216 = _S2450;
    }
    __syncthreads();
    bool _S2451 = warp_id_15 < 3U;
    uint _S2452 = __ballot_sync(_S2216, _S2451);
    if(_S2451)
    {
        bool _S2453 = lane_id_2 < 16U;
        uint _S2454 = __ballot_sync(_S2452, _S2453);
        if(_S2453)
        {
            float _S2455 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_15 * 16U];
            uint _S2456 = __ballot_sync(_S2452, true);
            v_2 = _S2455;
            _S2239 = _S2456;
        }
        else
        {
            uint _S2457 = __ballot_sync(_S2452, true);
            v_2 = 0.0f;
            _S2239 = _S2457;
        }
        float _S2458 = WaveActiveSum_0(v_2, _S2239);
        uint _S2459 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2460 = _S2458 != 0.0f;
            uint _S2461 = __ballot_sync(_S2239, true);
            _S2240 = _S2460;
            _S2239 = _S2461;
        }
        else
        {
            uint _S2462 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2462;
        }
        uint _S2463 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2464 = warp_id_15 == 0U;
            uint _S2465 = __ballot_sync(_S2463, _S2464);
            if(_S2464)
            {
                float _S2466 = atomicAdd(&((v_coeffs_31 + 4U)->x), _S2458);
                uint _S2467 = __ballot_sync(_S2463, true);
            }
            else
            {
                uint _S2468 = _S2463 & (~_S2465);
                bool _S2469 = warp_id_15 == 1U;
                uint _S2470 = __ballot_sync(_S2468, _S2469);
                if(_S2469)
                {
                    float _S2471 = atomicAdd(&((v_coeffs_31 + 4U)->y), _S2458);
                    uint _S2472 = __ballot_sync(_S2468, true);
                }
                else
                {
                    float _S2473 = atomicAdd(&((v_coeffs_31 + 4U)->z), _S2458);
                    uint _S2474 = __ballot_sync(_S2468, true);
                }
                uint _S2475 = __ballot_sync(_S2463, true);
            }
            uint _S2476 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2477 = __ballot_sync(_S2239, true);
        }
        uint _S2478 = __ballot_sync(_S2216, true);
        _S2216 = _S2478;
    }
    else
    {
        uint _S2479 = __ballot_sync(_S2216, true);
        _S2216 = _S2479;
    }
    float3  _S2480 = make_float3 (pSH6_21) * v_colors_31;
    float3  _S2481 = _S2480;
    bool _S2482 = (F32_isfinite((_S2480.x)));
    uint _S2483 = __ballot_sync(_S2216, _S2482);
    if(_S2482)
    {
        float _S2484 = _S2481.x;
        uint _S2485 = __ballot_sync(_S2216, true);
        v_2 = _S2484;
        _S2216 = _S2485;
    }
    else
    {
        uint _S2486 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2486;
    }
    *&((&_S2481)->x) = v_2;
    bool _S2487 = (F32_isfinite((_S2481.y)));
    uint _S2488 = __ballot_sync(_S2216, _S2487);
    if(_S2487)
    {
        float _S2489 = _S2481.y;
        uint _S2490 = __ballot_sync(_S2216, true);
        v_2 = _S2489;
        _S2216 = _S2490;
    }
    else
    {
        uint _S2491 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2491;
    }
    *&((&_S2481)->y) = v_2;
    bool _S2492 = (F32_isfinite((_S2481.z)));
    uint _S2493 = __ballot_sync(_S2216, _S2492);
    if(_S2492)
    {
        float _S2494 = _S2481.z;
        uint _S2495 = __ballot_sync(_S2216, true);
        v_2 = _S2494;
        _S2216 = _S2495;
    }
    else
    {
        uint _S2496 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2496;
    }
    *&((&_S2481)->z) = v_2;
    float _S2497 = WaveActiveSum_0(_S2481.x, _S2216);
    *&((&_S2481)->x) = _S2497;
    float _S2498 = WaveActiveSum_0(_S2481.y, _S2216);
    *&((&_S2481)->y) = _S2498;
    float _S2499 = WaveActiveSum_0(_S2481.z, _S2216);
    *&((&_S2481)->z) = _S2499;
    uint warp_id_16 = thread_id_3 / 32U;
    uint _S2500 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_16] = _S2481.x;
        (*&_sh_block_reduce_shared_0)[warp_id_16 + 16U] = _S2481.y;
        (*&_sh_block_reduce_shared_0)[warp_id_16 + 32U] = _S2481.z;
        uint _S2501 = __ballot_sync(_S2216, true);
        _S2216 = _S2501;
    }
    else
    {
        uint _S2502 = __ballot_sync(_S2216, true);
        _S2216 = _S2502;
    }
    __syncthreads();
    bool _S2503 = warp_id_16 < 3U;
    uint _S2504 = __ballot_sync(_S2216, _S2503);
    if(_S2503)
    {
        bool _S2505 = lane_id_2 < 16U;
        uint _S2506 = __ballot_sync(_S2504, _S2505);
        if(_S2505)
        {
            float _S2507 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_16 * 16U];
            uint _S2508 = __ballot_sync(_S2504, true);
            v_2 = _S2507;
            _S2239 = _S2508;
        }
        else
        {
            uint _S2509 = __ballot_sync(_S2504, true);
            v_2 = 0.0f;
            _S2239 = _S2509;
        }
        float _S2510 = WaveActiveSum_0(v_2, _S2239);
        uint _S2511 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2512 = _S2510 != 0.0f;
            uint _S2513 = __ballot_sync(_S2239, true);
            _S2240 = _S2512;
            _S2239 = _S2513;
        }
        else
        {
            uint _S2514 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2514;
        }
        uint _S2515 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2516 = warp_id_16 == 0U;
            uint _S2517 = __ballot_sync(_S2515, _S2516);
            if(_S2516)
            {
                float _S2518 = atomicAdd(&((v_coeffs_31 + 5U)->x), _S2510);
                uint _S2519 = __ballot_sync(_S2515, true);
            }
            else
            {
                uint _S2520 = _S2515 & (~_S2517);
                bool _S2521 = warp_id_16 == 1U;
                uint _S2522 = __ballot_sync(_S2520, _S2521);
                if(_S2521)
                {
                    float _S2523 = atomicAdd(&((v_coeffs_31 + 5U)->y), _S2510);
                    uint _S2524 = __ballot_sync(_S2520, true);
                }
                else
                {
                    float _S2525 = atomicAdd(&((v_coeffs_31 + 5U)->z), _S2510);
                    uint _S2526 = __ballot_sync(_S2520, true);
                }
                uint _S2527 = __ballot_sync(_S2515, true);
            }
            uint _S2528 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2529 = __ballot_sync(_S2239, true);
        }
        uint _S2530 = __ballot_sync(_S2216, true);
        _S2216 = _S2530;
    }
    else
    {
        uint _S2531 = __ballot_sync(_S2216, true);
        _S2216 = _S2531;
    }
    float3  _S2532 = make_float3 (pSH7_20) * v_colors_31;
    float3  _S2533 = _S2532;
    bool _S2534 = (F32_isfinite((_S2532.x)));
    uint _S2535 = __ballot_sync(_S2216, _S2534);
    if(_S2534)
    {
        float _S2536 = _S2533.x;
        uint _S2537 = __ballot_sync(_S2216, true);
        v_2 = _S2536;
        _S2216 = _S2537;
    }
    else
    {
        uint _S2538 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2538;
    }
    *&((&_S2533)->x) = v_2;
    bool _S2539 = (F32_isfinite((_S2533.y)));
    uint _S2540 = __ballot_sync(_S2216, _S2539);
    if(_S2539)
    {
        float _S2541 = _S2533.y;
        uint _S2542 = __ballot_sync(_S2216, true);
        v_2 = _S2541;
        _S2216 = _S2542;
    }
    else
    {
        uint _S2543 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2543;
    }
    *&((&_S2533)->y) = v_2;
    bool _S2544 = (F32_isfinite((_S2533.z)));
    uint _S2545 = __ballot_sync(_S2216, _S2544);
    if(_S2544)
    {
        float _S2546 = _S2533.z;
        uint _S2547 = __ballot_sync(_S2216, true);
        v_2 = _S2546;
        _S2216 = _S2547;
    }
    else
    {
        uint _S2548 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2548;
    }
    *&((&_S2533)->z) = v_2;
    float _S2549 = WaveActiveSum_0(_S2533.x, _S2216);
    *&((&_S2533)->x) = _S2549;
    float _S2550 = WaveActiveSum_0(_S2533.y, _S2216);
    *&((&_S2533)->y) = _S2550;
    float _S2551 = WaveActiveSum_0(_S2533.z, _S2216);
    *&((&_S2533)->z) = _S2551;
    uint warp_id_17 = thread_id_3 / 32U;
    uint _S2552 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_17] = _S2533.x;
        (*&_sh_block_reduce_shared_0)[warp_id_17 + 16U] = _S2533.y;
        (*&_sh_block_reduce_shared_0)[warp_id_17 + 32U] = _S2533.z;
        uint _S2553 = __ballot_sync(_S2216, true);
        _S2216 = _S2553;
    }
    else
    {
        uint _S2554 = __ballot_sync(_S2216, true);
        _S2216 = _S2554;
    }
    __syncthreads();
    bool _S2555 = warp_id_17 < 3U;
    uint _S2556 = __ballot_sync(_S2216, _S2555);
    if(_S2555)
    {
        bool _S2557 = lane_id_2 < 16U;
        uint _S2558 = __ballot_sync(_S2556, _S2557);
        if(_S2557)
        {
            float _S2559 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_17 * 16U];
            uint _S2560 = __ballot_sync(_S2556, true);
            v_2 = _S2559;
            _S2239 = _S2560;
        }
        else
        {
            uint _S2561 = __ballot_sync(_S2556, true);
            v_2 = 0.0f;
            _S2239 = _S2561;
        }
        float _S2562 = WaveActiveSum_0(v_2, _S2239);
        uint _S2563 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2564 = _S2562 != 0.0f;
            uint _S2565 = __ballot_sync(_S2239, true);
            _S2240 = _S2564;
            _S2239 = _S2565;
        }
        else
        {
            uint _S2566 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2566;
        }
        uint _S2567 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2568 = warp_id_17 == 0U;
            uint _S2569 = __ballot_sync(_S2567, _S2568);
            if(_S2568)
            {
                float _S2570 = atomicAdd(&((v_coeffs_31 + 6U)->x), _S2562);
                uint _S2571 = __ballot_sync(_S2567, true);
            }
            else
            {
                uint _S2572 = _S2567 & (~_S2569);
                bool _S2573 = warp_id_17 == 1U;
                uint _S2574 = __ballot_sync(_S2572, _S2573);
                if(_S2573)
                {
                    float _S2575 = atomicAdd(&((v_coeffs_31 + 6U)->y), _S2562);
                    uint _S2576 = __ballot_sync(_S2572, true);
                }
                else
                {
                    float _S2577 = atomicAdd(&((v_coeffs_31 + 6U)->z), _S2562);
                    uint _S2578 = __ballot_sync(_S2572, true);
                }
                uint _S2579 = __ballot_sync(_S2567, true);
            }
            uint _S2580 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2581 = __ballot_sync(_S2239, true);
        }
        uint _S2582 = __ballot_sync(_S2216, true);
        _S2216 = _S2582;
    }
    else
    {
        uint _S2583 = __ballot_sync(_S2216, true);
        _S2216 = _S2583;
    }
    float3  _S2584 = make_float3 (pSH8_20) * v_colors_31;
    float3  _S2585 = _S2584;
    bool _S2586 = (F32_isfinite((_S2584.x)));
    uint _S2587 = __ballot_sync(_S2216, _S2586);
    if(_S2586)
    {
        float _S2588 = _S2585.x;
        uint _S2589 = __ballot_sync(_S2216, true);
        v_2 = _S2588;
        _S2216 = _S2589;
    }
    else
    {
        uint _S2590 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2590;
    }
    *&((&_S2585)->x) = v_2;
    bool _S2591 = (F32_isfinite((_S2585.y)));
    uint _S2592 = __ballot_sync(_S2216, _S2591);
    if(_S2591)
    {
        float _S2593 = _S2585.y;
        uint _S2594 = __ballot_sync(_S2216, true);
        v_2 = _S2593;
        _S2216 = _S2594;
    }
    else
    {
        uint _S2595 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2595;
    }
    *&((&_S2585)->y) = v_2;
    bool _S2596 = (F32_isfinite((_S2585.z)));
    uint _S2597 = __ballot_sync(_S2216, _S2596);
    if(_S2596)
    {
        float _S2598 = _S2585.z;
        uint _S2599 = __ballot_sync(_S2216, true);
        v_2 = _S2598;
        _S2216 = _S2599;
    }
    else
    {
        uint _S2600 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2600;
    }
    *&((&_S2585)->z) = v_2;
    float _S2601 = WaveActiveSum_0(_S2585.x, _S2216);
    *&((&_S2585)->x) = _S2601;
    float _S2602 = WaveActiveSum_0(_S2585.y, _S2216);
    *&((&_S2585)->y) = _S2602;
    float _S2603 = WaveActiveSum_0(_S2585.z, _S2216);
    *&((&_S2585)->z) = _S2603;
    uint warp_id_18 = thread_id_3 / 32U;
    uint _S2604 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_18] = _S2585.x;
        (*&_sh_block_reduce_shared_0)[warp_id_18 + 16U] = _S2585.y;
        (*&_sh_block_reduce_shared_0)[warp_id_18 + 32U] = _S2585.z;
        uint _S2605 = __ballot_sync(_S2216, true);
        _S2216 = _S2605;
    }
    else
    {
        uint _S2606 = __ballot_sync(_S2216, true);
        _S2216 = _S2606;
    }
    __syncthreads();
    bool _S2607 = warp_id_18 < 3U;
    uint _S2608 = __ballot_sync(_S2216, _S2607);
    if(_S2607)
    {
        bool _S2609 = lane_id_2 < 16U;
        uint _S2610 = __ballot_sync(_S2608, _S2609);
        if(_S2609)
        {
            float _S2611 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_18 * 16U];
            uint _S2612 = __ballot_sync(_S2608, true);
            v_2 = _S2611;
            _S2239 = _S2612;
        }
        else
        {
            uint _S2613 = __ballot_sync(_S2608, true);
            v_2 = 0.0f;
            _S2239 = _S2613;
        }
        float _S2614 = WaveActiveSum_0(v_2, _S2239);
        uint _S2615 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2616 = _S2614 != 0.0f;
            uint _S2617 = __ballot_sync(_S2239, true);
            _S2240 = _S2616;
            _S2239 = _S2617;
        }
        else
        {
            uint _S2618 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2618;
        }
        uint _S2619 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2620 = warp_id_18 == 0U;
            uint _S2621 = __ballot_sync(_S2619, _S2620);
            if(_S2620)
            {
                float _S2622 = atomicAdd(&((v_coeffs_31 + 7U)->x), _S2614);
                uint _S2623 = __ballot_sync(_S2619, true);
            }
            else
            {
                uint _S2624 = _S2619 & (~_S2621);
                bool _S2625 = warp_id_18 == 1U;
                uint _S2626 = __ballot_sync(_S2624, _S2625);
                if(_S2625)
                {
                    float _S2627 = atomicAdd(&((v_coeffs_31 + 7U)->y), _S2614);
                    uint _S2628 = __ballot_sync(_S2624, true);
                }
                else
                {
                    float _S2629 = atomicAdd(&((v_coeffs_31 + 7U)->z), _S2614);
                    uint _S2630 = __ballot_sync(_S2624, true);
                }
                uint _S2631 = __ballot_sync(_S2619, true);
            }
            uint _S2632 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2633 = __ballot_sync(_S2239, true);
        }
        uint _S2634 = __ballot_sync(_S2216, true);
        _S2216 = _S2634;
    }
    else
    {
        uint _S2635 = __ballot_sync(_S2216, true);
        _S2216 = _S2635;
    }
    float fC1_y_10 = -2.0f * y_34;
    float fS1_x_10 = 2.0f * y_34;
    float pSH8_x_17 = 0.54627424478530884f * _S2375;
    float3  * _S2636 = coeffs_45 + int(3);
    float3  * _S2637 = coeffs_45 + int(7);
    float3  * _S2638 = coeffs_45 + int(6);
    float v_x_26 = _S2372 + dot_0(v_colors_31, make_float3 (0.54627424478530884f * fS1_x_10) * *_S2636 + make_float3 (pSH8_x_17) * *_S2637 + make_float3 (fTmp0B_25) * *_S2638);
    float3  * _S2639 = coeffs_45 + int(4);
    float v_y_26 = _S2373 + dot_0(v_colors_31, make_float3 (pSH8_x_17) * *_S2636 + make_float3 (0.54627424478530884f * fC1_y_10) * *_S2637 + make_float3 (fTmp0B_25) * *_S2639);
    float v_z_26 = _S2374 + dot_0(v_colors_31, make_float3 (1.89234936237335205f * z_33) * *(coeffs_45 + int(5)) + make_float3 (-1.09254848957061768f * x_36) * *_S2638 + make_float3 (-1.09254848957061768f * y_34) * *_S2639);
    float fTmp0C_15 = -2.28522896766662598f * z2_15 + 0.4570457935333252f;
    float fTmp1B_15 = 1.44530570507049561f * z_33;
    float pSH12_13 = z_33 * (1.86588168144226074f * z2_15 - 1.11952900886535645f);
    float pSH13_12 = fTmp0C_15 * x_36;
    float pSH11_12 = fTmp0C_15 * y_34;
    float pSH14_12 = fTmp1B_15 * fC1_15;
    float pSH10_12 = fTmp1B_15 * fS1_15;
    float pSH15_12 = -0.59004360437393188f * (x_36 * fC1_15 - y_34 * fS1_15);
    float3  _S2640 = make_float3 (-0.59004360437393188f * (x_36 * fS1_15 + y_34 * fC1_15)) * v_colors_31;
    float3  _S2641 = _S2640;
    bool _S2642 = (F32_isfinite((_S2640.x)));
    uint _S2643 = __ballot_sync(_S2216, _S2642);
    if(_S2642)
    {
        float _S2644 = _S2641.x;
        uint _S2645 = __ballot_sync(_S2216, true);
        v_2 = _S2644;
        _S2216 = _S2645;
    }
    else
    {
        uint _S2646 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2646;
    }
    *&((&_S2641)->x) = v_2;
    bool _S2647 = (F32_isfinite((_S2641.y)));
    uint _S2648 = __ballot_sync(_S2216, _S2647);
    if(_S2647)
    {
        float _S2649 = _S2641.y;
        uint _S2650 = __ballot_sync(_S2216, true);
        v_2 = _S2649;
        _S2216 = _S2650;
    }
    else
    {
        uint _S2651 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2651;
    }
    *&((&_S2641)->y) = v_2;
    bool _S2652 = (F32_isfinite((_S2641.z)));
    uint _S2653 = __ballot_sync(_S2216, _S2652);
    if(_S2652)
    {
        float _S2654 = _S2641.z;
        uint _S2655 = __ballot_sync(_S2216, true);
        v_2 = _S2654;
        _S2216 = _S2655;
    }
    else
    {
        uint _S2656 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2656;
    }
    *&((&_S2641)->z) = v_2;
    float _S2657 = WaveActiveSum_0(_S2641.x, _S2216);
    *&((&_S2641)->x) = _S2657;
    float _S2658 = WaveActiveSum_0(_S2641.y, _S2216);
    *&((&_S2641)->y) = _S2658;
    float _S2659 = WaveActiveSum_0(_S2641.z, _S2216);
    *&((&_S2641)->z) = _S2659;
    uint warp_id_19 = thread_id_3 / 32U;
    uint _S2660 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_19] = _S2641.x;
        (*&_sh_block_reduce_shared_0)[warp_id_19 + 16U] = _S2641.y;
        (*&_sh_block_reduce_shared_0)[warp_id_19 + 32U] = _S2641.z;
        uint _S2661 = __ballot_sync(_S2216, true);
        _S2216 = _S2661;
    }
    else
    {
        uint _S2662 = __ballot_sync(_S2216, true);
        _S2216 = _S2662;
    }
    __syncthreads();
    bool _S2663 = warp_id_19 < 3U;
    uint _S2664 = __ballot_sync(_S2216, _S2663);
    if(_S2663)
    {
        bool _S2665 = lane_id_2 < 16U;
        uint _S2666 = __ballot_sync(_S2664, _S2665);
        if(_S2665)
        {
            float _S2667 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_19 * 16U];
            uint _S2668 = __ballot_sync(_S2664, true);
            v_2 = _S2667;
            _S2239 = _S2668;
        }
        else
        {
            uint _S2669 = __ballot_sync(_S2664, true);
            v_2 = 0.0f;
            _S2239 = _S2669;
        }
        float _S2670 = WaveActiveSum_0(v_2, _S2239);
        uint _S2671 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2672 = _S2670 != 0.0f;
            uint _S2673 = __ballot_sync(_S2239, true);
            _S2240 = _S2672;
            _S2239 = _S2673;
        }
        else
        {
            uint _S2674 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2674;
        }
        uint _S2675 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2676 = warp_id_19 == 0U;
            uint _S2677 = __ballot_sync(_S2675, _S2676);
            if(_S2676)
            {
                float _S2678 = atomicAdd(&((v_coeffs_31 + 8U)->x), _S2670);
                uint _S2679 = __ballot_sync(_S2675, true);
            }
            else
            {
                uint _S2680 = _S2675 & (~_S2677);
                bool _S2681 = warp_id_19 == 1U;
                uint _S2682 = __ballot_sync(_S2680, _S2681);
                if(_S2681)
                {
                    float _S2683 = atomicAdd(&((v_coeffs_31 + 8U)->y), _S2670);
                    uint _S2684 = __ballot_sync(_S2680, true);
                }
                else
                {
                    float _S2685 = atomicAdd(&((v_coeffs_31 + 8U)->z), _S2670);
                    uint _S2686 = __ballot_sync(_S2680, true);
                }
                uint _S2687 = __ballot_sync(_S2675, true);
            }
            uint _S2688 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2689 = __ballot_sync(_S2239, true);
        }
        uint _S2690 = __ballot_sync(_S2216, true);
        _S2216 = _S2690;
    }
    else
    {
        uint _S2691 = __ballot_sync(_S2216, true);
        _S2216 = _S2691;
    }
    float3  _S2692 = make_float3 (pSH10_12) * v_colors_31;
    float3  _S2693 = _S2692;
    bool _S2694 = (F32_isfinite((_S2692.x)));
    uint _S2695 = __ballot_sync(_S2216, _S2694);
    if(_S2694)
    {
        float _S2696 = _S2693.x;
        uint _S2697 = __ballot_sync(_S2216, true);
        v_2 = _S2696;
        _S2216 = _S2697;
    }
    else
    {
        uint _S2698 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2698;
    }
    *&((&_S2693)->x) = v_2;
    bool _S2699 = (F32_isfinite((_S2693.y)));
    uint _S2700 = __ballot_sync(_S2216, _S2699);
    if(_S2699)
    {
        float _S2701 = _S2693.y;
        uint _S2702 = __ballot_sync(_S2216, true);
        v_2 = _S2701;
        _S2216 = _S2702;
    }
    else
    {
        uint _S2703 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2703;
    }
    *&((&_S2693)->y) = v_2;
    bool _S2704 = (F32_isfinite((_S2693.z)));
    uint _S2705 = __ballot_sync(_S2216, _S2704);
    if(_S2704)
    {
        float _S2706 = _S2693.z;
        uint _S2707 = __ballot_sync(_S2216, true);
        v_2 = _S2706;
        _S2216 = _S2707;
    }
    else
    {
        uint _S2708 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2708;
    }
    *&((&_S2693)->z) = v_2;
    float _S2709 = WaveActiveSum_0(_S2693.x, _S2216);
    *&((&_S2693)->x) = _S2709;
    float _S2710 = WaveActiveSum_0(_S2693.y, _S2216);
    *&((&_S2693)->y) = _S2710;
    float _S2711 = WaveActiveSum_0(_S2693.z, _S2216);
    *&((&_S2693)->z) = _S2711;
    uint warp_id_20 = thread_id_3 / 32U;
    uint _S2712 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_20] = _S2693.x;
        (*&_sh_block_reduce_shared_0)[warp_id_20 + 16U] = _S2693.y;
        (*&_sh_block_reduce_shared_0)[warp_id_20 + 32U] = _S2693.z;
        uint _S2713 = __ballot_sync(_S2216, true);
        _S2216 = _S2713;
    }
    else
    {
        uint _S2714 = __ballot_sync(_S2216, true);
        _S2216 = _S2714;
    }
    __syncthreads();
    bool _S2715 = warp_id_20 < 3U;
    uint _S2716 = __ballot_sync(_S2216, _S2715);
    if(_S2715)
    {
        bool _S2717 = lane_id_2 < 16U;
        uint _S2718 = __ballot_sync(_S2716, _S2717);
        if(_S2717)
        {
            float _S2719 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_20 * 16U];
            uint _S2720 = __ballot_sync(_S2716, true);
            v_2 = _S2719;
            _S2239 = _S2720;
        }
        else
        {
            uint _S2721 = __ballot_sync(_S2716, true);
            v_2 = 0.0f;
            _S2239 = _S2721;
        }
        float _S2722 = WaveActiveSum_0(v_2, _S2239);
        uint _S2723 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2724 = _S2722 != 0.0f;
            uint _S2725 = __ballot_sync(_S2239, true);
            _S2240 = _S2724;
            _S2239 = _S2725;
        }
        else
        {
            uint _S2726 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2726;
        }
        uint _S2727 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2728 = warp_id_20 == 0U;
            uint _S2729 = __ballot_sync(_S2727, _S2728);
            if(_S2728)
            {
                float _S2730 = atomicAdd(&((v_coeffs_31 + 9U)->x), _S2722);
                uint _S2731 = __ballot_sync(_S2727, true);
            }
            else
            {
                uint _S2732 = _S2727 & (~_S2729);
                bool _S2733 = warp_id_20 == 1U;
                uint _S2734 = __ballot_sync(_S2732, _S2733);
                if(_S2733)
                {
                    float _S2735 = atomicAdd(&((v_coeffs_31 + 9U)->y), _S2722);
                    uint _S2736 = __ballot_sync(_S2732, true);
                }
                else
                {
                    float _S2737 = atomicAdd(&((v_coeffs_31 + 9U)->z), _S2722);
                    uint _S2738 = __ballot_sync(_S2732, true);
                }
                uint _S2739 = __ballot_sync(_S2727, true);
            }
            uint _S2740 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2741 = __ballot_sync(_S2239, true);
        }
        uint _S2742 = __ballot_sync(_S2216, true);
        _S2216 = _S2742;
    }
    else
    {
        uint _S2743 = __ballot_sync(_S2216, true);
        _S2216 = _S2743;
    }
    float3  _S2744 = make_float3 (pSH11_12) * v_colors_31;
    float3  _S2745 = _S2744;
    bool _S2746 = (F32_isfinite((_S2744.x)));
    uint _S2747 = __ballot_sync(_S2216, _S2746);
    if(_S2746)
    {
        float _S2748 = _S2745.x;
        uint _S2749 = __ballot_sync(_S2216, true);
        v_2 = _S2748;
        _S2216 = _S2749;
    }
    else
    {
        uint _S2750 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2750;
    }
    *&((&_S2745)->x) = v_2;
    bool _S2751 = (F32_isfinite((_S2745.y)));
    uint _S2752 = __ballot_sync(_S2216, _S2751);
    if(_S2751)
    {
        float _S2753 = _S2745.y;
        uint _S2754 = __ballot_sync(_S2216, true);
        v_2 = _S2753;
        _S2216 = _S2754;
    }
    else
    {
        uint _S2755 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2755;
    }
    *&((&_S2745)->y) = v_2;
    bool _S2756 = (F32_isfinite((_S2745.z)));
    uint _S2757 = __ballot_sync(_S2216, _S2756);
    if(_S2756)
    {
        float _S2758 = _S2745.z;
        uint _S2759 = __ballot_sync(_S2216, true);
        v_2 = _S2758;
        _S2216 = _S2759;
    }
    else
    {
        uint _S2760 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2760;
    }
    *&((&_S2745)->z) = v_2;
    float _S2761 = WaveActiveSum_0(_S2745.x, _S2216);
    *&((&_S2745)->x) = _S2761;
    float _S2762 = WaveActiveSum_0(_S2745.y, _S2216);
    *&((&_S2745)->y) = _S2762;
    float _S2763 = WaveActiveSum_0(_S2745.z, _S2216);
    *&((&_S2745)->z) = _S2763;
    uint warp_id_21 = thread_id_3 / 32U;
    uint _S2764 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_21] = _S2745.x;
        (*&_sh_block_reduce_shared_0)[warp_id_21 + 16U] = _S2745.y;
        (*&_sh_block_reduce_shared_0)[warp_id_21 + 32U] = _S2745.z;
        uint _S2765 = __ballot_sync(_S2216, true);
        _S2216 = _S2765;
    }
    else
    {
        uint _S2766 = __ballot_sync(_S2216, true);
        _S2216 = _S2766;
    }
    __syncthreads();
    bool _S2767 = warp_id_21 < 3U;
    uint _S2768 = __ballot_sync(_S2216, _S2767);
    if(_S2767)
    {
        bool _S2769 = lane_id_2 < 16U;
        uint _S2770 = __ballot_sync(_S2768, _S2769);
        if(_S2769)
        {
            float _S2771 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_21 * 16U];
            uint _S2772 = __ballot_sync(_S2768, true);
            v_2 = _S2771;
            _S2239 = _S2772;
        }
        else
        {
            uint _S2773 = __ballot_sync(_S2768, true);
            v_2 = 0.0f;
            _S2239 = _S2773;
        }
        float _S2774 = WaveActiveSum_0(v_2, _S2239);
        uint _S2775 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2776 = _S2774 != 0.0f;
            uint _S2777 = __ballot_sync(_S2239, true);
            _S2240 = _S2776;
            _S2239 = _S2777;
        }
        else
        {
            uint _S2778 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2778;
        }
        uint _S2779 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2780 = warp_id_21 == 0U;
            uint _S2781 = __ballot_sync(_S2779, _S2780);
            if(_S2780)
            {
                float _S2782 = atomicAdd(&((v_coeffs_31 + 10U)->x), _S2774);
                uint _S2783 = __ballot_sync(_S2779, true);
            }
            else
            {
                uint _S2784 = _S2779 & (~_S2781);
                bool _S2785 = warp_id_21 == 1U;
                uint _S2786 = __ballot_sync(_S2784, _S2785);
                if(_S2785)
                {
                    float _S2787 = atomicAdd(&((v_coeffs_31 + 10U)->y), _S2774);
                    uint _S2788 = __ballot_sync(_S2784, true);
                }
                else
                {
                    float _S2789 = atomicAdd(&((v_coeffs_31 + 10U)->z), _S2774);
                    uint _S2790 = __ballot_sync(_S2784, true);
                }
                uint _S2791 = __ballot_sync(_S2779, true);
            }
            uint _S2792 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2793 = __ballot_sync(_S2239, true);
        }
        uint _S2794 = __ballot_sync(_S2216, true);
        _S2216 = _S2794;
    }
    else
    {
        uint _S2795 = __ballot_sync(_S2216, true);
        _S2216 = _S2795;
    }
    float3  _S2796 = make_float3 (pSH12_13) * v_colors_31;
    float3  _S2797 = _S2796;
    bool _S2798 = (F32_isfinite((_S2796.x)));
    uint _S2799 = __ballot_sync(_S2216, _S2798);
    if(_S2798)
    {
        float _S2800 = _S2797.x;
        uint _S2801 = __ballot_sync(_S2216, true);
        v_2 = _S2800;
        _S2216 = _S2801;
    }
    else
    {
        uint _S2802 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2802;
    }
    *&((&_S2797)->x) = v_2;
    bool _S2803 = (F32_isfinite((_S2797.y)));
    uint _S2804 = __ballot_sync(_S2216, _S2803);
    if(_S2803)
    {
        float _S2805 = _S2797.y;
        uint _S2806 = __ballot_sync(_S2216, true);
        v_2 = _S2805;
        _S2216 = _S2806;
    }
    else
    {
        uint _S2807 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2807;
    }
    *&((&_S2797)->y) = v_2;
    bool _S2808 = (F32_isfinite((_S2797.z)));
    uint _S2809 = __ballot_sync(_S2216, _S2808);
    if(_S2808)
    {
        float _S2810 = _S2797.z;
        uint _S2811 = __ballot_sync(_S2216, true);
        v_2 = _S2810;
        _S2216 = _S2811;
    }
    else
    {
        uint _S2812 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2812;
    }
    *&((&_S2797)->z) = v_2;
    float _S2813 = WaveActiveSum_0(_S2797.x, _S2216);
    *&((&_S2797)->x) = _S2813;
    float _S2814 = WaveActiveSum_0(_S2797.y, _S2216);
    *&((&_S2797)->y) = _S2814;
    float _S2815 = WaveActiveSum_0(_S2797.z, _S2216);
    *&((&_S2797)->z) = _S2815;
    uint warp_id_22 = thread_id_3 / 32U;
    uint _S2816 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_22] = _S2797.x;
        (*&_sh_block_reduce_shared_0)[warp_id_22 + 16U] = _S2797.y;
        (*&_sh_block_reduce_shared_0)[warp_id_22 + 32U] = _S2797.z;
        uint _S2817 = __ballot_sync(_S2216, true);
        _S2216 = _S2817;
    }
    else
    {
        uint _S2818 = __ballot_sync(_S2216, true);
        _S2216 = _S2818;
    }
    __syncthreads();
    bool _S2819 = warp_id_22 < 3U;
    uint _S2820 = __ballot_sync(_S2216, _S2819);
    if(_S2819)
    {
        bool _S2821 = lane_id_2 < 16U;
        uint _S2822 = __ballot_sync(_S2820, _S2821);
        if(_S2821)
        {
            float _S2823 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_22 * 16U];
            uint _S2824 = __ballot_sync(_S2820, true);
            v_2 = _S2823;
            _S2239 = _S2824;
        }
        else
        {
            uint _S2825 = __ballot_sync(_S2820, true);
            v_2 = 0.0f;
            _S2239 = _S2825;
        }
        float _S2826 = WaveActiveSum_0(v_2, _S2239);
        uint _S2827 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2828 = _S2826 != 0.0f;
            uint _S2829 = __ballot_sync(_S2239, true);
            _S2240 = _S2828;
            _S2239 = _S2829;
        }
        else
        {
            uint _S2830 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2830;
        }
        uint _S2831 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2832 = warp_id_22 == 0U;
            uint _S2833 = __ballot_sync(_S2831, _S2832);
            if(_S2832)
            {
                float _S2834 = atomicAdd(&((v_coeffs_31 + 11U)->x), _S2826);
                uint _S2835 = __ballot_sync(_S2831, true);
            }
            else
            {
                uint _S2836 = _S2831 & (~_S2833);
                bool _S2837 = warp_id_22 == 1U;
                uint _S2838 = __ballot_sync(_S2836, _S2837);
                if(_S2837)
                {
                    float _S2839 = atomicAdd(&((v_coeffs_31 + 11U)->y), _S2826);
                    uint _S2840 = __ballot_sync(_S2836, true);
                }
                else
                {
                    float _S2841 = atomicAdd(&((v_coeffs_31 + 11U)->z), _S2826);
                    uint _S2842 = __ballot_sync(_S2836, true);
                }
                uint _S2843 = __ballot_sync(_S2831, true);
            }
            uint _S2844 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2845 = __ballot_sync(_S2239, true);
        }
        uint _S2846 = __ballot_sync(_S2216, true);
        _S2216 = _S2846;
    }
    else
    {
        uint _S2847 = __ballot_sync(_S2216, true);
        _S2216 = _S2847;
    }
    float3  _S2848 = make_float3 (pSH13_12) * v_colors_31;
    float3  _S2849 = _S2848;
    bool _S2850 = (F32_isfinite((_S2848.x)));
    uint _S2851 = __ballot_sync(_S2216, _S2850);
    if(_S2850)
    {
        float _S2852 = _S2849.x;
        uint _S2853 = __ballot_sync(_S2216, true);
        v_2 = _S2852;
        _S2216 = _S2853;
    }
    else
    {
        uint _S2854 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2854;
    }
    *&((&_S2849)->x) = v_2;
    bool _S2855 = (F32_isfinite((_S2849.y)));
    uint _S2856 = __ballot_sync(_S2216, _S2855);
    if(_S2855)
    {
        float _S2857 = _S2849.y;
        uint _S2858 = __ballot_sync(_S2216, true);
        v_2 = _S2857;
        _S2216 = _S2858;
    }
    else
    {
        uint _S2859 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2859;
    }
    *&((&_S2849)->y) = v_2;
    bool _S2860 = (F32_isfinite((_S2849.z)));
    uint _S2861 = __ballot_sync(_S2216, _S2860);
    if(_S2860)
    {
        float _S2862 = _S2849.z;
        uint _S2863 = __ballot_sync(_S2216, true);
        v_2 = _S2862;
        _S2216 = _S2863;
    }
    else
    {
        uint _S2864 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2864;
    }
    *&((&_S2849)->z) = v_2;
    float _S2865 = WaveActiveSum_0(_S2849.x, _S2216);
    *&((&_S2849)->x) = _S2865;
    float _S2866 = WaveActiveSum_0(_S2849.y, _S2216);
    *&((&_S2849)->y) = _S2866;
    float _S2867 = WaveActiveSum_0(_S2849.z, _S2216);
    *&((&_S2849)->z) = _S2867;
    uint warp_id_23 = thread_id_3 / 32U;
    uint _S2868 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_23] = _S2849.x;
        (*&_sh_block_reduce_shared_0)[warp_id_23 + 16U] = _S2849.y;
        (*&_sh_block_reduce_shared_0)[warp_id_23 + 32U] = _S2849.z;
        uint _S2869 = __ballot_sync(_S2216, true);
        _S2216 = _S2869;
    }
    else
    {
        uint _S2870 = __ballot_sync(_S2216, true);
        _S2216 = _S2870;
    }
    __syncthreads();
    bool _S2871 = warp_id_23 < 3U;
    uint _S2872 = __ballot_sync(_S2216, _S2871);
    if(_S2871)
    {
        bool _S2873 = lane_id_2 < 16U;
        uint _S2874 = __ballot_sync(_S2872, _S2873);
        if(_S2873)
        {
            float _S2875 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_23 * 16U];
            uint _S2876 = __ballot_sync(_S2872, true);
            v_2 = _S2875;
            _S2239 = _S2876;
        }
        else
        {
            uint _S2877 = __ballot_sync(_S2872, true);
            v_2 = 0.0f;
            _S2239 = _S2877;
        }
        float _S2878 = WaveActiveSum_0(v_2, _S2239);
        uint _S2879 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2880 = _S2878 != 0.0f;
            uint _S2881 = __ballot_sync(_S2239, true);
            _S2240 = _S2880;
            _S2239 = _S2881;
        }
        else
        {
            uint _S2882 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2882;
        }
        uint _S2883 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2884 = warp_id_23 == 0U;
            uint _S2885 = __ballot_sync(_S2883, _S2884);
            if(_S2884)
            {
                float _S2886 = atomicAdd(&((v_coeffs_31 + 12U)->x), _S2878);
                uint _S2887 = __ballot_sync(_S2883, true);
            }
            else
            {
                uint _S2888 = _S2883 & (~_S2885);
                bool _S2889 = warp_id_23 == 1U;
                uint _S2890 = __ballot_sync(_S2888, _S2889);
                if(_S2889)
                {
                    float _S2891 = atomicAdd(&((v_coeffs_31 + 12U)->y), _S2878);
                    uint _S2892 = __ballot_sync(_S2888, true);
                }
                else
                {
                    float _S2893 = atomicAdd(&((v_coeffs_31 + 12U)->z), _S2878);
                    uint _S2894 = __ballot_sync(_S2888, true);
                }
                uint _S2895 = __ballot_sync(_S2883, true);
            }
            uint _S2896 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2897 = __ballot_sync(_S2239, true);
        }
        uint _S2898 = __ballot_sync(_S2216, true);
        _S2216 = _S2898;
    }
    else
    {
        uint _S2899 = __ballot_sync(_S2216, true);
        _S2216 = _S2899;
    }
    float3  _S2900 = make_float3 (pSH14_12) * v_colors_31;
    float3  _S2901 = _S2900;
    bool _S2902 = (F32_isfinite((_S2900.x)));
    uint _S2903 = __ballot_sync(_S2216, _S2902);
    if(_S2902)
    {
        float _S2904 = _S2901.x;
        uint _S2905 = __ballot_sync(_S2216, true);
        v_2 = _S2904;
        _S2216 = _S2905;
    }
    else
    {
        uint _S2906 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2906;
    }
    *&((&_S2901)->x) = v_2;
    bool _S2907 = (F32_isfinite((_S2901.y)));
    uint _S2908 = __ballot_sync(_S2216, _S2907);
    if(_S2907)
    {
        float _S2909 = _S2901.y;
        uint _S2910 = __ballot_sync(_S2216, true);
        v_2 = _S2909;
        _S2216 = _S2910;
    }
    else
    {
        uint _S2911 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2911;
    }
    *&((&_S2901)->y) = v_2;
    bool _S2912 = (F32_isfinite((_S2901.z)));
    uint _S2913 = __ballot_sync(_S2216, _S2912);
    if(_S2912)
    {
        float _S2914 = _S2901.z;
        uint _S2915 = __ballot_sync(_S2216, true);
        v_2 = _S2914;
        _S2216 = _S2915;
    }
    else
    {
        uint _S2916 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2916;
    }
    *&((&_S2901)->z) = v_2;
    float _S2917 = WaveActiveSum_0(_S2901.x, _S2216);
    *&((&_S2901)->x) = _S2917;
    float _S2918 = WaveActiveSum_0(_S2901.y, _S2216);
    *&((&_S2901)->y) = _S2918;
    float _S2919 = WaveActiveSum_0(_S2901.z, _S2216);
    *&((&_S2901)->z) = _S2919;
    uint warp_id_24 = thread_id_3 / 32U;
    uint _S2920 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_24] = _S2901.x;
        (*&_sh_block_reduce_shared_0)[warp_id_24 + 16U] = _S2901.y;
        (*&_sh_block_reduce_shared_0)[warp_id_24 + 32U] = _S2901.z;
        uint _S2921 = __ballot_sync(_S2216, true);
        _S2216 = _S2921;
    }
    else
    {
        uint _S2922 = __ballot_sync(_S2216, true);
        _S2216 = _S2922;
    }
    __syncthreads();
    bool _S2923 = warp_id_24 < 3U;
    uint _S2924 = __ballot_sync(_S2216, _S2923);
    if(_S2923)
    {
        bool _S2925 = lane_id_2 < 16U;
        uint _S2926 = __ballot_sync(_S2924, _S2925);
        if(_S2925)
        {
            float _S2927 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_24 * 16U];
            uint _S2928 = __ballot_sync(_S2924, true);
            v_2 = _S2927;
            _S2239 = _S2928;
        }
        else
        {
            uint _S2929 = __ballot_sync(_S2924, true);
            v_2 = 0.0f;
            _S2239 = _S2929;
        }
        float _S2930 = WaveActiveSum_0(v_2, _S2239);
        uint _S2931 = __ballot_sync(_S2239, _S2233);
        if(_S2233)
        {
            bool _S2932 = _S2930 != 0.0f;
            uint _S2933 = __ballot_sync(_S2239, true);
            _S2240 = _S2932;
            _S2239 = _S2933;
        }
        else
        {
            uint _S2934 = __ballot_sync(_S2239, true);
            _S2240 = false;
            _S2239 = _S2934;
        }
        uint _S2935 = __ballot_sync(_S2239, _S2240);
        if(_S2240)
        {
            bool _S2936 = warp_id_24 == 0U;
            uint _S2937 = __ballot_sync(_S2935, _S2936);
            if(_S2936)
            {
                float _S2938 = atomicAdd(&((v_coeffs_31 + 13U)->x), _S2930);
                uint _S2939 = __ballot_sync(_S2935, true);
            }
            else
            {
                uint _S2940 = _S2935 & (~_S2937);
                bool _S2941 = warp_id_24 == 1U;
                uint _S2942 = __ballot_sync(_S2940, _S2941);
                if(_S2941)
                {
                    float _S2943 = atomicAdd(&((v_coeffs_31 + 13U)->y), _S2930);
                    uint _S2944 = __ballot_sync(_S2940, true);
                }
                else
                {
                    float _S2945 = atomicAdd(&((v_coeffs_31 + 13U)->z), _S2930);
                    uint _S2946 = __ballot_sync(_S2940, true);
                }
                uint _S2947 = __ballot_sync(_S2935, true);
            }
            uint _S2948 = __ballot_sync(_S2239, true);
        }
        else
        {
            uint _S2949 = __ballot_sync(_S2239, true);
        }
        uint _S2950 = __ballot_sync(_S2216, true);
        _S2216 = _S2950;
    }
    else
    {
        uint _S2951 = __ballot_sync(_S2216, true);
        _S2216 = _S2951;
    }
    float3  _S2952 = make_float3 (pSH15_12) * v_colors_31;
    float3  _S2953 = _S2952;
    bool _S2954 = (F32_isfinite((_S2952.x)));
    uint _S2955 = __ballot_sync(_S2216, _S2954);
    if(_S2954)
    {
        float _S2956 = _S2953.x;
        uint _S2957 = __ballot_sync(_S2216, true);
        v_2 = _S2956;
        _S2216 = _S2957;
    }
    else
    {
        uint _S2958 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2958;
    }
    *&((&_S2953)->x) = v_2;
    bool _S2959 = (F32_isfinite((_S2953.y)));
    uint _S2960 = __ballot_sync(_S2216, _S2959);
    if(_S2959)
    {
        float _S2961 = _S2953.y;
        uint _S2962 = __ballot_sync(_S2216, true);
        v_2 = _S2961;
        _S2216 = _S2962;
    }
    else
    {
        uint _S2963 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2963;
    }
    *&((&_S2953)->y) = v_2;
    bool _S2964 = (F32_isfinite((_S2953.z)));
    uint _S2965 = __ballot_sync(_S2216, _S2964);
    if(_S2964)
    {
        float _S2966 = _S2953.z;
        uint _S2967 = __ballot_sync(_S2216, true);
        v_2 = _S2966;
        _S2216 = _S2967;
    }
    else
    {
        uint _S2968 = __ballot_sync(_S2216, true);
        v_2 = 0.0f;
        _S2216 = _S2968;
    }
    *&((&_S2953)->z) = v_2;
    float _S2969 = WaveActiveSum_0(_S2953.x, _S2216);
    *&((&_S2953)->x) = _S2969;
    float _S2970 = WaveActiveSum_0(_S2953.y, _S2216);
    *&((&_S2953)->y) = _S2970;
    float _S2971 = WaveActiveSum_0(_S2953.z, _S2216);
    *&((&_S2953)->z) = _S2971;
    uint warp_id_25 = thread_id_3 / 32U;
    uint _S2972 = __ballot_sync(_S2216, _S2233);
    if(_S2233)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_25] = _S2953.x;
        (*&_sh_block_reduce_shared_0)[warp_id_25 + 16U] = _S2953.y;
        (*&_sh_block_reduce_shared_0)[warp_id_25 + 32U] = _S2953.z;
        uint _S2973 = __ballot_sync(_S2216, true);
        _S2216 = _S2973;
    }
    else
    {
        uint _S2974 = __ballot_sync(_S2216, true);
        _S2216 = _S2974;
    }
    __syncthreads();
    bool _S2975 = warp_id_25 < 3U;
    uint _S2976 = __ballot_sync(_S2216, _S2975);
    if(_S2975)
    {
        bool _S2977 = lane_id_2 < 16U;
        uint _S2978 = __ballot_sync(_S2976, _S2977);
        if(_S2977)
        {
            float _S2979 = (*&_sh_block_reduce_shared_0)[lane_id_2 + warp_id_25 * 16U];
            uint _S2980 = __ballot_sync(_S2976, true);
            v_2 = _S2979;
            _S2216 = _S2980;
        }
        else
        {
            uint _S2981 = __ballot_sync(_S2976, true);
            v_2 = 0.0f;
            _S2216 = _S2981;
        }
        float _S2982 = WaveActiveSum_0(v_2, _S2216);
        uint _S2983 = __ballot_sync(_S2216, _S2233);
        if(_S2233)
        {
            _S2240 = _S2982 != 0.0f;
        }
        else
        {
            _S2240 = false;
        }
        if(_S2240)
        {
            if(warp_id_25 == 0U)
            {
                float _S2984 = atomicAdd(&((v_coeffs_31 + 14U)->x), _S2982);
            }
            else
            {
                if(warp_id_25 == 1U)
                {
                    float _S2985 = atomicAdd(&((v_coeffs_31 + 14U)->y), _S2982);
                }
                else
                {
                    float _S2986 = atomicAdd(&((v_coeffs_31 + 14U)->z), _S2982);
                }
            }
        }
    }
    float fTmp0C_z_10 = -4.57045793533325195f * z_33;
    float _S2987 = x_36 * _S2375;
    float _S2988 = y_34 * _S2375;
    float pSH14_x_10 = fTmp1B_15 * _S2375;
    float3  * _S2989 = coeffs_45 + int(8);
    float3  * _S2990 = coeffs_45 + int(14);
    float3  * _S2991 = coeffs_45 + int(9);
    float3  * _S2992 = coeffs_45 + int(13);
    float3  * _S2993 = coeffs_45 + int(12);
    float3  * _S2994 = coeffs_45 + int(10);
    float3  dir_n_24 = make_float3 (x_36, y_34, z_33);
    float3  v_dir_n_40 = make_float3 (v_x_26 + dot_0(v_colors_31, make_float3 (-0.59004360437393188f * (fS1_15 + x_36 * fS1_x_10 + _S2988)) * *_S2989 + make_float3 (-0.59004360437393188f * (fC1_15 + _S2987 - y_34 * fS1_x_10)) * *_S2990 + make_float3 (fTmp1B_15 * fS1_x_10) * *_S2991 + make_float3 (pSH14_x_10) * *_S2992 + make_float3 (fTmp0C_15) * *_S2993), v_y_26 + dot_0(v_colors_31, make_float3 (-0.59004360437393188f * (_S2987 + fC1_15 + y_34 * fC1_y_10)) * *_S2989 + make_float3 (-0.59004360437393188f * (x_36 * fC1_y_10 - fS1_15 - _S2988)) * *_S2990 + make_float3 (pSH14_x_10) * *_S2991 + make_float3 (fTmp1B_15 * fC1_y_10) * *_S2992 + make_float3 (fTmp0C_15) * *_S2994), v_z_26 + dot_0(v_colors_31, make_float3 (5.59764480590820312f * z2_15 - 1.11952900886535645f) * *(coeffs_45 + int(11)) + make_float3 (fTmp0C_z_10 * x_36) * *_S2993 + make_float3 (fTmp0C_z_10 * y_34) * *_S2994 + make_float3 (1.44530570507049561f * fC1_15) * *_S2992 + make_float3 (1.44530570507049561f * fS1_15) * *_S2991));
    *v_dir_11 = *v_dir_11 + (v_dir_n_40 - make_float3 (dot_0(v_dir_n_40, dir_n_24)) * dir_n_24) * make_float3 (inorm_8);
    return;
}

inline __device__ float3  sh4_to_color_dir(float3  dir_16, float3  coeff_dc_46, float3  * coeffs_46)
{
    float _S2995 = dir_16.x;
    float _S2996 = dir_16.y;
    float _S2997 = dir_16.z;
    float inv_norm_27 = (F32_rsqrt((_S2995 * _S2995 + _S2996 * _S2996 + _S2997 * _S2997)));
    float x_37 = _S2995 * inv_norm_27;
    float y_35 = _S2996 * inv_norm_27;
    float z_34 = _S2997 * inv_norm_27;
    float z2_16 = z_34 * z_34;
    float fTmp0B_26 = -1.09254848957061768f * z_34;
    float fC1_16 = x_37 * x_37 - y_35 * y_35;
    float fS1_16 = 2.0f * x_37 * y_35;
    float pSH6_22 = 0.94617468118667603f * z2_16 - 0.31539157032966614f;
    float fTmp0C_16 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_16 = 1.44530570507049561f * z_34;
    float fC2_6 = x_37 * fC1_16 - y_35 * fS1_16;
    float fS2_6 = x_37 * fS1_16 + y_35 * fC1_16;
    float pSH12_14 = z_34 * (1.86588168144226074f * z2_16 - 1.11952900886535645f);
    float fTmp0D_6 = z_34 * (-4.68332576751708984f * z2_16 + 2.00713968276977539f);
    float fTmp1C_6 = 3.31161141395568848f * z2_16 - 0.47308734059333801f;
    float fTmp2B_6 = -1.77013075351715088f * z_34;
    return make_float3 (0.282094806432724f) * coeff_dc_46 + make_float3 (0.48860251903533936f) * (make_float3 (- y_35) * *(coeffs_46 + int(0)) + make_float3 (z_34) * *(coeffs_46 + int(1)) - make_float3 (x_37) * *(coeffs_46 + int(2))) + (make_float3 (0.54627424478530884f * fS1_16) * *(coeffs_46 + int(3)) + make_float3 (fTmp0B_26 * y_35) * *(coeffs_46 + int(4)) + make_float3 (pSH6_22) * *(coeffs_46 + int(5)) + make_float3 (fTmp0B_26 * x_37) * *(coeffs_46 + int(6)) + make_float3 (0.54627424478530884f * fC1_16) * *(coeffs_46 + int(7))) + (make_float3 (-0.59004360437393188f * fS2_6) * *(coeffs_46 + int(8)) + make_float3 (fTmp1B_16 * fS1_16) * *(coeffs_46 + int(9)) + make_float3 (fTmp0C_16 * y_35) * *(coeffs_46 + int(10)) + make_float3 (pSH12_14) * *(coeffs_46 + int(11)) + make_float3 (fTmp0C_16 * x_37) * *(coeffs_46 + int(12)) + make_float3 (fTmp1B_16 * fC1_16) * *(coeffs_46 + int(13)) + make_float3 (-0.59004360437393188f * fC2_6) * *(coeffs_46 + int(14))) + (make_float3 (0.62583571672439575f * (x_37 * fS2_6 + y_35 * fC2_6)) * *(coeffs_46 + int(15)) + make_float3 (fTmp2B_6 * fS2_6) * *(coeffs_46 + int(16)) + make_float3 (fTmp1C_6 * fS1_16) * *(coeffs_46 + int(17)) + make_float3 (fTmp0D_6 * y_35) * *(coeffs_46 + int(18)) + make_float3 (1.9843134880065918f * z_34 * pSH12_14 - 1.00623059272766113f * pSH6_22) * *(coeffs_46 + int(19)) + make_float3 (fTmp0D_6 * x_37) * *(coeffs_46 + int(20)) + make_float3 (fTmp1C_6 * fC1_16) * *(coeffs_46 + int(21)) + make_float3 (fTmp2B_6 * fC2_6) * *(coeffs_46 + int(22)) + make_float3 (0.62583571672439575f * (x_37 * fC2_6 - y_35 * fS2_6)) * *(coeffs_46 + int(23)));
}

inline __device__ void sh4_to_color_dir_vjp_inplace(float3  dir_17, float3  coeff_dc_47, float3  * coeffs_47, float3  v_colors_32, float3  * v_coeff_dc_32, float3  * v_coeffs_32, float3  * v_dir_12)
{
    *v_coeff_dc_32 = *v_coeff_dc_32 + make_float3 (0.282094806432724f) * v_colors_32;
    float _S2998 = dir_17.x;
    float _S2999 = dir_17.y;
    float _S3000 = dir_17.z;
    float inorm_9 = (F32_rsqrt((_S2998 * _S2998 + _S2999 * _S2999 + _S3000 * _S3000)));
    float x_38 = _S2998 * inorm_9;
    float y_36 = _S2999 * inorm_9;
    float z_35 = _S3000 * inorm_9;
    float3  * _S3001 = v_coeffs_32 + int(0);
    *_S3001 = *_S3001 + make_float3 (-0.48860251903533936f * y_36) * v_colors_32;
    float3  * _S3002 = v_coeffs_32 + int(1);
    *_S3002 = *_S3002 + make_float3 (0.48860251903533936f * z_35) * v_colors_32;
    float3  * _S3003 = v_coeffs_32 + int(2);
    *_S3003 = *_S3003 + make_float3 (-0.48860251903533936f * x_38) * v_colors_32;
    float _S3004 = -0.48860251903533936f * dot_0(*(coeffs_47 + int(2)), v_colors_32);
    float _S3005 = -0.48860251903533936f * dot_0(*(coeffs_47 + int(0)), v_colors_32);
    float _S3006 = 0.48860251903533936f * dot_0(*(coeffs_47 + int(1)), v_colors_32);
    float z2_17 = z_35 * z_35;
    float fTmp0B_27 = -1.09254848957061768f * z_35;
    float fC1_17 = x_38 * x_38 - y_36 * y_36;
    float _S3007 = 2.0f * x_38;
    float fS1_17 = _S3007 * y_36;
    float pSH6_23 = 0.94617468118667603f * z2_17 - 0.31539157032966614f;
    float pSH7_21 = fTmp0B_27 * x_38;
    float pSH5_21 = fTmp0B_27 * y_36;
    float pSH8_21 = 0.54627424478530884f * fC1_17;
    float3  * _S3008 = v_coeffs_32 + int(3);
    *_S3008 = *_S3008 + make_float3 (0.54627424478530884f * fS1_17) * v_colors_32;
    float3  * _S3009 = v_coeffs_32 + int(4);
    *_S3009 = *_S3009 + make_float3 (pSH5_21) * v_colors_32;
    float3  * _S3010 = v_coeffs_32 + int(5);
    *_S3010 = *_S3010 + make_float3 (pSH6_23) * v_colors_32;
    float3  * _S3011 = v_coeffs_32 + int(6);
    *_S3011 = *_S3011 + make_float3 (pSH7_21) * v_colors_32;
    float3  * _S3012 = v_coeffs_32 + int(7);
    *_S3012 = *_S3012 + make_float3 (pSH8_21) * v_colors_32;
    float fC1_y_11 = -2.0f * y_36;
    float fS1_x_11 = 2.0f * y_36;
    float pSH6_z_8 = 1.89234936237335205f * z_35;
    float pSH8_x_18 = 0.54627424478530884f * _S3007;
    float3  * _S3013 = coeffs_47 + int(3);
    float3  * _S3014 = coeffs_47 + int(7);
    float3  * _S3015 = coeffs_47 + int(6);
    float v_x_27 = _S3004 + dot_0(v_colors_32, make_float3 (0.54627424478530884f * fS1_x_11) * *_S3013 + make_float3 (pSH8_x_18) * *_S3014 + make_float3 (fTmp0B_27) * *_S3015);
    float3  * _S3016 = coeffs_47 + int(4);
    float v_y_27 = _S3005 + dot_0(v_colors_32, make_float3 (pSH8_x_18) * *_S3013 + make_float3 (0.54627424478530884f * fC1_y_11) * *_S3014 + make_float3 (fTmp0B_27) * *_S3016);
    float v_z_27 = _S3006 + dot_0(v_colors_32, make_float3 (pSH6_z_8) * *(coeffs_47 + int(5)) + make_float3 (-1.09254848957061768f * x_38) * *_S3015 + make_float3 (-1.09254848957061768f * y_36) * *_S3016);
    float fTmp0C_17 = -2.28522896766662598f * z2_17 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_35;
    float fC2_7 = x_38 * fC1_17 - y_36 * fS1_17;
    float fS2_7 = x_38 * fS1_17 + y_36 * fC1_17;
    float pSH12_15 = z_35 * (1.86588168144226074f * z2_17 - 1.11952900886535645f);
    float pSH13_13 = fTmp0C_17 * x_38;
    float pSH11_13 = fTmp0C_17 * y_36;
    float pSH14_13 = fTmp1B_17 * fC1_17;
    float pSH10_13 = fTmp1B_17 * fS1_17;
    float pSH15_13 = -0.59004360437393188f * fC2_7;
    float3  * _S3017 = v_coeffs_32 + int(8);
    *_S3017 = *_S3017 + make_float3 (-0.59004360437393188f * fS2_7) * v_colors_32;
    float3  * _S3018 = v_coeffs_32 + int(9);
    *_S3018 = *_S3018 + make_float3 (pSH10_13) * v_colors_32;
    float3  * _S3019 = v_coeffs_32 + int(10);
    *_S3019 = *_S3019 + make_float3 (pSH11_13) * v_colors_32;
    float3  * _S3020 = v_coeffs_32 + int(11);
    *_S3020 = *_S3020 + make_float3 (pSH12_15) * v_colors_32;
    float3  * _S3021 = v_coeffs_32 + int(12);
    *_S3021 = *_S3021 + make_float3 (pSH13_13) * v_colors_32;
    float3  * _S3022 = v_coeffs_32 + int(13);
    *_S3022 = *_S3022 + make_float3 (pSH14_13) * v_colors_32;
    float3  * _S3023 = v_coeffs_32 + int(14);
    *_S3023 = *_S3023 + make_float3 (pSH15_13) * v_colors_32;
    float fTmp0C_z_11 = -4.57045793533325195f * z_35;
    float _S3024 = x_38 * _S3007;
    float fC2_x_4 = fC1_17 + _S3024 - y_36 * fS1_x_11;
    float _S3025 = y_36 * _S3007;
    float fC2_y_4 = x_38 * fC1_y_11 - fS1_17 - _S3025;
    float fS2_x_4 = fS1_17 + x_38 * fS1_x_11 + _S3025;
    float fS2_y_4 = _S3024 + fC1_17 + y_36 * fC1_y_11;
    float pSH12_z_6 = 5.59764480590820312f * z2_17 - 1.11952900886535645f;
    float pSH14_x_11 = fTmp1B_17 * _S3007;
    float3  * _S3026 = coeffs_47 + int(8);
    float3  * _S3027 = coeffs_47 + int(14);
    float3  * _S3028 = coeffs_47 + int(9);
    float3  * _S3029 = coeffs_47 + int(13);
    float3  * _S3030 = coeffs_47 + int(12);
    float v_x_28 = v_x_27 + dot_0(v_colors_32, make_float3 (-0.59004360437393188f * fS2_x_4) * *_S3026 + make_float3 (-0.59004360437393188f * fC2_x_4) * *_S3027 + make_float3 (fTmp1B_17 * fS1_x_11) * *_S3028 + make_float3 (pSH14_x_11) * *_S3029 + make_float3 (fTmp0C_17) * *_S3030);
    float3  * _S3031 = coeffs_47 + int(10);
    float v_y_28 = v_y_27 + dot_0(v_colors_32, make_float3 (-0.59004360437393188f * fS2_y_4) * *_S3026 + make_float3 (-0.59004360437393188f * fC2_y_4) * *_S3027 + make_float3 (pSH14_x_11) * *_S3028 + make_float3 (fTmp1B_17 * fC1_y_11) * *_S3029 + make_float3 (fTmp0C_17) * *_S3031);
    float v_z_28 = v_z_27 + dot_0(v_colors_32, make_float3 (pSH12_z_6) * *(coeffs_47 + int(11)) + make_float3 (fTmp0C_z_11 * x_38) * *_S3030 + make_float3 (fTmp0C_z_11 * y_36) * *_S3031 + make_float3 (1.44530570507049561f * fC1_17) * *_S3029 + make_float3 (1.44530570507049561f * fS1_17) * *_S3028);
    float fTmp0D_7 = z_35 * (-4.68332576751708984f * z2_17 + 2.00713968276977539f);
    float fTmp1C_7 = 3.31161141395568848f * z2_17 - 0.47308734059333801f;
    float fTmp2B_7 = -1.77013075351715088f * z_35;
    float pSH20_7 = 1.9843134880065918f * z_35 * pSH12_15 + -1.00623059272766113f * pSH6_23;
    float pSH21_5 = fTmp0D_7 * x_38;
    float pSH19_5 = fTmp0D_7 * y_36;
    float pSH22_5 = fTmp1C_7 * fC1_17;
    float pSH18_5 = fTmp1C_7 * fS1_17;
    float pSH23_5 = fTmp2B_7 * fC2_7;
    float pSH17_5 = fTmp2B_7 * fS2_7;
    float pSH24_5 = 0.62583571672439575f * (x_38 * fC2_7 - y_36 * fS2_7);
    float3  * _S3032 = v_coeffs_32 + int(15);
    *_S3032 = *_S3032 + make_float3 (0.62583571672439575f * (x_38 * fS2_7 + y_36 * fC2_7)) * v_colors_32;
    float3  * _S3033 = v_coeffs_32 + int(16);
    *_S3033 = *_S3033 + make_float3 (pSH17_5) * v_colors_32;
    float3  * _S3034 = v_coeffs_32 + int(17);
    *_S3034 = *_S3034 + make_float3 (pSH18_5) * v_colors_32;
    float3  * _S3035 = v_coeffs_32 + int(18);
    *_S3035 = *_S3035 + make_float3 (pSH19_5) * v_colors_32;
    float3  * _S3036 = v_coeffs_32 + int(19);
    *_S3036 = *_S3036 + make_float3 (pSH20_7) * v_colors_32;
    float3  * _S3037 = v_coeffs_32 + int(20);
    *_S3037 = *_S3037 + make_float3 (pSH21_5) * v_colors_32;
    float3  * _S3038 = v_coeffs_32 + int(21);
    *_S3038 = *_S3038 + make_float3 (pSH22_5) * v_colors_32;
    float3  * _S3039 = v_coeffs_32 + int(22);
    *_S3039 = *_S3039 + make_float3 (pSH23_5) * v_colors_32;
    float3  * _S3040 = v_coeffs_32 + int(23);
    *_S3040 = *_S3040 + make_float3 (pSH24_5) * v_colors_32;
    float fTmp0D_z_4 = -14.04997730255126953f * z2_17 + 2.00713968276977539f;
    float fTmp1C_z_4 = 6.62322282791137695f * z_35;
    float pSH22_x_4 = fTmp1C_7 * _S3007;
    float3  * _S3041 = coeffs_47 + int(15);
    float3  * _S3042 = coeffs_47 + int(23);
    float3  * _S3043 = coeffs_47 + int(16);
    float3  * _S3044 = coeffs_47 + int(22);
    float3  * _S3045 = coeffs_47 + int(17);
    float3  * _S3046 = coeffs_47 + int(21);
    float3  * _S3047 = coeffs_47 + int(20);
    float3  * _S3048 = coeffs_47 + int(18);
    float3  dir_n_25 = make_float3 (x_38, y_36, z_35);
    float3  v_dir_n_41 = make_float3 (v_x_28 + dot_0(v_colors_32, make_float3 (0.62583571672439575f * (fS2_7 + y_36 * fC2_x_4 + x_38 * fS2_x_4)) * *_S3041 + make_float3 (0.62583571672439575f * (fC2_7 + x_38 * fC2_x_4 - y_36 * fS2_x_4)) * *_S3042 + make_float3 (fTmp2B_7 * fS2_x_4) * *_S3043 + make_float3 (fTmp2B_7 * fC2_x_4) * *_S3044 + make_float3 (fTmp1C_7 * fS1_x_11) * *_S3045 + make_float3 (pSH22_x_4) * *_S3046 + make_float3 (fTmp0D_7) * *_S3047), v_y_28 + dot_0(v_colors_32, make_float3 (0.62583571672439575f * (x_38 * fS2_y_4 + fC2_7 + y_36 * fC2_y_4)) * *_S3041 + make_float3 (0.62583571672439575f * (x_38 * fC2_y_4 - fS2_7 - y_36 * fS2_y_4)) * *_S3042 + make_float3 (fTmp2B_7 * fS2_y_4) * *_S3043 + make_float3 (fTmp2B_7 * fC2_y_4) * *_S3044 + make_float3 (pSH22_x_4) * *_S3045 + make_float3 (fTmp1C_7 * fC1_y_11) * *_S3046 + make_float3 (fTmp0D_7) * *_S3048), v_z_28 + dot_0(v_colors_32, make_float3 (1.9843134880065918f * (pSH12_15 + z_35 * pSH12_z_6) + -1.00623059272766113f * pSH6_z_8) * *(coeffs_47 + int(19)) + make_float3 (fTmp0D_z_4 * x_38) * *_S3047 + make_float3 (fTmp0D_z_4 * y_36) * *_S3048 + make_float3 (fTmp1C_z_4 * fC1_17) * *_S3046 + make_float3 (fTmp1C_z_4 * fS1_17) * *_S3045 + make_float3 (-1.77013075351715088f * fC2_7) * *_S3044 + make_float3 (-1.77013075351715088f * fS2_7) * *_S3043));
    *v_dir_12 = *v_dir_12 + (v_dir_n_41 - make_float3 (dot_0(v_dir_n_41, dir_n_25)) * dir_n_25) * make_float3 (inorm_9);
    return;
}

inline __device__ void sh4_to_color_dir_vjp_atomic(float3  dir_18, float3  coeff_dc_48, float3  * coeffs_48, float3  v_colors_33, float3  * v_coeff_dc_33, float3  * v_coeffs_33, float3  * v_dir_13)
{
    *v_coeff_dc_33 = *v_coeff_dc_33 + make_float3 (0.282094806432724f) * v_colors_33;
    float _S3049 = dir_18.x;
    float _S3050 = dir_18.y;
    float _S3051 = dir_18.z;
    float inorm_10 = (F32_rsqrt((_S3049 * _S3049 + _S3050 * _S3050 + _S3051 * _S3051)));
    float x_39 = _S3049 * inorm_10;
    float y_37 = _S3050 * inorm_10;
    float z_36 = _S3051 * inorm_10;
    float3  temp_226 = make_float3 (-0.48860251903533936f * y_37) * v_colors_33;
    float _S3052 = dot_0(temp_226, temp_226);
    bool _S3053;
    if((F32_isfinite((_S3052))))
    {
        _S3053 = _S3052 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3054 = v_coeffs_33 + int(0);
        float _S3055 = atomicAdd(&(_S3054->x), temp_226.x);
        float _S3056 = atomicAdd(&(_S3054->y), temp_226.y);
        float _S3057 = atomicAdd(&(_S3054->z), temp_226.z);
    }
    float3  temp_227 = make_float3 (0.48860251903533936f * z_36) * v_colors_33;
    float _S3058 = dot_0(temp_227, temp_227);
    if((F32_isfinite((_S3058))))
    {
        _S3053 = _S3058 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3059 = v_coeffs_33 + int(1);
        float _S3060 = atomicAdd(&(_S3059->x), temp_227.x);
        float _S3061 = atomicAdd(&(_S3059->y), temp_227.y);
        float _S3062 = atomicAdd(&(_S3059->z), temp_227.z);
    }
    float3  temp_228 = make_float3 (-0.48860251903533936f * x_39) * v_colors_33;
    float _S3063 = dot_0(temp_228, temp_228);
    if((F32_isfinite((_S3063))))
    {
        _S3053 = _S3063 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3064 = v_coeffs_33 + int(2);
        float _S3065 = atomicAdd(&(_S3064->x), temp_228.x);
        float _S3066 = atomicAdd(&(_S3064->y), temp_228.y);
        float _S3067 = atomicAdd(&(_S3064->z), temp_228.z);
    }
    float _S3068 = -0.48860251903533936f * dot_0(*(coeffs_48 + int(2)), v_colors_33);
    float _S3069 = -0.48860251903533936f * dot_0(*(coeffs_48 + int(0)), v_colors_33);
    float _S3070 = 0.48860251903533936f * dot_0(*(coeffs_48 + int(1)), v_colors_33);
    float z2_18 = z_36 * z_36;
    float fTmp0B_28 = -1.09254848957061768f * z_36;
    float fC1_18 = x_39 * x_39 - y_37 * y_37;
    float _S3071 = 2.0f * x_39;
    float fS1_18 = _S3071 * y_37;
    float pSH6_24 = 0.94617468118667603f * z2_18 - 0.31539157032966614f;
    float pSH7_22 = fTmp0B_28 * x_39;
    float pSH5_22 = fTmp0B_28 * y_37;
    float pSH8_22 = 0.54627424478530884f * fC1_18;
    float3  temp_229 = make_float3 (0.54627424478530884f * fS1_18) * v_colors_33;
    float _S3072 = dot_0(temp_229, temp_229);
    if((F32_isfinite((_S3072))))
    {
        _S3053 = _S3072 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3073 = v_coeffs_33 + int(3);
        float _S3074 = atomicAdd(&(_S3073->x), temp_229.x);
        float _S3075 = atomicAdd(&(_S3073->y), temp_229.y);
        float _S3076 = atomicAdd(&(_S3073->z), temp_229.z);
    }
    float3  temp_230 = make_float3 (pSH5_22) * v_colors_33;
    float _S3077 = dot_0(temp_230, temp_230);
    if((F32_isfinite((_S3077))))
    {
        _S3053 = _S3077 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3078 = v_coeffs_33 + int(4);
        float _S3079 = atomicAdd(&(_S3078->x), temp_230.x);
        float _S3080 = atomicAdd(&(_S3078->y), temp_230.y);
        float _S3081 = atomicAdd(&(_S3078->z), temp_230.z);
    }
    float3  temp_231 = make_float3 (pSH6_24) * v_colors_33;
    float _S3082 = dot_0(temp_231, temp_231);
    if((F32_isfinite((_S3082))))
    {
        _S3053 = _S3082 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3083 = v_coeffs_33 + int(5);
        float _S3084 = atomicAdd(&(_S3083->x), temp_231.x);
        float _S3085 = atomicAdd(&(_S3083->y), temp_231.y);
        float _S3086 = atomicAdd(&(_S3083->z), temp_231.z);
    }
    float3  temp_232 = make_float3 (pSH7_22) * v_colors_33;
    float _S3087 = dot_0(temp_232, temp_232);
    if((F32_isfinite((_S3087))))
    {
        _S3053 = _S3087 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3088 = v_coeffs_33 + int(6);
        float _S3089 = atomicAdd(&(_S3088->x), temp_232.x);
        float _S3090 = atomicAdd(&(_S3088->y), temp_232.y);
        float _S3091 = atomicAdd(&(_S3088->z), temp_232.z);
    }
    float3  temp_233 = make_float3 (pSH8_22) * v_colors_33;
    float _S3092 = dot_0(temp_233, temp_233);
    if((F32_isfinite((_S3092))))
    {
        _S3053 = _S3092 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3093 = v_coeffs_33 + int(7);
        float _S3094 = atomicAdd(&(_S3093->x), temp_233.x);
        float _S3095 = atomicAdd(&(_S3093->y), temp_233.y);
        float _S3096 = atomicAdd(&(_S3093->z), temp_233.z);
    }
    float fC1_y_12 = -2.0f * y_37;
    float fS1_x_12 = 2.0f * y_37;
    float pSH6_z_9 = 1.89234936237335205f * z_36;
    float pSH8_x_19 = 0.54627424478530884f * _S3071;
    float3  * _S3097 = coeffs_48 + int(3);
    float3  * _S3098 = coeffs_48 + int(7);
    float3  * _S3099 = coeffs_48 + int(6);
    float v_x_29 = _S3068 + dot_0(v_colors_33, make_float3 (0.54627424478530884f * fS1_x_12) * *_S3097 + make_float3 (pSH8_x_19) * *_S3098 + make_float3 (fTmp0B_28) * *_S3099);
    float3  * _S3100 = coeffs_48 + int(4);
    float v_y_29 = _S3069 + dot_0(v_colors_33, make_float3 (pSH8_x_19) * *_S3097 + make_float3 (0.54627424478530884f * fC1_y_12) * *_S3098 + make_float3 (fTmp0B_28) * *_S3100);
    float v_z_29 = _S3070 + dot_0(v_colors_33, make_float3 (pSH6_z_9) * *(coeffs_48 + int(5)) + make_float3 (-1.09254848957061768f * x_39) * *_S3099 + make_float3 (-1.09254848957061768f * y_37) * *_S3100);
    float fTmp0C_18 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_36;
    float fC2_8 = x_39 * fC1_18 - y_37 * fS1_18;
    float fS2_8 = x_39 * fS1_18 + y_37 * fC1_18;
    float pSH12_16 = z_36 * (1.86588168144226074f * z2_18 - 1.11952900886535645f);
    float pSH13_14 = fTmp0C_18 * x_39;
    float pSH11_14 = fTmp0C_18 * y_37;
    float pSH14_14 = fTmp1B_18 * fC1_18;
    float pSH10_14 = fTmp1B_18 * fS1_18;
    float pSH15_14 = -0.59004360437393188f * fC2_8;
    float3  temp_234 = make_float3 (-0.59004360437393188f * fS2_8) * v_colors_33;
    float _S3101 = dot_0(temp_234, temp_234);
    if((F32_isfinite((_S3101))))
    {
        _S3053 = _S3101 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3102 = v_coeffs_33 + int(8);
        float _S3103 = atomicAdd(&(_S3102->x), temp_234.x);
        float _S3104 = atomicAdd(&(_S3102->y), temp_234.y);
        float _S3105 = atomicAdd(&(_S3102->z), temp_234.z);
    }
    float3  temp_235 = make_float3 (pSH10_14) * v_colors_33;
    float _S3106 = dot_0(temp_235, temp_235);
    if((F32_isfinite((_S3106))))
    {
        _S3053 = _S3106 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3107 = v_coeffs_33 + int(9);
        float _S3108 = atomicAdd(&(_S3107->x), temp_235.x);
        float _S3109 = atomicAdd(&(_S3107->y), temp_235.y);
        float _S3110 = atomicAdd(&(_S3107->z), temp_235.z);
    }
    float3  temp_236 = make_float3 (pSH11_14) * v_colors_33;
    float _S3111 = dot_0(temp_236, temp_236);
    if((F32_isfinite((_S3111))))
    {
        _S3053 = _S3111 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3112 = v_coeffs_33 + int(10);
        float _S3113 = atomicAdd(&(_S3112->x), temp_236.x);
        float _S3114 = atomicAdd(&(_S3112->y), temp_236.y);
        float _S3115 = atomicAdd(&(_S3112->z), temp_236.z);
    }
    float3  temp_237 = make_float3 (pSH12_16) * v_colors_33;
    float _S3116 = dot_0(temp_237, temp_237);
    if((F32_isfinite((_S3116))))
    {
        _S3053 = _S3116 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3117 = v_coeffs_33 + int(11);
        float _S3118 = atomicAdd(&(_S3117->x), temp_237.x);
        float _S3119 = atomicAdd(&(_S3117->y), temp_237.y);
        float _S3120 = atomicAdd(&(_S3117->z), temp_237.z);
    }
    float3  temp_238 = make_float3 (pSH13_14) * v_colors_33;
    float _S3121 = dot_0(temp_238, temp_238);
    if((F32_isfinite((_S3121))))
    {
        _S3053 = _S3121 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3122 = v_coeffs_33 + int(12);
        float _S3123 = atomicAdd(&(_S3122->x), temp_238.x);
        float _S3124 = atomicAdd(&(_S3122->y), temp_238.y);
        float _S3125 = atomicAdd(&(_S3122->z), temp_238.z);
    }
    float3  temp_239 = make_float3 (pSH14_14) * v_colors_33;
    float _S3126 = dot_0(temp_239, temp_239);
    if((F32_isfinite((_S3126))))
    {
        _S3053 = _S3126 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3127 = v_coeffs_33 + int(13);
        float _S3128 = atomicAdd(&(_S3127->x), temp_239.x);
        float _S3129 = atomicAdd(&(_S3127->y), temp_239.y);
        float _S3130 = atomicAdd(&(_S3127->z), temp_239.z);
    }
    float3  temp_240 = make_float3 (pSH15_14) * v_colors_33;
    float _S3131 = dot_0(temp_240, temp_240);
    if((F32_isfinite((_S3131))))
    {
        _S3053 = _S3131 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3132 = v_coeffs_33 + int(14);
        float _S3133 = atomicAdd(&(_S3132->x), temp_240.x);
        float _S3134 = atomicAdd(&(_S3132->y), temp_240.y);
        float _S3135 = atomicAdd(&(_S3132->z), temp_240.z);
    }
    float fTmp0C_z_12 = -4.57045793533325195f * z_36;
    float _S3136 = x_39 * _S3071;
    float fC2_x_5 = fC1_18 + _S3136 - y_37 * fS1_x_12;
    float _S3137 = y_37 * _S3071;
    float fC2_y_5 = x_39 * fC1_y_12 - fS1_18 - _S3137;
    float fS2_x_5 = fS1_18 + x_39 * fS1_x_12 + _S3137;
    float fS2_y_5 = _S3136 + fC1_18 + y_37 * fC1_y_12;
    float pSH12_z_7 = 5.59764480590820312f * z2_18 - 1.11952900886535645f;
    float pSH14_x_12 = fTmp1B_18 * _S3071;
    float3  * _S3138 = coeffs_48 + int(8);
    float3  * _S3139 = coeffs_48 + int(14);
    float3  * _S3140 = coeffs_48 + int(9);
    float3  * _S3141 = coeffs_48 + int(13);
    float3  * _S3142 = coeffs_48 + int(12);
    float v_x_30 = v_x_29 + dot_0(v_colors_33, make_float3 (-0.59004360437393188f * fS2_x_5) * *_S3138 + make_float3 (-0.59004360437393188f * fC2_x_5) * *_S3139 + make_float3 (fTmp1B_18 * fS1_x_12) * *_S3140 + make_float3 (pSH14_x_12) * *_S3141 + make_float3 (fTmp0C_18) * *_S3142);
    float3  * _S3143 = coeffs_48 + int(10);
    float v_y_30 = v_y_29 + dot_0(v_colors_33, make_float3 (-0.59004360437393188f * fS2_y_5) * *_S3138 + make_float3 (-0.59004360437393188f * fC2_y_5) * *_S3139 + make_float3 (pSH14_x_12) * *_S3140 + make_float3 (fTmp1B_18 * fC1_y_12) * *_S3141 + make_float3 (fTmp0C_18) * *_S3143);
    float v_z_30 = v_z_29 + dot_0(v_colors_33, make_float3 (pSH12_z_7) * *(coeffs_48 + int(11)) + make_float3 (fTmp0C_z_12 * x_39) * *_S3142 + make_float3 (fTmp0C_z_12 * y_37) * *_S3143 + make_float3 (1.44530570507049561f * fC1_18) * *_S3141 + make_float3 (1.44530570507049561f * fS1_18) * *_S3140);
    float fTmp0D_8 = z_36 * (-4.68332576751708984f * z2_18 + 2.00713968276977539f);
    float fTmp1C_8 = 3.31161141395568848f * z2_18 - 0.47308734059333801f;
    float fTmp2B_8 = -1.77013075351715088f * z_36;
    float pSH20_8 = 1.9843134880065918f * z_36 * pSH12_16 + -1.00623059272766113f * pSH6_24;
    float pSH21_6 = fTmp0D_8 * x_39;
    float pSH19_6 = fTmp0D_8 * y_37;
    float pSH22_6 = fTmp1C_8 * fC1_18;
    float pSH18_6 = fTmp1C_8 * fS1_18;
    float pSH23_6 = fTmp2B_8 * fC2_8;
    float pSH17_6 = fTmp2B_8 * fS2_8;
    float pSH24_6 = 0.62583571672439575f * (x_39 * fC2_8 - y_37 * fS2_8);
    float3  temp_241 = make_float3 (0.62583571672439575f * (x_39 * fS2_8 + y_37 * fC2_8)) * v_colors_33;
    float _S3144 = dot_0(temp_241, temp_241);
    if((F32_isfinite((_S3144))))
    {
        _S3053 = _S3144 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3145 = v_coeffs_33 + int(15);
        float _S3146 = atomicAdd(&(_S3145->x), temp_241.x);
        float _S3147 = atomicAdd(&(_S3145->y), temp_241.y);
        float _S3148 = atomicAdd(&(_S3145->z), temp_241.z);
    }
    float3  temp_242 = make_float3 (pSH17_6) * v_colors_33;
    float _S3149 = dot_0(temp_242, temp_242);
    if((F32_isfinite((_S3149))))
    {
        _S3053 = _S3149 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3150 = v_coeffs_33 + int(16);
        float _S3151 = atomicAdd(&(_S3150->x), temp_242.x);
        float _S3152 = atomicAdd(&(_S3150->y), temp_242.y);
        float _S3153 = atomicAdd(&(_S3150->z), temp_242.z);
    }
    float3  temp_243 = make_float3 (pSH18_6) * v_colors_33;
    float _S3154 = dot_0(temp_243, temp_243);
    if((F32_isfinite((_S3154))))
    {
        _S3053 = _S3154 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3155 = v_coeffs_33 + int(17);
        float _S3156 = atomicAdd(&(_S3155->x), temp_243.x);
        float _S3157 = atomicAdd(&(_S3155->y), temp_243.y);
        float _S3158 = atomicAdd(&(_S3155->z), temp_243.z);
    }
    float3  temp_244 = make_float3 (pSH19_6) * v_colors_33;
    float _S3159 = dot_0(temp_244, temp_244);
    if((F32_isfinite((_S3159))))
    {
        _S3053 = _S3159 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3160 = v_coeffs_33 + int(18);
        float _S3161 = atomicAdd(&(_S3160->x), temp_244.x);
        float _S3162 = atomicAdd(&(_S3160->y), temp_244.y);
        float _S3163 = atomicAdd(&(_S3160->z), temp_244.z);
    }
    float3  temp_245 = make_float3 (pSH20_8) * v_colors_33;
    float _S3164 = dot_0(temp_245, temp_245);
    if((F32_isfinite((_S3164))))
    {
        _S3053 = _S3164 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3165 = v_coeffs_33 + int(19);
        float _S3166 = atomicAdd(&(_S3165->x), temp_245.x);
        float _S3167 = atomicAdd(&(_S3165->y), temp_245.y);
        float _S3168 = atomicAdd(&(_S3165->z), temp_245.z);
    }
    float3  temp_246 = make_float3 (pSH21_6) * v_colors_33;
    float _S3169 = dot_0(temp_246, temp_246);
    if((F32_isfinite((_S3169))))
    {
        _S3053 = _S3169 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3170 = v_coeffs_33 + int(20);
        float _S3171 = atomicAdd(&(_S3170->x), temp_246.x);
        float _S3172 = atomicAdd(&(_S3170->y), temp_246.y);
        float _S3173 = atomicAdd(&(_S3170->z), temp_246.z);
    }
    float3  temp_247 = make_float3 (pSH22_6) * v_colors_33;
    float _S3174 = dot_0(temp_247, temp_247);
    if((F32_isfinite((_S3174))))
    {
        _S3053 = _S3174 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3175 = v_coeffs_33 + int(21);
        float _S3176 = atomicAdd(&(_S3175->x), temp_247.x);
        float _S3177 = atomicAdd(&(_S3175->y), temp_247.y);
        float _S3178 = atomicAdd(&(_S3175->z), temp_247.z);
    }
    float3  temp_248 = make_float3 (pSH23_6) * v_colors_33;
    float _S3179 = dot_0(temp_248, temp_248);
    if((F32_isfinite((_S3179))))
    {
        _S3053 = _S3179 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3180 = v_coeffs_33 + int(22);
        float _S3181 = atomicAdd(&(_S3180->x), temp_248.x);
        float _S3182 = atomicAdd(&(_S3180->y), temp_248.y);
        float _S3183 = atomicAdd(&(_S3180->z), temp_248.z);
    }
    float3  temp_249 = make_float3 (pSH24_6) * v_colors_33;
    float _S3184 = dot_0(temp_249, temp_249);
    if((F32_isfinite((_S3184))))
    {
        _S3053 = _S3184 != 0.0f;
    }
    else
    {
        _S3053 = false;
    }
    if(_S3053)
    {
        float3  * _S3185 = v_coeffs_33 + int(23);
        float _S3186 = atomicAdd(&(_S3185->x), temp_249.x);
        float _S3187 = atomicAdd(&(_S3185->y), temp_249.y);
        float _S3188 = atomicAdd(&(_S3185->z), temp_249.z);
    }
    float fTmp0D_z_5 = -14.04997730255126953f * z2_18 + 2.00713968276977539f;
    float fTmp1C_z_5 = 6.62322282791137695f * z_36;
    float pSH22_x_5 = fTmp1C_8 * _S3071;
    float3  * _S3189 = coeffs_48 + int(15);
    float3  * _S3190 = coeffs_48 + int(23);
    float3  * _S3191 = coeffs_48 + int(16);
    float3  * _S3192 = coeffs_48 + int(22);
    float3  * _S3193 = coeffs_48 + int(17);
    float3  * _S3194 = coeffs_48 + int(21);
    float3  * _S3195 = coeffs_48 + int(20);
    float3  * _S3196 = coeffs_48 + int(18);
    float3  dir_n_26 = make_float3 (x_39, y_37, z_36);
    float3  v_dir_n_42 = make_float3 (v_x_30 + dot_0(v_colors_33, make_float3 (0.62583571672439575f * (fS2_8 + y_37 * fC2_x_5 + x_39 * fS2_x_5)) * *_S3189 + make_float3 (0.62583571672439575f * (fC2_8 + x_39 * fC2_x_5 - y_37 * fS2_x_5)) * *_S3190 + make_float3 (fTmp2B_8 * fS2_x_5) * *_S3191 + make_float3 (fTmp2B_8 * fC2_x_5) * *_S3192 + make_float3 (fTmp1C_8 * fS1_x_12) * *_S3193 + make_float3 (pSH22_x_5) * *_S3194 + make_float3 (fTmp0D_8) * *_S3195), v_y_30 + dot_0(v_colors_33, make_float3 (0.62583571672439575f * (x_39 * fS2_y_5 + fC2_8 + y_37 * fC2_y_5)) * *_S3189 + make_float3 (0.62583571672439575f * (x_39 * fC2_y_5 - fS2_8 - y_37 * fS2_y_5)) * *_S3190 + make_float3 (fTmp2B_8 * fS2_y_5) * *_S3191 + make_float3 (fTmp2B_8 * fC2_y_5) * *_S3192 + make_float3 (pSH22_x_5) * *_S3193 + make_float3 (fTmp1C_8 * fC1_y_12) * *_S3194 + make_float3 (fTmp0D_8) * *_S3196), v_z_30 + dot_0(v_colors_33, make_float3 (1.9843134880065918f * (pSH12_16 + z_36 * pSH12_z_7) + -1.00623059272766113f * pSH6_z_9) * *(coeffs_48 + int(19)) + make_float3 (fTmp0D_z_5 * x_39) * *_S3195 + make_float3 (fTmp0D_z_5 * y_37) * *_S3196 + make_float3 (fTmp1C_z_5 * fC1_18) * *_S3194 + make_float3 (fTmp1C_z_5 * fS1_18) * *_S3193 + make_float3 (-1.77013075351715088f * fC2_8) * *_S3192 + make_float3 (-1.77013075351715088f * fS2_8) * *_S3191));
    *v_dir_13 = *v_dir_13 + (v_dir_n_42 - make_float3 (dot_0(v_dir_n_42, dir_n_26)) * dir_n_26) * make_float3 (inorm_10);
    return;
}

inline __device__ void sh4_to_color_dir_vjp_block_atomic(float3  dir_19, float3  coeff_dc_49, float3  * coeffs_49, float3  v_colors_34, float3  * v_coeff_dc_34, float3  * v_coeffs_34, float3  * v_dir_14, uint thread_id_4, uint _S3197)
{
    *v_coeff_dc_34 = *v_coeff_dc_34 + make_float3 (0.282094806432724f) * v_colors_34;
    float _S3198 = dir_19.x;
    float _S3199 = dir_19.y;
    float _S3200 = dir_19.z;
    float inorm_11 = (F32_rsqrt((_S3198 * _S3198 + _S3199 * _S3199 + _S3200 * _S3200)));
    float x_40 = _S3198 * inorm_11;
    float y_38 = _S3199 * inorm_11;
    float z_37 = _S3200 * inorm_11;
    float3  _S3201 = make_float3 (-0.48860251903533936f * y_38) * v_colors_34;
    float3  _S3202 = _S3201;
    bool _S3203 = (F32_isfinite((_S3201.x)));
    uint _S3204 = __ballot_sync(_S3197, _S3203);
    float v_3;
    uint _S3205;
    if(_S3203)
    {
        float _S3206 = _S3202.x;
        uint _S3207 = __ballot_sync(_S3197, true);
        v_3 = _S3206;
        _S3205 = _S3207;
    }
    else
    {
        uint _S3208 = __ballot_sync(_S3197, true);
        v_3 = 0.0f;
        _S3205 = _S3208;
    }
    *&((&_S3202)->x) = v_3;
    bool _S3209 = (F32_isfinite((_S3202.y)));
    uint _S3210 = __ballot_sync(_S3205, _S3209);
    if(_S3209)
    {
        float _S3211 = _S3202.y;
        uint _S3212 = __ballot_sync(_S3205, true);
        v_3 = _S3211;
        _S3205 = _S3212;
    }
    else
    {
        uint _S3213 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3213;
    }
    *&((&_S3202)->y) = v_3;
    bool _S3214 = (F32_isfinite((_S3202.z)));
    uint _S3215 = __ballot_sync(_S3205, _S3214);
    if(_S3214)
    {
        float _S3216 = _S3202.z;
        uint _S3217 = __ballot_sync(_S3205, true);
        v_3 = _S3216;
        _S3205 = _S3217;
    }
    else
    {
        uint _S3218 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3218;
    }
    *&((&_S3202)->z) = v_3;
    float _S3219 = WaveActiveSum_0(_S3202.x, _S3205);
    *&((&_S3202)->x) = _S3219;
    float _S3220 = WaveActiveSum_0(_S3202.y, _S3205);
    *&((&_S3202)->y) = _S3220;
    float _S3221 = WaveActiveSum_0(_S3202.z, _S3205);
    *&((&_S3202)->z) = _S3221;
    uint warp_id_26 = thread_id_4 / 32U;
    uint lane_id_3 = thread_id_4 % 32U;
    bool _S3222 = lane_id_3 == 0U;
    uint _S3223 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_26] = _S3202.x;
        (*&_sh_block_reduce_shared_0)[warp_id_26 + 16U] = _S3202.y;
        (*&_sh_block_reduce_shared_0)[warp_id_26 + 32U] = _S3202.z;
        uint _S3224 = __ballot_sync(_S3205, true);
        _S3205 = _S3224;
    }
    else
    {
        uint _S3225 = __ballot_sync(_S3205, true);
        _S3205 = _S3225;
    }
    __syncthreads();
    bool _S3226 = warp_id_26 < 3U;
    uint _S3227 = __ballot_sync(_S3205, _S3226);
    uint _S3228;
    bool _S3229;
    if(_S3226)
    {
        bool _S3230 = lane_id_3 < 16U;
        uint _S3231 = __ballot_sync(_S3227, _S3230);
        if(_S3230)
        {
            float _S3232 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_26 * 16U];
            uint _S3233 = __ballot_sync(_S3227, true);
            v_3 = _S3232;
            _S3228 = _S3233;
        }
        else
        {
            uint _S3234 = __ballot_sync(_S3227, true);
            v_3 = 0.0f;
            _S3228 = _S3234;
        }
        float _S3235 = WaveActiveSum_0(v_3, _S3228);
        uint _S3236 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3237 = _S3235 != 0.0f;
            uint _S3238 = __ballot_sync(_S3228, true);
            _S3229 = _S3237;
            _S3228 = _S3238;
        }
        else
        {
            uint _S3239 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3239;
        }
        uint _S3240 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3241 = warp_id_26 == 0U;
            uint _S3242 = __ballot_sync(_S3240, _S3241);
            if(_S3241)
            {
                float _S3243 = atomicAdd(&((v_coeffs_34 + 0U)->x), _S3235);
                uint _S3244 = __ballot_sync(_S3240, true);
            }
            else
            {
                uint _S3245 = _S3240 & (~_S3242);
                bool _S3246 = warp_id_26 == 1U;
                uint _S3247 = __ballot_sync(_S3245, _S3246);
                if(_S3246)
                {
                    float _S3248 = atomicAdd(&((v_coeffs_34 + 0U)->y), _S3235);
                    uint _S3249 = __ballot_sync(_S3245, true);
                }
                else
                {
                    float _S3250 = atomicAdd(&((v_coeffs_34 + 0U)->z), _S3235);
                    uint _S3251 = __ballot_sync(_S3245, true);
                }
                uint _S3252 = __ballot_sync(_S3240, true);
            }
            uint _S3253 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3254 = __ballot_sync(_S3228, true);
        }
        uint _S3255 = __ballot_sync(_S3205, true);
        _S3205 = _S3255;
    }
    else
    {
        uint _S3256 = __ballot_sync(_S3205, true);
        _S3205 = _S3256;
    }
    float3  _S3257 = make_float3 (0.48860251903533936f * z_37) * v_colors_34;
    float3  _S3258 = _S3257;
    bool _S3259 = (F32_isfinite((_S3257.x)));
    uint _S3260 = __ballot_sync(_S3205, _S3259);
    if(_S3259)
    {
        float _S3261 = _S3258.x;
        uint _S3262 = __ballot_sync(_S3205, true);
        v_3 = _S3261;
        _S3205 = _S3262;
    }
    else
    {
        uint _S3263 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3263;
    }
    *&((&_S3258)->x) = v_3;
    bool _S3264 = (F32_isfinite((_S3258.y)));
    uint _S3265 = __ballot_sync(_S3205, _S3264);
    if(_S3264)
    {
        float _S3266 = _S3258.y;
        uint _S3267 = __ballot_sync(_S3205, true);
        v_3 = _S3266;
        _S3205 = _S3267;
    }
    else
    {
        uint _S3268 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3268;
    }
    *&((&_S3258)->y) = v_3;
    bool _S3269 = (F32_isfinite((_S3258.z)));
    uint _S3270 = __ballot_sync(_S3205, _S3269);
    if(_S3269)
    {
        float _S3271 = _S3258.z;
        uint _S3272 = __ballot_sync(_S3205, true);
        v_3 = _S3271;
        _S3205 = _S3272;
    }
    else
    {
        uint _S3273 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3273;
    }
    *&((&_S3258)->z) = v_3;
    float _S3274 = WaveActiveSum_0(_S3258.x, _S3205);
    *&((&_S3258)->x) = _S3274;
    float _S3275 = WaveActiveSum_0(_S3258.y, _S3205);
    *&((&_S3258)->y) = _S3275;
    float _S3276 = WaveActiveSum_0(_S3258.z, _S3205);
    *&((&_S3258)->z) = _S3276;
    uint warp_id_27 = thread_id_4 / 32U;
    uint _S3277 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_27] = _S3258.x;
        (*&_sh_block_reduce_shared_0)[warp_id_27 + 16U] = _S3258.y;
        (*&_sh_block_reduce_shared_0)[warp_id_27 + 32U] = _S3258.z;
        uint _S3278 = __ballot_sync(_S3205, true);
        _S3205 = _S3278;
    }
    else
    {
        uint _S3279 = __ballot_sync(_S3205, true);
        _S3205 = _S3279;
    }
    __syncthreads();
    bool _S3280 = warp_id_27 < 3U;
    uint _S3281 = __ballot_sync(_S3205, _S3280);
    if(_S3280)
    {
        bool _S3282 = lane_id_3 < 16U;
        uint _S3283 = __ballot_sync(_S3281, _S3282);
        if(_S3282)
        {
            float _S3284 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_27 * 16U];
            uint _S3285 = __ballot_sync(_S3281, true);
            v_3 = _S3284;
            _S3228 = _S3285;
        }
        else
        {
            uint _S3286 = __ballot_sync(_S3281, true);
            v_3 = 0.0f;
            _S3228 = _S3286;
        }
        float _S3287 = WaveActiveSum_0(v_3, _S3228);
        uint _S3288 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3289 = _S3287 != 0.0f;
            uint _S3290 = __ballot_sync(_S3228, true);
            _S3229 = _S3289;
            _S3228 = _S3290;
        }
        else
        {
            uint _S3291 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3291;
        }
        uint _S3292 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3293 = warp_id_27 == 0U;
            uint _S3294 = __ballot_sync(_S3292, _S3293);
            if(_S3293)
            {
                float _S3295 = atomicAdd(&((v_coeffs_34 + 1U)->x), _S3287);
                uint _S3296 = __ballot_sync(_S3292, true);
            }
            else
            {
                uint _S3297 = _S3292 & (~_S3294);
                bool _S3298 = warp_id_27 == 1U;
                uint _S3299 = __ballot_sync(_S3297, _S3298);
                if(_S3298)
                {
                    float _S3300 = atomicAdd(&((v_coeffs_34 + 1U)->y), _S3287);
                    uint _S3301 = __ballot_sync(_S3297, true);
                }
                else
                {
                    float _S3302 = atomicAdd(&((v_coeffs_34 + 1U)->z), _S3287);
                    uint _S3303 = __ballot_sync(_S3297, true);
                }
                uint _S3304 = __ballot_sync(_S3292, true);
            }
            uint _S3305 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3306 = __ballot_sync(_S3228, true);
        }
        uint _S3307 = __ballot_sync(_S3205, true);
        _S3205 = _S3307;
    }
    else
    {
        uint _S3308 = __ballot_sync(_S3205, true);
        _S3205 = _S3308;
    }
    float3  _S3309 = make_float3 (-0.48860251903533936f * x_40) * v_colors_34;
    float3  _S3310 = _S3309;
    bool _S3311 = (F32_isfinite((_S3309.x)));
    uint _S3312 = __ballot_sync(_S3205, _S3311);
    if(_S3311)
    {
        float _S3313 = _S3310.x;
        uint _S3314 = __ballot_sync(_S3205, true);
        v_3 = _S3313;
        _S3205 = _S3314;
    }
    else
    {
        uint _S3315 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3315;
    }
    *&((&_S3310)->x) = v_3;
    bool _S3316 = (F32_isfinite((_S3310.y)));
    uint _S3317 = __ballot_sync(_S3205, _S3316);
    if(_S3316)
    {
        float _S3318 = _S3310.y;
        uint _S3319 = __ballot_sync(_S3205, true);
        v_3 = _S3318;
        _S3205 = _S3319;
    }
    else
    {
        uint _S3320 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3320;
    }
    *&((&_S3310)->y) = v_3;
    bool _S3321 = (F32_isfinite((_S3310.z)));
    uint _S3322 = __ballot_sync(_S3205, _S3321);
    if(_S3321)
    {
        float _S3323 = _S3310.z;
        uint _S3324 = __ballot_sync(_S3205, true);
        v_3 = _S3323;
        _S3205 = _S3324;
    }
    else
    {
        uint _S3325 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3325;
    }
    *&((&_S3310)->z) = v_3;
    float _S3326 = WaveActiveSum_0(_S3310.x, _S3205);
    *&((&_S3310)->x) = _S3326;
    float _S3327 = WaveActiveSum_0(_S3310.y, _S3205);
    *&((&_S3310)->y) = _S3327;
    float _S3328 = WaveActiveSum_0(_S3310.z, _S3205);
    *&((&_S3310)->z) = _S3328;
    uint warp_id_28 = thread_id_4 / 32U;
    uint _S3329 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_28] = _S3310.x;
        (*&_sh_block_reduce_shared_0)[warp_id_28 + 16U] = _S3310.y;
        (*&_sh_block_reduce_shared_0)[warp_id_28 + 32U] = _S3310.z;
        uint _S3330 = __ballot_sync(_S3205, true);
        _S3205 = _S3330;
    }
    else
    {
        uint _S3331 = __ballot_sync(_S3205, true);
        _S3205 = _S3331;
    }
    __syncthreads();
    bool _S3332 = warp_id_28 < 3U;
    uint _S3333 = __ballot_sync(_S3205, _S3332);
    if(_S3332)
    {
        bool _S3334 = lane_id_3 < 16U;
        uint _S3335 = __ballot_sync(_S3333, _S3334);
        if(_S3334)
        {
            float _S3336 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_28 * 16U];
            uint _S3337 = __ballot_sync(_S3333, true);
            v_3 = _S3336;
            _S3228 = _S3337;
        }
        else
        {
            uint _S3338 = __ballot_sync(_S3333, true);
            v_3 = 0.0f;
            _S3228 = _S3338;
        }
        float _S3339 = WaveActiveSum_0(v_3, _S3228);
        uint _S3340 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3341 = _S3339 != 0.0f;
            uint _S3342 = __ballot_sync(_S3228, true);
            _S3229 = _S3341;
            _S3228 = _S3342;
        }
        else
        {
            uint _S3343 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3343;
        }
        uint _S3344 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3345 = warp_id_28 == 0U;
            uint _S3346 = __ballot_sync(_S3344, _S3345);
            if(_S3345)
            {
                float _S3347 = atomicAdd(&((v_coeffs_34 + 2U)->x), _S3339);
                uint _S3348 = __ballot_sync(_S3344, true);
            }
            else
            {
                uint _S3349 = _S3344 & (~_S3346);
                bool _S3350 = warp_id_28 == 1U;
                uint _S3351 = __ballot_sync(_S3349, _S3350);
                if(_S3350)
                {
                    float _S3352 = atomicAdd(&((v_coeffs_34 + 2U)->y), _S3339);
                    uint _S3353 = __ballot_sync(_S3349, true);
                }
                else
                {
                    float _S3354 = atomicAdd(&((v_coeffs_34 + 2U)->z), _S3339);
                    uint _S3355 = __ballot_sync(_S3349, true);
                }
                uint _S3356 = __ballot_sync(_S3344, true);
            }
            uint _S3357 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3358 = __ballot_sync(_S3228, true);
        }
        uint _S3359 = __ballot_sync(_S3205, true);
        _S3205 = _S3359;
    }
    else
    {
        uint _S3360 = __ballot_sync(_S3205, true);
        _S3205 = _S3360;
    }
    float _S3361 = -0.48860251903533936f * dot_0(*(coeffs_49 + int(2)), v_colors_34);
    float _S3362 = -0.48860251903533936f * dot_0(*(coeffs_49 + int(0)), v_colors_34);
    float _S3363 = 0.48860251903533936f * dot_0(*(coeffs_49 + int(1)), v_colors_34);
    float z2_19 = z_37 * z_37;
    float fTmp0B_29 = -1.09254848957061768f * z_37;
    float fC1_19 = x_40 * x_40 - y_38 * y_38;
    float _S3364 = 2.0f * x_40;
    float fS1_19 = _S3364 * y_38;
    float pSH6_25 = 0.94617468118667603f * z2_19 - 0.31539157032966614f;
    float pSH7_23 = fTmp0B_29 * x_40;
    float pSH5_23 = fTmp0B_29 * y_38;
    float pSH8_23 = 0.54627424478530884f * fC1_19;
    float3  _S3365 = make_float3 (0.54627424478530884f * fS1_19) * v_colors_34;
    float3  _S3366 = _S3365;
    bool _S3367 = (F32_isfinite((_S3365.x)));
    uint _S3368 = __ballot_sync(_S3205, _S3367);
    if(_S3367)
    {
        float _S3369 = _S3366.x;
        uint _S3370 = __ballot_sync(_S3205, true);
        v_3 = _S3369;
        _S3205 = _S3370;
    }
    else
    {
        uint _S3371 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3371;
    }
    *&((&_S3366)->x) = v_3;
    bool _S3372 = (F32_isfinite((_S3366.y)));
    uint _S3373 = __ballot_sync(_S3205, _S3372);
    if(_S3372)
    {
        float _S3374 = _S3366.y;
        uint _S3375 = __ballot_sync(_S3205, true);
        v_3 = _S3374;
        _S3205 = _S3375;
    }
    else
    {
        uint _S3376 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3376;
    }
    *&((&_S3366)->y) = v_3;
    bool _S3377 = (F32_isfinite((_S3366.z)));
    uint _S3378 = __ballot_sync(_S3205, _S3377);
    if(_S3377)
    {
        float _S3379 = _S3366.z;
        uint _S3380 = __ballot_sync(_S3205, true);
        v_3 = _S3379;
        _S3205 = _S3380;
    }
    else
    {
        uint _S3381 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3381;
    }
    *&((&_S3366)->z) = v_3;
    float _S3382 = WaveActiveSum_0(_S3366.x, _S3205);
    *&((&_S3366)->x) = _S3382;
    float _S3383 = WaveActiveSum_0(_S3366.y, _S3205);
    *&((&_S3366)->y) = _S3383;
    float _S3384 = WaveActiveSum_0(_S3366.z, _S3205);
    *&((&_S3366)->z) = _S3384;
    uint warp_id_29 = thread_id_4 / 32U;
    uint _S3385 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_29] = _S3366.x;
        (*&_sh_block_reduce_shared_0)[warp_id_29 + 16U] = _S3366.y;
        (*&_sh_block_reduce_shared_0)[warp_id_29 + 32U] = _S3366.z;
        uint _S3386 = __ballot_sync(_S3205, true);
        _S3205 = _S3386;
    }
    else
    {
        uint _S3387 = __ballot_sync(_S3205, true);
        _S3205 = _S3387;
    }
    __syncthreads();
    bool _S3388 = warp_id_29 < 3U;
    uint _S3389 = __ballot_sync(_S3205, _S3388);
    if(_S3388)
    {
        bool _S3390 = lane_id_3 < 16U;
        uint _S3391 = __ballot_sync(_S3389, _S3390);
        if(_S3390)
        {
            float _S3392 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_29 * 16U];
            uint _S3393 = __ballot_sync(_S3389, true);
            v_3 = _S3392;
            _S3228 = _S3393;
        }
        else
        {
            uint _S3394 = __ballot_sync(_S3389, true);
            v_3 = 0.0f;
            _S3228 = _S3394;
        }
        float _S3395 = WaveActiveSum_0(v_3, _S3228);
        uint _S3396 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3397 = _S3395 != 0.0f;
            uint _S3398 = __ballot_sync(_S3228, true);
            _S3229 = _S3397;
            _S3228 = _S3398;
        }
        else
        {
            uint _S3399 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3399;
        }
        uint _S3400 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3401 = warp_id_29 == 0U;
            uint _S3402 = __ballot_sync(_S3400, _S3401);
            if(_S3401)
            {
                float _S3403 = atomicAdd(&((v_coeffs_34 + 3U)->x), _S3395);
                uint _S3404 = __ballot_sync(_S3400, true);
            }
            else
            {
                uint _S3405 = _S3400 & (~_S3402);
                bool _S3406 = warp_id_29 == 1U;
                uint _S3407 = __ballot_sync(_S3405, _S3406);
                if(_S3406)
                {
                    float _S3408 = atomicAdd(&((v_coeffs_34 + 3U)->y), _S3395);
                    uint _S3409 = __ballot_sync(_S3405, true);
                }
                else
                {
                    float _S3410 = atomicAdd(&((v_coeffs_34 + 3U)->z), _S3395);
                    uint _S3411 = __ballot_sync(_S3405, true);
                }
                uint _S3412 = __ballot_sync(_S3400, true);
            }
            uint _S3413 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3414 = __ballot_sync(_S3228, true);
        }
        uint _S3415 = __ballot_sync(_S3205, true);
        _S3205 = _S3415;
    }
    else
    {
        uint _S3416 = __ballot_sync(_S3205, true);
        _S3205 = _S3416;
    }
    float3  _S3417 = make_float3 (pSH5_23) * v_colors_34;
    float3  _S3418 = _S3417;
    bool _S3419 = (F32_isfinite((_S3417.x)));
    uint _S3420 = __ballot_sync(_S3205, _S3419);
    if(_S3419)
    {
        float _S3421 = _S3418.x;
        uint _S3422 = __ballot_sync(_S3205, true);
        v_3 = _S3421;
        _S3205 = _S3422;
    }
    else
    {
        uint _S3423 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3423;
    }
    *&((&_S3418)->x) = v_3;
    bool _S3424 = (F32_isfinite((_S3418.y)));
    uint _S3425 = __ballot_sync(_S3205, _S3424);
    if(_S3424)
    {
        float _S3426 = _S3418.y;
        uint _S3427 = __ballot_sync(_S3205, true);
        v_3 = _S3426;
        _S3205 = _S3427;
    }
    else
    {
        uint _S3428 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3428;
    }
    *&((&_S3418)->y) = v_3;
    bool _S3429 = (F32_isfinite((_S3418.z)));
    uint _S3430 = __ballot_sync(_S3205, _S3429);
    if(_S3429)
    {
        float _S3431 = _S3418.z;
        uint _S3432 = __ballot_sync(_S3205, true);
        v_3 = _S3431;
        _S3205 = _S3432;
    }
    else
    {
        uint _S3433 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3433;
    }
    *&((&_S3418)->z) = v_3;
    float _S3434 = WaveActiveSum_0(_S3418.x, _S3205);
    *&((&_S3418)->x) = _S3434;
    float _S3435 = WaveActiveSum_0(_S3418.y, _S3205);
    *&((&_S3418)->y) = _S3435;
    float _S3436 = WaveActiveSum_0(_S3418.z, _S3205);
    *&((&_S3418)->z) = _S3436;
    uint warp_id_30 = thread_id_4 / 32U;
    uint _S3437 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_30] = _S3418.x;
        (*&_sh_block_reduce_shared_0)[warp_id_30 + 16U] = _S3418.y;
        (*&_sh_block_reduce_shared_0)[warp_id_30 + 32U] = _S3418.z;
        uint _S3438 = __ballot_sync(_S3205, true);
        _S3205 = _S3438;
    }
    else
    {
        uint _S3439 = __ballot_sync(_S3205, true);
        _S3205 = _S3439;
    }
    __syncthreads();
    bool _S3440 = warp_id_30 < 3U;
    uint _S3441 = __ballot_sync(_S3205, _S3440);
    if(_S3440)
    {
        bool _S3442 = lane_id_3 < 16U;
        uint _S3443 = __ballot_sync(_S3441, _S3442);
        if(_S3442)
        {
            float _S3444 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_30 * 16U];
            uint _S3445 = __ballot_sync(_S3441, true);
            v_3 = _S3444;
            _S3228 = _S3445;
        }
        else
        {
            uint _S3446 = __ballot_sync(_S3441, true);
            v_3 = 0.0f;
            _S3228 = _S3446;
        }
        float _S3447 = WaveActiveSum_0(v_3, _S3228);
        uint _S3448 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3449 = _S3447 != 0.0f;
            uint _S3450 = __ballot_sync(_S3228, true);
            _S3229 = _S3449;
            _S3228 = _S3450;
        }
        else
        {
            uint _S3451 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3451;
        }
        uint _S3452 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3453 = warp_id_30 == 0U;
            uint _S3454 = __ballot_sync(_S3452, _S3453);
            if(_S3453)
            {
                float _S3455 = atomicAdd(&((v_coeffs_34 + 4U)->x), _S3447);
                uint _S3456 = __ballot_sync(_S3452, true);
            }
            else
            {
                uint _S3457 = _S3452 & (~_S3454);
                bool _S3458 = warp_id_30 == 1U;
                uint _S3459 = __ballot_sync(_S3457, _S3458);
                if(_S3458)
                {
                    float _S3460 = atomicAdd(&((v_coeffs_34 + 4U)->y), _S3447);
                    uint _S3461 = __ballot_sync(_S3457, true);
                }
                else
                {
                    float _S3462 = atomicAdd(&((v_coeffs_34 + 4U)->z), _S3447);
                    uint _S3463 = __ballot_sync(_S3457, true);
                }
                uint _S3464 = __ballot_sync(_S3452, true);
            }
            uint _S3465 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3466 = __ballot_sync(_S3228, true);
        }
        uint _S3467 = __ballot_sync(_S3205, true);
        _S3205 = _S3467;
    }
    else
    {
        uint _S3468 = __ballot_sync(_S3205, true);
        _S3205 = _S3468;
    }
    float3  _S3469 = make_float3 (pSH6_25) * v_colors_34;
    float3  _S3470 = _S3469;
    bool _S3471 = (F32_isfinite((_S3469.x)));
    uint _S3472 = __ballot_sync(_S3205, _S3471);
    if(_S3471)
    {
        float _S3473 = _S3470.x;
        uint _S3474 = __ballot_sync(_S3205, true);
        v_3 = _S3473;
        _S3205 = _S3474;
    }
    else
    {
        uint _S3475 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3475;
    }
    *&((&_S3470)->x) = v_3;
    bool _S3476 = (F32_isfinite((_S3470.y)));
    uint _S3477 = __ballot_sync(_S3205, _S3476);
    if(_S3476)
    {
        float _S3478 = _S3470.y;
        uint _S3479 = __ballot_sync(_S3205, true);
        v_3 = _S3478;
        _S3205 = _S3479;
    }
    else
    {
        uint _S3480 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3480;
    }
    *&((&_S3470)->y) = v_3;
    bool _S3481 = (F32_isfinite((_S3470.z)));
    uint _S3482 = __ballot_sync(_S3205, _S3481);
    if(_S3481)
    {
        float _S3483 = _S3470.z;
        uint _S3484 = __ballot_sync(_S3205, true);
        v_3 = _S3483;
        _S3205 = _S3484;
    }
    else
    {
        uint _S3485 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3485;
    }
    *&((&_S3470)->z) = v_3;
    float _S3486 = WaveActiveSum_0(_S3470.x, _S3205);
    *&((&_S3470)->x) = _S3486;
    float _S3487 = WaveActiveSum_0(_S3470.y, _S3205);
    *&((&_S3470)->y) = _S3487;
    float _S3488 = WaveActiveSum_0(_S3470.z, _S3205);
    *&((&_S3470)->z) = _S3488;
    uint warp_id_31 = thread_id_4 / 32U;
    uint _S3489 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_31] = _S3470.x;
        (*&_sh_block_reduce_shared_0)[warp_id_31 + 16U] = _S3470.y;
        (*&_sh_block_reduce_shared_0)[warp_id_31 + 32U] = _S3470.z;
        uint _S3490 = __ballot_sync(_S3205, true);
        _S3205 = _S3490;
    }
    else
    {
        uint _S3491 = __ballot_sync(_S3205, true);
        _S3205 = _S3491;
    }
    __syncthreads();
    bool _S3492 = warp_id_31 < 3U;
    uint _S3493 = __ballot_sync(_S3205, _S3492);
    if(_S3492)
    {
        bool _S3494 = lane_id_3 < 16U;
        uint _S3495 = __ballot_sync(_S3493, _S3494);
        if(_S3494)
        {
            float _S3496 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_31 * 16U];
            uint _S3497 = __ballot_sync(_S3493, true);
            v_3 = _S3496;
            _S3228 = _S3497;
        }
        else
        {
            uint _S3498 = __ballot_sync(_S3493, true);
            v_3 = 0.0f;
            _S3228 = _S3498;
        }
        float _S3499 = WaveActiveSum_0(v_3, _S3228);
        uint _S3500 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3501 = _S3499 != 0.0f;
            uint _S3502 = __ballot_sync(_S3228, true);
            _S3229 = _S3501;
            _S3228 = _S3502;
        }
        else
        {
            uint _S3503 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3503;
        }
        uint _S3504 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3505 = warp_id_31 == 0U;
            uint _S3506 = __ballot_sync(_S3504, _S3505);
            if(_S3505)
            {
                float _S3507 = atomicAdd(&((v_coeffs_34 + 5U)->x), _S3499);
                uint _S3508 = __ballot_sync(_S3504, true);
            }
            else
            {
                uint _S3509 = _S3504 & (~_S3506);
                bool _S3510 = warp_id_31 == 1U;
                uint _S3511 = __ballot_sync(_S3509, _S3510);
                if(_S3510)
                {
                    float _S3512 = atomicAdd(&((v_coeffs_34 + 5U)->y), _S3499);
                    uint _S3513 = __ballot_sync(_S3509, true);
                }
                else
                {
                    float _S3514 = atomicAdd(&((v_coeffs_34 + 5U)->z), _S3499);
                    uint _S3515 = __ballot_sync(_S3509, true);
                }
                uint _S3516 = __ballot_sync(_S3504, true);
            }
            uint _S3517 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3518 = __ballot_sync(_S3228, true);
        }
        uint _S3519 = __ballot_sync(_S3205, true);
        _S3205 = _S3519;
    }
    else
    {
        uint _S3520 = __ballot_sync(_S3205, true);
        _S3205 = _S3520;
    }
    float3  _S3521 = make_float3 (pSH7_23) * v_colors_34;
    float3  _S3522 = _S3521;
    bool _S3523 = (F32_isfinite((_S3521.x)));
    uint _S3524 = __ballot_sync(_S3205, _S3523);
    if(_S3523)
    {
        float _S3525 = _S3522.x;
        uint _S3526 = __ballot_sync(_S3205, true);
        v_3 = _S3525;
        _S3205 = _S3526;
    }
    else
    {
        uint _S3527 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3527;
    }
    *&((&_S3522)->x) = v_3;
    bool _S3528 = (F32_isfinite((_S3522.y)));
    uint _S3529 = __ballot_sync(_S3205, _S3528);
    if(_S3528)
    {
        float _S3530 = _S3522.y;
        uint _S3531 = __ballot_sync(_S3205, true);
        v_3 = _S3530;
        _S3205 = _S3531;
    }
    else
    {
        uint _S3532 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3532;
    }
    *&((&_S3522)->y) = v_3;
    bool _S3533 = (F32_isfinite((_S3522.z)));
    uint _S3534 = __ballot_sync(_S3205, _S3533);
    if(_S3533)
    {
        float _S3535 = _S3522.z;
        uint _S3536 = __ballot_sync(_S3205, true);
        v_3 = _S3535;
        _S3205 = _S3536;
    }
    else
    {
        uint _S3537 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3537;
    }
    *&((&_S3522)->z) = v_3;
    float _S3538 = WaveActiveSum_0(_S3522.x, _S3205);
    *&((&_S3522)->x) = _S3538;
    float _S3539 = WaveActiveSum_0(_S3522.y, _S3205);
    *&((&_S3522)->y) = _S3539;
    float _S3540 = WaveActiveSum_0(_S3522.z, _S3205);
    *&((&_S3522)->z) = _S3540;
    uint warp_id_32 = thread_id_4 / 32U;
    uint _S3541 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_32] = _S3522.x;
        (*&_sh_block_reduce_shared_0)[warp_id_32 + 16U] = _S3522.y;
        (*&_sh_block_reduce_shared_0)[warp_id_32 + 32U] = _S3522.z;
        uint _S3542 = __ballot_sync(_S3205, true);
        _S3205 = _S3542;
    }
    else
    {
        uint _S3543 = __ballot_sync(_S3205, true);
        _S3205 = _S3543;
    }
    __syncthreads();
    bool _S3544 = warp_id_32 < 3U;
    uint _S3545 = __ballot_sync(_S3205, _S3544);
    if(_S3544)
    {
        bool _S3546 = lane_id_3 < 16U;
        uint _S3547 = __ballot_sync(_S3545, _S3546);
        if(_S3546)
        {
            float _S3548 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_32 * 16U];
            uint _S3549 = __ballot_sync(_S3545, true);
            v_3 = _S3548;
            _S3228 = _S3549;
        }
        else
        {
            uint _S3550 = __ballot_sync(_S3545, true);
            v_3 = 0.0f;
            _S3228 = _S3550;
        }
        float _S3551 = WaveActiveSum_0(v_3, _S3228);
        uint _S3552 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3553 = _S3551 != 0.0f;
            uint _S3554 = __ballot_sync(_S3228, true);
            _S3229 = _S3553;
            _S3228 = _S3554;
        }
        else
        {
            uint _S3555 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3555;
        }
        uint _S3556 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3557 = warp_id_32 == 0U;
            uint _S3558 = __ballot_sync(_S3556, _S3557);
            if(_S3557)
            {
                float _S3559 = atomicAdd(&((v_coeffs_34 + 6U)->x), _S3551);
                uint _S3560 = __ballot_sync(_S3556, true);
            }
            else
            {
                uint _S3561 = _S3556 & (~_S3558);
                bool _S3562 = warp_id_32 == 1U;
                uint _S3563 = __ballot_sync(_S3561, _S3562);
                if(_S3562)
                {
                    float _S3564 = atomicAdd(&((v_coeffs_34 + 6U)->y), _S3551);
                    uint _S3565 = __ballot_sync(_S3561, true);
                }
                else
                {
                    float _S3566 = atomicAdd(&((v_coeffs_34 + 6U)->z), _S3551);
                    uint _S3567 = __ballot_sync(_S3561, true);
                }
                uint _S3568 = __ballot_sync(_S3556, true);
            }
            uint _S3569 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3570 = __ballot_sync(_S3228, true);
        }
        uint _S3571 = __ballot_sync(_S3205, true);
        _S3205 = _S3571;
    }
    else
    {
        uint _S3572 = __ballot_sync(_S3205, true);
        _S3205 = _S3572;
    }
    float3  _S3573 = make_float3 (pSH8_23) * v_colors_34;
    float3  _S3574 = _S3573;
    bool _S3575 = (F32_isfinite((_S3573.x)));
    uint _S3576 = __ballot_sync(_S3205, _S3575);
    if(_S3575)
    {
        float _S3577 = _S3574.x;
        uint _S3578 = __ballot_sync(_S3205, true);
        v_3 = _S3577;
        _S3205 = _S3578;
    }
    else
    {
        uint _S3579 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3579;
    }
    *&((&_S3574)->x) = v_3;
    bool _S3580 = (F32_isfinite((_S3574.y)));
    uint _S3581 = __ballot_sync(_S3205, _S3580);
    if(_S3580)
    {
        float _S3582 = _S3574.y;
        uint _S3583 = __ballot_sync(_S3205, true);
        v_3 = _S3582;
        _S3205 = _S3583;
    }
    else
    {
        uint _S3584 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3584;
    }
    *&((&_S3574)->y) = v_3;
    bool _S3585 = (F32_isfinite((_S3574.z)));
    uint _S3586 = __ballot_sync(_S3205, _S3585);
    if(_S3585)
    {
        float _S3587 = _S3574.z;
        uint _S3588 = __ballot_sync(_S3205, true);
        v_3 = _S3587;
        _S3205 = _S3588;
    }
    else
    {
        uint _S3589 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3589;
    }
    *&((&_S3574)->z) = v_3;
    float _S3590 = WaveActiveSum_0(_S3574.x, _S3205);
    *&((&_S3574)->x) = _S3590;
    float _S3591 = WaveActiveSum_0(_S3574.y, _S3205);
    *&((&_S3574)->y) = _S3591;
    float _S3592 = WaveActiveSum_0(_S3574.z, _S3205);
    *&((&_S3574)->z) = _S3592;
    uint warp_id_33 = thread_id_4 / 32U;
    uint _S3593 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_33] = _S3574.x;
        (*&_sh_block_reduce_shared_0)[warp_id_33 + 16U] = _S3574.y;
        (*&_sh_block_reduce_shared_0)[warp_id_33 + 32U] = _S3574.z;
        uint _S3594 = __ballot_sync(_S3205, true);
        _S3205 = _S3594;
    }
    else
    {
        uint _S3595 = __ballot_sync(_S3205, true);
        _S3205 = _S3595;
    }
    __syncthreads();
    bool _S3596 = warp_id_33 < 3U;
    uint _S3597 = __ballot_sync(_S3205, _S3596);
    if(_S3596)
    {
        bool _S3598 = lane_id_3 < 16U;
        uint _S3599 = __ballot_sync(_S3597, _S3598);
        if(_S3598)
        {
            float _S3600 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_33 * 16U];
            uint _S3601 = __ballot_sync(_S3597, true);
            v_3 = _S3600;
            _S3228 = _S3601;
        }
        else
        {
            uint _S3602 = __ballot_sync(_S3597, true);
            v_3 = 0.0f;
            _S3228 = _S3602;
        }
        float _S3603 = WaveActiveSum_0(v_3, _S3228);
        uint _S3604 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3605 = _S3603 != 0.0f;
            uint _S3606 = __ballot_sync(_S3228, true);
            _S3229 = _S3605;
            _S3228 = _S3606;
        }
        else
        {
            uint _S3607 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3607;
        }
        uint _S3608 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3609 = warp_id_33 == 0U;
            uint _S3610 = __ballot_sync(_S3608, _S3609);
            if(_S3609)
            {
                float _S3611 = atomicAdd(&((v_coeffs_34 + 7U)->x), _S3603);
                uint _S3612 = __ballot_sync(_S3608, true);
            }
            else
            {
                uint _S3613 = _S3608 & (~_S3610);
                bool _S3614 = warp_id_33 == 1U;
                uint _S3615 = __ballot_sync(_S3613, _S3614);
                if(_S3614)
                {
                    float _S3616 = atomicAdd(&((v_coeffs_34 + 7U)->y), _S3603);
                    uint _S3617 = __ballot_sync(_S3613, true);
                }
                else
                {
                    float _S3618 = atomicAdd(&((v_coeffs_34 + 7U)->z), _S3603);
                    uint _S3619 = __ballot_sync(_S3613, true);
                }
                uint _S3620 = __ballot_sync(_S3608, true);
            }
            uint _S3621 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3622 = __ballot_sync(_S3228, true);
        }
        uint _S3623 = __ballot_sync(_S3205, true);
        _S3205 = _S3623;
    }
    else
    {
        uint _S3624 = __ballot_sync(_S3205, true);
        _S3205 = _S3624;
    }
    float fC1_y_13 = -2.0f * y_38;
    float fS1_x_13 = 2.0f * y_38;
    float pSH6_z_10 = 1.89234936237335205f * z_37;
    float pSH8_x_20 = 0.54627424478530884f * _S3364;
    float3  * _S3625 = coeffs_49 + int(3);
    float3  * _S3626 = coeffs_49 + int(7);
    float3  * _S3627 = coeffs_49 + int(6);
    float v_x_31 = _S3361 + dot_0(v_colors_34, make_float3 (0.54627424478530884f * fS1_x_13) * *_S3625 + make_float3 (pSH8_x_20) * *_S3626 + make_float3 (fTmp0B_29) * *_S3627);
    float3  * _S3628 = coeffs_49 + int(4);
    float v_y_31 = _S3362 + dot_0(v_colors_34, make_float3 (pSH8_x_20) * *_S3625 + make_float3 (0.54627424478530884f * fC1_y_13) * *_S3626 + make_float3 (fTmp0B_29) * *_S3628);
    float v_z_31 = _S3363 + dot_0(v_colors_34, make_float3 (pSH6_z_10) * *(coeffs_49 + int(5)) + make_float3 (-1.09254848957061768f * x_40) * *_S3627 + make_float3 (-1.09254848957061768f * y_38) * *_S3628);
    float fTmp0C_19 = -2.28522896766662598f * z2_19 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_37;
    float fC2_9 = x_40 * fC1_19 - y_38 * fS1_19;
    float fS2_9 = x_40 * fS1_19 + y_38 * fC1_19;
    float pSH12_17 = z_37 * (1.86588168144226074f * z2_19 - 1.11952900886535645f);
    float pSH13_15 = fTmp0C_19 * x_40;
    float pSH11_15 = fTmp0C_19 * y_38;
    float pSH14_15 = fTmp1B_19 * fC1_19;
    float pSH10_15 = fTmp1B_19 * fS1_19;
    float pSH15_15 = -0.59004360437393188f * fC2_9;
    float3  _S3629 = make_float3 (-0.59004360437393188f * fS2_9) * v_colors_34;
    float3  _S3630 = _S3629;
    bool _S3631 = (F32_isfinite((_S3629.x)));
    uint _S3632 = __ballot_sync(_S3205, _S3631);
    if(_S3631)
    {
        float _S3633 = _S3630.x;
        uint _S3634 = __ballot_sync(_S3205, true);
        v_3 = _S3633;
        _S3205 = _S3634;
    }
    else
    {
        uint _S3635 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3635;
    }
    *&((&_S3630)->x) = v_3;
    bool _S3636 = (F32_isfinite((_S3630.y)));
    uint _S3637 = __ballot_sync(_S3205, _S3636);
    if(_S3636)
    {
        float _S3638 = _S3630.y;
        uint _S3639 = __ballot_sync(_S3205, true);
        v_3 = _S3638;
        _S3205 = _S3639;
    }
    else
    {
        uint _S3640 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3640;
    }
    *&((&_S3630)->y) = v_3;
    bool _S3641 = (F32_isfinite((_S3630.z)));
    uint _S3642 = __ballot_sync(_S3205, _S3641);
    if(_S3641)
    {
        float _S3643 = _S3630.z;
        uint _S3644 = __ballot_sync(_S3205, true);
        v_3 = _S3643;
        _S3205 = _S3644;
    }
    else
    {
        uint _S3645 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3645;
    }
    *&((&_S3630)->z) = v_3;
    float _S3646 = WaveActiveSum_0(_S3630.x, _S3205);
    *&((&_S3630)->x) = _S3646;
    float _S3647 = WaveActiveSum_0(_S3630.y, _S3205);
    *&((&_S3630)->y) = _S3647;
    float _S3648 = WaveActiveSum_0(_S3630.z, _S3205);
    *&((&_S3630)->z) = _S3648;
    uint warp_id_34 = thread_id_4 / 32U;
    uint _S3649 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_34] = _S3630.x;
        (*&_sh_block_reduce_shared_0)[warp_id_34 + 16U] = _S3630.y;
        (*&_sh_block_reduce_shared_0)[warp_id_34 + 32U] = _S3630.z;
        uint _S3650 = __ballot_sync(_S3205, true);
        _S3205 = _S3650;
    }
    else
    {
        uint _S3651 = __ballot_sync(_S3205, true);
        _S3205 = _S3651;
    }
    __syncthreads();
    bool _S3652 = warp_id_34 < 3U;
    uint _S3653 = __ballot_sync(_S3205, _S3652);
    if(_S3652)
    {
        bool _S3654 = lane_id_3 < 16U;
        uint _S3655 = __ballot_sync(_S3653, _S3654);
        if(_S3654)
        {
            float _S3656 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_34 * 16U];
            uint _S3657 = __ballot_sync(_S3653, true);
            v_3 = _S3656;
            _S3228 = _S3657;
        }
        else
        {
            uint _S3658 = __ballot_sync(_S3653, true);
            v_3 = 0.0f;
            _S3228 = _S3658;
        }
        float _S3659 = WaveActiveSum_0(v_3, _S3228);
        uint _S3660 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3661 = _S3659 != 0.0f;
            uint _S3662 = __ballot_sync(_S3228, true);
            _S3229 = _S3661;
            _S3228 = _S3662;
        }
        else
        {
            uint _S3663 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3663;
        }
        uint _S3664 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3665 = warp_id_34 == 0U;
            uint _S3666 = __ballot_sync(_S3664, _S3665);
            if(_S3665)
            {
                float _S3667 = atomicAdd(&((v_coeffs_34 + 8U)->x), _S3659);
                uint _S3668 = __ballot_sync(_S3664, true);
            }
            else
            {
                uint _S3669 = _S3664 & (~_S3666);
                bool _S3670 = warp_id_34 == 1U;
                uint _S3671 = __ballot_sync(_S3669, _S3670);
                if(_S3670)
                {
                    float _S3672 = atomicAdd(&((v_coeffs_34 + 8U)->y), _S3659);
                    uint _S3673 = __ballot_sync(_S3669, true);
                }
                else
                {
                    float _S3674 = atomicAdd(&((v_coeffs_34 + 8U)->z), _S3659);
                    uint _S3675 = __ballot_sync(_S3669, true);
                }
                uint _S3676 = __ballot_sync(_S3664, true);
            }
            uint _S3677 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3678 = __ballot_sync(_S3228, true);
        }
        uint _S3679 = __ballot_sync(_S3205, true);
        _S3205 = _S3679;
    }
    else
    {
        uint _S3680 = __ballot_sync(_S3205, true);
        _S3205 = _S3680;
    }
    float3  _S3681 = make_float3 (pSH10_15) * v_colors_34;
    float3  _S3682 = _S3681;
    bool _S3683 = (F32_isfinite((_S3681.x)));
    uint _S3684 = __ballot_sync(_S3205, _S3683);
    if(_S3683)
    {
        float _S3685 = _S3682.x;
        uint _S3686 = __ballot_sync(_S3205, true);
        v_3 = _S3685;
        _S3205 = _S3686;
    }
    else
    {
        uint _S3687 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3687;
    }
    *&((&_S3682)->x) = v_3;
    bool _S3688 = (F32_isfinite((_S3682.y)));
    uint _S3689 = __ballot_sync(_S3205, _S3688);
    if(_S3688)
    {
        float _S3690 = _S3682.y;
        uint _S3691 = __ballot_sync(_S3205, true);
        v_3 = _S3690;
        _S3205 = _S3691;
    }
    else
    {
        uint _S3692 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3692;
    }
    *&((&_S3682)->y) = v_3;
    bool _S3693 = (F32_isfinite((_S3682.z)));
    uint _S3694 = __ballot_sync(_S3205, _S3693);
    if(_S3693)
    {
        float _S3695 = _S3682.z;
        uint _S3696 = __ballot_sync(_S3205, true);
        v_3 = _S3695;
        _S3205 = _S3696;
    }
    else
    {
        uint _S3697 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3697;
    }
    *&((&_S3682)->z) = v_3;
    float _S3698 = WaveActiveSum_0(_S3682.x, _S3205);
    *&((&_S3682)->x) = _S3698;
    float _S3699 = WaveActiveSum_0(_S3682.y, _S3205);
    *&((&_S3682)->y) = _S3699;
    float _S3700 = WaveActiveSum_0(_S3682.z, _S3205);
    *&((&_S3682)->z) = _S3700;
    uint warp_id_35 = thread_id_4 / 32U;
    uint _S3701 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_35] = _S3682.x;
        (*&_sh_block_reduce_shared_0)[warp_id_35 + 16U] = _S3682.y;
        (*&_sh_block_reduce_shared_0)[warp_id_35 + 32U] = _S3682.z;
        uint _S3702 = __ballot_sync(_S3205, true);
        _S3205 = _S3702;
    }
    else
    {
        uint _S3703 = __ballot_sync(_S3205, true);
        _S3205 = _S3703;
    }
    __syncthreads();
    bool _S3704 = warp_id_35 < 3U;
    uint _S3705 = __ballot_sync(_S3205, _S3704);
    if(_S3704)
    {
        bool _S3706 = lane_id_3 < 16U;
        uint _S3707 = __ballot_sync(_S3705, _S3706);
        if(_S3706)
        {
            float _S3708 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_35 * 16U];
            uint _S3709 = __ballot_sync(_S3705, true);
            v_3 = _S3708;
            _S3228 = _S3709;
        }
        else
        {
            uint _S3710 = __ballot_sync(_S3705, true);
            v_3 = 0.0f;
            _S3228 = _S3710;
        }
        float _S3711 = WaveActiveSum_0(v_3, _S3228);
        uint _S3712 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3713 = _S3711 != 0.0f;
            uint _S3714 = __ballot_sync(_S3228, true);
            _S3229 = _S3713;
            _S3228 = _S3714;
        }
        else
        {
            uint _S3715 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3715;
        }
        uint _S3716 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3717 = warp_id_35 == 0U;
            uint _S3718 = __ballot_sync(_S3716, _S3717);
            if(_S3717)
            {
                float _S3719 = atomicAdd(&((v_coeffs_34 + 9U)->x), _S3711);
                uint _S3720 = __ballot_sync(_S3716, true);
            }
            else
            {
                uint _S3721 = _S3716 & (~_S3718);
                bool _S3722 = warp_id_35 == 1U;
                uint _S3723 = __ballot_sync(_S3721, _S3722);
                if(_S3722)
                {
                    float _S3724 = atomicAdd(&((v_coeffs_34 + 9U)->y), _S3711);
                    uint _S3725 = __ballot_sync(_S3721, true);
                }
                else
                {
                    float _S3726 = atomicAdd(&((v_coeffs_34 + 9U)->z), _S3711);
                    uint _S3727 = __ballot_sync(_S3721, true);
                }
                uint _S3728 = __ballot_sync(_S3716, true);
            }
            uint _S3729 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3730 = __ballot_sync(_S3228, true);
        }
        uint _S3731 = __ballot_sync(_S3205, true);
        _S3205 = _S3731;
    }
    else
    {
        uint _S3732 = __ballot_sync(_S3205, true);
        _S3205 = _S3732;
    }
    float3  _S3733 = make_float3 (pSH11_15) * v_colors_34;
    float3  _S3734 = _S3733;
    bool _S3735 = (F32_isfinite((_S3733.x)));
    uint _S3736 = __ballot_sync(_S3205, _S3735);
    if(_S3735)
    {
        float _S3737 = _S3734.x;
        uint _S3738 = __ballot_sync(_S3205, true);
        v_3 = _S3737;
        _S3205 = _S3738;
    }
    else
    {
        uint _S3739 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3739;
    }
    *&((&_S3734)->x) = v_3;
    bool _S3740 = (F32_isfinite((_S3734.y)));
    uint _S3741 = __ballot_sync(_S3205, _S3740);
    if(_S3740)
    {
        float _S3742 = _S3734.y;
        uint _S3743 = __ballot_sync(_S3205, true);
        v_3 = _S3742;
        _S3205 = _S3743;
    }
    else
    {
        uint _S3744 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3744;
    }
    *&((&_S3734)->y) = v_3;
    bool _S3745 = (F32_isfinite((_S3734.z)));
    uint _S3746 = __ballot_sync(_S3205, _S3745);
    if(_S3745)
    {
        float _S3747 = _S3734.z;
        uint _S3748 = __ballot_sync(_S3205, true);
        v_3 = _S3747;
        _S3205 = _S3748;
    }
    else
    {
        uint _S3749 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3749;
    }
    *&((&_S3734)->z) = v_3;
    float _S3750 = WaveActiveSum_0(_S3734.x, _S3205);
    *&((&_S3734)->x) = _S3750;
    float _S3751 = WaveActiveSum_0(_S3734.y, _S3205);
    *&((&_S3734)->y) = _S3751;
    float _S3752 = WaveActiveSum_0(_S3734.z, _S3205);
    *&((&_S3734)->z) = _S3752;
    uint warp_id_36 = thread_id_4 / 32U;
    uint _S3753 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_36] = _S3734.x;
        (*&_sh_block_reduce_shared_0)[warp_id_36 + 16U] = _S3734.y;
        (*&_sh_block_reduce_shared_0)[warp_id_36 + 32U] = _S3734.z;
        uint _S3754 = __ballot_sync(_S3205, true);
        _S3205 = _S3754;
    }
    else
    {
        uint _S3755 = __ballot_sync(_S3205, true);
        _S3205 = _S3755;
    }
    __syncthreads();
    bool _S3756 = warp_id_36 < 3U;
    uint _S3757 = __ballot_sync(_S3205, _S3756);
    if(_S3756)
    {
        bool _S3758 = lane_id_3 < 16U;
        uint _S3759 = __ballot_sync(_S3757, _S3758);
        if(_S3758)
        {
            float _S3760 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_36 * 16U];
            uint _S3761 = __ballot_sync(_S3757, true);
            v_3 = _S3760;
            _S3228 = _S3761;
        }
        else
        {
            uint _S3762 = __ballot_sync(_S3757, true);
            v_3 = 0.0f;
            _S3228 = _S3762;
        }
        float _S3763 = WaveActiveSum_0(v_3, _S3228);
        uint _S3764 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3765 = _S3763 != 0.0f;
            uint _S3766 = __ballot_sync(_S3228, true);
            _S3229 = _S3765;
            _S3228 = _S3766;
        }
        else
        {
            uint _S3767 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3767;
        }
        uint _S3768 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3769 = warp_id_36 == 0U;
            uint _S3770 = __ballot_sync(_S3768, _S3769);
            if(_S3769)
            {
                float _S3771 = atomicAdd(&((v_coeffs_34 + 10U)->x), _S3763);
                uint _S3772 = __ballot_sync(_S3768, true);
            }
            else
            {
                uint _S3773 = _S3768 & (~_S3770);
                bool _S3774 = warp_id_36 == 1U;
                uint _S3775 = __ballot_sync(_S3773, _S3774);
                if(_S3774)
                {
                    float _S3776 = atomicAdd(&((v_coeffs_34 + 10U)->y), _S3763);
                    uint _S3777 = __ballot_sync(_S3773, true);
                }
                else
                {
                    float _S3778 = atomicAdd(&((v_coeffs_34 + 10U)->z), _S3763);
                    uint _S3779 = __ballot_sync(_S3773, true);
                }
                uint _S3780 = __ballot_sync(_S3768, true);
            }
            uint _S3781 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3782 = __ballot_sync(_S3228, true);
        }
        uint _S3783 = __ballot_sync(_S3205, true);
        _S3205 = _S3783;
    }
    else
    {
        uint _S3784 = __ballot_sync(_S3205, true);
        _S3205 = _S3784;
    }
    float3  _S3785 = make_float3 (pSH12_17) * v_colors_34;
    float3  _S3786 = _S3785;
    bool _S3787 = (F32_isfinite((_S3785.x)));
    uint _S3788 = __ballot_sync(_S3205, _S3787);
    if(_S3787)
    {
        float _S3789 = _S3786.x;
        uint _S3790 = __ballot_sync(_S3205, true);
        v_3 = _S3789;
        _S3205 = _S3790;
    }
    else
    {
        uint _S3791 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3791;
    }
    *&((&_S3786)->x) = v_3;
    bool _S3792 = (F32_isfinite((_S3786.y)));
    uint _S3793 = __ballot_sync(_S3205, _S3792);
    if(_S3792)
    {
        float _S3794 = _S3786.y;
        uint _S3795 = __ballot_sync(_S3205, true);
        v_3 = _S3794;
        _S3205 = _S3795;
    }
    else
    {
        uint _S3796 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3796;
    }
    *&((&_S3786)->y) = v_3;
    bool _S3797 = (F32_isfinite((_S3786.z)));
    uint _S3798 = __ballot_sync(_S3205, _S3797);
    if(_S3797)
    {
        float _S3799 = _S3786.z;
        uint _S3800 = __ballot_sync(_S3205, true);
        v_3 = _S3799;
        _S3205 = _S3800;
    }
    else
    {
        uint _S3801 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3801;
    }
    *&((&_S3786)->z) = v_3;
    float _S3802 = WaveActiveSum_0(_S3786.x, _S3205);
    *&((&_S3786)->x) = _S3802;
    float _S3803 = WaveActiveSum_0(_S3786.y, _S3205);
    *&((&_S3786)->y) = _S3803;
    float _S3804 = WaveActiveSum_0(_S3786.z, _S3205);
    *&((&_S3786)->z) = _S3804;
    uint warp_id_37 = thread_id_4 / 32U;
    uint _S3805 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_37] = _S3786.x;
        (*&_sh_block_reduce_shared_0)[warp_id_37 + 16U] = _S3786.y;
        (*&_sh_block_reduce_shared_0)[warp_id_37 + 32U] = _S3786.z;
        uint _S3806 = __ballot_sync(_S3205, true);
        _S3205 = _S3806;
    }
    else
    {
        uint _S3807 = __ballot_sync(_S3205, true);
        _S3205 = _S3807;
    }
    __syncthreads();
    bool _S3808 = warp_id_37 < 3U;
    uint _S3809 = __ballot_sync(_S3205, _S3808);
    if(_S3808)
    {
        bool _S3810 = lane_id_3 < 16U;
        uint _S3811 = __ballot_sync(_S3809, _S3810);
        if(_S3810)
        {
            float _S3812 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_37 * 16U];
            uint _S3813 = __ballot_sync(_S3809, true);
            v_3 = _S3812;
            _S3228 = _S3813;
        }
        else
        {
            uint _S3814 = __ballot_sync(_S3809, true);
            v_3 = 0.0f;
            _S3228 = _S3814;
        }
        float _S3815 = WaveActiveSum_0(v_3, _S3228);
        uint _S3816 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3817 = _S3815 != 0.0f;
            uint _S3818 = __ballot_sync(_S3228, true);
            _S3229 = _S3817;
            _S3228 = _S3818;
        }
        else
        {
            uint _S3819 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3819;
        }
        uint _S3820 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3821 = warp_id_37 == 0U;
            uint _S3822 = __ballot_sync(_S3820, _S3821);
            if(_S3821)
            {
                float _S3823 = atomicAdd(&((v_coeffs_34 + 11U)->x), _S3815);
                uint _S3824 = __ballot_sync(_S3820, true);
            }
            else
            {
                uint _S3825 = _S3820 & (~_S3822);
                bool _S3826 = warp_id_37 == 1U;
                uint _S3827 = __ballot_sync(_S3825, _S3826);
                if(_S3826)
                {
                    float _S3828 = atomicAdd(&((v_coeffs_34 + 11U)->y), _S3815);
                    uint _S3829 = __ballot_sync(_S3825, true);
                }
                else
                {
                    float _S3830 = atomicAdd(&((v_coeffs_34 + 11U)->z), _S3815);
                    uint _S3831 = __ballot_sync(_S3825, true);
                }
                uint _S3832 = __ballot_sync(_S3820, true);
            }
            uint _S3833 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3834 = __ballot_sync(_S3228, true);
        }
        uint _S3835 = __ballot_sync(_S3205, true);
        _S3205 = _S3835;
    }
    else
    {
        uint _S3836 = __ballot_sync(_S3205, true);
        _S3205 = _S3836;
    }
    float3  _S3837 = make_float3 (pSH13_15) * v_colors_34;
    float3  _S3838 = _S3837;
    bool _S3839 = (F32_isfinite((_S3837.x)));
    uint _S3840 = __ballot_sync(_S3205, _S3839);
    if(_S3839)
    {
        float _S3841 = _S3838.x;
        uint _S3842 = __ballot_sync(_S3205, true);
        v_3 = _S3841;
        _S3205 = _S3842;
    }
    else
    {
        uint _S3843 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3843;
    }
    *&((&_S3838)->x) = v_3;
    bool _S3844 = (F32_isfinite((_S3838.y)));
    uint _S3845 = __ballot_sync(_S3205, _S3844);
    if(_S3844)
    {
        float _S3846 = _S3838.y;
        uint _S3847 = __ballot_sync(_S3205, true);
        v_3 = _S3846;
        _S3205 = _S3847;
    }
    else
    {
        uint _S3848 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3848;
    }
    *&((&_S3838)->y) = v_3;
    bool _S3849 = (F32_isfinite((_S3838.z)));
    uint _S3850 = __ballot_sync(_S3205, _S3849);
    if(_S3849)
    {
        float _S3851 = _S3838.z;
        uint _S3852 = __ballot_sync(_S3205, true);
        v_3 = _S3851;
        _S3205 = _S3852;
    }
    else
    {
        uint _S3853 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3853;
    }
    *&((&_S3838)->z) = v_3;
    float _S3854 = WaveActiveSum_0(_S3838.x, _S3205);
    *&((&_S3838)->x) = _S3854;
    float _S3855 = WaveActiveSum_0(_S3838.y, _S3205);
    *&((&_S3838)->y) = _S3855;
    float _S3856 = WaveActiveSum_0(_S3838.z, _S3205);
    *&((&_S3838)->z) = _S3856;
    uint warp_id_38 = thread_id_4 / 32U;
    uint _S3857 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_38] = _S3838.x;
        (*&_sh_block_reduce_shared_0)[warp_id_38 + 16U] = _S3838.y;
        (*&_sh_block_reduce_shared_0)[warp_id_38 + 32U] = _S3838.z;
        uint _S3858 = __ballot_sync(_S3205, true);
        _S3205 = _S3858;
    }
    else
    {
        uint _S3859 = __ballot_sync(_S3205, true);
        _S3205 = _S3859;
    }
    __syncthreads();
    bool _S3860 = warp_id_38 < 3U;
    uint _S3861 = __ballot_sync(_S3205, _S3860);
    if(_S3860)
    {
        bool _S3862 = lane_id_3 < 16U;
        uint _S3863 = __ballot_sync(_S3861, _S3862);
        if(_S3862)
        {
            float _S3864 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_38 * 16U];
            uint _S3865 = __ballot_sync(_S3861, true);
            v_3 = _S3864;
            _S3228 = _S3865;
        }
        else
        {
            uint _S3866 = __ballot_sync(_S3861, true);
            v_3 = 0.0f;
            _S3228 = _S3866;
        }
        float _S3867 = WaveActiveSum_0(v_3, _S3228);
        uint _S3868 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3869 = _S3867 != 0.0f;
            uint _S3870 = __ballot_sync(_S3228, true);
            _S3229 = _S3869;
            _S3228 = _S3870;
        }
        else
        {
            uint _S3871 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3871;
        }
        uint _S3872 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3873 = warp_id_38 == 0U;
            uint _S3874 = __ballot_sync(_S3872, _S3873);
            if(_S3873)
            {
                float _S3875 = atomicAdd(&((v_coeffs_34 + 12U)->x), _S3867);
                uint _S3876 = __ballot_sync(_S3872, true);
            }
            else
            {
                uint _S3877 = _S3872 & (~_S3874);
                bool _S3878 = warp_id_38 == 1U;
                uint _S3879 = __ballot_sync(_S3877, _S3878);
                if(_S3878)
                {
                    float _S3880 = atomicAdd(&((v_coeffs_34 + 12U)->y), _S3867);
                    uint _S3881 = __ballot_sync(_S3877, true);
                }
                else
                {
                    float _S3882 = atomicAdd(&((v_coeffs_34 + 12U)->z), _S3867);
                    uint _S3883 = __ballot_sync(_S3877, true);
                }
                uint _S3884 = __ballot_sync(_S3872, true);
            }
            uint _S3885 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3886 = __ballot_sync(_S3228, true);
        }
        uint _S3887 = __ballot_sync(_S3205, true);
        _S3205 = _S3887;
    }
    else
    {
        uint _S3888 = __ballot_sync(_S3205, true);
        _S3205 = _S3888;
    }
    float3  _S3889 = make_float3 (pSH14_15) * v_colors_34;
    float3  _S3890 = _S3889;
    bool _S3891 = (F32_isfinite((_S3889.x)));
    uint _S3892 = __ballot_sync(_S3205, _S3891);
    if(_S3891)
    {
        float _S3893 = _S3890.x;
        uint _S3894 = __ballot_sync(_S3205, true);
        v_3 = _S3893;
        _S3205 = _S3894;
    }
    else
    {
        uint _S3895 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3895;
    }
    *&((&_S3890)->x) = v_3;
    bool _S3896 = (F32_isfinite((_S3890.y)));
    uint _S3897 = __ballot_sync(_S3205, _S3896);
    if(_S3896)
    {
        float _S3898 = _S3890.y;
        uint _S3899 = __ballot_sync(_S3205, true);
        v_3 = _S3898;
        _S3205 = _S3899;
    }
    else
    {
        uint _S3900 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3900;
    }
    *&((&_S3890)->y) = v_3;
    bool _S3901 = (F32_isfinite((_S3890.z)));
    uint _S3902 = __ballot_sync(_S3205, _S3901);
    if(_S3901)
    {
        float _S3903 = _S3890.z;
        uint _S3904 = __ballot_sync(_S3205, true);
        v_3 = _S3903;
        _S3205 = _S3904;
    }
    else
    {
        uint _S3905 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3905;
    }
    *&((&_S3890)->z) = v_3;
    float _S3906 = WaveActiveSum_0(_S3890.x, _S3205);
    *&((&_S3890)->x) = _S3906;
    float _S3907 = WaveActiveSum_0(_S3890.y, _S3205);
    *&((&_S3890)->y) = _S3907;
    float _S3908 = WaveActiveSum_0(_S3890.z, _S3205);
    *&((&_S3890)->z) = _S3908;
    uint warp_id_39 = thread_id_4 / 32U;
    uint _S3909 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_39] = _S3890.x;
        (*&_sh_block_reduce_shared_0)[warp_id_39 + 16U] = _S3890.y;
        (*&_sh_block_reduce_shared_0)[warp_id_39 + 32U] = _S3890.z;
        uint _S3910 = __ballot_sync(_S3205, true);
        _S3205 = _S3910;
    }
    else
    {
        uint _S3911 = __ballot_sync(_S3205, true);
        _S3205 = _S3911;
    }
    __syncthreads();
    bool _S3912 = warp_id_39 < 3U;
    uint _S3913 = __ballot_sync(_S3205, _S3912);
    if(_S3912)
    {
        bool _S3914 = lane_id_3 < 16U;
        uint _S3915 = __ballot_sync(_S3913, _S3914);
        if(_S3914)
        {
            float _S3916 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_39 * 16U];
            uint _S3917 = __ballot_sync(_S3913, true);
            v_3 = _S3916;
            _S3228 = _S3917;
        }
        else
        {
            uint _S3918 = __ballot_sync(_S3913, true);
            v_3 = 0.0f;
            _S3228 = _S3918;
        }
        float _S3919 = WaveActiveSum_0(v_3, _S3228);
        uint _S3920 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3921 = _S3919 != 0.0f;
            uint _S3922 = __ballot_sync(_S3228, true);
            _S3229 = _S3921;
            _S3228 = _S3922;
        }
        else
        {
            uint _S3923 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3923;
        }
        uint _S3924 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3925 = warp_id_39 == 0U;
            uint _S3926 = __ballot_sync(_S3924, _S3925);
            if(_S3925)
            {
                float _S3927 = atomicAdd(&((v_coeffs_34 + 13U)->x), _S3919);
                uint _S3928 = __ballot_sync(_S3924, true);
            }
            else
            {
                uint _S3929 = _S3924 & (~_S3926);
                bool _S3930 = warp_id_39 == 1U;
                uint _S3931 = __ballot_sync(_S3929, _S3930);
                if(_S3930)
                {
                    float _S3932 = atomicAdd(&((v_coeffs_34 + 13U)->y), _S3919);
                    uint _S3933 = __ballot_sync(_S3929, true);
                }
                else
                {
                    float _S3934 = atomicAdd(&((v_coeffs_34 + 13U)->z), _S3919);
                    uint _S3935 = __ballot_sync(_S3929, true);
                }
                uint _S3936 = __ballot_sync(_S3924, true);
            }
            uint _S3937 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3938 = __ballot_sync(_S3228, true);
        }
        uint _S3939 = __ballot_sync(_S3205, true);
        _S3205 = _S3939;
    }
    else
    {
        uint _S3940 = __ballot_sync(_S3205, true);
        _S3205 = _S3940;
    }
    float3  _S3941 = make_float3 (pSH15_15) * v_colors_34;
    float3  _S3942 = _S3941;
    bool _S3943 = (F32_isfinite((_S3941.x)));
    uint _S3944 = __ballot_sync(_S3205, _S3943);
    if(_S3943)
    {
        float _S3945 = _S3942.x;
        uint _S3946 = __ballot_sync(_S3205, true);
        v_3 = _S3945;
        _S3205 = _S3946;
    }
    else
    {
        uint _S3947 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3947;
    }
    *&((&_S3942)->x) = v_3;
    bool _S3948 = (F32_isfinite((_S3942.y)));
    uint _S3949 = __ballot_sync(_S3205, _S3948);
    if(_S3948)
    {
        float _S3950 = _S3942.y;
        uint _S3951 = __ballot_sync(_S3205, true);
        v_3 = _S3950;
        _S3205 = _S3951;
    }
    else
    {
        uint _S3952 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3952;
    }
    *&((&_S3942)->y) = v_3;
    bool _S3953 = (F32_isfinite((_S3942.z)));
    uint _S3954 = __ballot_sync(_S3205, _S3953);
    if(_S3953)
    {
        float _S3955 = _S3942.z;
        uint _S3956 = __ballot_sync(_S3205, true);
        v_3 = _S3955;
        _S3205 = _S3956;
    }
    else
    {
        uint _S3957 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S3957;
    }
    *&((&_S3942)->z) = v_3;
    float _S3958 = WaveActiveSum_0(_S3942.x, _S3205);
    *&((&_S3942)->x) = _S3958;
    float _S3959 = WaveActiveSum_0(_S3942.y, _S3205);
    *&((&_S3942)->y) = _S3959;
    float _S3960 = WaveActiveSum_0(_S3942.z, _S3205);
    *&((&_S3942)->z) = _S3960;
    uint warp_id_40 = thread_id_4 / 32U;
    uint _S3961 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_40] = _S3942.x;
        (*&_sh_block_reduce_shared_0)[warp_id_40 + 16U] = _S3942.y;
        (*&_sh_block_reduce_shared_0)[warp_id_40 + 32U] = _S3942.z;
        uint _S3962 = __ballot_sync(_S3205, true);
        _S3205 = _S3962;
    }
    else
    {
        uint _S3963 = __ballot_sync(_S3205, true);
        _S3205 = _S3963;
    }
    __syncthreads();
    bool _S3964 = warp_id_40 < 3U;
    uint _S3965 = __ballot_sync(_S3205, _S3964);
    if(_S3964)
    {
        bool _S3966 = lane_id_3 < 16U;
        uint _S3967 = __ballot_sync(_S3965, _S3966);
        if(_S3966)
        {
            float _S3968 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_40 * 16U];
            uint _S3969 = __ballot_sync(_S3965, true);
            v_3 = _S3968;
            _S3228 = _S3969;
        }
        else
        {
            uint _S3970 = __ballot_sync(_S3965, true);
            v_3 = 0.0f;
            _S3228 = _S3970;
        }
        float _S3971 = WaveActiveSum_0(v_3, _S3228);
        uint _S3972 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S3973 = _S3971 != 0.0f;
            uint _S3974 = __ballot_sync(_S3228, true);
            _S3229 = _S3973;
            _S3228 = _S3974;
        }
        else
        {
            uint _S3975 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S3975;
        }
        uint _S3976 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S3977 = warp_id_40 == 0U;
            uint _S3978 = __ballot_sync(_S3976, _S3977);
            if(_S3977)
            {
                float _S3979 = atomicAdd(&((v_coeffs_34 + 14U)->x), _S3971);
                uint _S3980 = __ballot_sync(_S3976, true);
            }
            else
            {
                uint _S3981 = _S3976 & (~_S3978);
                bool _S3982 = warp_id_40 == 1U;
                uint _S3983 = __ballot_sync(_S3981, _S3982);
                if(_S3982)
                {
                    float _S3984 = atomicAdd(&((v_coeffs_34 + 14U)->y), _S3971);
                    uint _S3985 = __ballot_sync(_S3981, true);
                }
                else
                {
                    float _S3986 = atomicAdd(&((v_coeffs_34 + 14U)->z), _S3971);
                    uint _S3987 = __ballot_sync(_S3981, true);
                }
                uint _S3988 = __ballot_sync(_S3976, true);
            }
            uint _S3989 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S3990 = __ballot_sync(_S3228, true);
        }
        uint _S3991 = __ballot_sync(_S3205, true);
        _S3205 = _S3991;
    }
    else
    {
        uint _S3992 = __ballot_sync(_S3205, true);
        _S3205 = _S3992;
    }
    float fTmp0C_z_13 = -4.57045793533325195f * z_37;
    float _S3993 = x_40 * _S3364;
    float fC2_x_6 = fC1_19 + _S3993 - y_38 * fS1_x_13;
    float _S3994 = y_38 * _S3364;
    float fC2_y_6 = x_40 * fC1_y_13 - fS1_19 - _S3994;
    float fS2_x_6 = fS1_19 + x_40 * fS1_x_13 + _S3994;
    float fS2_y_6 = _S3993 + fC1_19 + y_38 * fC1_y_13;
    float pSH12_z_8 = 5.59764480590820312f * z2_19 - 1.11952900886535645f;
    float pSH14_x_13 = fTmp1B_19 * _S3364;
    float3  * _S3995 = coeffs_49 + int(8);
    float3  * _S3996 = coeffs_49 + int(14);
    float3  * _S3997 = coeffs_49 + int(9);
    float3  * _S3998 = coeffs_49 + int(13);
    float3  * _S3999 = coeffs_49 + int(12);
    float v_x_32 = v_x_31 + dot_0(v_colors_34, make_float3 (-0.59004360437393188f * fS2_x_6) * *_S3995 + make_float3 (-0.59004360437393188f * fC2_x_6) * *_S3996 + make_float3 (fTmp1B_19 * fS1_x_13) * *_S3997 + make_float3 (pSH14_x_13) * *_S3998 + make_float3 (fTmp0C_19) * *_S3999);
    float3  * _S4000 = coeffs_49 + int(10);
    float v_y_32 = v_y_31 + dot_0(v_colors_34, make_float3 (-0.59004360437393188f * fS2_y_6) * *_S3995 + make_float3 (-0.59004360437393188f * fC2_y_6) * *_S3996 + make_float3 (pSH14_x_13) * *_S3997 + make_float3 (fTmp1B_19 * fC1_y_13) * *_S3998 + make_float3 (fTmp0C_19) * *_S4000);
    float v_z_32 = v_z_31 + dot_0(v_colors_34, make_float3 (pSH12_z_8) * *(coeffs_49 + int(11)) + make_float3 (fTmp0C_z_13 * x_40) * *_S3999 + make_float3 (fTmp0C_z_13 * y_38) * *_S4000 + make_float3 (1.44530570507049561f * fC1_19) * *_S3998 + make_float3 (1.44530570507049561f * fS1_19) * *_S3997);
    float fTmp0D_9 = z_37 * (-4.68332576751708984f * z2_19 + 2.00713968276977539f);
    float fTmp1C_9 = 3.31161141395568848f * z2_19 - 0.47308734059333801f;
    float fTmp2B_9 = -1.77013075351715088f * z_37;
    float pSH20_9 = 1.9843134880065918f * z_37 * pSH12_17 + -1.00623059272766113f * pSH6_25;
    float pSH21_7 = fTmp0D_9 * x_40;
    float pSH19_7 = fTmp0D_9 * y_38;
    float pSH22_7 = fTmp1C_9 * fC1_19;
    float pSH18_7 = fTmp1C_9 * fS1_19;
    float pSH23_7 = fTmp2B_9 * fC2_9;
    float pSH17_7 = fTmp2B_9 * fS2_9;
    float pSH24_7 = 0.62583571672439575f * (x_40 * fC2_9 - y_38 * fS2_9);
    float3  _S4001 = make_float3 (0.62583571672439575f * (x_40 * fS2_9 + y_38 * fC2_9)) * v_colors_34;
    float3  _S4002 = _S4001;
    bool _S4003 = (F32_isfinite((_S4001.x)));
    uint _S4004 = __ballot_sync(_S3205, _S4003);
    if(_S4003)
    {
        float _S4005 = _S4002.x;
        uint _S4006 = __ballot_sync(_S3205, true);
        v_3 = _S4005;
        _S3205 = _S4006;
    }
    else
    {
        uint _S4007 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4007;
    }
    *&((&_S4002)->x) = v_3;
    bool _S4008 = (F32_isfinite((_S4002.y)));
    uint _S4009 = __ballot_sync(_S3205, _S4008);
    if(_S4008)
    {
        float _S4010 = _S4002.y;
        uint _S4011 = __ballot_sync(_S3205, true);
        v_3 = _S4010;
        _S3205 = _S4011;
    }
    else
    {
        uint _S4012 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4012;
    }
    *&((&_S4002)->y) = v_3;
    bool _S4013 = (F32_isfinite((_S4002.z)));
    uint _S4014 = __ballot_sync(_S3205, _S4013);
    if(_S4013)
    {
        float _S4015 = _S4002.z;
        uint _S4016 = __ballot_sync(_S3205, true);
        v_3 = _S4015;
        _S3205 = _S4016;
    }
    else
    {
        uint _S4017 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4017;
    }
    *&((&_S4002)->z) = v_3;
    float _S4018 = WaveActiveSum_0(_S4002.x, _S3205);
    *&((&_S4002)->x) = _S4018;
    float _S4019 = WaveActiveSum_0(_S4002.y, _S3205);
    *&((&_S4002)->y) = _S4019;
    float _S4020 = WaveActiveSum_0(_S4002.z, _S3205);
    *&((&_S4002)->z) = _S4020;
    uint warp_id_41 = thread_id_4 / 32U;
    uint _S4021 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_41] = _S4002.x;
        (*&_sh_block_reduce_shared_0)[warp_id_41 + 16U] = _S4002.y;
        (*&_sh_block_reduce_shared_0)[warp_id_41 + 32U] = _S4002.z;
        uint _S4022 = __ballot_sync(_S3205, true);
        _S3205 = _S4022;
    }
    else
    {
        uint _S4023 = __ballot_sync(_S3205, true);
        _S3205 = _S4023;
    }
    __syncthreads();
    bool _S4024 = warp_id_41 < 3U;
    uint _S4025 = __ballot_sync(_S3205, _S4024);
    if(_S4024)
    {
        bool _S4026 = lane_id_3 < 16U;
        uint _S4027 = __ballot_sync(_S4025, _S4026);
        if(_S4026)
        {
            float _S4028 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_41 * 16U];
            uint _S4029 = __ballot_sync(_S4025, true);
            v_3 = _S4028;
            _S3228 = _S4029;
        }
        else
        {
            uint _S4030 = __ballot_sync(_S4025, true);
            v_3 = 0.0f;
            _S3228 = _S4030;
        }
        float _S4031 = WaveActiveSum_0(v_3, _S3228);
        uint _S4032 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4033 = _S4031 != 0.0f;
            uint _S4034 = __ballot_sync(_S3228, true);
            _S3229 = _S4033;
            _S3228 = _S4034;
        }
        else
        {
            uint _S4035 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4035;
        }
        uint _S4036 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4037 = warp_id_41 == 0U;
            uint _S4038 = __ballot_sync(_S4036, _S4037);
            if(_S4037)
            {
                float _S4039 = atomicAdd(&((v_coeffs_34 + 15U)->x), _S4031);
                uint _S4040 = __ballot_sync(_S4036, true);
            }
            else
            {
                uint _S4041 = _S4036 & (~_S4038);
                bool _S4042 = warp_id_41 == 1U;
                uint _S4043 = __ballot_sync(_S4041, _S4042);
                if(_S4042)
                {
                    float _S4044 = atomicAdd(&((v_coeffs_34 + 15U)->y), _S4031);
                    uint _S4045 = __ballot_sync(_S4041, true);
                }
                else
                {
                    float _S4046 = atomicAdd(&((v_coeffs_34 + 15U)->z), _S4031);
                    uint _S4047 = __ballot_sync(_S4041, true);
                }
                uint _S4048 = __ballot_sync(_S4036, true);
            }
            uint _S4049 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4050 = __ballot_sync(_S3228, true);
        }
        uint _S4051 = __ballot_sync(_S3205, true);
        _S3205 = _S4051;
    }
    else
    {
        uint _S4052 = __ballot_sync(_S3205, true);
        _S3205 = _S4052;
    }
    float3  _S4053 = make_float3 (pSH17_7) * v_colors_34;
    float3  _S4054 = _S4053;
    bool _S4055 = (F32_isfinite((_S4053.x)));
    uint _S4056 = __ballot_sync(_S3205, _S4055);
    if(_S4055)
    {
        float _S4057 = _S4054.x;
        uint _S4058 = __ballot_sync(_S3205, true);
        v_3 = _S4057;
        _S3205 = _S4058;
    }
    else
    {
        uint _S4059 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4059;
    }
    *&((&_S4054)->x) = v_3;
    bool _S4060 = (F32_isfinite((_S4054.y)));
    uint _S4061 = __ballot_sync(_S3205, _S4060);
    if(_S4060)
    {
        float _S4062 = _S4054.y;
        uint _S4063 = __ballot_sync(_S3205, true);
        v_3 = _S4062;
        _S3205 = _S4063;
    }
    else
    {
        uint _S4064 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4064;
    }
    *&((&_S4054)->y) = v_3;
    bool _S4065 = (F32_isfinite((_S4054.z)));
    uint _S4066 = __ballot_sync(_S3205, _S4065);
    if(_S4065)
    {
        float _S4067 = _S4054.z;
        uint _S4068 = __ballot_sync(_S3205, true);
        v_3 = _S4067;
        _S3205 = _S4068;
    }
    else
    {
        uint _S4069 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4069;
    }
    *&((&_S4054)->z) = v_3;
    float _S4070 = WaveActiveSum_0(_S4054.x, _S3205);
    *&((&_S4054)->x) = _S4070;
    float _S4071 = WaveActiveSum_0(_S4054.y, _S3205);
    *&((&_S4054)->y) = _S4071;
    float _S4072 = WaveActiveSum_0(_S4054.z, _S3205);
    *&((&_S4054)->z) = _S4072;
    uint warp_id_42 = thread_id_4 / 32U;
    uint _S4073 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_42] = _S4054.x;
        (*&_sh_block_reduce_shared_0)[warp_id_42 + 16U] = _S4054.y;
        (*&_sh_block_reduce_shared_0)[warp_id_42 + 32U] = _S4054.z;
        uint _S4074 = __ballot_sync(_S3205, true);
        _S3205 = _S4074;
    }
    else
    {
        uint _S4075 = __ballot_sync(_S3205, true);
        _S3205 = _S4075;
    }
    __syncthreads();
    bool _S4076 = warp_id_42 < 3U;
    uint _S4077 = __ballot_sync(_S3205, _S4076);
    if(_S4076)
    {
        bool _S4078 = lane_id_3 < 16U;
        uint _S4079 = __ballot_sync(_S4077, _S4078);
        if(_S4078)
        {
            float _S4080 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_42 * 16U];
            uint _S4081 = __ballot_sync(_S4077, true);
            v_3 = _S4080;
            _S3228 = _S4081;
        }
        else
        {
            uint _S4082 = __ballot_sync(_S4077, true);
            v_3 = 0.0f;
            _S3228 = _S4082;
        }
        float _S4083 = WaveActiveSum_0(v_3, _S3228);
        uint _S4084 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4085 = _S4083 != 0.0f;
            uint _S4086 = __ballot_sync(_S3228, true);
            _S3229 = _S4085;
            _S3228 = _S4086;
        }
        else
        {
            uint _S4087 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4087;
        }
        uint _S4088 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4089 = warp_id_42 == 0U;
            uint _S4090 = __ballot_sync(_S4088, _S4089);
            if(_S4089)
            {
                float _S4091 = atomicAdd(&((v_coeffs_34 + 16U)->x), _S4083);
                uint _S4092 = __ballot_sync(_S4088, true);
            }
            else
            {
                uint _S4093 = _S4088 & (~_S4090);
                bool _S4094 = warp_id_42 == 1U;
                uint _S4095 = __ballot_sync(_S4093, _S4094);
                if(_S4094)
                {
                    float _S4096 = atomicAdd(&((v_coeffs_34 + 16U)->y), _S4083);
                    uint _S4097 = __ballot_sync(_S4093, true);
                }
                else
                {
                    float _S4098 = atomicAdd(&((v_coeffs_34 + 16U)->z), _S4083);
                    uint _S4099 = __ballot_sync(_S4093, true);
                }
                uint _S4100 = __ballot_sync(_S4088, true);
            }
            uint _S4101 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4102 = __ballot_sync(_S3228, true);
        }
        uint _S4103 = __ballot_sync(_S3205, true);
        _S3205 = _S4103;
    }
    else
    {
        uint _S4104 = __ballot_sync(_S3205, true);
        _S3205 = _S4104;
    }
    float3  _S4105 = make_float3 (pSH18_7) * v_colors_34;
    float3  _S4106 = _S4105;
    bool _S4107 = (F32_isfinite((_S4105.x)));
    uint _S4108 = __ballot_sync(_S3205, _S4107);
    if(_S4107)
    {
        float _S4109 = _S4106.x;
        uint _S4110 = __ballot_sync(_S3205, true);
        v_3 = _S4109;
        _S3205 = _S4110;
    }
    else
    {
        uint _S4111 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4111;
    }
    *&((&_S4106)->x) = v_3;
    bool _S4112 = (F32_isfinite((_S4106.y)));
    uint _S4113 = __ballot_sync(_S3205, _S4112);
    if(_S4112)
    {
        float _S4114 = _S4106.y;
        uint _S4115 = __ballot_sync(_S3205, true);
        v_3 = _S4114;
        _S3205 = _S4115;
    }
    else
    {
        uint _S4116 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4116;
    }
    *&((&_S4106)->y) = v_3;
    bool _S4117 = (F32_isfinite((_S4106.z)));
    uint _S4118 = __ballot_sync(_S3205, _S4117);
    if(_S4117)
    {
        float _S4119 = _S4106.z;
        uint _S4120 = __ballot_sync(_S3205, true);
        v_3 = _S4119;
        _S3205 = _S4120;
    }
    else
    {
        uint _S4121 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4121;
    }
    *&((&_S4106)->z) = v_3;
    float _S4122 = WaveActiveSum_0(_S4106.x, _S3205);
    *&((&_S4106)->x) = _S4122;
    float _S4123 = WaveActiveSum_0(_S4106.y, _S3205);
    *&((&_S4106)->y) = _S4123;
    float _S4124 = WaveActiveSum_0(_S4106.z, _S3205);
    *&((&_S4106)->z) = _S4124;
    uint warp_id_43 = thread_id_4 / 32U;
    uint _S4125 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_43] = _S4106.x;
        (*&_sh_block_reduce_shared_0)[warp_id_43 + 16U] = _S4106.y;
        (*&_sh_block_reduce_shared_0)[warp_id_43 + 32U] = _S4106.z;
        uint _S4126 = __ballot_sync(_S3205, true);
        _S3205 = _S4126;
    }
    else
    {
        uint _S4127 = __ballot_sync(_S3205, true);
        _S3205 = _S4127;
    }
    __syncthreads();
    bool _S4128 = warp_id_43 < 3U;
    uint _S4129 = __ballot_sync(_S3205, _S4128);
    if(_S4128)
    {
        bool _S4130 = lane_id_3 < 16U;
        uint _S4131 = __ballot_sync(_S4129, _S4130);
        if(_S4130)
        {
            float _S4132 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_43 * 16U];
            uint _S4133 = __ballot_sync(_S4129, true);
            v_3 = _S4132;
            _S3228 = _S4133;
        }
        else
        {
            uint _S4134 = __ballot_sync(_S4129, true);
            v_3 = 0.0f;
            _S3228 = _S4134;
        }
        float _S4135 = WaveActiveSum_0(v_3, _S3228);
        uint _S4136 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4137 = _S4135 != 0.0f;
            uint _S4138 = __ballot_sync(_S3228, true);
            _S3229 = _S4137;
            _S3228 = _S4138;
        }
        else
        {
            uint _S4139 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4139;
        }
        uint _S4140 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4141 = warp_id_43 == 0U;
            uint _S4142 = __ballot_sync(_S4140, _S4141);
            if(_S4141)
            {
                float _S4143 = atomicAdd(&((v_coeffs_34 + 17U)->x), _S4135);
                uint _S4144 = __ballot_sync(_S4140, true);
            }
            else
            {
                uint _S4145 = _S4140 & (~_S4142);
                bool _S4146 = warp_id_43 == 1U;
                uint _S4147 = __ballot_sync(_S4145, _S4146);
                if(_S4146)
                {
                    float _S4148 = atomicAdd(&((v_coeffs_34 + 17U)->y), _S4135);
                    uint _S4149 = __ballot_sync(_S4145, true);
                }
                else
                {
                    float _S4150 = atomicAdd(&((v_coeffs_34 + 17U)->z), _S4135);
                    uint _S4151 = __ballot_sync(_S4145, true);
                }
                uint _S4152 = __ballot_sync(_S4140, true);
            }
            uint _S4153 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4154 = __ballot_sync(_S3228, true);
        }
        uint _S4155 = __ballot_sync(_S3205, true);
        _S3205 = _S4155;
    }
    else
    {
        uint _S4156 = __ballot_sync(_S3205, true);
        _S3205 = _S4156;
    }
    float3  _S4157 = make_float3 (pSH19_7) * v_colors_34;
    float3  _S4158 = _S4157;
    bool _S4159 = (F32_isfinite((_S4157.x)));
    uint _S4160 = __ballot_sync(_S3205, _S4159);
    if(_S4159)
    {
        float _S4161 = _S4158.x;
        uint _S4162 = __ballot_sync(_S3205, true);
        v_3 = _S4161;
        _S3205 = _S4162;
    }
    else
    {
        uint _S4163 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4163;
    }
    *&((&_S4158)->x) = v_3;
    bool _S4164 = (F32_isfinite((_S4158.y)));
    uint _S4165 = __ballot_sync(_S3205, _S4164);
    if(_S4164)
    {
        float _S4166 = _S4158.y;
        uint _S4167 = __ballot_sync(_S3205, true);
        v_3 = _S4166;
        _S3205 = _S4167;
    }
    else
    {
        uint _S4168 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4168;
    }
    *&((&_S4158)->y) = v_3;
    bool _S4169 = (F32_isfinite((_S4158.z)));
    uint _S4170 = __ballot_sync(_S3205, _S4169);
    if(_S4169)
    {
        float _S4171 = _S4158.z;
        uint _S4172 = __ballot_sync(_S3205, true);
        v_3 = _S4171;
        _S3205 = _S4172;
    }
    else
    {
        uint _S4173 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4173;
    }
    *&((&_S4158)->z) = v_3;
    float _S4174 = WaveActiveSum_0(_S4158.x, _S3205);
    *&((&_S4158)->x) = _S4174;
    float _S4175 = WaveActiveSum_0(_S4158.y, _S3205);
    *&((&_S4158)->y) = _S4175;
    float _S4176 = WaveActiveSum_0(_S4158.z, _S3205);
    *&((&_S4158)->z) = _S4176;
    uint warp_id_44 = thread_id_4 / 32U;
    uint _S4177 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_44] = _S4158.x;
        (*&_sh_block_reduce_shared_0)[warp_id_44 + 16U] = _S4158.y;
        (*&_sh_block_reduce_shared_0)[warp_id_44 + 32U] = _S4158.z;
        uint _S4178 = __ballot_sync(_S3205, true);
        _S3205 = _S4178;
    }
    else
    {
        uint _S4179 = __ballot_sync(_S3205, true);
        _S3205 = _S4179;
    }
    __syncthreads();
    bool _S4180 = warp_id_44 < 3U;
    uint _S4181 = __ballot_sync(_S3205, _S4180);
    if(_S4180)
    {
        bool _S4182 = lane_id_3 < 16U;
        uint _S4183 = __ballot_sync(_S4181, _S4182);
        if(_S4182)
        {
            float _S4184 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_44 * 16U];
            uint _S4185 = __ballot_sync(_S4181, true);
            v_3 = _S4184;
            _S3228 = _S4185;
        }
        else
        {
            uint _S4186 = __ballot_sync(_S4181, true);
            v_3 = 0.0f;
            _S3228 = _S4186;
        }
        float _S4187 = WaveActiveSum_0(v_3, _S3228);
        uint _S4188 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4189 = _S4187 != 0.0f;
            uint _S4190 = __ballot_sync(_S3228, true);
            _S3229 = _S4189;
            _S3228 = _S4190;
        }
        else
        {
            uint _S4191 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4191;
        }
        uint _S4192 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4193 = warp_id_44 == 0U;
            uint _S4194 = __ballot_sync(_S4192, _S4193);
            if(_S4193)
            {
                float _S4195 = atomicAdd(&((v_coeffs_34 + 18U)->x), _S4187);
                uint _S4196 = __ballot_sync(_S4192, true);
            }
            else
            {
                uint _S4197 = _S4192 & (~_S4194);
                bool _S4198 = warp_id_44 == 1U;
                uint _S4199 = __ballot_sync(_S4197, _S4198);
                if(_S4198)
                {
                    float _S4200 = atomicAdd(&((v_coeffs_34 + 18U)->y), _S4187);
                    uint _S4201 = __ballot_sync(_S4197, true);
                }
                else
                {
                    float _S4202 = atomicAdd(&((v_coeffs_34 + 18U)->z), _S4187);
                    uint _S4203 = __ballot_sync(_S4197, true);
                }
                uint _S4204 = __ballot_sync(_S4192, true);
            }
            uint _S4205 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4206 = __ballot_sync(_S3228, true);
        }
        uint _S4207 = __ballot_sync(_S3205, true);
        _S3205 = _S4207;
    }
    else
    {
        uint _S4208 = __ballot_sync(_S3205, true);
        _S3205 = _S4208;
    }
    float3  _S4209 = make_float3 (pSH20_9) * v_colors_34;
    float3  _S4210 = _S4209;
    bool _S4211 = (F32_isfinite((_S4209.x)));
    uint _S4212 = __ballot_sync(_S3205, _S4211);
    if(_S4211)
    {
        float _S4213 = _S4210.x;
        uint _S4214 = __ballot_sync(_S3205, true);
        v_3 = _S4213;
        _S3205 = _S4214;
    }
    else
    {
        uint _S4215 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4215;
    }
    *&((&_S4210)->x) = v_3;
    bool _S4216 = (F32_isfinite((_S4210.y)));
    uint _S4217 = __ballot_sync(_S3205, _S4216);
    if(_S4216)
    {
        float _S4218 = _S4210.y;
        uint _S4219 = __ballot_sync(_S3205, true);
        v_3 = _S4218;
        _S3205 = _S4219;
    }
    else
    {
        uint _S4220 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4220;
    }
    *&((&_S4210)->y) = v_3;
    bool _S4221 = (F32_isfinite((_S4210.z)));
    uint _S4222 = __ballot_sync(_S3205, _S4221);
    if(_S4221)
    {
        float _S4223 = _S4210.z;
        uint _S4224 = __ballot_sync(_S3205, true);
        v_3 = _S4223;
        _S3205 = _S4224;
    }
    else
    {
        uint _S4225 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4225;
    }
    *&((&_S4210)->z) = v_3;
    float _S4226 = WaveActiveSum_0(_S4210.x, _S3205);
    *&((&_S4210)->x) = _S4226;
    float _S4227 = WaveActiveSum_0(_S4210.y, _S3205);
    *&((&_S4210)->y) = _S4227;
    float _S4228 = WaveActiveSum_0(_S4210.z, _S3205);
    *&((&_S4210)->z) = _S4228;
    uint warp_id_45 = thread_id_4 / 32U;
    uint _S4229 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_45] = _S4210.x;
        (*&_sh_block_reduce_shared_0)[warp_id_45 + 16U] = _S4210.y;
        (*&_sh_block_reduce_shared_0)[warp_id_45 + 32U] = _S4210.z;
        uint _S4230 = __ballot_sync(_S3205, true);
        _S3205 = _S4230;
    }
    else
    {
        uint _S4231 = __ballot_sync(_S3205, true);
        _S3205 = _S4231;
    }
    __syncthreads();
    bool _S4232 = warp_id_45 < 3U;
    uint _S4233 = __ballot_sync(_S3205, _S4232);
    if(_S4232)
    {
        bool _S4234 = lane_id_3 < 16U;
        uint _S4235 = __ballot_sync(_S4233, _S4234);
        if(_S4234)
        {
            float _S4236 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_45 * 16U];
            uint _S4237 = __ballot_sync(_S4233, true);
            v_3 = _S4236;
            _S3228 = _S4237;
        }
        else
        {
            uint _S4238 = __ballot_sync(_S4233, true);
            v_3 = 0.0f;
            _S3228 = _S4238;
        }
        float _S4239 = WaveActiveSum_0(v_3, _S3228);
        uint _S4240 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4241 = _S4239 != 0.0f;
            uint _S4242 = __ballot_sync(_S3228, true);
            _S3229 = _S4241;
            _S3228 = _S4242;
        }
        else
        {
            uint _S4243 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4243;
        }
        uint _S4244 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4245 = warp_id_45 == 0U;
            uint _S4246 = __ballot_sync(_S4244, _S4245);
            if(_S4245)
            {
                float _S4247 = atomicAdd(&((v_coeffs_34 + 19U)->x), _S4239);
                uint _S4248 = __ballot_sync(_S4244, true);
            }
            else
            {
                uint _S4249 = _S4244 & (~_S4246);
                bool _S4250 = warp_id_45 == 1U;
                uint _S4251 = __ballot_sync(_S4249, _S4250);
                if(_S4250)
                {
                    float _S4252 = atomicAdd(&((v_coeffs_34 + 19U)->y), _S4239);
                    uint _S4253 = __ballot_sync(_S4249, true);
                }
                else
                {
                    float _S4254 = atomicAdd(&((v_coeffs_34 + 19U)->z), _S4239);
                    uint _S4255 = __ballot_sync(_S4249, true);
                }
                uint _S4256 = __ballot_sync(_S4244, true);
            }
            uint _S4257 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4258 = __ballot_sync(_S3228, true);
        }
        uint _S4259 = __ballot_sync(_S3205, true);
        _S3205 = _S4259;
    }
    else
    {
        uint _S4260 = __ballot_sync(_S3205, true);
        _S3205 = _S4260;
    }
    float3  _S4261 = make_float3 (pSH21_7) * v_colors_34;
    float3  _S4262 = _S4261;
    bool _S4263 = (F32_isfinite((_S4261.x)));
    uint _S4264 = __ballot_sync(_S3205, _S4263);
    if(_S4263)
    {
        float _S4265 = _S4262.x;
        uint _S4266 = __ballot_sync(_S3205, true);
        v_3 = _S4265;
        _S3205 = _S4266;
    }
    else
    {
        uint _S4267 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4267;
    }
    *&((&_S4262)->x) = v_3;
    bool _S4268 = (F32_isfinite((_S4262.y)));
    uint _S4269 = __ballot_sync(_S3205, _S4268);
    if(_S4268)
    {
        float _S4270 = _S4262.y;
        uint _S4271 = __ballot_sync(_S3205, true);
        v_3 = _S4270;
        _S3205 = _S4271;
    }
    else
    {
        uint _S4272 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4272;
    }
    *&((&_S4262)->y) = v_3;
    bool _S4273 = (F32_isfinite((_S4262.z)));
    uint _S4274 = __ballot_sync(_S3205, _S4273);
    if(_S4273)
    {
        float _S4275 = _S4262.z;
        uint _S4276 = __ballot_sync(_S3205, true);
        v_3 = _S4275;
        _S3205 = _S4276;
    }
    else
    {
        uint _S4277 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4277;
    }
    *&((&_S4262)->z) = v_3;
    float _S4278 = WaveActiveSum_0(_S4262.x, _S3205);
    *&((&_S4262)->x) = _S4278;
    float _S4279 = WaveActiveSum_0(_S4262.y, _S3205);
    *&((&_S4262)->y) = _S4279;
    float _S4280 = WaveActiveSum_0(_S4262.z, _S3205);
    *&((&_S4262)->z) = _S4280;
    uint warp_id_46 = thread_id_4 / 32U;
    uint _S4281 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_46] = _S4262.x;
        (*&_sh_block_reduce_shared_0)[warp_id_46 + 16U] = _S4262.y;
        (*&_sh_block_reduce_shared_0)[warp_id_46 + 32U] = _S4262.z;
        uint _S4282 = __ballot_sync(_S3205, true);
        _S3205 = _S4282;
    }
    else
    {
        uint _S4283 = __ballot_sync(_S3205, true);
        _S3205 = _S4283;
    }
    __syncthreads();
    bool _S4284 = warp_id_46 < 3U;
    uint _S4285 = __ballot_sync(_S3205, _S4284);
    if(_S4284)
    {
        bool _S4286 = lane_id_3 < 16U;
        uint _S4287 = __ballot_sync(_S4285, _S4286);
        if(_S4286)
        {
            float _S4288 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_46 * 16U];
            uint _S4289 = __ballot_sync(_S4285, true);
            v_3 = _S4288;
            _S3228 = _S4289;
        }
        else
        {
            uint _S4290 = __ballot_sync(_S4285, true);
            v_3 = 0.0f;
            _S3228 = _S4290;
        }
        float _S4291 = WaveActiveSum_0(v_3, _S3228);
        uint _S4292 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4293 = _S4291 != 0.0f;
            uint _S4294 = __ballot_sync(_S3228, true);
            _S3229 = _S4293;
            _S3228 = _S4294;
        }
        else
        {
            uint _S4295 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4295;
        }
        uint _S4296 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4297 = warp_id_46 == 0U;
            uint _S4298 = __ballot_sync(_S4296, _S4297);
            if(_S4297)
            {
                float _S4299 = atomicAdd(&((v_coeffs_34 + 20U)->x), _S4291);
                uint _S4300 = __ballot_sync(_S4296, true);
            }
            else
            {
                uint _S4301 = _S4296 & (~_S4298);
                bool _S4302 = warp_id_46 == 1U;
                uint _S4303 = __ballot_sync(_S4301, _S4302);
                if(_S4302)
                {
                    float _S4304 = atomicAdd(&((v_coeffs_34 + 20U)->y), _S4291);
                    uint _S4305 = __ballot_sync(_S4301, true);
                }
                else
                {
                    float _S4306 = atomicAdd(&((v_coeffs_34 + 20U)->z), _S4291);
                    uint _S4307 = __ballot_sync(_S4301, true);
                }
                uint _S4308 = __ballot_sync(_S4296, true);
            }
            uint _S4309 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4310 = __ballot_sync(_S3228, true);
        }
        uint _S4311 = __ballot_sync(_S3205, true);
        _S3205 = _S4311;
    }
    else
    {
        uint _S4312 = __ballot_sync(_S3205, true);
        _S3205 = _S4312;
    }
    float3  _S4313 = make_float3 (pSH22_7) * v_colors_34;
    float3  _S4314 = _S4313;
    bool _S4315 = (F32_isfinite((_S4313.x)));
    uint _S4316 = __ballot_sync(_S3205, _S4315);
    if(_S4315)
    {
        float _S4317 = _S4314.x;
        uint _S4318 = __ballot_sync(_S3205, true);
        v_3 = _S4317;
        _S3205 = _S4318;
    }
    else
    {
        uint _S4319 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4319;
    }
    *&((&_S4314)->x) = v_3;
    bool _S4320 = (F32_isfinite((_S4314.y)));
    uint _S4321 = __ballot_sync(_S3205, _S4320);
    if(_S4320)
    {
        float _S4322 = _S4314.y;
        uint _S4323 = __ballot_sync(_S3205, true);
        v_3 = _S4322;
        _S3205 = _S4323;
    }
    else
    {
        uint _S4324 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4324;
    }
    *&((&_S4314)->y) = v_3;
    bool _S4325 = (F32_isfinite((_S4314.z)));
    uint _S4326 = __ballot_sync(_S3205, _S4325);
    if(_S4325)
    {
        float _S4327 = _S4314.z;
        uint _S4328 = __ballot_sync(_S3205, true);
        v_3 = _S4327;
        _S3205 = _S4328;
    }
    else
    {
        uint _S4329 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4329;
    }
    *&((&_S4314)->z) = v_3;
    float _S4330 = WaveActiveSum_0(_S4314.x, _S3205);
    *&((&_S4314)->x) = _S4330;
    float _S4331 = WaveActiveSum_0(_S4314.y, _S3205);
    *&((&_S4314)->y) = _S4331;
    float _S4332 = WaveActiveSum_0(_S4314.z, _S3205);
    *&((&_S4314)->z) = _S4332;
    uint warp_id_47 = thread_id_4 / 32U;
    uint _S4333 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_47] = _S4314.x;
        (*&_sh_block_reduce_shared_0)[warp_id_47 + 16U] = _S4314.y;
        (*&_sh_block_reduce_shared_0)[warp_id_47 + 32U] = _S4314.z;
        uint _S4334 = __ballot_sync(_S3205, true);
        _S3205 = _S4334;
    }
    else
    {
        uint _S4335 = __ballot_sync(_S3205, true);
        _S3205 = _S4335;
    }
    __syncthreads();
    bool _S4336 = warp_id_47 < 3U;
    uint _S4337 = __ballot_sync(_S3205, _S4336);
    if(_S4336)
    {
        bool _S4338 = lane_id_3 < 16U;
        uint _S4339 = __ballot_sync(_S4337, _S4338);
        if(_S4338)
        {
            float _S4340 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_47 * 16U];
            uint _S4341 = __ballot_sync(_S4337, true);
            v_3 = _S4340;
            _S3228 = _S4341;
        }
        else
        {
            uint _S4342 = __ballot_sync(_S4337, true);
            v_3 = 0.0f;
            _S3228 = _S4342;
        }
        float _S4343 = WaveActiveSum_0(v_3, _S3228);
        uint _S4344 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4345 = _S4343 != 0.0f;
            uint _S4346 = __ballot_sync(_S3228, true);
            _S3229 = _S4345;
            _S3228 = _S4346;
        }
        else
        {
            uint _S4347 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4347;
        }
        uint _S4348 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4349 = warp_id_47 == 0U;
            uint _S4350 = __ballot_sync(_S4348, _S4349);
            if(_S4349)
            {
                float _S4351 = atomicAdd(&((v_coeffs_34 + 21U)->x), _S4343);
                uint _S4352 = __ballot_sync(_S4348, true);
            }
            else
            {
                uint _S4353 = _S4348 & (~_S4350);
                bool _S4354 = warp_id_47 == 1U;
                uint _S4355 = __ballot_sync(_S4353, _S4354);
                if(_S4354)
                {
                    float _S4356 = atomicAdd(&((v_coeffs_34 + 21U)->y), _S4343);
                    uint _S4357 = __ballot_sync(_S4353, true);
                }
                else
                {
                    float _S4358 = atomicAdd(&((v_coeffs_34 + 21U)->z), _S4343);
                    uint _S4359 = __ballot_sync(_S4353, true);
                }
                uint _S4360 = __ballot_sync(_S4348, true);
            }
            uint _S4361 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4362 = __ballot_sync(_S3228, true);
        }
        uint _S4363 = __ballot_sync(_S3205, true);
        _S3205 = _S4363;
    }
    else
    {
        uint _S4364 = __ballot_sync(_S3205, true);
        _S3205 = _S4364;
    }
    float3  _S4365 = make_float3 (pSH23_7) * v_colors_34;
    float3  _S4366 = _S4365;
    bool _S4367 = (F32_isfinite((_S4365.x)));
    uint _S4368 = __ballot_sync(_S3205, _S4367);
    if(_S4367)
    {
        float _S4369 = _S4366.x;
        uint _S4370 = __ballot_sync(_S3205, true);
        v_3 = _S4369;
        _S3205 = _S4370;
    }
    else
    {
        uint _S4371 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4371;
    }
    *&((&_S4366)->x) = v_3;
    bool _S4372 = (F32_isfinite((_S4366.y)));
    uint _S4373 = __ballot_sync(_S3205, _S4372);
    if(_S4372)
    {
        float _S4374 = _S4366.y;
        uint _S4375 = __ballot_sync(_S3205, true);
        v_3 = _S4374;
        _S3205 = _S4375;
    }
    else
    {
        uint _S4376 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4376;
    }
    *&((&_S4366)->y) = v_3;
    bool _S4377 = (F32_isfinite((_S4366.z)));
    uint _S4378 = __ballot_sync(_S3205, _S4377);
    if(_S4377)
    {
        float _S4379 = _S4366.z;
        uint _S4380 = __ballot_sync(_S3205, true);
        v_3 = _S4379;
        _S3205 = _S4380;
    }
    else
    {
        uint _S4381 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4381;
    }
    *&((&_S4366)->z) = v_3;
    float _S4382 = WaveActiveSum_0(_S4366.x, _S3205);
    *&((&_S4366)->x) = _S4382;
    float _S4383 = WaveActiveSum_0(_S4366.y, _S3205);
    *&((&_S4366)->y) = _S4383;
    float _S4384 = WaveActiveSum_0(_S4366.z, _S3205);
    *&((&_S4366)->z) = _S4384;
    uint warp_id_48 = thread_id_4 / 32U;
    uint _S4385 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_48] = _S4366.x;
        (*&_sh_block_reduce_shared_0)[warp_id_48 + 16U] = _S4366.y;
        (*&_sh_block_reduce_shared_0)[warp_id_48 + 32U] = _S4366.z;
        uint _S4386 = __ballot_sync(_S3205, true);
        _S3205 = _S4386;
    }
    else
    {
        uint _S4387 = __ballot_sync(_S3205, true);
        _S3205 = _S4387;
    }
    __syncthreads();
    bool _S4388 = warp_id_48 < 3U;
    uint _S4389 = __ballot_sync(_S3205, _S4388);
    if(_S4388)
    {
        bool _S4390 = lane_id_3 < 16U;
        uint _S4391 = __ballot_sync(_S4389, _S4390);
        if(_S4390)
        {
            float _S4392 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_48 * 16U];
            uint _S4393 = __ballot_sync(_S4389, true);
            v_3 = _S4392;
            _S3228 = _S4393;
        }
        else
        {
            uint _S4394 = __ballot_sync(_S4389, true);
            v_3 = 0.0f;
            _S3228 = _S4394;
        }
        float _S4395 = WaveActiveSum_0(v_3, _S3228);
        uint _S4396 = __ballot_sync(_S3228, _S3222);
        if(_S3222)
        {
            bool _S4397 = _S4395 != 0.0f;
            uint _S4398 = __ballot_sync(_S3228, true);
            _S3229 = _S4397;
            _S3228 = _S4398;
        }
        else
        {
            uint _S4399 = __ballot_sync(_S3228, true);
            _S3229 = false;
            _S3228 = _S4399;
        }
        uint _S4400 = __ballot_sync(_S3228, _S3229);
        if(_S3229)
        {
            bool _S4401 = warp_id_48 == 0U;
            uint _S4402 = __ballot_sync(_S4400, _S4401);
            if(_S4401)
            {
                float _S4403 = atomicAdd(&((v_coeffs_34 + 22U)->x), _S4395);
                uint _S4404 = __ballot_sync(_S4400, true);
            }
            else
            {
                uint _S4405 = _S4400 & (~_S4402);
                bool _S4406 = warp_id_48 == 1U;
                uint _S4407 = __ballot_sync(_S4405, _S4406);
                if(_S4406)
                {
                    float _S4408 = atomicAdd(&((v_coeffs_34 + 22U)->y), _S4395);
                    uint _S4409 = __ballot_sync(_S4405, true);
                }
                else
                {
                    float _S4410 = atomicAdd(&((v_coeffs_34 + 22U)->z), _S4395);
                    uint _S4411 = __ballot_sync(_S4405, true);
                }
                uint _S4412 = __ballot_sync(_S4400, true);
            }
            uint _S4413 = __ballot_sync(_S3228, true);
        }
        else
        {
            uint _S4414 = __ballot_sync(_S3228, true);
        }
        uint _S4415 = __ballot_sync(_S3205, true);
        _S3205 = _S4415;
    }
    else
    {
        uint _S4416 = __ballot_sync(_S3205, true);
        _S3205 = _S4416;
    }
    float3  _S4417 = make_float3 (pSH24_7) * v_colors_34;
    float3  _S4418 = _S4417;
    bool _S4419 = (F32_isfinite((_S4417.x)));
    uint _S4420 = __ballot_sync(_S3205, _S4419);
    if(_S4419)
    {
        float _S4421 = _S4418.x;
        uint _S4422 = __ballot_sync(_S3205, true);
        v_3 = _S4421;
        _S3205 = _S4422;
    }
    else
    {
        uint _S4423 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4423;
    }
    *&((&_S4418)->x) = v_3;
    bool _S4424 = (F32_isfinite((_S4418.y)));
    uint _S4425 = __ballot_sync(_S3205, _S4424);
    if(_S4424)
    {
        float _S4426 = _S4418.y;
        uint _S4427 = __ballot_sync(_S3205, true);
        v_3 = _S4426;
        _S3205 = _S4427;
    }
    else
    {
        uint _S4428 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4428;
    }
    *&((&_S4418)->y) = v_3;
    bool _S4429 = (F32_isfinite((_S4418.z)));
    uint _S4430 = __ballot_sync(_S3205, _S4429);
    if(_S4429)
    {
        float _S4431 = _S4418.z;
        uint _S4432 = __ballot_sync(_S3205, true);
        v_3 = _S4431;
        _S3205 = _S4432;
    }
    else
    {
        uint _S4433 = __ballot_sync(_S3205, true);
        v_3 = 0.0f;
        _S3205 = _S4433;
    }
    *&((&_S4418)->z) = v_3;
    float _S4434 = WaveActiveSum_0(_S4418.x, _S3205);
    *&((&_S4418)->x) = _S4434;
    float _S4435 = WaveActiveSum_0(_S4418.y, _S3205);
    *&((&_S4418)->y) = _S4435;
    float _S4436 = WaveActiveSum_0(_S4418.z, _S3205);
    *&((&_S4418)->z) = _S4436;
    uint warp_id_49 = thread_id_4 / 32U;
    uint _S4437 = __ballot_sync(_S3205, _S3222);
    if(_S3222)
    {
        (*&_sh_block_reduce_shared_0)[warp_id_49] = _S4418.x;
        (*&_sh_block_reduce_shared_0)[warp_id_49 + 16U] = _S4418.y;
        (*&_sh_block_reduce_shared_0)[warp_id_49 + 32U] = _S4418.z;
        uint _S4438 = __ballot_sync(_S3205, true);
        _S3205 = _S4438;
    }
    else
    {
        uint _S4439 = __ballot_sync(_S3205, true);
        _S3205 = _S4439;
    }
    __syncthreads();
    bool _S4440 = warp_id_49 < 3U;
    uint _S4441 = __ballot_sync(_S3205, _S4440);
    if(_S4440)
    {
        bool _S4442 = lane_id_3 < 16U;
        uint _S4443 = __ballot_sync(_S4441, _S4442);
        if(_S4442)
        {
            float _S4444 = (*&_sh_block_reduce_shared_0)[lane_id_3 + warp_id_49 * 16U];
            uint _S4445 = __ballot_sync(_S4441, true);
            v_3 = _S4444;
            _S3205 = _S4445;
        }
        else
        {
            uint _S4446 = __ballot_sync(_S4441, true);
            v_3 = 0.0f;
            _S3205 = _S4446;
        }
        float _S4447 = WaveActiveSum_0(v_3, _S3205);
        uint _S4448 = __ballot_sync(_S3205, _S3222);
        if(_S3222)
        {
            _S3229 = _S4447 != 0.0f;
        }
        else
        {
            _S3229 = false;
        }
        if(_S3229)
        {
            if(warp_id_49 == 0U)
            {
                float _S4449 = atomicAdd(&((v_coeffs_34 + 23U)->x), _S4447);
            }
            else
            {
                if(warp_id_49 == 1U)
                {
                    float _S4450 = atomicAdd(&((v_coeffs_34 + 23U)->y), _S4447);
                }
                else
                {
                    float _S4451 = atomicAdd(&((v_coeffs_34 + 23U)->z), _S4447);
                }
            }
        }
    }
    float fTmp0D_z_6 = -14.04997730255126953f * z2_19 + 2.00713968276977539f;
    float fTmp1C_z_6 = 6.62322282791137695f * z_37;
    float pSH22_x_6 = fTmp1C_9 * _S3364;
    float3  * _S4452 = coeffs_49 + int(15);
    float3  * _S4453 = coeffs_49 + int(23);
    float3  * _S4454 = coeffs_49 + int(16);
    float3  * _S4455 = coeffs_49 + int(22);
    float3  * _S4456 = coeffs_49 + int(17);
    float3  * _S4457 = coeffs_49 + int(21);
    float3  * _S4458 = coeffs_49 + int(20);
    float3  * _S4459 = coeffs_49 + int(18);
    float3  dir_n_27 = make_float3 (x_40, y_38, z_37);
    float3  v_dir_n_43 = make_float3 (v_x_32 + dot_0(v_colors_34, make_float3 (0.62583571672439575f * (fS2_9 + y_38 * fC2_x_6 + x_40 * fS2_x_6)) * *_S4452 + make_float3 (0.62583571672439575f * (fC2_9 + x_40 * fC2_x_6 - y_38 * fS2_x_6)) * *_S4453 + make_float3 (fTmp2B_9 * fS2_x_6) * *_S4454 + make_float3 (fTmp2B_9 * fC2_x_6) * *_S4455 + make_float3 (fTmp1C_9 * fS1_x_13) * *_S4456 + make_float3 (pSH22_x_6) * *_S4457 + make_float3 (fTmp0D_9) * *_S4458), v_y_32 + dot_0(v_colors_34, make_float3 (0.62583571672439575f * (x_40 * fS2_y_6 + fC2_9 + y_38 * fC2_y_6)) * *_S4452 + make_float3 (0.62583571672439575f * (x_40 * fC2_y_6 - fS2_9 - y_38 * fS2_y_6)) * *_S4453 + make_float3 (fTmp2B_9 * fS2_y_6) * *_S4454 + make_float3 (fTmp2B_9 * fC2_y_6) * *_S4455 + make_float3 (pSH22_x_6) * *_S4456 + make_float3 (fTmp1C_9 * fC1_y_13) * *_S4457 + make_float3 (fTmp0D_9) * *_S4459), v_z_32 + dot_0(v_colors_34, make_float3 (1.9843134880065918f * (pSH12_17 + z_37 * pSH12_z_8) + -1.00623059272766113f * pSH6_z_10) * *(coeffs_49 + int(19)) + make_float3 (fTmp0D_z_6 * x_40) * *_S4458 + make_float3 (fTmp0D_z_6 * y_38) * *_S4459 + make_float3 (fTmp1C_z_6 * fC1_19) * *_S4457 + make_float3 (fTmp1C_z_6 * fS1_19) * *_S4456 + make_float3 (-1.77013075351715088f * fC2_9) * *_S4455 + make_float3 (-1.77013075351715088f * fS2_9) * *_S4454));
    *v_dir_14 = *v_dir_14 + (v_dir_n_43 - make_float3 (dot_0(v_dir_n_43, dir_n_27)) * dir_n_27) * make_float3 (inorm_11);
    return;
}

