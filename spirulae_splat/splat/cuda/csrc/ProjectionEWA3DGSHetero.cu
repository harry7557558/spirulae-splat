#include "ProjectionEWA3DGSHetero.cuh"

#include "helpers.cuh"
#include <algorithm>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>


template<typename T>
__device__ __forceinline__ unsigned upper_bound(
    const T *arr, unsigned n, T value
) {
    unsigned left = 0, right = n;
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] <= value)
            left = mid + 1;
        else
            right = mid;
    }
    return left;
}

__global__ void projection_ewa_3dgs_hetero_forward_kernel(
    const uint32_t C,
    const uint32_t nnz,
    const float *__restrict__ means,    // [N, 3]
    const float *__restrict__ quats,    // [N, 4]
    const float *__restrict__ scales,   // [N, 3]
    const float *__restrict__ opacities, // [N]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float *__restrict__ Ks,       // [C, 3, 3]
    const uint32_t image_width,  // TILE_SIZE
    const uint32_t image_height,  // TILE_SIZE
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    const gsplat::CameraModelType camera_model,
    const int32_t* __restrict__ intersection_count_map,  // [C+1]
    const int32_t* __restrict__ intersection_splat_id,  // [nnz]
    // outputs
    int64_t *__restrict__ camera_ids,    // [nnz]
    int64_t *__restrict__ gaussian_ids,  // [nnz]
    int32_t *__restrict__ radii,         // [nnz, 2]
    float *__restrict__ means2d,      // [nnz, 2]
    float *__restrict__ depths,       // [nnz]
    float *__restrict__ conics,       // [nnz, 3]
    float *__restrict__ compensations // [nnz] optional
) {
    int32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= nnz)
        return;
    int32_t camera_idx = upper_bound(intersection_count_map, C+1, thread_idx) - 1;
    int32_t splat_idx = intersection_splat_id[thread_idx];

    bool valid = true;

    // check if points are with camera near and far plane
    glm::vec3 mean_c;
    glm::mat3 R;
    if (valid) {
        // shift pointers to the current camera and gaussian
        means += splat_idx * 3;
        viewmats += camera_idx * 16;

        // glm is column-major but input is row-major
        R = glm::mat3(
            viewmats[0],
            viewmats[4],
            viewmats[8], // 1st column
            viewmats[1],
            viewmats[5],
            viewmats[9], // 2nd column
            viewmats[2],
            viewmats[6],
            viewmats[10] // 3rd column
        );
        glm::vec3 t = { viewmats[3], viewmats[7], viewmats[11] };

        // transform Gaussian center to camera space
        gsplat::posW2C(R, t, glm::make_vec3(means), mean_c);
        if (mean_c.z < near_plane || mean_c.z > far_plane) {
            valid = false;
        }
    }

    // check if the perspective projection is valid.
    glm::mat2 covar2d;
    glm::vec2 mean2d = glm::vec2(0);
    glm::mat2 covar2d_inv = glm::mat2(0);
    float compensation = 0.0f;
    float det;
    if (valid) {
        // transform Gaussian covariance to camera space
        glm::mat3 covar;
        
        quats += splat_idx * 4;
        scales += splat_idx * 3;
        gsplat::quat_scale_to_covar_preci(
            glm::make_vec4(quats), glm::make_vec3(scales), &covar, nullptr
        );

        glm::mat3 covar_c;
        gsplat::covarW2C(R, covar, covar_c);

        Ks += camera_idx * 9;
        switch (camera_model) {
        case gsplat::CameraModelType::PINHOLE: // perspective projection
            gsplat::persp_proj(
                mean_c,
                covar_c,
                Ks[0],
                Ks[4],
                Ks[2],
                Ks[5],
                covar2d,
                mean2d
            );
            break;
        case gsplat::CameraModelType::ORTHO: // orthographic projection
            gsplat::ortho_proj(
                mean_c,
                covar_c,
                Ks[0],
                Ks[4],
                Ks[2],
                Ks[5],
                covar2d,
                mean2d
            );
            break;
        case gsplat::CameraModelType::FISHEYE: // fisheye projection
            gsplat::fisheye_proj(
                mean_c,
                covar_c,
                Ks[0],
                Ks[4],
                Ks[2],
                Ks[5],
                covar2d,
                mean2d
            );
            break;
        }

        det = gsplat::add_blur(eps2d, covar2d, compensation);
        if (det <= 0.f) {
            valid = false;
        } else {
            // compute the inverse of the 2d covariance
            covar2d_inv = glm::inverse(covar2d);
        }
    }

    // check if the points are in the image region
    float radius_x, radius_y;
    if (valid) {
        float extend = 3.33f;
        if (opacities != nullptr) {
            float opacity = opacities[splat_idx];
            if (compensations != nullptr) {
                // we assume compensation term will be applied later on.
                opacity *= compensation;
            }
            if (opacity < ALPHA_THRESHOLD) {
                valid = false;
            }
            // Compute opacity-aware bounding box.
            // https://arxiv.org/pdf/2402.00525 Section B.2
            extend = min(extend, sqrt(2.0f * __logf(opacity / ALPHA_THRESHOLD)));
        }
        
        // compute tight rectangular bounding box (non differentiable)
        // https://arxiv.org/pdf/2402.00525
        radius_x = ceilf(extend * sqrtf(covar2d[0][0]));
        radius_y = ceilf(extend * sqrtf(covar2d[1][1]));
        
        if (radius_x <= radius_clip && radius_y <= radius_clip) {
            valid = false;
        }

        // mask out gaussians outside the image region
        if (mean2d.x + radius_x <= 0 || mean2d.x - radius_x >= image_width ||
            mean2d.y + radius_y <= 0 || mean2d.y - radius_y >= image_height) {
            valid = false;
        }
    }

    {
        // write to outputs
        camera_ids[thread_idx] = camera_idx;
        gaussian_ids[thread_idx] = splat_idx;
        radii[thread_idx * 2] = (int32_t)radius_x * int(valid);
        radii[thread_idx * 2 + 1] = (int32_t)radius_y * int(valid);
        means2d[thread_idx * 2] = mean2d.x;
        means2d[thread_idx * 2 + 1] = mean2d.y;
        depths[thread_idx] = valid ? mean_c.z : -0.0f;
        conics[thread_idx * 3] = covar2d_inv[0][0];
        conics[thread_idx * 3 + 1] = covar2d_inv[0][1];
        conics[thread_idx * 3 + 2] = covar2d_inv[1][1];
        if (compensations != nullptr) {
            compensations[thread_idx] = compensation;
        }
    }
}


std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // radii
    at::Tensor,  // means2d
    at::Tensor,  // depths
    at::Tensor,  // conics
    at::Tensor  // compensations
> projection_ewa_3dgs_hetero_forward_tensor(
    // inputs
    const at::Tensor means,                // [..., N, 3]
    const at::Tensor quats,  // [..., N, 4]
    const at::Tensor scales, // [..., N, 3]
    const std::optional<at::Tensor> opacities, // [..., N]
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    const bool calc_compensations,
    const gsplat::CameraModelType camera_model,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    uint32_t N = means.size(-2);          // number of gaussians
    uint32_t C = viewmats.size(-3);       // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = means.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor radii = at::empty({nnz, 2}, opt.dtype(at::kInt));
    at::Tensor means2d = at::empty({nnz, 2}, opt);
    at::Tensor depths = at::empty({nnz}, opt);
    at::Tensor conics = at::empty({nnz, 3}, opt);
    at::Tensor compensations;
    if (calc_compensations) {
        compensations = at::empty({nnz}, opt);
    }

    constexpr uint block = 256;
    projection_ewa_3dgs_hetero_forward_kernel<<<_CEIL_DIV(nnz, block), block>>>(
        C,
        nnz,
        means.data_ptr<float>(),
        quats.data_ptr<float>(),
        scales.data_ptr<float>(),
        opacities.has_value() ? opacities.value().data_ptr<float>() : nullptr,
        viewmats.data_ptr<float>(),
        Ks.data_ptr<float>(),
        image_width,
        image_height,
        eps2d,
        near_plane,
        far_plane,
        radius_clip,
        camera_model,
        intersection_count_map.data_ptr<int32_t>(),
        intersection_splat_id.data_ptr<int32_t>(),
        camera_ids.data_ptr<int64_t>(),
        gaussian_ids.data_ptr<int64_t>(),
        radii.data_ptr<int32_t>(),
        means2d.data_ptr<float>(),
        depths.data_ptr<float>(),
        conics.data_ptr<float>(),
        calc_compensations ? compensations.data_ptr<float>() : nullptr
    );

    return std::make_tuple(
        camera_ids,
        gaussian_ids,
        radii,
        means2d,
        depths,
        conics,
        compensations
    );
}




__global__ void projection_ewa_3dgs_hetero_backward_kernel(
    // fwd inputs
    const uint32_t C,
    const uint32_t nnz,
    const float *__restrict__ means,    // [N, 3]
    const float *__restrict__ quats,    // [N, 4]
    const float *__restrict__ scales,   // [N, 3]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float *__restrict__ Ks,       // [C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const gsplat::CameraModelType camera_model,
    // fwd outputs
    const int64_t *__restrict__ camera_ids,     // [nnz]
    const int64_t *__restrict__ gaussian_ids,   // [nnz]
    const float *__restrict__ conics,        // [nnz, 3]
    const float *__restrict__ compensations, // [nnz] optional
    // grad outputs
    const float *__restrict__ v_means2d,       // [nnz, 2]
    const float *__restrict__ v_depths,        // [nnz]
    const float *__restrict__ v_conics,        // [nnz, 3]
    const float *__restrict__ v_compensations, // [nnz] optional
    const bool sparse_grad, // whether the outputs are in COO format [nnz, ...]
    // grad inputs
    float *__restrict__ v_means,   // [N, 3] or [nnz, 3]
    float *__restrict__ v_quats,   // [N, 4] or [nnz, 4]
    float *__restrict__ v_scales,  // [N, 3] or [nnz, 3]
    float *__restrict__ v_viewmats // [C, 4, 4]
) {
    // parallelize over nnz.
    int32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= nnz)
        return;
    int32_t camera_idx = camera_ids[thread_idx];
    int32_t splat_idx = gaussian_ids[thread_idx];

    // shift pointers to the current camera and gaussian
    means += splat_idx * 3;
    viewmats += camera_idx * 16;
    Ks += camera_idx * 9;

    conics += thread_idx * 3;

    v_means2d += thread_idx * 2;
    v_depths += thread_idx;
    v_conics += thread_idx * 3;

    // vjp: compute the inverse of the 2d covariance
    glm::mat2 covar2d_inv(conics[0], conics[1], conics[1], conics[2]);
    glm::mat2 v_covar2d_inv(v_conics[0], v_conics[1] * .5f, v_conics[1] * .5f, v_conics[2]);
    glm::mat2 v_covar2d(0.f);
    gsplat::inverse_vjp(covar2d_inv, v_covar2d_inv, v_covar2d);

    if (v_compensations != nullptr) {
        // vjp: compensation term
        const float compensation = compensations[thread_idx];
        const float v_compensation = v_compensations[thread_idx];
        gsplat::add_blur_vjp(
            eps2d, covar2d_inv, compensation, v_compensation, v_covar2d
        );
    }

    // transform Gaussian to camera space
   glm:: mat3 R(
        viewmats[0],
        viewmats[4],
        viewmats[8], // 1st column
        viewmats[1],
        viewmats[5],
        viewmats[9], // 2nd column
        viewmats[2],
        viewmats[6],
        viewmats[10] // 3rd column
    );
    glm::vec3 t(viewmats[3], viewmats[7], viewmats[11]);
    glm::mat3 covar;
    glm::vec4 quat;
    glm::vec3 scale;
    {
        // compute it from quaternions and scales
        quat = glm::make_vec4(quats + splat_idx * 4);
        scale = glm::make_vec3(scales + splat_idx * 3);
        gsplat::quat_scale_to_covar_preci(quat, scale, &covar, nullptr);
    }
    glm::vec3 mean_c;
    gsplat::posW2C(R, t, glm::make_vec3(means), mean_c);
    glm::mat3 covar_c;
    gsplat::covarW2C(R, covar, covar_c);

    float fx = Ks[0], cx = Ks[2], fy = Ks[4], cy = Ks[5];
    glm::mat3 v_covar_c(0.f);
    glm::vec3 v_mean_c(0.f);
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        gsplat::persp_proj_vjp(
            mean_c,
            covar_c,
            fx,
            fy,
            cx,
            cy,
            v_covar2d,
            glm::make_vec2(v_means2d),
            v_mean_c,
            v_covar_c
        );
        break;
    case gsplat::CameraModelType::ORTHO: // orthographic projection
        gsplat::ortho_proj_vjp(
            mean_c,
            covar_c,
            fx,
            fy,
            cx,
            cy,
            v_covar2d,
            glm::make_vec2(v_means2d),
            v_mean_c,
            v_covar_c
        );
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        gsplat::fisheye_proj_vjp(
            mean_c,
            covar_c,
            fx,
            fy,
            cx,
            cy,
            v_covar2d,
            glm::make_vec2(v_means2d),
            v_mean_c,
            v_covar_c
        );
        break;
    }

    // add contribution from v_depths
    v_mean_c.z += v_depths[0];

    // vjp: transform Gaussian covariance to camera space
    glm::vec3 v_mean(0.f);
    glm::mat3 v_covar(0.f);
    glm::mat3 v_R(0.f);
    glm::vec3 v_t(0.f);
    gsplat::posW2C_VJP(R, t, glm::make_vec3(means), v_mean_c, v_R, v_t, v_mean);
    gsplat::covarW2C_VJP(R, covar, v_covar_c, v_R, v_covar);

    auto warp = cg::tiled_partition<32>(cg::this_thread_block());
    if (sparse_grad) {
        // write out results with sparse layout
        if (v_means != nullptr) {
            v_means += thread_idx * 3;
            #pragma unroll
            for (uint32_t i = 0; i < 3; i++) {
                v_means[i] = v_mean[i];
            }
        }
        {
            glm::mat3 rotmat = gsplat::quat_to_rotmat(quat);
            glm::vec4 v_quat(0.f);
            glm::vec3 v_scale(0.f);
            gsplat::quat_scale_to_covar_vjp(
                quat, scale, rotmat, v_covar, v_quat, v_scale
            );
            v_quats += thread_idx * 4;
            v_scales += thread_idx * 3;
            v_quats[0] = v_quat[0];
            v_quats[1] = v_quat[1];
            v_quats[2] = v_quat[2];
            v_quats[3] = v_quat[3];
            v_scales[0] = v_scale[0];
            v_scales[1] = v_scale[1];
            v_scales[2] = v_scale[2];
        }
    } else {
        // write out results with dense layout
        // #if __CUDA_ARCH__ >= 700
        // write out results with warp-level reduction
        auto warp_group_g = cg::labeled_partition(warp, splat_idx);
        if (v_means != nullptr) {
            warpSum(v_mean, warp_group_g);
            if (warp_group_g.thread_rank() == 0) {
                v_means += splat_idx * 3;
                #pragma unroll
                for (uint32_t i = 0; i < 3; i++) {
                    atomicAdd(v_means + i, v_mean[i]);
                }
            }
        }
        {
            // Directly output gradients w.r.t. the quaternion and scale
            glm::mat3 rotmat = gsplat::quat_to_rotmat(quat);
            glm::vec4 v_quat(0.f);
            glm::vec3 v_scale(0.f);
            gsplat::quat_scale_to_covar_vjp(
                quat, scale, rotmat, v_covar, v_quat, v_scale
            );
            warpSum(v_quat, warp_group_g);
            warpSum(v_scale, warp_group_g);
            if (warp_group_g.thread_rank() == 0) {
                v_quats += splat_idx * 4;
                v_scales += splat_idx * 3;
                atomicAdd(v_quats, v_quat[0]);
                atomicAdd(v_quats + 1, v_quat[1]);
                atomicAdd(v_quats + 2, v_quat[2]);
                atomicAdd(v_quats + 3, v_quat[3]);
                atomicAdd(v_scales, v_scale[0]);
                atomicAdd(v_scales + 1, v_scale[1]);
                atomicAdd(v_scales + 2, v_scale[2]);
            }
        }
    }
    // v_viewmats is always in dense layout
    if (v_viewmats != nullptr) {
        auto warp_group_c = cg::labeled_partition(warp, camera_idx);
        warpSum(v_R, warp_group_c);
        warpSum(v_t, warp_group_c);
        if (warp_group_c.thread_rank() == 0) {
            v_viewmats += camera_idx * 16;
            #pragma unroll
            for (uint32_t i = 0; i < 3; i++) { // rows
                #pragma unroll
                for (uint32_t j = 0; j < 3; j++) { // cols
                    atomicAdd(v_viewmats + i * 4 + j, v_R[j][i]);
                }
                atomicAdd(v_viewmats + i * 4 + 3, v_t[i]);
            }
        }
    }
}



std::tuple<
    at::Tensor,  // v_means
    at::Tensor,  // v_quats
    at::Tensor,  // v_scales
    at::Tensor  // v_viewmats
> projection_ewa_3dgs_hetero_backward_tensor(
    // fwd inputs
    const at::Tensor means, // [..., N, 3]
    const at::Tensor quats, // [..., N, 4]
    const at::Tensor scales, // [..., N, 3]
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const gsplat::CameraModelType camera_model,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor conics, // [nnz, 3]
    const std::optional<at::Tensor> compensations, // [nnz] optional
    // grad outputs
    const at::Tensor v_means2d, // [nnz, 2]
    const at::Tensor v_depths, // [nnz]
    const at::Tensor v_conics, // [nnz, 3]
    const std::optional<at::Tensor> v_compensations, // [nnz] optional
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    uint32_t N = means.size(-2);          // number of gaussians
    uint32_t C = viewmats.size(-3);       // number of cameras
    uint32_t nnz = camera_ids.size(0);

    auto opt = means.options();
    at::Tensor v_means, v_quats, v_scales, v_viewmats;
    if (sparse_grad) {
        v_means = at::zeros({nnz, 3}, opt);
        v_quats = at::zeros({nnz, 4}, opt);
        v_scales = at::zeros({nnz, 3}, opt);
    } else {
        v_means = at::zeros_like(means);
        v_quats = at::zeros_like(quats, opt);
        v_scales = at::zeros_like(scales, opt);
    }
    if (viewmats_requires_grad) {
        v_viewmats = at::zeros_like(viewmats, opt);
    }

    constexpr uint block = 256;
    projection_ewa_3dgs_hetero_backward_kernel<<<_CEIL_DIV(nnz, block), block>>>(
        C,
        nnz,
        means.data_ptr<float>(),
        quats.data_ptr<float>(),
        scales.data_ptr<float>(),
        viewmats.data_ptr<float>(),
        Ks.data_ptr<float>(),
        image_width,
        image_height,
        eps2d,
        camera_model,
        camera_ids.data_ptr<int64_t>(),
        gaussian_ids.data_ptr<int64_t>(),
        conics.data_ptr<float>(),
        compensations.has_value()
            ? compensations.value().data_ptr<float>()
            : nullptr,
        v_means2d.data_ptr<float>(),
        v_depths.data_ptr<float>(),
        v_conics.data_ptr<float>(),
        v_compensations.has_value()
            ? v_compensations.value().data_ptr<float>()
            : nullptr,
        sparse_grad,
        v_means.data_ptr<float>(),
        v_quats.data_ptr<float>(),
        v_scales.data_ptr<float>(),
        viewmats_requires_grad
            ? v_viewmats.data_ptr<float>()
            : nullptr
    );

    return std::make_tuple(v_means, v_quats, v_scales, v_viewmats);
}
