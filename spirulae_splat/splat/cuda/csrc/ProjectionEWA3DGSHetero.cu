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
    int32_t camera_idx = upper_bound(intersection_count_map, C+1, thread_idx) - 1;
    int32_t gauss_idx = intersection_splat_id[thread_idx];

    bool valid = (thread_idx < nnz);

    // check if points are with camera near and far plane
    glm::vec3 mean_c;
    glm::mat3 R;
    if (valid) {
        // shift pointers to the current camera and gaussian
        means += gauss_idx * 3;
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
        
        quats += gauss_idx * 4;
        scales += gauss_idx * 3;
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
            float opacity = opacities[gauss_idx];
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

    if (thread_idx < nnz) {
        // write to outputs
        camera_ids[thread_idx] = camera_idx;
        gaussian_ids[thread_idx] = gauss_idx;
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
    const at::optional<at::Tensor> opacities, // [..., N]
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
