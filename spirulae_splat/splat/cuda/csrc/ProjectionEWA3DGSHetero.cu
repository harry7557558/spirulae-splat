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

template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_hetero_forward_kernel(
    const uint32_t C,
    const uint32_t nnz,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float *__restrict__ Ks,       // [C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const int32_t* __restrict__ intersection_count_map,  // [C+1]
    const int32_t* __restrict__ intersection_splat_id,  // [nnz]
    // outputs
    int64_t *__restrict__ camera_ids,    // [nnz]
    int64_t *__restrict__ gaussian_ids,  // [nnz]
    int4 *__restrict__ aabbs,    // [nnz, 4]
    typename SplatPrimitive::Screen::Buffer splats_proj
) {
    int32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= nnz)
        return;
    int32_t camera_idx = upper_bound(intersection_count_map, C+1, thread_idx) - 1;
    int32_t splat_idx = intersection_splat_id[thread_idx];

    // Load camera
    viewmats += camera_idx * 16;
    Ks += camera_idx * 9;
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    typename SplatPrimitive::FwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
        near_plane, far_plane,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(camera_idx);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, splat_idx);

    // Projection
    int4 aabb;
    typename SplatPrimitive::Screen splat_proj;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_proj, aabb);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye(splat_world, cam, splat_proj, aabb);
        break;
    }

    // Save results
    // aabb.x = min(max(aabb.x, 0), tile_width-1);
    // aabb.y = min(max(aabb.y, 0), tile_height-1);
    // aabb.z = min(max(aabb.z, 0), tile_width-1);
    // aabb.w = min(max(aabb.w, 0), tile_height-1);
    int offset_x = (int)roundf(0.5f * image_width - cx);
    int offset_y = (int)roundf(0.5f * image_height - cy);
    aabb.x = min(max(aabb.x + offset_x, 0), image_width-1) - offset_x;
    aabb.y = min(max(aabb.y + offset_y, 0), image_height-1) - offset_y;
    aabb.z = min(max(aabb.z + offset_x, 0), image_width-1) - offset_x;
    aabb.w = min(max(aabb.w + offset_y, 0), image_height-1) - offset_y;
    camera_ids[thread_idx] = camera_idx;
    gaussian_ids[thread_idx] = splat_idx;
    aabbs[thread_idx] = aabb;
    splat_proj.saveParamsToBuffer(splats_proj, thread_idx);
}


// TODO: refactor

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    Vanilla3DGS::Screen::TensorTuple  // out splats
> projection_3dgs_hetero_forward_tensor(
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    Vanilla3DGS::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kInt));
    Vanilla3DGS::Screen::Tensor splats_proj =
        Vanilla3DGS::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, near_plane, far_plane, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb,
        splats_proj.tupleProjFwdPacked()
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    MipSplatting::Screen::TensorTuple  // out splats
> projection_mip_hetero_forward_tensor(
    // inputs
    const MipSplatting::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    MipSplatting::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kInt));
    MipSplatting::Screen::Tensor splats_proj =
        MipSplatting::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, near_plane, far_plane, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<MipSplatting, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<MipSplatting, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb,
        splats_proj.tupleProjFwdPacked()
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    Vanilla3DGUT::Screen::TensorTuple  // out splats
> projection_3dgut_hetero_forward_tensor(
    // inputs
    const Vanilla3DGUT::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    Vanilla3DGUT::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kInt));
    Vanilla3DGUT::Screen::Tensor splats_proj =
        Vanilla3DGUT::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, near_plane, far_plane, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<Vanilla3DGUT, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<Vanilla3DGUT, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb,
        splats_proj.tupleProjFwdPacked()
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    OpaqueTriangle::Screen::TensorTuple  // out splats
> projection_opaque_triangle_hetero_forward_tensor(
    // inputs
    const OpaqueTriangle::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    OpaqueTriangle::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kInt));
    OpaqueTriangle::Screen::Tensor splats_proj =
        OpaqueTriangle::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, near_plane, far_plane, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<OpaqueTriangle, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<OpaqueTriangle, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb,
        splats_proj.tupleProjFwdPacked()
    );
}




/*[AutoHeaderGeneratorExport]*/
template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_3dgs_hetero_backward_kernel(
    // fwd inputs
    const long C,
    const long N,
    const uint32_t nnz,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float *__restrict__ Ks,       // [C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    // fwd outputs
    const int64_t *__restrict__ camera_ids,     // [nnz]
    const int64_t *__restrict__ gaussian_ids,   // [nnz]
    const int4 *__restrict__ aabbs,          // [B, C, N, 4]
    // grad outputs
    typename SplatPrimitive::Screen::Buffer v_splats_proj,
    const bool sparse_grad, // whether the outputs are in COO format [nnz, ...]
    // grad inputs
    typename SplatPrimitive::World::Buffer v_splats_world,
    float *__restrict__ v_viewmats // [C, 4, 4]
) {
    // parallelize over nnz.
    int32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= nnz)
        return;

    int4 aabb = aabbs[thread_idx];
    aabb.x = min(max(aabb.x, 0), tile_width-1);
    aabb.y = min(max(aabb.y, 0), tile_height-1);
    aabb.z = min(max(aabb.z, 0), tile_width-1);
    aabb.w = min(max(aabb.w, 0), tile_height-1);
    if ((aabb.z-aabb.x) * (aabb.w-aabb.y) <= 0)
        return;

    int32_t camera_idx = camera_ids[thread_idx];
    int32_t splat_idx = gaussian_ids[thread_idx];

    // Load camera
    viewmats += camera_idx * 16;
    Ks += camera_idx * 9;
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    typename SplatPrimitive::BwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(camera_idx);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, splat_idx);
    typename SplatPrimitive::Screen v_splat_proj =
        SplatPrimitive::Screen::load(v_splats_proj, thread_idx);

    // Projection
    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp_vjp(splat_world, cam, v_splat_proj, v_splat_world, v_R, v_t);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye_vjp(splat_world, cam, v_splat_proj, v_splat_world, v_R, v_t);
        break;
    }
    
    // Save results
    auto warp = cg::tiled_partition<32>(cg::this_thread_block());
    if (sparse_grad) {
        v_splat_world.saveParamsToBuffer(v_splats_world, thread_idx);
    } else {
        v_splat_world.atomicAddGradientToBuffer(v_splats_world, splat_idx);
    }
    if (v_viewmats != nullptr) {
        auto warp_group_c = cg::labeled_partition(warp, camera_idx);
        warpSum(v_R[0], warp_group_c);
        warpSum(v_R[1], warp_group_c);
        warpSum(v_R[2], warp_group_c);
        warpSum(v_t, warp_group_c);
        if (warp_group_c.thread_rank() == 0) {
            v_viewmats += camera_idx * 16;
            #pragma unroll
            for (uint32_t i = 0; i < 3; i++) { // rows
                atomicAdd(v_viewmats + i * 4 + 0, v_R[i].x);
                atomicAdd(v_viewmats + i * 4 + 1, v_R[i].y);
                atomicAdd(v_viewmats + i * 4 + 2, v_R[i].z);
            }
            atomicAdd(v_viewmats + 0 * 4 + 3, v_t.x);
            atomicAdd(v_viewmats + 1 * 4 + 3, v_t.y);
            atomicAdd(v_viewmats + 2 * 4 + 3, v_t.z);
        }
    }
}



// TODO: refactor

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    Vanilla3DGS::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    Vanilla3DGS::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    Vanilla3DGS::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), sparse_grad, v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_hetero_backward_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const MipSplatting::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    MipSplatting::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    MipSplatting::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    MipSplatting::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), sparse_grad, v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<MipSplatting, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<MipSplatting, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    Vanilla3DGUT::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    Vanilla3DGUT::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    Vanilla3DGUT::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), sparse_grad, v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGUT, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGUT, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    OpaqueTriangle::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_opaque_triangle_hetero_backward_tensor(
    // fwd inputs
    const OpaqueTriangle::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const OpaqueTriangle::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    OpaqueTriangle::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    OpaqueTriangle::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    OpaqueTriangle::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), sparse_grad, v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<OpaqueTriangle, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<OpaqueTriangle, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}
