// Engine forward pass + debug forward.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include <stdexcept>


void forward_3dgs(
    std::string primitive,   // "3dgs", "mip", "3dgut"
    int sh_degree,
    bool packed
) {
    engine().primitive = primitive;
    engine().sh_degree = sh_degree;
    engine().packed = packed;

    // Build splats as DeviceTensorFloatND from typed device buffers.
    // For DeviceVector<float> (opacities), build a [N, 1] shape FND manually.
    DeviceTensorFloatND fnd_opac;
    {
        TorchTensorView opac_tv((uint64_t)engine().world.opacities.data_ptr(), 4,
            {engine().world.opacities.size(), 1LL, 1LL});
        fnd_opac = DeviceTensorFloatND(opac_tv);  // ndim=2, shape=[N, 1]
    }
    std::vector<DeviceTensorFloatND> in_splats = {
        DeviceTensorFloatND(engine().world.means),         // [N, 3]
        DeviceTensorFloatND(engine().world.quats),         // [N, 4]
        DeviceTensorFloatND(engine().world.scales),        // [N, 3]
        fnd_opac,                                          // [N, 1]
        DeviceTensorFloatND(engine().world.features_dc),   // [N, 3]
        DeviceTensorFloatND(engine().world.features_sh),   // [N, K, 3]
    };
    engine().fwd.splats_w = in_splats;

    DeviceTensorFloatND aabb_nd, depths_nd;

    // Allocate and zero radii buffer before projection (kernel uses atomicMax)
    engine().optim.radii.resize("eng.radii", engine().max_num_splats);
    engine().optim.radii.zero();

    // --- Projection ---
    if (packed) {
        DeviceVector<int32_t> cam_ids, gauss_ids;
        DeviceVector<float4> aabb_vec;
        DeviceVector<float> depths_vec;

        if (primitive == "3dgs") {
            auto [a, b, c, d, e] = projection_3dgs_packed_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            engine().fwd.splats_s = e;
        } else if (primitive == "mip") {
            auto [a, b, c, d, e] = projection_mip_packed_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            engine().fwd.splats_s = e;
        } else if (primitive == "3dgut") {
            auto [a, b, c, d, e] = projection_3dgut_packed_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            engine().fwd.splats_s = e;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        engine().fwd.camera_ids = cam_ids;
        engine().fwd.gaussian_ids = gauss_ids;
        engine().fwd.aabb = vec_to_2d_float4(aabb_vec);  // [nnz, 1] view for backward
        aabb_nd = DeviceTensorFloatND(aabb_vec);          // [nnz, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_vec);      // [nnz]   for intersect
    } else {
        DeviceTensor2D<float4> aabb_2d;
        DeviceTensor2D<float> depths_2d;  // [C, N] — sorting depths per camera

        if (primitive == "3dgs") {
            auto [a, b, c] = projection_3dgs_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            aabb_2d = a; depths_2d = b;
            engine().fwd.splats_s = c;
        } else if (primitive == "mip") {
            auto [a, b, c] = projection_mip_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            aabb_2d = a; depths_2d = b;
            engine().fwd.splats_s = c;
        } else if (primitive == "3dgut") {
            auto [a, b, c] = projection_3dgut_forward(
                engine().cur_num_splats, sh_degree, in_splats,
                _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
                (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
                engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
                engine().optim.radii);
            aabb_2d = a; depths_2d = b;
            engine().fwd.splats_s = c;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        engine().fwd.camera_ids = DeviceVector<int32_t>();
        engine().fwd.gaussian_ids = DeviceVector<int32_t>();
        engine().fwd.aabb = aabb_2d;                           // [C, N] for backward
        aabb_nd = DeviceTensorFloatND(aabb_2d);                // [C, N, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_2d);            // [C, N] -> numel=C*N for intersect
    }

    // --- Tile intersection (AABB mode) ---
    DeviceVector<int32_t>* image_ids_ptr = nullptr;
    if (packed && engine().fwd.camera_ids.data_ptr() != nullptr)
        image_ids_ptr = &engine().fwd.camera_ids;

    auto [isect_ids, flatten_ids, tile_offsets] = do_intersect_tile_generic(
        aabb_nd, depths_nd,
        nullptr, nullptr, nullptr,
        (uint32_t)engine().camera.num,
        _dv_tv(engine().camera.intrins),
        (uint32_t)engine().camera.width,
        (uint32_t)engine().camera.height,
        image_ids_ptr
    );

    engine().fwd.tile_offsets = tile_offsets;
    engine().fwd.flatten_ids = flatten_ids;

    // --- Rasterization forward ---
    RenderOutput::TensorTuple renders;
    DeviceTensor3D<float> render_Ts;
    DeviceTensor3D<int32_t> last_ids;

    if (primitive == "3dgs") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgs_fwd(
            in_splats, engine().fwd.splats_s, engine().fwd.gaussian_ids,
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "mip") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_mip_fwd(
            in_splats, engine().fwd.splats_s, engine().fwd.gaussian_ids,
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "3dgut") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgut_fwd(
            in_splats, engine().fwd.splats_s, engine().fwd.gaussian_ids,
            _dt2d_tv(engine().camera.viewmats), _dv_tv(engine().camera.intrins),
            engine().camera.model_str, _dt2d_tv(engine().camera.dist_coeffs),
            (uint32_t)engine().camera.width, (uint32_t)engine().camera.height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    }

    engine().fwd.renders = renders;
    engine().fwd.render_Ts = render_Ts;
    engine().fwd.last_ids = last_ids;

    // Results stay in pool — use engine_copy_render_to_host to fetch
}


void engine_debug_forward(
    TorchTensorView override_features_dc,
    int override_sh_degree,
    TorchTensorView out_rgb
) {
    if (std::get<0>(engine().fwd.renders).data_ptr() == nullptr)
        throw std::runtime_error("engine_debug_forward: forward_3dgs must be called first");

    // Swap in overrides (H->D copy for host tensor)
    DeviceVector<float3> saved_dc = engine().world.features_dc;
    int saved_sh = engine().sh_degree;
    if (std::get<0>(override_features_dc) != 0)
        engine().world.features_dc = _hv_to_dv<float3>("debug.features_dc", override_features_dc);
    if (override_sh_degree >= 0)
        engine().sh_degree = override_sh_degree;

    forward_3dgs(engine().primitive, engine().sh_degree, engine().packed);

    // Copy rgb result D->H
    auto& rgb = std::get<0>(engine().fwd.renders);
    if (rgb.data_ptr() && std::get<0>(out_rgb) != 0) {
        cudaMemcpy((void*)std::get<0>(out_rgb), rgb.data_ptr(),
                   rgb.numel() * sizeof(float3), cudaMemcpyDeviceToHost);
    }

    // Restore
    engine().world.features_dc = saved_dc;
    engine().sh_degree = saved_sh;
}
