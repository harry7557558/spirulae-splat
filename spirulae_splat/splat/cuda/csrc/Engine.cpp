#include "Engine.h"
#include "Tensor.h"

// TODO: write these into classes instead of using global variables

namespace Buffers {

/* 3DGS-like primitive in world-space */

int64_t num_splats;
int64_t max_num_splats;
int sh_degree;

enum class Primitive {
    Vanilla3DGS, MipSplatting, Vanilla3DGUT,
    Unknown,
} primitive;

DeviceTensor2D<float3> world_means;
DeviceTensor2D<float4> world_quats;
DeviceTensor2D<float3> world_scales;
DeviceTensor2D<float> world_opacities;
DeviceTensor2D<float3> world_features_dc;
DeviceTensor3D<float> world_features_sh;

/* Camera parameters */

int32_t num_cameras;
int32_t image_width;
int32_t image_height;
ssplat::CameraModelType camera_model;
DeviceTensor2D<float4> camera_viewmats;
DeviceVector<float4> camera_intrins;
DeviceTensor2D<float> camera_dist_coeffs;

/* 3DGS-like primitive in screen-space */

bool packed;

// DeviceTensor2D proj_means;
// DeviceTensor2D proj_quats;
// DeviceTensor2D proj_scales;
// DeviceTensor2D proj_opacities;
// DeviceTensor2D proj_features_dc;
// DeviceTensor2D proj_features_sh;

/* Forward outputs in image space */

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> renders;

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> renders2;

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> distortions;

inline DeviceTensor3D<float> render_Ts;
inline DeviceTensor3D<int32_t> render_last_ids;


/* Backward gradient in image space */

inline DeviceTensor3D<float> v_render_Ts;

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> v_renders;

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> v_renders2;

inline std::tuple<
    DeviceTensor3D<float3>,
    DeviceTensor3D<float>,
    DeviceTensor3D<float3>
> v_distortions;

/* Optimizer States */

DeviceTensor2D<float3> g1_world_means, g2_world_means;
DeviceTensor2D<float4> g1_world_quats, g2_world_quats;
DeviceTensor2D<float3> g1_world_scales, g2_world_scales;
DeviceTensor2D<float> g1_world_opacities, g2_world_opacities;
DeviceTensor2D<float3> g1_world_features_dc, g2_world_features_dc;
DeviceTensor3D<float> g1_world_features_sh, g2_world_features_sh;


}  // namespace Buffers


void setData3DGS(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    Buffers::num_splats = num_splats;

    Buffers::world_means = DeviceTensor2D<float3>(means);
    Buffers::world_quats = DeviceTensor2D<float4>(quats);
    Buffers::world_scales = DeviceTensor2D<float3>(scales);
    Buffers::world_opacities = DeviceTensor2D<float>(opacities);
    Buffers::world_features_dc = DeviceTensor2D<float3>(features_dc);
    Buffers::world_features_sh = DeviceTensor3D<float>(features_sh);

    Buffers::max_num_splats = std::get<2>(means)[0];
    if (std::get<2>(quats)[0] != Buffers::max_num_splats ||
        std::get<2>(scales)[0] != Buffers::max_num_splats ||
        std::get<2>(opacities)[0] != Buffers::max_num_splats ||
        std::get<2>(features_dc)[0] != Buffers::max_num_splats ||
        std::get<2>(features_sh)[0] != Buffers::max_num_splats)
        throw std::runtime_error("max_num_splats in quats/scales/opacities/features_dc/features_sh tensor is different from max_num_splats in means tensor");
    if (Buffers::max_num_splats < num_splats)
        throw std::runtime_error("max_num_splats in means tensor is smaller than num_splats");

}


void setCameraParams(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
) {
    Buffers::image_width = width;
    Buffers::image_height = height;
    Buffers::camera_model = cmt(camera_model);

    Buffers::camera_viewmats = DeviceTensor2D<float4>(viewmats);
    Buffers::camera_intrins = DeviceVector<float4>(intrins);
    Buffers::camera_dist_coeffs = DeviceTensor2D<float>(dist_coeffs);

    Buffers::num_cameras = std::get<2>(viewmats)[0];
    if (std::get<2>(intrins)[0] != Buffers::num_cameras ||
        std::get<2>(dist_coeffs)[0] != Buffers::num_cameras)
        throw std::runtime_error("num_cameras in intrins/dist_coeffs tensor is different from num_cameras in viewmats tensor");
}


void forward(
    std::string primitive,
    int sh_degree,
    bool packed
) {
    Buffers::primitive =
        primitive == "3dgs" ? Buffers::Primitive::Vanilla3DGS :
        primitive == "mip" ? Buffers::Primitive::MipSplatting :
        primitive == "3dgut" ? Buffers::Primitive::Vanilla3DGUT :
        Buffers::Primitive::Unknown;
    if (Buffers::primitive == Buffers::Primitive::Unknown)
        throw std::runtime_error("Unknown 3DGS primitive");

    Buffers::sh_degree = sh_degree;
    Buffers::packed = packed;

    if (Buffers::packed) {
        // projection_3dgut_packed_forward(
        //     Buffers::num_splats,
        //     Buffers::sh_degree,

        // )
    }
}

