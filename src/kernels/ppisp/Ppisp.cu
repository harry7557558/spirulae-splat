// Ppisp.cu -- PPISP (per-pixel image signal processing) forward, backward
// and regularization.
//
// Part of the PixelWise family -- see PixelWiseCommon.cuh.

#include "kernels/pixelwise/PixelWiseCommon.cuh"

// ================
// PPISP
// ================

static constexpr int kNumPPISPParams = 36;
static constexpr int kNumPPISPParamsRQS = 39;
static constexpr int kNumPPISPParamsNoCRF = 24;

enum class PPISPParamType : int {
    Original,
    RQS,
    NoCRF,
};

template<PPISPParamType param_type>
__global__ void ppisp_forward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [N_cam or B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    const int* __restrict__ cam_indices,  // [B], or nullptr -> identity
    TensorView<float, 4> out_image  // [B, H, W, C]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    int p_id = cam_indices ? cam_indices[bid] : (int)bid;

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[p_id * kNumParams + i];
    }

    float3 pixel = in_image.load3(bid, y, x);

    float3 out_pixel;
    if (param_type == PPISPParamType::Original)
        out_pixel = SlangPPISP::apply_ppisp(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params)
        );
    else if (param_type == PPISPParamType::RQS)
        out_pixel = SlangPPISP::apply_ppisp_rqs(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params)
        );
    else
        out_pixel = SlangPPISP::apply_ppisp_no_crf(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params)
        );

    out_image.store3(bid, y, x, out_pixel);
}

/*[AutoHeaderGeneratorExport]*/
void ppisp_forward(
    DeviceTensor3D<float3> in_image,    // [B, H, W, 3]
    TorchTensorView ppisp_params,       // [N_cam or B, PPISP_NUM_PARAMS]
    TorchTensorView intrins,            // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    std::string param_type,
    TorchTensorView cam_indices,        // [B] int32, or null -> identity (ppisp_params is [B,P])
    DeviceTensor3D<float3> out_image    // [B, H, W, 3]
) {
    long b = in_image.size<0>(), h = in_image.size<1>(), w = in_image.size<2>();
    const int* cam_idx_ptr = (std::get<0>(cam_indices) == 0) ?
        nullptr : (const int*)std::get<0>(cam_indices);
    if (param_type == "original" || param_type == "") {
        ppisp_forward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else if (param_type == "rqs") {
        ppisp_forward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else if (param_type == "no_crf") {
        ppisp_forward_kernel<PPISPParamType::NoCRF><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

template<PPISPParamType param_type>
__global__ void ppisp_backward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [N_cam or B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    const int* __restrict__ cam_indices,  // [B], or nullptr -> identity
    const TensorView<float, 4> v_out_image,  // [B, H, W, C]
    TensorView<float, 4> v_in_image,  // [B, H, W, C]
    float* __restrict__ v_ppisp_params  // [N_cam or B, PPISP_NUM_PARAMS]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    bool inside = (bid < B) && (gid < H*W);
    unsigned y = inside ? gid / W : 0u;
    unsigned x = inside ? gid % W : 0u;

    int p_id = (bid < B) ? (cam_indices ? cam_indices[bid] : (int)bid) : 0;

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
#if 0
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[p_id * kNumParams + i];
    }
#else
    __shared__ float params_shared[kNumParams];
    if (threadIdx.x < kNumParams) {  // assume blockDim.x >= kNumParams
        float value = ppisp_params[p_id * kNumParams + threadIdx.x];
        params_shared[threadIdx.x] = value;
    }
    __syncthreads();
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = params_shared[i];
        // int j = (i + threadIdx.x) % kNumParams;
        // params[j] = params_shared[j];
    }
#endif

    float3 v_pixel = make_float3(0.0f, 0.0f, 0.0f);
    FixedArray<float, kNumParams> v_params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++)
        v_params[i] = 0.0f;

    if (inside) {
        float3 pixel = in_image.load3(bid, y, x);
        float3 v_out_pixel = v_out_image.load3(bid, y, x);
        if (param_type == PPISPParamType::Original)
            SlangPPISP::apply_ppisp_vjp(
                pixel,
                make_float2((float)x, (float)y),
                make_float2(intrins[bid].z, intrins[bid].w),
                make_float2(actual_image_width, actual_image_height),
                *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
                v_out_pixel,
                &v_pixel,
                reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&v_params)
            );
        else if (param_type == PPISPParamType::RQS)
            SlangPPISP::apply_ppisp_rqs_vjp(
                pixel,
                make_float2((float)x, (float)y),
                make_float2(intrins[bid].z, intrins[bid].w),
                make_float2(actual_image_width, actual_image_height),
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
                v_out_pixel,
                &v_pixel,
                reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&v_params)
            );
        else
            SlangPPISP::apply_ppisp_no_crf_vjp(
                pixel,
                make_float2((float)x, (float)y),
                make_float2(intrins[bid].z, intrins[bid].w),
                make_float2(actual_image_width, actual_image_height),
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
                v_out_pixel,
                &v_pixel,
                reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&v_params)
            );

        v_in_image.store3(bid, y, x, v_pixel);
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        param = cg::reduce(warp, param, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && param != 0.0f)
            atomicAdd(&v_ppisp_params[p_id * kNumParams + i], param);
    }
}

/*[AutoHeaderGeneratorExport]*/
void ppisp_backward(
    DeviceTensor3D<float3> in_image,    // [B, H, W, 3]
    TorchTensorView ppisp_params,       // [N_cam or B, PPISP_NUM_PARAMS]
    TorchTensorView intrins,            // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    DeviceTensor3D<float3> v_out_image, // [B, H, W, 3]
    std::string param_type,
    TorchTensorView cam_indices,        // [B] int32, or null -> identity
    DeviceTensor3D<float3> v_in_image,  // [B, H, W, 3]
    TorchTensorView v_ppisp_params      // [N_cam or B, PPISP_NUM_PARAMS] (must be pre-zeroed)
) {
    long b = in_image.size<0>(), h = in_image.size<1>(), w = in_image.size<2>();
    const int* cam_idx_ptr = (std::get<0>(cam_indices) == 0) ?
        nullptr : (const int*)std::get<0>(cam_indices);
    if (param_type == "original" || param_type == "") {
        ppisp_backward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else if (param_type == "rqs") {
        ppisp_backward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else if (param_type == "no_crf") {
        ppisp_backward_kernel<PPISPParamType::NoCRF><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

template<PPISPParamType param_type>
__global__ void compute_raw_ppisp_regularization_forward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    float* __restrict__ raw_losses  // [B+1, RawPPISPRegLossIndex::length]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    bool inside = (bid < (unsigned)B);

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

    FixedArray<float, kNumRawLosses> losses;
    if (inside) {
        FixedArray<float, kNumParams> params;
        #pragma unroll
        for (int i = 0; i < kNumParams; i++) {
            params[i] = ppisp_params[bid * kNumParams + i];
        }

        if (param_type == PPISPParamType::Original)
            SlangPPISP::compute_raw_ppisp_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&losses)
            );
        else if (param_type == PPISPParamType::RQS)
            SlangPPISP::compute_raw_ppisp_rqs_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&losses)
            );
        else
            SlangPPISP::compute_raw_ppisp_no_crf_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&losses)
            );
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        float loss = (inside && isfinite(losses[i])) ? losses[i] : 0.0f;
        if (inside)
            raw_losses[bid*kNumRawLosses + i] = loss;
        loss = cg::reduce(warp, loss, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && loss != 0.0f)
            atomicAdd(&raw_losses[B*kNumRawLosses + i], loss);
    }
}

template<PPISPParamType param_type>
__global__ void compute_ppisp_regularization_forward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    float* __restrict__ losses_buffer  // [PPISPRegLossIndex::length]
) {
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

    FixedArray<float, kNumRawLosses> raw_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (int)PPISPRegLossIndex::length> losses;

    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_ppisp_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_ppisp_rqs_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );
    else
        SlangPPISP::compute_ppisp_no_crf_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );

    #pragma unroll
    for (int i = 0; i < (int)PPISPRegLossIndex::length; i++) {
        losses_buffer[i] = losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
void compute_ppsip_regularization_forward(
    TorchTensorView ppisp_params,       // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    std::string param_type,
    TorchTensorView losses,             // [PPISPRegLossIndex::length] (must be pre-zeroed)
    TorchTensorView raw_losses          // [B+1, RawPPISPRegLossIndex::length] (must be pre-zeroed)
) {
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = std::get<2>(ppisp_params)[0];

    if (param_type == "original" || param_type == "") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndex::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "no_crf") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::NoCRF>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::NoCRF>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndexNoCRF::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
}

template<PPISPParamType param_type>
__global__ void compute_raw_ppisp_regularization_backward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const float* __restrict__ v_raw_losses,  // [RawPPISPRegLossIndex::length]
    float* __restrict__ v_ppisp_params  // [B, PPISP_NUM_PARAMS]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= B)
        return;

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[bid * kNumParams + i];
    }

    FixedArray<float, kNumRawLosses> v_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        v_losses[i] = v_raw_losses[i];
    }

    FixedArray<float, kNumParams> v_params;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_raw_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&v_params)
        );
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_raw_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&v_params)
        );
    else
        SlangPPISP::compute_raw_ppisp_no_crf_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&v_params)
        );

    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        v_ppisp_params[bid * kNumParams + i] = param;
    }
}

template<PPISPParamType param_type>
__global__ void compute_ppisp_regularization_backward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    const float* __restrict__ v_losses_buffer,  // [PPISPRegLossIndex::length]
    float* __restrict__ v_raw_losses_buffer  // [RawPPISPRegLossIndex::length]
) {
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

    FixedArray<float, kNumRawLosses> raw_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (int)PPISPRegLossIndex::length> v_losses;
    #pragma unroll
    for (int i = 0; i < (int)PPISPRegLossIndex::length; i++) {
        v_losses[i] = v_losses_buffer[i];
    }

    FixedArray<float, kNumRawLosses> v_raw_losses;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_raw_losses)
        );
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_raw_losses)
        );
    else
        SlangPPISP::compute_ppisp_no_crf_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&v_raw_losses)
        );

    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        v_raw_losses_buffer[i] = v_raw_losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
void compute_ppsip_regularization_backward(
    TorchTensorView ppisp_params,       // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    TorchTensorView raw_losses,         // [B+1, RawPPISPRegLossIndex::length]
    TorchTensorView v_losses,           // [PPISPRegLossIndex::length]
    std::string param_type,
    TorchTensorView v_ppisp_params      // [B, PPISP_NUM_PARAMS] (must be pre-zeroed)
) {
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = std::get<2>(ppisp_params)[0];

    if (param_type == "original" || param_type == "") {
        // v_raw_losses is a small scratch buffer
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndex::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndex::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndex::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndexRQS::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndexRQS::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "no_crf") {
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndexNoCRF::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndexNoCRF::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::NoCRF>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndexNoCRF::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::NoCRF>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
}
