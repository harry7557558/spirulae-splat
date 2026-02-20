#include "common.cuh"

#define WARP_SIZE 32

template<uint32_t nx, uint32_t ny, uint32_t nz>
__global__ void dct3d_type1_ortho_kernel(
    float* __restrict__ data,
    uint32_t B,
    uint32_t actual_nx, uint32_t actual_ny, uint32_t actual_nz,
    uint32_t actual_sx, uint32_t actual_sy, uint32_t actual_sz
) {
    static_assert(nz <= WARP_SIZE && (nz & (nz - 1)) == 0);
    static_assert((nx * ny * nz) % WARP_SIZE == 0);

    uint32_t bid = blockIdx.x * WARP_SIZE + threadIdx.x;
    uint32_t zid = bid % nz; bid /= nz;
    uint32_t yid = bid % ny; bid /= ny;
    uint32_t xid = bid % nx; bid /= nx;
    bool inside = (xid < actual_nx && yid < actual_ny && zid < actual_nz);

    uint32_t offset = bid * actual_nx * actual_ny * actual_nz +
        xid * actual_sx + yid * actual_sy + zid * actual_sz;
    float s = (zid == 0 || zid == actual_nz-1) ? sqrtf(0.5f) : 1.0f;
    float x = inside ? s * data[offset] : 0.0f;

    float inv_nz_1 = 1.0f / (actual_nz - 1);

    // performance is almost entirely memory bound, brute force method suffice
    float y = 0.0f;
    #pragma unroll
    for (int n = 0; n < nz; ++n) {
        y += __shfl_sync(~0u, x, n, nz) * cosf(M_PIf * zid * n * inv_nz_1);
    }

    s = (zid == 0 || zid == actual_nz-1) ? sqrtf(inv_nz_1) : sqrtf(2.0f * inv_nz_1);
    if (inside) data[offset] = y * s;
}


template<uint32_t nx, uint32_t ny, uint32_t nz>
void launch_dct3d_type1_ortho_kernel(
    at::Tensor data,
    uint32_t actual_nx, uint32_t actual_ny, uint32_t actual_nz
) {
    uint32_t B = data.numel() / (actual_nx * actual_ny * actual_nz);

    if (nz > 1) {
        dct3d_type1_ortho_kernel<nx, ny, nz>
        <<<B * nx * ny / (WARP_SIZE / nz), WARP_SIZE>>>(
            data.data_ptr<float>(), B, actual_nx, actual_ny, actual_nz,
                actual_ny * actual_nz, actual_nz, 1u);
    }

    if (ny > 1) {
        dct3d_type1_ortho_kernel<nx, nz, ny>
        <<<B * nx * nz / (WARP_SIZE / ny), WARP_SIZE>>>(
            data.data_ptr<float>(), B, actual_nx, actual_nz, actual_ny,
                actual_ny * actual_nz, 1u, actual_nz);
    }

    if (nx > 1) {
        dct3d_type1_ortho_kernel<ny, nz, nx>
        <<<B * ny * nz / (WARP_SIZE / nx), WARP_SIZE>>>(
            data.data_ptr<float>(), B, actual_ny, actual_nz, actual_nx,
                actual_nz, 1u, actual_ny * actual_nz);
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor dct3d_type1_ortho_tensor(at::Tensor input) {
    DEVICE_GUARD(input);
    CHECK_INPUT(input);

    uint32_t actual_nx = input.size(-3);
    uint32_t actual_ny = input.size(-2);
    uint32_t actual_nz = input.size(-1);
    if (actual_nx > WARP_SIZE || actual_ny > WARP_SIZE || actual_nz > WARP_SIZE)
        throw std::runtime_error("Grid dimension must be <=" + std::to_string(WARP_SIZE));

    uint32_t nx = actual_nx, ny = actual_ny, nz = actual_nz;
    while (nx & (nx - 1)) ++nx;
    while (ny & (ny - 1)) ++ny;
    while (nz & (nz - 1)) ++nz;

    auto output = input.clone();

    // template instantiate selected kernels to balance runtime speed and compilation time
    #define _LAUNCH(cx, cy, cz) \
        if (nx <= cx && ny <= cy && nz <= cz) \
            launch_dct3d_type1_ortho_kernel<cx, cy, cz>(output, actual_nx, actual_ny, actual_nz)
    _LAUNCH(2, 4, 4);
    else _LAUNCH(4, 4, 4);
    else _LAUNCH(4, 8, 8);
    else _LAUNCH(8, 8, 8);
    else _LAUNCH(8, 16, 16);
    else _LAUNCH(16, 16, 16);
    else _LAUNCH(16, 32, 32);
    else _LAUNCH(32, 32, 32);
    #undef _LAUNCH

    return output;
}
