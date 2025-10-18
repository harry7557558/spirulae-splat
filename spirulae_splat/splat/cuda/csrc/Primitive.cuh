#pragma once


#include "common.cuh"

#include <tuple>


#ifdef __CUDACC__

inline __device__ float triangle_opac_approx(float2 p0, float2 p1, float2 p2, float2 p, float t) {
    float x_0 = p0.x, y_0 = p0.y, x_1 = p1.x, y_1 = p1.y, x_2 = p2.x, y_2 = p2.y, x = p.x, y = p.y;
    float v0=x_1-x_0, v1=y_1-y_0, v2=x_2-x_1, v3=y_2-y_1, v4=v0*v3, v5=v1*v2, v6=v4-v5, v7=v6>0.0f?1.0f:v6<0.0f?-1.0f:0.0f, v8=x-x_0, v9=y-y_0, v10=sqrt(v0*v0+v1*v1), v11=v0/v10, v12=v1/v10, v13=v8*v12, v14=v9*v11, v15=v13-v14, v16=v7*v15, v17=x-x_1, v18=y-y_1, v19=sqrt(v2*v2+v3*v3), v20=v2/v19, v21=v3/v19, v22=v17*v21, v23=v18*v20, v24=v22-v23, v25=v7*v24, v26=x-x_2, v27=y-y_2, v28=x_0-x_2, v29=y_0-y_2, v30=sqrt(v28*v28+v29*v29), v31=v28/v30, v32=v29/v30, v33=v26*v32, v34=v27*v31, v35=v33-v34, v36=v7*v35, v37=fmax(fmax(v16,v25),v36), v38=-v37, v39=v10+v19, v40=v39+v30, v41=v6/v40, v42=2.0f*v41, v43=1.0f-t, v44=v42*v43, v45=v44+t, v46=v38/v45, v47=v46+0.5f, v48=fmin(fmax(v47,0.0f),1.0f);
    return v48;
}

#endif


struct Vanilla3DGS {
    struct World;
    struct Screen;
};

struct Vanilla3DGS::World {

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means;
        at::Tensor quats;
        at::Tensor scales;
        at::Tensor opacities;

        Tensor(const TensorTuple& splats) {
            means = std::get<0>(splats);
            quats = std::get<1>(splats);
            scales = std::get<2>(splats);
            opacities = std::get<3>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means, quats, scales, opacities);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(means),
                at::zeros_like(quats),
                at::zeros_like(scales),
                at::zeros_like(opacities)
            ));
        }

        auto options() {
            return means.options();
        }
        long size() const {
            return quats.size(-2);
        }
        long batchSize() const {
            return quats.numel() / (4*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ means;
        float4* __restrict__ quats;
        float3* __restrict__ scales;
        float* __restrict__ opacities;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means);
            CHECK_INPUT(tensors.means);
            CHECK_INPUT(tensors.quats);
            CHECK_INPUT(tensors.scales);
            CHECK_INPUT(tensors.opacities);
            means = (float3*)tensors.means.data_ptr<float>();
            quats = (float4*)tensors.quats.data_ptr<float>();
            scales = (float3*)tensors.scales.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
        }
    };

    float3 mean;
    float4 quat;
    float3 scale;
    float opacity;

#ifdef __CUDACC__

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(mean, partition);
        warpSum(quat, partition);
        warpSum(scale, partition);
        warpSum(opacity, partition);
    }
    
    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.means + idx, mean);
        atomicAddFVec(buffer.quats + idx, quat);
        atomicAddFVec(buffer.scales + idx, scale);
        atomicAdd(buffer.opacities + idx, opacity);
    }

#endif  // #ifdef __CUDACC__
};

struct Vanilla3DGS::Screen {

    typedef std::tuple<at::Tensor, at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor means2d;
        at::Tensor conics;
        at::Tensor opacities;
        std::optional<at::Tensor> absgrad;

        Tensor(const TensorTuple& splats) {
            means2d = std::get<0>(splats);
            conics = std::get<1>(splats);
            opacities = std::get<2>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(means2d, conics, opacities);
        }

        Tensor zeros_like(bool absgrad) const {
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(means2d),
                at::zeros_like(conics),
                at::zeros_like(opacities)
            ));
            if (absgrad)
                result.absgrad = at::zeros_like(means2d);
            return result;
        }

        auto options() {
            return means2d.options();
        }
        bool isPacked() const {
            return means2d.dim() == 2;
        }
        long size() const {
            return means2d.size(-2);
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float2* __restrict__ means2d;  // [I, N, 2] or [nnz, 2]
        float3* __restrict__ conics;  // [I, N, 3] or [nnz, 3]
        float* __restrict__ opacities;  // [I, N] or [nnz]
        float2* __restrict__ absgrad;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.means2d);
            CHECK_INPUT(tensors.means2d);
            CHECK_INPUT(tensors.conics);
            CHECK_INPUT(tensors.opacities);
            means2d = (float2*)tensors.means2d.data_ptr<float>();
            conics = (float3*)tensors.conics.data_ptr<float>();
            opacities = tensors.opacities.data_ptr<float>();
            absgrad = tensors.absgrad.has_value() ?
                (float2*)tensors.absgrad.value().data_ptr<float>()
                : nullptr;
        }
    };

    float2 xy;
    float3 conic;
    float opac;
    float2 xy_abs;

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.means2d[idx],
            buffer.conics[idx],
            buffer.opacities[idx],
            // xy_abs is undefined
        };
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.0f,
            {0.f, 0.f}
        };
    }

    __device__ __forceinline__ void operator+=(const Screen &other) {
        xy += other.xy;
        conic += other.conic;
        opac += other.opac;
        xy_abs += fabs(other.xy);
    }

    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        float *v_conic_ptr = (float*)buffer.conics + 3 * idx;
        if (conic.x != 0.0f) atomicAdd(v_conic_ptr, conic.x);
        if (conic.y != 0.0f) atomicAdd(v_conic_ptr + 1, conic.y);
        if (conic.z != 0.0f) atomicAdd(v_conic_ptr + 2, conic.z);

        float *v_xy_ptr = (float*)buffer.means2d + 2 * idx;
        if (xy.x != 0.0f) atomicAdd(v_xy_ptr, xy.x);
        if (xy.y != 0.0f) atomicAdd(v_xy_ptr + 1, xy.y);

        if (buffer.absgrad != nullptr) {
            float *v_xy_abs_ptr = (float*)buffer.absgrad + 2 * idx;
            if (xy_abs.x != 0.0f) atomicAdd(v_xy_abs_ptr, xy_abs.x);
            if (xy_abs.y != 0.0f) atomicAdd(v_xy_abs_ptr + 1, xy_abs.y);
        }

        if (opac != 0.0f) atomicAdd((float*)buffer.opacities + idx, opac);
    }

    __device__ __forceinline__ float evaluate_alpha(float px, float py) {
        // float size = 0.2f * rsqrtf(conic.x*conic.z-0.25f*conic.y*conic.y);
        // float2 p0 = {xy.x + size, xy.y - size};
        // float2 p1 = {xy.x - size, xy.y - size};
        // float2 p2 = {xy.x, xy.y + size};
        // return 0.01f+0.98f*triangle_opac_approx(p0, p1, p2, {px,py}, 0.98f);
        float2 delta = {xy.x - px, xy.y - py};
        float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                conic.z * delta.y * delta.y) +
                        conic.y * delta.x * delta.y;
        float vis = __expf(-sigma);
        float alpha = min(0.999f, opac * vis);
        return sigma < 0.f ? 0.f : alpha;
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(float px, float py, float v_alpha) {
        float2 delta = {xy.x - px, xy.y - py};
        float sigma = 0.5f * (conic.x * delta.x * delta.x +
                                conic.z * delta.y * delta.y) +
                        conic.y * delta.x * delta.y;
        float vis = __expf(-sigma);
        float alpha = min(0.999f, opac * vis);

        Screen v_splat = Screen::zero();
        if (sigma >= 0.f && opac * vis <= 0.999f) {
            const float v_sigma = -opac * vis * v_alpha;
            v_splat.conic = {
                0.5f * v_sigma * delta.x * delta.x,
                v_sigma * delta.x * delta.y,
                0.5f * v_sigma * delta.y * delta.y
            };
            v_splat.xy = {
                v_sigma * (conic.x * delta.x + conic.y * delta.y),
                v_sigma * (conic.y * delta.x + conic.z * delta.y)
            };
            v_splat.opac = vis * v_alpha;
        }
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};
