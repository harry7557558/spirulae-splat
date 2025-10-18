#pragma once

#ifdef __CUDACC__
#define TensorView _Slang_TensorView
#include "generated/projection_opaque_triangle.cu"
#undef TensorView
#endif

#include "common.cuh"

#include <tuple>


#ifdef __CUDACC__

/*
https://harry7557558.github.io/spirulae/autodiff/

Language:
C/C++ (float)

Code function:
```
__forceinline__ __device__ float fun(float2 p0, float2 p1, float2 p2, float2 p, float t) {
    float x_0 = p0.x, y_0 = p0.y, x_1 = p1.x, y_1 = p1.y, x_2 = p2.x, y_2 = p2.y, x = p.x, y = p.y;
    float {%funbody%};
    return {%a%};
}
```
OR:
```
__forceinline__ __device__ void fun_vjp(float2 p0, float2 p1, float2 p2, float2 p, float t, float g_out, float2& g_p0, float2& g_p1, float2& g_p2, float& g_t) {
    float x_0 = p0.x, y_0 = p0.y, x_1 = p1.x, y_1 = p1.y, x_2 = p2.x, y_2 = p2.y, x = p.x, y = p.y;
    float {%funbody%};
    g_p0 = g_out * make_float2({%a;x_0%}, {%a;y_0%});
    g_p1 = g_out * make_float2({%a;x_1%}, {%a;y_1%});
    g_p2 = g_out * make_float2({%a;x_2%}, {%a;y_2%});
    g_t = g_out * {%a;t%};
}
```

Independent variable(s):
x_0,y_0,x_1,y_1,x_2,y_2,x,y,t

Dependent variable(s):
a

Expression:
```
v0=vec2(x0,y0)
v1=vec2(x1,y1)
v2=vec2(x2,y2)
e0=v1-v0
e1=v2-v1
e2=v0-v2
s_at(x)=clamp(x,0,1)
s_e=sign(det(e0,e1))
d_e0=s_e*det(vec2(x,y)-v0,normalize(e0))
d_e1=s_e*det(vec2(x,y)-v1,normalize(e1))
d_e2=s_e*det(vec2(x,y)-v2,normalize(e2))
d_v=max(d_e0,d_e1,d_e2)
d_max=det(e0,e1)/(|e0|+|e1|+|e2|)
a=s_at(-d_v/(2d_max*(1-t)+t)+0.5)
#a=sigmoid(-4d_v/(2d_max*(1-t)+t))
```
OR:
```
v0=vec2(x0,y0)
v1=vec2(x1,y1)
v2=vec2(x2,y2)
e0=v1-v0
e1=v2-v1
e2=v0-v2
s_at(x)=clamp(x,0,1)
d_v0(x,y)=|(vec2(x,y)-v0)-e0*s_at(dot(vec2(x,y)-v0,e0)/dot(e0,e0))|
d_v1(x,y)=|(vec2(x,y)-v1)-e1*s_at(dot(vec2(x,y)-v1,e1)/dot(e1,e1))|
d_v2(x,y)=|(vec2(x,y)-v2)-e2*s_at(dot(vec2(x,y)-v2,e2)/dot(e2,e2))|
s_v(x,y)=sign(max(det(vec2(x,y)-v0,e0),det(vec2(x,y)-v1,e1),det(vec2(x,y)-v2,e2)
d_v=s_v(x,y)*min(d_v0(x,y),d_v1(x,y),d_v2(x,y))
d_max=det(e0,e1)/(|e0|+|e1|+|e2|)
a=s_at(-d_v/(2d_max*(1-t)+t)+0.5)
#a=sigmoid(-4d_v/(2d_max*(1-t)+t))
```
*/

namespace _PrimitiveOpaqueTriangle {

__forceinline__ __device__ float triangle_opac_approx(float2 p0, float2 p1, float2 p2, float2 p, float t) {
    float x_0 = p0.x, y_0 = p0.y, x_1 = p1.x, y_1 = p1.y, x_2 = p2.x, y_2 = p2.y, x = p.x, y = p.y;
    float v0=x_1-x_0, v1=y_1-y_0, v2=x_2-x_1, v3=y_2-y_1, v4=v0*v3, v5=v1*v2, v6=v4-v5, v7=v6>0.0f?1.0f:v6<0.0f?-1.0f:0.0f, v8=x-x_0, v9=y-y_0, v10=sqrt(v0*v0+v1*v1), v11=v0/v10, v12=v1/v10, v13=v8*v12, v14=v9*v11, v15=v13-v14, v16=v7*v15, v17=x-x_1, v18=y-y_1, v19=sqrt(v2*v2+v3*v3), v20=v2/v19, v21=v3/v19, v22=v17*v21, v23=v18*v20, v24=v22-v23, v25=v7*v24, v26=x-x_2, v27=y-y_2, v28=x_0-x_2, v29=y_0-y_2, v30=sqrt(v28*v28+v29*v29), v31=v28/v30, v32=v29/v30, v33=v26*v32, v34=v27*v31, v35=v33-v34, v36=v7*v35, v37=fmax(fmax(v16,v25),v36), v38=-v37, v39=v10+v19, v40=v39+v30, v41=v6/v40, v42=2.0f*v41, v43=1.0f-t, v44=v42*v43, v45=v44+t, v46=v38/v45, v47=v46+0.5f, v48=fmin(fmax(v47,0.0f),1.0f);
    return v48;
}

__forceinline__ __device__ float triangle_opac_approx_vjp(float2 p0, float2 p1, float2 p2, float2 p, float t, float g_out, float2& g_p0, float2& g_p1, float2& g_p2, float& g_t) {
    float x_0 = p0.x, y_0 = p0.y, x_1 = p1.x, y_1 = p1.y, x_2 = p2.x, y_2 = p2.y, x = p.x, y = p.y;
    float v0=x_1-x_0, v1=y_1-y_0, v2=x_2-x_1, v3=y_2-y_1, v4=v0*v3, v5=v1*v2, v6=v4-v5, v7=v6>0.0f?1.0f:v6<0.0f?-1.0f:0.0f, v8=x-x_0, v9=y-y_0, v10=sqrt(v0*v0+v1*v1), v11=v0/v10, v12=v1/v10, v13=v8*v12, v14=v9*v11, v15=v13-v14, v16=v7*v15, v17=x-x_1, v18=y-y_1, v19=sqrt(v2*v2+v3*v3), v20=v2/v19, v21=v3/v19, v22=v17*v21, v23=v18*v20, v24=v22-v23, v25=v7*v24, v26=x-x_2, v27=y-y_2, v28=x_0-x_2, v29=y_0-y_2, v30=sqrt(v28*v28+v29*v29), v31=v28/v30, v32=v29/v30, v33=v26*v32, v34=v27*v31, v35=v33-v34, v36=v7*v35, v37=fmax(fmax(v16,v25),v36), v38=-v37, v39=v10+v19, v40=v39+v30, v41=v6/v40, v42=2.0f*v41, v43=1.0f-t, v44=v42*v43, v45=v44+t, v46=v38/v45, v47=v46+0.5f, v48=fmin(fmax(v47,0.0f),1.0f), v49=(-1.0f)*v3, v50=(-2.0f)*v0, v51=v0*v0, v52=v1*v1, v53=v51+v52, v54=sqrt(v53), v55=2.0f*v54, v56=v50/v55, v57=v1*v56, v58=-v57, v59=v54*v54, v60=v58/v59, v61=v1/v54, v62=(-1.0f)*v61, v63=v8*v60, v64=v62+v63, v65=v8*v61, v66=(-1.0f)*v54, v67=v0*v56, v68=v66-v67, v69=v68/v59, v70=v0/v54, v71=v9*v69, v72=v9*v70, v73=v64-v71, v74=v65-v72, v75=v7*v73, v76=v7*v74, v77=v2*v2, v78=v3*v3, v79=v77+v78, v80=sqrt(v79), v81=2.0f*v80, v82=v80*v80, v83=v3/v80, v84=v17*v83, v85=v2/v80, v86=v18*v85, v87=v84-v86, v88=v7*v87, v89=v76-v88, v90=v89>0.0f?v76:v88, v91=v90-v88, v92=v91>0.0f?v75:0.0f, v93=fmax(v76,v88), v94=2.0f*v28, v95=v28*v28, v96=v29*v29, v97=v95+v96, v98=sqrt(v97), v99=2.0f*v98, v100=v94/v99, v101=v29*v100, v102=-v101, v103=v98*v98, v104=v102/v103, v105=v29/v98, v106=v26*v104, v107=v26*v105, v108=v28*v100, v109=v98-v108, v110=v109/v103, v111=v28/v98, v112=v27*v110, v113=v27*v111, v114=v106-v112, v115=v107-v113, v116=v7*v114, v117=v7*v115, v118=v93-v117, v119=v118>0.0f?v93:v117, v120=v119-v117, v121=v120>0.0f?v92:v116, v122=fmax(v93,v117), v123=-v121, v124=-v122, v125=v54+v80, v126=v56+v100, v127=v125+v98, v128=v49*v127, v129=v6*v126, v130=v128-v129, v131=v127*v127, v132=v130/v131, v133=v6/v127, v134=2.0f*v132, v135=2.0f*v133, v136=v134*v43, v137=v135*v43, v138=v137+t, v139=v123*v138, v140=v124*v136, v141=v139-v140, v142=v138*v138, v143=v141/v142, v144=v124/v138, v145=v144+0.5f, v146=v145-1.0f, v147=v145>0.0f?v143:0.0f, v148=v146>0.0f?0.0f:v147, v149=(-1.0f)*v2, v150=-v149, v151=(-2.0f)*v1, v152=v151/v55, v153=v1*v152, v154=v66-v153, v155=v154/v59, v156=v8*v155, v157=v0*v152, v158=-v157, v159=v158/v59, v160=(-1.0f)*v70, v161=v9*v159, v162=v160+v161, v163=v156-v162, v164=v7*v163, v165=v91>0.0f?v164:0.0f, v166=2.0f*v29, v167=v166/v99, v168=v29*v167, v169=v98-v168, v170=v169/v103, v171=v26*v170, v172=v28*v167, v173=-v172, v174=v173/v103, v175=v27*v174, v176=v171-v175, v177=v7*v176, v178=v120>0.0f?v165:v177, v179=-v178, v180=v152+v167, v181=v150*v127, v182=v6*v180, v183=v181-v182, v184=v183/v131, v185=2.0f*v184, v186=v185*v43, v187=v179*v138, v188=v124*v186, v189=v187-v188, v190=v189/v142, v191=v145>0.0f?v190:0.0f, v192=v146>0.0f?0.0f:v191, v193=v1*(-1.0f), v194=v3-v193, v195=2.0f*v0, v196=v195/v55, v197=v1*v196, v198=-v197, v199=v198/v59, v200=v8*v199, v201=v0*v196, v202=v54-v201, v203=v202/v59, v204=v9*v203, v205=v200-v204, v206=v7*v205, v207=(-2.0f)*v2, v208=v207/v81, v209=v3*v208, v210=-v209, v211=v210/v82, v212=(-1.0f)*v83, v213=v17*v211, v214=v212+v213, v215=(-1.0f)*v80, v216=v2*v208, v217=v215-v216, v218=v217/v82, v219=v18*v218, v220=v214-v219, v221=v7*v220, v222=v91>0.0f?v206:v221, v223=v120>0.0f?v222:0.0f, v224=-v223, v225=v196+v208, v226=v194*v127, v227=v6*v225, v228=v226-v227, v229=v228/v131, v230=2.0f*v229, v231=v230*v43, v232=v224*v138, v233=v124*v231, v234=v232-v233, v235=v234/v142, v236=v145>0.0f?v235:0.0f, v237=v146>0.0f?0.0f:v236, v238=v0*(-1.0f), v239=v238-v2, v240=2.0f*v1, v241=v240/v55, v242=v1*v241, v243=v54-v242, v244=v243/v59, v245=v8*v244, v246=v0*v241, v247=-v246, v248=v247/v59, v249=v9*v248, v250=v245-v249, v251=v7*v250, v252=(-2.0f)*v3, v253=v252/v81, v254=v3*v253, v255=v215-v254, v256=v255/v82, v257=v17*v256, v258=v2*v253, v259=-v258, v260=v259/v82, v261=(-1.0f)*v85, v262=v18*v260, v263=v261+v262, v264=v257-v263, v265=v7*v264, v266=v91>0.0f?v251:v265, v267=v120>0.0f?v266:0.0f, v268=-v267, v269=v241+v253, v270=v239*v127, v271=v6*v269, v272=v270-v271, v273=v272/v131, v274=2.0f*v273, v275=v274*v43, v276=v268*v138, v277=v124*v275, v278=v276-v277, v279=v278/v142, v280=v145>0.0f?v279:0.0f, v281=v146>0.0f?0.0f:v280, v282=-v1, v283=2.0f*v2, v284=v283/v81, v285=v3*v284, v286=-v285, v287=v286/v82, v288=v17*v287, v289=v2*v284, v290=v80-v289, v291=v290/v82, v292=v18*v291, v293=v288-v292, v294=v7*v293, v295=v91>0.0f?0.0f:v294, v296=(-2.0f)*v28, v297=v296/v99, v298=v29*v297, v299=-v298, v300=v299/v103, v301=(-1.0f)*v105, v302=v26*v300, v303=v301+v302, v304=(-1.0f)*v98, v305=v28*v297, v306=v304-v305, v307=v306/v103, v308=v27*v307, v309=v303-v308, v310=v7*v309, v311=v120>0.0f?v295:v310, v312=-v311, v313=v284+v297, v314=v282*v127, v315=v6*v313, v316=v314-v315, v317=v316/v131, v318=2.0f*v317, v319=v318*v43, v320=v312*v138, v321=v124*v319, v322=v320-v321, v323=v322/v142, v324=v145>0.0f?v323:0.0f, v325=v146>0.0f?0.0f:v324, v326=2.0f*v3, v327=v326/v81, v328=v3*v327, v329=v80-v328, v330=v329/v82, v331=v17*v330, v332=v2*v327, v333=-v332, v334=v333/v82, v335=v18*v334, v336=v331-v335, v337=v7*v336, v338=v91>0.0f?0.0f:v337, v339=(-2.0f)*v29, v340=v339/v99, v341=v29*v340, v342=v304-v341, v343=v342/v103, v344=v26*v343, v345=v28*v340, v346=-v345, v347=v346/v103, v348=(-1.0f)*v111, v349=v27*v347, v350=v348+v349, v351=v344-v350, v352=v7*v351, v353=v120>0.0f?v338:v352, v354=-v353, v355=v327+v340, v356=v0*v127, v357=v6*v355, v358=v356-v357, v359=v358/v131, v360=2.0f*v359, v361=v360*v43, v362=v354*v138, v363=v124*v361, v364=v362-v363, v365=v364/v142, v366=v145>0.0f?v365:0.0f, v367=v146>0.0f?0.0f:v366, v368=v135*(-1.0f), v369=v368+1.0f, v370=v124*v369, v371=-v370, v372=v371/v142, v373=v145>0.0f?v372:0.0f, v374=v146>0.0f?0.0f:v373;
    g_p0 = g_out * make_float2(v148, v192);
    g_p1 = g_out * make_float2(v237, v281);
    g_p2 = g_out * make_float2(v325, v367);
    g_t = g_out * v374;
}

}

#endif


struct OpaqueTriangle {
    struct World;
    struct Screen;

#ifdef __CUDACC__

    struct FwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
        float near_plane, far_plane, radius_clip;
    };

    inline static __device__ void project_persp(
        World world, FwdProjCamera cam,
        Screen& screen, int2& radius, float& depth
    );

    struct BwdProjCamera {
        float3x3 R;
        float3 t;
        float fx, fy, cx, cy;
        uint width, height, antialiased;
    };

    inline static __device__ void project_persp_vjp(
        World world, BwdProjCamera cam,
        Screen v_screen, float v_depth,
        World& v_world, float3x3 &v_R, float3 &v_t
    );

#endif  // #ifdef __CUDACC__

};

struct OpaqueTriangle::World {

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor vertices;  // [..., 3, 3]
        at::Tensor hardness;

        Tensor(const TensorTuple& splats) {
            vertices = std::get<0>(splats);
            hardness = std::get<1>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(vertices, hardness);
        }

        Tensor zeros_like() const {
            return Tensor(std::make_tuple(
                at::zeros_like(vertices),
                hardness
            ));
        }

        auto options() {
            return vertices.options();
        }
        long size() const {
            return vertices.size(-3);
        }
        long batchSize() const {
            return vertices.numel() / (9*size());
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float3* __restrict__ vertices;
        float* __restrict__ hardness;

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.vertices);
            CHECK_INPUT(tensors.vertices);
            CHECK_INPUT(tensors.hardness);
            vertices = (float3*)tensors.vertices.data_ptr<float>();
            hardness = tensors.vertices.data_ptr<float>();
        }
    };

    float3 vert0;
    float3 vert1;
    float3 vert2;
    float hardness;

#ifdef __CUDACC__

    static __device__ World load(const Buffer &buffer, long idx) {
        return {
            buffer.vertices[3*idx+0],
            buffer.vertices[3*idx+1],
            buffer.vertices[3*idx+2],
            buffer.hardness[idx]
        };
    }

    static __device__ __forceinline__ World zero() {
        return {
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            {0.f, 0.f, 0.f},
            0.f
        };
    }

    template<typename Partition>
    __device__ void reduce(Partition& partition) {
        warpSum(vert0, partition);
        warpSum(vert1, partition);
        warpSum(vert2, partition);
        warpSum(hardness, partition);
    }
    
    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.vertices + 3*idx+0, vert0);
        atomicAddFVec(buffer.vertices + 3*idx+1, vert1);
        atomicAddFVec(buffer.vertices + 3*idx+2, vert2);
        atomicAdd(buffer.hardness + idx, hardness);
    }

#endif  // #ifdef __CUDACC__
};

struct OpaqueTriangle::Screen {

    typedef std::tuple<at::Tensor, at::Tensor> TensorTuple;

    struct Buffer;

    struct Tensor {
        at::Tensor vertices;  // [..., 3, 2]
        at::Tensor hardness;  // [..., 3, 2]
        std::optional<at::Tensor> absgrad;

        Tensor(const TensorTuple& splats) {
            vertices = std::get<0>(splats);
            hardness = std::get<1>(splats);
        }

        TensorTuple tuple() const {
            return std::make_tuple(vertices, hardness);
        }

        Tensor zeros_like(bool absgrad) const {
            Tensor result = Tensor(std::make_tuple(
                at::zeros_like(vertices),
                at::zeros_like(hardness)
            ));
            if (absgrad) {
                // TODO: batched tensor
                result.absgrad = at::zeros({size(), 2}, vertices.options());
            }
            return result;
        }

        auto options() {
            return vertices.options();
        }
        bool isPacked() const {
            return vertices.dim() == 3;
        }
        long size() const {
            return vertices.size(-3);
        }

        Buffer buffer() { return Buffer(*this); }
    };

    struct Buffer {
        float2* __restrict__ vertices;  // [I, N, 3, 2] or [nnz, 3, 2]
        float* __restrict__ hardness;  // [I, N] or [nnz]
        float2* __restrict__ absgrad;  // [I, N, 2] or [nnz, 2]

        Buffer(const Tensor& tensors) {
            DEVICE_GUARD(tensors.vertices);
            CHECK_INPUT(tensors.vertices);
            CHECK_INPUT(tensors.hardness);
            vertices = (float2*)tensors.vertices.data_ptr<float>();
            hardness = tensors.hardness.data_ptr<float>();
            absgrad = tensors.absgrad.has_value() ?
                (float2*)tensors.absgrad.value().data_ptr<float>()
                : nullptr;
        }
    };

    float2 vert0;
    float2 vert1;
    float2 vert2;
    float hardness;
    float2 absgrad;

#ifdef __CUDACC__

    static __device__ Screen load(const Buffer &buffer, long idx) {
        return {
            buffer.vertices[3*idx+0],
            buffer.vertices[3*idx+1],
            buffer.vertices[3*idx+2],
            buffer.hardness[idx]
            // absgrad is undefined
        };
    }

    static __device__ __forceinline__ Screen zero() {
        return {
            {0.f, 0.f},
            {0.f, 0.f},
            {0.f, 0.f},
            0.0f,
            {0.f, 0.f}
        };
    }

    __device__ __forceinline__ void operator+=(const Screen &other) {
        vert0 += other.vert0;
        vert1 += other.vert1;
        vert2 += other.vert2;
        hardness += other.hardness;
        absgrad += fabs(other.vert0 + other.vert1 + other.vert2) / 3.0f;
    }

    __device__ void atomicAddBuffer(Buffer &buffer, long idx) {
        atomicAddFVec(buffer.vertices + 3*idx+0, vert0);
        atomicAddFVec(buffer.vertices + 3*idx+1, vert1);
        atomicAddFVec(buffer.vertices + 3*idx+2, vert2);
        atomicAdd(buffer.hardness + idx, hardness);
        if (buffer.absgrad != nullptr)
            atomicAddFVec(buffer.absgrad + idx, absgrad);
    }

    __device__ __forceinline__ float evaluate_alpha(float px, float py) {
        return _PrimitiveOpaqueTriangle::triangle_opac_approx(
            vert0, vert1, vert2, {px, py}, hardness
        );
    }

    __device__ __forceinline__ Screen evaluate_alpha_vjp(float px, float py, float v_alpha) {
        Screen v_splat = Screen::zero();
        _PrimitiveOpaqueTriangle::triangle_opac_approx_vjp(
            vert0, vert1, vert2, {px, py}, hardness, v_alpha,
            v_splat.vert0, v_splat.vert1, v_splat.vert2, v_splat.hardness
        );
        return v_splat;
    }

#endif  // #ifdef __CUDACC__

};


#ifdef __CUDACC__

inline __device__ void OpaqueTriangle::project_persp(
    OpaqueTriangle::World world, OpaqueTriangle::FwdProjCamera cam,
    OpaqueTriangle::Screen& screen, int2& radius, float& depth
) {
    projection_opaque_triangle_persp(
        world.vert0, world.vert1, world.vert2, world.hardness,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.width, cam.height, cam.near_plane, cam.far_plane, cam.radius_clip,
        &radius, &depth, &screen.vert0, &screen.vert1, &screen.vert2, &screen.hardness
    );
}

inline __device__ void OpaqueTriangle::project_persp_vjp(
    OpaqueTriangle::World world, OpaqueTriangle::BwdProjCamera cam,
    OpaqueTriangle::Screen v_screen, float v_depth,
    OpaqueTriangle::World& v_world, float3x3 &v_R, float3 &v_t
) {
    projection_opaque_triangle_persp_vjp(
        world.vert0, world.vert1, world.vert2, world.hardness,
        cam.R, cam.t, cam.fx, cam.fy, cam.cx, cam.cy,
        cam.width, cam.height,
        v_depth, v_screen.vert0, v_screen.vert1, v_screen.vert2, v_screen.hardness,
        &v_world.vert0, &v_world.vert1, &v_world.vert2, &v_world.hardness,
        &v_R, &v_t
    );
}

#endif  // #ifdef __CUDACC__
