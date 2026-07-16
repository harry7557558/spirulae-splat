// Standalone WebGL viewer — performance-critical parsing / sorting / statistics.
//
// Compiled to WebAssembly (see CMakeLists.txt). Exposes a small C ABI consumed
// by js/wasm.js. Handles:
//   - 3DGS PLY parsing (binary LE + ascii), arbitrary property order, SH rest.
//   - Mesh PLY parsing (binary LE + ascii), triangulating polygon faces.
//   - Wavefront OBJ parsing (v/vt/vn, f with shared or separate indices).
//   - Depth sorting of splats (16-bit counting sort, back-to-front).
//   - Parameter histograms (opacities / scales / erank / rgb; mesh edge/area).
// GLTF/GLB are parsed in JavaScript (JSON + typed-array buffer views map
// directly), so they are intentionally absent here.
//
// This file is deliberately dependency-free (only the C++ standard library).

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <emscripten.h>

// Fast decimal parsing (FastFloat.h lives with the shared dataset parsers;
// the viewer CMake adds that directory to the include path). musl's
// strtod/strtof cost ~1.5us per call — tens of seconds over a 100+ MB OBJ
// or ascii PLY.
#include "FastFloat.h"

#define KEEP EMSCRIPTEN_KEEPALIVE extern "C"

// ---------------------------------------------------------------------------
// Result storage (module-global; one model at a time).
// ---------------------------------------------------------------------------

static const float SH_C0 = 0.28209479177387814f;

// ---- IEEE 754 half-float conversion. Everything except positions/opacity is
// stored as f16: it halves the memory of multi-GB models (a 20M-splat SH3
// model would not fit the 4GB wasm32 heap in f32) and uploads directly into
// RGBA16F/RGB16F textures without a driver-side conversion pass. ----
static inline uint16_t f2h(float f) {
    uint32_t x; memcpy(&x, &f, 4);
    uint32_t sign = (x >> 16) & 0x8000u;
    int32_t  exp  = (int32_t)((x >> 23) & 0xff) - 127 + 15;
    uint32_t man  = x & 0x7fffffu;
    if (exp <= 0) {                       // subnormal / underflow
        if (exp < -10) return (uint16_t)sign;
        man |= 0x800000u;
        uint32_t shift = (uint32_t)(14 - exp);
        uint16_t v = (uint16_t)(man >> shift);
        if ((man >> (shift - 1)) & 1u) v++;               // round to nearest
        return (uint16_t)(sign | v);
    }
    if (exp >= 31) return (uint16_t)(sign | 0x7c00u);     // inf / overflow
    uint16_t v = (uint16_t)(sign | ((uint32_t)exp << 10) | (man >> 13));
    if (man & 0x1000u) v++;               // round (carry into exponent is fine)
    return v;
}
static inline float h2f(uint16_t h) {
    uint32_t sign = ((uint32_t)h & 0x8000u) << 16;
    uint32_t exp  = (h >> 10) & 0x1fu;
    uint32_t man  = h & 0x3ffu;
    uint32_t x;
    if (exp == 0) {
        if (!man) x = sign;
        else {                             // subnormal
            exp = 127 - 15 + 1;
            while (!(man & 0x400u)) { man <<= 1; exp--; }
            man &= 0x3ffu;
            x = sign | (exp << 23) | (man << 13);
        }
    } else if (exp == 31) x = sign | 0x7f800000u | (man << 13);
    else x = sign | ((exp - 15 + 127) << 23) | (man << 13);
    float f; memcpy(&f, &x, 4); return f;
}

struct SplatData {
    uint32_t count = 0;
    int sh_degree = 0;            // 0..4
    int sh_rest = 0;              // rest coeffs per channel K (0,3,8,15,24)
    std::vector<float> posop;     // count*4 : x,y,z, opacity(sigmoid) — f32
    std::vector<uint16_t> quat;   // count*4 : x,y,z,w (normalized) — f16
    std::vector<uint16_t> scale;  // count*3 : exp(log_scale) world — f16
    std::vector<uint16_t> fdc;    // count*3 : raw f_dc — f16
    std::vector<uint16_t> sh;     // count*sh_rest*3 : [coeff][channel] — f16
    void clear() {
        count = 0; sh_degree = 0; sh_rest = 0;
        posop.clear(); quat.clear(); scale.clear(); fdc.clear(); sh.clear();
        posop.shrink_to_fit(); quat.shrink_to_fit(); scale.shrink_to_fit();
        fdc.shrink_to_fit(); sh.shrink_to_fit();
    }
};

struct MeshData {
    uint32_t nv = 0, nt = 0;
    std::vector<float> pos;      // nv*3
    std::vector<float> normal;   // nv*3 or empty
    std::vector<uint8_t> color;  // nv*3 or empty
    std::vector<float> uv;       // nv*2 or empty
    std::vector<uint32_t> idx;   // nt*3
    void clear() {
        nv = nt = 0;
        pos.clear(); normal.clear(); color.clear(); uv.clear(); idx.clear();
        pos.shrink_to_fit(); normal.shrink_to_fit(); color.shrink_to_fit();
        uv.shrink_to_fit(); idx.shrink_to_fit();
    }
};

static SplatData g_splat;
static MeshData  g_mesh;
static std::vector<uint32_t> g_order;   // sorted splat indices
static std::vector<float> g_depths;
static std::vector<float> g_bins;       // histogram output
static float g_minmax[2];
static int g_last_kind = 0;             // 0 none, 1 splat, 2 mesh

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

static inline float sigmoidf(float x) { return 1.0f / (1.0f + std::exp(-x)); }

// ---------------------------------------------------------------------------
// PLY parsing (shared header parser for splat + mesh)
// ---------------------------------------------------------------------------

enum PlyType { PT_I8, PT_U8, PT_I16, PT_U16, PT_I32, PT_U32, PT_F32, PT_F64, PT_NONE };

static int ply_type_size(int t) {
    switch (t) {
        case PT_I8: case PT_U8: return 1;
        case PT_I16: case PT_U16: return 2;
        case PT_I32: case PT_U32: case PT_F32: return 4;
        case PT_F64: return 8;
        default: return 0;
    }
}
static int ply_type_of(const std::string& s) {
    if (s=="char"||s=="int8") return PT_I8;
    if (s=="uchar"||s=="uint8") return PT_U8;
    if (s=="short"||s=="int16") return PT_I16;
    if (s=="ushort"||s=="uint16") return PT_U16;
    if (s=="int"||s=="int32") return PT_I32;
    if (s=="uint"||s=="uint32") return PT_U32;
    if (s=="float"||s=="float32") return PT_F32;
    if (s=="double"||s=="float64") return PT_F64;
    return PT_NONE;
}

struct PlyProp {
    std::string name;
    int type = PT_NONE;
    bool is_list = false;
    int count_type = PT_NONE;   // for lists
};
struct PlyElem {
    std::string name;
    uint64_t count = 0;
    std::vector<PlyProp> props;
};

struct PlyHeader {
    int format = 0;             // 0 ascii, 1 binary_le, 2 binary_be
    std::vector<PlyElem> elems;
    size_t data_offset = 0;
};

// read a little-endian scalar of given type as double
static inline double read_scalar_le(const uint8_t* p, int type) {
    switch (type) {
        case PT_I8:  return (double)(int8_t)p[0];
        case PT_U8:  return (double)p[0];
        case PT_I16: { int16_t v; memcpy(&v,p,2); return v; }
        case PT_U16: { uint16_t v; memcpy(&v,p,2); return v; }
        case PT_I32: { int32_t v; memcpy(&v,p,4); return v; }
        case PT_U32: { uint32_t v; memcpy(&v,p,4); return (double)v; }
        case PT_F32: { float v; memcpy(&v,p,4); return v; }
        case PT_F64: { double v; memcpy(&v,p,8); return v; }
    }
    return 0.0;
}

static bool parse_ply_header(const uint8_t* data, size_t len, PlyHeader& h) {
    // find "end_header\n"
    const char* marker = "end_header";
    size_t i = 0;
    // header is ascii text
    std::string text;
    size_t scan_limit = len < (1u<<20) ? len : (1u<<20);
    for (i = 0; i < scan_limit; i++) {
        text.push_back((char)data[i]);
        if (i >= 10 && text.size() >= 10 &&
            text.compare(text.size()-10, 10, marker) == 0) {
            // consume trailing newline
            size_t j = i + 1;
            if (j < len && data[j] == '\r') j++;
            if (j < len && data[j] == '\n') j++;
            h.data_offset = j;
            break;
        }
    }
    if (h.data_offset == 0) return false;

    // tokenize lines
    size_t pos = 0;
    auto next_line = [&](std::string& out)->bool{
        if (pos >= text.size()) return false;
        size_t e = text.find('\n', pos);
        if (e == std::string::npos) e = text.size();
        out = text.substr(pos, e - pos);
        if (!out.empty() && out.back()=='\r') out.pop_back();
        pos = e + 1;
        return true;
    };
    std::string line;
    if (!next_line(line)) return false; // "ply"
    while (next_line(line)) {
        if (line.rfind("format", 0) == 0) {
            if (line.find("ascii") != std::string::npos) h.format = 0;
            else if (line.find("binary_little_endian") != std::string::npos) h.format = 1;
            else if (line.find("binary_big_endian") != std::string::npos) h.format = 2;
        } else if (line.rfind("element", 0) == 0) {
            PlyElem e;
            char name[128]; unsigned long long cnt = 0;
            if (sscanf(line.c_str(), "element %127s %llu", name, &cnt) == 2) {
                e.name = name; e.count = cnt;
                h.elems.push_back(e);
            }
        } else if (line.rfind("property", 0) == 0) {
            if (h.elems.empty()) continue;
            PlyProp pr;
            char a[64], b[64], c[64];
            if (sscanf(line.c_str(), "property list %63s %63s %63s", a, b, c) == 3) {
                pr.is_list = true;
                pr.count_type = ply_type_of(a);
                pr.type = ply_type_of(b);
                pr.name = c;
            } else if (sscanf(line.c_str(), "property %63s %63s", a, b) == 2) {
                pr.type = ply_type_of(a);
                pr.name = b;
            } else continue;
            h.elems.back().props.push_back(pr);
        } else if (line.rfind("end_header", 0) == 0) {
            break;
        }
    }
    return !h.elems.empty();
}

// ---------- 3DGS splat PLY ----------

// Per-record parsing context for binary-LE splat vertices, shared by the
// whole-buffer and streaming paths. splat_ctx_init also (re)allocates the
// destination arrays (idempotent if already sized).
struct SplatRecCtx {
    int stride = 0;
    std::vector<int> off, type;                 // per-prop byte offset / type
    int ix=-1,iy=-1,iz=-1,io=-1;
    int is0=-1,is1=-1,is2=-1;
    int ir0=-1,ir1=-1,ir2=-1,ir3=-1;
    int id0=-1,id1=-1,id2=-1;
    int ired=-1,igrn=-1,iblu=-1;
    std::vector<int> irest;                     // f_rest_* in file order (3*Kf)
    int K=0, Kf=0;
    bool has_fdc=false, has_rgb=false;
    double rd(const uint8_t* rec, int k) const { return read_scalar_le(rec+off[k], type[k]); }
};

static bool splat_ctx_init(SplatRecCtx& ctx, const PlyElem& ve) {
    uint64_t N = ve.count;
    std::unordered_map<std::string,int> idx;
    for (size_t k=0;k<ve.props.size();k++){
        if (ve.props[k].is_list) return false;  // no list props in splat vertices
        idx[ve.props[k].name]=(int)k;
    }
    auto has=[&](const char*n){return idx.count(n)>0;};
    auto gi=[&](const char*n){auto it=idx.find(n);return it==idx.end()?-1:it->second;};
    ctx.off.resize(ve.props.size()); ctx.type.resize(ve.props.size());
    int stride=0;
    for(size_t k=0;k<ve.props.size();k++){
        ctx.off[k]=stride; ctx.type[k]=ve.props[k].type;
        stride+=ply_type_size(ve.props[k].type);
    }
    ctx.stride=stride;
    int rest=0; while(has(("f_rest_"+std::to_string(rest)).c_str())) rest++;
    ctx.Kf=rest/3;
    g_splat.count=(uint32_t)N;
    g_splat.sh_degree=(ctx.Kf>=24)?4:(ctx.Kf>=15)?3:(ctx.Kf>=8)?2:(ctx.Kf>=3)?1:0;
    static const int K_OF_DEG[5]={0,3,8,15,24};
    ctx.K=K_OF_DEG[g_splat.sh_degree];
    g_splat.sh_rest=ctx.K;
    g_splat.posop.resize(N*4); g_splat.quat.resize(N*4);
    g_splat.scale.resize(N*3); g_splat.fdc.resize(N*3);
    if(ctx.K>0) g_splat.sh.assign((size_t)N*ctx.K*3, 0);
    ctx.ix=gi("x"); ctx.iy=gi("y"); ctx.iz=gi("z");
    ctx.io=has("opacity")?gi("opacity"):gi("alpha");
    ctx.is0=gi("scale_0"); ctx.is1=gi("scale_1"); ctx.is2=gi("scale_2");
    ctx.ir0=gi("rot_0"); ctx.ir1=gi("rot_1"); ctx.ir2=gi("rot_2"); ctx.ir3=gi("rot_3");
    ctx.id0=gi("f_dc_0"); ctx.id1=gi("f_dc_1"); ctx.id2=gi("f_dc_2");
    ctx.ired=gi("red"); ctx.igrn=gi("green"); ctx.iblu=gi("blue");
    ctx.irest.resize(rest);
    for(int r=0;r<rest;r++) ctx.irest[r]=gi(("f_rest_"+std::to_string(r)).c_str());
    ctx.has_fdc=has("f_dc_0"); ctx.has_rgb=has("red")&&has("green")&&has("blue");
    return ctx.ix>=0 && ctx.iy>=0 && ctx.iz>=0 && stride>0;
}

static inline void splat_parse_record(const SplatRecCtx& c, const uint8_t* rec, uint64_t i){
    g_splat.posop[i*4+0]=(float)c.rd(rec,c.ix);
    g_splat.posop[i*4+1]=(float)c.rd(rec,c.iy);
    g_splat.posop[i*4+2]=(float)c.rd(rec,c.iz);
    g_splat.posop[i*4+3]= c.io>=0 ? sigmoidf((float)c.rd(rec,c.io)) : 1.0f;
    // quat stored wxyz (rot_0=w). Output xyzw normalized (f16).
    float qw= c.ir0>=0?(float)c.rd(rec,c.ir0):1.f;
    float qx= c.ir1>=0?(float)c.rd(rec,c.ir1):0.f;
    float qy= c.ir2>=0?(float)c.rd(rec,c.ir2):0.f;
    float qz= c.ir3>=0?(float)c.rd(rec,c.ir3):0.f;
    float qn=std::sqrt(qx*qx+qy*qy+qz*qz+qw*qw); if(qn<1e-12f)qn=1.f;
    g_splat.quat[i*4+0]=f2h(qx/qn); g_splat.quat[i*4+1]=f2h(qy/qn);
    g_splat.quat[i*4+2]=f2h(qz/qn); g_splat.quat[i*4+3]=f2h(qw/qn);
    g_splat.scale[i*3+0]= f2h(c.is0>=0?std::exp((float)c.rd(rec,c.is0)):1.f);
    g_splat.scale[i*3+1]= f2h(c.is1>=0?std::exp((float)c.rd(rec,c.is1)):1.f);
    g_splat.scale[i*3+2]= f2h(c.is2>=0?std::exp((float)c.rd(rec,c.is2)):1.f);
    if (c.has_fdc){
        g_splat.fdc[i*3+0]=f2h((float)c.rd(rec,c.id0));
        g_splat.fdc[i*3+1]=f2h((float)c.rd(rec,c.id1));
        g_splat.fdc[i*3+2]=f2h((float)c.rd(rec,c.id2));
    } else if (c.has_rgb){
        // convert display rgb (0..255) back into an f_dc so shader path is uniform
        g_splat.fdc[i*3+0]=f2h(((float)c.rd(rec,c.ired)/255.f-0.5f)/SH_C0);
        g_splat.fdc[i*3+1]=f2h(((float)c.rd(rec,c.igrn)/255.f-0.5f)/SH_C0);
        g_splat.fdc[i*3+2]=f2h(((float)c.rd(rec,c.iblu)/255.f-0.5f)/SH_C0);
    }
    if (c.K>0){
        // PLY rest is channel-major with the FILE's per-channel count:
        // [Kf red][Kf green][Kf blue]
        for (int ch=0;ch<3;ch++)
            for (int k=0;k<c.K;k++){
                int src=c.irest[ch*c.Kf+k];
                g_splat.sh[((size_t)i*c.K+k)*3+ch]= f2h(src>=0?(float)c.rd(rec,src):0.f);
            }
    }
}

static bool parse_splat_ply(const uint8_t* data, size_t len, const PlyHeader& h,
                            const PlyElem& ve) {
    uint64_t N = ve.count;
    g_splat.clear();
    g_splat.count = (uint32_t)N;
    g_splat.posop.resize(N*4);
    g_splat.quat.resize(N*4);
    g_splat.scale.resize(N*3);
    g_splat.fdc.resize(N*3);

    // Map property name -> index within the vertex record.
    std::unordered_map<std::string,int> idx;
    for (size_t k=0;k<ve.props.size();k++) idx[ve.props[k].name]=(int)k;
    auto has=[&](const char*n){return idx.count(n)>0;};
    auto gi=[&](const char*n){auto it=idx.find(n);return it==idx.end()?-1:it->second;};

    // count f_rest_*
    int rest = 0;
    while (has(("f_rest_"+std::to_string(rest)).c_str())) rest++;
    const int Kf = rest / 3;   // per channel, as stored in the file
    // degree from K: 3->1, 8->2, 15->3, 24->4
    g_splat.sh_degree = (Kf>=24)?4 : (Kf>=15)?3 : (Kf>=8)?2 : (Kf>=3)?1 : 0;
    static const int K_OF_DEG[5] = {0, 3, 8, 15, 24};
    int K = K_OF_DEG[g_splat.sh_degree];   // clamped to a complete degree
    g_splat.sh_rest = K;
    if (K>0) g_splat.sh.assign((size_t)N*K*3, 0);

    bool has_fdc = has("f_dc_0");
    bool has_rgb = has("red") && has("green") && has("blue");

    if (h.format == 1) {
        SplatRecCtx ctx;
        if (!splat_ctx_init(ctx, ve)) return false;
        const uint8_t* base = data + h.data_offset;
        if (h.data_offset + (uint64_t)ctx.stride*N > len) return false;
        for (uint64_t i=0;i<N;i++)
            splat_parse_record(ctx, base + (size_t)i*ctx.stride, i);
        return true;
    } else if (h.format == 0) {
        // ascii
        const char* p = (const char*)data + h.data_offset;
        const char* end = (const char*)data + len;
        std::vector<double> vals(ve.props.size());
        for (uint64_t i=0;i<N;i++){
            for (size_t k=0;k<ve.props.size();k++){
                char* np; vals[k]=fast_strtod(p,&np);
                if(np==p){ return i>0; } p=np;
                if(p>=end && k+1<ve.props.size()) return false;
            }
            auto v=[&](const char*n)->double{int j=gi(n);return j<0?0.0:vals[j];};
            g_splat.posop[i*4+0]=(float)v("x");
            g_splat.posop[i*4+1]=(float)v("y");
            g_splat.posop[i*4+2]=(float)v("z");
            g_splat.posop[i*4+3]= (has("opacity")||has("alpha")) ? sigmoidf((float)(has("opacity")?v("opacity"):v("alpha"))) : 1.0f;
            float qw=(float)v("rot_0"),qx=(float)v("rot_1"),qy=(float)v("rot_2"),qz=(float)v("rot_3");
            if(!has("rot_0")){qw=1;qx=qy=qz=0;}
            float qn=std::sqrt(qx*qx+qy*qy+qz*qz+qw*qw); if(qn<1e-12f)qn=1.f;
            g_splat.quat[i*4+0]=f2h(qx/qn);g_splat.quat[i*4+1]=f2h(qy/qn);g_splat.quat[i*4+2]=f2h(qz/qn);g_splat.quat[i*4+3]=f2h(qw/qn);
            g_splat.scale[i*3+0]=f2h(has("scale_0")?std::exp((float)v("scale_0")):1.f);
            g_splat.scale[i*3+1]=f2h(has("scale_1")?std::exp((float)v("scale_1")):1.f);
            g_splat.scale[i*3+2]=f2h(has("scale_2")?std::exp((float)v("scale_2")):1.f);
            if(has_fdc){g_splat.fdc[i*3+0]=f2h((float)v("f_dc_0"));g_splat.fdc[i*3+1]=f2h((float)v("f_dc_1"));g_splat.fdc[i*3+2]=f2h((float)v("f_dc_2"));}
            else if(has_rgb){g_splat.fdc[i*3+0]=f2h(((float)v("red")/255.f-0.5f)/SH_C0);g_splat.fdc[i*3+1]=f2h(((float)v("green")/255.f-0.5f)/SH_C0);g_splat.fdc[i*3+2]=f2h(((float)v("blue")/255.f-0.5f)/SH_C0);}
            if(K>0){for(int c=0;c<3;c++)for(int k=0;k<K;k++){g_splat.sh[((size_t)i*K+k)*3+c]=f2h((float)v(("f_rest_"+std::to_string(c*Kf+k)).c_str()));}}
        }
        return true;
    }
    return false; // big-endian unsupported
}

// ---------- mesh PLY ----------

// Per-record parsing context for binary-LE mesh vertices (shared by the
// whole-buffer and streaming paths). Allocates destination arrays.
struct MeshRecCtx {
    int stride = 0;
    std::vector<int> off, type;
    int ix=-1,iy=-1,iz=-1, inx=-1,iny=-1,inz=-1, ir=-1,ig=-1,ib=-1, iu=-1,iv=-1;
    bool has_n=false, has_c=false, has_uv=false;
    double rd(const uint8_t* rec, int k) const { return read_scalar_le(rec+off[k], type[k]); }
};

static bool mesh_ctx_init(MeshRecCtx& ctx, const PlyElem& ve) {
    uint64_t N = ve.count;
    std::unordered_map<std::string,int> idx;
    for (size_t k=0;k<ve.props.size();k++){
        if (ve.props[k].is_list) return false;  // vertex list props unsupported
        idx[ve.props[k].name]=(int)k;
    }
    auto has=[&](const char*n){return idx.count(n)>0;};
    auto gi=[&](const char*n){auto it=idx.find(n);return it==idx.end()?-1:it->second;};
    ctx.off.resize(ve.props.size()); ctx.type.resize(ve.props.size());
    int stride=0;
    for(size_t k=0;k<ve.props.size();k++){
        ctx.off[k]=stride; ctx.type[k]=ve.props[k].type;
        stride+=ply_type_size(ve.props[k].type);
    }
    ctx.stride=stride;
    ctx.has_n = has("nx")&&has("ny")&&has("nz");
    ctx.has_c = has("red")&&has("green")&&has("blue");
    ctx.has_uv= (has("u")&&has("v"))||(has("s")&&has("t"))||(has("texture_u")&&has("texture_v"));
    const char* un=has("u")?"u":has("s")?"s":"texture_u";
    const char* vn=has("v")?"v":has("t")?"t":"texture_v";
    g_mesh.nv=(uint32_t)N;
    g_mesh.pos.resize(N*3);
    if(ctx.has_n) g_mesh.normal.resize(N*3);
    if(ctx.has_c) g_mesh.color.resize(N*3);
    if(ctx.has_uv) g_mesh.uv.resize(N*2);
    ctx.ix=gi("x");ctx.iy=gi("y");ctx.iz=gi("z");
    ctx.inx=gi("nx");ctx.iny=gi("ny");ctx.inz=gi("nz");
    ctx.ir=gi("red");ctx.ig=gi("green");ctx.ib=gi("blue");
    ctx.iu=gi(un);ctx.iv=gi(vn);
    return ctx.ix>=0 && ctx.iy>=0 && ctx.iz>=0 && stride>0;
}

static inline void mesh_parse_vertex(const MeshRecCtx& c, const uint8_t* rec, uint64_t i){
    g_mesh.pos[i*3+0]=(float)c.rd(rec,c.ix);
    g_mesh.pos[i*3+1]=(float)c.rd(rec,c.iy);
    g_mesh.pos[i*3+2]=(float)c.rd(rec,c.iz);
    if(c.has_n){g_mesh.normal[i*3+0]=(float)c.rd(rec,c.inx);g_mesh.normal[i*3+1]=(float)c.rd(rec,c.iny);g_mesh.normal[i*3+2]=(float)c.rd(rec,c.inz);}
    if(c.has_c){g_mesh.color[i*3+0]=(uint8_t)c.rd(rec,c.ir);g_mesh.color[i*3+1]=(uint8_t)c.rd(rec,c.ig);g_mesh.color[i*3+2]=(uint8_t)c.rd(rec,c.ib);}
    if(c.has_uv){g_mesh.uv[i*2+0]=(float)c.rd(rec,c.iu);g_mesh.uv[i*2+1]=(float)c.rd(rec,c.iv);}
}

static bool parse_mesh_ply(const uint8_t* data, size_t len, const PlyHeader& h,
                           const PlyElem& ve, const PlyElem* fe) {
    uint64_t N = ve.count;
    g_mesh.clear();
    g_mesh.nv = (uint32_t)N;
    g_mesh.pos.resize(N*3);

    std::unordered_map<std::string,int> idx;
    for (size_t k=0;k<ve.props.size();k++) idx[ve.props[k].name]=(int)k;
    auto has=[&](const char*n){return idx.count(n)>0;};
    auto gi=[&](const char*n){auto it=idx.find(n);return it==idx.end()?-1:it->second;};
    bool has_n = has("nx")&&has("ny")&&has("nz");
    bool has_c = has("red")&&has("green")&&has("blue");
    bool has_uv= (has("u")&&has("v"))||(has("s")&&has("t"))||(has("texture_u")&&has("texture_v"));
    if(has_n) g_mesh.normal.resize(N*3);
    if(has_c) g_mesh.color.resize(N*3);
    if(has_uv) g_mesh.uv.resize(N*2);

    const char* un=has("u")?"u":has("s")?"s":"texture_u";
    const char* vn=has("v")?"v":has("t")?"t":"texture_v";

    if (h.format == 1) {
        MeshRecCtx ctx;
        if (!mesh_ctx_init(ctx, ve)) return false;
        const int stride = ctx.stride;
        const uint8_t* base=data+h.data_offset;
        if(h.data_offset+(uint64_t)stride*N>len) return false;
        for(uint64_t i=0;i<N;i++)
            mesh_parse_vertex(ctx, base+(size_t)i*stride, i);
        // faces
        if(fe && fe->count>0){
            // find list prop
            int lp=-1; for(size_t k=0;k<fe->props.size();k++) if(fe->props[k].is_list){lp=(int)k;break;}
            if(lp<0) return true;
            // face records may have leading/trailing scalar props; we assume single list prop
            const uint8_t* fp = base + (size_t)stride*N;
            const PlyProp& P=fe->props[lp];
            int cts=ply_type_size(P.count_type), es=ply_type_size(P.type);
            g_mesh.idx.reserve(fe->count*3);
            for(uint64_t f=0; f<fe->count; f++){
                if(fp+cts>data+len) break;
                int nverts=(int)read_scalar_le(fp,P.count_type); fp+=cts;
                if(fp+(size_t)nverts*es>data+len) break;
                // triangulate fan
                int v0=(int)read_scalar_le(fp,P.type);
                int prev=(int)read_scalar_le(fp+es,P.type);
                for(int t=2;t<nverts;t++){
                    int cur=(int)read_scalar_le(fp+(size_t)t*es,P.type);
                    g_mesh.idx.push_back((uint32_t)v0);
                    g_mesh.idx.push_back((uint32_t)prev);
                    g_mesh.idx.push_back((uint32_t)cur);
                    prev=cur;
                }
                fp+=(size_t)nverts*es;
            }
            g_mesh.nt=(uint32_t)(g_mesh.idx.size()/3);
        }
        return true;
    } else if (h.format == 0) {
        const char* p=(const char*)data+h.data_offset;
        const char* end=(const char*)data+len;
        std::vector<double> vals(ve.props.size());
        for(uint64_t i=0;i<N;i++){
            for(size_t k=0;k<ve.props.size();k++){char*np;vals[k]=fast_strtod(p,&np);if(np==p)return i>0;p=np;}
            auto v=[&](const char*n)->double{int j=gi(n);return j<0?0.0:vals[j];};
            g_mesh.pos[i*3+0]=(float)v("x");g_mesh.pos[i*3+1]=(float)v("y");g_mesh.pos[i*3+2]=(float)v("z");
            if(has_n){g_mesh.normal[i*3+0]=(float)v("nx");g_mesh.normal[i*3+1]=(float)v("ny");g_mesh.normal[i*3+2]=(float)v("nz");}
            if(has_c){g_mesh.color[i*3+0]=(uint8_t)v("red");g_mesh.color[i*3+1]=(uint8_t)v("green");g_mesh.color[i*3+2]=(uint8_t)v("blue");}
            if(has_uv){g_mesh.uv[i*2+0]=(float)v(un);g_mesh.uv[i*2+1]=(float)v(vn);}
        }
        if(fe&&fe->count>0){
            g_mesh.idx.reserve(fe->count*3);
            for(uint64_t f=0;f<fe->count;f++){
                char*np;long nverts=(long)fast_strtol(p,&np);if(np==p)break;p=np;
                std::vector<int> poly(nverts);
                for(long t=0;t<nverts;t++){poly[t]=(int)fast_strtol(p,&np);p=np;}
                for(long t=2;t<nverts;t++){g_mesh.idx.push_back(poly[0]);g_mesh.idx.push_back(poly[t-1]);g_mesh.idx.push_back(poly[t]);}
            }
            g_mesh.nt=(uint32_t)(g_mesh.idx.size()/3);
        }
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// OBJ parsing
// ---------------------------------------------------------------------------

// Flat open-addressing u64 -> u32 map for OBJ face-vertex dedup.
// std::unordered_map allocates a node per insert — for a 100+ MB OBJ that is
// tens of millions of mallocs and dominated the parse time.
struct FlatMap {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> vals;      // 0xFFFFFFFF = empty
    size_t mask = 0, count = 0;

    static inline uint64_t mix(uint64_t x) {   // splitmix64 finalizer
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
    }
    void init(size_t cap) {
        size_t n = 64;
        while (n < cap * 2) n <<= 1;
        keys.assign(n, 0);
        vals.assign(n, 0xFFFFFFFFu);
        mask = n - 1;
        count = 0;
    }
    void grow() {
        std::vector<uint64_t> ok = std::move(keys);
        std::vector<uint32_t> ov = std::move(vals);
        size_t n = (mask + 1) << 1;
        keys.assign(n, 0);
        vals.assign(n, 0xFFFFFFFFu);
        mask = n - 1;
        for (size_t i = 0; i < ov.size(); i++) {
            if (ov[i] == 0xFFFFFFFFu) continue;
            size_t j = mix(ok[i]) & mask;
            while (vals[j] != 0xFFFFFFFFu) j = (j + 1) & mask;
            keys[j] = ok[i]; vals[j] = ov[i];
        }
    }
    // returns existing value, or inserts `next` and returns it
    uint32_t find_or_insert(uint64_t key, uint32_t next, bool* inserted) {
        if ((count + 1) * 2 > mask + 1) grow();
        size_t j = mix(key) & mask;
        while (vals[j] != 0xFFFFFFFFu) {
            if (keys[j] == key) { *inserted = false; return vals[j]; }
            j = (j + 1) & mask;
        }
        keys[j] = key; vals[j] = next; count++;
        *inserted = true;
        return next;
    }
};

static bool parse_obj(const uint8_t* data, size_t len) {
    g_mesh.clear();
    std::vector<float> V, VT, VN;   // raw pools
    V.reserve(1<<16);
    // face vertex uniquification: key -> new index
    FlatMap vmap;
    vmap.init(len / 64 + 64);       // ~bytes per face line, avoids most rehashes
    std::vector<uint32_t> idx;
    std::vector<float> pos, nrm, uv;
    bool any_vt=false, any_vn=false;

    const char* p=(const char*)data;
    const char* end=(const char*)data+len;

    auto skipws=[&](){while(p<end&&(*p==' '||*p=='\t'))p++;};
    auto parsef=[&]()->float{skipws();char*np;float f=fast_strtof(p,&np);p=np;return f;};

    auto resolve=[&](long vi,long ti,long ni)->uint32_t{
        long nvp=(long)(V.size()/3);
        long nvt=(long)(VT.size()/2);
        long nvn=(long)(VN.size()/3);
        if(vi<0)vi=nvp+vi; else vi-=1;
        if(ti>0)ti-=1; else if(ti<0)ti=nvt+ti; else ti=-1;
        if(ni>0)ni-=1; else if(ni<0)ni=nvn+ni; else ni=-1;
        uint64_t key=((uint64_t)(uint32_t)vi)*2654435761u ^ ((uint64_t)(uint32_t)(ti+1)<<21) ^ ((uint64_t)(uint32_t)(ni+1)<<42);
        bool inserted;
        uint32_t ni_out=vmap.find_or_insert(key,(uint32_t)(pos.size()/3),&inserted);
        if(!inserted) return ni_out;
        if(vi>=0&&vi<nvp){pos.push_back(V[vi*3]);pos.push_back(V[vi*3+1]);pos.push_back(V[vi*3+2]);}
        else{pos.push_back(0);pos.push_back(0);pos.push_back(0);}
        if(ti>=0&&ti<nvt){uv.push_back(VT[ti*2]);uv.push_back(1.0f-VT[ti*2+1]);any_vt=true;}
        else{uv.push_back(0);uv.push_back(0);}
        if(ni>=0&&ni<nvn){nrm.push_back(VN[ni*3]);nrm.push_back(VN[ni*3+1]);nrm.push_back(VN[ni*3+2]);any_vn=true;}
        else{nrm.push_back(0);nrm.push_back(0);nrm.push_back(0);}
        return ni_out;
    };

    while(p<end){
        skipws();
        if(p>=end) break;
        char c0=*p;
        if(c0=='v'){
            char c1=(p+1<end)?p[1]:0;
            if(c1==' '||c1=='\t'){ p+=1; float x=parsef(),y=parsef(),z=parsef(); V.push_back(x);V.push_back(y);V.push_back(z); }
            else if(c1=='t'){ p+=2; float u=parsef(),v=parsef(); VT.push_back(u);VT.push_back(v); }
            else if(c1=='n'){ p+=2; float x=parsef(),y=parsef(),z=parsef(); VN.push_back(x);VN.push_back(y);VN.push_back(z); }
        } else if(c0=='f'&&(p+1<end)&&(p[1]==' '||p[1]=='\t')){
            p+=1;
            uint32_t face[64]; int fn=0;
            while(fn<64){
                skipws();
                if(p>=end||*p=='\n'||*p=='\r') break;
                // parse v[/vt][/vn]
                char* np;
                long vi=(long)fast_strtol(p,&np); if(np==p) break; p=np;
                long ti=0,ni=0;
                if(p<end&&*p=='/'){ p++;
                    if(p<end&&*p!='/'){ ti=(long)fast_strtol(p,&np); p=np; }
                    if(p<end&&*p=='/'){ p++; ni=(long)fast_strtol(p,&np); p=np; }
                }
                face[fn++]=resolve(vi,ti,ni);
            }
            for(int t=2;t<fn;t++){idx.push_back(face[0]);idx.push_back(face[t-1]);idx.push_back(face[t]);}
        }
        // skip to end of line
        while(p<end&&*p!='\n') p++;
        if(p<end) p++;
    }

    if(pos.empty()) return false;
    g_mesh.nv=(uint32_t)(pos.size()/3);
    g_mesh.nt=(uint32_t)(idx.size()/3);
    g_mesh.pos=std::move(pos);
    g_mesh.idx=std::move(idx);
    if(any_vn) g_mesh.normal=std::move(nrm);
    if(any_vt) g_mesh.uv=std::move(uv);
    return true;
}

// ---------------------------------------------------------------------------
// Public parse entry
// ---------------------------------------------------------------------------
// hint: 0 auto (by content), 1 obj
// returns kind: 0 error, 1 splat, 2 mesh
KEEP int ssv_parse(uint8_t* data, uint32_t len, int hint) {
    g_last_kind = 0;
    if (hint == 1) {
        if (parse_obj(data, len)) { g_last_kind = 2; return 2; }
        return 0;
    }
    // PLY?
    if (len >= 3 && data[0]=='p'&&data[1]=='l'&&data[2]=='y') {
        PlyHeader h;
        if (!parse_ply_header(data, len, h)) return 0;
        const PlyElem* ve=nullptr; const PlyElem* fe=nullptr;
        for (auto& e : h.elems){ if(e.name=="vertex")ve=&e; else if(e.name=="face")fe=&e; }
        if(!ve) return 0;
        bool is_mesh = (fe!=nullptr);
        // detect splat by scale_0/opacity even without faces
        std::unordered_map<std::string,int> vp; for(size_t k=0;k<ve->props.size();k++) vp[ve->props[k].name]=(int)k;
        bool looks_splat = vp.count("scale_0") && (vp.count("opacity")||vp.count("alpha")) && vp.count("rot_0");
        if (is_mesh && !(looks_splat && !fe)) {
            if(parse_mesh_ply(data,len,h,*ve,fe)){g_last_kind=2;return 2;}
            return 0;
        }
        if (looks_splat) { if(parse_splat_ply(data,len,h,*ve)){g_last_kind=1;return 1;} return 0; }
        if (fe){ if(parse_mesh_ply(data,len,h,*ve,fe)){g_last_kind=2;return 2;} }
        return 0;
    }
    return 0;
}

// Set mesh data directly from JS (used by GLTF/GLB path). JS provides typed
// arrays already; we copy into module heap so the same accessors work.
KEEP void ssv_set_mesh(uint32_t nv, uint32_t nt,
                       float* pos, float* nrm, uint8_t* col, float* uv, uint32_t* idx) {
    g_mesh.clear();
    g_mesh.nv=nv; g_mesh.nt=nt;
    g_mesh.pos.assign(pos, pos+(size_t)nv*3);
    if(nrm) g_mesh.normal.assign(nrm, nrm+(size_t)nv*3);
    if(col) g_mesh.color.assign(col, col+(size_t)nv*3);
    if(uv)  g_mesh.uv.assign(uv, uv+(size_t)nv*2);
    g_mesh.idx.assign(idx, idx+(size_t)nt*3);
    g_last_kind=2;
}

// ---------------------------------------------------------------------------
// Streaming parse (binary little-endian PLY only). For multi-GB files the
// whole-buffer path would need file + parsed arrays resident at once, which
// exceeds the wasm32 4GB heap; here JS feeds chunks and only a small carry
// buffer is retained. Usage: ssv_stream_begin(header bytes incl. end_header)
// -> 1 splat / 2 mesh / 3 unsupported (caller falls back) / 0 error; then
// ssv_stream_chunk() per chunk; ssv_stream_end() -> kind.
// ---------------------------------------------------------------------------
namespace stream {
    static int kind = 0;                  // 1 splat, 2 mesh
    static SplatRecCtx sctx;
    static MeshRecCtx mctx;
    static uint64_t vcount = 0, vdone = 0;
    static uint64_t fcount = 0, fdone = 0;
    static int f_cts = 0, f_es = 0;       // face list count/elem sizes
    static int f_ct = PT_NONE, f_et = PT_NONE;
    static bool ok = true;
    static std::vector<uint8_t> carry;    // unconsumed bytes across chunks
}

KEEP int ssv_stream_begin(uint8_t* hdr, uint32_t hlen) {
    using namespace stream;
    g_last_kind = 0; kind = 0; ok = true;
    vdone = fdone = 0; carry.clear();
    if (hlen < 3 || hdr[0]!='p'||hdr[1]!='l'||hdr[2]!='y') return 0;
    PlyHeader h;
    if (!parse_ply_header(hdr, hlen, h)) return 0;
    if (h.format != 1) return 3;          // ascii / big-endian: not streamable
    const PlyElem* ve=nullptr; const PlyElem* fe=nullptr;
    for (auto& e : h.elems){ if(e.name=="vertex")ve=&e; else if(e.name=="face")fe=&e; }
    if (!ve) return 0;
    std::unordered_map<std::string,int> vp;
    for(size_t k=0;k<ve->props.size();k++) vp[ve->props[k].name]=(int)k;
    bool looks_splat = vp.count("scale_0") && (vp.count("opacity")||vp.count("alpha")) && vp.count("rot_0");
    vcount = ve->count;
    if (looks_splat && !fe) {
        g_splat.clear();
        if (!splat_ctx_init(sctx, *ve)) return 0;
        fcount = 0;
        kind = 1;
    } else if (!fe && !looks_splat) {
        return 0;                          // plain point cloud: unsupported
    } else {
        if (fe) {
            int lp=-1; for(size_t k=0;k<fe->props.size();k++) if(fe->props[k].is_list){lp=(int)k;break;}
            if (lp<0) return 3;
            f_ct = fe->props[lp].count_type; f_et = fe->props[lp].type;
            f_cts = ply_type_size(f_ct); f_es = ply_type_size(f_et);
            fcount = fe->count;
        } else fcount = 0;
        g_mesh.clear();
        if (!mesh_ctx_init(mctx, *ve)) return 0;
        if (fcount) g_mesh.idx.reserve(fcount*3);
        kind = 2;
    }
    return kind;
}

KEEP int ssv_stream_chunk(uint8_t* data, uint32_t len) {
    using namespace stream;
    if (!kind || !ok) return -1;
    // append to carry (chunk sizes are small; one extra copy is fine)
    carry.insert(carry.end(), data, data+len);
    const uint8_t* p = carry.data();
    size_t avail = carry.size();
    size_t used = 0;
    const int stride = (kind==1) ? sctx.stride : mctx.stride;
    // vertices (fixed stride)
    while (vdone < vcount && avail - used >= (size_t)stride) {
        if (kind==1) splat_parse_record(sctx, p+used, vdone);
        else         mesh_parse_vertex(mctx, p+used, vdone);
        used += stride; vdone++;
    }
    // faces (variable-size list records)
    if (kind==2 && vdone==vcount) {
        while (fdone < fcount) {
            if (avail - used < (size_t)f_cts) break;
            int nverts=(int)read_scalar_le(p+used, f_ct);
            size_t recsz = (size_t)f_cts + (size_t)nverts*f_es;
            if (nverts < 0 || nverts > 1024) { ok=false; return -2; }
            if (avail - used < recsz) break;
            const uint8_t* fp = p+used+f_cts;
            if (nverts >= 3) {
                uint32_t v0=(uint32_t)read_scalar_le(fp, f_et);
                uint32_t prev=(uint32_t)read_scalar_le(fp+f_es, f_et);
                for(int t=2;t<nverts;t++){
                    uint32_t cur=(uint32_t)read_scalar_le(fp+(size_t)t*f_es, f_et);
                    g_mesh.idx.push_back(v0);
                    g_mesh.idx.push_back(prev);
                    g_mesh.idx.push_back(cur);
                    prev=cur;
                }
            }
            used += recsz; fdone++;
        }
    }
    carry.erase(carry.begin(), carry.begin()+used);
    return 0;
}

KEEP int ssv_stream_end() {
    using namespace stream;
    carry.clear(); carry.shrink_to_fit();
    if (!kind || !ok) { kind = 0; return 0; }
    if (vdone < vcount) { kind = 0; return 0; }   // truncated file
    if (kind==2) g_mesh.nt = (uint32_t)(g_mesh.idx.size()/3);
    g_last_kind = kind;
    int r = kind; kind = 0;
    return r;
}

// free the (large) SH rest array once it has been uploaded to the GPU;
// sorting / histograms / picking never read it
KEEP void ssv_free_splat_sh(){
    g_splat.sh.clear(); g_splat.sh.shrink_to_fit();
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------
KEEP int ssv_kind(){ return g_last_kind; }

KEEP uint32_t ssv_splat_count(){ return g_splat.count; }
KEEP int ssv_splat_sh_degree(){ return g_splat.sh_degree; }
KEEP int ssv_splat_sh_rest(){ return g_splat.sh_rest; }
KEEP float* ssv_splat_posop(){ return g_splat.posop.data(); }
// quat/scale/fdc/sh are stored as IEEE half floats (upload as HALF_FLOAT)
KEEP uint16_t* ssv_splat_quat(){ return g_splat.quat.data(); }
KEEP uint16_t* ssv_splat_scale(){ return g_splat.scale.data(); }
KEEP uint16_t* ssv_splat_fdc(){ return g_splat.fdc.data(); }
KEEP uint16_t* ssv_splat_sh(){ return g_splat.sh.data(); }

KEEP uint32_t ssv_mesh_nv(){ return g_mesh.nv; }
KEEP uint32_t ssv_mesh_nt(){ return g_mesh.nt; }
KEEP float* ssv_mesh_pos(){ return g_mesh.pos.data(); }
KEEP float* ssv_mesh_normal(){ return g_mesh.normal.empty()?nullptr:g_mesh.normal.data(); }
KEEP uint8_t* ssv_mesh_color(){ return g_mesh.color.empty()?nullptr:g_mesh.color.data(); }
KEEP float* ssv_mesh_uv(){ return g_mesh.uv.empty()?nullptr:g_mesh.uv.data(); }
KEEP uint32_t* ssv_mesh_idx(){ return g_mesh.idx.data(); }

KEEP void ssv_free_splat(){ g_splat.clear(); }
KEEP void ssv_free_mesh(){ g_mesh.clear(); }

// number of unique undirected edges (for statistics). For huge meshes the
// hash set would not fit in memory; fall back to the 2-manifold estimate
// E = 3T/2 (exact for closed triangle meshes).
KEEP uint32_t ssv_mesh_edge_count(){
    if(g_mesh.nt==0) return 0;
    if(g_mesh.nt > 20000000u) return (uint32_t)(((uint64_t)g_mesh.nt*3+1)/2);
    std::unordered_map<uint64_t,uint8_t> es;
    es.reserve(g_mesh.nt*3);
    auto add=[&](uint32_t a,uint32_t b){ if(a>b){uint32_t t=a;a=b;b=t;} es[((uint64_t)a<<32)|b]=1; };
    for(size_t f=0;f<g_mesh.nt;f++){
        uint32_t a=g_mesh.idx[f*3],b=g_mesh.idx[f*3+1],c=g_mesh.idx[f*3+2];
        add(a,b);add(b,c);add(c,a);
    }
    return (uint32_t)es.size();
}

// ---------------------------------------------------------------------------
// Depth sort (16-bit counting sort, back-to-front). Sorting depth follows
// get_sorting_depth in projection_utils.slang:
//   mode 0 (perspective/orthographic): planar depth dot(pos, dir) — fast,
//     camera position irrelevant.
//   mode 1 (equirectangular): distance from the camera (full-sphere view).
//   mode 2 (fisheye/equisolid): sqrt(sqrt(z^4 + (1/512)*(x^2+y^2)^2)) in
//     camera coordinates — continuous across the z=0 boundary so >180° views
//     with splats behind the camera sort correctly.
// (cx,cy,cz) = camera position, (dx,dy,dz) = unit forward, both in the
// model's native frame. Emits far -> near for "over" blending.
// ---------------------------------------------------------------------------
template<int mode>
static uint32_t* sort_positions(const float* P, uint32_t N,
                                float cx, float cy, float cz,
                                float dx, float dy, float dz){
    if(N==0){ g_order.clear(); g_depths.clear(); return nullptr; }
    if(g_order.size()!=N) g_order.resize(N);
    if(g_depths.size()!=N) g_depths.resize(N);
    auto depth=[&](uint32_t i)->float{
        float x=P[i*4], y=P[i*4+1], z=P[i*4+2];
        if constexpr(mode==0) return x*dx+y*dy+z*dz;
        float rx=x-cx, ry=y-cy, rz=z-cz;
        float r2=rx*rx+ry*ry+rz*rz;
        if constexpr(mode==1) return std::sqrt(r2);
        float zc=rx*dx+ry*dy+rz*dz;             // forward component
        float p2=r2-zc*zc; if(p2<0.f)p2=0.f;    // squared radial component
        float z2=zc*zc;
        return std::sqrt(std::sqrt(z2*z2 + p2*p2*(1.0f/512.0f)));
    };
    // find min/max depth
    float mn=1e30f, mx=-1e30f;
    for(uint32_t i=0;i<N;i++){
        float d=depth(i);
        if constexpr(mode==0)
            d = copysignf(std::sqrt(fabsf(d)), d);
        else
            d = std::sqrt(d);
        if(d<mn)mn=d; if(d>mx)mx=d;
        g_depths[i] = d;
    }
    float range=mx-mn; if(range<1e-12f)range=1e-12f;
    const int B=65536;
    static std::vector<uint32_t> cnt; cnt.assign(B+1,0);
    static std::vector<uint16_t> key; if(key.size()!=N) key.resize(N);
    float scale=(B-1)/range;
    for(uint32_t i=0;i<N;i++){
        float d=g_depths[i];
        int q=(int)((d-mn)*scale);
        if(q<0)q=0; if(q>=B)q=B-1;
        // reverse so that farthest (largest d) comes first
        int r=(B-1)-q;
        key[i]=(uint16_t)r;
        cnt[r+1]++;
    }
    for(int i=0;i<B;i++) cnt[i+1]+=cnt[i];
    for(uint32_t i=0;i<N;i++){
        uint16_t k=key[i];
        g_order[cnt[k]++]=i;
    }
    return g_order.data();
}

KEEP uint32_t* ssv_sort(int mode, float cx, float cy, float cz,
                        float dx, float dy, float dz){
    return (mode == 0 ? sort_positions<0> :
            mode == 1 ? sort_positions<1> :
            sort_positions<2>
        )(g_splat.posop.data(), g_splat.count,
                          cx,cy,cz, dx,dy,dz);
}

// External-positions sort, for the sort worker: the worker runs its own
// instance of this module holding only a copy of (x,y,z,opacity)*N, so the
// counting sort runs at native speed off the main thread.
static std::vector<float> g_wpos;

KEEP float* ssw_init(uint32_t n){
    g_wpos.resize((size_t)n*4);
    return g_wpos.data();
}

KEEP uint32_t* ssw_sort(int mode, float cx, float cy, float cz,
                        float dx, float dy, float dz){
    return (mode == 0 ? sort_positions<0> :
            mode == 1 ? sort_positions<1> :
            sort_positions<2>
        )(g_wpos.data(), (uint32_t)(g_wpos.size()/4),
                          cx,cy,cz, dx,dy,dz);
}

// ---------------------------------------------------------------------------
// Histograms
// ---------------------------------------------------------------------------
// param: 0 opacity, 1 log10 scale(geo mean), 2 erank, 3 r,4 g,5 b,
//        6 log10 sx,7 sy,8 sz, 9 anisotropy(log10 smax/smin)
// out_minmax[0..1] filled. bins written to g_bins (nbins).
static float splat_param(uint32_t i, int param){
    const uint16_t* S=g_splat.scale.data();
    switch(param){
        case 0: return g_splat.posop[i*4+3];
        case 1: { float a=h2f(S[i*3]),b=h2f(S[i*3+1]),c=h2f(S[i*3+2]); return std::log10(std::cbrt(a*b*c)+1e-20f);}
        case 2: {
            float sa=h2f(S[i*3]),sb=h2f(S[i*3+1]),sc=h2f(S[i*3+2]);
            float a=sa*sa, b=sb*sb, c=sc*sc;
            float s=a+b+c; if(s<1e-30f)return 1.f;
            float p0=a/s,p1=b/s,p2=c/s; float e=0;
            if(p0>1e-12f)e-=p0*std::log(p0); if(p1>1e-12f)e-=p1*std::log(p1); if(p2>1e-12f)e-=p2*std::log(p2);
            return std::exp(e);
        }
        case 3: return std::fmax(0.f,std::fmin(1.f,0.5f+SH_C0*h2f(g_splat.fdc[i*3+0])));
        case 4: return std::fmax(0.f,std::fmin(1.f,0.5f+SH_C0*h2f(g_splat.fdc[i*3+1])));
        case 5: return std::fmax(0.f,std::fmin(1.f,0.5f+SH_C0*h2f(g_splat.fdc[i*3+2])));
        case 6: return std::log10(h2f(S[i*3])+1e-20f);
        case 7: return std::log10(h2f(S[i*3+1])+1e-20f);
        case 8: return std::log10(h2f(S[i*3+2])+1e-20f);
        case 9: { float a=h2f(S[i*3]),b=h2f(S[i*3+1]),c=h2f(S[i*3+2]); float mx=std::fmax(a,std::fmax(b,c)),mn=std::fmin(a,std::fmin(b,c)); return std::log10(mx/(mn+1e-20f)+1e-20f);}
    }
    return 0;
}

KEEP float* ssv_splat_histogram(int param, int nbins){
    uint32_t N=g_splat.count;
    g_bins.assign(nbins,0.f);
    if(N==0){g_minmax[0]=0;g_minmax[1]=1;return g_bins.data();}
    float mn=1e30f,mx=-1e30f;
    for(uint32_t i=0;i<N;i++){float v=splat_param(i,param);if(std::isfinite(v)){if(v<mn)mn=v;if(v>mx)mx=v;}}
    if(!(mx>mn)){mx=mn+1.f;}
    g_minmax[0]=mn;g_minmax[1]=mx;
    float sc=nbins/(mx-mn);
    for(uint32_t i=0;i<N;i++){float v=splat_param(i,param);if(!std::isfinite(v))continue;int b=(int)((v-mn)*sc);if(b<0)b=0;if(b>=nbins)b=nbins-1;g_bins[b]+=1.f;}
    return g_bins.data();
}

// mesh param: 0 log10 edge length, 1 log10 tri area, 2 x,3 y,4 z coord
KEEP float* ssv_mesh_histogram(int param, int nbins){
    g_bins.assign(nbins,0.f);
    const float* P=g_mesh.pos.data();
    std::vector<float> vals;
    if(param==0){
        if (g_mesh.nt > 20000000u) {
            // huge mesh: skip dedup (hash set would not fit); interior edges
            // are counted twice but the distribution shape is unchanged
            vals.reserve((size_t)g_mesh.nt*3);
            auto edge=[&](uint32_t a,uint32_t b){float dx=P[a*3]-P[b*3],dy=P[a*3+1]-P[b*3+1],dz=P[a*3+2]-P[b*3+2];vals.push_back(std::log10(std::sqrt(dx*dx+dy*dy+dz*dz)+1e-20f));};
            for(size_t f=0;f<g_mesh.nt;f++){uint32_t a=g_mesh.idx[f*3],b=g_mesh.idx[f*3+1],c=g_mesh.idx[f*3+2];edge(a,b);edge(b,c);edge(c,a);}
        } else {
            std::unordered_map<uint64_t,uint8_t> es; es.reserve(g_mesh.nt*3);
            auto edge=[&](uint32_t a,uint32_t b){if(a>b){uint32_t t=a;a=b;b=t;}uint64_t k=((uint64_t)a<<32)|b;if(es.emplace(k,1).second){float dx=P[a*3]-P[b*3],dy=P[a*3+1]-P[b*3+1],dz=P[a*3+2]-P[b*3+2];vals.push_back(std::log10(std::sqrt(dx*dx+dy*dy+dz*dz)+1e-20f));}};
            for(size_t f=0;f<g_mesh.nt;f++){uint32_t a=g_mesh.idx[f*3],b=g_mesh.idx[f*3+1],c=g_mesh.idx[f*3+2];edge(a,b);edge(b,c);edge(c,a);}
        }
    } else if(param==1){
        vals.reserve(g_mesh.nt);
        for(size_t f=0;f<g_mesh.nt;f++){uint32_t a=g_mesh.idx[f*3],b=g_mesh.idx[f*3+1],c=g_mesh.idx[f*3+2];
            float ux=P[b*3]-P[a*3],uy=P[b*3+1]-P[a*3+1],uz=P[b*3+2]-P[a*3+2];
            float vx=P[c*3]-P[a*3],vy=P[c*3+1]-P[a*3+1],vz=P[c*3+2]-P[a*3+2];
            float cx=uy*vz-uz*vy,cy=uz*vx-ux*vz,cz=ux*vy-uy*vx;
            vals.push_back(std::log10(0.5f*std::sqrt(cx*cx+cy*cy+cz*cz)+1e-20f));}
    } else {
        int comp=param-2; vals.reserve(g_mesh.nv);
        for(uint32_t i=0;i<g_mesh.nv;i++) vals.push_back(P[i*3+comp]);
    }
    if(vals.empty()){g_minmax[0]=0;g_minmax[1]=1;return g_bins.data();}
    float mn=1e30f,mx=-1e30f; for(float v:vals){if(std::isfinite(v)){if(v<mn)mn=v;if(v>mx)mx=v;}}
    if(!(mx>mn))mx=mn+1.f; g_minmax[0]=mn;g_minmax[1]=mx;
    float sc=nbins/(mx-mn);
    for(float v:vals){if(!std::isfinite(v))continue;int b=(int)((v-mn)*sc);if(b<0)b=0;if(b>=nbins)b=nbins-1;g_bins[b]+=1.f;}
    return g_bins.data();
}
KEEP float* ssv_histogram_minmax(){ return g_minmax; }

// Bounding box of current model (min xyz, max xyz) -> 6 floats
static float g_bbox[6];
KEEP float* ssv_bbox(){
    float mn[3]={1e30f,1e30f,1e30f},mx[3]={-1e30f,-1e30f,-1e30f};
    if(g_last_kind==1){
        const float*P=g_splat.posop.data();
        for(uint32_t i=0;i<g_splat.count;i++)for(int k=0;k<3;k++){float v=P[i*4+k];if(v<mn[k])mn[k]=v;if(v>mx[k])mx[k]=v;}
    } else if(g_last_kind==2){
        const float*P=g_mesh.pos.data();
        for(uint32_t i=0;i<g_mesh.nv;i++)for(int k=0;k<3;k++){float v=P[i*3+k];if(v<mn[k])mn[k]=v;if(v>mx[k])mx[k]=v;}
    }
    for(int k=0;k<3;k++){g_bbox[k]=mn[k];g_bbox[k+3]=mx[k];}
    return g_bbox;
}

// Robust fitting sphere of the current model: per-axis median center and the
// median distance from it -> 4 floats (cx, cy, cz, median_dist). Trained
// splats and generated meshes often contain far-away outliers that blow up a
// bounding-box fit; medians ignore them. Subsamples uniformly for speed.
static float g_fitsph[4];
KEEP float* ssv_fit_sphere(){
    const float* P = nullptr; uint32_t N = 0, stride = 3;
    if (g_last_kind == 1) { P = g_splat.posop.data(); N = g_splat.count; stride = 4; }
    else if (g_last_kind == 2) { P = g_mesh.pos.data(); N = g_mesh.nv; stride = 3; }
    if (!P || N == 0) { g_fitsph[0]=g_fitsph[1]=g_fitsph[2]=0.f; g_fitsph[3]=1.f; return g_fitsph; }
    const uint32_t MAXS = 1u<<20;
    uint32_t step = (N + MAXS - 1) / MAXS; if (step < 1) step = 1;
    uint32_t M = (N + step - 1) / step;
    std::vector<float> tmp; tmp.reserve(M);
    for (int k = 0; k < 3; k++) {
        tmp.clear();
        for (uint32_t i = 0; i < N; i += step) tmp.push_back(P[(size_t)i*stride+k]);
        std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
        g_fitsph[k] = tmp[tmp.size()/2];
    }
    tmp.clear();
    for (uint32_t i = 0; i < N; i += step) {
        float dx=P[(size_t)i*stride]-g_fitsph[0], dy=P[(size_t)i*stride+1]-g_fitsph[1], dz=P[(size_t)i*stride+2]-g_fitsph[2];
        tmp.push_back(dx*dx+dy*dy+dz*dz);
    }
    std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
    g_fitsph[3] = std::sqrt(tmp[tmp.size()/2]);
    return g_fitsph;
}

// Brute-force ray/mesh intersection (Moller-Trumbore over all triangles) for
// double-click view centering. Ray in the model's native frame, dir need not
// be normalized (t is in units of |dir|). Returns nearest t > 0, or -1.
KEEP float ssv_raycast_mesh(float ox, float oy, float oz, float dx, float dy, float dz){
    if (g_mesh.nt == 0) return -1.f;
    const float* P = g_mesh.pos.data();
    const uint32_t* I = g_mesh.idx.data();
    float best = 3.4e38f;
    for (size_t f = 0; f < g_mesh.nt; f++){
        const float *a = P + (size_t)I[f*3]*3, *b = P + (size_t)I[f*3+1]*3, *c = P + (size_t)I[f*3+2]*3;
        float e1x=b[0]-a[0], e1y=b[1]-a[1], e1z=b[2]-a[2];
        float e2x=c[0]-a[0], e2y=c[1]-a[1], e2z=c[2]-a[2];
        float px=dy*e2z-dz*e2y, py=dz*e2x-dx*e2z, pz=dx*e2y-dy*e2x;
        float det=e1x*px+e1y*py+e1z*pz;
        if (std::fabs(det) < 1e-20f) continue;
        float inv=1.0f/det;
        float tx=ox-a[0], ty=oy-a[1], tz=oz-a[2];
        float u=(tx*px+ty*py+tz*pz)*inv;
        if (u<0.f||u>1.f) continue;
        float qx=ty*e1z-tz*e1y, qy=tz*e1x-tx*e1z, qz=tx*e1y-ty*e1x;
        float v=(dx*qx+dy*qy+dz*qz)*inv;
        if (v<0.f||u+v>1.f) continue;
        float t=(e2x*qx+e2y*qy+e2z*qz)*inv;
        if (t>1e-6f && t<best) best=t;
    }
    return best < 3.4e38f ? best : -1.f;
}

// allocate/free raw input buffer helpers (JS streams file bytes into here)
KEEP uint8_t* ssv_alloc(uint32_t n){ return (uint8_t*)malloc(n); }
KEEP void ssv_free(void* p){ free(p); }
