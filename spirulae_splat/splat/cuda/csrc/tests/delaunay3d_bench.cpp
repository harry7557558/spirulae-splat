// Standalone benchmark / correctness harness for Delaunay3D.
//
// Reads Gaussian means (x,y,z) from a 3DGS .ply file, builds the 3D Delaunay
// triangulation and reports wall-clock time and peak RSS.
//
// Build:
//   g++ -O3 -std=c++17 -pthread delaunay3d_bench.cpp ../Delaunay3D.cpp \
//       -o /tmp/delaunay3d_bench
// Run:
//   /tmp/delaunay3d_bench path/to/splat.ply [num_threads] [max_points]

#include "../Delaunay3D.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <sstream>
#include <thread>
#include <cmath>
#include <sys/resource.h>

struct PlyProp { std::string name; int size; };

static int type_size(const std::string& t) {
    if(t=="float"||t=="float32"||t=="int"||t=="int32"||t=="uint"||t=="uint32") return 4;
    if(t=="double"||t=="float64") return 8;
    if(t=="uchar"||t=="char"||t=="uint8"||t=="int8") return 1;
    if(t=="short"||t=="ushort"||t=="int16"||t=="uint16") return 2;
    return 4;
}

// Read x,y,z floats from a binary_little_endian ply.
static bool read_ply_means(const std::string& path, std::vector<double>& out,
                           size_t max_points) {
    std::ifstream f(path, std::ios::binary);
    if(!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); return false; }
    std::string line;
    size_t nverts = 0;
    std::vector<PlyProp> props;
    bool in_vertex = false;
    // header
    while(std::getline(f, line)) {
        if(!line.empty() && line.back()=='\r') line.pop_back();
        std::istringstream iss(line);
        std::string tok; iss >> tok;
        if(tok=="format") {
            std::string fmt; iss >> fmt;
            if(fmt!="binary_little_endian") {
                fprintf(stderr, "only binary_little_endian supported (%s)\n", fmt.c_str());
                return false;
            }
        } else if(tok=="element") {
            std::string ename; size_t n; iss >> ename >> n;
            in_vertex = (ename=="vertex");
            if(in_vertex) nverts = n;
        } else if(tok=="property" && in_vertex) {
            std::string t; iss >> t;
            if(t=="list") { fprintf(stderr,"list prop in vertex unsupported\n"); return false; }
            std::string pname; iss >> pname;
            props.push_back({pname, type_size(t)});
        } else if(tok=="end_header") {
            break;
        }
    }
    int off_x=-1, off_y=-1, off_z=-1, stride=0;
    for(auto& p : props) {
        if(p.name=="x") off_x=stride;
        else if(p.name=="y") off_y=stride;
        else if(p.name=="z") off_z=stride;
        stride += p.size;
    }
    if(off_x<0||off_y<0||off_z<0) { fprintf(stderr,"no x/y/z\n"); return false; }

    size_t n = nverts;
    if(max_points && n>max_points) n=max_points;
    out.resize(n*3);
    std::vector<char> buf(stride);
    for(size_t i=0;i<n;++i) {
        f.read(buf.data(), stride);
        if(!f) { fprintf(stderr,"unexpected EOF at %zu\n", i); return false; }
        float x,y,z;
        std::memcpy(&x, buf.data()+off_x, 4);
        std::memcpy(&y, buf.data()+off_y, 4);
        std::memcpy(&z, buf.data()+off_z, 4);
        out[3*i+0]=x; out[3*i+1]=y; out[3*i+2]=z;
    }
    fprintf(stderr, "read %zu / %zu vertices, stride=%d\n", n, nverts, stride);
    return true;
}

static long peak_rss_kb() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss; // KB on Linux
}

int main(int argc, char** argv) {
    if(argc < 2) {
        fprintf(stderr,"usage: %s file.ply [num_threads] [max_points]\n", argv[0]);
        return 1;
    }
    int num_threads = (argc>=3) ? atoi(argv[2]) : 0;
    size_t max_points = (argc>=4) ? (size_t)atoll(argv[3]) : 0;

    std::vector<double> pts;
    if(!read_ply_means(argv[1], pts, max_points)) return 1;
    int n = (int)(pts.size()/3);

    long rss_before = peak_rss_kb();
    auto t0 = std::chrono::high_resolution_clock::now();
    auto res = delaunay3d::compute_delaunay_3d(pts.data(), n, num_threads, true);
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1-t0).count();
    long rss_after = peak_rss_kb();

    double pts_mb = (double)pts.size()*sizeof(double)/1e6;
    printf("points=%d  threads=%d  tets=%d  time=%.3fs  rate=%.2fM pts/s\n",
           n, res.num_threads, res.nb_cells, secs, n/secs/1e6);
    printf("peak_RSS=%.1f MB  (delta over input load=%.1f MB)  input_pts=%.1f MB\n",
           rss_after/1024.0, (rss_after-rss_before)/1024.0, pts_mb);

    // light sanity check: vertex indices in range, Euler-ish ratio
    long bad=0;
    for(size_t i=0;i<res.cell_vertices.size();++i) {
        int v=res.cell_vertices[i];
        if(v<0||v>=n) ++bad;
    }
    printf("sanity: out_of_range_vertex_refs=%ld  tets/points=%.2f\n",
           bad, (double)res.nb_cells/n);
    return 0;
}
