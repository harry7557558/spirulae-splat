// The from-scratch GPU dense Cholesky, against a CPU reference on a random SPD
// system. One run per scalar configuration compiled into the binary.
//
//   sfm_cholesky_test [N] [--real float|double|df] [--device I]
//
// Prints PASS/FAIL and returns 0/1. See docs/testing.md.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "sfm/ba/Problem.h"
#include "sfm/ba/Solver.h"

// Validate the from-scratch GPU Cholesky against a CPU reference on a random
// SPD system of dimension n.
int selftestChol(uint32_t n, SolverOptions opt) {
    BAProblem P;  // empty problem, only n_dim used
    P.n_dim = n;
    BundleSolver solver(P, opt);
    solver.init();

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss;
    std::vector<double> A((size_t)n * n, 0.0), b(n);
    {
        std::vector<double> M((size_t)n * n);
        for (auto& v : M) v = gauss(rng);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = 0; j <= i; j++) {
                double s = 0;
                for (uint32_t k = 0; k < n; k++) s += M[(size_t)i * n + k] * M[(size_t)j * n + k];
                A[(size_t)i * n + j] = A[(size_t)j * n + i] = s + (i == j ? (double)n : 0.0);
            }
        for (auto& v : b) v = gauss(rng);
    }

    std::vector<double> packed((size_t)n * (n + 1) / 2);
    for (uint32_t i = 0; i < n; i++)
        for (uint32_t j = 0; j <= i; j++)
            packed[(size_t)i * (i + 1) / 2 + j] = A[(size_t)i * n + j];
    std::vector<uint8_t> tmp;
    packReals(tmp, packed.data(), packed.size(), opt.real);
    solver.ctx().upload(solver.bufS(), tmp.data(), tmp.size());
    packReals(tmp, b.data(), n, opt.real);
    solver.ctx().upload(solver.bufG(), tmp.data(), tmp.size());

    VkCommandBuffer cb = solver.ctx().begin();
    solver.recordCholesky(cb);
    solver.ctx().submit(cb);

    std::vector<uint8_t> raw(n * realSize(opt.real));
    solver.ctx().download(solver.bufG(), raw.data(), raw.size());
    std::vector<double> x;
    unpackReals(x, raw.data(), n, opt.real);

    // CPU reference solve (LLT via simple Cholesky)
    std::vector<double> L = A;
    for (uint32_t j = 0; j < n; j++) {
        for (uint32_t k = 0; k < j; k++)
            for (uint32_t i = j; i < n; i++)
                L[(size_t)i * n + j] -= L[(size_t)i * n + k] * L[(size_t)j * n + k];
        double d = std::sqrt(L[(size_t)j * n + j]);
        for (uint32_t i = j; i < n; i++) L[(size_t)i * n + j] /= d;
    }
    std::vector<double> y = b;
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) y[i] -= L[(size_t)i * n + j] * y[j];
        y[i] /= L[(size_t)i * n + i];
    }
    for (int i = (int)n - 1; i >= 0; i--) {
        for (uint32_t j = i + 1; j < n; j++) y[i] -= L[(size_t)j * n + i] * y[j];
        y[i] /= L[(size_t)i * n + i];
    }

    double maxRel = 0, maxAbs = 0;
    for (uint32_t i = 0; i < n; i++) {
        double e = std::fabs(x[i] - y[i]);
        maxAbs = std::max(maxAbs, e);
        maxRel = std::max(maxRel, e / std::max(1e-30, std::fabs(y[i])));
    }
    printf("selftest-chol n=%u real=%s: max abs err %.3e, max rel err %.3e\n", n,
           realCfgName(opt.real), maxAbs, maxRel);
    double tol = opt.real == RealCfg::F32 ? 1e-2 : 1e-7;
    printf("%s\n", maxRel < tol ? "PASS" : "FAIL");
    return maxRel < tol ? 0 : 1;
}

int main(int argc, char** argv) {
    uint32_t n = 500;
    SolverOptions opt;
    opt.verbose = false;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--real") {
            std::string s = next();
            opt.real = s == "float" ? RealCfg::F32 : s == "double" ? RealCfg::F64 : RealCfg::DF64;
        } else if (a == "--device") {
            opt.device = std::stoi(next());
        } else if (a[0] != '-') {
            n = (uint32_t)std::stoul(a);
        } else {
            fprintf(stderr, "unknown arg %s\n", a.c_str());
            return 2;
        }
    }
    return selftestChol(n, opt);
}
