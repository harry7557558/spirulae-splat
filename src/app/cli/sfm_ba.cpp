// `ssplat-sfm ba` -- the bundle adjuster driven directly on a BAL problem
// (Bundle Adjustment in the Large), for benchmarking and solver debugging.
// Its own translation unit because it shares nothing with the pipeline
// subcommands but the solver. See src/sfm/ba/README.md.
#include <chrono>
#include <cstring>
#include <fstream>
#include <random>

#include "sfm/ba/Problem.h"
#include "sfm/ba/Solver.h"

static void writePly(const char* path, const std::vector<double>& pts) {
    std::ofstream f(path, std::ios::binary);
    size_t n = pts.size() / 3;
    f << "ply\nformat binary_little_endian 1.0\nelement vertex " << n
      << "\nproperty float x\nproperty float y\nproperty float z\nend_header\n";
    for (size_t i = 0; i < n * 3; i++) {
        float v = (float)pts[i];
        f.write((const char*)&v, 4);
    }
}


int cmdBa(int argc, char** argv) {
    std::string file, ply_prefix, loss = "trivial", model = "snavely";
    SolverOptions opt;
    bool shared_intr = false;

    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--real") {
            std::string r = next();
            opt.real = r == "float" ? RealCfg::F32 : r == "df" ? RealCfg::DF64 : RealCfg::F64;
        } else if (a == "--loss") loss = next();
        else if (a == "--spv-path") opt.spv_path = next();
        else if (a == "--loss-param") opt.loss_param = std::stof(next());
        else if (a == "--model") model = next();
        else if (a == "--shared-intrinsics") shared_intr = true;
        else if (a == "--max-iters") opt.max_iters = std::stoi(next());
        else if (a == "--damping") opt.init_damping = std::stod(next());
        else if (a == "--rtol") opt.rtol = std::stod(next());
        else if (a == "--patience") opt.patience = std::stoi(next());
        else if (a == "--ply") ply_prefix = next();
        else if (a == "--solver") {
            std::string s = next();
            opt.solver = s == "dense" ? SolverSel::Dense
                       : s == "cg"    ? SolverSel::CG
                                      : SolverSel::Auto;
        } else if (a == "--vram-budget") opt.vram_budget_mb = std::stod(next());
        else if (a == "--cg-iters") opt.cg_max_iters = std::stoi(next());
        else if (a == "--cg-tol") opt.cg_tol = std::stod(next());
        else if (a == "--cg-fallback") {
            std::string s = next();
            opt.cg_fallback = s == "on"  ? CgFallback::On
                            : s == "off" ? CgFallback::Off
                                         : CgFallback::Auto;
        }
        else if (a == "--device") opt.device = std::stoi(next());
        else if (a == "--validate") opt.validate = true;
        else if (a == "--profile") opt.profile = true;
        else if (a == "--quiet") opt.verbose = false;
        else if (a[0] != '-') file = a;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }

    // The kernels are compiled into the binary; --spv-path overrides.
    opt.loss = loss;

    if (file.empty()) {
        fprintf(stderr,
            "usage: %s <bal_problem.txt> [--real float|double|df] [--loss trivial|huber|cauchy]\n"
            "  [--loss-param X] [--model snavely|snavely_f] [--shared-intrinsics]\n"
            "  [--max-iters N] [--damping X] [--rtol X] [--patience N] [--ply prefix]\n"
            "  [--solver auto|dense|cg] [--vram-budget MB] [--cg-iters N] [--cg-tol X]\n"
            "  [--cg-fallback auto|on|off]\n"
            "  [--device I] [--validate] [--quiet]\n"
            "  [--spv-path FILE]\n", "ssplat-sfm ba");
        return 1;
    }

    int model_id = -1;
    for (int m = 0; m < kNumModels; m++)
        if (model == kModels[m].name) model_id = m;
    if (model_id < 0) { fprintf(stderr, "unknown model %s\n", model.c_str()); return 1; }

    // both debug hooks pin the solver selection they need
    if (getenv("SSPLAT_SFM_DUMP_SG")) opt.solver = SolverSel::Dense;
    const bool cmp_step = getenv("SSPLAT_SFM_CMP_STEP") != nullptr;
    if (cmp_step) {
        opt.solver = SolverSel::CG;
        opt.cg_fallback = CgFallback::On;
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    BAProblem P = loadBAL(file, model_id, shared_intr);
    BundleSolver solver(P, opt);
    solver.init();
    auto t1 = std::chrono::high_resolution_clock::now();
    double t_pre = std::chrono::duration<double>(t1 - t0).count();

    if (!ply_prefix.empty()) writePly((ply_prefix + "_before.ply").c_str(), P.points);

    if (getenv("SSPLAT_SFM_DUMP_SG")) {
        solver.debugAssemble(atof(getenv("SSPLAT_SFM_DUMP_SG_LAMBDA") ?: "0.01"));
        std::vector<uint8_t> raw(solver.bufS().size);
        std::vector<double> v;
        solver.ctx().download(solver.bufS(), raw.data(), solver.bufS().size);
        unpackReals(v, raw.data(), (uint64_t)P.n_dim * (P.n_dim + 1) / 2, opt.real);
        std::string base = getenv("SSPLAT_SFM_DUMP_SG");
        std::ofstream fs(base + "_S.bin", std::ios::binary);
        fs.write((const char*)v.data(), v.size() * 8);
        solver.ctx().download(solver.bufG(), raw.data(), solver.bufG().size);
        unpackReals(v, raw.data(), P.n_dim, opt.real);
        std::ofstream fg(base + "_g.bin", std::ios::binary);
        fg.write((const char*)v.data(), v.size() * 8);
        return 0;
    }

    // SSPLAT_SFM_CMP_STEP: solve one assembly with both CG and dense, print the
    // step difference (use --cg-tol / --cg-iters to control CG accuracy)
    if (cmp_step) {
        double lam = atof(getenv("SSPLAT_SFM_CMP_STEP_LAMBDA") ?: "0.01");
        solver.debugCompareStep((float)lam);
        return 0;
    }

    solver.solve();
    solver.downloadParams();

    if (!ply_prefix.empty()) writePly((ply_prefix + "_after.ply").c_str(), P.points);

    const SolverStats& st = solver.stats();
    printf("real=%s loss=%s model=%s solver=%s%s\n", realCfgName(opt.real), loss.c_str(),
           model.c_str(), st.solver, shared_intr ? " shared-intrinsics" : "");
    printf("initial cost: %.6e\n", st.initial_cost);
    printf("final cost:   %.6e\n", st.final_cost);
    printf("iterations:   %d (%d accepted)\n", st.iterations, st.accepted);
    if (st.cg_solves)
        printf("cg: %.1f iters/solve avg, %d dense fallbacks\n",
               st.cg_iters_total / st.cg_solves, st.cg_fallbacks);
    printf("time: preprocess %.3f s, solve %.3f s, total %.3f s\n", t_pre, st.solve_seconds,
           t_pre + st.solve_seconds);
    printf("vram: %.1f MB\n", st.vram_mb);
    return 0;
}
