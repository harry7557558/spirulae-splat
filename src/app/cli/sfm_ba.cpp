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


// Printed by `ssplat-sfm ba --help` and by `ssplat-sfm help ba`. The solver's
// options are not in the SfmConfig table: nothing else drives them, and the
// pipeline's own bundle adjustment is configured through --ba-* on `map`.
void printBaHelp(FILE* out) {
    std::fprintf(out,
        "ssplat-sfm ba -- bundle-adjust a BAL problem (solver benchmark)\n"
        "\n"
        "Usage:\n"
        "  ssplat-sfm ba <BAL_PROBLEM.TXT> [options]\n"
        "\n"
        "Description:\n"
        "  Runs the GPU bundle adjuster directly on a problem in Bundle Adjustment in the\n"
        "  Large format, and reports cost, iterations, time and VRAM. This is how the\n"
        "  solver is benchmarked and debugged against a published reference; the pipeline\n"
        "  itself never reads BAL. See src/sfm/ba/README.md.\n"
        "\n"
        "Options:\n"
        "  --real {float|double|df}            [double]  scalar the solver works in; df is an\n"
        "                                      emulated double-float, for devices without fp64\n"
        "  --loss {trivial|huber|cauchy}       [trivial] robust loss\n"
        "  --loss-param X                      [1]       Huber delta / Cauchy c\n"
        "  --model {snavely|snavely_f}         [snavely] BAL camera model\n"
        "  --shared-intrinsics                 one intrinsics block for every camera\n"
        "  --solver {auto|dense|cg}            [auto]    dense Cholesky or implicit-Schur PCG\n"
        "  --max-iters N                       Levenberg-Marquardt iteration cap\n"
        "  --damping X                         initial LM damping\n"
        "  --rtol X                            relative cost improvement to stop below\n"
        "  --patience N                        steps below --rtol before stopping\n"
        "  --cg-iters N                        PCG iteration cap per LM step\n"
        "  --cg-tol X                          PCG relative residual tolerance\n"
        "  --cg-fallback {auto|on|off}         fall back to dense when PCG stalls\n"
        "  --vram-budget MB                    device memory the solver may use\n"
        "  --ply PREFIX                        write PREFIX_before.ply and PREFIX_after.ply\n"
        "  --device N                          Vulkan device index\n"
        "  --validate                          check the assembled system against a reference\n"
        "  --profile                           per-kernel timings\n"
        "  --quiet                             only the result lines\n"
        "  --spv-path FILE                     load the solver kernels from this SPIR-V file\n"
        "  -h, --help                          show this help and exit\n"
        "\n"
        "Examples:\n"
        "  ssplat-sfm ba problem-49-7776-pre.txt --real df --loss huber\n"
        "  ssplat-sfm ba problem-1778-993923-pre.txt --solver cg --vram-budget 4096\n"
        "\n"
        "A (--real, --loss) pair that was trimmed out of the build reports \"variant not\n"
        "built into this binary\"; see SSPLAT_SFM_REALS / SSPLAT_SFM_LOSSES.\n");
}

int cmdBa(int argc, char** argv) {
    std::string file, ply_prefix, loss = "trivial", model = "snavely";
    SolverOptions opt;
    bool shared_intr = false;

    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--help" || a == "-h") { printBaHelp(stdout); return 0; }
        else if (a == "--real") {
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
        else {
            fprintf(stderr, "ssplat-sfm ba: error: unknown option %s\n", a.c_str());
            fprintf(stderr, "Try 'ssplat-sfm ba --help' for more information.\n");
            return 1;
        }
    }

    // The kernels are compiled into the binary; --spv-path overrides.
    opt.loss = loss;

    if (file.empty()) {
        fprintf(stderr, "ssplat-sfm ba: error: a BAL problem file is required\n");
        fprintf(stderr, "Try 'ssplat-sfm ba --help' for more information.\n");
        return 1;
    }

    int model_id = -1;
    for (int m = 0; m < kNumModels; m++)
        if (model == kModels[m].name) model_id = m;
    if (model_id < 0) {
        fprintf(stderr, "ssplat-sfm ba: error: unknown --model %s\n", model.c_str());
        return 1;
    }

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
