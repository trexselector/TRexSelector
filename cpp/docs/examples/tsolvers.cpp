// Example program of the guide page "The T-solvers". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>
#include <tsolvers/ad/linear_model/afs_based.hpp>
#include <tsolvers/ad/linear_model/cd_based.hpp>
#include <tsolvers/ad/linear_model/lars_based.hpp>
#include <tsolvers/ad/linear_model/omp_based.hpp>

using trex::ml_methods::scaler_methods::LpNormScaler;
using trex::tsolvers::ad::SolverOptions;
namespace lars = trex::tsolvers::ad::linear_model::lars_based;
namespace omp = trex::tsolvers::ad::linear_model::omp_based;
namespace afs = trex::tsolvers::ad::linear_model::afs_based;
namespace cd = trex::tsolvers::ad::linear_model::cd_based;
//! [includes]

//! [helpers]
void print_indices(const char* label, const std::vector<std::size_t>& indices) {
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}

// Constructs a reader of the given solver type, walks its path to three dummies and prints
// the number of steps and the active predictors.
template <class Solver, class... Args>
void run_reader(const char* name, Args&... args) {
    Solver solver(args...);
    solver.executeStep(3);
    std::cout << std::left << std::setw(11) << name << std::setw(6) << solver.getNumSteps();
    print_indices(" ", solver.getActivePredictorIndices());
}

// Writes a snapshot of a path for the figure of the page: the dense coefficients, one row per
// step with the all-zero start first, and the code of the variable that enters at each step.
void write_snapshot(const lars::TLARS_Solver& solver, const std::string& tag) {
    const Eigen::MatrixXd path = solver.getBetaPathSparse().dense();   // (p + L) x (steps + 1)
    const std::string base = std::string(DOC_OUTPUT_DIR) + "/tsolvers.figure_" + tag;
    std::ofstream csv(base + ".csv");
    for (Eigen::Index s = 0; s < path.cols(); ++s) {
        for (Eigen::Index j = 0; j < path.rows(); ++j) csv << (j ? "," : "") << path(j, s);
        csv << '\n';
    }
    std::ofstream actions(base + ".actions.txt");
    for (const auto& step : solver.getActions()) actions << step.front() << '\n';
}
//! [helpers]

int main() {
    DocSections doc("tsolvers", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(31);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, p = 50, L = 50;
    const std::vector<std::size_t> true_support = {1, 8, 16, 32};

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta[1] = 2.0;
    beta[8] = -2.0;
    beta[16] = 1.5;
    beta[32] = -1.5;
    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::MatrixXd D = Eigen::MatrixXd::NullaryExpr(n, L, [&] { return normal(gen); });
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), n, L);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    //! [data]

    doc.begin("path");
    //! [path]
    lars::TLARS_Solver solver(X_map, D_map, y_map);
    solver.executeStep(/*T_stop=*/3);
    std::cout << "steps:             " << solver.getNumSteps() << '\n';
    print_indices("active predictors:", solver.getActivePredictorIndices());
    print_indices("active dummies:   ", solver.getActiveDummyIndices());
    //! [path]

    doc.begin("actives");
    //! [actives]
    print_indices("actives:          ", solver.getActives());
    std::cout << "dummy start index: " << solver.getDummyStartIndex() << '\n';
    std::cout << "first actions:    ";
    const auto& actions = solver.getActions();   // one list of signed codes per step
    for (std::size_t s = 0; s < 4 && s < actions.size(); ++s) {
        for (int code : actions[s]) std::cout << ' ' << code;
    }
    std::cout << '\n';
    //! [actives]

    doc.begin("beta");
    //! [beta]
    const auto path = solver.getBetaPathSparse();
    std::cout << "beta path: " << path.steps.size() << " steps of "
              << path.p_total << " coefficients\n";
    const Eigen::VectorXd beta_hat = solver.getBeta();   // the last step
    std::cout << std::fixed << std::setprecision(3) << "beta on the support:";
    for (std::size_t j : true_support) {
        std::cout << ' ' << beta_hat[static_cast<Eigen::Index>(j)];
    }
    std::cout << '\n';

    const std::vector<double>& lam = solver.getLambda();   // one value per step
    const std::size_t k = lam.size();                       // the others start at step zero
    const auto rss = solver.getRSS(), r2 = solver.getR2(), cp = solver.getCp();
    const auto dof = solver.getDoF();
    const auto tail = [k](const auto& v, std::size_t s) {
        return v[v.size() - k + s];
    };
    std::cout << "lambda, rss, r2, cp, dof per step:\n";
    for (std::size_t s = 0; s < 5 && s < k; ++s) {
        std::cout << lam[s] << "  " << tail(rss, s) << "  " << tail(r2, s) << "  "
                  << tail(cp, s) << "  " << tail(dof, s) << '\n';
    }
    solver.restore(X_map, D_map, y_map);   // give X and D back before the next solver
    //! [beta]

    doc.begin("figure");
    //! [figure]
    lars::TLARS_Solver s(X_map, D_map, y_map);
    s.executeStep(/*T_stop=*/1);
    write_snapshot(s, "T1");
    s.executeStep(/*T_stop=*/3);
    write_snapshot(s, "T3");
    s.restore(X_map, D_map, y_map);
    std::ofstream meta(std::string(DOC_OUTPUT_DIR) + "/tsolvers.figure.meta.txt");
    meta << "p " << p << "\ntrue_support";
    for (std::size_t j : true_support) meta << ' ' << j;
    meta << '\n';
    std::cout << "snapshots written, " << s.getNumSteps() << " steps after the continuation\n";
    //! [figure]

    doc.begin("options");
    //! [options]
    lars::TENET_Solver en(X_map, D_map, y_map, SolverOptions{.lambda2 = 0.5});
    en.executeStep(3);
    print_indices("active predictors:", en.getActivePredictorIndices());
    en.restore(X_map, D_map, y_map);
    //! [options]

    doc.begin("writers");
    //! [writers]
    const Eigen::MatrixXd X_before = X;
    const auto max_diff = [&] { return (X - X_before).cwiseAbs().maxCoeff(); };

    lars::TLASSO_Solver w(X_map, D_map, y_map);
    std::cout << std::scientific << std::setprecision(2);
    std::cout << "the scaled columns:  " << max_diff() << '\n';
    w.restore(X_map, D_map, y_map);
    std::cout << "the original values: " << max_diff() << '\n';
    //! [writers]

    doc.begin("readers");
    //! [readers]
    LpNormScaler sx(LpNormScaler::NormType::L2, /*center=*/true);
    sx.fit_transform_inplace(X_map);
    LpNormScaler sd(LpNormScaler::NormType::L2, /*center=*/true);
    sd.fit_transform_inplace(D_map);
    Eigen::VectorXd y_c = y.array() - y.mean();
    Eigen::Map<Eigen::VectorXd> yc_map(y_c.data(), n);

    SolverOptions reader_opts{.normalize = false, .intercept = false, .precentered = true};
    //! [readers]

    doc.begin("families");
    //! [families]
    std::cout << std::left << std::setw(11) << "solver" << std::setw(6) << "steps"
              << "  predictors\n";
    run_reader<lars::TLARS_Solver>("TLARS", X_map, D_map, yc_map, reader_opts);
    run_reader<lars::TLASSO_Solver>("TLASSO", X_map, D_map, yc_map, reader_opts);
    run_reader<lars::TSTEPWISE_Solver>("TSTEPWISE", X_map, D_map, yc_map, reader_opts);
    run_reader<lars::TSTAGEWISE_Solver>("TSTAGEWISE", X_map, D_map, yc_map, reader_opts);
    run_reader<lars::TENET_Solver>("TENET", X_map, D_map, yc_map, reader_opts);
    run_reader<lars::TENETAug_Solver>("TENETAug", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TOMP_Solver>("TOMP", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TOOLS_Solver>("TOOLS", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TGP_Solver>("TGP", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TACGP_Solver>("TACGP", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TMP_Solver>("TMP", X_map, D_map, yc_map, reader_opts);
    run_reader<omp::TNCGMP_Solver>("TNCGMP", X_map, D_map, yc_map, reader_opts);
    run_reader<afs::TAFS_Solver>("TAFS", X_map, D_map, yc_map, reader_opts);
    run_reader<cd::TCDLASSO_Solver>("TCDLASSO", X_map, D_map, yc_map, reader_opts);
    run_reader<cd::TCDENET_Solver>("TCDENET", X_map, D_map, yc_map, reader_opts);
    //! [families]

    doc.begin("structured");
    //! [structured]
    Eigen::VectorXi groups(p);
    for (Eigen::Index j = 0; j < p; ++j) groups[j] = static_cast<int>(j / 10);   // 5 groups
    Eigen::SparseMatrix<double> identity(p, p);
    identity.setIdentity();
    SolverOptions structured_opts{.normalize = false, .intercept = false,
                                  .precentered = true, .lambda2 = 0.5};

    run_reader<lars::TIENET_Solver>("TIENET", X_map, D_map, yc_map, groups, structured_opts);
    run_reader<cd::TCDTIKNET_Solver>("TCDTIKNET", X_map, D_map, yc_map, identity,
                                     structured_opts);
    //! [structured]

    doc.begin("save-load");
    //! [save-load]
    const auto work_dir = std::filesystem::temp_directory_path() / "trex_guide_tsolvers";
    std::filesystem::create_directories(work_dir);
    const std::string state_file = (work_dir / "tlars_state.bin").string();

    lars::TLARS_Solver a(X_map, D_map, yc_map, reader_opts);
    a.executeStep(2);
    a.save(state_file);
    lars::TLARS_Solver b = lars::TLARS_Solver::load(state_file, X_map, D_map);
    std::cout << "steps after load: " << b.getNumSteps() << '\n';
    b.executeStep(4);
    lars::TLARS_Solver c(X_map, D_map, yc_map, reader_opts);
    c.executeStep(4);
    std::cout << std::boolalpha << "same actives: " << (b.getActives() == c.getActives())
              << '\n';
    //! [save-load]

    doc.begin("memory");
    //! [memory]
    std::cout << "bytes held: " << a.memoryBytes() << '\n';
    //! [memory]

    doc.begin("cleanup");
    //! [cleanup]
    std::filesystem::remove_all(work_dir);
    //! [cleanup]

    return 0;
}
