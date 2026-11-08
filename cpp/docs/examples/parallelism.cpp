// Example program of the guide page "Parallelism". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>
#include <utils/openmp/utils_openmp.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::utils::openmp::get_max_threads;
using trex::utils::openmp::get_num_threads;
using trex::utils::openmp::set_num_threads;
namespace rates = trex::utils::eval::rates;
//! [includes]

int main() {
    DocSections doc("parallelism", DOC_OUTPUT_DIR);

    doc.begin("budget");
    //! [budget]
    // The two-thread rule of this guide, applied for the whole program
    const int old_threads = set_num_threads(std::min(2, get_max_threads()));
    std::cout << get_num_threads() << '\n';
    std::cout << get_max_threads() << '\n';
    //! [budget]

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(12);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, p = 150;
    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.head(5).setOnes();
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    //! [data]

    doc.begin("parallel");
    //! [parallel]
    TRexControlParameter serial_ctrl;
    serial_ctrl.K = 10;
    TRexSelector serial(X_map, y_map, 0.1, serial_ctrl, /*seed=*/42, /*verbose=*/false);
    serial.select();

    TRexControlParameter outer_ctrl = serial_ctrl;
    outer_ctrl.parallel_rnd_experiments = true;
    outer_ctrl.max_outer_threads = 2;
    TRexSelector parallel(X_map, y_map, 0.1, outer_ctrl, 42, false);
    parallel.select();

    TRexControlParameter inner_ctrl = serial_ctrl;
    inner_ctrl.max_inner_threads = 2;
    TRexSelector budget(X_map, y_map, 0.1, inner_ctrl, 42, false);
    budget.select();

    for (std::size_t j : serial.getSelectedIndices()) std::cout << j << ' ';
    std::cout << '\n' << std::boolalpha;
    std::cout << (serial.getSelectedIndices() == parallel.getSelectedIndices()) << '\n';
    std::cout << (serial.getSelectedIndices() == budget.getSelectedIndices()) << '\n';
    //! [parallel]

    doc.begin("mc");
    //! [mc]
    const int trials = 20;
    const std::vector<std::size_t> true_support = {0, 1, 2, 3, 4};
    std::vector<double> fdp(trials);

    #pragma omp parallel for num_threads(2) schedule(dynamic)
    for (int mc = 0; mc < trials; ++mc) {
        std::mt19937_64 gen_mc(static_cast<unsigned>(mc + 1));
        std::normal_distribution<double> normal_mc(0.0, 1.0);
        const auto draw = [&] { return normal_mc(gen_mc); };
        Eigen::MatrixXd X_mc = Eigen::MatrixXd::NullaryExpr(n, p, draw);
        Eigen::VectorXd y_mc = X_mc * beta + Eigen::VectorXd::NullaryExpr(n, draw);
        Eigen::Map<Eigen::MatrixXd> X_mc_map(X_mc.data(), n, p);
        Eigen::Map<Eigen::VectorXd> y_mc_map(y_mc.data(), n);

        TRexControlParameter ctrl;   // serial inside the worker, the default
        ctrl.K = 10;
        // seed -1: every trial draws its own dummies
        TRexSelector sel(X_mc_map, y_mc_map, 0.1, ctrl, /*seed=*/-1, false);
        sel.select();
        fdp[static_cast<std::size_t>(mc)] =
            rates::compute_fdp(sel.getSelectedIndices(), true_support);
    }

    double mean_fdp = 0.0;
    for (double f : fdp) mean_fdp += f / trials;
    std::cout << std::fixed << std::setprecision(3)
              << "mean FDP over " << trials << " trials: " << mean_fdp << '\n';
    //! [mc]

    doc.begin("restore-budget");
    //! [restore-budget]
    set_num_threads(old_threads);
    //! [restore-budget]

    return 0;
}
