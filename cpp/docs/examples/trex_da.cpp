// Example program of the guide page "T-Rex+DA: dependency-aware dummies". The regions
// between the //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_da/trex_da.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::trex_selector_methods::ad::trex_da::TRexDASelector;
using trex::trex_selector_methods::ad::trex_da::TRexDAControlParameter;
using trex::trex_selector_methods::ad::trex_da::DAMethod;
using trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;
namespace rates = trex::utils::eval::rates;
//! [includes]

//! [helpers]
void print_indices(const char* label, const std::vector<std::size_t>& indices) {
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}

const char* da_method_name(DAMethod method) {
    switch (method) {
        case DAMethod::AR1: return "AR1";
        case DAMethod::EQUI: return "EQUI";
        case DAMethod::BT: return "BT";
        case DAMethod::NN: return "NN";
        case DAMethod::PRIOR_GROUPS: return "PRIOR_GROUPS";
    }
    return "?";
}
//! [helpers]

int main() {
    DocSections doc("trex_da", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(346);
    std::normal_distribution<double> normal(0.0, 1.0);
    const auto draw = [&] { return normal(gen); };
    const Eigen::Index n = 100, p = 100;
    const double rho = 0.7;

    Eigen::MatrixXd X(n, p);
    X.col(0) = Eigen::VectorXd::NullaryExpr(n, draw);
    for (Eigen::Index j = 1; j < p; ++j) {
        X.col(j) = rho * X.col(j - 1)
                 + std::sqrt(1.0 - rho * rho) * Eigen::VectorXd::NullaryExpr(n, draw);
    }
    const std::vector<std::size_t> true_support = {9, 29, 49, 69, 89};
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (std::size_t j : true_support) beta[static_cast<Eigen::Index>(j)] = 2.0;
    Eigen::VectorXd y = X * beta + Eigen::VectorXd::NullaryExpr(n, draw);
    //! [data]

    doc.begin("compare");
    //! [compare]
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    TRexControlParameter ctrl;
    ctrl.K = 10;

    TRexSelector plain(X_map, y_map, 0.1, ctrl, /*seed=*/7, /*verbose=*/false);
    plain.select();

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::BT;
    da_ctrl.trex_ctrl = ctrl;
    TRexDASelector da(X_map, y_map, 0.1, da_ctrl, 7, false);
    da.select();

    const auto report = [&](const char* name, const TRexSelector& s) {
        const auto& idx = s.getSelectedIndices();
        std::cout << std::left << std::setw(10) << name << std::right
                  << std::setw(9) << idx.size()
                  << std::setw(7) << rates::compute_fdp(idx, true_support)
                  << std::setw(7) << rates::compute_tpp(idx, true_support) << '\n';
    };
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::left << std::setw(10) << "selector" << std::right
              << std::setw(9) << "selected" << std::setw(7) << "FDP"
              << std::setw(7) << "TPP\n";
    report("T-Rex", plain);
    report("T-Rex+DA", da);
    //! [compare]

    doc.begin("da-result");
    //! [da-result]
    const auto& da_result = da.getDAResult();
    std::cout << "method:     " << da_method_name(da_result.method) << '\n';
    std::cout << "rho_thresh: " << std::setprecision(4) << da_result.rho_thresh << '\n';
    std::cout << "rho_grid:   " << std::setprecision(3)
              << da_result.rho_grid.transpose() << '\n';
    //! [da-result]

    doc.begin("ar1");
    //! [ar1]
    TRexDAControlParameter ar1_ctrl;
    ar1_ctrl.method = DAMethod::AR1;
    ar1_ctrl.trex_ctrl = ctrl;
    TRexDASelector ar1(X_map, y_map, 0.1, ar1_ctrl, 7, false);
    ar1.select();
    print_indices("selected:", ar1.getSelectedIndices());
    std::cout << "cor_coef: " << std::setprecision(3) << ar1.getDAResult().cor_coef << '\n';
    //! [ar1]

    doc.begin("prior-groups");
    //! [prior-groups]
    std::vector<Eigen::Index> blocks(p);                       // ten blocks of ten variables
    for (Eigen::Index j = 0; j < p; ++j) blocks[j] = j / 10 + 1;

    TRexDAControlParameter pg_ctrl;
    pg_ctrl.prior_groups = {blocks};
    pg_ctrl.trex_ctrl = ctrl;
    TRexDASelector pg(X_map, y_map, 0.1, pg_ctrl, 7, false);
    pg.select();
    print_indices("selected:", pg.getSelectedIndices());
    std::cout << "method:   " << da_method_name(pg.getDAResult().method) << '\n';
    //! [prior-groups]

    doc.begin("greedy");
    //! [greedy]
    TRexControlParameter greedy_ctrl;
    greedy_ctrl.solver_type = SolverTypeForTRex::TOMP;
    greedy_ctrl.K = 10;
    greedy_ctrl.solver_params.exch_tie_alpha = 0.25;

    TRexDAControlParameter bt_ctrl;
    bt_ctrl.method = DAMethod::BT;
    bt_ctrl.trex_ctrl = greedy_ctrl;
    TRexDASelector greedy(X_map, y_map, 0.1, bt_ctrl, 7, false);
    greedy.select();
    print_indices("selected:", greedy.getSelectedIndices());
    //! [greedy]

    return 0;
}
