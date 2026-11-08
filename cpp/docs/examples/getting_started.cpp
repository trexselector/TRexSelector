// Example program of the guide page "Getting started". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;
namespace rates = trex::utils::eval::rates;
//! [includes]

//! [helpers]
void print_indices(const char* label, const std::vector<std::size_t>& indices) {
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}
//! [helpers]

int main() {
    DocSections doc("getting_started", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(123);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 50, p = 100, num_act = 3;

    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.head(num_act).setOnes();
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;

    std::vector<std::size_t> true_support;
    for (Eigen::Index j = 0; j < p; ++j)
        if (beta[j] != 0.0) true_support.push_back(static_cast<std::size_t>(j));
    print_indices("true support:", true_support);
    //! [data]

    doc.begin("first-run");
    //! [first-run]
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);

    TRexSelector sel(X_map, y_map, /*tFDR=*/0.1, TRexControlParameter{},
                     /*seed=*/1234, /*verbose=*/false);
    sel.select();
    print_indices("selected:", sel.getSelectedIndices());
    //! [first-run]

    doc.begin("metrics");
    //! [metrics]
    const std::vector<std::size_t>& selected = sel.getSelectedIndices();
    std::cout << "FDP: " << rates::compute_fdp(selected, true_support) << '\n';
    std::cout << "TPP: " << rates::compute_tpp(selected, true_support) << '\n';
    //! [metrics]

    doc.begin("result-fields");
    //! [result-fields]
    std::cout << "selected: " << sel.getNumSelected() << '\n';       // selected variables
    std::cout << "T_stop:   " << sel.getTStop() << '\n';             // stopping time
    std::cout << "L:        " << sel.getDummyMultiplierL() << '\n';  // L * p dummies
    std::cout << "v_thresh: " << sel.getVotingThreshold() << '\n';   // voting level
    std::cout << "K:        " << sel.getK() << '\n';                 // random experiments
    //! [result-fields]

    doc.begin("matrices");
    //! [matrices]
    const Eigen::MatrixXd& fdp_hat = sel.getFDPHatMat();   // stopping times x voting levels
    const Eigen::Index last = static_cast<Eigen::Index>(sel.getTStop()) - 1;
    std::cout << "FDP_hat_mat: " << fdp_hat.rows() << " x " << fdp_hat.cols() << '\n';
    std::cout << "voting grid: " << sel.getVotingGrid().head(6).transpose() << '\n';
    std::cout << "last row:    " << ((fdp_hat.row(last).array() * 1000.0).round() / 1000.0).matrix() << '\n';
    //! [matrices]

    doc.begin("solver");
    //! [solver]
    TRexControlParameter ctrl;
    ctrl.solver_type = SolverTypeForTRex::TLASSO;
    ctrl.K = 10;

    TRexSelector sel2(X_map, y_map, 0.1, ctrl, 1234, false);
    sel2.select();
    print_indices("selected:", sel2.getSelectedIndices());
    //! [solver]

    return 0;
}
