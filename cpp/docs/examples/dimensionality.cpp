// Example program of the guide page "FDR control and dimensionality". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
namespace rates = trex::utils::eval::rates;
//! [includes]

int main() {
    DocSections doc("dimensionality", DOC_OUTPUT_DIR);

    doc.begin("regimes");
    //! [regimes]
    std::mt19937_64 gen(7);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, num_act = 5;
    const int trials = 10;
    const std::vector<std::size_t> true_support = {0, 1, 2, 3, 4};

    for (Eigen::Index p : {50, 200, 1000}) {
        double fdr = 0.0, tpr = 0.0;
        for (int trial = 0; trial < trials; ++trial) {
            Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
            Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
            beta.head(num_act).setOnes();
            Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
            Eigen::VectorXd y = X * beta + noise;
            Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
            Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);

            TRexControlParameter ctrl;
            ctrl.K = 10;
            TRexSelector sel(X_map, y_map, 0.1, ctrl, 7 + trial, false);
            sel.select();
            fdr += rates::compute_fdp(sel.getSelectedIndices(), true_support) / trials;
            tpr += rates::compute_tpp(sel.getSelectedIndices(), true_support) / trials;
        }
        std::cout << std::fixed << std::setprecision(2)
                  << "p = " << std::setw(4) << p
                  << "  n/p = " << std::setw(5) << double(n) / double(p)
                  << "  FDR = " << fdr << "  TPR = " << tpr << '\n';
    }
    //! [regimes]

    return 0;
}
