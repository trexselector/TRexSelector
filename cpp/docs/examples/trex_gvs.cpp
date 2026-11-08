// Example program of the guide page "T-Rex GVS: grouped variable selection". The regions
// between the //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_gvs::TRexGVSSelector;
using trex::trex_selector_methods::ad::trex_gvs::TRexGVSControlParameter;
using trex::trex_selector_methods::ad::trex_gvs::GVSType;
using trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;
namespace rates = trex::utils::eval::rates;
//! [includes]

//! [helpers]
void print_sorted(const char* label, std::vector<std::size_t> indices) {
    std::sort(indices.begin(), indices.end());
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}
//! [helpers]

int main() {
    DocSections doc("trex_gvs", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(11);
    std::normal_distribution<double> normal(0.0, 1.0);
    const auto draw = [&] { return normal(gen); };
    const Eigen::Index n = 100, p = 100, group_size = 5;
    const Eigen::Index num_groups = p / group_size;

    Eigen::VectorXi groups(p);   // 0-based group label per variable
    for (Eigen::Index j = 0; j < p; ++j) groups[j] = static_cast<int>(j / group_size);

    Eigen::MatrixXd latent = Eigen::MatrixXd::NullaryExpr(n, num_groups, draw);
    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, draw);
    X *= std::sqrt(1.0 - 0.9 * 0.9);
    for (Eigen::Index j = 0; j < p; ++j) X.col(j) += 0.9 * latent.col(groups[j]);

    std::vector<std::size_t> true_support;   // the members of groups 0 and 1
    for (Eigen::Index j = 0; j < p; ++j) {
        if (groups[j] <= 1) true_support.push_back(static_cast<std::size_t>(j));
    }
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (std::size_t j : true_support) beta[static_cast<Eigen::Index>(j)] = 1.0;
    Eigen::VectorXd y = X * beta + Eigen::VectorXd::NullaryExpr(n, draw);
    //! [data]

    doc.begin("en");
    //! [en]
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);

    TRexGVSControlParameter gvs_ctrl;
    gvs_ctrl.gvs_type = GVSType::EN;
    gvs_ctrl.lambda_2 = 1.0;
    gvs_ctrl.trex_ctrl.K = 10;

    TRexGVSSelector gvs(X_map, y_map, 0.1, gvs_ctrl, /*seed=*/5, /*verbose=*/false);
    gvs.select();
    const auto& en_selected = gvs.getSelectedIndices();
    print_sorted("selected:", en_selected);
    std::cout << "FDP: " << rates::compute_fdp(en_selected, true_support) << '\n';
    std::cout << "TPP: " << rates::compute_tpp(en_selected, true_support) << '\n';
    //! [en]

    doc.begin("en-result");
    //! [en-result]
    const auto& res = gvs.getGVSResult();
    std::cout << "gvs_type:       " << (res.gvs_type == GVSType::EN ? "EN" : "IEN") << '\n';
    std::cout << "max_clusters:   " << res.max_clusters << '\n';
    std::cout << "hc_method_used: " << res.hc_method_used << '\n';
    std::cout << "lambda2_used:   " << res.lambda2_used << '\n';

    std::map<int, int> cluster_sizes;   // label -> number of members
    for (int label : res.groups_vec) ++cluster_sizes[label];
    std::cout << "first clusters:";
    int shown = 0;
    for (const auto& [label, size] : cluster_sizes) {
        if (shown++ == 5) break;
        std::cout << ' ' << label << ':' << size;
    }
    std::cout << '\n';
    //! [en-result]

    doc.begin("cv");
    //! [cv]
    TRexGVSControlParameter cv_ctrl;
    cv_ctrl.gvs_type = GVSType::EN;
    cv_ctrl.cv_n_lambda = 100;
    cv_ctrl.trex_ctrl.K = 10;

    TRexGVSSelector cv(X_map, y_map, 0.1, cv_ctrl, 5, false);
    cv.select();
    const auto& diag = cv.getLambda2CVDiagnostics();
    std::cout << std::boolalpha << std::fixed << std::setprecision(3);
    std::cout << "lambda2_used:      " << cv.getGVSResult().lambda2_used << '\n';
    std::cout << "ran:               " << diag.ran << '\n';
    std::cout << "seconds:           " << diag.seconds << '\n';
    std::cout << "n_profile_at_edge: " << diag.n_profile_at_edge << '\n';
    std::cout << "grid position:     " << diag.grid_index << " of " << diag.grid_size << '\n';
    std::cout << "boundary_pick:     " << diag.boundary_pick << '\n';
    //! [cv]

    doc.begin("ien");
    //! [ien]
    TRexGVSControlParameter ien_ctrl;
    ien_ctrl.gvs_type = GVSType::IEN;
    ien_ctrl.prior_groups = std::vector<Eigen::Index>(groups.begin(), groups.end());
    ien_ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TIENET;
    ien_ctrl.trex_ctrl.K = 10;

    TRexGVSSelector ien(X_map, y_map, 0.1, ien_ctrl, 5, false);
    ien.select();
    const auto& ien_selected = ien.getSelectedIndices();
    print_sorted("selected:", ien_selected);
    std::cout << std::defaultfloat << std::setprecision(4);
    std::cout << "lambda2_used: " << ien.getGVSResult().lambda2_used << '\n';
    std::cout << "FDP: " << rates::compute_fdp(ien_selected, true_support) << '\n';
    std::cout << "TPP: " << rates::compute_tpp(ien_selected, true_support) << '\n';
    //! [ien]

    return 0;
}
