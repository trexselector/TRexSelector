// Example program of the guide page "T-Rex Tikhonov". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_power_rule.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_tikhonov::TRexTikhonovSelector;
using trex::trex_selector_methods::ad::trex_tikhonov::TRexTikhonovControlParameter;
using trex::trex_selector_methods::ad::trex_tikhonov::TikPenaltyStructure;
using trex::trex_selector_methods::ad::trex_tikhonov::Lambda2PowerRuleConfig;
using trex::trex_selector_methods::ad::trex_tikhonov::tikhonov_lambda2_power_rule;
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
    DocSections doc("trex_tikhonov", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(2718);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 120, p = 60;

    std::vector<std::size_t> true_support;
    for (std::size_t j = 20; j < 30; ++j) true_support.push_back(j);
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.segment(20, 10).setOnes();

    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    //! [data]

    doc.begin("penalty");
    //! [penalty]
    // Gamma: first differences, (p - 1) x p
    std::vector<Eigen::Triplet<double>> entries;
    for (Eigen::Index i = 0; i + 1 < p; ++i) {
        entries.emplace_back(i, i, -1.0);
        entries.emplace_back(i, i + 1, 1.0);
    }
    Eigen::SparseMatrix<double> gamma(p - 1, p);
    gamma.setFromTriplets(entries.begin(), entries.end());

    Eigen::SparseMatrix<double> K = TRexTikhonovSelector::gammaToK(gamma);
    std::cout << "K: " << K.rows() << " x " << K.cols() << ", " << K.nonZeros()
              << " nonzeros\n";
    std::cout << Eigen::MatrixXd(K).topLeftCorner(4, 4) << '\n';
    //! [penalty]

    doc.begin("run");
    //! [run]
    TRexTikhonovControlParameter tik_ctrl;
    tik_ctrl.tikhonov_K = K;
    tik_ctrl.lambda_2 = 1.0;
    tik_ctrl.trex_ctrl.K = 10;

    TRexTikhonovSelector tik(X_map, y_map, /*tFDR=*/0.1, tik_ctrl,
                             /*seed=*/4, /*verbose=*/false);
    tik.select();
    const std::vector<std::size_t>& selected = tik.getSelectedIndices();
    print_indices("selected:", selected);
    std::cout << "FDP: " << rates::compute_fdp(selected, true_support) << '\n';
    std::cout << "TPP: " << rates::compute_tpp(selected, true_support) << '\n';
    std::cout << "lambda_2 used: " << tik.getLambda2Used() << '\n';
    //! [run]

    doc.begin("pilot");
    //! [pilot]
    TRexTikhonovControlParameter pilot_ctrl;
    pilot_ctrl.tikhonov_K = K;
    pilot_ctrl.pilot_grid = {0.0, 1.0, 10.0};
    pilot_ctrl.pilot_K = 4;
    pilot_ctrl.pilot_T = 2;
    pilot_ctrl.trex_ctrl.K = 10;

    TRexTikhonovSelector pilot(X_map, y_map, 0.1, pilot_ctrl, 4, false);
    pilot.select();
    const auto& pilot_result = pilot.getLambda2PilotResult();
    std::cout << "lambda_2 used: " << pilot.getLambda2Used() << '\n';
    std::cout << "grid:  ";
    for (double l2 : pilot_result.grid) std::cout << ' ' << l2;
    std::cout << "\ncounts:";
    for (std::size_t c : pilot_result.counts) std::cout << ' ' << c;
    std::cout << '\n';
    //! [pilot]

    doc.begin("power-rule");
    //! [power-rule]
    TRexTikhonovControlParameter rule_ctrl;
    rule_ctrl.tikhonov_K = K;
    rule_ctrl.trex_ctrl.K = 10;
    Lambda2PowerRuleConfig cfg;
    cfg.grid = {0.0, 1.0, 10.0};

    const auto profile = tikhonov_lambda2_power_rule(X_map, y_map, 0.1, rule_ctrl, 4, cfg);
    std::cout << "lambda_2: " << profile.lambda2 << '\n';
    std::cout << "counts:  ";
    for (std::size_t c : profile.counts()) std::cout << ' ' << c;
    std::cout << '\n';
    //! [power-rule]

    doc.begin("elastic-net");
    //! [elastic-net]
    Eigen::SparseMatrix<double> identity(p, p);
    identity.setIdentity();

    TRexTikhonovControlParameter en_ctrl;
    en_ctrl.tikhonov_K = identity;
    en_ctrl.lambda_2 = 1.0;
    en_ctrl.cluster_dummies = false;
    en_ctrl.fold_dummy_coupling = false;
    en_ctrl.trex_ctrl.K = 10;

    TRexTikhonovSelector en(X_map, y_map, 0.1, en_ctrl, 4, false);
    en.select();
    print_indices("selected:", en.getSelectedIndices());
    //! [elastic-net]

    doc.begin("group-centering");
    //! [group-centering]
    std::vector<Eigen::Index> groups(p);
    for (Eigen::Index j = 0; j < p; ++j) groups[j] = j / 10;   // six groups of ten

    TRexTikhonovControlParameter gc_ctrl;
    gc_ctrl.penalty_structure = TikPenaltyStructure::GROUP_CENTERING;
    gc_ctrl.prior_groups = groups;
    gc_ctrl.lambda_2 = 1.0;
    gc_ctrl.trex_ctrl.K = 10;

    TRexTikhonovSelector gc(X_map, y_map, 0.1, gc_ctrl, 4, false);
    gc.select();
    print_indices("selected:", gc.getSelectedIndices());
    //! [group-centering]

    doc.begin("solver");
    //! [solver]
    TRexTikhonovControlParameter pw_ctrl = tik_ctrl;
    pw_ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TIKNET;

    TRexTikhonovSelector pw(X_map, y_map, 0.1, pw_ctrl, 4, false);
    pw.select();
    print_indices("selected:", pw.getSelectedIndices());
    //! [solver]

    return 0;
}
