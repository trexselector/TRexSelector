// Example program of the guide page "The T-Rex selector". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <utils/datageneration/utils_dummygen.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::trex_selector_methods::ad::trex_core::LLoopStrategy;
using trex::trex_selector_methods::ad::utils::data_normalizer::ScalingMode;
namespace dummygen = trex::utils::datageneration::dummygen;
namespace rates = trex::utils::eval::rates;
//! [includes]

namespace {

void print_indices(const char* label, const std::vector<std::size_t>& indices) {
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}

}  // namespace

int main() {
    DocSections doc("trex_selector", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(123);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, p = 150;

    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.head(5).setOnes();
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;
    const std::vector<std::size_t> true_support = {0, 1, 2, 3, 4};
    //! [data]

    doc.begin("run");
    //! [run]
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);

    TRexSelector sel(X_map, y_map, /*tFDR=*/0.1, TRexControlParameter{},
                     /*seed=*/42, /*verbose=*/false);
    sel.select();
    print_indices("selected:", sel.getSelectedIndices());
    //! [run]

    doc.begin("voting-grid");
    //! [voting-grid]
    std::cout << "voting grid: " << sel.getVotingGrid().transpose() << '\n';
    //! [voting-grid]

    doc.begin("result-fields");
    //! [result-fields]
    std::cout << "T_stop:   " << sel.getTStop() << '\n';
    std::cout << "L:        " << sel.getDummyMultiplierL() << '\n';
    std::cout << "v_thresh: " << sel.getVotingThreshold() << '\n';
    //! [result-fields]

    doc.begin("matrices");
    //! [matrices]
    const Eigen::MatrixXd& fdp_hat = sel.getFDPHatMat();   // stopping times x voting levels
    const Eigen::MatrixXd& r_mat = sel.getRMat();          // selected variables per cell
    const Eigen::Index last = static_cast<Eigen::Index>(sel.getTStop()) - 1;

    std::cout << "FDP_hat_mat: " << fdp_hat.rows() << " x " << fdp_hat.cols() << '\n';
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "last row:    " << fdp_hat.row(last) << '\n';
    std::cout << std::setprecision(0);
    std::cout << "R_mat, last row: " << r_mat.row(last) << '\n';
    //! [matrices]

    doc.begin("phi");
    //! [phi]
    const Eigen::VectorXd phi = sel.getPhiMat().row(last);   // after the last stopping time
    std::cout << std::setprecision(2) << "phi of the true support:";
    for (std::size_t j : true_support) {
        std::cout << ' ' << phi[static_cast<Eigen::Index>(j)];
    }
    std::cout << '\n';

    int at_voting_level = 0, below_half = 0;
    for (Eigen::Index j = 5; j < p; ++j) {
        if (phi[j] >= sel.getVotingThreshold()) ++at_voting_level;
        if (phi[j] < 0.5) ++below_half;
    }
    std::cout << "null variables at the voting level: " << at_voting_level << '\n';
    std::cout << "null variables below 0.5: " << below_half << " of " << p - 5 << '\n';

    std::ofstream phi_csv(std::string(DOC_OUTPUT_DIR) + "/trex_selector.figure_phi.csv");
    for (Eigen::Index j = 0; j < p; ++j) {   // the data of the figure: variable, phi, true
        phi_csv << j << ',' << phi[j] << ',' << (j < 5 ? 1 : 0) << '\n';
    }
    std::ofstream phi_meta(std::string(DOC_OUTPUT_DIR) + "/trex_selector.figure_phi.meta.txt");
    phi_meta << "v_thresh " << sel.getVotingThreshold() << '\n';
    //! [phi]

    doc.begin("lloop");
    //! [lloop]
    struct Named { LLoopStrategy strategy; const char* name; };
    for (Named s : {Named{LLoopStrategy::HCONCAT, "HCONCAT"},
                    Named{LLoopStrategy::STANDARD, "STANDARD"},
                    Named{LLoopStrategy::SEEDED, "SEEDED"}}) {
        TRexControlParameter ctrl;
        ctrl.lloop_strategy = s.strategy;
        ctrl.K = 10;

        TRexSelector s_sel(X_map, y_map, 0.1, ctrl, 42, false);
        s_sel.select();
        std::cout << std::left << std::setw(10) << s.name
                  << " L = " << s_sel.getDummyMultiplierL()
                  << "  T_stop = " << s_sel.getTStop() << "  selected:";
        for (std::size_t j : s_sel.getSelectedIndices()) std::cout << ' ' << j;
        std::cout << '\n';
    }
    //! [lloop]

    doc.begin("dummies");
    //! [dummies]
    TRexControlParameter t_ctrl;
    t_ctrl.dummy_distribution = dummygen::Distribution::StudentT(3.0);
    t_ctrl.K = 10;

    TRexSelector t_sel(X_map, y_map, 0.1, t_ctrl, 42, false);
    t_sel.select();
    print_indices("selected:", t_sel.getSelectedIndices());
    //! [dummies]

    doc.begin("scaling");
    //! [scaling]
    // By hand what the selector does: centered unit-L2 columns, centered response
    Eigen::MatrixXd Xn = X.rowwise() - X.colwise().mean();
    Xn.array().rowwise() /= Xn.colwise().norm().array();
    Eigen::VectorXd yn = y.array() - y.mean();
    Eigen::Map<Eigen::MatrixXd> Xn_map(Xn.data(), n, p);
    Eigen::Map<Eigen::VectorXd> yn_map(yn.data(), n);

    TRexControlParameter n_ctrl;
    n_ctrl.scaling_mode = ScalingMode::NONE;
    TRexSelector n_sel(Xn_map, yn_map, 0.1, n_ctrl, 42, false);
    n_sel.select();
    std::cout << std::boolalpha << "same selection as under L2: "
              << (n_sel.getSelectedIndices() == sel.getSelectedIndices()) << '\n';
    //! [scaling]

    doc.begin("inplace");
    //! [inplace]
    const Eigen::MatrixXd X_before = X;
    const auto max_diff = [&] { return (X - X_before).cwiseAbs().maxCoeff(); };

    TRexSelector r_sel(X_map, y_map, 0.1, TRexControlParameter{}, 42, false);
    std::cout << std::scientific << std::setprecision(2);
    std::cout << "max |X - X_before| after construction: " << max_diff() << '\n';
    r_sel.restore();
    std::cout << "max |X - X_before| after restore():    " << max_diff() << '\n';
    //! [inplace]

    doc.begin("mc");
    //! [mc]
    const std::vector<double> targets = {0.1, 0.2, 0.3};
    const int MC = 20;
    std::ofstream mc_csv(std::string(DOC_OUTPUT_DIR) + "/trex_selector.figure_mc.csv");   // the figure
    for (std::size_t t = 0; t < targets.size(); ++t) {
        double fdr = 0.0, tpr = 0.0;
        for (int mc = 0; mc < MC; ++mc) {
            const int seed = 1000 * static_cast<int>(t + 1) + mc + 1;
            std::mt19937_64 gen_mc(static_cast<unsigned>(seed));
            std::normal_distribution<double> normal_mc(0.0, 1.0);
            const auto draw = [&] { return normal_mc(gen_mc); };
            Eigen::MatrixXd X_mc = Eigen::MatrixXd::NullaryExpr(n, p, draw);
            Eigen::VectorXd y_mc = X_mc * beta + Eigen::VectorXd::NullaryExpr(n, draw);
            Eigen::Map<Eigen::MatrixXd> X_mc_map(X_mc.data(), n, p);
            Eigen::Map<Eigen::VectorXd> y_mc_map(y_mc.data(), n);

            TRexControlParameter ctrl;
            ctrl.K = 10;
            TRexSelector s(X_mc_map, y_mc_map, targets[t], ctrl, seed, false);
            s.select();
            fdr += rates::compute_fdp(s.getSelectedIndices(), true_support) / MC;
            tpr += rates::compute_tpp(s.getSelectedIndices(), true_support) / MC;
        }
        std::cout << std::fixed << "tFDR " << std::setprecision(1) << targets[t]
                  << "   FDR " << std::setprecision(3) << fdr << "   TPR " << tpr << '\n';
        mc_csv << targets[t] << ',' << fdr << ',' << tpr << '\n';
    }
    //! [mc]

    return 0;
}
