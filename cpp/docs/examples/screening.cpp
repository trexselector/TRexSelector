// Example program of the guide page "Screening". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_screening/trex_biobank_screening.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>

using trex::trex_selector_methods::ad::trex_biobank_screening::BiobankScreenTRex;
using trex::trex_selector_methods::ad::trex_biobank_screening::BiobankScreenTRexControl;
using trex::trex_selector_methods::ad::trex_screening::ScreenTRexControlParameter;
using trex::trex_selector_methods::ad::trex_screening::ScreenTRexMethod;
using trex::trex_selector_methods::ad::trex_screening::ScreenTRexSelector;
//! [includes]

//! [helpers]
void print_indices(const char* label, const std::vector<std::size_t>& indices) {
    std::cout << label;
    for (std::size_t j : indices) std::cout << ' ' << j;
    std::cout << '\n';
}
//! [helpers]

int main() {
    DocSections doc("screening", DOC_OUTPUT_DIR);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(2024);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 100, p = 300;
    const std::vector<std::size_t> true_support = {2, 49, 119, 199, 279};

    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (std::size_t j : true_support) beta[static_cast<Eigen::Index>(j)] = 2.0;
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    //! [data]

    doc.begin("screen");
    //! [screen]
    ScreenTRexControlParameter ctrl;
    ctrl.trex_ctrl.K = 20;

    ScreenTRexSelector scr(X_map, y_map, ctrl, /*seed=*/3, /*verbose=*/false);
    scr.select();
    const auto& res = scr.getScreenResult();
    print_indices("selected:     ", scr.getSelectedIndices());
    std::cout << "estimated FDR: " << res.estimated_FDR << '\n';
    //! [screen]

    doc.begin("mats");
    //! [mats]
    std::cout << "beta_mat:    " << res.beta_mat.size() << '\n';
    std::cout << "dummy_betas: " << res.dummy_betas.size() << '\n';
    //! [mats]

    doc.begin("boot");
    //! [boot]
    ScreenTRexControlParameter boot_ctrl = ctrl;
    boot_ctrl.use_bootstrap_CI = true;
    boot_ctrl.R_boot = 200;

    ScreenTRexSelector boot(X_map, y_map, boot_ctrl, 3, false);
    boot.select();
    const auto& boot_res = boot.getScreenResult();
    std::cout << std::boolalpha;
    std::cout << "used_bootstrap:   " << boot_res.used_bootstrap << '\n';
    std::cout << "confidence_level: " << boot_res.confidence_level << '\n';
    std::cout << "estimated FDR:    " << boot_res.estimated_FDR << '\n';
    print_indices("selected:        ", boot.getSelectedIndices());
    //! [boot]

    doc.begin("da");
    //! [da]
    ScreenTRexControlParameter ar_ctrl = ctrl;
    ar_ctrl.trex_method = ScreenTRexMethod::TREX_DA_AR1;

    ScreenTRexSelector ar(X_map, y_map, ar_ctrl, 3, false);
    ar.select();
    std::cout << std::fixed << std::setprecision(3) << "estimated correlation: "
              << ar.getScreenResult().estimated_correlation << '\n';
    //! [da]

    doc.begin("biobank");
    //! [biobank]
    const auto draw = [&] { return normal(gen); };
    Eigen::MatrixXd Y(n, 3);
    Y.col(0) = y;
    Y.col(1) = X * (0.5 * beta) + Eigen::VectorXd::NullaryExpr(n, draw);
    Y.col(2) = Eigen::VectorXd::NullaryExpr(n, draw);
    const std::vector<std::string> phenotypes = {"strong", "weak", "null"};
    Eigen::Map<Eigen::MatrixXd> Y_map(Y.data(), n, 3);

    BiobankScreenTRexControl bio_ctrl;
    bio_ctrl.lower_bound_FDR = 0.05;
    bio_ctrl.upper_bound_FDR = 0.15;
    bio_ctrl.target_FDR_trex = 0.1;
    bio_ctrl.trex_screen_ctrl = ctrl;

    BiobankScreenTRex bio(X_map, Y_map, bio_ctrl, 3, false);
    const std::vector records = bio.screenPhenotypes();
    std::cout << "records: " << records.size() << '\n';
    //! [biobank]

    doc.begin("records");
    //! [records]
    std::cout << std::left << std::setw(10) << "phenotype" << std::setw(9) << "selected"
              << std::setw(15) << "estimated_FDR" << std::setw(28) << "method"
              << "fallback" << '\n';
    for (const auto& r : records) {
        std::cout << std::left << std::setw(10) << phenotypes[r.phenotype_index]
                  << std::setw(9) << r.selected_indices.size()
                  << std::setw(15) << r.estimated_FDR
                  << std::setw(28) << r.method_used
                  << (r.used_fallback_trex ? "true" : "false") << '\n';
    }
    print_indices("selected for 'strong':", records[0].selected_indices);
    //! [records]

    return 0;
}
