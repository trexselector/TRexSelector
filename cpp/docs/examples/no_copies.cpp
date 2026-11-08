// Example program of the guide page "Working without data copies". The regions between
// the //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlasso_solver.hpp>
#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>
#include <ml_methods/model_selection/ridge_cv_svd.hpp>

using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::tsolvers::ad::SolverOptions;
using trex::tsolvers::ad::linear_model::lars_based::TLARS_Solver;
using trex::tsolvers::ad::linear_model::lars_based::TLASSO_Solver;
using trex::ml_methods::scaler_methods::LpNormScaler;
using trex::ml_methods::model_selection::ridge_cv_svd;
//! [includes]

int main() {
    DocSections doc("no_copies", DOC_OUTPUT_DIR);

    doc.begin("double");
    //! [double]
    std::mt19937_64 gen(1);
    std::uniform_int_distribution<int> digits(0, 2);
    Eigen::MatrixXi Xi = Eigen::MatrixXi::NullaryExpr(60, 20, [&] { return digits(gen); });

    // The library binds double storage. An integer matrix is cast once, by the caller,
    // and the map shares the storage of the cast
    Eigen::MatrixXd Xf = Xi.cast<double>();
    Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), Xf.rows(), Xf.cols());
    std::cout << std::boolalpha << (Xf_map.data() == Xf.data()) << '\n';
    //! [double]

    doc.begin("selector");
    //! [selector]
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 80, p = 40;
    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = 2.0 * X.col(0) + noise;
    const Eigen::MatrixXd X_before = X;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), n);
    TRexControlParameter ctrl;
    ctrl.K = 10;

    TRexSelector sel(X_map, y_map, 0.1, ctrl, /*seed=*/2, /*verbose=*/false);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << X.colwise().norm().head(3) << '\n';   // the working design: unit columns
    sel.select();
    std::cout << std::scientific << std::setprecision(2);
    std::cout << (X - X_before).cwiseAbs().maxCoeff() << '\n';   // the originals are back
    //! [selector]

    doc.begin("writers");
    //! [writers]
    Eigen::MatrixXd D = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), n, p);

    TLARS_Solver w1(X_map, D_map, y_map);   // a writer: normalize and intercept on
    try {
        TLASSO_Solver w2(X_map, D_map, y_map);   // a second writer on claimed matrices
    } catch (const std::runtime_error& e) {
        std::cout << "error: " << std::string(e.what()).substr(0, 90) << " ...\n";
    }
    w1.restore(X_map, D_map, y_map);
    //! [writers]

    doc.begin("readers");
    //! [readers]
    LpNormScaler sx(LpNormScaler::NormType::L2, /*center=*/true);
    sx.fit_transform_inplace(X_map);
    LpNormScaler sd(LpNormScaler::NormType::L2, true);
    sd.fit_transform_inplace(D_map);
    Eigen::VectorXd y_c = y.array() - y.mean();
    Eigen::Map<Eigen::VectorXd> yc_map(y_c.data(), n);

    SolverOptions opts;
    opts.normalize = false;
    opts.intercept = false;
    opts.precentered = true;
    TLARS_Solver r1(X_map, D_map, yc_map, opts);
    TLASSO_Solver r2(X_map, D_map, yc_map, opts);
    r1.executeStep(/*T_stop=*/2);
    r2.executeStep(2);

    sx.inverse_transform_inplace(X_map);
    sd.inverse_transform_inplace(D_map);
    std::cout << (X - X_before).cwiseAbs().maxCoeff() << '\n';
    //! [readers]

    doc.begin("low-memory");
    //! [low-memory]
    ridge_cv_svd cv_a(/*low_memory=*/false);
    cv_a.fit(X, y, /*n_folds=*/5, /*n_lambda=*/30);
    ridge_cv_svd cv_b(true);
    cv_b.fit(X, y, 5, 30);
    std::cout << cv_a.cv_mse().isApprox(cv_b.cv_mse()) << '\n';
    //! [low-memory]

    doc.begin("bytes");
    //! [bytes]
    std::cout << r1.memoryBytes() << '\n';
    //! [bytes]

    return 0;
}
