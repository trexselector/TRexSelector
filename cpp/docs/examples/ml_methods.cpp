// Example program of the guide page "The ML methods". The regions between the
// //! [name] markers appear in the page, with the output they print next to them.
#include "doc_output.hpp"

//! [includes]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include <ml_methods/scaler_methods/z_score_scaler.hpp>
#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>
#include <ml_methods/ridge_regression/ridge.hpp>
#include <ml_methods/model_selection/ridge_cv_svd.hpp>
#include <ml_methods/model_selection/enet_cv_ccd.hpp>
#include <ml_methods/svd/svd.hpp>
#include <ml_methods/pca/pca.hpp>
#include <ml_methods/clustering/hierarchical/agglomerative/agglomerative_dispatcher.hpp>

using trex::ml_methods::scaler_methods::ZScoreScaler;
using trex::ml_methods::scaler_methods::LpNormScaler;
using trex::ml_methods::ridge::RidgeSolver;
using trex::ml_methods::model_selection::ridge_cv_svd;
using trex::ml_methods::model_selection::enet_gaussian;
using trex::ml_methods::model_selection::enet_cv_ccd;
using trex::ml_methods::svd::SVDSolver;
using trex::ml_methods::svd::SVDResult;
using trex::ml_methods::pca::PCA;
using trex::ml_methods::pca::PCAResult;
namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;
//! [includes]

int main() {
    DocSections doc("ml_methods", DOC_OUTPUT_DIR);
    std::cout << std::fixed << std::setprecision(3);

    doc.begin("data");
    //! [data]
    std::mt19937_64 gen(8);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Eigen::Index n = 120, p = 30;

    Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, p, [&] { return normal(gen); });
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta.head(5) << 2.0, -1.5, 0.0, 0.0, 1.0;
    Eigen::VectorXd noise = Eigen::VectorXd::NullaryExpr(n, [&] { return normal(gen); });
    Eigen::VectorXd y = X * beta + noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), n, p);   // the in-place methods take a map
    //! [data]

    doc.begin("zscore");
    //! [zscore]
    ZScoreScaler z;
    z.fit(X);
    std::cout << "centers: " << z.get_centers().head(3).transpose() << '\n';
    std::cout << "scales:  " << z.get_scales().head(3).transpose() << '\n';

    const Eigen::MatrixXd X_before = X;
    z.transform_inplace(X_map);
    Eigen::RowVectorXd means = X.colwise().mean();
    Eigen::ArrayXXd centered = (X.rowwise() - means).array();
    Eigen::RowVectorXd sds = (centered.square().colwise().sum() / (n - 1)).sqrt().matrix();
    std::cout << "largest |column mean| after: " << std::scientific
              << means.cwiseAbs().maxCoeff() << std::fixed << '\n';
    std::cout << "column sd after: " << sds.head(3) << '\n';

    z.inverse_transform_inplace(X_map);
    std::cout << "largest difference to the original: " << std::scientific
              << (X - X_before).cwiseAbs().maxCoeff() << std::fixed << '\n';
    //! [zscore]

    doc.begin("lpnorm");
    //! [lpnorm]
    LpNormScaler l2(LpNormScaler::NormType::L2, /*center=*/true);
    l2.fit_transform_inplace(X_map);
    std::cout << "column L2 norms: " << X.colwise().norm().head(3) << '\n';
    l2.inverse_transform_inplace(X_map);
    //! [lpnorm]

    doc.begin("ridge");
    //! [ridge]
    Eigen::VectorXd b_ridge = RidgeSolver::solve(X, y, /*lambda=*/1.0);
    std::cout << "ridge, first five coefficients: " << b_ridge.head(5).transpose() << '\n';

    ridge_cv_svd cv;
    cv.fit(X, y, /*n_folds=*/5, /*n_lambda=*/50);
    std::cout << std::setprecision(4) << "lambda_min: " << cv.cv_min()
              << "  lambda_1se: " << cv.cv_1se() << std::setprecision(3) << '\n';
    //! [ridge]

    doc.begin("ridge-curve");
    //! [ridge-curve]
    const Eigen::VectorXd& errors = cv.cv_mse();
    std::cout << "grid points: " << cv.lambdas().size() << '\n';
    std::cout << "CV error at lambda_min: " << errors(cv.index_min())
              << " (standard error " << cv.cv_sem()(cv.index_min()) << ")\n";
    std::cout << "CV error at lambda_1se: " << errors(cv.index_1se()) << '\n';
    std::cout << "largest CV error on the grid: " << errors.maxCoeff() << '\n';

    std::ofstream curve_csv(std::string(DOC_OUTPUT_DIR) + "/ml_methods.figure_ridge.csv");
    for (Eigen::Index i = 0; i < errors.size(); ++i) {   // the data of the figure: lambda, error
        curve_csv << cv.lambdas()[i] << ',' << errors[i] << '\n';
    }
    std::ofstream curve_meta(std::string(DOC_OUTPUT_DIR) + "/ml_methods.figure_ridge.meta.txt");
    curve_meta << "cv_min " << cv.cv_min() << "\ncv_1se " << cv.cv_1se() << '\n';
    //! [ridge-curve]

    doc.begin("enet");
    //! [enet]
    enet_gaussian en;
    en.fit(X, y, /*alpha=*/0.5, /*n_lambda=*/50);
    const Eigen::MatrixXd& B = en.coef();   // the path can stop before n_lambda points
    const Eigen::Index last = B.cols() - 1;
    std::cout << "coefficient path: " << B.rows() << " x " << B.cols() << '\n';
    std::cout << "grid points:      " << en.lambdas().size() << '\n';
    std::cout << std::boolalpha << "converged:        " << en.converged() << '\n';
    std::cout << "last point, first five coefficients: "
              << B.col(last).head(5).transpose() << '\n';
    std::cout << "predictions of three rows at the last point: "
              << en.predict(X.topRows(3)).col(last).transpose() << '\n';
    //! [enet]

    doc.begin("enet-cv");
    //! [enet-cv]
    enet_cv_ccd encv;
    encv.fit(X, y, /*alpha=*/0.5, /*n_folds=*/5, /*n_lambda=*/50);
    std::cout << std::setprecision(4) << "lambda_min: " << encv.cv_min()
              << "  lambda_1se: " << encv.cv_1se() << std::setprecision(3) << '\n';
    //! [enet-cv]

    doc.begin("svd");
    //! [svd]
    SVDSolver solver;
    SVDResult svd_res = solver.compute(X, /*M=*/3);
    std::cout << "U: " << svd_res.U.rows() << " x " << svd_res.U.cols()
              << "  V: " << svd_res.V.rows() << " x " << svd_res.V.cols() << '\n';
    std::cout << "S: " << svd_res.S.transpose() << '\n';
    //! [svd]

    doc.begin("pca");
    //! [pca]
    const Eigen::MatrixXd X_raw = X;
    PCA pca(X, /*M=*/3, /*center=*/true, /*normalize=*/true);
    PCAResult res = pca.fit();
    std::cout << "scores:   " << res.Z.rows() << " x " << res.Z.cols() << '\n';
    std::cout << "loadings: " << pca.getLoadings().rows() << " x "
              << pca.getLoadings().cols() << '\n';
    std::cout << "explained variance: " << pca.getExplainedVariance().transpose() << '\n';
    std::cout << "two rows transformed:\n" << pca.transform(X_raw.topRows(2)) << '\n';

    pca.restore(X);
    std::cout << "largest difference to the original: " << std::scientific
              << (X - X_raw).cwiseAbs().maxCoeff() << std::fixed << '\n';
    //! [pca]

    doc.begin("hac");
    //! [hac]
    // 60 observations in three groups of 20 with centers -3, 0 and 3, four features each
    std::mt19937_64 gen9(9);
    Eigen::MatrixXd Z(60, 4);
    for (Eigen::Index i = 0; i < 60; ++i) {
        const double center = (i < 20) ? -3.0 : (i < 40) ? 0.0 : 3.0;
        for (Eigen::Index j = 0; j < 4; ++j) Z(i, j) = center + normal(gen9);
    }

    // The engine clusters the columns of the matrix it is given
    const Eigen::MatrixXd Zt = Z.transpose();
    using Euclid = hac::DistancePolicy<Eigen::MatrixXd, hac::DistanceMetric::Euclidean>;
    constexpr auto ward = hac::LinkageMethod::Ward;
    std::vector<hac::MergeStep> merges =
        hac::AgglomerativeClustering::cluster<Eigen::MatrixXd, Euclid, ward>(
            Zt, /*use_mmap=*/false, /*verbose=*/false);
    std::cout << "merge steps: " << merges.size() << '\n';

    std::vector<Eigen::Index> labels = hac::DendrogramUtils::cut_tree(merges, 60, 3);
    Eigen::Matrix3i table = Eigen::Matrix3i::Zero();   // found clusters x true groups
    for (Eigen::Index i = 0; i < 60; ++i) table(labels[i], i / 20) += 1;
    std::cout << "found clusters (rows) against the true groups (columns):\n"
              << table << '\n';
    //! [hac]

    doc.begin("hac-vars");
    //! [hac-vars]
    // The design itself clusters its columns, the variables, by correlation distance
    using Corr = hac::DistancePolicy<Eigen::MatrixXd, hac::DistanceMetric::Correlation>;
    constexpr auto single = hac::LinkageMethod::Single;
    std::vector<hac::MergeStep> link_vars =
        hac::AgglomerativeClustering::cluster<Eigen::MatrixXd, Corr, single>(X, false, false);
    std::vector<Eigen::Index> var_labels = hac::DendrogramUtils::cut_tree(link_vars, p, 5);
    std::cout << "first six variables:";
    for (Eigen::Index j = 0; j < 6; ++j) std::cout << ' ' << var_labels[j];
    std::cout << '\n';
    //! [hac-vars]

    return 0;
}
