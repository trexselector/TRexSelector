// ========================================================================================
// test_trex_tikhonov.cpp
// ========================================================================================
/**
 * @file test_trex_tikhonov.cpp
 *
 * @brief Unit tests for TRexTikhonovSelector (general-Tikhonov T-Rex):
 * construction validation, the collapse onto the lasso base run at
 * lambda_2 = 0, the cluster-aware dummy path (gates, layer correlation,
 * prior-groups end-to-end), planted-group recovery with an informed K, and
 * the CV lambda_2 path.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// std includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// project trex includes
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_cluster_dummies.hpp>

// project utils includes
#include <utils/datageneration/utils_datagen.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::trex_selector_methods::ad::trex_tikhonov {

// Namespace aliases for convenience
using namespace trex::trex_selector_methods::ad::trex_tikhonov;
using namespace trex::trex_selector_methods::ad::trex_core;
using namespace trex::utils::datageneration::datagen;
namespace cdum = trex::trex_selector_methods::ad::utils::cluster_dummies;
namespace gvs  = trex::trex_selector_methods::ad::trex_gvs;

// ========================================================================================

namespace {

/** @brief Identity K (p x p). */
Eigen::SparseMatrix<double> identityK(Eigen::Index p) {
    Eigen::SparseMatrix<double> K(p, p);
    K.setIdentity();
    return K;
}

/** @brief Group-mean IEN matrix W = sum_m 1_m 1_m^T / p_m for contiguous
 *  groups of size `gsize`. */
Eigen::SparseMatrix<double> groupMeanK(Eigen::Index p, Eigen::Index gsize) {
    Eigen::SparseMatrix<double> W(p, p);
    std::vector<Eigen::Triplet<double>> trip;
    for (Eigen::Index i = 0; i < p; ++i)
        for (Eigen::Index j = 0; j < p; ++j)
            if (i / gsize == j / gsize)
                trip.emplace_back(i, j, 1.0 / static_cast<double>(gsize));
    W.setFromTriplets(trip.begin(), trip.end());
    return W;
}

/** @brief Correlated grouped design with the planted actives in the first
 *  two groups (mirrors the GVS test generator: within-group AR structure,
 *  additive signal on cols 0-5). */
struct GroupedData {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;

    GroupedData(Eigen::Index n, Eigen::Index p, double rho, double sigma,
                unsigned seed) {
        std::mt19937 rng(seed);
        std::normal_distribution<double> N01(0.0, 1.0);

        X.resize(n, p);
        for (Eigen::Index g = 0; g < p; g += 5) {
            Eigen::VectorXd z(n);
            for (Eigen::Index i = 0; i < n; ++i) z(i) = N01(rng);
            for (Eigen::Index j = g; j < std::min(g + 5, p); ++j) {
                for (Eigen::Index i = 0; i < n; ++i) {
                    X(i, j) = rho * z(i) + std::sqrt(1.0 - rho * rho) * N01(rng);
                }
            }
        }

        y.resize(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            double signal = 0.0;
            for (Eigen::Index j = 0; j < 3; ++j) { signal += X(i, j); }
            for (Eigen::Index j = 3; j < 6; ++j) { signal -= X(i, j); }
            y(i) = signal + sigma * N01(rng);
        }
    }
};

} // namespace

// ========================================================================================

class TRexTikhonovTest : public ::testing::Test {
protected:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::Map<Eigen::MatrixXd> X_map;
    Eigen::Map<Eigen::VectorXd> y_map;

    TRexTikhonovTest()
        : X(Eigen::MatrixXd::Random(30, 10)),
          y(Eigen::VectorXd::Random(30)),
          X_map(X.data(), X.rows(), X.cols()),
          y_map(y.data(), y.size()) {}
};


// ========================================================================================
// Construction validation
// ========================================================================================

/** @brief An empty or mis-dimensioned K, or a solver_type other than
 *  TIKNET, must throw at construction. */
TEST_F(TRexTikhonovTest, Validation_ThrowsOnBadKOrSolver) {

    // Default control: K is empty.
    TRexTikhonovControlParameter ctrl;
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                 std::invalid_argument);

    // Wrong dimensions (p = 10 here).
    ctrl.tikhonov_K = identityK(11);
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                 std::invalid_argument);

    // Correct K, wrong solver.
    ctrl.tikhonov_K = identityK(10);
    ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TLARS;
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                 std::invalid_argument);

    // Correct K: the default solver, TIKNET, constructs.
    EXPECT_EQ(TRexTikhonovControlParameter().trex_ctrl.solver_type, disp::SolverTypeForTRex::TIKNET);
    ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TIKNET;
    EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl));

    // The coordinate-descent solver is no longer accepted, and the message says so.
    const Eigen::MatrixXd X_before = X;
    ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDTIKNET;
    try {
        TRexTikhonovSelector refused(X_map, y_map, 0.1, ctrl);
        FAIL() << "TCDTIKNET was accepted";
    } catch (const std::invalid_argument& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("solver_type must be TIKNET"), std::string::npos) << what;
        EXPECT_NE(what.find("TCDTIKNET is no longer accepted"), std::string::npos) << what;
    }
    EXPECT_LT((X - X_before).cwiseAbs().maxCoeff(), 1e-12);   // a refused control leaves X alone
}


/** @brief The content of K is checked at construction, under the name of
 *  the selector: finite entries, symmetry, positive semi-definiteness. */
TEST_F(TRexTikhonovTest, Validation_ContentOfKIsCheckedAtConstruction) {
    const auto expect_refused = [&](const Eigen::SparseMatrix<double>& K) {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = K;
        ctrl.lambda_2   = 1.0;
        try {
            TRexTikhonovSelector sel(X_map, y_map, 0.1, ctrl, 1, false);
            ADD_FAILURE() << "the matrix was accepted";
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("TRexTikhonovSelector"), std::string::npos) << e.what();
        }
    };
    Eigen::SparseMatrix<double> K = identityK(10);
    K.coeffRef(0, 1) = 0.5;
    expect_refused(K);                                   // not symmetric
    K = identityK(10);
    K.coeffRef(9, 9) = -1.0;
    expect_refused(K);                                   // not positive semi-definite
    K = identityK(10);
    K.coeffRef(2, 2) = std::numeric_limits<double>::quiet_NaN();
    expect_refused(K);                                   // not finite

    // The three standard penalties pass.
    const std::vector<Eigen::Index> halves = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    for (const Eigen::SparseMatrix<double>& good :
         {identityK(10), groupMeanK(10, 5), TRexTikhonovSelector::buildGroupCenteringK(halves)}) {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = good;
        ctrl.lambda_2   = 1.0;
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl, 1, false));
    }
}


/** @brief lambda_2 is a finite number: a value >= 0 or the negative
 *  sentinel of the rules. NaN and infinity are refused at construction. */
TEST_F(TRexTikhonovTest, Validation_Lambda2MustBeFinite) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const double bad : {nan, inf, -inf}) {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = identityK(10);
        ctrl.lambda_2   = bad;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl, 1, false), std::invalid_argument)
            << "lambda_2 = " << bad;
    }
    for (const double good : {-1.0, 0.0, 2.5}) {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = identityK(10);
        ctrl.lambda_2   = good;
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl, 1, false)) << "lambda_2 = " << good;
    }
}


/** @brief The linkage and the sizes of the cross-validation are checked at
 *  construction when they will be used, and left alone when they will not. */
TEST_F(TRexTikhonovTest, Validation_LinkageAndCvSizesAreCheckedAtConstruction) {
    namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;
    const auto base = [] {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = identityK(10);
        ctrl.lambda_2   = 1.0;
        return ctrl;
    };
    const auto refused = [&](const TRexTikhonovControlParameter& ctrl) {
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl, 1, false), std::invalid_argument);
    };
    const auto accepted = [&](const TRexTikhonovControlParameter& ctrl) {
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl, 1, false));
    };

    // The linkage: used when the selector clusters the variables itself.
    for (const auto linkage : {hac::LinkageMethod::Ward, hac::LinkageMethod::Median,
                               hac::LinkageMethod::Centroid}) {
        auto ctrl = base();
        ctrl.hc_linkage = linkage;
        refused(ctrl);
        ctrl.prior_groups = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};   // supplied groups: no clustering
        accepted(ctrl);
        ctrl.prior_groups.clear();
        ctrl.cluster_dummies = false;                         // i.i.d. dummies: no clustering
        accepted(ctrl);
    }
    for (const auto linkage : {hac::LinkageMethod::Single, hac::LinkageMethod::Complete,
                               hac::LinkageMethod::Average, hac::LinkageMethod::WPGMA}) {
        auto ctrl = base();
        ctrl.hc_linkage = linkage;
        accepted(ctrl);
    }

    // The sizes of the cross-validation: used when a CV rule chooses lambda_2.
    auto cv = base();
    cv.lambda_2 = -1.0;
    cv.lambda2_method = TikLambda2Method::CV_1SE_IEN_CCD;
    accepted(cv);
    for (const int folds : {0, 1, 31}) {                      // n = 30
        auto ctrl = cv;
        ctrl.cv_n_folds = folds;
        refused(ctrl);
    }
    { auto ctrl = cv; ctrl.cv_n_lambda2 = 0; refused(ctrl); }
    { auto ctrl = cv; ctrl.cv_n_lambda1 = 0; refused(ctrl); }
    { auto ctrl = cv; ctrl.cv_n_lambda = 0; accepted(ctrl); }  // the SVD grid is not used
    cv.lambda2_method = TikLambda2Method::CV_MIN_TIK_SVD;
    { auto ctrl = cv; ctrl.cv_n_lambda = 0; refused(ctrl); }
    { auto ctrl = cv; ctrl.cv_n_lambda2 = 0; accepted(ctrl); } // the IEN grid is not used

    // A fixed lambda_2 and the pilot rule use none of them.
    for (const double lambda2 : {1.0, -1.0}) {
        auto ctrl = base();
        ctrl.lambda_2 = lambda2;                              // -1 with the default PILOT_KNEE
        ctrl.cv_n_folds = 0;
        ctrl.cv_n_lambda = 0;
        accepted(ctrl);
    }
}


/** @brief Under GROUP_CENTERING a supplied K has to be the centering matrix
 *  of the groups. Another matrix would let the cross-validation tune one
 *  penalty while the solver applies a second. */
TEST_F(TRexTikhonovTest, GroupCentering_SuppliedMatrixMustBeTheCenteringMatrix) {
    const std::vector<Eigen::Index> groups = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    const std::vector<Eigen::Index> other  = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
    const auto control = [&](const Eigen::SparseMatrix<double>& K) {
        TRexTikhonovControlParameter ctrl;
        ctrl.penalty_structure = TikPenaltyStructure::GROUP_CENTERING;
        ctrl.prior_groups      = groups;
        ctrl.tikhonov_K        = K;
        ctrl.lambda_2          = 1.0;
        ctrl.trex_ctrl.K       = 5;
        return ctrl;
    };
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, control(identityK(10)), 1, false),
                 std::invalid_argument);
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1,
                                      control(TRexTikhonovSelector::buildGroupCenteringK(other)), 1, false),
                 std::invalid_argument);
    EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, control(identityK(11)), 1, false),
                 std::invalid_argument);

    // The right matrix and no matrix give the same run.
    const auto run = [&](const TRexTikhonovControlParameter& ctrl) {
        Eigen::MatrixXd Xc = X;
        Eigen::VectorXd yc = y;
        Eigen::Map<Eigen::MatrixXd> Xm(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::VectorXd> ym(yc.data(), yc.size());
        TRexTikhonovSelector sel(Xm, ym, 0.3, ctrl, 7, false);
        sel.select();
        return sel.getSelectedIndices();
    };
    EXPECT_EQ(run(control(TRexTikhonovSelector::buildGroupCenteringK(groups))),
              run(control(Eigen::SparseMatrix<double>())));
}


/** @brief With cluster dummies the four strategies of the GVS family run.
 *  The permutations are refused, and run with i.i.d. dummies. */
TEST_F(TRexTikhonovTest, LLoopStrategies_WithAndWithoutClusterDummies) {
    const GroupedData d(/*n=*/60, /*p=*/20, /*rho=*/0.8, /*sigma=*/1.0, /*seed=*/9);
    const auto run = [&](LLoopStrategy strategy, bool cluster_dummies) {
        Eigen::MatrixXd Xc = d.X;
        Eigen::VectorXd yc = d.y;
        Eigen::Map<Eigen::MatrixXd> Xm(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::VectorXd> ym(yc.data(), yc.size());
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K      = identityK(Xc.cols());
        ctrl.lambda_2        = 1.0;
        ctrl.cluster_dummies = cluster_dummies;
        ctrl.trex_ctrl.K                    = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 3;
        ctrl.trex_ctrl.lloop_strategy       = strategy;
        TRexTikhonovSelector sel(Xm, ym, 0.2, ctrl, 4, false);
        sel.select();
        return sel.getSelectedIndices();
    };
    for (const LLoopStrategy strategy : {LLoopStrategy::STANDARD, LLoopStrategy::HCONCAT,
                                         LLoopStrategy::SKIPL, LLoopStrategy::SEEDED}) {
        EXPECT_NO_THROW(run(strategy, /*cluster_dummies=*/true));
    }
    EXPECT_EQ(run(LLoopStrategy::SEEDED, true), run(LLoopStrategy::STANDARD, true));
    for (const LLoopStrategy strategy : {LLoopStrategy::PERMUTATION, LLoopStrategy::PERMUTATION_SEEDED}) {
        EXPECT_THROW(run(strategy, /*cluster_dummies=*/true), std::invalid_argument);
        EXPECT_NO_THROW(run(strategy, /*cluster_dummies=*/false));
    }
}


/** @brief getDummyClusters() reports the clusters the dummies were drawn
 *  for: the supplied groups, or the clusters of the correlation clustering,
 *  with their number and their origin. */
TEST_F(TRexTikhonovTest, DummyClusters_AreReported) {
    namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;
    const GroupedData d(/*n=*/80, /*p=*/20, /*rho=*/0.9, /*sigma=*/1.0, /*seed=*/5);
    const Eigen::Index p = d.X.cols();
    const auto control = [&] {
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K  = identityK(p);
        ctrl.lambda_2    = 1.0;
        ctrl.trex_ctrl.K = 5;
        return ctrl;
    };
    const auto run = [&](const TRexTikhonovControlParameter& ctrl, bool select) {
        Eigen::MatrixXd X = d.X;
        Eigen::VectorXd y = d.y;
        Eigen::Map<Eigen::MatrixXd> Xm(X.data(), X.rows(), X.cols());
        Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
        TRexTikhonovSelector sel(Xm, ym, 0.2, ctrl, 3, false);
        if (select) { sel.select(); }
        return sel.getDummyClusters();
    };

    // Before a run the record is empty.
    const DummyClusters before = run(control(), /*select=*/false);
    EXPECT_EQ(before.groups_vec.size(), 0);
    EXPECT_EQ(before.max_clusters, 0U);
    EXPECT_TRUE(before.hc_method_used.empty());

    // Supplied groups are reported as they are.
    auto with_groups = control();
    for (Eigen::Index j = 0; j < p; ++j) { with_groups.prior_groups.push_back(j / 5); }
    const DummyClusters prior = run(with_groups, true);
    ASSERT_EQ(prior.groups_vec.size(), p);
    for (Eigen::Index j = 0; j < p; ++j) { EXPECT_EQ(prior.groups_vec(j), j / 5); }
    EXPECT_EQ(prior.max_clusters, 4U);
    EXPECT_EQ(prior.hc_method_used, "prior");

    // Found clusters: a partition into max_clusters non-empty clusters. With
    // single linkage no two variables of different clusters correlate above
    // corr_max, which is what the cut at the height 1 - corr_max means.
    Eigen::MatrixXd Xc = d.X.rowwise() - d.X.colwise().mean();
    Xc.colwise().normalize();
    const Eigen::MatrixXd corr = Xc.transpose() * Xc;
    for (const auto linkage : {hac::LinkageMethod::Single, hac::LinkageMethod::Complete}) {
        auto found_ctrl = control();
        found_ctrl.hc_linkage = linkage;
        found_ctrl.corr_max   = 0.5;
        const DummyClusters found = run(found_ctrl, true);
        ASSERT_EQ(found.groups_vec.size(), p);
        ASSERT_GE(found.max_clusters, 1U);
        std::vector<std::size_t> size(found.max_clusters, 0);
        for (Eigen::Index j = 0; j < p; ++j) {
            ASSERT_GE(found.groups_vec(j), 0);
            ASSERT_LT(static_cast<std::size_t>(found.groups_vec(j)), found.max_clusters);
            ++size[static_cast<std::size_t>(found.groups_vec(j))];
        }
        for (const std::size_t s : size) { EXPECT_GE(s, 1U); }
        if (linkage == hac::LinkageMethod::Single) {
            EXPECT_EQ(found.hc_method_used, "single");
            for (Eigen::Index a = 0; a < p; ++a) {
                for (Eigen::Index b = a + 1; b < p; ++b) {
                    if (found.groups_vec(a) != found.groups_vec(b)) {
                        EXPECT_LE(std::abs(corr(a, b)), 0.5 + 1e-9) << "variables " << a << " and " << b;
                    }
                }
            }
        } else {
            EXPECT_EQ(found.hc_method_used, "complete");
        }
    }

    // Without cluster dummies there are no clusters.
    auto iid = control();
    iid.cluster_dummies = false;
    const DummyClusters none = run(iid, true);
    EXPECT_EQ(none.groups_vec.size(), 0);
    EXPECT_EQ(none.max_clusters, 0U);
    EXPECT_TRUE(none.hc_method_used.empty());
}


/** @brief gammaToK builds K = Gamma^T Gamma (first-difference operator ->
 *  tridiagonal graph Laplacian). */
TEST_F(TRexTikhonovTest, GammaToK_FirstDifferenceGivesLaplacian) {
    const Eigen::Index p = 6;
    Eigen::SparseMatrix<double> gamma(p - 1, p);
    std::vector<Eigen::Triplet<double>> trip;
    for (Eigen::Index i = 0; i < p - 1; ++i) {
        trip.emplace_back(i, i, -1.0);
        trip.emplace_back(i, i + 1, 1.0);
    }
    gamma.setFromTriplets(trip.begin(), trip.end());

    const Eigen::MatrixXd K =
        Eigen::MatrixXd(TRexTikhonovSelector::gammaToK(gamma));
    const Eigen::MatrixXd K_ref =
        Eigen::MatrixXd(gamma.transpose()) * Eigen::MatrixXd(gamma);
    EXPECT_TRUE(K.isApprox(K_ref, 1e-14));
    EXPECT_DOUBLE_EQ(K(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(K(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(K(0, 1), -1.0);
}


// ========================================================================================
// lambda_2 = 0: the lasso
// ========================================================================================

/** @brief With lambda_2 = 0 the penalty vanishes and the solver runs its
 *  lasso path. With i.i.d. dummies (`cluster_dummies = false`) the selector
 *  then works on the dummies of the plain TRexSelector, so it selects what
 *  that selector selects with TLASSO. For lambda_2 > 0 there is no such
 *  collapse onto the elastic net, also not for K = I: in the FOLDED coupling
 *  a dummy shares the penalty of its variable. */
TEST_F(TRexTikhonovTest, Equivalence_ZeroLambda2MatchesLassoBaseRun) {
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};

    auto run_tik = [&]() {
        SyntheticData data(150, 40, support, coefs, 1.0, 42);
        auto Xd = data.getX();
        auto yd = data.getY();
        Eigen::Map<Eigen::MatrixXd> Xm(Xd.data(), Xd.rows(), Xd.cols());
        Eigen::Map<Eigen::VectorXd> ym(yd.data(), yd.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = identityK(40);
        ctrl.lambda_2   = 0.0;
        ctrl.cluster_dummies = false;   // i.i.d. dummies, as the base selector draws them
        ctrl.trex_ctrl.K = 20;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        TRexTikhonovSelector trex(Xm, ym, 0.3, ctrl, 42, false);
        return trex.select().selected_var;
    };

    auto run_lasso = [&]() {
        SyntheticData data(150, 40, support, coefs, 1.0, 42);
        auto Xd = data.getX();
        auto yd = data.getY();
        Eigen::Map<Eigen::MatrixXd> Xm(Xd.data(), Xd.rows(), Xd.cols());
        Eigen::Map<Eigen::VectorXd> ym(yd.data(), yd.size());

        TRexControlParameter ctrl;
        ctrl.solver_type = disp::SolverTypeForTRex::TLASSO;
        ctrl.K = 20;
        ctrl.max_dummy_multiplier = 2;
        TRexSelector trex(Xm, ym, 0.3, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const Eigen::VectorXi sel_tik   = run_tik();
    const Eigen::VectorXi sel_lasso = run_lasso();

    ASSERT_EQ(sel_tik.size(), sel_lasso.size());
    EXPECT_GT(sel_tik.sum(), 0) << "the run selected nothing";
    EXPECT_TRUE((sel_tik.array() == sel_lasso.array()).all())
        << "The Tikhonov run at lambda_2 = 0 diverged from the TLASSO base run.\n"
        << "  Tikhonov n_selected = " << sel_tik.sum() << "\n"
        << "  TLASSO   n_selected = " << sel_lasso.sum();
}


// ========================================================================================
// The default run on a chain Laplacian
// ========================================================================================

/** @brief The default control (homotopy solver, pilot rule over its default
 *  grid up to 3000) on an unscaled chain Laplacian selects the planted block
 *  (here all ten actives and two nulls), and the pick stays below the top
 *  of the grid. */
TEST_F(TRexTikhonovTest, DefaultRun_ChainLaplacianSelectsThePlantedBlock) {
    const Eigen::Index n = 120, p = 60;
    // Box-Muller on the raw engine output: the same data on every platform.
    std::mt19937 rng(2718);
    auto N01 = [&rng]() {
        const double u1 = (static_cast<double>(rng()) + 1.0) / 4294967297.0;
        const double u2 = (static_cast<double>(rng()) + 1.0) / 4294967297.0;
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    };
    Eigen::MatrixXd X(n, p);
    for (Eigen::Index j = 0; j < p; ++j) {
        for (Eigen::Index i = 0; i < n; ++i) { X(i, j) = N01(); }
    }
    Eigen::VectorXd y(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        y(i) = X.row(i).segment(20, 10).sum() + N01();   // ten equal actives, columns 20 to 29
    }
    Eigen::SparseMatrix<double> gamma(p - 1, p);
    std::vector<Eigen::Triplet<double>> trip;
    for (Eigen::Index i = 0; i < p - 1; ++i) {
        trip.emplace_back(i, i, -1.0);
        trip.emplace_back(i, i + 1, 1.0);
    }
    gamma.setFromTriplets(trip.begin(), trip.end());

    Eigen::Map<Eigen::MatrixXd> Xm(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> ym(y.data(), n);
    TRexTikhonovControlParameter ctrl;                 // lambda_2 < 0: the pilot rule
    ctrl.tikhonov_K = TRexTikhonovSelector::gammaToK(gamma);
    ctrl.trex_ctrl.K = 10;
    TRexTikhonovSelector trex(Xm, ym, 0.1, ctrl, 4, false);
    const Eigen::VectorXi selected = trex.select().selected_var;

    const int true_positives = selected.segment(20, 10).sum();
    const int false_positives = selected.sum() - true_positives;
    EXPECT_GE(true_positives, 8);
    EXPECT_LE(false_positives, 4);
    EXPECT_TRUE(trex.getLambda2PilotResult().ran);
    EXPECT_LT(trex.getLambda2Used(), 1000.0) << "the pilot picked the top of its grid";
}


// ========================================================================================
// End-to-end recovery with an informed K
// ========================================================================================

/** @brief Group-mean K on a correlated grouped design: the planted actives
 *  (cols 0-5) are recovered with few false selections, and the fixed
 *  lambda_2 is reported back verbatim. */
TEST_F(TRexTikhonovTest, EndToEnd_GroupMeanKRecoversPlanted) {
    const Eigen::Index n = 150, p = 40;
    GroupedData d(n, p, 0.75, 1.0, 7);
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexTikhonovControlParameter ctrl;
    ctrl.tikhonov_K = groupMeanK(p, 5);
    ctrl.lambda_2   = 1.0;
    ctrl.trex_ctrl.K = 5;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
    auto result = trex.select();

    ASSERT_EQ(result.selected_var.size(), p);
    EXPECT_DOUBLE_EQ(trex.getLambda2Used(), 1.0);
    for (Eigen::Index j = 0; j < result.Phi_prime.size(); ++j) {
        EXPECT_LE(result.Phi_prime(j), 1.0);
    }

    int hits = 0, false_sel = 0;
    for (Eigen::Index j = 0; j < p; ++j) {
        if (result.selected_var(j) == 1) {
            if (j < 6) { ++hits; } else { ++false_sel; }
        }
    }
    EXPECT_GE(hits, 4)
        << "Group-mean K failed to recover the planted actives.";
    EXPECT_LE(false_sel, 4)
        << "Group-mean K selected too many null variables.";
}


// ========================================================================================
// Auto lambda_2 (sparse-K CV tuners)
// ========================================================================================

/** @brief Analytic Tikhonov CV backbone on an n > p design: resolves a
 *  positive lambda_2 and completes selection. */
TEST_F(TRexTikhonovTest, AutoLambda2_TikSvdResolvesAndRuns) {
    const Eigen::Index n = 150, p = 40;
    GroupedData d(n, p, 0.75, 1.0, 7);
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexTikhonovControlParameter ctrl;
    ctrl.tikhonov_K     = groupMeanK(p, 5);
    ctrl.lambda_2       = -1.0;   // sentinel -> CV
    ctrl.lambda2_method = TikLambda2Method::CV_1SE_TIK_SVD;
    ctrl.cv_n_folds     = 5;
    ctrl.cv_n_lambda    = 50;
    ctrl.trex_ctrl.K = 3;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
    auto result = trex.select();

    ASSERT_EQ(result.selected_var.size(), p);
    EXPECT_GT(trex.getLambda2Used(), 0.0);
}


// ========================================================================================
// Cluster-aware dummies (default path, GVS parity)
// ========================================================================================

/** @brief Configurations the cluster-aware dummy path cannot serve must be
 *  rejected at construction: the PERMUTATION-family strategies (K
 *  row-permuted copies of one base cannot carry the per-cluster MVN
 *  convention), malformed prior_groups, and an out-of-range corr_max.
 *  Everything else — including SEEDED and memory mapping — is served
 *  natively by the DummyGenerator cluster mode. */
TEST_F(TRexTikhonovTest, ClusterDummies_GateThrowsOnUnsupportedConfig) {

    TRexTikhonovControlParameter base;
    base.tikhonov_K = identityK(10);

    {   // PERMUTATION cannot carry cluster-MVN dummies.
        auto ctrl = base;
        ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::PERMUTATION;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                     std::invalid_argument);
    }
    {   // PERMUTATION_SEEDED likewise.
        auto ctrl = base;
        ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::PERMUTATION_SEEDED;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                     std::invalid_argument);
    }
    {   // SEEDED and memory mapping are supported.
        auto ctrl = base;
        ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::SEEDED;
        ctrl.trex_ctrl.use_memory_mapping = true;
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl));
    }
    {   // prior_groups length mismatch (p = 10 here).
        auto ctrl = base;
        ctrl.prior_groups.assign(9, 0);
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                     std::invalid_argument);
    }
    {   // Non-contiguous cluster IDs (ID 1 missing).
        auto ctrl = base;
        ctrl.prior_groups.assign(10, 0);
        ctrl.prior_groups[5] = 2;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                     std::invalid_argument);
    }
    {   // corr_max out of range (HAC route).
        auto ctrl = base;
        ctrl.corr_max = 1.5;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl),
                     std::invalid_argument);
    }
    {   // With cluster_dummies off, the gated knobs are inert.
        auto ctrl = base;
        ctrl.cluster_dummies = false;
        ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::SEEDED;
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl));
    }
}


/** @brief The shared layer draw reproduces the per-cluster correlation of
 *  the working-scale X: strong within-cluster, near-zero across clusters. */
TEST_F(TRexTikhonovTest, ClusterDummies_LayerCarriesClusterCorrelation) {
    const Eigen::Index n = 200, p = 20;
    GroupedData d(n, p, 0.9, 1.0, 3);   // within-group corr = 0.81

    // Selector working scale: centered + unit-L2 columns.
    Eigen::MatrixXd Xs = d.X;
    Xs.rowwise() -= Xs.colwise().mean();
    for (Eigen::Index j = 0; j < p; ++j) { Xs.col(j) /= Xs.col(j).norm(); }
    Eigen::Map<Eigen::MatrixXd> Xm(Xs.data(), n, p);

    std::vector<std::vector<Eigen::Index>> clusters;
    for (Eigen::Index g = 0; g < p; g += 5) {
        std::vector<Eigen::Index> c;
        for (Eigen::Index j = g; j < g + 5; ++j) { c.push_back(j); }
        clusters.push_back(std::move(c));
    }

    const auto chol = cdum::buildClusterCholeskys(Xm, clusters, n, "test");
    std::mt19937 rng(123);
    Eigen::MatrixXd L = cdum::drawClusterDummyLayer(
        static_cast<std::size_t>(n), static_cast<std::size_t>(p),
        clusters, chol, rng);
    L.rowwise() -= L.colwise().mean();

    auto corr = [&L](Eigen::Index a, Eigen::Index b) {
        return L.col(a).dot(L.col(b)) / (L.col(a).norm() * L.col(b).norm());
    };

    double min_within = 1.0;
    double max_cross  = 0.0;
    for (Eigen::Index a = 0; a < p; ++a) {
        for (Eigen::Index b = a + 1; b < p; ++b) {
            const double c = std::abs(corr(a, b));
            if (a / 5 == b / 5) {
                min_within = std::min(min_within, c);
            } else {
                max_cross = std::max(max_cross, c);
            }
        }
    }
    EXPECT_GT(min_within, 0.6)
        << "Within-cluster dummy correlation too weak (design corr 0.81).";
    EXPECT_LT(max_cross, 0.35)
        << "Cross-cluster dummy correlation too strong (expected ~0).";
}


/** @brief Under cluster dummies the layer draws depend only on
 *  (seed, k, L): SEEDED (one D_k at a time, serial) produces content
 *  bit-identical to STANDARD at the same L, so both runs must select
 *  identically. */
TEST_F(TRexTikhonovTest, ClusterDummies_SeededMatchesStandard) {
    const Eigen::Index n = 150, p = 40;

    auto run_with = [&](LLoopStrategy strategy) {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = groupMeanK(p, 5);
        ctrl.lambda_2   = 1.0;
        ctrl.prior_groups.resize(static_cast<std::size_t>(p));
        for (Eigen::Index j = 0; j < p; ++j) {
            ctrl.prior_groups[static_cast<std::size_t>(j)] = j / 5;
        }
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        ctrl.trex_ctrl.lloop_strategy = strategy;

        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const Eigen::VectorXi sel_std = run_with(LLoopStrategy::STANDARD);
    const Eigen::VectorXi sel_od  = run_with(LLoopStrategy::SEEDED);

    ASSERT_EQ(sel_std.size(), sel_od.size());
    EXPECT_TRUE((sel_std.array() == sel_od.array()).all())
        << "SEEDED diverged from STANDARD under cluster dummies.\n"
        << "  STANDARD n_selected = " << sel_std.sum() << "\n"
        << "  SEEDED n_selected = " << sel_od.sum();
}


/** @brief The memory-mapped writers regenerate the cluster layers
 *  bit-exactly from the deterministic streams, so for every supported
 *  strategy the memory-mapped run must select identically to the
 *  in-memory run. */
TEST_F(TRexTikhonovTest, ClusterDummies_MemmapMatchesInMemory) {
    const Eigen::Index n = 150, p = 40;

    auto run_with = [&](LLoopStrategy strategy, bool use_mmap) {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K = groupMeanK(p, 5);
        ctrl.lambda_2   = 1.0;
        ctrl.prior_groups.resize(static_cast<std::size_t>(p));
        for (Eigen::Index j = 0; j < p; ++j) {
            ctrl.prior_groups[static_cast<std::size_t>(j)] = j / 5;
        }
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        ctrl.trex_ctrl.lloop_strategy = strategy;
        ctrl.trex_ctrl.use_memory_mapping = use_mmap;

        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const struct { LLoopStrategy s; const char* name; } strategies[] = {
        {LLoopStrategy::STANDARD, "STANDARD"},
        {LLoopStrategy::HCONCAT,  "HCONCAT"},
        {LLoopStrategy::SEEDED, "SEEDED"},
    };
    for (const auto& st : strategies) {
        const Eigen::VectorXi sel_mem  = run_with(st.s, false);
        const Eigen::VectorXi sel_mmap = run_with(st.s, true);
        ASSERT_EQ(sel_mem.size(), sel_mmap.size()) << st.name;
        EXPECT_TRUE((sel_mem.array() == sel_mmap.array()).all())
            << st.name << ": memory-mapped selection differs from the "
            << "in-memory selection.\n"
            << "  in-memory n_selected = " << sel_mem.sum() << "\n"
            << "  memmap    n_selected = " << sel_mmap.sum();
    }
}


/** @brief FOLDED dummy coupling with the group-mean K must reproduce
 *  TRexGVSSelector(IEN): identical prior groups, lambda_2 and seed give
 *  bit-identical cluster-dummy draws (shared makeLayerStream) and — with
 *  the folded K-coupling being the exact generalization of the GVS layered
 *  group convention — the identical solver problem, so the selected sets
 *  must coincide. This is the anchor of the FOLDED coupling. */
TEST_F(TRexTikhonovTest, Equivalence_FoldedGroupKMatchesGvsTcdienet) {
    const Eigen::Index n = 150, p = 40;

    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] = j / 5;
    }

    auto run_tik = [&]() {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K   = groupMeanK(p, 5);
        ctrl.lambda_2     = 1.0;
        ctrl.prior_groups = prior_groups;
        // The defaults under test: cluster dummies, FOLDED coupling, TIKNET.
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    auto run_gvs = [&]() {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        gvs::TRexGVSControlParameter ctrl;
        ctrl.gvs_type     = gvs::GVSType::IEN;
        ctrl.prior_groups = prior_groups;
        ctrl.lambda_2     = 1.0;
        ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDIENET;
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        ctrl.trex_ctrl.solver_params.cd_max_sweeps = 200000;  // the CD solvers converge on every draw of the dummies
        gvs::TRexGVSSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const Eigen::VectorXi sel_tik = run_tik();
    const Eigen::VectorXi sel_gvs = run_gvs();

    ASSERT_EQ(sel_tik.size(), sel_gvs.size());
    EXPECT_GT(sel_tik.sum(), 0);
    EXPECT_TRUE((sel_tik.array() == sel_gvs.array()).all())
        << "FOLDED group-mean-K Tikhonov run diverged from GVS-TCDIENET.\n"
        << "  Tikhonov n_selected = " << sel_tik.sum() << "\n"
        << "  GVS      n_selected = " << sel_gvs.sum();
}


/** @brief The FOLDED group-mean-K anchor again, on a partition that MIXES
 *  groups with SINGLETONS (the shape every whole-block demo DGP produces:
 *  a few blocks against a white singleton background). The two selectors
 *  must still draw identical cluster dummies and solve the identical
 *  problem — a singleton's K entry is 1 on the Tikhonov side and a
 *  one-member group row on the GVS side. */
TEST_F(TRexTikhonovTest, Equivalence_FoldedGroupKMatchesGvsTcdienet_Singletons) {
    const Eigen::Index n = 150, p = 40, grouped = 20, gsize = 5;

    // 4 groups of 5, then 20 singletons with their own ids.
    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] =
            j < grouped ? j / gsize : grouped / gsize + (j - grouped);
    }
    Eigen::SparseMatrix<double> K(p, p);
    {
        std::vector<Eigen::Triplet<double>> trip;
        for (Eigen::Index i = 0; i < p; ++i) {
            for (Eigen::Index j = 0; j < p; ++j) {
                if (prior_groups[static_cast<std::size_t>(i)] ==
                    prior_groups[static_cast<std::size_t>(j)]) {
                    const double w = (i < grouped) ? 1.0 / static_cast<double>(gsize) : 1.0;
                    trip.emplace_back(i, j, w);
                }
            }
        }
        K.setFromTriplets(trip.begin(), trip.end());
    }

    auto run_tik = [&]() {
        GroupedData d(n, p, 0.75, 1.0, 11);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
        TRexTikhonovControlParameter ctrl;
        ctrl.tikhonov_K   = K;
        ctrl.lambda_2     = 10.0;
        ctrl.prior_groups = prior_groups;
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };
    auto run_gvs = [&]() {
        GroupedData d(n, p, 0.75, 1.0, 11);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
        gvs::TRexGVSControlParameter ctrl;
        ctrl.gvs_type     = gvs::GVSType::IEN;
        ctrl.prior_groups = prior_groups;
        ctrl.lambda_2     = 10.0;
        ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDIENET;
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        ctrl.trex_ctrl.solver_params.cd_max_sweeps = 200000;  // the CD solvers converge on every draw of the dummies
        gvs::TRexGVSSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const Eigen::VectorXi sel_tik = run_tik();
    const Eigen::VectorXi sel_gvs = run_gvs();
    ASSERT_EQ(sel_tik.size(), sel_gvs.size());
    EXPECT_TRUE((sel_tik.array() == sel_gvs.array()).all())
        << "FOLDED group-mean-K Tikhonov run diverged from GVS-TCDIENET on a"
           " partition with singletons.\n"
        << "  Tikhonov n_selected = " << sel_tik.sum() << "\n"
        << "  GVS      n_selected = " << sel_gvs.sum();
}


/** @brief The FOLDED group-mean-K anchor at DEMO control (K = 20 random
 *  experiments, HCONCAT, T-loop stagnation stop, tFDR = 0.1) on the demo
 *  partition shape (4 groups of 5 against 80 singleton nulls), with the
 *  L-loop CAPPED at L = p (max_dummy_multiplier = 1): over 40 seeds and
 *  lambda_2 in {10, 0} the Tikhonov and GVS selectors select IDENTICAL
 *  sets. Together with the DISABLED test below this pins WHERE the two
 *  selectors part company: only when the L-loop extends the cluster
 *  dummies beyond L = p (found in the TRex_Tikhonov demos). */
TEST_F(TRexTikhonovTest, Equivalence_FoldedGroupKMatchesGvsTcdienet_DemoControl_LcappedAtP) {
    const Eigen::Index n = 150, p = 100, grouped = 20, gsize = 5;
    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] =
            j < grouped ? j / gsize : grouped / gsize + (j - grouped);
    }
    Eigen::SparseMatrix<double> K(p, p);
    {
        std::vector<Eigen::Triplet<double>> trip;
        for (Eigen::Index i = 0; i < p; ++i) {
            for (Eigen::Index j = 0; j < p; ++j) {
                if (prior_groups[static_cast<std::size_t>(i)] ==
                    prior_groups[static_cast<std::size_t>(j)]) {
                    const double w = (i < grouped) ? 1.0 / static_cast<double>(gsize) : 1.0;
                    trip.emplace_back(i, j, w);
                }
            }
        }
        K.setFromTriplets(trip.begin(), trip.end());
    }

    auto demo_control = [](TRexControlParameter& t, std::size_t multiplier) {
        t.K                        = 20;
        t.max_dummy_multiplier     = multiplier;
        t.use_max_T_stop           = true;
        t.lloop_strategy           = LLoopStrategy::HCONCAT;
        t.tloop_stagnation_stop    = true;
        t.tloop_max_stagnant_steps = 5;
        t.use_memory_mapping       = false;
        t.solver_params.cd_max_sweeps = 200000;  // the CD solvers converge on every draw of the dummies
    };

    int n_diff = 0;
    for (const double lambda2 : {10.0, 0.0}) {
        for (unsigned seed = 1; seed <= 40; ++seed) {
            auto run_tik = [&]() {
                GroupedData d(n, p, 0.75, 2.0, seed);
                Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
                Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
                TRexTikhonovControlParameter ctrl;
                ctrl.tikhonov_K   = K;
                ctrl.lambda_2     = lambda2;
                ctrl.prior_groups = prior_groups;
                demo_control(ctrl.trex_ctrl, 1);
                TRexTikhonovSelector trex(Xm, ym, 0.1, ctrl, 42, false);
                return trex.select().selected_var;
            };
            auto run_gvs = [&]() {
                GroupedData d(n, p, 0.75, 2.0, seed);
                Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
                Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
                gvs::TRexGVSControlParameter ctrl;
                ctrl.gvs_type     = gvs::GVSType::IEN;
                ctrl.prior_groups = prior_groups;
                ctrl.lambda_2     = lambda2;
                ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDIENET;
                demo_control(ctrl.trex_ctrl, 1);
                gvs::TRexGVSSelector trex(Xm, ym, 0.1, ctrl, 42, false);
                return trex.select().selected_var;
            };
            const Eigen::VectorXi sel_tik = run_tik();
            const Eigen::VectorXi sel_gvs = run_gvs();
            ASSERT_EQ(sel_tik.size(), sel_gvs.size());
            if (!(sel_tik.array() == sel_gvs.array()).all()) { ++n_diff; }
        }
    }
    EXPECT_EQ(n_diff, 0)
        << "Tikhonov and GVS-TCDIENET diverged on " << n_diff
        << " of 80 cases with the L-loop capped at L = p.";
}


/** @brief The same configuration with the L-loop allowed to grow
 *  (max_dummy_multiplier = 10) on the six seeds that used to diverge
 *  (6 of 160 over seeds 1-80 x lambda_2 {10, 0}; both directions,
 *  lambda_2 = 0 included). Root cause: the dummies and the
 *  per-experiment paths were identical; the classical runner accumulated
 *  Phi as K additions of 1/K, leaving a rounding residue that decided the
 *  strict `Phi > v` comparison at the ties the 1/K voting lattice makes
 *  routine (15 x 0.05 > 0.75), while the GVS K-loop divides exact counts
 *  once (0.75, not selected). One tie flips the L-loop stop, and the two
 *  selectors then stop at different L. Fixed by exact vote counts in
 *  ExperimentRunner::aggregateResults and a single-division voting grid;
 *  this test pins the equivalence beyond L = p. */
TEST_F(TRexTikhonovTest, Equivalence_FoldedGroupK_LGrowth) {
    const Eigen::Index n = 150, p = 100, grouped = 20, gsize = 5;
    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] =
            j < grouped ? j / gsize : grouped / gsize + (j - grouped);
    }
    Eigen::SparseMatrix<double> K(p, p);
    {
        std::vector<Eigen::Triplet<double>> trip;
        for (Eigen::Index i = 0; i < p; ++i) {
            for (Eigen::Index j = 0; j < p; ++j) {
                if (prior_groups[static_cast<std::size_t>(i)] ==
                    prior_groups[static_cast<std::size_t>(j)]) {
                    const double w = (i < grouped) ? 1.0 / static_cast<double>(gsize) : 1.0;
                    trip.emplace_back(i, j, w);
                }
            }
        }
        K.setFromTriplets(trip.begin(), trip.end());
    }
    auto demo_control = [](TRexControlParameter& t) {
        t.K                        = 20;
        t.max_dummy_multiplier     = 10;
        t.use_max_T_stop           = true;
        t.lloop_strategy           = LLoopStrategy::HCONCAT;
        t.tloop_stagnation_stop    = true;
        t.tloop_max_stagnant_steps = 5;
        t.use_memory_mapping       = false;
        // With up to 10 p dummies the paths saturate, where the CD solvers
        // need more sweeps than their default budget to converge.
        t.solver_params.cd_max_sweeps = 200000;
    };
    struct Case { double lambda2; unsigned seed; };
    const std::vector<Case> cases = {{10.0, 22}, {10.0, 51}, {10.0, 63},
                                     {0.0, 30},  {0.0, 45},  {0.0, 78}};
    for (const auto& c : cases) {
        GroupedData d(n, p, 0.75, 2.0, c.seed);
        Eigen::MatrixXd X1 = d.X; Eigen::VectorXd y1 = d.y;
        Eigen::MatrixXd X2 = d.X; Eigen::VectorXd y2 = d.y;
        Eigen::Map<Eigen::MatrixXd> Xm1(X1.data(), X1.rows(), X1.cols());
        Eigen::Map<Eigen::VectorXd> ym1(y1.data(), y1.size());
        Eigen::Map<Eigen::MatrixXd> Xm2(X2.data(), X2.rows(), X2.cols());
        Eigen::Map<Eigen::VectorXd> ym2(y2.data(), y2.size());

        TRexTikhonovControlParameter tc;
        tc.tikhonov_K = K; tc.lambda_2 = c.lambda2; tc.prior_groups = prior_groups;
        demo_control(tc.trex_ctrl);
        TRexTikhonovSelector tik(Xm1, ym1, 0.1, tc, 42, false);

        gvs::TRexGVSControlParameter gc;
        gc.gvs_type = gvs::GVSType::IEN; gc.prior_groups = prior_groups;
        gc.lambda_2 = c.lambda2;
        gc.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDIENET;
        demo_control(gc.trex_ctrl);
        gvs::TRexGVSSelector g(Xm2, ym2, 0.1, gc, 42, false);

        const Eigen::VectorXi sel_tik = tik.select().selected_var;
        const Eigen::VectorXi sel_gvs = g.select().selected_var;
        EXPECT_TRUE((sel_tik.array() == sel_gvs.array()).all())
            << "lambda_2 " << c.lambda2 << ", seed " << c.seed
            << ": Tikhonov n_selected = " << sel_tik.sum()
            << ", GVS n_selected = " << sel_gvs.sum();
    }
}


/** @brief GROUP_CENTERING with an ALL-SINGLETON partition — what a HAC cut
 *  above the within-group correlation produces (the near-white regime) —
 *  is legitimate, not an error: every block of K_GC = I - J/p_m vanishes,
 *  K is the zero matrix, and the path is the plain lasso on per-variable
 *  cluster dummies. Construction must therefore succeed, the run selects,
 *  and lambda_2 has no effect on it. */
TEST_F(TRexTikhonovTest, GroupCentering_AllSingletonGroupsIsZeroKLasso) {
    const Eigen::Index n = 150, p = 40;

    std::vector<Eigen::Index> singletons(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        singletons[static_cast<std::size_t>(j)] = j;
    }
    EXPECT_EQ(TRexTikhonovSelector::buildGroupCenteringK(singletons).nonZeros(), 0);

    auto run = [&](double lambda2) {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.penalty_structure = TikPenaltyStructure::GROUP_CENTERING;
        ctrl.prior_groups      = singletons;
        ctrl.lambda_2          = lambda2;
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    Eigen::VectorXi sel_10;
    EXPECT_NO_THROW(sel_10 = run(10.0));
    const Eigen::VectorXi sel_0 = run(0.0);

    ASSERT_EQ(sel_10.size(), p);
    EXPECT_GT(sel_10.sum(), 0) << "zero-K lasso path selected nothing";
    EXPECT_TRUE((sel_10.array() == sel_0.array()).all())
        << "lambda_2 changed the run although K is zero.\n"
        << "  lambda_2 = 10: n_selected = " << sel_10.sum() << "\n"
        << "  lambda_2 = 0:  n_selected = " << sel_0.sum();
}


/** @brief End-to-end with user-supplied prior groups: the cluster-aware
 *  default path recovers the planted actives with few false selections. */
TEST_F(TRexTikhonovTest, ClusterDummies_PriorGroupsEndToEndRecovers) {
    const Eigen::Index n = 150, p = 40;
    GroupedData d(n, p, 0.75, 1.0, 7);
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexTikhonovControlParameter ctrl;
    ctrl.tikhonov_K = groupMeanK(p, 5);
    ctrl.lambda_2   = 1.0;
    ctrl.prior_groups.resize(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        ctrl.prior_groups[static_cast<std::size_t>(j)] = j / 5;
    }
    ctrl.trex_ctrl.K = 5;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
    auto result = trex.select();

    ASSERT_EQ(result.selected_var.size(), p);

    int hits = 0, false_sel = 0;
    for (Eigen::Index j = 0; j < p; ++j) {
        if (result.selected_var(j) == 1) {
            if (j < 6) { ++hits; } else { ++false_sel; }
        }
    }
    EXPECT_GE(hits, 4)
        << "Cluster-aware dummies failed to recover the planted actives.";
    EXPECT_LE(false_sel, 4)
        << "Cluster-aware dummies selected too many null variables.";
}


/** @brief The profiled IEN-CCD tuner (recommended default) resolves a
 *  positive lambda_2 and completes selection. */
TEST_F(TRexTikhonovTest, AutoLambda2_IenCcdResolvesAndRuns) {
    const Eigen::Index n = 150, p = 40;
    GroupedData d(n, p, 0.75, 1.0, 7);
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexTikhonovControlParameter ctrl;
    ctrl.tikhonov_K     = groupMeanK(p, 5);
    ctrl.lambda_2       = -1.0;
    ctrl.lambda2_method = TikLambda2Method::CV_1SE_IEN_CCD;
    ctrl.cv_n_folds     = 5;
    ctrl.trex_ctrl.K = 3;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
    auto result = trex.select();

    ASSERT_EQ(result.selected_var.size(), p);
    EXPECT_GT(trex.getLambda2Used(), 0.0);
}


/** @brief GROUP_CENTERING, which builds K_GC itself, must reproduce the
 *  SPARSE_GENERAL run with the explicit K_GC matrix: identical data, seed,
 *  lambda_2 and cluster dummies give the identical solver problem. */
TEST_F(TRexTikhonovTest, Equivalence_GroupCenteringMatchesSparseCenteringMatrix) {
    const Eigen::Index n = 150, p = 40;

    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] = j / 5;
    }

    auto run_with = [&](TikPenaltyStructure structure) {
        GroupedData d(n, p, 0.75, 1.0, 7);
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

        TRexTikhonovControlParameter ctrl;
        ctrl.penalty_structure = structure;
        if (structure == TikPenaltyStructure::SPARSE_GENERAL) {
            ctrl.tikhonov_K =
                TRexTikhonovSelector::buildGroupCenteringK(prior_groups);
        }   // GROUP_CENTERING: K left empty -> built internally
        ctrl.lambda_2     = 2.0;
        ctrl.prior_groups = prior_groups;
        ctrl.trex_ctrl.K = 5;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
        return trex.select().selected_var;
    };

    const Eigen::VectorXi sel_sparse =
        run_with(TikPenaltyStructure::SPARSE_GENERAL);
    const Eigen::VectorXi sel_built =
        run_with(TikPenaltyStructure::GROUP_CENTERING);

    ASSERT_EQ(sel_built.size(), sel_sparse.size());
    EXPECT_GT(sel_built.sum(), 0);
    EXPECT_TRUE((sel_built.array() == sel_sparse.array()).all())
        << "GROUP_CENTERING diverged from SPARSE_GENERAL with the matrix K_GC.\n"
        << "  built  n_selected = " << sel_built.sum() << "\n"
        << "  sparse n_selected = " << sel_sparse.sum();
}


// ========================================================================================
// TIKNET solver routing (native sparse-K LARS)
// ========================================================================================

/** @brief TIKNET constructs with the FOLDED coupling. The decoupled dummy
 *  block is refused: it existed with TCDTIKNET only, which the selector no
 *  longer accepts. */
TEST_F(TRexTikhonovTest, SolverGate_TiknetFoldedOnly) {

    TRexTikhonovControlParameter base;
    base.tikhonov_K = identityK(10);

    {   // TIKNET constructs on the default (FOLDED) path.
        auto ctrl = base;
        ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TIKNET;
        EXPECT_NO_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl));
    }
    {   // The decoupled dummy block has no solver left.
        auto ctrl = base;
        ctrl.fold_dummy_coupling = false;
        try {
            TRexTikhonovSelector refused(X_map, y_map, 0.1, ctrl);
            FAIL() << "fold_dummy_coupling = false was accepted";
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("fold_dummy_coupling = false is no longer available"),
                      std::string::npos) << e.what();
        }
    }
    {   // With TCDTIKNET the solver is refused first.
        auto ctrl = base;
        ctrl.fold_dummy_coupling = false;
        ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TCDTIKNET;
        EXPECT_THROW(TRexTikhonovSelector(X_map, y_map, 0.1, ctrl), std::invalid_argument);
    }
}


/** @brief End-to-end with the native LARS path solver (TIKNET) on the
 *  K_GC penalty: the planted actives (cols 0-5, group-constant signs per
 *  group) are recovered with few false selections. */
TEST_F(TRexTikhonovTest, EndToEnd_TiknetKgcRecoversPlanted) {
    const Eigen::Index n = 150, p = 40;
    GroupedData d(n, p, 0.75, 1.0, 7);
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    std::vector<Eigen::Index> prior_groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) {
        prior_groups[static_cast<std::size_t>(j)] = j / 5;
    }

    TRexTikhonovControlParameter ctrl;
    ctrl.penalty_structure = TikPenaltyStructure::GROUP_CENTERING;
    ctrl.prior_groups      = prior_groups;   // defines K_GC + cluster dummies
    ctrl.lambda_2          = 2.0;
    ctrl.trex_ctrl.solver_type = disp::SolverTypeForTRex::TIKNET;
    ctrl.trex_ctrl.K = 5;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    TRexTikhonovSelector trex(Xm, ym, 0.2, ctrl, 42, false);
    auto result = trex.select();

    ASSERT_EQ(result.selected_var.size(), p);
    int hits = 0, false_sel = 0;
    for (Eigen::Index j = 0; j < p; ++j) {
        if (result.selected_var(j) == 1) {
            if (j < 6) { ++hits; } else { ++false_sel; }
        }
    }
    EXPECT_GE(hits, 4)
        << "TIKNET + K_GC failed to recover the planted actives.";
    EXPECT_LE(false_sel, 4)
        << "TIKNET + K_GC selected too many null variables.";
}


// ========================================================================================
} /* End of namespace trex::test::trex_selector_methods::ad::trex_tikhonov */
// ========================================================================================
