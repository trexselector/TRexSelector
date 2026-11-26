// ========================================================================================
// test_trex_screening.cpp
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// project trex includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_biobank_screening.hpp>

// project utils includes
#include <utils/datageneration/utils_datagen.hpp>


// ========================================================================================

// Embed into test namespace
namespace trex::test::trex_selector_methods::ad::trex_screening {

// Namespace aliases for convenience
using namespace trex::trex_selector_methods::ad::trex_screening;
using namespace trex::trex_selector_methods::ad::trex_core;
using namespace trex::utils::datageneration::datagen;
namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;
namespace disp = trex::trex_selector_methods::ad::utils::solver_dispatch;

// ========================================================================================

class TRexScreeningTest : public::testing::Test {
protected:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::Map<Eigen::MatrixXd> X_map;
    Eigen::Map<Eigen::VectorXd> y_map;

    TRexScreeningTest()
        : X(Eigen::MatrixXd::Random(50, 20)),
          y(Eigen::VectorXd::Random(50)),
          X_map(X.data(), X.rows(), X.cols()),
          y_map(y.data(), y.size()) {}
};


/** @brief Test of validation constraints */
TEST_F(TRexScreeningTest, Validation_ThrowsOnInvalidLLoopStrategy) {
    ScreenTRexControlParameter screen_params;

    // Supported: STANDARD, SEEDED, PERMUTATION, PERMUTATION_SEEDED.
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));

    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::SEEDED;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));

    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::PERMUTATION;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));

    screen_params.trex_ctrl.lloop_strategy =
        LLoopStrategy::PERMUTATION_SEEDED;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));

    // Unsupported: the L-loop growth strategies.
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::SKIPL;
    EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params),
                                    std::invalid_argument);

    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::HCONCAT;
    EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params),
                 std::invalid_argument);
}


/** @brief `rho_thr_DA` is checked at construction: [0, 1], and [0, 1) for
 *  TREX_DA_AR1, whose window has no lag left at a threshold of 1. */
TEST_F(TRexScreeningTest, Validation_ThrowsOnInvalidRhoThreshold) {
    ScreenTRexControlParameter screen_params;

    for (const double bad : {-0.1, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
        screen_params.rho_thr_DA = bad;
        for (const auto method : {ScreenTRexMethod::TREX, ScreenTRexMethod::TREX_DA_AR1,
                                  ScreenTRexMethod::TREX_DA_EQUI}) {
            screen_params.trex_method = method;
            EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
        }
    }

    screen_params.rho_thr_DA = 1.0;
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_AR1;
    EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_NN;
    EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_EQUI;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_BLOCK_EQUI;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));
}


/** @brief TREX_DA_BLOCK_EQUI takes between 1 and p blocks. The other
 *  variants do not read `n_blocks`. */
TEST_F(TRexScreeningTest, Validation_ThrowsOnInvalidBlockCount) {
    ScreenTRexControlParameter screen_params;
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_BLOCK_EQUI;

    for (const std::size_t bad : {std::size_t{0}, std::size_t{21}}) {  // p = 20
        screen_params.n_blocks = bad;
        EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
    }
    for (const std::size_t good : {std::size_t{1}, std::size_t{20}}) {
        screen_params.n_blocks = good;
        EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));
    }

    screen_params.n_blocks = 0;
    screen_params.trex_method = ScreenTRexMethod::TREX;
    EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));
}


/** @brief The bootstrap parameters are checked at construction, whether or
 *  not the bootstrap variant is on: `R_boot` at least 2 (the standard error
 *  needs two resamples) and `ci_grid_step` in (0, 1] (a step of 0 would never
 *  advance the grid over the confidence level). */
TEST_F(TRexScreeningTest, Validation_ThrowsOnInvalidBootstrapParameters) {
    for (const bool use_bootstrap : {false, true}) {
        ScreenTRexControlParameter screen_params;
        screen_params.use_bootstrap_CI = use_bootstrap;

        for (const std::size_t bad : {std::size_t{0}, std::size_t{1}}) {
            screen_params.R_boot = bad;
            EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
        }
        screen_params.R_boot = 2;
        EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));

        for (const double bad : {0.0, -0.001, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
            screen_params.ci_grid_step = bad;
            EXPECT_THROW(ScreenTRexSelector(X_map, y_map, screen_params), std::invalid_argument);
        }
        screen_params.ci_grid_step = 1.0;
        EXPECT_NO_THROW(ScreenTRexSelector(X_map, y_map, screen_params));
    }
}


/** @brief The estimated FDR is one over the number of selected variables, and
 *  1 when nothing is selected. A response orthogonal to the design and to the
 *  intercept leaves every variable behind the first dummy, so nothing is
 *  selected. */
TEST(TRexScreeningEstimateTest, Estimate_EmptySelectionReportsOne) {
    std::mt19937 rng(21);
    std::normal_distribution<double> N01(0.0, 1.0);
    const Eigen::Index n = 80;
    const Eigen::Index p = 40;
    Eigen::MatrixXd X(n, p);
    Eigen::VectorXd noise(n);
    for (Eigen::Index j = 0; j < p; ++j) {
        for (Eigen::Index i = 0; i < n; ++i) { X(i, j) = N01(rng); }
    }
    for (Eigen::Index i = 0; i < n; ++i) { noise(i) = N01(rng); }

    // Residual of the noise after the projection onto [1, X].
    Eigen::MatrixXd basis(n, p + 1);
    basis << Eigen::VectorXd::Ones(n), X;
    const Eigen::VectorXd y_orth =
        noise - basis * basis.colPivHouseholderQr().solve(noise);

    for (const bool use_bootstrap : {false, true}) {
        Eigen::MatrixXd X_run = X;
        Eigen::VectorXd y_run = y_orth;
        Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());

        ScreenTRexControlParameter screen_params;
        screen_params.use_bootstrap_CI = use_bootstrap;
        screen_params.R_boot = 100;

        ScreenTRexSelector selector(X_map, y_map, screen_params, 5, false);
        selector.select();
        const auto& result = selector.getScreenResult();

        EXPECT_EQ(result.selected_var.sum(), 0);
        EXPECT_DOUBLE_EQ(result.estimated_FDR, 1.0);
    }

    // With selections the estimate is one over their number.
    Eigen::MatrixXd X_run = X;
    Eigen::VectorXd y_signal = 5.0 * X.col(1) - 5.0 * X.col(9) + 5.0 * X.col(20) + 0.1 * noise;
    Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y_signal.data(), y_signal.size());
    ScreenTRexSelector selector(X_map, y_map, ScreenTRexControlParameter(), 5, false);
    selector.select();
    const auto& result = selector.getScreenResult();
    ASSERT_GT(result.selected_var.sum(), 0);
    EXPECT_DOUBLE_EQ(result.estimated_FDR, 1.0 / static_cast<double>(result.selected_var.sum()));
}


/** @brief TREX_DA_AR1 at the edges of its window: a zero threshold is the
 *  full window, and a negative coefficient gives the window of its magnitude. */
TEST(TRexScreeningWindowTest, Method_DA_AR1_WindowEdges) {
    const Eigen::Index n = 60;
    const Eigen::Index p = 30;
    std::srand(11);
    const Eigen::MatrixXd X0 = Eigen::MatrixXd::Random(n, p);
    const Eigen::VectorXd y0 =
        3.0 * X0.col(0) - 2.0 * X0.col(7) + 2.5 * X0.col(15) + 0.1 * Eigen::VectorXd::Random(n);

    auto run = [&](double cor_coef, double rho_thr_DA) {
        Eigen::MatrixXd X = X0;
        Eigen::VectorXd y = y0;
        Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

        ScreenTRexControlParameter screen_params;
        screen_params.trex_method = ScreenTRexMethod::TREX_DA_AR1;
        screen_params.cor_coef = cor_coef;
        screen_params.rho_thr_DA = rho_thr_DA;
        screen_params.trex_ctrl.K = 10;

        ScreenTRexSelector selector(X_map, y_map, screen_params, 7, false);
        return selector.select();
    };

    // The window follows |cor_coef|.
    const auto pos = run(0.5, 0.02);
    const auto neg = run(-0.5, 0.02);
    EXPECT_EQ(pos.selected_var, neg.selected_var);
    EXPECT_TRUE(pos.Phi_mat.isApprox(neg.Phi_mat));

    // A zero threshold is the full window: the run completes with finite
    // values and equals the run with a coefficient of magnitude 1.
    const auto full = run(0.5, 0.0);
    const auto unit = run(1.0, 0.02);
    EXPECT_TRUE(full.Phi_mat.allFinite());
    EXPECT_EQ(full.selected_var, unit.selected_var);
    EXPECT_TRUE(full.Phi_mat.isApprox(unit.Phi_mat));
}


/** @brief TREX_DA_NN compares a variable with the neighbors inside the AR(1)
 *  window whose sample correlation reaches the threshold. A zero threshold
 *  keeps every neighbor, so the variant equals TREX_DA_AR1 there. A variable
 *  without such a neighbor keeps its relative occurrence. */
TEST(TRexScreeningWindowTest, Method_DA_NN_GroupsAndEmptyGroup) {
    const Eigen::Index n = 60;
    const Eigen::Index p = 20;
    std::srand(13);

    auto run = [&](const Eigen::MatrixXd& X0, const Eigen::VectorXd& y0, ScreenTRexMethod method,
                   double cor_coef, double rho_thr_DA) {
        Eigen::MatrixXd X = X0;
        Eigen::VectorXd y = y0;
        Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

        ScreenTRexControlParameter screen_params;
        screen_params.trex_method = method;
        screen_params.cor_coef = cor_coef;
        screen_params.rho_thr_DA = rho_thr_DA;
        screen_params.trex_ctrl.K = 10;

        ScreenTRexSelector selector(X_map, y_map, screen_params, 7, false);
        return selector.select();
    };

    // 1. A zero threshold: the whole design is the window and every neighbor stays.
    const Eigen::MatrixXd X_any = Eigen::MatrixXd::Random(n, p);
    const Eigen::VectorXd y_any =
        3.0 * X_any.col(0) - 2.0 * X_any.col(7) + 2.5 * X_any.col(15) + 0.1 * Eigen::VectorXd::Random(n);
    const auto nn_all = run(X_any, y_any, ScreenTRexMethod::TREX_DA_NN, 0.5, 0.0);
    const auto ar1_all = run(X_any, y_any, ScreenTRexMethod::TREX_DA_AR1, 0.5, 0.0);
    EXPECT_EQ(nn_all.selected_var, ar1_all.selected_var);
    EXPECT_TRUE(nn_all.Phi_mat.isApprox(ar1_all.Phi_mat));

    // 2. Centered orthonormal columns: every sample correlation is 0, so no
    //    variable has a neighbor and the variant equals the plain screening.
    Eigen::MatrixXd centered = Eigen::MatrixXd::Random(n, p);
    centered.rowwise() -= centered.colwise().mean();
    const Eigen::MatrixXd X_orth =
        Eigen::HouseholderQR<Eigen::MatrixXd>(centered).householderQ() * Eigen::MatrixXd::Identity(n, p);
    const Eigen::VectorXd y_orth =
        3.0 * X_orth.col(0) - 2.0 * X_orth.col(7) + 2.5 * X_orth.col(15) + 0.01 * Eigen::VectorXd::Random(n);
    const auto nn_none = run(X_orth, y_orth, ScreenTRexMethod::TREX_DA_NN, 0.5, 0.5);
    const auto plain = run(X_orth, y_orth, ScreenTRexMethod::TREX, 0.5, 0.5);
    EXPECT_EQ(nn_none.selected_var, plain.selected_var);
    EXPECT_TRUE(nn_none.Phi_mat.isApprox(plain.Phi_mat));
    EXPECT_GT(nn_none.selected_var.sum(), 0)
        << "A penalty of one half for an empty group would leave nothing above the voting level.";
}


/** @brief A supplied coefficient outside [-1, 1] is replaced by the estimate,
 *  with a warning. A coefficient inside the range is used as supplied. */
TEST(TRexScreeningWindowTest, Method_DA_UnusableCoefficientIsEstimated) {
    std::srand(17);
    const Eigen::MatrixXd X0 = Eigen::MatrixXd::Random(60, 30);
    const Eigen::VectorXd y0 = 3.0 * X0.col(0) - 2.0 * X0.col(7) + 0.1 * Eigen::VectorXd::Random(60);

    auto coefficient = [&](ScreenTRexMethod method, double cor_coef, std::string& warning) {
        Eigen::MatrixXd X = X0;
        Eigen::VectorXd y = y0;
        Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

        ScreenTRexControlParameter screen_params;
        screen_params.trex_method = method;
        screen_params.cor_coef = cor_coef;
        screen_params.trex_ctrl.K = 5;

        ScreenTRexSelector selector(X_map, y_map, screen_params, 7, false);
        testing::internal::CaptureStderr();
        selector.select();
        warning = testing::internal::GetCapturedStderr();
        return selector.getScreenResult().estimated_correlation;
    };

    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto method : {ScreenTRexMethod::TREX_DA_AR1, ScreenTRexMethod::TREX_DA_EQUI,
                              ScreenTRexMethod::TREX_DA_BLOCK_EQUI, ScreenTRexMethod::TREX_DA_NN}) {
        std::string warning;
        const double estimated = coefficient(method, AUTO_ESTIMATE_CORRELATION, warning);
        EXPECT_TRUE(warning.empty()) << "The sentinel requests the estimate without a warning.";

        for (const double unusable : {-1.5, 1.5, -2.5, nan}) {
            EXPECT_DOUBLE_EQ(coefficient(method, unusable, warning), estimated);
            EXPECT_NE(warning.find("cor_coef"), std::string::npos);
        }

        EXPECT_DOUBLE_EQ(coefficient(method, -0.5, warning), -0.5);
        EXPECT_TRUE(warning.empty());
    }
}


/** @brief PERMUTATION and PERMUTATION_SEEDED key their base off the same
 *  resolved permutation_base_seed_, so for a fixed seed the stored and the
 *  stateless variants must produce bit-identical screening results. */
TEST_F(TRexScreeningTest, Equivalence_PermutationSeededMatchesPermutation) {

    auto run_with = [&](LLoopStrategy strategy) {
        Eigen::MatrixXd Xc = X;
        Eigen::VectorXd yc = y;
        Eigen::Map<Eigen::MatrixXd> Xm(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::VectorXd> ym(yc.data(), yc.size());

        ScreenTRexControlParameter screen_params;
        screen_params.trex_method = ScreenTRexMethod::TREX;
        screen_params.use_bootstrap_CI = false;
        screen_params.trex_ctrl.K = 5;
        screen_params.trex_ctrl.lloop_strategy = strategy;

        ScreenTRexSelector selector(Xm, ym, screen_params, /*seed=*/42);
        return selector.select();
    };

    const auto res_perm = run_with(LLoopStrategy::PERMUTATION);
    const auto res_pod  = run_with(LLoopStrategy::PERMUTATION_SEEDED);

    ASSERT_EQ(res_perm.selected_var.size(), res_pod.selected_var.size());
    EXPECT_TRUE((res_perm.selected_var.array()
                 == res_pod.selected_var.array()).all())
        << "PERMUTATION_SEEDED selection diverged from PERMUTATION.";
    ASSERT_EQ(res_perm.Phi_prime.size(), res_pod.Phi_prime.size());
    EXPECT_TRUE(res_perm.Phi_prime.isApprox(res_pod.Phi_prime, 1e-14))
        << "PERMUTATION_SEEDED Phi diverged from PERMUTATION.";
}


/** @brief SEEDED regenerates the identical draws STANDARD stores (both keyed
 *  on the first L step, screening runs at exactly L = p), one D_k at a time
 *  — the two must produce bit-identical screening results. */
TEST_F(TRexScreeningTest, Equivalence_SeededMatchesStandard) {

    auto run_with = [&](LLoopStrategy strategy) {
        Eigen::MatrixXd Xc = X;
        Eigen::VectorXd yc = y;
        Eigen::Map<Eigen::MatrixXd> Xm(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::VectorXd> ym(yc.data(), yc.size());

        ScreenTRexControlParameter screen_params;
        screen_params.trex_method = ScreenTRexMethod::TREX;
        screen_params.use_bootstrap_CI = false;
        screen_params.trex_ctrl.K = 5;
        screen_params.trex_ctrl.lloop_strategy = strategy;

        ScreenTRexSelector selector(Xm, ym, screen_params, /*seed=*/42);
        return selector.select();
    };

    const auto res_std = run_with(LLoopStrategy::STANDARD);
    const auto res_sd  = run_with(LLoopStrategy::SEEDED);

    ASSERT_EQ(res_std.selected_var.size(), res_sd.selected_var.size());
    EXPECT_TRUE((res_std.selected_var.array()
                 == res_sd.selected_var.array()).all())
        << "SEEDED selection diverged from STANDARD.";
    ASSERT_EQ(res_std.Phi_prime.size(), res_sd.Phi_prime.size());
    EXPECT_TRUE(res_std.Phi_prime.isApprox(res_sd.Phi_prime, 1e-14))
        << "SEEDED Phi diverged from STANDARD.";
}


/** @brief The getters of the base class hold the fixed values of a screening
 *         run, and the result holds the fields of the screening. */
TEST_F(TRexScreeningTest, Result_BaseGettersHoldTheFixedValuesOfAScreeningRun) {
    ScreenTRexControlParameter screen_params;
    screen_params.trex_ctrl.K = 5;

    ScreenTRexSelector selector(X_map, y_map, screen_params, 42, false);
    selector.select();
    const auto& res = selector.getScreenResult();
    const auto p = X.cols();
    const double n_selected = static_cast<double>(res.selected_var.sum());

    // One pass at T = 1 and L = p, without a target and without a calibration
    EXPECT_EQ(selector.getTStop(), 1u);
    EXPECT_EQ(selector.getDummyMultiplierL(), 1u);
    EXPECT_EQ(selector.getNumDummies(), static_cast<std::size_t>(p));
    EXPECT_DOUBLE_EQ(selector.getVotingThreshold(), 0.5);
    EXPECT_DOUBLE_EQ(selector.getTFDR(), 0.0);

    ASSERT_EQ(selector.getRMat().size(), 1);
    EXPECT_DOUBLE_EQ(selector.getRMat()(0, 0), n_selected);
    ASSERT_EQ(selector.getFDPHatMat().size(), 1);
    EXPECT_DOUBLE_EQ(selector.getFDPHatMat()(0, 0), 0.0);
    ASSERT_EQ(selector.getPhiMat().rows(), 1);
    EXPECT_TRUE(selector.getPhiPrime() == selector.getPhiMat().row(0).transpose());

    // The FDR estimate is a field of the result
    EXPECT_DOUBLE_EQ(res.estimated_FDR, 1.0 / std::max(1.0, n_selected));

    // The fields of the screening without bootstrap and without a filter
    EXPECT_FALSE(res.used_bootstrap);
    EXPECT_DOUBLE_EQ(res.confidence_level, 0.0);
    EXPECT_DOUBLE_EQ(res.estimated_correlation, 0.0);
    EXPECT_EQ(res.beta_mat.size(), p);
    EXPECT_GT(res.dummy_betas.size(), 0);
    EXPECT_TRUE((res.dummy_betas.array() != 0.0).all());

    // A supplied coefficient of a dependency-aware variant is reported as it is
    screen_params.trex_method = ScreenTRexMethod::TREX_DA_AR1;
    screen_params.cor_coef = 0.5;
    ScreenTRexSelector filtered(X_map, y_map, screen_params, 42, false);
    filtered.select();
    EXPECT_DOUBLE_EQ(filtered.getScreenResult().estimated_correlation, 0.5);
}


/** @brief Test execution of Ordinary Screening */
TEST_F(TRexScreeningTest, Execution_OrdinaryScreenTRexDoesNotThrow) {
    ScreenTRexControlParameter screen_params;
    screen_params.trex_method = ScreenTRexMethod::TREX;
    screen_params.use_bootstrap_CI = false;
    screen_params.trex_ctrl.K = 5; // Reduced K for quick unit testing
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;

    EXPECT_NO_THROW({
        ScreenTRexSelector selector(X_map, y_map, screen_params);
        auto result = selector.select();

        // Ensure static execution values
        EXPECT_EQ(result.T_stop, 1);

        // num_dummies is locked to p
        EXPECT_EQ(result.num_dummies, 20);
    });
}


/** @brief Test of Boostrap CI Screening */
TEST_F(TRexScreeningTest, Execution_BootstrapScreenTRexDoesNotThrow) {
    ScreenTRexControlParameter screen_params;
    screen_params.trex_method = ScreenTRexMethod::TREX;
    screen_params.use_bootstrap_CI = true;
    screen_params.R_boot = 50; // Small bootstrap resamples for fast testing
    screen_params.trex_ctrl.K = 5;
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;

    EXPECT_NO_THROW({
        ScreenTRexSelector selector(X_map, y_map, screen_params);

        auto result = selector.getScreenResult();
        selector.select();
    });
}

/// ========================================================================================
// Data Integrity
// ========================================================================================

/** @brief X is restored to its original values after select() returns (object still alive). */
TEST_F(TRexScreeningTest, DataIntegrity_XRestoredAfterSelect) {
    Eigen::MatrixXd X_copy = X;

    ScreenTRexControlParameter screen_params;
    screen_params.trex_ctrl.K = 3;
    screen_params.trex_ctrl.max_dummy_multiplier = 2;
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;

    ScreenTRexSelector selector(X_map, y_map, screen_params, 42, false);
    selector.select();

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was not restored after ScreenTRexSelector::select().";
}


/** @brief X is restored to its original values when the object is destroyed without
 *         calling select() (normalization happens in the constructor). */
TEST_F(TRexScreeningTest, DataIntegrity_XRestoredOnDestruction) {
    Eigen::MatrixXd X_copy = X;

    ScreenTRexControlParameter screen_params;
    screen_params.trex_ctrl.K = 3;
    screen_params.trex_ctrl.max_dummy_multiplier = 2;
    screen_params.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;

    {
        ScreenTRexSelector selector(X_map, y_map, screen_params, 42, false);
        // X is now normalized. Destructor fires here.
    }

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was not restored by ScreenTRexSelector destructor.";
}

// ========================================================================================
// Scaling-mode parity (DA-EQUI / DA-BLOCK-EQUI)
// ========================================================================================

/**
 * @brief Runs ScreenTRexSelector under a given scaling mode and returns the
 *        selected-variable indicator vector.
 *
 * estimateEquiCorrelation() / estimateBlockEquiCorrelation() accumulate
 * unit-norm columns (scale-mode agnostic, mirroring
 * TRexDASelector::estimateEquiCorrelation), so L2 and ZSCORE scaling must
 * yield bit-identical selections for the DA_EQUI and DA_BLOCK_EQUI variants.
 */
Eigen::VectorXi run_screen_trex_select(Eigen::Index n,
                                       Eigen::Index p,
                                       ScreenTRexMethod method,
                                       dn::ScalingMode scaling_mode)
{
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};

    SyntheticData data(n, p, support, coefs, 1.0, 42);

    auto X = data.getX();
    auto y = data.getY();

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter trex_ctrl;
    trex_ctrl.K = 3;
    trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;
    trex_ctrl.scaling_mode = scaling_mode;

    ScreenTRexControlParameter screen_ctrl;
    screen_ctrl.trex_method = method;
    screen_ctrl.trex_ctrl = trex_ctrl;

    ScreenTRexSelector selector(X_map, y_map, screen_ctrl, 42, false);
    auto result = selector.select();

    return result.selected_var;
}

/** @brief Screen-TRex DA-EQUI selections must be identical under L2 and ZSCORE scaling. */
TEST(TRexScreeningScalingTest, Method_DA_EQUI_ScalingModeParity) {
    Eigen::VectorXi sel_l2 =
        run_screen_trex_select(300, 100, ScreenTRexMethod::TREX_DA_EQUI, dn::ScalingMode::L2);
    Eigen::VectorXi sel_zscore =
        run_screen_trex_select(300, 100, ScreenTRexMethod::TREX_DA_EQUI, dn::ScalingMode::ZSCORE);

    ASSERT_EQ(sel_l2.size(), sel_zscore.size());
    EXPECT_TRUE((sel_l2.array() == sel_zscore.array()).all())
        << "L2 and ZSCORE scaling produced different Screen-TRex DA-EQUI selections.";
}

/** @brief Screen-TRex DA-BLOCK-EQUI selections must be identical under L2 and ZSCORE scaling. */
TEST(TRexScreeningScalingTest, Method_DA_BLOCK_EQUI_ScalingModeParity) {
    Eigen::VectorXi sel_l2 =
        run_screen_trex_select(300, 100, ScreenTRexMethod::TREX_DA_BLOCK_EQUI, dn::ScalingMode::L2);
    Eigen::VectorXi sel_zscore =
        run_screen_trex_select(300, 100, ScreenTRexMethod::TREX_DA_BLOCK_EQUI, dn::ScalingMode::ZSCORE);

    ASSERT_EQ(sel_l2.size(), sel_zscore.size());
    EXPECT_TRUE((sel_l2.array() == sel_zscore.array()).all())
        << "L2 and ZSCORE scaling produced different Screen-TRex DA-BLOCK-EQUI selections.";
}

/// ========================================================================================
} /* End of namespace trex::test::trex_selector_methods::ad::trex_screening */
