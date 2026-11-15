// ========================================================================================
// test_trex_da.cpp
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// project trex includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_da/trex_da.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_da_window.hpp>

// project utils includes
#include <utils/datageneration/utils_datagen.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::trex_selector_methods::ad::trex_da {

// Namespace aliases for convenience
using namespace trex::trex_selector_methods::ad::trex_da;
using namespace trex::trex_selector_methods::ad::trex_core;
using namespace trex::utils::datageneration::datagen;
namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;
namespace disp = trex::trex_selector_methods::ad::utils::solver_dispatch;
namespace dw = trex::trex_selector_methods::ad::utils::da_window;

// ========================================================================================

/**
 * @brief Helper to set up TRexDASelector boilerplate
 *
 * @param n Number of samples
 * @param p Number of features
 * @param method DA method to use
 * @param K Number of top variables to select
 */
void run_trex_da_test(Eigen::Index n,
                      Eigen::Index p,
                      DAMethod method,
                      std::size_t K = 3,
                      LLoopStrategy l_strategy = LLoopStrategy::HCONCAT,
                      disp::SolverTypeForTRex solver = disp::SolverTypeForTRex::TLARS)
{
    // Explicit vectors prevent template resolution errors in the constructor
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};

    SyntheticData data(n, p, support, coefs, 1.0, 42);

    auto X = data.getX();
    auto y = data.getY();

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter trex_ctrl;
    trex_ctrl.K = K;
    trex_ctrl.max_dummy_multiplier = 2;
    trex_ctrl.solver_type = solver;
    trex_ctrl.solver_params.cd_max_sweeps = 200000;  // the CD solvers converge on every draw of the dummies

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = method;
    da_ctrl.trex_ctrl = trex_ctrl;

    TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);

    // Test the select() execution safely
    auto result = trex.select();

    // Retrieve DA exact result extended format
    const auto& da_result = trex.getDAResult();

    // Verify returning correct struct with fields
    EXPECT_GE(da_result.selected_var.size(), 0);
    EXPECT_EQ(da_result.method, method);

    // Structural matrix sizing depending on configuration method
    if (method == DAMethod::AR1 || method == DAMethod::EQUI) {
        EXPECT_EQ(da_result.rho_grid.size(), 0);
    } else if (method == DAMethod::BT || method == DAMethod::NN) {
        EXPECT_GT(da_result.rho_grid.size(), 0);
    }
}

/** @brief Test TRexDASelector with AR1 method in low-dimensional setting */
TEST(TRexDATest, Method_AR1_LowDimensional) {
    run_trex_da_test(300, 100, DAMethod::AR1);
}


/** @brief Test TRexDASelector with AR1 method in high-dimensional setting */
TEST(TRexDATest, Method_AR1_HighDimensional) {
    run_trex_da_test(150, 300, DAMethod::AR1);
}


/** @brief Test TRexDASelector with BT method in low-dimensional setting */
TEST(TRexDATest, Method_BT_HighDimensional) {
    run_trex_da_test(150, 300, DAMethod::BT);
}


/** @brief Test TRexDASelector with NN method in low-dimensional setting */
TEST(TRexDATest, Method_NN_LowDimensional) {
    run_trex_da_test(300, 100, DAMethod::NN);
}

/** @brief Test TRexDASelector with EQUI method in low-dimensional setting */
TEST(TRexDATest, Method_EQUI_LowDimensional) {
    run_trex_da_test(300, 100, DAMethod::EQUI);
}


// ========================================================================================
// Scaling-mode parity
// ========================================================================================

/**
 * @brief Runs TRexDASelector under a given scaling mode and returns the
 *        selected-variable indicator vector.
 *
 * The DA correlation estimators for NN (pairwise Pearson correlation) and EQUI
 * (mean off-diagonal correlation from unit-norm columns) are normalization-
 * agnostic, so L2 and ZSCORE scaling must yield bit-identical selections.
 */
Eigen::VectorXi run_trex_da_select(Eigen::Index n,
                                   Eigen::Index p,
                                   DAMethod method,
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
    trex_ctrl.max_dummy_multiplier = 2;
    trex_ctrl.scaling_mode = scaling_mode;

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = method;
    da_ctrl.trex_ctrl = trex_ctrl;

    TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.select();

    return trex.getDAResult().selected_var;
}

/** @brief DA-NN selections must be identical under L2 and ZSCORE scaling. */
TEST(TRexDATest, Method_NN_ScalingModeParity) {
    Eigen::VectorXi sel_l2 =
        run_trex_da_select(300, 100, DAMethod::NN, dn::ScalingMode::L2);
    Eigen::VectorXi sel_zscore =
        run_trex_da_select(300, 100, DAMethod::NN, dn::ScalingMode::ZSCORE);

    ASSERT_EQ(sel_l2.size(), sel_zscore.size());
    EXPECT_TRUE((sel_l2.array() == sel_zscore.array()).all())
        << "L2 and ZSCORE scaling produced different DA-NN selections.";
}

/** @brief DA-EQUI selections must be identical under L2 and ZSCORE scaling. */
TEST(TRexDATest, Method_EQUI_ScalingModeParity) {
    Eigen::VectorXi sel_l2 =
        run_trex_da_select(300, 100, DAMethod::EQUI, dn::ScalingMode::L2);
    Eigen::VectorXi sel_zscore =
        run_trex_da_select(300, 100, DAMethod::EQUI, dn::ScalingMode::ZSCORE);

    ASSERT_EQ(sel_l2.size(), sel_zscore.size());
    EXPECT_TRUE((sel_l2.array() == sel_zscore.array()).all())
        << "L2 and ZSCORE scaling produced different DA-EQUI selections.";
}


/**
 * @brief Runs a DA-BT selection under a given BT candidate-cell policy and
 *        returns the selected-variable indicator vector.
 */
Eigen::VectorXi run_trex_da_bt_select(Eigen::Index n,
                                      Eigen::Index p,
                                      BTSelectionMode bt_mode)
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
    trex_ctrl.max_dummy_multiplier = 2;

    TRexDAControlParameter da_ctrl;
    da_ctrl.method            = DAMethod::BT;
    da_ctrl.bt_selection_mode = bt_mode;
    da_ctrl.trex_ctrl         = trex_ctrl;

    TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.select();

    return trex.getDAResult().selected_var;
}

/** @brief DA-BT bt_selection_mode: FeasibleOnly is the default, and both policies
 *  produce a valid binary indicator of length p. */
TEST(TRexDATest, Method_BT_SelectionModeField) {
    const Eigen::Index n = 150, p = 300;

    // The default candidate policy must be the FDR-controlling FeasibleOnly.
    TRexDAControlParameter da_default;
    EXPECT_EQ(da_default.bt_selection_mode, BTSelectionMode::FeasibleOnly);

    Eigen::VectorXi sel_feasible =
        run_trex_da_bt_select(n, p, BTSelectionMode::FeasibleOnly);
    Eigen::VectorXi sel_rfaithful =
        run_trex_da_bt_select(n, p, BTSelectionMode::RFaithful);

    // Both policies return a valid binary indicator of length p.
    ASSERT_EQ(sel_feasible.size(), p);
    ASSERT_EQ(sel_rfaithful.size(), p);
    EXPECT_GE(sel_feasible.minCoeff(), 0);
    EXPECT_LE(sel_feasible.maxCoeff(), 1);
    EXPECT_GE(sel_rfaithful.minCoeff(), 0);
    EXPECT_LE(sel_rfaithful.maxCoeff(), 1);
}


/** @brief Test TRexDASelector with AR1 method and SKIPL L loop strategy */
TEST(TRexDATest, Method_AR1_SKIPL) {
    run_trex_da_test(150, 300, DAMethod::AR1, 3, LLoopStrategy::SKIPL);
}


/** @brief DA with the CCD lasso solver (TCDLASSO): the deflation loop drives
 *  the shared dispatch, so the CD path solver must run end-to-end. */
TEST(TRexDATest, Param_Solvers_TCDLASSO) {
    run_trex_da_test(150, 300, DAMethod::AR1, 3, LLoopStrategy::HCONCAT,
                     disp::SolverTypeForTRex::TCDLASSO);
}


/** @brief DA with the CCD elastic-net solver (TCDENET). */
TEST(TRexDATest, Param_Solvers_TCDENET) {
    run_trex_da_test(150, 300, DAMethod::AR1, 3, LLoopStrategy::HCONCAT,
                     disp::SolverTypeForTRex::TCDENET);
}


/** @brief Test TRexDASelector validation with invalid parameters */
TEST(TRexDATest, Validation_ThrowsOnInvalidParameters) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(20, 10);
    Eigen::VectorXd y = Eigen::VectorXd::Random(20);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexDAControlParameter da_params;

    // 1. Invalid rho_thr_DA
    da_params.rho_thr_DA = 1.5; // Bounds are [0, 1]
    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    }, std::invalid_argument);

    // AR1 takes thresholds in [0, 1): at 1 the window holds no neighbor.
    da_params.method = DAMethod::AR1;
    da_params.rho_thr_DA = 1.0;
    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    }, std::invalid_argument);

    // EQUI accepts 1: the correction is then always skipped.
    da_params.method = DAMethod::EQUI;
    EXPECT_NO_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    });

    da_params.method = DAMethod::BT; // Reset to the default
    da_params.rho_thr_DA = 0.5;      // Reset to valid

    // 2. Prior Groups Length Mismatch
    da_params.prior_groups = {{1, 1, 1, 1, 1}}; // Incorrect length (p = 10)
    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    }, std::invalid_argument);

    da_params.prior_groups.clear(); // Reset

    // 3. Unsupported BT Linkage
    da_params.method = DAMethod::BT;
    // Supplying a linkage method not explicitly supported by our clustering logic
    da_params.hc_linkage =
        static_cast<trex::ml_methods::clustering::hierarchical::agglomerative::LinkageMethod>(999);
    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    }, std::invalid_argument);

    // 4. PRIOR_GROUPS chosen by hand, without groups
    da_params.hc_linkage =
        trex::ml_methods::clustering::hierarchical::agglomerative::LinkageMethod::Single;
    da_params.method = DAMethod::PRIOR_GROUPS;
    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
    }, std::invalid_argument);
}


/** @brief The matrices of the selector belong to the reference threshold of
 *  the grid, the selection and its counts to the selected threshold. */
TEST(TRexDATest, Result_GridMethodsReportTheReferenceThreshold) {
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};
    SyntheticData data(100, 40, support, coefs, 1.0, 42);

    auto X = data.getX();
    auto y = data.getY();
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    for (const DAMethod method : {DAMethod::BT, DAMethod::NN}) {
        TRexDAControlParameter da_ctrl;
        da_ctrl.method = method;
        da_ctrl.trex_ctrl.K = 5;

        TRexDASelector trex(X_map, y_map, 0.2, da_ctrl, 42, false);
        trex.select();
        const auto& res = trex.getDAResult();

        // Grid of 20 thresholds, reference element 15 (index 14)
        ASSERT_EQ(res.rho_grid.size(), 20);
        ASSERT_EQ(res.FDP_hat_array_BT.size(), 20u);
        const std::size_t ref = 14;
        EXPECT_TRUE(trex.getFDPHatMat() == res.FDP_hat_array_BT[ref]);
        EXPECT_TRUE(trex.getPhiMat() == res.Phi_array_BT[ref]);

        // The result holds the arrays, its two base matrices stay empty
        EXPECT_EQ(res.FDP_hat_mat.size(), 0);
        EXPECT_EQ(res.Phi_mat.size(), 0);

        // One row per stopping time, one less for the selection counts
        const auto T_stop = static_cast<Eigen::Index>(trex.getTStop());
        EXPECT_EQ(res.FDP_hat_array_BT[ref].rows(), T_stop);
        EXPECT_EQ(res.R_array_BT[ref].rows(), T_stop > 1 ? T_stop - 1 : T_stop);

        // The selection counts of the selector belong to the selected threshold
        if (res.rho_thresh >= 0.0) {
            bool matched = false;
            for (Eigen::Index r = 0; r < res.rho_grid.size(); ++r) {
                if (res.rho_grid(r) == res.rho_thresh &&
                    trex.getRMat() == res.R_array_BT[static_cast<std::size_t>(r)]) {
                    matched = true;
                }
            }
            EXPECT_TRUE(matched);
        }
    }
}


/** @brief Test confirming SKIPL executes without exceptions */
TEST(TRexDATest, Method_AR1_SKIPL_DoesNotThrow) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(50, 10);
    Eigen::VectorXd y = Eigen::VectorXd::Random(50);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter tc_params;
    tc_params.lloop_strategy = LLoopStrategy::SKIPL;

    TRexDAControlParameter da_params;
    da_params.method = DAMethod::AR1;
    da_params.trex_ctrl = tc_params;

    EXPECT_NO_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_params);
        auto result = selector.select();
        // Since SKIPL is used, verify initialized dummies
        EXPECT_EQ(result.T_stop >= 1, true);
    });
}


// ========================================================================================
// Data Integrity
// ========================================================================================

/** @brief X is restored to its original values after select() returns (object still alive). */
TEST(TRexDATest, DataIntegrity_XRestoredAfterSelect) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(30, 10);
    Eigen::VectorXd y = Eigen::VectorXd::Random(30);
    Eigen::MatrixXd X_copy = X;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter trex_ctrl;
    trex_ctrl.K = 3;
    trex_ctrl.max_dummy_multiplier = 2;

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::AR1;
    da_ctrl.trex_ctrl = trex_ctrl;

    TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.select();

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was not restored after TRexDASelector::select().";
}


/** @brief X is restored to its original values when the object is destroyed without
 *         calling select() (normalization happens in the constructor). */
TEST(TRexDATest, DataIntegrity_XRestoredOnDestruction) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(30, 10);
    Eigen::VectorXd y = Eigen::VectorXd::Random(30);
    Eigen::MatrixXd X_copy = X;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter trex_ctrl;
    trex_ctrl.K = 3;
    trex_ctrl.max_dummy_multiplier = 2;

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::AR1;
    da_ctrl.trex_ctrl = trex_ctrl;

    {
        TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);
        // X is now normalized. Destructor fires here.
    }

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was not restored by TRexDASelector destructor.";
}


// ========================================================================================
// Correlation estimation
// ========================================================================================

namespace {

/** @brief Probe exposing protected internals (estimators, DA setup) for testing. */
class TRexDASelectorProbe : public TRexDASelector {
public:
    using TRexDASelector::TRexDASelector;
    using TRexDASelector::estimateAR1Correlation;
    using TRexDASelector::setupDA;
    using TRexDASelector::da_setup_;
    using TRexDASelector::daCorrect;
};

} // namespace


// ========================================================================================
// Deflation
// ========================================================================================

/** @brief Every method deflates all steps with the penalty of the last step.
 *
 *  Phi^DA_t(j) = Phi_t(j) / (2 - min |Phi_T(j) - Phi_T(j')|) for t = 1, ..., T
 *  over the comparison set of j: the window for AR1, all other variables for
 *  EQUI, the group mates for the BT-style methods (a division by 2 for an
 *  empty group). One factor per variable keeps the increments non-negative.
 */
TEST(TRexDATest, Deflation_PenaltyOfTheLastStepForEveryMethod) {
    const Eigen::Index p = 4;
    const Eigen::Index T = 3;
    const Eigen::MatrixXd X = Eigen::MatrixXd::Random(30, p);
    const Eigen::VectorXd y = Eigen::VectorXd::Random(30);

    // Variable 3 stays at 0.5 while its neighbour 2 rises to 0.5.
    Eigen::MatrixXd phi(p, T);
    phi << 0.2, 0.6, 1.0,
           0.2, 0.3, 0.4,
           0.0, 0.1, 0.5,
           0.1, 0.5, 0.5;
    const Eigen::VectorXd phi_last = phi.col(T - 1);

    auto deflate = [&](DAMethod method, auto&& prepare) {
        Eigen::MatrixXd X_run = X;
        Eigen::VectorXd y_run = y;
        Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());

        TRexDAControlParameter da_ctrl;
        da_ctrl.method = method;
        TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);
        prepare(trex);
        return trex.daCorrect(phi, phi_last, static_cast<std::size_t>(T));
    };

    using Sets = std::vector<std::vector<Eigen::Index>>;
    auto expect_last_step = [&](const Eigen::MatrixXd& deflated, const Sets& sets,
                                const std::string& what) {
        ASSERT_EQ(deflated.rows(), p) << what;
        ASSERT_EQ(deflated.cols(), T) << what;
        for (Eigen::Index j = 0; j < p; ++j) {
            double min_diff = 0.0;  // an empty set divides by 2
            if (!sets[static_cast<std::size_t>(j)].empty()) {
                min_diff = 2.0;
                for (const auto k : sets[static_cast<std::size_t>(j)]) {
                    min_diff = std::min(min_diff, std::abs(phi_last(j) - phi_last(k)));
                }
            }
            for (Eigen::Index t = 0; t < T; ++t) {
                EXPECT_DOUBLE_EQ(deflated(j, t), phi(j, t) / (2.0 - min_diff))
                    << what << ", variable " << j << ", step " << t + 1;
            }
        }
        for (Eigen::Index t = 1; t < T; ++t) {
            EXPECT_GE((deflated.col(t) - deflated.col(t - 1)).minCoeff(), 0.0)
                << what << ": negative increment at step " << t + 1;
        }
    };

    // AR1 with a half-width of 1: the direct neighbours.
    const DACorrectionResult ar1 = deflate(DAMethod::AR1, [](TRexDASelectorProbe& trex) {
        trex.da_setup_.kap = 1;
    });
    expect_last_step(ar1.phi_T_mat, {{1}, {0, 2}, {1, 3}, {2}}, "AR1");
    EXPECT_TRUE(ar1.Phi == ar1.phi_T_mat.col(T - 1));

    // The penalty of step 2 alone would give 0.5 / 1.6 for variable 3, above
    // its value 0.5 / 2 of step 3: a negative increment.
    EXPECT_DOUBLE_EQ(ar1.phi_T_mat(3, 1), 0.5 / 2.0);
    EXPECT_DOUBLE_EQ(ar1.phi_T_mat(3, 2), 0.5 / 2.0);

    // AR1 with an empty window deflates nothing.
    const DACorrectionResult ar1_empty = deflate(DAMethod::AR1, [](TRexDASelectorProbe& trex) {
        trex.da_setup_.kap = 0;
    });
    EXPECT_TRUE(ar1_empty.phi_T_mat == phi);
    EXPECT_TRUE(ar1_empty.Phi == phi_last);

    // EQUI: all other variables.
    const DACorrectionResult equi = deflate(DAMethod::EQUI, [](TRexDASelectorProbe& trex) {
        trex.da_setup_.cor_coef = 0.5;
    });
    expect_last_step(equi.phi_T_mat, {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}}, "EQUI");
    EXPECT_TRUE(equi.Phi == equi.phi_T_mat.col(T - 1));

    // EQUI below the correlation threshold is skipped.
    const DACorrectionResult equi_off = deflate(DAMethod::EQUI, [](TRexDASelectorProbe& trex) {
        trex.da_setup_.cor_coef = 0.01;
    });
    EXPECT_TRUE(equi_off.phi_T_mat == phi);

    // BT-style with two cut levels: the pairs {0, 1} and {2, 3}, then singletons.
    const DACorrectionResult bt = deflate(DAMethod::BT, [](TRexDASelectorProbe& trex) {
        trex.da_setup_.use_BT_style = true;
        trex.da_setup_.rho_grid_len = 2;
        trex.da_setup_.gr_j_list = {{{1}, {}}, {{0}, {}}, {{3}, {}}, {{2}, {}}};
    });
    ASSERT_EQ(bt.phi_T_array_BT.size(), 2U);
    expect_last_step(bt.phi_T_array_BT[0], {{1}, {0}, {3}, {2}}, "BT, pair level");
    expect_last_step(bt.phi_T_array_BT[1], {{}, {}, {}, {}}, "BT, singleton level");
    EXPECT_TRUE(bt.Phi_BT.col(0) == bt.phi_T_array_BT[0].col(T - 1));
    EXPECT_TRUE(bt.Phi_BT.col(1) == bt.phi_T_array_BT[1].col(T - 1));
    EXPECT_TRUE(bt.phi_T_mat == phi) << "The BT-style methods keep the undeflated history.";
}


/** @brief The row-wise Yule-Walker AR(1) estimator recovers a known rho.
 *
 *  Tolerance context: the estimator approximates the R reference (per-row
 *  `arima(..., method = "ML")` averaged over rows); both are consistent, so
 *  on a long synthetic AR(1) design the estimate must sit close to the truth.
 */
TEST(TRexDATest, Estimation_AR1YuleWalkerRecoversKnownRho) {
    const Eigen::Index n = 100;
    const Eigen::Index p = 400;
    const double rho = 0.7;

    // Row-stationary AR(1) across the ordered predictors.
    std::mt19937 rng(123);
    std::normal_distribution<double> N01(0.0, 1.0);
    const double innov_sd = std::sqrt(1.0 - rho * rho);

    Eigen::MatrixXd X(n, p);
    for (Eigen::Index i = 0; i < n; ++i) {
        X(i, 0) = N01(rng);
        for (Eigen::Index t = 1; t < p; ++t) {
            X(i, t) = rho * X(i, t - 1) + innov_sd * N01(rng);
        }
    }
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);  // unused by the estimator

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexDAControlParameter da_ctrl;
    da_ctrl.method = DAMethod::AR1;

    // The constructor normalizes X in place; the estimator runs on the
    // normalized data, exactly as setupDA_AR1() would invoke it.
    TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);

    const double est = trex.estimateAR1Correlation();
    EXPECT_NEAR(est, rho, 0.05)
        << "Yule-Walker AR(1) estimate drifted from the true coefficient.";
}

/** @brief The AR(1) window half-width follows |cor_coef| and never drops to 0.
 *
 *  kap = ceil(log(rho_thr_DA) / log(|cor_coef|)): a negative coefficient gives
 *  the window of its magnitude, a zero coefficient the minimal window, and a
 *  zero threshold, a unit coefficient or an overlong window the full window.
 */
TEST(TRexDATest, Setup_AR1WindowFollowsCoefficientMagnitude) {
    const std::size_t p = 40;
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(30, static_cast<Eigen::Index>(p));
    Eigen::VectorXd y = Eigen::VectorXd::Random(30);

    auto half_width = [&](double cor_coef, double rho_thr_DA) {
        Eigen::MatrixXd X_run = X;
        Eigen::VectorXd y_run = y;
        Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());

        TRexDAControlParameter da_ctrl;
        da_ctrl.method = DAMethod::AR1;
        da_ctrl.cor_coef = cor_coef;
        da_ctrl.rho_thr_DA = rho_thr_DA;

        TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);
        trex.setupDA();
        return trex.da_setup_.kap;
    };

    // ceil(log(0.02) / log(0.5)) = ceil(5.64) = 6
    EXPECT_EQ(half_width(0.5, 0.02), 6U);
    EXPECT_EQ(half_width(-0.5, 0.02), 6U)
        << "A negative coefficient must give the window of its magnitude.";

    EXPECT_EQ(half_width(0.0, 0.02), 1U);        // minimal window
    EXPECT_EQ(half_width(0.5, 0.0), p);          // every lag reaches the threshold
    EXPECT_EQ(half_width(1.0, 0.02), p);         // perfect correlation
    EXPECT_EQ(half_width(-1.0, 0.02), p);
    EXPECT_EQ(half_width(0.999999, 0.02), p);    // wider than the design
    EXPECT_GE(half_width(0.5, 0.999), 1U);       // never an empty window
}

/** @brief A supplied coefficient outside [-1, 1] is replaced by the estimate,
 *  with a warning. A coefficient inside the range is used as supplied. */
TEST(TRexDATest, Setup_UnusableCoefficientIsEstimated) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(40, 30);
    Eigen::VectorXd y = Eigen::VectorXd::Random(40);

    auto coefficient = [&](DAMethod method, double cor_coef, std::string& warning) {
        Eigen::MatrixXd X_run = X;
        Eigen::VectorXd y_run = y;
        Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());

        TRexDAControlParameter da_ctrl;
        da_ctrl.method = method;
        da_ctrl.cor_coef = cor_coef;

        TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);
        testing::internal::CaptureStderr();
        trex.setupDA();
        warning = testing::internal::GetCapturedStderr();
        return trex.da_setup_.cor_coef;
    };

    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto method : {DAMethod::AR1, DAMethod::EQUI}) {
        std::string warning;
        const double estimated = coefficient(method, AUTO_ESTIMATE_CORRELATION, warning);
        EXPECT_TRUE(warning.empty()) << "The sentinel requests the estimate without a warning.";

        for (const double unusable : {-1.5, 1.5, -2.5, nan}) {
            EXPECT_DOUBLE_EQ(coefficient(method, unusable, warning), estimated);
            EXPECT_NE(warning.find("cor_coef"), std::string::npos);
        }

        for (const double usable : {-1.0, -0.5, 0.0, 1.0}) {
            EXPECT_DOUBLE_EQ(coefficient(method, usable, warning), usable);
            EXPECT_TRUE(warning.empty());
        }
    }
}


/** @brief The shared AR(1) half-width: formula, limits and rejected inputs. */
TEST(TRexDAWindowTest, AR1HalfWidth_FormulaLimitsAndRejections) {
    const std::size_t p = 40;

    // ceil(log(0.02) / log(0.5)) = 6, ceil(log(0.5) / log(0.7)) = 2
    EXPECT_EQ(dw::ar1HalfWidth(0.5, 0.02, p), 6U);
    EXPECT_EQ(dw::ar1HalfWidth(-0.5, 0.02, p), 6U);
    EXPECT_EQ(dw::ar1HalfWidth(0.7, 0.5, p), 2U);

    // Limits: at least 1, at most p
    EXPECT_EQ(dw::ar1HalfWidth(0.0, 0.02, p), 1U);
    EXPECT_EQ(dw::ar1HalfWidth(1e-300, 0.999, p), 1U);
    EXPECT_EQ(dw::ar1HalfWidth(0.5, 0.0, p), p);
    EXPECT_EQ(dw::ar1HalfWidth(1.0, 0.02, p), p);
    EXPECT_EQ(dw::ar1HalfWidth(-3.0, 0.02, p), p);
    EXPECT_EQ(dw::ar1HalfWidth(0.999999, 0.02, p), p);

    // Rejected: a threshold outside [0, 1) and inputs that are not numbers
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(dw::ar1HalfWidth(0.5, 1.0, p), std::invalid_argument);
    EXPECT_THROW(dw::ar1HalfWidth(0.5, -0.1, p), std::invalid_argument);
    EXPECT_THROW(dw::ar1HalfWidth(0.5, nan, p), std::invalid_argument);
    EXPECT_THROW(dw::ar1HalfWidth(nan, 0.02, p), std::invalid_argument);
}

// ========================================================================================
// Prior groups — constrained sub-clustering
// ========================================================================================

namespace {

/** @brief Group-factor design: within-group correlated columns, independent
 *  across groups. Column j of group g is a * f_g + b * noise.
 */
Eigen::MatrixXd make_group_factor_X(Eigen::Index n,
                                    Eigen::Index p,
                                    Eigen::Index group_size,
                                    double rho_within,
                                    unsigned seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    const double a = std::sqrt(rho_within);
    const double b = std::sqrt(1.0 - rho_within);

    Eigen::MatrixXd X(n, p);
    for (Eigen::Index i = 0; i < n; ++i) {
        double factor = 0.0;
        for (Eigen::Index j = 0; j < p; ++j) {
            if (j % group_size == 0) { factor = N01(rng); }
            X(i, j) = a * factor + b * N01(rng);
        }
    }
    return X;
}

} // namespace


/** @brief PRIOR_GROUPS builds a BT-style structure by sub-clustering WITHIN
 *  the prior groups: the rho grid ends at the conservative rho = 1 singleton
 *  anchor, no neighbourhood ever crosses a prior-group boundary, and the
 *  coarsest grid point recovers (at most) the whole prior group.
 */
TEST(TRexDATest, Method_PriorGroups_ConstrainedSubclusteringStructure) {
    const Eigen::Index n = 100, p = 40, group_size = 8;

    Eigen::MatrixXd X = make_group_factor_X(n, p, group_size, 0.5, 7);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexDAControlParameter da_ctrl;
    da_ctrl.prior_groups.resize(1);
    for (Eigen::Index j = 0; j < p; ++j) {
        da_ctrl.prior_groups[0].push_back(j / group_size);
    }

    TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.setupDA();

    const auto& setup = trex.da_setup_;
    ASSERT_TRUE(setup.use_BT_style);

    // Default grid length min(20, p) = 20; grid ascending with rho = 1 anchor.
    const auto grid_len = static_cast<Eigen::Index>(setup.rho_grid_len);
    ASSERT_EQ(grid_len, 20);
    ASSERT_EQ(setup.rho_grid.size(), grid_len);
    for (Eigen::Index r = 1; r < grid_len; ++r) {
        EXPECT_LE(setup.rho_grid(r - 1), setup.rho_grid(r));
    }
    EXPECT_DOUBLE_EQ(setup.rho_grid(grid_len - 1), 1.0);

    bool any_nonempty = false;
    for (Eigen::Index j = 0; j < p; ++j) {
        for (Eigen::Index r = 0; r < grid_len; ++r) {
            const auto& mates =
                setup.gr_j_list[static_cast<std::size_t>(j)]
                               [static_cast<std::size_t>(r)];
            if (!mates.empty()) { any_nonempty = true; }
            // Hard constraint: mates never leave the prior group.
            for (const auto k : mates) {
                EXPECT_EQ(k / group_size, j / group_size)
                    << "Neighbourhood crossed a prior-group boundary.";
                EXPECT_NE(k, j);
            }
            // Tightness: never more mates than the prior group offers.
            EXPECT_LE(mates.size(),
                      static_cast<std::size_t>(group_size - 1));
        }
        // rho = 1 anchor: all neighbourhoods empty (paper's conservative
        // Psi = 1/2 penalty applies to every variable).
        EXPECT_TRUE(setup.gr_j_list[static_cast<std::size_t>(j)]
                                   [static_cast<std::size_t>(grid_len - 1)]
                        .empty());
    }
    EXPECT_TRUE(any_nonempty)
        << "Sub-clustering produced no neighbourhoods at any grid point.";
}


/** @brief The constraint partition is the finest common refinement of ALL
 *  supplied levels: with level 0 = groups of 4 and level 1 = parity, the
 *  effective constraint groups have size 2.
 */
TEST(TRexDATest, Method_PriorGroups_FinestCommonRefinement) {
    const Eigen::Index n = 60, p = 16;

    Eigen::MatrixXd X = make_group_factor_X(n, p, 4, 0.6, 11);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexDAControlParameter da_ctrl;
    da_ctrl.prior_groups.resize(2);
    for (Eigen::Index j = 0; j < p; ++j) {
        da_ctrl.prior_groups[0].push_back(j / 4);  // groups of 4
        da_ctrl.prior_groups[1].push_back(j % 2);  // parity split
    }

    TRexDASelectorProbe trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.setupDA();

    const auto& setup = trex.da_setup_;
    for (Eigen::Index j = 0; j < p; ++j) {
        for (std::size_t r = 0; r < setup.rho_grid_len; ++r) {
            const auto& mates =
                setup.gr_j_list[static_cast<std::size_t>(j)][r];
            // Refinement group = {j, j xor 2} within each block of 4: the
            // only admissible mate shares the block AND the parity.
            EXPECT_LE(mates.size(), 1u);
            for (const auto k : mates) {
                EXPECT_EQ(k / 4, j / 4);
                EXPECT_EQ(k % 2, j % 2);
            }
        }
    }
}


/** @brief Full select() run under PRIOR_GROUPS: valid binary indicator and a
 *  data-driven rho grid (not the legacy one-point-per-level grid).
 */
TEST(TRexDATest, Method_PriorGroups_Select) {
    const Eigen::Index n = 150, p = 50, group_size = 5;

    Eigen::MatrixXd X = make_group_factor_X(n, p, group_size, 0.5, 3);
    // Signal on the first column of three distinct groups.
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    beta(0) = 5.0; beta(5) = -3.0; beta(10) = 2.0;
    std::mt19937 rng(99);
    std::normal_distribution<double> N01(0.0, 1.0);
    Eigen::VectorXd y = X * beta;
    for (Eigen::Index i = 0; i < n; ++i) { y(i) += N01(rng); }

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter trex_ctrl;
    trex_ctrl.K = 3;
    trex_ctrl.max_dummy_multiplier = 2;

    TRexDAControlParameter da_ctrl;
    da_ctrl.trex_ctrl = trex_ctrl;
    da_ctrl.prior_groups.resize(1);
    for (Eigen::Index j = 0; j < p; ++j) {
        da_ctrl.prior_groups[0].push_back(j / group_size);
    }

    TRexDASelector trex(X_map, y_map, 0.1, da_ctrl, 42, false);
    trex.select();

    const auto& da_result = trex.getDAResult();
    EXPECT_EQ(da_result.method, DAMethod::PRIOR_GROUPS);
    ASSERT_EQ(da_result.selected_var.size(), p);
    EXPECT_GE(da_result.selected_var.minCoeff(), 0);
    EXPECT_LE(da_result.selected_var.maxCoeff(), 1);

    // Data-driven rho grid: hc_grid_length points ending at the rho = 1
    // anchor — NOT the legacy grid with one point per prior-group level.
    ASSERT_EQ(da_result.rho_grid.size(), 20);
    EXPECT_DOUBLE_EQ(da_result.rho_grid(da_result.rho_grid.size() - 1), 1.0);
}


/** @brief Negative prior-group labels are rejected. */
TEST(TRexDATest, Method_PriorGroups_NegativeLabelThrows) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(20, 10);
    Eigen::VectorXd y = Eigen::VectorXd::Random(20);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexDAControlParameter da_ctrl;
    da_ctrl.prior_groups = {{0, 0, 1, 1, 2, 2, 3, 3, -1, 4}};

    EXPECT_THROW({
        TRexDASelector selector(X_map, y_map, 0.1, da_ctrl);
    }, std::invalid_argument);
}


// ========================================================================================
// Exchangeable-tie band: the default leaves the choice to the selector
// ========================================================================================

namespace {

/** @brief Ten clusters of four nearly equal columns, the signal in the first three clusters. */
struct ClusteredDesign {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;

    ClusteredDesign() {
        const Eigen::Index n = 120, clusters = 10, size = 4;
        std::mt19937 rng(11);
        std::normal_distribution<double> normal(0.0, 1.0);
        X.resize(n, clusters * size);
        y = Eigen::VectorXd::Zero(n);
        for (Eigen::Index c = 0; c < clusters; ++c) {
            Eigen::VectorXd base(n);
            for (Eigen::Index i = 0; i < n; ++i) { base(i) = normal(rng); }
            for (Eigen::Index j = 0; j < size; ++j) {
                for (Eigen::Index i = 0; i < n; ++i) { X(i, c * size + j) = base(i) + 0.1 * normal(rng); }
            }
            if (c < 3) { y += 3.0 * base; }
        }
        for (Eigen::Index i = 0; i < n; ++i) { y(i) += normal(rng); }
    }
};

struct BandRun {
    std::vector<std::size_t> selected;
    Eigen::MatrixXd phi;
};

/** @brief Two runs differ when the T-loops ended at different steps or an occurrence differs. */
bool differ(const BandRun& a, const BandRun& b) {
    if (a.phi.rows() != b.phi.rows() || a.phi.cols() != b.phi.cols()) { return true; }
    return a.phi.size() > 0 && (a.phi - b.phi).cwiseAbs().maxCoeff() > 0.0;
}

/** @brief One seeded run with TOMP. `alpha` unset keeps the default of the control. */
template <class Selector>
BandRun runWithBand(const ClusteredDesign& d, const double* alpha) {
    Eigen::MatrixXd Xc = d.X;
    Eigen::VectorXd yc = d.y;
    Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
    Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());
    TRexControlParameter trex_ctrl;
    trex_ctrl.K = 50;
    trex_ctrl.max_dummy_multiplier = 2;
    trex_ctrl.solver_type = disp::SolverTypeForTRex::TOMP;
    if (alpha != nullptr) { trex_ctrl.solver_params.exch_tie_alpha = *alpha; }
    if constexpr (std::is_same_v<Selector, TRexDASelector>) {
        TRexDAControlParameter da_ctrl;
        da_ctrl.method = DAMethod::BT;
        da_ctrl.trex_ctrl = trex_ctrl;
        TRexDASelector selector(X_map, y_map, 0.2, da_ctrl, 5, false);
        selector.select();
        return {selector.getSelectedIndices(), selector.getPhiMat()};
    } else {
        TRexSelector selector(X_map, y_map, 0.2, trex_ctrl, 5, false);
        selector.select();
        return {selector.getSelectedIndices(), selector.getPhiMat()};
    }
}

} // namespace


/** @brief A band left to the selector runs with the DA width, and 0 switches it off. */
TEST(TRexDATest, Band_LeftToTheSelectorRunsWithTheDAWidth) {
    EXPECT_EQ(TRexControlParameter().solver_params.exch_tie_alpha, -1.0);
    EXPECT_EQ(DA_EXCH_TIE_ALPHA, 0.25);

    const ClusteredDesign d;
    const double da_width = DA_EXCH_TIE_ALPHA, off = 0.0;
    const BandRun by_default = runWithBand<TRexDASelector>(d, nullptr);
    const BandRun with_width = runWithBand<TRexDASelector>(d, &da_width);
    const BandRun without = runWithBand<TRexDASelector>(d, &off);

    // The default is the run with the DA width
    EXPECT_EQ(by_default.selected, with_width.selected);
    EXPECT_FALSE(differ(by_default, with_width));

    // 0 gives another run: the band draws in ten clusters of 50 experiments
    EXPECT_TRUE(differ(without, by_default));
}


/** @brief The plain selector runs without the band unless the control sets it. */
TEST(TRexDATest, Band_LeftToTheSelectorIsOffForThePlainSelector) {
    const ClusteredDesign d;
    const double off = 0.0, on = DA_EXCH_TIE_ALPHA;
    const BandRun by_default = runWithBand<TRexSelector>(d, nullptr);
    const BandRun without = runWithBand<TRexSelector>(d, &off);
    const BandRun with_width = runWithBand<TRexSelector>(d, &on);

    EXPECT_EQ(by_default.selected, without.selected);
    EXPECT_FALSE(differ(by_default, without));
    EXPECT_TRUE(differ(with_width, by_default));
}


// ========================================================================================
} /* End of namespace trex::test::trex_selector_methods::ad::trex_da */
