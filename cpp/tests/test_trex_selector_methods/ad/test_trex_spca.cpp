// ========================================================================================
// test_trex_spca.cpp
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
#include <algorithm>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

// project trex includes
#include <trex_selector_methods/ad/trex_spca/trex_spca.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::trex_selector_methods::ad::trex_spca {

// Namespace alias for convenience
using namespace trex::trex_selector_methods::ad::trex_spca;
namespace tc = trex::trex_selector_methods::ad::trex_core;

// ========================================================================================

/**
 * @brief Test fixture providing a minimal (n=30, p=10) design matrix shared
 *        across all TRexSPCASelector tests.
 *
 * @details Dimensions are kept deliberately small so that the gtest_discover_tests
 *          post-build enumeration step (5 s limit) is never exceeded even under
 *          a heavily-loaded parallel build. A fixed seed ensures reproducibility.
 */
class TRexSPCATest : public ::testing::Test {
protected:
    static constexpr Eigen::Index n = 30;
    static constexpr Eigen::Index p = 10;
    static constexpr Eigen::Index M = 2;
    static constexpr int          seed = 42;

    Eigen::MatrixXd X;
    Eigen::Map<Eigen::MatrixXd> X_map;

    /** @brief Returns a TRexSPCAControl tuned for fast unit testing. */
    static TRexSPCAControl fast_ctrl(SPCAMode mode = SPCAMode::ACTIVE_SET) {
        TRexSPCAControl ctrl;
        ctrl.mode                           = mode;
        ctrl.lambda2_ridge_loadings         = 1e-6;
        ctrl.gvs_ctrl.trex_ctrl.K                    = 3;
        ctrl.gvs_ctrl.trex_ctrl.max_dummy_multiplier = 2;
        return ctrl;
    }

    TRexSPCATest()
        : X(Eigen::MatrixXd::Random(n, p)),
          X_map(X.data(), n, p) {}

    /** @brief Deterministic single-factor design engineered so the GVS active
     *         set of PC1 and the top-|A| ordinary-loading support diverge.
     *
     *  Ordinary loadings rank variables by covariance with z_1 (X^T z_1 is
     *  proportional to v_1), while the GVS sub-selector ranks by normalized
     *  correlation. The design therefore plants, next to the strong factor
     *  columns {0, 1, 2}:
     *    - columns {3, 4, 5}: tiny norm, near-unit correlation with the factor
     *      (GVS selects them, but their ordinary loadings are small), and
     *    - columns {6, 7, 8}: large norm, weak correlation (GVS rejects them,
     *      but their ordinary loadings can outrank columns {3, 4, 5}).
     *  The remaining columns are pure noise, so that the FDR-controlled active
     *  set is a proper subset of {1, ..., p} (with very small p the selector
     *  admits every column and the two supports trivially coincide).
     *  Uses its own RNG rather than Eigen::Random so the design is independent
     *  of test order. */
    static constexpr Eigen::Index n_f = 40;
    static constexpr Eigen::Index p_f = 20;

    static Eigen::MatrixXd plantedFactorX() {
        std::mt19937 rng(7u);
        std::normal_distribution<double> N01(0.0, 1.0);

        Eigen::VectorXd z(n_f);
        for (Eigen::Index i = 0; i < n_f; ++i) {
            z(i) = 3.0 * N01(rng);
        }

        Eigen::MatrixXd Xf(n_f, p_f);
        for (Eigen::Index j = 0; j < p_f; ++j) {
            for (Eigen::Index i = 0; i < n_f; ++i) {
                Xf(i, j) = N01(rng);
            }
        }
        for (Eigen::Index j = 0; j < 3; ++j) {
            Xf.col(j) = 0.9 * z + 0.3 * Xf.col(j);   // strong factor columns
        }
        for (Eigen::Index j = 3; j < 6; ++j) {
            Xf.col(j) = 0.05 * z + 0.02 * Xf.col(j); // tiny norm, corr ~ 1
        }
        for (Eigen::Index j = 6; j < 9; ++j) {
            Xf.col(j) = 0.05 * z + 1.5 * Xf.col(j);  // large norm, weak corr
        }
        return Xf;                                   // cols 9..19: pure noise
    }

    /** @brief Shared invariant check: after select(), Z must equal
     *         X_centered * V (scores computed over the full support of V, not
     *         merely the GVS active set), and every assembled loading column
     *         must be L2-normalized. */
    static void expectScoresMatchCenteredXTimesV(SPCAMode mode) {
        Eigen::MatrixXd Xf = plantedFactorX();
        Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), n_f, p_f);

        TRexSPCASelector spca(Xf_map, M, 0.1, fast_ctrl(mode), seed);
        auto result = spca.select();

        // The planted factor must be picked up, otherwise the invariant below
        // would only be checked on the trivial all-zero support.
        ASSERT_FALSE(result.active_sets[0].empty());

        const Eigen::RowVectorXd mu = Xf.colwise().mean();
        const Eigen::MatrixXd Xc = Xf.rowwise() - mu;
        const Eigen::MatrixXd Z_expected = Xc * result.V;

        EXPECT_LE((result.Z - Z_expected).norm(),
                  1e-10 * std::max(1.0, Z_expected.norm()))
            << "Z != X_centered * V for mode "
            << (mode == SPCAMode::ACTIVE_SET ? "ActiveSet" : "Thresholded");

        for (Eigen::Index m = 0; m < M; ++m) {
            const double norm_m = result.V.col(m).norm();
            if (norm_m > 0.0) {
                EXPECT_NEAR(norm_m, 1.0, 1e-12);
            }
        }
    }
};


// ========================================================================================
// Execution — does not throw
// ========================================================================================

/** @brief ActiveSet mode runs without throwing on valid inputs. */
TEST_F(TRexSPCATest, Execution_ActiveSetModeDoesNotThrow) {
    EXPECT_NO_THROW({
        TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(SPCAMode::ACTIVE_SET), seed);
        spca.select();
    });
}


/** @brief Thresholded mode runs without throwing on valid inputs. */
TEST_F(TRexSPCATest, Execution_ThresholdedModeDoesNotThrow) {
    EXPECT_NO_THROW({
        TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(SPCAMode::THRESHOLDED), seed);
        spca.select();
    });
}


// ========================================================================================
// Output dimensions
// ========================================================================================

/** @brief Result matrices have the expected shapes. */
TEST_F(TRexSPCATest, Dimensions_ResultMatricesCorrectShape) {
    TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
    auto result = spca.select();

    EXPECT_EQ(result.Z.rows(), n);
    EXPECT_EQ(result.Z.cols(), M);
    EXPECT_EQ(result.V.rows(), p);
    EXPECT_EQ(result.V.cols(), M);
    EXPECT_EQ(result.active_sets.size(), static_cast<std::size_t>(M));
    EXPECT_EQ(result.adjusted_ev.size(),  M);
    EXPECT_EQ(result.cumulative_ev.size(), M);
}


/** @brief Extracting only one sparse PC (M=1) produces correct shapes. */
TEST_F(TRexSPCATest, Dimensions_SingleComponentCorrectShape) {
    TRexSPCASelector spca(X_map, 1, 0.2, fast_ctrl(), seed);
    auto result = spca.select();

    EXPECT_EQ(result.Z.rows(), n);
    EXPECT_EQ(result.Z.cols(), 1);
    EXPECT_EQ(result.V.rows(), p);
    EXPECT_EQ(result.V.cols(), 1);
    EXPECT_EQ(result.active_sets.size(), 1u);
}


// ========================================================================================
// Active sets
// ========================================================================================

/** @brief All active-set indices are valid column indices of X. */
TEST_F(TRexSPCATest, ActiveSets_IndicesInBounds) {
    TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
    auto result = spca.select();

    for (const auto& as : result.active_sets) {
        for (Eigen::Index idx : as) {
            EXPECT_GE(idx, 0);
            EXPECT_LT(idx, p);
        }
    }
}


// ========================================================================================
// Explained variance
// ========================================================================================

/** @brief Adjusted explained variances are non-negative. */
TEST_F(TRexSPCATest, ExplainedVariance_AdjustedNonNegative) {
    TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
    auto result = spca.select();

    for (Eigen::Index k = 0; k < M; ++k) {
        EXPECT_GE(result.adjusted_ev(k), 0.0);
    }
}


/** @brief Cumulative explained variance is monotonically non-decreasing and in [0, 1]. */
TEST_F(TRexSPCATest, ExplainedVariance_CumulativeMonotoneInUnitInterval) {
    TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
    auto result = spca.select();

    double prev = 0.0;
    for (Eigen::Index k = 0; k < M; ++k) {
        double cum = result.cumulative_ev(k);
        EXPECT_GE(cum, prev - 1e-12); // non-decreasing (with numerical tolerance)
        EXPECT_GE(cum, 0.0);
        EXPECT_LE(cum, 1.0 + 1e-12);
        prev = cum;
    }
}


// ========================================================================================
// Score / loading consistency
// ========================================================================================

/** @brief ActiveSet mode: Z == X_centered * V and V columns are unit-norm. */
TEST_F(TRexSPCATest, Scores_MatchCenteredXTimesV_ActiveSet) {
    expectScoresMatchCenteredXTimesV(SPCAMode::ACTIVE_SET);
}


/** @brief Thresholded mode: Z == X_centered * V even when the thresholded
 *         support of V differs from the GVS active set (regression test for
 *         the score loop formerly iterating over the active set only). */
TEST_F(TRexSPCATest, Scores_MatchCenteredXTimesV_Thresholded) {
    expectScoresMatchCenteredXTimesV(SPCAMode::THRESHOLDED);
}


// ========================================================================================
// Data integrity
// ========================================================================================

/** @brief X is restored to its original values after select() returns. */
TEST_F(TRexSPCATest, DataIntegrity_XRestoredAfterSelect) {
    Eigen::MatrixXd X_copy = X;
    TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
    spca.select();

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was not restored after TRexSPCASelector::select().";
}


/** @brief The constructor does not mutate X — centering only happens inside select(). */
TEST_F(TRexSPCATest, DataIntegrity_ConstructorDoesNotMutateX) {
    Eigen::MatrixXd X_copy = X;
    { TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed); }

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was mutated by TRexSPCASelector constructor or destructor without select().";
}


/** @brief Destructor does not mutate X when select() was never called. */
TEST_F(TRexSPCATest, DataIntegrity_XUnchangedOnDestructionWithoutSelect) {
    Eigen::MatrixXd X_copy = X;

    {
        TRexSPCASelector spca(X_map, M, 0.2, fast_ctrl(), seed);
        // Destructor fires here; X_is_centered_ = false, so no inverse_transform.
    }

    EXPECT_TRUE(X.isApprox(X_copy, 1e-12))
        << "X was incorrectly modified by TRexSPCASelector destructor.";
}


// ========================================================================================
// GVS type recorded in result
// ========================================================================================

// ========================================================================================
// Exception handling
// ========================================================================================

/** @brief M > min(n, p) triggers std::invalid_argument in the constructor. */
TEST_F(TRexSPCATest, Exceptions_InvalidM) {
    EXPECT_THROW(
        (TRexSPCASelector(X_map, p + 1, 0.2, fast_ctrl(), seed)),
        std::invalid_argument
    );
}

/** @brief TRex-SPCA is defined on the elastic net: an IEN request is rejected
 *         at construction instead of being run as EN. */
TEST_F(TRexSPCATest, Exceptions_IENRejected) {
    TRexSPCAControl ctrl = fast_ctrl();
    ctrl.gvs_ctrl.gvs_type = tg::GVSType::IEN;
    EXPECT_THROW((TRexSPCASelector(X_map, M, 0.2, ctrl, seed)), std::invalid_argument);

    ctrl.gvs_ctrl.gvs_type = tg::GVSType::EN;
    EXPECT_NO_THROW((TRexSPCASelector(X_map, M, 0.2, ctrl, seed)));
}

/** @brief M below 1 is refused at construction, before any allocation. */
TEST_F(TRexSPCATest, Exceptions_MBelowOne) {
    EXPECT_THROW((TRexSPCASelector(X_map, 0, 0.2, fast_ctrl(), seed)), std::invalid_argument);
    EXPECT_THROW((TRexSPCASelector(X_map, -1, 0.2, fast_ctrl(), seed)), std::invalid_argument);
    EXPECT_NO_THROW((TRexSPCASelector(X_map, 1, 0.2, fast_ctrl(), seed)));
    EXPECT_NO_THROW((TRexSPCASelector(X_map, p, 0.2, fast_ctrl(), seed)));
}

/** @brief The target FDR is checked at construction: [0, 1] is accepted,
 *         values outside and NaN are refused. */
TEST_F(TRexSPCATest, Exceptions_TargetFDROutsideUnitInterval) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const double bad : {-0.1, 1.5, nan}) {
        EXPECT_THROW((TRexSPCASelector(X_map, M, bad, fast_ctrl(), seed)), std::invalid_argument)
            << "tFDR = " << bad;
    }
    for (const double good : {0.0, 0.2, 1.0}) {
        EXPECT_NO_THROW((TRexSPCASelector(X_map, M, good, fast_ctrl(), seed))) << "tFDR = " << good;
    }
}

/** @brief The rules of the informed elastic net for lambda_2 are refused at
 *         construction, and the four rules of the elastic net are accepted. */
TEST_F(TRexSPCATest, Exceptions_IENRulesRejectedAtConstruction) {
    using Rule = tg::GVSLambda2Method;
    for (const Rule rule : {Rule::CV_1SE_IEN_CCD, Rule::CV_MIN_IEN_CCD,
                            Rule::CV_1SE_TIK_SVD, Rule::CV_MIN_TIK_SVD}) {
        TRexSPCAControl ctrl = fast_ctrl();
        ctrl.gvs_ctrl.lambda2_method = rule;
        EXPECT_THROW((TRexSPCASelector(X_map, M, 0.2, ctrl, seed)), std::invalid_argument);
    }
    for (const Rule rule : {Rule::CV_1SE_CCD, Rule::CV_MIN_CCD,
                            Rule::CV_1SE_SVD, Rule::CV_MIN_SVD}) {
        TRexSPCAControl ctrl = fast_ctrl();
        ctrl.gvs_ctrl.lambda2_method = rule;
        EXPECT_NO_THROW((TRexSPCASelector(X_map, M, 0.2, ctrl, seed)));
    }
}

/** @brief The permutation strategies are refused at construction, and the
 *         four strategies of the T-Rex+GVS runs complete a selection. */
TEST_F(TRexSPCATest, LLoopStrategies_FourRunAndPermutationsAreRefused) {
    using Strategy = tc::LLoopStrategy;
    for (const Strategy strategy : {Strategy::PERMUTATION, Strategy::PERMUTATION_SEEDED}) {
        TRexSPCAControl ctrl = fast_ctrl();
        ctrl.gvs_ctrl.trex_ctrl.lloop_strategy = strategy;
        EXPECT_THROW((TRexSPCASelector(X_map, M, 0.2, ctrl, seed)), std::invalid_argument);
    }
    for (const Strategy strategy : {Strategy::STANDARD, Strategy::HCONCAT,
                                    Strategy::SKIPL, Strategy::SEEDED}) {
        TRexSPCAControl ctrl = fast_ctrl();
        ctrl.gvs_ctrl.trex_ctrl.lloop_strategy = strategy;
        TRexSPCASelector spca(X_map, M, 0.2, ctrl, seed);
        EXPECT_NO_THROW(spca.select());
        EXPECT_EQ(spca.getResult().active_sets.size(), static_cast<std::size_t>(M));
    }
}


// ========================================================================================
// Documented behaviour of the control and the result
// ========================================================================================

/** @brief The nested gvs_ctrl.en_solver and gvs_ctrl.trex_ctrl.solver_type
 *         have no effect: the selector replaces them by en_solver. */
TEST_F(TRexSPCATest, Control_NestedSolverFieldsHaveNoEffect) {
    using SolverType = decltype(tc::TRexControl::solver_type);

    auto run = [](const TRexSPCAControl& ctrl) {
        Eigen::MatrixXd Xf = plantedFactorX();
        Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), n_f, p_f);
        TRexSPCASelector spca(Xf_map, M, 0.1, ctrl, seed);
        return spca.select();
    };

    TRexSPCAControl plain = fast_ctrl();
    plain.gvs_ctrl.lambda_2 = 1.0;

    // A nested solver type outside the elastic-net family would end a
    // T-Rex+GVS run with an error.
    TRexSPCAControl nested = plain;
    nested.gvs_ctrl.en_solver = tg::ENSolverType::TENET_AUG;
    nested.gvs_ctrl.trex_ctrl.solver_type = SolverType::TLASSO;

    const TRexSPCAResult a = run(plain);
    const TRexSPCAResult b = run(nested);
    EXPECT_EQ(a.active_sets, b.active_sets);
    EXPECT_EQ(a.V, b.V);
}

/** @brief cumulative_ev is the running sum of adjusted_ev over the total
 *         variance of the centered X, a fraction and not a percentage. */
TEST_F(TRexSPCATest, ExplainedVariance_CumulativeIsFractionOfTotalVariance) {
    Eigen::MatrixXd Xf = plantedFactorX();
    const Eigen::RowVectorXd mu = Xf.colwise().mean();
    const double total_variance =
        (Xf.rowwise() - mu).squaredNorm() / static_cast<double>(n_f - 1);

    Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), n_f, p_f);
    TRexSPCASelector spca(Xf_map, M, 0.1, fast_ctrl(), seed);
    const TRexSPCAResult result = spca.select();

    double running = 0.0;
    for (Eigen::Index m = 0; m < M; ++m) {
        running += result.adjusted_ev(m);
        EXPECT_NEAR(result.cumulative_ev(m), running / total_variance, 1e-12);
        EXPECT_LE(result.cumulative_ev(m), 1.0 + 1e-12);
    }
}

/** @brief Thresholded mode keeps as many loadings as variables were selected. */
TEST_F(TRexSPCATest, Loadings_ThresholdedKeepsAsManyEntriesAsSelected) {
    Eigen::MatrixXd Xf = plantedFactorX();
    Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), n_f, p_f);
    TRexSPCASelector spca(Xf_map, M, 0.1, fast_ctrl(SPCAMode::THRESHOLDED), seed);
    const TRexSPCAResult result = spca.select();

    for (Eigen::Index m = 0; m < M; ++m) {
        const auto nonzero = (result.V.col(m).array() != 0.0).count();
        EXPECT_EQ(static_cast<std::size_t>(nonzero),
                  result.active_sets[static_cast<std::size_t>(m)].size());
    }
}

/** @brief A component without a selected variable has a zero column in V and
 *         Z and the explained variance 0, and every other loading has unit
 *         length. The rule holds for every component, whatever is selected.
 */
TEST_F(TRexSPCATest, Result_ComponentWithoutSelectionIsZero) {
    for (const double tFDR : {0.0, 0.05, 0.2}) {
        Eigen::MatrixXd Xf = plantedFactorX();
        Eigen::Map<Eigen::MatrixXd> Xf_map(Xf.data(), n_f, p_f);
        TRexSPCASelector spca(Xf_map, 3, tFDR, fast_ctrl(), seed);
        const TRexSPCAResult result = spca.select();

        double previous = 0.0;
        for (Eigen::Index m = 0; m < 3; ++m) {
            if (result.active_sets[static_cast<std::size_t>(m)].empty()) {
                EXPECT_EQ(result.V.col(m).norm(), 0.0) << "tFDR = " << tFDR;
                EXPECT_EQ(result.Z.col(m).norm(), 0.0) << "tFDR = " << tFDR;
                EXPECT_EQ(result.adjusted_ev(m), 0.0) << "tFDR = " << tFDR;
                EXPECT_EQ(result.cumulative_ev(m), previous) << "tFDR = " << tFDR;
            } else {
                EXPECT_NEAR(result.V.col(m).norm(), 1.0, 1e-12) << "tFDR = " << tFDR;
            }
            previous = result.cumulative_ev(m);
        }
    }
}

// ========================================================================================
} /* End of namespace trex::test::trex_selector_methods::ad::trex_spca */
