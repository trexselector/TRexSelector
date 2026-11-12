// =============================================================================
// test_trex_biobank_screening.cpp
// =============================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
#include <limits>
#include <random>
#include <stdexcept>

// project trex includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_biobank_screening.hpp>

// ===================================================================================

// Embed into test namespace
namespace trex::test::trex_selector_methods::ad::trex_biobank_screening {

// Namespace aliases for convenience
using namespace trex::trex_selector_methods::ad::trex_core;
using namespace trex::trex_selector_methods::ad::trex_screening;
using namespace trex::trex_selector_methods::ad::trex_biobank_screening;


// ===================================================================================


class TRexBiobankScreeningTest : public ::testing::Test {
protected:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::MatrixXd Y_matrix;

    Eigen::Map<Eigen::MatrixXd> X_map;
    Eigen::Map<Eigen::VectorXd> y_map;
    Eigen::Map<Eigen::MatrixXd> Y_matrix_map;

    TRexBiobankScreeningTest()
        : X(Eigen::MatrixXd::Random(50, 20)),
          y(Eigen::VectorXd::Random(50)),
          Y_matrix(Eigen::MatrixXd::Random(50, 3)),
          X_map(X.data(), X.rows(), X.cols()),
          y_map(y.data(), y.size()),
          Y_matrix_map(Y_matrix.data(), Y_matrix.rows(), Y_matrix.cols()) {}
};


/** @brief Ensure invalid phenotypic calls map to runtime_errors */
TEST_F(TRexBiobankScreeningTest, Validation_ThrowsOnMismatchedPhenotypeCalls) {
    BiobankScreenTRexControl control;

    // Single Phenotype constructor
    BiobankScreenTRex single_pheno_selector(X_map, y_map, control);
    EXPECT_THROW(single_pheno_selector.screenPhenotypes(), std::runtime_error);

    // Multiple Phenotype constructor
    BiobankScreenTRex multi_pheno_selector(X_map, Y_matrix_map, control);
    EXPECT_THROW(multi_pheno_selector.screenPhenotype(), std::runtime_error);
}


/** @brief Standard single execution works cleanly */
TEST_F(TRexBiobankScreeningTest, Execution_SinglePhenotypeDoesNotThrow) {
    BiobankScreenTRexControl control;
    control.trex_screen_ctrl.trex_method = ScreenTRexMethod::TREX;

    // As identified, we must override internal defaults to standard!
    control.trex_screen_ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;
    control.trex_screen_ctrl.trex_ctrl.K = 5;

    BiobankScreenTRex selector(X_map, y_map, control);
    EXPECT_NO_THROW({
        auto result = selector.screenPhenotype();
    });
}


namespace {

/** @brief Probe exposing the seed derivation for testing. */
class BiobankScreenTRexProbe : public BiobankScreenTRex {
public:
    using BiobankScreenTRex::BiobankScreenTRex;
    using BiobankScreenTRex::phenotypeSeed;
};

} // namespace


/** @brief Phenotype i runs with seed + i. Phenotype 0 keeps the seed, a
 *  negative seed is passed on, and the sum stays a non-negative int. */
TEST_F(TRexBiobankScreeningTest, Seed_PhenotypeSeedDerivation) {
    BiobankScreenTRexControl control;

    BiobankScreenTRexProbe seeded(X_map, Y_matrix_map, control, 11, false);
    EXPECT_EQ(seeded.phenotypeSeed(0), 11);
    EXPECT_EQ(seeded.phenotypeSeed(1), 12);
    EXPECT_EQ(seeded.phenotypeSeed(250), 261);

    BiobankScreenTRexProbe random(X_map, Y_matrix_map, control, -1, false);
    EXPECT_EQ(random.phenotypeSeed(0), -1);
    EXPECT_EQ(random.phenotypeSeed(7), -1);

    const int largest = std::numeric_limits<int>::max();
    BiobankScreenTRexProbe wrapped(X_map, Y_matrix_map, control, largest, false);
    EXPECT_EQ(wrapped.phenotypeSeed(0), largest);
    EXPECT_EQ(wrapped.phenotypeSeed(1), 0);
    EXPECT_EQ(wrapped.phenotypeSeed(3), 2);
}


/** @brief The phenotypes draw independent dummies, and each one is reproduced
 *  by a single-phenotype run with its seed. Twelve copies of one noise
 *  response make the dummy draw visible: the screening of such a response
 *  selects something in about half of the draws. */
TEST(TRexBiobankScreeningSeedTest, Execution_PhenotypesUseIndependentSeeds) {
    std::mt19937 rng(31);
    std::normal_distribution<double> N01(0.0, 1.0);
    const Eigen::Index n = 80;
    const Eigen::Index p = 40;
    const Eigen::Index copies = 12;
    Eigen::MatrixXd X(n, p);
    Eigen::VectorXd y(n);
    for (Eigen::Index j = 0; j < p; ++j) {
        for (Eigen::Index i = 0; i < n; ++i) { X(i, j) = N01(rng); }
    }
    for (Eigen::Index i = 0; i < n; ++i) { y(i) = N01(rng); }

    BiobankScreenTRexControl control;
    control.trex_screen_ctrl.R_boot = 100;
    control.trex_screen_ctrl.trex_ctrl.K = 10;
    const int seed = 11;

    auto single = [&](int run_seed) {
        Eigen::MatrixXd X_run = X;
        Eigen::VectorXd y_run = y;
        Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());
        BiobankScreenTRex selector(X_map, y_map, control, run_seed, false);
        return selector.screenPhenotype();
    };

    Eigen::MatrixXd X_run = X;
    Eigen::MatrixXd Y = y.replicate(1, copies);
    Eigen::Map<Eigen::MatrixXd> X_map(X_run.data(), X_run.rows(), X_run.cols());
    Eigen::Map<Eigen::MatrixXd> Y_map(Y.data(), Y.rows(), Y.cols());
    BiobankScreenTRex selector(X_map, Y_map, control, seed, false);
    const auto results = selector.screenPhenotypes();
    ASSERT_EQ(results.size(), static_cast<std::size_t>(copies));

    std::size_t differing = 0;
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto expected = single(seed + static_cast<int>(i));
        EXPECT_EQ(results[i].selected_indices, expected.selected_indices) << "phenotype " << i;
        EXPECT_EQ(results[i].selected_indices_screen_ordinary, expected.selected_indices_screen_ordinary);
        EXPECT_EQ(results[i].selected_indices_screen_bootstrap, expected.selected_indices_screen_bootstrap);
        EXPECT_DOUBLE_EQ(results[i].estimated_FDR_screen_ordinary, expected.estimated_FDR_screen_ordinary);
        EXPECT_DOUBLE_EQ(results[i].estimated_FDR_screen_bootstrap, expected.estimated_FDR_screen_bootstrap);
        EXPECT_EQ(results[i].method_used, expected.method_used);

        if (results[i].selected_indices_screen_ordinary != results[0].selected_indices_screen_ordinary) {
            ++differing;
        }
    }
    EXPECT_GT(differing, 0U)
        << "With one shared dummy draw the twelve copies would give twelve identical screening results.";
}


/** @brief An empty screening result reports an estimated FDR of 1, which no
 *  acceptance window with an upper bound below 1 holds. With a lower bound of
 *  0 the workflow therefore falls back to the T-Rex selector instead of
 *  accepting the empty result. */
TEST(TRexBiobankScreeningEstimateTest, Execution_EmptyScreeningFallsBack) {
    std::srand(23);
    const Eigen::Index n = 80;
    const Eigen::Index p = 40;
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, p);
    const Eigen::VectorXd noise = Eigen::VectorXd::Random(n);

    // Residual of the noise after the projection onto [1, X]: every variable
    // stays behind the first dummy, so the screening selects nothing.
    Eigen::MatrixXd basis(n, p + 1);
    basis << Eigen::VectorXd::Ones(n), X;
    Eigen::VectorXd y = noise - basis * basis.colPivHouseholderQr().solve(noise);

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    BiobankScreenTRexControl control;
    control.lower_bound_FDR = 0.0;
    control.upper_bound_FDR = 0.15;
    control.trex_screen_ctrl.R_boot = 100;
    control.trex_screen_ctrl.trex_ctrl.K = 5;

    BiobankScreenTRex selector(X_map, y_map, control, 5, false);
    const auto result = selector.screenPhenotype();

    EXPECT_TRUE(result.selected_indices_screen_ordinary.empty());
    EXPECT_DOUBLE_EQ(result.estimated_FDR_screen_ordinary, 1.0);
    EXPECT_DOUBLE_EQ(result.estimated_FDR_screen_bootstrap, 1.0);
    EXPECT_TRUE(result.used_fallback_trex);
    EXPECT_EQ(result.method_used, "T-Rex (fallback)");

    // A fallback result reports the target level, and the flag marks it.
    EXPECT_DOUBLE_EQ(result.estimated_FDR, control.target_FDR_trex);
}


/** @brief Execute fallback testing under impossible bounds */
TEST_F(TRexBiobankScreeningTest, Execution_FallbackLogicActivation) {
    BiobankScreenTRexControl control;
    // Force the decision logic to abandon Screen-TRex by setting impossible bounds:
    control.lower_bound_FDR = 0.99;
    control.upper_bound_FDR = 0.00;

    // Prevent the Screen-TRex components from crashing under SKIPL
    control.trex_screen_ctrl.trex_ctrl.lloop_strategy = LLoopStrategy::STANDARD;
    control.trex_screen_ctrl.trex_ctrl.K = 2; // Keep loop overhead low

    BiobankScreenTRex selector(X_map, y_map, control);
    auto result = selector.screenPhenotype();

    // Verify Algorithm 1 dumped us into the full T-Rex calculation tree
    EXPECT_EQ(result.method_used, "T-Rex (fallback)");
    EXPECT_TRUE(result.used_fallback_trex);
}

// ===================================================================================
} /* End of namespace trex::test::trex_selector_methods::ad::trex_biobank_screening */
