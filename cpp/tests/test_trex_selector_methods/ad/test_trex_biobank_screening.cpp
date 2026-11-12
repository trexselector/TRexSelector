// =============================================================================
// test_trex_biobank_screening.cpp
// =============================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
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
