// ========================================================================================
/**
 * @file test_trex_restore.cpp
 *
 * @brief TRexSelector::restore(): give X back before select().
 *
 * @details The selector normalizes X in place at construction and restores
 *          it when select() returns (and in the destructor). restore() hands
 *          X back right away for a selector that will not run; it is
 *          idempotent, a no-op after select(), and a restored selector
 *          refuses to run. The screening selector, which overrides select(),
 *          honours the same rule.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// project includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>
#include <tsolvers/ad/solver_utils/buffer_registry.hpp>
#include <utils/datageneration/utils_datagen.hpp>

// Eigen includes
#include <Eigen/Dense>

// ========================================================================================

namespace trex::test::trex_selector_methods::restore {

namespace datagen = trex::utils::datageneration::datagen;
namespace bufreg = trex::tsolvers::ad::solver_utils::buffer_registry;
using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::trex_selector_methods::ad::trex_gvs::GVSType;
using trex::trex_selector_methods::ad::trex_gvs::TRexGVSControlParameter;
using trex::trex_selector_methods::ad::trex_gvs::TRexGVSSelector;
using trex::trex_selector_methods::ad::trex_screening::ScreenTRexSelector;
using trex::trex_selector_methods::ad::trex_screening::ScreenTRexControlParameter;
using trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;

namespace {

struct Data {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Data() {
        const Eigen::Index n = 50, p = 40;
        std::vector<std::size_t> sup = {2, 9, 17};
        std::vector<double> coefs(sup.size(), 2.0);
        datagen::SyntheticData data(n, p, sup, coefs, /*snr=*/3.0, /*seed=*/4242,
                                    datagen::predictor_policy::Normal(),
                                    datagen::noisegen::noise_policy::Normal());
        X = data.getX();
        y = data.getY();
    }
};

TRexControlParameter smallControl() {
    TRexControlParameter ctrl;
    ctrl.K = 5;
    ctrl.max_dummy_multiplier = 2;
    ctrl.solver_type = SolverTypeForTRex::TLARS;
    return ctrl;
}

} // namespace

// ========================================================================================

/** @brief restore() before select() gives X back and releases the claim;
 *         it is idempotent; select() afterwards throws. */
TEST(TRexRestore, RestoreBeforeSelectGivesXBack) {
    Data d;
    const Eigen::MatrixXd X0 = d.X;
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexSelector sel(Xm, ym, /*tFDR=*/0.2, smallControl(), /*seed=*/3, /*verbose=*/false);
    EXPECT_FALSE(d.X.isApprox(X0, 1e-12));           // normalized in place
    EXPECT_TRUE(bufreg::isClaimed(d.X.data()));

    sel.restore();
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));
    EXPECT_FALSE(bufreg::isClaimed(d.X.data()));

    sel.restore();                                    // idempotent
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));

    EXPECT_THROW(sel.select(), std::runtime_error);   // a restored selector refuses to run
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));             // and left X alone

    // The buffer is free for a new selector, which runs normally
    TRexSelector fresh(Xm, ym, 0.2, smallControl(), 3, false);
    EXPECT_NO_THROW(fresh.select());
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));
}


/** @brief After select() the selector has restored X itself; restore() is
 *         a no-op and does not move the data. */
TEST(TRexRestore, RestoreAfterSelectIsNoOp) {
    Data d;
    const Eigen::MatrixXd X0 = d.X;
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexSelector sel(Xm, ym, 0.2, smallControl(), 3, false);
    sel.select();
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));
    const Eigen::MatrixXd X_after = d.X;
    sel.restore();
    EXPECT_EQ((d.X - X_after).cwiseAbs().maxCoeff(), 0.0);
}


/** @brief The screening selector overrides select() and honours the rule. */
TEST(TRexRestore, ScreeningSelectorHonoursRestore) {
    Data d;
    const Eigen::MatrixXd X0 = d.X;
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    ScreenTRexControlParameter sctrl;
    ScreenTRexSelector sel(Xm, ym, sctrl, /*seed=*/7);
    sel.restore();
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));
    EXPECT_THROW(sel.select(), std::runtime_error);
}



/** @brief A selector runs once. select() gives X back when it ends, so a
 *         second call throws instead of running on a design that is no
 *         longer standardized. The results of the run stay readable. */
TEST(TRexRestore, SecondSelectThrows) {
    Data d;
    const Eigen::MatrixXd X0 = d.X;
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexSelector sel(Xm, ym, 0.2, smallControl(), 3, false);
    sel.select();
    const std::vector<std::size_t> first = sel.getSelectedIndices();
    const Eigen::MatrixXd X_after = d.X;

    EXPECT_THROW(sel.select(), std::runtime_error);
    EXPECT_THROW(sel.select(), std::runtime_error);            // and every later call
    EXPECT_EQ(sel.getSelectedIndices(), first);                // the results stay
    EXPECT_EQ((d.X - X_after).cwiseAbs().maxCoeff(), 0.0);     // X is left alone
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));

    // A new selector on the same matrix repeats the run.
    TRexSelector again(Xm, ym, 0.2, smallControl(), 3, false);
    again.select();
    EXPECT_EQ(again.getSelectedIndices(), first);
}


/** @brief The rule holds under ScalingMode::NONE as well, where the selector
 *         writes nothing into X. */
TEST(TRexRestore, SecondSelectThrowsUnderScalingNone) {
    Data d;
    // A design the user standardized: centered columns of unit L2 norm.
    d.X.rowwise() -= d.X.colwise().mean();
    d.X.colwise().normalize();
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    TRexControlParameter ctrl = smallControl();
    ctrl.scaling_mode = trex::trex_selector_methods::ad::utils::data_normalizer::ScalingMode::NONE;
    TRexSelector sel(Xm, ym, 0.2, ctrl, 3, false);
    EXPECT_NO_THROW(sel.select());
    EXPECT_THROW(sel.select(), std::runtime_error);
}


/** @brief The screening selector has its own select() and follows the rule. */
TEST(TRexRestore, ScreeningSelectorRunsOnce) {
    Data d;
    const Eigen::MatrixXd X0 = d.X;
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());

    ScreenTRexControlParameter sctrl;
    ScreenTRexSelector sel(Xm, ym, sctrl, /*seed=*/7);
    EXPECT_NO_THROW(sel.select());
    EXPECT_THROW(sel.select(), std::runtime_error);
    EXPECT_TRUE(d.X.isApprox(X0, 1e-12));
}



/** @brief A select() that ends with an error gives X back at once and
 *         retires the selector. The error is an unconverged solve: the
 *         coordinate-descent solver of the informed elastic net at a ridge
 *         weight beyond its sweep budget (one group per variable). */
TEST(TRexRestore, SelectEndingWithAnErrorGivesXBack) {
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
    for (Eigen::Index i = 0; i < n; ++i) { y(i) = X.row(i).segment(20, 10).sum() + N01(); }
    const Eigen::MatrixXd X_before = X;

    Eigen::Map<Eigen::MatrixXd> Xm(X.data(), n, p);
    Eigen::Map<Eigen::VectorXd> ym(y.data(), n);
    TRexGVSControlParameter ctrl;
    ctrl.gvs_type = GVSType::IEN;
    ctrl.lambda_2 = 3000.0;
    ctrl.prior_groups.resize(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) { ctrl.prior_groups[static_cast<std::size_t>(j)] = j; }
    ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TCDIENET;
    ctrl.trex_ctrl.K = 10;
    TRexGVSSelector sel(Xm, ym, /*tFDR=*/0.1, ctrl, /*seed=*/4, /*verbose=*/false);
    EXPECT_GT((X - X_before).cwiseAbs().maxCoeff(), 1e-3);   // standardized in place

    try {
        sel.select();
        FAIL() << "select() returned from an unconverged path";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("sweep budget"), std::string::npos) << e.what();
    }
    EXPECT_LT((X - X_before).cwiseAbs().maxCoeff(), 1e-10);  // given back
    EXPECT_THROW(sel.select(), std::runtime_error);          // retired
}

} // namespace trex::test::trex_selector_methods::restore
