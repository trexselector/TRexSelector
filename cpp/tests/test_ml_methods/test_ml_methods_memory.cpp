// ==================================================================================
/**
 * @file test_ml_methods_memory.cpp
 *
 * @brief Peak-memory budget of the ML cores, in units of the design matrix
 *        (the allocation guard).
 *
 * @details Eigen allocates through malloc, not operator new, so a global
 *          operator-new counter cannot see it. The measure is the number of
 *          bytes that are live on the heap, read from the allocator
 *          (malloc_zone_statistics on macOS, mallinfo2 of glibc 2.33 or later
 *          on Linux). On every other platform this file holds one skipped
 *          test. A second thread reads the value without a pause while a
 *          call runs on a design of 32 MB, and the growth of its peak over
 *          the level before the call is the working memory of the call.
 *
 *          The peak depends on what the call allocates and not on how the
 *          allocator reuses freed blocks, so it is the same in every run.
 *          The sampler can miss a peak that is shorter than its sampling gap
 *          and can never report more than the peak. A busy machine therefore
 *          lowers the value at most and cannot make a test fail.
 *
 *          A full copy of X raises the value by one X. The budgets below are
 *          the CURRENT state of each core, tightened as copies are removed,
 *          and each lies below the measured peak plus one X, so a regression
 *          that adds a copy of X fails here. Set TREX_ML_MEMORY_VERBOSE=1 to
 *          print the measured growth per call.
 *
 *          Measured peaks in units of X, wide 500 x 8000, macOS with Eigen
 *          5.0 and Ubuntu with Eigen 3.4: ridge CV 2.84 and 4.41, elastic net
 *          1.09 and 1.09, elastic net CV 1.34 and 1.09, SVD 0.33 and 0.22,
 *          ridge 0.27 and 0.08, PCA 0.34 and 0.22, scaler 0.006 and 0.004,
 *          low-memory ridge CV 0.18 and 0.19, low-memory elastic net 0.093
 *          and 0.091, low-memory elastic net CV 0.094 and 0.091.
 */
// ==================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <thread>

// heap statistics of the platform
#if defined(__APPLE__)
#include <malloc/malloc.h>
#define TREX_TEST_HEAP_STATISTICS 1
#elif defined(__GLIBC__) && defined(__GLIBC_PREREQ)
#if __GLIBC_PREREQ(2, 33)
#include <malloc.h>
#define TREX_TEST_HEAP_STATISTICS 1
#endif
#endif

// Eigen includes
#include <Eigen/Dense>

// project ml_methods includes
#include <ml_methods/model_selection/enet_cv_ccd.hpp>
#include <ml_methods/model_selection/ridge_cv_svd.hpp>
#include <ml_methods/svd/svd.hpp>
#include <ml_methods/pca/pca.hpp>
#include <ml_methods/ridge_regression/ridge.hpp>
#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>

// ==================================================================================

namespace trex::test::ml_methods::memory {

using trex::ml_methods::model_selection::enet_cv_ccd;
using trex::ml_methods::model_selection::enet_gaussian;
using trex::ml_methods::model_selection::ridge_cv_svd;
using trex::ml_methods::pca::PCA;
using trex::ml_methods::ridge::RidgeSolver;
using trex::ml_methods::scaler_methods::LpNormScaler;
using trex::ml_methods::svd::SVDSolver;

#if !defined(TREX_TEST_HEAP_STATISTICS)

/** @brief The measure needs the heap statistics of macOS or of glibc 2.33 or later. */
TEST(MLMemory, SkippedWithoutHeapStatistics) {
    GTEST_SKIP() << "live heap bytes: malloc_zone_statistics (macOS) or mallinfo2 (glibc 2.33 or later) only";
}

#else

namespace {

/** @brief Bytes that are allocated on the heap at this moment. */
double live_heap_bytes() {
#if defined(__APPLE__)
    malloc_statistics_t st {};
    malloc_zone_statistics(nullptr, &st);
    return static_cast<double>(st.size_in_use);
#else
    const struct mallinfo2 mi = mallinfo2();
    return static_cast<double>(mi.uordblks) + static_cast<double>(mi.hblkhd);
#endif
}

Eigen::MatrixXd random_matrix(Eigen::Index n, Eigen::Index p, unsigned seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> dist(0.0, 1.0);
    return Eigen::MatrixXd::NullaryExpr(
        n, p, [&](Eigen::Index, Eigen::Index) { return dist(gen); });
}

/** @brief The wide regime: n = 500, p = 8000 -> X is 32 MB. */
constexpr Eigen::Index kN = 500;
constexpr Eigen::Index kP = 8000;

struct Design {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    double bytes;
    Design() : X(random_matrix(kN, kP, 5u)), y(X.col(0) * 2.0 - X.col(5) + random_matrix(kN, 1, 6u)),
               bytes(static_cast<double>(kN) * static_cast<double>(kP) * 8.0) {}
};

/**
 * @brief Runs f and returns the peak of the live heap bytes over the call,
 *        above the level before it, in units of X. Prints it if asked.
 *
 * @details A second thread reads the live heap bytes without a pause while f
 *          runs. It reports the peak or less, never more.
 */
template <class F>
double growth_in_X(const Design& d, const char* what, F&& f) {
    const double before = live_heap_bytes();
    double peak = before;
    std::atomic<bool> stop{false};
    std::thread sampler([&] {
        while (!stop.load(std::memory_order_relaxed)) {
            peak = std::max(peak, live_heap_bytes());
            std::this_thread::yield();
        }
    });
    try {
        f();
    } catch (...) {
        stop.store(true);
        sampler.join();
        throw;
    }
    stop.store(true);
    sampler.join();
    const double units = (peak - before) / d.bytes;
    if (const char* v = std::getenv("TREX_ML_MEMORY_VERBOSE"); v != nullptr && *v == '1') {
        std::cout << "[  MEMORY  ] " << what << ": peak growth " << units << " x X ("
                  << (peak - before) / 1e6 << " MB)\n";
    }
    return units;
}

} // namespace

// ==================================================================================
// Budgets: the current state of each core
// ==================================================================================

/** @brief RidgeCV: per fold X_train + standardized copy + Jacobi internals. */
TEST(MLMemory, RidgeCV) {
    Design d;
    const double g = growth_in_X(d, "ridge_cv_svd::fit", [&] {
        ridge_cv_svd cv;
        cv.fit(d.X, d.y, 5, 20, 1000.0, 1u);
    });
    // The Jacobi SVD keeps one internal copy of the fold under Eigen 5.0 (the
    // QR preconditioner) and two under Eigen 3.4 (plus its scaled copy):
    // measured 2.84 X on macOS / Eigen 5.0 and 4.41 X on Ubuntu / Eigen 3.4.
    const double budget = EIGEN_VERSION_AT_LEAST(5, 0, 0) ? 3.3 : 4.9;
    EXPECT_LT(g, budget) << "ridge CV allocated more than its budget of " << budget << " x X";
}

/** @brief ElasticNet: the centred member copy Xs_ plus p x K paths. */
TEST(MLMemory, ElasticNetPath) {
    Design d;
    const double g = growth_in_X(d, "enet_gaussian::fit", [&] {
        enet_gaussian en;
        en.fit(d.X, d.y, 1.0, 20, -1.0, true, true, false, 1000, 1e-5);
    });
    EXPECT_LT(g, 1.6) << "elastic net allocated more than its budget of 1.6 x X";
}

/** @brief RidgeCV, kernel engine: n x n kernels and column blocks only. */
TEST(MLMemory, RidgeCVLowMemory) {
    Design d;
    const double g = growth_in_X(d, "ridge_cv_svd::fit (low_memory)", [&] {
        ridge_cv_svd cv(/*low_memory=*/true);
        cv.fit(d.X, d.y, 5, 20, 1000.0, 1u);
    });
    // O(n^2) kernels plus the eigen solver's workspaces and the column-block
    // buffers: about 6 MB here (0.18 X of a 32 MB design), independent of p.
    EXPECT_LT(g, 0.5) << "low-memory ridge CV allocated more than its budget of 0.5 x X";
}

/** @brief ElasticNet, raw-X engine: no allocation of the size of X. */
TEST(MLMemory, ElasticNetPathLowMemory) {
    Design d;
    const double g = growth_in_X(d, "enet_gaussian::fit (low_memory)", [&] {
        enet_gaussian en(/*low_memory=*/true);
        en.fit(d.X, d.y, 1.0, 20, -1.0, true, true, false, 1000, 1e-5);
    });
    EXPECT_LT(g, 0.2) << "low-memory elastic net allocated more than its budget of 0.2 x X";
}

/** @brief ElasticNetCV, raw-X engine: masks and row-wise predictions only. */
TEST(MLMemory, ElasticNetCVLowMemory) {
    Design d;
    const double g = growth_in_X(d, "enet_cv_ccd::fit (low_memory)", [&] {
        enet_cv_ccd cv(/*low_memory=*/true);
        cv.fit(d.X, d.y, 0.5, 5, 20, -1.0, 2u, true, true, 1000, 1e-5);
    });
    // One fit needs about 0.1 X (the p x K paths, 0.093 X in ElasticNetPathLowMemory).
    // The six fits of a 5-fold CV run one after another and hold the paths of one
    // fit at a time: 0.094 X.
    EXPECT_LT(g, 0.3) << "low-memory elastic net CV allocated more than its budget of 0.3 x X";
}

/** @brief ElasticNetCV: full-data Xs_, then per fold X_tr + X_te + the fold's Xs_. */
TEST(MLMemory, ElasticNetCV) {
    Design d;
    const double g = growth_in_X(d, "enet_cv_ccd::fit", [&] {
        enet_cv_ccd cv;
        cv.fit(d.X, d.y, 0.5, 5, 20, -1.0, 2u, true, true, 1000, 1e-5);
    });
    EXPECT_LT(g, 1.6) << "elastic net CV allocated more than its budget of 1.6 x X";
}

/** @brief SVD, wide: the Gram path forms n x n and reads X blockwise. */
TEST(MLMemory, SVDWide) {
    Design d;
    const double g = growth_in_X(d, "SVDSolver::compute (wide, Gram path)", [&] {
        auto r = SVDSolver{}.compute(d.X, 5);
        volatile double sink = r.S.sum(); (void)sink;
    });
    EXPECT_LT(g, 0.5) << "wide SVD allocated more than its budget of 0.5 x X";
}

/** @brief Ridge, wide: the n x n kernel and its factorization only. */
TEST(MLMemory, RidgeWide) {
    Design d;
    const double g = growth_in_X(d, "RidgeSolver::solve (wide, dual)", [&] {
        Eigen::VectorXd b = RidgeSolver::solve(d.X, d.y, 1.0);
        volatile double sink = b.sum(); (void)sink;
    });
    EXPECT_LT(g, 0.5) << "wide ridge allocated more than its budget of 0.5 x X";
}

/** @brief PCA on a wide design: in place, then the Gram-path SVD. */
TEST(MLMemory, PCAWide) {
    Design d;
    const double g = growth_in_X(d, "PCA (wide)", [&] {
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        PCA pca(Xm, 5, true, true);
        auto r = pca.fit();
        volatile double sink = r.explained_variance.sum(); (void)sink;
        pca.restore(Xm);
    });
    EXPECT_LT(g, 0.5) << "wide PCA allocated more than its budget of 0.5 x X";
}

/** @brief Scalers: in place, no allocation of the size of X. */
TEST(MLMemory, Scaler) {
    Design d;
    const double g = growth_in_X(d, "LpNormScaler::fit_transform_inplace", [&] {
        Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
        LpNormScaler s;
        s.fit_transform_inplace(Xm);
        s.inverse_transform_inplace(Xm);
    });
    EXPECT_LT(g, 0.1) << "scaler allocated more than its budget of 0.1 x X";
}

#endif  // TREX_TEST_HEAP_STATISTICS

} // namespace trex::test::ml_methods::memory
