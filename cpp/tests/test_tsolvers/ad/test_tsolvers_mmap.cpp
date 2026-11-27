// ========================================================================================
// test_tsolvers_mmap.cpp
// ========================================================================================
/**
 * @file test_tsolvers_mmap.cpp
 *
 * @brief Unit tests verifying memory-mapping behavior interaction with solver preprocessing.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <filesystem>
#include <stdexcept>
#include <string>


// Eigen includes
#include <Eigen/Dense>

// tsolvers includes
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlasso_solver.hpp>
#include <tsolvers/ad/solver_options.hpp>

// utils includes
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::tsolvers::ad::linear_model {
namespace fs = std::filesystem;

using namespace trex::tsolvers::ad::linear_model::lars_based;
using namespace trex::utils::datageneration::datagen;

// ========================================================================================

/** @brief Test to verify that memory-mapped data is correctly mutated by the solver */
TEST(TSolverMmapTest, ValidateMapMutation) {
    SyntheticData data(
        100, 20, 20,
        {1, 2},
        {5.0, 2.0},
        1.0, 42, -1, -1,
        predictor_policy::Normal(),
        dummygen::DummyDistribution::Normal(),
        noisegen::noise_policy::Normal()
    );

    auto X = data.getX();
    auto D = data.getD();
    auto y = data.getY();

    // Modify to be un-centered
    X.array() += 10.0;

    Eigen::MatrixXd X_copy = X; // Keep copy for comparison

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), D.rows(), D.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TSolverLARS solver(X_map, D_map, y_map);

    // Initialization changes the original X memory (via Map) because of centering/scaling
    EXPECT_FALSE(X.isApprox(X_copy, 1e-10));
    EXPECT_NEAR(X.col(0).mean(), 0.0, 1e-10);
    EXPECT_NEAR(X.col(0).norm(), 1.0, 1e-10);
}

// ========================================================================================


/** @brief Test to verify that execution over disk-mapped data succeeds correctly */
TEST(TSolverMmapTest, ValidateDiskMappedDataExecution) {
    const std::string X_file = ::testing::TempDir() + "test_tsolvers_mmap_X.bin";
    const std::string D_file = ::testing::TempDir() + "test_tsolvers_mmap_D.bin";
    const std::string y_file = ::testing::TempDir() + "test_tsolvers_mmap_y.bin";

    // Clean up if files left over from previous aborted run
    if (fs::exists(X_file)) fs::remove(X_file);
    if (fs::exists(D_file)) fs::remove(D_file);
    if (fs::exists(y_file)) fs::remove(y_file);

    {
        // 1. Generate mapped data directly on disk using SyntheticDataMapped
        SyntheticDataMapped data(
            X_file, D_file, y_file,
            100, 20, 20,
            {1, 2}, {5.0, 2.0},
            1.0, 42
        );

        // 2. Retrieve Eigen::Map objects mapped over the file descriptors
        auto X_map = data.getX();
        auto D_map = data.getD();
        Eigen::Map<Eigen::VectorXd> y_map(data.getY().data(), data.rows());

        // 3. Feed the memory-mapped views to the solver
        TSolverLARS solver(X_map, D_map, y_map);
        solver.executeStep(5, true);

        // 4. Verify that execution over mapped memory succeeds correctly
        EXPECT_EQ(solver.getActiveDummyIndices().size(), 5);
        EXPECT_FALSE(solver.getActions().empty());
        EXPECT_GT(solver.getRSS().front(), solver.getRSS().back());

    } // the scoping block ensures SyntheticDataMapped instances destruct
      // and close their file handles before we attempt fs::remove()

    // 5. Clean up temporary files
    EXPECT_TRUE(fs::remove(X_file));
    EXPECT_TRUE(fs::remove(D_file));
    EXPECT_TRUE(fs::remove(y_file));
}



/** @brief A solver with normalize and intercept off only reads X and D, so
 *  both may be memory-mapped files opened ReadOnly: the run equals the
 *  in-memory run and leaves the files as they were. A writer is refused for
 *  such a file. */
TEST(TSolverMmapTest, ReaderRunsOnReadOnlyMemoryMappedFiles) {
    namespace mm = trex::utils::memmap;
    SyntheticData data(
        100, 20, 20,
        {1, 2},
        {5.0, 2.0},
        1.0, 42, -1, -1,
        predictor_policy::Normal(),
        dummygen::DummyDistribution::Normal(),
        noisegen::noise_policy::Normal()
    );
    // Centered columns of unit norm: the state a reader expects.
    auto standardize = [](Eigen::MatrixXd M) {
        M.rowwise() -= M.colwise().mean();
        M.colwise().normalize();
        return M;
    };
    const Eigen::MatrixXd X = standardize(data.getX());
    const Eigen::MatrixXd D = standardize(data.getD());
    Eigen::VectorXd y = data.getY();
    y.array() -= y.mean();

    const std::string X_file = (fs::path(::testing::TempDir()) / "tsolver_readonly_X.bin").string();
    const std::string D_file = (fs::path(::testing::TempDir()) / "tsolver_readonly_D.bin").string();
    fs::remove(X_file);
    fs::remove(D_file);
    {
        mm::MemoryMappedMatrix<double> X_writer(X_file, X.rows(), X.cols(), mm::AccessMode::ReadWrite);
        mm::MemoryMappedMatrix<double> D_writer(D_file, D.rows(), D.cols(), mm::AccessMode::ReadWrite);
        X_writer.getMap() = X;
        D_writer.getMap() = D;
    }

    trex::tsolvers::ad::TSolverOptions reader;
    reader.normalize = false;
    reader.intercept = false;
    trex::tsolvers::ad::TSolverOptions writer;  // normalize and intercept on
    EXPECT_FALSE(reader.writesDesign());
    EXPECT_TRUE(writer.writesDesign());

    // The in-memory reference
    Eigen::MatrixXd X_memory = X;
    Eigen::MatrixXd D_memory = D;
    Eigen::VectorXd y_memory = y;
    Eigen::Map<Eigen::MatrixXd> X_memory_map(X_memory.data(), X_memory.rows(), X_memory.cols());
    Eigen::Map<Eigen::MatrixXd> D_memory_map(D_memory.data(), D_memory.rows(), D_memory.cols());
    Eigen::Map<Eigen::VectorXd> y_memory_map(y_memory.data(), y_memory.size());
    TSolverLASSO reference(X_memory_map, D_memory_map, y_memory_map, reader);
    reference.executeStep(5, true);

    {
        mm::MemoryMappedMatrix<double> X_mapped(X_file, X.rows(), X.cols(), mm::AccessMode::ReadOnly);
        mm::MemoryMappedMatrix<double> D_mapped(D_file, D.rows(), D.cols(), mm::AccessMode::ReadOnly);

        EXPECT_THROW(X_mapped.getMapFor(writer.writesDesign()), std::runtime_error);

        auto X_map = X_mapped.getMapFor(reader.writesDesign());
        auto D_map = D_mapped.getMapFor(reader.writesDesign());
        Eigen::VectorXd y_run = y;
        Eigen::Map<Eigen::VectorXd> y_map(y_run.data(), y_run.size());
        {
            TSolverLASSO solver(X_map, D_map, y_map, reader);
            solver.executeStep(5, true);
            EXPECT_EQ(solver.getActions(), reference.getActions());
            EXPECT_TRUE(solver.getBeta(-1).isApprox(reference.getBeta(-1)));
        }

        // The files are as they were.
        EXPECT_TRUE(X_mapped.getReaderMap() == X);
        EXPECT_TRUE(D_mapped.getReaderMap() == D);
    }

    EXPECT_TRUE(fs::remove(X_file));
    EXPECT_TRUE(fs::remove(D_file));
}

// ========================================================================================
} /* End of namespace trex::test::tsolvers::ad::linear_model */
