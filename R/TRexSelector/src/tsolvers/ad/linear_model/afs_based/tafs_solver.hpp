// ===================================================================================
// tafs_solver.hpp
// ===================================================================================
#ifndef TSOLVERS_LINEAR_MODEL_AFS_BASED_TAFS_SOLVER_HPP
#define TSOLVERS_LINEAR_MODEL_AFS_BASED_TAFS_SOLVER_HPP
// ===================================================================================
/**
 * @file tafs_solver.hpp
 *
 * @brief Terminating Adaptive Forward Stepwise (T-AFS) solver.
 */
// ===================================================================================

// std includes
#include <string>

// Eigen includes
#include <Eigen/Dense>

// Cereal includes
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>

// Base class
#include <tsolvers/ad/tsolver_base.hpp>

// ===================================================================================

// Embedding into namespace trex::tsolvers::ad::linear_model::afs_based
namespace trex::tsolvers::ad::linear_model::afs_based {

// ===================================================================================

/**
 * @brief Terminating Adaptive Forward Stepwise (T-AFS) solver.
 * Implements and adapts the AFS algorithm by Zhang & Tibshirani (2026) into the
 * terminating solver framework. It interpolates between Forward Stepwise (rho = 1)
 * and LAR/LASSO (rho -> 0).
 *
 * @note The reported DoF is |active set| + intercept; with rho < 1 the
 * effective degrees of freedom of the shrunken fit are smaller, so
 * Cp values based on it penalize conservatively.
 *
 *  References:
 *  -----------
 *   I. Zhang and R. Tibshirani, "Adaptive Forward Stepwise: A Method for High Sparsity Regression," J. Mach.
 *     Learn. Res., vol. 27, no. 35, pp. 1-24, 2026.
 *   F. Scheidt and M. Muma, "Terminating Greedy Pursuits for High-Dimensional False Discovery Rate Control,"
 *     in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026, pp. 2511-2515.
 */
class TSolverAFS : public TSolverBase {
protected:
    /**
     * @brief Shrinkage step size in (0, 1]. Forward Stepwise (rho = 1) to LAR/LASSO
     * (rho -> 0). Set from TSolverOptions::rho_afs, default 0.3.
     */
    double rho_{0.3};

public:
    // ==========================================================================
    // Constructor/Destructor
    // ==========================================================================

    /**
     * @brief Construct a TSolverAFS object.
     *
     * @param X Map to the predictor matrix.
     * @param D Map to the dummy predictor matrix.
     * @param y Map to the response vector.
     * @param options Every other input (TSolverOptions); this solver reads
     *        rho_afs, the shrinkage parameter in (0, 1] (1 yields forward
     *        stepwise), and exch_tie_alpha / exch_tie_floor.
     */
    TSolverAFS(Eigen::Map<Eigen::MatrixXd>& X,
                Eigen::Map<Eigen::MatrixXd>& D,
                Eigen::Map<Eigen::VectorXd>& y,
                const TSolverOptions& options = TSolverOptions{});

    /** @brief Default constructor for serialization. */
    TSolverAFS();

    /** @brief Virtual destructor. */
    ~TSolverAFS() override = default;

    /** @brief Deleted copy constructor. */
    TSolverAFS(const TSolverAFS&) = delete;

    /** @brief Deleted copy assignment operator. */
    TSolverAFS& operator=(const TSolverAFS&) = delete;

    /** @brief Defaulted move operations. */
    TSolverAFS(TSolverAFS&&) = default;

    /** @brief Deleted move assignment operator. */
    TSolverAFS& operator=(TSolverAFS&&) = delete;

    // ==========================================================================
    // Main algorithm
    // ==========================================================================

    /**
     * @brief Execute T-AFS solution path until T_stop or early stopping criterion
     * is met.
     *
     * @param T_stop Maximum number of steps to execute. If 0, no limit is applied.
     * @param early_stop Whether to stop early if the stopping criterion is met.
     */
    void executeStep(std::size_t T_stop = 0, bool early_stop = true) override;

    /** @brief Return a string identifier for the solver type. */
    std::string solverTypeToString() const override {
        return "TAFS";
    }

    // ==========================================================================
    // Serialization
    // ==========================================================================

    /** @brief Friend class for cereal serialization. */
    friend class cereal::access;

    /**
     * @brief Serialize the TAFS solver object.
     *
     * @tparam Archive Cereal archive type for serialization.
     * @param archive Cereal archive object to serialize to.
     */
    template <class Archive>
    void serialize(Archive& archive) {
        archive(
            cereal::base_class<TSolverBase>(this),
            CEREAL_NVP(rho_)
        );
    }

    /**
     * @brief Save the current TAFS solver state to a file (binary serialization).
     *
     * @param filename Path to the output file.
     *
     * @note Numerical quantities are saved with a precision of 1e-12.
     */
    virtual void save(const std::string& filename) const;

    /**
     * @brief Load a TAFS solver state from a file (binary serialization).
     *
     * @param filename Path to the input file.
     * @param X Map to the predictor matrix.
     * @param D Map to the dummy predictor matrix.
     *
     * @return Deserialized TSolverAFS object connected to X and D.
     */
    static TSolverAFS load(const std::string& filename,
                            Eigen::Map<Eigen::MatrixXd>& X,
                            Eigen::Map<Eigen::MatrixXd>& D);

protected:
    // ==========================================================================
    // Internal Helpers
    // ==========================================================================

    /** @brief Updates the correlation vector against the current residual. */
    void updateCorrelations();
};

/** @brief Old name of TSolverAFS.
 *  @deprecated Use TSolverAFS. */
using TAFS_Solver [[deprecated("Use TSolverAFS")]] = TSolverAFS;

// ===================================================================================
} /* End of namespace trex::tsolvers::ad::linear_model::afs_based */
// ===================================================================================
#endif /* TSOLVERS_LINEAR_MODEL_AFS_BASED_TAFS_SOLVER_HPP */
