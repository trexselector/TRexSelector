// ===================================================================================
// tools_solver.hpp
// ===================================================================================
#ifndef TSOLVERS_LINEAR_MODEL_OMP_BASED_TOOLS_SOLVER_HPP
#define TSOLVERS_LINEAR_MODEL_OMP_BASED_TOOLS_SOLVER_HPP
// ===================================================================================
/**
 * @file tools_solver.hpp
 *
 * @brief Terminating Orthogonal Least Squares (T-OLS) solver.
 */
// ===================================================================================

// std includes
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// Cereal includes
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>

// Base class
#include <tsolvers/ad/linear_model/omp_based/tomp_solver.hpp>

// ===================================================================================

// Embedding into namespace trex::tsolvers::ad::linear_model::omp_based
namespace trex::tsolvers::ad::linear_model::omp_based {

// ===================================================================================

/**
 * @brief Terminating Orthogonal Least Squares (T-OLS) solver.
 *
 * @details Inherits the execution flow and state management from TSolverOMP, but
 * overrides the variable selection step to project candidate columns onto
 * the orthogonal complement of the active set prior to correlation calculation.
 *
 * The selection score of candidate j is |x_j^T r| / ||P_perp x_j||, with r the
 * least-squares residual on the active set (r is orthogonal to the active span,
 * so x_j^T r = (P_perp x_j)^T r). The squared projected norms of all candidates
 * are kept and downdated by (q^T x_j)^2 for every new orthonormal direction q of
 * the active span (q = X_A R^{-1} e_k, the k-th column of the thin Q factor),
 * inside the correlation pass: one O(n (p + L)) pass per step, as for T-OMP,
 * instead of re-orthogonalizing every candidate against the whole active set
 * (the modified Gram-Schmidt form of OLS, without an orthogonalized copy of
 * the dictionary).
 *
 *  References:
 *  -----------
 *   S. Chen, S. A. Billings, and W. Luo, "Orthogonal Least Squares Methods and Their Application to
 *     Non-Linear System Identification," Int. J. Control, vol. 50, no. 5, pp. 1873-1896, 1989,
 *     doi: 10.1080/00207178908953472.
 *   T. Blumensath and M. E. Davies, "On the Difference Between Orthogonal Matching Pursuit and Orthogonal
 *     Least Squares," unpublished note, draft of March 29, 2007.
 *   F. Scheidt and M. Muma, "Terminating Greedy Pursuits for High-Dimensional False Discovery Rate Control,"
 *     in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026, pp. 2511-2515.
 */
class TSolverOOLS : public tsolvers::ad::linear_model::omp_based::TSolverOMP {
public:
    // ==========================================================================
    // Constructor/Destructor
    // ==========================================================================

    /**
     * @brief Construct a new TSolverOOLS object.
     *
     * @param X Map to the predictor matrix (n x p).
     * @param D Map to the dummy predictor matrix (n x num_dummies).
     * @param y Map to the response vector (n x 1).
     * @param options Every other input (TSolverOptions): the preprocessing flags
     *        normalize (normalize the predictors), intercept (include an
     *        intercept term) and verbose (verbose output), the tolerance and
     *        the tie-breaking (exch_tie_alpha / exch_tie_floor, the
     *        exchangeable-tie band of the greedy solvers).
     */
    TSolverOOLS(Eigen::Map<Eigen::MatrixXd>& X,
                 Eigen::Map<Eigen::MatrixXd>& D,
                 Eigen::Map<Eigen::VectorXd>& y,
                 const TSolverOptions& options = TSolverOptions{});

    /** @brief Default constructor for serialization. */
    TSolverOOLS();

    /** @brief Virtual destructor. */
    ~TSolverOOLS() override = default;

    /** @brief Deleted copy operations. */
    TSolverOOLS(const TSolverOOLS&) = delete;

    /** @brief Deleted copy assignment operator. */
    TSolverOOLS& operator=(const TSolverOOLS&) = delete;

    /** @brief Defaulted move constructor. */
    TSolverOOLS(TSolverOOLS&&) = default;

    /** @brief Deleted move assignment operator. */
    TSolverOOLS& operator=(TSolverOOLS&&) = delete;

    /** @brief Get the string representation of the solver type. */
    std::string solverTypeToString() const override {
        return "TOOLS";
    }

    /** @brief TSolverOMP::memoryBytes() plus the projected norms of all columns. */
    std::size_t memoryBytes() const noexcept override;

    // ==========================================================================
    // Serialization
    // ==========================================================================

    /** @brief Friend class for serialization. */
    friend class cereal::access;

    /**
     * @brief Serialize the TSolverOOLS object.
     *
     * @tparam Archive The type of the archive.
     * @param archive The archive object.
     */
    template <class Archive>
    void serialize(Archive& archive) {
        archive(
            cereal::base_class<trex::tsolvers::ad::linear_model::omp_based::TSolverOMP>(this),
            CEREAL_NVP(proj_sq_), CEREAL_NVP(proj_rank_)
        );
    }

    /**
     * @brief Save the current TSolverOOLS object to a file (binary serialization).
     *
     * @param filename The name of the file to save the object to.
     */
    virtual void save(const std::string& filename) const override;

    /**
     * @brief Load a TSolverOOLS object from a file (binary deserialization).
     *
     * @param filename The name of the file to load the object from.
     * @param X Map to the predictor matrix (used for reconnecting after loading).
     * @param D Map to the dummy predictor matrix (used for reconnecting after loading).
     * @return A TSolverOOLS object loaded from the specified file.
     */
    static TSolverOOLS load(const std::string& filename,
                             Eigen::Map<Eigen::MatrixXd>& X,
                             Eigen::Map<Eigen::MatrixXd>& D);

protected:
    // ==========================================================================
    // OLS State
    // ==========================================================================

    /** @brief Squared norms ||P_perp x_j||^2 of all columns (unified index)
     *  projected onto the orthogonal complement of the active span; current
     *  for the inactive columns. */
    Eigen::VectorXd proj_sq_{};

    /** @brief Number of active positions already downdated into proj_sq_. */
    std::size_t proj_rank_{0};

    // ==========================================================================
    // Core OLS Selection
    // ==========================================================================

    /**
     * @brief Find inactive variables tied for the maximum OLS selection score
     * |x_j^T r| / ||P_perp x_j|| (candidates with a numerically vanishing
     * projection are excluded).
     */
    std::pair<double, std::vector<std::size_t>> findTiedMaxCorrelations() const override;

    /**
     * @brief Downdate proj_sq_ for the active positions added since the last
     * call, fused with the refresh of the correlations X^T r of the inactives.
     */
    void updateCorrelations() override;
};

/** @brief Old name of TSolverOOLS.
 *  @deprecated Use TSolverOOLS. */
using TOOLS_Solver [[deprecated("Use TSolverOOLS")]] = TSolverOOLS;

// ==================================================================================
} /* End of namespace trex::tsolvers::ad::linear_model::omp_based */
// ==================================================================================
#endif /* TSOLVERS_LINEAR_MODEL_OMP_BASED_TOOLS_SOLVER_HPP */
