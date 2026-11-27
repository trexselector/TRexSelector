// ===================================================================================
// tsolver_ncgmp.hpp
// ===================================================================================
#ifndef TSOLVERS_LINEAR_MODEL_OMP_BASED_TSOLVER_NCGMP_HPP
#define TSOLVERS_LINEAR_MODEL_OMP_BASED_TSOLVER_NCGMP_HPP
// ===================================================================================
/**
 * @file tsolver_ncgmp.hpp
 *
 * @brief Terminating Norm-Corrective Generalized Matching Pursuit (T-NCGMP) solver.
 *
 */
// ===================================================================================

// std includes
#include <string>

// Eigen includes
#include <Eigen/Dense>

// Cereal includes
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>

// tsolvers includes
#include <tsolvers/ad/tsolver_base.hpp>

// ===================================================================================

// Embedding into namespace trex::tsolvers::ad::linear_model::omp_based
namespace trex::tsolvers::ad::linear_model::omp_based {

// ===================================================================================

/** @brief Variants for the Norm-Corrective GMP solver — defined with the
 *  other option enums in tsolvers/tsolver_options.hpp and re-exported here. */
using trex::tsolvers::ad::NCGMPVariant;

/**
 * @brief Terminating Norm-Corrective Generalized Matching Pursuit (T-NCGMP) solver.
 *
 * @details Implements Algorithm 4 of Locatello et al. (2017).
 *
 * @note EQUIVALENCE FOR LINEAR MODELS: Algorithm 4 updates through the
 * quadratic upper bound b = x - (1/L) * grad f(x). For the least-squares
 * objective used here (L = 1, grad f = -r), Variant 0 (LineSearch) reduces
 * EXACTLY to Matching Pursuit (T-MP) and Variant 1 (FullyCorrective) reduces
 * EXACTLY to Orthogonal Matching Pursuit (T-OMP), stopping rules included
 * (at most effective_n_ active atoms, 8 x min(max_actives, effective_n_)
 * steps); both equivalences are pinned over the full path in the test
 * suite. The class is retained as
 * the GMP framing in anticipation of non-quadratic losses (e.g., the planned
 * logistic model), where the norm-corrective updates genuinely differ from
 * MP/OMP.
 *
 *  References:
 *  -----------
 *   F. Locatello, R. Khanna, M. Tschannen, and M. Jaggi, "A Unified Optimization View on Generalized Matching
 *     Pursuit and Frank-Wolfe," in Proc. 20th Int. Conf. Artif. Intell. Statist. (AISTATS), Proc. Mach.
 *     Learn. Res. vol. 54, 2017, pp. 860-868.
 *   F. Scheidt and M. Muma, "Terminating Greedy Pursuits for High-Dimensional False Discovery Rate Control,"
 *     in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026, pp. 2511-2515.
 */
class TSolverNCGMP : public TSolverBase {
protected:

    /** @brief Variant of the Norm-Corrective GMP solver. */
    NCGMPVariant variant_{NCGMPVariant::LINE_SEARCH};

public:
    // ==========================================================================
    // Constructor/Destructor
    // ==========================================================================

    /**
     * @brief Construct a new tncgmp solver object
     *
     * @param X Map to the predictor matrix (n x p).
     * @param D Map to the dummy predictor matrix (n x num_dummies).
     * @param y Map to the response vector (n x 1).
     * @param options Every other input (TSolverOptions); this solver reads
     *        ncgmp_variant (Line Search or Fully Corrective) and
     *        exch_tie_alpha / exch_tie_floor.
     */
    TSolverNCGMP(Eigen::Map<Eigen::MatrixXd>& X,
                  Eigen::Map<Eigen::MatrixXd>& D,
                  Eigen::Map<Eigen::VectorXd>& y,
                  const TSolverOptions& options = TSolverOptions{});

    /** @brief Default constructor for serialization. */
    TSolverNCGMP();

    /** @brief Virtual destructor. */
    ~TSolverNCGMP() override = default;

    /** @brief Deleted copy constructor. */
    TSolverNCGMP(const TSolverNCGMP&) = delete;

    /** @brief Deleted copy assignment operator. */
    TSolverNCGMP& operator=(const TSolverNCGMP&) = delete;

    /** @brief Defaulted move constructor. */
    TSolverNCGMP(TSolverNCGMP&&) = default;

    /** @brief Deleted move assignment operator. */
    TSolverNCGMP& operator=(TSolverNCGMP&&) = delete;

    // ==========================================================================
    // Main algorithm
    // ==========================================================================

    /**
     * @brief Execute a single step of the T-NCGMP algorithm.
     *
     * @param T_stop The maximum number of iterations to perform.
     * @param early_stop Whether to enable early stopping based on a convergence criterion.
     */
    void executeStep(std::size_t T_stop = 0, bool early_stop = true) override;

    /** @brief Get the string representation of the solver type. */
    std::string solverTypeToString() const override {
        return (variant_ == NCGMPVariant::LINE_SEARCH) ? "TNCGMP_LineSearch"
                                                      : "TNCGMP_FullyCorrective";
    }

    // ==========================================================================
    // Serialization
    // ==========================================================================

    /** @brief Friend class for serialization. */
    friend class cereal::access;

    /**
     * @brief Serialize the TSolverNCGMP object.
     *
     * @tparam Archive The type of the archive.
     * @param archive The archive object.
     */
    template <class Archive>
    void serialize(Archive& archive) {
        // Enum passes through a local int (reinterpret_cast aliasing is UB);
        // the assignment after archive() applies the loaded value and is a
        // no-op round-trip on save.
        int variant = static_cast<int>(variant_);
        archive(
            cereal::base_class<TSolverBase>(this),
            cereal::make_nvp("variant", variant)
        );
        variant_ = static_cast<NCGMPVariant>(variant);
    }

    /**
     * @brief Save the current TSolverNCGMP object to a file (binary serialization).
     *
     * @param filename The name of the file to save the object to.
     */
    virtual void save(const std::string& filename) const;

    /**
     * @brief Load a TSolverNCGMP object from a file (binary deserialization).
     *
     * @param filename The name of the file to load the object from.
     * @param X Map to the predictor matrix (required for reconstruction).
     * @param D Map to the dummy predictor matrix (required for reconstruction).
     * @return TSolverNCGMP The deserialized TSolverNCGMP object.
     */
    static TSolverNCGMP load(const std::string& filename,
                              Eigen::Map<Eigen::MatrixXd>& X,
                              Eigen::Map<Eigen::MatrixXd>& D);

protected:
    // ==========================================================================
    // Internal Helpers
    // ==========================================================================
    /** @brief Update the correlations between the residual and the predictors.
     *  Line Search (MP) refreshes all non-dropped variables to allow atom
     *  re-selection; Fully Corrective (OMP) refreshes inactive variables only. */
    void updateCorrelations();

    /** @brief Update the beta path and residuals. */
    void updateBetaPathAndResiduals();

    /** @brief Find tied max correlations. Line Search (MP) scans all non-dropped
     *  variables (re-selection); Fully Corrective (OMP) uses the base
     *  inactive-only scan. */
    std::pair<double, std::vector<std::size_t>>
        findTiedMaxCorrelations() const override;
};

/** @brief Old name of TSolverNCGMP.
 *  @deprecated Use TSolverNCGMP. */
using TNCGMP_Solver [[deprecated("Use TSolverNCGMP")]] = TSolverNCGMP;

// ===================================================================================
} /* End of namespace trex::tsolvers::ad::linear_model::omp_based */
// ===================================================================================
#endif /* TSOLVERS_LINEAR_MODEL_OMP_BASED_TSOLVER_NCGMP_HPP */
