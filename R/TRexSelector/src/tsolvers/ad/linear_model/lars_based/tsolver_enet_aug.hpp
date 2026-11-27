// ============================================================================
// tsolver_enet_aug.hpp
// ============================================================================
#ifndef TSOLVERS_LINEAR_MODEL_LARS_BASED_TSOLVER_ENET_AUG_HPP
#define TSOLVERS_LINEAR_MODEL_LARS_BASED_TSOLVER_ENET_AUG_HPP
// ============================================================================
/**
 * @file tsolver_enet_aug.hpp
 *
 * @brief T-ElasticNet (Augmented) solver: thin wrapper around TSolverLASSO
 * running on a data-augmented system that encodes the EN ridge penalty.
 *
 * @details
 *  Builds the augmented system once in the constructor:
 *    d1 = sqrt(lambda2),  d2 = 1/sqrt(1+lambda2)
 *
 *    X_aug = d2 * [ X    ]   ((n+p+L) × p)
 *                 [ d1·Ip]
 *                 [ 0_L  ]
 *
 *    D_aug = d2 * [ D    ]   ((n+p+L) × L)
 *                 [ 0_p  ]
 *                 [ d1·IL]
 *
 *    y_aug = [ y ; 0_{p+L} ]
 *
 *  Special case: when lambda2 = 0 the augmented rows are all zero and the
 *  inner TLASSO would see inflated effective_n = n+p+L instead of n, giving
 *  wrong max_actives.  The constructor therefore skips augmentation and passes
 *  the original X, D, y directly to the inner TLASSO (TENETAug collapses to
 *  TLASSO, consistent with TENET collapsing to TLASSO at lambda2 = 0).
 *
 *  All preprocessing (centering, L2-normalisation) and the LARS path are
 *  delegated to the inner TSolverLASSO via the normalize / intercept flags.
 *
 *  Coefficient recovery:
 *    getBetaPathSparse() returns the inner path with values / d2
 *    = beta * sqrt(1+lambda2)
 */
// ============================================================================

// std includes
#include <memory>
#include <string>

// Cereal includes
#include <cereal/cereal.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>

// Eigen includes
#include <Eigen/Dense>

// tsolvers includes
#include <tsolvers/ad/tsolver_base.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_lasso.hpp>

// ============================================================================

namespace trex::tsolvers::ad::linear_model::lars_based {

// ============================================================================

/**
 * @brief T-ElasticNet solver via data augmentation.
 *
 * @details Delegates all LARS logic to an inner TSolverLASSO running on
 * the augmented (n+p+L)-row system.  The augmentation encodes the ridge
 * penalty in the data itself, yielding the same Gram structure as the
 * T-Rex+GVS selector of the 1.x versions (trex(method = "trex+GVS")).
 *
 * Public interface: pass X, D, y, lambda2 and call executeStep() /
 * getBetaPathSparse().  In-memory warm-starting is supported via the inner
 * solver. The augmented copies are owned, so memoryBytes() grows with
 * 8 (n + p + L)(p + L + 1) bytes for lambda2 > 0 (8 n (p + L + 1) at 0).
 *
 *  References:
 *  -----------
 *   H. Zou and T. Hastie, "Regularization and Variable Selection via the Elastic Net," J. R. Stat. Soc. Ser.
 *     B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.
 *   B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least Angle Regression," Ann. Statist., vol. 32,
 *     no. 2, pp. 407-499, 2004.
 *   J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Grouped Variable Selection in
 *     High-Dimensional Linear Models Using the T-Knock Filter," in Proc. 30th Eur. Signal Process. Conf.
 *     (EUSIPCO), 2022, pp. 892-896.
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class TSolverENETAug : public TSolverBase {

    // ========================================================================
    // Owned augmented data
    // ========================================================================
    double   lambda2_{0.0};   ///< Ridge parameter
    double   d2_{1.0};        ///< 1 / sqrt(1 + lambda2) — beta-recovery scale

    Eigen::MatrixXd X_aug_owned_;  ///< Augmented, normalised X  ((n+p+L)×p)
    Eigen::MatrixXd D_aug_owned_;  ///< Augmented, normalised D  ((n+p+L)×L)
    Eigen::VectorXd y_aug_owned_;  ///< Augmented, centred y     (n+p+L)

    // Maps pointing into the above buffers (heap-allocated so that
    // their addresses are stable after the constructor runs)
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_aug_map_;
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> D_aug_map_;
    std::unique_ptr<Eigen::Map<Eigen::VectorXd>> y_aug_map_;

    // ========================================================================
    // Inner forward-selection solver
    // ========================================================================
    /** @brief Inner solver run on the augmented system. Held by base pointer so
     *  it can be either a TSolverLASSO (LARS-LASSO; can drop variables — the
     *  default, matching the augmented-Lasso EN recipe) or a plain TSolverLARS
     *  (pure LARS; never drops — matches R's trex(type="lar")). */
    std::unique_ptr<TSolverLARS> inner_;

    /** @brief If true, the inner solver is a pure-LARS TSolverLARS instead of
     *  the default TSolverLASSO. Diagnostic / R-parity option. */
    bool use_lars_inner_{false};

    // ========================================================================
    // Static helpers
    // ========================================================================

    /**
     * @brief Build the raw augmented X matrix (no normalisation).
     *
     * Constructs d2 · [X ; d1·Ip ; 0_L].  Centering and L2-normalisation
     * are delegated to the inner TSolverLASSO.
     *
     * @param X   Predictor matrix (n × p).
     * @param p   Number of original predictors.
     * @param L   Number of dummy predictors.
     * @param d1  sqrt(lambda2).
     * @param d2  1 / sqrt(1 + lambda2).
     * @return Raw augmented X, shape (n+p+L) × p.
     */
    static Eigen::MatrixXd buildXAug(
        const Eigen::Map<Eigen::MatrixXd>& X,
        std::size_t p,
        std::size_t L,
        double d1,
        double d2);

    /**
     * @brief Build the raw augmented D matrix (no normalisation).
     *
     * Constructs d2 · [D ; 0_p ; d1·IL].  Centering and L2-normalisation
     * are delegated to the inner TSolverLASSO.
     *
     * @param D   Dummy matrix (n × L).
     * @param p   Number of original predictors (X-ridge block rows).
     * @param L   Number of dummy predictors.
     * @param d1  sqrt(lambda2).
     * @param d2  1 / sqrt(1 + lambda2).
     * @return Raw augmented D, shape (n+p+L) × L.
     */
    static Eigen::MatrixXd buildDAug(
        const Eigen::Map<Eigen::MatrixXd>& D,
        std::size_t p,
        std::size_t L,
        double d1,
        double d2);

    /**
     * @brief Build the augmented response vector: y_aug = [y ; 0_{p+L}].
     *
     * Centering is delegated to the inner TSolverLASSO (intercept flag).
     *
     * @param y   Response vector (n).
     * @param pL  Number of augmentation rows (p + L).
     * @return Augmented y, shape (n + pL).
     */
    static Eigen::VectorXd buildYAug(
        const Eigen::Map<Eigen::VectorXd>& y,
        std::size_t pL);

    /** @brief Construct Eigen::Maps from the three owned buffers. */
    void buildMaps();

    /**
     * @brief Copy the inner solver's diagnostic state into the inherited
     * TSolverBase members so all base getters (getRSS, getActives,
     * getNumSteps, ...) report the inner solver's progress.
     *
     * Called after construction, after every executeStep(), and after
     * deserialization.
     */
    void syncDiagnosticsFromInner();

    /**
     * @brief Create the inner TSolverLASSO, forwarding preprocessing flags.
     * @param normalize  Passed to TSolverLASSO (L2-normalise columns).
     * @param intercept  Passed to TSolverLASSO (centre columns and y).
     * @param precentered Passed to TSolverLASSO (caller-centered data;
     *        false on the augmented system).
     * @param verbose    Passed to TSolverLASSO.
     * @param scaling_mode Passed to TSolverLASSO (L2 or z-score scaling).
     * @param tie_seed   Seed of the inner tie-break engine (TSolverOptions::tie_seed).
     * @param n_threads  Thread budget of the inner solver (TSolverOptions::n_threads).
     */
    void buildInnerSolver(bool normalize, bool intercept, bool precentered,
                          bool verbose, TSolverScalingMode scaling_mode, long long tie_seed,
                          int n_threads);

public:
    // ========================================================================
    // Constructors / destructor
    // ========================================================================

    /**
     * @brief Construct a TSolverENETAug.
     *
     * Builds the augmented system (d1/d2 formula), then creates an inner
     * TSolverLASSO that handles all centering and L2-normalisation.
     *
     * @param X          Predictor matrix (n × p).
     * @param D          Dummy matrix (n × L).
     * @param y          Response vector (n).
     * @param options    Every other input (TSolverOptions); this solver reads
     *                   lambda2, the ridge (L2) regularisation parameter
     *                   (≥ 0), the preprocessing flags forwarded to the inner
     *                   TSolverLASSO — normalize (columns of the augmented
     *                   matrices L2-normalised to unit norm, default true),
     *                   intercept (columns and response centred, default
     *                   true), verbose (default false) and scaling_mode (L2
     *                   or z-score column scaling, default L2) — and
     *                   tenet_aug_use_lars: if true, the inner
     *                   forward-selection solver is a pure-LARS TSolverLARS
     *                   (never drops variables) instead of the default
     *                   TSolverLASSO. Diagnostic / R-parity option.
     */
    TSolverENETAug(Eigen::Map<Eigen::MatrixXd>& X,
                    Eigen::Map<Eigen::MatrixXd>& D,
                    Eigen::Map<Eigen::VectorXd>& y,
                    const TSolverOptions& options = TSolverOptions{});

    /**
     * @brief Default constructor for Cereal deserialization.
     *
     * @note After deserialization the Maps and inner solver are rebuilt
     * automatically by the serialize() template (via buildMaps() +
     * buildInnerSolver()).
     */
    TSolverENETAug() = default;

    /** @brief Deleted copy constructor. */
    TSolverENETAug(const TSolverENETAug&) = delete;

    /** @brief Deleted copy assignment. */
    TSolverENETAug& operator=(const TSolverENETAug&) = delete;

    /** @brief Default move constructor. */
    TSolverENETAug(TSolverENETAug&&) = default;

    /** @brief Deleted move assignment. */
    TSolverENETAug& operator=(TSolverENETAug&&) = delete;

    /** @brief Virtual destructor. */
    ~TSolverENETAug() override = default;

    // ========================================================================
    // Core solver interface (TSolverBase pure virtuals)
    // ========================================================================

    /**
     * @brief Execute the T-LASSO solution path on the augmented system.
     *
     * Delegates to the inner TSolverLASSO; warm-starting is supported
     * (calling again with a larger T_stop continues the path).
     *
     * @param T_stop    Maximum number of dummy entries allowed in the active set.
     * @param early_stop If true, stop as soon as T_stop dummies have entered.
     */
    void executeStep(std::size_t T_stop = 0, bool early_stop = true) override;

    /** @brief Return the solver type as a string. */
    std::string solverTypeToString() const override { return "TENETAug"; }

    // ========================================================================
    // Solution extraction overrides
    // ========================================================================

    /**
     * @brief Return the sparse coefficient path scaled back to original units.
     *
     * The inner LASSO runs in the augmented space where each column is scaled
     * by d₂ = 1/sqrt(1+λ).  Dividing by d₂ (i.e. multiplying by sqrt(1+λ))
     * recovers coefficients in the same scale as the original predictors.
     */
    SparseBetaPath getBetaPathSparse() const override;

    /**
     * @brief Return coefficients at a given step in original units.
     * @param step Step index (-1 = last step).
     * @return Vector (p + L), original-scale coefficients.
     */
    Eigen::VectorXd getBeta(int step = -1) const override;

    /**
     * @brief Return the intercept at a given step (forwards to the inner solver).
     * @param step Step index (-1 = last step).
     */
    double getIntercept(int step = -1) const override;

    /** @brief Set the tolerance on the wrapper and on the inner solver. */
    void setTolerance(double tol) override { eps_ = tol; inner_->setTolerance(tol); }

    /** @brief Set the thread budget on the wrapper and on the inner solver
     *  (the inner solver owns every region). */
    void setNumThreads(int n_threads) override {
        TSolverBase::setNumThreads(n_threads);
        inner_->setNumThreads(n_threads);
    }

    /** @brief Set the tie-breaking seed on the inner solver (the wrapper's own
     *  random engine is unused). */
    void setTieSeed(uint32_t seed) override { inner_->setTieSeed(seed); }

    /** @brief The owned augmented X, D and y, the wrapper's diagnostic mirror,
     *  the inner solver object and view objects (heap allocations of the
     *  wrapper) and inner_->memoryBytes(); the inner X and D view the owned
     *  copies, so nothing is counted twice. */
    std::size_t memoryBytes() const noexcept override;

    // ========================================================================
    // Cereal serialization
    // ========================================================================

    friend class cereal::access;

    /**
     * @brief Serialize owned buffers, wrapper metadata, and inner solver state.
     *
     * @note X_ and D_ (non-owned TSolverBase pointers) are NOT serialized.
     * On load, the Maps are rebuilt from the owned buffers, the inner solver
     * is reconnected to them, and the wrapper's diagnostic mirror is refreshed.
     */
    template<class Archive>
    void serialize(Archive& archive) {
        archive(CEREAL_NVP(lambda2_),
                CEREAL_NVP(d2_),
                CEREAL_NVP(eps_),
                CEREAL_NVP(use_lars_inner_),
                CEREAL_NVP(p_original_), CEREAL_NVP(num_dummies_),
                CEREAL_NVP(dummy_start_idx_), CEREAL_NVP(effective_n_),
                CEREAL_NVP(normalize_), CEREAL_NVP(intercept_), CEREAL_NVP(precentered_),
                CEREAL_NVP(verbose_),
                CEREAL_NVP(X_aug_owned_),
                CEREAL_NVP(D_aug_owned_),
                CEREAL_NVP(y_aug_owned_));
        archive(CEREAL_NVP(inner_));

        // Rebuild maps and reconnect the inner solver ONLY on load (rebuilding
        // during save would destroy the Map objects the inner solver points
        // into). The rebuild must be unconditional on load: when deserializing
        // into an already-constructed solver (e.g. the R bindings' load path)
        // stale maps exist but point into buffers cereal may have reallocated,
        // and the freshly deserialized inner solver is always disconnected.
        if constexpr (Archive::is_loading::value) {
            if (X_aug_owned_.size() > 0) {
                buildMaps();
                if (inner_) {
                    inner_->reconnect(*X_aug_map_, *D_aug_map_);
                    syncDiagnosticsFromInner();
                }
            }
        }
    }

    /**
     * @brief Save the current TENETAug solver state to file (binary serialization).
     * @param filename Output file path.
     */
    void save(const std::string& filename) const;

    /**
     * @brief Load a TENETAug solver state from file (binary serialization).
     *
     * @details The wrapper owns its (augmented) data buffers, so unlike the
     * other solvers no reconnection to the passed maps happens; X and D are
     * only validated against the serialized dimensions for API parity.
     *
     * @param filename Path to input file.
     * @param X Map to the original feature matrix (dimension check only).
     * @param D Map to the original dummy matrix (dimension check only).
     * @return Deserialized TSolverENETAug.
     */
    static TSolverENETAug load(const std::string& filename,
                                Eigen::Map<Eigen::MatrixXd>& X,
                                Eigen::Map<Eigen::MatrixXd>& D);
};

/** @brief Old name of TSolverENETAug.
 *  @deprecated Use TSolverENETAug. */
using TENETAug_Solver [[deprecated("Use TSolverENETAug")]] = TSolverENETAug;

// ============================================================================
} // namespace trex::tsolvers::ad::linear_model::lars_based

// ============================================================================
// Cereal polymorphic registration
// ============================================================================
// The type is registered in tsolver_enet_aug.cpp, so its serialization code is
// compiled once. The macro below links that source into every program that
// includes this header.
#include <cereal/types/polymorphic.hpp>
CEREAL_FORCE_DYNAMIC_INIT(trex_tenet_aug_solver)

#endif /* TSOLVERS_LINEAR_MODEL_LARS_BASED_TSOLVER_ENET_AUG_HPP */
