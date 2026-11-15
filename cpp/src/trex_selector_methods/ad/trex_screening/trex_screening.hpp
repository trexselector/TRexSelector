// ===================================================================================
// trex_screening.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_SCREENING_HPP
#define TREX_SELECTOR_METHODS_TREX_SCREENING_HPP
// ===================================================================================
/**
 * @file trex_screening.hpp
 *
 * @brief T-Rex screening selector (T-Rex+Screening): fast FDR-controlled
 *        screening at T = 1, L = p.
 */
// ===================================================================================

// std includes
#include <cstddef>
#include <utility>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// TRex base
#include <trex_selector_methods/ad/trex_core/trex.hpp>

// ===================================================================================

namespace trex::trex_selector_methods::ad::trex_screening {

// Namespace aliases
namespace tc = trex::trex_selector_methods::ad::trex_core;
namespace er = trex::trex_selector_methods::ad::utils::experiment_runner;

// ===================================================================================
// Enums & Control Structures for T-Rex+Screening
// ===================================================================================

/**
 * @brief Variants of the screening: without or with a dependency-aware filter
 *        on the relative occurrences.
 */
enum class ScreenTRexMethod {
    /** @brief Screening without a dependency-aware filter. */
    TREX,

    /** @brief Screening with the dependency-aware filter for a first-order
     *  autoregressive dependency: a variable is compared with the neighbors
     *  inside the AR(1) window. */
    TREX_DA_AR1,

    /** @brief Screening with the dependency-aware filter for equi-correlated
     *  variables: a variable is compared with all other variables. */
    TREX_DA_EQUI,

    /** @brief Screening with the dependency-aware filter for block
     *  equi-correlated variables: the same comparison inside contiguous
     *  blocks. */
    TREX_DA_BLOCK_EQUI,

    /** @brief Screening with the dependency-aware filter of the nearest
     *  neighbors: a variable is compared with the neighbors inside the AR(1)
     *  window whose sample correlation reaches rho_thr_DA, and a variable
     *  without such a neighbor keeps its relative occurrence. */
    TREX_DA_NN
};


/** @brief Configuration for T-Rex+Screening specific parameters. */
struct ScreenTRexControlParameter {
    /** @brief Whether the confidence-based screening (bootstrap) makes the
     *  selection (default: false). The FDR estimate stays one over the number
     *  of selected variables. */
    bool use_bootstrap_CI = false;

    /** @brief Number of bootstrap resamples, at least 2 (default: 1000). */
    std::size_t R_boot = 1000;

    /** @brief Step of the grid over the confidence level, in (0, 1] (default: 0.001). */
    double ci_grid_step = 0.001;

    /** @brief Variant of the screening (default: TREX). See ScreenTRexMethod. */
    ScreenTRexMethod trex_method = ScreenTRexMethod::TREX;

    /** @brief Correlation threshold of the dependency-aware variants (default: 0.02).
     *  TREX_DA_AR1: the window around a variable reaches
     *  ceil(log(rho_thr_DA) / log(|cor_coef|)) columns to each side, up to
     *  the first lag whose correlation |cor_coef|^lag has fallen to the
     *  threshold, range [0, 1).
     *  TREX_DA_NN: the same window, of which the variables stay whose sample
     *  correlation with the variable reaches the threshold, range [0, 1).
     *  TREX_DA_EQUI and TREX_DA_BLOCK_EQUI: the filter is skipped when
     *  |cor_coef| does not exceed the threshold, range [0, 1].
     */
    double rho_thr_DA = 0.02;

    /** @brief Manually specified correlation coefficient.
     *  AUTO_ESTIMATE_CORRELATION (−2.0) triggers automatic estimation. The
     *  dependency-aware variants use the magnitude of a supplied coefficient.
     *  A supplied value outside [−1, 1] gives a warning and is replaced by
     *  the estimate.
     */
    double cor_coef = tc::AUTO_ESTIMATE_CORRELATION;

    /** @brief Number of blocks for block-equicorrelated DA variant (default: 5).
     *  TREX_DA_BLOCK_EQUI takes values between 1 and the number of variables.
     */
    std::size_t n_blocks = 5;

    /** @brief Base T-Rex algorithmic control parameters (nested).
     *  T-Rex+Screening accepts `lloop_strategy` in {STANDARD, SEEDED,
     *  PERMUTATION, PERMUTATION_SEEDED} (validated by
     *  `validateScreenTRexStrategy()`; HCONCAT / SKIPL describe L-loop
     *  growth, which screening does not have). The base class's own
     *  default (STANDARD) already satisfies this, so no override is
     *  needed. `max_dummy_multiplier` is pinned to
     *  1 internally (screening always runs at exactly L = p dummies),
     *  which also sizes memory-mapped D files at the width actually used.
     */
    tc::TRexControlParameter trex_ctrl;
};


/** @brief Result structure for T-Rex+Screening, extending SelectionResult.
 *
 * @details
 *  Adds T-Rex+Screening specific diagnostics on top of the standard result.
 */
struct ScreenTRexSelectionResult : public tc::TRexSelector::SelectionResult {
    /** @brief Estimated FDR = 1 / max(1, R) with R selected variables, so
     *  an empty selection reports 1. */
    double estimated_FDR = 0.0;

    /** @brief Whether the confidence-based screening (bootstrap) made the
     *  selection. */
    bool used_bootstrap = false;

    /** @brief Confidence level that the grid search of the bootstrap picked.
     *  0 without the bootstrap. */
    double confidence_level = 0.0;

    /** @brief Coefficients of the variables, averaged over the K experiments
     *  (one per variable, p x 1). */
    Eigen::VectorXd beta_mat;

    /** @brief Non-zero coefficients of the dummies, collected over the K
     *  experiments. They are the sample of the bootstrap. */
    Eigen::VectorXd dummy_betas;

    /** @brief Correlation coefficient of the dependency-aware variants: the
     *  supplied value, or the estimate from the data. 0 for TREX. */
    double estimated_correlation = 0.0;
};


// ===================================================================================
// ScreenTRexSelector
// ===================================================================================

/**
 * @brief T-Rex screening selector (T-Rex+Screening): fast FDR-controlled variable
 *        screening in one pass of the random experiments, until one dummy is
 *        included (T = 1), with as many dummies as variables (L = p).
 *
 * @details
 *  Overrides TRexSelector::select() to run K experiments at T=1, L=p with no
 *  T-loop or L-loop calibration.  Supports ordinary thresholding (Phi > 0.5),
 *  a bootstrap-CI variant, and dependency-aware filters on Phi for AR(1) and
 *  equi-correlated designs.
 *
 *  Getters of the base class after a screening run:
 *  ------------------------------------------------
 *  The screening has no target FDR and no calibration, so the getters of the
 *  base class hold fixed values: getTStop() is 1, getDummyMultiplierL() is 1,
 *  getNumDummies() is p, getVotingThreshold() is 0.5 and getTFDR() is 0.
 *  getRMat() is a 1 x 1 matrix with the number of selected variables,
 *  getPhiPrime() equals the relative occurrences, and getFDPHatMat() is a
 *  1 x 1 matrix with 0, because no FDP is estimated per voting level. The FDR
 *  estimate of the run is estimated_FDR of the result.
 *
 *  Supported dummy strategies: STANDARD, SEEDED, PERMUTATION, and
 *  PERMUTATION_SEEDED (single pass at L = p, T = 1 — the L-loop growth
 *  strategies HCONCAT / SKIPL have no meaning here). SEEDED is the
 *  memory-lean twin of STANDARD: the identical seeded draws, regenerated
 *  one D_k at a time in the serial runner loop and never stored.
 *  PERMUTATION_SEEDED is the stateless twin of PERMUTATION: bit-identical
 *  experiments (both key off the base class's resolved
 *  permutation_base_seed_), but the base matrix is re-derived from the
 *  seed instead of being stored for the selector's lifetime.
 *  Memory-mapped D storage is supported through the base class MemmapManager
 *  when TRexControlParameter::use_memory_mapping is set.
 *
 *  Variants:
 *  1. Ordinary T-Rex+Screening: selects variables where Phi > 0.5.
 *  2. Confidence-Based T-Rex+Screening: bootstrap CI on dummy betas.
 *  3. DA-AR1: adjusts Phi using a local AR(1) neighborhood window.
 *  4. DA-EQUI: adjusts Phi using a global equi-correlation comparison.
 *  5. DA-BLOCK-EQUI: the same comparison inside contiguous blocks.
 *  6. DA-NN: adjusts Phi using the neighbors inside the AR(1) window whose
 *     sample correlation reaches the threshold. A variable without such a
 *     neighbor keeps its Phi.
 *
 *  References:
 *  -----------
 *   J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Fast Screening of Large-Scale
 *     Genomics Biobanks," in Proc. 22nd IEEE Statist. Signal Process. Workshop (SSP), 2023, pp. 666-670.
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class ScreenTRexSelector : public tc::TRexSelector {

protected:
    // ==========================================================================
    // T-Rex+Screening Data Members
    // ==========================================================================

    /** @brief T-Rex+Screening control parameters (nests the base trex_ctrl). */
    ScreenTRexControlParameter trex_screen_ctrl_;

    /** @brief Sum of betas for original features across K experiments (p × 1). */
    Eigen::MatrixXd accumulated_beta_mat_;

    /** @brief Non-zero dummy betas collected across K experiments. */
    std::vector<double> collected_dummy_betas_;

    /** @brief Full T-Rex+Screening result from the last select() call. */
    ScreenTRexSelectionResult screen_result_;

public:
    // ==========================================================================
    // Constructor / Destructor
    // ==========================================================================

    /**
     * @brief Construct a ScreenTRexSelector.
     *
     * @param X             Feature matrix (n × p) — Eigen::Map, not copied.
     * @param y             Response vector (n × 1) — Eigen::Map, copied internally.
     * @param trex_screen_ctrl T-Rex+Screening specific parameters, nesting the base
     *                      algorithmic control (solver, dummy strategy, memmap)
     *                      as `trex_screen_ctrl.trex_ctrl`. Accepted
     *                      lloop_strategy: STANDARD, SEEDED, PERMUTATION,
     *                      PERMUTATION_SEEDED.
     * @param seed          Random seed. A negative value (default: -1) gives a run
     *                      that is not reproducible.
     * @param verbose       Enable verbose output.
     *
     * @throws std::invalid_argument on an unsupported lloop_strategy
     *         (HCONCAT / SKIPL).
     */
    ScreenTRexSelector(
        Eigen::Map<Eigen::MatrixXd>& X,
        Eigen::Map<Eigen::VectorXd>& y,
        ScreenTRexControlParameter   trex_screen_ctrl = ScreenTRexControlParameter(),
        int  seed    = -1,
        bool verbose = true
    );

    /** @brief Virtual destructor. */
    virtual ~ScreenTRexSelector() = default;

    /** @brief Deleted copy constructor. */
    ScreenTRexSelector(const ScreenTRexSelector&) = delete;

    /** @brief Deleted copy assignment. */
    ScreenTRexSelector& operator=(const ScreenTRexSelector&) = delete;

    /** @brief Deleted move constructor. */
    ScreenTRexSelector(ScreenTRexSelector&&) = delete;

    /** @brief Deleted move assignment. */
    ScreenTRexSelector& operator=(ScreenTRexSelector&&) = delete;

    // ==========================================================================
    // Main method
    // ==========================================================================

    /**
     * @brief Run T-Rex+Screening selection (T=1, L=p).
     *
     * @details Overrides the base class select().  Bypasses the T-loop and L-loop
     *          to run exactly one fixed-budget experiment set, then selects variables
     *          via ordinary thresholding or bootstrap CI on dummy betas.
     *
     * @return SelectionResult (sliced to base type). Use getScreenResult() for
     *         full T-Rex+Screening diagnostics.
     */
    virtual SelectionResult select() override;

    // ==========================================================================
    // Public Getter
    // ==========================================================================

    /**
     * @brief Retrieve the full T-Rex+Screening result from the last select() call.
     *
     * @return Const reference to ScreenTRexSelectionResult.
     */
    const ScreenTRexSelectionResult& getScreenResult() const;


    /**
     * @brief Apply bootstrap screening to an existing ordinary T-Rex+Screening result.
     *
     * @details Reuses Phi, beta_mat, and dummy_betas from a prior ordinary
     *          select() call to derive the confidence-based T-Rex+Screening result
     *          without re-running the K experiments.  This is used by the biobank
     *          orchestrator (Algorithm 1) to obtain both α̂ and α̂_C from a single
     *          set of K experiments.
     *
     * @param ordinary_result  Result from a prior ordinary select() call.
     *                         Must contain non-empty dummy_betas and beta_mat.
     *
     * @return ScreenTRexSelectionResult with bootstrap-based selection.
     */
    ScreenTRexSelectionResult applyBootstrapToOrdinaryResult(
        const ScreenTRexSelectionResult& ordinary_result
    );


protected:
    // =========================================================================
    // Core T-Rex+Screening Operations
    // =========================================================================

    /**
     * @brief Run K experiments at T=1, L=p and collect beta statistics.
     *
     * @details  Uses ExperimentRunner with collect_screen_betas=true.
     *           Populates accumulated_beta_mat_ and collected_dummy_betas_.
     *
     * @return ExperimentResults containing phi_T_mat and Phi.
     */
    er::ExperimentResults runScreenExperiments();


    /**
     * @brief Ordinary T-Rex+Screening: select variables where Phi > 0.5.
     *
     * @param Phi        Relative occurrences (p × 1).
     * @param result_out Result object to populate.
     */
    void performOrdinarySelection(
        const Eigen::VectorXd&    Phi,
        ScreenTRexSelectionResult& result_out
    ) const;


    /**
     * @brief Confidence-Based T-Rex+Screening: bootstrap CI on dummy betas.
     *
     * @details Grid-searches gamma ∈ [0, 1) for the smallest confidence level
     *          such that the number of selections ≤ ordinary T-Rex+Screening count.
     *
     * @param Phi        Relative occurrences (p × 1).
     * @param beta_means Averaged beta coefficients (p × 1).
     * @param result     Result object to populate with bootstrap diagnostics.
     */
    void performBootstrapSelection(
        const Eigen::VectorXd&    Phi,
        const Eigen::VectorXd&    beta_means,
        ScreenTRexSelectionResult& result
    );


    /**
     * @brief Compute estimated FDR = 1 / max(1, R).
     *
     * @param result Reference to result object (read selected_var, write estimated_FDR).
     */
    void computeEstimatedFDR(ScreenTRexSelectionResult& result) const;


    /**
     * @brief Calculate a bias-corrected normal bootstrap confidence interval.
     *
     * @details Matches R's boot::boot.ci(type = "norm"):
     *          center = 2·t₀ − mean(t*),  CI = center ± z·se(t*).
     *
     * @param boot_means       Bootstrap distribution of means (t*).
     * @param t0               Original statistic computed on the full sample.
     * @param confidence_level Confidence level in (0, 1).
     *
     * @return Pair {lower_bound, upper_bound}.
     */
    std::pair<double, double> calculateBootstrapCI(
        const std::vector<double>& boot_means,
        double t0,
        double confidence_level
    ) const;


    /**
     * @brief Inverse normal CDF (probit).
     *
     * @details Uses Boost Math erf_inv: z = √2 · erfinv(2p − 1).
     *
     * @param p Probability in (0, 1).
     * @return Corresponding z-score.
     */
    double probit(double p) const;


    // =========================================================================
    // Dependency-Aware (DA) Helpers
    // =========================================================================

    /**
     * @brief Apply Dependency-Aware Phi adjustment (AR1 or EQUI variant).
     *
     * @param Phi Relative occurrences (modified in-place).
     */
    void applyDependencyAwareAdjustment(Eigen::VectorXd& Phi);

    /**
     * @brief Estimate AR(1) correlation coefficient from X.
     *
     * @details Row-wise Yule-Walker estimator, averaged over n rows.
     *
     * @return Absolute average AR(1) autocorrelation.
     */
    double estimateAR1Correlation() const;

    /**
     * @brief Estimate equi-correlation coefficient from X.
     *
     * @details Uses variance of row sums: ρ̄ = (‖row_sums‖² − p) / (p(p−1)).
     *
     * @return Average off-diagonal correlation.
     */
    double estimateEquiCorrelation() const;

    /**
     * @brief Compute DA adjustment deltas for AR(1) strategy.
     *
     * @details Window size κ = ⌈log(ρ_thr) / log(ρ)⌉.  For each j, finds the
     *          minimum |Phi[j] − Phi[k]| in [j−κ, j+κ].
     *
     * @param Phi Relative occurrences (p × 1).
     * @param rho AR(1) correlation estimate.
     *
     * @return DA_delta adjustment vector (p × 1).
     */
    Eigen::VectorXd computeDA_AR1_Delta(
        const Eigen::VectorXd& Phi, double rho) const;

    /**
     * @brief Compute DA adjustment deltas for equi-correlation strategy.
     *
     * @details Exact O(p²) global minimum |Phi[j] − Phi[k]|.
     *          Skipped entirely when |ρ| ≤ ρ_thr_DA.
     *
     * @param Phi Relative occurrences (p × 1).
     * @param rho Equi-correlation estimate.
     *
     * @return DA_delta adjustment vector (p × 1).
     */
    Eigen::VectorXd computeDA_EQUI_Delta(
        const Eigen::VectorXd& Phi, double rho) const;

    /**
     * @brief Estimate within-block equi-correlation from X.
     *
     * @details Averages the equi-correlation estimate over each block.
     *
     * @return Average within-block correlation.
     */
    double estimateBlockEquiCorrelation() const;

    /**
     * @brief Compute DA adjustment deltas for block-equi-correlation strategy.
     *
     * @details Within-block neighbors only: for each j, finds
     *          min_{k ∈ block(j), k≠j} |Phi[j] − Phi[k]|.
     *
     * @param Phi Relative occurrences (p × 1).
     * @param rho Block equi-correlation estimate.
     *
     * @return DA_delta adjustment vector (p × 1).
     */
    Eigen::VectorXd computeDA_BLOCK_EQUI_Delta(
        const Eigen::VectorXd& Phi, double rho) const;

    /**
     * @brief Compute DA adjustment deltas for the nearest-neighbor strategy.
     *
     * @details The group of j holds the variables k in the AR(1) window
     *          [j−κ, j+κ], k≠j, whose sample correlation |corr(x_j, x_k)|
     *          reaches ρ_thr_DA. DA_delta[j] = 2 − min over the group of
     *          |Phi[j] − Phi[k]|. An empty group gives DA_delta[j] = 1, so
     *          the variable keeps its Phi. Cost O(n p κ).
     *
     * @param Phi Relative occurrences (p × 1).
     * @param rho AR(1) correlation estimate, which sets the window.
     *
     * @return DA_delta adjustment vector (p × 1).
     */
    Eigen::VectorXd computeDA_NN_Delta(
        const Eigen::VectorXd& Phi, double rho) const;


    // =========================================================================
    // Validation
    // =========================================================================

    /**
     * @brief Validate that lloop_strategy is one of STANDARD, SEEDED,
     *        PERMUTATION, PERMUTATION_SEEDED.
     *
     * @throws std::invalid_argument for unsupported strategies (HCONCAT /
     *         SKIPL, which describe L-loop growth).
     */
    void validateScreenTRexStrategy() const;

    /**
     * @brief Check the parameters of the dependency-aware variants.
     *
     * `rho_thr_DA` must be in [0, 1], and in [0, 1) for TREX_DA_AR1 and
     * TREX_DA_NN, whose window has no lag left at a threshold of 1.
     * `n_blocks` must be in [1, p] for TREX_DA_BLOCK_EQUI.
     *
     * @throws std::invalid_argument on a value outside its range.
     */
    void validateScreenDAParameters() const;

    /**
     * @brief Check the parameters of the bootstrap: `R_boot` at least 2 and
     *        `ci_grid_step` in (0, 1].
     *
     * The biobank screening applies the bootstrap to every phenotype, so the
     * check runs whether or not `use_bootstrap_CI` is set.
     *
     * @throws std::invalid_argument on a value outside its range.
     */
    void validateBootstrapParameters() const;

};

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_screening */
// ===================================================================================
#endif /* End of TREX_SELECTOR_METHODS_TREX_SCREENING_HPP */
