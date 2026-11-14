// ===================================================================================
// trex_da.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_DA_HPP
#define TREX_SELECTOR_METHODS_TREX_DA_HPP
// ===================================================================================
/**
 * @file trex_da.hpp
 *
 * @brief Dependency-aware T-Rex selector (T-Rex+DA).
 *
 * @details
 *  Extends the base class TRexSelector with dependency-aware deflation that accounts
 *  in special cases for the correlation structure among predictors.
 *
 *  Supported methods:
 *  -------------------
 *  - AR1:          AR(1) windowed deflation (ordered predictors).
 *  - EQUI:         Global equicorrelation deflation.
 *  - BT:           Binary-tree dendrogram from hierarchical clustering.
 *  - NN:           Nearest-neighbour threshold sweep on the correlation matrix.
 *  - PRIOR_GROUPS: User-supplied group structure at multiple hierarchy levels.
 *
 * Calibration Paths:
 *  - BT, NN, and PRIOR_GROUPS ("BT-style") perform 3D calibration over
 *    (T_stop, V, rho_grid);
 *  - AR1 and EQUI stay on the 2D (T_stop, V) grid of the base class.
 *
 * Penalty: all five methods deflate the whole occurrence history with the penalty
 * of the last step. This deviates from the papers, which write one penalty per
 * step (see the class TRexDASelector for the reason and the evidence).
 *
 * FDR limitation: BT, NN and PRIOR_GROUPS can exceed the target FDR at a very
 * high correlation within the groups (see the class TRexDASelector).
 */
// ===================================================================================

// std includes
#include <cstddef>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// TRex base class
#include <trex_selector_methods/ad/trex_core/trex.hpp>

// Hierarchical clustering
#include <ml_methods/clustering/hierarchical/agglomerative/agglomerative_dispatcher.hpp>
#include <ml_methods/clustering/hierarchical/agglomerative/agglomerative_types.hpp>
#include <ml_methods/clustering/hierarchical/agglomerative/dendrogram_utils.hpp>
#include <ml_methods/clustering/hierarchical/agglomerative/distance_policy.hpp>

// ===================================================================================

namespace trex::trex_selector_methods::ad::trex_da {

// Namespace aliases
namespace tc  = trex::trex_selector_methods::ad::trex_core;
namespace er  = trex::trex_selector_methods::ad::utils::experiment_runner;
namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;

// ===================================================================================
// Enums
// ===================================================================================

/**
 * @brief Dependency-aware correction methods.
 */
enum class DAMethod {
    /** @brief AR(1) windowed correction for ordered predictors. Made for a
     *  first-order autoregressive dependency across the columns. */
    AR1,

    /** @brief Equicorrelation correction (all-vs-all). Made for variables
     *  that all share one correlation. */
    EQUI,

    /** @brief Binary-tree dendrogram-based clustering. Made for disjoint
     *  groups of highly correlated variables. */
    BT,

    /** @brief Nearest-neighbour correlation threshold sweep. Made for
     *  overlapping groups of highly correlated variables. */
    NN,

    /** @brief User-supplied group labels at multiple hierarchy levels. The
     *  selector sets this method itself when prior_groups is supplied.
     *  Choosing it without prior_groups is an error of the constructor. */
    PRIOR_GROUPS
};


/**
 * @brief Candidate-cell policy for the variable selection of the methods BT,
 *        NN and PRIOR_GROUPS.
 *
 * @details
 *  In selectVariables_BT the winning selection count `val_max` is the largest
 *  count over FDP-feasible cells (FDP_hat <= tFDR). This flag controls which
 *  cells may then win the tie-break among those equal to `val_max`. AR1 and
 *  EQUI ignore it.
 */
enum class BTSelectionMode {
    /** @brief Only FDP-feasible cells (FDP_hat <= tFDR) may be selected. This
     *  controls the realized FDR and is the default. */
    FeasibleOnly,

    /** @brief Any cell whose count equals `val_max` may win, including
     *  FDP-infeasible ones. Reproduces the 1.x versions of the R package
     *  (`select_var_fun_DA_BT`, which scans the whole R_array), including
     *  their tendency to inflate the realized FDR at high within-group
     *  correlation. */
    RFaithful
};

// ===================================================================================
// Control Parameters
// ===================================================================================

/**
 * @brief DA-specific control parameters.
 */
struct TRexDAControlParameter {

    /** @brief Dependency-aware method (default: BT). Ignored when prior_groups is non-empty. */
    DAMethod method = DAMethod::BT;

    /** @brief Correlation coefficient for AR1/EQUI.
     *  AUTO_ESTIMATE_CORRELATION (-2.0) = auto-estimate from X. Both methods
     *  use the magnitude of a supplied coefficient. A supplied value outside
     *  [-1, 1] gives a warning and is replaced by the estimate.
     *
     *  @note The estimates: for EQUI, the mean of all pairwise Pearson
     *  correlations (as in the 1.x versions of the R package). For AR1, the
     *  magnitude of the lag-1 Yule-Walker autocorrelation of each row,
     *  averaged over the rows. The 1.x versions of the R package fit an
     *  AR(1) model per row via `arima(..., method="ML")` instead. The two
     *  estimators agree asymptotically and can differ slightly on short
     *  rows. Supply cor_coef to reproduce a result of those versions.
     */
    double cor_coef = tc::AUTO_ESTIMATE_CORRELATION;

    /** @brief Correlation threshold of the AR1 and EQUI methods (default: 0.02).
     *  AR1: the window around a variable reaches
     *  ceil(log(rho_thr_DA) / log(|cor_coef|)) columns to each side, up to
     *  the first lag whose correlation |cor_coef|^lag has fallen to the
     *  threshold, range [0, 1). EQUI: the correction is skipped when
     *  |cor_coef| does not exceed the threshold, range [0, 1].
     */
    double rho_thr_DA = 0.02;

    /** @brief Linkage of the hierarchical clustering (default: Single). Used
     *  by BT and by prior groups.
     *  Supported: Single, Complete, Average (UPGMA), WPGMA (McQuitty).
     *  Ward is rejected, because its update is valid for squared Euclidean
     *  distances only and the correlation distance 1 - |cor| is not Euclidean.
     */
    hac::LinkageMethod hc_linkage = hac::LinkageMethod::Single;

    /** @brief Number of correlation thresholds on the calibration grid of BT,
     *  NN and prior groups. 0 (default) takes min(20, p).
     */
    std::size_t hc_grid_length = 0;

    /** @brief Prior grouping constraints: one or more label vectors of length p
     *  (non-negative labels). When non-empty, the method field is ignored and
     *  the finest common refinement of all supplied levels becomes a hard
     *  merge CONSTRAINT on the dependency structure: hierarchical clustering
     *  (hc_linkage, correlation distance) is run WITHIN each constraint group,
     *  the rho grid is built from the pooled within-group dendrogram heights
     *  (subsampled to hc_grid_length, capped by the conservative rho = 1
     *  singleton anchor), and the BT-style 3D calibration selects the
     *  operating point. The nearest-partner deflation therefore always acts
     *  on tight data-driven neighbourhoods INSIDE the known groups — the
     *  groups themselves are never used directly as deflation neighbourhoods
     *  (that degenerates for group sizes > 2).
     */
    std::vector<std::vector<Eigen::Index>> prior_groups;

    /** @brief Candidate-cell policy for the selection of the methods BT, NN
     *  and PRIOR_GROUPS (default: FeasibleOnly). FeasibleOnly restricts the
     *  selection to FDP-feasible cells (FDR-controlling); RFaithful reproduces
     *  the 1.x versions of the R package by also admitting infeasible cells.
     *  Ignored by AR1 and EQUI. See BTSelectionMode.
     */
    BTSelectionMode bt_selection_mode = BTSelectionMode::FeasibleOnly;

    /** @brief Base T-Rex algorithmic control parameters (nested).
     *  DA imposes no restrictions on solver_type / lloop_strategy, so the
     *  base class's own defaults are used as-is.
     */
    tc::TRexControlParameter trex_ctrl;
};

// ===================================================================================
// DA Setup Result (internal)
// ===================================================================================

/**
 * @brief Internal structure holding the precomputed neighbourhood information.
 */
struct DASetupResult {
    /** @brief True for BT, NN, or PRIOR_GROUPS (3D calibration paths). */
    bool use_BT_style = false;

    /** @brief Group-mate indices: gr_j_list[j][rho_level] = vector of neighbour
     *  indices for variable j.
     */
    std::vector<std::vector<std::vector<Eigen::Index>>> gr_j_list;

    /** @brief Correlation-level grid (length rho_grid_len). Null-like (size 0) for AR1/EQUI. */
    Eigen::VectorXd rho_grid;

    /** @brief Length of the rho_grid (0 for AR1/EQUI). */
    std::size_t rho_grid_len = 0;

    /** @brief Reference index into rho_grid for the 75 % operating point. */
    std::size_t opt_point_BT = 0;

    /** @brief Estimated or user-supplied correlation coefficient. */
    double cor_coef = tc::AUTO_ESTIMATE_CORRELATION;

    /** @brief AR(1) window half-width (0 unless AR1). */
    std::size_t kap = 0;
};

// ===================================================================================
// DA Correction Result (internal)
// ===================================================================================

/**
 * @brief Output of the DA deflation step.
 */
struct DACorrectionResult {
    /** @brief phi_T_mat (p x T_stop). For AR1 and EQUI every column is deflated
     *  with the penalty of the last step. The undeflated copy for BT-style.
     */
    Eigen::MatrixXd phi_T_mat;

    /** @brief Deflated Phi (p x 1). */
    Eigen::VectorXd Phi;

    /** @brief BT-style only: deflated Phi per rho level (p x rho_grid_len). */
    Eigen::MatrixXd Phi_BT;

    /** @brief BT-style only: deflated phi_T_mat per rho level.
     *  phi_T_array_BT[rho] is a (p x T_stop) matrix, every column deflated
     *  with the penalty of the last step.
     */
    std::vector<Eigen::MatrixXd> phi_T_array_BT;
};


// ===================================================================================
// TRexDASelector
// ===================================================================================

/**
 * @brief Dependency-aware T-Rex selector (T-Rex+DA): FDR-controlled variable
 *        selection for dependent variables.
 *
 * @details
 *  Inherits from TRexSelector. The selector inserts a dependency-aware deflation
 *  step between the aggregation of the random experiments and the FDP estimation,
 *  through the hooks of the base class.
 *  For BT, NN and PRIOR_GROUPS the calibration is performed over a 3D grid
 *  (T_stop x V x rho_grid); for AR1/EQUI the standard 2D grid is used with modified
 *  phi_T / Phi values.
 *
 *  Reference threshold (BT, NN and PRIOR_GROUPS):
 *  ----------------------------------------------
 *  These methods evaluate every step on the whole grid of correlation thresholds
 *  (rho_grid). The L-loop and the T-loop read the FDP estimate at one of them, the
 *  reference threshold, which sits at three quarters of the grid: element 0.75 times
 *  its length counted from 1, a half rounded up (element 15 of the default 20). The
 *  1.x versions of the R package round a half to the even number and take one
 *  element lower for the grid lengths 6, 14, 22, 30 and so on. The selection then
 *  searches all cells of stopping time, voting level and threshold, and rho_thresh
 *  is the threshold of the cell it takes.
 *  getFDPHatMat(), getPhiMat() and getPhiPrime() therefore hold the values at the
 *  reference threshold, while getRMat(), the voting threshold and the selection
 *  belong to rho_thresh. The values at every threshold are in FDP_hat_array_BT,
 *  Phi_array_BT and R_array_BT of the result, whose FDP_hat_mat and Phi_mat are
 *  empty for these methods.
 *
 *  FDR limitation (BT, NN and PRIOR_GROUPS):
 *  -----------------------------------------
 *  At a very high correlation within the groups (about 0.9 and above) and with no
 *  independent ("white") predictors, the equiangular solver TLARS can exceed the
 *  target FDR, because the true signal and its collinear shadows become
 *  statistically exchangeable. This is inherent to the method, and it also appears
 *  in the 1.x versions of the R package. Mitigate with the default
 *  BTSelectionMode::FeasibleOnly and by preferring a greedy solver (TOMP) in designs
 *  with high correlation.
 *
 *  Penalty of the last step (a deliberate deviation from the papers):
 *  ------------------------------------------------------------------
 *  The dependency-aware relative occurrences of variable j at step t of the T-loop are
 *  Phi^DA_t(j) = Psi(j) * Phi_t(j), with the penalty factor
 *  Psi(j) = 1 / (2 - min |Phi(j) - Phi(j')|) over the variables j' that j is compared with
 *  (its group, its window or all other variables), and Psi(j) = 1/2 for an empty group.
 *
 *  Read literally, the papers below recompute the factor at every step, Psi_t from Phi_t.
 *  This class takes the factor of the last step T of the evaluated path, Psi_T from Phi_T,
 *  and applies it to every step t = 1, ..., T. It does so for all five methods.
 *
 *  Reason: the FDP estimator sums over the increments Phi^DA_t(j) - Phi^DA_(t-1)(j) and is
 *  a count of expected false discoveries only while no increment is negative. One factor
 *  per variable keeps the non-negative increments of Phi_t non-negative. A factor per step
 *  does not, because the factor of a variable changes from step to step. The estimator then
 *  leaves its range, Phi' exceeds 1 and is capped, the variable counts as a certain true
 *  discovery and the FDP estimate is too small.
 *
 *  Evidence: a Monte Carlo comparison of both constructions for all methods (1,000 runs per
 *  setting, block, AR(1) and independent designs, T-LARS, November 2026).
 *  - FDR: the per-step penalty exceeded a target of 0.1 with NN, 0.15 on an independent
 *    design and 0.13 on a block design, where the last-step penalty gave 0.00 and 0.03.
 *    No setting had the last-step penalty above the target and the per-step penalty below.
 *  - Power: the per-step penalty selected more variables only where its estimator was out
 *    of range, and the gain vanished once its history was made monotone. On AR(1) designs
 *    the last-step penalty had the higher TPR, by 0.01 to 0.02.
 *  - Range: Phi' was capped in 1 % to 27 % of the evaluations with the per-step penalty
 *    and in at most 0.01 % with the last-step penalty.
 *
 *  AR1 and EQUI used the per-step penalty up to C++ 1.0.0, R 2.0.1 and Python 1.0.0, as
 *  the 1.x versions of the R package do. BT, NN and PRIOR_GROUPS have always used the
 *  last step, in those versions as well. Do not move a method to a penalty per step, or
 *  add a clamp to Phi', without repeating that comparison (HISTORY.md, 2026-11-13).
 *
 *  References:
 *  -----------
 *   J. Machkour, M. Muma, and D. P. Palomar, "High-Dimensional False Discovery Rate Control for Dependent
 *     Variables," Signal Process., vol. 234, Art. no. 109990, 2025, doi: 10.1016/j.sigpro.2025.109990.
 *   J. Machkour, D. P. Palomar, and M. Muma, "FDR-Controlled Portfolio Optimization for Sparse Financial
 *     Index Tracking," in Proc. IEEE Int. Conf. Acoust., Speech, Signal Process. (ICASSP), 2025, pp. 1-5,
 *     doi: 10.1109/ICASSP49660.2025.10890629.
 *   J. Machkour, M. Muma, and D. P. Palomar, "FDR-Controlled Sparse Index Tracking with Autoregressive Stock
 *     Dependency Models," in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2024, pp. 2662-2666.
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class TRexDASelector : public tc::TRexSelector {

public:

    // ============================================================
    // Result Structure
    // ============================================================

    /**
     * @brief Extended result of T-Rex+DA, adding the record of the dependency
     *        estimation.
     */
    struct DASelectionResult : public SelectionResult {
        /** @brief Correlation threshold of the selected cell of the calibration
         *  grid (BT, NN and PRIOR_GROUPS). -2 (AUTO_ESTIMATE_CORRELATION) for
         *  AR1 and EQUI and when no cell meets the target.
         */
        double rho_thresh = tc::AUTO_ESTIMATE_CORRELATION;

        /** @brief Grid of correlation thresholds of BT, NN and PRIOR_GROUPS:
         *  hc_grid_length ascending values. BT and PRIOR_GROUPS take them from
         *  the merge heights of the dendrogram (1 - height) and end at 1, NN
         *  spaces them equally on [0, 1]. Empty for AR1 and EQUI.
         */
        Eigen::VectorXd rho_grid;

        /** @brief Dependency-aware method of the run. PRIOR_GROUPS when prior
         *  groups were supplied.
         */
        DAMethod method = DAMethod::BT;

        /** @brief Correlation coefficient of the methods AR1 and EQUI: the
         *  supplied value, or the estimate from the data when -2 requested it.
         *  The other methods do not use it and return the value of the control.
         */
        double cor_coef = tc::AUTO_ESTIMATE_CORRELATION;

        /** @brief FDP estimates of the methods BT, NN and PRIOR_GROUPS: one
         *  matrix per correlation threshold (rho_grid), with one row per
         *  stopping time (T_stop) and one column per voting level
         *  (voting_grid). Empty for AR1 and EQUI.
         */
        std::vector<Eigen::MatrixXd> FDP_hat_array_BT;

        /** @brief Dependency-aware relative occurrences of the methods BT, NN
         *  and PRIOR_GROUPS: one matrix per correlation threshold (rho_grid),
         *  with one row per stopping time (T_stop) and one column per variable.
         *  Empty for AR1 and EQUI.
         */
        std::vector<Eigen::MatrixXd> Phi_array_BT;

        /** @brief Numbers of selected variables of the methods BT, NN and
         *  PRIOR_GROUPS: one matrix per correlation threshold (rho_grid), with
         *  one row per stopping time that can be selected (T_stop - 1) and one
         *  column per voting level (voting_grid). Empty for AR1 and EQUI.
         */
        std::vector<Eigen::MatrixXd> R_array_BT;

        /** @brief Virtual destructor. */
        ~DASelectionResult() override = default;
    };


    // ============================================================
    // Constructor / Destructor
    // ============================================================

    /**
     * @brief Construct a TRexDASelector.
     *
     * @param X             Feature matrix (n x p) — Eigen::Map, not copied.
     * @param y             Response vector (n x 1) — Eigen::Map, copied internally.
     * @param tFDR          Target FDR level (default: 0.1).
     * @param trex_da_ctrl  DA-specific control parameters, nesting the base
     *                      T-Rex algorithmic control parameters as
     *                      `trex_da_ctrl.trex_ctrl`.
     * @param seed          Random seed. A negative value (default: -1) gives a run
     *                      that is not reproducible.
     * @param verbose       Enable verbose output (default: true).
     */
    TRexDASelector(
        Eigen::Map<Eigen::MatrixXd>& X,
        Eigen::Map<Eigen::VectorXd>& y,
        double tFDR = 0.1,
        TRexDAControlParameter trex_da_ctrl = TRexDAControlParameter(),
        int  seed    = -1,
        bool verbose = true
    );

    /** @brief Virtual destructor. */
    ~TRexDASelector() override = default;

    /** @brief Deleted copy constructor. */
    TRexDASelector(const TRexDASelector&) = delete;

    /** @brief Deleted copy assignment operator. */
    TRexDASelector& operator=(const TRexDASelector&) = delete;

    /** @brief Deleted move constructor. */
    TRexDASelector(TRexDASelector&&) = delete;

    /** @brief Deleted move assignment operator. */
    TRexDASelector& operator=(TRexDASelector&&) = delete;


    // ============================================================
    // Main Method
    // ============================================================

    /**
     * @brief Retrieve the full DA result from the last select() call.
     *
     * @details `select()` itself is inherited from `TRexSelector`; DA plugs
     * in via the variant hooks `onSelectBegin()` (dependency-structure
     * setup), `evaluateStep()` (deflation), the accumulator overrides, and
     * `finalizeSelectionResult()` (DA selection + result widening).
     *
     * @return Const reference to DASelectionResult.
     */
    const DASelectionResult& getDAResult() const noexcept { return da_result_; }


protected:

    // ============================================================
    // DA Data Members
    // ============================================================

    /** @brief DA-specific control parameters (nests the base trex_ctrl). */
    TRexDAControlParameter trex_da_ctrl_;

    /** @brief Precomputed neighborhood / correlation structure. */
    DASetupResult da_setup_;

    /** @brief Full DA result from the last select() call. */
    DASelectionResult da_result_;


    // ============================================================
    // Per-step side state (populated by evaluateStep override)
    // ============================================================

    /** @brief BT-style snapshot from the most recent evaluateStep call:
     *  deflated Phi per rho level (p x rho_grid_len). Empty for non-BT.
     */
    Eigen::MatrixXd Phi_BT_last_;

    /** @brief BT-style snapshot from the most recent evaluateStep call:
     *  Phi_prime per rho level (p x rho_grid_len). Empty for non-BT.
     */
    Eigen::MatrixXd Phi_prime_BT_last_;

    /** @brief BT-style snapshot from the most recent evaluateStep call:
     *  FDP estimates per rho level (v_len x rho_grid_len). Empty for non-BT.
     */
    Eigen::MatrixXd FDP_hat_BT_last_;

    /** @brief Canonical Phi_prime from the most recent evaluateStep call.
     *  For BT-style this is the column at opt_point_BT; for non-BT this is
     *  the standard 1D Phi_prime computed on the deflated phi_T_mat.
     */
    Eigen::VectorXd Phi_prime_last_;

    // ============================================================
    // BT-style T-loop accumulators (rho_grid_len matrices each)
    // ============================================================

    /** @brief BT-style accumulator: phi_acc_bt_[r] = (capacity x p) of
     *  deflated Phi rows over the T-loop. Allocated only when use_BT_style.
     */
    std::vector<Eigen::MatrixXd> phi_acc_bt_;

    /** @brief BT-style accumulator: fdp_acc_bt_[r] = (capacity x v_len) of
     *  FDP rows over the T-loop. Allocated only when use_BT_style.
     */
    std::vector<Eigen::MatrixXd> fdp_acc_bt_;


    // ============================================================
    // Hook overrides — lifecycle, per-step body, accumulator management
    // ============================================================

    /**
     * @brief Override: build the dependency structure before the L-loop.
     *
     * @details Runs `setupDA()` (clustering / correlation estimation) and
     *  prints the method banner. Called by the inherited `select()` after
     *  the voting-grid setup.
     */
    void onSelectBegin() override;

    /**
     * @brief Override: DA variable selection, result widening, and
     *        `Phi_prime_` correction.
     *
     * @details Replaces the base 2D selection with `selectVariables_BT()`
     *  for BT-style methods (the base `selectedVariables()` result computed
     *  on the canonical opt-rho channel is discarded), re-runs the base 2D
     *  selection semantics for AR1/EQUI on the deflated accumulators, and
     *  overwrites the raw `Phi_prime_` the base pipeline computed with the
     *  deflated one from the last `evaluateStep()`. Populates `da_result_`
     *  and returns the base-sliced result; cleanup and X denormalization
     *  remain with the inherited `select()`.
     */
    SelectionResult finalizeSelectionResult(SelectionResult&& base_result) override;

    /**
     * @brief Override: run experiment, apply DA deflation, compute FDP.
     *
     * @details
     *  BT-style: stashes full 2D matrices (Phi_BT_last_,
     *  Phi_prime_BT_last_, FDP_hat_BT_last_) and writes the column at
     *  `opt_point_BT` into the canonical view (Phi, Phi_prime, FDP_hat).
     *  Non-BT: writes deflated 1D channel into the view.
     *
     *  Always caches `Phi_prime_last_` for later result assembly.
     */
    StepView evaluateStep(std::size_t num_dummies,
                          std::size_t T_stop,
                          bool use_warm_start,
                          er::ExperimentStrategy strategy,
                          std::size_t seed_factor,
                          std::size_t existing_on_disk) override;

    /** @brief Override: also allocate BT accumulators when use_BT_style. */
    void initTAccumulators(std::size_t capacity, std::size_t v_len) override;

    /** @brief Override: also grow BT accumulators in lock-step. */
    void growTAccumulators(std::size_t old_capacity,
                           std::size_t new_capacity) override;

    /** @brief Override: also append BT rows from Phi_BT_last_ /
     *  FDP_hat_BT_last_ for every rho level. */
    void appendTRow(Eigen::Index row_idx, const StepView& view) override;

    /** @brief Override: also trim BT accumulators. */
    void trimTAccumulators(Eigen::Index n_rows) override;


    // ============================================================
    // DA Setup — Neighbourhood / Correlation Dispatch
    // ============================================================

    /**
     * @brief Build the dependency structure from X (once, before the loops).
     *
     * Dispatches to:
     *  - setupDA_PriorGroups() if prior_groups is non-empty.
     *  - setupDA_BT()          for DAMethod::BT.
     *  - setupDA_NN()          for DAMethod::NN.
     *  - setupDA_AR1()         for DAMethod::AR1.
     *  - setupDA_EQUI()        for DAMethod::EQUI.
     */
    void setupDA();

    /** @brief DA-Setup 1: constrained sub-clustering within user-supplied
     *  prior groups. The finest common refinement of all prior_groups levels
     *  constrains the dendrogram (no merges across group boundaries); HAC
     *  runs per group, the rho grid is pooled from the within-group merge
     *  heights (plus the rho = 1 singleton anchor), and gr_j_list holds the
     *  per-rho tight neighbourhoods obtained by cutting each group's
     *  dendrogram — never the raw prior groups themselves.
     */
    void setupDA_PriorGroups();

    /** @brief DA-Setup 2: dendrogram-based clustering on |cor(X)|. */
    void setupDA_BT();

    /** @brief DA-Setup 3: nearest-neighbour correlation threshold sweep. */
    void setupDA_NN();

    /** @brief DA-Setup 4: AR(1) model — estimate ρ and compute window κ. */
    void setupDA_AR1();

    /** @brief DA-Setup 5: equicorrelation — estimate ρ̄. */
    void setupDA_EQUI();


    // ============================================================
    // DA Correction
    // ============================================================

    /**
     * @brief Apply the dependency-aware deflation to experiment results.
     *
     * @details The penalty is computed once from `Phi`, the relative occurrences
     *  of the last step, and deflates every column of `phi_T_mat` (see the
     *  class documentation, "Penalty of the last step").
     *
     * @param phi_T_mat Relative occurrence matrix (p x T_stop).
     * @param Phi       Relative occurrence vector at T_stop (p x 1), the last
     *                  column of `phi_T_mat`.
     * @param T_stop    Current stopping time.
     *
     * @return DACorrectionResult with deflated values.
     */
    DACorrectionResult daCorrect(
        const Eigen::MatrixXd& phi_T_mat,
        const Eigen::VectorXd& Phi,
        std::size_t T_stop) const;


    // ============================================================
    // FDP Computation for BT-style
    // ============================================================

    /**
     * @brief Compute Phi_prime and FDP_hat across the rho_grid for BT-style methods.
     *
     * @param da_corr   DA correction result (contains phi_T_array_BT, Phi_BT).
     * @param num_dummies Number of dummy variables used.
     *
     * @return Pair of (Phi_prime_BT, FDP_hat_BT) where each is a matrix (p x rho_grid_len)
     *         and (V_len x rho_grid_len), respectively.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> computeFDP_BT(
        const DACorrectionResult& da_corr,
        std::size_t num_dummies) const;


    // ============================================================
    // Variable Selection for BT-style
    // ============================================================

    /**
     * @brief Select variables over the 3D grid (T, V, rho_grid) for BT-style methods.
     *
     * @param FDP_hat_array_BT  FDP estimates: one (T_stop x V_len) matrix per rho level.
     * @param Phi_array_BT      Deflated Phi: one (T_stop x p) matrix per rho level.
     * @param voting_grid       Voting grid V.
     * @param rho_grid          Correlation-level grid.
     * @param T_stop            Stopping time.
     *
     * @return DASelectionResult with selected variables, thresholds, and R_array.
     */
    DASelectionResult selectVariables_BT(
        const std::vector<Eigen::MatrixXd>& FDP_hat_array_BT,
        const std::vector<Eigen::MatrixXd>& Phi_array_BT,
        const Eigen::VectorXd& voting_grid,
        const Eigen::VectorXd& rho_grid,
        std::size_t T_stop) const;


    // ============================================================
    // Correlation Estimation Helpers
    // ============================================================

    /** @brief Whether the coefficient is estimated: `cor_coef` is the
     *  sentinel AUTO_ESTIMATE_CORRELATION, or it lies outside [-1, 1], which
     *  also prints a warning.
     */
    bool estimatesCorrelation() const;

    /** @brief Estimate AR(1) coefficient from X via row-wise Yule-Walker. */
    double estimateAR1Correlation() const;

    /** @brief Estimate mean pairwise correlation from X. */
    double estimateEquiCorrelation() const;


    // ============================================================
    // Validation
    // ============================================================

    /** @brief Validate DA-specific parameters. */
    void validateDAParameters() const;


    // ============================================================
    // Clustering Dispatch Helper
    // ============================================================

    /**
     * @brief Run hierarchical clustering with the configured linkage method.
     *
     * @return Sorted merge matrix (p-1 merge steps).
     */
    std::vector<hac::MergeStep> runClustering() const;

    /**
     * @brief Run hierarchical clustering on a column subset of X with the
     *  configured linkage method (correlation distance).
     *
     * @param cols Global column indices forming the subset (size m >= 2).
     *
     * @return Sorted merge matrix (m-1 merge steps) over LOCAL indices
     *  0..m-1 in the order of @p cols.
     */
    std::vector<hac::MergeStep> clusterColumns(
        const std::vector<Eigen::Index>& cols) const;

// -----------------------------------------------------------------------
}; /* End of class TRexDASelector */

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_da */
// ===================================================================================

#endif /* TREX_SELECTOR_METHODS_TREX_DA_HPP */
