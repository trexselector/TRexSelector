// ===================================================================================
// trex.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_CORE_HPP
#define TREX_SELECTOR_METHODS_TREX_CORE_HPP
// ===================================================================================
/**
 * @file trex.hpp
 *
 * @brief T-Rex Selector for FDR-controlled variable selection.
 */
// ===================================================================================

// std includes
#include <cstddef>
#include <exception>
#include <functional>
#include <optional>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// OpenMP compatibility layer
#include <utils/openmp/utils_openmp.hpp>

// TRex includes
#include <trex_selector_methods/ad/trex_utils/trex_data_normalizer.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_solver_dispatch.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_warm_start_manager.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_memmap_manager.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_dummy_generator.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_experiment_runner.hpp>

// Dummy generation distribution type
#include <utils/datageneration/utils_dummygen.hpp>

// ===================================================================================

// Embedded into namespace trex::trex_selector_methods::ad::trex_core
namespace trex::trex_selector_methods::ad::trex_core {

// Namespace alias
namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;
namespace dg = trex::trex_selector_methods::ad::utils::dummy_generator;
namespace wsm = trex::trex_selector_methods::ad::utils::warm_start_manager;
namespace mm = trex::trex_selector_methods::ad::utils::memmap_manager;
namespace disp = trex::trex_selector_methods::ad::utils::solver_dispatch;
namespace er = trex::trex_selector_methods::ad::utils::experiment_runner;
namespace dummygen = trex::utils::datageneration::dummygen;

// ===================================================================================
// Constants
// ===================================================================================

/**
 * @brief Sentinel value for "estimate this correlation coefficient from X".
 *
 * @details A value outside [−1, 1], which no correlation coefficient takes.
 */
inline constexpr double AUTO_ESTIMATE_CORRELATION = -2.0;

/**
 * @brief Whether a solver is greedy (non-equiangular forward selection).
 *
 * @details Greedy solvers tend to run into a noise trap in the T-loop: they
 *          keep admitting variables until the maximum number of dummies is
 *          selected instead of terminating. The path solvers are free of this
 *          behavior: TLARS, TLASSO, TENET, TENET_AUG, TIENET, TIENET_AUG
 *          and TIKNET. Every other solver
 *          counts as greedy: TSTEPWISE, TSTAGEWISE, TOMP, TGP, TACGP, TMP,
 *          TNCGMP, TOOLS and TAFS. Used to resolve the AUTO default of
 *          TRexControl::tloop_stagnation_stop and for the
 *          corresponding validation warning.
 */
inline constexpr bool isGreedySolver(
    trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex st) noexcept
{
    namespace sdx = trex::trex_selector_methods::ad::utils::solver_dispatch;
    switch (st) {
        case sdx::SolverTypeForTRex::TLARS:
        case sdx::SolverTypeForTRex::TLASSO:
        case sdx::SolverTypeForTRex::TENET:
        case sdx::SolverTypeForTRex::TENET_AUG:
        case sdx::SolverTypeForTRex::TIENET:
        case sdx::SolverTypeForTRex::TIENET_AUG:
        case sdx::SolverTypeForTRex::TIKNET:
            return false;
        default:
            return true;
    }
}

// ===================================================================================
// Enums & Control Structures
// ===================================================================================

/**
 * @brief L-loop calibration strategies.
 *
 * @details Two orthogonal axes:
 *  - Dummy source: fresh INDEPENDENT draws per experiment (SKIPL / STANDARD /
 *    HCONCAT / SEEDED) vs one shared base matrix whose rows are permuted per
 *    experiment (PERMUTATION / PERMUTATION_SEEDED).
 *  - Storage: STORED strategies keep matrices in the DummyGenerator for the
 *    whole run; the SEEDED strategies re-derive everything from the seed at
 *    each step (zero persistent state). SEEDED is the storage-free twin of
 *    STANDARD: the same (k, L step) keys, so the two draw bit-identical
 *    dummies for the same seed. PERMUTATION_SEEDED is the storage-free twin
 *    of PERMUTATION (prefix-stable base, deriveBlockSeed64 l_tag = 0).
 */
enum class LLoopStrategy {
    SKIPL,                 // Stored, one-shot: fixed dummies = max_dummy_multiplier * p.
    STANDARD,              // Stored, fresh dummies at each L-loop iteration (default).
    HCONCAT,               // Stored, the matrices grow by p columns per L step (faster, same FDR).
    PERMUTATION,           // Stored base dummy matrix + deterministic row permutations per k.
    PERMUTATION_SEEDED,  // Seed-derived base + row permutations per k; nothing stored.
    SEEDED               // STANDARD's draws re-derived per k and L step; nothing stored.
};


/**
 * @brief Algorithmic control parameters for T-Rex Selector.
 */
struct TRexControl {

    // ====================================
    // Experiment Configuration
    // ====================================

    /** @brief Number of random experiments (default: 20, at least 2). */
    std::size_t K = 20;

    /** @brief Maximum number of dummy variables as a multiple of p (default: 10). */
    std::size_t max_dummy_multiplier = 10;

    /** @brief If true, limit T_stop to ceiling(n/2) (default: true). If false,
     *  the T-loop can run up to the number of dummies. */
    bool use_max_T_stop = true;

    /** @brief Reference voting level of the L-loop (default: 0.75, in [0, 1]).
     *  The L-loop raises L until the FDP estimate at the largest voting level
     *  that does not exceed opt_threshold is at or below tFDR. */
    double opt_threshold = 0.75;

    // ===================================
    // Dummy Generation
    // ===================================

    /** @brief Distribution for generating dummies (default: Normal). */
    dummygen::DummyDistribution dummy_distribution = dummygen::DummyDistribution::Normal();

    // ===================================
    // L-Loop Strategy
    // ===================================

    /** @brief Strategy for L-loop calibration (default: STANDARD).
     *
     *  STANDARD draws a fresh set of dummies at every L-loop iteration, as
     *  the 1.x versions do (random_experiments / add_dummies_GVS). HCONCAT
     *  keeps the dummies of the earlier iterations and appends p columns per
     *  iteration, which is faster. The two gave the same FDR in simulations
     *  and can select differently. LLoopStrategy lists the other strategies.
     */
    LLoopStrategy lloop_strategy = LLoopStrategy::STANDARD;

    // ==================================
    // T-Loop Early Stopping
    // ==================================

    /** @brief T-loop early stopping when the selection stagnates.
     *
     *  Tri-state: std::nullopt (default, "auto") resolves by solver family in
     *  the TRexSelector constructor. It is off for the path solvers (TLARS,
     *  TLASSO, TENET, TENET_AUG, TIENET, TIENET_AUG and TIKNET), which
     *  terminate on their own, as in the
     *  1.x versions and the T-Rex paper, neither of which has a stagnation
     *  stop. It is on for the greedy solvers (TSTEPWISE, TSTAGEWISE, TOMP,
     *  TGP, TACGP, TMP, TNCGMP, TOOLS, TAFS), which otherwise run into a
     *  noise trap and select until the maximum number of dummies is reached.
     *  Set an explicit true/false to override.
     */
    std::optional<bool> tloop_stagnation_stop = std::nullopt;

    /** @brief Number of T steps in a row without a growing selection that
     *  ends the T-loop under the stagnation stop (default: 5, at least 1).
     *  The selection is counted at the highest voting level whose FDP
     *  estimate does not exceed tFDR.
     */
    std::size_t tloop_max_stagnant_steps = 5;

    // =================================
    // Threads: one owner of the cores at a time
    // =================================

    /** @brief Run the K random experiments in parallel (default: false, the
     *  experiments run in sequence). Inside a parallel K-loop every solver
     *  runs serially; nested teams are never formed. The lean L-loop
     *  strategies (SEEDED, PERMUTATION, PERMUTATION_SEEDED) always run their
     *  experiments in sequence. */
    bool parallel_rnd_experiments = false;

    /** @brief Threads of the parallel K-loop (default: 1, at least 1); the
     *  count is capped at K. Ignored while parallel_rnd_experiments is
     *  false. */
    int max_outer_threads = 1;

    /** @brief Thread budget of the solvers' own regions, passed as
     *  TSolverOptions::n_threads (default: 1, serial; at least 1). Live when
     *  the experiments run in sequence; inside a parallel K-loop the solvers
     *  stay serial whatever the budget says. */
    int max_inner_threads = 1;

    /** @brief Validate the two thread options (a value below 1 throws) and
     *  cap the outer count at K. The selector constructor calls it on its
     *  copy of the control. */
    void resolveThreadOptions() {
        auto check = [](int value, const char* name) {
            if (value < 1) {
                throw std::invalid_argument(
                    std::string(name) + " must be >= 1. Got: " + std::to_string(value));
            }
        };
        check(max_outer_threads, "max_outer_threads");
        check(max_inner_threads, "max_inner_threads");
        if (K > 0) {
            max_outer_threads = std::min(max_outer_threads, static_cast<int>(K));
        }
    }

    // ===================================
    // Memory-Mapping
    // ===================================

    /** @brief Keep the dummy matrices in memory-mapped files instead of in
     *  memory (default: false). */
    bool use_memory_mapping = false;

    // ===================================
    // Solver Configuration
    // ===================================

    /** @brief Solver type (default: TLARS). */
    disp::SolverTypeForTRex solver_type = disp::SolverTypeForTRex::TLARS;

    /** @brief Bundled hyperparameters for all solver algorithms. */
    disp::TSolverParams solver_params{};

    // ===================================
    // Data Scaling
    // ===================================

    /** @brief Working-design convention. ScalingMode::L2 (default) centers X
     *  and y and divides each column of X by its L2 norm; ScalingMode::ZSCORE
     *  centers and divides by the sample SD (n-1) to unit variance (glmnet
     *  `standardize=TRUE` convention). For the plain LASSO path the two are
     *  equivalent (uniform column rescaling), but they differ for the
     *  Elastic-Net ridge term — see TRexGVSSelector. ScalingMode::NONE leaves
     *  X and y exactly as supplied (no centering, no scaling): the caller
     *  pre-centers and pre-scales as the model demands, or supplies an
     *  intercept-free dictionary (steering, Vandermonde, Fourier atoms). The
     *  dummies keep the unit-L2 standard, so a NONE design should be
     *  unit-L2-commensurate; a design pre-normalized with
     *  data_normalizer::centerAndL2NormalizeX / centerY selects under NONE
     *  exactly what the raw design selects under L2. Non-finite entries and
     *  dead columns are rejected in every mode (std::invalid_argument). L2
     *  and ZSCORE standardize X in place and restore it, so a memory-mapped
     *  X needs AccessMode::ReadWrite. NONE writes nothing into X, so
     *  AccessMode::ReadOnly serves (data_normalizer::writesDesign()). */
    dn::ScalingMode scaling_mode = dn::ScalingMode::L2;

};

/** @brief Old name of TRexControl.
 *  @deprecated Use TRexControl. */
using TRexControlParameter [[deprecated("Use TRexControl")]] = TRexControl;


// ===================================================================================
// TRexSelector
// ===================================================================================

/** @brief Result from variable selection. */
struct TRexResult {

    /** @brief Binary support vector (0/1, length p). */
    Eigen::VectorXi selected_var;

    /** @brief Voting threshold applied to Phi to obtain the selection. */
    double v_thresh;

    /** @brief Number of selected variables per stopping time (rows) and
     *  voting level (columns), (T_stop - 1) × v_len. The stopping time
     *  T_stop is not part of it. One row for T_stop = 1. */
    Eigen::MatrixXd R_mat;

    /** @brief Stopping time T_stop: the number of included dummies at
     *  which the T-loop stopped, so the random experiments ran until
     *  T_stop dummies were included. For T_stop > 1 it is not the
     *  stopping time of the selection: the selection is taken from the
     *  stopping times 1 to T_stop - 1, mostly from T_stop - 1 and
     *  sometimes from an earlier one (see select()). */
    std::size_t T_stop;

    /** @brief Total number of dummy variables used (L · p). */
    std::size_t num_dummies;

    /** @brief Dummy multiplier L determined by the L-loop. */
    std::size_t dummy_factor_L;

    /** @brief Estimated FDP per stopping time (rows) and voting level
     *  (columns), T_stop × v_len. */
    Eigen::MatrixXd FDP_hat_mat;

    /** @brief Relative occurrences of the variables over the K random
     *  experiments per stopping time (rows), T_stop × p. */
    Eigen::MatrixXd Phi_mat;

    /** @brief Calibrated relative occurrences Phi' that enter the FDP
     *  estimate (p × 1), at most 1 and possibly negative. */
    Eigen::VectorXd Phi_prime;

    /** @brief Voting levels 0.5, 0.5 + 1/K, ... below 1, closed by
     *  1 - eps. */
    Eigen::VectorXd voting_grid;

    virtual ~TRexResult() = default;
};

/**
 * @brief T-Rex Selector for FDR-controlled variable selection.
 *
 * @details
 *  Base class implementing the T-Rex algorithm.
 *  Supports in-memory and memory-mapped data via Eigen::Map.
 *  Base class of the selectors
 *  - T-Rex+DA (TRexDASelector)
 *  - T-Rex+GVS (TRexGVSSelector)
 *  - T-Rex+Screening (TRexScreeningSelector)
 *  - T-Rex+Tikhonov (TRexTikhonovSelector)
 *
 *  T-Rex+SPCA (TRexSPCASelector) and the biobank screening (TRexBiobankScreeningSelector) run
 *  selectors of these classes.
 *
 *  Key design:
 *    - Does not own X (p >> n, too large to copy).
 *    - Normalizes X in-place and restores it, copies y and centers the copy.
 *    - Delegates dummy generation, solver dispatch, warm-start management,
 *      memory mapping, and experiment coordination to helper modules.
 *    - select() returns TRexResult with all key outputs.
 *
 *  References:
 *  -----------
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class TRexSelector {
protected:

    // ============================================================
    // Data Members — Input Data
    // ============================================================

    /** @brief Non-owning pointer to design matrix X (n x p). */
    Eigen::Map<Eigen::MatrixXd>* X_;

    /** @brief Owned response vector y (n x 1). */
    Eigen::VectorXd y_{};

    /** @brief Number of observations. */
    const std::size_t n_;

    /** @brief Number of features. */
    const std::size_t p_;

    // ============================================================
    // Data Members — Configuration
    // ============================================================

    /** @brief Target False Discovery Rate. */
    const double tFDR_;

    /** @brief Algorithmic control parameters. */
    const TRexControl trex_ctrl_;

    /** @brief Random seed. */
    const int seed_;

    /** @brief Verbose output flag. */
    const bool verbose_;

    /** @brief Numerical tolerance.
     *  @note  Initialized to 0 here; the constructor always overwrites this
     *         with `solver_params.tol` (default 1e-6) or, when tol <= 0,
     *         with `std::numeric_limits<double>::epsilon()`. */
    double eps_{0.0};

    // ============================================================
    // Data Members — Normalization (via trex_data_normalizer)
    // ============================================================

    /** @brief Normalization parameters for X and y. */
    dn::NormalizationParams norm_params_;

    /** @brief Tracks whether X is currently in normalized state.
     *  It guards against double denormalization (select() + destructor).
     */
    bool X_is_normalized_{false};

    /** @brief True after restore() or a failed run gave X back: the
     *  selector's normalized view is gone and a run must refuse. */
    bool released_{false};

    /** @brief True from the start of the first run. A selector runs once. */
    bool has_run_{false};

    /** @brief Entry of a run. Throws std::runtime_error if X was given back
     *  before (restore(), a failed run) or if the selector has run already,
     *  and marks the selector as run. Every run calls it first: select() of
     *  the base and of the variants, and a variant's own run. */
    void beginRun();

    /** @brief Scope guard of the select() implementations: an exception that
     *  leaves the scope gives X back (original scale, claim released) and
     *  retires the selector, so the caller's design never stays standardized
     *  after an error. */
    class GiveXBackOnError {
    public:
        explicit GiveXBackOnError(TRexSelector& selector) noexcept
            : selector_(selector), pending_(std::uncaught_exceptions()) {}
        ~GiveXBackOnError() {
            if (std::uncaught_exceptions() > pending_) { selector_.giveXBackAfterError(); }
        }
        GiveXBackOnError(const GiveXBackOnError&) = delete;
        GiveXBackOnError& operator=(const GiveXBackOnError&) = delete;

    private:
        TRexSelector& selector_;
        int pending_;
    };

    /** @brief What GiveXBackOnError runs: restore X and retire the selector.
     *  Never throws, the destructor stays the last line of defense. */
    void giveXBackAfterError() noexcept;

    /** @brief Restore X to its original scale and release the shared-buffer
     *  claim. Idempotent (no-op if X is already denormalized). Subclasses that
     *  restore X at the end of their own select() override must call this
     *  instead of denormalizing directly, so the buffer claim is released. */
    void denormalizeAndReleaseX();

    // ============================================================
    // Data Members — Results
    // ============================================================

    /** @brief Binary support vector (0/1). */
    Eigen::VectorXi selected_var_;

    /** @brief Indices of selected variables. */
    std::vector<std::size_t> selected_indices_;

    /** @brief Stopping time at which the T-loop stopped (TRexResult::T_stop). */
    std::size_t T_stop_{0};

    /** @brief Number of dummy variables used. */
    std::size_t num_dummies_{0};

    /** @brief Voting grid used. */
    Eigen::VectorXd voting_grid_;

    /** @brief Voting threshold. */
    double voting_threshold_{0.0};

    /** @brief Estimated FDP matrix. */
    Eigen::MatrixXd FDP_hat_mat_;

    /** @brief Relative occurrence matrix. */
    Eigen::MatrixXd Phi_mat_;

    /** @brief Numbers of selected variables (see TRexResult::R_mat). */
    Eigen::MatrixXd R_mat_;

    /** @brief Adjusted relative occurrences (p x 1). */
    Eigen::VectorXd Phi_prime_;

    /** @brief Optimization point for dummy calibration. */
    std::size_t opt_point_;

    /** @brief Deflated Phi (length p) from the last L-loop `evaluateStep`.
     *
     *  Captured at the end of the L-loop so that the T-loop seed row (T = 1)
     *  carries the *deflated* Phi rather than the raw experiment Phi. For the
     *  base (non-deflating) selector this equals `exp_results.Phi`; for
     *  dependency-aware variants it is the deflated calibration-level Phi,
     *  as the seed row `Phi_mat <- matrix(Phi, nrow = 1)` of the 1.x versions. */
    Eigen::VectorXd lloop_last_Phi_;

    /** @brief Current dummy multiplier in L-loop. */
    std::size_t dummy_multiplier_LL_{0};

    /** @brief Cached base seed for the PERMUTATION L-loop strategy.
     *
     *  Resolved once in the constructor so that `prepareDummiesForLStep`
     *  observes the same value across every L-iteration regardless of
     *  whether `seed_ >= 0` (deterministic) or `seed_ < 0` (random_device).
     *  64-bit to stay in the collision-free seeding pipeline (mix_seed64).
     */
    std::uint64_t permutation_base_seed_{0};

    // ============================================================
    // Helper Module Members
    // ============================================================

    /** @brief Dummy generation (all strategies). */
    dg::DummyGenerator dummy_gen_;

    /** @brief Warm-start management (in-memory or serialized). */
    wsm::WarmStartManager warm_start_mgr_;

    /** @brief Memory-mapped D file management (nullptr if not used). */
    std::unique_ptr<mm::MemmapManager> memmap_mgr_;

    /** @brief Experiment coordination (K-experiment loop). */
    std::unique_ptr<er::ExperimentRunner> experiment_runner_;


public:

    // ============================================================
    // Result Structure
    // ============================================================

    /** @brief Old name of TRexResult.
     *  @deprecated Use TRexResult. */
    using SelectionResult [[deprecated("Use TRexResult")]] = TRexResult;


    // ============================================================
    // Constructor / Destructor
    // ============================================================

    /**
     * @brief Construct a TRexSelector object.
     *
     * @param X Design matrix (n x p) — Eigen::Map, not copied.
     * @param y Response vector (n x 1) — Eigen::Map, copied internally.
     * @param tFDR Target False Discovery Rate in [0, 1] (default: 0.1).
     * @param trex_control Algorithmic control parameters.
     * @param seed Random seed (default: -1). A seed >= 0 makes the run
     *        deterministic and reproducible: the same seed draws the same
     *        dummies. A negative seed draws them from the hardware entropy
     *        source, and the run is not reproducible bit by bit. A Monte Carlo
     *        study gives every trial a seed of its own, or uses -1.
     * @param verbose Enable verbose output (default: true).
     *
     * @throws std::invalid_argument for an argument or a control value out of
     *         range, a response of the wrong length, a non-finite entry of X
     *         or y, or a dead column of X.
     * @throws std::runtime_error if another live selector or solver holds X.
     */
    TRexSelector(Eigen::Map<Eigen::MatrixXd>& X,
                 Eigen::Map<Eigen::VectorXd>& y,
                 double tFDR = 0.1,
                 TRexControl trex_control = TRexControl(),
                 int seed = -1,
                 bool verbose = true);

    /** @brief Virtual destructor. */
    virtual ~TRexSelector();

    /** @brief Delete copy constructor */
    TRexSelector(const TRexSelector&) = delete;

    /** @brief Delete copy assignment operator. */
    TRexSelector& operator=(const TRexSelector&) = delete;

    /** @brief Delete move constructor. */
    TRexSelector(TRexSelector&&) = delete;

    /** @brief Delete move assignment operator. */
    TRexSelector& operator=(TRexSelector&&) = delete;


    // ============================================================
    // Main Method
    // ============================================================

    /**
     * @brief Run T-Rex selection algorithm.
     *
     * @details A selector runs once. select() gives X back when it ends, so
     *  a second call throws std::runtime_error instead of running on a design
     *  that is no longer standardized: construct a new selector to run
     *  again. The results of the run stay readable.
     *
     *  The T-loop raises the stopping time until the FDP estimate at the
     *  highest voting level exceeds tFDR, or until its bound or the
     *  stagnation stop ends it, and T_stop is where it stopped.
     *  The selection is the cell with the largest number of selected
     *  variables among the stopping times 1 to T_stop - 1 (the stopping time
     *  1 for T_stop = 1) and the voting levels whose FDP estimate does not
     *  exceed tFDR. Among equal counts the highest voting level wins, then
     *  the latest stopping time. The stopping time of the selection is
     *  therefore below T_stop, mostly by one and sometimes by more.
     *
     *  A select() that ends with an exception gives X back at once (original
     *  scale, claim released) and retires the selector as well.
     *
     * @return TRexResult containing selected variables and
     *         diagnostics.
     */
    virtual TRexResult select();

    /**
     * @brief Give X back before select(): undo the in-place normalization
     *        applied at construction and release the buffer claim.
     *
     * @details The constructor normalizes X in place; select() restores it
     *          on return and the destructor restores it as a safety net. A
     *          caller who constructs a selector and decides not to run it can
     *          hand X back right away with restore() instead of waiting for
     *          the destructor (in R and Python: the garbage collector).
     *          Idempotent; a no-op after select(). A restored selector cannot
     *          run: select() throws, construct a new one.
     */
    void restore();


    // ============================================================
    // Public Getters — Configuration
    // ============================================================

    /** @brief Get target FDR level. */
    double getTFDR() const noexcept { return tFDR_; }

    /** @brief Get number of random experiments K. */
    std::size_t getK() const noexcept { return trex_ctrl_.K; }

    /** @brief Get the cap of the dummy multiplier (max_dummy_multiplier): the
     *  L-loop adds dummies up to this multiple of p. The value is a multiplier,
     *  not a number of dummies. */
    std::size_t getMaxNumDummies() const noexcept { return trex_ctrl_.max_dummy_multiplier; }

    /** @brief Check whether the bound on the stopping time is on
     *  (use_max_T_stop): the T-loop then ends at ceiling(n/2) at the latest.
     *  The value is a flag, not a stopping time. */
    bool getMaxTStop() const noexcept { return trex_ctrl_.use_max_T_stop; }

    /** @brief Check if stagnation check is enabled (AUTO is resolved by
     *  solver family in the constructor, so this is always definite). */
    bool getStagnationCheck() const noexcept {
        return trex_ctrl_.tloop_stagnation_stop.value_or(false);
    }

    // ============================================================
    // Public Getters — Normalization
    // ============================================================

    /** @brief Get means of the design matrix columns (zeros under ScalingMode::NONE). */
    const Eigen::VectorXd& getXMeans() const noexcept { return norm_params_.X_means; }

    /** @brief Get the scale divisors of the design matrix columns: the L2
     *  norms of the centered columns under ScalingMode::L2, the sample
     *  standard deviations under ZSCORE, ones under NONE. */
    const Eigen::VectorXd& getXL2Norms() const noexcept { return norm_params_.X_l2norms; }

    /** @brief Get mean of the response vector (0 under ScalingMode::NONE). */
    double getYMean() const noexcept { return norm_params_.y_mean; }

    // ============================================================
    // Public Getters — Results
    // ============================================================

    /** @brief Get selected variables as a binary vector. */
    const Eigen::VectorXi& getSelectedVar() const noexcept { return selected_var_; }

    /** @brief Get indices of selected variables. */
    const std::vector<std::size_t>& getSelectedIndices() const noexcept {
        return selected_indices_; }

    /** @brief Get number of selected variables. */
    std::size_t getNumSelected() const noexcept { return selected_indices_.size(); }

    /** @brief Number of dummies used by the last select() (the L-loop's calibration). */
    std::size_t getNumDummies() const noexcept { return num_dummies_; }

    /** @brief Get the stopping time T_stop at which the T-loop stopped. The
     *  selection is taken from an earlier one (TRexResult::T_stop). */
    std::size_t getTStop() const noexcept { return T_stop_; }

    /** @brief Get voting grid vector V = [0.5, 0.5+1/K, ..., 1-eps]. */
    const Eigen::VectorXd& getVotingGrid() const noexcept { return voting_grid_; }

    /** @brief Get voting threshold. */
    double getVotingThreshold() const noexcept { return voting_threshold_; }

    /** @brief Get FDP hat matrix. */
    const Eigen::MatrixXd& getFDPHatMat() const noexcept { return FDP_hat_mat_; }

    /** @brief Get Phi matrix. */
    const Eigen::MatrixXd& getPhiMat() const noexcept { return Phi_mat_; }

    /** @brief Get the numbers of selected variables (TRexResult::R_mat). */
    const Eigen::MatrixXd& getRMat() const noexcept { return R_mat_; }

    /** @brief Get Phi prime vector. */
    const Eigen::VectorXd& getPhiPrime() const noexcept { return Phi_prime_; }

    /** @brief Get dummy multiplier L. */
    std::size_t getDummyMultiplierL() const noexcept { return dummy_multiplier_LL_; }

protected:

    // ============================================================
    // Core T-Rex Computation — Voting & FDP
    // ============================================================

    /** @brief Generate voting grid V = [0.5, 0.5 + 1/K, ..., 1-eps] for the
     *  configured number of random experiments K. */
    Eigen::VectorXd createVotingGrid() const;

    /** @brief The same grid for an arbitrary number of experiments K (a
     *  sub-run with fewer experiments, e.g. a lambda_2 pilot, aggregates
     *  its Phi at resolution 1/K and needs the matching grid). */
    Eigen::VectorXd createVotingGrid(std::size_t K) const;

    /**
     * @brief Set optimization point.
     *
     * @param v_len Length of the voting grid.
     * @param threshold Voting threshold to determine optimization point (default: 0.75).
     */
    void setPercentOptPoint(std::size_t v_len, double threshold = 0.75);

    /**
     * @brief Compute Phi_prime (calibration parameter).
     *
     * @param phi_T_mat T x p matrix of relative occurrences.
     * @param Phi p x 1 vector of relative occurrences at T_stop.
     * @param num_dummies Number of dummy variables used (for scaling).
     *
     * @return Phi-prime vector (p x 1) for use in FDP estimation.
     */
    Eigen::VectorXd computePhiPrime(
        const Eigen::MatrixXd& phi_T_mat,
        const Eigen::VectorXd& Phi,
        std::size_t num_dummies) const;

    /**
     * @brief Compute FDP estimates for voting grid.
     *
     * @param voting_grid Voting grid vector.
     * @param Phi p x 1 vector of relative occurrences at T_stop.
     * @param Phi_prime p x 1 vector of calibrated occurrences.
     * @return FDP estimates for each voting threshold.
     */
    Eigen::VectorXd computeFDPHat(
        const Eigen::VectorXd& voting_grid,
        const Eigen::VectorXd& Phi,
        const Eigen::VectorXd& Phi_prime) const;

    /**
     * @brief Select variables with FDR control.
     *
     * @param FDP_hat_mat Matrix of FDP estimates.
     * @param Phi_mat Matrix of relative occurrences.
     * @param voting_grid Voting grid vector.
     * @param T_stop Stopping time.
     * @param drop_last_row Drop the last T-row when T_stop > 1 (default, the
     *        convention of the 1.x versions for the main run). After an exit
     *        by the FDP rule that row is the step that broke the rule. After
     *        an exit at T_max or by the stagnation stop it is a feasible step
     *        and is dropped all the same, which mirrors the 1.x versions. A
     *        run over a FIXED T budget, such as a lambda_2 pilot, passes
     *        false: none of its rows is a boundary step.
     *
     * @return Selection result containing selected variables and related information.
     */
    TRexResult selectedVariables(
        const Eigen::MatrixXd& FDP_hat_mat,
        const Eigen::MatrixXd& Phi_mat,
        const Eigen::VectorXd& voting_grid,
        std::size_t T_stop,
        bool drop_last_row = true) const;

    /**
     * @brief The T-Rex selection rule applied to ONE runner result: builds
     *        the (T x p) Phi and (T x v) FDP-hat accumulators row by row
     *        from `phi_T_mat` (T = its number of columns) with
     *        computePhiPrime() / computeFDPHat(), then selectedVariables().
     *
     * @details The path to T dummies has every shorter path as a prefix, so
     *   `phi_T_mat.leftCols(t + 1)` is what a run stopped at t + 1 dummies
     *   would return, and row t of the accumulators equals the main T-loop's
     *   row t. Used by the lambda_2 pilot of the Tikhonov selector; a GVS
     *   pilot would call it the same way.
     *
     * @param phi_T_mat   p x T relative-occurrence matrix of one runner call.
     * @param num_dummies Number of dummies that run used (Phi-prime scaling).
     * @param voting_grid The voting grid matching the run's number of experiments.
     * @param drop_last_row Forwarded to selectedVariables().
     */
    TRexResult selectFromPhiT(
        const Eigen::MatrixXd& phi_T_mat,
        std::size_t num_dummies,
        const Eigen::VectorXd& voting_grid,
        bool drop_last_row) const;


    // ============================================================
    // Per-Iteration Step + Accumulator Hooks (overridable)
    // ============================================================

    /**
     * @brief Per-L-iteration dummy-preparation context.
     *
     * @details
     *  Passed to `prepareDummiesForLStep()` once per L-loop iteration. The
     *  hook is expected to populate / refresh whatever per-experiment dummy
     *  state the subsequent `evaluateStep()` call will consume, and to write
     *  back `existing_on_disk` for HCONCAT memmap accounting.
     */
    struct LStepContext {
        /** Current dummy multiplier (1-based; equals max_dummy_multiplier
         *  for the SKIPL one-shot path). */
        std::size_t L_iter = 1;

        /** Total number of dummy columns this iteration (`L_iter * p_`). */
        std::size_t num_dummies = 0;

        /** Active L-loop strategy. */
        LLoopStrategy strategy = LLoopStrategy::STANDARD;

        /** Output: number of dummy columns already persisted on disk.
         *  Used by the HCONCAT memmap path to skip rewriting prefix
         *  columns; set to 0 by all other strategies. */
        std::size_t existing_on_disk = 0;
    };

    /**
     * @brief Lifecycle hook fired once at the start of `select()`, after
     *        voting-grid setup and before the L-loop.
     *
     * @details
     *  Default base implementation: no-op. Subclasses may override to perform
     *  one-shot setup (e.g. precomputing cluster structure) that the L-loop
     *  / T-loop drivers depend on.
     */
    virtual void onSelectBegin();

    /**
     * @brief Per-L-iteration dummy-preparation hook.
     *
     * @details
     *  Default base implementation: dispatches on `ctx.strategy` and performs
     *  the canonical T-Rex dummy work via `dummy_gen_` (and invalidates
     *  `warm_start_mgr_` when fresh dummy data has been produced). Subclasses
     *  that bypass `dummy_gen_` (e.g. cluster-aware MVN draws in GVS) may
     *  override this hook to populate their own per-experiment dummy state
     *  while still inheriting the L-loop control flow.
     *
     *  Subclasses MUST set `ctx.existing_on_disk` to the appropriate value
     *  (typically 0 unless they explicitly implement an HCONCAT memmap
     *  variant).
     *
     * @param ctx Per-iteration context (input + output).
     */
    virtual void prepareDummiesForLStep(LStepContext& ctx);

    /**
     * @brief Final result post-processing hook.
     *
     * @details
     *  Called from `select()` after `selectedVariables()` and `storeResults()`,
     *  but before warm-start / memmap cleanup and before X is denormalized.
     *  Subclasses may override to widen the returned object (e.g. into a
     *  derived result struct) and to perform variant-specific cleanup of
     *  any extra resources allocated in `onSelectBegin()` or in the L-loop.
     *
     *  Default base implementation: returns `base_result` unchanged.
     *
     * @param base_result The base TRexResult produced by the standard pipeline.
     * @return The (possibly widened) TRexResult to return from `select()`.
     */
    virtual TRexResult finalizeSelectionResult(TRexResult&& base_result);

    /**
     * @brief Bundle of per-step outputs shared by the L-loop and T-loop drivers.
     *
     * @details
     *  Subclasses (e.g. dependency-aware variants) override `evaluateStep` to
     *  inject extra processing — typically a deflation step — between the
     *  experiment runner and the FDP computation.
     *  They populate `FDP_hat` with the canonical 1D channel that the base
     *  drivers rely on for the L-loop reference-point guard, the T-loop
     *  edge guard, and stagnation tracking. Side data (e.g. per-rho 2D matrices)
     *  is stashed in subclass-private member state.
     */
    struct StepView {
        /** @brief Raw experiment-runner output (Phi, phi_T_mat, etc.). */
        er::ExperimentResults exp_results;

        /** @brief Canonical FDP estimates for the voting grid (length v_len). */
        Eigen::VectorXd FDP_hat;

        /** @brief Canonical Phi vector (length p). For DA-BT variants this is
         *  the deflated Phi at the calibration rho-level. */
        Eigen::VectorXd Phi;

        /** @brief Canonical Phi_prime vector (length p). */
        Eigen::VectorXd Phi_prime;
    };

    /**
     * @brief Run one experiment step and produce the canonical FDP/Phi view.
     *
     * @details
     *  Default base implementation: runs the experiment runner, computes
     *  `Phi_prime = computePhiPrime(...)` and `FDP_hat = computeFDPHat(...)`.
     *  Subclasses may override to insert a deflation step.
     *
     *  Subclasses MUST populate `view.FDP_hat`, `view.Phi`, and
     *  `view.Phi_prime` with vectors that the L-loop reference-point guard,
     *  the T-loop edge guard, and the base accumulator append logic can use.
     *
     * @param num_dummies        Number of dummies for this step.
     * @param T_stop             Stopping time for this step.
     * @param use_warm_start     Pass-through to the experiment runner.
     * @param strategy           Experiment strategy.
     * @param seed_factor        Seed factor for the runner config.
     * @param existing_on_disk   `existing_cols_on_disk` for HCONCAT memmap.
     * @return Populated StepView.
     */
    virtual StepView evaluateStep(
        std::size_t num_dummies,
        std::size_t T_stop,
        bool use_warm_start,
        er::ExperimentStrategy strategy,
        std::size_t seed_factor,
        std::size_t existing_on_disk);

    /**
     * @brief Initialise T-loop accumulators with the given capacity.
     *
     * @details
     *  Default base: zero-initialises `Phi_mat_` (capacity x p) and
     *  `FDP_hat_mat_` (capacity x v_len). Subclasses may extend to allocate
     *  their own per-rho accumulators.
     */
    virtual void initTAccumulators(std::size_t capacity, std::size_t v_len);

    /**
     * @brief Grow T-loop accumulators to a larger capacity (capacity-doubling).
     *
     * @details
     *  Default base: `conservativeResize` on `Phi_mat_` / `FDP_hat_mat_` and
     *  zero the newly added bottom rows.
     */
    virtual void growTAccumulators(std::size_t old_capacity,
                                   std::size_t new_capacity);

    /**
     * @brief Append the current step's row at the given index.
     *
     * @details
     *  Default base: writes `view.Phi` and `view.FDP_hat` into row `row_idx`
     *  of `Phi_mat_` and `FDP_hat_mat_`.
     */
    virtual void appendTRow(Eigen::Index row_idx, const StepView& view);

    /**
     * @brief Trim T-loop accumulators to the actually populated rows.
     *
     * @details Default base: `conservativeResize` to `n_rows` rows.
     */
    virtual void trimTAccumulators(Eigen::Index n_rows);


    // ============================================================
    // L-Loop Calibration Strategies
    // ============================================================

    /**
     * @brief Dispatch L-loop strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoop(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);

    /**
     * @brief Per-iteration policy context for the unified L-loop driver.
     *
     * Bundles the strategy enum, a flag telling the driver whether to forward
     * the current `dummy_multiplier_LL_` as `seed_factor` (vs. fixed 0), and
     * the strategy-specific dummy-generation callable.
     */
    struct LLoopPolicyContext {
        /** Experiment strategy to forward to ExperimentRunnerConfig. */
        er::ExperimentStrategy strategy;
        /** If true, pass `dummy_multiplier_LL_` as `seed_factor`; else 0. */
        bool seed_factor_uses_LL;
        /**
         * @brief Strategy-specific per-iteration dummy work.
         * @param LL Current value of `dummy_multiplier_LL_` (1-based).
         * @param num_dummies `LL * p_` for this iteration.
         * @param existing_on_disk_out Output: columns already on memmap disk
         *        (used by HCONCAT for the runner config); set to 0 if N/A.
         *
         * Contract: the callable must perform any required calls on
         * `dummy_gen_` (e.g. `generateAndStore`, `expandStored`,
         * `storeBaseDummies`) and invoke `warm_start_mgr_.invalidate()`
         * whenever it mutates dummy data that previously cached solvers
         * depend on.
         */
        std::function<void(std::size_t LL,
                           std::size_t num_dummies,
                           std::size_t& existing_on_disk_out)> policy;
    };

    /**
     * @brief Unified L-loop calibration driver shared by the five looping
     *        strategies (STANDARD, HCONCAT, PERMUTATION, SEEDED,
     *        PERMUTATION_SEEDED).
     *
     * Implements the common while-loop skeleton:
     *   - guard `dummy_multiplier_LL_ <= max_dummy_multiplier &&
     *           (FDP_hat empty || FDP_hat(opt_point_) > tFDR_)`,
     *   - compute `num_dummies_ = dummy_multiplier_LL_ * p_`,
     *   - delegate strategy-specific dummy work to `ctx.policy`,
     *   - build runner config + run experiments,
     *   - update `Phi_prime` and `FDP_hat` via `computePhiPrime` /
     *     `computeFDPHat`,
     *   - increment `dummy_multiplier_LL_`.
     *
     * @note SKIPL is intentionally NOT routed through this driver: it has no
     *       loop, performs a single shot with `num_dummies_ =
     *       max_dummy_multiplier * p_`, and never increments
     *       `dummy_multiplier_LL_`.
     *
     * @param ctx Strategy policy context.
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoopCalibration(const LLoopPolicyContext& ctx,
                             Eigen::VectorXd& FDP_hat,
                             er::ExperimentResults& exp_results);

    /**
     * @brief Run L-loop calibration using the SKIPL strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoopCalibration_SKIPL(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);

    /**
     * @brief Run L-loop calibration using the STANDARD strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoopCalibration_Standard(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);

    /**
     * @brief Run L-loop calibration using the HCONCAT strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results  Current experiment results (updated in-place).
     */
    void runLLoopCalibration_HCONCAT(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);

    /**
     * @brief Run L-loop calibration using the PERMUTATION strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoopCalibration_Permutation(Eigen::VectorXd& FDP_hat,
                                         er::ExperimentResults& exp_results);

    /**
     * @brief Run L-loop calibration using the SEEDED or the PERMUTATION_SEEDED strategy.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runLLoopCalibration_Seeded(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);


    // ============================================================
    // T-Loop Calibration
    // ============================================================

    /** @brief Run T-loop to determine stopping threshold T_stop.
     *
     * @param FDP_hat Current FDP estimates (updated in-place).
     * @param exp_results Current experiment results (updated in-place).
     */
    void runTLoop(Eigen::VectorXd& FDP_hat, er::ExperimentResults& exp_results);


    // ============================================================
    // Utilities
    // ============================================================

    /** @brief Validate input data dimensions. */
    void validateTRexParameters() const;

    /**
    * @brief Print progress message (if verbose).
    *
    * @param message Message to print.
    */
    void printProgress(const std::string& message) const;

    /**  @brief Update selected_indices_ from selected_var_. */
    void updateSelectedIndices();

    /**
     * @brief Store results after selection.
     *
     * @param result Selection result to store.
     */
    void storeResults(const TRexResult& result);

    /**
     * @brief Build an ExperimentRunnerConfig from current TRex state.
     *
     * @param num_dummies Number of dummies for this iteration.
     * @param T_stop Stopping threshold.
     * @param use_warm_start Whether to use warm start.
     * @param strategy Experiment strategy.
     * @param seed_factor L-loop iteration multiplier for seed mixing (default: 0).
     * @param existing_cols_on_disk Columns already on disk for HCONCAT memmap path (default: 0).
     *
     * @note Virtual so selector variants can post-process the config
     *       (TRexTikhonovSelector injects its Tikhonov matrix and the
     *       CV-resolved lambda_2).
     *
     * @return Populated ExperimentRunnerConfig.
     */
    virtual er::ExperimentRunnerConfig buildRunnerConfig(
        std::size_t num_dummies,
        std::size_t T_stop,
        bool use_warm_start,
        er::ExperimentStrategy strategy,
        std::size_t seed_factor = 0,
        std::size_t existing_cols_on_disk = 0) const;

// -----------------------------------------------------------------------
}; /* End of class TRexSelector */

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_core */
// ===================================================================================
#endif /* TREX_SELECTOR_METHODS_TREX_CORE_HPP */
