// ===================================================================================
// trex_solver_dispatch.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_UTILS_SOLVER_DISPATCH_HPP
#define TREX_SELECTOR_METHODS_UTILS_SOLVER_DISPATCH_HPP
// ===================================================================================
/**
 * @file trex_solver_dispatch.hpp
 *
 * @brief Type-safe solver dispatch for TRex.
 *
 * @details
 *
 *  Design:
 *  -------
 *    1.  SolverConfig — bundles the parameters every solver needs (data refs,
 *        flags, optionals like lambda2).
 *
 *    2.  dispatchSolver<TSolver>() - a function template that encapsulates the
 *        full create-or-warm-start → execute → extract → optionally-save
 *        lifecycle. A single template handles all solver types uniformly.
 *
 *    3.  dispatchByType() — the one remaining switch/case that maps the
 *        runtime SolverTypeForTRex enum to a compile-time concrete type and
 *        forwards to dispatchSolver<TSolver>.
 *
 *  Adding a new solver:
 *    1.  Write the solver class (inheriting TSolver_Base) with the standard
 *        constructor, default constructor, save(), and static load().
 *    2.  Add one `case` line in dispatchByType().
 */
// ===================================================================================

// std includes
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Solver includes
#include <tsolvers/tsolvers.hpp>

// ===================================================================================

// Embedded into namespace trex::trex_selector_methods::ad::utils::solver_dispatch
namespace trex::trex_selector_methods::ad::utils::solver_dispatch {

// ===================================================================================

/** @brief Sparse beta path handed from the solvers to the T-Rex aggregation
 *  (per-step support + values; O(sum_t |support_t|) memory instead of the
 *  former dense (p + num_dummies) x steps matrix). */
using SparseBetaPath = trex::tsolvers::ad::SparseBetaPath;

// ===================================================================================


// ===================================================================================
// Solver type enum for TRex dispatch
// ===================================================================================

/**
 * @brief Solver types available for the T-Rex Selector.
 */
enum class SolverTypeForTRex {

    // LARS-based solvers
    // -------------------------------
    TLARS,      // Terminating LARS
    TLASSO,     // Terminating LASSO
    TENET,      // Terminating ENET
    TENET_AUG,  // Terminating ENET via augmented LASSO (GVS only)
    TIENET,     // Terminating Informed ENET, native pathwise (GVS only)
    TIENET_AUG, // Terminating Informed ENET via augmented LASSO (GVS only)
    TIKNET,     // Terminating Tikhonov Net, native pathwise (sparse-K, FOLDED)
    TSTEPWISE,  // Terminating Stepwise
    TSTAGEWISE, // Terminating Stagewise
    // -------------------------------

    // OMP-based solvers
    // -------------------------------
    TOMP,       // Terminating OMP
    TGP,        // Terminating Gradient Pursuit
    TACGP,      // Terminating Approximate Conjugate Gradient Pursuit
    TMP,        // Terminating Matching Pursuit
    TNCGMP,     // Terminating Norm-Corrected Generalized Matching Pursuit
    TOOLS,      // Terminating Orthogonal Least Squares
    // -------------------------------

    // AFS-based solvers
    // --------------------------------
    TAFS
    // --------------------------------
};


// ===================================================================================
// Encapsulated Solver Hyperparameters
// ===================================================================================
/**
 * @brief Bundles all specific hyperparameters required by individual solver algorithms.
 *
 * @details Adding a new parameter here automatically makes it available to the
 * dispatcher without touching intermediate routing layers.
 */
struct SolverHyperparameters {
    /** @brief L2 penalty for ENET-type solvers (ignored by non-using solvers). */
    double lambda2 = 0.1;   //
    /** @brief Shrinkage step size for AFS solver in (0, 1] (ignored by non-using solvers). */
    double rho_afs = 0.3;   //
    /** @brief Variant selector for NCGMP solver
     * (0: Line Search/MP, 1: Fully Corrective/OMP) (ignored by non-using solvers). */
    int ncgmp_variant = 0;  //
    /** @brief Numerical tolerance of the selector: the threshold of its own
     *  comparisons and the tolerance of the solvers in its random experiments.
     *  A stand-alone solver defaults to machine epsilon (SolverOptions::tol).
     *  The two tolerances serve different purposes and differ on purpose. */
    double tol = 1e-6;      //
    /** @brief Exchangeable-tie band width of the greedy solvers (TOMP,
     * TOOLS, TGP, TACGP, TMP, TNCGMP and TAFS), in units of the pairwise
     * ranking-noise sd. When > 0, statistically indistinguishable top
     * candidates within highly correlated clusters are picked uniformly at
     * random per step, restoring the within-experiment occurrence spread
     * that the DA deflation's FDR control relies on.
     * -1 (default): the selector decides. The T-Rex+DA selector runs with
     * 0.25, every other selector without the band. 0 = off, also under
     * T-Rex+DA. The path solvers ignore the option. */
    double exch_tie_alpha = -1.0;
    /** @brief Minimum |correlation| for exchangeable-tie candidates in
     * (0, 1) (ignored unless exch_tie_alpha > 0). */
    double exch_tie_floor = 0.5;

    /** @brief TENET_AUG: run the augmented system with a pure-LARS inner
     * path that never drops variables (true) instead of the default LASSO
     * inner path that can (false). Same knob as the GVS
     * tenet_aug_use_lars; ignored by every other solver. */
    bool tenet_aug_use_lars = false;
};


// ===================================================================================
// SolverConfig — everything a solver needs to be created or warm-started
// ===================================================================================

/**
 * @brief Bundles all parameters required to construct or warm-start any TSolver.
 *
 * @note References to X, D, y must outlive the SolverConfig (non-owning).
 */
struct SolverConfig {

    // --- Data references (non-owning, must outlive config) ---

    /** @brief Reference to the original design matrix X (n x p). */
    Eigen::Map<Eigen::MatrixXd>& X;

    /** @brief Reference to the dummy design matrix D (n x num_dummies). */
    Eigen::Map<Eigen::MatrixXd>& D;

    /** @brief Reference to the response vector y (n x 1). */
    Eigen::Map<Eigen::VectorXd>& y;

    // --- Solver construction flags ---

    /** @brief Whether the solver should normalize columns internally.
     *  In TRex context always false (TRex normalizes X/D beforehand).
     */
    bool normalize = false;

    /** @brief Whether the solver should fit an intercept internally.
     *  In TRex context always false (TRex centers y beforehand).
     */
    bool intercept = false;

    /** @brief The caller centered X, D and y (SolverOptions::precentered):
     *  the solver caps its active set at n - 1. In TRex context true for
     *  every ScalingMode except NONE.
     */
    bool precentered = false;

    /** @brief Enable verbose output from the solver. */
    bool verbose   = false;

    /** @brief Internal scaling mode applied by the solver to its own columns.
     *  In TRex context X/D are pre-normalized, so this is usually L2; set to
     *  ZSCORE only when the caller wants the solver to z-score internally
     *  (e.g. augmented EN paths that must keep the ridge block consistent).
     */
    trex::tsolvers::ad::ScalingMode scaling_mode = trex::tsolvers::ad::ScalingMode::L2;

    // --- Execution parameters ---

    /** @brief Early stopping threshold (number of active dummies). */
    std::size_t T_stop = 0;

    /** @brief If true, stop when count_active_dummies >= T_stop. */
    bool early_stop = true;

    // --- Warm-start / serialization ---

    /** @brief If true, attempt to load from solver_file before executing. */
    bool use_warm_start = false;

    /**
     * @brief Path to the serialized solver state file. Empty string means "no serialization".
     */
    std::string solver_file{};

    // --- Solver-specific extras ---
    /** @brief Encapsulated solver hyperparameters for t-solvers. */
    SolverHyperparameters hyperparams{};

    /** @brief Deterministic tie-break seed for the solver's internal RNG
     *  (pruneTiedDummies shuffle). < 0 (the T-Rex default, an unseeded
     *  selector): std::random_device seeding (nondeterministic). Applied only on fresh
     *  construction — warm-started / deserialized solvers keep their RNG
     *  state so the shuffle stream continues instead of restarting. */
    long long tie_seed = -1;

    /** @brief Thread budget of the solver's own OpenMP regions
     *  (SolverOptions::n_threads): 1 serial, at least 1. Applied at
     *  construction and re-applied to warm-started and deserialized solvers,
     *  since the budget is not serialized. Inside a parallel K-loop the
     *  solver keeps its regions serial whatever the budget says. */
    int n_threads = 1;

    // --- In-memory warm start ---

    /** @brief Retained solver to continue in-memory (nullptr = none).
     *  When set, dispatchSolver() skips construction entirely and resumes
     *  this solver via executeStep(T_stop, early_stop). The solver holds
     *  non-owning views into X and its dummy matrix from the run that
     *  created it; cfg.X / cfg.D are ignored on this path. */
    tsolvers::ad::TSolver_Base* warm_solver = nullptr;

    /** @brief Retention sink for the freshly constructed solver (nullptr = none).
     *  When set, dispatchSolver() heap-allocates the solver and moves it here
     *  after execution so the caller can retain it for the next T-step.
     *  The caller must keep the X and D buffers alive for as long as the
     *  solver is retained. */
    std::unique_ptr<tsolvers::ad::TSolver_Base>* retain_sink = nullptr;

    // --- General-Tikhonov penalty (TIKNET dispatch) ---

    /** @brief Non-owning pointer to a general Tikhonov matrix K = Gamma^T
     *  Gamma (p x p, PSD; nullptr = none). Required by the general-K
     *  solver TIKNET (used by TRexTikhonovSelector).
     *  Must outlive the dispatch call; on serialized warm starts K travels
     *  inside the solver state instead. */
    const Eigen::SparseMatrix<double>* tikhonov_K = nullptr;

    /** @brief tikhonov_K was checked by the caller (finite, symmetric,
     *  positive semi-definite), so the solver skips its own check. */
    bool tikhonov_K_validated = false;

};


// Aliases for clean reading in the if constexpr block
namespace lars = trex::tsolvers::ad::linear_model::lars_based;
namespace omp  = trex::tsolvers::ad::linear_model::omp_based;
namespace afs  = trex::tsolvers::ad::linear_model::afs_based;


/**
 * @brief Translate a SolverConfig into the solver-level SolverOptions.
 *
 * @details The T-Rex layer keeps its own SolverHyperparameters (lambda2 0.1
 *          and rho_afs 0.3 as the solvers' own defaults, its own tol 1e-6).
 *          This is the single place where they are mapped onto the solvers'
 *          construction-time options.
 */
inline trex::tsolvers::ad::SolverOptions toSolverOptions(const SolverConfig& cfg) {
    trex::tsolvers::ad::SolverOptions o;
    const SolverHyperparameters& hp = cfg.hyperparams;
    o.normalize      = cfg.normalize;
    o.intercept      = cfg.intercept;
    o.precentered    = cfg.precentered;
    o.verbose        = cfg.verbose;
    o.scaling_mode   = cfg.scaling_mode;
    o.tol            = hp.tol;
    o.tie_seed       = cfg.tie_seed;
    o.exch_tie_alpha = hp.exch_tie_alpha;
    o.exch_tie_floor = hp.exch_tie_floor;
    o.lambda2        = hp.lambda2;
    o.rho_afs        = hp.rho_afs;
    o.ncgmp_variant  = static_cast<trex::tsolvers::ad::NCGMPVariant>(hp.ncgmp_variant);
    o.tenet_aug_use_lars = hp.tenet_aug_use_lars;
    o.n_threads      = cfg.n_threads;
    o.tikhonov_validated = cfg.tikhonov_K_validated;
    return o;
}


/**
 * @brief Construct a solver on the heap from the config.
 *
 * @details Every solver takes (X, D, y, options). The sparse-K solver TIKNET
 *          takes its sparse Tikhonov matrix in addition, which is data and
 *          not an option. The options come from toSolverOptions(). The
 *          tolerance, the tie-break seed and the exchangeable-tie
 *          configuration are applied at construction, so no
 *          post-construction setter is needed on this path (the warm-start
 *          paths still apply them to a restored solver).
 */
template <typename TSolver>
std::unique_ptr<TSolver> makeSolverForConfig(const SolverConfig& cfg) {
    const trex::tsolvers::ad::SolverOptions options = toSolverOptions(cfg);

    if constexpr (std::is_same_v<TSolver, lars::TIKNET_Solver>) {
        // TIKNET implements the FOLDED dummy convention: dummy copy t
        // couples as its base (t - p) mod p.
        if (cfg.tikhonov_K == nullptr) {
            throw std::invalid_argument(
                "trex_solver_dispatch: the sparse-K solver requires a "
                "Tikhonov matrix (SolverConfig::tikhonov_K).");
        }
        return std::make_unique<TSolver>(
            cfg.X, cfg.D, cfg.y, *cfg.tikhonov_K, options);
    } else {
        return std::make_unique<TSolver>(cfg.X, cfg.D, cfg.y, options);
    }
}


/**
 * @brief Execute a single solver lifecycle: construct (or warm-start), run, extract path,
 *        and optionally serialize.
 *
 * @details Three paths, checked in order:
 *          1. In-memory warm start (cfg.warm_solver != nullptr): the retained
 *             solver is resumed via executeStep(T_stop, early_stop); no
 *             construction, no serialization.
 *          2. Serialized warm start (cfg.use_warm_start && solver_file set):
 *             the solver is deserialized from disk, reconnected to cfg.X /
 *             cfg.D, continued, and re-saved.
 *          3. Fresh construction via makeSolverForConfig(). If a solver_file
 *             is set, the state is saved for later resumption; if
 *             cfg.retain_sink is set, ownership of the solver is handed to
 *             the caller for in-memory retention.
 *
 * @tparam TSolver Concrete solver type (must satisfy the TSolver_Base interface).
 *
 * @param cfg Populated SolverConfig with data references and execution flags.
 *
 * @return Sparse beta path ((p + num_dummies)-space support per step).
 */
template <typename TSolver>
SparseBetaPath dispatchSolver(const SolverConfig& cfg) {

    // 1. In-memory warm start: resume the retained solver. executeStep()
    //    continues the existing path up to the new T_stop; getBetaPathSparse()
    //    returns the full accumulated path.
    if (cfg.warm_solver != nullptr) {
        cfg.warm_solver->setTolerance(cfg.hyperparams.tol);
        cfg.warm_solver->setExchangeableTie(cfg.hyperparams.exch_tie_alpha,
                                            cfg.hyperparams.exch_tie_floor);
        cfg.warm_solver->setNumThreads(cfg.n_threads);
        cfg.warm_solver->executeStep(cfg.T_stop, cfg.early_stop);
        return cfg.warm_solver->getBetaPathSparse();
    }

    // 2. Serialized warm start: restore from disk, continue, re-save.
    if (cfg.use_warm_start && !cfg.solver_file.empty()) {
        TSolver solver = TSolver::load(cfg.solver_file, cfg.X, cfg.D);
        solver.setTolerance(cfg.hyperparams.tol);
        solver.setExchangeableTie(cfg.hyperparams.exch_tie_alpha,
                                  cfg.hyperparams.exch_tie_floor);
        solver.setNumThreads(cfg.n_threads);   // not serialized: a loaded solver starts serial
        solver.executeStep(cfg.T_stop, cfg.early_stop);
        SparseBetaPath path = solver.getBetaPathSparse();
        solver.save(cfg.solver_file);
        return path;
    }

    // 3. Fresh construction: tolerance, tie seed and exchangeable ties are
    //    part of the SolverOptions (toSolverOptions()).
    std::unique_ptr<TSolver> solver = makeSolverForConfig<TSolver>(cfg);
    solver->executeStep(cfg.T_stop, cfg.early_stop);
    SparseBetaPath path = solver->getBetaPathSparse();
    if (!cfg.solver_file.empty()) solver->save(cfg.solver_file);
    if (cfg.retain_sink != nullptr) *cfg.retain_sink = std::move(solver);
    return path;
}


/**
 * @brief Map a runtime SolverTypeForTRex enum to a concrete solver type and dispatch.
 *
 * @details The sole switch/case in the dispatch layer. Each case instantiates the
 *          corresponding solver type via dispatchSolver<TSolver>(cfg).
 *          To add a new solver: write the class, add one case here.
 *
 * @param solver_type  Runtime enum value identifying the desired solver.
 * @param cfg          Populated SolverConfig with data references and execution flags.
 *
 * @return Sparse beta path ((p + num_dummies)-space support per step).
 *
 * @throws std::invalid_argument if solver_type is not a recognized enumerator.
 */
inline SparseBetaPath dispatchByType(SolverTypeForTRex solver_type, const SolverConfig& cfg) {
    switch (solver_type) {

        // LARS family
        case SolverTypeForTRex::TLARS:      return dispatchSolver<lars::TLARS_Solver>(cfg);
        case SolverTypeForTRex::TLASSO:     return dispatchSolver<lars::TLASSO_Solver>(cfg);
        case SolverTypeForTRex::TSTEPWISE:  return dispatchSolver<lars::TSTEPWISE_Solver>(cfg);
        case SolverTypeForTRex::TENET:      return dispatchSolver<lars::TENET_Solver>(cfg);
        case SolverTypeForTRex::TENET_AUG: return dispatchSolver<lars::TENETAug_Solver>(cfg);
        case SolverTypeForTRex::TIENET:     throw std::invalid_argument(
                "trex_solver_dispatch: TIENET is a GVS-only solver (needs a "
                "group assignment); use TRexGVSSelector with gvs_type=IEN.");
        case SolverTypeForTRex::TIENET_AUG: throw std::invalid_argument(
                "trex_solver_dispatch: TIENET_AUG is a GVS-only solver; "
                "use TRexGVSSelector with gvs_type=IEN.");
        case SolverTypeForTRex::TIKNET:
            // Dispatchable only with a Tikhonov matrix (TRexTikhonovSelector
            // sets cfg.tikhonov_K; with penalty_structure GROUP_CENTERING
            // the selector materializes K_GC, which this sparse path
            // consumes directly).
            if (cfg.tikhonov_K != nullptr) {
                return dispatchSolver<lars::TIKNET_Solver>(cfg);
            }
            throw std::invalid_argument(
                "trex_solver_dispatch: TIKNET needs a Tikhonov matrix "
                "(TRexTikhonovSelector / SolverConfig::tikhonov_K).");
        case SolverTypeForTRex::TSTAGEWISE: return dispatchSolver<lars::TSTAGEWISE_Solver>(cfg);

        // OMP family
        case SolverTypeForTRex::TOMP:       return dispatchSolver<omp::TOMP_Solver>(cfg);
        case SolverTypeForTRex::TGP:        return dispatchSolver<omp::TGP_Solver>(cfg);
        case SolverTypeForTRex::TACGP:      return dispatchSolver<omp::TACGP_Solver>(cfg);
        case SolverTypeForTRex::TMP:        return dispatchSolver<omp::TMP_Solver>(cfg);
        case SolverTypeForTRex::TNCGMP:     return dispatchSolver<omp::TNCGMP_Solver>(cfg);
        case SolverTypeForTRex::TOOLS:      return dispatchSolver<omp::TOOLS_Solver>(cfg);

        // AFS family
        case SolverTypeForTRex::TAFS:       return dispatchSolver<afs::TAFS_Solver>(cfg);

        default:
            throw std::invalid_argument("trex_solver_dispatch: unsupported SolverTypeForTRex");
    }
}

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::utils::solver_dispatch */
// ===================================================================================
#endif /* End of TREX_SELECTOR_METHODS_UTILS_SOLVER_DISPATCH_HPP */
