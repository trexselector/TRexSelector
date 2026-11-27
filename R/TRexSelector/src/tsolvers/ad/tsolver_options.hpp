// ===================================================================================
// tsolver_options.hpp
// ===================================================================================
/**
 * @file tsolver_options.hpp
 *
 * @brief TSolverOptions — every T-solver input except the data X, D and y.
 *
 * @details One structure carries the preprocessing flags, the numerical
 *  tolerance, the tie-breaking configuration and every family
 *  hyperparameter, so a solver is constructed as
 *  `Solver(X, D, y, options)` (plus the structural inputs that are data,
 *  not options: the group vector of the IEN solvers, the Tikhonov matrix of
 *  the sparse-K solvers). Each family reads the fields it needs and ignores
 *  the rest; the per-field comments say which solvers consume a field. The
 *  defaults are the solver-level defaults (what a direct C++ caller got
 *  before); the T-Rex layer maps its own `TSolverParams` onto this
 *  structure in the dispatch (trex_solver_dispatch.hpp) and keeps its own
 *  defaults there.
 */
// ===================================================================================
#ifndef TREX_TSOLVERS_TSOLVER_OPTIONS_HPP
#define TREX_TSOLVERS_TSOLVER_OPTIONS_HPP

// std includes
#include <cstddef>
#include <limits>

// tsolvers includes
#include <tsolvers/ad/solver_utils/solver_preprocessing.hpp>

// ===================================================================================

namespace trex::tsolvers::ad {

// ===================================================================================

/** @brief Column-scaling convention used by all solvers (unit-L2 or z-score). */
using TSolverScalingMode = solver_utils::preprocessing::TSolverScalingMode;

/** @brief Old name of TSolverScalingMode.
 *  @deprecated Use TSolverScalingMode. */
using ScalingMode [[deprecated("Use TSolverScalingMode")]] = TSolverScalingMode;

/**
 * @brief Variants of the Norm-Corrective GMP solver (TNCGMP).
 */
enum class NCGMPVariant : int {
    LINE_SEARCH = 0,       ///< Standard Matching Pursuit (1D projection)
    FULLY_CORRECTIVE = 1,  ///< Orthogonal Matching Pursuit (least squares over the active set)

    // old names, kept for one release
    LineSearch [[deprecated("Use LINE_SEARCH")]] = LINE_SEARCH,
    FullyCorrective [[deprecated("Use FULLY_CORRECTIVE")]] = FULLY_CORRECTIVE
};

// ===================================================================================

/**
 * @brief Every T-solver input except X, D and y.
 *
 * @note Construction-time only. The execution parameters (T_stop,
 *  early_stop) belong to executeStep(); serialization keeps its own record
 *  of the fields it needs.
 */
struct TSolverOptions {

    // --- Preprocessing (every solver) -------------------------------------

    /** @brief Scale each column of X (and D) to unit L2 norm (or unit SD,
     *  see scaling_mode). Off when the caller pre-normalized (the T-Rex
     *  selectors). Default true. */
    bool normalize = true;

    /** @brief Center X, D and y (an intercept is fitted). Off when the caller
     *  pre-centered. Default true. */
    bool intercept = true;

    /** @brief Whether a solver with these options writes into X and D: it
     *  preprocesses them in place when normalize or intercept is on (a
     *  writer). With both off it only reads them (a reader). */
    bool writesDesign() const noexcept { return normalize || intercept; }

    /** @brief The caller centered X, D and y (intercept = false on centered
     *  data, as in the T-Rex selectors). The centered columns span at most
     *  n - 1 dimensions, so the solver counts the intercept as with
     *  intercept = true: the active set is capped at n - 1 and the degrees
     *  of freedom and Cp include the intercept. The data are not touched.
     *  Redundant when intercept is true. Default false. */
    bool precentered = false;

    /** @brief Column-scaling convention applied when normalize is true:
     *  L2 (unit L2 norm, default) or ZSCORE (sample SD). */
    TSolverScalingMode scaling_mode = TSolverScalingMode::L2;

    /** @brief Print status and diagnostics during execution. Default false. */
    bool verbose = false;

    // --- Threads (every solver) -------------------------------------------

    /** @brief Thread budget of the solver's own OpenMP regions, the passes
     *  over the design (correlation refreshes, the step-length search, the
     *  T-OLS downdate, the standardization and its restore). 1 (default):
     *  serial; at least 1. A region runs in parallel only when the budget
     *  is above 1, the work of the pass exceeds
     *  TSolverBase::kParallelWorkThreshold entries and the solver is not
     *  inside a parallel region of the caller. Not serialized: a loaded
     *  solver starts serial, setNumThreads() gives it a budget. */
    int n_threads = 1;

    // --- Numerics and tie-breaking (every solver) -------------------------

    /** @brief Numerical tolerance of the solver's internal comparisons
     *  (collinearity, step roots, ...). Default: machine epsilon. A selector
     *  sets its own tolerance on the solvers it runs (1e-6 by default). */
    double tol = std::numeric_limits<double>::epsilon();

    /** @brief Seed of the solver's tie-break engine (pruneTiedDummies
     *  shuffle; the exchangeable-tie draw). Default 0: a fixed seed, so a
     *  solver run is reproducible also where exact ties occur. < 0: seeded
     *  from std::random_device, i.e. nondeterministic. */
    long long tie_seed = 0;

    /** @brief Exchangeable-tie band width of the greedy solvers (TOMP, TOOLS,
     *  TGP, TACGP, TMP, TNCGMP and TAFS), in units of the pairwise
     *  ranking-noise sd: statistically indistinguishable top candidates
     *  within highly correlated clusters are picked uniformly at random per
     *  step. 0 (default): off. The path solvers ignore it. */
    double exch_tie_alpha = 0.0;

    /** @brief Minimum |correlation| for exchangeable-tie candidates in
     *  [0, 1). Default 0.5. Ignored unless exch_tie_alpha > 0. */
    double exch_tie_floor = 0.5;

    // --- Family hyperparameters --------------------------------------------

    /** @brief L2 (ridge / group-ridge / Tikhonov) penalty weight, >= 0.
     *  Read by TENET, TENETAug, TIENET, TIENETAug and TIKNET. Default 0.1,
     *  as inside the selectors. With 0 there is no quadratic penalty and the
     *  elastic nets reduce to their lasso. */
    double lambda2 = 0.1;

    /** @brief Shrinkage step of TAFS in (0, 1]. Default 0.3, as inside the
     *  selectors; 1 is forward stepwise. Read by TAFS only. */
    double rho_afs = 0.3;

    /** @brief NCGMP variant. Read by TNCGMP only. Default LineSearch. */
    NCGMPVariant ncgmp_variant = NCGMPVariant::LINE_SEARCH;

    /** @brief TENETAug / TIENETAug: run the augmented system with a pure-LARS
     *  inner path that never drops variables (true) instead of the default
     *  LASSO inner path (false). */
    bool tenet_aug_use_lars = false;

    /** @brief TIKNET: the caller has checked the Tikhonov matrix (finite
     *  entries, symmetric, positive semi-definite), and the solver skips its
     *  own check. For a caller that builds many solvers on one matrix, as
     *  the T-Rex selector does. Default false. */
    bool tikhonov_validated = false;
};

/** @brief Old name of TSolverOptions.
 *  @deprecated Use TSolverOptions. */
using SolverOptions [[deprecated("Use TSolverOptions")]] = TSolverOptions;

// ===================================================================================
} /* End of namespace trex::tsolvers::ad */

#endif /* TREX_TSOLVERS_TSOLVER_OPTIONS_HPP */
