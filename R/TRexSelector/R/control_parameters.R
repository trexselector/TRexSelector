# =============================================================================
# Control parameter constructor functions for TRex selector variants.
#
# These mirror the Cpp struct hierarchy:
#   trex_control()        -> TRexControlParameter
#   trex_da_control()     -> TRexDAControlParameter
#   trex_gvs_control()    -> GVS control struct
#   trex_screen_control() -> ScreenTRexControlParameter
#   trex_biobank_control()-> BiobankScreenTRex FDR bounds
# =============================================================================

#' @name trex_control
#'
#' @title TRex Control Parameters
#'
#' @description Build a control list for TRex selector algorithms.
#'
#' @note **Seed recommendation for Monte Carlo FDR simulations:** always pass
#'   \code{seed = -1} (the default) to \code{TRexSelector$new()} or any other
#'   selector constructor when running repeated trials. A fixed integer seed
#'   causes every trial to draw K per-experiment dummy seeds as
#'   \code{mix_seed(fixed, 0)}, \code{mix_seed(fixed, 1)}, \ldots — deterministic but
#'   structurally correlated across trials — which can inflate the empirical FDR
#'   in a sweep. With \code{seed = -1} each of the K experiments independently
#'   draws hardware entropy via \code{std::random_device}, matching the
#'   behaviour of the standalone C++ demo.
#'
#' @param solver Solver algorithm: "TLARS", "TLASSO", "TENET", "TENET_AUG",
#'   "TSTAGEWISE", "TSTEPWISE", "TOMP", "TGP", "TACGP", "TMP", "TNCGMP",
#'   "TOOLS", "TAFS", "TCDLASSO", "TCDENET" (default: "TLARS"). "TENET_AUG"
#'   is the row-augmented formulation of the TENET elastic net (proved
#'   equivalent for \code{lambda2 >= 0}) and is accepted wherever TENET
#'   runs; it solves with the standard LASSO inner path unless the optional
#'   \code{tenet_aug_use_lars} is set. The coordinate-descent (CD) family
#'   (TCDLASSO/TCDENET) records exact penalized-lasso minimizers per dummy
#'   crossing. The IEN-family names "TIENET", "TIENET_AUG", and "TCDIENET"
#'   are accepted by the parser but only dispatchable inside
#'   \code{TRexGVSSelector} with \code{gvs_type = "IEN"}; the
#'   Tikhonov-family names "TCDTIKNET" and "TIKNET" need a
#'   penalty matrix and are dispatchable only inside
#'   \code{TRexTikhonovSelector}. The plain selector rejects both families
#'   with an informative error. The coordinate-descent names "TCDLASSO",
#'   "TCDENET", "TCDIENET" and "TCDTIKNET" are deprecated and will be removed
#'   in TRexSelector 2.1.0: a selector that receives one of them raises a
#'   deprecation warning. Use "TLASSO", "TENET", "TIENET" and "TIKNET", the
#'   homotopy solvers of the same problems.
#' @param K Number of random experiments (default: 20).
#' @param tloop_stagnation_stop Whether to stop the T-loop on stagnation. NULL
#'   (default) leaves the choice to the C++ core, the same rule on every
#'   surface: FALSE for the path solvers (TLARS, TLASSO, TENET, TENET_AUG,
#'   TIENET, TIENET_AUG, TIKNET and the CD family TCDLASSO, TCDENET, TCDIENET,
#'   TCDTIKNET), TRUE for the greedy ones (the OMP family, TSTEPWISE,
#'   TSTAGEWISE, TAFS, TNCGMP, TOOLS).
#' @param tloop_max_stagnant_steps Steps before stagnation trigger (default: 5).
#' @param parallel_rnd_experiments Run the K random experiments in parallel
#'   (default: \code{FALSE}, the experiments run in sequence). Inside a
#'   parallel run every solver runs serially; the lean L-loop strategies
#'   (\code{"SEEDED"}, \code{"PERMUTATION"}, \code{"PERMUTATION_SEEDED"})
#'   always run their experiments in sequence. Was \code{use_openmp}.
#' @param use_memory_mapping Use out-of-core memory mapping (default: FALSE).
#' @param max_outer_threads Threads of the parallel K-loop (default: 1; a
#'   whole number of at least 1), capped at \code{K}. Ignored while
#'   \code{parallel_rnd_experiments} is \code{FALSE}.
#' @param max_inner_threads Thread budget of the solvers' own passes over
#'   the design (default: 1, serial; a whole number of at least 1). Live
#'   when the experiments run in sequence; inside a parallel run the
#'   solvers stay serial. A pass splits only above about a million entries
#'   of the design, so small problems run serially whatever the budget.
#' @param max_dummy_multiplier Cap of the dummy multiplier: the L-loop adds dummies up to
#'   this many times the number of variables (default: 10).
#' @param use_max_T_stop Whether the T-loop ends at \code{ceiling(n / 2)} at the latest
#'   (default: TRUE). Without the bound it can run up to the number of dummies.
#' @param opt_threshold Optimization threshold (default: 0.75).
#' @param lloop_strategy L-loop strategy: "SKIPL", "STANDARD", "HCONCAT",
#'   "PERMUTATION", "PERMUTATION_SEEDED", or "SEEDED" (default: "STANDARD",
#'   as in C++ and Python).
#'   The SEEDED variants re-derive dummies from the seed at every step and
#'   store nothing: "SEEDED" is the storage-free twin of "STANDARD" (fresh
#'   dummies at every L step, bit-identical to "STANDARD" for the same seed)
#'   and "PERMUTATION_SEEDED" the storage-free twin of "PERMUTATION".
#' @param tol Numerical tolerance of the selector (default: 1e-6): the threshold of
#'   its own comparisons and the tolerance of the solvers in its random experiments.
#'   A stand-alone solver defaults to machine epsilon (\code{\link{tsolver_options}}).
#'   The two tolerances serve different purposes and differ on purpose.
#' @param lambda2 Elastic-net L2 penalty (TENET hyperparameter, default: 0.1).
#' @param rho_afs Shrinkage step of TAFS in (0, 1] (default: 0.3); 1 is forward stepwise.
#' @param ncgmp_variant NCGMP variant: 0 = line_search, 1 = fully_corrective (default: 0).
#' @param exch_tie_alpha Exchangeable-tie band width for greedy solvers
#'   (TOMP/TAFS), in units of the pairwise ranking-noise sd. When > 0,
#'   statistically indistinguishable top candidates within highly correlated
#'   clusters are picked uniformly at random per step, restoring the
#'   within-experiment occurrence spread that the dependency-aware (DA)
#'   deflation's FDR control relies on. Recommended for greedy solvers under
#'   trex+DA: 0.25. \code{0} (default) turns the feature off (exact legacy
#'   behavior; ignored by path solvers).
#' @param exch_tie_floor Minimum absolute correlation for exchangeable-tie
#'   candidates, in (0, 1) (default: 0.5). Ignored unless
#'   \code{exch_tie_alpha > 0}.
#' @param tenet_aug_use_lars Logical, for \code{solver = "TENET_AUG"} in the
#'   plain/DA/screening selectors. The default \code{FALSE} runs the
#'   STANDARD LASSO inner path (variables can be dropped); \code{TRUE}
#'   opts into the pure-LARS inner path that never drops variables. Inside
#'   \code{TRexGVSSelector} use the \code{gvs_control} knob of the same
#'   name instead. Ignored by every other solver.
#' @param cd_lambda_rel_tol Relative lambda step tolerance of the CCD family's
#'   crossing search (default: 1e-6, as in C++ and Python). Values <= 0 keep
#'   the solver's internal default. Ignored by the non-CCD solvers.
#'   Deprecated: it will be removed with the coordinate-descent family in
#'   TRexSelector 2.1.0.
#' @param cd_tol_certify KKT certification tolerance of the CCD family. Must
#'   be set together with \code{cd_tol_probe}; the default -1 keeps the
#'   solver's own tolerances. Ignored by the non-CCD solvers. Deprecated: it
#'   will be removed with the coordinate-descent family in TRexSelector 2.1.0.
#' @param cd_tol_probe Probing tolerance of the CCD family's coordinate
#'   sweeps. Must be set together with \code{cd_tol_certify}; the default -1
#'   keeps the solver's own tolerances. Ignored by the non-CCD solvers.
#'   Deprecated: it will be removed with the coordinate-descent family in
#'   TRexSelector 2.1.0.
#' @param cd_gram_cap Column cap for the CCD family's cached Gram matrix;
#'   above it the solver falls back to the naive updates (default: 0 = keep
#'   the solver's internal default). Ignored by the non-CCD solvers.
#'   Deprecated: it will be removed with the coordinate-descent family in
#'   TRexSelector 2.1.0.
#' @param cd_max_sweeps Maximum coordinate sweeps per CCD subproblem
#'   (default: 0 = keep the solver's internal default of 2000). A solve that
#'   stops at this budget far from a solution (KKT violation above 0.001
#'   lambda_1) raises an error instead of continuing on an unconverged path.
#'   Ignored by the non-CCD solvers. Deprecated: it will be removed with the
#'   coordinate-descent family in TRexSelector 2.1.0.
#' @param cd_certified_localization Logical (default \code{TRUE}): the CD
#'   family localizes every dummy crossing with certificates — a duality
#'   gap and coordinate statuses at every bracket endpoint, closed-form
#'   continuation between them, the record placed at the knot.
#'   \code{FALSE} selects the LEGACY probe bisection with the
#'   support-nesting stop (the pre-2026-08 behaviour, kept as the A/B
#'   baseline). Ignored by the non-CD solvers. Deprecated: it will be removed
#'   with the coordinate-descent family in TRexSelector 2.1.0.
#' @param dummy_distribution Dummy variable distribution, created with
#'   \code{\link{trex_dummy_distribution}} (default: standard Normal).
#' @param scaling_mode Working-design convention, one of \code{"L2"}
#'   (default), \code{"ZSCORE"} or \code{"NONE"}. \code{"L2"} centers
#'   \code{X} and \code{y} inside the selector and divides each column of
#'   \code{X} by its L2 norm; \code{"ZSCORE"} centers and divides by the
#'   sample SD (glmnet \code{standardize = TRUE} convention); the dummies
#'   follow the same convention. \code{"NONE"} takes \code{X} and \code{y}
#'   exactly as supplied (no centering, no scaling): pre-center and
#'   pre-scale as the model demands (\code{\link{trex_normalize_design}}
#'   reproduces the selector's own preprocessing), or supply an
#'   intercept-free dictionary (steering vectors, Vandermonde, Fourier
#'   atoms) whose columns must not be centered. The dummies keep the
#'   unit-L2 standard, so a \code{"NONE"} design should be
#'   unit-L2-commensurate (a warning names columns whose L2 norm leaves
#'   \[0.25, 4\]). Non-finite entries and dead columns (constant under
#'   \code{"L2"} / \code{"ZSCORE"}, all-zero under \code{"NONE"}) are
#'   rejected in every mode with the column named. \code{"L2"} and
#'   \code{"ZSCORE"} standardize the design in place and restore it, so a
#'   memory-mapped design must be open in \code{"readwrite"} mode; under
#'   \code{"NONE"} nothing is written and \code{"readonly"} mode suffices.
#'
#' @return A named list suitable for passing to TRex selector constructors.
#'
#' @examples
#' trex_control()
#' trex_control(solver = "TLASSO", K = 10)
#' trex_control(scaling_mode = "NONE")  # user-managed centering and scaling
#'
#' @family T-Rex selectors and controls
#' @export
trex_control <- function(solver = "TLARS",
                         K = 20,
                         tloop_stagnation_stop = NULL,
                         tloop_max_stagnant_steps = 5,
                         parallel_rnd_experiments = FALSE,
                         use_memory_mapping = FALSE,
                         max_outer_threads = 1,
                         max_inner_threads = 1,
                         max_dummy_multiplier = 10,
                         use_max_T_stop = TRUE,
                         opt_threshold = 0.75,
                         lloop_strategy = "STANDARD",
                         tol = 1e-6,
                         lambda2 = 0.1,
                         rho_afs = 0.3,
                         ncgmp_variant = 0,
                         exch_tie_alpha = 0,
                         exch_tie_floor = 0.5,
                         tenet_aug_use_lars = FALSE,
                         cd_lambda_rel_tol = 1e-6,
                         cd_tol_certify = -1,
                         cd_tol_probe = -1,
                         cd_gram_cap = 0,
                         cd_max_sweeps = 0,
                         cd_certified_localization = TRUE,
                         dummy_distribution = trex_dummy_distribution(),
                         scaling_mode = c("L2", "ZSCORE", "NONE")) {

  if (max_dummy_multiplier < 1) {
    stop("max_dummy_multiplier must be >= 1. Got: ", max_dummy_multiplier)
  }
  if (!is.character(scaling_mode) || length(scaling_mode) < 1 ||
      !(scaling_mode[1] %in% c("L2", "ZSCORE", "NONE"))) {
    stop("scaling_mode must be one of \"L2\", \"ZSCORE\", \"NONE\". Got: ",
         paste(scaling_mode, collapse = ", "))
  }
  scaling_mode <- scaling_mode[1]
  if (tloop_max_stagnant_steps < 0) {
    stop("tloop_max_stagnant_steps must be >= 0. Got: ", tloop_max_stagnant_steps)
  }
  if (max_outer_threads < 1) {
    stop("max_outer_threads must be >= 1. Got: ", max_outer_threads)
  }
  if (max_inner_threads < 1) {
    stop("max_inner_threads must be >= 1. Got: ", max_inner_threads)
  }
  # The C++ side applies the certify/probe pair only when BOTH are > 0;
  # fail loudly here instead of silently ignoring a half-set pair.
  if (xor(cd_tol_certify > 0, cd_tol_probe > 0)) {
    stop("cd_tol_certify and cd_tol_probe must be set together (both > 0), ",
         "or both left at their defaults. Got: cd_tol_certify = ",
         cd_tol_certify, ", cd_tol_probe = ", cd_tol_probe)
  }
  if (cd_gram_cap < 0) {
    stop("cd_gram_cap must be >= 0. Got: ", cd_gram_cap)
  }
  if (cd_max_sweeps < 0) {
    stop("cd_max_sweeps must be >= 0. Got: ", cd_max_sweeps)
  }

  # NULL stays NULL: the C++ core resolves it by its own solver rule, the one
  # source of the default for the three surfaces
  if (!is.null(tloop_stagnation_stop)) {
    tloop_stagnation_stop <- as.logical(tloop_stagnation_stop)
  }

  list(
    solver                  = solver,
    K                       = K,
    max_dummy_multiplier    = max_dummy_multiplier,
    use_max_T_stop          = use_max_T_stop,
    opt_threshold           = opt_threshold,
    lloop_strategy          = lloop_strategy,
    tloop_stagnation_stop   = tloop_stagnation_stop,
    tloop_max_stagnant_steps = tloop_max_stagnant_steps,
    parallel_rnd_experiments = as.logical(parallel_rnd_experiments),
    use_memory_mapping      = use_memory_mapping,
    max_outer_threads       = max_outer_threads,
    max_inner_threads       = max_inner_threads,
    tol                     = tol,
    lambda2                 = lambda2,
    rho_afs                 = rho_afs,
    ncgmp_variant           = ncgmp_variant,
    exch_tie_alpha          = exch_tie_alpha,
    exch_tie_floor          = exch_tie_floor,
    tenet_aug_use_lars      = tenet_aug_use_lars,
    cd_lambda_rel_tol       = cd_lambda_rel_tol,
    cd_tol_certify          = cd_tol_certify,
    cd_tol_probe            = cd_tol_probe,
    cd_gram_cap             = cd_gram_cap,
    cd_max_sweeps           = cd_max_sweeps,
    cd_certified_localization = as.logical(cd_certified_localization),
    dummy_distribution      = dummy_distribution,
    scaling_mode            = scaling_mode
  )
}


#' @name trex_da_control
#'
#' @title TRex DA Control Parameters
#'
#' @description Build a control list for the dependency-aware correction of TRexDASelector.
#'
#' @param da_method Dependency model (default: "BT"). "BT", the binary tree, is made
#'   for disjoint groups of highly correlated variables. "NN", nearest neighbours, is
#'   made for overlapping groups of highly correlated variables and sweeps the
#'   correlation thresholds of a rho grid of length \code{hc_grid_length}. "AR1" is
#'   made for a first-order autoregressive dependency across the columns. "EQUI" is
#'   made for variables that all share one correlation.
#' @param cor_coef Correlation coefficient for AR1/EQUI. -2.0 triggers auto-estimation
#'   (default: -2.0). The estimate is the mean lag-1 autocorrelation of the rows for AR1
#'   and the mean pairwise correlation of the variables for EQUI. Both methods use the
#'   magnitude of a supplied coefficient. Any other value outside \code{[-1, 1]} gives a
#'   warning and is replaced by the estimate.
#' @param rho_thr_DA Correlation threshold of the AR1 and EQUI methods (default: 0.02).
#'   AR1: the window around a variable reaches
#'   \code{ceiling(log(rho_thr_DA) / log(abs(cor_coef)))} columns to each side, up to the
#'   first lag whose correlation \code{abs(cor_coef)^lag} has fallen to the threshold,
#'   range \code{[0, 1)}. EQUI: the correction is skipped when \code{abs(cor_coef)} does
#'   not exceed the threshold, range \code{[0, 1]}.
#' @param hc_linkage Hierarchical clustering linkage: "Single", "Complete", "Average",
#'   "WPGMA" (default: "Single"). Used by BT and by prior groups. "Ward" is rejected,
#'   because its update is valid for squared Euclidean distances only and the correlation
#'   distance \code{1 - abs(cor)} is not Euclidean.
#' @param hc_grid_length Number of correlation thresholds on the calibration grid of BT,
#'   NN and prior groups. \code{0} (default) takes \code{min(20, p)}.
#' @param prior_groups Optional prior grouping constraints: a list of integer
#'   vectors, each of length \code{p}, giving non-negative group labels per
#'   variable. When supplied (non-\code{NULL}, non-empty), the prior-groups
#'   dependency-aware path is used and \code{da_method} is ignored: the finest
#'   common refinement of all supplied levels acts as a hard merge constraint,
#'   hierarchical clustering (\code{hc_linkage}, correlation distance) runs
#'   within each constraint group, and the calibration rho grid (length
#'   \code{hc_grid_length}) is built from the pooled within-group dendrogram
#'   heights plus the conservative rho = 1 singleton anchor. The deflation
#'   therefore acts on tight data-driven neighbourhoods inside the known
#'   groups, never on the raw groups themselves. (default: \code{NULL}).
#' @param bt_selection_mode Which cells of the calibration grid may win the selection
#'   of the methods BT, NN and prior groups. \code{"FeasibleOnly"} (default) admits
#'   only cells whose FDP estimate is at or below the target, which controls the
#'   realized FDR. \code{"RFaithful"} admits every cell with the largest selection
#'   count, as the 1.x versions of the package do, and can exceed the target at
#'   high within-group correlation. Ignored by AR1 and EQUI.
#'
#' @return A named list for use in TRexDASelector.
#'
#' @examples
#' trex_da_control()
#' trex_da_control(da_method = "BT", hc_linkage = "Average")
#' # Prior-groups path (constraints over p = 4 variables):
#' trex_da_control(prior_groups = list(c(1, 1, 2, 2)))
#'
#' @family T-Rex selectors and controls
#' @export
trex_da_control <- function(da_method = "BT",
                            cor_coef = -2.0,
                            rho_thr_DA = 0.02,
                            hc_linkage = "Single",
                            hc_grid_length = 0,
                            prior_groups = NULL,
                            bt_selection_mode = "FeasibleOnly") {
  bt_selection_mode <- match.arg(bt_selection_mode, c("FeasibleOnly", "RFaithful"))
  list(
    da_method         = da_method,
    cor_coef          = cor_coef,
    rho_thr_DA        = rho_thr_DA,
    hc_linkage        = hc_linkage,
    hc_grid_length    = hc_grid_length,
    prior_groups      = prior_groups,
    bt_selection_mode = bt_selection_mode
  )
}


#' @name trex_gvs_control
#'
#' @title TRex GVS Control Parameters
#'
#' @description Build a control list for group variable selection augmentation policy in
#'   TRexGVSSelector. The solver is determined automatically from \code{gvs_type}:
#'   "EN" selects the T-Elastic-Net solver via \code{en_solver} — TENET
#'   (default), with TENET_AUG and TCDENET as optional variants; "IEN" honors
#'   an IEN-family \code{solver} field in \code{\link{trex_control}} —
#'   TIENET (default, native pathwise), with TIENET_AUG and TCDIENET as
#'   optional variants. For \code{gvs_type = "EN"} the \code{solver} field
#'   is ignored. TCDIENET is not finally developed: its sweeps grow with
#'   \code{lambda_2} and on saturating designs, and a solve that does not
#'   converge within \code{cd_max_sweeps} raises an error. TCDENET and
#'   TCDIENET are deprecated and will be removed with the coordinate-descent
#'   family in TRexSelector 2.1.0: the selector raises a deprecation warning
#'   for them. Use TENET and TIENET.
#'
#' @param gvs_type Augmentation policy: "EN" (Elastic Net) or "IEN" (Informed Elastic Net)
#'   (default: "EN").
#' @param corr_max Maximum pairwise correlation for automatic cluster formation (default: 0.5).
#' @param hc_linkage Hierarchical clustering linkage: "Single", "Complete", "Average",
#'   "WPGMA" (default: "Single").
#' @param lambda_2 L2 penalty weight for group structure, in LARS units.
#'   A negative value (default \code{-1}) triggers automatic determination via
#'   the method selected by \code{lambda2_method}; \code{0} is the degenerate
#'   pure T-LASSO case (no ridge); a positive value is used as-is.
#' @param lambda2_method Auto-determination strategy for \code{lambda_2} when
#'   \code{lambda_2 < 0}. \code{NULL} (default) resolves by track:
#'   \code{"CV_1SE_CCD"} for \code{gvs_type = "EN"} and
#'   \code{"CV_1SE_IEN_CCD"} for \code{gvs_type = "IEN"} (the IEN-geometry
#'   tuner consumes the group structure; the 1SE rule is preferred because
#'   the MIN rule tends to pin to the grid floor). The C++ core applies the
#'   same rule at construction on every surface, so an explicit
#'   \code{"CV_1SE_CCD"} on the IEN track also runs the IEN tuner. Explicit choices:
#'   EN-shaped ridge CV \code{"CV_1SE_CCD"}, \code{"CV_MIN_CCD"},
#'   \code{"CV_1SE_SVD"}, \code{"CV_MIN_SVD"} (allowed on both tracks);
#'   IEN-geometry \code{"CV_1SE_IEN_CCD"}, \code{"CV_MIN_IEN_CCD"},
#'   \code{"CV_1SE_TIK_SVD"}, \code{"CV_MIN_TIK_SVD"} (require
#'   \code{gvs_type = "IEN"} — they consume the group structure).
#' @param en_solver Elastic Net solver variant for \code{gvs_type = "EN"}:
#'   \code{"TENET"} (default, Gram-based), \code{"TENET_AUG"}
#'   (row-augmented LASSO), or \code{"TCDENET"} (coordinate descent; exact
#'   fixed-lambda1 elastic-net minimizers per dummy crossing; deprecated,
#'   it will be removed in TRexSelector 2.1.0, use \code{"TENET"}). All solve
#'   the same EN
#'   problem and select the same variables up to path-vs-minimizer
#'   differences; the switch exists to validate solver equivalence,
#'   \strong{not} to choose a LARS vs LASSO path (for that, see
#'   \code{tenet_aug_use_lars}).
#' @param tenet_aug_use_lars Logical. When \code{en_solver = "TENET_AUG"}, run
#'   the augmented system with a pure-LARS inner path that never drops variables
#'   (\code{TRUE}) instead of the default LASSO inner path that can
#'   (\code{FALSE}, default). This is the LARS-vs-LASSO distinction. Ignored for
#'   \code{en_solver = "TENET"}.
#' @param groups Integer vector of manually specified group assignments (1-based, default: NULL).
#' @param group_labels Optional character vector with a name per group of \code{groups}, in
#'   the order of the group numbers. Its length must equal the number of groups. The names
#'   are returned in the field \code{group_labels} of the selector (default: NULL).
#' @param cv_n_folds Number of folds for the k-fold cross-validation used by the
#'   CV \code{lambda2_method} paths (default: 10, matching \code{cv.glmnet}).
#'   Ignored when \code{lambda_2 >= 0}.
#' @param cv_n_lambda Number of points on the geometric lambda grid searched by
#'   the CV auto-computation of \code{lambda_2} (default: 1000). Ignored when
#'   \code{lambda_2 >= 0}.
#' @param cv_n_lambda1 IEN-geometry coordinate-descent tuners
#'   (\code{"CV_1SE_IEN_CCD"}, \code{"CV_MIN_IEN_CCD"}) only: points on the
#'   profiled lambda_1 path per lambda_2 value (default: 40, the shortest
#'   path that still reproduces the 100-point picks; the cost is U-shaped
#'   in this number). Ignored by the other rules and when
#'   \code{lambda_2 >= 0}.
#' @param cv_n_lambda2 IEN-geometry coordinate-descent tuners only: points
#'   on the lambda_2 grid (default: 100). Ignored by the other rules and
#'   when \code{lambda_2 >= 0}.
#' @param cv_dfmax_ratio IEN-geometry coordinate-descent tuners only: support
#'   cap as a fraction of the smallest fold's training size; each lambda_1
#'   path stops once more than \code{floor(cv_dfmax_ratio * n_train)}
#'   coefficients are nonzero (default: 0.8). A value \code{<= 0} disables
#'   the cap. Larger values allow larger supports at a higher tuner cost.
#'   Must be finite. Ignored by the other rules and when
#'   \code{lambda_2 >= 0}.
#' @param cv_seed Seed for the deterministic fold-permutation RNG used by the CV
#'   \code{lambda2_method} paths. \code{-1} (default) derives it from the T-Rex
#'   \code{seed} (deterministic when \code{seed >= 0}, hardware entropy when
#'   \code{seed < 0}); a value \code{>= 0} overrides \code{seed} for CV folds
#'   only, giving per-trial fold reproducibility.
#'
#' @return A named list for use in TRexGVSSelector.
#'
#' @examples
#' trex_gvs_control()
#' trex_gvs_control(gvs_type = "IEN", corr_max = 0.3)
#'
#' @family T-Rex selectors and controls
#' @export
trex_gvs_control <- function(gvs_type = "EN",
                             corr_max = 0.5,
                             hc_linkage = "Single",
                             lambda_2 = -1.0,
                             lambda2_method = NULL,
                             en_solver = "TENET",
                             tenet_aug_use_lars = FALSE,
                             groups = NULL,
                             cv_n_folds = 10,
                             cv_n_lambda = 1000,
                             cv_n_lambda1 = 40,
                             cv_n_lambda2 = 100,
                             cv_dfmax_ratio = 0.8,
                             cv_seed = -1,
                             group_labels = NULL) {
  # NULL resolves by track: the EN-shaped ridge CV for EN, the IEN-geometry
  # tuner (group-consuming, no p/2 conversion) for IEN. 1SE in both cases —
  # the MIN rule tends to pin to the grid floor.
  if (is.null(lambda2_method)) {
    lambda2_method <- if (gvs_type == "IEN") "CV_1SE_IEN_CCD" else "CV_1SE_CCD"
  }
  lambda2_method <- match.arg(lambda2_method,
                              c("CV_1SE_CCD", "CV_MIN_CCD",
                                "CV_1SE_SVD", "CV_MIN_SVD",
                                "CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD",
                                "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD"))
  ien_methods <- c("CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD",
                   "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD")
  if (gvs_type == "EN" && lambda2_method %in% ien_methods) {
    stop("lambda2_method = \"", lambda2_method, "\" is an IEN-geometry ",
         "tuner (it consumes the group structure) and requires ",
         "gvs_type = \"IEN\".")
  }
  en_solver <- match.arg(en_solver, c("TENET", "TENET_AUG", "TCDENET"))
  dfmax_ok <- is.numeric(cv_dfmax_ratio) && length(cv_dfmax_ratio) == 1 &&
    is.finite(cv_dfmax_ratio)
  if (!dfmax_ok) {
    stop("cv_dfmax_ratio must be a single finite number (<= 0 disables the ",
         "support cap). Got: ", format(cv_dfmax_ratio))
  }
  list(
    gvs_type           = gvs_type,
    corr_max           = corr_max,
    hc_linkage         = hc_linkage,
    lambda_2           = lambda_2,
    lambda2_method     = lambda2_method,
    en_solver          = en_solver,
    tenet_aug_use_lars = as.logical(tenet_aug_use_lars),
    groups             = groups,
    group_labels       = if (is.null(group_labels)) NULL else as.character(group_labels),
    cv_n_folds         = as.integer(cv_n_folds),
    cv_n_lambda        = as.integer(cv_n_lambda),
    cv_n_lambda1       = as.integer(cv_n_lambda1),
    cv_n_lambda2       = as.integer(cv_n_lambda2),
    cv_dfmax_ratio     = as.numeric(cv_dfmax_ratio),
    cv_seed            = as.integer(cv_seed)
  )
}


#' @name trex_screen_control
#'
#' @title TRex Screen Control Parameters
#'
#' @description Build a control list for the screening step in TRexScreeningSelector and
#'   TRexBiobankScreeningSelector.
#'
#' @param use_bootstrap_CI Whether to use bootstrap confidence intervals for screening
#'   (default: FALSE).
#' @param R_boot Number of bootstrap resamples, at least 2 (default: 1000).
#' @param ci_grid_step Step of the grid over the confidence level, in \code{(0, 1]}
#'   (default: 0.001).
#' @param trex_method Screening method: "TREX", "TREX_DA_AR1", "TREX_DA_EQUI",
#'   "TREX_DA_BLOCK_EQUI", "TREX_DA_NN" (default: "TREX"). "TREX_DA_NN" compares a variable
#'   with the neighbors inside the AR(1) window whose sample correlation reaches
#'   \code{rho_thr_DA}, and a variable without such a neighbor keeps its relative occurrence.
#' @param cor_coef Correlation coefficient of the dependency-aware variants. -2.0 triggers
#'   auto-estimation (default: -2.0). The variants use the magnitude of a supplied coefficient.
#'   Any other value outside \code{[-1, 1]} gives a warning and is replaced by the estimate.
#' @param rho_thr_DA Correlation threshold of the dependency-aware variants (default: 0.02).
#'   "TREX_DA_AR1": the window around a variable holds the lags whose correlation
#'   \code{abs(cor_coef)^lag} reaches the threshold, range \code{[0, 1)}. "TREX_DA_NN": the
#'   same window, of which the variables stay whose sample correlation with the variable
#'   reaches the threshold, range \code{[0, 1)}. "TREX_DA_EQUI" and
#'   "TREX_DA_BLOCK_EQUI": the correction is skipped when \code{abs(cor_coef)} does not exceed
#'   the threshold, range \code{[0, 1]}.
#' @param n_blocks Number of blocks of "TREX_DA_BLOCK_EQUI", between 1 and the number of
#'   variables (default: 5).
#'
#' @return A named list for use in TRexScreeningSelector or TRexBiobankScreeningSelector.
#'
#' @examples
#' trex_screen_control()
#' trex_screen_control(use_bootstrap_CI = TRUE, R_boot = 500)
#'
#' @family T-Rex selectors and controls
#' @export
trex_screen_control <- function(use_bootstrap_CI = FALSE,
                                R_boot = 1000,
                                ci_grid_step = 0.001,
                                trex_method = "TREX",
                                cor_coef = -2.0,
                                rho_thr_DA = 0.02,
                                n_blocks = 5) {
  list(
    use_bootstrap_CI = use_bootstrap_CI,
    R_boot           = R_boot,
    ci_grid_step     = ci_grid_step,
    trex_method      = trex_method,
    cor_coef         = cor_coef,
    rho_thr_DA       = rho_thr_DA,
    n_blocks         = n_blocks
  )
}


#' @name trex_biobank_control
#'
#' @title TRex Biobank Control Parameters
#'
#' @description Build a control list for the FDR acceptance bounds in
#'   TRexBiobankScreeningSelector. The target FDR level (\code{tFDR}) is a separate
#'   top-level parameter of the selector constructor.
#'
#' @param lower_bound_FDR Lower bound of the acceptable estimated FDR range (default: 0.05).
#' @param upper_bound_FDR Upper bound of the acceptable estimated FDR range (default: 0.15).
#'
#' @return A named list for use in TRexBiobankScreeningSelector.
#'
#' @examples
#' trex_biobank_control()
#' trex_biobank_control(lower_bound_FDR = 0.02, upper_bound_FDR = 0.10)
#'
#' @family T-Rex selectors and controls
#' @export
trex_biobank_control <- function(lower_bound_FDR = 0.05,
                                 upper_bound_FDR = 0.15) {
  list(
    lower_bound_FDR = lower_bound_FDR,
    upper_bound_FDR = upper_bound_FDR
  )
}


#' @name trex_dummy_distribution
#'
#' @title TRex Dummy Variable Distribution
#'
#' @description Build a distribution specification for the dummy variables
#'   injected by the T-Rex Selector. Pass the result to the
#'   \code{dummy_distribution} argument of \code{\link{trex_control}}.
#'
#' @param type Character string naming the distribution.  One of:
#'   \code{"Normal"}, \code{"Uniform"}, \code{"Rademacher"},
#'   \code{"StudentT"}, \code{"Laplace"}, \code{"Gumbel"},
#'   \code{"Triangle"}, \code{"UniformSphere"},
#'   \code{"Mammen"}, \code{"ConstrainedSparseRademacher"}, \code{"Logistic"}.
#'   Defaults to \code{"Normal"}.
#' @param a   Scale / upper-bound parameter (used by \code{"Uniform"}).
#' @param df  Degrees of freedom (used by \code{"StudentT"}).
#' @param location  Location parameter (used by \code{"Laplace"},
#'   \code{"Gumbel"}, \code{"Logistic"}).
#' @param scale  Scale parameter (used by \code{"Laplace"},
#'   \code{"Gumbel"}, \code{"Logistic"}).
#' @param b  Upper vertex of the triangle (used by \code{"Triangle"}).
#' @param c  Right boundary of the triangle (used by \code{"Triangle"}).
#' @param dim  Dimension of the sphere (used by \code{"UniformSphere"}).
#' @param s  Sparsity parameter (used by \code{"ConstrainedSparseRademacher"}).
#'
#' @return A named list describing the dummy distribution, for use in
#'   \code{\link{trex_control}}.
#'
#' @references
#' J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
#'   High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
#'   Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
#'
#' E. Mammen, "Bootstrap and Wild Bootstrap for High Dimensional Linear Models," Ann. Statist., vol. 21,
#'   no. 1, pp. 255-285, 1993, doi: 10.1214/aos/1176349025.
#'
#' @examples
#' # Standard Normal (default)
#' trex_dummy_distribution()
#'
#' # Student-t with 3 degrees of freedom
#' trex_dummy_distribution("StudentT", df = 3)
#'
#' @family T-Rex selectors and controls
#' @export
trex_dummy_distribution <- function(type     = "Normal",
                                    a        = NULL,
                                    df       = NULL,
                                    location = NULL,
                                    scale    = NULL,
                                    b        = NULL,
                                    c        = NULL,
                                    dim      = NULL,
                                    s        = NULL) {
  valid_types <- c("Normal", "Uniform", "Rademacher", "StudentT",
                   "Laplace", "Gumbel", "Triangle",
                   "UniformSphere", "Mammen", "ConstrainedSparseRademacher",
                   "Logistic")
  if (!type %in% valid_types) {
    stop("Unknown distribution type: '", type, "'. ",
         "Must be one of: ", paste(valid_types, collapse = ", "), ".")
  }
  dist <- list(type = type)
  if (!is.null(a))        dist$a        <- a
  if (!is.null(df))       dist$df       <- df
  if (!is.null(location)) dist$location <- location
  if (!is.null(scale))    dist$scale    <- scale
  if (!is.null(b))        dist$b        <- b
  if (!is.null(c))        dist$c        <- c
  if (!is.null(dim))      dist$dim      <- dim
  if (!is.null(s))        dist$s        <- s
  dist
}


#' @name trex_tikhonov_control
#'
#' @title TRex Tikhonov Control Parameters
#'
#' @description Build a control list for \code{\link{TRexTikhonovSelector}},
#'   the general-Tikhonov T-Rex: a penalty matrix
#'   \code{K = t(Gamma) \%*\% Gamma} drives the Tikhonov-net solvers
#'   (TIKNET by default; TCDTIKNET through the \code{solver} field
#'   of \code{\link{trex_control}}, deprecated and removed in TRexSelector
#'   2.1.0). \code{K = I} gives the ridge penalty of the
#'   elastic net, a group-mean K reproduces the GVS-IEN penalty, and a
#'   banded or Laplacian K encodes smoothness. With
#'   \code{penalty_structure = "GROUP_CENTERING"} the selector builds K
#'   itself from \code{prior_groups} (block-diagonal \code{I - J / p_m}:
#'   the members of a group share a coefficient).
#'
#' @param tikhonov_K Penalty matrix K (p x p): a base \code{matrix} or a
#'   sparse \code{Matrix::dgCMatrix}. Build it from a regularization
#'   operator Gamma with \code{\link{trex_gamma_to_k}}. Required for
#'   \code{penalty_structure = "SPARSE_GENERAL"}. Optional under
#'   \code{"GROUP_CENTERING"}: left out, the selector builds the centering
#'   matrix from \code{prior_groups}. A supplied matrix is used as it is
#'   and has to be that centering matrix.
#' @param penalty_structure \code{"SPARSE_GENERAL"} (default): the
#'   user-supplied K. \code{"GROUP_CENTERING"}: K is the block-diagonal
#'   centering penalty over the groups of \code{prior_groups}, which this
#'   structure requires. An all-singleton partition is legitimate there: K
#'   vanishes and the path is the plain lasso on per-variable cluster
#'   dummies. With \code{solver = "TCDTIKNET"} this path runs through the
#'   group-centering hooks of TCDIENET, so its warnings (sweep budget,
#'   certification fallback) carry the tag \code{[TCDIENET]}.
#' @param lambda_2 Penalty weight, solver scale. A negative value (default
#'   \code{-1}) triggers auto-determination via \code{lambda2_method};
#'   \code{0} is the degenerate pure-lasso case; a positive value is used
#'   as-is.
#' @param lambda2_method Auto-determination strategy when
#'   \code{lambda_2 < 0}. \code{"PILOT_KNEE"} (default) reads lambda_2 off
#'   a few random experiments per candidate of \code{pilot_grid}, by the
#'   selector's own principle: the smallest lambda_2 whose pilot selection
#'   count is within \code{pilot_knee_tolerance} of the largest (the knee
#'   of the profile). The pilot record is available as \code{lambda2_pilot}
#'   on the selector. The k-fold CV rules consume K and return lambda_2 on
#'   the solver scale (no p/2 conversion): \code{"CV_1SE_IEN_CCD"}
#'   (well-posed in any p/n regime), \code{"CV_MIN_IEN_CCD"},
#'   \code{"CV_1SE_TIK_SVD"}, \code{"CV_MIN_TIK_SVD"}. The knee of the full
#'   selection profile, one selection per grid value, is
#'   \code{\link{trex_tikhonov_lambda2_power_rule}}.
#' @param cv_n_folds Folds for the CV paths (default: 10).
#' @param cv_n_lambda Points on the CV lambda grid of the SVD rules
#'   (default: 100).
#' @param cv_n_lambda1 IEN-geometry coordinate-descent rules: points on the
#'   profiled lambda_1 path per lambda_2 value (default: 40).
#' @param cv_n_lambda2 IEN-geometry coordinate-descent rules: points on the
#'   lambda_2 grid (default: 100).
#' @param cv_dfmax_ratio IEN-geometry coordinate-descent rules: support cap
#'   as a fraction of the smallest fold's training size; each lambda_1 path
#'   stops once more than \code{floor(cv_dfmax_ratio * n_train)}
#'   coefficients are nonzero (default: 0.8). A value \code{<= 0} disables
#'   the cap. Larger values allow larger supports at a higher tuner cost.
#'   Must be finite.
#' @param cv_seed Seed for the CV fold permutation; \code{-1} (default)
#'   derives it from the selector seed, \code{>= 0} overrides it for CV
#'   folds only.
#' @param pilot_grid Pilot rule: lambda_2 candidates on the working column
#'   scale (default \code{c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)}),
#'   sorted and de-duplicated internally. Keep 0 on it so that "no
#'   penalty" is an admissible answer.
#' @param pilot_K Pilot rule: pilot experiments per grid point (default:
#'   10; at least 2, capped at \code{K} of \code{\link{trex_control}}).
#' @param pilot_T Pilot rule: stopping time of the pilot paths (default:
#'   5; at least 1, capped at \code{min(L_0, ceiling(n / 2))}).
#' @param pilot_L_multiplier Pilot rule: pilot dummies
#'   \code{L_0 = pilot_L_multiplier * p} (default: 1; at most
#'   \code{max_dummy_multiplier}).
#' @param pilot_knee_tolerance Pilot rule: selections below the largest
#'   pilot count that still count as "on the plateau" (default: 1; 0
#'   reads the largest count, ties toward the smallest lambda_2).
#' @param pilot_statistic Pilot rule, the count the pick is read off:
#'   \code{"MAJORITY"} (default; originals ahead of the T-th dummy in more
#'   than half of the pilot experiments), \code{"CALIBRATED"} (the T-Rex
#'   selection rule applied to the pilot data), or \code{"GROUP_MAJORITY"}
#'   (the majority count restricted to the groups of \code{prior_groups}
#'   in which at least \code{pilot_group_fraction} of the members are
#'   ahead; needs \code{prior_groups}).
#' @param pilot_group_fraction \code{"GROUP_MAJORITY"}: the share of a
#'   group's members that must be ahead, in (0, 1] (default: 0.5).
#' @param pilot_shared_dummies Logical (default \code{FALSE}): the pilot
#'   draws its dummies from a reserved seed stream, independent of the
#'   main run's. \code{TRUE} reuses the main strategy's first dummy block
#'   (not available under the PERMUTATION strategies).
#' @param fold_dummy_coupling Logical (default \code{TRUE}): dummy copies
#'   join their originating variable's K-coupling (FOLDED — the exact
#'   generalization of the GVS layered convention; dummy/null
#'   exchangeability holds for arbitrary K). \code{FALSE} restores the
#'   decoupled unit-ridge dummy block (INDEPENDENT_RIDGE) — a validation
#'   escape hatch that preserves the \code{K = I} == TCDENET collapse but
#'   deforms exchangeability whenever K has off-diagonals. Only the
#'   TCDTIKNET solver offers it, so it needs \code{solver = "TCDTIKNET"} in
#'   \code{\link{trex_control}}; the default solver TIKNET raises an error.
#'   Deprecated: it will be removed with TCDTIKNET and the coordinate-descent
#'   family in TRexSelector 2.1.0.
#' @param cluster_dummies Logical (default \code{TRUE}): draw dummies as
#'   per-cluster colored MVN layers (GVS parity) on the clustering below.
#'   \code{FALSE} restores canonical i.i.d. dummies — validation runs only;
#'   with an informed K, i.i.d. dummies deform the dummy/null
#'   exchangeability the FDR calibration relies on.
#' @param prior_groups Integer vector of 1-based cluster IDs per variable
#'   (length p) for the dummy clustering. \code{NULL} (default): clusters
#'   are discovered by correlation-HAC with \code{corr_max} /
#'   \code{hc_linkage}. Ignored when \code{cluster_dummies = FALSE}.
#' @param corr_max Maximum cross-cluster |correlation| for the HAC dummy
#'   clustering (default: 0.5). Ignored when \code{prior_groups} is set or
#'   \code{cluster_dummies = FALSE}.
#' @param hc_linkage HAC linkage: "Single" (default), "Complete",
#'   "Average", "WPGMA". Ignored under \code{prior_groups} or
#'   \code{cluster_dummies = FALSE}.
#'
#' @return A named list for use in TRexTikhonovSelector.
#'
#' @examples
#' trex_tikhonov_control(tikhonov_K = diag(10))
#' trex_tikhonov_control(penalty_structure = "GROUP_CENTERING",
#'                       prior_groups = rep(1:2, each = 5),
#'                       pilot_grid = c(0, 1, 10, 100))
#' trex_tikhonov_control(tikhonov_K = diag(10),
#'                       lambda2_method = "CV_1SE_IEN_CCD")
#'
#' @family T-Rex selectors and controls
#' @export
trex_tikhonov_control <- function(tikhonov_K = NULL,
                                  penalty_structure = "SPARSE_GENERAL",
                                  lambda_2 = -1.0,
                                  lambda2_method = "PILOT_KNEE",
                                  cv_n_folds = 10,
                                  cv_n_lambda = 100,
                                  cv_n_lambda1 = 40,
                                  cv_n_lambda2 = 100,
                                  cv_dfmax_ratio = 0.8,
                                  cv_seed = -1,
                                  pilot_grid = c(0, 1, 3, 10, 30, 100,
                                                 300, 1000, 3000),
                                  pilot_K = 10,
                                  pilot_T = 5,
                                  pilot_L_multiplier = 1,
                                  pilot_knee_tolerance = 1,
                                  pilot_statistic = "MAJORITY",
                                  pilot_group_fraction = 0.5,
                                  pilot_shared_dummies = FALSE,
                                  fold_dummy_coupling = TRUE,
                                  cluster_dummies = TRUE,
                                  prior_groups = NULL,
                                  corr_max = 0.5,
                                  hc_linkage = "Single") {
  penalty_structure <- match.arg(penalty_structure,
                                 c("SPARSE_GENERAL", "GROUP_CENTERING"))
  lambda2_method <- match.arg(lambda2_method,
                              c("PILOT_KNEE",
                                "CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD",
                                "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD"))
  pilot_statistic <- match.arg(pilot_statistic,
                               c("MAJORITY", "CALIBRATED", "GROUP_MAJORITY"))
  hc_linkage <- match.arg(hc_linkage,
                          c("Single", "Complete", "Average", "WPGMA"))
  if (penalty_structure == "SPARSE_GENERAL" && is.null(tikhonov_K)) {
    stop("tikhonov_K is required (p x p penalty matrix; use trex_gamma_to_k() ",
         "to build K from a regularization operator Gamma) unless ",
         "penalty_structure = \"GROUP_CENTERING\".")
  }
  dfmax_ok <- is.numeric(cv_dfmax_ratio) && length(cv_dfmax_ratio) == 1 &&
    is.finite(cv_dfmax_ratio)
  if (!dfmax_ok) {
    stop("cv_dfmax_ratio must be a single finite number (<= 0 disables the ",
         "support cap). Got: ", format(cv_dfmax_ratio))
  }
  if (!is.numeric(pilot_grid) || length(pilot_grid) == 0 ||
      any(!is.finite(pilot_grid)) || any(pilot_grid < 0)) {
    stop("pilot_grid must be a non-empty vector of finite values >= 0.")
  }
  if (pilot_K < 2) {
    stop("pilot_K must be >= 2. Got: ", pilot_K)
  }
  if (pilot_T < 1) {
    stop("pilot_T must be >= 1. Got: ", pilot_T)
  }
  if (pilot_L_multiplier < 1) {
    stop("pilot_L_multiplier must be >= 1. Got: ", pilot_L_multiplier)
  }
  if (pilot_knee_tolerance < 0) {
    stop("pilot_knee_tolerance must be >= 0. Got: ", pilot_knee_tolerance)
  }
  if (!(pilot_group_fraction > 0 && pilot_group_fraction <= 1)) {
    stop("pilot_group_fraction must be in (0, 1]. Got: ", pilot_group_fraction)
  }
  list(
    tikhonov_K           = tikhonov_K,
    penalty_structure    = penalty_structure,
    lambda_2             = lambda_2,
    lambda2_method       = lambda2_method,
    cv_n_folds           = as.integer(cv_n_folds),
    cv_n_lambda          = as.integer(cv_n_lambda),
    cv_n_lambda1         = as.integer(cv_n_lambda1),
    cv_n_lambda2         = as.integer(cv_n_lambda2),
    cv_dfmax_ratio       = as.numeric(cv_dfmax_ratio),
    cv_seed              = as.integer(cv_seed),
    pilot_grid           = as.numeric(pilot_grid),
    pilot_K              = as.integer(pilot_K),
    pilot_T              = as.integer(pilot_T),
    pilot_L_multiplier   = as.integer(pilot_L_multiplier),
    pilot_knee_tolerance = as.numeric(pilot_knee_tolerance),
    pilot_statistic      = pilot_statistic,
    pilot_group_fraction = as.numeric(pilot_group_fraction),
    pilot_shared_dummies = as.logical(pilot_shared_dummies),
    fold_dummy_coupling  = fold_dummy_coupling,
    cluster_dummies      = cluster_dummies,
    prior_groups         = prior_groups,
    corr_max             = corr_max,
    hc_linkage           = hc_linkage
  )
}
