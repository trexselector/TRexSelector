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
#'   \code{TRexGVSSelector} with \code{gvs_type = "IEN"}; "TIKNET" needs a
#'   penalty matrix and is dispatchable only inside
#'   \code{TRexTikhonovSelector}. The plain selector rejects both families
#'   with an informative error. The coordinate-descent names "TCDLASSO",
#'   "TCDENET" and "TCDIENET" are deprecated and will be removed in
#'   TRexSelector 2.1.0: a selector that receives one of them raises a
#'   deprecation warning. "TCDTIKNET" is deprecated as well, and no selector
#'   accepts it any more. Use "TLASSO", "TENET", "TIENET" and "TIKNET", the
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


# Group numbers of an R argument: 1 to M without gaps, one whole number per
# variable. The core counts from 0, so its own messages would not fit R.
.check_group_numbers <- function(groups, name) {
  whole <- is.numeric(groups) && length(groups) > 0 && !anyNA(groups) &&
    all(groups == round(groups))
  if (!whole || min(groups) < 1 || length(unique(groups)) != max(groups)) {
    stop(name, " must number the groups 1 to M without gaps, with one whole ",
         "number per variable.")
  }
  invisible(NULL)
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
#' @param hc_linkage Hierarchical clustering linkage: "single", "complete", "average",
#'   "wpgma" (default: "single"). Used by BT and by prior groups. "ward" is rejected,
#'   because its update is valid for squared Euclidean distances only and the correlation
#'   distance \code{1 - abs(cor)} is not Euclidean. The capitalized names of
#'   TRexSelector 2.0.1 ("Single") are accepted as well.
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
#' trex_da_control(da_method = "BT", hc_linkage = "average")
#' # Prior-groups path (constraints over p = 4 variables):
#' trex_da_control(prior_groups = list(c(1, 1, 2, 2)))
#'
#' @family T-Rex selectors and controls
#' @export
trex_da_control <- function(da_method = "BT",
                            cor_coef = -2.0,
                            rho_thr_DA = 0.02,
                            hc_linkage = "single",
                            hc_grid_length = 0,
                            prior_groups = NULL,
                            bt_selection_mode = "FeasibleOnly") {
  bt_selection_mode <- match.arg(bt_selection_mode, c("FeasibleOnly", "RFaithful"))
  list(
    da_method         = da_method,
    cor_coef          = cor_coef,
    rho_thr_DA        = rho_thr_DA,
    hc_linkage        = tolower(hc_linkage),
    hc_grid_length    = hc_grid_length,
    prior_groups      = prior_groups,
    bt_selection_mode = bt_selection_mode
  )
}


#' @name trex_gvs_control
#'
#' @title TRex GVS Control Parameters
#'
#' @description Build a control list for the grouped variable selection of
#'   \code{\link{TRexGVSSelector}}. The solver follows \code{gvs_type}. For
#'   "EN" \code{en_solver} chooses it: TENET (default), TENET_AUG or the
#'   deprecated TCDENET. For "IEN" the \code{solver} of
#'   \code{\link{trex_control}} names it, TIENET (native pathwise),
#'   TIENET_AUG or the deprecated TCDIENET, so the IEN track needs a base
#'   control with one of them. The selector refuses a base control whose
#'   solver does not fit the track. TCDIENET is not finally developed: its
#'   sweeps grow with \code{lambda_2} and on saturating designs, and a solve
#'   that does not converge within \code{cd_max_sweeps} raises an error.
#'   TCDENET and TCDIENET are deprecated and will be removed with the
#'   coordinate-descent family in TRexSelector 2.1.0: the selector raises a
#'   deprecation warning for them. Use TENET and TIENET.
#'
#' @param gvs_type Penalty of the grouped selection: "EN" (elastic net) or "IEN" (informed
#'   elastic net) (default: "EN").
#' @param corr_max Correlation threshold of the automatic clustering (default: 0.5). The
#'   dendrogram is cut at the height \code{1 - corr_max}. With "single" linkage (default)
#'   no two variables of different clusters correlate above \code{corr_max}. With
#'   "complete" linkage every two variables of one cluster correlate at least at
#'   \code{corr_max}, and variables of different clusters can exceed it.
#' @param hc_linkage Hierarchical clustering linkage: "single", "complete", "average",
#'   "wpgma" (default: "single"). Used when no \code{groups} are supplied. "ward" is
#'   rejected, because its update is valid for squared Euclidean distances only and the
#'   correlation distance \code{1 - abs(cor)} is not Euclidean. The capitalized names of
#'   TRexSelector 2.0.1 ("Single") are accepted as well.
#' @param lambda_2 L2 penalty weight for group structure.
#'   A negative value (default \code{-1}) triggers automatic determination via
#'   the method selected by \code{lambda2_method}; \code{0} is the degenerate
#'   pure T-LASSO case (no ridge); a positive value is used as-is. A supplied value
#'   refers to the scale of the standardized columns. Under
#'   \code{scaling_mode = "ZSCORE"} the columns have the squared norm \code{n - 1}, so the
#'   same shrinkage needs the value times \code{n - 1}.
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
#' @param tenet_aug_use_lars Logical. For the augmented solvers, \code{en_solver =
#'   "TENET_AUG"} or \code{gvs_type = "IEN"} with the solver \code{"TIENET_AUG"}: run
#'   the augmented system with a pure-LARS inner path that never drops variables
#'   (\code{TRUE}) instead of the default LASSO inner path that can
#'   (\code{FALSE}, default). This is the LARS-vs-LASSO distinction. Ignored for
#'   \code{"TENET"} and \code{"TIENET"}.
#' @param groups Integer vector of manually specified group assignments, one entry per
#'   variable, numbered 1 to M without gaps (default: NULL).
#' @param group_labels Optional character vector with a name per group of \code{groups}, in
#'   the order of the group numbers. Its length must equal the number of groups. The names
#'   are returned in the field \code{group_labels} of the selector (default: NULL).
#' @param cv_n_folds Number of folds for the k-fold cross-validation used by the
#'   CV \code{lambda2_method} paths (default: 10, matching \code{cv.glmnet}).
#'   Ignored when \code{lambda_2 >= 0}.
#' @param cv_n_lambda Number of points on the \code{lambda_2} grid of the rules
#'   \code{"CV_1SE_SVD"}, \code{"CV_MIN_SVD"}, \code{"CV_1SE_CCD"}, \code{"CV_MIN_CCD"},
#'   \code{"CV_1SE_TIK_SVD"} and \code{"CV_MIN_TIK_SVD"} (default: 1000). The rules
#'   \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"} use \code{cv_n_lambda2}.
#'   Ignored when \code{lambda_2 >= 0}.
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
                             hc_linkage = "single",
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
  if (!is.null(groups)) {
    .check_group_numbers(groups, "groups")
  }
  dfmax_ok <- is.numeric(cv_dfmax_ratio) && length(cv_dfmax_ratio) == 1 &&
    is.finite(cv_dfmax_ratio)
  if (!dfmax_ok) {
    stop("cv_dfmax_ratio must be a single finite number (<= 0 disables the ",
         "support cap). Got: ", format(cv_dfmax_ratio))
  }
  list(
    gvs_type           = gvs_type,
    corr_max           = corr_max,
    hc_linkage         = tolower(hc_linkage),
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
#' @param use_bootstrap_CI Whether the confidence-based screening (bootstrap) makes the
#'   selection (default: FALSE). The FDR estimate stays one over the number of selected
#'   variables.
#' @param R_boot Number of bootstrap resamples, at least 2 (default: 1000).
#' @param ci_grid_step Step of the grid over the confidence level, in \code{(0, 1]}
#'   (default: 0.001).
#' @param trex_method Variant of the screening (default: "TREX"). "TREX" screens without a
#'   dependency-aware filter. The other variants apply a dependency-aware filter to the
#'   relative occurrences. "TREX_DA_AR1", for a first-order autoregressive dependency,
#'   compares a variable with the neighbors inside the AR(1) window. "TREX_DA_EQUI", for
#'   equi-correlated variables, compares a variable with all other variables.
#'   "TREX_DA_BLOCK_EQUI", for block equi-correlated variables, makes the same comparison
#'   inside contiguous blocks. "TREX_DA_NN" compares a variable with the neighbors inside the
#'   AR(1) window whose sample correlation reaches \code{rho_thr_DA}, and a variable without
#'   such a neighbor keeps its relative occurrence.
#' @param cor_coef Correlation coefficient of the dependency-aware variants. -2.0 triggers
#'   auto-estimation (default: -2.0). The variants use the magnitude of a supplied coefficient.
#'   Any other value outside \code{[-1, 1]} gives a warning and is replaced by the estimate.
#' @param rho_thr_DA Correlation threshold of the dependency-aware variants (default: 0.02).
#'   "TREX_DA_AR1": the window around a variable reaches
#'   \code{ceiling(log(rho_thr_DA) / log(abs(cor_coef)))} columns to each side, up to the
#'   first lag whose correlation \code{abs(cor_coef)^lag} has fallen to the threshold, range
#'   \code{[0, 1)}. "TREX_DA_NN": the same window, of which the variables stay whose sample
#'   correlation with the variable reaches the threshold, range \code{[0, 1)}.
#'   "TREX_DA_EQUI" and "TREX_DA_BLOCK_EQUI": the filter is skipped when
#'   \code{abs(cor_coef)} does not exceed the threshold, range \code{[0, 1]}.
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
#' @description Build a control list for TRexBiobankScreeningSelector: the accepted range of
#'   the FDR estimate of the screening and the target FDR level of the fallback T-Rex
#'   selector.
#'
#' @param lower_bound_FDR Lower bound of the accepted range of the FDR estimate of the
#'   screening (default: 0.05). An estimate below it suggests too little power.
#' @param upper_bound_FDR Upper bound of the accepted range of the FDR estimate of the
#'   screening (default: 0.15). An estimate above it suggests too many false discoveries.
#' @param target_FDR_trex Target FDR level of the fallback T-Rex selector (default: 0.1).
#'
#' @return A named list for use in TRexBiobankScreeningSelector.
#'
#' @examples
#' trex_biobank_control()
#' trex_biobank_control(lower_bound_FDR = 0.02, upper_bound_FDR = 0.10)
#' trex_biobank_control(target_FDR_trex = 0.2)
#'
#' @family T-Rex selectors and controls
#' @export
trex_biobank_control <- function(lower_bound_FDR = 0.05,
                                 upper_bound_FDR = 0.15,
                                 target_FDR_trex = 0.1) {
  list(
    lower_bound_FDR = lower_bound_FDR,
    upper_bound_FDR = upper_bound_FDR,
    target_FDR_trex = target_FDR_trex
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
#' @description Build the control list of \code{\link{TRexTikhonovSelector}}, the
#'   T-Rex selector with a general Tikhonov penalty (T-Rex+Tikhonov). The list
#'   carries the penalty, the rule that chooses its weight \code{lambda_2}, and
#'   how the dummies are drawn. The help of the class describes the method.
#'
#' @param tikhonov_K Penalty matrix K (p x p), symmetric and positive
#'   semi-definite: a base \code{matrix} or a sparse \code{Matrix::dgCMatrix}.
#'   \code{\link{trex_gamma_to_k}} builds it from a regularization operator
#'   Gamma. Required for \code{penalty_structure = "SPARSE_GENERAL"}, where the
#'   selector checks its content once at construction (finite entries,
#'   symmetry, positive semi-definiteness). Optional under
#'   \code{"GROUP_CENTERING"}: left out, the selector builds the centering
#'   matrix from \code{prior_groups}. A supplied matrix is compared with that
#'   centering matrix, and another matrix is an error.
#' @param penalty_structure Structural form of the penalty.
#'   \code{"SPARSE_GENERAL"} (default): the supplied \code{tikhonov_K} is the
#'   penalty. \code{"GROUP_CENTERING"}: K is the block-diagonal centering
#'   matrix \code{I - J / p_m} over the groups of \code{prior_groups}, which
#'   this structure requires, so the members of a group share a coefficient
#'   (\code{\link{trex_build_group_centering_k}} returns the matrix). A
#'   partition into single variables is allowed: K is then zero and the path
#'   is the lasso path.
#' @param lambda_2 Weight of the penalty. A negative value (default: -1) lets
#'   \code{lambda2_method} choose it, \code{0} is no penalty (the lasso path of
#'   the solver), and a positive value is used as it is. The value must be
#'   finite. A supplied value refers to the scale of the standardized columns.
#'   Under \code{scaling_mode = "ZSCORE"} the columns have the squared norm
#'   \code{n - 1}, so the same shrinkage needs the value times \code{n - 1}.
#'   The rules report their pick on this scale.
#' @param lambda2_method Rule that chooses \code{lambda_2} when
#'   \code{lambda_2 < 0}. \code{"PILOT_KNEE"} (default), the pilot rule, reads
#'   it off a few random experiments per value of \code{pilot_grid}: the
#'   smallest lambda_2 whose pilot selection count is within
#'   \code{pilot_knee_tolerance} of the largest, the knee of the profile. The
#'   field \code{lambda2_pilot} of the selector holds its record. The
#'   cross-validation rules tune on the normalized data before the random
#'   experiments. \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"} profile
#'   the lambda_1 path out per lambda_2 candidate and are well-posed for every
#'   p and n. \code{"CV_1SE_TIK_SVD"} and \code{"CV_MIN_TIK_SVD"} are an
#'   analytic ridge cross-validation (lambda_1 = 0), exact for n > p and
#'   refused when the unpenalized null space interpolates the training folds.
#'   The knee of the full selection profile, one selection per grid value, is
#'   \code{\link{trex_tikhonov_lambda2_power_rule}}.
#' @param cv_n_folds Folds of the cross-validation rules, 2 to n (default: 10).
#' @param cv_n_lambda \code{"CV_1SE_TIK_SVD"} and \code{"CV_MIN_TIK_SVD"}:
#'   points on the lambda_2 grid (default: 100).
#' @param cv_n_lambda1 \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"}:
#'   points on the profiled lambda_1 path per lambda_2 value (default: 40).
#' @param cv_n_lambda2 \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"}:
#'   points on the lambda_2 grid (default: 100).
#' @param cv_dfmax_ratio \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"}:
#'   support cap as a fraction of the smallest fold's training size. Each
#'   lambda_1 path stops once more than
#'   \code{floor(cv_dfmax_ratio * n_train)} coefficients are nonzero (default:
#'   0.8). A value \code{<= 0} disables the cap. Larger values allow larger
#'   supports at a higher cost. Must be finite.
#' @param cv_seed Seed of the fold permutation of the cross-validation rules.
#'   \code{-1} (default) derives it from the seed of the selector, a value
#'   \code{>= 0} replaces it for the folds only.
#' @param pilot_grid Pilot rule: lambda_2 candidates, stated for columns of
#'   unit norm, sorted and de-duplicated internally. The default is the grid
#'   of the C++ core, \code{\link{trex_default_lambda2_grid}()}: 0, 1, 3, 10,
#'   30, 100, 300, 1000, 3000. Under \code{scaling_mode = "ZSCORE"} the pilot
#'   runs at these values times \code{n - 1}, as the cross-validation rules
#'   return their pick, so one grid means the same shrinkage under every
#'   scaling mode. The pick is reported on the working scale and can be
#'   supplied as \code{lambda_2}. Keep 0 on the grid so that "no penalty" is
#'   an admissible answer: a penalty that misfits the truth loses power as
#'   lambda_2 grows, and the rule then switches it off.
#' @param pilot_K Pilot rule: pilot experiments per grid value (default: 10;
#'   at least 2, capped at \code{K} of \code{\link{trex_control}}). The
#'   majority count has the resolution \code{1 / pilot_K}.
#' @param pilot_T Pilot rule: stopping time of the pilot paths (default: 5;
#'   at least 1, capped at \code{min(L_0, ceiling(n / 2))}, as the T-loop caps
#'   its stopping time). The pick reads the count at this stopping time.
#' @param pilot_L_multiplier Pilot rule: pilot dummies
#'   \code{L_0 = pilot_L_multiplier * p} (default: 1; at most
#'   \code{max_dummy_multiplier} of \code{\link{trex_control}}).
#' @param pilot_knee_tolerance Pilot rule: selections below the largest
#'   pilot count that still count as "on the plateau" (default: 1). \code{0}
#'   reads the largest count, ties toward the smallest lambda_2.
#' @param pilot_statistic Pilot rule: the count the pick is read off, per grid
#'   value. \code{"MAJORITY"} (default) is the number of variables ahead of
#'   the last pilot dummy (the \code{pilot_T}-th) in more than half of the
#'   pilot experiments, the count of the T-Rex selector at the voting level
#'   1/2. It is zero under the global null. \code{"CALIBRATED"} is the T-Rex
#'   selection rule applied to the pilot data, the largest selection whose
#'   estimated FDP stays within the target. It is coarser at a small
#'   \code{pilot_K}. \code{"GROUP_MAJORITY"} is the majority count restricted
#'   to the groups of \code{prior_groups} in which at least
#'   \code{pilot_group_fraction} of the members are ahead. A group that enters
#'   as a whole counts with its members, a group that enters as a minority
#'   counts nothing. It needs \code{prior_groups}. The record
#'   \code{lambda2_pilot} holds all three counts.
#' @param pilot_group_fraction \code{"GROUP_MAJORITY"}: the share of a
#'   group's members that must be ahead for the group to count, in (0, 1]
#'   (default: 0.5).
#' @param pilot_shared_dummies Logical (default \code{FALSE}): the pilot
#'   draws its dummies from a seed stream of its own, so the dummies of the
#'   main run do not depend on the pick, and the main run equals a run at the
#'   picked lambda_2 with the same seed. \code{TRUE}: the pilot uses the first
#'   dummy block of the main strategy (not with the PERMUTATION strategies).
#' @param fold_dummy_coupling Logical (default \code{TRUE}): in the penalty
#'   every dummy takes the place of the variable it was drawn for (folded
#'   coupling), so dummies and null variables stay exchangeable for every K.
#'   \code{FALSE}, a ridge block of the dummies decoupled from K, is an
#'   error: it existed with the solver \code{"TCDTIKNET"} only, which the
#'   selector no longer accepts. Deprecated: the argument will be removed
#'   with the coordinate-descent family in TRexSelector 2.1.0.
#' @param cluster_dummies Logical (default \code{TRUE}): draw the dummies per
#'   cluster of correlated variables, as in \code{\link{TRexGVSSelector}}.
#'   \code{FALSE} gives i.i.d. dummies, for validation runs only: with an
#'   informed K they deform the exchangeability of dummies and null variables
#'   that the FDR calibration relies on.
#' @param prior_groups Groups of the variables: an integer vector with the
#'   group of each variable, numbered 1 to M without gaps (length p). They
#'   define the clusters of the dummies (with \code{cluster_dummies}), the
#'   penalty under \code{"GROUP_CENTERING"} and the count
#'   \code{"GROUP_MAJORITY"} of the pilot rule. \code{NULL} (default): the
#'   clusters of the dummies are found by clustering the correlations with
#'   \code{corr_max} and \code{hc_linkage}.
#' @param corr_max Correlation threshold of the automatic clustering, between
#'   0 and 1 (default: 0.5). The dendrogram is cut at the height
#'   \code{1 - corr_max}. With "single" linkage (default) no two variables of
#'   different clusters correlate above \code{corr_max}. With "complete"
#'   linkage every two variables of one cluster correlate at least at
#'   \code{corr_max}, and variables of different clusters can exceed it.
#'   Ignored when \code{prior_groups} is set or \code{cluster_dummies = FALSE}.
#' @param hc_linkage Hierarchical clustering linkage: "single", "complete",
#'   "average", "wpgma" (default: "single"). "ward" is rejected, because its
#'   update is valid for squared Euclidean distances only and the correlation
#'   distance \code{1 - abs(cor)} is not Euclidean. Ignored when
#'   \code{prior_groups} is set or \code{cluster_dummies = FALSE}. The
#'   capitalized names of TRexSelector 2.0.1 ("Single") are accepted as well.
#'
#' @return A named list for use in \code{\link{TRexTikhonovSelector}}.
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
                                  pilot_grid = trex_default_lambda2_grid(),
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
                                  hc_linkage = "single") {
  penalty_structure <- match.arg(penalty_structure,
                                 c("SPARSE_GENERAL", "GROUP_CENTERING"))
  lambda2_method <- match.arg(lambda2_method,
                              c("PILOT_KNEE",
                                "CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD",
                                "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD"))
  pilot_statistic <- match.arg(pilot_statistic,
                               c("MAJORITY", "CALIBRATED", "GROUP_MAJORITY"))
  hc_linkage <- match.arg(tolower(hc_linkage),
                          c("single", "complete", "average", "wpgma"))
  if (!is.null(prior_groups)) {
    .check_group_numbers(prior_groups, "prior_groups")
  }
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
