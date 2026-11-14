# The solver of a base control that is passed to the Tikhonov selector. NULL
# stands for the base control of the C++ class, whose solver is TIKNET.
.check_tikhonov_solver <- function(control) {
  if (is.null(control)) {
    return(invisible(NULL))
  }
  solver <- control$solver
  if (!identical(solver, "TIKNET")) {
    stop("control$solver must be \"TIKNET\" for the Tikhonov selector. Got: \"",
         paste(format(solver), collapse = ", "), "\". ",
         if (identical(solver, "TCDTIKNET")) {
           "The deprecated solver \"TCDTIKNET\" is no longer accepted by the selector."
         } else {
           paste0("trex_control() names \"TLARS\" by default, so pass ",
                  "trex_control(solver = \"TIKNET\", ...), or leave control out ",
                  "to run the default control.")
         })
  }
  invisible(NULL)
}

# One group number per variable, checked with the name of the R argument.
.check_tikhonov_groups_length <- function(tik_control, X) {
  groups <- tik_control$prior_groups
  if (!is.null(groups) && length(groups) != ncol(X)) {
    stop("tik_control$prior_groups must have one entry per variable (p = ", ncol(X),
         "). Got: ", length(groups), ".")
  }
  invisible(NULL)
}


#' @name TRexTikhonovSelector
#'
#' @title T-Rex Selector with a General Tikhonov Penalty (T-Rex+Tikhonov)
#'
#' @description FDR-controlled variable selection with the penalty
#'   \code{(lambda_2 / 2) * t(beta) \%*\% K \%*\% beta}, whose positive
#'   semi-definite matrix K carries the structure among the variables.
#'
#' @section The penalty:
#' \code{K = t(Gamma) \%*\% Gamma} (p x p) for a regularization operator Gamma,
#' see \code{\link{trex_gamma_to_k}}. Special cases:
#' \itemize{
#'   \item \code{K = I}: the ridge part of the elastic net.
#'   \item The group-mean matrix, with the block \code{J / p_m} for every
#'     group m of size \code{p_m} (J is the matrix of ones): the penalty of the
#'     informed elastic net of \code{\link{TRexGVSSelector}}.
#'   \item The centering matrix, with the block \code{I - J / p_m} for every
#'     group (\code{penalty_structure = "GROUP_CENTERING"},
#'     \code{\link{trex_build_group_centering_k}}): the members of a group
#'     share a coefficient.
#'   \item A banded matrix or a graph Laplacian: neighbouring coefficients are
#'     smooth.
#' }
#'
#' @section The solver:
#' \code{"TIKNET"}, the homotopy solver, is the one solver of the selector. It
#' is the solver of the base control that applies when \code{control} is left
#' out, and a control that is passed has to name it. Every other solver is
#' refused, also \code{"TCDTIKNET"}, the deprecated coordinate-descent solver,
#' which the selector accepted up to TRexSelector 2.0.1.
#'
#' @section The dummies:
#' By default they are drawn per cluster of correlated variables, as in
#' \code{\link{TRexGVSSelector}} (\code{cluster_dummies}): the clusters are
#' \code{prior_groups}, or come from a hierarchical clustering of the
#' correlations (\code{corr_max}, \code{hc_linkage}). The field
#' \code{dummy_clusters} reports them. With an informed K, i.i.d. dummies would
#' be less competitive than correlated null variables, which deforms the
#' exchangeability the FDR calibration relies on. In the penalty every dummy
#' takes the place of the variable it was drawn for (folded coupling,
#' \code{fold_dummy_coupling}), so dummies and null variables stay exchangeable
#' for every K.
#'
#' @section L-loop strategies:
#' With cluster dummies the strategies \code{"STANDARD"}, \code{"HCONCAT"},
#' \code{"SKIPL"} and \code{"SEEDED"} of \code{\link{trex_control}()} are
#' supported, and \code{"SEEDED"} equals \code{"STANDARD"}.
#' \code{"PERMUTATION"} and \code{"PERMUTATION_SEEDED"} are refused, because a
#' row permutation destroys the covariance of the dummies within a cluster.
#' With \code{cluster_dummies = FALSE} the two run as well.
#'
#' @section lambda_2:
#' A supplied value is used as it is. Otherwise \code{lambda2_method} of
#' \code{\link{trex_tikhonov_control}()} chooses it before the L-loop, on the
#' normalized data:
#' \itemize{
#'   \item \code{"PILOT_KNEE"} (default), the pilot rule. For every value of
#'     \code{pilot_grid}, \code{pilot_K} random experiments run with
#'     \code{L_0 = p} dummies to the stopping time \code{pilot_T}, on one set
#'     of dummies, and the variables ahead of the dummies in more than half of
#'     the experiments are counted (\code{pilot_statistic}). The pick is the
#'     knee of this profile, the smallest lambda_2 whose count is within
#'     \code{pilot_knee_tolerance} of the largest. It costs about the length of
#'     the grid times \code{pilot_K} paths, a fraction of one selection per
#'     grid value. The pilot dummies come from a seed stream of their own, so
#'     the main run is a T-Rex run at a fixed lambda_2 whose dummies do not
#'     depend on the pick. The pick itself depends on y and lies outside the
#'     FDR guarantee of a fixed lambda_2. In simulations the pick cost nothing
#'     under the global null and at most 0.02 of true positive rate against
#'     the best fixed value under the alternatives. The record is the field
#'     \code{lambda2_pilot}, and
#'     \code{\link{trex_tikhonov_lambda2_power_rule}()} runs the same rule on
#'     full selections.
#'   \item \code{"CV_1SE_IEN_CCD"} and \code{"CV_MIN_IEN_CCD"}:
#'     cross-validation with the lambda_1 path profiled out per lambda_2
#'     candidate, well-posed for every p and n.
#'   \item \code{"CV_1SE_TIK_SVD"} and \code{"CV_MIN_TIK_SVD"}: analytic ridge
#'     cross-validation (lambda_1 = 0), exact for n > p and refused when the
#'     unpenalized null space interpolates the training folds.
#' }
#' The scale. A supplied \code{lambda_2} refers to the scale of the
#' standardized columns: under \code{scaling_mode = "ZSCORE"} the columns have
#' the squared norm \code{n - 1}, so the same shrinkage needs the value times
#' \code{n - 1}. \code{pilot_grid} is stated for columns of unit norm, and the
#' pilot rule and the cross-validation rules apply the factor \code{n - 1}
#' themselves. Every rule reports its pick on the scale of a supplied value
#' (the field \code{lambda2_used}).
#'
#' @section A selector runs once:
#' \code{select()} and \code{run_pilot_only()} exclude each other, and neither
#' runs a second time. Construct a new selector to run again.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' A. N. Tikhonov, "Solution of Incorrectly Formulated Problems and the Regularization Method," Soviet Math.
#'   Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl. Akad. Nauk SSSR 151 (1963), no. 3,
#'   501-504).
#'
#' Y. Park, L. Reichel, G. Rodriguez, and X. Yu, "Parameter Determination for Tikhonov Regularization Problems
#'   in General Form," J. Comput. Appl. Math., vol. 343, pp. 12-25, 2018, doi: 10.1016/j.cam.2018.04.049.
#'
#' A. J. Smola and R. Kondor, "Kernels and Regularization on Graphs," in Learning Theory and Kernel Machines
#'   (COLT/Kernel 2003), Lecture Notes in Artificial Intelligence vol. 2777, 2003, pp. 144-158.
#'
#' F. R. K. Chung, Spectral Graph Theory. CBMS Regional Conference Series in Mathematics, no. 92. American
#'   Mathematical Society, 1997.
#'
#' A. E. Hoerl and R. W. Kennard, "Ridge Regression: Biased Estimation for Nonorthogonal Problems,"
#'   Technometrics, vol. 42, no. 1, pp. 80-86, 2000 (Reprint of the 1970 article in the Special 40th
#'   Anniversary Issue).
#'
#' J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
#'   High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
#'   Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
#'
#' @examples
#' \donttest{
#' data(Gauss_data)
#' X <- Gauss_data$X
#' y <- Gauss_data$y
#' ctrl <- trex_tikhonov_control(tikhonov_K = diag(ncol(X)), lambda_2 = 0.5)
#' sel <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, verbose = FALSE,
#'                                 tik_control = ctrl)
#' sel$select()
#' sel$selected_indices
#' }
#' @family T-Rex selectors and controls
#' @export
TRexTikhonovSelector <- R6::R6Class("TRexTikhonovSelector",
  inherit = TRexSelector,
  public = list(
    #' @description Create a new TRexTikhonovSelector object.
    #'
    #' @param X Feature matrix in memory, or an \code{mmap_matrix}, the handle
    #'   of a memory-mapped file. Both are bound without a copy. The selector
    #'   standardizes a memory-mapped design in place and restores it, which
    #'   needs \code{"readwrite"} mode. Under \code{scaling_mode = "NONE"} it
    #'   writes nothing and \code{"readonly"} mode suffices.
    #' @param y Response vector.
    #' @param tFDR Target FDR level between 0 and 1 (default: 0.1).
    #' @param seed Random seed. A negative value (default: -1) gives a run that
    #'   is not reproducible.
    #' @param verbose Whether to print progress (default: TRUE).
    #' @param tik_control A Tikhonov control list from
    #'   \code{\link{trex_tikhonov_control}()}. Required, because it carries the
    #'   penalty matrix \code{tikhonov_K} or the penalty structure.
    #' @param control A control list from \code{\link{trex_control}()}, or
    #'   \code{NULL} (default). With \code{NULL} the base control of the C++
    #'   class applies: the defaults of \code{trex_control()} with the solver
    #'   \code{"TIKNET"}, the pathwise homotopy solver. A control that is
    #'   passed is taken as it is, so its \code{solver} must be
    #'   \code{"TIKNET"}: \code{trex_control(solver = "TIKNET", ...)}. Every
    #'   other solver is an error, also \code{"TLARS"}, the default of
    #'   \code{trex_control()}, and \code{"TCDTIKNET"}, the deprecated
    #'   coordinate-descent solver, which the selector accepted up to
    #'   TRexSelector 2.0.1.
    initialize = function(X, y,
                          tFDR = 0.1,
                          seed = -1,
                          verbose = TRUE,
                          tik_control,
                          control = NULL) {

      if (missing(tik_control)) {
        stop("tik_control is required (build it with trex_tikhonov_control(); ",
             "it carries the penalty matrix K or the penalty structure).")
      }
      .check_tikhonov_solver(control)
      .check_tikhonov_groups_length(tik_control, X)

      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_tikhonov_mmap_create(
          get_ptr(X), y, tFDR, tik_control, control, seed, verbose
        )
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_tikhonov_create(
          X, y, tFDR, tik_control, control, seed, verbose
        )
      }
    },

    #' @description Run the T-Rex+Tikhonov selection. A selector runs once: a
    #'   second call is an error, and so is a call after
    #'   \code{run_pilot_only()}, see \code{\link{TRexSelector}}.
    select = function() {
      trex_tikhonov_select(private$ptr)
      invisible(self)
    },

    #' @description Run only the lambda_2 pilot, without the selection. The
    #'   pilot builds the cluster model, runs the pilot experiments and picks
    #'   lambda_2, and \code{X} is given back. It is for a caller who wants
    #'   the pick alone. It needs \code{lambda_2 < 0} and
    #'   \code{lambda2_method = "PILOT_KNEE"} in \code{tik_control} and is an
    #'   error otherwise. A selector runs once: after \code{run_pilot_only()}
    #'   a call of \code{select()} or of \code{run_pilot_only()} is an error,
    #'   and so is \code{run_pilot_only()} after \code{select()}.
    #' @return The record of the pilot, as the field \code{lambda2_pilot}
    #'   holds it afterwards.
    run_pilot_only = function() {
      trex_tikhonov_run_pilot_only(private$ptr)
    }
  ),
  active = list(
    #' @field lambda2_used The lambda_2 the last run used: the supplied value,
    #'   or the pick of the pilot rule or of a cross-validation rule, on the
    #'   scale of a supplied value. \code{-1} before a run.
    lambda2_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_used(private$ptr)
    },

    #' @field lambda2_cv_diagnostics Health record of the rule that chose
    #'   lambda_2 in the last run: a list with \code{ran} (a cross-validation
    #'   rule or the pilot ran, \code{FALSE} when the value was supplied),
    #'   \code{seconds}, \code{n_profile_at_edge} (lambda_2 values whose
    #'   profiled lambda_1 optimum sat on a truncation edge; 0 is healthy),
    #'   \code{support_cap} (-1 = disabled), \code{dfmax_ratio} (the
    #'   \code{cv_dfmax_ratio} the IEN rules were called with; 0 otherwise),
    #'   \code{low_dimensional}, and the fields of the boundary-pick guard,
    #'   \code{grid_index} and \code{grid_size}, \code{at_grid_edge},
    #'   \code{profile_at_edge}, \code{profile_at_support_cap} and
    #'   \code{boundary_pick} (see \code{\link{TRexGVSSelector}}). The
    #'   cross-validation rules warn once on a boundary pick and return it
    #'   unchanged. The pilot rule leaves the guard fields at their defaults.
    lambda2_cv_diagnostics = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_cv_diagnostics(private$ptr)
    },

    #' @field lambda2_pilot Record of the lambda_2 pilot of the last
    #'   \code{select()} or \code{run_pilot_only()}, a list:
    #'   \describe{
    #'     \item{\code{ran}}{A pilot ran. \code{FALSE}: lambda_2 came from the
    #'       caller or from a cross-validation rule.}
    #'     \item{\code{grid}}{The lambda_2 values the pilot ran at, ascending
    #'       and without duplicates: \code{pilot_grid} times \code{grid_scale}.}
    #'     \item{\code{grid_scale}}{The factor between \code{pilot_grid} and
    #'       \code{grid}: \code{n - 1} under \code{scaling_mode = "ZSCORE"}, 1
    #'       otherwise.}
    #'     \item{\code{counts}}{The count the pick was read off, per grid
    #'       value: the count of \code{pilot_statistic}.}
    #'     \item{\code{curve_majority}}{The majority count per grid value (row)
    #'       and pilot stopping time 1 to \code{T_pilot} (column).}
    #'     \item{\code{curve_soft}}{The mean number of variables ahead of the
    #'       t-th dummy, per grid value (row) and t (column). A diagnostic.}
    #'     \item{\code{calibrated}}{The calibrated pilot selection count per
    #'       grid value.}
    #'     \item{\code{group_majority}}{The group-majority count per grid
    #'       value: the members ahead, summed over the groups that reach
    #'       \code{pilot_group_fraction}. Zeros without \code{prior_groups}.}
    #'     \item{\code{groups_qualifying}}{The number of groups that reach
    #'       \code{pilot_group_fraction}, per grid value. Zeros without
    #'       \code{prior_groups}.}
    #'     \item{\code{probe_seconds}}{Wall-clock seconds per grid value.}
    #'     \item{\code{index}}{Index of the pick into \code{grid}, counted
    #'       from 1.}
    #'     \item{\code{lambda2}}{The picked lambda_2, \code{grid[index]}, on
    #'       the working scale.}
    #'     \item{\code{max_selected}}{The largest count on the grid.}
    #'     \item{\code{n_ties}}{Grid values within the tolerance of the largest
    #'       count.}
    #'     \item{\code{at_lower_end}}{The pick sits on the first grid value.}
    #'     \item{\code{at_upper_end}}{The pick sits on the last grid value:
    #'       widen the grid before trusting it.}
    #'     \item{\code{K_pilot}}{Pilot experiments per grid value that ran:
    #'       \code{pilot_K}, capped at \code{K}.}
    #'     \item{\code{T_pilot}}{Stopping time of the pilot paths that was
    #'       used: \code{pilot_T}, capped at \code{min(L_0, ceiling(n / 2))}.}
    #'     \item{\code{L_0}}{Number of pilot dummies, \code{pilot_L_multiplier}
    #'       times p.}
    #'     \item{\code{dummy_tag}}{Seed tag under which the pilot dummies were
    #'       drawn.}
    #'     \item{\code{seconds}}{Wall-clock seconds of the whole pilot.}
    #'   }
    lambda2_pilot = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_pilot_result(private$ptr)
    },

    #' @field dummy_clusters The clusters the dummies of the last run were
    #'   drawn for, a list with the names of the fields of
    #'   \code{\link{TRexGVSSelector}}: \code{groups} (the cluster of each
    #'   variable, numbered from 1), \code{max_clusters} (their number M) and
    #'   \code{hc_method_used} (\code{"prior"} when \code{prior_groups} was
    #'   supplied, else the linkage of the clustering). Before a run and with
    #'   \code{cluster_dummies = FALSE} the record is empty: no groups, 0 and
    #'   \code{""}.
    dummy_clusters = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_dummy_clusters(private$ptr)
    }
  )
)


#' @name trex_tikhonov_lambda2_power_rule
#'
#' @title The knee of the full lambda_2 profile of the T-Rex+Tikhonov selector
#'
#' @description Runs one FDR-controlled selection of
#'   \code{\link{TRexTikhonovSelector}} per value of \code{grid} with the same
#'   seed, so the probes share their dummies and differ by the penalty weight
#'   alone, and picks the smallest lambda_2 whose selection is within
#'   \code{knee_tolerance} selections of the largest (the knee of the profile).
#'   This is the reference rule behind the default
#'   \code{lambda2_method = "PILOT_KNEE"} of
#'   \code{\link{trex_tikhonov_control}}, which reads the same knee off a few
#'   random experiments per grid value at a fraction of the cost. The full
#'   profile costs one selection per grid value and returns every probe's
#'   selection. The pick depends on y and lies outside the fixed-lambda_2 FDR
#'   guarantee that each probe carries.
#'   \code{\link{trex_gvs_lambda2_power_rule}} is the same rule on the
#'   T-Rex+GVS selector.
#'
#' @param X Feature matrix (n x p) in memory, or an \code{mmap_matrix}, the
#'   handle of a memory-mapped file. Every probe normalizes it in place and
#'   hands it back, so \code{X} is unchanged on return and never copied. A
#'   memory-mapped design needs \code{"readwrite"} mode for that, and
#'   \code{"readonly"} mode suffices under \code{scaling_mode = "NONE"}.
#' @param y Response vector (length n).
#' @param tFDR Target FDR level (default: 0.1).
#' @param tik_control A Tikhonov control list from
#'   \code{\link{trex_tikhonov_control}()}. Its \code{lambda_2} is replaced
#'   per probe, everything else is shared by every probe.
#' @param control A control list from \code{\link{trex_control}()}, or
#'   \code{NULL} (default), taken as by \code{\link{TRexTikhonovSelector}}:
#'   with \code{NULL} the base control of the C++ class applies, and a
#'   control that is passed must name the solver \code{"TIKNET"}.
#' @param seed Selector seed, the same for every probe (default: 0). A
#'   negative seed would draw fresh dummies per probe and is refused.
#' @param grid lambda_2 candidates, stated for columns of unit norm, sorted
#'   and de-duplicated. The default is the grid of the pilot rule,
#'   \code{\link{trex_default_lambda2_grid}()}: 0, 1, 3, 10, 30, 100, 300,
#'   1000, 3000. Under \code{scaling_mode = "ZSCORE"} the probes run at
#'   these values times \code{n - 1}, as the pilot rule does. Keep 0 on the
#'   grid so that "no penalty" is an admissible answer.
#' @param knee_tolerance Selections below the largest count that still
#'   count as "on the plateau" (default: 1). \code{0} reads the largest
#'   selection, ties toward the smallest lambda_2.
#' @param verbose Whether the probes print their progress (default:
#'   \code{FALSE}).
#'
#' @return A list:
#'   \describe{
#'     \item{\code{lambda2}}{The pick, on the working scale of the runs.
#'       Supplied as \code{lambda_2}, it reproduces the selection at the pick.}
#'     \item{\code{index}}{Index of the pick into \code{grid}, counted from 1.}
#'     \item{\code{selected}}{The selection at the pick (indices from 1,
#'       ascending).}
#'     \item{\code{grid}}{The values the probes ran at, ascending and without
#'       duplicates: the supplied grid times \code{grid_scale}.}
#'     \item{\code{grid_scale}}{The factor between the supplied grid and
#'       \code{grid}: \code{n - 1} under \code{scaling_mode = "ZSCORE"}, 1
#'       otherwise.}
#'     \item{\code{counts}}{The selection sizes per grid value.}
#'     \item{\code{probes}}{One list per grid value with \code{lambda2},
#'       \code{selected}, \code{T_stop} (the stopping time of the run) and
#'       \code{seconds}.}
#'     \item{\code{max_selected}}{The largest selection on the grid.}
#'     \item{\code{n_ties}}{Grid values within the tolerance of the largest
#'       selection.}
#'     \item{\code{at_lower_end}}{The pick sits on the first grid value.}
#'     \item{\code{at_upper_end}}{The pick sits on the last grid value, a
#'       boundary solution: widen the grid before trusting it.}
#'     \item{\code{seconds}}{Wall-clock seconds of the whole profile.}
#'   }
#'
#' @examples
#' \donttest{
#' data(Gauss_data)
#' X <- Gauss_data$X
#' y <- Gauss_data$y
#' ctrl <- trex_tikhonov_control(tikhonov_K = diag(ncol(X)))
#' prof <- trex_tikhonov_lambda2_power_rule(X, y, tFDR = 0.1, tik_control = ctrl,
#'                                          grid = c(0, 1, 10), seed = 1)
#' prof$lambda2
#' prof$counts
#' }
#' @family T-Rex selectors and controls
#' @export
trex_tikhonov_lambda2_power_rule <- function(X, y,
                                             tFDR = 0.1,
                                             tik_control,
                                             control = NULL,
                                             seed = 0,
                                             grid = trex_default_lambda2_grid(),
                                             knee_tolerance = 1,
                                             verbose = FALSE) {
  if (missing(tik_control)) {
    stop("tik_control is required (build it with trex_tikhonov_control(); ",
         "it carries the penalty matrix K or the penalty structure).")
  }
  .check_power_rule_arguments(seed, grid, knee_tolerance)
  .check_tikhonov_solver(control)
  .check_tikhonov_groups_length(tik_control, X)
  y <- .as_double_vector(y)
  if (inherits(X, "MemoryMappedMatrix")) {
    return(trex_tikhonov_lambda2_power_rule_mmap_impl(
      get_ptr(X), y, tFDR, tik_control, control, as.integer(seed),
      as.numeric(grid), as.numeric(knee_tolerance), verbose
    ))
  }
  X <- .as_double_matrix(X)
  trex_tikhonov_lambda2_power_rule_impl(X, y, tFDR, tik_control, control,
                                        as.integer(seed), as.numeric(grid),
                                        as.numeric(knee_tolerance), verbose)
}


# The arguments the two lambda_2 power rules share.
.check_power_rule_arguments <- function(seed, grid, knee_tolerance) {
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed < 0) {
    stop("seed must be a single non-negative integer: the probes share ",
         "their dummies. Got: ", format(seed))
  }
  if (!is.numeric(grid) || length(grid) == 0 || any(!is.finite(grid)) ||
      any(grid < 0)) {
    stop("grid must be a non-empty vector of finite values >= 0.")
  }
  if (!is.numeric(knee_tolerance) || length(knee_tolerance) != 1 ||
      !is.finite(knee_tolerance) || knee_tolerance < 0) {
    stop("knee_tolerance must be a single number >= 0. Got: ",
         format(knee_tolerance))
  }
  invisible(NULL)
}


#' @name trex_default_lambda2_grid
#'
#' @title The default lambda_2 grid of the pilot rule and of the power rules
#'
#' @description Returns the grid that the C++ core holds once: 0, 1, 3, 10, 30,
#'   100, 300, 1000, 3000. It is the default of \code{pilot_grid} in
#'   \code{\link{trex_tikhonov_control}} and of \code{grid} in
#'   \code{\link{trex_tikhonov_lambda2_power_rule}} and
#'   \code{\link{trex_gvs_lambda2_power_rule}}. The values are stated for
#'   columns of unit norm.
#'
#' @return A numeric vector, ascending.
#'
#' @examples
#' trex_default_lambda2_grid()
#'
#' @family T-Rex selectors and controls
#' @export
trex_default_lambda2_grid <- function() {
  trex_default_lambda2_grid_impl()
}


#' @name trex_gamma_to_k
#'
#' @title Build a Tikhonov penalty matrix K from an operator Gamma
#'
#' @description Computes \code{K = t(Gamma) \%*\% Gamma} (pruned and
#'   compressed) for use as \code{tikhonov_K} in
#'   \code{\link{trex_tikhonov_control}}. For example, a first-difference
#'   Gamma yields the smoothness-encoding graph Laplacian of a chain.
#'
#' @param gamma Regularization operator (m x p): a base \code{matrix} or a
#'   sparse \code{Matrix::dgCMatrix}.
#'
#' @return A sparse \code{Matrix::dgCMatrix} (p x p). Requires the Matrix
#'   package (Suggests) for the sparse return class.
#'
#' @examples
#' # First-difference operator over 5 variables -> chain Laplacian
#' if (requireNamespace("Matrix", quietly = TRUE)) {
#'   p <- 5
#'   gamma <- diff(diag(p))
#'   K <- trex_gamma_to_k(gamma)
#' }
#'
#' @family T-Rex selectors and controls
#' @export
trex_gamma_to_k <- function(gamma) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("trex_gamma_to_k() returns a sparse Matrix::dgCMatrix and needs the ",
         "Matrix package installed.")
  }
  trex_tikhonov_gamma_to_k(gamma)
}


#' @name trex_build_group_centering_k
#'
#' @title Build the centering penalty matrix of a group assignment
#'
#' @description Computes the block-diagonal matrix
#'   \code{K = blockdiag(I - J / p_m)} over the groups, with \code{p_m} the
#'   size of group m and J the matrix of ones: the penalty under which the
#'   members of a group share a coefficient. It is the matrix that
#'   \code{penalty_structure = "GROUP_CENTERING"} of
#'   \code{\link{trex_tikhonov_control}} stands for, so a selector with that
#'   structure needs the groups alone. A group with one member contributes
#'   nothing.
#'
#' @param prior_groups Integer vector with the group of each variable,
#'   numbered 1 to M without gaps (length p).
#'
#' @return A sparse \code{Matrix::dgCMatrix} (p x p). Requires the Matrix
#'   package (Suggests) for the sparse return class.
#'
#' @examples
#' if (requireNamespace("Matrix", quietly = TRUE)) {
#'   K <- trex_build_group_centering_k(rep(1:2, each = 3))
#' }
#'
#' @family T-Rex selectors and controls
#' @export
trex_build_group_centering_k <- function(prior_groups) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("trex_build_group_centering_k() returns a sparse Matrix::dgCMatrix and ",
         "needs the Matrix package installed.")
  }
  .check_group_numbers(prior_groups, "prior_groups")
  trex_tikhonov_group_centering_k(as.integer(prior_groups))
}
