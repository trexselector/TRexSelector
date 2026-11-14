# The solver of a base control that is passed to the Tikhonov selector. NULL
# stands for the base control of the C++ class, whose solver is TIKNET.
.check_tikhonov_solver <- function(control) {
  if (is.null(control)) {
    return(invisible(NULL))
  }
  solver <- control$solver
  if (!is.character(solver) || length(solver) != 1L ||
      !(solver %in% c("TIKNET", "TCDTIKNET"))) {
    stop("control$solver must be \"TIKNET\" for the Tikhonov selector. Got: \"",
         paste(format(solver), collapse = ", "), "\". trex_control() names ",
         "\"TLARS\" by default, so pass trex_control(solver = \"TIKNET\", ...), ",
         "or leave control out to run the default control.")
  }
  .warn_cd_solver_name(solver)
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
#' @title T-Rex Selector with a general Tikhonov penalty
#'
#' @description FDR-controlled variable selection with a Tikhonov penalty
#'   matrix \code{K = t(Gamma) \%*\% Gamma} driving the Tikhonov-net solvers
#'   (TIKNET by default, or the deprecated TCDTIKNET). \code{K = I} collapses to
#'   the plain elastic net; a group-mean K reproduces the GVS-IEN penalty;
#'   banded / Laplacian K encode smoothness-informed selection;
#'   \code{penalty_structure = "GROUP_CENTERING"} builds the centering
#'   penalty from the dummy clustering. lambda_2 comes from the caller, a
#'   CV tuner, or the pilot rule (see \code{\link{trex_tikhonov_control}}).
#'   Dummies are cluster-aware and FOLDED-coupled by default.
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
    #' @param tFDR Target FDR level (default: 0.1).
    #' @param seed Random seed (default: -1).
    #' @param verbose Whether to print progress (default: TRUE).
    #' @param tik_control A Tikhonov control list from
    #'   \code{\link{trex_tikhonov_control}()} (required — it carries the
    #'   penalty matrix K).
    #' @param control A control list from \code{\link{trex_control}()}, or
    #'   \code{NULL} (default). With \code{NULL} the base control of the C++
    #'   class applies: the defaults of \code{trex_control()} with the solver
    #'   \code{"TIKNET"}, the pathwise homotopy solver. A control that is
    #'   passed is taken as it is, so its \code{solver} must be
    #'   \code{"TIKNET"}: \code{trex_control(solver = "TIKNET", ...)}. Every
    #'   other solver is an error, also \code{"TLARS"}, the default of
    #'   \code{trex_control()}. \code{"TCDTIKNET"}, the coordinate-descent
    #'   solver, is deprecated and will be removed in TRexSelector 2.1.0: the
    #'   selector raises a deprecation warning for it. It is not finally
    #'   developed: once \code{lambda_2} times K is large it can stop at its
    #'   sweep budget (\code{cd_max_sweeps}) before it has converged, and
    #'   \code{select()} then stops with an error instead of a selection from
    #'   an unconverged path.
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

    #' @description Run the Tikhonov T-Rex algorithm. A selector runs once: a second call is an
    #'   error, see \code{\link{TRexSelector}}.
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
    #' @field lambda2_used The lambda_2 actually used by the last
    #'   \code{select()} run (user-supplied or CV-resolved; solver scale).
    #'   \code{-1} before \code{select()}.
    lambda2_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_used(private$ptr)
    },

    #' @field lambda2_cv_diagnostics Health record of the lambda_2 tuner of
    #'   the last \code{select()}: a list with \code{ran} (a tuner or the
    #'   pilot ran), \code{seconds}, \code{n_profile_at_edge} (lambda_2
    #'   values whose profiled lambda_1 optimum sat on a truncation edge; 0
    #'   is healthy), \code{support_cap} (-1 = disabled),
    #'   \code{dfmax_ratio} (IEN tuner: the \code{cv_dfmax_ratio} it was
    #'   called with; 0 otherwise), \code{low_dimensional}, and the
    #'   boundary-pick guard fields \code{grid_index} / \code{grid_size},
    #'   \code{at_grid_edge}, \code{profile_at_edge},
    #'   \code{profile_at_support_cap} and \code{boundary_pick} (see
    #'   \code{TRexGVSSelector}; the CV methods warn once on a boundary pick
    #'   and return it unchanged; the pilot rule leaves them at their
    #'   defaults).
    lambda2_cv_diagnostics = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_cv_diagnostics(private$ptr)
    },

    #' @field lambda2_pilot Record of the lambda_2 pilot rule of the last
    #'   \code{select()} (\code{ran == FALSE} when no pilot ran): the
    #'   \code{grid} of the values the pilot ran at (\code{pilot_grid} times
    #'   \code{grid_scale}, which is \code{n - 1} under
    #'   \code{scaling_mode = "ZSCORE"} and 1 otherwise), the \code{counts}
    #'   the pick was read off, the
    #'   per-grid-point curves over the pilot stopping times
    #'   (\code{curve_majority}, \code{curve_soft}; G x T matrices), the
    #'   \code{calibrated} and \code{group_majority} counts, the 1-based
    #'   \code{index} of the pick, the picked \code{lambda2},
    #'   \code{max_selected}, \code{n_ties}, the boundary flags
    #'   \code{at_lower_end} / \code{at_upper_end} (widen the grid before
    #'   trusting a boundary pick), the budget actually used (\code{K_pilot},
    #'   \code{T_pilot}, \code{L_0}), the seed tag under which the pilot
    #'   dummies were drawn (\code{dummy_tag}) and the timings.
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
#' @title The knee of the full lambda_2 profile of the Tikhonov T-Rex
#'
#' @description Runs one FDR-controlled Tikhonov T-Rex selection per value of
#'   \code{grid} with the same seed, so the probes share their dummies and
#'   differ by the penalty weight alone, and picks the smallest lambda_2 whose
#'   selection is within \code{knee_tolerance} selections of the largest (the
#'   knee of the profile). This is the reference rule behind the default
#'   \code{lambda2_method = "PILOT_KNEE"} of
#'   \code{\link{trex_tikhonov_control}}, which reads the same knee off a few
#'   random experiments per grid value at a fraction of the cost; the full
#'   profile costs one selection per grid value and returns every probe's
#'   selection. The pick depends on y and lies outside the fixed-lambda_2 FDR
#'   guarantee that each probe carries.
#'
#' @param X Feature matrix (n x p) in memory, or an \code{mmap_matrix}, the
#'   handle of a memory-mapped file. Every probe normalizes it in place and
#'   hands it back, so \code{X} is unchanged on return and never copied. A
#'   memory-mapped design needs \code{"readwrite"} mode for that, and
#'   \code{"readonly"} mode suffices under \code{scaling_mode = "NONE"}.
#' @param y Response vector (length n).
#' @param tFDR Target FDR level (default: 0.1).
#' @param tik_control A Tikhonov control list from
#'   \code{\link{trex_tikhonov_control}()}; its \code{lambda_2} is
#'   overwritten per probe, everything else is shared by every probe.
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
#'   count as "on the plateau" (default: 1); 0 reads the largest selection,
#'   ties toward the smallest lambda_2.
#' @param verbose Whether the probes print their progress (default:
#'   \code{FALSE}).
#'
#' @return A list: \code{lambda2} (the pick, on the working scale: supplied as
#'   \code{lambda_2} it reproduces the selection), \code{index} (1-based, into
#'   \code{grid}), \code{selected} (the 1-based selection at the pick),
#'   \code{grid} (the values the probes ran at, ascending and de-duplicated),
#'   \code{grid_scale} (the factor between the supplied grid and \code{grid},
#'   \code{n - 1} under \code{"ZSCORE"} and 1 otherwise), \code{counts} (selection sizes
#'   per grid value), \code{probes} (one list per grid value with
#'   \code{lambda2}, \code{selected}, \code{T_stop} and \code{seconds}),
#'   \code{max_selected}, \code{n_ties} (grid values within the tolerance
#'   of the maximum), the boundary flags \code{at_lower_end} /
#'   \code{at_upper_end} (widen the grid before trusting a boundary pick)
#'   and \code{seconds}.
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
