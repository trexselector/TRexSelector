# The solver of a GVS run, checked in the words of R. NULL stands for the base
# control of the C++ class, whose solver is "TENET". Nothing is changed: EN
# takes its solver from en_solver, IEN needs it named in the base control.
.check_gvs_solver <- function(gvs_control, control) {
  solver <- if (is.null(control)) "TENET" else control$solver
  got <- if (is.null(control)) {
    "No control was passed."
  } else {
    paste0("Got: \"", paste(format(solver), collapse = ", "), "\".")
  }
  named <- is.character(solver) && length(solver) == 1L
  if (gvs_control$gvs_type == "IEN") {
    if (!named || !(solver %in% c("TIENET", "TIENET_AUG", "TCDIENET"))) {
      stop("gvs_type = \"IEN\" needs a base control that names its solver, ",
           "control = trex_control(solver = \"TIENET\", ...) or ",
           "\"TIENET_AUG\". ", got)
    }
    .warn_cd_solver_name(solver)
  } else {
    fits <- unique(c("TENET", gvs_control$en_solver))
    if (!named || !(solver %in% fits)) {
      stop("control$solver must be ", paste0("\"", fits, "\"", collapse = " or "),
           " for gvs_type = \"EN\". ", got, " The elastic-net solver is chosen ",
           "with en_solver of trex_gvs_control(). trex_control() names ",
           "\"TLARS\" by default, so pass trex_control(solver = \"TENET\", ...), ",
           "or leave control out to run the default control.")
    }
    .warn_cd_solver_name(gvs_control$en_solver, "en_solver")
  }
  invisible(NULL)
}

# One group entry per variable, checked with the name of the R argument.
.check_gvs_groups_length <- function(gvs_control, X) {
  if (!is.null(gvs_control$groups) && length(gvs_control$groups) != ncol(X)) {
    stop("gvs_control$groups must have one entry per variable (p = ", ncol(X),
         "). Got: ", length(gvs_control$groups), ".")
  }
  invisible(NULL)
}


#' @name TRexGVSSelector
#'
#' @title T-Rex Selector for Grouped Variable Selection (T-Rex+GVS)
#'
#' @description FDR-controlled variable selection for groups of highly correlated variables.
#'
#' @section L-loop strategies:
#' The selector supports the strategies \code{"STANDARD"}, \code{"HCONCAT"}, \code{"SKIPL"} and
#' \code{"SEEDED"} of \code{\link{trex_control}()}, and \code{"SEEDED"} equals \code{"STANDARD"}
#' here. \code{"PERMUTATION"} and \code{"PERMUTATION_SEEDED"} are refused, because a row permutation
#' destroys the covariance of the dummies within a cluster.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Grouped Variable Selection in
#'   High-Dimensional Linear Models Using the T-Knock Filter," in Proc. 30th Eur. Signal Process. Conf.
#'   (EUSIPCO), 2022, pp. 892-896.
#'
#' J. Machkour, M. Muma, and D. P. Palomar, "The Informed Elastic Net for Fast Grouped Variable Selection and
#'   FDR Control in Genomics Research," in Proc. IEEE 9th Int. Workshop Comput. Adv. Multi-Sensor Adapt.
#'   Process. (CAMSAP), 2023, pp. 466-470.
#'
#' H. Zou and T. Hastie, "Regularization and Variable Selection via the Elastic Net," J. R. Stat. Soc. Ser. B
#'   Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.
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
#' sel <- TRexGVSSelector$new(X, y, tFDR = 0.1, verbose = FALSE,
#'                            gvs_control = trex_gvs_control(gvs_type = "EN"))
#' sel$select()
#' sel$selected_indices
#' }
#' @family T-Rex selectors and controls
#' @export
TRexGVSSelector <- R6::R6Class("TRexGVSSelector",
  inherit = TRexSelector,
  public = list(
    #' @description Create a new TRexGVSSelector object.
    #'
    #' @param X Feature matrix in memory, or an \code{mmap_matrix}, the handle
    #'   of a memory-mapped file. Both are bound without a copy. The selector
    #'   standardizes a memory-mapped design in place and restores it, which
    #'   needs \code{"readwrite"} mode. Under \code{scaling_mode = "NONE"} it
    #'   writes nothing and \code{"readonly"} mode suffices.
    #' @param y Response vector.
    #' @param tFDR Target FDR level (default: 0.1).
    #' @param seed Random seed. A negative value (default: -1) gives a run that is not
    #'   reproducible.
    #' @param verbose Whether to print progress (default: TRUE).
    #' @param gvs_control A GVS control list from \code{\link{trex_gvs_control}()} (default:
    #'   \code{trex_gvs_control()}). For \code{gvs_type = "EN"} its \code{en_solver}
    #'   chooses the solver, "TENET" by default. "TCDENET" and "TCDIENET" are
    #'   deprecated and will be removed in TRexSelector 2.1.0: the selector raises a
    #'   deprecation warning for them.
    #' @param control A control list from \code{\link{trex_control}()}, or \code{NULL}
    #'   (default). With \code{NULL} the base control of the C++ class applies: the
    #'   defaults of \code{trex_control()} with the solver "TENET". A control that is
    #'   passed is taken as it is, so it has to name a solver of the track. For
    #'   \code{gvs_type = "EN"} its \code{solver} must be "TENET" or the solver of
    #'   \code{en_solver}. For \code{gvs_type = "IEN"} it must be "TIENET" (native
    #'   pathwise) or "TIENET_AUG" (row-augmented), so the IEN track needs a base
    #'   control: \code{trex_control(solver = "TIENET", ...)}. Every other solver is an
    #'   error, also "TLARS", the default of \code{trex_control()}.
    initialize = function(X, y,
                          tFDR = 0.1,
                          seed = -1,
                          verbose = TRUE,
                          gvs_control = trex_gvs_control(),
                          control = NULL) {

      .check_gvs_solver(gvs_control, control)

      # lambda_2 == 0 is the degenerate no-ridge case, NOT automatic CV
      # (the sentinel for that is lambda_2 < 0). A user passing exactly 0
      # in the belief it triggers CV silently gets pure T-LASSO.
      if (isTRUE(gvs_control$lambda_2 == 0)) {
        warning("gvs_control$lambda_2 = 0 is the degenerate no-ridge case ",
                "(pure T-LASSO), not automatic cross-validation. ",
                "Use lambda_2 < 0 (e.g. -1) to auto-determine lambda_2 via CV.")
      }

      .check_gvs_groups_length(gvs_control, X)

      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- .r_positions(trex_gvs_mmap_create(
          get_ptr(X), y, tFDR, gvs_control, control, seed, verbose
        ))
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- .r_positions(trex_gvs_create(
          X, y, tFDR, gvs_control, control, seed, verbose
        ))
      }
    },

    #' @description Run the GVS algorithm. A selector runs once: a second call is an
    #'   error, see \code{\link{TRexSelector}}.
    select = function() {
      trex_gvs_select(private$ptr)
      invisible(self)
    }
  ),
  active = list(
    #' @field lambda2_used \code{lambda_2} of the last \code{select()}: the supplied value,
    #'   or the pick of the cross-validation rule \code{lambda2_method}.
    lambda2_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_lambda2_used(private$ptr)
    },

    #' @field lambda2_cv_diagnostics Health record of the lambda_2 tuner of
    #'   the last \code{select()}: a list with \code{ran}, \code{seconds},
    #'   \code{n_profile_at_edge} (lambda_2 values whose profiled lambda_1
    #'   optimum sat on a truncation edge; 0 is healthy), \code{support_cap}
    #'   (-1 = disabled), \code{dfmax_ratio} (IEN tuner: the
    #'   \code{cv_dfmax_ratio} it was called with; 0 when no IEN-geometry
    #'   tuner ran), \code{low_dimensional}, and the boundary-pick guard
    #'   fields: \code{grid_index} / \code{grid_size} (0-based position of
    #'   the returned lambda_2 on the tuner's grid), \code{at_grid_edge} (the
    #'   pick is the first or last grid point, so the grid boundary rather
    #'   than a CV optimum), \code{profile_at_edge} (IEN tuner: the profiled
    #'   lambda_1 optimum of the chosen lambda_2 sat on the truncation edge),
    #'   \code{profile_at_support_cap} (that edge is the support cap; raise
    #'   \code{cv_dfmax_ratio} to move it) and \code{boundary_pick}
    #'   (\code{at_grid_edge} or \code{profile_at_edge}; the selector warned
    #'   once and returned the pick unchanged).
    lambda2_cv_diagnostics = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_lambda2_cv_diagnostics(private$ptr)
    },

    #' @field gvs_type Penalty of the grouped selection that was applied, \code{"EN"} or
    #'   \code{"IEN"}.
    gvs_type = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_gvs_type(private$ptr)
    },

    #' @field max_clusters Number of groups M, supplied or found by the clustering.
    max_clusters = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_max_clusters(private$ptr)
    },

    #' @field groups Integer vector of length \code{p} giving the group of each
    #'   variable, numbered from 1: the supplied groups, or the assignment of the
    #'   clustering.
    groups = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_groups(private$ptr)
    },

    #' @field group_labels Optional human-readable cluster names (length
    #'   \code{max_clusters}), or an empty character vector when none were
    #'   supplied.
    group_labels = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_gvs_get_group_labels(private$ptr)
    },

    #' @field hc_method_used Origin of the groups: \code{"prior"} when they were
    #'   supplied, else the linkage of the clustering, \code{"single"},
    #'   \code{"complete"}, \code{"average"} or \code{"wpgma"}.
    hc_method_used = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      trex_gvs_get_hc_method_used(private$ptr)
    }
  )
)


#' @name trex_gvs_lambda2_power_rule
#'
#' @title The knee of the full lambda_2 profile of the T-Rex+GVS selector
#'
#' @description Runs one FDR-controlled selection of
#'   \code{\link{TRexGVSSelector}} per value of \code{grid} for
#'   \code{lambda_2}, the ridge part of the elastic net or of the informed
#'   elastic net, with the same seed, so the probes share their dummies and
#'   differ by the penalty weight alone. The pick is the smallest lambda_2
#'   whose selection is within \code{knee_tolerance} selections of the largest
#'   (the knee of the profile). The result holds every probe's selection. The
#'   pick depends on y and lies outside the fixed-lambda_2 FDR guarantee that
#'   each probe carries. \code{\link{trex_tikhonov_lambda2_power_rule}} is the
#'   same rule on the Tikhonov selector.
#'
#' @param X Feature matrix (n x p) in memory, or an \code{mmap_matrix}, the
#'   handle of a memory-mapped file. Every probe normalizes it in place and
#'   hands it back, so \code{X} is unchanged on return and never copied. A
#'   memory-mapped design needs \code{"readwrite"} mode for that, and
#'   \code{"readonly"} mode suffices under \code{scaling_mode = "NONE"}.
#' @param y Response vector (length n).
#' @param tFDR Target FDR level (default: 0.1).
#' @param gvs_control A GVS control list from \code{\link{trex_gvs_control}()}
#'   (default: \code{trex_gvs_control()}); its \code{lambda_2} is replaced
#'   per probe, everything else is shared by every probe.
#' @param control A control list from \code{\link{trex_control}()}, or
#'   \code{NULL} (default), taken as by \code{\link{TRexGVSSelector}}: with
#'   \code{NULL} the base control of the C++ class applies, and a control
#'   that is passed has to name a solver of the track.
#' @param seed Selector seed, the same for every probe (default: 0). A
#'   negative seed would draw fresh dummies per probe and is refused.
#' @param grid lambda_2 candidates, stated for columns of unit norm, sorted
#'   and de-duplicated. The default is
#'   \code{\link{trex_default_lambda2_grid}()}: 0, 1, 3, 10, 30, 100, 300,
#'   1000, 3000. Under \code{scaling_mode = "ZSCORE"} the probes run at these
#'   values times \code{n - 1}. The value 0 is the T-Rex with the lasso.
#' @param knee_tolerance Selections below the largest count that still
#'   count as "on the plateau" (default: 1); 0 reads the largest selection,
#'   ties toward the smallest lambda_2.
#' @param verbose Whether the probes print their progress (default:
#'   \code{FALSE}).
#'
#' @return A list with the fields of
#'   \code{\link{trex_tikhonov_lambda2_power_rule}}: \code{lambda2} (the pick,
#'   on the working scale: supplied as \code{lambda_2} it reproduces the
#'   selection), \code{index} (1-based, into \code{grid}), \code{selected}
#'   (the 1-based selection at the pick), \code{grid} (the values the probes
#'   ran at, ascending and de-duplicated), \code{grid_scale} (the factor
#'   between the supplied grid and \code{grid}, \code{n - 1} under
#'   \code{"ZSCORE"} and 1 otherwise), \code{counts} (selection sizes per
#'   grid value), \code{probes} (one list per grid value with
#'   \code{lambda2}, \code{selected}, \code{T_stop} and \code{seconds}),
#'   \code{max_selected}, \code{n_ties}, the boundary flags
#'   \code{at_lower_end} and \code{at_upper_end} (the pick sits on the first
#'   or on the last grid value, and a pick on the last one asks for a wider
#'   grid) and \code{seconds}.
#'
#' @examples
#' \donttest{
#' data(Gauss_data)
#' X <- Gauss_data$X
#' y <- Gauss_data$y
#' prof <- trex_gvs_lambda2_power_rule(X, y, tFDR = 0.1, seed = 1,
#'                                     grid = c(0, 1, 10))
#' prof$lambda2
#' prof$counts
#' }
#' @family T-Rex selectors and controls
#' @export
trex_gvs_lambda2_power_rule <- function(X, y,
                                        tFDR = 0.1,
                                        gvs_control = trex_gvs_control(),
                                        control = NULL,
                                        seed = 0,
                                        grid = trex_default_lambda2_grid(),
                                        knee_tolerance = 1,
                                        verbose = FALSE) {
  .check_power_rule_arguments(seed, grid, knee_tolerance)
  .check_gvs_solver(gvs_control, control)
  .check_gvs_groups_length(gvs_control, X)
  y <- .as_double_vector(y)
  if (inherits(X, "MemoryMappedMatrix")) {
    return(.r_positions(trex_gvs_lambda2_power_rule_mmap_impl(
      get_ptr(X), y, tFDR, gvs_control, control, as.integer(seed),
      as.numeric(grid), as.numeric(knee_tolerance), verbose
    )))
  }
  X <- .as_double_matrix(X)
  .r_positions(trex_gvs_lambda2_power_rule_impl(X, y, tFDR, gvs_control, control,
                                   as.integer(seed), as.numeric(grid),
                                   as.numeric(knee_tolerance), verbose))
}
