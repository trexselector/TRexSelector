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
    #'   \code{trex_gvs_control()}). For \code{gvs_type = "EN"} the solver is
    #'   derived from \code{en_solver} (TENET, TENET_AUG, or TCDENET) and the
    #'   \code{solver} field in \code{control} is ignored (a warning is issued
    #'   if it is set to a non-matching solver). For \code{gvs_type = "IEN"}
    #'   the \code{solver} field selects the IEN solver: "TIENET"
    #'   (default, native pathwise), "TIENET_AUG", or "TCDIENET" (not finally
    #'   developed, see \code{\link{trex_gvs_control}}). "TCDENET" and
    #'   "TCDIENET" are deprecated and will be removed in TRexSelector 2.1.0:
    #'   the selector raises a deprecation warning for them.
    #' @param control A control list from \code{\link{trex_control}()} (default:
    #'   \code{trex_control()}). For \code{gvs_type = "EN"} the \code{solver}
    #'   field is overridden; for "IEN" it selects the IEN-family solver.
    initialize = function(X, y,
                          tFDR = 0.1,
                          seed = -1,
                          verbose = TRUE,
                          gvs_control = trex_gvs_control(),
                          control = trex_control()) {

      # Solver selection mirrors the C++ design:
      #   EN  -> derived from en_solver (TENET / TENET_AUG / TCDENET);
      #          control$solver is ignored (warn when it was explicitly set
      #          to a non-matching value; "TLARS" is the trex_control()
      #          default and cannot be distinguished from an explicit
      #          choice, so it never warns).
      #   IEN -> control$solver IS the IEN-solver choice: "TIENET" (native
      #          pathwise, default), "TIENET_AUG" (row-augmented), or
      #          "TCDIENET" (coordinate descent). The default "TLARS" maps
      #          to "TIENET"; any other value warns and falls back.
      if (gvs_control$gvs_type == "IEN") {
        ien_family <- c("TIENET", "TIENET_AUG", "TCDIENET")
        if (!(control$solver %in% c("TLARS", ien_family))) {
          warning("control$solver = \"", control$solver,
                  "\" is not an IEN-family solver; gvs_type = \"IEN\" ",
                  "accepts \"TIENET\", \"TIENET_AUG\", or \"TCDIENET\". ",
                  "Using \"TIENET\".")
          control$solver <- "TLARS"  # glue defaults it to TIENET
        }
        .warn_cd_solver_name(control$solver)
      } else {
        derived_solver <- if (identical(gvs_control$en_solver, "TENET_AUG")) {
          "TENET_AUG"
        } else if (identical(gvs_control$en_solver, "TCDENET")) {
          "TCDENET"
        } else {
          "TENET"
        }
        if (!(control$solver %in% c("TLARS", derived_solver))) {
          warning("control$solver = \"", control$solver,
                  "\" is ignored by TRexGVSSelector: gvs_type = \"EN\" ",
                  "derives solver \"", derived_solver,
                  "\" from en_solver.")
        }
        # The C++ glue derives solver_type from en_solver regardless.
        control$solver <- "TLARS"
        .warn_cd_solver_name(gvs_control$en_solver, "en_solver")
      }

      # lambda_2 == 0 is the degenerate no-ridge case, NOT automatic CV
      # (the sentinel for that is lambda_2 < 0). A user passing exactly 0
      # in the belief it triggers CV silently gets pure T-LASSO.
      if (isTRUE(gvs_control$lambda_2 == 0)) {
        warning("gvs_control$lambda_2 = 0 is the degenerate no-ridge case ",
                "(pure T-LASSO), not automatic cross-validation. ",
                "Use lambda_2 < 0 (e.g. -1) to auto-determine lambda_2 via CV.")
      }

      # One group entry per variable, checked here with the name of the R argument.
      if (!is.null(gvs_control$groups) && length(gvs_control$groups) != ncol(X)) {
        stop("gvs_control$groups must have one entry per variable (p = ", ncol(X),
             "). Got: ", length(gvs_control$groups), ".")
      }

      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_gvs_mmap_create(
          get_ptr(X), y, tFDR, gvs_control, control, seed, verbose
        )
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_gvs_create(
          X, y, tFDR, gvs_control, control, seed, verbose
        )
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
