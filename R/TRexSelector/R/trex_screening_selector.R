# =========================================================================
# trex_screening_selector.R - TRexScreeningSelector class definition
# =========================================================================

#' @name TRexScreeningSelector
#'
#' @title T-Rex Screening Selector (T-Rex+Screening)
#'
#' @description Fast FDR-controlled variable screening in one pass of the random experiments, until
#'   one dummy is included (T = 1), with as many dummies as variables (L = p).
#'
#' @section Fields inherited from TRexSelector:
#' The screening has no target FDR and no calibration, so the fields of \code{TRexSelector} hold fixed
#' values after a run: \code{T_stop} is 1, \code{L} is 1, \code{num_dummies} is \code{p},
#' \code{v_thresh} is 0.5 and \code{tFDR} is 0. \code{r_mat} is a 1 x 1 matrix with the number of
#' selected variables, \code{phi_prime} equals the relative occurrences, and \code{fdp_hat_mat} is a
#' 1 x 1 matrix with 0, because no FDP is estimated per voting level. The FDR estimate of the run is
#' \code{estimated_FDR}.
#'
#' @section L-loop strategies:
#' The screening takes the strategies \code{"STANDARD"} (default), \code{"SEEDED"},
#' \code{"PERMUTATION"} and \code{"PERMUTATION_SEEDED"} of \code{\link{trex_control}()}.
#' \code{"SEEDED"} gives the results of \code{"STANDARD"} and keeps one dummy matrix in memory at a
#' time. \code{"PERMUTATION_SEEDED"} gives the results of \code{"PERMUTATION"} without storing the base
#' dummy matrix. \code{"HCONCAT"} and \code{"SKIPL"} are refused, because they describe the growth of
#' the dummy matrix over the L-loop, and the screening is one pass at L = p.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Fast Screening of Large-Scale
#'   Genomics Biobanks," in Proc. 22nd IEEE Statist. Signal Process. Workshop (SSP), 2023, pp. 666-670.
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
#' sel <- TRexScreeningSelector$new(X, y, verbose = FALSE)
#' sel$select()
#' sel$selected_indices
#' }
#'
#' @family T-Rex selectors and controls
#' @export
TRexScreeningSelector <- R6::R6Class("TRexScreeningSelector",
  inherit = TRexSelector,
  public = list(
    #' @description Create a new TRexScreeningSelector object.
    #'
    #' @param X Feature matrix in memory, or an \code{mmap_matrix}, the handle
    #'   of a memory-mapped file. Both are bound without a copy. The selector
    #'   standardizes a memory-mapped design in place and restores it, which
    #'   needs \code{"readwrite"} mode. Under \code{scaling_mode = "NONE"} it
    #'   writes nothing and \code{"readonly"} mode suffices.
    #' @param y Response vector.
    #' @param seed Random seed. A negative value (default: -1) gives a run that is not
    #'   reproducible.
    #' @param verbose Whether to print progress (default: TRUE).
    #' @param screen_control A screen control list from \code{\link{trex_screen_control}()}
    #'   (default: \code{trex_screen_control()}).
    #' @param control A control list from \code{\link{trex_control}()} (default:
    #'   \code{trex_control()}). See the section on the L-loop strategies for the values of
    #'   \code{lloop_strategy} that the screening takes.
    initialize = function(X, y,
                          seed = -1,
                          verbose = TRUE,
                          screen_control = trex_screen_control(),
                          control = trex_control()) {

      .warn_cd_solver_name(control$solver)

      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_screening_mmap_create(
          get_ptr(X), y, screen_control, control, seed, verbose
        )
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_screening_create(
          X, y, screen_control, control, seed, verbose
        )
      }
    },

    #' @description Run the Screening algorithm.
    select = function() {
      trex_screening_select(private$ptr)
      invisible(self)
    }
  ),
  active = list(
    #' @field estimated_FDR The estimated FDR, one over the number of selected variables. An
    #'   empty selection reports 1.
    estimated_FDR = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_estimated_fdr(private$ptr)
    },

    #' @field used_bootstrap Whether the confidence-based screening (bootstrap) made the
    #'   selection.
    used_bootstrap = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_used_bootstrap(private$ptr)
    },

    #' @field confidence_level Confidence level that the grid search of the bootstrap picked.
    #'   0 without the bootstrap.
    confidence_level = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_confidence_level(private$ptr)
    },

    #' @field estimated_correlation Correlation coefficient of the dependency-aware variants:
    #'   the supplied value, or the estimate from the data. 0 for "TREX".
    estimated_correlation = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_estimated_correlation(private$ptr)
    },

    #' @field beta_mat Coefficients of the variables, averaged over the K experiments (one
    #'   per variable, a p x 1 matrix).
    beta_mat = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_beta_mat(private$ptr)
    },

    #' @field dummy_betas Non-zero coefficients of the dummies, collected over the K
    #'   experiments. They are the sample of the bootstrap.
    dummy_betas = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_screening_get_dummy_betas(private$ptr)
    }
  )
)
