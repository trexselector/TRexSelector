#' @name TRexDASelector
#'
#' @title Dependency-Aware T-Rex Selector (TRexDA)
#'
#' @description FDR-controlled variable selection algorithm for dependency scenarios.
#'
#' @section The penalty of the last step:
#' The dependency-aware relative occurrences of variable \eqn{j} at step \eqn{t} of the T-loop are
#' \deqn{\Phi^{DA}_t(j) = \Psi_T(j) \, \Phi_t(j), \qquad
#'   \Psi_T(j) = \frac{1}{2 - \min_{j'} |\Phi_T(j) - \Phi_T(j')|}, \qquad t = 1, \dots, T,}{
#'   Phi^DA_t(j) = Psi_T(j) * Phi_t(j),  Psi_T(j) = 1 / (2 - min_j' |Phi_T(j) - Phi_T(j')|),  t = 1, ..., T,}
#' where \eqn{T} is the last step of the evaluated path and \eqn{j'} runs over the variables that
#' \eqn{j} is compared with (its group, its window or all other variables). A variable without a group
#' gets \eqn{\Psi_T(j) = 1/2}.
#'
#' This is a deliberate deviation from the papers below. Read literally, they recompute the penalty at
#' every step, \eqn{\Psi_t(j)} from \eqn{\Phi_t}. The selector takes the penalty of the last step for
#' all steps, with all five methods. The reason is the FDP estimator. It sums over the increments
#' \eqn{\Phi^{DA}_t(j) - \Phi^{DA}_{t-1}(j)} and counts expected false discoveries only while no
#' increment is negative. One penalty per variable keeps the increments non-negative. A penalty per
#' step does not, because the penalty of a variable changes from step to step, and the FDP estimate
#' then becomes too small.
#'
#' A Monte Carlo comparison of the two constructions (1,000 runs per setting on block, AR(1) and
#' independent designs) supports the choice. With the method \code{"NN"} the per-step penalty exceeded
#' a target FDR of 0.1, with 0.15 on an independent design and 0.13 on a block design, where the
#' last-step penalty gave 0.00 and 0.03. No setting had the last-step penalty above the target and the
#' per-step penalty below it. The per-step penalty selected more variables only where its estimate was
#' invalid, and on AR(1) designs the last-step penalty had the higher true positive rate.
#'
#' The methods \code{"AR1"} and \code{"EQUI"} used the per-step penalty up to TRexSelector 2.0.1, so
#' their selections can differ slightly from earlier versions. \code{"BT"}, \code{"NN"} and prior
#' groups have always used the penalty of the last step.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' J. Machkour, M. Muma, and D. P. Palomar, "High-Dimensional False Discovery Rate Control for Dependent
#'   Variables," Signal Process., vol. 234, Art. no. 109990, 2025, doi: 10.1016/j.sigpro.2025.109990.
#'
#' J. Machkour, D. P. Palomar, and M. Muma, "FDR-Controlled Portfolio Optimization for Sparse Financial Index
#'   Tracking," in Proc. IEEE Int. Conf. Acoust., Speech, Signal Process. (ICASSP), 2025, pp. 1-5,
#'   doi: 10.1109/ICASSP49660.2025.10890629.
#'
#' J. Machkour, M. Muma, and D. P. Palomar, "FDR-Controlled Sparse Index Tracking with Autoregressive Stock
#'   Dependency Models," in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2024, pp. 2662-2666.
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
#' sel <- TRexDASelector$new(X, y, tFDR = 0.1, verbose = FALSE,
#'                           da_control = trex_da_control(da_method = "BT"))
#' sel$select()
#' sel$selected_indices
#' }
#'
#' @family T-Rex selectors and controls
#' @export
TRexDASelector <- R6::R6Class("TRexDASelector",
  inherit = TRexSelector,
  public = list(
    #' @description Create a new TRexDASelector object.
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
    #' @param da_control A DA control list from \code{\link{trex_da_control}()} (default:
    #'   \code{trex_da_control()}).
    #' @param control A control list from \code{\link{trex_control}()} (default:
    #'   \code{trex_control()}).
    initialize = function(X, y,
                          tFDR = 0.1,
                          seed = -1,
                          verbose = TRUE,
                          da_control = trex_da_control(),
                          control = trex_control()) {

      .warn_cd_solver_name(control$solver)

      # No need to call super$initialize() because the Cpp pointer constructor is different.
      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_da_mmap_create(
          get_ptr(X), y, tFDR, da_control, control, seed, verbose
        )
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_da_create(
          X, y, tFDR, da_control, control, seed, verbose
        )
      }
    },

    #' @description Run the Dependency-Aware TRexSelector algorithm.
    select = function() {
      trex_da_select(private$ptr)
      invisible(self)
    }
  ),
  active = list(
    #' @field rho_thresh Optimal correlation threshold selected by the DA process.
    rho_thresh = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_rho_thresh(private$ptr)
    },
    #' @field rho_grid The grid of correlation thresholds evaluated.
    rho_grid = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_rho_grid(private$ptr)
    },

    #' @field da_method_used The integer code (or name mapped later) of the DA method used.
    da_method_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_method_string(private$ptr)
    },

    #' @field FDP_hat_array_BT A 3D array of FDP estimates for BT method:
    #' (Variables x Thresholds x K).
    FDP_hat_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array (T_stop x p x R) built once in C++; an empty result
      # is a length-0 numeric vector.
      trex_da_get_fdp_hat_array_bt(private$ptr)
    },

    #' @field Phi_array_BT A 3D array of variables selected for BT method.
    Phi_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array (T_stop x p x R) built once in C++; an empty result
      # is a length-0 numeric vector.
      trex_da_get_phi_array_bt(private$ptr)
    },

    #' @field R_array_BT A 3D array mapping original to clustered features for BT method.
    R_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array (T_stop x p x R) built once in C++; an empty result
      # is a length-0 numeric vector.
      trex_da_get_r_array_bt(private$ptr)
    },

    #' @field cor_coef The correlation coefficient of the dependency-aware
    #'   dummy generation (AR1 / EQUI methods).
    cor_coef = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      trex_da_get_cor_coef(private$ptr)
    }
  )
)
