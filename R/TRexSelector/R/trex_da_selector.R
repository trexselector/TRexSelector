#' @name TRexDASelector
#'
#' @title Dependency-Aware T-Rex Selector (T-Rex+DA)
#'
#' @description FDR-controlled variable selection for dependent variables.
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
#' @section The reference threshold:
#' The methods \code{"BT"}, \code{"NN"} and prior groups evaluate every step on the whole grid of
#' correlation thresholds (\code{rho_grid}). The L-loop and the T-loop read the FDP estimate at one of
#' them, the reference threshold, which sits at three quarters of the grid: element 0.75 times its
#' length, a half rounded up (element 15 of the default 20). The 1.x versions of the package round a
#' half to the even number and take one element lower for the grid lengths 6, 14, 22, 30 and so on.
#' The selection then searches all cells of stopping time, voting level and threshold, and
#' \code{rho_thresh} is the threshold of the cell it takes. The inherited fields \code{fdp_hat_mat},
#' \code{phi_mat} and \code{phi_prime} therefore hold the values at the reference threshold, while
#' \code{r_mat}, \code{v_thresh} and the selection belong to \code{rho_thresh}. The values at every
#' threshold are in \code{FDP_hat_array_BT}, \code{Phi_array_BT} and \code{R_array_BT}.
#'
#' @section FDR limitation of BT, NN and prior groups:
#' At a very high correlation within the groups (about 0.9 and above) and with no independent
#' ("white") predictors, the equiangular solver \code{"TLARS"} can exceed the target FDR, because
#' the true signal and its collinear shadows become statistically exchangeable. This is inherent to
#' the method, and it also appears in the 1.x versions of the package. Mitigate with the default
#' \code{bt_selection_mode = "FeasibleOnly"} and by preferring a greedy solver (\code{"TOMP"}) in
#' designs with high correlation.
#'
#' @section Greedy solvers and the exchangeable-tie band:
#' A greedy solver (\code{"TOMP"}, \code{"TOOLS"}, \code{"TGP"}, \code{"TACGP"}, \code{"TMP"},
#' \code{"TNCGMP"}, \code{"TAFS"}) picks the top candidate of a correlated cluster deterministically,
#' which removes the spread of occurrences within the cluster that the penalty relies on. The selector
#' therefore runs a greedy solver with the exchangeable-tie band at 0.25 when
#' \code{trex_control(exch_tie_alpha = )} is left at its default of -1. A value of 0 switches the
#' band off, and a positive value sets its width. The path solvers ignore the band.
#'
#' Evidence: a Monte Carlo study on AR(1) designs (1,000 runs per setting, target FDR 0.2, method
#' \code{"AR1"}, November 2026). At a correlation of 0.9 each of the seven greedy solvers had an
#' FDR of 0.21 to 0.27 without the band and of 0.13 to 0.17 with it, at a true positive rate lower
#' by 0.07 to 0.13. At correlations up to 0.8 the target held without the band as well. The study
#' shows the control up to a correlation of 0.9. At 0.95 the target no longer held (500 runs per
#' setting), with an FDR of up to 0.23 for the path solver \code{"TLARS"} and of up to 0.30 for
#' \code{"TOMP"} with the band.
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
#' # Dependent variables: an AR(1) process across the columns, five of them active
#' set.seed(321)
#' n <- 100
#' p <- 100
#' rho <- 0.7
#' X <- matrix(0, n, p)
#' X[, 1] <- rnorm(n)
#' for (j in 2:p) X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
#' true_support <- c(10, 30, 50, 70, 90)
#' beta <- numeric(p)
#' beta[true_support] <- 2
#' y <- as.numeric(X %*% beta + rnorm(n))
#'
#' sel <- TRexDASelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
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
    #' @param seed Random seed. A negative value (default: -1) gives a run that is not
    #'   reproducible.
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

      # No need to call super$initialize() because the Cpp pointer constructor is different.
      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- .r_positions(trex_da_mmap_create(
          get_ptr(X), y, tFDR, da_control, control, seed, verbose
        ))
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- .r_positions(trex_da_create(
          X, y, tFDR, da_control, control, seed, verbose
        ))
      }
    },

    #' @description Run the Dependency-Aware TRexSelector algorithm.
    #'   A selector runs once: a second call is an error, see \code{\link{TRexSelector}}.
    select = function() {
      trex_da_select(private$ptr)
      invisible(self)
    }
  ),
  active = list(
    #' @field rho_thresh Correlation threshold of the selected cell of the calibration grid
    #'   (BT, NN and prior groups). -2 for AR1 and EQUI and when no cell meets the target.
    rho_thresh = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_rho_thresh(private$ptr)
    },
    #' @field rho_grid Grid of correlation thresholds of BT, NN and prior groups:
    #'   \code{hc_grid_length} ascending values. BT and prior groups take them from the merge
    #'   heights of the dendrogram (1 - height) and end at 1, NN spaces them equally on
    #'   \code{[0, 1]}. Empty for AR1 and EQUI.
    rho_grid = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_rho_grid(private$ptr)
    },

    #' @field da_method_used The dependency-aware method of the last \code{select()}:
    #'   \code{"BT"}, \code{"NN"}, \code{"AR1"}, \code{"EQUI"}, or \code{"PRIOR_GROUPS"} when
    #'   prior groups were supplied.
    da_method_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_da_get_method_string(private$ptr)
    },

    #' @field FDP_hat_array_BT FDP estimates of the methods BT, NN and prior groups: an array
    #'   with one row per stopping time (\code{T_stop}), one column per voting level
    #'   (\code{voting_grid}) and one slice per correlation threshold (\code{rho_grid}).
    #'   Empty for AR1 and EQUI.
    FDP_hat_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array (T_stop x V x rho levels) built once in C++; an empty
      # result is a length-0 numeric vector.
      trex_da_get_fdp_hat_array_bt(private$ptr)
    },

    #' @field Phi_array_BT Dependency-aware relative occurrences of the methods BT, NN and
    #'   prior groups: an array with one row per stopping time (\code{T_stop}), one column
    #'   per variable and one slice per correlation threshold (\code{rho_grid}). Empty for
    #'   AR1 and EQUI.
    Phi_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array (T_stop x p x rho levels) built once in C++; an empty
      # result is a length-0 numeric vector.
      trex_da_get_phi_array_bt(private$ptr)
    },

    #' @field R_array_BT Numbers of selected variables of the methods BT, NN and prior
    #'   groups: an array with one row per stopping time that can be selected
    #'   (\code{T_stop - 1}), one column per voting level (\code{voting_grid}) and one
    #'   slice per correlation threshold (\code{rho_grid}). Empty for AR1 and EQUI.
    R_array_BT = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      # A 3-d array ((T_stop - 1) x V x rho levels) built once in C++; an
      # empty result is a length-0 numeric vector.
      trex_da_get_r_array_bt(private$ptr)
    },

    #' @field cor_coef Correlation coefficient of the methods AR1 and EQUI: the supplied
    #'   value, or the estimate from the data when -2 requested it. The other methods do not
    #'   use it and return the value of the control.
    cor_coef = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      trex_da_get_cor_coef(private$ptr)
    }
  )
)
