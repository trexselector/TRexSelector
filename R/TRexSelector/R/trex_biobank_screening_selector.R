# =============================================================================
# TRexBiobankScreeningSelector.R - Biobank Screening T-Rex Selector
# =============================================================================

#' @name TRexBiobankScreeningSelector
#'
#' @title Biobank Screening with the T-Rex Selector
#'
#' @description Per phenotype the fast screening (T-Rex+Screening), and the T-Rex selector as
#'   fallback when no FDR estimate of the screening lies in the accepted range.
#'
#' @section Decision rule per phenotype:
#' The screening gives two FDR estimates, the one of its ordinary variant and the one of its
#' bootstrap variant. A variant is accepted when its estimate lies in
#' \code{[lower_bound_FDR, upper_bound_FDR]} of \code{\link{trex_biobank_control}()}. Of two
#' accepted variants the one with the larger estimate is taken, and the bootstrap variant on a tie.
#' With no accepted variant the T-Rex selector runs at \code{target_FDR_trex}.
#'
#' Both estimates come from one set of K experiments. The screening runs once, and the bootstrap
#' variant is derived from its coefficients.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Fast Screening of Large-Scale
#'   Genomics Biobanks," in Proc. 22nd IEEE Statist. Signal Process. Workshop (SSP), 2023, pp. 666-670.
#'
#' F. Scheidt, J. Machkour, and M. Muma, "Solving FDR-Controlled Sparse Regression Problems with Five Million
#'   Variables on a Laptop," in Proc. IEEE 9th Int. Workshop Comput. Adv. Multi-Sensor Adapt. Process.
#'   (CAMSAP), 2023, pp. 116-120.
#'
#' J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
#'   High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
#'   Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
#'
#' @examples
#' \donttest{
#' data(Gauss_data)
#' X <- Gauss_data$X
#' # A vector is one phenotype and gives one record, a matrix gives a list of records
#' y <- drop(Gauss_data$y)
#' ctrl <- trex_biobank_control(target_FDR_trex = 0.1)
#' sel <- TRexBiobankScreeningSelector$new(X, y, biobank_control = ctrl)
#' res <- sel$select()
#' res$selected_indices
#' }
#'
#' @family T-Rex selectors and controls
#' @export
TRexBiobankScreeningSelector <- R6::R6Class("TRexBiobankScreeningSelector",
  public = list(
    #' @description Create a new TRexBiobankScreeningSelector object.
    #' @param X Feature matrix (n x p) in memory, or an \code{mmap_matrix},
    #'   the handle of a memory-mapped file. Both are bound without a copy.
    #'   The inner selectors standardize a memory-mapped design in place per
    #'   call and restore it, which needs \code{"readwrite"} mode. Under
    #'   \code{scaling_mode = "NONE"} they write nothing and
    #'   \code{"readonly"} mode suffices.
    #' @param Y Response vector (1D) or matrix (2D). The matrix is in memory
    #'   or an \code{mmap_matrix}, the handle of a memory-mapped file, and is
    #'   bound without a copy. It is only read, so \code{"readonly"} mode
    #'   suffices for it.
    #' @param tFDR Target FDR level of the fallback T-Rex selector. \code{NULL} (default)
    #'   takes \code{target_FDR_trex} of \code{biobank_control}. A supplied value replaces it,
    #'   as in TRexSelector 2.0.1.
    #' @param seed Random seed. A negative value (default: -1) gives a run that is not
    #'   reproducible. Phenotype \code{k} (the \code{k}-th column of
    #'   \code{Y}) runs with \code{seed + k - 1}, so the phenotypes draw independent dummies and
    #'   a single phenotype is reproduced by a single-phenotype run with that seed.
    #' @param verbose Whether to print progress (default: FALSE).
    #' @param biobank_control A biobank control list from \code{\link{trex_biobank_control}()}
    #'   (default: \code{trex_biobank_control()}).
    #' @param screen_control A screen control list from \code{\link{trex_screen_control}()}
    #'   (default: \code{trex_screen_control()}). Its \code{use_bootstrap_CI} has no effect
    #'   here: the biobank screening computes the ordinary and the bootstrap variant for every
    #'   phenotype.
    #' @param control A control list from \code{\link{trex_control}()} (default:
    #'   \code{trex_control()}). The screening takes the values \code{"STANDARD"} (default),
    #'   \code{"SEEDED"}, \code{"PERMUTATION"} and \code{"PERMUTATION_SEEDED"} of
    #'   \code{lloop_strategy} and refuses \code{"HCONCAT"} and \code{"SKIPL"}, as
    #'   \code{\link{TRexScreeningSelector}} does.
    initialize = function(X, Y,
                          tFDR = NULL,
                          seed = -1,
                          verbose = FALSE,
                          biobank_control = trex_biobank_control(),
                          screen_control = trex_screen_control(),
                          control = trex_control()) {

      # A supplied tFDR replaces the target of the control, as in 2.0.1.
      biobank_control_list <- biobank_control
      if (!is.null(tFDR)) biobank_control_list$target_FDR_trex <- tFDR

      y_mmap <- inherits(Y, "mmap_matrix")
      if (!is.numeric(Y) && !y_mmap) {
        stop("Response evaluation Y must be a numeric vector, a numeric matrix or an mmap_matrix.")
      }

      # A memory-mapped X is bound directly. The inner selectors normalize
      # it in place per call and restore it ("readwrite"), except under
      # scaling_mode = "NONE", where they only read it. An R matrix is bound
      # zero-copy in double storage.
      is_mmap <- inherits(X, "mmap_matrix")
      if (!is_mmap) X <- .as_double_matrix(X)
      n_x <- if (is_mmap) dim(X)[1] else nrow(X)

      if (is.vector(Y) || (is.array(Y) && length(dim(Y)) == 1)) {
        if (n_x != length(Y)) stop("Dimension mismatch: nrow(X) != length(Y)")
        Y <- .as_double_vector(Y)
        private$refs <- list(X = X, Y = Y)
        private$is_multi <- FALSE
        private$ptr <- if (is_mmap) {
          .r_positions(trex_biobank_screening_1d_mmap_create(
            X, Y, biobank_control_list, screen_control, control, seed, verbose))
        } else {
          .r_positions(trex_biobank_screening_1d_create(
            X, Y, biobank_control_list, screen_control, control, seed, verbose))
        }
      } else if (is.matrix(Y) || y_mmap) {
        n_y <- if (y_mmap) dim(Y)[1] else nrow(Y)
        if (n_x != n_y) stop("Dimension mismatch: nrow(X) != nrow(Y)")
        if (!y_mmap) Y <- .as_double_matrix(Y)
        private$refs <- list(X = X, Y = Y)
        private$is_multi <- TRUE
        private$ptr <- if (is_mmap) {
          .r_positions(trex_biobank_screening_2d_mmap_create(
            X, Y, biobank_control_list, screen_control, control, seed, verbose))
        } else {
          .r_positions(trex_biobank_screening_2d_create(
            X, Y, biobank_control_list, screen_control, control, seed, verbose))
        }
      } else {
        stop("Response evaluation Y must be a numeric vector, a numeric matrix or an mmap_matrix.")
      }
    },

    #' @description Run the Biobank Screening algorithm.
    #'
    #' @return For a vector \code{Y}, a single per-phenotype record: a named list
    #'   with \code{phenotype_index}, \code{selected_indices}, \code{estimated_FDR},
    #'   \code{method_used} (\code{"Screen-TRex (ordinary)"}, \code{"Screen-TRex (bootstrap-CI)"}
    #'   or \code{"T-Rex (fallback)"}), \code{estimated_FDR_screen_ordinary},
    #'   \code{estimated_FDR_screen_bootstrap}, \code{selected_indices_screen_ordinary},
    #'   \code{selected_indices_screen_bootstrap}, and \code{used_fallback_trex}
    #'   (indices are 1-based). \code{estimated_FDR} is the estimate of the accepted screening
    #'   variant, or, for a fallback result, the target level \code{target_FDR_trex}, at which the T-Rex
    #'   selector controls the FDR. An estimate can equal the target, so
    #'   \code{used_fallback_trex} tells the two cases apart. The two screening estimates and
    #'   the two screening index vectors are stored for every phenotype, a fallback result
    #'   included. For a matrix \code{Y}, a plain list of such records,
    #'   one per column, so \code{res[[i]]} has the same shape as the vector case.
    #'   This mirrors the Python binding, which returns one result object (vector)
    #'   or a list of them (matrix). A tidy summary table is one line away, e.g.
    #'   \code{do.call(rbind, lapply(res, function(r)}
    #'   \code{as.data.frame(r[c("phenotype_index", "estimated_FDR", "method_used", "used_fallback_trex")])))}.
    select = function() {
      if (private$is_multi) {
        return(.r_positions(trex_biobank_screening_screen_phenotypes(private$ptr)))
      } else {
        return(.r_positions(trex_biobank_screening_screen_phenotype(private$ptr)))
      }
    }
  ),
  private = list(
    ptr = NULL,
    refs = NULL,
    is_multi = FALSE
  )
)
