# =============================================================================
# trex_spca_selector.R - TRexSPCASelector class and control parameter constructor
# =============================================================================

#' @name trex_spca_control
#'
#' @title T-Rex+SPCA Control Parameters
#'
#' @description Build the control list of the T-Rex sparse PCA selector (T-Rex+SPCA),
#'   \code{\link{TRexSPCASelector}}. The list holds the parameters below. Every other parameter
#'   of the T-Rex+GVS run of each component (the clustering, the cross-validation grid, the
#'   L-loop strategy, the scaling mode, the threads) runs at its default, see
#'   \code{\link{trex_gvs_control}()} and \code{\link{trex_control}()}.
#'
#' @param mode Assembly of the sparse loading vectors. \code{"ActiveSet"} (default) is the
#'   ridge regression of the component on the selected variables, normalized to unit length
#'   (the T-Rex PCA of the paper). \code{"Thresholded"} keeps the largest ordinary PCA
#'   loadings in absolute value, as many as variables were selected, normalized to unit
#'   length (the T-Rex Thresholded PCA of the paper).
#' @param lambda2_ridge_loadings Ridge penalty of the loading assembly in
#'   \code{"ActiveSet"} mode (default: \code{1e-6}). Not used in \code{"Thresholded"} mode.
#' @param gvs_type Penalty of the T-Rex+GVS runs. Only \code{"EN"} (default, elastic net) is
#'   accepted: the method is defined on the elastic net, and \code{"IEN"} is refused.
#' @param lambda_2 Ridge part of the elastic net of the T-Rex+GVS runs. A negative value
#'   (default: \code{-1}) lets the rule \code{lambda2_method} choose it, \code{0} is the
#'   T-LASSO without the ridge part, and a positive value is used as it is. A supplied
#'   value refers to the scale of the standardized columns.
#' @param lambda2_method Rule that chooses \code{lambda_2} when \code{lambda_2 < 0}:
#'   \code{"CV_1SE_CCD"} (default), \code{"CV_MIN_CCD"}, \code{"CV_1SE_SVD"} or
#'   \code{"CV_MIN_SVD"}. The rules of the informed elastic net are refused.
#' @param en_solver Elastic net solver of the T-Rex+GVS runs. \code{"TENET_AUG"} (default) is
#'   the elastic net as a LASSO on row-augmented data, the form of the paper.
#'   \code{"TENET"} is the elastic net on the Gram matrix. Both solve the same elastic net
#'   problem.
#' @param tenet_aug_use_lars Logical. For \code{en_solver = "TENET_AUG"}: run the augmented
#'   system with a pure-LARS inner path that never drops variables (\code{TRUE}) instead of
#'   the LASSO inner path that can (\code{FALSE}, default). Ignored for \code{"TENET"}.
#' @param K Number of random experiments of each T-Rex+GVS run (default: 20).
#' @param max_dummy_multiplier Cap of the dummy multiplier of each T-Rex+GVS run: the L-loop
#'   adds dummies up to this many times the number of variables (default: 10).
#'
#' @return A named list for passing to \code{TRexSPCASelector$new()}.
#'
#' @examples
#' trex_spca_control()
#' trex_spca_control(mode = "Thresholded", K = 10)
#' trex_spca_control(lambda_2 = 0.01)
#'
#' @family T-Rex selectors and controls
#' @export
trex_spca_control <- function(mode                   = "ActiveSet",
                              lambda2_ridge_loadings  = 1e-6,
                              gvs_type               = "EN",
                              lambda_2               = -1.0,
                              lambda2_method          = "CV_1SE_CCD",
                              en_solver              = "TENET_AUG",
                              tenet_aug_use_lars     = FALSE,
                              K                      = 20L,
                              max_dummy_multiplier   = 10L) {
  mode          <- match.arg(mode,          c("ActiveSet", "Thresholded"))
  # T-Rex+SPCA is defined on the elastic net, so an IEN request is an error
  if (!identical(gvs_type, "EN"))
    stop("trex_spca_control: gvs_type must be \"EN\". T-Rex+SPCA is ",
         "defined on the elastic net, and the informed elastic net (IEN) is ",
         "not supported. Got: ", format(gvs_type))
  lambda2_method <- match.arg(lambda2_method,
                              c("CV_1SE_CCD", "CV_MIN_CCD",
                                "CV_1SE_SVD", "CV_MIN_SVD"))
  en_solver     <- match.arg(en_solver, c("TENET_AUG", "TENET"))
  stopifnot(
    is.numeric(lambda2_ridge_loadings), lambda2_ridge_loadings >= 0,
    is.numeric(lambda_2),
    is.numeric(K),                      K >= 1L,
    is.numeric(max_dummy_multiplier),   max_dummy_multiplier >= 1L
  )
  list(
    mode                  = mode,
    lambda2_ridge_loadings = as.numeric(lambda2_ridge_loadings),
    gvs_type              = gvs_type,
    lambda_2              = as.numeric(lambda_2),
    lambda2_method        = lambda2_method,
    en_solver             = en_solver,
    tenet_aug_use_lars    = as.logical(tenet_aug_use_lars),
    K                     = as.integer(K),
    max_dummy_multiplier  = as.integer(max_dummy_multiplier)
  )
}


#' @name TRexSPCASelector
#'
#' @title T-Rex Sparse PCA Selector (T-Rex+SPCA)
#'
#' @description Sparse principal components whose non-zero loadings are selected with FDR
#'   control. Each ordinary principal component is the response of a T-Rex+GVS run with the
#'   elastic net.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item \code{X} is centered and its first \code{M} ordinary principal components are
#'     computed (a PCA of the covariance, the columns are not scaled).
#'   \item For each component, a T-Rex+GVS run with the elastic net and the component as
#'     response selects its variables at the target FDR.
#'   \item The sparse loading vector is assembled from the selected set (see \code{mode} of
#'     \code{\link{trex_spca_control}()}), and the scores of the sparse component are the
#'     centered \code{X} times the loading vector.
#'   \item The explained variance is adjusted for components that are not orthogonal, by
#'     the QR decomposition of the scores.
#' }
#'
#' A component without a selected variable has a zero column in \code{V} and \code{Z} and
#' the explained variance 0.
#'
#' @importFrom R6 R6Class
#'
#' @references
#' J. Machkour, A. Breloy, M. Muma, D. P. Palomar, and F. Pascal, "Sparse PCA with False Discovery Rate
#'   Controlled Variable Selection," in Proc. IEEE 49th Int. Conf. Acoust., Speech, Signal Process. (ICASSP),
#'   2024, pp. 9716-9720.
#'
#' J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
#'   High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
#'   Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100; p <- 30
#' X <- matrix(rnorm(n * p), n, p)
#' sel <- TRexSPCASelector$new(X, M = 2, tFDR = 0.2,
#'                              control = trex_spca_control(K = 5))
#' sel$select()
#' sel$adjusted_ev
#' sel$active_sets
#' }
#'
#' @family T-Rex selectors and controls
#' @export
TRexSPCASelector <- R6::R6Class("TRexSPCASelector",
  private = list(
    X_ref   = NULL,
    M_      = NULL,
    tFDR_   = NULL,
    seed_   = NULL,
    verbose_ = NULL,
    ctrl    = NULL,
    result_ = NULL
  ),
  public = list(

    #' @description Create a new TRexSPCASelector.
    #'
    #' @param X Design matrix (n x p) in memory, or an \code{mmap_matrix} opened in
    #'   \code{"readwrite"} mode. Both are bound without a copy. \code{X} is centered in place
    #'   during \code{select()} and restored on return.
    #' @param M Number of sparse principal components, 1 to \code{min(n, p)}.
    #' @param tFDR Target FDR level of each component, in \code{[0, 1]}.
    #' @param control A control list from \code{\link{trex_spca_control}()}
    #'   (default: \code{trex_spca_control()}).
    #' @param seed Random seed. A negative value (default: -1) gives a run that is not
    #'   reproducible. Each component runs with its own seed derived from it.
    #' @param verbose Whether to print one line per component (default: \code{FALSE}). The
    #'   T-Rex+GVS runs stay quiet.
    initialize = function(X, M, tFDR,
                          control = trex_spca_control(),
                          seed    = -1L,
                          verbose = FALSE) {
      # A memory-mapped X is centered in the file during select() and
      # restored ("readwrite"); an R matrix is bound zero-copy.
      is_mmap <- inherits(X, "mmap_matrix")
      if (!is_mmap && !is.matrix(X)) X <- as.matrix(X)
      stopifnot(
        is_mmap || is.numeric(X),
        is.numeric(M),    length(M) == 1L,
        is.numeric(tFDR), length(tFDR) == 1L,
        is.numeric(control$lambda2_ridge_loadings),
        is.logical(verbose)
      )
      # The bounds of the C++ constructor, checked here so that the error comes
      # at construction: the C++ object is built in select().
      if (is.na(M) || M < 1 || M > min(dim(X))) {
        stop("M must be a whole number from 1 to min(n, p) = ", min(dim(X)), ". Got: ", M)
      }
      if (is.na(tFDR) || tFDR < 0 || tFDR > 1) {
        stop("tFDR must be in [0, 1]. Got: ", tFDR)
      }
      # lambda_2 == 0 is the degenerate no-ridge case, NOT automatic CV
      # (the sentinel for that is lambda_2 < 0). Warn so a user who passed 0
      # expecting CV is not silently handed pure T-LASSO.
      if (isTRUE(control$lambda_2 == 0)) {
        warning("control$lambda_2 = 0 is the degenerate no-ridge case ",
                "(pure T-LASSO), not automatic cross-validation. ",
                "Use lambda_2 < 0 (e.g. -1) to auto-determine lambda_2 via CV.")
      }
      # Double storage so the zero-copy select() view aliases this buffer
      # instead of a silent integer->double copy (see .as_double_matrix).
      private$X_ref    <- if (is_mmap) X else .as_double_matrix(X)
      private$M_       <- as.integer(M)
      private$tFDR_    <- as.numeric(tFDR)
      private$seed_    <- as.integer(seed)
      private$verbose_ <- as.logical(verbose)
      private$ctrl     <- control
    },

    #' @description Run the T-Rex sparse PCA selection.
    #'
    #' @return Invisible \code{self} for method chaining. The results are in the fields
    #'   \code{Z}, \code{V}, \code{active_sets}, \code{adjusted_ev} and \code{cumulative_ev}.
    select = function() {
      selector <- if (inherits(private$X_ref, "mmap_matrix")) trex_spca_select_mmap
                  else trex_spca_select
      private$result_ <- .r_positions(selector(
        private$X_ref,
        private$M_,
        private$tFDR_,
        private$ctrl,
        private$seed_,
        private$verbose_
      ))
      invisible(self)
    }
  ),
  active = list(

    #' @field Z Scores of the sparse components (n x M): column m is the centered \code{X}
    #'   times column m of \code{V}.
    Z = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      private$result_$Z
    },

    #' @field V Sparse loadings (p x M), each column of unit L2 norm. In \code{"ActiveSet"}
    #'   mode the non-zero entries of column m are the variables of \code{active_sets[[m]]}.
    #'   In \code{"Thresholded"} mode they are the largest ordinary PCA loadings, as many as
    #'   \code{active_sets[[m]]} has entries, which can be other variables.
    V = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      private$result_$V
    },

    #' @field active_sets Selected variables per component, a list of M integer vectors
    #'   (1-based), the result of the T-Rex+GVS runs. The FDR control holds for these sets. In
    #'   \code{"Thresholded"} mode they can differ from the non-zero entries of \code{V}.
    active_sets = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      private$result_$active_sets
    },

    #' @field adjusted_ev Adjusted explained variance per component (length M):
    #'   \code{r_mm^2 / (n - 1)}, with \code{r_mm} the m-th diagonal element of R in the QR
    #'   decomposition \code{Z = QR}. The adjustment accounts for sparse components that are
    #'   not orthogonal.
    adjusted_ev = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      private$result_$adjusted_ev
    },

    #' @field cumulative_ev Cumulative adjusted explained variance as a fraction of the total
    #'   variance of the centered \code{X} (length M), between 0 and 1. It is not the PEV of
    #'   Definition 1 of the paper, whose denominator (signal EV plus mixed EV) needs the true
    #'   support.
    cumulative_ev = function(value) {
      if (!missing(value)) stop("Field is read-only.")
      private$result_$cumulative_ev
    }
  )
)
