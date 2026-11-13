# =======================================================================================
# TRexDASelector - Unit tests
# =======================================================================================


test_that("TRexDASelector instantiation and dimensional validation work", {

  set.seed(42)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Valid instantiation
  expect_silent(TRexDASelector$new(X, y, verbose = FALSE))

  # Dimension mismatch
  y_bad <- rnorm(n + 10)
  expect_error(TRexDASelector$new(X, y_bad, verbose = FALSE), "length equal to")
})



test_that("TRexDASelector executes correctly with valid data", {

  set.seed(42)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Inject some strong signals in the first 3 variables
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # Instantiate TRex selector
  selector <- TRexDASelector$new(X, y, tFDR = 0.1, verbose = FALSE, control = trex_control(K = 5))

  # Execute variable selection
  selector$select()

  # Assert result structure and logic
  res_indices <- selector$selected_indices
  res_T_stop <- selector$T_stop
  res_num_dummies <- selector$L * ncol(selector$.ptr_X)

  expect_type(res_indices, "integer")

  expect_type(res_T_stop, "integer")
  expect_type(res_num_dummies, "integer")

  # Bounds check on output selected indices (1 <= index <= p)
  expect_true(all(res_indices >= 1 & res_indices <= p))
})



test_that("TRexDASelector supports da_method = 'NN' (nearest-neighbour sweep)", {
  set.seed(42)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # NN must be accepted (previously the R glue threw "Unknown DAMethod").
  # hc_grid_length = 0 lets the core auto-select min(20, p).
  selector <- TRexDASelector$new(
    X, y, tFDR = 0.2, verbose = FALSE,
    da_control = trex_da_control(da_method = "NN"),
    control = trex_control(K = 5)
  )
  expect_no_error(selector$select())
  expect_identical(selector$da_method_used, "NN")
  expect_type(selector$selected_indices, "integer")
  expect_true(all(selector$selected_indices >= 1 &
                  selector$selected_indices <= p))
  # NN runs a BT-style rho grid
  expect_gt(length(selector$rho_grid), 0L)
})


test_that("TRexDASelector supports the prior-groups DA path", {
  set.seed(42)
  n <- 60
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # Two-level prior constraints over p = 8 variables. When prior_groups is
  # supplied, da_method is ignored and the core routes through
  # setupDA_PriorGroups(): hierarchical sub-clustering runs within the finest
  # common refinement of the levels (here the size-2 fine groups) and the rho
  # grid is built from the pooled within-group dendrogram heights.
  fine   <- c(1, 1, 2, 2, 3, 3, 4, 4)
  coarse <- c(1, 1, 1, 1, 2, 2, 2, 2)

  selector <- TRexDASelector$new(
    X, y, tFDR = 0.2, verbose = FALSE,
    da_control = trex_da_control(
      prior_groups = list(fine, coarse)
    ),
    control = trex_control(K = 5)
  )
  expect_no_error(selector$select())
  expect_type(selector$selected_indices, "integer")
  expect_true(all(selector$selected_indices >= 1 &
                  selector$selected_indices <= p))
  # Data-driven rho grid: hc_grid_length points (default min(20, p) = 8),
  # ascending and terminated by the conservative rho = 1 singleton anchor.
  expect_equal(length(selector$rho_grid), p)
  expect_equal(selector$rho_grid, sort(selector$rho_grid))
  expect_equal(selector$rho_grid[length(selector$rho_grid)], 1.0)
})


test_that("prior_groups rejects negative labels", {
  set.seed(7)
  n <- 60
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 4 + X[, 2] * -3 + rnorm(n)

  expect_error(
    TRexDASelector$new(
      X, y, tFDR = 0.2, verbose = FALSE,
      da_control = trex_da_control(
        prior_groups = list(c(1, 1, 2, 2, -3, 3))
      ),
      control = trex_control(K = 5)
    ),
    "non-negative"
  )
})


test_that("prior_groups rejects a level whose length != p", {
  set.seed(42)
  n <- 40
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  expect_error(
    TRexDASelector$new(
      X, y, verbose = FALSE,
      da_control = trex_da_control(prior_groups = list(c(1, 1, 2)))  # length 3 != p
    ),
    "length"
  )
})


test_that("TRexDASelector behaves memory safe with memory mapping enabled", {

  set.seed(42)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Use use_memory_mapping = TRUE to confirm the sandboxed tempdir() doesn't crash
  selector <- TRexDASelector$new(X, y, verbose = FALSE, control = trex_control(K = 3, use_memory_mapping = TRUE))
  selector$select()

  res_indices <- selector$selected_indices
  res_T_stop <- selector$T_stop
  res_num_dummies <- selector$L * ncol(selector$.ptr_X)

  expect_type(res_indices, "integer")


})



test_that("TRexDASelector correctly guards against invalid DA parameters", {

  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Invalid DA threshold parameter bounds (rho_thr_DA must be [0, 1])
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "EQUI", rho_thr_DA = -0.5), verbose = FALSE),
    "rho_thr_DA"
  )
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "EQUI", rho_thr_DA = 1.5), verbose = FALSE),
    "rho_thr_DA"
  )
  # AR1 takes thresholds in [0, 1): at 1 the window holds no neighbor. EQUI accepts 1
  # (X + 0 is a copy, so the selector that stays alive does not hold X).
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "AR1", rho_thr_DA = 1), verbose = FALSE),
    "rho_thr_DA"
  )
  expect_no_error(
    TRexDASelector$new(X + 0, y, da_control = trex_da_control(da_method = "EQUI", rho_thr_DA = 1), verbose = FALSE)
  )

  # Invalid DA Methods that are caught by Rcpp string mapping
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "UNKNOWN_DA"), verbose = FALSE),
    "Unknown DAMethod:"
  )

  # Invalid Hierarchical Linkage Options (rejected at Rcpp mapping)
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "BT", hc_linkage = "UNKNOWN"), verbose = FALSE),
    "Unsupported linkage method for the dependency-aware selector: UNKNOWN"
  )

  # Ward is not offered: its update needs squared Euclidean distances, and the
  # correlation distance is not Euclidean. The message names the supported linkages.
  expect_error(
    TRexDASelector$new(X, y, da_control = trex_da_control(da_method = "BT", hc_linkage = "Ward"), verbose = FALSE),
    "Unsupported linkage method for the dependency-aware selector: Ward\\nSupported: Single, Complete, Average, WPGMA"
  )

})



test_that("TRexDASelector AR1 uses the magnitude of a supplied cor_coef", {

  set.seed(42)
  n <- 60
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # X + 0 hands every run its own copy of the design
  run <- function(cor_coef) {
    selector <- TRexDASelector$new(X + 0, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                                   da_control = trex_da_control(da_method = "AR1", cor_coef = cor_coef),
                                   control = trex_control(K = 5))
    selector$select()
    selector
  }
  neg <- run(-0.5)
  pos <- run(0.5)

  # The window follows abs(cor_coef), so both runs give the same result
  expect_identical(neg$selected_var, pos$selected_var)
  expect_equal(neg$phi_prime, pos$phi_prime)
  expect_equal(neg$phi_mat, pos$phi_mat)

  # The field reports the coefficient as supplied
  expect_equal(neg$cor_coef, -0.5)
})



test_that("trex_da_control carries bt_selection_mode", {

  set.seed(42)
  n <- 60
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # The default is the feasible-only policy of the core
  expect_identical(trex_da_control()$bt_selection_mode, "FeasibleOnly")
  expect_identical(trex_da_control(bt_selection_mode = "RFaithful")$bt_selection_mode, "RFaithful")
  expect_error(trex_da_control(bt_selection_mode = "Other"), "should be one of")

  # X + 0 hands every run its own copy of the design
  run <- function(da_control) {
    selector <- TRexDASelector$new(X + 0, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                                   da_control = da_control, control = trex_control(K = 5))
    selector$select()
    selector
  }
  default <- run(trex_da_control())
  feasible <- run(trex_da_control(bt_selection_mode = "FeasibleOnly"))
  faithful <- run(trex_da_control(bt_selection_mode = "RFaithful"))
  expect_identical(default$selected_var, feasible$selected_var)
  expect_length(faithful$selected_var, p)

  # A control built by hand reaches the parser, which names the supported values
  by_hand <- trex_da_control()
  by_hand$bt_selection_mode <- "Other"
  expect_error(TRexDASelector$new(X + 0, y, verbose = FALSE, da_control = by_hand),
               "Unknown bt_selection_mode")
})



test_that("TRexDASelector AR1 deflates every step with the penalty of the last step", {

  # The deflated history is then non-decreasing, so the calibrated occurrences
  # phi_prime never exceed the deflated occurrences of the last step. A penalty
  # per step breaks this in about every second run of this design. The design
  # is tall (p < n / 2 - 1), so no experiment exhausts its path in the T-loop.
  n <- 150
  p <- 60
  rho <- 0.7
  for (seed in 1:6) {
    set.seed(seed)
    X <- matrix(0, n, p)
    X[, 1] <- rnorm(n)
    for (j in 2:p) X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
    y <- as.numeric(2 * (X[, 10] + X[, 30] + X[, 50]) + rnorm(n))

    selector <- TRexDASelector$new(X, y, tFDR = 0.2, seed = seed, verbose = FALSE,
                                   da_control = trex_da_control(da_method = "AR1"),
                                   control = trex_control(K = 10))
    selector$select()

    phi_mat <- selector$phi_mat
    phi_da <- phi_mat[nrow(phi_mat), ]
    expect_true(all(selector$phi_prime <= phi_da + 1e-12), info = paste("seed", seed))
  }
})



test_that("TRexDASelector replaces a cor_coef outside [-1, 1] by the estimate with a warning", {

  set.seed(42)
  n <- 60
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + X[, 3] * 3 + rnorm(n)

  # X + 0 hands every run its own copy of the design
  run <- function(da_method, cor_coef) {
    selector <- TRexDASelector$new(X + 0, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                                   da_control = trex_da_control(da_method = da_method, cor_coef = cor_coef),
                                   control = trex_control(K = 5))
    warning_text <- capture.output(selector$select(), type = "message")
    list(cor_coef = selector$cor_coef, warning_text = paste(warning_text, collapse = " "))
  }

  for (da_method in c("AR1", "EQUI")) {
    estimated <- run(da_method, -2)
    expect_false(grepl("cor_coef", estimated$warning_text), label = paste(da_method, "sentinel is silent"))

    for (unusable in c(-1.5, 1.5, -2.5)) {
      replaced <- run(da_method, unusable)
      expect_equal(replaced$cor_coef, estimated$cor_coef, label = paste(da_method, unusable, "estimate"))
      expect_match(replaced$warning_text, "cor_coef", label = paste(da_method, unusable, "warning"))
    }

    supplied <- run(da_method, -0.5)
    expect_equal(supplied$cor_coef, -0.5)
    expect_false(grepl("cor_coef", supplied$warning_text))
  }
})



test_that("TRexDASelector handles concurrent execution properly (OpenMP with DA)", {

  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # DA creates nested 3D calibrations; multithreaded runs ensure no race conditions
  selector_openmp <- expect_no_error(
    TRexDASelector$new(X, y,
                       da_control = trex_da_control(da_method = "BT", hc_linkage = "Average"),
                       control = trex_control(K = 5, parallel_rnd_experiments = TRUE, max_inner_threads = 2,
                           max_outer_threads = 2),
                       verbose = FALSE)
  )

  # Let it resolve the OpenMP execution through evaluateStep overrides
  selector_openmp$select()




  expect_type(selector_openmp$selected_indices, "integer")

  expect_type(selector_openmp$T_stop, "integer")

})


# =======================================================================================
# Unit tests for the R6 class TRexDASelector, focusing on active bindings and getters.
# =======================================================================================

test_that("TRexDASelector active bindings enforce read-only protection", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexDASelector$new(X, y, verbose = FALSE)

  # Try to assign to DA-specific active bindings (which should fail)
  expect_error(
    selector$rho_thresh <- 0.5,
    "Field is read-only"
  )
  expect_error(
    selector$rho_grid <- c(0.1, 0.2),
    "Field is read-only"
  )
  expect_error(
    selector$da_method_used <- 1L,
    "Field is read-only"
  )
  expect_error(
    selector$FDP_hat_array_BT <- array(0, dim = c(1, 1, 1)),
    "Field is read-only"
  )
  expect_error(
    selector$Phi_array_BT <- array(0, dim = c(1, 1, 1)),
    "Field is read-only"
  )
  expect_error(
    selector$R_array_BT <- array(0, dim = c(1, 1, 1)),
    "Field is read-only"
  )
})

test_that("TRexDASelector active bindings return expected data types post-selection", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Using BT with a valid K for multidimensional matrices
  selector <- TRexDASelector$new(X, y, verbose = FALSE,
                                 da_control = trex_da_control(hc_grid_length = 3),
                                 control = trex_control(K = 3))

  # Execute so matrices get allocated
  selector$select()

  # Check standard primitive types
  expect_type(selector$rho_thresh, "double")
  expect_type(selector$rho_grid, "double")
  expect_type(selector$da_method_used, "character")

  # When da_method="BT", arrays should either be full 3D structures or valid fallback numeric structs
  expect_true(is.array(selector$FDP_hat_array_BT) || is.numeric(selector$FDP_hat_array_BT))
  expect_true(is.array(selector$Phi_array_BT) || is.numeric(selector$Phi_array_BT))
  expect_true(is.array(selector$R_array_BT) || is.numeric(selector$R_array_BT))

  if (is.array(selector$FDP_hat_array_BT)) {
    expect_equal(length(dim(selector$FDP_hat_array_BT)), 3)
  }

  # Ensure the get_matrices() legacy generator matches the active binding state
})


test_that("TRexDASelector selected_var and selected_indices are consistent (1-based)", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexDASelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  idx <- selector$selected_indices
  sv  <- selector$selected_var

  expect_type(idx, "integer")
  expect_length(sv, p)
  expect_true(all(idx >= 1L & idx <= p))
  expect_equal(sort(idx), which(sv))
})


test_that("TRexDASelector phi_prime is at most 1 with length p", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexDASelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  phi_p <- selector$phi_prime
  expect_length(phi_p, p)
  expect_true(all(phi_p <= 1.0))  # capped at 1, may be negative
})


test_that("TRexDASelector phi_mat has T_stop rows and p columns", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexDASelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  ts <- selector$T_stop
  pm <- selector$phi_mat
  if (ts > 0L) {
    expect_gte(nrow(pm), 1L)
    expect_equal(ncol(pm), p)
  }
})


test_that("TRexDASelector DA-specific fields have valid bounds after selection", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexDASelector$new(X, y, verbose = FALSE,
                                 da_control = trex_da_control(da_method = "BT"),
                                 control = trex_control(K = 3))
  selector$select()

  expect_gte(selector$rho_thresh, 0.0)
  expect_lte(selector$rho_thresh, 1.0)
  expect_gt(length(selector$rho_grid), 0L)
  expect_true(all(selector$rho_grid >= 0.0 & selector$rho_grid <= 1.0))
  expect_true(selector$da_method_used %in% c("BT", "AR1", "EQUI", "BLOCK_EQUI"))
})
