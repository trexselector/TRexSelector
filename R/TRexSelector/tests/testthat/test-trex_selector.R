# =======================================================================================
# Unit tests for the R6 class TRexSelector, focusing on active bindings and getters.
# =======================================================================================

test_that("TRexSelector active bindings enforce read-only protection", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE)

  # Try to assign to an active binding (which should fail)
  expect_error(
    selector$v_thresh <- 0.5,
    "v_thresh is read-only"
  )
  expect_error(
    selector$selected_var <- c(TRUE, FALSE),
    "selected_var is read-only"
  )
  expect_error(
    selector$fdp_hat_mat <- matrix(0),
    "fdp_hat_mat is read-only"
  )
  expect_error(
    selector$.ptr <- NULL,
    ".ptr is read-only"
  )
})


test_that("TRexSelector active bindings return expected data types post-selection", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))

  # Before execution, properties may understandably be empty or matrix of 0
  # We invoke the algorithm directly
  selector$select()

  # Check vector state returns
  expect_type(selector$selected_indices, "integer")
  expect_type(selector$selected_var, "logical")
  expect_length(selector$selected_var, p)

  # Check numeric values
  expect_type(selector$v_thresh, "double")
  expect_type(selector$T_stop, "integer")
  expect_type(selector$L, "integer")
  expect_type(selector$y_mean, "double")

  # Check vector distributions
  expect_type(selector$x_means, "double")
  expect_length(selector$x_means, p)

  # Check the large matrix returns
  expect_true(is.matrix(selector$fdp_hat_mat))
  expect_true(is.matrix(selector$phi_mat))

})


test_that("TRexSelector selected_indices are 1-based and in [1, p]", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  idx <- selector$selected_indices
  expect_type(idx, "integer")
  expect_true(all(idx >= 1L & idx <= p))
})


test_that("TRexSelector selected_var and selected_indices are consistent", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  idx <- selector$selected_indices
  sv  <- selector$selected_var

  # selected_var has one entry per predictor
  expect_length(sv, p)
  # selected_indices == which(selected_var) (both 1-based)
  expect_equal(sort(idx), which(sv))
})


test_that("TRexSelector phi_prime is at most 1 with length p", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  phi_p <- selector$phi_prime
  expect_length(phi_p, p)
  expect_true(all(phi_p <= 1.0))  # capped at 1, may be negative
})


test_that("TRexSelector phi_mat has T_stop rows and p columns", {
  n <- 50
  p <- 10
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  selector <- TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 3))
  selector$select()

  ts <- selector$T_stop
  pm <- selector$phi_mat
  if (ts > 0L) {
    expect_equal(nrow(pm), ts)
    expect_equal(ncol(pm), p)
  }
})


test_that("TRexSelector multiple solver methods return 1-based indices consistent with selected_var", {
  n <- 100
  p <- 30
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 5] * 4.0 + X[, 15] * (-2.5) + rnorm(n)

  for (method in c("TLARS", "TLASSO", "TOMP")) {
    sel <- TRexSelector$new(X, y, tFDR = 0.1, verbose = FALSE,
                            control = trex_control(K = 20, solver = method))
    sel$select()

    idx <- sel$selected_indices
    sv  <- sel$selected_var

    expect_type(idx, "integer")
    expect_length(sv, p)

    expect_true(all(idx >= 1L & idx <= p), label = paste(method, "1-based in [1,p]"))
    expect_equal(sort(idx), which(sv), label = paste(method, "consistency"))
  }
})


test_that("exchangeable tie-breaking is plumbed through for greedy solvers", {
  # Control defaults: the band is left to the selector (-1), floor 0.5.
  ctrl <- trex_control()
  expect_equal(ctrl$exch_tie_alpha, -1)
  expect_equal(ctrl$exch_tie_floor, 0.5)

  # A correlated cluster so the exchangeable tie set is non-trivial.
  n <- 100
  p <- 12
  set.seed(42)
  base <- rnorm(n)
  X <- matrix(rnorm(n * p), n, p)
  for (j in 1:4) X[, j] <- base + 0.3 * rnorm(n)
  y <- base * 3.0 + rnorm(n)

  # Seeded selector (dummies and tie-breaks) and an unconditional
  # expectation, so the expectation count of the suite does not depend on
  # whether this run happened to select anything.
  for (method in c("TOMP", "TAFS")) {
    sel <- TRexSelector$new(
      X, y, tFDR = 0.2, seed = 1, verbose = FALSE,
      control = trex_control(K = 5, solver = method,
                             exch_tie_alpha = 0.25, exch_tie_floor = 0.5)
    )
    sel$select()
    idx <- sel$selected_indices
    expect_true(all(idx >= 1L & idx <= p),
                label = paste(method, "tie-broken run in [1,p]"))
  }
})


# =============================================================================
# restore(): give X back before select()
# =============================================================================

test_that("TRexSelector$restore() gives X back before select(); select() then errors", {
  set.seed(51)
  n <- 50; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  X0 <- X + 0
  sel <- TRexSelector$new(X, y, tFDR = 0.2, verbose = FALSE, control = trex_control(K = 5))
  expect_false(isTRUE(all.equal(X, X0)))          # normalized in place at construction
  sel$restore()
  expect_equal(X, X0, tolerance = 1e-12)
  sel$restore()                                    # idempotent
  expect_equal(X, X0, tolerance = 1e-12)
  expect_error(sel$select(), "restore")
  # A fresh selector on the same matrix runs and restores X itself
  sel2 <- TRexSelector$new(X, y, tFDR = 0.2, verbose = FALSE, control = trex_control(K = 5))
  sel2$select()
  expect_equal(X, X0, tolerance = 1e-12)
  sel2$restore()                                   # no-op after select()
  expect_equal(X, X0, tolerance = 1e-12)
})

test_that("the derived selectors inherit restore()", {
  set.seed(52)
  n <- 40; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 2] * 2 + rnorm(n)
  X0 <- X + 0
  sel <- TRexScreeningSelector$new(X, y, seed = 3, verbose = FALSE)
  sel$restore()
  expect_equal(X, X0, tolerance = 1e-12)
  expect_error(sel$select(), "restore")
})

test_that("a selector runs once: a second select() is an error and the results stay", {
  set.seed(54)
  n <- 50; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  X0 <- X + 0
  sel <- TRexSelector$new(X, y, tFDR = 0.2, seed = 3, verbose = FALSE, control = trex_control(K = 5))
  sel$select()
  first <- sel$selected_indices
  expect_error(sel$select(), "has run already")
  expect_error(sel$select(), "has run already")     # and every later call
  expect_identical(sel$selected_indices, first)
  expect_equal(X, X0, tolerance = 1e-12)
  # A new selector on the same matrix repeats the run
  again <- TRexSelector$new(X, y, tFDR = 0.2, seed = 3, verbose = FALSE, control = trex_control(K = 5))
  again$select()
  expect_identical(again$selected_indices, first)
})

test_that("the derived selectors run once as well", {
  set.seed(55)
  n <- 60; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 2] * 2 + rnorm(n)
  X0 <- X + 0
  makers <- list(
    screening = function() TRexScreeningSelector$new(X, y, seed = 3, verbose = FALSE),
    da = function() TRexDASelector$new(X, y, tFDR = 0.2, seed = 3, verbose = FALSE,
                                       control = trex_control(K = 5)),
    gvs = function() TRexGVSSelector$new(X, y, tFDR = 0.2, seed = 3, verbose = FALSE,
                                         gvs_control = trex_gvs_control(lambda_2 = 1),
                                         control = trex_control(solver = "TENET", K = 5)),
    tikhonov = function() TRexTikhonovSelector$new(X, y, tFDR = 0.2, seed = 3, verbose = FALSE,
                                                   tik_control = trex_tikhonov_control(
                                                     tikhonov_K = diag(p), lambda_2 = 1),
                                                   control = trex_control(solver = "TIKNET", K = 5))
  )
  for (name in names(makers)) {
    sel <- makers[[name]]()
    sel$select()
    expect_error(sel$select(), "has run already", info = name)
    expect_equal(X, X0, tolerance = 1e-12, info = name)
  }
})

test_that("a target FDR that is not a number is refused", {
  set.seed(56)
  X <- matrix(rnorm(40 * 10), 40, 10)
  y <- X[, 1] + rnorm(40)
  expect_error(TRexSelector$new(X, y, tFDR = NaN, verbose = FALSE), "tFDR must be in \\[0, 1\\]")
})

test_that("the selector control getters and num_dummies are exposed", {
  set.seed(53)
  n <- 40; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  sel <- TRexSelector$new(X, y, tFDR = 0.2, verbose = FALSE,
                          control = trex_control(K = 5, max_dummy_multiplier = 2))
  expect_equal(sel$tFDR, 0.2)
  expect_identical(sel$K, 5L)
  expect_identical(sel$max_num_dummies, 2L)
  expect_type(sel$max_T_stop, "logical")
  expect_type(sel$stagnation_check, "logical")
  sel$select()
  expect_identical(sel$num_selected, length(sel$selected_indices))
  expect_true(sel$num_dummies >= p)
})

test_that("SKIPL with a large budget warns only without the stagnation stop", {
  set.seed(42)
  X <- matrix(rnorm(50 * 10), 50, 10)
  y <- rnorm(50)
  warns <- function(strategy, multiplier, stop = NULL) {
    ctrl <- trex_control(solver = "TLARS", K = 3, lloop_strategy = strategy,
                         max_dummy_multiplier = multiplier, tloop_stagnation_stop = stop)
    out <- utils::capture.output(
      TRexSelector$new(X, y, verbose = FALSE, control = ctrl)$select(),
      type = "message"
    )
    any(grepl("lloop_strategy SKIPL", out, fixed = TRUE))
  }
  expect_true(warns("SKIPL", 10))
  expect_false(warns("SKIPL", 5))
  expect_false(warns("SKIPL", 20, TRUE))
  expect_false(warns("STANDARD", 20))
})


# =======================================================================================
# Solvers of the plain selector
# =======================================================================================

# Planted design: strong signal on the variables 1:5.
.planted_data <- function(n = 100, p = 20, data_seed = 42) {
  set.seed(data_seed)
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(5, 5), rep(0, p - 5))
  y <- as.vector(X %*% beta + rnorm(n))
  list(X = X, y = y, support = 1:5, p = p)
}


test_that("GVS-only solvers are rejected by the plain selector with guidance", {
  d <- .planted_data(n = 50, p = 10)

  run_with <- function(s) {
    selector <- TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 42,
                                 verbose = FALSE,
                                 control = trex_control(solver = s, K = 3))
    selector$select()
  }
  expect_error(run_with("TIENET"), "GVS-only solver")
  expect_error(run_with("TIENET_AUG"), "GVS-only solver")
  # TIKNET needs a penalty matrix: TRexTikhonovSelector only.
  expect_error(run_with("TIKNET"), "Tikhonov penalty matrix")
})


test_that("TENET_AUG runs in the plain selector and matches TENET", {
  d <- .planted_data()

  run_with <- function(s) {
    selector <- TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 42,
                                 verbose = FALSE,
                                 control = trex_control(solver = s, K = 5))
    selector$select()
    selector$selected_var
  }

  sel_tenet <- run_with("TENET")
  sel_aug   <- run_with("TENET_AUG")

  # TENET_AUG is the row-augmented formulation of the same elastic-net
  # problem (proved equivalent for lambda2 >= 0): identical data and seed
  # must select identically.
  expect_identical(sel_aug, sel_tenet)
  expect_true(all(d$support %in% which(sel_aug)))
})


test_that("an error from inside the random experiments reaches R, and the design comes back", {
  # The trigger is an option that the solver checks and the selector does not: TAFS refuses a
  # shrinkage step of 0 when it is constructed inside the experiments. select() then ends with
  # the error of the solver, the design is given back at once and the selector is retired.
  d <- .planted_data(n = 50, p = 10)
  X <- d$X
  X0 <- X + 0

  sel <- TRexSelector$new(X, d$y, tFDR = 0.2, seed = 4, verbose = FALSE,
                          control = trex_control(solver = "TAFS", rho_afs = 0, K = 5))
  expect_false(isTRUE(all.equal(X, X0)))          # standardized in place
  expect_error(sel$select(), "TAFS_Solver: rho must be in the range (0, 1].", fixed = TRUE)
  expect_equal(X, X0, tolerance = 1e-12)          # given back
  expect_error(sel$select(), "can no longer run") # retired
})
