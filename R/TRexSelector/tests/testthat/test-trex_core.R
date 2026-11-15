# =======================================================================================
# TRexSelector core - Unit tests
# =======================================================================================

test_that("TRexSelector correctly guards against invalid backend parameters", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # tFDR bounds constraints ([0, 1])
  expect_error(
    TRexSelector$new(X, y, tFDR = -0.1, verbose = FALSE),
    "tFDR must be in"
  )
  expect_error(
    TRexSelector$new(X, y, tFDR = 1.5, verbose = FALSE),
    "tFDR must be in"
  )


  # K value checks
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(K = 1)),
    "K must be >= 2"
  )


  # Invalid selector solver enum matching
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(solver = "INVALID_METHOD")),
    "Unknown TRex selector solver:"
  )


  # Invalid LLoop strategy enum matching
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE,
                     control = trex_control(solver = "SCREEN", lloop_strategy = "INVALID")),
    "Unknown LLoopStrategy:"
  )

  # max_dummy_multiplier
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(max_dummy_multiplier = 0)),
    "max_dummy_multiplier must be >= 1. Got: 0"
  )
})


test_that("TRexSelector handles concurrent execution properly (OpenMP)", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Ensure the object can be created and bound appropriately with concurrency settings
  selector_openmp <- expect_no_error(
    TRexSelector$new(X, y, verbose = FALSE,
                     control = trex_control(parallel_rnd_experiments = TRUE, max_outer_threads = 2))
  )

  # Run the algorithm targeting out threaded branches in C++ and falls back cleanly
  # if openMP is disabled.
  res_openmp <- expect_no_error(selector_openmp$select())

  expect_type(res_openmp, "environment")
  expect_true("selected_indices" %in% names(res_openmp))
})


test_that("a seeded selector selects the same variables whatever the thread options", {
  set.seed(7)
  n <- 80
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] * 3 - X[, 2] * 2 + rnorm(n))
  run <- function(ctrl) {
    sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 11, verbose = FALSE, control = ctrl)
    sel$select()
    sel$selected_indices
  }
  serial   <- run(trex_control(K = 4))
  parallel <- run(trex_control(K = 4, parallel_rnd_experiments = TRUE, max_outer_threads = 2))
  budget   <- run(trex_control(K = 4, max_inner_threads = 2))
  both     <- run(trex_control(K = 4, parallel_rnd_experiments = TRUE, max_outer_threads = 2,
                               max_inner_threads = 2))
  expect_identical(parallel, serial)
  expect_identical(budget, serial)
  expect_identical(both, serial)

  # use_openmp was renamed parallel_rnd_experiments
  expect_error(trex_control(use_openmp = TRUE), "unused argument")
})


test_that("TRexSelector supports memory mapping workflow correctly", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Activate the memmap_manager_ branch deeply embedded in the backend utilities
  selector_mmap <- expect_no_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(use_memory_mapping = TRUE))
  )

  # Validate the full data stream workflow holds mapping pointers open effectively
  res_mmap <- expect_no_error(selector_mmap$select())

  expect_type(res_mmap, "environment")
  expect_true("selected_indices" %in% names(res_mmap))
})


test_that("TRexSelector correctly guards against edge-case hyperparameter inputs", {

  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # opt_threshold bounds constraints ([0, 1])
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(opt_threshold = -0.5)),
    "opt_threshold"
  )
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(opt_threshold = 1.5)),
    "opt_threshold"
  )

  # tloop_max_stagnant_steps is at least 1
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(tloop_max_stagnant_steps = -1)),
    "tloop_max_stagnant_steps"
  )
  expect_error(trex_control(tloop_max_stagnant_steps = 0), "tloop_max_stagnant_steps must be >= 1. Got: 0",
               fixed = TRUE)
  # a hand-built control reaches the check of the core
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = list(tloop_max_stagnant_steps = 0)),
    "tloop_max_stagnant_steps must be >= 1. Got: 0", fixed = TRUE
  )
  expect_silent(trex_control(tloop_max_stagnant_steps = 1))

  # The thread options are whole numbers of at least 1 (there is no "every core" value)
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(max_outer_threads = -1)),
    "max_outer_threads"
  )
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(max_inner_threads = -1)),
    "max_inner_threads"
  )
  expect_error(trex_control(max_outer_threads = 0), "max_outer_threads")
  expect_error(trex_control(max_inner_threads = 0), "max_inner_threads")
  ctrl_bad <- trex_control()
  ctrl_bad$max_outer_threads <- 0   # bypass the R check: the backend rejects it too
  expect_error(TRexSelector$new(X, y, verbose = FALSE, control = ctrl_bad), "max_outer_threads")

  # lambda2 (ENET constraint) and rho_afs (AFS constraint) bounds
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(lambda2 = -0.1)),
    "lambda2"
  )
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(rho_afs = -0.5)),
    "rho_afs"
  )

  # ncgmp_variant valid types (assumed boolean / specific int)
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(ncgmp_variant = 5)),
    "ncgmp_variant"
  )

  # Explicitly unsupported solver variants
  expect_error(
    TRexSelector$new(X, y, verbose = FALSE, control = trex_control(solver = "UNKNOWN_SOLVER")),
    "Unknown TRex selector solver"
  )
})

# =======================================================================================
# trex_dummy_distribution() — structural tests
# =======================================================================================

test_that("trex_dummy_distribution returns correct list structure for Normal (default)", {
  d <- trex_dummy_distribution()
  expect_type(d, "list")
  expect_equal(d$type, "Normal")
  expect_null(d$df)
})

test_that("trex_dummy_distribution returns correct list structure for StudentT", {
  d <- trex_dummy_distribution("StudentT", df = 3)
  expect_equal(d$type, "StudentT")
  expect_equal(d$df, 3)
})

test_that("trex_dummy_distribution returns correct list structure for Uniform", {
  d <- trex_dummy_distribution("Uniform", a = 2.0)
  expect_equal(d$type, "Uniform")
  expect_equal(d$a, 2.0)
})

test_that("trex_dummy_distribution rejects unknown type", {
  expect_error(trex_dummy_distribution("GaussianMixture"), "Unknown distribution type")
})

test_that("trex_dummy_distribution passes the Mammen probabilities on", {
  d <- trex_dummy_distribution("Mammen", p1 = 0.7, p2 = 0.3)
  expect_equal(d$type, "Mammen")
  expect_equal(d$p1, 0.7)
  expect_equal(d$p2, 0.3)
  expect_null(trex_dummy_distribution("Mammen")$p1)
})

test_that("trex_control includes dummy_distribution field", {
  ctrl <- trex_control()
  expect_true("dummy_distribution" %in% names(ctrl))
  expect_equal(ctrl$dummy_distribution$type, "Normal")
})

test_that("trex_control passes custom dummy_distribution through", {
  ctrl <- trex_control(dummy_distribution = trex_dummy_distribution("Rademacher"))
  expect_equal(ctrl$dummy_distribution$type, "Rademacher")
})

# =======================================================================================
# trex_dummy_distribution() — integration tests (TRex actually runs). Each takes
# about 10 ms, so they run everywhere: no skip_on_cran(), which also skipped
# them under plain test_dir() / R CMD check (NOT_CRAN unset).
# =======================================================================================

test_that("TRexSelector runs with Rademacher dummy distribution", {
  set.seed(1)
  n <- 60; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + X[, 2] + rnorm(n)
  sel <- TRexSelector$new(X, y, tFDR = 0.2, verbose = FALSE,
                          control = trex_control(
                            K = 5,
                            dummy_distribution = trex_dummy_distribution("Rademacher")
                          ))
  result <- sel$select()
  expect_type(result$selected_var, "logical")
  expect_length(result$selected_var, p)
})

test_that("the Mammen probabilities reach the dummy draw", {
  set.seed(3)
  n <- 60; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + X[, 2] + rnorm(n)
  run <- function(...) {
    sel <- TRexSelector$new(X + 0, y, tFDR = 0.2, seed = 11, verbose = FALSE,
                            control = trex_control(
                              K = 5,
                              dummy_distribution = trex_dummy_distribution("Mammen", ...)
                            ))
    sel$select()
    sel$phi_mat
  }
  golden <- (sqrt(5) + 1) / (2 * sqrt(5))
  # the defaults of the core are the golden-ratio probabilities
  expect_identical(run(), run(p1 = golden, p2 = 1 - golden))
  # other probabilities draw other dummies
  expect_false(identical(run(), run(p1 = 0.5, p2 = 0.5)))
})

test_that("the sparse Rademacher dummies run with the default s = 0.1", {
  set.seed(4)
  n <- 60; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + X[, 2] + rnorm(n)
  run <- function(...) {
    sel <- TRexSelector$new(X + 0, y, tFDR = 0.2, seed = 11, verbose = FALSE,
                            control = trex_control(
                              K = 5,
                              dummy_distribution = trex_dummy_distribution("ConstrainedSparseRademacher", ...)
                            ))
    sel$select()
    sel$phi_mat
  }
  expect_identical(run(), run(s = 0.1))
})

# =======================================================================================
# The errors of the data normalizer count rows and columns from 1 in R.
# =======================================================================================

test_that("a constant column is named by its R position", {
  set.seed(5)
  n <- 40; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  X[, 5] <- 1
  y <- rnorm(n)
  msg <- tryCatch(TRexSelector$new(X + 0, y, verbose = FALSE), error = conditionMessage)
  expect_match(msg, "column 5 of X is constant", fixed = TRUE)
  expect_false(grepl("0-based", msg, fixed = TRUE))
  expect_match(msg, "scaling_mode = \"NONE\"", fixed = TRUE)
  expect_error(trex_normalize_design(X, y), "column 5 of X is constant", fixed = TRUE)
  X[, 5] <- 0
  expect_error(TRexSelector$new(X + 0, y, verbose = FALSE, control = trex_control(scaling_mode = "NONE")),
               "column 5 of X has L2 norm", fixed = TRUE)
})

test_that("an entry that is not finite is named by its R position", {
  set.seed(6)
  n <- 40; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  X[3, 7] <- NA
  y <- rnorm(n)
  msg <- tryCatch(TRexSelector$new(X + 0, y, verbose = FALSE), error = conditionMessage)
  expect_match(msg, "X[3, 7] is not finite", fixed = TRUE)
  expect_match(msg, "impute column 7 before", fixed = TRUE)
  expect_false(grepl("0-based", msg, fixed = TRUE))
  expect_error(TRexDASelector$new(X + 0, y, verbose = FALSE), "X[3, 7] is not finite", fixed = TRUE)
})

test_that("an error without a position passes unchanged", {
  set.seed(7)
  X <- matrix(rnorm(40 * 12), 40, 12)
  expect_error(TRexSelector$new(X, rnorm(40), tFDR = 2, verbose = FALSE), "Target FDR tFDR must be in [0, 1]",
               fixed = TRUE)
})

test_that("TRexSelector runs with StudentT dummy distribution (df=3)", {
  set.seed(2)
  n <- 60; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + X[, 2] + rnorm(n)
  sel <- TRexSelector$new(X, y, tFDR = 0.2, verbose = FALSE,
                          control = trex_control(
                            K = 5,
                            dummy_distribution = trex_dummy_distribution("StudentT", df = 3)
                          ))
  result <- sel$select()
  expect_type(result$selected_var, "logical")
  expect_length(result$selected_var, p)
})
