# =======================================================================================
# TRexBiobankScreeningSelector - Unit tests
# =======================================================================================

test_that("TRexBiobankScreeningSelector instantiation and data types work", {

  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)

  # Valid 1D
  y_vec <- rnorm(n)
  expect_silent(TRexBiobankScreeningSelector$new(X, y_vec, verbose = FALSE))

  # Valid 2D
  Y_mat <- matrix(rnorm(n * 3), n, 3)
  expect_silent(TRexBiobankScreeningSelector$new(X, Y_mat, verbose = FALSE))

  # Invalid type for Y
  expect_error(
    TRexBiobankScreeningSelector$new(X, "invalid", verbose = FALSE),
    "Response evaluation Y must be a numeric vector, a numeric matrix or an mmap_matrix"
  )

  # Dimension mismatch 1D
  y_bad <- rnorm(n + 10)
  expect_error(TRexBiobankScreeningSelector$new(X, y_bad, verbose = FALSE))


  # Dimension mismatch 2D
  Y_bad <- matrix(rnorm((n + 1) * 3), n + 1, 3)
  expect_error(TRexBiobankScreeningSelector$new(X, Y_bad, verbose = FALSE))
})



test_that("TRexBiobankScreeningSelector screens single phenotype correctly (1D)", {

  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)

  # Inject specific signal
  y <- X[, 1] * 5 + X[, 2] * -4 + rnorm(n)

  selector <- TRexBiobankScreeningSelector$new(X, y, verbose = FALSE)
  res <- selector$select()

  # Check object structure
  expect_type(res, "list")
  expect_true("phenotype_index" %in% names(res))
  expect_true("selected_indices" %in% names(res))
  expect_true("estimated_FDR" %in% names(res))
  expect_true("method_used" %in% names(res))
  expect_true("used_fallback_trex" %in% names(res))

  # R indexing base check
  expect_equal(res$phenotype_index, 1)
})



test_that("TRexBiobankScreeningSelector screens multiple phenotypes correctly (2D)", {

  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 3), n, 3)

  # Inject signals for some phenotypes
  Y[, 1] <- X[, 1] * 5 + rnorm(n)
  Y[, 2] <- X[, 2] * 4 + X[, 3] * -3 + rnorm(n)

  selector <- TRexBiobankScreeningSelector$new(X, Y, verbose = FALSE)
  res <- selector$select()

  # 2-D returns a plain list of per-phenotype records: res[[i]] has the same
  # shape as the 1-D return (mirrors Python's list[BiobankScreenTRexResult]).
  expect_type(res, "list")
  expect_equal(length(res), 3)

  record_fields <- c("phenotype_index", "selected_indices", "estimated_FDR",
                     "method_used", "estimated_FDR_screen_ordinary",
                     "estimated_FDR_screen_bootstrap",
                     "selected_indices_screen_ordinary",
                     "selected_indices_screen_bootstrap", "used_fallback_trex")

  for (i in seq_len(3)) {
    rec <- res[[i]]
    expect_type(rec, "list")
    expect_true(all(record_fields %in% names(rec)),
                label = paste("phenotype", i, "record fields"))
    expect_equal(rec$phenotype_index, i)
  }

  # The former $statistics data.frame is a one-line projection of the records.
  stats <- do.call(rbind, lapply(res, function(r) {
    as.data.frame(r[c("phenotype_index", "estimated_FDR",
                      "method_used", "used_fallback_trex")])
  }))
  expect_s3_class(stats, "data.frame")
  expect_equal(nrow(stats), 3)
})



test_that("TRexBiobankScreeningSelector correctly triggers fallback T-Rex decision logic", {

  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)

  # Strong noisy signal to trigger an FDR outside constraints
  y <- X[, 1] * 5 + rnorm(n)

  # Force fallback by setting an impossible FDR acceptance range for Screen-TRex
  selector <- TRexBiobankScreeningSelector$new(
    X, y,
    biobank_control = trex_biobank_control(lower_bound_FDR = 0.0001, upper_bound_FDR = 0.0002),
    verbose = FALSE
  )
  res <- selector$select()

  # With such extreme bounds, it should reject Screen-TRex evaluation and use the full T-Rex
  expect_true(res$used_fallback_trex)

  # A fallback result reports the target level (tFDR, default 0.1), the flag marks it, and the
  # screening estimates stay in their own fields
  expect_equal(res$estimated_FDR, 0.1)
  expect_equal(res$method_used, "T-Rex (fallback)")
  expect_true(res$estimated_FDR_screen_ordinary > 0 && res$estimated_FDR_screen_ordinary <= 1)
  expect_true(res$estimated_FDR_screen_bootstrap > 0 && res$estimated_FDR_screen_bootstrap <= 1)
})



test_that("TRexBiobankScreeningSelector protects against invalid framework strings", {

  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Invalid solver fallback
  expect_error(
    TRexBiobankScreeningSelector$new(X, y,
                                     control = trex_control(solver = "UNKNOWN_SOLVER"),
                                     verbose = FALSE),
    "Unknown TRex selector solver"
  )

  # Invalid core screening method parameter mapping
  expect_error(
    TRexBiobankScreeningSelector$new(X, y,
                                     screen_control = trex_screen_control(
                                       trex_method = "UNKNOWN_VARIANT"
                                     ),
                                     verbose = FALSE),
    "Unknown ScreenTRexMethod"
  )
})


test_that("TRexBiobankScreeningSelector 1D result has 1-based indices in [1, p]", {
  set.seed(42)
  n <- 50
  p <- 15
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 5 + X[, 2] * -4 + rnorm(n)

  selector <- TRexBiobankScreeningSelector$new(X, y, verbose = FALSE)
  res <- selector$select()

  # phenotype_index is 1 (single phenotype, 1-based)
  expect_equal(res$phenotype_index, 1L)

  # selected_indices are 1-based integers in [1, p]
  idx <- res$selected_indices
  expect_type(idx, "integer")
  expect_true(all(idx >= 1L & idx <= p))

  # sub-method index fields are present and also 1-based in [1, p]
  expect_true("selected_indices_screen_ordinary" %in% names(res))
  expect_true("selected_indices_screen_bootstrap" %in% names(res))

  for (field in c("selected_indices_screen_ordinary", "selected_indices_screen_bootstrap")) {
    sub_idx <- res[[field]]
    expect_type(sub_idx, "integer")
    expect_true(all(sub_idx >= 1L & sub_idx <= p), label = paste(field, "1-based in [1,p]"))
  }

  # estimated_FDR is a non-negative double
  expect_type(res$estimated_FDR, "double")
  expect_true(res$estimated_FDR >= 0.0)
})


test_that("TRexBiobankScreeningSelector 2D result has 1-based phenotype_index and per-phenotype indices", {

  set.seed(42)
  n <- 50
  p <- 15
  q <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  Y[, 1] <- X[, 1] * 5 + rnorm(n)
  Y[, 2] <- X[, 2] * 4 + X[, 3] * -3 + rnorm(n)

  selector <- TRexBiobankScreeningSelector$new(X, Y, verbose = FALSE)
  res <- selector$select()

  # 2-D returns a list of q per-phenotype records.
  expect_equal(length(res), q)

  # phenotype_index is 1-based: res[[i]]$phenotype_index == i
  expect_equal(vapply(res, function(r) r$phenotype_index, numeric(1)),
               as.numeric(seq_len(q)))

  # each record's index vectors are 1-based integers in [1, p]
  index_fields <- c("selected_indices", "selected_indices_screen_ordinary",
                    "selected_indices_screen_bootstrap")
  for (i in seq_len(q)) {
    rec <- res[[i]]
    for (field in index_fields) {
      idx <- rec[[field]]
      expect_type(idx, "integer")
      expect_true(all(idx >= 1L & idx <= p),
                  label = paste(field, "phenotype", i, "1-based in [1,p]"))
    }
  }
})


# =============================================================================
# Memory-mapped X: mapped directly instead of
# being read whole; results equal the in-memory path.
# =============================================================================

test_that("TRexBiobankScreeningSelector on a readwrite mmap_matrix equals the in-memory run", {
  set.seed(71)
  n <- 60; p <- 40
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 3 - X[, 5] * 2 + rnorm(n)
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f), add = TRUE)
  rw <- mem_matrix_from(X, f)
  X0 <- X + 0

  a <- TRexBiobankScreeningSelector$new(X, y, tFDR = 0.2, seed = 9, verbose = FALSE); a$select()
  b <- TRexBiobankScreeningSelector$new(rw, y, tFDR = 0.2, seed = 9, verbose = FALSE); b$select()
  expect_identical(b$selected_indices, a$selected_indices)
  expect_equal(as.matrix(rw), X0, tolerance = 1e-12)
})


test_that("TRexBiobankScreeningSelector takes a memory-mapped Y and only reads it", {
  set.seed(73)
  n <- 80; p <- 40; q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- sapply(seq_len(q), function(k) 4 * X[, k] - 3 * X[, k + 10] + rnorm(n))
  fx <- tempfile(fileext = ".bin")
  fy <- tempfile(fileext = ".bin")
  on.exit(unlink(c(fx, fy)), add = TRUE)
  mx <- mem_matrix_from(X, fx)
  my <- mem_matrix_from(Y, fy)
  ro <- mem_matrix(fy, rows = n, cols = q, mode = "readonly")
  md5 <- unname(tools::md5sum(fy))
  X0 <- X + 0

  run <- function(design, response) {
    TRexBiobankScreeningSelector$new(design, response, tFDR = 0.2, seed = 9, verbose = FALSE,
                                     screen_control = trex_screen_control(R_boot = 100),
                                     control = trex_control(K = 10))$select()
  }
  expected <- run(X + 0, Y)
  expect_length(expected, q)

  # Y in a file, in either mode, with X in memory and with X in a file
  for (res in list(run(X + 0, my), run(X + 0, ro), run(mx, ro))) {
    expect_length(res, q)
    for (k in seq_len(q)) {
      expect_identical(res[[k]]$selected_indices, expected[[k]]$selected_indices)
      expect_identical(res[[k]]$method_used, expected[[k]]$method_used)
    }
  }

  # The file of Y is as it was, and the file of X is restored
  expect_identical(as.matrix(ro), Y)
  expect_identical(unname(tools::md5sum(fy)), md5)
  expect_equal(as.matrix(mx), X0, tolerance = 1e-12)

  expect_error(TRexBiobankScreeningSelector$new(X + 0, "Y"), "numeric matrix or an mmap_matrix")
  short <- mem_matrix_from(Y[-1, ], tempfile(fileext = ".bin"))
  expect_error(TRexBiobankScreeningSelector$new(X + 0, short), "Dimension mismatch")
})


test_that("TRexBiobankScreeningSelector gives every phenotype its own seed", {

  set.seed(42)
  n <- 80
  p <- 40
  copies <- 12
  X <- matrix(rnorm(n * p), n, p)
  y_noise <- rnorm(n)

  # Twelve copies of one noise response make the dummy draw visible: the screening of such
  # a response selects something in about half of the draws
  Y <- matrix(y_noise, n, copies)
  screen_control <- trex_screen_control(R_boot = 100)
  multi <- TRexBiobankScreeningSelector$new(X + 0, Y, seed = 11, verbose = FALSE,
                                            screen_control = screen_control,
                                            control = trex_control(K = 10))$select()
  expect_length(multi, copies)

  # Phenotype k runs with seed + k - 1, and a single-phenotype run with that seed reproduces it
  for (k in seq_len(copies)) {
    single <- TRexBiobankScreeningSelector$new(X + 0, y_noise, seed = 11 + k - 1, verbose = FALSE,
                                               screen_control = screen_control,
                                               control = trex_control(K = 10))$select()
    for (field in c("selected_indices", "selected_indices_screen_ordinary",
                    "selected_indices_screen_bootstrap", "estimated_FDR_screen_ordinary",
                    "estimated_FDR_screen_bootstrap", "method_used")) {
      expect_equal(multi[[k]][[field]], single[[field]], label = paste("phenotype", k, field))
    }
  }

  # With one shared dummy draw the twelve copies would give twelve identical results
  first <- multi[[1]]$selected_indices_screen_ordinary
  differing <- sum(!vapply(multi, function(r) identical(r$selected_indices_screen_ordinary, first), logical(1)))
  expect_gt(differing, 0)
})
