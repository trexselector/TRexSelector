# =======================================================================================
# TRexTikhonovSelector - Unit tests (general-Tikhonov T-Rex through R)
# =======================================================================================

# Grouped design shared by the anchor tests: blocks of 5 correlated
# columns, signal on the first three variables of block 1.
#
# The dimensions mirror the C++ FOLDED anchor (n = 150, p = 40,
# rho = 0.75) DELIBERATELY: the FOLDED sparse-K penalty and GVS's GROUP
# hooks are mathematically identical but accumulate in different
# floating-point orders (w += lambda2*delta*K.col vs sigma_m += delta,
# multiplied once), so selections can differ at ULP-level margins. On a
# knife-edge design (n = 80, rho = 0.85) the GVS anchor flipped 2
# borderline variables; keep the margins wide here.
.make_tik_data <- function(n = 150, p = 40, rho = 0.75, data_seed = 42) {
  set.seed(data_seed)
  n_groups <- p / 5
  z <- matrix(rnorm(n * n_groups), n, n_groups)
  X <- matrix(0, n, p)
  for (g in 1:n_groups) {
    for (j in 1:5) {
      X[, (g - 1) * 5 + j] <- rho * z[, g] + sqrt(1 - rho^2) * rnorm(n)
    }
  }
  y <- as.vector(rowSums(X[, 1:3]) * 3 + rnorm(n))
  groups <- rep(1:n_groups, each = 5)
  list(X = X, y = y, groups = groups, p = p)
}

# Group-mean K = sum_m 1_m 1_m^T / p_m for equally-sized groups of 5.
.group_mean_k <- function(p, group_size = 5) {
  K <- matrix(0, p, p)
  for (g in seq_len(p / group_size)) {
    idx <- ((g - 1) * group_size + 1):(g * group_size)
    K[idx, idx] <- 1 / group_size
  }
  K
}

# A base control for the Tikhonov selector: a control that is passed names the solver.
.tik_base <- function(...) trex_control(solver = "TIKNET", ...)


test_that("trex_tikhonov_control validates its arguments", {
  expect_error(trex_tikhonov_control(), "tikhonov_K is required")
  expect_error(
    trex_tikhonov_control(tikhonov_K = diag(4), lambda2_method = "CV_1SE_CCD"),
    "arg"  # match.arg rejects the EN-shaped name: only IEN-geometry tuners
  )
  ctrl <- trex_tikhonov_control(tikhonov_K = diag(4))
  expect_identical(ctrl$lambda2_method, "PILOT_KNEE")
  expect_error(
    trex_tikhonov_control(tikhonov_K = diag(4), lambda2_method = "PILOT_MAX"),
    "arg"  # the argmax variant is gone; the knee with tolerance 0 reads it
  )
  expect_identical(ctrl$penalty_structure, "SPARSE_GENERAL")
  expect_true(ctrl$cluster_dummies)
  expect_identical(ctrl$cv_n_lambda1, 40L)
  expect_identical(ctrl$cv_n_lambda2, 100L)
  expect_identical(ctrl$cv_dfmax_ratio, 0.8)
  wide <- trex_tikhonov_control(tikhonov_K = diag(4), cv_dfmax_ratio = 1)
  expect_identical(wide$cv_dfmax_ratio, 1)
  expect_error(
    trex_tikhonov_control(tikhonov_K = diag(4), cv_dfmax_ratio = NaN),
    "cv_dfmax_ratio"
  )

  # GROUP_CENTERING builds K itself: no tikhonov_K needed.
  gc <- trex_tikhonov_control(penalty_structure = "GROUP_CENTERING",
                              lambda2_method = "PILOT_KNEE")
  expect_null(gc$tikhonov_K)
  expect_identical(gc$lambda2_method, "PILOT_KNEE")
  expect_identical(gc$pilot_statistic, "MAJORITY")

  # Pilot gates.
  expect_error(trex_tikhonov_control(tikhonov_K = diag(4), pilot_K = 1),
               "pilot_K must be >= 2")
  expect_error(trex_tikhonov_control(tikhonov_K = diag(4), pilot_T = 0),
               "pilot_T must be >= 1")
  expect_error(trex_tikhonov_control(tikhonov_K = diag(4), pilot_grid = c(1, -1)),
               "pilot_grid")
  expect_error(trex_tikhonov_control(tikhonov_K = diag(4),
                                     pilot_group_fraction = 1.5),
               "pilot_group_fraction")
})


test_that("construction gates translate as R errors", {
  d <- .make_tik_data()

  # Mis-dimensioned K.
  expect_error(
    TRexTikhonovSelector$new(
      d$X, d$y, verbose = FALSE,
      tik_control = trex_tikhonov_control(tikhonov_K = diag(d$p - 1))),
    "[Kk]"
  )

  # PERMUTATION-family strategies cannot carry cluster-MVN dummies.
  expect_error(
    TRexTikhonovSelector$new(
      d$X, d$y, verbose = FALSE,
      tik_control = trex_tikhonov_control(tikhonov_K = diag(d$p)),
      control = .tik_base(lloop_strategy = "PERMUTATION")),
    "PERMUTATION"
  )

  # Malformed prior_groups (length mismatch).
  expect_error(
    TRexTikhonovSelector$new(
      d$X, d$y, verbose = FALSE,
      tik_control = trex_tikhonov_control(tikhonov_K = diag(d$p),
                                          prior_groups = rep(1, d$p - 1))),
    "prior_groups"
  )
})


test_that("a base control that is passed has to name the solver TIKNET", {
  d <- .make_tik_data()
  X0 <- d$X + 0
  tik <- trex_tikhonov_control(tikhonov_K = diag(d$p), lambda_2 = 0.5)
  run <- function(control) {
    sel <- TRexTikhonovSelector$new(d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
                                    tik_control = tik, control = control)
    sel$select()
    sel$selected_indices
  }

  # Every solver outside the family is an error, in the words of R. TLARS is
  # one of them, although it is the default of trex_control().
  for (solver in c("TLARS", "TLASSO", "TENET", "TOMP", "TIENET")) {
    expect_error(run(trex_control(solver = solver)),
                 "control\\$solver must be \"TIKNET\" for the Tikhonov selector. Got: \"",
                 info = solver)
    expect_error(run(trex_control(solver = solver)), solver, fixed = TRUE, info = solver)
  }
  expect_error(run(trex_control()), "leave control out")
  expect_error(run(trex_control(K = 5)), "TLARS")
  # No selector was built, so the design is untouched
  expect_identical(d$X, X0)

  # Without a base control the one of the C++ class applies, whose solver is
  # TIKNET: the run equals the one with the defaults of trex_control() and TIKNET
  expect_identical(run(NULL), run(trex_control(solver = "TIKNET")))
  expect_gt(length(run(NULL)), 0)

  # The power rule takes the base control the same way
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y, tik_control = tik, seed = 1,
                                                grid = c(0, 1), control = trex_control()),
               "control\\$solver must be \"TIKNET\"")
  a <- trex_tikhonov_lambda2_power_rule(d$X, d$y, tFDR = 0.2, tik_control = tik, seed = 42,
                                        grid = c(0, 1))
  b <- trex_tikhonov_lambda2_power_rule(d$X, d$y, tFDR = 0.2, tik_control = tik, seed = 42,
                                        grid = c(0, 1), control = trex_control(solver = "TIKNET"))
  expect_identical(a$selected, b$selected)
  expect_identical(a$counts, b$counts)
  expect_equal(d$X, X0, tolerance = 1e-12)
})


test_that("GROUP_CENTERING with one group per variable is the zero-K path", {
  d <- .make_tik_data()

  # A HAC cut above the within-group correlation yields singletons: K_GC
  # vanishes and the path is the plain lasso on per-variable cluster
  # dummies. The run selects, and lambda_2 has no effect on it.
  run <- function(lambda_2) {
    sel <- TRexTikhonovSelector$new(
      d$X + 0, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
      tik_control = trex_tikhonov_control(
        penalty_structure = "GROUP_CENTERING",
        prior_groups = seq_len(d$p), lambda_2 = lambda_2),
      control = .tik_base(K = 5, max_dummy_multiplier = 2))
    sel$select()
    sel$selected_var
  }

  expect_gt(sum(run(10)), 0)
  expect_identical(run(10), run(0))
})


test_that("the pilot rule picks lambda_2 from its grid and records the profile", {
  d <- .make_tik_data()

  sel <- TRexTikhonovSelector$new(
    d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
    tik_control = trex_tikhonov_control(
      tikhonov_K = .group_mean_k(d$p), prior_groups = d$groups,
      lambda_2 = -1, lambda2_method = "PILOT_KNEE",
      pilot_grid = c(10, 0, 1), pilot_K = 5, pilot_T = 3),
    control = .tik_base(K = 5, max_dummy_multiplier = 2))
  sel$select()

  pilot <- sel$lambda2_pilot
  expect_true(pilot$ran)
  expect_identical(pilot$grid, c(0, 1, 10))        # sorted, de-duplicated
  expect_length(pilot$counts, 3)
  expect_equal(dim(pilot$curve_majority), c(3, 3))  # G x T_pilot
  expect_equal(dim(pilot$curve_soft), c(3, 3))
  expect_identical(pilot$K_pilot, 5)
  expect_identical(pilot$T_pilot, 3)
  expect_identical(pilot$L_0, d$p)
  expect_true(pilot$index %in% 1:3)
  expect_identical(pilot$grid[pilot$index], pilot$lambda2)
  expect_identical(sel$lambda2_used, pilot$lambda2)
  expect_identical(pilot$counts[pilot$index], pilot$max_selected)
  # The pilot draws its dummies under a seed tag of its own
  expect_identical(pilot$dummy_tag, 2^30)
  # The pilot is timed as the lambda_2 resolution step.
  expect_true(sel$lambda2_cv_diagnostics$ran)

  # With shared dummies the pilot takes the first block of the main strategy
  shared <- TRexTikhonovSelector$new(
    d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
    tik_control = trex_tikhonov_control(
      tikhonov_K = .group_mean_k(d$p), prior_groups = d$groups,
      pilot_grid = c(0, 1), pilot_K = 5, pilot_T = 3, pilot_shared_dummies = TRUE),
    control = .tik_base(K = 5, max_dummy_multiplier = 2))
  shared$select()
  expect_identical(shared$lambda2_pilot$dummy_tag, 1)
})


test_that("prior-groups end-to-end recovers the planted signal", {
  d <- .make_tik_data()

  sel <- TRexTikhonovSelector$new(
    d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
    tik_control = trex_tikhonov_control(
      tikhonov_K = .group_mean_k(d$p), lambda_2 = 1.0,
      prior_groups = d$groups),
    control = .tik_base(K = 5, max_dummy_multiplier = 2))
  sel$select()

  idx <- sel$selected_indices
  expect_true(length(intersect(idx, 1:3)) >= 2)
  expect_identical(sel$lambda2_used, 1.0)
  # A caller-supplied lambda_2: neither a tuner nor a pilot ran.
  expect_false(sel$lambda2_cv_diagnostics$ran)
  expect_false(sel$lambda2_pilot$ran)
})


test_that("lambda_2 auto-resolves via the IEN-geometry tuner", {
  d <- .make_tik_data()

  sel <- TRexTikhonovSelector$new(
    d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
    tik_control = trex_tikhonov_control(
      tikhonov_K = .group_mean_k(d$p), lambda_2 = -1,
      prior_groups = d$groups, lambda2_method = "CV_1SE_IEN_CCD",
      cv_n_lambda1 = 20, cv_n_lambda2 = 10, cv_dfmax_ratio = 0.9),
    control = .tik_base(K = 5, max_dummy_multiplier = 2))
  sel$select()

  expect_gt(sel$lambda2_used, 0)
  diag <- sel$lambda2_cv_diagnostics
  expect_true(diag$ran)
  expect_gte(diag$seconds, 0)
  expect_gte(diag$n_profile_at_edge, 0)
  expect_false(sel$lambda2_pilot$ran)
  # The support-cap ratio reached the tuner and is recorded as passed.
  expect_identical(diag$dfmax_ratio, 0.9)
  expect_true(diag$profile_at_edge || !diag$profile_at_support_cap)

  # Boundary-pick guard fields.
  expect_identical(diag$grid_size, 10L)
  expect_true(diag$grid_index >= 0L && diag$grid_index < 10L)
  expect_identical(diag$boundary_pick, diag$at_grid_edge || diag$profile_at_edge)
})


test_that("trex_gamma_to_k builds the chain Laplacian from a first-difference Gamma", {
  skip_if_not_installed("Matrix")
  p <- 5
  gamma <- diff(diag(p))
  K <- trex_gamma_to_k(gamma)
  Kd <- as.matrix(K)
  expect_equal(dim(Kd), c(p, p))
  expect_equal(Kd, t(Kd))
  expect_equal(diag(Kd), c(1, 2, 2, 2, 1))
  expect_equal(Kd[1, 2], -1)
  # A sparse dgCMatrix Gamma is accepted too.
  K2 <- trex_gamma_to_k(Matrix::Matrix(gamma, sparse = TRUE))
  expect_equal(as.matrix(K2), Kd)
})


test_that("the full-profile knee runs one selection per grid value and keeps the contract", {
  d <- .make_tik_data()
  X_before <- d$X + 0  # a private copy: the probes hand X back unchanged
  tik <- trex_tikhonov_control(tikhonov_K = .group_mean_k(d$p),
                               prior_groups = d$groups)
  ctrl <- .tik_base(K = 5, max_dummy_multiplier = 2)

  prof <- trex_tikhonov_lambda2_power_rule(d$X, d$y, tFDR = 0.2, tik_control = tik,
                                           control = ctrl, seed = 42,
                                           grid = c(10, 0, 1, 10))
  expect_identical(prof$grid, c(0, 1, 10))          # sorted, de-duplicated
  expect_length(prof$probes, 3)
  expect_identical(prof$counts, vapply(prof$probes, function(p) length(p$selected), 1))
  expect_true(prof$index %in% 1:3)
  expect_identical(prof$lambda2, prof$grid[prof$index])
  expect_identical(prof$selected, prof$probes[[prof$index]]$selected)
  expect_identical(prof$max_selected, max(prof$counts))
  expect_true(prof$seconds > 0)
  expect_equal(d$X, X_before)

  # The knee: the smallest grid value within one selection of the largest.
  knee <- which(prof$counts >= max(prof$counts) - 1)[1]
  expect_identical(prof$index, knee)

  # Tolerance 0 reads the largest selection, ties toward the smallest lambda_2.
  argmax <- trex_tikhonov_lambda2_power_rule(d$X, d$y, tFDR = 0.2, tik_control = tik,
                                             control = ctrl, seed = 42,
                                             grid = c(0, 1, 10), knee_tolerance = 0)
  expect_identical(argmax$index, which.max(prof$counts))

  # The contract: a fixed-lambda_2 run at the pick with the same seed selects
  # the same variables.
  tik_fixed <- tik
  tik_fixed$lambda_2 <- prof$lambda2
  sel <- TRexTikhonovSelector$new(d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
                                  tik_control = tik_fixed, control = ctrl)
  sel$select()
  expect_identical(sort(sel$selected_indices), prof$selected)

  # Gates.
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y, tik_control = tik, seed = -1),
               "seed")
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y, tik_control = tik, grid = c(-1, 1)),
               "grid")
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y, tik_control = tik,
                                                knee_tolerance = -1),
               "knee_tolerance")
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y), "tik_control")
})


test_that("the default run on a chain Laplacian selects the planted block", {
  # The default control: the homotopy solver and the pilot rule over its default grid up to
  # 3000. The pick stays below the top of the grid.
  set.seed(2718)
  n <- 120
  p <- 60
  support <- 21:30
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, support]) + rnorm(n))
  K <- crossprod(diff(diag(p)))

  sel <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                                  tik_control = trex_tikhonov_control(tikhonov_K = K),
                                  control = .tik_base(K = 10))
  sel$select()

  expect_gte(sum(support %in% sel$selected_indices), 8)
  expect_lte(length(setdiff(sel$selected_indices, support)), 4)
  expect_true(sel$lambda2_pilot$ran)
  expect_lt(sel$lambda2_used, 1000)   # not the top of the pilot grid
})


# =============================================================================
# Checks at construction
# =============================================================================

test_that("the Tikhonov selector checks K, lambda_2 and the CV sizes at construction", {
  set.seed(31)
  n <- 40; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  X0 <- X + 0
  build <- function(...) {
    TRexTikhonovSelector$new(X, y, tFDR = 0.2, seed = 2, verbose = FALSE,
                             tik_control = trex_tikhonov_control(...),
                             control = .tik_base(K = 5))
  }
  # The content of K, under the name of the selector
  K_asym <- diag(p); K_asym[1, 2] <- 0.5
  K_neg <- diag(p); K_neg[p, p] <- -1
  expect_error(build(tikhonov_K = K_asym, lambda_2 = 1), "TRexTikhonovSelector.*not symmetric")
  expect_error(build(tikhonov_K = K_neg, lambda_2 = 1), "TRexTikhonovSelector.*positive semi-definite")
  # lambda_2 must be finite
  for (bad in c(NaN, Inf, -Inf)) {
    expect_error(build(tikhonov_K = diag(p), lambda_2 = bad), "lambda_2 must be finite")
  }
  # The folds of a cross-validation rule
  expect_error(build(tikhonov_K = diag(p), lambda2_method = "CV_1SE_IEN_CCD", cv_n_folds = 1),
               "cv_n_folds must be in \\[2, n\\]")
  expect_error(build(tikhonov_K = diag(p), lambda2_method = "CV_1SE_IEN_CCD", cv_n_folds = n + 1),
               "cv_n_folds must be in \\[2, n\\]")
  # A fixed lambda_2 does not use the folds
  expect_no_error(build(tikhonov_K = diag(p), lambda_2 = 1, cv_n_folds = 1)$restore())
  # A refused control leaves X alone
  expect_equal(X, X0, tolerance = 1e-12)
})

test_that("GROUP_CENTERING refuses a supplied K that is not the centering matrix", {
  set.seed(32)
  n <- 40; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  groups <- rep(1:2, each = 5)
  run <- function(K) {
    sel <- TRexTikhonovSelector$new(X, y, tFDR = 0.3, seed = 2, verbose = FALSE,
                                    tik_control = trex_tikhonov_control(
                                      tikhonov_K = K, penalty_structure = "GROUP_CENTERING",
                                      prior_groups = groups, lambda_2 = 1),
                                    control = .tik_base(K = 5))
    sel$select()
    sel$selected_indices
  }
  expect_error(run(diag(p)), "must be the centering matrix")
  # The right matrix and no matrix give the same run
  K_gc <- matrix(0, p, p)
  for (g in 1:2) {
    idx <- which(groups == g)
    K_gc[idx, idx] <- diag(5) - 1 / 5
  }
  expect_identical(run(K_gc), run(NULL))
})

test_that("the power rules and the pilot rule have one default grid, the one of C++", {
  grid <- trex_default_lambda2_grid()
  expect_identical(grid, c(0, 1, 3, 10, 30, 100, 300, 1000, 3000))
  expect_identical(trex_tikhonov_control(tikhonov_K = diag(2))$pilot_grid, grid)
  # R holds no copy of the numbers: the three defaults are the call of the function
  call <- "trex_default_lambda2_grid()"
  expect_identical(deparse(formals(trex_tikhonov_control)$pilot_grid), call)
  expect_identical(deparse(formals(trex_tikhonov_lambda2_power_rule)$grid), call)
  expect_identical(deparse(formals(trex_gvs_lambda2_power_rule)$grid), call)
})

test_that("a lambda_2 grid means the same shrinkage under ZSCORE as under L2", {
  d <- .make_tik_data()
  n <- nrow(d$X)
  grid <- c(0, 1, 10)
  pilot <- function(mode) {
    sel <- TRexTikhonovSelector$new(d$X + 0, d$y, tFDR = 0.2, seed = 5, verbose = FALSE,
                                    tik_control = trex_tikhonov_control(
                                      penalty_structure = "GROUP_CENTERING", prior_groups = d$groups,
                                      pilot_grid = grid, pilot_K = 5),
                                    control = .tik_base(K = 10, scaling_mode = mode))
    sel$select()
    sel
  }
  l2 <- pilot("L2")
  z <- pilot("ZSCORE")
  expect_identical(l2$lambda2_pilot$grid_scale, 1)
  expect_identical(z$lambda2_pilot$grid_scale, n - 1)
  expect_equal(l2$lambda2_pilot$grid, grid)
  expect_equal(z$lambda2_pilot$grid, grid * (n - 1))
  expect_identical(z$lambda2_pilot$counts, l2$lambda2_pilot$counts)
  expect_identical(z$lambda2_pilot$index, l2$lambda2_pilot$index)
  expect_equal(z$lambda2_used, l2$lambda2_used * (n - 1))
  expect_identical(z$selected_indices, l2$selected_indices)

  rule <- function(mode) {
    trex_tikhonov_lambda2_power_rule(
      d$X + 0, d$y, tFDR = 0.2, seed = 5, grid = grid,
      tik_control = trex_tikhonov_control(penalty_structure = "GROUP_CENTERING",
                                          prior_groups = d$groups),
      control = .tik_base(K = 10, scaling_mode = mode))
  }
  a <- rule("L2")
  b <- rule("ZSCORE")
  expect_identical(a$grid_scale, 1)
  expect_identical(b$grid_scale, n - 1)
  expect_equal(b$grid, grid * (n - 1))
  expect_identical(b$counts, a$counts)
  expect_identical(b$selected, a$selected)
})

# =============================================================================
# What the C++ class has: the pilot alone, the clusters, the centering matrix
# =============================================================================

# The pilot rule on the grouped design with a small pilot
.pilot_selector <- function(d, X = d$X + 0, ...) {
  TRexTikhonovSelector$new(X, d$y, tFDR = 0.2, seed = 5, verbose = FALSE,
                           tik_control = trex_tikhonov_control(
                             penalty_structure = "GROUP_CENTERING", prior_groups = d$groups,
                             pilot_grid = c(0, 1, 10), pilot_K = 5, ...),
                           control = .tik_base(K = 10))
}

test_that("run_pilot_only() gives the pilot of a selection and spends the selector", {
  d <- .make_tik_data()
  full <- .pilot_selector(d)
  full$select()

  X <- d$X + 0
  only <- .pilot_selector(d, X = X)
  record <- only$run_pilot_only()
  expect_true(record$ran)
  expect_identical(record$grid, c(0, 1, 10))
  expect_identical(record$counts, full$lambda2_pilot$counts)
  expect_identical(record$lambda2, full$lambda2_pilot$lambda2)
  expect_identical(only$lambda2_used, full$lambda2_used)
  expect_identical(only$lambda2_pilot, record)     # the field holds the record
  expect_equal(X, d$X, tolerance = 1e-12)          # X is given back

  # A selector runs once
  expect_error(only$select(), "has run already")
  expect_error(only$run_pilot_only(), "has run already")
  expect_error(full$run_pilot_only(), "has run already")
})

test_that("run_pilot_only() needs the pilot rule", {
  d <- .make_tik_data()
  fixed <- .pilot_selector(d, lambda_2 = 1)
  expect_error(fixed$run_pilot_only(), "PILOT_KNEE")
  # The refused call did not spend the selector
  expect_no_error(fixed$select())
  expect_identical(fixed$lambda2_used, 1)
})

test_that("the clusters of the dummies are reported", {
  d <- .make_tik_data()
  empty <- list(groups = integer(0), max_clusters = 0L, hc_method_used = "")

  before <- .pilot_selector(d, lambda_2 = 1)
  expect_identical(before$dummy_clusters, empty)
  before$restore()

  prior <- .pilot_selector(d, lambda_2 = 1)
  prior$select()
  expect_identical(prior$dummy_clusters,
                   list(groups = as.integer(d$groups), max_clusters = 8L, hc_method_used = "prior"))

  found <- function(...) {
    sel <- TRexTikhonovSelector$new(d$X + 0, d$y, tFDR = 0.2, seed = 5, verbose = FALSE,
                                    tik_control = trex_tikhonov_control(
                                      tikhonov_K = diag(d$p), lambda_2 = 1, ...),
                                    control = .tik_base(K = 5))
    sel$select()
    sel$dummy_clusters
  }
  single <- found()
  expect_identical(single$hc_method_used, "single")
  expect_length(single$groups, d$p)
  expect_identical(sort(unique(single$groups)), seq_len(single$max_clusters))
  expect_identical(found(hc_linkage = "complete")$hc_method_used, "complete")
  expect_identical(found(cluster_dummies = FALSE), empty)

  # The pilot alone reports them as well
  only <- .pilot_selector(d)
  only$run_pilot_only()
  expect_identical(only$dummy_clusters$hc_method_used, "prior")
  expect_error(only$dummy_clusters <- empty, "read-only")
})

test_that("prior_groups is checked in the words of R", {
  d <- .make_tik_data()
  numbers <- "prior_groups must number the groups 1 to M without gaps"
  for (bad in list(c(0, 0, 1, 1), c(1, 1, 3, 3), c(1, 1.5, 2, 2), c(1, NA, 2, 2), "a")) {
    expect_error(trex_tikhonov_control(tikhonov_K = diag(4), prior_groups = bad), numbers)
  }
  expect_no_error(trex_tikhonov_control(tikhonov_K = diag(4), prior_groups = c(2, 1, 1, 2)))

  short <- trex_tikhonov_control(tikhonov_K = diag(d$p), prior_groups = rep(1, d$p - 1), lambda_2 = 1)
  len <- "tik_control\\$prior_groups must have one entry per variable \\(p = 40\\). Got: 39."
  expect_error(TRexTikhonovSelector$new(d$X, d$y, verbose = FALSE, tik_control = short), len)
  expect_error(trex_tikhonov_lambda2_power_rule(d$X, d$y, tik_control = short, grid = c(0, 1)), len)
})

test_that("trex_build_group_centering_k builds the matrix of GROUP_CENTERING", {
  skip_if_not_installed("Matrix")
  K <- trex_build_group_centering_k(c(1, 1, 1, 2, 2, 3))
  expect_s4_class(K, "dgCMatrix")
  ref <- matrix(0, 6, 6)
  ref[1:3, 1:3] <- diag(3) - 1 / 3
  ref[4:5, 4:5] <- diag(2) - 1 / 2
  expect_equal(unname(as.matrix(K)), ref)             # the singleton contributes nothing
  # The order of the variables is free
  expect_equal(unname(as.matrix(trex_build_group_centering_k(c(2, 1, 2)))),
               rbind(c(0.5, 0, -0.5), 0, c(-0.5, 0, 0.5)))
  numbers <- "prior_groups must number the groups 1 to M without gaps"
  expect_error(trex_build_group_centering_k(c(0, 1)), numbers)
  expect_error(trex_build_group_centering_k(c(1, 3)), numbers)

  # GROUP_CENTERING takes the built matrix and runs as without one
  d <- .make_tik_data()
  run <- function(K) {
    sel <- TRexTikhonovSelector$new(d$X + 0, d$y, tFDR = 0.2, seed = 5, verbose = FALSE,
                                    tik_control = trex_tikhonov_control(
                                      tikhonov_K = K, penalty_structure = "GROUP_CENTERING",
                                      prior_groups = d$groups, lambda_2 = 1),
                                    control = .tik_base(K = 5))
    sel$select()
    sel$selected_indices
  }
  expect_identical(run(trex_build_group_centering_k(d$groups)), run(NULL))
})

test_that("the power rule takes a memory-mapped design", {
  d <- .make_tik_data()
  n <- nrow(d$X)
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f), add = TRUE)
  mem_matrix_from(d$X, f)
  tik <- trex_tikhonov_control(tikhonov_K = .group_mean_k(d$p), prior_groups = d$groups)
  rule <- function(X, ...) {
    trex_tikhonov_lambda2_power_rule(X, d$y, tFDR = 0.2, tik_control = tik, seed = 42,
                                     grid = c(0, 1, 10), control = .tik_base(K = 5, ...))
  }
  in_memory <- rule(d$X + 0)
  on_file <- rule(mem_matrix(f, rows = n, cols = d$p, mode = "readwrite"))
  for (field in c("lambda2", "index", "selected", "grid", "counts")) {
    expect_identical(on_file[[field]], in_memory[[field]], info = field)
  }
  # Every probe hands the design back, so the file holds it as before
  expect_equal(as.matrix(mem_matrix(f, rows = n, cols = d$p, mode = "readonly")), d$X,
               tolerance = 1e-12)

  # A read-only file is refused where the probes write, and serves under NONE
  ro <- mem_matrix(f, rows = n, cols = d$p, mode = "readonly")
  expect_error(rule(ro), "ReadOnly")
  pre <- trex_normalize_design(d$X, d$y)
  g <- tempfile(fileext = ".bin")
  on.exit(unlink(g), add = TRUE)
  mem_matrix_from(pre$X, g)
  md5 <- unname(tools::md5sum(g))
  none <- function(X) {
    trex_tikhonov_lambda2_power_rule(X, pre$y, tFDR = 0.2, tik_control = tik, seed = 42,
                                     grid = c(0, 1, 10),
                                     control = .tik_base(K = 5, scaling_mode = "NONE"))
  }
  expect_identical(none(mem_matrix(g, rows = n, cols = d$p, mode = "readonly"))$selected,
                   none(pre$X + 0)$selected)
  expect_identical(unname(tools::md5sum(g)), md5)     # nothing was written
})
