# =======================================================================================
# CCD solver family (TCDLASSO / TCDENET / TCDIENET) through the R selector surface
# =======================================================================================

# Shared planted design: strong signal on variables 1:5.
.make_cd_data <- function(n = 100, p = 20, data_seed = 42) {
  set.seed(data_seed)
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(5, 5), rep(0, p - 5))
  y <- as.vector(X %*% beta + rnorm(n))
  list(X = X, y = y, support = 1:5, p = p)
}


test_that("trex_control validates the CD knobs", {

  # certify/probe must be set together (the C++ side applies the pair only
  # when both are > 0 — a half-set pair must fail loudly, not silently).
  expect_error(
    trex_control(cd_tol_certify = 1e-6),
    "cd_tol_certify and cd_tol_probe must be set together"
  )
  expect_error(
    trex_control(cd_tol_probe = 1e-6),
    "cd_tol_certify and cd_tol_probe must be set together"
  )
  expect_silent(trex_control(cd_tol_certify = 1e-6, cd_tol_probe = 1e-7))

  expect_error(trex_control(cd_gram_cap = -1), "cd_gram_cap must be >= 0")
  expect_error(trex_control(cd_max_sweeps = -1), "cd_max_sweeps must be >= 0")

  # The stagnation-stop default is NULL on the R side and resolved by the C++
  # core's solver rule (isGreedySolver), the one source for the three
  # surfaces; an explicit value is carried as given.
  for (s in c("TCDLASSO", "TCDENET", "TCDIENET", "TCDTIKNET", "TIKNET",
              "TENET_AUG", "TOMP")) {
    expect_null(trex_control(solver = s)$tloop_stagnation_stop)
  }
  expect_false(trex_control(solver = "TOMP", tloop_stagnation_stop = FALSE)$tloop_stagnation_stop)
  expect_true(trex_control(tloop_stagnation_stop = TRUE)$tloop_stagnation_stop)

  # Certified crossing localization is the CD default; the flag is carried.
  expect_true(trex_control()$cd_certified_localization)
  expect_false(trex_control(cd_certified_localization = FALSE)$cd_certified_localization)
})


test_that("TCDLASSO and TCDENET run end-to-end and recover the planted support", {
  d <- .make_cd_data()

  for (s in c("TCDLASSO", "TCDENET")) {
    selector <- quiet_cd(TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 42,
                                          verbose = FALSE,
                                          control = trex_control(solver = s, K = 5,
                                                                 cd_max_sweeps = cd_test_sweeps)))
    selector$select()

    expect_length(selector$selected_var, d$p)
    idx <- selector$selected_indices
    expect_true(all(d$support %in% idx),
                info = paste0(s, ": planted support not recovered. Selected: ",
                              paste(idx, collapse = ", ")))
    expect_lte(length(setdiff(idx, d$support)), 3)
  }
})


test_that("CD selection is seed-reproducible", {
  d <- .make_cd_data()

  run_once <- function() {
    selector <- TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 7,
                                 verbose = FALSE,
                                 control = trex_control(solver = "TCDENET", K = 5,
                                                        cd_max_sweeps = cd_test_sweeps))
    selector$select()
    selector$selected_var
  }
  expect_identical(quiet_cd(run_once()), quiet_cd(run_once()))
})


test_that("CD knobs pass through the control parser", {
  d <- .make_cd_data()

  ctrl <- trex_control(solver = "TCDENET", K = 3,
                       cd_lambda_rel_tol = 1e-4,
                       cd_tol_certify = 1e-6, cd_tol_probe = 1e-7,
                       cd_gram_cap = 100, cd_max_sweeps = cd_test_sweeps)
  selector <- quiet_cd(TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 42,
                                        verbose = FALSE, control = ctrl))
  selector$select()
  expect_true(all(d$support %in% selector$selected_indices))
})


test_that("GVS-only solvers are rejected by the plain selector with guidance", {
  d <- .make_cd_data(n = 50, p = 10)

  run_with <- function(s) {
    selector <- TRexSelector$new(d$X, d$y, tFDR = 0.2, seed = 42,
                                 verbose = FALSE,
                                 control = trex_control(solver = s, K = 3))
    selector$select()
  }
  expect_error(quiet_cd(run_with("TCDIENET")), "needs a group assignment")
  expect_error(run_with("TIENET"), "GVS-only solver")
  expect_error(run_with("TIENET_AUG"), "GVS-only solver")
  # The Tikhonov family needs a penalty matrix: TRexTikhonovSelector only.
  expect_error(quiet_cd(run_with("TCDTIKNET")), "Tikhonov penalty matrix")
  expect_error(run_with("TIKNET"), "Tikhonov penalty matrix")
})


test_that("both CD crossing-localization modes run and recover the support", {
  d <- .make_cd_data()

  for (certified in c(TRUE, FALSE)) {
    selector <- quiet_cd(TRexSelector$new(
      d$X, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
      control = trex_control(solver = "TCDLASSO", K = 5, cd_max_sweeps = cd_test_sweeps,
                             cd_certified_localization = certified)))
    selector$select()
    expect_true(all(d$support %in% selector$selected_indices),
                info = paste("cd_certified_localization =", certified))
  }
})


test_that("TENET_AUG runs in the plain selector and matches TENET", {
  d <- .make_cd_data()

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


test_that("GVS runs the CCD family: en_solver = TCDENET and IEN solver choices", {
  set.seed(42)
  n <- 80
  p <- 20
  # Grouped design: 4 blocks of 5 correlated columns, signal on block 1.
  z <- matrix(rnorm(n * 4), n, 4)
  X <- matrix(0, n, p)
  for (g in 1:4) {
    for (j in 1:5) {
      X[, (g - 1) * 5 + j] <- 0.85 * z[, g] + sqrt(1 - 0.85^2) * rnorm(n)
    }
  }
  y <- as.vector(rowSums(X[, 1:3]) * 3 + rnorm(n))
  # A live selector claims and normalizes its X buffer in place
  # (shared-buffer claim guard), so every constructor gets its own copy.
  fresh <- function() X + 0

  # EN track, CCD elastic net.
  sel_en <- quiet_cd(TRexGVSSelector$new(
    fresh(), y, tFDR = 0.2, seed = 42, verbose = FALSE,
    gvs_control = trex_gvs_control(gvs_type = "EN", en_solver = "TCDENET",
                                   lambda_2 = 1.0),
    control = trex_control(K = 5, cd_max_sweeps = cd_test_sweeps)
  ))
  sel_en$select()
  expect_length(sel_en$selected_var, p)

  # IEN track, CCD informed elastic net via control$solver.
  sel_ien <- quiet_cd(TRexGVSSelector$new(
    fresh(), y, tFDR = 0.2, seed = 42, verbose = FALSE,
    gvs_control = trex_gvs_control(gvs_type = "IEN", lambda_2 = 1.0),
    control = trex_control(solver = "TCDIENET", K = 5, cd_max_sweeps = cd_test_sweeps)
  ))
  sel_ien$select()
  expect_length(sel_ien$selected_var, p)

  # EN track warns when control$solver conflicts with the en_solver axis.
  expect_warning(
    TRexGVSSelector$new(
      fresh(), y, tFDR = 0.2, seed = 42, verbose = FALSE,
      gvs_control = trex_gvs_control(gvs_type = "EN", lambda_2 = 1.0),
      control = trex_control(solver = "TCDENET", K = 5)
    ),
    "derives solver"
  )

  # IEN track warns on a non-IEN-family control$solver and falls back.
  expect_warning(
    TRexGVSSelector$new(
      fresh(), y, tFDR = 0.2, seed = 42, verbose = FALSE,
      gvs_control = trex_gvs_control(gvs_type = "IEN", lambda_2 = 1.0),
      control = trex_control(solver = "TOMP", K = 5)
    ),
    "not an IEN-family solver"
  )
})


# =======================================================================================
# Deprecation of the family
# =======================================================================================

.cd_name_message <- function(argument, old, new) {
  paste0(argument, " = \"", old, "\" is deprecated and will be removed in TRexSelector 2.1.0. ",
         "Use \"", new, "\", the homotopy solver of the same problem.")
}


test_that("a selector warns once for a coordinate-descent solver name", {
  d <- .make_cd_data(n = 50, p = 10)

  # One warning at construction that names the replacement, none from the experiments
  w <- cd_deprecations(
    sel <- TRexSelector$new(d$X + 0, d$y, tFDR = 0.2, seed = 42, verbose = FALSE,
                            control = trex_control(solver = "TCDLASSO", K = 3, cd_max_sweeps = cd_test_sweeps)))
  expect_identical(w, .cd_name_message("solver", "TCDLASSO", "TLASSO"))
  expect_length(cd_deprecations(sel$select()), 0)

  # Every selector that takes the base control
  makers <- list(
    function(ctrl) TRexSelector$new(d$X + 0, d$y, verbose = FALSE, control = ctrl),
    function(ctrl) TRexDASelector$new(d$X + 0, d$y, verbose = FALSE, control = ctrl),
    function(ctrl) TRexScreeningSelector$new(d$X + 0, d$y, verbose = FALSE, control = ctrl),
    function(ctrl) TRexBiobankScreeningSelector$new(d$X + 0, d$y, verbose = FALSE, control = ctrl)
  )
  for (make in makers) {
    expect_identical(cd_deprecations(make(trex_control(solver = "TCDENET", K = 3))),
                     .cd_name_message("solver", "TCDENET", "TENET"))
    expect_length(cd_deprecations(make(trex_control(solver = "TENET", K = 3))), 0)
    expect_length(cd_deprecations(make(trex_control(K = 3))), 0)
  }
})


test_that("the grouped and the Tikhonov selector warn for their coordinate-descent solvers", {
  d <- .make_cd_data(n = 50, p = 10)
  groups <- rep(1:5, each = 2)
  tik <- trex_tikhonov_control(tikhonov_K = diag(d$p), lambda_2 = 0.5)

  # Elastic-net track: the name comes from en_solver
  w <- cd_deprecations(TRexGVSSelector$new(
    d$X + 0, d$y, verbose = FALSE,
    gvs_control = trex_gvs_control(gvs_type = "EN", en_solver = "TCDENET", lambda_2 = 1)))
  expect_identical(w, .cd_name_message("en_solver", "TCDENET", "TENET"))

  # Informed elastic-net track: the name comes from the base control
  w <- cd_deprecations(TRexGVSSelector$new(
    d$X + 0, d$y, verbose = FALSE,
    gvs_control = trex_gvs_control(gvs_type = "IEN", lambda_2 = 1, groups = groups),
    control = trex_control(solver = "TCDIENET")))
  expect_identical(w, .cd_name_message("solver", "TCDIENET", "TIENET"))

  w <- cd_deprecations(TRexTikhonovSelector$new(
    d$X + 0, d$y, verbose = FALSE, tik_control = tik,
    control = trex_control(solver = "TCDTIKNET")))
  expect_identical(w, .cd_name_message("solver", "TCDTIKNET", "TIKNET"))

  w <- cd_deprecations(trex_tikhonov_lambda2_power_rule(
    d$X + 0, d$y, tik_control = tik, seed = 1, grid = c(0, 1),
    control = trex_control(solver = "TCDTIKNET", K = 3, cd_max_sweeps = cd_test_sweeps)))
  expect_identical(w, .cd_name_message("solver", "TCDTIKNET", "TIKNET"))

  # The defaults raise none
  expect_length(cd_deprecations(TRexGVSSelector$new(
    d$X + 0, d$y, verbose = FALSE, gvs_control = trex_gvs_control(lambda_2 = 1))), 0)
  expect_length(cd_deprecations(TRexTikhonovSelector$new(
    d$X + 0, d$y, verbose = FALSE, tik_control = tik)), 0)
})


test_that("constructing a coordinate-descent solver class warns and names the replacement", {
  set.seed(3)
  n <- 20
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  D <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  opts <- tsolver_options(lambda2 = 0.1)
  groups <- c(1, 1, 2, 2)

  makers <- list(
    list("TCDLASSO", "TLASSO", function() TCDLASSO_Solver$new(X, D, y)),
    list("TCDENET", "TENET", function() TCDENET_Solver$new(X, D, y, options = opts)),
    list("TCDIENET", "TIENET", function() TCDIENET_Solver$new(X, D, y, groups = groups, options = opts)),
    list("TCDTIKNET", "TIKNET", function() TCDTIKNET_Solver$new(X, D, y, tikhonov_K = diag(p), options = opts))
  )
  for (m in makers) {
    w <- cd_deprecations(solver <- m[[3]]())
    solver$restore()
    expect_identical(w, paste0(m[[1]], "_Solver is deprecated and will be removed in TRexSelector 2.1.0. ",
                               "Use ", m[[2]], "_Solver, the homotopy solver of the same problem."))
  }

  # The condition is the one .Deprecated() raises
  w <- tryCatch(TCDLASSO_Solver$new(X, D, y), deprecatedWarning = function(w) w)
  expect_s3_class(w, "deprecatedWarning")
  expect_identical(w$old, "TCDLASSO_Solver")
  expect_identical(w$new, "TLASSO_Solver")
  expect_identical(w$package, "TRexSelector")

  # The homotopy solvers raise none
  replacements <- list(
    function() TLASSO_Solver$new(X, D, y),
    function() TENET_Solver$new(X, D, y, options = opts),
    function() TIENET_Solver$new(X, D, y, groups = groups, options = opts),
    function() TIKNET_Solver$new(X, D, y, tikhonov_K = diag(p), options = opts)
  )
  for (make in replacements) {
    expect_length(cd_deprecations(solver <- make()), 0)
    solver$restore()
  }
})
