# =======================================================================================
# Helpers for the tests that run the deprecated coordinate-descent solver family
# =======================================================================================

# Sweep budget of the selector runs of the family. At the default budget such a run can end with the
# convergence error on some draws of the dummies, and the dummies differ per platform.
cd_test_sweeps <- 2e5

# Runs code without the deprecation warning of the family.
quiet_cd <- function(expr) suppressWarnings(expr, classes = "deprecatedWarning")

# Messages of the deprecation warnings that code raises.
cd_deprecations <- function(expr) {
  seen <- character()
  withCallingHandlers(expr, deprecatedWarning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  seen
}
