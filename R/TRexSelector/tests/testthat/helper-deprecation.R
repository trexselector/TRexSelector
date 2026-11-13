# =======================================================================================
# Helpers for the tests that run the deprecated coordinate-descent solver family
# =======================================================================================

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
