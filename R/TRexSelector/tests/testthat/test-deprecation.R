# =======================================================================================
# The deprecation warning of the package
# =======================================================================================

# The function is internal, and a test run against the attached package sees exported names alone.
.warn_deprecated <- TRexSelector:::.warn_deprecated

test_that("the deprecation warning carries its class, the old name and the replacement", {
  w <- tryCatch(.warn_deprecated("old_name", "new_name"), warning = function(w) w)
  expect_s3_class(w, "deprecatedWarning")
  expect_identical(conditionMessage(w), "old_name is deprecated. Use new_name.")
  expect_identical(w$old, "old_name")
  expect_identical(w$new, "new_name")
  expect_identical(w$package, "TRexSelector")
  expect_null(conditionCall(w))
})


test_that("the deprecation warning is silenced by its class alone", {
  expect_no_warning(suppressWarnings(.warn_deprecated("old_name", "new_name"),
                                     classes = "deprecatedWarning"))
  expect_warning(suppressWarnings(.warn_deprecated("old_name", "new_name"),
                                  classes = "someOtherWarning"),
                 class = "deprecatedWarning")
})
