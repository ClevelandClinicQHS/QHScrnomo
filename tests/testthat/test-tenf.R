# Test that I get the number of original rows
test_that("The number of predictions equals development row count", {
  expect_equal(length(tenf.crr(prostate.crr, time = 60, fold = 2)), nrow(prostate.dat))
})

# Ensure I get a warning with multiple time points
test_that("A warning is generated for supplying multiple time points", {
  expect_warning(tenf.crr(prostate.crr, time = c(50, 60), fold = 2))
})

# Ensure I get a warning when supplying a time point with a linear predictor
test_that("A warning is generated for requesting predictions on the linear predictor scale but also giving a time point", {
  expect_warning(tenf.crr(prostate.crr, time = 60, fold = 2, lps = TRUE))
})

# Ensure I get an error when specifying too large of a time point
test_that("An error is generated for supplying a time point too large", {
  expect_error(tenf.crr(prostate.crr, time = Inf))
})

# Ensure I get an error when specifying too small of a time point
test_that("An error is generated for supplying a time point too small", {
  expect_error(tenf.crr(prostate.crr, time = -Inf))
})
