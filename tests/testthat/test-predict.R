# Test that I get the number of original rows
test_that("The number of predictions equals development row count when no newdata is supplied", {
  expect_equal(length(predict(prostate.crr, time = 60)), nrow(prostate.dat))
})

# Ensure I get a warning with multiple time points
test_that("A warning is generated for supplying multiple time points", {
  expect_warning(predict(prostate.crr, time = c(50, 60)))
})

# Ensure I get an error when specifying too large of a time point
test_that("An error is generated for supplying a time point too large", {
  expect_error(predict(prostate.crr, time = Inf))
})

# Ensure I get an error when specifying too small of a time point
test_that("An error is generated for supplying a time point too small", {
  expect_error(predict(prostate.crr, time = -Inf))
})
