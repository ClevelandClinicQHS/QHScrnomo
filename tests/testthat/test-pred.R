# Ensure we get the correct number of rows
test_that("We get a single row when no groups are used", {
  expect_equal(nrow(pred.ci(cum, tm1 = 60, failcode = 1)), 1)
})
test_that("We get a row per group when a grouping variable is used", {
  expect_equal(nrow(pred.ci(cum2, tm1 = 60, failcode = 1)), length(unique(prostate.dat$TX)))
})

# Ensure I get a warning with multiple time points
test_that("A warning is generated for supplying multiple time points", {
  expect_warning(pred.ci(cum, tm1 = c(50, 60)))
})

# Ensure I get an error when specifying too large of a time point
test_that("An error is generated for supplying a time point too large", {
  expect_error(pred.ci(cum, tm1 = Inf))
})

# Ensure I get an error when specifying too small of a time point
test_that("An error is generated for supplying a time point too small", {
  expect_error(pred.ci(cum, tm1 = -1))
})
