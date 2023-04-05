# Test for the assigned class of output
test_that("The model is a cmprsk object", {
  expect_equal(inherits(prostate.crr, "cmprsk"), TRUE)
})

# Test that the model was refit with the design
test_that("The design (x) was added", {
  expect_equal("x" %in% names(prostate.crr$cph.f), TRUE)
})

# Test that the supplied statuses must be present
test_that("The supplied 'failcode' value must be in the status column", {
  expect_error(crr.fit(prostate.f, cencode = 0, failcode = 5))
})
test_that("The supplied 'cencode' value must be in the status column", {
  expect_error(crr.fit(prostate.f, cencode = -1, failcode = 2))
})
