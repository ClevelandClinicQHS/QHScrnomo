# Errors
test_that("An error is thrown when a time point is not specified", {
  expect_error(nomogram.crr(prostate.crr))
})
test_that("An error is thrown when trying to pass a rms::cph object", {
  expect_error(nomogram.crr(prostate.crr$cph.f, failtime = 50))
})

# Return value
test_that("A nomogram object is returned", {
  expect_equal(class(nomogram.crr(prostate.crr, failtime = 50)), "nomogram")
})
