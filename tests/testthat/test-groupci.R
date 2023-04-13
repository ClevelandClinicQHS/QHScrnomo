# Errors for erroneous inputs
test_that("An error is thrown for not supplying a time point of interest", {
  expect_error(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD))
})
test_that("An error is thrown for specifying invalid event codes", {
  expect_error(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, u = 120, cencode = 10))
})
test_that("An error is thrown for having unequal input lengths", {
  expect_error(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD[-1], u = 120))
})
test_that("An error is thrown when too many groups are being defined (lack of events in groups)", {
  expect_error(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, u = 20, g = 100))
})

# Warnings
test_that("A warning is generated for supplying multiple time points", {
  expect_warning(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, u = c(50, 120)))
})
test_that("A warning is generated for supplying missing data", {
  expect_warning(groupci(ifelse(prostate.dat$PSA > 50, NA_real_, prostate.dat$PSA), prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, u = 20))
})

# Output
test_that("The result is 1 row per group", {
  expect_equal(nrow(groupci(prostate.dat$PSA, prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, u = 20, g = 5)), 5)
})

