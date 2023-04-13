# Errors
test_that("An error throws for entering an unselectable option", {
  expect_error(cindex(prostate.dat$PSA, prostate.dat$EVENT_DOD, type = "linear"))
})
test_that("An error throws for having inputs of unequal length", {
  expect_error(cindex(prostate.dat$PSA, prostate.dat$EVENT_DOD[-1]))
})
test_that("An error throws for having inputs of unequal length (survival/crr mode)", {
  expect_error(cindex(prostate.dat$PSA, prostate.dat$EVENT_DOD, prostate.dat$TIME_EVENT[-1]))
})
test_that("An error is thrown for specifying invalid event code", {
  expect_error(cindex(prostate.dat$PSA, prostate.dat$EVENT_DOD, prostate.dat$TIME_EVENT, failcode = 10))
})
