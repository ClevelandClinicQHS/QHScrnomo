# Ensure I get a warning with no output
test_that("A warning is generated for not generating any output", {
  expect_warning(sas.cmprsk(prostate.crr, baseonly = TRUE))
})

# Ensure I get an error when specifying too large of a time point
test_that("An error is generated for supplying a time point too large", {
  expect_error(sas.cmprsk(prostate.crr, time = Inf))
})

# Ensure I get an error when specifying too small of a time point
test_that("An error is generated for supplying a time point too small", {
  expect_error(sas.cmprsk(prostate.crr, time = 0))
})

# There is a NULL object returned
test_that("The output doesn't return any object, it just prints to the console.", {
  expect_equal(class(sas.cmprsk(prostate.crr, time = 10)), "NULL")
})
