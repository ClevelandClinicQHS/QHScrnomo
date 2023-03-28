# Test for the assigned class of output
test_that("The model is a cmprsk object", {
  expect_equal(methods::is(mod_crr, "cmprsk"), TRUE)
})

# Test that the model was refit with the design
test_that("The design (x) was added", {
  expect_equal("x" %in% names(mod_crr$cph.f), TRUE)
})
