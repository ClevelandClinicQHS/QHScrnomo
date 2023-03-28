# Make sure that it returns an rms summary
test_that("It is an summary.rms object", {
  expect_equal(methods::is(summary(mod_crr)), "summary.rms")
})

# Ensure we can pass additional arguments
test_that("We can pass additional summary.rms arguments", {
  expect_equal(methods::is(summary(mod_crr, conf.int = 0.75, antilog = FALSE)), "summary.rms")
})
