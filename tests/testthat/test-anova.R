# Make sure that it returns an rms anova
test_that("It is an anova.rms object", {
  expect_equal(methods::is(anova(prostate.crr)), "anova.rms")
})
