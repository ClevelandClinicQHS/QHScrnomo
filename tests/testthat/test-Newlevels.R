# Ensure the levels make their way to new level set
test_that("The set factor levels are the same in the new object", {
  expect_equal(all.equal(c('Treatment 1','Treatment 2', 'Treatment 3'),  Newlevels(prostate.crr, levels = list(TX = c('Treatment 1','Treatment 2', 'Treatment 3')))$cph.f$Design$parms$TX), TRUE)
})
