# Ensure the label makes its way to new label set
test_that("The set label exists in the new object", {
  expect_equal("Treatment options" %in% Newlabels(prostate.crr, labels = c(TX = "Treatment options"))$cph.f$Design$label, TRUE)
})
