# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(QHScrnomo)

# Run a basic model
dd <- datadist(prostate.dat)
options(datadist = "dd")
mod_cph <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX, data = prostate.dat)
mod_crr <- crr.fit(mod_cph, cencode = 0, failcode = 1)

test_check("QHScrnomo")
