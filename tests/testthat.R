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
prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
                    BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
                    RACE_AA, data = prostate.dat,
                  x = FALSE, y = TRUE, surv = TRUE, time.inc = 144)
prostate.crr <- crr.fit(prostate.f, cencode = 0, failcode = 1)

# Fit some cumulative incidence curves
cum <- cmprsk::cuminc(prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, cencode = 0)
cum2 <- cmprsk::cuminc(prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, prostate.dat$TX, cencode = 0)

test_check("QHScrnomo")
