## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE, 
  message = FALSE,
  comment = "#>"
)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
library(QHScrnomo)
sessionInfo()

## -----------------------------------------------------------------------------
data("prostate.dat")
str(prostate.dat)

## -----------------------------------------------------------------------------
dd <- datadist(prostate.dat)
options(datadist = "dd")
prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
           BX_GLSN_CAT + CLIN_STG + rcs(AGE,3) +
           RACE_AA, data = prostate.dat,
           x = TRUE, y= TRUE, surv=TRUE, time.inc = 144)
prostate.crr <- crr.fit(prostate.f, cencode = 0, failcode = 1)
summary(prostate.crr)

## ----fig.width=7, fig.height=5------------------------------------------------
prostate.g <- Newlabels(prostate.crr,
                        c(TX = 'Treatment options', 
                          BX_GLSN_CAT = 'Biopsy Gleason Score Sum',
                          CLIN_STG = 'Clinical stage'))
nomogram.crr(prostate.g,
             failtime = 120,
             lp=FALSE,
             xfrac=0.65,
             fun.at = seq(0.2, 0.45, 0.05),
             funlabel = "Predicted 10-year cumulative incidence")

## -----------------------------------------------------------------------------
# output a math formula
sas.cmprsk(prostate.crr, time = 120)

## -----------------------------------------------------------------------------
# anova table
anova(prostate.crr)

## -----------------------------------------------------------------------------
# prediction from 10-fold cross validation
prostate.dat$preds.tenf <- tenf.crr(prostate.crr, time=120, fold = 10)
str(prostate.dat$preds.tenf)

## -----------------------------------------------------------------------------
## calculate the CRR version of concordance index
with(prostate.dat, cindex(preds.tenf,
                          ftime = TIME_EVENT,
                          fstatus =EVENT_DOD, type = "crr"))["cindex"]

## ----fig.width=5, fig.height=5------------------------------------------------
## generate the calibration curve for predicted 10-year cancer
## specific mortality
with(prostate.dat,
     groupci(
         preds.tenf, 
         ftime = TIME_EVENT,
         fstatus =EVENT_DOD, g = 5, u = 120,
         xlab = "Nomogram predicted 10-year cancerspecific mortality",
         ylab = "Observed predicted 10-year cancerspecific mortality")
)

