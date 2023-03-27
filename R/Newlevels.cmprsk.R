##'
##' ##' Change levels of categorical variable
##'
##' This method function was written for competing risks regression
##' model for facilitating to change the levels of categorical
##' predictors when construct a nomogram. It is used for the generic
##' function \code{\link{Newlevels}}
##' @title Change levels of categorical variable for a model fit
##' @param fit a model fit
##' @param levels  a list of named vectors specifying new level labels
##' for categorical predictors. This will override parms as well as
##' datadist information (if available) that were stored with the fit.
##' @param ... other arguments
##' @return returns a new model fit object with the levels adjusted.
##' @export 
##'
##' @keywords attributes
##'
##' @examples
##' 
##' data(prostate.dat)
##' dd <- datadist(prostate.dat)
##' options(datadist = "dd")
##' prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
##'            BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
##'            RACE_AA, data = prostate.dat,
##'            x = TRUE, y= TRUE, surv=TRUE,time.inc = 144)
##' prostate.crr <- crr.fit(prostate.f,cencode = 0,failcode = 1)
##' prostate.g <- Newlevels(prostate.crr, 
##'     list(TX=c('Treatment 1','Treatment 2', 'Treatment 3')))
##'
Newlevels.cmprsk <- function(fit, levels, ...) {
    cph.f <- fit$cph.f
    fit$cph.f <- Newlevels(cph.f, levels = levels, ... = ...)
    fit
}
