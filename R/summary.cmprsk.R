##' generate summary information
##'
##' summarize a crr regression model
##'
##' @title summary for competing risks regression
##' @usage \S3method{summary}{cmprsk}(object,\dots)
##' @param object a crr model object from function 
##' \code{\link[QHScrnomo]{crr.fit}}
##' @param ... other parameters
##' @return  a matrix
##' @author changhong
##' @export
##'
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
##' summary(prostate.crr)
##' 
##' 

summary.cmprsk <- function(object, ...) {
    if ("cmprsk" %nin% class(object)) 
        stop("the model is not fitter by crr.fit!!")
    cph.f <- object$cph.f
    coe.name <- names(cph.f$coefficients)
    cph.f$coefficients <- object$coef
    cph.f$var <- object$var
    summary(cph.f, ...)
}
