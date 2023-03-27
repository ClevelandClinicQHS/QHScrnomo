##' generate anova table for crr
##'
##' generate anova table for competing risks regression model
##' @title anova table for competing risks regression
##' @usage \S3method{anova}{cmprsk}(object, \dots)
##' @param object a competing risks regression model object
##' built from funciton \code{\link[QHScrnomo]{crr.fit}}
##' @param ... other arguments
##' @return anova table in matrix
##' @importFrom stats anova
##' @export
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
##' ## anova test
##' anova(prostate.crr)
##'
anova.cmprsk <- function(object, ...) {
    if ("cmprsk" %nin% class(object)) {
        stop("the model is not fitter by crr.fit!!")
    }
    cph.f <- object$cph.f
    coe.name <- names(cph.f$coefficients)
    cph.f$coefficients <- object$coef
    cph.f$var <- object$var
    anova(cph.f)
}
