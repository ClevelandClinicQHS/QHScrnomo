##' Summary of a Competing Risks Regression Model
##'
##' Uses the \code{\link[rms]{summary.rms}} method to construct a summary for the competing risks regression model fit from \code{\link[QHScrnomo]{crr.fit}}.
##'
##' @usage \S3method{summary}{cmprsk}(object, \dots)
##'
##' @param object A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param ... Other arguments for \code{\link[rms]{summary.rms}}
##'
##' @return A \code{\link[rms]{summary.rms}} matrix
##'
##' @note This function requires that the \code{\link{rms}} package is attached
##' @author Changhong Yu. Department of
##' Quantitative Health Sciences, Cleveland Clinic
##' @seealso \code{\link[QHScrnomo]{crr.fit}} \code{\link[rms]{summary.rms}}
##'
##' @import rms
##' @export
##' @examples
##' dd <- datadist(prostate.dat)
##' options(datadist = "dd")
##' prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
##'            BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
##'            RACE_AA, data = prostate.dat,
##'            x = TRUE, y = TRUE, surv = TRUE,time.inc = 144)
##' prostate.crr <- crr.fit(prostate.f, cencode = 0, failcode = 1)
##' summary(prostate.crr)
##'
summary.cmprsk <-
  function(
      object, # An object fit from crr.fit
      ...
    ) {

    # Check for fit
    if(missing(object))
      stop("Please supply a model fit from crr.fit")

    # Check if the object was fit from crr.fit
    if(!inherits(object, "cmprsk"))
      stop("The object is not a 'cmprsk' object fit from crr.fit")

    # Extract the original rms::cph model
    cph.f <- object$cph.f

    # Set the coefficients and covariance matrix that of the crr.fit model
    cph.f$coefficients <- object$coef
    cph.f$var <- object$var

    # Return the summary for the now 'rms' object
    summary(cph.f, ...)

  }
