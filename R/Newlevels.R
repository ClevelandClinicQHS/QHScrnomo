##' Change the Level Labels of Categorical Predictors for a Competing Risks Regression Model
##'
##' Uses the \code{\link[rms]{Newlevels}} function to change the labels predictors when constructing a nomogram.
##'
##' @usage \S3method{Newlevels}{cmprsk}(fit, levels, \dots)
##'
##' @param fit A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param levels A list of named vectors specifying the new level labels for categorical predictors.
##' @param ... Other arguments for \code{\link[rms]{Newlevels}}
##'
##' @return A new \code{\link[QHScrnomo]{crr.fit}} object with adjusted labels on the factor levels
##'
##' @note This will override \code{parms} and \code{\link[rms]{datadist}} information that were stored with the fit.
##' @author Changhong Yu. Department of
##' Quantitative Health Sciences, Cleveland Clinic
##' @seealso \code{\link[QHScrnomo]{Newlabels.cmprsk}}
##' @keywords attributes
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
##' prostate.g <- Newlevels(prostate.crr, list(TX = c('Treatment 1', 'Treatment 2', 'Treatment 3')))
##'
Newlevels.cmprsk <-
  function(
    fit, # An object fit from crr.fit
    levels, # Named list of character vector to rename categorical variables
    ...
  ) {

    # Check for fit
    if(missing(fit))
      stop("Please supply a model fit from crr.fit")

    # Check for levels
    if(missing(levels))
      stop("Please supply a named list of character vectors with labels for categorical factor levels.")

    # Check if the object was fit from crr.fit
    if(!inherits(fit, "cmprsk"))
      stop("The object is not a 'cmprsk' object fit from crr.fit")

    # Extract the original rms::cph model
    cph.f <- fit$cph.f

    # Set the new labels
    fit$cph.f <- Newlevels(cph.f, levels = levels, ... = ...)

    # Return the new fit
    fit

  }
