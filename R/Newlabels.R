##' Change the Predictor Labels for a Competing Risks Regression Model
##'
##' Uses the \code{\link[rms]{Newlabels}} function to change the labels predictors when constructing a nomogram.
##'
##' @usage \S3method{Newlabels}{cmprsk}(fit, labels, \dots)
##'
##' @param fit A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param labels A character vector specifying the new labels for variables in a fit.
##' @param ... Other arguments for \code{\link[rms]{Newlabels}}
##'
##' @details To give new labels for all variables, you can specify labels of the
##' form \code{labels = c("Age in Years","Cholesterol")}, where the list of new labels
##' is assumed to be the length of all main effect-type variables in the fit
##' and in their original order in the model formula. You may specify a named
##' vector to give new labels in any order for a subset of the variables,
##' e.g., \code{labels = c(age = "Age in Years", chol = "Cholesterol")}.
##'
##' @return A new \code{\link[QHScrnomo]{crr.fit}} object with adjusted labels
##'
##' @author Changhong Yu. Department of
##' Quantitative Health Sciences, Cleveland Clinic
##' @seealso \code{\link[QHScrnomo]{Newlevels.cmprsk}}
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
##' prostate.g <- Newlabels(
##' prostate.crr,
##' c(
##'   TX = 'Treatment options',
##'   BX_GLSN_CAT = 'Biopsy Gleason Score Sum',
##'   CLIN_STG = 'Clinical stage'
##' )
##' )
##'
Newlabels.cmprsk <-
  function(
    fit, # An object fit from crr.fit
    labels, # Character vector of labels to rename variables
    ...
  ) {

    # Check for fit
    if(missing(fit))
      stop("Please supply a model fit from crr.fit")

    # Check for labels
    if(missing(labels))
      stop("Please supply a character vector of labels for model variables.")

    # Check if the object was fit from crr.fit
    if(!inherits(fit, "cmprsk"))
      stop("The object is not a 'cmprsk' object fit from crr.fit")

    # Extract the original rms::cph model
    cph.f <- fit$cph.f

    # Set the new labels
    fit$cph.f <- Newlabels(cph.f, labels = labels, ... = ...)

    # Return the new fit
    fit

  }
