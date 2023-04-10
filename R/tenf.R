##' Obtain K-Fold Cross-Validated Predictions
##'
##' Computes "out-of-sample" predictions by K-fold cross-validation for each observation in the modeling data set from a \code{\link[QHScrnomo]{crr.fit}} object.
##'
##' @param fit A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param time A single time point to calculate the failure probability
##' @param lps Should the linear predictor be returned instead of the failure probability? Defaults to \code{FALSE}.
##' @param fold The number of folds. Defaults to \code{10}.
##' @param trace Should the progress of cross-validation be printed to the console? Defaults to \code{TRUE}.
##'
##' @return A vector of failure probabilities at the specified time point (or linear predictors if \code{lps=TRUE}) with length equal to the number of rows in the original data set.
##'
##' @author Changhong Yu, Michael Kattan, Ph.D \cr Department of Quantitative
##'   Health Sciences\cr Cleveland Clinic\cr
##' @seealso \code{\link[QHScrnomo]{crr.fit}}, \code{\link[cmprsk]{crr}}
##' @keywords models survival
##'
##' @export
##'
##' @examples
##' dd <- datadist(prostate.dat)
##' options(datadist = "dd")
##' prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
##'            BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
##'            RACE_AA, data = prostate.dat,
##'            x = TRUE, y = TRUE, surv = TRUE,time.inc = 144)
##' prostate.crr <- crr.fit(prostate.f, cencode = 0, failcode = 1)
##' tenf.crr(prostate.crr, time = 120, fold = 3)
##'
tenf.crr <-
  function(
      fit, # An object fit from crr.fit
      time = NULL, # Time point for validating at
      lps = FALSE, # Should the cross-validated values be returned on the linear predictor scale?
      fold = 10, # Number of folds
      trace = TRUE # Should the progress be traced?
    ) {

    # Check for fit
    if(missing(fit))
      stop("Please supply a model fit from crr.fit")

    # Check if the object was fit from crr.fit
    if(!inherits(fit, "cmprsk"))
      stop("The fit is not a 'cmprsk' object fit from crr.fit")

    # Check inputs when a probability is requested
    if(!lps) {

      # Check for a time point
      if(is.null(time))
        stop("Please specify the time point for computing the failure probability.")

      # Check for multiple times
      if(length(time) > 1) {

        # Set to the first one
        time <- time[1]

        # Give a warning
        warning(paste0("Multiple time points supplied, but only one can be used. Defaulting to 'time=", time, "'"))

      }

      # Check for too large of a time point
      if(time > max(fit$uftime))
        stop(paste0("'time=", time, "' is larger than the maximum observed failure time ", max(fit$uftime), ". Please select a smaller time value."))

      # Check for too small of a time point
      if(time < min(fit$uftime))
        stop(paste0("'time=", time, "' is smaller than the minimum observed failure time ", min(fit$uftime), ". Please select a larger time value."))

    } else {

      # Check for time-dependent terms
      if(length(fit$tfs) > 1.)
        stop("There are time-dependent terms present. Cannot return linear-predictor. Please rerun with 'lps=FALSE' and specify a 'time' value.")

      # Check for a time value
      if(!is.null(time))
        warning("Since 'lps=TRUE', the time argument was disregarded.")

    }

    # Assign the fit
    assign("fit", fit)

    # Extract the model-fit data
    thedata <- fit$cphdat

    # Extract dimensions
    nobs <- nrow(thedata)
    nc <- ncol(thedata)

    # Set up some parameters
    pred <- rep(NA, length = nobs)
    cencode <- fit[["cencode"]]
    failcode <- fit[["failcode"]]

    # Create the folds
    folds <- seq_len(fold)
    rand <- sample(rep(folds, length.out = nobs), nobs, replace = FALSE)

    # Iterate the folds
    for(i in folds) {

      # Trace the folds (if requested)
      if(trace)
        cat(i, " ")

      # Split the data
      train.dat <- thedata[rand != i, ]
      test.dat <- thedata[rand == i, ]

      # Fit the model
      newfit <-
        cmprsk::crr(
          ftime = train.dat[, 1],
          fstatus = train.dat[, 2],
          cov1 = as.matrix(train.dat[,seq(3, nc)]),
          cencode = cencode,
          failcode = failcode
        )

      # Obtain the predictions
      pred[rand == i] <-
        pred3.crr(
          z = newfit,
          cov1 = as.matrix(test.dat[,seq(3, nc)]),
          time = time,
          lps = lps
        )

    }

    # Fill in NA's for observations missing from the original fit
    pred2 <- rep(NA, length = length(fit$subst))
    pred2[fit$subst] <- pred

    # Return the vector
    pred2

  }
