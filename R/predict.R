##' Calculate the Failure Time Probability from a Competing Risks Regression Model
##'
##' Computes the predicted probability of the event of interest at a specified time point for a competing risks regression model fit by \code{\link[QHScrnomo]{crr.fit}}. This function is adapted from \code{\link[cmprsk]{predict.crr}}.
##'
##' @usage \S3method{predict}{cmprsk}(object, newdata = NULL, time, lps, \dots)
##'
##' @param object A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param newdata A \code{data.frame} for prediction containing values of covariates in the model. If missing, the model development dataset (\code{object$cphdat}) is used.
##' @param time A single time point to calculate the failure probability
##' @param lps Should the linear predictor be returned instead of the failure probability? Defaults to \code{FALSE}.
##' @param ... Additional arguments such as \code{cov2} as in \code{\link[cmprsk]{crr}}
##'
##' @return A vector of failure probabilities for the specified time point (or linear predictors) with length equal to the number of rows in \code{newdata}
##' @author Michael W. Kattan, Ph.D. and Changhong Yu.\cr Department of
##' Quantitative Health Sciences, Cleveland Clinic
##' @references \code{ Fine JP and Gray RJ (1999)} A proportional hazards model
##' for the subdistribution of a competing risk.  \code{JASA} 94:496-509.
##' @seealso \code{\link[QHScrnomo]{crr.fit}}, \code{\link[cmprsk]{predict.crr}}
##' @keywords survival datagen
##'
##' @import rms
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
##' prostate.dat$pred.60 <- predict(prostate.crr, time = 60)
##'
predict.cmprsk <-
  function(
    object, # An object fit from crr.fit
    newdata = NULL, # Data frame to get predictions on
    time, # The time for evaluating probability
    lps = FALSE, # Should the linear predictor be returned instead of a probability?
    ...
  ) {

    # Check inputs when a probability is requested
    if(!lps) {

      # Check for a time point
      if(missing(time))
        stop("Please specify the time point for computing the failure probability.")

      # Check for multiple times
      if(length(time) > 1) {

        # Set to the first one
        time <- time[1]

        # Give a warning
        warning(paste0("Multiple time points supplied, but only one can be used. Defaulting to 'time=", time, "'"))

      }

      # Check for too large of a time point
      if(time > max(object$uftime))
        stop(paste0("'time=", time, "' is larger than the maximum observed failure time ", max(object$uftime), ". Please select a smaller time value."))

      # Check for too small of a time point
      if(time < min(object$uftime))
        stop(paste0("'time=", time, "' is smaller than the minimum observed failure time ", min(object$uftime), ". Please select a larger time value."))
    }

    # Get the design matrix to evaluate the model on
    if(is.null(newdata)) {

      # Set to the observed modeling data set (removing the event time/status columns)
      cov1 <- as.matrix(object$cphdat[,-c(1,2)])

    } else {

      # Otherwise transform the input data set to conform for model evaluation
      cov1 <- rms::predictrms(object$cph.f, newdata = newdata, type = "x")

    }

    # Count the model parameters
    np <- length(object$coef)

    ## Check for time-dependent covariates

    # No time dependent covariates
    if(length(object$tfs) <= 1.) {

      # If there is only a single observation to predict for
      if(length(object$coef) == length(cov1)) {

        # Compute the linear predictor
        lp <- sum(cov1 * object$coef)

        # Compute the cumulative subdistribution hazard
        lhat <- cumsum(exp(lp) * object$bfitj)

      # Otherwise we need to do matrix operations
      } else {

        # Set the matrix
        cov1 <- as.matrix(cov1)

        # Make the placeholders for linear predictor and subdistribution hazard (failure times X observations)
        lp <- matrix(0., nrow = length(object$uftime), ncol = nrow(cov1))
        lhat <- lp

        # Iterate the observations
        for(j in seq_len(nrow(cov1))) {

          # Compute the linear predictor
          lp[, j] <- sum(cov1[j, ] * object$coef)

          # Compute the cumulative subdistribution hazard
          lhat[, j] <- cumsum(exp(lp[, j]) * object$bfitj)

        }

        # The linear predictor is the same at all time points, so just take the first row
        lp <- lp[1., ]

      }

    # Time dependent covariates (assumes a 'cov2' argument has been submitted as in predict.crr)
    } else {

      # Check for only time functions (note: object$tfs = tf(object$uftime) where 'tf' is the time function as in cmprsk::crr)
      if(length(object$coef) == ncol(as.matrix(object$tfs))) {

        # Check for a single observation (assumes 'cov2' was submitted)
        if(length(object$coef) == length(cov2)) {

          # Compute the cumulative subdistribution hazard for that observation
          lhat <- cumsum(exp(object$tfs %*% c(cov2 * object$coef)) * object$bfitj)

        # Otherwise get estimates for collection of observations
        } else {

          # Convert to a matrix (generally should be already)
          cov2 <- as.matrix(cov2)

          # Make a placeholder for values (failure times by observations)
          lhat <- matrix(0., nrow = length(object$uftime), ncol = nrow(cov1))

          # Iterate to get the cumulative subdistribution hazard for each observation
          for(j in seq_len(nrow(cov2)))
            lhat[, j] <- cumsum(exp(object$tfs %*% c(cov2[j, ] * object$coef)) * object$bfitj)

        }

      # Otherwise, if there is a mix of time functions and fixed covariates
      } else {

        # Check for a single observation
        if(length(object$coef) == length(cov1) + length(cov2)) {

          # Compute the cumulative subdistribution hazard for that observation
          lhat <- cumsum(exp(sum(cov1 * object$coef[seq_len(length(cov1))]) + object$tfs %*% c(cov2, object$coef[seq(np - length(cov2) + 1., np)])) * object$bfitj)

        # Other compute for all observations
        } else {

          # Set matrices
          cov1 <- as.matrix(cov1)
          cov2 <- as.matrix(cov2)

          # Create placeholder for values (failure times by observations)
          lhat <- matrix(0., nrow = length(object$uftime), ncol = nrow(cov1))

          # Iterate to get the cumulative subdistribution hazard for each observation
          for(j in seq_len(nrow(cov1)))
            lhat[, j] <- cumsum(exp(sum(cov1[j, ] * object$coef[seq_len(ncol(cov1))]) + object$tfs %*% c(cov2[j, ] * object$coef[seq(np - ncol(cov2) + 1., np)])) * object$bfitj)

        }

      }

    }

    # Finally, if the linear predictor was requested, return that
    if(lps) {

      # Return the vector
      lp

    # Otherwise, return the failure probability
    } else {

      # Compute the failure probability (at a given time point (row) the CIF for each observation (columns))
      lhat <- cbind(object$uftime, 1. - exp(-lhat))

      # Find a cutoff for the requested time point (closest before the timepoint)
      lhat <- lhat[lhat[, 1.] <= time + 1e-10, ]

      # Take the last row (time point) in the reduced set; remove the time column
      lhat <- lhat[dim(lhat)[1], -1.]

      # Return the vector
      lhat

    }

  }
