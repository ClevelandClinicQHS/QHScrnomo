#' Fit A Competing Risks Regression Model
#'
#' Fits a competing risks regression model using the \code{\link[cmprsk]{crr}} function from an existing \code{\link[rms]{cph}}
#' object which can then be used to construct a nomogram.
#'
#' @param fit A Cox proportional hazards regression model constructed from \code{\link[rms]{cph}} (by Frank Harrell)
#' @param cencode The value of the status column that indicates a censored observation
#' @param failcode The value of the status column that indicates an event of interest
#'
#' @return Returns a list of class \code{cmprsk}, with components:
##' \item{coef }{the estimated regression coefficients}
##' \item{loglik }{log pseudo-liklihood evaluated at coef}
##' \item{lscore }{derivitives of the log pseudo-likelihood evaluated at coef}
##' \item{inf}{-second derivatives of the log pseudo-likelihood}
##' \item{var}{estimated variance covariance matrix of coef}
##' \item{res}{matrix of residuals giving the contribution to each score
##' (columns) at each unique failure time (rows)}
##'  \item{uftime}{vector of unique failure times}
##' \item{bfitj}{jumps in the Breslow-type estimate of the underlying
##' sub-distribution cumulative hazard (used by predict.crr())}
##' \item{tfs}{the tfs matrix (output of  tf(), if used)}
##' \item{converged}{TRUE if the iterative algorithm converged.}
##' \item{cencode }{the value of the status indicator that indicates
##' a censored observation}
##' \item{failcode}{the value of the status indicator that indicates an
##' event of interest}
##' \item{cph.f}{regular survival model fitted by cph which is saved for
##'   function \code{\link[QHScrnomo]{nomogram.crr}} to adjust lp for
##'   competing risks}
##' \item{cphdat}{data used for cph model, where all
##'   predictors are represented in numeric format, which is used by function
##'   \code{\link[QHScrnomo]{tenf.crr}} to do ten fold cross-validation}
##' @note This function requires that the \code{\link{rms}} and \code{\link{cmprsk}} libraries are
##'   attached.
##' @author Michael W. Kattan, Ph.D. and Changhong Yu. Department of
##'   Quantitative Health Sciences, Cleveland Clinic
##' @importFrom cmprsk crr
##' @export
##' @seealso \code{\link[rms]{cph}} \code{\link[cmprsk]{crr}}
##'   \code{\link[QHScrnomo]{nomogram.crr}}
##' @references Michael W. Kattan, Glenn Heller and Murray F. Brennan (2003). A
##'   competing-risks nomogram for sarcoma-specific death following local
##'   recurrence. Statistics in Medicine. \code{Stat Med}. 2003;22:3515-3525.
#'
##' @examples
##' 1
##'
##' @keywords survival multivariate
##'
crr.fit <-
  function(
    fit, # An rms::cph object
    cencode = 0, # Censoring indicator
    failcode = 1 # Event code of interest
  ) {

    # Record the call
    call <- match.call()

    # Check for fit
    if(missing(fit))
      stop("Please supply a model fit from rms::cph")

    # Check for the right object
    if(!methods::is(fit, "cph"))
      stop("The model fit must be from rms::cph")

    # Check if the design was part of the the original model call (update if not)
    if(!("x" %in% names(fit)))
      fit <- stats::update(fit, x = TRUE)

    # Retrieve the data set from the original model that is loaded into memory
    thedata <- get(paste(as.list(fit$call)$data))

    # Check if a subset of data was used in the original fit
    if(is.null(as.list(fit$call)$subset)) {

      # No subset --> make a vector inclusive of all rows
      subst <- rep(TRUE, nrow(thedata))

    } else {

      # Execute the subsetting expression on the data set (produces vector of T/F) (e.g., Produces: with(data, AGE>50); Equivalent: data$AGE>50)
      subst <- with(thedata, eval(as.list(fit$call)$subset))

      # Set any NA values to FALSE as well
      subst[!subst | is.na(subst)] <- FALSE

    }

    # Extract the modeling data set
    data <- thedata[subst,]

    # Retrieve the name of the event time column in the original data
    timevar <- as.character(stats::formula(fit$call)[[2]])[2]
    timevar <- ifelse(is.null(timevar), 1, timevar)

    # Retrieve the name of the event column/status call
    x <- as.character(stats::formula(fit$call)[[2]])[3]

    # Search for event binary declaration in the call
    event_declared <- regexpr("[^[:alnum:]._]", x)

    # Remove event declaration from the event call (e.g., 'status == 1'; want 'status')
    statvar <- substring(x, 1, ifelse(!event_declared, 100000, event_declared - 1))

    # Fit the competing risk model
    newfit <-
      cmprsk::crr(
        ftime = data[, timevar],
        fstatus = data[, statvar],
        cov1 = as.matrix(fit$x[subst,]),
        cencode = cencode,
        failcode = failcode
      )

    # Generate a design matrix
    design.m <- rms::predictrms(fit, newdata = data, type = "x")

    # Set some attributes
    attr(newfit$coef, "names") <- names(fit$coefficients)
    newfit$cencode <- cencode # Used for tenf.crr
    newfit$failcode <- failcode
    newfit$cph.f <- fit # Used for nomogram.crr
    newfit$cphdat <- data.frame(thedata[subst, c(timevar, statvar)], design.m)
    newfit$timevar <- timevar
    newfit$statvar <- statvar
    newfit$subst <- subst

    # Check for large jumps in Breslow-estimator
    if(mean(newfit$bfitj) > 1e+100)
      stop("Infinite 'bfitj' causes by predictors with large integers!")

    # Set the class
    oldClass(newfit) <- c("cmprsk", "crr")

    # Return the new fit
    return(newfit)

  }
