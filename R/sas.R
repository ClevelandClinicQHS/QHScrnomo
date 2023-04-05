##' Generate a Prediction Equation for a Competing Risks Regression Model
##'
##' Uses \code{\link[rms]{Function}} to generate and print the linear predictor to the console or file which can be hard-coded into a function for evaluation on the scale of the original data. A time point can optionally be specified to retrieve the base sub-cumulative rate at that time point (i.e, the failure probability when all covariate values are 0).
##'
##' @param f A model fit by \code{\link[QHScrnomo]{crr.fit}}
##' @param time A single time point to calculate the failure probability
##' @param baseonly Should we only display the failure probability at the specified \code{time}? Defaults to \code{FALSE}.
##' @param file An optional connection or character string naming the file to print to.
##' @param append Only used if the \code{file} argument is specified. If \code{TRUE}, the output will be appended to \code{file}, otherwise it will overwritten.
##'
##' @return A printed equation using \code{\link{cat}} (invisible \code{NULL})
##'
##' @author Changhong Yu. Department of
##' Quantitative Health Sciences, Cleveland Clinic
##' @seealso \code{\link[QHScrnomo]{crr.fit}} \code{\link[rms]{sascode}} \code{\link[rms]{Function}}
##'
##' @import rms
##' @importFrom Hmisc Function
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
##' sas.cmprsk(prostate.crr, time = 60)
##'
sas.cmprsk <-
  function(
      f, # An object fit from crr.fit
      time = NA,
      baseonly = FALSE,
      file = "",
      append = FALSE
    ) {

    # Get the call
    call <- match.call()
    call <- as.list(call)

    # Check for fit
    if(missing(f))
      stop("Please supply a model fit from crr.fit")

    # Check if the object was fit from crr.fit
    if(!inherits(f, "cmprsk"))
      stop("'f' is not a 'cmprsk' object fit from crr.fit")

    # Check for a time value
    if(!is.na(time)) {

      # Check if the input time is too small
      if(time <= 0)
        stop(paste0("'time=", time, "' is not allowed. The minimum observed failure time ", min(f$uftime), ". Please select a larger time value."))

      # Check if the input time is too large
      if(time > max(f$uftime))
        stop(paste0("'time=", time, "' is larger than the maximum observed failure time ", max(f$uftime), ". Please select a smaller time value."))

      # Compute the failure probability when all covariates are 0
      lhat <- cbind(f$uftime, 1 - exp(-cumsum(f$bfitj)))

      # Keep times up to the desired time point
      lhat <- lhat[lhat[, 1.] <= time + 1e-10, ]

      # Extract the failure probability at the nearest time
      lhat <- lhat[dim(lhat)[1.], 2.]

      # If there's a file, print to that; otherwise, print to console
      if(file != "") {
        cat("Base is ", lhat, "\n", file = file, append = append)
      } else {
        cat("Base is ", lhat, "\n")
      }

    }

    # Check if the full equation should be printed
    if(!baseonly) {

      # Extract the original rms::cph model
      cph.f <- f$cph.f

      # Set the coefficients and covariance matrix that of the crr.fit model
      cph.f$coefficients <- f$coef
      cph.f$var <- f$var

      # Get the equation (from the Function.rms generic)
      funout <- as.list(Function(cph.f))
      funout <- funout[[length(funout)]]

      # Get the length
      fn <- length(funout)

      # Check for multiple equations (unclear exactly when this condition would hold, but assuming)
      if(fn > 2) {

        # Iterate elements
        for(j in 2:(fn - 1)) {

          # Separate lines into vector
          trans <- deparse(funout[[j]])

          # Replace characters
          trans <-
            Hmisc::sedit(
              text = trans,
              from = c("pmax", "pmin", "<-", "==", "^"),
              to = c("max", "min", "=", "=", "**"),
              wild.literal = TRUE
            )
          cat(trans, "\n")
        }

      }

      # Extract the equation; separate into vector
      out <- (funout)[[fn]]
      out2 <- deparse(out)

      # Identify and remove the intercept term (center of cph model; unapplicable to crr)
      pos <- gregexpr("[+-]", out2[1])[[1]]
      pos2 <- pos[pos > 1][1]
      out2[1] <-
        ifelse(
          substring(out2[1], pos2, pos2) == "+",
          substring(out2[1], pos2 + 1),
          substring(out2[1], pos2)
        )
      out3 <- out2

      # Replace characters in equation
      out4 <-
        Hmisc::sedit(
          text = out3,
          from = c("pmax", "pmin", "<-", "==", "^"),
          to = c("max", "min", "=", "=", "**"),
          wild.literal = TRUE
        )

      # If there's a file, print to that; otherwise, print to console
      if(file != "") {
        cat(out4, sep = "\n", file = file, append = TRUE)
      } else {
        cat(out4, sep = "\n")
      }

    }

    # Give a warning for no output being rendered
    if(baseonly & is.na(time))
      warning("No output produced. Please specify a 'time' value or set 'baseonly = FALSE")

  }
