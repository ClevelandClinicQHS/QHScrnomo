##' Concordance Index Calculation (C-Index)
##'
##' Computes the concordance index for a predictor as a discrimination metric for binary, time-to-event, and competing risks outcomes.
##'
##' @param prob A risk score (typically a probability giving the risk of event failure)
##' @param fstatus The event status
##' @param ftime The event times. Applies when the \code{type} argument is \code{"survival"} or \code{"crr"}
##' @param type The outcome type: \code{"logistic"} for binary, \code{"survival"} for ordinary time-to-event, and \code{"crr"} for competing risks outcomes. Defaults to \code{"crr"}.
##' @param failcode The value of \code{fstatus} that indicates the event of interest. Defaults to \code{1}.
##' @param cencode The censoring event code. Defaults to \code{0}.
##' @param tol Error tolerance (not used)
##'
##' @return A named vector with following elements:
##' \item{N}{Total number of observations in the input data}
##' \item{n}{Number of observations used for calculation}
##' \item{usable}{Total number of usable pairs}
##' \item{concordant}{Number of concordant pairs}
##' \item{cindex}{The concordance index: number of concordant pairs divided by the total number of usable pairs}
##'
##' @author Changhong Yu, Michael Kattan, Brian Wells, Amy Nowacki.
##' @keywords semiparametric regression
##'
##' @useDynLib QHScrnomo , cindexCrr, .registration = TRUE
##' @useDynLib QHScrnomo , cindexLog, .registration = TRUE
##' @useDynLib QHScrnomo , cindexSurv, .registration = TRUE
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
##'
##' # Cross-validated predictions
##' prostate.dat$preds.cv.prostate.crr.120 <- tenf.crr(prostate.crr, time = 120, fold = 2)
##'
##' ## calculate the competing-risks version of concordance index
##' with(prostate.dat, cindex(preds.cv.prostate.crr.120,
##'                           ftime = TIME_EVENT,
##'                           fstatus =EVENT_DOD, type = "crr"))["cindex"]
##'
cindex <-
  function(
      prob, # A risk score
      fstatus, # Event status
      ftime, # Event time
      type = "crr", # Type of validation
      failcode = 1, # Event status of interest
      cencode = 0, # Censoring code
      tol = 1e-20
    ) {

    # Extract the validation type
    type <- match.arg(type, c("logistic", "survival", "crr"))

    # Check for lengths; throw error if not equal
    unequal_length <- length(prob) != length(fstatus)
    if(unequal_length)
      stop("'prob' and 'fstatus' are of different lengths.")

    # Check that the supplied code is found in the input
    if(!(failcode %in% unique(fstatus)))
      stop(paste0("The supplied 'failcode=", failcode, "' was not found in 'fstatus'"))

    # Find NA cases
    isna <- is.na(prob) + is.na(fstatus)

    # Check validation type
    if(type == "logistic") {

      # Set parameters
      n <- sum(isna == 0)
      prob <- prob[isna == 0]
      fstatus <- fstatus[isna == 0]
      fstatus <- ifelse(fstatus %in% failcode, 1, 0)

      # Run C function
      out <-
        .C(
          "cindexLog",
          prob = as.double(prob),
          fstatus = as.integer(fstatus),
          n = as.integer(n),
          npair = integer(2),
          cindex = double(1),
          PACKAGE = "QHScrnomo"
        )

    } else {

      # Check for additional equal lengths
      if(length(fstatus) != length(ftime))
        stop("'fstatus' and 'ftime' are of different lengths.")

      # Add additional check for NA cases
      isna <- isna + is.na(ftime)

      # Set parameters
      n <- sum(isna == 0)
      prob <- prob[isna == 0]
      fstatus <- fstatus[isna == 0]

      # Reorder the variables
      ftorder <- order(ftime)
      prob <- prob[ftorder]
      fstatus <- fstatus[ftorder]
      ftime <- ftime[ftorder]

      # Check for type of survival
      if(type == "survival") {

        # Set the status
        fstatus <- ifelse(fstatus %in% failcode, 1, 0)

        # Set the C function name
        c_func <- "cindexSurv"

      } else {

        # Set the status
        fstatus <- ifelse(fstatus %in% failcode, 1, ifelse(fstatus %in% cencode, 0, 2))

        # Set the C function name
        c_func <- "cindexCrr"

      }

      # Run the C function
      out <-
        .C(
          c_func,
          prob = as.double(prob),
          fstatus = as.integer(fstatus),
          ftime = as.double(ftime),
          n = as.integer(n),
          npair = integer(2),
          cindex = double(1),
          PACKAGE = "QHScrnomo"
        )

      # Extract components if needed
      if(type == "survival")
        out <- out[4:6]

    }

    # Gather and return results
    c(
      N = length(prob),
      n = n,
      usable = out$npair[1],
      concordant = out$npair[2],
      cindex = out$cindex
    )

  }
