##' Extract Cumulative Incidence Estimates at a Specified Time Point
##'
##' Extracts the cumulative incidence estimates from a \code{\link[cmprsk]{cuminc}} object for the cause of interest at a specified time point into a \code{\link{data.frame}}.
##'
##' @param cum A \code{\link[cmprsk]{cuminc}} object
##' @param tm1 A single time point to return the cumulative incidence at
##' @param failcode The value of the status column that indicates the event of interest
##'
##' @return A \code{\link{data.frame}} with 3 columns:
##' \item{Group}{The group name. If the \code{group} argument was used in the \code{\link[cmprsk]{cuminc}} fit, there will be one row per group. Otherwise this is non-informative.}
##' \item{CI.Prob}{The cumulative incidence probability at the desired time point}
##' \item{CI.Var}{The estimated variance of the cumulative incidence estimate}
##'
##' @author Michael W. Kattan, Ph.D. and Changhong Yu.\cr Department of Quantitative Health Sciences, Cleveland Clinic
##' @seealso \code{\link[cmprsk]{cuminc}}
##'
##' @export
##'
##' @examples
##'  cum <- cmprsk::cuminc(prostate.dat$TIME_EVENT, prostate.dat$EVENT_DOD, cencode = 0)
##'  pred.ci(cum, 60, failcode = 1)
##'
pred.ci <-
  function(
      cum, # A cmprsk::cuminc object
      tm1, # Time point to evaluate cumulative incidence
      failcode = 1 # The desired cause to return the cumulative incidence for
    ) {

    # Check for object
    if(missing(cum))
      stop("Please supply a cmprsk::cuminc object")

    # Check if the object was fit from cmprsk::cuminc
    if(!inherits(cum, "cuminc"))
      stop("The object is not a 'cuminc' object fit from cmprsk::cuminc")

    # Check for a time point
    if(missing(tm1))
      stop("Please specify the time point to extract the cumulative incidence rate at.")

    # Check for multiple times
    if(length(tm1) > 1) {

      # Set to the first one
      tm1 <- tm1[1]

      # Give a warning
      warning(paste0("Multiple time points supplied, but only one can be used. Defaulting to 'tm1=", tm1, "'"))

    }

    # Extract object names
    nms <- names(cum)

    # Parse the event status codes
    statuscode <- substring(nms, regexpr(" ", nms) + 1)

    # Find the names of the list elements associated with the cause of interest
    subgrp <- nms[regexpr(failcode, statuscode) != -1]

    # Check for too large of a time point (larger than the )
    max_time <- max(vapply(subgrp, function(x) max(cum[[x]]$time), 0))
    if(tm1 > max_time)
      stop(paste0("'tm1=", tm1, "' is larger than the maximum observed failure time ", max_time, ". Please select a smaller time value."))

    # Check for too small of a time point
    if(tm1 < 0)
      stop(paste0("The time point cannot be negative. Please select a time value greater than 0."))

    # Extract the estimated cumulative incidence and variance at the closest observed failure time for each group
    temp_cuminc <-

      # Iterate the list
      sapply(

        # For each group element of the desired subset
        subgrp,

        # Apply this function
        function(i) {

          # Extract this element
          this_cum <- cum[[i]]

          # Extract the cumulative incidence curves into a matrix for the current element
          lhat <- cbind(this_cum$time, this_cum$est, this_cum$var)

          # Keep rows up to the desired time point
          lhat <- lhat[lhat[,1] < tm1, ]

          # Now only keep the row closest to the desired time point
          lhat <- lhat[nrow(lhat), ]

          lhat

        },

        # Convert to a matrix
        simplify = TRUE

      )

    # Remove the cause from the group name(s)
    colnames(temp_cuminc) <- substring(colnames(temp_cuminc), 1, regexpr(paste0("[ ]", failcode, "$"), colnames(temp_cuminc)) - 1)

    # Transpose the matrix so each group is along the rows
    temp_cuminc <- t(temp_cuminc)

    # Make a data frame
    data.frame(Group = rownames(temp_cuminc), CI.Prob = temp_cuminc[,2], CI.Var = temp_cuminc[,3], row.names = NULL)

  }
