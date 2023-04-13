##' Assess Calibration for a Competing Risks Endpoint
##'
##' Uses \code{\link[cmprsk]{cuminc}} to estimate the cumulative incidence at a given time point within subgroups of a continuous variable (often predicted failure probabilities from a \code{\link[QHScrnomo]{crr.fit}} model).
##'
##' @param x A numeric variable to assess calibration for
##' @param ftime The event time variable. See \code{\link[cmprsk]{cuminc}}.
##' @param fstatus The event status variable. See \code{\link[cmprsk]{cuminc}}.
##' @param u A single time point to assess calibration at
##' @param cencode The censoring event code. See \code{\link[cmprsk]{cuminc}}.
##' @param failcode The value of \code{fstatus} that indicates the event of interest
##' @param ci Should the failure probability be assessed? Defaults to \code{TRUE}. If \code{FALSE}, the event-free probability is assessed.
##' @param m Minimum number of observations in each group. See \code{\link[Hmisc]{cut2}}.
##' @param g Number of quantile groups. See \code{\link[Hmisc]{cut2}}.
##' @param cuts Actual cut points to use for \code{x}. See \code{\link[Hmisc]{cut2}}.
##' @param pl Should the calibration curve be plotted? Defaults to \code{TRUE}.
##' @param conf.int Confidence limit on error bars. Defaults to 0.95. Set to \code{FALSE} to suppress.
##' @param xlab The x-axis label. Uses \code{\link[Hmisc]{label}} or name of calling argument if not specified.
##' @param ylab The y-axis label. Uses a default label is none specified.
##' @param xlim The x-axis limits. Defaults to c(0, 1).
##' @param ylim The y-axis limits. Defaults to c(0, 1).
##' @param lty Line type for connecting estimates and error bars
##' @param add Defaults to \code{FALSE}. Set to \code{TRUE} to add to an existing plot.
##' @param cex.subtitle Character size for subtitle (default 0.7). Defaults to \code{FALSE} to suppress.
##' @param ab Should a reference line be added? See \code{\link[graphics]{abline}}.
##' @param a The intercept for the reference line. See \code{\link[graphics]{abline}}.
##' @param b The slope for the reference line. See \code{\link[graphics]{abline}}.
##' @param ... Other arguments passed to \code{\link[graphics]{lines}} and \code{\link[Hmisc]{errbar}}.
##'
##' @details To divide \code{x}, the function first looks for \code{cuts}, then \code{g}, then \code{m}.
##'
##' @return A matrix with a row for each group of \code{x}:
##' \item{x}{Mean value of \code{x}}
##' \item{n}{Number of observations}
##' \item{events}{Number of events (of type \code{failcode})}
##' \item{ci}{Estimated cumulative incidence (or event-free probability if \code{ci=FALSE})}
##' \item{std.err}{Estimated standard error for the \code{ci} value}
##' If \code{pl=TRUE}, a calibration plot is also displayed.
##'
##' @author Changhong Yu, Michael Kattan, Ph.D \cr Department of Quantitative Health Sciences\cr Cleveland Clinic\cr
##' @seealso \code{\link[cmprsk]{cuminc}} \code{\link[QHScrnomo]{pred.ci}} \code{\link[Hmisc]{cut2}}
##' @keywords survival nonparametric
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
##'
##' # Cross-validated predictions
##' prostate.dat$preds.cv.prostate.crr.120 <- tenf.crr(prostate.crr, time = 120, fold = 3)
##'
##' with(prostate.dat,
##'      groupci(preds.cv.prostate.crr.120, ftime = TIME_EVENT,
##'              fstatus =EVENT_DOD, g = 5, u = 120,
##'              xlab = "Nomogram predicted 10-year cancerspecific mortality",
##'              ylab = "Observed predicted 10-year cancerspecific mortality")
##' )
##'
groupci <-
  function(
      x,
      ftime,
      fstatus,
      u,
      cencode = 0,
      failcode = 1,
      ci = TRUE,
      m = 50,
      g = NULL,
      cuts = NULL,
      pl = TRUE,
      conf.int = 0.95,
      xlab = NULL,
      ylab = NULL,
      xlim = c(0, 1),
      ylim = c(0, 1),
      lty = 1,
      add = FALSE,
      cex.subtitle = FALSE,
      ab = TRUE,
      a = 0,
      b = 1,
      ...
    ) {

    # Check for missing required inputs
    if(missing(x))
      stop("Please supply a continuous variable 'x'.")
    if(missing(ftime))
      stop("Please supply a vector for 'ftime' indicating the follow-up time.")
    if(missing(fstatus))
      stop("Please supply a vector for 'fstatus' indicating the event status.")
    if(missing(u))
      stop("Please supply a time point of interest 'u' to evaluate calibration at.")

    # Verify that the inputs are the same length
    input_lengths <- c(length(x), length(ftime), length(fstatus))
    if(length(unique(input_lengths)) > 1)
      stop(paste0("'x', 'ftime', 'fstatus' must be the same lengths. The current lengths are ", paste(input_lengths, collapse = ', '), ", respectively."))

    # Check that the supplied codes are found in the inputs
    unique_statuses <- unique(fstatus)
    if(!(cencode %in% unique_statuses))
      stop(paste0("The supplied 'cencode=", cencode, "' was not found in 'fstatus'"))
    if(!(failcode %in% unique_statuses))
      stop(paste0("The supplied 'failcode=", failcode, "' was not found in 'fstatus'"))

    # Check for multiple times
    if(length(u) > 1) {

      # Set to the first one
      u <- u[1]

      # Give a warning
      warning(paste0("Multiple time points supplied, but only one can be used. Defaulting to 'u=", u, "'"))

    }

    # Set the x-label
    if(is.null(xlab)) {

      # Extract the label
      xlab <- Hmisc::label(x)

      # Check for no value, and set to variable name
      if(xlab == "")
        xlab <- as.character(sys.call()[2])

    }

    # Indicate missingness of the vectors
    s <- !(is.na(x) | is.na(ftime) | is.na(fstatus))

    # Index the subset of complete cases (if necessary)
    if(sum(!s) > 0) {

      # Extract the subsets
      x <- x[s]
      ftime <- ftime[s]
      fstatus <- fstatus[s]

      # Give a warning to user
      warning(paste0(sum(!s), " observations removed due to missingness."))

    }

    # Set implied zeros
    x[abs(x) < 1e-10] <- 0

    # Extract/set the units
    unit <- attr(ftime, "units")
    if(is.null(unit) || unit == "")
      unit <- "Day"

    # Set the groups in order of precedence depending on user input
    if(!is.null(cuts)) {

      # User-supplied cut points
      q <- Hmisc::cut2(x, cuts = cuts)

    } else if(!is.null(g)) {

      # Number of quantile (equally-sized) groups
      q <- Hmisc::cut2(x, g = g)

    } else {

      # Minimum number of observations per group
      q <- Hmisc::cut2(x, m = m)

    }

    # Remove the factor class
    q <- unclass(q)

    # Get the group count
    g <- length(levels(q))

    # Set some parameters
    cmi <- double(g)
    std.err <- pred <- cmi
    numobs <- events <- integer(g)
    e <- fstatus

    # Iterate the groups
    for(i in seq_len(g)) {

      # Find members of the current group
      s <- q == i

      # Group observation and event volumes
      nobs <- sum(s)
      ne <- sum(e[s] == failcode)

      # Check for events
      if(ne == 0)
        stop(paste0("There are no events of type 'failcode=", failcode, "' in group ", i, ". Please regroup 'x'."))

      # Check for observations
      if(nobs == 0) {

        # Set to null values
        numobs[i] <- 0
        events[i] <- 0
        pred[i] <- NA
        cmi[i] <- NA
        std.err[i] <- NA

      # If there are observations and events, do some work
      } else {

        # Observed group average
        pred[i] <- mean(x[s], na.rm = TRUE)

        # Actual group cumulative incidence
        f <- cmprsk::cuminc(ftime[s], fstatus[s], cencode = cencode)
        cumci <- pred.ci(f, u, failcode = failcode)
        cmi[i] <- cumci[["CI.Prob"]]

        # Change to survival probability (if requested)
        if(!ci)
          cmi[i] <- 1 - cmi[i]

        # Set additional results
        std.err[i] <- sqrt(cumci[["CI.Var"]])
        numobs[i] <- nobs
        events[i] <- ne

        # Find the largest failure time
        fnm <- names(f)
        statuscode <- substring(fnm, regexpr(" ", fnm) + 1)
        tt <- f[[fnm[regexpr(failcode, statuscode) != -1]]]$time
        n <- length(tt)

        # Check time point
        if(u > tt[n] + 1e-06) {

          # Set to nulls
          cmi[i] <- NA
          std.err[i] <- NA

          # Give a warning
          warning(paste0("The maximum failure time for group ", i, " is ", tt[n], " which is less than the requested time point 'u=", u, "'. Cannot evaluate."))

        }

      }

    }

    # Place into a matrix
    z <-
      cbind(
        x = pred,
        n = numobs,
        events = events,
        ci = cmi,
        std.err = std.err
      )

    # Check if a plot should be produced
    if(pl) {

      # Set y variable
      y <- cmi

      # Get the y-axis label
      if(is.null(ylab))
        ylab <- paste0("Observed ", format(u), "-", unit, " competing risks probability")

      # Check to create a new plot
      if(!add)
        graphics::plot(x = pred, y = y, xlab = xlab, ylab = ylab, type = "n", xlim = xlim, ylim = ylim)

      # Add lines to the plot
      graphics::lines(pred, y, lty = lty, ...)

      # Check to add an ab line
      if(ab)
        graphics::abline(a, b)

      # Check to compute confidence limits
      if(conf.int) {

        # Critical value
        zcrit <- stats::qnorm((conf.int + 1) / 2)

        # Set limits
        low <- cmi - zcrit*std.err
        hi <- cmi + zcrit*std.err

        # Set implied boundaries
        low[low < 0] <- 0
        hi[hi > 1] <- 1

        # Add to plot
        Hmisc::errbar(pred, y, hi, low, add = TRUE, lty = lty, ...)

      }

      # Check for subtitle
      if(!is.logical(cex.subtitle)) {

        # Compute average observation count
        nn <- sum(numobs, na.rm = TRUE)
        mm <- round(nn / g)

        # Add subtitle to graph
        graphics::title(sub = paste0("n=", nn, " d=", sum(events, na.rm = TRUE), ",\n\tavg. ", mm, " patients per group"), adj = 0, cex = cex.subtitle)

      }
    }

    # Return the data matrix
    z

  }
