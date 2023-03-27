##' Cumulative Incidence Estimates vs. a Continuous Variable
##'
##'
##' Function to divide a continuous variable \code{x} (e.g. age, or predicted
##' cumulative incidence at time \code{u} created by
##' \code{\link[QHScrnomo]{predict.cmprsk}} into \code{g} quantile groups, get
##' cumulative incidence estimates at time \code{u} (a scaler), and to return a
##' matrix with columns \code{x}=mean \code{x} in quantile, \code{n}=number of
##' subjects, \code{events}=no. events, and \code{ci}= cumulattive incidence at
##' time \code{u}, \code{std.err} = standard error. Instead of supplying
##' \code{g}, the user can supply the minimum number of subjects to have in the
##' quantile group (\code{m}, default=50). If \code{cuts} is given (e.g.
##' \code{cuts=c(0,.1,.2,\dots{},.9,.1)}), it overrides \code{m} and \code{g}.
##'
##' @title make calibration curve for competing risks endpoint
##' @param x a continuous variable
##' @param ftime vector of follow-up time
##' @param fstatus vector of failure status
##' @param cencode value indicating cencering.
##' @param failcode value indicating event of interest
##' @param ci logical flag to output event free probability if setting
##'   \code{FALSE}
##' @param m desired minimum number of observations in a group
##' @param g number of quantile groups
##' @param cuts actual cuts in \code{x}, e.g. \code{c(0,1,2)} to use [0,1),
##'   [1,2].
##' @param u time for which to estimate cumulative incidence
##' @param pl TRUE to plot results
##' @param conf.int defaults to \code{.95} for 0.95 confidence bars.  Set to
##'   \code{FALSE} to suppress bars
##' @param xlab if \code{pl=TRUE}, is x-axis label.  Default is \code{label(x)}
##'   or name of calling argument
##' @param ylab if \code{pl=TRUE}, is y-axis label.  Default is constructed
##'   from \code{u} and time \code{units}
##' @param xlim range of x axis
##' @param ylim range of y axis
##' @param lty line time for primary line connecting estimates
##' @param add set to \code{TRUE} if adding to an existing plot
##' @param cex.subtitle character size for subtitle. Default is \code{.7}.  Use
##'   \code{FALSE} to suppress subtitle.
##' @param ab \code{TRUE} to add a 45 degree line
##' @param \dots plotting parameters to pass to the plot and errbar functions
##' @return matrix with columns named \code{x} (mean predictor value in
##'   interval), \code{n} (sample size in interval), \code{events} (number of
##'   events in interval), \code{ci} (cumulative incidence estimate),
##'   \code{std.err} (standard error of cumulative incidence)
##' @note This function is adapted from Harrell's function.
##' @author Changhong Yu, Michael Kattan, Ph.D \cr Department of Quantitative
##'   Health Sciences\cr Cleveland Clinic\cr
##' @importFrom Hmisc label cut2 errbar
##' @importFrom cmprsk cuminc
##' @export
##' @seealso \code{\link[cmprsk]{cuminc}},\code{\link[QHScrnomo]{pred.ci}}
##' @keywords survival nonparametric
##'
##' @examples
##' 
##' \donttest{
##' data(prostate.dat)
##' dd <- datadist(prostate.dat)
##' options(datadist = "dd")
##' prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
##'            BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
##'            RACE_AA, data = prostate.dat,
##'            x = TRUE, y= TRUE, surv=TRUE,time.inc = 144)
##' prostate.crr <- crr.fit(prostate.f,cencode = 0,failcode = 1)
##'
##' ## ten fold cross validation
##' prostate.dat$preds.tenf.cv.prostate.crr.120 <-
##'                                        tenf.crr(prostate.crr,time = 120)
##'
##' with(prostate.dat,
##'      groupci(preds.tenf.cv.prostate.crr.120 , ftime = TIME_EVENT,
##'              fstatus =EVENT_DOD, g = 5, u = 120,
##'              xlab = "Nomogram predicted 10-year cancerspecific mortality",
##'              ylab = "Observed predicted 10-year cancerspecific mortality")
##' )
##' }
##' 
##'
groupci <-
    function(
        x,
        ftime,
        fstatus,
        cencode = 0,
        failcode = 1,
        ci = TRUE,
        m = 50,
        g,
        cuts,
        u,
        pl = TRUE,
        conf.int = 0.95,
        xlab,
        ylab,
        xlim = c(0, 1),
        ylim = c(0, 1),
        lty = 1,
        add = FALSE,
        cex.subtitle = FALSE,
        ab = TRUE,
        ...) {
        if (missing(u)) {
            stop("u (time point) must be given")
        }
        if (missing(xlab)) {
            xlab <- label(x)
        }
        if (xlab == "") {
            xlab <- as.character(sys.call()[2])
        }
        s <- !(is.na(x) | is.na(ftime) | is.na(fstatus))
        x <- x[s]
        ftime <- ftime[s]
        fstatus <- fstatus[s]
        x[abs(x) < 1e-10] <- 0
        e <- fstatus
        if (length(ftime) != length(fstatus) |
            length(ftime) != length(x)) {
            stop("lengths of x and Srv must match")
        }
        unit <- attr(ftime, "units")
        if (is.null(unit) || unit == "") {
            unit <- "Day"
        }
        if (!missing(cuts)) {
            q <- cut2(x, cuts)
        } else if (!missing(g)) {
            q <- cut2(x, g = g)
        } else {
            q <- cut2(x, m = m)
        }
        q <- oldUnclass(q)
        g <- length(levels(q))
        cmi <- single(g)
        std.err <- pred <- cmi
        numobs <- events <- integer(g)
        for (i in seq_len(g)) {
            s <- q == i
            nobs <- sum(s)
            ne <- sum(e[s] == failcode)
            if (ne == 0) {
                stop(
                    "regroup x as there is no event of interest in group ",
                    i, "\n"
                )
            }
            if (nobs == 0) {
                numobs[i] <- 0
                events[i] <- 0
                pred[i] <- NA
                cmi[i] <- NA
                std.err[i] <- NA
            } else {
                pred[i] <- mean(x[s], na.rm = TRUE)
                f <- cuminc(ftime[s], fstatus[s], cencode = cencode)
                cumci <- pred.ci(f, u, failcode = failcode)
                cmi[i] <- cumci[["CI.Prob"]]
                if (!ci) {
                    pred[i] <- 1 - pred[i]
                    cmi[i] <- 1 - cmi[i]
                }
                std.err[i] <- sqrt(cumci[["CI.Var"]])
                numobs[i] <- nobs
                events[i] <- ne
                fnm <- names(f)
                statuscode <-
                    substring(fnm, regexpr(" ", fnm) + 1)
                # to accomodate faicode = 1
                tt <- f[[fnm[regexpr(failcode, statuscode) != -1]]]$time
                n <- length(tt)
                if (u > tt[n] + 1e-06) {
                    cat(
                        "group ",
                        i,
                        ":",
                        "maxftime== ",
                        tt[n],
                        "less than u(",
                        u,
                        ")\n"
                    )
                    cmi[i] <- NA
                    std.err[i] <- NA
                }
            }
        }
        z <- cbind(
            x = pred,
            n = numobs,
            events = events,
            ci = cmi,
            std.err = std.err
        )
        if (pl) {
            y <- cmi
            if (conf.int) {
                zcrit <- stats::qnorm((conf.int + 1) / 2)
                low <- cmi - zcrit * std.err
                hi <- cmi + zcrit * std.err
                # changhong edit ???
                low[low < 0] <- 0
                hi[hi > 1] <- 1
            }
            if (missing(ylab)) {
                ylab <- paste(
                    "Competing Risks",
                    format(u),
                    "-",
                    unit,
                    " Survival",
                    sep = ""
                )
            }
            if (!add) {
                plot(
                    pred,
                    y,
                    xlab = xlab,
                    ylab = ylab,
                    type = "n",
                    xlim = xlim,
                    ylim = ylim
                )
            }
            graphics::lines(pred, y, lty = lty, ...)
            if (ab) {
                graphics::abline(0, 1)
            }
            if (conf.int) {
                errbar(pred, y, hi, low, add = TRUE, lty = lty, ...)
            }
            if (!is.logical(cex.subtitle)) {
                nn <- sum(numobs, na.rm = TRUE)
                mm <- round(nn / g)
                graphics::title(
                    sub = paste(
                        "n=",
                        nn,
                        " d=",
                        sum(events, na.rm = TRUE),
                        ",\n       avg. ",
                        mm,
                        " patients per group",
                        sep = ""
                    ),
                    adj = 0,
                    cex = cex.subtitle
                )
            }
        }
        z
    }
