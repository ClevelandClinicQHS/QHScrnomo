##' Calculate Cumulative Incidence
##'
##'
##' Extract cumulative incidence and its variance from a object generated from
##' function \code{\link[cmprsk]{cuminc}}.
##'
##' @title Calculate Cumulative Incidence
##' @param cum a object from function \code{\link[cmprsk]{cuminc}}
##' @param tm1 expected failure time
##' @param failcode value indicating the event of interest
##' @return a data frame with 3 columns.  \item{column 3: }{Group name.}
##'   \item{column 2: }{Cumulative Incidence Probability.} \item{column 3:
##'   }{Variance}
##' @author Michael W. Kattan, Ph.D. and Changhong Yu.\cr Department of
##'   Quantitative Health Sciences, Cleveland Clinic
##' @export
##' @seealso \code{\link[cmprsk]{cuminc}}
##' @examples
##'
##'  data(prostate.dat) # get demo data set
##'  cum <- cuminc(prostate.dat$TIME_EVENT,prostate.dat$EVENT_DOD,
##'                cencode = 0)
##'  # calculate the expected cumulative incidence by 5 year for death from
##'  # prostate cancer
##'  # Here, code for cause A is 'DOA'.
##'  pred.ci(cum,60,failcode = 1)
##'
##' @keywords survival datagen
##'
pred.ci <- function(cum, tm1, failcode) {
    nms <- names(cum)
    statuscode <-
        substring(nms, regexpr(" ", nms) + 1) # to accomodate faicode = 1
    subgrp <- nms[regexpr(failcode, statuscode) != -1]
    if (tm1 > max(vapply(subgrp, function(x) {
        max(cum[[x]]$time)
    }, 0))) {
        stop("expected failure time is too large !!!")
    }
    outmatrix <- data.frame(matrix(ncol = 3, nrow = 100))
    j <- 1
    for (i in nms) {
        ln <- nchar(i)
        fc <- substring(i, regexpr("[ ]", i) + 1, ln)
        if (fc == as.character(failcode)) {
            lhat <- cbind(cum[[i]]$time, cum[[i]]$est, cum[[i]]$var)
            lhat <-
                lhat[max((seq(1, length(cum[[i]]$time)))[cum[[i]]$time <=
                    tm1 + 1e-10]), c(2, 3)]
            outmatrix[j, c(2, 3)] <- lhat
        }
        outmatrix[j, 1] <- substring(i, 1, regexpr("[ ]", i) -
            1)
        j <- j + 1
    }
    outpt <- stats::na.omit(outmatrix)
    names(outpt) <- c("Group", "CI.Prob", "CI.Var")
    outpt
}
