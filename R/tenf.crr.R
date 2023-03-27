##' Ten fold cross validation for crr endpoint
##'
##'
##' Do cross validation on a competing risk regression model.
##' @title Ten fold cross validation for competing risks regression
##' @param fit a competing risks regression model fittd by function
##'   \code{\link[QHScrnomo]{crr.fit}}.
##' @param time the expected time point.
##' @param lps logical flag. If true, values of predicted X beta will be output
##'   instead of cumulative incidence
##' @param fold number of fold. the default is 10 fold cross validation.
##' @return A vector of predicted values of cumulative incidence or X beta for
##'   each observation.
##' @note Before the function is called, packages 'Hmisc', 'rms' and 'cmprsk'
##' should be loaded as the function will call some funcitons in these packages.
##' @importFrom cmprsk crr
##' @author Changhong Yu, Michael Kattan, Ph.D \cr Department of Quantitative
##'   Health Sciences\cr Cleveland Clinic\cr
##' @export
##' @seealso \code{\link[QHScrnomo]{crr.fit}}, \code{\link[cmprsk]{crr}}
##'
##' @keywords models survival
##'
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
##' prostate.dat$preds.tenf<-
##'     tenf.crr(prostate.crr,time = 120)
##' }
##'
tenf.crr <- function(fit, time = NA, lps = FALSE, fold = 10) {
    if (is.na(time)) {
        stop("Specify an expected survival time!!")
    }
    if (time > max(fit$uftime)) {
        stop(
            "the expected survival time is greater than the",
            "largest event time. Please pick a smaller time.")
    }
    if (time < min(fit$uftime)) {
        stop(
            "the expected survival time is less than the smallest event time.",
            "Please pick a greater time.")
    }
    # oldwarn <- options("warn")[[1]]
    # options(warn = -1)
    # on.exit(options(warn = oldwarn))
    assign("fit", fit)
    thedata <- fit$cphdat
    nobs <- nrow(thedata)
    nc <- ncol(thedata)
    pred <- rep(NA, length = nobs)
    rand <- sample(rep(seq(1,fold), length.out = nobs), nobs, replace = FALSE)
    cencode <- fit[["cencode"]]
    failcode <- fit[["failcode"]]
    for (i in sort(unique(rand))) {
        cat(i, " ")
        train.dat <- thedata[rand != i, ]
        test.dat <- thedata[rand == i, ]
        newfit <- crr(train.dat[, 1], train.dat[, 2],
            as.matrix(train.dat[, seq(3,nc)]),
            cencode = cencode, failcode = failcode
        )
        pred[rand == i] <- pred3.crr(newfit,
            as.matrix(test.dat[, seq(3,nc)]),
            time = time, lps = lps
        )
    }
    # add NA to rows that were excluded from analysis because of not
    # meeting the inclusion criterion in subset argument
    pred2 <- rep(NA, length = length(fit$subst))
    pred2[fit$subst] <- pred
    pred2
}
