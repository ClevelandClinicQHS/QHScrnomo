##' Fit Competing Risks Regression Model
##' Fits a competing risks regression model from an existing Cox proportional
##' hazards object and allows a nomogram to be constructed from the
##' competing risks regression object.
##'
##' This function uses the \code{\link[cmprsk]{crr}} function in
##' the \code{cmprsk} package to construct a competing risk regression object.
##'
##' @param fit a Cox proportional hazards regression model constructed from
##'   \code{\link[rms]{cph}} in rms library (by Frank Harrell)
##' @param cencode the value of the status indicator that indicates a censored
##'   observation
##' @param failcode the value of the status indicator that indicates an event
##'   of interest
##' @return Returns a list of class cmprsk, with components:
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
##'
##' @note This function requires that the rms and cmprsk libraries are
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
##' ## anova test
##' anova(prostate.crr)
##' ## hazards ratio
##' summary(prostate.crr)
##'
##' ## ten fold cross validation
##' prostate.dat$preds.tenf.cv.prostate.crr.120 <-
##'                                        tenf.crr(prostate.crr,time = 120)
##'
##' ## make a CRR nomogram
##' nomogram.crr(prostate.crr,failtime = 120,lp=FALSE,
##' funlabel = "Predicted 10-year cumulative incidence")
##'
##' ## calculate the CRR version of concordance index
##' with(prostate.dat, cindex(preds.tenf.cv.prostate.crr.120 ,
##'                           ftime = TIME_EVENT,
##'                           fstatus =EVENT_DOD, type = "crr"))["cindex"]
##'
##' ## generate the calibration curve for predicted 10-year cancer
##' ## specific mortality
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
##' @keywords survival multivariate
##'
crr.fit <- function(fit,
                    cencode = 0,
                    failcode = 1) {
    call <- match.call()
    fit <- fit
    if (!match("x", names(fit))) {
        fit <- stats::update(fit, x = TRUE)
    }
    thedata <- get(paste(as.list(fit$call)$data)) # original data
    if (is.null(as.list(fit$call)$subset)) {
        subst <- rep(TRUE, nrow(thedata))
    } else {
        subst <- with(thedata, eval(as.list(fit$call)$subset))
        # exclude subjects with missing subset from model building
        subst[!subst | is.na(subst)] <- FALSE
    }
    # browser()
    data <-
        thedata[subst, ] # dataset used for model fitting
    timevar <- as.character(stats::formula(fit$call)[[2]])[2]
    timevar <- ifelse(is.null(timevar), 1, timevar)
    #assign("timevar", timevar, 1)
    x <- as.character(stats::formula(fit$call)[[2]])[3]
    statvar <-
        substring(x, 1, ifelse(
            regexpr("[^[:alnum:]._]", x) == -1,
            100000,
            regexpr("[^[:alnum:]._]", x) - 1
        ))
    statvar <- ifelse(is.null(statvar), 1, statvar)
    #assign("statvar", statvar, 1)
    # sub <- as.numeric(row.names(data.frame(fit$x)))
    #sub <- row.names(data.frame(fit$x))
    newfit <- crr(
        data[, timevar],
        data[, statvar],
        as.matrix(fit$x[subst,]),
        cencode = cencode,
        failcode = failcode
    )
    attr(newfit$coef, "names") <- names(fit$coefficients)
    newfit$cencode <- cencode # used for tenf.crr
    newfit$failcode <- failcode
    design.m <- rms::predictrms(fit, newdat = data, type = "x")
    newfit$cph.f <- fit # used for nomogram.crr
    newfit$cphdat <-
        data.frame(thedata[subst, c(timevar, statvar)], design.m)
    # newfit$call <- call
    newfit$timevar <- timevar
    newfit$statvar <- statvar
    newfit$subst <- subst
    if (mean(newfit$bfitj) > 1e+100) {
        stop(
            "infinite bfitj caused by predictors with",
            "large integers !"
        )
    }
    oldClass(newfit) <- c("cmprsk", "crr")
    return(newfit)
}
