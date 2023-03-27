##' Calculate concordance index
##'
##' to calculate the discrimination metric, concordance index for binary,
##' time-to event and competing risks outcomes
##'
##' @title Concordance index calculation
##' @param prob predicted risk of failure event, either probability or
##' risk score
##' @param fstatus failure(event) variable
##' @param ftime follow-up time variable for survival or competing risks
##' predictions
##' @param type type of regression models corresponding to different type
##' of outcomes.
##' 'logis' is the default value for binary outcome, 'surv' for ordinary
##' survival outcome
##' and 'crr' for competing risks outcome.
##' @param failcode coding for failure(event). 1 is the default value.
##' @param cencode coding for censoring. 0 is the defaul
##' @param tol error tolerance. the default value is 1e-20.
##' @return a vector of returned values.
##' \item{N}{the total number of observations in the input data}
##' \item{n}{the nonmissing number of observations that was used f
##' or calculation}
##' \item{usable}{the total number of usable pairs.}
##' \item{oncordant}{the number of concordant pairs}
##' \item{cindex}{the concordance index that equal to the number of
##' concordant pairs
##' divided by the total number of usable pairs.}
##' @author Changhong Yu, Michael Kattan, Brian Wells, Amy Nowacki.
##' @export
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
##' ## calculate the CRR version of concordance index
##' with(prostate.dat, cindex(preds.tenf.cv.prostate.crr.120 ,
##'                           ftime = TIME_EVENT,
##'                           fstatus =EVENT_DOD, type = "crr"))["cindex"]
##' }
##'
##' @keywords semiparametric regression
##' @useDynLib QHScrnomo , cindexCrr, .registration = TRUE
##' @useDynLib QHScrnomo , cindexLog, .registration = TRUE
##' @useDynLib QHScrnomo , cindexSurv, .registration = TRUE
##'  
##'

cindex <-
    function(
        prob, fstatus, ftime, type = "crr", failcode = 1, 
        cencode = 0, tol = 1e-20) {
        type <-
            match.arg(
                type, c("logistic", "survival", "crr"),
                several.ok = FALSE
            )
        if (
            all(regexpr(toupper(type), toupper(c(
                "logistic", "survival", "crr"
            ))) == -1)) {
            stop("type should be one of 'logistic','survival' or 'crr' !!!")
        }
        
        if (!is.na(pmatch("LOG", toupper(type)))) {
            if ((N <- length(prob)) != length(fstatus)) {
                stop(
                    "event variable has different length",
                    "from the predicted risk variable!"
                )
            }
            isna <- is.na(prob) + is.na(fstatus)
            n <- sum(isna == 0)
            prob <- prob[isna == 0]
            fstatus <- fstatus[isna == 0]
            fstatus <- ifelse(fstatus %in% failcode, 1, 0)
            
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
            if (!is.na(pmatch("SURV", toupper(type)))) {
                if (((N <- length(prob)) != length(fstatus)) |
                    (length(fstatus) != length(ftime))) {
                    stop(
                        "event variable has different length from",
                        "the predicted risk variable!"
                    )
                }
                isna <- is.na(prob) + is.na(fstatus) + is.na(ftime)
                n <- sum(isna == 0)
                prob <- prob[isna == 0]
                fstatus <- fstatus[isna == 0]
                ftime <- ftime[isna == 0]
                fstatus <- ifelse(fstatus %in% failcode, 1, 0)
                # sort the follow-up time in ascending order
                ftorder <- order(ftime)
                prob <- prob[ftorder]
                fstatus <- fstatus[ftorder]
                ftime <- ftime[ftorder]
                out <-
                    .C(
                        "cindexSurv",
                        prob = as.double(prob),
                        fstatus = as.integer(fstatus),
                        ftime = as.double(ftime),
                        n = as.integer(n),
                        npair = integer(2),
                        cindex = double(1),
                        PACKAGE = "QHScrnomo"
                    )[4:6]
            } else {
                if (((N <- length(prob)) != length(fstatus)) |
                    (length(fstatus) != length(ftime))) {
                    stop(
                        "event variable has different length from",
                        "the predicted risk variable!"
                    )
                }
                isna <- is.na(prob) + is.na(fstatus) + is.na(ftime)
                n <- sum(isna == 0)
                prob <- prob[isna == 0]
                fstatus <- fstatus[isna == 0]
                ftime <- ftime[isna == 0]
                fstatus <- ifelse(
                    fstatus %in% failcode, 1,
                    ifelse(fstatus %in% cencode, 0, 2)
                )
                # sort the follow-up time in ascending order
                ftorder <- order(ftime)
                prob <- prob[ftorder]
                fstatus <- fstatus[ftorder]
                ftime <- ftime[ftorder]
                out <-
                    .C(
                        "cindexCrr",
                        prob = as.double(prob),
                        fstatus = as.integer(fstatus),
                        ftime = as.double(ftime),
                        n = as.integer(n),
                        npair = integer(2),
                        cindex = double(1),
                        PACKAGE = "QHScrnomo"
                    )
            }
        }
        # browser()
        out <-
            c(
                N = N,
                n = n,
                usable = out$npair[1],
                concordant = out$npair[2],
                cindex = out$cindex
            )
        return(out)
    }
