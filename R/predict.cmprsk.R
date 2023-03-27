##' Calculate Predicted Competing Risks Probability.
##'
##'
##' Calculate predicted probabilities for a competing risks regression model,
##' which is fitted by function \code{\link[QHScrnomo]{crr.fit}}.
##'
##' @title Calculate Predicted Competing Risks Probability
##' @usage \S3method{predict}{cmprsk}(object, newdata = NULL, time, 
##' lps = FALSE, \dots)
##' @param object a saved crr model fit crated by function
##' \code{crr.fit}
##' @param newdata data frame for prediction. Each row of the data
##' frame   contains values of covariates that are required in the crr
##' model. If   missing, the original data set that was used to
##' develop the crr model   will be used for prediction.
##' @param time expected time point for evaluating the competing risks
##' probability.
##' @param lps set \code{TRUE} to return linear predictor values
##' instead of   failure probabilities.
##' @param ... other arguments
##' @return A vector with the length equal to the number of rows in the data
##'   frame, which was used to make prediction. Each element corresponds to a
##'   predicted failure probability at the expected time point.
##' @note This function is adapted from function
##'   \code{\link[cmprsk]{predict.crr}} in package \code{cmprsk}.
##' @author Michael W. Kattan, Ph.D. and Changhong Yu.\cr Department of
##'   Quantitative Health Sciences, Cleveland Clinic
##' @export
##' @seealso \code{\link[QHScrnomo]{crr.fit}},
##'   \code{\link[cmprsk]{predict.crr}}
##' @references \code{ Fine JP and Gray RJ (1999)} A proportional hazards model
##'   for the subdistribution of a competing risk.  \code{JASA} 94:496-509.
##' @keywords survival datagen
##'
##' @examples 
##' 
##' data(prostate.dat)
##' dd <- datadist(prostate.dat)
##' options(datadist = "dd")
##' prostate.f <- cph(Surv(TIME_EVENT,EVENT_DOD == 1) ~ TX  + rcs(PSA,3) +
##'            BX_GLSN_CAT +  CLIN_STG + rcs(AGE,3) +
##'            RACE_AA, data = prostate.dat,
##'            x = TRUE, y= TRUE, surv=TRUE,time.inc = 144)
##' prostate.crr <- crr.fit(prostate.f,cencode = 0,failcode = 1)
##' prostate.dat$pred.60 <- predict(prostate.crr, time=60)
##' 
##' 
predict.cmprsk <- function(object, newdata = NULL, time, lps = FALSE, ...){
    if(!lps){
        if(missing(time)) stop("specify expected time point.")
        if(time > max(object$uftime)) stop("select a  smaller time.")
        if(time < min(object$uftime)) stop("select a larger time.")
    }
    if(is.null(newdata)) cov1 <- as.matrix(object$cphdat[,-c(1,2)]) else
        cov1 <- predictrms(object$cph.f,newdata,type="x")
    np <- length(object$coef)
    if(length(object$tfs) <= 1.) {
        if(length(object$coef) == length(cov1)) {
            lhat <- cumsum(exp(sum(cov1 * object$coef)) * object$bfitj)
            lp <- sum(cov1 * object$coef)
        }
        else {
            cov1 <- as.matrix(cov1)
            lhat <- matrix(0., nrow = length(object$uftime), ncol = nrow(
                cov1))
            lp <- matrix(0., nrow = length(object$uftime), ncol = nrow(
                cov1))
            for(j in seq_len(nrow(cov1))) {
                lhat[, j] <- cumsum(exp(sum(cov1[j,  ] * object$
                                                coef)) * object$bfitj)
                lp[, j] <- sum(cov1[j,  ] * object$coef)
            }
            lp <- lp[1.,  ]
        }
    }
    else {
        if(length(object$coef) == ncol(as.matrix(object$tfs))) {
            if(length(object$coef) == length(cov2))
                lhat <- cumsum(
                    exp(object$tfs %*% c(cov2 * object$coef)) * object$bfitj)
            else {
                cov2 <- as.matrix(cov2)
                lhat <- matrix(
                    0., nrow = length(object$uftime), ncol = nrow(cov1))
                for(j in seq_len(nrow(cov2)))
                    lhat[, j] <- cumsum(exp(object$tfs %*% c(
                        cov2[j,  ] * object$coef)) * object$bfitj)
            }
        }
        else {
            if(length(object$coef) == length(cov1) + length(cov2))
                lhat <- cumsum(exp(sum(cov1 * object$coef[seq_len(length(
                    cov1))]) + object$tfs %*% c(cov2 * object$coef[
                        seq(np - length(cov2) + 1.,np)])) * object$
                        bfitj)
            else {
                cov1 <- as.matrix(cov1)
                cov2 <- as.matrix(cov2)
                lhat <- matrix(
                    0., nrow = length(object$uftime), ncol = nrow(cov1))
                for(j in seq_len(nrow(cov1)))
                    lhat[, j] <- cumsum(exp(sum(cov1[j,
                    ] * object$coef[seq_len(ncol(cov1))]) + object$tfs %*% c(
                        cov2[j,  ] * object$coef[seq(
                            (np - ncol(cov2) + 1.), np)])) * object$bfitj)
            }
        }
    }
    
    if(lps)
        lp
    else {
        lhat <- cbind(object$uftime, 1. - exp( - lhat))
        lhat <- lhat[lhat[, 1.] <= time + 1e-10,  ]
        lhat <- lhat[dim(lhat)[1.], -1.]
        lhat
    }
}
