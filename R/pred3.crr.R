##' Predict cumulative incidence used internally
##' an internal function.
##'
##' Internally used only.
##' @title Predict cumulative incidence used internally
##' @param z  the fitter crr model
##' @param cov1 covarite matrix 1
##' @param cov2 covariate matrix 2
##' @param time time point at which the prediction will make
##' @param lps logical flag. if the liner predictor be generated.
##' @return A list. See \code{\link[cmprsk]{crr}} for details.
##' @note an internal function called by \code{\link{crr.fit}}.
##' @author changhong
##' @references NULL
##'
##'
pred3.crr <- function(z, cov1, cov2, time, lps = FALSE) {
    np <- length(z$coef)
    if (length(z$tfs) <= 1.) {
        if (length(z$coef) == length(cov1)) {
            lhat <- cumsum(exp(sum(cov1 * z$coef)) * z$bfitj)
            lp <- sum(cov1 * z$coef)
        } else {
            cov1 <- as.matrix(cov1)
            lhat <- matrix(0., nrow = length(z$uftime), ncol = nrow(cov1))
            lp <- matrix(0., nrow = length(z$uftime), ncol = nrow(cov1))
            for (j in 1.:nrow(cov1)) {
                lhat[, j] <- cumsum(exp(sum(cov1[j, ] * z$coef)) * z$bfitj)
                lp[, j] <- sum(cov1[j, ] * z$coef)
            }
            lp <- lp[1., ]
        }
    } else {
        if (length(z$coef) == ncol(as.matrix(z$tfs))) {
            if (length(z$coef) == length(cov2)) {
                lhat <- cumsum(exp(z$tfs %*% c(cov2 * z$coef)) * z$bfitj)
            } else {
                cov2 <- as.matrix(cov2)
                lhat <- matrix(
                    0.,
                    nrow = length(z$uftime),
                    ncol = nrow(cov1)
                )
                for (j in 1.:nrow(cov2)) {
                    lhat[, j] <- cumsum(exp(z$tfs %*% c(
                        cov2[j, ] * z$coef
                    )) * z$bfitj)
                }
            }
        } else {
            if (length(z$coef) == length(cov1) + length(cov2)) {
                lhat <- cumsum(exp(sum(cov1 * z$coef[1.:length(
                    cov1
                )]) + z$tfs %*% c(cov2 * z$coef[
                    (np - length(cov2) + 1.):np
                ])) * z$
                    bfitj)
            } else {
                cov1 <- as.matrix(cov1)
                cov2 <- as.matrix(cov2)
                lhat <- matrix(
                    0.,
                    nrow = length(z$uftime),
                    ncol = nrow(cov1)
                )
                for (j in 1.:nrow(cov1)) {
                    lhat[, j] <- 
                        cumsum(
                            exp(sum(
                                cov1[j, ] * z$coef[1.:ncol(cov1)]) + z$tfs %*% 
                                    c(cov2[j, ] * z$coef[seq(
                                        (np - ncol(cov2) + 1.), np)])) * 
                                z$bfitj)
                }
            }
        }
    }
    lhat <- cbind(z$uftime, 1. - exp(-lhat))
    lhat <- lhat[lhat[, 1.] <= time + 1e-10, ]
    lhat <- lhat[dim(lhat)[1.], -1.]
    if (lps) {
        lp
    } else {
        lhat
    }
}
