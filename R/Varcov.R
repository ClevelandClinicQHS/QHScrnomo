# Function to return variance-covariance matrix, optionally deleting
# rows and columns corresponding to parameters such as scale parameters
# in parametric survival models

Varcov <- function(object, ...) UseMethod("Varcov")

Varcov.lrm <- function(object, regcoef.only = FALSE, ...) {
    Varcov.default(object, regcoef.only, ...)
} # for fastbw etc.
Varcov.ols <- function(object, regcoef.only = FALSE, ...) {
    Varcov.default(object, regcoef.only, ...)
}
Varcov.cph <- function(object, regcoef.only = FALSE, ...) {
    Varcov.default(object, regcoef.only, ...)
}
Varcov.psm <- function(object, regcoef.only = FALSE, ...) {
    Varcov.default(object, regcoef.only, ...)
}

Varcov.default <- function(object, regcoef.only = FALSE, ...) {
    vc <- object$Varcov
    if (length(vc)) {
        if (regcoef.only) {
            return(object$var)
        } else {
            return(vc(object, which = "var"))
        }
    }
    cov <- object$var
    if (is.null(cov)) {
        stop("object fit does not have variance-covariance matrix")
    }
    if (regcoef.only) {
        p <- length(object$coefficients) # 14Sep00
        cov <- cov[seq(1, p), seq(1, p), drop = FALSE]
    }
    cov
}

Varcov.lm <- function(object, ...) {
    cof <- object$coefficients
    rinv <- solve(object$R, diag(length(cof)))
    cov <- rinv %*% t(rinv)
    cov <- sum(object$residuals^2) * cov / object$df.residual
    nm <- names(cof)
    dimnames(cov) <- list(nm, nm)
    cov
}

Varcov.glm <- function(object, ...) {
    s <- stats::summary.glm(object)
    s$cov.unscaled * s$dispersion
}

Varcov.fit.mult.impute <- function(object, ...) object$var

Varcov.multinom <- function(object, ...) stats::vcov(object)
