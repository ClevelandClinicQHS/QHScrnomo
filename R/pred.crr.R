pred.crr <- function(
    f.crr,
    newdata = NULL,
    time,
    lps = FALSE) {
    if (!lps) {
        if (missing(time)) {
            stop("specify expected time point.")
        }
        if (time > max(f.crr$uftime)) {
            stop("select a larger time.")
        }
        if (time < min(f.crr$uftime)) {
            stop("select a smaller time.")
        }
    }
    if (is.null(newdata)) {
        cov1 <- as.matrix(f.crr$cphdat[, -c(1,2,3)])
    } else {
        cov1 <- predictrms(f.crr$cph.f, newdata, type = "x")
    }
    np <- length(f.crr$coef)
    if (length(f.crr$tfs) <= 1.) {
        if (length(f.crr$coef) == length(cov1)) {
            lhat <- cumsum(exp(sum(cov1 * f.crr$coef)) * f.crr$bfitj)
            lp <- sum(cov1 * f.crr$coef)
        } else {
            cov1 <- as.matrix(cov1)
            lhat <- matrix(
                0.,
                nrow = length(f.crr$uftime),
                ncol = nrow(cov1)
            )
            lp <- matrix(
                0.,
                nrow = length(f.crr$uftime),
                ncol = nrow(cov1)
            )
            for (j in seq_len(nrow(cov1))) {
                lhat[, j] <- cumsum(
                    exp(sum(cov1[j, ] * f.crr$coef)) * f.crr$bfitj)
                lp[, j] <- sum(cov1[j, ] * f.crr$coef)
            }
            lp <- lp[1., ]
        }
    } else {
        if (length(f.crr$coef) == ncol(as.matrix(f.crr$tfs))) {
            if (length(f.crr$coef) == length(cov2)) {
                lhat <- cumsum(
                    exp(f.crr$tfs %*% c(cov2 * f.crr$coef)) * f.crr$bfitj)
            } else {
                cov2 <- as.matrix(cov2)
                lhat <- matrix(
                    0.,
                    nrow = length(f.crr$uftime),
                    ncol = nrow(cov1)
                )
                for (j in seq_len(nrow(cov2))) {
                    lhat[, j] <- cumsum(
                        exp(f.crr$tfs %*% c(cov2[j, ] * f.crr$coef)) * 
                            f.crr$bfitj)
                }
            }
        } else {
            if (length(f.crr$coef) == length(cov1) + length(cov2)) {
                lhat <-
                    cumsum(exp(
                        sum(cov1 * f.crr$coef[seq_len(length(cov1))]) +
                            f.crr$tfs %*% c(
                                cov2 * f.crr$coef[
                                    seq((np - length(cov2) + 1.),np)])
                    ) * f.crr$bfitj)
            } else {
                cov1 <- as.matrix(cov1)
                cov2 <- as.matrix(cov2)
                lhat <-
                    matrix(
                        0.,
                        nrow = length(f.crr$uftime),
                        ncol = nrow(cov1)
                    )
                for (j in seq_len(nrow(cov1))) {
                    lhat[, j] <-
                        cumsum(exp(
                            sum(cov1[j, ] * f.crr$coef[seq_len(ncol(cov1))]) +
                                f.crr$tfs %*% c(
                                    cov2[j, ] * 
                                        f.crr$coef[
                                            seq((np - ncol(cov2) + 1.),np)])
                        ) * f.crr$bfitj)
                }
            }
        }
    }

    if (lps) {
        lp
    } else {
        lhat <- cbind(f.crr$uftime, 1. - exp(-lhat))
        lhat <- lhat[lhat[, 1.] <= time + 1e-10, ]
        lhat <- lhat[dim(lhat)[1.], -1.]
        lhat
    }
}
