##' @importFrom Hmisc Function sedit
sas.crr <- function(f, time = NA, baseonly = FALSE, file = "", append = FALSE) {
    call <- match.call()
    call <- as.list(call)
    if (!any(regexpr("crr", attr(f, "class")) != -1)) {
        stop(call$f, " should be a model fitted by function crr.fit !")
    }
    # removed fill=F from command line
    if (!is.na(time)) {
        lhat <- cbind(f$uftime, exp(-cumsum(f$bfitj)))
        if (time > max(lhat[, 1])) stop("choose a smaller predicting time")
        lhat <- lhat[lhat[, 1.] <= time + 1e-10, ]
        lhat <- lhat[dim(lhat)[1.], 2.]
        if (file != "") {
            cat("Base is ", lhat, "\n", file = file, append = append)
        } else {
            cat("Base is ", lhat, "\n")
        }
    }
    if (!baseonly) {
        cph.f <- f$cph.f
        cph.f$coefficients <- f$coef
        cph.f$var <- f$var
        funout <- as.list(Function(cph.f))
        n <- length(funout)
        funout <- funout[[n]]
        fn <- length(funout)
        if (fn > 2) {
            for (j in 2:(fn - 1)) {
                trans <- deparse(funout[[j]])
                trans <- sedit(
                    trans,
                    from = c("pmax", "pmin", "<-", "==", "^"),
                    to = c("max", "min", "=", "=", "**"), wild.literal = TRUE
                )
                cat(trans, "\n")
            }
        }
        out <- (funout)[[fn]]
        out2 <- deparse(out)
        pos <- gregexpr("[+-]", out2[1])[[1]]
        pos2 <- pos[pos > 1][1]
        # remove the intercept that is the center of cph model
        # but unapplicable for crr
        out2[1] <- ifelse(
            substring(out2[1], pos2, pos2) == "+",
            substring(out2[1], pos2 + 1),
            substring(out2[1], pos2)
        )
        out3 <- out2
        out4 <- sedit(
            out3,
            from = c("pmax", "pmin", "<-", "==", "^"),
            to = c("max", "min", "=", "=", "**"), wild.literal = TRUE
        )
        if (file == "") {
            cat(out4, sep = "\n")
        } else {
            cat(out4, sep = "\n", file = file, append = TRUE)
        }
    }
}
