##' Generate an equation to calculate X beta from a crr model fit
##' If specify a time point, the function also generates the subcumulative rate
##' at the time point.
##'
##' f should be fitted by the function \code{\link[QHScrnomo]{crr.fit}}
##'
##' @title generate prediction equation for a competing risks regression models
##' @param f a model fit from the competing risks regression.
##' @param time time point
##' @param baseonly logical variable. If true, only base survival probability
##'   will be printed.
##' @param file A connection, or a character string naming the file to print to
##' @param append logical. Only used if the argument file is the name of file.
##'   If TRUE output will be appended to file; otherwise, it will overwrite the
##'   content of file
##' @return \item{out }{a character vector that can be output as a formula by
##'   function cat} \item{Rout }{same as out except replacing
##'   \code{"max"},\code{"min"},\code{"="} and \code{"**"} with
##'   \code{"pmax"},\code{"pmin"}, \code{"=="} and
##'   \eqn{\mbox{\textasciicircum}}{"^"} respectively so that the formula can
##'   be pasted to R session and compute \code{X} beta directly without any
##'   further modification}
##' @author changhong
##' @importFrom Hmisc Function sedit
##' @export
##' @seealso
##'   \code{\link[rms]{sascode}},\code{\link[rms]{Function}},
##'   \code{\link[QHScrnomo]{crr.fit}}
##' @keywords regression survival
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
##' sas.cmprsk(prostate.crr, time = 60)
##'
sas.cmprsk <- function(
    f, time = NA, baseonly = FALSE,
    file = "", append = FALSE) {
    call <- match.call()
    call <- as.list(call)
    if (!any(regexpr("cmprsk", attr(f, "class")) != -1)) {
        stop(
            call$f, " should be a model fitted by function crr.fit !")
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
                    to = c("max", "min", "=", "=", "**"),
                    wild.literal = TRUE
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
