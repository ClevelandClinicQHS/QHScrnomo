oldUnclass <- function (x)  unclass(x)

##
value.chk <- rms:::value.chk

##
axisf <- rms:::axisf

##
is.category <- function (x)
    length(attr(x, "levels")) > 0 && mode(x) == "numeric"

##
getOldDesign <- function(fit) {
    at <- attr(fit$terms,'Design')
    if(is.null(at))
        stop('fit was not created by a Design library fitting function')
    at
}