## 
Getlim <- function(at, allow.null = FALSE, need.all = TRUE) {
    nam <- at$name[at$assume != "interaction"]
    limits <- at$limits
    values <- at$values
    
    XDATADIST <- .Options$datadist
    X <- lims <- vals <- NULL
    if (!is.null(XDATADIST) && exists(XDATADIST)) {
        X <- eval(as.name(XDATADIST))
        # X <- if (.R.) {
        #     eval(as.name(XDATADIST))
        # } else {
        #     eval(as.name(XDATADIST), local = FALSE)
        # } # 27May99  9Apr02
        lims <- X$limits
        if (is.null(lims)) {
            stop(
                "options(datadist=", XDATADIST,
                ") not created with datadist"
            )
        }
        vals <- X$values
    }
    
    if ((length(X) + length(limits)) == 0) {
        if (allow.null) {
            lims <- list()
            for (nn in nam) lims[[nn]] <- rep(NA, 7)
            lims <- structure(
                lims,
                class = "data.frame",
                row.names = c(
                    "Low:effect", "Adjust to", "High:effect", "Low:prediction",
                    "High:prediction", "Low", "High"
                )
            )
            return(list(limits = lims, values = values))
        }
        stop("no datadist in effect now or during model fit")
    }
    
    na <- if (length(limits)) {
        vapply(limits, function(x) all(is.na(x)), FALSE)
    } else {
        rep(TRUE, length(nam))
    }
    if (length(lims) && any(na)) {
        for (n in nam[na]) { # if() assumes NA stored in fit
            # for missing vars
            z <- limits[[n]]
            u <- if (match(n, names(lims), 0) > 0) lims[[n]] else NULL
            # This requires exact name match, not substring match
            if (is.null(u)) {
                if (need.all) {
                    stop(
                        "variable", n,
                        "does not have limits defined in fit or with datadist"
                    )
                } else {
                    limits[[n]] <- rep(NA, 7)
                } # Added 28 Jul 94
            } else {
                limits[[n]] <- u
            }
        }
    }
    limits <- structure(
        limits,
        class = "data.frame",
        row.names = c(
            "Low:effect", "Adjust to", "High:effect", "Low:prediction",
            "High:prediction", "Low", "High"
        )
    )
    
    if (length(vals)) {
        values <- c(
            values,
            vals[
                match(names(vals), nam, 0) > 0 & 
                    match(names(vals), names(values), 0) == 0]
        )
    } # add in values from datadist corresponding to vars in model
    # not already defined for model
    
    list(limits = limits, values = values)
}

# Function to return limits for an individual variable, given an object
# created by Getlim

Getlimi <- function(name, Limval, need.all = TRUE) {
    lim <- if (match(name, names(Limval$limits), 0) > 0) {
        Limval$limits[[name]]
    } else {
        NULL
    }
    if (is.null(Limval) || is.null(lim) || all(is.na(lim))) {
        if (need.all) {
            stop(
                "no limits defined by datadist for variable",
                name
            )
        }
        return(rep(NA, 7))
    }
    lim
}


.R. <- TRUE
