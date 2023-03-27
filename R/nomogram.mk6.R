##' Draw a Nomogram from a regression model
##'
##' a modified version of nomogram in \code{rms} package
##'
##' @title Draw a Nomogram with modified function \code{Nomogram}
##' @usage nomogram.mk6(fit, ..., adj.to, lp = TRUE, lp.at,
##' lplabel = "Linear Predictor",
##' fun, fun.at, fun.lp.at, funlabel = "Predicted Value", fun.side,
##' interact = NULL, intercept = 1, conf.int = FALSE,
##' col.conf = c(1, 12), conf.space = c(0.08, 0.2),
##' conf.lp = c("representative", "all", "none"), est.all = TRUE,
##' abbrev = FALSE, minlength = 4, maxscale = 100, nint = 10,
##' label.every = 1, force.label = FALSE, xfrac = 0.35, cex.axis = 0.85,
##' cex.var = 1, col.grid = NULL, vnames = c("labels", "names"),
##' varname.label = TRUE, varname.label.sep = "=", ia.space = 0.7,
##' tck = NA, tcl = -0.25, lmgp = 0.4, omit = NULL, naxes, 
##' points.label = "Points",
##' total.points.label = "Total Points", total.sep.page = FALSE,
##' total.fun, verbose = FALSE, cap.labels = FALSE, total.min,
##' total.max, survtime, mikeomit = NULL)
##'
##' @param fit a regression model fit that was created with \code{library(rms)}
##' in effect, and (usually) with \code{options(datadist = "object.name")} in 
##' effect.
##' @param ...     settings of variables to use in constructing axes. 
##' If \code{datadist} was in effect, the default is to use 
##' \code{pretty(total range, nint)} for continuous variables, and the class 
##' levels for discrete ones.
##' For \code{legend.nomabbrev}, \code{\dots} specifies optional parameters to 
##' pass to \code{legend}.  Common ones are \code{bty = "n"} to suppress 
##' drawing the box.  You may want to specify a non-proportionally spaced font
##' (e.g., courier) number if abbreviations are more than one letter long.
##' This will make the abbreviation definitions line up (e.g., specify
##' \code{font = 2}, the default for courier).  Ignored for \code{print}.
##' @param adj.to   If you didn't define \code{datadist} for all predictors,
##' you will have to define adjustment settings for the undefined ones, e.g.
##' \code{adj.to= list(age = 50, sex = "female")}.
##' @param lp  Set to \code{FALSE} to suppress creation of an axis for scoring
##' \eqn{X\beta}{X beta}.
##' @param lp.at If \code{lp=TRUE}, \code{lp.at} may specify a vector of 
##' settings of \eqn{X\beta}{X beta}.
##' Default is to use \code{pretty(range of linear predictors, nint)}.
##' @param lplabel label for linear predictor axis.
##' Default is \code{"Linear Predictor"}.
##' @param fun on another axis.  If more than one transformation is plotted, 
##' put them in a list, e.g. \code{list(function(x) x/2, function(x) 2*x)}.
##' Any function values equal to \code{NA} will be ignored.
##' @param fun.at function values to label on axis.  Default \code{fun} 
##' evaluated at \code{lp.at}.   If more than one \code{fun} was specified, 
##' using a vector for \code{fun.at} will cause all functions to be evaluated 
##' at the same argument values.  To use different values, specify a list of 
##' vectors for \code{fun.at}, with elements corresponding to the different 
##' functions (lists of vectors also applies to \code{fun.lp.at} and 
##' \code{fun.side}).
##' @param fun.lp.at If you want to
##' evaluate one of the functions at a different set of linear predictor
##' values than may have been used in constructing the linear predictor axis,
##' specify a vector or list of vectors
##' of linear predictor values at which to evaluate the function.  This is
##' especially useful for discrete functions.  The presence of this attribute
##' also does away with the need for \code{nomogram} to compute numerical 
##' approximations of the inverse of the function.  It also allows the 
##' user-supplied function to return \code{factor} objects, which is useful 
##' when e.g. a single tick mark position actually represents a range.
##' If the \code{fun.lp.at} parameter is present, the \code{fun.at}
##' vector for that function is ignored.
##' @param funlabel label for \code{fun} axis.  If more than one function was 
##' given but funlabel is of length one, it will be duplicated as needed. 
##' If \code{fun} is a list of functions for which you specified names (see the
##' final example below), these names will be used as labels.
##' @param fun.side a vector or list of vectors of \code{side} parameters for 
##' the \code{axis} function for labeling function values.
##' Values may be 1 to position a tick mark label below the axis (the default),
##' or 3 for above the axis.  If for example an axis has 5 tick mark labels
##' and the second and third will run into each other, specify
##' \code{fun.side=c(1,1,3,1,1)} (assuming only one function is specified as 
##' \code{fun}).
##' @param interact When a continuous variable interacts with a discrete one, 
##' axes are constructed so that the continuous variable moves within the axis,
##' and separate axes represent levels of interacting factors.  For 
##' interactions between two continuous variables, all but the axis variable 
##' must have discrete levels defined in \code{interact}.
##' For discrete interacting factors, you may specify levels to use in
##' constructing the multiple axes.  For continuous interacting factors,
##' you must do this.  Examples: \code{interact = list(age = seq(10,70,by=10),
##' treat = c("A","B","D"))}.
##' @param intercept for models such as the ordinal logistic model with 
##' multiple intercepts, specifies which one to use in evaluating the linear 
##' predictor.
##' @param conf.int confidence levels to display for each scoring.
##' Default is \code{FALSE} to display no confidence limits.
##' Setting \code{conf.int} to \code{TRUE} is the same as
##' setting it to \code{c(0.7, 0.9)},
##' with the line segment between the 0.7 and 0.9 levels shaded using
##' gray scale.
##' @param col.conf colors corresponding to \code{conf.int}.
##' Use fractions for gray scale(for UNIX S-PLUS).
##' @param conf.space a 2-element vector with the vertical range within which 
##' to draw confidence bars, in units of 1=spacing between main bars.  Four 
##' heights are used within this range (8 for the linear predictor if more than
##' 16 unique values were evaluated), cycling them among separate confidence
##' intervals to reduce overlapping
##' @param conf.lp default is \code{"representative"} to group all linear 
##' predictors evaluated into deciles, and to show, for the linear predictor 
##' confidence intervals, only the mean linear predictor within the deciles 
##' along with the median standard error within the deciles.  Set 
##' \code{conf.lp = "none"} to suppress confidence limits for the linear 
##' predictors, and to \code{"all"} to show all confidence limits.
##' @param est.all To plot axes for only the subset of variables named in
##'  \code{\dots}, set \code{est.all = FALSE}.  Note: This option only works 
##' when zero has a special meaning for the variables that are omitted from the
##' graph
##' @param abbrev Set to \code{TRUE} to use the \code{\link{abbreviate}} 
##' function to abbreviate levels of categorical factors, both for labeling 
##' tick marks and for axis titles. If you only want to abbreviate certain 
##' predictor variables, set \code{abbrev} to a vector of character strings 
##' containing their names.
##' @param minlength \code{\link{abbreviate}} function. If you set 
##' \code{minlength = 1}, the letters of the alphabet are used to label tick 
##' marks for categorical predictors, and all letters are drawn no matter how 
##' close together they are.  For labeling axes (interaction settings), 
##' \code{minlength = 1} causes \code{minlength = 4} to be used
##' @param maxscale default maximum point score is 100
##' @param nint number of intervals to label for axes representing continuous 
##' variables.
##' See \code{\link{pretty}}.
##' @param label.every Specify \code{label.every = i} to label on every 
##' \code{i}th tick mark
##' @param force.label set to \code{TRUE} to force every tick mark intended to 
##' be labeled to have a label plotted (whether the labels run into each other 
##' or not)
##' @param xfrac fraction of horizontal plot to set aside for axis titles
##' @param cex.axis character size for tick mark labels
##' @param cex.var character size for axis titles (variable names)
##' @param col.grid If left unspecified, no vertical reference lines are drawn.
##' Specify a vector of length one (to use the same color for both minor and 
##' major reference lines) or two (corresponding to the color for the major and
##' minor divisions, respectively) containing colors, to cause vertical 
##' reference lines to the top points scale to be drawn.  For R, a good choice 
##' is \code{col.grid = gray(c(0.8, 0.95))}.
##' @param vnames By default, variable labels are used to label axes.  Set
##' \code{vnames = "names"} to instead use variable names.
##' @param varname.label In constructing axis titles for interactions, the 
##' default is to add \code{(interacting.varname = level)} on the right. 
##' Specify \code{varname.label = FALSE} to instead use \code{"(level)"}.
##' @param varname.label.sep If \code{varname.label = TRUE}, you can change the
##' separator to something other than
##' \code{=} by specifying this parameter.
##' @param ia.space When multiple axes are draw for levels of interacting 
##' factors, the default is to group combinations related to a main effect.  
##' This is done by spacing the axes for the second to last of these
##' within a group only
##' 0.7 (by default) of the way down as compared with normal space of 1 unit.
##' @param tck see \code{tck} under \code{\link{par}}
##' @param tcl length of tick marks in nomogram
##' @param lmgp spacing between numeric axis labels and axis 
##' (see \code{\link{par}} for \code{mgp})
##' @param omit vector of character strings containing names of variables for 
##' which to suppress drawing axes.  Default is to show all variables.
##' @param naxes  maximum number of axes to allow on one plot.  If the nomogram
##' requires more than one \dQuote{page}, the \dQuote{Points} axis will be 
##' repeated at the top of each page when necessary.
##' @param points.label a character string giving the axis label for the points
##' scale
##' @param total.points.label  a character string giving the axis label for 
##' the total points scale
##' @param total.sep.page set to \code{TRUE} to force the total points and 
##' later axes to be placed on a separate page
##' @param total.fun  a user-provided function that will be executed before the
##' total points axis is drawn.  Default is not toe xecute a function.  
##' This is useful e.g. when \code{total.sep.page = TRUE} and you wish to use
##' \code{locator} to find the coordinates for positioning an abbreviation 
##' legend before it's too late and a new page is started (i.e., 
##' \code{total.fun = function() print(locator(1))}).
##' @param verbose set to \code{TRUE} to get printed output detailing how tick 
##' marks are chosen and labeled for function axes.  This is useful in seeing 
##' how certain linear predictor values cannot be solved for using inverse 
##' linear interpolation on the (requested linear predictor values, function 
##' values at these lp values).  When this happens you will see \code{NA}s in 
##' the verbose output, and the corresponding tick marks will not appear in 
##' the nomogram.
##' @param cap.labels logical: should the factor labels have their first
##' letter capitalized?
##' @param total.min the minimum point for the total point axis
##' @param total.max the maxmum point for the total point axis
##' @param survtime specified survival time for the predicted survival 
##' probability
##' @param mikeomit a modified version of \code{omit}
##' @return a nomogram object
##' @importFrom Hmisc oPar setParNro strgraphwrap sedit capitalize cut2
##' @export
##' @note internal use only. please reference to \code{\link[rms]{nomogram}} 
##' details.
##' @references
##' Banks J: Nomograms. Encylopedia of Statistical Sciences, Vol 6.
##' Editors: S Kotz and NL Johnson.  New York: Wiley; 1985.
##'
##' Lubsen J, Pool J, van der Does, E: A practical device for the application
##' of a diagnostic or prognostic function.  Meth. Inform. Med. 17:127--129;
##' 1978.
##'
##' Wikipedia: Nomogram, \url{https://en.wikipedia.org/wiki/Nomogram}.
##'
##' @seealso  \code{\link{rms}}, 
##' \code{\link{plot.summary.rms}},
##' \code{\link{axis}}, \code{\link{pretty}}, \code{\link{approx}}
##'
##' @keywords graphics
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
##' ## make a cph nomogram
##' nomogram.mk6(prostate.f, survtime=120, lp=FALSE,
##' funlabel = "Predicted 10-year cumulative incidence")
##' 
nomogram.mk6 <-
    function(
        fit, ..., adj.to, lp = TRUE, lp.at,
        lplabel = "Linear Predictor",
        fun, fun.at, fun.lp.at, funlabel = "Predicted Value", fun.side,
        interact = NULL, intercept = 1, conf.int = FALSE,
        col.conf = c(1, 12), conf.space = c(0.08, 0.2),
        conf.lp = c("representative", "all", "none"), est.all = TRUE,
        abbrev = FALSE, minlength = 4, maxscale = 100, nint = 10,
        label.every = 1, force.label = FALSE, xfrac = 0.35, 
        cex.axis = 0.85, cex.var = 1, col.grid = NULL, 
        vnames = c("labels", "names"),varname.label = TRUE, 
        varname.label.sep = "=", ia.space = 0.7,tck = NA, tcl = -0.25, 
        lmgp = 0.4, omit = NULL, naxes, points.label = "Points",
        total.points.label = "Total Points", total.sep.page = FALSE,
        total.fun, verbose = FALSE, cap.labels = FALSE, total.min,
        total.max, survtime, mikeomit = NULL) {
        if (!missing(survtime)) {
            # lp <- FALSE
            surv <- Survival(fit)
            fun <- function(x) surv(survtime, x)
            assign("surv", surv)
            assign("survtime", survtime)
        }
        if (any(stats::na.omit(match(attr(fit, "class"), "lrm"))) & 
            missing(fun)) {
            # lp <- FALSE
            fun <- function(x) 1 / (1 + exp(-x))
        }
        ols <- any(stats::na.omit(match(attr(fit, "class"), "ols")))
        # if (!missing(total.min) && !missing(total.max) && lp && !ols)
        # stop("Total points control combined with lp is not setup yet.",
        #       "Pick one.")
        
        conf.lp <- match.arg(conf.lp)
        vnames <- match.arg(vnames)
        abb <- (is.logical(abbrev) && abbrev) || is.character(abbrev)
        if (is.logical(conf.int) && conf.int) {
            conf.int <- c(0.7, 0.9)
        }
        if (!is.logical(conf.int) && (length(conf.int) != length(col.conf))) {
            stop("conf.int and col.conf must have same length")
        }
        oldpar <- oPar()
        mgp <- oldpar$mgp
        mar <- oldpar$mar
        graphics::par(mgp = c(mgp[1], lmgp, mgp[3]), mar = c(
            mar[1], 1.1, mar[3],
            mar[4]
        ))
        on.exit(setParNro(oldpar))
        tck2 <- tck / 2
        tcl2 <- tcl / 2
        tck3 <- tck / 3
        tcl3 <- tcl / 3
        se <- FALSE
        if (any(conf.int > 0)) {
            se <- TRUE
            zcrit <- stats::qnorm((conf.int + 1) / 2)
            bar <- function(x, y, zcrit, se, col.conf, nlev = 4) {
                y <- rep(seq(y[1], y[2], length = nlev),length.out = length(x))
                for (j in seq(1,length(x))) {
                    xj <- x[j]
                    yj <- y[j]
                    W <- c(0, zcrit) * se[j]
                    for (i in seq(1,length(zcrit))) {
                        graphics::segments(
                            xj - W[i + 1], yj, xj - W[i], yj,
                            col = col.conf[i], lwd = 1
                        )
                        graphics::segments(
                            xj + W[i + 1], yj, xj + W[i], yj,
                            col = col.conf[i], lwd = 1
                        )
                    }
                }
            }
        }
        nfun <- if (missing(fun)) {
            0
        } else if (is.list(fun)) {
            length(fun)
        } else {
            1
        }
        if (nfun > 1 && length(funlabel) == 1) {
            funlabel <- rep(funlabel, nfun)
        }
        if (nfun > 0 && is.list(fun) && length(names(fun))) {
            funlabel <- names(fun)
        }
        if (!missing(fun.at)) {
            if (!is.list(fun.at)) {
                fun.at <- rep(list(fun.at), nfun)
            }
        }
        if (!missing(fun.lp.at)) {
            if (!is.list(fun.lp.at)) {
                fun.lp.at <- rep(list(fun.lp.at), nfun)
            }
        }
        if (!missing(fun.side)) {
            if (!is.list(fun.side)) {
                fun.side <- rep(list(fun.side), nfun)
            }
            if (any(!(unlist(fun.side) %in% c(1, 3)))) {
                stop("fun.side must contain only the numbers 1 and 3")
            }
        }
        at <- fit$Design
        if (!length(at)) {
            at <- getOldDesign(fit)
        }
        assume <- at$assume.code
        if (any(assume == 10)) {
            stop("does not currently work with matrix factors in model")
        }
        name <- at$name
        names(assume) <- name
        parms <- at$parms
        label <- if (vnames == "labels") {
            at$label
        } else {
            name
        }
        if (any(d <- duplicated(name))) {
            stop("duplicated variable names:", name[d])
        }
        label <- name
        if (vnames == "labels") {
            label <- at$label
            if (any(d <- duplicated(label))) {
                stop("duplicated variable labels:", label[d])
            }
        }
        ia <- at$interactions
        factors <- list(...)
        nf <- length(factors)
        which <- if (est.all) {
            seq(1,length(assume))[assume != 8]
        } else {
            seq(1,length(assume))[assume != 8 & assume != 9]
        }
        if (nf > 0) {
            jw <- charmatch(names(factors), name, 0)
            if (any(jw == 0)) {
                stop(
                    "factor name(s) not in the design:", 
                    names(factors)[jw == 0])
            }
            if (!est.all) {
                which <- jw
            }
        }
        Limval <- Getlim(at, allow.null = TRUE, need.all = FALSE)
        values <- Limval$values
        lims <- Limval$limits[c(6, 2, 7), , drop = FALSE]
        lims <- oldUnclass(lims)
        for (i in seq(1,length(lims))) {
            if (is.factor(lims[[i]])) {
                lims[[i]] <- as.character(lims[[i]])
            }
        }
        attr(lims, "class") <- "data.frame"
        ucat <- rep(FALSE, length(assume))
        names(ucat) <- name
        for (i in seq(1,length(assume))[assume != 5 & assume < 8]) {
            ucat[i] <- !is.null(V <- values[[name[i]]])
            if (ucat[i]) {
                parms[[name[i]]] <- V
            }
        }
        discrete <- assume == 5 | assume == 8 | ucat
        names(discrete) <- name
        nrp <- num.intercepts(fit)
        Intercept <- if (nrp > 0) {
            fit$coefficients[intercept]
        } else if (!is.null(fit$center)) {
            -fit$center
        } else {
            0
        }
        intercept.offset <- if (nrp < 2) {
            0
        } else {
            fit$coefficients[intercept] - fit$coefficients[1]
        }
        settings <- list()
        for (i in which[assume[which] < 9]) {
            ni <- name[i]
            z <- factors[[ni]]
            lz <- length(z)
            if (lz < 2) {
                settings[[ni]] <- value.chk(
                    at, i, NA, -nint, Limval,
                    type.range = "full"
                )
            } else if (lz > 0 && any(is.na(z))) {
                stop("may not specify NA as a variable value")
            }
            if (lz == 1) {
                lims[2, i] <- z
            } else if (lz > 1) {
                settings[[ni]] <- z
                if (is.null(lims[[ni]]) || is.na(lims[2, ni])) {
                    lims[[ni]] <- c(NA, z[1], NA)
                    stop(
                        "adjustment values for ", ni, 
                        " not defined in datadist;",
                        "taken to be first value specified (",
                        z[1], ")",
                    )
                }
            }
        }
        adj <- lims[2, , drop = FALSE]
        if (!missing(adj.to)) {
            for (nn in names(adj.to)) adj[[nn]] <- adj.to[[nn]]
        }
        isna <- vapply(adj, is.na, FALSE)
        if (any(isna)) {
            stop(
                "adjustment values not defined here or with datadist for",
                name[assume != 9][isna]
            )
        }
        num.lines <- 0
        space.used <- 0
        entities <- 0
        set <- list()
        nset <- character(0)
        iset <- 0
        start <- len <- NULL
        end <- 0
        main.effects <- which[assume[which] < 8]
        if (any(assume == 9)) {
            main.effects <- main.effects[
                order(10 * discrete[main.effects] + (
                    name[main.effects] %in% names(interact)))]
        }
        rel <- related.predictors(at)
        already.done <- structure(rep(FALSE, length(name)), names = name)
        for (i in main.effects) {
            nam <- name[i]
            ## changhong add
            if (already.done[nam] || (nam %in% omit) || (nam %in% mikeomit)) {
                next
            }
            r <- if (length(rel[[nam]])) {
                sort(rel[[nam]])
            } else {
                NULL
            }
            if (length(r) == 0) {
                num.lines <- num.lines + 1
                space.used <- space.used + 1
                entities <- entities + 1
                x <- list()
                x[[nam]] <- settings[[nam]]
                iset <- iset + 1
                attr(x, "info") <- list(
                    predictor = nam, effect.name = nam,
                    type = "main"
                )
                set[[iset]] <- x
                nset <- c(nset, label[i])
                start <- c(start, end + 1)
                n <- length(settings[[nam]])
                len <- c(len, n)
                end <- end + n
                NULL
            } else {
                namo <- name[r]
                s <- !(name[r] %in% names(interact))
                if (any(s)) {
                    if (is.null(interact)) {
                        interact <- list()
                    }
                    for (j in r[s]) {
                        nj <- name[j]
                        if (discrete[j]) {
                            interact[[nj]] <- parms[[nj]]
                        }
                        NULL
                    }
                    s <- !(name[r] %in% names(interact))
                }
                if (any(s)) {
                    stop(
                        "factors not defined in interact=list(...):",
                        name[r[s]]
                    )
                }
                combo <- expand.grid(interact[namo])
                oldClass(combo) <- NULL
                acombo <- combo
                if (abb) {
                    for (n in if (is.character(abbrev)) {
                        abbrev
                    } else {
                        names(acombo)
                    }) {
                        if (discrete[n]) {
                            acombo[[n]] <- abbreviate(
                                parms[[n]], 
                                minlength = if (minlength ==1) {
                                    4
                                } else {
                                    minlength
                                })[combo[[n]]]
                        }
                    }
                }
                for (n in names(combo)) {
                    if (is.factor(combo[[n]])) {
                        combo[[n]] <- as.character(combo[[n]])
                        acombo[[n]] <- as.character(acombo[[n]])
                        NULL
                    }
                }
                entities <- entities + 1
                already.done[namo] <- TRUE
                for (k in seq(1,length(combo[[1]]))) {
                    num.lines <- num.lines + 1
                    space.used <- space.used + if (k == 1) {
                        1
                    } else {
                        ia.space
                    }
                    x <- list()
                    x[[nam]] <- settings[[nam]]
                    for (nm in namo) x[[nm]] <- combo[[nm]][k]
                    iset <- iset + 1
                    set.name <- paste(nam, " (", sep = "")
                    for (j in seq(1,length(acombo))) {
                        set.name <- paste(set.name, if (varname.label) {
                            paste(namo[j], varname.label.sep, sep = "")
                        } else {
                            ""
                        }, format(acombo[[j]][k]), sep = "")
                        if (j < length(acombo)) {
                            set.name <- paste(set.name, " ", sep = "")
                        }
                    }
                    set.name <- paste(set.name, ")", sep = "")
                    ia.names <- NULL
                    for (j in r) {
                        ia.names <- c(ia.names, name[interactions.containing(
                            at,
                            j
                        )])
                    }
                    ia.names <- unique(ia.names)
                    attr(x, "info") <- list(predictor = nam, effect.name = c(
                        nam,
                        namo[assume[namo] != 8], ia.names
                    ), 
                    type = if (k == 1) {
                        "first"
                    } else {
                        "continuation"
                    })
                    set[[iset]] <- x
                    nset <- c(nset, set.name)
                    start <- c(start, end + 1)
                    n <- length(settings[[nam]])
                    len <- c(len, n)
                    end <- end + n
                    NULL
                }
                NULL
            }
        }
        xadj <- oldUnclass(Design.levels(adj, at))
        for (k in seq(1,length(xadj))) xadj[[k]] <- rep(xadj[[k]], sum(len))
        j <- 0
        for (S in set) {
            j <- j + 1
            ns <- names(S)
            nam <- names(S)
            for (k in seq(1,length(nam))) {
                xadj[[nam[k]]][start[j]:(start[j] + len[j] - 1)] <- S[[k]]
                NULL
            }
        }
        xadj <- structure(
            xadj, class = "data.frame", 
            row.names = as.character(seq(1,sum(len))))
        xx <- predictDesign(
            fit,
            newdata = xadj, type = "terms",
            center.terms = FALSE, se.fit = FALSE, kint = intercept
        )
        if (any(is.infinite(xx))) {
            stop(
                "variable limits and transformations are such that an", 
                "infinite axis value has resulted.\nRe-run specifying",
                "your own limits to variables.")
        }
        if (se) {
            xse <- predictDesign(
                fit,
                newdata = xadj, se.fit = TRUE,
                kint = intercept
            )
        }
        R <- matrix(NA, nrow = 2, ncol = length(main.effects), dimnames = list(
            NULL,
            name[main.effects]
        ))
        R[1, ] <- 1e+30
        R[2, ] <- -1e+30
        for (i in seq(1,num.lines)) {
            is <- start[i]
            ie <- is + len[i] - 1
            s <- set[[i]]
            setinfo <- attr(s, "info")
            nam <- setinfo$effect.name
            xt <- xx[is:ie, nam]
            if (length(nam) > 1) {
                xt <- apply(xt, 1, sum)
            }
            set[[i]]$Xbeta <- xt
            r <- range(xt)
            pname <- setinfo$predictor
            R[1, pname] <- min(R[1, pname], r[1])
            R[2, pname] <- max(R[2, pname], r[2])
            if (se) {
                set[[i]]$Xbeta.whole <- xse$linear.predictors[is:ie]
                set[[i]]$se.fit <- xse$se.fit[is:ie]
                NULL
            }
            NULL
        }
        R <- R[, R[1, ] < 1e+30, drop = FALSE]
        sc <- maxscale / max(R[2, ] - R[1, ])
        xl <- -xfrac * maxscale
        Intercept <- Intercept + sum(R[1, ])
        if (missing(naxes)) {
            naxes <- if (total.sep.page) {
                max(space.used + 1, nfun + lp + 1)
            } else {
                space.used + 1 + nfun + lp + 1
            }
        }
        Format <- function(x) {
            f <- character(l <- length(x))
            for (i in seq(1,l)) f[i] <- format(x[i])
            f
        }
        newpage <- function(
            naxes, xl, maxscale, cex.var, nint, space.used,
            col.grid, cex.axis, tck, tck2, tcl, tcl2, 
            label.every, force.label, points = TRUE, 
            points.label = "Points", usr) {
            y <- naxes - 1
            plot(
                0, 0,
                xlim = c(xl, maxscale), ylim = c(0, y), type = "n",
                axes = FALSE, xlab = "", ylab = ""
            )
            if (!missing(usr)) {
                graphics::par(usr = usr)
            }
            if (!points) {
                return(y + 1)
            }
            ax <- c(0, maxscale)
            graphics::text(xl, y, points.label, adj = 0, cex = cex.var)
            x <- pretty(ax, n = nint)
            dif <- x[2] - x[1]
            x2 <- seq((x[1] + x[2]) / 2, max(x), by = dif)
            x2 <- sort(c(x2 - dif / 4, x2, x2 + dif / 4))
            if (length(col.grid)) {
                graphics::segments(
                    x, y, x, y - space.used,
                    col = col.grid[1],
                    lwd = 1
                )
                graphics::segments(
                    x2, y, x2, y - space.used,
                    col = col.grid[-1],
                    lwd = 1
                )
            }
            axisf(
                3,
                at = x, pos = y, cex = cex.axis, tck = tck,
                tcl = tcl, label.every = label.every, 
                force.label = force.label, padj = 0
            )
            axisf(
                3,
                at = x2, labels = FALSE, pos = y, tck = tck2,
                tcl = tcl2, cex = cex.axis
            )
            y
        }
        y <- newpage(
            naxes, xl, maxscale, cex.var, nint, space.used,
            col.grid, cex.axis, tck, tck2, tcl, tcl2,
            label.every = label.every,
            force.label = force.label, points.label = points.label
        )
        i <- 0
        Abbrev <- list()
        for (S in set) {
            i <- i + 1
            setinfo <- attr(S, "info")
            type <- setinfo$type
            y <- y - (if (type == "continuation") {
                ia.space
            } else {
                1
            })
            if (y < -0.05) {
                y <- newpage(
                    naxes, xl, maxscale, cex.var, nint,
                    space.used, col.grid, cex.axis, tck, tck2, tcl,
                    tcl2,
                    label.every = label.every, force.label = force.label,
                    points.label = points.label
                ) - (if (type == "continuation") {
                    ia.space
                } else {
                    1
                })
            }
            graphics::text(
                xl, y, paste(
                    strgraphwrap(
                        nset[[i]], 
                        abs(xl), cex = cex.var),
                    collapse = "\n"
                ), adj = 0, cex = cex.var)
            x <- S[[1]]
            nam <- names(S)[1]
            fx <- if (is.character(x)) {
                x
            } else {
                sedit(Format(x), " ", "")
            }
            if (abb && discrete[nam] && 
                (is.logical(abbrev) || nam %in% abbrev)) {
                old.text <- fx
                fx <- if (abb && minlength == 1) {
                    letters[seq(1,length(fx))]
                } else {
                    abbreviate(fx, minlength = minlength)
                }
                Abbrev[[nam]] <- list(abbrev = fx, full = old.text)
            }
            j <- match(nam, name, 0)
            #if (any(j == 0)) stop("program logic error 1")
            is <- start[i]
            ie <- is + len[i] - 1
            xt <- (S$Xbeta - R[1, nam]) * sc
            set[[i]]$points <- xt
            r <- rle(xt)
            if (any(r$length > 1)) {
                is <- 1
                for (j in r$length) {
                    ie <- is + j - 1
                    if (j > 1) {
                        fx[ie] <- if (discrete[nam] || ie < length(xt)) {
                            paste(fx[is], "-", fx[ie], sep = "")
                        } else {
                            paste(fx[is], "+", sep = "")
                        }
                        fx[is:(ie - 1)] <- ""
                        xt[is:(ie - 1)] <- NA
                    }
                    is <- ie + 1
                }
                fx <- fx[!is.na(xt)]
                xt <- xt[!is.na(xt)]
            }
            side <- c(1, 3)
            padj <- c(1, 0)
            new.mgp <- vector(mode = "list", 2)
            new.mgp[[2]] <- c(0, lmgp, 0)
            new.mgp[[1]] <- new.mgp[[2]] - c(0, 0.6, 0)
            ch <- if (length(xt) > 2) {
                c(FALSE, FALSE, diff(diff(xt) > 0) != 0)
            } else {
                rep(FALSE, length(xt))
            }
            if (discrete[nam] && length(xt) > 1) {
                j <- order(xt)
                graphics::lines(range(xt), rep(y, 2))
                for (k in c(1,2)) {
                    is <- j[seq(k, length(j), by = 2)]
                    new.labs <- if (cap.labels) {
                        capitalize(fx[is])
                    } else {
                        fx[is]
                    }
                    axisf(
                        side[k],
                        at = xt[is], labels = new.labs,
                        pos = y, cex = cex.axis, tck = tck, tcl = tcl,
                        force.label = force.label || 
                            (abb && minlength == 1 && (
                                is.logical(abbrev) || nam %in% abbrev)),
                        disc = TRUE, mgp = new.mgp[[k]], padj = padj[k]
                    )
                    if (se) {
                        bar(
                            xt[is], if (k == 1) {
                                y - conf.space - 0.32
                            } else {
                                y + conf.space + 0.32
                            }, zcrit, sc * S$se.fit[is],
                            col.conf
                        )
                    }
                }
            } else if (!any(ch)) {
                axisf(
                    1,
                    at = xt, labels = fx, pos = y, cex = cex.axis,
                    tck = tck, tcl = tcl, mgp = new.mgp[[1]], 
                    label.every = label.every, force.label = force.label, 
                    disc = discrete[nam], padj = padj[1]
                )
                if (se) {
                    bar(
                        xt, y + conf.space, zcrit, sc * S$se.fit,
                        col.conf
                    )
                }
            } else {
                graphics::lines(range(xt), rep(y, 2))
                j <- seq(1,length(ch))[ch]
                if (max(j) < length(ch)) {
                    j <- c(j, length(ch) + 1)
                }
                flag <- 1
                is <- 1
                for (k in j) {
                    ie <- k - 1
                    axisf(
                        side[flag],
                        at = xt[is:ie], labels = fx[is:ie],
                        pos = y, cex = cex.axis, tck = tck, tcl = tcl,
                        label.every = label.every, force.label = force.label,
                        mgp = new.mgp[[flag]], disc = discrete[nam],
                        padj = padj[flag]
                    )
                    if (se) {
                        bar(
                            xt[is:ie], if (side[flag] == 1) {
                                y - conf.space - 0.32
                            } else {
                                y + conf.space + 0.32
                            }, zcrit, sc * S$se.fit[is:ie],
                            col.conf
                        )
                    }
                    flag <- if (flag == 2) {
                        1
                    } else {
                        2
                    }
                    is <- ie + 1
                }
            }
        }
        # browser()
        if (missing(lp.at)) {
            xb <- fit$linear.predictors
            if (!length(xb)) {
                xb <- fit$fitted.values
            }
            if (!length(xb)) {
                xb <- fit$fitted
            }
            if (!length(xb)) {
                stop(
                    "lp.at not given and fit did not store", 
                    "linear.predictors or fitted.values")
            }
            if (nrp > 1) {
                xb <- xb + intercept.offset
            }
            lp.at <- pretty(range(xb), n = nint)
        }
        sum.max <- if (entities == 1) {
            maxscale
        } else {
            max(maxscale, sc * max(lp.at - Intercept))
        }
        x <- pretty(c(0, sum.max), n = nint)
        new.max <- max(x)
        xl.old <- xl
        xl <- -xfrac * new.max
        u <- graphics::par()$usr
        if (!missing(total.fun)) {
            total.fun()
        }
        usr <- c(xl * u[1] / xl.old, new.max * u[2] / maxscale, u[3:4])
        graphics::par(usr = usr)
        x.double <- seq(x[1], new.max, by = (x[2] - x[1]) / 5)
        y <- y - 1
        if (y < -0.05 || total.sep.page) {
            y <- newpage(
                naxes, xl, maxscale, cex.var, nint, space.used,
                col.grid, cex.axis, tck, tck2, tcl, tcl2,
                label.every = label.every,
                force.label = force.label, points = FALSE, usr = usr
            ) -
                1
        }
        graphics::text(xl, y, total.points.label, adj = 0, cex = cex.var)
        # changhong edit here
        if (missing(total.max) || missing(total.min)) {
            axisf(
                1,
                at = x, pos = y, cex = cex.axis, tck = tck, tcl = tcl,
                label.every = label.every, force.label = force.label,
                mgp = c(0, lmgp - 0.6, 0), padj = 1
            )
            axisf(
                1,
                at = x.double, labels = FALSE, pos = y, tck = tck2,
                tcl = tcl2, cex = cex.axis
            )
        } else {
            axisf(
                1,
                at = x, pos = y, 
                labels = round(seq(
                    total.min,total.max,
                    by = (total.max - total.min) / (length(x) - 1)
                )),
                cex = cex.axis, tck = tck, tcl = tcl,
                label.every = label.every, force.label = force.label,
                mgp = c(0, lmgp - 0.6, 0), padj = 1
            )
            axisf(
                1,
                at = x.double, labels = FALSE, pos = y, tck = tck2,
                tcl = tcl2, cex = cex.axis
            )
        }
        iset <- iset + 1
        nset <- c(nset, "total.points")
        set[[iset]] <- list(x = x, y = y)
        # browser()
        if (lp) {
            x2 <- seq(lp.at[1], max(lp.at), by = (lp.at[2] - lp.at[1]) / 2)
            scaled.x <- (lp.at - Intercept) * sc
            scaled.x2 <- (x2 - Intercept) * sc
            y <- y - 1
            if (y < -0.05) {
                y <- newpage(
                    naxes, xl, maxscale, cex.var, nint,
                    space.used, col.grid, cex.axis, tck, tck2, tcl,
                    tcl2,
                    label.every = label.every, force.label = force.label,
                    points = FALSE, usr = usr
                ) - 1
            }
            graphics::text(xl, y, lplabel, adj = 0, cex = cex.var)
            # changhong edit here on 05/12/2010
            if (missing(total.max) || missing(total.min)) {
                axisf(
                    1,
                    at = scaled.x, labels = Format(lp.at), pos = y,
                    cex = cex.axis, tck = tck, tcl = tcl, 
                    label.every = label.every,
                    force.label = force.label, mgp = c(
                        0, lmgp - 0.6,
                        0
                    ), padj = 1
                )
            } else {
                axisf(
                    1,
                    at = ((scaled.x - total.min) * (
                        max(x) - min(x))) / (total.max - total.min), 
                    labels = Format(lp.at), pos = y,
                    cex = cex.axis, tck = tck, tcl = tcl, 
                    label.every = label.every,
                    force.label = force.label, mgp = c(
                        0, lmgp - 0.6,
                        0
                    ), padj = 1
                )
            }
            axisf(
                1,
                at = scaled.x2, labels = FALSE, tck = tck2,
                tcl = tcl2, pos = y, cex = cex.axis
            )
            iset <- iset + 1
            nset <- c(nset, "lp")
            set[[iset]] <- list(x = scaled.x, y = y, x.real = lp.at)
            if (se && conf.lp != "none") {
                xxb <- NULL
                xse <- NULL
                for (S in set) {
                    xxb <- c(xxb, S$Xbeta.whole)
                    xse <- c(xse, S$se.fit)
                }
                i <- order(xxb)
                if (length(xxb) < 16 | conf.lp == "representative") {
                    nlev <- 4
                    w <- 1
                } else {
                    nlev <- 8
                    w <- 2
                }
                if (conf.lp == "representative") {
                    deciles <- cut2(xxb[i], g = 10)
                    mean.xxb <- tapply(xxb[i], deciles, mean)
                    median.se <- tapply(xse[i], deciles, stats::median)
                    bar((mean.xxb - Intercept) * sc, y + c(
                        conf.space[1],
                        conf.space[1] + w * diff(conf.space)
                    ), zcrit,
                    sc * median.se, col.conf,
                    nlev = nlev
                    )
                } else {
                    bar((xxb[i] - Intercept) * sc, y + c(
                        conf.space[1],
                        conf.space[1] + w * diff(conf.space)
                    ), zcrit,
                    sc * xse[i], col.conf,
                    nlev = nlev
                    )
                }
            }
        }
        # browser()
        if (nfun > 0) {
            if (!is.list(fun)) {
                fun <- list(fun)
            }
            i <- 0
            for (func in fun) {
                i <- i + 1
                if (!missing(fun.lp.at)) {
                    xseq <- fun.lp.at[[i]]
                    fat <- func(xseq)
                    w <- xseq
                } else {
                    if (missing(fun.at)) {
                        fat <- pretty(func(range(lp.at)), n = nint)
                    } else {
                        fat <- fun.at[[i]]
                    }
                    if (verbose) {
                        cat(
                            "Function", i, 
                            "values at which to place tick marks:\n")
                        print(fat)
                    }
                    xseq <- seq(min(lp.at), max(lp.at), length = 1000)
                    fu <- func(xseq)
                    s <- !is.na(fu)
                    w <- stats::approx(fu[s], xseq[s], fat, ties = mean)$y
                    if (verbose) {
                        cat("Estimated inverse function values (lp):\n")
                        print(w)
                    }
                }
                s <- !(is.na(w) | is.na(fat))
                w <- w[s]
                fat <- fat[s]
                fat.orig <- fat
                fat <- if (is.category(fat)) {
                    as.character(fat)
                } else {
                    Format(fat)
                }
                scaled <- (w - Intercept) * sc
                y <- y - 1
                if (y < -0.05) {
                    y <- newpage(
                        naxes, xl, maxscale, cex.var, nint,
                        space.used, col.grid, cex.axis, tck, tck2,
                        tcl, tcl2,
                        label.every = label.every, force.label = force.label,
                        points = FALSE, usr = usr
                    ) - 1
                }
                graphics::text(xl, y, funlabel[i], adj = 0, cex = cex.var)
                sides <- if (missing(fun.side)) {
                    rep(1, length(fat))
                } else {
                    (fun.side[[i]])[s]
                }
                if (length(sides) != length(fat)) {
                    stop(
                        "fun.side vector not same length as",
                        "fun.at or fun.lp.at")
                }
                
                # changhong edit here
                for (jj in seq(1,length(fat))) {
                    if (missing(total.max) || missing(total.min)) {
                        graphics::axis(
                            sides[jj],
                            at = scaled[jj], label = fat[jj],
                            pos = y, cex.axis = cex.axis, tck = tck, tcl = tcl,
                            mgp = if (sides[jj] == 1) {
                                c(0, lmgp - 0.6, 0)
                            } else {
                                c(0, lmgp, 0)
                            }, padj = if (sides[jj] == 1) {
                                1
                            } else {
                                0
                            }
                        )
                        graphics::lines(range(scaled), rep(y, 2))
                    } else {
                        graphics::axis(
                            sides[jj],
                            at = ((
                                scaled[jj] - total.min) * (
                                    max(x) - min(x))) / 
                                (total.max - total.min), 
                            label = fat[jj],
                            pos = y, cex.axis = cex.axis, 
                            tck = tck, tcl = tcl,
                            mgp = if (sides[jj] == 1) {
                                c(0, lmgp - 0.6, 0)
                            } else {
                                c(0, lmgp, 0)
                            }, padj = if (sides[jj] == 1) {
                                1
                            } else {
                                0
                            }
                        )
                        graphics::lines(((
                            range(scaled) - total.min) * 
                                (max(x) - min(x))) / 
                                (total.max - total.min), rep(y, 2))
                    }
                }
                
                iset <- iset + 1
                nset <- c(nset, funlabel[i])
                set[[iset]] <- list(x = scaled, y = y, x.real = fat.orig)
            }
        }
        names(set) <- nset
        set$abbrev <- Abbrev
        oldClass(set) <- "nomogram"
        invisible(set)
    }
