##'  Draws a partial nomogram that can be used to manually obtain predicted
##'  values from a regression model that was fitted with \code{rms} in effect.
##'
##' The nomogram does not have lines representing sums, but it has a reference
##' line for reading scoring points (default range 0--100).  Once the reader
##' manually totals the points, the predicted values can be read at the bottom.
##' Non-monotonic transformations of continuous variables are handled (scales
##' wrap around), as are transformations which have flat sections (tick marks
##' are labeled with ranges).
##' @title Draw a Nomogram
##' @param fit a regression model fit that was created with \code{library(rms)}
##' in effect, and (usually) with \code{options(datadist = "object.name")} in
##' effect.
##' @param ... settings of variables to use in constructing axes.
##' If \code{datadist} was in effect, the default is to use
##' \code{pretty(total range, nint)} for continuous variables, and the class
##' levels for discrete ones.
##' For \code{legend.nomabbrev}, \code{\dots} specifies optional parameters
##' to pass to \code{legend}.  Common ones are \code{bty = "n"} to suppress
##' drawing the box. You may want to specify a non-proportionally spaced font
##' (e.g., courier) number if abbreviations are more than one letter long.
##' This will make the abbreviation definitions line up (e.g., specify
##'  \code{font = 2}, the default for courier).  Ignored for \code{print}.
##' @note internal use only
##' @author   Frank Harrell\cr
##' Department of Biostatistics\cr
##' Vanderbilt University\cr
##' \email{f.harrell@@vanderbilt.edu}

## nomogram <- function(fit, ...) UseMethod("nomogram")

##'
##' Draw a Competing Risks Nomogram
##'
##'
##' Draws a partial nomogram adjusting for competing risks for a cox ph
##' survival model.
##'
##' @title Draw nomogram for competing risks regression models
##' @param fit a competing risks regression model fit that was created with
##'   function \code{\link[QHScrnomo]{crr.fit}}.
##' @param failtime the expected failure time for calculating cumalative
##'   incidence.
##' @param ci logical flag to output cumulative incidence or event free
##'   probability if setting \code{FALSE}.
##' @param ...  settings of variables to use in constructing axes.
##' If datadist was in effect, the default is to use pretty(total range, nint)
##' for continuous variables, and the class levels for discrete ones.
##' For legend.nomabbrev, ... specifies optional parameters to pass to legend.
##' Common ones are bty = "n" to suppress drawing the box. You may want to
##' specify a non-proportionally spaced font (e.g., courier) number if
##' abbreviations are more than one letter long. This will make the abbreviation
##' definitions line up (e.g., specify font = 2, the default for courier).
##' Ignored for print.
##' @param adj.to If you didn't define \code{datadist} for all predictors, you
##'   will have to define adjustment settings for the undefined ones, e.g.
##'   \code{adj.to=list(age=50, sex="female")}.
##' @param lp Set to \code{FALSE} to suppress creation of an axis for scoring
##'   \eqn{X\beta}{X beta}.
##' @param lp.at If \code{lp=TRUE}, \code{lp.at} may specify a vector of
##'   settings of \eqn{X\beta}{X beta}. Default is to use \code{pretty(range of
##'   linear predictors, nint)}.
##' @param lplabel label for linear predictor axis.  Default is \code{"Linear
##'   Predictor"}.
##' @param fun.at function values to label on axis.  Default \code{fun}
##'   evaluated at \code{lp.at}.  If more than one \code{fun} was specified,
##'   using a vector for \code{fun.at} will cause all functions to be evaluated
##'   at the same argument values.  To use different values, specify a list of
##'   vectors for \code{fun.at}, with elements corresponding to the different
##'   functions (lists of vectors also applies to \code{fun.lp.at} and
##'   \code{fun.side}).
##' @param fun.lp.at If you want to evaluate one of the functions at a
##'   different set of linear predictor values than may have been used in
##'   constructing the linear predictor axis, specify a vector or list of
##'   vectors of linear predictor values at which to evaluate the function.
##'   This is especially useful for discrete functions.  The presence of this
##'   attribute also does away with the need for \code{nomogram} to compute
##'   numerical approximations of the inverse of the function.  It also allows
##'   the user-supplied function to return \code{factor} objects, which is
##'   useful when e.g. a single tick mark position actually represents a range.
##'   If the \code{fun.lp.at} parameter is present, the \code{fun.at} vector
##'   for that function is ignored.
##' @param funlabel label for \code{fun} axis.  If more than one function was
##'   given but funlabel is of length one, it will be duplicated as needed.  If
##'   \code{fun} is a list of functions for which you specified names (see the
##'   final example below), these names will be used as labels.
##' @param fun.side a vector or list of vectors of \code{side} parameters for
##'   the \code{axis} function for labeling function values. Values may be 1 to
##'   position a tick mark label below the axis (the default), or 3 for above
##'   the axis.  If for example an axis has 5 tick mark labels and the second
##'   and third will run into each other, specify \code{fun.side=c(1,1,3,1,1)}
##'   (assuming only one function is specified as \code{fun}).
##' @param interact When a continuous variable interacts with a discrete one,
##'   axes are constructed so that the continuous variable moves within the
##'   axis, and separate axes represent levels of interacting factors.  For
##'   interactions between two continuous variables, all but the axis variable
##'   must have discrete levels defined in \code{interact}.  For discrete
##'   interacting factors, you may specify levels to use in constructing the
##'   multiple axes.  For continuous interacting factors, you must do this.
##'   Examples: \code{interact=list(age=seq(10,70,by=10),
##'   treat=c("A","B","D"))}.
##' @param intercept for models such as the ordinal logistic model with
##'   multiple intercepts, specifies which one to use in evaluating the linear
##'   predictor.
##' @param conf.int confidence levels to display for each scoring.  Default is
##'   \code{FALSE} to display no confidence limits.  Setting \code{conf.int} to
##'   \code{TRUE} is the same as setting it to \code{c(0.7, 0.9)}, with the
##'   line segment between the 0.7 and 0.9 levels shaded using gray scale.
##' @param col.conf colors corresponding to \code{conf.int}.  Use fractions for
##'   gray scale (for UNIX S-PLUS).
##' @param conf.space a 2-element vector with the vertical range within which
##'   to draw confidence bars, in units of 1=spacing between main bars.  Four
##'   heights are used within this range (8 for the linear predictor if more
##'   than 16 unique values were evaluated), cycling them among separate
##'   confidence intervals to reduce overlapping.
##' @param conf.lp default is \code{"representative"} to group all linear
##'   predictors evaluated into deciles, and to show, for the linear predictor
##'   confidence intervals, only the mean linear predictor within the deciles
##'   along with the median standard error within the deciles.  Set
##'   \code{conf.lp="none"} to suppress confidence limits for the linear
##'   predictors, and to \code{"all"} to show all confidence limits.
##' @param est.all To plot axes for only the subset of variables named in
##'   \code{\dots{}}, set \code{est.all=FALSE}.  Note: This option only works
##'   when zero has a special meaning for the variables that are omitted from
##'   the graph.
##' @param abbrev Set to \code{TRUE} to use the \code{abbreviate} function to
##'   abbreviate levels of categorical factors, both for labeling tick marks
##'   and for axis titles. If you only want to abbreviate certain predictor
##'   variables, set \code{abbrev} to a vector of character strings containing
##'   their names.
##' @param minlength applies if \code{abbrev=TRUE}.  Is the minimum
##'   abbreviation length passed to the \code{abbreviate} function.  If you set
##'   \code{minlength=1}, the letters of the alphabet are used to label tick
##'   marks for categorical predictors, and all letters are drawn no matter how
##'   close together they are.  For labeling axes (interaction settings),
##'   \code{minlength=1} causes \code{minlength=4} to be used.
##' @param maxscale default maximum point score is 100
##' @param nint number of intervals to label for axes representing continuous
##'   variables. See \code{pretty}.
##' @param label.every Specify \code{label.every=i} to label on every
##'   \code{i}th tick mark.
##' @param force.label set to \code{TRUE} to force every tick mark intended to
##'   be labeled to have a label plotted (whether the labels run into each
##'   other or not)
##' @param xfrac fraction of horizontal plot to set aside for axis titles
##' @param cex.axis character size for tick mark labels
##' @param cex.var character size for axis titles (variable names)
##' @param col.grid If \code{col.grid=1}, no gray scale is used, but an
##'   ordinary line is drawn.  If \code{0<col.grid<1}, a \code{col} (gray
##'   scale) of \code{col.grid} is used to draw vertical reference lines for
##'   major axis divisions and \code{col.grid/2} for minor divisions. The
##'   default is \code{col.grid=FALSE}, i.e., reference lines are omitted.
##'   Specifying \code{col.grid=TRUE} is the same as specifying a gray scale
##'   level of \code{col.grid=.2} (5 for Windows S-PLUS).
##' @param vnames By default, variable labels are used to label axes.  Set
##'   \code{vnames="names"} to instead use variable names.
##' @param varname.label In constructing axis titles for interactions, the
##'   default is to add \code{"(interacting.varname=level)"} on the right.
##'   Specify \code{varname.label=FALSE} to instead use \code{"(level)"}.
##' @param varname.label.sep If \code{varname.label=TRUE}, you can change the
##'   separator to something other than \code{=} by specifying this parameter.
##' @param ia.space When multiple axes are draw for levels of interacting
##'   factors, the default is to group combinations related to a main effect.
##'   This is done by spacing the axes for the second to last of these within a
##'   group only 0.7 (by default) of the way down as compared with normal space
##'   of 1 unit.
##' @param tck see \code{tck} under \code{par}
##' @param lmgp spacing between numeric axis labels and axis (see \code{par}
##'   for \code{mgp})
##' @param omit vector of character strings containing names of variables for
##'   which to suppress drawing axes.  Default is to show all variables.
##' @param naxes maximum number of axes to allow on one plot.  If the nomogram
##'   requires more than one "page", the "Points" axis will be repeated at the
##'   top of each page when necessary.
##' @param points.label a character string giving the axis label for the points
##'   scale
##' @param total.points.label a character string giving the axis label for the
##'   total points scale
##' @param total.sep.page set to \code{TRUE} to force the total points and
##'   later axes to be placed on a separate page
##' @param total.fun a user-provided function that will be executed before the
##'   total points axis is drawn.  Default is not to execute a function.  This
##'   is useful e.g. when \code{total.sep.page=TRUE} and you wish to use
##'   \code{locator} to find the coordinates for positioning an abbreviation
##'   legend before it's too late and a new page is started (i.e.,
##'   \code{total.fun=function()print(locator(1))}).
##' @param verbose set to \code{TRUE} to get printed output detailing how tick
##'   marks are chosen and labeled for function axes.  This is useful in seeing
##'   how certain linear predictor values cannot be solved for using inverse
##'   linear interpolation on the (requested linear predictor values, function
##'   values at these lp values).  When this happens you will see \code{NA}s in
##'   the \code{verbose} output, and the corresponding tick marks will not
##'   appear in the nomogram.
##' @param total.min Setting the minimal value in the total point axis on the
##'   nomogram.
##' @param total.max Setting the maximal value in the total point axis.
##' @param mikeomit The predictor variables specified by their names here will
##'   not be shown in the nomogram.  The predicted outcome based on this
##'   reduced nomogram would be the same as if users were using the full
##'   version of the nomogram by entering the some values for the predictors
##'   remaining in the reduced nomogram but adjusted values for the hiden
##'   predictors so that 0 points will be achieved from these hiden predictor
##'   variables in the full nomogram.
##' @return a list of class \code{"nomogram"} that contains information used in
##'   plotting the axes. Please see \code{\link[rms]{nomogram}} for details.
##' @author Changhong Yu, Michael Kattan, Ph.D \cr Department of Quantitative
##'   Health Sciences\cr Cleveland Clinic\cr
##' @export
##' @seealso \code{\link[rms]{nomogram}}, \code{\link[QHScrnomo]{crr.fit}},
##'   \code{\link[QHScrnomo]{pred2.crr}}, \code{\link[QHScrnomo]{nomo2.crr}}
##' @references Banks J: Nomograms. Encylopedia of Statistical Sciences, Vol 6.
##'   Editors: S Kotz and NL Johnson.  New York: Wiley; 1985.
##'
##' Lubsen J, Pool J, van der Does, E: A practical device for the application
##'   of a diagnostic or prognostic function.  Meth. Inform. Med. 17:127--129;
##'   1978.
##'
##' Wikipedia: Nomogram, \url{https://en.wikipedia.org/wiki/Nomogram}.
##'
##' Michael W. Kattan, Glenn Heller and Murray F. Brennan (2003). A
##'   competing-risks nomogram \cr for sarcoma-specific death following local
##'   recurrence. Statistics in Medicine. \code{Stat Med}. 2003;22:3515-3525.
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
##' ## make a CRR nomogram
##' nomogram.crr(prostate.crr,failtime = 120,lp=FALSE,
##' funlabel = "Predicted 10-year cumulative incidence")
##' 
##' @keywords models regression hplot
##'
nomogram.crr <-
    function(
        fit,
        failtime = NULL,
        ci = TRUE,
        ...,
        adj.to,
        lp = TRUE,
        lp.at,
        lplabel = "Linear Predictor",
        fun.at,
        fun.lp.at,
        funlabel = "Predicted Value",
        fun.side,
        interact = NULL,
        intercept = 1,
        conf.int = FALSE,
        col.conf = c(1, 12),
        conf.space = c(0.08, 0.2),
        conf.lp = c(
            "representative",
            "all", "none"
        ),
        est.all = TRUE,
        abbrev = FALSE,
        minlength = 4,
        maxscale = 100,
        nint = 10,
        label.every = 1,
        force.label = FALSE,
        xfrac = 0.35,
        cex.axis = 0.85,
        cex.var = 1,
        col.grid = FALSE,
        vnames = c("labels", "names"),
        varname.label = TRUE,
        varname.label.sep = "=",
        ia.space = 0.7,
        tck = -0.009,
        lmgp = 0.4,
        omit = NULL,
        naxes,
        points.label = "Points",
        total.points.label = "Total Points",
        total.sep.page = FALSE,
        total.fun,
        verbose = FALSE,
        total.min,
        total.max,
        mikeomit = NULL) {
        nt <- length(failtime)
        if (nt == 0) {
            stop("Specify the expected failure time!!")
        }
        
        cph.f <- fit$cph.f
        #cph.f <- ifelse(is.null(cph.f), 1, cph.f)
        #fit <- ifelse(is.null(fit), 1, fit)
        #assign("cph.f", cph.f, 1)
        #assign("fit", fit, 1)
        nm <- names(cph.f$coefficients)
        cph.f$coefficients <- fit$coef
        names(cph.f$coefficients) <- nm
        func <- mapply(list, rep(NA, nt))
        if (ci) {
            for (i in seq_len(nt)) {
                expr <- parse(
                    text = paste(
                        "function(x) nomo2.crr(x + cph.f$center, fit,time = ",
                        failtime[i],
                        ")",
                        sep = ""
                    )
                )
                func[[i]] <- eval(expr)
            }
        } else {
            for (i in seq_len(nt)) {
                expr <-
                    parse(
                        text = paste(
                            "function(x) 1 - ", 
                            "nomo2.crr(x + cph.f$center, fit,time = ",
                            failtime[i],
                            ")",
                            sep = ""
                        )
                    )
                func[[i]] <- eval(expr)
            }
        }
        nomogram.mk6(
            cph.f,
            ...,
            adj.to = adj.to,
            lp = lp,
            lp.at = lp.at,
            lplabel = lplabel,
            fun = func,
            fun.at = fun.at,
            fun.lp.at = fun.lp.at,
            funlabel = funlabel,
            fun.side = fun.side,
            interact = interact,
            intercept = intercept,
            conf.int = conf.int,
            col.conf = col.conf,
            conf.space = conf.space,
            conf.lp = conf.lp,
            est.all = est.all,
            abbrev = abbrev,
            minlength = minlength,
            maxscale = maxscale,
            nint = nint,
            label.every = label.every,
            force.label = force.label,
            xfrac = xfrac,
            cex.axis = cex.axis,
            cex.var = cex.var,
            col.grid = col.grid,
            vnames = vnames,
            varname.label = varname.label,
            varname.label.sep = varname.label.sep,
            ia.space = ia.space,
            tck = tck,
            lmgp = lmgp,
            naxes = naxes,
            points.label = points.label,
            total.points.label = total.points.label,
            total.sep.page = total.sep.page,
            total.fun = total.fun,
            verbose = verbose,
            total.min = total.min,
            total.max = total.max,
            mikeomit = mikeomit
        )
    }
