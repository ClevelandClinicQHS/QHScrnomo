predictDesign <-
  function(
    fit, newdata = NULL,
    type = c(
      "lp", "x", "data.frame", "terms", "adjto", "adjto.data.frame",
      "model.frame"
    ),
    se.fit = FALSE, conf.int = FALSE,
    conf.type = c("mean", "individual"),
    incl.non.slopes = NULL, non.slopes = NULL, kint = 1,
    na.action = Hmisc::na.keep, expand.na = TRUE, center.terms = TRUE, ...) {
    type <- match.arg(type)
    conf.type <- match.arg(conf.type)
    ## R does not preserve missing(x):   31jul02
    mnon.slopes <- missing(non.slopes) || !length(non.slopes)
    # was missing( ) 6jan04


    at <- fit$Design
    if (!length(at)) at <- getOldDesign(fit)
    assume <- at$assume.code
    Limval <- Getlim(at, allow.null = TRUE, need.all = FALSE)
    Values <- Limval$values
    non.ia <- assume != 9
    non.strat <- assume != 8
    f <- sum(non.ia)
    nstrata <- sum(assume == 8)
    somex <- any(non.strat)
    rnam <- NULL
    cox <- inherits(fit, "cph") ||
      (length(fit$fitFunction) && any(fit$fitFunction == "cph"))
    ## 14Nov00 22May01
    naa <- fit$na.action
    if (!expand.na) naresid <- function(a, b) b
    # don't really call naresid if drop NAs

    parms <- at$parms
    name <- at$name
    coeff <- fit$coefficients
    nrp <- Hmisc::num.intercepts(fit)

    if (mnon.slopes) {
      non.slopes <- rep(0, nrp)
      non.slopes[kint] <- 1 # 13Sep94
    } else if (nrp > 0 & length(non.slopes) != nrp) {
      stop("length of non.slopes incompatible with fit")
    }

    int.pres <- nrp > 0 # was !(cox|lrm)
    if (somex) cov <- Varcov(fit, regcoef.only = TRUE) # remove scale params
    # if(missing(incl.non.slopes))   6jan04
    if (missing(incl.non.slopes) || !length(incl.non.slopes)) {
      incl.non.slopes <- !mnon.slopes | (!missing(kint)) |
        int.pres | type != "x"
    }
    ## added 12Feb93   !missing() added 18Feb93, 2nd one 13Sep94
    int.pres <- int.pres & incl.non.slopes

    assign <- fit$assign

    nama <- names(assign)[1]
    asso <- 1 * (nama == "Intercept" | nama == "(Intercept)")

    Center <- if (cox) fit$center else 0

    oldopts <- options(
      contrasts = c(factor = "contr.treatment", ordered = "contr.poly"),
      Design.attr = at
    )

    ## 20Nov00   In SV4 options(two lists) causes problems
    on.exit({
      options(contrasts = oldopts$contrasts)
      options(Design.attr = NULL)
    })

    ## Terms <- delete.response(terms(attr(fit$terms,"formula"),
    ## specials="strat"))

    Terms <- if (.R.) {
      stats::delete.response(
        stats::terms(stats::formula(fit), specials = "strat"))
    } else {
      stats::delete.response(stats::terms(fit$terms, specials = "strat"))
    } ## 17Apr02  30may02
    attr(Terms, "response") <- 0
    attr(Terms, "intercept") <- 1 # was int.pres 12Feb93
    ## Need intercept whenever design matrix is generated to get
    ## current list of dummy variables for factor variables

    stra <- attr(Terms, "specials")$strat

    offset <- 0 # used if no newdata

    Terms.ns <- Terms
    if (length(stra)) {
      Terms.ns <- Terms[-stra]
      # Uses [.terms function. 3.0 did not add 1!
      ## was [1-stra], changed 7June94

      ## For some reason attr(...) <- pmin(attr(...)) changed a detail
      ## in factors attribute in R but R and SV4 don't seem to need this
      ## anyway   1may02
      if (!.R.) {
        tfac <- attr(Terms.ns, "factors")
        if (length(tfac) && any(tfac > 1)) {
          attr(Terms.ns, "factors") <- pmin(tfac, 1)
        }
      }
    }

    if (conf.int) {
      vconstant <- 0
      if (conf.type == "individual") {
        vconstant <- fit$stats["Sigma"]^2
        if (is.na(vconstant)) {
          stop(
            'conf.type="individual" requires",
                        "that fit be from ols')
        }
      }
      zcrit <- if (length(idf <- fit$df.residual)) {
        stats::qt((1 + conf.int) / 2, idf)
      } else {
        stats::qnorm((1 + conf.int) / 2)
      }
    }

    if (type != "adjto" & type != "adjto.data.frame") {
      X <- NULL
      if (missing(newdata) || !length(newdata)) {
        if (type == "lp" && length(fit$linear.predictors)) {
          LP <- naresid(naa, fit$linear.predictors)
          # changed 8June94
          if (kint > 1) LP <- LP - fit$coef[1] + fit$coef[kint]
          # added 13Sep94
          if (length(stra<-attr(fit$linear.predictors, "strata"))) {
            attr(LP, "strata") <- naresid(naa, stra)
          }
          if (!se.fit && !conf.int) {
            return(LP)
          } ## 7Mar99
          else if (length(fit$se.fit)) {
            if (kint > 1) {
              stop(
                "se.fit is retrieved from the fit",
                "but it corresponded to kint=1"
              )
            }
            retlist <- list(
              linear.predictors = LP,
              se.fit = naresid(naa, fit$se.fit)
            )
            if (conf.int) {
              plminus <- zcrit * sqrt(
                retlist$se.fit^2 + vconstant
              )
              retlist$lower <- LP - plminus
              retlist$upper <- LP + plminus
            }
            return(retlist)
          }
        } else if (type == "x") {
          return(structure(
            naresid(naa, fit$x),
            strata = if (length(stra <- attr(fit$x, "strata"))) {
              naresid(naa, stra)
            } else {
              NULL
            }
          ))
        }
        X <- fit$x
        rnam <- dimnames(X)[[1]]
        if (!any(names(fit) == "x")) X <- NULL # fit$x can get fit$xref
        if (!length(X)) {
          stop("newdata not given and fit did not store x")
        }
      }
      if (!length(X)) {
        if (!is.data.frame(newdata)) {
          if (is.list(newdata)) {
            loc <- if (!length(names(newdata))) {
              seq(1, f)
            } else {
              name[assume != 9]
            }
            new <- matrix(
              if (.R.) double(1) else single(1),
              nrow = length(newdata[[1]]),
              ncol = length(newdata)
            )
            for (j in seq(1, ncol(new))) {
              new[, j] <- newdata[[loc[j]]]
            }
            newdata <- new
          }
          if (!is.matrix(newdata)) {
            newdata <- matrix(newdata, ncol = f)
          }
          if (ncol(newdata) != f) {
            stop("# columns in newdata not= # factors in design")
          }
          X <- list()
          k <- 0
          ii <- 0
          for (i in (seq(1, length(assume)))[non.ia]) {
            ii <- ii + 1
            xi <- newdata[, ii]
            as <- assume[i]
            allna <- all(is.na(xi))
            ## if(as!=10 && allna) xi <- at$limits[3,ii]
            if (as == 5 | as == 8) {
              xi <- as.integer(xi)
              levels(xi) <- parms[[name[i]]]
              oldClass(xi) <- "factor"
            } else if (as == 10) {
              if (i == 1) {
                ifact <- 1
              } else {
                ifact <- 1 + sum(assume[seq(1, (i - 1))] != 8)
              }
              ## Accounts for assign not being output
              ## for strata factors
              ncols <- length(assign[[ifact + asso]])
              if (allna) {
                xi <- matrix(
                  if (.R.) double(1) else single(1),
                  nrow = length(xi), ncol = ncols
                )
                for (j in seq(1, ncol(xi))) {
                  xi[, j] <- parms[[name[i]]][j]
                }
              } else {
                xi <- matrix(
                  if (.R.) xi else as.single(xi),
                  nrow = length(xi), ncol = ncols
                )
              }
            }
            ## Duplicate single value for all parts of matrix
            k <- k + 1
            X[[k]] <- xi
          }
          names(X) <- name[non.ia]
          attr(X, "row.names") <- as.character(seq(1, nrow(newdata)))
          oldClass(X) <- "data.frame"
          newdata <- X
          ## Note: data.frame() converts matrix variables to
          ## individual variables
          if (type == "data.frame") {
            return(newdata)
          }
        } else {
          ## Need to convert any factors to have all levels in
          ## original fit Otherwise, wrong dummy variables will
          ## be generated by model.matrix
          nm <- names(newdata)
          for (i in seq(1, ncol(newdata))) {
            j <- match(nm[i], name)
            if (!is.na(j)) {
              asj <- assume[j]
              w <- newdata[, i]
              V <- NULL
              if (asj == 5 | asj == 7 | asj == 8 |
                  (name[j] %in% names(Values) &&
                   length(V <- Values[[name[j]]]) &&
                   is.character(V))) {
                if (length(Pa <- parms[[name[j]]])) V <- Pa
                # added 8Apr94
                ## if(is.null(V)) V <- parms[[name[j]]]
                # subtracted 8Apr94
                newdata[, i] <- factor(w, V)
                ## Handles user specifying numeric values
                ## without quotes, that are levels
                ww <- is.na(newdata[, i])&!is.na(unclass(w))
                if (any(ww)) {
                  cat(
                    "Error in predictDesign: Values in",
                    names(newdata)[i],
                    "not in", V, ":\n"
                  )
                  print(as.character(w[ww]), quote = FALSE)
                  stop()
                }
              }
            }
          }
        }

        newdata <- addOffset4ModelFrame(Terms, newdata) ## 23nov02
        X <- stats::model.frame(
          Terms,
          newdata, na.action = na.action, ...)
        if (type == "model.frame") {
          return(X)
        }
        naa <- attr(X, "na.action")
        rnam <- row.names(X)

        offs <- attr(Terms, "offset")
        if (!length(offs)) {
          offset <- rep(0, length(rnam))
        } else {
          offset <- X[[offs]]
        }

        ## if(ncol(X) != sum(non.ia))
        ##     stop("improperly formed model frame")
        strata <- list()
        nst <- 0
        ii <- 0 ## 23nov02
        for (i in setdiff(seq(1, ncol(X)), offs)) {
          ## setdiff() was 1:ncol(X) 23nov02
          ii <- ii + 1
          xi <- X[[i]]
          asi <- attr(xi, "assume.code")
          as <- assume[ii] ## was i 23nov02
          if (!length(asi) && as == 7) {
            attr(X[, i], "contrasts") <-
              attr(scored(xi, name = name[ii]), "contrasts")
            ## was i 23nov02
            if (length(xi) == 1) {
              stop(
                "a bug in model.matrix can produce",
                "incorrect results\n",
                "when only one observation is being",
                "predicted for an ordered variable"
              )
            }
          }

          if (as == 8) {
            nst <- nst + 1
            strata[[nst]] <- paste(
              name[ii], "=",
              parms[[name[ii]]][X[, i]],
              sep = ""
            )
            ## was name[i] 23nov02
          }
        }
        if (!somex) {
          X <- NULL
        } else if (int.pres && nrp == 1) {
          X <- stats::model.matrix(Terms.ns, X)
        } # nrp Jan94
        else {
          X <- stats::model.matrix(Terms.ns, X)[, -1, drop = FALSE]
        } # 12Feb93
        if (nstrata > 0) {
          names(strata) <- paste("S", seq(1, nstrata), sep = "")
          strata <- factor(interaction(
            as.data.frame(strata),
            drop = TRUE
          ),
          levels = fit$strata
          )
        }
      } else {
        strata <- attr(X, "strata")
      }

      added.col <- if (
        incl.non.slopes & (nrp > 1 | !int.pres)) nrp else 0
      # nrp>1 Jan94
      ## & !scale.pres removed from following statement 20Feb93
      if (incl.non.slopes & nrp > 0 & somex & added.col > 0) {
        xx <- matrix(
          if (.R.) double(1) else single(1),
          nrow = nrow(X), ncol = added.col
        )
        for (j in seq(1, nrp)) xx[, j] <- non.slopes[j]
      } else {
        xx <- NULL
      }
    }

    ## For models with multiple intercepts, delete elements of
    ## covariance matrix containing unused intercepts
    elements.to.delete <- 9999
    if (somex && nrp > 1) {
      i <- seq(1, nrp)[non.slopes == 0]
      cov <- cov[-i, -i, drop = FALSE]
      elements.to.delete <- i
    }

    if (type == "adjto" | type == "adjto.data.frame" |
        (center.terms && type == "terms") |
        (cox & (se.fit | conf.int))) {
      ## Form design matrix for adjust-to values
      adjto <- list()
      ii <- 0
      for (i in seq(1, length(assume))[non.ia]) {
        ii <- ii + 1
        xi <- Getlimi(name[i], Limval, need.all = TRUE)[2]
        # was =F  5Feb94
        if (assume[i] == 5 | assume[i] == 8) {
          xi <- factor(xi, parms[[name[i]]])
        } else if (assume[i] == 7) {
          xi <- scored(xi, name = name[i])
        } else if (assume[i] == 10) {
          xi <- matrix(parms[[name[i]]], nrow = 1)
        }
        # matrx col medians
        adjto[[ii]] <- xi
      }
      names(adjto) <- name[non.ia]
      ##   adjto <- data.frame(adjto,check.names=FALSE)
      ##   data.frame will take matrix factors and convert into
      ##   individual vars
      attr(adjto, "row.names") <- "1"
      oldClass(adjto) <- "data.frame"
      if (type == "adjto.data.frame") {
        return(adjto)
      }
      adjto <- addOffset4ModelFrame(Terms, adjto) ## 23nov02
      adjto <- stats::model.frame(Terms, adjto)
      adjto <- if (int.pres) {
        stats::model.matrix(Terms.ns, adjto)
      } else {
        stats::model.matrix(Terms.ns, adjto)[, -1, drop = FALSE]
      } # -1 added 12Feb93
      ## added drop=FALSE 27feb03
      if (type == "adjto") {
        k <- (if (int.pres) {
          seq(1, length(coeff))
        } else {
          seq(nrp + 1, length(coeff))
        })
        if (is.matrix(adjto)) {
          dimnames(adjto) <- list(
            dimnames(adjto)[[1]],
            names(coeff)[k]
          )
        } else {
          names(adjto) <- names(coeff)[k]
        }
        return(adjto)
      }
    }

    if (length(xx) & type != "terms" & incl.non.slopes) {
      X <- cbind(xx, X)
      dimnames(X) <- list(rnam, names(coeff))
      if ((center.terms && type == "terms") |
          (cox & (se.fit | conf.int))) {
        adjto <- c(xx[1, ], adjto)
      } # 12Feb93
    } else if (somex) {
      dimnames(X) <-
        ## list(rnam,names(coeff)[
        ## (nrp+1-(int.pres & incl.non.slopes)):length(coeff)])
        list(rnam, names(coeff)[
          seq(1 + length(coeff) - ncol(X), length(coeff))
        ])
    } # 22Jun95


    storage.mode(X) <- "double"
    if (type == "x") {
      return(
        structure(
          naresid(naa, X),
          strata = if (nstrata > 0) naresid(naa, strata) else NULL,
          offset = if (length(offs)) naresid(naa, offset) else NULL,
          na.action = if (expand.na) NULL else naa
        )
      )
    }

    if (type == "lp") {
      if (somex) {
        ## if( ) 28apr02
        if (any(elements.to.delete == 9999)) {
          cof <- coeff
        } else {
          cof <- coeff[-elements.to.delete]
          X <- X[, -elements.to.delete, drop = FALSE]
        }
        xb <- Hmisc::matxv(X, cof) + offset - Center
        names(xb) <- rnam
        if (!.R.) storage.mode(xb) <- "single"
      } else {
        if (!length(offs)) xb <- NULL else xb <- offset
      }
      xb <- naresid(naa, xb)
      if (nstrata > 0) attr(xb, "strata") <- naresid(naa, strata)
      if ((se.fit | conf.int) & somex) {
        if (cox) X <- sweep(X, 2, adjto) # Center columns
        se <- drop(sqrt(((X %*% cov) * X) %*% rep(1, ncol(X))))
        names(se) <- rnam
        if (!.R.) storage.mode(se) <- "single"
        retlist <- structure(list(
          linear.predictors = xb,
          se.fit = naresid(naa, se)
        ),
        na.action = if (expand.na) NULL else naa
        )
        if (conf.int) {
          plminus <- zcrit * sqrt(retlist$se.fit^2 + vconstant)
          retlist$lower <- xb - plminus
          retlist$upper <- xb + plminus
        }
        return(retlist)
      } else {
        return(structure(xb, na.action = if (expand.na) NULL else naa))
      }
    }

    if (type == "terms") {
      if (!somex) {
        stop(
          'type="terms" may not be given',
          "unless covariables present"
        )
      }
      fitted <- array(
        0, c(nrow(X), sum(non.strat)),
        list(rnam, name[non.strat])
      )
      if (se.fit) se <- fitted
      j <- 0
      if (center.terms) {
        ## 31jul02: lrm and perhaps others put out fit$x without
        ## column of intercepts but model has intercept
        if (ncol(adjto) != ncol(X)) {
          if (dimnames(adjto)[[2]][1] %in% c(
            "Intercept",
            "(Intercept)"
          ) &&
          !(dimnames(X)[[2]][1] %in% c(
            "Intercept",
            "(Intercept)"
          ))) {
            adjto <- adjto[, -1, drop = FALSE]
          }
        }
        # if (ncol(adjto) != ncol(X)) stop("program logic error")
        X <- sweep(X, 2, adjto) # center columns
      }
      # PROBLEM: adjto = c(Intercept=1, sexmale=0); no 1s col in f$x
      num.intercepts.not.in.X <- length(coeff) - ncol(X) # 23Jan95
      for (i in seq(1, length(assume))[non.strat]) {
        j <- j + 1
        k <- assign[[j + asso]] # ; m <- k+int.pres
        ko <- k - num.intercepts.not.in.X # 23Jun95
        fitted[, j] <- Hmisc::matxv(X[, ko, drop = FALSE], coeff[k])
        ## was X[,m], coeff[nrp+k]
        if (se.fit) {
          se[, j] <- (((
            X[, ko, drop = FALSE] %*%
              cov[ko, ko, drop = FALSE]) *
              X[, ko, drop = FALSE]) %*% rep(1, length(ko)))^.5
        }
      }
      if (!.R.) storage.mode(fitted) <- "single"
      fitted <- structure(
        naresid(naa, fitted),
        strata = if (nstrata == 0) {
          NULL
        } else {
          naresid(naa, strata)
        }
      )
      if (se.fit) {
        if (!.R.) storage.mode(se) <- "single"
        return(structure(list(
          fitted = fitted,
          se.fit = naresid(naa, se)
        ),
        na.action = if (expand.na) NULL else naa
        ))
      } else {
        return(structure(
          fitted,
          na.action = if (expand.na) NULL else naa
        ))
      }
    }
  }

addOffset4ModelFrame <- function(Terms, newdata, offset = 0) {
  offs <- attr(Terms, "offset")
  if (!length(offs)) {
    return(newdata)
  }
  ##  offsetVarname <- all.names(attr(Terms,'variables')[offs+1])[1] 12mar04
  offsetVarname <- setdiff(
    all.names(attr(Terms, "variables")[offs + 1]),
    "offset"
  )
  if (length(offsetVarname) > 1) {
    stop(
      "More than one offset variable, only first used:",
      offsetVarname
    )
    offsetVarname <- offsetVarname[1]
  }
  ##  offsetVarname <- offsetVarname[offsetVarname != 'offset']
  if (!(offsetVarname %in% names(newdata))) {
    newdata[[offsetVarname]] <- rep(offset, length = nrow(newdata))
    stop(
      "offset variable set to",
      format(offset)
    )
  }
  newdata
}
