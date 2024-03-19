# QHScrnomo (development version)

# QHScrnomo 3.0.1

* Fixed issue in `crr.fit` caused by extracting column when input data set is a `tibble` (#1)

# QHScrnomo 3.0.0

## Major changes

* Lazily-load internal package data set (i.e., `QHScrnomo::prostate.dat`) instead of using the `data()` function
* Changed dependency structure to only use the minimal necessary (see DESCRIPTION)
  + `rms` is now the only package that `QHScrnomo` "Depends" on (i.e., it will be loaded when this package is loaded)
  + `cmprsk` and `Hmisc` have been moved to the "Imports" field and are accessed internally with the `pkg::func` format. Therefore, users should _not_ expect these packages to load with `QHScrnomo` (though they are required to be installed).
* Changed the output of `sas.cmprsk` to print the base _failure_ probability (`1 - exp(-cumsum(f$bfitj)`) when the `time` argument is used instead of the _survival_ probability (`exp(-cumsum(f$bfitj)`)
* Removed the `pred.crr` function from the package. This function was not previously exported or used, and was highly duplicative of `predict.cmprsk`. Nevertheless, it is no longer in the package.
* Removed .Rd files for functions that are not exported (`addOffset4ModelFrame`, `nomo2.crr`, `pred2.crr`, `pred3.crr`, `predictDesign`)
* Changed the default argument for `failcode` in `pred.ci` to `failcode=1` to match `cmprsk::crr` and `QHScrnomo::crr.fit`
* Set the default argument for `time` in `tenf.crr` and `predict.crr` to `time=NULL`; also set default arguments for `g`, `cuts`, `xlab`, `ylab` arguments in `groupci`
* Added the `trace` argument to `tenf.crr` to _optionally_ trace the process in the console. It still does it by default, but previously it always traced it with no option to turn it off.
* In `groupci`, the `ci` argument allows the user to control whether we work on the scale of _survival_ or _failure_ (default) probabilities. When the argument was set to `FALSE`, this would apply a `1-value` operation to _both_ the within-group cumulative incidence estimate (to convert to survival) _and_ the group-level mean (of `x`). However, `x` is supposed to be an arbitrary continuous variable that we can assess calibration for, not necessarily a probability. Therefore, the change made in this version is to _only_ apply the `1-value` transformation to the within-group cumulative incidence and _not_ the group level mean for `x`. 
  + To retain the same behavior when a probability is entered for `x`, simply pass `1-x` into the function call. 
  + Users can also use the **new** `a` and `b` arguments (for the intercept and slope, respectively), as well as `ab`, `xlim`, and `ylim` arguments to adjust the graph accordingly. 

## Minor improvements or bug fixes

* Added more comprehensive error-handling throughout with more informative messaging:
  + Checking for the appropriate input objects in (almost) all functions
  + Checking for the appropriate time point inputs in `predict.cmprsk`, `sas.cmprsk`, `pred.ci`, `tenf.crr`, `groupci`, `nomogram.crr`
  + Other miscellaneous warnings and/or errors (see unit tests)
* Added a unit testing suite via `testthat`
* Changed, edited, and/or rewrote most (if not all) function documentation text (.Rd files) to be more descriptive and provide more context
* Removed definition and usage of the `oldUnclass` utility (internal) function. It's definition was `function(x) unclass(x)`, so we just use `unclass(x)` instead.
* In `groupci`, the user is now able to supply `xlab = ""` as an argument. Also removed the usage of `single` in the function definition in favor of `double` as the former is only to be used in the context of `.C/.Fortran` (per documentation), which is not the case here.
* In `nomogram.crr`, changed some arguments:
  + Removed the `NULL` default for `failtime` to make it more clear that this is a required argument
* Majority of internal/utility functions are used in `nomogram.crr` and `nomogram.mk6`. Moved short utility functions to be housed in `R/nomogram.R` directly (`nomo2.crr`, `is.category`, `getOldDesign`, `Design.levels`, `value.chk`, `axisf`)
  + Still maintaining separate files in `R/` for `Getlim`, `predictDesign` (contains `addOffset4ModelFrame`), `Varcov`

## Historical changes

Adding historical change-tracking here that was included in the previous version's NEWS file (2.2.0) (in reverse order), but in a Markdown-friendly format:

* 2011-05-23 - fix a bug in predict.cmprsk
* 2011-04-19 - change the random sampling scheme in tenf.crr so that the number of patients in each fold is more balanced
* 2011-04-08 - fix a bug in tenf.crr
* 2011-04-13 - add Newlabels.cmprsk and Newlevels.cmprsk
* 2011-03-01 - add cindex, summary.crr and anova.crr
* 2009-06-03 - fixed a bug in sas.crr
* 2008-12-31 - fixed a bug in crr.fit
* 2008-10-16 - edit crr.fit and tenf.crr to deal with subset options in cph
* 2008-08-04 - modify function 'pred.ci' to accomodate when failcode == 1
* 2008-04-02 - add 'tenf.crr'
* 2008-04-02 - add a check to 'crr.fit' if bfitj are very large(infinite), stop giving a error message to check if there is a predictor with large integers( i.e, year of surgery)
* 2008-02-21 - add function 'pred.crr','pred.ci' and 'groupci'
* 2008-03-20 - set version to 0.1-3; add flag 'ci' to 'nomogram.crr'
* 2008-03-17 - capitalize package name as QHScrnomo
* 2008-03-14 - add nomogram.crr to correct lp and set verion number to 0.1-2; add function zzz.R
* 2008-03-11 - start to construct the first version 0.1-1
