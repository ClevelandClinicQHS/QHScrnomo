#' Prostate cancer data set
#'
#' This is an artificial prostate cancer dataset used for illustrating the usage of functions in \code{QHScrnomo}
#'
#' @format ## `prostate.dat`
#' A data frame with 2000 observations on the following 9 variables:
#' \describe{
#'   \item{UNIQID}{patient ID}
#'   \item{TX}{Treatment options of prostate cancera with levels \code{EBRT}, \code{PI}, \code{RP}}
#'   \item{PSA}{Pre-treatment PSA levels}
#'   \item{BX_GLSN_CAT}{Biopsy Gleason Score Sum. a factor with levels \code{1} for 2-6 \code{2} for 7 and \code{3} for 8-10}
#'   \item{CLIN_STG}{Clinical stage with levels \code{T1}, \code{T2}, \code{T3}}
#'   \item{AGE}{Age at treatment date}
#'   \item{RACE_AA}{patient ethnicity.a factor with levels \code{0} for other and \code{1} for African American}
#'   \item{TIME_EVENT}{follow up time in months}
#'   \item{EVENT_DOD}{followup status, 0 - censored, 1 - died of prostate cancer, 2 - died of other causes}
#' }
#' @details This is a simulated dataset
"prostate.dat"
