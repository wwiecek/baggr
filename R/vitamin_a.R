#' Vitamin A supplementation and child mortality
#'
#' Summary-level data from the Imdad et al. (2022) meta-analysis of vitamin A
#' supplementation and all-cause child mortality in children aged six months to
#' five years. Effects are log risk ratios.
#'
#' The source CSV contains one row for Lin 2008 with `tau = 0` and `se = 0`.
#' That row is excluded here because the object is intended to be used directly
#' with [baggr()] and a zero standard error is not a valid summary-data input.
#'
#' @format A data frame with 18 rows and 3 columns:
#' \describe{
#'   \item{group}{Study label.}
#'   \item{tau}{Estimated log risk ratio for all-cause child mortality.}
#'   \item{se}{Standard error of the estimated log risk ratio.}
#' }
#' @references Imdad, A., Mayo-Wilson, E., Herzer, K., & Bhutta, Z. A. (2022).
#' Vitamin A supplementation for preventing morbidity and mortality in children
#' from six months to five years of age. Cochrane Database of Systematic
#' Reviews, 2022(3), CD008524.
"vitamin_a"
