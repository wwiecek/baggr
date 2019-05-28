#' baggr - a package for Bayesian meta-analysis
#'
#' This is _baggr_ (pronounced as _bagger_ or _badger_), a Bayesian meta-analysis package for R using [Stan](https://mc-stan.org/). _Baggr_ is intended to be user-friendly and transparent so that it's easier to understand the models you are building and criticise them. The current version is a stable prototype of a tool that's in active development so we are counting on your feedback.
#'
#' _Baggr_ provides a suite of Bayesian aggregation models for both summary statistics and full data sets to synthesise evidence collected from different groups, contexts or time periods.
#'
#' @section Getting help:
#'
#' This is only a sumple package help file.
#' For documentation of the main function for conducting analyses see [baggr].
#' For description of models, data types and priors available in the package,
#' try the built-in [vignette](baggr.html).
#'
#' @docType package
#' @name baggr_package
#' @useDynLib baggr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#' @author Witold Wiecek, Rachael Meager

NULL
