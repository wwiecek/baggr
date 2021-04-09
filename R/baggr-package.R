#' baggr - a package for Bayesian meta-analysis
#'
#' This is _baggr_ (pronounced as _bagger_ or _badger_), a Bayesian meta-analysis
#' package for R that uses [Stan](https://mc-stan.org/) to fit the models.
#' _Baggr_ is intended to be user-friendly and transparent so that
#' it's easier to understand the models you are building and criticise them.
#'
#' _Baggr_ package provides a suite of models that work with both summary data and full data sets,
#' to synthesise evidence collected from different groups, contexts or time periods.
#' The [baggr] command automatically detects the data type and, by default, fits a partial
#' pooling model (which you may know as
#' [random effects models](https://stats.stackexchange.com/questions/4700/what-is-the-difference-between-fixed-effect-random-effect-and-mixed-effect-mode))
#' with weakly informative priors by calling [Stan](https://mc-stan.org/) to carry
#' out Bayesian inference. Modelling of variances or quantiles, standardisation and
#' transformation of data are also possible.
#'
#'
#' @section Getting help:
#'
#' This is only a simple package help file.
#' For documentation of the main function for conducting analyses see [baggr].
#' For description of models, data types and priors available in the package,
#' try the built-in vignette (`vignette("baggr")`).
#'
#' @docType package
#' @name baggr-package
#' @useDynLib baggr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
NULL
