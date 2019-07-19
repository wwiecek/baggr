
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *baggr*

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/wwiecek/baggr.svg?branch=master)](https://travis-ci.org/wwiecek/baggr)
<!-- badges: end -->

This is *baggr* (pronounced as *bagger*, not *badger*), a Bayesian
meta-analysis package for R using [Stan](https://mc-stan.org/). *Baggr*
is intended to be user-friendly and transparent so that it’s easier to
understand the models you are building and criticise them.

The current version (v0.1, June 2019) is a stable prototype of a tool
that’s in active development so we are counting on your feedback.

*Baggr* provides a suite of Bayesian aggregation models that work with
both summary data and full data sets, to synthesise evidence collected
from different groups, contexts or time periods.

The `baggr()` command automatically detects the data type and, by
default, fits Rubin’s (1981) partial pooling model with weakly
informative priors by calling Stan under the hood to carry out Bayesian
inference. Modelling of variances or quantiles, standardisation and
transformation of data is also possible.

## Installation

Before starting, `baggr` will not work if you don’t have RStan. In that
case, please follow [the installation instructions for
RStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

The easiest way to install `baggr` is directly from GitHub; you need to
have `devtools` package installed in R.

``` r
# compilation of models should take 5-10 minutes
devtools::install_github("wwiecek/baggr", 
                         build_vignettes = TRUE,
                         build_opts = c("--no-resave-data", "--no-manual"))
```

## Basic use case

`baggr` is designed to work well with both individual-level (“full”) and
aggregate/summary (“group”) data on treatment effect. In basic cases
only the summary information on treatment effects (such as means and
their standard errors) is needed. Data are always specified in a single
data frame and the same `baggr()` function is used for all models.

For the “standard” cases of modelling means, the appropriate model is
detected from the shape of data.

``` r
library(baggr)
df_pooled <- data.frame("tau" = c(28,8,-3,7,-1,1,18,12),
                        "se"  = c(15,10,16,11,9,11,10,18))
bgfit <- baggr(df_pooled, pooling = "partial")
```

The user has access to the underlying `stanfit` object and full joint
posterior distribution at all times. The user can also specify the model
type from several choices, the pooling type (`"none"`, `"partial"` or
`"full"`), and certain aspects of the priors.

Flexible plotting methods are included, together with automatic
comparison of multiple models (e.g. comparing no, partial and full
pooling) through `baggr_compare()` command. Various statistics can be
calculated: in particular, the `loocv()` command automatically performs
leave-one-group-out cross-validation, allowing us to compare and select
models.

Try `vignette('baggr')` for overview of these functions example of
meta-analysis workflow with `baggr`.

## Features and use cases

Current list of main features:

  - Hierarchical models for continuous outcomes
  - Meta-analysis specfic summaries and plots
  - Compatibility with `rstan` and `bayesplot` features
  - Both full and aggregate data sets can be used
  - Modelling of quantiles and SE’s
  - Modelling of log-normal data
  - Automatic standardisation of variables
  - Automatic choice of priors
  - Automatic calculation of pooling metrics
  - Leave-one-out cross-validation
