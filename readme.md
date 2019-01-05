
<!-- README.md is generated from README.Rmd. Please edit that file -->
baggr
=====

Overview
--------

`baggr` (pronounced as *bagger* or *badger*) is a Bayesian meta-analysis package for R using Stan. `baggr` is intended to be user-friendly and transparent so that it's easier to understand the models you are building and criticise them. The current version (Dec 2018) is a stable prototype of a tool that's in active development.

To install `baggr` from github source use the `devtools` package (installation requires `rstan`, `Rcpp` and RTools):

``` r
devtools::install_github("wwiecek/baggr")
```

Basic use case
--------------

`baggr` is designed to work well with both individual-level ("full") and aggregate ("group") data on treatment effect. In basic cases only summary information on treatment effects (such as means and their standard errors) is needed. Data are always specified in a single data frame and the same `baggr()` function is used for all models. For "standard" cases of modelling means the appropriate model is detected from the shape of data:

``` r
library(baggr)
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3), 
                        "se" = rep(1, 8))
bgfit <- baggr(df_pooled, pooling = "full")
```

Flexible plotting methods are included, together with automatic comparison of multiple models (e.g. no, partial and full pooling).

Features and use cases
----------------------

Current list of main features:

-   Hierarchical models for continuous outcomes
-   Either full or aggregate data can be used
-   Modelling of log-normal data
-   Automatic standardisation of variables
-   Automatic choice of priors
-   Modelling of quantiles and SE's
-   Automatic calculation of pooling metrics
-   Leave-one-out cross-validation
-   Meta-analysis specfic summaries and plots

Try `vignette('baggr')` for example of meta-analysis workflow with `baggr`.
