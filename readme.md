
<!-- README.md is generated from README.Rmd. Please edit that file -->
baggr
=====

-- This is a prototype version of a tool that's in development; we plan to release 1.0 version in October 2018 --

Overview
--------

`baggr` (pronounced as *bagger* or *badger*) is a Bayesian meta-analysis package for R using Stan. `baggr` is intended to be user-friendly and transparent so that it's easier to understand the models you are building and criticise them.

Basic use case
--------------

`baggr` is designed to work well with both individual-level ("full") and summarised ("group") data on treatment effect. In most basic case only summary information on treatment effects (such as means and their standard errors) is needed. Data are always specified in a single data frame and the same `baggr()` function is used for all models. For "standard" cases of modelling means the appropriate model is detected from shape of data:

``` r
library(baggr)
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3), 
                        "se" = rep(1, 8))
bgfit <- baggr(df_pooled, pooling = "full")
```

Flexible plotting methods are included, together with comparison of multiple models.

Realistic use cases
-------------------

Try `vignette('baggr')` for example of meta-analysis workflow with `baggr`.
