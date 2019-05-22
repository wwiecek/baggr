
<!-- README.md is generated from README.Rmd. Please edit that file -->

# baggr

`baggr` (pronounced as *bagger* or *badger*) is a Bayesian meta-analysis
package for R using [Stan](https://mc-stan.org/). `baggr` is intended to
be user-friendly and transparent so that it’s easier to understand the
models you are building and criticise them. The current version is a
stable prototype of a tool that’s in active development so we are
counting on your feedback.

## Installation

Before starting, `baggr` will not work if you don’t have RStan. In that
case please follow [the installation instructions for
RStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

To install `baggr` directly from GitHub you can use the `devtools`
package:

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
data frame and the same `baggr()` function is used for all models. For
the “standard” cases of modelling means, the appropriate model is
detected from the shape of data:

``` r
library(baggr)
df_pooled <- data.frame("tau" = c(28,8,-3,7,-1,1,18,12),
                        "se"  = c(15,10,16,11,9,11,10,18))
bgfit <- baggr(df_pooled, pooling = "full")
```

Flexible plotting methods are included, together with automatic
comparison of multiple models (e.g. comparing no, partial and full
pooling). Try `vignette('baggr')` for example of meta-analysis workflow
with `baggr`.

## Features and use cases

Current list of main features:

  - Hierarchical models for continuous outcomes
  - Meta-analysis specfic summaries and plots
  - Compatibility with RStan and `bayesplot` features
  - Both full and aggregate data sets can be used
  - Modelling of quantiles and SE’s
  - Modelling of log-normal data
  - Automatic standardisation of variables
  - Automatic choice of priors
  - Automatic calculation of pooling metrics
  - Leave-one-out cross-validation
