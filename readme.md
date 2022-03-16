# baggr: Bayesian aggregator for R, v0.6 (March 2022)


<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-last-release/baggr?color=green)](http://cran.r-project.org/package=baggr)
[![](https://cranlogs.r-pkg.org/badges/baggr)](https://cran.rstudio.com/web/packages/baggr/index.html)
[![Codecov\_test\_coverage](https://codecov.io/gh/wwiecek/baggr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/wwiecek/baggr?branch=master)
<!-- badges: end -->

This is *baggr*, an [R package](https://www.r-project.org/) for Bayesian
meta-analysis using [Stan](https://mc-stan.org/). *Baggr* is intended to
be user-friendly and transparent so that it’s easier to understand the
models you are building and criticise them.

*Baggr* provides a suite of models that work with both summary data and
full data sets, to synthesise evidence collected from different groups,
contexts or time periods. The `baggr()` command automatically detects
the data type and, by default, fits a partial pooling model (which you
may know as [random effects
models](https://stats.stackexchange.com/questions/4700/what-is-the-difference-between-fixed-effect-random-effect-and-mixed-effect-mode))
with weakly informative priors by calling [Stan](https://mc-stan.org/)
to carry out Bayesian inference. Modelling of variances or quantiles,
standardisation and transformation of data is also possible.

The current version is a stable version of a tool that’s in active
development so we are counting on your feedback.

## Installation

Before starting, please follow [the installation instructions for
RStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started),
which is responsible for Bayesian inference in *baggr*. If you don’t
have Stan, it’s worth following the instructions step-by-step.

The package itself is available on CRAN:

``` r
install.packages("baggr")
```

You can also install the most up-to-date version of `baggr` directly
from GitHub; this is what we recommend, but to do that you will need the
`remotes` package:

``` r
#installation this way may take 5-15 minutes
remotes::install_github("wwiecek/baggr", 
                        ref = "devel", #if problems try changing to ref = "master"
                        build_vignettes = TRUE, quiet = TRUE,
                        build_opts = c("--no-resave-data", "--no-manual"))
```

Most common issue in installing *baggr* is with updating other packages.
Try updating your packages (and ensure R is at least version 4) before
trying the `remotes` command.

## Basic use case

`baggr` is designed to work well with both individual-level (“full”) and
aggregate/summary (“group”) data on treatment effect. In basic cases,
only the summary information on treatment effects (such as means and
their standard errors) is needed. Data are always specified in a single
input data frame and the same `baggr()` function is used for different
models.

For the “standard” cases of modelling means, the appropriate model is
detected from the shape of data.

``` r
library(baggr)
df_pooled <- data.frame("tau" = c(28,8,-3,7,-1,1,18,12),
                        "se"  = c(15,10,16,11,9,11,10,18))
bg <- baggr(df_pooled, pooling = "partial")
```

You can specify the model type from several choices, the pooling type
(`"none"`, `"partial"` or `"full"`), and certain aspects of the priors,
as well as other options for data preparation, prediction and more. You
can access the underlying `stanfit` object through `bg$fit`.

Flexible plotting methods are included, together with an automatic
comparison of multiple models (e.g. comparing no, partial and full
pooling) through `baggr_compare()` command. Various statistics can be
calculated: in particular, `pooling()` for pooling metrics and `loocv()`
for leave-one-group-out cross-validation, allowing us to then compare
and select models via `loo_compare()`. Forest plots and plots of
treatment effects are available.

Try `vignette('baggr')` for an overview of these functions and an
example of meta-analysis workflow with `baggr`. If working with binary
data, try `vignette("baggr_binary")`. Compiled vignettes are available
[on CRAN](https://cran.r-project.org/web/packages/baggr/index.html).

## Current and future releases

Included in baggr v0.6 (2022):

  - Meta-analysis of continuous and binary outcomes
  - Both full and aggregate data sets can be used
  - Summaries and plots specific to meta-analysis, typical diagnostic
    plots
  - Meta-regression / fixed effects modelling
  - Compatibility with `rstan` and `bayesplot` features
  - Automatic choice of priors or “plain-text” specification of priors
  - Calculation of pooling/heterogeneity metrics
  - Cross-validation (leave-one-group-out)
  - Prior and posterior predictive distributions

Check [NEWS.md] for more information on recent changes to the package.
