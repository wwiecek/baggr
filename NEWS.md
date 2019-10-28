# baggr 0.2.0

### New features

* Users can now define priors in `baggr()` using a syntax similar to `rstanarm`.
  Extra priors are available
* `baggr()` outputs prior predictive distributions; they can be examined using
  `baggr_compare` and `effect_plot`, `effect_draw` -- 2 new functions
* More types of model comparisons are possible
* LOO CV function has been reworked
* Full pooling and no pooling models have been reworked to avoid divergent 
  transitions.

# baggr 0.1.0

First package version for CRAN. 
