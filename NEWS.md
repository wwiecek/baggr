# baggr 0.3.0

### New features

* Binary data models for both summary and individual-level data.
* New vignette for working with binary data; see `vignette("baggr_binary")`.
* Expanded helper functions (esp. `prepare_ma`), esp. for prepping binary data.
* Added forest plot functionality for all types of models.
* Various outputs can now be transformed (main use case is `exp`, but any transform is allowed).
* Reworked vignette sections for pooling and cross-validation.
* Pooling statistics are now calculated for the whole model and better documented.
* More consistent theming, similar to bayesplot (thanks to Brice Green)
* Comparison of leave-one-out cross-validations with `loo_compare` (thanks to Brice Green)
  
### Bug fixes

* Re-enabled missing Cauchy priors


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
