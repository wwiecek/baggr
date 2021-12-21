# baggr 0.6.10-0.6.14 (Sept-Dec 2021)

* You can add numerical values to `plot.baggr_compare` and `baggr_plot` graphics
  (a la forest plot)
* You don't need to convert summary data to individual-level data before running
  `model="logit"`, call to `baggr()` should detect it automatically now
* `pooling()` includes extra metrics, including study weights calculation
  (and better documentation)

Misc:

* You can plot only hyperparameter values (without group-specific estimates) in
  baggr_compare now
* Removed an unnecessary dependency on quantreg
* Rare event corrections (`prepare_ma()`) can now be applied either to particular 
  studies or all data
* Clearer prompts about priors and pooling in control arms when working with
  individual-level data models.

Bug fixes:

* LOO CV printing errors fixed
* LOO CV with full pooling and binary outcomes now works again after being broken in 0.6
* Individual-level Rubin  model with covariates was broken in recent 0.6 
  releases, now fixed 
* Fixed a calculation of default beta prior
* No more confusing warnings about setting `prior_control` for `"logit"` model.
* `binary_to_individual` with non-integer number of events warns user and throws 
  an error now
* Confusing results in `baggr_binary` vignette (rare events section)
* Fixes crashes for elpd calculations with unusual binary input data



# baggr 0.6.5-0.6.9 (June-August 2021)

* Mu & tau models now also print correlations between effects, via a new
  function `mutau_cor`
* You can now change type of visual comparison (`baggr_compare`) on the fly 
  (between `"effects"` and `"groups"`). Printing comparisons also returns posterior
  predictive draws.
* Upgraded forest plots to work with `forestplot` 2.0

Minor bug fixes:

* Fixed errors that could happen when using multiple factor covariates, or
  various covariate models with `loocv()`
* Fixed a bug with reporting wrong SD's for effect in the v0.6 `mutau` model
  when using `plot.baggr_compare`
* Fixed ordering of groups in `baggr_compare()`
* Various small changes to reduce amount of persistent messages 
  triggered by normal user behaviour.
* Fixed a bug where priors for meta-regressions were set even though there were
  no covariates.



# baggr 0.6.3-0.6.4 (May 2021)

* Various documentation fixes for re-submission of v0.6 to CRAN 
  (first one since v0.4).
* Added `summary` option for `effect_draw`.
* Factor covariates will work (better) now.
* Removed some non-essential code for faster compilation on CRAN.


# baggr 0.6.2 (April 2021)

* New `"mutau_full"` model is a generalisation of the `"mutau"` model into individual-level data.
  The idea is similar as for the recent `"rubin_full"` changes, see version 0.6.0.
  
* I also reparameterised the `mutau` model. It should be faster and have fewer divergent 
  transition warnings.Some of the code around the mu and tau model has also been 
  rewritten on the back end.


On the back end the package now follows the rstantools recommended way of compiling models.
The user experience should be exactly the same, but this may avoid some problems 
when installing the package from GitHub or otherwise compiling it locally.


# baggr 0.6.0 (February 2021)

### New features

* Spike and slab model can be called via `model="sslab"`. See `?baggr` for basics of
  working with this type of a model. A vignette will be added soon.
* Rubin model with full data is now called via `model="rubin_full"` rather than `"full"`. 
  Old syntax will still work, however. Made some documentation and code improvements
  around this issue.
* Leave-one-out cross-validation works for `model="rubin_full"` now. It works the same 
  way as for `model="logit"`. See `?baggr` for more information on how to use it.
* It's now possible to use `model="rubin"` with the same inputs as `model="mutau"`.
  Some data columns are removed automatically in that case.

For v0.6 we added more generic code around plotting, printing, grabbing treatment effects
etc. While there are no differences on the front-end, this means that for the next
versions we will be able to consider some new models and have more homogeneous syntax
for all models.

### Bugs

* Fixed a few issues with formatting data for individual-level data models.
* Fixed a major bug with distributions of baselines in the `rubin_full` (`full`) model.
* Fixed glitchy display for some `baggr_compare` plots.


# baggr 0.5.0 (June 2020)

### New features

* Fixed and random effects for `baggr` models now have their own separate functions,
  `fixed_effects` and `random_effects`, in addition to `group_effects`
* LOO CV works for the logistic model (as does general cross-validation).
* Vignette for binary data analysis has been rewritten in parts.
* L'Abbe plots for binary data, see `labbe()`.
* There is now more automatic conversion between summary-level and individual-level
  data for binary data (e.g. you can run `baggr()` with summary data and `model="logit"`
  for automatic conversion)
* For logistic model, priors can be specified for rates of events in the control arm,
  see arguments `prior_control` and `prior_control_sd` in `baggr()`
* There are experimental features for working with models of quantiles. 
  We advise against fitting such models using the package until these features
  have been fully tested and documented.
  

### Bug fixes

* Fixed some issues with printing of coefficients in meta-regressions,
  where wrong values were given for some models.




# baggr 0.4.0 (February 2020)

### New features

* Covariates can now be used in all baggr() models: in "rubin" model they give meta-regression
  (group-level covariates), while in "full" and "logit" models they can be used for "regular"
  regression (individual-level covariates)
* Priors for covariates are set through the argument prior_beta
* You can work with regression coefficients for covariates 
    + you can access and summarise coefficients through fixed_effects(),
    + you will also see them when printing baggr objects; 
    + when using forest_plot() you can request `show = "covariates"`
* Prototype of pp_check() function now works for Rubin model (thanks to Brice Green)
  you can apply it to generate new datasets according to posterior distribution of treatment effect
  and contrast them with the observed quantities as part of model checking
* baggr_compare() function now has standard output which you can print() or plot(), 
  thanks to Brice Green
* Vignettes and documentation were updated to better describe binary data analysis
* We now give more warnings when plugging in stupid inputs.
  
### Bug fixes

* Messages for setting priors were accidentally given when e.g. running full pooling models
* All models were re-written to standardise our approach and syntax. 
  + "Full" model might now behave differently.
  + "Mutau" model will be re-written and generalised for next release.
  + LOO CV is also disabled for some models. Prompts will be given.



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
