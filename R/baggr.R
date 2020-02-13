#' Bayesian aggregate treatment effects model
#'
#' Bayesian inference on parameters of an average treatment effects model
#' that's appropriate to the supplied
#' individual- or group-level data, using Hamiltonian Monte Carlo in Stan.
#' (For overall package help file see [baggr-package])
#'
#'
#' @importFrom rstan summary
#' @importFrom rstan sampling
#'
#' @param data data frame with summary or individual level data to meta-analyse
#' @param model if \code{NULL}, detected automatically from input data
#'              otherwise choose from
#'              \code{"rubin"}, \code{"mutau"}, \code{"individual"}, \code{"quantiles"}
#'              (see Details).
#' @param pooling Type of pooling;
#'                choose from \code{"none"}, \code{"partial"} (default) and \code{"full"}.
#'                If you are not familiar with the terms, consult the vignette;
#'                "partial" can be understood as random effects and "full" as fixed effects
#' @param effect Label for effect. Will default to "mean" in most cases, "log OR" in logistic model,
#'               quantiles in `quantiles` model etc.
#'               These labels are used in various print and plot outputs.
#'               Comparable models (e.g. in [baggr_compare]) should have same `effect`.
#' @param covariates Character vector with column names in `data`. The corresponding columns are used as
#'                   covariates (fixed effects) in the meta-regression model (in case of aggregate data).
#'                   In the case of individual level data the model does not differentiate between group-level
#'                   variables (same values of the covariate for all rows related to a given group) and
#'                   individual-level covariates.
#' @param prior_hypermean prior distribution for hypermean; you can use "plain text" notation like
#'              `prior_hypermean=normal(0,100)` or `uniform(-10, 10)`.
#'              See Details:Priors below for more possible specifications.
#'              If unspecified, the priors will be derived automatically based on data
#'              (and printed out in the console).
#' @param prior_hypersd  prior for hyper-standard deviation, used
#'                       by Rubin and `"mutau"`` models;
#'                       same rules apply as for `_hypermean`;
#' @param prior_hypercor prior for hypercorrelation matrix, used by the `"mutau"` model
#' @param prior_beta prior for regression coefficients if `covariates` are specified; will default to
#'                       experimental normal(0, 10^2) distribution
#' @param prior alternative way to specify all priors as a named list with `hypermean`,
#'              `hypersd`, `hypercor`, `beta`, analogous to `prior_` arguments above,
#'              e.g. `prior = list(hypermean = normal(0,10), beta = uniform(-50, 50))`
#' @param ppd       logical; use prior predictive distribution? (_p.p.d._) Default is no.
#'                  If `ppd=TRUE`, Stan model will sample from the prior distributions
#'                  and ignore `data` in inference. However, `data` argument might still
#'                  be used to infer the correct model and to set the default priors.
#' @param outcome   character; column name in (individual-level)
#'                  \code{data} with outcome variable values
#' @param group     character; column name in \code{data} with grouping factor;
#'                  it's necessary for individual-level data, for summarised data
#'                  it will be used as labels for groups when displaying results
#' @param treatment character; column name in (individual-level) \code{data} with treatment factor;
#' @param quantiles if \code{model = "quantiles"}, a vector indicating which quantiles of data to use
#'                  (with values between 0 and 1)
#' @param test_data data for cross-validation; NULL for no validation, otherwise a data frame
#'                  with the same columns as `data` argument
#' @param warn print an additional warning if Rhat exceeds 1.05
#' @param ... extra options passed to Stan function, e.g. \code{control = list(adapt_delta = 0.99)},
#'            number of iterations etc.
#' @return `baggr` class structure: a list including Stan model fit
#'          alongside input data, pooling metrics, various model properties.
#'          If test data is used, mean value of -2*lpd is reported as `mean_lpd`
#'
#'
#' @details
#'
#' Running `baggr` requires 1/ data preparation, 2/ choice of model, 3/ choice of priors.
#' All three are discussed in depth in the package vignette (`vignette("baggr")`).
#'
#' __Data.__ For aggregate data models you need a data frame with columns
#' `tau` and `se` or `tau`, `mu`, `se.tau`, `se.mu`.
#' An additional column can be used to provide labels for each group
#' (by default column `group` is used if available, but this can be
#' customised -- see the example below).
#' For individual level data three columns are needed: outcome, treatment, group. These
#' are identified by using the `outcome`, `treatment` and `group` arguments.
#'
#' When working with individual-level data,
#' many data preparation steps (summarising, standardisation etc.)
#' can be done through a helper function [prepare_ma].
#' Using it will also automatically format data inputs to be
#' work with `baggr()`.
#'
#' __Models.__ Available models are:
#'
#' * for the __continuous variable__ means:
#'   `"rubin"` model for average treatment effect, `"mutau"`
#'   version which takes into account means of control groups, `"full"`,
#'   which works with individual-level data
#' * for __continuous variable quantiles__: `"quantiles"`` model
#'   (see Meager, 2019 in references)
#' * for __binary data__: `"logit"` model can be used on individual-level data;
#'   you can also analyse continuous statistics such as
#'   log odds ratios and logs risk ratios using the models listed above;
#'   see `vignette("baggr_binary")` for tutorial with examples
#'
#'  If no model is specified, the function tries to infer the appropriate
#'  model automatically.
#'  Additionally, the user must specify type of pooling.
#'  The default is always partial pooling.
#'
#' __Priors.__ It is optional to specify priors yourself,
#' as the package will try propose an appropriate
#' prior for the input data if you do not pass a `prior` argument.
#' To set the priors yourself, use `prior_` arguments. For specifying many priors at once
#' (or re-using between models), a single `prior = list(...)` argument can be used instead.
#' Appropriate examples are given in `vignette("baggr")`.
#'
#' __Outputs.__ Standard functions for analysing the `baggr` object are
#'
#' * [treatment_effect] for distribution of hyperparameters
#' * [group_effects] for distributions of group-specific parameters
#' * [fixed_effects] for coefficients in meta-regression
#' * [effect_draw] and [effect_plot] for posterior predictive distributions
#' * [baggr_compare] for comparing multiple `baggr` models
#' * [loocv] for cross-validation
#'
#'
#' @author Witold Wiecek, Rachael Meager
#'
#' @examples
#' df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
#'                         "se" = rep(1, 8),
#'                         "state" = datasets::state.name[1:8])
#' baggr(df_pooled) #baggr automatically detects the input data
#' # correct labels, different pooling & passing some options to Stan
#' baggr(df_pooled, group = "state", pooling = "full", iter = 500)
#' #change the priors:
#' baggr(df_pooled, prior_hypermean = normal(5,5))
#'
#' # "mu & tau" model, using a built-in dataset
#' # prepare_ma() can summarise individual-level data
#' \donttest{
#' microcredit_summary_data <- prepare_ma(microcredit_simplified,
#'                                        outcome = "consumerdurables")
#' baggr(microcredit_summary_data, model = "mutau",
#'       pooling = "partial", prior_hypercor = lkj(1),
#'       prior_hypersd = normal(0,10),
#'       prior_hypermean = multinormal(c(0,0),matrix(c(10,3,3,10),2,2)))
#' }
#'
#'
#' @importFrom rstan summary
#' @importFrom rstan sampling
#' @export

baggr <- function(data, model = NULL, pooling = "partial",
                  effect = NULL,
                  covariates = c(),
                  prior_hypermean = NULL, prior_hypersd = NULL, prior_hypercor=NULL,
                  prior_beta = NULL,
                  # log = FALSE, cfb = FALSE, standardise = FALSE,
                  # baseline = NULL,
                  prior = NULL, ppd = FALSE,
                  test_data = NULL, quantiles = seq(.05, .95, .1),
                  outcome = "outcome", group = "group", treatment = "treatment",
                  warn = TRUE, ...) {

  # check that it is data.frame of at least 1 row
  if(!inherits(data, "data.frame") || nrow(data) == 1)
    stop("data argument must be a data.frame of >1 rows")

  # For now we recommend that users format their data before passing to baggr()
  # data <- prepare_ma(data,
  #                    standardise = standardise, log = log,
  #                    summarise = FALSE, cfb = cfb,
  #                    treatment=treatment, group=group,
  #                    outcome=outcome, baseline=baseline)
  attr(data, "outcome") <- outcome
  attr(data, "group") <- group
  attr(data, "treatment") <- treatment

  stan_data <- convert_inputs(data,
                              model,
                              covariates = covariates,
                              quantiles = quantiles,
                              outcome = outcome,
                              group = group,
                              treatment = treatment,
                              test_data = test_data)
  # model might've been chosen automatically (if NULL)
  # within convert_inptuts(), otherwise it's unchanged
  model <- attr(stan_data, "model")

  # remember number of groups:
  n_groups <- attr(stan_data, "n_groups")

  # labels for what the effect parameters represent:


  if(model == "quantiles"){
    if(is.null(effect))
      effect <- paste0(100*quantiles, "% quantile mean")
    else if(length(effect) == 1)
      effect <- paste0(100*quantiles, "% quantile on ", effect)
    else if(length(length(effect) != length(quantiles)))
      stop("'effect' must be of length 1 or same as number of quantiles")
  }
  if(model == "logit"){
    if(is.null(effect))
      effect <- "logOR"
  }
  # In all other cases we set it to mean
  if(is.null(effect))
    effect <- "mean"

  # pooling type
  if(pooling %in% c("none", "partial", "full")) {
    stan_data[["pooling_type"]] <- switch(pooling,
                                          "none" = 0,
                                          "partial" = 1,
                                          "full" = 2)
  } else {
    stop('Wrong pooling parameter; choose from c("none", "partial", "full")')
  }

  # Prior settings:
  if(is.null(prior))
    prior <- list(hypermean = prior_hypermean,
                  hypercor  = prior_hypercor,
                  hypersd   = prior_hypersd,
                  beta      = prior_beta)
  else {
    if(!is.null(prior_hypermean) || !is.null(prior_beta) ||
       !is.null(prior_hypercor)  || !is.null(prior_hypersd))
      message("Both 'prior' and 'prior_' arguments specified. Using 'prior' only.")
    if(class(prior) != "list" ||
       !all(names(prior) %in% c('hypermean', 'hypercor', 'hypersd', 'beta')))
      stop(paste("Prior argument must be a list with names",
                 "'hypermean', 'hypercor', 'hypersd', 'beta'"))
  }
  # If extracting prior from another model, we need to do a swapsie switcheroo:
  stan_args <- list(...)
  if("formatted_prior" %in% names(stan_args)){
    formatted_prior <- stan_args$formatted_prior
    stan_args$formatted_prior <- NULL
  } else { # extract priors from inputs & fill in missing priors
    formatted_prior <- prepare_prior(prior, data, stan_data, model,
                                     pooling, covariates, quantiles = quantiles)
  }
  for(nm in names(formatted_prior))
    stan_data[[nm]] <- formatted_prior[[nm]]

  stan_args$object <- stanmodels[[model]]

  if(ppd){
    if(pooling == "none")
      stop("Can't use prior predictive distribution with no pooling.")
    stan_data <- remove_data_for_prior_pred(stan_data)
  }
  stan_args$data <- stan_data

  fit <- do.call(rstan::sampling, stan_args)

  result <- list(
    "data" = data,
    "inputs" = stan_data,
    "user_prior" = prior,
    "formatted_prior" = formatted_prior,
    "n_groups" = n_groups,
    "n_parameters" = ifelse(model == "quantiles", length(quantiles), 1),
    "effects" = effect,
    "covariates" = covariates,
    "pooling" = pooling,
    "fit" = fit,
    "model" = model
  )

  class(result) <- c("baggr")

  attr(result, "ppd") <- ppd

  if(grepl("individual", attr(stan_data, "data_type")))
    result$summary_data <- prepare_ma(data,
                                      effect = ifelse(model == "logit", "logOR", "mean"),
                                      group = attr(data, "group"),
                                      treatment = attr(data, "treatment"),
                                      outcome = attr(data, "outcome"))

  if(model == "quantiles")
    result[["quantiles"]]    <- quantiles
  if(!ppd){
    result[["pooling_metric"]] <- pooling(result)
    if(!is.null(test_data)){
      result[["test_data"]]    <- test_data
      result[["mean_lpd"]]     <- -2*mean(rstan::extract(fit, "logpd[1]")[[1]])
    }
  }

  # Check convergence
  rhat <- rstan::summary(fit)$summary[,"Rhat"]
  rhat <- rhat[!is.nan(rhat)] #drop some nonsensical parameters
  if(warn && any(rhat > 1.05))
    warning(paste0("Rhat statistic for ", sum(rhat > 1.05),
                   " parameters exceeded 1.05, with maximum equal to ",
                   round(max(rhat),2), ". This suggests lack of convergence.",
                   "\n No further warning will be issued.",
                   "\n Stan model saved as $fit in the returned object. \n"))

  return(result)
}

check_if_baggr <- function(bg) {
  if(!inherits(bg, "baggr"))
    stop("Object of class 'baggr' required.")
}

remove_data_for_prior_pred <- function(data) {
  scalars_to0 <- c("K", "N", "Nc")
  vectors_to_remove <- c("theta_hat_k", "se_theta_k",
                         "y", "treatment", "site")
  matrices_to_remove <- c("X")
  for(nm in scalars_to0)
    if(!is.null(data[[nm]]))
      data[[nm]] <- 0

  for(nm in vectors_to_remove)
    if(!is.null(data[[nm]]))
      data[[nm]] <- array(0, dim = c(0))

  for(nm in matrices_to_remove)
    if(!is.null(data[[nm]]))
      data[[nm]] <- array(0, dim = c(0,0))

  data
}
