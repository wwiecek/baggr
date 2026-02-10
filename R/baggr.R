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
#' @param data data frame with summary or individual level data to meta-analyse;
#'             see Details section for how to format your data
#' @param model if \code{NULL}, detected automatically from input data
#'              otherwise choose from `rubin`, `rubin_full`, `logit`, `mutau` (see Details).
#' @param pooling Type of pooling;
#'                choose from \code{"none"}, \code{"partial"} (default) and \code{"full"}.
#'                If you are not familiar with the terms, consult the vignette;
#'                "partial" can be understood as random effects and "full" as fixed effects
#' @param pooling_control Pooling for group-specific control mean terms in models using
#'                        individual-level data. Typically we use either `"none"` or `"partial"`,
#'                        but if you want to remove the group-specific intercept altogether,
#'                        set this to `"remove"`.
#' @param effect_label How to label the effect(s). These labels are used in various print and plot outputs.
#'                     Will default to `"mean"` in most models, `"log OR"` in logistic model etc.
#'                     If you plan on comparing models (see [baggr_compare]), use the same labels.
#' @param covariates Character vector with column names in `data`. The corresponding columns are used as
#'                   covariates (fixed effects) in the meta-regression model (in case of aggregate data).
#'                   In the case of individual level data the model does not differentiate between group-level
#'                   variables (same values of the covariate for all rows related to a given group) and
#'                   individual-level covariates.
#' @param prior_hypermean prior distribution for hypermean; you can use "plain text" notation like
#'              `prior_hypermean=normal(0,100)` or `uniform(-10, 10)`.
#'              See _Details:Priors_ section below for more possible specifications.
#'              If unspecified, the priors will be derived automatically based on data
#'              (and printed out in the console).
#' @param prior_hypersd  prior for hyper-standard deviation, used
#'                       by Rubin and `"mutau"` models;
#'                       same rules apply as for `_hypermean`;
#' @param prior_hypercor prior for hypercorrelation matrix, used by the `"mutau"` model
#' @param prior_beta prior for regression coefficients if `covariates` are specified; will default to
#'                       experimental normal(0, 10^2) distribution
#' @param prior_cluster priors for SDs of cluster random effects in each study
#'                      (i.e. assuming normal(0, sigma_k^2), with different sigma in each `group`)
#' @param prior_control prior for the mean in the control arm (baseline), currently
#'                      used in `"logit"` model only;
#'                      if `pooling_control = "partial"`, the prior is hyperprior
#'                      for all baselines, if `"none"`,
#'                      then it is an independent prior for all baselines
#' @param prior_selection prior for the log of relative publication probability (see __Selection__ section below and argument `selection`)
#' @param prior_control_sd prior for the SD in the control arm (baseline), currently
#'                         used in `"logit"` model only;
#'                         this can only be used if `pooling_control = "partial"`
#' @param prior_sigma prior for error terms in linear regression models (`"rubin_full"` or `"mutau_full"`)
#' @param prior alternative way to specify all of the priors above as a single named list, with `hypermean`,
#'              `hypersd`, `hypercor`, `beta` etc., with names dropping `prior_` prefix,
#'              e.g. `prior = list(hypermean = normal(0,10), beta = uniform(-5, 5))`
#' @param ppd       logical; use prior predictive distribution? (_p.p.d._)
#'                  If `ppd=TRUE`, Stan model will sample from the prior distribution(s)
#'                  and ignore `data` in inference. However, `data` argument might still
#'                  be used to infer the correct model (if `model=NULL`) and to set the
#'                  default priors, therefore you must specify it.
#' @param outcome   column name in \code{data} (used in individual-level only) with outcome variable values
#' @param group     column name in \code{data} with grouping factor;
#'                  it's necessary for individual-level data, for summarised data
#'                  it will be used as labels for groups when displaying results
#' @param treatment column name in (individual-level) \code{data} with treatment factor;
#' @param cluster   optional; column name in (individual-level) \code{data}; if defined,
#'                  random cluster effects will be fitted in each study
#' @param selection optional vector of z-values which implements symmetrical selection (relative probability of publication) model;
#'                  most commonly you'd set this to 1.96 (p=0.05) of `c(1.96, 2.58)` (p=0.05 and p=0.01)
#' @param quantiles if \code{model = "quantiles"}, a vector indicating which quantiles of data to use
#'                  (with values between 0 and 1)
#' @param test_data data for cross-validation; NULL for no validation, otherwise a data frame
#'                  with the same columns as `data` argument. See "Cross-validation" section below.
#' @param silent Whether to silence messages about prior settings and about other automatic behaviour.
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
#' Below we briefly discuss 1/ data preparation, 2/ choice of model, 3/ choice of priors.
#' All three are discussed in more depth in the package vignette, `vignette("baggr")`.
#'
#' __Data.__ For aggregate data models you need a data frame with columns
#' `tau` and `se` (Rubin model) or `tau`, `mu`, `se.tau`, `se.mu` ("mu & tau" model).
#' An additional column can be used to provide labels for each group
#' (by default column `group` is used if available, but this can be
#' customised -- see the example below).
#' For individual level data three columns are needed: outcome, treatment, group. These
#' are identified by using the `outcome`, `treatment` and `group` arguments.
#'
#' Many data preparation steps can be done through a helper function [prepare_ma].
#' It can convert individual to summary-level data, calculate
#' odds/risk ratios (with/without corrections) in binary data, standardise variables and more.
#' Using it will automatically format data inputs to work with `baggr()`.
#'
#'
#' __Models.__ Available models are:
#'
#' * for the __continuous variable__ means:
#'   `"rubin"` model for average treatment effect (using summary data), `"mutau"`
#'   version which takes into account means of control groups (also using summary data),
#'   `"rubin_full"`,  which is the same model as `"rubin"` but works with individual-level data
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
#'
#' __Covariates.__
#' Both aggregate and individual-level data can include extra columns, given by `covariates` argument
#' (specified as a character vector of column names) to be used in regression models.
#'  We also refer to impact of these covariates as _fixed effects_.
#'
#' Two types of covariates may be present in your data:
#'
#' * In `"rubin"` and `"mutau"` models, covariates that __change according to group unit__.
#'   In that case, the model accounting
#'   for the group covariates is a
#'   [meta-regression](https://handbook-5-1.cochrane.org/chapter_9/9_6_4_meta_regression.htm)
#'   model. It can be modelled on summary-level data.
#' * In `"logit"` and `"rubin_full"` models, covariates that __change according to individual unit__.
#'   Then, such a model is often referred to as a mixed model. It has to be
#'   fitted to individual-level data. Note that meta-regression is a special
#'   case of a mixed model for individual-level data.
#'
#'
#' __Priors.__ It is optional to specify priors yourself,
#' as the package will try propose an appropriate
#' prior for the input data if you do not pass a `prior` argument.
#' To set the priors yourself, use `prior_` arguments. For specifying many priors at once
#' (or re-using between models), a single `prior = list(...)` argument can be used instead.
#' Meaning of the prior parameters may slightly change from model to model.
#' Details and examples are given in `vignette("baggr")`.
#' Setting `ppd=TRUE` can be used to obtain prior predictive distributions,
#' which is useful for understanding the prior assumptions,
#' especially useful in conjunction with [effect_plot]. You can also [baggr_compare]
#' different priors by setting `baggr_compare(..., compare="prior")`.
#'
#' __Cross-validation.__ When `test_data` are specified, an extra parameter, the
#' log predictive density, will be returned by the model.
#' (The fitted model itself is the same regardless of whether there are `test_data`.)
#' To understand this parameter, see documentation of [loocv], a function that
#' can be used to assess out of sample prediction of the model using all available data.
#' If using individual-level data model, `test_data` should only include treatment arms
#' of the groups of interest. (This is because in cross-validation we are not typically
#' interested in the model's ability to fit heterogeneity in control arms, but
#' only heterogeneity in treatment arms.)
#' For using aggregate level data, there is no such restriction.
#'
#' __Selection model.__ If the `selection` argument is not `NULL`, `baggr()` fits a symmetric
#' selection-on-z-values model (currently only in the `"rubin"` summary-data model).
#' The values in `selection` are cut-points on |z| = |tau / se|; for example,
#' `selection = c(1.96, 2.58)` gives three intervals, `[0, 1.96)`, `[1.96, 2.58)`,
#' and `[2.58, Inf)`. Note how inequality works: 1.96 is "already" treated as more significant than 1.95.
#' (You should also consider defining it to more significant digits in large datasets.)
#' Each interval has its own relative publication probability
#' (weight), with the highest-|z| interval normalised to 1, and these weights
#' are estimated jointly with the usual Rubin parameters.
#'
#' The selection model assumes that publication depends only on |z| (not on the
#' sign or other study features). Inference can be very sensitive to the choice of
#' cut-points and priors on the selection weights, which you should set manually.
#' For more complex cases you should consider using other methods.
#'
#' __Outputs.__ By default, some outputs are printed. There is also a
#' plot method for _baggr_ objects which you can access via [baggr_plot] (or simply `plot()`).
#' Other standard functions for working with `baggr` object are
#'
#' * [treatment_effect] for distribution of hyperparameters
#' * [group_effects] for distributions of group-specific parameters (alias: [study_effects], we use the two interchangeably)
#' * [fixed_effects] for coefficients in (meta-)regression
#' * [effect_draw] and [effect_plot] for posterior predictive distributions
#' * [baggr_compare] for comparing multiple `baggr` models
#' * [loocv] for cross-validation
#'
#' @examples
#' df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
#'                         "se" = rep(1, 8),
#'                         "state" = datasets::state.name[1:8])
#' baggr(df_pooled) #baggr automatically detects the input data
#' # same model, but with correct labels,
#' # different pooling & passing some options to Stan
#' baggr(df_pooled, group = "state", pooling = "full", iter = 500)
#' # model with non-default (and very informative) priors
#'
#' baggr(df_pooled, prior_hypersd = normal(0, 2))
#'
#' \donttest{
#' # "mu & tau" model, using a built-in dataset
#' # prepare_ma() can summarise individual-level data
#' ms <- microcredit_simplified
#' microcredit_summary_data <- prepare_ma(ms, outcome = "consumption")
#' baggr(microcredit_summary_data, model = "mutau",
#'       iter = 500, #this is just for illustration -- don't set it this low normally!
#'       pooling = "partial", prior_hypercor = lkj(1),
#'       prior_hypersd = normal(0,10),
#'       prior_hypermean = multinormal(c(0,0),matrix(c(10,3,3,10),2,2)))
#' }
#'
#'
#' @importFrom rstan summary
#' @importFrom rstan sampling
#' @export

baggr <- function(data,
                  model = NULL,
                  pooling = c("partial", "none", "full"),
                  effect_label = NULL,
                  covariates = c(),
                  prior_hypermean = NULL, prior_hypersd = NULL, prior_hypercor=NULL,
                  prior_beta = NULL, prior_cluster = NULL,
                  prior_control = NULL, prior_control_sd = NULL,
                  prior_selection = NULL,
                  prior_sigma = NULL,
                  # log = FALSE, cfb = FALSE, standardise = FALSE,
                  # baseline = NULL,
                  prior = NULL, ppd = FALSE,
                  pooling_control = c("none", "partial", "remove"),
                  test_data = NULL, quantiles = seq(.05, .95, .1),
                  outcome = "outcome", group = "group", treatment = "treatment",
                  cluster = NULL,
                  selection = NULL,
                  silent = FALSE, warn = TRUE, ...) {

  # check that it is data.frame of at least 1 row
  # if(!inherits(data, "data.frame") || nrow(data) == 1)
    # stop("data argument must be a data.frame of >1 rows")

  # Match arguments
  pooling <- match.arg(pooling)
  pooling_control <- match.arg(pooling_control)

  # Rename effect_label to effect (old naming convention in code)
  effect <- effect_label

  # For now we recommend that users format their data before passing to baggr()
  # data <- prepare_ma(data,
  #                    standardise = standardise, log = log,
  #                    summarise = FALSE, cfb = cfb,
  #                    treatment=treatment, group=group,
  #                    outcome=outcome, baseline=baseline)

  if(!is.null(model) && (model == "full")){
    message("From v0.6 model 'full' is named 'rubin_full'. Please update your code in the future.")
    model <- "rubin_full"
  }
  if(!is.null(model) && model == "mutau_sum") {
    message("Experimental: using mu & tau model with sum specification (tau = mu + theta)")
    model <- "mutau"
    cumsum_mutau <- 0
  } else {
    cumsum_mutau <- 1
  }

  stan_data <- convert_inputs(data,
                              model,
                              effect,
                              covariates = covariates,
                              quantiles = quantiles,
                              outcome = outcome,
                              group = group,
                              treatment = treatment,
                              cluster = cluster,
                              selection = selection,
                              test_data = test_data,
                              silent = silent)

  # model might've been chosen automatically (if NULL)
  # within convert_inptuts(), otherwise it's unchanged
  model <- attr(stan_data, "model")
  # data might also change if Rubin model requested but mutau
  # type inputs supplied
  data <- attr(stan_data, "data")
  # Lastly, effect changes if binary data were passed
  # and effect arg was not specified
  effect <- attr(stan_data, "effect")
  attr(data, "outcome") <- outcome
  attr(data, "group") <- group
  attr(data, "treatment") <- treatment

  # remember number of groups:
  n_groups <- attr(stan_data, "n_groups")

  # labels for what the effect parameters represent:


  if(model == "quantiles"){
    if(is.null(effect))
      effect <- paste0(100*quantiles, "% quantile mean")
    else if(length(effect) == 1)
      effect <- paste0(100*quantiles, "% quantile on ", effect)
    else if(length(length(effect) != length(quantiles)))
      stop("For quantile models, 'effect' must be of length 1",
           "or same as number of quantiles")
  } else if(model == "logit"){
    if(is.null(effect))
      effect <- "logOR"
  } else if(model == "sslab") {
    if(is.null(effect))
      effect <- c("Location of negative tail",
                  "Location of positive tail",
                  "Scale of negative tail",
                  "Scale of positive tail",
                  "LogOR on being negative",
                  "LogOR on being equal to 0")
    else if(length(effect) != 6)
      stop("For spike & slab models, 'effect' must be of length 6")
  }
  # In all other cases we set it to mean
  if(is.null(effect))
    effect <- "mean"

  # For models with multiple REs (treatment is factor), we now expand it:
  if((model %in% c("logit", "rubin_full")) && (attr(stan_data, "n_re") > 1))
    effect <- paste(effect, attr(stan_data, "treatment_levels"))

  # Number of TE parameters
  # (in the future this can be built into the models):
  n_parameters <- length(effect)


  # Pooling type:
  stan_data[["pooling_type"]] <- switch(pooling,
                                        "none" = 0,
                                        "partial" = 1,
                                        "full" = 2)
  if(model %in% c("logit", "rubin_full", "mutau_full"))
    stan_data[["pooling_baseline"]] <- switch(pooling_control,
                                              "none" = 0,
                                              "partial" = 1,
                                              "remove" = 2)
  if(model == "mutau_full"){
    # For now this model only used for joint prior
    stan_data[["joint_prior_mean"]] <- 1
    stan_data[["joint_prior_variance"]] <- 1
  }
  if(model == "mutau")
    stan_data[["cumsum"]] <- cumsum_mutau



  # Prior settings:
  if(is.null(prior)){
    prior <- list(hypermean = prior_hypermean,
                  hypercor  = prior_hypercor,
                  hypersd   = prior_hypersd,
                  beta      = prior_beta,
                  cluster   = prior_cluster,
                  control   = prior_control,
                  control_sd= prior_control_sd,
                  selection = prior_selection
                  )
  } else {
    if(!is.null(prior_hypermean) || !is.null(prior_beta) || is.null(prior_cluster) ||
       !is.null(prior_control)   || !is.null(prior_control_sd) ||
       !is.null(prior_hypercor)  || !is.null(prior_hypersd))
      message("Both 'prior$' and 'prior_' arguments specified. Using 'prior' only.")
    if(!inherits(prior, "list") ||
       !all(names(prior) %in% c('hypermean', 'hypercor', 'hypersd', 'cluster',
                                'beta', 'control', 'control_sd')))
      warning(paste("Only names used in the prior argument are:",
                    "'hypermean', 'hypercor', 'hypersd', 'cluster',
                    'beta', 'control', 'control_sd'"))
  }

  # If extracting prior from another model, we need to do this switcheroo:
  stan_args <- list(...)
  if("formatted_prior" %in% names(stan_args)){
    formatted_prior <- stan_args$formatted_prior
    stan_args$formatted_prior <- NULL
    prior_dist <- NULL
  } else { # extract priors from inputs & fill in missing priors
    auto_prior <- prepare_prior(prior, data, stan_data, model,
                                     pooling, covariates, selection,
                                     silent = silent)
    formatted_prior <- auto_prior$stan
    prior_dist      <- auto_prior$dist
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


  # ENTIRELY TEMPORARY CHANGE TO HOW THE FULL RUBIN/LOGIT MODEL OPERATES WITH P==1
  if(model %in% c("rubin_full", "logit")) {
    # browser()
    # stan_args$data$P <- 1
    # stan_args$data$treatment <- array(stan_args$data$treatment, c(stan_args$data$N, 1))
    # stan_args$data$test_treatment <- array(stan_args$data$test_treatment, c(stan_args$data$test_N, 1))

    # wip:
    stan_args$data$test_treatment      <- array(stan_args$data$test_treatment,
                                                c(stan_args$data$N_test, attr(stan_data, "n_re")))

    # stan_args$data$prior_hypermean_fam <- array(stan_args$data$prior_hypermean_fam, 1)
    # stan_args$data$prior_hypersd_fam   <- array(stan_args$data$prior_hypersd_fam, 1)
    # stan_args$data$prior_hypermean_val <- list(stan_args$data$prior_hypermean_val)
    # stan_args$data$prior_hypersd_val   <- list(stan_args$data$prior_hypersd_val)
  }

  # If data are small, we use Rubin model and no covariates, you can default refresh to 0
  if(is.null(stan_args[["refresh"]]) && (model == "rubin") && length(covariates) == 0)
    stan_args[["refresh"]] <- 0

  # Cluster-level random effects are latent nuisance parameters in these full IPD models.
  # Exclude them from monitored outputs unless user explicitly controls monitored pars.
  if(model %in% c("rubin_full", "logit") &&
     is.null(stan_args[["pars"]]) &&
     is.null(stan_args[["include"]])) {
    stan_args[["pars"]] <- "eta_cluster"
    stan_args[["include"]] <- FALSE
  }

  # SAMPLING IS HERE
  fit <- do.call(rstan::sampling, stan_args)



  result <- list(
    "data" = data,
    "inputs" = stan_data,
    "user_prior" = prior,
    "prior_dist" = prior_dist,
    "formatted_prior" = formatted_prior,
    "n_groups" = n_groups,
    "n_parameters" = n_parameters,
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
                                      # rare_event_correction = 0.5,
                                      effect = ifelse(model == "logit", "logOR", "mean"),
                                      group = attr(data, "group"),
                                      treatment = attr(data, "treatment"),
                                      outcome = attr(data, "outcome"),
                                      pooling = TRUE)

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
    warning(paste0("\nRhat statistic for ", sum(rhat > 1.05),
                   " parameters exceeded 1.05, with maximum equal to ",
                   round(max(rhat),2),
                   ".\n This suggests lack of convergence.",
                   "\n No further warning will be issued.",
                   "\n Stan model saved as $fit in the returned object. \n"))

  return(result)
}

check_if_baggr <- function(bg) {
  if(!inherits(bg, "baggr"))
    stop("Object of class 'baggr' required.")
}

remove_data_for_prior_pred <- function(data) {
  scalars_to0 <- c("K", "N", "Nc",
                   # specific to sslab (generalise this)
                   "N_neg", "N_pos")
  objects_to_remove <- c("theta_hat_k", "se_theta_k",
                         "y", "site",
                         # specific to sslab
                         "treatment_neg", "treatment_pos", "cat",
                         "site_neg", "site_pos", "y_neg", "y_pos",
                         "X")

  # Specific to quantiles:
  matrices_to_rescale <- c("y_0", "y_1")
  arrays_to_rescale <- c("Sigma_y_k_0", "Sigma_y_k_1")
  matrices_remove_rows <- c("x")

  # Treatment can be matrix or vector, so:
  if(is.vector(data[["treatment"]]))
    objects_to_remove <- c(objects_to_remove, "treatment")
  if(is.matrix(data[["treatment"]]))
    matrices_remove_rows <- c(matrices_remove_rows, "treatment")

  for(nm in scalars_to0)
    if(!is.null(data[[nm]]))
      data[[nm]] <- 0

  for(nm in objects_to_remove){
    if(!is.null(data[[nm]])){
      # Make sure we're replacing it with the right data type
      cdim <- ifelse(is.null(dim(data[[nm]])), 1, length(dim(data[[nm]])))
      data[[nm]] <- array(0, dim = rep(0, cdim))
    }
  }

  # for(nm in matrices_to_remove)
  #   if(!is.null(data[[nm]]))
  #     data[[nm]] <- array(0, dim = c(0,0))

  for(nm in matrices_remove_rows)
    if(!is.null(data[[nm]]))
      data[[nm]] <- array(0, dim = c(0,ncol(data[[nm]])))

  # This is for quantiles model where you have K x Nq inputs
  # i.e. sites x Nquantiles
  # For PPD we will need 0 x Nq
  for(nm in matrices_to_rescale)
    if(!is.null(data[[nm]]))
      data[[nm]] <- array(0, dim = c(0, data$Nq))
  for(nm in arrays_to_rescale)
    if(!is.null(data[[nm]]))
      data[[nm]] <- array(0, dim = c(0, data$Nq, data$Nq))

  data
}
