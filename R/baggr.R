#' Bayesian aggregate treatment effects model
#'
#' Estimate parameters of an ATE model
#' that's appropriate to the supplied
#' individual- or group-level data.
#' (For overall package help file see \code{?baggr_package})
#'
#' @param data data frame with summary or individual level data to meta-analyse
#' @param model if \code{NULL}, detected automatically from input data
#'              otherwise choose from \code{"rubin"}, \code{"mutau"}, \code{"individual"}
#' @param pooling choose from \code{none}, \code{partial} (default) and \code{full}
#' @param prior list of prior arguments passed directly to each model (see Details)
#' @param joint_prior If \code{TRUE}, \code{mu} and \code{tau} will have joint distribution.
#'                    If \code{FALSE}, they have independent priors. Ignored if no control
#'                    (\code{mu}) data exists.
#' @param outcome   character; column name in (individual-level) \code{data} with outcome variable values
#' @param group     character; column name in \code{data} with grouping factor;
#'                  it's necessary for individual-level data, for summarised data
#'                  it will be used as labels for groups when displaying results
#' @param treatment character; column name in (individual-level) \code{data} with treatment factor;
#' @param quantiles if \code{model = "quantiles"}, a vector indicating which quantiles of data to use
#'                  (with values between 0 and 1)
#' @param test_data data for cross-validation; NULL for no validation, otherwise a data frame
#'                  with the same columns as `data` argument (see \code{\link[baggr]{loocv}} for automation)
#' @param warn print warning if Rhat exceeds 1.05
#' @param ... extra options passed to Stan function, e.g. \code{control = list(adapt_delta = 0.99)},
#'            number of iterations etc.
#' @return `baggr` class structure: list with Stan model fit embedded inside it,
#'          alongside input data, pooling metrics, various model properties
#'
#' @details
#'
#' Running `baggr` requires 1/ data preparation, 2/ choice of model, 3/ choice of priors.
#' All three are discussed in depth in [the package vignette](baggr.html).
#'
#' __Data.__ For aggregate data models you need a data frame with columns `tau` and `se` or `tau`, `mu`, `se.tau`, `se.mu`.
#' For individual level data three columns are needed: outcome, treatment, group. These
#' are identified by using the `outcome`, `treatment` and `group` arguments.
#'
#' Most of data preparation steps (summarising, standardisation etc.) can be done
#' automatically through \code{\link[baggr]{prepare_ma}}.
#' While the preparation step is optional it will also automatically format data inputs to be
#' immediately recognisable by `baggr()`.
#'
#' __Models.__ Available models are:
#'
#' * for the means: `"rubin"`` model for average treatment effect, `"mutau"` version which takes into account means in control group, `"full"`` model which reduces to "mu and tau" (if no covariates are used)
#' * "quantiles" model is also available (see Meager, 2019 and [vignette](baggr.html) for details)
#'
#'  If no model is specified, the function tries to infer the appropriate model automatically.
#'
#' __Priors.__ It is optional to specify priors yourself, as the package will try propose an appropriate
#' prior for the input data if `prior=NULL`.
#' To priors yourself, please refer to the list in the [vignette](baggr.html).
#'
#'
#'
#' @author Witold Wiecek, Rachael Meager
#'
#' @examples
#' df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
#' "se" = rep(1, 8),
#' "state" = datasets::state.name[1:8])
#' baggr(df_pooled) #automatically detects the input data
#' # correct labels & passing some options to Stan
#' baggr(df_pooled, group = "state", iter = 200)
#'
#' @export

baggr <- function(data, model = NULL, prior = NULL, pooling = "partial",
                  # log = FALSE, cfb = FALSE, standardise = FALSE,
                  # baseline = NULL,
                  joint_prior = TRUE,
                  test_data = NULL, quantiles = seq(.05, .95, .1),
                  outcome = "outcome", group = "group", treatment = "treatment",
                  warn = TRUE, ...) {

  # For now we recommend that users format their data before passing to baggr()
  # data <- prepare_ma(data,
  #                    standardise = standardise, log = log,
  #                    summarise = FALSE, cfb = cfb,
  #                    treatment=treatment, group=group,
  #                    outcome=outcome, baseline=baseline)

  stan_data <- convert_inputs(data,
                              model,
                              quantiles = quantiles,
                              outcome = outcome,
                              group = group,
                              treatment = treatment,
                              test_data = test_data)
  # model might've been chosen automatically (if NULL)
  # within convert_inptuts(), otherwise it's unchanged
  model <- attr(stan_data, "model")


  # choice whether the parameters have a joint prior or not
  # for now fixed
  if(model != "rubin")
    stan_data[["joint"]] <- joint_prior

  # remember number of groups:
  n_groups <- attr(stan_data, "n_groups")

  # labels for what the effect parameters represent:
  if(model == "quantiles")
    effects <- paste0(100*quantiles, "% quantile mean")
  else
    effects <- "mean"

  # pooling type
  if(pooling %in% c("none", "partial", "full")) {
    # if(model %in% c("rubin", "mutau")) {
    # switch? separate scripts for each pooling type? third way?
    stan_data[["pooling_type"]] <- switch(pooling,
                                          "none" = 0,
                                          "partial" = 1,
                                          "full" = 2)
  } else {
    stop('Wrong pooling parameter; choose from c("none", "partial", "full")')
  }

  # default priors
  if(is.null(prior)) {
    prior <- auto_prior(data, stan_data, model,
                        outcome = outcome, quantiles = quantiles)

  } else {
    # !!!check for allowed priors here!!!
  }
  for(nm in names(prior))
    stan_data[[nm]] <- prior[[nm]]

  fit <- rstan::sampling(stanmodels[[model]], data = stan_data, ...)

  result <- list(
    "data" = data,
    "inputs" = stan_data,
    "prior" = prior,
    "n_groups" = n_groups,
    "n_parameters" = ifelse(model == "quantiles", length(quantiles), 1),
    "effects" = effects,
    "pooling" = pooling,
    "fit" = fit,
    "model" = model
  )

  class(result) <- c("baggr")

  result[["pooling_metric"]] <- pooling(result)
  if(model == "quantiles")
    result[["quantiles"]] <- quantiles
  if(!is.null(test_data)){
    result[["test_data"]] <- test_data
    result[["mean_lpd"]] <- mean(rstan::extract(fit, "logpd")[[1]])
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
