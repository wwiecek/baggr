#' Bayesian aggregate treatment effects model
#'
#' Estimate parameters of an ATE model
#' that's appropriate to the supplied
#' individual- or group-level data.
#' (For overall package help file see \code{?baggr_package})
#'
#' @param data data frame with summary or individual level data to meta-analyse
#' @param model if \code{NULL}, detected automatically from input data
#'              otherwise choose from \code{rubin}, \code{mutau}, \code{individual}
#' @param prior list of prior arguments passed directly to each model (see Details)
#' @param pooling choose from \code{none}, \code{partial} (default) and \code{full}
#' @param joint_prior If \code{TRUE}, \code{mu} and \code{tau} will have joint distribution.
#'                    If \code{FALSE}, they have independent priors. Ignored if no control
#'                    (\code{mu}) data exists.
#' @param standardise logical; determines if data inputs are standardised
#' @param outcome   character; column name in (individual-level) \code{data} with outcome variable values
#' @param grouping  character; column name in \code{data} with grouping factor;
#'                  it's necessary for individual-level data, for summarised data
#'                  it will be used as labels for groups when displaying results
#' @param treatment character; column name in (individual-level) \code{data} with treatment factor;
#' @param quantiles if \code{model = "quantiles"}, a vector indicating which quantiles of data to use
#'                  (with values between 0 and 1)
#' @param test_data data for cross-validation; NULL for no validation, otherwise a data frame
#'                  with the same columns as `data` argument (see \code{\link[baggr]{loocv}} for automation)
#' @param ... extra options passed to Stan function, e.g. \code{control = list(adapt_delta = 0.99)},
#'            number of iterations etc.
#' @return `baggr` class structure: list with Stan model fit embedded inside it,
#'          alongside input data, pooling metrics, various model properties
#'
#' @details
#' This part of documentation is in development.
#'
#' @author Witold Wiecek, Rachael Meager
#'
#' @examples
#' df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
#' "se" = rep(1, 8),
#' "state" = datasets::state.name[1:8])
#' baggr(df_pooled) #automatically detects the input data
#' # correct labels & passing some options to Stan
#' baggr(df_pooled, grouping = "state", iter = 200)
#'
#' @importFrom broom tidyMCMC
#' @export

baggr <- function(data, model = NULL, prior = NULL, pooling = "partial",
                  joint_prior = TRUE, standardise = FALSE,
                  test_data = NULL,
                  quantiles = seq(.05, .95, .1),
                  outcome = "outcome", grouping = "site", treatment = "treatment", ...) {

  stan_data <- convert_inputs(data, model,
                              quantiles = quantiles,
                              outcome = outcome,
                              grouping = grouping,
                              treatment = treatment,
                              standardise = standardise,
                              test_data = test_data)
  # model might've been chosen automatically (if NULL)
  # within convert_inptuts(), otherwise it's unchanged
  model <- attr(stan_data, "model")

  # choice whether the parameters have a joint prior or not
  # for now fixed
  if(model != "rubin")
    stan_data[["joint"]] <- joint_prior

  # remember number of sites:
  n_sites <- attr(stan_data, "n_sites")

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
    "n_sites" = n_sites,
    "effects" = effects,
    "pooling" = pooling,
    "fit" = fit,
    "model" = model
  )

  class(result) <- c("baggr")

  result[["pooling_metric"]] <- pooling(result)
  if(model == "quantiles")
    result[["quantiles"]] <- quantiles

  # Check convergence
  rhat <- rstan::summary(fit)$summary[,"Rhat"]
  rhat <- rhat[!is.nan(rhat)] #drop some nonsensical parameters
  if(any(rhat > 1.05))
    cat(crayon::red(paste0("Rhat statistic for ", sum(rhat > 1.05),
                           " parameters exceeded 1.05, with maximum equal to ",
                           round(max(rhat),2), ". This suggests lack of convergence.",
                           "\n No further warning will be issued. \n")))

  return(result)
}
