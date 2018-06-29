#' Bayesian aggregate treatment effects model
#'
#' Estimate parameters in an ATE model
#' that's appropriate to the supplied
#' individual- or group-level data.
#'
#' @param data data frame with columns 'outcome', 'treatment' and 'site'
#' @param model choose from 'rubin', ...
#' @param pooling choose from `none`, `partial` (default) and `full`
#' @param model choose from 'rubin', ...
#' @param prior WIP
#' @param warnings if `TRUE` (default), Stan warnings will be displayed
#'
#' @return `stanfit` object with appropriate parameters
#'
#' @details
#' WIP
#' @author Witold Wiecek
#'
#' @import rstan
#' @export

baggr <- function(data,
                  model = NULL,
                  prior = NULL,
                  pooling = "partial",
                  warnings = TRUE,
                  outcome = "outcome", grouping = "site", treatment = "treatment",
                  ...) {



  stan_data <- convert_inputs(data, model,
                              outcome = outcome,
                              grouping = grouping,
                              treatment = treatment)
  # model might've been chosen automatically
  # when we prepared inputs, take note:
  model <- attr(stan_data, "model")
  stan_model  <- get_model(model)

  # remember number of sites:
  n_sites <- attr(stan_data, "n_sites")



  # default priors
  if(pooling %in% c("none", "partial", "full")) {
    if(model == "rubin") {
      # switch? separate scripts for each pooling type? third way?
      stan_data[["pooling_type"]] <- switch(pooling,
                                            "none" = 0,
                                            "partial" = 1,
                                            "full" = 2)
    }
  } else {
    stop('Wrong pooling parameter; choose from c("none", "partial", "full")')
  }

  if(is.null(prior)) {
    message("Automatically setting prior values.")
    if(model == "rubin") {
      stan_data[["prior_max"]] <- 5*var(data$tau)
      message(paste0("sigma_tau ~ Uniform(0, ", round(stan_data[["prior_max"]], 2), ")"))
    }
    if(model == "joint") {
      # empirical variance of outcome:
      vhat <- var(stan_data[["y"]])
      message(paste(
        "Prior variance set to 5 times the observed variance in outcome."))
      stan_data[["P"]] <- 2
      stan_data[["mutau_prior_mean"]]  <- rep(0, stan_data$P)
      stan_data[["mutau_prior_sigma"]] <- (vhat)*diag(stan_data$P)
    }
  } else {
    # check for allowed priors here
    for(nm in names(prior))
      stan_data[[nm]] <- prior[[nm]]
  }
  if(warnings)
    fit <- rstan::sampling(stan_model, data = stan_data, ...)
  else
    fit <- suppressWarnings(rstan::sampling(stan_model, data = stan_data, ...))

  result <- list(
    "data" = data,
    "inputs" = stan_data,
    "n_sites" = n_sites,
    "pooling" = pooling,
    "fit" = fit,
    "model" = model
  )

  class(result) <- c("baggr")

  result[["pooling_metric"]] <- pooling(result)

  return(result)
}
