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
#' @return `baggr` class list with Stan model fit embedded inside it,
#'          alongside input data, pooling metrics, various model properties
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
                  joint_prior = TRUE,
                  warnings = TRUE,
                  standardise = FALSE,
                  outcome = "outcome", grouping = "site", treatment = "treatment",
                  ...) {



  stan_data <- convert_inputs(data, model,
                              outcome = outcome,
                              grouping = grouping,
                              treatment = treatment,
                              standardise = standardise)
  # model might've been chosen automatically
  # when we prepared inputs, take note:
  model <- attr(stan_data, "model")
  stan_model  <- get_model(model)

  # choice whether the parameters have a joint prior or not
  # for now fixed
  if(model != "rubin")
    stan_data[["joint"]] <- joint_prior

  # remember number of sites:
  n_sites <- attr(stan_data, "n_sites")



  # default priors
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

  if(is.null(prior)) {
    message("Automatically setting prior values.")
    if(model %in% c("rubin")) {
      stan_data[["prior_upper_sigma_tau"]] <- 5*var(data$tau)
      message(paste0("sigma_tau ~ Uniform(0, ",
                     round(stan_data[["prior_upper_sigma_tau"]], 2), ")"))
      stan_data[["prior_tau_mean"]] <- 0
      stan_data[["prior_tau_scale"]] <- 1000
      message(paste0("tau ~ Normal(0, 1000)"))
    }
    if(model == "mutau") {
      # Remember, first row is always mu (baseline), second row is tau (effect)
      stan_data[["prior_upper_sigma_tau"]] <- c(5*var(data$mu), 5*var(data$tau))
      message(paste0("sigma_mu ~ Uniform(0, ",
                     round(stan_data[["prior_upper_sigma_tau"]][1], 2), ")"))
      message(paste0("sigma_tau ~ Uniform(0, ",
                     round(stan_data[["prior_upper_sigma_tau"]][2], 2), ")"))
      stan_data[["prior_tau_mean"]] <- c(0,0)
      stan_data[["prior_tau_scale"]] <- 1000*diag(2)
      message(paste0("(mu, tau) ~ Normal([0,0], (1000^2)*Id_2)"))
    }
    if(model == "joint") {
      # empirical variance of outcome:
      vhat <- var(stan_data[["y"]])
      message(paste("Prior variance set to 5 times the observed variance in outcome."))
      stan_data[["P"]] <- 2
      stan_data[["mutau_prior_mean"]]  <- rep(0, stan_data$P)
      stan_data[["mutau_prior_sigma"]] <- (vhat)*diag(stan_data$P)
    }
  } else {
    # !!!check for allowed priors here!!!
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
