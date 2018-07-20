#' Bayesian aggregate treatment effects model
#'
#' Estimate parameters of an ATE model
#' that's appropriate to the supplied
#' individual- or group-level data.
#'
#' @param data data frame with columns 'outcome', 'treatment' and 'site'
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
#' @param ... extra options passed to Stan function, e.g. \code{control = list(adapt_delta = 0.99)},
#'            number of iterations etc.
#' @param warnings if `TRUE` (default), Stan warnings will be displayed
#' @return `baggr` class structure: list with Stan model fit embedded inside it,
#'          alongside input data, pooling metrics, various model properties
#'
#' @details
#' This part of documentation is in development.
#'
#' @author Witold Wiecek
#' @examples
#' load_baggr_models()
#' df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
#' "se" = rep(1, 8),
#' "state" = datasets::state.name[1:8])
#' baggr(df_pooled) #automatically detects the input data
#' # correct labels
#' baggr(df_pooled, grouping = "state")
#' # Stan options
#' baggr(df_pooled, iter = 200)
#' @export

baggr <- function(data,
                  model = NULL,
                  prior = NULL,
                  pooling = "partial",
                  joint_prior = TRUE,
                  standardise = FALSE,
                  outcome = "outcome", grouping = "site", treatment = "treatment",
                  warnings = TRUE,
                  ...) {



  stan_data <- convert_inputs(data, model,
                              outcome = outcome,
                              grouping = grouping,
                              treatment = treatment,
                              standardise = standardise)
  # model might've been chosen automatically
  # when we prepared inputs, take note:
  model <- attr(stan_data, "model")

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
    if(model == "full") {
      # empirical variance of outcome:
      vhat <- var(stan_data[["y"]])
      message(paste("Prior variance set to 5 times the observed variance in outcome."))
      stan_data[["P"]] <- 2
      stan_data[["mutau_prior_mean"]]  <- rep(0, stan_data$P)
      stan_data[["mutau_prior_sigma"]] <- 5*vhat*diag(stan_data$P)
    }
  } else {
    # !!!check for allowed priors here!!!
    for(nm in names(prior))
      stan_data[[nm]] <- prior[[nm]]
  }

  # if(warnings)
    # fit <- rstan::sampling(stan_model, data = stan_data, ...)
  fit <- rstan::sampling(stanmodels[[model]], data = stan_data, ...)
  # else
    # fit <- suppressWarnings(rstan::sampling(stan_model, data = stan_data, ...))

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
