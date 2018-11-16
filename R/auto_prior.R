# This automatically sets priors for all baggr models
# That is inferred from data

auto_prior <- function(data, stan_data, model, outcome = "outcome",
                       quantiles = c()) {
  prior_list <- list()

  message("Automatically setting prior values:")

  if(model %in% c("rubin")) {
    prior_list[["prior_upper_sigma_tau"]] <- 10*var(data$tau)
    message(paste0("* sigma_tau ~ Uniform(0, ",
                   round(prior_list[["prior_upper_sigma_tau"]], 2), ")"))
    prior_list[["prior_tau_mean"]] <- 0
    prior_list[["prior_tau_scale"]] <- 1000
    message(paste0("* tau ~ Normal(0, 1000)"))
  }
  if(model == "mutau") {
    # Remember, first row is always mu (baseline), second row is tau (effect)
    prior_list[["prior_upper_sigma_tau"]] <- c(10*sd(data$mu), 10*sd(data$tau))^2
    message(paste0("* sigma_mu ~ Uniform(0, ",
                   round(sqrt(prior_list[["prior_upper_sigma_tau"]][1]), 2), ")"))
    message(paste0("* sigma_tau ~ Uniform(0, ",
                   round(sqrt(prior_list[["prior_upper_sigma_tau"]][2]), 2), ")"))
    prior_list[["prior_tau_mean"]] <- c(0,0)
    prior_list[["prior_tau_scale"]] <- 1000*diag(2)
    message(paste0("mean priors: (mu, tau) ~ Normal([0,0], (1000^2)*Id_2)"))
  }
  if(model == "full") {
    # empirical variance of outcome:
    vhat <- var(data[[outcome]]) #this may give trouble, look out!
    message(paste("* Prior variance set to 5 times the observed variance in outcome."))
    prior_list[["P"]] <- 2
    prior_list[["mutau_prior_mean"]]  <- rep(0, prior_list$P)
    prior_list[["mutau_prior_sigma"]] <- 10*vhat*diag(prior_list$P)
  }
  if(model == "quantiles") {
    prior_list[["prior_dispersion_on_beta_0"]] <- 1000*diag(stan_data$N)
    prior_list[["prior_dispersion_on_beta_1"]] <- 1000*diag(stan_data$N)
    message("prior_dispersion_on_beta = 1000*diag(stan_data$N)")
  }

  return(prior_list)
}
