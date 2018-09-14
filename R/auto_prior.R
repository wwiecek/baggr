# This automatically sets priors for all baggr models
# That is inferred from data

auto_prior <- function(data, model, outcome = "outcome") {
  stan_data <- list()

  message("Automatically setting prior values:")

  if(model %in% c("rubin")) {
    stan_data[["prior_upper_sigma_tau"]] <- 5*var(data$tau)
    message(paste0("* sigma_tau ~ Uniform(0, ",
                   round(stan_data[["prior_upper_sigma_tau"]], 2), ")"))
    stan_data[["prior_tau_mean"]] <- 0
    stan_data[["prior_tau_scale"]] <- 1000
    message(paste0("* tau ~ Normal(0, 1000)"))
  }
  if(model == "mutau") {
    # Remember, first row is always mu (baseline), second row is tau (effect)
    stan_data[["prior_upper_sigma_tau"]] <- c(5*var(data$mu), 5*var(data$tau))
    message(paste0("* sigma_mu ~ Uniform(0, ",
                   round(stan_data[["prior_upper_sigma_tau"]][1], 2), ")"))
    message(paste0("* sigma_tau ~ Uniform(0, ",
                   round(stan_data[["prior_upper_sigma_tau"]][2], 2), ")"))
    stan_data[["prior_tau_mean"]] <- c(0,0)
    stan_data[["prior_tau_scale"]] <- 1000*diag(2)
    message(paste0("(mu, tau) ~ Normal([0,0], (1000^2)*Id_2)"))
  }
  if(model == "full") {
    # empirical variance of outcome:
    vhat <- var(data[[outcome]]) #this may give trouble, look out!
    message(paste("* Prior variance set to 5 times the observed variance in outcome."))
    stan_data[["P"]] <- 2
    stan_data[["mutau_prior_mean"]]  <- rep(0, stan_data$P)
    stan_data[["mutau_prior_sigma"]] <- 5*vhat*diag(stan_data$P)
  }

  return(stan_data)
}
