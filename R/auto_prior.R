# Extracts priors specified by the user;
# Then, if any priors are missing, sets priors for all baggr models

prepare_prior <- function(prior, data, stan_data, model,
                          quantiles = c()) {
  if(missing(prior))
    prior <- list()

  prior_list <- list()

  # message("Automatically setting prior values:")
  if(model %in% c("rubin")) {
    # Hyper-SD
    if(is.null(prior$hypersd)){
      prior_list <- set_prior_val(prior_list, "prior_hypersd", uniform(0, 10*sd(data$tau)))
      if(nrow(data) < 5)
        message(paste("/Dataset has only", nrow(data),
                       "rows -- consider setting variance prior manually./"))

      message(paste0("* Set hypersd: sigma_tau ~ Uniform(0, ",
                     format(10*sd(data$tau), digits = 2), ")"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypersd", prior$hypersd)
    }

    # Hypermean
    if(is.null(prior$hypermean)){
      prior_list <- set_prior_val(prior_list, "prior_hypermean", normal(0, 100))
      message(paste0("* tau ~ Normal(0, 100^2)"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypermean", prior$hypermean)
    }

    check_eligible_priors(prior_list,
                          list("hypersd" = c("normal", "uniform"),
                               "hypermean" = c("normal", "uniform", "cauchy")))
  }

  if(model == "mutau") {
    # Remember, first row is always mu (baseline), second row is tau (effect)

    # Hypermean
    if(is.null(prior$hypermean)){
      prior_list <- set_prior_val(prior_list, "prior_hypermean",
                                  multinormal(c(0,0), 10000*diag(2)))
      message(paste0("* hypermean (mu, tau) ~ Normal([0,0], (1000^2)*Id_2)"))
      if(nrow(data) < 5)
        message(paste("/Dataset has only", nrow(data),
                      "rows -- consider setting variance prior manually./"))
    } else {
      if(prior$hypermean$dist == "normal")
        prior$hypermean <- multinormal(rep(prior$hypermean$values[1], 2),
                                       (prior$hypermean$values[2]^2)*diag(2))
      if(prior$hypermean$dimension != 2)
        stop("Prior for mu & tau model must have 2 dimensions.")
      prior_list <- set_prior_val(prior_list, "prior_hypermean", prior$hypermean)
    }

    # Hypervariance
    if(is.null(prior$hypervar)){
      prior_list <- set_prior_val(prior_list, "prior_hypervar", cauchy(0,10))
      message(paste0("* hypervariance (mu, tau) ~ Cauchy(0,10)"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypervar", prior$hypervar)
    }

    # Hypercorrelation (Only LKJ enabled for now)
    if(is.null(prior$hypercor)){
      prior_list$prior_hypercor_fam <- 4
      prior_list$prior_hypercor_val <- 3
      message(paste0("* hypercorrelation (mu, tau) ~ LKJ(shape=3)"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypercor", prior$hypercor)
    }
    # Make sure prior_list$prior_hypercor_val is an array:
    prior_list$prior_hypercor_val <- array(prior_list$prior_hypercor_val, dim=c(1))

    check_eligible_priors(prior_list,
                          list("hypervar" = c("cauchy", "normal", "uniform"),
                               "hypercor" = c("lkj"),
                               "hypermean" = c("multinormal")))
  }
  if(model == "full") {
    # empirical variance of outcome:
    vhat <- var(stan_data$y) #this may give trouble, look out!
    message(paste0("SD of treatment effect is Uniform(0, ",
                   format(10*sqrt(vhat), digits=2),
                   "); 10*(observed outcome SD)"))
    prior_list[["P"]] <- 2
    prior_list[["mutau_prior_mean"]]  <- rep(0, prior_list$P)
    prior_list[["mutau_prior_sigma"]] <- 100*vhat*diag(prior_list$P)
  }
  if(model == "quantiles") {
    prior_list[["prior_dispersion_on_beta_0"]] <- 1000*diag(stan_data$N)
    prior_list[["prior_dispersion_on_beta_1"]] <- 1000*diag(stan_data$N)
    message("prior_dispersion_on_beta = 1000*diag(stan_data$N)")
  }

  return(prior_list)
}

check_eligible_priors <- function(prior, eligible) {
  for(i in seq_along(eligible)){
    if(!(paste0("prior_", names(eligible)[i], "_fam") %in% names(prior)))
      stop(paste("Prior needed for", names(eligible)[i]))

    allowed_dist <- prior_dist_fam[eligible[[i]]]

    if(!any(prior[[paste0("prior_", names(eligible)[i], "_fam")]] == allowed_dist))
      stop("Prior for ", names(eligible)[i], " must be one of ", eligible[i])
  }
}
