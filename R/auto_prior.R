#' Prepare prior values for Stan models in baggr
#'
#' This is an internal function called by [baggr]. You can use it for debugging
#' or to run modified models.
#' It extracts and prepares priors passed by the user.
#' Then, if any necessary priors are missing, it sets them automatically
#' and notifies user about these automatic choices.
#'
#' @param prior `prior` argument passed from [baggr] call
#' @param data  `data` another argument in [baggr]
#' @param stan_data list of inputs that will be used by sampler
#'                  this is already pre-obtained through [convert_inputs]
#' @param model same as in [baggr]
#' @param pooling same as in [baggr]
#' @param covariates same as in [baggr]
#' @param quantiles  same as in [baggr]
#' @param silent same as in [baggr]
#'
#' @return A named list with prior values that can be appended to `stan_data`
#'         and passed to a Stan model.
#'

prepare_prior <- function(prior, data, stan_data, model, pooling, covariates,
                          quantiles = c(), silent = FALSE) {

  if(missing(prior))
    prior <- list()

  prior_list <- list()

  # Swap data for summary data...
  if(model %in% c("full", "logit")) {
    pma_data <- data.frame(outcome = stan_data$y,
                           group = stan_data$site,
                           treatment = stan_data$treatment)
    if(model == "logit")
      data <- prepare_ma(pma_data, effect = "logOR")
    if(model == "full")
      data <- prepare_ma(pma_data, effect = "mean")
  }
  # ...then proceed as you would in Rubin model

  if(model %in% c("rubin", "logit", "full")) {
    # Hypermean
    if(is.null(prior$hypermean)){
      val <- 10*max(abs(data$tau))
      prior_list <- set_prior_val(prior_list, "prior_hypermean", normal(0, val))
      priorname <- ifelse(pooling == "none", "mean in each group", "hypermean")

      effect_on <- ifelse(model == "logit", "on log OR ", "")

      if(!silent) {
        message(paste0("Setting prior for ", priorname, " according to max effect ",
                       effect_on, "(",
                       format(val/10, digits = 2), "):"))
        message(paste0("* tau ~ Normal(0, (10*",
                       format(val/10, digits = 2), ")^2)"))
      }
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypermean", prior$hypermean)
    }

    # Hyper-SD
    if(is.null(prior$hypersd)){
      if(pooling == "partial") {
        prior_list <- set_prior_val(prior_list, "prior_hypersd", uniform(0, 10*sd(data$tau)))
        if(nrow(data) < 5)
          message(paste("/Dataset has only", nrow(data),
                        "groups -- consider setting variance prior manually./"))
        if(!silent) {
          message(paste0("Setting hyper-SD prior using 10 times the naive SD across sites (",
                         format(10*sd(data$tau), digits = 2), ")"))
          message(paste0("* sigma_tau ~ Uniform(0, ",
                         format(10*sd(data$tau), digits = 2), ")"))
        }
      } else {
        prior_list <- set_prior_val(prior_list, "prior_hypersd", uniform(0, 0))
      }
    } else {
      if(pooling == "full")
        message("Prior for hyper-SD set, but pooling is set to full. Ignoring SD prior.")
      prior_list <- set_prior_val(prior_list, "prior_hypersd", prior$hypersd)
    }

    # Controls
    if(model == "logit"){
      # Remember that at this stage data is summary data, so we can estimate control props
      prop_ctrl <- data$c / (data$c + data$d)
      # val <- max(abs(qlogis(prop_ctrl)))
      # But I wouldn't guess from data

      if(max(prop_ctrl) > .999 | min(prop_ctrl) < .001)
        message("Baseline proportion of events is very low or very common. Consider manually setting the prior.")

      prior_list <- set_prior_val(prior_list, "prior_hbasesd", normal(0, 10))

      if(is.null(prior$control)){
        prior_list <- set_prior_val(prior_list, "prior_hbasemean", normal(0, 10))

        if(!silent)
          message(paste0("* log odds of event rate in untreated: mean ~ normal(0, 10^2)"))
          if(stan_data$pooling_baseline != 0)
            message(paste0("sd ~ normal(0, 10)"))

      } else {
        prior_list <- set_prior_val(prior_list, "prior_hbasemean", prior$control)
      }
    }
    check_eligible_priors(prior_list,
                          list("hypersd"   = c("normal", "uniform", "cauchy"),
                               "hypermean" = c("normal", "uniform", "cauchy")))
  }

  if(model == "mutau") {
    # Remember, first row is always mu (baseline), second row is tau (effect)

    # Hypermean
    if(is.null(prior$hypermean)){
      val1 <- 100*max(abs(data$mu))
      val2 <- 100*max(abs(data$tau))
      # Behaviour for joint prior-type behaviour:
      prior_list <- set_prior_val(prior_list, "prior_hypermean",
                                  multinormal(c(0,0), c(val1, val2)*diag(2)))
      if(!silent) {

        message("Set hypermean prior according to max effect:")
        message(paste0("* hypermean (mu, tau) ~ Normal([0,0], [",
                       format(val1, digits = 2), ", ",
                       format(val2, digits = 2), "]*Id_2)"))
      }
      if(nrow(data) < 5)
        message(paste("/Dataset has only", nrow(data),
                      "groups -- consider setting variance prior manually./"))
    } else {
      if(prior$hypermean$dist == "normal")
        prior$hypermean <- multinormal(rep(prior$hypermean$values[1], 2),
                                       (prior$hypermean$values[2]^2)*diag(2))
      if(prior$hypermean$dimension != 2)
        stop("Prior for mu & tau model must have 2 dimensions.")
      prior_list <- set_prior_val(prior_list, "prior_hypermean", prior$hypermean)
    }

    # Hyper-SD
    if(is.null(prior$hypersd)){
      val <- max(10*sd(data$mu), 10*sd(data$tau))
      prior_list <- set_prior_val(prior_list, "prior_hypersd", cauchy(0,val))

      if(!silent) {
        message(paste0("Set hyper-SD prior using 10 times the naive SD across sites (",
                       format(val, digits = 2), ")"))
        message(paste0("* hyper-SD (mu, tau) ~ Cauchy(0,",
                       format(val, digits = 2), ") (i.i.d.)"))
      }
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypersd", prior$hypersd)
    }

    # Hypercorrelation (Only LKJ enabled for now)
    if(is.null(prior$hypercor)){
      prior_list$prior_hypercor_fam <- 4
      prior_list$prior_hypercor_val <- 3
      if(!silent) {
        message(paste0("* hypercorrelation (mu, tau) ~ LKJ(shape=3)"))
      }
    } else {
      prior_list <- set_prior_val(prior_list, "prior_hypercor", prior$hypercor)
    }
    # Make sure prior_list$prior_hypercor_val is an array:
    prior_list$prior_hypercor_val <- array(prior_list$prior_hypercor_val, dim=c(1))

    check_eligible_priors(prior_list,
                          list("hypersd"   = c("cauchy", "normal", "uniform"),
                               "hypercor"  = c("lkj"),
                               "hypermean" = c("multinormal")))
  }
  # if(model == "full") {
  #   # empirical variance of outcome:
  #   vhat <- var(stan_data$y) #this may give trouble, look out!
  #   message(paste0("SD of treatment effect is Uniform(0, ",
  #                  format(10*sqrt(vhat), digits=2),
  #                  "); 10*(observed outcome SD)"))
  #   prior_list[["joint"]] <- 1
  #   prior_list[["P"]] <- 2
  #   prior_list[["mutau_prior_mean"]]  <- rep(0, prior_list$P)
  #   prior_list[["mutau_prior_sigma"]] <- 100*vhat*diag(prior_list$P)
  # }

  if(model == "quantiles") {
    prior_list[["prior_dispersion_on_beta_0"]] <- 1000*diag(stan_data$N)
    prior_list[["prior_dispersion_on_beta_1"]] <- 1000*diag(stan_data$N)
    if(!silent) {
      message("prior_dispersion_on_beta = 1000*diag(stan_data$N)")
    }
  }

  # Setting covariates prior
  if(length(covariates) > 0) {
    if(is.null(prior$beta)){
      val <- max(10*sd(data$mu), 10*sd(data$tau))
      prior_list <- set_prior_val(prior_list, "prior_beta", normal(0, 10))
      message(paste0("Set beta prior (on covariates in regression)to N(0, 10^2)",
                     " -- purely experimental, use with caution"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_beta", prior$beta)
    }
  } else {
    prior_list <- set_prior_val(prior_list, "prior_beta", uniform(0,0))
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
