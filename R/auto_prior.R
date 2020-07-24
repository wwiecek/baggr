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
      data <- prepare_ma(pma_data, effect = "logOR", rare_event_correction = 0.1)
    if(model == "full")
      data <- prepare_ma(pma_data, effect = "mean")
  }
  # ...then proceed as you would in Rubin model

  if(model %in% c("rubin", "logit", "full", "sslab")) {
    # For each model we need a list of what priors can be specified and
    # what families of distributions are allowed for them
    # (this second part should be changed to just specifying dimensionality/type)
    priors_spec <- list(
      "rubin" = c("hypermean" = "real", "hypersd" = "positive_real"),
      "full"  = c("hypermean" = "real", "hypersd" = "positive_real"),
      "logit"  = c("hypermean" = "real", "hypersd" = "positive_real",
                   "control" = "real", "control_sd" = "positive_real"),
      "sslab"  = c("hypermean" = "real", "hypersd" = "positive_real",
                   "control" = "real", "control_sd" = "positive_real",
                   "control_sigma"  = "real",
                   "control_sigma_sd"  = "real",
                   "sigma"  = "real",
                   "sigma_sd"  = "real",
                   "kappa"  = "real",
                   "kappa_sd"  = "real")

    )

    # "mutau" = list("hypermean"        = list(allowed = c("multinormal"),
    #                                          default = function(data) normal(0, 10*max(abs(data$tau)))),
    #                "hypersd"          = list(allowed = c("cauchy", "normal", "uniform"),
    #                                          default = function(data) normal(0, 10*max(abs(data$tau)))),
    #                "hypercor"         = list(allowed = c("lkj"),
    #                                          default = function(data) normal(0, 10*max(abs(data$tau))))

    for(current_prior in names(priors_spec[[model]])) {
      if(is.null(prior[[current_prior]])) {

        if(model != "sslab")
          default_prior_dist <- switch(current_prior,
                                       "hypermean"  = normal(0, 10*max(abs(data$tau))),
                                       "hypersd"    = uniform(0, 10*sd(data$tau)),
                                       "control"    = normal(0, 10),
                                       "control_sd" = normal(0, 10)
          )

        if(model == "sslab") {

          sd_y_neg_ctrl <- sd(stan_data$y_neg[stan_data$treatment_neg == 0])
          sd_y_neg_trt  <- sd(stan_data$y_neg[stan_data$treatment_neg == 1])
          sd_y_pos_ctrl <- sd(stan_data$y_pos[stan_data$treatment_pos == 0])
          sd_y_pos_trt  <- sd(stan_data$y_pos[stan_data$treatment_pos == 1])
          max_hypermean <- max(c(sd_y_neg_trt, sd_y_pos_trt))
          max_control   <- max(c(sd_y_neg_ctrl, sd_y_pos_ctrl))

          default_prior_dist <- switch(current_prior,
                                       # Prior settings from Rachael, 17/07/2020
                                       "hypermean"        = normal(0, 10*max_hypermean),
                                       "hypersd"          = normal(0, 2.5),
                                       "control"          = normal(0, 10*max_control),
                                       "control_sd"       = normal(0, 10),
                                       "control_sigma"    = normal(0, 10),
                                       "control_sigma_sd" = normal(0, 10),
                                       "sigma"            = normal(0, 10),
                                       "sigma_sd"         = normal(0, 10),
                                       "kappa"            = normal(0, 2.5),
                                       "kappa_sd"         = normal(0, 2.5)
          )
        }

        prior_list <- set_prior_val(prior_list,
                                    paste0("prior_", current_prior),
                                    default_prior_dist)

      } else {
        prior_list <- set_prior_val(prior_list,
                                    paste0("prior_", current_prior),
                                    prior[[current_prior]])
      }

      # Print out priors
      if(!silent) {
        if(is.null(prior[[current_prior]])) {
          if(current_prior == "hypermean") {
            priorname <- ifelse(pooling == "none", "mean in each group", "hypermean")
            effect_on <- switch(model,
                                "logit" = "on log OR ",
                                "sslab" = "in either of the tails",
                                "")
            message(paste0("Setting prior for ", priorname, " using 10 times the max effect ",
                           effect_on, ":"))
            message(paste0("* tau ~ ", print_dist(default_prior_dist)))
          } else if(current_prior == "hypersd") {
            if(nrow(data) < 5)
              message(paste("/Dataset has only", nrow(data),
                            "groups -- consider setting variance prior manually./"))
            if(model != "sslab")
              message(paste0("Setting hyper-SD prior using 10 times the naive SD across sites"))
            message(paste0("* sigma_tau ~ ", print_dist(default_prior_dist)))

          } else if(current_prior == "control" && model == "logit") {
            prop_ctrl <- data$c / (data$c + data$d)
            if(max(prop_ctrl) > .999 | min(prop_ctrl) < .001)
              message("Baseline proportion of events is very low or very common.",
                      "Consider manually setting prior_control.")
            message(paste0("* log odds of event rate in untreated: mean ~ ",
                           print_dist(default_prior_dist)))
          } else if(current_prior == "control_sd" && model == "logit") {
            if(stan_data$pooling_baseline != 0)
              message(paste0("* log odds of event rate in untreated: sd ~ ",
                             print_dist(default_prior_dist)))
          } else {
            message(paste0("* ", current_prior, " ~ ", print_dist(default_prior_dist)))
          }
        } else {
          if(current_prior == "hypersd" && pooling != "partial")
            message("Prior for hyper-SD set, but pooling is not partial. Ignoring.")
          if(current_prior == "control_sd" && model == "logit" && stan_data$pooling_baseline != 0)
            message("SD hyperparameter for control groups defined, but there is no pooling. Ignoring.")
        }
      }
    }

    # Check eligibility
    check_eligible_priors(prior_list, priors_spec[[model]])
  }

  # This is a special case for now.
  if(model == "mutau") {
    # Remember, first row is always mu (baseline), second row is tau (effect)
    # THIS WILL BE REVIDES TO USE control/control_sd type of specification

    # Hypermean
    if(is.null(prior$hypermean)){
      val1 <- 100*max(abs(data$mu))
      val2 <- 100*max(abs(data$tau))
      # Behaviour for joint prior-type behaviour:
      prior_list <- set_prior_val(prior_list, "prior_hypermean",
                                  multinormal(c(0,0), c(val1, val2)*diag(2)))
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

  }

  if(model == "quantiles") {
    prior_list[["prior_dispersion_on_beta_0"]] <- 1000*diag(stan_data$Nq)
    prior_list[["prior_dispersion_on_beta_1"]] <- 1000*diag(stan_data$Nq)
    if(!silent) {
      message("Sigma (hypervariance-covariance) = 1000*I_{Nquantiles}")
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



check_eligible_priors <- function(prior, spec) {
  for(nm in names(spec)){
    if(!(paste0("prior_", nm, "_fam") %in% names(prior)))
      stop(paste("Prior needed for", nm))

    allowed_dist <- available_priors[[spec[[nm]]]]

    if(!any(prior_dist_fam[allowed_dist] == prior[[paste0("prior_", nm, "_fam")]]))
      stop("Prior for ", nm, " must be one of ", allowed_dist)
  }
}
