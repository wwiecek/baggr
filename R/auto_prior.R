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

  if(model %in% c("mutau", "mutau_full"))
    re_dim <- 2
  else
    re_dim <- 1

  # Swap data for summary data...
  if(model %in% c("rubin_full", "logit", "mutau_full")) {
    pma_data <- data.frame(outcome = stan_data$y,
                           group = stan_data$site,
                           treatment = stan_data$treatment)
    if(model == "logit"){
      data <- prepare_ma(pma_data, effect = "logOR", rare_event_correction = 0.1)
      # Must add mu for automatic definition of baseline prior
      p <- data$c/data$n2
      data$mu <- log(p/(1-p))
    }
    if(model %in% c("mutau_full", "rubin_full"))
      data <- prepare_ma(pma_data, effect = "mean")
  }
  # ...then proceed as you would in Rubin model

  if(model %in% c("rubin", "mutau", "logit", "rubin_full", "sslab", "mutau_full")) {
    # For each model we need a list of what priors can be specified and
    # what families of distributions are allowed for them
    # (this second part should be changed to just specifying dimensionality/type)
    # This must match names.R
    priors_spec <- list(
      "rubin" = c("hypermean" = "real", "hypersd" = "positive_real"),
      "rubin_full"  = c("hypermean" = "real", "hypersd" = "positive_real",
                        "control" = "real", "control_sd" = "positive_real"),
      "mutau" = c("hypermean" = "real_2",
                  "hypercor" = "corr",
                  "hypersd" = "positive_real"),
      "mutau_full" = c("hypermean" = "real_2",
                       "hypercor" = "corr",
                       "hypersd" = "positive_real",
                       "control_sd" = "positive_real"),
      "logit"  = c("hypermean" = "real", "hypersd" = "positive_real",
                   "control" = "real", "control_sd" = "positive_real"),
      "sslab"  = c("hypermean" = "real", "hypersd" = "positive_real",
                   "control" = "real", "control_sd" = "positive_real",
                   "scale_control"  = "real",
                   "scale_control_sd"  = "positive_real",
                   "scale"  = "real",
                   "scale_sd"  = "positive_real",
                   "kappa"  = "real",
                   "kappa_sd"  = "positive_real")
    )

    # "mutau" = list("hypermean" = list(allowed = c("multinormal"),
    #                                   default = function(data) normal(0, 10*max(abs(data$tau)))),
    #                "hypersd"   = list(allowed = c("cauchy", "normal", "uniform"),
    #                                   default = function(data) normal(0, 10*max(abs(data$tau)))),
    #                "hypercor"  = list(allowed = c("lkj"),
    #                                   default = function(data) normal(0, 10*max(abs(data$tau))))

    for(current_prior in names(priors_spec[[model]])) {
      if(is.null(prior[[current_prior]])) {
        if(model != "sslab")
          default_prior_dist <- switch(current_prior,
                                       "hypermean"  = if(re_dim == 1)
                                         normal(0, 10*max(abs(data$tau)))
                                       else
                                         multinormal(c(0,0),
                                                     c(100*max(abs(data$mu)),
                                                       100*max(abs(data$tau)))*diag(2)),
                                       "hypersd"    = uniform(0, 10*sd(data$tau)),
                                       "hypercor"   = lkj(3),
                                       "control"    = normal(0, 10*max(abs(data$mu))),
                                       "control_sd" = uniform(0, 10*sd(data$mu))
          )

        if(model == "sslab") {
          data_pos <- prepare_ma(data.frame(outcome = stan_data$y_pos,
                                            group = stan_data$site_pos,
                                            treatment = stan_data$treatment_pos),
                                 effect = "mean")
          data_neg <- prepare_ma(data.frame(outcome = stan_data$y_neg,
                                            group = stan_data$site_neg,
                                            treatment = stan_data$treatment_neg),
                                 effect = "mean")
          max_hypermean <- max(c(abs(data_pos$tau), abs(data_neg$tau)))
          max_control   <- max(c(abs(data_pos$mu),  abs(data_neg$mu)))
          sd_hypermean  <- 10*max(c(sd(data_pos$mu),  sd(data_neg$mu)))
          sd_control    <- 10*max(c(sd(data_pos$tau), sd(data_neg$tau)))

          default_prior_dist <- switch(current_prior,
                                       # Prior settings from Rachael, 17/07/2020
                                       "hypermean"        = normal(0, 10*max_hypermean),
                                       "hypersd"          = normal(0, 10*sd_hypermean),
                                       "control"          = normal(0, 10*max_control),
                                       "control_sd"       = normal(0, 10*sd_control),
                                       "scale_control"    = normal(0, 10),
                                       "scale_control_sd" = normal(0, 10),
                                       "scale"            = normal(0, 5),
                                       "scale_sd"         = normal(0, 5),
                                       "kappa"            = normal(0, 5),
                                       "kappa_sd"         = normal(0, 2.5)
          )
        }

        prior_list <- set_prior_val(prior_list,
                                    paste0("prior_", current_prior),
                                    default_prior_dist)

      } else {

        # Expand Normal() to Normal_p() when P > 1
        if(re_dim > 1 && prior[[current_prior]][["dist"]] == "normal")
          prior$hypermean <- multinormal(rep(prior[[current_prior]]$values[1], 2),
                                         (prior[[current_prior]]$values[2]^2)*diag(2))

        prior_list <- set_prior_val(prior_list,
                                    paste0("prior_", current_prior),
                                    prior[[current_prior]])
      }

      # Print out priors
      if(!silent) {
        special_name <- ""

        if(is.null(prior[[current_prior]])) {

          # 1) Go through various automated prompts
          if(current_prior == "hypermean") {
            priorname <- ifelse(pooling == "none", "mean in each group", "hypermean")
            effect_on <- switch(model,
                                "logit" = "on log OR ",
                                "sslab" = "in either of the tails",
                                "")
            message(paste0("Setting prior for ", priorname,
                           " using 10 times the max effect ",
                           effect_on, ":"))
            if(model != "sslab")
              special_name <- "tau"
          } else if(current_prior == "hypersd") {
            if(nrow(data) < 5)
              message(paste("/Dataset has only", nrow(data),
                            "groups -- consider setting variance prior manually./"))
            if(model != "sslab"){
              message("Setting hyper-SD prior using 10 times the naive SD across sites")
              special_name <- "sigma_tau"
            }
          } else if(current_prior == "control" && model == "logit") {
            prop_ctrl <- data$c / (data$c + data$d)
            if(max(prop_ctrl) > .999 | min(prop_ctrl) < .001)
              message("Baseline proportion of events is very close to 0 or 1.",
                      "Consider manually setting prior_control.")
            special_name <- "log odds of event rate in untreated: mean"
          } else if(current_prior == "control_sd" && model == "logit") {
            if(stan_data$pooling_baseline != 0)
              special_name <- "log odds of event rate in untreated: sd"
          }

          # 2) Print the prior:

          if(special_name == "")
            message(paste0("* ", current_prior, " ~ ",
                           print_dist(default_prior_dist)))
          else
            message(paste0("* ", current_prior, " [", special_name, "] ~ ",
                           print_dist(default_prior_dist)))

        } else {
          if(current_prior == "hypersd" &&
             pooling != "partial")
            message("Prior for hyper-SD set, but pooling is not partial. Ignoring.")
          if(current_prior == "control_sd" &&
             model == "logit" &&
             stan_data$pooling_baseline != 0)
            message("SD hyperparameter for control groups defined,",
                    "but there is no pooling. Ignoring.")
        }
      }
    }

    # Check eligibility
    check_eligible_priors(prior_list, priors_spec[[model]])
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
      val <- max(10*unlist(lapply(data[, covariates, drop = FALSE], function(x)
        sd(as.numeric(as.factor(x))))))
      prior_list <- set_prior_val(prior_list, "prior_beta", normal(0, val))
      message(
        paste0("Setting prior for covariates in regression to normal,",
               " with SD equal to 10*(highest SD among covariates):\n",
               "* beta ~ ", print_dist(normal(0, val)),
               # Normal(0, ", format(val, digits = 2), "^2)",
               " -- purely experimental, use with caution"))
    } else {
      prior_list <- set_prior_val(prior_list, "prior_beta", prior$beta)
    }
  } else {
    prior_list <- set_prior_val(prior_list, "prior_beta", uniform(0, 1))
  }

  return(prior_list)
}



check_eligible_priors <- function(prior, spec) {
  for(nm in names(spec)){
    if(!(paste0("prior_", nm, "_fam") %in% names(prior)))
      stop(paste("Prior needed for", nm))

    allowed_dist <- available_priors[[spec[[nm]]]]

    if(!any(prior_dist_fam[allowed_dist] == prior[[paste0("prior_", nm, "_fam")]]))
      stop("Prior for ", nm, " must be one of the following: ",
           paste(allowed_dist, collapse = ", "))

    if(spec[[nm]] == "real_2"){
      if(length(prior[[paste0("prior_", nm, "_mean")]]) != 2)
        stop("Prior for ", nm, " must have the dimension equal to 2")
    }
  }
}
