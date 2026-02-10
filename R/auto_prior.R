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
#' @param selection  same as in [baggr]
#' @param silent same as in [baggr]
#'
#' @return A named list with prior values that can be appended to `stan_data`
#'         and passed to a Stan model.
#'

prepare_prior <- function(prior,
                          data,
                          stan_data,
                          model,
                          pooling = c("partial", "none", "full"),
                          covariates = c(),
                          selection = NULL,
                          silent = FALSE) {

  if(missing(prior))
    prior <- list()

  # lists for saving the results of this procedure
  prior_list <- list() # List of priors formatted for Stan input (including some superfluous items)
  prior_dist_list <- list() # List of priors that were set in default format e.g. "normal(0,2.5)"
  # prior_dist_list used for presentation/readibility; prior_list used for running models

  pooling <- match.arg(pooling)

  # Determine number of random effects in the model, which will determine dimension of some priors
  if(model %in% c("rubin_full", "logit")) re_dim <- attr(stan_data, "n_re") else re_dim <- 1

  # In models that use individual-level data, swap `data` argument for summary data
  if(model %in% c("rubin_full", "logit", "mutau_full")) {
    trt_col <- attr(stan_data, "columns")[["treatment"]]
    pma_data <- data.frame(outcome = stan_data$y,
                           group = stan_data$site,
                           treatment = attr(stan_data, "data")[[trt_col]]
    )

    # if(model %in% c("rubin_full", "logit") && re_dim > 1)
    # stop("For a model with multiple interventions you must set prior manually.")
    if(model == "logit"){
      data <- prepare_ma(pma_data, effect = "logOR", rare_event_correction = 0.1)
      # Must add mu for automatic definition of baseline prior
      p <- data$c/data$n2
      data$mu <- log(p/(1-p))
    }
    if(model %in% c("mutau_full", "rubin_full")) {
      data <- prepare_ma(pma_data, effect = "mean")
      # mu and tau columns could be NA because in some studies you may not
      # have all treatments present; for purpose of deriving prior, we need
      # to ignore them:
      data <- data[!is.na(data$mu) & !is.na(data$tau),]
    }
  }
  # ...then proceed as you would in Rubin model

  if(!(model %in% c("rubin", "mutau", "logit", "rubin_full", "sslab", "mutau_full")))
    stop(paste0("Cannot prepare priors for model of type ", model))



  # Now set the priors in several steps:
  # (1) Setting model-specific priors (and checking they conform to spec)

  # For each model we need a list of what priors can be specified and
  # what families of distributions are allowed for them
  # (this second part should be changed to just specifying dimensionality/type)
  # This must match names.R
  priors_spec <- list(
    "rubin" = c("hypermean" = "real", "hypersd" = "positive_real"),
    "rubin_full"  = c("hypermean" = "real", "hypersd" = "positive_real",
                      "control" = "real", "control_sd" = "positive_real",
                      "sigma" = "positive_real"),
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
                 "scale_control"  = "real", "scale_control_sd"  = "positive_real",
                 "scale"  = "real", "scale_sd"  = "positive_real",
                 "kappa"  = "real", "kappa_sd"  = "positive_real")
  )

  # Now loop over all priors that need to be specified for a given model
  for(current_prior in names(priors_spec[[model]])) {
    if(is.null(prior[[current_prior]])) {
      if(model != "sslab"){
        default_prior_dist <- switch(
          current_prior,
          "hypermean"  =
            if(model %in% c("mutau", "mutau_full"))
              multinormal(c(0,0),
                          c(100*max(abs(data$mu)),
                            100*max(abs(data$tau)))*diag(2))
          else
            normal(0, 10*max(abs(data$tau))),
          # else
          # replicate(re_dim, normal(0, 10*max(abs(data$tau))), simplify = FALSE),
          "sigma"      = uniform(0, 10*max(c(sqrt(data$n.mu)*data$se.mu,
                                             sqrt(data$n.tau*data$se.tau)))),
          "hypersd"    = if(attr(stan_data, "n_groups") == 1)
            normal(0, 1) #this is not used because stop() will trigger below
          else
            uniform(0, 10*sd(data$tau)),
          "hypercor"   = lkj(3),
          "control"    = normal(0, 10*max(abs(data$mu))),
          "control_sd" = if(attr(stan_data, "n_groups") == 1)
            normal(0, 1) #this is not used because stop() will trigger below
          else
            uniform(0, 10*sd(data$mu))
        )
      }
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

      dist_to_set <- default_prior_dist

    } else {
      dist_to_set <- prior[[current_prior]]
    }


    convert_to_array <-
      model %in% c("logit", "rubin_full") &&
      current_prior %in% c("hypermean", "hypersd")

    # SET THE PRIOR VALUE HERE
    prior_list <- set_prior_val(prior_list,
                                paste0("prior_", current_prior),
                                dist_to_set,
                                p = if(convert_to_array) re_dim else 1,
                                to_array = convert_to_array)
    # Also save it in the "dist" format for future reference
    prior_dist_list[[current_prior]] <- dist_to_set

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
          if(attr(stan_data, "n_groups") == 1)
            stop("You must specify hyper-SD prior manually when data has only one row.")
          if(attr(stan_data, "n_groups") < 5)
            message(paste("/Dataset has only", attr(stan_data, "n_groups"),
                          "groups -- consider setting variance prior manually./"))
          if(model != "sslab"){
            message("Setting hyper-SD prior using 10 times the naive SD across sites")
            special_name <- "sigma_tau"
          }
        } else if(current_prior == "control") {
          if(model == "logit"){
            prop_ctrl <- data$c / (data$c + data$d)
            if(max(prop_ctrl) > .999 | min(prop_ctrl) < .001)
              message("Baseline proportion of events is very close to 0 or 1.",
                      "Consider manually setting prior_control.")
            special_name <- "log odds of event rate in untreated: mean"
          }
          if(model %in% c("rubin_full", "mutau_full")){
            if(stan_data$pooling_baseline == 1)
              special_name <- "mu (hyperparameter)"
            else if(stan_data$pooling_baseline == 0)
              special_name <- "independent prior on control means"
            else
              special_name <- "DNP"
          }
        } else if(current_prior == "control_sd") {
          if((attr(stan_data, "n_groups") == 1) && (stan_data$pooling_baseline == 1))
            stop("You must specify SD in baseline rates (control_sd) prior ",
                 "manually when data has only one row.")

          if(stan_data$pooling_baseline == 1){
            if(model == "logit")
              special_name <- "log odds of event rate in untreated: sd"
            if(model %in% c("rubin_full", "mutau_full")){
              special_name <- "sigma_mu (hyperparameter)"
            }
          }else
            special_name <- "DNP" #do not print, this is not used!
        } else if(current_prior == "sigma") {
          special_name <- "error term in linear regresion"
        }

        # 2) Print the prior:

        if(special_name == "")
          message(paste0("* ", current_prior, " ~ ",
                         print_dist(default_prior_dist)))
        else if(special_name != "DNP")
          message(paste0("* ", current_prior, " [", special_name, "] ~ ",
                         print_dist(default_prior_dist)))

      } else {
        if(current_prior == "hypersd" &&
           pooling != "partial")
          message("Prior for hyper-SD set, but pooling is not partial. Ignoring.")
        if(current_prior == "control_sd" &&
           model %in% c("logit", "rubin_full", "mutau_full") &&
           stan_data$pooling_baseline != 1)
          message("SD hyperparameter for control groups defined,",
                  "but there is no pooling. Ignoring it.")
      }
    }
  }

  # Check eligibility
  check_eligible_priors(prior_list, priors_spec[[model]])

  # (2) Setting fixed effects/beta/covariates prior
  if(length(covariates) > 0) {
    if(is.null(prior$beta)){
      val <- 10*max(apply(stan_data$X, 2, sd))
      prior$beta <- normal(0, val)
      prior_list <- set_prior_val(prior_list, "prior_beta", normal(0, val))
      message(
        paste0("Setting prior for covariates in regression to normal,",
               " with SD equal to 10*(highest SD among covariates):\n",
               "* beta ~ ", print_dist(normal(0, val)),
               "---purely experimental, we recommend you set this manually"))
    }
    prior_dist_list[["beta"]] <- prior$beta
  } else {
    prior$beta <- uniform(0,1)
  }
  prior_list <- set_prior_val(prior_list, "prior_beta", prior$beta)

  # (3) Setting clustering prior
  if(!is.null(stan_data$clustered) && (stan_data$clustered == 1)) {
    if(is.null(prior$cluster)){
      val <- 10*sd(abs(tapply(stan_data$y, stan_data$cluster, mean)))
      prior$cluster <- normal(0, val)
      message(
        paste0("Setting prior for cluster SD to normal,",
               " with SD equal to 10*(SD in all cluster-level means):\n",
               "* cluster ~ ", print_dist(normal(0, val)),
               "---purely experimental, we recommend you set this manually"))
    }
    prior_dist_list[["cluster"]] <- prior$cluster
  } else {
    prior$cluster <- uniform(0,1)
  }
  prior_list <- set_prior_val(prior_list, "prior_cluster", prior$cluster)

  # (4) Setting the selection prior (for relative publication Pr on log scale)
  if(model == "rubin" && !is.null(selection)){
    if(is.null(prior$selection)){
      prior$selection <- normal(0, 2)
      message("* log(relative publication probability) ~ N(0, 2); consider setting this manually.")
    }
    prior_dist_list[["selection"]] <- prior$selection
  } else {
    prior$selection <- uniform(0, 1)
  }
  prior_list <- set_prior_val(prior_list, "prior_sel", prior$selection)

  return(list(
    stan = prior_list,
    dist = prior_dist_list))
}



check_eligible_priors <- function(prior, spec) {
  for(nm in names(spec)){
    if(!(paste0("prior_", nm, "_fam") %in% names(prior)))
      stop(paste("Prior needed for", nm))

    allowed_dist <- available_priors[[spec[[nm]]]]

    # Usually priors family is length 1 but sometimes it can be a vector, so loop:
    for(cud in prior[[paste0("prior_", nm, "_fam")]]){
      if(!any(prior_dist_fam[allowed_dist] == cud))
        stop("Prior for ", nm, " must be one of the following: ",
             paste(allowed_dist, collapse = ", "))
    }
    if(spec[[nm]] == "real_2"){
      if(length(prior[[paste0("prior_", nm, "_mean")]]) != 2)
        stop("Prior for ", nm, " must have the dimension equal to 2")
    }
  }
}
