
#' Average treatment effect in a baggr model
#'
#' @param bg a [baggr] model
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `group_effects` from other
#'                  plotting or printing functions)
#' @return A list with 2 vectors (corresponding to MCMC samples)
#'         `tau` (mean effect) and `sigma_tau` (SD)
#' @export
#' @importFrom rstan extract


treatment_effect <- function(bg, transform = NULL) {
  check_if_baggr(bg)

  if(bg$pooling == "none"){
    message("There is no treatment effect estimated when pooling = 'none'.")
    return(list(tau = as.numeric(NA), sigma_tau = as.numeric(NA)))
  }
  if(bg$model %in% c("rubin", "mutau", "logit")) {
    tau <- rstan::extract(bg$fit, pars="tau")[[1]]
    if(bg$model %in% c("rubin", "logit"))
      tau <- c(tau)
    if(bg$model == "mutau")
      tau <- tau[,1,2]

    if(bg$pooling == "partial"){
      sigma_tau <- rstan::extract(bg$fit, pars="sigma_tau")[[1]]
      if(bg$model %in% c("rubin", "logit"))
        sigma_tau <- c(sigma_tau)
      if(bg$model == "mutau")
        sigma_tau <- sqrt(sigma_tau[,1,2,2])
    }
    if(bg$pooling == "full")
      sigma_tau <- 0 #same dim as tau, but by convention set to 0
  } else if(bg$model == "full") {
    tau <- as.matrix(bg$fit)[,"mutau[2]"]
    sigma_tau <- as.matrix(bg$fit)[,"sigma_mutau[2,2]"]
    # in model with correlation, we have Var(), not SD():
    sigma_tau <- sqrt(sigma_tau)
  } else if(bg$model == "quantiles") {
    # In this case we have N columns = N quantiles
    tau <- as.matrix(bg$fit, "beta_1")
    # only take diagonals:
    sigma_tau <- t(apply(rstan::extract(bg$fit, "Sigma_1")[[1]], 1, diag))
    # in model with correlation, we have Var(), not SD()
    sigma_tau <- sqrt(sigma_tau)
  }

  if(!is.null(transform)){
    tau <- do.call(transform, list(tau))
    sigma_tau <- NA # by convention we set it to NA so that people don't convert
                    # and then do operations on it by accident
  }
  return(list(tau = tau, sigma_tau = sigma_tau))
}



#' Make posterior draws for treatment effect
#'
#' This function takes the samples of hyperparameters of a `baggr` model
#' (commonly hypermean tau and hyper-SD sigma_tau) and simulates values of
#' new realisations of tau (a mean effect in some unobserved group).
#'
#' @param x A `baggr` class object.
#' @param n How many values to draw? The default is the same
#'          as number of samples in the model (default is 2,000).
#' @return A vector of possible values of the treatment effect.
#' @export
#'
effect_draw <- function(x, n) {
  check_if_baggr(x)
  # Draw sigma and tau
  te <- do.call(cbind, treatment_effect(x))
  if(!missing(n))
    te <- te[sample(nrow(te), n, replace = T),]
  new_tau <- apply(te, 1, function(x) {
    if(any(is.na(x)))
      return(NA)
    else
      return(rnorm(1, x[1], x[2]))
  })
  new_tau
}

#' Plot posterior distribution for treatment effect
#'
#' This function plots the [effect_draw] for one or more baggr objects.
#'
#' @param ... Object(s) of class `baggr`. If there is more than one,
#'            the names of objects will be used as a plot legend (see example).
#' @return A ggplot.
#' @import bayesplot
#' @export
#' @seealso [baggr_compare] can be used as a shortcut for `effect_plot` with argument
#'          `compare = "effects"`
#' @examples
#'
#' # A single effects plot
#' bg1 <- baggr(schools, prior_hypersd = uniform(0, 20))
#' effect_plot(bg1)
#'
#' # Compare how posterior depends on the prior choice
#' bg2 <- baggr(schools, prior_hypersd = normal(0, 5))
#' effect_plot("Uniform prior on SD"=bg1,
#'             "Normal prior on SD"=bg2)
#'
#' # Compare the priors themselves (ppd=T)
#' bg1_ppd <- baggr(schools, prior_hypersd = uniform(0, 20), ppd=TRUE)
#' bg2_ppd <- baggr(schools, prior_hypersd = normal(0, 5), ppd=TRUE)
#' effect_plot("Uniform prior on SD"=bg1_ppd,
#'             "Normal prior on SD"=bg2_ppd)
#'
effect_plot <- function(...) {
  l <- list(...)

  caption <- "Possible treatment effect"
  if(!all(unlist(lapply(l, inherits, "baggr"))))
    stop("Effects plots can only be drawn for baggr class objects")
  if(all(unlist(lapply(l, attr, "ppd"))))
    caption <- "Possible treatment effect (prior predictive)"
  if(is.null(names(l))){
    if(length(l) > 1)
      message("Automatically naming models; please use named arguments to override.")
    names(l) <- paste("Model", 1:length(l))
  }

  # Check effects and prepare X label
  if(any(unlist(lapply(l, function(x) length(x$effects))) > 1))
    stop("Effect_plot is only possible for models with 1-dimensional treatment effects")
  effects <- paste("Effect on", unique(unlist(lapply(l, function(x) x$effects))))
  if(length(effects) > 1)
    stop("All models must have same effects")

  l <- lapply(l, effect_draw)
  df <- data.frame()
  for(i in seq_along(l))
    df <- rbind(df, data.frame("model"=names(l)[i],
                               "value" = l[[i]]))
  single_model_flag <- (length(l) == 1)
  model <- value <- NULL
  ggplot(df, aes(value, group = model, fill = model)) +
    baggr_theme_get() +
    geom_density(alpha = .25) +
    ggtitle(caption) +
    xlab(effects) +
    {if(single_model_flag) theme(legend.position = "none")}
}
