
#' Average treatment effect in a baggr model
#'
#' @param bg a [baggr] model
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param interval uncertainty interval width (numeric between 0 and 1), if summarising
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `treatment_effect` from other
#'                  plotting or printing functions)
#' @return A list with 2 vectors (corresponding to MCMC samples)
#'         `tau` (mean effect) and `sigma_tau` (SD). If `summary=TRUE`,
#'         both vectors are summarised as mean and lower/upper bounds according to
#'         `interval`
#' @export
#' @importFrom rstan extract


treatment_effect <- function(bg, summary = FALSE,
                             transform = NULL, interval = .95) {
  check_if_baggr(bg)

  if(bg$pooling == "none"){
    message("There is no treatment effect estimated when pooling = 'none'.")
    return(list(tau = as.numeric(NA), sigma_tau = as.numeric(NA)))
  }
  if(bg$model %in% c("rubin", "mutau", "mutau_full", "logit", "rubin_full")) {
    tau <- rstan::extract(bg$fit, pars="mu")[[1]]
    if(bg$model %in% c("rubin", "logit"))
      tau <- c(tau)
    if(bg$model %in% c("mutau", "mutau_full"))
      tau <- tau[,1,2]

    if(bg$pooling == "partial"){
      sigma_tau <- rstan::extract(bg$fit, pars="tau")[[1]]
      if(bg$model %in% c("rubin", "logit"))
        sigma_tau <- c(sigma_tau)
      if(bg$model %in% c("mutau", "mutau_full"))
        sigma_tau <- sqrt(sigma_tau[,1,2,2])
    }
    if(bg$pooling == "full")
      sigma_tau <- 0 #same dim as tau, but by convention set to 0
  } else if(bg$model == "quantiles") {
    # In this case we have N columns = N quantiles
    tau <- as.matrix(bg$fit, "beta_1")
    # only take diagonals:
    sigma_tau <- t(apply(rstan::extract(bg$fit, "Sigma_1")[[1]], 1, diag))
    # in model with correlation, we have Var(), not SD()
    sigma_tau <- sqrt(sigma_tau)
  } else if(bg$model == "sslab") {
    mean_params <- c("tau[1]", "tau[2]",
                     "sigma_TE[1]", "sigma_TE[2]",
                     "kappa[1,1,2]", "kappa[1,2,2]")
    sigma_params <- c("hypersd_tau[1]", "hypersd_tau[2]",
                      "hypersd_sigma_TE[1]", "hypersd_sigma_TE[2]",
                      "hypersd_kappa[1,1,2]", "hypersd_kappa[1,2,2]")
    tau <- as.matrix(bg$fit, mean_params)
    sigma_tau <- as.matrix(bg$fit, sigma_params)

  } else {
    stop("Can't calculate treatment effect for this model.")
  }

  if(length(bg$effects) > 1)
    colnames(tau) <- colnames(sigma_tau) <- bg$effects

  if(!is.null(transform)){
    tau <- do.call(transform, list(tau))
    sigma_tau <- NA # by convention we set it to NA so that people don't transform
    # and then do operations on it by accident
  }
  if(summary) {
    tau <- mint(tau, int=interval, median=TRUE, sd = TRUE)
    sigma_tau <- mint(sigma_tau, int=interval, median=TRUE, sd = TRUE)
  }

  return(list(tau = tau, sigma_tau = sigma_tau))
}



#' Make predictive draws from baggr model
#'
#' This function takes the samples of hyperparameters from a [baggr] model
#' (typically hypermean and hyper-SD, which you can see using [treatment_effect]) and draws values of
#' new realisations of treatment effect, i.e. an additional draw from the "population of studies".
#' This can be used for both prior and posterior draws, depending on [baggr] model.
#'
#' @param x A `baggr` class object.
#' @param transform a transformation (an R function) to apply to the result of a draw.
#' @param summary logical; if TRUE returns summary statistics rather than samples from the distribution;
#' @param interval uncertainty interval width (numeric between 0 and 1), if `summary=TRUE`
#' @param n How many values to draw? The default is as long as the number of samples
#'          in the `baggr` object (see _Details_).
#'
#' @return A vector (with `n` values) for models with one treatment effect parameter,
#'         a matrix (`n` rows and same number of columns as number of parameters) otherwise.
#' @export
#'
#' @seealso [treatment_effect] returns samples of hypermean and hyper-SD
#'          which are used by this function
#'
#' @details
#' The predictive distribution can be used to "combine" heterogeneity between treatment effects and
#' uncertainty in the mean treatment effect. This is useful both in understanding impact of
#' heterogeneity (see Riley et al, 2011, for a simple introduction) and for study design e.g.
#' as priors in analysis of future data (since the draws can be seen as an expected treatment effect
#' in a hypothetical study).
#'
#' The default number of samples is the same as what is returned by Stan model implemented in [baggr],
#' (depending on such options as `iter`, `chains`, `thin`). If `n` is larger than what is available
#' in Stan model, we draw values with replacement. This is not recommended and warning is printed in
#' these cases.
#'
#' Under default settings in [baggr], a _posterior_ predictive distribution is obtained. But
#' `effect_draw` can also be used for _prior_ predictive distributions when
#' setting `ppd=T` in [baggr]. The two outputs work exactly the same way.
#'
#' @references
#' Riley, Richard D., Julian P. T. Higgins, and Jonathan J. Deeks.
#' "Interpretation of Random Effects Meta-Analyses".
#' _BMJ 342 (10 February 2011)._ <https://doi.org/10.1136/bmj.d549>.
#'
effect_draw <- function(x, n, transform=NULL, summary = FALSE, interval = .95) {
  check_if_baggr(x)

  te <- treatment_effect(x)

  # Resize trt effects to the demanded size by making extra draws
  neffects <- length(x$effects)
  if(!missing(n)){
    if(neffects > 1){
      if(n > nrow(te$tau))
        warning("Making more effect draws than there are available samples in Stan object.",
                "Consider running baggr() with higher iter= setting")
      rows <- sample(nrow(te$tau), n, replace = TRUE)
      te$tau   <- te$tau[rows,]
      te$sigma_tau <- te$sigma_tau[rows,]
    }
    if(neffects == 1){
      if(n > length(te$tau))
        warning("Making more effect draws than there are available samples in Stan object.",
                "Consider running baggr() with higher iter= setting")
      rows <- sample(length(te$tau), n, replace = TRUE)
      te$tau   <- te$tau[rows]
      te$sigma_tau <- te$sigma_tau[rows]
    }
  }

  # Make draws using normal distribution (for now it's the only option)
  new_tau <- rnorm(length(te$tau), c(te$tau), c(te$sigma_tau))

  if(neffects > 1){
    new_tau <- matrix(new_tau, nrow(te$tau), ncol(te$tau))
    colnames(new_tau) <- colnames(te$tau)
  }
  if(!is.null(transform))
    new_tau <- do.call(transform, list(new_tau))

  if(summary) {
    new_tau <- mint(new_tau, int=interval, median=TRUE, sd = TRUE)
  }

  new_tau
}



#' Plot predictive draws from baggr model
#'
#' This function plots values from [effect_draw], the predictive distribution
#' (under default settings, _posterior_ predictive),
#' for one or more `baggr` objects.
#'
#' @param ... Object(s) of class [baggr]. If there is more than one,
#'            a comparison will be plotted and  names of objects
#'            will be used as a plot legend (see examples).
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `group_effects` from other
#'                  plotting or printing functions)
#' @return A `ggplot` object.
#' @import bayesplot
#' @export
#' @seealso [effect_draw] documents the process of drawing values;
#'          [baggr_compare] can be used as a shortcut for `effect_plot` with argument
#'          `compare = "effects"`
#'
#' @details
#' Under default settings in [baggr] posterior predictive is obtained. But
#' `effect_plot` can also be used for _prior_ predictive distributions when
#' setting `ppd=T` in [baggr]. The two outputs work exactly the same, but
#' labels will change to indicate this difference.
#'
#' @examples
#'
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
effect_plot <- function(..., transform=NULL) {
  l <- list(...)

  caption <- list(
    title = "Posterior distribution for possible treatment effect",
    subtitle = ""
  )
  if(!all(unlist(lapply(l, inherits, "baggr"))))
    stop("Effects plots can only be drawn for baggr class objects")
  if(all(unlist(lapply(l, attr, "ppd"))))
    caption <- list(
      title = "Prior distribution for possible treatment effect",
      subtitle = "No data used, only sampling from priors"
    )
  if(is.null(names(l))){
    if(length(l) > 1)
      message("Automatically naming models; please use named arguments to override.")
    names(l) <- paste("Model", 1:length(l))
  }

  # Check effects and prepare X label
  # if(any(unlist(lapply(l, function(x) length(x$effects))) > 1))
  # stop("Effect_plot is only possible for models with 1-dimensional treatment effects")
  # effects <- paste("Treatment effect on",
  # unique(unlist(lapply(l, function(x) x$effects))))
  effects <- unique(unlist(lapply(l, function(x) x$effects)))
  n_parameters <- unique(unlist(lapply(l, function(x) x$n_parameters)))
  if(length(n_parameters) != 1)
    stop("All models must have the same number of parameters")
  if(length(effects) > n_parameters)
    stop("All models must have same effects")



  # dim <-
  l <- lapply(l, effect_draw, transform=transform)

  df <- data.frame()
  for(i in seq_along(l)){
    if(n_parameters == 1)
      df <- rbind(df, data.frame("model"=names(l)[i],
                                 "value" = l[[i]]))
    else{
      # Melt with base R
      for(nm in colnames(l[[i]]))
        df <- rbind(df, data.frame("model"=names(l)[i],
                                   "variable" = nm,
                                   "value" = l[[i]][,nm]))
    }
  }

  single_model_flag <- (length(l) == 1)
  model <- value <- NULL
  ggplot(df, aes(value, group = model, fill = model)) +
    baggr_theme_get() +
    geom_density(alpha = .25) +
    ggtitle(label = caption$title,
            subtitle = caption$subtitle) +
    {if(n_parameters == 1) xlab(effects)} +
    {if(single_model_flag) theme(legend.position = "none")} +
    {if(n_parameters > 1) facet_wrap(~variable, ncol = 3, scales = "free")}
}
