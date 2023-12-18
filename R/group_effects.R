#' Extract baggr study/group effects
#'
#' Given a baggr object, returns the raw MCMC draws of the posterior for
#' each group's effect or a summary of these draws. (We use "group" and "study" interchangeably.)
#' If there are no covariates in the model, this effect is a single random variable.
#' If there are covariates, the group effect is a sum of effect of covariates (fixed effects)
#' and the study-specific random variable (random effects).
#' This is an internal function currently used as a helper for plotting and
#' printing of results.
#'
#' @param bg baggr object
#' @param summary logical; if `TRUE` returns summary statistics as explained below.
#' @param interval uncertainty interval width (numeric between 0 and 1), if summarising
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `group_effects` from other
#'                  plotting or printing functions)
#' @param random_only logical; for meta-regression models, should [fixed_effects] be included in the
#'                    returned group effect?
#' @param rename_int logical; if `TRUE` then rather than returning `median`, `lci` and `uci`
#'                   columns they are renamed to e.g. `50%`, `2.5%`, `97.5%`; this only
#'                   works if `summary=TRUE`
#' @return Either an array with MCMC samples (if `summary = FALSE`)
#'         or a summary of these samples (if `summary = TRUE`).
#'         For arrays the three dimensions are: N samples, N groups and N effects
#'         (equal to 1 for the basic models).
#' @examples
#' fit1 <- baggr(schools)
#' group_effects(fit1, summary = TRUE, interval = 0.5)
#' @details If `summary = TRUE`, the returned object contains, for each study
#' or group, the following 5 values:
#' the posterior medians, the lower and upper bounds of the
#' uncertainty intervals using the central posterior credible interval
#' of width specified in the argument `interval`, the posterior mean, and
#' the posterior standard deviation.
#'
#' @seealso [fixed_effects] for effects of covariates on outcome. To extract random effects
#'          when covariates are present, you can use either [random_effects] or, equivalently,
#'          `group_effects(random_only=TRUE)`.
#'
#' @export
group_effects <- function(bg, summary = FALSE, transform = NULL, interval = .95,
                          random_only = FALSE,
                          rename_int = FALSE) {
  check_if_baggr(bg)

  # Grab group labels
  par_names <- group_names(bg)

  # Grab effect names
  effect_names <- bg$effects



  # m <- as.matrix(bg$fit)
  if(attr(bg , "ppd"))
    stop("There are no group effects in prior predictive distribution baggr objects.")

  if(bg$pooling == "full"){
    tau <- treatment_effect(bg)[["tau"]]
    k <- attr(bg$inputs, "n_groups")
    m <- replicate(k, tau)
    if(length(dim(m)) == 3)
      m <- aperm(m, c(1, 3, 2))

  } else {
    # choose correct columns for the given models:
    if(bg$model %in% c("rubin", "mutau", "mutau_full", "logit", "rubin_full")) {
      #replace by extract:
      # m <- m[, grepl("^tau_k", colnames(m))]
      m <- rstan::extract(bg$fit, pars = "theta_k")[[1]]
      # drop mu if model has mu (baseline/control value)
      if(bg$model %in% c("mutau", "mutau_full"))
        m <- m[,1,2,]

      # If dealing with a meta-regression model, we automatically add effect of covariates
      # unless user requests random_only
      if(bg$model == "rubin" && !random_only && !is.null(bg$covariates)) {
        m_fe <- fixed_effects(bg) %*% t(bg$inputs$X)
        m <- m + m_fe
      }

    } else if(bg$model == "quantiles") {
      # In this case we have 3D array, last dim is quantiles
      m <- rstan::extract(bg$fit, pars = "beta_1_k")[[1]]
    } else if(bg$model == "sslab") {
      m1 <- rstan::extract(bg$fit, "tau_k")[[1]]
      m2 <- rstan::extract(bg$fit, "sigma_TE_k")[[1]]
      m3 <- rstan::extract(bg$fit, "kappa_k")[[1]][,,1:2,2]
      m <- array(c(m1, m2, m3), c(nrow(m1), ncol(m1), length(effect_names)))
    } else {
      stop("Can't calculate treatment effect for this model.")
    }
  }
  # for consistency with quantiles, except we have 1 parameter only
  if(length(dim(m)) == 2)
    m <- array(m, dim = c(dim(m), 1))

  # Assing correct dimnames
  dimnames(m)[[2]] <- par_names
  dimnames(m)[[3]] <- effect_names

  # will summarise if requested:
  if(summary) {
    intval <- c((1-interval)/2, .5, 1 - (1-interval)/2)
    m <- apply(m, c(2,3), function(x) c(quantile(x, intval), mean(x), sd(x)))
    if(is.null(dimnames(m)[[2]]))
      dimnames(m)[[2]] <- 1:nrow(m)
    if(rename_int){
      dimnames(m)[[1]] <- c(paste0(as.character(100*intval), "%"), "mean", "sd")
    }else
      dimnames(m)[[1]] <- c("lci", "median", "uci", "mean", "sd")
    m <- aperm(m, c(2,1,3))
  }

  if(!is.null(transform))
    m <- do.call(transform, list(m))

  return(m)
}

#' Extract only random effects from a baggr model
#'
#' This function is a shortcut for `group_effects(random_only=TRUE, ...)`
#'
#' @param ... arguments passed to [group_effects]
#' @export
random_effects <- function(...) {
  group_effects(random_only = TRUE, ...)
}

#' `study_effects` is just an alias for `group_effects`
#' @rdname group_effects
#' @export
study_effects <- group_effects
