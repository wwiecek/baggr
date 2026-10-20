
#' Average treatment effects in a baggr model
#'
#' The most general `treatment_effect` displays
#' both hypermean and hyperSD (as a list of length 2),
#' whereas `hypermean` and `hypersd` can be used as shorthands.
#' The `transform` argument applies to posterior draws of model parameters:
#' for example, `hypermean(bg, transform = exp)` exponentiates only the
#' posterior distribution of the mean treatment effect. This is useful when a
#' model was fit on a transformed scale, such as a log risk ratio or log odds
#' ratio scale. It is different from [effect_draw()], which first draws from
#' the posterior predictive distribution of a new treatment effect, including
#' heterogeneity or uncertainty around the mean, and then applies the
#' transformation. Because a transformed hyper-SD is not generally
#' interpretable, `treatment_effect()` returns `sigma = NA` when
#' `transform` is supplied; this is also what `hypersd(bg, transform = ...)`
#' returns.
#'
#' With `summary = TRUE`, summaries are numeric matrices with rows for
#' the quantities being summarised and columns `mean`, `median`, `sd`,
#' `lower`, `upper`, `width`, in that order. The bounds form a central
#' interval; `width` is its requested probability mass, even when estimates
#' are unavailable. For multiple effects, the result is a rows by statistics
#' by effects array, with effect names in the third dimension. A single-effect
#' matrix stores its effect name in the `effect` attribute. Row names are
#' preserved. Aggregate summaries have one row; predictions conditional on
#' `newdata` have one row per observation.
#'
#' @param bg a [baggr::baggr()] model
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param prob Probability mass of the central uncertainty interval, on the 0--1 scale.
#' @param interval Compatibility alias for `prob`; supply only one of these arguments.
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `treatment_effect` from other
#'                  plotting or printing functions)
#' @param rename_int logical; with `summary = TRUE`, rename only `lower`
#'   and `upper` to endpoint percentages and omit `width`. `median` retains
#'   its name. Defaults to `FALSE`.
#' @param message logical; use to disable messages prompted by using with
#'                no pooling models
#' @describeIn treatment_effect A list with 2 components (posterior draws)
#'         `theta` (mean effect) and `sigma` (SD). If `summary=TRUE`,
#'         each component follows the summary layout described here.
#' @export
#' @importFrom rstan extract
#' @examples
#' \donttest{
#' # Fit on a log scale, then exponentiate the hypermean draws to obtain ratios.
#' log_rr <- data.frame(yi = log(c(0.7, 0.9, 1.1, 1.3)),
#'                      sei = c(0.2, 0.25, 0.3, 0.35))
#' bg_log <- baggr(log_rr, pooling = "partial", iter = 500, refresh = 0)
#'
#' hypermean(bg_log, transform = exp)
#' }


treatment_effect <- function(bg, summary = FALSE,
                             transform = NULL, prob = .95,
                             message = TRUE, interval = NULL, rename_int = FALSE) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
  check_if_baggr(bg)

  if(bg$pooling == "none"){
    if(message)
      message("There is no treatment effect estimated when pooling = 'none'.")
    if(summary) {
      unavailable <- summary_draws(array(NA_real_, c(1, 1, length(bg$effects))),
                                    prob, bg$effects, rename_int = rename_int)
      return(list(theta = unavailable, sigma = unavailable))
    }
    return(list(theta = as.numeric(NA), sigma = as.numeric(NA)))
  }
  if(bg$model %in% c("rubin", "vivs", "vivs_full", "logit", "rubin_full")) {
    theta <- rstan::extract(bg$fit, pars="mu")[[1]]
    if(bg$model %in% c("vivs", "vivs_full"))
      theta <- theta[,1,2] #in this joint model the second element is trt effect, first is baseline
    if(is.matrix(theta) && ncol(theta) == 1) theta <- c(theta) #convert to vector for 1D cases

    if(bg$pooling == "partial"){
      if(bg$model %in% c("vivs", "vivs_full"))
        sigma <- rstan::extract(bg$fit, pars="hypersd[1,2]")[[1]]
      else
        sigma <- rstan::extract(bg$fit, pars="tau")[[1]]
      if(is.matrix(sigma) && ncol(sigma) == 1) sigma <- c(sigma)
    }
    if(bg$pooling == "full")
      sigma <- rep(0, length(bg$effects)) #same dim as theta, but by convention set to 0
  } else if(bg$model == "quantiles") {
    # In this case we have N columns = N quantiles
    theta <- as.matrix(bg$fit, "beta_1")
    # only take diagonals:
    sigma <- t(apply(rstan::extract(bg$fit, "Sigma_1")[[1]], 1, diag))
    # in model with correlation, we have Var(), not SD()
    sigma <- sqrt(sigma)
  } else if(bg$model == "sslab") {
    mean_params <- c("tau[1]", "tau[2]",
                     "sigma_TE[1]", "sigma_TE[2]",
                     "kappa[1,1,2]", "kappa[1,2,2]")
    sigma_params <- c("hypersd_tau[1]", "hypersd_tau[2]",
                      "hypersd_sigma_TE[1]", "hypersd_sigma_TE[2]",
                      "hypersd_kappa[1,1,2]", "hypersd_kappa[1,2,2]")
    theta <- as.matrix(bg$fit, mean_params)
    sigma <- as.matrix(bg$fit, sigma_params)

  } else {
    stop("Can't calculate treatment effect for this model.")
  }
  if(length(bg$effects) > 1){
    colnames(theta) <- bg$effects
    if(bg$pooling != "full") colnames(sigma) <- bg$effects
  }
  if(!is.null(transform)){
    theta <- do.call(transform, list(theta))
    sigma <- NA_real_ # by convention we set it to NA so that people don't transform
    # and then do operations on it by accident
  }
  if(summary) {
    theta <- summary_draws(theta, prob, bg$effects, aggregate = TRUE,
                          rename_int = rename_int)
    if(!is.null(transform) || bg$pooling == "full")
      sigma <- array(if(is.null(transform)) 0 else NA_real_,
                         c(1, 1, length(bg$effects)))
    sigma <- summary_draws(sigma, prob, bg$effects, aggregate = TRUE,
                                rename_int = rename_int)
  }

  return(list(theta = theta, sigma = sigma))
}



#' @describeIn treatment_effect The hypermean of a `baggr` model, shorthand for `treatment_effect(x, s=T)[[1]]`
#' @export
hypermean <- function(bg, transform = NULL, prob = 0.95, message = FALSE,
                      summary = TRUE, interval = NULL, rename_int = FALSE) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
  t <- treatment_effect(bg,summary=summary,transform=transform,prob=prob,message=FALSE,
                        rename_int=rename_int)
  if(message){
    cat(paste0("Hypermean of ",bg$effects," with ",prob*100,"% interval:\n"))
    print(t[[1]])
  } else
    t[[1]]
}

#' @describeIn treatment_effect The hyper-SD of a `baggr` model, shorthand for `treatment_effect(x, s=T)[[2]]`
#' @export
hypersd <- function(bg, transform = NULL, prob = 0.95, message = FALSE,
                    summary = TRUE, interval = NULL, rename_int = FALSE) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
  t <- treatment_effect(bg,summary=summary,transform=transform,prob=prob,message=FALSE,
                        rename_int=rename_int)
  if(message){
    cat(paste0("Hyper SD of ",bg$effects," with ",prob*100,"% interval:\n"))
    print(t[[2]])
  } else
    t[[2]]
}


#' Correlation between baseline means and treatment effects in a baggr model
#'
#' For varying-intercept, varying-slope models, `vivs_cor()` returns the
#' correlation between control means and treatment effects across studies.
#' `mutau_cor()` is retained as an alias.
#'
#' @param bg A varying-intercept, varying-slope model (`"vivs"` or `"vivs_full"`).
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param prob Probability mass of the central uncertainty interval, on the 0--1 scale.
#' @param interval Compatibility alias for `prob`; supply only one of these arguments.
#' @return Raw draws as a vector, or a one-row numeric summary matrix with
#'   columns `mean`, `median`, `sd`, `lower`, `upper`, `width` and an
#'   `effect` attribute of `"correlation"`. Correlation summaries are
#'   unavailable for models with no or full pooling; their `width` remains known.
#' @export
vivs_cor <- function(bg, summary = FALSE, prob = 0.95, interval = NULL) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
  if(summary && bg$pooling %in% c("none", "full"))
    return(summary_draws(NA_real_, prob, "correlation"))
  # m <- matrix(apply(as.matrix(bg$fit, "L_Omega"), 2, mean), 2, 2)
  # h <- apply(as.matrix(bg$fit, "hypersd"), 2, mean)
  # (m %*% t(m))[2,1]
  # diag(h) %*% (m %*% t(m)) %*% diag(h)

  m <- as.matrix(bg$fit, "L_Omega")
  # we want entry (2,1) in LL^T which is simply (1,1)*(2,1)
  # in a lower-triangular matrix; and (1,1) should be == 1 everywhere
  if(!all(m[,1] == 1))
    warning("Error with correlation matrix in the joint baseline/effect model. Inspect L_Omega parameters.")
  if(summary)
    return(summary_draws(m[,2], prob, "correlation"))
  else
    return(m[,2])

}

#' @rdname vivs_cor
#' @export
mutau_cor <- vivs_cor
