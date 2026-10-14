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
#' @param bg baggr object
#' @param summary logical; if `TRUE` returns summary statistics as explained below.
#' @param prob Probability mass of the central uncertainty interval, on the 0--1 scale.
#' @param interval Compatibility alias for `prob`; supply only one of these arguments.
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `group_effects` from other
#'                  plotting or printing functions)
#' @param random_only logical; for meta-regression models, should [fixed_effects] be included in the
#'                    returned group effect?
#' @param rename_int logical; with `summary = TRUE`, rename only `lower`
#'   and `upper` to endpoint percentages and omit `width`. `median` retains
#'   its name. Defaults to `FALSE`.
#' @return Either an array with MCMC samples (if `summary = FALSE`)
#'         or a summary of these samples (if `summary = TRUE`).
#'         For arrays the three dimensions are: N samples, N groups and N effects
#'         (equal to 1 for the basic models).
#' @examples
#' fit1 <- baggr(schools)
#' group_effects(fit1, summary = TRUE, prob = 0.5)
#' @seealso [fixed_effects] for effects of covariates on outcome. To extract random effects
#'          when covariates are present, you can use either [random_effects] or, equivalently,
#'          `group_effects(random_only=TRUE)`.
#'
#' @export
group_effects <- function(bg, summary = FALSE, transform = NULL, prob = .95,
                          random_only = FALSE,
                          rename_int = FALSE, interval = NULL) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
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
        m <- array(m[,1,2,,drop = FALSE], c(dim(m)[1], dim(m)[4], 1))
      # Some full-data models keep a singleton effects dimension.
      # Reduce to 2D so fixed-effect draws can be added by group.
      if(length(dim(m)) == 3 && dim(m)[3] == 1)
        m <- matrix(m, nrow = dim(m)[1], ncol = dim(m)[2])

      # If dealing with a meta-regression model, we automatically add effect of covariates
      # unless user requests random_only
      if(!random_only) {
        x_group <- group_effects_covariate_matrix(bg, n_groups = ncol(m))
        if(!is.null(x_group) && ncol(x_group) > 0) {
          m_fe <- fixed_effects(bg) %*% t(x_group)
          m <- m + m_fe
        }
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

  if(!is.null(transform))
    m <- do.call(transform, list(m))

  # will summarise if requested:
  if(summary) {
    m <- summary_draws(m, prob, effect_names, rename_int = rename_int)
  }

  return(m)
}


# Group-level design matrix for adding fixed effects to random group effects
# in meta-regression models.
group_effects_covariate_matrix <- function(bg, n_groups) {
  if(length(bg$covariates) == 0)
    return(NULL)

  covariate_coding <- attr(bg$inputs, "covariate_coding")
  if(is.null(covariate_coding) || length(covariate_coding) == 0)
    return(NULL)

  # Summary-data models already have one row per group.
  if(bg$model == "rubin") {
    x_group <- bg$inputs$X
    if(is.null(dim(x_group)))
      return(NULL)
    if(nrow(x_group) != n_groups)
      stop("Unable to align covariate matrix and group effects.")
    return(x_group)
  }

  if(is.null(bg$summary_data))
    return(NULL)

  # For individual-data models, summary_data includes only meta-regression
  # covariates that are fixed within studies.
  mr_covariates <- intersect(bg$covariates, names(bg$summary_data))
  if(length(mr_covariates) == 0)
    return(NULL)

  mm <- model.matrix(
    as.formula(paste("~", paste(mr_covariates, collapse = "+"))),
    data = bg$summary_data[, mr_covariates, drop = FALSE]
  )

  x_group <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  x_group <- x_group[, intersect(covariate_coding, colnames(x_group)), drop = FALSE]

  if(length(covariate_coding) > 0) {
    missing_cols <- setdiff(covariate_coding, colnames(x_group))
    if(length(missing_cols) > 0) {
      zeros <- matrix(0, nrow = nrow(x_group), ncol = length(missing_cols))
      colnames(zeros) <- missing_cols
      x_group <- cbind(x_group, zeros)
    }
    x_group <- x_group[, covariate_coding, drop = FALSE]
  }

  if(nrow(x_group) != n_groups)
    stop("Unable to align group-level covariates and group effects.")

  x_group
}


# Grabbing labels for group/study names
group_names <- function(bg) {
  if(is.null(attr(bg$inputs, "group_label")))
    par_names <- paste0("Group ", 1:attr(bg$inputs, "n_groups"))
  else
    par_names <- attr(bg$inputs, "group_label")
  par_names
}

#' Extract only random (treatment) effects from a baggr model
#'
#' This function is a shortcut for `group_effects(random_only=TRUE, ...)`. Note that this
#' is different to cluster random effects in individual-level data: by random effects we
#' mean the random component of the group-wide effect
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
