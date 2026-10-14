
#' Effects of covariates on outcome in baggr models
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
#' @param summary logical; if `TRUE` returns summary statistic instead of all MCMC samples
#' @param transform a transformation (R function) to apply to the result;
#'                  (this is commonly used when calling from other
#'                  plotting or printing functions)
#' @param prob Probability mass of the central uncertainty interval, on the 0--1 scale.
#' @param interval Compatibility alias for `prob`; supply only one of these arguments.
#' @return With `summary = TRUE`, one row per covariate coefficient in the
#'   summary layout described here. Otherwise a matrix: columns are covariate coefficients and rows are draws from the posterior distribution.
#'         Number of rows depends on iterations in the MCMC (i.e. `x` in baggr(..., iter = x`)
#' @seealso [treatment_effect] for overall treatment effect across groups,
#'          [group_effects] for effects within each group,
#'          [effect_draw] and [effect_plot] for predicted treatment effect in new group
#'          (which you can condition on fixed effects using new data argument)
#' @export
#' @importFrom rstan extract


fixed_effects <- function(bg, summary = FALSE,
                             transform = NULL, prob = .95, interval = NULL) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
  check_if_baggr(bg)
  if(length(bg$covariates) == 0) {
    if(summary)
      return(summary_draws(array(0, c(1, 0, length(bg$effects))),
                           prob, bg$effects))
    return(array(0, dim = c(0,0,1)))
  }

  beta <- rstan::extract(bg$fit, pars="beta")[[1]]
  # colnames(beta) <- bg$covariates
  colnames(beta) <- attr(bg$inputs, "covariate_coding")

  if(!is.null(transform))
    beta <- do.call(transform, list(beta))

  if(summary) {
    beta <- array(beta, c(nrow(beta), ncol(beta), length(bg$effects)),
                  dimnames = list(NULL, colnames(beta), bg$effects))
    beta <- summary_draws(beta, prob, bg$effects)
  }
  return(beta)
}
