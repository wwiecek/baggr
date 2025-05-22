
#' Effects of covariates on outcome in baggr models
#'
#' @param bg a [baggr] model
#' @param summary logical; if `TRUE` returns summary statistic instead of all MCMC samples
#' @param transform a transformation (R function) to apply to the result;
#'                  (this is commonly used when calling from other
#'                  plotting or printing functions)
#' @param interval uncertainty interval width (numeric between 0 and 1), if `summary=TRUE`
#' @return A matrix: columns are covariate coefficients and rows are draws from the posterior distribution.
#'         Number of rows depends on iterations in the MCMC (i.e. `x` in baggr(..., iter = x`)
#' @seealso [treatment_effect] for overall treatment effect across groups,
#'          [group_effects] for effects within each group,
#'          [effect_draw] and [effect_plot] for predicted treatment effect in new group
#'          (which you can condition on fixed effects using new data argument)
#' @export
#' @importFrom rstan extract


fixed_effects <- function(bg, summary = FALSE,
                             transform = NULL, interval = .95) {
  check_if_baggr(bg)
  if(length(bg$covariates) == 0)
    return(array(0, dim = c(0,0,1)))

  beta <- rstan::extract(bg$fit, pars="beta")[[1]]
  # colnames(beta) <- bg$covariates
  colnames(beta) <- attr(bg$inputs, "covariate_coding")

  if(!is.null(transform))
    beta <- do.call(transform, list(beta))

  if(summary){
    beta <- t(apply(beta, 2, mint, int=interval, median=TRUE, sd = TRUE))
    colnames(beta) <- c("lci", "mean", "uci", "median", "sd")
    beta <- replicate(1, beta) #by the convention where 3rd dimension is effect
  }
  return(beta)
}
