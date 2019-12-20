#' Pooling metrics for baggr
#'
#' Compute heterogeneity/pooling statistics (for groups or whole model)
#' given a [baggr] meta-analysis model.
#' The statistics are the pooling metric by Gelman & Pardoe (2006) or its
#' complement, the _I-squared_ statistic.
#'
#' @param bg output of a baggr() function
#' @param type Calculation is done for each of the `"groups"`
#'            (default) or for `"total"` hypereffect(s).
#'             See details for how calculation is done.
#' @param summary logical; if `FALSE` a whole vector of pooling values is returned,
#'                otherwise only the means and intervals
#'
#' @details
#' Pooling statistic describes the extent to which group-level estimates of treatment
#' effect are "pooled" (or pulled!) toward average treatment effect in the meta-analysis model.
#' If `pooling = "none"` or "full" in [baggr], then the returned values are always 0 or 1, respectively.
#' If `pooling = "partial"`, the value is somewhere between 0 and 1.
#'
#' **Formulae for the calculations below are provided in main package vignette.** See `vignette("baggr").`
#'
#' #' **Estimate of pooling in a group**: this is the calculation done by `pooling()`
#' if `type = "groups"` (default).
#'
#' In a partial pooling model (see [baggr]), group _k_ (e.g. study) has a treatment effect
#' estimate, with some SE around the real treatment effect (TE).
#' Each TE itself is distributed with mean and variance.
#'
#' The quantity of interest is ratio of variability in \eqn{\tau} to total variability.
#' By convention, we subtract it from 1, to obtain a _pooling metric_ _p_.
#'
#' \deqn{p = 1 - (\sigma(\tau)^2 / (\sigma_(\tau)^2 + se_k^2))}
#'
#' * If \eqn{p < 0.5}, that means the variation across studies is higher than variation within studies.
#' * Values close to 1 indicate nearly full pooling. Variation across studies dominates.
#' * Values close to 0 -- no pooling. Variation within studies dominates.
#'
#' Note that, since \eqn{\sigma_{\tau}^2} is a Bayesian parameter (rather than a single fixed value)
#' _p_ is also a parameter. It is typical for _p_ to have very high dispersion, as in many cases we
#' cannot precisely estimate \eqn{\sigma_{\tau}}. To obtain the whole distribution of_p_
#' (rather than summarised values), set `summary=FALSE`.
#'
#' **Overall pooling (in the model)**
#'
#' Typically it is a single measure of heterogeneity that is of interest to researchers.
#' This is calculated by setting `type = "total"` or simply writing `heterogeneity(mymodel)`
#'
#' In many contexts, i.e. medical statistics, it is typical to report _1-P_, called \eqn{I^2}
#' (see Higgins _et al_, 2003). Higher values of _I-squared_ indicate higher heterogeneity.
#' Von Hippel (2015) provides useful details for _I-squared_ calulations.
#'
#' Same as for group-specific estimates, _P_ is a Bayesian parameter and its dispersion can be high.
#'
#'
#' **Relationship to R-squared statistic**
#'
#' See Gelman & Pardoe (2006) Section 1.1 for a short explanation of how \eqn{R^2}
#' statistic relates to the pooling metric.
#'
#'
#' @seealso
#' [metafor::rma.uni] for frequentist I^2 calculations of heterogeneity statistics
#'
#'
#' @references
#' Gelman, Andrew, and Iain Pardoe.
#' "Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models."
#' _Technometrics 48, no. 2 (May 2006): 241-51_. <https://doi.org/10.1198/004017005000000517>.
#'
#' Higgins, Julian P T, Simon G Thompson, Jonathan J Deeks, and Douglas G Altman.
#' "Measuring Inconsistency in Meta-Analyses."
#' _British Medical Journal 327, no. 7414 (September 6, 2003): 557-60._
#'
#' Hippel, Paul T von. "The Heterogeneity Statistic I2 Can Be Biased in Small Meta-Analyses."
#' _BMC Medical Research Methodology 15 (April 14, 2015)._ <https://doi.org/10.1186/s12874-015-0024-z>.
#'
#' @return Matrix with mean and intervals for chosen pooling metric,
#'         each row corresponding to one meta-analysis group.
#' @export
#'

pooling <- function(bg,
                    type = c("groups", "total"),
                    summary = TRUE) {
  type <- match.arg(type, c("groups", "total"))

  # we have to rig it for no pooling cases
  # because sigma_tau parameter might be meaningless then
  if(bg$pooling == "none")
    return(array(0, c(3, bg$n_groups, bg$n_parameters)))
  if(bg$pooling == "full")
    return(array(1, c(3, bg$n_groups, bg$n_parameters)))

  # we'll replace by switch() in the future
  if(bg$model %in% c("rubin", "mutau", "logit", "full")) {

    sigma_tau <- treatment_effect(bg)$sigma_tau

    # Grab the appropriate SE:
    if(bg$model == "mutau")
      sigma_k <- bg$data$se.tau
    if(bg$model == "rubin")
      sigma_k <- bg$data$se
    if(bg$model == "logit")
      sigma_k <- prepare_ma(bg$data, effect = "logOR")$se
    if(bg$model == "full")
      sigma_k <- group_effects(bg, summary = TRUE)[, "sd", 1]

    if(type == "groups")
      ret <- sapply(sigma_k, function(se) se^2 / (se^2 + sigma_tau^2))
    if(type == "total")
      ret <- replicate(1, mean(sigma_k^2) / (mean(sigma_k^2) + sigma_tau^2))
    ret <- replicate(1, ret) #third dim is always N parameters, by convention

  } else if(bg$model == "quantiles") {
    # compared to individual-level, here we are dealing with 1 more dimension
    # which is number of quantiles; so if sigma_tau above is sigma of trt effect
    # here it is N effects on N quantiles etc.

    # discard everything off-diagonal, we only care about variances here:
    sigma_tau <- treatment_effect(bg)$sigma_tau
    for(i in 1:dim(sigma_tau)[2])
      sigma_tau[,i,1] <- sigma_tau[,i,i]
    sigma_tau <- sigma_tau[,,1] #rows are samples, columns are quantiles

    # now SE's of study treatment effects: this is our input data!
    # Sigma_y_k_1 is now vcov of trt effect (used to be control + effect, mind)
    sigma_k <- t(apply(bg$inputs$Sigma_y_k_1, 1, diag)) #rows are studies, cols are quantiles

    # output is an array with dim = samples, studies, parameters
    # so let's convert to these dimensions:
    # for sigma_k the values are fixed over samples
    sigma_k <- array(rep(sigma_k, each = nrow(sigma_tau)),
                     dim = c(nrow(sigma_tau), nrow(sigma_k), ncol(sigma_k)))
    # for sigma_tau it doesn't change with studies, so:
    sigma_tau <- replicate(dim(sigma_k)[2], sigma_tau)
    # but this has (parameters, studies) rather than (studies, parameters) so
    sigma_tau <- aperm(sigma_tau, c(1, 3, 2))

    # now we can just do the operation on arrays and preserve dimensions:
    ret <- sigma_k / (sigma_k + sigma_tau)

    if(type == "total"){
      warning("Total pooling not implemented for quantiles model")
      ret <- 0
    }
  }

  if(summary)
    ret <- apply(ret, c(2,3), mint)

  return(ret)

}

#' @rdname pooling
#' @export
heterogeneity <- function(bg, summary = TRUE)
  pooling(bg, type = "total", summary = summary)
