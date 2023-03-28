#' Pooling metrics and related statistics for baggr
#'
#' Compute statistics relating to
#' `pooling` in a given [baggr] meta-analysis model returns statistics, for
#' either the entire model or individual groups, such as
#' pooling statistic by Gelman & Pardoe (2006), _I-squared_, _H-squared_, or study weights;
#' `heterogeneity` is a shorthand for `pooling(type = "total")`
#' `weights` is shorthand for `pooling(metric = "weights")`
#'
#' @param bg a [baggr] model
#' @param metric `"pooling"` for Gelman & Pardoe statistic _P_,
#'               `"isq"` for I-squared statistic (_1-P_, Higgins & Thompson, 2002)
#'               `"hsq"` for H squared statistic (_1/P_, ibid.);
#'               `"weights"` for study weights;
#'               also see _Details_
#' @param type In `pooling` calculation is done for each of the `"groups"`
#'            (default) or for `"total"` hypereffect(s).
#' @param summary logical; if `FALSE` a whole vector of pooling values is returned,
#'                otherwise only the means and intervals
#'
#' @details
#'
#' Pooling statistic (Gelman & Pardoe, 2006) describes the extent to which
#' group-level estimates of treatment
#' effect are "pooled" toward average treatment effect in the meta-analysis model.
#' If `pooling = "none"` or `"full"` (which you specify when calling [baggr]),
#' then the values are always 0 or 1, respectively.
#' If `pooling = "partial"`, the value is somewhere between 0 and 1.
#' We can distinguish between pooling of individual groups and overall pooling in
#' the model.
#'
#' In many contexts, i.e. medical statistics, it is typical to report _1-P_, called \eqn{I^2}
#' (see Higgins and Thompson, 2002; sometimes another statistic, \eqn{H^2 = 1 / P},
#' is used).
#' Higher values of _I-squared_ indicate higher heterogeneity;
#' Von Hippel (2015) provides useful details for _I-squared_ calculations (and some
#' issues related to it, especially in frequentist models).
#' See Gelman & Pardoe (2006) Section 1.1 for a short explanation of how \eqn{R^2}
#' statistic relates to the pooling metric.
#'
#' @section Group pooling:
#'
#' This is the calculation done by `pooling()` if `type = "groups"` (default).
#' In a partial pooling model (see [baggr] and above), group _k_ (e.g. study) has
#' standard error of treatment effect estimate, \eqn{se_k}.
#' The treatment effect (across _k_ groups) is variable across groups, with
#' hyper-SD parameter \eqn{\sigma_(\tau)}.
#'
#' The quantity of interest is ratio of variation in treatment effects to the
#' total variation.
#' By convention, we subtract it from 1, to obtain a _pooling metric_ _P_.
#'
#' \deqn{p = 1 - (\sigma_(\tau)^2 / (\sigma_(\tau)^2 + se_k^2))}
#'
#' * If \eqn{p < 0.5}, the variation across studies is higher than variation within studies.
#' * Values close to 1 indicate nearly full pooling. Variation across studies dominates.
#' * Values close to 0 indicate no pooling. Variation within studies dominates.
#'
#' Note that, since \eqn{\sigma_{\tau}^2} is a Bayesian parameter (rather than a
#' single fixed value),
#' _P_ is also a parameter. It is typical for _P_ to have very high dispersion,
#' as in many cases we
#' cannot precisely estimate \eqn{\sigma_{\tau}}. To obtain samples from the distribution
#' of _P_ (rather than summarised values), set `summary=FALSE`.
#'
#'
#' @section Study weights:
#'
#' Contributions of each group (e.g. each study) to the mean meta-analysis estimate
#' can be calculated by calculating for each study *w_k* the inverse of sum of group-specific
#' SE squared and between-study variation.
#' To obtain weights, this vector (across all studies) has to be normalised to 1, i.e.
#' *w_k/sum(w_k)* for each _k_.
#'
#' SE is typically treated as a fixed quantity
#' (and usually reported on the reported point estimate),
#' but between-study variance is a model parameter,
#' hence the weights themselves are also random variables.
#'
#'
#'
#'
#' @section Overall pooling in the model:
#'
#' Typically researchers want to report a single measure from the model,
#' relating to heterogeneity across groups.
#' This is calculated by either `pooling(mymodel, type = "total")` or simply
#' `heterogeneity(mymodel)`
#'
#' Formulae for the calculations below are provided in main package vignette and
#' almost analogous to the group calculation above, but using mean variance across
#' all studies. In other words, pooling _P_ is simply ratio of the expected within-study
#' variance term to total variance.
#'
#' The typical study variance is calculated following Eqn. (1) and (9)
#' in Higgins and Thompson (see References). We use this formulation
#' to make our pooling and I^2 comparable with other meta-analysis implementations,
#' but users should be aware that this is only one possibility for calculating
#' that "typical" within-study variance.
#'
#' Same as for group-specific estimates, _P_ is a Bayesian parameter and its
#' dispersion can be high.
#'
#'
#' @section Value:
#'
#' Matrix with mean and intervals for chosen pooling metric,
#' each row corresponding to one meta-analysis group.
#'
#' @references
#' Gelman, Andrew, and Iain Pardoe.
#' "Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models."
#' _Technometrics 48, no. 2 (May 2006): 241-51_.
#'
#' Higgins, Julian P. T., and Simon G. Thompson.
#' "Quantifying Heterogeneity in a Meta-Analysis."
#' _Statistics in Medicine, vol. 21, no. 11, June 2002, pp. 1539-58_.
#'
#' Hippel, Paul T von. "The Heterogeneity Statistic I2 Can Be Biased in Small Meta-Analyses."
#' _BMC Medical Research Methodology 15 (April 14, 2015)._
#'
#'
#' @export
#'

pooling <- function(bg,
                    metric = c("pooling", "isq", "hsq", "weights"),
                    type = c("groups", "total"),
                    summary = TRUE) {
  type <- match.arg(type, c("groups", "total"))
  metric <- match.arg(metric, c("pooling", "isq", "hsq", "weights"))

  # we have to rig it for no pooling cases
  # because sigma_tau parameter might be meaningless then
  if(bg$pooling == "none")
    return(array(0, c(3, bg$n_groups, bg$n_parameters)))
  if(bg$pooling == "full")
    return(array(1, c(3, bg$n_groups, bg$n_parameters)))

  # we'll replace by switch() in the future
  # if(bg$model %in% c("rubin", "mutau", "logit", "rubin_full")) {
  if(bg$n_parameters == 1) {
    # 1-dimensional vector of S values (S=N samples)
    sigma_tau <- hypersd(bg,message=FALSE)

    # Grab the appropriate SE (k values)
    sigma_k <- switch(bg$model,
                      "mutau"       = bg$data$se.tau,
                      "rubin"       = bg$data$se,
                      # May need to apply rare event correction for
                      # the pooling models to work correctly
                      "logit"       = apply_cont_corr(bg$summary_data, 0.25, "single",
                                                      add_or = TRUE, pooling = TRUE)$se,
                      # These are SDs after pooling, don't use this (here for tests)
                      # "rubin_full"  = group_effects(bg, summary = TRUE)[, "sd", 1],
                      "rubin_full"  = bg$summary_data$se.tau,
                      "mutau_full"  = bg$summary_data$se.tau
                      )


    if(type == "groups")
      ret <- sapply(sigma_k, function(se) se^2 / (se^2 + sigma_tau^2))
    if(type == "total"){
      # nu <- mean(sigma_k^2)
      w <- 1/(sigma_k^2)
      nu <- (length(sigma_k)-1)*sum(w) / (sum(w)^2 - sum(w^2))
      ret <- replicate(1, nu / (nu + sigma_tau^2))
    }
    if(metric == "weights" && type == "groups"){
      precisions <- sapply(sigma_k, function(se) 1 / (se^2 + sigma_tau^2))
      ret <- t(apply(precisions, 1, function(x) x/sum(x)))
    }
    if(metric == "weights" && type == "total")
      stop("Weights can be calculated only for type = 'groups'")

    ret <- replicate(1, ret) #third dim is always N parameters, by convention,
                             #so we set it to 1 here

  } else if(bg$model == "quantiles") {
    # compared to individual-level, here we are dealing with 1 more dimension
    # which is number of quantiles; so if sigma_tau above is sigma of trt effect
    # here it is N effects on N quantiles etc.

    # Nq-dimensional sigma_tau
    sigma_tau <- hypersd(bg,message=FALSE)
    # for(i in 1:dim(sigma_tau)[2])
      # sigma_tau[,i,1] <- sigma_tau[,i,i]
    # sigma_tau <- sigma_tau[,,1] #rows are samples, columns are quantiles

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
    ret <- sigma_k^2 / (sigma_k^2 + sigma_tau^2)

    if(type == "total"){
      warning("Total pooling not implemented for quantiles model")
      return(array(NA, c(3, bg$n_groups, bg$n_parameters)))
    }
  } else if(bg$model == "sslab") {
    if(type == "total"){
      warning("Total pooling not implemented for spike & slab model")
      return(array(NA, c(3, bg$n_groups, bg$n_parameters)))
    }

    # Must calculate pooling for each effect vector?
    te <- hypersd(bg,message=FALSE)
    warning("In this version of baggr calculations of pooling for spike & slab may be wrong")
    warning("Please contact package authors if youy are using this feature")
    ge <- group_effects(bg, summary = T)[,"sd",]
    sigma_tau <- aperm(replicate(5, te), c(1,3,2))
    sigma_k   <- aperm(replicate(dim(sigma_tau)[1], ge), c(3,1,2))
    ret <- sigma_k^2 / (sigma_k^2 + sigma_tau^2)

  } else {
    stop("Cannot calculate pooling metrics.")
  }

  if(metric == "isq") ret <- 1-ret
  if(metric == "hsq") ret <- 1/ret
  if(metric == "h")   ret <- sqrt(1/ret)


  if(summary)
    ret <- apply(ret, c(2,3), mint)

  # if(bg$n_parameters == 1)
  #   ret <- ret[1,,]

  return(ret)

}

#' @rdname pooling
#' @export
heterogeneity <- function(
  bg,
  metric = c("pooling", "isq", "hsq", "weights"),
  summary = TRUE)
  pooling(bg, metric, type = "total", summary = summary)

#' @rdname pooling
#' @param object [baggr] model for which to calculate group (study) weights
#' @param ... Unused, please ignore.
#' @export
#' @importFrom stats weights
weights.baggr <- function(object, ...)
  pooling(object, metric = "weights", type = "groups", ...)

