#' Pooling metrics for baggr
#'
#' Compute pooling metrics given a baggr meta-analysis model
#'
#' @param bg output of a baggr() function
#' @param metric for now we use \code{"gelman-hill"}, other options will be added in the future
#' @param summary logical; if FALSE a whole vector of pooling statistics is returned,
#'                otherwise only the means and intervals
#'
#' @details
#' Pooling statistic describes the extent to which group-level estimates of treatment
#' effect are shrunk toward average treatment effect in the meta-analysis model.
#'
#' @return Matrix with mean and intervals for chosen pooling metric,
#'         each row corresponding to one meta-analysis group.
#' @author Witold Wiecek, Rachael Meager


pooling <- function(bg, metric = "gelman-hill", summary = TRUE) {
  # we have to rig it for no pooling cases
  # because sigma_tau parameter might be meaningless then
  if(bg$pooling == "none")
    return(array(0, c(3, bg$n_groups, bg$n_parameters)))
  if(bg$pooling == "full")
    return(array(1, c(3, bg$n_groups, bg$n_parameters)))

  # we'll replace by switch() in the future
  if(bg$model == "rubin" || bg$model == "mutau") {

    if(bg$model == "mutau"){
      sigma_k <- bg$data$se.tau
      sigma_tau <- treatment_effect(bg)$sigma_tau
    }
    if(bg$model == "rubin"){
      sigma_k <- bg$data$se
      sigma_tau <- treatment_effect(bg)$sigma_tau
    }
    # for summary level data cases:
    ret <- sapply(sigma_k, function(se) se^2 / (se^2 + sigma_tau^2))
    ret <- replicate(1, ret) #third dim is always N parameters, by convention

  } else if(bg$model == "full") {
    # note that we use a point estimate for sigma_k, so it's not fully Bayesian
    # but in the end we usually report a point estimate so this is acceptable for now
    sigma_k <- group_effects(bg, summary = TRUE)[, "sd", 1]
    sigma_tau <- treatment_effect(bg)$sigma_tau
    # we add 3rd dimension to allow more than 1 parameter (in the future), e.g. Var
    ret <- array(0, dim = c(length(sigma_tau), bg$n_groups, 1))
    ret[,,1] <- sapply(sigma_k, function(sigma) sigma^2 / (sigma^2 + sigma_tau^2))

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


  }

  if(summary)
    ret <- apply(ret, c(2,3), mint)
  return(ret)

}
