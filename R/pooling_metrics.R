# this pooling metric function should be made more general so that it's

#' Pooling metrics
#'
#' Compute pooling metrics (of a few different varieties)
#' given a baggr meta-analysis model
#'
#' @param bg output of a baggr() function
#' @param metric \code{"gelman-hill"} or \code{WIP}
#'
#' @details
#' Pooling statistic describes the extent to which group-level estimates of treatment effect are shrunk toward
#' average treatment effect in the meta-analysis model.
#'
#' Different measures of pooling can be estimated and choosing the right one depends on the research questions.
#' The default (\code{"gelman-hill"}) "pooling factor" statistic by Gelman & Hill (2007) [REFERENCES ARE WIP] is
#' \deqn{\omega(\tau_k) = \frac{WIP}}
#'
#'
#' @return Matrix with mean and intervals for chosen pooling metric, each row corresponding to one meta-analysis group.
#' @author Witold Wiecek, Rachael Meager


pooling <- function(bg, metric = "gelman-hill") {
  if(bg$model == "quantiles")
    return(0)

  # we have to rig it for no pooling cases
  # because sigma_tau parameter might be meaningless then
  if(bg$pooling == "none")
    return(matrix(0, bg$n_sites, 3))
  if(bg$pooling == "full")
    return(matrix(1, bg$n_sites, 3))

  # we'll replace by switch() in the future
  if(bg$model == "rubin") {
    sigma_tau <- rstan::extract(bg$fit, pars = "sigma_tau")[[1]]
    sigma_k <- bg$data$se
  } else if(bg$model == "mutau") {
    sigma_tau <- rstan::extract(bg$fit, pars = "sigma_tau")[[1]]
    sigma_k <- bg$data$se.tau
  } else if(bg$model == "full") {
    #now it will be vectors, not a single value
    m <- as.matrix(bg$fit)
    #replace by extract:
    sigma_k <- m[, grepl("^sigma_y_k", colnames(m))]
    sigma_tau <- rstan::extract(bg$fit, pars = "sigma_mutau[2,2]")[[1]]
    ret <- t(apply(sigma_k, 2, function(sigma) {
      x <- sigma^2 / (sigma^2 + sigma_tau^2)
      quantile(x, c(.025, .5, .975))
    }))
    return(ret)
  } else if(bg$model == "quantiles") {
    # compared to individual-level, here we are dealing with 1 more dimension
    # which is number of quantiles; so if sigma_tau above is sigma of trt effect
    # here it is N effects on N quantiles etc.

    # discard everything off-diagonal, we only care about variances here:
    sigma_tau <- rstan::extract(bg$fit, "Sigma_1")[[1]]
    for(i in 1:dim(sigma_tau)[2])
      sigma_tau[,i,1] <- sigma_tau[,i,i]
    sigma_tau <- sigma_tau[,,1] #rows are samples, columns are

  }

  # for summary level data cases:
  ret <- t(sapply(sigma_k, function(se) {
    x <- se^2 / (se^2 + sigma_tau^2)
    quantile(x, c(.025, .5, .975))
  }))
  return(ret)
}
