# this pooling metric function should be made more general so that it's
# 1 function to prepare data for pooling calculations
# ...followed by 1 function that calculates (different methods of calc.)
# the later function will be more generalisable to results not necc. from baggr()
pooling <- function(baggr, interval = FALSE) {

  # replace by switch() in the future
  if(baggr$model == "rubin") {
    sigma_tau <- rstan::extract(baggr$fit, pars = "sigma_tau")[[1]]
    sigma_k <- baggr$data$se
  } else if(baggr$model == "mutau") {
    sigma_tau <- rstan::extract(baggr$fit, pars = "sigma_tau")[[1]]
    sigma_k <- baggr$data$se.tau
  } else if(baggr$model == "joint") {
    #now it will be vectors, not a single value
    m <- as.matrix(baggr$fit)
    sigma_k <- m[, grepl("^sigma_y_k", colnames(m))]
    sigma_tau <- rstan::extract(baggr$fit, pars = "sigma_mutau[2,2]")[[1]]
    ret <- t(apply(sigma_k, 2, function(sigma) {
      x <- sigma^2 / (sigma^2 + sigma_tau^2)
      quantile(x, c(.025, .5, .975))
    }))
    return(ret)
  }
  ret <- t(sapply(sigma_k, function(se) {
    x <- se^2 / (se^2 + sigma_tau^2)
    quantile(x, c(.025, .5, .975))
  }))
  return(ret)
}
