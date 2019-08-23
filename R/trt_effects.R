
#' Average treatment effect in a baggr model
#'
#' @param bg a [baggr] model
#' @return A list with 2 vectors (corresponding to MCMC samples)
#'         `tau` (mean effect) and `sigma_tau` (SD)
#' @export
#' @importFrom rstan extract

treatment_effect <- function(bg) {
  if(class(bg) != "baggr")
    stop("treatment_effect requires a baggr object")
  if(bg$pooling == "none"){
    message("There is no treatment effect estimated when pooling = 'none'.")
    return(NULL)
  }
  if(bg$model %in% c("rubin", "mutau")) {
    tau <- rstan::extract(bg$fit, pars="tau")[[1]]
    if(bg$model == "rubin")
      tau <- c(tau)
    if(bg$model == "mutau")
      tau <- tau[,1,2]

    if(bg$pooling == "partial"){
      sigma_tau <- rstan::extract(bg$fit, pars="sigma_tau")[[1]]
      if(bg$model == "rubin")
        sigma_tau <- c(sigma_tau)
      if(bg$model == "mutau")
        sigma_tau <- sqrt(sigma_tau[,1,2,2])
    }
    if(bg$pooling == "full")
      sigma_tau <- 0 #same dim as tau, but by convention set to 0
  } else if(bg$model == "full") {
    tau <- as.matrix(bg$fit)[,"mutau[2]"]
    sigma_tau <- as.matrix(bg$fit)[,"sigma_mutau[2,2]"]
    # in model with correlation, we have Var(), not SD():
    sigma_tau <- sqrt(sigma_tau)
  } else if(bg$model == "quantiles") {
    # In this case we have N columns = N quantiles
    tau <- as.matrix(bg$fit, "beta_1")
    # only take diagonals:
    sigma_tau <- t(apply(rstan::extract(bg$fit, "Sigma_1")[[1]], 1, diag))
    # in model with correlation, we have Var(), not SD()
    sigma_tau <- sqrt(sigma_tau)
  }

  return(list(tau = tau, sigma_tau = sigma_tau))
}
