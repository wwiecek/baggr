# given a baggr model
# extract list of 2 vectors: tau and sigma_tau
# these are not summarised but full sample

treatment_effect <- function(bg) {
  if(bg$pooling == "none"){
    message("There is no treatment effect estimated when pooling = 'none'.")
    return(NULL)
  }
  if(bg$model %in% c("rubin", "mutau")) {
    tau <- rstan::extract(bg$fit, pars="tau")[[1]]
    sigma_tau <- rstan::extract(bg$fit, pars="sigma_tau")[[1]]
    if(bg$model == "mutau") {
      tau <- tau[,2]
      # in model with correlation, we have Var(), not SD()
      sigma_tau <- sqrt(sigma_tau[,2,2])
    }
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
