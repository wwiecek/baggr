# given a baggr model
# extract list of 2 vectors: tau and sigma_tau

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
      sigma_tau <- sigma_tau[,2,2]
    }
  } else if(bg$model == "joint") {
    tau <- as.matrix(bg$fit)[,"mutau[2]"]
    sigma_tau <- as.matrix(bg$fit)[,"sigma_mutau[2,2]"]
  }
  return(list(tau = tau, sigma_tau = sigma_tau))
}
