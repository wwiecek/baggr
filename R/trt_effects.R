treatment_effect <- function(bg) {
  if(bg$model == "rubin") {
    tau <- as.matrix(bg$fit)[,"tau"]
    sigma_tau <- as.matrix(bg$fit)[,"sigma_tau"]
  } else if(bg$model == "joint") {
    tau <- as.matrix(bg$fit)[,"mutau[2]"]
    sigma_tau <- as.matrix(bg$fit)[,"sigma_mutau[2,2]"]
  }

  return(list(tau = tau, sigma_tau = sigma_tau))
}
