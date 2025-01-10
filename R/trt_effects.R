
#' Average treatment effects in a baggr model
#'
#' The most general `treatment_effect` displays
#' both hypermean and hyperSD (as a list of length 2),
#' whereas `hypermean` and `hypersd` can be used as shorthands.
#'
#' @param bg a [baggr] model
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param interval uncertainty interval width (numeric between 0 and 1), if summarising
#' @param transform a transformation to apply to the result, should be an R function;
#'                  (this is commonly used when calling `treatment_effect` from other
#'                  plotting or printing functions)
#' @param message logical; use to disable messages prompted by using with
#'                no pooling models
#' @describeIn treatment_effect A list with 2 vectors (corresponding to MCMC samples)
#'         `tau` (mean effect) and `sigma_tau` (SD). If `summary=TRUE`,
#'         both vectors are summarised as mean and lower/upper bounds according to
#'         `interval`
#' @export
#' @importFrom rstan extract


treatment_effect <- function(bg, summary = FALSE,
                             transform = NULL, interval = .95,
                             message = TRUE) {
  check_if_baggr(bg)

  if(bg$pooling == "none"){
    if(message)
      message("There is no treatment effect estimated when pooling = 'none'.")
    return(list(tau = as.numeric(NA), sigma_tau = as.numeric(NA)))
  }
  if(bg$model %in% c("rubin", "mutau", "mutau_full", "logit", "rubin_full")) {
    tau <- rstan::extract(bg$fit, pars="mu")[[1]]
    if(bg$model %in% c("mutau", "mutau_full"))
      tau <- tau[,1,2] #in this mutau case the second element is trt effect, first is baseline
    if(is.matrix(tau) && ncol(tau) == 1) tau <- c(tau) #convert to vector for 1D cases

    if(bg$pooling == "partial"){
      if(bg$model %in% c("mutau", "mutau_full"))
        sigma_tau <- rstan::extract(bg$fit, pars="hypersd[1,2]")[[1]]
      else
        sigma_tau <- rstan::extract(bg$fit, pars="tau")[[1]]
      if(is.matrix(tau) && ncol(sigma_tau) == 1) sigma_tau <- c(sigma_tau)
    }
    if(bg$pooling == "full")
      sigma_tau <- rep(0, length(bg$effects)) #same dim as tau, but by convention set to 0
  } else if(bg$model == "quantiles") {
    # In this case we have N columns = N quantiles
    tau <- as.matrix(bg$fit, "beta_1")
    # only take diagonals:
    sigma_tau <- t(apply(rstan::extract(bg$fit, "Sigma_1")[[1]], 1, diag))
    # in model with correlation, we have Var(), not SD()
    sigma_tau <- sqrt(sigma_tau)
  } else if(bg$model == "sslab") {
    mean_params <- c("tau[1]", "tau[2]",
                     "sigma_TE[1]", "sigma_TE[2]",
                     "kappa[1,1,2]", "kappa[1,2,2]")
    sigma_params <- c("hypersd_tau[1]", "hypersd_tau[2]",
                      "hypersd_sigma_TE[1]", "hypersd_sigma_TE[2]",
                      "hypersd_kappa[1,1,2]", "hypersd_kappa[1,2,2]")
    tau <- as.matrix(bg$fit, mean_params)
    sigma_tau <- as.matrix(bg$fit, sigma_params)

  } else {
    stop("Can't calculate treatment effect for this model.")
  }
  if(length(bg$effects) > 1){
    colnames(tau) <- bg$effects
    if(bg$pooling != "full") colnames(sigma_tau) <- bg$effects
  }
  if(!is.null(transform)){
    tau <- do.call(transform, list(tau))
    sigma_tau <- NA # by convention we set it to NA so that people don't transform
    # and then do operations on it by accident
  }
  if(summary) {
    tau <- mint(tau, int=interval, median=TRUE, sd = TRUE)
    sigma_tau <- mint(sigma_tau, int=interval, median=TRUE, sd = TRUE)
  }

  return(list(tau = tau, sigma_tau = sigma_tau))
}



#' @describeIn treatment_effect The hypermean of a `baggr` model, shorthand for `treatment_effect(x, s=T)[[1]]`
#' @export
hypermean <- function(bg,transform=NULL,interval = 0.95,message=FALSE, summary=TRUE) {
  t <- treatment_effect(bg,summary=summary,transform=transform,interval=interval,message=FALSE)
  if(message){
    cat(paste0("Hypermean of ",bg$effects," with ",interval*100,"% interval:\n"))
    print(t[[1]])
  } else
    t[[1]]
}

#' @describeIn treatment_effect The hyper-SD of a `baggr` model, shorthand for `treatment_effect(x, s=T)[[2]]`
#' @export
hypersd <- function(bg,transform=NULL,interval = 0.95,message=FALSE, summary=TRUE) {
  t <- treatment_effect(bg,summary=summary,transform=transform,interval=interval,message=FALSE)
  if(message){
    cat(paste0("Hyper SD of ",bg$effects," with ",interval*100,"% interval:\n"))
    print(t[[2]])
  } else
    t[[2]]
}


#' Correlation between mu and tau in a baggr model
#'
#' @param bg  a [baggr] model where `model = "mutau"`
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param interval uncertainty interval width (numeric between 0 and 1), if summarising
#' @return a vector of values
#' @export
mutau_cor <- function(bg,
                      summary = FALSE,
                      interval = 0.95) {
  # m <- matrix(apply(as.matrix(bg$fit, "L_Omega"), 2, mean), 2, 2)
  # h <- apply(as.matrix(bg$fit, "hypersd"), 2, mean)
  # (m %*% t(m))[2,1]
  # diag(h) %*% (m %*% t(m)) %*% diag(h)

  m <- as.matrix(bg$fit, "L_Omega")
  # we want entry (2,1) in LL^T which is simply (1,1)*(2,1)
  # in a lower-triangular matrix; and (1,1) should be == 1 everywhere
  if(!all(m[,1] == 1))
    warning("Error with correlation matrix in the mu&tau model. Inspect L_Omega parameters.")
  if(summary)
    return(mint(m[,2], int=interval, median=TRUE, sd = TRUE))
  else
    return(m[,2])

}

