#' Given a model return treatment effects
#' (Used as a helper for plotting and printing of results)

study_effects <- function(bg, interval = FALSE) {
  if(class(bg) != "baggr")
    stop("study_effects only works with 'baggr' class objects")

  m <- as.matrix(bg$fit)

  if(bg$pooling == "full"){
    # tau <- m[,"tau"]
    tau <- treatment_effect(bg)[["tau"]] #for consistency we have a separate function for this
    k <- attr(bg$inputs, "n_sites")
    m <- matrix(tau, nrow(m), k, byrow = F)
    # return(NULL) #for now we don't compute them at all, they're all the same as avg effect
  } else{
    # choose correct columns for the given models:
    if(bg$model %in% c("rubin", "mutau")) {
      #replace by extract:
      # m <- m[, grepl("^tau_k", colnames(m))]
      m <- rstan::extract(bg$fit, pars = "tau_k")[[1]]
      # drop mu if model has mu (baseline/control value)
      if(bg$model == "mutau")
        m <- m[,,2]
    } else if(bg$model == "joint") {
      m <- rstan::extract(bg$fit, pars = "mutau_k")[[1]][,,2]
      #replace by extract:
      # m <- m[, grepl("^mutau_k", colnames(m))]
      # m <- m[, grepl(",2]", colnames(m))]
    }
  }
  par_names <- attr(bg$inputs, "site_label")
  if(!is.null(par_names))
    colnames(m) <- par_names
  else
    colnames(m) <- paste0("Site ", 1:attr(bg$inputs, "n_sites"))

  # will summarise if requested:
  if(interval) {
    intval <- c((1-interval)/2, .5, 1 - (1-interval)/2)
    m <- t(apply(m, 2, function(x) quantile(x, intval)))
    if(is.null(rownames(m)))
      rownames(m) <- 1:nrow(m)
    m <- data.frame(rownames(m), m)
    colnames(m) <- c("parameter", "lci", "median", "uci")
  }

  return(m)
}
