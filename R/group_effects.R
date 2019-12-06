#' Extract baggr study effects
#'
#' Given a baggr object, returns the raw MCMC draws of the posterior for
#' each group's effect, or a summary of these draws. This is an internal
#' function currently used as a helper for plotting and printing of results.
#'
#' @param bg baggr object
#' @param summary logical; if TRUE returns summary statistics as explained below.
#' @param interval uncertainty interval width (numeric between 0 and 1), if summarising
#'
#' @return Either a matrix with MCMC samples (if summary = FALSE)
#'         or a summary of these samples (if summary = TRUE).
#' @examples
#' fit1 <- baggr(schools)
#' group_effects(fit1, summary = TRUE, interval = 0.5)
#' @details If summary = TRUE, the returned object contains for each study
#' or group: the posterior medians, the lower and upper bounds of the
#' uncertainty intervals using the central posterior credible interval
#' of width specified in the argument "interval", the posterior mean, and
#' the posterior standard deviation.
#'
#' @export

group_effects <- function(bg, summary = FALSE, transform = NULL, interval = .95) {
  check_if_baggr(bg)

  # m <- as.matrix(bg$fit)
  if(attr(bg , "ppd"))
    stop("There are no group effects in prior predictive distribution baggr objects.")

  if(bg$pooling == "full"){
    tau <- treatment_effect(bg)[["tau"]]
    k <- attr(bg$inputs, "n_groups")
    m <- replicate(k, tau)
    if(length(dim(m)) == 3)
      m <- aperm(m, c(1, 3, 2))

  } else {
    # choose correct columns for the given models:
    if(bg$model %in% c("rubin", "mutau", "logit")) {
      #replace by extract:
      # m <- m[, grepl("^tau_k", colnames(m))]
      m <- rstan::extract(bg$fit, pars = "tau_k")[[1]]
      # drop mu if model has mu (baseline/control value)
      if(bg$model == "mutau")
        m <- m[,,2]
    } else if(bg$model == "full") {
      m <- rstan::extract(bg$fit, pars = "mutau_k")[[1]][,,2]
    } else if(bg$model == "quantiles") {
      # In this case we have 3D array, last dim is quantiles
      m <- rstan::extract(bg$fit, pars = "beta_1_k")[[1]]
    }
  }
  # for consistency with quantiles, except we have 1 parameter only
  if(length(dim(m)) == 2)
    m <- array(m, dim = c(dim(m), 1))

  par_names <- attr(bg$inputs, "group_label")

  if(!is.null(par_names))
    dimnames(m)[[2]] <- par_names
  else
    dimnames(m)[[2]] <- paste0("Group ", 1:attr(bg$inputs, "n_groups"))

  # will summarise if requested:
  if(summary) {
    intval <- c((1-interval)/2, .5, 1 - (1-interval)/2)
    m <- apply(m, c(2,3), function(x) c(quantile(x, intval), mean(x), sd(x)))
    if(is.null(dimnames(m)[[2]]))
      dimnames(m)[[2]] <- 1:nrow(m)
    dimnames(m)[[1]] <- c("lci", "median", "uci", "mean", "sd")
    m <- aperm(m, c(2,1,3))
  }

  if(!is.null(transform))
    m <- do.call(transform, list(m))

  return(m)
}
