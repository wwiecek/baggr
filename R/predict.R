#' Predict method for baggr objects
#' @param object model to predict from
#' @param newdata optional, new data to predict observations from
#' @param allow_new_levels whether to allow the model to make predictions
#' about unobserved groups. Without additional group-level information
#' the model will use the unconditional, pooled estimate.
#' @param nsamples Number of samples to draw from the posterior.
#' Cannot exceed the number of samples in the fitted model.
#' @param ... other arguments to pass to predict function
#'            (currently not used)
#' @export
predict.baggr <- function(object, nsamples,
                          newdata = NULL,
                          allow_new_levels = TRUE,
                          ...) {
  switch(object$model,
         rubin = predict_rubin(object,
                               nsamples = nsamples,
                               newdata = newdata,
                               allow_new_levels = allow_new_levels),
         predict_unknown(object))
}



#' Predict method for model that is unknown or not implemented
#' @param x baggr model to generate predictions from
predict_unknown <- function(x) {
  stop("The ", x$model, " model",
       " does not yet have a posterior prediction method",
       " implemented. If you would like this to be implemented,",
       " open a github issue here: \n",
       "https://github.com/wwiecek/baggr/issues")
}

#' Make model matrix for the rubin data
#' @param x model to get data from
#' @param newdata new data to use with model
#' @param allow_new_levels whether to allow for unobserved groups
rubin_data <- function(x, newdata = NULL, allow_new_levels = TRUE) {
  check_if_baggr(x)
  # if(x$model != "rubin")
    # stop("Model must be type Rubin.")

  group_label <- attr(x$inputs, "group_label")
  group_num   <- 1:attr(x$inputs, "n_groups")
  group_col   <- names(which(sapply(x$data,
                                    function(x, labels) any(labels %in% x),
                                    labels = group_label)))
  dat <- x$data
  other_cols <- setdiff(colnames(dat), group_col)
  if(is.null(newdata)) {
    dat <- x$data
  } else {
    dat <- newdata
  }
  dat[,group_col] <- factor(dat[,group_col],
                            levels = as.character(group_label))

  if(allow_new_levels != T)
    if(any(is.na(dat[,group_col])))
      stop("Data contains new levels. If this behavior is desired,",
           "set allow_new_levels to TRUE.")

  predmat <- matrix(nrow = nrow(dat), ncol = x$n_groups)
  for(i in 1:ncol(predmat)) {
    lvl <- as.integer(dat[i,group_col])
    predmat[lvl,i] <- 1
  }
  predmat[which(is.na(predmat))] <- 0
  cbind(1, predmat)
}

#' Predict function for the rubin model
#' @importFrom rstan extract
#' @param x model to predict from
#' @param nsamples number of samples to predict
#' @param newdata new data to predict, defaults to NULL
#' @param allow_new_levels allow the predictive of new, unobserved groups
predict_rubin <- function(x,
                          nsamples,
                          newdata = NULL,
                          allow_new_levels = TRUE) {
  if(missing(nsamples)){
    nsamples <- get_n_samples(x)
  }

  pred_data <- rubin_data(x, newdata)
  se <- sapply(x$data$se, rep, times = nsamples)

  eta <- rstan::extract(x$fit, c("eta"))[[1]][1:nsamples,]
  tau <- treatment_effect(x)$tau[1:nsamples]

  pred_means <- pred_data %*% t(cbind(tau, eta))
  epsilon <- rnorm(length(se),
                   mean = 0, se)
  pp_dist <- pred_means + epsilon
  t(pp_dist)
}

#' Predict function for the mu & tau model
#' @importFrom rstan extract
#' @param x model to predict from
#' @param nsamples number of samples to predict
#' @param newdata new data to predict, defaults to NULL
#' @param allow_new_levels allow the predictive of new, unobserved groups
predict_mutau <- function(x,
                          nsamples,
                          newdata = NULL,
                          allow_new_levels = TRUE) {
  if(missing(nsamples))
    nsamples <- get_n_samples(x)

  stop_not_implemented()
}

#' Predict function for the quantiles model
#' @importFrom rstan extract
#' @param x model to predict from
#' @param nsamples number of samples to predict
#' @param newdata new data to predict, defaults to NULL
#' @param allow_new_levels allow the predictive of new, unobserved groups
predict_quantiles <- function(x,
                              nsamples,
                              newdata = NULL,
                              allow_new_levels = TRUE){
  if(missing(nsamples))
    nsamples <- get_n_samples(x)

  stop_not_implemented()
}

#' Posterior predictive checks for baggr model
#'
#' Performs posterior predictive checks with the
#' \pkg{bayesplot} package.
#'
#' @param x Model to check
#' @param type type of pp_check. For a list see
#'             \pkg{\link[bayesplot:available_ppc]{here}}.
#' @param nsamples number of samples to compare
#' @aliases pp_check
#'
#' @details For a detailed explanation of each of the ppc functions,
#' see the \code{\link[bayesplot:PPC-overview]{PPC}}
#' documentation of the `bayesplot`
#' package.
#'
#' @import bayesplot
#' @importFrom utils getFromNamespace
#' @export

pp_check.baggr <- function(x, type = "dens_overlay", nsamples = 40) {
  pp_fun <- utils::getFromNamespace(paste0("ppc_",type),ns = "bayesplot")
  col <- switch(x$model,
                rubin = "tau",
                mutau = "tau",
                quantiles = stop_not_implemented(),
                full = stop_not_implemented()
  )
  y <- x$data[,col]
  yrep <- predict(x, nsamples = nsamples)
  pp_fun(y, yrep)
}

# Helper functions -----
#' Stop with informative error
stop_not_implemented <- function() {
  stop("Method not implemented.")
}

#' Extract number of samples from a baggr object
#' @param x baggr fit to get samples from
#' @details Checks for number of iterations and
#' number of Markov chains, returns maximum number
#' of valid samples
get_n_samples <- function(x) {
  check_if_baggr(x)
  iter <- attr(x$fit, "stan_args")[[1]]$iter - attr(x$fit, "stan_args")[[1]]$warmup
  chains <- max(sapply(x$fit@stan_args,function(x) x$chain_id))
  iter * chains
}
