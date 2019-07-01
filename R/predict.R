#' Predict method for baggr objects
#' @param x model to predict from
#' @param ... other arguments to pass to predict
#' @export
#' @details This is currently only implemented for the baggr model
#' but could also be extended to the others. Currently a WIP, please
#' report bugs.
predict.baggr <- function(x, newdata = NULL,
                          allow_new_levels = T, nsamples, ...) {
  switch(x$model,
         rubin = predict_rubin(x, newdata = newdata,
                               allow_new_levels = allow_new_levels,
                               nsamples = nsamples, ...),
         quantiles = predict_quantiles(x, newdata = newdata,
                                       allow_new_levels = allow_new_levels,
                                       nsamples = nsamples,  ...))
}

#' Make model matrix for the rubin data
#' @importFrom reshape2 dcast
#' @param x model to get data from
#' @param newdata new data to use with model
#' @param allow_new_levels whether to allow for unobserved groups
rubin_data <- function(x, newdata = NULL, allow_new_levels = T) {
  if(x$model != "rubin") {
    stop("Model must be type Rubin.")
  }

  group_label <- attr(x$inputs, "group_label")
  group_num <- 1:attr(x$inputs, "n_groups")
  group_col <- names(which(sapply(x$data,
                            function(x, labels) any(labels %in% x),
                            labels = group_label)))

  dat <- x$data

  other_cols <- setdiff(colnames(dat), group_col)

  if(is.null(newdata)) {
    dat <- x$data
  } else {
    dat <- newdata
  }

  dat[,group_col] <-
        factor(dat[,group_col],
               levels = as.character(group_label))


  if(allow_new_levels != T) {
    if(any(is.na(dat[,group_col]))) stop("Data contains new levels. If this behavior is desired, set allow_new_levels to TRUE.")
  }

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
#' @param newdata new data to predict, defaults to NULL
#' @param allow_new_levels allow the predictive of new, unobserved groups
#' @param nsamples number of samples to predict
predict_rubin <- function(x,
                          newdata = NULL,
                          allow_new_levels = T,
                          nsamples,
                          ...) {
  pred_data <- rubin_data(x, newdata)

  sigmas <- sapply(x$data$se, rep, times = nsamples)

  params <- rstan::extract(x$fit, c("tau","sigma_tau","tau_k","eta"))

  tau_k <- params$tau_k[1:nsamples,]
  eta <- params$eta[1:nsamples,]

  tau <- params$tau[1:nsamples]
  sigma_tau <- params$sigma_tau[1:nsamples]

  means <- cbind(tau, eta)

  pred_means <- pred_data %*% t(means)

  epsilon <- rnorm(length(sigmas),
                   mean = 0, sigmas)

  pp_dist <- pred_means + epsilon

  t(pp_dist)

}

#' Stop with informative error
stop_not_implemented <- function() {
  stop("Method not implemented.")
}

#' Get Y from various models
#' @param x model to get y
#' @details grabs the relevant Y variable
#' for use with posterior or prior predictive checks.
get_y <- function(x) {
  switch(x$model,
         rubin = "tau",
         mutau = "tau",
         quantiles = stop_not_implemented(),
         full = stop_not_implemented()
         )
}

#' Posterior predictive checks for baggr model
#' @import bayesplot
#' @param x Model to check
#' @param type type of pp_check. For a list see \link{bayesplot::available_ppc}
#' @param nsamples number of samples to compare
#' @export
pp_check.baggr <- function(x, type = "dens_overlay", nsamples = 40) {
  pp_fun <- getFromNamespace(paste0("ppc_",type),ns = "bayesplot")
  model_type <- x$model
  y <- x$data[,get_y(x)]
  yrep <- predict(x, nsamples = nsamples)
  pp_fun(y, yrep)
}
