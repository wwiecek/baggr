#' Leave one out cross-validation for \code{baggr} models
#'
#' Performs leave-one-out cross-validation on a \code{baggr} model at the group level.
#' This function automatically runs K `baggr` models, leaving out one group at a time,
#' and then calculating log predictive density for that group (see Gelman _et al_ 2014).
#' The main output is -2 times the log predictive density averaged over the K models,
#' which corresponds to the Watanabe-Akaike Information Criterion.
#' This function takes in the same arguments as `baggr()`, plus an option
#' (`return_models`) for whether to return all the models or just the summary statistics.
#'
#'
#' @param data Input data frame - same as for [baggr] function.
#' @param return_models logical; if FALSE, summary statistics will be returned and the
#'                      models discarded;
#'                      if TRUE, a list of models will be returned alongside summaries
#' @param ... Additional arguments passed to [baggr].
#' @return log predictive density value, an object of class `baggr_cv`;
#' full model, prior values and _lpd_ of each model are also returned.
#' These can be examined by using `attributes()` function.
#'
#' @details
#' This function can be used to understand how any one group affects
#' the overall result, as well as how well the model predicts
#' the omitted group. Because this function runs K models in total,
#' it is recommended to set `mc.cores` option before running `loocv`,
#' e.g. `options(mc.cores = 4)`.
#' Even with this option enabled, this function often has a long
#' runtime even for simple examples.
#' The main output is -2 times the log predictive density averaged
#' over `K`` models, which corresponds to the Watanabe-Akaike
#' Information Criterion. A WAIC value closer to zero (i.e. a smaller
#' number in magnitude) means a better fit.
#' For more information on cross-validation see
#' [this overview article](http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf)
#'
#' For running more computation-intensive models, consider setting the `mc.cores`
#' option before running `loocv`, e.g. `options(mc.cores = 4)` (by default `baggr` runs
#' 4 MCMC chains in parallel).
#' As a default, rstan runs "silently" (`refresh=0`). To see sampling progress, please
#' set e.g. `loocv(data, refresh = 500)`.
#'
#' @examples
#' \donttest{
#' # even simple examples may take a while
#' cv <- loocv(schools, pooling = "partial")
#' print(cv)      # returns the lpd value
#' attributes(cv) # more information is included in the object
#' }
#'
#' @author Witold Wiecek
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'

loocv <- function(data, return_models = FALSE, ...) {
  # Set prior, if not specified by the user
  args <- list(...)
  # Stop if pooling == "none"
  if(!is.null(args[["pooling"]]) && (args[["pooling"]] == "none"))
    stop("For a model with no pooling LOO CV doesn't exist.")

  # Set parallel up...
  # if(is.null(getOption("mc.cores"))){
  #   cat(paste0(
  #     "loocv() temporarily set options(mc.cores = parallel::detectCores()) \n"))
  #   temp_cores <- TRUE
  #   options(mc.cores = parallel::detectCores())
  # } else {
  #   temp_cores <- FALSE
  # }
  temp_cores <- FALSE

  # Model with all of data:
  full_fit <- try(baggr(data, refresh = 0, ...))
  if(class(full_fit) == "try-error")
    stop("Inference failed for the model with all data")

  # Prepare the arguments
  args[["data"]] <- data
  args[["model"]] <- full_fit$model
  if(!("prior" %in% names(args))) {
    message("(Prior distributions taken from the model with all data. See $prior.)")
    args[["prior"]] <- full_fit$prior
  }

  # Determine number and names of groups
  if(args[["model"]] %in% c("full", "quantiles")) {
    if(!is.null(args[["group"]]))
      group_col <- args[["group"]]
    else
      group_col <- "group"
    group_names <- unique(data[[group_col]])
    K <- length(group_names)
  } else {
    K <- nrow(data)
  }
  if(K > 50){
    message(paste("Large number of groups - note that the
                   model will be repeated", K, "times"))
  } else {
    message(paste0("Repeating baggr() for ", K, " separate models"))
  }

  # LOO CV models
  cat("\n")
  pb <- utils::txtProgressBar(style = 3)
  kfits <- lapply(as.list(1:K), function(i) {
    if(args[["model"]] %in% c("full", "quantiles")) {
      args$data      <- data[data[[group_col]] != group_names[i], ]
      args$test_data <- data[data[[group_col]] == group_names[i], ]
    } else {
      args$data      <- data[-i,]
      args$test_data <- data[i,]
    }

    utils::setTxtProgressBar(pb, (i-1)/K)

    # Run baggr models:
    if(is.null(args[["refresh"]]))
      args[["refresh"]] <- 0
    res <- do.call(baggr, args)

    # Sanitized version:
    # res <- try(do.call(baggr, args))
    # if(class(res) == "try-error") {
    #   stop(paste0("Inference failed for model number ", i, " out of ", K))
    # } else {
    #   return(res)
    # }
  })
  utils::setTxtProgressBar(pb, 1)
  close(pb)

  tau_estimate <-
    lapply(kfits, function(x) mean(treatment_effect(x)[[1]]))
  sd_estimate <-
    lapply(kfits, function(x) mean(treatment_effect(x)[[2]]))
  loglik <-
    lapply(kfits, function(x) apply(as.matrix(x$fit, "logpd"), 2, mean))

  out <- list(
    lpd = -2*sum(unlist(loglik)),
    df  = data.frame(
      "tau" = unlist(tau_estimate),
      "sigma_tau" = unlist(sd_estimate),
      "lpd" = unlist(loglik)),
    full_model = full_fit,
    prior = args[["prior"]]
  )
  if(return_models)
    out$models <- kfits
  class(out) <- "baggr_cv"

  if(temp_cores)
    options(mc.cores = NULL)

  return(out)
}

#' @export
print.baggr_cv <- function(x, ...) {
  print(paste("log predictive density = ", format(x$lpd)))
}
