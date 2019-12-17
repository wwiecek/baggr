#' Leave one group out cross-validation for \code{baggr} models
#'
#' Performs exact leave-one-group-out cross-validation on a baggr model.
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
#'
#' The values returned by `loocv()` can be used to understand how any
#' one group affects the overall result, as well as how well the model
#' predicts the omitted group.
#'
#' This function automatically runs K baggr models, leaving out one group at a time,
#' and then calculates expected log predictive density (ELPD) for
#' that group (see Gelman et al 2013). The main output is the cross-validation
#' information criterion, or -2 times the ELPD averaged over 'K' models.
#' This is related to, and often approximated by, the Watanabe-Akaike
#' Information Criterion. A value closer to zero (i.e. a smaller number in magnitude)
#' means a better fit. For more information on cross-validation see
#' [this overview article](http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf)
#'
#' For running more computation-intensive models, consider setting the mc.cores option before running loocv, e.g. options(mc.cores = 4)
#' (by default baggr runs 4 MCMC chains in parallel).
#' As a default, rstan runs "silently" (refresh=0).
#' To see sampling progress, please set e.g. loocv(data, refresh = 500).
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
#' @references Gelman, Andrew, et al. Bayesian data analysis. Chapman and Hall/CRC, 2013.
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
  if(full_fit$model == "full")
    stop("LOO CV is not implemented for full data models yet.")
  if(!("prior" %in% names(args))) {
    message("(Prior distributions taken from the model with all data. See $prior.)")
    args[["formatted_prior"]] <- full_fit$formatted_prior
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
    lapply(kfits, function(x) apply(as.matrix(x$fit, "logpd[1]"), 2, mean))

  elpds <- unlist(loglik)

  out <- list(
    se = sqrt(length(elpds) * var(elpds)),
    elpd = sum(unlist(elpds)),
    looic = -2*sum(unlist(elpds)),
    df  = data.frame(
      "tau" = unlist(tau_estimate),
      "sigma_tau" = unlist(sd_estimate),
      "lpd" = unlist(loglik)),
    full_model = full_fit,
    prior = args[["prior"]],
    K = K,
    pointwise = elpds
  )
  if(return_models)
    out$models <- kfits
  class(out) <- "baggr_cv"

  if(temp_cores)
    options(mc.cores = NULL)

  return(out)
}

#' Check if something is a baggr_cv object
#' @param x object to check
is.baggr_cv <- function(x) {
  inherits(x, "baggr_cv")
}

#' Compare fitted models on loo
#' @param x An object of class "baggr_cv" or a list of such objects.
#' @param ... Additional objects of class "baggr_cv"
#' @importFrom loo loo_compare
#' @export loo_compare
#' @examples
#' # 2 models with more/less informative priors
#' cv_1 <- loocv(schools, model = "rubin", pooling = "partial")
#' cv_2 <- loocv(schools, model = "rubin", pooling = "partial", prior_hypermean = normal(0, 5), prior_hypersd = cauchy(0,4))
#' loo_compare(cv_1, cv_2)
#' @export
loo_compare <- function(x, ...) {
  UseMethod("loo_compare")
}

#' @aliases loo_compare
#' @export
loo_compare.baggr_cv <- function(x, ...) {
  if (is.baggr_cv(x)) {
    dots <- list(...)
    loos <- c(list(x), dots)
  } else {
    if (!is.list(x) || !length(x)) {
      stop("'x' must be a list if not a 'loo' object.")
    }
    loos <- x
  }
  if (!all(sapply(loos, is.baggr_cv))) {
    stop("All inputs should have class 'baggr_cv'.")
  }
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.")
  }
  Ns <- sapply(loos, function(x) nrow(x$df))
  if (!all(Ns == Ns[1L])) {
    stop("Not all models have the same number of data points.")
  }

  elpds <- lapply(loos, function(x) x$pointwise)
  comp <- Reduce(cbind, elpds)

  diffs <- list()

  for(i in 2:ncol(comp)) {
    diffname <- paste0("Model_1",
                       " - ",
                       "Model", i)

    rawdiffs <- comp[,1] - comp[,i]

    diffs[[i-1]] <-
      matrix(nrow = 1, ncol = 2,
             c(sum(rawdiffs),
                     sqrt(length(rawdiffs))*sd(rawdiffs)),
               dimnames = list(diffname, c("ELPD", "ELPD SE")))
  }

  diffs <- Reduce(rbind, diffs)


  class(diffs) <- c("compare_baggr_cv", class(comp))
  diffs
}

#' Print baggr_cv comparisons
#' @param x baggr_cv comparison to print
#' @param digits number of digits to print
#' @importFrom testthat capture_output
#' @importFrom crayon bold
#' @export
print.compare_baggr_cv <- function(x, digits = 3, ...) {

  mat <- as.matrix(x)

  class(mat) <- "matrix"

  cat(
    crayon::bold(paste0("Comparison of cross-validation\n")),
    "\n",
    testthat::capture_output(print(signif(mat, digits = digits)))
  )

}

#' Print baggr cv objects nicely
#' @param x baggr_cv object to print
#' @importFrom testthat capture_output
#' @importFrom crayon bold
#' @export
print.baggr_cv <- function(x, digits = 3, ...) {

  mat <- matrix(nrow = 2, ncol = 2)

  mat[1,] <- c(x$elpd, x$se)
  mat[2,] <- c(x$looic, -2*x$se)

  colnames(mat) <- c("Estimate", "Standard Error")
  rownames(mat) <- c("elpd", "looic")

  cat(
    crayon::bold(paste0("Based on ", x$K, "-fold cross-validation\n")),
    "\n",
    testthat::capture_output(print(signif(mat, digits = digits)))
  )

}
