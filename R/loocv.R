#' Leave one out cross-validation for \code{baggr} models
#'
#' @param data Input data frame - same as for \code{baggr()} function.
#' @param return_models logical; if FALSE, summary statistics will be returned and the
#'                      models discarded;
#'                      if TRUE, a list of models will be returned alongside summaries
#' @param ... Additional arguments passed to \code{baggr()}.
#' @return A data.frame of results and list of baggr models
#'
#' @details
#' For each model we calculate log predictive density l.p.d. (TBC)
#'
#' @examples
#' #even simple examples may take a long moment
#' loocv(schools, pooling = "partial")
#' loocv(schools, pooling = "full")
#'
#' @author Witold Wiecek, Rachael Meager
#' @references Gelman A, Hwang J, Vehtari A.
#'             Understanding predictive information criteria for Bayesian models.
#'             Statistics and Computing. 2014 Nov 24(6):997-1016.
#'             \href{http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf}{(PDF link)}
#'
#' @export
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'

loocv <- function(data, return_models = FALSE, ...) {
  K <- nrow(data)
  if(K > 50)
    message(paste("Large number of groups - note that the
                   model will be repeated", K, "times"))

  # Set prior, if not specified by the user
  args <- list(...)
  # Stop if pooling == "none"
  if(!is.null(args[["pooling"]]) && (args[["pooling"]] == "none"))
    stop("For a model with no pooling LOO CV doesn't exist.")


  # Model with all of data:
  full_fit <- baggr(data, ...)

  # Prepare the arguments
  args[["data"]] <- data
  args[["model"]] <- full_fit$model
  if(!("prior" %in% names(args))) {
    message("(Prior distributions taken from the model with all data. See $prior.)")
    args[["prior"]] <- full_fit$prior
  }

  # LOO CV models
  cat("\n")
  pb <- utils::txtProgressBar(style = 3)
  # should sink() Stan print()'s?
  kfits <- lapply(as.list(1:K), function(i) {
    # baggr(data = data[-i,], test_data = data[i,], ...)
    args$data <- data[-i,]
    args$test_data <- data[i,]
    utils::setTxtProgressBar(pb, (i-1)/K)
    do.call(baggr, args)
  })
  close(pb)

  tau_estimate <-
    lapply(kfits, function(x) apply(as.matrix(x$fit, "tau"), 2, mean))
  loglik <-
    lapply(kfits, function(x) apply(as.matrix(x$fit, "logpd"), 2, mean))

  out <- structure(
    -2*sum(unlist(loglik)),
    df = data.frame(
      "tau" = unlist(tau_estimate),
      "lpd" = unlist(loglik)),
    full_model = full_fit,
    prior = args[["prior"]]
  )
  if(return_models)
    out$models <- kfits

  class(out) <- "baggr_cv"
  return(out)
}

print.baggr_cv <- function(x, ...) {
  attributes(x) <- NULL
  print(x)
}
