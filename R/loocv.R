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
#' @seealso [loo_compare] for comparison of many LOO CV results; you can print and plot
#'          output via [plot.baggr_cv] and [print.baggr_cv]
#' @details
#'
#' The values returned by `loocv()` can be used to understand how excluding
#' any one group affects the overall result, as well as how well the model
#' predicts the omitted group. LOO-CV approaches are a good general practice
#' for comparing Bayesian models, not only in meta-analysis.
#'
#' To learn about cross-validation see Gelman et al 2014.
#'
#' This function automatically runs _K_ baggr models, where _K_ is number of groups (e.g. studies),
#' leaving out one group at a time. For each run, it calculates
#' _expected log predictive density_ (ELPD) for that group (see Gelman et al 2013).
#' (In the logistic model, where the proportion in control group is unknown, each of
#' the groups is divided into data for controls, which is kept for estimation, and data for
#' treated units, which is not used for estimation but only for calculating predictive density.
#' This is akin to fixing the baseline risk and only trying to infer the odds ratio.)
#'
#' The main output is the cross-validation
#' information criterion, or -2 times the ELPD summed over _K_ models.
#' (We sum the terms as we are working with logarithms.)
#' This is related to, and often approximated by, the Watanabe-Akaike
#' Information Criterion. When comparing models, smaller values mean
#' a better fit.
#'
#' For running more computation-intensive models, consider setting the
#' `mc.cores` option before running loocv, e.g. `options(mc.cores = 4)`
#' (by default baggr runs 4 MCMC chains in parallel).
#' As a default, rstan runs "silently" (`refresh=0`).
#' To see sampling progress, please set e.g. `loocv(data, refresh = 500)`.
#'
#' @examples
#' \dontrun{
#' # even simple examples may take a while
#' cv <- loocv(schools, pooling = "partial")
#' print(cv)      # returns the lpd value
#' attributes(cv) # more information is included in the object
#' }
#'
#' @author Witold Wiecek
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @references
#' Gelman, Andrew, Jessica Hwang, and Aki Vehtari.
#' 'Understanding Predictive Information Criteria for Bayesian Models.'
#' Statistics and Computing 24, no. 6 (November 2014): 997–1016.
#' @export

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
  # temp_cores <- FALSE

  # Model with all of data:
  full_fit <- tryCatch(
    baggr(data, refresh = 0, ...),
    error = function(e) {
      stop(
        "Inference failed for the model with all data, please try fitting outside of loocv().\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )})

  # Prepare the arguments
  args[["data"]] <- data
  args[["model"]] <- full_fit$model
  # if(full_fit$model == "rubin_full")
  # stop("LOO CV is not implemented for full data models yet.")
  if(!("prior" %in% names(args))) {
    message("(Prior distributions taken from the model with all data. See $prior.)")
    args[["formatted_prior"]] <- full_fit$formatted_prior
  }

  # Determine number and names of groups
  if(args[["model"]] %in% c("rubin_full", "quantiles", "logit")) {
    #For individual-level data models we need to calculate N groups
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

    # Partitioning data into the fit data and LOO
    if(args[["model"]] %in% c("quantiles")) {
      args$data      <- data[data[[group_col]] != group_names[i], ]
      args$test_data <- data[data[[group_col]] == group_names[i], ]
    } else if(args[["model"]] %in% c("rubin_full", "logit")) {
      trt_column <- ifelse(is.null(args[["treatment"]]), "treatment", args[["treatment"]])
      # for these models we have to use control group in leave out study to calculate baseline
      args$data      <- data[data[[group_col]] != group_names[i] | data[[trt_column]] == 0, ]
      # ...and then use treatment arm of the left out study
      args$test_data <- data[data[[group_col]] == group_names[i] & data[[trt_column]] == 1, ]
    } else {
      args$data      <- data[-i,]
      args$test_data <- data[i,]
    }
    # Done partitioning data.

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

  # if(temp_cores)
  # options(mc.cores = NULL)

  return(out)
}

#' Check if something is a baggr_cv object
#' @param x object to check
is.baggr_cv <- function(x) {
  inherits(x, "baggr_cv")
}

#' Compare LOO CV models
#'
#' Given multiple [loocv] outputs, calculate differences in their expected log
#' predictive density.
#'
#' @param ... A series of `baggr_cv` objects passed as arguments, with a minimum of 2
#'          arguments required for comparison. `baggr_cv` objects can be created via the
#'          [loocv] function. In instances where more than 2 arguments are passed, the
#'          first model will be compared sequentially to all other provided models.
#'          Arguments can be passed with names (see example below).
#' @return  Returns a series of comparisons in order of the arguments provided as Model 1 - Model N for
#'          N loocv objects provided. Model 1 corresponds to the first object passed and
#'          Model N corresponds to the Nth object passed.
#' @export loo_compare
#' @seealso [loocv] for fitting LOO CV objects and explanation of the procedure;
#' loo package by Vehtari et al (available on CRAN) for a more comprehensive approach
#' @examples
#' \dontrun{
#' # 2 models with more/less informative priors -- this will take a while to run
#' cv_1 <- loocv(schools, model = "rubin", pooling = "partial")
#' cv_2 <- loocv(schools, model = "rubin", pooling = "partial",
#'               prior_hypermean = normal(0, 5), prior_hypersd = cauchy(0,2.5))
#' loo_compare("Default prior"=cv_1,"Alternative prior"=cv_2)
#' }
#' @export
loo_compare <- function(...) {
  UseMethod("loo_compare")
}

#' @aliases loo_compare
#' @export
loo_compare.baggr_cv <- function(...) {
  l <- list(...)
  if(all(unlist(lapply(l, class)) == "baggr_cv")) {
    if(is.null(names(l)))
      names(l) <- paste("Model", 1:length(l))
    if(length(unique(names(l))) != length(names(l)))
      stop("You must use unique model names")
    loos <- l
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
    diffname <- paste0(names(loos)[1], " - ", names(loos)[i])
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
#' @param ... additional arguments for s3 consistency
#' @importFrom crayon bold
#' @export
print.compare_baggr_cv <- function(x, digits = 3, ...) {
  mat <- as.matrix(x)
  class(mat) <- "matrix"
  cat(crayon::bold(paste0("Comparison of cross-validation\n\n")))
  print(signif(mat, digits = digits))
  cat("\n")
  cat("Positive ELPD indicates the reference group is preferred.")
}

#' Print baggr cv objects nicely
#' @param x `baggr_cv` object obtained from [loocv] to print
#' @param digits number of digits to print
#' @param ... Unused, ignore
#' @importFrom crayon bold
#' @export
print.baggr_cv <- function(x, digits = 3, ...) {

  mat <- matrix(nrow = 2, ncol = 2)

  mat[1,] <- c(x$elpd, x$se)
  mat[2,] <- c(x$looic, abs(-2*x$se))

  colnames(mat) <- c("Estimate", "Standard Error")
  rownames(mat) <- c("elpd", "looic")

  cat("LOO estimate based on", crayon::bold(paste0(x$K, "-fold cross-validation\n")))
  cat("\n")
  print(signif(mat, digits = digits))

}



#' Plotting method for results of baggr LOO analyses
#'
#' @param x output from [loocv] that has `return_models = TRUE`
#' @param y Unused, ignore
#' @param ... Unused, ignore
#' @param add_values logical; if `TRUE`, values of _elpd_ are printed next to each
#'                   study
#' @return `ggplot2` plot in similar style to [baggr_compare] default plots
#' @export
#'
plot.baggr_cv <- function(x, y, ..., add_values = TRUE){
  loo_model <- x
  input_data <- loo_model$full_model$summary_data
  if(is.null(input_data))
    input_data <- loo_model$full_model$data

  if(is.null(loo_model$models))
    stop("To plot, the loocv() output must include models (return_models = TRUE).")

  mm1 <- do.call(rbind,
                 lapply(loo_model$models, function(x) treatment_effect(x, summary = T)$tau))
  df1 <- data.frame(
    setNames( as.data.frame(mm1[,c(1,4,3)]) , c("lci", "median", "uci")),
    model = "LOOCV estimate",
    elpd = loo_model$pointwise,
    group = input_data$group)
  df2 <- data.frame(lci = input_data$tau - 1.96*input_data$se,
                    median = input_data$tau,
                    uci = input_data$tau + 1.96*input_data$se,
                    group = input_data$group,
                    elpd = NA,
                    model = "Test data (out-of-sample)")
  df <- rbind(df1, df2)
  pl <- single_comp_plot(df)

  fmti <- function(x, digits = values_digits) {
    format(round(x, digits), nsmall = digits)
  }

  # add_values
  if(add_values){
    group <- median <- lci <- uci <- model <- elpd <- NULL
    values_digits <- 2
    values_size <- 2.5

    pl <- pl +
      ggplot2::geom_text(
        aes(label = fmti(elpd, values_digits),
            hjust = 1,
            y = 1.1*max(uci)),
        position = ggplot2::position_dodge(width = .5),
        size = values_size) +
      ggplot2::theme(strip.text.x = ggplot2::element_blank()) +
      ggplot2::coord_flip(clip = "off")
  }

  pl
}
