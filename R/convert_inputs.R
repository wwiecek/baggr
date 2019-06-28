#' @title Convert inputs for baggr models
#'
#' Converts data to Stan inputs, checks intergrity of data
#' and suggests default model if needed.
#'
#' @param data `data.frame`` with desired modelling input
#' @param model valid model name used by baggr;
#'              see [baggr()] for allowed models
#'              if `model = NULL`, this function will try to find appropriate model
#'              automatically
#' @param quantiles vector of quantiles to use (only applicable if `model = "quantiles"`)
#' @param group name of the column with grouping variable
#' @param outcome name of column with outcome variable
#' @param treatment name of column with treatment variable
#' @param test_data same format as `data` argument, gets left aside for
#'                  testing purposes (see [baggr()] and [loocv()])
#' @return R structure that's appropriate for use by [baggr()] Stan models;
#'         `group_label`, `model` and `n_groups` are incuded as attributes
#'         and are necessary for [baggr()] to work correctly
#' @details Typically this function is only called within [baggr()] and you do
#'          not need to use it yourself. It can be useful to understand inputs
#'          or to run models which you modified yourself.
#'
#'          For quantile models, [summarise_quantiles_data()] is called from within
#'          this function.
#'
#' @author Witold Wiecek, Rachael Meager
#' @export

convert_inputs <- function(data,
                           model,
                           quantiles,
                           group  = "group",
                           outcome   = "outcome",
                           treatment = "treatment",
                           test_data = NULL) {

  # check what kind of data is required for the model & what's available
  model_data_types <- c("rubin" = "pool_noctrl_narrow",
                        "mutau" = "pool_wide",
                        "full" = "individual",
                        "quantiles" = "individual") #for now no quantiles model from summary level data
  available_data <- detect_input_type(data, group)

  # if(available_data == "unknown")
  # stop("Cannot automatically determine type of input data.")
  # let's assume data is individual-level
  # if we can't determine it
  # because it may have custom columns
  if(available_data == "unknown")
    available_data <- "individual" #in future can call it 'inferred ind.'

  if(available_data == "individual")
    check_columns(data, outcome, group, treatment)


  if(is.null(model)) {
    message("Attempting to infer the correct model for data.")
    # we take FIRST MODEL THAT SUITS OUR DATA!
    model <- names(model_data_types)[which(model_data_types == available_data)[1]]
    message(paste0("Chosen model ", model))
  } else {
    if(!(model %in% names(model_data_types)))
      stop("Unrecognised model, can't format data.")
  }

  required_data <- model_data_types[[model]]


  if(required_data != available_data)
    stop(paste(
      "Data provided is of type", available_data,
      "and the model requires", required_data))
  #for now this means no automatic conversion of individual->pooled

  # individual level data -----
  if(required_data == "individual"){
    # # check correctness of inputs:
    # (This check moved up now.)
    # if(is.null(data[[group]]))
    #   stop("No 'group' column in data.")
    # if(is.null(data[[outcome]]))
    #   stop("No outcome column in data.")
    # if(is.null(data[[treatment]]))
    #   stop("No treatment column in data.")

    groups <- as.factor(as.character(data[[group]]))
    group_numeric <- as.numeric(groups)
    group_label <- levels(groups)

    if(model == "full")
      out <- list(
        K = max(group_numeric),
        N = nrow(data),
        P = 2, #will be dynamic
        y = data[[outcome]],
        ITT = data[[treatment]],
        site = group_numeric
      )
    if(model == "quantiles"){
      if((any(quantiles < 0)) ||
         (any(quantiles > 1)))
        stop("quantiles must be between 0 and 1")
      if(length(quantiles) < 2)
        stop("cannot model less then 2 quantiles")
      data[[group]] <- group_numeric
      out <- summarise_quantiles_data(data, quantiles,
                                      outcome, group, treatment)
      message("Data have been automatically summarised for quantiles model.")

      # Fix for R 3.5.1. on Windows
      # https://stackoverflow.com/questions/51343022/
      out$temp <- out[["y_0"]]
      out$y_0 <- NULL
      out[["y_0"]] <- out$temp
      out$temp <- NULL
    }

    # Cross-validation:
    if(is.null(test_data)){
      out$K_test <- 0
      out$test_tau_hat_k <- array(0, dim = 0)
      out$test_se_k <- array(0, dim = 0)
      out$test_y_0 <- array(0, dim = c(0, ncol(out$y_0)))
      out$test_y_1 <- array(0, dim = c(0, ncol(out$y_0)))
      out$test_Sigma_y_k_0 <- array(0, dim = c(0, ncol(out$y_0), ncol(out$y_0)))
      out$test_Sigma_y_k_1 <- array(0, dim = c(0, ncol(out$y_0), ncol(out$y_0)))
    } else {
      out_test <- summarise_quantiles_data(test_data, quantiles,
                                      outcome, group, treatment)
      out$K_test <- out_test$K #reminder: K is number of sites, N is number of quantiles
      out$test_y_0 <- out_test$y_0
      out$test_y_1 <- out_test$y_1
      out$test_Sigma_y_k_0 <- out_test$Sigma_y_k_0
      out$test_Sigma_y_k_1 <- out_test$Sigma_y_k_1
    }
  }

  # summary data: treatment effect only -----
  if(required_data == "pool_noctrl_narrow"){
    group_label <- data[[group]]
    out <- list(
      K = nrow(data),
      tau_hat_k = data[["tau"]],
      se_tau_k = data[["se"]]
    )
    if(is.null(test_data)){
      out$K_test <- 0
      out$test_tau_hat_k <- array(0, dim = 0)
      out$test_se_k <- array(0, dim = 0)
    } else {
      if(is.null(test_data[["tau"]]) ||
         is.null(test_data[["se"]]))
        stop("Test data must be of the same format as input data")
      out$K_test <- nrow(test_data)
      # remember that for 1-dim cases we need to pass array()
      out$test_tau_hat_k <- array(test_data[["tau"]], dim = c(nrow(test_data)))
      out$test_se_k <- array(test_data[["se"]], dim = c(nrow(test_data)))
    }
  }


  # summary data: baseline & treatment effect -----
  if(required_data == "pool_wide"){
    group_label <- data[[group]]
    nr <- nrow(data)
    out <- list(
      K = nr,
      P = 2, #fixed for this case
      # Remember, first row is always mu (baseline), second row is tau (effect)
      # (Has to be consistent against ordering of prior values.)
      tau_hat_k = matrix(c(data[["mu"]], data[["tau"]]), 2, nr, byrow = T),
      se_tau_k = matrix(c(data[["se.mu"]], data[["se.tau"]]), 2, nr, byrow = T)
    )
    if(is.null(test_data)){
      out$K_test <- 0
      out$test_tau_hat_k <- array(0, dim = c(2,0))
      out$test_se_k <- array(0, dim = c(2,0))
    } else {
      if(is.null(test_data[["mu"]]) ||
         is.null(test_data[["tau"]]) ||
         is.null(test_data[["se.mu"]]) ||
         is.null(test_data[["se.tau"]]))
        stop("Test data must be of the same format as input data")
      out$K_test <- nrow(test_data)
      out$test_tau_hat_k <- matrix(c(test_data[["mu"]], test_data[["tau"]]),
                                   2, nrow(test_data), byrow = T)
      out$test_se_k <- matrix(c(test_data[["se.mu"]], test_data[["se.tau"]]),
                              2, nrow(test_data), byrow = T)
    }
  }

  na_cols <- unlist(lapply(out, function(x) any(is.na(x))))
  if(any(na_cols))
    stop(paste0("baggr() does not allow NA values in inputs (see vectors ",
                paste(names(out)[na_cols], collapse = ", "), ")"))


  return(structure(
    out,
    group_label = group_label,
    n_groups = out[["K"]],
    model = model))

}
