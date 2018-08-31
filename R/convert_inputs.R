#' @title Convert inputs for baggr Stan models
#'
#' @description
#' Allows conversions from full to reduced (summary) data and
#' in summary between long (control and treatment on separate rows)
#' and wide data (control and treatment in separate columns).
#'
#' @param data data.frame with desired modelling input
#' @param model valid model name used by baggr;
#'              see \code{?baggr} for allowed models
#' @param grouping name of the column with grouping variable
#' @param outcome name of column with outcome variable
#' @param treatment name of column with treatment variable
#' @param standardise logical; whether to standardise data when converting
#' @return data.frame of class baggr_data that baggr() uses
#' @details
#' The conversions will typically happen automatically when data is fed to baggr()
#' function. This function can be used to explicitly convert from full to reduced
#' data without analysing it in a model.
#' @author Witold Wiecek
#' @export

convert_inputs <- function(data,
                           model,
                           grouping  = "site",
                           outcome   = "outcome",
                           treatment = "treatment",
                           standardise = FALSE,
                           test_data = NULL) {

  # check what kind of data is required for the model & what's available
  model_data_types <- c("rubin" = "pool_noctrl_narrow",
                        "mutau" = "pool_wide",
                        "full" = "individual")
  available_data <- detect_input_type(data, grouping)

  if(available_data == "individual") {
    # stop if variables are not available
    if(is.null(data[[grouping]]))
      stop("No grouping column in data.")
    if(is.null(data[[outcome]]))
      stop("No outcome column in data.")
    if(is.null(data[[treatment]]))
      stop("No treatment column in data.")
  }

  # if(available_data == "unknown")
    # stop("Cannot automatically determine type of input data.")
  # let's assume data is individual-level
  # if we can't determine it
  # because it may have custom columns
  if(available_data == "unknown")
    available_data <- "individual" #in future call it 'inferred ind.'

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
    # check correctness of inputs:
    if(is.null(data[[grouping]]))
      stop("No 'site' column in data.")
    if(is.null(data[[outcome]]))
      stop("No outcome column in data.")
    if(is.null(data[[treatment]]))
      stop("No treatment column in data.")

    site_numeric <- as.numeric(as.factor(as.character(data[[grouping]])))
    site_label <- unique(as.character(data[[grouping]]))

    if(standardise)
      data[[outcome]] <- as.vector(scale(data[[outcome]]))

    out <- list(
      K = max(site_numeric),
      N = nrow(data),
      P = 2, #will be dynamic
      y = data[[outcome]],
      ITT = data[[treatment]],
      site = site_numeric
    )
  }

  # summary data: treatment effect only -----
  if(required_data == "pool_noctrl_narrow"){
    site_label <- data[[grouping]]
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
      out$test_tau_hat_k <- test_data[["tau"]]
      out$test_se_k <- test_data[["se"]]
    }
  }


  # summary data: baseline & treatment effect -----
  if(required_data == "pool_wide"){
    site_label <- data[[grouping]]
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
      out$test_tau_hat_k <- array(0, dim = 0)
      out$test_se_tau_k <- array(0, dim = 0)
    } else {
      if(is.null(test_data[["tau"]]) ||
         is.null(test_data[["se"]]))
        stop("Test data must be of the same format as input data")
      out$K_test <- nrow(test_data)
      out$test_tau_hat_k <- test_data[["tau"]]
      out$test_se_tau_k <- test_data[["se"]]
    }
  }


  return(structure(
    out,
    standardised = standardise,
    site_label = site_label,
    n_sites = out[["K"]],
    model = model))

}
