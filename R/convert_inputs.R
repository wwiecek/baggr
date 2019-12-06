#' Convert inputs for baggr models
#'
#' Converts data to Stan inputs, checks integrity of data
#' and suggests default model if needed. Typically used
#' automatically by [baggr] but useful for debugging.
#'
#' @param data `data.frame`` with desired modelling input
#' @param model valid model name used by baggr;
#'              see [baggr] for allowed models
#'              if `model = NULL`, this function will try to find appropriate model
#'              automatically
#' @param quantiles vector of quantiles to use (only applicable if `model = "quantiles"`)
#' @param group name of the column with grouping variable
#' @param outcome name of column with outcome variable (designated as string)
#' @param treatment name of column with treatment variable
#' @param test_data same format as `data` argument, gets left aside for
#'                  testing purposes (see [baggr])
#' @return R structure that's appropriate for use by [baggr] Stan models;
#'         `group_label`, `model` and `n_groups` are included as attributes
#'         and are necessary for [baggr] to work correctly
#' @details Typically this function is only called within [baggr] and you do
#'          not need to use it yourself. It can be useful to understand inputs
#'          or to run models which you modified yourself.
#'
#'
#' @author Witold Wiecek
#' @examples
#' # simple meta-analysis example:
#' convert_inputs(schools, "rubin")
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
                        "logit" = "individual_binary",
                        "full" = "individual",
                        #for now no quantiles model from summary level data
                        "quantiles" = "individual")
  data_type_names <- c("pool_noctrl_narrow" = "Aggregate (effects only)",
                       "pool_wide" = "Aggregate (control and effects)",
                       "individual" = "Individual-level with continuous outcome",
                       "individual_binary" = "Individual-level with binary outcome")

  available_data <- detect_input_type(data, group, treatment, outcome)

  if(!is.null(test_data)){
    available_data_test <- detect_input_type(test_data, group, treatment, outcome)
    if(available_data != available_data_test)
      stop("'test_data' is of type ", available_data_test, " and 'data' is of type ", available_data)
  }
  # if(available_data == "unknown")
  # stop("Cannot automatically determine type of input data.")
  # let's assume data is individual-level
  # if we can't determine it
  # because it may have custom columns
  if(available_data == "unknown")
    available_data <- "individual" #in future can call it 'inferred ind.'

  if(grepl("individual", available_data)){
    check_columns(data, outcome, group, treatment)
    if(!is.null(test_data))
      check_columns(data, outcome, group, treatment)
  }
  if(is.null(model)) {
    message("Attempting to infer the correct model for data.")
    # we take FIRST MODEL THAT SUITS OUR DATA!
    model <- names(model_data_types)[which(model_data_types == available_data)[1]]
    message(paste0("Chosen model ", model))
  } else {
    if(!(model %in% names(model_data_types)))
      stop("Unrecognised model, can't format data.")
  }

  # Convert mutau data to Rubin model data if requested
  # (For now disabled as it won't work with functions that then
  #  re-use the input data)
  # if(model == "rubin" && available_data == "pool_wide"){
  #   data$se <- data$se.tau
  #   data$se.tau <- data$se.mu <- data$mu <- NULL
  #   available_data <- "pool_noctrl_narrow"
  # }

  required_data <- model_data_types[[model]]


  if(required_data != available_data)
    stop(paste(
      "Data provided is of type", data_type_names[available_data],
      "and the model requires", data_type_names[required_data]))
  #for now this means no automatic conversion of individual->pooled

  # individual level data -----
  if(grepl("individual", required_data)) {

    groups <- as.factor(as.character(data[[group]]))
    group_numeric <- as.numeric(groups)
    group_label <- levels(groups)

    if(!is.null(test_data)) {
      groups_test <- as.factor(as.character(test_data[[group]]))
      group_numeric_test <- as.numeric(groups_test)
      group_label_test <- levels(groups_test)
      if(any(group_label_test %in% group_label))
        message(
          "Test data has some groups that have same labels as groups in data. ",
          "For cross-validation they will be treated as 'new' groups.")
    }

    if(model %in% c("full", "logit")){
      out <- list(
        K = max(group_numeric),
        N = nrow(data),
        P = 2, #will be dynamic
        y = data[[outcome]],
        treatment = data[[treatment]],
        site = group_numeric
      )
    }
    if(model == "logit") {
      if(is.null(test_data)) {
        out$N_test <- 0
        out$K_test <- 0
        out$test_y <- array(0, dim = 0)
        out$test_site <- array(0, dim = 0)
        out$test_treatment <- array(0, dim = 0)
      } else {
        out$N_test <- nrow(test_data)
        out$K_test <- max(group_numeric_test)
        out$test_y <- test_data[[outcome]]
        out$test_treatment <- test_data[[treatment]]
        out$test_site <- group_numeric_test
      }
    }
    if(model == "quantiles"){
      if((any(quantiles < 0)) ||
         (any(quantiles > 1)))
        stop("quantiles must be between 0 and 1")
      if(length(quantiles) < 2)
        stop("cannot model less then 2 quantiles")
      data[[group]] <- group_numeric
      # out <- summarise_quantiles_data(data, quantiles,
      # outcome, group, treatment)
      message("Data have been automatically summarised for quantiles model.")

      # Fix for R 3.5.1. on Windows
      # https://stackoverflow.com/questions/51343022/
      out$temp <- out[["y_0"]]
      out$y_0 <- NULL
      out[["y_0"]] <- out$temp
      out$temp <- NULL

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
        # Disabled until summarise_quantiles_data() gets included
        # in the release again.
        # out_test <- summarise_quantiles_data(test_data, quantiles,
        # outcome, group, treatment)
        # out$K_test <- out_test$K #reminder: K is number of sites,
        # N is number of quantiles
        # out$test_y_0 <- out_test$y_0
        # out$test_y_1 <- out_test$y_1
        # out$test_Sigma_y_k_0 <- out_test$Sigma_y_k_0
        # out$test_Sigma_y_k_1 <- out_test$Sigma_y_k_1
      }
    }
  }

  # summary data: treatment effect only -----
  if(required_data == "pool_noctrl_narrow"){
    group_label <- data[[group]]
    if(is.null(data[[group]]) && (group != "group"))
      warning(paste0("Column '", group,
                     "' does not exist in data. No labels will be added."))
    check_columns_numeric(data[,c("tau", "se")])
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
    group_label <- data[[group]]
    if(is.null(data[[group]]) && (group != "group"))
      warning(paste0("Column '", group,
                     "' does not exist in data. No labels will be added."))
    check_columns_numeric(data[,c("tau", "se.tau", "mu", "se.mu")])
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
    data_type = available_data,
    group_label = group_label,
    n_groups = out[["K"]],
    model = model))

}
