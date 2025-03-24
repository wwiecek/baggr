#' Convert inputs for baggr models
#'
#' Converts data to a list of inputs suitable for Stan models,
#' checks integrity of data and suggests the appropriate default model
#' if needed. Typically all of this is
#' done automatically by [baggr], so __this function is included only for debugging__
#' or running (custom) models "by hand".
#'
#' @param data `data.frame`` with desired modelling input
#' @param model valid model name used by baggr;
#'              see [baggr] for allowed models
#'              if `model = NULL`, this function will try to find appropriate model
#'              automatically
#' @param covariates Character vector with column names in `data`.
#'                   The corresponding columns are used as
#'                   covariates (fixed effects) in the meta-regression model.
#' @param effect Only matters for binary data, use `logOR`, `logRR`, or `RD`. Otherwise ignore.
#'                See [prepare_ma] for details.
#' @param quantiles vector of quantiles to use (only applicable if `model = "quantiles"`)
#' @param group name of the column with grouping variable
#' @param outcome name of column with outcome variable (designated as string)
#' @param treatment name of column with treatment variable
#' @param test_data same format as `data` argument, gets left aside for
#'                  testing purposes (see [baggr])
#' @param silent Whether to print messages when evaluated
#' @return R structure that's appropriate for use by [baggr] Stan models;
#'         `group_label`, `model`, `effect` and `n_groups` are included as attributes
#'         and are necessary for [baggr] to work correctly
#' @details Typically this function is only called within [baggr] and you do
#'          not need to use it yourself. It can be useful to understand inputs
#'          or to run models which you modified yourself.
#'
#'
#' @author Witold Wiecek
#' @examples
#' # simple meta-analysis example,
#' # this is the formatted input for Stan models in baggr():
#' convert_inputs(schools, "rubin")
#' @export

convert_inputs <- function(data,
                           model,
                           quantiles,
                           effect = NULL,
                           group  = "group",
                           outcome   = "outcome",
                           treatment = "treatment",
                           covariates = c(),
                           test_data = NULL,
                           silent = FALSE) {

  # If lazy users forgot to define their group column,
  # check if the first column is usable
  # group <- find_group_column(data, group)

  # Step 1: check what data are available (with some conversions) -----

  available_data <- detect_input_type(data, group, treatment, outcome)

  if(!is.null(test_data)){
    available_data_test <- detect_input_type(test_data, group, treatment, outcome)
    if(available_data != available_data_test)
      stop("'test_data' is of type ", available_data_test,
           " and 'data' is of type ", available_data)
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
    # model <- names(model_data_types)[which(model_data_types == available_data)[1]]
    model <- data_type_default_model[[available_data]]
    if(!silent)
      message(paste0("Automatically chose ", crayon::bold(model_names[model]),
                     " based on input data."))
  } else {
    if(!(model %in% names(model_data_types)))
      stop("Unrecognised model, can't format data.")
  }

  # Convert mutau data to Rubin model data if requested
  if(model == "rubin" && available_data == "pool_wide"){
    test_data$se <- test_data$se.tau
    test_data$se.tau <- test_data$se.mu <- test_data$mu <- NULL
    data$se <- data$se.tau
    data$se.tau <- data$se.mu <- data$mu <- NULL
    available_data <- "pool_noctrl_narrow"
  }



  # Step 2: check what data are required by the model and match -----
  required_data <- model_data_types[[model]]

  if(required_data == "individual_binary" && available_data == "pool_binary") {
    data <- binary_to_individual(data, group, covariates, FALSE)
    available_data <- "individual_binary"
    message("Data were automatically converted from summary to individual-level.")
  }


  if(model == "rubin" && available_data == "pool_binary"){
    if(is.null(effect) || !(effect %in% c("logOR", "logRR", "RD"))) {
        message('Automatically summarising binary data with logOR.
              In baggr() set effect to one of "logOR", "logRR", "RD".
              Alternatively, use ?prepare_ma to do this manually before running.')
        effect <- "logOR"
    }
    data <- prepare_ma(data, effect = effect, group = group)
    group <- "group"
    available_data <- required_data
  }

  if(required_data != available_data)
    stop(paste(
      "Data provided is of type", data_type_names[available_data],
      "and the model requires", data_type_names[required_data]))
  #for now this means no automatic conversion of individual->pooled


  # Step 3: conversions of data -----

  # 3.1. individual level data
  if(grepl("individual", required_data)) {
    groups        <- factor(as.character(data[[group]]),
                            levels = unique(data[[group]]))
    group_numeric <- as.numeric(groups)
    group_label   <- levels(groups)

    if(!is.null(test_data)) {
      groups_test <- as.factor(as.character(test_data[[group]]))
      group_numeric_test <- as.numeric(groups_test)
      group_label_test <- levels(groups_test)
      if(any(group_label_test %in% group_label) && model != "logit")
        message(
          "Test data has some groups that have same labels as groups in data. ",
          "For cross-validation they will be treated as 'new' groups.")
      if((!all(group_label_test %in% group_label) && model == "logit") || !all(test_data[[treatment]] == 1))
        message(
          "Test data for ", model, " model should include treated units only. ",
          "Baselines for all these groups should be included in data argument.")
    }

    if(model %in% c("rubin_full", "mutau_full", "logit")){
      out <- list(
        K = max(group_numeric),
        N = nrow(data),
        P = 2, #will be dynamic
        y = data[[outcome]],
        treatment = data[[treatment]],
        site = group_numeric
      )
      # }
      # if(model == "logit") {
      if(is.null(test_data)) {
        out$N_test <- 0
        out$K_test <- 0
        out$test_y <- array(0, dim = 0)
        out$test_site <- array(0, dim = 0)
        out$test_treatment <- array(0, dim = 0)
        if(model %in% c("rubin_full", "mutau_full"))
          out$test_sigma_y_k <- array(0, dim = 0)

      } else {
        out$N_test <- nrow(test_data)
        out$K_test <- max(group_numeric_test)
        out$test_y <- test_data[[outcome]]
        out$test_treatment <- test_data[[treatment]]
        out$test_site <- group_numeric_test
        # calculate SEs in each test group
        if(model %in% c("rubin_full", "mutau_full")){
          se_in_each_group <- sapply(
            1:max(group_numeric_test), function(i) {
              n <- sum(group_numeric_test == i)
              sd(test_data[[outcome]][group_numeric_test == i])/sqrt(n)
            })
          if(any(is.na(se_in_each_group)))
            stop("Cannot calculate SE in groups in test data. Each out-of-sample ",
                 "group must be of size at least 2.")
          out$test_sigma_y_k <- array(se_in_each_group, dim = max(group_numeric_test))
        }
      }
    }
    if(model == "quantiles"){
      # This is currently disabled, together with the quantiles model

      # if((any(quantiles < 0)) ||
      #    (any(quantiles > 1)))
      #   stop("quantiles must be between 0 and 1")
      # if(length(quantiles) < 2)
      #   stop("cannot model less than 2 quantiles")
      # data[[group]] <- group_numeric
      # out <- summarise_quantiles_data(data, quantiles, outcome, group, treatment)
      # message("Data have been automatically summarised for quantiles model.")
      #
      # # Fix for R 3.5.1. on Windows
      # # https://stackoverflow.com/questions/51343022/
      # out$temp <- out[["y_0"]]
      # # out$y_0 <- NULL
      # out[["y_0"]] <- out$temp
      # out$temp <- NULL
      #
      # # Cross-validation:
      # if(is.null(test_data)){
      #   out$K_test <- 0
      #   out$test_theta_hat_k <- array(0, dim = 0)
      #   out$test_se_theta_k <- array(0, dim = 0)
      #   out$test_y_0 <- array(0, dim = c(0, ncol(out$y_0)))
      #   out$test_y_1 <- array(0, dim = c(0, ncol(out$y_0)))
      #   out$test_Sigma_y_k_0 <- array(0, dim = c(0, ncol(out$y_0), ncol(out$y_0)))
      #   out$test_Sigma_y_k_1 <- array(0, dim = c(0, ncol(out$y_0), ncol(out$y_0)))
      # } else {
      #   out_test <- summarise_quantiles_data(test_data, quantiles,
      #                                        outcome, group, treatment)
      #   out$K_test <- out_test$K #reminder: K is number of sites,
      #   # N is number of quantiles
      #   out$test_y_0 <- out_test$y_0
      #   out$test_y_1 <- out_test$y_1
      #   out$test_Sigma_y_k_0 <- out_test$Sigma_y_k_0
      #   out$test_Sigma_y_k_1 <- out_test$Sigma_y_k_1
      # }
    }

    if(model == "sslab") {
      # Generic code for dividing observations into positive, negative and == 0 components
      cat <- ifelse(data[[outcome]] < 0, 1, ifelse(data[[outcome]] == 0, 2, 3))
      out <- list(
        K = max(group_numeric),
        N = nrow(data),
        M = 3,
        P = 2, #covariates are taken care of later in this function
        x = array(cbind(1, data[[treatment]]), c(nrow(data), 2)),
        N_neg = sum(cat == 1),
        N_pos = sum(cat == 3),
        y_neg = -1*data[[outcome]][cat == 1],
        site_neg = group_numeric[cat == 1],
        y_pos = data[[outcome]][cat == 3],
        site_pos = group_numeric[cat == 3],
        cat = cat,
        treatment_pos = data[[treatment]][cat == 3],
        treatment_neg = data[[treatment]][cat == 1],
        site = group_numeric
      )

      if(is.null(test_data)){
        # This will have to be done for v0.8 release
        out_test <- list()
      } else {
        out_test <- list()
      }

      for(nm in names(out_test))
        out[[nm]] <- out_test[[nm]]
    }
  }

  # 3.2. summary data: treatment effect only -----
  if(required_data == "pool_noctrl_narrow"){
    group_label <- data[[group]]
    if(is.null(data[[group]]) && (group != "group"))
      warning(paste0("Column '", group,
                     "' does not exist in data. No labels will be added."))
    check_columns_numeric(data[,c("tau", "se")])
    out <- list(
      K = nrow(data),
      theta_hat_k = data[["tau"]],
      se_theta_k = data[["se"]]
    )
    if(is.null(test_data)){
      out$K_test <- 0
      out$test_theta_hat_k <- array(0, dim = 0)
      out$test_se_theta_k <- array(0, dim = 0)
    } else {
      if(is.null(test_data[["tau"]]) ||
         is.null(test_data[["se"]]))
        stop("Test data must be of the same format as input data")
      out$K_test <- nrow(test_data)
      # remember that for 1-dim cases we need to pass array()
      out$test_theta_hat_k <- array(test_data[["tau"]], dim = c(nrow(test_data)))
      out$test_se_theta_k <- array(test_data[["se"]], dim = c(nrow(test_data)))
    }
  }


  # 3.3. summary data: baseline & treatment effect -----
  if(required_data == "pool_wide"){
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
      theta_hat_k = matrix(c(data[["mu"]], data[["tau"]]), 2, nr, byrow = TRUE),
      se_theta_k = matrix(c(data[["se.mu"]], data[["se.tau"]]), 2, nr, byrow = TRUE)
    )
    if(is.null(test_data)){
      out$K_test <- as.integer(0)
      out$test_theta_hat_k <- array(0, dim = c(2,0))
      out$test_se_theta_k <- array(0, dim = c(2,0))
    } else {
      if(is.null(test_data[["mu"]]) ||
         is.null(test_data[["tau"]]) ||
         is.null(test_data[["se.mu"]]) ||
         is.null(test_data[["se.tau"]]))
        stop("Test data must be of the same format as input data")
      out$K_test <- nrow(test_data)
      out$test_theta_hat_k <- matrix(c(test_data[["mu"]], test_data[["tau"]]),
                                     2, nrow(test_data), byrow = TRUE)
      out$test_se_theta_k <- matrix(c(test_data[["se.mu"]], test_data[["se.tau"]]),
                                    2, nrow(test_data), byrow = TRUE)
    }
  }

  # 4. Include covariates ------
  # if(required_data != "individual") {
  if(length(covariates) > 0) {
    if(model == "quantiles")
      stop("Quantiles model cannot regress on covariates.")

    if(!all(covariates %in% names(data)))
      stop(paste0("Covariates ",
                  paste(covariates[!(covariates %in% names(data))], collapse=","),
                  " are not columns in input data"))

    for(cov in covariates)
      if(any(is.na(data[[cov]])))
        stop("NA values present in covariates")

    # Test_data preparation
    cov_bind <- data[,covariates, drop = FALSE]
    if(!is.null(test_data)){
      cov_bind <- try(rbind(
        cov_bind,
        test_data[,covariates, drop = FALSE]))
      if(inherits(cov_bind, "try-error"))
        stop("Cannot bind data and test_data. Ensure that all ",
             "covariates are present and same levels are used.")
    }
    cov_bind$tau <- 0
    cov_bind[] <- lapply(cov_bind, function(x) if(is.character(x)) factor(x) else x)
    cov_mm <- model.matrix(as.formula(
      paste("tau ~", paste(covariates, collapse="+"))),
      data=cov_bind)
    out$X <- cov_mm[1:nrow(data), 2:ncol(cov_mm), drop = FALSE]
    out$Nc <- ncol(out$X)
    if(!is.null(test_data))
      out$X_test <- cov_mm[(nrow(data)+1):nrow(cov_mm), 2:ncol(cov_mm), drop = FALSE]
    else
      out$X_test <- array(0, dim=c(0, out$Nc))

    covariate_coding <- colnames(out$X)
    covariate_levels <- lapply(cov_bind, levels)
    covariate_levels[["tau"]] <- NULL

  } else {
    covariate_coding <- c()
    covariate_levels <- c()
    out$Nc <- 0
    if(model != "quantiles"){
      out$X <- array(0, dim=c(nrow(data), 0))
      out$X_test <- array(0, dim=c(ifelse(is.null(test_data), 0, nrow(test_data)), 0))
    }
  }

  na_cols <- unlist(lapply(out, function(x) any(is.na(x))))
  if(any(na_cols))
    stop(paste0("baggr() does not allow NA values in inputs (see vectors ",
                paste(names(out)[na_cols], collapse = ", "), ")"))

  # When using data frames with 1 rows, we need to explicitly define them
  # as arrays before passing from R to Stan. So here I check if any of them
  # are length 1 and change them to arrays
  numeric_to_array_c <- c("theta_hat_k", "se_theta_k")
  for(nm in numeric_to_array_c) {
    if(length(out[[nm]]) == 1)
      out[[nm]] <- array(out[[nm]], dim = 1)
  }

  return(structure(
    out,
    data_type = available_data,
    data = data,
    covariate_coding = covariate_coding,
    covariate_levels = covariate_levels,
    group_label = group_label,
    n_groups = out[["K"]],
    model = model,
    effect = effect))
}
