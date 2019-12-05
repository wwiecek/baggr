#' @title Convert from individual to summary data in meta-analyses
#'
#' @description Allows one-way conversion from full to summary data.
#'              Input must be pre-formatted appropriately.
#'
#' @param data data.frame of individual-level observations
#'             with columns for outcome (numeric), treatment (values 0 and 1) and
#'             group (numeric, character or factor);
#'             column names can be user-defined (see below)
#' @param log logical; log-transform the outcome variable?
#' @param rare_event_correction This correction is used when working with
#'             binary data only. The value of correction is added to all arms
#'             in trials where some arms had 0 events.
#'             Using corrections may bias results but is the only alternative to
#'             avoid infinite values.
#' @param cfb logical; calculate change from baseline? If yes, the outcome
#'            variable is taken as a difference between values in `outcome` and
#'            `baseline` columns
#' @param summarise logical; convert to aggregate level data?
#' @param group name of the column with grouping variable
#' @param outcome name of column with outcome variable
#' @param treatment name of column with treatment variable
#' @param baseline name of column with baseline variable
#'
#' @return data.frame with columns \code{mu}, \code{se.mu},
#'         \code{tau} and \code{se.tau}
#'
#' @details
#' The conversions done by this function are not typically needed and may happen automatically
#' when data is fed to [baggr]. However, this function can be used to explicitly
#' convert from full to reduced (summarised) data without analysing it in any model.
#' It can be useful for examining your data.
#'
#' If multiple operations are performed, they are taken in this order:
#' 1) conversion to log scale,
#' 2) calculating change from baseline,
#' 3) summarising data.
#'
#' @author Witold Wiecek
#' @seealso [convert_inputs] for how data is converted into Stan inputs;
#' @export
#' @import stats
#' @import dplyr
#'

prepare_ma <- function(data, #standardise = NULL,
                       effect = c("mean", "logOR", "logRR"),
                       rare_event_correction = 0.25,
                       log = FALSE, cfb = FALSE, summarise = TRUE,
                       treatment="treatment",
                       baseline = NULL,
                       group="group",
                       outcome="outcome") {
  if(grepl("pool", detect_input_type(data, group, treatment)))
    stop("Data must be individual-level to use prepare_ma.")
  check_columns(data, outcome, group, treatment)


  # Input checks and prep
  data <- data[,c(treatment, group, outcome, baseline)]
  if(is.null(baseline))
    names(data) <- c("treatment", "group", "outcome")
  else
    names(data) <- c("treatment", "group", "outcome", "baseline")

  if(any(!stats::complete.cases(data))){
    if(summarise)
      warning("NA values present in data - they were dropped when summarising")
    else
      warning("NA values present in data")
  }

  effect <- match.arg(effect, c("mean", "logOR", "logRR"))
  if(effect %in% c("logOR", "logRR")) {
    if(!is_binary(data$outcome))
      stop("Outcome column is not binary (only 0 and 1 values allowed).")
  }

  # 1. transform data (for now only log)
  if(log) {
    data$outcome  <- log(data$outcome)
    if(!is.null(baseline))
      data$baseline  <- log(data$baseline)
  }

  # 2. Change from baseline
  if(cfb){
    if(is.null(baseline))
      stop("Define baseline column to calculate change from baseline")
    data$outcome <- data$outcome - data$baseline
    data$baseline <- NULL
  }

  # 3. Standardise
  standardise <- NULL #until standardisation is fixed, we cauterise this
  if(!is.null(standardise)) {
    # Whole sample
    if(standardise == "all")
      data$outcome <- (data$outcome - mean(data$outcome)) / sd(data$outcome)
    else if(standardise == "bsl"){
      if(cfb)
        stop("Can't standardise by baseline value if outcome is change to baseline.")
      data$outcome <- (data$outcome - mean(data$baseline)) / sd(data$baseline)
    }
    # In each group separately
    else if(standardise %in% c("group", "by group", "by_group")){
      agg <- stats::aggregate(outcome ~ group,
                              function(x) {c(mean=mean(x), sd=sd(x))},
                              data = data)
      means <- agg$outcome[,"mean"]
      sds <- agg$outcome[,"sd"]
      names(means) <- names(sds) <- agg$group
      data$outcome <- (data$outcome - means[data$group]) / sds[data$group]
    }
    else if(standardise %in% c("group_bsl")){
      if(cfb)
        stop("Can't standardise by baseline value if outcome is change to baseline.")
      if(is.null(baseline))
        stop("Can't standardise by baseline value if baseline=NULL.")
      agg <- stats::aggregate(baseline ~ group,
                              function(x) {c(mean=mean(x), sd=sd(x))},
                              data = data)
      means <- agg$baseline[,"mean"]
      sds <- agg$baseline[,"sd"]
      names(means) <- names(sds) <- agg$group
      data$outcome <- (data$outcome - means[data$group]) / sds[data$group]
    } else if(standardise %in% c("control")) {
      mean_ctrl <- mean(data$outcome[data$treatment == 0])
      sd_ctrl <- sd(data$outcome[data$treatment == 0])
      data$outcome <- (data$outcome -  mean_ctrl)/ sd_ctrl

    } else if(standardise %in% c("group_control")) {
      agg <- stats::aggregate(outcome ~ group,
                              function(x) {c(mean=mean(x), sd=sd(x))},
                              data = data[data$treatment == 0,])
      means <- agg$outcome[,"mean"]
      sds <- agg$outcome[,"sd"]
      names(means) <- names(sds) <- agg$group
      data$outcome <- (data$outcome - means[data$group]) / sds[data$group]

    } else
      stop("Wrong standardise argument, expecting 'all', 'bsl', 'group' or 'group_bsl'")

  }

  # 4. Summarising
  if(summarise){
    if(effect == "mean") {
      magg  <- stats::aggregate(outcome ~ treatment + group,
                                mean, data = data)
      seagg <- stats::aggregate(outcome ~ treatment + group,
                                function(x) sd(x)/sqrt(length(x)), data = data)
      mwide <- stats::reshape(data = magg, timevar = "treatment",
                              idvar = "group", direction = "wide")
      sewide <- stats::reshape(data = seagg, timevar = "treatment",
                               idvar = "group", direction = "wide")
      out <- data.frame(group = mwide$group,
                        mu = mwide$outcome.0,
                        tau = mwide$outcome.1 - mwide$outcome.0,
                        se.mu = sewide$outcome.0,
                        se.tau = sqrt(sewide$outcome.0^2 + sewide$outcome.1^2))
    }

    # Prepare event counts for binary data models
    # (including rare event corrections)
    v <- rare_event_correction
    binary_data_table <- dplyr::group_by(data, group) %>%
      dplyr::summarise(a = sum(outcome[treatment == 1]),
                       n1 = sum(treatment == 1),
                       c = sum(outcome[treatment == 0]),
                       n2 = sum(treatment == 0)) %>%
      dplyr::mutate(b = n1 - a,
                    d = n2 - c) %>%
    # Rare-event correction
      dplyr::mutate(rare = (a == 0 | b == 0 | c == 0 | d == 0)) %>%
      dplyr::mutate(a = v*rare + a,
                    b = v*rare + b,
                    c = v*rare + c,
                    d = v*rare + d) %>%
      dplyr::select(-rare)


    if(effect == "logRR") {
      out <- dplyr::mutate(binary_data_table,
                           tau = log((a/(a+b))/(c/(c+d))),
                           se  = sqrt(1/a + 1/c - 1/(a+b) - 1/(c+d)))
    }
    if(effect == "logOR") {
      out <- dplyr::mutate(binary_data_table,
                           tau = log((a*d)/(b*c)),
                           se  = sqrt(1/a + 1/b + 1/c + 1/d))
    }
  } else {
    out <- data
  }

  out
}

