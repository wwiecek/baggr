#' @title Convert from individual to summary data in meta-analyses
#'
#' @description Allows for one-way conversion from full to summary data or for calculation of effects for binary data.
#'              Input must be pre-formatted appropriately.
#'
#' @param data __either__ a data.frame of individual-level observations
#'             with columns for outcome (numeric), treatment (values 0 and 1) and
#'             group (numeric, character or factor); __or__, a data frame with binary data
#'             (must have columns `a`, `c`, `b`/`n1`, `d`/`n2`).
#' @param effect what effect to calculate? a `mean` (and SE) of outcome in groups or
#'               (for binary data) `logOR` (odds ratio), `logRR` (risk ratio);
#' @param log logical; log-transform the outcome variable?
#' @param rare_event_correction This correction is used when working with
#'             binary data (effect `logOR` or `logRR`)
#'             The value of correction is added to all cells in
#'             either some or all rows (groups), depending on `correction_type`.
#'             Using corrections may bias results but is the only alternative to
#'             avoid infinite values.
#' @param correction_type If `"single"` then rare event correction is only applied to the
#'             particular rows that have 0 cells, if `"all"`, then to all studies
#' @param cfb logical; calculate change from baseline? If yes, the outcome
#'            variable is taken as a difference between values in `outcome` and
#'            `baseline` columns
#' @param summarise logical; `TRUE` by default, but you can disable it to obtain
#'                  converted (e.g. logged) data with columns renamed
#' @param group name of the column with grouping variable
#' @param outcome name of column with outcome variable
#' @param treatment name of column with treatment variable
#' @param baseline name of column with baseline variable
#'
#' @return
#' * If you `summarise`: a data.frame with columns for `group`, `tau` and `se.tau`
#'         (for `effect = "mean"`, also baseline means; for `"logRR"` or `"logOR"` also
#'         `a`, `b`, `c`, `d`, which correspond to typical contingency table notation, that is:
#'         `a` = events in exposed; `b` = no events in exposed, `c` = events in unexposed,
#'         `d` = no events in unexposed).
#' * If you do not summarise data, individual level data will be returned, but
#'   some columns may be renamed or transformed (see the arguments above).
#'
#' @details
#' The conversions done by this function are not typically needed and may happen automatically
#' when `data` is given to [baggr]. However, this function can be used to explicitly
#' convert from full to reduced (summarised) data without analysing it in any model.
#' It can be useful for examining your data and generating summary tables.
#'
#' If multiple operations are performed, they are taken in this order:
#' 1) conversion to log scale,
#' 2) calculating change from baseline,
#' 3) summarising data (using appropriate `effect`)
#'
#' @author Witold Wiecek
#' @seealso [convert_inputs] for how any type of data is (internally) converted into
#'          a list of Stan inputs; vignette `baggr_binary` for more details about
#'          rare event corrections
#' @export
#' @import stats
#'
#' @examples
#'
#' # Example of working with binary outcomes data
#' # Make up some individual-level data first:
#' df_rare <- data.frame(group = paste("Study", LETTERS[1:5]),
#'                       a = c(0, 2, 1, 3, 1), c = c(2, 2, 3, 3, 5),
#'                       n1i = c(120, 300, 110, 250, 95),
#'                       n2i = c(120, 300, 110, 250, 95))
#' df_rare_ind <- binary_to_individual(df_rare)
#' # Calculate ORs; default rare event correction will be applied
#' prepare_ma(df_rare_ind, effect = "logOR")
#' # Add 0.5 to all rows
#' prepare_ma(df_rare_ind, effect = "logOR",
#'            correction_type = "all",
#'            rare_event_correction = 0.5)

prepare_ma <- function(data, #standardise = NULL,
                       effect = c("mean", "logOR", "logRR"),
                       rare_event_correction = 0.25,
                       correction_type = c("single", "all"),
                       log = FALSE, cfb = FALSE, summarise = TRUE,
                       treatment="treatment",
                       baseline = NULL,
                       group="group",
                       outcome="outcome") {

  effect <- match.arg(effect)
  correction_type <- match.arg(correction_type)

  if(grepl("pool|unknown", detect_input_type(data, group, treatment, outcome))){
    if(effect %in% c("logOR", "logRR")){
      check_columns_binary(data)
      data <- binary_to_individual(data, group)
      group <- "group"
    }else
      stop("Data must be individual-level (if summarising) or binary data (if converting), see ?prepare_ma")
  }

  check_columns(data, outcome, group, treatment, stop.for.na = FALSE)


  # Input checks and prep
  data <- data[,c(treatment, group, outcome, baseline)]
  if(is.null(baseline))
    names(data) <- c("treatment", "group", "outcome")
  else
    names(data) <- c("treatment", "group", "outcome", "baseline")

  if(any(!stats::complete.cases(data))){
    if(summarise)
      warning("NA values present in data - they may be dropped when summarising")
    else
      check_columns(data, outcome, group, treatment, stop.for.na = TRUE)
  }

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
      magg   <- stats::aggregate(outcome ~ treatment + group,
                                 mean, data = data)
      seagg  <- stats::aggregate(outcome ~ treatment + group,
                                 function(x) sd(x)/sqrt(length(x)), data = data)
      mwide  <- stats::reshape(data = magg, timevar = "treatment",
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
    if(effect %in% c("logOR", "logRR")) {
      v <- rare_event_correction
      binary_data_table <-
        do.call(rbind, by(data, list(data$group), function(x) {
          with(x,
               data.frame(
                 group = unique(as.character(group)),
                 a     = sum(outcome[treatment == 1]),
                 n1    = sum(treatment == 1),
                 c     = sum(outcome[treatment == 0]),
                 n2    = sum(treatment == 0),
                 b     = sum(treatment == 1) - sum(outcome[treatment == 1]),
                 d     = sum(treatment == 0) - sum(outcome[treatment == 0])))
        }))

      rare <- with(binary_data_table, (a == 0 | b == 0 | c == 0 | d == 0))

      if(sum(rare) == 1 && rare_event_correction == 0.25)
        message("Applied default rare event correction (0.25) in 1 study")
      if(sum(rare) > 1 && rare_event_correction == 0.25)
        message("Applied default rare event corrections (0.25) in ", sum(rare), " studies")

      if(correction_type == "single")
        cc_value <- v*rare
      if(correction_type == "all")
        cc_value <- rare_event_correction

      binary_data_table$a  <- cc_value + binary_data_table$a
      binary_data_table$b  <- cc_value + binary_data_table$b
      binary_data_table$c  <- cc_value + binary_data_table$c
      binary_data_table$d  <- cc_value + binary_data_table$d
      binary_data_table$n1 <- cc_value + binary_data_table$n1
      binary_data_table$n2 <- cc_value + binary_data_table$n2


      out <- binary_data_table
      rownames(out) <- NULL

      if(effect == "logRR") {
        out$tau <- with(out, log((a/(a+b))/(c/(c+d))))
        out$se  <- with(out, sqrt(1/a + 1/c - 1/(a+b) - 1/(c+d)))
      }
      if(effect == "logOR") {
        out$tau <- with(out, log((a*d)/(b*c)))
        out$se  <- with(out, sqrt(1/a + 1/b + 1/c + 1/d))
      }

    }
  } else {
    out <- data
  }

  out
}

