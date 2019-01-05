#' @title Convert from individual to summary data in meta-analyses
#' @description Allows only one-way conversion from full to summary data.
#'              Input must be pre-formatted appropriately.
#'
#' @param data data.frame of individual-level observations
#'             with columns \code{outcome} (numeric),
#'             \code{treatment} (values 0 and 1) and \code{group} (numeric, character or factor)
#' @param standardise logical; if TRUE, values of outcome are standardised within each group
#' @return data.frame with columns \code{mu}, \code{se.mu}, \code{tau} and \code{se.tau}
#' @details
#' The conversions will typically happen automatically when data is fed to baggr()
#' function. This function can be used to explicitly convert from full to reduced
#' data without analysing it in any model.
#' @author Witold Wiecek, Rachael Meager
#' @seealso \code{\link{convert_inputs}}
#' @export
#' @import stats

summarise_ma <- function(data, standardise = FALSE) {
  if(any(!stats::complete.cases(data[,c("treatment", "group", "outcome")])))
    warning("NA values present in data - they were dropped when summarising")

  if(standardise) {
    agg <- stats::aggregate(outcome ~ group, function(x) {c(mean=mean(x), sd=sd(x))}, data = data)
    means <- agg$outcome[,"mean"]
    sds <- agg$outcome[,"sd"]
    names(means) <- names(sds) <- agg$group

    data$outcome <- (data$outcome - means[data$group]) / sds[data$group]
  }

  magg  <- stats::aggregate(outcome ~ treatment + group,
                     mean, data = data)
  seagg <- stats::aggregate(outcome ~ treatment + group,
                     function(x) sd(x)/sqrt(length(x)), data = data)
  mwide <- stats::reshape(data = magg, timevar = "treatment",
                          idvar = "group", direction = "wide")
  sewide <- stats::reshape(data = seagg, timevar = "treatment",
                           idvar = "group", direction = "wide")
  data.frame(group = mwide$group,
             mu = mwide$outcome.0,
             tau = mwide$outcome.1 - mwide$outcome.0,
             se.mu = sewide$outcome.0,
             se.tau = sqrt(sewide$outcome.0^2 + sewide$outcome.1^2))
}

