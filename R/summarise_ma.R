#' @title Convert from individual to summary data in meta-analyses
#' @description Allows only one-way conversion from full to summary data.
#'              Input must be pre-formatted appropriately.
#'
#' @param data data.frame of individual-level observations
#'             with columns \code{outcome} (numeric),
#'             \code{treatment} (values 0 and 1) and \code{grouping}
#' @return data.frame with columns \code{mu}, \code{se.mu}, \code{tau} and \code{se.tau}
#' @details
#' The conversions will typically happen automatically when data is fed to baggr()
#' function. This function can be used to explicitly convert from full to reduced
#' data without analysing it in any model.
#' @author Witold Wiecek, Rachael Meager
#' @seealso \code{\link{convert_inputs}}
#' @export

summarise_ma <- function(data) {
  magg  <- aggregate(outcome ~ treatment + grouping,
                     mean, data = data)
  seagg <- aggregate(outcome ~ treatment + grouping,
                     function(x) sd(x)/sqrt(length(x)), data = data)
  mwide <- stats::reshape(data = magg, timevar = "treatment",
                          idvar = "grouping", direction = "wide")
  sewide <- stats::reshape(data = seagg, timevar = "treatment",
                           idvar = "grouping", direction = "wide")
  data.frame(grouping = mwide$grouping,
             mu = mwide$outcome.0,
             tau = mwide$outcome.1 - mwide$outcome.0,
             se.mu = sewide$outcome.0,
             se.tau = sqrt(sewide$outcome.0^2 + sewide$outcome.1^2))
}

