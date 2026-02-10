#' Relative publication probabilities of a baggr model
#'
#' Extracts the posterior draws (or summaries) of the selection-model parameters
#' that describe the relative publication probabilities across |z|-intervals.
#'
#' @param bg a [baggr::baggr()] model fitted with a selection-enabled likelihood
#'   (currently `model = "rubin"`).
#' @param summary logical; if `TRUE`, returns summary statistics (mean and
#'   the uncertainty interval defined by `interval`).
#' @param interval uncertainty interval width (numeric between 0 and 1),
#'   used only when `summary = TRUE`.
#' @return If `summary = FALSE`, a matrix of posterior draws where each column
#'   corresponds to a |z|-interval. Otherwise, a matrix with one row per interval
#'   and columns containing the mean and bounds of the uncertainty interval.
#' @export
#' @importFrom rstan extract
selection <- function(bg, summary = TRUE, interval = 0.95) {
  check_if_baggr(bg)

  if(!bg$model %in% c("rubin"))
    stop("Selection parameters are only available for selection models.")

  if(is.null(bg$inputs$M) || bg$inputs$M <= 0)
    stop("The fitted model does not define any selection intervals.")

  omega_draws <- rstan::extract(bg$fit, pars = "omega")[[1]]
  if(is.null(dim(omega_draws)))
    omega_draws <- matrix(omega_draws, ncol = bg$inputs$M)
  if(!is.matrix(omega_draws))
    omega_draws <- as.matrix(omega_draws)

  if(summary)
    return(mint(omega_draws, int = interval, median = TRUE, sd = TRUE))

  omega_draws
}
