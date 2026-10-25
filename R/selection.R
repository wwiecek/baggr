#' Relative publication probabilities of a baggr model
#'
#' Extracts the posterior draws (or summaries) of the selection-model parameters
#' that describe the relative publication probabilities across |z|-intervals.
#' In normal-normal hierarchical (Rubin) models, `z = yi / sei`; for example,
#' `selection = 1.96` estimates the relative publication probability below the
#' conventional two-sided 5% normal threshold compared with the implicit
#' highest-|z| interval.
#'
#' `selection()` returns draws or summaries of the selection weights only. These
#' are relative publication probabilities, not treatment effects or corrected
#' treatment effects. The highest-|z| interval is the reference category and is
#' fixed to 1, so it is not returned.
#'
#' With `summary = TRUE`, summaries are numeric matrices with rows for
#' the quantities being summarised and columns `mean`, `median`, `sd`,
#' `lower`, `upper`, `width`, in that order. The bounds form a central
#' interval; `width` is its requested probability mass, even when estimates
#' are unavailable. For multiple effects, the result is a rows by statistics
#' by effects array, with effect names in the third dimension. A single-effect
#' matrix stores its effect name in the `effect` attribute. Row names are
#' preserved. Aggregate summaries have one row; predictions conditional on
#' `newdata` have one row per observation.
#'
#' @param bg a [baggr::baggr()] model fitted with a selection-enabled likelihood
#'   (currently `model = "rubin"`).
#' @param summary logical; if `TRUE`, returns summary statistics (mean and
#'   the uncertainty interval defined by `prob`).
#' @param prob Probability mass of the central uncertainty interval, on the
#'   0--1 scale; used only when `summary = TRUE`.
#' @param interval Compatibility alias for `prob`; supply only one argument.
#' @return If `summary = FALSE`, a matrix of posterior draws where each column
#'   corresponds to one non-reference |z|-interval. Otherwise, a matrix with one
#'   row per interval in the summary layout described here, labelled with
#'   `effect = "selection"`. The final highest-|z| interval is fixed to relative
#'   probability 1 and is therefore not returned.
#' @export
#' @importFrom rstan extract
selection <- function(bg, summary = TRUE, prob = 0.95, interval = NULL) {
  prob <- interval_prob(prob, interval, missing(prob), missing(interval))
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
    return(summary_draws(omega_draws, prob, "selection"))

  omega_draws
}
