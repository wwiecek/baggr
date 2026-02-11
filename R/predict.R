#' Predict treatment effects from baggr model
#'
#' Alias for [effect_draw], returning draws from the predictive distribution
#' for a `baggr` model.
#'
#' @param object A `baggr` class object.
#' @param ... Additional arguments passed to [effect_draw].
#' @return See [effect_draw].
#' @export
predict.baggr <- function(object, ...) {
  effect_draw(object, ...)
}
