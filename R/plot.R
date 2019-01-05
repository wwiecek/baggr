#' Generic plot for baggr package
#'
#' Using generic `plot()` on `baggr` output invokes \code{\link{baggr::baggr_plot}} visual.
#' See therein for customisation options. Note that plot output is `ggplot2` object.`
#'
#' @param x object of class `baggr`
#' @param ... optional arguments, see `baggr_plot`
#' @return ggplot2 object from `baggr_plot`
#' @export
#' @author Witold Wiecek

plot.baggr <- function(x, ...) {
  l <- list(...)
  l[["bg"]] <- x
  do.call(baggr_plot, l)
}
