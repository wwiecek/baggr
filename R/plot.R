#' Generic plot for baggr package
#'
#' Using generic `plot()` on `baggr` output invokes \code{\link{baggr_plot}} visual.
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

#' Set, get, and replace themes for baggr plots
#' @rdname baggr_theme_set
#' @importFrom ggplot2 theme_grey
#' @importFrom bayesplot theme_default
#' @param new New theme to use for all baggr plots
#' @description These functions get, set, and modify the ggplot2 themes
#' of the baggr plots. By default, these match the bayesplot's package
#' theme to make it easier to form cohesive graphs across this package
#' and others.
#' @return The get method returns the current theme, but all of the
#' others invisibly return the old theme.
#' @export
baggr_theme_set <- function(new = bayesplot::theme_default()) {
  bayesplot::bayesplot_theme_set(new)
}

#' @rdname baggr_theme_set
#' @importFrom bayesplot bayesplot_theme_get
baggr_theme_get <- function() {
  bayesplot::bayesplot_theme_get()
}

#' @rdname baggr_theme_set
#' @importFrom bayesplot bayesplot_theme_update
#' @export
baggr_theme_update <- function(...) {
  bayesplot::bayesplot_theme_update(...)
}

#' @rdname baggr_theme_set
#' @export
#' @importFrom ggplot2 %+replace%
#' @importFrom bayesplot bayesplot_theme_replace
baggr_theme_replace <- function(...) {
  bayesplot::bayesplot_theme_replace(...)
}
