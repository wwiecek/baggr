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
#' @param ... A named list of theme settings
#' @description These functions get, set, and modify the ggplot2 themes
#' of the baggr plots. \code{baggr_theme_get()} returns a ggplot2 theme function for
#' adding themes to a plot. \code{baggr_theme_set()} assigns a new theme
#' for all plots of baggr objects. \code{baggr_theme_update()} edits a specific
#' theme element for the current theme while holding the theme's
#' other aspects constant. \code{baggr_theme_replace()} is used for
#' wholesale replacing aspects of a plot's theme (see [ggplot2::theme_replace()]).
#'
#' @details
#'
#' Under the hood, many of the visualizations rely on the
#' bayesplot package, and thus these leverage the [bayesplot::bayesplot_theme_get()]
#' functions. By default, these match the bayesplot's package
#' theme to make it easier to form cohesive graphs across this package
#' and others. The trickiest of these to use is \code{baggr_theme_replace};
#' 9 times out of 10 you want baggr_theme_update.
#'
#' @examples
#'
#' # make plot look like default ggplots
#'
#' library(ggplot2)
#'
#' fit <- baggr(schools)
#' baggr_theme_set(theme_grey())
#' baggr_plot(fit)
#'
#' # use baggr_theme_get to return theme elements for current theme
#' qplot(mtcars$mpg) + baggr_theme_get()
#'
#' # update specific aspect of theme you are interested in
#' library(extrafont)
#' baggr_theme_update(text = element_text(family = "Wingdings"))
#'
#' # undo that silliness
#' baggr_theme_update(text = element_text(family = "Calibri"))
#'
#' # update and replace are similar, but replace overwrites the
#' # whole element, update just edits the aspect of the element
#' # that you give it
#' # this will error:
#' # baggr_theme_replace(text = element_text(family = "Times"))
#' # baggr_plot(fit)
#' # because it deleted everything else to do with text elements
#'
#' @return The get method returns the current theme, but all of the
#' others invisibly return the old theme.
#' @seealso [bayesplot::bayesplot_theme_get]
#' @export
baggr_theme_set <- function(new = bayesplot::theme_default()) {
  bayesplot::bayesplot_theme_set(new)
}

#' @rdname baggr_theme_set
#' @importFrom bayesplot bayesplot_theme_get
#' @export
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
