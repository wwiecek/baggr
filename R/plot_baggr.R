#' Plotting method in baggr package
#'
#' @param bg object of class \code{baggr}
#' @param style one of 'areas', 'intervals' or (WIP)
#' @param ... extra arguments passed to bayesplot functions
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' @import bayesplot

plot.baggr <- function(bg, style = "areas", ...) {
  m <- study_effects(bg)
  if(!(style %in% c("areas", "intervals")))
    stop('plot "style" argument must be one of: "areas", "intervals"')
  p <- switch(style,
              "areas" = bayesplot::mcmc_areas(m, ...),
              "intervals" = bayesplot::mcmc_intervals(m, ...))
  p + ggplot2::labs(x = "Effect size")
}
