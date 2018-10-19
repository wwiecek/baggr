#' Plotting method in baggr package
#'
#' Extracts study effects from the  \code{baggr} model and sends them to
#' one of \code{bayesplot} package plotting functions (with minimal embelishments).
#'
#' @param bg object of class \code{baggr}
#' @param style one of \code{areas}, \code{intervals}
#' @param prob Probability mass for the inner interval in visualisation
#' @param prob_outer Probability mass for the outer interval in visualisation
#' @param vline logical; show vertical line through 0 in the plot
#' @param order logical; sort groups by magnitude of treatment effect?
#' @param ... extra arguments passed to bayesplot functions
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' @import bayesplot
#' @author Witold Wiecek, Rachael Meager
#'

plot.baggr <- function(bg, style = "intervals", prob = 0.5, prob_outer = 0.95,
                       vline = TRUE, order = TRUE, ...) {
  m <- study_effects(bg)
  if(!(style %in% c("areas", "intervals")))
    stop('plot "style" argument must be one of: "areas", "intervals"')

  ret_list <- lapply(as.list(1:dim(m)[3]), function(i) {
    if(order)
      m[,,i] <- m[,order(apply(m[,,i], 2, mean)),i]
    p <- switch(style,
                "areas"     = bayesplot::mcmc_areas(m[,,i],     prob = prob, prob_outer = prob_outer, ...),
                "intervals" = bayesplot::mcmc_intervals(m[,,i], prob = prob, prob_outer = prob_outer, ...))
    p + ggplot2::labs(x = paste("Effect size:", bg$effects[i])) +
    {if(vline) geom_vline(xintercept = 0, lty = "dashed")}
  })

  if(length(ret_list) == 1)
    return(ret_list[[1]])
  else
    return(ret_list)
}
