#' Plotting method in baggr package
#'
#' Extracts study effects from the  \code{baggr} model and sends them to
#' one of \code{bayesplot} package plotting functions.
#'
#' @param bg object of class \code{baggr}
#' @param hyper logical; show hypereffect as the last row of the plot?
#' @param style one of \code{areas}, \code{intervals}
#' @param prob Probability mass for the inner interval in visualisation
#' @param prob_outer Probability mass for the outer interval in visualisation
#' @param vline logical; show vertical line through 0 in the plot?
#' @param order logical; sort groups by magnitude of treatment effect?
#' @param ... extra arguments to pass to the `bayesplot` functions
#'
#' @return ggplot2 object
#'
#' @examples
#' fit <- baggr(schools, pooling = "none")
#' plot(fit)
#' plot(fit, style = "areas", order = FALSE)
#'
#' @export
#' @import ggplot2
#' @import bayesplot
#'
#' @author Witold Wiecek, Rachael Meager
#' @seealso [bayesplot::MCMC-intervals] for more information about _bayesplot_ functionality;
#'          [forest_plot] for a typical meta-analysis alternative; [effect_plot] for plotting
#'          treatment effects for a new group

baggr_plot <- function(bg, hyper=FALSE,
                       style = "intervals",
                       prob = 0.5, prob_outer = 0.95,
                       vline = TRUE, order = TRUE, ...) {
  if(attr(bg, "ppd")){
    message("Baggr model is prior predictive; returning effect_plot().")
    return(effect_plot(bg))
  }
  m <- group_effects(bg)
  if(hyper)
    te <- treatment_effect(bg)$tau
  effect_labels <- bg$effects

  if(!(style %in% c("areas", "intervals")))
    stop('plot "style" argument must be one of: "areas", "intervals"')

  ret_list <- lapply(as.list(1:dim(m)[3]), function(i) {
    if(order)
      mat_to_plot <- m[,order(apply(m[,,i], 2, mean)),i] #assigning to m[,,i] wouldn't reorder dimnames
    else
      mat_to_plot <- m[,,i]
    if(hyper){
      mat_to_plot <- cbind(mat_to_plot, te)
      colnames(mat_to_plot)[ncol(mat_to_plot)] <- effect_labels[i]
    }
    p <- switch(style,
                "areas"     = bayesplot::mcmc_areas(mat_to_plot, prob = prob, prob_outer = prob_outer, ...),
                "intervals" = bayesplot::mcmc_intervals(mat_to_plot, prob = prob, prob_outer = prob_outer, ...))
    p +
      ggplot2::labs(x = paste("Effect size:", bg$effects[i])) +
      baggr_theme_get() +
      {if(vline) geom_vline(xintercept = 0, lty = "dashed")} +
      {if(hyper & style == "intervals") geom_hline(yintercept = 1.5)}
  })

  if(length(ret_list) == 1)
    return(ret_list[[1]])
  else
    return(ret_list)
}
