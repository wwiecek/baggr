#' Plotting method in baggr package
#'
#' Extracts study effects from the  \code{baggr} model and sends them to
#' one of \code{bayesplot} package plotting functions.
#'
#' @param bg object of class \code{baggr}
#' @param hyper logical; show hypereffect as the last row of the plot?
#' @param style either `"intervals"` or `"areas"`
#' @param transform a function (e.g. `exp()`, `log()`) to apply to the
#'                  values of group (and hyper, if `hyper=TRUE`) effects
#'                  before plotting; when working with effects that are on
#'                  log scale, exponent transform is used automatically, you can
#'                  plot on log scale by setting `transform = identity`
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
#' @author Witold Wiecek; the visual style is based on _bayesplot_ package
#' @seealso [bayesplot::MCMC-intervals] for more information about _bayesplot_ functionality;
#'          [forest_plot] for a typical meta-analysis alternative; [effect_plot] for plotting
#'          treatment effects for a new group

baggr_plot <- function(bg, hyper=FALSE,
                       style = "intervals",
                       transform = NULL,
                       prob = 0.5, prob_outer = 0.95,
                       vline = FALSE, order = TRUE, ...) {
  if(attr(bg, "ppd")){
    message("Baggr model is prior predictive; returning effect_plot().")
    return(effect_plot(bg))
  }
  m <- group_effects(bg, transform = transform)
  effect_labels <- bg$effects

  if(!(style %in% c("areas", "intervals")))
    stop('plot "style" argument must be one of: "areas", "intervals"')

  if(!is.null(transform))
    vline_value <- do.call(transform, list(0))
  else
    vline_value <- 0

  if(bg$n_groups == 1)
    stop("You can only plot meta-analyses with more than 1 group.")

  ret_list <- lapply(as.list(1:dim(m)[3]), function(i) {
    if(order)
      #assigning to m[,,i] wouldn't reorder dimnames
      mat_to_plot <- m[,order(apply(m[,,i], 2, mean)),i]
    else
      mat_to_plot <- m[,,i]
    if(hyper && bg$pooling != "none"){
      ate <- hypermean(bg, transform = transform,message=FALSE) 
      if(length(bg$effects) > 1) #ATE is a matrix
        mat_to_plot <- cbind(mat_to_plot, ate[,i])
      else #ATE is a vector
        mat_to_plot <- cbind(mat_to_plot, ate)
      colnames(mat_to_plot)[ncol(mat_to_plot)] <- "Hypermean"
    }

    p <- switch(style,
                "areas"     = bayesplot::mcmc_areas(mat_to_plot, prob = prob,
                                                    prob_outer = prob_outer, ...),
                "intervals" = bayesplot::mcmc_intervals(mat_to_plot, prob = prob,
                                                        prob_outer = prob_outer, ...))
    p +
      ggplot2::labs(x = paste("Effect on", bg$effects[i])) +
      baggr_theme_get() +
      {if(hyper & style == "intervals") geom_hline(yintercept = 1.5)} +
      {if(vline) geom_vline(xintercept = vline_value, lty = "dashed")}
  })

  if(length(ret_list) == 1)
    return(ret_list[[1]])
  else
    return(ret_list)
}
