#' Plotting method in baggr package
#'
#' Extracts study effects from the  \code{baggr} model and plots them,
#' possibly next to the hypereffect estimate.
#'
#' @param bg object of class \code{baggr}
#' @param hyper logical; show hypereffect as the last row of the plot?
#'              alternatively you can pass colour for the hypermean row,
#'              e.g. `hyper = "red"`
#' @param style `"forest_plot"` imitates the visual style of forest plots
#'              and also prints means and intervals next to each row;
#'              `"intervals"` (default) or `"areas"` use package `bayesplot` styles
#' @param transform a function (e.g. `exp()`, `log()`) to apply to the
#'                  values of group (and hyper, if `hyper=TRUE`) effects
#'                  before plotting; when working with effects that are on
#'                  log scale, exponent transform is used automatically, you can
#'                  plot on log scale by setting `transform = identity`
#' @param prob Probability mass for the inner interval in visualisation
#' @param prob_outer Probability mass for the outer interval in visualisation
#' @param vline logical; show vertical line through 0 in the plot?
#' @param order logical; sort groups by magnitude of treatment effect?
#' @param values_outer logical; use the interval corresponding to `prob_outer` when `style = "forest_plot"`?
#'                     if not, the "inner" interval (`prob`) is used
#' @param values_size size of the text values in the plot when `style = "forest_plot"`
#' @param values_digits number of significant digits to use when `style = "forest_plot"`
#' @param ... extra arguments to pass to the `bayesplot` functions
#'
#' @return ggplot2 object
#'
#' @examples
#' fit <- baggr(schools, pooling = "none")
#' plot(fit, hyper = "red")
#' plot(fit, style = "areas", order = FALSE)
#' plot(fit, style = "forest_plot", order = FALSE)
#'
#' @export
#' @import ggplot2
#' @import bayesplot
#'
#' @author Witold Wiecek; the visual style is based on _bayesplot_ package
#' @seealso [bayesplot::MCMC-intervals] for more information about _bayesplot_ functionality;
#'          [forest_plot] for a typical meta-analysis alternative (which you can imitate using `style = "forest_plot"`);
#'          [effect_plot] for plotting treatment effects for a new group

baggr_plot <- function(bg, hyper=FALSE,
                       style = c("intervals", "areas", "forest_plot"),
                       transform = NULL,
                       prob = 0.5, prob_outer = 0.95,
                       vline = TRUE, order = TRUE,
                       values_outer = TRUE,
                       values_size = 4,
                       values_digits = 1,
                       ...) {
  if(attr(bg, "ppd")){
    message("Baggr model is prior predictive; returning effect_plot().")
    return(effect_plot(bg))
  }
  if(is.character(hyper)){
    hyper_colour <- hyper
    hyper <- TRUE
  } else
    hyper_colour <- NULL

  m <- group_effects(bg, transform = transform)
  effect_labels <- bg$effects

  style <- match.arg(style, c("intervals", "areas", "forest_plot"))

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
      ate <- treatment_effect(bg, transform = transform)$tau
      if(length(bg$effects) > 1) #ATE is a matrix
        mat_to_plot <- cbind(mat_to_plot, ate[,i])
      else #ATE is a vector
        mat_to_plot <- cbind(mat_to_plot, ate)
      colnames(mat_to_plot)[ncol(mat_to_plot)] <- "Hypermean"
    }
    if(!style=="forest_plot"){
      p <- switch(style,
                  "areas"     = bayesplot::mcmc_areas(mat_to_plot, prob = prob,
                                                      prob_outer = prob_outer, ...),
                  "intervals" = bayesplot::mcmc_intervals(mat_to_plot, prob = prob,
                                                          prob_outer = prob_outer, ...))
      p <- p +
        ggplot2::labs(x = paste("Effect on", bg$effects[i])) +
        baggr_theme_get() +
        {if(hyper & style == "intervals") geom_hline(yintercept = 1.5)} +
        {if(vline) geom_vline(xintercept = vline_value, lty = "dashed")}

    }
    if(style=="forest_plot"){
      if(values_digits<3)
        buffer = (values_size*8)/3.5
      else
        buffer = (values_size*values_digits*3)/4

      parameter <- ll <- l <- m <- h <- hh <- NULL

      data <- bayesplot::mcmc_intervals_data(mat_to_plot, prob = prob,
                                             prob_outer = prob_outer, ...)

      if(values_outer){
        lci_col <- "ll"; uci_col <- "hh"
        value_text <- paste3_formatter(data$m, data$ll, data$hh, values_digits)
      }else{
        lci_col <- "l"; uci_col <- "h"
        value_text <- paste3_formatter(data$m, data$l, data$h, values_digits)
      }

      p <- ggplot2::ggplot(data, aes(x = m, y = parameter)) +
              ggplot2::scale_y_discrete(limits = rev) +
              ggplot2::geom_point(size=5,shape=15) +
              ggplot2::geom_errorbarh(aes(y = parameter,xmin = ll, xmax = hh),
                                      linewidth=0.5, height=0.15,inherit.aes=FALSE) +
              ggplot2::geom_errorbarh(aes(y = parameter,xmin = l, xmax = h),
                                      linewidth=2,height=0,inherit.aes=FALSE) +
              ggplot2::labs(x = paste("Effect on", bg$effects[i])) +
              ggplot2::theme(panel.background = element_rect(fill = "white", colour = "white"),
                legend.position="none",
                axis.title.y=element_blank(),
                axis.title.x=element_text(size=values_size*3+2),
                axis.text.y=element_text(size=values_size*3),
                axis.text.x=element_text(size=values_size*3),
                axis.ticks.y = element_blank(),
                plot.margin = unit(c(1,buffer,1,1), "lines")) +
              ggplot2::geom_text(aes(x=max(hh),hjust=-0.25, label=value_text), size=values_size) +
              ggplot2::coord_cartesian(clip = 'off') +
              {if(hyper) geom_hline(yintercept = 1.5)} +
              {if(vline) geom_vline(xintercept = vline_value, lty = "dashed")}
    }

    if(!is.null(hyper_colour) && style != "areas") {
      p <- add_color_to_plot(p, c("Hypermean" = hyper_colour))
    }

    p
  })



  if(length(ret_list) == 1)
    return(ret_list[[1]])
  else
    return(ret_list)
}


#' Add colors to baggr plots
#'
#' @param p    A ggplot object to add colors to
#' @param what A named vector, e.g. `c(Hypermean = "red", "Group A" = "green")`.
#' @importFrom ggplotify as.ggplot

add_color_to_plot <- function(p,what){
  target <- names(what)
  colors <- what
  where <- sapply(target, function(x) which(p$data == x))
  pb <- ggplot2::ggplot_build(p)
  for(i in 1:length(what)){
    pb$data[[2]][where[[i]],'colour'] <- colors[[i]]
    pb$data[[3]][where[[i]],'colour'] <- colors[[i]]
    pb$data[[4]][where[[i]],'colour'] <- colors[[i]]
    pb$data <- pb$data
  }
  pb <- ggplot2::ggplot_gtable(pb)
  return(ggplotify::as.ggplot(pb))
}
