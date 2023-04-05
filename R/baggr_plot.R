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
#' @param add_color_to  A list of objects and colors to be added to the plots.
#'                      Lists are constructed such that objects and colors are given
#'                      in pairs. For example, a list to turn School A red and School B
#'                      blue would be constructed as `list('School A','red','School B','blue')`.
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
#' @import ggplotify
#'
#' @author Witold Wiecek; the visual style is based on _bayesplot_ package
#' @seealso [bayesplot::MCMC-intervals] for more information about _bayesplot_ functionality;
#'          [forest_plot] for a typical meta-analysis alternative; [effect_plot] for plotting
#'          treatment effects for a new group

baggr_plot <- function(bg, hyper=FALSE,
                       style = "intervals",
                       transform = NULL,
                       prob = 0.5, prob_outer = 0.95,
                       vline = FALSE, order = TRUE, 
                       add_color_to=NULL,...) {
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
      ate <- treatment_effect(bg, transform = transform)$tau
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
      {if(vline) geom_vline(xintercept = vline_value, lty = "dashed")} 
    
    if(length(add_color_to)>=1){ 
      p <- p + {if(hyper & style == "intervals") geom_hline(yintercept = 1.5)}
      p <- add_color(p,add_color_to,style) 
      
    } else {
      p <- p + {if(hyper & style == "intervals") geom_hline(yintercept = 1.5)}
      
    }
      

  })

  if(length(ret_list) == 1)
    return(ret_list[[1]])
  else
    return(ret_list)
}


#' Add colors to baggr plots
#'
#' @param p             A ggplot object to add colors to
#' @param add_color_to  A list of objects and colors to be added to the plots.
#'                      Lists are constructed such that objects and colors are given
#'                      in pairs. For example, a list to turn School A red and School B
#'                      blue would be constructed as `list('School A','red','School B','blue')`.
#'                      CSS hex numbers are also supported. 
#' @param style         The style of the plot. Currently, only plots of style `intervals` are 
#'                      supported.
#' @examples
#' bg <- baggr(schools)
#' colors <- list('School A','red','School B','blue')
#' plot(bg,add_color_to=colors)

add_color <- function(p,add_color_to,style){
  if(style == 'areas'){
    stop('Adding color is currently only supported for plots in interval style. Please change the `style` argument to `intervals` and try again.')
  }
  if(style == "intervals"){
    obj_names <- add_color_to[c(TRUE,FALSE)]
    objects <- NULL
    for(name in obj_names){
      objects[length(objects)+1] <- which(p$data==name)
    }
    colors <- add_color_to[c(FALSE,TRUE)]
    pb <- ggplot2::ggplot_build(p)
    for(obj in 1:length(objects)){
      pb$data[[2]][objects[obj],'colour'] <- colors[obj]
      pb$data[[3]][objects[obj],'colour'] <- colors[obj] 
      pb$data[[4]][objects[obj],'colour'] <- colors[obj] 
      pb$data <- pb$data
    }
    pb <- ggplot2::ggplot_gtable(pb)
    return(ggplotify::as.ggplot(pb))
  }
}

