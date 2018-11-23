#' Compare a(ny) number of baggr models side by side
#'
#' @param ... Either any number of objects of class `baggr`
#'            (you can name your objects, see example below)
#'            or the same arguments you'd pass to baggr()
#'            function, but with `pooling = ...` omitted.
#'            In the latter case 3 models will be run, with
#'            pooling set to `none`, `partial` and `full`.
#' @param style What kind of plot to display - see options for `plot.baggr`
#' @param arrange If `single`, generate a single plot, if `grid`
#'                display multiple plots side-by-side.
#' @return ggplot is plotted and an extra comparison printed.
#'         Returned value is a list that contains
#'         the models, the plot and printed summaries. (WIP)
#' @author Witold Wiecek
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr bind_rows
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
#' @author Witold Wiecek, Rachael Meager

baggr_compare <- function(...,
                          compare = "group",
                          style   = "areas",
                          arrange = "single") {
  l <- list(...)
  if(length(l) == 0)
    stop("Must provide baggr models or model specification.")
  if(all(unlist(lapply(l, class)) == "baggr")) {
    models <- l
  } else {
    if("pooling" %in% names(l))
      stop("Can't run the model comparison with pooling setting
            already set to a particular value.")

    # This is where the models are run:
    models <- lapply(list("none", "partial", "full"), function(pool){
      try(do.call(baggr, c(l, "pooling" = pool)))
    })
    names(models) <- c("none", "partial", "full")
  }
  if(arrange == "grid") {
    plots <- lapply(models, plot, style = style)
    gridExtra::grid.arrange(grobs = plots, ncol = length(plots))
  }
  if(arrange == "single") {
    # Note: pipe operators are not used here for compatibility and
    # to reduce dependencies but i still went for dplyr as it's so much simpler
    df_groups <- dplyr::bind_rows(
      lapply(models, function(x) {
        # will need to be modified for quantiles models case:
        m <- as.data.frame(study_effects(x, summary = TRUE)[,,1])
        m$parameter <- rownames(m)
        m
      }), .id = "model")
    df_groups$parameter <- factor(df_groups$parameter,
                                  levels = unique(df_groups$parameter[order(df_groups$median)]))

    # now treatment effect:
    # df_trt <- dplyr::bind_rows(
    #   lapply(models, function(x) {
    #     m <- mint(treatment_effect(x)$tau)
    #     m$parameter <- rownames(m)
    #     m
    #   }), .id = "model")

    # df <- rbind(df_groups, df_trt)
    df <- df_groups
    comparison_plot <- ggplot2::ggplot(df, aes(x = parameter, y = median, ymin = lci, ymax = uci,
                            group = interaction(model),
                            color = model)) +
      geom_point(size = 2, position=position_dodge(width=0.5)) +
      # geom_jitter(size = 2) +
      geom_errorbar(size = 1.2, width = 0, alpha = .5, position=position_dodge(width=0.5)) +
      coord_flip() +
      labs(x = "", y = "Treatment effect (95% interval)",
           title = "Comparison of treatment effects by group")
  }

  # render the plot (for now not saved)
  plot(comparison_plot)

  return(models)
}
