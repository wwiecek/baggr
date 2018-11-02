#' Compare a(ny) number of baggr models side by side
#'
#' @param ... Either any number of objects of class `baggr`
#'            (you can name your objects, see example below)
#'            or the same arguments you'd pass to baggr()
#'            function, but with `pooling = ...` omitted.
#'            In the latter case 3 models will be run, with
#'            pooling set to `none`, `partial` and `full`.
#' @param compare can be 'groups' (compare different groups in each model),
#'                or 'ate' (compare aggregate treatment effects only)
#' @param style What kind of plot to display - see options for `plot.baggr`
#' @param arrange If `single`, generate a single plot, if `grid`
#'                display multiple plots side-by-side.
#' @return ggplot is plotted and an extra comparison printed.
#'         Returned value is a list that contains
#'         the models, the plot and printed summaries. (WIP)
#' @author Witold Wiecek
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
#' @author Witold Wiecek, Rachael Meager

baggr_compare <- function(...,
                          compare = "group",
                          style = "areas",
                          arrange = "single") {
  l <- list(...)
  if(length(l) == 0)
    stop("Provide baggr models or model specification.")
  if(all(unlist(lapply(l, class)) == "baggr")) {
    models <- l
  } else {
    if("pooling" %in% names(l))
      stop("Can't run the model comparison with pooling setting
            already set to a particular value.")
    models <- lapply(list("none", "partial"), function(pool){
      try(do.call(baggr, c(l, "pooling" = pool)))
    })
    names(models) <- c("none", "partial")
  }
  if(arrange == "grid") {
    plots <- lapply(models, plot, style = style)
    gridExtra::grid.arrange(grobs = plots, ncol = length(plots))
  }
  if(arrange == "single") {
    # Note: pipe operators are not used here for compatibility and
    # to reduce dependencies
    df <- dplyr::bind_rows(
      lapply(models, function(x) as.data.frame(study_effects(x, summary = TRUE))),
      .id = "model")
    df <- dplyr::mutate(df, parameter = factor(parameter, levels = unique(parameter[order(median)])))
    ggplot2::ggplot(df, aes(x = parameter, y = median, ymin = lci, ymax = uci,
                   group = interaction(model),
                   color = model)) +
      geom_point(size = 2, position=position_dodge(width=0.5)) +
      # geom_jitter(size = 2) +
      geom_errorbar(size = 1.2, width = 0, alpha = .5, position=position_dodge(width=0.5)) +
      coord_flip() +
      labs(x = "", y = "Treatment effect (95% interval)",
           title = "Comparison of treatment effects by group")
  }
}

