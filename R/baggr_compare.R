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
#' @import ggplot2
#' @author Witold Wiecek, Rachael Meager

baggr_compare <- function(...,
                          style   = "areas",
                          arrange = "single") {
  l <- list(...)
  if(length(l) == 0)
    stop("Must provide baggr models or model specification.")
  if(all(unlist(lapply(l, class)) == "baggr")) {
    return_models_flag <- 0
    if(is.null(names(l)))
      names(l) <- paste("Model", 1:length(l))
    models <- l
  } else {
    if("pooling" %in% names(l))
      stop("Can't run the model comparison with pooling setting
            already set to a particular value.")
    return_models_flag <- 1
    # This is where the models are run:
    models <- lapply(list("none", "partial", "full"), function(pool){
      try(do.call(baggr, c(l, "pooling" = pool)))
    })
    names(models) <- c("none", "partial", "full")
  }

  effect_names <- lapply(models, function(x) x$effects)
  # quite a mouthful:
  if(!all(unlist(lapply(effect_names, function(x) all.equal(effect_names[[1]], x))) == 1))
    stop("Models must have the same effects to be comparable")
  effect_names <- effect_names[[1]]

  if(arrange == "grid") {
    plots <- lapply(models, baggr_plot, style = style, order = FALSE)
    grid_width <- length(plots)
    # if each plots element contains multiple plots (like with quantiles):
    if(class(plots[[1]])[1] == "list")
      plots <- unlist(plots, recursive = FALSE)
    gridExtra::grid.arrange(grobs = plots, ncol = grid_width)
  }
  if(arrange == "single") {
    plots <- lapply(as.list(1:(length(effect_names))), function(i) {
      # Note: pipe operators are dplyr not used here for compatibility
      ll <- lapply(models, function(x) {
        # will need to be modified for quantiles models case:
        m <- as.data.frame(study_effects(x, summary = TRUE)[,,i])
        m$group <- rownames(m)
        m
      })
      df_groups <- data.frame()
      for(j in 1:length(ll))
        df_groups <- rbind(df_groups,
                           data.frame(model = names(ll)[j], ll[[j]]))
      df_groups$group <- factor(df_groups$group,
                                    levels = unique(df_groups$group[order(df_groups$median)]))

      # now treatment effect:
      # df_trt <- dplyr::bind_rows(
      #   lapply(models, function(x) {
      #     m <- mint(treatment_effect(x)$tau)
      #     m$group <- rownames(m)
      #     m
      #   }), .id = "model")

      # df <- rbind(df_groups, df_trt)
      df <- df_groups

      # Refer to global variables outside of ggplot context to pass CMD CHECK, see:
      # https://stackoverflow.com/questions/9439256/
      #   how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
      lci <- uci <- model <- group <- NULL

      comparison_plot <- ggplot2::ggplot(df, aes(x = group, y = median, ymin = lci, ymax = uci,
                                                 group = interaction(model),
                                                 color = model)) +
        # geom_jitter(size = 2) +
        geom_errorbar(size = 1.2, width = 0, position=position_dodge(width=0.5)) +
        geom_point(size = 2, stroke = 1.5, fill = "white", position=position_dodge(width=0.5), pch = 21) +
        coord_flip() +
        labs(x = "", y = "Treatment effect (95% interval)",
             title = effect_names[i]) +
        theme(legend.position="top")
      # plot(comparison_plot)
      return(comparison_plot)
    })
  }

  if(length(plots) == 1)
    plots <- plots[[1]]

  if(return_models_flag)
    return(list(plot = plots, models = models))
  else
    return(plots)
}
