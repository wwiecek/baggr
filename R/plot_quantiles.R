#' plot quantiles
#'
#' Plot results for baggr quantile models. Displays results facetted per group.
#' Results are `ggplot2` plots and can be modified.
#'
#' @param fit an object of class `baggr`
#' @param ncol number of columns for the plot; defaults to half of number of groups
#' @param hline logical; plots a line through 0
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' bg <- baggr(microcredit_simplified, model = "quantiles",
#'             quantiles = c(0.25, 0.50, 0.75),
#'             iter = 1000, refresh = 0,
#'             outcome = "consumption")
#' #vanilla plot
#' plot_quantiles(bg)[[1]]
#' plot_quantiles(bg, hline = TRUE)[[2]] +
#'   ggplot2::coord_cartesian(ylim = c(-2, 5)) +
#'   ggplot2::ggtitle("Works like a ggplot2 plot!")
#' }
#'
#' @import ggplot2

plot_quantiles <- function(fit, ncol, hline = TRUE) {
  if(!("baggr" %in% class(fit)) || (fit$model != "quantiles"))
    stop("fit must be a baggr 'quantiles' model")

  if(missing(ncol))
    ncol <- round((length(fit$quantiles) + 1)/2)

  # Group-specific treatment effect
  ste <- group_effects(fit, summary = TRUE)
  ate <- treatment_effect(fit, summary = TRUE)$tau
  df <- df_ate <- data.frame()
  for(i in seq_along(fit$quantiles)) {
    df2 <- as.data.frame(effect_slice(ste, i)[,c("lower", "mean", "upper"),
                                               drop = FALSE])
    df2$group <- rownames(df2)
    rownames(df2) <- NULL
    df2$quantile <- as.character(fit$quantiles[i])
    df <- rbind(df, df2)
    df_ate <- rbind(df_ate, as.data.frame(effect_slice(ate, i)[,
      c("lower", "mean", "upper"), drop = FALSE]))
  }

  # Average treatment effect
  df_ate$quantile <- as.character(fit$quantiles)
  df_ate$group <- "All countries"

  lower <- upper <- group <- NULL #check fix

  gg_groups <-
    ggplot2::ggplot(df, aes(y = mean, x=quantile, ymin = lower, ymax = upper)) +
    # theme_bw() +
    {if(hline) geom_hline(yintercept = 0, lty = "dashed")} +
    geom_line(aes(group = group)) +
    geom_errorbar(width = 0) +
    geom_point(size = 2, stroke = 1.5, pch = 21, fill = "white") +
    scale_x_discrete(labels = paste0(100*fit$quantiles, "%")) +
    ylab("mean treatment effect (95% interval)") +
    facet_wrap(~group, ncol = ncol) +
    theme(panel.grid.major.x = element_blank())

  gg_trt <-
    ggplot2::ggplot(df_ate, aes(y = mean, x=quantile, ymin = lower, ymax = upper)) +
    # theme_bw() +
    {if(hline) geom_hline(yintercept = 0, lty = "dashed")} +
    geom_line(aes(group = group)) +
    geom_errorbar(width = 0) +
    geom_point(size = 2, stroke = 1.5, pch = 21, fill = "white") +
    scale_x_discrete(labels = paste0(100*fit$quantiles, "%")) +
    ylab("mean treatment effect (95% interval)") +
    theme(panel.grid.major.x = element_blank())

  return(list(gg_trt, gg_groups))
}
