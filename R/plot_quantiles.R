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
#' \donttest{
#' bg <- baggr(microcredit_simplified, model = "quantiles",
#'             quantiles = c(0.25, 0.50, 0.75), iter = 1000,
#'             outcome = "consumption")
#' #vanilla plot
#' plot_quantiles(bg)
#'
#' plot_quantiles(bg, hline = TRUE) +
#'   ggplot2::coord_cartesian(ylim = c(-2, 2)) +
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
  ste <- apply(group_effects(fit), c(2,3),
               function(x) c(quantile(x, .025), "mean" = mean(x), quantile(x, .975)))
  dimnames(ste)[[3]] <- fit$quantiles
  df <- data.frame()
  for(i in 1:dim(ste)[3]) {
    df2 <- as.data.frame(t(ste[,,i]))
    df2$country <- rownames(df2)
    rownames(df2) <- NULL
    df2$quantile <- dimnames(ste)[[3]][i]
    df <- rbind(df, df2)
  }
  colnames(df) <- c("lci", "mean", "uci", "group", "quantile")

  # Average treatment effect
  ate <- apply(treatment_effect(fit)$tau, 2,
               function(x) c(quantile(x, .025), "mean" = mean(x), quantile(x, .975)))
  df_ate <- data.frame(t(ate))
  names(df_ate) <- c("lci", "mean", "uci")
  df_ate$quantile <- as.character(fit$quantiles)
  df_ate$group <- "All countries"

  lci <- uci <- group <- NULL #check fix

  gg_groups <-
    ggplot2::ggplot(df, aes(y = mean, x=quantile, ymin = lci, ymax = uci)) +
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
    ggplot2::ggplot(df_ate, aes(y = mean, x=quantile, ymin = lci, ymax = uci)) +
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
