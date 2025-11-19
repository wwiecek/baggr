#' Funnel plots for baggr models
#'
#' @param bg a [baggr()] model
#' @param show whether to plot raw study-level inputs (`"inputs"`) or
#'   posterior summaries (`"posterior"`)
#' @param level confidence level for reference lines
#' @param label logical: add study/group labels?
#'
#' @return A ggplot funnel plot for the supplied model
#' @export
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' bg <- baggr(schools, iter = 500, refresh = 0)
#' funnel(bg, label = TRUE)
#'
funnel <- function(bg,
                   show = c("inputs", "posterior"),
                   level = 0.95,
                   label = FALSE) {
  stopifnot(inherits(bg, "baggr"))
  if(length(bg$effects) != 1)
    stop("Funnel plot currently defined for 1-dimensional effects.")
  show <- match.arg(show)

  # raw study-level data (tau/se/group)
  if(bg$model %in% c("rubin_full", "logit")) {
    studies <- bg$summary_data
  } else {
    studies <- bg$data
  }
  if(bg$model %in% c("rubin_full", "mutau"))
    studies$se <- studies$se.tau
  if(is.null(studies$group))
    studies$group <- group_names(bg)

  raw_df <- studies[, c("group", "tau", "se")]
  names(raw_df) <- c("group", "effect", "se")

  # posterior summaries (mean & SD) for partial/full pooling outputs
  post_df <- as.data.frame(group_effects(bg, summary = TRUE,
                                         interval = level)[,,1])
  post_df$group <- group_names(bg)
  post_df <- post_df[, c("group", "mean", "sd")]
  names(post_df) <- c("group", "effect", "se")

  plot_df <- if(show == "inputs") raw_df else post_df

  # hyper-mean + heterogeneity
  if(bg$pooling == "none")
    stop("Need a pooled model to plot average treatment effect for the funnel.")
  hypermean <- hypermean(bg)[["mean"]]
  hetero <- if(bg$pooling == "partial") hypersd(bg)[["mean"]] else 0

  crit <- stats::qnorm(1 - (1 - level)/2)
  se_grid <- seq(0, max(plot_df$se), length.out = 200)
  fan <- data.frame(
    se = se_grid,
    lower = hypermean - crit * sqrt(se_grid^2 + hetero^2),
    upper = hypermean + crit * sqrt(se_grid^2 + hetero^2)
  )

  effect <- group <- head <- lower <- se <- tail <- upper <- NULL

  ggplot2::ggplot(plot_df, ggplot2::aes(x = effect, y = se)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = fan,
                       ggplot2::aes(x = lower, y = se),
                       linetype = "dashed") +
    ggplot2::geom_line(data = fan,
                       ggplot2::aes(x = upper, y = se),
                       linetype = "dashed") +
    ggplot2::geom_vline(xintercept = hypermean, linewidth = 0.4) +
    {if(label) ggrepel::geom_text_repel(ggplot2::aes(label = group),
                                        min.segment.length = 0)} +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = bg$effects,
                  y = "Standard error") +
    baggr_theme_get()
}
