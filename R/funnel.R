#' Funnel plots for baggr models
#'
#' Funnel plots are diagnostic visuals for meta-analyses: they plot study
#' standard errors (or precisions) on the vertical axis against estimated
#' treatment effects on the horizontal axis. Under the assumption of no
#' publication bias or heterogeneity-driven asymmetry, the points should
#' form a roughly symmetric inverted funnel around the pooled effect line.
#'
#' @param bg a [baggr()] model
#' @param show whether to plot raw study-level inputs (`"inputs"`) or
#'   posterior summaries (`"posterior"`)
#' @param level confidence level for reference lines
#' @param label logical: add study/group labels?
#' @param covariate optional name of a column in the input data to use for
#'   colouring points; ignored (with a warning) if the column is absent
#'
#' @details For baggr models that include covariates (meta-regressions), the
#'   plotted study effects already incorporate the covariate adjustments used
#'   in the model.
#'
#' @return A ggplot funnel plot for the supplied model
#' @export
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' bg <- baggr(schools, iter = 500, refresh = 0)
#' funnel_plot(bg, label = TRUE)
#'
funnel_plot <- function(bg,
                        show = c("inputs", "posterior"),
                        level = 0.95,
                        label = FALSE,
                        covariate = NULL) {
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

  covariate_values <- NULL
  if(!is.null(covariate)) {
    if(!covariate %in% names(studies)) {
      warning("Covariate column ", covariate, " not found; ignoring `covariate` argument.")
      covariate <- NULL
    } else {
      covariate_values <- studies[[covariate]]
    }
  }

  raw_df <- studies[, c("group", "tau", "se")]
  names(raw_df) <- c("group", "effect", "se")
  if(!is.null(covariate_values))
    raw_df$covariate_value <- covariate_values

  # posterior summaries (mean & SD) for partial/full pooling outputs
  post_df <- as.data.frame(group_effects(bg, summary = TRUE,
                                         interval = level)[,,1])
  post_df$group <- group_names(bg)
  post_df <- post_df[, c("group", "mean", "sd")]
  names(post_df) <- c("group", "effect", "se")
  if(!is.null(covariate_values))
    post_df$covariate_value <- covariate_values

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

  effect <- group <- head <- lower <- se <- tail <- upper <- covariate_value <- NULL

  aes_args <- list(x = quote(effect), y = quote(se))
  if(!is.null(covariate))
    aes_args$colour <- quote(covariate_value)

  ggplot2::ggplot(plot_df, do.call(ggplot2::aes, aes_args)) +
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
