
#' L'Abbe plot for binary data
#'
#' This plot shows relationship between proportions of events in control and treatment groups in binary data.
#'
#' @param data a data frame with binary data
#'             (must have columns `a`, `c`, `b`/`n1`, `d`/`n2`)
#' @param group a character string specifying group names (e.g. study names), used for labels;
#' @param plot_model if `TRUE`, then odds ratios and risk ratios [baggr] models are estimated (using default
#'                   settings) and their mean estimates of effects are plotted as lines
#' @param labels if `TRUE`, names from the `group` column are displayed
#' @param shade_se if `"none"`, nothing is plotted, if `"or"` or `"rr"`, a shaded area corresponding to
#'                 inverse of effect's (OR or RR) SE is added to each data point; the default is `"rr"`
#' @seealso `vignette("baggr_binary")` for an illustrative example
#' @return A `ggplot` object
#' @export
#' @import ggplot2
#'
labbe <- function(data, group = "group", plot_model = FALSE, labels = TRUE,
                  shade_se = c("rr", "or", "none")) {
  # Global bindings fix:
  risk_trt <- risk_control <- se_size <- p_trt <- p_bsl <- a <-b<-c<-d <- model <- NULL

  if(is.null(data[[group]])){
    data$group <- 1:nrow(data)
    labels <- FALSE
  }

  shade_se <- match.arg(shade_se, c("rr", "or", "none"))

  # Bring data to a common format (this will guarantee we have a,b,c,d cols)
  data_or <- prepare_ma(data, group=group, effect = "logOR", rare_event_correction = 0)
  data_rr  <- prepare_ma(data_or, effect = "logRR", rare_event_correction = 0)

  ggdata <-  data_or
  ggdata$risk_trt <- ggdata$a/(ggdata$a+ggdata$b)
  ggdata$risk_control <- ggdata$c/(ggdata$c+ggdata$d)

  if(shade_se == "or") ggdata$se_size <- 1 / data_or$se
  if(shade_se == "rr") ggdata$se_size <- 1 / data_rr$se

  if(plot_model) {
    model_or <- baggr(data_or, refresh=0, silent = TRUE)

    model_rr <- baggr(data_rr, refresh=0, silent = TRUE)

    # Apply trt effects to control pr events to get the right line
    te_or <- hypermean(model_or) #treatment_effect(model_or, summary = TRUE)$tau
    te_rr <- hypermean(model_rr) #treatment_effect(model_rr, summary = TRUE)$tau
    p_bsl <- seq(0,1,length=1000)
    p_trt_rr <- exp(te_rr[["mean"]])*p_bsl
    odds_trt <- exp(te_or[["mean"]])*(p_bsl/(1-p_bsl))
    p_trt_or <- odds_trt / (1 + odds_trt)
    p_lines <- rbind(
      data.frame(p_bsl, "p_trt" = p_trt_rr, model = "RR"),
      data.frame(p_bsl, "p_trt" = p_trt_or, model = "OR"))
    # p_trt_or <- te_or[["mean"]]
  }

  lim_min <- min(c(ggdata$risk_trt*0.9, ggdata$risk_control*0.9))
  lim_max <- max(c(ggdata$risk_trt*1.1, ggdata$risk_control*1.1))
  lim_max <- min(c(lim_max, 1))


  ggplot(ggdata, aes(x = risk_control, y = risk_trt)) +
    geom_point() +
    {if(shade_se != "none") geom_point(aes(size = se_size), alpha = .25)} +
    # ggforce::geom_circle(aes(x0 = risk_control, y0 = risk_trt, r=se_trt)) +
    {if(labels) geom_text(aes(label = group, hjust=-0.1, vjust=-0.1))} +
    {if(plot_model) geom_line(data = p_lines, aes(x=p_bsl, y=p_trt, lty = model))} +
    {if(plot_model) scale_linetype_manual(values = c("RR" = "dotted", "OR" = "dashed"))} +
    coord_cartesian(xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max)) +
    geom_abline(slope = 1, intercept = 0) +
    {if(shade_se != "none") scale_size_continuous(range = c(1,20))} +
    theme(legend.position = "bottom") +
    guides(size = "none") +
    labs(y = "Proportion of events in treatment arm",
         x = "Proportion of events in control arm") +
    {if(shade_se != "none")
      labs(subtitle = paste0("Shaded area is proportional in radius to 1/SE of treatment effect (log ",
                             toupper(shade_se), ")"))}


}
