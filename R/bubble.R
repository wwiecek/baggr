#' Bubble plots for meta-regression models
#'
#' @param bg a [baggr()] model using summary-level data, with covariates
#' @param covariate one of the covariates present in the model
#' @param fit logical: show mean model prediction? (slope is mean estimate of [fixed_effects()], intercept is [hypermean()])
#' @param label logical: label study/group names?
#'
#' @return A simple bubble plot in `ggplot` style.
#' Dot sizes are proportional to inverse of variance of each study (more precise studies are larger).
#' @seealso [labbe()] for an exploratory plot of binary data in similar style
#' @export
#'
bubble <- function(bg, covariate, fit=TRUE, label=TRUE) {
  if(!inherits(bg, "baggr"))
    stop("bg must be a baggr object")
  if(is.null(bg$covariates))
    stop("bg must include covariates")
  if(!(covariate %in% bg$covariate))
    stop("requested covariate not specified in the baggr object")

  if(!(bg$model %in% c("rubin", "mutau")))
    stop("At the moment bubble plots are only available for summary-level data models")

  data <- bg$data
  data$dotsize <- 1/(data$se^2)
  data$mean <- as.data.frame(group_effects(bg, summary = TRUE)[,,1])$mean
  group <- NULL
  data$.covariate <- data[[covariate]]
  fe <- fixed_effects(bg, summary=TRUE)[,"mean",1]
  te <- treatment_effect(bg,summary=TRUE)$tau[["mean"]]
  if(is.null(data$group))
    data$group <- group_names(bg)
  data$.group <- data$group

  ggplot2::ggplot(data, aes(x = .covariate,y=mean)) +
    ggplot2::geom_point(aes(size=dotsize)) +
    ggplot2::guides(size = "none") +
    ggplot2::xlab(covariate) + ggplot2::ylab(bg$effects) +
    # {if(pred) ggplot2::geom_smooth(method="lm", level = interval, col="black")} +
    {if(prediction) ggplot2::geom_abline(slope = fe, intercept = te, lty = "dashed")} +
    {if(label) ggrepel::geom_text_repel(aes(x = .covariate,
                                            y=mean,
                                            label=.group),
                                        box.padding=1)}

}
