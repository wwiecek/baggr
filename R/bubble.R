#' Draw a bubble plot for a baggr object
#'
#' @param bg          A [baggr] object with at least one covariate
#' @param covariate  Character: a name of one of the covariates present in `bg`
#' @param pred  TRUE/FALSE; print the linear trend for the variable?
#' @param label   A logical value that determines if group labels are added to the plot.
#' @return A bubble plot (`ggplot`)
#' @examples
#'
#' schools$a <- rnorm(8)
#' bg <- baggr(schools, cov = "a")
#' bubble(bg, "a")
#'
#' @details
#'
#' The bubble plot shows the observed outcomes of individual
#' studies plotted against study-level covariates. The size of each point reflects
#' the precision of the study (inverse variance).
#'
#' @export
#' @import ggrepel
bubble <- function(bg, covariate, pred=TRUE, label=TRUE) {
  if(!inherits(bg, "baggr"))
    stop("bg must be a baggr object")
  if(is.null(bg$covariates))
    stop("bg must include covariates")
  if(!(covariate %in% bg$covariate))
    stop("requested covariate not specified in the baggr object")

  .covariate <- dotsize <- NULL

  data <- bg$data
  data$dotsize <- 1/(data$se^2)
  data$mean <- as.data.frame(group_effects(bg, summary = TRUE)[,,1])$mean
  group <- NULL
  data$.covariate <- data[[covariate]]
  fe <- fixed_effects(bg, summary=TRUE)[,"mean",1]
  te <- treatment_effect(bg,summary=TRUE)$tau[["mean"]]

  ggplot2::ggplot(data, aes(x = .covariate,y=mean)) +
    ggplot2::geom_point(aes(size=dotsize)) +
    ggplot2::guides(size = "none") +
    ggplot2::xlab(covariate) + ggplot2::ylab(bg$effects) +
    # {if(pred) ggplot2::geom_smooth(method="lm", level = interval, col="black")} +
    {if(pred) ggplot2::geom_abline(slope = fe, intercept = te, lty = "dashed")} +
    {if(label) ggrepel::geom_text_repel(aes(x = .covariate,
                                            y=mean,
                                            label=group),
                                            box.padding=1)}

}
