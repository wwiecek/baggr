#' Draw a forest plot for a baggr model
#'
#' The forest plot functionality in _baggr_ is a simple interface for
#' calling the [forestplot] function. By default the forest plot
#' displays raw (unpooled) estimates for groups and the treatment effect
#' estimate underneath. This behaviour can be modified to display pooled
#' group estimates.
#'
#' @param bg a [baggr] class object
#' @param show if `inputs`, then plotted points and lines correspond to raw inputs;
#'             if `posterior` -- to posterior distribution; you can also plot `both`
#' @param print which values to print next to the plot: `inputs` or `posterior`?
#' @param prob width of the intervals (lines) for the plot
#' @param digits number of digits to display when printing out mean and SD
#'        in the plot
#' @param ... other arguments passed to [forestplot]
#'
#' @seealso [forestplot] function and its associated vignette for examples;
#'          [effect_plot] and [baggr_plot] for non-forest plots of baggr results
#'
#' @examples
#' bg <- baggr(schools, iter = 500)
#' forest_plot(bg)
#' forest_plot(bg, show = "posterior", print = "inputs", digits = 2)
#'
#' @export
#' @import forestplot

forest_plot <- function(bg, show = c("inputs", "posterior", "both"),
                        print = show,
                        prob = .95,
                        digits = 3,
                        ...) {
  if(!inherits(bg, "baggr"))
    stop("forest_plot can only be used with baggr objects")
  if(length(bg$effects) > 1)
    stop("forest_plot only works with 1-dimensional effects")

  # Get the summary-level data with tau and se columns
  if(bg$model %in% c("full", "logit"))
    ge_raw <- bg$summary_data
  else
    ge_raw <- bg$data
  if(bg$model %in% c("full", "mutau"))
    ge_raw$se <- ge_raw$se.tau

  ge_raw$mean  <- ge_raw$tau
  ge_raw$sd    <- ge_raw$se
  ge_raw$lci   <- ge_raw$tau - 1.96*ge_raw$se
  ge_raw$uci   <- ge_raw$tau + 1.96*ge_raw$se
  ge_posterior <- as.data.frame(
    group_effects(bg, summary = TRUE, interval = prob)[,,1])

  show <- match.arg(show, c("inputs", "posterior", "both"))

  if(show == "inputs")
    ge <- ge_raw
  else if(show == "posterior")
    ge <- ge_posterior
  else if(show == "both")
    ge <- list(mean = cbind(ge_raw$mean, ge_posterior$mean),
               lci  = cbind(ge_raw$lci,  ge_posterior$lci),
               uci  = cbind(ge_raw$uci,  ge_posterior$uci))
  te <- treatment_effect(bg)$tau
  if(print == "posterior")
    ge_printed <- paste0(format(ge_posterior$mean, digits = digits), " (",
                         format(ge_posterior$sd, digits = digits), ")")
  else
    ge_printed <- paste0(format(ge_raw$mean, digits = digits), " (",
                         format(ge_raw$sd, digits = digits), ")")

  te_printed <- paste0(format(mean(te), digits = digits), " (",
                       format(sd(te), digits = digits), ")")

  fp_labels <- c("Group", as.character(attr(bg$inputs, "group_label")), NA, "Hypermean")
  fp_printed_vals <- c("Mean (SD)", ge_printed, NA, te_printed)
  fp_text <- matrix(c(fp_labels, fp_printed_vals), nrow(ge_raw) + 3, 2, byrow = FALSE)

  l <- list(...)
  if(!("hrzl_lines" %in% names(l)))
    l[["hrzl_lines"]] <- TRUE
  if(!("is.summary" %in% names(l)))
    l[["is.summary"]] <- c(TRUE, rep(FALSE, nrow(ge_raw)), FALSE, FALSE)
  if(!("labeltext" %in% names(l)))
    l[["labeltext"]] <- fp_text
  l[["mean"]]  <- rbind(NA, as.matrix(ge$mean),NA, mint(te)[2])
  l[["lower"]] <- rbind(NA, as.matrix(ge$lci), NA, mint(te)[1])
  l[["upper"]] <- rbind(NA, as.matrix(ge$uci), NA, mint(te)[3])
  if(show == "both"){
    l[["fn.ci_norm"]]  <- c(fpDrawCircleCI, fpDrawNormalCI)
    l[["legend"]]      <- c("Input", "Estimate")
    l[["boxsize"]]     <- .2
  }
  if(!("xlab" %in% names(l)) && bg$effects != "mean")
    l[["xlab"]] <- bg$effects


  do.call(forestplot, l)
}
