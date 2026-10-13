#' Draw a forest plot for a baggr model
#'
#' The forest plot functionality in _baggr_ is a simple interface for
#' calling \link[forestplot]{forestplot} By default the forest plot
#' displays raw (unpooled) estimates for groups and the treatment effect
#' estimate underneath. This behaviour can be modified to display pooled
#' group estimates.
#'
#' @param bg a [baggr::baggr()] class object
#' @param show if `"inputs"`, then plotted points and lines
#'             correspond to raw inputs for each group;
#'             if `"posterior"` -- to posterior distribution;
#'             you can also plot `"both"` inputs and posteriors;
#'             if `"covariates"`, then fixed effect coefficients are plotted
#' @param print which values to print next to the plot: values of `"inputs"`
#'              or `"posterior"` means?
#'              (if `show="covariates"`, it must be `"posterior"`);
#'              the default for `show="both"` prints posterior values.
#' @param prob Probability mass of the central intervals, on the 0--1 scale.
#' @param digits number of digits to display when printing out mean and SD
#'        in the plot
#' @param ... other arguments passed to \link[forestplot]{forestplot}
#'
#' @seealso \link[forestplot]{forestplot} function and its vignette for examples;
#'          [effect_plot] and [baggr_plot] for non-forest plots of baggr results
#'
#' @examples
#' bg <- baggr(schools, iter = 500)
#' forest_plot(bg)
#' forest_plot(bg, show = "posterior", print = "inputs", digits = 2)
#'
#' @export
#' @import forestplot

forest_plot <- function(bg,
                        show = c("inputs", "posterior", "both", "covariates"),
                        print = show,
                        prob = .95,
                        digits = 3,
                        ...) {
  prob <- interval_prob(prob)
  if(!inherits(bg, "baggr"))
    stop("forest_plot can only be used with baggr objects")
  if(length(bg$effects) > 1)
    stop("forest_plot only works with 1-dimensional effects")

  # Get the summary-level data with tau and se columns
  if(bg$model %in% c("rubin_full", "logit"))
    ge_raw <- bg$summary_data
  else
    ge_raw <- bg$data
  if(bg$model %in% c("rubin_full", "mutau"))
    ge_raw$se <- ge_raw$se.tau

  ge_raw$mean  <- ge_raw$tau
  ge_raw$sd    <- ge_raw$se
  ge_raw$lci   <- ge_raw$tau - stats::qnorm((1 + prob)/2)*ge_raw$se
  ge_raw$uci   <- ge_raw$tau + stats::qnorm((1 + prob)/2)*ge_raw$se
  ge_posterior <- as.data.frame(
    effect_slice(group_effects(bg, summary = TRUE, prob = prob), 1))
  fe_posterior <- as.data.frame(
    effect_slice(fixed_effects(bg, summary = TRUE, prob = prob), 1))

  show <- match.arg(show, c("inputs", "posterior", "both", "covariates"))

  if(show == "inputs"){
    ge <- ge_raw
  }else if(show == "posterior"){
    ge <- ge_posterior
  }else if(show == "both"){
    ge <- list(mean = cbind(ge_raw$mean, ge_posterior$mean),
               sd   = cbind(ge_raw$sd,   ge_posterior$sd),
               lci  = cbind(ge_raw$lci,  ge_posterior$lci),
               uci  = cbind(ge_raw$uci,  ge_posterior$uci))
  } else if(show == "covariates"){
    if(nrow(fe_posterior) == 0)
      stop("No covariates to include in forest plot.")
    ge <- fe_posterior
  }

  # The default for both series displays posterior values, as before.
  if(print %in% c("both", "covariates"))
    print <- "posterior"
  print <- match.arg(print, c("inputs", "posterior"))
  if(show == "covariates" && print != "posterior")
    stop("Covariate plots require print = 'posterior'")
  printed <- if(show == "covariates") fe_posterior else
    if(print == "inputs") ge_raw else ge_posterior
  ge_printed <- paste0(format(printed$mean, digits = digits), " (",
                       format(printed$sd, digits = digits), ")")

  if(show == "covariates")
    n_ge_rows <- nrow(fe_posterior)
  else
    n_ge_rows <- nrow(ge_raw)

  te <- treatment_effect(bg)$tau

  te_printed <- paste0(format(mean(te), digits = digits), " (",
                       format(sd(te), digits = digits), ")")


  if(show == "covariates")
    fp_labels <- c("Covariate", bg$covariates, NA, "Treatment effect (hypermean)")
  else
    fp_labels <- c("Group mean treatment effect", as.character(group_names(bg)), NA, "Hypermean treatment effect")

  fp_printed_vals <- c("Mean (SD)", ge_printed, NA, te_printed)
  fp_text <- matrix(c(fp_labels, fp_printed_vals), n_ge_rows + 3, 2, byrow = FALSE)

  l <- list(...)
  if(!("hrzl_lines" %in% names(l)))
    l[["hrzl_lines"]] <- TRUE
  if(!("is.summary" %in% names(l)))
    l[["is.summary"]] <- c(TRUE, rep(FALSE, n_ge_rows), FALSE, FALSE)
  if(!("labeltext" %in% names(l)))
    l[["labeltext"]] <- fp_text

  # Fix to make X ticks bigger
  # (solving the world hunger this ain't)
  custom_xticks <- pretty(c(ge$lci, ge$uci), n = 8)
  # custom_xticks[which.min(abs(custom_xticks))] <- 0
  if(!("xticks" %in% names(l)))
    l[["xticks"]] <- round(sort(c(0, custom_xticks)))
  l[["txt_gp"]]  <- fpTxtGp(ticks = grid::gpar(cex = 0.85))

  l[["mean"]]  <- rbind(NA, as.matrix(ge$mean),NA, mint(te, prob = prob)["mean"])
  l[["lower"]] <- rbind(NA, as.matrix(ge$lci), NA, mint(te, prob = prob)[interval_names(prob)[1]])
  l[["upper"]] <- rbind(NA, as.matrix(ge$uci), NA, mint(te, prob = prob)[interval_names(prob)[2]])
  if(show == "both"){
    l[["fn.ci_norm"]]  <- c(fpDrawCircleCI, fpDrawNormalCI)
    l[["legend"]]      <- c("Input", "Estimate")
    l[["boxsize"]]     <- .15
  }
  if(show == "covariates"){
    l[["boxsize"]]     <- .15
  }
  if(!("xlab" %in% names(l)) && bg$effects != "mean")
    l[["xlab"]] <- bg$effects


  do.call(forestplot::forestplot, l)
}
