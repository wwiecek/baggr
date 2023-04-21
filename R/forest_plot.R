#' Draw a forest plot for a baggr model
#'
#' The forest plot functionality in _baggr_ is a simple interface for
#' calling the [forestplot] function. By default the forest plot
#' displays raw (unpooled) estimates for groups and the treatment effect
#' estimate underneath. This behaviour can be modified to display pooled
#' group estimates.
#'
#' @param bg a [baggr] class object
#' @param show if `"inputs"`, then plotted points and lines
#'             correspond to raw inputs for each group;
#'             if `"posterior"` -- to posterior distribution;
#'             you can also plot `"both"` inputs and posteriors;
#'             if `"covariates"`, then fixed effect coefficients are plotted
#' @param print which values to print next to the plot: values of `"inputs"`
#'              or `"posterior"` means?
#'              (if `show="covariates"`, it must be `"posterior"`)
#' @param prob width of the intervals (lines) for the plot
#' @param digits number of digits to display when printing out mean and SD
#'        in the plot
#' @param ... other arguments passed to [forestplot]
#'
#' @seealso [forestplot] function and its vignette for examples;
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
  ge_raw$lci   <- ge_raw$tau - 1.96*ge_raw$se
  ge_raw$uci   <- ge_raw$tau + 1.96*ge_raw$se
  ge_posterior <- as.data.frame(
    group_effects(bg, summary = TRUE, interval = prob)[,,1])
  fe_posterior <- as.data.frame(
    fixed_effects(bg, summary = TRUE, interval = prob)[,,1])

  show <- match.arg(show, c("inputs", "posterior", "both", "covariates"))

  if(show == "inputs"){
    ge <- ge_raw
    ge_printed <- paste0(format(ge_raw$mean, digits = digits), " (",
                         format(ge_raw$sd, digits = digits), ")")
  }else if(show == "posterior"){
    ge <- ge_posterior
    ge_printed <- paste0(format(ge_posterior$mean, digits = digits), " (",
                         format(ge_posterior$sd, digits = digits), ")")
  }else if(show == "both"){
    ge <- list(mean = cbind(ge_raw$mean, ge_posterior$mean),
               sd   = cbind(ge_raw$sd,   ge_posterior$sd),
               lci  = cbind(ge_raw$lci,  ge_posterior$lci),
               uci  = cbind(ge_raw$uci,  ge_posterior$uci))
    ge_printed <- paste0(format(ge_posterior$mean, digits = digits), " (",
                         format(ge_posterior$sd, digits = digits), ")")
  } else if(show == "covariates"){
    if(nrow(fe_posterior) == 0)
      stop("No covariates to include in forest plot.")
    ge <- fe_posterior
    ge_printed <- paste0(format(fe_posterior$mean, digits = digits), " (",
                         format(fe_posterior$sd, digits = digits), ")")
  }

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
    fp_labels <- c("Group mean treatment effect", as.character(attr(bg$inputs, "group_label")),
                   NA, "Hypermean treatment effect")

  fp_printed_vals <- c("Mean (SD)", ge_printed, NA, te_printed)
  fp_text <- matrix(c(fp_labels, fp_printed_vals), n_ge_rows + 3, 2, byrow = FALSE)

  l <- list(...)
  if(!("hrzl_lines" %in% names(l)))
    l[["hrzl_lines"]] <- TRUE
  if(!("is.summary" %in% names(l)))
    l[["is.summary"]] <- c(TRUE, rep(FALSE, nrow(ge_raw)), FALSE, FALSE)
  if(!("labeltext" %in% names(l)))
    l[["labeltext"]] <- fp_text

  # Fix to make X ticks bigger
  # (solving the world hunger this ain't)
  custom_xticks <- pretty(c(ge$lci, ge$uci), n = 8)
  # custom_xticks[which.min(abs(custom_xticks))] <- 0
  if(!("xticks" %in% names(l)))
    l[["xticks"]] <- round(sort(c(0, custom_xticks)))
  l[["txt_gp"]]  <- fpTxtGp(ticks = grid::gpar(cex = 0.85))

  l[["mean"]]  <- rbind(NA, as.matrix(ge$mean),NA, mint(te)[2])
  l[["lower"]] <- rbind(NA, as.matrix(ge$lci), NA, mint(te)[1])
  l[["upper"]] <- rbind(NA, as.matrix(ge$uci), NA, mint(te)[3])
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


  do.call(forestplot, l)
}

#' Draw a bubble plot using ggplot syntax 
#'
#' A bubble plot shows the inverse variance (1/SE^2) as a function of the mean treatment effect
#' for given covariates. The size of the bubble corresponds to the magnitude of the inverse
#' variance. 
#'
#' @param bg          A [baggr] object.
#' @param covariates  Covariates for the [baggr] data. 
#' @param regression  A logical value that sets an overlying linear regression and 
#'                    confidence interval on or off.
#' @param interval    The confidence interval used in the linear regression, 
#'                    given as a value between 0 and 1.
#' @param add_label   A logical value that determines if group labels are added to the plot. 
#' @return A bubble plot of the inverse variance of the mean treatment effect as a function  
#'         of the given covariates.  
#' @examples
#' 
#' covariates <- rnorm(8)
#' bg <- baggr(schools)
#' bubble_plot(bg,covariates)
#'
#' @details
#' 
#' The bubble plot shows the observed outcomes (e.g. mean or logOR) of individual
#' studies plotted against quantitative covariates. The size of the points reflect
#' the wieght the studies have in the meta analysis, measured by inverse variance. 
#' A regression and 95% confidence interval shows the predicted mean treatment effect
#' as a function of the covariate.
#'
#' @export
#' @import ggrepel
bubble_plot <- function(bg, covariates, regression=TRUE, interval=0.95,add_label=TRUE) {
  data <- bg$data
  data$se <- 1/(data$se^2)
  data$mean <- as.data.frame(group_effects(bg, summary = TRUE, interval = interval)[,,1])$mean

  group <- NULL
  
  ggplot2::ggplot(data, aes(x = covariates,y=mean)) + 
  ggplot2::geom_point(aes(size=data$se)) +
  ggplot2::scale_size_continuous(name = "1/SE^2") + 
  ggplot2::xlab("Covariates") + ggplot2::ylab(bg$effects) +
  {if(regression) ggplot2::geom_smooth(method="lm",level=interval, col="black")} +
  {if(add_label) ggrepel::geom_text_repel(aes(x = covariates,y=mean,label=data$group),box.padding=1)}
  
}
