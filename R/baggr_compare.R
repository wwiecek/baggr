#' (Run and) compare multiple baggr models
#'
#' @description Compare multiple [baggr] models by either
#' providing multiple already existing models as (named) arguments or
#' passing parameters necessary to run a [baggr] model.
#'
#' @param ... Either some (at least 1) objects of class `baggr`
#'            (you should name your objects, see the example below)
#'            or the same arguments you'd pass to [baggr].
#'            In the latter case you must specify `what` to compare.
#' @param what  One of `"pooling"` (comparison between no, partial and
#'              full pooling) or `"prior"` (comparison between prior and
#'              posterior predictive). If pre-existing baggr models are
#'              passed to `...`, this argument is ignored.
#' @param compare When plotting, choose between comparison of `"groups"`
#'                (default), `"hyperpars"` (to omit group-specific estimates)
#'                or (predicted) `"effects"`.
#'                The `"groups"` option is not available when `what = "prior"`.
#' @param transform a function (e.g. exp(), log()) to apply to
#'                  the the sample of group (and hyper, if `hyper=TRUE`)
#'                  effects before plotting; when working with
#'                  effects that are on log scale,
#'                  exponent transform is used automatically,
#'                  you can plot on log scale by setting
#'                  transform = identity
#' @param prob  Width of uncertainty interval (defaults to 95%)
#' @param plot logical; calls [plot.baggr_compare] when running `baggr_compare`
#' @return an object of class `baggr_compare`
#' @seealso [plot.baggr_compare] and [print.baggr_compare]
#'          for working with results of this function
#' @author Witold Wiecek, Brice Green
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#' @export
#' @details If you pass parameters to the function you must specify
#' what kind of comparison you want, either `"pooling"`, which
#' will run fully/partially/un-pooled models and then compare them,
#' or `"prior"` which will generate estimates without the data
#' and compare them to the model with the full data. For more
#' details see [baggr], specifically the `ppd` argument.
#' @examples \donttest{
#' # Most basic comparison between no, partial and full pooling
#' # (This will run the models)
#' # run model with just prior and then full data for comparison
#' # with the same arguments that are passed to baggr
#' prior_comparison <-
#'   baggr_compare(schools,
#'                 model = 'rubin',
#'                 #this is just for illustration -- don't set it this low normally!
#'                 iter = 500,
#'                 prior_hypermean = normal(0, 3),
#'                 prior_hypersd = normal(0,2),
#'                 prior_hypercor = lkj(2),
#'                 what = "prior")
#' # print the aggregated treatment effects
#' prior_comparison
#' # plot the comparison of the two distributions
#' plot(prior_comparison)
#' # Now compare different types of pooling for the same model
#' pooling_comparison <-
#'   baggr_compare(schools,
#'                 model = 'rubin',
#'                 #this is just for illustration -- don't set it this low normally!
#'                 iter = 500,
#'                 prior_hypermean = normal(0, 3),
#'                 prior_hypersd = normal(0,2),
#'                 prior_hypercor = lkj(2),
#'                 what = "pooling",
#'                 # You can automatically plot:
#'                 plot = TRUE)
#' # Compare existing models (you don't have to, but best to name them):
#' bg1 <- baggr(schools, pooling = "partial")
#' bg2 <- baggr(schools, pooling = "full")
#' baggr_compare("Partial pooling model" = bg1, "Full pooling" = bg2)
#'
#' #' ...or simply draw from prior predictive dist (note ppd=T)
#' bg1 <- baggr(schools, ppd=TRUE)
#' bg2 <- baggr(schools, prior_hypermean = normal(0, 5), ppd=TRUE)
#' baggr_compare("Prior A, p.p.d."=bg1,
#'               "Prior B p.p.d."=bg2,
#'               compare = "effects")
#'
#' # Compare how posterior predictive effect varies with e.g. choice of prior
#' bg1 <- baggr(schools, prior_hypersd = uniform(0, 20))
#' bg2 <- baggr(schools, prior_hypersd = normal(0, 5))
#' baggr_compare("Uniform prior on SD"=bg1,
#'                    "Normal prior on SD"=bg2,
#'                    compare = "effects", plot = TRUE)
#'
#' # Models don't have to be identical. Compare different subsets of input data:
#' bg1_small <- baggr(schools[1:6,], pooling = "partial")
#' baggr_compare("8 schools model" = bg1, "First 6 schools" = bg1_small,
#'               plot = TRUE)
#' }

baggr_compare <- function(...,
                          what    = "pooling",
                          compare = c("groups", "hyperpars", "effects"),
                          transform = NULL,
                          prob = 0.95,
                          plot = FALSE) {
  l <- list(...)
  if(length(l) == 0)
    stop("Must provide baggr models or model specification.")
  compare <- match.arg(compare, c("groups", "hyperpars", "effects"))

  if(all(unlist(lapply(l, class)) == "baggr_cv")) {
    message("LOO CV models used instead of baggr models.",
    "Using full models.")
    l <- lapply(l, function(x) x$full_model)
  }

  if(all(unlist(lapply(l, class)) == "baggr")) {
    # return_models_flag <- 0
    if(is.null(names(l)))
      names(l) <- paste("Model", 1:length(l))
    if(length(unique(names(l))) != length(names(l)))
      stop("You must use unique model names")
    models <- l
  } else {
    # return_models_flag <- 1
    if(what == "pooling"){
      if("pooling" %in% names(l))
        stop("Can't run the model comparison with pooling setting",
             "already set to a particular value.")
      models <- lapply(list("none", "partial", "full"), function(pool){
        # message to display progress
        # message(paste0("Sampling for model with pooling set to ", pool))

        # suppress baggr/rstan output
        model <- do.call(baggr, c(l, "pooling" = pool,
                                  "silent" = TRUE,
                                  "refresh" = 0))

        # return model
        model
      })
      names(models) <- c("none", "partial", "full")
    } else if(what == "prior") {
      if("ppd" %in% names(l))
        stop("Can't run the model comparison with ppd setting",
             "already set to a particular value.")
      models <- lapply(list(TRUE, FALSE), function(ppdv){
        check_which <- ifelse(ppdv, "just the prior", "prior and full data")
        message(paste0("Sampling for model with ", check_which, "."))
        model <- do.call(baggr, c(l, "ppd" = ppdv,
                                  "silent" = TRUE,
                                  "refresh" = 0))
        model
      })
      names(models) <- c("Prior", "Posterior")
      compare <- "effects"
    }
  }

  # For PPD models always switch comparison to "effects"
  if(any(unlist(lapply(models, attr, "ppd")))) {
    if(compare != "effects")
      message("Models to compare are PPD -- switching to compare='effects'")
    compare <- "effects"
  }

  effect_names <- lapply(models, function(x) x$effects)
  # quite a mouthful:
  if(!all(
    unlist(
      lapply(
        effect_names,
        function(x) all.equal(effect_names[[1]], x))) == 1)
  )
    stop("Models must have the same effects to be comparable")

  effect_names <- effect_names[[1]]



  # Check if, for group comparisons, covariates are the same
  cov_effects <- NULL
  covs <- lapply(models, function(x) x$covariates)
  if(length(unique(covs)) > 1)
    message("Compared models use different covariates. ",
            "Hyperparameters are not directly comparable.")
  # else if(length(effect_names) == 1 && length(unique(unlist(covs))) > 0) {
  if(length(effect_names) == 1 && length(unique(unlist(covs))) > 0) {
    # Grab covariates, if there are any:
    cov_effects <- lapply(models, function(x) {
      est <- fixed_effects(x, transform = transform,
                           summary = TRUE)
      if(dim(est)[1] == 0)
        # return(c(NA, NA))
        return(NULL)
      else
        return(
          data.frame(
            model = "",
            covariate = rownames(est),
            mean = est[, "mean", 1],
            sd = est[, "sd", 1],
            lci = est[, "lci", 1],
            uci = est[, "uci", 1]
          )
        )
    })
    cov_effects <- cov_effects[!unlist(lapply(cov_effects, is.null))]
    for(i in 1:length(cov_effects)){
      cov_effects[[i]]$model <- names(models)[i]
      # rownames(cov_effects[[i]]) <- NULL
    }
    cov_effects <- do.call(rbind, cov_effects)
    rownames(cov_effects) <- NULL
  }

  # Return treatment effects, hyperSDs, predicted effects
  mean_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      est <- treatment_effect(x, transform = transform,
                              interval = prob,
                              summary = TRUE, message = FALSE)$tau
      if(is.matrix(est)) {
        if(nrow(est) == 1) est <- est[1,]
      }
      est
    })))
  sd_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      est <- treatment_effect(x, transform = transform,
                              interval = prob,
                              summary = TRUE, message = FALSE)$sigma_tau
      if(is.matrix(est)) {
        if(nrow(est) == 1) est <- est[1,]
      }
      est
    })))
  posteriorpd_trt <- do.call(rbind, (
    lapply(models, function(x) {
      est <- effect_draw(x, transform = transform,
                         interval = prob,
                         summary = TRUE)
      if(is.matrix(est)) {
        if(nrow(est) == 1) est <- est[1,]
      }
      est
    })))

  bgc <- structure(
    list(
      models = models,
      mean_trt = mean_trt_effects,
      sd_trt = sd_trt_effects,
      posteriorpd_trt = posteriorpd_trt,
      covariates = cov_effects,
      compare = compare,
      effect_names = effect_names,
      transform = transform,
      prob = prob
    ),
    class = "baggr_compare")

  if(plot)
    print(plot(bgc))

  bgc
}

#' Print method for baggr_compare models
#' @param x baggr_compare model
#' @param digits number of significant digits for effect estimates
#' @param ... other parameters passed to print
#' @export
print.baggr_compare <- function(x, digits, ...){
  cat("\nMean treatment effects:\n")
  print(signif(x$mean_trt, digits = digits))
  cat("\nSD for treatment effects:\n")
  print(signif(x$sd_trt, digits = digits))
  cat("\nPosterior predictive effects:\n")
  print(signif(x$posteriorpd_trt, digits = digits))

  if(!is.null(x$covariates)){
    cat("\nMean (SD) for covariates:\n")
    mcov <- x$covariates
    mcov$meansd <- paste0(signif(mcov$mean, digits = digits), " (",
                          signif(mcov$sd, digits = digits), ")",
                          sep = "")
    mcov <- mcov[, c("model", "covariate", "meansd")]
    ycov <- reshape(mcov, idvar = "model", timevar = "covariate", direction = "wide")
    colnames(ycov) <- gsub("meansd.", "", colnames(ycov))
    print(ycov, row.names = FALSE)
  }

  if(!is.null(x$transform))
    cat(paste0("\nTransform: ", deparse(x$transform)))
}

#' Plot method for baggr_compare models
#' @description Allows plots that compare multiple baggr models
#' that were passed for comparison purposes to baggr compare or
#' run automatically by baggr_compare
#' @param x baggr_compare model to plot
#' @param compare When plotting, choose between comparison of `"groups"`
#'                (default), `"hyperpars"` (to omit group-specific estimates)
#'                or (predicted) `"effects"`.
#'                The `"groups"` option is not available when `what = "prior"`.
#' @param grid_models If `FALSE` (default), generate a single comparison plot;
#'                if `TRUE`, display each model (using individual [baggr_plot]'s)
#'                side-by-side.
#' @param grid_parameters if `TRUE`, uses `ggplot`-style facetting when plotting models
#'                        with many parameters (especially `"quantiles"`, `"sslab"`);
#'                        if `FALSE`, returns separate plot for each parameter
#' @param style What kind of plot to display (if `grid_models = TRUE`),
#'              passed to the `style` argument in [baggr_plot].
#' @param prob  Width of uncertainty interval (defaults to 95%)
#' @param hyper Whether to plot pooled treatment effect
#'              in addition to group treatment effects when `compare = "groups"`
#' @param transform a function (e.g. exp(), log())
#'                  to apply to the values of group (and hyper, if hyper=TRUE)
#'                  effects before plotting
#' @param order Whether to sort by median treatment effect by group.
#'              If yes, medians from the model with largest range of estimates
#'              are used for sorting.
#'              If not, groups are shown alphabetically.
#' @param add_values logical; if TRUE, values will be printed next to the plot,
#'                   in a style that's similar to what is done for forest plots
#' @param values_digits number of significant digits to use when printing values,
#' @param values_size size of font for the values, if `add_values == TRUE`
#' @param vline logical; show vertical line through 0 in the plot?
#' @param ... ignored for now, may be used in the future
#' @export
plot.baggr_compare <- function(x,
                               compare = x$compare,
                               style   = "areas",
                               grid_models = FALSE,
                               grid_parameters = TRUE,
                               prob = x$prob,
                               hyper = TRUE,
                               transform = NULL,
                               order = F,
                               vline = FALSE,
                               add_values = FALSE,
                               values_digits = 2,
                               values_size = 2,
                               ...) {

  # Refer to global variables outside of ggplot context to pass CMD CHECK, see:
  # https://stackoverflow.com/questions/9439256/
  #   how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  lci <- uci <- model <- group <- NULL

  compare <- match.arg(compare, c("groups", "hyperpars", "effects"))

  if(is.null(transform))
    transform <- x$transform

  models <- x$models
  effect_names <- x$effect_names

  if(grid_models) {
    if(length(effect_names) > 1)
      stop("Cannot plot models with more than 1 effect in a grid")
    plots <- lapply(models, baggr_plot,
                    style = style,
                    order = FALSE,
                    transform = transform,
                    vline = vline)
    grid_width <- length(plots)
    # if each plots element contains multiple plots (like with quantiles):
    if(class(plots[[1]])[1] == "list")
      plots <- unlist(plots, recursive = FALSE)

    for(i in 1:length(models))
      plots[[i]] <- plots[[i]] + ggplot2::ggtitle(names(models)[i])

    return(gridExtra::grid.arrange(grobs = plots, ncol = grid_width))
  }

  # This is non-grid, baggr_compare-specific code:

  if(compare %in% c("groups", "hyperpars")) {

    # Create input data frames for ggplots
    plot_dfs <- lapply(as.list(1:(length(effect_names))), function(i) {

      # Note: pipe operators or dplyr not used here anymore
      # to reduce dependencies
      effects_list <- lapply(models, function(cmodel) {
        m <- as.data.frame(group_effects(cmodel,
                                         interval = prob,
                                         summary = TRUE,
                                         transform = transform)[,,i])
        m$group <- rownames(m)

        if(cmodel$pooling != "none" && hyper) {
          if(length(effect_names) == 1)
            hyper_treat <- treatment_effect(cmodel,
                                            transform = transform,
                                            message=FALSE)$tau
          if(length(effect_names) > 1)
            hyper_treat <- treatment_effect(cmodel,
                                            transform = transform,
                                            message=FALSE)$tau[,i]
          hyper_effects <- data.frame(
            lci = quantile(hyper_treat, (1 - prob)/2),
            median = quantile(hyper_treat, 0.5),
            uci = quantile(hyper_treat, 1 - (1 - prob)/2),
            mean = mean(hyper_treat),
            sd = sd(hyper_treat),
            group = "Pooled Estimate"
          )
          if(compare == "hyperpars")
            m <- hyper_effects
          if(compare == "groups")
            m <- rbind(hyper_effects, m)
        }
        return(m)
      })

      df_groups <- data.frame()
      for(j in 1:length(effects_list))
        df_groups <- rbind(
          df_groups,
          data.frame(model = names(effects_list)[j], effects_list[[j]])
        )

      # Find the ordering for groups
      if(order == T) {
        ord <- get_order(df_groups, hyper)
      } else if(hyper) {
        ord <- c(setdiff(unique(df_groups$group),
                         "Pooled Estimate"),
                 "Pooled Estimate")
      } else {
        ord <- unique(df_groups$group)
      }

      # Rearrange ordering of groups/hypereffects
      df_groups$group <- factor(df_groups$group,
                                levels = rev(ord) # because of flipped coordinates
      )

      df_groups
    })

    names(plot_dfs) <- effect_names

    # This is an unpleasant way of doing map2...
    # Sorry to everyone.
    if(grid_parameters){
      big_df <- data.frame()
      for(i in 1:length(plot_dfs))
        big_df <- rbind(big_df,
                        data.frame(parameter = effect_names[i], plot_dfs[[i]]))

      plots <- single_comp_plot(big_df, "", grid = T,
                                ylab = paste0("Treatment effect (",
                                              as.character(100*prob), "% interval)"),
                                add_values = add_values,
                                values_digits = values_digits,
                                values_size = values_size)
    } else {
      plots <- lapply(as.list(1:(length(effect_names))), function(i) {
        single_comp_plot(plot_dfs[[i]], effect_names[i], grid = F,
                         ylab = paste0("Treatment effect (",
                                       as.character(100*prob), "% interval)"),
                         add_values = add_values,
                         values_digits = values_digits,
                         values_size = values_size)
      })
    }
  }

  if(compare == "effects"){
    arglist <- models
    arglist[["transform"]] <- transform
    plots <- do.call(effect_plot, arglist) + ggplot2::labs(fill = NULL)
  }

  # return the plots
  plots
}



#' Plot single comparison ggplot in `baggr_compare` style
#'
#' @param df data.frame with columns `group`, `median`, `lci`, `uci`,
#'           `model` (character or factor listing compared models) and,
#'           optionally, `parameter` (character or factor with name of parameter)
#' @param title `ggtitle` argument passed to ggplot
#' @param legend `legend.position`  argument passed to ggplot
#' @param ylab Y axis label
#' @param grid logical; if `TRUE`, facets the plot by values in the `parameter` column
#' @param points you can optionally specify a (`numeric`) column that has values of points
#'               to be plotted next to intervals
#' @param add_values logical; if `TRUE`, values will be printed next to the plot,
#'                   in a style that's similar to what is done for forest plots
#' @param values_digits number of significant digits to use when printing values,
#' @param values_size size of font for the values, if `add_values == TRUE`
#' @return a `ggplot2` object
#' @export
single_comp_plot <- function(df, title="", legend = "top", ylab = "", grid = F,
                             points = FALSE,
                             add_values = FALSE, values_digits = 2, values_size = 2.5) {

  group <- median <- lci <- uci <- model <- NULL

  # A bit of code to avoid printing pointless colour legend
  if(is.null(df[["model"]])) df$model <- ""
  if(length(unique(df[["model"]])) == 1)
    legend <- "none"

  pl <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = median,
                                         ymin = lci, ymax = uci,
                                         group = interaction(model),
                                         color = model)) +
    {if(grid) ggplot2::facet_wrap( ~ parameter, ncol = 3)} +
    ggplot2::geom_errorbar(linewidth = 1.2, width = 0,
                           position = ggplot2::position_dodge(width = 0.5)
    ) +
    # use position_dodge2 to reverse order of points to match legend
    ggplot2::geom_point(size = 2, stroke = 1.5, fill = "white",
                        position = ggplot2::position_dodge2(width=0.5, reverse=TRUE),
                        pch = 21) +
    { if(points != 0 && !is.null(df[[points]]) && is.numeric(df[[points]]))
      ggplot2::geom_point(aes(x = group, y = .data[[points]]), color = "black") } +
    { if(!add_values) ggplot2::coord_flip() } +
    ggplot2::labs(x = "",
                  y = ylab) +
    {if(title != "") ggplot2::ggtitle(paste0("Effect of treatment on ", title)) } +
    baggr_theme_get() +
    ggplot2::theme(legend.position=legend)

  if(add_values)
    pl <- add_plot_values(pl, values_digits, values_size)

  pl
}


add_plot_values <- function(pl, values_digits = 2, values_size = 2.5) {

  fmti <- function(x, digits = values_digits) {
    format(round(x, digits), nsmall = digits)
  }

  group <- median <- lci <- uci <- model <- NULL

  pl <- pl +
    ggplot2::geom_text(
      aes(label = fmti(lci, values_digits),
          hjust = 1,
          y = 1.1*max(uci)),
      position = ggplot2::position_dodge(width = .5),
      size = values_size) +
    ggplot2::geom_text(
      aes(label = fmti(mean, values_digits),
          hjust = 1,
          y = 1.2*max(uci)),
      position = ggplot2::position_dodge(width = .5),
      size = values_size) +
    ggplot2::geom_text(
      aes(label = fmti(uci, values_digits),
          hjust = 1,
          y = 1.3*max(uci)),
      position = ggplot2::position_dodge(width = .5),
      size = values_size) +
    ggplot2::theme(strip.text.x = ggplot2::element_blank()) +
    ggplot2::coord_flip(clip = "off")

  pl
}



#' Separate out ordering so we can test directly
#' @param df_groups data.frame of group effects used in [plot.baggr_compare]
#' @param hyper show parameter estimate? same as in [plot.baggr_compare]
#' @details Given a set of effects measured by models, identifies the
#' model which has the biggest range of estimates and ranks groups
#' by those estimates, returning the order

get_order <- function(df_groups, hyper) {
  model_spread <- sapply(
    split(df_groups, df_groups$model), function(x) max(x$median) - min(x$median)
  )
  max_spread_model <- names(model_spread[which.max(model_spread)])
  rank_data <- df_groups[which(df_groups$model == max_spread_model),]
  lvls <- setdiff(as.character(rank_data[order(rank_data$median),]$group),
                  "Pooled Estimate")
  if(hyper)
    ord <- c(rev(lvls), "Pooled Estimate")
  else
    ord <- rev(lvls)
  ord
}
