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
#'                (default) or (hyper-) `"effects"`. The former is not available
#'                when `what = "prior"`.
#' @param transform a function (e.g. exp(), log()) to apply to
#'                  the values of group (and hyper, if hyper=TRUE)
#'                  effects before plotting; when working with
#'                  effects that are on log scale,
#'                  exponent transform is used automatically,
#'                  you can plot on log scale by setting
#'                  transform = identity
#' @return an object of class `baggr_compare`
#' @seealso [plot.baggr_compare] and [print.baggr_compare] for working with results of this function
#' @author Witold Wiecek, Brice Green
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#' @export
#' @details If you pass parameters to the function you must specify
#' what kind of comparison you want, either "pooling" which
#' will run fully/partially/un-pooled models and compare them
#' or "prior" which will generate estimates without the data
#' and compare them to the model with the full data. For more
#' details see [baggr](baggr), specifically the PPD argument.
#' @examples \donttest{
#' # Most basic comparison between no, partial and full pooling
#' # (This will run the models)
#'
#' # run model with just prior and then full data for comparison
#' # with the same arguments that are passed to baggr
#' prior_comparison <-
#'     baggr_compare(schools,
#'                   model = 'rubin',
#'                   prior_hypermean = normal(0, 3),
#'                   prior_hypersd = normal(0,2),
#'                   prior_hypercor = lkj(2),
#'                   what = "prior")
#'
#' # print the aggregated treatment effects
#' prior_comparison
#'
#' # plot the comparison of the two distributions
#' plot(prior_comparison)
#'
#' # Now compare different types of pooling for the same model
#' pooling_comparison <-
#'    baggr_compare(schools,
#'                  model = 'rubin',
#'                  prior_hypermean = normal(0, 3),
#'                  prior_hypersd = normal(0,2),
#'                  prior_hypercor = lkj(2),
#'                  what = "pooling")
#'
#' # plot this comparison
#' plot(pooling_comparison)
#'
#' # Compare existing models:
#' bg1 <- baggr(schools, pooling = "partial")
#' bg2 <- baggr(schools, pooling = "full")
#' baggr_compare("Partial pooling model" = bg1, "Full pooling" = bg2,
#'               arrange = "grid")
#'
#' #' ...or simply draw prior predictive dist (note ppd=T)
#' bg1 <- baggr(schools, ppd=T)
#' bg2 <- baggr(schools, prior_hypermean = normal(0, 5), ppd=T)
#' baggr_compare("Prior A, p.p.d."=bg1,
#'               "Prior B p.p.d."=bg2,
#'               compare = "effects")
#'
#' # Compare posterior effects as a function of priors (note ppd=F)
#' bg1 <- baggr(schools, prior_hypersd = uniform(0, 20))
#' bg2 <- baggr(schools, prior_hypersd = normal(0, 5))
#' baggr_compare("Uniform prior on SD"=bg1,
#'               "Normal prior on SD"=bg2,
#'               compare = "effects")

#' # You can also compare different subsets of input data
#' bg1_small <- baggr(schools[1:6,], pooling = "partial")
#' baggr_compare("8 schools model" = bg1, "First 6 schools" = bg1_small)
#' }

baggr_compare <- function(...,
                          what    = "pooling",
                          compare = "groups",
                          transform = NULL) {
  l <- list(...)
  if(length(l) == 0)
    stop("Must provide baggr models or model specification.")
  if(!compare %in% c("groups","effects")){
    stop("'compare' argument must be set to either 'groups' or 'effects'.")
  }
  if(all(unlist(lapply(l, class)) == "baggr")) {
    # return_models_flag <- 0
    if(is.null(names(l)))
      names(l) <- paste("Model", 1:length(l))
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
                                  "silence_messages" = T,
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
                                  "silence_messages" = T,
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

  # Return treatment effects
  mean_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      mint(treatment_effect(x, transform = transform)$tau)
    })))
  sd_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      mint(treatment_effect(x, transform = transform)$sigma_tau)
    })))


  structure(list(
              models = models,
              mean_trt = mean_trt_effects,
              sd_trt = sd_trt_effects,
              compare = compare,
              effect_names = effect_names,
              transform = deparse(substitute(transform))),
            class = "baggr_compare")
}

#' Print method for baggr_compare models
#' @param x baggr_compare model
#' @param digits number of significant digits for effect estimates
#' @param ... other parameters passed to print
#' @export
print.baggr_compare <- function(x, digits, ...){
  cat("Mean treatment effects:\n")
  print(signif(x$mean_trt, digits = digits))
  cat("\n")
  cat("SD for treatment effects:\n")
  print(signif(x$sd_trt, digits = digits))
  cat("\n")
  cat(paste0("Transform: ", x$transform))
}

#' Plot method for baggr_compare models
#' @description Allows plots that compare multiple baggr models
#' that were passed for comparison purposes to baggr compare or
#' run automatically by baggr_compare
#' @param x baggr_compare model to plot
#' @param arrange If `"single"` (default), generate a single comparison plot;
#'                if `"grid"`, display multiple plots side-by-side.
#' @param style What kind of plot to display (if `arrange = "grid"`),
#'              passed to the `style` argument in [baggr_plot].
#' @param interval probability level used for display of posterior interval
#' @param hyper Whether to plot pooled treatment effect
#' in addition to group treatment effects
#' @param transform a function (e.g. exp(), log())
#' to apply to the values of group (and hyper, if hyper=TRUE)
#' effects before plotting; when working with effects that are on log scale, exponent transform is used automatically,
#' you can plot on log scale by setting transform = identity
#' @param order Whether to order by median treatment effect by group. If not, this
#' sorts group alphabetically. The pooled estimate is always listed first, when applicable.
#' @param ... ignored for now, may be used in the future
#' @export
plot.baggr_compare <- function(x,
                               style   = "areas",
                               arrange = "single",
                               interval = 0.95,
                               hyper = T,
                               transform = NULL,
                               order = F,
                               ...) {

  models <- x$models
  compare <- x$compare
  effect_names <- x$effect_names

  if(arrange == "grid") {
    plots <- lapply(models, baggr_plot,
                    style = style,
                    order = FALSE,
                    transform = transform)
    grid_width <- length(plots)
    # if each plots element contains multiple plots (like with quantiles):
    if(class(plots[[1]])[1] == "list")
      plots <- unlist(plots, recursive = FALSE)
    gridExtra::grid.arrange(grobs = plots, ncol = grid_width)
  }
  if(arrange == "single") {
    if(compare == "groups") {
      plots <- lapply(as.list(1:(length(effect_names))), function(i) {
        # Note: pipe operators are dplyr not used here for compatibility
        ll <- lapply(models, function(x) {
          # will need to be modified for quantiles models case:
          if(x$pooling != "none"){

            hyper_treat <- treatment_effect(x)[[i]]
            hyper_effects <- data.frame(
              lci = quantile(hyper_treat, (1 - interval)/2),
              median = quantile(hyper_treat, 0.5),
              uci = quantile(hyper_treat, 1 - (1 - interval)/2),
              mean = mean(hyper_treat),
              sd = sd(hyper_treat),
              group = "Pooled Estimate"
            )
            m <- as.data.frame(group_effects(x, interval = interval,
                                             summary = TRUE,
                                             transform = transform)[,,i])
            m$group <- rownames(m)
            m <- rbind(hyper_effects, m)
            m
          } else {
            m <- as.data.frame(group_effects(x,
                                             interval = interval,
                                             summary = TRUE,
                                             transform = transform)[,,i])
            m$group <- rownames(m)
            m
          }

        })
        df_groups <- data.frame()
        for(j in 1:length(ll))
          df_groups <- rbind(df_groups,
                             data.frame(model = names(ll)[j], ll[[j]]))

          if(order == T) {

            model_spread <- sapply(
              split(df_groups, df_groups$model), function(x) max(x$median) - min(x$median)
            )

            max_spread_model <- names(model_spread[which.max(model_spread)])

            rank_data <- df_groups[which(df_groups$model == max_spread_model),]

            lvls <- rank_data[order(rank_data$median),]$group

            ord <- c(
              "Pooled Estimate",
              rev(lvls)
            )

        } else {
          ord <- c(
            "Pooled Estimate",
            setdiff(sort(unique(df_groups$group)),
                    "Pooled Estimate")
          )
        }

        df_groups$group <- factor(df_groups$group,
                                  levels = rev(ord) # because of flipped coordinates
                                  )

        # df <- rbind(df_groups, df_trt)
        df <- df_groups

        # Refer to global variables outside of ggplot context to pass CMD CHECK, see:
        # https://stackoverflow.com/questions/9439256/
        #   how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
        lci <- uci <- model <- group <- NULL

        comparison_plot <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = median,
                                                            ymin = lci, ymax = uci,
                                                            group = interaction(model),
                                                            color = model)) +
          # geom_jitter(size = 2) +
          ggplot2::geom_errorbar(size = 1.2, width = 0,
                                 position = ggplot2::position_dodge(width=0.5)) +
          ggplot2::geom_point(size = 2, stroke = 1.5, fill = "white",
                              position = ggplot2::position_dodge(width=0.5),
                              pch = 21) +
          ggplot2::coord_flip() +
          ggplot2::labs(x = "", y = "Treatment effect (95% interval)",
                        title = paste0(
                          "Pooled and group-specific effect ",
                          "of treatment on ",
                          effect_names[i],
                          " outcome."
                          )
                        ) +
          baggr_theme_get() +
          ggplot2::theme(legend.position="top")
        return(comparison_plot)
      })
    } else if(compare == "effects"){
      plots <- do.call(effect_plot, models) +
        labs(fill = NULL)
    } else {
      stop("Argument compare = must be 'effects' or 'groups'.")
    }
  }

  if("ggplot" %in% class(plots)){
    return(plots)
  } else {
    structure(plots, class = "plot_list")
  }
}

#' Print list of baggr plots
#' @param x list of plots to print
#' @param ... ignored for now
#' @details prints plots in a loop, internal use only
print.plot_list <- function(x) {
  if(length(x) == 1) {
    print(x[[1]])
  } else {
    for(i in 1:length(plots)) {
      print(plots[[i]])
    }
  }
}


