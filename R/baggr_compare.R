#' (Run and) compare multiple baggr models
#'
#' @param ... Either a number (at least 1) of objects of class `baggr`
#'            (you should name your objects, see the example below)
#'            or the same arguments you'd pass to baggr(),
#'            in which you must case specify `what` to compare.
#' @param what  One of `"pooling"` (comparison between no, partial and
#'              full pooling) or `"prior"` (comparison between prior and
#'              posterior predictive). If pre-existing baggr models are
#'              passed to `...`, this argument is ignored.
#' @param compare When plotting, choose between comparison of `"groups"`
#'                (default) or (hyper-) `"effects"`. The former is not available
#'                when `what = "prior"`.
#' @param arrange If `"single"` (default), generate a single comparison plot;
#'                if `"grid"`, display multiple plots side-by-side.
#' @param style What kind of plot to display (if `arrange = "grid"`),
#'              passed to the `style` argument in [baggr_plot].
#' @return a `ggplot` is rendered and/or returned
#' @author Witold Wiecek
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' # Most basic comparison between no, partial and full pooling
#' # (This will run the models)
#' baggr_compare(schools)
#'
#' # Compare prior vs posterior
#' baggr_compare(schools, what = "prior")
#'
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
                          style   = "areas",
                          arrange = "single") {
  l <- list(...)
  if(length(l) == 0)
    stop("Must provide baggr models or model specification.")
  if(all(unlist(lapply(l, class)) == "baggr")) {
    return_models_flag <- 0
    if(is.null(names(l)))
      names(l) <- paste("Model", 1:length(l))
    models <- l
  } else {
    return_models_flag <- 1
    if(what == "pooling"){
      if("pooling" %in% names(l))
        stop("Can't run the model comparison with pooling setting",
             "already set to a particular value.")
      models <- lapply(list("none", "partial", "full"), function(pool){
        try(do.call(baggr, c(l, "pooling" = pool)))
      })
      names(models) <- c("none", "partial", "full")
    } else if(what == "prior") {
      if("ppd" %in% names(l))
        stop("Can't run the model comparison with ppd setting",
             "already set to a particular value.")
      models <- lapply(list(TRUE, FALSE), function(ppdv){
        try(do.call(baggr, c(l, "ppd" = ppdv)))
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
  if(!all(unlist(lapply(effect_names, function(x) all.equal(effect_names[[1]], x))) == 1))
    stop("Models must have the same effects to be comparable")
  effect_names <- effect_names[[1]]

  if(arrange == "grid") {
    plots <- lapply(models, baggr_plot, style = style, order = FALSE)
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
          m <- as.data.frame(group_effects(x, summary = TRUE)[,,i])
          m$group <- rownames(m)
          m
        })
        df_groups <- data.frame()
        for(j in 1:length(ll))
          df_groups <- rbind(df_groups,
                             data.frame(model = names(ll)[j], ll[[j]]))
        df_groups$group <- factor(df_groups$group,
                                  levels = unique(df_groups$group[order(df_groups$median)]))

        # df <- rbind(df_groups, df_trt)
        df <- df_groups

        # Refer to global variables outside of ggplot context to pass CMD CHECK, see:
        # https://stackoverflow.com/questions/9439256/
        #   how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
        lci <- uci <- model <- group <- NULL

        comparison_plot <- ggplot2::ggplot(df, aes(x = group, y = median, ymin = lci, ymax = uci,
                                                   group = interaction(model),
                                                   color = model)) +
          # geom_jitter(size = 2) +
          geom_errorbar(size = 1.2, width = 0, position=position_dodge(width=0.5)) +
          geom_point(size = 2, stroke = 1.5, fill = "white", position=position_dodge(width=0.5), pch = 21) +
          coord_flip() +
          labs(x = "", y = "Treatment effect (95% interval)",
               title = effect_names[i]) +
          theme(legend.position="top")
        # plot(comparison_plot)
        return(comparison_plot)
      })
    } else if(compare == "effects"){
      plots <- do.call(effect_plot, models)
    } else {
      stop("Argument compare= must be 'effects' or 'groups'.")
    }
  }

  if(length(plots) == 1)
    plots <- plots[[1]]

  # Print treatment effects
  mean_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      mint(treatment_effect(x)$tau)
    })))
  sd_trt_effects <- do.call(rbind, (
    lapply(models, function(x) {
      mint(treatment_effect(x)$sigma_tau)
    })))

  cat("Mean treatment effects:\n")
  print(mean_trt_effects)
  cat("\n")
  cat("SD for treatment effects:\n")
  print(sd_trt_effects)


  if(return_models_flag)
    return(list(plot = plots, models = models))
  else
    return(plots)
}
