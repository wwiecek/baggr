
#' S3 print method for objects of class `baggr` (model fits)
#'
#' This prints a concise summary of the main [baggr::baggr()] model features.
#' More info is included in the summary of the model and its attributes.
#' When the model supports publication selection (currently `"rubin"`), the
#' printout also reports the posterior summaries of the relative publication
#' probability across the |z|-intervals defined in the selection model.
#'
#' @param x object of class `baggr`
#' @param exponent if `TRUE`, results (for means) are converted to exp scale
#' @param digits Number of significant digits to print.
#' @param prob Probability mass of the central interval (defaults to 0.95).
#' @param group logical; print group effects? If unspecified,
#'              they are printed only if
#'              less than 20 groups are present
#' @param fixed logical: print fixed effects?
#' @param ... currently unused by this package: further arguments passed
#'            to or from other methods (\code{print}  requirement)
#' @importFrom stats sd var median quantile
#' @importFrom crayon bold
#' @importFrom crayon red
#' @importFrom crayon blue
#' @importFrom utils head tail
#' @export
#' @method print baggr
#'

print.baggr <- function(x,
                        exponent=FALSE,
                        digits = 2,
                        prob = 0.95,
                        group,
                        fixed = TRUE, ...) {
  prob <- interval_prob(prob)
  ppd <- attr(x, "ppd")
  is_selection_model <- x$model %in% c("rubin") &&
    !is.null(x$inputs$M) && x$inputs$M > 0 && !is.null(x$inputs$c)

  # Announce model type
  if(ppd) {
    cat("Model type: Prior predictive draws for",
        crayon::bold(model_names[x$model]), "\n")
  } else {
    if(!is_selection_model)
      cat("Model type:", crayon::bold(model_names[x$model]), "\n")
    else
      cat("Model type:", crayon::bold(model_names[x$model]), "and with selection on |z|\n")

    cat("Pooling of effects:", crayon::bold(x$pooling), "\n")
  }

  cat("\n")

  if(length(x$effects) == 1)
    cat(crayon::bold(paste0("Aggregate treatment effect (on ",
                            x$effects, "), ",
                            x$n_groups," groups:\n")))
  else
    cat(crayon::bold(paste0("Aggregate treatment effects, ",
                            x$n_groups," groups:\n")))

  if(x$pooling == "none") {
    cat("No treatment effect estimated as pooling = 'none'.\n\n")
  } else {
    # Means:
    if(exponent){
      te <- treatment_effect(x, transform = base::exp)
      de <- effect_draw(x, transform = base::exp)
    }else{
      te <- treatment_effect(x)
      de <- effect_draw(x)
    }
    if(exponent)
      cat("Exponent of hypermean (exp(theta))")
    else
      cat("Hypermean (theta)")

    intervaltxt <- paste0("with ", as.character(100*prob), "% interval")

    #trim=T avoids whitespace in place of minus sign
    if(length(x$effects) == 1){
      tau <- format(interval_stats(mint(te$tau, prob = prob), prob),
                         digits = digits, trim = TRUE)
      sigma_tau <- format(interval_stats(mint(te$sigma_tau, prob = prob), prob),
                         digits = digits, trim = TRUE)
      posteriorpred <- format(interval_stats(mint(de,  prob = prob), prob),
                         digits = digits, trim = TRUE)
      cat(" = ", tau["mean"], intervaltxt, tau["lower"], "to", tau["upper"], "\n")
      if(x$pooling == "partial"){
        if(!exponent)
          cat("Hyper-SD (sigma_theta) =", sigma_tau["mean"], intervaltxt,
            sigma_tau["lower"], "to", sigma_tau["upper"], "\n")
        cat("Posterior predictive effect =", posteriorpred["mean"], intervaltxt,
            posteriorpred["lower"], "to", posteriorpred["upper"], "\n")
        tot_pool <- format(interval_stats(effect_slice(heterogeneity(x, prob = prob), 1)[1,], prob),
                           digits = digits, trim = TRUE)
        if(!ppd)
          cat("Total pooling (1 - I^2) =", tot_pool["mean"], intervaltxt,
              tot_pool["lower"], "to", tot_pool["upper"], "\n")
      }

      if(x$model == "mutau" && x$pooling != "full") {
        mutauc <- format(interval_stats(effect_slice(mutau_cor(x, summary = TRUE, prob = prob), 1)[1,], prob),
                         digits = digits, trim = TRUE)
        cat("Correlation of treatment effect and baseline =", mutauc["mean"], intervaltxt,
            mutauc["lower"], "to", mutauc["upper"], "\n")
      }

    } else {
      cat("\n")
      tau <- interval_stats(
        mint(te$tau, prob = prob, median = TRUE, sd = TRUE), prob)
      sigma_tau <- interval_stats(
        mint(te$sigma_tau, prob = prob, median = TRUE, sd = TRUE), prob)

      if(x$model == "quantiles"){
        rownames(tau) <- paste0(100*x$quantiles, "% quantile")
        if(!exponent)
          rownames(sigma_tau) <- paste0(100*x$quantiles, "% quantile")
      } else if(x$model == "sslab") {
        rownames(tau) <- x$effects
        if(!exponent)
          rownames(sigma_tau) <- x$effects
      }

      print(interval_table(tau, prob, include_sd = !exponent), digits = digits)
      if(!exponent && x$pooling == "partial"){
        cat(crayon::bold("\nSD of treatement effects:"))
        print(interval_table(sigma_tau, prob), digits = digits)


      }
    }
  }
  if(x$pooling == "full")
    cat("(SD(theta) undefined.)\n")
  cat("\n")

  if(is_selection_model) {
    omega_summary <- selection(x, prob = prob)
    if(is.null(dim(omega_summary)))
      omega_summary <- matrix(omega_summary, nrow = 1, dimnames = list(NULL, names(omega_summary)))
    cuts <- x$inputs$c
    fmt_cut <- function(val) formatC(val, format = "f", digits = 2, drop0trailing = TRUE)
    lower_bounds <- c(0, utils::head(cuts, -1))
    upper_bounds <- cuts
    omega_summary <- interval_stats(omega_summary, prob)
    intervaltxt_selection <- paste0("with ", as.character(100*prob), "% interval")
    cat("Publication probability relative to |z| in (", fmt_cut(utils::tail(cuts, 1)),
        ", Inf):\n", sep = "")
    for(i in seq_len(nrow(omega_summary))) {
      mean_val <- format(omega_summary[i, "mean"], digits = digits, trim = TRUE)
      lower_val <- format(omega_summary[i, "lower"], digits = digits, trim = TRUE)
      upper_val <- format(omega_summary[i, "upper"], digits = digits, trim = TRUE)
      cat("* |z| in (", fmt_cut(lower_bounds[i]), ", ", fmt_cut(upper_bounds[i]), "] = ",
          mean_val, " ", intervaltxt_selection, " ", lower_val, " to ",
          upper_val, "\n", sep = "")
    }
    cat("\n")
  }

  # If this is just drawing from prior, stop here
  if(ppd)
    return(invisible(x))

  # Group effects
  # Check if we should print them:
  group_warning_flag <- FALSE
  if(missing(group)){
    group_warning_flag <- TRUE
    group <- ifelse(x$n_groups > 20, FALSE, TRUE)
  }



  if(group){ #Print groups
    if(x$pooling != "full") {
      # study_eff_tab <- apply(group_effects(x), c(2,3),
      # function(x) c("mean" = mean(x), "sd" = sd(x)))
      pooling_tab <- pooling(x, summary = TRUE, prob = prob)
      if(exponent)
        study_eff_tab <- group_effects(x, summary = TRUE,
                                       transform=exp, prob = prob,
                                       rename_int = FALSE)
      else
        study_eff_tab <- group_effects(x, summary = TRUE,
                                       prob = prob,
                                       rename_int = FALSE)
      if(length(x$effects) == 1 && x$n_groups > 1)
        cat(paste0("Group-specific treatment effects"))

      for(i in 1:length(x$effects)){
        if(length(x$effects) > 1 || x$n_groups == 1){
          if(x$n_groups == 1)
            cat(paste0("Group-specific treatment effect on ", x$effects[i]))
          else
            cat(paste0("Treatment effects on ", x$effects[i]))
        }

        if(exponent)
          cat(" (converted to exp scale):\n")
        else
          cat(":\n")
        # Keep the study dimension even for a single study.
        tab <- effect_slice(study_eff_tab, i)
        tab <- interval_table(tab, prob, include_sd = !exponent)
        tab <- cbind(tab, pooling = effect_slice(pooling_tab, i)[,"mean"])
        print(tab, digits = digits)
      }
      cat("\n")
    }
  } else if(group_warning_flag) { #No printing of groups
    cat("Group effects omitted, as number of groups is > 20.",
        "\nUse print.baggr() with group = TRUE to print them.\n")
  }

  if(fixed && length(x$covariates) > 0) {
    # fe <- fixed_effects(x, summary = TRUE, prob = prob)

    if(exponent)
      fixed_eff_tab <- fixed_effects(x, summary = TRUE, transform=exp, prob = prob)
    else
      fixed_eff_tab <- fixed_effects(x, summary = TRUE, prob = prob)

    for(i in 1:length(x$effects)){
      cat(paste0("Covariate (fixed) effects on ", x$effects[i]))
      if(exponent)
        cat(" (converted to exp scale):\n")
      else
        cat(":\n")
      tab <- effect_slice(fixed_eff_tab, i)
      tab <- interval_table(tab, prob, include_sd = !exponent)
      print(tab, digits = digits)
      cat("\n")
    }
  }

  if(!is.null(x[["mean_lpd"]]))
    cat("Cross-validation result: -2 * mean log predictive density =",
        crayon::bold(format(x$mean_lpd)), "\n")

  invisible(x)
}
