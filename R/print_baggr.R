
#' S3 print method for objects of class `baggr` (model fits)
#'
#' This \code{print} method for a very concise summary of main model features.
#' More info is included in the summary of the model and its attributes.
#'
#' @param x object of class `baggr`
#' @param exponent if `TRUE`, results (for means) are converted to exp scale
#' @param digits Number of significant digits to print.
#' @param group logical; print group effects? If unspecified,
#'              they are printed if only
#'              less than 20 groups are present
#' @param ... currently unused by this package: further arguments passed
#'            to or from other methods (\code{print}  requirement)
#' @importFrom stats sd var median quantile
#' @importFrom crayon bold
#' @importFrom crayon red
#' @importFrom crayon blue
#' @export
#' @method print baggr
#'

print.baggr <- function(x, exponent=FALSE, digits = 2, group, ...) {
  ppd <- attr(x, "ppd")

  # Announce model type
  if(ppd) {
    cat("Model type: Prior predictive draws for",
        crayon::bold(model_names[x$model]), "\n")
  } else {
    cat("Model type:", crayon::bold(model_names[x$model]), "\n")
    cat("Pooling of effects:", crayon::bold(x$pooling), "\n")
  }

  cat("\n")

  if(length(x$effects) == 1)
    cat(crayon::bold(paste0("Aggregate treatment effect (on ", x$effects, "):\n")))
  else
    cat(crayon::bold(paste0("Aggregate treatment effect:\n")))


  if(x$pooling == "none") {
    cat("No treatment effect estimated as pooling = 'none'.\n\n")
  } else {
    # Means:
    if(exponent)
      te <- treatment_effect(x, transform = base::exp)
    else
      te <- treatment_effect(x)
    #trim=T avoids whitespace in place of minus sign
    if(x$model != "quantiles"){
      tau       <- format(mint(te[[1]]), digits = digits, trim = T)
      sigma_tau <- format(mint(te[[2]]), digits = digits, trim = T)
      if(exponent)
        cat("Exponent of hypermean (exp(tau))")
      else
        cat("Hypermean (tau)")
      cat(" = ", tau[2], "with 95% interval", tau[1], "to", tau[3], "\n")
      if(x$pooling == "partial" && !exponent){
        cat("Hyper-SD (sigma_tau) =", sigma_tau[2], "with 95% interval",
            sigma_tau[1], "to", sigma_tau[3], "\n")
        tot_pool <- format(heterogeneity(x)[,,1], digits = digits, trim = T)
        cat("Total pooling (1 - I^2) =", tot_pool[2], "with 95% interval",
            tot_pool[1], "to", tot_pool[3], "\n")
      }
    } else { #quantiles
      tau <- mint(te[[1]])
      sigma_tau <- mint(te[[2]])
      rownames(tau) <- rownames(sigma_tau) <- paste0(100*x$quantiles, "% quantile")
      print(tau, digits = digits)
      if(x$pooling == "partial"){
        cat(crayon::bold("\nSD of treatement effects:"))
        print(sigma_tau, digits = digits)
      }
    }
  }
  if(x$pooling == "full")
    cat("(SD(tau) undefined.)\n")
  cat("\n")

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
      pooling_tab <- pooling(x, summary = TRUE)
      if(exponent)
        study_eff_tab <- group_effects(x, summary = TRUE, transform=exp)
      else
        study_eff_tab <- group_effects(x, summary = TRUE)

      for(i in 1:dim(study_eff_tab)[3]){
        cat(paste0("Treatment effects on ", x$effects[i]))
        if(exponent){
          cat(" (converted to exp scale):\n")
          tab <- cbind(study_eff_tab[,c("mean", "lci", "uci"),i],
                       pooling = pooling_tab[2,,i])
        } else{
          cat(":\n")
          tab <- cbind(study_eff_tab[,c("mean", "sd"),i],
                       pooling = pooling_tab[2,,i])
        }
        print(tab, digits = digits)
      }
      cat("\n")
    }
  } else if(group_warning_flag) { #No printing of groups
    cat("Group effects omitted, as number of groups is > 20.\n",
        "Use print.baggr() with group = TRUE to print them.\n")
  }

  if(!is.null(x[["mean_lpd"]]))
    cat("Cross-validation result: mean lpd =",
        crayon::bold(format(x$mean_lpd)), "\n")

  invisible(x)
}
