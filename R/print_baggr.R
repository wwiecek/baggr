
#' S3 print method for objects of class `baggr` (model fits)
#'
#' This \code{print} method for a very concise summary of main model features.
#' More info is included in the summary of the model and its attributes.
#'
#' @param x object of class `baggr`
#' @param ... currently unused by this pacakge: further arguments passed
#'            to or from other methods (\code{print}  requirement)
#' @importFrom stats sd var median quantile
#' @importFrom crayon bold
#' @importFrom crayon red
#' @importFrom crayon blue
#' @export
#' @method print baggr
#'

print.baggr <- function(x, ...) {
  # cat(crayon::red("---this Baggr printing module is under construction---\n\n"))
  cat("Model type:", crayon::bold(model_names[x$model]), "\n")
  cat("Pooling of effects:", crayon::bold(x$pooling), "\n")
  cat("\n")

  cat(crayon::bold("Aggregate treatment effect:\n"))
  if(x$pooling == "none") {
    cat("No treatment effect estimated as pooling = 'none'.\n\n")
  } else {
    # Means:
    te <- treatment_effect(x)
    #trim=T avoids whitespace in place of minus sign
    if(x$model != "quantiles"){
      tau       <- format(mint(te[[1]]), digits = 2, trim = T)
      sigma_tau <- format(mint(te[[2]]), digits = 2, trim = T)
      cat("Mean(tau) =", tau[2], "with 95% interval", tau[1], "to", tau[3], "\n")
      if(x$pooling == "partial")
        cat("SD(tau) =", sigma_tau[2], "with 95% interval", sigma_tau[1], "to", sigma_tau[3], "\n")
    } else { #quantiles
      tau <- mint(te[[1]])
      sigma_tau <- mint(te[[2]])
      rownames(tau) <- rownames(sigma_tau) <- paste0(100*x$quantiles, "% quantile")
      print(tau, digits = 2)
      if(x$pooling == "partial"){
        cat(crayon::bold("\nSD of treatement effects:"))
        print(sigma_tau, digits = 2)
      }
    }
  }
  if(x$pooling == "full")
    cat("(SD(tau) undefined.)\n")
  cat("\n")

  if(x$pooling != "full") {
    # study_eff_tab <- apply(study_effects(x), c(2,3),
                             # function(x) c("mean" = mean(x), "sd" = sd(x)))
    study_eff_tab <- study_effects(x, summary = TRUE)
    # attach pooling metric:
    pooling_tab <- pooling(x, summary = TRUE)
    for(i in 1:dim(study_eff_tab)[3]){
      cat(paste0("Treatment effects on ", x$effects[i] , ":\n"))
      tab <- cbind(study_eff_tab[,c("mean", "sd"),i], pooling = pooling_tab[2,,i])
      print(tab, digits = 2)
    }
    cat("\n")
  }

  if(!is.null(x[["mean_lpd"]]))
    cat("Cross-validation result: mean lpd =", crayon::bold(format(x$mean_lpd)), "\n")

  invisible(x)
}
