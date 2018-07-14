
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
#' @details
#' TBC
#'

print.baggr <- function(x, ...) {
  cat(crayon::red("---this Baggr printing module is under construction---\n\n"))
  cat("Model type:", crayon::bold(x$model), "\n")
  cat("Pooling of effects:", crayon::bold(x$pooling), "\n")
  cat("\n")


  cat("Aggregate treatment effect:\n")
  if(x$pooling == "none") {
    cat("No treatment effect estimated as pooling = 'none'.\n")
  } else {
    te <- treatment_effect(x)
    tau <- te[[1]]; sigma_tau <- te[[2]]
    cat("Mean(tau) = ", round(mean(tau), 2), "; 95% interval", round(quantile(tau, .025),2), "to", round(quantile(tau, .975), 2))
    cat("\n")
    if(x$pooling == "partial")
      cat("SD(tau) = ", round(mean(sigma_tau), 2), "; 95% interval",
          round(quantile(sigma_tau, .025), 2), "to", round(quantile(sigma_tau, .975), 2), "\n")
    if(x$pooling == "full")
      cat("(SD(tau) undefined.)\n")
    cat("\n")
  }
  if(x$pooling != "full") {
    cat("Study effects:\n")
    study_eff_tab <- t(apply(study_effects(x), 2,
                             function(x) c("mean" = mean(x), "sd" = sd(x))))
    names(dimnames(study_eff_tab))[1] <- ""
    # attach pooling metric:
    study_eff_tab <- cbind(study_eff_tab, pooling(x)[,2])

    colnames(study_eff_tab) <- c("mean", "SD", "pooling")
    print(study_eff_tab, digits = 2)
    cat("\n")
  }


  invisible(x)
}

