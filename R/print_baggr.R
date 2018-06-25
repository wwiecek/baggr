
#' S3 print method for objects of class `baggr` (model fits)
#'
#' This \code{print} method for a very concise summary of main model features.
#' More info is included in the summary of the model and its attributes.
#'
#' @param bg object of class `baggr`
#' @importFrom crayon bold
#' @importFrom crayon red
#' @importFrom crayon blue
#' @export
#' @method print baggr
#' @details
#' TBC
#'

print.baggr <- function(bg) {
  cat(crayon::red("---this Baggr printing module is under construction---\n\n"))
  cat("Model type:", crayon::bold(bg$model), "\n")
  cat("Pooling of effects:", crayon::bold(bg$pooling), "\n")
  cat("\n")


  cat("Aggregate treatment effect:\n")
  te <- treatment_effect(bg)
  tau <- te[[1]]; sigma_tau <- te[[2]]
  cat("Mean(tau) = ", round(mean(tau), 2), "; 95% interval", round(quantile(tau, .025),2), "to", round(quantile(tau, .975), 2))
  cat("\n")
  cat("SD(tau) = ", round(mean(sigma_tau), 2), "; 95% interval", round(quantile(sigma_tau, .025), 2), "to", round(quantile(sigma_tau, .975), 2))
  cat("\n\n")

  if(bg$pooling != "full") {
    cat("Study effects:\n")
    study_eff_tab <- t(apply(study_effects(bg), 2,
                             function(x) c("mean" = mean(x), "sd" = sd(x))))
    names(dimnames(study_eff_tab))[1] <- ""
    # attach pooling metric:
    study_eff_tab <- cbind(study_eff_tab, pooling(bg)[,2])

    colnames(study_eff_tab) <- c("mean", "SD", "pooling")
    print(study_eff_tab, digits = 2)
    cat("\n")
  }


  invisible(bg)
}

