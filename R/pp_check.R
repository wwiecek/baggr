#' Posterior predictive checks for baggr model
#'
#' Performs posterior predictive checks with the
#' \pkg{bayesplot} package.
#'
#' @param x Model to check
#' @param type type of pp_check. For a list see
#'             \pkg{\link[bayesplot:available_ppc]{here}}.
#' @param nsamples number of samples to compare
#' @aliases pp_check
#'
#' @details For a detailed explanation of each of the ppc functions,
#' see the \code{\link[bayesplot:PPC-overview]{PPC}}
#' documentation of the `bayesplot`
#' package.
#'
#' @import bayesplot
#' @importFrom utils getFromNamespace
#' @export

pp_check.baggr <- function(x, type = "dens_overlay", nsamples = 40) {
  pp_fun <- utils::getFromNamespace(paste0("ppc_",type),ns = "bayesplot")
  col <- switch(x$model,
                rubin = "tau",
                mutau = "tau",
                quantiles = stop_not_implemented(),
                full = stop_not_implemented()
  )
  y <- x$data[,col]
  yrep <- predict(x, nsamples = nsamples)
  pp_fun(y, yrep)
}
