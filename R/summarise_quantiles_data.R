#' summarise_quantiles_data
#'
#' Given individual level data, return list of summary statistics
#' Of quantile means and Sigma's, as well as K, N.
#' This operation is done automatically inside [baggr] function
#' but it can be used for exploring input data.
#'
#' @param data      a data.frame with grouping, outcome, treatment variables
#' @param quantiles a vector of quantiles (between 0 and 1)
#' @param outcome   character: column name in `data` for outcome variable
#' @param group     character: column name in `data` for grouping variable
#' @param treatment character: column name in `data` for treatment variable
#' @param means_only logical; if true, covariances and extra outputs are not
#'                   returned, but only means (`y_0` and `y_1`)
#' @details
#' This function is intended for data exploration outside of [baggr]
#' function. When calling [baggr], individual-level data should be supplied --
#' they will be summarised automatically.
#' (See [baggr] for documentation of how columns should be formatted.)
#'
#' Estimates are obtained externally via `quantreg::rq` function of the
#' `quantreg` package. This implementation is experimental and will change in the
#' future versions of the package.
#'
#' @return
#' A list with mean estimates for control and treatment (`y_0` and `y_1`) and
#' variance-covariance matrices (`Sigma_y_k_0` and `Sigma_y_k_1`).
#'
#' @examples
#' \donttest{
#' summarise_quantiles_data(microcredit_simplified, c(.2, .4, .6),
#'                          outcome = "consumption")
#' }
#' @importFrom quantreg rq
#' @import stats
#' @export

summarise_quantiles_data <- function(data, quantiles,
                                     outcome   = "outcome",
                                     group  = "group",
                                     treatment = "treatment",
                                     means_only = FALSE) {
  N <- length(quantiles)
  K <- length(unique(data[[group]]))
  if(!all(quantiles[-1] - quantiles[-length(quantiles)] > 0))
    stop("Vector of quantiles should be strictly monotonic")
  if(N < 2)
    stop("At least 2 quantiles needed")
  # if(K < 2)
  # stop("Needs at least 2 groups to summarise data for")

  check_columns(data, outcome, group, treatment)

  # Calculate means and SE's of our quantiles via quantreg::qr()
  # Not very elegant & very slow.
  groups <- unique(data[[group]])

  data <- data[,c(group, outcome, treatment)]
  names(data) <- c("group", "outcome", "treatment")

  # for(group in groups) {
  #   qr <- quantreg::rq(data[[outcome]][data[[group]] == group] ~ data[[treatment]][data[[group]] == group],
  #                tau = quantiles)
  #   y_0 <- qr$coef[1,]
  #   y_1 <- qr$coef[2,]
  # }
  calc <- stats::aggregate(outcome ~ group + treatment, function(x) {
    qr <- quantreg::rq(x ~ 1, tau = quantiles)
    qs <- summary(qr, se = "iid")
    unlist(lapply(qs, function(coef) {
      c("mean" = coef$coef[1,1],
        "se"   = coef$coef[1,2])
    }))
  }, data = data)

  # Proceed carefully, as the ordering of columns is crucial here!
  y_0    <- calc$outcome[calc$treatment == 0, 2*(1:length(quantiles)) - 1, drop = F]
  y_0_se <- calc$outcome[calc$treatment == 0, 2*(1:length(quantiles)), drop = F]
  y_1    <- calc$outcome[calc$treatment == 1, 2*(1:length(quantiles)) - 1, drop = F]
  y_1_se <- calc$outcome[calc$treatment == 1, 2*(1:length(quantiles)), drop = F]
  colnames(y_0) <- colnames(y_1) <- quantiles
  rownames(y_0) <- calc$group[calc$treatment == 0]
  rownames(y_1) <- calc$group[calc$treatment == 1]

  # Obtain Sigma's from SE's (RM's code rewritten in R)
  Sigma_y_k_0 <- array(0, c(K, N, N))
  Sigma_y_k_1 <- array(0, c(K, N, N))
  density_eval_k_0 <- matrix(0, K, N)
  density_eval_k_1 <- matrix(0, K, N)
  for(k in 1:K){
    for (i in 1:N){
      density_eval_k_0[k,i] <- sqrt(quantiles[i]*(1-quantiles[i])/(y_0_se[k,i]^2))
      density_eval_k_1[k,i] <- sqrt(quantiles[i] *(1-quantiles[i])/(y_1_se[k,i]^2))
      Sigma_y_k_0[k,i,i] <- (y_0_se[k,i]^2)
      Sigma_y_k_1[k,i,i] <- (y_1_se[k,i]^2)
    }
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma_y_k_0[k,i,j] <- quantiles[i]*(1-quantiles[j])/(density_eval_k_0[k,i]*density_eval_k_0[k,j])
        Sigma_y_k_1[k,i,j] <- quantiles[i]*(1-quantiles[j])/(density_eval_k_1[k,i]*density_eval_k_1[k,j])
        Sigma_y_k_0[k,j,i] <- Sigma_y_k_0[k,i,j]
        Sigma_y_k_1[k,j,i] <- Sigma_y_k_1[k,i,j]
      }
    }
  }

  # Output:
  if(means_only)
    return(list(
      y_0 = y_0,
      y_1 = y_1
    ))
  else
    return(list(y_0 = y_0,
                y_1 = y_1,
                Sigma_y_k_0 = Sigma_y_k_0,
                Sigma_y_k_1 = Sigma_y_k_1,
                quantiles = quantiles,
                K = K, Nq = N
    ))
}
