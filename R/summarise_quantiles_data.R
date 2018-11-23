#' summarise_quantiles_data
#' Given individual level data, return list of summary statistics
#' Of quantile means and Sigma's, as well as K, N
#' Already formatted as a valid input into a Stan model of quantiles
# summarise_quantiles_data(microcredit, c(.2, .4, .6), "study", "profit")
#' @importFrom quantreg rq
#' @export

summarise_quantiles_data <- function(data, quantiles,
                                     group  = "group",
                                     outcome   = "outcome",
                                     treatment = "treatment") {
  N <- length(quantiles)
  K <- length(unique(data[[group]]))
  if(!all(quantiles[-1] - quantiles[-length(quantiles)] > 0))
    stop("Vector of quantiles should be strictly monotonic")

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
  calc <- aggregate(outcome ~ group + treatment, function(x) {
    qr <- quantreg::rq(x ~ 1, tau = quantiles)
    qs <- summary(qr, se = "iid")
    unlist(lapply(qs, function(coef) {
      c("mean" = coef$coef[1,1],
        "se"   = coef$coef[1,2])
    }))
  }, data = data)

  # Proceed carefully, as the ordering of columns is crucial here!
  y_0    <- calc$outcome[calc$treatment == 0, 2*(1:length(quantiles)) - 1]
  y_0_se <- calc$outcome[calc$treatment == 0, 2*(1:length(quantiles))]
  y_1    <- calc$outcome[calc$treatment == 1, 2*(1:length(quantiles)) - 1]
  y_1_se <- calc$outcome[calc$treatment == 1, 2*(1:length(quantiles))]
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
  list(y_0 = y_0,
       y_1 = y_1,
       Sigma_y_k_0 = Sigma_y_k_0,
       Sigma_y_k_1 = Sigma_y_k_1,
       quantiles = quantiles,
       K = K, N = N
  )
}
