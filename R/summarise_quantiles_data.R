#' summarise_quantiles_data
#' Given individual level data, return list of summary statistics
#' Of quantile means and Sigma's, as well as K, N
#' Already formatted as a valid input into a Stan model of quantiles
# summarise_quantiles_data(microcredit, c(.2, .4, .6), "study", "profit")
#' @importFrom quantreg rq
#'
summarise_quantiles_data <- function(data, quantiles,
                                     grouping  = "site",
                                     outcome   = "outcome",
                                     treatment = "treatment") {
  N <- length(quantiles)
  K <- length(unique(data[[grouping]]))

  # Calculate means and SE's of our quantiles via quantreg::qr()
  # Not very elegant & very slow.
  # for(group in unique(data[[grouping]])) {
  #   qr <- quantreg::rq(data[[outcome]][data[[grouping]] == group] ~ data[[treatment]][data[[grouping]] == group],
  #                tau = quantiles_list)
  #   y_0 <- qr$coef[1,]
  #   y_1 <- qr$coef[2,]
  # }
  i <- 0 #there's no shame in a loop
  groups <- unique(data[[grouping]])
  y_0 <- matrix(NA, length(groups), length(quantiles))
  y_1 <- matrix(NA, length(groups), length(quantiles))
  y_0_se <- matrix(NA, length(groups), length(quantiles))
  y_1_se <- matrix(NA, length(groups), length(quantiles))
  for(i in 1:length(groups)) {
    gr <- groups[i]
    a <- quantreg::rq(data[[outcome]][data[[grouping]] == gr] ~
                        data[[treatment]][data[[grouping]] == gr],
                      tau = quantiles)
    for(j in 1:length(quantiles)) {
      y_0[i,j] <- summary(a, se = "iid")[[j]]$coef[1,1]
      y_1[i,j] <- summary(a, se = "iid")[[j]]$coef[2,1]
      y_0_se[i,j] <- summary(a, se = "iid")[[j]]$coef[1,2]
      y_1_se[i,j] <- summary(a, se = "iid")[[j]]$coef[2,2]
    }
  }
  # Obtain Sigma's from data
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
