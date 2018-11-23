data {
  int<lower=1> N; // number of quantiles
  int<lower=1> K; // number of sites
  //means:
  vector[N] y_0[K]; // quantile of control group
  vector[N] y_1[K]; // quantile of treatment group
  //variances:
  cov_matrix[N] Sigma_y_k_0[K];
  cov_matrix[N] Sigma_y_k_1[K];
  //priors:
  cov_matrix[N] prior_dispersion_on_beta_0; // prior dispersion on beta_0
  cov_matrix[N] prior_dispersion_on_beta_1; // prior dispersion on beta_1

  int pooling_type; //0 if none, 1 if partial, 2 if full
}
transformed data {
  int K_pooled; // number of modelled sites if we take into account pooling
  // int N_pooled;
  if(pooling_type == 2){
    K_pooled = 0;
    // N_pooled = 0;
  }
  if(pooling_type != 2){
    K_pooled = K;
    // N_pooled = N;
  }
}
parameters {
  positive_ordered[N] beta_0; // true (level 2) DIFFERENCE in treatment and control quantiles
  positive_ordered[N] beta_0_k[K_pooled]; // intermediate treatment EFFECT

  //means:
  // positive_ordered[N] treatment_quantiles;
  positive_ordered[N] treatment_quantiles;
  positive_ordered[N] treatment_quantiles_k[K_pooled];

  //variances:
  cholesky_factor_corr[N] L_Omega_0; // cholesky factor of correlation matrix
  cholesky_factor_corr[N] L_Omega_1; // cholesky factor of correlation matrix
  vector<lower=0>[N] tau_0; // scale
  vector<lower=0>[N] tau_1; // scale
}
transformed parameters {
  cov_matrix[N] Sigma_0;
  cov_matrix[N] Sigma_1;
  // positive_ordered[N] beta_1;
  vector[N] beta_1;
  vector[N] beta_1_k[K_pooled];
  beta_1 =  treatment_quantiles - beta_0;
  if(K_pooled > 0)
    for(k in 1:K_pooled)
      beta_1_k[k] =  treatment_quantiles_k[k] - beta_0_k[k];
  Sigma_0 = diag_pre_multiply(tau_0,L_Omega_0) * diag_pre_multiply(tau_0,L_Omega_0)' ;
  Sigma_1 = diag_pre_multiply(tau_1,L_Omega_1) * diag_pre_multiply(tau_1,L_Omega_1)' ;
}

model {
  //mean priors:

  //variance priors:
  tau_0 ~ cauchy(0,20);
  tau_1 ~ cauchy(0,20);
  L_Omega_0 ~ lkj_corr_cholesky(1);
  L_Omega_1 ~ lkj_corr_cholesky(1);

  //distribution of betas (level 2 and level 1):
  if(pooling_type == 0) {
    // beta_0 and beta_1 allowed to wander around (but not too much),
    // we will discard them afterwards
    beta_0 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_0);
    beta_1 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_1);
    for(k in 1:K_pooled){ //no level 2 there
      beta_0_k[k] ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_0);
      beta_1_k[k] ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_1);
    }
  } else {
    beta_0 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_0);
    beta_1 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_1);
    if(K_pooled > 0) { //if pooling_type == 2 K_pooled is 0 anyway
      for(k in 1:K_pooled){
        beta_0_k[k] ~ multi_normal(beta_0, Sigma_0);
        beta_1_k[k] ~ multi_normal(beta_1, Sigma_1);
      }
    }
  }

  // model for observed means:
  if(pooling_type != 2) {
    for(k in 1:K){
      y_0[k] ~ multi_normal(beta_0_k[k], Sigma_y_k_0[k]);
      y_1[k] ~ multi_normal(treatment_quantiles_k[k], Sigma_y_k_1[k]);
    }
  } else {
   for(k in 1:K){
      y_0[k] ~ multi_normal(beta_0, Sigma_y_k_0[k]);
      y_1[k] ~ multi_normal(treatment_quantiles, Sigma_y_k_1[k]);
    }
  }
}
