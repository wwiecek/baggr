data {
  int<lower=1> N; // number of quantiles
  int<lower=1> K; // number of sites
  //means:
  vector[N] y_0[K]; // quantile of control group
  vector[N] y_1[K]; // quantile of treatment effect
                    // (treatment group quantiles = y_0 + y_1)
  //variances:
  cov_matrix[N] Sigma_y_k_0[K];
  cov_matrix[N] Sigma_y_k_1[K];
  //priors:
  cov_matrix[N] prior_dispersion_on_beta_0; // prior dispersion on beta_0
  cov_matrix[N] prior_dispersion_on_beta_1; // prior dispersion on beta_1
}
parameters {
  //means:
  positive_ordered[N] beta_0;
  positive_ordered[N] beta_0_k[K];
  vector[N] beta_1; // true parent DIFFERENCE in treatment and control quantiles
  vector[N] beta_1_k[K]; // intermediate treatment EFFECT
  //variances:
  cholesky_factor_corr[N] L_Omega_0; // cholesky factor of correlation matrix
  cholesky_factor_corr[N] L_Omega_1; // cholesky factor of correlation matrix
  vector<lower=0>[N] tau_0; // scale
  vector<lower=0>[N] tau_1; // scale
}
transformed parameters {
  // positive_ordered[N]  treatment_quantiles[K]; // disabled for 18 Oct 2018
  positive_ordered[N]  treatment_quantiles_k[K]; // true treatment quantiles
  cov_matrix[N] Sigma_0;
  cov_matrix[N] Sigma_1;

  for(k in 1:K)
   treatment_quantiles_k[k] = beta_1_k[k] + beta_0_k[k];
  Sigma_0 = diag_pre_multiply(tau_0,L_Omega_0) * diag_pre_multiply(tau_0,L_Omega_0)' ;
  Sigma_1 = diag_pre_multiply(tau_1,L_Omega_1) * diag_pre_multiply(tau_1,L_Omega_1)' ;
}
model {
  //mean priors:
  beta_0 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_0);
  beta_1 ~ multi_normal(rep_vector(0, N), prior_dispersion_on_beta_1);
  //variance priors:
  tau_0 ~ cauchy(0,20);
  tau_1 ~ cauchy(0,20);
  L_Omega_0 ~ lkj_corr_cholesky(1);
  L_Omega_1 ~ lkj_corr_cholesky(1);
  //data model:
  for(k in 1:K){
    beta_0_k[k] ~ multi_normal(beta_0, Sigma_0);
    beta_1_k[k] ~ multi_normal(beta_1, Sigma_1);
    y_0[k] ~ multi_normal(beta_0_k[k], Sigma_y_k_0[k]);
    y_1[k] ~ multi_normal(beta_1_k[k], Sigma_y_k_1[k]);
  }
}
// generated quantities{
//   vector[N] posterior_predictive_beta_0_k;
//   vector[N] posterior_predictive_treatment_quantiles;
//   vector[N] posterior_predictive_beta_1_k;
//   posterior_predictive_beta_0_k = multi_normal_rng(beta_0, Sigma_0);
//   posterior_predictive_beta_1_k = multi_normal_rng(beta_1, Sigma_1);
//   posterior_predictive_treatment_quantiles = posterior_predictive_beta_0_k + posterior_predictive_beta_1_k;
// }
