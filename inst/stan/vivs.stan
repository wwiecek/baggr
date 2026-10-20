functions {
#include /functions/prior_increment.stan
}

data {
  int<lower=0> K; // number of sites
  int<lower=2> P; // number of parameters (1 or 2)
  array[P,K] real yi; // estimated treatment effects
  array[P,K] real<lower=0> sei; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int<lower=0, upper=1> cumsum; //if 1, use baseline and effect separately;
                                //if 0, add baseline to the effect mean

  // priors:
  int prior_hypermean_fam;
  vector[P] prior_hypermean_mean;
  matrix<lower=0>[P, P] prior_hypermean_scale;
  int prior_hypersd_fam;
  vector[3] prior_hypersd_val;
  int prior_hypercor_fam; //only LKJ allowed for now...
  real prior_hypercor_val;

  //cross-validation variables:
  int<lower=0> K_test; // number of sites
  array[P,K_test] real test_yi; // estimated treatment effects
  array[P,K_test] real<lower=0> test_sei; // s.e. of effect estimates

}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
}

parameters {
  array[pooling_type != 0] vector[P] theta;
  array[pooling_type == 1] cholesky_factor_corr[P] L_Omega;
  array[pooling_type == 1] vector<lower=0>[P] sigma;    //  scale
  array[pooling_type != 2] matrix[P,K] eta;
}

transformed parameters {
  array[pooling_type != 2] matrix[P,K] theta_k;
  array[pooling_type == 1] matrix[P,P] L_Sigma;

  if(pooling_type == 0)
    theta_k[1] = eta[1];
  if(pooling_type == 1){
    // We use Cholesky decomposition per Stan manual recommendation
    // diag_pre_multiply is diag(sigma) %*% L_Omega in R code
    // L %*% t(L) is correlation matrix C,
    // diag(sigma) %*% C %*% diag(sigma) is var-covariance
    // it follows that pre-multiplying gives a Cholesky factor of variance-cov.
    L_Sigma[1] = diag_pre_multiply(sigma[1], L_Omega[1]);
    if(cumsum==1)
      theta_k[1] = rep_matrix(theta[1], K) + L_Sigma[1] * eta[1];
    else
      theta_k[1] = rep_matrix(cumulative_sum(theta[1]), K) + L_Sigma[1] * eta[1];
  }
}

model {

  // priors: hypermean
  // if(pooling_type == 0 && K > 0)
    // for (k in 1:K)
      // theta_k[k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  if(pooling_type != 0) {
    if(prior_hypermean_fam == 3)
      theta[1] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  } else {
    for(k in 1:K)
      eta[1][,k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  }

  //priors variance/correlation
  if(pooling_type == 1) {
    to_vector(eta[1]) ~ std_normal();

    for(p in 1:P)
      sigma[1][p] ~ realprior(prior_hypersd_fam, prior_hypersd_val);

    //for Omega only LKJ allowed for now
    L_Omega[1] ~ lkj_corr_cholesky(prior_hypercor_val);
  }

  if(pooling_type == 1 && K > 0) {
    for (k in 1:K) {
      // theta_k[k] ~ multi_normal(theta[1], L_Sigma[1]);
      for(p in 1:P)
        yi[p,k] ~ normal(theta_k[1][p,k], sei[p,k]);
    }
  }
  if(pooling_type == 2 && K > 0)
    for (k in 1:K)
      for(p in 1:P)
        yi[p,k] ~ normal(theta[1][p], sei[p,k]);
}

generated quantities {
  array[K_test > 0] real logpd;

  if (K_test > 0) {
    logpd[1] = 0;

    for (k in 1:K_test) {
      vector[P] y = to_vector(test_yi[, k]);
      vector[P] se = to_vector(test_sei[, k]);

      if (pooling_type == 1) {
        matrix[P,P] Sigma = multiply_lower_tri_self_transpose(L_Sigma[1]);
        matrix[P,P] V = Sigma + diag_matrix(square(se));
        logpd[1] += multi_normal_lpdf(y | theta[1], V);
      }

      if (pooling_type == 2) {
        matrix[P,P] V = diag_matrix(square(se));
        logpd[1] += multi_normal_lpdf(y | theta[1], V);
      }
    }
  }
}

