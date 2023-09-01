functions {
#include /functions/prior_increment.stan
}

data {
  int<lower=0> K; // number of sites
  int<lower=2> P; // number of parameters (1 or 2)
  array[P,K] real theta_hat_k; // estimated treatment effects
  array[P,K] real<lower=0> se_theta_k; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int<lower=0, upper=1> cumsum; //if 1, tau = control + trt effect,
                                //if 0, tau (TE) does not get added to mu ("control")

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
  array[P,K_test] real test_theta_hat_k; // estimated treatment effects
  array[P,K_test] real<lower=0> test_se_theta_k; // s.e. of effect estimates

}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
}

parameters {
  array[pooling_type != 0] vector[P] mu;
  array[pooling_type == 1] cholesky_factor_corr[P] L_Omega;
  array[pooling_type == 1] vector<lower=0>[P] hypersd;    //  scale
  array[pooling_type != 2] matrix[P,K] eta;
}

transformed parameters {
  array[pooling_type != 2] matrix[P,K] theta_k;
  array[pooling_type == 1] matrix[P,P] tau;

  if(pooling_type == 0)
    theta_k[1] = eta[1];
  if(pooling_type == 1){
    // We use Cholesky decomposition per Stan manual recommendation
    // diag_pre_multiply is diag(hypersd) %*% L_Omega in R code
    // L %*% t(L) is correlation matrix C,
    // diag(hypersd) %*% C %*% diag(hypersd) is var-covariance
    // it follows that pre-multiplying gives a Cholesky factor of variance-cov.
    tau[1] = diag_pre_multiply(hypersd[1], L_Omega[1]);
    if(cumsum==1)
      theta_k[1] = rep_matrix(mu[1], K) + tau[1] * eta[1];
    else
      theta_k[1] = rep_matrix(cumulative_sum(mu[1]), K) + tau[1] * eta[1];
  }
}

model {

  // priors: hypermean
  // if(pooling_type == 0 && K > 0)
    // for (k in 1:K)
      // theta_k[k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  if(pooling_type != 0) {
    if(prior_hypermean_fam == 3)
      mu[1] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  } else {
    for(k in 1:K)
      eta[1][,k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  }

  //priors variance/correlation
  if(pooling_type == 1) {
    to_vector(eta[1]) ~ std_normal();

    for(p in 1:P)
      hypersd[1][p] ~ realprior(prior_hypersd_fam, prior_hypersd_val);

    //for Omega only LKJ allowed for now
    L_Omega[1] ~ lkj_corr_cholesky(prior_hypercor_val);
  }

  if(pooling_type == 1 && K > 0) {
    for (k in 1:K) {
      // theta_k[k] ~ multi_normal(mu[1], tau[1]);
      for(p in 1:P)
        theta_hat_k[p,k] ~ normal(theta_k[1][p,k], se_theta_k[p,k]);
    }
  }
  if(pooling_type == 2 && K > 0)
    for (k in 1:K)
      for(p in 1:P)
        theta_hat_k[p,k] ~ normal(mu[1][p], se_theta_k[p,k]);
}

generated quantities {
  array[K_test > 0] real logpd;
  if(K_test > 0) {
    logpd[1] = 0;
    for(k in 1:K_test){
      for(p in 1:P) {
        if(pooling_type == 1)
        logpd[1] += normal_lpdf(test_theta_hat_k[p,k] | mu[1],
                                sqrt(tau[1][p,p]^2 + test_se_theta_k[p,k]^2));
        if(pooling_type == 2)
        logpd[1] += normal_lpdf(test_theta_hat_k[p,k] | mu[1],
                                sqrt(test_se_theta_k[p,k]^2));
      }
    }
  }
}
