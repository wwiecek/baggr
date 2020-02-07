functions {
#include /functions/prior_increment.stan
}

data {
  //controls
  int pooling_type; //0 if none, 1 if partial, 2 if full

  //data
  int<lower=0> K; // number of groups
  vector[K] theta_hat_k;
  vector<lower=0>[K] se_theta_k;
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[K,Nc] X;  //covariate values (design matrix for FE)

  //priors
  int prior_hypermean_fam;
  int prior_hypersd_fam;
  vector[3] prior_hypermean_val;
  vector[3] prior_hypersd_val;

  //test data (cross-validation)
  int<lower=0> K_test;
  vector[K_test] test_theta_hat_k;
  vector<lower=0>[K_test] test_se_theta_k;
  matrix[K_test,Nc] X_test;  //covariate values (design matrix for FE)
}

transformed data {
  int K_pooled; // number of modelled sites if we take pooling into account
  if(pooling_type == 2)
    K_pooled = 0;
  if(pooling_type != 2)
    K_pooled = K;
}

parameters {
  real mu[pooling_type != 0? 1: 0];
  real<lower=0> tau[pooling_type == 1? 1: 0];
  vector[K_pooled] eta;
  vector[Nc] beta;
}
transformed parameters {
  vector[K_pooled] theta_k;
  for(k in 1:K_pooled){
    if(pooling_type == 0)
      theta_k[k] = eta[k];
    if(pooling_type == 1)
      theta_k[k] = mu[1] + eta[k]*tau[1];
  }
}
model {
  vector[K] fe_k;
  if(K > 0){
    if(Nc == 0)
      fe_k = rep_vector(0.0, K);
    else
      fe_k = X*beta;
  }

  //hypermean priors:
  if(pooling_type > 0)
    target += prior_increment_vec(prior_hypermean_fam, mu[1], prior_hypermean_val);
  else{
    for(k in 1:K)
      target += prior_increment_vec(prior_hypermean_fam, eta[k], prior_hypermean_val);
  }

  //hyper-SD priors:
  if(pooling_type == 1)
    target += prior_increment_vec(prior_hypersd_fam, tau[1], prior_hypersd_val);

  //fixed effect coefficients
  beta ~ normal(0, 10);

  //likelihood (block evaluated only if there are data, i.e. K>0)
  if(K > 0) {
    if(pooling_type == 1)
        eta ~ normal(0,1);
    if(pooling_type != 2)
        theta_hat_k ~ normal(theta_k + fe_k, se_theta_k);
    if(pooling_type == 2)
        theta_hat_k ~ normal(rep_vector(mu[1], K) + fe_k, se_theta_k);
  }
}

generated quantities {
  real logpd[K_test > 0? 1: 0];
  vector[K_test] fe_k_test;
  if(K_test > 0){
    if(Nc == 0)
      fe_k_test = rep_vector(0.0, K_test);
    else
      fe_k_test = X_test*beta;
    logpd[1] = 0;
    for(k in 1:K_test){
      if(pooling_type == 1)
        logpd[1] += normal_lpdf(test_theta_hat_k[k] | mu[1] + fe_k_test, sqrt(tau[1]^2 + test_se_theta_k[k]^2));
      if(pooling_type == 2)
        logpd[1] += normal_lpdf(test_theta_hat_k[k] | mu[1] + fe_k_test, sqrt(test_se_theta_k[k]^2));
    }
  }
}
