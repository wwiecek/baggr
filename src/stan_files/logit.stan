functions {
#include /functions/prior_increment.stan
}

data {
  int<lower=0> K;  // number of sites
  int<lower=0> N;  // total number of observations
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[N,Nc] X;  //covariate values (design matrix for FE)
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int<lower=0,upper=1> y[N];
  int<lower=0,upper=K> site[N];
  vector<lower=0,upper=1>[N] treatment;

  //priors
  int prior_hypermean_fam;
  int prior_hypersd_fam;
  vector[3] prior_hypermean_val;
  vector[3] prior_hypersd_val;

  //cross-validation variables:
  int<lower=0> N_test;
  int<lower=0> K_test;
  real test_y[N_test];
  int<lower=0, upper=K> test_site[N_test];
  int<lower=0, upper=1> test_treatment[N_test];
}
transformed data {
  int K_pooled; // number of modelled sites if we take pooling into account
  if(pooling_type == 2)
    K_pooled = 0;
  if(pooling_type != 2)
    K_pooled = K;
}
parameters {
  real baseline[K];
  real mu[pooling_type != 0? 1: 0];
  real<lower=0> tau[pooling_type == 1? 1: 0];
  real eta[K_pooled];
  vector[Nc] beta;
}
transformed parameters {
  real theta_k[K_pooled];
  for(k in 1:K_pooled){
    if(pooling_type == 0)
      theta_k[k] = eta[k];
    if(pooling_type == 1)
      theta_k[k] = mu[1] + eta[k]*tau[1];
  }
}
model {
  vector[N] fe;
  if(N > 0){
    if(Nc == 0)
      fe = rep_vector(0.0, N);
    else
      fe = X*beta;
  }

  baseline ~ normal(0, 10);

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

  if(pooling_type == 1)
    eta ~ normal(0,1);

  for(i in 1:N){
    if(pooling_type < 2)
      y[i] ~ bernoulli_logit(baseline[site[i]] + theta_k[site[i]] * treatment[i] + fe[i]);
    if(pooling_type == 2)
      y[i] ~ bernoulli_logit(baseline[site[i]] + mu[1] * treatment[i] + fe[i]);
  }
}

/*
generated_quantities {
  real logpd[K_test > 0? 1: 0];
  // vector[K_test] fe_k_test;
  if(K_test > 0){
    // if(Nc == 0)
      // fe_k_test = rep_vector(0.0, K_test);
    // else
      // fe_k_test = X_test*beta;
    logpd[1] = 0;
    for(k in 1:K_test){
      if(pooling_type == 1)
        logpd[1] += bernoulli_logit_lpmf(test_theta_hat_k[k] | baseline[site[i]] + mu[1] * treatment[i] + fe[i]);
      // if(pooling_type == 2)
        // logpd[1] += bernoulli_logit_lpmf(test_theta_hat_k[k] | baseline[site[i]] + mu[1] * treatment[i] + fe[i]);
    }
  }
}
*/
