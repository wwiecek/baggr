functions {
#include /functions/prior_increment.stan
}

data {
  int<lower=0> K;  // number of sites
  int<lower=0> N;  // total number of observations
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[N,Nc] X;  //covariate values (design matrix for FE)
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int pooling_baseline; //pooling for proportions in control arm;
                        //0 if none, 1 if partial
  int<lower=0,upper=1> y[N];
  int<lower=0,upper=K> site[N];
  vector<lower=0,upper=1>[N] treatment;

  //priors for baseline parameters
  int prior_hbasemean_fam;
  int prior_hbasesd_fam;
  vector[3] prior_hbasemean_val;
  vector[3] prior_hbasesd_val;
  //priors for effects::
  int prior_hypermean_fam;
  int prior_hypersd_fam;
  int prior_beta_fam;
  vector[3] prior_hypermean_val;
  vector[3] prior_hypersd_val;
  vector[3] prior_beta_val;

  //cross-validation variables:
  int<lower=0> N_test;
  int<lower=0> K_test;
  matrix[N_test, Nc] X_test;
  int<lower=0,upper=1> test_y[N_test];
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
  // vector[K] baseline_k;
  real mu_baseline[pooling_baseline != 0? 1: 0];
  real<lower=0> tau_baseline[pooling_baseline != 0? 1: 0];
  real mu[pooling_type != 0? 1: 0];
  real<lower=0> tau[pooling_type == 1? 1: 0];
  vector[K_pooled] eta;
  vector[K] eta_baseline;
  vector[Nc] beta;
}
transformed parameters {
  vector[K_pooled] theta_k;
  vector[K] baseline_k;

  if(pooling_type == 0)
    theta_k = eta;
  else if(pooling_type == 1)
    theta_k = rep_vector(mu[1], K_pooled) + tau[1]*eta;

  if(pooling_baseline == 0)
    baseline_k = eta_baseline;
  else if(pooling_baseline == 1)
    baseline_k = rep_vector(mu_baseline[1], K) + tau_baseline[1]*eta_baseline;
}
model {

  //controls/baselines (hyper)priors
  if(pooling_baseline == 0)
    target += prior_increment_vec(prior_hbasemean_fam, eta_baseline, prior_hbasemean_val);
  if(pooling_baseline == 1){
    eta_baseline ~ normal(0,1);
    target += prior_increment_real(prior_hbasemean_fam, mu_baseline[1], prior_hbasemean_val);
    target += prior_increment_real(prior_hbasesd_fam, tau_baseline[1], prior_hbasesd_val);
  }

  //hypermean priors:
  if(pooling_type > 0)
    target += prior_increment_real(prior_hypermean_fam, mu[1], prior_hypermean_val);
  else{
    for(k in 1:K)
      target += prior_increment_real(prior_hypermean_fam, eta[k], prior_hypermean_val);
  }

  //hyper-SD priors:
  if(pooling_type == 1)
    target += prior_increment_real(prior_hypersd_fam, tau[1], prior_hypersd_val);

  //fixed effect coefficients
  target += prior_increment_vec(prior_beta_fam, beta, prior_beta_val);

  if(pooling_type == 1)
    eta ~ normal(0,1);


  // Branching logic to account for fixed-effects yes/no, full pooling yes/no
  if(pooling_type < 2){
    if(Nc == 0){
      for(i in 1:N)
        y[i] ~ bernoulli_logit(baseline_k[site[i]] + theta_k[site[i]] * treatment[i]);
    }else{
      for(i in 1:N)
        y[i] ~ bernoulli_logit(baseline_k[site[i]] + theta_k[site[i]] * treatment[i] + X[i,]*beta);
    }
  }
  if(pooling_type == 2) {
    if(Nc == 0){
      for(i in 1:N)
        y[i] ~ bernoulli_logit(baseline_k[site[i]] + mu[1] * treatment[i]);
    }else{
      for(i in 1:N)
        y[i] ~ bernoulli_logit(baseline_k[site[i]] + mu[1] * treatment[i] + X[i,]*beta);
    }
  }
}


generated quantities {
  real logpd[K_test > 0? 1: 0];
  real theta_k_test[K_test];
  if(K_test > 0){
    logpd[1] = 0;
    if(pooling_type != 0){
      for(k in 1:K_test)
      theta_k_test[k] = normal_rng(mu[1], tau[1]);
      // This will only work if we predict for baselines which are already estimated
      if(Nc == 0){
        for(i in 1:N_test)
          logpd[1] += bernoulli_logit_lpmf(test_y[i] | baseline_k[test_site[i]] +
                                           theta_k_test[test_site[i]] * test_treatment[i]);
      } else {
        for(i in 1:N_test)
          logpd[1] += bernoulli_logit_lpmf(test_y[i] | baseline_k[test_site[i]] +
                                           theta_k_test[test_site[i]] * test_treatment[i] + X_test*beta);
      }
    }
  }
}

