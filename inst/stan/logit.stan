functions {
#include /functions/prior_increment.stan
}

data {
  // SHARED ACROSS FULL MODELS:
  int<lower=0> N;  // total number of observations
  int<lower=0> K;  // number of sites
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[N,Nc] X;  //covariate values (design matrix for FE)
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int pooling_baseline; //pooling for proportions in control arm;
                        //0 if none, 1 if partial, else no bsl (==0)
  int<lower=0,upper=K> site[N];
  vector<lower=0,upper=1>[N] treatment;

  //priors for baseline parameters
  int prior_control_fam;
  int prior_control_sd_fam;
  vector[3] prior_control_val;
  vector[3] prior_control_sd_val;
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
  int<lower=0, upper=K> test_site[N_test];
  vector<lower=0, upper=1>[N_test] test_treatment;

  //LOGIT-specific:
  int<lower=0,upper=1> y[N];
  int<lower=0,upper=1> test_y[N_test];
}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
  int K_bsl_pooled = pooling_baseline == 2 ? 0 : K;
}

parameters {
  // SHARED ACROSS FULL MODELS:
  real mu_baseline[pooling_baseline == 1];
  real mu[pooling_type != 0];
  real<lower=0> tau_baseline[pooling_baseline == 1];
  real<lower=0> tau[pooling_type == 1];
  vector[K_pooled] eta;
  vector[K_bsl_pooled] eta_baseline;
  vector[Nc] beta;
}

transformed parameters {
  vector[K_pooled] theta_k;
  vector[K] baseline_k;

  if(pooling_type == 0) {
    theta_k = eta;
  } else if(pooling_type == 1) {
    theta_k = mu[1] + tau[1] * eta;
  }

  if(pooling_baseline == 0) {
    baseline_k = eta_baseline;
  } else if(pooling_baseline == 1) {
    baseline_k = mu_baseline[1] + tau_baseline[1]*eta_baseline;
  } else
    baseline_k = rep_vector(0.0, K);
}

model {
  //fixed effects:
  vector[N] fe;
  if(N > 0){
    if(Nc == 0){
      fe = rep_vector(0.0, N);
    } else {
      fe = X*beta;
      beta ~ vecprior(prior_beta_fam, prior_beta_val);
    }
  }

  //controls/baselines (hyper)priors
  if(pooling_baseline == 0){
    eta_baseline ~ vecprior(prior_control_fam, prior_control_val);
  } else if(pooling_baseline == 1) {
    eta_baseline ~ normal(0,1);
    mu_baseline[1]  ~ realprior(prior_control_fam, prior_control_val);
    tau_baseline[1] ~ realprior(prior_control_sd_fam, prior_control_sd_val);
  }

  //priors for trt effect and Bernoulli likelihood
  if(pooling_type == 0){
    eta ~ vecprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ bernoulli_logit_lpmf(baseline_k[site] + theta_k[site] .* treatment + fe);
  } else if(pooling_type == 1){
    eta ~ normal(0,1);
    tau[1] ~ realprior(prior_hypersd_fam, prior_hypersd_val);
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ bernoulli_logit_lpmf(baseline_k[site] + theta_k[site] .* treatment + fe);
  } else {
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ bernoulli_logit_lpmf(baseline_k[site] + mu[1] * treatment + fe);
  }
}

generated quantities {
  real logpd[K_test > 0];
  vector[pooling_type == 1? K_test: 0] theta_k_test;
  vector[N_test] fe_test;
  if(K_test > 0){
    if(Nc == 0)
      fe_test = rep_vector(0.0, N_test);
    else
      fe_test = X_test*beta;
    logpd[1] = 0;
    if(pooling_type == 1){
      for(k in 1:K_test)
        theta_k_test[k] = normal_rng(mu[1], tau[1]);
      // This will only work if we predict for baselines which are already estimated
        logpd[1] += bernoulli_logit_lpmf(test_y | baseline_k[test_site] +
                                         theta_k_test[test_site] .* test_treatment + fe_test);
    }
    if(pooling_type == 2){
        logpd[1] += bernoulli_logit_lpmf(test_y | baseline_k[test_site] +
                                         mu[1] * test_treatment + fe_test);
    }
  }
}

