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
  array[N] int<lower=0,upper=K> site;
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
  //priors for regression model:
  int prior_sigma_fam;
  vector[3] prior_sigma_val;

  //cross-validation variables:
  int<lower=0> N_test;
  int<lower=0> K_test;
  matrix[N_test, Nc] X_test;
  array[N_test] int<lower=0, upper=K> test_site;
  array[N_test] int<lower=0, upper=1> test_treatment;

  // NORMAL specific:
  array[N] real y;
  array[N_test] real test_y;
  array[K_test] real test_sigma_y_k;

}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
  int K_bsl_pooled = pooling_baseline == 2 ? 0 : K;
}
parameters {
  // SHARED ACROSS FULL MODELS:
  array[pooling_baseline == 1] real mu_baseline;
  array[pooling_type != 0] real mu;
  array[pooling_baseline == 1] real<lower=0> tau_baseline;
  array[pooling_type == 1] real<lower=0> tau;
  vector[K_pooled] eta;
  vector[K_bsl_pooled] eta_baseline;
  vector[Nc] beta;

  // NORMAL specific:
  vector<lower=0>[K] sigma_y_k;
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
  else
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
  if(pooling_baseline == 0)
    eta_baseline ~ vecprior(prior_control_fam, prior_control_val);
  else if(pooling_baseline == 1){
    eta_baseline ~ normal(0,1);
    mu_baseline[1]  ~ realprior(prior_control_fam, prior_control_val);
    tau_baseline[1] ~ realprior(prior_control_sd_fam, prior_control_sd_val);
  }

  //priors for trt effect and normal likelihood
  if(pooling_type == 0){
    eta ~ vecprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ normal(baseline_k[site] + theta_k[site] .* treatment + fe, sigma_y_k[site]);
  } else if(pooling_type == 1){
    eta ~ normal(0,1);
    tau[1] ~ realprior(prior_hypersd_fam, prior_hypersd_val);
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ normal(baseline_k[site] + theta_k[site] .* treatment + fe, sigma_y_k[site]);
  } else {
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    y ~ normal(baseline_k[site] + mu[1] * treatment + fe, sigma_y_k[site]);
  }
  // (normal model only:) error term priors
  sigma_y_k ~ vecprior(prior_sigma_fam, prior_sigma_val);
}
generated quantities {
  // to do this, we must first (outside of Stan) calculate SEs in each test group,
  // i.e. test_sigma_y_k
  array[K_test > 0] real logpd;
  vector[N_test] fe_test;
  if(K_test > 0){
    if(Nc == 0)
      fe_test = rep_vector(0.0, N_test);
    else
      fe_test = X_test*beta;
    logpd[1] = 0;
    for(i in 1:N_test){
      if(pooling_type == 1)
        logpd[1] += normal_lpdf(test_y[i] | baseline_k[test_site[i]] + mu[1] * test_treatment[i] + fe_test[i],
                                sqrt(tau[1]^2 + test_sigma_y_k[test_site[i]]^2));
      if(pooling_type == 2)
        logpd[1] += normal_lpdf(test_y[i] | baseline_k[test_site[i]] + mu[1] * test_treatment[i] + fe_test[i],
                                sqrt(test_sigma_y_k[test_site[i]]^2));
    }
  }
}
