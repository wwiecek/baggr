data {
  int<lower=0> K;  // number of sites
  int<lower=0> N;  // total number of observations
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int<lower=0,upper=1> y[N];
  int<lower=0,upper=K> site[N];
  vector<lower=0,upper=1>[N] treatment;

  //priors (proof of concept)
  //0 = uniform, 1 = normal
  int prior_hypermean_fam;
  int prior_hypersd_fam;
  real prior_hypermean_val[2];
  real prior_hypersd_val[2];

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
  real mu[K];
  real tau[pooling_type != 0? 1: 0];
  real<lower=0> sigma_tau[pooling_type == 1? 1: 0];
  real eta[K_pooled];
}
transformed parameters {
  real tau_k[K_pooled];
  for(k in 1:K_pooled){
    if(pooling_type == 0)
      tau_k[k] = eta[k];
    if(pooling_type == 1)
      tau_k[k] = tau[1] + eta[k]*sigma_tau[1];
  }
}
model {

  mu ~ normal(0, 10);

  if(pooling_type > 0) {
    if(prior_hypermean_fam == 0)
    tau ~ uniform(prior_hypermean_val[1], prior_hypermean_val[2]);
    if(prior_hypermean_fam == 1)
    tau ~ normal(prior_hypermean_val[1], prior_hypermean_val[2]);
    if(prior_hypermean_fam == 2)
    tau ~ cauchy(prior_hypermean_val[1], prior_hypermean_val[2]);
  } else {
    if(prior_hypermean_fam == 0)
    eta ~ uniform(prior_hypermean_val[1], prior_hypermean_val[2]);
    if(prior_hypermean_fam == 1)
    eta ~ normal(prior_hypermean_val[1], prior_hypermean_val[2]);
    if(prior_hypermean_fam == 2)
    eta ~ cauchy(prior_hypermean_val[1], prior_hypermean_val[2]);
  }

  //hyper-SD priors:
  if(pooling_type == 1){
    if(prior_hypersd_fam == 0)
    target += uniform_lpdf(sigma_tau |
    prior_hypersd_val[1], prior_hypersd_val[2]);
    if(prior_hypersd_fam == 1)
    target += normal_lpdf(sigma_tau |
    prior_hypersd_val[1], prior_hypersd_val[2]);
    if(prior_hypersd_fam == 2)
    target += cauchy_lpdf(sigma_tau |
    prior_hypersd_val[1], prior_hypersd_val[2]);
  }

  if(pooling_type == 1)
    eta ~ normal(0,1);

  for(i in 1:N){
    if(pooling_type < 2)
      y[i] ~ bernoulli_logit(mu[site[i]] + tau_k[site[i]] * treatment[i]);
    if(pooling_type == 2)
      y[i] ~ bernoulli_logit(mu[site[i]] + tau[1] * treatment[i]);
  }
}
