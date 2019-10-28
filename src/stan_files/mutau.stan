data {
  int<lower=2> K; // number of sites
  int<lower=2> P; // number of parameters (1 or 2)
  real tau_hat_k[P,K]; // estimated treatment effects
  real<lower=0> se_tau_k[P,K]; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full

  // priors:
  int prior_hypermean_fam;
  vector[P] prior_hypermean_mean;
  matrix<lower=0>[P, P] prior_hypermean_scale;
  int prior_hypersd_fam;
  real prior_hypersd_val[2];
  int prior_hypercor_fam; //only LKJ allowed for now...
  real prior_hypercor_val[1];

  //cross-validation variables:
  int<lower=0> K_test; // number of sites
  real test_tau_hat_k[P,K_test]; // estimated treatment effects
  real<lower=0> test_se_k[P,K_test]; // s.e. of effect estimates

}

transformed data {
  int K_pooled; // number of modelled sites if we take into account pooling
  if(pooling_type == 2)
  K_pooled = 0;
  if(pooling_type != 2)
  K_pooled = K;
}

parameters {
  vector[P] tau[pooling_type != 0? 1: 0];
  vector[P] tau_k[K_pooled];
  corr_matrix[P] Omega[pooling_type == 1? 1: 0];        //  correlation
  vector<lower=0>[P] theta[pooling_type == 1? 1: 0];    //  scale
}

transformed parameters {
  matrix[P,P] sigma_tau[pooling_type == 1? 1: 0];
  if(pooling_type == 1)
  sigma_tau[1] = quad_form_diag(Omega[1],theta[1]);
}

model {

  // priors: hypermean
  if(pooling_type == 0)
    for (k in 1:K)
      tau_k[k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  if(pooling_type != 0) {
    if(prior_hypermean_fam == 3) //only LKJ allowed at the moment
      tau[1] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  }

  //priors variance/correlation
  if(pooling_type == 1) {
    if(prior_hypersd_fam == 0)
      theta[1] ~ uniform(prior_hypersd_val[1], prior_hypersd_val[2]);
    if(prior_hypersd_fam == 1)
      theta[1] ~ normal(prior_hypersd_val[1], prior_hypersd_val[2]);
    if(prior_hypersd_fam == 2)
      theta[1] ~ cauchy(prior_hypersd_val[1], prior_hypersd_val[2]);

    //for Omega only LKJ allowed for now
    Omega[1] ~ lkj_corr(prior_hypercor_val[1]);
  }

  if(pooling_type == 1) {
    for (k in 1:K) {
      tau_k[k] ~ multi_normal(tau[1], sigma_tau[1]);
      for(p in 1:P)
        tau_hat_k[p,k] ~ normal(tau_k[k,p], se_tau_k[p,k]);
    }
  }

  if(pooling_type == 2)
    for (k in 1:K)
      for(p in 1:P)
        tau_hat_k[p,k] ~ normal(tau[1][p], se_tau_k[p,k]);
}

generated quantities {
  real logpd = 0;
  if(K_test > 0)
  for(k in 1:K_test){
    for(p in 1:P) {
      //sigma_tau[p,p] is questionable!
      if(pooling_type == 1)
      logpd += normal_lpdf(test_tau_hat_k[p,k] | tau[1], sqrt(sigma_tau[1][p,p]^2 + test_se_k[p,k]^2));
      if(pooling_type == 2)
      logpd += normal_lpdf(test_tau_hat_k[p,k] | tau[1], sqrt(test_se_k[p,k]^2));
    }}
}
