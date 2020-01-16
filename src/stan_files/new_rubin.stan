functions {
#include /functions/prior_increment.stan
}

data {
  //controls
  int<lower=1> P; // number of parameters
  int pooling_type; //0 if none, 1 if partial, 2 if full

  //data
  int<lower=0> K; // number of groups
  vector[K] theta_hat_k[P];
  vector<lower=0>[K] se_theta_k[P];

  //priors
  int prior_hypermean_fam[P];
  int prior_hypersd_fam[P];
  vector[3] prior_hypermean_val[P];
  vector[3] prior_hypersd_val[P];
  // these are used only if P > 1
  int prior_hypercor_fam;
  real prior_hypercor_val[1];

  //test data (cross-validation)
  int<lower=0> K_test;
  matrix[K,P] test_theta_hat_k;
  matrix<lower=0>[K,P] test_se_theta_k[K];
}

transformed data {
  int K_pooled; // number of modelled sites if we take pooling into account
  if(pooling_type == 2)
    K_pooled = 0;
  if(pooling_type != 2)
    K_pooled = K;
}

parameters {
  vector[P] mu[pooling_type != 0? 1: 0];
  vector<lower=0>[P] sigma_tau[pooling_type == 1? 1: 0];
  matrix[P,K_pooled] eta;
  corr_matrix[P] Omega[pooling_type == 1? 1: 0];
}
transformed parameters {
  vector[K_pooled] theta_k[P];
  for(p in 1:P){
    for(k in 1:K_pooled){
      if(pooling_type == 0)
        theta_k[p,k] = eta[p,k];
      if(pooling_type == 1)
        theta_k[p,k] = mu[1,p] + eta[p,k]*sigma_tau[1,p];
    }
  }
}
model {
  for(p in 1:P){
    //hypermean priors:
    if(pooling_type > 0)
      target += prior_increment_vec(prior_hypermean_fam[p], mu[1,p], prior_hypermean_val[p]);
    else{
      for(k in 1:K)
        target += prior_increment_vec(prior_hypermean_fam[p], eta[p,k], prior_hypermean_val[p]);
    }
    //hyper-SD priors:
    if(pooling_type == 1)
      target += prior_increment_vec(prior_hypersd_fam[p], sigma_tau[1,p], prior_hypersd_val[p]);
  }

  //likelihood (block evaluated only if there are data, i.e. K>0)
  if(K > 0) {
    if(pooling_type == 1){
      for(k in 1:K)
        eta[,k] ~ multi_normal(rep_vector(0, P), Omega[1]);
    if(pooling_type != 2)
      for(p in 1:P)
        theta_hat_k[p] ~ normal(theta_k[p], se_theta_k[p]);
    }
    if(pooling_type == 2){
      for(p in 1:P)
        theta_hat_k[p] ~ normal(mu[1,p], se_theta_k[p]);
    }
  }
}
generated quantities {
}
