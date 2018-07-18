data {
  int<lower=0> K; // number of sites
  // int<lower=0> K_pooled; // number of sites if we take into account pooling
  real tau_hat_k[K]; // estimated treatment effects
  real<lower=0> se_tau_k[K]; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full
  real prior_upper_sigma_tau;
  real prior_tau_mean;
  real prior_tau_scale;
}
transformed data {
  int K_pooled; // number of modelled sites if we take into account pooling
  if(pooling_type == 2)
    K_pooled = 0;
  if(pooling_type != 2)
    K_pooled = K;
}
parameters {
  real tau;
  real<lower=0> sigma_tau;
  real tau_k[K_pooled];
}
model {
  sigma_tau ~ uniform(0, prior_upper_sigma_tau);

    if(pooling_type == 0){
      // we don't care about the other parameter, let it wander!
      // but in that case we need to give a sensible prior for individual tau's and their SE's?
      tau ~ normal(0, 1); //wander tau but avoid divergent transitions
      tau_k ~ normal(prior_tau_mean, prior_tau_scale); //tau_k 'take on' tau's prior
      tau_hat_k ~ normal(tau_k, se_tau_k);
    }
    if(pooling_type == 1){
      tau ~ normal(prior_tau_mean, prior_tau_scale);//
      tau_k ~ normal(tau, sigma_tau);
      tau_hat_k ~ normal(tau_k, se_tau_k);
    }
    if(pooling_type == 2){
      tau ~ normal(prior_tau_mean, prior_tau_scale);
      tau_hat_k ~ normal(tau, se_tau_k);
    }
}

generated quantities{
  // real predicted_tau_k;
  // predicted_tau_k = normal_rng(tau, sigma_tau);
}
