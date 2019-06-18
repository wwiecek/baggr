// parameterization using normal(0,1)
// for Rubin's model

data {
  int<lower=0> K; // number of sites
  // int<lower=0> K_pooled; // number of sites if we take into account pooling
  real tau_hat_k[K]; // estimated treatment effects
  real<lower=0> se_tau_k[K]; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full
  real prior_upper_sigma_tau;
  real prior_tau_mean;
  real prior_tau_scale;

  //cross-validation variables:
  int<lower=0> K_test; // number of sites
  real test_tau_hat_k[K_test]; // estimated treatment effects
  real<lower=0> test_se_k[K_test]; // s.e. of effect estimates

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
  real eta[K_pooled];
}
transformed parameters {
  real tau_k[K_pooled];
  for(k in 1:K_pooled){
    if(pooling_type == 0)
      tau_k[k] = eta[k];
    if(pooling_type == 1)
      tau_k[k] = tau + eta[k]*sigma_tau;
  }
}
model {
  sigma_tau ~ uniform(0, prior_upper_sigma_tau);
  if(pooling_type == 0){
    tau ~ normal(0, 1); //wander tau but avoid divergent transitions
    eta ~ normal(prior_tau_mean, prior_tau_scale);
    tau_hat_k ~ normal(tau_k, se_tau_k);
  }
  if(pooling_type == 1){
    tau ~ normal(prior_tau_mean, prior_tau_scale);
    eta ~ normal(0,1);
    tau_hat_k ~ normal(tau_k, se_tau_k);
  }
  if(pooling_type == 2){
    tau ~ normal(prior_tau_mean, prior_tau_scale);
    tau_hat_k ~ normal(tau, se_tau_k);
  }
}

generated quantities {
  real logpd = 0;
  if(K_test > 0){
    for(k in 1:K_test){
      if(pooling_type == 1)
        logpd += normal_lpdf(test_tau_hat_k[k] | tau, sqrt(sigma_tau^2 + test_se_k[k]^2));
      if(pooling_type == 2)
        logpd += normal_lpdf(test_tau_hat_k[k] | tau, sqrt(test_se_k[k]^2));
    }
  }
}
