data {
  int<lower=2> K; // number of sites
  int<lower=2> P; // number of parameters (1 or 2)
  real tau_hat_k[P,K]; // estimated treatment effects
  real<lower=0> se_tau_k[P,K]; // s.e. of effect estimates
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int<lower=0, upper=1> joint; //is the distribution on parameters (mu and tau) joint?
  // real prior_upper_sigma_tau[P];
  vector[P] prior_tau_mean;
  vector[P] prior_upper_sigma_tau;
  matrix<lower=0>[P,P] prior_tau_scale;

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

  if(pooling_type != 0) { //hyperparam only if there's pooling
  if(joint == 1)
  tau[1] ~ multi_normal(prior_tau_mean, prior_tau_scale);
  if(joint == 0) {
    for(p in 1:P)
    tau[1][p] ~ normal(prior_tau_mean[p], prior_tau_scale[p,p]);
  }
  }
  if(pooling_type == 0) {
    //tau_k's 'take over' tau's prior distribution
    for (k in 1:K)
      tau_k[k] ~ multi_normal(prior_tau_mean, prior_tau_scale);
  }

  if(pooling_type == 1) {
    // parameter variance priors
    theta[1] ~ cauchy(0,10);
    Omega[1] ~ lkj_corr(3); // pushes towards independence
    for (k in 1:K) {
      tau_k[k] ~ multi_normal(tau[1], sigma_tau[1]);
      for(p in 1:P)
        tau_hat_k[p,k] ~ normal(tau_k[k,p], se_tau_k[p,k]);
    }
  }

  if(pooling_type == 2) {
    for (k in 1:K) {
      for(p in 1:P)
      tau_hat_k[p,k] ~ normal(tau[1][p], se_tau_k[p,k]);
    }
  }
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
