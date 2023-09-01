functions {
#include /functions/prior_increment.stan
}

data {
  //controls
  int pooling_type; //0 if none, 1 if partial, 2 if full

  //data
  int<lower=0> K; // number of groups
  vector[K] theta_hat_k;
  vector<lower=0>[K] se_theta_k;
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[K,Nc] X;  //covariate values (design matrix for FE)

  //priors
  int prior_hypermean_fam;
  int prior_hypersd_fam;
  int prior_beta_fam;
  vector[3] prior_hypermean_val;
  vector[3] prior_hypersd_val;
  vector[3] prior_beta_val;

  //test data (cross-validation)
  int<lower=0> K_test;
  vector[K_test] test_theta_hat_k;
  vector<lower=0>[K_test] test_se_theta_k;
  matrix[K_test,Nc] X_test;  //covariate values (design matrix for FE)
}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
  vector[K_test] test_var_theta_k = test_se_theta_k .* test_se_theta_k;
}

parameters {
  array[pooling_type != 0] real mu;
  array[pooling_type == 1] real<lower=0> tau;
  vector[K_pooled] eta;
  vector[Nc] beta;
}
transformed parameters {
  /* if there is no pooling then eta's assume role of study means
     this is done to avoid defining yet another parameter but rather
     recycle something that already exists */
  vector[K_pooled] theta_k;
  if(pooling_type == 0)
    theta_k = eta;
  else if(pooling_type == 1)
    theta_k = mu[1] + eta*tau[1];
}
model {
  vector[K] fe_k;

  if(K > 0){
    if(Nc == 0){
      fe_k = rep_vector(0.0, K);
    }else{
      fe_k = X*beta;
      beta ~ vecprior(prior_beta_fam, prior_beta_val);
    }
  }
  if(pooling_type == 0){
    eta ~ vecprior(prior_hypermean_fam, prior_hypermean_val);
    theta_hat_k ~ normal(theta_k + fe_k, se_theta_k);
  }else if(pooling_type == 1){
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    tau[1] ~ realprior(prior_hypersd_fam, prior_hypersd_val);
    eta ~ normal(0,1);
    theta_hat_k ~ normal(theta_k + fe_k, se_theta_k);
  }else{
    mu[1] ~ realprior(prior_hypermean_fam, prior_hypermean_val);
    theta_hat_k ~ normal(mu[1] + fe_k, se_theta_k);
  }

}

generated quantities {
  array[K_test > 0] real logpd;
  vector[K_test] fe_k_test;
  if(K_test > 0){
    if(Nc == 0)
      fe_k_test = rep_vector(0.0, K_test);
    else
      fe_k_test = X_test*beta;
    logpd[1] = 0;
    for(k in 1:K_test){
      if(pooling_type == 1)
        logpd[1] += normal_lpdf(test_theta_hat_k[k] | mu[1] + fe_k_test, sqrt(tau[1]^2 + test_var_theta_k));
      if(pooling_type == 2)
        logpd[1] += normal_lpdf(test_theta_hat_k[k] | mu[1] + fe_k_test, test_se_theta_k);
    }
  }
}

