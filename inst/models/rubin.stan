  data {
    int<lower=0> K; // number of sites
    real tau_hat_k[K]; // estimated treatment effects
    real<lower=0> se_k[K]; // s.e. of effect estimates
    int pooling_type; //0 if none, 1 if partial, 2 if full
    real prior_max;
  }
  parameters {
    real tau;
    real<lower=0> sigma_tau;
    real tau_k[K];
  }
  transformed parameters {

  }
  model {
    if(pooling_type == 0)
      tau_hat_k ~ normal(tau_k, se_k); //third level normal

    sigma_tau ~ uniform(0,prior_max);//

    if(pooling_type == 1){
      tau ~ normal(0,1000);//
      tau_k ~ normal(tau, sigma_tau); // second level normal
      tau_hat_k ~ normal(tau_k, se_k); //third level normal
    }


    if(pooling_type == 2){
      tau ~ normal(0,1000);//
      tau_k ~ normal(tau, 0); // second level normal
      tau_hat_k ~ normal(tau_k, se_k); //third level normal
    }


  }
generated quantities{
  real predicted_tau_k;
  predicted_tau_k = normal_rng(tau, sigma_tau);
}
