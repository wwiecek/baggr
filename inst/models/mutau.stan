

data {
  int<lower=0> K;  // number of sites
  real mu_k_hat[K] ; // means
  real tau_k_hat[K] ; // slopes
  real se_mu_k[K] ; // ses for mus
  real se_tau_k[K] ; // ses for taus

}


parameters {
  real tau;//
  real mu;//
  vector[K] tau_k;//
  vector[K] mu_k;//
  real<lower=0> sigma_tau;//
  real<lower=0> sigma_mu;//

}
transformed parameters {

}

model {
     //  let me try with bounded uniform
  sigma_tau ~ uniform(0,100000);
  sigma_mu ~ uniform(0,100000);


  //  I am hard coding the priors on the hyperparameters here
  tau ~ normal(0,1000);
  mu ~ normal(0,1000); // one could later insert a more realistic mean but it doesnt matter

  tau_k ~ normal(tau, sigma_tau);

  mu_k ~ normal(mu, sigma_mu);

  for (k in 1:K) {

    mu_k_hat[k] ~ normal(mu_k[k],se_mu_k[k]);
    tau_k_hat[k] ~ normal(tau_k[k],se_tau_k[k]);
  }


}

generated quantities{
  real predicted_tau_k;
  real predicted_mu_k;
  predicted_tau_k = normal_rng(tau, sigma_tau);
  predicted_mu_k = normal_rng(mu, sigma_mu);
}
