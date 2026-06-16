functions {
#include /functions/prior_increment.stan
}

data {
  // SHARED ACROSS FULL MODELS:
  int<lower=0> N;   // total number of observations
  int<lower=0> P;   // number of treatment columns
  int<lower=0> K;   // number of sites
  int<lower=0> Nc;  //number of covariates (fixed effects)
  int<lower=0> Ncluster;  //number of clusters
  matrix[N,Nc] X;   //covariate values (design matrix for FE)
  int clustered;    //indicator for cluster REs within trials
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int pooling_baseline; //pooling for proportions in control arm;
                        //0 if none, 1 if partial, else no bsl (==0)
  array[N] int<lower=0,upper=K> site;
  array[clustered == 1 ? N : 0] int<lower=0> cluster;
  matrix<lower=0,upper=1>[N,P] treatment;

  //PRIORS
  //baseline parameters:
  int prior_control_fam;
  int prior_control_sd_fam;
  vector[3] prior_control_val;
  vector[3] prior_control_sd_val;
  // effects:
  int prior_hypermean_fam[P];
  int prior_hypersd_fam[P];
  vector[3] prior_hypermean_val[P];
  vector[3] prior_hypersd_val[P];
  // coefficients in regression:
  int prior_beta_fam;
  vector[3] prior_beta_val;
  // clusters:
  int prior_cluster_fam;
  vector[3] prior_cluster_val;


  // cross-validation variables:
  int<lower=0> N_test;
  int<lower=0> K_test;
  matrix[N_test, Nc] X_test;
  array[N_test] int<lower=0, upper=K> test_site;
  matrix<lower=0,upper=1>[N_test,P] test_treatment;

  //LOGIT-specific:
  array[N] int<lower=0,upper=1> y;
  array[N_test] int<lower=0,upper=1> test_y;
}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
  int K_bsl_pooled = pooling_baseline == 2 ? 0 : K;
}
parameters {
  vector[pooling_type != 0 ? P : 0] mu;
  vector<lower=0>[pooling_type == 1 ? P : 0] tau;
  array[pooling_baseline == 1] real mu_baseline;
  array[pooling_baseline == 1] real<lower=0> tau_baseline;
  vector[K_pooled * P] eta;
  vector[K_bsl_pooled] eta_baseline;
  vector[Nc] beta;
  vector<lower=0>[clustered == 1 ? K : 0] sigma_cluster;
  vector[clustered == 1 ? Ncluster : 0] eta_cluster;
}

transformed parameters {
  matrix[K_pooled,P] theta_k;
  vector[K] baseline_k;
  matrix[K_pooled, P] eta_matrix;
  eta_matrix = to_matrix(eta, K_pooled, P);

  if(pooling_type == 0)
    theta_k = eta_matrix;
  else if(pooling_type == 1)
    theta_k = rep_matrix(mu', K_pooled) + eta_matrix .* rep_matrix(tau', K_pooled);

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

  //cluster REs
  vector[N] re = rep_vector(0, N);
  if(clustered) {
    eta_cluster ~ normal(0, 1);
    sigma_cluster ~ vecprior(prior_cluster_fam, prior_cluster_val);
    re = eta_cluster[cluster] .* sigma_cluster[site];
  }

  //controls/baselines (hyper)priors
  if(pooling_baseline == 0)
    eta_baseline ~ vecprior(prior_control_fam, prior_control_val);
  else if(pooling_baseline == 1){
    eta_baseline ~ normal(0,1);
    mu_baseline[1]  ~ realprior(prior_control_fam, prior_control_val);
    tau_baseline[1] ~ realprior(prior_control_sd_fam, prior_control_sd_val);
  }

  //priors for trt effect and Bernoulli-logit link likelihood
  if(pooling_type == 0){
    //each column of eta_matrix should use a different "P-specific" distribution
    for (p in 1:P)
      eta_matrix[:, p] ~ vecprior(prior_hypermean_fam[p], prior_hypermean_val[p]);
    y ~ bernoulli_logit(baseline_k[site] + rows_dot_product(theta_k[site,], treatment) + fe + re);
  } else if(pooling_type == 1) {
    eta ~ normal(0,1);
    for(p in 1:P){
      tau[p] ~ realprior(prior_hypersd_fam[p], prior_hypersd_val[p]);
      mu[p]  ~ realprior(prior_hypermean_fam[p], prior_hypermean_val[p]);
    }
    y   ~ bernoulli_logit(baseline_k[site] + rows_dot_product(theta_k[site,], treatment) + fe + re);
  } else {
    for(p in 1:P)
      mu[p] ~ realprior(prior_hypermean_fam[p], prior_hypermean_val[p]);
    y ~ bernoulli_logit(baseline_k[site] + treatment * mu + fe + re);
  }

}

generated quantities {
  array[K_test > 0] real logpd;
  if(K_test > 0) logpd[1] = 0;
  if (pooling_type != 0) {
    for (i in 1:N_test) {
      real mu_pred;

      if (pooling_type == 1) {
        matrix[K,P] theta_k_test;
        // Generate random effects for each site and treatment
        for (k in 1:K) {
          for (p in 1:P) {
            theta_k_test[k,p] = normal_rng(mu[p], tau[p]);
          }
        }
        mu_pred = baseline_k[test_site[i]]
                  + dot_product(theta_k_test[test_site[i],], test_treatment[i,])
                  + (Nc > 0 ? dot_product(X_test[i,], beta) : 0.0);
      } else {
        mu_pred = baseline_k[test_site[i]]
                  + dot_product(test_treatment[i,], mu)
                  + (Nc > 0 ? dot_product(X_test[i,], beta) : 0.0);
      }
      // Accumulate log predictive density
      logpd[1] += bernoulli_logit_lpmf(test_y[i] | mu_pred);
    }
  }
}

