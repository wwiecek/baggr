
functions {
#include /functions/prior_increment.stan
}

data {
  // actual data inputs
  int M; // number of categories
  int N; // number of observations
  int P; // dimentionality of beta parameter
  int K; // number of sites
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int cat[N]; // category indicator
  int N_neg; // number of obvs in negative tail
  int N_pos; // number of obvs in positive tail
  vector[P] x[N]; // covariates
  int site[N]; // site indicator
  vector[N_neg] treatment_neg; // treatment in neg tail
  vector[N_pos] treatment_pos; // treatment in pos tail
  int site_neg[N_neg]; // site indicator in neg tail
  int site_pos[N_pos]; // site indicator in pos tail
  real y_neg[N_neg]; // the negative tail
  real y_pos[N_pos]; // the positive tail

  // now the prior inputs
  int prior_hypermu_fam;
  int prior_hypertau_fam;
  int prior_hypersigmacontrol_fam;
  int prior_hypersigmaTE_fam;
  int prior_hyperbeta_fam;
  int prior_hypersd_mu_fam;
  int prior_hypersd_tau_fam;
  int prior_hypersd_sigmacontrol_fam;
  int prior_hypersd_sigmaTE_fam;
  int prior_hypersd_beta_fam;
  vector[3] prior_hypermu_val;
  vector[3] prior_hypertau_val;
  vector[3] prior_hypersigmacontrol_val;
  vector[3] prior_hypersigmaTE_val;
  vector[3] prior_hyperbeta_val;
  vector[3] prior_hypersd_mu_val;
  vector[3] prior_hypersd_tau_val;
  vector[3] prior_hypersd_sigmacontrol_val;
  vector[3] prior_hypersd_sigmaTE_val;
  vector[3] prior_hypersd_beta_val;
}

transformed data {
  int K_pooled = (pooling_type == 2? 0: K); // number of modelled sites if we take pooling into account
}

parameters {
  real mu[pooling_type != 0? 2: 0];
  real tau[pooling_type != 0? 2: 0];
  real<lower=0> hypersd_mu[pooling_type == 1? 2: 0];
  real<lower=0> hypersd_tau[pooling_type == 1? 2: 0];
  real sigma_control[pooling_type != 0? 2: 0];
  real sigma_TE[pooling_type != 0? 2: 0];
  real<lower=0> hypersd_sigma_control[pooling_type == 1? 2: 0];
  real<lower=0> hypersd_sigma_TE[pooling_type == 1? 2: 0];
  matrix[K_pooled,2] eta_mu_k;
  matrix[K_pooled,2] eta_tau_k;
  matrix[K_pooled,2] eta_sigma_control_k;
  matrix[K_pooled,2] eta_sigma_TE_k;
  matrix[M-1,P] beta[pooling_type != 0? 1: 0]; // the parent parameters minus the Mth category
  matrix<lower=0>[M,P] hypersd_beta[pooling_type == 1? 1: 0]; // the set of M*P parent variances (not a covariance matrix)
  matrix[M,P] beta_k_raw[K_pooled]; // the hierarchical increments
}

transformed parameters{
  matrix[M,P] beta_k[K_pooled];
  matrix[K_pooled,2] mu_k;
  matrix[K_pooled,2] tau_k;
  matrix[K_pooled,2] sigma_control_k;
  matrix[K_pooled,2] sigma_TE_k;

  if(pooling_type == 0){
    mu_k = eta_mu_k;
    tau_k = eta_tau_k;
    sigma_control_k = eta_sigma_control_k;
    sigma_TE_k = eta_sigma_TE_k;
    beta_k = beta_k_raw;
  }

  if(pooling_type == 1){
    for (i in 1:2){
      mu_k[,i] = mu[i] + hypersd_mu[i]*eta_mu_k[,i];
      tau_k[,i] = tau[i] + hypersd_tau[i]*eta_tau_k[,i];
      sigma_control_k[,i] = sigma_control[i] + hypersd_sigma_control[i]*eta_sigma_control_k[,i];
      sigma_TE_k[,i] = sigma_TE[i] + hypersd_sigma_TE[i]*eta_sigma_TE_k[,i];
    }
    for (k in 1:K_pooled)
      //IS THIS BIT CORRECT? (mean Beta = 0 on the last row, but then group-effect is added)
      beta_k[k] = append_row(beta[1],rep_row_vector(0, P)) + hypersd_beta[1] .* beta_k_raw[k];
  }
}

model {

  // PRIORS
  if(pooling_type==0){
    for (m in 1:M)
      for (k in 1:K)
        target += prior_increment_vec(prior_hyperbeta_fam, beta_k_raw[k,m]', prior_hyperbeta_val);

    for (k in 1:K){ // should have the HYPERPARAMETER'S PRIORS WHEN YOU FIX PRIORS
      for (i in 1:2){
        target += prior_increment_real(prior_hypermu_fam, eta_mu_k[k,i], prior_hypermu_val);
        target += prior_increment_real(prior_hypertau_fam, eta_tau_k[k,i], prior_hypertau_val);
        target += prior_increment_real(prior_hypersigmacontrol_fam, eta_sigma_control_k[k,i], prior_hypersigmacontrol_val);
        target += prior_increment_real(prior_hypersigmaTE_fam, eta_sigma_TE_k[k,i], prior_hypersigmaTE_val);
      }
    }
  } // closes the pooling = 0 case

  if(pooling_type == 1){
    for (i in 1:2){
      target += prior_increment_real(prior_hypermu_fam, mu[i], prior_hypermu_val);
      target += prior_increment_real(prior_hypertau_fam, tau[i], prior_hypertau_val);
      target += prior_increment_real(prior_hypersd_mu_fam, hypersd_mu[i], prior_hypersd_mu_val);
      target += prior_increment_real(prior_hypersd_tau_fam, hypersd_tau[i], prior_hypersd_tau_val);
      target += prior_increment_real(prior_hypersigmacontrol_fam, sigma_control[i], prior_hypersigmacontrol_val);
      target += prior_increment_real(prior_hypersd_sigmacontrol_fam, hypersd_sigma_control[i], prior_hypersd_sigmacontrol_val);
      target += prior_increment_real(prior_hypersigmaTE_fam, sigma_TE[i], prior_hypersigmaTE_val);
      target += prior_increment_real(prior_hypersd_sigmaTE_fam, hypersd_sigma_TE[i], prior_hypersd_sigmaTE_val);
    } // closes the i loop
    target += prior_increment_vec(prior_hyperbeta_fam, to_vector(beta[1]) , prior_hyperbeta_val);
    // WW: HYPERSD_beta matrix here is converted to vector and then iid priors given on each element of this matrix
    //     is this ok?
    target += prior_increment_vec(prior_hypersd_beta_fam, to_vector(hypersd_beta[1]) , prior_hypersd_beta_val);
  } // closes the pooling = 1 case

  if(pooling_type ==2){
    target += prior_increment_vec(prior_hyperbeta_fam, to_vector(beta[1]) , prior_hyperbeta_val);
    for (i in 1:2){
      target += prior_increment_real(prior_hypermu_fam, mu[i], prior_hypermu_val);
      target += prior_increment_real(prior_hypertau_fam, tau[i], prior_hypertau_val);
      target += prior_increment_real(prior_hypersigmacontrol_fam, sigma_control[i], prior_hypersigmacontrol_val);
      target += prior_increment_real(prior_hypersigmaTE_fam, sigma_TE[i], prior_hypersigmaTE_val);
    } // closes the for loop indexed by i
  } // closes the pooling = 2 case



  // LIKELIHOOD
  if(N > 0){
    //Likelihood: 1/ hierarchy
    if(pooling_type==1){
      for (k in 1:K){
        for (m in 1:M){
          beta_k_raw[k,m] ~ normal(0,1);
        }
        eta_mu_k[k] ~ normal(0,1);
        eta_tau_k[k] ~ normal(0,1);
        eta_sigma_control_k[k] ~ normal(0,1);
        eta_sigma_TE_k[k] ~ normal(0,1);
      }
    } // closes the if pooling = 1 statement


    // Likelihood: 2/ categorical logit
    // All pooling types need the data level but split up as follows
    if(pooling_type < 2){
      for (n in 1:N)
        cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);
    } else if(pooling_type == 2){
      for (n in 1:N)
        cat[n] ~ categorical_logit(append_row(beta[1],rep_row_vector(0, P)) * x[n]);
    }

    //Likelihood: 3/ log-normal components
    if(pooling_type < 2) {
      y_pos ~ lognormal(mu_k[site_pos,2] + tau_k[site_pos,2].*treatment_pos,
                           exp(sigma_control_k[site_pos,2] + sigma_TE_k[site_pos,2].*treatment_pos)) ;
      y_neg ~ lognormal(mu_k[site_neg,1] + tau_k[site_neg,1].*treatment_neg,
                           exp(sigma_control_k[site_neg,1] + sigma_TE_k[site_neg,1].*treatment_neg));
    } else if(pooling_type == 2) {
      y_pos ~ lognormal(mu[2] + tau[2]*treatment_pos,
                        exp(sigma_control[2] + sigma_TE[2]*treatment_pos));
      y_neg ~ lognormal(mu[1] + tau[1]*treatment_neg,
                        exp(sigma_control[1] + sigma_TE[1]*treatment_neg));
    }
  }
}
