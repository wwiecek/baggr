functions {
#include /functions/prior_increment.stan
}

data {
  int<lower=0> N;  // total number of observations
  int<lower=0> K;  // number of sites
  int<lower=0> Nc; //number of covariates (fixed effects)
  matrix[N,Nc] X;  //covariate values (design matrix for FE)
  int pooling_type; //0 if none, 1 if partial, 2 if full
  int pooling_baseline; //pooling for proportions in control arm;
                        //0 if none, 1 if partial
  int joint_prior_mean;     //0 if no, 1 if yes
  int joint_prior_variance; //0 if no, 1 if yes
  int<lower=0,upper=K> site[N];
  vector<lower=0,upper=1>[N] treatment;
  int<lower=0> P;  //number of parameters -- CURRENTLY IT HAS TO BE 2, more will fail

  // priors:

  // joint prior for the mean
  int prior_hypermean_fam;
  vector[P] prior_hypermean_mean;
  matrix<lower=0>[P, P] prior_hypermean_scale;

  //priors for independent mean parameters (if joint_prior == 0)
  // int prior_control_fam;
  // vector[3] prior_control_val;
  // vector[3] prior_hypermean_val;

  //dispersions for effect and baselines:
  int prior_hypersd_fam;
  vector[3] prior_hypersd_val;
  int prior_control_sd_fam;
  vector[3] prior_control_sd_val;

  //correlations (only LKJ allowed for now)
  int prior_hypercor_fam;
  real prior_hypercor_val;

  //fixed effects (same dist. for all betas)
  int prior_beta_fam;
  vector[3] prior_beta_val;

    //cross-validation variables:
  int<lower=0> N_test;
  int<lower=0> K_test;
  matrix[N_test, Nc] X_test;
  int<lower=0, upper=K> test_site[N_test];
  int<lower=0, upper=1> test_treatment[N_test];

  // NORMAL specific:
  vector[N] y;
  vector[N_test] test_y;
  vector[K_test] test_sigma_y_k;
}

transformed data {
  int K_pooled = pooling_type == 2 ? 0 : K;
}

parameters {
  // corr_matrix[P] Omega;        // prior correlation
  cholesky_factor_corr[P] L_Omega[joint_prior_variance == 1 && pooling_type == 1];
  vector<lower=0>[P] hypersd[pooling_type == 1];        // prior scale
  vector[P] mu[pooling_type != 0];
  matrix[P,K] eta[pooling_type != 2];
  vector[Nc] beta;
  // NORMAL specific:
  vector<lower=0>[K] sigma_y_k;
}
transformed parameters {
  matrix[P,K] theta_k[pooling_type != 2];
  matrix[P,P] tau[pooling_type == 1];
  if(pooling_type == 0)
    theta_k[1] = eta[1];
  if(pooling_type == 1){
    if(joint_prior_variance)
      tau[1] = diag_pre_multiply(hypersd[1], L_Omega[1]);
    else
      tau[1] = diag_matrix(hypersd[1]);
    theta_k[1] = rep_matrix(mu[1], K) + tau[1] * eta[1];
  }
}
model {
  vector[N] y_mean;
  if(N > 0){

    // Impact of covariates:
    if(Nc == 0)
      y_mean = rep_vector(0.0, N);
    else{
      y_mean = X*beta;
      target += prior_increment_vec(prior_beta_fam, beta, prior_beta_val);
    }

    // Impact of baseline:
    if(pooling_baseline != 2)
      y_mean += to_vector(theta_k[1][1, site]);
    else
      y_mean += mu[1][1];

    // Impact of treatment:
    if(pooling_type != 2)
      y_mean += treatment .* to_vector(theta_k[1][2, site]);
    else
      y_mean += treatment * mu[1][2];
  }

  /* PRIORS */
  if(pooling_type > 0){
    if(joint_prior_mean)
      mu[1] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
    // else {
    //   target += prior_increment_real(prior_control_fam,   mu[1][1], prior_control_val);
    //   target += prior_increment_real(prior_hypermean_fam, mu[1][2], prior_hypermean_val);
    // }
  } else {
    for(k in 1:K)
      eta[1][,k] ~ multi_normal(prior_hypermean_mean, prior_hypermean_scale);
  }

  if(pooling_type == 1) {
    to_vector(eta[1]) ~ std_normal();
    target += prior_increment_real(prior_hypersd_fam,    hypersd[1][1], prior_hypersd_val);
    target += prior_increment_real(prior_control_sd_fam, hypersd[1][2], prior_control_sd_val);
    if(joint_prior_variance)
      L_Omega[1] ~ lkj_corr_cholesky(prior_hypercor_val);
  }


  /* LIKELIHOOD */
  y ~ normal(y_mean, sigma_y_k[site]);
}

generated quantities {
  // to do this, we must first (outside of Stan) calculate SEs in each test group,
  // i.e. test_sigma_y_k
  real logpd[K_test > 0];
  vector[N_test] fe_test;
  if(K_test > 0){
    if(Nc == 0)
      fe_test = rep_vector(0.0, N_test);
    else
      fe_test = X_test*beta;
    logpd[1] = 0;
    for(i in 1:N_test){
      if(pooling_type == 1)
        logpd[1] += normal_lpdf(test_y[i] | mu[1][1] + mu[1][2] * test_treatment[i] + fe_test[i],
                                sqrt(hypersd[1][1]^2 + hypersd[1][2]^2 + test_sigma_y_k[test_site[i]]^2));
      if(pooling_type == 2)
        logpd[1] += normal_lpdf(test_y[i] | mu[1][1] + mu[1][2] * test_treatment[i] + fe_test[i],
                                sqrt(test_sigma_y_k[test_site[i]]^2));
    }
  }
}

