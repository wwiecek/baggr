functions {
  real normal_ss_log(int N, real y_sq_sum, vector xy_sum,
                     matrix xx_sum, vector beta, real sigma) {
    real beta_xy;
    real lp;
    beta_xy = dot_product(xy_sum, beta);
    lp = -.5*(y_sq_sum - 2*beta_xy + sum(beta * beta' .* xx_sum))/sigma^2;
    lp = lp - .5*N*log(sigma^2);
    return lp;
  }
}

data {
  int<lower=0> K;  // number of sites
  int<lower=0> N;  // total number of observations
  int<lower=0> P;  // dimensionality of parameter vector which  is
                   // jointly distributed - here, it is 2 dimensional
  real y[N];       // outcome variable of interest
  int ITT[N];      // intention to treat indicator
  int site[N];     // factor variable to split them out into K sites
  matrix[P,P] mutau_prior_sigma;
  vector[P] mutau_prior_mean;
}
transformed data {
  int N_k[K];           // number of observations from site K
  real y_sq_sum[K];     // sum_i y_{ki}^2
  vector[P] xy_sum[K];  // sum_i y_ki [1, ITT_{ki}]
  matrix[P,P] xx_sum[K];// sum_i [1, ITT_{ki}] [1, ITT_{ki}]'
  int s;
  vector[P] x;
  // initialize everything to zero
  N_k = rep_array(0, K);
  y_sq_sum = rep_array(0.0, K);
  xy_sum = rep_array(rep_vector(0.0, P), K);
  xx_sum = rep_array(rep_matrix(0.0, P, P), K);
  // x[1] is always 1
  x[1] = 1.0;
  for (n in 1:N) {
    s = site[n];
    x[2] = ITT[n];
    N_k[s] = N_k[s] + 1;
    y_sq_sum[s] = y_sq_sum[s] + y[n]^2;
    xy_sum[s] = xy_sum[s] + y[n]*x;
    xx_sum[s] = xx_sum[s] + x*x';
  }
}

parameters {
  vector[P] mutau;
  matrix[K,P] mutau_k;

  real<lower=0> sigma_y_k[K];
  corr_matrix[P] Omega;        //  correlation
  vector<lower=0>[P] theta;    //  scale
}
transformed parameters {
  matrix[P,P] sigma_mutau;
  sigma_mutau = quad_form_diag(Omega,theta);
}

model {
  //data variance priors
  sigma_y_k ~ uniform(0,100000);
  // sigma_y_k ~ inv_gamma(0.1,10);

  // parameter variance priors
  theta ~ cauchy(0,10);
  
  // theta ~ normal(0,100);
  Omega ~ lkj_corr(3);

  // hyperparameter priors
  mutau ~ multi_normal(mutau_prior_mean, mutau_prior_sigma);

  for (k in 1:K) {
    mutau_k[k] ~ multi_normal(mutau, sigma_mutau);
    target += normal_ss_log(N_k[k], y_sq_sum[k], xy_sum[k], 
                            xx_sum[k], mutau_k[k]', sigma_y_k[k]);
  }
}
generated quantities{
  vector[2] predicted_mutau_k;
  real signal_noise_ratio_mu;
  real signal_noise_ratio_tau;
  signal_noise_ratio_mu = mutau[1]/sqrt(sigma_mutau[1,1]);
  signal_noise_ratio_tau = mutau[2]/sqrt(sigma_mutau[2,2]);
  predicted_mutau_k = multi_normal_rng(mutau, sigma_mutau);

}
