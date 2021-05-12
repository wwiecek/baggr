library(rstan)
set.seed(1990)

convert_logit_2_prob <- function(vector_betas){
  K <- length(vector_betas)
  vector_betas <- vector_betas - vector_betas[K]
  probs <- rep(NA,K)
  denominator_excluding_K <- (1+sum(exp((vector_betas[1:(K-1)]))))
  probs[K] <- 1/denominator_excluding_K
  for(k in 1:(K-1)){
    probs[k] <- exp(vector_betas[k])/denominator_excluding_K
  }

  return(probs)
}


get_posterior_intervals_function <- function(vector_draws){
  posterior_interval <- quantile(vector_draws, c(0.025, 0.25, 0.5, 0.75,0.975))
  return(posterior_interval)
  print(posterior_interval)
}


### GENERATE TEST DATA

mu <- c(1,2)
tau <- c(1,0)
sd_mu <- c(1,2)
sd_tau <- c(0.5,1)
sigma_control <- c(1,1)
sd_sigma_control <- c(1,1)
sigma_TE <- c(1,1)
sd_sigma_TE <- c(0.5,0.5)
K <- 5
N <- K*300
x <- data.frame( rep(1,N), rbinom(N,1,0.5) )
site <- rep(seq(1:K),each = N/K)
beta <- matrix(c(0.5,0,0,0.5,-0.5,0), 3,2)
sigma <- matrix(rep(0.5,6),3,2)
control_cat_probs <- convert_logit_2_prob(beta[,1]) # functions as a sanity check
treatment_cat_probs <- convert_logit_2_prob(beta[,2]) # sanity!
stan_gen_data <- list(N = N, K = K, x = x, site = site, beta = beta[-3,], sigma = sigma,
                      mu = mu, tau = tau, sd_mu = sd_mu, sd_tau = sd_tau, sigma_control = sigma_control,
                      sd_sigma_control = sd_sigma_control, sigma_TE = sigma_TE, sd_sigma_TE = sd_sigma_TE)

model_code <- "
data {
int N; // number of observations
int K; // number of sites
vector[2] x[N]; // covariates, but i only allow a constant and a binary treatment indicator
int site[N]; // site indicator
real mu[2];
real tau[2];
real<lower=0> sd_mu[2];
real<lower=0> sd_tau[2];
real sigma_control[2];
real sigma_TE[2];
real<lower=0> sd_sigma_control[2];
real<lower=0> sd_sigma_TE[2];
matrix[2,2] beta; // the parent parameters minus the omitted beta category
matrix<lower=0>[3,2] sigma; // the set of logit parent variances (not a covariance matrix)

}

parameters {
matrix[K,2] mu_k;
matrix[K,2] tau_k;
matrix[K,2] sigma_control_k;
matrix[K,2] sigma_TE_k;
matrix[3,2] beta_k_raw[K]; // the hierarchical increments

}
transformed parameters{
matrix[3,2] beta_full;
matrix[3,2] beta_k[K];
beta_full = append_row(beta,rep_row_vector(0, 2));
for (m in 1:3){
 for (k in 1:K){
  for (p in 1:2){
   beta_k[k,m,p] = beta_full[m,p] + sigma[m,p]*beta_k_raw[k,m,p];
}}}
}
model {
  for (m in 1:3){
    for (k in 1:K){
      beta_k_raw[k,m] ~ normal(0,1);
  }}
  for (k in 1:K){
    for (i in 1:2){
      mu_k[k,i] ~ normal(mu[i],sd_mu[i]);
      tau_k[k,i] ~ normal(tau[i], sd_tau[i]);
      sigma_control_k[k,i] ~ normal(sigma_control[i],sd_sigma_control[i]);
      sigma_TE_k[k,i] ~ normal(sigma_TE[i], sd_sigma_TE[i]);
  }
}
}
generated quantities{
int cat[N]; // category indicator
real y[N];//
for (n in 1:N){
  cat[n] = categorical_logit_rng(beta_k[site[n]] * x[n]);
  if(cat[n]==1) y[n] = -1*lognormal_rng(mu_k[site[n],1] + tau_k[site[n],1]*x[n,2], exp(sigma_control_k[site[n],1] + sigma_TE_k[site[n],1]*x[n,2]));
  else if(cat[n] == 2) y[n] = 0;
  else y[n] = lognormal_rng( mu_k[site[n],2] + tau_k[site[n],2]*x[n,2], exp(sigma_control_k[site[n],2] + sigma_TE_k[site[n],2]*x[n,2]));
  }

}
"
sm <- stan_model(model_code = model_code)
generated_data_out <- sampling(sm,
                           data = stan_gen_data, iter = 1000, chains = 1)

tau_k_draws <- extract(generated_data_out, pars = "tau_k", permuted = FALSE, inc_warmup = FALSE)
single_tau_k_draw <- tau_k_draws[500,1,]
cat_draws <- extract(generated_data_out, pars = "cat", permuted = FALSE, inc_warmup = FALSE)
single_cat_draw <- cat_draws[500,1,]
y_draws <- extract(generated_data_out, pars = "y", permuted = FALSE, inc_warmup = FALSE)
single_y_draw <- y_draws[500,1,]
data_out <- data.frame(single_cat_draw, single_y_draw, x, site)

data <- data.frame(single_y_draw, x, site)
colnames(data) <- c("profit","constant","treatment", "site")
N <- length(data$profit) # number of draws from the distribution / dimensionality of data vector
M <- 3 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X
X <- cbind(rep(1,N), data$treatment)
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$profit[i] < 0){ cat[i] <- 1  }
  else if(identical(data$profit[i],0)){cat[i] <- 2}
  else{cat[i] <- 3}
}

# data <- data.frame(data, cat)
# saveRDS(data, file = "models_spike/spike-slab-data.rds")
data_spike <- data
usethis::use_data(data_spike, overwrite = TRUE)
