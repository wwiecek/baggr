# baggr models compilation, most basic version
# involves external storage

model_precompilation <- function(
  installation_path = "~/baggr/models_stan/")
{
  # library(rstan)


  wd <- getwd()
  # this is risky, if foo crashes before we revert to old wd we're stuck
  # but for now let's say it's OK...
  setwd(installation_path)

  lf <- list.files()
  baggr_models <- list()
  if(length(lf) == 0)
    stop("No models to compile. Stopping.")

  for(nm in lf) {
    if(grepl(".stan", nm)) {
      sm <- rstan::stan_model(stanc_ret = rstan::stanc(nm))
      stan_inputs <- list(K = 2, tau_hat_k = c(0,0), se_k = c(1,1))
      assign(gsub(".stan", "", nm),
             # sampling(sm, data = stan_inputs, iter = 100))
             sm)
      baggr_models[[gsub(".stan", "", nm)]] <- get(gsub(".stan", "", nm))
    }}

  rm(stan_inputs, sm)
  save(baggr_models, file = "../precompiled_models.Rdata")
  setwd(wd)
}
