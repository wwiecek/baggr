# for compilation of individual models use something like this:
recompile_one <- function(model){
  baggr_models[[model]] <-
    rstan::stan_model(stanc_ret = rstan::stanc(
      paste0("inst/models/", model, ".stan")))
  save(baggr_models, file = "~/baggr/precompiled_models.Rdata")
}
