# fancy display names in the package
model_names <- c("mutau" = "Aggregate data (with control group)",
                 "mutau_full" = "Individual-level data with control and treatment",
                 "rubin" = "Rubin model with aggregate data",
                 "rubin_full" = "Rubin model with individual-level data",
                 "quantiles" = "Effects on quantiles",
                 "logit" = "Logistic model (individual-level data)",
                 "sslab" = "Spike & slab model (3 components)")

model_data_types <- c("rubin" = "pool_noctrl_narrow",
                      "mutau" = "pool_wide",
                      "logit" = "individual_binary",
                      "rubin_full"  = "individual",
                      "mutau_full"  = "individual",
                      "sslab" = "individual",
                      #for now no quantiles model from summary level data
                      "quantiles" = "individual")

data_type_names <- c("pool_noctrl_narrow" = "Aggregate (effects only)",
                     "pool_wide" = "Aggregate (control and effects)",
                     "individual" = "Individual-level with continuous outcome",
                     "individual_binary" = "Individual-level with binary outcome")

available_priors <- list(
  "real" = c("normal", "uniform", "cauchy"),
  "positive_real" = c("normal", "uniform", "cauchy"),
  "real_2" = c("multinormal"),
  "corr" = c("lkj")
)
