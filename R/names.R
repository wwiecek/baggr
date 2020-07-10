# fancy display names in the package
model_names <- c("mutau" = "Aggregate data (with control group)",
                 "rubin" = "Rubin model with aggregate data",
                 "full" = "Rubin model with individual-level data",
                 "quantiles" = "Effects on quantiles",
                 "logit" = "Logistic model (individual-level data)",
                 "sslab" = "Spike & slab model (3 components)")

model_data_types <- c("rubin" = "pool_noctrl_narrow",
                      "mutau" = "pool_wide",
                      "logit" = "individual_binary",
                      "full"  = "individual",
                      "sslab" = "individual",
                      #for now no quantiles model from summary level data
                      "quantiles" = "individual")

data_type_names <- c("pool_noctrl_narrow" = "Aggregate (effects only)",
                     "pool_wide" = "Aggregate (control and effects)",
                     "individual" = "Individual-level with continuous outcome",
                     "individual_binary" = "Individual-level with binary outcome")
