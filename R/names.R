# fancy display names in the package
model_names <- c("mutau" = "Aggregate data (with control group)",
                 "mutau_full" = "Individual-level data with control and treatment",
                 "rubin" = "Rubin model with aggregate data",
                 "rubin_full" = "Rubin model with individual-level data",
                 "quantiles" = "Effects on quantiles",
                 "logit" = "Logistic model",
                 "sslab" = "Spike & slab model (3 components)")

data_type_names <- c("pool_noctrl_narrow" = "Aggregate (effects only)",
                     "pool_wide" = "Aggregate (control and effects)",
                     "pool_binary" = "Aggregate (event counts)",
                     "individual" = "Individual-level with continuous outcome",
                     "individual_binary" = "Individual-level with binary outcome")


# Default settings for matching models with input data types (convert_inputs)
model_data_types <- c("rubin" = "pool_noctrl_narrow",
                      "mutau" = "pool_wide",
                      "logit" = "individual_binary",
                      "rubin_full"  = "individual",
                      "mutau_full"  = "individual",
                      "sslab" = "individual",
                      #for now no quantiles model from summary level data
                      "quantiles" = "individual")

data_type_default_model <- c("pool_noctrl_narrow" = "rubin",
                             "pool_wide" = "mutau",
                             "pool_binary" = "logit",
                             "individual" = "rubin_full",
                             "individual_binary" = "logit")


# Allowed priors for different types of variables
available_priors <- list(
  "real" = c("normal", "uniform", "cauchy", "lognormal", "student_t"),
  "positive_real" = c("normal", "uniform", "cauchy", "lognormal", "student_t"),
  "real_2" = c("multinormal"),
  "corr" = c("lkj")
)

prior_dist_fam <- c("uniform" = 0,
                    "normal" = 1,
                    "cauchy" = 2,
                    "multinormal" = 3,
                    "lkj" = 4,
                    "lognormal" = 5,
                    "student_t" = 6)
