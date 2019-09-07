# brms testing script
library(brms)
library(baggr)

fit <- brm(tau | se(se) ~ 1 + (1 | group),
           data = schools,
           control = list(adapt_delta = 0.95),
           prior = c(set_prior("normal(0,100)", class = "Intercept"),
                     set_prior("uniform(0, 104.44)", class = "sd")),
           file = "misc/brms_kfold_test")

brms_kfold <- kfold(fit, group = "group")

baggr_kfold <- loocv(schools, control = list(adapt_delta = 0.99),
                     iter = 5000)


