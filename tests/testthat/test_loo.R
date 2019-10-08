context("test loo in baggr")

# brms testing script
# library(brms)
library(baggr)

set.seed(1999)

# fit <- brm(tau | se(se) ~ 1 + (1 | group),
#            data = schools,
#            control = list(adapt_delta = 0.95),
#            prior = c(set_prior("normal(0,100)", class = "Intercept"),
#                      set_prior("uniform(0, 104.44)", class = "sd")))
#
# brms_kfold <- kfold(fit, group = "group")

# output from brms kfold
brms_kfold <- list(estimates = structure(c(-30.9625079222176, NA,
                                           61.9250158444352,
                                           1.00614837109732, NA,
                                           2.01229674219465),
                                         .Dim = 3:2,
                                         .Dimnames = list(
                               c("elpd_kfold", "p_kfold", "kfoldic"),
                               c("Estimate", "SE"))),
                   pointwise = structure(c(-4.60651920550889,
                                           -3.46670445367466,
                                           -4.0169751024801,
                                           -3.52425367649299,
                                           -3.82156646558605,
                                           -3.68970464199433,
                                           -3.92034704268037,
                                           -3.91643733380024),
                                         .Dim = c(8L, 1L),
                                         .Dimnames = list(
                                           NULL, "elpd_kfold")))

baggr_kfold <- loocv(schools, control = list(adapt_delta = 0.99),
                     iter = 5000)

tol = 0.5
test_that(desc = "baggr and brms are at least close", {
  # cross-validation scores
  expect_lt(brms_kfold$estimates[1,1] - baggr_kfold$elpd, 1)


})
