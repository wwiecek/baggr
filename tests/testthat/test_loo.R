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

baggr_kfold <- expect_warning(loocv(schools,
                                    # control = list(adapt_delta = 0.9),
                                    iter = 5000))

# baggr_ranef <- group_effects(baggr_fit, summary = T)[]
test_that(desc = "baggr and brms are at least close", {
  # cross-validation scores
  expect_lt(brms_kfold$estimates[1,1] - baggr_kfold$elpd, 1)
})

# should be 0
comp <- loo_compare(baggr_kfold, baggr_kfold)

# test various things about the loo comparison method
test_that(desc = "loo_compare method works", {
  expect_error(loo_compare(baggr_kfold, brms_kfold))
  expect_error(loo_compare(list(baggr_kfold, brms_kfold)))
  expect_equal(comp[,1], 0)
  expect_equal(comp[,2], 0)
  expect_is(comp, "compare_baggr_cv")
  capture_output(print(comp))
})



