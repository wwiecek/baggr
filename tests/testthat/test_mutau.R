context("baggr() calls with mu and tau model")
library(baggr)
library(testthat)

# prepare inputs ----------------------------------------------------------
set.seed(1990)

# pooled, with equal SE's!
df_mutau <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                       "se.tau" = rep(1, 8),
                       "mu" = rnorm(8),
                       "se.mu" = rep(1, 8),
                       "state" = datasets::state.name[1:8])

#
# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  # model, data or pooling mismatch
  expect_error(baggr(df_mutau, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_mutau, pooling = "nune"), "Wrong pooling")

  # NA or NULL inputs
  df_na <- df_mutau; df_na$mu[1] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_mutau; df_na$se.mu[2] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_mutau; df_na$se.mu <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_mutau; df_na$mu <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_mutau; df_na$mu <- as.character(df_na$mu)
  expect_error(baggr(df_na),"are not numeric")

  expect_warning(baggr(df_mutau, group = "state1000", iter = 50, refresh = 0),
                 "No labels will be added.")
  expect_identical(names(convert_inputs(df_mutau, "mutau")),
                   c("K", "P", "theta_hat_k", "se_theta_k",
                     "K_test", "test_theta_hat_k", "test_se_theta_k", "Nc", "X", "X_test"))
})

bg5_n <- expect_warning(baggr(df_mutau, pooling = "none", group = "state",
               iter = 200, chains = 2, refresh = 0))
bg5_p <- expect_warning(baggr(df_mutau, pooling = "partial", group = "state",
               iter = 200, chains = 2, refresh = 0))
bg5_f <- expect_warning(baggr(df_mutau, pooling = "full", group = "state",
               iter = 200, chains = 2, refresh = 0))

test_that("Different pooling methods work for mu tau model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
})

test_that("Extra args to Stan passed via ... work well", {
  expect_equal(nrow(as.matrix(bg5_p$fit)), 200) #right dimension means right iter
  expect_error(baggr(df_mutau, rubbish = 41))
})

test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, 8)
  expect_equal(bg5_p$effects, "mean")
  expect_equal(bg5_p$model, "mutau")
  expect_is(bg5_p$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_is(bg5_n$data, "data.frame")
  expect_is(bg5_p$data, "data.frame")
  expect_is(bg5_f$data, "data.frame")
})

test_that("Pooling metrics", {
  # all pooling metric are the same as SE's are the same
  expect_equal(length(unique(bg5_p$pooling_metric[1,,1])), 1) #expect_length()
  expect_equal(length(unique(bg5_p$pooling_metric[2,,1])), 1)
  expect_equal(length(unique(bg5_p$pooling_metric[3,,1])), 1)
  # all pooling stats are 0 if no pooling
  expect_equal(unique(as.numeric(bg5_n$pooling_metric)), 0)
  # full pooling means 1's everywhere
  expect_equal(unique(as.numeric(bg5_f$pooling_metric)), 1)

  pp <- pooling(bg5_p)
  expect_is(pp, "array")
  expect_gt(min(pp), 0)
  expect_lt(max(pp), 1)
  # since all SEs are the same, pooling should be the same for all sites
  capture_output(pp)
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 200 iter
  expect_equal(length(unique(pp[2,,1])), 1)
  expect_equal(as.numeric(pp[2,1,1]), .7, tolerance = .2)
})


test_that("Calculation of effects works", {
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")
  expect_length(treatment_effect(bg5_p, summary = TRUE)$tau, 5)
  expect_length(treatment_effect(bg5_p, summary = TRUE)$sigma_tau, 5)

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(200, 8 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting works", {
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  expect_is(forest_plot(bg5_n), "vpPath")
  expect_is(forest_plot(bg5_p), "vpPath")
  expect_is(forest_plot(bg5_f), "vpPath")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")
})

test_that("Test data can be used in the mu tau model", {
  bg_lpd <- expect_warning(baggr(df_mutau[1:6,], test_data = df_mutau[7:8,],
                  iter = 2000, chains = 2, refresh = 0))
  expect_is(bg_lpd, "baggr")
  # make sure that we have 6 sites, not 8:
  expect_equal(dim(group_effects(bg_lpd)), c(2000, 6, 1))
  # make sure it's not 0 but something sensible
  expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd[1]")[[1]]), -13, tolerance = 1)

  # wrong test_data
  df_na <- df_mutau[7:8,]; df_na$tau <- NULL
  expect_error(baggr(df_mutau[1:6,], test_data = df_na))
})


# test helpers -----

test_that("Extracting treatment/study effects works", {
  expect_error(treatment_effect(df_mutau))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "numeric")
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")
})

comp_mt <- baggr_compare(
  bg5_p, bg5_f
)

test_that("baggr comparison method works for mu-tau models", {

  expect_is(comp_mt, "baggr_compare")
  expect_output(print(comp_mt))
  expect_gt(length(comp_mt), 0)

  expect_is(plot(comp_mt), "gg")
  expect_is(plot(comp_mt, grid_models = TRUE), "gtable")

})

