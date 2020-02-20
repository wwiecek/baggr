context("baggr() calls with logistic regression model")
library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

df_binary <- data.frame(treatment = rbinom(400, 1, .5),
                      group = rep(paste("Trial", LETTERS[1:5]), each = 80))
df_binary$outcome <- ifelse(df_binary$treatment, rbinom(400, 1, .3), rbinom(400, 1, .15))


# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  expect_error(baggr(df_binary, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_binary, pooling = "nune"), "Wrong pooling")

  # test_that("Converting inputs works correctly") more explicitly
  expect_identical(names(convert_inputs(df_binary, "logit")),
                   c("K", "N", "P", "y", "treatment", "site", "N_test", "K_test",
                     "test_y", "test_site", "test_treatment", "Nc", "X", "X_test"))
})

bg5_n <- expect_warning(baggr(df_binary, "logit", pooling = "none",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))
bg5_p <- expect_warning(baggr(df_binary, "logit", pooling = "partial",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))
bg5_f <- expect_warning(baggr(df_binary, "logit", pooling = "full",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))
bg5_ppd <- expect_warning(baggr(df_binary, "logit", ppd = T,
                                iter = 200, chains = 2, refresh = 0,
                                show_messages = F))

test_that("Different pooling methods work for Rubin model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
})

test_that("Extra args to Stan passed via ... work well", {
  expect_equal(nrow(as.matrix(bg5_p$fit)), 200) #right dimension means right iter
  expect_error(baggr(df_binary, rubbish = 41))
})

test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, 5)
  expect_equal(bg5_p$effects, "logOR")
  expect_equal(bg5_p$model, "logit")
  expect_is(bg5_p$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_is(bg5_n$data, "data.frame")
  expect_is(bg5_p$data, "data.frame")
  expect_is(bg5_f$data, "data.frame")
})

test_that("Pooling metrics", {
  # all pooling stats are 0 if no pooling
  expect_equal(unique(as.numeric(bg5_n$pooling_metric)), 0)
  # full pooling means 1's everywhere
  expect_equal(unique(as.numeric(bg5_f$pooling_metric)), 1)

  # pp <- pooling(bg5_p)
  # expect_is(pp, "array")
  # expect_gt(min(pp), 0)
  # expect_lt(max(pp), 1)
  # expect_identical(bg5_p$pooling_metric, pooling(bg5_p))

  # since all SEs are the same, pooling should be the same for all sites
  # capture_output(print(pp))
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 200 iter
  # expect_equal(length(unique(pp[2,,1])), 1)
  # expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)
})


test_that("Calculation of effects works", {
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")
  expect_length(treatment_effect(bg5_p, summary = T)$tau, 5)
  expect_length(treatment_effect(bg5_p, summary = T)$sigma_tau, 5)

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(200, 5 , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(200, 5 , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(200, 5 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting works", {
  expect_is(plot(bg5_ppd), "gg")
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, transform = exp), "gg")
  expect_is(plot(bg5_p, transform = exp, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  expect_is(forest_plot(bg5_n), "vpPath")
  expect_is(forest_plot(bg5_p), "vpPath")
  expect_is(forest_plot(bg5_f), "vpPath")
  expect_is(forest_plot(bg5_f, graph.pos = 1), "vpPath")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")
})

test_that("printing works", {
  capture_output(print(bg5_n))
  capture_output(print(bg5_p))
  capture_output(print(bg5_f))
  capture_output(print(bg5_p, exponent = TRUE))
})

# test_that("Test data can be used in the Rubin model", {
#   bg_lpd <- expect_warning(baggr(df_binary[1:6,], test_data = df_binary[7:8,],
#                                  iter = 500, refresh = 0))
#   expect_is(bg_lpd, "baggr")
#   # make sure that we have 6 sites, not 8:
#   expect_equal(dim(group_effects(bg_lpd)), c(1000, 6, 1))
#   # make sure it's not 0
#   expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd[1]")[[1]]), -3.6, tolerance = 1)
#
#   # wrong test_data
#   df_na <- df_binary[7:8,]; df_na$tau <- NULL
#   expect_error(baggr(df_binary[1:6,], test_data = df_na), "must be of the same format as input")
# })


# test helpers -----

test_that("Extracting treatment/study effects works", {
  expect_error(treatment_effect(df_binary))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "numeric") #this might change to accommodate more dim's
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")

  # Drawing values of tau:
  expect_error(effect_draw(cars))
  expect_is(effect_draw(bg5_p), "numeric")
  expect_length(effect_draw(bg5_p), 200)
  expect_length(effect_draw(bg5_p,7), 7)
  expect_identical(effect_draw(bg5_n), NA)

  # Plotting tau:
  expect_is(effect_plot(bg5_p), "gg")
  expect_is(effect_plot(bg5_p, bg5_f), "gg")
  expect_is(effect_plot("Model A" = bg5_p, "Model B" = bg5_f), "gg")
  # Crashes when passing nonsense
  expect_error(effect_plot(cars), "baggr class")
  expect_error(effect_plot(cars, cars, bg5_f), "baggr class")
})



# covariates ------
sa <- df_binary
sa$a <- rnorm(nrow(df_binary))
sa$b <- rnorm(nrow(df_binary))
sb <- sa
sb$b <- NULL
bg_cov <- baggr(sa, covariates = c("a", "b"), iter = 200, refresh = 0)

test_that("Model with covariates works fine", {
  expect_is(bg_cov, "baggr")
  expect_error(baggr(sa, covariates = c("made_up_covariates")), "are not columns")
  expect_error(baggr(sa, covariates = c("a", "b", "made_up_covariates")))
  expect_length(bg5_p$covariates, 0)
  expect_length(bg_cov$covariates, 2)
  expect_null(bg_cov$mean_lpd)

  # Fixed effects extraction
  expect_is(fixed_effects(bg_cov), "matrix")
  expect_is(fixed_effects(bg_cov, transform = exp), "matrix")
  expect_equal(dim(fixed_effects(bg_cov, summary = TRUE)), c(2,5,1))
  expect_equal(dim(fixed_effects(bg_cov, summary = FALSE))[2], 2)
})



# tests for helper functions -----

test_that("baggr_compare basic cases work with logit models", {
  # If I pass nothing
  expect_error(baggr_compare(), "Must provide baggr models")
  # pooling
  expect_error(baggr_compare(schools, pooling = "full"))
  # if I pass rubbish
  expect_error(baggr_compare(cars))
  # if I pass list of rubbish
  expect_error(baggr_compare("Fit 1" = cars, "Fit 2" = cars))
  # try to make nonexistant comparison:
  expect_error(baggr_compare(bg5_p, bg5_n, bg5_f, compare = "sreffects"))
  # Run models from baggr_compare:
  bgcomp <- expect_warning(baggr_compare(schools,
                                         iter = 200, refresh = 0))
  expect_is(bgcomp, "baggr_compare")
  # Compare prior vs posterior:
  bgcomp <- expect_warning(baggr_compare(schools, iter = 200,
                                         what = "prior", refresh = 0))
  expect_is(bgcomp, "baggr_compare")
  # Compare existing models:
  bgcomp2 <- plot(baggr_compare(bg5_p, bg5_n, bg5_f), arrange = "single")
  # bgcomp3 <- baggr_compare(bg5_p, bg5_n, bg5_f, arrange = "grid")
  expect_is(bgcomp2, "plot_list")
  expect_is(bgcomp2[[1]], "gg")

})

test_that("loocv", {
  # Rubbish model
  # expect_error(loocv(schools, model = "mutau"))
  # Can't do pooling none
  expect_error(loocv(schools, pooling = "none"))

  loo_model <- expect_warning(loocv(schools, return_models = TRUE, iter = 200, refresh = 0))
  expect_is(loo_model, "baggr_cv")
  capture_output(print(loo_model))
})
