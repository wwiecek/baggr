context("baggr() calls with logistic regression model")
library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

df_binary <- data.frame(treatment = rbinom(400, 1, .5),
                        group = rep(paste("Trial", LETTERS[1:5]), each = 80))
df_binary$outcome <- ifelse(df_binary$treatment, rbinom(400, 1, .3), rbinom(400, 1, .15))
df_summ <- prepare_ma(df_binary, "logOR")

# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  expect_error(baggr(df_binary, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_binary, pooling = "nune"), "should be one of")

  # test_that("Converting inputs works correctly") more explicitly
  expect_identical(names(convert_inputs(df_binary, "logit")),
                   c("K", "N", "P", "y", "treatment", "site", "N_test", "K_test",
                     "test_y", "test_site", "test_treatment", "Nc", "X", "X_test"))
})

bg5_n <- expect_warning(baggr(df_binary, "logit", pooling = "none",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
bg5_p <- expect_warning(baggr(df_binary, "logit", pooling = "partial",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
bg5_f <- expect_warning(baggr(df_binary, "logit", pooling = "full",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
bg5_ppd <- expect_warning(baggr(df_binary, "logit", ppd = TRUE,
                                iter = 150, chains = 2, refresh = 0,
                                show_messages = F))
bg5_summarydt <- expect_warning(baggr(df_summ,
                                      iter = 150, chains = 2, refresh = 0,
                                      show_messages = F))

test_that("Different pooling methods work for Rubin model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
  expect_is(bg5_summarydt, "baggr")
})

test_that("Extra args to Stan passed via ... work well", {
  expect_equal(nrow(as.matrix(bg5_p$fit)), 150) #right dimension means right iter
  expect_error(baggr(df_binary, rubbish = 41))
})

# Run without summarising data first
bg5_or <- expect_warning(baggr(df_summ[,c("group", "a", "n1", "c", "n2")],
                     iter = 150, chains = 2, refresh = 0,
                     show_messages = F))
# Same but now manually set the effect to RR
bg5_rr <- expect_warning(
  baggr(df_summ[,c("group", "a", "n1", "c", "n2")], effect = "logRR",
                iter = 150, chains = 2, refresh = 0,
                show_messages = F))

test_that("We can run models without prepare_ma'ing data", {
  expect_is(bg5_or, "baggr")
  expect_is(bg5_rr, "baggr")
  expect_equal(bg5_or$effects, "logOR")
  expect_equal(bg5_rr$effects, "logRR")
})

test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, 5)
  expect_equal(bg5_p$effects, "logOR")
  expect_equal(bg5_p$model, "logit")
  expect_equal(bg5_summarydt$model, "rubin")
  expect_equal(bg5_summarydt$pooling, "partial")
  expect_is(bg5_p$fit, "stanfit")
  expect_is(bg5_summarydt$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_is(bg5_n$data, "data.frame")
  expect_is(bg5_p$data, "data.frame")
  expect_is(bg5_f$data, "data.frame")
  expect_is(bg5_n$summary_data, "data.frame")
  expect_is(bg5_p$summary_data, "data.frame")
  expect_is(bg5_f$summary_data, "data.frame")
  expect_null(bg5_summarydt$summary_data)
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
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 150 iter
  # expect_equal(length(unique(pp[2,,1])), 1)
  # expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)
})

test_that("extra pooling stats work", {
  # Extra pooling checks
  # Calculation of I^2 and H^2
  i2 <- pooling(bg5_p, metric = "isq")
  expect_is(i2, "array")
  expect_gte(min(i2), 0)
  expect_lte(max(i2), 1)
  h2 <- pooling(bg5_p, metric = "hsq")
  expect_is(h2, "array")
  expect_gte(min(h2), 1)
  h <- pooling(bg5_p, metric = "h")
  expect_is(h, "array")
  expect_gte(min(h), 1)
  # Calculation of weights makes sense
  wt <- weights(bg5_p)
  expect_is(wt, "array")
  expect_equal(dim(wt), c(3,5,1))
  expect_equal(sum(wt[2,,1]), 1)
  expect_lte(sum(wt[1,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[3,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[1,,1]), 0)
  wt2 <- pooling(bg5_p, metric = "weights")
  expect_identical(wt, wt2)
})

test_that("Calculation of effects works", {
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")
  expect_length(hypermean(bg5_p,message=FALSE), 5) 
  expect_length(hypersd(bg5_p,message=FALSE), 5)

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(150, 5 , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(150, 5 , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(150, 5 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting and printing works", {
  expect_is(plot(bg5_ppd), "gg")
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, transform = exp), "gg")
  expect_is(plot(bg5_p, transform = exp, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  # This has been changed in forestplot 2.0:
  expect_is(forest_plot(bg5_n), "gforge_forestplot")
  expect_is(forest_plot(bg5_p), "gforge_forestplot")
  expect_is(forest_plot(bg5_f), "gforge_forestplot")
  expect_is(forest_plot(bg5_f, graph.pos = 1), "gforge_forestplot")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")

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
  expect_is(effect_draw(bg5_p, transform = exp), "numeric")
  expect_length(effect_draw(bg5_p), 150)
  expect_length(effect_draw(bg5_p,7), 7)

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


test_that("Model with covariates works fine", {
  bg_cov <- expect_warning(
    baggr(sa, covariates = c("a", "b"), iter = 150, chains = 1, refresh = 0))

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


# covariates, group-level ------
df_cov <- df_summ
df_cov$somecov <- rnorm(5)

test_that("Model with covariates works fine", {
  bg_cov2 <- expect_warning(
    baggr(df_cov, covariates = c("somecov"), iter = 150, chains = 1, refresh = 0))

  expect_is(bg_cov2, "baggr")
  expect_error(baggr(df_cov, covariates = c("made_up_covariates")), "are not columns")
  expect_length(bg_cov2$covariates, 1)
  expect_null(bg_cov2$mean_lpd)

  # Fixed effects extraction
  expect_is(fixed_effects(bg_cov2), "matrix")
  expect_is(fixed_effects(bg_cov2, transform = exp), "matrix")
  expect_equal(dim(fixed_effects(bg_cov2, summary = TRUE)), c(1,5,1))
  expect_equal(dim(fixed_effects(bg_cov2, summary = FALSE))[2], 1)
})



# tests for helper functions -----

test_that("baggr_compare basic cases work with logit models", {
  # try to make nonexistant comparison:
  expect_error(baggr_compare(bg5_p, bg5_n, bg5_f, compare = "sreffects"))
  # Compare existing models:
  expect_is(plot(baggr_compare(bg5_p, bg5_n, bg5_f)), "gg")
})

test_that("loocv", {
  # Rubbish model
  # expect_error(loocv(schools, model = "mutau"))
  # Can't do pooling none
  expect_error(loocv(df_binary, pooling = "none"))

  skip_on_cran()

  loo_model <- expect_warning(loocv(df_binary, model = "logit",
                                    return_models = TRUE, iter = 150, chains = 1, refresh = 0))
  expect_is(loo_model, "baggr_cv")
  capture_output(print(loo_model))
  expect_is(plot(loo_model), "gg")

  loo_full <- expect_warning(loocv(df_binary, model = "logit", pooling = "full",
                                   return_models = TRUE, iter = 150, chains = 1, refresh = 0))
  expect_is(loo_full, "baggr_cv")
  capture_output(print(loo_full))
  expect_is(plot(loo_full, add_values = FALSE), "gg")

  looc <- loo_compare(loo_model, loo_full)
  expect_is(looc, "compare_baggr_cv")
  capture_output(looc)
})

comp_pl <- expect_warning(baggr_compare(
  df_binary, model = "logit", iter = 150, what = "pooling"
))

comp_pr <- expect_warning(baggr_compare(
  df_binary, model = "logit", iter = 150, what = "prior"
))

comp_existmodels <- baggr_compare(bg5_p, bg5_f)

test_that("baggr comparison method works for Full model", {

  expect_is(comp_pl, "baggr_compare")
  expect_is(comp_pr, "baggr_compare")
  expect_is(comp_existmodels, "baggr_compare")

  expect_is(testthat::capture_output(print(comp_pl)), "character")
  expect_is(testthat::capture_output(print(comp_pr)), "character")
  expect_is(testthat::capture_output(print(comp_existmodels)), "character")

  expect_gt(length(comp_pl), 0)
  expect_gt(length(comp_pr), 0)
  expect_gt(length(comp_existmodels), 0)

  expect_is(plot(comp_pl), "gg")
  expect_is(plot(comp_pl, grid_models = TRUE), "gtable")
  expect_is(plot(comp_pr), "gg")
  expect_is(plot(comp_pr, grid_models = TRUE), "gtable")
  expect_is(plot(comp_existmodels), "gg")
  expect_is(plot(comp_existmodels, grid_models = TRUE), "gtable")

})



# Setting control pooling, control priors -----

test_that("Prior specifications for baselines work", {

  skip_on_cran()

  bg1 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
                              pooling_control = "partial",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
  bg2 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
                              pooling_control = "partial", prior_control = normal(0, 5),
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
  bg3 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
                              pooling_control = "partial", prior_control = normal(0, 5), prior_control_sd = uniform(0, 1),
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))

  expect_is(bg1, "baggr")
  expect_is(bg2, "baggr")
  expect_is(bg3, "baggr")

  expect_error(baggr(df_binary, "logit", pooling = "none",
                     pooling_control = "partial",
                     prior_control = multinormal(c(0,0), diag(2)),
                     iter = 150, chains = 2, refresh = 0,
                     show_messages = F))
})
