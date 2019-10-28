context("baggr() calls with Rubin model")
library(baggr)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

# pooled, with equal SE's!
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])

# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  expect_error(baggr("Text"), "data argument")
  expect_error(baggr(5), "data argument")
  expect_error(baggr(df_pooled[1,]), "data argument")

  # model or pooling type doesn't exist
  expect_error(baggr(df_pooled, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_pooled, pooling = "nune"), "Wrong pooling")

  # NA or NULL inputs
  df_na <- df_pooled; df_na$tau[1] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_pooled; df_na$se[2] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_pooled; df_na$se <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_pooled; df_na$tau <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_pooled; df_na$tau <- as.character(df_na$tau)
  expect_error(baggr(df_na),"are not numeric")


  # test_that("Converting inputs works correctly") more explicitly
  expect_identical(names(convert_inputs(df_pooled, "rubin")),
                   c("K", "tau_hat_k", "se_tau_k", "K_test", "test_tau_hat_k", "test_se_k"))

  expect_warning(baggr(df_pooled, group = "state1000", iter = 50, refresh = 0),
                 "No labels will be added.")

})


# There will always be a divergent transition / ESS warning produced by Stan
# at iter = 200.
bg5_n <- expect_warning(baggr(df_pooled, "rubin", pooling = "none", group = "state",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))
bg5_p <- expect_warning(baggr(df_pooled, "rubin", pooling = "partial", group = "state",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))
bg5_f <- expect_warning(baggr(df_pooled, "rubin", pooling = "full", group = "state",
                              iter = 200, chains = 2, refresh = 0,
                              show_messages = F))

test_that("Different pooling methods work for Rubin model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
})

test_that("Extra args to Stan passed via ... work well", {
  expect_equal(nrow(as.matrix(bg5_p$fit)), 200) #right dimension means right iter
  expect_error(baggr(df_pooled, rubbish = 41))
})

test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, 8)
  expect_equal(bg5_p$effects, "mean")
  expect_equal(bg5_p$model, "rubin")
  expect_is(bg5_p$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_identical(bg5_n$data, df_pooled)
  expect_identical(bg5_p$data, df_pooled)
  expect_identical(bg5_f$data, df_pooled)
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
  expect_identical(bg5_p$pooling_metric, pooling(bg5_p))

  # since all SEs are the same, pooling should be the same for all sites
  capture_output(print(pp))
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 200 iter
  expect_equal(length(unique(pp[2,,1])), 1)
  expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)
})


test_that("Calculation of effects works", {
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(200, 8 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting works", {
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")
})

test_that("printing works", {
  capture_output(print(bg5_n))
  capture_output(print(bg5_p))
  capture_output(print(bg5_f))
})

test_that("Test data can be used in the Rubin model", {
  bg_lpd <- expect_warning(baggr(df_pooled[1:6,], test_data = df_pooled[7:8,],
                                 iter = 500, refresh = 0))
  expect_is(bg_lpd, "baggr")
  # make sure that we have 6 sites, not 8:
  expect_equal(dim(group_effects(bg_lpd)), c(1000, 6, 1))
  # make sure it's not 0
  expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd")[[1]]), -3.6, tolerance = 1)

  # wrong test_data
  df_na <- df_pooled[7:8,]; df_na$tau <- NULL
  expect_error(baggr(df_pooled[1:6,], test_data = df_na), "must be of the same format as input")
})


# test helpers -----

test_that("Extracting treatment/study effects works", {
  expect_error(treatment_effect(df_pooled))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "numeric") #this might change to accommodate more dim's
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")
})



# to-do list for tests -----

# show_model()
# For this we need a pre-commit hook to copy models from src/ to inst/models
# check_columns()
# convert_inputs(): if(required_data != available_data)
# detect_input_type()
# mint()
# prepare_ma()
# print_baggr()
# group_effects() (above)

# v0.2
# baggr_compare()
# loocv()
# plot_quantiles()
# summarise_quantiles_data()
# mutau, full, qunatiles in
#   baggr, study effects, trt effects, convert_inputs, pooling_metrics
# baggr_plot() with multiple effects



# tests for helper functions -----

test_that("baggr_compare", {
  # If I pass nothing
  expect_error(baggr_compare(), "Must provide baggr models")
  # pooling
  expect_error(baggr_compare(schools, pooling = "full"))
  # if I pass rubbish
  expect_error(baggr_compare(cars))
  # if I pass list of rubbish
  expect_error(baggr_compare("Fit 1" = cars, "Fit 2" = cars))
  # Run models from baggr_compare:
  bgcomp <- baggr_compare(schools, refresh = 0)
  expect_is(bgcomp, "list")
  # Compare existing models:
  bgcomp2 <- baggr_compare(bg5_p, bg5_n, bg5_f, arrange = "single")
  bgcomp3 <- baggr_compare(bg5_p, bg5_n, bg5_f, arrange = "grid")
  expect_is(bgcomp2, "gg")
  expect_is(bgcomp3, "list")
  expect_is(bgcomp3[[1]], "gg")

})

test_that("loocv", {
  # Rubbish model
  expect_error(loocv(schools, model = "mutau"))
  # Can't do pooling none
  expect_error(loocv(schools, pooling = "none"))

  loo_model <- expect_warning(loocv(schools, return_models = TRUE, iter = 200, refresh = 0))
  expect_is(loo_model, "baggr_cv")
  capture_output(print(loo_model))
})





# 8 schools correctness test -----

bg_s <- baggr(schools, refresh = 0, iter = 2000, control = list(adapt_delta = .99))

test_that("The default 8 schools result is close to the result in BDA", {
  expect_equal(mean(treatment_effect(bg_s)$tau), 8, tolerance = .25)
  expect_equal(mean(treatment_effect(bg_s)$sigma_tau), 7, tolerance = .25)
})


# Bangert-Drowns correctness check -----

# dt_bd <- metafor::dat.bangertdrowns2004 %>%
# dplyr::mutate(study = paste(author, year), tau = yi, se = sqrt(vi))
# write.table(dt_bd, file = "inst/tests/bangertdrowns2004.csv")
dt_bd <- read.table(system.file("tests", "bangertdrowns2004.csv", package = "baggr"))
bg_bd <- baggr(dt_bd, group = "study", model = "rubin", iter = 4000,
               prior_hypersd = uniform(0, 10), refresh = 0)
# RE results from metafor:
# metafor::rma(yi, vi, data=dt2)
test_that("Bangert-Drowns meta-analysis result is close to metafor output", {
  expect_equal(mean(treatment_effect(bg_bd)$tau), 0.22, tolerance = .01)
  expect_equal(as.numeric(quantile(treatment_effect(bg_bd)$tau, .025)),
               0.13, tolerance = .015)
  expect_equal(as.numeric(quantile(treatment_effect(bg_bd)$tau, .975)),
               0.31, tolerance = .015)
})
