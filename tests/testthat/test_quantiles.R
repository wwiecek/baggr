context("baggr() calls with quantiles model")
library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

# This example is based on Rachael Meager's testing code
K <- 7
N0 <- 50*K
N1 <- 51*K
sites0 <- rep(seq(1:7), N0/K)
sites1 <- rep(seq(1:7), N1/K)
y0 <- rnorm(N0,50 + sites0, 10)
y1 <- rnorm(N1,60 + sites1, 20)
chosen_quantiles <- c(0.2,0.4,0.5)
Q <- length(chosen_quantiles)
N1 <- length(y1)
N0 <- length(y0)
df_quantiles <- rbind(
  data.frame(outcome = y0, treatment = 0, group = sites0),
  data.frame(outcome = y1, treatment = 1, group = sites1)
)
# df_quantiles %>% dplyr::group_by(group, treatment) %>%
#   dplyr::summarise(m=quantile(outcome, .2)) %>% tidyr::spread(treatment, m)

# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  expect_error(baggr(df_quantiles[1,]), "data argument")
  expect_error(baggr(df_quantiles, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_quantiles, pooling = "nune"), "Wrong pooling")

  # NA or NULL inputs
  df_na <- df_quantiles; df_na$outcome[1] <- NA
  expect_error(baggr(df_na),"NA")
  df_na <- df_quantiles; df_na$treatment[2] <- NA
  expect_error(baggr(df_na),"NA")
  df_na <- df_quantiles; df_na$group[5] <- NA
  expect_error(baggr(df_na),"NA")
  df_na <- df_quantiles; df_na$group <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_quantiles; df_na$outcome <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_quantiles; df_na$treatment <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_quantiles; df_na$outcome <- as.character(df_na$outcome)
  expect_error(baggr(df_na),"numeric")
  expect_error(baggr(df_quantiles, group = "state1000"), "no column")

  expect_error(expect_warning(
    baggr(df_quantiles, model = "quantiles", rubbish = 41)), "unknown")

  # Improper quantiles spec:
  expect_error(convert_inputs(df_quantiles, "quantiles"), "quantiles")
  expect_error(convert_inputs(df_quantiles, "quantiles", quantiles = c(.5)), "less than 2")
  expect_error(baggr(df_quantiles, "quantiles", quantiles = c(.5)), "less than 2")
  expect_error(convert_inputs(df_quantiles, "quantiles", quantiles = c(.5, 1.1)), "must be between 0 and 1")
  expect_error(baggr(df_quantiles, "quantiles", quantiles = c(.5, 1.1)), "must be between 0 and 1")

  # test_that("Converting inputs works correctly") more explicitly
  ci <- suppressWarnings(convert_inputs(df_quantiles, "quantiles", quantiles = c(.2, .4, .5)))
  expect_length(ci, 15)
  expect_identical(names(ci)[1:4],
                   c("y_0", "y_1", "Sigma_y_k_0", "Sigma_y_k_1"))


})


test_that("All basic quantile models tests", {

  # skip("Test")
  skip_on_cran()

  # There will always be divergent transitions / ESS warning produced by Stan
  # at iter = 200.
  bg5_n <- expect_warning(baggr(df_quantiles, "quantiles", pooling = "none",
                                quantiles = chosen_quantiles,
                                iter = 200, chains = 2, refresh = 0,
                                show_messages = F))
  bg5_p <- expect_warning(baggr(df_quantiles, "quantiles", pooling = "partial",
                                quantiles = chosen_quantiles,
                                iter = 200, chains = 2, refresh = 0,
                                show_messages = F))
  bg5_f <- expect_warning(baggr(df_quantiles, "quantiles", pooling = "full",
                                quantiles = chosen_quantiles,
                                iter = 200, chains = 2, refresh = 0,
                                show_messages = F))
  bg5_ppd <- expect_warning(baggr(df_quantiles, "quantiles", ppd = TRUE,
                                  quantiles = chosen_quantiles,
                                  iter = 200, chains = 2, refresh = 0,
                                  show_messages = F))
  bg5_labels <- expect_warning(baggr(df_quantiles, "quantiles", pooling = "partial",
                                     quantiles = chosen_quantiles,
                                     effect = "Special", iter = 200, chains = 2, refresh = 0))

  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
  expect_is(bg5_labels, "baggr")

  expect_equal(nrow(as.matrix(bg5_p$fit)), 200) #right dimension means right iter


  # Various attr of baggr object are correct
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, length(chosen_quantiles))
  expect_equal(bg5_p$n_groups, 7)
  expect_equal(bg5_p$effects, c("20% quantile mean", "40% quantile mean", "50% quantile mean"))
  expect_equal(bg5_p$model, "quantiles")
  expect_equal(bg5_labels$effect, c("20% quantile on Special",
                                    "40% quantile on Special",
                                    "50% quantile on Special"))
  expect_error(
    suppressWarnings(baggr(df_quantiles, "quantiles", pooling = "partial",
                           quantiles = chosen_quantiles, effect = c("L1", "L2"))),
    "length")
  expect_is(bg5_p$fit, "stanfit")

  # Data are available in baggr object
  expect_is(bg5_n$data, "data.frame")
  expect_is(bg5_p$data, "data.frame")
  expect_is(bg5_f$data, "data.frame")
  expect_identical(bg5_f$data, bg5_p$data)

  # Pooling metrics
  # all pooling metric are the same as SE's are the same
  expect_equal(dim(bg5_p$pooling_metric), c(3,7,3))
  expect_equal(dim(bg5_f$pooling_metric), c(3,7,3))
  expect_equal(dim(bg5_n$pooling_metric), c(3,7,3))
  expect_equal(length(unique(bg5_p$pooling_metric[1,,1])), 7) #expect_length()
  expect_equal(length(unique(bg5_p$pooling_metric[2,,1])), 7)
  expect_equal(length(unique(bg5_p$pooling_metric[3,,1])), 7)
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
  expect_equal(length(unique(pp[2,,1])), 7)
  expect_gt(as.numeric(pp[2,1,1]), .5)


  # Calculation of effects works
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(200, 7, 3)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(200, 7, 3)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(200, 7, 3)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))


  # Plotting works
  expect_is(plot(bg5_n), "list")
  expect_is(plot(bg5_p, order = TRUE), "list")
  expect_is(plot(bg5_f, order = FALSE), "list")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")

  # printing works
  expect_output(print(bg5_n))
  expect_output(print(bg5_p))
  expect_output(print(bg5_p, group = FALSE))
  capture_output(expect_error(print(bg5_p, group = "abc"), "logical"))
  expect_output(print(bg5_f))
  expect_output(print(bg5_ppd))

  # Try this:
  expect_output(print(bg5_p, exponent = TRUE))
  expect_output(print(bg5_p, digits = 3))

  #Forest plots for quantiles model
  # expect_is(forest_plot(bg5_n), "vpPath")
  # expect_is(forest_plot(bg5_p), "vpPath")
  # expect_is(forest_plot(bg5_p, show = "posterior"), "vpPath")
  # expect_is(forest_plot(bg5_p, show = "both"), "vpPath")
  # expect_is(forest_plot(bg5_f), "vpPath")
  # expect_is(forest_plot(bg5_f, graph.pos = 1), "vpPath")
  # expect_error(forest_plot(cars), "baggr objects")
  # expect_error(forest_plot(bg5_p, show = "abc"), "should be one of")



  # Extracting treatment/study effects works
  expect_error(treatment_effect(df_quantiles))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "matrix") #this might change to accommodate more dim's
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")

  # Drawing values of tau:
  # expect_error(effect_draw(cars))
  # expect_is(effect_draw(bg5_p), "numeric")
  # expect_length(effect_draw(bg5_p), 200)
  # expect_length(effect_draw(bg5_p,7), 7)
  # expect_identical(effect_draw(bg5_n), NA)

  # Plotting tau:
  # expect_is(effect_plot(bg5_p), "gg")
  # expect_is(effect_plot(bg5_p, bg5_f), "gg")
  # expect_is(effect_plot("Model A" = bg5_p, "Model B" = bg5_f), "gg")
  # Crashes when passing nonsense
  # expect_error(effect_plot(cars), "baggr class")
  # expect_error(effect_plot(cars, cars, bg5_f), "baggr class")


  # test_that("baggr_compare basic cases work with quantiles", {
  #  # try to make nonexistant comparison:
  #   expect_error(baggr_compare(bg5_p, bg5_n, bg5_f, compare = "sreffects"),
  #                "argument must be set")
  #  # Compare prior vs posterior:
  #   bgcomp <- expect_warning(baggr_compare(df_quantiles, iter = 200, model = "quantiles",
  #                                          what = "prior", refresh = 0))
  #   expect_is(bgcomp, "baggr_compare")
  #   # Compare existing models:
  #   bgcomp2 <- baggr_compare(bg5_p, bg5_n, bg5_f)
  #   expect_is(bgcomp2, "baggr_compare")
  # })
  #
  # test_that("loocv", {
  #   # Rubbish model
  #   expect_error(loocv(schools, model = "rubbish"))
  #   # Can't do pooling none
  #   expect_error(loocv(schools, pooling = "none"))
  #
  #   loo_model <- expect_warning(loocv(schools, return_models = TRUE, iter = 200, refresh = 0))
  #   expect_is(loo_model, "baggr_cv")
  #   capture_output(print(loo_model))
  # })

  # Plot quantiles
  expect_error(plot_quantiles(cars), "a baggr")
  expect_is(plot_quantiles(bg5_labels), "list")
  expect_is(plot_quantiles(bg5_labels, ncol  = 3), "list")

})



# Test data -----

# This should work:


test_that("Test data can be used in the quantiles model", {

  df_notest <- df_quantiles[df_quantiles$group < 6,]
  df_test <- df_quantiles[df_quantiles$group >= 3,]

  skip_on_cran()

  bg_lpd <- expect_warning(baggr(df_notest, test_data = df_test, model = "quantiles",
                                 quantiles = chosen_quantiles, iter = 200, refresh = 0))


  # Wrong data type:
  expect_error(baggr(data = df_quantiles, test_data = cars), "is of type")

  # NA values or wrong cols:
  df2 <- df_quantiles; df2$group[1] <- NA
  expect_error(
    suppressWarnings(baggr(data = df_quantiles, test_data = df2, model = "quantiles")),
    "NA")


  expect_is(bg_lpd, "baggr")
  # make sure that we have 5 sites, not 7:
  expect_equal(dim(group_effects(bg_lpd)), c(400, 5, 3))
  # make sure it's not 0
  expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd[1]")[[1]]), -4000, tolerance = 500)

})


# covariates ------
test_that("Model with covariates works fine", {
  df2 <- df_quantiles
  df2$x <- rnorm(nrow(df2))
  # expect_error(baggr(df2, covariates = c("made_up_covariates"), model = "quantiles"), "made_up_covariates")
  expect_error(
    suppressWarnings(baggr(df2, covariates = c("x"), model = "quantiles")),
    "Quantiles model cannot regress on covariates.")
})




# Helper: summarise quantiles data -----

test_that("Summarising quantiles works", {
  expect_error(summarise_quantiles_data(df_quantiles,  quantiles = c(.1, .5, .1)))
  expect_error(summarise_quantiles_data(df_quantiles,  quantiles = c(.1, .5, .5)))
  expect_error(summarise_quantiles_data(df_quantiles,  quantiles = c(.1)))
  expect_length(
    suppressWarnings(summarise_quantiles_data(df_quantiles, quantiles = c(.1, .5, .9), means_only = TRUE)), 2)
})

