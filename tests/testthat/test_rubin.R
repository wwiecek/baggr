context("baggr() calls with Rubin model")
library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

# pooled, with equal SE's!
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])

# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  expect_error(baggr("Text"), "it's not a data.frame")
  expect_error(baggr(5), "it's not a data.frame")
  expect_error(baggr(df_pooled[1,]), "specify hyper-SD prior manually")

  # model or pooling type doesn't exist
  expect_error(baggr(df_pooled, "made_up_model"), "Unrecognised model")
  expect_error(baggr(df_pooled, pooling = "nune"), "one of")

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
                   c("K", "theta_hat_k", "se_theta_k", "K_test",
                     "test_theta_hat_k", "test_se_theta_k", "Nc", "X", "X_test",
                     "M", "c"))

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
bg5_ppd <- expect_warning(baggr(df_pooled, "rubin", ppd = TRUE,
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
  expect_identical(bg5_p$pooling_metric, pooling(bg5_p))

  # since all SEs are the same, pooling should be the same for all sites
  capture_output(print(pp))
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 200 iter
  expect_equal(length(unique(pp[2,,1])), 1)
  expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)

  # Calculation of I^2 and H^2
  i2 <- pooling(bg5_p, metric = "isq")
  expect_is(i2, "array")
  expect_gte(min(i2), 0)
  expect_lte(max(i2), 1)
  h2 <- pooling(bg5_p, metric = "hsq")
  expect_is(h2, "array")
  expect_gt(min(h2), 1)

  # Calculation of weights makes sense
  wt <- weights(bg5_p)
  expect_is(wt, "array")
  expect_equal(dim(wt), c(3,8,1))
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

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(200, 8 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting works", {
  expect_is(plot(bg5_ppd), "gg")
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_p, style = "forest"), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "one of")
})


test_that("printing works", {
  capture_output(print(bg5_n))
  capture_output(print(bg5_p))
  capture_output(print(bg5_p, group = FALSE))
  expect_error(print(bg5_p, group = "abc"), "logical")
  capture_output(print(bg5_f))
  capture_output(print(bg5_ppd))
})

test_that("Forest plots for Rubin model", {
  expect_is(forest_plot(bg5_n), "gforge_forestplot")
  expect_is(forest_plot(bg5_p), "gforge_forestplot")
  expect_is(forest_plot(bg5_p, show = "posterior"), "gforge_forestplot")
  expect_is(forest_plot(bg5_p, show = "both"), "gforge_forestplot")
  expect_is(forest_plot(bg5_f), "gforge_forestplot")
  expect_is(forest_plot(bg5_f, graph.pos = 1), "gforge_forestplot")
  expect_error(forest_plot(cars), "baggr objects")
  expect_error(forest_plot(bg5_p, show = "abc"), "one of")
})
test_that("Test data can be used in the Rubin model", {
  # Wrong data type:
  expect_error(baggr(data = df_pooled, test_data = cars), "is of type")

  # NA values or wrong cols:
  df2 <- df_pooled; df2$tau[1] <- NA
  expect_error(baggr(data = df_pooled, test_data = df2), "NA")
  df_na <- df_pooled[7:8,]; df_na$tau <- NULL
  expect_error(baggr(df_pooled[1:6,], test_data = df_na), "is of type")

  # This should work:
  bg_lpd <- expect_warning(baggr(df_pooled[1:6,], test_data = df_pooled[7:8,],
                                 iter = 500, refresh = 0))
  expect_is(bg_lpd, "baggr")
  # make sure that we have 6 sites, not 8:
  expect_equal(dim(group_effects(bg_lpd)), c(1000, 6, 1))
  # make sure it's not 0
  expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd[1]")[[1]]), -3.6, tolerance = 1)

})


# covariates ------
sa <- schools
sa$a <- rnorm(8)
sa$b <- rnorm(8)
sa$f <- as.factor(c(rep("Yes", 4), rep("No", 4)))
sa$f2 <- as.factor(c(rep("Yes", 3), rep("No", 2), rep("Maybe", 3)))
sa$bad <- c(rnorm(7), NA)

sb <- sa
sb$b <- NULL
bg_cov <- expect_warning(
  baggr(sa, covariates = c("a", "b"), iter = 200, refresh = 0))
bg_cov_factor <- expect_warning(
  baggr(sa, covariates = c("f"), iter = 200, refresh = 0))
bg_cov_factor2 <- expect_warning(
  baggr(sa, covariates = c("f2"), iter = 200, refresh = 0))
expect_identical(attr(bg_cov_factor$inputs, "covariate_coding"), c("fYes"))
expect_identical(attr(bg_cov_factor2$inputs, "covariate_coding"), c("f2No", "f2Yes"))
bg_cov_test <- expect_warning(
  baggr(sa, covariates = c("a"), test_data = sb, iter = 200, refresh = 0))
bg_cov_prior1 <- expect_warning(
  baggr(sa, covariates = c("a", "b"),
        iter = 200, refresh = 0, prior_beta = normal(0, 3)))
bg_cov_prior2 <- expect_warning(
  baggr(sa, covariates = c("a", "b"),
        iter = 200, refresh = 0, prior = list("beta" = uniform(-5, 5))))

test_that("Model with covariates works fine", {
  expect_is(bg_cov, "baggr")
  expect_equal(bg_cov$formatted_prior$prior_beta_fam, 1)
  expect_equal(bg_cov$formatted_prior$prior_beta_val[1], 0)
  expect_gt(bg_cov$formatted_prior$prior_beta_val[2], 0)
  expect_equal(bg_cov$formatted_prior$prior_beta_val[3], 0)
  expect_error(baggr(sa, covariates = c("made_up_covariates")))
  expect_error(baggr(sa, covariates = c("a", "b", "made_up_covariates")))
  expect_error(baggr(sa, covariates = c("bad")), "NA")
  expect_length(bg5_p$covariates, 0)
  expect_length(bg_cov$covariates, 2)
  expect_null(bg_cov$mean_lpd)

  # Fixed effects extraction
  expect_is(fixed_effects(bg_cov), "matrix")
  expect_is(fixed_effects(bg_cov, transform = exp), "matrix")
  expect_equal(dim(fixed_effects(bg_cov, summary = TRUE)), c(2,5,1))
  expect_equal(dim(fixed_effects(bg_cov, summary = FALSE))[2], 2)

  # Bubble plots
  p1 <- bubble(bg_cov, "a")
  p2 <- bubble(bg_cov_factor, "f")
  p1
  p2
  expect_is(p1, "gg")
  expect_is(p2, "gg")

  # covariates and test_data
  expect_error(baggr(sa, covariates = c("a", "b"), test_data = sb), "Cannot bind")
  expect_error(baggr(sb, model = "rubin", test_data=sa[1:2,], covariates = c("b")), "are not columns")
  expect_is(bg_cov_test$mean_lpd, "numeric")
  expect_length(bg_cov_test$covariates, 1)

  # Setting priors for covariates manually works
  expect_is(bg_cov_prior1, "baggr")
  expect_is(bg_cov_prior2, "baggr")
  expect_equal(bg_cov_prior1$formatted_prior$prior_beta_fam, 1)
  expect_equal(bg_cov_prior1$formatted_prior$prior_beta_val, c(0,3,0))
  expect_equal(bg_cov_prior2$formatted_prior$prior_beta_fam, 0)
  expect_equal(bg_cov_prior2$formatted_prior$prior_beta_val, c(-5,5,0))
})

# test helpers -----

test_that("Extracting treatment/study effects works", {
  expect_error(treatment_effect(df_pooled))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "numeric") #this might change to accommodate more dim's
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")
  expect_length(treatment_effect(bg5_p, summary = TRUE)$tau, 5)

  # Drawing values of tau:
  expect_error(effect_draw(cars))
  expect_is(effect_draw(bg5_p), "numeric")
  expect_length(effect_draw(bg5_p), 200)
  expect_length(effect_draw(bg5_p,7), 7)
  eds1 <- effect_draw(bg5_p, summary = T)
  eds2 <- effect_draw(bg5_p, summary = T, interval = .5)
  expect_length(eds1, 5)
  expect_gt(eds2[1], eds1[1]) #narrower interval
  expect_warning(effect_draw(bg5_p, 1e05), "more effect draws than there are available samples")

  # Plotting tau:
  expect_is(effect_plot(bg5_p), "gg")
  expect_is(effect_plot(bg5_p, bg5_f), "gg")
  expect_is(effect_plot("Model A" = bg5_p, "Model B" = bg5_f), "gg")
  # Crashes when passing nonsense
  expect_error(effect_plot(cars), "baggr class")
  expect_error(effect_plot(cars, cars, bg5_f), "baggr class")
})


test_that("baggr_compare basic cases work with Rubin", {
  # If I pass nothing
  expect_error(baggr_compare(), "Must provide baggr models")
  # pooling
  expect_error(baggr_compare(schools, pooling = "full"))
  # if I pass rubbish
  expect_error(baggr_compare(cars))
  # if I pass list of rubbish
  expect_error(baggr_compare("Fit 1" = cars, "Fit 2" = cars))
  # try to make nonexistant comparison:
  expect_error(baggr_compare(bg5_p, bg5_n, bg5_f, compare = "sreffects"), "one of")
  # Run models from baggr_compare:
  bgcomp <- expect_warning(baggr_compare(schools,
                                         iter = 200, refresh = 0))
  expect_is(bgcomp, "baggr_compare")
  # Compare prior vs posterior:
  bgcomp <- expect_warning(baggr_compare(schools, iter = 200,
                                         what = "prior", refresh = 0))
  expect_is(bgcomp, "baggr_compare")
  # Compare existing models:
  expect_error(baggr_compare("Name" = bg5_p, "Name" = bg5_n), "unique model names")

  bgcomp2 <- baggr_compare(bg5_p, bg5_n, bg5_f)
  expect_is(bgcomp2, "baggr_compare")
  p1 <- plot(bgcomp2)
  p2 <- plot(bgcomp2, add_values = TRUE)
  expect_is(p1, "gg")
  expect_is(p2, "gg")
})

test_that("loocv", {
  skip_on_cran()

  # Rubbish model
  expect_error(loocv(schools, model = "rubbish"), "Inference failed")
  # Can't do pooling none
  expect_error(loocv(schools, pooling = "none"))

  loo_model <- expect_warning(loocv(schools, return_models = TRUE, iter = 200, refresh = 0))
  expect_is(loo_model, "baggr_cv")
  capture_output(print(loo_model))
})





# 8 schools correctness test -----


test_that("The default 8 schools result is close to the result in BDA", {
  skip_on_cran()
  bg_s <- baggr(schools, refresh = 0, iter = 2000, control = list(adapt_delta = .99))
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
               0.13, tolerance = .02)
  expect_equal(as.numeric(quantile(treatment_effect(bg_bd)$tau, .975)),
               0.31, tolerance = .02)
})


comp_rbpl <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 200, what = "pooling"
))

comp_rbpr <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 200, what = "prior"
))

test_that("baggr comparison method works for Rubin model", {

  expect_is(comp_rbpr, "baggr_compare")
  expect_is(comp_rbpl, "baggr_compare")

  expect_is(testthat::capture_output(print(comp_rbpl)), "character")
  expect_is(testthat::capture_output(print(comp_rbpr)), "character")

  expect_gt(length(comp_rbpl), 0)
  expect_gt(length(comp_rbpr), 0)

  expect_is(plot(comp_rbpl), "gg")

  expect_is(plot(comp_rbpr), "gg")

  expect_is(plot(comp_rbpl, grid_models = TRUE), "gtable")

  expect_is(plot(comp_rbpr, grid_models = TRUE), "gtable")
})

test_that("Plot quantiles", {
  expect_error(plot_quantiles(bg5_p))
})

test_that("You can run Rubin model with mutau type inputs", {
  df_mutau <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                         "se.tau" = rep(1, 8),
                         "mu" = rnorm(8),
                         "se.mu" = rep(1, 8),
                         "state" = datasets::state.name[1:8])
  bg <- expect_warning(baggr(df_mutau, model = "rubin", iter = 20, refresh = 0))
  expect_is(bg, "baggr")
  expect_equal(bg$model, "rubin")
  bg <- expect_warning(baggr(df_mutau[1:6,], test_data = df_mutau[7:8,],
                             model = "rubin", iter = 20, refresh = 0))
  expect_is(bg, "baggr")
  expect_gt(bg$mean_lpd, 0)
})
