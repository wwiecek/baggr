context("baggr() calls with IPD version of Rubin model")
library(baggr)

skip_on_cran()

set.seed(1990)

# Generate 8 schools like IPD data
schools_ipd <- data.frame()
N <- c(rep(10, 4), rep(20, 4))
bsl <- vector(length=8)
for(i in 1:8){
  bsl[i] <- rnorm(1, 0, 5)

  x <- rnorm(N[i])
  x <- (x-mean(x))/sd(x)
  x <- x*schools$se[i]*sqrt(N[i])/1.41 + schools$tau[i]

  y <- rnorm(N[i])
  y <- (y-mean(y))/sd(y)
  y <- y*schools$se[i]*sqrt(N[i])/1.41

  schools_ipd <- rbind(schools_ipd,
                       data.frame(group = schools$group[i], outcome = bsl[i] + x, treatment = 1),
                       # This is just so that we don't trip off prepare_ma:
                       data.frame(group = schools$group[i], outcome = bsl[i] + y, treatment = 0))
}



bg_n <- expect_warning(baggr(schools_ipd, pooling = "none", iter = 150, refresh=0))
bg_p <- expect_warning(baggr(schools_ipd, pooling = "partial", iter = 150, refresh=0))
bg_f <- expect_warning(baggr(schools_ipd, pooling = "full", iter = 150, refresh=0))


test_that("Different pooling methods work for the rubin_full model", {
  expect_is(bg_n, "baggr")
  expect_is(bg_p, "baggr")
  expect_is(bg_f, "baggr")
})

test_that("Basic operations on rubin_full model", {
  expect_error(baggr(schools_ipd, rubbish = 41))
  expect_is(pooling(bg_p)[,,1], "matrix")
  expect_is(plot(bg_p), "gg")
  expect_is(effect_plot(bg_p), "gg")
  expect_is(funnel_plot(bg_p), "gg")
  expect_is(forest_plot(bg_p), "gforge_forestplot")
  bgc <- try(baggr_compare(bg_n, bg_p, bg_f))
  expect_is(bgc, "baggr_compare")

  # No pooling gives sensible results:
  expect_equal(
    as.numeric(group_effects(bg_n,s=T)[,"mean",1]),
    schools$tau, tolerance = 1)

})

test_that("extra pooling stats work", {
  # Extra pooling checks
  # Calculation of I^2 and H^2
  i2 <- pooling(bg_p, metric = "isq")
  expect_is(i2, "array")
  expect_gte(min(i2), 0)
  expect_lte(max(i2), 1)
  h2 <- pooling(bg_p, metric = "hsq")
  expect_is(h2, "array")
  expect_gte(min(h2), 1)
  # Calculation of weights makes sense
  wt <- weights(bg_p)
  expect_is(wt, "array")
  expect_equal(dim(wt), c(3,8,1))
  expect_equal(sum(wt[2,,1]), 1)
  expect_lte(sum(wt[1,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[3,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[1,,1]), 0)
  wt2 <- pooling(bg_p, metric = "weights")
  expect_identical(wt, wt2)
})

test_that("rubin_full model crashes with nonsense inputs", {
  expect_error(baggr(schools_ipd, outcome = 2), "Arguments")
  expect_error(baggr(schools_ipd, group = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = "wrong"), "no column")
  schools_ipd2 <- schools_ipd; schools_ipd2$treatment <- as.character(schools_ipd$treatment)
  expect_error(baggr(schools_ipd2), "has to be numeric")
  schools_ipd2 <- schools_ipd; schools_ipd2$outcome <- as.character(schools_ipd$outcome)
  expect_error(baggr(schools_ipd2), "has to be numeric")

})

test_that("Using old syntax (model = full) still works", {
  bg <- expect_warning(expect_message(baggr(schools_ipd,
                                            model = "full",
                                            pooling = "none", iter = 10, refresh=0)))
  expect_is(bg, "baggr")
})


comp_flpl <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 150, what = "pooling"
))

comp_flpr <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 150, what = "prior"
))

comp_flmdls <- baggr_compare(bg_f, bg_p)

test_that("baggr comparison method works for rubin_full model", {

  expect_is(comp_flpl, "baggr_compare")
  expect_is(comp_flpr, "baggr_compare")
  expect_is(comp_flmdls, "baggr_compare")

  expect_is(testthat::capture_output(print(comp_flpl)), "character")
  expect_is(testthat::capture_output(print(comp_flpr)), "character")
  expect_is(testthat::capture_output(print(comp_flmdls)), "character")

  expect_gt(length(comp_flpl), 0)
  expect_gt(length(comp_flpr), 0)
  expect_gt(length(comp_flmdls), 0)

  expect_is(plot(comp_flpl), "gg")

  expect_is(plot(comp_flpr), "ggplot")

  expect_is(plot(comp_flmdls), "gg")

  expect_is(plot(comp_flpl, grid_models = TRUE), "gtable")

  expect_is(plot(comp_flpr, grid_models = TRUE), "gtable")

  expect_is(plot(comp_flmdls, grid_models = TRUE), "gtable")
})

test_that("rubin_full cross-validation works", {
  # Run it first with test data that includes baseline, this will gen a message:
  bg <- expect_warning(expect_message(
    baggr(subset(schools_ipd, group != "School A"), iter = 20, refresh = 0,
          test_data = subset(schools_ipd, group == "School A")),
    "Baselines for all these groups"
  )
  )

  # Now repeat without bsl data
  bg <- expect_warning(baggr(subset(schools_ipd, group != "School A"), iter = 20, refresh = 0,
                             test_data = subset(schools_ipd, group == "School A" & treatment == 1)))
  expect_is(bg, "baggr")
  expect_gt(bg$mean_lpd, 0)
})

test_that("rubin_full CV handles numeric/character groups and K_test is correct", {
  run_cv_case <- function(df, held_out_groups) {
    train_data <- subset(df, !(group %in% held_out_groups) | treatment == 0)
    test_data <- subset(df, group %in% held_out_groups & treatment == 1)

    bg <- expect_warning(
      baggr(train_data, iter = 20, chains = 1, refresh = 0, test_data = test_data)
    )

    expect_is(bg, "baggr")
    expect_true(is.finite(bg$mean_lpd))
    expect_equal(bg$inputs$K_test, length(unique(test_data$group)))
    expect_equal(length(unique(bg$inputs$test_site)), bg$inputs$K_test)
  }

  schools_ipd_num <- schools_ipd
  schools_ipd_num$group <- as.numeric(factor(schools_ipd_num$group, levels = unique(schools_ipd_num$group)))

  run_cv_case(schools_ipd, "School A")
  run_cv_case(schools_ipd, c("School A", "School C", "School E"))
  run_cv_case(schools_ipd_num, 1)
  run_cv_case(schools_ipd_num, c(1, 3, 5))
})



# covariates ------
sa <- schools_ipd
sa$a <- rnorm(nrow(schools_ipd))
sa$b <- rnorm(nrow(schools_ipd))
sa$c <- sample(c("A", "B", "C"), nrow(schools_ipd), replace = T)

sb <- sa
sb$b <- NULL


test_that("Model with covariates works fine", {
  bg_cov <- expect_warning(
    baggr(sa, covariates = c("a", "b", "c"), iter = 150, chains = 1, refresh = 0))

  expect_is(bg_cov, "baggr")
  expect_error(baggr(sa, covariates = c("made_up_covariates")), "are not columns")
  expect_error(baggr(sa, covariates = c("a", "b", "made_up_covariates")), "are not columns")
  expect_length(bg_p$covariates, 0)
  expect_length(bg_cov$covariates, 3)
  expect_null(bg_cov$mean_lpd)

  # Fixed effects extraction
  expect_is(fixed_effects(bg_cov), "matrix")
  expect_is(fixed_effects(bg_cov, transform = exp), "matrix")
  expect_equal(dim(fixed_effects(bg_cov, summary = TRUE)), c(4,5,1))
  expect_equal(dim(fixed_effects(bg_cov, summary = FALSE))[2], 4)
})



test_that("Fixed within-study covariates are added to summary_data", {
  sa_fixed <- schools_ipd
  sa_fixed$site_cov <- as.numeric(factor(sa_fixed$group))
  bg_cov_fixed <- expect_warning(
    baggr(sa_fixed, covariates = c("site_cov"), iter = 150, chains = 1, refresh = 0)
  )

  expect_true("site_cov" %in% names(bg_cov_fixed$summary_data))
  expected <- sa_fixed$site_cov[match(bg_cov_fixed$summary_data$group, sa_fixed$group)]
  expect_equal(bg_cov_fixed$summary_data$site_cov, expected)
})


test_that("Within-study varying covariates are reported as not meta-regression", {
  sa_vary <- schools_ipd
  sa_vary$ind_cov <- rnorm(nrow(sa_vary))

  bg_cov_vary <- expect_message(
    suppressWarnings(
      baggr(sa_vary, covariates = c("ind_cov"),
            iter = 150, chains = 1, refresh = 0, warn = FALSE)
    ),
    "Covariate ind_cov varies within studies. Model fitting will work but is not a meta-regression."
  )
  expect_is(bg_cov_vary, "baggr")
})
bg_pr <- expect_warning(baggr(schools_ipd,
                              pooling = "partial",
                              pooling_control = "remove",
                              iter = 150,
                              refresh=0))
bg_pn <- expect_warning(baggr(schools_ipd,
                              pooling = "partial",
                              pooling_control = "none",
                              iter = 150,
                              refresh=0))
test_that("You can change pooling on values in control group", {
  expect_is(bg_pr, "baggr")
  expect_is(bg_pn, "baggr")
  bsl_k <- apply(rstan::extract(bg_pr$fit, "baseline_k")[[1]], 2, mean)
  expect_length(bsl_k, 8)
  expect_equal(bsl_k, rep(0, 8))
  bsl_k <- apply(rstan::extract(bg_pn$fit, "baseline_k")[[1]], 2, mean)
  expect_length(bsl_k, 8)
  expect_equal(bsl_k, bsl, tolerance = 1)
})
