context("baggr() calls with multi-arm IPD setup")
library(baggr)

skip_on_cran()

set.seed(1990)

N  <- 1000
df <- data.frame(
  treatment = factor(sample(c("A", "B", "C"), N, replace = TRUE), levels = c("A", "B", "C")),
  group = sample(paste("Study", 1:10), N, replace = TRUE)
)
df$cl <- sample(1:10, N, replace = TRUE)
df$outcome_cont <- rnorm(N) +
  (df$treatment == "A") * 0.2 +
  (df$treatment == "B") * 0.4 +
  (df$treatment == "C") * 0.6
df$outcome_bin  <- 1 * (df$outcome_cont > 0.2)

bg_n <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "none", iter = 150, refresh = 0))
bg_p <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "partial", iter = 150, refresh = 0))
bg_f <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "full", iter = 150, refresh = 0))

bg_rubin_cluster <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "partial",
                                       iter = 150, refresh = 0, cluster = "cl"))

bg_logit <- expect_warning(baggr(df, outcome = "outcome_bin", model = "logit",
                                 pooling = "partial", iter = 150, refresh = 0, cluster = "cl"))

test_that("multi-arm setup produces baggr objects", {
  expect_s3_class(bg_n, "baggr")
  expect_s3_class(bg_p, "baggr")
  expect_s3_class(bg_f, "baggr")
  expect_s3_class(bg_rubin_cluster, "baggr")
  expect_s3_class(bg_logit, "baggr")
})

test_that("multi-arm outputs are correctly formed for rubin_full", {
  expect_equal(bg_p$n_groups, 10)
  expect_equal(bg_p$n_parameters, 2)
  expect_equal(length(bg_p$effects), 2)
  expect_equal(bg_p$effects, c("mean B", "mean C"))

  ge <- group_effects(bg_p, summary = TRUE)
  expect_true(is.array(ge))
  expect_equal(dim(ge), c(10, 5, 2))

  te <- treatment_effect(bg_p, summary= TRUE)
  expect_equal(dim(te$sigma), c(2,5))
  expect_equal(dim(te$sigma_tau), c(2,5))

  pm <- pooling(bg_p)
  expect_true(is.array(pm))
  expect_equal(dim(pm), c(3, 10, 2))

  wt <- weights(bg_p)
  expect_true(is.array(wt))
  expect_equal(dim(wt), c(3, 10, 2))
  expect_equal(sum(wt[2, , 1]), 1, tolerance = 1e-6)
  expect_equal(sum(wt[2, , 2]), 1, tolerance = 1e-6)

  bsl_k <- apply(rstan::extract(bg_p$fit, "baseline_k")[[1]], 2, mean)
  expect_length(bsl_k, 10)
  expect_false(anyNA(bsl_k))
})

test_that("multi-arm plotting and comparison helpers work", {
  plotbg <- plot(bg_p)
  expect_length(plot(bg_p), 2)
  expect_s3_class(plotbg[[1]], "gg")
  expect_s3_class(plotbg[[2]], "gg")
  expect_s3_class(effect_plot(bg_p), "gg")
  expect_s3_class(funnel_plot(bg_p), "gg")
  expect_s3_class(forest_plot(bg_p), "gforge_forestplot")

  bgc <- baggr_compare(bg_n, bg_p, bg_f)
  expect_s3_class(bgc, "baggr_compare")
  expect_s3_class(plot(bgc), "gg")
})

test_that("logit multi-arm model keeps core outputs and drops eta_cluster draws", {
  expect_equal(bg_logit$n_groups, 10)
  expect_equal(bg_logit$n_parameters, 2)
  expect_equal(bg_logit$effects, c("logOR B", "logOR C"))

  par_names <- names(rstan::extract(bg_logit$fit))
  expect_false("eta_cluster" %in% par_names)

  ge <- group_effects(bg_logit)
  expect_true(is.array(ge))
  expect_equal(dim(ge), c(10, 5, 2))

  te <- treatment_effect(bg_logit)
  expect_true(is.array(te))
  expect_equal(dim(te), c(5, 2))

  pm <- pooling(bg_logit)
  expect_true(is.array(pm))
  expect_equal(dim(pm), c(3, 10, 2))
})

test_that("clustered rubin_full model also drops eta_cluster draws", {
  expect_equal(bg_rubin_cluster$n_groups, 10)
  expect_equal(bg_rubin_cluster$n_parameters, 2)

  par_names <- names(rstan::extract(bg_rubin_cluster$fit))
  expect_false("eta_cluster" %in% par_names)

  ge <- group_effects(bg_rubin_cluster)
  expect_true(is.array(ge))
  expect_equal(dim(ge), c(10, 5, 2))
})
