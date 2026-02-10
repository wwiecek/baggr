context("group_effects meta-regression covariates in full-data models")
library(baggr)

set.seed(2026)

make_continuous_ipd <- function() {
  k <- 6
  n_per_group <- 24
  group <- rep(seq_len(k), each = n_per_group)
  treatment <- rbinom(k * n_per_group, 1, 0.5)

  x_fixed_values <- seq(-1, 1, length.out = k)
  x_fixed <- x_fixed_values[group]
  x_varying <- rnorm(k * n_per_group)

  baseline <- rep(rnorm(k, 0, 0.3), each = n_per_group)
  tau_group <- rep(rnorm(k, 0.4, 0.2), each = n_per_group)
  outcome <- baseline + treatment * tau_group + 0.7 * x_fixed + 0.5 * x_varying + rnorm(k * n_per_group, 0, 0.6)

  data.frame(
    group = as.character(group),
    treatment = treatment,
    outcome = outcome,
    x_fixed = x_fixed,
    x_varying = x_varying
  )
}

make_binary_ipd <- function() {
  k <- 6
  n_per_group <- 40
  group <- rep(seq_len(k), each = n_per_group)
  treatment <- rbinom(k * n_per_group, 1, 0.5)

  x_fixed_values <- seq(-1, 1, length.out = k)
  x_fixed <- x_fixed_values[group]
  x_varying <- rnorm(k * n_per_group)

  alpha_group <- rep(rnorm(k, -1.8, 0.3), each = n_per_group)
  tau_group <- rep(rnorm(k, 0.8, 0.2), each = n_per_group)
  lp <- alpha_group + treatment * tau_group + 0.7 * x_fixed + 0.5 * x_varying
  outcome <- rbinom(k * n_per_group, 1, plogis(lp))

  data.frame(
    group = as.character(group),
    treatment = treatment,
    outcome = outcome,
    x_fixed = x_fixed,
    x_varying = x_varying
  )
}

check_group_effect_decomposition <- function(bg) {
  ge_all <- group_effects(bg, random_only = FALSE)[,,1]
  ge_random <- group_effects(bg, random_only = TRUE)[,,1]
  x_group <- baggr:::group_effects_covariate_matrix(bg, n_groups = ncol(ge_random))
  fe_component <- fixed_effects(bg) %*% t(x_group)

  expect_equal(ge_all - ge_random, fe_component, tolerance = 1e-10)
  expect_true("x_fixed" %in% colnames(bg$summary_data))
  expect_false("x_varying" %in% colnames(bg$summary_data))
}

test_that("rubin_full and mutau_full add study-level fixed effects in group_effects", {
  data_cont <- make_continuous_ipd()

  fit_rubin_full <- expect_warning(
    baggr(
      data_cont,
      model = "rubin_full",
      covariates = c("x_fixed", "x_varying"),
      pooling = "partial",
      iter = 120,
      chains = 2,
      refresh = 0,
      show_messages = FALSE
    )
  )

  fit_mutau_full <- expect_warning(
    baggr(
      data_cont,
      model = "mutau_full",
      covariates = c("x_fixed", "x_varying"),
      pooling = "partial",
      iter = 120,
      chains = 2,
      refresh = 0,
      show_messages = FALSE
    )
  )

  check_group_effect_decomposition(fit_rubin_full)
  check_group_effect_decomposition(fit_mutau_full)
})

test_that("logit model adds study-level fixed effects in group_effects", {
  data_bin <- make_binary_ipd()

  fit_logit <- expect_warning(
    baggr(
      data_bin,
      model = "logit",
      covariates = c("x_fixed", "x_varying"),
      pooling = "partial",
      iter = 120,
      chains = 2,
      refresh = 0,
      show_messages = FALSE
    )
  )

  check_group_effect_decomposition(fit_logit)
})
