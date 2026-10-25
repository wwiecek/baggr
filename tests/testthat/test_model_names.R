context("model names and posterior parameter notation")

model_summary <- data.frame(group = LETTERS[1:4], yi = c(1, 2, -1, 3),
                            sei = c(.4, .6, .5, .3),
                            y0i = c(8, 12, 10, 11), se0i = rep(.3, 4))
model_individual <- data.frame(group = rep(LETTERS[1:4], each = 12),
                               treatment = rep(rep(0:1, each = 6), 4),
                               outcome = rep(1:6, 8) +
                                 rep(c(0, 1, 2, 4, 3, 2, 4, 7), each = 6))
model_aliases <- list(rubin = c("rubin", "nnhm"),
                      vivs = c("vivs", "mutau", "bivariate"),
                      vivs_full = c("vivs_full", "mutau_full", "bivariate_full"))

test_that("aliases normalise before conversion and prior construction", {
  for(model in names(model_aliases)) {
    data <- if(model == "vivs_full") model_individual else model_summary
    if(model == "rubin") data <- data[, c("group", "yi", "sei")]
    input <- convert_inputs(data, model, silent = TRUE)
    prior <- prepare_prior(list(), data, input, model, silent = TRUE)
    for(alias in model_aliases[[model]]) {
      converted <- convert_inputs(data, alias, silent = TRUE)
      expect_identical(converted, input)
      expect_identical(attr(converted, "model"), model)
      expect_identical(prepare_prior(list(), data, input, alias, silent = TRUE),
                       prior)
    }
  }
  expect_identical(attr(convert_inputs(model_summary, NULL, silent = TRUE),
                        "model"), "vivs")
  expect_message(convert_inputs(model_summary[, 1:3], NULL),
                 "Normal-normal hierarchical \\(Rubin\\)")
  expect_message(convert_inputs(model_summary, NULL),
                 "Varying-intercept, varying-slope")
})

test_that("aliases fit the same model with unchanged posterior draws", {
  for(model in names(model_aliases)) {
    data <- if(model == "vivs_full") model_individual else model_summary
    reference <- NULL
    reference_draws <- NULL
    for(alias in model_aliases[[model]]) {
      fit <- suppressWarnings(baggr(data, model = alias, iter = 100,
        warmup = 50, chains = 1, seed = 1909, refresh = 0,
        silent = TRUE, warn = FALSE))
      expect_true(all(c("theta", "sigma", "theta_k") %in%
                        fit$fit@sim$pars_oi))
      expect_false(any(c("mu", "tau", "hypersd") %in% fit$fit@sim$pars_oi))
      expect_identical(fit$fit@stanmodel@model_name, model)
      posterior <- treatment_effect(fit, summary = TRUE)
      draws <- rstan::extract(fit$fit, permuted = FALSE)
      if(is.null(reference_draws)) reference_draws <- draws else
        expect_identical(draws, reference_draws, info = paste(model, alias))
      expect_identical(fit$model, model)
      expect_named(posterior, c("theta", "sigma"))
      expect_identical(posterior$theta, hypermean(fit))
      expect_identical(posterior$sigma, hypersd(fit))
      if(is.null(reference)) reference <- posterior else
        expect_equal(posterior, reference, tolerance = 1e-12,
                     info = paste(model, alias))
      if(model != "rubin")
        expect_identical(vivs_cor(fit), mutau_cor(fit))
    }
  }
})

test_that("LOO summaries use theta and sigma for canonical models and aliases", {
  references <- list()
  for(model in c("rubin", "nnhm", "vivs", "mutau", "bivariate")) {
    canonical <- canonical_model(model)
    data <- if(canonical == "rubin") model_summary[, 1:3] else model_summary
    cv <- suppressWarnings(suppressMessages(loocv(data, model = model,
      iter = 100, warmup = 50, chains = 1, seed = 1909,
      silent = TRUE, warn = FALSE)))
    expect_named(cv$df, c("theta", "sigma", "lpd"))
    expect_identical(cv$full_model$model, canonical)
    # Permuting identical draws can change rounding in their means.
    if(is.null(references[[canonical]])) references[[canonical]] <- cv$df else
      expect_equal(cv$df, references[[canonical]], tolerance = 1e-12,
                   info = paste(model, canonical))
  }
})

test_that("individual VIVS inputs keep treatment vectors for held-out data", {
  train <- model_individual[model_individual$group != "D", ]
  held_out <- model_individual[model_individual$group == "D", ]
  for(alias in model_aliases$vivs_full) {
    empty <- convert_inputs(train, alias, silent = TRUE)
    expect_identical(dim(empty$test_treatment), 0L)
    input <- convert_inputs(train, alias, test_data = held_out, silent = TRUE)
    expect_identical(dim(input$test_treatment), nrow(held_out))
    expect_equal(as.vector(input$test_treatment), held_out$treatment)
    fit <- suppressWarnings(baggr(train, model = alias, test_data = held_out,
      iter = 100, warmup = 50, chains = 1, seed = 1909, refresh = 0,
      silent = TRUE, warn = FALSE))
    expect_true(is.finite(fit$mean_lpd))
  }
})

test_that("Rubin and logistic individual models retain treatment matrices", {
  for(model in c("rubin_full", "logit")) {
    data <- model_individual
    if(model == "logit") data$outcome <- as.integer(data$outcome %% 2 == 0)
    empty <- convert_inputs(data, model, silent = TRUE)
    expect_identical(dim(empty$treatment), c(nrow(data), 1L))
    expect_identical(dim(empty$test_treatment), c(0L, 1L))
    for(n in c(2L, 3L)) {
      held_out <- data[which(data$treatment == 1)[seq_len(n)], ]
      input <- suppressMessages(convert_inputs(data, model,
        test_data = held_out, silent = TRUE))
      expect_identical(dim(input$treatment), c(nrow(data), 1L))
      expect_identical(dim(input$test_treatment), c(n, 1L))
      expect_equal(as.vector(input$test_treatment), held_out$treatment)
    }
    # Multiple contrasts still require matrices even without held-out rows.
    data$treatment <- factor(rep(c("control", "A", "B"), length.out = nrow(data)),
                             levels = c("control", "A", "B"))
    input <- convert_inputs(data, model, silent = TRUE)
    expect_identical(dim(input$treatment), c(nrow(data), 2L))
    expect_identical(dim(input$test_treatment), c(0L, 2L))
    expect_equal(unname(input$treatment),
                 unname(model.matrix(~ treatment, data)[, -1, drop = FALSE]))
  }
})
