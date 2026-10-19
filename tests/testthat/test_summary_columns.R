context("canonical summary-data columns")

summary_rubin <- data.frame(group = LETTERS[1:6], yi = c(2, 1, -1, 3, 0, 2),
                            sei = c(.3, .5, .4, .6, .5, .4))
summary_mutau <- transform(summary_rubin, y0i = 10:15, se0i = rep(.2, 6))
legacy_rubin <- summary_rubin
names(legacy_rubin)[2:3] <- c("tau", "se")
legacy_mutau <- summary_mutau
names(legacy_mutau)[2:5] <- c("tau", "se.tau", "mu", "se.mu")

# Stored data retain source columns; compare all actual Stan inputs separately.
stan_values <- function(x) {
  attributes(x) <- list(names = names(x))
  x
}

test_that("summary conventions yield the same Stan inputs", {
  expect_identical(normalise_summary_columns(summary_rubin)$data, summary_rubin)
  for(model in c("rubin", "vivs")) {
    canonical <- if(model == "rubin") summary_rubin else summary_mutau
    legacy <- if(model == "rubin") legacy_rubin else legacy_mutau
    expect_message(old <- convert_inputs(legacy, model), "Starting in baggr 0.9")
    new <- convert_inputs(canonical, model)
    expect_identical(stan_values(old), stan_values(new))
    expect_identical(attr(convert_inputs(canonical, NULL, silent = TRUE), "model"),
                     model)
    expect_message(convert_inputs(legacy, model, silent = TRUE), NA)
    expect_message(convert_inputs(legacy, model, yi = "tau",
      sei = if(model == "rubin") "se" else "se.tau"), NA)
  }
  expect_identical(detect_input_type(transform(summary_rubin, treatment = 1)),
                   "pool_noctrl_narrow")
  inp <- convert_inputs(summary_mutau, "vivs")
  expect_equal(inp$theta_hat_k[1, ], summary_mutau$y0i)
  expect_equal(inp$theta_hat_k[2, ], summary_mutau$yi)
  expect_equal(inp$se_theta_k[1, ], summary_mutau$se0i)
  expect_equal(inp$se_theta_k[2, ], summary_mutau$sei)
  expect_equal(stan_values(convert_inputs(summary_mutau, "rubin")),
               stan_values(convert_inputs(summary_rubin, "rubin")))
})

test_that("exact aliases and custom strings work", {
  aliases <- summary_rubin
  names(aliases)[2:3] <- c("y", "se")
  expect_message(alias <- convert_inputs(aliases, "rubin"), "preferred names")
  expect_equal(alias$theta_hat_k, summary_rubin$yi)
  wide <- summary_mutau
  names(wide)[2:5] <- c("y", "se", "y0", "se0")
  expect_message(alias <- convert_inputs(wide, "vivs"), "preferred names")
  expect_equal(alias$theta_hat_k[1, ], summary_mutau$y0i)
  custom <- summary_mutau
  names(custom)[2:5] <- c("effect", "effect_se", "baseline", "baseline_se")
  inp <- convert_inputs(custom, "vivs", yi = "effect", sei = "effect_se",
                         y0i = "baseline", se0i = "baseline_se")
  expect_identical(stan_values(inp),
                   stan_values(convert_inputs(summary_mutau, "vivs")))
  for(name in c("se_effect", "se_tau", "stderr")) {
    wrong <- summary_rubin
    names(wrong)[3] <- name
    expect_error(convert_inputs(wrong, "rubin"))
  }
})

test_that("conflicting conventions require an explicit mapping", {
  both <- transform(summary_rubin, tau = yi, se = sei, y = yi)
  expect_message(normalise_summary_columns(both), NA)
  expect_equal(normalise_summary_columns(both)$data$yi, summary_rubin$yi)
  both$tau[1] <- 99
  expect_error(normalise_summary_columns(both), "Conflicting summary columns")
  expect_equal(normalise_summary_columns(both, yi = "tau", sei = "se")$data$yi,
               both$tau)
  both$tau <- both$yi
  both$y[1] <- 99
  expect_error(normalise_summary_columns(both), "Conflicting summary columns")
})

test_that("mappings and held-out columns are validated early", {
  for(mapping in list(1:6, TRUE, c("yi", "sei"), NA_character_, ""))
    expect_error(convert_inputs(summary_rubin, "rubin", yi = mapping,
                                sei = "sei"), "single column name")
  expect_error(convert_inputs(summary_rubin, "rubin", yi = "yi"), "both")
  expect_error(convert_inputs(summary_rubin, "rubin", sei = "sei"), "both")
  expect_error(convert_inputs(summary_mutau, "vivs", y0i = "y0i",
                              se0i = "se0i"), "together")
  expect_error(convert_inputs(summary_rubin, "rubin", yi = "absent",
                              sei = "sei"), "Missing mapped summary columns")
  for(model in c("rubin", "vivs")) {
    df <- if(model == "rubin") legacy_rubin else legacy_mutau
    inp <- convert_inputs(df[-1, ], model, test_data = df[1, ], silent = TRUE)
    expect_equal(as.numeric(inp$test_theta_hat_k),
                 if(model == "rubin") summary_rubin$yi[1] else
                   c(summary_mutau$y0i[1], summary_mutau$yi[1]))
    expect_error(convert_inputs(df[-1, ], model,
      test_data = df[1, setdiff(names(df), "tau")], silent = TRUE), "test_data")
  }
})

test_that("slab column defaults apply only to summary data", {
  df <- summary_rubin
  names(df)[1] <- "slab"
  expect_message(inp <- convert_inputs(df, "rubin"), NA)
  expect_identical(attr(inp, "group_label"), df$slab)
  expect_message(inp <- convert_inputs(df, "rubin", group = "group"), NA)
  expect_null(attr(inp, "group_label"))
  df$group <- letters[1:6]
  expect_identical(attr(convert_inputs(df, "rubin"), "group_label"), df$group)
  ipd <- data.frame(outcome = 1:12, treatment = rep(0:1, 6),
                    group = rep(1:3, each = 4), slab = "unused")
  expect_message(inp <- convert_inputs(ipd, "rubin_full"), NA)
  expect_identical(attr(inp, "data"), ipd)
})

test_that("prepare_ma emits canonical columns with unchanged calculations", {
  ipd <- data.frame(group = rep(LETTERS[1:2], each = 8),
                    treatment = rep(rep(0:1, each = 4), 2),
                    outcome = c(1:4, 3:6, 2:5, 5:8))
  prepared <- prepare_ma(ipd)
  expect_true(all(c("y0i", "yi", "se0i", "sei", "n0i", "ni") %in%
                    names(prepared)))
  expect_false(any(c("tau", "mu", "se.mu", "se.tau") %in% names(prepared)))
  expect_equal(prepared$y0i, c(2.5, 3.5))
  expect_equal(prepared$yi, c(2, 3))
  expect_equal(prepared$sei, rep(sqrt(2)*sd(1:4)/2, 2))
  counts <- data.frame(group = LETTERS[1:2], a = c(10, 20), b = c(90, 80),
                        c = c(15, 25), d = c(85, 75))
  for(effect in c("logOR", "logRR", "RD")) {
    prepared <- prepare_ma(counts, effect = effect)
    expect_true(all(c("a", "b", "c", "d", "n1", "n2", "yi", "sei") %in%
                      names(prepared)))
    expected <- with(counts, switch(effect,
      logOR = log(a*d/(b*c)), logRR = log((a/(a+b))/(c/(c+d))),
      RD = a/(a+b) - c/(c+d)))
    expect_equal(prepared$yi, expected)
    expected_se <- with(counts, switch(effect,
      logOR = sqrt(1/a + 1/b + 1/c + 1/d),
      logRR = sqrt(1/a + 1/c - 1/(a+b) - 1/(c+d)),
      RD = sqrt((a/(a+b))*(b/(a+b))/(a+b) +
                (c/(c+d))*(d/(c+d))/(c+d))))
    expect_equal(prepared$sei, expected_se)
  }
})

test_that("canonical and legacy fits have identical priors and posterior summaries", {
  for(model in c("rubin", "vivs")) {
    canonical <- if(model == "rubin") summary_rubin else summary_mutau
    legacy <- if(model == "rubin") legacy_rubin else legacy_mutau
    new <- suppressWarnings(baggr(canonical, iter = 100, warmup = 50,
      chains = 1, seed = 235, refresh = 0, silent = TRUE, warn = FALSE))
    old <- suppressWarnings(baggr(legacy, iter = 100, warmup = 50,
      chains = 1, seed = 235, refresh = 0, silent = TRUE, warn = FALSE))
    expect_identical(new$model, model)
    expect_identical(new$formatted_prior, old$formatted_prior)
    expect_equal(treatment_effect(new, summary = TRUE),
                 treatment_effect(old, summary = TRUE), tolerance = 1e-12)
    expect_equal(pooling(new), pooling(old), tolerance = 1e-12)
    expect_s3_class(funnel_plot(new), "ggplot")
    expect_s3_class(forest_plot(new, show = "inputs"), "gforge_forestplot")
    custom <- canonical
    names(custom)[names(custom) == "yi"] <- "effect"
    names(custom)[names(custom) == "sei"] <- "effect_se"
    mapped <- suppressWarnings(baggr(custom, model = model,
      yi = "effect", sei = "effect_se", iter = 100, warmup = 50,
      chains = 1, seed = 235, refresh = 0, silent = TRUE, warn = FALSE))
    expect_identical(mapped$formatted_prior, new$formatted_prior)
    expect_equal(treatment_effect(mapped, summary = TRUE),
                 treatment_effect(new, summary = TRUE), tolerance = 1e-12)
  }
})

test_that("baggr consumes slab and stores canonical held-out data", {
  sampling_args <- NULL
  original_sampling <- rstan::sampling
  local_mocked_bindings(sampling = function(object, ...) {
    sampling_args <<- list(...)
    original_sampling(object, ...)
  }, .package = "rstan")
  df <- legacy_rubin
  names(df)[1] <- "study_label"
  messages <- capture_messages(bg <- suppressWarnings(baggr(
    df[-1, ], test_data = df[1, ], slab = "study_label", iter = 100,
    warmup = 50, chains = 1, seed = 235, refresh = 0, warn = FALSE)))
  expect_false(any(grepl("Using `slab`", messages, fixed = TRUE)))
  expect_false("slab" %in% names(sampling_args))
  expect_identical(group_names(bg), df$study_label[-1])
  expect_equal(bg$test_data$yi, df$tau[1])
  expect_equal(bg$test_data$sei, df$se[1])
  expect_error(baggr(df, slab = "absent"), "must name a column")
  df$group <- letters[1:6]
  messages <- capture_messages(bg <- suppressWarnings(baggr(
    df, group = "group", slab = "study_label", iter = 100,
    warmup = 50, chains = 1, seed = 235, refresh = 0, warn = FALSE)))
  expect_true(any(grepl("using `group`", messages, fixed = TRUE)))
  expect_identical(group_names(bg), df$group)
})

test_that("wide duplicate conventions are checked without changing priority", {
  both <- transform(summary_mutau, mu = y0i, se.mu = se0i,
                     tau = yi, se.tau = sei, y = yi, se = sei)
  expect_message(normalise_summary_columns(both), NA)
  both$mu[1] <- 99
  expect_error(normalise_summary_columns(both), "Conflicting summary columns")
  both$mu <- both$y0i
  both$se.tau[1] <- 99
  expect_error(normalise_summary_columns(both), "Conflicting summary columns")
  expect_message(normalise_summary_columns(both, yi = "yi", sei = "sei",
                                           y0i = "y0i", se0i = "se0i"), NA)
  test <- both[1, ]
  both$se.tau <- both$sei
  expect_error(normalise_summary_columns(both, test), "Conflicting summary columns")
})

test_that("explicit mappings remain authoritative in downstream plots", {
  df <- transform(summary_rubin, effect = yi + 1, effect_se = sei*2,
                   tau = yi - 1, se = sei/2, y = yi + 2)
  bg <- suppressWarnings(baggr(df, yi = "effect", sei = "effect_se",
    iter = 100, warmup = 50, chains = 1, seed = 235, refresh = 0,
    silent = TRUE, warn = FALSE))
  expect_equal(bg$data$yi, df$effect)
  expect_equal(funnel_plot(bg)$data$effect, df$effect)
  expect_equal(funnel_plot(bg)$data$se, df$effect_se)
  expect_s3_class(forest_plot(bg, show = "inputs"), "gforge_forestplot")
  expect_true(all(is.finite(pooling(bg))))

  df <- summary_rubin
  names(df)[1] <- "slab"
  messages <- capture_messages(bg <- suppressWarnings(baggr(df,
    iter = 100, warmup = 50, chains = 1, seed = 235, refresh = 0, warn = FALSE)))
  expect_false(any(grepl("Using `slab`", messages, fixed = TRUE)))
  expect_identical(group_names(bg), df$slab)
})
