context("coverage for selection and validation helpers")

library(baggr)

test_that("selection() validates inputs before extracting draws", {
  expect_error(selection(list()), "Object of class 'baggr' required")

  fake_bg <- structure(list(model = "rubin_full", inputs = list(M = 2)), class = "baggr")
  expect_error(selection(fake_bg), "only available for selection models")

  fake_bg$model <- "rubin"
  fake_bg$inputs$M <- 0
  expect_error(selection(fake_bg), "does not define any selection intervals")
})

test_that("numeric selection input maps to default symmetric selection data", {
  inp <- convert_inputs(schools, "rubin", selection = c(1.96, 2.58))

  expect_equal(inp$M, 2L)
  expect_equal(as.numeric(inp$c), c(1.96, 2.58))
  expect_equal(inp$symmetric, 1L)
  expect_equal(as.integer(inp$possible_selection), rep(1L, nrow(schools)))
})

test_that("list selection input controls z, symmetry and possible flags", {
  possible <- c(1, 0, rep(1, nrow(schools) - 2))
  inp <- convert_inputs(
    schools, "rubin",
    selection = list(z = c(1.64, 1.96), symmetrical = FALSE, possible = possible)
  )

  expect_equal(inp$M, 2L)
  expect_equal(as.numeric(inp$c), c(1.64, 1.96))
  expect_equal(inp$symmetric, 0L)
  expect_equal(as.integer(inp$possible_selection), possible)
})

test_that("numeric shorthand is equivalent to explicit all-possible selection", {
  shorthand <- convert_inputs(schools, "rubin", selection = c(1.96, 2.58))
  explicit <- convert_inputs(
    schools, "rubin",
    selection = list(z = c(1.96, 2.58), symmetrical = TRUE,
                     possible = rep(1, nrow(schools)))
  )

  expect_identical(shorthand$M, explicit$M)
  expect_equal(shorthand$c, explicit$c)
  expect_identical(shorthand$symmetric, explicit$symmetric)
  expect_identical(shorthand$possible_selection, explicit$possible_selection)
})

test_that("selection input is validated", {
  expect_error(
    convert_inputs(schools, "rubin", selection = list(z = 1.96, possible = rep(1, 8))),
    "named elements z, symmetrical and possible"
  )
  expect_error(
    convert_inputs(schools, "rubin",
                   selection = list(z = c(1.96, 1.64), symmetrical = TRUE,
                                    possible = rep(1, 8))),
    "strictly increasing"
  )
  expect_error(
    convert_inputs(schools, "rubin",
                   selection = list(z = 0, symmetrical = TRUE, possible = rep(1, 8))),
    "positive finite"
  )
  expect_error(
    convert_inputs(schools, "rubin",
                   selection = list(z = 1.96, symmetrical = 1, possible = rep(1, 8))),
    "TRUE or FALSE"
  )
  expect_error(
    convert_inputs(schools, "rubin",
                   selection = list(z = 1.96, symmetrical = TRUE, possible = rep(1, 7))),
    "one value per study"
  )
  expect_error(
    convert_inputs(schools, "rubin",
                   selection = list(z = 1.96, symmetrical = TRUE,
                                    possible = c(1, 2, rep(1, 6)))),
    "0/1 or logical"
  )
  expect_error(
    convert_inputs(schools, "mutau", selection = 1.96),
    "only for model = 'rubin'"
  )
})

test_that("possible selection flag is used by Stan selection likelihood", {
  stan_file <- system.file("stan/functions/selection.stan", package = "baggr")
  if(!nzchar(stan_file))
    stan_file <- test_path("../../inst/stan/functions/selection.stan")
  stan_code <- paste(readLines(stan_file), collapse = "\n")

  expect_match(stan_code, "possible_selection == 0", fixed = TRUE)
  expect_match(stan_code, "if (symmetric)\n      z = fabs(z);", fixed = TRUE)
  expect_match(stan_code, "a = symmetric ? 0.0 : negative_infinity();",
               fixed = TRUE)
  expect_match(stan_code, "w_obs = (j <= M) ? omega[j] : 1.0;", fixed = TRUE)
})

test_that("baggr() accepts list selection input", {
  fit <- suppressMessages(suppressWarnings(
    baggr(
      schools, selection = list(z = 1.96, symmetrical = TRUE, possible = rep(1, 8)),
      iter = 50, chains = 1, refresh = 0, seed = 123
    )
  ))

  expect_s3_class(fit, "baggr")
  expect_equal(fit$inputs$M, 1L)
  expect_equal(as.numeric(fit$inputs$c), 1.96)
  expect_equal(fit$inputs$symmetric, 1L)
  expect_equal(as.integer(fit$inputs$possible_selection), rep(1L, 8))

  omega_draws <- selection(fit, summary = FALSE)
  omega_summary <- selection(fit)

  expect_equal(ncol(omega_draws), 1L)
  expect_equal(nrow(omega_summary), 1L)
  expect_true(all(c("mean", "sd") %in% colnames(omega_summary)))
})

test_that("selection can be disabled study-wise in prior predictive fits", {
  fit <- suppressMessages(suppressWarnings(
    baggr(
      schools,
      selection = list(z = 1.96, symmetrical = TRUE, possible = rep(0, 8)),
      ppd = TRUE, iter = 50, chains = 1, refresh = 0, seed = 123
    )
  ))

  expect_s3_class(fit, "baggr")
  expect_equal(fit$inputs$K, 0L)
  expect_equal(length(fit$inputs$possible_selection), 0L)
})

test_that("check_columns_ipd() validates columns and argument types", {
  dat <- data.frame(
    outcome = c(1.2, 0.8, 1.5),
    group = c("a", "a", "b"),
    treatment = c(0, 1, 0),
    stringsAsFactors = FALSE
  )

  expect_error(
    baggr:::check_columns_ipd(dat, outcome = 1, group = "group", treatment = "treatment"),
    'must be of type "character"'
  )
  expect_error(
    baggr:::check_columns_ipd(dat, outcome = "missing", group = "group", treatment = "treatment"),
    "There's no column 'missing'"
  )

  dat_bad <- dat
  dat_bad$outcome <- as.character(dat_bad$outcome)
  expect_error(
    baggr:::check_columns_ipd(dat_bad, outcome = "outcome", group = "group", treatment = "treatment"),
    "Outcome variable in baggr has to be numeric"
  )

  dat_bad <- dat
  dat_bad$treatment <- c("a", "b", "a")
  expect_error(
    baggr:::check_columns_ipd(dat_bad, outcome = "outcome", group = "group", treatment = "treatment"),
    "Treatment variable in baggr has to be numeric or a factor"
  )

  dat_bad <- dat
  dat_bad$group[2] <- NA
  expect_error(
    baggr:::check_columns_ipd(dat_bad, outcome = "outcome", group = "group", treatment = "treatment"),
    "Some of group values are NA"
  )

  dat_bad <- dat
  dat_bad$treatment <- c(0, 2, 0)
  expect_error(
    baggr:::check_columns_ipd(dat_bad, outcome = "outcome", group = "group", treatment = "treatment", trt_binary = TRUE),
    "Treatment column has to have values 0 or 1"
  )

  dat_ok <- dat
  expect_null(
    baggr:::check_columns_ipd(dat_ok, outcome = "outcome", group = "group", treatment = "treatment", trt_binary = FALSE)
  )
})

test_that("loo_compare and print methods work with synthetic baggr_cv objects", {
  mk_cv <- function(pointwise) {
    structure(list(
      df = data.frame(lpd = pointwise),
      pointwise = pointwise,
      elpd = sum(pointwise),
      se = 0.1,
      looic = -2 * sum(pointwise),
      K = length(pointwise)
    ), class = "baggr_cv")
  }

  cv1 <- mk_cv(c(-1.0, -1.1, -0.9))
  cv2 <- mk_cv(c(-1.2, -1.0, -0.8))

  comp <- loo_compare(Reference = cv1, Alternative = cv2)
  expect_s3_class(comp, "compare_baggr_cv")
  expect_equal(dim(comp), c(1, 2))

  unnamed_comp <- loo_compare(cv1, cv2)
  expect_true(any(grepl("Model 1 - Model 2", rownames(unnamed_comp))))

  expect_error(loo_compare(cv1), "requires at least two models")
  expect_error(loo_compare(A = cv1, A = cv2), "unique model names")

  cv3 <- mk_cv(c(-1.0, -1.1))
  expect_error(loo_compare(cv1, cv3), "same number of data points")

  expect_type(testthat::capture_output(print(comp)), "character")
  expect_type(testthat::capture_output(print(cv1)), "character")
})
