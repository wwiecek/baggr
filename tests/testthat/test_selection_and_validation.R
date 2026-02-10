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
