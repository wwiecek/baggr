context("Additional coverage for internal helper functions")

library(baggr)

test_that("detect_input_type() identifies supported formats", {
  df_narrow <- data.frame(group = letters[1:3], treatment = c(0, 1, 1),
                          tau = c(0.1, 0.2, 0.3), se = c(0.01, 0.02, 0.03))
  expect_identical(baggr:::detect_input_type(df_narrow), "pool_narrow")

  df_noctrl <- df_narrow[, c("group", "tau", "se")]
  expect_identical(baggr:::detect_input_type(df_noctrl), "pool_noctrl_narrow")

  df_binary <- data.frame(group = LETTERS[1:2], a = c(3, 4), c = c(1, 2),
                          n1 = c(10, 10), n2 = c(12, 14))
  expect_identical(baggr:::detect_input_type(df_binary), "pool_binary")

  df_wide <- data.frame(group = LETTERS[1:3],
                        tau = c(0.1, 0.2, 0.3),
                        mu = c(1, 2, 3),
                        se.mu = c(0.1, 0.1, 0.1),
                        se.tau = c(0.2, 0.2, 0.2))
  expect_identical(baggr:::detect_input_type(df_wide), "pool_wide")

  df_ipd_binary <- data.frame(group = rep(LETTERS[1:2], each = 2),
                              treatment = c(0, 1, 0, 1),
                              outcome = c(0, 1, 0, 1))
  expect_identical(baggr:::detect_input_type(df_ipd_binary), "individual_binary")

  df_ipd <- transform(df_ipd_binary, outcome = c(0.1, 0.2, 0.3, 0.4))
  expect_identical(baggr:::detect_input_type(df_ipd), "individual")

  expect_identical(baggr:::detect_input_type(data.frame(x = 1:3)), "unknown")
  expect_error(baggr:::detect_input_type(list(x = 1:3)), "not a data.frame")
})

test_that("column checks and binary checks behave as expected", {
  expect_identical(baggr:::check_columns_binary(data.frame(a = 1, c = 1, n1 = 2, n2 = 2), stop = FALSE), 1)
  expect_identical(baggr:::check_columns_binary(data.frame(a = 1, c = 1), stop = FALSE), 0)
  expect_error(baggr:::check_columns_binary(data.frame(a = 1, c = 1)), "Binary data")

  expect_identical(baggr:::is.binary(c(0, 1, 1, 0)), 1)
  expect_identical(baggr:::is.binary(c(0, 0, 0), both = TRUE), 0)
  expect_identical(baggr:::is.binary(c(0, 1, NA), both = TRUE), 1)
  expect_identical(baggr:::is.binary(c(0, 2, 1)), 0)
})

test_that("find_group_column() defaults to first character/factor column", {
  dat <- data.frame(id = 1:3, study = letters[1:3], y = rnorm(3))
  expect_identical(baggr:::find_group_column(dat, "study"), "study")
  expect_message(
    expect_identical(baggr:::find_group_column(dat, "group"), "id"),
    "No grouping column found"
  )
  expect_error(baggr:::find_group_column(list(x = 1), "group"), "not a data.frame")
})

test_that("prior helper utilities cover key branches", {
  expect_error(baggr:::check_scalar(c(1, 2)), "scalar")
  expect_error(baggr:::check_scalar("a"), "numeric")

  expect_identical(baggr:::print_dist(normal(0, 2)), "normal(0, 2^2)")
  expect_identical(baggr:::print_dist(multinormal(c(0, 1), diag(2))), "multinormal(...)")

  base <- list()
  prior <- normal(0, 1)
  out <- baggr:::set_prior_val(base, "hypermean", prior, p = 2)
  expect_identical(out$hypermean_fam, rep(out$hypermean_fam[1], 2))
  expect_equal(dim(out$hypermean_val), c(2, 3))

  out_arr <- baggr:::set_prior_val(base, "hypermean", prior, to_array = TRUE)
  expect_equal(dim(out_arr$hypermean_fam), 1)
  expect_equal(dim(out_arr$hypermean_val), c(1, 3))

  expect_error(
    baggr:::set_prior_val(base, "hypercor", lkj(1), p = 2),
    "can't be 'replicated'"
  )
})
