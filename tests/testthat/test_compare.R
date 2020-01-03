testthat::context("test calls to baggr_compare()")

library(baggr)
set.seed(1990)

# pre-set models
schools_pooled <- baggr(schools,
                        model = "rubin",
                        pooling = "full")

schools_partial_pooled <- baggr(schools,
                                model = "rubin",
                                pooling = "partial")

preset_comparison <- baggr_compare(schools_pooled,
                                   schools_partial_pooled)

params_comparison <- baggr_compare(schools,
                                   model = 'rubin',
                                   prior_hypermean = normal(0, 3),
                                   prior_hypersd = normal(0,2),
                                   prior_hypercor = lkj(2),
                                   what = "pooling")

microcredit_summary_data <- prepare_ma(microcredit_simplified,
                                       outcome = "consumerdurables")

mutau_comparison <-
  baggr_compare(data = microcredit_summary_data,
                model = "mutau",
      what = "prior", prior_hypercor = lkj(1),
      prior_hypersd = normal(0,10),
      prior_hypermean = multinormal(c(0,0),matrix(c(10,3,3,10),2,2)))



test_that("baggr_compare works", {
  expect_invisible(baggr:::make_silent(tmp <- 1:10))
  expect_is(preset_comparison, "baggr_compare")
  expect_is(params_comparison, "baggr_compare")
  expect_is(mutau_comparison, "baggr_compare")
  expect_is(testthat::capture_output(print(params_comparison)), "character")
  expect_gt(length(preset_comparison), 0)
})

rm(tmp)

# preset comparison
pooling_plt <- plot(preset_comparison)

# prior comparison
prior_plot <- plot(mutau_comparison)

# simple tests that show plots exist
test_that("plotting works", {
  expect_is(pooling_plt, "plot_list")
  expect_is(pooling_plt[[1]], "ggplot2")
  expect_is(prior_plot, "ggplot")
})

