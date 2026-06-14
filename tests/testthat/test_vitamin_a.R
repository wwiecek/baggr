context("vitamin A data and publication bias vignette helpers")

library(baggr)

test_that("vitamin_a is analysis-ready summary data", {
  expect_s3_class(vitamin_a, "data.frame")
  expect_equal(names(vitamin_a), c("group", "tau", "se"))
  expect_equal(nrow(vitamin_a), 18)
  expect_false("Lin 2008" %in% vitamin_a$group)
  expect_true(all(is.finite(vitamin_a$tau)))
  expect_true(all(is.finite(vitamin_a$se)))
  expect_true(all(vitamin_a$se > 0))
})

test_that("funnel() delegates to funnel_plot()", {
  fit <- suppressMessages(suppressWarnings(
    baggr(vitamin_a, iter = 50, chains = 1, refresh = 0, seed = 123,
          effect = "log risk ratio")
  ))

  expect_s3_class(funnel(fit), "ggplot")
  expect_s3_class(funnel_plot(fit), "ggplot")
})
