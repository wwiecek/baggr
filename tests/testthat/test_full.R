context("baggr() calls with full model")
library(baggr)
set.seed(1990)

ms <- microcredit_simplified[sample(nrow(microcredit_simplified), 1000), ]
ms$outcome <- log(ms$consumerdurables + 1)
bg_n <- expect_warning(baggr(ms, pooling = "none", iter = 200))
bg_p <- expect_warning(baggr(ms, pooling = "partial", iter = 200))
bg_f <- expect_warning(baggr(ms, pooling = "full", iter = 200))


test_that("Different pooling methods work for the full model", {
  expect_is(bg_n, "baggr")
  expect_is(bg_p, "baggr")
  expect_is(bg_f, "baggr")
})

test_that("Basic operations on full data model", {
  expect_error(baggr(ms, rubbish = 41))
  expect_is(pooling(bg_p)[,,1], "matrix")
  expect_is(plot(bg_p), "gg")
  bgc <- baggr_compare(bg_n, bg_p, bg_f)
  expect_is(bgc, "gg")
  expect_error(loocv(ms))
})

