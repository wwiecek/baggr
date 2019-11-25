context("baggr() calls with full model")
library(baggr)
set.seed(1990)

ms <- microcredit_simplified[sample(nrow(microcredit_simplified), 1000), ]
ms$outcome <- log(ms$consumerdurables + 1)
bg_n <- expect_warning(baggr(ms, pooling = "none", iter = 200, refresh=0))
bg_p <- expect_warning(baggr(ms, pooling = "partial", iter = 200, refresh=0))
bg_f <- expect_warning(baggr(ms, pooling = "full", iter = 200, refresh=0))


test_that("Different pooling methods work for the full model", {
  expect_is(bg_n, "baggr")
  expect_is(bg_p, "baggr")
  expect_is(bg_f, "baggr")
})

test_that("Basic operations on full data model", {
  expect_error(baggr(ms, rubbish = 41))
  expect_is(pooling(bg_p)[,,1], "matrix")
  expect_is(plot(bg_p), "gg")
  expect_is(effect_plot(bg_p), "vpPath")
  expect_is(forest_plot(bg_p), "vpPath")
  bgc <- try(baggr_compare(bg_n, bg_p, bg_f))
  expect_is(bgc, "gg")
  expect_error(loocv(ms))
})

test_that("Full model crashes with nonsense inputs", {
  expect_error(baggr(ms, outcome = 2), "Arguments")
  expect_error(baggr(ms, group = 2), "Arguments")
  expect_error(baggr(ms, treatment = 2), "Arguments")
  expect_error(baggr(ms, treatment = "wrong"), "no column")
  ms2 <- ms; ms2$treatment <- as.character(ms$treatment)
  expect_error(baggr(ms2), "has to be numeric")
  ms2 <- ms; ms2$outcome <- as.character(ms$outcome)
  expect_error(baggr(ms2), "has to be numeric")

})
