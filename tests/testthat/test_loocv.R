context("Cross-validation methods")
library(baggr)


# tests ----------------------------------------------------------
test_that("LOO CV runs for all models", {
  #Rubin models:
  bgcv1 <- loocv(schools, pooling = "full", iter = 100)
  bgcv2 <- loocv(schools, pooling = "partial", iter = 100)
  expect_is(bgcv1, "baggr_cv")
  expect_is(bgcv2, "baggr_cv")
})

