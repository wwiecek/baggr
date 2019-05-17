context("Cross-validation methods")
library(baggr)

# assuming test are for now always executed on WW's machine or similar
# we could enable this:
# options(mc.cores = 4)

# tests ----------------------------------------------------------
test_that("LOO CV runs for all models", {
  #Rubin models:
  bgcv1 <- loocv(schools[1:4,], pooling = "full", iter = 100)
  bgcv2 <- loocv(schools[1:4,], pooling = "partial", iter = 100)
  expect_is(bgcv1, "baggr_cv")
  expect_is(bgcv2, "baggr_cv")
})

