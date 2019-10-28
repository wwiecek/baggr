context("Prior and posterior predictions")
library(baggr)
library(testthat)
set.seed(11241)

df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])
df_mutau <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                       "se.tau" = rep(1, 8),
                       "mu" = rnorm(8),
                       "se.mu" = rep(1, 8),
                       "state" = datasets::state.name[1:8])

bg_ppd <- expect_warning(baggr(schools, iter = 200, refresh = 0, ppd = TRUE))

test_that("Basic ppd usage", {
  expect_error(baggr(schools, pooling = "none", ppd = TRUE))
  expect_is(bg_ppd, "baggr")
  expect_true(attr(bg_ppd, "ppd"))
  capture_output(bg_ppd) #printing
  expect_is(treatment_effect(bg_ppd), "list") #extracting TE

})
