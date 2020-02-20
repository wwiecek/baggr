context("Specifying priors for baggr models")
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

test_that("Wrong prior specifications crash baggr()", {
  expect_error(baggr(df_mutau, prior_hypermean = 2))
  expect_error(baggr(df_mutau, prior_hypermean = normal(2, -100)))
  expect_error(baggr(df_mutau, prior_hypermean = lkj(4)))
  expect_error(baggr(df_mutau, prior_hypermean = list(dist = "not_a_dist", 0, 10)))
  expect_error(baggr(df_mutau, prior_hypermean = list(dist = "not_a_dist")))
  expect_error(baggr(df_mutau, prior_hypermean = list(dist = "normal", a = 5, b = 6)))
})

test_that("Prior specification via different arguments", {
  custom_prior <- list(hypermean = normal(0, 10), hypersd = uniform(0, 20))
  bg_prior1 <- expect_warning(baggr(df_pooled, "rubin",
                                    iter = 200, chains = 2, refresh = 0, seed = 1990,
                                    prior_hypermean = normal(0, 2),
                                    prior = custom_prior)) #custom_prior OVERRIDES hypermean!
  bg_prior2 <- expect_warning(baggr(df_pooled, "rubin",
                                    iter = 200, chains = 2, refresh = 0, seed = 1990,
                                    prior_hypermean = normal(0,10),
                                    prior_hypersd = uniform(0,20)))
  bg_prior3 <- expect_warning(baggr(df_pooled, "rubin",
                                    iter = 200, chains = 2, refresh = 0, seed = 1990,
                                    formatted_prior = bg_prior2$formatted_prior))
  # Same result (given the same seed):
  te3 <- sort(treatment_effect(bg_prior3)[[1]])
  te2 <- sort(treatment_effect(bg_prior2)[[1]])
  te1 <- sort(treatment_effect(bg_prior1)[[1]])
  expect_identical(te1, te2)
  expect_identical(te3, te2)

  # Wrong names in the list
  expect_error(baggr(df_pooled, prior = list(hypermeann = normal(0,5))),
               "Prior argument")
  expect_error(baggr(df_pooled, prior = list(hypermeann = normal(0,5),
                                             hypermean = normal(0,5))),
               "Prior argument")
})

test_that("All possible prior dist's work", {
  expect_is(normal(0, 10), "list")
  expect_is(cauchy(0, 10), "list")
  expect_is(uniform(0, 10), "list")
  expect_is(multinormal(c(0,0), diag(2)), "list")
  expect_is(lkj(5), "list")

  expect_error(multinormal(0, 10))
  expect_error(normal(c(0,0), diag(2)))
  expect_error(cauchy(0, 5, 8))

  expect_error(normal("0", -1))
  expect_error(normal(0, -1))
  expect_error(uniform(0, -1))
  expect_error(cauchy(0, -1))
  expect_error(multinormal(c(0,0), matrix(c(-1,0,-1,1),2,2)), "positive")
  expect_error(multinormal(c(0,0,0), diag(2)), "dimensions")
})

test_that("Different priors for mutau model", {
  bg1 <- expect_warning(baggr(df_mutau, prior_hypermean = normal(0, 5),
                              iter = 200, chains = 2, refresh = 0))
  bg2 <- expect_warning(baggr(df_mutau, prior_hypercor = lkj(4),
                              iter = 200, chains = 2, refresh = 0))
  bg3 <- expect_warning(baggr(df_mutau, prior_hypersd = normal(0, 5),
                              iter = 200, chains = 2, refresh = 0))
  expect_error(baggr(df_mutau, prior_hypermean = multinormal(c(0,0,0), diag(3))))
  expect_error(baggr(df_mutau, prior_hypercor  = multinormal(c(0,0), diag(2))))
  expect_is(bg1, "baggr")
  expect_is(bg2, "baggr")
  expect_is(bg3, "baggr")
})


test_that("Prior vs posterior and PPD comparisons work", {
  # Invalid comparison for prior vs posterior:
  expect_error(baggr_compare(schools, ppd = TRUE, what = "prior"))

  # Typical PPD objects:
  bg_ppd1 <- expect_warning(baggr(schools, ppd = T, refresh = 0, iter = 200))
  bg_ppd2 <- expect_warning(baggr(schools, ppd = T, prior_hypermean = normal(0,10), refresh = 0, iter = 200))
  expect_is(bg_ppd1, "baggr")
  expect_is(bg_ppd2, "baggr")
  # Regular comparison (don't have to say compare = "groups")
  bgc <- baggr_compare(bg_ppd1, bg_ppd2)
  expect_is(bgc, "baggr_compare")

  # Prior vs posterior
  bgc2 <- expect_warning(baggr_compare(schools, what = "prior", refresh = 0, iter = 200))
  expect_is(bgc2, "baggr_compare")

  # Effect plot of PPD:
  gg <- effect_plot(bg_ppd1)
  expect_identical(gg$labels$title, "Prior distribution for pooled treatment effect")


})
