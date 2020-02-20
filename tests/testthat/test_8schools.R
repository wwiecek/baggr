# tests of 8 schools
context("test 8 schools model in baggr")

# brms testing script
# library(brms)
library(baggr)

set.seed(1999)


# effect estimates from brms
# ran for 100,000 iterations
brms_fixef <- 7.88433

brms_ranef <- c(`School A` = 3.50503394056709, `School B` = -0.0197197298800334,
                `School C` = -1.79659389944439, `School D` = -0.251061954670167,
                `School E` = -2.7779755236811, `School F` = -1.76165581026491,
                `School G` = 2.76042774120207, `School H` = 0.535207950933801
)

brms_totalef <- brms_ranef + brms_fixef

# baggr estimates, run each time
baggr_fit <- baggr(schools, control = list(adapt_delta = 0.99),
                   iter = 5000, refresh = 0)
baggr_groupef <- group_effects(baggr_fit, summary = T)[]
baggr_fixef   <- mean(baggr::treatment_effect(baggr_fit)$tau)

# tolerance for differences
tol <- 0.25

test_that(desc = "baggr and brms are at least close for schools estimates", {
  # estimates within tolerance of each other
  for(i in 1:nrow(baggr_groupef)){

    expect_lt(abs(brms_totalef[i] - baggr_groupef[i,"mean",]),
              tol)
  }

  expect_lt(abs(baggr_fixef - brms_fixef), tol)

})
