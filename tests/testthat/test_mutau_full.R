library(baggr)

set.seed(1990)

test_that("prepare_ma works for synthetic mutau_full-style IPD", {
  df <- data.frame()

  K <- 10
  N <- c(rep(20, 5), rep(20, 5))
  bsl <- c(0,0,0,0,0,1,1,1,1,1)
  se <- rep(1, 10)
  tau <- c(1,1,1,1,1,0,0,0,0,0)

  for(i in 1:K){
    df <- rbind(df,
                data.frame(group = i, outcome = bsl[i] + rnorm(N[i], 0, se[i]), treatment = 0),
                data.frame(group = i, outcome = bsl[i] + tau[i] + rnorm(N[i], 0, se[i]), treatment = 1))
  }

  ma <- prepare_ma(df)
  expect_s3_class(ma, "data.frame")
  expect_true(all(c("group", "mu", "tau", "se.mu", "se.tau", "n.mu", "n.tau") %in% names(ma)))
  expect_equal(nrow(ma), K)
})
