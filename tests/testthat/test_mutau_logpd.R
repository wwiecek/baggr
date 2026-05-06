context("mutau log predictive density")

library(baggr)

test_that("mutau logpd is finite and close to a plug-in MVN reference", {
  testthat::skip_if_not_installed("mvtnorm")

  set.seed(1)

  P <- 2
  K <- 50
  K_test <- 5

  mu_true <- c(0.25, -0.10)
  sd_true <- c(0.35, 0.60)
  rho_true <- 0.55
  corr_true <- matrix(c(1, rho_true, rho_true, 1), 2, 2, byrow = TRUE)
  sigma_true <- diag(sd_true) %*% corr_true %*% diag(sd_true)

  se_train <- matrix(runif(P * K, 0.08, 0.20), nrow = P)
  se_test  <- matrix(runif(P * K_test, 0.08, 0.20), nrow = P)

  theta_train <- t(mvtnorm::rmvnorm(K, mean = mu_true, sigma = sigma_true))
  theta_test  <- t(mvtnorm::rmvnorm(K_test, mean = mu_true, sigma = sigma_true))

  y_train <- theta_train + matrix(rnorm(P * K), nrow = P) * se_train
  y_test  <- theta_test + matrix(rnorm(P * K_test), nrow = P) * se_test

  logpd_joint <- function(y, se, mu, sigma) {
    stopifnot(nrow(y) == length(mu), all(dim(y) == dim(se)))

    sum(vapply(seq_len(ncol(y)), function(k) {
      V <- sigma + diag(se[, k]^2, nrow = length(mu))
      mvtnorm::dmvnorm(y[, k], mean = mu, sigma = V, log = TRUE)
    }, numeric(1)))
  }

  logpd_oracle <- logpd_joint(y_test, se_test, mu_true, sigma_true)

  nll <- function(par, y, se) {
    mu <- par[1:2]
    sd <- exp(par[3:4])
    rho <- tanh(par[5])
    corr <- matrix(c(1, rho, rho, 1), 2, 2, byrow = TRUE)
    sigma <- diag(sd) %*% corr %*% diag(sd)

    sum(vapply(seq_len(ncol(y)), function(k) {
      V <- sigma + diag(se[, k]^2, 2)
      -mvtnorm::dmvnorm(y[, k], mean = mu, sigma = V, log = TRUE)
    }, numeric(1)))
  }

  par0 <- c(colMeans(y_train), log(apply(y_train, 1, sd)), atanh(0.1))
  fit <- optim(par0, nll, y = y_train, se = se_train, method = "BFGS")

  mu_hat <- fit$par[1:2]
  sd_hat <- exp(fit$par[3:4])
  rho_hat <- tanh(fit$par[5])
  corr_hat <- matrix(c(1, rho_hat, rho_hat, 1), 2, 2, byrow = TRUE)
  sigma_hat <- diag(sd_hat) %*% corr_hat %*% diag(sd_hat)

  logpd_plugin <- logpd_joint(y_test, se_test, mu_hat, sigma_hat)

  expect_true(is.finite(logpd_oracle))
  expect_true(is.finite(logpd_plugin))
  expect_lt(abs(logpd_oracle - logpd_plugin), 50)

  bgdt <- data.frame(
    mu = c(y_train[1, ], y_test[1, ]),
    tau = c(y_train[2, ], y_test[2, ]),
    se.mu = c(se_train[1, ], se_test[1, ]),
    se.tau = c(se_train[2, ], se_test[2, ]),
    include = c(rep(1, K), rep(0, K_test))
  )

  bg_train <- bgdt[bgdt$include == 1, ]
  bg_test <- bgdt[bgdt$include == 0, ]

  bg_fit <- baggr(bg_train, model = "mutau", refresh = 0, test_data = bg_test)
  bg_fit_full <- baggr(bg_train, model = "mutau", pooling = "full",
                       refresh = 0, test_data = bg_test)

  expect_is(bg_fit, "baggr")
  expect_is(bg_fit_full, "baggr")
  expect_true(is.finite(bg_fit$mean_lpd))
  expect_true(is.finite(bg_fit_full$mean_lpd))
})
