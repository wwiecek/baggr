library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

# pooled, with equal SE's!
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])


bg_onerow_p <- NULL
bg_onerow_n <- NULL
bg_onerow_f <- NULL
bg_onerow_binary <- NULL

setup({
  bg_onerow_p <<- expect_warning(baggr(df_pooled[1,], pooling = "partial", group = "state",
                                       iter = 200, chains = 2, refresh = 0,
                                       show_messages = F, prior_hypersd = normal(0,1)))
  bg_onerow_n <<- expect_warning(baggr(df_pooled[1,], pooling = "none", group = "state",
                                       iter = 200, chains = 2, refresh = 0,
                                       show_messages = F, prior_hypersd = normal(0,1)))
  bg_onerow_f <<- expect_warning(baggr(df_pooled[1,], pooling = "full", group = "state",
                                       iter = 200, chains = 2, refresh = 0,
                                       show_messages = F, prior_hypersd = normal(0,1)))
  bg_onerow_binary <<- expect_warning(
    baggr(yusuf[1,], model = "logit",
          prior_hypersd = normal(0,1),
          prior_control_sd = normal(0, 1),
          pooling_control = "partial",
          iter = 200, chains = 2, refresh = 0)
  )
})

test_that("The thing runs", {

  expect_s3_class(bg_onerow_f, "baggr")
  expect_s3_class(bg_onerow_p, "baggr")
  expect_s3_class(bg_onerow_n, "baggr")
  expect_s3_class(bg_onerow_binary, "baggr")

  expect_error(baggr(df_pooled[1,], pooling = "partial", group = "state",
                     iter = 200, chains = 2, refresh = 0,
                     show_messages = F), "specify hyper-SD prior")
  expect_error(baggr(yusuf[1,], model = "logit",
                     prior_hypersd = normal(0,1),
                     # prior_control_sd = normal(0, 1),
                     pooling_control = "partial"),
               "You must specify SD in baseline rates")
  expect_error(plot(bg_onerow_p), "You can only plot meta-analyses with more than 1 group.")
  baggr_compare(bg_onerow_p, bg_onerow_f)

  gg1 <- baggr_compare(bg_onerow_p, bg_onerow_f, compare = "effects") %>% plot
  gg2 <- baggr_compare(bg_onerow_p, bg_onerow_f, compare = "hyperpars") %>% plot
  expect_s3_class(gg1, "gg")
  expect_s3_class(gg2, "gg")
})




