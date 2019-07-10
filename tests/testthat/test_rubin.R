context("baggr() calls with Rubin model")
library(baggr)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

# pooled, with equal SE's!
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])

# tests ----------------------------------------------------------
test_that("Error messages for wrong inputs are in place", {
  # model doesn't exist
  expect_error(baggr(df_pooled, "made_up_model"), "Unrecognised model")

  # NA or NULL inputs
  df_na <- df_pooled; df_na$tau[1] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_pooled; df_na$se[2] <- NA
  expect_error(baggr(df_na),"NA values")
  df_na <- df_pooled; df_na$se <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_pooled; df_na$tau <- NULL
  expect_error(baggr(df_na),"no column")
  df_na <- df_pooled; df_na$tau <- as.character(df_na$tau)
  expect_error(baggr(df_na),"are not numeric")

  # test_that("Converting inputs works correctly") more explicitly
  expect_identical(names(convert_inputs(df_pooled, "rubin")),
                   c("K", "tau_hat_k", "se_tau_k", "K_test", "test_tau_hat_k", "test_se_k"))

  expect_warning(baggr(df_pooled, group = "state1000", iter = 50), "No labels will be added.")

})



bg5_n <- baggr(df_pooled, "rubin", pooling = "none", group = "state",
               iter = 200, chains = 2)
bg5_p <- baggr(df_pooled, "rubin", pooling = "partial", group = "state",
               iter = 200, chains = 2)
bg5_f <- baggr(df_pooled, "rubin", pooling = "full", group = "state",
               iter = 200, chains = 2)

test_that("Extra args to Stan passed via ... work well", {
  expect_equal(nrow(as.matrix(bg5_p$fit)), 200) #right dimension means right iter
  expect_error(baggr(df_pooled, rubbish = 41))
})

test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, 8)
  expect_equal(bg5_p$effects, "mean")
  expect_equal(bg5_p$model, "rubin")
  expect_is(bg5_p$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_identical(bg5_n$data, df_pooled)
  expect_identical(bg5_p$data, df_pooled)
  expect_identical(bg5_f$data, df_pooled)
})

test_that("Different pooling methods work for Rubin model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
})

test_that("Pooling metrics", {
  # all pooling metric are the same as SE's are the same
  expect_equal(length(unique(bg5_p$pooling_metric[1,,1])), 1) #expect_length()
  expect_equal(length(unique(bg5_p$pooling_metric[2,,1])), 1)
  expect_equal(length(unique(bg5_p$pooling_metric[3,,1])), 1)
  # all pooling stats are 0 if no pooling
  expect_equal(unique(as.numeric(bg5_n$pooling_metric)), 0)
  # full pooling means 1's everywhere
  expect_equal(unique(as.numeric(bg5_f$pooling_metric)), 1)

  pp <- pooling(bg5_p)
  expect_is(pp, "array")
  expect_gt(min(pp), 0)
  expect_lt(max(pp), 1)
  # since all SEs are the same, pooling should be the same for all sites
  print(pp)
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 200 iter
  expect_equal(length(unique(pp[2,,1])), 1)
  expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)
})


test_that("Calculation of effects works", {
  expect_is(study_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")

  expect_identical(dim(study_effects(bg5_n)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(study_effects(bg5_p)), as.integer(c(200, 8 , 1)))
  expect_identical(dim(study_effects(bg5_f)), as.integer(c(200, 8 , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting works", {
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p), "gg")
  expect_is(plot(bg5_f), "gg")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")
})


# 8 schools test -----

bg_s <- baggr(schools)

test_that("The default 8 schools result is close to the result in BDA", {
  expect_equal(mean(treatment_effect(bg_s)$tau), 8, tolerance = .25)
  expect_equal(mean(treatment_effect(bg_s)$sigma_tau), 7, tolerance = .25)

})
