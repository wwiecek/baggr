context("baggr helper functions")
library(baggr)
set.seed(1990)

test_that("prepare_ma()", {
  expect_error(prepare_ma(schools), "individual-level")
  expect_error(prepare_ma(microcredit_simplified), "no column")
  expect_warning(prepare_ma(microcredit, outcome = "consumption"))
  expect_warning(prepare_ma(microcredit[!is.na(microcredit$treatment),],
                          outcome = "consumption"))
  expect_error(prepare_ma(microcredit_simplified, outcome = "consumerdurables",
                          effect = "logRR"), "not binary")
  expect_error(prepare_ma(microcredit_simplified, outcome = "consumerdurables",
                          effect = "logOR"), "not binary")

  pm <- prepare_ma(microcredit_simplified, outcome = "consumerdurables")
  expect_is(pm, "data.frame")
  expect_equal(dim(pm), c(4,5))
  expect_identical(names(pm), c("group", "mu", "tau", "se.mu", "se.tau"))

  pm <- prepare_ma(microcredit_simplified, outcome = "consumerdurables", summarise = FALSE)
  expect_identical(dim(pm), dim(microcredit_simplified))

  mc2 <- microcredit_simplified
  names(mc2)[1] <- "study"
  expect_error(prepare_ma(mc2, outcome = "consumerdurables"), "no column")
  expect_is(prepare_ma(mc2, group = "study", outcome = "consumerdurables"), "data.frame")

  # prepare_ma for binary data
  df_pat2 <- data.frame(treatment = rbinom(900, 1, .5),
                        group = rep(paste("Trial", LETTERS[1:10]), each = 90)) %>%
    mutate(outcome = ifelse(treatment, rbinom(900, 1, .3), rbinom(900, 1, .15)))
  expect_is(prepare_ma(df_pat2, effect = "logOR"), "data.frame")
  expect_is(prepare_ma(df_pat2, effect = "logRR"), "data.frame")

})

test_that("convert_inputs()", {
  # Rubin model
  expect_is(convert_inputs(schools, "rubin"), "list")
  expect_error(convert_inputs(schools, "mutau"))
  expect_is(convert_inputs(schools, "rubin", test_data = schools[7:8,]), "list")
})

