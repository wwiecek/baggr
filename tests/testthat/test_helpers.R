context("baggr helper functions")
library(baggr)
set.seed(1990)

test_that("prepare_ma()", {
  expect_error(prepare_ma(schools), "individual-level")
  expect_error(prepare_ma(microcredit_simplified), "no column")
  expect_error(prepare_ma(microcredit, outcome = "consumption"), "treatment values are NA")
  expect_error(prepare_ma(microcredit[!is.na(microcredit$treatment),], outcome = "consumption"),
               "outcome values are NA")

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

})

test_that("convert_inputs()", {
  # Rubin model
  expect_is(convert_inputs(schools, "rubin"), "list")
  expect_error(convert_inputs(schools, "mutau"))
  expect_is(convert_inputs(schools, "rubin", test_data = schools[7:8,]), "list")
})

