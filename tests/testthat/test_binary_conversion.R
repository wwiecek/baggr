library(baggr)
library(testthat)

binary_counts <- data.frame(group = c("A", "B"), a = c(10, 20),
                            n1 = c(100, 120), c = c(15, 25),
                            n2 = c(100, 120))

test_that("automatic binary conversion reports each selected measure once", {
  for(effect in c("logOR", "logRR", "RD")) {
    messages <- capture_messages(
      inputs <- convert_inputs(binary_counts, model = "rubin", effect = effect))
    expect_identical(messages, paste0(
      "Automatically summarising binary data with ", effect, ".\n"))
    expect_identical(attr(inputs, "effect"), effect)
    expect_equal(attr(inputs, "data"), prepare_ma(binary_counts, effect = effect))
  }
})

test_that("default binary conversion reports logOR once", {
  for(effect in list(NULL, "custom label")) {
    messages <- capture_messages(
      inputs <- convert_inputs(binary_counts, model = "rubin", effect = effect))
    expect_identical(messages,
                     "Automatically summarising binary data with logOR.\n")
    expect_identical(attr(inputs, "effect"), "logOR")
  }
})

test_that("silent binary conversion and prepared estimates emit no message", {
  for(effect in list(NULL, "logOR", "logRR", "RD")) {
    expect_message(convert_inputs(binary_counts, model = "rubin",
                                  effect = effect, silent = TRUE), NA)
  }
  prepared <- prepare_ma(binary_counts, effect = "logRR")
  expect_message(convert_inputs(prepared, model = "rubin",
                                effect = "custom label"), NA)
})
