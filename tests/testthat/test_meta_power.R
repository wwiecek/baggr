context("meta_power")

library(baggr)

test_that("meta_power works with data frame input", {
  out <- meta_power(schools, n = 7, print_plot = FALSE)

  expect_true(is.list(out))
  expect_true(all(c("plot", "values") %in% names(out)))
  expect_equal(nrow(out$values$grid_wide), 49)
  expect_true("power_RE" %in% names(out$values$grid_wide))
  expect_false("power_FE" %in% names(out$values$grid_wide))
})

test_that("meta_power accepts baggr object and adds highlighted point", {
  set.seed(1999)
  bg <- baggr(schools,
              chains = 1,
              iter = 300,
              warmup = 150,
              refresh = 0,
              control = list(adapt_delta = 0.95))

  out <- meta_power(bg, n = 6, contours = c(), print_plot = FALSE)

  expect_true(inherits(out$plot, "ggplot"))

  built <- ggplot2::ggplot_build(out$plot)
  has_point_layer <- any(vapply(
    built$data,
    function(d) {
      ("shape" %in% names(d)) && any(d$shape == 21)
    },
    logical(1)
  ))

  expect_true(has_point_layer)
})
