library(baggr)

test_that("LOO plot labels match estimates, intervals and comparison style", {
  summaries <- rbind(c(-5, -6, 2, -12, -1, .95),
                     c(-3, -4, 1, -8, -2, .95))
  colnames(summaries) <- c("mean", "median", "sd", "lower", "upper", "width")
  local_mocked_bindings(
    treatment_effect = function(bg, summary) {
      tau <- summaries[bg,,drop = FALSE]
      attr(tau, "effect") <- "mean"
      list(theta = tau)
    },
    .package = "baggr"
  )
  cv <- structure(list(
    full_model = list(effects = "mean",
                      inputs = structure(list(), n_groups = 2,
                                         group_label = c("A", "B")),
                      data = data.frame(group = c("A", "B"),
                                       yi = c(-10, -2), sei = c(1, 0.5))),
    models = list(1, 2), pointwise = c(-30, -40)
  ), class = "baggr_cv")

  pl <- plot(cv)
  built <- ggplot2::ggplot_build(pl)
  labels <- built$data[[3]]
  expect_equal(nrow(pl$data), 4)
  expect_setequal(pl$data$group, c("A", "B"))
  expect_equal(pl$data$lower[1:2], c(-12, -8))
  expect_equal(pl$data$median[1:2], c(-6, -4))
  expect_equal(pl$data$upper[1:2], c(-1, -2))
  expect_equal(pl$data$lower[3:4], c(-11.96, -2.98))
  expect_equal(pl$data$upper[3:4], c(-8.04, -1.02))
  expect_equal(labels$label,
               baggr:::paste3_formatter(pl$data$median, pl$data$lower,
                                        pl$data$upper, 1)[c(1, 3, 2, 4)])
  expect_equal(labels$x, sort(built$data[[2]]$x))
  expect_equal(labels$size, rep(4, 4))
  expect_true(all(labels$hjust < 0))
  expect_equal(pl$theme$plot.margin,
               single_comp_plot(pl$data, add_values = TRUE)$theme$plot.margin)

  custom <- ggplot2::ggplot_build(plot(cv, values_digits = 3, values_size = 5))
  expect_equal(custom$data[[3]]$size, rep(5, 4))
  expect_equal(custom$data[[3]]$label,
               baggr:::paste3_formatter(pl$data$median, pl$data$lower,
                                        pl$data$upper, 3)[c(1, 3, 2, 4)])
  expect_length(plot(cv, add_values = FALSE)$layers, 2)

  cv$full_model$summary_data <- cv$full_model$data
  cv$full_model$data <- NULL
  expect_equal(plot(cv)$data, pl$data)
  cv$models <- NULL
  expect_error(plot(cv), "return_models = TRUE", fixed = TRUE)
})
