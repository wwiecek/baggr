context("baggr() calls with full model")
library(baggr)
set.seed(1990)

# Generate 8 schools like IPD data
schools_ipd <- data.frame()
N <- c(rep(10, 4), rep(20, 4))
for(i in 1:8){
  x <- rnorm(N[i])
  x <- (x-mean(x))/sd(x)
  x <- x*schools$se[i]*sqrt(N[i])/1.41 + schools$tau[i]

  y <- rnorm(N[i])
  y <- (y-mean(y))/sd(y)
  y <- y*schools$se[i]*sqrt(N[i])/1.41

  schools_ipd <- rbind(schools_ipd,
                       data.frame(group = schools$group[i], outcome = x, treatment = 1),
                       # This is just so that we don't trip off prepare_ma:
                       data.frame(group = schools$group[i], outcome = y, treatment = 0))
}


bg_n <- expect_warning(baggr(schools_ipd, pooling = "none", iter = 150, refresh=0))
bg_p <- expect_warning(baggr(schools_ipd, pooling = "partial", iter = 150, refresh=0))
bg_f <- expect_warning(baggr(schools_ipd, pooling = "full", iter = 150, refresh=0))


test_that("Different pooling methods work for the full model", {
  expect_is(bg_n, "baggr")
  expect_is(bg_p, "baggr")
  expect_is(bg_f, "baggr")
})

test_that("Basic operations on full data model", {
  expect_error(baggr(schools_ipd, rubbish = 41))
  expect_is(pooling(bg_p)[,,1], "matrix")
  expect_is(plot(bg_p), "gg")
  expect_is(effect_plot(bg_p), "gg")
  expect_is(forest_plot(bg_p), "vpPath")
  bgc <- try(baggr_compare(bg_n, bg_p, bg_f))
  expect_is(bgc, "baggr_compare")
})

test_that("Full model crashes with nonsense inputs", {
  expect_error(baggr(schools_ipd, outcome = 2), "Arguments")
  expect_error(baggr(schools_ipd, group = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = "wrong"), "no column")
  schools_ipd2 <- schools_ipd; schools_ipd2$treatment <- as.character(schools_ipd$treatment)
  expect_error(baggr(schools_ipd2), "has to be numeric")
  schools_ipd2 <- schools_ipd; schools_ipd2$outcome <- as.character(schools_ipd$outcome)
  expect_error(baggr(schools_ipd2), "has to be numeric")

})

test_that("Using old syntax works still", {
  bg <- expect_message(baggr(schools_ipd,
                       model = "full",
                       pooling = "none", iter = 10, refresh=0))
  expect_is(bg, "baggr")
})


comp_flpl <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 150, what = "pooling"
))

comp_flpr <- expect_warning(baggr_compare(
  schools, model = "rubin", iter = 150, what = "prior"
))

comp_flmdls <- baggr_compare(bg_f, bg_p)

test_that("baggr comparison method works for Full model", {

  expect_is(comp_flpl, "baggr_compare")
  expect_is(comp_flpr, "baggr_compare")
  expect_is(comp_flmdls, "baggr_compare")

  expect_is(testthat::capture_output(print(comp_flpl)), "character")
  expect_is(testthat::capture_output(print(comp_flpr)), "character")
  expect_is(testthat::capture_output(print(comp_flmdls)), "character")

  expect_gt(length(comp_flpl), 0)
  expect_gt(length(comp_flpr), 0)
  expect_gt(length(comp_flmdls), 0)

  expect_is(plot(comp_flpl), "gg")

  expect_is(plot(comp_flpr), "ggplot")

  expect_is(plot(comp_flmdls), "gg")

  expect_is(plot(comp_flpl, grid_models = TRUE), "gtable")

  expect_is(plot(comp_flpr, grid_models = TRUE), "gtable")

  expect_is(plot(comp_flmdls, grid_models = TRUE), "gtable")
})

